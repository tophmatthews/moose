//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADLPSViscoPlasticityStressUpdate.h"

#include "libmesh/utility.h"

registerADMooseObject("TensorMechanicsApp", ADLPSViscoPlasticityStressUpdate);

defineADValidParams(
    ADLPSViscoPlasticityStressUpdate, ADStressUpdateBase, params.addClassDescription("");
    params += validParams<ADSingleVariableReturnMappingSolution<RESIDUAL>>();
    params.addParam<Real>("max_inelastic_increment",
                          1.0e-4,
                          "The maximum inelastic strain increment allowed in a time step");
    params.addParam<std::string>(
        "effective_inelastic_strain_name",
        "effective_creep_strain",
        "Name of the material property that stores the effective inelastic strain");
    params.addRequiredParam<std::vector<MaterialPropertyName>>(
        "coefficients", "Leading coefficient in power-law equation");
    params.addRequiredParam<std::vector<Real>>(
        "powers", "Exponent on effective stress in power-law equation");
    params.addParam<bool>("verbose", false, "Flag to output verbose information");
    params.addRangeCheckedParam<Real>(
        "initial_porosity", 0.0, "initial_porosity>0.0", "Initial porosity");

    params.addParamNamesToGroup("effective_inelastic_strain_name", "Advanced"););

template <ComputeStage compute_stage>
ADLPSViscoPlasticityStressUpdate<compute_stage>::ADLPSViscoPlasticityStressUpdate(
    const InputParameters & parameters)
  : ADStressUpdateBase<compute_stage>(parameters),
    ADSingleVariableReturnMappingSolution<compute_stage>(parameters),
    _effective_inelastic_strain(adDeclareADProperty<Real>(
        _base_name + adGetParam<std::string>("effective_inelastic_strain_name"))),
    _effective_inelastic_strain_old(adGetMaterialPropertyOld<Real>(
        _base_name + adGetParam<std::string>("effective_inelastic_strain_name"))),
    _creep_strain(adDeclareADProperty<RankTwoTensor>(_base_name + "creep_strain")),
    _creep_strain_old(adGetMaterialPropertyOld<RankTwoTensor>(_base_name + "creep_strain")),
    _porosity(adDeclareADProperty<Real>("porosity")),
    _porosity_old(adGetMaterialPropertyOld<Real>("porosity")),
    _max_inelastic_increment(adGetParam<Real>("max_inelastic_increment")),
    _powers(adGetParam<std::vector<Real>>("powers")),
    _num_models(_powers.size()),
    _coefficients(_num_models),
    _gauge_stresses(_num_models),
    _initial_porosity(adGetParam<Real>("initial_porosity")),
    _verbose(adGetParam<bool>("verbose")),
    _hydro_stress(0.0),
    _identity_two(RankTwoTensor::initIdentity),
    _dhydro_stress_dsigma(_identity_two / 3.0),
    _derivative(0.0),
    _current_n(0.0)
{
  std::vector<MaterialPropertyName> coeff_names =
      adGetParam<std::vector<MaterialPropertyName>>("coefficients");
  if (_num_models != coeff_names.size())
    paramError(
        "powers", "In ", _name, ": Number of powers and number of coefficients must be the same!");

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _coefficients[i] = &adGetADMaterialProperty<Real>(coeff_names[i]);
    _gauge_stresses[i] = &adDeclareADProperty<Real>(
        _base_name + "gauge_stress_n" + Moose::stringify(_powers[i]) + "_i" + Moose::stringify(i));
  }
}

template <ComputeStage compute_stage>
void
ADLPSViscoPlasticityStressUpdate<compute_stage>::initQpStatefulProperties()
{
  _effective_inelastic_strain[_qp] = 0.0;
  _porosity[_qp] = _initial_porosity;
  _creep_strain[_qp].zero();
}

template <ComputeStage compute_stage>
void
ADLPSViscoPlasticityStressUpdate<compute_stage>::propagateQpStatefulProperties()
{
  _effective_inelastic_strain[_qp] = _effective_inelastic_strain_old[_qp];
  _porosity[_qp] = _porosity_old[_qp];
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

template <ComputeStage compute_stage>
void
ADLPSViscoPlasticityStressUpdate<compute_stage>::updateState(
    ADRankTwoTensor & strain_increment,
    ADRankTwoTensor & inelastic_strain_increment,
    const ADRankTwoTensor & /*rotation_increment*/,
    ADRankTwoTensor & stress,
    const RankTwoTensor & /*stress_old*/,
    const ADRankFourTensor & elasticity_tensor,
    const RankTwoTensor & elastic_strain_old)
{
  // Compute initial hydrostatic stress and porosity
  _hydro_stress = stress.trace() / 3.0;
  _porosity[_qp] = _porosity_old[_qp];
  _effective_inelastic_strain[_qp] = _effective_inelastic_strain_old[_qp];
  inelastic_strain_increment.zero();

  // Compute swelling increment due to gas diffusion
  const ADReal scalar_swelling_increment = computeSwellingIncrement(_hydro_stress);

  // Compute new intermediate stress and scalar stresses
  if (scalar_swelling_increment)
  {
    inelastic_strain_increment += scalar_swelling_increment * _identity_two;
    strain_increment -= scalar_swelling_increment * _identity_two;
    stress = elasticity_tensor * (elastic_strain_old + strain_increment);
    _hydro_stress = stress.trace() / 3.0;
    _effective_inelastic_strain[_qp] += scalar_swelling_increment;
    _porosity[_qp] = (1.0 - _porosity_old[_qp]) * scalar_swelling_increment + _porosity_old[_qp];
  }

  // Compute intermediate equivalent stress
  const ADRankTwoTensor dev_stress = stress.deviatoric();
  const ADReal dev_stress_squared = dev_stress.doubleContraction(dev_stress);
  const ADReal equiv_stress = MooseUtils::absoluteFuzzyEqual(dev_stress_squared, 0.0)
                                  ? 0.0
                                  : std::sqrt(3.0 / 2.0 * dev_stress_squared);

  computeStressInitialize(equiv_stress, elasticity_tensor);

  if (equiv_stress || _hydro_stress)
  {
    // Initialize sum of all derivatives of residual with respect to each gauge stress
    ADReal sum_dpsi_dgauge(0);

    ADRankTwoTensor creep_strain_increment;
    for (unsigned int i = 0; i < _num_models; ++i)
    {
      ADReal dpsi_dgauge(0); // Initalize derivative of residual with respect to the gauge stress
      computeCreepStrainIncrementN((*_gauge_stresses[i])[_qp],
                         dpsi_dgauge,
                         creep_strain_increment,
                         equiv_stress,
                         dev_stress,
                         stress,
                         i);
      sum_dpsi_dgauge += dpsi_dgauge;
    }

    inelastic_strain_increment += creep_strain_increment;
    strain_increment -= creep_strain_increment;
    _effective_inelastic_strain[_qp] += sum_dpsi_dgauge * _dt;
    stress = elasticity_tensor * (elastic_strain_old + strain_increment);
    _creep_strain[_qp] = _creep_strain_old[_qp] + inelastic_strain_increment;

    _porosity[_qp] =
        (1.0 - _porosity_old[_qp]) * inelastic_strain_increment.tr() + _porosity_old[_qp];
  }

  if (_verbose)
    Moose::out << "current porosity: " << MetaPhysicL::raw_value(_porosity[_qp])
               << " old porosity: " << _porosity_old[_qp] << std::endl;

  if (_porosity[_qp] < 0.0 || _porosity[_qp] > 1.0)
  {
    std::stringstream exception_meessage;
    exception_meessage << "In " << _name << ": porosity (" << _porosity[_qp]
                       << ") is less than zero or greater than one. Cutting timestep";
    throw MooseException(exception_meessage.str());
  }

  computeStressFinalize(inelastic_strain_increment);
}

template <ComputeStage compute_stage>
Real
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeReferenceResidual(
    const ADReal & /*effective_trial_stress*/, const ADReal & scalar_effective_inelastic_strain)
{
  return MetaPhysicL::raw_value(scalar_effective_inelastic_strain);
}

template <ComputeStage compute_stage>
Real
ADLPSViscoPlasticityStressUpdate<compute_stage>::maximumPermissibleValue(
    const ADReal & effective_trial_stress) const
{
  return MetaPhysicL::raw_value(effective_trial_stress) * 1.0e6; // todo check
}

template <ComputeStage compute_stage>
Real
ADLPSViscoPlasticityStressUpdate<compute_stage>::minimumPermissibleValue(
    const ADReal & effective_trial_stress) const
{
  return MetaPhysicL::raw_value(effective_trial_stress) / 1.0e6; // todo check
}

template <ComputeStage compute_stage>
Real
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeTimeStepLimit()
{
  const Real scalar_inelastic_strain_incr =
      MetaPhysicL::raw_value(_effective_inelastic_strain[_qp]) -
      _effective_inelastic_strain_old[_qp];

  if (MooseUtils::absoluteFuzzyEqual(scalar_inelastic_strain_incr, 0.0))
    return std::numeric_limits<Real>::max();

  return _dt * _max_inelastic_increment / scalar_inelastic_strain_incr;
}

template <ComputeStage compute_stage>
void
ADLPSViscoPlasticityStressUpdate<compute_stage>::outputIterationSummary(
    std::stringstream * iter_output, const unsigned int total_it)
{
  if (iter_output)
  {
    *iter_output << "At element " << _current_elem->id() << " _qp=" << _qp << " Coordinates "
                 << _q_point[_qp] << " block=" << _current_elem->subdomain_id() << '\n';
  }
  ADSingleVariableReturnMappingSolution<compute_stage>::outputIterationSummary(iter_output,
                                                                               total_it);
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeResidual(const ADReal & equiv_stress,
                                                                 const ADReal & trial_gauge)
{
  const ADReal M = computeM(_hydro_stress, trial_gauge);
  const ADReal h = computeH(_current_n, M);
  const Real n_mp = (_current_n - 1.0) / (_current_n + 1.0);
  const ADReal res = Utility::pow<2>(equiv_stress / trial_gauge) + _porosity[_qp] * (h + n_mp / h) -
                     1.0 - n_mp * Utility::pow<2>(_porosity[_qp]);

  const ADReal dM_dgauge_stress = -M / trial_gauge;
  const ADReal dh_dM = computeH(_current_n, M, true);
  const ADReal dres_dh = _porosity[_qp] * (1.0 - n_mp / Utility::pow<2>(h));
  const ADReal dres_dgauge_stress_left =
      -2.0 * Utility::pow<2>(trial_gauge) / Utility::pow<3>(trial_gauge);
  _derivative = dres_dgauge_stress_left + dres_dh * dh_dM * dM_dgauge_stress;

  if (_verbose)
  {
    Moose::out << "in computeResidual:\n"
               << "  pos: " << _q_point[_qp] << " hydro: " << MetaPhysicL::raw_value(_hydro_stress)
               << " equiv: " << MetaPhysicL::raw_value(equiv_stress)
               << " gauge: " << MetaPhysicL::raw_value(trial_gauge)
               << " M: " << MetaPhysicL::raw_value(M) << " h: " << MetaPhysicL::raw_value(h)
               << std::endl;
    Moose::out << "  res: " << MetaPhysicL::raw_value(res)
               << "  deriv: " << MetaPhysicL::raw_value(_derivative) << std::endl;
  }

  return res;
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeH(const Real n,
                                                          const ADReal & M,
                                                          const bool derivative)
{
  if (derivative)
    return (n + 1.0) / n * std::pow(M, 1.0 / n) *
           std::pow(1.0 + std::pow(M, (n + 1.0) / n) / n, n - 1.0);
  return std::pow(1.0 + std::pow(M, (n + 1.0) / n) / n, n);
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeM(const ADReal & hydro_stress,
                                                          const ADReal & gauge_stress)
{
  return 3.0 / 2.0 * std::abs(hydro_stress) / gauge_stress;
}

template <ComputeStage compute_stage>
ADRankTwoTensor
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeDGaugeDSigma(
    const ADReal & gauge_stress,
    const ADReal & equiv_stress,
    const ADRankTwoTensor & dev_stress,
    const ADRankTwoTensor & stress,
    const Real n)
{
  const ADReal M = computeM(_hydro_stress, gauge_stress);
  const ADReal h = computeH(n, M);

  ADReal dres_dhydro_stress = 0.0;
  if (_hydro_stress)
  {
    const ADReal dres_dh = _porosity[_qp] * (1.0 - (n - 1.0) / (n + 1.0) / Utility::pow<2>(h));
    const ADReal dh_dM = computeH(n, M, true);
    const ADReal dM_dhydro_stress = M / _hydro_stress;
    dres_dhydro_stress = dres_dh * dh_dM * dM_dhydro_stress;
  }

  ADReal dres_dequiv_stress = 0.0;
  ADRankTwoTensor dequiv_stress_dsigma;
  if (equiv_stress)
  {
    dres_dequiv_stress = 2.0 * equiv_stress / Utility::pow<2>(gauge_stress);
    dequiv_stress_dsigma = 3.0 / 2.0 * dev_stress / equiv_stress;
  }

  const ADRankTwoTensor dres_dsigma =
      dres_dhydro_stress * _dhydro_stress_dsigma + dres_dequiv_stress * dequiv_stress_dsigma;
  const ADRankTwoTensor dgauge_dsigma =
      gauge_stress * dres_dsigma / dres_dsigma.doubleContraction(stress);

  return dgauge_dsigma;
}

template <ComputeStage compute_stage>
void
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeCreepStrainIncrementN(
    ADReal & gauge_stress,
    ADReal & dpsi_dgauge,
    ADRankTwoTensor & creep_strain_increment,
    const ADReal & equiv_stress,
    const ADRankTwoTensor & dev_stress,
    const ADRankTwoTensor & stress,
    const unsigned int i)
{
  // Set model parameters for evaluations
  _current_n = _powers[i];

  // Run non-linear solve for gauge stress
  if (_hydro_stress && _porosity[_qp])
    returnMappingSolve(equiv_stress, gauge_stress, _console);
  else if (equiv_stress)
    gauge_stress =
        equiv_stress / std::sqrt(1.0 + (_current_n - 1.0) / (_current_n + 1.0) * _porosity[_qp]);
  else
    mooseError(
        "In ",
        _name,
        ": No equivalent and hydro stress is zero. Something is wrong in computeCreepStrainIncrementN");

  if (gauge_stress < 0.0)
    mooseError("In ",
               _name,
               ": Gauge stress is less than zero. Something is wrong with the inner Newton solve");

  dpsi_dgauge = (*_coefficients[i])[_qp] * std::pow(gauge_stress, _current_n);

  creep_strain_increment +=
      _dt * dpsi_dgauge *
      computeDGaugeDSigma(gauge_stress, equiv_stress, dev_stress, stress, _current_n);
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::initialGuess(const ADReal & effective_trial_stress)
{
  if (effective_trial_stress)
    return effective_trial_stress;
  else if (_hydro_stress)
    return std::abs(_hydro_stress);
  else
    mooseError("negative initial guess");
}

adBaseClass(ADLPSViscoPlasticityStressUpdate);
