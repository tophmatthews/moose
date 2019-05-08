//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADLPSViscoPlasticityStressUpdate.h"

#include "MooseMesh.h"
#include "ElasticityTensorTools.h"
#include "libmesh/utility.h"
#include "RankTwoScalarTools.h"

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
    params.addRequiredCoupledVar("temperature", "Coupled temperature");
    params.addRequiredRangeCheckedParam<Real>("coefficient",
                                              "coefficient>=0.0",
                                              "Leading coefficient in power-law equation");
    params.addRequiredRangeCheckedParam<Real>("n_exponent",
                                              "n_exponent>0.0",
                                              "Exponent on effective stress in power-law equation");
    params.addRequiredRangeCheckedParam<Real>("activation_energy",
                                              "activation_energy>0",
                                              "Activation energy");
    params.addParam<Real>("gas_constant", 8.3143, "Universal gas constant");
    params.addParam<bool>("verbose", false, "Flag to output verbose information");
    params.addRangeCheckedParam<Real>(
        "initial_porosity", 0.0, "initial_porosity>=0.0", "Initial porosity");
    params.addRangeCheckedParam<unsigned int>(
        "max_its",
        50,
        "max_its>0 & max_its<10000",
        "Maximum number of bi-section iterations for gauge stress calculation");
    params.addParam<bool>("use_newton", true, "Flag to use newton method");
    params.addParam<Real>("abs_tol", 1e-6, "absolute_tolerance");

    params.addParamNamesToGroup("effective_inelastic_strain_name", "Advanced"););

template <ComputeStage compute_stage>
ADLPSViscoPlasticityStressUpdate<compute_stage>::ADLPSViscoPlasticityStressUpdate(
    const InputParameters & parameters)
  : ADStressUpdateBase<compute_stage>(parameters),
    ADSingleVariableReturnMappingSolution<compute_stage>(parameters),
    _three_shear_modulus(1.0),
    _effective_inelastic_strain(adDeclareADProperty<Real>(
        _base_name + adGetParam<std::string>("effective_inelastic_strain_name"))),
    _effective_inelastic_strain_old(adGetMaterialPropertyOld<Real>(
        _base_name + adGetParam<std::string>("effective_inelastic_strain_name"))),
    _porosity(adDeclareADProperty<Real>("porosity")),
    _porosity_old(adGetMaterialPropertyOld<Real>("porosity")),
    _max_inelastic_increment(adGetParam<Real>("max_inelastic_increment")),
    _temperature(adCoupledValue("temperature")),
    _coefficient(adGetParam<Real>("coefficient")),
    _n(adGetParam<Real>("n_exponent")),
    _activation_energy(adGetParam<Real>("activation_energy")),
    _gas_constant(adGetParam<Real>("gas_constant")),
    _derivative(0.0),
    _creep_strain(adDeclareADProperty<RankTwoTensor>(_base_name + "creep_strain")),
    _creep_strain_old(adGetMaterialPropertyOld<RankTwoTensor>(_base_name + "creep_strain")),
    _verbose(adGetParam<bool>("verbose")),
    _initial_porosity(adGetParam<Real>("initial_porosity")),
    _max_its(adGetParam<unsigned int>("max_its")),
    _identity_two(RankTwoTensor::initIdentity),
    _dhydro_stress_dsigma(_identity_two / 3.0),
    _gauge_stress(adDeclareADProperty<Real>("gauge_stress")),
    _newton(adGetParam<bool>("use_newton")),
    _abs_tol(adGetParam<Real>("abs_tol"))
{
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
ADLPSViscoPlasticityStressUpdate<compute_stage>::propagateQpStatefulPropertiesRadialReturn()
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
  inelastic_strain_increment.zero();
  ADReal dpsi_dgauge = 0.0;

  ADRankTwoTensor dev_stress = stress.deviatoric();
  ADReal dev_trial_stress_squared = dev_stress.doubleContraction(dev_stress);
  ADReal equiv_stress = MooseUtils::absoluteFuzzyEqual(dev_trial_stress_squared, 0.0)
                            ? 0.0
                            : std::sqrt(3.0 / 2.0 * dev_trial_stress_squared);

  if (equiv_stress)
  {
    // Hydrostatic stress
    const ADReal hydro_stress = stress.trace() / 3.0;

    // Gauge trial stress that will be solved for
    _gauge_stress[_qp] = equiv_stress;

    computeGaugeStress(hydro_stress, equiv_stress, _gauge_stress[_qp]);

    // power law constants
    const Real sigma_0 = 1.0;
    const ADReal epsilon_0 =
        _coefficient * std::exp(-_activation_energy / _gas_constant / _temperature[_qp]);

    if (_gauge_stress[_qp] < 0.0)
      mooseError("In ", _name, ": internal error. gauge stress is negative");
    dpsi_dgauge = epsilon_0 / sigma_0 * std::pow(_gauge_stress[_qp], _n);

    const ADReal M = 2.0 / 3.0 * std::abs(hydro_stress / _gauge_stress[_qp]);
    const ADReal h = computeH(M);

    const ADReal dres_dh =
        2.0 * _porosity[_qp] * (1.0 - (_n - 1.0) / (_n + 1.0) / Utility::pow<2>(h));
    const ADReal dh_dM = computeH(M, true);
    const ADReal dM_dhydro_stress = M / hydro_stress;
    const ADReal dres_dhydro_stress = dres_dh * dh_dM * dM_dhydro_stress;

    const ADReal dres_dequiv_stress = 2.0 * equiv_stress / Utility::pow<2>(_gauge_stress[_qp]);
    const ADRankTwoTensor dequiv_stress_dsigma =
        1.5 / equiv_stress * (stress - hydro_stress * _identity_two);

    const ADRankTwoTensor dres_dsigma =
        dres_dhydro_stress * _dhydro_stress_dsigma + dres_dequiv_stress * dequiv_stress_dsigma;
    const ADRankTwoTensor dgauge_dsigma =
        _gauge_stress[_qp] * dres_dsigma / dres_dsigma.doubleContraction(stress);

    const ADRankTwoTensor strain_rate = dpsi_dgauge * dgauge_dsigma;

    if (dpsi_dgauge)
      inelastic_strain_increment = strain_rate * _dt;
    strain_increment -= inelastic_strain_increment;
    _effective_inelastic_strain[_qp] = _effective_inelastic_strain_old[_qp] + dpsi_dgauge;
  }

  stress = elasticity_tensor * (elastic_strain_old + strain_increment);

  _porosity[_qp] = (1.0 - _porosity_old[_qp]) * strain_increment.tr() + _porosity_old[_qp];
  _creep_strain[_qp] = _creep_strain_old[_qp] + strain_increment;

  if (_verbose)
    Moose::out << "new porosity: " << _porosity[_qp] << " old porosity: " << _porosity_old[_qp]
               << std::endl;
}

template <ComputeStage compute_stage>
void
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeGaugeStress(const ADReal & hydro_stress,
                                                                    const ADReal & equiv_stress,
                                                                    ADReal & gauge_stress)
{
  if (!gauge_stress)
    return;

  if (_newton)
  {
    _hydro_stress = hydro_stress;
    returnMappingSolve(equiv_stress, gauge_stress, _console);
    return;
  }

  computeBisectionSolve(hydro_stress, equiv_stress, gauge_stress);
}

template <ComputeStage compute_stage>
void
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeBisectionSolve(
    const ADReal & hydro_stress, const ADReal & equiv_stress, ADReal & gauge_stress)
{
  const ADReal initial_R = computeGaugeStressResidual(hydro_stress, equiv_stress, gauge_stress);

  if (_verbose)
    Moose::out << "initial_R: " << MetaPhysicL::raw_value(initial_R)
               << " gauge_stress: " << MetaPhysicL::raw_value(gauge_stress) << std::endl;

  ADReal lower_gauge_bound = gauge_stress;
  ADReal upper_gauge_bound = gauge_stress;

  Real factor = 1.1;
  if (initial_R < 0)
    factor = 1.0 / factor;

  ADReal bounding_R;
  ADReal trial_gauge_bound = gauge_stress;
  do
  {
    trial_gauge_bound *= factor;
    bounding_R = computeGaugeStressResidual(hydro_stress, equiv_stress, trial_gauge_bound);
  } while (std::abs(MetaPhysicL::raw_value(trial_gauge_bound)) < 1.0e20 &&
           MetaPhysicL::raw_value(bounding_R) * MetaPhysicL::raw_value(initial_R) > 0);

  if (MetaPhysicL::raw_value(initial_R) < 0)
    upper_gauge_bound = trial_gauge_bound;
  else
    lower_gauge_bound = trial_gauge_bound;

  if (_verbose)
    Moose::out << "Bounding calcualations: gauge_stress: " << MetaPhysicL::raw_value(gauge_stress)
               << " initial_R: " << MetaPhysicL::raw_value(initial_R)
               << " lower_gauge_bound: " << MetaPhysicL::raw_value(lower_gauge_bound)
               << " lower bound R: "
               << MetaPhysicL::raw_value(
                      computeGaugeStressResidual(hydro_stress, equiv_stress, lower_gauge_bound))
               << " upper_gauge_bound: " << MetaPhysicL::raw_value(upper_gauge_bound)
               << " upper bound R: "
               << MetaPhysicL::raw_value(
                      computeGaugeStressResidual(hydro_stress, equiv_stress, upper_gauge_bound))
               << std::endl;

  if (trial_gauge_bound >= 1.0e20)
  {
    mooseWarning(
        "In ",
        _name,
        ": internal error. Second residual for gauge stress calculation could not be found");
    throw MooseException("cutting timestep");
  }

  ADReal trial_gauge;
  ADReal R;
  unsigned int its = 0;
  do
  {
    ++its;

    trial_gauge = (upper_gauge_bound + lower_gauge_bound) / 2.0;
    R = computeGaugeStressResidual(hydro_stress, equiv_stress, trial_gauge);

    if (R < 0.0)
      lower_gauge_bound = trial_gauge;
    else if (R > 0.0)
      upper_gauge_bound = trial_gauge;
    else
      mooseError("In ",
                 _name,
                 ": internal error. Something is wrong with the bi-section calculation for the "
                 "gauge stress");

    if (_verbose)
      Moose::out << "   it: " << its << " equiv_stress: " << MetaPhysicL::raw_value(equiv_stress)
                 << " trial_gauge: " << MetaPhysicL::raw_value(trial_gauge)
                 << " porosity: " << MetaPhysicL::raw_value(_porosity[_qp]) << " R "
                 << MetaPhysicL::raw_value(R)
                 << " upper_gauge_bound: " << MetaPhysicL::raw_value(upper_gauge_bound)
                 << " lower_gauge_bound: " << MetaPhysicL::raw_value(lower_gauge_bound)
                 << std::endl;

  } while (its < _max_its && std::abs(MetaPhysicL::raw_value(R)) > _abs_tol);

  if (its >= _max_its)
    mooseError("Unable to converge on gauge stress");

  gauge_stress = trial_gauge;
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeGaugeStressResidual(
    const ADReal & hydro_stress, const ADReal & equiv_stress, const ADReal & gauge_stress)
{
  const ADReal M = 2.0 / 3.0 * std::abs(hydro_stress / gauge_stress);
  const ADReal h = computeH(M);
  const Real n_mp = (_n - 1.0) / (_n + 1.0);
  if (_verbose)
  {
    Moose::out << "in computeGaugeStressResidual:\n hydro: " << MetaPhysicL::raw_value(hydro_stress)
               << " equiv: " << MetaPhysicL::raw_value(equiv_stress)
               << " gauge: " << MetaPhysicL::raw_value(gauge_stress) << " M: " << M
               << " h: " << MetaPhysicL::raw_value(h) << std::endl;
    Moose::out << " res: "
               << MetaPhysicL::raw_value(Utility::pow<2>(equiv_stress / gauge_stress) +
                                         2.0 * _porosity[_qp] * (h + n_mp / h) - 1.0 -
                                         n_mp * Utility::pow<2>(_porosity[_qp]))
               << std::endl;
  }

  return Utility::pow<2>(equiv_stress / gauge_stress) + 2.0 * _porosity[_qp] * (h + n_mp / h) -
         1.0 - n_mp * Utility::pow<2>(_porosity[_qp]);
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeGaugeStressDerivative(
    const ADReal & hydro_stress, const ADReal & equiv_stress, const ADReal & gauge_stress)
{
  const ADReal M = 2.0 / 3.0 * std::abs(hydro_stress / gauge_stress);
  const ADReal dM_dgauge_stress = -M / gauge_stress;

  const ADReal h = computeH(M);
  const ADReal dh_dM = computeH(M, true);

  const ADReal dres_dh =
      2.0 * _porosity[_qp] * (1.0 - (_n - 1.0) / (_n + 1.0) / Utility::pow<2>(h));

  const ADReal dres_dgauge_stress_left =
      -2.0 * Utility::pow<2>(equiv_stress) / Utility::pow<3>(gauge_stress);

  return dres_dgauge_stress_left + dres_dh * dh_dM * dM_dgauge_stress;
}

template <ComputeStage compute_stage>
Real
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeReferenceResidual(
    const ADReal & effective_trial_stress, const ADReal & /*scalar_effective_inelastic_strain*/)
{
  return MetaPhysicL::raw_value(effective_trial_stress);
}

template <ComputeStage compute_stage>
Real
ADLPSViscoPlasticityStressUpdate<compute_stage>::maximumPermissibleValue(
    const ADReal & effective_trial_stress) const
{
  return MetaPhysicL::raw_value(effective_trial_stress) * 1.0e6;
}

template <ComputeStage compute_stage>
Real
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeTimeStepLimit()
{
  Real scalar_inelastic_strain_incr = MetaPhysicL::raw_value(_effective_inelastic_strain[_qp]) -
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
  return computeGaugeStressResidual(_hydro_stress, equiv_stress, trial_gauge);
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeDerivative(const ADReal & equiv_stress,
                                                                   const ADReal & trial_gauge)
{
  return computeGaugeStressDerivative(_hydro_stress, equiv_stress, trial_gauge);
}

template <ComputeStage compute_stage>
ADReal
ADLPSViscoPlasticityStressUpdate<compute_stage>::computeH(const ADReal & M, const bool derivative)
{
  if (derivative)
    return (_n + 1.0) / _n * std::pow(M, 1.0 / _n) *
           std::pow(1.0 + std::pow(M, (_n + 1.0) / _n) / _n, _n - 1.0);
  return std::pow(1.0 + std::pow(M, (_n + 1.0) / _n) / _n, _n);
}
