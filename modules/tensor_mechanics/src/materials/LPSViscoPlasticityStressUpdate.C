//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LPSViscoPlasticityStressUpdate.h"

#include "MooseMesh.h"
#include "ElasticityTensorTools.h"
#include "libmesh/utility.h"
#include "RankTwoScalarTools.h"

registerMooseObject("TensorMechanicsApp", LPSViscoPlasticityStressUpdate);

template <>
InputParameters
validParams<LPSViscoPlasticityStressUpdate>()
{
  InputParameters params = validParams<StressUpdateBase>();
  params += validParams<SingleVariableReturnMappingSolution>();
  params.addParam<Real>("max_inelastic_increment",
                        1.0e-4,
                        "The maximum inelastic strain increment allowed in a time step");
  params.addParam<std::string>(
      "effective_inelastic_strain_name",
      "effective_creep_strain",
      "Name of the material property that stores the effective inelastic strain");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "coefficients", "Leading coefficient in power-law equation");
  params.addRequiredParam<std::vector<Real>>("powers",
                                             "Exponent on effective stress in power-law equation");
  params.addParam<bool>("verbose", false, "Flag to output verbose information");
  params.addRangeCheckedParam<Real>(
      "initial_porosity", 0.0, "initial_porosity>0.0", "Initial porosity");

  params.addParamNamesToGroup("effective_inelastic_strain_name", "Advanced");
  return params;
}

LPSViscoPlasticityStressUpdate::LPSViscoPlasticityStressUpdate(const InputParameters & parameters)
  : StressUpdateBase(parameters),
    SingleVariableReturnMappingSolution(parameters),
    _effective_inelastic_strain(declareProperty<Real>(
        _base_name + getParam<std::string>("effective_inelastic_strain_name"))),
    _effective_inelastic_strain_old(getMaterialPropertyOld<Real>(
        _base_name + getParam<std::string>("effective_inelastic_strain_name"))),
    _creep_strain(declareProperty<RankTwoTensor>(_base_name + "creep_strain")),
    _creep_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "creep_strain")),
    _porosity(declareProperty<Real>("porosity")),
    _porosity_old(getMaterialPropertyOld<Real>("porosity")),
    _max_inelastic_increment(getParam<Real>("max_inelastic_increment")),
    _powers(getParam<std::vector<Real>>("powers")),
    _num_models(_powers.size()),
    _coefficients(_num_models),
    _gauge_stresses(_num_models),
    _initial_porosity(getParam<Real>("initial_porosity")),
    _verbose(getParam<bool>("verbose")),
    _hydro_stress(0.0),
    _identity_two(RankTwoTensor::initIdentity),
    _dhydro_stress_dsigma(_identity_two / 3.0),
    _derivative(0.0),
    _current_n(0.0),
    _current_activation_energy(0.0),
    _current_coefficient(0.0)
{
  std::vector<MaterialPropertyName> coeff_names =
      getParam<std::vector<MaterialPropertyName>>("coefficients");
  if (_num_models != coeff_names.size())
    paramError(
        "powers", "In ", _name, ": Number of powers and number of coefficients must be the same!");
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _coefficients[i] = &adGetADMaterialProperty<Real>(coeff_names[i]);
    _gauge_stresses[i] = &declareProperty<Real>(
        _base_name + "gauge_stress_n" + Moose::stringify(_powers[i]) + "_i" + Moose::stringify(i));
  }
}

void
LPSViscoPlasticityStressUpdate::initQpStatefulProperties()
{
  _effective_inelastic_strain[_qp] = 0.0;
  _porosity[_qp] = _initial_porosity;
  _creep_strain[_qp].zero();
}

void
LPSViscoPlasticityStressUpdate::propagateQpStatefulPropertiesRadialReturn()
{
  _effective_inelastic_strain[_qp] = _effective_inelastic_strain_old[_qp];
  _porosity[_qp] = _porosity_old[_qp];
  _creep_strain[_qp] = _creep_strain_old[_qp];
}

void
LPSViscoPlasticityStressUpdate::updateState(RankTwoTensor & strain_increment,
                                            RankTwoTensor & inelastic_strain_increment,
                                            const RankTwoTensor & /*rotation_increment*/,
                                            RankTwoTensor & stress,
                                            const RankTwoTensor & /*stress_old*/,
                                            const RankFourTensor & elasticity_tensor,
                                            const RankTwoTensor & elastic_strain_old,
                                            bool /*compute_full_tangent_operator*/,
                                            RankFourTensor & /*tangent_operator*/)
{
  inelastic_strain_increment.zero();

  RankTwoTensor dev_stress = stress.deviatoric();
  Real dev_trial_stress_squared = dev_stress.doubleContraction(dev_stress);
  Real equiv_stress = MooseUtils::absoluteFuzzyEqual(dev_trial_stress_squared, 0.0)
                          ? 0.0
                          : std::sqrt(3.0 / 2.0 * dev_trial_stress_squared);

  if (equiv_stress)
  {
    _hydro_stress = stress.trace() / 3.0;

    RankTwoTensor strain_rate;
    strain_rate.zero();
    Real dpsi_dgauge(0);
    Real sum_dpsi_dgauge(0);
    for (unsigned int i = 0; i < _num_models; ++i)
    {
      computeNStrainRate(
          (*_gauge_stresses[i])[_qp], dpsi_dgauge, strain_rate, equiv_stress, dev_stress, stress, i);
      sum_dpsi_dgauge += dpsi_dgauge;
    }

    inelastic_strain_increment = strain_rate * _dt;
    strain_increment -= inelastic_strain_increment;
    _effective_inelastic_strain[_qp] = _effective_inelastic_strain_old[_qp] + sum_dpsi_dgauge * _dt;
  }

  stress = elasticity_tensor * (elastic_strain_old + strain_increment);

  _porosity[_qp] =
      (1.0 - _porosity_old[_qp]) * inelastic_strain_increment.tr() + _porosity_old[_qp];

  if (_verbose)
    Moose::out << "new porosity: " << (_porosity[_qp]) << " old porosity: " << _porosity_old[_qp]
               << std::endl;

  if (_porosity[_qp] < 0.0 || _porosity[_qp] > 1.0)
  {
    std::stringstream exception_meessage;
    exception_meessage << "In " << _name << ": porosity (" << _porosity[_qp]
                       << ") is less than zero or greater than one. Cutting timestep";
    throw MooseException(exception_meessage.str());
  }

  _creep_strain[_qp] = _creep_strain_old[_qp] + strain_increment;
}

Real
LPSViscoPlasticityStressUpdate::computeReferenceResidual(
    const Real /*effective_trial_stress*/, const Real scalar_effective_inelastic_strain)
{
  return (scalar_effective_inelastic_strain);
}

Real
LPSViscoPlasticityStressUpdate::maximumPermissibleValue(const Real effective_trial_stress) const
{
  return (effective_trial_stress)*1.0e6;
}

Real
LPSViscoPlasticityStressUpdate::minimumPermissibleValue(const Real effective_trial_stress) const
{
  return (effective_trial_stress) / 1.0e6;
}

Real
LPSViscoPlasticityStressUpdate::computeTimeStepLimit()
{
  // Real scalar_inelastic_strain_incr = (_porosity[_qp]) -
  // _porosity_old[_qp];
  const Real scalar_inelastic_strain_incr =
      (_effective_inelastic_strain[_qp]) - _effective_inelastic_strain_old[_qp];
  if (MooseUtils::absoluteFuzzyEqual(scalar_inelastic_strain_incr, 0.0))
    return std::numeric_limits<Real>::max();

  return _dt * _max_inelastic_increment / scalar_inelastic_strain_incr;
}

void
LPSViscoPlasticityStressUpdate::outputIterationSummary(std::stringstream * iter_output,
                                                       const unsigned int total_it)
{
  if (iter_output)
  {
    *iter_output << "At element " << _current_elem->id() << " _qp=" << _qp << " Coordinates "
                 << _q_point[_qp] << " block=" << _current_elem->subdomain_id() << '\n';
  }
  SingleVariableReturnMappingSolution::outputIterationSummary(iter_output, total_it);
}

Real
LPSViscoPlasticityStressUpdate::computeResidual(const Real equiv_stress, const Real trial_gauge)
{
  const Real M = computeM(_hydro_stress, trial_gauge);
  const Real h = computeH(_current_n, M);
  const Real n_mp = (_current_n - 1.0) / (_current_n + 1.0);
  const Real res = Utility::pow<2>(equiv_stress / trial_gauge) +
                   2.0 * _porosity_old[_qp] * (h + n_mp / h) - 1.0 -
                   n_mp * Utility::pow<2>(_porosity_old[_qp]);

  const Real dM_dgauge_stress = computeM(_hydro_stress, trial_gauge, 2);
  const Real dh_dM = computeH(_current_n, M, true);
  const Real dres_dh = 2.0 * _porosity_old[_qp] * (1.0 - n_mp / Utility::pow<2>(h));
  const Real dres_dgauge_stress_left =
      -2.0 * Utility::pow<2>(trial_gauge) / Utility::pow<3>(trial_gauge);
  _derivative = dres_dgauge_stress_left + dres_dh * dh_dM * dM_dgauge_stress;

  if (_verbose)
  {
    Moose::out << "in computeResidual:\n"
               << "  pos: " << _q_point[_qp] << " hydro: " << (_hydro_stress)
               << " equiv: " << (equiv_stress) << " gauge: " << (trial_gauge) << " M: " << (M)
               << " h: " << (h) << std::endl;
    Moose::out << "  res: " << (res) << "  deriv: " << (_derivative) << std::endl;
  }

  return res;
}

Real
LPSViscoPlasticityStressUpdate::computeDerivative(const Real /*equiv_stress*/,
                                                  const Real /*trial_gauge*/)
{
  return _derivative;
}

Real
LPSViscoPlasticityStressUpdate::computeH(const Real n, const Real & M, const bool derivative)
{
  if (derivative)
    return (n + 1.0) / n * std::pow(M, 1.0 / n) *
           std::pow(1.0 + std::pow(M, (n + 1.0) / n) / n, n - 1.0);
  return std::pow(1.0 + std::pow(M, (n + 1.0) / n) / n, n);
}

Real
LPSViscoPlasticityStressUpdate::computeM(const Real & hydro_stress,
                                         const Real & gauge_stress,
                                         const unsigned int derivative)
{
  if (derivative == 0)
    return 3.0 / 2.0 * std::abs(hydro_stress) / gauge_stress;
  else if (derivative == 1)
    return 3.0 / 2.0 * std::abs(hydro_stress) / gauge_stress / hydro_stress;
  else if (derivative == 2)
    return -3.0 / 2.0 * std::abs(hydro_stress) / gauge_stress / gauge_stress;
  else
    mooseError("In ", _name, ": Internal error in LPSViscoPlasticityStressUpdate::computeM");
}

RankTwoTensor
LPSViscoPlasticityStressUpdate::computeDGaugeDSigma(const Real & gauge_stress,
                                                    const Real & equiv_stress,
                                                    const RankTwoTensor & dev_stress,
                                                    const RankTwoTensor & stress,
                                                    const Real n)
{
  const Real M = computeM(_hydro_stress, gauge_stress);
  const Real h = computeH(n, M);

  const Real dres_dh = _porosity_old[_qp] * (1.0 - (n - 1.0) / (n + 1.0) / Utility::pow<2>(h));
  const Real dh_dM = computeH(n, M, true);
  const Real dM_dhydro_stress = computeM(_hydro_stress, gauge_stress, 1);
  const Real dres_dhydro_stress = dres_dh * dh_dM * dM_dhydro_stress;

  const Real dres_dequiv_stress = 2.0 * equiv_stress / Utility::pow<2>(gauge_stress);
  const RankTwoTensor dequiv_stress_dsigma = 3.0 / 2.0 * dev_stress / equiv_stress;

  const RankTwoTensor dres_dsigma =
      dres_dhydro_stress * _dhydro_stress_dsigma + dres_dequiv_stress * dequiv_stress_dsigma;
  const RankTwoTensor dgauge_dsigma =
      gauge_stress * dres_dsigma / dres_dsigma.doubleContraction(stress);

  return dgauge_dsigma;
}

void
LPSViscoPlasticityStressUpdate::computeNStrainRate(Real & gauge_stress,
                                                   Real & dpsi_dgauge,
                                                   RankTwoTensor & strain_rate,
                                                   const Real & equiv_stress,
                                                   const RankTwoTensor & dev_stress,
                                                   const RankTwoTensor & stress,
                                                   const unsigned int i)
{
  // Set model parameters for evaluations
  _current_n = _powers[i];

  // Run non-linear solve for gauge stress
  returnMappingSolve(equiv_stress, gauge_stress, _console);

  if (gauge_stress < 0.0)
    mooseError("In ",
               _name,
               ": Gauge stress is less than zero. Something is wrong with the inner Newton solve");

  dpsi_dgauge = (*_coefficients[i])[_qp] * std::pow(gauge_stress, _current_n);

  strain_rate +=
      dpsi_dgauge * computeDGaugeDSigma(gauge_stress, equiv_stress, dev_stress, stress, _current_n);
}
