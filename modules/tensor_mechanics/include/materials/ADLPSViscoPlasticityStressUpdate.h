//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADStressUpdateBase.h"
#include "ADSingleVariableReturnMappingSolution.h"

#define usingLPSViscoPlasticityStressUpdateMembers                                                 \
  usingStressUpdateBaseMembers;                                                                    \
  usingSingleVariableReturnMappingSolutionMembers;                                                 \
  using ADLPSViscoPlasticityStressUpdate<compute_stage>::_three_shear_modulus;                     \
  using ADLPSViscoPlasticityStressUpdate<compute_stage>::propagateQpStatefulPropertiesRadialReturn

// Forward declarations
template <ComputeStage>
class ADLPSViscoPlasticityStressUpdate;

declareADValidParams(ADLPSViscoPlasticityStressUpdate);

/**
 * ADLPSViscoPlasticityStressUpdate computes the radial return stress increment for
 * an isotropic elastic-viscoplasticity model after interating on the difference
 * between new and old trial stress increments.  This radial return mapping class
 * acts as a base class for the radial return creep and plasticity classes / combinations.
 * The stress increment computed by ADLPSViscoPlasticityStressUpdate is used by
 * ComputeMultipleInelasticStress which computes the elastic stress for finite
 * strains.  This return mapping class is acceptable for finite strains but not
 * total strains.
 * This class is based on the Elasto-viscoplasticity algorithm in F. Dunne and N.
 * Petrinic's Introduction to Computational Plasticity (2004) Oxford University Press.
 */
template <ComputeStage compute_stage>
class ADLPSViscoPlasticityStressUpdate : public ADStressUpdateBase<compute_stage>,
                                         public ADSingleVariableReturnMappingSolution<compute_stage>
{
public:
  ADLPSViscoPlasticityStressUpdate(const InputParameters & parameters);

  /**
   * A radial return (J2) mapping method is performed with return mapping
   * iterations.
   * @param strain_increment              Sum of elastic and inelastic strain increments
   * @param inelastic_strain_increment    Inelastic strain increment calculated by this class
   * @param rotation increment            Not used by this class
   * @param stress                    New trial stress from pure elastic calculation
   * @param stress_old                    Old state of stress
   * @param elasticity_tensor             Rank 4 C_{ijkl}, must be isotropic
   * @param elastic_strain_old            Old state of total elastic strain
   */
  virtual void updateState(ADRankTwoTensor & strain_increment,
                           ADRankTwoTensor & inelastic_strain_increment,
                           const ADRankTwoTensor & rotation_increment,
                           ADRankTwoTensor & stress,
                           const RankTwoTensor & stress_old,
                           const ADRankFourTensor & elasticity_tensor,
                           const RankTwoTensor & elastic_strain_old) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & scalar_effective_inelastic_strain) override;

  virtual Real maximumPermissibleValue(const ADReal & effective_trial_stress) const override;
  virtual Real minimumPermissibleValue(const ADReal & effective_trial_stress) const override;

  /**
   * Compute the limiting value of the time step for this material
   * @return Limiting time step
   */
  virtual Real computeTimeStepLimit() override;

  /**
   * Does the model require the elasticity tensor to be isotropic?
   */
  bool requiresIsotropicTensor() override { return true; }

protected:
  virtual void initQpStatefulProperties() override;

  /**
   * Propagate the properties pertaining to this intermediate class.  This
   * is intended to be called from propagateQpStatefulProperties() in
   * classes that inherit from this one.
   * This is intentionally named uniquely because almost all models that derive
   * from this class have their own stateful properties, and this forces them
   * to define their own implementations of propagateQpStatefulProperties().
   */
  void propagateQpStatefulPropertiesRadialReturn();

  /**
   * Calculate the derivative of the strain increment with respect to the updated stress.
   * @param effective_trial_stress Effective trial stress
   * @param scalar                 Inelastic strain increment magnitude being solved for
   */
  virtual Real computeStressDerivative(const ADReal & /*effective_trial_stress*/,
                                       const ADReal & /*scalar*/)
  {
    return 0.0;
  }

  /**
   * Compute an initial guess for the value of the scalar. For some cases, an
   * intellegent starting point can provide enhanced robustness in the Newton
   * iterations. This is also an opportunity for classes that derive from this
   * to perform initialization tasks.
   * @param effective_trial_stress Effective trial stress
   */
  virtual ADReal initialGuess(const ADReal & effective_trial_stress) override
  {
    return effective_trial_stress;
  }

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & scalar) override;
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & scalar) override;

  void outputIterationSummary(std::stringstream * iter_output,
                              const unsigned int total_it) override;

  ADReal computeM(const ADReal & hydro_stress,
                  const ADReal & gauge_stress,
                  const unsigned int derivative = 0);

  ADReal computeH(const Real n, const ADReal & M, const bool derivative = false);

  ADReal computePowerLawConstant();

  ADRankTwoTensor computeDGaugeDSigma(const ADReal & gauge_stress,
                                      const ADReal & equiv_stress,
                                      const ADRankTwoTensor & dev_stress,
                                      const ADRankTwoTensor & stress,
                                      const Real n);

  void computeNStrainRate(ADReal & gauge_stress,
                          ADReal & dpsi_dgauge,
                          ADRankTwoTensor & strain_rate,
                          const ADReal & equiv_stress,
                          const ADRankTwoTensor & dev_stress,
                          const ADRankTwoTensor & stress,
                          const Real n);

  ///@{ Effective inelastic strain material property
  ADMaterialProperty(Real) & _effective_inelastic_strain;
  const MaterialProperty<Real> & _effective_inelastic_strain_old;
  ///@}

  /// Gauge stress
  ADMaterialProperty(Real) & _gauge_stress;

  ///@{ Creep strain material property
  ADMaterialProperty(RankTwoTensor) & _creep_strain;
  const MaterialProperty<RankTwoTensor> & _creep_strain_old;
  ///@}

  ///@{ Porosity material property
  ADMaterialProperty(Real) & _porosity;
  const MaterialProperty<Real> & _porosity_old;
  ///@}

  /// Max increment for inelastic strain
  Real _max_inelastic_increment;

  /// Temperature variable value
  const ADVariableValue & _temperature;

  /// Leading coefficient
  const Real _coefficient;

  /// Exponent on the effective stress
  const Real _n;

  /// Activation energy for exp term
  const Real _activation_energy;

  /// Gas constant for exp term
  const Real _gas_constant;

  /// Initial porosity to setup stateful materials
  const Real _initial_porosity;

  /// Flag to enable verbose output
  const bool _verbose;

  /// Container for hydrostatic stress
  ADReal _hydro_stress;

  /// Rank two identity tensor
  const RankTwoTensor _identity_two;

  /// Derivative of hydrostatic stress with respect to the stress tensor
  const RankTwoTensor _dhydro_stress_dsigma;

  /// Container for _derivative
  ADReal _derivative;

  /// Container for current value of n
  Real _current_n;

  usingStressUpdateBaseMembers;
  usingSingleVariableReturnMappingSolutionMembers;
};
