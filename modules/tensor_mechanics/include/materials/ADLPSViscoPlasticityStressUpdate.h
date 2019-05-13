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
  using ADLPSViscoPlasticityStressUpdate<compute_stage>::_porosity;                                \
  using ADLPSViscoPlasticityStressUpdate<compute_stage>::_porosity_old

template <ComputeStage>
class ADLPSViscoPlasticityStressUpdate;

declareADValidParams(ADLPSViscoPlasticityStressUpdate);

template <ComputeStage compute_stage>
class ADLPSViscoPlasticityStressUpdate : public ADStressUpdateBase<compute_stage>,
                                         public ADSingleVariableReturnMappingSolution<compute_stage>
{
public:
  ADLPSViscoPlasticityStressUpdate(const InputParameters & parameters);

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

  virtual void propagateQpStatefulProperties() override;

  /**
   * Perform any necessary initialization before return mapping iterations
   * @param effective_trial_stress Effective trial stress
   * @param elasticityTensor     Elasticity tensor
   */
  virtual void computeStressInitialize(const ADReal & /*effective_trial_stress*/,
                                       const ADRankFourTensor & /*elasticity_tensor*/)
  {
  }

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
  virtual ADReal initialGuess(const ADReal & effective_trial_stress) override;

  /**
   * Perform any necessary steps to finalize state after return mapping iterations
   * @param inelasticStrainIncrement Inelastic strain increment
   */
  virtual void computeStressFinalize(const ADRankTwoTensor & /*plastic_strain_increment*/) {}

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & scalar) override;
  virtual ADReal computeDerivative(const ADReal & /*effective_trial_stress*/,
                                   const ADReal & /*scalar*/) override
  {
    return _derivative;
  }

  void outputIterationSummary(std::stringstream * iter_output,
                              const unsigned int total_it) override;

  virtual ADReal computeSwellingIncrement(const ADReal & /*hydrostatic_stress*/) { return 0.0; }

  ADReal computeM(const ADReal & hydro_stress, const ADReal & gauge_stress);

  ADReal computeH(const Real n, const ADReal & M, const bool derivative = false);

  ADRankTwoTensor computeDGaugeDSigma(const ADReal & gauge_stress,
                                      const ADReal & equiv_stress,
                                      const ADRankTwoTensor & dev_stress,
                                      const ADRankTwoTensor & stress,
                                      const Real n);

  void computeCreepStrainIncrementN(ADReal & gauge_stress,
                          ADReal & dpsi_dgauge,
                          ADRankTwoTensor & creep_strain_increment,
                          const ADReal & equiv_stress,
                          const ADRankTwoTensor & dev_stress,
                          const ADRankTwoTensor & stress,
                          const unsigned int i);

  ///@{ Effective inelastic strain material property
  ADMaterialProperty(Real) & _effective_inelastic_strain;
  const MaterialProperty<Real> & _effective_inelastic_strain_old;
  ///@}

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

  /// Exponent on the effective stress
  const std::vector<Real> _powers;

  /// Number of total models
  const unsigned int _num_models;

  /// Leading coefficient
  std::vector<const ADMaterialProperty(Real) *> _coefficients;

  /// Gauge stress
  std::vector<ADMaterialProperty(Real) *> _gauge_stresses;

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

  /// Containers for current model parameters
  Real _current_n;

  usingStressUpdateBaseMembers;
  usingSingleVariableReturnMappingSolutionMembers;
};
