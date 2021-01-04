//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVMatAdvection.h"

class INSFVVelocityVariable;
class INSFVPressureVariable;

/**
 * An advection kernel that implements interpolation schemes specific to Navier-Stokes flow
 * physics
 */
class INSFVMomentumAdvection : public FVMatAdvection
{
public:
  static InputParameters validParams();
  INSFVMomentumAdvection(const InputParameters & params);

protected:
  /**
   * interpolation overload for the velocity
   */
  void interpolate(Moose::FV::InterpMethod m,
                   ADRealVectorValue & interp_v,
                   const ADRealVectorValue & elem_v,
                   const ADRealVectorValue & neighbor_v);

  ADReal computeQpResidual() override;

  void residualSetup() override final { clearRCCoeffs(); }
  void jacobianSetup() override final { clearRCCoeffs(); }

  /// The dynamic viscosity on the FaceInfo elem
  const ADMaterialProperty<Real> & _mu_elem;

  /// The dynamic viscosity on the FaceInfo neighbor
  const ADMaterialProperty<Real> & _mu_neighbor;

  /**
   * Returns the Rhie-Chow 'a' coefficient for the requested element \p elem
   * @param elem The elem to get the Rhie-Chow coefficient for
   * @param mu The dynamic viscosity
   */
  const VectorValue<ADReal> & rcCoeff(const Elem & elem, const ADReal & mu) const;

  /**
   * method for computing the Rhie-Chow 'a' coefficients for the given elem \p elem
   * @param elem The elem to compute the Rhie-Chow coefficient for
   * @param mu The dynamic viscosity
   */
  VectorValue<ADReal> coeffCalculator(const Elem & elem, const ADReal & mu) const;

  /**
   * Clear the RC 'a' coefficient cache
   */
  void clearRCCoeffs();

  bool skipForBoundary(const FaceInfo & fi) const override;

  /// pressure variable
  const INSFVPressureVariable * const _p_var;
  /// x-velocity
  const INSFVVelocityVariable * const _u_var;
  /// y-velocity
  const INSFVVelocityVariable * const _v_var;
  /// z-velocity
  const INSFVVelocityVariable * const _w_var;

  /// Density
  const Real & _rho;

  /// the dimension of the simulation
  const unsigned int _dim;

  /// The interpolation method to use for the velocity
  Moose::FV::InterpMethod _velocity_interp_method;

  /// Boundary IDs with no slip walls
  std::set<BoundaryID> _no_slip_wall_boundaries;

  /// Boundary IDs with slip walls
  std::set<BoundaryID> _slip_wall_boundaries;

  /// Flow Boundary IDs
  std::set<BoundaryID> _flow_boundaries;

private:
  /// A map from elements to the 'a' coefficients used in the Rhie-Chow interpolation. The size of
  /// the vector is equal to the number of threads in the simulation. We maintain a map from
  /// MooseApp pointer to RC coefficients in order to support MultiApp simulations
  static std::unordered_map<const MooseApp *,
                            std::vector<std::unordered_map<const Elem *, VectorValue<ADReal>>>>
      _rc_a_coeffs;
};
