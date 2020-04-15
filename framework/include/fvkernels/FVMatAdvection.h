#pragma once

#include "FVFluxKernel.h"

class FVMatAdvection : public FVFluxKernel
{
public:
  static InputParameters validParams();
  FVMatAdvection(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<RealVectorValue> & _vel_left;
  const ADMaterialProperty<RealVectorValue> & _vel_right;
};
