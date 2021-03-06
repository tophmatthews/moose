//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVMatAdvectionOutflowBC.h"

registerADMooseObject("MooseApp", FVMatAdvectionOutflowBC);

InputParameters
FVMatAdvectionOutflowBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params.addRequiredParam<MaterialPropertyName>("vel", "advection velocity");
  params.addParam<MaterialPropertyName>(
      "advected_quantity",
      "An optional parameter for specifying an advected quantity from a material property. If this "
      "is not specified, then the advected quantity will simply be the variable that this object "
      "is acting on");
  return params;
}

FVMatAdvectionOutflowBC::FVMatAdvectionOutflowBC(const InputParameters & params)
  : FVFluxBC(params), _vel(getADMaterialProperty<RealVectorValue>("vel"))
{
  if (isParamValid("advected_quantity"))
    _adv_quant = &getADMaterialProperty<Real>("advected_quantity").get();
  else
    _adv_quant = &_u;
}

ADReal
FVMatAdvectionOutflowBC::computeQpResidual()
{
  return _normal * _vel[_qp] * (*_adv_quant)[_qp];
}
