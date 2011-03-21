#ifndef SIDEINTEGRAL_H_
#define SIDEINTEGRAL_H_

#include "SidePostprocessor.h"

//Forward Declarations
class SideIntegral;

template<>
InputParameters validParams<SideIntegral>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class SideIntegral : public SidePostprocessor
{
public:
  SideIntegral(const std::string & name, InputParameters parameters);
  
  virtual void initialize();
  virtual void execute();
  virtual Real getValue();

protected:
  virtual Real computeQpIntegral();

  Real _integral_value;
};
 
#endif
