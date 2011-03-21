#ifndef POSTPROCESSORDATA_H_
#define POSTPROCESSORDATA_H_

#include <map>
//MOOSE includes
#include "Moose.h"


class PostprocessorData
{
public:
  PostprocessorData();

  void init(const std::string & name);

  bool hasPostprocessor(const std::string & name);

  PostprocessorValue & getPostprocessorValue(const std::string & name);

  void storeValue(const std::string & name, PostprocessorValue value);
  
protected:
  std::map<std::string, PostprocessorValue> _values;
};

#endif //POSTPROCESSORDATA_H_
