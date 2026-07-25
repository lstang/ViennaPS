#pragma once

#include "ConstantDiffusion.hpp"

#include <string>

namespace viennaps {

/// NitrogenDiffusion — N transport (ConstantDiffusion specialization).
template <class NumericType>
class NitrogenDiffusion : public ConstantDiffusion<NumericType> {
public:
  NitrogenDiffusion() : ConstantDiffusion<NumericType>("Nitrogen") {
    this->setName("NitrogenDiffusion");
    this->setDiffusivity(NumericType(0.05), NumericType(3.0));
  }
};

} // namespace viennaps
