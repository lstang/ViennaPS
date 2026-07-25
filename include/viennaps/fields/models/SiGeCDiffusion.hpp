#pragma once

/// SiGeCDiffusion — carbon I-trapping in SiGe context (extends Phase 4 Carbon).

#include "CarbonDiffusion.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class SiGeCDiffusion : public CarbonDiffusion<NumericType> {
public:
  SiGeCDiffusion() { this->setName("SiGeCDiffusion"); }

  void setGeFraction(NumericType x) { x_Ge_ = x; }
  NumericType geFraction() const { return x_Ge_; }

  /// TED enhancement factor: free I reduced by carbon trapping.
  static NumericType tedFactor(NumericType C_I, NumericType C_I_eq,
                               NumericType C_C, NumericType trapStrength) {
    const NumericType freeI =
        C_I / (NumericType(1) + trapStrength * std::max(C_C, NumericType(0)));
    return freeI / std::max(C_I_eq, NumericType(1));
  }

private:
  NumericType x_Ge_ = NumericType(0.2);
};

} // namespace viennaps
