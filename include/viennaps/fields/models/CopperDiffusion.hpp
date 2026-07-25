#pragma once

/// CopperDiffusion — fast Cu diffuser with ion-pairing enhancement in
/// doped regions: D = D0 * (1 + beta * C_dopant / ni)

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class CopperDiffusion : public DiffusionModel<NumericType> {
public:
  CopperDiffusion() { this->setName("CopperDiffusion"); }

  void setD0(NumericType D0) { D0_ = D0; }
  void setIonPairing(NumericType beta) { beta_ = beta; }

  NumericType getDiffusivity(NumericType C_dopant, NumericType T) const {
    IntrinsicCarrier<NumericType> ic;
    const NumericType ni = ic.ni(T, "Si");
    return D0_ * (NumericType(1) +
                  beta_ * C_dopant / std::max(ni, NumericType(1)));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {"Copper"};
  }

private:
  NumericType D0_ = NumericType(1e-5); // Cu is fast
  NumericType beta_ = NumericType(1);
};

} // namespace viennaps
