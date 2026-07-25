#pragma once

/// ChargedPairDiffusion — pair diffusion with Fermi-level dependent D_pair.

#include "PairDiffusion.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
#include <string>

namespace viennaps {

template <class NumericType>
class ChargedPairDiffusion : public PairDiffusion<NumericType> {
public:
  ChargedPairDiffusion(std::string dopant = "Boron",
                       std::string interstitial = "Interstitial")
      : PairDiffusion<NumericType>(std::move(dopant),
                                   std::move(interstitial)) {
    this->setName("ChargedPairDiffusion");
  }

  void setFermiEnhancement(NumericType alpha, NumericType ni) {
    alpha_ = alpha;
    ni_ = ni;
  }

  NumericType getDiffusivity(NumericType C_I, NumericType C_dopant,
                             NumericType T) const {
    NumericType ni = ni_;
    if (ni <= NumericType(0)) {
      ni = IntrinsicCarrier<NumericType>{}.ni(T, "Si");
    }
    const NumericType n = std::max(C_dopant, ni);
    const NumericType fermiFactor =
        NumericType(1) + alpha_ * n / std::max(ni, NumericType(1));
    return PairDiffusion<NumericType>::getDiffusivity(C_I, T) * fermiFactor;
  }

private:
  NumericType alpha_ = NumericType(1);
  NumericType ni_ = NumericType(0);
};

} // namespace viennaps
