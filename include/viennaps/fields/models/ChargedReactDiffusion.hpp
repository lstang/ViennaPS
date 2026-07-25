#pragma once

/// ChargedReactDiffusion — I+V recombination with charge-enhanced rate.
/// k = k0 * (1 + gamma * n/ni)

#include "ReactDiffusion.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ChargedReactDiffusion : public ReactDiffusion<NumericType> {
public:
  ChargedReactDiffusion(std::string interstitial = "Interstitial",
                        std::string vacancy = "Vacancy")
      : ReactDiffusion<NumericType>(std::move(interstitial),
                                    std::move(vacancy)) {
    this->setName("ChargedReactDiffusion");
  }

  void setChargeEnhancement(NumericType gamma, NumericType ni) {
    gamma_ = gamma;
    ni_ = ni;
  }

  void setDopantConcentration(NumericType C) { C_dopant_ = C; }

  NumericType effectiveRate(NumericType T) const {
    NumericType ni = ni_;
    if (ni <= NumericType(0)) {
      ni = IntrinsicCarrier<NumericType>{}.ni(T, "Si");
    }
    const NumericType n = std::max(C_dopant_, ni);
    return this->recombinationRate() *
           (NumericType(1) + gamma_ * n / std::max(ni, NumericType(1)));
  }

  /// Host reaction step using charge-enhanced k.
  void applyChargedReactionStep(std::vector<NumericType> &I,
                                std::vector<NumericType> &V, NumericType dt,
                                NumericType T) const {
    const NumericType kEff = effectiveRate(T);
    const std::size_t n = std::min(I.size(), V.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType r = kEff * I[i] * V[i] * dt;
      I[i] = std::max(NumericType(0), I[i] - r);
      V[i] = std::max(NumericType(0), V[i] - r);
    }
  }

private:
  NumericType gamma_ = NumericType(1);
  NumericType ni_ = NumericType(0);
  NumericType C_dopant_ = NumericType(0);
};

} // namespace viennaps
