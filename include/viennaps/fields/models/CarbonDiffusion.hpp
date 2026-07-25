#pragma once

/// CarbonDiffusion — C traps interstitials (C + I → CI), suppressing TED.

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class CarbonDiffusion : public DiffusionModel<NumericType> {
public:
  CarbonDiffusion() { this->setName("CarbonDiffusion"); }

  void setTrapRate(NumericType kf) { kf_ = kf; }

  void applyTrapStep(std::vector<NumericType> &C, std::vector<NumericType> &I,
                     std::vector<NumericType> &CI, NumericType dt) const {
    const std::size_t n = std::min({C.size(), I.size(), CI.size()});
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType r = kf_ * C[i] * I[i] * dt;
      C[i] = std::max(NumericType(0), C[i] - r);
      I[i] = std::max(NumericType(0), I[i] - r);
      CI[i] += r;
    }
  }

  int numSpecies() const override { return 3; }
  std::vector<std::string> speciesNames() const override {
    return {"Carbon", "Interstitial", "CarbonInterstitial"};
  }

private:
  NumericType kf_ = NumericType(1e-18);
};

} // namespace viennaps
