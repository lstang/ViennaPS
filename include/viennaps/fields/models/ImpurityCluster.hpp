#pragma once

/// ImpurityCluster (BIC) — B + I ⇌ BIC
/// dC_BIC/dt = k_f * C_B * C_I - k_r * C_BIC

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ImpurityCluster : public DiffusionModel<NumericType> {
public:
  ImpurityCluster(std::string bic = "BIC", std::string boron = "Boron",
                  std::string interstitial = "Interstitial")
      : bic_(std::move(bic)), B_(std::move(boron)),
        I_(std::move(interstitial)) {
    this->setName("ImpurityCluster");
  }

  void setRates(NumericType kf, NumericType kr) {
    kf_ = kf;
    kr_ = kr;
  }

  void applyReactionStep(std::vector<NumericType> &B,
                         std::vector<NumericType> &I,
                         std::vector<NumericType> &BIC,
                         NumericType dt) const {
    const std::size_t n = std::min({B.size(), I.size(), BIC.size()});
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType form = kf_ * B[i] * I[i];
      const NumericType diss = kr_ * BIC[i];
      const NumericType dC = (form - diss) * dt;
      BIC[i] = std::max(NumericType(0), BIC[i] + dC);
      B[i] = std::max(NumericType(0), B[i] - dC);
      I[i] = std::max(NumericType(0), I[i] - dC);
    }
  }

  int numSpecies() const override { return 3; }
  std::vector<std::string> speciesNames() const override {
    return {bic_, B_, I_};
  }

private:
  std::string bic_, B_, I_;
  NumericType kf_ = NumericType(1e-18);
  NumericType kr_ = NumericType(1e-3);
};

} // namespace viennaps
