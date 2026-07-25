#pragma once

/// Cluster311 — {311} interstitial cluster growth / dissociation.
///
/// dC_311/dt = k_f * C_I^n - k_r * C_311
/// dC_I/dt  -= n * (k_f * C_I^n - k_r * C_311)   (atom conservation)

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class Cluster311 : public DiffusionModel<NumericType> {
public:
  Cluster311(std::string cluster = "311", std::string interstitial = "Interstitial")
      : cluster_(std::move(cluster)), I_(std::move(interstitial)) {
    this->setName("Cluster311");
  }

  void setRates(NumericType kf, NumericType kr, int n = 2) {
    kf_ = kf;
    kr_ = kr;
    n_ = n;
  }

  void applyReactionStep(std::vector<NumericType> &I,
                         std::vector<NumericType> &C311,
                         NumericType dt) const {
    const std::size_t n = std::min(I.size(), C311.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType growth =
          kf_ * std::pow(std::max(I[i], NumericType(0)), n_);
      const NumericType dissoc = kr_ * C311[i];
      const NumericType dC = (growth - dissoc) * dt;
      C311[i] = std::max(NumericType(0), C311[i] + dC);
      I[i] = std::max(NumericType(0),
                      I[i] - static_cast<NumericType>(n_) * dC);
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {cluster_, I_};
  }

  NumericType kf() const { return kf_; }
  NumericType kr() const { return kr_; }

private:
  std::string cluster_, I_;
  NumericType kf_ = NumericType(1e-20);
  NumericType kr_ = NumericType(1e-3);
  int n_ = 2;
};

} // namespace viennaps
