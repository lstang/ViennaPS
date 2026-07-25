#pragma once

/// DislocationLoop — loop growth from interstitial supersaturation.
/// dC_loop/dt = k * max(C_I/C_I_eq - 1, 0)^p

#include "../DiffusionModel.hpp"
#include "../PointDefectEquilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class DislocationLoop : public DiffusionModel<NumericType> {
public:
  explicit DislocationLoop(std::string loop = "DislocationLoop",
                           std::string interstitial = "Interstitial")
      : loop_(std::move(loop)), I_(std::move(interstitial)) {
    this->setName("DislocationLoop");
  }

  void setGrowth(NumericType k, NumericType p = NumericType(1)) {
    k_ = k;
    p_ = p;
  }

  void setCIEq(NumericType Ceq) { Ceq_ = Ceq; }

  void applyReactionStep(std::vector<NumericType> &I,
                         std::vector<NumericType> &loop,
                         NumericType dt) const {
    const std::size_t n = std::min(I.size(), loop.size());
    const NumericType Ceq =
        (Ceq_ > NumericType(0))
            ? Ceq_
            : PointDefectEquilibrium<NumericType>{}.C_I_eq(this->T_, "Si");
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType super =
          std::max(NumericType(0), I[i] / std::max(Ceq, NumericType(1)) -
                                       NumericType(1));
      const NumericType dL = k_ * std::pow(super, p_) * dt;
      loop[i] += dL;
      // Consume interstitials proportionally (1:1 atomistic mapping).
      I[i] = std::max(NumericType(0), I[i] - dL);
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {loop_, I_};
  }

private:
  std::string loop_, I_;
  NumericType k_ = NumericType(1e12);
  NumericType p_ = NumericType(1);
  NumericType Ceq_ = NumericType(0); // 0 => lookup
};

} // namespace viennaps
