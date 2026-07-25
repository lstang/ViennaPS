#pragma once

/// GrainModel — average grain radius evolution for polysilicon.
/// dR/dt = k * exp(-Ea/(kB*T)) / R^n

#include <algorithm>
#include <cmath>

namespace viennaps {

template <class NumericType>
class GrainModel {
public:
  void setGrowthParameters(NumericType k, NumericType Ea_eV,
                           NumericType n = NumericType(1)) {
    k_ = k;
    Ea_ = Ea_eV;
    n_ = n;
  }

  void setRadius(NumericType R) { R_ = std::max(R, NumericType(1e-12)); }
  NumericType radius() const { return R_; }

  /// Advance grain growth over dt at temperature T [K].
  void advance(NumericType T, NumericType dt) {
    if (T <= NumericType(0) || dt <= NumericType(0))
      return;
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    const NumericType rate =
        k_ * std::exp(-Ea_ / (kB * T)) /
        std::pow(std::max(R_, NumericType(1e-12)), n_);
    R_ += rate * dt;
  }

  /// Grain-boundary volume fraction (simple 2D estimate): f_gb ~ δ / R
  NumericType grainBoundaryFraction(NumericType boundaryWidth) const {
    return std::min(NumericType(1),
                    boundaryWidth / std::max(R_, NumericType(1e-12)));
  }

private:
  NumericType R_ = NumericType(1e-5); // cm, ~0.1 um default
  NumericType k_ = NumericType(1e-8);
  NumericType Ea_ = NumericType(2.5);
  NumericType n_ = NumericType(1);
};

} // namespace viennaps
