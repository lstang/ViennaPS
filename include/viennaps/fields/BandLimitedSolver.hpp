#pragma once

/// BandLimitedSolver - Restrict field solves to an active band (3D scaling).
/// Marks profile bins / mesh elements inside [yMin,yMax] and optionally
/// solves diffusion only on that band (production performance path).

#include "PhysicsField.hpp"
#include "AmgclSolver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace viennaps {

template <class NumericType>
class BandLimitedSolver {
public:
  void setBand(NumericType yMin, NumericType yMax) {
    yMin_ = std::min(yMin, yMax);
    yMax_ = std::max(yMin, yMax);
  }

  void setDimension(int d) { dim_ = d; }
  int getDimension() const { return dim_; }

  /// Initialize 2D or 3D mesh on the field (MFEM when available).
  void initMesh(PhysicsField<NumericType>& field, int n = 16) {
    if (dim_ >= 3) {
      field.initMesh3D(0, 1, 0, 1, 0, 1, n, n, std::max(4, n / 2));
    } else {
      field.initMeshFromBounds(0, 1, 0, 1, n, n);
    }
    std::cout << "[BandLimitedSolver] init mesh dim=" << dim_ << " n=" << n << "\n";
  }

  /// Count profile bins inside the active band (normalized depth 0..1 as y).
  int countActiveBins(const PhysicsField<NumericType>& field) const {
    auto n = field.getProfileSize();
    int active = 0;
    for (std::size_t i = 0; i < n; ++i) {
      NumericType y = (static_cast<NumericType>(i) + NumericType(0.5)) /
                      static_cast<NumericType>(n);
      if (y >= yMin_ && y <= yMax_) active++;
    }
    return active;
  }

  /// Diffuse only inside the band using amgcl/CG; outside bins frozen.
  AmgclSolveResult diffuseSpeciesBand(PhysicsField<NumericType>& field,
                                      const std::string& species,
                                      NumericType Ddt) {
    auto* prof = field.getProfileMutable(species);
    if (!prof || prof->empty()) return {};

    // Extract band
    std::vector<std::size_t> idx;
    std::vector<NumericType> band;
    for (std::size_t i = 0; i < prof->size(); ++i) {
      NumericType y = (static_cast<NumericType>(i) + NumericType(0.5)) /
                      static_cast<NumericType>(prof->size());
      if (y >= yMin_ && y <= yMax_) {
        idx.push_back(i);
        band.push_back((*prof)[i]);
      }
    }
    if (band.size() < 2) {
      std::cout << "[BandLimitedSolver] band too small\n";
      return {};
    }
    auto res = solveProfileDiffusion(band, Ddt);
    if (res.ok) {
      for (std::size_t k = 0; k < idx.size(); ++k)
        (*prof)[idx[k]] = band[k];
      field.refreshDose(species);
    }
    lastActive_ = static_cast<int>(band.size());
    lastTotal_ = static_cast<int>(prof->size());
    std::cout << "[BandLimitedSolver] diffuse " << species << " active="
              << lastActive_ << "/" << lastTotal_ << " dim=" << dim_ << "\n";
    return res;
  }

  int getLastActive() const { return lastActive_; }
  int getLastTotal() const { return lastTotal_; }

private:
  NumericType yMin_ = NumericType(0.2);
  NumericType yMax_ = NumericType(0.7);
  int dim_ = 2;
  int lastActive_ = 0;
  int lastTotal_ = 0;
};

} // namespace viennaps
