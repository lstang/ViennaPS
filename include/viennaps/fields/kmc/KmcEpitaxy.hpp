#pragma once

/// KMC lattice epitaxy helpers (Phase 8 skeleton).

#include "KmcLattice.hpp"

#include <algorithm>
#include <cmath>

namespace viennaps {

struct KmcSurfaceEvent {
  enum class Kind { Adsorb, Desorb, Diffuse, Twin } kind = Kind::Adsorb;
  double rate = 0.0;
  int i = 0, j = 0, k = 0;
};

class KmcEpitaxyModel {
public:
  void setGrowthRate(double r) { growthRate_ = r; }
  void setGeFraction(double x) { xGe_ = std::clamp(x, 0.0, 1.0); }

  /// Planar growth: deposit on topmost empty site in each column.
  int planarGrow(KmcLattice &lat, int speciesCode, int layers) {
    int deposited = 0;
    for (int layer = 0; layer < layers; ++layer) {
      for (int j = 0; j < lat.ny(); ++j)
        for (int i = 0; i < lat.nx(); ++i) {
          for (int k = 0; k < lat.nz(); ++k) {
            auto &s = lat.at(i, j, k);
            if (!s.occupied) {
              s.occupied = true;
              s.species = speciesCode;
              ++deposited;
              break;
            }
          }
        }
    }
    return deposited;
  }

  /// Coordination-based attachment probability ~ z/zmax.
  double attachmentProbability(int coordination, int zmax = 4) const {
    return std::min(1.0, static_cast<double>(coordination) /
                             static_cast<double>(std::max(zmax, 1)));
  }

  /// Ge mole-fraction dependent growth rate factor.
  double geGrowthFactor() const {
    // SiGe growth often slows mildly with Ge fraction (simple linear model).
    return 1.0 - 0.3 * xGe_;
  }

  double effectiveGrowthRate() const {
    return growthRate_ * geGrowthFactor();
  }

private:
  double growthRate_ = 1.0;
  double xGe_ = 0.0;
};

class KmcVisibility {
public:
  /// Simple z-buffer shadowing: site visible if no occupied site above.
  static bool isVisible(const KmcLattice &lat, int i, int j, int k) {
    for (int kk = k + 1; kk < lat.nz(); ++kk)
      if (lat.at(i, j, kk).occupied)
        return false;
    return true;
  }
};

} // namespace viennaps
