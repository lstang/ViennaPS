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
  void setSurfaceSegregation(double s) { surfSeg_ = std::clamp(s, 0.0, 1.0); }

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

  /// Coordination-based growth: only deposit if neighbor count ≥ minCoord.
  int coordinationGrow(KmcLattice &lat, int speciesCode, int minCoord = 2) {
    int deposited = 0;
    for (int k = 0; k < lat.nz(); ++k)
      for (int j = 0; j < lat.ny(); ++j)
        for (int i = 0; i < lat.nx(); ++i) {
          auto &s = lat.at(i, j, k);
          if (s.occupied)
            continue;
          int z = 0;
          const int nbr[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                                 {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
          for (auto &d : nbr) {
            const int i1 = i + d[0], j1 = j + d[1], k1 = k + d[2];
            if (i1 < 0 || j1 < 0 || k1 < 0 || i1 >= lat.nx() || j1 >= lat.ny() ||
                k1 >= lat.nz())
              continue;
            if (lat.at(i1, j1, k1).occupied)
              ++z;
          }
          if (z >= minCoord &&
              attachmentProbability(z) * geGrowthFactor() > 0.1) {
            s.occupied = true;
            s.species = speciesCode;
            ++deposited;
          }
        }
    return deposited;
  }

  /// Surface segregation: swap subsurface Ge-rich (species 5) with surface Si.
  int surfaceSegregate(KmcLattice &lat, int geCode = 5, int siCode = 2) {
    int swaps = 0;
    if (surfSeg_ <= 0.0)
      return 0;
    for (int j = 0; j < lat.ny(); ++j)
      for (int i = 0; i < lat.nx(); ++i) {
        int kTop = -1;
        for (int k = lat.nz() - 1; k >= 0; --k) {
          if (lat.at(i, j, k).occupied) {
            kTop = k;
            break;
          }
        }
        if (kTop <= 0)
          continue;
        auto &surf = lat.at(i, j, kTop);
        auto &sub = lat.at(i, j, kTop - 1);
        if (sub.species == geCode && surf.species == siCode) {
          std::swap(surf.species, sub.species);
          ++swaps;
        }
      }
    return swaps;
  }

  /// Twin defect: mark a surface site with twin code (6) at low probability.
  int formTwin(KmcLattice &lat, int twinCode = 6) {
    int n = 0;
    for (int j = 0; j < lat.ny(); ++j)
      for (int i = 0; i < lat.nx(); ++i) {
        for (int k = lat.nz() - 1; k >= 0; --k) {
          auto &s = lat.at(i, j, k);
          if (s.occupied) {
            if ((i + j + k) % 7 == 0) {
              s.species = twinCode;
              ++n;
            }
            break;
          }
        }
      }
    return n;
  }

  /// Coordination-based attachment probability ~ z/zmax.
  double attachmentProbability(int coordination, int zmax = 4) const {
    return std::min(1.0, static_cast<double>(coordination) /
                             static_cast<double>(std::max(zmax, 1)));
  }

  /// Ge mole-fraction dependent growth rate factor.
  double geGrowthFactor() const {
    return 1.0 - 0.3 * xGe_;
  }

  double effectiveGrowthRate() const {
    return growthRate_ * geGrowthFactor();
  }

private:
  double growthRate_ = 1.0;
  double xGe_ = 0.0;
  double surfSeg_ = 0.0;
};

} // namespace viennaps
