#pragma once

/// GeometryFieldCoupler - Level-set / material-map ↔ PhysicsField mesh coupling.
///
/// Provides region marking without requiring a live ViennaLS Domain:
///   material ids: 1=Si, 2=SiO2, 3=Mask, 0=Ambient
///
/// Full Domain integration: call markFromDepthInterfaces with oxide/mask depths
/// extracted from level-set zero contours.

#include "PhysicsField.hpp"

#include <functional>
#include <iostream>
#include <memory>
#include <cmath>
#include <utility>
#include <vector>

namespace viennaps {

template <class NumericType>
class GeometryFieldCoupler {
public:
  GeometryFieldCoupler() = default;

  explicit GeometryFieldCoupler(std::shared_ptr<PhysicsField<NumericType>> field)
      : field_(std::move(field)) {}

  void setPhysicsField(std::shared_ptr<PhysicsField<NumericType>> f) { field_ = std::move(f); }

  /// Mark a simple LOCOS-like stack:
  ///   y > surfaceY            -> Ambient (0)
  ///   oxideBottom < y <= surface -> Oxide (2) if |x| > maskHalfWidth else Mask (3) near surface
  ///   y <= oxideBottom        -> Si (1)
  void markLocosLike(NumericType surfaceY = NumericType(0.5),
                     NumericType oxideBottom = NumericType(0.3),
                     NumericType maskHalfWidth = NumericType(0.25),
                     NumericType domainHalfWidth = NumericType(0.5)) {
    if (!field_) return;

    field_->initMeshFromBounds(-static_cast<double>(domainHalfWidth),
                               static_cast<double>(domainHalfWidth), 0.0, 1.0, 32, 32);

    field_->markMeshRegions(
        [=](double x, double y, double /*z*/) -> int {
          if (y > static_cast<double>(surfaceY)) return 0; // ambient
          if (y > static_cast<double>(oxideBottom)) {
            if (std::abs(x) < static_cast<double>(maskHalfWidth) &&
                y > static_cast<double>(surfaceY) - 0.05)
              return 3; // nitride mask
            return 2;   // oxide
          }
          return 1; // silicon
        });

    lastSurfaceY_ = surfaceY;
    lastOxideBottom_ = oxideBottom;
    std::cout << "[GeometryFieldCoupler] LOCOS-like mark: surfaceY=" << surfaceY
              << " oxideBottom=" << oxideBottom << " maskHalfW=" << maskHalfWidth
              << "\n";
  }

  /// Mark from depth interfaces only (1D depth stack, x ignored).
  /// depths are normalized 0..1 from top (surface) toward bulk.
  void markFromDepthInterfaces(NumericType oxideThicknessNorm = NumericType(0.15),
                               NumericType maskThicknessNorm = NumericType(0.05)) {
    if (!field_) return;
    field_->markMeshRegions([=](double /*x*/, double y, double /*z*/) -> int {
      // y is treated as depth from top (0=surface, 1=bulk) in profile mode
      if (y < static_cast<double>(maskThicknessNorm)) return 3; // mask / surface film
      if (y < static_cast<double>(maskThicknessNorm + oxideThicknessNorm)) return 2;
      return 1; // Si bulk
    });
    std::cout << "[GeometryFieldCoupler] Depth-stack mark oxideNorm="
              << oxideThicknessNorm << " maskNorm=" << maskThicknessNorm << "\n";
  }

  /// Warm-start field remap after geometry movement (oxide growth delta).
  void remapAfterGrowth(NumericType oxideDeltaNorm) {
    if (!field_) return;
    lastOxideBottom_ = std::max(NumericType(0), lastOxideBottom_ - oxideDeltaNorm);
    lastSurfaceY_ = std::min(NumericType(1), lastSurfaceY_ + oxideDeltaNorm * NumericType(0.5));
    field_->remapAfterGeometryUpdate();
    std::cout << "[GeometryFieldCoupler] Remap after growth delta=" << oxideDeltaNorm
              << " newOxideBottom=" << lastOxideBottom_ << "\n";
  }

  /// Live material-map coupling: mark from a dense material-id sampler
  /// (Domain MaterialMap / cell set query). Signature: matId = f(x,y,z).
  /// This is the Domain-facing path without requiring ViennaLS headers here.
  int markFromMaterialMap(
      const std::function<int(double x, double y, double z)>& materialAt,
      double x0 = -0.5, double x1 = 0.5, double y0 = 0.0, double y1 = 1.0,
      int nx = 32, int ny = 32) {
    if (!field_ || !materialAt) return 0;
    field_->initMeshFromBounds(x0, x1, y0, y1, nx, ny);
    field_->markMeshRegions(materialAt);
    // Count unique materials in profile MaterialID
    int nSi = 0, nOx = 0, nMask = 0, nAmb = 0;
    auto n = field_->getProfileSize();
    for (std::size_t i = 0; i < n; ++i) {
      int m = field_->getMaterialAtNormalizedDepth(
          static_cast<NumericType>(i) / static_cast<NumericType>(std::max<std::size_t>(1, n)));
      if (m == 1) nSi++;
      else if (m == 2) nOx++;
      else if (m == 3) nMask++;
      else nAmb++;
    }
    lastMapCounts_ = nSi + nOx + nMask + nAmb;
    std::cout << "[GeometryFieldCoupler] markFromMaterialMap bins: Si=" << nSi
              << " Ox=" << nOx << " Mask=" << nMask << " Amb=" << nAmb << "\n";
    return lastMapCounts_;
  }

  /// Convenience: Domain-like LS stack inferred from material layers at depths.
  /// materials ordered from surface to bulk, each with thickness in normalized y.
  int markFromLayerStack(
      const std::vector<std::pair<int, NumericType>>& layersFromSurface) {
    if (!field_) return 0;
    return markFromMaterialMap([=](double /*x*/, double y, double /*z*/) -> int {
      NumericType depth = static_cast<NumericType>(y);
      NumericType acc = 0;
      for (const auto& [mat, th] : layersFromSurface) {
        acc += th;
        if (depth < static_cast<double>(acc)) return mat;
      }
      return 1; // default Si bulk
    });
  }

  NumericType getLastSurfaceY() const { return lastSurfaceY_; }
  NumericType getLastOxideBottom() const { return lastOxideBottom_; }
  int getLastMapCounts() const { return lastMapCounts_; }

private:
  std::shared_ptr<PhysicsField<NumericType>> field_;
  NumericType lastSurfaceY_ = NumericType(0.5);
  NumericType lastOxideBottom_ = NumericType(0.3);
  int lastMapCounts_ = 0;
};

} // namespace viennaps
