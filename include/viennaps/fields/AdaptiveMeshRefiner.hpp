#pragma once

/// Adaptive mesh refinement helpers (Phase 11 skeleton).

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

struct RefinementBox {
  double x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  bool contains(double x, double y) const {
    return x >= x0 && x <= x1 && y >= y0 && y <= y1;
  }
};

class MeshQualityEstimator {
public:
  /// Relative difference criterion between two scalar fields.
  static std::vector<double> relativeDifference(const std::vector<double> &a,
                                                const std::vector<double> &b) {
    const std::size_t n = std::min(a.size(), b.size());
    std::vector<double> e(n);
    for (std::size_t i = 0; i < n; ++i) {
      const double denom = std::max({std::abs(a[i]), std::abs(b[i]), 1e-30});
      e[i] = std::abs(a[i] - b[i]) / denom;
    }
    return e;
  }

  /// Simple gradient magnitude on a 1D grid.
  static std::vector<double> gradient1D(const std::vector<double> &u,
                                        double dx) {
    std::vector<double> g(u.size(), 0.0);
    if (u.size() < 2 || dx <= 0)
      return g;
    g[0] = std::abs(u[1] - u[0]) / dx;
    for (std::size_t i = 1; i + 1 < u.size(); ++i)
      g[i] = std::abs(u[i + 1] - u[i - 1]) / (2 * dx);
    g.back() = std::abs(u.back() - u[u.size() - 2]) / dx;
    return g;
  }

  static std::vector<double> logarithmic(const std::vector<double> &u) {
    std::vector<double> out(u.size());
    for (std::size_t i = 0; i < u.size(); ++i)
      out[i] = std::log10(std::max(u[i], 1e-30));
    return out;
  }

  static std::vector<double> asinh(const std::vector<double> &u) {
    std::vector<double> out(u.size());
    for (std::size_t i = 0; i < u.size(); ++i)
      out[i] = std::asinh(u[i]);
    return out;
  }
};

class AdaptiveMeshRefiner {
public:
  void setMaxLevel(int l) { maxLevel_ = l; }
  void setThreshold(double t) { threshold_ = t; }
  int maxLevel() const { return maxLevel_; }

  /// Mark elements (or bins) where indicator > threshold.
  static std::vector<int> mark(const std::vector<double> &indicator,
                               double threshold) {
    std::vector<int> marks;
    for (std::size_t i = 0; i < indicator.size(); ++i)
      if (indicator[i] > threshold)
        marks.push_back(static_cast<int>(i));
    return marks;
  }

  /// Uniform scale factor for mesh spacing.
  static double uniformScale(double h, double factor) {
    return h / std::max(factor, 1e-12);
  }

#ifdef VIENNAPS_HAS_MFEM
  /// Static box refinement: return element indices whose centers are in box.
  // Non-const mesh: MFEM GetElementCenter is non-const.
  static std::vector<int> markBox(mfem::Mesh &mesh, const RefinementBox &box) {
    std::vector<int> ids;
    for (int e = 0; e < mesh.GetNE(); ++e) {
      mfem::Vector c;
      mesh.GetElementCenter(e, c);
      if (box.contains(c(0), c.Size() > 1 ? c(1) : 0.0))
        ids.push_back(e);
    }
    return ids;
  }
#endif

private:
  int maxLevel_ = 3;
  double threshold_ = 0.1;
};

} // namespace viennaps
