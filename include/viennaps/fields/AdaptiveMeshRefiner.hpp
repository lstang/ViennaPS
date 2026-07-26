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

  /// Local refinement via MFEM GeneralRefinement (list of element ids).
  static int refineMarked(mfem::Mesh &mesh, const std::vector<int> &elemIds) {
    if (elemIds.empty())
      return 0;
    mfem::Array<int> el_to_refine;
    for (int id : elemIds)
      if (id >= 0 && id < mesh.GetNE())
        el_to_refine.Append(id);
    if (el_to_refine.Size() == 0)
      return 0;
    const int n = el_to_refine.Size();
    mesh.GeneralRefinement(el_to_refine);
    return n;
  }

  /// Refine mesh and prolongate all GridFunctions to the new space (MFEM native).
  /// Returns number of elements refined.
  static int refineMarkedWithProlongation(mfem::Mesh &mesh,
                                          mfem::FiniteElementSpace &fes,
                                          const std::vector<int> &elemIds,
                                          std::vector<mfem::GridFunction *> &gfs) {
    if (elemIds.empty())
      return 0;
    mfem::Array<int> el_to_refine;
    for (int id : elemIds)
      if (id >= 0 && id < mesh.GetNE())
        el_to_refine.Append(id);
    if (el_to_refine.Size() == 0)
      return 0;

    mesh.GeneralRefinement(el_to_refine);
    fes.Update(true);
    for (auto *gf : gfs)
      gf->Update();
    fes.UpdatesFinished();
    return static_cast<int>(el_to_refine.Size());
  }

  /// Gradient-based marking: elements with |∇u| ≥ fraction * max|∇u|.
  static std::vector<int>
  markByGradient(mfem::Mesh &mesh, const mfem::GridFunction &u,
                 double fractionOfMax = 0.5) {
    std::vector<double> g(static_cast<std::size_t>(mesh.GetNE()), 0.0);
    double gmax = 0.0;
    for (int e = 0; e < mesh.GetNE(); ++e) {
      mfem::ElementTransformation *T = mesh.GetElementTransformation(e);
      const mfem::IntegrationRule *ir =
          &mfem::IntRules.Get(mesh.GetElementBaseGeometry(e), 2);
      double local = 0.0;
      for (int i = 0; i < ir->GetNPoints(); ++i) {
        const mfem::IntegrationPoint &ip = ir->IntPoint(i);
        T->SetIntPoint(&ip);
        mfem::Vector grad;
        u.GetGradient(*T, grad);
        local = std::max(local, grad.Norml2());
      }
      g[static_cast<std::size_t>(e)] = local;
      gmax = std::max(gmax, local);
    }
    const double thr = fractionOfMax * gmax;
    std::vector<int> marks;
    for (int e = 0; e < mesh.GetNE(); ++e)
      if (g[static_cast<std::size_t>(e)] >= thr && thr > 0.0)
        marks.push_back(e);
    return marks;
  }

  /// Project GridFunction from old space onto new space after mesh refine.
  /// Caller must rebuild FiniteElementSpace on the refined mesh first; this
  /// helper does nodal injection when sizes match, else zero-fills.
  static void transferField(const mfem::GridFunction &oldGf,
                            mfem::GridFunction &newGf) {
    if (oldGf.Size() == newGf.Size()) {
      newGf = oldGf;
      return;
    }
    newGf = 0.0;
    const int n = std::min(oldGf.Size(), newGf.Size());
    for (int i = 0; i < n; ++i)
      newGf(i) = oldGf(i);
  }

  /// ZZ-style element indicator: |∇u − average(∇u)| proxy via local |∇u|.
  static std::vector<double> zzIndicator(mfem::Mesh &mesh,
                                         const mfem::GridFunction &u) {
    std::vector<double> ind(static_cast<std::size_t>(mesh.GetNE()), 0.0);
    for (int e = 0; e < mesh.GetNE(); ++e) {
      mfem::ElementTransformation *T = mesh.GetElementTransformation(e);
      const mfem::IntegrationRule *ir =
          &mfem::IntRules.Get(mesh.GetElementBaseGeometry(e), 2);
      double gsum = 0.0;
      int nq = 0;
      for (int i = 0; i < ir->GetNPoints(); ++i) {
        const mfem::IntegrationPoint &ip = ir->IntPoint(i);
        T->SetIntPoint(&ip);
        mfem::Vector grad;
        u.GetGradient(*T, grad);
        gsum += grad.Norml2();
        ++nq;
      }
      ind[static_cast<std::size_t>(e)] = (nq > 0) ? gsum / nq : 0.0;
    }
    return ind;
  }

  /// Threshold refiner: mark elements with indicator ≥ threshold.
  static std::vector<int> thresholdRefine(const std::vector<double> &indicator,
                                          double threshold) {
    return mark(indicator, threshold);
  }

  /// Coarsen proxy: mark elements with indicator ≤ threshold for derefinement
  /// candidate list (caller applies NCMesh derefines if available).
  static std::vector<int> markDerefine(const std::vector<double> &indicator,
                                       double threshold) {
    std::vector<int> marks;
    for (std::size_t i = 0; i < indicator.size(); ++i)
      if (indicator[i] <= threshold)
        marks.push_back(static_cast<int>(i));
    return marks;
  }

  /// Attempt uniform coarsening by removing finest NC level when mesh allows;
  /// returns number of elements before−after (0 if unsupported).
  static int tryDerefine(mfem::Mesh &mesh) {
    const int ne0 = mesh.GetNE();
    // Serial MFEM Cartesian meshes have no NC derefines; no-op returns 0.
    // Hook kept so runtime AMR can call a single API.
    (void)mesh;
    return 0 * ne0;
  }
#endif

private:
  int maxLevel_ = 3;
  double threshold_ = 0.1;
};

} // namespace viennaps
