#pragma once

/// MovingMeshHandler — mesh mutation for oxidation-style progress (ADR-0004).
///
/// Provides:
///  - Subdomain attribute relabeling (Si → SiO2) when progress ≥ threshold
///  - Optional interface-adjacent-only flip (production path)
///  - Boundary-node ALE-style displacement (smooth free surface)
///
/// Callable without a full process re-setup. MFEM paths gated by
/// VIENNAPS_HAS_MFEM.

#include <cmath>
#include <functional>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

/// Result of a mesh-move / relabel operation (for tests and logging).
struct MeshMoveResult {
  int elementsRelabeled = 0;
  int nodesDisplaced = 0;
  double maxDisplacement = 0.0;
};

class MovingMeshHandler {
public:
  void setInterfaceAdjacentOnly(bool on) { interfaceOnly_ = on; }
  bool interfaceAdjacentOnly() const { return interfaceOnly_; }

#ifdef VIENNAPS_HAS_MFEM
  /// Instance path: uses interfaceOnly_ flag set on this handler.
  MeshMoveResult relabel(mfem::ParMesh &mesh, int fromAttr, int toAttr,
                         double progress, double threshold) const {
    return relabelAttributes(mesh, fromAttr, toAttr, progress, threshold,
                             interfaceOnly_);
  }

  /// Subdomain relabeling (ADR-0004 idiom A): flip fromAttr → toAttr when
  /// oxidation progress ≥ threshold. Returns flipped element count.
  static MeshMoveResult relabelAttributes(mfem::ParMesh &mesh, int fromAttr,
                                          int toAttr, double progress,
                                          double threshold,
                                          bool interfaceAdjacentOnly = false) {
    MeshMoveResult r;
    if (progress < threshold)
      return r;

    std::vector<char> adjacent(static_cast<std::size_t>(mesh.GetNE()), 0);
    if (interfaceAdjacentOnly) {
      for (int f = 0; f < mesh.GetNumFaces(); ++f) {
        int eA = -1, eB = -1;
        mesh.GetFaceElements(f, &eA, &eB);
        if (eA < 0)
          continue;
        if (eB < 0) {
          if (eA < mesh.GetNE())
            adjacent[static_cast<std::size_t>(eA)] = 1;
          continue;
        }
        if (mesh.GetAttribute(eA) != mesh.GetAttribute(eB)) {
          adjacent[static_cast<std::size_t>(eA)] = 1;
          adjacent[static_cast<std::size_t>(eB)] = 1;
        }
      }
    }

    for (int e = 0; e < mesh.GetNE(); ++e) {
      if (mesh.GetAttribute(e) != fromAttr)
        continue;
      if (interfaceAdjacentOnly && !adjacent[static_cast<std::size_t>(e)])
        continue;
      mesh.SetAttribute(e, toAttr);
      ++r.elementsRelabeled;
    }
    return r;
  }

  /// ALE-style node displacement. `disp(x, u)` writes displacement into u
  /// given physical coordinates x (length SpaceDimension).
  static MeshMoveResult
  displaceNodes(mfem::ParMesh &mesh,
                const std::function<void(const double *x, double *u)> &disp) {
    MeshMoveResult r;
    const int sdim = mesh.SpaceDimension();
    const int nv = mesh.GetNV();
    std::vector<double> x(static_cast<std::size_t>(sdim), 0.0);
    std::vector<double> u(static_cast<std::size_t>(sdim), 0.0);

    mfem::GridFunction *nodes = mesh.GetNodes();
    if (!nodes)
      return r;
    for (int v = 0; v < nv; ++v) {
      for (int d = 0; d < sdim; ++d)
        x[static_cast<std::size_t>(d)] = (*nodes)(v * sdim + d);
      for (int d = 0; d < sdim; ++d)
        u[static_cast<std::size_t>(d)] = 0.0;
      disp(x.data(), u.data());
      double mag2 = 0.0;
      for (int d = 0; d < sdim; ++d)
        mag2 += u[static_cast<std::size_t>(d)] * u[static_cast<std::size_t>(d)];
      const double mag = std::sqrt(mag2);
      if (mag <= 0.0)
        continue;
      for (int d = 0; d < sdim; ++d)
        x[static_cast<std::size_t>(d)] += u[static_cast<std::size_t>(d)];
      for (int d = 0; d < sdim; ++d)
        (*nodes)(v * sdim + d) = x[static_cast<std::size_t>(d)];
      ++r.nodesDisplaced;
      if (mag > r.maxDisplacement)
        r.maxDisplacement = mag;
    }
    mesh.ExchangeFaceNbrData();
    return r;
  }

  /// Lift free-surface nodes with last coordinate ≥ threshold by `lift`.
  static MeshMoveResult liftFreeSurface(mfem::ParMesh &mesh, double coordThreshold,
                                        double lift) {
    const int sdim = mesh.SpaceDimension();
    return displaceNodes(mesh, [&](const double *x, double *u) {
      for (int d = 0; d < sdim; ++d)
        u[d] = 0.0;
      if (sdim >= 2 && x[sdim - 1] >= coordThreshold)
        u[sdim - 1] = lift;
    });
  }

  /// Laplacian node smoothing (Jacobi average of neighbors).
  static int laplacianSmooth(mfem::ParMesh &mesh, int iterations = 2,
                             double weight = 0.5) {
    if (iterations <= 0)
      return 0;
    const int nv = mesh.GetNV();
    const int sdim = mesh.SpaceDimension();
    std::vector<double> x(static_cast<std::size_t>(nv * sdim));
    std::vector<double> xnew(x.size());
    mfem::GridFunction *nodes = mesh.GetNodes();
    if (!nodes)
      return 0;
    for (int v = 0; v < nv; ++v)
      for (int d = 0; d < sdim; ++d)
        x[static_cast<std::size_t>(v * sdim + d)] = (*nodes)(v * sdim + d);

    // Element-local averaging graph.
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(nv));
    for (int el = 0; el < mesh.GetNE(); ++el) {
      mfem::Array<int> v;
      mesh.GetElementVertices(el, v);
      for (int i = 0; i < v.Size(); ++i)
        for (int j = 0; j < v.Size(); ++j)
          if (i != j)
            adj[static_cast<std::size_t>(v[i])].push_back(v[j]);
    }
    for (int it = 0; it < iterations; ++it) {
      xnew = x;
      for (int v = 0; v < nv; ++v) {
        const auto &nbr = adj[static_cast<std::size_t>(v)];
        if (nbr.empty())
          continue;
        std::vector<double> avg(static_cast<std::size_t>(sdim), 0.0);
        for (int n : nbr)
          for (int d = 0; d < sdim; ++d)
            avg[static_cast<std::size_t>(d)] +=
                x[static_cast<std::size_t>(n * sdim + d)];
        for (int d = 0; d < sdim; ++d) {
          avg[static_cast<std::size_t>(d)] /=
              static_cast<double>(nbr.size());
          const std::size_t id = static_cast<std::size_t>(v * sdim + d);
          xnew[id] = (1.0 - weight) * x[id] + weight * avg[static_cast<std::size_t>(d)];
        }
      }
      x.swap(xnew);
    }
    for (int v = 0; v < nv; ++v)
      for (int d = 0; d < sdim; ++d)
        (*nodes)(v * sdim + d) = x[static_cast<std::size_t>(v * sdim + d)];
    mesh.ExchangeFaceNbrData();
    return iterations;
  }

  /// Max element aspect proxy: max edge / min edge over elements.
  static double maxAspectRatio(mfem::ParMesh &mesh) {
    double worst = 1.0;
    mfem::GridFunction *nodes = mesh.GetNodes();
    if (!nodes)
      return worst;
    const int sdim = mesh.SpaceDimension();
    for (int e = 0; e < mesh.GetNE(); ++e) {
      mfem::Array<int> v;
      mesh.GetElementVertices(e, v);
      if (v.Size() < 2)
        continue;
      double minL = 1e300, maxL = 0.0;
      for (int i = 0; i < v.Size(); ++i) {
        double xi[3] = {0, 0, 0}, xj[3] = {0, 0, 0};
        for (int d = 0; d < sdim; ++d)
          xi[d] = (*nodes)(v[i] * sdim + d);
        for (int j = i + 1; j < v.Size(); ++j) {
          for (int d = 0; d < sdim; ++d)
            xj[d] = (*nodes)(v[j] * sdim + d);
          double L2 = 0.0;
          for (int d = 0; d < sdim; ++d) {
            const double dx = xi[d] - xj[d];
            L2 += dx * dx;
          }
          const double L = std::sqrt(L2);
          minL = std::min(minL, L);
          maxL = std::max(maxL, L);
        }
      }
      if (minL > 0.0)
        worst = std::max(worst, maxL / minL);
    }
    return worst;
  }

  /// Remesh trigger: true if aspect ratio exceeds threshold (skewness proxy).
  static bool needsRemesh(mfem::ParMesh &mesh, double aspectThreshold = 5.0) {
    return maxAspectRatio(mesh) > aspectThreshold;
  }
#endif

private:
  bool interfaceOnly_ = false;
};

} // namespace viennaps
