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
  MeshMoveResult relabel(mfem::Mesh &mesh, int fromAttr, int toAttr,
                         double progress, double threshold) const {
    return relabelAttributes(mesh, fromAttr, toAttr, progress, threshold,
                             interfaceOnly_);
  }

  /// Subdomain relabeling (ADR-0004 idiom A): flip fromAttr → toAttr when
  /// oxidation progress ≥ threshold. Returns flipped element count.
  static MeshMoveResult relabelAttributes(mfem::Mesh &mesh, int fromAttr,
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
  displaceNodes(mfem::Mesh &mesh,
                const std::function<void(const double *x, double *u)> &disp) {
    MeshMoveResult r;
    const int sdim = mesh.SpaceDimension();
    const int nv = mesh.GetNV();
    std::vector<double> x(static_cast<std::size_t>(sdim), 0.0);
    std::vector<double> u(static_cast<std::size_t>(sdim), 0.0);

    for (int v = 0; v < nv; ++v) {
      mesh.GetNode(v, x.data());
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
      mesh.SetNode(v, x.data());
      ++r.nodesDisplaced;
      if (mag > r.maxDisplacement)
        r.maxDisplacement = mag;
    }
    mesh.NodesUpdated();
    return r;
  }

  /// Lift free-surface nodes with last coordinate ≥ threshold by `lift`.
  static MeshMoveResult liftFreeSurface(mfem::Mesh &mesh, double coordThreshold,
                                        double lift) {
    const int sdim = mesh.SpaceDimension();
    return displaceNodes(mesh, [&](const double *x, double *u) {
      for (int d = 0; d < sdim; ++d)
        u[d] = 0.0;
      if (sdim >= 2 && x[sdim - 1] >= coordThreshold)
        u[sdim - 1] = lift;
    });
  }
#endif

private:
  bool interfaceOnly_ = false;
};

} // namespace viennaps
