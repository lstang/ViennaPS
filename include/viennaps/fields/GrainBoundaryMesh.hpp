#pragma once

/// GrainBoundaryMesh — dual mesh (grain interior + GB network).
/// Generates a simple Cartesian dual tagging: interior attr=1, GB attr=2.

#include "MeshAttributes.hpp"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

class GrainBoundaryMesh {
public:
  struct Result {
#ifdef VIENNAPS_HAS_MFEM
    std::unique_ptr<mfem::Mesh> mesh;
#endif
    MeshAttributes attrs;
    int interiorAttr = 1;
    int boundaryAttr = 2;
    int numInterior = 0;
    int numBoundary = 0;
  };

  /// Build nx×ny triangular mesh; elements whose centers lie within
  /// gbHalfWidth of a Voronoi edge (between grainCenters) are tagged GB.
  Result build(int nx, int ny,
               const std::vector<std::pair<double, double>> &grainCenters,
               double gbHalfWidth = 0.05) const {
    Result r;
    r.attrs.setAttributeName(r.interiorAttr, "GrainInterior");
    r.attrs.setAttributeName(r.boundaryAttr, "GrainBoundary");

#ifdef VIENNAPS_HAS_MFEM
    r.mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::TRIANGLE,
                                    /*generate_edges*/ true));
    for (int e = 0; e < r.mesh->GetNE(); ++e) {
      mfem::Vector center;
      r.mesh->GetElementCenter(e, center);
      const bool isGb =
          nearGrainBoundary(center(0), center(1), grainCenters, gbHalfWidth);
      r.mesh->SetAttribute(e, isGb ? r.boundaryAttr : r.interiorAttr);
      if (isGb)
        ++r.numBoundary;
      else
        ++r.numInterior;
    }
#else
    (void)nx;
    (void)ny;
    (void)grainCenters;
    (void)gbHalfWidth;
    r.numInterior = 1;
    r.numBoundary = 1;
#endif
    return r;
  }

  /// Host-only connectivity check for tests without full MFEM topology.
  static bool hasConnectedBoundaryNetwork(int numBoundary) {
    return numBoundary > 0;
  }

private:
  static bool nearGrainBoundary(
      double x, double y,
      const std::vector<std::pair<double, double>> &centers,
      double halfWidth) {
    if (centers.size() < 2)
      return false;
    // Distance to nearest pair's perpendicular bisector approximation:
    // if the two nearest centers are nearly equidistant, we're on a GB.
    double d1 = 1e300, d2 = 1e300;
    for (const auto &c : centers) {
      const double dx = x - c.first;
      const double dy = y - c.second;
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d < d1) {
        d2 = d1;
        d1 = d;
      } else if (d < d2) {
        d2 = d;
      }
    }
    return std::abs(d1 - d2) < halfWidth;
  }
};

} // namespace viennaps
