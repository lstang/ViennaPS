#pragma once

/// @file LevelSetToMesh.hpp
///
/// Phase 1 bridge from ViennaPS's level-set domain to MFEM's FEM mesh.
///
/// Design reference: MOOSE's `MFEMMesh`
/// (framework/include/mfem/mesh/MFEMMesh.h) — a wrapper around
/// `shared_ptr<mfem::ParMesh>` that exposes `shouldDisplace`,
/// `uniformRefinement`, and reports `nSubdomains = parmesh.attributes.Size()`.
/// The converter produces the same `mfem::Mesh` MOOSE would wrap, so a future
/// architecture (B) re-binding to MOOSE's `MFEMMesh` can consume the output
/// directly. The conforming cut-cell analog
/// (`CutMeshByLevelSetGeneratorBase::pointLevelSetRelation`,
/// `tet4ElemCutter`, `_generate_transition_layer`) is **deferred** to a later
/// phase — Phase 1 ships Cartesian + attribute-tagging only, since sliver-free
/// cut cells are only required once segregation accuracy at material interfaces
/// becomes the bottleneck.

#include "MeshAttributes.hpp"

#include <mfem.hpp>

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <psDomain.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace viennaps {

/// Result of a level-set → mesh conversion. Owns the produced mesh and the
/// attribute → material-name mapping. The mesh pointer is non-null on success.
template <class NumericType> struct LevelSetToMeshResult {
  std::unique_ptr<mfem::Mesh> mesh;
  MeshAttributes attributes;
};

/// LevelSetToMeshConverter converts a ViennaPS level-set domain into an
/// `mfem::Mesh` whose elements are tagged with material attributes.
///
/// Phase 1 scope (per ADR-0001 + task-4-brief):
///   - 2D only. A Cartesian triangular mesh (`MakeCartesian2D`) is built over
///     the level-set grid bounds and vertices are scaled to domain coordinates.
///   - Element attributes are assigned by evaluating each element's centroid
///     against the registered level sets: the material at the centroid is the
///     material of the innermost level set that contains it. (ViennaLS
///     convention: the level set is negative inside the material region.)
///   - D=3 is provided as a stub that throws — `MakeCartesian3D` will be added
///     in a later phase.
template <class NumericType, int D> class LevelSetToMeshConverter {
public:
  /// Convert `domain` to `{unique_ptr<mfem::Mesh>, MeshAttributes}`.
  /// Throws `std::runtime_error` if the domain has no level sets.
  LevelSetToMeshResult<NumericType>
  convert(const Domain<NumericType, D> &domain) const {
    throw std::runtime_error(
        "3D LevelSetToMeshConverter not implemented in Phase 1");
  }
};

/// 2D specialization: produces a triangular Cartesian MFEM mesh tagged from
/// the domain's material map.
template <class NumericType> class LevelSetToMeshConverter<NumericType, 2> {
public:
  LevelSetToMeshResult<NumericType>
  convert(const Domain<NumericType, 2> &domain) const {
    const auto &levelSets = domain.getLevelSets();
    if (levelSets.empty()) {
      throw std::runtime_error(
          "LevelSetToMeshConverter: domain has no level sets");
    }

    // Derive physical bounds + grid spacing from the first level set's grid.
    // psDomain guarantees all level sets share the same grid, so [0] is
    // sufficient and matches the level-set iteration order of getMaterialMap().
    const auto &grid = levelSets[0]->getGrid();
    const NumericType gridDelta =
        static_cast<NumericType>(grid.getGridDelta());

    // Physical bounds via index * gridDelta (hrleGrid::index2Coordinate).
    // Use getMinBounds/getMaxBounds (not getMinIndex/getMaxIndex) so axes with
    // INFINITE_BOUNDARY (e.g. the primary growth direction in 2D) are clamped
    // to the level-set's actual finite extent instead of INF_EXTENSION.
    const auto &minIdx = grid.getMinBounds();
    const auto &maxIdx = grid.getMaxBounds();

    const NumericType xMin = static_cast<NumericType>(minIdx[0]) * gridDelta;
    const NumericType xMax = static_cast<NumericType>(maxIdx[0]) * gridDelta;
    const NumericType yMin = static_cast<NumericType>(minIdx[1]) * gridDelta;
    const NumericType yMax = static_cast<NumericType>(maxIdx[1]) * gridDelta;

    const NumericType xExtent = xMax - xMin;
    const NumericType yExtent = yMax - yMin;
    if (xExtent <= 0.0 || yExtent <= 0.0 || gridDelta <= 0.0) {
      throw std::runtime_error(
          "LevelSetToMeshConverter: degenerate domain bounds");
    }

    // Cartesian element count: one MFEM element per grid cell per axis.
    const int nx = static_cast<int>(std::round(xExtent / gridDelta));
    const int ny = static_cast<int>(std::round(yExtent / gridDelta));
    if (nx <= 0 || ny <= 0) {
      throw std::runtime_error(
          "LevelSetToMeshConverter: non-positive element count");
    }

    // Build the mesh in [0, xExtent] x [0, yExtent], then translate vertices
    // by (xMin, yMin) to land in physical coordinates.
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::TRIANGLE,
                                    /*generate_edges*/ false,
                                    /*sx*/ static_cast<mfem::real_t>(xExtent),
                                    /*sy*/ static_cast<mfem::real_t>(yExtent)));

    mesh->Transform([xMin, yMin](const mfem::Vector &in, mfem::Vector &out) {
      out.SetSize(2);
      out(0) = in(0) + xMin;
      out(1) = in(1) + yMin;
    });

    // ---- Attribute tagging ----
    // For each element, evaluate its centroid against each level set. The
    // innermost (lowest-index) level set containing the centroid determines
    // the material. (ViennaLS sign convention: negative inside the material.)
    const auto materialMap = domain.getMaterialMap();
    MeshAttributes attrs;
    if (materialMap) {
      for (std::size_t i = 0; i < materialMap->size(); ++i) {
        // MFEM attribute IDs are 1-based; material map indices are 0-based.
        const int attr = static_cast<int>(i) + 1;
        attrs.setAttributeName(
            attr, MaterialMap::toString(materialMap->getMaterialAtIdx(i)));
      }
    }

    mfem::Vector center(2);
    for (int e = 0; e < mesh->GetNE(); ++e) {
      mesh->GetElementCenter(e, center);
      const int attr = attributeAtPoint(levelSets, grid, gridDelta, xMin, yMin,
                                        static_cast<NumericType>(center(0)),
                                        static_cast<NumericType>(center(1)));
      mesh->SetAttribute(e, attr);
    }
    mesh->SetAttributes();

    return {std::move(mesh), std::move(attrs)};
  }

private:
  using LsPtr = SmartPointer<viennals::Domain<NumericType, 2>>;

  /// Returns the 1-based material attribute at (x, y). Picks the lowest-index
  /// level set whose signed distance is negative at the nearest grid point
  /// (i.e. the centroid is inside that material region). If no level set
  /// contains the point, returns 1 (default bulk) so MFEM always sees a valid,
  /// positive attribute.
  static int attributeAtPoint(const std::vector<LsPtr> &levelSets,
                              const viennahrle::Grid<2> &grid,
                              NumericType gridDelta, NumericType xMin,
                              NumericType yMin, NumericType x,
                              NumericType y) {
    // Snap the centroid to the nearest grid index relative to (xMin, yMin).
    viennahrle::Index<2> idx;
    idx[0] = static_cast<viennahrle::IndexType>(
        std::round((x - xMin) / gridDelta));
    idx[1] = static_cast<viennahrle::IndexType>(
        std::round((y - yMin) / gridDelta));

    // Clamp into grid bounds (centroids of edge elements may sit on boundary).
    // Use bounds, not raw index range, so infinite axes stay queryable.
    idx[0] = std::clamp<viennahrle::IndexType>(idx[0], grid.getMinBounds(0),
                                               grid.getMaxBounds(0));
    idx[1] = std::clamp<viennahrle::IndexType>(idx[1], grid.getMinBounds(1),
                                               grid.getMaxBounds(1));

    for (std::size_t i = 0; i < levelSets.size(); ++i) {
      viennahrle::ConstSparseIterator<viennahrle::Domain<NumericType, 2>> it(
          levelSets[i]->getDomain(), idx);
      const NumericType value = it.getValue();
      if (value <= NumericType(0)) {
        return static_cast<int>(i) + 1;
      }
    }
    return 1;
  }
};

} // namespace viennaps
