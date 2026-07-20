#pragma once

/// PhysicsField - Container for all bulk physical fields (dopants, defects, damage, stress).
///
/// Current stage: functional stub supporting total-dose accumulation + per-species
/// 1D profile storage. This is sufficient for early 2D implant + simple diffusion demos.
///
/// Next stages:
///   - Own or reference an MFEM mesh + GridFunction/ParGridFunction per species (when VIENNAPS_HAS_MFEM)
///   - Couple to level-set geometry (band-limited active region)
///   - Hand off to SUNDIALS time integrators + kernel-assembled residuals
///   - Support vector/tensor fields for stress

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <iostream>
#include <numeric>
#include <algorithm>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#include <functional>
#endif

namespace viennaps {

template <class NumericType>
class PhysicsField {
public:
  PhysicsField() = default;

  void addSpecies(const std::string& name) {
    if (species_.find(name) == species_.end()) {
      int id = static_cast<int>(species_.size());
      species_[name] = id;
      // Allocate a simple profile storage for current stub mode
      profiles_[name] = std::vector<NumericType>(128, NumericType(0));
      std::cout << "[PhysicsField] Added species: " << name << std::endl;
    }
  }

  void injectImplantProfile(const std::string& species, const std::vector<NumericType>& profile) {
    addSpecies(species);
    totalDose_[species] += std::accumulate(profile.begin(), profile.end(), NumericType(0));

    // Store a simple resampled version into our internal profile (demo only)
    auto& target = profiles_[species];
    if (!profile.empty()) {
      for (size_t i = 0; i < target.size(); ++i) {
        size_t src = i * profile.size() / target.size();
        if (src < profile.size()) {
          target[i] += profile[src];
        }
      }
    }

#ifdef VIENNAPS_HAS_MFEM
    // Project 1D depth profile onto MFEM GridFunction (assume depth ~ y for 2D demo mesh)
    ensureGridFunction(species);
    if (gridFunctions_.count(species) && mesh_ && !profile.empty()) {
      mfem::GridFunction& gf = *gridFunctions_[species];
      NumericType peak = 0;
      for (auto v : profile) if (v > peak) peak = v;
      // Simple: set GF to a function of y (depth)
      mfem::FunctionCoefficient depthProf([&](const mfem::Vector& x) {
        // x[1] is y, map to profile
        double normY = (x.Size() > 1 ? x[1] : 0.0) / 2.0; // rough scale
        size_t idx = std::min(profile.size()-1, size_t(normY * (profile.size()-1)));
        return profile[idx];
      });
      gf.ProjectCoefficient(depthProf);
    }
#endif

    std::cout << "[PhysicsField] Injected " << species
              << "  totalDose=" << totalDose_[species] << std::endl;
  }

  NumericType getTotalDose(const std::string& species) const {
    auto it = totalDose_.find(species);
    return it != totalDose_.end() ? it->second : NumericType(0);
  }

  // Return a reference to the current (stub) depth profile for a species.
  // In full MFEM mode this will be replaced by projection/interpolation from GridFunction.
  const std::vector<NumericType>& getProfile(const std::string& species) const {
    static const std::vector<NumericType> empty;
    auto it = profiles_.find(species);
    return (it != profiles_.end()) ? it->second : empty;
  }

  // Very simple "get concentration at normalized depth" (0..1) using the profile.
  // Used by stub kernels.
  NumericType getConcentration(const std::string& species, NumericType normalizedDepth01) const {
#ifdef VIENNAPS_HAS_MFEM
    if (gridFunctions_.count(species) && mesh_) {
      // Evaluate at a point along depth (y)
      mfem::Vector pt(mesh_->Dimension());
      pt = 0.0;
      if (pt.Size() > 1) pt[1] = normalizedDepth01 * 2.0; // rough
      return (*gridFunctions_.at(species)).GetValue(pt);
    }
#endif
    auto it = profiles_.find(species);
    if (it == profiles_.end() || it->second.empty()) return NumericType(0);
    size_t idx = static_cast<size_t>(std::clamp(normalizedDepth01, NumericType(0), NumericType(0.999)) *
                                     (it->second.size() - 1));
    return it->second[idx];
  }

#ifdef VIENNAPS_HAS_MFEM
  mfem::GridFunction* getGridFunction(const std::string& species) {
    auto it = gridFunctions_.find(species);
    return (it != gridFunctions_.end()) ? it->second.get() : nullptr;
  }
  mfem::FiniteElementSpace* getFESpace() { return fespace_.get(); }
  mfem::Mesh* getMesh() { return mesh_.get(); }
#endif

  // Simple evolve stub - will be replaced by SUNDIALS + kernel system
  void evolve(NumericType dt) {
    std::cout << "[PhysicsField] Evolving all fields for dt=" << dt << " (stub)\n";
    // In current stub we only support the BasicDiffusion hack via inject.
    // Real evolution happens inside concrete PhysicsKernels.
  }

  // Future: allow external code (kernels/adapter) to directly scale a profile (demo helper)
  void scaleProfile(const std::string& species, NumericType factor) {
    auto it = profiles_.find(species);
    if (it != profiles_.end()) {
      for (auto& v : it->second) v *= factor;
    }
    auto td = totalDose_.find(species);
    if (td != totalDose_.end()) td->second *= factor;

#ifdef VIENNAPS_HAS_MFEM
    if (gridFunctions_.count(species)) {
      *gridFunctions_[species] *= factor;
    }
#endif
  }

  // Remap fields after geometry update (warm start for new mesh or LS change)
  void remapAfterGeometryUpdate() {
#ifdef VIENNAPS_HAS_MFEM
    if (mesh_) {
      std::cout << "[PhysicsField] Remapping fields after geometry update (warm-start stub)\n";
      // In full: interpolate old GF to new mesh, or project from LS
    }
#endif
  }

  // TODO full:
  //   MFEM: Mesh, FiniteElementSpace, map<string, GridFunction> or ParGridFunction
  //   SUNDIALS: N_Vector wrappers around the dofs
  //   Geometry sync: mark active elements from level-set bands

private:
  std::map<std::string, int> species_;
  std::map<std::string, NumericType> totalDose_;
  std::map<std::string, std::vector<NumericType>> profiles_;

#ifdef VIENNAPS_HAS_MFEM
  // MFEM backing (structured Cartesian mesh preferred for initial LS coupling)
  std::unique_ptr<mfem::Mesh> mesh_;
  std::unique_ptr<mfem::H1_FECollection> fec_;
  std::unique_ptr<mfem::FiniteElementSpace> fespace_;
  std::map<std::string, std::unique_ptr<mfem::GridFunction>> gridFunctions_;

  void ensureMFEMMesh(int dim = 2, int nx = 64, int ny = 64) {
    if (mesh_) return;
    if (dim == 2) {
      mesh_ = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL));
    } else {
      mesh_ = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(nx, ny, std::max(4, nx/2), mfem::Element::HEXAHEDRON));
    }
    fec_ = std::make_unique<mfem::H1_FECollection>(1, mesh_->Dimension());
    fespace_ = std::make_unique<mfem::FiniteElementSpace>(mesh_.get(), fec_.get());
    std::cout << "[PhysicsField] Created MFEM structured mesh (dim=" << dim << "), dofs=" << fespace_->GetTrueVSize() << std::endl;
  }

  void ensureGridFunction(const std::string& species) {
    ensureMFEMMesh();
    if (gridFunctions_.find(species) == gridFunctions_.end()) {
      auto gf = std::make_unique<mfem::GridFunction>(fespace_.get());
      *gf = 0.0;
      gridFunctions_[species] = std::move(gf);
    }
  }

  // Initialize a structured mesh covering the given bounds (for coupling to LS domain)
  void initMFEMMeshFromBounds(double xMin, double xMax, double yMin, double yMax, int nx = 64, int ny = 64) {
    // Simple 2D for now; extend for 3D later
    mesh_ = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, false, xMax-xMin, yMax-yMin));
    // Shift to origin if needed (MakeCartesian starts at 0,0)
    mfem::Vector shift(2); shift[0] = xMin; shift[1] = yMin;
    mesh_->MoveNodes(shift);  // translate
    fec_ = std::make_unique<mfem::H1_FECollection>(1, 2);
    fespace_ = std::make_unique<mfem::FiniteElementSpace>(mesh_.get(), fec_.get());
    std::cout << "[PhysicsField] MFEM mesh initialized from bounds, dofs=" << fespace_->GetTrueVSize() << std::endl;
  }

  // Mark regions for the background mesh using a classifier function (geometry coupling).
  // classifier(point) -> material id
  void markMeshRegions(std::function<int(const mfem::Vector&)> classifier) {
    if (!mesh_) return;
    int nMarked = 0;
    for (int i = 0; i < mesh_->GetNE(); ++i) {
      mfem::Element* el = mesh_->GetElement(i);
      mfem::Array<int> verts;
      el->GetVertices(verts);
      mfem::Vector center(mesh_->Dimension());
      center = 0.0;
      for (int v = 0; v < verts.Size(); ++v) {
        mfem::Vector pt;
        mesh_->GetVertex(verts[v], pt);
        center += pt;
      }
      center /= verts.Size();
      int mat = classifier(center);
      mesh_->SetAttribute(i, mat);
      nMarked++;
    }
    std::cout << "[PhysicsField] Marked " << nMarked << " mesh elements using classifier.\n";
  }

  mfem::GridFunction* getMaterialIDGridFunction() {
    if (!mesh_ || !fespace_) return nullptr;
    static const std::string matIDName = "MaterialID";
    if (gridFunctions_.find(matIDName) == gridFunctions_.end()) {
      auto gf = std::make_unique<mfem::GridFunction>(fespace_.get());
      *gf = 0.0;
      gridFunctions_[matIDName] = std::move(gf);
    }
    return gridFunctions_[matIDName].get();
  }
#endif
};

} // namespace viennaps