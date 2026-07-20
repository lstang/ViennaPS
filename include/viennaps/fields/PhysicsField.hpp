#pragma once

/// PhysicsField - Container for bulk physical fields (dopants, defects, damage, stress).
///
/// Supports:
///   - Per-species 1D depth profiles (always available; used by CVODE packed state)
///   - Optional MFEM GridFunction storage when VIENNAPS_HAS_MFEM is defined
///   - packState / unpackState for multi-species SUNDIALS residuals
///   - Mesh region marking via classifier (geometry coupling skeleton)

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <functional>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
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
      speciesOrder_.push_back(name);
      profiles_[name] = std::vector<NumericType>(profileSize_, NumericType(0));
      totalDose_[name] = NumericType(0);
      std::cout << "[PhysicsField] Added species: " << name << std::endl;
    }
  }

  void setProfileSize(std::size_t n) {
    if (n < 4) n = 4;
    profileSize_ = n;
    for (auto& [name, prof] : profiles_) {
      prof.resize(profileSize_, NumericType(0));
    }
  }

  std::size_t getProfileSize() const { return profileSize_; }

  const std::vector<std::string>& getSpeciesOrder() const { return speciesOrder_; }

  void injectImplantProfile(const std::string& species,
                            const std::vector<NumericType>& profile) {
    addSpecies(species);
    totalDose_[species] += std::accumulate(profile.begin(), profile.end(), NumericType(0));

    auto& target = profiles_[species];
    if (!profile.empty()) {
      for (size_t i = 0; i < target.size(); ++i) {
        size_t src = i * profile.size() / target.size();
        if (src < profile.size()) {
          target[i] += profile[src];
        }
      }
    }
    clampProfileNonNegative(species);

#ifdef VIENNAPS_HAS_MFEM
    ensureGridFunction(species);
    if (gridFunctions_.count(species) && mesh_ && !profile.empty()) {
      mfem::GridFunction& gf = *gridFunctions_[species];
      mfem::FunctionCoefficient depthProf([&](const mfem::Vector& x) {
        double normY = (x.Size() > 1 ? x[1] : 0.0) / 2.0;
        size_t idx = std::min(profile.size() - 1, size_t(std::max(0.0, normY) * (profile.size() - 1)));
        return static_cast<double>(profile[idx]);
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

  const std::vector<NumericType>& getProfile(const std::string& species) const {
    static const std::vector<NumericType> empty;
    auto it = profiles_.find(species);
    return (it != profiles_.end()) ? it->second : empty;
  }

  /// Writable profile access (kernels / CVODE unpack).
  std::vector<NumericType>* getProfileMutable(const std::string& species) {
    auto it = profiles_.find(species);
    return (it != profiles_.end()) ? &it->second : nullptr;
  }

  NumericType getConcentration(const std::string& species,
                               NumericType normalizedDepth01) const {
#ifdef VIENNAPS_HAS_MFEM
    if (gridFunctions_.count(species) && mesh_) {
      mfem::Vector pt(mesh_->Dimension());
      pt = 0.0;
      if (pt.Size() > 1) pt[1] = static_cast<double>(normalizedDepth01) * 2.0;
      return static_cast<NumericType>((*gridFunctions_.at(species)).GetValue(pt));
    }
#endif
    auto it = profiles_.find(species);
    if (it == profiles_.end() || it->second.empty()) return NumericType(0);
    size_t idx = static_cast<size_t>(
        std::clamp(normalizedDepth01, NumericType(0), NumericType(0.999)) *
        (it->second.size() - 1));
    return it->second[idx];
  }

  void setConcentration(const std::string& species, std::size_t index, NumericType value) {
    addSpecies(species);
    auto& p = profiles_[species];
    if (index >= p.size()) return;
    p[index] = std::max(NumericType(0), value);
    recomputeTotalDose(species);
  }

  /// Pack all species profiles into a single state vector for SUNDIALS.
  /// Layout: speciesOrder_[0][0..N-1], speciesOrder_[1][0..N-1], ...
  std::vector<NumericType> packState() const {
    std::vector<NumericType> y;
    y.reserve(speciesOrder_.size() * profileSize_);
    for (const auto& name : speciesOrder_) {
      auto it = profiles_.find(name);
      if (it == profiles_.end()) {
        y.insert(y.end(), profileSize_, NumericType(0));
      } else {
        y.insert(y.end(), it->second.begin(), it->second.end());
      }
    }
    return y;
  }

  void unpackState(const std::vector<NumericType>& y) {
    std::size_t offset = 0;
    for (const auto& name : speciesOrder_) {
      addSpecies(name);
      auto& p = profiles_[name];
      p.resize(profileSize_, NumericType(0));
      for (std::size_t i = 0; i < profileSize_ && offset < y.size(); ++i, ++offset) {
        p[i] = std::max(NumericType(0), y[offset]);
      }
      recomputeTotalDose(name);
#ifdef VIENNAPS_HAS_MFEM
      ensureGridFunction(name);
      if (gridFunctions_.count(name)) {
        // Project 1D profile onto GF (depth ~ y)
        auto& gf = *gridFunctions_[name];
        const auto& prof = p;
        mfem::FunctionCoefficient depthProf([&](const mfem::Vector& x) {
          double normY = (x.Size() > 1 ? x[1] : 0.0) / 2.0;
          size_t idx = std::min(prof.size() - 1,
                                size_t(std::max(0.0, normY) * (prof.size() - 1)));
          return static_cast<double>(prof[idx]);
        });
        gf.ProjectCoefficient(depthProf);
      }
#endif
    }
  }

  std::size_t getStateSize() const { return speciesOrder_.size() * profileSize_; }

  /// Offset of species block in packed state, or npos if missing.
  std::size_t getSpeciesOffset(const std::string& species) const {
    for (std::size_t i = 0; i < speciesOrder_.size(); ++i) {
      if (speciesOrder_[i] == species) return i * profileSize_;
    }
    return static_cast<std::size_t>(-1);
  }

  bool hasSpecies(const std::string& species) const {
    return species_.find(species) != species_.end();
  }

#ifdef VIENNAPS_HAS_MFEM
  mfem::GridFunction* getGridFunction(const std::string& species) {
    auto it = gridFunctions_.find(species);
    return (it != gridFunctions_.end()) ? it->second.get() : nullptr;
  }
  mfem::FiniteElementSpace* getFESpace() { return fespace_.get(); }
  mfem::Mesh* getMesh() { return mesh_.get(); }
#endif

  void evolve(NumericType dt) {
    std::cout << "[PhysicsField] Evolving all fields for dt=" << dt << " (stub)\n";
  }

  void scaleProfile(const std::string& species, NumericType factor) {
    auto it = profiles_.find(species);
    if (it != profiles_.end()) {
      factor = std::max(NumericType(0), factor);
      for (auto& v : it->second) {
        v = std::max(NumericType(0), v * factor);
      }
    }
    recomputeTotalDose(species);

#ifdef VIENNAPS_HAS_MFEM
    if (gridFunctions_.count(species)) {
      *gridFunctions_[species] *= static_cast<double>(factor);
    }
#endif
  }

  /// Recompute totalDose from the stored profile (after direct profile mutation).
  void refreshDose(const std::string& species) { recomputeTotalDose(species); }

  void remapAfterGeometryUpdate() {
#ifdef VIENNAPS_HAS_MFEM
    if (mesh_) {
      std::cout << "[PhysicsField] Remapping fields after geometry update (warm-start)\n";
    }
#endif
    std::cout << "[PhysicsField] Remap after geometry update (profiles retained)\n";
  }

  /// Initialize / ensure MFEM structured mesh (public for geometry coupling).
  void initMeshFromBounds(double xMin, double xMax, double yMin, double yMax,
                          int nx = 32, int ny = 32) {
#ifdef VIENNAPS_HAS_MFEM
    initMFEMMeshFromBounds(xMin, xMax, yMin, yMax, nx, ny);
#else
    (void)xMin; (void)xMax; (void)yMin; (void)yMax; (void)nx; (void)ny;
    std::cout << "[PhysicsField] initMeshFromBounds: MFEM not available, profile-only mode\n";
#endif
  }

  /// Mark mesh elements by material id using a point classifier (geometry coupling).
  /// Without MFEM, stores a 1D material id profile for the depth axis.
  void markMeshRegions(const std::function<int(double x, double y, double z)>& classifier) {
#ifdef VIENNAPS_HAS_MFEM
    ensureMFEMMesh();
    markMeshRegionsMFEM([&](const mfem::Vector& pt) {
      double x = pt.Size() > 0 ? pt[0] : 0.0;
      double y = pt.Size() > 1 ? pt[1] : 0.0;
      double z = pt.Size() > 2 ? pt[2] : 0.0;
      return classifier(x, y, z);
    });
#else
    // Profile-only material map along depth
    addSpecies("MaterialID");
    auto& mid = profiles_["MaterialID"];
    for (std::size_t i = 0; i < mid.size(); ++i) {
      double y = static_cast<double>(i) / static_cast<double>(mid.size());
      mid[i] = static_cast<NumericType>(classifier(0.0, y, 0.0));
    }
    recomputeTotalDose("MaterialID");
    std::cout << "[PhysicsField] Marked " << mid.size()
              << " profile bins with material classifier (no MFEM)\n";
#endif
  }

  int getMaterialAtNormalizedDepth(NumericType normalizedDepth01) const {
    auto it = profiles_.find("MaterialID");
    if (it == profiles_.end() || it->second.empty()) return 0;
    size_t idx = static_cast<size_t>(
        std::clamp(normalizedDepth01, NumericType(0), NumericType(0.999)) *
        (it->second.size() - 1));
    return static_cast<int>(it->second[idx]);
  }

private:
  std::map<std::string, int> species_;
  std::vector<std::string> speciesOrder_;
  std::map<std::string, NumericType> totalDose_;
  std::map<std::string, std::vector<NumericType>> profiles_;
  std::size_t profileSize_ = 128;

  void recomputeTotalDose(const std::string& species) {
    auto it = profiles_.find(species);
    if (it == profiles_.end()) {
      totalDose_[species] = NumericType(0);
      return;
    }
    totalDose_[species] =
        std::accumulate(it->second.begin(), it->second.end(), NumericType(0));
  }

  void clampProfileNonNegative(const std::string& species) {
    auto it = profiles_.find(species);
    if (it == profiles_.end()) return;
    for (auto& v : it->second) v = std::max(NumericType(0), v);
    recomputeTotalDose(species);
  }

#ifdef VIENNAPS_HAS_MFEM
  std::unique_ptr<mfem::Mesh> mesh_;
  std::unique_ptr<mfem::H1_FECollection> fec_;
  std::unique_ptr<mfem::FiniteElementSpace> fespace_;
  std::map<std::string, std::unique_ptr<mfem::GridFunction>> gridFunctions_;

  void ensureMFEMMesh(int dim = 2, int nx = 64, int ny = 64) {
    if (mesh_) return;
    if (dim == 2) {
      mesh_ = std::make_unique<mfem::Mesh>(
          mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL));
    } else {
      mesh_ = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
          nx, ny, std::max(4, nx / 2), mfem::Element::HEXAHEDRON));
    }
    fec_ = std::make_unique<mfem::H1_FECollection>(1, mesh_->Dimension());
    fespace_ = std::make_unique<mfem::FiniteElementSpace>(mesh_.get(), fec_.get());
    std::cout << "[PhysicsField] Created MFEM structured mesh (dim=" << dim
              << "), dofs=" << fespace_->GetTrueVSize() << std::endl;
  }

  void ensureGridFunction(const std::string& species) {
    ensureMFEMMesh();
    if (gridFunctions_.find(species) == gridFunctions_.end()) {
      auto gf = std::make_unique<mfem::GridFunction>(fespace_.get());
      *gf = 0.0;
      gridFunctions_[species] = std::move(gf);
    }
  }

  void initMFEMMeshFromBounds(double xMin, double xMax, double yMin, double yMax,
                              int nx = 64, int ny = 64) {
    mesh_ = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(
        nx, ny, mfem::Element::QUADRILATERAL, false, xMax - xMin, yMax - yMin));
    mfem::Vector shift(2);
    shift[0] = xMin;
    shift[1] = yMin;
    mesh_->Transform([&](const mfem::Vector& x, mfem::Vector& p) {
      p = x;
      p[0] += xMin;
      p[1] += yMin;
    });
    fec_ = std::make_unique<mfem::H1_FECollection>(1, 2);
    fespace_ = std::make_unique<mfem::FiniteElementSpace>(mesh_.get(), fec_.get());
    std::cout << "[PhysicsField] MFEM mesh initialized from bounds, dofs="
              << fespace_->GetTrueVSize() << std::endl;
  }

  void markMeshRegionsMFEM(std::function<int(const mfem::Vector&)> classifier) {
    if (!mesh_) return;
    int nMarked = 0;
    for (int i = 0; i < mesh_->GetNE(); ++i) {
      mfem::Element* el = mesh_->GetElement(i);
      mfem::Array<int> verts;
      el->GetVertices(verts);
      mfem::Vector center(mesh_->Dimension());
      center = 0.0;
      for (int v = 0; v < verts.Size(); ++v) {
        const double* coords = mesh_->GetVertex(verts[v]);
        for (int d = 0; d < mesh_->Dimension(); ++d)
          center[d] += coords[d];
      }
      center /= verts.Size();
      int mat = classifier(center);
      mesh_->SetAttribute(i, mat);
      nMarked++;
    }
    // Also store MaterialID profile for depth coupling
    addSpecies("MaterialID");
    auto& mid = profiles_["MaterialID"];
    for (std::size_t i = 0; i < mid.size(); ++i) {
      double y = static_cast<double>(i) / static_cast<double>(mid.size());
      mfem::Vector pt(mesh_->Dimension());
      pt = 0.0;
      if (pt.Size() > 1) pt[1] = y * 2.0;
      mid[i] = static_cast<NumericType>(classifier(pt));
    }
    recomputeTotalDose("MaterialID");
    std::cout << "[PhysicsField] Marked " << nMarked
              << " mesh elements using classifier.\n";
  }
#endif
};

} // namespace viennaps
