#pragma once

#include "DiffusionModel.hpp"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <algorithm>

namespace viennaps {

template <class NumericType>
class DiffusionPhysics {
public:
  // ---- Species registration (MOOSE PhysicsBase::saveSolverVariableName) ----
  void addSpecies(const std::string& name) {
    if (!hasSpecies(name)) {
      species_.push_back(name);
      solverVariables_.insert(name);
    }
  }

  bool hasSpecies(const std::string& name) const {
    return std::find(species_.begin(), species_.end(), name)
           != species_.end();
  }

  /// MOOSE PhysicsBase::variableExists analog (verified PhysicsBase.h:152).
  bool variableExists(const std::string& name) const {
    return solverVariables_.find(name) != solverVariables_.end();
  }

  int numSpecies() const { return static_cast<int>(species_.size()); }
  const std::vector<std::string>& speciesNames() const { return species_; }

  // ---- Model registration ----
  void addModel(std::shared_ptr<DiffusionModel<NumericType>> m) {
    models_.push_back(m);
  }
  int numModels() const { return static_cast<int>(models_.size()); }
  const std::vector<std::shared_ptr<DiffusionModel<NumericType>>>&
  models() const { return models_; }

  // ---- Composition gatekeeper (MOOSE PhysicsBase::shouldCreateTimeDerivative,
  //      verified PhysicsBase.h:237). Returns true only the first time a
  //      species' time derivative is requested within the current solve,
  //      so composing FermiDiffusion + CddDiffusion on the same species
  //      does not produce a double dC/dt.
  bool shouldCreateTimeDerivative(const std::string& species,
                                  const DiffusionModel<NumericType>& model) {
    (void)model;  // identity not used in this minimal form; MOOSE tracks by physics ptr
    if (timeDerivativeClaimed_.find(species) != timeDerivativeClaimed_.end())
      return false;
    timeDerivativeClaimed_.insert(species);
    return true;
  }

  /// Clear the per-species time-derivative claim set. Called by the engine
  /// at the start of each solve() so composing models can re-claim dC/dt
  /// on the same species across multiple solves on the same physics object.
  /// Mirrors MOOSE PhysicsBase semantics where the gatekeeper is
  /// per-add-kernel (in our case, per-solve), not per-physics-lifetime.
  /// Without this, a second solve() on the same physics would deny the
  /// mass matrix for every species (timeDerivativeClaimed_ still holds
  /// them from the first solve) and silently fall back to identity mass.
  void resetTimeDerivativeClaims() { timeDerivativeClaimed_.clear(); }

  // ---- Per-species BC list (MOOSE MultiSpeciesDiffusionPhysicsBase pattern,
  //      verified: std::vector<std::vector<BoundaryName>> _neumann_boundaries).
  //      Outer key = species, inner = that species' BCs. Replaces the previous
  //      flat std::vector<BCSpec> which scales badly when species have
  //      divergent BC sets.
  struct BCSpec {
    std::string boundary;
    std::string type;  // "neumann", "dirichlet", "segregation", "robin"
    NumericType value;
  };

  void addNeumannBC(const std::string& sp, const std::string& bnd,
                    NumericType flux) {
    bcs_[sp].push_back({bnd, "neumann", flux});
  }

  void addDirichletBC(const std::string& sp, const std::string& bnd,
                      NumericType val) {
    bcs_[sp].push_back({bnd, "dirichlet", val});
  }

  /// Robin (dose-loss) BC: -D dC/dn = h * C. `value` is the transfer
  /// coefficient h [length/time]. Engine adds ∫_Γ h u v to the weak form.
  void addRobinBC(const std::string& sp, const std::string& bnd,
                  NumericType h) {
    bcs_[sp].push_back({bnd, "robin", h});
  }

  /// Per-species BC list. Empty vector if species has no BCs registered.
  const std::vector<BCSpec>& boundaryConditions(const std::string& sp) const {
    static const std::vector<BCSpec> empty;
    auto it = bcs_.find(sp);
    return it == bcs_.end() ? empty : it->second;
  }

  /// All BCs across all species (flat view, for backward compatibility with
  /// engines that loop species-outer). Each entry is tagged with its species.
  struct TaggedBCSpec {
    std::string species;
    BCSpec bc;
  };
  std::vector<TaggedBCSpec> allBoundaryConditions() const {
    std::vector<TaggedBCSpec> out;
    for (const auto& [sp, list] : bcs_)
      for (const auto& bc : list) out.push_back({sp, bc});
    return out;
  }

  void setTemperature(NumericType T) { T_ = T; }
  NumericType temperature() const { return T_; }

private:
  std::vector<std::string> species_;
  std::set<std::string> solverVariables_;        // PhysicsBase::saveSolverVariableName
  std::set<std::string> timeDerivativeClaimed_;  // PhysicsBase::shouldCreateTimeDerivative
  std::vector<std::shared_ptr<DiffusionModel<NumericType>>> models_;
  std::map<std::string, std::vector<BCSpec>> bcs_;  // per-species BC list
  NumericType T_ = NumericType(1273.15);
};

} // namespace viennaps
