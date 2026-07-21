#pragma once

#include "MeshAttributes.hpp"

#include <string>
#include <vector>
#include <map>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

/// Abstract base for all diffusion physics models.
///
/// Design reference: MOOSE `Kernel`. Each concrete model operates on one
/// species at a time, located by named GridFunction lookup rather than by
/// manual offset arithmetic into a packed vector. Concrete subclasses
/// declare the species they own via `numSpecies()` / `speciesNames()` and
/// contribute to the global system by overriding the MFEM-gated
/// `assembleStiffness` / `assembleReaction` / `assembleMass` hooks.
template <class NumericType>
class DiffusionModel {
public:
  virtual ~DiffusionModel() = default;

  virtual void setup(const MeshAttributes& attrs, NumericType T) {
    attrs_ = &attrs;
    T_ = T;
  }

#ifdef VIENNAPS_HAS_MFEM
  /// Contribute to stiffness matrix K for this species.
  /// speciesGF: this species' GridFunction (for concentration-dependent D)
  /// allSpecies: map of all species GridFunctions (for coupled models)
  virtual void assembleStiffness(
      mfem::BilinearForm& K,
      const mfem::GridFunction& speciesGF,
      const std::map<std::string, mfem::GridFunction*>& allSpecies,
      const mfem::GridFunction* temp) const {}

  /// Contribute to nonlinear reaction RHS R for this species.
  virtual void assembleReaction(
      mfem::LinearForm& R,
      const mfem::GridFunction& speciesGF,
      const std::map<std::string, mfem::GridFunction*>& allSpecies,
      const mfem::GridFunction* temp) const {}

  /// Contribute to mass matrix M.
  virtual void assembleMass(mfem::BilinearForm& M) const {}
#endif

  virtual int numSpecies() const = 0;
  virtual std::vector<std::string> speciesNames() const = 0;
  virtual std::vector<int> applicableAttributes() const { return {}; }

  void setName(const std::string& n) { name_ = n; }
  const std::string& getName() const { return name_; }

protected:
  const MeshAttributes* attrs_ = nullptr;
  NumericType T_ = NumericType(1273.15);
  std::string name_;
};

} // namespace viennaps
