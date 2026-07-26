#pragma once

/// KernelTerm — one physics contribution (MOOSE Kernel pattern).
///
/// Each term targets a single species and contributes residual / stiffness /
/// mass. DiffusionModel subclasses (e.g. CddDiffusion) compose many terms.

#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

class KernelTerm {
public:
  virtual ~KernelTerm() = default;

#ifdef VIENNAPS_HAS_MFEM
  virtual void
  assembleResidual(mfem::ParLinearForm &R,
                   const std::map<std::string, mfem::ParGridFunction *> &species,
                   const mfem::ParGridFunction *temp) const {
    (void)R;
    (void)species;
    (void)temp;
  }

  virtual void
  assembleStiffness(mfem::ParBilinearForm &K,
                    const std::map<std::string, mfem::ParGridFunction *> &species,
                    const mfem::ParGridFunction *temp) const {
    (void)K;
    (void)species;
    (void)temp;
  }

  virtual void assembleMass(mfem::ParBilinearForm &M) const { (void)M; }
#endif

  /// Species whose residual/stiffness this term contributes to.
  virtual std::string targetSpecies() const = 0;

  void setName(const std::string &n) { name_ = n; }
  const std::string &getName() const { return name_; }

protected:
  std::string name_;
};

/// Equilibrium (secondary) species derived from a primary after each step.
/// Mirrors MOOSE CoupledBEEquilibriumSub / AuxKernel pattern: not a primary
/// unknown; evaluated after Newton as AuxVariable = f(primary).
template <class NumericType>
class EquilibriumSpeciesAuxKernel {
public:
  EquilibriumSpeciesAuxKernel(std::string primary, std::string equilibrium,
                              NumericType Keq)
      : primary_(std::move(primary)), equilibrium_(std::move(equilibrium)),
        Keq_(Keq) {}

  const std::string &primarySpecies() const { return primary_; }
  const std::string &equilibriumSpecies() const { return equilibrium_; }
  NumericType equilibriumConstant() const { return Keq_; }
  void setEquilibriumConstant(NumericType K) { Keq_ = K; }

  /// Evaluate C_eq = K_eq * C_primary (element-wise / dof-wise).
  void evaluate(const std::vector<NumericType> &primary,
                std::vector<NumericType> &equilibrium) const {
    equilibrium.resize(primary.size());
    for (std::size_t i = 0; i < primary.size(); ++i) {
      equilibrium[i] = Keq_ * primary[i];
    }
  }

#ifdef VIENNAPS_HAS_MFEM
  void evaluate(const mfem::GridFunction &primary,
                mfem::GridFunction &equilibrium) const {
    // Nodal collocation: C_eq_i = Keq * C_primary_i
    for (int i = 0; i < primary.Size(); ++i) {
      equilibrium(i) = static_cast<double>(Keq_) * primary(i);
    }
  }
#endif

private:
  std::string primary_;
  std::string equilibrium_;
  NumericType Keq_ = NumericType(1);
};

} // namespace viennaps
