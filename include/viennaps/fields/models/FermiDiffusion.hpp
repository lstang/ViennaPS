#pragma once

/// FermiDiffusion — extrinsic (Fermi) diffusion model.
///
/// D_eff(C,T) = D_i(T) * (1 + alpha * n/ni)  with n ≈ max(C, ni) for n-type.
///
/// Jacobian strategy (b1): supplies analytic dD/dC via
/// assembleStiffnessJacobian + FermiDdCCoef (MOOSE MatDiffusion +
/// DerivativeMaterialInterface pattern). See Phase 2 Task 2 plan.

#include "../DiffusivityMaterial.hpp"
#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class FermiDiffusion : public DiffusionModel<NumericType> {
public:
  explicit FermiDiffusion(const std::string &species = "Boron") {
    this->setName("FermiDiffusion(" + species + ")");
    species_ = species;
  }

  void setSpecies(const std::string &s) { species_ = s; }

  /// Intrinsic diffusivity D_i and enhancement factor alpha.
  void setDiffusivity(NumericType D_i, NumericType alpha = NumericType(1)) {
    D_i_ = D_i;
    alpha_ = alpha;
  }

  void setIntrinsicCarrierConcentration(NumericType ni) { ni_ = ni; }

  /// Use IntrinsicCarrier to set ni from material/temperature.
  void setNiFromMaterial(const std::string &material = "Si") {
    IntrinsicCarrier<NumericType> ic;
    ni_ = ic.ni(this->T_, material);
  }

  void setup(const MeshAttributes &attrs, NumericType T) override {
    DiffusionModel<NumericType>::setup(attrs, T);
    // Refresh ni at the simulation temperature if still at default sentinel.
    if (ni_ <= NumericType(0)) {
      setNiFromMaterial("Si");
    }
  }

  NumericType getDiffusivity(NumericType C, NumericType T) const {
    // T is accepted for API symmetry; D_i_ is treated as already evaluated
    // at T (caller sets D_i via setDiffusivity). ni_ is fixed after setup.
    (void)T;
    return FermiDiffusivity<NumericType>::D(C, D_i_, alpha_, ni_);
  }

  NumericType evalDdC(NumericType C, NumericType T) const {
    (void)T;
    return FermiDiffusivity<NumericType>::dDdC(C, D_i_, alpha_, ni_);
  }

  NumericType getDi() const { return D_i_; }
  NumericType getAlpha() const { return alpha_; }
  NumericType getNi() const { return ni_; }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }
  std::vector<int> applicableAttributes() const override { return attrs_; }
  void setApplicableAttributes(std::vector<int> a) { attrs_ = std::move(a); }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> & /*allSpecies*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    // Owned coefficient: DiffusionIntegrator holds Coefficient& across
    // Assemble() (same lifetime pattern as ConstantDiffusion).
    stiffCoef_ = std::make_unique<FermiDCoef>(
        static_cast<double>(D_i_), static_cast<double>(alpha_),
        static_cast<double>(ni_), static_cast<double>(this->T_));
    stiffCoef_->SetConcentrationField(&speciesGF);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleStiffnessJacobian(
      mfem::MixedBilinearForm &dKdC,
      const mfem::ParGridFunction &speciesGF) const override {
    ddcCoef_ = std::make_unique<FermiDdCCoef>(static_cast<double>(D_i_),
                                              static_cast<double>(alpha_),
                                              static_cast<double>(ni_));
    ddcCoef_->SetConcentrationField(&speciesGF);
    dKdC.AddDomainIntegrator(new mfem::MixedGradGradIntegrator(*ddcCoef_));
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  NumericType D_i_ = NumericType(1e-13);
  NumericType alpha_ = NumericType(1);
  NumericType ni_ = NumericType(0); // 0 => refresh from IntrinsicCarrier in setup
  std::vector<int> attrs_;
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<FermiDCoef> stiffCoef_;
  mutable std::unique_ptr<FermiDdCCoef> ddcCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
