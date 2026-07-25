#pragma once

/// ChargedFermiDiffusion — charge-state-decomposed Fermi diffusion.
///
/// D_eff = sum_z D^z * f^z(n, p, T)
/// where charge-state fractions are
///   f^z ∝ exp(-z * (E_F - E_i) / kT)
/// normalized so sum_z f^z = 1.
///
/// E_F is estimated from charge neutrality for n-type doping:
///   n ≈ max(C, ni), p = ni^2/n, and (E_F - E_i)/kT ≈ ln(n/ni).
///
/// Default charge states z ∈ {-1, 0, +1} with equal intrinsic D^z = D0/3
/// so that D_eff → D0 at C = ni (intrinsic limit).

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ChargedFermiDiffusion : public DiffusionModel<NumericType> {
public:
  explicit ChargedFermiDiffusion(const std::string &species = "Boron") {
    this->setName("ChargedFermiDiffusion(" + species + ")");
    species_ = species;
    // Default charge states z ∈ {-1, 0, +1} with unequal D^z so that
    // shifting f^z with Fermi level changes D_eff (equal D^z would make
    // D independent of charge-state fractions). Relative weights applied
    // in setDiffusivity via D0_.
    z_ = {-1, 0, 1};
    relDz_ = {NumericType(0.1), NumericType(1), NumericType(5)};
    customDz_ = false;
    setDiffusivity(D0_);
  }

  void setSpecies(const std::string &s) { species_ = s; }

  /// Set reference diffusivity D0. Absolute D^z = relDz_i * D0 unless the
  /// user supplied absolute D^z via setChargeStates.
  void setDiffusivity(NumericType D0) {
    D0_ = D0;
    if (!customDz_) {
      Dz_.resize(relDz_.size());
      for (std::size_t i = 0; i < relDz_.size(); ++i) {
        Dz_[i] = relDz_[i] * D0_;
      }
    }
  }

  /// Absolute charge-state diffusivities (marks Dz as user-custom).
  void setChargeStates(std::vector<int> z, std::vector<NumericType> Dz) {
    z_ = std::move(z);
    Dz_ = std::move(Dz);
    customDz_ = true;
    if (z_.size() != Dz_.size()) {
      const std::size_t n = std::min(z_.size(), Dz_.size());
      z_.resize(n);
      Dz_.resize(n);
    }
  }

  void setIntrinsicCarrierConcentration(NumericType ni) { ni_ = ni; }

  void setup(const MeshAttributes &attrs, NumericType T) override {
    DiffusionModel<NumericType>::setup(attrs, T);
    if (ni_ <= NumericType(0)) {
      IntrinsicCarrier<NumericType> ic;
      ni_ = ic.ni(this->T_, "Si");
    }
  }

  /// Effective diffusivity at concentration C and temperature T.
  NumericType getDiffusivity(NumericType C, NumericType T) const {
    const NumericType ni = (ni_ > 0) ? ni_ : IntrinsicCarrier<NumericType>{}.ni(T, "Si");
    const NumericType n = std::max(C, ni);
    // (E_F - E_i)/kT ≈ ln(n/ni) for n-type non-degenerate approx.
    const NumericType eta =
        (ni > 0 && n > 0) ? std::log(n / ni) : NumericType(0);

    NumericType sumW = NumericType(0);
    NumericType sumD = NumericType(0);
    for (std::size_t i = 0; i < z_.size(); ++i) {
      // f^z ∝ exp(-z * eta); z is integer charge state.
      const NumericType w = std::exp(-static_cast<NumericType>(z_[i]) * eta);
      sumW += w;
      sumD += Dz_[i] * w;
    }
    if (sumW <= NumericType(0)) {
      return D0_;
    }
    return sumD / sumW;
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }
  std::vector<int> applicableAttributes() const override { return attrs_; }
  void setApplicableAttributes(std::vector<int> a) { attrs_ = std::move(a); }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    // Picard: lag D from a representative field value (mean concentration).
    // Full quadrature-dependent D is available via getDiffusivity; a dedicated
    // Coefficient can replace this when engine Newton is wired (Phase 2+).
    const double Cmean = speciesGF.Size() > 0
                             ? speciesGF.Sum() / speciesGF.Size()
                             : 0.0;
    const double Deff =
        static_cast<double>(getDiffusivity(static_cast<NumericType>(Cmean),
                                           this->T_));
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(Deff);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-13);
  NumericType ni_ = NumericType(0);
  std::vector<int> z_;
  std::vector<NumericType> Dz_;
  std::vector<NumericType> relDz_; ///< relative weights; Dz = relDz * D0
  bool customDz_ = false;
  std::vector<int> attrs_;
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
