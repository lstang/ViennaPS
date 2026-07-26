#pragma once

/// ChargedEquilibriumDiffusion — charge-state equilibrium diffusivity.
/// D_eff = Σ_z D^z f^z  with f^z ∝ exp(−z · η), η = (E_F − E_i)/kT ≈ ln(n/ni)
/// (SProcess 4.194 / same math as ChargedFermiDiffusion).

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ChargedEquilibriumDiffusion : public DiffusionModel<NumericType> {
public:
  explicit ChargedEquilibriumDiffusion(std::string species = "Boron")
      : species_(std::move(species)) {
    this->setName("ChargedEquilibriumDiffusion(" + species_ + ")");
    z_ = {-1, 0, 1};
    // Unequal D^z so extrinsic shift changes D_eff (same defaults as ChargedFermi).
    relDz_ = {NumericType(0.1), NumericType(1), NumericType(5)};
    setD0(D0_);
  }

  void setD0(NumericType D0) {
    D0_ = D0;
    Dz_.resize(relDz_.size());
    for (std::size_t i = 0; i < relDz_.size(); ++i)
      Dz_[i] = relDz_[i] * D0_;
  }

  /// Absolute D^z per charge state (optional).
  void setChargeStates(std::vector<int> z, std::vector<NumericType> Dz) {
    z_ = std::move(z);
    Dz_ = std::move(Dz);
    if (z_.size() != Dz_.size()) {
      const std::size_t n = std::min(z_.size(), Dz_.size());
      z_.resize(n);
      Dz_.resize(n);
    }
  }

  void setNi(NumericType ni) { ni_ = ni; }

  /// Legacy API (ignored for formula; kept so existing setters compile).
  void setAlpha(NumericType) {}

  NumericType getDiffusivity(NumericType C, NumericType T) const {
    NumericType ni = ni_;
    if (ni <= NumericType(0)) {
      IntrinsicCarrier<NumericType> ic;
      ni = ic.ni(T, "Si");
    }
    const NumericType n = std::max(C, ni);
    const NumericType eta =
        (ni > 0 && n > 0) ? std::log(n / ni) : NumericType(0);

    NumericType sumW = NumericType(0);
    NumericType sumD = NumericType(0);
    for (std::size_t i = 0; i < z_.size(); ++i) {
      const NumericType w =
          std::exp(-static_cast<NumericType>(z_[i]) * eta);
      sumW += w;
      sumD += Dz_[i] * w;
    }
    if (sumW <= NumericType(0))
      return D0_;
    return sumD / sumW;
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  class EqDCoef : public mfem::Coefficient {
  public:
    EqDCoef(const ChargedEquilibriumDiffusion *m, double T)
        : m_(m), T_(T), conc_(nullptr) {}
    void SetField(const mfem::ParGridFunction *c) { conc_ = c; }
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double C = 0.0;
      if (conc_)
        C = conc_->GetValue(T, ip);
      return static_cast<double>(
          m_->getDiffusivity(static_cast<NumericType>(C),
                             static_cast<NumericType>(T_)));
    }

  private:
    const ChargedEquilibriumDiffusion *m_;
    double T_;
    const mfem::ParGridFunction *conc_;
  };

  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> & /*allSpecies*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    stiffCoef_ =
        std::make_unique<EqDCoef>(this, static_cast<double>(this->T_));
    stiffCoef_->SetField(&speciesGF);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-14);
  NumericType ni_ = NumericType(0);
  std::vector<int> z_;
  std::vector<NumericType> Dz_;
  std::vector<NumericType> relDz_;
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<EqDCoef> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
