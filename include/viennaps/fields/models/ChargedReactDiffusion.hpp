#pragma once

/// ChargedReactDiffusion — I+V recombination with charge-enhanced rate.
/// k_eff = k0 * (1 + gamma * n/ni), n ≈ max(C_dopant, ni)
///
/// FEM: overrides assembleReaction with ProductCoef scaled by k_eff.
/// When a dopant GridFunction is registered (setDopantSpecies), n is
/// evaluated at quadrature points; otherwise uses setDopantConcentration.

#include "ReactDiffusion.hpp"
#include "../IntrinsicCarrier.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class ChargedReactDiffusion : public ReactDiffusion<NumericType> {
public:
  ChargedReactDiffusion(std::string interstitial = "Interstitial",
                        std::string vacancy = "Vacancy")
      : ReactDiffusion<NumericType>(std::move(interstitial),
                                    std::move(vacancy)) {
    this->setName("ChargedReactDiffusion");
  }

  void setChargeEnhancement(NumericType gamma, NumericType ni) {
    gamma_ = gamma;
    ni_ = ni;
  }

  void setDopantConcentration(NumericType C) { C_dopant_ = C; }

  /// Optional coupled dopant field for QP-local n(x).
  void setDopantSpecies(std::string name) { dopantSpecies_ = std::move(name); }

  NumericType effectiveRate(NumericType T,
                            NumericType C_local = NumericType(-1)) const {
    NumericType ni = ni_;
    if (ni <= NumericType(0)) {
      ni = IntrinsicCarrier<NumericType>{}.ni(T, "Si");
    }
    const NumericType C =
        (C_local >= NumericType(0)) ? C_local : C_dopant_;
    const NumericType n = std::max(C, ni);
    return this->recombinationRate() *
           (NumericType(1) + gamma_ * n / std::max(ni, NumericType(1)));
  }

  /// Host reaction step using charge-enhanced k (uniform C_dopant_).
  void applyChargedReactionStep(std::vector<NumericType> &I,
                                std::vector<NumericType> &V, NumericType dt,
                                NumericType T) const {
    const NumericType kEff = effectiveRate(T);
    const std::size_t n = std::min(I.size(), V.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType r = kEff * I[i] * V[i] * dt;
      I[i] = std::max(NumericType(0), I[i] - r);
      V[i] = std::max(NumericType(0), V[i] - r);
    }
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    auto itI = allSpecies.find(this->interstitialName());
    auto itV = allSpecies.find(this->vacancyName());
    if (itI == allSpecies.end() || itV == allSpecies.end() || !itI->second ||
        !itV->second)
      return;

    const mfem::ParGridFunction *dopantGF = nullptr;
    if (!dopantSpecies_.empty()) {
      auto itD = allSpecies.find(dopantSpecies_);
      if (itD != allSpecies.end())
        dopantGF = itD->second;
    }

    chargedCoefs_.clear();
    chargedCoefs_.push_back(std::make_unique<ChargedProductCoef>(
        *itI->second, *itV->second, dopantGF, this,
        static_cast<double>(this->T_)));
    R.AddDomainIntegrator(
        new mfem::DomainLFIntegrator(*chargedCoefs_.back()));
    (void)speciesGF;
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  /// R = -k_eff(n) * C_I * C_V at each QP.
  class ChargedProductCoef : public mfem::Coefficient {
  public:
    ChargedProductCoef(const mfem::ParGridFunction &I, const mfem::ParGridFunction &V,
                       const mfem::ParGridFunction *dopant,
                       const ChargedReactDiffusion *model, double T)
        : I_(&I), V_(&V), dopant_(dopant), model_(model), T_(T) {}

    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double cI = I_->GetValue(T, ip);
      const double cV = V_->GetValue(T, ip);
      double C = -1.0;
      if (dopant_)
        C = dopant_->GetValue(T, ip);
      const double kEff = static_cast<double>(model_->effectiveRate(
          static_cast<NumericType>(T_), static_cast<NumericType>(C)));
      return -kEff * cI * cV;
    }

  private:
    const mfem::ParGridFunction *I_;
    const mfem::ParGridFunction *V_;
    const mfem::ParGridFunction *dopant_;
    const ChargedReactDiffusion *model_;
    double T_;
  };
  mutable std::vector<std::unique_ptr<ChargedProductCoef>> chargedCoefs_;
#endif

  NumericType gamma_ = NumericType(1);
  NumericType ni_ = NumericType(0);
  NumericType C_dopant_ = NumericType(0);
  std::string dopantSpecies_;
};

} // namespace viennaps
