#pragma once

/// ChargedPairDiffusion — pair diffusion with Fermi-level dependent D_pair.
/// D_eff = D_pair * (C_I / C_I_eq) * (1 + alpha * n/ni)
///
/// FEM: QP-local coefficient; n from the dopant GridFunction at each IP.

#include "PairDiffusion.hpp"
#include "../IntrinsicCarrier.hpp"
#include "../PointDefectEquilibrium.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>

namespace viennaps {

template <class NumericType>
class ChargedPairDiffusion : public PairDiffusion<NumericType> {
public:
  ChargedPairDiffusion(std::string dopant = "Boron",
                       std::string interstitial = "Interstitial")
      : PairDiffusion<NumericType>(std::move(dopant),
                                   std::move(interstitial)) {
    this->setName("ChargedPairDiffusion");
  }

  void setFermiEnhancement(NumericType alpha, NumericType ni) {
    alpha_ = alpha;
    ni_ = ni;
  }

  NumericType getDiffusivity(NumericType C_I, NumericType C_dopant,
                             NumericType T) const {
    NumericType ni = ni_;
    if (ni <= NumericType(0)) {
      ni = IntrinsicCarrier<NumericType>{}.ni(T, "Si");
    }
    const NumericType n = std::max(C_dopant, ni);
    const NumericType fermiFactor =
        NumericType(1) + alpha_ * n / std::max(ni, NumericType(1));
    return PairDiffusion<NumericType>::getDiffusivity(C_I, T) * fermiFactor;
  }

#ifdef VIENNAPS_HAS_MFEM
  class ChargedPairDCoef : public mfem::Coefficient {
  public:
    ChargedPairDCoef(const mfem::GridFunction *CI,
                     const mfem::GridFunction *Cdop, double Dpair, double Ceq,
                     double alpha, double ni)
        : CI_(CI), Cdop_(Cdop), Dpair_(Dpair), Ceq_(std::max(Ceq, 1.0)),
          alpha_(alpha), ni_(std::max(ni, 1.0)) {}

    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double cI = 0.0;
      if (CI_)
        cI = std::max(0.0, CI_->GetValue(T, ip));
      double C = 0.0;
      if (Cdop_)
        C = std::max(0.0, Cdop_->GetValue(T, ip));
      const double n = std::max(C, ni_);
      const double fermi = 1.0 + alpha_ * n / ni_;
      return Dpair_ * (cI / Ceq_) * fermi;
    }

  private:
    const mfem::GridFunction *CI_;
    const mfem::GridFunction *Cdop_;
    double Dpair_, Ceq_, alpha_, ni_;
  };

  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    const mfem::GridFunction *CI = nullptr;
    auto it = allSpecies.find(this->interstitialName());
    if (it != allSpecies.end())
      CI = it->second;

    double Ceq = static_cast<double>(this->ciEqOverride());
    if (Ceq <= 0.0) {
      Ceq = static_cast<double>(
          PointDefectEquilibrium<NumericType>{}.C_I_eq(this->T_, "Si"));
    }

    double ni = static_cast<double>(ni_);
    if (ni <= 0.0)
      ni = static_cast<double>(
          IntrinsicCarrier<NumericType>{}.ni(this->T_, "Si"));

    chargedPairCoef_ = std::make_unique<ChargedPairDCoef>(
        CI, &speciesGF, static_cast<double>(this->pairDiffusivity()), Ceq,
        static_cast<double>(alpha_), ni);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*chargedPairCoef_));

    if (this->supgEnabled() && CI && CI->Size() > 0) {
      // Reuse parent SUPG path via temporary PairDiffusion-style artificial D.
      const double cMean = CI->Sum() / CI->Size();
      const double Deff =
          static_cast<double>(this->pairDiffusivity()) * cMean / Ceq;
      const double vmag = Deff / std::max(static_cast<double>(this->hmin_), 1e-12);
      SupgAdvectionTerm supg(this->dopantName(), vmag, 0.0,
                             static_cast<double>(this->hmin_));
      supg.assembleStiffness(K, allSpecies, nullptr);
    }
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    PairDiffusion<NumericType>::assembleMass(M);
  }
#endif

private:
  NumericType alpha_ = NumericType(1);
  NumericType ni_ = NumericType(0);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<ChargedPairDCoef> chargedPairCoef_;
#endif
};

} // namespace viennaps
