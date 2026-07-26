#pragma once

/// MobileImpurity — general mobile impurity with ion-pairing and optional
/// drift (Nernst–Planck): J = −D (∇C + (q/kT) z C E).
///
/// The second template parameter `SpeciesTag` provides compile-time defaults
/// (name, D0, Ea, charge, pairRate).  `GenericImpurityTag` preserves the
/// existing runtime-string behaviour so all existing call sites compile
/// unchanged.  See MobileImpurityTags.hpp for the tag definitions.
///
/// CopperDiffusion is a thin specialization of this class (Phase 4 Task 7–8).

#include "../DiffusionModel.hpp"
#include "../IntrinsicCarrier.hpp"
#include "MobileImpurityTags.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

/// Boltzmann constant (eV/K) — shared default for drift terms.
inline constexpr double kB_eV = 8.617333262145e-5;

template <class NumericType, class SpeciesTag = GenericImpurityTag>
class MobileImpurity : public DiffusionModel<NumericType> {
public:
  explicit MobileImpurity(std::string species = SpeciesTag::name,
                          std::string pairSpecies = "")
      : species_(std::move(species)), pairSpecies_(std::move(pairSpecies)),
        D0_(SpeciesTag::D0), z_(SpeciesTag::charge) {
    this->setName("MobileImpurity(" + species_ + ")");
  }

  void setD0(NumericType D0) { D0_ = D0; }
  void setIonPairing(NumericType beta) { beta_ = beta; }
  void setDopantConcentration(NumericType C) { C_dopant_ = C; }
  void setDopantSpecies(std::string s) { dopantSpecies_ = std::move(s); }

  /// Drift: charge state z and electric field E (V/cm components).
  void setChargeState(NumericType z) { z_ = z; }
  void setElectricField(NumericType Ex, NumericType Ey = 0,
                        NumericType Ez = 0) {
    Ex_ = Ex;
    Ey_ = Ey;
    Ez_ = Ez;
    driftEnabled_ = (Ex != 0 || Ey != 0 || Ez != 0) && (z_ != 0);
  }
  void setDriftEnabled(bool on) { driftEnabled_ = on; }

  /// Pairing reaction mobile + acceptor ⇌ pair (immobile).
  void setPairingRates(NumericType kPair, NumericType kDiss) {
    kPair_ = kPair;
    kDiss_ = kDiss;
  }
  void setAcceptorSpecies(std::string s) { acceptorSpecies_ = std::move(s); }

  NumericType getDiffusivity(NumericType C_dopant, NumericType T) const {
    if (beta_ == NumericType(0))
      return D0_;
    IntrinsicCarrier<NumericType> ic;
    const NumericType ni = ic.ni(T, "Si");
    return D0_ * (NumericType(1) +
                  beta_ * C_dopant / std::max(ni, NumericType(1)));
  }

  /// Nernst factor beta_NP = (q/kT)*z ≈ z/(kB*T) with E in V/cm, T in K.
  NumericType nernstBeta(NumericType T) const {
    if (T <= NumericType(0))
      return NumericType(0);
    return z_ / (static_cast<NumericType>(kB_eV) * T);
  }

  int numSpecies() const override {
    return pairSpecies_.empty() ? 1 : 2;
  }
  std::vector<std::string> speciesNames() const override {
    if (pairSpecies_.empty())
      return {species_};
    return {species_, pairSpecies_};
  }

#ifdef VIENNAPS_HAS_MFEM
  class ImpDCoef : public mfem::Coefficient {
  public:
    ImpDCoef(const MobileImpurity *m, double T, const mfem::GridFunction *dop,
             double Cconst)
        : m_(m), T_(T), dop_(dop), Cconst_(Cconst) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double C = Cconst_;
      if (dop_)
        C = dop_->GetValue(T, ip);
      return static_cast<double>(
          m_->getDiffusivity(static_cast<NumericType>(C),
                             static_cast<NumericType>(T_)));
    }

  private:
    const MobileImpurity *m_;
    double T_;
    const mfem::GridFunction *dop_;
    double Cconst_;
  };

  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    // Immobile pair species: no diffusion/drift.
    if (!pairSpecies_.empty()) {
      auto itP = allSpecies.find(pairSpecies_);
      if (itP != allSpecies.end() && itP->second == &speciesGF) {
        zeroCoef_ = std::make_unique<mfem::ConstantCoefficient>(0.0);
        K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*zeroCoef_));
        return;
      }
    }

    const mfem::GridFunction *dop = nullptr;
    if (!dopantSpecies_.empty()) {
      auto it = allSpecies.find(dopantSpecies_);
      if (it != allSpecies.end())
        dop = it->second;
    }
    dCoef_ = std::make_unique<ImpDCoef>(
        this, static_cast<double>(this->T_), dop,
        static_cast<double>(C_dopant_));
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*dCoef_));

    // Nernst-Planck drift: J = -D(grad C + (q/kT) z C E). Continuity gives
    //   dC/dt = div(D grad C) - div(C v_drift),  v_drift = (D q/kT) z E.
    // Weak form (zero boundary flux):
    //   (dC/dt,v) + (D grad C, grad v) - (C v_drift, grad v) = 0.
    // Using the IBP identity (a.grad u, v) = -(a u, grad v), the drift term
    // -(C v_drift, grad v) = (a grad u, v) requires a = +v_drift = +D bNP E.
    // (Earlier code had a = -D bNP E, which inverted the drift direction;
    //  caught by the copper-drift-direction test.)
    if (driftEnabled_) {
      const double Drep = static_cast<double>(
          getDiffusivity(C_dopant_, this->T_));
      const double bNP =
          static_cast<double>(nernstBeta(this->T_));
      // a = +D * bNP * E  (positive z, positive E => drift along +E).
      mfem::Vector a(3);
      a = 0.0;
      a(0) = Drep * bNP * static_cast<double>(Ex_);
      a(1) = Drep * bNP * static_cast<double>(Ey_);
      if (a.Size() > 2)
        a(2) = Drep * bNP * static_cast<double>(Ez_);
      // ConvectionIntegrator assembles (a . grad u, v); engine form is
      // M du/dt + K u = 0, so this contributes +a to the advection velocity
      // (Cu+ with z>0, E>0 moves +x, as physically expected).
      const int sdim =
          K.FESpace() ? K.FESpace()->GetMesh()->SpaceDimension() : 2;
      mfem::Vector aUse(sdim);
      for (int d = 0; d < sdim; ++d)
        aUse(d) = (d < a.Size()) ? a(d) : 0.0;
      driftVel_ = std::make_unique<mfem::VectorConstantCoefficient>(aUse);
      K.AddDomainIntegrator(new mfem::ConvectionIntegrator(*driftVel_));
    }
    (void)speciesGF;
  }

  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    if (pairSpecies_.empty() || kPair_ == NumericType(0))
      return;
    auto itM = allSpecies.find(species_);
    auto itP = allSpecies.find(pairSpecies_);
    if (itM == allSpecies.end() || itP == allSpecies.end() || !itM->second ||
        !itP->second)
      return;
    const mfem::GridFunction *acc = nullptr;
    if (!acceptorSpecies_.empty()) {
      auto itA = allSpecies.find(acceptorSpecies_);
      if (itA != allSpecies.end())
        acc = itA->second;
    }
    double scale = 0.0;
    if (&speciesGF == itP->second)
      scale = 1.0;
    else if (&speciesGF == itM->second)
      scale = -1.0;
    else
      return;
    pairCoefs_.clear();
    pairCoefs_.push_back(std::make_unique<PairRateCoef>(
        *itM->second, *itP->second, acc, static_cast<double>(kPair_),
        static_cast<double>(kDiss_), scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*pairCoefs_.back()));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

  /// Enable immobile pair species name (numSpecies becomes 2).
  void setPairSpecies(std::string s) { pairSpecies_ = std::move(s); }

protected:
  std::string species_;
  std::string pairSpecies_;
  std::string dopantSpecies_;
  std::string acceptorSpecies_;
  NumericType D0_ = NumericType(1e-10);
  NumericType beta_ = NumericType(0);
  NumericType C_dopant_ = NumericType(0);
  NumericType z_ = NumericType(1);
  NumericType Ex_ = NumericType(0), Ey_ = NumericType(0), Ez_ = NumericType(0);
  bool driftEnabled_ = false;
  NumericType kPair_ = NumericType(0);
  NumericType kDiss_ = NumericType(0);

#ifdef VIENNAPS_HAS_MFEM
  class PairRateCoef : public mfem::Coefficient {
  public:
    PairRateCoef(const mfem::GridFunction &mob, const mfem::GridFunction &pair,
                 const mfem::GridFunction *acc, double kf, double kr,
                 double scale)
        : mob_(&mob), pair_(&pair), acc_(acc), kf_(kf), kr_(kr),
          scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      const double Cm = std::max(0.0, mob_->GetValue(T, ip));
      const double Cp = std::max(0.0, pair_->GetValue(T, ip));
      const double Ca = acc_ ? std::max(0.0, acc_->GetValue(T, ip)) : 1.0;
      return scale_ * (kf_ * Cm * Ca - kr_ * Cp);
    }

  private:
    const mfem::GridFunction *mob_;
    const mfem::GridFunction *pair_;
    const mfem::GridFunction *acc_;
    double kf_, kr_, scale_;
  };
  mutable std::unique_ptr<ImpDCoef> dCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> zeroCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
  mutable std::unique_ptr<mfem::VectorConstantCoefficient> driftVel_;
  mutable std::vector<std::unique_ptr<PairRateCoef>> pairCoefs_;
#endif
};

} // namespace viennaps
