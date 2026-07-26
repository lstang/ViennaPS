#pragma once

/// Flash/Laser anneal models (Phase 9 skeleton).

#include "../DiffusionModel.hpp"
#include "ConstantDiffusion.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class HeatTransfer : public DiffusionModel<NumericType> {
public:
  HeatTransfer() { this->setName("HeatTransfer"); }
  void setThermalDiffusivity(NumericType alpha) { alpha_ = alpha; }
  void setSource(NumericType q) { q_ = q; }
  NumericType thermalDiffusivity() const { return alpha_; }

  /// Latent-heat coupling (SProcess eq. 213 term ρ·L·∂φ/∂t). The melt
  /// fraction φ is read from a registered species GridFunction named
  /// `meltSpecies_` (default "MeltFraction", matching MeltDiffusion's
  /// convention). `previousPhi_` is updated by the orchestrator between
  /// steps via `setPreviousPhi()`. When φ is rising (melting), the term
  /// is negative (absorbs heat); when falling (solidifying), positive.
  void setMeltSpecies(std::string s) { meltSpecies_ = std::move(s); }
  void setLatentHeat(NumericType rhoL) { rhoL_ = rhoL; }
  void setPreviousPhi(const mfem::ParGridFunction *prev) { previousPhi_ = prev; }

  /// 1D explicit heat step on T profile.
  void step(std::vector<NumericType> &T, NumericType dx, NumericType dt) const {
    if (T.size() < 3)
      return;
    const NumericType a = alpha_ * dt / (dx * dx);
    std::vector<NumericType> n = T;
    for (std::size_t i = 1; i + 1 < T.size(); ++i) {
      n[i] = T[i] + a * (T[i - 1] - NumericType(2) * T[i] + T[i + 1]) +
             q_ * dt;
    }
    T.swap(n);
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {"Temperature"};
  }

#ifdef VIENNAPS_HAS_MFEM
  /// FEM heat equation stiffness: α ∇T·∇v (+ optional volumetric source in R).
  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::ParGridFunction *> & /*allSpecies*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(
        static_cast<double>(alpha_));
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }
  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    // Volumetric heat source Q (constant).
    if (q_ != NumericType(0)) {
      srcCoef_ = std::make_unique<mfem::ConstantCoefficient>(
          static_cast<double>(q_));
      R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*srcCoef_));
    }
    // Latent-heat coupling ρ·L·∂φ/∂t (SProcess eq. 213). Read current φ
    // from the registered MeltFraction species; ∂φ/∂t ≈ (φ_curr - φ_prev)/dt.
    // The previous-step φ pointer is set by the orchestrator via
    // setPreviousPhi(); dt is set via setDt(). When φ is rising (melting),
    // ∂φ/∂t > 0 and the term is negative (absorbs latent heat).
    if (rhoL_ != NumericType(0) && dt_ > NumericType(0) && previousPhi_) {
      auto it = allSpecies.find(meltSpecies_);
      if (it != allSpecies.end() && it->second) {
        latentCoef_ = std::make_unique<LatentHeatCoef>(
            it->second, previousPhi_, static_cast<double>(rhoL_),
            static_cast<double>(dt_));
        R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*latentCoef_));
      }
    }
  }
  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
  /// Set the current time step (for ∂φ/∂t finite difference).
  void setDt(NumericType dt) { dt_ = dt; }
#endif

private:
  NumericType alpha_ = NumericType(0.8); // cm^2/s order Si thermal
  NumericType q_ = NumericType(0);
  NumericType rhoL_ = NumericType(0); // latent heat coefficient ρ·L
  NumericType dt_ = NumericType(0);   // current dt for ∂φ/∂t
  std::string meltSpecies_ = "MeltFraction";
  const mfem::ParGridFunction *previousPhi_ = nullptr;
#ifdef VIENNAPS_HAS_MFEM
  class LatentHeatCoef : public mfem::Coefficient {
  public:
    LatentHeatCoef(const mfem::ParGridFunction *phiCurr,
                   const mfem::ParGridFunction *phiPrev, double rhoL, double dt)
        : phiCurr_(phiCurr), phiPrev_(phiPrev), rhoL_(rhoL), dt_(dt) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      if (!phiCurr_ || !phiPrev_ || dt_ <= 0.0)
        return 0.0;
      const double cur = phiCurr_->GetValue(T, ip);
      const double prev = phiPrev_->GetValue(T, ip);
      const double dphiDt = (cur - prev) / dt_;
      // ∂φ/∂t > 0 (melting) -> term is negative (absorbs heat).
      return -rhoL_ * dphiDt;
    }

  private:
    const mfem::ParGridFunction *phiCurr_;
    const mfem::ParGridFunction *phiPrev_;
    double rhoL_, dt_;
  };
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> srcCoef_;
  mutable std::unique_ptr<LatentHeatCoef> latentCoef_;
#endif
};

template <class NumericType>
class LaserIntensity {
public:
  void setPeak(NumericType I0) { I0_ = I0; }
  void setAbsorption(NumericType alpha) { alpha_ = alpha; }

  /// Beer's law: I(z) = I0 * exp(-alpha*z)
  NumericType intensity(NumericType z) const {
    return I0_ * std::exp(-alpha_ * std::max(z, NumericType(0)));
  }

private:
  NumericType I0_ = NumericType(1e5);
  NumericType alpha_ = NumericType(1e4); // 1/cm
};

/// Allen–Cahn phase-field term: dφ/dt = -L * (df/dφ - κ ∇²φ)
template <class NumericType>
class AllenCahnTerm {
public:
  void setMobility(NumericType L) { L_ = L; }
  void setGradientEnergy(NumericType kappa) { kappa_ = kappa; }

  /// Double-well df/dφ = 2φ(1-φ)(1-2φ)
  static NumericType df_dphi(NumericType phi) {
    return NumericType(2) * phi * (NumericType(1) - phi) *
           (NumericType(1) - NumericType(2) * phi);
  }

  NumericType mobility() const { return L_; }
  NumericType kappa() const { return kappa_; }

private:
  NumericType L_ = NumericType(1);
  NumericType kappa_ = NumericType(1e-12);
};

template <class NumericType>
class MeltingPhaseField {
public:
  void setMeltingPoint(NumericType Tm) { Tm_ = Tm; }

  /// Equilibrium order parameter: 0 solid, 1 liquid.
  NumericType phiEq(NumericType T) const {
    return (T >= Tm_) ? NumericType(1) : NumericType(0);
  }

  void relax(std::vector<NumericType> &phi, const std::vector<NumericType> &T,
             NumericType rate, NumericType dt) const {
    const std::size_t n = std::min(phi.size(), T.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType target = phiEq(T[i]);
      phi[i] += rate * (target - phi[i]) * dt;
      phi[i] = std::min(NumericType(1), std::max(NumericType(0), phi[i]));
    }
  }

private:
  NumericType Tm_ = NumericType(1687); // Si melting K
};

template <class NumericType>
class MeltDiffusion : public DiffusionModel<NumericType> {
public:
  explicit MeltDiffusion(std::string species = "Boron",
                         std::string meltFrac = "MeltFraction")
      : species_(std::move(species)), melt_(std::move(meltFrac)) {
    this->setName("MeltDiffusion(" + species_ + ")");
  }

  void setSolidD(NumericType D) { Ds_ = D; }
  void setLiquidD(NumericType D) { Dl_ = D; }
  void setMeltSpecies(std::string m) { melt_ = std::move(m); }

  NumericType getDiffusivity(NumericType phi) const {
    return Ds_ * (NumericType(1) - phi) + Dl_ * phi;
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  class MeltDCoef : public mfem::Coefficient {
  public:
    MeltDCoef(const MeltDiffusion *m, const mfem::ParGridFunction *phi)
        : m_(m), phi_(phi) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double p = 0.0;
      if (phi_)
        p = std::min(1.0, std::max(0.0, phi_->GetValue(T, ip)));
      return static_cast<double>(
          m_->getDiffusivity(static_cast<NumericType>(p)));
    }

  private:
    const MeltDiffusion *m_;
    const mfem::GridFunction *phi_;
  };

  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    const mfem::ParGridFunction *phi = nullptr;
    auto it = allSpecies.find(melt_);
    if (it != allSpecies.end())
      phi = it->second;
    meltCoef_ = std::make_unique<MeltDCoef>(this, phi);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*meltCoef_));
  }
  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  std::string melt_;
  NumericType Ds_ = NumericType(1e-12);
  NumericType Dl_ = NumericType(1e-4);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<MeltDCoef> meltCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

/// FEM Allen-Cahn melting phase field (SProcess eq. 22575):
/// dφ/dt = -L * (df/dφ - κ ∇²φ),  f(φ,T) = (φ²-1)²/4 - λ(T-Tm)φ
/// df/dφ = φ³ - φ - λ(T-Tm)
/// Assembled as: M dφ/dt + L*κ*K*φ = -L*R(φ,T)
/// where K is the diffusion stiffness and R is the reaction (bulk driving force).
template <class NumericType>
class MeltingPhaseFieldFEM : public DiffusionModel<NumericType> {
public:
  MeltingPhaseFieldFEM() {
    this->setName("MeltingPhaseFieldFEM");
  }
  void setMobility(NumericType L) { L_ = L; }
  void setGradientEnergy(NumericType kappa) { kappa_ = kappa; }
  void setMeltingPoint(NumericType Tm) { Tm_ = Tm; }
  void setCoupling(NumericType lambda) { lambda_ = lambda; }
  void setTemperatureSpecies(std::string s) { tempSpecies_ = std::move(s); }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {"MeltFraction"};
  }

#ifdef VIENNAPS_HAS_MFEM
  /// Stiffness: L*κ * ∇φ·∇v (gradient energy term).
  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::ParGridFunction *> & /*allSpecies*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    const double lk = static_cast<double>(L_ * kappa_);
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(lk);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  /// Reaction: -L * (φ³ - φ - λ(T-Tm)) as a QP-local source.
  /// Reads φ from the species GF and T from the temperature species.
  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    const mfem::ParGridFunction *Tgf = nullptr;
    if (!tempSpecies_.empty()) {
      auto it = allSpecies.find(tempSpecies_);
      if (it != allSpecies.end())
        Tgf = it->second;
    }
    reactCoef_ = std::make_unique<AllenCahnReactionCoef>(
        &speciesGF, Tgf, static_cast<double>(L_),
        static_cast<double>(lambda_), static_cast<double>(Tm_));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*reactCoef_));
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  NumericType L_ = NumericType(1.0);
  NumericType kappa_ = NumericType(1e-4);
  NumericType Tm_ = NumericType(1687);
  NumericType lambda_ = NumericType(1.0);
  std::string tempSpecies_ = "Temperature";
#ifdef VIENNAPS_HAS_MFEM
  class AllenCahnReactionCoef : public mfem::Coefficient {
  public:
    AllenCahnReactionCoef(const mfem::ParGridFunction *phi,
                          const mfem::ParGridFunction *T, double L, double lambda,
                          double Tm)
        : phi_(phi), T_(T), L_(L), lambda_(lambda), Tm_(Tm) {}
    double Eval(mfem::ElementTransformation &tr,
                const mfem::IntegrationPoint &ip) override {
      double p = 0.0;
      if (phi_)
        p = std::min(1.0, std::max(-1.0, phi_->GetValue(tr, ip)));
      double T = Tm_;
      if (T_)
        T = T_->GetValue(tr, ip);
      // df/dφ = φ³ - φ - λ(T-Tm). Reaction = -L * df/dφ.
      return -L_ * (p * p * p - p - lambda_ * (T - Tm_));
    }

  private:
    const mfem::ParGridFunction *phi_;
    const mfem::ParGridFunction *T_;
    double L_, lambda_, Tm_;
  };
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<AllenCahnReactionCoef> reactCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

/// FEM Allen-Cahn crystallinity phase field (SProcess eq. 23069):
/// df/dφ = φ³ - φ - λ*v_SPER(T)*φ, driven by SPER velocity.
template <class NumericType>
class CrystallinityPhaseFieldFEM : public DiffusionModel<NumericType> {
public:
  CrystallinityPhaseFieldFEM() {
    this->setName("CrystallinityPhaseFieldFEM");
  }
  void setMobility(NumericType L) { L_ = L; }
  void setGradientEnergy(NumericType kappa) { kappa_ = kappa; }
  void setSperVelocity(NumericType v0, NumericType Ea) {
    v0_ = v0;
    Ea_ = Ea;
  }
  void setCoupling(NumericType lambda) { lambda_ = lambda; }
  void setTemperatureSpecies(std::string s) { tempSpecies_ = std::move(s); }

  NumericType sperVelocity(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return v0_;
    return v0_ * std::exp(-Ea_ / (kB * T));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {"Crystallinity"};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::ParGridFunction *> & /*allSpecies*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    const double lk = static_cast<double>(L_ * kappa_);
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(lk);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    const mfem::ParGridFunction *Tgf = nullptr;
    if (!tempSpecies_.empty()) {
      auto it = allSpecies.find(tempSpecies_);
      if (it != allSpecies.end())
        Tgf = it->second;
    }
    reactCoef_ = std::make_unique<CrystReactionCoef>(
        &speciesGF, Tgf, static_cast<double>(L_),
        static_cast<double>(lambda_), static_cast<double>(v0_),
        static_cast<double>(Ea_));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*reactCoef_));
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  NumericType L_ = NumericType(1.0);
  NumericType kappa_ = NumericType(1e-4);
  NumericType lambda_ = NumericType(1.0);
  NumericType v0_ = NumericType(1e-6);
  NumericType Ea_ = NumericType(2.7);
  std::string tempSpecies_ = "Temperature";
#ifdef VIENNAPS_HAS_MFEM
  class CrystReactionCoef : public mfem::Coefficient {
  public:
    CrystReactionCoef(const mfem::ParGridFunction *phi,
                      const mfem::ParGridFunction *T, double L, double lambda,
                      double v0, double Ea)
        : phi_(phi), T_(T), L_(L), lambda_(lambda), v0_(v0), Ea_(Ea) {}
    double Eval(mfem::ElementTransformation &tr,
                const mfem::IntegrationPoint &ip) override {
      double p = 0.0;
      if (phi_)
        p = std::min(1.0, std::max(-1.0, phi_->GetValue(tr, ip)));
      double T = 1000.0;
      if (T_)
        T = T_->GetValue(tr, ip);
      const double kB = 8.617333262145e-5;
      const double vSper =
          (T > 0) ? v0_ * std::exp(-Ea_ / (kB * T)) : v0_;
      // df/dφ = φ³ - φ - λ*v_SPER(T)*φ. Reaction = -L * df/dφ.
      return -L_ * (p * p * p - p - lambda_ * vSper * p);
    }

  private:
    const mfem::ParGridFunction *phi_;
    const mfem::ParGridFunction *T_;
    double L_, lambda_, v0_, Ea_;
  };
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<CrystReactionCoef> reactCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

template <class NumericType>
class FlashLaserAnneal {
public:
  void setPulse(NumericType Tpeak, NumericType duration) {
    Tpeak_ = Tpeak;
    duration_ = duration;
  }

  NumericType effectiveTimeAt(NumericType Tref) const {
    (void)Tref;
    return duration_;
  }

  NumericType peakTemperature() const { return Tpeak_; }
  NumericType duration() const { return duration_; }

  /// Orchestrate heat → melt phase → effective melt diffusivity history.
  struct Result {
    std::vector<NumericType> T;
    std::vector<NumericType> phi;
    std::vector<NumericType> Deff;
  };

  Result runPulse(std::vector<NumericType> T0, NumericType dx,
                  NumericType dt, int nSteps) {
    Result r;
    r.T = std::move(T0);
    r.phi.assign(r.T.size(), NumericType(0));
    r.Deff.reserve(static_cast<std::size_t>(nSteps));
    // Deposit laser energy near surface (index 0).
    if (!r.T.empty())
      r.T[0] += Tpeak_ * NumericType(0.1);
    for (int s = 0; s < nSteps; ++s) {
      heat.step(r.T, dx, dt);
      melt.relax(r.phi, r.T, NumericType(5), dt);
      NumericType Dmean = NumericType(0);
      for (auto p : r.phi)
        Dmean += meltDiff.getDiffusivity(p);
      Dmean /= static_cast<NumericType>(std::max<std::size_t>(1, r.phi.size()));
      r.Deff.push_back(Dmean);
    }
    return r;
  }

  HeatTransfer<NumericType> heat;
  LaserIntensity<NumericType> laser;
  MeltingPhaseField<NumericType> melt;
  MeltDiffusion<NumericType> meltDiff;

private:
  NumericType Tpeak_ = NumericType(1500);
  NumericType duration_ = NumericType(1e-3);
};

} // namespace viennaps
