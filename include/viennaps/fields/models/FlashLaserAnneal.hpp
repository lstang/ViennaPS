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
      mfem::BilinearForm &K, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(
        static_cast<double>(alpha_));
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    if (q_ == NumericType(0))
      return;
    srcCoef_ = std::make_unique<mfem::ConstantCoefficient>(
        static_cast<double>(q_));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*srcCoef_));
  }
  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  NumericType alpha_ = NumericType(0.8); // cm^2/s order Si thermal
  NumericType q_ = NumericType(0);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> srcCoef_;
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
  explicit MeltDiffusion(std::string species = "Boron")
      : species_(std::move(species)) {
    this->setName("MeltDiffusion(" + species_ + ")");
  }

  void setSolidD(NumericType D) { Ds_ = D; }
  void setLiquidD(NumericType D) { Dl_ = D; }

  NumericType getDiffusivity(NumericType phi) const {
    // phi=0 solid, phi=1 liquid
    return Ds_ * (NumericType(1) - phi) + Dl_ * phi;
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

private:
  std::string species_;
  NumericType Ds_ = NumericType(1e-12);
  NumericType Dl_ = NumericType(1e-4);
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
