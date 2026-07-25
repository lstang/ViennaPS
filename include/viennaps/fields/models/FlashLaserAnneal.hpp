#pragma once

/// Flash/Laser anneal models (Phase 9 skeleton).

#include "../DiffusionModel.hpp"
#include "ConstantDiffusion.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class HeatTransfer {
public:
  void setThermalDiffusivity(NumericType alpha) { alpha_ = alpha; }
  void setSource(NumericType q) { q_ = q; }

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

private:
  NumericType alpha_ = NumericType(0.8); // cm^2/s order Si thermal
  NumericType q_ = NumericType(0);
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

  /// Effective anneal: raise T for pulse, return time-averaged D scale.
  NumericType effectiveTimeAt(NumericType Tref) const {
    (void)Tref;
    return duration_;
  }

  NumericType peakTemperature() const { return Tpeak_; }
  NumericType duration() const { return duration_; }

  HeatTransfer<NumericType> heat;
  LaserIntensity<NumericType> laser;
  MeltingPhaseField<NumericType> melt;
  MeltDiffusion<NumericType> meltDiff;

private:
  NumericType Tpeak_ = NumericType(1500);
  NumericType duration_ = NumericType(1e-3);
};

} // namespace viennaps
