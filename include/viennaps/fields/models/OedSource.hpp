#pragma once

/// OedSource — oxidation-enhanced diffusion interstitial injection.
/// Flux of I at the moving Si/SiO2 interface: Γ_I = θ * v_ox
/// where v_ox = dx_ox/dt and θ is the injection efficiency.

#include "../DiffusionModel.hpp"

#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class OedSource : public DiffusionModel<NumericType> {
public:
  explicit OedSource(std::string interstitial = "Interstitial")
      : I_(std::move(interstitial)) {
    this->setName("OedSource");
  }

  void setInjectionEfficiency(NumericType theta) { theta_ = theta; }
  void setOxidationRate(NumericType v_ox) { v_ox_ = v_ox; }

  /// Interstitial injection flux [cm^-2 s^-1] at the interface.
  NumericType injectionFlux() const { return theta_ * v_ox_; }

  /// Inject into a 1D host profile near the interface (test helper).
  /// Adds flux * dt / dx to the interface bin.
  void applyInjection(std::vector<NumericType> &C_I, int interfaceBin,
                      NumericType dx, NumericType dt) const {
    if (interfaceBin < 0 ||
        interfaceBin >= static_cast<int>(C_I.size()) || dx <= 0)
      return;
    C_I[static_cast<std::size_t>(interfaceBin)] +=
        injectionFlux() * dt / dx;
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override { return {I_}; }

private:
  std::string I_;
  NumericType theta_ = NumericType(0.01); // typical OED efficiency
  NumericType v_ox_ = NumericType(0);     // cm/s
};

} // namespace viennaps
