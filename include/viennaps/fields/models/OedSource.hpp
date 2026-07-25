#pragma once

/// OedSource — oxidation-enhanced diffusion interstitial injection.
/// Flux Γ_I = θ * v_ox at the moving Si/SiO2 interface (ADR-0004).

#include "../DiffusionModel.hpp"
#include "../DiffusionPhysics.hpp"

#include <string>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

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
  void setBoundary(std::string b) { boundary_ = std::move(b); }

  NumericType injectionFlux() const { return theta_ * v_ox_; }

  /// Register as Neumann flux BC on the interface boundary attribute.
  void registerWith(DiffusionPhysics<NumericType> &physics) const {
    physics.addNeumannBC(I_, boundary_, injectionFlux());
  }

  void applyInjection(std::vector<NumericType> &C_I, int interfaceBin,
                      NumericType dx, NumericType dt) const {
    if (interfaceBin < 0 ||
        interfaceBin >= static_cast<int>(C_I.size()) || dx <= 0)
      return;
    C_I[static_cast<std::size_t>(interfaceBin)] +=
        injectionFlux() * dt / dx;
  }

  /// Advance interface attribute tag for subdomain-relabeling idiom (A).
  /// Elements with material `fromAttr` whose "oxidation progress" exceeds
  /// threshold flip to `toAttr` (Si → SiO2).
  static int relabelOxidized(
#ifdef VIENNAPS_HAS_MFEM
      mfem::Mesh &mesh,
#else
      int /*dummy*/,
#endif
      int fromAttr, int toAttr, NumericType progress,
      NumericType threshold) {
    int flipped = 0;
#ifdef VIENNAPS_HAS_MFEM
    if (progress < threshold)
      return 0;
    for (int e = 0; e < mesh.GetNE(); ++e) {
      if (mesh.GetAttribute(e) == fromAttr) {
        // Simple global flip when progress crosses threshold (unit test).
        // Production: flip only elements adjacent to the interface.
        mesh.SetAttribute(e, toAttr);
        ++flipped;
      }
    }
#else
    (void)fromAttr;
    (void)toAttr;
    (void)progress;
    (void)threshold;
#endif
    return flipped;
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override { return {I_}; }

private:
  std::string I_;
  std::string boundary_ = "all";
  NumericType theta_ = NumericType(0.01);
  NumericType v_ox_ = NumericType(0);
};

} // namespace viennaps
