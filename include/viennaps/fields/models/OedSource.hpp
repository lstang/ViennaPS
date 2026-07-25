#pragma once

/// OedSource — oxidation-enhanced diffusion interstitial injection.
/// Flux Γ_I = θ * v_ox at the moving Si/SiO2 interface (ADR-0004).

#include "../DiffusionModel.hpp"
#include "../DiffusionPhysics.hpp"
#include "../MovingMeshHandler.hpp"

#include <map>
#include <memory>
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
      NumericType threshold, bool interfaceAdjacentOnly = true) {
#ifdef VIENNAPS_HAS_MFEM
    // Delegate to MovingMeshHandler (ADR-0004); default interface-adjacent.
    auto r = MovingMeshHandler::relabelAttributes(
        mesh, fromAttr, toAttr, static_cast<double>(progress),
        static_cast<double>(threshold), interfaceAdjacentOnly);
    return r.elementsRelabeled;
#else
    (void)fromAttr;
    (void)toAttr;
    (void)progress;
    (void)threshold;
    (void)interfaceAdjacentOnly;
    return 0;
#endif
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override { return {I_}; }

#ifdef VIENNAPS_HAS_MFEM
  /// Volumetric OED injection proxy when interface BC attrs are unavailable:
  /// R += Γ_I as a uniform DomainLFIntegrator (unit test / 1D-like path).
  /// Production path prefers registerWith() Neumann on the moving boundary.
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    const double flux = static_cast<double>(injectionFlux());
    if (flux == 0.0)
      return;
    srcCoef_ = std::make_unique<mfem::ConstantCoefficient>(flux);
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*srcCoef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string I_;
  std::string boundary_ = "all";
  NumericType theta_ = NumericType(0.01);
  NumericType v_ox_ = NumericType(0);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> srcCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
