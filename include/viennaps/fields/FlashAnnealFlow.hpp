#pragma once

/// FlashAnnealFlow — end-to-end flash/laser anneal orchestration on
/// DiffusionEngine (GAP_ANALYSIS §4.10): surface heat pulse → melting
/// phase field (latent heat ρ·L·∂φ/∂t coupled) → melt-enhanced dopant
/// diffusion + solidification trapping + SPER crystallinity. All species
/// and models run in ONE multi-species solve per segment; between
/// segments the MeltFraction copy is refreshed for the ∂φ/∂t terms.

#include "DiffusionEngine.hpp"
#include "DiffusionPhysics.hpp"
#include "models/FlashLaserAnneal.hpp"

#include <memory>

namespace viennaps {

template <class NumericType> class FlashAnnealFlow {
public:
  FlashAnnealFlow() { buildPhysics(); }

  // ---- Configuration ------------------------------------------------------
  void setPulse(NumericType Tpeak, NumericType duration) {
    Tpeak_ = Tpeak;
    duration_ = duration;
  }
  void setLaserAbsorption(NumericType alpha) { laserAlpha_ = alpha; }
  void setLatentHeat(NumericType rhoL) { heat_->setLatentHeat(rhoL); }
  void setDopantDiffusivities(NumericType Ds, NumericType Dl) {
    meltDiff_->setSolidD(Ds);
    meltDiff_->setLiquidD(Dl);
  }
  void setTrappingStrength(NumericType r) { trap_->setTrappingStrength(r); }
  void setMeltParameters(NumericType Tm, NumericType L, NumericType kappa,
                         NumericType lambda) {
    meltFem_->setMeltingPoint(Tm);
    meltFem_->setMobility(L);
    meltFem_->setGradientEnergy(kappa);
    meltFem_->setCoupling(lambda);
  }

  /// Seed the Temperature species with a Beer's-law surface pulse:
  /// T(z) = T0 + I0 * exp(-alpha * z) along the dof axis, projected with
  /// projectIntegralPreserving and rescaled to the analytical mean
  /// (assumes a unit-area domain, matching the test meshes).
  void seedLaserPulse(DiffusionEngine<NumericType, 2> &engine, NumericType T0,
                      NumericType I0) const {
    const int n = 64;
    std::vector<NumericType> samples(static_cast<std::size_t>(n),
                                     NumericType(0));
    for (int i = 0; i < n; ++i) {
      const double z = static_cast<double>(i) / static_cast<double>(n - 1);
      samples[static_cast<std::size_t>(i)] =
          T0 + I0 * std::exp(-static_cast<double>(laserAlpha_) * z);
    }
    // Analytical integral of T(z) over the unit depth (per unit area).
    const double alpha = static_cast<double>(laserAlpha_);
    const NumericType dose = static_cast<NumericType>(
        static_cast<double>(T0) + I0 * (1.0 - std::exp(-alpha)) / alpha);
    engine.projectIntegralPreserving("Temperature", samples, dose);
  }

  /// Run `nSteps` segments of `dt` each. Between segments, the
  /// previous-φ copy is refreshed and handed to the latent-heat and
  /// trapping models via setPreviousPhi/setDt.
  void apply(DiffusionEngine<NumericType, 2> &engine, NumericType dt,
             int nSteps) {
    // Ensure all 4 species are initialized in the engine if not already set.
    if (engine.getSolution("Temperature").Size() == 0)
      engine.initializeSpecies("Temperature", NumericType(300));
    if (engine.getSolution("MeltFraction").Size() == 0)
      engine.initializeSpecies("MeltFraction", NumericType(0));
    if (engine.getSolution("Crystallinity").Size() == 0)
      engine.initializeSpecies("Crystallinity", NumericType(1));
    if (engine.getSolution("Dopant").Size() == 0)
      engine.initializeSpecies("Dopant", NumericType(0));

    NumericType t = NumericType(0);
    for (int s = 0; s < nSteps; ++s) {
#ifdef VIENNAPS_HAS_MFEM
      if (!prevPhi_) {
        prevPhi_ = std::make_unique<mfem::ParGridFunction>(
            engine.getSolution("MeltFraction"));
      } else {
        *prevPhi_ = engine.getSolution("MeltFraction");
      }
      heat_->setPreviousPhi(prevPhi_.get());
      heat_->setDt(dt);
      trap_->setPreviousPhi(prevPhi_.get());
      trap_->setDt(dt);
#endif
      engine.solve(t, t + dt, dt);
      t += dt;
    }
  }

  DiffusionPhysics<NumericType> &physics() { return physics_; }

private:
  void buildPhysics() {
    physics_.addSpecies("Temperature");
    physics_.addSpecies("MeltFraction");
    physics_.addSpecies("Crystallinity");
    physics_.addSpecies("Dopant");

    heat_ = std::make_shared<HeatTransfer<NumericType>>();
    heat_->setThermalDiffusivity(NumericType(0.8));
    meltFem_ = std::make_shared<MeltingPhaseFieldFEM<NumericType>>();
    meltFem_->setMeltingPoint(NumericType(1687));
    meltFem_->setTemperatureSpecies("Temperature");
    meltDiff_ = std::make_shared<MeltDiffusion<NumericType>>("Dopant");
    cryst_ = std::make_shared<CrystallinityPhaseFieldFEM<NumericType>>();
    cryst_->setTemperatureSpecies("Temperature");
    trap_ = std::make_shared<SolidificationTrapping<NumericType>>("Dopant");

    physics_.addModel(heat_);
    physics_.addModel(meltFem_);
    physics_.addModel(meltDiff_);
    physics_.addModel(cryst_);
    physics_.addModel(trap_);
    physics_.setTemperature(NumericType(300));
  }

  DiffusionPhysics<NumericType> physics_;
  std::shared_ptr<HeatTransfer<NumericType>> heat_;
  std::shared_ptr<MeltingPhaseFieldFEM<NumericType>> meltFem_;
  std::shared_ptr<MeltDiffusion<NumericType>> meltDiff_;
  std::shared_ptr<CrystallinityPhaseFieldFEM<NumericType>> cryst_;
  std::shared_ptr<SolidificationTrapping<NumericType>> trap_;
#ifdef VIENNAPS_HAS_MFEM
  std::unique_ptr<mfem::ParGridFunction> prevPhi_;
#endif
  NumericType Tpeak_ = NumericType(1500);
  NumericType duration_ = NumericType(1e-3);
  NumericType laserAlpha_ = NumericType(1e4);
};

} // namespace viennaps
