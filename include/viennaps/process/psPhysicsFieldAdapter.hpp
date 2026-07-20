#pragma once

/// PhysicsFieldAdapter - Thin coupling layer between the new unified PhysicsField
/// (dopants, defects, stress) and existing high-quality models (especially Oxidation).
///
/// Design principle: NEVER modify the internals of psOxidation / existing models.
/// Instead, the adapter:
///   - Extracts current field state (dose, profiles, stress) before oxidation step
///   - Injects OED defects, dopant-modified rates, stress info into the existing model (via public params when available)
///   - Pulls back newly generated interstitials / vacancies / stress updates after oxidation
///
/// This preserves all the validated LOCOS, mask bending, viscous flow, stress, and Deal-Grove logic.

#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"
#include "psDomain.hpp"

namespace viennaps {

template <class NumericType, int D>
class PhysicsFieldAdapter {
public:
  PhysicsFieldAdapter(std::shared_ptr<PhysicsField<NumericType>> field,
                      std::shared_ptr<MaterialPropertySystem<NumericType>> mat = nullptr)
      : field_(field), material_(mat) {}

  void setPhysicsField(std::shared_ptr<PhysicsField<NumericType>> f) { field_ = f; }
  void setMaterialProperties(std::shared_ptr<MaterialPropertySystem<NumericType>> m) { material_ = m; }

  // Called before invoking the existing Oxidation model.
  // Here we would configure dopant-dependent oxidation rate, supply stress, etc.
  void applyToOxidation(Domain<NumericType, D>& domain) {
    if (!field_) return;

    std::cout << "[PhysicsFieldAdapter] Pre-oxidation sync:\n"
              << "  Dopant total dose = " << field_->getTotalDose("Dopant") << "\n"
              << "  Interstitials     = " << field_->getTotalDose("Interstitial") << "\n";

    // Placeholder for real coupling:
    // - Read current dopant profile near the Si/SiO2 interface
    // - Modify oxidation parameters (e.g. linear/parabolic rates via dopant effect)
    // - Pass hydrostatic stress from field to oxidation stress solver
    //
    // Example (when Oxidation exposes setters):
    //   auto& oxParams = ... ;
    //   NumericType surfB = field_->getConcentration("Boron", 0.0);
    //   oxParams.setDopantEnhancedFactor( computeEnhancement(surfB) );
  }

  // Called after the existing Oxidation step.
  // Pull newly injected defects (OED) back into the unified field.
  void updateFromOxidation(Domain<NumericType, D>& domain) {
    if (!field_) return;

    // In real implementation the oxidation model (or its velocity/stress output)
    // would tell us how many interstitials were injected.
    // For now we demonstrate by adding a small OED source.
    field_->addSpecies("Interstitial");
    field_->addSpecies("OxidationDefects");

    // Simulate a small OED injection (orders of magnitude will be calibrated)
    std::vector<NumericType> oed(32, static_cast<NumericType>(1e11));
    field_->injectImplantProfile("Interstitial", oed);
    field_->injectImplantProfile("OxidationDefects", oed);

    std::cout << "[PhysicsFieldAdapter] Post-oxidation: OED defects injected into unified field.\n";
  }

  // Future: stress feedback round-trip
  void syncStressToOxidation() { /* TODO */ }
  void syncStressFromOxidation() { /* TODO */ }

private:
  std::shared_ptr<PhysicsField<NumericType>> field_;
  std::shared_ptr<MaterialPropertySystem<NumericType>> material_;
};

} // namespace viennaps