#pragma once

/// PhysicsFieldAdapter - Thin coupling between PhysicsField and existing Oxidation.
///
/// Does NOT modify psOxidation internals. Provides:
///  - Pre-oxidation hooks: dopant-enhanced rate factor, hydrostatic stress readout
///  - Post-oxidation OED defect injection into unified fields
///  - Field-only APIs for unit tests without a full Domain/level-set setup

#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace viennaps {

// Forward declare to avoid pulling ViennaLS headers into multiphysics-only tests.
template <class NumericType, int D>
class Domain;

template <class NumericType, int D>
class PhysicsFieldAdapter {
public:
  PhysicsFieldAdapter(std::shared_ptr<PhysicsField<NumericType>> field,
                      std::shared_ptr<MaterialPropertySystem<NumericType>> mat = nullptr)
      : field_(std::move(field)), material_(std::move(mat)) {}

  void setPhysicsField(std::shared_ptr<PhysicsField<NumericType>> f) { field_ = std::move(f); }
  void setMaterialProperties(std::shared_ptr<MaterialPropertySystem<NumericType>> m) {
    material_ = std::move(m);
  }

  void setOEDDosePerStep(NumericType dose) { oedDose_ = dose; }
  void setDopantEnhancementScale(NumericType s) { dopantScale_ = s; }

  /// Dopant-dependent oxidation enhancement factor (linear in surface-ish dopant dose).
  /// Returns ~1.0 with no dopant; grows gently with Dopant total dose.
  NumericType getDopantEnhancedOxidationFactor() const {
    if (!field_) return NumericType(1);
    NumericType dose = field_->getTotalDose("Dopant");
    // Reference dose 1e13 cm^-2 -> +10% max enhancement in demo units
    NumericType enh = NumericType(1) + dopantScale_ * dose / NumericType(1e14);
    return std::min(NumericType(2), std::max(NumericType(1), enh));
  }

  /// Hydrostatic stress available for oxidation stress coupling (Pa-like units).
  NumericType getHydrostaticStressForOxidation() const {
    if (!field_) return NumericType(0);
    return field_->getTotalDose("HydrostaticStress");
  }

  /// Inject OED-style interstitial/vacancy defects after oxidation growth.
  /// oxideThicknessDelta: relative growth proxy (dimensionless or um); scales OED dose.
  void injectOEDDefects(NumericType oxideThicknessDelta = NumericType(1)) {
    if (!field_) return;
    field_->addSpecies("Interstitial");
    field_->addSpecies("Vacancy");
    field_->addSpecies("OxidationDefects");

    NumericType dose = oedDose_ * std::max(NumericType(0), oxideThicknessDelta);
    // Dopant can slightly increase OED (TED-like coupling)
    NumericType enh = getDopantEnhancedOxidationFactor();
    dose *= enh;

    // Dose-conserving inject via size-1 total (not hardcoded profile length)
    field_->addDose("Interstitial", dose);
    // Vacancies slightly less than I for net interstitial injection (classic OED)
    field_->addDose("Vacancy", dose * NumericType(0.7));
    field_->addDose("OxidationDefects", dose);

    lastOEDDose_ = dose;
    std::cout << "[PhysicsFieldAdapter] OED inject dose=" << dose
              << " (oxideDelta=" << oxideThicknessDelta
              << ", dopantEnh=" << enh << ")\n";
  }

  /// Called before existing Oxidation model: log coupling state and cache factors.
  void applyToOxidation(Domain<NumericType, D>& /*domain*/) {
    if (!field_) return;
    preOxDopantFactor_ = getDopantEnhancedOxidationFactor();
    preOxStress_ = getHydrostaticStressForOxidation();
    std::cout << "[PhysicsFieldAdapter] Pre-oxidation sync:\n"
              << "  Dopant total dose = " << field_->getTotalDose("Dopant") << "\n"
              << "  Interstitials     = " << field_->getTotalDose("Interstitial") << "\n"
              << "  Dopant ox factor  = " << preOxDopantFactor_ << "\n"
              << "  HydrostaticStress = " << preOxStress_ << "\n";
  }

  /// Called after existing Oxidation: inject OED defects into unified field.
  void updateFromOxidation(Domain<NumericType, D>& /*domain*/) {
    if (!field_) return;
    // Use a default growth proxy; real path would read oxide thickness from domain
    injectOEDDefects(NumericType(1));
    std::cout << "[PhysicsFieldAdapter] Post-oxidation: OED defects injected into unified field.\n";
  }

  /// Field-only apply/update for tests without Domain construction.
  void applyToOxidationFieldOnly() {
    if (!field_) return;
    preOxDopantFactor_ = getDopantEnhancedOxidationFactor();
    preOxStress_ = getHydrostaticStressForOxidation();
    std::cout << "[PhysicsFieldAdapter] Pre-oxidation (field-only) dopantFactor="
              << preOxDopantFactor_ << " stress=" << preOxStress_ << "\n";
  }

  void updateFromOxidationFieldOnly(NumericType oxideThicknessDelta = NumericType(1)) {
    injectOEDDefects(oxideThicknessDelta);
  }

  void syncStressToOxidation() {
    preOxStress_ = getHydrostaticStressForOxidation();
    std::cout << "[PhysicsFieldAdapter] syncStressToOxidation stress=" << preOxStress_ << "\n";
  }

  void syncStressFromOxidation(NumericType residualStress = NumericType(0)) {
    if (!field_) return;
    if (residualStress != 0) {
      // Add residual onto existing hydrostatic channel (dose-conserving)
      field_->addDose("HydrostaticStress", residualStress);
    }
    std::cout << "[PhysicsFieldAdapter] syncStressFromOxidation residual="
              << residualStress << "\n";
  }

  NumericType getLastOEDDose() const { return lastOEDDose_; }
  NumericType getCachedDopantFactor() const { return preOxDopantFactor_; }
  NumericType getCachedStress() const { return preOxStress_; }

private:
  std::shared_ptr<PhysicsField<NumericType>> field_;
  std::shared_ptr<MaterialPropertySystem<NumericType>> material_;
  NumericType oedDose_ = static_cast<NumericType>(1e11);
  NumericType dopantScale_ = static_cast<NumericType>(0.1);
  NumericType lastOEDDose_ = 0;
  NumericType preOxDopantFactor_ = 1;
  NumericType preOxStress_ = 0;
};

} // namespace viennaps
