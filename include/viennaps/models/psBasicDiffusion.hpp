#pragma once

/// BasicDiffusion - Early stub diffusion process model.
/// Demonstrates use of the unified PhysicsField + MaterialPropertySystem.
///
/// In the full implementation this will be replaced / complemented by:
///   - Multiple composable PhysicsKernels (Fermi, ChargedReact, Pair, ... )
///   - SUNDIALS time integration of the coupled nonlinear system
///   - Full MFEM spatial discretization

#include "psProcessModel.hpp"
#include "psDomain.hpp"
#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"

#include <iostream>

namespace viennaps {

template <class NumericType, int D>
class BasicDiffusion : public ProcessModelCPU<NumericType, D> {
public:
  BasicDiffusion() {
    this->setName("BasicDiffusion");
  }

  void setPhysicsField(std::shared_ptr<PhysicsField<NumericType>> field) {
    field_ = field;
  }

  void setMaterialProperties(std::shared_ptr<MaterialPropertySystem<NumericType>> mat) {
    material_ = mat;
  }

  void setTime(NumericType t) { time_ = t; }
  void setTemperature(NumericType T_K) { temperature_ = T_K; }

  void apply(Domain<NumericType, D>& domain) override {
    if (!field_) {
      field_ = domain.getPhysicsField();
    }
    if (!field_) {
      std::cout << "[BasicDiffusion] No physics field in domain, skipping.\n";
      return;
    }

    NumericType T = (temperature_ > 0) ? temperature_ : static_cast<NumericType>(1273.15);

    std::cout << "[BasicDiffusion] Evolving " << time_ << " s at T=" << T << " K (stub)\n";

    if (material_) {
      auto D = material_->getDiffusivity("Boron", "Si", T);
      std::cout << "  Using material diffusivity (B in Si) ~ " << D << " cm^2/s\n";
    }

    // Very crude "diffusion" demo: reduce total dose a little (as if some out-diffusion or clustering)
    auto dose = field_->getTotalDose("Dopant");
    if (dose > 0) {
      NumericType newDose = dose * static_cast<NumericType>(0.98);
      field_->scaleProfile("Dopant", newDose / dose);
      std::cout << "  Dopant dose: " << dose << " -> " << newDose << "\n";
    }

    // Evolve defects a tiny bit too
    auto iDose = field_->getTotalDose("Interstitial");
    if (iDose > 0) {
      field_->scaleProfile("Interstitial", static_cast<NumericType>(0.95));
    }

    field_->evolve(time_);
  }

private:
  std::shared_ptr<PhysicsField<NumericType>> field_;
  std::shared_ptr<MaterialPropertySystem<NumericType>> material_;
  NumericType time_ = 1.0;
  NumericType temperature_ = 0; // K, 0 means use default in material system
};

} // namespace viennaps