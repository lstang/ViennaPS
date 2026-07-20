#pragma once

/// DiffusionKernel - Transport kernel for dopant/defect diffusion.
/// Uses MaterialPropertySystem for D(T, conc, stress).
/// For now: explicit update on stub profiles (or MFEM GF in future).
/// Later: assembles MFEM DiffusionIntegrator + mass, driven by SUNDIALS.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <cmath>

namespace viennaps {

template <class NumericType>
class DiffusionKernel : public PhysicsKernel<NumericType> {
public:
  DiffusionKernel(const std::string& species, NumericType temperatureK = 1273.15)
      : species_(species), T_(temperatureK) {
    this->setName("Diffusion(" + species + ")");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }

  void setup() override {
    if (this->material_) {
      D_ = this->material_->getDiffusivity(species_, "Si", T_);
    } else {
      D_ = NumericType(1e-14); // fallback
    }
    std::cout << "[DiffusionKernel] " << species_ << " D=" << D_ << " cm^2/s at T=" << T_ << " K\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;

    // Very simple explicit "diffusion" demo on the stub profile:
    // In practice this will be replaced by MFEM + SUNDIALS time step.
    // Here we just apply a decay/spread factor proportional to D*dt (toy model for visibility).

    auto totalBefore = this->field_->getTotalDose(species_);
    if (totalBefore <= 0) return;

    // Effective reduction factor (simulates some loss + spreading)
    NumericType alpha = D_ * dt * NumericType(1e10); // tune scale for demo
    NumericType factor = NumericType(1) / (NumericType(1) + alpha * 0.1);
    if (factor < NumericType(0.7)) factor = NumericType(0.7);

    this->field_->scaleProfile(species_, factor);

    // TODO (MFEM path): when PhysicsField exposes GridFunction, we will
    // assemble a simple mass-matrix update or call SUNDIALS here.

    std::cout << "[DiffusionKernel] " << species_ << " evolved dt=" << dt
              << "  factor=" << factor
              << "  total " << totalBefore << " -> " << this->field_->getTotalDose(species_) << "\n";
  }

  // Contribute rate to RHS for SUNDIALS: ydot = -alpha(D) * y  (collected by CVODE)
  void addToRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) &&
        myIndex < static_cast<int>(y.size())) {
      // Rate scaled from diffusivity so real D feeds the residual
      NumericType alpha = (D_ > 0) ? (D_ * NumericType(1e10)) : NumericType(0.01);
      if (alpha > NumericType(1)) alpha = NumericType(1); // clamp for stability in demo
      ydot[myIndex] += -alpha * y[myIndex];
    }
  }

private:
  std::string species_;
  NumericType T_ = 1273.15;
  NumericType D_ = 0;
};

} // namespace viennaps
