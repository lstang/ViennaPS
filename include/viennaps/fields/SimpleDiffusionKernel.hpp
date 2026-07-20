#pragma once

/// SimpleDiffusionKernel - A minimal concrete PhysicsKernel example.
/// For demonstration and early testing only.
///
/// It uses the MaterialPropertySystem to look up diffusivity and then
/// performs a trivial explicit update on the stub profiles in PhysicsField.
/// 
/// Real kernels will:
///   - Assemble MFEM forms (diffusion operator + reactions)
///   - Be registered with a time integrator (SUNDIALS)
///   - Support multiple species + coupling (I/V recombination, clustering, OED sources, stress)

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <vector>
#include <cmath>

namespace viennaps {

template <class NumericType>
class SimpleDiffusionKernel : public PhysicsKernel<NumericType> {
public:
  SimpleDiffusionKernel(const std::string& species = "Dopant")
      : species_(species) {
    this->setName("SimpleDiffusion(" + species + ")");
  }

  void setup() override {
    if (this->material_) {
      // Pre-query a reference D at a typical temp (will be per-step in real)
      refD_ = this->material_->getDiffusivity(species_, "Si", 1273.15);
    } else {
      refD_ = NumericType(1e-15); // fallback
    }
    std::cout << "[SimpleDiffusionKernel] setup for " << species_
              << "  refD=" << refD_ << std::endl;
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;

    // Extremely simplified explicit "diffusion" on the stub profile:
    // We just scale the profile slightly according to D*dt (toy model)
    // In reality we would do finite-volume / FEM step here or via MFEM operators.

    const auto& prof = this->field_->getProfile(species_);
    if (prof.empty()) return;

    // Compute a fake "loss" or spread factor
    NumericType factor = NumericType(1) - std::min(NumericType(0.05),
                           refD_ * dt * NumericType(1e12)); // tune for demo visibility

    if (factor < 0.5) factor = 0.5;

    // Apply to the profile via the field helper
    // (we reuse scale because the stub profile is internal)
    // A better kernel would directly mutate a copy or use field API extensions.

    // For demo we call a simple rescale on total + profile
    auto totalBefore = this->field_->getTotalDose(species_);
    this->field_->scaleProfile(species_, factor);

    std::cout << "[SimpleDiffusionKernel] " << species_
              << " evolved dt=" << dt << "  factor=" << factor
              << "  total " << totalBefore << " -> " << this->field_->getTotalDose(species_) << std::endl;
  }

private:
  std::string species_;
  NumericType refD_ = NumericType(1e-15);
};

} // namespace viennaps
