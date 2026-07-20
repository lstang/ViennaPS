#pragma once

/// StressKernel / ViscoelasticStressKernel
/// Computes stress fields (hydrostatic, etc.) from:
/// - Growth / mismatch (oxidation, deposition, etc.)
/// - Elastic response
/// - Viscoelastic relaxation (Maxwell model)
///
/// Outputs to PhysicsField as "HydrostaticStress", "GrowthStress" etc.
/// Feeds back to diffusion (stress-dependent D) and oxidation adapter.
/// Uses MFEM when available for proper elasticity solve (future full impl).
/// For now: average + relaxation on stub + MFEM GridFunction storage.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <cmath>

namespace viennaps {

template <class NumericType>
class ViscoelasticStressKernel : public PhysicsKernel<NumericType> {
public:
  ViscoelasticStressKernel(NumericType temperatureK = 1273.15)
      : T_(temperatureK) {
    this->setName("ViscoelasticStress");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }

  void setup() override {
    if (this->material_) {
      E_ = this->material_->getYoungModulus("Si", T_);
      nu_ = this->material_->getPoissonRatio("Si", T_);
      eta_ = this->material_->getViscosity("SiO2", T_);
      alpha_ = this->material_->getCTE("Si", T_);
    } else {
      E_ = 130.0; nu_ = 0.28; eta_ = 1e13; alpha_ = 2.6e-6;
    }
    // lame or bulk for hydrostatic
    bulk_ = E_ / (3.0 * (1.0 - 2.0 * nu_));
    std::cout << "[ViscoelasticStressKernel] E=" << E_ << " GPa, nu=" << nu_
              << ", eta=" << eta_ << ", bulk=" << bulk_ << " at T=" << T_ << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;

    // Get or create stress fields
    this->field_->addSpecies("HydrostaticStress");
    this->field_->addSpecies("GrowthStress");

    // Growth stress
    NumericType growthStress = this->material_ ? this->material_->getProperty("Si", "GrowthStress", T_) : 300.0e6;
    NumericType currentGrowth = this->field_->getTotalDose("GrowthStress");
    NumericType newGrowth = currentGrowth + growthStress * 0.01 * dt;

    std::vector<NumericType> gs(1, newGrowth);
    this->field_->injectImplantProfile("GrowthStress", gs);

    // Viscoelastic relaxation on top
    NumericType G = E_ / (2.0 * (1.0 + nu_));
    NumericType tau = (eta_ > 0) ? (eta_ / G) * 1e-9 : 1e10;
    NumericType currentStress = this->field_->getTotalDose("HydrostaticStress");

    NumericType relaxed = currentStress * std::exp( -dt / std::max(tau, NumericType(1e-6)) );

    NumericType dT = 0;
    NumericType thermal = bulk_ * 3.0 * alpha_ * dT * 1e9;

    NumericType totalStress = relaxed + thermal + newGrowth * 0.1;

    std::vector<NumericType> hs(1, totalStress);
    this->field_->injectImplantProfile("HydrostaticStress", hs);

    std::cout << "[ViscoelasticStressKernel] stress evolved dt=" << dt
              << "  hydrostatic=" << totalStress / 1e6 << " MPa (relaxed)\n";
  }

  // Contribute relaxation/growth rate to stress "state"
  virtual void addToRHS(NumericType t, const std::vector<NumericType>& y, std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) && myIndex < static_cast<int>(y.size())) {
      NumericType G = E_ / (2.0 * (1.0 + nu_));
      NumericType tau = (eta_ > 0) ? (eta_ / G) * 1e-9 : 1e10;
      NumericType rate = - y[myIndex] / std::max(tau, NumericType(1e-6));
      ydot[myIndex] += rate;
    }
  }

private:
  NumericType T_ = 1273.15;
  NumericType E_ = 130.0, nu_ = 0.28, eta_ = 1e13, alpha_ = 2.6e-6;
  NumericType bulk_ = 0;
};

} // namespace viennaps
