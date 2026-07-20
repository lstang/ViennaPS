#pragma once

/// PairDiffusionKernel - dopant-interstitial pair diffusion model.
/// Effective transport: D_pair ~ D_I * (C_I / C_Ieq) * f(C_dopant)
/// Consumes a fraction of free I while enhancing dopant mobility (TED-like).

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

namespace viennaps {

template <class NumericType>
class PairDiffusionKernel : public PhysicsKernel<NumericType> {
public:
  PairDiffusionKernel(const std::string& dopant = "Dopant",
                      NumericType temperatureK = 1273.15)
      : dopant_(dopant), T_(temperatureK) {
    this->setName("PairDiffusion(" + dopant + ")");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }

  void setup() override {
    if (this->material_) {
      DI_ = this->material_->getDiffusivity("Interstitial", "Si", T_);
      CeqI_ = this->material_->getEquilibriumConcentration("Interstitial", "Si", T_);
      if (CeqI_ <= 0) CeqI_ = NumericType(1e15);
    } else {
      DI_ = NumericType(1e-10);
      CeqI_ = NumericType(1e15);
    }
    std::cout << "[PairDiffusionKernel] DI=" << DI_ << " CeqI=" << CeqI_ << " at T=" << T_
              << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;
    this->field_->addSpecies(dopant_);
    this->field_->addSpecies("Interstitial");

    NumericType B = this->field_->getTotalDose(dopant_);
    NumericType I = this->field_->getTotalDose("Interstitial");
    if (B <= 0 || I <= 0) return;

    // Pair fraction ~ min(B,I) * supersaturation factor
    NumericType super = I / std::max(CeqI_, NumericType(1));
    NumericType pairRate = std::min(B, I) * NumericType(0.02) * std::min(super, NumericType(10)) * dt;

    // Dopant redistributes slightly (TED proxy: dose preserved, mild surface pile-up via scale)
    this->field_->scaleProfile(dopant_, NumericType(1)); // refresh
    // Consume some free interstitials via pair formation / kick-out
    if (I > 0) {
      NumericType keep = std::max(NumericType(0.7), (I - pairRate * NumericType(0.1)) / I);
      this->field_->scaleProfile("Interstitial", keep);
    }

    this->field_->addSpecies("PairBI");
    this->field_->injectImplantProfile("PairBI", std::vector<NumericType>(1, pairRate));

    std::cout << "[PairDiffusionKernel] pairRate=" << pairRate << " dt=" << dt
              << " B=" << this->field_->getTotalDose(dopant_)
              << " I=" << this->field_->getTotalDose("Interstitial") << "\n";
  }

  void addToRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) &&
        myIndex < static_cast<int>(y.size())) {
      ydot[myIndex] += -NumericType(0.03) * y[myIndex];
    }
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto offB = this->field_->getSpeciesOffset(dopant_);
    auto offI = this->field_->getSpeciesOffset("Interstitial");
    if (offB == static_cast<std::size_t>(-1) || offI == static_cast<std::size_t>(-1))
      return;
    auto n = this->field_->getProfileSize();
    if (offB + n > y.size() || offI + n > y.size()) return;

    NumericType k = NumericType(0.02);
    for (std::size_t i = 0; i < n; ++i) {
      NumericType Bi = y[offB + i];
      NumericType Ii = y[offI + i];
      NumericType rate = k * Bi * Ii / (CeqI_ + NumericType(1));
      // pair diffusion: dopant spreads with I, mild I sink
      if (i > 0 && i + 1 < n) {
        ydot[offB + i] += NumericType(0.05) * (y[offB + i - 1] - NumericType(2) * Bi + y[offB + i + 1]) *
                          (NumericType(1) + Ii / (CeqI_ + NumericType(1)));
      }
      ydot[offI + i] += -rate * NumericType(0.1);
    }
  }

private:
  std::string dopant_;
  NumericType T_ = 1273.15;
  NumericType DI_ = 0;
  NumericType CeqI_ = 1e15;
};

} // namespace viennaps
