#pragma once

/// ChargedReactKernel - charged point-defect reaction model (I-V recomb + charge).
/// Captures basic charged-defect kinetics used in ATHENA/SProcess:
///   I + V -> 0  (recombination)
///   rate enhanced by local supersaturation product C_I * C_V

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

namespace viennaps {

template <class NumericType>
class ChargedReactKernel : public PhysicsKernel<NumericType> {
public:
  explicit ChargedReactKernel(NumericType temperatureK = 1273.15)
      : T_(temperatureK) {
    this->setName("ChargedReact");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }
  void setRecombinationRate(NumericType k) { kIV_ = k; }

  void setup() override {
    if (this->material_) {
      CeqI_ = this->material_->getEquilibriumConcentration("Interstitial", "Si", T_);
      CeqV_ = this->material_->getEquilibriumConcentration("Vacancy", "Si", T_);
    }
    if (CeqI_ <= 0) CeqI_ = NumericType(1e15);
    if (CeqV_ <= 0) CeqV_ = NumericType(1e15);
    std::cout << "[ChargedReactKernel] kIV=" << kIV_ << " CeqI=" << CeqI_
              << " CeqV=" << CeqV_ << " at T=" << T_ << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;
    this->field_->addSpecies("Interstitial");
    this->field_->addSpecies("Vacancy");

    NumericType I = this->field_->getTotalDose("Interstitial");
    NumericType V = this->field_->getTotalDose("Vacancy");
    if (I <= 0 || V <= 0) return;

    NumericType recomb = kIV_ * std::min(I, V) * dt;
    recomb = std::min(recomb, std::min(I, V) * NumericType(0.5));

    if (I > 0) this->field_->scaleProfile("Interstitial", (I - recomb) / I);
    if (V > 0) this->field_->scaleProfile("Vacancy", (V - recomb) / V);

    this->field_->addDose("RecombinedIV", recomb);

    std::cout << "[ChargedReactKernel] recomb=" << recomb << " dt=" << dt
              << " I=" << this->field_->getTotalDose("Interstitial")
              << " V=" << this->field_->getTotalDose("Vacancy") << "\n";
  }

  void addToRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) &&
        myIndex < static_cast<int>(y.size())) {
      ydot[myIndex] += -kIV_ * y[myIndex];
    }
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto offI = this->field_->getSpeciesOffset("Interstitial");
    auto offV = this->field_->getSpeciesOffset("Vacancy");
    if (offI == static_cast<std::size_t>(-1) || offV == static_cast<std::size_t>(-1))
      return;
    auto n = this->field_->getProfileSize();
    if (offI + n > y.size() || offV + n > y.size()) return;

    for (std::size_t i = 0; i < n; ++i) {
      NumericType Ii = std::max(NumericType(0), y[offI + i]);
      NumericType Vi = std::max(NumericType(0), y[offV + i]);
      NumericType rate = kIV_ * Ii * Vi / (CeqI_ * CeqV_ + NumericType(1));
      // Cap rate for CVODE stability
      rate = std::min(rate, std::min(Ii, Vi) * NumericType(0.5));
      ydot[offI + i] += -rate;
      ydot[offV + i] += -rate;
    }
  }

private:
  NumericType T_ = 1273.15;
  NumericType kIV_ = NumericType(0.1);
  NumericType CeqI_ = 1e15;
  NumericType CeqV_ = 1e15;
};

} // namespace viennaps
