#pragma once

/// SPERKernel - Solid Phase Epitaxial Regrowth of amorphous Si.
/// Consumes AmorphousFraction / AmorphousDepth, restores crystalline Si,
/// releases trapped dopants, and emits residual end-of-range defects.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace viennaps {

template <class NumericType>
class SPERKernel : public PhysicsKernel<NumericType> {
public:
  explicit SPERKernel(NumericType temperatureK = 873.15) // ~600 C typical SPER
      : T_(temperatureK) {
    this->setName("SPER");
  }

  void setTemperature(NumericType T) { T_ = T; }

  /// Orientation factor: (100)=1, (110)~0.7, (111)~0.5 typical SPER ratio.
  void setOrientation(const std::string &ori) {
    if (ori == "110")
      orientFactor_ = NumericType(0.7);
    else if (ori == "111")
      orientFactor_ = NumericType(0.5);
    else
      orientFactor_ = NumericType(1.0);
  }
  NumericType orientationFactor() const { return orientFactor_; }

  void setup() override {
    // v = v0 exp(-Ea/kT)  (cm/s → normalized units)
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    NumericType v0 = NumericType(1e7);
    NumericType Ea = NumericType(2.7);
    if (this->material_) {
      // Prefer ParameterDatabase-fed MaterialPropertySystem keys if present
      NumericType pv = this->material_->getProperty("Si", "SPER_Velocity", T_);
      if (pv > 0) {
        // interpret as already Arrhenius-evaluated velocity proxy
        velocity_ = pv;
      } else {
        velocity_ = v0 * std::exp(-Ea / (kB * std::max(T_, NumericType(1))));
      }
    } else {
      velocity_ = v0 * std::exp(-Ea / (kB * std::max(T_, NumericType(1))));
    }
    velocity_ *= orientFactor_;
    // Normalize for demo evolve steps
    velocity_ = std::min(NumericType(1), std::max(NumericType(1e-6), velocity_ * NumericType(1e-6)));
    std::cout << "[SPERKernel] T=" << T_ << " velocity=" << velocity_
              << " orient=" << orientFactor_ << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;
    this->field_->addSpecies("AmorphousFraction");
    this->field_->addSpecies("CrystallineFraction");
    this->field_->addSpecies("EOR_Defects");

    NumericType a0 = this->field_->getTotalDose("AmorphousFraction");
    if (a0 <= 0) {
      // seed tiny amorphous if none (so tests can still exercise)
      std::cout << "[SPERKernel] no amorphous dose; skip\n";
      return;
    }

    // Regrow fraction of amorphous layer
    NumericType regrow = std::min(a0, velocity_ * dt * a0);
    regrow = std::min(regrow, a0 * NumericType(0.95));
    if (a0 > 0)
      this->field_->scaleProfile("AmorphousFraction", (a0 - regrow) / a0);

    this->field_->addDose("CrystallineFraction", regrow);
    // End-of-range defects ~ 10% of regrown amorphous
    this->field_->addDose("EOR_Defects", regrow * NumericType(0.1));
    // Mild dopant reactivation (increase free Dopant slightly)
    if (this->field_->hasSpecies("Dopant")) {
      NumericType B = this->field_->getTotalDose("Dopant");
      this->field_->addDose("Dopant", regrow * NumericType(0.01) * std::max(B, NumericType(1)) /
                                          std::max(a0, NumericType(1)));
    }

    lastRegrown_ = regrow;
    std::cout << "[SPERKernel] regrown=" << regrow << " a0=" << a0
              << " EOR=" << this->field_->getTotalDose("EOR_Defects") << "\n";
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto off = this->field_->getSpeciesOffset("AmorphousFraction");
    if (off == static_cast<std::size_t>(-1)) return;
    auto n = this->field_->getProfileSize();
    if (off + n > y.size()) return;
    for (std::size_t i = 0; i < n; ++i)
      ydot[off + i] += -velocity_ * y[off + i];
  }

  NumericType getLastRegrown() const { return lastRegrown_; }

private:
  NumericType T_ = 873.15;
  NumericType velocity_ = NumericType(0.1);
  NumericType orientFactor_ = NumericType(1);
  NumericType lastRegrown_ = 0;
};

} // namespace viennaps
