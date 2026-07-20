#pragma once

/// psSilicidation - Metal-silicon reaction forming silicide (TiSi2, NiSi, ...).
/// Consumes Si / metal dose, grows silicide thickness, injects stress and
/// optional dopant segregation into PhysicsField.

#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"
#include "fields/ParameterDatabase.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class SilicidationModel {
public:
  SilicidationModel(std::string silicide = "NiSi", NumericType temperatureK = 773.15)
      : silicide_(std::move(silicide)), T_(temperatureK) {}

  void setSilicide(const std::string& s) { silicide_ = s; }
  void setTemperature(NumericType T) { T_ = T; }
  void setMetalDose(NumericType d) { metalDose_ = d; }
  void setTime(NumericType t) { time_ = t; }

  void setParameterDatabase(std::shared_ptr<ParameterDatabase<NumericType>> db) {
    db_ = std::move(db);
  }

  /// Advance silicidation; returns grown thickness (normalized um proxy).
  NumericType evolve(PhysicsField<NumericType>& field) {
    NumericType k0 = NumericType(1e-4);
    NumericType Ea = NumericType(1.2);
    if (db_) {
      NumericType gr = db_->get(silicide_, "GrowthRate", T_);
      if (gr > 0) k0 = gr;
    }
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    NumericType rate = k0 * std::exp(-Ea / (kB * std::max(T_, NumericType(1))));
    // Normalize for demo
    rate = std::min(NumericType(0.5), std::max(NumericType(1e-6), rate * NumericType(1e3)));

    NumericType thickness = rate * time_ * std::min(NumericType(1), metalDose_ / NumericType(1e15));
    thickness = std::min(thickness, NumericType(0.2)); // cap

    field.addDose("SilicideThickness", thickness);
    field.addDose("SilicideMaterial",
                  silicide_ == "NiSi" ? NumericType(1) : NumericType(2));
    // Consume metal / Si proxy
    field.addDose("ConsumedMetal", metalDose_ * thickness);
    field.addDose("ConsumedSi", metalDose_ * thickness * NumericType(1.5));
    // Growth stress
    field.addDose("SilicideStress", thickness * NumericType(1e8));
    // Dopant segregation into silicide (sink)
    if (field.hasSpecies("Dopant")) {
      NumericType B = field.getTotalDose("Dopant");
      NumericType seg = B * thickness * NumericType(0.1);
      if (B > 0) field.scaleProfile("Dopant", std::max(NumericType(0.5), (B - seg) / B));
      field.addDose("SilicideDopant", seg);
    }

    lastThickness_ = thickness;
    std::cout << "[SilicidationModel] " << silicide_ << " T=" << T_
              << " thickness=" << thickness << " stress="
              << field.getTotalDose("SilicideStress") << "\n";
    return thickness;
  }

  NumericType getLastThickness() const { return lastThickness_; }

private:
  std::string silicide_ = "NiSi";
  NumericType T_ = 773.15;
  NumericType metalDose_ = 1e15;
  NumericType time_ = 60;
  NumericType lastThickness_ = 0;
  std::shared_ptr<ParameterDatabase<NumericType>> db_;
};

} // namespace viennaps
