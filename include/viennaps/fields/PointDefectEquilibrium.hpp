#pragma once

/// PointDefectEquilibrium — C_I^eq(T), C_V^eq(T) for Si and other materials.
///
/// Arrhenius form: C_eq = C0 * exp(-E_f / (kB * T))
/// Defaults are order-of-magnitude TCAD values suitable for TED/OED tests.

#include "MaterialPropertySystem.hpp"

#include <cmath>
#include <string>

namespace viennaps {

template <class NumericType>
class PointDefectEquilibrium {
public:
  PointDefectEquilibrium() = default;

  explicit PointDefectEquilibrium(MaterialPropertySystem<NumericType> mps)
      : mps_(std::move(mps)) {}

  void setMaterialPropertySystem(MaterialPropertySystem<NumericType> mps) {
    mps_ = std::move(mps);
  }

  /// Equilibrium interstitial concentration [cm^-3].
  NumericType C_I_eq(NumericType T_K,
                     const std::string &material = "Si") const {
    return equilibrium("Interstitial", T_K, material);
  }

  /// Equilibrium vacancy concentration [cm^-3].
  NumericType C_V_eq(NumericType T_K,
                     const std::string &material = "Si") const {
    return equilibrium("Vacancy", T_K, material);
  }

  /// Generic species equilibrium via MaterialPropertySystem key "{species}_Ceq"
  /// or Arrhenius pair "{species}_Ceq0" / "{species}_Ef".
  NumericType equilibrium(const std::string &species, NumericType T_K,
                          const std::string &material = "Si") const {
    // Prefer full Arrhenius if prefactor + formation energy registered.
    const NumericType C0 =
        mps_.getProperty(material, species + "_Ceq0", T_K);
    const NumericType Ef =
        mps_.getProperty(material, species + "_Ef", T_K);
    if (C0 > NumericType(0) && Ef > NumericType(0) && T_K > NumericType(0)) {
      const NumericType kB =
          static_cast<NumericType>(8.617333262145e-5);
      return C0 * std::exp(-Ef / (kB * T_K));
    }
    // Fallback: constant Ceq property (Phase 1 defaults ~1e15).
    return mps_.getEquilibriumConcentration(species, material, T_K);
  }

private:
  MaterialPropertySystem<NumericType> mps_{};
};

} // namespace viennaps
