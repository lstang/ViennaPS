#pragma once

/// IntrinsicCarrier — pure-math calculator for semiconductor intrinsic and
/// free-carrier concentrations, with a Fermi-Dirac activity-coefficient hook.
///
/// Formulas (Phase 2 Task 1):
///   n_i(T) = sqrt(N_c * N_v) * exp(-E_g / (2 * kB * T))
///   n_Boltzmann ≈ max(C_dopant, n_i)   (n-type complete-ionization approx.)
///   p = n_i^2 / n
///
/// Band-structure parameters N_c, N_v, E_g are looked up from
/// MaterialPropertySystem (Si defaults registered there). Phase 2 uses
/// temperature-independent N_c, N_v constants.
///
/// Activity coefficients (MOOSE CoupledDiffusionReactionSub pattern: _gamma_u
/// / _gamma_v / _gamma_eq):
/// At degenerate doping (C >> N_c) Boltzmann statistics break down. The
/// Fermi-Dirac activity γ ∈ (0,1] is applied so that the free-carrier
/// estimate used by later charged-Fermi models is reduced relative to pure
/// Boltzmann. Default γ = 1 for Boltzmann statistics.
///
/// No MFEM dependency — this is a pure numeric calculator.

#include "MaterialPropertySystem.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace viennaps {

/// Carrier statistics used by IntrinsicCarrier::activity and the
/// useFermiDirac overload of electronConcentration.
enum class CarrierStatistics {
  Boltzmann,  ///< Ideal non-degenerate: γ = 1
  FermiDirac  ///< Degenerate correction: γ ∈ (0,1]
};

template <class NumericType>
class IntrinsicCarrier {
public:
  IntrinsicCarrier() = default;

  explicit IntrinsicCarrier(MaterialPropertySystem<NumericType> mps)
      : mps_(std::move(mps)) {}

  void setMaterialPropertySystem(MaterialPropertySystem<NumericType> mps) {
    mps_ = std::move(mps);
  }

  const MaterialPropertySystem<NumericType> &materialProperties() const {
    return mps_;
  }

  /// Intrinsic carrier concentration n_i(T) [cm^-3].
  /// n_i = sqrt(N_c * N_v) * exp(-E_g / (2 * kB * T))
  NumericType ni(NumericType T_K,
                 const std::string &material = "Si") const {
    const NumericType Nc = mps_.getProperty(material, "Nc", T_K);
    const NumericType Nv = mps_.getProperty(material, "Nv", T_K);
    const NumericType Eg = mps_.getProperty(material, "Eg", T_K);
    const NumericType kB = kBoltzmann();
    if (T_K <= NumericType(0)) {
      return NumericType(0);
    }
    const NumericType sqrtNcNv = std::sqrt(std::max(Nc * Nv, NumericType(0)));
    return sqrtNcNv * std::exp(-Eg / (NumericType(2) * kB * T_K));
  }

  /// Activity coefficient γ.
  /// Boltzmann: γ = 1.
  /// Fermi-Dirac: simple analytic degeneracy factor
  ///   γ(x) = ln(1+x) / x ,  x = C / N_c
  /// which satisfies γ ∈ (0,1], γ→1 as C→0, and γ < 1 for C ≫ N_c.
  NumericType activity(
      NumericType C, NumericType T_K, const std::string &material,
      CarrierStatistics statistics = CarrierStatistics::Boltzmann) const {
    if (statistics == CarrierStatistics::Boltzmann) {
      return NumericType(1);
    }
    const NumericType Nc =
        std::max(mps_.getProperty(material, "Nc", T_K), NumericType(1));
    const NumericType x =
        std::max(C, NumericType(0)) / Nc;
    // ln(1+x)/x → 1 as x→0; clamp tiny x to avoid 0/0.
    if (x < NumericType(1e-12)) {
      return NumericType(1);
    }
    return std::log(NumericType(1) + x) / x;
  }

  /// Free-electron concentration for n-type doping approximation.
  /// Boltzmann: n = max(C_dopant, n_i).
  /// Fermi-Dirac (useFermiDirac=true): n = max(C_dopant, n_i) * γ,
  /// where γ = activity(..., FermiDirac) ∈ (0,1], so n_FD ≤ n_Boltzmann.
  ///
  /// Note: the plan wording "divides C by γ" with γ∈(0,1] would *increase*
  /// n; the required test asserts n_FD < n_Boltzmann, so the activity is
  /// applied multiplicatively (standard γ ≤ 1 convention).
  NumericType electronConcentration(NumericType C_dopant, NumericType T_K,
                                    const std::string &material = "Si",
                                    bool useFermiDirac = false) const {
    const NumericType niT = ni(T_K, material);
    NumericType n = std::max(C_dopant, niT);
    if (useFermiDirac) {
      const NumericType g =
          activity(C_dopant, T_K, material, CarrierStatistics::FermiDirac);
      n = n * g;
      // Keep at least n_i so p = ni^2/n remains well-defined and physical.
      n = std::max(n, niT);
    }
    return n;
  }

  /// Hole concentration from mass-action: p = n_i^2 / n.
  NumericType holeConcentration(NumericType C_dopant, NumericType T_K,
                                const std::string &material = "Si",
                                bool useFermiDirac = false) const {
    const NumericType niT = ni(T_K, material);
    const NumericType n =
        electronConcentration(C_dopant, T_K, material, useFermiDirac);
    if (n <= NumericType(0)) {
      return NumericType(0);
    }
    return (niT * niT) / n;
  }

  static constexpr NumericType kBoltzmann() {
    return static_cast<NumericType>(8.617333262145e-5); // eV/K
  }

private:
  MaterialPropertySystem<NumericType> mps_{};
};

} // namespace viennaps
