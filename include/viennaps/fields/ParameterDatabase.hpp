#pragma once

/// ParameterDatabase - Manual-style material parameter DB with inheritance
/// and like-materials interpolation (ATHENA/SProcess parity goal).
///
/// Layers:
///   1. Base material defaults (Si, SiO2, Si3N4, PolySi, SiGe, ...)
///   2. Species-specific Arrhenius keys (Boron_D, Interstitial_Ceq, ...)
///   3. Inheritance: child material falls back to parent
///   4. Like-materials blend: weighted interpolation between two materials

#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

namespace viennaps {

template <class NumericType>
class ParameterDatabase {
public:
  ParameterDatabase() { registerBuiltins(); }

  void setParent(const std::string& child, const std::string& parent) {
    parents_[child] = parent;
  }

  void setArrhenius(const std::string& material, const std::string& key,
                    NumericType D0, NumericType Ea_eV) {
    mats_.setArrhenius(material, key, D0, Ea_eV);
  }

  void setProperty(const std::string& material, const std::string& key,
                   NumericType value) {
    mats_.setProperty(material, key, value);
  }

  /// Query with inheritance chain: material → parent → ... → default.
  NumericType get(const std::string& material, const std::string& key,
                  NumericType T_K = 1273.15, NumericType conc = 0,
                  NumericType stress = 0) const {
    // Direct material query (MaterialPropertySystem already falls back to "default")
    NumericType v = mats_.getProperty(material, key, T_K, conc, stress);
    NumericType defV = mats_.getProperty("default", key, T_K, conc, stress);
    if (v != defV || material == "default") return v;

    // Walk parents when child only hit the default fallback
    std::string m = material;
    for (int guard = 0; guard < 8; ++guard) {
      auto pit = parents_.find(m);
      if (pit == parents_.end()) break;
      m = pit->second;
      v = mats_.getProperty(m, key, T_K, conc, stress);
      if (v != defV || m == "default") return v;
    }
    return defV;
  }

  /// Like-materials interpolation: blend = (1-w)*A + w*B
  NumericType blend(const std::string& matA, const std::string& matB,
                    const std::string& key, NumericType weightB,
                    NumericType T_K = 1273.15) const {
    weightB = std::clamp(weightB, NumericType(0), NumericType(1));
    NumericType a = get(matA, key, T_K);
    NumericType b = get(matB, key, T_K);
    return (NumericType(1) - weightB) * a + weightB * b;
  }

  MaterialPropertySystem<NumericType>& materialSystem() { return mats_; }
  const MaterialPropertySystem<NumericType>& materialSystem() const { return mats_; }

  void printSummary() const {
    std::cout << "[ParameterDatabase] parents:\n";
    for (const auto& [c, p] : parents_)
      std::cout << "  " << c << " -> " << p << "\n";
    mats_.printSummary();
  }

  int parentCount() const { return static_cast<int>(parents_.size()); }

private:
  MaterialPropertySystem<NumericType> mats_;
  std::map<std::string, std::string> parents_;

  void registerBuiltins() {
    // Inheritance tree
    setParent("SiGe", "Si");
    setParent("PolySi", "Si");
    setParent("DopedOxide", "SiO2");
    setParent("ThermalOxide", "SiO2");
    setParent("Nitride", "Si3N4");

    // Extra materials
    setArrhenius("Si3N4", "Boron_D", 1e-4, 3.8);
    setProperty("Si3N4", "YoungModulus", 300.0);
    setProperty("Si3N4", "PoissonRatio", 0.25);
    setProperty("Si3N4", "CTE", 3.0e-6);

    setArrhenius("PolySi", "Boron_D", 1.2, 3.4);
    setProperty("PolySi", "YoungModulus", 160.0);

    setArrhenius("SiGe", "Boron_D", 1.5, 3.3);
    setProperty("SiGe", "YoungModulus", 120.0);
    setProperty("SiGe", "CTE", 4.0e-6);

    setArrhenius("DopedOxide", "Boron_D", 5e-3, 3.2);

    // Silicide kinetics placeholders
    setArrhenius("TiSi2", "GrowthRate", 1e-4, 1.5);
    setProperty("TiSi2", "YoungModulus", 250.0);
    setArrhenius("NiSi", "GrowthRate", 2e-4, 1.2);

    // SPER / amorph
    setArrhenius("Si", "SPER_Velocity", 1e7, 2.7);
    setProperty("Si", "AmorphThreshold", 1e22);

    // --- Defect / cluster / interface params (Phase 3 full-depth) ---
    // Order-of-magnitude TCAD defaults; SProcess 311 / recombination chapters
    // use calibrated tables — override via setProperty/setArrhenius for device work.
    // Point-defect equilibria (also registered on MaterialPropertySystem defaults).
    setProperty("Si", "Interstitial_Ceq0", NumericType(2.9e24));
    setProperty("Si", "Interstitial_Ef", NumericType(3.46));
    setProperty("Si", "Vacancy_Ceq0", NumericType(1.4e23));
    setProperty("Si", "Vacancy_Ef", NumericType(2.6));
    // Bulk I+V recombination rate scale [cm^3/s] at process T (test-friendly).
    setProperty("Si", "IV_Recombination_k", NumericType(1e-15));
    // {311} formation / dissociation (kf on C_I^n, kr linear).
    setProperty("Si", "Cluster311_kf", NumericType(1e-20));
    setProperty("Si", "Cluster311_kr", NumericType(1e-3));
    setProperty("Si", "Cluster311_n", NumericType(2));
    // Vacancy clusters
    setProperty("Si", "VacancyCluster_kf", NumericType(1e-20));
    setProperty("Si", "VacancyCluster_kr", NumericType(1e-3));
    // BIC / impurity-interstitial clusters
    setProperty("Si", "BIC_kf", NumericType(1e-18));
    setProperty("Si", "BIC_kr", NumericType(1e-3));
    // Dislocation loop growth from I supersaturation
    setProperty("Si", "DislocationLoop_k", NumericType(1e12));
    setProperty("Si", "DislocationLoop_p", NumericType(1));
    // Segregation coefficient m = C_mat2/C_mat1 at Si/SiO2 (B-like default)
    setProperty("Si", "Segregation_m_B_SiO2", NumericType(0.1));
    setProperty("SiO2", "Segregation_m_B_Si", NumericType(10.0));
    // Dose-loss / surface transfer velocity h [cm/s] placeholder
    setProperty("Si", "DoseLoss_h_Boron", NumericType(1e-4));
    // Pair diffusivity scale for TED (cm^2/s) when not taken from Arrhenius D
    setProperty("Si", "Pair_D_Boron", NumericType(1e-13));
  }
};

} // namespace viennaps
