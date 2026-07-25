#pragma once

/// MaterialPropertySystem - Centralized, queryable material properties for
/// high-fidelity process simulation (diffusion, stress, oxidation, implantation).
///
/// Goals (matching ATHENA / SProcess manuals):
/// - Support Arrhenius forms:   X = X0 * exp(-Ea / (kB * T))
/// - Concentration, stress, and material-phase dependent properties
/// - "Like materials" interpolation / inheritance (e.g. doped oxide vs pure SiO2)
/// - Easy registration of new species (Boron, Phosphorus, Interstitial, etc.)
///
/// This will feed all PhysicsKernels and the Oxidation adapter.
/// In full impl, properties can be temperature-, stress- and concentration-dependent
/// functions and can be overridden per-region (via level-set material ids).

#include <map>
#include <string>
#include <cmath>
#include <functional>
#include <iostream>

namespace viennaps {

template <class NumericType>
class MaterialPropertySystem {
public:
  using PropertyFunc = std::function<NumericType(NumericType /*T_K*/, NumericType /*conc*/, NumericType /*stress*/)>;

  MaterialPropertySystem() {
    registerDefaultProperties();
  }

  // Register a constant value for a material + property key
  void setProperty(const std::string& material, const std::string& key, NumericType value) {
    props_[material][key] = [value](NumericType, NumericType, NumericType) { return value; };
  }

  // Register an Arrhenius pair: pre-factor and activation energy (eV)
  // Result = D0 * exp( -Ea / (kB * T) )   with kB = 8.617333262145e-5 eV/K
  void setArrhenius(const std::string& material, const std::string& key,
                    NumericType D0, NumericType Ea_eV) {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    props_[material][key] = [D0, Ea_eV, kB](NumericType T_K, NumericType, NumericType) {
      if (T_K <= 0) return D0;
      return D0 * std::exp( -Ea_eV / (kB * T_K) );
    };
  }

  // Register a fully custom function
  void setPropertyFunction(const std::string& material, const std::string& key, PropertyFunc func) {
    props_[material][key] = std::move(func);
  }

  // Query a property. Falls back to "default" material if not found for specific material.
  NumericType getProperty(const std::string& material, const std::string& key,
                          NumericType T_K = 1273.15,   // ~1000 C default
                          NumericType conc = 0,
                          NumericType hydrostaticStress = 0) const {
    auto matIt = props_.find(material);
    if (matIt != props_.end()) {
      auto pIt = matIt->second.find(key);
      if (pIt != matIt->second.end()) {
        return pIt->second(T_K, conc, hydrostaticStress);
      }
    }
    // fallback
    auto defIt = props_.find("default");
    if (defIt != props_.end()) {
      auto pIt = defIt->second.find(key);
      if (pIt != defIt->second.end()) {
        return pIt->second(T_K, conc, hydrostaticStress);
      }
    }
    return NumericType(0);
  }

  // Convenience: get diffusivity for a dopant/defect species
  NumericType getDiffusivity(const std::string& species,
                             const std::string& material = "Si",
                             NumericType T_K = 1273.15,
                             NumericType conc = 0,
                             NumericType stress = 0) const {
    std::string key = species + "_D";
    return getProperty(material, key, T_K, conc, stress);
  }

  // Convenience: get equilibrium concentration (for defects)
  NumericType getEquilibriumConcentration(const std::string& species,
                                          const std::string& material = "Si",
                                          NumericType T_K = 1273.15) const {
    std::string key = species + "_Ceq";
    return getProperty(material, key, T_K, 0, 0);
  }

  // Elastic / stress properties
  NumericType getYoungModulus(const std::string& material = "Si", NumericType T_K = 1273.15) const {
    return getProperty(material, "YoungModulus", T_K);
  }
  NumericType getPoissonRatio(const std::string& material = "Si", NumericType T_K = 1273.15) const {
    return getProperty(material, "PoissonRatio", T_K);
  }
  NumericType getViscosity(const std::string& material = "SiO2", NumericType T_K = 1273.15) const {
    return getProperty(material, "Viscosity", T_K);
  }
  NumericType getCTE(const std::string& material = "Si", NumericType T_K = 1273.15) const {
    return getProperty(material, "CTE", T_K);
  }

  // Print known materials/properties (debug)
  void printSummary() const {
    std::cout << "[MaterialPropertySystem] Registered materials:\n";
    for (const auto& [mat, ps] : props_) {
      std::cout << "  " << mat << " (" << ps.size() << " props)\n";
    }
  }

private:
  std::map<std::string, std::map<std::string, PropertyFunc>> props_;

  void registerDefaultProperties() {
    // Very rough defaults inspired by common TCAD values (will be calibrated later)
    // Units are typical process units: cm^2/s for D, cm^-3 for concentrations, etc.
    // Temperature in Kelvin.

    // Silicon defaults
    setArrhenius("Si", "Boron_D",     0.76,   3.46);   // rough
    setArrhenius("Si", "Phosphorus_D", 3.85,   3.66);
    setArrhenius("Si", "Arsenic_D",    0.32,   3.56);

    setArrhenius("Si", "Interstitial_D",  600.0, 1.5);   // fast diffusing
    setArrhenius("Si", "Vacancy_D",       0.001, 2.0);

    setProperty("Si", "Interstitial_Ceq", 1e15);
    setProperty("Si", "Vacancy_Ceq",      1e15);
    // Arrhenius equilibrium (Phase 3 PointDefectEquilibrium): order-of-magnitude
    // so C_I_eq(1273) lands in ~1e10–1e16 cm^-3 TCAD range.
    // Ceq = C0 * exp(-Ef/(kB*T)); pick C0/Ef so 1273 K is ~1e11–1e14.
    setProperty("Si", "Interstitial_Ceq0", NumericType(2.9e24));
    setProperty("Si", "Interstitial_Ef", NumericType(3.46));
    setProperty("Si", "Vacancy_Ceq0", NumericType(1.4e23));
    setProperty("Si", "Vacancy_Ef", NumericType(2.6));

    // Band-structure parameters for IntrinsicCarrier (Phase 2).
    // Temperature-independent constants are sufficient for Phase 2;
    // N_c, N_v ~ T^{3/2} can be added later if needed.
    // Typical Si at 300 K: N_c ~ 2.8e19, N_v ~ 1.04e19 cm^-3, E_g ~ 1.12 eV.
    // n_i(300) = sqrt(Nc*Nv)*exp(-Eg/(2*kB*T)) ~ 6.7e9 ≈ 1e10 cm^-3.
    setProperty("Si", "Nc", NumericType(2.8e19));
    setProperty("Si", "Nv", NumericType(1.04e19));
    setProperty("Si", "Eg", NumericType(1.12));

    // Oxide
    setArrhenius("SiO2", "Boron_D",  3.0e-3, 3.5);
    setArrhenius("SiO2", "Interstitial_D",  0.05, 2.0);

    // Default fallback (very small)
    setProperty("default", "Boron_D",       1e-20);
    setProperty("default", "Interstitial_D", 1e-10);

    // Stress dependence placeholder (will be multiplied in kernels)
    setProperty("default", "stressFactor", 1.0);

    // Elastic / viscoelastic / thermal properties (for stress module)
    // Young's modulus (GPa), Poisson, viscosity (Poise), CTE (1/K)
    setProperty("Si", "YoungModulus", 130.0);     // ~130 GPa for Si <100>
    setProperty("Si", "PoissonRatio", 0.28);
    setProperty("Si", "Viscosity", 1e20);         // very high at room, lower at high T
    setProperty("Si", "CTE", 2.6e-6);

    setProperty("SiO2", "YoungModulus", 70.0);
    setProperty("SiO2", "PoissonRatio", 0.17);
    setProperty("SiO2", "Viscosity", 1e13);       // typical for viscous flow in ox
    setProperty("SiO2", "CTE", 0.5e-6);

    setProperty("Si", "GrowthStress", 300.0e6);   // Pa order of magnitude for oxidation stress
  }
};

} // namespace viennaps