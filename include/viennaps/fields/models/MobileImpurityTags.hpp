#pragma once

/// Tag structs for compile-time species identification in MobileImpurity.
/// Each tag carries compile-time constants: name, D0, Ea, charge, default
/// pairing rate.  GenericImpurityTag preserves the existing runtime-string
/// behaviour so all existing call sites compile unchanged.

#include <string>

namespace viennaps {

struct GenericImpurityTag {
  static constexpr const char *name = "Impurity";
  static constexpr double D0 = 1e-10;
  static constexpr double Ea = 0.0;
  static constexpr double charge = 1.0;
  static constexpr double pairRate = 0.0;
};

struct CopperTag {
  static constexpr const char *name = "Copper";
  static constexpr double D0 = 1e-5;
  static constexpr double Ea = 0.0;
  static constexpr double charge = 1.0;
  static constexpr double pairRate = 0.0;
};

struct SodiumTag {
  static constexpr const char *name = "Sodium";
  static constexpr double D0 = 1e-8;
  static constexpr double Ea = 0.0;
  static constexpr double charge = 1.0;
  static constexpr double pairRate = 0.0;
};

struct IronTag {
  static constexpr const char *name = "Iron";
  static constexpr double D0 = 1e-9;
  static constexpr double Ea = 0.0;
  static constexpr double charge = 2.0;
  static constexpr double pairRate = 0.0;
};

} // namespace viennaps