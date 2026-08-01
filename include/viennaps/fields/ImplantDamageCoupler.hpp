#pragma once

/// ImplantDamageCoupler — bridge MCBca implant damage into DiffusionEngine
/// initial conditions (dopant, interstitial, vacancy) with dose-conserving
/// 1D depth-profile projection. TED physics factories are added in Task 4.
///
/// Stateless: reuses DiffusionEngine::projectIntegralPreserving (the
/// IntegralPreservingFunctionIC-style hook) so it adds no solver or
/// assembly code. GAP_ANALYSIS §4.3: commercial tools feed implant damage
/// into the anneal step automatically; this is that bridge.

#include "../models/psMCBcaImplant.hpp"
#include "DiffusionEngine.hpp"

#include <numeric>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType> class ImplantDamageCoupler {
public:
  /// Seed dopant + defect species from a BCA result onto the engine's
  /// species fields. Each 1D depth profile is projected with
  /// projectIntegralPreserving, rescaling the field integral to the
  /// profile's total dose (dose conservation per species).
  static void seedFromBca(DiffusionEngine<NumericType, 2> &engine,
                          const MCBcaResult<NumericType> &result) {
    seedSpecies(engine, "Dopant", result.dopant);
    seedSpecies(engine, "Interstitial", result.interstitial);
    seedSpecies(engine, "Vacancy", result.vacancy);
  }

  /// Project one 1D depth profile with dose conservation.
  static void seedSpecies(DiffusionEngine<NumericType, 2> &engine,
                          const std::string &name,
                          const std::vector<NumericType> &profile) {
    if (profile.empty())
      return;
    const NumericType dose =
        std::accumulate(profile.begin(), profile.end(), NumericType(0));
    engine.projectIntegralPreserving(name, profile, dose);
  }
};

} // namespace viennaps
