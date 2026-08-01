#pragma once

/// ImplantDamageCoupler — bridge MCBca implant damage into DiffusionEngine
/// initial conditions (dopant, interstitial, vacancy) with dose-conserving
/// 1D depth-profile projection. TED physics factories are added in Task 4.
///
/// Stateless: reuses DiffusionEngine::projectIntegralPreserving (the
/// IntegralPreservingFunctionIC-style hook) so it adds no solver or
/// assembly code. GAP_ANALYSIS §4.3: commercial tools feed implant damage
/// into the anneal step automatically; this is that bridge.

#include "DiffusionEngine.hpp"
#include "models/ConstantDiffusion.hpp"
#include "models/PairDiffusion.hpp"
#include "models/psMCBcaImplant.hpp"

#include <memory>
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

  /// TED pair-diffusion factory: D_eff = D_pair * C_I / C_I_eq.
  static std::shared_ptr<PairDiffusion<NumericType>>
  makeTedPair(const std::string &dopant, const std::string &interstitial,
              NumericType D_pair, NumericType C_Ieq) {
    auto m =
        std::make_shared<PairDiffusion<NumericType>>(dopant, interstitial);
    m->setPairDiffusivity(D_pair);
    m->setCIEq(C_Ieq);
    return m;
  }

  /// Defect transport factory (Arrhenius constant diffusivity).
  static std::shared_ptr<ConstantDiffusion<NumericType>>
  makeDefectTransport(const std::string &species, NumericType D0,
                      NumericType Ea_eV) {
    auto m = std::make_shared<ConstantDiffusion<NumericType>>(species);
    m->setDiffusivity(D0, Ea_eV);
    return m;
  }
};

} // namespace viennaps
