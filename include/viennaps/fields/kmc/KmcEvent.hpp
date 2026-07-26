#pragma once

/// KmcEvent — event types and rates for BKL KMC.

#include <cmath>

namespace viennaps {

enum class KmcEventType {
  Hop,
  Recombine,
  Cluster,
  Dissociate,
  Deposit,   ///< Epitaxial surface attachment (Arrhenius, BKL event)
  Desorb,    ///< Surface desorption (Arrhenius, BKL event)
  Twin       ///< Twin-defect formation on {111} (Arrhenius, BKL event)
};

struct KmcEvent {
  KmcEventType type = KmcEventType::Hop;
  int i0 = 0, j0 = 0, k0 = 0;
  int i1 = 0, j1 = 0, k1 = 0;
  double rate = 0.0;
};

struct KmcParameters {
  double T = 1273.0;
  double hopPreFactor = 1e13; // Hz
  double hopBarrier = 0.7;    // eV
  double recombPreFactor = 1e12;
  double recombBarrier = 0.3;
  double clusterPreFactor = 1e12;
  double clusterBarrier = 0.5;
  double dissocPreFactor = 1e12;
  double dissocBarrier = 1.2; // binding + migration

  // Epitaxy rates (Phase 8 plan, Arrhenius form). Surface attachment,
  // desorption, and twin-defect formation as stochastic BKL events.
  double attachPreFactor = 1e6;   // Hz (gas-flux-limited)
  double attachBarrier = 0.5;     // eV (H-desorption-limited below ~600C)
  double desorbPreFactor = 1e10;
  double desorbBarrier = 1.5;     // eV
  double twinPreFactor = 1e8;
  double twinBarrier = 1.0;       // eV (stacking-fault formation)

  static constexpr double kB = 8.617333262145e-5;

  double hopRate() const {
    return hopPreFactor * std::exp(-hopBarrier / (kB * T));
  }
  double recombRate() const {
    return recombPreFactor * std::exp(-recombBarrier / (kB * T));
  }
  double clusterRate() const {
    return clusterPreFactor * std::exp(-clusterBarrier / (kB * T));
  }
  double dissocRate() const {
    return dissocPreFactor * std::exp(-dissocBarrier / (kB * T));
  }
  double attachRate() const {
    return attachPreFactor * std::exp(-attachBarrier / (kB * T));
  }
  double desorbRate() const {
    return desorbPreFactor * std::exp(-desorbBarrier / (kB * T));
  }
  double twinRate() const {
    return twinPreFactor * std::exp(-twinBarrier / (kB * T));
  }
};

// Species codes: 0 empty/Si, 1 interstitial, 2 vacancy, 3 {311}-like cluster,
// 4 Si (deposited), 5 Ge, 6 dopant, 7 twin-defect marker, 8 amorphous.
enum KmcSpeciesCode : int {
  KmcEmpty = 0,
  KmcInterstitial = 1,
  KmcVacancy = 2,
  KmcCluster311 = 3,
  KmcSi = 4,
  KmcGe = 5,
  KmcDopant = 6,
  KmcTwin = 7,
  KmcAmorphous = 8
};

} // namespace viennaps
