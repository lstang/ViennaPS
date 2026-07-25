#pragma once

/// KmcEvent — event types and rates for BKL KMC.

#include <cmath>

namespace viennaps {

enum class KmcEventType {
  Hop,
  Recombine,
  Cluster,
  Dissociate,
  Deposit,
  Desorb
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
};

// Species codes: 0 empty/Si, 1 interstitial, 2 vacancy, 3 {311}-like cluster
enum KmcSpeciesCode : int {
  KmcEmpty = 0,
  KmcInterstitial = 1,
  KmcVacancy = 2,
  KmcCluster311 = 3
};

} // namespace viennaps
