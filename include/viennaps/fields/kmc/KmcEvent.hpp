#pragma once

/// KmcEvent — event types and rates for BKL KMC.

#include <cmath>
#include <string>

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
  double hopBarrier = 0.7;   // eV
  double recombPreFactor = 1e12;
  double recombBarrier = 0.3;

  double hopRate() const {
    const double kB = 8.617333262145e-5;
    return hopPreFactor * std::exp(-hopBarrier / (kB * T));
  }

  double recombRate() const {
    const double kB = 8.617333262145e-5;
    return recombPreFactor * std::exp(-recombBarrier / (kB * T));
  }
};

} // namespace viennaps
