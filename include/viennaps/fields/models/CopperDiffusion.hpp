#pragma once

/// CopperDiffusion — MobileImpurity specialization for Cu (Phase 4 Task 7).
///
/// D = D0 * (1 + beta * C_dopant / ni) with optional Nernst–Planck drift
/// J = −D (∇C + (q/kT) z C E) and Cu + acceptor ⇌ CuA pairing.

#include "MobileImpurity.hpp"

#include <string>

namespace viennaps {

template <class NumericType>
class CopperDiffusion : public MobileImpurity<NumericType> {
public:
  CopperDiffusion() : MobileImpurity<NumericType>("Copper") {
    this->setName("CopperDiffusion");
    this->setD0(NumericType(1e-5)); // Cu is fast
    this->setIonPairing(NumericType(1));
    this->setChargeState(NumericType(1)); // Cu+
    // Pairing species enabled when setPairingRates + enablePairSpecies.
  }

  /// Enable Cu + acceptor ⇌ CuA (adds second species "CopperPair").
  void enablePairSpecies(const std::string &pair = "CopperPair",
                         const std::string &acceptor = "Boron") {
    this->pairSpecies_ = pair;
    this->setAcceptorSpecies(acceptor);
  }
};

} // namespace viennaps
