// testKmcEpitaxy.cpp — LKMC epitaxy production-hardening tests
// (GAP_ANALYSIS §4.9): coordination-scaled attachment, SiGe composition,
// z-buffer visibility, desorption balance, rare twins.
#include <fields/kmc/KmcAtomisticEngine.hpp>
#include <fields/kmc/KmcEpitaxy.hpp>
#include <fields/kmc/KmcEvent.hpp>
#include <fields/kmc/KmcLattice.hpp>
#include <vcTestAsserts.hpp>

#include <algorithm>
#include <cstdio>
#include <iostream>

namespace viennacore {

using namespace viennaps;

template <class NumericType, int D> void RunTest() {
  // Seed a full (8x8) Si floor at k=0.
  auto seedFloor = [](KmcLattice &lat) {
    for (int i = 0; i < lat.nx(); ++i)
      for (int j = 0; j < lat.ny(); ++j) {
        auto &s = lat.at(i, j, 0);
        s.occupied = true;
        s.species = KmcSi;
      }
  };
  auto makeEngine = [&](unsigned seed) {
    KmcLattice lat;
    lat.resize(8, 8, 8);
    seedFloor(lat);
    KmcAtomisticEngine<NumericType> engine(seed);
    engine.setLattice(lat);
    engine.setEpitaxyEnabled(true);
    engine.setDiamondNeighbors(true);
    return engine;
  };

  // --- Test 1: bottom-up growth with coordination scaling.
  {
    std::cout << "[kmc-epi] Test 1 start" << std::endl;
    auto engine = makeEngine(42);
    KmcParameters p;
    p.desorbPreFactor = 0; // isolate attachment
    p.twinPreFactor = 0;
    engine.setParameters(p);
    engine.run(64); // one layer worth of sites
    const KmcLattice &out = engine.lattice();
    int l1 = 0, l2 = 0;
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 8; ++j) {
        if (out.at(i, j, 1).occupied)
          ++l1;
        if (out.at(i, j, 2).occupied)
          ++l2;
      }
    std::cout << "[kmc-epi] deposits=" << engine.depositCount() << " l1=" << l1
              << " l2=" << l2 << "\n";
    VC_TEST_ASSERT(engine.depositCount() == 64); // no desorb/twin
    VC_TEST_ASSERT(l1 > 0);
    VC_TEST_ASSERT(l1 >= l2); // lower layers fill before upper ones
  }

  // --- Test 2: Ge mole fraction selects KmcGe species.
  {
    std::cout << "[kmc-epi] Test 2 start" << std::endl;
    auto runGe = [&](double xGe) {
      auto engine = makeEngine(7);
      KmcParameters p;
      p.desorbPreFactor = 0;
      p.twinPreFactor = 0;
      engine.setParameters(p);
      engine.setGeFraction(xGe);
      engine.run(200);
      return engine.lattice().countSpecies(KmcGe);
    };
    VC_TEST_ASSERT(runGe(0.0) == 0);
    VC_TEST_ASSERT(runGe(0.5) > 0); // deterministic seed → stable count
  }

  // --- Test 3: visibility — a floating atom shadows its column.
  {
    std::cout << "[kmc-epi] Test 3 start" << std::endl;
    auto engine = makeEngine(11);
    KmcParameters p;
    p.desorbPreFactor = 0;
    p.twinPreFactor = 0;
    engine.setParameters(p);
    engine.setVisibilityEnabled(true);
    // Floater above column (0,0): blocks its line of sight from k=2 up,
    // so the (0,0,1) deposit site is shadowed and must not grow.
    auto &f = engine.lattice().at(0, 0, 3);
    f.occupied = true;
    f.species = KmcSi;
    engine.run(200);
    const KmcLattice &out = engine.lattice();
    VC_TEST_ASSERT(!out.at(0, 0, 1).occupied); // shadowed column does not grow
    VC_TEST_ASSERT(out.countSpecies(KmcSi) > 64); // other columns still grow
  }

  // --- Test 4: minCoordination gate — c >= 5 impossible (diamond).
  {
    std::cout << "[kmc-epi] Test 4 start" << std::endl;
    auto engine = makeEngine(13);
    KmcParameters p;
    p.desorbPreFactor = 0;
    p.twinPreFactor = 0;
    engine.setParameters(p);
    engine.setMinCoordination(5); // diamond max coordination is 4
    engine.run(100);
    VC_TEST_ASSERT(engine.depositCount() == 0);
  }

  // --- Test 5: desorption balances attachment at high desorb rate.
  {
    std::cout << "[kmc-epi] Test 5 start" << std::endl;
    auto engine = makeEngine(17);
    KmcParameters p;
    p.attachPreFactor = 1e6;  // ~1.05e4 Hz at 1273 K
    p.desorbPreFactor = 1e13; // ~1.1e7 Hz — dominates attachment
    p.twinPreFactor = 0;
    engine.setParameters(p);
    engine.run(2000);
    std::cout << "[kmc-epi] steps=" << engine.steps()
              << " deposit=" << engine.depositCount()
              << " desorb=" << engine.desorbCount() << std::endl;
    VC_TEST_ASSERT(engine.steps() > 0);
    VC_TEST_ASSERT(engine.desorbCount() > engine.depositCount());
  }

  // --- Test 6: twins are a rare minority.
  {
    std::cout << "[kmc-epi] Test 6 start" << std::endl;
    auto engine = makeEngine(3);
    KmcParameters p;
    p.desorbPreFactor = 0;
    p.twinPreFactor = 1e6;
    engine.setParameters(p);
    engine.run(500);
    std::cout << "[kmc-epi] deposit=" << engine.depositCount()
              << " twin=" << engine.twinCount() << std::endl;
    VC_TEST_ASSERT(engine.depositCount() > 0);
    VC_TEST_ASSERT(engine.twinCount() < engine.depositCount() * 0.15);
  }

  // --- Test 7: KmcEpitaxyModel rate-based parity.
  {
    std::cout << "[kmc-epi] Test 7 start" << std::endl;
    KmcLattice lat;
    lat.resize(8, 8, 8);
    seedFloor(lat);
    KmcEpitaxyModel model;
    model.setGrowthRate(1.0);
    model.setGeFraction(0.0);
    const int deposits = model.runRateBased(lat, 64);
    std::cout << "[kmc-epi] model.runRateBased deposits=" << deposits
              << std::endl;
    VC_TEST_ASSERT(deposits > 0);
    VC_TEST_ASSERT(lat.countSpecies(KmcSi) > 64);
  }
}

} // namespace viennacore

int main() {
  std::setvbuf(stdout, NULL, _IONBF, 0);
  VC_RUN_ALL_TESTS
}
