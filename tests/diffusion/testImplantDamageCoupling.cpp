// testImplantDamageCoupling.cpp — MCBca implant damage → DiffusionEngine
// initial-condition coupling (GAP_ANALYSIS §4.3).
#include <fields/DiffusionEngine.hpp>
#include <fields/DiffusionPhysics.hpp>
#include <fields/ImplantDamageCoupler.hpp>
#include <models/psMCBcaImplant.hpp>
#include <vcTestAsserts.hpp>

#include <cmath>
#include <memory>
#include <numeric>

namespace viennacore {

using namespace viennaps;

#ifdef VIENNAPS_HAS_MFEM

// Test 1: seeding from a deterministic BCA run conserves dose. The coupler
// projects each 1D depth profile onto the species GridFunction and rescales
// the field integral to the profile's total dose.
void TestImplantDamageSeeding() {
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");
  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

  MCBcaEngine<double> bca;
  bca.params().energyKeV = 50;
  bca.params().dose = 1e13;
  bca.params().nIons = 64;
  bca.params().nDepth = 64;
  bca.params().seed = 42;
  const MCBcaResult<double> result = bca.run();

  const double dopantDose =
      std::accumulate(result.dopant.begin(), result.dopant.end(), 0.0);

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);
  DiffusionPhysics<double> physics;
  physics.addSpecies("Dopant");
  physics.addSpecies("Interstitial");
  physics.addSpecies("Vacancy");
  engine.setPhysics(physics);

  ImplantDamageCoupler<double>::seedFromBca(engine, result);

  const double seeded = engine.getIntegral("Dopant");
  const double iDose = engine.getIntegral("Interstitial");
  const double vDose = engine.getIntegral("Vacancy");
  std::cout << "[implant-coupling] target=" << dopantDose
            << " seeded=" << seeded << " I=" << iDose << " V=" << vDose
            << "\n";
  VC_TEST_ASSERT(std::abs(seeded - dopantDose) /
                     std::max(dopantDose, 1.0) <
                 1e-3);
  VC_TEST_ASSERT(iDose > 0.0);
  VC_TEST_ASSERT(vDose > 0.0);
  VC_TEST_ASSERT(result.dopant.size() ==
                 static_cast<std::size_t>(bca.params().nDepth));
}

#endif // VIENNAPS_HAS_MFEM

} // namespace viennacore

int main() {
#ifdef VIENNAPS_HAS_MFEM
  try {
    viennacore::TestImplantDamageSeeding();
    std::cout << "All implant damage coupling tests passed.\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "TEST EXCEPTION: " << ex.what() << "\n";
    return 1;
  }
#else
  std::cout << "MFEM not available, skipping implant damage coupling tests.\n";
  return 0;
#endif
}
