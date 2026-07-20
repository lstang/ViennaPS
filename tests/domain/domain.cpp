#define VIENNAPS_HAS_SUNDIALS 1
#include <vcTestAsserts.hpp>

// Multiphysics Track 1 / Phase 1 smoke without full ViennaLS geometry includes.
#include "psVersion.hpp"
#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"
#include "fields/DiffusionKernel.hpp"
#include "fields/FermiDiffusionKernel.hpp"
#include "fields/PairDiffusionKernel.hpp"
#include "fields/ChargedReactKernel.hpp"
#include "fields/StressKernel.hpp"
#include "fields/DefectClusterKernel.hpp"
#include "fields/SundialsTimeIntegrator.hpp"
#include "fields/GeometryFieldCoupler.hpp"
#include "process/psPhysicsFieldAdapter.hpp"
#include "ProcessOrchestrator.hpp"

#include "fields/PhysicsKernel.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace viennacore {

template <class NumericType, int D>
void RunTest() {
  using Field = viennaps::PhysicsField<NumericType>;
  using MatSys = viennaps::MaterialPropertySystem<NumericType>;

  auto field = std::make_shared<Field>();
  field->setProfileSize(32); // keep CVODE dense state modest
  auto mats = std::make_shared<MatSys>();
  mats->setArrhenius("Si", "Dopant_D", 0.1, 2.5);
  mats->setArrhenius("Si", "Interstitial_D", 10.0, 1.8);
  mats->setProperty("Si", "YoungModulus", 130.0);
  mats->setProperty("Si", "PoissonRatio", 0.28);
  mats->setProperty("Si", "GrowthStress", 300.0e6);

  // Implant-like seed
  std::vector<NumericType> profile(32, static_cast<NumericType>(1e12 / 32));
  field->injectImplantProfile("Dopant", profile);
  field->injectImplantProfile("Interstitial", profile);
  field->injectImplantProfile("Vacancy",
                              std::vector<NumericType>(32, static_cast<NumericType>(5e11 / 32)));
  VC_TEST_ASSERT(field->getTotalDose("Dopant") > 0);
  VC_TEST_ASSERT(field->getStateSize() == field->getSpeciesOrder().size() * 32);

  // Pack / unpack round-trip
  {
    auto packed = field->packState();
    VC_TEST_ASSERT(packed.size() == field->getStateSize());
    auto dose0 = field->getTotalDose("Dopant");
    field->unpackState(packed);
    VC_TEST_ASSERT(std::abs(field->getTotalDose("Dopant") - dose0) < dose0 * NumericType(1e-6) + NumericType(1));
    std::cout << "[state-check] pack/unpack OK size=" << packed.size() << "\n";
  }

  // Geometry coupling (depth stack material map)
  {
    viennaps::GeometryFieldCoupler<NumericType> geo(field);
    geo.markFromDepthInterfaces(NumericType(0.15), NumericType(0.05));
    int matSurf = field->getMaterialAtNormalizedDepth(NumericType(0.02));
    int matBulk = field->getMaterialAtNormalizedDepth(NumericType(0.8));
    VC_TEST_ASSERT(matSurf == 3 || matSurf == 2); // mask or oxide near surface
    VC_TEST_ASSERT(matBulk == 1);                 // Si bulk
    std::cout << "[geo-check] surfaceMat=" << matSurf << " bulkMat=" << matBulk << "\n";
  }

  // Diffusion + Fermi + Pair + ChargedReact + stress
  {
    viennaps::DiffusionKernel<NumericType> k("Dopant", 1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }
  {
    viennaps::FermiDiffusionKernel<NumericType> k("Dopant", 1273.15, 1e18);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }
  {
    viennaps::PairDiffusionKernel<NumericType> k("Dopant", 1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("PairBI") > 0);
    std::cout << "[pair-check] PairBI=" << field->getTotalDose("PairBI") << "\n";
  }
  {
    auto I0 = field->getTotalDose("Interstitial");
    auto V0 = field->getTotalDose("Vacancy");
    viennaps::ChargedReactKernel<NumericType> k(1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("Interstitial") <= I0 + NumericType(1));
    VC_TEST_ASSERT(field->getTotalDose("Vacancy") <= V0 + NumericType(1));
    VC_TEST_ASSERT(field->getTotalDose("Vacancy") >= 0);
    std::cout << "[charged-check] I0=" << I0 << " V0=" << V0
              << " I=" << field->getTotalDose("Interstitial")
              << " V=" << field->getTotalDose("Vacancy") << "\n";
  }
  {
    viennaps::ViscoelasticStressKernel<NumericType> k(1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }
  {
    viennaps::ElasticStressKernel<NumericType> k(1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setMismatchStrain(NumericType(0.002));
    k.setup();
    k.evolve(30.0);
    VC_TEST_ASSERT(field->getTotalDose("ElasticStress") > 0);
    std::cout << "[elastic-check] ElasticStress=" << field->getTotalDose("ElasticStress")
              << "\n";
  }

  // --- Clustering models (recomb / 311 / bic / loop) ---
  {
    // Fresh non-negative I/V for recomb
    field->injectImplantProfile("Interstitial", profile);
    field->injectImplantProfile("Vacancy",
                                std::vector<NumericType>(32, static_cast<NumericType>(5e11 / 32)));

    auto I0 = field->getTotalDose("Interstitial");
    auto V0 = field->getTotalDose("Vacancy");
    VC_TEST_ASSERT(I0 > 0 && V0 > 0);

    viennaps::DefectClusterKernel<NumericType> recomb(1273.15, "recomb");
    recomb.setPhysicsField(field);
    recomb.setMaterialProperties(mats);
    recomb.setup();
    recomb.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("RecombinedIV") > 0);
    VC_TEST_ASSERT(field->getTotalDose("Vacancy") >= 0);
    VC_TEST_ASSERT(field->getTotalDose("Interstitial") >= 0);
    std::cout << "[cluster-check] recomb: I0=" << I0 << " V0=" << V0
              << " RecombinedIV=" << field->getTotalDose("RecombinedIV") << "\n";

    field->injectImplantProfile("Interstitial", profile);
    viennaps::DefectClusterKernel<NumericType> c311(1273.15, "311");
    c311.setPhysicsField(field);
    c311.setMaterialProperties(mats);
    c311.setup();
    c311.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("Cluster311") > 0);
    std::cout << "[cluster-check] 311: Cluster311=" << field->getTotalDose("Cluster311")
              << "\n";

    field->injectImplantProfile("Dopant", profile);
    field->injectImplantProfile("Interstitial", profile);
    viennaps::DefectClusterKernel<NumericType> bic(1273.15, "bic");
    bic.setPhysicsField(field);
    bic.setMaterialProperties(mats);
    bic.setup();
    bic.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("BIC") > 0);
    std::cout << "[cluster-check] bic: BIC=" << field->getTotalDose("BIC") << "\n";

    field->injectImplantProfile("Interstitial", profile);
    viennaps::DefectClusterKernel<NumericType> loop(1273.15, "loop");
    loop.setPhysicsField(field);
    loop.setMaterialProperties(mats);
    loop.setup();
    loop.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("DislocationLoop") > 0);
    std::cout << "[cluster-check] loop: DislocationLoop="
              << field->getTotalDose("DislocationLoop") << "\n";
  }

  // Sundials with multi-kernel on packed field dofs
  {
    // Use a compact field for CVODE
    auto f2 = std::make_shared<Field>();
    f2->setProfileSize(16);
    f2->injectImplantProfile("Dopant", std::vector<NumericType>(16, NumericType(1e12 / 16)));
    f2->injectImplantProfile("Interstitial",
                             std::vector<NumericType>(16, NumericType(1e12 / 16)));
    f2->injectImplantProfile("Vacancy", std::vector<NumericType>(16, NumericType(5e11 / 16)));

    viennaps::SundialsTimeIntegrator<NumericType> integ;
    integ.setPhysicsField(f2);
    integ.setUseFieldState(true);
    auto kd = std::make_shared<viennaps::DiffusionKernel<NumericType>>("Dopant", 1273.15);
    auto kf = std::make_shared<viennaps::FermiDiffusionKernel<NumericType>>("Interstitial", 1273.15);
    auto kc = std::make_shared<viennaps::ChargedReactKernel<NumericType>>(1273.15);
    auto ks = std::make_shared<viennaps::ViscoelasticStressKernel<NumericType>>(1273.15);
    for (auto k : std::vector<std::shared_ptr<viennaps::PhysicsKernel<NumericType>>>{kd, kf, kc, ks}) {
      k->setPhysicsField(f2);
      k->setMaterialProperties(mats);
      integ.addKernel(k);
    }
    auto doseBefore = f2->getTotalDose("Dopant");
    integ.evolve(0.0, 60.0, 20.0);
    auto doseAfter = f2->getTotalDose("Dopant");
    VC_TEST_ASSERT(doseAfter > 0);
    std::cout << "[cvode-field-check] Dopant dose " << doseBefore << " -> " << doseAfter
              << " stateSize=" << f2->getStateSize() << "\n";
  }

  // --- Oxidation adapter (field-only path for OED / dopant / stress hooks) ---
  {
    viennaps::PhysicsFieldAdapter<NumericType, 2> adapter(field, mats);
    adapter.setOEDDosePerStep(static_cast<NumericType>(1e11));
    adapter.applyToOxidationFieldOnly();
    auto factor = adapter.getDopantEnhancedOxidationFactor();
    auto stress = adapter.getHydrostaticStressForOxidation();
    VC_TEST_ASSERT(factor >= 1);
    std::cout << "[adapter-check] pre: dopantFactor=" << factor << " stress=" << stress
              << "\n";

    auto I_before = field->getTotalDose("Interstitial");
    adapter.updateFromOxidationFieldOnly(static_cast<NumericType>(1.5));
    auto I_after = field->getTotalDose("Interstitial");
    auto oed = field->getTotalDose("OxidationDefects");
    VC_TEST_ASSERT(I_after > I_before);
    VC_TEST_ASSERT(oed > 0);
    VC_TEST_ASSERT(adapter.getLastOEDDose() > 0);
    std::cout << "[adapter-check] OED: I_before=" << I_before << " I_after=" << I_after
              << " OxidationDefects=" << oed
              << " lastOEDDose=" << adapter.getLastOEDDose() << "\n";

    adapter.syncStressToOxidation();
    adapter.syncStressFromOxidation(static_cast<NumericType>(1e6));
  }

  // Orchestrator multi-step
  {
    auto f3 = std::make_shared<Field>();
    f3->setProfileSize(16);
    viennaps::ProcessOrchestrator<NumericType, 2> orch;
    orch.runMultiStepExample(f3, mats, NumericType(100));
    VC_TEST_ASSERT(f3->getTotalDose("Dopant") > 0);
    std::cout << "[orchestrator-check] multi-step OK Dopant=" << f3->getTotalDose("Dopant")
              << " OED-ish I=" << f3->getTotalDose("Interstitial") << "\n";
  }

  std::cout << "[domain test] Track 1 physics smoke passed (field-dofs CVODE + pair/charged + "
               "elastic + geo + clustering + adapter + orchestrator).\n";
  std::cout << "[domain test] Multiphysics validation markers: state-check geo-check "
               "pair-check charged-check elastic-check cluster-check cvode-field-check "
               "adapter-check orchestrator-check OK\n";
}

} // namespace viennacore

int main() {
  std::cout << "Running ViennaPS version: " << viennaps::version << std::endl;
  std::cout << "Major: " << viennaps::versionMajor
            << ", Minor: " << viennaps::versionMinor
            << ", Patch: " << viennaps::versionPatch << std::endl;
  VC_RUN_ALL_TESTS
}
