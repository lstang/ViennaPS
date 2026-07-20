#define VIENNAPS_HAS_SUNDIALS 1
#include <vcTestAsserts.hpp>

// Multipysics Track 1 / Phase 1 smoke without full ViennaLS geometry includes.
#include "psVersion.hpp"
#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"
#include "fields/DiffusionKernel.hpp"
#include "fields/FermiDiffusionKernel.hpp"
#include "fields/StressKernel.hpp"
#include "fields/DefectClusterKernel.hpp"
#include "fields/SundialsTimeIntegrator.hpp"
#include "process/psPhysicsFieldAdapter.hpp"

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
  auto mats = std::make_shared<MatSys>();
  mats->setArrhenius("Si", "Dopant_D", 0.1, 2.5);
  mats->setArrhenius("Si", "Interstitial_D", 10.0, 1.8);
  mats->setProperty("Si", "YoungModulus", 130.0);
  mats->setProperty("Si", "PoissonRatio", 0.28);
  mats->setProperty("Si", "GrowthStress", 300.0e6);

  // Implant-like seed
  std::vector<NumericType> profile(64, static_cast<NumericType>(1e12 / 64));
  field->injectImplantProfile("Dopant", profile);
  field->injectImplantProfile("Interstitial", profile);
  field->injectImplantProfile("Vacancy",
                              std::vector<NumericType>(64, static_cast<NumericType>(5e11 / 64)));
  VC_TEST_ASSERT(field->getTotalDose("Dopant") > 0);

  // Diffusion + Fermi + stress
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
    viennaps::ViscoelasticStressKernel<NumericType> k(1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }

  // --- Clustering models (recomb / 311 / bic / loop) ---
  {
    auto I0 = field->getTotalDose("Interstitial");
    auto V0 = field->getTotalDose("Vacancy");

    viennaps::DefectClusterKernel<NumericType> recomb(1273.15, "recomb");
    recomb.setPhysicsField(field);
    recomb.setMaterialProperties(mats);
    recomb.setup();
    recomb.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("RecombinedIV") > 0);
    std::cout << "[cluster-check] recomb: I0=" << I0 << " V0=" << V0
              << " RecombinedIV=" << field->getTotalDose("RecombinedIV") << "\n";

    // Reset interstitial for 311
    field->injectImplantProfile("Interstitial", profile);
    viennaps::DefectClusterKernel<NumericType> c311(1273.15, "311");
    c311.setPhysicsField(field);
    c311.setMaterialProperties(mats);
    c311.setup();
    c311.evolve(10.0);
    VC_TEST_ASSERT(field->getTotalDose("Cluster311") > 0);
    std::cout << "[cluster-check] 311: Cluster311=" << field->getTotalDose("Cluster311")
              << "\n";

    // BIC needs dopant + I
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

  // Sundials with multi-kernel (CVODE collects RHS via addToRHS)
  {
    viennaps::SundialsTimeIntegrator<NumericType> integ;
    integ.setPhysicsField(field);
    auto kd = std::make_shared<viennaps::DiffusionKernel<NumericType>>("Dopant", 1273.15);
    auto kf = std::make_shared<viennaps::FermiDiffusionKernel<NumericType>>("Interstitial", 1273.15);
    auto ks = std::make_shared<viennaps::ViscoelasticStressKernel<NumericType>>(1273.15);
    kd->setPhysicsField(field);
    kd->setMaterialProperties(mats);
    kf->setPhysicsField(field);
    kf->setMaterialProperties(mats);
    ks->setPhysicsField(field);
    ks->setMaterialProperties(mats);
    integ.addKernel(kd);
    integ.addKernel(kf);
    integ.addKernel(ks);
    integ.evolve(0.0, 90.0, 30.0);
  }

  // --- Oxidation adapter (field-only path for OED / dopant / stress hooks) ---
  {
    // Use D=2 for adapter template; field-only APIs do not need real Domain
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

  std::cout << "[domain test] Track 1 physics smoke passed (kernels + stress coupling + "
               "clustering + Sundials + Field).\n";
  std::cout << "[domain test] Multiphysics validation markers: cluster-check adapter-check "
               "OK\n";
}

} // namespace viennacore

int main() {
  std::cout << "Running ViennaPS version: " << viennaps::version << std::endl;
  std::cout << "Major: " << viennaps::versionMajor
            << ", Minor: " << viennaps::versionMinor
            << ", Patch: " << viennaps::versionPatch << std::endl;
  VC_RUN_ALL_TESTS
}
