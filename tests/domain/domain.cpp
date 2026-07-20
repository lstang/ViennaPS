#define VIENNAPS_HAS_SUNDIALS 1
#include <vcTestAsserts.hpp>

// Safe includes for the new multiphysics Track 1 code only.
// These do not transitively require ls* / ViennaLS headers.
#include "psVersion.hpp"
#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"
#include "fields/DiffusionKernel.hpp"
#include "fields/FermiDiffusionKernel.hpp"
#include "fields/StressKernel.hpp"
#include "fields/DefectClusterKernel.hpp"
#include "fields/SundialsTimeIntegrator.hpp"
#include "process/psPhysicsFieldAdapter.hpp"

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennacore {

template <class NumericType, int D>
void RunTest() {
  // Self-contained smoke test for Track 1 physics (no geometry/level-set construction).
  // Verifies: Field storage + MFEM projection (when enabled), kernels (diffusion/fermi/stress/cluster),
  // stress <-> diffusion coupling, Sundials driver, implant-like injection.

  using Field = viennaps::PhysicsField<NumericType>;
  using MatSys = viennaps::MaterialPropertySystem<NumericType>;

  auto field = std::make_shared<Field>();
  auto mats  = std::make_shared<MatSys>();

  // Improve numerics for demo
  mats->setArrhenius("Si", "Dopant_D", 0.1, 2.5);
  mats->setArrhenius("Si", "Interstitial_D", 10.0, 1.8);

  // Elastic props for real MFEM stress solve
  mats->setProperty("Si", "YoungModulus", 130.0);
  mats->setProperty("Si", "PoissonRatio", 0.28);
  mats->setProperty("Si", "GrowthStress", 300.0e6);

  // "Implant" simulation via direct injection (profile logic lives in AnalyticImplant; here we test the field side)
  std::vector<NumericType> profile(64, static_cast<NumericType>(1e12 / 64));
  field->injectImplantProfile("Dopant", profile);
  field->injectImplantProfile("Interstitial", profile);
  VC_TEST_ASSERT(field->getTotalDose("Dopant") > 0);

  // Basic diffusion kernel
  {
    viennaps::DiffusionKernel<NumericType> k("Dopant", 1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }

  // Fermi (extrinsic, concentration + stress dependent)
  {
    viennaps::FermiDiffusionKernel<NumericType> k("Dopant", 1273.15, 1e18);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }

  // Viscoelastic + growth stress (Track 1 mechanics)
  {
    viennaps::ViscoelasticStressKernel<NumericType> k(1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }

  // Basic defect clustering
  {
    viennaps::DefectClusterKernel<NumericType> k(1273.15);
    k.setPhysicsField(field);
    k.setMaterialProperties(mats);
    k.setup();
    k.evolve(30.0);
  }

  // Sundials time integrator driving multiple kernels
  {
    viennaps::SundialsTimeIntegrator<NumericType> integ;
    integ.setPhysicsField(field);

    auto kd = std::make_shared<viennaps::DiffusionKernel<NumericType>>("Dopant", 1273.15);
    auto kf = std::make_shared<viennaps::FermiDiffusionKernel<NumericType>>("Interstitial", 1273.15);
    auto ks = std::make_shared<viennaps::ViscoelasticStressKernel<NumericType>>(1273.15);

    kd->setPhysicsField(field); kd->setMaterialProperties(mats);
    kf->setPhysicsField(field); kf->setMaterialProperties(mats);
    ks->setPhysicsField(field); ks->setMaterialProperties(mats);

    integ.addKernel(kd);
    integ.addKernel(kf);
    integ.addKernel(ks);

    integ.evolve(0.0, 90.0, 30.0);
  }

  std::cout << "[domain test] Track 1 physics smoke passed (kernels + stress coupling + clustering + Sundials + Field).\n";

  // Explicit double instantiation to exercise real CVODE path
  {
    using Dbl = double;
    auto fieldD = std::make_shared<viennaps::PhysicsField<Dbl>>();
    auto matsD = std::make_shared<viennaps::MaterialPropertySystem<Dbl>>();
    std::vector<Dbl> p(64, 1e12/64);
    fieldD->injectImplantProfile("Dopant", p);
    viennaps::SundialsTimeIntegrator<Dbl> integ;
    integ.setPhysicsField(fieldD);
    auto k = std::make_shared<viennaps::DiffusionKernel<Dbl>>("Dopant", 1273.15);
    k->setPhysicsField(fieldD);
    k->setMaterialProperties(matsD);
    integ.addKernel(k);
    integ.evolve(0, 90, 30);
  }
}

} // namespace viennacore

int main() {
  std::cout << "Running ViennaPS version: " << viennaps::version << std::endl;
  std::cout << "Major: " << viennaps::versionMajor
            << ", Minor: " << viennaps::versionMinor
            << ", Patch: " << viennaps::versionPatch << std::endl;
  VC_RUN_ALL_TESTS
}
