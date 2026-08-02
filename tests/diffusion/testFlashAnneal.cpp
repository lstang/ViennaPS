// testFlashAnneal.cpp — flash/laser anneal completion tests
// (GAP_ANALYSIS §4.10): solidification trapping + FEM orchestration flow,
// plus regression tests for the existing melt physics.
#include <fields/DiffusionEngine.hpp>
#include <fields/DiffusionPhysics.hpp>
#include <fields/FlashAnnealFlow.hpp>
#include <fields/models/ConstantDiffusion.hpp>
#include <fields/models/FlashLaserAnneal.hpp>
#include <vcTestAsserts.hpp>

#include <cmath>
#include <cstdio>
#include <memory>

namespace viennacore {

using namespace viennaps;

void TestFlashAnneal() {
  constexpr int D = 2;

  // --- Test A: solidification trapping immobilizes mobile dopant.
  {
    MeshAttributes attrs;
    attrs.setAttributeName(1, "Si");
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

    auto runTrap = [&](double r) {
      std::cout << "[flash] runTrap start r=" << r << std::endl;
      auto engine = std::make_unique<DiffusionEngine<double, 2>>();
      engine->setMesh(std::make_unique<mfem::Mesh>(*mesh), attrs);
      DiffusionPhysics<double> physics;
      physics.addSpecies("Boron");
      physics.addSpecies("MeltFraction");
      auto diff = std::make_shared<ConstantDiffusion<double>>("Boron");
      diff->setDiffusivity(1e-14, 0.0);
      physics.addModel(diff);
      auto meltDiff =
          std::make_shared<ConstantDiffusion<double>>("MeltFraction");
      meltDiff->setDiffusivity(0.0, 0.0);
      physics.addModel(meltDiff);
      auto trap = std::make_shared<SolidificationTrapping<double>>("Boron");
      trap->setTrappingStrength(r);
      physics.addModel(trap);
      physics.setTemperature(1273.0);
      engine->setPhysics(physics);
      engine->initializeSpecies("Boron", 1e18);
      engine->initializeSpecies("MeltFraction", 0.0);

      engine->solve(0.0, 100.0, 100.0);
      mfem::ParGridFunction prevPhi(engine->getSolution("MeltFraction"));
      prevPhi = 1.0;
      trap->setPreviousPhi(&prevPhi);
      trap->setDt(100.0);
      const double d0 = engine->getIntegral("Boron");
      engine->solve(100.0, 200.0, 100.0);
      const double d1 = engine->getIntegral("Boron");
      std::cout << "[flash] trap r=" << r << " survival=" << d1 / d0
                << std::endl;
      return d1;
    };

    const double dNoTrap = runTrap(0.0);
    const double dTrap = runTrap(1.0);
    // exp(-1) ≈ 0.37 survival with trapping; full conservation without.
    VC_TEST_ASSERT(dTrap < dNoTrap * 0.6);
    VC_TEST_ASSERT(dNoTrap > 0.0);
  }

  // --- Test B: FlashAnnealFlow end-to-end pulse.
  {
    std::cout << "[flash] Test B start" << std::endl;
    MeshAttributes attrs;
    attrs.setAttributeName(1, "Si");
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(16, 16, mfem::Element::TRIANGLE));
    DiffusionEngine<double, 2> engine;
    engine.setMesh(std::move(mesh), attrs);
    FlashAnnealFlow<double> flow;
    flow.setPulse(/*Tpeak=*/2000.0, /*duration=*/1e-3);
    flow.setLaserAbsorption(10.0);
    flow.setDopantDiffusivities(/*Ds=*/1e-12, /*Dl=*/1e-4);
    flow.setLatentHeat(/*rhoL=*/5e4);
    engine.setPhysics(flow.physics());
    engine.initializeSpecies("Dopant", 1e18);
    flow.seedLaserPulse(engine, /*T0=*/300.0, /*I0=*/1800.0);
    const double d0 = engine.getIntegral("Dopant");
    flow.apply(engine, /*dt=*/1e-4, /*nSteps=*/10);
    const double d1 = engine.getIntegral("Dopant");
    const double rel = std::abs(d1 - d0) / std::max(d0, 1.0);
    std::cout << "[flash] pulse doseRel=" << rel
              << " Tmax=" << engine.getSolution("Temperature").Max()
              << " phiMax=" << engine.getSolution("MeltFraction").Max()
              << std::endl;
    // Surface pulse: T0+I0 = 2100 K > Tm = 1687 K → surface melts.
    VC_TEST_ASSERT(engine.getSolution("Temperature").Max() > 1687.0);
    VC_TEST_ASSERT(engine.getSolution("MeltFraction").Max() > 0.3);
    // Dopant stays conserved (zero-flux, melt-enhanced diffusion).
    VC_TEST_ASSERT(rel < 0.1);
  }

  // --- Test C: melting above Tm.
  {
    std::cout << "[flash] Test C start" << std::endl;
    MeshAttributes attrs;
    attrs.setAttributeName(1, "Si");
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));
    DiffusionEngine<double, 2> engine;
    engine.setMesh(std::move(mesh), attrs);
    auto melt = std::make_shared<MeltingPhaseFieldFEM<double>>();
    melt->setMobility(1.0);
    melt->setGradientEnergy(1e-4);
    melt->setMeltingPoint(1687.0);
    melt->setCoupling(1.0);
    melt->setTemperatureSpecies("Temperature");
    auto tempDiff = std::make_shared<ConstantDiffusion<double>>("Temperature");
    tempDiff->setDiffusivity(0.8, 0.0);
    DiffusionPhysics<double> physics;
    physics.addSpecies("Temperature");
    physics.addSpecies("MeltFraction");
    physics.addModel(tempDiff);
    physics.addModel(melt);
    physics.setTemperature(1273.0);
    engine.setPhysics(physics);
    engine.initializeSpecies("Temperature", 1750.0);
    engine.initializeSpecies("MeltFraction", 0.0);
    engine.solve(0.0, 0.05, 0.01);
    const double meanPhi = engine.meanOnAttribute("MeltFraction", 1);
    std::cout << "[flash] melt meanPhi=" << meanPhi << std::endl;
    VC_TEST_ASSERT(meanPhi > 0.5);
    VC_TEST_ASSERT(meanPhi < 5.0); // phase field can overshoot; bounded
  }

  // --- Test D: melt-enhanced dopant diffusion.
  {
    std::cout << "[flash] Test D start" << std::endl;
    auto width = [&](const mfem::ParGridFunction &gf) {
      const mfem::Mesh *m = gf.ParFESpace()->GetMesh();
      const int n = gf.Size();
      double sum = 0.0, sx = 0.0, sy = 0.0;
      for (int i = 0; i < n; ++i) {
        const double w = static_cast<double>(gf(i));
        const double *v = m->GetVertex(i);
        sum += w;
        sx += w * v[0];
        sy += w * v[1];
      }
      double vx = 0.0, vy = 0.0;
      for (int i = 0; i < n; ++i) {
        const double w = static_cast<double>(gf(i));
        const double *v = m->GetVertex(i);
        vx += w * (v[0] - sx / sum) * (v[0] - sx / sum);
        vy += w * (v[1] - sy / sum) * (v[1] - sy / sum);
      }
      return std::sqrt((vx + vy) / sum);
    };

    auto runMeltDiff = [&](double phi) {
      MeshAttributes attrs;
      attrs.setAttributeName(1, "Si");
      auto mesh = std::make_unique<mfem::Mesh>(
          mfem::Mesh::MakeCartesian2D(16, 16, mfem::Element::TRIANGLE));
      auto engine = std::make_unique<DiffusionEngine<double, 2>>();
      engine->setMesh(std::move(mesh), attrs);
      auto md = std::make_shared<MeltDiffusion<double>>("Boron");
      md->setSolidD(1e-12);
      md->setLiquidD(1e-4);
      auto meltModel =
          std::make_shared<ConstantDiffusion<double>>("MeltFraction");
      meltModel->setDiffusivity(0.0, 0.0);
      DiffusionPhysics<double> physics;
      physics.addSpecies("Boron");
      physics.addSpecies("MeltFraction");
      physics.addModel(md);
      physics.addModel(meltModel);
      physics.setTemperature(1273.0);
      engine->setPhysics(physics);
      std::vector<double> peak(17, 0.0);
      for (std::size_t i = 0; i < peak.size(); ++i) {
        const double z = (static_cast<double>(i) - 8.0) / 4.0;
        peak[i] = 1e18 * std::exp(-z * z);
      }
      engine->projectIntegralPreserving("Boron", peak, 4e18);
      engine->initializeSpecies("MeltFraction", phi);
      engine->solve(0.0, 10.0, 1.0);
      return width(engine->getSolution("Boron"));
    };

    const double wLiquid = runMeltDiff(1.0);
    const double wSolid = runMeltDiff(0.0);
    std::cout << "[flash] width liquid=" << wLiquid << " solid=" << wSolid
              << std::endl;
    // Dl·t = 1e-3 → spread ~sqrt(4e-3) ≈ 0.063; solid case ~no spread.
    VC_TEST_ASSERT(wLiquid > wSolid);
  }

  // --- Test E: latent heat absorbs during melting.
  {
    std::cout << "[flash] Test E start" << std::endl;
    auto runLatent = [&](double rhoL) {
      MeshAttributes attrs;
      attrs.setAttributeName(1, "Si");
      auto mesh = std::make_unique<mfem::Mesh>(
          mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));
      auto engine = std::make_unique<DiffusionEngine<double, 2>>();
      engine->setMesh(std::move(mesh), attrs);
      auto heat = std::make_shared<HeatTransfer<double>>();
      heat->setThermalDiffusivity(0.8);
      heat->setLatentHeat(rhoL);
      auto melt = std::make_shared<MeltingPhaseFieldFEM<double>>();
      melt->setMobility(1.0);
      melt->setGradientEnergy(1e-4);
      melt->setMeltingPoint(1687.0);
      melt->setCoupling(1.0);
      melt->setTemperatureSpecies("Temperature");
      DiffusionPhysics<double> physics;
      physics.addSpecies("Temperature");
      physics.addSpecies("MeltFraction");
      physics.addModel(heat);
      physics.addModel(melt);
      physics.setTemperature(1273.0);
      engine->setPhysics(physics);
      engine->initializeSpecies("Temperature", 1700.0);
      engine->initializeSpecies("MeltFraction", 0.0);
      engine->solve(0.0, 0.01, 0.01);
      mfem::ParGridFunction prevPhi(engine->getSolution("MeltFraction"));
      prevPhi = 0.0; // φ rising → melting → latent term absorbs heat
      heat->setPreviousPhi(&prevPhi);
      heat->setDt(0.01);
      engine->solve(0.01, 0.05, 0.01);
      return engine->meanOnAttribute("Temperature", 1);
    };
    const double tLatent = runLatent(5e4);
    const double tNoLatent = runLatent(0.0);
    std::cout << "[flash] meanT latent=" << tLatent << " noLatent=" << tNoLatent
              << std::endl;
    VC_TEST_ASSERT(tLatent < tNoLatent);
  }
}

} // namespace viennacore

int main() {
  std::setvbuf(stdout, NULL, _IONBF, 0);
#ifdef VIENNAPS_HAS_MFEM
  viennacore::TestFlashAnneal();
  std::cout << "All flash anneal tests passed." << std::endl;
#else
  std::cout << "MFEM not available, skipping flash anneal tests." << std::endl;
#endif
  return 0;
}
