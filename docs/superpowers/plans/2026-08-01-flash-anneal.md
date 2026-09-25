# Flash/Laser Anneal Completion Implementation Plan

> For agentic workers: REQUIRED SUB-SKILL — `test-driven-development` for every task
> (write the failing test first, verify it fails, implement, verify it passes, commit);
> `verification-before-completion` before claiming any task done; execute this plan
> via `executing-plans` (inline) or `subagent-driven-development` (subagents), one
> task at a time. Do not skip the format step (Task 4).

## Goal

Close the medium-term gap "flash/laser anneal completion" from
`docs/GAP_ANALYSIS.md` §4.10: the physics skeleton exists
(`HeatTransfer`, `LaserIntensity`, `MeltingPhaseFieldFEM`,
`CrystallinityPhaseFieldFEM`, `MeltDiffusion`, 1D `FlashLaserAnneal::runPulse`),
but there is **no process-level FEM orchestration** — no multi-species
`DiffusionEngine` flow that couples heat → melting (latent heat ρ·L·∂φ/∂t) →
melt-enhanced dopant diffusion → resolidification trapping. This plan adds:

1. `SolidificationTrapping` — a dopant reaction model: at the resolidification
   front (φ falling) mobile dopant is immobilized into the regrown lattice
   (`R = −r·max(0,−∂φ/∂t)·C`, ∂φ/∂t via the previous-φ copy pattern already
   used by `HeatTransfer`).
2. `FlashAnnealFlow` — an orchestrator that registers
   `{Temperature, MeltFraction, Crystallinity, Dopant}` + the existing models
   into one `DiffusionPhysics`, seeds a Beer's-law surface temperature pulse,
   and runs sequential `DiffusionEngine::solve` segments, refreshing the
   previous-φ copy between segments for the latent-heat/trapping terms.
3. Regression tests that lock the existing physics (melting above Tm,
   melt-enhanced diffusion, latent-heat absorption) plus red-green tests for
   the two new pieces.

## Architecture

```
fields/models/FlashLaserAnneal.hpp (edit)
  + SolidificationTrapping<NumericType> : DiffusionModel<NumericType>
      species: {Dopant}; reads MeltFraction + previous-φ copy (setPreviousPhi/setDt)
      assembleReaction: R = -r * max(0, -(φ_cur - φ_prev)/dt) * C

fields/FlashAnnealFlow.hpp (NEW)
  FlashAnnealFlow<NumericType>
    physics(): DiffusionPhysics<NumericType>&  — {Temperature, MeltFraction,
               Crystallinity, Dopant} + HeatTransfer + MeltingPhaseFieldFEM +
               MeltDiffusion + CrystallinityPhaseFieldFEM + SolidificationTrapping
    setPulse(Tpeak, duration) / setLaserAbsorption(alpha) / setLatentHeat(rhoL) /
    setDopantDiffusivities(Ds, Dl) / setTrappingStrength(r) / setMeltParameters(...)
    seedLaserPulse(engine, T0, I0)  — T(z) = T0 + I0·exp(-αz) via
               projectIntegralPreserving (unit-area domain; dose = mean)
    apply(engine, dt, nSteps)       — per segment: refresh prevPhi_ ← MeltFraction,
               setPreviousPhi/setDt on heat+trap models, engine.solve(t, t+dt, dt)
```

## Tech Stack

- C++20, header-only; MFEM-backed `DiffusionEngine<double, 2>`; tests in
  `tests/diffusion/` (same conventions as `testDiffusion.cpp`).

## Global Constraints

- C++20; LLVM clang-format (2-space, 80-col).
- `FlashAnnealFlow` lives in `fields/` (it includes `DiffusionEngine.hpp`);
  `SolidificationTrapping` lives in `fields/models/FlashLaserAnneal.hpp` next
  to the other melt models (it is a pure `DiffusionModel`). The existing 1D
  `FlashLaserAnneal::runPulse` API is untouched.
- The `setPreviousPhi(const mfem::ParGridFunction*)` pattern is copied
  verbatim from `HeatTransfer` (raw pointer, orchestrator-managed lifetime —
  the flow owns the copy).
- `prevPhi` copies: allocate with `mfem::ParGridFunction(engine.fes())` then
  copy-assign from `engine.getSolution("MeltFraction")` (same space). A
  default-constructed `ParGridFunction` cannot be copy-assigned from a
  different-size field.
- Tests: new `tests/diffusion/testFlashAnneal.cpp` + CMake registration;
  build commands per AGENTS.md; commit style `feat(fields): ...`.

---

## Task 1 — Failing tests: trapping + end-to-end flow

**Files**
- `tests/diffusion/testFlashAnneal.cpp` (new)
- `tests/diffusion/CMakeLists.txt` (edit)

**Interfaces**
- Consumes (new, test-only): `SolidificationTrapping` (ctor `(species)` +
  `setTrappingStrength(r)` + `setMeltSpecies` + `setPreviousPhi` + `setDt`),
  `FlashAnnealFlow` (`physics()`, `setPulse`, `setDopantDiffusivities`,
  `setLatentHeat`, `seedLaserPulse`, `apply`).
- Consumes (existing): `MeltingPhaseFieldFEM`, `MeltDiffusion`, `HeatTransfer`
  (as in `FlashLaserAnneal.hpp`), `DiffusionEngine` ICs
  (`initializeSpecies`, `projectIntegralPreserving`), `engine.fes()`,
  `engine.getSolution(name).Max()`, `meanOnAttribute`, `getIntegral`.
- Produces: `tests/diffusion/testFlashAnneal` ctest target.

**Steps**

1. Create `tests/diffusion/testFlashAnneal.cpp`:

```cpp
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
#include <memory>

namespace viennacore {

using namespace viennaps;

template <class NumericType> void RunTest() {
  constexpr int D = 2;

  // --- Test A (red before Task 2): solidification trapping immobilizes
  // mobile dopant. With previous φ = 1 and current φ = 0, ∂φ/∂t < 0 and
  // the trapping reaction removes C at rate r·(1/dt)·C; over one segment
  // of duration dt the surviving fraction is ~exp(-r).
  {
    MeshAttributes attrs;
    attrs.setAttributeName(1, "Si");
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

    auto runTrap = [&](double r) {
      auto engine = std::make_unique<DiffusionEngine<double, 2>>();
      engine->setMesh(std::make_unique<mfem::Mesh>(*mesh), attrs);
      DiffusionPhysics<double> physics;
      physics.addSpecies("Boron");
      physics.addSpecies("MeltFraction");
      auto diff = std::make_shared<ConstantDiffusion<double>>("Boron");
      diff->setDiffusivity(1e-14, 0.0);
      physics.addModel(diff);
      auto trap = std::make_shared<SolidificationTrapping<double>>("Boron");
      trap->setTrappingStrength(r);
      physics.addModel(trap);
      physics.setTemperature(1273.0);
      engine->setPhysics(physics);
      engine->initializeSpecies("Boron", 1e18);
      engine->initializeSpecies("MeltFraction", 0.0);
      // Previous-φ field: fully molten before the segment.
      mfem::ParGridFunction prevPhi(engine->fes());
      prevPhi = 1.0;
      trap->setPreviousPhi(&prevPhi);
      trap->setDt(100.0);
      const double d0 = engine->getIntegral("Boron");
      engine->solve(0.0, 100.0, 100.0);
      const double d1 = engine->getIntegral("Boron");
      std::cout << "[flash] trap r=" << r << " survival=" << d1 / d0 << "\n";
      return d1;
    };

    const double dNoTrap = runTrap(0.0);
    const double dTrap = runTrap(1.0);
    // exp(-1) ≈ 0.37 survival with trapping; full conservation without.
    VC_TEST_ASSERT(dTrap < dNoTrap * 0.6);
    VC_TEST_ASSERT(dNoTrap > 0.0);
  }

  // --- Test B (red before Task 3): FlashAnnealFlow end-to-end pulse.
  {
    MeshAttributes attrs;
    attrs.setAttributeName(1, "Si");
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(16, 16, mfem::Element::TRIANGLE));
    DiffusionEngine<double, 2> engine;
    engine.setMesh(std::move(mesh), attrs);
    FlashAnnealFlow<double> flow;
    flow.setPulse(/*Tpeak=*/2000.0, /*duration=*/1e-3);
    flow.setDopantDiffusivities(/*Ds=*/1e-12, /*Dl=*/1e-4);
    flow.setLatentHeat(/*rhoL=*/5e4);
    engine.setPhysics(flow.physics());
    engine.initializeSpecies("Dopant", 1e18);
    flow.seedLaserPulse(engine, /*T0=*/300.0, /*I0=*/1400.0);
    const double d0 = engine.getIntegral("Dopant");
    flow.apply(engine, /*dt=*/1e-4, /*nSteps=*/10);
    const double d1 = engine.getIntegral("Dopant");
    const double rel = std::abs(d1 - d0) / std::max(d0, 1.0);
    std::cout << "[flash] pulse doseRel=" << rel
              << " Tmax=" << engine.getSolution("Temperature").Max()
              << " phiMax=" << engine.getSolution("MeltFraction").Max()
              << "\n";
    // Surface pulse: T0+I0 = 1700 K > Tm = 1687 K → surface melts.
    VC_TEST_ASSERT(engine.getSolution("Temperature").Max() > 1687.0);
    VC_TEST_ASSERT(engine.getSolution("MeltFraction").Max() > 0.5);
    // Dopant stays conserved (zero-flux, melt-enhanced diffusion).
    VC_TEST_ASSERT(rel < 0.1);
  }

  // --- Test C (regression, green already): melting above Tm.
  {
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
    DiffusionPhysics<double> physics;
    physics.addSpecies("Temperature");
    physics.addSpecies("MeltFraction");
    physics.addModel(melt);
    physics.setTemperature(1273.0);
    engine.setPhysics(physics);
    engine.initializeSpecies("Temperature", 1750.0);
    engine.initializeSpecies("MeltFraction", 0.0);
    engine.solve(0.0, 0.05, 0.01);
    const double meanPhi = engine.meanOnAttribute("MeltFraction", 1);
    std::cout << "[flash] melt meanPhi=" << meanPhi << "\n";
    VC_TEST_ASSERT(meanPhi > 0.5);
    VC_TEST_ASSERT(meanPhi < 2.0); // phase field can overshoot; bounded
  }

  // --- Test D (regression, green already): melt-enhanced dopant diffusion.
  {
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
      DiffusionPhysics<double> physics;
      physics.addSpecies("Boron");
      physics.addSpecies("MeltFraction");
      physics.addModel(md);
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
      return width(*engine->getSolution("Boron"));
    };

    const double wLiquid = runMeltDiff(1.0);
    const double wSolid = runMeltDiff(0.0);
    std::cout << "[flash] width liquid=" << wLiquid << " solid=" << wSolid
              << "\n";
    // Dl·t = 1e-3 → spread ~sqrt(4e-3) ≈ 0.063; solid case ~no spread.
    VC_TEST_ASSERT(wLiquid > wSolid * 2.0);
  }

  // --- Test E (regression, green already): latent heat absorbs during
  // melting — with ρ·L > 0 the final mean temperature is lower.
  {
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
      mfem::ParGridFunction prevPhi(engine->fes());
      prevPhi = 0.0; // φ rising → melting → latent term absorbs heat
      heat->setPreviousPhi(&prevPhi);
      heat->setDt(0.01);
      engine->solve(0.0, 0.05, 0.01);
      return engine->meanOnAttribute("Temperature", 1);
    };
    const double tLatent = runLatent(5e4);
    const double tNoLatent = runLatent(0.0);
    std::cout << "[flash] meanT latent=" << tLatent
              << " noLatent=" << tNoLatent << "\n";
    VC_TEST_ASSERT(tLatent < tNoLatent);
  }
}

} // namespace viennacore

int main() { VC_RUN_ALL_TESTS }
```

2. Append to `tests/diffusion/CMakeLists.txt`:

```cmake
project(testFlashAnneal LANGUAGES CXX)
viennaps_add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
add_dependencies(ViennaPS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
```

3. Build and run — expect **FAIL** (red): `SolidificationTrapping` and
   `FlashAnnealFlow` do not exist. Tests C/D/E (existing physics) are
   expected to pass once the file compiles.

```powershell
cmake --build build --config Release --target testFlashAnneal --parallel
# expected: error C2039: 'SolidificationTrapping': is not a member of 'viennaps'
#           error C2065: 'FlashAnnealFlow': undeclared identifier
```

**Commit**: none yet.

## Task 2 — Implement SolidificationTrapping

**Files**
- `include/viennaps/fields/models/FlashLaserAnneal.hpp` (edit — add the class
  before `FlashLaserAnneal`)

**Interfaces**
- Produces (namespace `viennaps`):
  - `SolidificationTrapping<NumericType> : DiffusionModel<NumericType>` —
    ctor `explicit SolidificationTrapping(std::string species = "Boron")`;
    `setTrappingStrength(NumericType r)` (0 = disabled);
    `setMeltSpecies(std::string)` (default `"MeltFraction"`);
    `setPreviousPhi(const mfem::ParGridFunction *prev)`;
    `setDt(NumericType dt)`; `numSpecies() = 1`;
    `speciesNames() = {species}`; MFEM `assembleReaction` adding
    `R = −r·max(0,−(φ_cur−φ_prev)/dt)·C` as a `DomainLFIntegrator`.
- Consumes: the `DiffusionModel` base hooks, the `allSpecies` map (current
  MeltFraction), the previous-φ copy.

**Steps**

1. Add to `FlashLaserAnneal.hpp`, after `CrystallinityPhaseFieldFEM`:

```cpp
/// Solidification trapping: at the resolidification front (φ falling),
/// mobile dopant is immobilized into the regrown lattice.
/// R = -r * max(0, -∂φ/∂t) * C, with ∂φ/∂t = (φ_cur - φ_prev)/dt read
/// from the registered MeltFraction species and the orchestrator-supplied
/// previous-φ copy (same pattern as HeatTransfer::LatentHeatCoef).
/// Approximate model (documented in GAP_ANALYSIS §4.10): the immobilized
/// fraction leaves the mobile dopant field; r = 0 disables trapping.
template <class NumericType>
class SolidificationTrapping : public DiffusionModel<NumericType> {
public:
  explicit SolidificationTrapping(std::string species = "Boron")
      : species_(std::move(species)) {
    this->setName("SolidificationTrapping(" + species_ + ")");
  }
  void setTrappingStrength(NumericType r) { r_ = r; }
  void setMeltSpecies(std::string m) { melt_ = std::move(m); }
  void setPreviousPhi(const mfem::ParGridFunction *prev) {
    previousPhi_ = prev;
  }
  void setDt(NumericType dt) { dt_ = dt; }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    if (r_ == NumericType(0) || !previousPhi_ || dt_ <= NumericType(0))
      return;
    const mfem::ParGridFunction *phiCurr = nullptr;
    auto it = allSpecies.find(melt_);
    if (it != allSpecies.end())
      phiCurr = it->second;
    trapCoef_ = std::make_unique<TrappingCoef>(
        &speciesGF, phiCurr, previousPhi_, static_cast<double>(r_),
        static_cast<double>(dt_));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*trapCoef_));
  }
#endif

private:
  std::string species_;
  std::string melt_ = "MeltFraction";
  NumericType r_ = NumericType(0);
  NumericType dt_ = NumericType(0);
  const mfem::ParGridFunction *previousPhi_ = nullptr;
#ifdef VIENNAPS_HAS_MFEM
  class TrappingCoef : public mfem::Coefficient {
  public:
    TrappingCoef(const mfem::ParGridFunction *C,
                 const mfem::ParGridFunction *phiCurr,
                 const mfem::ParGridFunction *phiPrev, double r, double dt)
        : C_(C), phiCurr_(phiCurr), phiPrev_(phiPrev), r_(r), dt_(dt) {}
    double Eval(mfem::ElementTransformation &tr,
                const mfem::IntegrationPoint &ip) override {
      if (!C_ || !phiCurr_ || !phiPrev_ || dt_ <= 0.0)
        return 0.0;
      const double c = std::max(0.0, C_->GetValue(tr, ip));
      const double cur = phiCurr_->GetValue(tr, ip);
      const double prev = phiPrev_->GetValue(tr, ip);
      const double dphiDt = (cur - prev) / dt_;
      // Only solidifying cells (∂φ/∂t < 0) trap; strength r.
      return -r_ * std::max(0.0, -dphiDt) * c;
    }

  private:
    const mfem::ParGridFunction *C_;
    const mfem::ParGridFunction *phiCurr_;
    const mfem::ParGridFunction *phiPrev_;
    double r_, dt_;
  };
  mutable std::unique_ptr<TrappingCoef> trapCoef_;
#endif
};
```

2. Build and run — expect Test A **PASS** (green); Test B still **FAIL**
   (FlashAnnealFlow missing):

```powershell
cmake --build build --config Release --target testFlashAnneal --parallel
ctest -R testFlashAnneal --test-dir build -C Release --output-on-failure
# expected: "[flash] trap r=1 survival≈0.37" ; Test B fails to compile
```

**Commit**: `feat(fields): solidification dopant trapping model (SolidificationTrapping)`

## Task 3 — Implement FlashAnnealFlow

**Files**
- `include/viennaps/fields/FlashAnnealFlow.hpp` (new)

**Interfaces**
- Produces (namespace `viennaps`):
  - `FlashAnnealFlow<NumericType>` — `physics()` →
    `DiffusionPhysics<NumericType>&`; setters `setPulse(Tpeak, duration)`,
    `setLaserAbsorption(alpha)`, `setLatentHeat(rhoL)`,
    `setDopantDiffusivities(Ds, Dl)`, `setTrappingStrength(r)`,
    `setMeltParameters(Tm, L, kappa, lambda)`;
    `void seedLaserPulse(DiffusionEngine<NumericType, 2> &engine,
    NumericType T0, NumericType I0) const`;
    `void apply(DiffusionEngine<NumericType, 2> &engine, NumericType dt,
    int nSteps)`.
- Consumes: `DiffusionEngine` (setPhysics/solve/initializeSpecies/
  projectIntegralPreserving/getSolution/fes), all models in
  `FlashLaserAnneal.hpp`.

**Steps**

1. Create `include/viennaps/fields/FlashAnnealFlow.hpp`:

```cpp
#pragma once

/// FlashAnnealFlow — end-to-end flash/laser anneal orchestration on
/// DiffusionEngine (GAP_ANALYSIS §4.10): surface heat pulse → melting
/// phase field (latent heat ρ·L·∂φ/∂t coupled) → melt-enhanced dopant
/// diffusion + solidification trapping + SPER crystallinity. All species
/// and models run in ONE multi-species solve per segment; between
/// segments the MeltFraction copy is refreshed for the ∂φ/∂t terms.

#include "DiffusionEngine.hpp"
#include "DiffusionPhysics.hpp"
#include "models/FlashLaserAnneal.hpp"

#include <memory>

namespace viennaps {

template <class NumericType> class FlashAnnealFlow {
public:
  FlashAnnealFlow() { buildPhysics(); }

  // ---- Configuration ------------------------------------------------------
  void setPulse(NumericType Tpeak, NumericType duration) {
    Tpeak_ = Tpeak;
    duration_ = duration;
  }
  void setLaserAbsorption(NumericType alpha) { laserAlpha_ = alpha; }
  void setLatentHeat(NumericType rhoL) {
    heat_->setLatentHeat(rhoL);
  }
  void setDopantDiffusivities(NumericType Ds, NumericType Dl) {
    meltDiff_->setSolidD(Ds);
    meltDiff_->setLiquidD(Dl);
  }
  void setTrappingStrength(NumericType r) {
    trap_->setTrappingStrength(r);
  }
  void setMeltParameters(NumericType Tm, NumericType L, NumericType kappa,
                         NumericType lambda) {
    meltFem_->setMeltingPoint(Tm);
    meltFem_->setMobility(L);
    meltFem_->setGradientEnergy(kappa);
    meltFem_->setCoupling(lambda);
  }

  /// Seed the Temperature species with a Beer's-law surface pulse:
  /// T(z) = T0 + I0 * exp(-alpha * z) along the dof axis, projected with
  /// projectIntegralPreserving and rescaled to the analytical mean
  /// (assumes a unit-area domain, matching the test meshes).
  void seedLaserPulse(DiffusionEngine<NumericType, 2> &engine,
                      NumericType T0, NumericType I0) const {
    const int n = 64;
    std::vector<NumericType> samples(static_cast<std::size_t>(n),
                                     NumericType(0));
    for (int i = 0; i < n; ++i) {
      const double z = static_cast<double>(i) / static_cast<double>(n - 1);
      samples[static_cast<std::size_t>(i)] =
          T0 + I0 * std::exp(-static_cast<double>(laserAlpha_) * z);
    }
    // Analytical integral of T(z) over the unit depth (per unit area).
    const double alpha = static_cast<double>(laserAlpha_);
    const NumericType dose = static_cast<NumericType>(
        static_cast<double>(T0) + I0 * (1.0 - std::exp(-alpha)) / alpha);
    engine.projectIntegralPreserving("Temperature", samples, dose);
  }

  /// Run `nSteps` segments of `dt` each. Between segments, the
  /// previous-φ copy is refreshed and handed to the latent-heat and
  /// trapping models via setPreviousPhi/setDt.
  void apply(DiffusionEngine<NumericType, 2> &engine, NumericType dt,
             int nSteps) {
    NumericType t = NumericType(0);
    for (int s = 0; s < nSteps; ++s) {
      if (!prevPhi_) {
        prevPhi_ = std::make_unique<mfem::ParGridFunction>(engine.fes());
      }
      *prevPhi_ = engine.getSolution("MeltFraction");
      heat_->setPreviousPhi(prevPhi_.get());
      heat_->setDt(dt);
      trap_->setPreviousPhi(prevPhi_.get());
      trap_->setDt(dt);
      engine.solve(t, t + dt, dt);
      t += dt;
    }
  }

  DiffusionPhysics<NumericType> &physics() { return physics_; }

private:
  void buildPhysics() {
    physics_.addSpecies("Temperature");
    physics_.addSpecies("MeltFraction");
    physics_.addSpecies("Crystallinity");
    physics_.addSpecies("Dopant");

    heat_ = std::make_shared<HeatTransfer<NumericType>>();
    heat_->setThermalDiffusivity(NumericType(0.8));
    meltFem_ = std::make_shared<MeltingPhaseFieldFEM<NumericType>>();
    meltFem_->setMeltingPoint(NumericType(1687));
    meltFem_->setTemperatureSpecies("Temperature");
    meltDiff_ = std::make_shared<MeltDiffusion<NumericType>>("Dopant");
    cryst_ = std::make_shared<CrystallinityPhaseFieldFEM<NumericType>>();
    cryst_->setTemperatureSpecies("Temperature");
    trap_ = std::make_shared<SolidificationTrapping<NumericType>>("Dopant");

    physics_.addModel(heat_);
    physics_.addModel(meltFem_);
    physics_.addModel(meltDiff_);
    physics_.addModel(cryst_);
    physics_.addModel(trap_);
    physics_.setTemperature(NumericType(300));
  }

  DiffusionPhysics<NumericType> physics_;
  std::shared_ptr<HeatTransfer<NumericType>> heat_;
  std::shared_ptr<MeltingPhaseFieldFEM<NumericType>> meltFem_;
  std::shared_ptr<MeltDiffusion<NumericType>> meltDiff_;
  std::shared_ptr<CrystallinityPhaseFieldFEM<NumericType>> cryst_;
  std::shared_ptr<SolidificationTrapping<NumericType>> trap_;
  std::unique_ptr<mfem::ParGridFunction> prevPhi_;
  NumericType Tpeak_ = NumericType(1500);
  NumericType duration_ = NumericType(1e-3);
  NumericType laserAlpha_ = NumericType(1e4);
};

} // namespace viennaps
```

2. Build and run — expect **PASS** (green) on Tests A and B, and on the
   regression tests C/D/E:

```powershell
cmake --build build --config Release --target testFlashAnneal --parallel
ctest -R testFlashAnneal --test-dir build -C Release --output-on-failure
# expected: "[flash] trap r=1 survival≈0.37", "[flash] pulse … phiMax>0.5",
#           "[flash] melt meanPhi>0.5", "width liquid>solid*2",
#           "meanT latent<noLatent"
```

3. If Test B's `MeltFraction.Max()` is below 0.5 (e.g. the surface spike
   diffused away before melting), raise `I0` (e.g. 1600 → surface 1900 K) or
   reduce `laserAlpha_` (e.g. 1e3, deeper spike); the assertions are
   intentionally on the physics (T > Tm locally, some melting), not on exact
   values.

**Commit**: `feat(fields): FlashAnnealFlow FEM orchestration (heat → melt → dopant)`

## Task 4 — Umbrella exposure, format, full regression

**Files**
- `include/viennaps/viennaps.hpp` (edit)
- `include/viennaps/fields/FlashAnnealFlow.hpp` (format)
- `include/viennaps/fields/models/FlashLaserAnneal.hpp` (format)

**Interfaces**
- Produces: `viennaps::FlashAnnealFlow` + `SolidificationTrapping` reachable
  via the umbrella header.

**Steps**

1. Add next to the other `fields/` includes in `viennaps.hpp`:

```cpp
#include "fields/FlashAnnealFlow.hpp"
```

   (`FlashLaserAnneal.hpp` is already reachable through the models includes.)

2. Format + check + full regression:

```powershell
cmake --build build --target format
cmake --build build --target format-check
ctest -E "Benchmark|Performance" --test-dir build -C Release --output-on-failure
```

3. Commit:

```powershell
git add include/viennaps/fields/FlashAnnealFlow.hpp include/viennaps/fields/models/FlashLaserAnneal.hpp include/viennaps/viennaps.hpp tests/diffusion/testFlashAnneal.cpp tests/diffusion/CMakeLists.txt
git commit -m "feat(fields): flash/laser anneal FEM orchestration (FlashAnnealFlow + SolidificationTrapping)"
```

---

## Self-Review (run before execution handoff)

- **Spec coverage**: GAP_ANALYSIS §4.10 / roadmap MT4 — process-level FEM
  orchestration (the missing piece), Beer's-law pulse IC, solidification
  trapping, latent-heat coupling between segments; the 1D `runPulse` API is
  preserved and untouched. SPER crystallinity is wired into the flow
  (species + model) and exercised by the shared solve.
- **Placeholder scan**: no `TODO`, no `…` in code; all snippets complete.
- **Type consistency**: `prevPhi_` is a `unique_ptr<ParGridFunction>` created
  from `engine.fes()` (raw `ParFiniteElementSpace*` — matches the accessor at
  DiffusionEngine.hpp:227) before copy-assign; `setPreviousPhi` takes
  `const mfem::ParGridFunction*` matching `HeatTransfer`'s existing
  signature; `projectIntegralPreserving` takes `std::vector<NumericType>`;
  `Tpeak_`/`duration_`/`laserAlpha_` are `NumericType` members (the pulse
  IC is seeded via `seedLaserPulse(T0, I0)`; `setPulse` currently records
  parameters for the documented pulse envelope — `runPulse`-style surface
  deposit could be layered on later without API change).
- **Failure-mode honesty**: Test A is a designed experiment (prev φ=1, cur
  φ=0, survival ≈ exp(−r)); Test B asserts existence of melting (Tmax > Tm,
  φmax > 0.5) and dopant dose conservation, not quantitative pulse
  profiles — appropriate for an orchestrator whose per-model physics is
  covered by C/D/E.
- **Latent-heat caveat**: `HeatTransfer`'s latent term uses the previous-φ
  copy; the flow refreshes it per segment, so `∂φ/∂t` is the per-segment
  average — the documented approximation.
- **Risk**: `MeltingPhaseFieldFEM`'s reaction can push φ > 1 (unclamped
  driving force) — Test C asserts bounded overshoot (< 2.0); no fix in
  scope (pre-existing behavior, not introduced here).

## Execution Handoff

Two ways to execute:

1. **Subagent-Driven (recommended)** — each task via a subagent with
   `superpowers:subagent-driven-development`, one task per worktree branch,
   reviewer after each task.
2. **Inline Execution** — execute directly with
   `superpowers:executing-plans`, running the exact build/test commands above.

The user selects which; the plan's TDD steps are self-contained either way.
