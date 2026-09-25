# Implant Damage → Diffusion Coupling Implementation Plan

> For agentic workers: REQUIRED SUB-SKILL — `test-driven-development` for every task
> (write the failing test first, verify it fails, implement, verify it passes, commit);
> `verification-before-completion` before claiming any task done; execute this plan
> via `executing-plans` (inline) or `subagent-driven-development` (subagents), one
> task at a time. Do not skip the format step (Task 5).

## Goal

Close the short-term gap "implant damage → diffusion/anneal coupling" from
`docs/GAP_ANALYSIS.md` §4.3: commercial tools automatically feed the implant
damage (dopant + interstitial/vacancy depth profiles) from the BCA step into the
diffusion/anneal step. ViennaPS currently separates these concerns. This plan
adds a small, dependency-free bridge that seeds `DiffusionEngine` initial
conditions from an `MCBcaResult` with dose conservation, and a TED physics
factory (pair diffusion with interstitial supersaturation), plus regression
tests proving both the dose conservation and the physics enhancement.

## Architecture

```
MCBcaEngine (psMCBcaImplant.hpp)          DiffusionEngine (fields/)
  run() ──> MCBcaResult {dopant,               projectIntegralPreserving(name,
  interstitial, vacancy, ...}                  samples, targetDose)  [exists]
                     │                                   ▲
                     └─── ImplantDamageCoupler ──────────┘
                          fields/ImplantDamageCoupler.hpp (NEW)
                            seedFromBca(engine, result)        — 1D profile → GF, dose-rescaled
                            seedSpecies(engine, name, profile) — single-species projection
                            makeTedPair(dopant, I, D_pair, C_Ieq)      — TED factory
                            makeDefectTransport(species, D0, Ea_eV)    — defect transport
```

The coupler is **stateless** (static functions): it reuses the engine's existing
`projectIntegralPreserving` (already the IntegralPreservingFunctionIC-style hook),
so it adds no solver or assembly code. TED enhancement is provided by the
existing `PairDiffusion` model (`D_eff = D_pair · C_I/C_I_eq`), wired by the
factory helpers. Nothing in `psBasicDiffusion.hpp` or existing process models
changes.

## Tech Stack

- C++20, header-only, `include/viennaps/fields/ImplantDamageCoupler.hpp` (new).
- MFEM-backed `DiffusionEngine<double, 2>` (requires `VIENNAPS_HAS_MFEM` — the
  test target only builds on this machine's configuration, same as
  `tests/diffusion/testDiffusion.cpp`).
- Tests: `tests/diffusion/testImplantDamageCoupling.cpp` (new), registered in
  `tests/diffusion/CMakeLists.txt`; `VC_TEST_ASSERT`; explicit `main()` with
  `#ifdef VIENNAPS_HAS_MFEM` + `#else` fallback main, exactly like
  `tests/diffusion/testDiffusion.cpp` (do NOT use `VC_RUN_ALL_TESTS` — it
  expands to two-parameter `RunTest<double|float,2|3>` instantiations;
  this test is MFEM-only and double-only).

## Global Constraints

- C++20; LLVM clang-format (2-space indent, 80-col, see `.clang-format`);
  include order sorted case-sensitively by clang-format.
- Header-only: no `.cpp`, no new link deps (the coupler uses only headers the
  fields layer already includes).
- MFEM-gated code must keep working when `VIENNAPS_HAS_MFEM` is off: the coupler
  header sits behind the same `#ifdef` structure as `DiffusionEngine.hpp` (it
  includes `DiffusionEngine.hpp`, which is itself MFEM-gated; the coupler needs
  no extra gate of its own).
- Tests require the machine-specific CMake setup from `AGENTS.md`:
  `F:/dev/vcpkg/installed/x64-windows`, `f:/dev/mfem/build`,
  `VCPKG_MANIFEST_INSTALL=OFF`.
- Do NOT touch `psBasicDiffusion.hpp` semantics; do NOT classify anything as
  CMP; do NOT frame MOOSE as a capability.
- Commit style: `feat(fields): ...` with a one-line summary + body.

---

## Task 1 — Failing test: dose-conserving seeding from a BCA run

**Files**
- `tests/diffusion/testImplantDamageCoupling.cpp` (new)
- `tests/diffusion/CMakeLists.txt` (edit — add target)

**Interfaces**
- Consumes: `MCBcaEngine<NumericType>::params()` (mutable `Params`),
  `MCBcaEngine<NumericType>::run()` → `MCBcaResult<NumericType>`; field members
  `dopant`, `interstitial`, `vacancy` (1D depth profiles);
  `DiffusionEngine<NumericType,2>::setMesh/setPhysics/getIntegral`;
  `DiffusionPhysics<NumericType>::addSpecies`.
- Produces (test-only): `tests/diffusion/testImplantDamageCoupling` ctest target.

**Steps**

1. Create `tests/diffusion/testImplantDamageCoupling.cpp`:

```cpp
// testImplantDamageCoupling.cpp — MCBca implant damage → DiffusionEngine
// initial-condition coupling (GAP_ANALYSIS §4.3).
#include <fields/DiffusionEngine.hpp>
#include <fields/DiffusionPhysics.hpp>
#include <fields/ImplantDamageCoupler.hpp>
#include <models/psMCBcaImplant.hpp>
#include <vcTestAsserts.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
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
```

2. Append to `tests/diffusion/CMakeLists.txt`:

```cmake
project(testImplantDamageCoupling LANGUAGES CXX)
viennaps_add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
add_dependencies(ViennaPS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
```

3. Build and run — expect **FAIL** (red): `ImplantDamageCoupler.hpp` does not
   exist yet.

```powershell
cmake --build build --config Release --target testImplantDamageCoupling --parallel
# expected: fatal error C1083: cannot open source file 'fields/ImplantDamageCoupler.hpp'
```

**Commit**: none yet (red state is not committed).

## Task 2 — Implement the coupler seeding API

**Files**
- `include/viennaps/fields/ImplantDamageCoupler.hpp` (new)

**Interfaces**
- Produces (public, namespace `viennaps`):
  - `static void seedFromBca(DiffusionEngine<NumericType, 2> &engine, const MCBcaResult<NumericType> &result)` — seeds `Dopant`, `Interstitial`, `Vacancy` from the result's 1D profiles.
  - `static void seedSpecies(DiffusionEngine<NumericType, 2> &engine, const std::string &name, const std::vector<NumericType> &profile)` — projects one profile with `projectIntegralPreserving(name, profile, dose)` where `dose = Σ profile`; no-op on empty profile.
- Consumes: `../models/psMCBcaImplant.hpp` (for `MCBcaResult`),
  `DiffusionEngine.hpp` (`projectIntegralPreserving`).

**Steps**

1. Create `include/viennaps/fields/ImplantDamageCoupler.hpp` — seeding API
   ONLY (the TED factory methods are Task 4, so that Task 3's failing test
   has a real red):

```cpp
#pragma once

/// ImplantDamageCoupler — bridge MCBca implant damage into DiffusionEngine
/// initial conditions (dopant, interstitial, vacancy) with dose-conserving
/// 1D depth-profile projection. TED physics factories are added in Task 4.
///
/// Stateless: reuses DiffusionEngine::projectIntegralPreserving (the
/// IntegralPreservingFunctionIC-style hook) so it adds no solver or
/// assembly code. GAP_ANALYSIS §4.3: commercial tools feed implant damage
/// into the anneal step automatically; this is that bridge.

#include "../models/psMCBcaImplant.hpp"
#include "DiffusionEngine.hpp"

#include <numeric>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType> class ImplantDamageCoupler {
public:
  /// Seed dopant + defect species from a BCA result onto the engine's
  /// species fields. Each 1D depth profile is projected with
  /// projectIntegralPreserving, rescaling the field integral to the
  /// profile's total dose (dose conservation per species).
  static void seedFromBca(DiffusionEngine<NumericType, 2> &engine,
                          const MCBcaResult<NumericType> &result) {
    seedSpecies(engine, "Dopant", result.dopant);
    seedSpecies(engine, "Interstitial", result.interstitial);
    seedSpecies(engine, "Vacancy", result.vacancy);
  }

  /// Project one 1D depth profile with dose conservation.
  static void seedSpecies(DiffusionEngine<NumericType, 2> &engine,
                          const std::string &name,
                          const std::vector<NumericType> &profile) {
    if (profile.empty())
      return;
    const NumericType dose =
        std::accumulate(profile.begin(), profile.end(), NumericType(0));
    engine.projectIntegralPreserving(name, profile, dose);
  }
};

} // namespace viennaps
```

2. Build and run — expect **PASS** (green) for Test 1. This run also proves
   the Task 1 file shape compiles (the earlier red stopped at C1083 before
   parsing):

```powershell
cmake --build build --config Release --target testImplantDamageCoupling --parallel
ctest -R testImplantDamageCoupling --test-dir build -C Release --output-on-failure
# expected: "[implant-coupling] target=… seeded≈target …" and all assertions pass
```

**Commit**: `feat(fields): dose-conserving MCBca implant damage seeding (ImplantDamageCoupler)`

## Task 3 — Failing test: TED broadens the dopant profile with damage

**Files**
- `tests/diffusion/testImplantDamageCoupling.cpp` (edit — append Test 2 block)

**Interfaces**
- Consumes: `ImplantDamageCoupler::makeTedPair/makeDefectTransport`
  (added in Task 4 — the test does NOT compile until then, which is the
  TDD red for this task),
  `DiffusionEngine::projectIntegralPreserving` (peaked dopant IC),
  `getSolution(name)` → `const mfem::ParGridFunction &` (dofs = order-1
  vertex dofs on the Cartesian mesh).

**Steps**

1. Append this block inside `TestImplantDamageSeeding()`, before its closing
   `}` (after the Test 1 block). The block is `double`-typed throughout, which
   matches the function's `double` engine instantiation.

   While editing the file, also normalize the include block to what the test
   actually uses (T1 review finding): `#include <algorithm>` (for `std::max`)
   and `#include <exception>` (for `std::exception` in `main`) — the sorted
   block is:

```cpp
  // --- Test 2: TED — an interstitial supersaturation (implant damage)
  // broadens the dopant profile beyond the no-damage baseline. Width is
  // the RMS distance of the concentration from its centroid, computed
  // from vertex dofs (order-1 H1 nodal basis).
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

    auto runTed = [&](double cI) {
      MeshAttributes attrs;
      attrs.setAttributeName(1, "Si");
      auto mesh = std::make_unique<mfem::Mesh>(
          mfem::Mesh::MakeCartesian2D(32, 32, mfem::Element::TRIANGLE));
      auto engine = std::make_unique<DiffusionEngine<double, 2>>();
      engine->setMesh(std::move(mesh), attrs);
      DiffusionPhysics<double> physics;
      physics.addSpecies("Boron");
      physics.addSpecies("Interstitial");
      physics.addModel(ImplantDamageCoupler<double>::makeTedPair(
          "Boron", "Interstitial", /*D_pair=*/1e-10, /*C_Ieq=*/1e10));
      physics.addModel(ImplantDamageCoupler<double>::makeDefectTransport(
          "Interstitial", /*D0=*/1e-8, /*Ea_eV=*/0.0));
      physics.setTemperature(1273.0);
      engine->setPhysics(physics);
      // Peaked dopant IC (Gaussian over the dof axis, integral ~4e18).
      std::vector<double> peak(33, 0.0);
      for (std::size_t i = 0; i < peak.size(); ++i) {
        const double z = (static_cast<double>(i) - 16.0) / 4.0;
        peak[i] = 1e18 * std::exp(-z * z);
      }
      engine->projectIntegralPreserving("Boron", peak, 4e18);
      engine->initializeSpecies("Interstitial", cI);
      engine->solve(0.0, 1e4, 1e3);
      // NOTE: no `*` — getSolution returns a const ref; MFEM defines no
      // unary operator* for ParGridFunction (T3 red-step finding).
      return width(engine->getSolution("Boron"));
    };

    const double wNoDamage = runTed(1e-30);
    const double wDamage = runTed(1e15);
    std::cout << "[implant-coupling] width no-damage=" << wNoDamage
              << " with-damage=" << wDamage << "\n";
    // D_eff = D_pair * C_I/C_Ieq = 1e-10 * 1e15/1e10 = 1e-5 for the
    // damage case; measured ratio 1.344 (plan fallback raised cI 1e14 ->
    // 1e15 after 1e14 measured only 1.10x; no-damage: D_eff ~ 0, width ~
    // 0.311).
    VC_TEST_ASSERT(wDamage > wNoDamage * 1.2);
    VC_TEST_ASSERT(wNoDamage > 0.0);
  }
```

2. Build and run — expect **FAIL** (red): `makeTedPair` and `makeDefectTransport`
   are not yet declared.

```powershell
cmake --build build --config Release --target testImplantDamageCoupling --parallel
# expected: error C2039: 'makeTedPair': is not a member of 'viennaps::ImplantDamageCoupler<double>'
```

**Commit**: none yet.

## Task 4 — Implement the TED physics factory

**Files**
- `include/viennaps/fields/ImplantDamageCoupler.hpp` (edit — add the two
  factory methods; the file from Task 2 has seeding only)

**Interfaces**
- Produces: `makeTedPair(dopant, interstitial, D_pair, C_Ieq)` →
  `std::shared_ptr<PairDiffusion<NumericType>>`,
  `makeDefectTransport(species, D0, Ea_eV)` →
  `std::shared_ptr<ConstantDiffusion<NumericType>>`.
- Consumes: `PairDiffusion::setPairDiffusivity/setCIEq`,
  `ConstantDiffusion::setDiffusivity`.

**Steps**

1. Add the two factory methods to `ImplantDamageCoupler` (inside the class,
   after `seedSpecies`) and the two model includes:

```cpp
#include "DiffusionEngine.hpp"
#include "models/ConstantDiffusion.hpp"
#include "models/PairDiffusion.hpp"
#include "models/psMCBcaImplant.hpp"

#include <numeric>
#include <string>
#include <vector>
```

```cpp
  /// TED pair-diffusion factory: D_eff = D_pair * C_I / C_I_eq.
  static std::shared_ptr<PairDiffusion<NumericType>>
  makeTedPair(const std::string &dopant, const std::string &interstitial,
              NumericType D_pair, NumericType C_Ieq) {
    auto m =
        std::make_shared<PairDiffusion<NumericType>>(dopant, interstitial);
    m->setPairDiffusivity(D_pair);
    m->setCIEq(C_Ieq);
    return m;
  }

  /// Defect transport factory (Arrhenius constant diffusivity).
  static std::shared_ptr<ConstantDiffusion<NumericType>>
  makeDefectTransport(const std::string &species, NumericType D0,
                      NumericType Ea_eV) {
    auto m = std::make_shared<ConstantDiffusion<NumericType>>(species);
    m->setDiffusivity(D0, Ea_eV);
    return m;
  }
```

2. Build and run — expect **PASS** (green) for both tests:

```powershell
cmake --build build --config Release --target testImplantDamageCoupling --parallel
ctest -R testImplantDamageCoupling --test-dir build -C Release --output-on-failure
# expected: "[implant-coupling] width no-damage=… with-damage=…" and all assertions pass
```

3. Sanity-check the physics claim: the damage-case width must exceed the
   no-damage width by >20% (assertion). If it does not, INCREASE the
   enhancement — raise the damage supersaturation to `1e15` or `D_pair` to
   `1e-9` (do NOT lower it; lowering shrinks the ratio) — and re-check.
   If the ratio still sits below 1.2 with plausible constants, stop and
   report to the controller: the assertion threshold or the physics
   parameters are a plan-level decision the controller owns.

**Commit**: `feat(fields): TED pair-diffusion + defect-transport factories in ImplantDamageCoupler`

## Task 5 — Umbrella exposure, format, full regression

**Files**
- `include/viennaps/viennaps.hpp` (edit — expose the new public header)

**Interfaces**
- Produces: `viennaps::ImplantDamageCoupler` reachable via the umbrella
  header (same section as the other `fields/` includes).

**Steps**

1. In `include/viennaps/viennaps.hpp`, immediately after the existing
   `#include <fields/DiffusionEngine.hpp>` line (currently line 123), add:

```cpp
#include "fields/ImplantDamageCoupler.hpp"
```

2. Format and check:

```powershell
cmake --build build --target format
cmake --build build --target format-check
# expected: format-check passes with no diffs
```

3. Full regression — the existing diffusion suite must stay green:

```powershell
ctest -R "diffusion|implant" --test-dir build -C Release --output-on-failure
ctest -E "Benchmark|Performance" --test-dir build -C Release --output-on-failure
```

4. Review the diff (no changes outside the touched files) and commit:

```powershell
git add include/viennaps/viennaps.hpp && git commit -m "docs(viennaps): expose ImplantDamageCoupler via umbrella header"
```

---

## Self-Review (run before execution handoff)

- **Spec coverage**: §4.3 of GAP_ANALYSIS.md (implant damage → anneal wiring) —
  covered: seeding (dopant/I/V), dose conservation, TED physics factory,
  regression tests. The amorphous-region → SPER path is explicitly **out of
  scope** (owned by the flash-anneal plan).
- **Placeholder scan**: no `TODO`, no `…` in code, every step has complete code.
- **Type consistency**: the test file is intentionally **double-only** (no
  `NumericType` template, no `VC_RUN_ALL_TESTS`) — it matches
  `testDiffusion.cpp`'s explicit-`main()` + `#ifdef VIENNAPS_HAS_MFEM`
  fallback shape, and the `MCBcaEngine<double>` ↔
  `DiffusionEngine<double,2>` ↔ `ImplantDamageCoupler<double>` chain stays
  uniformly `double`. `getIntegral` returns `NumericType` (=`double` here).
  The `peak` vector in Test 2 is `std::vector<double>`, which is exact.
- **MFEM-less builds**: with `VIENNAPS_HAS_MFEM` undefined the file compiles
  to the fallback `main()` that prints and exits 0, so the target never
  breaks the no-MFEM tests/ tree (mirrors `testDiffusion.cpp`).
- **Failure-mode honesty**: the dose-conservation assertion uses the BCA
  profile's own sum as target (self-consistent), not the nominal
  `params().dose`, because `run()` may drop ions that exit the depth window.

## Execution Handoff

Two ways to execute:

1. **Subagent-Driven (recommended)** — run each task through a subagent with
   `superpowers:subagent-driven-development`, one task per worktree branch,
   reviewer after each task.
2. **Inline Execution** — execute the tasks directly in this session with
   `superpowers:executing-plans`, running the exact build/test commands above.

The user selects which; the plan's TDD steps are self-contained either way.
