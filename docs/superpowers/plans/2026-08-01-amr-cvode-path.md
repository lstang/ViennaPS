# Runtime AMR on the CVODE Path Implementation Plan

> For agentic workers: REQUIRED SUB-SKILL — `test-driven-development` for every task
> (write the failing test first, verify it fails, implement, verify it passes, commit);
> `verification-before-completion` before claiming any task done; execute this plan
> via `executing-plans` (inline) or `subagent-driven-development` (subagents), one
> task at a time. Do not skip the format step (Task 4).

## Goal

Close the short-term gap "runtime AMR on the CVODE path + refinebox boxes" from
`docs/GAP_ANALYSIS.md` §4.7: today `DiffusionEngine::solve` silently forces
fixed-step implicit Euler whenever runtime AMR is enabled
(`if (forceImplicitEuler_ || runtimeAmr_)`), because refinement was only
implemented on the implicit-Euler path. This plan makes the CVODE (BDF) path
perform runtime AMR via a **checkpoint–restart loop**: integrate a segment,
refine + prolong the species fields, rebuild the packed state, systems,
RHS operator and a fresh `CVODESolver`, and continue from the current time.

The public API is unchanged (`setRuntimeAmrBox` etc. keep their signatures);
a new read-only getter `usedImplicitEulerPath()` makes the integrator choice
observable (and is the red-green hook for the test).

## Architecture

```
solve() dispatch (today)                     solve() dispatch (after)
  if (forceImplicitEuler_ || runtimeAmr_)      if (forceImplicitEuler_)
      └─> solveImplicitEuler  (+ WARNING)          └─> solveImplicitEuler
  else                                            else
      └─> solveCVODE  (no AMR)                        └─> solveCVODE  (AMR via
                                                        checkpoint–restart)

solveCVODE (after):
  loop { cvode->Step(state, t, dt);
    if (runtimeAmr_ && (++amrStepCounter_ % amrEvery_ == 0)) {
      1) sync  state ──> species GridFunctions      (GFs are STALE mid-loop;
                                                     state is the live copy)
      2) refineBetweenSteps()  ──> GeneralRefinement + EnsureNCMesh +
                                    prolongation, fes_ update
      3) repack GFs ──> state at the new ndof/totalSize
      4) systems = assembleAllSpecies();
         buildIntegrator();   // NEW DiffusionRHSOperator + NEW CVODESolver
                              // (a solver's operator binds at Init; ReInit
                              // cannot swap operators → fresh solver per
                              // segment)
    } }
```

Critical ordering (root-cause note): `refineBetweenSteps()` operates on the
species GridFunctions; mid-integration those GFs are **stale** (the live
solution is the local packed `state` Vector). The state→GF sync MUST happen
before the refine call, or the prolongation would propagate old data.

## Tech Stack

- C++20, MFEM + SUNDIALS build (`MFEM_USE_SUNDIALS` defined in the MFEM at
  `f:/dev/mfem/build`), header `include/viennaps/fields/DiffusionEngine.hpp`
  (edit), tests in `tests/diffusion/testDiffusion.cpp` (edit).
- Test observability: new getter `bool usedImplicitEulerPath() const` +
  member `usedImplicitEuler_` (set per `solve()` dispatch).

## Global Constraints

- C++20; LLVM clang-format (2-space, 80-col).
- The implicit-Euler AMR path must remain byte-for-byte behavior-identical
  (all existing AMR tests keep passing unchanged).
- The existing `DiffusionRHSOperator` (iterative_mode = false inner CG,
  essential-dof zeroing in `Mult`, `SUNImplicitSetup` caching of
  `(M + gamma·K)`) is untouched — the plan only changes how operator/solver
  instances are built and rebuilt.
- Keep the `state`-before-solver declaration order (destruction order:
  `~CVODESolver` may touch the state vector).
- Tests require the AGENTS.md machine-specific CMake setup; commit style
  `feat(fields): ...`.

---

## Task 1 — Failing tests + integrator observability

**Files**
- `include/viennaps/fields/DiffusionEngine.hpp` (edit — getter + member +
  set it in the two dispatch branches)
- `tests/diffusion/testDiffusion.cpp` (edit — append the three AMR-CVODE
  test blocks at the end of `RunTest()`, before its closing brace)

**Interfaces**
- Produces: `bool usedImplicitEulerPath() const` — true iff the last
  `solve()` ran `solveImplicitEuler` (forced by flag, by non-SUNDIALS build,
  or — before this plan — by AMR fallback).
- Consumes: existing `setRuntimeAmrBox(x0, x1, y0, y1, everyNSteps)`,
  `runtimeAmrRefineCount()`, `lastAmrMarkCount()`, `getIntegral(name)`.

**Steps**

1. In `DiffusionEngine.hpp`:
   - Add the public getter next to `runtimeAmrRefineCount()` (public section):

```cpp
  /// True if the last solve() used the fixed-step implicit-Euler path
  /// (forced via setForceImplicitEuler, or a build without SUNDIALS).
  /// Runtime AMR used to force this fallback; after the CVODE-AMR plan it
  /// no longer does.
  bool usedImplicitEulerPath() const { return usedImplicitEuler_; }
```

   - Add the member next to `forceImplicitEuler_` (private state):

```cpp
  bool usedImplicitEuler_ = false;
```

   - In `solve()`, inside `#ifdef MFEM_USE_SUNDIALS`, set the flag in the
     existing dispatch (this makes the test red before the real change):

```cpp
    if (forceImplicitEuler_ || runtimeAmr_) {
      usedImplicitEuler_ = true;
      // (existing warning + solveImplicitEuler call unchanged)
    } else {
      usedImplicitEuler_ = false;
      solveCVODE(tStart, tEnd, dtMax);
    }
```

   - In the `#else` (no-SUNDIALS) branch, set `usedImplicitEuler_ = true;`.

2. Append to `tests/diffusion/testDiffusion.cpp` (end of `RunTest()`):

```cpp
  // --- Runtime AMR on the CVODE path: refinement must run WITHOUT the
  // implicit-Euler fallback, and dose must stay conserved.
  {
    MeshAttributes attrs;
    attrs.setAttributeName(1, "Si");
    auto mesh = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));
    DiffusionEngine<double, 2> engine;
    engine.setMesh(std::move(mesh), attrs);
    auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
    model->setDiffusivity(1e-4, 0.0);
    DiffusionPhysics<double> physics;
    physics.addSpecies("Boron");
    physics.addModel(model);
    physics.setTemperature(1273.0);
    engine.setPhysics(physics);
    engine.initializeSpecies("Boron", 1e18);
    engine.setRuntimeAmrBox(0.25, 0.75, 0.25, 0.75, /*everyNSteps=*/2);
    const double d0 = engine.getIntegral("Boron");
    engine.solve(0.0, 0.2, 0.05);
    const double d1 = engine.getIntegral("Boron");
    const double rel = std::abs(d1 - d0) / std::max(d0, 1.0);
    std::cout << "[amr-cvode] implicitPath=" << engine.usedImplicitEulerPath()
              << " refineCount=" << engine.runtimeAmrRefineCount()
              << " marks=" << engine.lastAmrMarkCount() << " doseRel=" << rel
              << "\n";
    VC_TEST_ASSERT(!engine.usedImplicitEulerPath());
    VC_TEST_ASSERT(engine.runtimeAmrRefineCount() > 0);
    VC_TEST_ASSERT(engine.lastAmrMarkCount() > 0);
    VC_TEST_ASSERT(rel < 0.05);
  }

  // --- AMR CVODE vs AMR implicit-Euler: same physics, comparable dose
  // (both integrators solve the same ODE; CVODE is adaptive, Euler is
  // fixed-step, so allow 5% agreement).
  {
    MeshAttributes attrs;
    attrs.setAttributeName(1, "Si");
    auto runAmr = [&](bool forceEuler) {
      auto mesh = std::make_unique<mfem::Mesh>(
          mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));
      DiffusionEngine<double, 2> engine;
      engine.setMesh(std::move(mesh), attrs);
      auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
      model->setDiffusivity(1e-4, 0.0);
      DiffusionPhysics<double> physics;
      physics.addSpecies("Boron");
      physics.addModel(model);
      physics.setTemperature(1273.0);
      engine.setPhysics(physics);
      engine.initializeSpecies("Boron", 1e18);
      engine.setRuntimeAmrBox(0.25, 0.75, 0.25, 0.75, /*everyNSteps=*/2);
      engine.setForceImplicitEuler(forceEuler);
      engine.solve(0.0, 0.2, 0.05);
      return engine.getIntegral("Boron");
    };
    const double dCvode = runAmr(false);
    const double dEuler = runAmr(true);
    const double rel = std::abs(dCvode - dEuler) / std::max(dEuler, 1.0);
    std::cout << "[amr-cvode] cvode=" << dCvode << " euler=" << dEuler
              << " rel=" << rel << "\n";
    VC_TEST_ASSERT(rel < 0.05);
  }
```

3. Build and run — expect **FAIL** (red) on the first block's first assert:
   AMR still forces the implicit path, so `usedImplicitEulerPath()` returns
   true.

```powershell
cmake --build build --config Release --target testDiffusion --parallel
ctest -R testDiffusion --test-dir build -C Release --output-on-failure
# expected: "[amr-cvode] implicitPath=1 …" and VC_TEST_ASSERT(!usedImplicitEulerPath()) fails
```

**Commit**: none yet.

## Task 2 — Dispatch: AMR no longer forces implicit Euler

**Files**
- `include/viennaps/fields/DiffusionEngine.hpp` (edit — `solve()` dispatch)

**Interfaces**
- Produces: `solve()` routes AMR-enabled runs to `solveCVODE` when
  `forceImplicitEuler_` is false.
- Consumes: unchanged.

**Steps**

1. Replace the dispatch in `solve()`:

```cpp
#ifdef MFEM_USE_SUNDIALS
    if (forceImplicitEuler_) {
      usedImplicitEuler_ = true;
      solveImplicitEuler(tStart, tEnd, dtMax);
    } else {
      usedImplicitEuler_ = false;
      solveCVODE(tStart, tEnd, dtMax);
    }
#else
    usedImplicitEuler_ = true;
    solveImplicitEuler(tStart, tEnd, dtMax);
#endif
```

   Remove the `warnedAmrCvode` block and the "Runtime AMR … ONLY on the
   implicit-Euler path" comment paragraphs; replace with a note that AMR is
   handled inside `solveCVODE` via checkpoint–restart.

2. Build and run — the test **will still fail to compile** until Task 3 makes
   `solveCVODE` AMR-capable? No: `solveCVODE` compiles today; the test's
   `runtimeAmrRefineCount() > 0` assert will fail (CVODE path currently never
   refines). Verify red:

```powershell
cmake --build build --config Release --target testDiffusion --parallel
ctest -R testDiffusion --test-dir build -C Release --output-on-failure
# expected: implicitPath=0 now, but refineCount=0 → assert fails (still red)
```

**Commit**: none yet.

## Task 3 — Checkpoint–restart AMR loop in solveCVODE

**Files**
- `include/viennaps/fields/DiffusionEngine.hpp` (edit — `solveCVODE`)

**Interfaces**
- Consumes: `refineBetweenSteps()` (existing; operates on `mesh_`/`fes_`/
  species GFs, returns true if any element was refined, clears
  implicit-path caches), `assembleAllSpecies()`, `resolveBoundaryMasks(name)`,
  `DiffusionRHSOperator(engine, systems, names, ndof, nSpecies, totalSize,
  bdrMasks)`, `mfem::CVODESolver(CV_BDF)::Init/SetSStolerances/SetMaxStep/
  UseMFEMLinearSolver/Step`.
- Produces: `solveCVODE` executes runtime AMR between segments; `state`,
  `systems`, operator and solver are rebuilt at the refined resolution.

**Steps**

1. In `solveCVODE`, replace the declarations

```cpp
    DiffusionRHSOperator op(*this, systems, names, ndof, nSpecies, totalSize,
                            bdrMasks);

    mfem::CVODESolver cvode(CV_BDF);
    cvode.Init(op);
    cvode.SetSStolerances(/*reltol*/ 1e-6, /*abstol*/ 1e5);
    cvode.SetMaxStep(static_cast<double>(dtMax));
    cvode.UseMFEMLinearSolver();
```

   with (note: `totalSize` must become non-const; keep the
   `state`-declared-first destruction order):

```cpp
    int totalSize = ndof * nSpecies;
    // Declare `state` BEFORE the solver holders so destruction happens in
    // the right order (~CVODESolver may touch the state vector).
    mfem::Vector state(totalSize);
    std::unique_ptr<DiffusionRHSOperator> op;
    std::unique_ptr<mfem::CVODESolver> cvode;

    // Build (or rebuild, after an AMR segment) the RHS operator and a
    // fresh CVODE solver. A solver's operator is bound at Init, so a
    // refined mesh requires a NEW solver; ReInit cannot swap operators.
    auto buildIntegrator = [&]() {
      op = std::make_unique<DiffusionRHSOperator>(
          *this, systems, names, ndof, nSpecies, totalSize, bdrMasks);
      cvode = std::make_unique<mfem::CVODESolver>(CV_BDF);
      cvode->Init(*op);
      cvode->SetSStolerances(/*reltol*/ 1e-6, /*abstol*/ 1e5);
      cvode->SetMaxStep(static_cast<double>(dtMax));
      cvode->UseMFEMLinearSolver();
    };
    buildIntegrator();
```

2. Replace the integration loop

```cpp
    while (t < tFinal) {
      const double targetTime = std::min(t + dt, tFinal);
      dt = targetTime - t;
      cvode.SetMaxStep(tFinal - t);
      cvode.Step(state, t, dt);
      if (t >= tFinal)
        break;
    }
```

   with:

```cpp
    while (t < tFinal) {
      const double targetTime = std::min(t + dt, tFinal);
      dt = targetTime - t;
      cvode->SetMaxStep(tFinal - t);
      cvode->Step(state, t, dt);

      // Runtime AMR checkpoint on the CVODE path (segment restart):
      // 1) sync the live packed state into the species GridFunctions —
      //    mid-integration those GFs are STALE, and refineBetweenSteps()
      //    prolongs whatever is in them, so this sync MUST come first.
      // 2) refine + prolong, 3) repack at the new resolution,
      // 4) rebuild systems + operator + a fresh CVODESolver, restart BDF
      //    from the current time.
      if (runtimeAmr_ && (++amrStepCounter_ % amrEvery_ == 0)) {
        for (int s = 0; s < nSpecies; ++s) {
          auto &gf = *allSpecies_[names[s]];
          mfem::Vector block(state.GetData() + s * ndof, ndof);
          gf = block;
        }
        if (refineBetweenSteps()) {
          ++amrRefineCount_;
          ndof = fes_->GetVSize();
          totalSize = ndof * nSpecies;
          state.SetSize(totalSize);
          for (int s = 0; s < nSpecies; ++s) {
            auto &gf = *allSpecies_[names[s]];
            mfem::Vector block(state.GetData() + s * ndof, ndof);
            block = gf;
          }
          systems = assembleAllSpecies();
          buildIntegrator();
        }
      }

      if (t >= tFinal)
        break;
    }
```

   The final unpack loop already uses the (possibly updated) `ndof`, which is
   correct because `state` was repacked at the final refinement's resolution.

3. Build and run — expect **PASS** (green) on both new blocks AND all
   pre-existing AMR (implicit-path) tests:

```powershell
cmake --build build --config Release --target testDiffusion --parallel
ctest -R testDiffusion --test-dir build -C Release --output-on-failure
# expected: "[amr-cvode] implicitPath=0 refineCount>0 marks>0 doseRel<0.05"
#           and all existing diffusion tests pass
```

4. Check `git diff include/viennaps/fields/DiffusionEngine.hpp` — only
   `solve()`, `solveCVODE()`, the new getter/member, and comments changed.

**Commit**: `feat(fields): runtime AMR on the CVODE path via checkpoint-restart (DiffusionEngine)`

## Task 4 — Format, full regression, docs touch-up

**Files**
- `include/viennaps/fields/DiffusionEngine.hpp` (format)
- `docs/REFINEMENT_REPORT.md` (edit — drop the "AMR only on implicit-Euler"
  caveat, if present)
- `docs/GAP_ANALYSIS.md` (edit — §4.7/§5: mark "AMR on CVODE path" as
  addressed; refinebox boxes remain future work unless
  `setRuntimeAmrBox`-style region boxes are deemed sufficient — they are the
  requested shape)

**Interfaces**
- Consumes: none new.

**Steps**

1. Format + check:

```powershell
cmake --build build --target format
cmake --build build --target format-check
```

2. Full regression:

```powershell
ctest -E "Benchmark|Performance" --test-dir build -C Release --output-on-failure
```

3. Update `docs/REFINEMENT_REPORT.md` and `docs/GAP_ANALYSIS.md` (§4.7, §5
   roadmap ST3): CVODE-path AMR now implemented; note the segment-restart
   behavior (BDF history resets per segment) as a documented characteristic.
4. Commit:

```powershell
git add include/viennaps/fields/DiffusionEngine.hpp tests/diffusion/testDiffusion.cpp docs/REFINEMENT_REPORT.md docs/GAP_ANALYSIS.md
git commit -m "feat(fields): runtime AMR on CVODE path (checkpoint-restart segments)"
```

---

## Self-Review (run before execution handoff)

- **Spec coverage**: §4.7 + roadmap ST3 — AMR now runs on CVODE; the
  "refinebox-style user boxes" gap is covered by the existing
  `setRuntimeAmrBox(x0,x1,y0,y1,everyNSteps)` API (asserted in tests);
  the remaining "refinebox" vocabulary maps to that API.
- **Placeholder scan**: no `TODO`/`…` in code; all snippets complete.
- **Type consistency**: `ndof`/`totalSize` stay `int`; `state.GetData()` +
  `s * ndof` matches the operator's packed layout at every point; the
  operator is rebuilt with the SAME `bdrMasks` (BCs are resolution-
  independent markers — `GetEssentialTrueDofs` re-resolves against the new
  `fes_` inside the operator ctor); `buildIntegrator` captures by reference,
  so `totalSize`/`ndof`/`systems` updates are visible.
- **State-sync ordering**: the state→GF sync precedes `refineBetweenSteps()`
  (advisory-confirmed root-cause requirement); repack precedes
  `assembleAllSpecies()` (systems assemble against the refined `fes_` and
  read the updated GFs only through the operator, which is built after).
- **Risk**: a refine at a segment boundary discards CVODE's BDF history
  (order/error estimates reset). Acceptable: refinement is rare, and the
  implicit path has the same discrete-step granularity.
- **Boundary condition edge**: `cvode->SetMaxStep(tFinal - t)` per segment
  keeps the final segment aligned with `tFinal`, preserving the existing
  chunking semantics.

## Execution Handoff

Two ways to execute:

1. **Subagent-Driven (recommended)** — each task via a subagent with
   `superpowers:subagent-driven-development`, one task per worktree branch,
   reviewer after each task.
2. **Inline Execution** — execute directly with
   `superpowers:executing-plans`, running the exact build/test commands above.

The user selects which; the plan's TDD steps are self-contained either way.
