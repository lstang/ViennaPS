# Diffusion Engine — Phase 1 Follow-ups

Tracked follow-up items emerging from Phase 1 implementation review.
These are NOT Phase 1 blockers (all Phase 1 tests pass); they are
known scope limits or improvements that Phase 2+ should pick up.

## Important (will bite Phase 2 tasks if not addressed first)

### F1. Dirichlet/Neumann BCs stored but not applied in engine

**Where:** `include/viennaps/fields/DiffusionPhysics.hpp` defines
`addNeumannBC` / `addDirichletBC` / `boundaryConditions` /
`allBoundaryConditions` (lines ~70-98). `DiffusionEngine` never reads
them. `solveImplicitEuler` keeps `essTdof` empty with the comment
"No essential BC in Phase 1 (zero-flux Neumann is natural)."

**Why it's OK for Phase 1:** every Phase 1 test (`TestDiffusionEngineAssembly`,
`TestDoseConservation`) uses zero-flux Neumann (natural) BCs, which is the
default when no essential BCs are applied. Dose conservation validated to
1.28e-16.

**Phase 2 trigger:** the first Dirichlet/segregation/Robin BC test (Phase 2
Task 5 segregation, Phase 4 Task 3 DoseLossBC) will silently behave as
Neumann-zero unless this is fixed.

**Fix:** in `solveImplicitEuler` and `solveCVODE` (the linear-system
construction sites), iterate `physics_->allBoundaryConditions()`:
- For each `BCSpec{type="dirichlet", boundary, value}`: collect boundary
  attributes matching `boundary` from `MeshAttributes`, add their dofs to
  `essTdof`, and apply
  `gf.ProjectBdrCoefficient(...)` or use `EliminateEssentialBC` on the
  system matrix before the solve.
- For each `BCSpec{type="neumann", boundary, value}`: add a
  `mfem::BoundaryLFIntegrator` with the flux coefficient to the RHS
  `LinearForm` (skip if value is zero — natural BC).
- For `"segregation"` / `"robin"`: Phase 2/4 will add the two-sided
  `InterfaceReaction`-style integrator (see Phase 2 Task 5 plan).

**Test coverage to add:** a Dirichlet-BC test (fix C=1e18 on one boundary,
verify the steady-state solution is constant 1e18) and a non-zero Neumann
test (constant flux in, verify dose increases at the expected rate).

### F2. `shouldCreateTimeDerivative` not re-entrant across `solve()` calls

**Where:** `DiffusionPhysics.hpp:53-56` claims the species into
`timeDerivativeClaimed_` and never clears it. The engine calls this in
`assembleAllSpecies()`.

**Why it's OK for Phase 1:** every Phase 1 test creates a fresh physics
object per test, so the gatekeeper is hit once per species per test run.

**Phase 2 trigger:** any test that calls `engine.solve()` twice on the
same physics object (e.g. a multi-step anneal with different temperatures)
will hit the bug: every species gets denied the mass matrix on the second
solve → silent fallback to identity mass → wrong time scale.

**Fix (pick one, document the choice):**
- (A) Clear `timeDerivativeClaimed_` at the top of `assembleAllSpecies()`
  — simplest, matches "gatekeeper is per-solve-call".
- (B) Make the gatekeeper take a "session id" or call scope — more
  general, overkill for now.
- (C) Change semantics so the *first model* on a species claims the
  derivative, evaluated per-solve, with no persistent claim.

Recommended: **(A)**. The MOOSE `PhysicsBase::shouldCreateTimeDerivative`
semantics it mirrors are per-add-kernel, which in our case is per-solve.

### F3. No multi-species integration test

**Where:** both engine tests (`TestDiffusionEngineAssembly`,
`TestDoseConservation`) use exactly 1 species.

**Why it's OK for Phase 1:** `ConstantDiffusion` is a 1-species model;
multi-species is a Phase 3 (CDD) concern.

**Phase 2 trigger:** the species-outer loop (`DiffusionEngine.hpp:269`),
`allSpecies_` map plumbing, packed-block CVODE layout
(`DiffusionEngine.hpp:397-402`), and `modelTargetsSpecies` filter are all
multi-species code paths exercised only by inspection. The first Fermi+Cdd
composition test (Phase 3 Task 10) will be the first real multi-species
exercise.

**Fix:** add a 2-species smoke test in Phase 2 Task 1 — register two
`ConstantDiffusion` models for two species, initialize each, solve, verify
each species conserves its own dose independently. Cheap test, big safety
net.

## Medium (performance / cleanliness)

### F4. Implicit-Euler rebuilds mass matrix every step

**Where:** `solveImplicitEuler` (DiffusionEngine.hpp ~line 538)
re-assembles `Aform` (MassIntegrator + dt*K) every time step, even when
`dt` is constant.

**Why it's OK for Phase 1:** meshes are tiny (4x4, 8x8); rebuild cost is
negligible relative to the linear solve. Also, rebuilding is *correct*
behavior when D is concentration-dependent (Phase 2+ nonlinear) because
K changes each step.

**Phase 2+ fix:** cache `(M + dt K)` and only rebuild when `dt` changes
OR when a model's `assembleStiffness` produces a different K (Picard
update for nonlinear D). Add a `bool dirtySystem_` flag.

### F5. MFEM built without SUNDIALS / MPI / HYPRE

**Where:** `f:/dev/mfem/build` MFEM 4.9.1 was built with neither
`MFEM_USE_SUNDIALS`, `MFEM_USE_MPI`, nor `MFEM_USE_HYPRE`. The CVODE and
HypreBoomerAMG code paths in `DiffusionEngine.hpp` are gated by
`MFEM_USE_SUNDIALS` / `MFEM_USE_MPI` and compile cleanly, but the runtime
falls back to implicit-Euler + DSmoother-preconditioned CG.

**Why it's OK for Phase 1:** the fallback is explicitly in-spec
("Falls back to implicit Euler if SUNDIALS unavailable" — Phase 1 Task 5
brief). Dose conservation validated.

**Phase 2+ fix:** rebuild MFEM with
`MFEM_USE_MPI=ON MFEM_USE_SUNDIALS=ON MFEM_USE_HYPRE=ON` (HYPRE + SUNDIALS
libs are already in vcpkg). After rebuild, the CVODE + HypreBoomerAMG
paths activate automatically (no engine code changes needed) and Phase 2
tests should verify the production solver path works.

## Documentation / scope

### F6. Conforming cut-cell mesh deferred (not just a code comment)

**Where:** `include/viennaps/fields/LevelSetToMesh.hpp` — Phase 1 ships
Cartesian + attribute tagging only.

**What's deferred (per Phase 1 Task 4 scope note):**
- Conforming cut-cell mesh generation at material boundaries (the MOOSE
  `framework/include/meshgenerators/CutMeshByLevelSetGenerator.h` analog).
  MOOSE's generator takes a FunctionParser string for the level set;
  ViennaPS level sets are discrete grids, so direct reuse requires a
  ViennaPS-side `pointLevelSetRelation` adapter (a small bit of new code).
- Transition-layer option (`_generate_transition_layer`) — the standard
  technique to avoid sliver elements at material interfaces. Relevant when
  segregation accuracy at sliver elements becomes the bottleneck.
- `MFEMCutTransitionSubMesh` (`framework/include/mfem/submeshes/`)
  integration for cut-region labeling — relevant to Phase 11 Task 9
  (AMR during moving boundary) and Phase 2 Task 5 (segregation at
  Si/SiO2 interface).

**Trigger to revisit:** Phase 2 Task 5 segregation, or whenever the
flat Cartesian mesh's attribute-tagging produces sliver elements that
degrade segregation accuracy at material interfaces.

### F7. Plan filename inconsistency (cosmetic)

The plan body referred to `tests/diffusion/testDiffusionEngine.cpp` while
the CMakeLists uses `project(testDiffusion)` + `${PROJECT_NAME}.cpp`,
which makes the actual file `testDiffusion.cpp`. Commit `9003343` fixed
the plan; tracked here so future plan edits don't regress it.
