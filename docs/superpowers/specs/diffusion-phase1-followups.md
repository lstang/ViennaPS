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

**Status (2026-07-22): RESOLVED.** MFEM rebuilt with `MFEM_USE_MPI=ON`,
`MFEM_USE_SUNDIALS=ON`, `MFEM_USE_HYPRE=ON` (HYPRE 2.23.0). ViennaPS
configures and builds cleanly against the new MFEM. See F8 and F9 for
follow-ups that emerged from the rebuild.

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

## Post-MFEM-rebuild follow-ups (2026-07-22)

After MFEM was rebuilt with `MFEM_USE_MPI`/`MFEM_USE_SUNDIALS`/`MFEM_USE_HYPRE`,
two new issues surfaced during verification.

### F8. Convert DiffusionEngine to parallel MFEM objects (HypreBoomerAMG)

**Where:** `include/viennaps/fields/DiffusionEngine.hpp` —
`makeMassSolver()` uses `DSmoother` unconditionally.

**Why it's needed:** `mfem::HypreBoomerAMG::SetOperator` (hypre.cpp:5384)
asserts `new Operator must be a HypreParMatrix` — it requires a parallel
distributed matrix as input. Our engine builds serial `SparseMatrix`
operators on a serial `FiniteElementSpace` (FES from `setMesh` is
`FiniteElementSpace`, not `ParFiniteElementSpace`). `MFEM_USE_MPI` being
defined is necessary but not sufficient: the engine itself must construct
`ParFiniteElementSpace` + `ParBilinearForm` + `ParGridFunction` and run
under `MPI_COMM_WORLD` (rank 1 is fine for serial testing).

For Phase 1's 4x4 and 8x8 test meshes (25-81 DOFs) DSmoother-preconditioned
CG converges in <10 iterations, so the parallel path offers no benefit.
Phase 2+'s larger meshes (>1000 DOFs) will benefit substantially.

**Fix:** convert `DiffusionEngine` to use parallel MFEM objects:
- `mesh_` becomes `std::unique_ptr<mfem::ParMesh>` (construct from a serial
  `mfem::Mesh` via `ParMesh(MPI_COMM_WORLD, mesh)`).
- `fes_` becomes `std::unique_ptr<mfem::ParFiniteElementSpace>`.
- `species_` stores `std::unique_ptr<mfem::ParGridFunction>` per species.
- `assembleAllSpecies` builds `ParBilinearForm` per species.
- `makeMassSolver` switches on `MFEM_USE_MPI`: returns `HypreBoomerAMG`-
  preconditioned `CGSolver` for `ParBilinearForm`; keeps `DSmoother` for
  serial fallback.
- The CVODE `DiffusionRHSOperator` works against `HypreParMatrix` (CVODE
  itself is rank-0 only in `SundialsNVector`, but the linear solve via
  `HypreBoomerAMG` is parallel).
- `getIntegral` reduces across MPI ranks (use `ParGridFunction::ComputeL2Norm`
  or similar; the L1 integral requires `MPI_Allreduce`).

**Test:** introduce a 32x32 or 64x64 mesh test where DSmoother-preconditioned
CG takes >100 iterations; verify HypreBoomerAMG converges in <20. Run
single-rank first; multi-rank testing is a separate concern.

### F9. CVODE BDF path is flaky (~20% failure rate)

**Where:** `include/viennaps/fields/DiffusionEngine.hpp` — `solveCVODE()`
via `DiffusionRHSOperator` (implements `Mult`, `SUNImplicitSetup`,
`SUNImplicitSolve`).

**Symptom:** when the engine runs CVODE (env var `VIENNAPS_USE_CVODE=1`
set), the test suite fails ~20% of the time with `Verification failed:
(IsFinite(nom)) is false: nom = -nan` in `CGSolver::Mult`. The crash
happens inside CVODE's Newton iteration; adding `std::cerr` debug prints
to `Mult` masks the issue (suggests timing-sensitive uninitialized read).

**Verified working:** when CVODE does complete (most runs), results are
correct and BETTER than the implicit-Euler fallback:
- Dose conservation: 1.28e-16 (same as implicit-Euler)
- Dirichlet BC: interior reaches 8.35e17 of clamped 1e18 (vs 7.6e17 with
  implicit-Euler — CVODE's smaller adaptive steps give better accuracy)
- Neumann BC: rel_diff 0.4% (vs 4% with implicit-Euler — 10x better)

**Likely cause:** either (a) a SUNDIALS internal SUNContext race, (b)
uninitialized memory in CVODE's predictor, or (c) the operator's
`implicitSolvers_` map being read during a step where setup hasn't yet
run. The SUNImplicitSolve code returns -1 if setup hasn't run, which
should be safe; the issue is more likely inside SUNDIALS itself or in
the operator↔CVODE bridge.

**Workaround (current):** `solve()` selects the integrator at runtime via
the `VIENNAPS_USE_CVODE` env var. Default is implicit-Euler (rock-solid:
8/8 runs pass). CVODE can be opted into per-run for experimentation.

**Fix candidate:** rebuild SUNDIALS from source against the same MSVC
toolchain + flags as MFEM, then run under AddressSanitizer to localize
the uninitialized read. May also be a SUNDIALS 7.x SUNContext thread-
safety issue (we're single-threaded, but SUNContext is sometimes shared
across N_Vectors).
