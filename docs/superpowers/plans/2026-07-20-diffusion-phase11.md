# Diffusion Engine Phase 11: Adaptive Mesh Refinement

**Goal:** Add adaptive mesh refinement (AMR) during diffusion: refine/coarsen based on solution gradients, dose error, and user-specified criteria. Uses MFEM's built-in AMR.

**Depends on:** Phase 1 (DiffusionEngine, LevelSetToMesh)

## File Structure

| File | Responsibility |
|------|---------------|
| `AdaptiveMeshRefiner.hpp` | AMR controller: criteria evaluation, refine/coarsen |
| `RefinementBox.hpp` | User-specified static refinement region |
| `MeshQualityEstimator.hpp` | Skewness, Jacobian, aspect ratio checks |

## Tasks

### Task 1: RefinementBox - Static Refinement
- `RefinementBox<NumericType, D>` - user-specified rectangular region with max element size. Before diffusion, refine mesh inside box to specified resolution. Uses MFEM `Mesh::GeneralRefinement` with attribute-based region selection
- Test: 8x8 mesh, refine box in corner -> verify corner elements smaller than rest
- Commit: `"feat: add RefinementBox for static mesh refinement"`

### Task 2: MeshQualityEstimator
- `MeshQualityEstimator` - computes per-element quality metrics: skewness, Jacobian determinant, aspect ratio. Identifies elements needing refinement or remeshing
- Test: regular mesh -> all quality ~ 1.0. Distorted mesh -> low quality elements flagged
- Commit: `"feat: add MeshQualityEstimator for element quality assessment"`

### Task 3: Relative Difference Refinement Criterion
- `RelativeDifferenceCriterion<NumericType>` - refine elements where |C_new - C_old| / max(|C|, eps) > threshold. Coarsen where below. Tracks solution change between time steps.

**MOOSE analog (verified):** `framework/include/markers/ValueThresholdMarker.h` driven by an aux-variable storing |C_new - C_old| — the standard MOOSE idiom for time-step-difference refinement. Implementing this criterion = maintain an aux `dC` GridFunction updated after each solve, then threshold.

- Test: sharp profile → refine near front. Uniform profile → no refinement.
- Commit: `"feat: add relative difference AMR criterion (MOOSE ValueThresholdMarker pattern)"`

### Task 4: Gradient Refinement Criterion
- `GradientCriterion<NumericType>` - refine elements where |grad(C)| > threshold. Captures steep concentration gradients (implant profiles, junctions).

**MOOSE references (verified — prefer the indicator-driven path):**
- `framework/include/indicators/GradientJumpIndicator.h` — `InternalSideIndicator` subclass that integrates gradient jumps across internal sides. **This is the closest built-in to `GradientCriterion` and is theoretically grounded** (jump indicators have a posteriori error estimates for diffusion). Prefer implementing the criterion as an *indicator* feeding Task 7's `MFEMRefinementMarker`-style refiner, not as a free-standing criterion.
- `framework/include/indicators/LaplacianJumpIndicator.h` — jump in Laplacian (higher-order analog for smooth solutions).

**Action:** expose `GradientCriterion` either as a free-standing criterion OR (preferred) as an indicator that produces a per-element error estimate consumed by Task 7's `AdaptiveMeshRefiner`. The indicator path integrates naturally with `MFEML2ZienkiewiczZhuIndicator` (Task 7) for combined error estimation.

- Test: error function profile → refine near steepest gradient. **Add:** compare element marking against analytic `|grad C|` profile within 5%.
- Commit: `"feat: add gradient-based AMR criterion as indicator (MOOSE GradientJumpIndicator pattern)"`

### Task 5: Local Dose Error Criterion
- `LocalDoseErrorCriterion<NumericType>` - refine elements where integral of |C| change exceeds threshold. Ensures dose accuracy in critical regions.

**MOOSE analog (verified):** `framework/include/markers/ErrorFractionMarker.h` driven by an `ElementIntegralIndicator` (verified at `framework/include/indicators/ElementIntegralIndicator.h`). The MOOSE pattern: indicator computes per-element integral, marker picks the top fraction. Implement the same way — `LocalDoseErrorCriterion` produces a per-element dose-error indicator; Task 7's `AdaptiveMeshRefiner` thresholds it via `ErrorFractionMarker`-style logic.

- Test: high-dose region → refined. **Add:** total dose error (sum over elements) decreases monotonically as threshold tightens.
- Commit: `"feat: add local dose error AMR criterion as indicator (MOOSE ErrorFractionMarker + ElementIntegralIndicator pattern)"`

### Task 6: Logarithmic and Asinh Criteria

- `LogarithmicCriterion` - |log(C_new/C_old)| > threshold. Better for concentrations spanning many orders of magnitude.
- `AsinhCriterion` - inverse hyperbolic sine difference. Handles sign changes and wide dynamic range.

**No MOOSE analog — these are genuinely new.** Dopant concentrations span 10+ orders of magnitude (1e10 substrate to 1e20 implant), so linear `RelativeDifferenceCriterion` (Task 3) either over-refines the substrate (where tiny absolute changes look large relatively) or under-refines the implant tail. `LogarithmicCriterion` and `AsinhCriterion` are the correct response. Well-justified additions; cite this rationale in the commit.

- Test: C from 1e10 to 1e20 → log criterion refines where orders of magnitude change. **Add:** on the same profile, `RelativeDifferenceCriterion` either over-refines near 1e10 OR under-refines near 1e15 — log criterion avoids both.
- Commit: `"feat: add logarithmic and asinh AMR criteria (no MOOSE analog; justified for wide-dynamic-range concentration fields)"`

### Task 7: AdaptiveMeshRefiner Controller

- `AdaptiveMeshRefiner<NumericType, D>` - orchestrates AMR. Collects active criteria + an optional error estimator, evaluates per-element, marks for refine/coarsen, calls MFEM refinement, transfers solution. Configurable: max h-level, max p-level, min element size, coarsen threshold, **refinement mode (h-only / p-only / hp)**, **combination policy (max-error-wins default, OR, AND)**.

**MOOSE MFEM reference (PRIMARY — verified):** `framework/include/mfem/markers/MFEMRefinementMarker.h` already wraps `mfem::ThresholdRefiner` with exactly the knobs we need:
- `_error_threshold` (real_t) — element error above this gets refined
- `_max_h_level`, `_max_p_level` (unsigned) — independent caps for h and p
- `_rebalance` (bool) — repartition the mesh after h-refinement (parallel load balance)
- `hRefine()` / `pRefine()` — separate methods; **p-refinement is MFEM's hallmark and the plan must not omit it.** For diffusion with smooth solutions away from junctions, p-refinement is dramatically more efficient than h-refinement (exponential vs algebraic convergence).

**MOOSE MFEM reference for error estimation (verified):** `framework/include/mfem/indicators/MFEML2ZienkiewiczZhuIndicator.h` — wraps the Zienkiewicz-Zhu flux-projection estimator (MFEM example 6p). This is a true residual-based error estimator, not a heuristic on solution values. The current plan's heuristic criteria (`GradientCriterion`, `RelativeDifferenceCriterion`) should be **additions to** this, not replacements for it.

**MOOSE framework reference for combination policy (verified):** `framework/include/markers/ComboMarker.h` combines multiple markers; the convention is max-error-wins (each marker returns REFINE/DONT_REFINE/COARSEN, the most aggressive action wins). The plan must pick a policy explicitly — currently it says "combined marking" without defining the combination logic.

**MOOSE framework reference for interface protection (verified):** `framework/include/markers/BoundaryPreservedMarker.h` prevents coarsening from eating a tracked boundary. Relevant for the Si/SiO2 interface during AMR cycles — without it, coarsening can remove the very interface elements that segregation (Phase 2 Task 5) needs.

**Two projection cases — DO NOT CONFLATE (verified from `framework/src/transfers/MultiAppProjectionTransfer.C`):**
1. **Same mesh after AMR** → use MFEM's native prolongation/restriction operators: `FiniteElementSpace::Update()` + `GridFunction::SetUpdateOperator()`. O(N), no solve.
2. **Cross-mesh** (KMC bins → FEM, master ↔ sub-app, or remesh after oxidation in Task 9) → L2 projection solve: build mass matrix `M_ij = ∫φ_iφ_j`, RHS `b_i = ∫source_value*φ_i`, solve `Mx=b`. This is exactly `MultiAppProjectionTransfer::assembleL2()` (verified lines 125-186).

The previous plan said "transfer solution via L2 projection" for both cases — wrong for case 1 (slow, unnecessary solve).

- [ ] **Step 1: Write failing test** - `TestAdaptiveMeshRefinerHP()`:
  - (a) Smooth Gaussian-bump concentration on 4x4 mesh, run refiner with `mode=hp`, assert final error (L2 ZZ estimator) decays faster with hp than h-only at same N_DoF.
  - (b) Multiple criteria active (gradient + ZZ estimator) with `policy=max_wins`; assert no element is refined twice and the union of refined elements equals union of each criterion's set.
  - (c) After refinement, `FiniteElementSpace::Update()` is called and `GridFunction` prolongation preserves total dose within 0.1% (NOT via L2 solve — via native update operator).
  - (d) `BoundaryPreservedMarker` analog: coarsening never removes an element adjacent to the Si/SiO2 interface.

- [ ] **Step 2: Run to verify failure** -> FAIL

- [ ] **Step 3: Implement AdaptiveMeshRefiner** wrapping `mfem::ThresholdRefiner` (mirror `MFEMRefinementMarker.h`):
```cpp
enum class RefinementMode { HOnly, POnly, HP };
enum class CombinationPolicy { MaxWins, Or, And };  // mirrors MOOSE ComboMarker

template <class NumericType, int D>
class AdaptiveMeshRefiner {
public:
  // Mirror MOOSE MFEMRefinementMarker knobs (verified header):
  void setErrorThreshold(NumericType e) { errorThreshold_ = e; }
  void setMaxHLevel(unsigned n) { maxH_ = n; }
  void setMaxPLevel(unsigned n) { maxP_ = n; }      // <-- new; absent in original plan
  void setRebalance(bool b) { rebalance_ = b; }     // parallel load balance
  void setMode(RefinementMode m) { mode_ = m; }     // <-- new
  void setCombinationPolicy(CombinationPolicy p) { policy_ = p; }  // <-- new

  // Error estimator (MFEML2ZienkiewiczZhuIndicator analog) — optional but
  // strongly recommended; this is the only path with theoretical guarantees.
  // Heuristic criteria (Gradient-, RelativeDifference-, etc.) are added via
  // addCriterion() and combined per policy_.
  void setErrorEstimator(std::shared_ptr<mfem::ErrorEstimator> ee) {
    estimator_ = ee;
  }

  void addCriterion(std::shared_ptr<AmrCriterion<NumericType, D>> c) {
    criteria_.push_back(c);
  }

  // Interface-protection marker (BoundaryPreservedMarker analog). Elements
  // on this attribute list are never coarsened.
  void preserveAttribute(int attr) { preservedAttrs_.insert(attr); }

  /// Apply one AMR cycle. Returns true if mesh changed.
  bool refine(mfem::ParMesh& mesh, mfem::ParFiniteElementSpace& fes,
              mfem::ParGridFunction& sol) {
    // 1. Compute per-element errors: estimator_->ComputeEstimates(error_, fes, sol)
    //    OR if no estimator, evaluate heuristic criteria_ and combine per policy_.
    // 2. Build mfem::ThresholdRefiner; set maximum hypotheses per mode_:
    //      HOnly  -> thresholdRefiner.SetTotalErrorFraction / max_h_level
    //      POnly  -> pRefine() only
    //      HP     -> hRefine() then pRefine() (mirror MFEMRefinementMarker.h:33-36)
    // 3. Apply BoundaryPreservedMarker analog: un-mark any element in
    //    preservedAttrs_ from the coarsen set.
    // 4. mesh.Refine(refiner) -> fes.Update() -> sol.SetUpdateOperator()
    //    (CASE 1: native prolongation, NOT L2 solve)
    // 5. If rebalance_ and parallel: mesh.Rebalance()
  }

private:
  NumericType errorThreshold_ = 1e-3;
  unsigned maxH_ = 4, maxP_ = 2;
  bool rebalance_ = false;
  RefinementMode mode_ = RefinementMode::HP;     // default to HP — MFEM's strength
  CombinationPolicy policy_ = CombinationPolicy::MaxWins;
  std::shared_ptr<mfem::ErrorEstimator> estimator_;
  std::vector<std::shared_ptr<AmrCriterion<NumericType, D>>> criteria_;
  std::set<int> preservedAttrs_;                  // BoundaryPreservedMarker analog
};
```

- [ ] **Step 4: Run to verify pass** -> PASS

- [ ] **Step 5: Commit** - `"feat: add AdaptiveMeshRefiner with hp-refinement, ZZ error estimator, and ComboMarker-style policy (MOOSE MFEMRefinementMarker pattern)"`

### Task 8: AMR During Diffusion

- Integrate AMR into `DiffusionEngine::solve()`. Evaluate criteria, refine/coarsen, transfer solution, continue.

**MOOSE framework reference (verified — distinguishes two AMR timings):** `framework/include/base/Adaptivity.h` (lines 50-159) separates two distinct AMR timings that the original plan conflated into "every N time steps":
1. **`initial_steps`** — pre-solve AMR. Run several refinement cycles at t=0 (before the first time step) to resolve the initial condition. Essential for sharp implant profiles where the IC itself triggers refinement.
2. **`cycles_per_step`** — per-step AMR cycles. Number of refine/coarsen cycles per time step. Usually 1, but for stiff nonlinear diffusion (clustering in Phase 3) the within-cycle recomputation matters.
3. **`recompute_markers_during_cycles`** — whether markers are re-evaluated within a multi-cycle step (true) or only at the start (false). True is more accurate but slower.

**Action:** expose three knobs on `DiffusionEngine::solve()`:
```cpp
engine.setAMRInitialSteps(unsigned n);      // Adaptivity::init initial_steps
engine.setAMRCyclesPerStep(unsigned n);     // Adaptivity::init cycles_per_step (default 1)
engine.setAMRInterval(unsigned n);          // run AMR every N steps (default 1); perf knob
engine.setAMRRecomputeMarkers(bool b);      // recompute_markers_during_cycles
```
Default: `initial_steps=2`, `cycles_per_step=1`, `interval=1`, `recompute=true`. Document each.

- Test: implant profile (sharp) → AMR refines near junction, coarse away. Solution accuracy maintained. **Add:** with `initial_steps=4` vs `initial_steps=0`, the first-step dose error is <0.1% vs >5% — demonstrating why `initial_steps` matters.
- Commit: `"feat: integrate AMR into DiffusionEngine with initial_steps / cycles_per_step / interval knobs (MOOSE Adaptivity.h pattern)"`
- Commit: `"feat: integrate AMR into DiffusionEngine time stepping"`

### Task 9: AMR During Moving Boundary

- For oxidation: after mesh deformation, check quality. If degraded, trigger remesh + AMR. Refine near moving interface.

**MOOSE pattern (verified three-class composition):** This is the canonical MOOSE recipe for "AMR while a mesh is also being deformed/transferred," composed of three pieces:
1. **`modules/level_set/include/base/LevelSetProblem.h`** — overrides `adaptMesh()` to fire MultiAppTransfers at `EXEC_ADAPT_MESH`. The ViennaPS analog: have `DiffusionEngine::solve()` hook AMR into the same loop that applies oxide-driven mesh motion, with an explicit "after mesh moved" callback.
2. **`modules/level_set/include/transfers/LevelSetMeshRefinementTransfer.h`** — pushes the refinement marker from the master problem to the sub-app; `initialSetup()` configures sub-app adaptivity once, `execute()` toggles adaptivity only during `EXEC_ADAPT_MESH`. Cite this for any sub-app-based coupling.
3. **`framework/src/transfers/MultiAppProjectionTransfer.C` (lines 125-186)** — the L2-projection algorithm (CASE 2 from Task 7) used to move the dopant field from the deformed mesh to the new AMR'd mesh. This is the right transfer when mesh topology changes — native prolongation (CASE 1) does NOT apply here.

**Solution transfer is L2 projection (CASE 2), not native prolongation.** Task 7's CASE 1 only works when the mesh is the same object before/after AMR. Here, oxidation changes the geometry, so the field must be projected across two distinct meshes via `Mx=b` solve.

- Test: oxidation step -> mesh deforms -> quality drops -> remesh with AMR near interface, dose conserved within 1% after L2 projection transfer
- Commit: `"feat: add AMR for moving boundary problems (LevelSetProblem + MultiAppProjectionTransfer pattern)"`

### Task 10: Interface-Aligned Refinement
- Automatic refinement at material interfaces (where element attribute changes between neighbors). Densifies mesh near Si/SiO2, Si/Si3N4 boundaries for accurate segregation
- Test: 2-material mesh -> verify finer elements near interface
- Commit: `"feat: add interface-aligned mesh refinement"`

### Task 11: Uniform Mesh Scaling
- Global mesh refinement/coarsening by factor. Useful for mesh convergence studies. `engine.refineGlobally(factor)`, `engine.coarsenGlobally(factor)`
- Test: 8x8 -> refine 2x -> 16x16. Verify solution converges
- Commit: `"feat: add uniform mesh scaling for convergence studies"`

### Task 12: Mesh Convergence Test
- Solve same diffusion problem on progressively finer meshes. Verify solution (1D profile, dose) converges. Richardson extrapolation for error estimate
- Test: 4x4, 8x8, 16x16, 32x32 -> profiles converge. Dose error decreases
- Commit: `"test: add mesh convergence validation test"`
