# Code Review: Diffusion Engine Phases 3–4 + 6–11 (751dc336..c475863b)

**Date:** 2026-07-25
**Scope:** 4 commits, 23 files, ~1638 insertions. "Full-depth FEM" implementations covering Phases 3 (charged/cluster models), 4 (OED/impurities), 6–10 (SiGe/KMC/PDE bridge), and 11 mesh infrastructure (MovingMeshHandler + SolutionTransfer).
**Base:** 751dc336 (fix: address diffusion engine code review findings)
**Head:** c475863b (feat(diffusion): MovingMeshHandler + SolutionTransfer full-depth)

**Commits in range:**
- `0f88a0a` — Phase 3 full-depth: ChargedFermi (quadrature-point D), ChargedPair, ChargedReact, Cluster311, VacancyCluster, ImpurityCluster (BIC), DislocationLoop
- `e13fba3` — Phase 4 full-depth: CarbonDiffusion (C_s-I trapping), ChargedEquilibriumDiffusion, CopperDiffusion (drift+pairing), MobileImpurity (template), OedSource
- `332bc6b` — Phases 6–10 deepen: SiGeDiffusion FEM assembly, KmcAtomisticEngine recombination/clustering events, PdeApi bridge
- `c475863` — MovingMeshHandler + SolutionTransfer (Phase 11 mesh infra)

**Reviewed against:**
- `docs/superpowers/plans/2026-07-20-diffusion-phase{3,4,6,7,10,11}.md` (recently expanded with manual equations, MOOSE citations, per-task acceptance criteria)
- `docs/superpowers/specs/2026-07-20-diffusion-parity-design.md`
- `docs/superpowers/reviews/2026-07-25-spec-gap-analysis.md` (12 gap-row closure claims verified)

---

## Strengths

- **Phase 3 cluster/charged models are genuinely well-implemented.** `ChargedFermiDiffusion::ChargedFermiDCoef` (`include/viennaps/fields/models/ChargedFermiDiffusion.hpp:114-135`) correctly evaluates `getDiffusivity(C, T)` per quadrature point via `conc_->GetValue(T, ip)`, closing the mean-Picard gap the right way. `Cluster311`, `VacancyCluster`, `ImpurityCluster`, `DislocationLoop` all follow a consistent, correct FEM reaction-assembly pattern: per-species `scale` selected by GridFunction pointer identity, residual assembled into `LinearForm` via a custom `Coefficient` that reads all coupled species at QPs, with `std::max(0, ...)` floors on concentrations before they enter rate laws (`Cluster311.hpp:104`, `DislocationLoop.hpp:106`).

- **Engine integration is clean and correct.** `DiffusionEngine::assembleAllSpecies` (`DiffusionEngine.hpp:446-487`) calls all four hooks (`assembleStiffness`/`assembleMass`/`assembleReaction`/`finalizeReaction`) in the right order, applies the `shouldCreateTimeDerivative` gatekeeper (Phase 1 Task 3.5) to prevent double `dC/dt`, and falls back to an identity mass when no model claims the time derivative (`DiffusionEngine.hpp:499-503`). The Robin dose-loss fold-into-K is thoughtful.

- **Coefficient ownership is safe.** Every model stores its QP coefficients in `mutable std::unique_ptr` members cleared at the top of each `assembleReaction`/`assembleStiffness` call (e.g. `ChargedReactDiffusion.hpp:85`, `Cluster311.hpp:80`). This keeps the coefficient alive across `LinearForm::Assemble()` — a subtle MFEM lifetime bug class that was correctly avoided.

- **`SolutionTransfer::integrate` and dose bookkeeping** (`SolutionTransfer.hpp:35-58, 107-122`) genuinely compute the L2-dual dose and enforce the 0.1% conservation invariant with a meaningful `TransferResult::ok` flag. The `meshVolume` fallback for zero-target-dose (`SolutionTransfer.hpp:112-114`) is a nice edge-case guard.

- **`KernelTerms.hpp` MFEM/non-MFEM split** is correct — the two `SupgAdvectionTerm` definitions sit in `#ifdef/#else` branches (`KernelTerms.hpp:198, 302`), so there's no ODR violation; the non-MFEM stubs let unit tests compile.

---

## Gap-Row Closure Verification

The gap analysis (`docs/superpowers/reviews/2026-07-25-spec-gap-analysis.md`) marks 12 rows as recently closed by this code. Each was verified against the actual code (not taken at face value):

| # | Gap Row | Verdict | Evidence |
|---|---|---|---|
| 1 | ChargedFermi quadrature-point D | ✅ **Confirmed closed** | `ChargedFermiDCoef::Eval` reads `conc_->GetValue(T, ip)` per QP (`ChargedFermiDiffusion.hpp:121-129`) |
| 2 | ChargedPair FEM assembly | ✅ **Confirmed closed** | `ChargedPairDCoef` + `assembleStiffness` (`ChargedPairDiffusion.hpp:47-109`) |
| 3 | ChargedReact FEM assembly | ✅ **Confirmed closed** | `ChargedProductCoef` evaluates `effectiveRate` per QP (`ChargedReactDiffusion.hpp:97-123`) |
| 4 | {311}/VC/BIC/DislocationLoop stiffness | ✅ **Confirmed closed** | All four have `assembleReaction` + `assembleMass`; `DislocationLoop` uses supersaturation `max(C_I/C_eq - 1, 0)^p` (`DislocationLoop.hpp:104-108`) |
| 5 | ChargedEquilibrium FEM | 🟡 **Partially closed** | FEM assembly exists but the **physics formula is wrong**: code uses empirical `D0*(1+alpha*n/ni)` (`ChargedEquilibriumDiffusion.hpp:36`) instead of the plan-mandated Fermi charge-state partition `D_eff = D0·g_0/Σg_z` with `g_z = exp((E_F−E_z)/kT)` (Phase 4 Task 4 manual eq.). FEM scaffolding landed; the actual SProcess 4.194 model did not. |
| 6 | Carbon I-trapping FEM | 🟡 **Partially closed** | FEM residual exists (`CarbonDiffusion.hpp:70-94`), but the model implements **only the forward trap** `r = kf*C*I`. The plan's detailed-balance reverse rate `k_r = k_f·C*_I` (SProcess 4.193) is entirely absent — no `kr_` member, no dissociation term. The complex `CI` only ever grows; it cannot reach the equilibrium `C_sI_eq = C_s·C_I/C*_I` the plan requires. |
| 7 | Copper drift+pairing FEM (key novelty) | ❌ **FALSE CLOSURE** | `CopperDiffusion.hpp` implements only `D = D0*(1+beta*C_dopant/ni)` — ion-pairing-enhanced scalar D. There is **no drift term** `J = −D·(∇C + (q/kT)·z·C·E)`, no electric field `E = −∇φ`, no `MixedScalarVectorGradientIntegrator`, no pairing reaction `k_pair·C_Cu·C_A`, no immobile `C_CuA` species (`numSpecies()==1`). The Phase 4 plan explicitly calls this "the first drift term — the one architectural novelty" and requires it to land in K. It did not. The test even disables the only physics present (`cu->setIonPairing(0.0)`, testDiffusion.cpp `TestPhase4FullDepthFem`). |
| 8 | SiGe FEM assembly | 🟡 **Partially closed** | FEM stiffness assembly now exists (`SiGeDiffusion.hpp:66-73`, was 1D explicit), but uses plain Arrhenius `D0·exp(−Ea/kT)`. The plan's defect-mediated `D_inter = D_V*·(C_V/C_V*) + D_I*·(C_I/C_I*)` (SProcess 26790–26850) is not implemented. The wrong physics is being assembled. |
| 9 | KMC recombination/clustering events | 🟡 **Partially closed / false closure on sub-claims** | Recombination genuinely implemented (`KmcAtomisticEngine.hpp:70-87, 159-167`). But: (a) **Event tree is still O(N)** — `step()` rebuilds the full `events` vector every call and does a linear scan to select (`KmcAtomisticEngine.hpp:42-116`); the "O(log N) binary heap" gap row remains false-closed. (b) **Clustering/Dissociation/Pair events not implemented** — enum lists them (`KmcEvent.hpp:10-17`) but engine never creates them; only Hop+Recombine. (c) **Diamond lattice is fake** — `setDiamondNeighbors(true)` just adds 4 body diagonals to a cubic grid (`KmcAtomisticEngine.hpp:135-143`), not 2-interpenetrating FCC sublattices. (d) **No IDW smoothing + L2 projection** in `KmcDeatomize::deatomize` (`KmcAtomisticEngine.hpp:206-221`) — it bins per cell. |
| 10 | PDE API engine integration | 🟡 **Partially closed** | `PdeEquation::applyTo` bridges into `DiffusionPhysics` (`PdeApi.hpp:126-132`), so "PdeEquation cannot drive engine" is no longer literally true. But `buildModels` silently drops `ReactionPdeTerm` (the `dynamic_cast` returns null and skips, `PdeApi.hpp:91-93`); PdeIC is never applied (stored but `applyTo` ignores it); `ResultsExtractor` still operates on raw vectors, not MFEM GridFunctions. |
| 11 | MovingMeshHandler exists with deformation/remesh | 🟡 **Partially closed** | File now exists (was Missing). `relabelAttributes` (idiom A) and `displaceNodes`/`liftFreeSurface` (idiom C) are correctly implemented (`MovingMeshHandler.hpp:38-119`). But the **remesh trigger** (skewness/Jacobian quality → remesh), Laplacian smoothing, and 3D moving boundary — all listed in gap-analysis Section 3 — are absent. The instance member `interfaceOnly_` is dead (static methods don't read it). |
| 12 | SolutionTransfer L2 projection | 🟡 **Partially closed (defensible deviation)** | Cross-mesh transfer with dose conservation genuinely works (`SolutionTransfer.hpp:62-123`, verified by tests at 0.1%). But the implementation uses `target.ProjectCoefficient(SourceSampleCoef)` with `smesh->FindPoints(...)` — a **point-sampling projection hybrid**, NOT the `MultiAppProjectionTransfer::assembleL2()` mass-matrix-solve `Mx=b` pattern the plan (Phase 11 Task 7 CASE 2, Phase 7 Task 5) explicitly cites. `FindPoints` with `warn=false` silently returns `-1` for unmapped points near moving boundaries (`SolutionTransfer.hpp:93-94`), losing mass that the post-rescale only approximately restores. Functional, but a real accuracy/robustness regression vs the cited algorithm on non-matching meshes. |

**Score: 4 confirmed, 7 partially closed, 1 false closure (plus 3 additional false-closure sub-claims inside row 9).**

---

## Issues

### Critical (Must Fix)

**C1. CopperDiffusion drift term is false-closure (gap-row 7).**
- **File:** `include/viennaps/fields/models/CopperDiffusion.hpp:28-33`
- **What's wrong:** The plan (Phase 4 Task 7) makes this the architectural centerpiece — "the first model in the engine with a genuine drift contribution" — with manual equation `J_Cu = −D·(∇C + (q/kT)zCE)`, and self-review note: "the drift term must land in the same stiffness matrix K so the HypreBoomerAMG preconditioner sees it." The implementation is ion-pairing-enhanced scalar D, with no drift, no E-field, no pairing reaction, no immobile pair species.
- **Why it matters:** This was the single architectural novelty the plan called out as propagating a constraint into `assembleStiffness`'s signature. Marking it closed blocks the actual work and ships a model that looks like Cu drift to callers but isn't.
- **How to fix:** Implement the drift integrator (MFEM `MixedScalarVectorGradientIntegrator` or hand-assembled `(q/kT)·z·E·φ_j·φ_i` into K), add the `C_CuA` pair species (`numSpecies()==2`), and add the `k_pair`/`k_diss` reaction. Until then, mark gap-row 7 🟡, not ✅.

**C2. KMC event tree is still O(N) linear scan (gap-row 9 sub-claim).**
- **File:** `include/viennaps/fields/kmc/KmcAtomisticEngine.hpp:42-116`
- **What's wrong:** The plan (Phase 7 Task 4) is explicit: "binary-heap event tree for O(log N) selection (was O(N) linear scan)." The code rebuilds the entire `events` vector every step and does `for (const auto &e : events) acc += e.rate; if (acc >= thresh)`. This is O(N·neighbors) per step — worse than the original, not better.
- **Why it matters:** For production KMC (10⁶+ atoms), this is the difference between seconds and hours. The gap row claims the architectural fix landed; it didn't.
- **How to fix:** Use `std::priority_queue` or a binary-heap-with-indices indexed by site/event, with partial rebuild of only the affected site's events after each apply.

**C3. KMC clustering/dissociation/pairing events not implemented (gap-row 9 sub-claim).**
- **File:** `include/viennaps/fields/kmc/KmcAtomisticEngine.hpp` (entire engine)
- **What's wrong:** The plan (Task 2) requires 5 event classes; enum lists them (`KmcEvent.hpp:10-17`) but engine only constructs `Hop` and `Recombine`. The `{311}` cluster formation — the canonical KMC TED mechanism — is absent.
- **Why it matters:** Without clustering, KMC cannot model the "{311} formation → I release → B TED" sequence that is SProcess Ch.5's canonical validation scenario (Phase 7 Task 9).
- **How to fix:** Add `ClusterEvent` construction in `step()` when I meets I, and `DissociateEvent` for bound clusters; both with Arrhenius rates `ν_0·exp(−(E_m+ΔE)/kT)` and binding-energy `ΔE`.

### Important (Should Fix)

**I1. CarbonDiffusion has no reverse rate (gap-row 6 partial).**
- **File:** `include/viennaps/fields/models/CarbonDiffusion.hpp:33-42, 104-118`
- **What's wrong:** Only forward trap `r = kf*C*I` is implemented. The plan's detailed balance `k_r = k_f·C*_I` (so `C_sI_eq = C_s·C_I/C*_I`) requires a `kr_` member and `−kr·CI` term in the residual. Without it, the complex grows monotonically and cannot equilibrate — the model is structurally unable to reproduce the carbon-suppresses-TED equilibrium that the Phase 4 integration test asserts.
- **How to fix:** Add `kr_`, extend `TrapCoef::Eval` to `scale*(kf*C*I − kr*CI)`, and wire `k_r = k_f * PointDefectEquilibrium::C_I_eq(T, "Si")` by default. The Phase 3 cluster models already show the correct pattern (`Cluster311.hpp:104` does `kf*pow(cI,n) − kr*c311`).

**I2. ChargedEquilibriumDiffusion uses the wrong formula (gap-row 5 partial).**
- **File:** `include/viennaps/fields/models/ChargedEquilibriumDiffusion.hpp:29-37`
- **What's wrong:** Implements `D0*(1+alpha*n/ni)`. The plan (Task 4, SProcess 4.194) requires `D_eff = D0·g_0/Σg_z` with `g_z = exp((E_F−E_z)/kT)` over charge states. The code currently looks like a stub-grade linear fit masquerading as the equilibrium model.
- **How to fix:** Port the `Σ_z Dz_z · exp(−z·eta) / Σ_z exp(−z·eta)` pattern from `ChargedFermiDiffusion::getDiffusivity` (`ChargedFermiDiffusion.hpp:83-103`) — that model already does this correctly; ChargedEquilibrium should reuse the same math.

**I3. SiGeDiffusion uses plain Arrhenius, not defect-mediated (gap-row 8 partial).**
- **File:** `include/viennaps/fields/models/SiGeDiffusion.hpp:29-34, 66-73`
- **What's wrong:** Plan (Task 2) requires `D_inter = D_V*·(C_V/C_V*) + D_I*·(C_I/C_I*)`. The code uses `D0*exp(−Ea/kT)`. A passing FEM test on this validates only "diffusion happens," not the SProcess model.
- **How to fix:** Add `C_I`/`C_V` GridFunction lookups (as `PairDiffusion` does for `C_I`), compute `D_inter` per QP.

**I4. MobileImpurity does not generalize CopperDiffusion (DRY failure).**
- **Files:** `include/viennaps/fields/models/CopperDiffusion.hpp` + `MobileImpurity.hpp`
- **What's wrong:** Plan (Task 8) specifies `MobileImpurity<NumericType, SpeciesTag>` with `using CopperDiffusion = MobileImpurity<double, CopperTag>`. The implementation is a non-templated `MobileImpurity<NumericType>` that duplicates CopperDiffusion line-for-line (~50 lines, same `D0*(1+beta*C/ni)` formula, same `ImpDCoef`/`CuDCoef` classes with different names). Two parallel implementations now exist where the plan wanted one.
- **How to fix:** Make `MobileImpurity` a `SpeciesTag` template (or a parameter struct), and `CopperDiffusion` a typedef. Then both share one code path.

**I5. SolutionTransfer deviates from the cited L2-projection algorithm.**
- **File:** `include/viennaps/fields/SolutionTransfer.hpp:79-104`
- **What's wrong:** Plan (Phase 11 Task 7 CASE 2, Phase 7 Task 5) cites `MultiAppProjectionTransfer::assembleL2()` (build `M_ij`, build `b_i`, solve `Mx=b`). The code uses `target.ProjectCoefficient(SourceSampleCoef)` with `FindPoints(pts, el, ips, warn=false)`. Points the locator can't find return 0 silently (`SolutionTransfer.hpp:93-94`), which can lose mass near moving boundaries; the post-rescale only roughly compensates.
- **Why it matters:** The mass-matrix L2 approach has no point-location failure mode and is the algorithm the plan cited. The deviation is defensible but introduces a real edge-case regression.
- **How to fix:** Either (a) implement the cited `Mx=b` solve with a `GridFunctionCoefficient(source)` RHS integrator on the target mesh (no FindPoints needed), or (b) document why point-sampling is preferred and add a hard error when `FindPoints` returns `-1` for any QP (fail-loud, not silent-zero).

**I6. Tests verify "runs + dose conserved + monotonic," not physics.**
- **File:** `tests/diffusion/testDiffusion.cpp` (`TestPhase3FullDepthFem`, `TestPhase4FullDepthFem`, etc.)
- **What's wrong:** The plan's quantitative pass criteria are almost entirely absent. Examples:
  - Carbon test checks `CI > 0`, `I1 < I0` — never checks `ΔC_s + ΔC_sI ≈ 0` (Phase 4 Task 5 requirement).
  - Copper test sets `setIonPairing(0.0)` — explicitly disabling the only physics — then asserts dose conservation. The drift-up-gradient test (Task 7 Step 1) is absent.
  - Cluster311 test checks `C1 > 1.0` — never checks dissociation at long time (Task 6 requirement).
  - No test checks `C_sI_eq = C_s·C_I/C*_I` equilibrium, `C^gb/C^gi → m` segregation, or `{311}` Ostwald ripening.
  - Dose tolerance is 5% (`< 0.05`) throughout; the plan specifies 0.1% in several places.
- **Why it matters:** "Runs without crashing" tests give false confidence. Several false closures above (C1, C2, C3) would have been caught by physics-verification tests.
- **How to fix:** Add tests with known analytic equilibria: initialize C_s + I, run to steady state, assert `CI/(C_s·C_I) → 1/C*_I`. Add an O(N) vs O(log N) KMC scaling test (double the lattice, assert near-constant per-step cost). Add a Cu drift test with a known E-field.

### Minor (Nice to Have)

**M1. `MovingMeshHandler::interfaceOnly_` is dead state.** The class has one instance member (`MovingMeshHandler.hpp:33, 123`) that no static method reads. Either remove it, or make the methods non-static and honor the flag (cleaner API: callers do `handler.setInterfaceAdjacentOnly(true); handler.relabel(...)` instead of passing the bool per call).

**M2. `OedSource::relabelOxidized` flips globally, not interface-adjacent.** `OedSource.hpp:63-69` flips every `fromAttr` element when progress crosses threshold — the comment admits "Production: flip only elements adjacent to the interface." Use `MovingMeshHandler::relabelAttributes(mesh, ..., interfaceAdjacentOnly=true)` instead.

**M3. `SolutionTransfer::FindPoints` signature may be MFEM-version-sensitive.** The `FindPoints(pts, el, ips, warn)` form (`SolutionTransfer.hpp:92`) is the older MFEM C API; MFEM 4.4+ moved toward `FindPointsGFCoef` / returning a struct. Not verifiable without the vendored MFEM version, but worth confirming against `f:/dev/mfem/build` per AGENTS.md.

**M4. Duplicated `kB` constant.** The Boltzmann constant `8.617333262145e-5` is re-typed in `SiGeDiffusion.hpp:30,126`, `PdeApi.hpp:191,219`, `KmcEvent.hpp:34,39`. Extract to a single `constexpr` in a shared header.

**M5. Phase 4 plan cites SProcess eq. 4.193/4.194/4.195 body line ranges** (19364–19564, 19733–19910, 19869+). These are credible but I could not independently verify against the manual (not in the review checkout). Implementer should double-check the reverse-rate `k_r = k_f·C*_I` form against the actual manual text before the I1 fix — the plan asserts it, but I1's fix depends on it being correct.

---

## Recommendations

1. **Reopen gap-analysis rows 5, 6, 7, 8, 9 (sub-claims), 10, 11, 12 as 🟡** to reflect actual state. The current ✅ marks overstate closure and will mislead downstream phases (e.g., Phase 5+ builds that consume CopperDiffusion drift, Phase 9 OED workflows that assume full KMC clustering).

2. **Add a "physics-formula" column to the gap analysis**, separate from the "FEM assembly" column. Several rows here got FEM scaffolding without the cited physics; tracking these separately would prevent the recurring "stub formula behind a real-looking assembly" pattern.

3. **Before any more model headers land, write the equilibrium test first.** The plan specifies quantitative equilibria (`C_sI_eq = C_s·C_I/C*_I`, `C^gb/C^gi → m`, `D_inter` increasing with `x_Ge`). These are cheap to write and would have caught C1/C2/C3/I1/I2/I3 immediately.

4. **Consolidate MobileImpurity + CopperDiffusion** before adding Na/K/Fe tags — otherwise the duplication compounds 5×.

5. **Consider porting `SolutionTransfer` to the cited `Mx=b` form** (or adding it as a second path). The point-sampling approach is fine for refinement-only topology changes but is the wrong tool for oxidation-driven remesh where boundary points genuinely don't exist on the source mesh.

---

## Assessment

**Ready to merge?** **No — with fixes.**

**Reasoning:** Phase 3 (rows 1–4) is genuinely well-done and merge-worthy. But Phase 4 row 7 (CopperDiffusion drift) is a false closure of the plan's single most-called-out architectural novelty, rows 5/6/8 ship FEM scaffolding with stub physics behind it, and the KMC engine's O(N) event tree + missing clustering events (row 9) contradict the plan's explicit closure claims. The test suite consistently verifies "runs + dose conserved" rather than the plan's quantitative physics criteria, so it would not catch these regressions. The honest path is to land Phase 3 + the Phase 11 mesh infrastructure now, reopen the over-stated gap rows, and rework Copper/Carbon/ChargedEquilibrium/SiGe/KMC against the actual manual equations before marking them closed.

---

## Relevant file paths (all absolute)

- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\ChargedFermiDiffusion.hpp` (row 1 — good)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\ChargedPairDiffusion.hpp` (row 2 — good)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\ChargedReactDiffusion.hpp` (row 3 — good)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\Cluster311.hpp`, `VacancyCluster.hpp`, `ImpurityCluster.hpp`, `DislocationLoop.hpp` (row 4 — good)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\ChargedEquilibriumDiffusion.hpp` (row 5 — wrong formula)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\CarbonDiffusion.hpp` (row 6 — missing `k_r`)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\CopperDiffusion.hpp` (row 7 — **false closure**)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\MobileImpurity.hpp` (DRY failure)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\models\SiGeDiffusion.hpp` (row 8 — plain Arrhenius, not defect-mediated)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\kmc\KmcAtomisticEngine.hpp` (row 9 — O(N) scan, only Hop+Recombine)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\PdeApi.hpp` (row 10 — reaction terms dropped, ICs not applied)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\MovingMeshHandler.hpp` (row 11 — no remesh trigger)
- `D:\dev\ViennaPS_mod\include\viennaps\fields\SolutionTransfer.hpp` (row 12 — point-sampling, not cited `Mx=b`)
- `D:\dev\ViennaPS_mod\tests\diffusion\testDiffusion.cpp` (tests verify "runs", not physics)
