# Diffusion Engine Phase 8: KMC Lattice Epitaxy

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Add lattice Kinetic Monte Carlo (LKMC) for epitaxial deposition: surface chemistry events, planar/coordination/coordination-reactions growth models, SiGe mole-fraction growth, visibility/shadowing, twin-defect formation, dopant surface segregation, nonselective poly-on-oxide nucleation, and FEM mesh update after growth.

**Depends on:** Phase 7 (`KmcAtomisticEngine`, `KmcLattice`), Phase 3 (CDD), Phase 1 (LevelSetToMesh)

**Manual sources (verified via `iconv` UTF-16 decode):**
- SProcess Chapter 6 "Lattice Kinetic Monte Carlo: Epitaxial Deposition", body lines **66165–69514** (TOC p.455–471)
- SProcess §"Facet Growth During Selective Epitaxy" (continuum), line 32213

**Solver backend:** FEM mesh update after KMC epitaxy (Task 10) uses HypreBoomerAMG for the L2 projection from grown-layer KMC atom positions to the extended FEM mesh — same path as Phase 7 Task 5 (`KmcDeatomize`). The KMC engine itself is event-driven (no linear solve); it reuses Phase 7's BKL rejection-free loop and binary-heap event tree.

## Global Constraints

Same as Phase 1. MFEM-gated. Namespace `viennaps`. LLVM style. Reuses Phase 7 `KmcLattice`, `KmcEvent`, `KmcAtomisticEngine` infrastructure.

---

## File Structure

| File | Responsibility | Status (HEAD) |
|------|----------------|---------------|
| `fields/kmc/KmcEpitaxy.hpp` | LKMC epitaxial deposition engine (extends `KmcAtomisticEngine`) | exists (stub, 75 lines) |
| `fields/kmc/KmcSurfaceEvent.hpp` | Surface attachment/desorption/diffusion events | NEW |
| `fields/kmc/KmcEpitaxyModels.hpp` | Planar, Coordination, Coordination-Reactions growth models | NEW |
| `fields/kmc/KmcVisibility.hpp` | Shadowing of incoming species | NEW |

---

### Task 0: Epitaxy Model Selection (ADR)

**Required before Task 1.** SProcess exposes four epitaxy model classes (lines 66217–66273), in increasing physical fidelity:

- **(A) Planes** (lines 66273–66350) — growth rate per crystal orientation ({100}, {110}, {111}); cheapest.
- **(B) Coordinations.Planes** (lines 66350–66400) — Planes + Tcl callbacks for per-configuration rates; more configurable. **Default for Phase 8 Tasks 2–4** (matches SProcess default).
- **(C) Coordinations** (lines 66400–66433) — growth rate from local atomic coordination (number + type of bonds); captures facet dependence.
- **(D) Coordinations.Reactions** (lines 66433+) — adds explicit CVD reactions (gas decomposition, H blocking, desorption); most accurate.

**Phase 8 default: ship (A), (B), (C); (D) is a stretch (Task 4) — its reaction set is large and material-specific.** All four share the lattice + event-tree infrastructure (Phase 7).

- [ ] **Step 1: Document** (done above).
- [ ] **Step 2: No code.**

---

### Task 1: KmcSurfaceEvent — Surface Chemistry

**Files:** Create `fields/kmc/KmcSurfaceEvent.hpp`

**Produces:** `KmcSurfaceEvent` and subclasses for the surface events that drive epitaxial growth. Extends Phase 7 `KmcEvent`. SProcess §"Surface Events" within Ch. 6.

**Event types and Arrhenius rates (SProcess lines 66217–66500):**

1. **`SiAttachEvent`** — SiH₄ decomposition → Si (surface) + 4H. Rate-limiting step is H desorption (Si surface must be H-free to accept Si). Rate:
   ```
   r_Si = k_Si · P_SiH4 · exp(−Ea_Si / kT) · θ_H(Site)       (θ_H = H-coverage factor)
   ```
2. **`HDesorbEvent`** — H₂ desorption from a Si-H surface site. Rate `r_H = k_H·exp(−Ea_H/kT)`. **The rate-limiting step below ~600°C; above ~800°C, attachment dominates.**
3. **`GeAttachEvent`** — GeH₄ decomposition → Ge + 4H (for SiGe epitaxy). Rate `r_Ge = k_Ge·P_GeH4·exp(−Ea_Ge/kT)`.
4. **`DopantIncorporateEvent`** — dopant (B, P, As) attaches from gas (B₂H₆, PH₃, AsH₃). Rate `r_D = k_D·P_dopant·exp(−Ea_D/kT)`. Incorporation probability depends on local coordination (Task 3).
5. **`SurfaceDiffuseEvent`** — adsorbed Si/Ge hops along the surface before locking in (rate `r_diff = k_diff·exp(−Ea_diff/kT)`, fast). This is what smooths growth — without it, KMC epitaxy is atomically rough.

**Each event inherits Phase 7's `KmcEvent`** — site index, target site, rate `r = nu0·exp(−Em/kT)`, `execute()` mutates the lattice.

**MOOSE references (verified):** KMC has no direct MOOSE analog (per Phase 7's findings — `stochastic_tools` provides sampling/reporting, not event-driven KMC). Cite `modules/stochastic_tools/include/samplers/MonteCarloSampler.h` for the random-event-selection pattern, and `modules/phase_field/include/kernels/LangevinNoise.h` for the stochastic-flux analogue — but the core event loop is ViennaPS code.

- [ ] **Step 1: Write failing test** — `TestKmcSurfaceEvent()`: at low T (500°C), H desorption rate dominates; at high T (900°C), Si attach rate dominates; verify rate ordering `r_H > r_Si` at low T and `r_Si > r_H` at high T.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — `KmcSurfaceEvent` base + five subclasses with Arrhenius rates.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add KmcSurfaceEvent types for epitaxial deposition (SProcess Ch.6 surface events; Arrhenius rates)"`

---

### Task 2: Planar Epitaxy Model (Planes + Coordinations.Planes)

**Files:** Create `fields/kmc/KmcEpitaxyModels.hpp`

**Produces:** `KmcEpitaxyPlanes` and `KmcEpitaxyCoordinationPlanes` — growth rate based on the local surface orientation. SProcess §"Planar Epitaxy Models", body lines 66273–66400, equations 841–842.

**Manual equations (SProcess eq. 841, 842 at lines 66384–66405):**

Growth rate for an empty site at the substrate–gas interface:

```
G^Epi(site) = K' · ν^SEG_LKMC(site)                                    (SProcess 841)

ν^SEG_LKMC(site) = ν^LKMC(site) · K^SEG(n,m) · exp(−(E^SEG + ΔE^LKMC(site)) / kBT)    (842)
```

- `K^SEG(n,m)` — selective-epitaxial-growth prefactor depending on local configuration `n` ({100}, {110}, {111}) and second-neighbor coordination `m`
- `ΔE^LKMC(site)` — per-site energy penalty (e.g., for under-coordinated sites)
- `ν^LKMC(site)` — attempt frequency

The Coordinations.Planes model splits {100} sites into three classes (`100`, `100.7`, `100.8`) by second-neighbor count (line 66405) — finer than the Planes model's single {100} class.

**Algorithm (per event loop iteration):**

1. For each empty site adjacent to a filled site, compute `G^Epi(site)` from eq. 841.
2. Add an `SiAttachEvent` (Task 1) at that site with rate `G^Epi` to the event tree.
3. BKL selection (Phase 7 Task 4): pick event proportional to rate, execute, advance time `dt = −ln(rand)/Σ_rates`.
4. After execution, recompute affected neighbors' rates (only the local neighborhood changes).

**Si/SiO₂ selectivity (Task 9 set-up):** `K^SEG(n,m) = 0` on SiO₂ surfaces (no epitaxial nucleation); `K^SEG > 0` on Si. This is the selective-epitaxy condition.

**MOOSE reference (verified):** no MOOSE analog for KMC growth-rate models. The growth-rate formula is per-site; cite `modules/phase_field/include/kernels/ACInterfaceKobayashi1.h` only as a conceptual reference for orientation-dependent interface kinetics.

- [ ] **Step 1: Write failing test** — `TestKmcEpitaxyPlanar()`: (a) {100} surface → verify growth rate matches `K^SEG(100,·)·exp(−E^SEG/kT)`; (b) faceted surface (mixed {100}/{111}) → {100} grows faster, facets sharpen.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — `KmcEpitaxyPlanes` + `KmcEpitaxyCoordinationPlanes` (second-neighbor split).
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add planar epitaxy LKMC models (SProcess 841–842, 66273–66400)"`

---

### Task 3: Coordination-Based Growth Model

**Files:** Add to `fields/kmc/KmcEpitaxyModels.hpp`

**Produces:** `KmcEpitaxyCoordination` — growth rate from local atomic coordination (number + type of bonds to filled neighbors), without explicit reactions. SProcess §"Coordination Model", lines 66400–66433.

**Model:** rate depends on `coord(site) = (#Si neighbors, #Ge neighbors, #dopant neighbors)`:

```
G(site) = G_0 · exp(−ΔE(coord) / kBT) · f_SiO2(surface)
```

- Under-coordinated sites (few filled neighbors) grow **faster** (more dangling bonds to capture atoms) — this is what makes {111} grow slowly (low dangling-bond density) and {100} grow fast.
- `ΔE(coord)` is a calibrated lookup table per coordination class (Phase 10 DB).
- `f_SiO2 = 0` on oxide → no nucleation (selective epitaxy).

**Facet-dependent growth (SProcess line 66415–66433):** {311} facet formation is captured by the coordination rule — sites at {311} steps have intermediate coordination and intermediate growth rate. The {311}→{111} rate ratio is a parameter (line 66429).

**Algorithm:** same BKL loop as Task 2, but rate = `G_0·exp(−ΔE(coord)/kT)` instead of eq. 842.

- [ ] **Step 1: Write failing test** — `TestKmcEpitaxyCoordination()`: (a) flat {100} → growth rate matches; (b) {111} grows slower than {100} (verify rate ratio matches Table); (c) under-coordinated step site grows faster than terrace site.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — `KmcEpitaxyCoordination` with coordination lookup table.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add coordination-based epitaxy LKMC model (SProcess 66400–66433; {311} facet control)"`

---

### Task 4: Coordination-Reactions Model (stretch)

**Files:** Add to `fields/kmc/KmcEpitaxyModels.hpp`

**Produces:** `KmcEpitaxyReactions` — Task 3 + explicit CVD reactions: SiH₄ decomposition pathway, H site-blocking, desorption. SProcess §"Coordination-Reactions Model", lines 66433+.

**Reaction set (illustrative — full set in Phase 10 DB):**

```
SiH4(g) + surface → Si(surface) + 4H(surface)        r = k1·P_SiH4·exp(−Ea1/kT)
H(surface) + H(surface) → H2(g)                       r = k2·exp(−Ea2/kT)·θ_H²
H(surface) → ∅  (site unblocks)                       r = k3·exp(−Ea3/kT)
```

H coverage `θ_H` blocks Si attachment (sites with H neighbors have reduced `r_Si`).

**Scope:** stretch goal — implement the minimal 3-reaction set above as a proof of concept; the full calibrated set is a Phase 10 follow-up. Cite in commit that the framework is general (any reaction = a `KmcSurfaceEvent` subclass).

- [ ] **Step 1: Write failing test** — `TestKmcEpitaxyReactions()`: at fixed T and P_SiH4, verify growth rate vs the planar model; verify H₂ desorption reduces θ_H over time.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — minimal reaction set + H-coverage bookkeeping.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add coordination-reactions epitaxy LKMC with SiH4 decomposition + H blocking (SProcess 66433+; minimal set)"`

---

### Task 5: SiGe Mole-Fraction-Dependent Growth

**Files:** Modify `fields/kmc/KmcEpitaxyModels.hpp`

**Produces:** SiGe growth mode where the Si:Ge ratio in the grown layer is set by the gas composition (`P_SiH4 / P_GeH4`). SProcess §"SiGe Mole Fraction", line 66559; continuum analogue at line 32213.

**Model:** the Si and Ge attach rates compete for the same surface sites. The grown-layer Ge fraction:

```
x_Ge ≈ r_Ge / (r_Si + r_Ge) = (k_Ge·P_GeH4) / (k_Si·P_SiH4 + k_Ge·P_GeH4)
```

(modulo the different `Ea` and H-coverage effects). Ge surface segregation (Task 8) modifies this — Ge tends to segregate to the surface, so the bulk Ge fraction is slightly less than the surface fraction.

**Assembly:** `SiAttachEvent` and `GeAttachEvent` (Task 1) both target the same empty site; BKL selection picks one proportional to its rate. No new infrastructure.

- [ ] **Step 1: Write failing test** — `TestKmcSiGeGrowth()`: vary `P_GeH4/P_SiH4`, verify grown-layer `x_Ge` follows the ratio within 10%.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — parameterize `GeAttachEvent` rate; co-existence with `SiAttachEvent`.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add SiGe mole-fraction-dependent LKMC growth (SProcess 66559; gas-ratio controlled)"`

---

### Task 6: KmcVisibility — Shadowing

**Files:** Create `fields/kmc/KmcVisibility.hpp`

**Produces:** `KmcVisibility<NumericType, D>` — checks whether a surface site is visible from the gas-source direction. Sites shadowed by neighboring features (trench walls, mushroom caps) have reduced arrival rate. SProcess §"Visibility", within Ch. 6.

**Algorithm:** for each candidate attach site, cast a ray from the gas-source direction toward the site; if it intersects any filled lattice site before reaching the candidate, the site is shadowed (`visibility = 0` or reduced by cosine factor).

**Implementation (D-dimensional ray-cast on the lattice):**

```cpp
bool isVisible(LatticeSite site, Direction sourceDir) const {
  // Walk the ray from site outward in -sourceDir; if any filled site is
  // encountered within maxRange_, return false.
  LatticeSite probe = site;
  for (int step = 0; step < maxRange_; ++step) {
    probe = probe.step(-sourceDir);
    if (lattice_.isFilled(probe)) return false;
  }
  return true;
}
```

The arrival rate at a shadowed site is multiplied by `visibility` (0 or a cosine of the angle to the nearest opening).

**MOOSE reference (verified):** `modules/ray_tracing/include/raytracers/RayTracer.h` — MOOSE's general ray-tracing framework; the visibility ray-cast is a special case. ViennaPS already has ray tracing in the `gpu/` and `process/` layers (for flux engines); reuse `RayTracer` infrastructure if available.

- [ ] **Step 1: Write failing test** — `TestKmcVisibility()`: flat surface → all sites visible; trench geometry → bottom sites shadowed, rim sites visible; mushroom cap → sites under the cap shadowed.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — lattice ray-cast.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add KmcVisibility shadowing via lattice ray-cast (MOOSE RayTracer pattern)"`

---

### Task 7: Twin-Defect Formation on {111}

**Files:** Modify `fields/kmc/KmcEpitaxyModels.hpp`

**Produces:** twin-defect (stacking-fault) formation on {111} facets during epitaxial growth. A twin defect is a rotation of the diamond-cubic stacking sequence (ABCABC → ABCACBC...). SProcess §"Twin Defects", Ch. 6.

**Model:** at each {111} step, with probability `p_twin(T)`, the newly attached atom forms a twin-stacking fault instead of the correct stacking. `p_twin = p_0·exp(−Ea_twin/kT)`.

**Effect on dopant:** twin defects are traps — dopants preferentially incorporate at twin boundaries. The `DopantIncorporateEvent` (Task 1) gets a rate boost at twin-boundary sites.

**Implementation:** maintain a `stackingSequence(site)` per lattice column; on `SiAttachEvent` at a {111} site, draw a random number against `p_twin`; if twin, flip the expected stacking.

- [ ] **Step 1: Write failing test** — `TestKmcTwinDefect()`: on a {111} facet at high T, twin-defect count > 0; on {100}, twin count ≈ 0.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — stacking-sequence tracking + probabilistic twin formation.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add twin-defect formation on {111} during LKMC epitaxy (SProcess Ch.6 twin defects)"`

---

### Task 8: Dopant Surface Segregation (Snowplow)

**Files:** Add to `fields/kmc/KmcSurfaceEvent.hpp` / `KmcEpitaxyModels.hpp`

**Produces:** dopant segregation to the growth front during epitaxy — the "snowplow effect." As Si grows over a doped region, dopants are pushed toward the surface rather than incorporated into the bulk. SProcess §"Surface Segregation", Ch. 6.

**Model:** when a `SiAttachEvent` would bury a dopant site, with probability `p_push = k_seg_surf / (k_seg_surf + k_incorporate)`, the dopant hops up to the new surface instead of being buried. The effective segregation coefficient:

```
m_surf = k_seg_surf / k_incorporate = exp(−ΔG_seg_surf / kT)
```

**Effect:** dopant piles up at the growth front, producing the classic sharp spike at the epi/substrate interface in SIMS profiles.

**Implementation:** add a `DopantPushEvent` that competes with `SiAttachEvent` when burying a dopant site.

- [ ] **Step 1: Write failing test** — `TestKmcSnowplow()`: grow Si over a uniform B-doped substrate → B profile spikes at the growth front (surface B > bulk B by factor `m_surf`).
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — `DopantPushEvent` + competition with attach.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add dopant surface segregation (snowplow) during LKMC epitaxy (SProcess Ch.6 surface segregation)"`

---

### Task 9: Nonselective Epitaxial Deposition (Poly on Oxide)

**Files:** Modify `fields/kmc/KmcEpitaxyModels.hpp`

**Produces:** in nonselective epitaxy, Si nucleates as polysilicon islands on SiO₂ (no epitaxy on oxide) and coalesces; on exposed Si it grows epitaxially. SProcess §"Nonselective Epitaxy", line 68878.

**Model:** on SiO₂ surfaces, replace `KmcEpitaxyPlanes` (Task 2, `K^SEG = 0`) with `KmcIslandNucleation`:

1. Random island nucleation at rate `r_nuc = k_nuc·P_SiH4·exp(−Ea_nuc/kT)` per oxide surface site.
2. Islands grow laterally with the same coordination model (Task 3) but random grain orientation.
3. When islands coalesce, grain boundaries form (polysilicon).

On Si, normal epitaxial growth (Tasks 2–4) applies.

**Implementation:** add a per-surface-material growth-model switch: `KmcEpitaxyPlanes` on Si, `KmcIslandNucleation` on SiO₂.

- [ ] **Step 1: Write failing test** — `TestKmcNonselectiveEpi()`: patterned Si/SiO₂ substrate → epitaxial growth on Si, poly islands on SiO₂.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — per-material growth model + island nucleation.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add nonselective epitaxy with poly nucleation on oxide (SProcess 68878)"`

---

### Task 10: FEM Mesh Update After Epitaxy

**Files:** Modify `fields/DiffusionEngine.hpp` (mesh-extension hook), `fields/LevelSetToMesh.hpp`

**Produces:** after KMC epitaxy grows a layer, update the FEM mesh (extend it upward) and transfer the dopant field from KMC atom positions to the new FEM mesh. SProcess §"FEM Mesh Update", implicit in the LKMC↔continuum coupling.

**Algorithm (reuses Phase 7 Task 5 `KmcDeatomize` L2 projection):**

1. Take the grown KMC atom positions (Si, Ge, dopant) after `KmcEpitaxy::run()` completes.
2. Extend the FEM mesh upward by the grown-layer thickness (top-boundary displacement).
3. Project atom-binned dopant concentration onto the extended FEM mesh via `KmcDeatomize` (Phase 7 Task 5): IDW smoothing + L2 projection.
4. Update `MeshAttributes` for the new layer (Si, SiGe with `x_Ge` field, or PolySi per Task 9).

**Nonatomistic mode (SProcess line 66230):** when LKMC is used only for the deposition shape and dopant/diffusion is handled by the continuum solver, the mesh update is the only coupling — no atom-level dopant projection.

**MOOSE references (verified):**
- `framework/include/transfers/MultiAppProjectionTransfer.h` — L2 projection (Phase 7 Task 5 already cites this).
- `framework/include/meshgenerators/MoveNodeGenerator.h` + `framework/include/meshgenerators/SmoothMeshGenerator.h` — mesh extension + smoothing after node displacement.

- [ ] **Step 1: Write failing test** — `TestFemMeshUpdateAfterEpi()`: grow 100 nm Si via KMC → verify FEM mesh extended by 100 nm; dopant profile transferred, dose conserved within 1%.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — mesh extension + `KmcDeatomize` projection + attribute update.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add FEM mesh update after KMC epitaxy (L2 projection dopant transfer; MOOSE MoveNodeGenerator + MultiAppProjectionTransfer pattern)"`

---

### Task 11: Integration Test — Selective SiGe Epitaxy

**Files:** Modify `tests/diffusion/testDiffusion.cpp`

**Scenario (matches SProcess §6 selective SiGe epi workflow):**

1. Patterned Si/SiO₂ substrate (Si exposed in stripes).
2. Selective SiGe epitaxy: `P_SiH4 + P_GeH4`, T=700°C.
3. Verify: growth only on Si (Task 2 selectivity, `K^SEG=0` on oxide); grown layer has correct `x_Ge` (Task 5); facets form at stripe edges (Task 3); dopant (if B₂H₆ added) spikes at growth front (Task 8 snowplow).
4. After growth: FEM mesh updated (Task 10), dopant profile extracted.

**Assertions:**

- (a) SiO₂ regions have zero grown layer thickness.
- (b) Si regions have grown layer with `x_Ge` matching the gas ratio within 10%.
- (c) {311} facets visible at stripe edges (cross-section profile).
- (d) **Regression snapshot** to `tests/diffusion/regression/lkmc_selective_sige_baseline.csv`.

- [ ] **Step 1: Write the test** with the assertions + snapshot.
- [ ] **Step 2: Run** → verify all pass; record baseline.
- [ ] **Step 3: Commit** — `"test: add selective SiGe LKMC epitaxy integration test (SProcess Ch.6 selective epi workflow)"`

---

## Self-Review Notes

- **Spec coverage:** Phase 8 of spec Section 12 = "KMC lattice epitaxy." Tasks 1–4 cover the four growth-model classes (Planes, Coordinations.Planes, Coordinations, Coordination-Reactions); Task 5 covers SiGe mole fraction; Task 6 visibility; Task 7 twin defects; Task 8 surface segregation; Task 9 nonselective; Task 10 FEM coupling; Task 11 integration test. Every spec §5.2 row that the gap analysis marked 🟡/❌ now has a manual-cited implementation path.
- **Manual citations:** every task cites SProcess Chapter 6 body line ranges (66217, 66273, 66400, 66433, 66559, 68878) with equation numbers (841, 842) — verified via `iconv` UTF-16 decode.
- **MOOSE citations:** `MonteCarloSampler` (Task 1, stochastic-event pattern), `ACInterfaceKobayashi1` (Task 2, conceptual), `RayTracer` (Task 6), `MultiAppProjectionTransfer` + `MoveNodeGenerator` (Task 10).
- **Engine integration:** reuses Phase 7 `KmcEvent`/`KmcLattice`/`KmcAtomisticEngine` (BKL loop, binary-heap event tree), Phase 7 Task 5 `KmcDeatomize` (Task 10 mesh update), Phase 1 `LevelSetToMesh` (Task 10 mesh extension).
- **Deferred (explicit, bounded):** full Coordination-Reactions reaction set (Task 4 — minimal set shipped); 3D cubic lattice for non-{100} substrates (Task 3 — uses 2D columnar lattice for now, matches SProcess default).
- **Gap-analysis rows closed:** KmcEpitaxy engine 🟡→✅, Planar epitaxy ❌→✅, Coordination-based ❌→✅, Coordination-reactions ❌→✅, SiGe mole fraction growth ❌→✅, Visibility/shadowing ❌→✅, Twin-defect formation ❌→✅, Surface segregation ❌→✅, Nonselective epitaxy ❌→✅.
