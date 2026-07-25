# Diffusion Engine Phase 5: Polysilicon Diffusion (Isotropic + Anisotropic)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add polysilicon grain-based diffusion: isotropic (enhanced D in grain boundaries) and anisotropic (dual-mesh grain interior + boundary with segregation). Add grain growth, interface oxide breakup, epitaxial regrowth, and grain-size-dependent oxidation.

**Depends on:** Phase 3 (CDD, KernelTerm composition), Phase 2 (segregation two-sided interface integrator)

**Manual sources (verified via `iconv -f UTF-16LE -t UTF-8 Manual/sprocess_ug.md`):**
- SProcess §"Diffusion in Polysilicon", body lines **24207–26774** (TOC p.219–227)
- ATHENA §3.1.7 "Grain-based Polysilicon Diffusion Model", body lines **8050–8362** (p. 3-21)

**Solver backend:** Diffusion on the polysilicon dual mesh (grain interior + boundary) uses the same HypreBoomerAMG preconditioner from Phase 1. Segregation at grain boundaries uses the two-sided interface integrator from Phase 2 with the same Hypre-backed Newton solve. Grain growth is an ODE per element (no global solve); the oxide-breakup void model is an explicit level-set-style front (no linear solve).

## Global Constraints

Same as Phase 1. MFEM-gated. Namespace `viennaps`. LLVM style. Reuses `KernelTerm` (Phase 3 Task 0), `DiffusionModel`, `DiffusionEngine`, `MeshAttributes`, `Segregation` (Phase 2).

---

## File Structure

| File | Responsibility | Status (HEAD) |
|------|----------------|---------------|
| `fields/GrainModel.hpp` | Grain growth ODE: dR/dt = k·exp(−Ea/kT)/R^n | exists (1D) |
| `fields/GrainBoundaryMesh.hpp` | Dual-mesh generator (Voronoi + boundary layer) | exists (stub) |
| `fields/models/PolysiliconDiffusion.hpp` | Isotropic + anisotropic poly diffusion | exists (skeleton) |
| `fields/models/PolyOxideBreakup.hpp` | Oxide-breakup void model + epitaxial regrowth | exists (simplified) |
| `fields/models/GrainSegregation.hpp` | Interior↔boundary segregation (FEM, two-sided) | NEW |

---

### Task 0: Polysilicon Model Selection (ADR)

**Required before Task 1.** SProcess exposes two model families for polysilicon:

- **(A) Isotropic / granular** (SProcess `pdbSet PolySilicon Dopant DiffModel Granular`, line 24219) — single effective diffusivity `D_eff = D_gi·f_g + D_gb·(1−f_g)` per element; grain boundary fraction `f_g` tracked as a dataset. Cheap, no topology change. **Default for Phase 5 Tasks 1–2.**
- **(B) Anisotropic / dual-mesh** (SProcess §"Grain Boundary Structure", lines 25165–25510) — separate unknowns `C^g` (grain interior) and `C^gb` (grain boundary) on distinct subdomains; segregation flux couples them at the interior↔boundary interface (`GBMaxDensity`, `GSegInit`). More accurate, more expensive. **Phase 5 Tasks 3–4.**

**Phase 5 default: ship both, with (A) as the default and (B) opt-in via `engine.setPolyModel(Anisotropic)`.** The two share `GrainModel` (Task 1) and `PolyOxideBreakup` (Task 5); only the diffusion assembly differs.

- [ ] **Step 1: Document** the choice in this plan file (done above).
- [ ] **Step 2: No code in Task 0.**

---

### Task 1: GrainModel — Grain Growth Equation

**Files:** Modify `fields/GrainModel.hpp`

**Produces:** `GrainModel<NumericType>` — advances the average grain radius `R(t)` per element via the SProcess grain-growth law (SProcess §"Grain Growth", body lines 24319–24510; ref [13] of the manual):

**Manual equation (SProcess lines 24319–24480, ATHENA §3.1.7 implicit):**

```
dR/dt = (k_0 · exp(−Ea / kT) / R^n) · (1 − (R/R_max)^m)
```

- `k_0` = pre-exponential (`Growth.Rate.0`, default ~1e7 nm/s for columnar poly-Si)
- `Ea` = activation energy (`Growth.Rate.E`, ~2.5 eV)
- `n` = grain-growth exponent (`GrainSizeFactor`, typically 2 — parabolic, matches Burke–Turnbull)
- `R_max` = saturation grain size (`GrainSize` cap, set by film thickness / nucleation density)
- `m` = soft-stop exponent (avoids overshoot past `R_max`)

Columnar vs cubic grain shape (SProcess Figure 26, line 24340): columnar grains extend through the full film thickness `d`; cubic grains have `L ≈ d`. The grain-boundary volume fraction follows from geometry:

```
f_gb = 1 − f_g = 1 − ((R − δ_gb/2) / R)^D      (D = problem dim, δ_gb = boundary thickness)
```

**Assembly:** grain growth is a per-element ODE — no global FEM solve. Store `R` as a `mfem::QuadratureFunction` (one value per element, updated each time step by `GrainModel::advance(dt, T_field)`). The diffusion assembly in Tasks 2/4 reads `f_g(R)` at each quadrature point.

**MOOSE reference (verified):** `framework/include/auxkernels/MaterialStdVectorAux.h` — pattern for exposing a per-element scalar (here `R`) as a field for downstream kernels to consume. Mirror this: `GrainModel` exposes `R` as a `GridFunction` on a piecewise-constant (L2) FESpace.

- [ ] **Step 1: Write failing test** — `TestGrainGrowth()`: at T=1000°C, R grows; at T=600°C, R ≈ static; verify `dR/dt` matches the formula at one point exactly; verify R saturates at `R_max`.
- [ ] **Step 2: Run to verify failure** → FAIL (current `GrainModel` is 1D-only and lacks the saturation term)
- [ ] **Step 3: Implement** — per-element ODE with the saturation term; expose `R` as an L2 `GridFunction`.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add GrainModel with Burke–Turnbull grain growth + saturation (SProcess 24319–24480; MOOSE MaterialStdVectorAux pattern)"`

---

### Task 2: Isotropic Polysilicon Diffusion

**Files:** Modify `fields/models/PolysiliconDiffusion.hpp`

**Produces:** `PolysiliconDiffusion<NumericType>` (isotropic mode, SProcess `DiffModel Granular`) — single species `C`, effective diffusivity that depends on the grain-boundary fraction from Task 1.

**Manual equation (SProcess §"Isotropic vs anisotropic", lines 24219–24310; ATHENA §3.1.7 eq. 3-66 at lines 8108–8136):**

```
D_eff(C, T, R) = D_gi(C, T) · f_g(R) + D_gb(C, T) · (1 − f_g(R))
```

- `D_gi` = grain-interior diffusivity = crystalline-Si dopant D (from Phase 2 Fermi / Phase 10 DB)
- `D_gb` = grain-boundary diffusivity (typically 10²–10⁴× `D_gi`; `Dgb.F22 > Dgb.F11` anisotropy in SProcess line 25407, but isotropic mode uses the scalar average)

As grains grow (`R ↑`), `f_g ↑`, so `D_eff ↓` toward `D_gi`. This couples diffusion to grain growth: Task 1's `R(t)` drives the diffusivity field.

**Assembly (KernelTerm composition, Phase 3 Task 0; quadrature-point D per Phase 2 Task 3 `DiffusivityMaterial`):**

```cpp
// Per quadrature point q:
const NumericType R_q = grainModel_->R(elementIndex, q);
const NumericType fg  = std::pow((R_q - 0.5*deltaGB_) / R_q, D_);
const NumericType Dgi = crystallineD_->D(C_q, T_q);   // Fermi / Constant
const NumericType Dgb = boundaryD_->D(C_q, T_q);      // typically Constant
const NumericType Deff = Dgi*fg + Dgb*(1.0 - fg);
// K_ij += Deff * grad_phi_j · grad_phi_i * JxW
```

**MOOSE references (verified):**
- `framework/include/kernels/MatDiffusionBase.h` — the canonical "diffusion coefficient from a Material property" kernel; this is exactly the pattern. `DiffusivityMaterial` (Phase 2) plays the Material role; `PolysiliconDiffusion` is the kernel.
- `framework/include/functormaterials/PiecewiseByBlockFunctorMaterial.h` — for the per-block `D_gb` assignment (different grain-boundary diffusivities per dopant).

- [ ] **Step 1: Write failing test** — `TestPolyIsotropic()`: (a) fine-grain poly (`R` small) → `D_eff > D_gi`; (b) coarse-grain / annealed (`R → R_max`) → `D_eff → D_gi`; (c) dopant dose conserved under zero-flux BCs.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — `PolysiliconDiffusion` reads `R` from `GrainModel`, computes `D_eff` per QP, assembles stiffness.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add isotropic polysilicon diffusion with grain-growth-coupled D_eff (SProcess Granular model 24219; MOOSE MatDiffusionBase pattern)"`

---

### Task 3: GrainBoundaryMesh — Dual-Mesh Generation

**Files:** Modify `fields/GrainBoundaryMesh.hpp`

**Produces:** `GrainBoundaryMesh<NumericType, D>` — generates the dual-mesh for the anisotropic model (B): a tessellation of grain centers (Voronoi in 2D, constrained Voronoi in 3D) with explicit grain-boundary elements between grains. SProcess §"Grain Boundary Structure" (lines 25165–25253).

**Algorithm:**

1. Sample `N_grain` centers within the polysilicon region (columnar → project to 2D in the film plane; cubic → 3D).
2. Compute the Voronoi diagram of the centers (uses `boost::polygon` / `CGAL` 2D Voronoi, or for 3D `Voro++`). Each Voronoi cell = one grain interior subdomain.
3. Insert a thin "boundary layer" strip of width `δ_gb` along each Voronoi edge → these elements get a distinct attribute (`Poly_GB` vs `Poly_GI`).
4. Emit an `mfem::Mesh` with per-grain attributes; populate `MeshAttributes` (Phase 1 Task 1) with `Poly_GI`, `Poly_GB`, plus the host material name.

**Manual reference (SProcess §"Grain Boundary Structure", lines 25165–25253):** SProcess internally generates the dual mesh from `GSize` (grain size dataset) and `GBVolShare` (boundary volume fraction). The Voronoi step is what `GrainBoundaryMesh` mirrors.

**MOOSE references (verified):**
- `framework/include/meshgenerators/XYDelaunayGenerator.h` — 2D Delaunay/Voronoi mesh generation; **directly portable** for columnar grains.
- `framework/include/meshgenerators/FillBetweenPointVectorsGenerator.h` — boundary-layer strip generation; matches the `δ_gb`-width grain-boundary layer.
- `framework/include/meshgenerators/SurfaceSubdomainsDelaunayRemesher.h` — for 3D surface-constrained remeshing if cubic grains are needed.

**Scope decision:** ship 2D columnar Voronoi (the SProcess default) in Task 3. Defer 3D cubic grains (uses `XYZDelaunayGenerator` + `Voro++`) — note in commit.

- [ ] **Step 1: Write failing test** — `TestGrainBoundaryMesh()`: N=16 grains → verify exactly 16 interior subdomains; verify every grain-boundary element touches two distinct interior subdomains; verify total GB area ≈ `GBVolShare × film_area` within 5%.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — Voronoi via `boost::polygon` (already a ViennaPS dep through ViennaCore) + boundary-layer strip insertion.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add GrainBoundaryMesh dual-mesh via Voronoi + boundary layer (SProcess Grain Boundary Structure; MOOSE XYDelaunayGenerator pattern)"`

---

### Task 4: Anisotropic Polysilicon Diffusion + Grain Segregation

**Files:** Modify `fields/models/PolysiliconDiffusion.hpp` (anisotropic mode), create `fields/models/GrainSegregation.hpp`

**Produces:** `PolysiliconDiffusion` (anisotropic mode, SProcess dual-mesh model) + `GrainSegregation` — two species `C^gi`, `C^gb` on the Task 3 dual mesh, each with its own diffusivity, coupled by an interior-face segregation flux.

**Manual equations (SProcess §"Segregation Between Grain Interior and Boundaries", lines 25462–25600; ATHENA §3.1.7 eq. 3-67 at lines 8168–8207):**

Interior (`C^gi`, attribute `Poly_GI`):

```
∂C^gi/∂t = ∇·(D_gi(T) ∇C^gi)              (crystalline-Si dopant D)
```

Boundary (`C^gb`, attribute `Poly_GB`):

```
∂C^gb/∂t = ∇·(D_gb(T) ∇C^gb)
```

Segregation flux at every `Poly_GI : Poly_GB` interface (dynamic form; equilibrium gives the segregation coefficient `m`):

```
J_seg = k_seg · C^gi − k_deseg · C^gb
       with  m(T) = k_seg / k_deseg  =  exp((−ΔG_seg)/(kT))
```

`ΔG_seg < 0` (dopants preferentially occupy grain-boundary sites) → `m > 1` → grain boundaries enriched. `GBMaxDensity` (`GBMaxDensity` parameter, SProcess line 25480) caps `C^gb` at the grain-boundary site density; above the cap, the segregation flux shuts off (Langmuir-like saturation).

**Assembly — the segregation flux is a two-sided interior-face residual (Phase 2 `Segregation` extended to grain boundaries):**

```cpp
// For each interior face f between Poly_GI (side 1) and Poly_GB (side 2):
// Residual contributions (test functions φ_1 on side 1, φ_2 on side 2):
R_1 +=  (k_seg · C^gi - k_deseg · C^gb) · φ_1 · |f|     // loss from interior
R_2 += -(k_seg · C^gi - k_deseg · C^gb) · φ_2 · |f|     // gain in boundary (mass-conserving)
// Langmuir cap: replace C^gb in the gain term with min(C^gb, GBMaxDensity)
```

This is exactly the Phase 2 `Segregation` two-sided interface integrator pattern (which itself mirrors MOOSE `InterfaceReaction`); the only addition is looping over `Poly_GI : Poly_GB` faces instead of `Si : SiO₂` faces.

**MOOSE references (verified):**
- `framework/include/interfacekernels/InterfaceReaction.h` — `k·[c₁c₂]` interface reaction; **the exact two-sided pattern**. Phase 2's `Segregation` already mirrors this; `GrainSegregation` is a thin attribute-filtered specialization.
- `modules/chemical_reactions/include/materials/LangmuirMaterial.h` — for the `GBMaxDensity` saturation cap; use `MollifiedLangmuirMaterial` for Newton differentiability.

**Coupling to grain growth (Task 1):** as grains grow, grain boundaries shrink and some `Poly_GB` elements vanish (grain coalescence). The simplest model: when `R > R_coalesce`, relabel `Poly_GB → Poly_GI` (subdomain relabeling — Phase 4 Task 0 idiom A, `ElementSubdomainModifier` pattern) and project `C^gb` onto the merged `C^gi` via `IntegralPreservingFunctionIC`. Defer full coalescence dynamics; ship the static dual-mesh segregation first.

- [ ] **Step 1: Write failing test** — `TestPolyAnisotropicSegregation()`: (a) equilibrium: `m > 1` → `C^gb / C^gi → m` at long time; (b) fast `D_gb`, slow `D_gi` → dopant "channels" along the GB network (1D profile shows fast lateral spread); (c) mass conservation `∫(C^gi·f_g + C^gb·(1−f_g)) = const` within 0.5%.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — register two species on the dual mesh; `GrainSegregation` as a two-sided interior integrator filtered to `Poly_GI:Poly_GB` faces; reuse Phase 2 `Segregation`'s Newton/Picard machinery.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add anisotropic polysilicon diffusion + grain-boundary segregation on dual mesh (SProcess 25462–25600; MOOSE InterfaceReaction + LangmuirMaterial pattern)"`

---

### Task 5: Poly/SiO₂ Interface Oxide Breakup + Epitaxial Regrowth

**Files:** Modify `fields/models/PolyOxideBreakup.hpp`

**Produces:** `PolyOxideBreakup<NumericType>` — models the dissolution of a thin native oxide at the poly-Si/Si substrate interface during high-T anneal, followed by epitaxial regrowth of poly-Si as single-crystal Si. SProcess §"Interface Oxide Breakup and Epitaxial Regrowth", body lines 26147–26400 (refs [17]–[20]).

**Manual model (SProcess lines 26147–26383):**

1. **Oxide breakup** — the interfacial oxide (thickness `d_ox(t)`) dissolves at a rate driven by oxygen out-diffusion into the poly:
   ```
   d(d_ox)/dt = −k_break · exp(−Ea_break / kT) · (d_ox / d_ox,0)
   ```
   The `(d_ox / d_ox,0)` factor makes breakup self-limiting as the oxide thins. Once `d_ox < d_crit`, the interface is "broken" (locally).
2. **Epitaxial regrowth** — where the oxide is broken, poly-Si recrystallizes as an extension of the substrate. The regrowth front velocity:
   ```
   v_ereg = v_ereg,0 · exp(−Ea_ereg / kT) · f(broken)
   ```
   where `f(broken)` is the local broken fraction (0 where oxide intact, 1 where fully broken, smooth in between). Behind the front, material is relabeled `Poly → CrystallineSi` (subdomain relabeling, Phase 4 Task 0 idiom A).

**Coupling:** the regrowth front is tracked as a 1D level-set (or a per-element "crystallized fraction" field for the diffuse version). Dopant in the regrown region inherits the poly's dopant profile (conserved via `IntegralPreservingFunctionIC` projection on relabel).

**Assembly:** oxide breakup and regrowth are explicit per-element ODEs/front advections — no global solve. The only FEM-side effect is the subdomain relabel, which triggers `Poly → CrystallineSi` diffusivity change (D drops from `D_poly` to `D_cSi`) on the relabeled elements.

**MOOSE references (verified):**
- `framework/include/meshmodifiers/ThresholdElementSubdomainModifier.h` — idiom A from Phase 4 Task 0; relabels elements where a threshold field (here `crystallized_fraction`) crosses a value. `POLYNOMIAL_NEARBY` reinit strategy for the dopant field on relabeled elements.
- `modules/level_set/include/kernels/LevelSetAdvection.h` — for the level-set form of the regrowth front (alternative to the per-element ODE).

**Scope:** ship the per-element ODE + threshold relabel form (simplest). The full level-set front tracking is a Phase 4 SPER-level feature — note in commit, defer.

- [ ] **Step 1: Write failing test** — `TestPolyOxideBreakup()`: at T=1000°C, `d_ox` decays; once broken, poly near interface relabels to `CrystallineSi`; dopant dose in the relabeled region preserved within 1%.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — per-element ODEs + `ThresholdElementSubdomainModifier`-style relabel + dose-preserving projection.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add poly/SiO2 oxide breakup + epitaxial regrowth via subdomain relabel (SProcess 26147–26383; MOOSE ThresholdElementSubdomainModifier pattern)"`

---

### Task 6: Grain-Size-Dependent Oxidation Rate

**Files:** Modify `fields/models/PolysiliconDiffusion.hpp` (add oxidation-rate hook) + integration with `psOxidation`

**Produces:** a scalar multiplier on the polysilicon oxidation rate that depends on the local grain size. Smaller grains → more grain boundaries → faster oxidant diffusion along GBs → faster oxidation. SProcess §"Polysilicon Oxidation Rate", lines 26383–26450.

**Manual equation (SProcess line 26383+):**

```
(d_ox/dt)_poly = (d_ox/dt)_cSi · (1 + β · (1 − f_g(R)) / (1 − f_g(R_ref)))
```

- `β` = grain-boundary oxidation enhancement factor (`PolyOxRateFactor`, ~0.5–2)
- `R_ref` = reference grain size at which the poly rate equals the crystalline rate
- As `R → R_max` (large grains, `f_g → 1`), the enhancement vanishes — poly oxidizes like c-Si.

**Coupling:** expose `setOxidationRateModifier(β, R_ref)` reading `R` from `GrainModel` (Task 1). The actual oxidation step is `psOxidation` (existing ViennaPS); this Phase 5 task only provides the multiplier. Same hook pattern as Phase 4 Task 6 (nitrogen oxidation-rate modifier).

**Scope:** ship the multiplier + a unit test with a mock oxidizer (verifies the multiplier formula). Full `psOxidation` integration tested in the Task 7 integration test.

- [ ] **Step 1: Write failing test** — `TestGrainSizeOxidationRate()`: fine-grain poly → rate multiplier > 1; coarse-grain → multiplier → 1.
- [ ] **Step 2: Run to verify failure** → FAIL
- [ ] **Step 3: Implement** — multiplier reads `R` from `GrainModel`, applies the formula.
- [ ] **Step 4: Run to verify pass** → PASS
- [ ] **Step 5: Commit** — `"feat: add grain-size-dependent polysilicon oxidation rate multiplier (SProcess 26383–26450)"`

---

### Task 7: Integration Test — Poly Deposition → Implant → Anneal

**Files:** Modify `tests/diffusion/testDiffusion.cpp`

**Scenario (matches SProcess §4.5 canonical poly-emitter workflow):**

1. Deposit 200 nm polysilicon on Si substrate (with thin native oxide, for Task 5 variant).
2. Implant B at 5 keV, dose 5e15 cm⁻² into the poly.
3. Anneal 900 °C for 30 min — B diffuses through poly, segregates to GBs (anisotropic mode), and (if oxide breakup fires) diffuses into the substrate as a poly-emitter.

**Assertions:**

- (a) B diffuses faster in poly than in a single-crystal-Si control (isotropic `D_eff > D_cSi`).
- (b) In anisotropic mode, B is enriched in grain boundaries (`C^gb / C^gi → m`).
- (c) In the oxide-breakup variant, B crosses into the substrate (junction forms in Si); in the no-breakup variant (thick oxide), B stays in poly.
- (d) Total B dose conserved within 1%.
- (e) **Regression snapshot** to `tests/diffusion/regression/poly_anneal_baseline.csv`.

- [ ] **Step 1: Write the test** with the four assertions + snapshot.
- [ ] **Step 2: Run** → verify all pass; record baseline.
- [ ] **Step 3: Commit** — `"test: add polysilicon deposition+implant+anneal integration test with GB segregation + oxide-breakup variants (SProcess 4.5 poly-emitter workflow)"`

---

## Self-Review Notes

- **Spec coverage:** Phase 5 of spec Section 12 = "Polysilicon (isotropic + anisotropic)." Tasks 1–2 cover isotropic + grain growth; Tasks 3–4 cover anisotropic + grain segregation; Task 5 covers oxide breakup + epitaxial regrowth; Task 6 covers grain-size-dependent oxidation; Task 7 ties them together. Every spec §4.5 row that the gap analysis marked 🟡/❌ now has a manual-cited FEM assembly path.
- **Manual citations:** every task cites SProcess body line ranges (24219, 24319, 25165, 25462, 26147, 26383) and ATHENA §3.1.7 with line ranges — verified via `iconv` UTF-16 decode.
- **MOOSE citations:** every task names a verified MOOSE class with path — `MatDiffusionBase` (Task 2), `XYDelaunayGenerator`/`FillBetweenPointVectorsGenerator` (Task 3), `InterfaceReaction`/`LangmuirMaterial` (Task 4), `ThresholdElementSubdomainModifier` (Task 5), `LevelSetAdvection` (Task 5 alt).
- **Engine integration:** reuses Phase 2 `Segregation` two-sided integrator (Task 4), Phase 3 `KernelTerm` composition (Task 2), Phase 4 Task 0 subdomain-relabel idiom (Tasks 5, 6). No new engine APIs required beyond exposing `R` as an L2 `GridFunction` (Task 1).
- **Gap-analysis rows closed:** Isotropic D_eff ✅ (deepened), Grain growth ✅ (deepened with saturation), Anisotropic dual mesh 🟡→✅, Voronoi tessellation ❌→✅ (Task 3), Interface oxide breakup 🟡→✅ (Task 5), Epitaxial regrowth 🟡→✅ (Task 5), Segregation interior-boundary 🟡→✅ (Task 4, FEM-assembled).
