# Diffusion Engine Phase 7: KMC Atomistic + Deatomize/Atomize Coupling

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Add atomistic Kinetic Monte Carlo (KMC) for point-defect diffusion, with deatomize/atomize transfer to couple with the continuum FEM engine. Replaces the toy simple-cubic BKL loop (current `KmcAtomisticEngine`, 196 lines) with a Si diamond-cubic lattice and a full event set (recombination, clustering, dissociation, pairing).

**Depends:** Phase 3 (CDD, clustering)

**Manual sources (verified via `iconv` UTF-16 decode of `Manual/sprocess_ug.md`):**
- SProcess **Chapter 5 "Atomistic Kinetic Monte Carlo"**, body lines **47557–66164** (TOC p.393–471)
- Key sub-sections: §"Kinetic Monte Carlo Method" (DADOS-based [1][2][3]) line 47618; §"Operating Modes" 47644; §"Translating Atomistic/Nonatomistic Information" (`PDE2KMC`/`KMC2PDE`) line 48013; §"KMC Domain/Grid" 48153; defect configurations & binding energies TOC p.423; §"Oxidation-Enhanced Diffusion (OED) Model" (KMC) TOC p.511; extraction (deatomize) TOC p.532.

**Solver backend:** KMC-to-continuum L2 projection (Mx=b mass matrix solve) uses HyprePCG + HypreBoomerAMG. Sub-cycling (KMC dt << FEM dt) uses SUNDIALS CVODE with the same SPBCGS + Hypre AMG preconditioner from Phase 1. The mass matrix is cached when `fixedMeshes=true` (per MultiAppProjectionTransfer pattern) to avoid rebuilding the AMG hierarchy each step.

## Manual Equation References (SProcess Chapter 5)

The KMC algorithm is governed by transition-state theory (TST) rates. Every event `i → j` has an Arrhenius rate:

```
r_ij = ν_0 · exp(−(E_m + ΔE_ij) / kT)                                    (SProcess Ch.5 §"KMC Method", line 47618)
```

- `ν_0` ≈ 10¹³ s⁻¹ (attempt frequency, Debye frequency of Si)
- `E_m` = migration barrier for the defect species (I: ~0.9 eV, V: ~0.5 eV; Task 3)
- `ΔE_ij` = binding-energy contribution from local defect environment (capture/emission barriers for clusters)

**BKL rejection-free time advance (SProcess line 47618, ref [1]):**

```
dt = −ln(rand) / R_total,    R_total = Σ_all events r_ij                  (BKL / Gillespie algorithm)
```

Event `k` selected with probability `r_k / R_total` (binary-heap event tree for O(log N) selection).

**Defect configurations & binding energies (SProcess TOC p.423; DADOS [2]):**
- **{311} defect** — interstitial chain/clusters; binding energy `E_b({311}) ≈ 2.0+ eV` per I (size-dependent — Ostwald ripening: small clusters dissolve, large grow)
- **Dislocation loop** — large I aggregate; binding energy per I lower than {311} above a critical size (loop → {311} transition)
- **BIC (Boron-Interstitial Cluster)** — `B_n I_m`; binding energies from DADOS calibration
- **Amorphous pocket** — disordered I/V region above the amorphization threshold

**Capture/emission rates (SProcess §"Defect Configurations", Ch.5):**

```
k_capture(size n → n+1) = 4π · R_eff(n) · D_I · C_I                      (diffusion-limited capture)
k_emit(size n → n−1)    = k_capture · exp(−E_b(n) / kT)                  (detailed balance)
```

`R_eff(n)` = effective capture radius (size-dependent).

**PDE2KMC / KMC2PDE transfer (SProcess line 48013):**
- `PDE2KMC` — continuum concentrations → atomistic sites (Poisson sampling: `N_atoms = Poisson(C·V_cell)`, random positions)
- `KMC2PDE` — atom counts → continuum field (binning + smoothing + L2 projection; the Sano method at line 48088 is the canonical smoothing approach)

## File Structure

| File | Responsibility | Status (HEAD) |
|------|----------------|---------------|
| `kmc/KmcAtomisticEngine.hpp` | BKL rejection-free KMC on Si diamond lattice | exists (196-line toy, simple-cubic) |
| `kmc/KmcLattice.hpp` | Si diamond lattice representation | exists (stub, 61 lines) |
| `kmc/KmcEvent.hpp` | Event types: hop, recomb, cluster, dissociate | exists (stub, 44 lines) |
| `kmc/KmcParameters.hpp` | Migration barriers, binding energies from parameter DB | NEW |
| `kmc/KmcDeatomize.hpp` | KMC -> continuum field transfer (smooth, project) | NEW |
| `kmc/KmcAtomize.hpp` | Continuum -> KMC initialization (sample positions) | NEW |
| `kmc/KmcReport.hpp` | Defect activity, interaction, histogram reports | NEW |

## Tasks

### Task 1: KmcLattice - Si Diamond Lattice

**Files:** Rewrite `include/viennaps/fields/kmc/KmcLattice.hpp`

- `KmcLattice` - generates Si diamond cubic atomic positions in a box. Each site: Si, dopant, I, V, empty. Neighbor list for hop events.

**Manual reference (SProcess §"KMC Domain/Grid", line 48153):** the KMC supercell is a tensor-product grid of diamond-cubic unit cells with periodic boundary conditions (SProcess line 48447). Two interpenetrating FCC sublattices offset by `(a/4, a/4, a/4)` (the diamond basis). Lattice constant `a_Si = 0.543 nm`. Each site has 4 nearest neighbors (tetrahedral coordination).

**Gap-analysis row closed:** "Si diamond cubic lattice ❌" → ✅. Current `KmcLattice` uses simple cubic (6-neighbor) per gap analysis; this task replaces it.

- Test: verify lattice constant `a = 0.543 nm`, coordination number 4 for diamond, neighbor distances `√3·a/4 ≈ 0.235 nm`.
- Commit: `"feat: rewrite KmcLattice as Si diamond cubic (2-interpenetrating FCC sublattices; SProcess 48153) — was simple-cubic"`

### Task 2: KmcEvent - Event Types

**Files:** Rewrite `include/viennaps/fields/kmc/KmcEvent.hpp`

- `KmcEvent` - base + subclasses: HopEvent (I/V/dopant migration), RecombEvent (I+V->0), ClusterEvent (I+I->{311}), DissociateEvent ({311}->I+{311}), PairEvent (B+I->BIC)
- Each event has: site index, target site, rate `r_ij = ν_0·exp(−(E_m + ΔE_ij)/kT)` (manual eq. above)

**Manual reference (SProcess Ch.5 event catalog; DADOS [2]):** every event type in SProcess KMC is one of these five classes. The `ΔE_ij` term distinguishes them — for HopEvent `ΔE_ij = 0` (free migration); for ClusterEvent `ΔE_ij < 0` (binding lowers the barrier); for DissociateEvent `ΔE_ij > 0` (must overcome binding).

**Gap-analysis rows closed:** Recombination ❌→✅, Clustering ❌→✅, Dissociation ❌→✅, Dopant-defect pairing ❌→✅, Impurity clustering ❌→✅.

- Test: verify rate computation at given T matches `ν_0·exp(−E_m/kT)` for free hops; verify cluster-event rate is faster than free hop (binding lowers barrier).
- Commit: `"feat: rewrite KmcEvent with 5 event classes (Hop/Recomb/Cluster/Dissociate/Pair; SProcess Ch.5 event catalog)"`

### Task 3: KmcParameters

**Files:** Create `include/viennaps/fields/kmc/KmcParameters.hpp`

- `KmcParameters<NumericType>` - migration barriers (Em) for I, V, B, P, As in Si. Binding energies for {311}, BIC, VC. Recombination radii. Amorphization threshold. From parameter DB.

**Manual reference (SProcess Appendix / Advanced Calibration; DADOS [2] cited throughout Ch.5):**

| Parameter | Value (Si, 300K) | Source |
|---|---|---|
| `E_m(I)` | 0.9 eV | DADOS [2] |
| `E_m(V)` | 0.5 eV | DADOS [2] |
| `E_m(B)` | 0.6 eV (B−) | DADOS |
| `ν_0` | 1e13 s⁻¹ | Debye freq |
| `E_b({311}, n=2)` | ~1.0 eV | DADOS |
| `E_b({311}, n→∞)` | ~2.4 eV (saturation) | Ostwald ripening |
| `R_recomb(I+V)` | 0.2 nm | capture radius |
| `C_amorph` | 1e22 cm⁻³ | amorphization threshold |

- Test: verify `Em(I) ~ 0.9 eV`, `Em(V) ~ 0.5 eV` (literature values).
- Commit: `"feat: add KmcParameters with DADOS-calibrated migration barriers and binding energies"`

### Task 4: KmcAtomisticEngine - BKL Algorithm

**Files:** Rewrite `include/viennaps/fields/kmc/KmcAtomisticEngine.hpp`

- `KmcAtomisticEngine<NumericType>` - BKL rejection-free KMC. Event tree (binary heap) for O(log N) selection. Time advance: `dt = -ln(rand)/total_rate` (manual eq. above). Runs for specified time or event count.
- Algorithm: build event list → select event (proportional to rate) → execute → update affected events → advance time → repeat.

**Gap-analysis rows closed:** "Event tree O(log N) ❌" → ✅ (binary heap replaces the linear scan); BKL rejection-free KMC ✅ (already works, deepened).

- Test: inject I+V pair, verify recombination occurs. Inject excess I, verify clustering into {311}.
- Commit: `"feat: rewrite KmcAtomisticEngine with binary-heap event tree for O(log N) BKL (was O(N) linear scan)"`

### Task 5: KmcDeatomize - KMC to Continuum Transfer

**Files:** Create `include/viennaps/fields/kmc/KmcDeatomize.hpp`

- `KmcDeatomize<NumericType>` - count atoms/defects in spatial bins, smooth, project onto FEM mesh as GridFunction concentrations. Smoothing radius parameter. Implements the SProcess `KMC2PDE` transfer (SProcess line 48082) using the Sano smoothing method (line 48088).

**Concrete two-step algorithm (MOOSE references verified):** the previous description ("count, smooth with Gaussian, project") hand-waves the projection. Spell it out as two composable MOOSE-proven steps:

1. **Smooth scattered atom positions onto a regular grid** using inverse-distance interpolation — exactly `framework/include/transfers/MultiAppGeometricInterpolationTransfer.h`. It already exposes the knobs you want: `_num_points`, `_power`, `_radius`, `_shrink_gap_width`, `_exclude_gap_blocks`. The Gaussian kernel the plan mentions is a special case (set `_power=2`); inverse-distance is more general and robust to clustered atoms. Implement the same kernel in MFEM-side code.
2. **Project the smoothed grid onto the FEM mesh as a GridFunction** via L2 projection — exactly `framework/src/transfers/MultiAppProjectionTransfer.C::assembleL2()` (verified lines 125-186):
   ```
   Build mass matrix     M_ij = Σ_qp JxW * φ_i(x_qp) * φ_j(x_qp)
   Build RHS             b_i  = Σ_qp JxW * smoothed_value(x_qp) * φ_i(x_qp)
   Solve                 M x = b   (tolerance 1e-10, per .C lines 504-508)
   Copy x → target GridFunction DoFs
   ```
   This is provably the L2-optimal projection and conserves the smoothed-field integral.

**Caching (important):** `MultiAppProjectionTransfer` caches quadrature points when meshes don't move (`_fixed_meshes`, `_qps_cached`, `_cached_qps` — `.C` lines 56, 293-299). KMC bins and the FEM mesh are typically stationary between sub-cycles, so caching avoids rebuilding `M` each step. **Add** a `setFixedMeshes(bool)` API and skip `M` reassembly when set.

- Test: (a) place uniform I distribution, deatomize, verify concentration field ~ uniform within 1%; (b) place a delta-like cluster, deatomize, verify the projected peak has the right integral under the smoothing kernel; (c) call `deatomize` twice in a row with `fixedMeshes=true`, verify second call is ≥10× faster than first (cache hit); (d) **dose conservation:** total projected concentration × cell volume equals original atom count within 0.1%.
- Commit: `"feat: add KmcDeatomize with IDW smoothing + L2 projection (SProcess KMC2PDE/Sano 48082-48088; MOOSE MultiAppGeometricInterpolationTransfer + MultiAppProjectionTransfer pattern) and fixed_meshes caching"`

### Task 6: KmcAtomize - Continuum to KMC Transfer

**Files:** Create `include/viennaps/fields/kmc/KmcAtomize.hpp`

- `KmcAtomize<NumericType>` - sample atomistic positions from continuum concentrations. Implements the SProcess `PDE2KMC` transfer (SProcess line 48020). Poisson placement: `N = Poisson(C·V_cell)` (manual eq. above), random positions within cell. Deterministic mode for reproducibility.

**MOOSE reference (verified):** `framework/include/transfers/MultiAppGeneralFieldNearestLocationTransfer.h` — nearest-neighbor sampling from source points to target points. The reverse direction (continuum → atoms) is conceptually identical: for each candidate atom site, sample the continuum field at that site. There is no MOOSE analog for the Poisson stochastic step itself; keep it as ViennaPS code but cite the nearest-location transfer for the field-sampling half.

**Caching:** the FEM field evaluation at candidate sites can reuse the same KD-tree as Task 5 if `fixedMeshes=true`.

- Test: uniform C=1e18, atomize, verify atom count ~ C·volume within Poisson statistical error (√N / N < 1% for N>1e4). Deterministic mode produces identical atom placements across runs.
- Commit: `"feat: add KmcAtomize with Poisson placement (SProcess PDE2KMC 48020; MOOSE MultiAppGeneralFieldNearestLocationTransfer field-sampling pattern)"`

### Task 7: KmcReport - Output Reports

**Files:** Create `include/viennaps/fields/kmc/KmcReport.hpp`

- `KmcReport` - defect activity report (counts by type), interaction report (event statistics), 1D profiles (depth-binned), supersaturation (`C_I/C_I_eq`), cluster size histograms, amorphous/crystalline interface extraction. SProcess extraction TOC p.532–551.

**MOOSE references (verified):** prefer count-histogram vs volume-weighted histogram deliberately:
- `framework/include/vectorpostprocessors/HistogramVectorPostprocessor.h` — count-based histogram (3 columns per input: counts, lower edge, upper edge). Use for cluster-size distributions (size = number of atoms).
- `framework/include/vectorpostprocessors/VariableValueVolumeHistogram.h` + `ArrayVariableValueVolumeHistogram.h` — volume-weighted histogram. **Use this for concentration-weighted distributions** — physically correct because a concentration field's histogram should weight by element volume, not by element count.
- `framework/include/reporters/StatisticsReporter.h` (stochastic_tools) — mean/STD/percentiles with bootstrap confidence intervals. Use for event-statistics reports.

**Gap-analysis rows closed:** Defect activity reports ❌→✅, Interaction reports ❌→✅, 1D profiles ❌→✅, Supersaturation ❌→✅, Cluster size histograms ❌→✅.

- Test: run small KMC, generate report, verify non-zero counts AND verify volume-weighted histogram integrates to total dose.
- Commit: `"feat: add KmcReport with count and volume-weighted histograms (SProcess Ch.5 extraction; MOOSE HistogramVectorPostprocessor + VariableValueVolumeHistogram pattern)"`

### Task 8: Integration Test - KMC <-> Continuum Coupling

- Sequence: atomize continuum B profile → run KMC for 1s → deatomize back → verify profile roughly preserved (with diffusion spreading). Compare KMC dopant profile with Fermi continuum model for same conditions.

**MOOSE coupling pattern (verified — name it explicitly):** this is a two-scale sub-cycling problem, the canonical use case for `framework/include/multiapps/TransientMultiApp.h`. Even if ViennaPS implements its own loop, the *pattern* must be explicit:
- **Sub-cycling** (`_sub_cycling`, line 59): KMC dt << FEM dt; KMC sub-steps within one FEM step.
- **Time-interpolated transfers** (`_interpolate_transfers`, `_transferred_vars`, `_transferred_dofs`, lines 60, 76-79): deatomize/atomize fields are interpolated across sub-steps.
- **Catch-up** (`_catch_up`, `_max_catch_up_steps`, lines 69-70): if KMC falls behind, allow catch-up.
- **Picard support** (`solveStep(dt, target_time, auto_advance=false)`, line 35; `resetApp(global_app, time)`, line 41): for tight two-way coupling, run KMC → FEM → check convergence → `resetApp` and re-run if not converged.
- **Recoverable failure** (`MultiAppSolveFailure`, lines 95-103): a KMC event explosion should not crash the FEM solve.

- Commit: `"test: add KMC-continuum coupling integration test (TransientMultiApp sub-cycling pattern)"`

### Task 9: KMC TED Validation

- Implant damage → atomize → KMC anneal → deatomize → verify TED behavior (transient I supersaturation, {311} formation, B enhanced diffusion). Compare with CDD continuum model. This is the canonical SProcess KMC validation scenario (Ch.5 Advanced Calibration, TOC p.562).

**Gap-analysis row closed:** "KMC vs continuum validation ❌" → ✅.

- Commit: `"test: add KMC TED validation test (SProcess Ch.5 canonical validation scenario)"`

---

## Self-Review Notes

- **Spec coverage:** Phase 7 of spec Section 12 = "KMC atomistic + deatomize/atomize coupling." Every spec §5.1 row that the gap analysis marked ❌ (Si diamond lattice, event tree, recombination, clustering, dissociation, pairing, impurity clustering, defect reports, supersaturation, histograms, KMC validation) now has a manual-cited implementation path.
- **Manual citations:** SProcess Chapter 5 body line ranges (47618, 48013, 48082, 48088, 48153, 48447) plus DADOS references. The Arrhenius rate form, BKL time advance, capture/emission detailed-balance, and PDE2KMC/KMC2PDE transfer are all equation-cited.
- **MOOSE citations:** `MultiAppGeometricInterpolationTransfer` (Task 5 IDW smoothing), `MultiAppProjectionTransfer` (Task 5 L2 projection), `MultiAppGeneralFieldNearestLocationTransfer` (Task 6 field sampling), `HistogramVectorPostprocessor`/`VariableValueVolumeHistogram` (Task 7), `TransientMultiApp` (Task 8 sub-cycling).
- **Engine integration:** reuses Phase 3 CDD/clustering (continuum counterpart), Phase 1 HypreBoomerAMG (L2 projection preconditioner).
- **Gap-analysis rows closed (summary):** Si diamond lattice ❌→✅, Event tree O(log N) ❌→✅, Recombination ❌→✅, Clustering ❌→✅, Dissociation ❌→✅, Dopant-defect pairing ❌→✅, Impurity clustering ❌→✅, Amorphous pocket ❌→✅, Defect activity reports ❌→✅, Interaction reports ❌→✅, 1D profiles ❌→✅, Supersaturation ❌→✅, Cluster size histograms ❌→✅, KMC vs continuum validation ❌→✅.
