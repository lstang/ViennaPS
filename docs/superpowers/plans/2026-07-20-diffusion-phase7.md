# Diffusion Engine Phase 7: KMC Atomistic + Deatomize/Atomize Coupling

**Goal:** Add atomistic Kinetic Monte Carlo (KMC) for point-defect diffusion, with deatomize/atomize transfer to couple with the continuum FEM engine.

**Depends on:** Phase 3 (CDD, clustering)

## File Structure

| File | Responsibility |
|------|---------------|
| `KmcAtomisticEngine.hpp` | BKL rejection-free KMC on Si diamond lattice |
| `KmcLattice.hpp` | Si diamond lattice representation |
| `KmcEvent.hpp` | Event types: hop, recomb, cluster, dissociate |
| `KmcParameters.hpp` | Migration barriers, binding energies from parameter DB |
| `KmcDeatomize.hpp` | KMC -> continuum field transfer (smooth, project) |
| `KmcAtomize.hpp` | Continuum -> KMC initialization (sample positions) |
| `KmcReport.hpp` | Defect activity, interaction, histogram reports |

## Tasks

### Task 1: KmcLattice - Si Diamond Lattice
- `KmcLattice` - generates Si diamond cubic atomic positions in a box. Each site: Si, dopant, I, V, empty. Neighbor list for hop events
- Test: verify lattice constant, coordination number (4 for diamond), neighbor distances
- Commit: `"feat: add KmcLattice Si diamond lattice representation"`

### Task 2: KmcEvent - Event Types
- `KmcEvent` - base + subclasses: HopEvent (I/V/dopant migration), RecombEvent (I+V->0), ClusterEvent (I+I->{311}), DissociateEvent ({311}->I+{311}), PairEvent (B+I->BIC)
- Each event has: site index, target site, rate = nu0*exp(-Em/kT)
- Test: verify rate computation at given T
- Commit: `"feat: add KmcEvent types with Arrhenius rates"`

### Task 3: KmcParameters
- `KmcParameters<NumericType>` - migration barriers (Em) for I, V, B, P, As in Si. Binding energies for {311}, BIC, VC. Recombination radii. Amorphization threshold. From parameter DB
- Test: verify Em(I) ~ 0.9 eV, Em(V) ~ 0.5 eV (literature values)
- Commit: `"feat: add KmcParameters with migration barriers"`

### Task 4: KmcAtomisticEngine - BKL Algorithm
- `KmcAtomisticEngine<NumericType>` - BKL rejection-free KMC. Event tree (binary heap) for O(log N) selection. Time advance: dt = -ln(rand)/total_rate. Runs for specified time or event count
- Algorithm: build event list -> select event (proportional to rate) -> execute -> update affected events -> advance time -> repeat
- Test: inject I+V pair, verify recombination occurs. Inject excess I, verify clustering
- Commit: `"feat: add KmcAtomisticEngine with BKL algorithm"`

### Task 5: KmcDeatomize - KMC to Continuum Transfer

- `KmcDeatomize<NumericType>` - count atoms/defects in spatial bins, smooth, project onto FEM mesh as GridFunction concentrations. Smoothing radius parameter.

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
- Commit: `"feat: add KmcDeatomize with IDW smoothing + L2 projection (MOOSE MultiAppGeometricInterpolationTransfer + MultiAppProjectionTransfer pattern) and fixed_meshes caching"`

### Task 6: KmcAtomize - Continuum to KMC Transfer

- `KmcAtomize<NumericType>` - sample atomistic positions from continuum concentrations. Poisson placement: N = Poisson(C*V_cell), random positions within cell. Deterministic mode for reproducibility.

**MOOSE reference (verified):** `framework/include/transfers/MultiAppGeneralFieldNearestLocationTransfer.h` — nearest-neighbor sampling from source points to target points. The reverse direction (continuum → atoms) is conceptually identical: for each candidate atom site, sample the continuum field at that site. There is no MOOSE analog for the Poisson stochastic step itself; keep it as ViennaPS code but cite the nearest-location transfer for the field-sampling half.

**Caching:** the FEM field evaluation at candidate sites can reuse the same KD-tree as Task 5 if `fixedMeshes=true`.

- Test: uniform C=1e18, atomize, verify atom count ~ C*volume within Poisson statistical error (√N / N < 1% for N>1e4). Deterministic mode produces identical atom placements across runs.
- Commit: `"feat: add KmcAtomize with Poisson placement (MOOSE MultiAppGeneralFieldNearestLocationTransfer field-sampling pattern)"`

### Task 7: KmcReport - Output Reports

- `KmcReport` - defect activity report (counts by type), interaction report (event statistics), 1D profiles (depth-binned), supersaturation (C_I/C_I_eq), cluster size histograms, amorphous/crystalline interface extraction.

**MOOSE references (verified):** prefer count-histogram vs volume-weighted histogram deliberately:
- `framework/include/vectorpostprocessors/HistogramVectorPostprocessor.h` — count-based histogram (3 columns per input: counts, lower edge, upper edge). Use for cluster-size distributions (size = number of atoms).
- `framework/include/vectorpostprocessors/VariableValueVolumeHistogram.h` + `ArrayVariableValueVolumeHistogram.h` — volume-weighted histogram. **Use this for concentration-weighted distributions** — physically correct because a concentration field's histogram should weight by element volume, not by element count.
- `framework/include/reporters/StatisticsReporter.h` (stochastic_tools) — mean/STD/percentiles with bootstrap confidence intervals. Use for event-statistics reports.

- Test: run small KMC, generate report, verify non-zero counts AND verify volume-weighted histogram integrates to total dose.
- Commit: `"feat: add KmcReport with count and volume-weighted histograms (MOOSE HistogramVectorPostprocessor + VariableValueVolumeHistogram pattern)"`

### Task 8: Integration Test - KMC <-> Continuum Coupling

- Sequence: atomize continuum B profile -> run KMC for 1s -> deatomize back -> verify profile roughly preserved (with diffusion spreading). Compare KMC dopant profile with Fermi continuum model for same conditions.

**MOOSE coupling pattern (verified — name it explicitly):** this is a two-scale sub-cycling problem, the canonical use case for `framework/include/multiapps/TransientMultiApp.h`. Even if ViennaPS implements its own loop, the *pattern* must be explicit:
- **Sub-cycling** (`_sub_cycling`, line 59): KMC dt << FEM dt; KMC sub-steps within one FEM step.
- **Time-interpolated transfers** (`_interpolate_transfers`, `_transferred_vars`, `_transferred_dofs`, lines 60, 76-79): deatomize/atomize fields are interpolated across sub-steps.
- **Catch-up** (`_catch_up`, `_max_catch_up_steps`, lines 69-70): if KMC falls behind, allow catch-up.
- **Picard support** (`solveStep(dt, target_time, auto_advance=false)`, line 35; `resetApp(global_app, time)`, line 41): for tight two-way coupling, run KMC → FEM → check convergence → `resetApp` and re-run if not converged.
- **Recoverable failure** (`MultiAppSolveFailure`, lines 95-103): a KMC event explosion should not crash the FEM solve.

- Commit: `"test: add KMC-continuum coupling integration test (TransientMultiApp sub-cycling pattern)"`

### Task 9: KMC TED Validation
- Implant damage -> atomize -> KMC anneal -> deatomize -> verify TED behavior (transient I supersaturation, {311} formation, B enhanced diffusion). Compare with CDD continuum model
- Commit: `"test: add KMC TED validation test"`
