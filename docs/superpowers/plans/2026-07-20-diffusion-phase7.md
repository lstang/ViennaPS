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
- `KmcDeatomize<NumericType>` - count atoms/defects in spatial bins, smooth (Gaussian kernel), project onto FEM mesh as GridFunction concentrations. Smoothing radius parameter
- Test: place uniform I distribution, deatomize, verify concentration field ~ uniform
- Commit: `"feat: add KmcDeatomize transfer from atomistic to continuum"`

### Task 6: KmcAtomize - Continuum to KMC Transfer
- `KmcAtomize<NumericType>` - sample atomistic positions from continuum concentrations. Poisson placement: N = Poisson(C*V_cell), random positions within cell. Deterministic mode for reproducibility
- Test: uniform C=1e18, atomize, verify atom count ~ C*volume
- Commit: `"feat: add KmcAtomize transfer from continuum to atomistic"`

### Task 7: KmcReport - Output Reports
- `KmcReport` - defect activity report (counts by type), interaction report (event statistics), 1D profiles (depth-binned), supersaturation (C_I/C_I_eq), cluster size histograms, amorphous/crystalline interface extraction
- Test: run small KMC, generate report, verify non-zero counts
- Commit: `"feat: add KmcReport for defect activity and profiles"`

### Task 8: Integration Test - KMC <-> Continuum Coupling
- Sequence: atomize continuum B profile -> run KMC for 1s -> deatomize back -> verify profile roughly preserved (with diffusion spreading). Compare KMC dopant profile with Fermi continuum model for same conditions
- Commit: `"test: add KMC-continuum coupling integration test"`

### Task 9: KMC TED Validation
- Implant damage -> atomize -> KMC anneal -> deatomize -> verify TED behavior (transient I supersaturation, {311} formation, B enhanced diffusion). Compare with CDD continuum model
- Commit: `"test: add KMC TED validation test"`
