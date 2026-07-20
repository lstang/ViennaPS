# Diffusion Engine Phase 8: KMC Lattice Epitaxy

**Goal:** Add lattice Kinetic Monte Carlo for epitaxial deposition with coordination-based growth models, SiGe mole fraction, visibility, twin defects, surface segregation.

**Depends on:** Phase 7 (KmcAtomisticEngine, KmcLattice)

## File Structure

| File | Responsibility |
|------|---------------|
| `KmcLatticeEngine.hpp` | Lattice KMC epitaxial deposition engine |
| `KmcSurfaceEvent.hpp` | Surface attachment, desorption, diffusion events |
| `KmcEpitaxyModels.hpp` | Planar, coordination, coordination-reactions models |
| `KmcVisibility.hpp` | Shadowing of incoming species by surface features |

## Tasks

### Task 1: KmcSurfaceEvent - Surface Chemistry
- `KmcSurfaceEvent` - Si attachment (SiH4 decomposition -> Si + 4H), H desorption, dopant incorporation, Ge attachment. Rates depend on surface temperature and gas partial pressures
- Test: verify Si attachment rate increases with T (H desorption limited)
- Commit: `"feat: add KmcSurfaceEvent for epitaxial deposition"`

### Task 2: Planar Epitaxy Model
- `KmcEpitaxyPlanar` - growth rate model: v = k0*exp(-Ea/kT)*P_SiH4^alpha. Dopant-dependent growth rate factor. Simple rate-based deposition
- Test: verify growth rate matches Arrhenius at given T
- Commit: `"feat: add planar epitaxy KMC model"`

### Task 3: Coordination-Based Growth Model
- `KmcEpitaxyCoordination` - growth rate depends on local atomic coordination (number of bonded neighbors). Captures facet-dependent growth. Under-coordinated sites grow faster
- Test: verify (111) surface grows slower than (100) (fewer dangling bonds)
- Commit: `"feat: add coordination-based epitaxy KMC model"`

### Task 4: Coordination-Reactions Model
- `KmcEpitaxyReactions` - detailed surface chemistry with explicit reaction rates. SiH4 decomposition, H desorption, Si attachment, dopant incorporation, Ge attachment. Multiple reaction pathways
- Test: verify growth rate with different gas compositions (SiH4 + GeH4)
- Commit: `"feat: add coordination-reactions epitaxy KMC model"`

### Task 5: SiGe Mole Fraction-Dependent Growth
- Growth rate depends on Ge fraction. Si_{1-x}Ge_x growth rate != Si growth rate. Ge surface segregation during growth
- Test: verify growth rate changes with GeH4/SiH4 ratio
- Commit: `"feat: add SiGe mole fraction-dependent epitaxy"`

### Task 6: KmcVisibility - Shadowing
- `KmcVisibility` - checks if a surface site is visible from the gas source direction. Sites shadowed by neighboring features have reduced arrival rate. Ray casting on lattice
- Test: trench geometry -> bottom of trench has reduced visibility vs top
- Commit: `"feat: add KmcVisibility shadowing for epitaxy"`

### Task 7: Twin-Defect Formation
- Twin defects form on {111} facets during epitaxy. Stacking fault creation probability. Affects dopant incorporation at defects
- Test: verify twin defect formation on (111) growth
- Commit: `"feat: add twin-defect formation in epitaxy KMC"`

### Task 8: Surface Segregation
- Dopant surface segregation during growth: dopant pushed to surface (snowplow effect). Segregation coefficient at growth front
- Test: growing Si with B doping -> B accumulates at surface
- Commit: `"feat: add dopant surface segregation during epitaxy"`

### Task 9: Nonselective Epitaxial Deposition
- Nucleation of polysilicon on oxide during nonselective epitaxy. Island formation on oxide, coalescence, transition to epitaxial on exposed Si
- Test: verify poly nucleation on SiO2, epitaxial on Si
- Commit: `"feat: add nonselective epitaxial deposition of polysilicon"`

### Task 10: FEM Mesh Update During/After Epitaxy
- After KMC epitaxy, update FEM mesh: add new elements for grown layer. Transfer dopant fields from KMC to FEM. Nonatomistic mode for continuum coupling
- Test: grow 100nm Si -> verify FEM mesh extended, dopant profile transferred
- Commit: `"feat: add FEM mesh update after KMC epitaxy"`

### Task 11: Integration Test - Selective Epitaxy
- Patterned Si/SiO2 substrate -> selective SiGe epitaxy on Si -> verify growth only on Si (not SiO2), Ge profile, dopant incorporation
- Commit: `"test: add selective epitaxy KMC integration test"`
