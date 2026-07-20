# Diffusion Functionality Parity Design Spec

**Date:** 2026-07-20  
**Goal:** Full parity with Silvaco ATHENA and SProcess diffusion capabilities  
**Scope:** All diffusion sub-models described in the manuals (continuum, KMC, flash/laser anneal, polysilicon, SiGe, III-V, PDE API)

## 1. Context

ViennaPS_mod is a local fork of ViennaPS (header-only C++ library for semiconductor process simulation). It currently has toy-scale 1D profile diffusion kernels (`DiffusionKernel`, `FermiDiffusionKernel`, `PairDiffusionKernel`, `ChargedReactKernel`, `DefectClusterKernel`, `SPERKernel`) that operate on 1D depth profiles with explicit time stepping. These are placeholders with simplified physics.

The Silvaco ATHENA User's Manual (Chapter 3.1-3.2) and SProcess User Guide (Chapters 4-7) describe full 2D/3D FEM-based diffusion with 10+ transport models, point-defect clustering, segregation, polysilicon grain diffusion, SiGe/SiGeC interdiffusion, III-V compound diffusion, flash/laser anneal, atomistic KMC, lattice KMC epitaxy, and Alagator PDE specification.

This spec designs a full FEM-based diffusion engine to replace the 1D placeholder kernels.

## 2. Architecture

### 2.1 Design Decision: Monolithic FEM Engine

A single MFEM-based `DiffusionEngine` handles all continuum diffusion models. Each model contributes stiffness matrix entries (diffusion terms) and nonlinear RHS contributions (reaction/clustering terms) to the global FEM system. KMC runs as a separate engine with a deatomize/atomize transfer interface.

**Rationale:** The models are coupled through shared unknowns (dopant + interstitial + vacancy + clusters in one system). A per-kernel RHS callback interface cannot cleanly express the coupled stiffness assembly. SProcess/ATHENA internally use a monolithic FEM solver.

### 2.2 System Architecture

```
DiffusionEngine
  +-- MeshModule (LevelSet -> MFEM, moving boundary)
  +-- ModelRegistry (Fermi, CDD, cluster, segregation, ...)
  +-- FEMAssembly (global system: M*du/dt = -K(u)*u + R(u))
  +-- TimeIntegrator (SUNDIALS CVODE/BDF)

KmcEngine (separate)
  +-- KmcAtomisticEngine (point-defect KMC, Chapter 5)
  +-- KmcLatticeEngine (epitaxy KMC, Chapter 6)
  +-- DeatomizeTransfer (KMC <-> continuum coupling)

ResultsExtractor (1D cuts, dose, sheet resistance, fitting)
PdeAPI (C++ function-object PDE specification)
ParameterDatabase (calibrated values from manuals)
```

### 2.3 Key Classes

```cpp
DiffusionEngine<NumericType, D>
  - owns: MFEM::Mesh, GridFunction vector u, ModelRegistry
  - assemble(): builds M, K(u), R(u) from all registered models
  - solve(t_start, t_end, dt_max): SUNDIALS integration
  - setMeshFromLevelSet(domain): generates FEM mesh from level-set
  - updateMesh(domain): moving-boundary mesh deformation/remesh

DiffusionModel (abstract base)
  - assembleStiffness(BilinearForm& K, const Vector& u, attrs)
  - assembleReaction(Vector& R, const Vector& u, attrs)
  - assembleMass(BilinearForm& M, attrs)
  - numSpecies(), speciesNames()

KmcEngine<NumericType>
  - runAtomistic(params): atomistic KMC simulation
  - runLattice(params): lattice KMC epitaxy
  - deatomize(field): KMC -> continuum field transfer
  - atomize(field): continuum -> KMC initialization
```

## 3. Mesh Generation (LevelSet -> MFEM)

### 3.1 LevelSetToMeshConverter

- **Boundary extraction**: Marching squares (2D) / marching cubes (3D) on the level-set grid. Each level-set layer boundary becomes mesh edges/faces.
- **Triangulation**: Delaunay (2D) / Delaunay tetrahedralization (3D) of interior, constrained by extracted boundaries.
- **Material attribution**: Each mesh element inherits material ID from the level-set grid cell at its centroid. Stored as MFEM element attribute.
- **Interface-aligned refinement**: Automatic mesh densification near material boundaries.
- **User refinement boxes**: `RefinementBox(xMin, yMin, xMax, yMax, maxElementSize)` equivalent to SProcess `refinebox`.

### 3.2 MovingMeshHandler

- **Mesh deformation**: Nodes move with boundary displacement. Interior nodes smoothed via Laplacian smoothing. Matches SProcess `MovingMesh` algorithm.
- **Remesh trigger**: When max element skewness exceeds 0.8 or min Jacobian < 0, regenerate mesh.
- **Grid spacing control**: Min/max element size, max growth rate, boundary layer spacing.
- **3D moving boundary**: `MovingMesh3D` algorithm with node relaxation and quality-based remeshing.

### 3.3 SolutionTransfer

- **Mesh regeneration**: GridFunctions projected onto new mesh via MFEM L2 projection.
- **Mesh deformation**: GridFunction stays on same mesh; nodes moved. Values preserved at DoF locations.
- **Data interpolation**: Near boundaries, special interpolation avoids artificial dose loss/gain.

### 3.4 Adaptive Mesh Refinement

- **Static**: User-specified refinement boxes before diffusion step.
- **Adaptive**: During diffusion, refine/coarsen based on relative difference, gradient, local dose error criteria. Uses MFEM `Mesh::GeneralRefinement` and `Mesh::Derefine`.

## 4. Continuum Diffusion Models

### 4.1 Unknown Vector Layout

The global unknown vector `u` packs all species: `[C_Boron, C_Phosphorus, ..., C_I, C_V, C_311, C_VC, C_BIC, C_loop, ...]`. Each species occupies a block of DoFs. System size is `nNodes * nSpecies`.

### 4.2 Transport Models

| Model | Formula | Manual Ref |
|-------|---------|------------|
| Constant | D = D0 * exp(-Ea/kT) | SProcess 4.191 |
| Fermi | D = D_i*(1 + alpha*n/ni) + D_v*(1 + beta*p/ni) | SProcess 4.190, ATHENA 3.1.2 |
| ChargedFermi | D = sum_z D^z * f^z(n, p, T) with explicit charge states | SProcess 4.188 |
| Pair | D = D_pair * (C_I / C_I_eq) | SProcess 4.187 |
| ChargedPair | Pair + Fermi-level coupling | SProcess 4.185 |
| React | I+V recombination: dC_I/dt -= k*C_I*C_V | SProcess 4.183 |
| ChargedReact | React + charge-state-dependent rates | SProcess 4.176 |
| NeutralReact | Neutral defect reactions | SProcess 4.192 |
| CDD | Full coupled: dopant + I + V with pairing, recomb, clustering | SProcess 4.175, ATHENA 3.2.1 |
| ChargedEquilibrium | Equilibrium concentrations with charge | SProcess 4.194 |
| Carbon | Carbon-specific diffusion and I trapping | SProcess 4.193 |
| Nitrogen | Nitrogen diffusion model | SProcess 4.193 |
| Copper | Cu diffusion + ion-pairing | SProcess 4.195 |

### 4.3 Cluster and Deactivation Models

| Model | Unknowns | Reaction | Manual Ref |
|-------|----------|----------|------------|
| {311} clusters | C_311 | dC_311/dt = k_f*C_I^n - k_r*C_311 | ATHENA 3.2.3, SProcess 5 |
| Vacancy clusters (VC) | C_VC | dC_VC/dt = k_f*C_V^m - k_r*C_VC | ATHENA 3.2.4 |
| Impurity clustering (DDC) | C_BIC, C_cluster | Dopant + I/V -> cluster | ATHENA 3.2.5 |
| Dislocation loops | C_loop | Loop growth from I supersaturation | ATHENA 3.2, SProcess 5 |
| Solid solubility | (modifies active C) | C_active = min(C, C_ss(T)); excess -> cluster | ATHENA 3.2.2 |

### 4.4 Interface Physics

| Model | Implementation | Manual Ref |
|-------|---------------|------------|
| Segregation | BC at material interface: C_2 = m(T)*C_1 | SProcess 4.173, ATHENA 3.1.3 |
| OED | Oxidation injects I at Si/SiO2 interface: flux BC on C_I | SProcess 9, ATHENA 3.1.5 |
| TED | Initial condition: C_I(x,0) from implant damage profile | SProcess 5, ATHENA 3.2 |
| Dose loss | Surface BC: -D*dC/dn = h*C (evaporation) | SProcess 4.216 |
| Grain boundary segregation | C_gb = m*C_interior (polysilicon) | SProcess 4.232 |

### 4.5 Polysilicon Diffusion (SProcess 4.224-239)

**Isotropic model:** D_eff = D_bulk + D_gb * f_gb (grain boundary fraction). Grain growth: dR/dt = k * exp(-Ea/kT) / R^n. Segregation between grain interior and boundary.

**Anisotropic model:** Dual mesh (grain interior + grain boundary). Different diffusivities. Voronoi tessellation of grain centers. Grain size model with surface nucleation. Interface oxide breakup and epitaxial regrowth.

### 4.6 SiGe/SiGeC Diffusion (SProcess 4.241-251)

- SiGe interdiffusion: D_SiGe(x_Ge, T) with concentration-dependent terms
- Bandgap model: E_g(x_Ge, strain) affects n_i and Fermi-level-dependent D
- Boron diffusion: enhanced/suppressed by Ge via bandgap and point-defect parameters
- Carbon suppression: C traps interstitials, suppressing B transient diffusion
- Strain effects on dopant activation and point-defect equilibrium
- Ge-B pairing model with cluster initialization

### 4.7 III-V Compound Semiconductor Diffusion (SProcess 4.251-252)

- Material conversion: Convert Si mesh to GaAs/InP with new properties
- Species-specific D in III-V materials
- Different I/V equilibrium in compound semiconductors

### 4.8 SPER (SProcess 4.197)

- Interface velocity: v = v0 * exp(-Ea/kT), orientation-dependent
- Dopant release: trapped dopants released at crystallization front
- Defect emission: end-of-range defects as interstitial clusters
- Coupling: SPER front tracked as level-set; FEM mesh updated as front advances

### 4.9 Flash/Laser Anneal (SProcess 4.201-223)

**Thermal module:** Heat transfer equation (rho*c_p*dT/dt = div(k*grad(T)) + Q). Energy implantation, transfer matrix method (optical absorption), FDTD for sub-wavelength features.

**Phase field:** Melting phase field (liquid/solid), crystallinity phase field (amorphous/crystalline). Coupled with SPER.

**Dopant diffusion in melt:** Liquid-phase D (orders of magnitude faster). Solute transport in liquid Si. Resolidification trapping.

**Intensity models:** Gaussian (flash), table lookup, user-specified (PDE API), scanning laser with control parameters.

### 4.10 FEM Assembly Pattern

```cpp
class DiffusionModel {
public:
  virtual void assembleStiffness(mfem::BilinearForm& K,
                                  const mfem::Vector& u,
                                  const MeshAttributes& attrs) const = 0;
  virtual void assembleReaction(mfem::Vector& R,
                                const mfem::Vector& u,
                                const MeshAttributes& attrs) const = 0;
  virtual void assembleMass(mfem::BilinearForm& M,
                            const MeshAttributes& attrs) const = 0;
  virtual int numSpecies() const = 0;
  virtual std::vector<std::string> speciesNames() const = 0;
  virtual std::vector<int> applicableAttributes() const = 0;
};
```

## 5. KMC Engine

### 5.1 Atomistic KMC (SProcess Chapter 5)

**KmcAtomisticEngine** - point-defect diffusion via kinetic Monte Carlo.

**Lattice:** Si diamond cubic lattice. Each site holds Si atom, dopant, I, V, or empty.

**Events:** Hop (I/V/dopant migration, rate = nu0*exp(-Em/kT)), recombination (I+V->0), clustering (I+I->{311}), dissociation, amorphous pocket formation, dopant-defect pairing (B+I->BIC), impurity clustering.

**Algorithm:** BKL rejection-free KMC with event tree for O(log N) selection.

**Coupling:** `deatomize(field)` converts atomistic counts to continuum concentrations (smoothed, projected onto FEM mesh). `atomize(field)` samples atomistic positions from continuum concentrations.

**Output:** Defect activity reports, interaction reports, 1D profiles, supersaturation, amorphous/crystalline interfaces, cluster size histograms.

### 5.2 Lattice KMC Epitaxy (SProcess Chapter 6)

**KmcLatticeEngine** - epitaxial deposition via lattice KMC.

**Models:** Planar epitaxy (growth rate), coordination-based (atomic bonding), coordination-reactions (detailed surface chemistry with explicit reaction rates).

**Physics:** SiGe mole fraction-dependent growth, visibility/shadowing, twin-defect formation on {111}, surface segregation of dopants, secondary reactions, dangling bonds via shared occupancy, nonselective epitaxial deposition of polysilicon on oxide, dopant activation/clustering during epitaxy.

**Coupling:** Updates FEM mesh during/after epitaxy. Nonatomistic mode for continuum coupling.

## 6. PDE API (Alagator Equivalent)

C++ function-object based PDE specification:

```cpp
auto dopantEqn = PdeEquation("Dopant")
  .diffusion([D](auto T, auto C) { return D; })
  .reaction([](auto C) { return -k * C * C; })
  .boundary("surface", PdeBC::evaporation(h))
  .boundary("interface", PdeBC::segregation(m, "Oxide"));
engine.addEquation(dopantEqn);
```

- **PdeTerm**: base class for diffusion, reaction, growth, boundary condition terms
- **PdeEquation**: collection of terms for one species (one block in global system)
- **PdeBC**: Dirichlet, Neumann, segregation, flux boundary conditions
- **Built-in terms**: DiffusionTerm, ReactionTerm, GrowthTerm, SegregationBC, FluxBC, DirichletBC
- **Modifying built-ins**: `engine.addEquationTerm("Dopant", extraTerm)`, `engine.subEquationTerm("Dopant", term)`

## 7. Parameter Database

Extend existing `ParameterDatabase` with calibrated values from manuals:

- **Species parameters**: D0, Ea for each dopant (B, P, As, Sb, In, Ga) in each material (Si, SiO2, Si3N4, PolySi, SiGe)
- **Point-defect parameters**: I/V equilibrium concentrations, diffusivities, formation energies
- **Cluster parameters**: {311} binding energy, loop growth rate, BIC binding
- **Segregation coefficients**: m(T) for each dopant/material-interface pair
- **Oxidation parameters**: linear/parabolic rate constants, orientation factors
- **Inheritance**: SiGe inherits from Si with Ge-dependent modifications
- **Like-materials blend**: Si_{1-x}Ge_x interpolates between Si and Ge parameters

## 8. Results Extraction

- **1D data cuts**: `extract1D(field, startPoint, endPoint, nPoints)` - linear interpolation through FEM solution
- **Dose calculation**: `integrateField(field, region)` - area integral of concentration
- **Level crossings**: `findLevel(field, threshold, direction)` - depth where C = threshold
- **Sheet resistance**: `calcSheetResistance(field, mobility)` - integral of 1/(q*mu*n) dx
- **Fitting**: FitArrhenius, FitLine, FitPearson, FitPearsonFloor - least-squares utilities
- **Output**: VTK (existing), TDR-equivalent (HDF5 via MFEM), CSV

## 9. Testing and Validation

- **Unit tests**: Each model tested against analytical solutions (constant-D diffusion, recomb decay, segregation equilibrium)
- **Integration tests**: Multi-step sequences (implant -> diffuse -> oxidize) compared against manual examples
- **Regression tests**: Snapshot of 1D profiles for known conditions (e.g., 50keV B implant, 1000C 30s anneal)
- **KMC validation**: Compare KMC dopant profile with continuum Fermi model for same conditions
- **Mesh convergence**: Verify solution converges as mesh refines
- **Test registration**: New `tests/diffusion/` directory with CTest entries

## 10. File Structure

New files under `include/viennaps/`:

```
fields/
  DiffusionEngine.hpp          # Main FEM diffusion engine
  DiffusionModel.hpp           # Abstract base for all models
  LevelSetToMesh.hpp           # Level-set -> MFEM mesh converter
  MovingMeshHandler.hpp        # Moving boundary mesh management
  SolutionTransfer.hpp         # Field transfer between meshes
  PdeEquation.hpp              # PDE API (C++ function objects)
  PdeTerm.hpp                  # PDE term base + built-in terms
  PdeBC.hpp                    # Boundary condition types
  KmcAtomisticEngine.hpp       # Atomistic KMC
  KmcLatticeEngine.hpp         # Lattice KMC epitaxy
  KmcDeatomize.hpp             # KMC <-> continuum transfer
  ResultsExtractor.hpp         # 1D cuts, dose, fitting
  models/
    ConstantDiffusion.hpp
    FermiDiffusion.hpp
    ChargedFermiDiffusion.hpp
    PairDiffusion.hpp
    ChargedPairDiffusion.hpp
    ReactDiffusion.hpp
    ChargedReactDiffusion.hpp
    NeutralReactDiffusion.hpp
    CddDiffusion.hpp
    ChargedEquilibriumDiffusion.hpp
    Cluster311.hpp
    VacancyCluster.hpp
    ImpurityCluster.hpp
    DislocationLoop.hpp
    SolidSolubility.hpp
    Segregation.hpp
    OedSource.hpp
    PolysiliconDiffusion.hpp
    SiGeDiffusion.hpp
    FlashLaserAnneal.hpp
    HeatTransfer.hpp
    MeltingPhaseField.hpp
```

## 11. Backward Compatibility

- Existing 1D profile kernels (`DiffusionKernel`, `FermiDiffusionKernel`, etc.) remain as deprecated fallbacks when MFEM is not available (`#ifndef VIENNAPS_HAS_MFEM`).
- `ProcessOrchestrator` is extended to optionally use the new `DiffusionEngine` when MFEM is present.
- Existing `PhysicsField` 1D profile API remains for backward compatibility with existing tests.

## 11.1. Environment

Prebuilt MFEM is at `f:/dev/mfem/build` (Release) and `f:/dev/mfem/build_debug` (Debug). This resolves the MSVC CRT mismatch: link Release app with Release MFEM (`build`), Debug app with Debug MFEM (`build_debug`). vcpkg installed deps at `f:/dev/vcpkg/installed/x64-windows/`.

## 12. Implementation Phasing

Although this spec covers the full scope, implementation will be phased:

| Phase | Content | Dependencies |
|-------|---------|-------------|
| 1 | Mesh generation (LevelSet->MFEM) + Constant diffusion + SUNDIALS coupling | None |
| 2 | Fermi + ChargedFermi + segregation + solid solubility | Phase 1 |
| 3 | CDD + React/ChargedReact + Pair/ChargedPair + clustering ({311}, VC, BIC, loops) | Phase 2 |
| 4 | OED + TED + dose loss + interface physics | Phase 3 |
| 5 | Polysilicon (isotropic + anisotropic) | Phase 3 |
| 6 | SiGe/SiGeC + III-V | Phase 3 |
| 7 | KMC atomistic + deatomize/atomize coupling | Phase 3 |
| 8 | KMC lattice epitaxy | Phase 7 |
| 9 | Flash/laser anneal (heat transfer + phase field + melt diffusion) | Phase 3 |
| 10 | PDE API + results extraction + parameter database calibration | Phase 2 |
| 11 | Adaptive mesh refinement | Phase 1 |
