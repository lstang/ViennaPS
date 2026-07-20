# Diffusion Engine Phase 5: Polysilicon Diffusion (Isotropic + Anisotropic)

**Goal:** Add polysilicon grain-based diffusion: isotropic (enhanced D in grain boundaries) and anisotropic (dual-mesh grain interior + boundary with segregation).

**Depends on:** Phase 3 (CDD)

## File Structure

| File | Responsibility |
|------|---------------|
| `models/PolysiliconDiffusion.hpp` | Isotropic + anisotropic poly diffusion |
| `GrainModel.hpp` | Grain growth equation, Voronoi tessellation |
| `GrainBoundaryMesh.hpp` | Dual-mesh for anisotropic model |

## Tasks

### Task 1: GrainModel - Grain Growth
- `GrainModel<NumericType>` - grain growth equation: dR/dt = k*exp(-Ea/kT)/R^n. Tracks average grain radius over time
- Test: verify grain grows at high T, static at low T
- Commit: `"feat: add GrainModel grain growth equation"`

### Task 2: Isotropic Polysilicon Diffusion
- `PolysiliconDiffusion<NumericType>` (isotropic mode) - D_eff = D_bulk + D_gb * f_gb(R). Enhanced diffusivity from grain boundary fraction. f_gb decreases as grains grow
- Test: verify D_poly > D_single_crystal. Verify D decreases as grains grow
- Commit: `"feat: add isotropic polysilicon diffusion"`

### Task 3: GrainBoundaryMesh - Dual Mesh
- `GrainBoundaryMesh` - generates Voronoi tessellation of grain centers, creates dual mesh (interior + boundary elements). Uses MFEM mesh with distinct attributes for grain interior vs boundary
- Test: verify dual mesh has correct topology, boundary elements form connected network
- Commit: `"feat: add GrainBoundaryMesh dual-mesh generation"`

### Task 4: Anisotropic Polysilicon Diffusion
- `PolysiliconDiffusion<NumericType>` (anisotropic mode) - different D in grain interior vs boundary. Segregation at grain/boundary interface: C_gb = m * C_interior. Uses dual mesh
- Test: verify dopant segregates to grain boundaries. Verify faster diffusion along GB network
- Commit: `"feat: add anisotropic polysilicon diffusion with grain boundary segregation"`

### Task 5: Interface Oxide Breakup
- Model for poly-Si/SiO2 interface oxide breakup and epitaxial regrowth during anneal. Oxide dissolves, poly recrystallizes as epitaxial Si
- Test: verify oxide thickness decreases, poly converts to single-crystal near interface
- Commit: `"feat: add poly/SiO2 interface oxide breakup model"`

### Task 6: Grain-Size-Dependent Oxidation Rate
- Poly oxidation rate depends on grain size (smaller grains -> faster oxidation via GB diffusion of oxidant)
- Test: verify fine-grain poly oxidizes faster than coarse-grain
- Commit: `"feat: add grain-size-dependent polysilicon oxidation"`

### Task 7: Integration Test
- Deposit poly -> implant B -> anneal -> verify B diffuses faster in poly than single-crystal Si, segregates to grain boundaries
- Commit: `"test: add polysilicon diffusion integration test"`
