# Diffusion Engine Phase 9: Flash/Laser Anneal

**Goal:** Add flash and laser anneal simulation: heat transfer, melting phase field, crystallinity phase field, dopant diffusion in melt, and intensity models. Phase fields use Allen-Cahn equation (MOOSE `ADAllenCahn` + `ACInterface` pattern) for non-conserved order parameters.

**Depends on:** Phase 3 (CDD)

## File Structure

| File | Responsibility |
|------|---------------|
| `models/AllenCahnTerm.hpp` | Allen-Cahn kernel: dη/dt = L*(κ∇²η - df/dη) (MOOSE ADAllenCahn pattern) |
| `models/ACInterfaceTerm.hpp` | Gradient energy term: κ∇η·∇test (MOOSE ACInterface pattern) |
| `models/HeatTransfer.hpp` | Thermal FEM: rho*cp*dT/dt = div(k*grad(T)) + Q (MOOSE ADHeatConduction pattern) |
| `models/LaserIntensity.hpp` | Gaussian, table lookup, scanning laser intensity |
| `models/TransferMatrix.hpp` | Optical absorption in multilayer (thin film optics) |
| `models/MeltingPhaseField.hpp` | Liquid/solid Allen-Cahn: f(η,T) = (η²-1)²/4 - λ(T-T_m)*η |
| `models/CrystallinityPhaseField.hpp` | Amorphous/crystalline Allen-Cahn: f(η,T) driven by SPER velocity |
| `models/MeltDiffusion.hpp` | Liquid-phase dopant diffusion (very fast D, activated where η_melt>0.5) |
| `models/FlashLaserAnneal.hpp` | Orchestrator: thermal + phase + diffusion |

## Tasks

### Task 1: HeatTransfer Model
- `HeatTransfer<NumericType>` - FEM thermal solver. rho*c_p*dT/dt = div(k*grad(T)) + Q(x,t). Temperature-dependent k, c_p. Uses same MFEM mesh as diffusion. Adds T as a species in the unknown vector
- Test: uniform Q -> T rises uniformly. Q=0, boundary T=1000 -> steady-state gradient
- Commit: `"feat: add HeatTransfer FEM thermal solver"`

### Task 2: TransferMatrix - Optical Absorption
- `TransferMatrix<NumericType>` - thin film optics for multilayer stack. Computes absorbed energy Q(x) from laser intensity, layer refractive indices, absorption coefficients. Fresnel equations at each interface
- Test: single Si layer, known absorption coefficient -> verify exponential decay of Q with depth
- Commit: `"feat: add TransferMatrix optical absorption model"`

### Task 3: LaserIntensity Models
- `LaserIntensity<NumericType>` - Gaussian (flash): I(r,t) = I0*exp(-r^2/2sigma^2)*pulse(t). Scanning laser: I(x,t) = I0*exp(-(x-v*t)^2/2sigma^2). Table lookup: user-specified I(x,t). User-specified: via PDE API callback
- Test: Gaussian peak at center, scanning laser moves with v*t
- Commit: `"feat: add LaserIntensity models (Gaussian, scanning, table)"`

### Task 3.5: AllenCahnTerm + ACInterfaceTerm (MOOSE Phase Field Pattern)
- Create `models/AllenCahnTerm.hpp` + `models/ACInterfaceTerm.hpp`
- **MOOSE reference:** `ADAllenCahn::computeDFDOP()` returns `_dFdEta[_qp]` (bulk driving force). `ACInterface::computeQpResidual()` = `_grad_u * kappaNablaLPsi()` (gradient energy).
- `AllenCahnTerm`: bulk driving force -L*df/dη. `ACInterfaceTerm`: gradient energy L*κ*∇η·∇test. Combined: dη/dt = L*(κ∇²η - df/dη)
- Test: double-well free energy f=(η²-1)²/4. Starting η=0.1 -> evolves to η=1. Starting η=-0.1 -> evolves to η=-1
- Commit: `"feat: add AllenCahnTerm + ACInterfaceTerm (MOOSE ADAllenCahn + ACInterface pattern)"`

### Task 4: MeltingPhaseField (Allen-Cahn Equation)
- `MeltingPhaseField<NumericType>` - Allen-Cahn equation for liquid/solid order parameter η_m in [-1,1] (-1=solid, 1=liquid)
- Free energy: f(η,T) = (η²-1)²/4 - λ(T-T_melt)*η (double-well + thermal driving force)
- df/dη = η³ - η - λ(T-T_melt). When T>T_melt, driving force pushes η->1 (liquid)
- Latent heat release couples back to heat equation via L*∂η/∂t term in HeatTransfer
- Composed using AllenCahnTerm + ACInterfaceTerm + TimeDerivativeTerm
- Test: T > T_melt -> η_m -> 1 (liquid). T < T_melt -> η_m -> -1 (solid). Sharp interface at T_melt
- Commit: `"feat: add MeltingPhaseField using Allen-Cahn equation"`

### Task 5: CrystallinityPhaseField (Allen-Cahn Equation)
- `CrystallinityPhaseField<NumericType>` - Allen-Cahn for amorphous/crystalline η_c in [-1,1] (-1=amorphous, 1=crystalline)
- Free energy: f(η,T) = (η²-1)²/4 - λ*v_SPER(T)*η. Driving force from SPER velocity
- df/dη = η³ - η - λ*v_SPER(T). SPER velocity v = v0*exp(-Ea/kT)
- During melt (η_m>0): η_c resets to 0 (liquid has no crystallinity) via coupling term
- Coupled with existing SPERKernel for interface velocity
- Test: amorphous layer -> anneal -> η_c advances from crystalline seed. Melt -> η_c resets
- Commit: `"feat: add CrystallinityPhaseField using Allen-Cahn equation"`

### Task 6: MeltDiffusion
- `MeltDiffusion<NumericType>` - liquid-phase dopant diffusion. D_liquid >> D_solid (orders of magnitude). Activated where η_m > 0 (liquid phase field from Allen-Cahn). D_eff = D_solid + (D_liquid - D_solid) * (1+η_m)/2. Solute transport in liquid Si. Resolidification traps dopant at solidification front
- Test: melt region -> dopant diffuses much faster than solid. Resolidify -> trap
- Commit: `"feat: add MeltDiffusion for liquid-phase dopant transport"`

### Task 7: FlashLaserAnneal Orchestrator
- `FlashLaserAnneal<NumericType>` - ties together: laser intensity -> Q(x,t) -> heat transfer -> T(x,t) -> melting phase field -> crystallinity phase field -> dopant diffusion (solid + melt). Multi-physics coupled solve
- Test: flash pulse -> T rises -> melts -> dopant redistributes -> cools -> resolidifies
- Commit: `"feat: add FlashLaserAnneal multi-physics orchestrator"`

### Task 8: FDTD for Sub-Wavelength Features
- Optional: FDTD simulation of electromagnetic field for sub-wavelength structures. Computes Q(x,y,z) from Maxwell's equations. Used when feature size ~ wavelength
- Test: simple 1D layered -> verify matches transfer matrix
- Commit: `"feat: add FDTD for sub-wavelength laser absorption"`

### Task 9: Saving Thermal Profile
- Save T(x,t) profile for use in subsequent diffusion steps. Load saved profile for non-thermal-coupled diffusion
- Test: save -> load -> verify identical
- Commit: `"feat: add thermal profile save/load"`

### Task 10: Integration Test - Flash Anneal
- Implant B -> flash anneal (melt) -> verify: (1) T reaches T_melt, (2) dopant redistributes in melt, (3) resolidification traps dopant, (4) final profile is shallower/broader than solid-state anneal
- Commit: `"test: add flash anneal integration test"`

### Task 11: Integration Test - Laser Scan
- Scanning laser across surface -> verify T profile follows laser, melting zone moves with laser, dopant redistribution follows melt zone
- Commit: `"test: add scanning laser anneal integration test"`
