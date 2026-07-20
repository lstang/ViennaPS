# Diffusion Engine Phase 9: Flash/Laser Anneal

**Goal:** Add flash and laser anneal simulation: heat transfer, melting phase field, crystallinity phase field, dopant diffusion in melt, and intensity models.

**Depends on:** Phase 3 (CDD)

## File Structure

| File | Responsibility |
|------|---------------|
| `models/HeatTransfer.hpp` | Thermal FEM: rho*cp*dT/dt = div(k*grad(T)) + Q |
| `models/LaserIntensity.hpp` | Gaussian, table lookup, scanning laser intensity |
| `models/TransferMatrix.hpp` | Optical absorption in multilayer (thin film optics) |
| `models/MeltingPhaseField.hpp` | Liquid/solid phase field tracking |
| `models/CrystallinityPhaseField.hpp` | Amorphous/crystalline tracking, SPER coupling |
| `models/MeltDiffusion.hpp` | Liquid-phase dopant diffusion (very fast D) |
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

### Task 4: MeltingPhaseField
- `MeltingPhaseField<NumericType>` - phase field phi_m in [0,1] (0=solid, 1=liquid). Evolves via: dphi_m/dt = M * (driving_force). Driving force = T - T_melt. Latent heat release couples back to heat equation
- Test: T > T_melt -> phi_m -> 1. T < T_melt -> phi_m -> 0. Sharp interface at T_melt
- Commit: `"feat: add MeltingPhaseField for liquid/solid tracking"`

### Task 5: CrystallinityPhaseField
- `CrystallinityPhaseField<NumericType>` - phase field phi_c in [0,1] (0=amorphous, 1=crystalline). SPER velocity drives crystallization. Coupled with SPERKernel. During melt: phi_c resets (liquid has no crystallinity)
- Test: amorphous layer -> anneal -> phi_c advances from crystalline seed. Melt -> phi_c = 0
- Commit: `"feat: add CrystallinityPhaseField for SPER tracking"`

### Task 6: MeltDiffusion
- `MeltDiffusion<NumericType>` - liquid-phase dopant diffusion. D_liquid >> D_solid (orders of magnitude). Activated where phi_m > 0.5. Solute transport in liquid Si. Resolidification traps dopant at solidification front
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
