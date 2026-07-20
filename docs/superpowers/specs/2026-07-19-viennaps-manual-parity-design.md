# Design: ViennaPS Functionality Parity with ATHENA and Sentaurus Process Manuals

**Date**: 2026-07-19  
**Author**: Grok (in collaboration with user)  
**Status**: Design validated through incremental sections; ready for implementation planning.

## Executive Summary

ViennaPS is a modern, high-performance C++ library for level-set based topography simulation with strong Monte Carlo flux capabilities and an advanced thermal oxidation model (including stress, mask bending, and LOCOS). However, it currently lacks the high-fidelity bulk process physics present in commercial TCAD tools such as Silvaco ATHENA and Synopsys Sentaurus Process.

The goal is to reach **functional parity** in the core areas of:
- Ion implantation (analytic + Monte Carlo)
- High-fidelity diffusion (multiple transport models, point defects, clustering, OED, stress effects, SPER, etc.)
- Coupled stress/mechanics
- Full process step orchestration with moving boundaries

We will achieve this by introducing a **unified multi-physics field layer** on top of the existing level-set geometry engine, while **preserving the existing advanced oxidation implementation** through well-defined adapter interfaces.

**Selected Path**: Path 1 (MFEM + amgcl + SUNDIALS unified framework) with high-fidelity doping core (C) and broad process completeness (B).

**Key Libraries**:
- MFEM (vcpkg): FEM discretization for bulk fields
- amgcl (3rdparty): AMG solvers/preconditioners
- SUNDIALS (vcpkg): Stiff time integration
- Eigen3 (vcpkg): Auxiliary linear algebra
- Existing Vienna tools + level-set geometry remain the foundation

## Selected Approaches and Trade-offs

**Overall Scope**: B (broader process completeness) + C (high-fidelity doping core)

**Phasing Strategy**: B (architecture-first + parallel development)

**Dimension Strategy**: B (2D/3D parallel from architecture start)

**Framework Style**: B (unified multi-physics framework)

**Oxidation Handling**: A (preserve existing implementation via interfaces)

**Recommended Library Stack**: MFEM-centric with amgcl + SUNDIALS (Path 1)

**Trade-offs**:
- Path 1 was chosen over lighter custom solvers for long-term maintainability, accuracy, and extensibility.
- Preserving existing oxidation (A) reduces risk and re-use of validated code (LOCOS, mask bending, stress) at the cost of a clean adapter layer.
- Starting with high-fidelity implant + diffusion + stress in parallel tracks allows early demonstration of value.

## Overall Architecture

ViennaPS will retain its core strengths:
- Level-set geometry evolution (ViennaLS)
- Monte Carlo flux calculation (ViennaRay)
- Existing high-quality Oxidation model (diffusion + viscous flow + stress + mask bending)

A new **Unified Physics Field Layer** will be introduced for bulk quantities:
- Dopant concentrations
- Point defects (I, V)
- Clusters and damage
- Stress/strain fields

**Layers**:
1. Geometry Layer (level-set, unchanged core)
2. Unified Physics Field Layer (new, MFEM-based)
3. Process Model Layer (new kernels + adapters)
4. Coupling & Orchestration Layer

The existing Oxidation model will be kept intact and integrated via a thin `OxidationFieldAdapter`.

## Unified Physics Field Layer

Inspired by MOOSE's Kernel + Material pattern (highly relevant reference from f:\dev\moose), we introduce:

- `PhysicsKernel` abstraction: each physical mechanism (diffusion term, reaction, stress divergence, etc.) is a separate kernel contributing to residual and Jacobian.
- `MaterialPropertySystem`: centralized management of temperature-, concentration-, and stress-dependent parameters (matching manual-style parameter inheritance and like-materials).
- Fields managed via MFEM `GridFunction`s on a background mesh (structured or adaptive).
- amgcl for efficient AMG solves.
- SUNDIALS for time integration of stiff systems.

Kernels are composable, allowing different diffusion models (Fermi, ChargedReact, Pair, etc.) and cluster models to be enabled independently.

## Geometry and Bulk Field Coupling

- Three level-sets continue to drive geometry: φ_Si (reaction interface), φ_amb (free surface), φ_mask (nitride).
- Background mesh for MFEM fields is marked dynamically using level-set information (oxide band, mask interior, etc.).
- Interface conditions are extracted from level-sets and applied to MFEM fields (Robin/Neumann for oxidant consumption, defect injection, stress traction).
- Constrained velocity for ambient interface based on mask (re-using existing logic).
- Data mapping and warm-starting between level-set interfaces and MFEM dofs.
- Boolean clipping and interior filling operations remain.

This allows realistic LOCOS bird's beak, mask bending, and lateral oxidant diffusion while fields evolve on the volume.

## Implantation Module

Supports both analytic and Monte Carlo paths, outputting directly to the unified fields.

**AnalyticImplantKernel**:
- Pearson / Dual-Pearson / Gaussian distributions
- Support for common table formats (Dios, Tasch, etc.)
- Damage models (Hobler-style initial I/V)
- Tilt/rotation, multilayer, screening layers

**MonteCarloImplantKernel**:
- BCA-based (crystal/amorphous)
- Damage accumulation
- Outputs dopant + defect + damage fields

Integration with existing GDS masks and level-set geometry for realistic masking and shadowing. Results seed the diffusion initial conditions.

## Diffusion Models and Solver

Pluggable kernels for high fidelity:

**Transport Kernels**:
- Fermi, ChargedReact, Pair, Constant, etc.
- Stress-dependent diffusivity and solubility

**Defect & Cluster Kernels**:
- Point defect evolution (I, V)
- Multiple cluster models (BIC, ChargedCluster, 311, Loop, FVCluster, Transient, etc.)
- Amorphous pocket and SPER handling

**Coupling Kernels**:
- OED sources
- Implant damage annealing
- Stress coupling

Solved using MFEM weak forms + amgcl + SUNDIALS time integration. Material system provides all parameters.

## Stress and Mechanics Module

Extends the existing oxidation stress logic via the adapter while adding general capabilities.

Kernels:
- Viscoelastic (Maxwell/SLS — re-using current oxidation implementation)
- Elastic (anisotropic)
- Plastic (incremental)
- Growth / mismatch stress

Outputs stress fields that feed back into diffusion kernels and the oxidation adapter. Supports thermal/lattice mismatch and rebalancing after etch/deposition.

## Time Stepping and Coupled Solver Orchestration

- SUNDIALS (CVODE/IDA) as primary time integrator for stiff diffusion + clustering + stress systems.
- Residual/Jacobian assembly from all active PhysicsKernels.
- Geometry (level-set) advancement driven by velocity fields (constrained for masks).
- Existing Oxidation called via adapter after or interleaved with field solves.
- Warm-starting of fields, stress history, and solver states.
- Combined CFL + error control for stability.

## Complete Process Step Loop and Integration with Existing Oxidation

A typical step follows:
1. Optional implantation (populates fields)
2. Field solve (diffusion + stress) via SUNDIALS + kernels + adapter
3. Geometry advancement (level-set)
4. Call existing Oxidation model via thin adapter (OED, defect injection, stress feedback)
5. Remap / warm-start fields to new geometry
6. Update material properties

The adapter exposes current dopant/defect/stress fields to the existing Oxidation and pulls back generated defects and updated stress without modifying its internal algorithms.

## Parallelism and 2D/3D Strategy

- Dimension-agnostic design from day one (templates + runtime dimension).
- MFEM assembly/solves use built-in OpenMP/CUDA.
- amgcl backends (builtin OpenMP, CUDA) for preconditioning.
- SUNDIALS parallel vector support.
- Level-set advection re-uses existing GPU/OpenMP paths.
- Background mesh marking limited to active bands for efficiency.
- 2D first for validation, 3D structures implemented in parallel.

## Implementation Phasing and Parallel Tracks

**Phase 0: Architecture Foundation** (serial)
- Unified field layer, kernel/material system, geometry coupling, adapter.
- Basic 2D/3D support.

**Parallel Track 1 (start mid-Phase 0)**: High-fidelity core
- Implantation (analytic + MC)
- Diffusion kernels + major cluster models
- Stress/mechanics kernels
- SUNDIALS integration

**Parallel Track 2**: Oxidation coupling
- Full adapter for OED, rates, defect injection, stress feedback

**Phase 1**: Orchestration & 3D
- Complete step loops
- 3D validation and performance tuning

**Phase 2**: Extensions
- Silicidation, advanced topography, lithography support, additional materials

## Dependencies

- MFEM, SUNDIALS, Eigen3 (via vcpkg)
- amgcl (already in 3rdparty)
- Existing ViennaPS dependencies (ViennaLS, ViennaRay, VTK, etc.)

## Validation Strategy

- Unit tests for individual kernels
- Regression against current oxidation behavior via adapter
- Comparison to manual examples (analytic profiles, LOCOS bird's beak with doping, stress effects, clustering)
- 2D first, then 3D

## Risks and Mitigations

- Integration precision with existing oxidation → thin adapter + dual validation
- Performance of MFEM on LS-driven problems → amgcl + band-limited solves + warm starts
- Parameter calibration for high-fidelity models → implement framework first, calibrate iteratively against manuals

This design provides a clean, extensible path to commercial-grade process physics while protecting ViennaPS's existing strengths in topography and oxidation.

---

**Next Steps** (after your review):
- Full spec self-review
- User review of this document
- Transition to detailed implementation plan (via writing-plans skill)

Please review the document above and let me know if anything needs adjustment before we commit it and move forward.