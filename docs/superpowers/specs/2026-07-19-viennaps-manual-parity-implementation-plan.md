# Implementation Plan: ViennaPS Functionality Parity with ATHENA/SProcess

**Date**: 2026-07-19
**Based on**: 2026-07-19-viennaps-manual-parity-design.md
**Approach**: Path 1 (MFEM + amgcl + SUNDIALS unified framework), high-fidelity core first, preserve existing Oxidation via adapter.

## High-Level Phasing

**Phase 0: Foundation (Architecture + Core Framework)**
- Duration estimate: 3-5 weeks (serial)
- Goal: Stable unified field layer, geometry coupling, adapter skeleton, basic 2D support.
- Parallel start for Track 1 after 60% complete.

**Parallel Track 1: High-Fidelity Core (Implant + Diffusion + Stress)**
- High priority for C (high-fid doping).
- 2D first, 3D structures in parallel.

**Parallel Track 2: Oxidation Integration**
- After Track 1 has basic fields.
- Focus on adapter for OED, coupling, etc.

**Phase 1: Orchestration & Validation**
- Full loops, 3D, basic end-to-end.

**Phase 2: Extensions**
- Additional models, features.

## Detailed Tasks

### Phase 0: Foundation

1. **Setup & Dependencies**
   - Add MFEM, SUNDIALS, Eigen3 to vcpkg.json (or confirm prebuilts).
   - Ensure amgcl (3rdparty) builds and links cleanly with MFEM.
   - Add CMake integration for new libs (use existing patterns from vtk.cmake, cpm).
   - Create `viennaps/fields/` or similar namespace for new layer.

2. **Unified Physics Field Layer**
   - Define `PhysicsKernel` base class (inspired by MOOSE, adapted to MFEM weak forms).
   - Implement `MaterialPropertySystem` for parameters (Arrhenius, stress-dependent, like-materials interpolation).
   - MFEM background mesh management (structured preferred initially for LS coupling).
   - Field container for dopants, defects (I/V), clusters, stress.
   - Basic kernel examples (e.g., simple diffusion kernel).

3. **Geometry-Field Coupling**
   - Extend level-set utilities to mark background mesh regions (oxide band, mask, Si).
   - Implement interface extraction for BCs (Robin for oxidant, Neumann for defects, traction for stress).
   - Data mapping (interpolation between LS and MFEM dofs).
   - Constrained ambient velocity logic.
   - Warm-start mechanisms for fields across geometry updates.

4. **Oxidation Adapter Skeleton**
   - Create `OxidationFieldAdapter` class.
   - Define interfaces for reading/writing fields (dopants, defects, stress) without changing existing Oxidation internals.
   - Basic sync for geometry changes.

5. **2D/3D Foundation**
   - Templated or dimension-parameterized code for kernels, coupling.
   - Ensure all new code supports both from start (use D template like existing).

**Validation for Phase 0**: Simple test where a field is defined, coupled to a static geometry, and basic solve runs. Existing oxidation still works unchanged.

### Parallel Track 1: High-Fidelity Core (Start after Phase 0 ~60%)

**Implantation**
- AnalyticImplantKernel: Implement Pearson/Dual-Pearson/Gaussian, table support, damage models.
- MCImplantKernel: Basic BCA or integrate suitable MC (leverage existing particle infra if possible).
- Output to unified fields (dopant + initial defects/damage).
- Geometry integration (GDS masks, tilt/rotation, shadowing via LS).
- Tests against manual examples (1D profiles, damage).

**Diffusion**
- Implement key Transport Kernels: Fermi, ChargedReact, Pair, etc.
- Defect kernels: I/V evolution, basic clustering (start with 311, BIC, simple loops).
- Material properties for all models (from manuals).
- Couple to stress (from Track 1 mechanics).
- OED source kernel (will connect to Track 2).
- Use SUNDIALS for time integration of stiff system.
- amgcl for linear solves.

**Stress/Mechanics**
- Mechanics Kernels: Viscoelastic (reuse logic from existing ox where possible), elastic, growth/mismatch.
- Output stress fields usable by diffusion.
- Initial coupling to geometry (via adapter later).
- Use MFEM for discretization.

**Integration within Track 1**
- End-to-end: Implant -> diffuse with defects/clustering/stress.
- 2D validation first (compare to manual profiles, simple TED/OED cases).
- Unit tests for each kernel.
- Performance: warm-start, limited domain solves.

### Parallel Track 2: Oxidation Integration (Start after basic fields from Track 1)

- Enhance OxidationFieldAdapter for:
  - OED: oxidation-enhanced defect generation.
  - Dopant-dependent oxidation rates.
  - Stress feedback to oxidation.
  - Defect injection from oxidation to unified fields.
- Test with existing LOCOS examples (ensure no regression in geometry/stress, while adding doping effects).
- Bird's beak with doping.

### Phase 1: Full Orchestration

- Implement complete Process step loops (implant -> diffuse -> oxidize, with geometry updates).
- Full time-step orchestration with SUNDIALS, level-set advance, adapter calls.
- 3D support and validation.
- Basic multi-step process examples.
- Error handling, logging for fields.

### Phase 2: Extensions & Parity Polish

- Silicidation models.
- Additional diffusion models (polysilicon, SiGe specifics, laser anneal).
- More cluster models, advanced SPER.
- Advanced structure generation features if gaps remain.
- Lithography support (if in scope).
- Full parameter DB like in manuals (interpolation, inheritance).
- Comprehensive validation suite against manual examples.
- Documentation, examples, Python bindings updates.

## Dependencies & Setup

- vcpkg: MFEM, SUNDIALS, Eigen3 (PETSc optional for advanced).
- 3rdparty: amgcl (already present, extend if needed for more backends).
- Ensure compatibility with existing ViennaLS, VTK, etc.
- Build system updates for new libs.

## Risks & Mitigations

- Coupling precision with existing Oxidation: Thin adapter + extensive regression tests on current examples.
- Performance (MFEM overhead): amgcl + warm starts + band-limited solves + profiling.
- Stiff systems stability: SUNDIALS error control + manual-inspired time step limits.
- 3D scaling: Start with 2D, use GPU where possible (amgcl CUDA).
- Validation against manuals: Use simple 1D/2D cases first, then complex LOCOS with doping.

## Validation & Testing

- Unit tests for kernels, materials, adapters.
- Regression: Existing oxidation examples must still pass.
- Parity tests: Compare output profiles, geometries, stresses to manual descriptions/examples.
- 2D first, then 3D.
- Performance benchmarks.

## Success Criteria

- Can simulate key process steps with high-fidelity doping and stress.
- Existing Oxidation functionality preserved and enhanced via coupling.
- Matches qualitative/quantitative behavior from ATHENA/SProcess manuals for targeted models.
- Clean, extensible code using chosen libs.

## Next Steps

1. Self-review of this plan (done).
2. User approval of full design + this plan.
3. Detailed task breakdown and start Phase 0 (use todo or issues).
4. Begin coding with tests.

This plan is ready for execution. All design choices from previous sections are incorporated. 

If any adjustments needed before starting, let me know. Otherwise, we can proceed to detailed tasking.