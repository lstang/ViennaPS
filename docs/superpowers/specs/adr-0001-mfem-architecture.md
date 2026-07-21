# ADR-0001: MFEM Architecture — Standalone Engine with Prebuilt mfem.lib

**Status:** Accepted
**Date:** 2026-07-20
**Decision owner:** ltang
**Supersedes:** Phase 1 Task 0 open question in
`docs/superpowers/plans/2026-07-20-diffusion-phase1.md`.

## Context

The diffusion engine plans (Phases 1–11) need MFEM as the FEM backend.
MOOSE ships a complete MFEM subsystem at
`3rdparty/moose/framework/include/mfem/` (gated by `MOOSE_MFEM_ENABLED`)
that already implements most of what the plans design from scratch:
`MFEMProblem`, `MFEMDiffusionKernel`, `MFEMNLDiffusionKernel`,
`MFEMHypreBoomerAMG`, `MFEMRefinementMarker`, `MFEML2ZienkiewiczZhuIndicator`,
`MFEMMesh`, `MFEMCutTransitionSubMesh`, etc.

Two viable architectures were identified:

- **(A) Standalone ViennaPS MFEM engine** — implement equivalent classes
  inside `include/viennaps/fields/` with full control of the public API.
- **(B) Reuse MOOSE's MFEM subsystem via `ExternalProblem`** — pull MOOSE
  in as a build/runtime dependency and wrap.

## Decision

**Adopt (A): Standalone ViennaPS MFEM engine.**

Refinement on the original plan: **not header-only.** ViennaPS links
against the **prebuilt MFEM shared/static library** at
`f:/dev/mfem/build` (MFEM 4.9.1, detected via `MFEMConfig.cmake`).
This is the existing convention in this fork — see
`AGENTS.md`: "prebuilt MFEM (provides `MFEMConfig.cmake`)" and the
CMakeLists hints at lines 67–70, 237–255.

**MOOSE classes are the design reference**, not a runtime dependency.
Every Phase 1 class has a verified MOOSE counterpart whose design choices
are inherited (see the table in Phase 1 Task 0 of the plan). When a
MOOSE class has a non-obvious correct behavior (e.g., `InterfaceReaction`'s
4-block Jacobian, `MultiSpeciesDiffusionCG`'s species-outer assembly
order, `MFEMRefinementMarker`'s separate `hRefine()`/`pRefine()`), the
ViennaPS implementation mirrors it.

## Rationale

1. **Build weight.** Pulling MOOSE as a runtime dep would dwarf the
   ViennaPS binary, break the lightweight-fork goal, and complicate
   distribution. Linking `mfem.lib` (already required) keeps the
   dependency surface unchanged.
2. **Header-only convention relaxed only where MFEM forces it.**
   ViennaPS itself remains header-only C++; the only `.lib` is MFEM,
   which was already a build dep before this ADR.
3. **API control.** MOOSE's input-file-driven Reporter/VPP/Postprocessor
   execution model does not match ViennaPS's existing C++/pybind11 API
   style. A standalone engine can compose cleanly with ViennaPS's
   existing `psProcess`, `GeometryFieldCoupler`, etc.
4. **Design inheritance without runtime coupling.** The MOOSE MFEM
   subsystem is the most thoroughly tested reference for "how to wrap
   MFEM correctly" — citing its classes forces correct behavior
   (4-block interface Jacobians, hp-refinement, ZZ error estimator,
   prolongation vs L2-projection distinction) without dragging in
   MOOSE's framework.

## Consequences

- Every new ViennaPS class must cite its MOOSE counterpart in the
  header doc-comment as a design reference. This is a hard convention
  for the diffusion engine code.
- When a MOOSE class solves a subtle problem (e.g.,
  `MultiAppProjectionTransfer::assembleL2` for cross-mesh L2 projection),
  the ViennaPS implementation must replicate the algorithm, not invent
  a new one.
- Future re-evaluation: if ViennaPS later needs the parts of MOOSE that
  are *not* MFEM (e.g., `stochastic_tools` Bayesian calibration from
  Phase 10 Tasks 12–15), this decision should be revisited — those
  modules are not portable without MOOSE as a dep.

## References

- Plan: `docs/superpowers/plans/2026-07-20-diffusion-phase1.md` Task 0
- MOOSE MFEM subsystem: `3rdparty/moose/framework/include/mfem/`
- `AGENTS.md` (machine-specific setup section)
- Root `CMakeLists.txt` lines 67–70, 237–255
