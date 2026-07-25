# ADR-0004: Moving Si/SiO2 Interface for OED

**Status:** Accepted  
**Date:** 2026-07-25  
**Phase:** 4 Task 0

## Context

Oxidation-enhanced diffusion (OED) injects interstitials at a **moving** Si/SiO2
boundary. The Phase 4 plan requires choosing one of three MOOSE-verified idioms.

## Decision

**Default idiom: (A) Element subdomain relabeling**

Mirror MOOSE `ElementSubdomainModifierBase` / `ThresholdElementSubdomainModifier`:

- Elements flip attribute Si → SiO2 as oxidation progress crosses a threshold.
- Interface boundary is the face set between the two attributes.
- Field reinit on newly oxidized elements: polynomial neighbor extrapolation
  (MOOSE `POLYNOMIAL_NEARBY`) — Phase 4 implements a simplified global flip
  in `OedSource::relabelOxidized` for unit tests; production should flip only
  interface-adjacent elements.

**Secondary idiom: (C) Boundary-node ALE** for the free oxide surface when
smooth surface motion is needed without topology change
(`INSADDisplaceBoundaryBC` pattern).

**Not default: (B) MFEM ALE mesh displace** — reserved for small growth without
element inversion; harder to combine with sharp material change.

## Consequences

- `OedSource::registerWith` applies Neumann flux on the current interface
  boundary attribute string.
- `OedSource::assembleReaction` provides a uniform volumetric flux proxy when
  mesh boundary attrs are not yet wired (engine tests).
- Full ALE + L2 solution transfer across remesh is Phase 11 / MovingMeshHandler
  follow-up (gap analysis high priority).

## References

- MOOSE `framework/include/meshmodifiers/ElementSubdomainModifierBase.h`
- Phase 4 plan Task 0
- SProcess UG: Emulated Oxidation-Enhanced Diffusion
