# ADR-0004: Moving-Interface Idiom for Si/SiO2 Oxidation

**Status:** Accepted  
**Date:** 2026-07-25  
**Context:** Phase 4 Task 0 — OED interstitial injection tracks a moving Si/SiO2 boundary.

## Decision

**Default idiom: (A) subdomain relabeling**, mirroring MOOSE
`ElementSubdomainModifierBase` / `ThresholdElementSubdomainModifier`.

- Elements flip Si → SiO2 as oxide grows.
- Moving boundary is tracked via material-attribute pairs (attribute map on
  `MeshAttributes`), not a fixed Neumann segment.
- Field reinitialization on newly oxidized elements uses a
  nearest-neighbor / polynomial-nearby strategy (MOOSE
  `ReinitStrategy::POLYNOMIAL_NEARBY`).

**Secondary:** free-surface ALE (idiom C, `INSADDisplaceBoundaryBC` pattern)
may be layered later for smooth oxide free-surface motion. Idiom B (full
MFEM mesh displace) is reserved for small smooth growth without topology
change.

## Consequences

- `OedSource` queries the current interface faces each step from mesh
  attributes (not a static boundary id).
- Full mesh topology updates remain out of scope for the Phase 4 unit-test
  skeleton; host-side flux APIs are provided first.
