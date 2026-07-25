# Architecture Decision Record: Standalone MFEM Engine vs. MOOSE MFEM Subsystem

**Decision:** We choose to implement a Standalone ViennaPS MFEM engine (Architecture A) instead of directly wrapping MOOSE's MFEM subsystem as a runtime dependency.

**Rationale:** This maintains ViennaPS as a lightweight, header-only library under `include/viennaps/` without coupling its public API to MOOSE's input-file-driven or Postprocessor paradigms. We will, however, extensively reference and mirror MOOSE's battle-tested internal patterns (e.g., separating `DiffusionPhysics` from `DiffusionEngine`, using per-species GridFunctions, and leveraging HypreBoomerAMG preconditioning) to inherit the architectural benefits of a proven parallel FEM system while preserving standalone independence.