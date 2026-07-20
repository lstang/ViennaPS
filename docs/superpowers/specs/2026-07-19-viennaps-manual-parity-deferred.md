# Deferred items (manual parity plan)

Items from `2026-07-19-viennaps-manual-parity-implementation-plan.md` not finished in the multiphysics foundation work, and intentionally deferred:

## Deferred
- Full MC BCA implant (crystal channeling, full cascade fidelity) — sampling + Hobler damage stub exists
- Full MFEM vector elasticity (ElasticityIntegrator + traction BCs + stress tensor recovery)
- Production-scale SUNDIALS on full MFEM GridFunction dofs with amgcl preconditioning (profile-packed CVODE path exists)
- Live ViennaLS Domain material maps → MFEM region marking (GeometryFieldCoupler depth/LOCOS-like path exists)
- LOCOS bird’s-beak + doping parity regression vs ATHENA/SProcess manuals
- 3D production validation and performance tuning (amgcl CUDA, band-limited solves)
- Phase 2: silicidation, lithography, advanced SPER, full parameter DB, Python bindings polish
- Debug builds against Release-prebuilt mfem (known MD/MDd runtime mismatch)

## Done in this goal (reference)
- Multiphysics foundation commit
- Named cluster models: recomb, 311, bic, loop (+ combined)
- Oxidation adapter OED/dopant/stress field-only and Domain hooks
- Domain multiphysics smoke + cluster/adapter validation markers
- **PhysicsField packState/unpackState multi-species CVODE state**
- **SundialsTimeIntegrator field-dofs mode (addToFieldRHS)**
- **PairDiffusionKernel + ChargedReactKernel**
- **ElasticStressKernel (Hooke proxy + optional MFEM mass energy)**
- **GeometryFieldCoupler (depth stack + LOCOS-like mark + remap)**
- **MC implant nuclear/electronic sampling + Hobler I/V + amorph proxy**
- **ProcessOrchestrator multi-kernel field-dof loop**
