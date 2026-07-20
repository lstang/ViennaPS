# Deferred items (manual parity plan)

Items from `2026-07-19-viennaps-manual-parity-implementation-plan.md` not finished in the current multiphysics commit goal, and intentionally deferred:

## Deferred
- Full MC BCA implant (crystal channeling, damage cascade fidelity) beyond sampling stub
- Full MFEM elasticity assembly (ElasticityIntegrator + traction BCs + stress recovery tensor fields)
- Full SUNDIALS residual on MFEM GridFunction dofs (state = concatenated species dofs)
- Real level-set ↔ MFEM region marking from Domain material maps / interface BCs
- LOCOS bird’s-beak + doping parity regression vs ATHENA/SProcess manuals
- 3D production validation and performance tuning (amgcl CUDA, band-limited solves)
- Phase 2: silicidation, lithography, advanced SPER, full parameter DB, Python bindings polish
- Debug builds against Release-prebuilt mfem (known MD/MDd runtime mismatch)

## Done in this goal (reference)
- Multiphysics foundation commit
- Named cluster models: recomb, 311, bic, loop (+ combined)
- Oxidation adapter OED/dopant/stress field-only and Domain hooks
- Domain multiphysics smoke + cluster/adapter validation markers
