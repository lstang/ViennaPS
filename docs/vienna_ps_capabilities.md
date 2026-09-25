# ViennaPS Current Capabilities Inventory

Generated from `include/viennaps/*` directories.

## Models

### psAnalyticImplant.hpp
- Classes: NumericType, AnalyticImplant
- Brief: psAnalyticImplant - Analytic ion implantation model for ViennaPS

### psBasicDiffusion.hpp
- Classes: NumericType, BasicDiffusion
- Brief: BasicDiffusion - Early stub diffusion process model.

### psCF4O2Etching.hpp
- Classes: CF4O2SurfaceModel, CF4O2Ion, CF4O2Etchant, CF4O2Oxygen, CF4O2Polymer

### psCF4O2Parameters.hpp
- Classes: 

### psCMP.hpp
- Classes: NumericType, CmpVelocityField, CMP
- Brief: CMP — chemical-mechanical planarization (Preston law) process model.

### psCSVFileProcess.hpp
- Classes: VelocityFieldFromFile, CSVFileProcess

### psDirectionalProcess.hpp
- Classes: NumericType, NumericType, DirectionalVelocityField, NumericType, DirectionalVelocityFieldSimple

### psFaradayCageEtching.hpp
- Classes: PeriodicSource, FaradayCageEtching, FaradayCageEtching

### psFluorocarbonEtching.hpp
- Classes: FluorocarbonSurfaceModel, FluorocarbonIon, FluorocarbonNeutral, FluorocarbonEtching

### psGeometricDistributionModels.hpp
- Classes: SphereDistribution, BoxDistribution, CustomSphereDistribution, T, TrenchDistribution

### psHBrO2Etching.hpp
- Classes: HBrO2Etching, HBrO2Etching

### psIonBeamEtching.hpp
- Classes: IBESurfaceModel, IBEIonWithRedeposition, IonBeamEtching, IonBeamEtching

### psIonBeamParameters.hpp
- Classes: 

### psIonModelUtil.hpp
- Classes: NumericType, NumericType, NumericType

### psIsotropicProcess.hpp
- Classes: NumericType, IsotropicVelocityField, IsotropicProcess

### psLithography.hpp
- Classes: NumericType, LithographyModel
- Brief: psLithography - Simple aerial-image / threshold lithography for mask generation.

### psMCBcaImplant.hpp
- Classes: NumericType, NumericType, MCBcaEngine, NumericType, MCBcaImplant
- Brief: psMCBcaImplant - Binary Collision Approximation Monte Carlo implant.

### psMultiParticleProcess.hpp
- Classes: MultiParticleSurfaceModel, IonParticle, DiffuseParticle, MultiParticleProcess, MultiParticleProcess

### psNeutralTransport.hpp
- Classes: NeutralTransportSurfaceModel, NeutralTransportParticle, NeutralTransport, NeutralTransport

### psOxidation.hpp
- Classes: OxidantType, SiliconOrientation, NumericType, Oxidation
- Brief: psOxidation — ViennaPS process model for thermal silicon oxidation.

### psOxideRegrowth.hpp
- Classes: NumericType, SelectiveEtchingVelocityField, NumericType, RedepositionVelocityField, T

### psPipelineParameters.hpp
- Classes: 

### psPlasmaEtching.hpp
- Classes: PlasmaEtchingSurfaceModel, PlasmaEtchingIon, PlasmaEtchingNeutral

### psPlasmaEtchingParameters.hpp
- Classes: 

### psSF6C4F8Etching.hpp
- Classes: SF6C4F8Etching, SF6C4F8Etching

### psSF6O2Etching.hpp
- Classes: SF6O2Etching, SF6O2Etching

### psSelectiveEpitaxy.hpp
- Classes: NumericType, EpitaxyVelocityField, SelectiveEpitaxy

### psSilicidation.hpp
- Classes: NumericType, SilicidationModel
- Brief: psSilicidation - Metal-silicon reaction forming silicide (TiSi2, NiSi, ...).

### psSingleParticleALD.hpp
- Classes: SingleParticleALDSurfaceModel, SingleParticleALDParticle, SingleParticleALD, SingleParticleALD

### psSingleParticleProcess.hpp
- Classes: SingleParticleSurfaceModel, SingleParticle, SingleParticleProcess, SingleParticleProcess

### psTEOSDeposition.hpp
- Classes: NumericType, SingleTEOSSurfaceModel, NumericType, MultiTEOSSurfaceModel, NumericType

### psTEOSPECVD.hpp
- Classes: NumericType, PECVDSurfaceModel, Ion, NumericType, TEOSPECVD

### psWetEtching.hpp
- Classes: NumericType, WetEtchingVelocityField, WetEtching

## Process / Strategy

### psALPStrategy.hpp
- Classes: ALPStrategy

### psAdvectionCallback.hpp
- Classes: AdvectionCallback

### psAdvectionHandler.hpp
- Classes: AdvectionHandler

### psAnalyticProcessStrategy.hpp
- Classes: AnalyticProcessStrategy

### psCPUDiskEngine.hpp
- Classes: CPUDiskEngine

### psCPUTriangleEngine.hpp
- Classes: CPUTriangleEngine

### psCallbackOnlyStrategy.hpp
- Classes: CallbackOnlyStrategy

### psCoverageManager.hpp
- Classes: CoverageManager

### psDesorptionSource.hpp
- Classes: DesorptionSource

### psFluxEngine.hpp
- Classes: FluxEngine

### psFluxProcessStrategy.hpp
- Classes: FluxProcessStrategy

### psGPUDiskEngine.hpp
- Classes: GPUDiskEngine

### psGPULineEngine.hpp
- Classes: GPULineEngine, AT, BT

### psGPUTriangleEngine.hpp
- Classes: GPUTriangleEngine

### psGeometricModel.hpp
- Classes: GeometricModel

### psGeometricProcessStrategy.hpp
- Classes: GeometricProcessStrategy

### psOxidationStrategy.hpp
- Classes: OxidationStrategy

### psPhysicsFieldAdapter.hpp
- Classes: Domain, NumericType, PhysicsFieldAdapter
- Brief: PhysicsFieldAdapter - Thin coupling between PhysicsField and existing Oxidation.

### psProcess.hpp
- Classes: Process

### psProcessContext.hpp
- Classes: ProcessResult

### psProcessModel.hpp
- Classes: ProcessModelBase, ProcessModelCPU, ProcessModelGPU

### psProcessParams.hpp
- Classes: ProcessParams

### psProcessStrategy.hpp
- Classes: ProcessStrategy

### psSurfaceDiffusion.hpp
- Classes: NumericType, NumericType, NeighborSearch, NumericType, SurfaceDiffusionStencil

### psSurfaceModel.hpp
- Classes: SurfaceModel

### psTranslationField.hpp
- Classes: TranslationField

### psVelocityField.hpp
- Classes: VelocityField, DefaultVelocityField

## Fields / Multiphysics

### AdaptiveMeshRefiner.hpp
- Classes: MeshQualityEstimator, AdaptiveMeshRefiner
- Brief: Adaptive mesh refinement helpers (Phase 11 skeleton).

### AmgclSolver.hpp
- Classes: NumericType
- Brief: AmgclSolver - Thin wrapper around amgcl for sparse CSR systems.

### BandLimitedSolver.hpp
- Classes: NumericType, BandLimitedSolver
- Brief: BandLimitedSolver - Restrict field solves to an active band (3D scaling).

### BandgapModel.hpp
- Classes: NumericType, BandgapModel
- Brief: BandgapModel — E_g(x_Ge, T, strain) for SiGe and related materials.

### ChargedReactKernel.hpp
- Classes: NumericType, ChargedReactKernel
- Brief: ChargedReactKernel - charged point-defect reaction model (I-V recomb + charge).

### DefectClusterKernel.hpp
- Classes: NumericType, DefectClusterKernel
- Brief: DefectClusterKernel - Point defect clustering models for high-fid diffusion.

### DiffusionEngine.hpp
- Classes: NumericType, DiffusionEngine, DiffusionRHSOperator, body, NumericType
- Brief: @file DiffusionEngine.hpp

### DiffusionKernel.hpp
- Classes: NumericType, DiffusionKernel
- Brief: DiffusionKernel - Transport kernel for dopant/defect diffusion.

### DiffusionModel.hpp
- Classes: NumericType, DiffusionModel

### DiffusionPhysics.hpp
- Classes: NumericType, DiffusionPhysics

### DiffusivityMaterial.hpp
- Classes: NumericType, FermiDCoef, FermiDdCCoef
- Brief: DiffusivityMaterial — concentration-dependent diffusivity coefficients

### FermiDiffusionKernel.hpp
- Classes: NumericType, FermiDiffusionKernel
- Brief: FermiDiffusionKernel - Fermi (extrinsic) diffusion model.

### FlashAnnealFlow.hpp
- Classes: NumericType, FlashAnnealFlow
- Brief: FlashAnnealFlow — end-to-end flash/laser anneal orchestration on

### GeometryFieldCoupler.hpp
- Classes: NumericType, GeometryFieldCoupler
- Brief: GeometryFieldCoupler - Level-set / material-map ↔ PhysicsField mesh coupling.

### GrainBoundaryMesh.hpp
- Classes: GrainBoundaryMesh
- Brief: GrainBoundaryMesh — dual mesh (grain interior + GB network).

### GrainModel.hpp
- Classes: NumericType, GrainModel
- Brief: GrainModel — average grain radius evolution for polysilicon.

### ImplantDamageCoupler.hpp
- Classes: NumericType, ImplantDamageCoupler
- Brief: ImplantDamageCoupler — bridge MCBca implant damage into DiffusionEngine

### IntrinsicCarrier.hpp
- Classes: CarrierStatistics, NumericType, IntrinsicCarrier
- Brief: IntrinsicCarrier — pure-math calculator for semiconductor intrinsic and

### KernelTerm.hpp
- Classes: KernelTerm, NumericType, EquilibriumSpeciesAuxKernel
- Brief: KernelTerm — one physics contribution (MOOSE Kernel pattern).

### KernelTerms.hpp
- Classes: DiffusionTerm, ReactionTerm, ScaledGFCoef, CoupledForceTerm, ScaledGFCoef
- Brief: Built-in KernelTerm implementations (MOOSE MatDiffusion / Reaction /

### LevelSetToMesh.hpp
- Classes: NumericType, LevelSetToMeshConverter, NumericType, LevelSetToMeshConverter, NumericType
- Brief: @file LevelSetToMesh.hpp

### LocosDopingValidator.hpp
- Classes: NumericType, NumericType, LocosDopingValidator
- Brief: LocosDopingValidator - Synthetic LOCOS bird's-beak + doping regression.

### MaterialConverter.hpp
- Classes: MaterialConverter
- Brief: MaterialConverter — re-tag mesh materials (e.g. Si → GaAs).

### MaterialPropertySystem.hpp
- Classes: NumericType, MaterialPropertySystem
- Brief: MaterialPropertySystem - Centralized, queryable material properties for

### MeshAttributes.hpp
- Classes: MeshAttributes

### MfemElasticityKernel.hpp
- Classes: NumericType, MfemElasticityKernel
- Brief: MfemElasticityKernel - Vector isotropic elasticity with MFEM assembly.

### MovingMeshHandler.hpp
- Classes: MovingMeshHandler
- Brief: MovingMeshHandler — mesh mutation for oxidation-style progress (ADR-0004).

### PairDiffusionKernel.hpp
- Classes: NumericType, PairDiffusionKernel
- Brief: PairDiffusionKernel - dopant-interstitial pair diffusion model.

### ParameterDatabase.hpp
- Classes: NumericType, ParameterDatabase
- Brief: ParameterDatabase - Manual-style material parameter DB with inheritance

### PdeApi.hpp
- Classes: Type, PdeTerm, DiffusionPdeTerm, ReactionPdeTerm, PdeEquation
- Brief: PDE API — composable equation terms + results extraction (Phase 10).

### PhysicsField.hpp
- Classes: NumericType, PhysicsField
- Brief: PhysicsField - Container for bulk physical fields (dopants, defects, damage, stress).

### PhysicsKernel.hpp
- Classes: NumericType, PhysicsField, NumericType, MaterialPropertySystem, NumericType
- Brief: PhysicsKernel base for the unified multiphysics framework.

### PointDefectEquilibrium.hpp
- Classes: NumericType, PointDefectEquilibrium
- Brief: PointDefectEquilibrium — C_I^eq(T), C_V^eq(T) for Si and other materials.

### SPERKernel.hpp
- Classes: NumericType, SPERKernel
- Brief: SPERKernel - Solid Phase Epitaxial Regrowth of amorphous Si.

### SimpleDiffusionKernel.hpp
- Classes: NumericType, SimpleDiffusionKernel
- Brief: SimpleDiffusionKernel - A minimal concrete PhysicsKernel example.

### SolutionTransfer.hpp
- Classes: SolutionTransfer, SourceSampleCoef
- Brief: SolutionTransfer — field transfer between meshes with dose preservation.

### StressKernel.hpp
- Classes: NumericType, ViscoelasticStressKernel, NumericType, ElasticStressKernel
- Brief: Stress kernels for multiphysics process simulation.

### SundialsTimeIntegrator.hpp
- Classes: NumericType, NumericType, SundialsTimeIntegrator
- Brief: SundialsTimeIntegrator - SUNDIALS (CVODE) time integration for multiphysics fields.

## GDS / Layout

### psGDSGeometry.hpp
- Classes: NumericType, GDSGeometry

### psGDSMaskProximity.hpp
- Classes: GDSMaskProximity

### psGDSReader.hpp
- Classes: reads, GDSReader

### psGDSUtils.hpp
- Classes: ElementType, RecordNumbers, T, T, T

## Geometries

### psGeometryFactory.hpp
- Classes: NumericType, GeometryFactory

### psMakeFin.hpp
- Classes: NumericType, MakeFin

### psMakeHole.hpp
- Classes: HoleShape, NumericType, MakeHole

### psMakePlane.hpp
- Classes: provides, NumericType, MakePlane
- Brief: This class provides a simple way to create a plane in a level set. It can be

### psMakeStack.hpp
- Classes: NumericType, MakeStack

### psMakeTrench.hpp
- Classes: NumericType, MakeTrench

## Materials

### psBuiltInMaterial.hpp
- Classes: MaterialCategory, BuiltInMaterial

### psMaterial.hpp
- Classes: Material, Kind

### psMaterialMap.hpp
- Classes: that, and, MaterialMap, T, T

### psMaterialRegistry.hpp
- Classes: MaterialRegistry

### psMaterialValueMap.hpp
- Classes: T, MaterialValueMap, MapLike, Iterator, Phase

## Compact

### psCSVDataSource.hpp
- Classes: CSVDataSource

### psCSVReader.hpp
- Classes: for, NumericType, CSVReader

### psCSVWriter.hpp
- Classes: Iterator, template, NumericType, CSVWriter

### psDataScaler.hpp
- Classes: for, DataScaler, StandardScaler, MedianDistanceScaler

### psDataSource.hpp
- Classes: DataSource, has, has

### psNearestNeighborsInterpolation.hpp
- Classes: NearestNeighborsInterpolation

### psRectilinearGridInterpolation.hpp
- Classes: RectilinearGridInterpolation

### psValueEstimator.hpp
- Classes: ValueEstimator
