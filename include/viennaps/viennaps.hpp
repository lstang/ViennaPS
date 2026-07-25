#pragma once

#include <psConstants.hpp>
#include <psDomain.hpp>
#include <psExtrude.hpp>
#include <psReader.hpp>
#include <psSlice.hpp>
#include <psUnits.hpp>
#include <psWriter.hpp>

#include <psPlanarize.hpp>

#include <geometries/psMakeFin.hpp>
#include <geometries/psMakeHole.hpp>
#include <geometries/psMakePlane.hpp>
#include <geometries/psMakeStack.hpp>
#include <geometries/psMakeTrench.hpp>

#include <gds/psGDSReader.hpp>

#include <process/psProcess.hpp>
#include <process/psSurfaceDiffusion.hpp>

#include <models/psCF4O2Etching.hpp>
#include <models/psCSVFileProcess.hpp>
#include <models/psDirectionalProcess.hpp>
#include <models/psFaradayCageEtching.hpp>
#include <models/psFluorocarbonEtching.hpp>
#include <models/psGeometricDistributionModels.hpp>
#include <models/psHBrO2Etching.hpp>
#include <models/psIonBeamEtching.hpp>
#include <models/psIsotropicProcess.hpp>
#include <models/psMultiParticleProcess.hpp>
#include <models/psNeutralTransport.hpp>
#include <models/psOxidation.hpp>
#include <models/psOxideRegrowth.hpp>
#include <models/psSF6C4F8Etching.hpp>
#include <models/psSF6O2Etching.hpp>
#include <models/psSelectiveEpitaxy.hpp>
#include <models/psSingleParticleALD.hpp>
#include <models/psSingleParticleProcess.hpp>
#include <models/psTEOSDeposition.hpp>
#include <models/psTEOSPECVD.hpp>
#include <models/psWetEtching.hpp>
#include <models/psAnalyticImplant.hpp>

#include <fields/PhysicsField.hpp>
#include <fields/MaterialPropertySystem.hpp>
#include <fields/PhysicsKernel.hpp>
#include <fields/SimpleDiffusionKernel.hpp>
#include <fields/DiffusionKernel.hpp>
#include <fields/FermiDiffusionKernel.hpp>
#include <fields/PairDiffusionKernel.hpp>
#include <fields/ChargedReactKernel.hpp>
#include <fields/SundialsTimeIntegrator.hpp>
#include <fields/StressKernel.hpp>
#include <fields/MfemElasticityKernel.hpp>
#include <fields/DefectClusterKernel.hpp>
#include <fields/GeometryFieldCoupler.hpp>
#include <fields/AmgclSolver.hpp>
#include <fields/BandLimitedSolver.hpp>
#include <fields/ParameterDatabase.hpp>
#include <fields/SPERKernel.hpp>
#include <fields/LocosDopingValidator.hpp>
#include <ProcessOrchestrator.hpp>
#include <process/psPhysicsFieldAdapter.hpp>
#include <models/psBasicDiffusion.hpp>
#include <models/psMCBcaImplant.hpp>
#include <models/psSilicidation.hpp>
#include <models/psLithography.hpp>

#include <fields/MeshAttributes.hpp>
#include <fields/DiffusionModel.hpp>
#include <fields/DiffusionPhysics.hpp>
#include <fields/IntrinsicCarrier.hpp>
#include <fields/DiffusivityMaterial.hpp>
#include <fields/KernelTerm.hpp>
#include <fields/KernelTerms.hpp>
#include <fields/PointDefectEquilibrium.hpp>
#include <fields/models/ConstantDiffusion.hpp>
#include <fields/models/FermiDiffusion.hpp>
#include <fields/models/ChargedFermiDiffusion.hpp>
#include <fields/models/SolidSolubility.hpp>
#include <fields/models/Segregation.hpp>
#include <fields/models/ReactDiffusion.hpp>
#include <fields/models/ChargedReactDiffusion.hpp>
#include <fields/models/PairDiffusion.hpp>
#include <fields/models/ChargedPairDiffusion.hpp>
#include <fields/models/NeutralReactDiffusion.hpp>
#include <fields/models/Cluster311.hpp>
#include <fields/models/VacancyCluster.hpp>
#include <fields/models/ImpurityCluster.hpp>
#include <fields/models/DislocationLoop.hpp>
#include <fields/models/CddDiffusion.hpp>
#ifdef VIENNAPS_HAS_MFEM
#include <fields/LevelSetToMesh.hpp>
#include <fields/DiffusionEngine.hpp>
#endif

// These macros might be defined on some systems (MSCV), undefine them to avoid
// conflicts
#ifdef ERROR
#undef ERROR
#endif
#ifdef WARNING
#undef WARNING
#endif
#ifdef INFO
#undef INFO
#endif
#ifdef DEBUG
#undef DEBUG
#endif
