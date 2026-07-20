#pragma once

/// ProcessOrchestrator - Phase 1 full multi-step process orchestration.
/// Sequence:
///   1. Implant (profile seed)
///   2. Diffuse (Sundials field-dofs + transport/defect/stress kernels)
///   3. Oxidize (PhysicsFieldAdapter OED path)
///   4. Geometry mark + remap (GeometryFieldCoupler)

#include "fields/PhysicsField.hpp"
#include "fields/MaterialPropertySystem.hpp"
#include "fields/PhysicsKernel.hpp"
#include "fields/SundialsTimeIntegrator.hpp"
#include "fields/DiffusionKernel.hpp"
#include "fields/FermiDiffusionKernel.hpp"
#include "fields/PairDiffusionKernel.hpp"
#include "fields/ChargedReactKernel.hpp"
#include "fields/StressKernel.hpp"
#include "fields/DefectClusterKernel.hpp"
#include "fields/GeometryFieldCoupler.hpp"
#include "process/psPhysicsFieldAdapter.hpp"

#include <iostream>
#include <memory>
#include <vector>

namespace viennaps {

template <class NumericType, int D>
class ProcessOrchestrator {
public:
  ProcessOrchestrator() = default;

  void runMultiStepExample(std::shared_ptr<PhysicsField<NumericType>> field,
                           std::shared_ptr<MaterialPropertySystem<NumericType>> mats,
                           NumericType totalTime = 300.0) {
    std::cout << "[ProcessOrchestrator] Starting full multi-step sequence (Phase 1)\n";

    // Step 1: Implant
    {
      std::vector<NumericType> implantProfile(64, static_cast<NumericType>(1e13 / 64));
      field->injectImplantProfile("Dopant", implantProfile);
      field->injectImplantProfile("Interstitial", implantProfile);
      field->injectImplantProfile("Vacancy",
                                  std::vector<NumericType>(64, static_cast<NumericType>(5e12 / 64)));
      std::cout << "  [Step 1] Implant completed\n";
    }

    // Step 2: Diffuse with Sundials on packed field state
    {
      SundialsTimeIntegrator<NumericType> integ;
      integ.setPhysicsField(field);
      integ.setUseFieldState(true);

      auto add = [&](std::shared_ptr<PhysicsKernel<NumericType>> k) {
        k->setPhysicsField(field);
        k->setMaterialProperties(mats);
        integ.addKernel(std::move(k));
      };

      add(std::make_shared<DiffusionKernel<NumericType>>("Dopant", 1273.15));
      add(std::make_shared<FermiDiffusionKernel<NumericType>>("Dopant", 1273.15));
      add(std::make_shared<PairDiffusionKernel<NumericType>>("Dopant", 1273.15));
      add(std::make_shared<ChargedReactKernel<NumericType>>(1273.15));
      add(std::make_shared<ViscoelasticStressKernel<NumericType>>(1273.15));
      add(std::make_shared<ElasticStressKernel<NumericType>>(1273.15));
      add(std::make_shared<DefectClusterKernel<NumericType>>(1273.15));

      integ.evolve(0.0, totalTime * 0.6, totalTime * 0.1);
      std::cout << "  [Step 2] Diffuse (Sundials field-dofs + kernels) completed\n";
    }

    // Step 3: Oxidize via adapter
    {
      PhysicsFieldAdapter<NumericType, D> adapter(field, mats);
      adapter.applyToOxidationFieldOnly();
      adapter.updateFromOxidationFieldOnly(NumericType(1));
      std::cout << "  [Step 3] Oxidize via adapter completed (OED + stress feedback)\n";
    }

    // Step 4: Geometry update + field remap
    {
      GeometryFieldCoupler<NumericType> geo(field);
      geo.markFromDepthInterfaces(NumericType(0.12), NumericType(0.04));
      geo.remapAfterGrowth(NumericType(0.02));
      field->remapAfterGeometryUpdate();
      std::cout << "  [Step 4] Geometry update + field remap completed\n";
    }

    std::cout << "[ProcessOrchestrator] Full multi-step sequence finished.\n";
  }
};

} // namespace viennaps
