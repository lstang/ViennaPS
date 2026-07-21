#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <vcTestAsserts.hpp>
#include <fields/MeshAttributes.hpp>
#include <fields/DiffusionModel.hpp>
#include <fields/models/ConstantDiffusion.hpp>
#include <fields/DiffusionPhysics.hpp>
#include <fields/LevelSetToMesh.hpp>
#include <geometries/psMakePlane.hpp>
#include <psDomain.hpp>

#ifdef VIENNAPS_HAS_MFEM
using namespace viennaps;

void TestMeshAttributes() {
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");
  attrs.setAttributeName(2, "SiO2");
  VC_TEST_ASSERT(attrs.materialName(1) == "Si");
  VC_TEST_ASSERT(attrs.materialName(2) == "SiO2");
  VC_TEST_ASSERT(attrs.isMaterial(1, "Si"));
  VC_TEST_ASSERT(!attrs.isMaterial(1, "SiO2"));
  VC_TEST_ASSERT(attrs.numMaterials() == 2);
  VC_TEST_ASSERT(attrs.attributeOf("Si") == 1);
  VC_TEST_ASSERT(attrs.hasMaterial("SiO2"));
}

void TestDiffusionModelInterface() {
  struct TestModel : public DiffusionModel<double> {
    int numSpecies() const override { return 1; }
    std::vector<std::string> speciesNames() const override {
      return {"TestSpecies"};
    }
    std::vector<int> applicableAttributes() const override { return {1}; }
  };
  TestModel m;
  VC_TEST_ASSERT(m.numSpecies() == 1);
  VC_TEST_ASSERT(m.speciesNames()[0] == "TestSpecies");
}

void TestConstantDiffusion() {
  ConstantDiffusion<double> model("Boron");
  model.setDiffusivity(1e-13, 3.46);
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");
  model.setup(attrs, 1273.15);
  VC_TEST_ASSERT(model.numSpecies() == 1);
  VC_TEST_ASSERT(model.speciesNames()[0] == "Boron");
  double kB = 8.617333262145e-5;
  double expectedD = 1e-13 * std::exp(-3.46 / (kB * 1273.15));
  VC_TEST_ASSERT(std::abs(model.getDiffusivity() - expectedD) / expectedD < 1e-6);
}

// ---- Minimal stubs for Phase 2/3 models referenced by the composition test ----
// These exist ONLY so TestDiffusionPhysicsComposition can run in Phase 1.
// Real FermiDiffusion (Phase 2) and CddDiffusion (Phase 3 Task 10) will
// replace them. Each stub claims the same single species ("Boron") so the
// gatekeeper test exercises the double-dC/dt prevention path.
template <class NumericType>
class FermiDiffusionStub : public DiffusionModel<NumericType> {
public:
  explicit FermiDiffusionStub(const std::string& species = "Boron") {
    this->setName("FermiDiffusion");
    species_ = species;
  }
  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override { return {species_}; }

private:
  std::string species_;
};

template <class NumericType>
class CddDiffusionStub : public DiffusionModel<NumericType> {
public:
  CddDiffusionStub() { this->setName("CddDiffusion"); }
  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override { return {"Boron"}; }
};

void TestDiffusionPhysics() {
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addSpecies("Interstitial");
  VC_TEST_ASSERT(physics.numSpecies() == 2);
  VC_TEST_ASSERT(physics.hasSpecies("Boron"));

  // Register a model
  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-13, 3.46);
  physics.addModel(model);

  // BC specification
  physics.addNeumannBC("Boron", "surface", 0.0);  // zero flux
  physics.addDirichletBC("Boron", "bottom", 1e18);

  VC_TEST_ASSERT(physics.numModels() == 1);

  // Per-species BC lookup (not flat list)
  const auto& boronBCs = physics.boundaryConditions("Boron");
  VC_TEST_ASSERT(boronBCs.size() == 2);
  const auto& iBCs = physics.boundaryConditions("Interstitial");
  VC_TEST_ASSERT(iBCs.empty());  // Interstitial has no BCs yet
}

void TestDiffusionPhysicsComposition() {
  // Verify the composition gatekeepers: two models on the same species must
  // not create a double time-derivative. Mirrors MOOSE PhysicsBase::
  // shouldCreateTimeDerivative (framework/include/physics/PhysicsBase.h:237).
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addSpecies("Interstitial");

  // Both models touch Boron; both naively want a time derivative on it.
  // The physics must guarantee exactly ONE time derivative per species.
  auto fermi = std::make_shared<FermiDiffusionStub<double>>("Boron");
  auto cdd   = std::make_shared<CddDiffusionStub<double>>();  // composes PairTerm on Boron
  physics.addModel(fermi);
  physics.addModel(cdd);

  VC_TEST_ASSERT(physics.shouldCreateTimeDerivative("Boron", *fermi));
  // Second model on the same species must be denied the time derivative —
  // it can still contribute stiffness/reaction terms, but not dC/dt.
  VC_TEST_ASSERT(!physics.shouldCreateTimeDerivative("Boron", *cdd));
  VC_TEST_ASSERT(physics.shouldCreateTimeDerivative("Interstitial", *cdd));

  // species-existence gatekeeper
  VC_TEST_ASSERT(physics.variableExists("Boron"));
  VC_TEST_ASSERT(!physics.variableExists("Arsenic"));
}

void TestLevelSetToMesh2D() {
  // Build a 2D level-set domain with one material via MakePlane.
  // MakePlane produces a substrate extending in the +y direction centered on
  // the origin; it is the lightest geometry generator available and yields a
  // domain whose first level set's grid provides the bounds + gridDelta the
  // converter reads.
  auto domain = viennaps::Domain<double, 2>::New();
  viennaps::MakePlane<double, 2>(domain, /*gridDelta*/ 0.5, /*xExtent*/ 8.0,
                                 /*yExtent*/ 0.0, /*baseHeight*/ 0.0,
                                 /*periodic*/ false,
                                 /*material*/ viennaps::Material::Si)
      .apply();
  VC_TEST_ASSERT(domain->getLevelSets().size() == 1);
  VC_TEST_ASSERT(domain->getMaterialMap());

  viennaps::LevelSetToMeshConverter<double, 2> converter;
  auto [mesh, attrs] = converter.convert(*domain);

  VC_TEST_ASSERT(mesh != nullptr);
  VC_TEST_ASSERT(mesh->GetNV() > 0);
  VC_TEST_ASSERT(mesh->GetNE() > 0);

  // Element attributes must form a non-empty set (materials were tagged).
  std::set<int> attributes;
  for (int i = 0; i < mesh->GetNE(); ++i) {
    attributes.insert(mesh->GetAttribute(i));
  }
  VC_TEST_ASSERT(!attributes.empty());
  VC_TEST_ASSERT(attrs.numMaterials() >= 1);

  std::cout << "[level-set-to-mesh-check] nv=" << mesh->GetNV()
            << " ne=" << mesh->GetNE()
            << " attributes=" << attributes.size()
            << " materials=" << attrs.numMaterials() << "\n";
}

int main() {
  TestMeshAttributes();
  TestDiffusionModelInterface();
  TestConstantDiffusion();
  TestDiffusionPhysics();
  TestDiffusionPhysicsComposition();
  TestLevelSetToMesh2D();
  std::cout << "All diffusion tests passed.\n";
  return 0;
}
#else
int main() {
  std::cout << "MFEM not available, skipping.\n";
  return 0;
}
#endif
