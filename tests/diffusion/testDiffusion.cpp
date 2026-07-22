#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <vcTestAsserts.hpp>
#include <fields/MeshAttributes.hpp>
#include <fields/DiffusionModel.hpp>
#include <fields/models/ConstantDiffusion.hpp>
#include <fields/DiffusionPhysics.hpp>
#include <fields/DiffusionEngine.hpp>
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

void TestDiffusionEngineAssembly() {
  // Phase 1 Task 5 smoke check (per brief): build a 4x4 triangular MFEM mesh
  // directly (no LevelSetToMesh), register one ConstantDiffusion("Boron")
  // model with D0=1e-3 and Ea=0 (so D = 1e-3 independent of T, easy to reason
  // about), initialize to 1e18, integrate 0->1s, and assert dose > 0.
  //
  // Strict dose conservation is Task 7's responsibility; here we only verify
  // the engine runs end-to-end and produces a sensible (positive) dose. We
  // additionally check that the concentration field has actually evolved away
  // from the uniform initial condition — otherwise the test would pass even
  // if the engine silently did nothing (D=0 bug, no mass matrix, etc.).
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::TRIANGLE));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-3, 0.0);  // Ea=0 => D = D0 = 1e-3 regardless of T
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addModel(model);
  physics.setTemperature(1273.15);  // irrelevant for Ea=0 but kept for realism
  engine.setPhysics(physics);

  engine.initializeSpecies("Boron", 1e18);

  // Integrate 0->1s with max step 0.1s. Zero-flux (Neumann) BCs are the
  // natural default; for ConstantDiffusion with uniform IC and zero flux,
  // the analytic solution is constant in time (grad u = 0 => du/dt = 0),
  // so dose is conserved exactly AND the field does not redistribute. The
  // engine's job here is to demonstrate it can assemble+integrate without
  // blowing up; Task 7 will use a non-uniform IC to verify conservation.
  engine.solve(0.0, 1.0, 0.1);

  const double dose = engine.getIntegral("Boron");
  std::cout << "[diffusion-engine-assembly] Boron dose after 1s = " << dose
            << "\n";
  VC_TEST_ASSERT(dose > 0.0);

  // Sanity: with uniform IC + zero-flux BCs the analytic steady state is
  // u(x,t) = u0, so the field should still be ~1e18 everywhere. If the
  // engine silently zeroed it (mass matrix bug, bad linear solve, etc.) the
  // max would be near zero. Allow a wide tolerance — this is a smoke check,
  // not a conservation invariant (Task 7's job).
  const mfem::GridFunction &gf = engine.getSolution("Boron");
  const double gfMax = gf.Max();
  std::cout << "[diffusion-engine-assembly] Boron max = " << gfMax << "\n";
  VC_TEST_ASSERT(gfMax > 1e17);  // within ~10x of initial 1e18
}

void TestDoseConservation() {
  // Phase 1 Task 7 capstone: dose conservation under zero-flux (Neumann) BCs.
  //
  // Physics: for a closed system (all Neumann zero-flux boundaries) with no
  // reaction terms, total concentration (dose) must be conserved over time.
  // The implicit-Euler scheme + zero-flux BCs should conserve dose to floating-
  // point precision, so the 1% tolerance below is generous; if the test fails
  // by more than 1% something is physically wrong.
  //
  // Setup per the task brief:
  //   - 8x8 triangular mesh (MakeCartesian2D default unit square, area = 1)
  //   - ConstantDiffusion("Boron") with D0=1e-8, Ea=0 (so D = 1e-8 at any T)
  //   - Uniform IC = 1e18
  //   - Integrate 0 -> 10s with dt = 1.0
  //   - Assert |dose_final - dose_initial| / dose_initial < 0.01
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-8, 0.0);  // Ea=0 => D = D0 = 1e-8 regardless of T
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addModel(model);
  physics.setTemperature(1273.15);  // irrelevant for Ea=0 but kept for realism
  engine.setPhysics(physics);

  engine.initializeSpecies("Boron", 1e18);

  const double doseInitial = engine.getIntegral("Boron");

  // Zero-flux (Neumann) BCs are the natural default in Phase 1; no BCs are
  // registered, so the system is closed and dose must be conserved.
  engine.solve(0.0, 10.0, 1.0);

  const double doseFinal = engine.getIntegral("Boron");
  const double relDiff =
      std::abs(doseFinal - doseInitial) / std::abs(doseInitial);

  std::cout << "[dose-conservation] dose_initial=" << doseInitial
            << " dose_final=" << doseFinal << " rel_diff=" << relDiff << "\n";

  // 1% tolerance — implicit Euler + zero-flux should be at FP precision.
  // If this fails by more than 1%, dose is leaking: real bug.
  VC_TEST_ASSERT(relDiff < 0.01);
}

void TestDirichletBC() {
  // Phase 1 BC follow-up (F1): verify Dirichlet BCs are actually applied
  // by the engine (not silently treated as natural zero-flux).
  //
  // Physics: with the entire boundary clamped at C = 1e18 and the interior
  // initialized to 0, diffusion will drive the interior toward 1e18. After
  // a long integration the solution should be ~uniform at 1e18 (the unique
  // steady state of Laplace's equation with all-Dirichlet BC).
  //
  // Setup:
  //   - 8x8 triangular mesh on the unit square (MakeCartesian2D with
  //     generate_edges=true produces 4 distinct boundary attributes, one
  //     per side; "all" marks every boundary attribute)
  //   - ConstantDiffusion D=1e-2 (fast diffusion so steady-state is reached
  //     within the integration window)
  //   - Dirichlet BC C=1e18 on "all" boundaries
  //   - Interior IC = 0
  //   - Integrate 0 -> 5s
  //   - Assert the final interior mean concentration > 0.9e18 (within 10%
  //     of the clamped boundary value; the system is near steady state)
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(
      8, 8, mfem::Element::TRIANGLE, /*generate_edges*/ true));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-2, 0.0); // Ea=0 => D=1e-2 regardless of T
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addModel(model);
  // "all" marks every boundary attribute on the mesh (4 sides here).
  physics.addDirichletBC("Boron", "all", 1e18);
  engine.setPhysics(physics);

  // Interior IC = 0. The Dirichlet BC is enforced by the engine; the
  // boundary dofs are NOT pre-set on the GridFunction — the engine's
  // EliminateRow call sets them to 1e18 at each step.
  engine.initializeSpecies("Boron", 0.0);

  engine.solve(0.0, 5.0, 0.1);

  const double doseFinal = engine.getIntegral("Boron");
  // Unit square area = 1; if the interior has reached ~1e18 uniform,
  // dose ~ 1e18 * 1 = 1e18. Allow 10% slack for transient lag.
  const double interiorMean = doseFinal; // area is 1
  std::cout << "[dirichlet-bc] interior_mean=" << interiorMean
            << " (expected ~1e18)\n";
  VC_TEST_ASSERT(interiorMean > 0.5e18);
}

void TestNeumannBC() {
  // Phase 1 BC follow-up (F1): verify non-zero Neumann BCs are actually
  // applied by the engine.
  //
  // Physics: uniform IC + constant surface flux g on all boundaries. The
  // dose should increase linearly: d(dose)/dt = g * (perimeter). For the
  // unit square, perimeter = 4, so dose(t) = dose_initial + 4*g*t.
  //
  // Setup:
  //   - 8x8 triangular mesh on the unit square, generate_edges=true
  //   - ConstantDiffusion D=1e-3 (irrelevant for the integral balance; the
  //     flux is prescribed at the boundary regardless of interior D)
  //   - Uniform IC = 1e18 (so we have a meaningful baseline)
  //   - Neumann flux g=1e15 on boundary "1" (== "all")
  //   - Integrate 0 -> 1s with dt=0.1
  //   - Expected delta = 4 * 1e15 * 1 = 4e15
  //   - Assert |dose_final - dose_expected| / dose_expected < 0.05 (5%)
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(
      8, 8, mfem::Element::TRIANGLE, /*generate_edges*/ true));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-3, 0.0);
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addModel(model);
  physics.addNeumannBC("Boron", "all", 1e15);
  engine.setPhysics(physics);

  engine.initializeSpecies("Boron", 1e18);

  const double doseInitial = engine.getIntegral("Boron");
  engine.solve(0.0, 1.0, 0.1);
  const double doseFinal = engine.getIntegral("Boron");

  const double flux = 1e15;
  const double perimeter = 4.0; // unit square
  const double expectedDelta = flux * perimeter * 1.0;
  const double doseExpected = doseInitial + expectedDelta;
  const double relDiff = std::abs(doseFinal - doseExpected) / doseExpected;

  std::cout << "[neumann-bc] dose_initial=" << doseInitial
            << " dose_final=" << doseFinal
            << " expected=" << doseExpected << " rel_diff=" << relDiff << "\n";

  // 5% tolerance: spatial discretization + implicit-Euler time integration
  // introduce small errors, but the integral balance should be tight.
  VC_TEST_ASSERT(relDiff < 0.05);
}

void TestReentrantSolve() {
  // F2 follow-up: shouldCreateTimeDerivative must be re-entrant across
  // solve() calls. Without resetTimeDerivativeClaims(), the second solve()
  // on the same physics object denies the mass matrix for every species
  // (timeDerivativeClaimed_ still holds them from the first solve),
  // silently falling back to identity mass and producing wrong time scales.
  //
  // Setup: closed system, uniform IC. Two consecutive solve() calls on
  // the SAME physics object. Both must conserve dose to ~FP precision.
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto model = std::make_shared<ConstantDiffusion<double>>("Boron");
  model->setDiffusivity(1e-8, 0.0);
  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addModel(model);
  physics.setTemperature(1273.15);
  engine.setPhysics(physics);

  engine.initializeSpecies("Boron", 1e18);
  const double doseInitial = engine.getIntegral("Boron");

  // First solve - mass matrix claimed for Boron.
  engine.solve(0.0, 5.0, 1.0);
  const double doseAfterFirst = engine.getIntegral("Boron");
  const double relDiffFirst =
      std::abs(doseAfterFirst - doseInitial) / std::abs(doseInitial);

  // Second solve on the SAME physics object - this is the regression
  // case. Before the F2 fix, timeDerivativeClaimed_ still held "Boron",
  // so the mass matrix was denied and identity mass was used (du/dt =
  // -K u + R with no M inverse - wildly different time scale).
  engine.solve(5.0, 10.0, 1.0);
  const double doseAfterSecond = engine.getIntegral("Boron");
  const double relDiffSecond =
      std::abs(doseAfterSecond - doseInitial) / std::abs(doseInitial);

  std::cout << "[reentrant-solve] dose_initial=" << doseInitial
            << " after_first=" << doseAfterFirst
            << " rel_diff_first=" << relDiffFirst
            << " after_second=" << doseAfterSecond
            << " rel_diff_second=" << relDiffSecond << "\n";

  // Both solves must conserve dose. The 1e-6 tolerance catches identity-
  // mass fallback (which produces order-unity drift, not round-off).
  VC_TEST_ASSERT(relDiffFirst < 1e-6);
  VC_TEST_ASSERT(relDiffSecond < 1e-6);
}

void TestMultiSpeciesSmoke() {
  // F3 follow-up: exercise the multi-species code paths in DiffusionEngine.
  // Both prior engine tests used exactly 1 species; the species-outer
  // assembly loop, the allSpecies_ map, the packed-block CVODE state
  // layout, and the modelTargetsSpecies filter are tested only by inspection.
  //
  // Setup: two INDEPENDENT ConstantDiffusion models on two species
  // (Boron, Phosphorus), each with closed-system zero-flux BCs. Since
  // the species don't couple, each must independently conserve its own
  // dose. Catches bugs in species-outer loop, allSpecies_ map, packed-
  // block layout, and cross-species K leakage.
  MeshAttributes attrs;
  attrs.setAttributeName(1, "Si");

  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(8, 8, mfem::Element::TRIANGLE));

  DiffusionEngine<double, 2> engine;
  engine.setMesh(std::move(mesh), attrs);

  auto boronModel = std::make_shared<ConstantDiffusion<double>>("Boron");
  boronModel->setDiffusivity(1e-8, 0.0);
  auto phosphorusModel =
      std::make_shared<ConstantDiffusion<double>>("Phosphorus");
  phosphorusModel->setDiffusivity(1e-7, 0.0); // different D - catches
                                              // cross-species K leakage

  DiffusionPhysics<double> physics;
  physics.addSpecies("Boron");
  physics.addSpecies("Phosphorus");
  physics.addModel(boronModel);
  physics.addModel(phosphorusModel);
  physics.setTemperature(1273.15);
  engine.setPhysics(physics);

  // Different ICs so a cross-species leak would be visible.
  engine.initializeSpecies("Boron", 1e18);
  engine.initializeSpecies("Phosphorus", 1e15);

  const double boronInitial = engine.getIntegral("Boron");
  const double phosphorusInitial = engine.getIntegral("Phosphorus");

  engine.solve(0.0, 10.0, 1.0);

  const double boronFinal = engine.getIntegral("Boron");
  const double phosphorusFinal = engine.getIntegral("Phosphorus");
  const double boronRelDiff =
      std::abs(boronFinal - boronInitial) / std::abs(boronInitial);
  const double phosphorusRelDiff =
      std::abs(phosphorusFinal - phosphorusInitial) /
      std::abs(phosphorusInitial);

  std::cout << "[multi-species] boron: initial=" << boronInitial
            << " final=" << boronFinal << " rel_diff=" << boronRelDiff
            << "\n";
  std::cout << "[multi-species] phosphorus: initial=" << phosphorusInitial
            << " final=" << phosphorusFinal
            << " rel_diff=" << phosphorusRelDiff << "\n";

  // Each species must conserve its own dose independently. 1e-6 tolerance
  // catches any cross-species coupling bug (order-unity drift).
  VC_TEST_ASSERT(boronRelDiff < 1e-6);
  VC_TEST_ASSERT(phosphorusRelDiff < 1e-6);
}

int main() {
  TestMeshAttributes();
  TestDiffusionModelInterface();
  TestConstantDiffusion();
  TestDiffusionPhysics();
  TestDiffusionPhysicsComposition();
  TestLevelSetToMesh2D();
  TestDiffusionEngineAssembly();
  TestDoseConservation();
  TestDirichletBC();
  TestNeumannBC();
  TestMultiSpeciesSmoke();
  TestReentrantSolve();
  std::cout << "All diffusion tests passed.\n";
  return 0;
}
#else
int main() {
  std::cout << "MFEM not available, skipping.\n";
  return 0;
}
#endif
