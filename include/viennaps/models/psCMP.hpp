#pragma once

/// CMP — chemical-mechanical planarization (Preston law) process model.
///
/// Removal rate: V = K_p * P * v_rel * s(material) * f_pattern(h), where
/// K_p is the Preston coefficient, P the downforce, v_rel the pad-wafer
/// relative velocity, s(material) a per-material selectivity (0 for
/// non-polish materials and hard stops), and
///   f_pattern(h) = clamp(1 + alpha * (h - h_ref) / L_p, 0.1, 2.0)
/// the pattern-density term: protrusions (h > h_ref) polish faster,
/// depressions slower (planarization). GAP_ANALYSIS §4.4 — this is polish
/// physics, distinct from the geometric plane cut in psPlanarize.hpp.

#include "../materials/psMaterialMap.hpp"
#include "../process/psProcessModel.hpp"

#include <vcVectorType.hpp>

#include <algorithm>
#include <vector>

namespace viennaps {

using namespace viennacore;

namespace impl {

template <class NumericType, int D>
class CmpVelocityField : public VelocityField<NumericType, D> {
  const NumericType prestonRate_;
  const NumericType alpha_;
  const NumericType planarizationLength_;
  const NumericType refHeight_;
  const std::vector<std::pair<Material, NumericType>> &materials_;
  const std::vector<Material> &hardStop_;

public:
  CmpVelocityField(
      NumericType prestonRate, NumericType alpha,
      NumericType planarizationLength, NumericType refHeight,
      const std::vector<std::pair<Material, NumericType>> &materials,
      const std::vector<Material> &hardStop)
      : prestonRate_(prestonRate), alpha_(alpha),
        planarizationLength_(planarizationLength), refHeight_(refHeight),
        materials_(materials), hardStop_(hardStop) {}

  NumericType getScalarVelocity(const Vec3D<NumericType> &coordinate,
                                int material, const Vec3D<NumericType> &nv,
                                unsigned long pointID) override {
    // Hard-stop materials never polish.
    for (const auto &m : hardStop_) {
      if (MaterialMap::isMaterial(material, m))
        return NumericType(0);
    }
    // Only configured polish materials are removed; selectivity scales
    // the base rate per material (e.g. nitride vs oxide).
    NumericType selectivity = NumericType(0);
    for (const auto &pm : materials_) {
      if (MaterialMap::isMaterial(material, pm.first)) {
        selectivity = pm.second;
        break;
      }
    }
    if (selectivity <= NumericType(0))
      return NumericType(0);

    // Preston base rate with pattern-density modulation. Protrusions
    // (h > h_ref) polish faster, depressions slower; clamped so the
    // field stays positive and bounded.
    const NumericType h = coordinate[D - 1];
    NumericType patternFactor = NumericType(1);
    if (planarizationLength_ > NumericType(0)) {
      patternFactor += alpha_ * (h - refHeight_) / planarizationLength_;
      patternFactor =
          std::clamp(patternFactor, NumericType(0.1), NumericType(2.0));
    }
    return prestonRate_ * selectivity * patternFactor;
  }
};

} // namespace impl

// Model for a chemical-mechanical planarization (CMP) process.
template <typename NumericType, int D>
class CMP : public ProcessModelCPU<NumericType, D> {
public:
  CMP() { initialize(); }

  // --- Preston parameters -------------------------------------------------
  void setPressure(NumericType p) {
    pressure_ = p;
    refreshRate();
  }
  void setRelativeVelocity(NumericType v) {
    velocity_ = v;
    refreshRate();
  }
  void setPrestonCoefficient(NumericType k) {
    prestonK_ = k;
    refreshRate();
  }

  // --- Pattern-density parameters -----------------------------------------
  void setPatternDensity(NumericType alpha, NumericType planarizationLength) {
    alpha_ = alpha;
    planarizationLength_ = planarizationLength;
    refreshField();
  }
  void setReferenceHeight(NumericType hRef) {
    refHeight_ = hRef;
    refreshField();
  }

  // --- Materials ----------------------------------------------------------
  void addPolishingMaterial(Material m, NumericType selectivity = 1.) {
    polishMaterials_.emplace_back(m, selectivity);
    refreshField();
  }
  void setHardStopMaterials(std::vector<Material> stops) {
    hardStop_ = std::move(stops);
    refreshField();
  }

private:
  void initialize() {
    // default surface model (pure advection, no surface reactions)
    auto surfModel = SmartPointer<SurfaceModel<NumericType>>::New();

    this->setSurfaceModel(surfModel);
    this->setProcessName("CMP");
    refreshRate();
    refreshField();

    processMetaData["PrestonCoefficient"] = {prestonK_};
    processMetaData["Pressure"] = {pressure_};
    processMetaData["RelativeVelocity"] = {velocity_};
    processMetaData["ReferenceHeight"] = {refHeight_};
    for (const auto &m : polishMaterials_) {
      processMetaData[MaterialMap::toString(m.first) + " Selectivity"] =
          std::vector<double>{m.second};
    }
  }

  void refreshRate() { prestonRate_ = prestonK_ * pressure_ * velocity_; }

  void refreshField() {
    auto velField = SmartPointer<impl::CmpVelocityField<NumericType, D>>::New(
        prestonRate_, alpha_, planarizationLength_, refHeight_,
        polishMaterials_, hardStop_);
    this->setVelocityField(velField);
  }

  NumericType pressure_ = 5.;
  NumericType velocity_ = 1.;
  NumericType prestonK_ = 1e-4;
  NumericType alpha_ = 0.5;
  NumericType planarizationLength_ = 1.;
  NumericType refHeight_ = 0.;
  NumericType prestonRate_ = 0.;
  std::vector<std::pair<Material, NumericType>> polishMaterials_;
  std::vector<Material> hardStop_;
  using ProcessModelCPU<NumericType, D>::processMetaData;
};

PS_PRECOMPILE_PRECISION_DIMENSION(CMP)

} // namespace viennaps
