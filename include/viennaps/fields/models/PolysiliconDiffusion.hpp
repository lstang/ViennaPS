#pragma once

/// PolysiliconDiffusion — isotropic (D_bulk + D_gb*f_gb) and anisotropic
/// (attribute-dependent D on dual mesh interior/GB).

#include "../DiffusionModel.hpp"
#include "../GrainBoundaryMesh.hpp"
#include "../GrainModel.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class PolysiliconDiffusion : public DiffusionModel<NumericType> {
public:
  enum class Mode { Isotropic, Anisotropic };

  explicit PolysiliconDiffusion(std::string species = "Boron")
      : species_(std::move(species)) {
    this->setName("PolysiliconDiffusion(" + species_ + ")");
  }

  void setMode(Mode m) { mode_ = m; }
  Mode mode() const { return mode_; }

  void setDiffusivities(NumericType D_bulk, NumericType D_gb) {
    D_bulk_ = D_bulk;
    D_gb_ = D_gb;
  }

  void setGrainModel(GrainModel<NumericType> g) { grains_ = std::move(g); }
  GrainModel<NumericType> &grains() { return grains_; }
  const GrainModel<NumericType> &grains() const { return grains_; }

  void setBoundaryWidth(NumericType w) { boundaryWidth_ = w; }
  void setSegregationCoefficient(NumericType m) { m_seg_ = m; }
  void setAttributes(int interiorAttr, int boundaryAttr) {
    interiorAttr_ = interiorAttr;
    boundaryAttr_ = boundaryAttr;
  }

  NumericType getIsotropicDiffusivity() const {
    const NumericType f = grains_.grainBoundaryFraction(boundaryWidth_);
    return D_bulk_ + D_gb_ * f;
  }

  NumericType getInteriorDiffusivity() const { return D_bulk_; }
  NumericType getBoundaryDiffusivity() const { return D_gb_; }
  NumericType segregationCoefficient() const { return m_seg_; }

  void applySegregationStep(NumericType &C_int, NumericType &C_gb,
                            NumericType k0, NumericType dt) const {
    const NumericType kf = m_seg_ * k0;
    const NumericType kb = k0;
    const NumericType rate = kf * C_int - kb * C_gb;
    C_int -= rate * dt;
    C_gb += rate * dt;
  }

  /// Dual-species names when GB segregation is active as FEM residual.
  void enableGbSegregationSpecies(std::string gbSpecies = "Boron_GB") {
    gbSpecies_ = std::move(gbSpecies);
  }

  int numSpecies() const override {
    return gbSpecies_.empty() ? 1 : 2;
  }
  std::vector<std::string> speciesNames() const override {
    if (gbSpecies_.empty())
      return {species_};
    return {species_, gbSpecies_};
  }

#ifdef VIENNAPS_HAS_MFEM
  /// PWConst coefficient: D_bulk on interior attr, D_gb on boundary attr.
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    if (mode_ == Mode::Isotropic) {
      stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(
          static_cast<double>(getIsotropicDiffusivity()));
      K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
      return;
    }
    // Anisotropic: PWConst by element attribute (1-based in MFEM).
    // Build maxAttr from attrs_ if available.
    int maxAttr = std::max(interiorAttr_, boundaryAttr_);
    pwValues_.assign(static_cast<std::size_t>(maxAttr),
                     static_cast<double>(D_bulk_));
    if (interiorAttr_ >= 1 &&
        interiorAttr_ <= static_cast<int>(pwValues_.size()))
      pwValues_[static_cast<std::size_t>(interiorAttr_ - 1)] =
          static_cast<double>(D_bulk_);
    if (boundaryAttr_ >= 1 &&
        boundaryAttr_ <= static_cast<int>(pwValues_.size()))
      pwValues_[static_cast<std::size_t>(boundaryAttr_ - 1)] =
          static_cast<double>(D_gb_);
    pwCoef_ = std::make_unique<mfem::PWConstCoefficient>(
        mfem::Vector(pwValues_.data(), static_cast<int>(pwValues_.size())));
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*pwCoef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }

  /// FEM residual for GI↔GB exchange when dual species are registered.
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    if (gbSpecies_.empty())
      return;
    auto itI = allSpecies.find(species_);
    auto itG = allSpecies.find(gbSpecies_);
    if (itI == allSpecies.end() || itG == allSpecies.end() || !itI->second ||
        !itG->second)
      return;
    double scale = 0.0;
    if (&speciesGF == itG->second)
      scale = 1.0;
    else if (&speciesGF == itI->second)
      scale = -1.0;
    else
      return;
    const double kf = static_cast<double>(m_seg_ * k0_fem_);
    const double kb = static_cast<double>(k0_fem_);
    segCoefs_.clear();
    segCoefs_.push_back(std::make_unique<SegCoef>(
        *itI->second, *itG->second, kf, kb, scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*segCoefs_.back()));
  }
#endif

  void setFemSegregationRate(NumericType k0) { k0_fem_ = k0; }

private:
  std::string species_;
  std::string gbSpecies_;
  NumericType k0_fem_ = NumericType(1e-3);
  Mode mode_ = Mode::Isotropic;
  NumericType D_bulk_ = NumericType(1e-14);
  NumericType D_gb_ = NumericType(1e-10);
  NumericType boundaryWidth_ = NumericType(1e-7);
  NumericType m_seg_ = NumericType(10);
  int interiorAttr_ = 1;
  int boundaryAttr_ = 2;
  GrainModel<NumericType> grains_;
#ifdef VIENNAPS_HAS_MFEM
  class SegCoef : public mfem::Coefficient {
  public:
    SegCoef(const mfem::GridFunction &gi, const mfem::GridFunction &gb,
            double kf, double kb, double scale)
        : gi_(&gi), gb_(&gb), kf_(kf), kb_(kb), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return scale_ * (kf_ * gi_->GetValue(T, ip) - kb_ * gb_->GetValue(T, ip));
    }

  private:
    const mfem::GridFunction *gi_;
    const mfem::GridFunction *gb_;
    double kf_, kb_, scale_;
  };
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
  mutable std::unique_ptr<mfem::PWConstCoefficient> pwCoef_;
  mutable std::vector<double> pwValues_;
  mutable std::vector<std::unique_ptr<SegCoef>> segCoefs_;
#endif
};

template <class NumericType>
class PolyOxideBreakup {
public:
  void setRates(NumericType k_dissolve, NumericType k_epi) {
    k_dissolve_ = k_dissolve;
    k_epi_ = k_epi;
  }
  void setOxideThickness(NumericType t) {
    t_ox_ = std::max(t, NumericType(0));
  }
  void setEpiFraction(NumericType f) {
    f_epi_ = std::min(NumericType(1), std::max(NumericType(0), f));
  }
  NumericType oxideThickness() const { return t_ox_; }
  NumericType epiFraction() const { return f_epi_; }
  void advance(NumericType T, NumericType dt) {
    (void)T;
    t_ox_ = std::max(NumericType(0), t_ox_ - k_dissolve_ * dt);
    if (t_ox_ <= NumericType(0))
      f_epi_ = std::min(NumericType(1), f_epi_ + k_epi_ * dt);
  }

private:
  NumericType t_ox_ = NumericType(1e-7);
  NumericType f_epi_ = NumericType(0);
  NumericType k_dissolve_ = NumericType(1e-8);
  NumericType k_epi_ = NumericType(0.1);
};

template <class NumericType>
class PolyOxidationRate {
public:
  void setBaseRate(NumericType v0) { v0_ = v0; }
  void setGrainSensitivity(NumericType alpha) { alpha_ = alpha; }
  NumericType rate(NumericType R) const {
    const NumericType Rref = NumericType(1e-5);
    return v0_ * (NumericType(1) +
                  alpha_ * Rref / std::max(R, NumericType(1e-12)));
  }

private:
  NumericType v0_ = NumericType(1e-9);
  NumericType alpha_ = NumericType(1);
};

} // namespace viennaps
