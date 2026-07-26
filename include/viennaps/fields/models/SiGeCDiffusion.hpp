#pragma once

/// SiGeCDiffusion — carbon I-trapping in SiGe with FEM residual + Ge factor.

#include "CarbonDiffusion.hpp"
#include "../BandgapModel.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class SiGeCDiffusion : public CarbonDiffusion<NumericType> {
public:
  SiGeCDiffusion() { this->setName("SiGeCDiffusion"); }

  void setGeFraction(NumericType x) { x_Ge_ = x; }
  NumericType geFraction() const { return x_Ge_; }

  /// TED enhancement factor: free I reduced by carbon trapping.
  static NumericType tedFactor(NumericType C_I, NumericType C_I_eq,
                               NumericType C_C, NumericType trapStrength) {
    const NumericType freeI =
        C_I / (NumericType(1) + trapStrength * std::max(C_C, NumericType(0)));
    return freeI / std::max(C_I_eq, NumericType(1));
  }

  /// Ge-enhanced carbon diffusivity scale (ni ratio proxy).
  NumericType carbonDiffusivityScale(NumericType T) const {
    BandgapModel<NumericType> bg;
    return bg.niRatioToSi(x_Ge_, T);
  }

private:
  NumericType x_Ge_ = NumericType(0.2);
};

/// Ge-B pairing as FEM residual model (B + Ge ⇌ pair).
template <class NumericType>
class GeBPairingModel : public DiffusionModel<NumericType> {
public:
  GeBPairingModel(std::string boron = "Boron", std::string germanium = "Germanium",
                  std::string pair = "GeBPair")
      : B_(std::move(boron)), Ge_(std::move(germanium)), P_(std::move(pair)) {
    this->setName("GeBPairingModel");
  }

  void setRates(NumericType kf, NumericType kr) {
    kf_ = kf;
    kr_ = kr;
  }

  int numSpecies() const override { return 3; }
  std::vector<std::string> speciesNames() const override {
    return {B_, Ge_, P_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &speciesGF,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    auto itB = allSpecies.find(B_);
    auto itG = allSpecies.find(Ge_);
    auto itP = allSpecies.find(P_);
    if (itB == allSpecies.end() || itG == allSpecies.end() ||
        itP == allSpecies.end() || !itB->second || !itG->second || !itP->second)
      return;
    double scale = 0.0;
    if (&speciesGF == itP->second)
      scale = 1.0;
    else if (&speciesGF == itB->second || &speciesGF == itG->second)
      scale = -1.0;
    else
      return;
    coefs_.clear();
    coefs_.push_back(std::make_unique<PairCoef>(
        *itB->second, *itG->second, *itP->second, static_cast<double>(kf_),
        static_cast<double>(kr_), scale));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coefs_.back()));
  }
  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
#ifdef VIENNAPS_HAS_MFEM
  class PairCoef : public mfem::Coefficient {
  public:
    PairCoef(const mfem::ParGridFunction &B, const mfem::ParGridFunction &G,
             const mfem::ParGridFunction &P, double kf, double kr, double scale)
        : B_(&B), G_(&G), P_(&P), kf_(kf), kr_(kr), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return scale_ * (kf_ * B_->GetValue(T, ip) * G_->GetValue(T, ip) -
                       kr_ * P_->GetValue(T, ip));
    }

  private:
    const mfem::ParGridFunction *B_, *G_, *P_;
    double kf_, kr_, scale_;
  };
  mutable std::vector<std::unique_ptr<PairCoef>> coefs_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
  std::string B_, Ge_, P_;
  NumericType kf_ = NumericType(1e-20);
  NumericType kr_ = NumericType(1e-3);
};

/// Strain-modified constant-D FEM model: D = D0 * exp(-alpha*eps/kT).
template <class NumericType>
class StrainDiffusionModel : public DiffusionModel<NumericType> {
public:
  explicit StrainDiffusionModel(std::string species = "Boron")
      : species_(std::move(species)) {
    this->setName("StrainDiffusionModel(" + species_ + ")");
  }
  void setD0(NumericType D0) { D0_ = D0; }
  void setStrain(NumericType eps) { eps_ = eps; }
  void setAlpha(NumericType a) { alpha_ = a; }

  NumericType getDiffusivity(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return D0_;
    return D0_ * std::exp(-alpha_ * eps_ / (kB * T));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::ParGridFunction *> & /*allSpecies*/,
      const mfem::ParGridFunction * /*temp*/) const override {
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(
        static_cast<double>(getDiffusivity(this->T_)));
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }
  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  NumericType D0_ = NumericType(1e-13);
  NumericType eps_ = NumericType(0);
  NumericType alpha_ = NumericType(0.5);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
