#pragma once

/// SiGeDiffusion — Ge interdiffusion (defect-mediated) + B D(x_Ge).
/// D_inter = D_V* · (C_V/C_V*) + D_I* · (C_I/C_I*)  (SProcess-style)
/// Falls back to Arrhenius D0·exp(−Ea/kT) when defect fields are absent.

#include "../BandgapModel.hpp"
#include "../DiffusionModel.hpp"
#include "../PointDefectEquilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class SiGeDiffusion : public DiffusionModel<NumericType> {
public:
  SiGeDiffusion() { this->setName("SiGeDiffusion"); }

  void setGeDiffusivity(NumericType D0, NumericType Ea) {
    D0_Ge_ = D0;
    Ea_Ge_ = Ea;
  }

  void setDefectMediated(NumericType D_Istar, NumericType D_Vstar) {
    D_Istar_ = D_Istar;
    D_Vstar_ = D_Vstar;
    useDefectMediated_ = true;
  }

  void setBoronBaseDiffusivity(NumericType D0) { D0_B_ = D0; }
  void setInterstitialName(std::string s) { I_ = std::move(s); }
  void setVacancyName(std::string s) { V_ = std::move(s); }

  NumericType geDiffusivity(NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return D0_Ge_;
    return D0_Ge_ * std::exp(-Ea_Ge_ / (kB * T));
  }

  /// Defect-mediated interdiffusivity at given C_I, C_V, T.
  NumericType interdiffusivity(NumericType C_I, NumericType C_V,
                               NumericType T) const {
    if (!useDefectMediated_)
      return geDiffusivity(T);
    PointDefectEquilibrium<NumericType> pde;
    const NumericType CIeq = std::max(pde.C_I_eq(T, "Si"), NumericType(1));
    const NumericType CVeq = std::max(pde.C_V_eq(T, "Si"), NumericType(1));
    return D_Istar_ * (C_I / CIeq) + D_Vstar_ * (C_V / CVeq);
  }

  NumericType boronDiffusivity(NumericType x_Ge, NumericType T) const {
    BandgapModel<NumericType> bg;
    const NumericType ratio = bg.niRatioToSi(x_Ge, T);
    return D0_B_ * ratio;
  }

  void applyIntermixStep(std::vector<NumericType> &xGe, NumericType T,
                         NumericType dx, NumericType dt) const {
    if (xGe.size() < 3 || dx <= 0)
      return;
    const NumericType D = geDiffusivity(T);
    const NumericType alpha = D * dt / (dx * dx);
    std::vector<NumericType> next = xGe;
    for (std::size_t i = 1; i + 1 < xGe.size(); ++i) {
      next[i] = xGe[i] + alpha * (xGe[i - 1] - NumericType(2) * xGe[i] +
                                  xGe[i + 1]);
      next[i] = std::min(NumericType(1), std::max(NumericType(0), next[i]));
    }
    xGe.swap(next);
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {"Germanium"};
  }

#ifdef VIENNAPS_HAS_MFEM
  class DefectDCoef : public mfem::Coefficient {
  public:
    DefectDCoef(const SiGeDiffusion *m, double T, const mfem::GridFunction *CI,
                const mfem::GridFunction *CV)
        : m_(m), T_(T), CI_(CI), CV_(CV) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double cI = 0.0, cV = 0.0;
      if (CI_)
        cI = std::max(0.0, CI_->GetValue(T, ip));
      if (CV_)
        cV = std::max(0.0, CV_->GetValue(T, ip));
      if (CI_ || CV_) {
        return static_cast<double>(
            m_->interdiffusivity(static_cast<NumericType>(cI),
                                 static_cast<NumericType>(cV),
                                 static_cast<NumericType>(T_)));
      }
      return static_cast<double>(
          m_->geDiffusivity(static_cast<NumericType>(T_)));
    }

  private:
    const SiGeDiffusion *m_;
    double T_;
    const mfem::GridFunction *CI_;
    const mfem::GridFunction *CV_;
  };

  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    const mfem::GridFunction *CI = nullptr;
    const mfem::GridFunction *CV = nullptr;
    auto itI = allSpecies.find(I_);
    auto itV = allSpecies.find(V_);
    if (itI != allSpecies.end())
      CI = itI->second;
    if (itV != allSpecies.end())
      CV = itV->second;

    if (useDefectMediated_ && (CI || CV)) {
      defectCoef_ = std::make_unique<DefectDCoef>(
          this, static_cast<double>(this->T_), CI, CV);
      K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*defectCoef_));
    } else {
      stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(
          static_cast<double>(geDiffusivity(this->T_)));
      K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
    }
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  NumericType D0_Ge_ = NumericType(1e-3);
  NumericType Ea_Ge_ = NumericType(4.0);
  NumericType D0_B_ = NumericType(1e-13);
  NumericType D_Istar_ = NumericType(1e-13);
  NumericType D_Vstar_ = NumericType(1e-14);
  bool useDefectMediated_ = false;
  std::string I_ = "Interstitial";
  std::string V_ = "Vacancy";
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<DefectDCoef> defectCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

template <class NumericType>
class GeBPairing {
public:
  void setRates(NumericType kf, NumericType kr) {
    kf_ = kf;
    kr_ = kr;
  }

  void applyStep(std::vector<NumericType> &B, std::vector<NumericType> &Ge,
                 std::vector<NumericType> &pair, NumericType dt) const {
    const std::size_t n = std::min({B.size(), Ge.size(), pair.size()});
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType form = kf_ * B[i] * Ge[i];
      const NumericType diss = kr_ * pair[i];
      const NumericType d = (form - diss) * dt;
      pair[i] = std::max(NumericType(0), pair[i] + d);
      B[i] = std::max(NumericType(0), B[i] - d);
      Ge[i] = std::max(NumericType(0), Ge[i] - d);
    }
  }

private:
  NumericType kf_ = NumericType(1e-20);
  NumericType kr_ = NumericType(1e-3);
};

template <class NumericType>
class StrainDiffusionModifier {
public:
  void setAlpha(NumericType a) { alpha_ = a; }

  NumericType modifyD(NumericType D0, NumericType strain,
                      NumericType T) const {
    const NumericType kB = static_cast<NumericType>(8.617333262145e-5);
    if (T <= 0)
      return D0;
    return D0 * std::exp(-alpha_ * strain / (kB * T));
  }

private:
  NumericType alpha_ = NumericType(0.5);
};

} // namespace viennaps
