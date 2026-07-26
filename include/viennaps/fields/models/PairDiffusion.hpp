#pragma once

/// PairDiffusion — dopant–interstitial pair (TED): D_eff = D_pair * C_I/C_I_eq
/// Optional SUPG artificial diffusion for sharp C_I fronts.

#include "../DiffusionModel.hpp"
#include "../KernelTerms.hpp"
#include "../PointDefectEquilibrium.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class PairDiffusion : public DiffusionModel<NumericType> {
public:
  PairDiffusion(std::string dopant = "Boron",
                std::string interstitial = "Interstitial")
      : dopant_(std::move(dopant)), I_(std::move(interstitial)) {
    this->setName("PairDiffusion(" + dopant_ + ")");
  }

  void setPairDiffusivity(NumericType D_pair) { D_pair_ = D_pair; }
  void setCIEq(NumericType Ceq) { Ceq_ = Ceq; }
  void setSupg(bool on, NumericType hmin = NumericType(0.1)) {
    useSupg_ = on;
    hmin_ = hmin;
  }
  bool supgEnabled() const { return useSupg_; }

  NumericType getDiffusivity(NumericType C_I, NumericType T) const {
    const NumericType Ceq =
        (Ceq_ > NumericType(0))
            ? Ceq_
            : PointDefectEquilibrium<NumericType>{}.C_I_eq(T, "Si");
    const NumericType ratio = C_I / std::max(Ceq, NumericType(1));
    return D_pair_ * std::max(ratio, NumericType(0));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {dopant_};
  }

#ifdef VIENNAPS_HAS_MFEM
  /// Coefficient: D(x) = D_pair * C_I(x) / C_I_eq evaluated at QPs.
  class PairDCoef : public mfem::Coefficient {
  public:
    PairDCoef(const mfem::ParGridFunction *CI, double Dpair, double Ceq)
        : CI_(CI), Dpair_(Dpair), Ceq_(std::max(Ceq, 1.0)) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      double cI = 0.0;
      if (CI_)
        cI = std::max(0.0, CI_->GetValue(T, ip));
      return Dpair_ * (cI / Ceq_);
    }

  private:
    const mfem::ParGridFunction *CI_;
    double Dpair_, Ceq_;
  };

  void assembleStiffness(
      mfem::ParBilinearForm &K, const mfem::ParGridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::ParGridFunction *> &allSpecies,
      const mfem::ParGridFunction * /*temp*/) const override {
    const mfem::ParGridFunction *CI = nullptr;
    auto it = allSpecies.find(I_);
    if (it != allSpecies.end())
      CI = it->second;

    const double Ceq = static_cast<double>(
        (Ceq_ > NumericType(0))
            ? Ceq_
            : PointDefectEquilibrium<NumericType>{}.C_I_eq(this->T_, "Si"));

    pairCoef_ = std::make_unique<PairDCoef>(CI, static_cast<double>(D_pair_),
                                            Ceq);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*pairCoef_));

    if (useSupg_ && CI && CI->Size() > 0) {
      // Estimate drift from mean |grad C_I| * D_pair / Ceq as advection scale.
      // Skeleton SUPG: isotropic artificial diffusion tau*|v|^2 with
      // |v| ~ D_pair * (mean C_I)/Ceq / h  (order-of-magnitude).
      const double cMean = CI->Sum() / CI->Size();
      const double Deff = static_cast<double>(D_pair_) * cMean / Ceq;
      const double vmag = Deff / std::max(static_cast<double>(hmin_), 1e-12);
      SupgAdvectionTerm supg(dopant_, vmag, 0.0, static_cast<double>(hmin_));
      supg.assembleStiffness(K, allSpecies, nullptr);
    }
  }

  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

  /// Accessors for derived charged-pair FEM assembly.
  const std::string &interstitialName() const { return I_; }
  const std::string &dopantName() const { return dopant_; }
  NumericType pairDiffusivity() const { return D_pair_; }
  NumericType ciEqOverride() const { return Ceq_; }

protected:
  std::string dopant_, I_;
  NumericType D_pair_ = NumericType(1e-13);
  NumericType Ceq_ = NumericType(0);
  bool useSupg_ = false;
  NumericType hmin_ = NumericType(0.1);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<PairDCoef> pairCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
