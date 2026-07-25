#pragma once

/// PairDiffusion — dopant–interstitial pair (TED): D_eff = D_pair * C_I/C_I_eq

#include "../DiffusionModel.hpp"
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

  NumericType getDiffusivity(NumericType C_I, NumericType T) const {
    const NumericType Ceq =
        (Ceq_ > NumericType(0))
            ? Ceq_
            : PointDefectEquilibrium<NumericType>{}.C_I_eq(T, "Si");
    const NumericType ratio =
        C_I / std::max(Ceq, NumericType(1));
    return D_pair_ * std::max(ratio, NumericType(0));
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {dopant_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction & /*speciesGF*/,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction * /*temp*/) const override {
    double C_I_mean = 0.0;
    auto it = allSpecies.find(I_);
    if (it != allSpecies.end() && it->second && it->second->Size() > 0) {
      C_I_mean = it->second->Sum() / it->second->Size();
    }
    const double Deff = static_cast<double>(
        getDiffusivity(static_cast<NumericType>(C_I_mean), this->T_));
    stiffCoef_ = std::make_unique<mfem::ConstantCoefficient>(Deff);
    K.AddDomainIntegrator(new mfem::DiffusionIntegrator(*stiffCoef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string dopant_, I_;
  NumericType D_pair_ = NumericType(1e-13);
  NumericType Ceq_ = NumericType(0);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> stiffCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
