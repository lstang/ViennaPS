#pragma once

/// LinearReactionDiffusion — FEM model for PDE API ReactionPdeTerm.
/// Residual R = −k · C (decay) or +k · C when k < 0 (source).

#include "../DiffusionModel.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class LinearReactionDiffusion : public DiffusionModel<NumericType> {
public:
  LinearReactionDiffusion(std::string species, NumericType k)
      : species_(std::move(species)), k_(k) {
    this->setName("LinearReactionDiffusion(" + species_ + ")");
  }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> & /*allSpecies*/,
      const mfem::GridFunction * /*temp*/) const override {
    coef_ = std::make_unique<ScaledGF>(speciesGF, -static_cast<double>(k_));
    R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coef_));
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::string species_;
  NumericType k_ = NumericType(0);
#ifdef VIENNAPS_HAS_MFEM
  class ScaledGF : public mfem::Coefficient {
  public:
    ScaledGF(const mfem::GridFunction &gf, double s) : gf_(&gf), s_(s) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return s_ * gf_->GetValue(T, ip);
    }

  private:
    const mfem::GridFunction *gf_;
    double s_;
  };
  mutable std::unique_ptr<ScaledGF> coef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
