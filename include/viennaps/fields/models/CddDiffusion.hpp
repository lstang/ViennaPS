#pragma once

/// CddDiffusion — Classical Dopant Diffusion as composable KernelTerms.
///
/// Composes pair-enhanced dopant diffusion, I/V transport, recombination,
/// and optional cluster bookkeeping. Uses shouldCreateTimeDerivative gate
/// when registered with DiffusionPhysics to avoid double dC/dt.

#include "../DiffusionModel.hpp"
#include "../KernelTerm.hpp"
#include "../KernelTerms.hpp"
#include "PairDiffusion.hpp"
#include "ReactDiffusion.hpp"
#include "Cluster311.hpp"
#include "ImpurityCluster.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class CddDiffusion : public DiffusionModel<NumericType> {
public:
  CddDiffusion() {
    this->setName("CddDiffusion");
    // Default term composition (TED-capable skeleton).
    terms_.push_back(
        std::make_shared<DiffusionTerm>("Interstitial", 1e-10));
    terms_.push_back(std::make_shared<DiffusionTerm>("Vacancy", 1e-12));
    terms_.push_back(std::make_shared<TimeDerivativeTerm>("Boron"));
    terms_.push_back(std::make_shared<TimeDerivativeTerm>("Interstitial"));
    terms_.push_back(std::make_shared<TimeDerivativeTerm>("Vacancy"));
  }

  void addTerm(std::shared_ptr<KernelTerm> term) {
    terms_.push_back(std::move(term));
  }

  const std::vector<std::shared_ptr<KernelTerm>> &terms() const {
    return terms_;
  }

  void setPairDiffusivity(NumericType D_pair) { D_pair_ = D_pair; }
  void setRecombinationRate(NumericType k) { k_recomb_ = k; }

  /// Host-side TED smoke: enhanced D when C_I >> C_I_eq.
  NumericType effectiveDopantDiffusivity(NumericType C_I, NumericType C_I_eq,
                                         NumericType T) const {
    (void)T;
    PairDiffusion<NumericType> pair("Boron", "Interstitial");
    pair.setPairDiffusivity(D_pair_);
    pair.setCIEq(C_I_eq);
    return pair.getDiffusivity(C_I, T);
  }

  int numSpecies() const override {
    // Unique species from terms + defaults.
    return static_cast<int>(speciesNames().size());
  }

  std::vector<std::string> speciesNames() const override {
    std::vector<std::string> names;
    auto add = [&](const std::string &s) {
      if (std::find(names.begin(), names.end(), s) == names.end())
        names.push_back(s);
    };
    add("Boron");
    add("Interstitial");
    add("Vacancy");
    for (const auto &t : terms_) {
      if (t)
        add(t->targetSpecies());
    }
    return names;
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction *temp) const override {
    // Pair-enhanced dopant diffusion when assembling Boron.
    auto itB = allSpecies.find("Boron");
    if (itB != allSpecies.end() && itB->second == &speciesGF) {
      PairDiffusion<NumericType> pair("Boron", "Interstitial");
      pair.setPairDiffusivity(D_pair_);
      pair.setup(*this->attrs_, this->T_);
      pair.assembleStiffness(K, speciesGF, allSpecies, temp);
      return;
    }
    for (const auto &t : terms_) {
      if (!t)
        continue;
      auto it = allSpecies.find(t->targetSpecies());
      if (it != allSpecies.end() && it->second == &speciesGF) {
        t->assembleStiffness(K, allSpecies, temp);
      }
    }
  }

  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction *temp) const override {
    ReactDiffusion<NumericType> react("Interstitial", "Vacancy");
    react.setRecombinationRate(k_recomb_);
    react.assembleReaction(R, speciesGF, allSpecies, temp);
    for (const auto &t : terms_) {
      if (!t)
        continue;
      auto it = allSpecies.find(t->targetSpecies());
      if (it != allSpecies.end() && it->second == &speciesGF) {
        t->assembleResidual(R, allSpecies, temp);
      }
    }
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    // Mass contributed by TimeDerivativeTerm entries when selected by engine
    // gatekeeper. Default unit mass for any species this model owns.
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

private:
  std::vector<std::shared_ptr<KernelTerm>> terms_;
  NumericType D_pair_ = NumericType(1e-13);
  NumericType k_recomb_ = NumericType(1e-15);
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
