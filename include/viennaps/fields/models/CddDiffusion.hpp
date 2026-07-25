#pragma once

/// CddDiffusion — Classical Dopant Diffusion as composable KernelTerms +
/// reaction network specification (MOOSE ReactionNetworkPhysicsBase pattern).

#include "../DiffusionModel.hpp"
#include "../KernelTerm.hpp"
#include "../KernelTerms.hpp"
#include "Cluster311.hpp"
#include "ImpurityCluster.hpp"
#include "PairDiffusion.hpp"
#include "ReactDiffusion.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace viennaps {

template <class NumericType>
class CddDiffusion : public DiffusionModel<NumericType> {
public:
  /// One reaction: reactants -> products with rate k (and optional Keq).
  struct ReactionSpec {
    std::vector<std::string> reactants;
    std::vector<std::string> products;
    NumericType rate = NumericType(0);
    NumericType Keq = NumericType(0); // 0 = irreversible
  };

  CddDiffusion() {
    this->setName("CddDiffusion");
    buildDefaultTerms();
  }

  void addTerm(std::shared_ptr<KernelTerm> term) {
    terms_.push_back(std::move(term));
  }

  void addReaction(ReactionSpec r) {
    reactions_.push_back(std::move(r));
  }

  const std::vector<std::shared_ptr<KernelTerm>> &terms() const {
    return terms_;
  }
  const std::vector<ReactionSpec> &reactions() const { return reactions_; }

  void setPairDiffusivity(NumericType D_pair) { D_pair_ = D_pair; }
  void setRecombinationRate(NumericType k) { k_recomb_ = k; }
  void setCIEq(NumericType ceq) { Ceq_ = ceq; }
  void setSupg(bool on) { useSupg_ = on; }
  void enableClusters(bool on) { clusters_ = on; }

  NumericType effectiveDopantDiffusivity(NumericType C_I, NumericType C_I_eq,
                                         NumericType T) const {
    PairDiffusion<NumericType> pair("Boron", "Interstitial");
    pair.setPairDiffusivity(D_pair_);
    pair.setCIEq(C_I_eq);
    return pair.getDiffusivity(C_I, T);
  }

  int numSpecies() const override {
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
    if (clusters_) {
      add("311");
      add("BIC");
    }
    for (const auto &t : terms_)
      if (t)
        add(t->targetSpecies());
    for (const auto &r : reactions_) {
      for (const auto &s : r.reactants)
        add(s);
      for (const auto &s : r.products)
        add(s);
    }
    return names;
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleStiffness(
      mfem::BilinearForm &K, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction *temp) const override {
    auto itB = allSpecies.find("Boron");
    if (itB != allSpecies.end() && itB->second == &speciesGF) {
      PairDiffusion<NumericType> pair("Boron", "Interstitial");
      pair.setPairDiffusivity(D_pair_);
      if (Ceq_ > NumericType(0))
        pair.setCIEq(Ceq_);
      pair.setSupg(useSupg_);
      if (this->attrs_)
        pair.setup(*this->attrs_, this->T_);
      else
        pair.setup(MeshAttributes{}, this->T_);
      pair.assembleStiffness(K, speciesGF, allSpecies, temp);
      return;
    }
    for (const auto &t : terms_) {
      if (!t)
        continue;
      auto it = allSpecies.find(t->targetSpecies());
      if (it != allSpecies.end() && it->second == &speciesGF)
        t->assembleStiffness(K, allSpecies, temp);
    }
  }

  void assembleReaction(
      mfem::LinearForm &R, const mfem::GridFunction &speciesGF,
      const std::map<std::string, mfem::GridFunction *> &allSpecies,
      const mfem::GridFunction *temp) const override {
    ReactDiffusion<NumericType> react("Interstitial", "Vacancy");
    react.setRecombinationRate(k_recomb_);
    react.assembleReaction(R, speciesGF, allSpecies, temp);

    // Host-specified reactions as product coefficients (bimolecular sink/source).
    for (const auto &rx : reactions_) {
      if (rx.reactants.size() == 2) {
        auto itA = allSpecies.find(rx.reactants[0]);
        auto itB = allSpecies.find(rx.reactants[1]);
        if (itA == allSpecies.end() || itB == allSpecies.end() ||
            !itA->second || !itB->second)
          continue;
        // If assembling a reactant: sink -k*A*B; if a product: source +k*A*B.
        double scale = 0.0;
        for (const auto &p : rx.products)
          if (allSpecies.count(p) && allSpecies.at(p) == &speciesGF)
            scale = static_cast<double>(rx.rate);
        for (const auto &r : rx.reactants)
          if (allSpecies.count(r) && allSpecies.at(r) == &speciesGF)
            scale = -static_cast<double>(rx.rate);
        if (scale == 0.0)
          continue;
        prodCoef_ = std::make_unique<ProductCoef>(*itA->second, *itB->second,
                                                  scale);
        R.AddDomainIntegrator(new mfem::DomainLFIntegrator(*prodCoef_));
      }
    }

    for (const auto &t : terms_) {
      if (!t)
        continue;
      auto it = allSpecies.find(t->targetSpecies());
      if (it != allSpecies.end() && it->second == &speciesGF)
        t->assembleResidual(R, allSpecies, temp);
    }
  }

  void assembleMass(mfem::BilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }
#endif

  /// TED host sequence helper: enhanced then relaxed D_eff timeline.
  static std::vector<NumericType>
  tedTimeline(NumericType D_pair, NumericType C_I0, NumericType C_I_eq,
              NumericType k_recomb, int nSteps, NumericType dt) {
    std::vector<NumericType> Deff;
    NumericType CI = C_I0;
    NumericType CV = C_I0; // damage
    for (int i = 0; i < nSteps; ++i) {
      Deff.push_back(D_pair * CI / std::max(C_I_eq, NumericType(1)));
      const NumericType r = k_recomb * CI * CV * dt;
      CI = std::max(NumericType(0), CI - r);
      CV = std::max(NumericType(0), CV - r);
    }
    return Deff;
  }

private:
  void buildDefaultTerms() {
    terms_.clear();
    terms_.push_back(std::make_shared<DiffusionTerm>("Interstitial", 1e-10));
    terms_.push_back(std::make_shared<DiffusionTerm>("Vacancy", 1e-12));
    terms_.push_back(std::make_shared<TimeDerivativeTerm>("Boron"));
    terms_.push_back(std::make_shared<TimeDerivativeTerm>("Interstitial"));
    terms_.push_back(std::make_shared<TimeDerivativeTerm>("Vacancy"));
  }

#ifdef VIENNAPS_HAS_MFEM
  class ProductCoef : public mfem::Coefficient {
  public:
    ProductCoef(const mfem::GridFunction &a, const mfem::GridFunction &b,
                double scale)
        : a_(&a), b_(&b), scale_(scale) {}
    double Eval(mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override {
      return scale_ * a_->GetValue(T, ip) * b_->GetValue(T, ip);
    }

  private:
    const mfem::GridFunction *a_;
    const mfem::GridFunction *b_;
    double scale_;
  };
  mutable std::unique_ptr<ProductCoef> prodCoef_;
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif

  std::vector<std::shared_ptr<KernelTerm>> terms_;
  std::vector<ReactionSpec> reactions_;
  NumericType D_pair_ = NumericType(1e-13);
  NumericType k_recomb_ = NumericType(1e-15);
  NumericType Ceq_ = NumericType(0);
  bool useSupg_ = false;
  bool clusters_ = false;
};

} // namespace viennaps
