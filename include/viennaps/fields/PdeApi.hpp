#pragma once

/// PDE API — composable equation terms + results extraction (Phase 10).
/// Engine bridge: PdeEquation::buildModels + applyBCs → DiffusionPhysics.

#include "DiffusionModel.hpp"
#include "DiffusionPhysics.hpp"
#include "models/ConstantDiffusion.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

struct PdeIC {
  std::string species;
  double value = 0.0;
};

struct PdeBC {
  enum class Type { Dirichlet, Neumann, Robin } type = Type::Neumann;
  std::string species;
  std::string boundary = "all";
  double value = 0.0;
};

class PdeTerm {
public:
  virtual ~PdeTerm() = default;
  virtual std::string name() const = 0;
  virtual std::string targetSpecies() const = 0;
};

class DiffusionPdeTerm : public PdeTerm {
public:
  DiffusionPdeTerm(std::string species, double D)
      : species_(std::move(species)), D_(D) {}
  std::string name() const override { return "Diffusion"; }
  std::string targetSpecies() const override { return species_; }
  double D() const { return D_; }

private:
  std::string species_;
  double D_;
};

class ReactionPdeTerm : public PdeTerm {
public:
  ReactionPdeTerm(std::string species, double k)
      : species_(std::move(species)), k_(k) {}
  std::string name() const override { return "Reaction"; }
  std::string targetSpecies() const override { return species_; }
  double k() const { return k_; }

private:
  std::string species_;
  double k_;
};

class PdeEquation {
public:
  void addTerm(std::shared_ptr<PdeTerm> t) { terms_.push_back(std::move(t)); }
  void addIC(PdeIC ic) { ics_.push_back(std::move(ic)); }
  void addBC(PdeBC bc) { bcs_.push_back(std::move(bc)); }

  const std::vector<std::shared_ptr<PdeTerm>> &terms() const { return terms_; }
  const std::vector<PdeIC> &ics() const { return ics_; }
  const std::vector<PdeBC> &bcs() const { return bcs_; }

  std::vector<std::string> species() const {
    std::vector<std::string> s;
    for (const auto &t : terms_) {
      if (!t)
        continue;
      if (std::find(s.begin(), s.end(), t->targetSpecies()) == s.end())
        s.push_back(t->targetSpecies());
    }
    return s;
  }

  /// Build ConstantDiffusion models from DiffusionPdeTerm entries (Ea=0).
  template <class NumericType>
  std::vector<std::shared_ptr<DiffusionModel<NumericType>>>
  buildModels() const {
    std::vector<std::shared_ptr<DiffusionModel<NumericType>>> out;
    for (const auto &t : terms_) {
      auto *diff = dynamic_cast<const DiffusionPdeTerm *>(t.get());
      if (!diff)
        continue;
      auto m = std::make_shared<ConstantDiffusion<NumericType>>(
          diff->targetSpecies());
      m->setDiffusivity(static_cast<NumericType>(diff->D()),
                        NumericType(0));
      out.push_back(std::move(m));
    }
    return out;
  }

  /// Register PdeBC entries on a DiffusionPhysics instance.
  template <class NumericType>
  void applyBCs(DiffusionPhysics<NumericType> &physics) const {
    for (const auto &bc : bcs_) {
      switch (bc.type) {
      case PdeBC::Type::Dirichlet:
        physics.addDirichletBC(bc.species, bc.boundary,
                               static_cast<NumericType>(bc.value));
        break;
      case PdeBC::Type::Neumann:
        physics.addNeumannBC(bc.species, bc.boundary,
                             static_cast<NumericType>(bc.value));
        break;
      case PdeBC::Type::Robin:
        physics.addRobinBC(bc.species, bc.boundary,
                           static_cast<NumericType>(bc.value));
        break;
      }
    }
  }

  /// Convenience: add species + models + BCs to physics from this equation.
  template <class NumericType>
  void applyTo(DiffusionPhysics<NumericType> &physics) const {
    for (const auto &sp : species())
      physics.addSpecies(sp);
    for (auto &m : buildModels<NumericType>())
      physics.addModel(m);
    applyBCs(physics);
  }

private:
  std::vector<std::shared_ptr<PdeTerm>> terms_;
  std::vector<PdeIC> ics_;
  std::vector<PdeBC> bcs_;
};

class ResultsExtractor {
public:
  /// 1D cut along samples: return values at indices.
  static std::vector<double> cut1D(const std::vector<double> &field,
                                   const std::vector<std::size_t> &indices) {
    std::vector<double> out;
    out.reserve(indices.size());
    for (auto i : indices)
      if (i < field.size())
        out.push_back(field[i]);
    return out;
  }

  static double dose(const std::vector<double> &field, double dx) {
    double s = 0.0;
    for (double v : field)
      s += v;
    return s * dx;
  }

  /// First index where field crosses level (linear search).
  static int levelCrossing(const std::vector<double> &field, double level) {
    for (std::size_t i = 1; i < field.size(); ++i) {
      if ((field[i - 1] - level) * (field[i] - level) <= 0.0)
        return static_cast<int>(i);
    }
    return -1;
  }

  /// Sheet resistance proxy: Rs ∝ 1 / ∫ μ(C) C dx  (μ~const → 1/dose).
  static double sheetResistanceProxy(const std::vector<double> &field,
                                     double dx, double mu = 1.0) {
    const double d = dose(field, dx);
    if (d <= 0.0 || mu <= 0.0)
      return 1e300;
    return 1.0 / (mu * d);
  }
};

/// Simple Arrhenius parameter table for calibrated D0/Ea.
struct CalibratedParameters {
  std::map<std::string, std::pair<double, double>> dopantD0Ea;

  void setDopant(const std::string &name, double D0, double Ea) {
    dopantD0Ea[name] = {D0, Ea};
  }

  double diffusivity(const std::string &name, double T) const {
    auto it = dopantD0Ea.find(name);
    if (it == dopantD0Ea.end() || T <= 0)
      return 0.0;
    const double kB = 8.617333262145e-5;
    return it->second.first * std::exp(-it->second.second / (kB * T));
  }

  /// Blend two parameter sets (e.g. inheritance Si → SiGe).
  static CalibratedParameters blend(const CalibratedParameters &a,
                                    const CalibratedParameters &b,
                                    double w) {
    CalibratedParameters out = a;
    for (const auto &kv : b.dopantD0Ea) {
      auto it = out.dopantD0Ea.find(kv.first);
      if (it == out.dopantD0Ea.end())
        out.dopantD0Ea[kv.first] = kv.second;
      else {
        it->second.first =
            (1.0 - w) * it->second.first + w * kv.second.first;
        it->second.second =
            (1.0 - w) * it->second.second + w * kv.second.second;
      }
    }
    return out;
  }
};

/// Least-squares fit of D0 from (T, D) samples with fixed Ea (log-linear).
struct FittingUtilities {
  static double fitD0FixedEa(const std::vector<double> &T,
                             const std::vector<double> &D, double Ea) {
    const double kB = 8.617333262145e-5;
    double num = 0.0, den = 0.0;
    const std::size_t n = std::min(T.size(), D.size());
    for (std::size_t i = 0; i < n; ++i) {
      if (T[i] <= 0 || D[i] <= 0)
        continue;
      const double y = std::log(D[i]) + Ea / (kB * T[i]);
      num += y;
      den += 1.0;
    }
    if (den <= 0)
      return 0.0;
    return std::exp(num / den);
  }
};

/// Build a minimal PdeEquation for constant-D diffusion of one species.
inline PdeEquation makeConstantDiffusionEquation(const std::string &species,
                                                 double D, double C0) {
  PdeEquation eq;
  eq.addTerm(std::make_shared<DiffusionPdeTerm>(species, D));
  eq.addIC({species, C0});
  eq.addBC({PdeBC::Type::Neumann, species, "all", 0.0});
  return eq;
}

} // namespace viennaps

