#pragma once

/// PDE API — composable equation terms + results extraction (Phase 10 skeleton).

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
};

} // namespace viennaps
