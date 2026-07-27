#pragma once

/// PDE API — composable equation terms + results extraction (Phase 10).
/// Engine bridge: PdeEquation::buildModels + applyBCs → DiffusionPhysics.

#include "DiffusionModel.hpp"
#include "DiffusionPhysics.hpp"
#include "models/ConstantDiffusion.hpp"
#include "models/LinearReactionDiffusion.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

struct PdeIC {
  std::string species;
  double value = 0.0;
};

struct PdeBC {
  enum class Type {
    Dirichlet,
    Neumann,
    Robin,
    Flux,
    Segregation
  } type = Type::Neumann;
  std::string species;
  std::string boundary = "all";
  double value = 0.0;
  /// Segregation: partner species on the other side of the interface.
  std::string partnerSpecies;
  double segregationM = 1.0;
};

struct FluxBC {
  std::string species;
  std::string boundary = "all";
  double flux = 0.0;
};

struct SegregationBC {
  std::string speciesA;
  std::string speciesB;
  double m = 1.0;
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

  /// Build models from DiffusionPdeTerm and ReactionPdeTerm.
  template <class NumericType>
  std::vector<std::shared_ptr<DiffusionModel<NumericType>>>
  buildModels() const {
    std::vector<std::shared_ptr<DiffusionModel<NumericType>>> out;
    for (const auto &t : terms_) {
      if (auto *diff = dynamic_cast<const DiffusionPdeTerm *>(t.get())) {
        auto m = std::make_shared<ConstantDiffusion<NumericType>>(
            diff->targetSpecies());
        m->setDiffusivity(static_cast<NumericType>(diff->D()), NumericType(0));
        out.push_back(std::move(m));
      } else if (auto *rx = dynamic_cast<const ReactionPdeTerm *>(t.get())) {
        out.push_back(std::make_shared<LinearReactionDiffusion<NumericType>>(
            rx->targetSpecies(), static_cast<NumericType>(rx->k())));
      }
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
      case PdeBC::Type::Flux:
        physics.addNeumannBC(bc.species, bc.boundary,
                             static_cast<NumericType>(bc.value));
        break;
      case PdeBC::Type::Robin:
        physics.addRobinBC(bc.species, bc.boundary,
                           static_cast<NumericType>(bc.value));
        break;
      case PdeBC::Type::Segregation:
        // Segregation is applied as operator-split by engine when models
        // register; BC value stores m for host tools.
        (void)bc.segregationM;
        break;
      }
    }
  }

  /// Return IC map species → value for engine.initializeSpecies.
  std::map<std::string, double> initialConditions() const {
    std::map<std::string, double> m;
    for (const auto &ic : ics_)
      m[ic.species] = ic.value;
    return m;
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

  /// Convenience overload: also apply ICs to the engine when provided.
  /// Solves the prior review's M1 finding ("applyTo does not call applyICs;
  /// user must remember separately"). The 1-arg overload above remains for
  /// backward compatibility.
  template <class NumericType, class Engine>
  void applyTo(DiffusionPhysics<NumericType> &physics, Engine *engine) const {
    applyTo(physics);
    if (engine)
      applyICs(*engine);
  }

  /// Apply ICs to an engine (initializeSpecies for each stored PdeIC).
  template <class Engine> void applyICs(Engine &engine) const {
    for (const auto &ic : ics_)
      engine.initializeSpecies(ic.species, ic.value);
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

#ifdef VIENNAPS_HAS_MFEM
  /// Sample GridFunction along a line from (x0,y0) to (x1,y1) with n points.
  static std::vector<double> cut1D(const mfem::GridFunction &gf, double x0,
                                   double y0, double x1, double y1, int n) {
    std::vector<double> out;
    if (!gf.FESpace() || !gf.FESpace()->GetMesh() || n <= 0)
      return out;
    mfem::Mesh *mesh = gf.FESpace()->GetMesh();
    out.resize(static_cast<std::size_t>(n), 0.0);
    mfem::DenseMatrix pts(mesh->SpaceDimension(), n);
    for (int i = 0; i < n; ++i) {
      const double t = (n == 1) ? 0.0 : static_cast<double>(i) / (n - 1);
      pts(0, i) = x0 + t * (x1 - x0);
      if (mesh->SpaceDimension() > 1)
        pts(1, i) = y0 + t * (y1 - y0);
      for (int d = 2; d < mesh->SpaceDimension(); ++d)
        pts(d, i) = 0.0;
    }
    mfem::Array<int> el(n);
    mfem::Array<mfem::IntegrationPoint> ips(n);
    mesh->FindPoints(pts, el, ips, /*warn=*/false);
    for (int i = 0; i < n; ++i) {
      if (el[i] < 0)
        continue;
      out[static_cast<std::size_t>(i)] =
          gf.GetValue(*mesh->GetElementTransformation(el[i]), ips[i]);
    }
    return out;
  }

  static double dose(const mfem::GridFunction &gf) {
    mfem::ConstantCoefficient one(1.0);
    mfem::LinearForm mass(const_cast<mfem::FiniteElementSpace *>(gf.FESpace()));
    mass.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    mass.Assemble();
    return gf * mass;
  }
#endif

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

  /// Write depth profile CSV: columns depth,value.
  static bool writeCSV(const std::string &path,
                       const std::vector<double> &depth,
                       const std::vector<double> &value) {
    std::ofstream os(path);
    if (!os)
      return false;
    os << "depth,value\n";
    const std::size_t n = std::min(depth.size(), value.size());
    for (std::size_t i = 0; i < n; ++i)
      os << depth[i] << "," << value[i] << "\n";
    return true;
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
                                    const CalibratedParameters &b, double w) {
    CalibratedParameters out = a;
    for (const auto &kv : b.dopantD0Ea) {
      auto it = out.dopantD0Ea.find(kv.first);
      if (it == out.dopantD0Ea.end())
        out.dopantD0Ea[kv.first] = kv.second;
      else {
        it->second.first = (1.0 - w) * it->second.first + w * kv.second.first;
        it->second.second =
            (1.0 - w) * it->second.second + w * kv.second.second;
      }
    }
    return out;
  }
};

/// Least-squares fit utilities.
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

  /// Linear fit y = a + b x → returns {a, b}.
  static std::pair<double, double> fitLine(const std::vector<double> &x,
                                           const std::vector<double> &y) {
    const std::size_t n = std::min(x.size(), y.size());
    if (n < 2)
      return {0.0, 0.0};
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < n; ++i) {
      sx += x[i];
      sy += y[i];
      sxx += x[i] * x[i];
      sxy += x[i] * y[i];
    }
    const double nd = static_cast<double>(n);
    const double den = nd * sxx - sx * sx;
    if (std::abs(den) < 1e-30)
      return {sy / nd, 0.0};
    const double b = (nd * sxy - sx * sy) / den;
    const double a = (sy - b * sx) / nd;
    return {a, b};
  }

  /// Fit Pearson-IV moments proxy: returns {Rp, deltaRp, skew} from depth C.
  static std::array<double, 3> fitPearson(const std::vector<double> &z,
                                          const std::vector<double> &C) {
    const std::size_t n = std::min(z.size(), C.size());
    double dose = 0.0;
    for (std::size_t i = 0; i < n; ++i)
      dose += std::max(0.0, C[i]);
    if (dose <= 0.0 || n < 2)
      return {0.0, 0.0, 0.0};
    double m1 = 0.0;
    for (std::size_t i = 0; i < n; ++i)
      m1 += z[i] * std::max(0.0, C[i]);
    m1 /= dose;
    double m2 = 0.0, m3 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double dz = z[i] - m1;
      const double w = std::max(0.0, C[i]);
      m2 += w * dz * dz;
      m3 += w * dz * dz * dz;
    }
    m2 /= dose;
    m3 /= dose;
    const double sig = std::sqrt(std::max(m2, 0.0));
    const double skew = (sig > 0) ? m3 / (sig * sig * sig) : 0.0;
    return {m1, sig, skew};
  }

  /// Pearson with floor: clamp C to floor before fit.
  static std::array<double, 3> fitPearsonFloor(const std::vector<double> &z,
                                               const std::vector<double> &C,
                                               double floor) {
    std::vector<double> Cf = C;
    for (auto &v : Cf)
      v = std::max(v, floor);
    return fitPearson(z, Cf);
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
