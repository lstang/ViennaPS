#pragma once

/// KmcAtomisticEngine — simple BKL (n-fold way) KMC loop.

#include "KmcEvent.hpp"
#include "KmcLattice.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace viennaps {

class KmcAtomisticEngine {
public:
  explicit KmcAtomisticEngine(unsigned seed = 42) : rng_(seed) {}

  void setLattice(KmcLattice lattice) { lattice_ = std::move(lattice); }
  KmcLattice &lattice() { return lattice_; }
  const KmcLattice &lattice() const { return lattice_; }

  void setParameters(KmcParameters p) { params_ = std::move(p); }
  const KmcParameters &parameters() const { return params_; }

  double time() const { return time_; }
  int steps() const { return steps_; }

  /// One BKL step: build hop events for occupied sites, pick, advance time.
  bool step() {
    std::vector<KmcEvent> events;
    events.reserve(lattice_.size() * 6);
    const double rHop = params_.hopRate();
    for (int k = 0; k < lattice_.nz(); ++k)
      for (int j = 0; j < lattice_.ny(); ++j)
        for (int i = 0; i < lattice_.nx(); ++i) {
          const auto &s = lattice_.at(i, j, k);
          if (!s.occupied)
            continue;
          // 6-neighbor hops on cubic mesh (Si diamond approx for skeleton).
          const int nbr[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                                 {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
          for (auto &d : nbr) {
            const int i1 = i + d[0], j1 = j + d[1], k1 = k + d[2];
            if (i1 < 0 || j1 < 0 || k1 < 0 || i1 >= lattice_.nx() ||
                j1 >= lattice_.ny() || k1 >= lattice_.nz())
              continue;
            if (lattice_.at(i1, j1, k1).occupied)
              continue;
            KmcEvent ev;
            ev.type = KmcEventType::Hop;
            ev.i0 = i;
            ev.j0 = j;
            ev.k0 = k;
            ev.i1 = i1;
            ev.j1 = j1;
            ev.k1 = k1;
            ev.rate = rHop;
            events.push_back(ev);
          }
        }
    if (events.empty())
      return false;

    double Rtot = 0.0;
    for (const auto &e : events)
      Rtot += e.rate;
    if (Rtot <= 0.0)
      return false;

    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double u1 = std::max(U(rng_), 1e-16);
    const double u2 = U(rng_);
    time_ += -std::log(u1) / Rtot;
    double thresh = u2 * Rtot;
    double acc = 0.0;
    const KmcEvent *chosen = &events.back();
    for (const auto &e : events) {
      acc += e.rate;
      if (acc >= thresh) {
        chosen = &e;
        break;
      }
    }
    apply(*chosen);
    ++steps_;
    return true;
  }

  void run(int maxSteps) {
    for (int i = 0; i < maxSteps; ++i) {
      if (!step())
        break;
    }
  }

private:
  void apply(const KmcEvent &e) {
    if (e.type == KmcEventType::Hop) {
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      auto &b = lattice_.at(e.i1, e.j1, e.k1);
      b.species = a.species;
      b.occupied = true;
      a.occupied = false;
      a.species = 0;
    }
  }

  KmcLattice lattice_;
  KmcParameters params_;
  double time_ = 0.0;
  int steps_ = 0;
  std::mt19937 rng_;
};

/// Transfer continuum field → discrete KMC occupations (Poisson sampling).
class KmcAtomize {
public:
  static void atomize(KmcLattice &lat, const std::vector<double> &conc,
                      double volumePerSite, int speciesCode,
                      std::mt19937 &rng) {
    std::size_t idx = 0;
    for (int k = 0; k < lat.nz(); ++k)
      for (int j = 0; j < lat.ny(); ++j)
        for (int i = 0; i < lat.nx(); ++i) {
          if (idx >= conc.size())
            return;
          const double expected = conc[idx++] * volumePerSite;
          std::poisson_distribution<int> pois(
              std::max(0.0, std::min(expected, 20.0)));
          const int n = pois(rng);
          if (n > 0) {
            auto &s = lat.at(i, j, k);
            s.occupied = true;
            s.species = speciesCode;
          }
        }
  }
};

/// Transfer KMC occupations → continuum concentration.
class KmcDeatomize {
public:
  static void deatomize(const KmcLattice &lat, std::vector<double> &conc,
                        double volumePerSite, int speciesCode) {
    conc.assign(lat.size(), 0.0);
    std::size_t idx = 0;
    for (int k = 0; k < lat.nz(); ++k)
      for (int j = 0; j < lat.ny(); ++j)
        for (int i = 0; i < lat.nx(); ++i) {
          const auto &s = lat.at(i, j, k);
          conc[idx++] =
              (s.occupied && s.species == speciesCode) ? (1.0 / volumePerSite)
                                                       : 0.0;
        }
  }
};

struct KmcReport {
  double time = 0;
  int steps = 0;
  int hopCount = 0;
  int occupied = 0;

  static KmcReport fromEngine(const KmcAtomisticEngine &eng) {
    KmcReport r;
    r.time = eng.time();
    r.steps = eng.steps();
    r.occupied = 0;
    for (int k = 0; k < eng.lattice().nz(); ++k)
      for (int j = 0; j < eng.lattice().ny(); ++j)
        for (int i = 0; i < eng.lattice().nx(); ++i)
          if (eng.lattice().at(i, j, k).occupied)
            ++r.occupied;
    return r;
  }
};

/// Continuum ↔ KMC coupling facade for TED validation loops.
class KmcContinuumCoupler {
public:
  /// Run KMC hops then deatomize back to continuum concentration.
  static std::vector<double>
  hopAndDeatomize(KmcLattice lat, KmcParameters params, int speciesCode,
                  double volumePerSite, int steps, unsigned seed = 1) {
    KmcAtomisticEngine eng(seed);
    eng.setLattice(std::move(lat));
    eng.setParameters(params);
    eng.run(steps);
    std::vector<double> conc;
    KmcDeatomize::deatomize(eng.lattice(), conc, volumePerSite, speciesCode);
    return conc;
  }
};

} // namespace viennaps
