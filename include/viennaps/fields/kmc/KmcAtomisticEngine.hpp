#pragma once

/// BKL kinetic Monte Carlo atomistic engine (SProcess Ch. 5 partial;
/// Phase 7 skeleton + Phase 8 LKMC epitaxy production rates).

#include "KmcEvent.hpp"
#include "KmcLattice.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace viennaps {

template <class Real> class FenwickTree {
public:
  explicit FenwickTree(std::size_t n = 0) { init(n); }

  void init(std::size_t n) {
    n_ = n;
    tree_.assign(n_ + 1, Real(0));
  }

  std::size_t size() const { return n_; }

  void add(std::size_t idx, Real val) {
    for (std::size_t i = idx + 1; i <= n_; i += (i & -i))
      tree_[i] += val;
  }

  Real total() const {
    Real sum = 0;
    for (std::size_t i = n_; i > 0; i -= (i & -i))
      sum += tree_[i];
    return sum;
  }

  /// 0-indexed lower bound: returns smallest 0-based index whose prefix sum >= val.
  std::size_t lower_bound(Real val) const {
    if (n_ == 0)
      return 0;
    std::size_t idx = 0;
    std::size_t mask = 1;
    while (mask <= n_)
      mask <<= 1;
    mask >>= 1;

    for (; mask > 0; mask >>= 1) {
      const std::size_t nextIdx = idx + mask;
      if (nextIdx <= n_ && val > tree_[nextIdx]) {
        idx = nextIdx;
        val -= tree_[idx];
      }
    }
    return std::min(idx, n_ - 1);
  }

private:
  std::size_t n_ = 0;
  std::vector<Real> tree_;
};

template <class NumericType> class KmcAtomisticEngine {
public:
  KmcAtomisticEngine() = default;
  explicit KmcAtomisticEngine(unsigned seed) : rng_(seed) {}

  void setLattice(KmcLattice lat) {
    lattice_ = std::move(lat);
    fenwickInitialized_ = false;
  }
  const KmcLattice &lattice() const { return lattice_; }
  KmcLattice &lattice() { return lattice_; }

  void setParameters(KmcParameters p) {
    params_ = p;
    fenwickInitialized_ = false;
  }
  const KmcParameters &parameters() const { return params_; }

  /// Ge mole fraction for Deposit events: the deposited species is KmcGe
  /// with probability xGe, else growthSpecies_; the deposit rate is scaled
  /// by the Ge growth factor (1 - 0.3*xGe) like KmcEpitaxyModel.
  void setGeFraction(double x) { xGe_ = std::clamp(x, 0.0, 1.0); }
  double geFraction() const { return xGe_; }
  /// Deposit candidates need at least this many occupied neighbors.
  void setMinCoordination(int c) { minCoord_ = std::max(0, c); }
  /// Require deposit sites to pass KmcVisibility::isVisible (z-buffer).
  void setVisibilityEnabled(bool on) { visibility_ = on; }
  bool diamondNeighbors() const { return diamond_; }

  void setDiamondNeighbors(bool on) { diamond_ = on; }
  void setRecombinationEnabled(bool on) { recomb_ = on; }
  void setClusteringEnabled(bool on) { cluster_ = on; }
  void setEpitaxyEnabled(bool on) { epitaxy_ = on; }
  void setGrowthSpecies(int code) { growthSpecies_ = code; }

  double time() const { return time_; }
  int steps() const { return steps_; }
  int recombCount() const { return recombCount_; }
  int clusterCount() const { return clusterCount_; }
  int dissocCount() const { return dissocCount_; }
  int depositCount() const { return depositCount_; }
  int desorbCount() const { return desorbCount_; }
  int twinCount() const { return twinCount_; }

  /// Build per-site event lists on first call; thereafter pick via Fenwick
  /// tree and incrementally update affected neighborhoods.
  bool step() {
    if (!fenwickInitialized_) {
      buildAllSiteEvents();
      fenwickInitialized_ = true;
    }

    const double totalRate = fenwick_.total();
    if (totalRate <= 0.0)
      return false;

    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double u1 = std::max(U(rng_), 1e-16);
    const double u2 = U(rng_);
    time_ += -std::log(u1) / totalRate;
    const double thresh = u2 * totalRate;
    std::size_t siteIdx =
        std::min(fenwick_.lower_bound(thresh), fenwick_.size() - 1);

    const auto &siteEvents = eventsBySite_[siteIdx];
    double r = U(rng_) * siteTotalRates_[siteIdx];
    const KmcEvent *chosen = nullptr;
    for (const auto &ev : siteEvents) {
      r -= ev.rate;
      if (r <= 0.0) {
        chosen = &ev;
        break;
      }
    }
    if (!chosen && !siteEvents.empty() && siteTotalRates_[siteIdx] > 0.0)
      chosen = &siteEvents.back();

    if (chosen) {
      apply(*chosen);
      rebuildAffectedSites(*chosen);
    }
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
  template <class Fn>
  void forEachNeighbor(int i, int j, int k, Fn &&fn) const {
    if (diamond_) {
      // Diamond A-B only: neighbors are opposite sublattice via body diagonals
      // (2-FCC approximation on a cubic grid).
      static const int diag[4][3] = {
          {1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
      for (const auto &d : diag) {
        const int i1 = i + d[0], j1 = j + d[1], k1 = k + d[2];
        if (inBounds(i1, j1, k1))
          fn(i1, j1, k1);
      }
      return;
    }
    static const int cubic[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                                    {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (const auto &d : cubic) {
      const int i1 = i + d[0], j1 = j + d[1], k1 = k + d[2];
      if (inBounds(i1, j1, k1))
        fn(i1, j1, k1);
    }
  }

  bool inBounds(int i, int j, int k) const {
    return i >= 0 && j >= 0 && k >= 0 && i < lattice_.nx() &&
           j < lattice_.ny() && k < lattice_.nz();
  }

  std::size_t linearIndex(int i, int j, int k) const {
    return static_cast<std::size_t>(
        (k * lattice_.ny() + j) * lattice_.nx() + i);
  }

  int coordinatedNeighbors(int i, int j, int k) const {
    int c = 0;
    forEachNeighbor(i, j, k, [&](int i1, int j1, int k1) {
      if (lattice_.at(i1, j1, k1).occupied)
        ++c;
    });
    return c;
  }

  void buildAllSiteEvents() {
    const std::size_t n = static_cast<std::size_t>(lattice_.size());
    eventsBySite_.assign(n, std::vector<KmcEvent>());
    siteTotalRates_.assign(n, 0.0);
    fenwick_.init(n);
    for (std::size_t idx = 0; idx < n; ++idx)
      buildSiteEvents(idx);
  }

  void buildSiteEvents(std::size_t siteIdx) {
    const int nx = lattice_.nx(), ny = lattice_.ny(), nz = lattice_.nz();
    const int k = static_cast<int>(siteIdx / (ny * nx));
    const int j = static_cast<int>((siteIdx % (ny * nx)) / nx);
    const int i = static_cast<int>(siteIdx % nx);
    const auto &s = lattice_.at(i, j, k);

    double totalRate = 0.0;
    auto &evList = eventsBySite_[siteIdx];
    evList.clear();

    const double rHop = params_.hopRate();
    const double rRec = params_.recombRate();
    const double rCl = params_.clusterRate();
    const double rDi = params_.dissocRate();

    if (s.occupied) {
      if (cluster_ && s.species == KmcCluster311) {
        forEachNeighbor(i, j, k, [&](int i1, int j1, int k1) {
          if (lattice_.at(i1, j1, k1).occupied)
            return;
          KmcEvent ev;
          ev.type = KmcEventType::Dissociate;
          ev.i0 = i;  ev.j0 = j;  ev.k0 = k;
          ev.i1 = i1; ev.j1 = j1; ev.k1 = k1;
          ev.rate = rDi;
          evList.push_back(ev);
          totalRate += rDi;
        });
      }

      forEachNeighbor(i, j, k, [&](int i1, int j1, int k1) {
        const auto &n = lattice_.at(i1, j1, k1);
        if (!n.occupied) {
          if (s.species == KmcInterstitial || s.species == KmcVacancy) {
            KmcEvent ev;
            ev.type = KmcEventType::Hop;
            ev.i0 = i;  ev.j0 = j;  ev.k0 = k;
            ev.i1 = i1; ev.j1 = j1; ev.k1 = k1;
            ev.rate = rHop;
            evList.push_back(ev);
            totalRate += rHop;
          }
          return;
        }
        if (recomb_ && s.species == KmcInterstitial &&
            n.species == KmcVacancy) {
          KmcEvent ev;
          ev.type = KmcEventType::Recombine;
          ev.i0 = i;  ev.j0 = j;  ev.k0 = k;
          ev.i1 = i1; ev.j1 = j1; ev.k1 = k1;
          ev.rate = rRec;
          evList.push_back(ev);
          totalRate += rRec;
        }
        if (cluster_ && s.species == KmcInterstitial &&
            n.species == KmcInterstitial) {
          const std::size_t a = siteIdx;
          const std::size_t b = linearIndex(i1, j1, k1);
          if (a < b) {
            KmcEvent ev;
            ev.type = KmcEventType::Cluster;
            ev.i0 = i;  ev.j0 = j;  ev.k0 = k;
            ev.i1 = i1; ev.j1 = j1; ev.k1 = k1;
            ev.rate = rCl;
            evList.push_back(ev);
            totalRate += rCl;
          }
        }
      });
    }

    if (epitaxy_) {
      // Generalized surface events. Deposit: any empty site with >= minCoord
      // occupied neighbors (and, optionally, clear line of sight), rate
      // scaled by coordination fraction c/c_max and the Ge growth factor.
      // On a perfect crystal this reduces to the old column-top rule with
      // full-coordination rates (c == c_max → factor 1).
      if (!s.occupied) {
        const int c = coordinatedNeighbors(i, j, k);
        if (c >= minCoord_ &&
            (!visibility_ || KmcVisibility::isVisible(lattice_, i, j, k))) {
          const double maxCoord = diamond_ ? 4.0 : 6.0;
          const double coordFactor = static_cast<double>(c) / maxCoord;
          KmcEvent ev;
          ev.type = KmcEventType::Deposit;
          ev.i0 = i;
          ev.j0 = j;
          ev.k0 = k;
          ev.rate =
              params_.attachRate() * coordFactor * (1.0 - 0.3 * xGe_);
          evList.push_back(ev);
          totalRate += ev.rate;
        }
      } else {
        int kTop = -1;
        for (int kk = nz - 1; kk >= 0; --kk) {
          if (lattice_.at(i, j, kk).occupied) {
            kTop = kk;
            break;
          }
        }
        if (k == kTop) {
          const int c = coordinatedNeighbors(i, j, k);
          const double maxCoord = diamond_ ? 4.0 : 6.0;
          const double coordFactor = static_cast<double>(c) / maxCoord;
          KmcEvent ev;
          ev.type = KmcEventType::Desorb;
          ev.i0 = i;  ev.j0 = j;  ev.k0 = kTop;
          // Weakly bonded atoms desorb faster (fewer bonds → higher rate).
          ev.rate = params_.desorbRate() * (1.0 - 0.5 * coordFactor);
          evList.push_back(ev);
          totalRate += ev.rate;
          if (s.species != KmcTwin) {
            KmcEvent ev2;
            ev2.type = KmcEventType::Twin;
            ev2.i0 = i;  ev2.j0 = j;  ev2.k0 = kTop;
            // Twins form preferentially on well-coordinated {111}-like sites.
            ev2.rate = params_.twinRate() * (c >= 3 ? 1.0 : 0.2);
            evList.push_back(ev2);
            totalRate += ev2.rate;
          }
        }
      }
    }

    double oldTotal = siteTotalRates_[siteIdx];
    siteTotalRates_[siteIdx] = totalRate;
    fenwick_.add(siteIdx, totalRate - oldTotal);
  }

  void rebuildAffectedSites(const KmcEvent &e) {
    std::vector<std::size_t> sites;
    const std::size_t idx0 = linearIndex(e.i0, e.j0, e.k0);
    const std::size_t idx1 = linearIndex(e.i1, e.j1, e.k1);
    sites.push_back(idx0);
    if (idx1 != idx0)
      sites.push_back(idx1);

    auto addNeighbors = [&](int ci, int cj, int ck) {
      forEachNeighbor(ci, cj, ck, [&](int ni, int nj, int nk) {
        sites.push_back(linearIndex(ni, nj, nk));
      });
      // In epitaxy mode, changing any site in column (ci, cj) alters kTop for that column,
      // so rebuild all sites in column (ci, cj).
      if (epitaxy_) {
        for (int kk = 0; kk < lattice_.nz(); ++kk) {
          sites.push_back(linearIndex(ci, cj, kk));
        }
      }
    };

    addNeighbors(e.i0, e.j0, e.k0);
    if (idx1 != idx0)
      addNeighbors(e.i1, e.j1, e.k1);

    std::sort(sites.begin(), sites.end());
    sites.erase(std::unique(sites.begin(), sites.end()), sites.end());
    for (auto idx : sites)
      buildSiteEvents(idx);
  }

  void apply(const KmcEvent &e) {
    if (e.type == KmcEventType::Hop) {
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      auto &b = lattice_.at(e.i1, e.j1, e.k1);
      b.species = a.species;
      b.occupied = true;
      a.occupied = false;
      a.species = 0;
    } else if (e.type == KmcEventType::Recombine) {
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      auto &b = lattice_.at(e.i1, e.j1, e.k1);
      a.occupied = false;
      a.species = 0;
      b.occupied = false;
      b.species = 0;
      ++recombCount_;
    } else if (e.type == KmcEventType::Cluster) {
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      auto &b = lattice_.at(e.i1, e.j1, e.k1);
      // Keep cluster on lower-index site; free the other.
      a.species = KmcCluster311;
      a.occupied = true;
      b.occupied = false;
      b.species = 0;
      ++clusterCount_;
    } else if (e.type == KmcEventType::Dissociate) {
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      auto &b = lattice_.at(e.i1, e.j1, e.k1);
      a.species = KmcInterstitial;
      a.occupied = true;
      b.species = KmcInterstitial;
      b.occupied = true;
      ++dissocCount_;
    } else if (e.type == KmcEventType::Deposit) {
      // Epitaxial surface attachment: fill the empty site with the growth
      // species (Si) or Ge with probability xGe (SiGe composition).
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      a.occupied = true;
      std::uniform_real_distribution<double> U(0.0, 1.0);
      a.species = (U(rng_) < xGe_) ? KmcGe : growthSpecies_;
      ++depositCount_;
    } else if (e.type == KmcEventType::Desorb) {
      // Surface desorption: remove the surface atom.
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      a.occupied = false;
      a.species = 0;
      ++desorbCount_;
    } else if (e.type == KmcEventType::Twin) {
      // Twin-defect formation: mark the surface site as a twin (stacking fault).
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      a.species = KmcTwin;
      ++twinCount_;
    }
  }

  KmcLattice lattice_;
  KmcParameters params_;
  std::vector<std::vector<KmcEvent>> eventsBySite_;
  std::vector<double> siteTotalRates_;
  viennaps::FenwickTree<double> fenwick_;
  bool fenwickInitialized_ = false;
  double time_ = 0.0;
  int steps_ = 0;
  int recombCount_ = 0;
  int clusterCount_ = 0;
  int dissocCount_ = 0;
  int depositCount_ = 0;
  int desorbCount_ = 0;
  int twinCount_ = 0;
  bool diamond_ = false;
  bool recomb_ = true;
  bool cluster_ = true;
  bool epitaxy_ = false;
  int growthSpecies_ = KmcSi;
  double xGe_ = 0.0;
  int minCoord_ = 1;
  bool visibility_ = false;
  std::mt19937 rng_;
};

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
          if (pois(rng) > 0) {
            auto &s = lat.at(i, j, k);
            s.occupied = true;
            s.species = speciesCode;
          }
        }
  }
};

} // namespace viennaps
