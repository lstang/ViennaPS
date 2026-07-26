#pragma once

/// KmcAtomisticEngine — BKL KMC with Hop, Recombine, Cluster, Dissociate.
///
/// Event selection uses Fenwick tree over per-site total rates for O(log N)
/// site selection, with O(1)-per-affected-site incremental rebuild after each
/// event. Species codes: 0 empty, 1 interstitial, 2 vacancy, 3 {311}-like.

#include "KmcEvent.hpp"
#include "KmcLattice.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace viennaps {

// Fenwick Tree (Binary Indexed Tree) for O(log N) prefix sums and updates.
template <class Real>
class FenwickTree {
public:
  explicit FenwickTree(std::size_t n = 0) { init(n); }

  void init(std::size_t n) {
    n_ = n;
    bit_.assign(n + 1, Real(0));
  }

  void add(std::size_t i, double delta) {
    for (std::size_t idx = i + 1; idx <= n_; idx += idx & -idx)
      bit_[idx] += delta;
  }

  double sum(std::size_t i) const {
    double res = 0;
    for (std::size_t idx = i + 1; idx > 0; idx -= idx & -idx)
      res += bit_[idx];
    return res;
  }

  double total() const { return sum(n_ - 1); }

  std::size_t lower_bound(double target) const {
    if (target <= 0)
      return 0;
    std::size_t idx = 0;
    double s = 0;
    std::size_t bit = 1;
    while (bit << 1 <= n_)
      bit <<= 1;
    for (; bit > 0; bit >>= 1) {
      std::size_t next = idx + bit;
      if (next <= n_ && s + bit_[next] < target) {
        idx = next;
        s += bit_[next];
      }
    }
    return std::min(idx, n_ - 1);
  }

  std::size_t size() const { return n_; }

private:
  std::size_t n_ = 0;
  std::vector<double> bit_;
};

class KmcAtomisticEngine {
public:
  explicit KmcAtomisticEngine(unsigned seed = 42) : rng_(seed) {}

  void setLattice(KmcLattice lattice) { lattice_ = std::move(lattice); }
  KmcLattice &lattice() { return lattice_; }
  const KmcLattice &lattice() const { return lattice_; }

  void setParameters(KmcParameters p) { params_ = std::move(p); }
  const KmcParameters &parameters() const { return params_; }

  void setDiamondNeighbors(bool on) { diamond_ = on; }
  bool diamondNeighbors() const { return diamond_; }

  void setRecombinationEnabled(bool on) { recomb_ = on; }
  bool recombinationEnabled() const { return recomb_; }

  void setClusteringEnabled(bool on) { cluster_ = on; }
  bool clusteringEnabled() const { return cluster_; }

  /// Enable epitaxial surface events (Deposit/Desorb/Twin) as BKL events.
  void setEpitaxyEnabled(bool on) { epitaxy_ = on; }
  bool epitaxyEnabled() const { return epitaxy_; }
  /// Set the growth species code for Deposit events (default KmcSi=4).
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
    if (!chosen && !siteEvents.empty())
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
          const int a = static_cast<int>(siteIdx);
          const int b = (k1 * ny + j1) * nx + i1;
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
      int kTop = -1;
      for (int kk = nz - 1; kk >= 0; --kk) {
        if (lattice_.at(i, j, kk).occupied) {
          kTop = kk;
          break;
        }
      }
      if (kTop >= 0 && k == kTop + 1 && kTop + 1 < nz &&
          !lattice_.at(i, j, kTop + 1).occupied) {
        KmcEvent ev;
        ev.type = KmcEventType::Deposit;
        ev.i0 = i;     ev.j0 = j;     ev.k0 = kTop + 1;
        ev.rate = params_.attachRate();
        evList.push_back(ev);
        totalRate += ev.rate;
      }
      if (k == kTop && s.occupied) {
        KmcEvent ev;
        ev.type = KmcEventType::Desorb;
        ev.i0 = i;  ev.j0 = j;  ev.k0 = kTop;
        ev.rate = params_.desorbRate();
        evList.push_back(ev);
        totalRate += ev.rate;
        KmcEvent ev2;
        ev2.type = KmcEventType::Twin;
        ev2.i0 = i;  ev2.j0 = j;  ev2.k0 = kTop;
        ev2.rate = params_.twinRate();
        evList.push_back(ev2);
        totalRate += ev2.rate;
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
      // Epitaxial surface attachment: fill the empty site with growth species.
      auto &a = lattice_.at(e.i0, e.j0, e.k0);
      a.occupied = true;
      a.species = growthSpecies_;
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
  FenwickTree<double> fenwick_;
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

  /// Inverse-distance weighted smooth of discrete occupations onto a 1D/flat
  /// continuum grid of length `nOut` (depth bins along k).
  static std::vector<double> deatomizeIDW(const KmcLattice &lat,
                                          int speciesCode, int nOut,
                                          double power = 2.0) {
    std::vector<double> out(static_cast<std::size_t>(std::max(nOut, 1)), 0.0);
    if (nOut <= 0 || lat.nz() <= 0)
      return out;
    std::vector<double> wsum(out.size(), 0.0);
    for (int k = 0; k < lat.nz(); ++k) {
      int count = 0;
      for (int j = 0; j < lat.ny(); ++j)
        for (int i = 0; i < lat.nx(); ++i)
          if (lat.at(i, j, k).occupied &&
              lat.at(i, j, k).species == speciesCode)
            ++count;
      const double zk = static_cast<double>(k) / std::max(lat.nz() - 1, 1);
      for (int b = 0; b < nOut; ++b) {
        const double zb =
            static_cast<double>(b) / std::max(nOut - 1, 1);
        const double d = std::abs(zk - zb) + 1e-9;
        const double w = 1.0 / std::pow(d, power);
        out[static_cast<std::size_t>(b)] += w * count;
        wsum[static_cast<std::size_t>(b)] += w;
      }
    }
    for (std::size_t i = 0; i < out.size(); ++i)
      if (wsum[i] > 0)
        out[i] /= wsum[i];
    return out;
  }
};

/// Amorphous pocket: mark a cubic region as amorphous species (code 7).
class KmcAmorphousPocket {
public:
  static int implant(KmcLattice &lat, int i0, int j0, int k0, int r,
                     int amorphCode = 7) {
    int n = 0;
    for (int k = k0 - r; k <= k0 + r; ++k)
      for (int j = j0 - r; j <= j0 + r; ++j)
        for (int i = i0 - r; i <= i0 + r; ++i) {
          if (i < 0 || j < 0 || k < 0 || i >= lat.nx() || j >= lat.ny() ||
              k >= lat.nz())
            continue;
          if ((i - i0) * (i - i0) + (j - j0) * (j - j0) + (k - k0) * (k - k0) >
              r * r)
            continue;
          auto &s = lat.at(i, j, k);
          s.occupied = true;
          s.species = amorphCode;
          ++n;
        }
    return n;
  }
};

struct KmcReport {
  double time = 0;
  int steps = 0;
  int hopCount = 0;
  int occupied = 0;
  int recombCount = 0;
  int clusterCount = 0;
  int dissocCount = 0;
  double supersaturationI = 0.0;
  std::vector<int> depthProfileI;
  std::map<int, int> clusterSizeHistogram; // size→count (size=1 free I, 2=cluster)

  static KmcReport fromEngine(const KmcAtomisticEngine &eng,
                              double C_I_eq = 1e12) {
    KmcReport r;
    r.time = eng.time();
    r.steps = eng.steps();
    r.recombCount = eng.recombCount();
    r.clusterCount = eng.clusterCount();
    r.dissocCount = eng.dissocCount();
    r.occupied = 0;
    int nI = 0, nCl = 0;
    for (int k = 0; k < eng.lattice().nz(); ++k)
      for (int j = 0; j < eng.lattice().ny(); ++j)
        for (int i = 0; i < eng.lattice().nx(); ++i) {
          const auto &s = eng.lattice().at(i, j, k);
          if (!s.occupied)
            continue;
          ++r.occupied;
          if (s.species == 1)
            ++nI;
          if (s.species == 3)
            ++nCl;
        }
    r.depthProfileI = eng.lattice().profile1D(1);
    r.clusterSizeHistogram[1] = nI;
    r.clusterSizeHistogram[2] = nCl;
    // Supersaturation proxy: free-I count / (C_I_eq * volume units).
    const double volSites =
        static_cast<double>(eng.lattice().size());
    r.supersaturationI =
        (C_I_eq > 0 && volSites > 0)
            ? static_cast<double>(nI) / (C_I_eq * volSites * 1e-24 + 1e-30)
            : 0.0;
    return r;
  }
};

class KmcContinuumCoupler {
public:
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