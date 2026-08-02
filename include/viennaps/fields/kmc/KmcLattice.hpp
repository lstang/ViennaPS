#pragma once

/// KmcLattice — Si diamond lattice sites for atomistic KMC.

#include <cstddef>
#include <cstdint>
#include <vector>

namespace viennaps {

struct KmcSite {
  double x = 0, y = 0, z = 0;
  int species = 0; // 0 empty/Si, >0 dopant/defect codes
  bool occupied = false;
};

class KmcLattice {
public:
  void resize(int nx, int ny, int nz, double a0 = 5.43e-8) {
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    a0_ = a0;
    sites_.assign(static_cast<std::size_t>(nx * ny * nz), KmcSite{});
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
          auto &s = at(i, j, k);
          s.x = i * a0_;
          s.y = j * a0_;
          s.z = k * a0_;
        }
  }

  int nx() const { return nx_; }
  int ny() const { return ny_; }
  int nz() const { return nz_; }
  std::size_t size() const { return sites_.size(); }

  KmcSite &at(int i, int j, int k) {
    return sites_[static_cast<std::size_t>((k * ny_ + j) * nx_ + i)];
  }
  const KmcSite &at(int i, int j, int k) const {
    return sites_[static_cast<std::size_t>((k * ny_ + j) * nx_ + i)];
  }

  int countSpecies(int code) const {
    int n = 0;
    for (const auto &s : sites_)
      if (s.occupied && s.species == code)
        ++n;
    return n;
  }

  /// Mark sites as diamond sublattice A (parity even) or B (parity odd).
  /// Neighbors between A–B approximate the two FCC sublattices of diamond.
  int sublattice(int i, int j, int k) const {
    return ((i + j + k) & 1); // 0 = A, 1 = B
  }

  /// 1D depth profile of species count per k-plane.
  std::vector<int> profile1D(int speciesCode) const {
    std::vector<int> p(static_cast<std::size_t>(nz_), 0);
    for (int k = 0; k < nz_; ++k)
      for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i)
          if (at(i, j, k).occupied && at(i, j, k).species == speciesCode)
            ++p[static_cast<std::size_t>(k)];
    return p;
  }

private:
  int nx_ = 0, ny_ = 0, nz_ = 0;
  double a0_ = 5.43e-8;
  std::vector<KmcSite> sites_;
};

class KmcVisibility {
public:
  /// Simple z-buffer shadowing: site visible if no occupied site above.
  static bool isVisible(const KmcLattice &lat, int i, int j, int k) {
    for (int kk = k + 1; kk < lat.nz(); ++kk)
      if (lat.at(i, j, kk).occupied)
        return false;
    return true;
  }
};

} // namespace viennaps
