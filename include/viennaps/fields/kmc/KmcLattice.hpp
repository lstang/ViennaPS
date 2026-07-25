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

private:
  int nx_ = 0, ny_ = 0, nz_ = 0;
  double a0_ = 5.43e-8;
  std::vector<KmcSite> sites_;
};

} // namespace viennaps
