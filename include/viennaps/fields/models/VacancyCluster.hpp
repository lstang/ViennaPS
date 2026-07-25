#pragma once

/// VacancyCluster — vacancy cluster (VC) growth / dissociation.
/// dC_VC/dt = k_f * C_V^m - k_r * C_VC

#include "../DiffusionModel.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class VacancyCluster : public DiffusionModel<NumericType> {
public:
  VacancyCluster(std::string cluster = "VC",
                 std::string vacancy = "Vacancy")
      : cluster_(std::move(cluster)), V_(std::move(vacancy)) {
    this->setName("VacancyCluster");
  }

  void setRates(NumericType kf, NumericType kr, int m = 2) {
    kf_ = kf;
    kr_ = kr;
    m_ = m;
  }

  void applyReactionStep(std::vector<NumericType> &V,
                         std::vector<NumericType> &VC,
                         NumericType dt) const {
    const std::size_t n = std::min(V.size(), VC.size());
    for (std::size_t i = 0; i < n; ++i) {
      const NumericType growth =
          kf_ * std::pow(std::max(V[i], NumericType(0)), m_);
      const NumericType dissoc = kr_ * VC[i];
      const NumericType dC = (growth - dissoc) * dt;
      VC[i] = std::max(NumericType(0), VC[i] + dC);
      V[i] = std::max(NumericType(0),
                      V[i] - static_cast<NumericType>(m_) * dC);
    }
  }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {cluster_, V_};
  }

private:
  std::string cluster_, V_;
  NumericType kf_ = NumericType(1e-20);
  NumericType kr_ = NumericType(1e-3);
  int m_ = 2;
};

} // namespace viennaps
