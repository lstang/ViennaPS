#pragma once

/// DiffusionKernel - Transport kernel for dopant/defect diffusion.
/// Field-packed RHS: 1D discrete Laplacian on the species profile segment.
/// Explicit evolve: mild profile smoothing + total-dose conserving scale.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

namespace viennaps {

template <class NumericType>
class DiffusionKernel : public PhysicsKernel<NumericType> {
public:
  DiffusionKernel(const std::string& species, NumericType temperatureK = 1273.15)
      : species_(species), T_(temperatureK) {
    this->setName("Diffusion(" + species + ")");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }

  void setup() override {
    if (this->material_) {
      D_ = this->material_->getDiffusivity(species_, "Si", T_);
    } else {
      D_ = NumericType(1e-14);
    }
    std::cout << "[DiffusionKernel] " << species_ << " D=" << D_ << " cm^2/s at T=" << T_
              << " K\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;
    auto* prof = this->field_->getProfileMutable(species_);
    if (!prof || prof->empty()) return;

    auto totalBefore = this->field_->getTotalDose(species_);
    if (totalBefore <= 0) return;

    // One explicit smoothing pass (conserves total approximately)
    std::vector<NumericType> next = *prof;
    NumericType alpha = std::min(NumericType(0.25), D_ * dt * NumericType(1e12));
    for (std::size_t i = 1; i + 1 < prof->size(); ++i) {
      next[i] = (*prof)[i] + alpha * ((*prof)[i - 1] - NumericType(2) * (*prof)[i] + (*prof)[i + 1]);
      next[i] = std::max(NumericType(0), next[i]);
    }
    *prof = next;

    // Keep dose ~ stable under pure diffusion
    NumericType sum = 0;
    for (auto v : *prof) sum += v;
    if (sum > 0 && totalBefore > 0) {
      NumericType scale = totalBefore / sum;
      for (auto& v : *prof) v *= scale;
    }
    this->field_->refreshDose(species_);

    std::cout << "[DiffusionKernel] " << species_ << " evolved dt=" << dt
              << "  total " << totalBefore << " -> " << this->field_->getTotalDose(species_)
              << "\n";
  }

  void addToRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) &&
        myIndex < static_cast<int>(y.size())) {
      NumericType alpha = (D_ > 0) ? (D_ * NumericType(1e10)) : NumericType(0.01);
      if (alpha > NumericType(1)) alpha = NumericType(1);
      ydot[myIndex] += -alpha * y[myIndex];
    }
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto off = this->field_->getSpeciesOffset(species_);
    if (off == static_cast<std::size_t>(-1)) return;
    auto n = this->field_->getProfileSize();
    if (off + n > y.size() || off + n > ydot.size()) return;

    // Dimensionless diffusion rate for packed residual (stable for CVODE demo)
    NumericType Dscale = std::min(NumericType(0.5),
                                  std::max(NumericType(1e-4), D_ * NumericType(1e12)));
    for (std::size_t i = 0; i < n; ++i) {
      NumericType left = (i > 0) ? y[off + i - 1] : y[off + i];
      NumericType right = (i + 1 < n) ? y[off + i + 1] : y[off + i];
      NumericType center = y[off + i];
      ydot[off + i] += Dscale * (left - NumericType(2) * center + right);
    }
  }

private:
  std::string species_;
  NumericType T_ = 1273.15;
  NumericType D_ = 0;
};

} // namespace viennaps
