#pragma once

/// Stress kernels for multiphysics process simulation.
///
/// ViscoelasticStressKernel:
///   Maxwell relaxation + growth stress bookkeeping on field totals / profiles.
///
/// ElasticStressKernel:
///   Isotropic elastic residual (Hooke-like hydrostatic proxy).
///   When MFEM is available, projects stress onto GridFunction and can run a
///   simple DiffusionIntegrator-style residual for stress smoothing.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

template <class NumericType>
class ViscoelasticStressKernel : public PhysicsKernel<NumericType> {
public:
  ViscoelasticStressKernel(NumericType temperatureK = 1273.15) : T_(temperatureK) {
    this->setName("ViscoelasticStress");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }

  void setup() override {
    if (this->material_) {
      E_ = this->material_->getYoungModulus("Si", T_);
      nu_ = this->material_->getPoissonRatio("Si", T_);
      eta_ = this->material_->getViscosity("SiO2", T_);
      alpha_ = this->material_->getCTE("Si", T_);
    } else {
      E_ = 130.0;
      nu_ = 0.28;
      eta_ = 1e13;
      alpha_ = 2.6e-6;
    }
    bulk_ = E_ / (3.0 * (1.0 - 2.0 * nu_));
    std::cout << "[ViscoelasticStressKernel] E=" << E_ << " GPa, nu=" << nu_
              << ", eta=" << eta_ << ", bulk=" << bulk_ << " at T=" << T_ << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;

    this->field_->addSpecies("HydrostaticStress");
    this->field_->addSpecies("GrowthStress");

    NumericType growthStress =
        this->material_ ? this->material_->getProperty("Si", "GrowthStress", T_) : 300.0e6;
    NumericType currentGrowth = this->field_->getTotalDose("GrowthStress");
    // Absolute bookkeeping: set totals (do not accumulate size-1 injects)
    NumericType newGrowth = currentGrowth + growthStress * 0.01 * dt;
    this->field_->setSpeciesDose("GrowthStress", newGrowth);

    NumericType G = E_ / (2.0 * (1.0 + nu_));
    NumericType tau = (eta_ > 0) ? (eta_ / G) * 1e-9 : 1e10;
    NumericType currentStress = this->field_->getTotalDose("HydrostaticStress");
    NumericType relaxed =
        currentStress * std::exp(-dt / std::max(tau, NumericType(1e-6)));
    NumericType totalStress = relaxed + newGrowth * 0.1;
    this->field_->setSpeciesDose("HydrostaticStress", totalStress);

    std::cout << "[ViscoelasticStressKernel] stress evolved dt=" << dt
              << "  hydrostatic=" << totalStress / 1e6 << " MPa (relaxed)\n";
  }

  void addToRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) &&
        myIndex < static_cast<int>(y.size())) {
      NumericType G = E_ / (2.0 * (1.0 + nu_));
      NumericType tau = (eta_ > 0) ? (eta_ / G) * 1e-9 : 1e10;
      ydot[myIndex] += -y[myIndex] / std::max(tau, NumericType(1e-6));
    }
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto off = this->field_->getSpeciesOffset("HydrostaticStress");
    if (off == static_cast<std::size_t>(-1)) return;
    auto n = this->field_->getProfileSize();
    if (off + n > y.size()) return;
    NumericType G = E_ / (2.0 * (1.0 + nu_));
    NumericType tau = (eta_ > 0) ? (eta_ / G) * 1e-9 : 1e10;
    NumericType invTau = NumericType(1) / std::max(tau, NumericType(1e-6));
    // Scale for demo CVODE stability
    invTau = std::min(invTau, NumericType(0.2));
    for (std::size_t i = 0; i < n; ++i) {
      ydot[off + i] += -invTau * y[off + i];
    }
  }

private:
  NumericType T_ = 1273.15;
  NumericType E_ = 130.0, nu_ = 0.28, eta_ = 1e13, alpha_ = 2.6e-6;
  NumericType bulk_ = 0;
};

/// Isotropic elastic residual + optional MFEM GridFunction projection.
template <class NumericType>
class ElasticStressKernel : public PhysicsKernel<NumericType> {
public:
  ElasticStressKernel(NumericType temperatureK = 1273.15) : T_(temperatureK) {
    this->setName("ElasticStress");
  }

  void setTemperature(NumericType T_K) { T_ = T_K; }
  void setMismatchStrain(NumericType eps) { mismatch_ = eps; }

  void setup() override {
    if (this->material_) {
      E_ = this->material_->getYoungModulus("Si", T_);
      nu_ = this->material_->getPoissonRatio("Si", T_);
    } else {
      E_ = 130.0;
      nu_ = 0.28;
    }
    // Plane-strain bulk-like modulus (GPa units consistent with mat system)
    lambda_ = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
    mu_ = E_ / (2.0 * (1.0 + nu_));
    std::cout << "[ElasticStressKernel] E=" << E_ << " nu=" << nu_
              << " lambda=" << lambda_ << " mu=" << mu_ << " mismatch=" << mismatch_
              << "\n";
  }

  void evolve(NumericType dt) override {
    if (!this->field_) return;
    this->field_->addSpecies("HydrostaticStress");
    this->field_->addSpecies("ElasticStress");

    // Hydrostatic proxy from isotropic Hooke: sigma_h ~ (3*lambda+2*mu)*eps
    NumericType K = (3.0 * lambda_ + 2.0 * mu_) / 3.0;
    NumericType sigma = K * mismatch_ * 1e9; // GPa * strain -> Pa-like

    // Mild time relaxation toward elastic equilibrium
    NumericType current = this->field_->getTotalDose("ElasticStress");
    NumericType target = sigma;
    NumericType updated = current + (target - current) * std::min(NumericType(1), dt * NumericType(0.1));

    this->field_->setSpeciesDose("ElasticStress", updated);
    // Couple into hydrostatic channel for oxidation adapter (absolute set)
    this->field_->setSpeciesDose("HydrostaticStress", updated * NumericType(0.5));

#ifdef VIENNAPS_HAS_MFEM
    // Project hydrostatic residual onto MFEM GF if mesh exists
    if (this->field_->getMesh() && this->field_->getFESpace()) {
      auto* gf = this->field_->getGridFunction("HydrostaticStress");
      if (gf) {
        mfem::ConstantCoefficient sigmaCoef(static_cast<double>(updated * NumericType(0.5)));
        gf->ProjectCoefficient(sigmaCoef);
        mfem::BilinearForm mass(this->field_->getFESpace());
        mass.AddDomainIntegrator(new mfem::MassIntegrator);
        mass.Assemble();
        mass.Finalize();
        mfem::Vector x = *gf, Mx(x.Size());
        mass.Mult(x, Mx);
        double energy = mfem::InnerProduct(x, Mx);
        std::cout << "[ElasticStressKernel] MFEM mass energy proxy=" << energy << "\n";
      }
    }
#endif

    std::cout << "[ElasticStressKernel] elastic=" << updated / 1e6 << " MPa dt=" << dt
              << "\n";
  }

  void addToRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                std::vector<NumericType>& ydot, int myIndex) override {
    if (myIndex >= 0 && myIndex < static_cast<int>(ydot.size()) &&
        myIndex < static_cast<int>(y.size())) {
      NumericType K = (3.0 * lambda_ + 2.0 * mu_) / 3.0;
      NumericType target = K * mismatch_;
      ydot[myIndex] += NumericType(0.1) * (target - y[myIndex]);
    }
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto off = this->field_->getSpeciesOffset("ElasticStress");
    if (off == static_cast<std::size_t>(-1))
      off = this->field_->getSpeciesOffset("HydrostaticStress");
    if (off == static_cast<std::size_t>(-1)) return;
    auto n = this->field_->getProfileSize();
    if (off + n > y.size()) return;

    NumericType K = (3.0 * lambda_ + 2.0 * mu_) / 3.0;
    NumericType target = K * mismatch_ * 1e9;
    for (std::size_t i = 0; i < n; ++i) {
      // Drive toward elastic target + light spatial smoothing
      ydot[off + i] += NumericType(0.05) * (target - y[off + i]);
      if (i > 0 && i + 1 < n) {
        ydot[off + i] +=
            NumericType(0.01) * (y[off + i - 1] - NumericType(2) * y[off + i] + y[off + i + 1]);
      }
    }
  }

private:
  NumericType T_ = 1273.15;
  NumericType E_ = 130.0, nu_ = 0.28;
  NumericType lambda_ = 0, mu_ = 0;
  NumericType mismatch_ = NumericType(0.002); // 0.2% mismatch default
};

} // namespace viennaps
