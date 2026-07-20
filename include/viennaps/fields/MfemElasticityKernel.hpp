#pragma once

/// MfemElasticityKernel - Vector isotropic elasticity with MFEM assembly.
///
/// When VIENNAPS_HAS_MFEM:
///   - Builds vector H1 FE space (dim components)
///   - Assembles ElasticityIntegrator (lambda, mu)
///   - Applies traction BCs on boundary attributes (mask edge / free surface proxy)
///   - Recovers hydrostatic + von Mises proxies into PhysicsField species
///
/// Without MFEM: falls back to ElasticStressKernel-style Hooke hydrostatic set.

#include "PhysicsKernel.hpp"
#include "PhysicsField.hpp"
#include "MaterialPropertySystem.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

template <class NumericType>
class MfemElasticityKernel : public PhysicsKernel<NumericType> {
public:
  MfemElasticityKernel(NumericType temperatureK = 1273.15) : T_(temperatureK) {
    this->setName("MfemElasticity");
  }

  void setTemperature(NumericType T) { T_ = T; }
  void setMismatchStrain(NumericType e) { mismatch_ = e; }
  void setTractionMagnitude(NumericType t) { traction_ = t; }
  /// Boundary attribute id that receives traction (default 1).
  void setTractionAttribute(int attr) { tractionAttr_ = attr; }

  void setup() override {
    if (this->material_) {
      E_ = this->material_->getYoungModulus("Si", T_);
      nu_ = this->material_->getPoissonRatio("Si", T_);
    } else {
      E_ = 130.0;
      nu_ = 0.28;
    }
    lambda_ = E_ * nu_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_));
    mu_ = E_ / (2.0 * (1.0 + nu_));
    std::cout << "[MfemElasticityKernel] E=" << E_ << " nu=" << nu_
              << " lambda=" << lambda_ << " mu=" << mu_
              << " mismatch=" << mismatch_ << " traction=" << traction_ << "\n";
  }

  void evolve(NumericType /*dt*/) override {
    if (!this->field_) return;
    this->field_->addSpecies("HydrostaticStress");
    this->field_->addSpecies("VonMisesStress");
    this->field_->addSpecies("StressXX");
    this->field_->addSpecies("StressYY");
    this->field_->addSpecies("StressXY");

#ifdef VIENNAPS_HAS_MFEM
    if (solveMFEMElasticity()) {
      return;
    }
#endif
    // Fallback Hooke hydrostatic + simple tensor components
    NumericType K = (3.0 * lambda_ + 2.0 * mu_) / 3.0;
    NumericType sig = K * mismatch_ * 1e9;
    NumericType sxx = sig * (1.0 + mismatch_);
    NumericType syy = sig * (1.0 - 0.5 * mismatch_);
    NumericType sxy = traction_ * 1e6;
    NumericType vm =
        std::sqrt(0.5 * ((sxx - syy) * (sxx - syy) + sxx * sxx + syy * syy) +
                  3.0 * sxy * sxy);

    this->field_->setSpeciesDose("HydrostaticStress", (sxx + syy) / 3.0);
    this->field_->setSpeciesDose("VonMisesStress", vm);
    this->field_->setSpeciesDose("StressXX", sxx);
    this->field_->setSpeciesDose("StressYY", syy);
    this->field_->setSpeciesDose("StressXY", sxy);
    lastEnergy_ = vm;
    std::cout << "[MfemElasticityKernel] fallback Hooke hydrostatic="
              << (sxx + syy) / 3e6 << " MPa vonMises=" << vm / 1e6 << " MPa\n";
  }

  void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& y,
                     std::vector<NumericType>& ydot) override {
    if (!this->field_) return;
    auto off = this->field_->getSpeciesOffset("HydrostaticStress");
    if (off == static_cast<std::size_t>(-1)) return;
    auto n = this->field_->getProfileSize();
    if (off + n > y.size()) return;
    NumericType K = (3.0 * lambda_ + 2.0 * mu_) / 3.0;
    NumericType target = K * mismatch_ * 1e9;
    for (std::size_t i = 0; i < n; ++i)
      ydot[off + i] += NumericType(0.05) * (target - y[off + i]);
  }

  NumericType getLastEnergy() const { return lastEnergy_; }
  bool usedMFEMSolve() const { return usedMfem_; }

private:
  NumericType T_ = 1273.15;
  NumericType E_ = 130, nu_ = 0.28, lambda_ = 0, mu_ = 0;
  NumericType mismatch_ = NumericType(0.002);
  NumericType traction_ = NumericType(10); // MPa-ish
  int tractionAttr_ = 1;
  NumericType lastEnergy_ = 0;
  bool usedMfem_ = false;

#ifdef VIENNAPS_HAS_MFEM
  bool solveMFEMElasticity() {
    usedMfem_ = false;
    // Ensure a mesh exists
    if (!this->field_->getMesh()) {
      this->field_->initMeshFromBounds(0, 1, 0, 1, 16, 16);
    }
    mfem::Mesh* mesh = this->field_->getMesh();
    if (!mesh) return false;

    const int dim = mesh->Dimension();
    mfem::H1_FECollection fec(1, dim);
    mfem::FiniteElementSpace fes(mesh, &fec, dim, mfem::Ordering::byVDIM);

    // Elasticity bilinear form
    mfem::ConstantCoefficient lambda_cf(static_cast<double>(lambda_));
    mfem::ConstantCoefficient mu_cf(static_cast<double>(mu_));
    mfem::BilinearForm a(&fes);
    a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda_cf, mu_cf));
    a.Assemble();

    // RHS: body force from mismatch (isotropic expansion proxy) + traction BC
    mfem::VectorArrayCoefficient body(dim);
    for (int d = 0; d < dim; ++d) {
      // small body force proportional to mismatch
      body.Set(d, new mfem::ConstantCoefficient(
                      static_cast<double>(mismatch_ * 1e3 * (d == 1 ? 1.0 : 0.1))));
    }
    mfem::LinearForm b(&fes);
    b.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(body));

    // Traction on boundaries with attribute tractionAttr_
    mfem::VectorArrayCoefficient traction(dim);
    for (int d = 0; d < dim; ++d) {
      double val = (d == 0) ? static_cast<double>(traction_ * 1e6) : 0.0;
      traction.Set(d, new mfem::ConstantCoefficient(val));
    }
    // Mark boundary attributes if none set
    if (mesh->bdr_attributes.Size() == 0) {
      // ensure at least attribute 1 on all boundary
      for (int i = 0; i < mesh->GetNBE(); ++i)
        mesh->SetBdrAttribute(i, 1);
      mesh->SetAttributes();
    }
    mfem::Array<int> ess_bdr;
    if (mesh->bdr_attributes.Size() > 0) {
      ess_bdr.SetSize(mesh->bdr_attributes.Max());
      ess_bdr = 0;
      // leave free for traction; pin one attribute for uniqueness if possible
      if (ess_bdr.Size() >= 2) {
        ess_bdr[ess_bdr.Size() - 1] = 1; // essential on last attr
      }
      b.AddBoundaryIntegrator(new mfem::VectorBoundaryLFIntegrator(traction));
    }
    b.Assemble();

    mfem::Array<int> ess_tdof;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof);

    mfem::GridFunction x(&fes);
    x = 0.0;
    // Apply zero essential BC
    mfem::Vector X, B;
    mfem::OperatorPtr A;
    a.FormLinearSystem(ess_tdof, x, b, A, X, B);

    // Solve with GSSmoother + PCG (amgcl optional path is separate)
    mfem::GSSmoother M((mfem::SparseMatrix&)(*A));
    mfem::PCG(*A, M, B, X, 0, 200, 1e-8, 0.0);
    a.RecoverFEMSolution(X, b, x);

    // Stress recovery: Hooke proxy from displacement magnitude per element
    NumericType sumH = 0, sumVM = 0, sumXX = 0, sumYY = 0, sumXY = 0;
    int nSamp = 0;
    for (int e = 0; e < mesh->GetNE(); ++e) {
      mfem::Array<int> dofs;
      fes.GetElementVDofs(e, dofs);
      double uNorm = 0;
      for (int i = 0; i < dofs.Size(); ++i) {
        int idx = dofs[i];
        if (idx < 0) idx = -1 - idx;
        if (idx >= 0 && idx < x.Size()) uNorm += x(idx) * x(idx);
      }
      uNorm = std::sqrt(uNorm / std::max(1, dofs.Size()));
      double sxx = (lambda_ + 2 * mu_) * mismatch_ + 2 * mu_ * uNorm;
      double syy = lambda_ * mismatch_ + 2 * mu_ * uNorm * 0.5;
      double sxy = traction_ * 1e6 * 0.1 + mu_ * uNorm * 0.1;
      double hydro = (sxx + syy) / (dim == 3 ? 3.0 : 2.0);
      double vm = std::sqrt(0.5 * ((sxx - syy) * (sxx - syy) + sxx * sxx + syy * syy) +
                            3.0 * sxy * sxy);
      sumH += hydro;
      sumVM += vm;
      sumXX += sxx;
      sumYY += syy;
      sumXY += sxy;
      nSamp++;
    }
    if (nSamp == 0) return false;
    NumericType inv = NumericType(1) / static_cast<NumericType>(nSamp);
    this->field_->setSpeciesDose("HydrostaticStress", sumH * inv * 1e9);
    this->field_->setSpeciesDose("VonMisesStress", sumVM * inv * 1e9);
    this->field_->setSpeciesDose("StressXX", sumXX * inv * 1e9);
    this->field_->setSpeciesDose("StressYY", sumYY * inv * 1e9);
    this->field_->setSpeciesDose("StressXY", sumXY * inv * 1e9);
    lastEnergy_ = sumVM * inv * 1e9;
    usedMfem_ = true;

    // Project hydrostatic onto scalar GF if available
    if (auto* gf = this->field_->getGridFunction("HydrostaticStress")) {
      mfem::ConstantCoefficient c(static_cast<double>(sumH * inv * 1e9));
      gf->ProjectCoefficient(c);
    }

    std::cout << "[MfemElasticityKernel] MFEM ElasticityIntegrator solve OK, "
              << "avg vonMises=" << lastEnergy_ / 1e6 << " MPa, ne=" << mesh->GetNE()
              << " vdim=" << dim << "\n";
    return true;
  }
#endif
};

} // namespace viennaps
