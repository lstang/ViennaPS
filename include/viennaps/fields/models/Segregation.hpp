#pragma once

/// Segregation — two-sided multi-material interface exchange (MOOSE
/// InterfaceReaction). Dual-species representation: H1 cannot hold C2/C1≠1
/// on one continuous field, so mat1/mat2 concentrations are separate species
/// coupled on interior faces between the two attributes.

#include "../DiffusionModel.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

template <class NumericType>
class SegregationCondition {
public:
  void setSegregationCoefficient(NumericType m, NumericType k0) {
    m_ = m;
    k0_ = k0;
    kf_ = m * k0;
    kb_ = k0;
  }

  NumericType segregationCoefficient() const { return m_; }
  NumericType kf() const { return kf_; }
  NumericType kb() const { return kb_; }
  NumericType equilibriumRatio() const {
    return (kb_ != NumericType(0)) ? (kf_ / kb_) : m_;
  }
  NumericType rate(NumericType C1, NumericType C2) const {
    return kf_ * C1 - kb_ * C2;
  }

#ifdef VIENNAPS_HAS_MFEM
  void assembleElementSide(
      mfem::Vector &R_elem, mfem::DenseMatrix &K_ee, mfem::DenseMatrix &K_en,
      const mfem::Vector &C1_dofs, const mfem::Vector &C2_dofs,
      const mfem::Vector &shape_e, const mfem::Vector &shape_n,
      double weight) const {
    const double C1 = mfem::InnerProduct(C1_dofs, shape_e);
    const double C2 = mfem::InnerProduct(C2_dofs, shape_n);
    const double r =
        static_cast<double>(kf_) * C1 - static_cast<double>(kb_) * C2;
    for (int i = 0; i < shape_e.Size(); ++i) {
      R_elem(i) += weight * shape_e(i) * r;
      for (int j = 0; j < shape_e.Size(); ++j)
        K_ee(i, j) +=
            weight * shape_e(i) * static_cast<double>(kf_) * shape_e(j);
      for (int j = 0; j < shape_n.Size(); ++j)
        K_en(i, j) +=
            weight * shape_e(i) * (-static_cast<double>(kb_)) * shape_n(j);
    }
  }

  void assembleNeighborSide(
      mfem::Vector &R_nbr, mfem::DenseMatrix &K_nn, mfem::DenseMatrix &K_ne,
      const mfem::Vector &C1_dofs, const mfem::Vector &C2_dofs,
      const mfem::Vector &shape_e, const mfem::Vector &shape_n,
      double weight) const {
    const double C1 = mfem::InnerProduct(C1_dofs, shape_e);
    const double C2 = mfem::InnerProduct(C2_dofs, shape_n);
    const double r =
        static_cast<double>(kf_) * C1 - static_cast<double>(kb_) * C2;
    for (int i = 0; i < shape_n.Size(); ++i) {
      R_nbr(i) += -weight * shape_n(i) * r;
      for (int j = 0; j < shape_n.Size(); ++j)
        K_nn(i, j) +=
            weight * shape_n(i) * static_cast<double>(kb_) * shape_n(j);
      for (int j = 0; j < shape_e.Size(); ++j)
        K_ne(i, j) +=
            -weight * shape_n(i) * static_cast<double>(kf_) * shape_e(j);
    }
  }

  /// Mass-conserving interface residual for engine form M du/dt = -K u + R.
  /// When r = kf*C1 - kb*C2 > 0, mass leaves mat1 and enters mat2:
  ///   R_mat1 gets -r * w * phi,  R_mat2 gets +r * w * phi  on the face.
  /// side=1 assembles into R for mat1 species; side=2 for mat2.
  void assembleInterfaceResidual(
      mfem::ParLinearForm &R, const mfem::ParGridFunction &C1,
      const mfem::ParGridFunction &C2, mfem::ParMesh &mesh, int attr1, int attr2,
      int side) const {
    // MFEM FESpace() is const; face/vdof APIs need non-const Mesh/FES.
    mfem::ParFiniteElementSpace *fes =
        dynamic_cast<mfem::ParFiniteElementSpace *>(
            const_cast<mfem::FiniteElementSpace *>(C1.FESpace()));
    if (!fes || C2.FESpace() != C1.FESpace())
      return;

    for (int f = 0; f < mesh.GetNumFaces(); ++f) {
      int eA = -1, eB = -1;
      mesh.GetFaceElements(f, &eA, &eB);
      if (eA < 0 || eB < 0)
        continue;

      const int aA = mesh.GetAttribute(eA);
      const int aB = mesh.GetAttribute(eB);
      int e1 = -1, e2 = -1;
      if (aA == attr1 && aB == attr2) {
        e1 = eA;
        e2 = eB;
      } else if (aA == attr2 && aB == attr1) {
        e1 = eB;
        e2 = eA;
      } else {
        continue;
      }

      mfem::Array<int> vdofs1, vdofs2;
      fes->GetElementVDofs(e1, vdofs1);
      fes->GetElementVDofs(e2, vdofs2);
      if (vdofs1.Size() == 0 || vdofs2.Size() == 0)
        continue;

      // Element-average concentrations (robust H1 dual-field coupling).
      double c1sum = 0.0, c2sum = 0.0;
      int n1 = 0, n2 = 0;
      for (int i = 0; i < vdofs1.Size(); ++i) {
        const int id = vdofs1[i];
        if (id < 0)
          continue;
        c1sum += C1(id);
        ++n1;
      }
      for (int i = 0; i < vdofs2.Size(); ++i) {
        const int id = vdofs2[i];
        if (id < 0)
          continue;
        c2sum += C2(id);
        ++n2;
      }
      if (n1 == 0 || n2 == 0)
        continue;
      const double C1v = c1sum / n1;
      const double C2v = c2sum / n2;
      const double r =
          static_cast<double>(kf_) * C1v - static_cast<double>(kb_) * C2v;

      // Face measure (length in 2D, area in 3D).
      mfem::FaceElementTransformations *FT =
          mesh.GetFaceElementTransformations(f);
      if (!FT)
        continue;
      const mfem::IntegrationRule &ir =
          mfem::IntRules.Get(FT->GetGeometryType(), 2);
      double faceMeas = 0.0;
      for (int i = 0; i < ir.GetNPoints(); ++i) {
        const mfem::IntegrationPoint &ip = ir.IntPoint(i);
        FT->SetAllIntPoints(&ip);
        faceMeas += ip.weight * FT->Weight();
      }

      // Distribute ±r*faceMeas uniformly over element dofs (mass conserving).
      const double load1 = -r * faceMeas; // leave mat1
      const double load2 = +r * faceMeas; // enter mat2
      if (side == 1) {
        const double per = load1 / static_cast<double>(n1);
        for (int i = 0; i < vdofs1.Size(); ++i) {
          const int id = vdofs1[i];
          if (id >= 0)
            R[id] += per;
        }
      } else if (side == 2) {
        const double per = load2 / static_cast<double>(n2);
        for (int i = 0; i < vdofs2.Size(); ++i) {
          const int id = vdofs2[i];
          if (id >= 0)
            R[id] += per;
        }
      }
    }
  }
#endif

private:
  NumericType m_ = NumericType(0.1);
  NumericType k0_ = NumericType(1e-3);
  NumericType kf_ = NumericType(1e-4);
  NumericType kb_ = NumericType(1e-3);
};

/// Dual-species segregation model for DiffusionEngine.
template <class NumericType>
class Segregation : public DiffusionModel<NumericType> {
public:
  Segregation(std::string speciesMat1 = "Boron_Si",
              std::string speciesMat2 = "Boron_Ox", int attrMat1 = 1,
              int attrMat2 = 2)
      : sp1_(std::move(speciesMat1)), sp2_(std::move(speciesMat2)),
        attr1_(attrMat1), attr2_(attrMat2) {
    this->setName("Segregation(" + sp1_ + "," + sp2_ + ")");
  }

  void setSegregationCoefficient(NumericType m, NumericType k0) {
    cond_.setSegregationCoefficient(m, k0);
  }

  SegregationCondition<NumericType> &condition() { return cond_; }
  const SegregationCondition<NumericType> &condition() const { return cond_; }

  int numSpecies() const override { return 2; }
  std::vector<std::string> speciesNames() const override {
    return {sp1_, sp2_};
  }
  std::vector<int> applicableAttributes() const override {
    return {attr1_, attr2_};
  }

  const std::string &speciesMat1() const { return sp1_; }
  const std::string &speciesMat2() const { return sp2_; }

#ifdef VIENNAPS_HAS_MFEM
  void assembleMass(mfem::ParBilinearForm &M) const override {
    massCoef_ = std::make_unique<mfem::ConstantCoefficient>(1.0);
    M.AddDomainIntegrator(new mfem::MassIntegrator(*massCoef_));
  }

  /// Stable per-step operator-split exchange inside the engine solve loop.
  /// Well-mixed two-compartment model with analytic exponential relaxation
  /// toward C2/C1 = m, conserving total dose M = ∫C1 + ∫C2.
  void applyOperatorSplitStep(
      mfem::ParGridFunction &C1, mfem::ParGridFunction &C2, mfem::ParMesh &mesh,
      double dose1, double dose2, double dt) const {
    if (dt <= 0.0)
      return;
    const double M = dose1 + dose2;
    if (M <= 0.0)
      return;

    // Volumes of each material (unit density → measure of elements).
    double V1 = 0.0, V2 = 0.0, Aiface = 0.0;
    for (int e = 0; e < mesh.GetNE(); ++e) {
      mfem::ElementTransformation *T = mesh.GetElementTransformation(e);
      const mfem::IntegrationRule &ir =
          mfem::IntRules.Get(T->GetGeometryType(), 2);
      double vol = 0.0;
      for (int i = 0; i < ir.GetNPoints(); ++i) {
        const mfem::IntegrationPoint &ip = ir.IntPoint(i);
        T->SetIntPoint(&ip);
        vol += ip.weight * T->Weight();
      }
      if (mesh.GetAttribute(e) == attr1_)
        V1 += vol;
      else if (mesh.GetAttribute(e) == attr2_)
        V2 += vol;
    }
    if (V1 <= 0.0 || V2 <= 0.0)
      return;

    for (int f = 0; f < mesh.GetNumFaces(); ++f) {
      int eA = -1, eB = -1;
      mesh.GetFaceElements(f, &eA, &eB);
      if (eA < 0 || eB < 0)
        continue;
      const int aA = mesh.GetAttribute(eA);
      const int aB = mesh.GetAttribute(eB);
      if (!((aA == attr1_ && aB == attr2_) || (aA == attr2_ && aB == attr1_)))
        continue;
      mfem::FaceElementTransformations *FT =
          mesh.GetFaceElementTransformations(f);
      if (!FT)
        continue;
      const mfem::IntegrationRule &ir =
          mfem::IntRules.Get(FT->GetGeometryType(), 2);
      for (int i = 0; i < ir.GetNPoints(); ++i) {
        const mfem::IntegrationPoint &ip = ir.IntPoint(i);
        FT->SetAllIntPoints(&ip);
        Aiface += ip.weight * FT->Weight();
      }
    }
    if (Aiface <= 0.0)
      Aiface = 1.0; // fallback

    const double m = static_cast<double>(cond_.segregationCoefficient());
    const double kf = static_cast<double>(cond_.kf());
    const double kb = static_cast<double>(cond_.kb());
    // C1_eq = M / (V1 + m V2), C2_eq = m C1_eq
    const double C1_eq = M / (V1 + m * V2);
    const double C2_eq = m * C1_eq;
    const double C1_0 = dose1 / V1;
    const double C2_0 = dose2 / V2;
    // Linearization rate about equilibrium for two-compartment exchange.
    const double lambda = Aiface * (kf / V1 + kb / V2);
    const double alpha = 1.0 - std::exp(-lambda * dt);
    const double C1n = C1_0 + alpha * (C1_eq - C1_0);
    const double C2n = C2_0 + alpha * (C2_eq - C2_0);

    // Restrict each species to its material: zero elsewhere, uniform on-side.
    C1 = 0.0;
    C2 = 0.0;
    mfem::ParFiniteElementSpace *fes =
        dynamic_cast<mfem::ParFiniteElementSpace *>(
            const_cast<mfem::FiniteElementSpace *>(C1.FESpace()));
    mfem::Array<int> vdofs;
    for (int e = 0; e < mesh.GetNE(); ++e) {
      fes->GetElementVDofs(e, vdofs);
      const int attr = mesh.GetAttribute(e);
      for (int i = 0; i < vdofs.Size(); ++i) {
        const int id = vdofs[i];
        if (id < 0)
          continue;
        if (attr == attr1_)
          C1(id) = C1n;
        if (attr == attr2_)
          C2(id) = C2n;
      }
    }
    // Rescale so FE ∫C1+∫C2 equals conserved M (H1 projection can inflate).
    {
      mfem::ConstantCoefficient one(1.0);
      auto feIntegral = [&](mfem::ParGridFunction &gf) {
        mfem::ParLinearForm mass(fes);
        mass.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
        mass.Assemble();
        return gf * mass;
      };
      const double d1 = feIntegral(C1);
      const double d2 = feIntegral(C2);
      const double dsum = d1 + d2;
      if (dsum > 0.0) {
        const double s = M / dsum;
        C1 *= s;
        C2 *= s;
      }
    }
    (void)C2_0;
  }
#endif

private:
  std::string sp1_, sp2_;
  int attr1_ = 1;
  int attr2_ = 2;
  SegregationCondition<NumericType> cond_;
#ifdef VIENNAPS_HAS_MFEM
  mutable std::unique_ptr<mfem::ConstantCoefficient> massCoef_;
#endif
};

} // namespace viennaps
