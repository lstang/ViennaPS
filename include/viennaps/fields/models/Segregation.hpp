#pragma once

/// Segregation — two-sided Si/SiO2 (or general mat1/mat2) interface condition.
///
/// rate = kf * C_1 - kb * C_2
/// equilibrium: C_2 / C_1 = kf / kb = m(T)
///
/// Mirrors MOOSE InterfaceReaction (framework/interfacekernels/InterfaceReaction):
///   Element residual:   +test      * rate
///   Neighbor residual:  -test_nbr  * rate
/// with the full 4-block Jacobian.
///
/// Implementation option (Phase 2 Task 5): **two integrators** — one per
/// material side on the shared interior interface faces. Dose conservation
/// across the interface is the acceptance criterion.

#include "../DiffusionModel.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace viennaps {

/// Standalone two-sided segregation condition (usable without DiffusionModel).
template <class NumericType>
class SegregationCondition {
public:
  void setSegregationCoefficient(NumericType m, NumericType k0) {
    // kf = m * k0 (mat1 → mat2), kb = k0 (mat2 → mat1)
    // equilibrium: C2/C1 = kf/kb = m
    m_ = m;
    k0_ = k0;
    kf_ = m * k0;
    kb_ = k0;
  }

  NumericType segregationCoefficient() const { return m_; }
  NumericType kf() const { return kf_; }
  NumericType kb() const { return kb_; }

  /// Equilibrium interface ratio C2/C1.
  NumericType equilibriumRatio() const {
    return (kb_ != NumericType(0)) ? (kf_ / kb_) : m_;
  }

  /// Instantaneous transfer rate (positive = flux mat1 → mat2).
  NumericType rate(NumericType C1, NumericType C2) const {
    return kf_ * C1 - kb_ * C2;
  }

#ifdef VIENNAPS_HAS_MFEM
  /// Element (mat1) residual contribution at a quadrature point:
  ///   R_elem += w * test * rate
  /// and 2 of the 4 Jacobian blocks (EE, EN).
  void assembleElementSide(
      mfem::Vector &R_elem, mfem::DenseMatrix &K_ee, mfem::DenseMatrix &K_en,
      const mfem::Vector &C1_dofs, const mfem::Vector &C2_dofs,
      const mfem::Vector &shape_e, const mfem::Vector &shape_n,
      double weight) const {
    const double C1 = mfem::InnerProduct(C1_dofs, shape_e);
    const double C2 = mfem::InnerProduct(C2_dofs, shape_n);
    const double r = static_cast<double>(kf_) * C1 - static_cast<double>(kb_) * C2;
    // R_elem_i += w * phi_i * r
    for (int i = 0; i < shape_e.Size(); ++i) {
      R_elem(i) += weight * shape_e(i) * r;
      for (int j = 0; j < shape_e.Size(); ++j) {
        K_ee(i, j) += weight * shape_e(i) * static_cast<double>(kf_) * shape_e(j);
      }
      for (int j = 0; j < shape_n.Size(); ++j) {
        K_en(i, j) +=
            weight * shape_e(i) * (-static_cast<double>(kb_)) * shape_n(j);
      }
    }
  }

  /// Neighbor (mat2) residual — SIGN FLIP of element residual.
  void assembleNeighborSide(
      mfem::Vector &R_nbr, mfem::DenseMatrix &K_nn, mfem::DenseMatrix &K_ne,
      const mfem::Vector &C1_dofs, const mfem::Vector &C2_dofs,
      const mfem::Vector &shape_e, const mfem::Vector &shape_n,
      double weight) const {
    const double C1 = mfem::InnerProduct(C1_dofs, shape_e);
    const double C2 = mfem::InnerProduct(C2_dofs, shape_n);
    const double r = static_cast<double>(kf_) * C1 - static_cast<double>(kb_) * C2;
    // R_nbr_i += -w * phi_n_i * r
    for (int i = 0; i < shape_n.Size(); ++i) {
      R_nbr(i) += -weight * shape_n(i) * r;
      for (int j = 0; j < shape_n.Size(); ++j) {
        // dR_nbr/dC2 = -test_n * (-kb) * phi_n = +kb * test_n * phi_n
        // MOOSE: NeighborNeighbor = -test_n * (-kb) * phi_n = + kb * ...
        K_nn(i, j) +=
            weight * shape_n(i) * static_cast<double>(kb_) * shape_n(j);
      }
      for (int j = 0; j < shape_e.Size(); ++j) {
        // dR_nbr/dC1 = -test_n * kf * phi_e
        K_ne(i, j) +=
            -weight * shape_n(i) * static_cast<double>(kf_) * shape_e(j);
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

/// DiffusionModel wrapper that exposes segregation as a named model.
/// Full mesh wiring (interior face assembly) is done by the engine or by
/// the Phase 2 segregation integration test via SegregationCondition
/// directly.
template <class NumericType>
class Segregation : public DiffusionModel<NumericType> {
public:
  explicit Segregation(const std::string &species = "Boron",
                       int attrMat1 = 1, int attrMat2 = 2) {
    this->setName("Segregation(" + species + ")");
    species_ = species;
    attr1_ = attrMat1;
    attr2_ = attrMat2;
  }

  void setSegregationCoefficient(NumericType m, NumericType k0) {
    cond_.setSegregationCoefficient(m, k0);
  }

  SegregationCondition<NumericType> &condition() { return cond_; }
  const SegregationCondition<NumericType> &condition() const { return cond_; }

  int numSpecies() const override { return 1; }
  std::vector<std::string> speciesNames() const override {
    return {species_};
  }
  std::vector<int> applicableAttributes() const override {
    return {attr1_, attr2_};
  }

private:
  std::string species_;
  int attr1_ = 1;
  int attr2_ = 2;
  SegregationCondition<NumericType> cond_;
};

} // namespace viennaps
