#pragma once

/// PhysicsKernel base for the unified multiphysics framework.
/// Inspired by MOOSE kernels: each kernel contributes residual / Jacobian terms.
///
/// Two RHS styles:
///   1. addToRHS(t, y, ydot, myIndex)  — legacy per-kernel scalar state
///   2. addToFieldRHS(t, y, ydot)      — full packed multi-species profile state
///
/// SundialsTimeIntegrator prefers (2) when the field has packable species state.

#include <string>
#include <memory>
#include <vector>

namespace viennaps {

template <class NumericType>
class PhysicsField;

template <class NumericType>
class MaterialPropertySystem;

template <class NumericType>
class PhysicsKernel {
public:
  PhysicsKernel() = default;
  virtual ~PhysicsKernel() = default;

  virtual void setName(const std::string& n) { name_ = n; }
  const std::string& getName() const { return name_; }

  virtual void setPhysicsField(std::shared_ptr<PhysicsField<NumericType>> field) {
    field_ = field;
  }
  virtual void setMaterialProperties(std::shared_ptr<MaterialPropertySystem<NumericType>> mat) {
    material_ = mat;
  }

  virtual void setup() {}
  virtual void computeResidual() {}
  virtual void computeJacobian() {}

  virtual void evolve(NumericType /*dt*/) {}

  /// Legacy: one state entry per kernel.
  virtual void addToRHS(NumericType /*t*/, const std::vector<NumericType>& /*y*/,
                        std::vector<NumericType>& /*ydot*/, int /*myIndex*/) {}

  /// Field-packed residual: y / ydot are concatenations of all species profiles.
  /// Default no-op; diffusion / reaction kernels override.
  virtual void addToFieldRHS(NumericType /*t*/, const std::vector<NumericType>& /*y*/,
                             std::vector<NumericType>& /*ydot*/) {}

  virtual void onGeometryUpdate() {}

protected:
  std::string name_;
  std::shared_ptr<PhysicsField<NumericType>> field_;
  std::shared_ptr<MaterialPropertySystem<NumericType>> material_;
};

} // namespace viennaps
