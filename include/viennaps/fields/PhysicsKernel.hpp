#pragma once

/// PhysicsKernel base for the unified multiphysics framework.
/// Inspired by MOOSE kernels: each kernel contributes to the residual/Jacobian
/// of a coupled system (dopant diffusion, defect reactions, stress, etc.).
///
/// In full implementation:
/// - Kernels assemble into MFEM weak forms (Bilinear/Linear forms)
/// - SUNDIALS (CVODE/IDA) drives time integration of the stiff nonlinear system
/// - amgcl provides AMG preconditioned solves
///
/// Currently provides a working stub interface for early integration with
/// AnalyticImplant + BasicDiffusion and the Oxidation adapter.

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

  // Associate with the unified field and material property system
  virtual void setPhysicsField(std::shared_ptr<PhysicsField<NumericType>> field) {
    field_ = field;
  }
  virtual void setMaterialProperties(std::shared_ptr<MaterialPropertySystem<NumericType>> mat) {
    material_ = mat;
  }

  // Called once before time stepping (allocate temporaries, setup BCs, etc.)
  virtual void setup() {}

  // Contribute to residual (for nonlinear solve / time integrator)
  virtual void computeResidual() {}

  // Contribute to Jacobian (for implicit methods)
  virtual void computeJacobian() {}

  // Simple explicit-style evolve for early demos and stub models
  // Full version will be driven by SUNDIALS through the kernel assembly
  virtual void evolve(NumericType dt) {
    // Derived kernels override to advance their species/reactions
  }

  // Contribute to RHS for SUNDIALS/CVODE: ydot[myIndex] += f(t, y) from this kernel
  // In the integrator, y is a flat vector (one entry per active kernel for demo).
  // Real version: kernels would read/write MFEM dofs for their species.
  virtual void addToRHS(NumericType t, const std::vector<NumericType>& y, std::vector<NumericType>& ydot, int myIndex) {
    // Default: no contribution. Override in concrete kernels.
  }

  // Hook for geometry change (level-set moved) — remap or invalidate cached data
  virtual void onGeometryUpdate() {}

protected:
  std::string name_;
  std::shared_ptr<PhysicsField<NumericType>> field_;
  std::shared_ptr<MaterialPropertySystem<NumericType>> material_;
};

} // namespace viennaps