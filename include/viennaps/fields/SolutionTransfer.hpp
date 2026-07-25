#pragma once

/// SolutionTransfer — field transfer between meshes with dose preservation.
///
/// Path (MFEM):
///  1. Project source field onto target FE space by evaluating the source
///     GridFunction at physical points of the target mesh (FindPoints).
///  2. Rescale so ∫C_target = ∫C_source (IntegralPreservingFunctionIC pattern;
///     same mass dual as TedInitializer::projectToGridFunction).
///
/// Dose = (1, C)_L2 via DomainLFIntegrator(ConstantCoefficient(1)).

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef VIENNAPS_HAS_MFEM
#include <mfem.hpp>
#endif

namespace viennaps {

struct TransferResult {
  double doseSource = 0.0;
  double doseTargetBeforeScale = 0.0;
  double doseTargetAfterScale = 0.0;
  double relativeDoseError = 0.0;
  bool ok = false;
};

class SolutionTransfer {
public:
#ifdef VIENNAPS_HAS_MFEM
  /// L2 dual pairing with the constant-1 functional: ∫ C dV.
  static double integrate(const mfem::GridFunction &gf) {
    if (gf.Size() == 0 || !gf.FESpace())
      return 0.0;
    mfem::ConstantCoefficient one(1.0);
    mfem::LinearForm mass(const_cast<mfem::FiniteElementSpace *>(gf.FESpace()));
    mass.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    mass.Assemble();
    return gf * mass;
  }

  static double meshVolume(mfem::Mesh &mesh) {
    double vol = 0.0;
    for (int e = 0; e < mesh.GetNE(); ++e) {
      mfem::ElementTransformation *T = mesh.GetElementTransformation(e);
      const mfem::IntegrationRule *ir =
          &mfem::IntRules.Get(mesh.GetElementBaseGeometry(e), 2);
      for (int i = 0; i < ir->GetNPoints(); ++i) {
        const mfem::IntegrationPoint &ip = ir->IntPoint(i);
        T->SetIntPoint(&ip);
        vol += ip.weight * T->Weight();
      }
    }
    return vol;
  }

  /// Transfer source → target (possibly different meshes). After sampling,
  /// rescale so dose matches source when source dose > 0.
  static TransferResult
  transferIntegralPreserving(const mfem::GridFunction &source,
                             mfem::GridFunction &target) {
    TransferResult r;
    r.doseSource = integrate(source);
    if (target.Size() == 0 || !target.FESpace() || !source.FESpace()) {
      r.ok = false;
      return r;
    }
    mfem::Mesh *tmesh = target.FESpace()->GetMesh();
    mfem::Mesh *smesh = source.FESpace()->GetMesh();
    if (!tmesh || !smesh) {
      r.ok = false;
      return r;
    }

    // Coefficient: at each target QP, sample source via FindPoints.
    class SourceSampleCoef : public mfem::Coefficient {
    public:
      SourceSampleCoef(const mfem::GridFunction *src, mfem::Mesh *smesh)
          : src_(src), smesh_(smesh) {}
      double Eval(mfem::ElementTransformation &T,
                  const mfem::IntegrationPoint &ip) override {
        mfem::Vector x;
        T.Transform(ip, x);
        mfem::DenseMatrix pts(x.Size(), 1);
        for (int d = 0; d < x.Size(); ++d)
          pts(d, 0) = x(d);
        mfem::Array<int> el(1);
        mfem::Array<mfem::IntegrationPoint> ips(1);
        smesh_->FindPoints(pts, el, ips, /*warn=*/false);
        if (el[0] < 0)
          return 0.0;
        return src_->GetValue(*smesh_->GetElementTransformation(el[0]), ips[0]);
      }

    private:
      const mfem::GridFunction *src_;
      mfem::Mesh *smesh_;
    };

    SourceSampleCoef coef(&source, smesh);
    target.ProjectCoefficient(coef);

    r.doseTargetBeforeScale = integrate(target);
    if (r.doseSource > 0.0) {
      const double now = r.doseTargetBeforeScale;
      if (now > 0.0) {
        target *= static_cast<mfem::real_t>(r.doseSource / now);
      } else {
        const double vol = meshVolume(*tmesh);
        if (vol > 0.0)
          target = static_cast<mfem::real_t>(r.doseSource / vol);
      }
    }
    r.doseTargetAfterScale = integrate(target);
    const double denom = std::max(std::abs(r.doseSource), 1e-30);
    r.relativeDoseError =
        std::abs(r.doseTargetAfterScale - r.doseSource) / denom;
    r.ok = (r.relativeDoseError <= 1e-3); // 0.1%
    return r;
  }

  /// Same-size dof copy then optional rescale (identical FE spaces).
  static TransferResult
  copyIntegralPreserving(const mfem::GridFunction &source,
                         mfem::GridFunction &target) {
    if (target.Size() != source.Size())
      return transferIntegralPreserving(source, target);
    TransferResult r;
    r.doseSource = integrate(source);
    target = source;
    r.doseTargetBeforeScale = integrate(target);
    r.doseTargetAfterScale = r.doseTargetBeforeScale;
    const double denom = std::max(std::abs(r.doseSource), 1e-30);
    r.relativeDoseError =
        std::abs(r.doseTargetAfterScale - r.doseSource) / denom;
    r.ok = (r.relativeDoseError <= 1e-3);
    return r;
  }
#endif
};

} // namespace viennaps
