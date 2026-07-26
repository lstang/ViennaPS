#pragma once

/// SolutionTransfer — field transfer between meshes with dose preservation.
///
/// Paths (MFEM):
///  A) L2 projection (preferred): build mass M on the target space, RHS
///     b_i = ∫ source(x) φ_i dV via FindPoints-sampled coefficient, solve
///     M x = b (MultiAppProjectionTransfer::assembleL2 pattern).
///  B) Pointwise ProjectCoefficient fallback + integral-preserving rescale.
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
  int unmappedQuadraturePoints = 0;
  bool ok = false;
};

class SolutionTransfer {
public:
#ifdef VIENNAPS_HAS_MFEM
  /// L2 dual pairing with the constant-1 functional: ∫ C dV.
  static double integrate(const mfem::ParGridFunction &gf) {
    if (gf.Size() == 0 || !gf.ParFESpace())
      return 0.0;
    mfem::ConstantCoefficient one(1.0);
    mfem::ParLinearForm mass(gf.ParFESpace());
    mass.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    mass.Assemble();
    return gf * mass;
  }

  static double meshVolume(mfem::ParMesh &mesh) {
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

  /// Sample source GF at a physical point (FindPoints).
  class SourceSampleCoef : public mfem::Coefficient {
  public:
    SourceSampleCoef(const mfem::ParGridFunction *src, mfem::ParMesh *smesh)
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
      if (el[0] < 0) {
        ++misses_;
        return 0.0;
      }
      return src_->GetValue(*smesh_->GetElementTransformation(el[0]), ips[0]);
    }
    int misses() const { return misses_; }

  private:
    const mfem::ParGridFunction *src_;
    mfem::ParMesh *smesh_;
    mutable int misses_ = 0;
  };

  /// L2 projection: solve M x = b with b_i = ∫ s(x) φ_i, then dose rescale.
  static TransferResult
  transferL2(const mfem::ParGridFunction &source, mfem::ParGridFunction &target) {
    TransferResult r;
    r.doseSource = integrate(source);
    if (target.Size() == 0 || !target.ParFESpace() || !source.ParFESpace()) {
      r.ok = false;
      return r;
    }
    mfem::ParMesh *smesh = source.ParFESpace()->GetParMesh();
    mfem::ParFiniteElementSpace *fes = target.ParFESpace();
    if (!smesh || !fes) {
      r.ok = false;
      return r;
    }

    SourceSampleCoef coef(&source, smesh);
    mfem::ConstantCoefficient one(1.0);
    mfem::ParBilinearForm M(fes);
    M.AddDomainIntegrator(new mfem::MassIntegrator(one));
    M.Assemble();
    M.Finalize();

    mfem::ParLinearForm b(fes);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));
    b.Assemble();
    r.unmappedQuadraturePoints = coef.misses();

    mfem::HypreParMatrix *Mh = M.ParallelAssemble();
    mfem::HypreBoomerAMG amg(*Mh);
    amg.SetPrintLevel(0);
    target = 0.0;
    mfem::HyprePCG pcg(*Mh);
    pcg.SetPreconditioner(amg);
    pcg.SetTol(1e-12);
    pcg.SetMaxIter(400);
    pcg.SetPrintLevel(0);
    pcg.Mult(b, target);
    delete Mh;

    r.doseTargetBeforeScale = integrate(target);
    if (r.doseSource > 0.0) {
      const double now = r.doseTargetBeforeScale;
      if (now > 0.0)
        target *= static_cast<mfem::real_t>(r.doseSource / now);
      else {
        const double vol = meshVolume(*fes->GetParMesh());
        if (vol > 0.0)
          target = static_cast<mfem::real_t>(r.doseSource / vol);
      }
    }
    r.doseTargetAfterScale = integrate(target);
    const double denom = std::max(std::abs(r.doseSource), 1e-30);
    r.relativeDoseError =
        std::abs(r.doseTargetAfterScale - r.doseSource) / denom;
    r.ok = (r.relativeDoseError <= 1e-3);
    return r;
  }

  /// Default transfer: L2 Mx=b path with integral-preserving rescale.
  static TransferResult
  transferIntegralPreserving(const mfem::ParGridFunction &source,
                             mfem::ParGridFunction &target) {
    return transferL2(source, target);
  }

  /// Same-size dof copy then optional rescale (identical FE spaces).
  static TransferResult
  copyIntegralPreserving(const mfem::ParGridFunction &source,
                         mfem::ParGridFunction &target) {
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
