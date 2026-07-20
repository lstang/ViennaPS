#pragma once

/// AmgclSolver - Thin wrapper around amgcl for sparse CSR systems.
/// Used for production-scale linear solves (diffusion residual, elasticity).
/// Header-only; requires 3rdparty/amgcl on include path (CMake already sets it).

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

// amgcl's make_solver path pulls boost::property_tree. Enable only when both exist.
#if defined(__has_include)
#if __has_include(<boost/property_tree/ptree.hpp>) && __has_include(<amgcl/make_solver.hpp>)
#define VIENNAPS_HAS_AMGCL 1
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/backend/builtin.hpp>
#else
#define VIENNAPS_HAS_AMGCL 0
#endif
#else
#define VIENNAPS_HAS_AMGCL 0
#endif

namespace viennaps {

/// Simple CSR matrix (0-based).
struct CsrMatrix {
  std::vector<int> ptr;    // size n+1
  std::vector<int> col;    // size nnz
  std::vector<double> val; // size nnz
  int n = 0;
};

/// Build a 1D Laplacian CSR of size n with Dirichlet ends (for smoke tests).
inline CsrMatrix makeLaplacian1D(int n, double diagScale = 2.0) {
  CsrMatrix A;
  A.n = n;
  A.ptr.reserve(n + 1);
  A.ptr.push_back(0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      A.col.push_back(i - 1);
      A.val.push_back(-1.0);
    }
    A.col.push_back(i);
    A.val.push_back(diagScale);
    if (i + 1 < n) {
      A.col.push_back(i + 1);
      A.val.push_back(-1.0);
    }
    A.ptr.push_back(static_cast<int>(A.col.size()));
  }
  return A;
}

struct AmgclSolveResult {
  bool ok = false;
  int iters = 0;
  double residual = 0;
  bool usedAmgcl = false;
};

/// Solve A x = b. Uses amgcl when available; otherwise Jacobi-preconditioned CG.
inline AmgclSolveResult solveCSR(const CsrMatrix& A, const std::vector<double>& b,
                                 std::vector<double>& x, double tol = 1e-8,
                                 int maxIters = 500) {
  AmgclSolveResult res;
  if (A.n <= 0 || static_cast<int>(b.size()) != A.n) return res;
  x.assign(A.n, 0.0);

#if VIENNAPS_HAS_AMGCL
  try {
    typedef amgcl::backend::builtin<double> Backend;
    typedef amgcl::make_solver<
        amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation,
                   amgcl::relaxation::spai0>,
        amgcl::solver::bicgstab<Backend>>
        Solver;

    auto backend_prm = Backend::params();
    Solver::params prm;
    prm.solver.tol = tol;
    prm.solver.maxiter = maxIters;

    auto sys = std::make_tuple(A.n, A.ptr, A.col, A.val);
    Solver solve(sys, prm, backend_prm);

    std::vector<double> rhs = b;
    auto [iters, error] = solve(rhs, x);
    res.ok = true;
    res.iters = static_cast<int>(iters);
    res.residual = error;
    res.usedAmgcl = true;
    std::cout << "[AmgclSolver] amgcl BiCGSTAB iters=" << res.iters
              << " residual=" << res.residual << " n=" << A.n << "\n";
    return res;
  } catch (const std::exception& ex) {
    std::cout << "[AmgclSolver] amgcl failed (" << ex.what()
              << "), falling back to Jacobi-CG\n";
  }
#endif

  // Jacobi-preconditioned CG fallback
  auto matvec = [&](const std::vector<double>& v, std::vector<double>& Av) {
    Av.assign(A.n, 0.0);
    for (int i = 0; i < A.n; ++i) {
      for (int j = A.ptr[i]; j < A.ptr[i + 1]; ++j)
        Av[i] += A.val[j] * v[A.col[j]];
    }
  };
  std::vector<double> diag(A.n, 1.0);
  for (int i = 0; i < A.n; ++i) {
    for (int j = A.ptr[i]; j < A.ptr[i + 1]; ++j)
      if (A.col[j] == i) diag[i] = A.val[j];
    if (std::abs(diag[i]) < 1e-30) diag[i] = 1.0;
  }

  std::vector<double> r(A.n), z(A.n), p(A.n), Ap(A.n);
  matvec(x, Ap);
  for (int i = 0; i < A.n; ++i) r[i] = b[i] - Ap[i];
  for (int i = 0; i < A.n; ++i) z[i] = r[i] / diag[i];
  p = z;
  double rz = 0;
  for (int i = 0; i < A.n; ++i) rz += r[i] * z[i];
  double bnorm = 0;
  for (double bi : b) bnorm += bi * bi;
  bnorm = std::sqrt(std::max(bnorm, 1e-30));

  int it = 0;
  double rel = 1.0;
  for (; it < maxIters; ++it) {
    matvec(p, Ap);
    double pAp = 0;
    for (int i = 0; i < A.n; ++i) pAp += p[i] * Ap[i];
    if (std::abs(pAp) < 1e-30) break;
    double alpha = rz / pAp;
    for (int i = 0; i < A.n; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
    }
    double rnorm = 0;
    for (double ri : r) rnorm += ri * ri;
    rel = std::sqrt(rnorm) / bnorm;
    if (rel < tol) break;
    for (int i = 0; i < A.n; ++i) z[i] = r[i] / diag[i];
    double rzNew = 0;
    for (int i = 0; i < A.n; ++i) rzNew += r[i] * z[i];
    double beta = rzNew / std::max(rz, 1e-30);
    for (int i = 0; i < A.n; ++i) p[i] = z[i] + beta * p[i];
    rz = rzNew;
  }
  res.ok = rel < tol * 10 || it > 0;
  res.iters = it;
  res.residual = rel;
  res.usedAmgcl = false;
  std::cout << "[AmgclSolver] Jacobi-CG iters=" << res.iters
            << " residual=" << res.residual << " n=" << A.n << "\n";
  return res;
}

/// Convenience: solve diffusion mass+stiffness style system on packed profile.
template <class NumericType>
inline AmgclSolveResult solveProfileDiffusion(std::vector<NumericType>& profile,
                                              NumericType Ddt, double tol = 1e-8) {
  const int n = static_cast<int>(profile.size());
  if (n < 2) return {};
  // (M + dt D K) c_new = M c_old  with lumped mass, 1D K Laplacian
  CsrMatrix A = makeLaplacian1D(n, 2.0);
  // Scale off-diagonals by -Ddt, diagonal by 1+2*Ddt
  double a = static_cast<double>(std::max(NumericType(0), Ddt));
  for (int i = 0; i < n; ++i) {
    for (int j = A.ptr[i]; j < A.ptr[i + 1]; ++j) {
      if (A.col[j] == i)
        A.val[j] = 1.0 + 2.0 * a;
      else
        A.val[j] = -a;
    }
  }
  std::vector<double> b(n), x(n, 0.0);
  for (int i = 0; i < n; ++i) b[i] = static_cast<double>(profile[i]);
  auto res = solveCSR(A, b, x, tol, 400);
  if (res.ok) {
    for (int i = 0; i < n; ++i)
      profile[i] = static_cast<NumericType>(std::max(0.0, x[i]));
  }
  return res;
}

} // namespace viennaps
