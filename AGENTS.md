# AGENTS.md

This is a **local fork** of [ViennaPS](https://github.com/ViennaTools/ViennaPS) (header-only C++ library for semiconductor process/topography simulation). It adds a multiphysics field layer (MFEM, SUNDIALS, amgcl) not present upstream.

## Machine-specific setup (critical)

The root `CMakeLists.txt` contains **hardcoded local paths** that must exist for a full build:
- `f:/dev/mfem/build` — prebuilt MFEM (provides `MFEMConfig.cmake`); enables `VIENNAPS_HAS_MFEM`
- `F:/dev/vcpkg/installed/x64-windows` — vcpkg installed deps (VTK, Embree, SUNDIALS, Eigen3, zlib)
- `3rdparty/amgcl/` — vendored locally; **gitignored**, must be present on disk

`VCPKG_MANIFEST_INSTALL` is forced `OFF` by default — deps are expected pre-installed, not built from manifest.

If MFEM/SUNDIALS are not found, configure still succeeds but multiphysics features are silently disabled (warnings printed). Check CMake status lines for `Found MFEM` / `Found SUNDIALS`.

MFEM CRT note (MSVC): link Release app with Release MFEM (MD). Debug (MDd) against Release MFEM is unsupported.

## Build commands (Windows / MSVC)

```powershell
# Configure (tests on)
cmake -B build -DVIENNAPS_BUILD_TESTS=ON -DCMAKE_TOOLCHAIN_FILE=F:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -DVCPKG_INSTALLED_DIR=F:/dev/vcpkg/installed -DVCPKG_MANIFEST_INSTALL=OFF

# Build
cmake --build build --config Release --parallel

# Examples / GPU / Python (optional)
cmake -B build -DVIENNAPS_BUILD_EXAMPLES=ON -DVIENNAPS_USE_GPU=ON -DVIENNAPS_BUILD_PYTHON=ON
```

Key CMake options: `VIENNAPS_BUILD_TESTS`, `VIENNAPS_BUILD_EXAMPLES`, `VIENNAPS_BUILD_PYTHON`, `VIENNAPS_USE_GPU`, `VIENNAPS_PRECOMPILE_HEADERS` (shared lib), `VIENNAPS_STATIC_BUILD`, `VIENNAPS_ENABLE_SANITIZER`, `VIENNAPS_ENABLE_CLANG_TIDY`.

## Tests

```powershell
ctest -E "Benchmark|Performance" --test-dir build -C Release --output-on-failure
```

Tests live in `tests/<name>/` — each is a single `.cpp` linked against `ViennaPS` and registered via `add_test`. Benchmarks are excluded by default. To run one test: `ctest -R <testName> --test-dir build`.

## Format (required before PRs)

```powershell
cmake -B build
cmake --build build --target format       # apply
cmake --build build --target format-check # CI checks
```

Style: LLVM-based clang-format — **2-space indent, 80-col limit, no tabs, Attach braces, pointer right-aligned**. See `.clang-format`.

## Architecture

- **Header-only**: all C++ code is in `include/viennaps/`. The `lib/` dir contains precompiled template specializations (only built when `VIENNAPS_PRECOMPILE_HEADERS=ON`).
- **`include/viennaps/models/`** — physics process models (etching, deposition, oxidation, etc.)
- **`include/viennaps/fields/`** — multiphysics field layer (fork addition): `PhysicsField`, `MfemElasticityKernel`, `SundialsTimeIntegrator`, `DiffusionKernel`, `StressKernel`, `AmgclSolver`, etc. Gated by `VIENNAPS_HAS_MFEM` / `VIENNAPS_HAS_SUNDIALS`.
- **`include/viennaps/process/`** — process execution: `psProcess`, flux engines, process strategies, GPU/CPU engines.
- **`gpu/`** — OptiX/CUDA kernels for ray tracing (experimental, requires `VIENNAPS_USE_GPU=ON`).
- **`python/`** — pybind11 bindings; built via scikit-build-core (`pip install .`).
- ViennaTools deps (ViennaCore, ViennaLS, ViennaRay, ViennaCS) are fetched via CPM.cmake at configure time.

## Python bindings

```powershell
pip install .                                    # basic
python python/scripts/install_ViennaPS.py        # with GPU support
```

Python module is `viennaps` (import as `import viennaps as vps`). Dimension defaults to 2D; set via `vps.setDimension(3)`.

## Conventions

- C++20 required.
- Follow LLVM coding guidelines.
- Include order is sorted case-sensitively by clang-format.
- `viennaps.hpp` is the umbrella header — new public headers should be added there.
