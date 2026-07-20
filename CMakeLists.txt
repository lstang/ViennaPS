cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(
  ViennaPS
  LANGUAGES CXX C
  VERSION 4.6.1)

# --------------------------------------------------------------------------------------------------------
# Include guards
# --------------------------------------------------------------------------------------------------------

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# --------------------------------------------------------------------------------------------------------
# Library switches
# --------------------------------------------------------------------------------------------------------

option(VIENNAPS_PRECOMPILE_HEADERS "Enable precompiled headers" OFF)
option(VIENNAPS_STATIC_BUILD "Build ViennaPS as static library" OFF)
option(VIENNAPS_ENABLE_SANITIZER "Enable sanitizers" OFF)
option(VIENNAPS_ENABLE_CLANG_TIDY "Run clang-tidy during build" OFF)
option(VIENNAPS_VTK_RENDERING "Enable VTK rendering support" ON)

option(VIENNAPS_USE_GPU "Enable GPU support" OFF)

option(VIENNAPS_BUILD_EXAMPLES "Build examples" OFF)
option(VIENNAPS_BUILD_TESTS "Build tests" OFF)

option(VIENNAPS_BUILD_PYTHON "Build python bindings" OFF)
option(VIENNAPS_PACKAGE_PYTHON "Build python bindings with intent to publish wheel" OFF)

# --------------------------------------------------------------------------------------------------------
# Library options
# --------------------------------------------------------------------------------------------------------

set(VIENNAPS_LOOKUP_DIRS
    ""
    CACHE STRING "Directories to account for when searching installed dependencies (VTK or embree)")

list(APPEND CMAKE_PREFIX_PATH ${VIENNAPS_LOOKUP_DIRS})

# Force skip vcpkg manifest install when using custom prebuilts (like your mfem in f:/dev/mfem/build)
# This prevents "vcpkg install failed" when baseline doesn't have all ports or for custom setups.
if(NOT DEFINED VCPKG_MANIFEST_INSTALL)
  set(VCPKG_MANIFEST_INSTALL OFF CACHE BOOL "Skip vcpkg manifest install for prebuilt/custom deps" FORCE)
endif()

# --------------------------------------------------------------------------------------------------------
# Early hints for custom / vcpkg dependencies (especially MFEM at user-specified location)
# --------------------------------------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH
     "f:/dev/mfem/build"
     "F:/dev/vcpkg/installed/x64-windows")

# Help transitive deps from VTK/HDF5 (zlib etc.) when using vcpkg toolchain
if(NOT ZLIB_ROOT)
  set(ZLIB_ROOT "F:/dev/vcpkg/installed/x64-windows" CACHE PATH "Hint for ZLIB used by VTK/HDF5")
endif()
if(NOT ZLIB_DIR)
  set(ZLIB_DIR "F:/dev/vcpkg/installed/x64-windows/share/zlib" CACHE PATH "")
endif()

# Force early discovery of custom MFEM (build tree)
if(EXISTS "f:/dev/mfem/build/MFEMConfig.cmake" OR EXISTS "F:/dev/mfem/build/MFEMConfig.cmake")
  set(mfem_DIR "f:/dev/mfem/build" CACHE PATH "MFEM build/install tree containing MFEMConfig.cmake" FORCE)
  set(MFEM_DIR "f:/dev/mfem/build" CACHE PATH "" FORCE)
endif()

# --------------------------------------------------------------------------------------------------------
# Global CMake Configuration
# └ We depend on the vtk python package for our python build and re-use
#   their existing vtk libraries, which are not versioned. Because we (possibly) pull in VTK ourselves
#   during build time though, we disable the soname versioning so that we can easily set the rpath to
#   link against the vtk-python libs
# --------------------------------------------------------------------------------------------------------

if(VIENNAPS_PACKAGE_PYTHON)
  message(
    STATUS
      "[ViennaPS] Using remote packages due to 'VIENNAPS_PACKAGE_PYTHON', the build will take a long time!"
  )

  # CIBUILDWHEEL is absolutely cursed and won't work when this is turned on.
  # I don't know why...
  if(NOT VIENNAPS_IS_CI)
    set(CPM_DOWNLOAD_ALL ON)
  endif()

  set(VIENNAPS_BUILD_PYTHON ON)
endif()

if(VIENNAPS_BUILD_PYTHON)
  set(CMAKE_PLATFORM_NO_VERSIONED_SONAME ON)
endif()

if(VIENNAPS_STATIC_BUILD)
  message(STATUS "[ViennaPS] Enabling precompiled headers for static build")
  set(VIENNAPS_PRECOMPILE_HEADERS ON)
  set(BUILD_SHARED_LIBS OFF)
endif()

if(MSVC)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
    set(OpenMP_RUNTIME_MSVC
        "llvm"
        CACHE STRING "MSVC OpenMP runtime selected for FindOpenMP" FORCE)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp:llvm")
  endif()

  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:/bigobj>")
endif()

# Reduce wheel size when packaging Python wheels.
# These flags are set before dependencies are added, so bundled dependencies use them too.
if(MSVC AND VIENNAPS_PACKAGE_PYTHON)
  string(APPEND CMAKE_C_FLAGS_RELEASE " /O2 /Gy /Gw")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " /O2 /Gy /Gw")
  string(APPEND CMAKE_MODULE_LINKER_FLAGS_RELEASE " /OPT:REF /OPT:ICF")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS_RELEASE " /OPT:REF /OPT:ICF")
  string(APPEND CMAKE_EXE_LINKER_FLAGS_RELEASE " /OPT:REF /OPT:ICF")
endif()

if((CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang") AND VIENNAPS_PACKAGE_PYTHON)
  string(APPEND CMAKE_C_FLAGS_RELEASE " -O3 -DNDEBUG -ffunction-sections -fdata-sections")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O3 -DNDEBUG -ffunction-sections -fdata-sections")
  if(APPLE)
    string(APPEND CMAKE_MODULE_LINKER_FLAGS_RELEASE " -Wl,-dead_strip")
    string(APPEND CMAKE_SHARED_LINKER_FLAGS_RELEASE " -Wl,-dead_strip")
  else()
    string(APPEND CMAKE_MODULE_LINKER_FLAGS_RELEASE " -Wl,--gc-sections")
    string(APPEND CMAKE_SHARED_LINKER_FLAGS_RELEASE " -Wl,--gc-sections")
  endif()
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Library
# --------------------------------------------------------------------------------------------------------

add_library(${PROJECT_NAME} INTERFACE)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES CXX_STANDARD 20
             CXX_EXTENSIONS OFF
             CXX_STANDARD_REQUIRED ON
             WINDOWS_EXPORT_ALL_SYMBOLS ON)

# --------------------------------------------------------------------------------------------------------
# Include directories
# --------------------------------------------------------------------------------------------------------

# Generate version header
configure_file("${PROJECT_SOURCE_DIR}/cmake/psVersion.hpp.in"
               "${PROJECT_SOURCE_DIR}/include/viennaps/psVersion.hpp" @ONLY)

target_include_directories(
  ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/viennaps>
                            $<INSTALL_INTERFACE:include/viennaps-${PROJECT_VERSION}>)

# --------------------------------------------------------------------------------------------------------
# Setup Dependencies
# --------------------------------------------------------------------------------------------------------
include("cmake/cpm.cmake")
include("cmake/vtk.cmake")

CPMAddPackage(
  NAME ViennaCore
  VERSION 2.2.1
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaCore"
  EXCLUDE_FROM_ALL ${VIENNAPS_BUILD_PYTHON}
  OPTIONS "VIENNACORE_FORMAT_EXCLUDE build/" "VIENNACORE_USE_GPU ${VIENNAPS_USE_GPU}")

CPMAddPackage(
  NAME PackageProject
  VERSION 1.13.0
  GIT_REPOSITORY "https://github.com/TheLartians/PackageProject.cmake"
  EXCLUDE_FROM_ALL ${VIENNAPS_BUILD_PYTHON})

CPMAddPackage(
  NAME ViennaRay
  VERSION 4.3.1
  GIT_TAG benchmark
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaRay"
  EXCLUDE_FROM_ALL ${VIENNAPS_BUILD_PYTHON}
  OPTIONS "VIENNARAY_USE_GPU ${VIENNAPS_USE_GPU}")

# --------------------------------------------------------------------------------------------------------
# Pre-resolve common transitive dependencies that VTK/HDF5 pull in (ZLIB etc.)
# This helps when vcpkg toolchain + manifest are used together with custom locations.
# --------------------------------------------------------------------------------------------------------
find_package(ZLIB QUIET)
if(NOT ZLIB_FOUND)
  set(ZLIB_ROOT "F:/dev/vcpkg/installed/x64-windows" CACHE PATH "" FORCE)
  set(ZLIB_INCLUDE_DIR "F:/dev/vcpkg/installed/x64-windows/include" CACHE PATH "" FORCE)
  if(EXISTS "F:/dev/vcpkg/installed/x64-windows/lib/zlib.lib")
    set(ZLIB_LIBRARY "F:/dev/vcpkg/installed/x64-windows/lib/zlib.lib" CACHE FILEPATH "" FORCE)
  elseif(EXISTS "F:/dev/vcpkg/installed/x64-windows/lib/zlibd.lib")
    set(ZLIB_LIBRARY "F:/dev/vcpkg/installed/x64-windows/lib/zlibd.lib" CACHE FILEPATH "" FORCE)
  endif()
  find_package(ZLIB QUIET)
endif()
if(ZLIB_FOUND)
  message(STATUS "[ViennaPS] ZLIB resolved: ${ZLIB_LIBRARIES}")
endif()

viennaps_prepare_vtk()

if(VIENNAPS_BUILD_PYTHON OR NOT VTK_FOUND)
  set(CMAKE_SKIP_INSTALL_RULES TRUE)
endif()

CPMAddPackage(
  NAME ViennaLS
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaLS"
  GIT_TAG v5.8.3
  EXCLUDE_FROM_ALL ${VIENNAPS_BUILD_PYTHON}
  OPTIONS "VIENNALS_PRECOMPILE_HEADERS ${VIENNAPS_PRECOMPILE_HEADERS}"
          "VIENNALS_VTK_RENDERING ${VIENNAPS_VTK_RENDERING}" "VIENNALS_USE_GPU ${VIENNAPS_USE_GPU}")

CPMAddPackage(
  NAME ViennaCS
  VERSION 2.0.1
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaCS"
  EXCLUDE_FROM_ALL ${VIENNAPS_BUILD_PYTHON})

target_link_libraries(${PROJECT_NAME} INTERFACE ViennaCore ViennaLS ViennaRay ViennaCS)

# --------------------------------------------------------------------------------------------------------
# Multipysics Dependencies (for unified field layer, high-fidelity diffusion/stress)
# --------------------------------------------------------------------------------------------------------
# MFEM for FEM discretization of bulk fields (dopants, defects, stress)
# We use a custom/prebuilt MFEM at f:/dev/mfem/build (not from vcpkg)
set(_mfem_candidate "f:/dev/mfem/build")
if(EXISTS "${_mfem_candidate}/MFEMConfig.cmake")
  set(mfem_DIR "${_mfem_candidate}" CACHE PATH "Location of MFEMConfig.cmake (build or install tree)" FORCE)
  set(MFEM_DIR "${_mfem_candidate}" CACHE PATH "" FORCE)
endif()
if(NOT mfem_DIR AND NOT MFEM_DIR)
  set(mfem_DIR "f:/dev/mfem/build" CACHE PATH "Location of MFEMConfig.cmake (build or install tree)")
endif()

# Prioritize the custom MFEM location
find_package(mfem CONFIG QUIET
             HINTS "${mfem_DIR}" "${MFEM_DIR}" "${_mfem_candidate}"
             PATHS "${mfem_DIR}" "${MFEM_DIR}" "${_mfem_candidate}"
             NO_DEFAULT_PATH)
if(mfem_FOUND)
  message(STATUS "[ViennaPS] Found MFEM ${MFEM_VERSION} for multiphysics fields at ${mfem_DIR}")
  message(STATUS "           MFEM_INCLUDE_DIRS=${MFEM_INCLUDE_DIRS}")

  # If MFEM was built with CUDA, try to locate toolkit but do not fail configure
  # (our current multiphysics code is host-only; full CUDA kernels later)
  if(MFEM_USE_CUDA)
    message(STATUS "[ViennaPS] MFEM built with CUDA support — locating CUDAToolkit (optional for now)")
    if(NOT DEFINED CMAKE_CUDA_COMPILER AND EXISTS "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.4/bin/nvcc.exe")
      set(CMAKE_CUDA_COMPILER "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.4/bin/nvcc.exe" CACHE FILEPATH "")
    endif()
    find_package(CUDAToolkit QUIET)
    if(CUDAToolkit_FOUND)
      message(STATUS "           CUDAToolkit version ${CUDAToolkit_VERSION}")
    else()
      message(WARNING "MFEM wants CUDA but CUDAToolkit not found — will link mfem but CUDA kernels disabled for now")
    endif()
    # Do NOT call enable_language(CUDA) here to avoid hard failure if no full CUDA setup for host compile.
  endif()

  if(TARGET mfem)
    target_link_libraries(${PROJECT_NAME} INTERFACE mfem)
  elseif(DEFINED MFEM_LIBRARY AND MFEM_LIBRARY)
    target_link_libraries(${PROJECT_NAME} INTERFACE ${MFEM_LIBRARY})
  else()
    # Fallback: some build trees expose it via imported target after include
    include("${mfem_DIR}/MFEMTargets.cmake" OPTIONAL)
    if(TARGET mfem)
      target_link_libraries(${PROJECT_NAME} INTERFACE mfem)
    endif()
  endif()

  if(MFEM_INCLUDE_DIRS)
    target_include_directories(${PROJECT_NAME} INTERFACE ${MFEM_INCLUDE_DIRS})
  elseif(MFEM_INCLUDE_DIR)
    target_include_directories(${PROJECT_NAME} INTERFACE ${MFEM_INCLUDE_DIR})
  endif()

  target_compile_definitions(${PROJECT_NAME} INTERFACE VIENNAPS_HAS_MFEM)

  # Debug vs Release CRT mismatch (MSVC MD/MDd) when linking a Release-prebuilt mfem.
  # Option: build a matching Debug mfem, or force Release for multiphysics targets.
  option(VIENNAPS_MFEM_REQUIRE_MATCHING_CRT
         "Warn if CMAKE_BUILD_TYPE/config may mismatch prebuilt MFEM CRT" ON)
  if(VIENNAPS_MFEM_REQUIRE_MATCHING_CRT AND MSVC)
    message(STATUS
      "[ViennaPS] MFEM linked from ${mfem_DIR}. "
      "Use matching CRT: Release app ↔ Release mfem (MD). "
      "Debug (MDd) against Release mfem is unsupported — build Debug mfem or use Release.")
    target_compile_definitions(${PROJECT_NAME} INTERFACE VIENNAPS_MFEM_CRT_NOTE=1)
  endif()
else()
  message(WARNING "[ViennaPS] MFEM not found (tried mfem_DIR=${mfem_DIR}). Multipysics features limited. Build MFEM or point -Dmfem_DIR=... to f:/dev/mfem/build")
endif()

# SUNDIALS for stiff time integration (diffusion, clustering, coupled systems)
find_package(SUNDIALS CONFIG QUIET)
if(SUNDIALS_FOUND)
  message(STATUS "[ViennaPS] Found SUNDIALS for time integration")
  target_link_libraries(${PROJECT_NAME} INTERFACE
    SUNDIALS::cvode
    SUNDIALS::ida
    SUNDIALS::nvecserial
    SUNDIALS::sunmatrixdense
    SUNDIALS::sunlinsoldense
    SUNDIALS::core)
  target_compile_definitions(${PROJECT_NAME} INTERFACE VIENNAPS_HAS_SUNDIALS)
else()
  message(WARNING "[ViennaPS] SUNDIALS not found - time stepping for stiff systems limited.")
endif()

# Eigen3 for auxiliary linear algebra (used by amgcl adapters etc.)
find_package(Eigen3 CONFIG QUIET)
if(Eigen3_FOUND)
  message(STATUS "[ViennaPS] Found Eigen3")
  target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
endif()

# amgcl (local in 3rdparty) for AMG solvers/preconditioners
# Use BUILD_INTERFACE to avoid "prefixed in the source directory" error for INTERFACE targets
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/3rdparty/amgcl>
)

# --------------------------------------------------------------------------------------------------------
# Setup Sanitizer or Clang-Tidy
# --------------------------------------------------------------------------------------------------------

if(VIENNAPS_ENABLE_SANITIZER)
  if(VIENNAPS_USE_GPU)
    message(WARNING "[ViennaPS] Sanitizers cannot be used with GPU support enabled!")
  else()
    viennacore_enable_sanitizer()
    message(STATUS "[ViennaPS] Enabled Sanitizers")
  endif()
endif()

if(VIENNAPS_ENABLE_CLANG_TIDY)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(NOT CLANG_TIDY_EXE)
    message(FATAL_ERROR "clang-tidy not found")
  endif()

  set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}" "--quiet" "--warnings-as-errors=*")
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Shared/Static Library
# --------------------------------------------------------------------------------------------------------

if(VIENNAPS_PRECOMPILE_HEADERS)
  set(VIENNAPS_LINKAGE SHARED)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

  if(VIENNAPS_STATIC_BUILD)
    set(VIENNAPS_LINKAGE STATIC)
  endif()

  set(LIB_NAME "viennaps")
  add_library(${LIB_NAME} ${VIENNAPS_LINKAGE})

  target_link_libraries(${LIB_NAME} PRIVATE ${PROJECT_NAME})
  target_compile_definitions(${LIB_NAME} PRIVATE VIENNAPS_USE_PRECOMPILED=1
                                                 VIENNACORE_COMPILE_SHARED_LIB=1)
  target_sources(
    ${LIB_NAME} PRIVATE "lib/specMain.cpp" "lib/specGeometries.cpp" "lib/specModelsEmulation.cpp"
                        "lib/specModelsSimulation.cpp" "lib/specProcess.cpp")

  set_target_properties(
    ${LIB_NAME}
    PROPERTIES VERSION ${PROJECT_VERSION} # full x.y.z appears in filename
               SOVERSION ${PROJECT_VERSION_MAJOR} # major version for SONAME
               WINDOWS_EXPORT_ALL_SYMBOLS ON)

  target_link_libraries(${PROJECT_NAME} INTERFACE ${LIB_NAME})
  target_compile_definitions(${PROJECT_NAME} INTERFACE VIENNAPS_USE_PRECOMPILED)

  install(
    TARGETS ${LIB_NAME}
    EXPORT ViennaPSTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup GPU
# --------------------------------------------------------------------------------------------------------

if(VIENNAPS_USE_GPU)
  message(STATUS "[ViennaPS] Enabled GPU Support")
  add_subdirectory(gpu)
else()
  message(STATUS "[ViennaPS] Disabled GPU Support")
  set(VIENNAPS_USE_GPU OFF)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Examples
# --------------------------------------------------------------------------------------------------------

include("cmake/ViennaPSUtils.cmake")

if(VIENNAPS_BUILD_EXAMPLES)
  message(STATUS "[ViennaPS] Building Examples")
  add_subdirectory(examples)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Tests
# --------------------------------------------------------------------------------------------------------

if(VIENNAPS_BUILD_TESTS)
  message(STATUS "[ViennaPS] Building Tests")

  enable_testing()
  add_subdirectory(tests)
endif()

# --------------------------------------------------------------------------------------------------------
# Setup Python Bindings
# --------------------------------------------------------------------------------------------------------

if(VIENNAPS_BUILD_PYTHON)
  message(STATUS "[ViennaPS] Building Python Bindings")
  add_subdirectory(python)
  return()
endif()

# --------------------------------------------------------------------------------------------------------
# Install Target
# --------------------------------------------------------------------------------------------------------

packageProject(
  NAME ${PROJECT_NAME} NAMESPACE ViennaTools
  VERSION ${PROJECT_VERSION}
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/viennaps
  INCLUDE_DESTINATION include/viennaps-${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  DEPENDENCIES "ViennaCore;ViennaLS;ViennaRay;ViennaCS")
