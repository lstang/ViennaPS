# --------------------------------------------------------------------------------------------------------
# Dependency source patches.
#
# ViennaRay 4.5.0 applies a `reduction` clause to the *split* `#pragma omp parallel` region in
# include/viennaray/rayTraceKernel.hpp (the loop is a bare `#pragma omp for` inside it, not a fused
# `parallel for`). With MSVC's LLVM OpenMP runtime (`/openmp:llvm` -> libomp140) the reduction
# finalization then mis-identifies every thread of the team, so the process prints
#   OMP: Error #132: Thread identifier invalid.
# and aborts. On this fork that kills six tests, all of which run that kernel:
#   cpuDiskEngine, fluxEngines, multiParticleProcess, process, rngSeed, singleParticleProcess
#
# The patch replaces the reduction with per-thread counters that are aggregated after the region
# (identical semantics, no atomics). Verified against ViennaRay 4.5.0: fluxEngines exits 3 before,
# 0 after; full ctest suite 44/44 with it.
#
# Remove this call (and cmake/patches/viennaray-4.5.0-omp-reduction.patch) once ViennaRay carries
# the fix upstream.
# --------------------------------------------------------------------------------------------------------
function(viennaps_patch_viennaray_openmp_reduction VIENNARAY_SOURCE_DIR)
  set(_viennaray_kernel
      "${VIENNARAY_SOURCE_DIR}/include/viennaray/rayTraceKernel.hpp")

  if(NOT EXISTS "${_viennaray_kernel}")
    message(WARNING "[ViennaPS] ViennaRay OpenMP patch skipped: ${_viennaray_kernel} not found")
    return()
  endif()

  file(READ "${_viennaray_kernel}" _viennaray_contents)
  string(FIND "${_viennaray_contents}" "threadCounters" _already_patched)
  if(NOT _already_patched EQUAL -1)
    message(STATUS "[ViennaPS] ViennaRay OpenMP reduction patch already applied")
    return()
  endif()

  find_package(Git QUIET)
  if(NOT GIT_EXECUTABLE)
    message(WARNING "[ViennaPS] ViennaRay OpenMP patch skipped: git executable not found")
    return()
  endif()

  set(_viennaray_patch
      "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/patches/viennaray-4.5.0-omp-reduction.patch")
  message(STATUS "[ViennaPS] Patching ViennaRay OpenMP reduction region")

  execute_process(
    COMMAND "${GIT_EXECUTABLE}" apply --ignore-whitespace "${_viennaray_patch}"
    WORKING_DIRECTORY "${VIENNARAY_SOURCE_DIR}"
    RESULT_VARIABLE _viennaray_apply_result
    ERROR_VARIABLE _viennaray_apply_error)

  if(_viennaray_apply_result EQUAL 0)
    message(
      STATUS
        "[ViennaPS] Applied ViennaRay OpenMP reduction patch (MSVC libomp 'OMP: Error #132')")
  else()
    message(
      WARNING
        "[ViennaPS] Could not apply the ViennaRay OpenMP reduction patch.\n"
        "          Without it, MSVC builds abort six ray-tracing tests with 'OMP: Error #132: Thread identifier invalid'.\n"
        "          Patch: ${_viennaray_patch}\n"
        "          ViennaRay source may have changed - update the patch. git said: ${_viennaray_apply_error}")
  endif()
endfunction()
