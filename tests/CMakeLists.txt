# Standalone or subdirectory unit-test build for GCTA helpers (no MKL / no gcta64).
# Standalone:
#   cmake -S tests -B build-unit && cmake --build build-unit && ctest --test-dir build-unit --output-on-failure
#
# From repo root (after enable_testing + add_subdirectory(tests)):
#   uses parent Eigen3::Eigen when available.

cmake_minimum_required(VERSION 3.16)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(gcta_cojo_unit_tests LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 17)
    set(GCTA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
    set(GCTA_TESTS_STANDALONE ON)
else()
    set(GCTA_ROOT "${CMAKE_SOURCE_DIR}")
    set(GCTA_TESTS_STANDALONE OFF)
endif()

include(CTest)
enable_testing()

# GoogleTest via CPM when standalone; reuse parent CPM if present.
if(GCTA_TESTS_STANDALONE)
    list(APPEND CMAKE_MODULE_PATH "${GCTA_ROOT}/cmake")
    if(EXISTS "${GCTA_ROOT}/cmake/CPM.cmake")
        include("${GCTA_ROOT}/cmake/CPM.cmake")
    else()
        include(FetchContent)
        function(CPMAddPackage)
            # minimal fallback unused
        endfunction()
    endif()
    set(FETCHCONTENT_QUIET OFF)
    CPMAddPackage("gh:eigen-mirror/eigen#3.4.0")
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION 1.15.2
        OPTIONS "INSTALL_GTEST OFF" "BUILD_GMOCK OFF")
else()
    include("${GCTA_ROOT}/cmake/CPM.cmake")
    CPMAddPackage(
        NAME googletest
        GITHUB_REPOSITORY google/googletest
        VERSION 1.15.2
        OPTIONS "INSTALL_GTEST OFF" "BUILD_GMOCK OFF")
endif()

add_library(gcta_test_support STATIC
    "${GCTA_ROOT}/main/CommFunc.cpp"
    "${GCTA_ROOT}/main/StrFunc.cpp"
    "${GCTA_ROOT}/main/StatFunc.cpp"
    "${GCTA_ROOT}/main/dcdflib.cpp"
    "${GCTA_ROOT}/src/Logger.cpp")
target_include_directories(gcta_test_support PUBLIC
    "${GCTA_ROOT}/main"
    "${GCTA_ROOT}/include")
target_link_libraries(gcta_test_support PUBLIC Eigen3::Eigen)

add_executable(cojo_unit_tests
    unit/test_commfunc.cpp
    unit/test_strfunc.cpp
    unit/test_cojo_algebra.cpp
    unit/test_statfunc.cpp)
target_link_libraries(cojo_unit_tests PRIVATE gcta_test_support GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(cojo_unit_tests
    DISCOVERY_TIMEOUT 30
    PROPERTIES LABELS "unit")

# Tier A golden CLI tests — skip cleanly until fixtures exist (--allow-empty).
set(GOLDEN_SCRIPT "${GCTA_ROOT}/scripts/ci/golden_test.sh")
set(GOLDEN_TESTS_DIR "${GCTA_ROOT}/tests")
if(TARGET gcta64)
    set(_gcta_for_golden "$<TARGET_FILE:gcta64>")
else()
    set(_gcta_for_golden "gcta64")
endif()
add_test(
    NAME cojo_golden
    COMMAND bash "${GOLDEN_SCRIPT}"
        --gcta "${_gcta_for_golden}"
        --tests-dir "${GOLDEN_TESTS_DIR}")
set_tests_properties(cojo_golden PROPERTIES
    LABELS "golden"
    WORKING_DIRECTORY "${GCTA_ROOT}")
