# =============================================================================
#  Bounded Distortion Parametrization in the Space of Metrics
#  (Chien, Levi, Weber - SIGGRAPH Asia 2016)
#
#  C++ (CGAL + Eigen + Polyscope) front end, MATLAB (Engine API) numerical back end.
#
#  All C++ dependencies are downloaded into <source>/external on the first
#  configure (FetchContent):
#     Eigen 3.4.0          header only
#     Boost 1.86 headers   header only (required by CGAL)
#     CGAL 6.0.1           header only
#     Polyscope 2.3.0      built from source (brings GLFW, GLAD, ImGui, GLM)
#
#  MATLAB: the installation is taken from the environment variable MATLAB_64_DIR
#  (the MATLAB root folder, e.g. C:\Program Files\MATLAB\R2024b; a subfolder
#  such as <root>\bin\win64 is accepted too):
#     set MATLAB_64_DIR=C:\Program Files\MATLAB\R2024b
#     cmake -S . -B build
# =============================================================================
cmake_minimum_required(VERSION 3.20)
project(BoundedDistortionParametrization LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(BDP_USE_MATLAB    "Build the MATLAB Engine API bridge"          ON)
option(BDP_USE_POLYSCOPE "Build the Polyscope viewer"                  ON)

if(POLICY CMP0135)            # timestamps of extracted archives = extraction time
  cmake_policy(SET CMP0135 NEW)
endif()

# -----------------------------------------------------------------------------
#  Dependencies, downloaded into the local working directory
# -----------------------------------------------------------------------------
set(BDP_DEPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external" CACHE PATH
    "Download/build directory of the third-party dependencies")
include(FetchContent)
# normal variable: overrides FetchContent's cached default (<build>/_deps)
set(FETCHCONTENT_BASE_DIR "${BDP_DEPS_DIR}")
set(FETCHCONTENT_QUIET OFF)

# Header-only packages: SOURCE_SUBDIR points to a folder without CMakeLists.txt,
# so FetchContent_MakeAvailable only downloads/extracts them (no add_subdirectory,
# no tests/docs of the upstream projects are configured).
FetchContent_Declare(eigen
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
  SOURCE_SUBDIR _header_only_)
FetchContent_Declare(boost_headers
  URL https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-b2-nodocs.tar.xz
  SOURCE_SUBDIR _header_only_)
FetchContent_Declare(cgal
  URL https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1-library.tar.xz
  SOURCE_SUBDIR _header_only_)
FetchContent_MakeAvailable(eigen boost_headers cgal)

add_library(bdp_eigen INTERFACE)
target_include_directories(bdp_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")

add_library(bdp_cgal INTERFACE)
target_include_directories(bdp_cgal SYSTEM INTERFACE "${cgal_SOURCE_DIR}/include" "${boost_headers_SOURCE_DIR}")
target_compile_definitions(bdp_cgal INTERFACE
  CGAL_DISABLE_GMP=1          # only the double-precision Simple_cartesian kernel is used
  CGAL_NO_GMP=1
  BOOST_ALL_NO_LIB=1)         # header-only Boost, no auto-linking on MSVC

if(BDP_USE_POLYSCOPE)
  FetchContent_Declare(polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        v2.3.0
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(polyscope)
endif()

# -----------------------------------------------------------------------------
#  MATLAB Engine API for C++
# -----------------------------------------------------------------------------
if(BDP_USE_MATLAB)
  # The MATLAB installation is located through the environment variable
  # MATLAB_64_DIR. It may point to the MATLAB root or to a folder inside it
  # (e.g. <root>/bin/win64); we walk up until <root>/extern/include/MatlabEngine.hpp
  # is found and hand that root to CMake's FindMatlab as Matlab_ROOT_DIR.
  if(NOT DEFINED ENV{MATLAB_64_DIR} OR "$ENV{MATLAB_64_DIR}" STREQUAL "")
    message(FATAL_ERROR
      "BDP_USE_MATLAB is ON but the environment variable MATLAB_64_DIR is not set.\n"
      "Set it to your MATLAB installation, e.g.\n"
      "  Windows : set MATLAB_64_DIR=C:\\Program Files\\MATLAB\\R2024b\n"
      "  Linux   : export MATLAB_64_DIR=/usr/local/MATLAB/R2024b\n"
      "then configure again in that shell (or use -DBDP_USE_MATLAB=OFF).")
  endif()
  file(TO_CMAKE_PATH "$ENV{MATLAB_64_DIR}" BDP_MATLAB_ENV_DIR)
  string(REGEX REPLACE "/+$" "" BDP_MATLAB_ENV_DIR "${BDP_MATLAB_ENV_DIR}")
  set(BDP_MATLAB_ROOT "")
  set(_bdp_candidate "${BDP_MATLAB_ENV_DIR}")
  foreach(_bdp_level RANGE 0 3)
    if(EXISTS "${_bdp_candidate}/extern/include/MatlabEngine.hpp")
      set(BDP_MATLAB_ROOT "${_bdp_candidate}")
      break()
    endif()
    get_filename_component(_bdp_candidate "${_bdp_candidate}" DIRECTORY)
  endforeach()
  if(NOT BDP_MATLAB_ROOT)
    message(FATAL_ERROR
      "MATLAB_64_DIR = '$ENV{MATLAB_64_DIR}' does not contain a MATLAB installation with the "
      "C++ Engine API (extern/include/MatlabEngine.hpp, MATLAB R2018b or newer).")
  endif()
  # FORCE: the environment variable is the single source of truth for the location
  set(Matlab_ROOT_DIR "${BDP_MATLAB_ROOT}" CACHE PATH "MATLAB root (from MATLAB_64_DIR)" FORCE)
  message(STATUS "MATLAB_64_DIR    : $ENV{MATLAB_64_DIR}")

  # Provides Matlab_INCLUDE_DIRS, Matlab_ENGINE_LIBRARY (libMatlabEngine) and
  # Matlab_DATAARRAY_LIBRARY (libMatlabDataArray).
  find_package(Matlab REQUIRED COMPONENTS ENGINE_LIBRARY DATAARRAY_LIBRARY)
  message(STATUS "MATLAB root      : ${Matlab_ROOT_DIR}")
  message(STATUS "MATLAB engine lib: ${Matlab_ENGINE_LIBRARY}")
endif()

# -----------------------------------------------------------------------------
#  Executable
# -----------------------------------------------------------------------------
set(BDP_SOURCES
  src/main.cpp
  src/Mesh.cpp
  src/MetricData.cpp
  src/Geometry.cpp
  src/Cones.cpp
  src/Layout.cpp
  src/MatlabBridge.cpp)
if(BDP_USE_POLYSCOPE)
  list(APPEND BDP_SOURCES src/Viewer.cpp)
endif()

add_executable(bdp ${BDP_SOURCES})
target_include_directories(bdp PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(bdp PRIVATE bdp_eigen bdp_cgal)
# Location of bdp_solve.m and friends (added to the MATLAB path at start-up)
target_compile_definitions(bdp PRIVATE BDP_MATLAB_DIR="${CMAKE_CURRENT_SOURCE_DIR}/matlab")

if(BDP_USE_POLYSCOPE)
  target_link_libraries(bdp PRIVATE polyscope)
  target_compile_definitions(bdp PRIVATE BDP_HAVE_POLYSCOPE)
endif()

if(BDP_USE_MATLAB)
  target_include_directories(bdp PRIVATE ${Matlab_INCLUDE_DIRS})
  target_link_libraries(bdp PRIVATE ${Matlab_ENGINE_LIBRARY} ${Matlab_DATAARRAY_LIBRARY})
  target_compile_definitions(bdp PRIVATE BDP_HAVE_MATLAB)
  find_package(Threads REQUIRED)
  target_link_libraries(bdp PRIVATE Threads::Threads)
endif()

if(MSVC)
  target_compile_options(bdp PRIVATE /bigobj /permissive- /Zc:__cplusplus /EHsc)
  target_compile_definitions(bdp PRIVATE NOMINMAX _USE_MATH_DEFINES)
endif()

# -----------------------------------------------------------------------------
#  Run helpers: the MATLAB Engine DLLs/shared objects must be on the search path
#  at run time (extern/bin/<arch> and bin/<arch> of the MATLAB installation).
# -----------------------------------------------------------------------------
if(BDP_USE_MATLAB)
  if(WIN32)
    set(BDP_MATLAB_RUNTIME "${Matlab_ROOT_DIR}/extern/bin/win64;${Matlab_ROOT_DIR}/bin/win64")
    file(TO_NATIVE_PATH "${BDP_MATLAB_RUNTIME}" BDP_MATLAB_RUNTIME_NATIVE)
    # run_bdp.bat next to bdp.exe (one per configuration)
    file(GENERATE OUTPUT "$<TARGET_FILE_DIR:bdp>/run_bdp.bat"
         CONTENT "@echo off\r\nset \"PATH=${BDP_MATLAB_RUNTIME_NATIVE};%PATH%\"\r\n\"%~dp0bdp.exe\" %*\r\n")
    # used by "Start Debugging" in Visual Studio
    set_target_properties(bdp PROPERTIES VS_DEBUGGER_ENVIRONMENT "PATH=${BDP_MATLAB_RUNTIME_NATIVE};%PATH%")
  else()
    if(APPLE)
      set(BDP_ARCH maci64)
      if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        set(BDP_ARCH maca64)
      endif()
      set(BDP_LDVAR DYLD_LIBRARY_PATH)
    else()
      set(BDP_ARCH glnxa64)
      set(BDP_LDVAR LD_LIBRARY_PATH)
    endif()
    file(GENERATE OUTPUT "$<TARGET_FILE_DIR:bdp>/run_bdp.sh"
         CONTENT "#!/bin/sh\nexport ${BDP_LDVAR}=\"${Matlab_ROOT_DIR}/extern/bin/${BDP_ARCH}:${Matlab_ROOT_DIR}/bin/${BDP_ARCH}:\$${BDP_LDVAR}\"\nexec \"$(dirname \"$0\")/bdp\" \"$@\"\n"
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()
endif()
