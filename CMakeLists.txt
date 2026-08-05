cmake_minimum_required(VERSION 3.15)
project(PySCF_CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Find BLAS/LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find libcint
set(CINT_LIBRARY "/usr/local/lib/libcint.so" CACHE FILEPATH "libcint library")
set(CINT_INCLUDE_DIR "/usr/local/include" CACHE PATH "libcint include directory")

if(NOT EXISTS ${CINT_LIBRARY})
    message(WARNING "libcint not found at ${CINT_LIBRARY}, integral calculations will be approximate")
    set(CINT_LIBRARY "")
endif()

# Create the library
add_library(pyscf_core STATIC
    src/molecule.cpp
    src/gto/gto_integrals.cpp
    src/dft/lebedev_grid.cpp
    src/dft/b3lyp_xc.cpp
    src/dft/scf.cpp
    src/geom/optimizer.cpp
    src/bond/mayer_bond.cpp
    src/utils/basis_parser.cpp
)

target_include_directories(pyscf_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(CINT_INCLUDE_DIR)
    target_include_directories(pyscf_core PRIVATE ${CINT_INCLUDE_DIR})
endif()

target_link_libraries(pyscf_core PUBLIC
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    m
    pthread
)

if(CINT_LIBRARY)
    target_link_libraries(pyscf_core PUBLIC ${CINT_LIBRARY})
endif()

# Create the main executable
add_executable(pyscf_demo examples/demo.cpp)
target_link_libraries(pyscf_demo PRIVATE pyscf_core)

# Create H2O demo
add_executable(pyscf_h2o_demo examples/demo_h2o.cpp)
target_link_libraries(pyscf_h2o_demo PRIVATE pyscf_core)

# Create geometry optimization demo
add_executable(pyscf_optimize_demo examples/demo_optimize.cpp)
target_link_libraries(pyscf_optimize_demo PRIVATE pyscf_core)

# Create bond analysis demo
add_executable(pyscf_bond_demo examples/demo_bond.cpp)
target_link_libraries(pyscf_bond_demo PRIVATE pyscf_core)

# Installation
install(TARGETS pyscf_core pyscf_demo pyscf_h2o_demo pyscf_optimize_demo pyscf_bond_demo
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY include/
    DESTINATION include
)

install(DIRECTORY basis/
    DESTINATION share/pyscf_cpp/basis
)
