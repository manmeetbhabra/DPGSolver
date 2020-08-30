#find_package(MPI REQUIRED)
#include_directories(${MPI_C_INCLUDE_PATH})

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

find_package(GSL 2.4 REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

set(RUN_PETSC_TESTS FALSE)
find_package(PETSc 3.8 REQUIRED)
include_directories(${PETSC_INCLUDES})

set(RUN_SLEPC_TESTS FALSE)
find_package(SLEPc 3.8 REQUIRED)
include_directories(${SLEPC_INCLUDE_DIRS})

find_package(NLPQLP REQUIRED)

# Currently including all directories here so that `#include`s need not specify relative paths.
include_directories(
	allocators
	containers
	definitions
	element
	general
	intrusive
	simulation
	templates
	simulation/computational_elements
	simulation/geometry
	simulation/mesh
	simulation/operators
	simulation/optimization
	simulation/solvers
	simulation/solvers/adaptation
	simulation/solvers/dg
	simulation/solvers/dpg
	simulation/solvers/opg
	simulation/solvers/numerical_flux
	simulation/test_case
	simulation/test_case/boundary
	simulation/test_case/compute_error
	simulation/test_case/flux
	simulation/test_case/solution
	simulation/test_case/solution/advection
	simulation/test_case/solution/diffusion
	simulation/test_case/solution/euler
	simulation/test_case/solution/navier_stokes
	simulation/test_case/solution/burgers_inviscid
	simulation/test_case/solution/restart
	simulation/visualization
	testing
	testing/integration
	)

configure_file_h_in(${PROJECT_SOURCE_DIR}/src/ definitions_core.h.in)

add_include_directories()

add_subdirectory(containers)
add_subdirectory(element)
add_subdirectory(general)
add_subdirectory(intrusive)
add_subdirectory(simulation)
add_subdirectory(testing)

set	(SOURCE
	 core.c
	)

set	(LIBS_DEPEND
	)

set(LIB Core)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
