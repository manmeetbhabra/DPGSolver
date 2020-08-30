add_subdirectory(mesh)
add_subdirectory(operators)

set	(SOURCE
	 simulation.c

	 computational_elements/computational_elements.c
	 computational_elements/face.c
	 computational_elements/volume.c

	 geometry/element_geometry.c
	 geometry/geometry.c
	 geometry/geometry_blended.c
	 geometry/geometry_parametric.c
	 geometry/geometry_surface.c
	 geometry/geometry_normals.c

	 optimization/optimize.c
	 optimization/optimization_case.c
	 optimization/optimizer_line_search.c
	 optimization/optimizer_NLPQLP.c
	 optimization/functionals.c
	 optimization/adjoint.c
	 optimization/sensitivities.c
	 optimization/gradient.c
	 optimization/output_progress.c

	 solvers/compute_rlhs.c
	 solvers/compute_face_rlhs.c
	 solvers/compute_volume_rlhs.c
	 solvers/element_solver.c
	 solvers/face_solver.c
	 solvers/solve.c
	 solvers/solve_explicit.c
	 solvers/solve_implicit.c
	 solvers/volume_solver.c

	 solvers/adaptation/adaptation.c
	 solvers/adaptation/face_solver_adaptive.c
	 solvers/adaptation/element_adaptation.c
	 solvers/adaptation/volume_solver_adaptive.c

	 solvers/dg/compute_face_rlhs_dg.c
	 solvers/dg/compute_grad_coef_dg.c
	 solvers/dg/compute_source_rlhs_dg.c
	 solvers/dg/compute_volume_rlhs_dg.c
	 solvers/dg/element_solver_dg.c
	 solvers/dg/face_solver_dg.c
	 solvers/dg/solve_dg.c
	 solvers/dg/volume_solver_dg.c

	 solvers/dpg/compute_all_rlhs_dpg.c
	 solvers/dpg/face_solver_dpg.c
	 solvers/dpg/element_solver_dpg.c
	 solvers/dpg/solve_dpg.c
	 solvers/dpg/volume_solver_dpg.c

	 solvers/opg/compute_volume_rlhs_opg.c
	 solvers/opg/compute_face_rlhs_opg.c
	 solvers/opg/face_solver_opg.c
	 solvers/opg/element_solver_opg.c
	 solvers/opg/penalty_opg.c
	 solvers/opg/penalty_opg_advection.c
	 solvers/opg/solve_opg.c
	 solvers/opg/volume_solver_opg.c

	 solvers/numerical_flux/numerical_flux.c
	 solvers/numerical_flux/numerical_flux_advection.c
	 solvers/numerical_flux/numerical_flux_diffusion.c
	 solvers/numerical_flux/numerical_flux_euler.c
	 solvers/numerical_flux/numerical_flux_navier_stokes.c
	 solvers/numerical_flux/numerical_flux_burgers_inviscid.c

	 test_case/test_case.c
	 test_case/test_case_c.c

	 test_case/boundary/boundary.c
	 test_case/boundary/boundary_advection.c
	 test_case/boundary/boundary_diffusion.c
	 test_case/boundary/boundary_euler.c
	 test_case/boundary/boundary_navier_stokes.c

	 test_case/compute_error/compute_error.c
	 test_case/compute_error/compute_error_advection.c
	 test_case/compute_error/compute_error_diffusion.c
	 test_case/compute_error/compute_error_euler.c
	 test_case/compute_error/compute_error_navier_stokes.c
	 test_case/compute_error/compute_error_burgers_inviscid.c

	 test_case/flux/flux.c
	 test_case/flux/flux_advection.c
	 test_case/flux/flux_diffusion.c
	 test_case/flux/flux_euler.c
	 test_case/flux/flux_navier_stokes.c
	 test_case/flux/flux_burgers_inviscid.c

	 test_case/solution/element_solution.c
	 test_case/solution/solution.c
	 test_case/solution/advection/solution_advection.c
	 test_case/solution/advection/solution_advection_default.c
	 test_case/solution/advection/free_stream_advection/solution_free_stream_advection.c
	 test_case/solution/advection/peterson/solution_peterson.c
	 test_case/solution/advection/vortex_advection/solution_vortex_advection.c
	 test_case/solution/advection/hyperbolic_tan/solution_hyperbolic_tan.c
	 test_case/solution/diffusion/solution_diffusion.c
	 test_case/solution/diffusion/default_steady/solution_diffusion_default_steady.c
	 test_case/solution/euler/solution_euler.c
	 test_case/solution/euler/periodic_vortex/solution_periodic_vortex.c
	 test_case/solution/euler/supersonic_vortex/solution_supersonic_vortex.c
	 test_case/solution/euler/free_stream/solution_free_stream.c
	 test_case/solution/navier_stokes/solution_navier_stokes.c
	 test_case/solution/navier_stokes/taylor_couette/solution_taylor_couette.c
	 test_case/solution/burgers_inviscid/solution_burgers_inviscid.c
	 test_case/solution/burgers_inviscid/trigonometric/solution_trigonometric.c
	 test_case/solution/restart/approximate_nearest_neighbor.c
	 test_case/solution/restart/inverse_mapping.c
	 test_case/solution/restart/restart.c
	 test_case/solution/restart/restart_writers.c

	 visualization/element_plotting.c
	 visualization/visualization.c
	)

set	(LIBS_DEPEND
	 ${MPI_C_LIBRARIES}
	 ${PETSC_LIBRARIES}
	 ${NLPQLP_LIBRARIES}
	 Containers
	 Intrusive
	 Mesh
	)

set(LIB Simulation)
add_library(${LIB} STATIC ${SOURCE})
target_link_libraries(${LIB} ${LIBS_DEPEND})
