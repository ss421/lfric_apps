!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the adjoint miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module adjoint_test_driver_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use field_mod,                  only : field_type
  use geometric_constants_mod,    only : get_coordinates, &
                                         get_panel_id
  use log_mod,                    only : log_event, LOG_LEVEL_INFO
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection

  implicit none

  private

  public run

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief   Runs adjoint tests.
  !> @details Runs algorithm layer adjoint tests.
  subroutine run()

    ! PSyAD generated tests
    use gen_adj_kernel_tests_mod,                   only : run_gen_adj_kernel_tests

    ! Handwritten tests
    use adjt_interpolation_alg_mod,                 only : adjt_interp_w3wth_to_w2_alg, &
                                                           adjt_interp_w2_to_w3wth_alg
    use adjt_convert_cart2sphere_vector_alg_mod,    only : adjt_convert_cart2sphere_vector_alg
    use adjt_sci_convert_hdiv_field_alg_mod,        only : adjt_sci_convert_hdiv_field_alg

    implicit none

    type(field_type),             pointer :: chi(:) => null()
    type(field_type),             pointer :: panel_id => null()
    type(mesh_type),              pointer :: mesh => null()

    mesh => mesh_collection%get_mesh( prime_mesh_name )
    chi => get_coordinates( mesh%get_id() )
    panel_id => get_panel_id( mesh%get_id() )

    call log_event( "TESTING generated adjoint kernels", LOG_LEVEL_INFO )
    call run_gen_adj_kernel_tests( mesh, chi, panel_id )
    call log_event( "TESTING handwritten adjoints", LOG_LEVEL_INFO )
    call adjt_convert_cart2sphere_vector_alg( mesh )
    call adjt_sci_convert_hdiv_field_alg( mesh, chi, panel_id )
    call adjt_interp_w3wth_to_w2_alg( mesh )
    call adjt_interp_w2_to_w3wth_alg( mesh )
    call log_event( "TESTING COMPLETE", LOG_LEVEL_INFO )

  end subroutine run

end module adjoint_test_driver_mod
