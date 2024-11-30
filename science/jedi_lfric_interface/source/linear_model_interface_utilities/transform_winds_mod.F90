!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief  Wind transform code used to interpolate between vector and scalar
!>        winds
!>
!>@detail This code will be replaced by code in #3734 when it is ready.
!>
module transform_winds_mod

  use constants_mod,                 only: i_def
  use field_collection_mod,          only: field_collection_type
  use field_mod,                     only: field_type, field_proxy_type

  implicit none

  private
  public :: wind_scalar_to_vector
  public :: adj_wind_scalar_to_vector
  public :: wind_vector_to_scalar
  public :: adj_wind_vector_to_scalar

  contains

  !> @brief   Interpolate the cell-centre scalar winds to the edge based vector wind
  !>
  !> @param[in,out] fields  A field collection that contains the wind fields
  subroutine wind_scalar_to_vector( fields )

    use interpolation_alg_mod, only : interp_w3wth_to_w2_alg

    implicit none

    type( field_collection_type ), intent(inout) :: fields

    ! Local
    type( field_type ), pointer :: u_in_w3
    type( field_type ), pointer :: v_in_w3
    type( field_type ), pointer :: w_in_wth
    type( field_type ), pointer :: u_in_w2

    nullify( u_in_w3, v_in_w3, w_in_wth, u_in_w2 )

    ! Get the fields
    call fields%get_field('u_in_w3', u_in_w3)
    call fields%get_field('v_in_w3', v_in_w3)
    call fields%get_field('w_in_wth', w_in_wth)
    call fields%get_field('u', u_in_w2)

    ! Interpolate cell centred zonal/meridional winds
    call interp_w3wth_to_w2_alg(u_in_w2, u_in_w3, v_in_w3, w_in_wth)

  end subroutine wind_scalar_to_vector

  !> @brief  Adjoint of interpolation of the cell-centred scalar winds to the edge-based vector wind
  !>
  !> @param[in,out] fields  A field collection that contains the wind fields
  subroutine adj_wind_scalar_to_vector( fields )

    use adj_interpolation_alg_mod, only : adj_interp_w3wth_to_w2_alg

    implicit none

    type( field_collection_type ), intent(inout) :: fields

    ! Local
    type( field_type ), pointer :: u_in_w3
    type( field_type ), pointer :: v_in_w3
    type( field_type ), pointer :: w_in_wth
    type( field_type ), pointer :: u_in_w2

    nullify( u_in_w3, v_in_w3, w_in_wth, u_in_w2 )

    ! Get the fields
    call fields%get_field('u_in_w3', u_in_w3)
    call fields%get_field('v_in_w3', v_in_w3)
    call fields%get_field('w_in_wth', w_in_wth)
    call fields%get_field('u', u_in_w2)

    ! Adjoint-interpolate cell centred zonal/meridional winds
    call adj_interp_w3wth_to_w2_alg(u_in_w2, u_in_w3, v_in_w3, w_in_wth)

  end subroutine adj_wind_scalar_to_vector

  !> @brief  Interpolate the edge based vector wind to the cell centre scalar winds
  !>
  !> @param[in,out] fields  A field collection that contains the wind fields
  subroutine wind_vector_to_scalar( fields )

    use interpolation_alg_mod, only : interp_w2_to_w3wth_alg

    implicit none

    type( field_collection_type ), intent(inout) :: fields

    ! Local
    type( field_type ),          pointer :: u_in_w3
    type( field_type ),          pointer :: v_in_w3
    type( field_type ),          pointer :: w_in_wth
    type( field_type ),          pointer :: u_in_w2

    nullify( u_in_w3, v_in_w3, w_in_wth, u_in_w2 )

    ! Get the fields from the prognostic_fields collection
    call fields%get_field('u_in_w3', u_in_w3)
    call fields%get_field('v_in_w3', v_in_w3)
    call fields%get_field('w_in_wth', w_in_wth)
    call fields%get_field('u', u_in_w2)

    call interp_w2_to_w3wth_alg(u_in_w2, u_in_w3, v_in_w3, w_in_wth)

  end subroutine wind_vector_to_scalar

  !> @brief  Adjoint of interpolation of the edge-based vector wind to the cell-centred scalar winds
  !>
  !> @param[in,out] fields  A field collection that contains the wind fields
  subroutine adj_wind_vector_to_scalar( fields )

    use adj_interpolation_alg_mod, only : adj_interp_w2_to_w3wth_alg

    implicit none

    type( field_collection_type ), intent(inout) :: fields

    ! Local
    type( field_type ),          pointer :: u_in_w3
    type( field_type ),          pointer :: v_in_w3
    type( field_type ),          pointer :: w_in_wth
    type( field_type ),          pointer :: u_in_w2

    nullify( u_in_w3, v_in_w3, w_in_wth, u_in_w2 )

    ! Get the fields from the prognostic_fields collection
    call fields%get_field('u_in_w3', u_in_w3)
    call fields%get_field('v_in_w3', v_in_w3)
    call fields%get_field('w_in_wth', w_in_wth)
    call fields%get_field('u', u_in_w2)

    call adj_interp_w2_to_w3wth_alg(u_in_w2, u_in_w3, v_in_w3, w_in_wth)

  end subroutine adj_wind_vector_to_scalar

end module transform_winds_mod
