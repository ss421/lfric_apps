!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Standalone adjoint of invoke_convert_cart2sphere_vector subroutine.
module invoke_adj_convert_cart2sphere_vector_mod

  use constants_mod,       only: r_def, PI
  use field_mod,           only: field_type, field_proxy_type
  use coord_transform_mod, only: sphere2cart_vector

  implicit none

  contains

  !> @brief Applies the adjoint of invoke_convert_cart2sphere_vector
  !> @param[in,out] field  A bundle of 3 fields to be transformed
  !> @param[in]     coords A bundle of 3 fields that includes the
  !>                       field coordinates
  subroutine invoke_adj_convert_cart2sphere_vector(field, coords)

    implicit none

    type(field_type), intent(inout) :: field(3)
    type(field_type), intent(in)    :: coords(3)

    ! Local
    type(field_proxy_type) :: f_p(3), x_p(3)
    integer                :: i, df, undf

    do i = 1,3
      f_p(i) = field(i)%get_proxy()
      x_p(i) = coords(i)%get_proxy()
    end do

    undf = f_p(1)%vspace%get_last_dof_annexed()

!Please see PSyclone issues #1351 regarding this implementation
!$omp parallel default(none)                                                   &
!$omp private(df)                                                              &
!$omp shared(undf,f_p,x_p)
!$omp do schedule(static)
    do df = undf, 1, -1
      call adj_cart2sphere_scalar( x_p(1)%data(df), x_p(2)%data(df), &
                                   x_p(3)%data(df), f_p(1)%data(df), &
                                   f_p(2)%data(df), f_p(3)%data(df) )
    end do
!$omp end do
!$omp end parallel

    call f_p(1)%set_dirty()
    call f_p(2)%set_dirty()
    call f_p(3)%set_dirty()

  end subroutine invoke_adj_convert_cart2sphere_vector

  !> @brief Applies the adjoint of cart2sphere_scalar
  !> @param[in]    x x component of location in Cartesian coodinates
  !> @param[in]    y y component of location in Cartesian coodinates
  !> @param[in]    z z component of location in Cartesian coodinates
  !> @param[inout] u u component of flux vector in spherical
  !>                 coordinates, out in Cartesian coordinates
  !> @param[inout] v v component of flux vector in spherical
  !>                 coordinates, out in Cartesian coordinates
  !> @param[inout] w w component of flux vector in spherical
  !>                 coordinates, out in Cartesian coordinates
  subroutine adj_cart2sphere_scalar(x, y, z, u, v, w)

    implicit none

    real(kind=r_def), intent(in)    :: x, y, z
    real(kind=r_def), intent(inout) :: u, v, w

    ! Local
    real(kind=r_def) :: t, r, phi, c1, c2
    real(kind=r_def) :: u_initial, v_initial, w_initial

    t = x**2 + y**2
    r = sqrt(t + z**2)
    phi = 0.5_r_def * PI - acos(z / r)

    ! Note: adjoint formed by taking original subroutine and adjointing
    ! equations by hand, rather than writing line-by-line adjoint
    c1 = r * cos(phi) / t
    c2 = z / (r * sqrt(t))

    u_initial = u
    v_initial = v
    w_initial = w

    u = u_initial * (-y * c1) + v_initial * (-x * c2)     + w_initial * (x / r)
    v = u_initial * (x * c1)  + v_initial * (-y * c2)     + w_initial * (y / r)
    w =                         v_initial * (sqrt(t) / r) + w_initial * (z / r)

  end subroutine adj_cart2sphere_scalar

end module invoke_adj_convert_cart2sphere_vector_mod
