!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the horizontal mass flux for FFSL using a first-order scheme
!!        as the update part of the cheap update.
!> @details This kernel computes the flux in the x and y directions using a
!!          constant reconstruction. This is multiplied by the fractional
!!          velocity to give the flux. As this kernel is only used for the
!!          update part of the cheap update it is assumed that |CFL| < 1, and
!!          thus this kernel must only be called with a stencil extent of 1.
!!          This kernel is used for the FFSL horizontal transport scheme.
!!
!> @note This kernel only works when field is a W3 field at lowest order since
!!       it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing the
!!       relevant dofmaps.

module ffsl_flux_xy_first_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    CELL_COLUMN, GH_WRITE,     &
                                    GH_READ, GH_SCALAR,        &
                                    STENCIL, X1D, GH_INTEGER,  &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    Y1D
  use constants_mod,         only : i_def, r_tran
  use fs_continuity_mod,     only : W3
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : W, E, N, S

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_flux_xy_first_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),          & ! field_for_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(Y1D)),          & ! field_for_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! frac_dry_flux
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ew
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ns
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     )                         & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_flux_xy_first_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: ffsl_flux_xy_first_code

contains

  !> @brief Compute the horizontal fluxes for FFSL using a first order scheme.
  !> @param[in]     nlayers             Number of layers
  !> @param[in,out] flux                The output first order flux
  !> @param[in]     field_for_x         Field to use in evaluating x-flux
  !> @param[in]     stencil_size_x      Local length of field_for_x W3 stencil
  !> @param[in]     stencil_map_x       Dofmap for the field_for_x stencil
  !> @param[in]     field_for_y         Field to use in evaluating x-flux
  !> @param[in]     stencil_size_y      Local length of field_for_y W3 stencil
  !> @param[in]     stencil_map_y       Dofmap for the field_for_y stencil
  !> @param[in]     frac_dry_flux       Fractional part of the dry flux or wind
  !> @param[in]     face_selector_ew    2D field indicating which W/E faces to
  !!                                    loop over for this column
  !> @param[in]     face_selector_ns    2D field indicating which N/S faces to
  !!                                    loop over for this column
  !> @param[in]     dt                  Time step
  !> @param[in]     ndf_w2h             Num of DoFs for W2h per cell
  !> @param[in]     undf_w2h            Num of DoFs for W2h in this partition
  !> @param[in]     map_w2h             Map for W2h
  !> @param[in]     ndf_w3              Num of DoFs for W3 per cell
  !> @param[in]     undf_w3             Num of DoFs for W3 in this partition
  !> @param[in]     map_w3              Map for W3
  !> @param[in]     ndf_w3_2d           Num of DoFs for 2D W3 per cell
  !> @param[in]     undf_w3_2d          Num of DoFs for this partition for 2D W3
  !> @param[in]     map_w3_2d           Map for 2D W3
  subroutine ffsl_flux_xy_first_code( nlayers,             &
                                      flux,                &
                                      field_for_x,         &
                                      stencil_size_x,      &
                                      stencil_map_x,       &
                                      field_for_y,         &
                                      stencil_size_y,      &
                                      stencil_map_y,       &
                                      frac_dry_flux,       &
                                      face_selector_ew,    &
                                      face_selector_ns,    &
                                      dt,                  &
                                      ndf_w2h,             &
                                      undf_w2h,            &
                                      map_w2h,             &
                                      ndf_w3,              &
                                      undf_w3,             &
                                      map_w3,              &
                                      ndf_w3_2d,           &
                                      undf_w3_2d,          &
                                      map_w3_2d )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_w3_2d
    integer(kind=i_def), intent(in) :: stencil_size_x
    integer(kind=i_def), intent(in) :: stencil_size_y

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_w3,stencil_size_x)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_w3,stencil_size_y)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: field_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: frac_dry_flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dt
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)

    ! Integers
    integer(kind=i_def) :: local_dofs_x(2), local_dofs_y(2)
    integer(kind=i_def) :: dof_iterator
    integer(kind=i_def) :: stencil_idx
    integer(kind=i_def) :: sign_displacement
    integer(kind=i_def) :: rel_idx
    integer(kind=i_def) :: dof_offset, sign_offset
    integer(kind=i_def) :: k, idx_3d
    integer(kind=i_def) :: lam_edge_size

    ! x-direction
    local_dofs_x = (/ W, E /)
    local_dofs_y = (/ S, N /)

    ! The stencil should be of size 3 (i.e. extent = 1) therefore we
    ! are at the edge of a LAM if the stencil size is less than 3
    lam_edge_size = 3_i_def

    ! = X Calculation ======================================================== !

    if ( lam_edge_size > stencil_size_x ) then

      ! At edge of LAM, so set output to zero ----------------------------------
      do k = 0, nlayers - 1
        do dof_iterator = 1, face_selector_ew(map_w3_2d(1))
          flux( map_w2h(local_dofs_x(dof_iterator)) + k ) = 0.0_r_tran
        end do
      end do

    else

      ! Not at edge of LAM so compute fluxes -----------------------------------

      ! Loop over the x direction dofs to compute flux at each dof
      do dof_iterator = 1, face_selector_ew(map_w3_2d(1))

        ! Pull out index to avoid multiple indirections
        idx_3d = map_w2h(local_dofs_x(dof_iterator))

        ! Set a local offset, dependent on the face we are looping over
        select case (local_dofs_x(dof_iterator))
        case ( W, S )
          dof_offset = 0
        case ( E, N )
          dof_offset = 1
        end select

        ! Loop over vertical levels
        do k = 0, nlayers - 1

          ! Flux calculation for a single face =================================

          ! Pull out dry flux / dep point, the sign determines the upwind cell
          ! as we assume |CFL| < 1
          sign_displacement = int(sign(1.0_r_tran, frac_dry_flux(idx_3d + k)))

          ! Set an offset for the stencil index, based on dep point sign
          sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

          ! The local index of the departure cell
          rel_idx = dof_offset + sign_offset - 1

          ! Fractional part ====================================================
          ! Extract reconstruction data

          ! Determine the index in the stencil from rel_idx
          ! e.g. for extent 1:
          ! Relative idx is   | -1 |  0 |  1 |
          ! Stencil has order |  2 |  1 |  3 |
          stencil_idx = 1 + ABS(rel_idx) + (1 - SIGN(1, -rel_idx))/2

          ! Assign flux ========================================================
          flux(idx_3d + k) = (field_for_x(stencil_map_x(1,stencil_idx) + k) &
                              * frac_dry_flux(idx_3d + k) ) / dt

        end do ! vertical levels k

      end do ! dof_iterator
    end if

    ! = Y Calculation ======================================================== !

    if ( lam_edge_size > stencil_size_y ) then

      ! At edge of LAM, so set output to zero ----------------------------------
      do k = 0, nlayers - 1
        do dof_iterator = 1, face_selector_ns(map_w3_2d(1))
          flux( map_w2h(local_dofs_y(dof_iterator)) + k ) = 0.0_r_tran
        end do
      end do

    else

      ! Not at edge of LAM so compute fluxes -----------------------------------

      ! Loop over the y direction dofs to compute flux at each dof
      do dof_iterator = 1, face_selector_ns(map_w3_2d(1))

        ! Pull out index to avoid multiple indirections
        idx_3d = map_w2h(local_dofs_y(dof_iterator))

        ! Set a local offset, dependent on the face we are looping over
        select case (local_dofs_y(dof_iterator))
        case ( W, S )
          dof_offset = 0
        case ( E, N )
          dof_offset = 1
        end select

        ! Loop over vertical levels
        do k = 0, nlayers - 1

          ! Flux calculation for a single face =================================

          ! Pull out dry flux / dep point, the sign determines the upwind cell
          ! as we assume |CFL| < 1
          ! NB: minus signs are needed because the Y1D stencil runs from S to N
          ! but winds blowing S to N have negative sign
          sign_displacement = int(sign(1.0_r_tran, -frac_dry_flux(idx_3d + k)))

          ! Set an offset for the stencil index, based on dep point sign
          sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

          ! The local index of the departure cell
          rel_idx = dof_offset + sign_offset - 1

          ! Fractional part ====================================================
          ! Extract reconstruction data

          ! Determine the index in the stencil from rel_idx
          ! e.g. for extent 1:
          ! Relative idx is   | -1 |  0 |  1 |
          ! Stencil has order |  2 |  1 |  3 |
          stencil_idx = 1 + ABS(rel_idx) + (1 - SIGN(1, -rel_idx))/2

          ! Assign flux ========================================================
          flux(idx_3d + k) = (field_for_y(stencil_map_y(1,stencil_idx) + k) &
                              * frac_dry_flux(idx_3d + k) ) / dt

        end do ! vertical levels k

      end do ! dof_iterator

    end if

  end subroutine ffsl_flux_xy_first_code

end module ffsl_flux_xy_first_kernel_mod
