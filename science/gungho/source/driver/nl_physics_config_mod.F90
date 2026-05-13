!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Global flag controlling whether NL physics is active.
!> @details Provides a single module-level logical that all gungho driver
!>          modules use to conditionally execute UM physics code paths.
!>          Set to .false. by default; set to .true. to activate full NL
!>          physics (e.g. for an lfric_atm configuration).
module nl_physics_config_mod

  implicit none

  logical, public, save :: use_nl_physics = .false.

end module nl_physics_config_mod
