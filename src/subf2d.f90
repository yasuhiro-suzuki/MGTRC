!> @file subf2d.f90
!------------------------------------------------------------------------------
!
! SUBROUITNE: subf2d
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> @brief
!>
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
SUBROUTINE subf2d (t, x,    &
  &                xp, iout &
  &               )

  USE kind_spec
  USE cylindrical_coord_mod, ONLY : mgval1

  IMPLICIT NONE

!Arguments
  INTEGER, INTENT(OUT) :: iout
  REAL(DP), INTENT(IN) :: t,   &
    &                     x(2)
  REAL(DP), INTENT(OUT) :: xp(2)
!Local variables
  REAL(DP) :: r,   &
    &         z,   &
    &         phi, &
    &         br,  &
    &         bp,  &
    &         bz,  &
    &         bb


  iout =  0
  r    =  x(1)
  z    =  x(2)
  phi  =  t 

!----------------------------------------
  CALL mgval1(r, phi, z, br, bp, bz, bb)
!----------------------------------------

  IF(bp == 0.0_DP)THEN
    xp   =  0.0_DP
    iout =  1
    RETURN
  END IF

  xp(1) =  r * br / bp
  xp(2) =  r * bz / bp


  RETURN
END SUBROUTINE subf2d
