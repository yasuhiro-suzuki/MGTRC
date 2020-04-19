!=limiter_mod.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
! SUBROUTINE for limiter
!
!==Reference
!
! None
!
!==Error Handlings
!
! None
!
!==Known Bugs
!
! None
!
!==Note
!
! None
!
!==TODO
!
! None
!

MODULE limiter_mod

  USE kind_spec
  USE param1,   ONLY : pi, &
    &                  pi2

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: llimiter,       &
    &        lcheck_limiter
  INTEGER, PARAMETER :: ntheta_lim =  500
  INTEGER :: nphi_lim
  CHARACTER(LEN=10) :: limiter_model = '2015a'

  INTEGER :: mtor
  INTEGER, ALLOCATABLE :: ntheta(:)
  REAL(DP) :: pi2m, &
    &         dphi
  REAL(DP), ALLOCATABLE :: phi_lim(:), &
    &                      xx(:,:),    &
    &                      yy(:,:),    &
    &                      zz(:,:),    &
    &                      rr(:,:)

  PUBLIC :: llimiter,         &
    &       lcheck_limiter,   &
    &       limiter_model,    &
    &       nphi_lim,         &
    &       phi_lim,          &
    &       read_limiter,     &
    &       free_mem_limiter, &
    &       check_limiter

CONTAINS

  SUBROUTINE make_mem_limiter

    IMPLICIT NONE


    ALLOCATE(phi_lim(nphi_lim), ntheta(nphi_lim), xx(ntheta_lim,nphi_lim), yy(ntheta_lim,nphi_lim), zz(ntheta_lim,nphi_lim), rr(ntheta_lim,nphi_lim))


  END SUBROUTINE make_mem_limiter

  SUBROUTINE free_mem_limiter

    IMPLICIT NONE


    DEALLOCATE(phi_lim, ntheta, xx, yy, zz, rr)


  END SUBROUTINE free_mem_limiter

  SUBROUTINE read_limiter

    IMPLICIT NONE

    !SELECT CASE(TRIM(limiter_model))
    !  CASE("2015a")
    !    CALL read_limiter_2015a
    !  CASE DEFAULT
        CALL read_limiter_2015a
    !END SELECT


  END SUBROUTINE read_limiter

  SUBROUTINE read_limiter_2015a

    IMPLICIT NONE
!Local variables
    INTEGER :: i, &
      &        j


    READ(35,*) mtor
    READ(35,*) nphi_lim

    pi2m =  pi2  / mtor

    CALL make_mem_limiter

    DO j=1,nphi_lim
      READ(35,*) phi_lim(j)
      READ(35,*) ntheta(j)
      IF(ntheta(j) > 500) STOP ' ntheta must be smaller than 500!'
      DO i=1,ntheta(j)
        READ(35,*) rr(i,j), zz(i,j)
      END DO
    END DO

    phi_lim(:) =  pi2 * phi_lim(:) / 360

    WHERE(phi_lim == 0.0_DP) phi_lim =  phi_lim + 1.0E-04_DP
 

  END SUBROUTINE read_limiter_2015a

  SUBROUTINE check_limiter (r, z, iphi, & !(in)
    &                       iout        & !(out)
    &                      )

    IMPLICIT NONE
!Arguments
    INTEGER, INTENT(IN) :: iphi
    INTEGER, INTENT(OUT) :: iout
    REAL(DP), INTENT(IN) :: r, &
      &                     z
!Local variables
    INTEGER :: i
    REAL(DP) :: theta0, &
      &         theta,  &
      &         theta1, &
      &         dtheta


    theta  =  0.0_DP
    theta0 =  ATAN2(zz(1,iphi) - z, rr(1,iphi) - r)
    IF(theta0 < 0.0_DP) theta0 =  theta0 + pi2
    DO i=2,ntheta(iphi)
      theta1 =  ATAN2(zz(i,iphi) - z, rr(i,iphi) - r)
      IF(theta1 < 0.0_DP) theta1 =  theta1 + pi2
      dtheta =  theta1 - theta0
      IF(dtheta >  pi) dtheta =  dtheta - pi2
      IF(dtheta < -pi) dtheta =  dtheta + pi2
      theta  =  theta + dtheta
      theta0 =  theta1
    END DO
    theta =  theta / pi2
    IF(theta >= 0.0_DP)THEN
      theta =  theta + 1.0e-06_DP
    ELSE
      theta =  theta - 1.0e-06_DP
    END IF

    iout =  0
    IF(INT(theta) == 0) iout =  1

 
  END SUBROUTINE check_limiter

END MODULE limiter_mod
