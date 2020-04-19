!=mgltrc.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
!==Reference
!
!==Error Handlings
!
!==Known Bugs
!
!==Note
!
!==TODO
!
SUBROUTINE mgltrc

  USE fline_mod,    ONLY : mode,             &
    &                      lpoin,            &
    &                      lflxqnt
  USE vessel_mod,   ONLY : lvessel,          &
    &                      read_vessel,      &
    &                      free_mem_vessel
  USE limiter_mod,  ONLY : llimiter,         &
    &                      read_limiter,     &
    &                      free_mem_limiter
  USE divertor_mod, ONLY : ldivertor,        &
    &                      read_divertor,    &
    &                      free_mem_divertor

  IMPLICIT NONE


  IF(lvessel)   CALL read_vessel

  IF(llimiter)  CALL read_limiter

  IF(ldivertor) CALL read_divertor

  IF(lpoin)THEN

    SELECT CASE(TRIM(mode))
      CASE('2D', '2d')
        CALL trace2d
      CASE('3D', '3d')
        CALL trace3d
      CASE('RTHETA', 'rtheta')
        CALL trace2d_rtheta
      CASE('DIV_TRACE', 'div_trace')
        CALL div_trace3d
      CASE('VMEC', 'vmec')
        CALL mgvmec
      CASE DEFAULT
        CALL trace2d
    END SELECT

  END IF

  IF(lflxqnt)   CALL cal_flxqnt

  IF(lvessel)   CALL free_mem_vessel

  IF(llimiter)  CALL free_mem_limiter

  IF(ldivertor) CALL free_mem_divertor


END SUBROUTINE mgltrc
