!=main.f90
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
PROGRAM MGTRC

  USE param1,                ONLY : ver_info
  USE file_mod,              ONLY : lfile
  USE fline_mod,             ONLY : lfigout
  USE cylindrical_coord_mod, ONLY : magset,         &
    &                               free_mem_field
  USE plot_mod,              ONLY : init_plot,      &
    &                               end_plot

  IMPLICIT NONE


  ver_info = '1.0.0'

  PRINT *
  PRINT *, ' MGTRC --- MAGNETIC FIELD LINE TRACING CODE, VERSION ', ver_info, ' ---'
  PRINT *

  CALL mgcpu(0, 'start of execution  ')

  CALL mgcpu(0, 'start of execution  ')
  CALL readin

  IF(lfile)THEN
    CALL file_open
  END IF

  IF(lfigout) CALL init_plot

  CALL mgcpu(0, 'end of readin       ')
  CALL magset
  CALL mgcpu(0, 'end of magset       ')
  !CALL magadj
  CALL mgaxis
  CALL mgcpu(0, 'end of mgaxis       ')
  CALL mgltrc
  CALL mgcpu(0, 'end of mgltrc       ')
  CALL free_mem_field
  CALL mgcpu(0, 'end of free_mem     ')
  CALL mgcpu(1, 'end of execution    ')

  IF(lfigout) CALL end_plot

  IF(lfile)THEN
    CALL file_close
  END IF


END PROGRAM MGTRC
