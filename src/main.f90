! Copyright (c) 2020, Yasuhiro Suzuki
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright notice,
!   this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
! * Neither the name of the <organization> nor the names of its contributors
!   may be used to endorse or promote products derived from this software
!   without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!------------------------------------------------------------------------------
!> @file main.f90
!------------------------------------------------------------------------------
! MGTRC Code
!------------------------------------------------------------------------------
!
! PROGRAM: main.f90
!
!> @author
!> Yasuhiro Suzuki, National Institute for Fusion Science
!
! DESCRIPTION:
!> @brief
!> Main driver of code.
!
! REVISION HISTORY:
!> @date 19 Apr 2020
!
!> @version Initial Version
!
!------------------------------------------------------------------------------
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
