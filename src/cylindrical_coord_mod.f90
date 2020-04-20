!> @file cylindrical_coord_mod.f90
!------------------------------------------------------------------------------
!
! MODULE: cylindrical_coord_mod
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
MODULE cylindrical_coord_mod

  USE kind_spec
  USE param1,        ONLY : pi2
  USE spline_mod,    ONLY : l2d,    &
    &                       l3d,    &
    &                       f2d,    &
    &                       f3d,    &
    &                       nx2d,   &
    &                       ny2d,   &
    &                       nx3d,   &
    &                       ny3d,   &
    &                       nz3d,   &
    &                       splin2, &
    &                       spl2df, &
    &                       splin3, &
    &                       spl3df, &
    &                       spl3dd
#ifdef NETCDF
  USE netcdf
#endif
#ifdef MPI
  USE mpi_param_mod, ONLY : myrank
  USE mpi
#endif

  IMPLICIT NONE

  PRIVATE

!
! A parameter FILE_FORMAT selects the Fortran binary, netCDF, and HDF5
!

  CHARACTER(LEN=6) :: file_format = 'bin'

!
! A parameter MAG_FILE specifies the file name of the magnetic field
!

  CHARACTER(LEN=200) :: mag_file

!----------------------------------------------------------------------------
! A parameter MAG_FORM prescribes the table format of the magnetic field.
! A parameter VERSION prescribes the version for the old legacy format.
!
! 1. mgrid:   file format used in the MAKEGRID code
! 2. mgo:     file format used in KMAG/KMAG2 codes
! 3. mag:     file format for the equilibrium field used in HINT/HINT2 codes
! 3.1. ver2:  legacy HINT2 format
! 4. vac:     file format for the vacuum field used in HINT/HINT2 codes
! 4.1. ver2:  legacy HINT2 format
! 5. movie:   file format used in the MIPS code
! 6. mips_eq: file format for the equilibrium field used in the MIPS code
!----------------------------------------------------------------------------

  CHARACTER(LEN=20) :: mag_form = '', &
    &                  version  = ''

!
! A parameter LFLUX switch on/off to read the flux distribution.
!

  LOGICAL :: lflux

!
! A parameter FLX_FILE specifies the file name of the magnetic field
!

  CHARACTER(LEN=200) :: flx_file

!
! A parameter FLX_FORM prescribes the table format of the flux distribution
!
! 1. xss:    file format used in the SMAP code
! 1. flx:    file format used in the HINT code
! 2. eqdsk:  file format used in the EFIT code
! 3. eqdata: file format used in the TOPICS
!

  CHARACTER(LEN=20) :: flx_form = 'xss'

!
! A parameter LVACOUT switch on/off to output magnetic field
!

  LOGICAL :: lvacout

!

  LOGICAL :: lsymmetry

! 

  INTEGER ::  mtor,                    & ! toroidal field period
    &         nr0b,                    & ! grid number along R-direction
    &         nt0b,                    & ! grid number along phi-direction
    &         nz0b,                    & ! grid number along Z-direction
!
    &         ipfcoil(500) =  0,       & ! index of PF coil to specify the magnetic axis
!                                        ! CAUTION! array size is 500 (max)
!
    &         kstep        =  99999,   & ! time steps of MIPS
    &         igrid(4)                   ! grid number of MIPS
  REAL(DP) :: bmax         =  5.0_DP,  & ! torelance of maximum B field
    &         sedge        =  0.98_DP, & ! the edge toroidal flux (default)
    &         pi2m,                    & ! one toroidal field period [rad]
    &         rmaxb,                   & ! major radius of the inner side of computational DOmain [m]
    &         rminb,                   & ! major radius of the outer side of computational DOmain [m]
    &         zmaxb,                   & ! height of the lower side of computational DOmain [m]
    &         zminb,                   & ! height of the higher side of the computational DOmain [m]
    &         delrb,                   & ! grid size along R-direction [m]
    &         delzb,                   & ! grid size along Z-direction [m]
    &         badjust      =  1.0_DP,  & ! normalization factor of magnetic field
    &         bnorm        =  3.0_DP,  & ! normalization factor of magnetic field
    &         cj(500)      =  0.0_DP,  & ! work array to store coil current in each coils [A]
!                                        ! CAUTION! array size is 500 (max)
    &         extcur(500)  =  0.0_DP,  & ! work array to store coil current in each coils [A]
!                                        ! #NOTE# compatibility for MAKEGRID of VMEC
!                                        ! CAUTION! array size is 500 (max)
    &         cturn(500)   =  1.0_DP,  & ! work array to store coil current in each coils [A]
!                                        ! CAUTION! array size is 500 (max)
    &         cfact(500)   =  1.0_DP,  & ! work array to store coil current in each coils [A]
!                                        ! CAUTION! array size is 500 (max)
    &         cpfcoil(500) =  1.0_DP,  &
    &         mbound(4)    =  0.0_DP     ! boundary of MIPS
!                                        ! mbound(1) = Rmin [m]
!                                        ! mbound(2) = Rmax [m]
!                                        ! mbound(3) = Zmin [m]
!                                        ! mbound(4) = Zmax [m]
  REAL(DP), ALLOCATABLE :: rg(:), &
    &                      zg(:)


  NAMELIST /nlinp_coil_dat/ file_format, &
    &                       mag_file,    &
    &                       flx_file,    &
    &                       mag_form,    &
    &                       flx_form,    &
    &                       version,     &
    &                       lflux,       &
    &                       lsymmetry,   &
    &                       lvacout,     &
    &                       cj,          &
    &                       extcur,      &
    &                       cturn,       &
    &                       cfact,       &
    &                       sedge,       &
    &                       badjust,     &
    &                       ipfcoil,     &
    &                       cpfcoil,     &
    &                       bnorm,       &
    &                       igrid,       &
    &                       mbound,      &
    &                       kstep

  PUBLIC :: file_format,    &
    &       mag_file,       &
    &       flx_file,       &
    &       mag_form,       &
    &       flx_form,       &
    &       version,        &
    &       lflux,          &
    &       lsymmetry,      &
    &       lvacout,        &
    &       cj,             &
    &       extcur,         &
    &       cturn,          &
    &       cfact,          &
    &       sedge,          &
    &       badjust,        &
    &       ipfcoil,        &
    &       cpfcoil,        &
    &       bnorm,          &
    &       igrid,          &
    &       mbound,         &
    &       kstep,          &
    &       pi2m,           &
    &       mtor,           &
    &       nr0b,           &
    &       nt0b,           &
    &       nz0b,           &
    &       rminb,          &
    &       rmaxb,          &
    &       zminb,          &
    &       zmaxb,          &
    &       delrb,          &
    &       delzb,          &
    &       rg,             &
    &       zg,             &
    &       mag_file_open,  &
    &       mag_file_close, &
    &       free_mem_field, &
    &       magset,         &
    &       mgval1,         &
    &       mgval2,         &
    &       mgval3,         &
    &       nlinp_coil_dat

CONTAINS


  SUBROUTINE mag_file_open

    IMPLICIT NONE


    OPEN(25, FILE=mag_file, FORM='unformatted', STATUS='old')

    IF(lflux) OPEN(26, FILE=flx_file, FORM='unformatted', STATUS='old')


  END SUBROUTINE mag_file_open

  SUBROUTINE mag_file_close

    IMPLICIT NONE


    CLOSE(25)

    IF(lflux) CLOSE(26)


  END SUBROUTINE mag_file_close

  SUBROUTINE read_field

    IMPLICIT NONE


    l3d =  4
    IF(lflux)THEN
      l3d =  5
    END IF

    SELECT CASE(TRIM(mag_form))
      CASE('mgrid', "MGRID")
#ifdef NETCDF
        SELECT CASE(TRIM(file_format))
          CASE('bin',"BIN")
#endif
            CALL read_mgrid_bin
#ifdef NETCDF
          CASE('netcdf',"NETCDF")
            CALL read_mgrid_nc
          END SELECT
#endif
      CASE('mgo', "MGO")
        CALL read_mgo
      CASE('mag', "MAG")
#ifdef NETCDF
        SELECT CASE(TRIM(file_format))
          CASE('bin',"BIN")
#endif
            CALL read_hint_eq
#ifdef NETCDF
          CASE('netcdf',"NETCDF")
            CALL read_hint_eq
          END SELECT
#endif
      CASE('vac', "VAC")
#ifdef NETCDF
        SELECT CASE(TRIM(file_format))
          CASE('bin',"BIN")
#endif
            CALL read_hint_vac
#ifdef NETCDF
          CASE('netcdf',"NETCDF")
            CALL read_hint_vac
          END SELECT
#endif
      CASE('mips', "MIPS")
        CALL read_mips
      CASE('mips_eq', "MIPS_EQ")
        CALL read_mips_eq
      CASE DEFAULT
        CALL read_mgrid_bin
    END SELECT

    IF(lflux)THEN
      SELECT CASE(TRIM(flx_form))
        CASE('xss')
          CALL read_xss
        CASE('flx')
#ifdef NETCDF
          SELECT CASE(TRIM(file_format))
            CASE('bin',"BIN")
#endif
              CALL read_flx
#ifdef NETCDF
            CASE('netcdf',"NETCDF")
              CALL read_flx
          END SELECT
#endif
        CASE('eqdsk')
          CALL read_eqdsk
        CASE('eqdata')
          CALL read_eqdata
        CASE DEFAULT
          CALL read_xss
      END SELECT
    END IF


  END SUBROUTINE read_field

  SUBROUTINE read_mgrid_bin

    IMPLICIT NONE

    INTEGER, SAVE :: isave =  1
    INTEGER :: nextcur, &
      &        i,       &
      &        j,       &
      &        k,       &
      &        n
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=30), ALLOCATABLE :: curlabel(:)
    CHARACTER(LEN=100) :: fmt
    LOGICAL :: lstyle2000
#ifdef MPI
!For MPI
    INTEGER :: icode, &
      &        ierr
#endif

    REWIND 25

    READ(25) nr0b, nz0b, nt0b, mtor, nextcur
    READ(25) rminb, zminb, rmaxb, zmaxb

    IF(nextcur < 0) lstyle2000 = .true.
    nextcur = ABS(nextcur)

#ifdef MPI
    IF(nextcur > 500)THEN
      PRINT *, ' nextcur > 500'
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
    END IF
#else
    IF(nextcur > 500) STOP ' nextcur > 500'
#endif

    ALLOCATE(curlabel(nextcur))

    IF(isave == 1)THEN
      ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))
      isave =  0
    END IF

    f3d(:,:,:,:) =  0.0_DP

    READ(25) (curlabel(i), i=1,nextcur)

    ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

    fmt = '(I7, ES17.6, I7, ES12.4)'

#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *, ' Reading coil data'
    PRINT *

    PRINT *
    PRINT *, ' coil No    c_I [A/T]      trun  c_factor '
    PRINT *, '--------------------------------------------'

    IF((SUM(cj(:)) == 0.0_DP) .AND. (SUM(extcur(:)) /= 0.0_DP)) cj(:) =  extcur(:)

    DO n=1,nextcur
      PRINT fmt, n, cj(n), INT(cturn(n)), cfact(n)
    END DO

    PRINT *
#ifdef MPI
    END IF
#endif

    fmt = '(A10, I3, A15, A5, F12.6, A6)'

    loop010 : DO n=1,nextcur 
#ifdef MPI
      IF(myrank == 0)THEN
#endif
      PRINT fmt, ' coil NO: ', n, TRIM(curlabel(n)), ' c_I ', cj(n) * cfact(n) * cturn(n) * badjust * 1.0E-06_DP, ' [MA]'
#ifdef MPI
      END IF
#endif
      IF(lstyle2000)THEN
        READ(25) br, bp, bz
      ELSE
        READ(25) (((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)
      END IF
      f3d(1,1:nr0b,1:nz0b,1:nt0b) =  f3d(1,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * br(1:nr0b,1:nz0b,1:nt0b) !< B_R
      f3d(2,1:nr0b,1:nz0b,1:nt0b) =  f3d(2,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * bp(1:nr0b,1:nz0b,1:nt0b) !< B_phi
      f3d(3,1:nr0b,1:nz0b,1:nt0b) =  f3d(3,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * bz(1:nr0b,1:nz0b,1:nt0b) !< B_Z
    END DO loop010

    IF(lvacout)THEN
      WRITE(35) nr0b, nz0b, nt0b, mtor
      WRITE(35) rminb, zminb, rmaxb, zmaxb
      WRITE(35) (((f3d(1,i,j,k), f3d(2,i,j,k), f3d(3,i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)
    END IF 


    DEALLOCATE(br, bp, bz)
    DEALLOCATE(curlabel)


  END SUBROUTINE read_mgrid_bin

#ifdef NETCDF
  SUBROUTINE read_mgrid_nc

    IMPLICIT NONE

    INTEGER, SAVE :: isave =  1
    INTEGER :: nextcur, &
      &        n
    REAL(DP), ALLOCATABLE :: raw_coil_current(:), &
      &                      br(:,:,:),           &
      &                      bp(:,:,:),           &
      &                      bz(:,:,:)
    CHARACTER(LEN=30), ALLOCATABLE :: curlabel(:)
    CHARACTER(LEN=100) :: fmt
!For netCDF
    INTEGER :: ncid,          &
      &        ir_varid,      &
      &        jz_varid,      &
      &        kp_varid,      &
      &        mtor_varid,    &
      &        rmin_varid,    &
      &        rmax_varid,    &
      &        zmin_varid,    &
      &        zmax_varid,    &
      &        nextcur_varid, &
      &        coilgrp_varid, &
      &        coilcur_varid, &
      &        br_varid,      &
      &        bp_varid,      &
      &        bz_varid,      &
      &        x_dimid,       &
      &        y_dimid,       &
      &        z_dimid,       &
      &        dim3ids(3)
    REAL(DP), ALLOCATABLE :: w3d(:,:,:)
    CHARACTER(LEN=100) :: vname
#ifdef MPI
!For MPI
    INTEGER :: icode, &
      &        ierr
#endif


    CALL check(NF90_OPEN(mag_file, NF90_NOWRITE, ncid))

    CALL check(NF90_INQ_VARID(ncid, 'ir',      ir_varid))
    CALL check(NF90_INQ_VARID(ncid, 'jz',      jz_varid))
    CALL check(NF90_INQ_VARID(ncid, 'kp',      kp_varid))
    CALL check(NF90_INQ_VARID(ncid, 'nfp',     mtor_varid))
    CALL check(NF90_INQ_VARID(ncid, 'nextcur', nextcur_varid))

    CALL check(NF90_GET_VAR(ncid, ir_varid,      nr0b))
    CALL check(NF90_GET_VAR(ncid, jz_varid,      nz0b))
    CALL check(NF90_GET_VAR(ncid, kp_varid,      nt0b))
    CALL check(NF90_GET_VAR(ncid, mtor_varid,    mtor))
    CALL check(NF90_GET_VAR(ncid, nextcur_varid, nextcur))

    CALL check(NF90_INQ_VARID(ncid, 'rmin',  rmin_varid))
    CALL check(NF90_INQ_VARID(ncid, 'rmax',  rmax_varid))
    CALL check(NF90_INQ_VARID(ncid, 'zmin',  zmin_varid))
    CALL check(NF90_INQ_VARID(ncid, 'zmax',  zmax_varid))

    CALL check(NF90_GET_VAR(ncid, rmin_varid, rminb))
    CALL check(NF90_GET_VAR(ncid, rmax_varid, rmaxb))
    CALL check(NF90_GET_VAR(ncid, zmin_varid, zminb))
    CALL check(NF90_GET_VAR(ncid, zmax_varid, zmaxb))

#ifdef MPI
    IF(nextcur > 500)THEN
      PRINT *, ' nextcur > 500'
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
    END IF
#else
    IF(nextcur > 500) STOP ' nextcur > 500'
#endif

    ALLOCATE(curlabel(nextcur), raw_coil_current(nextcur))

    CALL check(NF90_INQ_VARID(ncid, 'coil_group',    coilgrp_varid))
    CALL check(NF90_INQ_VARID(ncid, 'raw_coil_cur',  coilcur_varid))

    CALL check(NF90_GET_VAR(ncid, coilgrp_varid, curlabel))
    CALL check(NF90_GET_VAR(ncid, coilcur_varid, raw_coil_current))

    IF(isave == 1)THEN
      ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))
      isave =  0
    END IF

    f3d(:,:,:,:) =  0.0_DP

    ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

    fmt = '(I7, ES17.6, I7, ES12.4)'

#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *, ' Reading coil data'
    PRINT *

    PRINT *
    PRINT *, ' coil No    c_I [A/T]      trun  c_factor '
    PRINT *, '--------------------------------------------'

    IF((SUM(cj(:)) == 0.0_DP) .AND. (SUM(extcur(:)) /= 0.0_DP)) cj(:) =  extcur(:)

    DO n=1,nextcur
      PRINT fmt, n, cj(n), INT(cturn(n)), cfact(n)
    END DO

    PRINT *
#ifdef MPI
    END IF
#endif

    fmt = '(A10, I3, A15, A5, F12.6, A6)'

    loop010 : DO n=1,nextcur
#ifdef MPI
      IF(myrank == 0)THEN
#endif
      PRINT fmt, ' coil NO: ', n, TRIM(curlabel(n)), ' c_I ', cj(n) * cfact(n) * cturn(n) * badjust * 1.0E-06_DP, ' [MA]'
#ifdef MPI
      END IF
#endif

      WRITE(vname, '(I3.3)') n
      vname =  'br_' // vname

      CALL check(NF90_INQ_VARID(ncid, TRIM(vname),  br_varid))

      CALL check(NF90_GET_VAR(ncid, br_varid, br))

      WRITE(vname, '(I3.3)') n
      vname =  'bp_' // vname

      CALL check(NF90_INQ_VARID(ncid, TRIM(vname),  bp_varid))

      CALL check(NF90_GET_VAR(ncid, bp_varid, bp))

      WRITE(vname, '(I3.3)') n
      vname =  'bz_' // vname

      CALL check(NF90_INQ_VARID(ncid, TRIM(vname),  bz_varid))

      CALL check(NF90_GET_VAR(ncid, bz_varid, bz))

      f3d(1,1:nr0b,1:nz0b,1:nt0b) =  f3d(1,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * br(1:nr0b,1:nz0b,1:nt0b) !< B_R
      f3d(2,1:nr0b,1:nz0b,1:nt0b) =  f3d(2,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * bp(1:nr0b,1:nz0b,1:nt0b) !< B_phi
      f3d(3,1:nr0b,1:nz0b,1:nt0b) =  f3d(3,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * bz(1:nr0b,1:nz0b,1:nt0b) !< B_Z
    END DO loop010

    CALL check(NF90_CLOSE(ncid))

    IF(lvacout)THEN

      CALL check(NF90_CREATE('vacfile.nc', NF90_CLOBBER, ncid))

      CALL check(NF90_DEF_VAR(ncid, "mtor",  NF90_INT,   mtor_varid))

      CALL check(NF90_DEF_VAR(ncid, "rminb", NF90_REAL8, rmin_varid))
      CALL check(NF90_DEF_VAR(ncid, "rmaxb", NF90_REAL8, rmax_varid))
      CALL check(NF90_DEF_VAR(ncid, "zminb", NF90_REAL8, zmin_varid))
      CALL check(NF90_DEF_VAR(ncid, "zmaxb", NF90_REAL8, zmax_varid))

      CALL check(NF90_DEF_DIM(ncid, "R",   nr0b,   x_dimid))
      CALL check(NF90_DEF_DIM(ncid, "Z",   nz0b,   y_dimid))
      CALL check(NF90_DEF_DIM(ncid, "phi", nt0b, z_dimid))

      dim3ids =  (/x_dimid, y_dimid, z_dimid/)

      CALL check(NF90_DEF_VAR(ncid, "B_R",   NF90_REAL8, dim3ids, br_varid))
      CALL check(NF90_DEF_VAR(ncid, "B_phi", NF90_REAL8, dim3ids, bp_varid))
      CALL check(NF90_DEF_VAR(ncid, "B_Z",   NF90_REAL8, dim3ids, bz_varid))

      CALL check(NF90_ENDDEF(ncid))

      CALL check(NF90_PUT_VAR(ncid, mtor_varid, mtor))

      CALL check(NF90_PUT_VAR(ncid, rmin_varid, rminb))
      CALL check(NF90_PUT_VAR(ncid, rmax_varid, rmaxb))
      CALL check(NF90_PUT_VAR(ncid, zmin_varid, zminb))
      CALL check(NF90_PUT_VAR(ncid, zmax_varid, zmaxb))

      ALLOCATE(w3d(nr0b,nz0b,nt0b))

      w3d(1:nr0b,1:nz0b,1:nt0b) =  f3d(1,1:nr0b,1:nz0b,1:nt0b)

      CALL check(NF90_PUT_VAR(ncid, br_varid, w3d))

      w3d(1:nr0b,1:nz0b,1:nt0b) =  f3d(2,1:nr0b,1:nz0b,1:nt0b)

      CALL check(NF90_PUT_VAR(ncid, bp_varid, w3d))

      w3d(1:nr0b,1:nz0b,1:nt0b) =  f3d(3,1:nr0b,1:nz0b,1:nt0b)

      CALL check(NF90_PUT_VAR(ncid, bz_varid, w3d))

      CALL check(NF90_CLOSE(ncid))

      DEALLOCATE(w3d)

    END IF


    DEALLOCATE(br, bp, bz)
    DEALLOCATE(curlabel, raw_coil_current)


  END SUBROUTINE read_mgrid_nc
#endif

  SUBROUTINE read_mgo

    IMPLICIT NONE

    INTEGER :: i, &
      &        j, &
      &        k
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=100) :: fmt


    READ(25) nr0b, nz0b, nt0b, mtor
    READ(25) rminb, zminb, rmaxb, zmaxb


    ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

    f3d(:,:,:,:) =  0.0_DP

    ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

#ifdef MPI
    IF(myrank == 0)THEN
#endif

    PRINT *, ' Reading magnetic field data as MGO format'
    PRINT *

#ifdef MPI
    END IF
#endif

    loop010 : DO k=1,nt0b
      READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
    END DO loop010

    f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b) !< B_R
    f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b) !< B_phi
    f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b) !< B_Z

    DEALLOCATE(br, bp, bz)


  END SUBROUTINE read_mgo

  SUBROUTINE read_hint_eq

    IMPLICIT NONE

    INTEGER :: nt0bh,   &
      &        nt1b,    &
      &        kstep,   &
      &        i,       &
      &        j,       &
      &        jq,      &
      &        k,       &
      &        itemp(4)
    REAL(DP) :: time,   &
      &         kpitch, &
      &         dt,     &
      &         temp
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=100) :: fmt


    SELECT CASE(TRIM(version))

      CASE('ver2')

        READ(25) time
        READ(25) nr0b, nz0b, nt1b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb, kpitch

        IF(kpitch == 0.5_DP)THEN
          PRINT *
          PRINT *,' !!! kpitch=0.5 !!!: assumed helical symmetry.'
          PRINT *
          nt0bh =  nt1b - 4
          nt0b  =  2 * (nt1b - 5)
        ELSE
          nt0b  =  nt1b - 5
          mtor  =  mtor / kpitch
        END IF

        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt1b), bp(nr0b,nz0b,nt1b), bz(nr0b,nz0b,nt1b))

#ifdef MPI
        IF(myrank == 0)THEN
#endif

        PRINT *, ' Reading magnetic field data as HINT2 format'
        PRINT *

#ifdef MPI
        END IF
#endif

        REWIND 25

 1      CONTINUE

        READ(25,END=2) time
        READ(25) (itemp(i), i=1,4)
        READ(25) temp, temp, temp, temp, temp
        READ(25) temp

        fmt =  '(A7,ES12.4)' 

#ifdef MPI
        IF(myrank == 0)THEN
#endif

        PRINT fmt, ' time= ', time

#ifdef MPI
        END IF
#endif
 
        loop010 : DO k=1,nt1b
          READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
          READ(25) ((temp, temp, temp, i=1,nr0b), j=1,nz0b)
          READ(25) ((temp, i=1,nr0b), j=1,nz0b)
        END DO loop010

        GOTO 1

 2      CONTINUE

        IF(kpitch == 0.5_DP)THEN
          f3d(1,1:nr0b,1:nz0b,1:nt0bh) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0bh) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0bh) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
          DO k=1,nt0bh-1
            DO j=1,nz0b
              jq =  nz0b - j + 1
              DO i=1,nr0b
                f3d(1,i,j,nt0b+2-k) = -f3d(1,i,jq,k)
                f3d(2,i,j,nt0b+2-k) =  f3d(2,i,jq,k)
                f3d(3,i,j,nt0b+2-k) =  f3d(3,i,jq,k)
              END DO
            END DO
          END DO
        ELSE
          f3d(1,1:nr0b,1:nz0b,1:nt0b+1) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0b+1) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0b+1) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
        END IF

      CASE DEFAULT

        READ(25) kstep
        READ(25) time
        READ(25) nr0b, nz0b, nt0b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb


        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

#ifdef MPI
        IF(myrank == 0)THEN
#endif

        PRINT *, ' Reading magnetic field data as HINT format'
        PRINT *

#ifdef MPI
        END IF
#endif

        REWIND 25

 3      CONTINUE

        READ(25,END=4) kstep
        READ(25) time
        READ(25) (itemp(i), i=1,4)
        READ(25) temp, temp, temp, temp

        fmt =  '(A7,ES12.4)'

#ifdef MPI
        IF(myrank == 0)THEN
#endif

        PRINT fmt, ' time= ', time

#ifdef MPI
        END IF
#endif

        READ(25) (((br(i,j,k), bp(i,j,k), bz(i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)
        READ(25) (((temp, temp, temp, i=1,nr0b), j=1,nz0b), k=1,nt0b)
        READ(25) (((temp, i=1,nr0b), j=1,nz0b), k=1,nt0b)

        GOTO 3

 4      CONTINUE

        f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b)
        f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b)
        f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b)

    END SELECT

    DEALLOCATE(br, bp, bz)


  END SUBROUTINE read_hint_eq

  SUBROUTINE read_hint_vac

    IMPLICIT NONE

    INTEGER :: nt0bh, &
      &        nt1b,  &
      &        i,     &
      &        j,     &
      &        jq,    &
      &        k
    REAL(DP) :: kpitch
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=100) :: fmt


    SELECT CASE(TRIM(version))

      CASE('ver2')

        READ(25) nr0b, nz0b, nt1b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb, kpitch

        IF(kpitch == 0.5_DP)THEN
#ifdef MPI
          IF(myrank == 0)THEN
#endif
          PRINT *
          PRINT *,' !!! kpitch=0.5 !!!: assumed helical symmetry.'
          PRINT *
#ifdef MPI
          END IF
#endif
          nt0bh =  nt1b - 4
          nt0b  =  2 * (nt1b - 5)
        ELSE
          nt0b  =  nt1b - 5
          mtor  =  mtor / kpitch
        END IF

        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt1b), bp(nr0b,nz0b,nt1b), bz(nr0b,nz0b,nt1b))

#ifdef MPI
        IF(myrank == 0)THEN
#endif

        PRINT *, ' Reading magnetic field data as legacy HINT2 format'
        PRINT *

#ifdef MPI
        END IF
#endif

        loop010 : DO k=1,nt1b
          READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
        END DO loop010

        IF(kpitch == 0.5_DP)THEN
          f3d(1,1:nr0b,1:nz0b,1:nt0bh) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0bh) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0bh) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
          DO k=1,nt0bh-1
            DO j=1,nz0b
              jq =  nz0b - j + 1
              DO i=1,nr0b
                f3d(1,i,j,nt0b+2-k) = -f3d(1,i,jq,k)
                f3d(2,i,j,nt0b+2-k) =  f3d(2,i,jq,k)
                f3d(3,i,j,nt0b+2-k) =  f3d(3,i,jq,k)
              END DO
            END DO
          END DO
        ELSE
          f3d(1,1:nr0b,1:nz0b,1:nt0b+1) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0b+1) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0b+1) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
        END IF

      CASE DEFAULT

        READ(25) nr0b, nz0b, nt0b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb

        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

#ifdef MPI
        IF(myrank == 0)THEN
#endif

        PRINT *, ' Reading vacuum magnetic field data as HINT format'
        PRINT *

#ifdef MPI
        END IF
#endif

        READ(25) (((br(i,j,k), bp(i,j,k), bz(i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)

        f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b)
        f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b)
        f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b)

    END SELECT

    DEALLOCATE(br, bp, bz)


  END SUBROUTINE read_hint_vac

  SUBROUTINE read_mips_eq

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        jtemp, &
      &        ktemp
    REAL(DP) :: pminb,  &
      &         pmaxb,  &
      &         dtempr, &
      &         dtempz, &
      &         dtempt
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:), &
      &                      p(:,:,:)


    mtor =  igrid(4)
    IF(mtor <= 0) mtor =  1

    nr0b =  igrid(1)
    nz0b =  igrid(2)
    nt0b =  igrid(3) - 4

    ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

    f3d(:,:,:,:) =  0.0_DP

    ALLOCATE(br(nr0b,nz0b,nt0b+4), bp(nr0b,nz0b,nt0b+4), bz(nr0b,nz0b,nt0b+4), p(nr0b,nz0b,nt0b+4))

#ifdef MPI
    IF(myrank == 0)THEN
#endif

    PRINT *, ' Reading magnetic field data as MIPS format'
    PRINT *

#ifdef MPI
    END IF
#endif

    READ(25) itemp, jtemp, ktemp,    &
      &      rminb, rmaxb,           &
      &      zminb, zmaxb,           &
      &      pminb, pmaxb,           &
      &      dtempr, dtempz, dtempt, &
      &      br, bz, bp, p

    f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,3:nt0b+2)
    f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,3:nt0b+2)
    f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,3:nt0b+2)

    DEALLOCATE(br, bp, bz, p)


  END SUBROUTINE read_mips_eq

  SUBROUTINE read_mips

    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(DP) :: temp, bmax
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)


    bmax =  2.333413159206456_DP
    !bmax =  2.18529336349835_DP
    !bmax =  2.891266292494260_DP

    READ(25) nr0b, nz0b, nt0b, mtor
    READ(25) rminb, zminb, rmaxb, zmaxb

    ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

    f3d(:,:,:,:) =  0.0_DP

    ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

#ifdef MPI
    IF(myrank == 0)THEN
#endif

    PRINT *, ' Reading vacuum magnetic field data as HINT format'
    PRINT *

#ifdef MPI
    END IF
#endif

    DO k=1,nt0b
      READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
      READ(25) ((temp, temp, temp, i=1,nr0b), j=1,nz0b)
      READ(25) ((temp, temp, i=1,nr0b), j=1,nz0b)
    END DO

    f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b) * bmax
    f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b) * bmax
    f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b) * bmax


    DEALLOCATE(br, bp, bz)


    RETURN
  END SUBROUTINE read_mips

  SUBROUTINE read_xss

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        jtemp, &
      &        ktemp, &
      &        mtemp, &
      &        i,     &
      &        j,     &
      &        k
    REAL(DP) :: r1temp, &
      &         r2temp, &
      &         z1temp, &
      &         z2temp
    REAL(DP), ALLOCATABLE :: ss(:,:,:)
    CHARACTER(LEN=100) :: fmt
#ifdef MPI
!For MPI
    INTEGER :: icode, &
      &        ierr
#endif


#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE READ_XSS                                    '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *
#ifdef MPI
    END IF
#endif

    READ(26) itemp, jtemp, ktemp, mtemp
    READ(26) r1temp, z1temp, r2temp, z2temp

    IF(itemp /= nr0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  itemp =  ', itemp,  '  nr0b  = ', nr0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(jtemp /= nz0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  jtemp =  ', jtemp,  '  nz0b  = ', nz0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(ktemp /= nt0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  ktemp =  ', ktemp,  '  nt0b  = ', nt0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(mtemp /= mtor)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  mtemp =  ', mtemp,  '  mtor  = ', mtor
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(r1temp /= rminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  r1temp = ', r1temp, '  rminb = ', rminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(r2temp /= rmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  r2temp = ', r2temp, '  rmaxb = ', rmaxb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(z1temp /= zminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  z1temp = ', z1temp, '  zminb = ', zminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(z2temp /= zmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_xss: parameter error  z2temp = ', z2temp, '  zmaxb = ', zmaxb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    ALLOCATE(ss(nr0b,nz0b,nt0b))

    ss(:,:,:) =  0.0_DP

#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *, ' Reading normalized toroidal flux as XSS format'
    PRINT *
#ifdef MPI
    END IF
#endif

    loop010 : DO k=1,nt0b
      READ(26) ((ss(i,j,k), i=1,nr0b), j=1,nz0b)
    END DO loop010

    f3d(5,1:nr0b,1:nz0b,1:nt0b) =  ss(1:nr0b,1:nz0b,1:nt0b)


    DEALLOCATE(ss)


  END SUBROUTINE read_xss

  SUBROUTINE read_flx

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        jtemp, &
      &        ktemp, &
      &        mtemp, &
      &        i,     &
      &        j,     &
      &        k
    REAL(DP) :: r1temp, &
      &         r2temp, &
      &         z1temp, &
      &         z2temp
    REAL(DP), ALLOCATABLE :: ss(:,:,:)
    CHARACTER(LEN=100) :: fmt
#ifdef MPI
!For MPI
    INTEGER :: icode, &
      &        ierr
#endif


#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE READ_FLX                                    '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *
#ifdef MPI
    END IF
#endif

    READ(26) itemp, jtemp, ktemp, mtemp
    READ(26) r1temp, z1temp, r2temp, z2temp

    IF(itemp /= nr0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  itemp =  ', itemp,  '  nr0b  = ', nr0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(jtemp /= nz0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  jtemp =  ', jtemp,  '  nz0b  = ', nz0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(ktemp /= nt0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  ktemp =  ', ktemp,  '  nt0b  = ', nt0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(mtemp /= mtor)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  mtemp =  ', mtemp,  '  mtor  = ', mtor
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(r1temp /= rminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  r1temp = ', r1temp, '  rminb = ', rminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(r2temp /= rmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  r2temp = ', r2temp, '  rmaxb = ', rmaxb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(z1temp /= zminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  z1temp = ', z1temp, '  zminb = ', zminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(z2temp /= zmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_flx: parameter error  z2temp = ', z2temp, '  zmaxb = ', zmaxb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    ALLOCATE(ss(nr0b,nz0b,nt0b))

    ss(:,:,:) =  0.0_DP

#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *, ' Reading normalized toroidal flux as HINT format'
    PRINT *
#ifdef MPI
    END IF
#endif

    READ(26) (((ss(i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)

    f3d(5,1:nr0b,1:nz0b,1:nt0b) =  ss(1:nr0b,1:nz0b,1:nt0b)


    DEALLOCATE(ss)


  END SUBROUTINE read_flx

  SUBROUTINE read_eqdsk

    IMPLICIT NONE

    INTEGER :: nr, &
      &        nz, &
      &        ns, &
      &        i,  &
      &        j,  &
      &        k
    REAL(DP) :: rdim,   &
      &         zdim,   &
      &         rleft,  &
      &         zmid,   &
      &         psi0,   &
      &         psia,   &
      &         rmin,   &
      &         rmax,   &
      &         zmin,   &
      &         zmax,   &
      &         r,      &
      &         z,      &
      &         dr,     &
      &         dz,     &
      &         eps,    &
      &         temp,   &
      &         x2d(2), &
      &         w2d(1)
    REAL(DP), ALLOCATABLE :: psig(:,:)
    CHARACTER(LEN=10) :: header(6)
    CHARACTER(LEN=100) :: fmt
#ifdef MPI
!For MPI
    INTEGER :: icode, &
      &        ierr
#endif


#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE READ_EQDSK                                  '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *
#ifdef MPI
    END IF
#endif

    READ(26,'( 6A8, 3I4 )') (header(i), i=1,6), j, nr, nz

    ns =  nr

    READ(26,'( 5E16.9 )') rdim, zdim, temp, rleft, zmid
    READ(26,'( 5E16.9 )') temp, temp, psi0, psia, temp
    READ(26,'( 5E16.9 )') temp, temp, temp, temp, temp
    READ(26,'( 5E16.9 )') temp, temp, temp, temp, temp

    rmin =  rleft
    rmax =  rleft + rdim
    zmin =  zmid - 0.5_DP * zdim
    zmax =  zmid + 0.5_DP * zdim

#ifdef MPI
    IF(myrank == 0)THEN
#endif

    PRINT *
    PRINT *,  '   COORDINATES SYSTEMS '
    PRINT *,  '  --------------------------------------------------------'
    PRINT *,  '       RMIN          RMAX          ZMIN          ZMAX '

    fmt = '(2X,4F14.9)'

    PRINT fmt, rmin, rmax, zmin, zmax

    PRINT *
    PRINT *,  '   RESOLUSIONS'
    PRINT *,  '  ------------------------'

    PRINT *,  '      NR      NZ      NS '

    fmt = '(2X,3I8)'

    PRINT fmt, nr, nz, ns

    PRINT *
    PRINT *,  '   EQUILIBRIUM PARAMETRS '
    PRINT *,  '  ----------------------------'
    PRINT *,  '       PSI0          PSIA '

    fmt = '(2x,2F14.9)'

    PRINT fmt, psi0, psia

    PRINT *
    PRINT *

#ifdef MPI
    END IF
#endif

    ALLOCATE(psig(nr,nz))

    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') ((psig(i,j), i=1,nr), j=1,nz)


    IF(rmin /= rminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdsk: parameter error  rleft       =', rmin, '  rminb = ', rminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(rmax /= rmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdsk: parameter error  rleft+rdim  =', rmax, '  rmaxb = ', rmaxb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(zmin /= zminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdsk: parameter error  zmid-zdim/2 =', zmin, '  zminb = ', zminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(zmax /= zmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdsk: parameter error  zmid+zdim/2 =', zmax, '  zminb = ', zminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *, ' Reading normalized poloidal flux as EQDSK format'
    PRINT *
#ifdef MPI
    END IF
#endif

    !psig(:,:) =  psig(:,:) - psia
    !psig(:,:) =  psig(:,:) / (psia - psi0) + 1.0_DP
    psig(:,:) = 1.0_DP - (psig(:,:) - psia) / (psi0 - psia)

    IF((nr /= nr0b) .OR. (nz /= nz0b))THEN

#ifdef MPI
      IF(myrank == 0)THEN
#endif
      PRINT *
      PRINT *
      PRINT *
#ifdef MPI
      END IF
#endif

      l2d =  1

      ALLOCATE(rg(nr), zg(nz), f2d(l2d,-2:nr+3,-2:nz+3))

      f2d(1,1:nr,1:nz) =  psig(1:nr,1:nz)

      dr = (rmaxb - rminb) / (nr - 1)
      dz = (zmaxb - zminb) / (nz - 1)

      DO i=1,nr
        rg(i) =  rminb + dr * (i - 1)
      END DO

      DO j=1,nz
        zg(j) =  zminb + dz * (j - 1)
      END DO

#ifdef MPI
      IF(myrank == 0)THEN
#endif
      DO j=1,nz
        DO i=1,nr
          WRITE(71,'(10ES15.7)') rg(i), zg(j), psig(i,j)
        END DO
        WRITE(71,*)
      END DO
#ifdef MPI
      END IF
#endif

      DO j=1,nz
        DO i=1,3
          r =  rg(1)  - dr * (4 - i)
          CALL polint(rg(1:4),     f2d(1,1:4,j),     4, r, f2d(1,i-3,j),  eps)
          r =  rg(nr) + dr * i
          CALL polint(rg(nr-3:nr), f2d(1,nr-3:nr,j), 4, r, f2d(1,nr+i,j), eps)
        END DO
      END DO

      DO i=-2,nr+3
        DO j=1,3
          z =  zg(1)  - dz * (4 - j)
          CALL polint(zg(1:4),     f2d(1,i,1:4),     4, z, f2d(1,i,j-3),  eps)
          z =  zg(nz) + dz * j
          CALL polint(zg(nz-3:nz), f2d(1,i,nz-3:nz), 4, z, f2d(1,i,nz+j), eps)
        END DO
      END DO

      nx2d =  nr
      ny2d =  nz

      CALL splin2(rminb, rmaxb, zminb, zmaxb)

      dr = (rmaxb - rminb) / (nr0b - 1)
      dz = (zmaxb - zminb) / (nz0b - 1)

#ifdef MPI
      IF(myrank == 0)THEN
#endif
      DO j=1,nz0b
        DO i=1,nr0b
          x2d(1) =  rminb + dr * (i - 1)
          x2d(2) =  zminb + dz * (j - 1)
          CALL spl2df(x2d, w2d)
          f3d(5,i,j,k) =  w2d(1)
          WRITE(72,'(10ES15.7)') x2d(1), x2d(2), w2d(1)
        END DO
        WRITE(72,*)
      END DO
#ifdef MPI
      END IF
#endif

      DO k=1,nt0b
        DO j=1,nz0b
          DO i=1,nr0b
            x2d(1) =  rminb + dr * (i - 1)
            x2d(2) =  zminb + dz * (j - 1)
            CALL spl2df(x2d, w2d)
            f3d(5,i,j,k) =  w2d(1)
          END DO
        END DO
      END DO

      DEALLOCATE(rg, zg, f2d)
 
    ELSE

      DO k=1,nt0b
        f3d(5,1:nr0b,1:nz0b,k) =  psig(1:nr0b,1:nz0b)
      END DO

    END IF

    DEALLOCATE(psig)


  END SUBROUTINE read_eqdsk

  SUBROUTINE read_eqdata

    IMPLICIT NONE
    INTEGER :: irdm,      &
      &        izdm2,     &
      &        irzdm2,    &
      &        ivdm,      &
      &        i,         &
      &        j,         &
      &        k,         &
      &        itemp(10)
    REAL(DP) :: saxis, &
      &         temp
    REAL(DP), ALLOCATABLE :: rr(:),    &
      &                      zz(:),    &
      &                      psi(:),   &
      &                      psig(:,:)
    CHARACTER(LEN=10) :: header(6)
#ifdef MPI
!For MPI
    INTEGER :: icode, &
      &        ierr
#endif


#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE READ_EQDATA                                 '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *
#ifdef MPI
    END IF
#endif

    READ(26,*) itemp(1), temp
    READ(26,*) irdm, izdm2, irzdm2, itemp(1), itemp(2), ivdm, itemp(3), itemp(4), itemp(5)

    ALLOCATE(rr(irdm), zz(izdm2), psi(irzdm2))

    READ(26,*) (psi(i), temp, i=1,irzdm2)
    READ(26,*) (rr(i), i=1,irdm)
    READ(26,*) (temp, i=1,irdm)
    READ(26,*) (zz(i), i=1,izdm2)

    READ(26,*) (temp, i=1,ivdm)
    READ(26,*) (temp, temp, temp, i=1,ivdm)

    READ(26,*) temp, temp, temp, saxis, temp, temp, temp, temp, temp, temp, temp, temp, temp, temp, temp

    ALLOCATE(psig(irdm,izdm2))

    k =  0
    DO j=1,izdm2
      DO i=1,irdm
        k         =  k + 1
        psig(i,j) =  psi(k)
      END DO
    END DO

    IF(irdm /= nr0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdata: parameter error  irdm      = ', irdm,     '  nr0b  = ', nr0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(izdm2 /= nz0b)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdata: parameter error  izdm2     = ', izdm2,    '  nz0b  = ', nz0b
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(rr(1) /= rminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdata: parameter error  rr(1)     =', rr(1),     '  rminb = ', rminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(rr(irdm) /= rmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdata: parameter error  rr(irdm)  =', rr(irdm),  '  rmaxb = ', rmaxb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(zz(1) /= zminb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdata: parameter error  zz(1)     =', zz(1),     '  zminb = ', zminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

    IF(zz(izdm2) /= zmaxb)THEN
#ifdef MPI
      IF(myrank == 0) &
#endif
      PRINT *, 'read_eqdata: parameter error  zz(izdm2) =', zz(izdm2), '  zminb = ', zminb
#ifdef MPI
      CALL MPI_ABORT(MPI_COMM_WORLD, icode, ierr)
#else
      STOP
#endif
    END IF

#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *, ' Reading normalized poloidal flux as EQDATA format'
    PRINT *
#ifdef MPI
    END IF
#endif

    psig(:,:) =  psig(:,:) - saxis
    psig(:,:) =  psig(:,:) / ABS(saxis)

    DO k=1,nt0b
      f3d(5,1:nr0b,1:nz0b,k) =  psig(1:nr0b,1:nz0b)
    END DO

    DEALLOCATE(rr, zz, psi, psig)


  END SUBROUTINE read_eqdata

  SUBROUTINE free_mem_field

    IMPLICIT NONE


    DEALLOCATE(f3d)


  END SUBROUTINE free_mem_field

  SUBROUTINE magset

    IMPLICIT NONE

    INTEGER :: i,  &
      &        j,  &
      &        jq, &
      &        k,  &
      &        l
    REAL(DP) :: r,  &
      &         z,  &
      &         bb, &
      &         eps
    CHARACTER(LEN=100) :: fmt


#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE MAGSET                                      '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *
#ifdef MPI
    END IF
#endif

    CALL read_field

    ALLOCATE(rg(nr0b), zg(nz0b))

    pi2m  =  pi2 / mtor

    nx3d =  nr0b
    ny3d =  nz0b
    nz3d =  nt0b + 1

    delrb = (rmaxb - rminb) / (nr0b - 1)
    delzb = (zmaxb - zminb) / (nz0b - 1)

    DO i=1,nr0b
      rg(i) =  rminb + delrb * (i - 1)
    END DO

    DO j=1,nz0b
      zg(j) =  zminb + delzb * (j - 1)
    END DO
 
    CALL splin3(rminb, rmaxb, zminb, zmaxb, 0.0_DP, pi2m)

    f3d(4,:,:,:) =  SQRT(f3d(1,:,:,:)**2 + f3d(2,:,:,:)**2 + f3d(3,:,:,:)**2)

    bmax =  5.0_DP

    DO k=1,nt0b
      DO j=1,nz0b
        DO i=1,nr0b
          bb =  f3d(4,i,j,k)
          IF(bb > bmax)THEN
            f3d(1,i,j,k) =  bmax * f3d(1,i,j,k) / bb
            f3d(2,i,j,k) =  bmax * f3d(2,i,j,k) / bb
            f3d(3,i,j,k) =  bmax * f3d(3,i,j,k) / bb
          END IF
        END DO
      END DO
    END DO

    DO k=1,nt0b
      DO j=1,nz0b
        DO i=1,3
          DO l=1,l3d
            r =  rminb - delrb * (4 - i)
            CALL polint(rg(1:4),         f3d(l,1:4,j,k),         4, r, f3d(l,i-3,j,k),    eps)
            r =  rmaxb + delrb * i
            CALL polint(rg(nr0b-3:nr0b), f3d(l,nr0b-3:nr0b,j,k), 4, r, f3d(l,nr0b+i,j,k), eps)
          END DO
        END DO
      END DO
    END DO

    DO k=1,nt0b
      DO i=-2,nr0b+3
        DO j=1,3
          DO l=1,l3d
            z =  zminb - delzb * (4 - j)
            CALL polint(zg(1:4),         f3d(l,i,1:4,k),         4, z, f3d(l,i,j-3,k),    eps)
            z =  zmaxb + delzb * j
            CALL polint(zg(nz0b-3:nz0b), f3d(l,i,nz0b-3:nz0b,k), 4, z, f3d(l,i,nz0b+j,k), eps)
          END DO
        END DO
      END DO
    END DO

    IF(nt0b == 1)THEN
      f3d(:,:,:,nt0b+1) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+2) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+3) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+4) =  f3d(:,:,:,1)

      f3d(:,:,:,-2)     =  f3d(:,:,:,1)
      f3d(:,:,:,-1)     =  f3d(:,:,:,1)
      f3d(:,:,:,0)      =  f3d(:,:,:,1)
    ELSE
      IF(lsymmetry)THEN

        IF(mod(nt0b,2) /= 0)THEN
#ifdef MPI
          IF(myrank == 0)THEN
#endif
          PRINT *
          PRINT *, ' parameter nt0b is NOT even number!!!'
          PRINT *
#ifdef MPI
          END IF
#endif
          STOP
        END IF

        f3d(:,:,:,1)        =  0.5_DP * (f3d(:,:,:,2)      + f3d(:,:,:,nt0b))
        f3d(:,:,:,nt0b/2+1) =  0.5_DP * (f3d(:,:,:,nt0b/2) + f3d(:,:,:,nt0b/2+2))
      END IF 

      f3d(:,:,:,nt0b+1) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+2) =  f3d(:,:,:,2)
      f3d(:,:,:,nt0b+3) =  f3d(:,:,:,3)
      f3d(:,:,:,nt0b+4) =  f3d(:,:,:,4)

      f3d(:,:,:,-2)     =  f3d(:,:,:,nt0b-2)
      f3d(:,:,:,-1)     =  f3d(:,:,:,nt0b-1)
      f3d(:,:,:,0)      =  f3d(:,:,:,nt0b)
    END IF

#ifdef MPI
    IF(myrank == 0)THEN
#endif
    PRINT *
    PRINT *, ' COORDINATES SYSTEMS '
    PRINT *, '------------------------------------------'
    PRINT *, '   M    Rmin     Rmax     Zmin     Zmax'

    fmt = '(I5, 4F9.4)'

    PRINT fmt, mtor, rminb, rmaxb, zminb, zmaxb

    PRINT *
    PRINT *, ' RESOLUSIONS'
    PRINT *, '----------------------------------------------------'
    PRINT *, '   nr0b    nt0b    nz0b   delrb     dphi    delzb'

    fmt = '(3I8, 4F9.4)'
    PRINT fmt, nr0b, nt0b, nz0b, delrb, pi2m / nt0b, delzb
#ifdef MPI
    END IF
#endif

    DEALLOCATE(rg, zg)


  END SUBROUTINE magset

  SUBROUTINE mgval1 (r, phi, z,     & ! (in)
    &                br, bp, bz, bb & ! (out)
    &               )

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: r,   & ! major radius in computational region [m]
      &                      phi, & ! toroidal angle in computational region [rad]
      &                      z      ! height in computational region [m]
    REAL(DP), INTENT(OUT) :: br, & ! BR component on (R,phi,Z
      &                      bp, & ! Bphi component on (R,phi,Z)
      &                      bz, & ! BZ component on (R,phi,Z)
      &                      bb    ! B strength on (R,phi,Z) [T]
    INTEGER :: iphi
    REAL(DP) :: phi1,   &
      &         xd(3),  &
      &         w0(l3d)


    iphi =  phi / pi2m
    phi1 =  phi - pi2m * iphi
    IF(phi1 < 0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)  phi1 =  phi1 - pi2m       

    xd(1) =  r
    xd(2) =  z
    xd(3) =  phi1

    CALL spl3df(xd(:), w0(:))

    br =  w0(1)
    bp =  w0(2)
    bz =  w0(3)
    bb =  w0(4)


    RETURN
  END SUBROUTINE mgval1

  SUBROUTINE mgval2 (r, phi, z,          & ! (in)
    &                b, dbdr, dbdp, dbdz & ! (out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: r,         & ! major radius in computational region [m]
      &                      phi,       & ! toroidal angle in computational region [rad]
      &                      z            ! height in computational region [m]
    REAL(DP), INTENT(OUT) :: b(l3d),    & ! B vector (1:BR,2:BT,3:BZ,4:B)
      &                      dbdr(l3d), & ! R-derivative of B vector (1:BR,2:BZ,3:BT,4:B)
      &                      dbdp(l3d), & ! phi-derivative of B vector (1:BR,2:BZ,3:BT,4:B)
      &                      dbdz(l3d)    ! Z-derivative of B vector (1:BR,2:BZ,3:BT,4:B)
!Local variables
    INTEGER :: iphi
    REAL(DP) :: phi1,    &
      &         xd(3),   &
      &         w0(l3d), &
      &         wx(l3d), &
      &         wy(l3d), &
      &         wz(l3d)


    iphi =  phi / pi2m
    phi1 =  phi - pi2m * iphi
    IF(phi1 < 0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)  phi1 =  phi1 - pi2m      

    xd(1) =  r
    xd(2) =  z
    xd(3) =  phi1

    CALL spl3dd(xd(:), w0(:), wx(:), wy(:), wz(:))

    b(:)    =  w0(:)
    dbdr(:) =  wx(:)
    dbdz(:) =  wy(:)
    dbdp(:) =  wz(:)


    RETURN
  END SUBROUTINE mgval2

  SUBROUTINE mgval3 (r, phi, z, & ! (in)
    &                s          & ! (out)
    &               )

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: r,   & ! major radius in computational region [m]
      &                      phi, & ! toroidal angle in computational region [rad]
      &                      z      ! height in computational region [m]
    REAL(DP), INTENT(OUT) :: s     ! flux label

    INTEGER :: iphi
    REAL(DP) :: phi1,  &
      &         xd(3),  &
      &         w0(l3d)


    iphi =  phi / pi2m
    phi1 =  phi - pi2m * iphi
    IF(phi1 < 0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)  phi1 =  phi1 - pi2m

    xd(1) =  r
    xd(2) =  z
    xd(3) =  phi1

    CALL spl3df(xd(:), w0(:))

    s =  w0(5)


    RETURN
  END SUBROUTINE mgval3

  SUBROUTINE polint (xa, ya, n, x, & ! (in)
    &                y, dy         & ! (out)
    &               )

    IMPLICIT NONE

!Arguments
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x,     &
      &                     xa(n), &
      &                     ya(n)
    REAL(DP), INTENT(OUT) :: y,  &
    &                        dy
!Local variables
    INTEGER :: m,  &
      &        ns, &
      &        i
    REAL(DP) :: den,  &
      &         dif,  &
      &         dift, &
      &         ho,   &
      &         hp,   &
      &         w,    &
      &         c(n), &
      &         d(n)
  

    ns  =  1
    dif =  ABS(x - xa(1))

    loop100 : DO i=1,n
      dift =  ABS(x - xa(i))
      IF(dift < dif)THEN
        ns  =  i
        dif =  dift
      END IF
      c(i) =  ya(i)
      d(i) =  ya(i)
    END DO loop100

    y  =  ya(ns)
    ns =  ns - 1

    loop200 : DO m=1,n-1
      loop210 : DO i=1,n-m
        ho  =  xa(i)   - x
        hp  =  xa(i+m) - x
        w   =  c(i+1)  - d(i)
        den =  ho      - hp
        IF(den == 0.0_DP) STOP 'failure in polint'
        den  =  w  / den
        d(i) =  hp * den
        c(i) =  ho * den
      END DO loop210
      IF(2 * ns < n - m)THEN
        dy =  c(ns+1)
      ELSE
        dy =  d(ns)
        ns =  ns - 1
      END IF
      y =  y + dy
    END DO loop200


    RETURN
  END SUBROUTINE polint

#ifdef NETCDF
  SUBROUTINE check(status)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: status


    IF(STATUS /= NF90_NOERR)THEN
      PRINT *, TRIM(NF90_STRERROR(status))
      STOP "Stopped"
    END IF


  END SUBROUTINE check
#endif


END MODULE cylindrical_coord_mod
