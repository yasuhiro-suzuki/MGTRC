!=cal_avr2d.f90
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
SUBROUTINE cal_flxqnt

  USE kind_spec
  USE param1,                ONLY : pi,                &
    &                               pi2,               &
    &                               c12,               &
    &                               c13,               &
    &                               c14
  USE fline_mod,             ONLY : def_start,         &
    &                               lfigout,           &
    &                               ltext,             &
    &                               mr,                &
    &                               mcirc,             &
    &                               nstep,             &
    &                               rstart,            &
    &                               zstart,            &
    &                               pstart,            &
    &                               sstart,            &
    &                               rout,              &
    &                               zout,              &
    &                               sout,              &
    &                               drflx,             &
    &                               dzflx
  USE axis_mod,              ONLY : raxis,             &
    &                               zaxis
  USE cylindrical_coord_mod, ONLY : pi2m,              &
    &                               mtor,              &
    &                               rmaxb,             &
    &                               mgval1,            &
    &                               mgval3
  USE ode_mod,               ONLY : odeint,            &
    &                               odeadm5,           &
    &                               dopri5
  USE plot_mod,              ONLY : init_plot_profile, &
    &                               end_plot_profile,  &
    &                               plot_profile_reff

  IMPLICIT NONE

  INTEGER :: ntheta,     &
    &        ntheta1,    &
    &        m_turn,     &
    &        m_turn_old, &
    &        i_index,    &
    &        i,          &
    &        j,          &
    &        js,         &
    &        k,          &
    &        n
  REAL(DP), PARAMETER :: eps =  1.0E-03_DP
  REAL(DP) :: s,       &
    &         ds,      &
    &         s1,      &
    &         s2,      &
    &         tf1,     &
    &         tf2,     &
    &         x,       &
    &         y,       &
    &         z,       &
    &         r,       &
    &         phi,     &
    &         r0,      &
    &         z0,      &
    &         r1,      &
    &         z1,      &
    &         br,      &
    &         bp,      &
    &         bz,      &
    &         bb,      &
    &         dl,      &
    &         theta,   &
    &         theta0,  &
    &         theta1,  &
    &         dtheta,  &
    &         dtheta2, &
    &         dradphi, &
    &         dzadphi, &
    &         drdphi,  &
    &         dzdphi,  &
    &         th0,     &
    &         th1,     &
    &         th2,     &
    &         th01,    &
    &         th02,    &
    &         th12,    &
    &         ath,     &
    &         bth,     &
    &         cth,     &
    &         rmin1,   &
    &         rmax1,   &
    &         zmin1,   &
    &         zmax1,   &
    &         rmin2,   &
    &         rmax2,   &
    &         zmin2,   &
    &         zmax2,   &
    &         elng1,   &
    &         elng2,   &
    &         x1,      &
    &         y1,      &
    &         x2,      &
    &         y2,      &
    &         x3,      &
    &         y3,      &
    &         x4,      &
    &         y4,      &
    &         b1,      &
    &         b2,      &
    &         b3,      &
    &         b4
  REAL(DP), ALLOCATABLE :: rav1(:),       &
    &                      rav2(:),       &
    &                      phit(:),       &
    &                      psip(:),       &
    &                      iota1(:),      &
    &                      iota2(:),      &
    &                      vp(:),         &
    &                      volume(:),     &
    &                      lc(:),         &
    &                      elong(:),      &
    &                      axisr(:),      &
    &                      axisz(:),      &
    &                      wrk1(:),       &
    &                      rk(:),         &
    &                      thetak(:),     &
    &                      ri(:),         &
    &                      thetai(:),     &
    &                      rsurfk(:,:),   &
    &                      zsurfk(:,:),   &
    &                      rsurfi(:,:,:), &
    &                      zsurfi(:,:,:)
  CHARACTER(LEN=100) :: fmt
!
  REAL(DP) :: rmin, &
    &         rmax, &
    &         rmid
  REAL(DP), ALLOCATABLE :: ss(:)
!for odeint
  INTEGER, PARAMETER :: ln =  2
  INTEGER :: lnx
  INTEGER, ALLOCATABLE :: icount(:), &
    &                     iout(:)
  REAL(DP) :: h, &
    &         x0
  REAL(DP), ALLOCATABLE :: f0(:,:), &
    &                      f(:,:,:)
  EXTERNAL :: subf2d
!for check vessel
  INTEGER :: ivessel


  PRINT *
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *, '          SUBROUTINE cal_flxqnt                                  '
  PRINT *, ' ----------------------------------------------------------------'
  PRINT *

  x0      =  pi2m * pstart
  h       =  pi2m / nstep
  ntheta  =  mtor * mcirc
  lnx     =  nstep * ntheta

  ntheta1 =  361

  fmt = '(A11, I12)'

  PRINT fmt, ' nstep   = ', nstep
  PRINT fmt, ' ntheta  = ', ntheta
  PRINT fmt, ' ntheta1 = ', ntheta1

  fmt = '(A11, ES12.4)'

  PRINT fmt, ' h       = ', h
  PRINT *
  PRINT *

  ALLOCATE(f0(0:mr,ln), f(0:1,ln,lnx), icount(0:mr), iout(0:mr))
  ALLOCATE(rav1(0:mr), rav2(0:mr), ss(0:mr), phit(0:mr), psip(0:mr), iota1(0:mr), iota2(0:mr), vp(0:mr), volume(0:mr), lc(0:mr), elong(0:mr))
  ALLOCATE(axisr(nstep), axisz(nstep), wrk1(lnx), rk(-1:ntheta+2), thetak(-1:ntheta+2), ri(ntheta1), thetai(ntheta1), rsurfk(nstep,ntheta), zsurfk(nstep,ntheta), rsurfi(mr,nstep,ntheta), zsurfi(mr,nstep,ntheta))

  icount(:) =  0
  iout(:)   =  0

  rav1(:)   =  0.0_DP
  rav2(:)   =  0.0_DP
  ss(:)     =  0.0_DP
  phit(:)   =  0.0_DP
  psip(:)   =  0.0_DP
  iota1(:)  =  0.0_DP
  iota2(:)  =  0.0_DP
  vp(:)     =  0.0_DP
  volume(:) =  0.0_DP
  lc(:)     =  0.0_DP
  elong(:)  =  0.0_DP

  dtheta =  pi2 / (ntheta1 - 1)
  DO i=1,ntheta1
    thetai(i) =  dtheta * (i - 1)
  END DO

  SELECT CASE(TRIM(def_start))

    CASE('rtheta')

      IF(sstart > sout)THEN
        sstart =  0.0_DP
        sout   =  1.0_DP
      END IF

      IF(sstart < 0) sstart =  0.0_DP
      IF(sout   > 1) sout   =  1.0_DP

      ds = (sout - sstart) / mr
      DO js=0,mr
        s      =  sstart + ds * js
        ss(js) =  s**2
      END DO

      IF(ss(mr) == 1.0_DP) ss(mr) =  ss(mr) - eps

      f0(0,1) =  raxis
      f0(0,2) =  zaxis

      DO js=1,mr

        rmin   =  raxis
        rmax   =  rmaxb

        DO i=1,100
          rmid = (rmin + rmax) / 2

          CALL mgval3(rmid, x0, zaxis, s)

          IF(s < ss(js))THEN
            rmin =  rmid
          ELSE
            rmax =  rmid
          END IF
        END DO

        f0(js,1) =  rmid
        f0(js,2) =  zaxis

        !IF(js == mr) f0(js,1) =  f0(js,1) - eps

      END DO

    CASE('file', 'FILE')

      STOP 'DEF_START:FILE mode does not work in CAL_FLXQNT'

    CASE DEFAULT

      f0(0,1) =  raxis
      f0(0,2) =  zaxis

      IF(rstart == 0.0_DP) rstart =  raxis
      IF(zstart == 0.0_DP) zstart =  zaxis

      IF(rout /= 0.0_DP) drflx = (rout - rstart) / mr
      IF(zout /= 0.0_DP) drflx = (zout - zstart) / mr

      DO js=1,mr
        f0(js,1) =  rstart + drflx * js
        f0(js,2) =  zstart + dzflx * js
      END DO

  END SELECT

  !CALL odeint(x0, f0(0,:), ln, h, lnx, dopri5, subf2d, f(0,:,:), icount(0), iout(0))
  CALL odeadm5(x0, f0(0,:), ln, h, lnx, dopri5, subf2d, f(0,:,:), icount(0), iout(0))

  IF(iout(0) == 1) STOP ' iout = 1!!!!'

  DO i=1,nstep
    axisr(i) =  f(0,1,i)
    axisz(i) =  f(0,2,i)
  END DO

  PRINT *
  PRINT *, ' flux surface quantities '
  PRINT *

  fmt = '(A)'

  PRINT fmt, '   mr  iout   icount     R        Z       <r>(2)      iota(1)      iota(2)         Vp       elongation'
  PRINT fmt, '                                        (int:rdl/B)  (dth/dphi)   (Rbp/rbt)    (int:dl/B)'

  fmt = '(2I5, I10, 2F9.5, 20ES13.5)'

  DO js=1,mr

    !CALL odeint(x0, f0(js,:), ln, h, lnx, dopri5, subf2d, f(1,:,:), icount(js), iout(js))
    CALL odeadm5(x0, f0(js,:), ln, h, lnx, dopri5, subf2d, f(1,:,:), icount(js), iout(js))

    DO i=1,icount(js)
      phi      =  x0 + h * (i - 1)
      r1       =  f(1,1,i)
      z1       =  f(1,2,i)
      r0       =  f(0,1,i)
      z0       =  f(0,2,i)
      r        =  r1 - r0
      z        =  z1 - z0
      wrk1(i)  =  ATAN2(z, r)
      CALL mgval1(r1, phi, z1, br, bp, bz, bb)
      dl       =  bb * r1 * h / bp
      vp(js)   =  vp(js) + ABS(dl) / bb
      lc(js)   =  lc(js) + ABS(dl)
      rav2(js) =  rav2(js) + SQRT(r**2 + z**2) * ABS(dl) / bb
    END DO

    rav2(js) =  rav2(js) / vp(js)
    vp(js)   =  vp(js) / mcirc

    theta      =  0.0_DP
    m_turn_old =  0
    DO i=2,icount(js)
      dtheta =  wrk1(i) - wrk1(i-1)
      IF(dtheta < 0.0_DP)THEN
        dtheta2 =  dtheta + pi2
        IF(ABS(dtheta2) < ABS(dtheta)) dtheta =  dtheta2
      ELSE
        dtheta2 =  dtheta - pi2
        IF(ABS(dtheta2) < ABS(dtheta)) dtheta =  dtheta2
      END IF
      theta =  theta + dtheta
      IF(theta > 0.0_DP)THEN
        m_turn =  theta / pi2
        IF((m_turn /= 0) .AND. (m_turn /= m_turn_old))THEN
          m_turn_old =  m_turn
          i_index    =  i
        END IF
      ELSE
        m_turn =  theta / pi2
        IF((m_turn /= 0) .AND. (m_turn /= m_turn_old))THEN
          m_turn_old =  m_turn
          i_index    =  i
        END IF
      END IF
    END DO

    theta     =  pi2 * m_turn_old
    phi       =  h * i_index
    iota1(js) =  theta / phi

    theta =  0.0_DP
    DO i=2,icount(js)
      phi     =  x0 + h * (i - 1)
      r1      =  f(1,1,i)
      z1      =  f(1,2,i)
      r0      =  f(0,1,i)
      z0      =  f(0,2,i)
      r       =  r1 - r0
      z       =  z1 - z0
      CALL mgval1(r0, phi, z0, br, bp, bz, bb)
      dradphi = r0 * br / bp
      dzadphi = r0 * bz / bp
      CALL mgval1(r1, phi, z1, br, bp, bz, bb)
      drdphi = r1 * br / bp
      dzdphi = r1 * bz / bp
      dtheta = r / (r**2 + z**2) * (dzdphi - dzadphi) - z / (r**2 + z**2) * (drdphi - dradphi)
      theta =  theta + dtheta * h
    END DO

    phi       =  h * icount(js)
    iota2(js) =  theta / phi

    IF(iout(js) == 0)THEN

      DO i=1,ntheta
        DO n=1,nstep
          k           = (i - 1) * nstep + n
          rsurfk(n,i) =  f(1,1,k)
          zsurfk(n,i) =  f(1,2,k)
        END DO
      END DO

      DO n=1,nstep

        CALL order(ntheta, rsurfk(n,:), zsurfk(n,:), axisr(n), axisz(n))

        DO i=1,ntheta
          r     =  rsurfk(n,i) - axisr(n)
          z     =  zsurfk(n,i) - axisz(n)
          rk(i) =  SQRT(r**2 + z**2)
          theta =  ATAN2(z, r)
          IF(theta < 0.0_DP) theta =  pi2   + theta
          IF(theta > pi2)    theta =  theta - pi2
          thetak(i) =  theta
        END DO

        rk(-1)           =  rk(ntheta-1)
        rk(0)            =  rk(ntheta)
        rk(ntheta+1)     =  rk(1)
        rk(ntheta+2)     =  rk(2)
        thetak(-1)       =  thetak(ntheta-1) - pi2
        thetak(0)        =  thetak(ntheta) - pi2
        thetak(ntheta+1) =  thetak(1) + pi2
        thetak(ntheta+2) =  thetak(2) + pi2

        DO i=1,ntheta1
          loop_check : DO j=-1,ntheta+2
            IF(thetak(j) > thetai(i))THEN
              dtheta =  thetak(j) - thetak(j-1)
              IF(dtheta <= 0.2_DP)THEN
                ri(i) =  rk(j-1) + (rk(j) - rk(j-1)) * (thetai(i) - thetak(j-1)) / (thetak(j) - thetak(j-1))
              ELSE
                th0   =  thetai(i)   - thetak(j-1)
                th1   =  thetai(i)   - thetak(j)
                th2   =  thetai(i)   - thetak(j+1)
                th01  =  thetak(j-1) - thetak(j)
                th02  =  thetak(j-1) - thetak(j+1)
                th12  =  thetak(j)   - thetak(j+1)
                ath   = (th1 * th2) / (th01  * th02)
                bth   = (th0 * th2) / (-th01 * th12)
                cth   = (th0 * th1) / (th02  * th12)
                ri(i) =  ath * rk(j-1) + bth * rk(j) + cth * rk(j+1)
              END IF
              EXIT loop_check
            END IF
          END DO loop_check
        END DO

        ri(ntheta) = ri(1)

        DO i=1,ntheta1
          rsurfi(js,n,i) =  ri(i) * COS(thetai(i)) + axisr(n)
          zsurfi(js,n,i) =  ri(i) * SIN(thetai(i)) + axisz(n)
        END DO

      END DO

      rmin1     =  MINVAL(rsurfi(js,1,:))
      rmax1     =  MAXVAL(rsurfi(js,1,:))
      zmin1     =  MINVAL(zsurfi(js,1,:))
      zmax1     =  MAXVAL(zsurfi(js,1,:))

      rmin2     =  MINVAL(rsurfi(js,nstep/2+1,:))
      rmax2     =  MAXVAL(rsurfi(js,nstep/2+1,:))
      zmin2     =  MINVAL(zsurfi(js,nstep/2+1,:))
      zmax2     =  MAXVAL(zsurfi(js,nstep/2+1,:))

      elng1     = (zmax1 - zmin1) / (rmax1 - rmin1)
      elng2     = (zmax2 - zmin2) / (rmax2 - rmin2)

      elong(js) = (elng1 + elng2) / 2

    END IF

    PRINT fmt, js, iout(js), icount(js), f0(js,1), f0(js,2), rav2(js), iota1(js), iota2(js), vp(js), elong(js)

  END DO

  IF(iout(1) == 0)THEN
    s2  =  0.0_DP
    tf2 =  0.0_DP
    DO n=1,nstep
      phi =  x0 + h * (n - 1)
      x1  =  axisr(n)
      y1  =  axisz(n)
      CALL mgval1(x1, phi, y1, br, b1, bz, bb)
      s1  =  0.0_DP
      tf1 =  0.0_DP
      DO i=2,ntheta1
        x2 =  rsurfi(1,n,i-1)
        y2 =  zsurfi(1,n,i-1)
        x3 =  rsurfi(1,n,i)
        y3 =  zsurfi(1,n,i)
        ds =  c12 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
        CALL mgval1(x2, phi, y2, br, b2, bz, bb)
        CALL mgval1(x3, phi, y3, br, b3, bz, bb)
        s1  =  s1  + ds
        tf1 =  tf1 + ds * c13 * (b1 + b2 + b3)
      END DO
      s2  =  s2  + s1
      tf2 =  tf2 + tf1
    END DO
    ss(1)   =  s2  / nstep
    phit(1) =  tf2 / nstep
  END IF

  IF(mr >= 2)THEN 
    DO js=2,mr
      IF(iout(js) == 0)THEN
        s2  =  0.0_DP
        tf2 =  0.0_DP
        DO n=1,nstep
          phi =  x0 + h * (n - 1)
          s1  =  0.0_DP
          tf1 =  0.0_DP
          DO i=2,ntheta1
            x1 =  rsurfi(js-1,n,i-1)
            y1 =  zsurfi(js-1,n,i-1)
            x2 =  rsurfi(js,  n,i-1)
            y2 =  zsurfi(js,  n,i-1)
            x3 =  rsurfi(js,  n,i)
            y3 =  zsurfi(js,  n,i)
            x4 =  rsurfi(js-1,n,i)
            y4 =  zsurfi(js-1,n,i)
            ds =  c12 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) + (x3 - x1) * (y4 - y1) - (x4 - x1) * (y3 - y1))
            CALL mgval1(x1, phi, y1, br, b1, bz, bb)
            CALL mgval1(x2, phi, y2, br, b2, bz, bb)
            CALL mgval1(x3, phi, y3, br, b3, bz, bb)
            CALL mgval1(x4, phi, y4, br, b4, bz, bb)
            s1  =  s1  + ds
            tf1 =  tf1 + ds * c14 * (b1 + b2 + b3 + b4)
          END DO
          s2  =  s2  + s1
          tf2 =  tf2 + tf1
        END DO
        ss(js)   =  ss(js-1)   + s2  / nstep
        phit(js) =  phit(js-1) + tf2 / nstep
      END IF
    END DO
  END IF

  rav1(:)  =  SQRT(ABS(ss(:)) / pi)

  iota1(0) =  3 * iota1(1) - 3 * iota1(2) + iota1(3)
  iota2(0) =  3 * iota2(1) - 3 * iota2(2) + iota2(3)
  vp(0)    =  3 * vp(1)    - 3 * vp(2)    + vp(3)
  lc(0)    =  3 * lc(1)    - 3 * lc(2)    + lc(3)

  DO js=2,mr
    psip(js)   =  psip(js-1)   + 0.5_DP * (iota1(js-1) + iota1(js)) * (phit(js) - phit(js-1))
    volume(js) =  volume(js-1) + 0.5_DP * (vp(js-1)    + vp(js)) * ABS(phit(js) - phit(js-1))
  END DO

  DO js=2,mr-1
    psip(js)   = (psip(js+1) + psip(js)) / 2
    volume(js) = (volume(js+1) + volume(js)) / 2
  END DO

  psip(0)    =  0.0_DP
  volume(0)  =  0.0_DP

  psip(mr)   =  3 * psip(mr-1)   - 3 * psip(mr-2)   + psip(mr-3)
  volume(mr) =  3 * volume(mr-1) - 3 * volume(mr-2) + volume(mr-3)

  DO js=1,mr
    IF(iout(js) == 1)THEN
      psip(js)   =  0.0_DP
      volume(js) =  0.0_DP
    END IF
  END DO
 
  PRINT *
  PRINT *, ' flux surface quantities '
  PRINT *

  fmt = '(A)'

  PRINT fmt, '   mr  iout     <s>        <r>(1)       <r>(2)         Phi          Psi       iota(1)      iota(2)         Vp         Volume        well          Lc'
  PRINT fmt, '                         (from <s>)   (int:rdl/B)  (from <s>)   (int:iota1)  (dth/dphi)   (Rbp/rbt)    (int:dl/B)    (int:Vp)    (1-vp/vp0)    (int:dl)'

  fmt = '(2I5, 20ES13.5)'

  loop0060 : DO js=0,mr
    PRINT fmt, js, iout(js), ss(js), rav1(js), rav2(js), phit(js), psip(js), iota1(js), iota2(js), vp(js), volume(js), (vp(0)-vp(js))/vp(0), lc(js)
  END DO loop0060

  IF(ltext)THEN
    DO js=0,mr
      WRITE(60,'(20(ES15.7,A1))') f0(js,1), ',', f0(js,2), ',', rav1(js), ',', phit(js), ',', iota1(js), ',', iota2(js), ',', (vp(0) - vp(js)) / vp(0), ',', lc(js), ','
    END DO
  END IF

  IF(lfigout)THEN

    CALL init_plot_profile(mr)
    CALL plot_profile_reff(mr, rav2, rav2, iota1, vp, volume, lc)
    CALL end_plot_profile

  END IF

  DEALLOCATE(axisr, axisz, wrk1, rsurfk, zsurfk, rsurfi, zsurfi)
  DEALLOCATE(rav1, rav2, ss, phit, psip, iota1, iota2, vp, volume, lc)
  DEALLOCATE(f0, f)


END SUBROUTINE cal_flxqnt
