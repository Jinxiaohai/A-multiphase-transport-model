!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!....................amptsub.f
!.....this file contains 4 sections:
!.....1. ART subroutines;
!.....2. ART functions;
!.....3. ART block data;
!.....4. subprocesses borrowed from other codes.
!.....5. the previous artana.f
!.....6. the previous zpcsub.f
!.....7. subroutine getnp
!.....Note that Parts1-4 are the previous artsub.f
!
!=======================================================================
!.....subroutine to set up ART parameters and analysis files
!.....before looping different events
Subroutine artset
!
  Parameter (amu=0.9383, nxymax=10001)
  Double Precision dpcoal, drcoal, ecritl
  Integer zta, zpr
  Common /gg/dx, dy, dz, dpx, dpy, dpz
!lin-10/03/03
!     "SAVE   " (without argument) is used for most subroutines and functions,
!     this is important for the success when using "f77" to compile:
!c      SAVE /gg/
  Common /zz/zta, zpr
!c      SAVE /zz/
  Common /run/num
!c      SAVE /RUN/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
!c      SAVE /INPUT2/
  Common /input3/plab, elab, zeropt, b0, bi, bm, dencut, cycbox
!c      SAVE /INPUT3/
  Common /imulst/iperts
!c      SAVE /imulst/
  Common /coal/dpcoal, drcoal, ecritl
  Common /anim/nevent, isoft, isflag, izpc
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  Common /xyembed/nxyjet, xyjet(nxymax, 2)
  Save
!lin-10/03/03  ecritl: local energy density below which a parton
!     will freeze out (in GeV/fm^3), for improvements on string melting,
!     not used in this version of AMPT:
!lin-4/2008
!      data ecritl/1.d0/
  ecritl = 1.D0
!
!     combine ART initialization into ampt.ini:
!     (Note that the following values are relics from the old ART structure)
!.....input parameter file
!      OPEN(13, FILE = 'art1.ini', STATUS = 'UNKNOWN')
!      READ (13, *) MASSTA, ZTA
  massta = 1
  zta = 1
!      write(12,*) massta, zta, ' massta, zta'
!      READ (13, *) MASSPR, ZPR
  masspr = 1
  zpr = 1
!      write(12,*) masspr, zpr, ' masspr, zpr'
!      READ (13, *) PLAB, IPLAB
  plab = 14.6
  iplab = 2
!      write(12,*) plab, iplab, ' plab, iplab'
  If (iplab==2) Then
    elab = sqrt(plab**2+amu**2) - amu
  Else
    elab = plab
  End If
  elab = elab*1000.
!      READ (13, *) ZEROPT
  zeropt = 0.
!      write(12,*) zeropt, ' zeropt'
!lin-10/03/03 ISEED was used as a seed for random number inside ART,
!     not used in AMPT:
  iseed = 700721
!     0/1: (Normal or Perturbative) multistrange partice production.
!     Perturbative option is disabled for now:
  iperts = 0
!      READ (13, *) MANYB, B0, BI, BM
!     2/04/00 MANYB MUST BE SET TO 1 !
!     in order to skip impact parameter setting by ART, then B0 has no effect.
  manyb = 1
  b0 = 1
  bi = 0
  bm = 0
!      write(12,*) manyb, b0, bi, bm, ' manyb, b0, bi, bm'
!      READ (13, *) ISEED
!      write(12,*) iseed, ' iseed'
!      READ (13, *) DT
!      write(12,*) dt, ' dt'
!      READ (13, *) NTMAX
!      write(12,*) ntmax, ' ntmax'
!      READ (13, *) ICOLL
  icoll = -1
!      write(12,*) icoll, ' icoll'
!      READ (13, *) NUM
!     2/11/03 run events without test particles for now:
  num = 1
!      write(12,*) num, ' num'
!      READ (13, *) INSYS
  insys = 1
!      write(12,*) insys, ' insys'
!      READ (13, *) IPOT
  ipot = 3
!      write(12,*) ipot, ' ipot'
!      READ (13, *) MODE
  mode = 0
  If (icoll==-1) ipot = 0
!      write(12,*) mode, ' mode'
!      READ (13, *) DX, DY, DZ
  dx = 2.73
  dy = 2.73
  dz = 2.73
!      write(12,*) dx,dy,dz,' dx,dy,dz'
!      READ (13, *) DPX, DPY, DPZ
  dpx = 0.6
  dpy = 0.6
  dpz = 0.6
!      write(12,*) dpx,dpy,dpz,' dpx,dpy,dpz'
!      READ (13, *) IAVOID
  iavoid = 1
!      write(12,*) iavoid, ' iavoid'
!      READ (13, *) IMOMEN
  imomen = 1
!      write(12,*) imomen, ' imomen'
  If (icoll==-1) imomen = 3
!      READ (13, *) NFREQ
  nfreq = 10
!      write(12,*) nfreq, ' nfreq'
!      READ (13, *) ICFLOW
  icflow = 0
!      write(12,*) ICFLOW, ' ICFLOW'
!      READ (13, *) ICRHO
  icrho = 0
!      write(12,*) ICRHO, ' ICRHO'
!      READ (13, *) ICOU
  icou = 0
!      write(12,*)icou, ' icou'
! kaon potential control parameter
! KMUL IS A MULTIPLIER TO THE STANDARD K-N SCATTERING LENGTH
!      READ (13, *) KPOTEN, KMUL
  kpoten = 0
  kmul = 1
!      write(12,*)kpoten,kmul, ' kpoten, kmul'
! mean field control parameter FOR BARYONS
! no mean filed is used for baryons if their
! local density is higher than dencut.
!      READ (13, *) DENCUT
  dencut = 15
!      write(12,*)dencut, ' dencut'
! test reactions in a box of side-length cycbox
! input cycbox
!      READ (13, *) CYCBOX
  cycbox = 0
!      write(12,*) cycbox, ' cycbox'
!
!lin-5b/2008
!      if(ioscar.eq.2) then
  If (ioscar==2 .Or. ioscar==3) Then
    Open (92, File='ana/parton-initial-afterPropagation.dat', Status='UNKNOWN')
  End If
  If (ioscar==3) Then
!lin-6/2009 write out full parton collision history:
    Open (95, File='ana/parton-collisionsHistory.dat', Status='UNKNOWN')
!lin-6/2009 write out initial minijet information:
    Open (96, File='ana/minijet-initial-beforePropagation.dat', Status='UNKNOWN')
!lin-6/2009 write out parton info after coalescence:
    If (isoft==4 .Or. isoft==5) Then
      Open (85, File='ana/parton-after-coalescence.dat', Status='UNKNOWN')
    End If
  End If
!lin-6/2009 write out initial transverse positions of initial nucleons:
  Open (94, File='ana/npart-xy.dat', Status='UNKNOWN')
!
!lin-8/2009 In case that random positions are used to embed high-Pt jets:
  If (iembed==3 .Or. iembed==4) Then
    Open (97, File='embed-jet-xy.txt', Status='UNKNOWN')
    Read (97, *) nxyjet
!     Save positions in array to reuse when embedding more jet pairs
!     than the number of entries in the position file:
    If (nevent>nxyjet) Then
      If (nxyjet>nxymax) Then
        Print *, 'Too many lines in embed-jet-xy.txt:             increase value of the parameter nxymax'
        Stop
      Else If (nxyjet<=0) Then
        Print *, 'Check number of entries in embed-jet-xy.txt'
        Stop
      End If
      Do ixy = 1, nxyjet
        Read (97, *) xyjet(ixy, 1), xyjet(ixy, 2)
      End Do
    End If
  End If

  Return
End Subroutine artset

!-----------------------------------------------------------------------

!.....subroutine to initialize cascade.

Subroutine arini

!.....before invoking ARINI:
!.....IAPAR2(1), IAINT2(1) must be set.
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
!c      SAVE /ARPRNT/
  Save

!test off for resonance (phi, K*) studies:
!      OPEN (89, FILE = 'ana/decay_rec.dat', STATUS = 'UNKNOWN')

  iflg = iapar2(1)
  Goto (200, 200, 300) iflg

!.....error choice of initialization
  Print *, 'IAPAR2(1) must be 1, 2, or 3'
  Stop

!.....to use default initial conditions generated by the cascade,
!.....or to read in initial conditions.
  200 Return

!.....to generate formation time and the position at formation time from
!.....read-in initial conditions with an averaged formation proper time.
  300 Call arini1
!.....ordering the particle label according to increasing order of
!.....formation time.
  Call artord
  Return

End Subroutine arini

!-----------------------------------------------------------------------

!.....subroutine to generate formation time and position at formation time
!.....from read-in initial conditions with an averaged formation proper
!.....time.

Subroutine arini1

!.....before invoking ARINI1:
!.....ARPAR1(1), IAINT2(1) must be set:
  Parameter (maxstr=150001)
  Double Precision smearp, smearh

  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
!c      SAVE /ARPRNT/
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
!c      SAVE /ARPRC/
  Common /smearz/smearp, smearh
!c      SAVE /smearz/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /anim/nevent, isoft, isflag, izpc
!c      SAVE /anim/
  Common /nzpc/nattzp
!c      SAVE /nzpc/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Common /para8/idpert, npertd, idxsec
  Save

!lin-5/2008 for perturbatively-produced hadrons (currently only deuterons):
  Open (91, File='ana/deuteron_processes.dat', Status='UNKNOWN')
  If (idpert==1 .Or. idpert==2) Then
    Open (90, File='ana/ampt_pert.dat', Status='UNKNOWN')
  End If
!.....generate formation time and position at formation time.
  tau0 = arpar1(1)
  np = iaint2(1)
!lin-7/10/01     initial positions already given for hadrons
!     formed from partons inside ZPC (from string melting):
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
!lin-8/2015 fixed a bug that may skip "dpertp(I)=1." in addhad and
!     cause the first few events to be missing in ampt.dat
!     (mostly for low-multiplicity events such as PP collisions):
!         if(NP.le.nattzp) return
    If (np>nattzp) Then

      Do i = nattzp + 1, np
!lin-9/2012 determine rapidity more generally
!     to prevent overflow when Pt~=0 and E=|Pz|:
!            IF (ABS(PZAR(I)) .GE. PEAR(I)) THEN
!               PRINT *, ' IN ARINI1'
!               PRINT *, 'ABS(PZ) .GE. EE for particle ', I
!               PRINT *, ' FLAV = ', ITYPAR(I), ' PX = ', PXAR(I),
!     &              ' PY = ', PYAR(I)
!               PRINT *, ' PZ = ', PZAR(I), ' EE = ', PEAR(I)
!               PRINT *, ' XM = ', XMAR(I)
!               RAP = 1000000.0
!               GOTO 50
!            END IF
!c            RAP=0.5*LOG((PEAR(I)+PZAR(I))/(PEAR(I)-PZAR(I)))
!            RAP=0.5*LOG((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
! 50         CONTINUE
        If ((xmar(i)**2+pxar(i)**2+pyar(i)**2)>0.) Then
          rap = asinh(pzar(i)/sqrt(xmar(i)**2+pxar(i)**2+pyar(i)**2))
        Else
          Print *, ' IN ARINI1 mt=0'
          rap = 1000000.0*sign(1., pzar(i))
        End If

        vx = pxar(i)/pear(i)
        vy = pyar(i)/pear(i)
        ftar(i) = tau0*cosh(rap)
        gxar(i) = gxar(i) + vx*ftar(i)
        gyar(i) = gyar(i) + vy*ftar(i)
        gzar(i) = tau0*sinh(rap)
!lin-5/2009 No formation time for spectator projectile or target nucleons:
        If (pxar(i)==0 .And. pyar(i)==0 .And. (itypar(i)==2112 .Or. itypar(i)==2212)) Then
!lin-2/2013 for spectator target nucleons in LAB frame:
!     1           .and.(PEAR(I)*2/HINT1(1)).gt.0.99
          If ((pear(i)/hint1(6)>0.99 .And. pear(i)/hint1(6)<1.01) .Or. (pear(i)/hint1(7)>0.99 .And. pear(i)/hint1(7)<1.01)) Then
!
            taui = 1.E-20
            ftar(i) = taui*cosh(rap)
            gzar(i) = taui*sinh(rap)
          End If
        End If
      End Do
!lin-8/2015:
    End If
!lin-7/10/01-end
!lin-3/2009 cleanup of program flow:
  Else
    Do i = 1, np
!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZAR(I)) .GE. PEAR(I)) THEN
!               PRINT *, ' IN ARINI1'
!               PRINT *, 'ABS(PZ) .GE. EE for particle ', I
!               PRINT *, ' FLAV = ', ITYPAR(I), ' PX = ', PXAR(I),
!     &              ' PY = ', PYAR(I)
!               PRINT *, ' PZ = ', PZAR(I), ' EE = ', PEAR(I)
!               PRINT *, ' XM = ', XMAR(I)
!               RAP = 1000000.0
!               GOTO 100
!c               STOP
!            END IF
! 100        CONTINUE
!            RAP=0.5*LOG((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
      If ((xmar(i)**2+pxar(i)**2+pyar(i)**2)>0.) Then
        rap = asinh(pzar(i)/sqrt(xmar(i)**2+pxar(i)**2+pyar(i)**2))
      Else
        Print *, ' IN ARINI1 mt=0'
        rap = 1000000.0*sign(1., pzar(i))
      End If

      vx = pxar(i)/pear(i)
      vy = pyar(i)/pear(i)
!.....give initial formation time shift
      taui = ftar(i) + tau0
      ftar(i) = taui*cosh(rap)
      gxar(i) = gxar(i) + vx*tau0*cosh(rap)
      gyar(i) = gyar(i) + vy*tau0*cosh(rap)
!     4/25/03: hadron z-position upon formation determined the same way as x,y:
      gzar(i) = taui*sinh(rap)
!     the old prescription:
!            GZAR(I) = GZAR(I) + TAU0 * SINH(RAP)
      zsmear = sngl(smearh)*(2.*ranart(nseed)-1.)
      gzar(i) = gzar(i) + zsmear
!bz1/28/99end
!     10/05/01 no formation time for spectator projectile or target nucleons:
      If (pxar(i)==0 .And. pyar(i)==0 .And. (itypar(i)==2112 .Or. itypar(i)==2212)) Then
!lin-2/2013 for spectator target nucleons in LAB frame:
!     1           .and.(PEAR(I)*2/HINT1(1)).gt.0.99
        If ((pear(i)/hint1(6)>0.99 .And. pear(i)/hint1(6)<1.01) .Or. (pear(i)/hint1(7)>0.99 .And. pear(i)/hint1(7)<1.01)) Then
!
!lin-5/2008:
!               TAUI=0.00001
          taui = 1.E-20
          ftar(i) = taui*cosh(rap)
          gzar(i) = taui*sinh(rap) + zsmear
        End If
      End If
    End Do
!lin-3/2009 cleanup of program flow:
  End If

!lin-3/2009 Add initial hadrons before the hadron cascade starts:
  Call addhad

  Return
End Subroutine arini1

!-----------------------------------------------------------------------

!.....subroutine to order particle labels according to increasing
!.....formation time

Subroutine artord

!.....before invoking ARTORD:
!.....IAINT2(1) must be set:
  Parameter (maxstr=150001, maxr=1)
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
!c      SAVE /ARPRNT/
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
!c      SAVE /ARPRC/
!lin-3/2009 Take care of particle weights when user inserts initial hadrons:
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Dimension dptemp(maxstr)
!
  Dimension ityp0(maxstr), gx0(maxstr), gy0(maxstr), gz0(maxstr), ft0(maxstr), px0(maxstr), py0(maxstr), pz0(maxstr), ee0(maxstr), xm0(maxstr)
  Dimension indx(maxstr)
  External arindx
  Save
!
  npar = 0
  np = iaint2(1)
  Do i = 1, np
    ityp0(i) = itypar(i)
    gx0(i) = gxar(i)
    gy0(i) = gyar(i)
    gz0(i) = gzar(i)
    ft0(i) = ftar(i)
    px0(i) = pxar(i)
    py0(i) = pyar(i)
    pz0(i) = pzar(i)
    ee0(i) = pear(i)
    xm0(i) = xmar(i)
!lin-3/2009:
    dptemp(i) = dpertp(i)
  End Do
  Call arindx(maxstr, np, ft0, indx)
  Do i = 1, np
!bz12/3/98
!         IF (ITYP0(INDX(I)) .EQ. 211) THEN
!         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 321) THEN
!         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 2212 .OR.
!     &      ITYP0(INDX(I)) .EQ. 2112 .OR. ITYP0(INDX(I)) .EQ. -211 .OR.
!     &      ITYP0(INDX(I)) .EQ. 111) THEN
!         IF (ITYP0(INDX(I)) .EQ. 211 .OR. ITYP0(INDX(I)) .EQ. 2212 .OR.
!     &      ITYP0(INDX(I)) .EQ. 2112) THEN
    npar = npar + 1
!         ITYPAR(I) = ITYP0(INDX(I))
!         GXAR(I) = GX0(INDX(I))
!         GYAR(I) = GY0(INDX(I))
!         GZAR(I) = GZ0(INDX(I))
!         FTAR(I) = FT0(INDX(I))
!         PXAR(I) = PX0(INDX(I))
!         PYAR(I) = PY0(INDX(I))
!         PZAR(I) = PZ0(INDX(I))
!         PEAR(I) = EE0(INDX(I))
!         XMAR(I) = XM0(INDX(I))
    itypar(npar) = ityp0(indx(i))
    gxar(npar) = gx0(indx(i))
    gyar(npar) = gy0(indx(i))
    gzar(npar) = gz0(indx(i))
    ftar(npar) = ft0(indx(i))
    pxar(npar) = px0(indx(i))
    pyar(npar) = py0(indx(i))
    pzar(npar) = pz0(indx(i))
    pear(npar) = ee0(indx(i))
    xmar(npar) = xm0(indx(i))
!lin-3/2009:
    dpertp(npar) = dptemp(indx(i))
!         END IF
!bz12/3/98end
  End Do
  iaint2(1) = npar
!
  Return
End Subroutine artord

!-----------------------------------------------------------------------

!.....subroutine to copy individually generated particle record into
!.....particle record for many test particle runs.

Subroutine arini2(k)

  Parameter (maxstr=150001, maxr=1)
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
!c      SAVE /ARPRNT/
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
!c      SAVE /ARPRC/
  Common /arerc1/multi1(maxr)
!c      SAVE /ARERC1/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
!c      SAVE /ARPRC1/
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
!c      SAVE /tdecay/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
!c      SAVE /INPUT2/
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
!c      SAVE /RNDF77/
  Save

  multi1(k) = iaint2(1)
  Do i = 1, multi1(k)
    ityp1(i, k) = itypar(i)
    gx1(i, k) = gxar(i)
    gy1(i, k) = gyar(i)
    gz1(i, k) = gzar(i)
    ft1(i, k) = ftar(i)
    px1(i, k) = pxar(i)
    py1(i, k) = pyar(i)
    pz1(i, k) = pzar(i)
    ee1(i, k) = pear(i)
    xm1(i, k) = xmar(i)
!lin-3/2009 hadron weights are initialized in addhad():
!lin-5/2008 all hadrons not perturbatively-produced have the weight of 1:
!         dpp1(I,K)=1.
    dpp1(i, k) = dpertp(i)
  End Do

!     initialize final time of each particle to ntmax*dt except for
!     decay daughters, which have values given by tfdcy() and >(ntmax*dt):
  Do ip = 1, maxstr
    tfdcy(ip) = ntmax*dt
    tft(ip) = ntmax*dt
  End Do
!
  Do irun = 1, maxr
    Do ip = 1, maxstr
      tfdpi(ip, irun) = ntmax*dt
    End Do
  End Do

  Return
End Subroutine arini2

!=======================================================================

!.....function to convert PDG flavor code into ART flavor code.

Function iarflv(ipdg)

  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

!.....anti-Delta-
  If (ipdg==-1114) Then
    iarflv = -6
    Return
  End If

!.....anti-Delta0
  If (ipdg==-2114) Then
    iarflv = -7
    Return
  End If

!.....anti-Delta+
  If (ipdg==-2214) Then
    iarflv = -8
    Return
  End If

!.....anti-Delta++
  If (ipdg==-2224) Then
    iarflv = -9
    Return
  End If

!bzdbg2/23/99
!.....anti-proton
  If (ipdg==-2212) Then
    iarflv = -1
    Return
  End If

!.....anti-neutron
  If (ipdg==-2112) Then
    iarflv = -2
    Return
  End If
!bzdbg2/23/99end

!.....eta
  If (ipdg==221) Then
    iarflv = 0
    Return
  End If

!.....proton
  If (ipdg==2212) Then
    iarflv = 1
    Return
  End If

!.....neutron
  If (ipdg==2112) Then
    iarflv = 2
    Return
  End If

!.....pi-
  If (ipdg==-211) Then
    iarflv = 3
    Return
  End If

!.....pi0
  If (ipdg==111) Then
    iarflv = 4
    Return
  End If

!.....pi+
  If (ipdg==211) Then
    iarflv = 5
    Return
  End If

!.....Delta-
  If (ipdg==1114) Then
    iarflv = 6
    Return
  End If

!.....Delta0
  If (ipdg==2114) Then
    iarflv = 7
    Return
  End If

!.....Delta+
  If (ipdg==2214) Then
    iarflv = 8
    Return
  End If

!.....Delta++
  If (ipdg==2224) Then
    iarflv = 9
    Return
  End If

!.....Lambda
  If (ipdg==3122) Then
    iarflv = 14
    Return
  End If

!.....Lambda-bar
  If (ipdg==-3122) Then
    iarflv = -14
    Return
  End If

!.....Sigma-
  If (ipdg==3112) Then
    iarflv = 15
    Return
  End If

!.....Sigma-bar
  If (ipdg==-3112) Then
    iarflv = -15
    Return
  End If

!.....Sigma0
  If (ipdg==3212) Then
    iarflv = 16
    Return
  End If

!.....Sigma0-bar
  If (ipdg==-3212) Then
    iarflv = -16
    Return
  End If

!.....Sigma+
  If (ipdg==3222) Then
    iarflv = 17
    Return
  End If

!.....Sigma+ -bar
  If (ipdg==-3222) Then
    iarflv = -17
    Return
  End If

!.....K-
  If (ipdg==-321) Then
    iarflv = 21
    Return
  End If

!.....K+
  If (ipdg==321) Then
    iarflv = 23
    Return
  End If

!.....temporary entry for K0
  If (ipdg==311) Then
    iarflv = 23
    Return
  End If

!.....temporary entry for K0bar
  If (ipdg==-311) Then
    iarflv = 21
    Return
  End If

!.....temporary entry for K0S and K0L
  If (ipdg==310 .Or. ipdg==130) Then
    r = ranart(nseed)
    If (r>0.5) Then
      iarflv = 23
    Else
      iarflv = 21
    End If
    Return
  End If

!.....rho-
  If (ipdg==-213) Then
    iarflv = 25
    Return
  End If

!.....rho0
  If (ipdg==113) Then
    iarflv = 26
    Return
  End If

!.....rho+
  If (ipdg==213) Then
    iarflv = 27
    Return
  End If

!.....omega
  If (ipdg==223) Then
    iarflv = 28
    Return
  End If

!.....phi
  If (ipdg==333) Then
    iarflv = 29
    Return
  End If

!.....K*+
  If (ipdg==323) Then
    iarflv = 30
    Return
  End If
!.....K*-
  If (ipdg==-323) Then
    iarflv = -30
    Return
  End If
!.....temporary entry for K*0
  If (ipdg==313) Then
    iarflv = 30
    Return
  End If
!.....temporary entry for K*0bar
  If (ipdg==-313) Then
    iarflv = -30
    Return
  End If

!...... eta-prime
  If (ipdg==331) Then
    iarflv = 31
    Return
  End If

!...... a1
!     IF (IPDG .EQ. 777) THEN
!        IARFLV = 32
!        RETURN
!     END IF

!... cascade-
  If (ipdg==3312) Then
    iarflv = 40
    Return
  End If

!... cascade+ (bar)
  If (ipdg==-3312) Then
    iarflv = -40
    Return
  End If

!... cascade0
  If (ipdg==3322) Then
    iarflv = 41
    Return
  End If

!... cascade0 -bar
  If (ipdg==-3322) Then
    iarflv = -41
    Return
  End If

!... Omega-
  If (ipdg==3334) Then
    iarflv = 45
    Return
  End If

!... Omega+ (bar)
  If (ipdg==-3334) Then
    iarflv = -45
    Return
  End If

!... Di-Omega
  If (ipdg==6666) Then
    iarflv = 44
    Return
  End If
! sp06/05/01 end

!lin-3/2009 keep the same ID numbers in case there are initial deuterons:
  If (ipdg==42 .Or. ipdg==-42) Then
    iarflv = ipdg
    Return
  End If

!.....other
  iarflv = ipdg + 10000

  Return
End Function iarflv

!-----------------------------------------------------------------------

!.....function to convert ART flavor code into PDG flavor code.

Function invflv(iart)

  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

!.....anti-Delta-
  If (iart==-6) Then
    invflv = -1114
    Return
  End If

!.....anti-Delta0
  If (iart==-7) Then
    invflv = -2114
    Return
  End If

!.....anti-Delta+
  If (iart==-8) Then
    invflv = -2214
    Return
  End If

!.....anti-Delta++
  If (iart==-9) Then
    invflv = -2224
    Return
  End If

!bzdbg2/23/99
!.....anti-proton
  If (iart==-1) Then
    invflv = -2212
    Return
  End If

!.....anti-neutron
  If (iart==-2) Then
    invflv = -2112
    Return
  End If
!bzdbg2/23/99end

!.....eta
  If (iart==0) Then
    invflv = 221
    Return
  End If

!.....proton
  If (iart==1) Then
    invflv = 2212
    Return
  End If

!.....neutron
  If (iart==2) Then
    invflv = 2112
    Return
  End If

!.....pi-
  If (iart==3) Then
    invflv = -211
    Return
  End If

!.....pi0
  If (iart==4) Then
    invflv = 111
    Return
  End If

!.....pi+
  If (iart==5) Then
    invflv = 211
    Return
  End If

!.....Delta-
  If (iart==6) Then
    invflv = 1114
    Return
  End If

!.....Delta0
  If (iart==7) Then
    invflv = 2114
    Return
  End If

!.....Delta+
  If (iart==8) Then
    invflv = 2214
    Return
  End If

!.....Delta++
  If (iart==9) Then
    invflv = 2224
    Return
  End If

!c.....N*(1440), N*(1535) temporary entry
!      IF (IART .GE. 10 .AND. IART .LE.13) THEN
!         INVFLV = 0
!         RETURN
!      END IF

!.....Lambda
  If (iart==14) Then
    invflv = 3122
    Return
  End If
!.....Lambda-bar
  If (iart==-14) Then
    invflv = -3122
    Return
  End If

!bz3/12/99
!.....temporary entry for Sigma's
!      IF (IART .EQ. 15) THEN
!         R = RANART(NSEED)
!         IF (R .GT. 2. / 3.) THEN
!            INVFLV = 3112
!         ELSE IF (R .GT. 1./ 3. .AND. R .LE. 2. / 3.) THEN
!            INVFLV = 3212
!         ELSE
!            INVFLV = 3222
!         END IF
!         RETURN
!      END IF

!.....Sigma-
  If (iart==15) Then
    invflv = 3112
    Return
  End If

!.....Sigma- bar
  If (iart==-15) Then
    invflv = -3112
    Return
  End If

!.....Sigma0
  If (iart==16) Then
    invflv = 3212
    Return
  End If

!.....Sigma0 -bar
  If (iart==-16) Then
    invflv = -3212
    Return
  End If

!.....Sigma+
  If (iart==17) Then
    invflv = 3222
    Return
  End If

!.....Sigma+ -bar
  If (iart==-17) Then
    invflv = -3222
    Return
  End If

!lin-2/23/03 K0S and K0L are generated at the last timestep:
!.....temporary entry for K- and K0bar
  If (iart==21) Then
!         R = RANART(NSEED)
!         IF (R .GT. 0.5) THEN
    invflv = -321
!         ELSE
!            INVFLV = -311
!            R = RANART(NSEED)
!            IF (R .GT. 0.5) THEN
!               INVFLV = 310
!            ELSE
!               INVFLV = 130
!            END IF
!         END IF
    Return
  End If

!.....temporary entry for K+ and K0
  If (iart==23) Then
!         R = RANART(NSEED)
!         IF (R .GT. 0.5) THEN
    invflv = 321
!         ELSE
!            INVFLV = 311
!            R = RANART(NSEED)
!            IF (R .GT. 0.5) THEN
!               INVFLV = 310
!            ELSE
!               INVFLV = 130
!            END IF
!         END IF
    Return
  End If

!.....K0Long:
  If (iart==22) Then
    invflv = 130
    Return
  End If
!.....K0Short:
  If (iart==24) Then
    invflv = 310
    Return
  End If

!.....rho-
  If (iart==25) Then
    invflv = -213
    Return
  End If

!.....rho0
  If (iart==26) Then
    invflv = 113
    Return
  End If

!.....rho+
  If (iart==27) Then
    invflv = 213
    Return
  End If

!.....omega
  If (iart==28) Then
    invflv = 223
    Return
  End If

!.....phi
  If (iart==29) Then
    invflv = 333
    Return
  End If

!.....temporary entry for K*+ and K*0
  If (iart==30) Then
    invflv = 323
    If (ranart(nseed)>0.5) invflv = 313
    Return
  End If

!.....temporary entry for K*- and K*0bar
  If (iart==-30) Then
    invflv = -323
    If (ranart(nseed)>0.5) invflv = -313
    Return
  End If

!... eta-prime (bar)
  If (iart==31) Then
    invflv = 331
    Return
  End If

!... a1
  If (iart==32) Then
    invflv = 777
    Return
  End If

!... cascade-
  If (iart==40) Then
    invflv = 3312
    Return
  End If

!... cascade+ (bar)
  If (iart==-40) Then
    invflv = -3312
    Return
  End If

!... cascade0
  If (iart==41) Then
    invflv = 3322
    Return
  End If

!... cascade0 -bar
  If (iart==-41) Then
    invflv = -3322
    Return
  End If

!... Omega-
  If (iart==45) Then
    invflv = 3334
    Return
  End If

!... Omega+ (bar)
  If (iart==-45) Then
    invflv = -3334
    Return
  End If

!... Di-Omega
  If (iart==44) Then
    invflv = 6666
    Return
  End If
! sp 12/19/00 end

!lin-5/2008 deuteron ID numbers in ART and ampt.dat:
  If (iart==42) Then
    invflv = 42
    Return
  Else If (iart==-42) Then
    invflv = -42
    Return
  End If
!
!.....other
  invflv = iart - 10000

  Return
End Function invflv

!=======================================================================

Block Data ardata

  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
!c      SAVE /ARPRNT/
  Save
  Data arpar1/1.19, 99*0.0/
  Data iapar2/3, 49*0/
  Data arint1/100*0.0/
  Data iaint2/50*0/

End Block Data ardata

!=======================================================================

!.....Routine borrowed from ZPC.
!.....double precision  is modified to real*4.

!bz1/29/99
!      subroutine index1(n, m, arrin, indx)
Subroutine arindx(n, m, arrin, indx)
!bz1/29/99end
!     indexes the first m elements of ARRIN of length n, i.e., outputs INDX
!     such that ARRIN(INDEX(J)) is in ascending order for J=1,...,m

!      implicit real*4 (a-h, o-z)

  Dimension arrin(n), indx(n)
  Save
  Do j = 1, m
    indx(j) = j
  End Do
  l = m/2 + 1
  ir = m
  10 Continue
  If (l>1) Then
    l = l - 1
    indxt = indx(l)
    q = arrin(indxt)
  Else
    indxt = indx(ir)
    q = arrin(indxt)
    indx(ir) = indx(1)
    ir = ir - 1
    If (ir==1) Then
      indx(1) = indxt
      Return
    End If
  End If
  i = l
  j = l + l
  20 If (j<=ir) Then
    If (j<ir) Then
      If (arrin(indx(j))<arrin(indx(j+1))) j = j + 1
    End If
    If (q<arrin(indx(j))) Then
      indx(i) = indx(j)
      i = j
      j = j + j
    Else
      j = ir + 1
    End If
    Goto 20
  End If
  indx(i) = indxt
  Goto 10

End Subroutine arindx

!-----------------------------------------------------------------------

!.....extracted from G. Song's ART expasion including K- interactions
!.....file `NEWKAON.FOR'

!     5/01/03 send iblock value into art1f.f, necessary for resonance studies:
!        subroutine newka(icase,irun,iseed,dt,nt,ictrl,i1,i2,
!     &                                   srt,pcx,pcy,pcz)
Subroutine newka(icase, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, iblock)
  Parameter (maxstr=150001, maxr=1)
  Parameter (aka=0.498)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
  Common /nn/nnn
!c      SAVE /NN/
  Common /run/num
!c      SAVE /RUN/
  Common /pa/rpion(3, maxstr, maxr)
!c      SAVE /PA/
  Common /pb/ppion(3, maxstr, maxr)
!c      SAVE /PB/
  Common /pc/epion(maxstr, maxr)
!c      SAVE /PC/
  Common /pd/lpion(maxstr, maxr)
!c      SAVE /PD/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
  Logical lb1bn, lb2bn, lb1mn, lb2mn
!bz3/7/99 neutralk
!        logical lb1bn1, lb2bayon1, lb1bn0, lb2bn0
  Logical lb1bn1, lb2bn1, lb1bn0, lb2bn0
!bz3/7/99 neutralk end
  Logical lb1mn0, lb2mn0, lb1mn1, lb2mn1
  Logical lb1mn2, lb2mn2
  icase = -1
!        icase: flag for the type of reaction that is going to happen.
!        icase=-1,  no desired reaction, return to main program.
!              1,  NN,ND,DD
!              2,  PI+N, PI+D
!              3,  K(-) absorption.
  nchrg = -100
!        nchrg: Net charges of the two incoming particles.
  ictrl = 1
  lb1 = lb(i1)
  lb2 = lb(i2)
  em1 = e(i1)
  em2 = e(i2)
  lb1bn = lb1 == 1 .Or. lb1 == 2 .Or. (lb1>5 .And. lb1<=13)
  lb2bn = lb2 == 1 .Or. lb2 == 2 .Or. (lb2>5 .And. lb2<=13)
  lb1bn0 = lb1 == 2 .Or. lb1 == 7 .Or. lb1 == 10 .Or. lb1 == 12
  lb2bn0 = lb2 == 2 .Or. lb2 == 7 .Or. lb2 == 10 .Or. lb2 == 12
  lb1bn1 = lb1 == 1 .Or. lb1 == 8 .Or. lb1 == 11 .Or. lb1 == 13
  lb2bn1 = lb2 == 1 .Or. lb2 == 8 .Or. lb2 == 11 .Or. lb2 == 13
  lb1mn = em1 < 0.2 .Or. lb1 == 0 .Or. (lb1>=25 .And. lb1<=29)
  lb2mn = em2 < 0.2 .Or. lb2 == 0 .Or. (lb2>=25 .And. lb2<=29)
  lb1mn0 = lb1 == 0 .Or. lb1 == 4 .Or. lb1 == 26 .Or. lb1 == 28 .Or. lb1 == 29
  lb2mn0 = lb2 == 0 .Or. lb2 == 4 .Or. lb2 == 26 .Or. lb2 == 28 .Or. lb2 == 29
  lb1mn1 = lb1 == 5 .Or. lb1 == 27
  lb2mn1 = lb2 == 5 .Or. lb2 == 27
  lb1mn2 = lb1 == 3 .Or. lb1 == 25
  lb2mn2 = lb2 == 3 .Or. lb2 == 25

!        1. consider N+N, N+Resonance, R + R reactions
  If (lb1bn .And. lb2bn) Then
!     NN,ND,DD:
    icase = 1
!     total cross section
    sig = 40.
    If (lb1==9 .And. lb2==9) Then
      nchrg = 4
    End If
    If ((lb1bn1 .And. lb2==9) .Or. (lb2bn1 .And. lb1==9)) Then
      nchrg = 3
    End If
    If ((lb1bn0 .And. lb2==9) .Or. (lb2bn0 .And. lb1==9) .Or. (lb1bn1 .And. lb2bn1)) Then
      nchrg = 2
    End If
    If ((lb1bn1 .And. lb2bn0) .Or. (lb1==6 .And. lb2==9) .Or. (lb2bn1 .And. lb1bn0) .Or. (lb2==6 .And. lb1==9)) Then
      nchrg = 1
    End If
    If ((lb1bn0 .And. lb2bn0) .Or. (lb1bn1 .And. lb2==6) .Or. (lb2bn1 .And. lb1==6)) Then
      nchrg = 0
    End If
    If ((lb1bn0 .And. lb2==6) .Or. (lb2bn0 .And. lb1==6)) Then
      nchrg = -1
    End If
    If (lb1==6 .And. lb2==6) Then
      nchrg = -2
    End If
!     brsig = x2kaon_no_isospin(srt)
    If (nchrg>=-1 .And. nchrg<=2) Then
!     K,Kbar prduction x sect.
      brsig = x2kaon(srt)
    Else
      brsig = 0.0
!                if(nchrg.eq.-2.or.nchrg.eq.3) then
!                   brsig = x2kaon(srt+0.938-1.232)
!                else
!     nchrg=4
!                   brsig = x2kaon(srt+2.*(0.938-1.232))
!                endif
    End If

!bz3/7/99 neutralk
    brsig = 2.0*brsig
!bz3/7/99 neutralk end

  End If

!        2. consider PI(meson:eta,omega,rho,phi) + N(N*,D)
  If ((lb1bn .And. lb2mn) .Or. (lb2bn .And. lb1mn)) Then
!     PN,PD
    icase = 2
    sig = 20.
    sigma0 = pinsg0(srt)
    brsig = 0.0
    If ((lb1bn1 .And. lb2mn0) .Or. (lb2bn1 .And. lb1mn0) .Or. (lb1bn0 .And. lb2mn1) .Or. (lb2bn0 .And. lb1mn1) .Or. (lb1==9 .And. lb2mn2) .Or. (lb2==9 .And. lb1mn2)) Then
      nchrg = 1
!bz3/2/99/song
!                if(lb1bn1.or.lb2bn1) brsig=2.0*sigma0
!                if(lb1bn0.or.lb2bn0) brsig=0.5*sigma0
      If (lb1bn1 .Or. lb2bn1) brsig = 0.5*sigma0
      If (lb1bn0 .Or. lb2bn0) brsig = 2.0*sigma0
!bz3/2/99/song end
!                if(lb1.eq.9.or.lb2.eq.9) brsig=1.5*sigma0
    End If
    If ((lb1bn0 .And. lb2mn0) .Or. (lb2bn0 .And. lb1mn0) .Or. (lb1bn1 .And. lb2mn2) .Or. (lb2bn1 .And. lb1mn2) .Or. (lb1==6 .And. lb2mn1) .Or. (lb2==6 .And. lb1mn1)) Then
      nchrg = 0
      If (lb1bn1 .Or. lb2bn1) Then
!bz3/2/99/song
!                  brsig=1.5*sigma0
        brsig = 3.0*sigma0
!bz3/2/99/song end
!bz3/11/99/song
!                  ratiok = 1./3.
        ratiok = 2./3.
!bz3/11/99/song end

!                  ratiok: the ratio of channels: ->nK+k- vs. -> pK0K-
      End If
      If (lb1bn0 .Or. lb2bn0) Then
        brsig = 2.5*sigma0
!bz3/2/99/song
!                  ratiok = 0.8
        ratiok = 0.2
!bz3/2/99/song end
      End If
!                if(lb1.eq.6.or.lb2.eq.6) then
!     lb=6 : D-
!                  brsig=1.5*sigma0
!                  ratiok = 0.5
!                endif
    End If
    If ((lb1bn0 .And. lb2mn2) .Or. (lb2bn0 .And. lb1mn2) .Or. (lb1==6 .And. lb2mn0) .Or. (lb2==6 .And. lb1mn0)) Then
      nchrg = -1
      If (lb1bn0 .Or. lb2bn0) brsig = sigma0
!                if(lb1.eq.6.or.lb2.eq.6) brsig=sigma0
    End If
!          if((lb1.eq.6.and.lb2mn2).or.(lb2.eq.6.and.lb1mn2))then
!                nchrg=-2
!          endif
!          if((lb1bn1.and.lb2mn1).or.(lb2bn1.and.lb1mn1)
!    &           .or.(lb1.eq.9.and.lb2mn0).or.(lb2.eq.9.and.lb1mn0)) then
!                nchrg=2
!          endif

!bz3/11/99 neutralk
    If ((lb1==6 .And. lb2mn2) .Or. (lb2==6 .And. lb1mn2)) Then
      nchrg = -2
    End If
!bz3/11/99 neutralk
!bz3/8/99 neutralk
    If ((lb1bn1 .And. lb2mn1) .Or. (lb2bn1 .And. lb1mn1) .Or. (lb1==9 .And. lb2mn0) .Or. (lb2==9 .And. lb1mn0)) Then
      nchrg = 2
    End If
!bz3/8/99 neutralk end

!bz3/7/99 neutralk
    If (nchrg>=-2 .And. nchrg<=2) Then
      brsig = 3.0*sigma0
    End If
!bz3/7/99 neutralk end

  End If

!        3. consider K- + N(N*,D) absorption.
!        if((lb1bn.and.lb2.eq.21).OR.(lb2bn.and.lb1.eq.21)) then
  If ((lb1bn .And. (lb2==21 .Or. lb2==-30)) .Or. (lb2bn .And. (lb1==21 .Or. lb1==-30))) Then
!          bmass=em1+em2-aka
    bmass = 0.938
    If (srt<=(bmass+aka)) Then
!bz3/2/99
!                write(100,*)'--lb1,lb2,em1,em2,srt',lb1,lb2,em1,em2,srt
!bz3/2/99end
      pkaon = 0.
    Else
      pkaon = sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
    End If
    sig = 0.
    If (lb1==1 .Or. lb2==1 .Or. lb1==8 .Or. lb2==8 .Or. lb1==11 .Or. lb2==11 .Or. lb1==13 .Or. lb2==13) Then
!          K- + (D+,N*+)p ->
      nchrg = 0
      sigela = akpel(pkaon)
      sigsgm = 3.*akpsgm(pkaon)
      sig = sigela + sigsgm + akplam(pkaon)
    End If
    If (lb1==2 .Or. lb2==2 .Or. lb1==7 .Or. lb2==7 .Or. lb1==10 .Or. lb2==10 .Or. lb1==12 .Or. lb2==12) Then
!          K- + (D0, N*0)n ->
      nchrg = -1
      sigela = aknel(pkaon)
      sigsgm = 2.*aknsgm(pkaon)
      sig = sigela + sigsgm + aknlam(pkaon)
    End If
    If (lb1==6 .Or. lb2==6) Then
!     K- + D-
      nchrg = -2
      sigela = aknel(pkaon)
      sigsgm = aknsgm(pkaon)
      sig = sigela + sigsgm
    End If
    If (lb1==9 .Or. lb2==9) Then
!     K- + D++
      nchrg = 1
      sigela = akpel(pkaon)
      sigsgm = 2.*akpsgm(pkaon)
      sig = sigela + sigsgm + akplam(pkaon)
    End If

!bz3/8/99 neutralk
    sigela = 0.5*(akpel(pkaon)+aknel(pkaon))
    sigsgm = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
    sig = sigela + sigsgm + akplam(pkaon)
!bz3/8/99 neutralk end

    If (sig>1.E-7) Then
!     K(-) + N reactions
      icase = 3
      brel = sigela/sig
      brsgm = sigsgm/sig
!              branch_lambda=akNlam(pkaon)/sig
      brsig = sig
    End If
  End If

!        4. meson + hyperon -> K- + N
!        if(((lb1.ge.14.and.lb1.le.17).and.lb2mn).OR.
!     &     ((lb2.ge.14.and.lb2.le.17).and.lb1mn)) then
  If (((lb1>=14 .And. lb1<=17) .And. (lb2>=3 .And. lb2<=5)) .Or. ((lb2>=14 .And. lb2<=17) .And. (lb1>=3 .And. lb1<=5))) Then
!        first classify the reactions due to total charge.
    nchrg = -100
    If ((lb1==15 .And. (lb2==3 .Or. lb2==25)) .Or. (lb2==15 .And. (lb1==3 .Or. lb1==25))) Then
      nchrg = -2
!     D-
      bmass = 1.232
    End If
    If ((lb1==15 .And. lb2mn0) .Or. (lb2==15 .And. lb1mn0) .Or. ((lb1==14 .Or. lb1==16) .And. (lb2==3 .Or. lb2==25)) .Or. ((lb2==14 .Or. lb2==16) .And. (lb1==3 .Or. lb1==25))) Then
      nchrg = -1
!     n
      bmass = 0.938
    End If
    If ((lb1==15 .And. (lb2==5 .Or. lb2==27)) .Or. (lb2==15 .And. (lb1==5 .Or. lb1==27)) .Or. (lb1==17 .And. (lb2==3 .Or. lb2==25)) .Or. (lb2==17 .And. (lb1==3 .Or. lb1==25)) .Or. ((lb1==14 .Or. lb1==16) .And. lb2mn0) .Or. ((lb2==14 .Or. lb2==16) .And. lb1mn0)) Then
      nchrg = 0
!     p
      bmass = 0.938
    End If
    If ((lb1==17 .And. lb2mn0) .Or. (lb2==17 .And. lb1mn0) .Or. ((lb1==14 .Or. lb1==16) .And. (lb2==5 .Or. lb2==27)) .Or. ((lb2==14 .Or. lb2==16) .And. (lb1==5 .Or. lb1==27))) Then
      nchrg = 1
!     D++
      bmass = 1.232
    End If
    sig = 0.
    If (nchrg/=-100 .And. srt>(aka+bmass)) Then
!     PI+sigma or PI + Lambda => Kbar + N reactions
      icase = 4
!             pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
      pkaon = sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
!     lambda + Pi
      If (lb1==14 .Or. lb2==14) Then
        If (nchrg>=0) sigma0 = akplam(pkaon)
        If (nchrg<0) sigma0 = aknlam(pkaon)
!     sigma + pi
      Else
!     K-p or K-D++
        If (nchrg>=0) sigma0 = akpsgm(pkaon)
!     K-n or K-D-
        If (nchrg<0) sigma0 = aknsgm(pkaon)

!bz3/8/99 neutralk
        sigma0 = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
!bz3/8/99 neutralk end

      End If
      sig = (srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
!bz3/8/99 neutralk
!     if(nchrg.eq.-2.or.nchrg.eq.1) sig=2.*sig K-D++, K-D-
!     K0barD++, K-D-
      If (nchrg==-2 .Or. nchrg==2) sig = 2.*sig

!bz3/8/99 neutralk end

!             the factor 2 comes from spin of delta, which is 3/2
!             detailed balance. copy from Page 423 of N.P. A614 1997

!bz3/8/99 neutralk
      If (lb1==14 .Or. lb2==14) Then
        sig = 4.0/3.0*sig
      Else If (nchrg==-2 .Or. nchrg==2) Then
        sig = 8.0/9.0*sig
      Else
        sig = 4.0/9.0*sig
      End If
!bz3/8/99 neutralk end
      brsig = sig
      If (sig<1.E-7) sig = 1.E-7
    End If
!sp05/07/01
! comment icase=4 statement below if only inelastic
!     PI+L/Si => Kbar + N  OR ELASTIC SCATTERING
    icase = 4
    brsig = sig
!     elastic xsecn of 10mb
    sigela = 10.
    sig = sig + sigela
    brel = sigela/sig
!c          brsig = sig
!sp05/07/01 end
  End If
!
!        if(em2.lt.0.2.and.em1.lt.0.2) then
!     PI + PI
!             icase=5
!     assumed PI PI total x section.
!              sig=50.
!     Mk + Mkbar
!              s0=aka+aka
!              brsig = 0.
!              if(srt.gt.s0) brsig = 2.7*(1.-s0**2/srt**2)**0.76
!              x section for PIPI->KKbar   PRC43 (1991) 1881
!        endif
  If (icase==-1) Then
    ictrl = -1
    Return
  End If
  px1cm = pcx
  py1cm = pcy
  pz1cm = pcz
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
!        ec=3.59709
!        if((e(i1).ge.1.).and.(e(i2).ge.1.)) ec = 4.75

  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
  If (ic==-1) Then
!     no anti-kaon production
    ictrl = -1
!           in=in+1
!           write(60,*)'--------------distance-----',in
    Return
  End If

!lin-10/24/02 set to 0: ik,ik0-3,il,im,im3-4,in,inpion,ipipi,
!     sgsum,sgsum1,sgsum3:
  ik = 0
  ik0 = 0
  ik1 = 0
  ik2 = 0
  ik3 = 0
  il = 0
  im = 0
  im3 = 0
  im4 = 0
  in = 0
  inpion = 0
  ipipi = 0
  sgsum = 0.
  sgsum1 = 0.
  sgsum3 = 0.
  If (icase==1) Then
    ik = ik + 1
    If (srt>2.8639) Then
      ik0 = ik0 + 1
      If (em1<1.0 .And. em2<1.0) Then
        ik1 = ik1 + 1
        sgsum1 = sgsum1 + brsig
!                        ratio_1=sgsum1/ik1/40.
      End If
      If (em1>1.0 .And. em2>1.0) Then
        ik3 = ik3 + 1
        sgsum3 = sgsum3 + brsig
!                        ratio_3=sgsum3/ik3/40.
      End If
      If (em1>1.0 .And. em2<1.0) ik2 = ik2 + 1
      If (em1<1.0 .And. em2>1.0) ik2 = ik2 + 1
      sgsum = sgsum + brsig
!                ratio=sgsum/ik0/40.
    End If
  End If
  If (icase==2) inpion = inpion + 1
  If (icase==5) ipipi = ipipi + 1
!        write(62,*)'ik1,ik2,ik3',ik1,ik2,ik3,ratio_1,ratio_3,ratio
!        write(62,*)'inpion,ipipi',inpion,ipipi
  If (ranart(nseed)>(brsig/sig)) Then
!     no anti-kaon production
    ictrl = -1
    Return
  End If
  il = il + 1
!        kaons could be created now.
  If (icase==1) Then
    in = in + 1
!          write(60,*)'------in,s2kaon,sig=',in,brsig,sig,lb1,lb2
    Call nnkaon(irun, iseed, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg)
  End If
  If (icase==2) Then
    im = im + 1
!          call npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
!     &              pcx,pcy,pcz,nchrg,ratiok)
    Call npik(irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg, ratiok, iblock)
  End If
!
  If (icase==3) Then
    im3 = im3 + 1
!          write(63,*)'im3,lb1,lb2,pkaon',im3,lb1,lb2,pkaon
!          write(63,*)'sig,el,sigma',sig,brel,brsgm
!          write(63,*)'srt,pcx,pcy,pcz,em1,em2',srt,pcx,pcy,pcz,em1,em2
    Call kaonn(brel, brsgm, irun, iseed, dt, nt, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg)
!         this subroutine format is diff. since three final states are possible
  End If
!

  If (icase==4) Then
    im4 = im4 + 1
!          write(64,*)'im4,sigma0,branch,sig=',im4,sigma0,brsig,sig
!          write(64,*)'lb1,lb2,em1,em2,pkaon=',lb1,lb2,em1,em2,pkaon

!sp06/07/01
    If (ranart(nseed)<brel) Then
      ielstc = 1
    Else
      ielstc = 0
    End If
!          call Pihypn(ielstc,irun,iseed,dt,nt,ictrl,i1,i2,srt,
!     &                   pcx,pcy,pcz,nchrg)
    Call pihypn(ielstc, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg, iblock)

!sp06/07/01 end
  End If
!        if(icase.eq.5) then
!          im5=im5+1
!          write(65,*)'---im5,s2kaon,sig=',im5,brsig,sig
!          call pipikaon(irun,iseed,dt,nt,ictrl,i1,i2,srt,pcx,pcy,pcz)
!        endif
!bz3/2/99
!        write(101,*)lb1,lb2,lb(i1),lb(i2)
!        write(101,*)em1,em2,e(i1),e(i2),srt
!bz3/2/99end

  Return
End Subroutine newka

!*****************************************
! for pp-->pp + kaon + anti-kaon
!      real*4 function X2kaon(srt)
Real Function x2kaon(srt)
  Save
!  This function contains the experimental total pp->pp+K(+)K(-) Xsections    *
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in mb                                    *
!                                                                             *
!*****************************************
!     minimum c.m.s. energy to create 2 kaon. = 2*(mp+mk)
  smin = 2.8639
  x2kaon = 0.0000001
  If (srt<smin) Return
  sigma1 = 2.8
  sigma2 = 7.7
  sigma3 = 3.9
  x = srt**2/smin**2 + 0.0000001
  f1 = (1.+1./sqrt(x))*alog(x) - 4.*(1.-1./sqrt(x))
  f2 = 1. - (1./sqrt(x))*(1.+alog(sqrt(x)))
  f3 = ((x-1.)/x**2)**3.5
  x2kaon = (1.-1./x)**3*(sigma1*f1+sigma2*f2) + sigma3*f3
  Return
End Function x2kaon

Real Function pinsg0(srt)
  Save
! cross section in mb for PI- + P -> P + K0 + K-
!     Mn + 2* Mk
  srt0 = 0.938 + 2.*0.498
  If (srt<srt0) Then
    pinsg0 = 0.0
    Return
  End If
  ratio = srt0**2/srt**2
  pinsg0 = 1.121*(1.-ratio)**1.86*ratio**2
  Return
End Function pinsg0

Real Function aknel(pkaon)
  Save
!cross section in mb for K- + N reactions.
!        the following data come from PRC 41 (1701)
!        sigma1: K(-) + neutron elastic
  If (pkaon<0.5 .Or. pkaon>=4.0) sigma1 = 0.
  If (pkaon>=0.5 .And. pkaon<1.0) sigma1 = 20.*pkaon**2.74
  If (pkaon>=1.0 .And. pkaon<4.0) sigma1 = 20.*pkaon**(-1.8)
  aknel = sigma1
  Return
End Function aknel

Real Function akpel(pkaon)
  Save
!cross section in mb for K- + N reactions.
!        the following data come from PRC 41 (1701)
!        sigma2: K(-) + proton elastic
  If (pkaon<0.25 .Or. pkaon>=4.0) sigma2 = 0.
  If (pkaon>=0.25 .And. pkaon<4.0) sigma2 = 13.*pkaon**(-0.9)
  akpel = sigma2
  Return
End Function akpel

Real Function aknsgm(pkaon)
  Save
!cross section in mb for K- + N reactions.
!        sigma2: x section for K- + n -> sigma0 + PI-
  If (pkaon<0.5 .Or. pkaon>=6.0) sigma2 = 0.
  If (pkaon>=0.5 .And. pkaon<1.0) sigma2 = 1.2*pkaon**(-1.3)
  If (pkaon>=1.0 .And. pkaon<6.0) sigma2 = 1.2*pkaon**(-2.3)
  aknsgm = sigma2
  Return
End Function aknsgm

Real Function akpsgm(pkaon)
  Save
!cross section in mb for K- + N reactions.
!        sigma1: x section for K- + p -> sigma0 + PI0
  If (pkaon<0.2 .Or. pkaon>=1.5) sigma1 = 0.
  If (pkaon>=0.2 .And. pkaon<1.5) sigma1 = 0.6*pkaon**(-1.8)
  akpsgm = sigma1
  Return
End Function akpsgm

Real Function akplam(pkaon)
  Save
!cross section in mb for K- + N reactions.
!        sigma: x section for K- + p -> lambda + PI0
  p = pkaon
  If (pkaon<0.2 .Or. pkaon>=10.0) sigma = 0.
  If (pkaon>=0.2 .And. pkaon<0.9) sigma = 50.*p**2 - 67.*p + 24.
  If (pkaon>=0.9 .And. pkaon<10.0) sigma = 3.0*pkaon**(-2.6)
  akplam = sigma
  Return
End Function akplam

Real Function aknlam(pkaon)
  Save
!cross section in mb for K- + N reactions.
  aknlam = akplam(pkaon)
  Return
End Function aknlam

! GQ Li parametrization (without resonance)
Real Function aknpsg(pkaon)
  Save
!cross section in mb for K- + N reactions.
!       sigma1: x section for K- + p/n -> sigma0 + PI0
  If (pkaon<=0.345) Then
    sigma1 = 0.624*pkaon**(-1.83)
  Else
    sigma1 = 0.7*pkaon**(-2.09)
  End If
  aknpsg = sigma1
  Return
End Function aknpsg

!-----------------------------------------------------------------------

!.....extracted from G. Song's ART expasion including K- interactions
!.....file `NEWNNK.FOR'

Subroutine nnkaon(irun, iseed, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg)
!        <pt>=0.27+0.037*log(srt) was changed to 0.632 + ... on Aug. 14, 1997
!     CANCELED also alpha=1 changed to alpha=3 to decrease the leadng effect.
  Parameter (maxstr=150001, maxr=1)
  Parameter (aka=0.498)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
  Common /nn/nnn
!c      SAVE /NN/
  Common /run/num
!c      SAVE /RUN/
  Common /pa/rpion(3, maxstr, maxr)
!c      SAVE /PA/
  Common /pb/ppion(3, maxstr, maxr)
!c      SAVE /PB/
  Common /pc/epion(maxstr, maxr)
!c      SAVE /PC/
  Common /pd/lpion(maxstr, maxr)
!c      SAVE /PD/
  Dimension px(4), py(4), pz(4)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
!      dm1=e(i1)
!      dm2=e(i2)
  dm3 = 0.938
  dm4 = 0.938
!     10/24/02 initialize n to 0:
  n = 0

!bz3/11/99 neutralk
!        if(nchrg.eq.-2.or.nchrg.ge.3) dm3=1.232
!        if(nchrg.eq.4) dm4=1.232
  If (nchrg<=-1 .Or. nchrg>=3) dm3 = 1.232
  If (nchrg==-2 .Or. nchrg==4) dm4 = 1.232
!bz3/11/99 neutralk end
  iblock = 0
  Call fstate(iseed, srt, dm3, dm4, px, py, pz, iflag)
  If (iflag<0) Then
!           write(60,*)'------------final state fail-------',n
!     no anti-kaon production
    ictrl = -1
    n = n + 1
    Return
  End If
  iblock = 12
! Rotate the momenta of particles in the cms of I1 & I2
! px(1), py(1), pz(1): momentum of I1
! px(2), py(2), pz(2): momentum of I2
! px(3), py(3), pz(3): momentum of anti-kaon
! px(4), py(4), pz(4): momentum of kaon


!     10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = px(1)
  pyrota = py(1)
  pzrota = pz(1)
!        call rotate(pcx,pcy,pcz,px(1),py(1),pz(1))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(1) = pxrota
  py(1) = pyrota
  pz(1) = pzrota
!
  pxrota = px(2)
  pyrota = py(2)
  pzrota = pz(2)
!        call rotate(pcx,pcy,pcz,px(2),py(2),pz(2))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(2) = pxrota
  py(2) = pyrota
  pz(2) = pzrota
!
  pxrota = px(3)
  pyrota = py(3)
  pzrota = pz(3)
!        call rotate(pcx,pcy,pcz,px(3),py(3),pz(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(3) = pxrota
  py(3) = pyrota
  pz(3) = pzrota
!
  pxrota = px(4)
  pyrota = py(4)
  pzrota = pz(4)
!        call rotate(pcx,pcy,pcz,px(4),py(4),pz(4))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  px(4) = pxrota
  py(4) = pyrota
  pz(4) = pzrota

  nnn = nnn + 2
!     K+
  lpion(nnn, irun) = 23
  If (nchrg==-1 .Or. nchrg==-2) Then
!        To keep charge conservation. D-n->nnK0K-, D-D- -> nD-K0K-

!bz3/7/99 neutralk
!           lpion(nnn,irun)=24 ! K0
!bz3/7/99 neutralk end

  End If
!     aka: rest mass of K
  epion(nnn, irun) = aka
!     K-
  lpion(nnn-1, irun) = 21
!     aka: rest mass of K
  epion(nnn-1, irun) = aka
! Find the momenta of particles in the final state in the nucleus_nucleus
! cms frame.   Lorentz transformation into lab frame.
  e1cm = sqrt(dm3**2+px(1)**2+py(1)**2+pz(1)**2)
  p1beta = px(1)*betax + py(1)*betay + pz(1)*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px(1)
  pt2i1 = betay*transf + py(1)
  pt3i1 = betaz*transf + pz(1)
  eti1 = dm3
!        lb1   = lb(i1)
  lb1 = 2
  If (nchrg>=-2 .And. nchrg<=1) lb1 = 2

!bz3/7/99 neutralk
  If (nchrg==-2 .Or. nchrg==-1) Then
    lb1 = 6
  End If
!bz3/7/99 neutralk end

!bz3/11/99 neutralk
!        if(nchrg.eq.2.or.nchrg.eq.3) lb1=1
!        if(nchrg.eq.4) lb1=9
  If (nchrg==1 .Or. nchrg==2) lb1 = 1
  If (nchrg==3 .Or. nchrg==4) lb1 = 9
!bz3/11/99 neutralk end

! For second nulceon, same
  e2cm = sqrt(dm4**2+px(2)**2+py(2)**2+pz(2)**2)
  p2beta = px(2)*betax + py(2)*betay + pz(2)*betaz
  transf = gamma*(gamma*p2beta/(gamma+1)+e2cm)
  pt1i2 = betax*transf + px(2)
  pt2i2 = betay*transf + py(2)
  pt3i2 = betaz*transf + pz(2)
  eti2 = dm4
!        lb2   = lb(i2)
  lb2 = 2

!bz3/11/99 neutralk
!        if(nchrg.eq.-1.or.nchrg.eq.0) lb2=2
!        if(nchrg.eq. 2.or.nchrg.eq.1) lb2=1
!        if(nchrg.eq. 4.or.nchrg.eq.3) lb2=9
!        if(nchrg.eq.-2) lb2=6
  If (nchrg>=-1 .Or. nchrg<=1) lb2 = 2
  If (nchrg==2 .Or. nchrg==3) lb2 = 1
  If (nchrg==4) lb2 = 9
  If (nchrg==-2) lb2 = 6
!bz3/11/99 neutralk end

!        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lb1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lb2

!                px1 = p(1,i1)
!                py1 = p(2,i1)
!                pz1 = p(3,i1)
!                em1 = e(i1)
!                id(i1) = 2
!                id(i2) = 2
!                id1 = id(i1)
!                iblock = 101  ! K(+)K(-) production
! Get anti-kaons' momenta and coordinates in nucleus-nucleus cms. frame.
  epcmk = sqrt(epion(nnn-1,irun)**2+px(3)**2+py(3)**2+pz(3)**2)
  betak = px(3)*betax + py(3)*betay + pz(3)*betaz
  transf = gamma*(gamma*betak/(gamma+1.)+epcmk)
  ppion(1, nnn-1, irun) = betax*transf + px(3)
  ppion(2, nnn-1, irun) = betay*transf + py(3)
  ppion(3, nnn-1, irun) = betaz*transf + pz(3)
  rpion(1, nnn-1, irun) = r(1, i1)
  rpion(2, nnn-1, irun) = r(2, i1)
  rpion(3, nnn-1, irun) = r(3, i1)
!lin-5/2008:
  dppion(nnn-1, irun) = dpertp(i1)*dpertp(i2)
! Same thing for kaon **************************************
  epcmak = sqrt(epion(nnn,irun)**2+px(4)**2+py(4)**2+pz(4)**2)
  betaak = px(4)*betax + py(4)*betay + pz(4)*betaz
  transf = gamma*(gamma*betaak/(gamma+1.)+epcmak)
  ppion(1, nnn, irun) = betax*transf + px(4)
  ppion(2, nnn, irun) = betay*transf + py(4)
  ppion(3, nnn, irun) = betaz*transf + pz(4)
  rpion(1, nnn, irun) = r(1, i2)
  rpion(2, nnn, irun) = r(2, i2)
  rpion(3, nnn, irun) = r(3, i2)
!lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  Return
End Subroutine nnkaon

Subroutine lorntz(ilo, b, pi, pj)
!       It uses to perform Lorentz (or inverse Lorentz) transformation
  Dimension pi(4), pj(4), b(3)
  Save
!       dimension db(3)
  bb = b(1)*b(1) + b(2)*b(2) + b(3)*b(3)
  deno3 = sqrt(1.-bb)
  If (deno3==0.) deno3 = 1.E-10
  gam = 1./deno3
  ga = gam*gam/(gam+1.)
  If (ilo==1) Goto 100
!       Lorentz transformation
  pib = pi(1)*b(1) + pi(2)*b(2) + pi(3)*b(3)
  pjb = pj(1)*b(1) + pj(2)*b(2) + pj(3)*b(3)
!       drb=drd(1)*b(1)+drd(2)*b(2)+drd(3)*b(3)
!       drdb=db(1)*b(1)+db(2)*b(2)+db(3)*b(3)
  Do i = 1, 3
    pi(i) = pi(i) + b(i)*(ga*pib-gam*pi(4))
    pj(i) = pj(i) + b(i)*(ga*pjb-gam*pj(4))
!       drd(i)=drd(i)+b(i)*ga*drb
!       db(i)=db(i)+b(i)*ga*drdb
  End Do
  pi(4) = gam*(pi(4)-pib)
  pj(4) = gam*(pj(4)-pjb)
  Return
  100 Continue
!       inverse Lorentz transformation
  pib = pi(1)*b(1) + pi(2)*b(2) + pi(3)*b(3)
  pjb = pj(1)*b(1) + pj(2)*b(2) + pj(3)*b(3)
  Do i = 1, 3
    pi(i) = pi(i) + b(i)*(ga*pib+gam*pi(4))
    pj(i) = pj(i) + b(i)*(ga*pjb+gam*pj(4))
  End Do
  pi(4) = gam*(pi(4)+pib)
  pj(4) = gam*(pj(4)+pjb)
  Return
End Subroutine lorntz



!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine fstate(iseed, srt, dm3, dm4, px, py, pz, iflag)
!        function: decide final momentum for N,N,K(+),and K(-)
  Dimension px(4), py(4), pz(4), pe(4)
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  iflag = -1
!        iflag=-1: fail to find momenta
!             = 1: success
  pio = 3.1415926
  aka = 0.498
!        v=0.43
!        w=-0.84
!        b=3.78
!        c=0.47
!        d=3.60
!        fmax=1.056
!        gmax=1.+c

  icount = 0
  ekmax = (srt-dm3-dm4)/2.
  If (ekmax<=aka) Return
  pkmax = sqrt(ekmax**2-aka**2)

  If (dm3<=0.0 .Or. dm4<=0.0) Then
    Write (1, *) 'error: minus mass!!!'
    Return
  End If

!        after we have the momenta for both nucleus, we sample the
!        transverse momentum for K-.
!        dsigma/dpt**2 = exp(-4.145*pt**2) obtained by fitting data on
!        page 72, fig 23i.
  50 Continue
  icount = icount + 1
  If (icount>10) Return
  ptkmi2 = -1./4.145*alog(ranart(nseed))
  ptkm = sqrt(ptkmi2)
  3 v1 = ranart(nseed)
  v2 = ranart(nseed)
  rsq = v1**2 + v2**2
  If (rsq>=1.0 .Or. rsq<=0.) Goto 3
  fac = sqrt(-2.*alog(rsq)/rsq)
  guass = v1*fac
  If (guass>=5.) Goto 3
  xstar = guass/5.
  pzkm = pkmax*xstar
  ekm = sqrt(aka**2+pzkm**2+ptkm**2)
  If (ranart(nseed)>aka/ekm) Goto 50
  bbb = ranart(nseed)
  px(3) = ptkm*cos(2.*pio*bbb)
  py(3) = ptkm*sin(2.*pio*bbb)
  If (ranart(nseed)>0.5) pzkm = -1.*pzkm
  pz(3) = pzkm
  pe(3) = ekm
  150 ptkpl2 = -1./3.68*alog(ranart(nseed))
  ptkp = sqrt(ptkpl2)
  13 v1 = ranart(nseed)
  v2 = ranart(nseed)
  rsq = v1**2 + v2**2
  If (rsq>=1.0 .Or. rsq<=0.) Goto 13
  fac = sqrt(-2.*alog(rsq)/rsq)
  guass = v1*fac
  If (guass>=3.25) Goto 13
  xstar = guass/3.25
  pzkp = pkmax*xstar
  ekp = sqrt(aka**2+pzkp**2+ptkp**2)
  If (ranart(nseed)>aka/ekp) Goto 150
  bbb = ranart(nseed)
  px(4) = ptkp*cos(2.*pio*bbb)
  py(4) = ptkp*sin(2.*pio*bbb)
  If (ranart(nseed)>0.5) pzkp = -1.*pzkp
  pz(4) = pzkp
  pe(4) = ekp

  resten = srt - pe(3) - pe(4)
  restpz = -pz(3) - pz(4)
!     resample
  If (resten<=abs(restpz)) Goto 50
  restms = sqrt(resten**2-restpz**2)
!     resample
  If (restms<(dm3+dm4)) Goto 50
  ptp2 = -1./2.76*alog(ranart(nseed))
  ptp = sqrt(ptp2)
  bbb = ranart(nseed)
  px(2) = ptp*cos(2.*pio*bbb)
  py(2) = ptp*sin(2.*pio*bbb)
  px(1) = -1.*(px(4)+px(3)+px(2))
  py(1) = -1.*(py(4)+py(3)+py(2))
!     transverse mass for K-
  rmt3 = sqrt(dm3**2+px(1)**2+py(1)**2)
!     transverse mass for K+
  rmt4 = sqrt(dm4**2+px(2)**2+py(2)**2)
  If (restms<(rmt3+rmt4)) Goto 50
!        else: sampling success!
  pzcms = sqrt((restms**2-(rmt3+rmt4)**2)*(restms**2-(rmt3-rmt4)**2))/2./restms
  If (ranart(nseed)>0.5) Then
    pz(1) = pzcms
    pz(2) = -pzcms
  Else
    pz(1) = -pzcms
    pz(2) = pzcms
  End If
  beta = restpz/resten
  gama = 1./sqrt(1.-beta**2)
  pz(1) = pz(1)*gama + beta*gama*sqrt(rmt3**2+pz(1)**2)
  pz(2) = pz(2)*gama + beta*gama*sqrt(rmt4**2+pz(2)**2)
  pe(1) = sqrt(rmt3**2+pz(1)**2)
  pe(2) = sqrt(rmt4**2+pz(2)**2)

  iflag = 1
  Return
End Subroutine fstate

!-----------------------------------------------------------------------

!.....extracted from G. Song's ART expasion including K- interactions
!.....file `NPIK.FOR'

!***************************************
!        subroutine npik(irun,iseed,dt,nt,ictrl,i1,i2,srt,
!     &                  pcx,pcy,pcz,nchrg,ratiok)
Subroutine npik(irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg, ratiok, iblock)
!
! Process: PI + N -> K(-) + ANYTHING
! 1.  PI- + P -> P + K0 + K-
! 2.  PI+ + N -> P + K+ + K-
! 3.  PI0 + P -> P + K+ + K-
! 4.  PI0 + N -> P + K0 + K-
! 5.  PI0 + N -> N + K+ + K-
! 6.  PI- + P -> N + K+ + K-
! 7.  PI- + N -> N + K0 + K-
! NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
!***************************************
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (aka=0.498)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
  Common /nn/nnn
!c      SAVE /NN/
  Common /run/num
!c      SAVE /RUN/
  Common /pa/rpion(3, maxstr, maxr)
!c      SAVE /PA/
  Common /pb/ppion(3, maxstr, maxr)
!c      SAVE /PB/
  Common /pc/epion(maxstr, maxr)
!c      SAVE /PC/
  Common /pd/lpion(maxstr, maxr)
!c      SAVE /PD/
  Dimension bb(3), p1(4), p2(4), p3(4), px(4), py(4), pz(4)
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
  px1cm = pcx
  py1cm = pcy
  pz1cm = pcz
  ictrl = 1
  lb1 = lb(i1)
  lb2 = lb(i2)
  k1 = i1
  k2 = i2
!        k1 must be bayron. k2 be meson. If not, exchange.
  If (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13)) Then
    k1 = i2
    k2 = i1
  End If
!bz3/8/99 neutralk
!bz10/12/99
!        LB(I1) = 1 + 2 * RANART(NSEED)
!        LB(I2) = 23
  lb(k1) = 1 + int(2*ranart(nseed))
  lb(k2) = 23
!       pkmax=sqrt((srt**2-(aka+0.938+aka)**2)*(srt**2-(aka+0.938-aka)**2))
!     &           /2./srt
  pkmax = sqrt((srt**2-(aka+0.938+aka)**2)*(srt**2-(aka+0.938-aka)**2))/2./srt
  pk = ranart(nseed)*pkmax
!-----------------------------------------------------
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1.-css**2)
  fai = 2*3.1415926*ranart(nseed)
  p3(1) = pk*sss*cos(fai)
  p3(2) = pk*sss*sin(fai)
  p3(3) = pk*css
  eip = srt - sqrt(aka**2+pk**2)
  rmnp = sqrt(eip**2-pk**2)
  Do i = 1, 3
    bb(i) = -1.*p3(i)/eip
  End Do
!        bb: velocity of the other two particles as a whole.
  pznp = sqrt((rmnp**2-(aka+0.938)**2)*(rmnp**2-(0.938-aka)**2))/2./rmnp
!-----------------------------------------------------
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1.-css**2)
  fai = 2*3.1415926*ranart(nseed)
  p1(1) = pznp*sss*cos(fai)
  p1(2) = pznp*sss*sin(fai)
  p1(3) = pznp*css
  p1(4) = sqrt(0.938**2+pznp**2)
  p2(4) = sqrt(aka**2+pznp**2)
  Do i = 1, 3
    p2(i) = -1.*p1(i)
  End Do
!        p1,p2: the momenta of the two particles in their cms
!        p1: momentum of N or P
!        p2: momentum of anti_kaon
!        p3: momentum of K0 or K+
  ilo = 1
!        write(61,*)'--------p1,p2',p1,p2
!        write(61,*)'--------bb',bb
  Call lorntz(ilo, bb, p1, p2)
!******* Checking *************
!        pxsum = p1(1)+p2(1)+p3(1)
!        pysum = p1(2)+p2(2)+p3(2)
!        pzsum = p1(3)+p2(3)+p3(3)
!        pesum = p1(4)+p2(4)+sqrt(p3(1)**2+p3(2)**2+p3(3)**2+aka**2)-srt
!        write(61,*)'---p1,pxsum',p1,pxsum
!        write(61,*)'---p2,pysum',p2,pysum
!        write(61,*)'---p3,pzsum',p3,pzsum
!        write(61,*)'---pesum',pesum
!***********************************

! Rotate the momenta of particles in the cms of I1 & I2
! px(1), py(1), pz(1): momentum of I1
! px(2), py(2), pz(2): momentum of I2
! px(3), py(3), pz(3): momentum of anti-kaon

!     10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = p1(1)
  pyrota = p1(2)
  pzrota = p1(3)
!        call rotate(pcx,pcy,pcz,p1(1),p1(2),p1(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p1(1) = pxrota
  p1(2) = pyrota
  p1(3) = pzrota
!
  pxrota = p2(1)
  pyrota = p2(2)
  pzrota = p2(3)
!        call rotate(pcx,pcy,pcz,p2(1),p2(2),p2(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p2(1) = pxrota
  p2(2) = pyrota
  p2(3) = pzrota
!
  pxrota = p3(1)
  pyrota = p3(2)
  pzrota = p3(3)
!        call rotate(pcx,pcy,pcz,p3(1),p3(2),p3(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p3(1) = pxrota
  p3(2) = pyrota
  p3(3) = pzrota

  nnn = nnn + 1
!     K(-)
  lpion(nnn, irun) = 21
!     aka: rest mass of K
  epion(nnn, irun) = aka
! Find the momenta of particles in the final state in the nucleus_nucleus
! cms frame.   Lorentz transformation into lab frame.
  e1cm = sqrt(0.938**2+p1(1)**2+p1(2)**2+p1(3)**2)
  p1beta = p1(1)*betax + p1(2)*betay + p1(3)*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + p1(1)
  pt2i1 = betay*transf + p1(2)
  pt3i1 = betaz*transf + p1(3)
  eti1 = 0.938
  lb1 = lb(k1)

! For second nulceon, same
  e2cm = sqrt(aka**2+p3(1)**2+p3(2)**2+p3(3)**2)
  p2beta = p3(1)*betax + p3(2)*betay + p3(3)*betaz
  transf = gamma*(gamma*p2beta/(gamma+1)+e2cm)
  pt1i2 = betax*transf + p3(1)
  pt2i2 = betay*transf + p3(2)
  pt3i2 = betaz*transf + p3(3)
  eti2 = aka
  lb2 = lb(k2)

!        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
!       k1 stand for nucleon, k2 stand for kaon. lpion stand for Kbar.
  p(1, k1) = pt1i1
  p(2, k1) = pt2i1
  p(3, k1) = pt3i1
  e(k1) = eti1
  lb(k1) = lb1
  p(1, k2) = pt1i2
  p(2, k2) = pt2i2
  p(3, k2) = pt3i2
  e(k2) = eti2
  lb(k2) = lb2

!                px1 = p(1,i1)
!                py1 = p(2,i1)
!                pz1 = p(3,i1)
!                em1 = e(i1)
!                id(i1) = 2
!                id(i2) = 2
!                id1 = id(i1)
!     K(+)K(-) production
  iblock = 101
! Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
!  p2:  momentum of anti-kaon.
!        epcmk = sqrt(epion(nnn,irun)**2 + p2(1)**2 + p2(2)**2 + p2(3)**2)
  epcmk = sqrt(epion(nnn,irun)**2+p2(1)**2+p2(2)**2+p2(3)**2)
  betak = p2(1)*betax + p2(2)*betay + p2(3)*betaz
  transf = gamma*(gamma*betak/(gamma+1.)+epcmk)
  ppion(1, nnn, irun) = betax*transf + p2(1)
  ppion(2, nnn, irun) = betay*transf + p2(2)
  ppion(3, nnn, irun) = betaz*transf + p2(3)
!lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
!bz3/2/99
!        write(400,*)'2 ', ppion(1,nnn,irun), ppion(2,nnn,irun),
!     &                    ppion(3,nnn,irun), dt*nt, srt
!bz3/2/99end
!        write(420,*)ppion(1,nnn,irun), ppion(2,nnn,irun),
!     &                    ppion(3,nnn,irun), dt*nt, srt
  k = i2
  If (lb(i1)==1 .Or. lb(i1)==2) k = i1
  rpion(1, nnn, irun) = r(1, k)
  rpion(2, nnn, irun) = r(2, k)
  rpion(3, nnn, irun) = r(3, k)
  Return
End Subroutine npik

!-----------------------------------------------------------------------

!.....extracted from G. Song's ART expasion including K- interactions
!.....file `PIHYPN.FOR'

!*****************************************
Subroutine pihypn(ielstc, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, nchrg, iblock)
!
! Process: PI + sigma(or Lambda) -> Kbar + N
! NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
!*****************************************

! NOTE: for PI + Hyperon: the produced kaons have mass 0.498
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (aka=0.498)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
  Common /nn/nnn
!c      SAVE /NN/
  Common /run/num
!c      SAVE /RUN/
  Common /pa/rpion(3, maxstr, maxr)
!c      SAVE /PA/
  Common /pb/ppion(3, maxstr, maxr)
!c      SAVE /PB/
  Common /pc/epion(maxstr, maxr)
!c      SAVE /PC/
  Common /pd/lpion(maxstr, maxr)
!c      SAVE /PD/
  Dimension p1(4), p2(4)
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
  px1cm = pcx
  py1cm = pcy
  pz1cm = pcz
  ictrl = 1
!sp06/07/01
  If (ielstc==1) Then
!    L/Si + meson -> L/Si + meson
    k1 = i1
    k2 = i2
    dm3 = e(k1)
    dm4 = e(k2)
    iblock = 10
  Else
    iblock = 12
!sp06/07/01 end
!        PI + Sigma(or Lambda) -> Kbar + N
    k1 = i1
    k2 = i2
!        k1 must be bayron! So if I1 is PI, exchange k1 & k2.
    If (lb(i1)<14 .Or. lb(i1)>17) Then
      k1 = i2
      k2 = i1
    End If
!bz3/8/99 neutralk
    lb(k1) = 1 + int(2*ranart(nseed))
    If (nchrg==-2) lb(k1) = 6
!     if(nchrg.eq.-1) lb(k1)=2
!     if(nchrg.eq. 0) lb(k1)=1
!     if(nchrg.eq. 1) lb(k1)=9
    If (nchrg==2) lb(k1) = 9
!bz3/8/99 neutralk end

!     K-
    lb(k2) = 21
    dm3 = 0.938
    If (nchrg==-2 .Or. nchrg==1) dm3 = 1.232
    dm4 = aka
!        dm3,dm4: the mass of final state particles.
  End If

!*******Now, antikaon will be created.
!        call antikaon_fstate(iseed,srt,dm1,dm2,dm3,dm4,px,py,pz,icou1)
!        pkmax: the maximum momentum of anti-kaon
  pkmax = sqrt((srt**2-(dm3+dm4)**2)*(srt**2-(dm3-dm4)**2))/2./srt
  pk = pkmax
!-----------------------------------------------------
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1.-css**2)
  fai = 2*3.1415926*ranart(nseed)
  p1(1) = pk*sss*cos(fai)
  p1(2) = pk*sss*sin(fai)
  p1(3) = pk*css
  Do i = 1, 3
    p2(i) = -1.*p1(i)
  End Do
!        p1,p2: the momenta of the two particles in their cms
!        p1: momentum of kaon
!        p2: momentum of Kbar

! Rotate the momenta of particles in the cms of I1 & I2
!lin-10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = p1(1)
  pyrota = p1(2)
  pzrota = p1(3)
!        call rotate(pcx,pcy,pcz,p1(1),p1(2),p1(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p1(1) = pxrota
  p1(2) = pyrota
  p1(3) = pzrota
!
  pxrota = p2(1)
  pyrota = p2(2)
  pzrota = p2(3)
!        call rotate(pcx,pcy,pcz,p2(1),p2(2),p2(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p2(1) = pxrota
  p2(2) = pyrota
  p2(3) = pzrota
!lin-10/28/02-end

! Find the momenta of particles in the final state in the nucleus_nucleus
! cms frame.   Lorentz transformation into lab frame.
  e1cm = sqrt(dm3**2+p1(1)**2+p1(2)**2+p1(3)**2)
  p1beta = p1(1)*betax + p1(2)*betay + p1(3)*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + p1(1)
  pt2i1 = betay*transf + p1(2)
  pt3i1 = betaz*transf + p1(3)
  eti1 = dm3
  lb1 = lb(k1)

! For second kaon, same
  e2cm = sqrt(dm4**2+p2(1)**2+p2(2)**2+p2(3)**2)
  p2beta = p2(1)*betax + p2(2)*betay + p2(3)*betaz
  transf = gamma*(gamma*p2beta/(gamma+1)+e2cm)
  pt1i2 = betax*transf + p2(1)
  pt2i2 = betay*transf + p2(2)
  pt3i2 = betaz*transf + p2(3)
!bz3/2/99
!        write(400,*)'3 ', pt1i2, pt2i2, pt3i2, dt*nt, srt
!bz3/2/99end
!        write(430,*)pt1i2, pt2i2, pt3i2, dt*nt, srt
  eti2 = dm4
  lb2 = lb(k2)

!        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
!        k1=i1
!        k2=i2
!       k1 stand for nucleon, k2 stand for kaon.
  p(1, k1) = pt1i1
  p(2, k1) = pt2i1
  p(3, k1) = pt3i1
  e(k1) = eti1
  lb(k1) = lb1
  p(1, k2) = pt1i2
  p(2, k2) = pt2i2
  p(3, k2) = pt3i2
  e(k2) = eti2
  lb(k2) = lb2

!c                iblock = 101  ! K(+)K(-) production
! Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
  Return
End Subroutine pihypn

!-----------------------------------------------------------------------

!.....extracted from G. Song's ART expasion including K- interactions
!.....file `KAONN.FOR'

!***************************************
Subroutine kaonn(brel, brsgm, irun, iseed, dt, nt, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, nchrg)
!
! Process: PI + sigma(or Lambda) <- Kbar + N
! NOTE: the mass of K is assumed to be same as K0. ie. 0.498 NOT 0.494
!***************************************
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
  Common /nn/nnn
!c      SAVE /NN/
  Common /run/num
!c      SAVE /RUN/
  Common /pa/rpion(3, maxstr, maxr)
!c      SAVE /PA/
  Common /pb/ppion(3, maxstr, maxr)
!c      SAVE /PB/
  Common /pc/epion(maxstr, maxr)
!c      SAVE /PC/
  Common /pd/lpion(maxstr, maxr)
!c      SAVE /PD/
  Dimension p1(4), p2(4)
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
  px1cm = pcx
  py1cm = pcy
  pz1cm = pcz
  ictrl = 1
!        ratio: used for isospin decision.
  k1 = i1
  k2 = i2
!        k1 must be bayron! So if I1 is Kaon, exchange k1 & k2.
  If (e(i1)<0.5 .And. e(i1)>0.01) Then
    k1 = i2
    k2 = i1
  End If
!** note: for print out only *******************************
!     record kaon's mass
  eee = e(k2)
!** end **************
  rrr = ranart(nseed)
  If (rrr<brel) Then
!       Kbar + N -> Kbar + N
    lb1 = lb(k1)
    lb2 = lb(k2)
    em1 = e(k1)
    em2 = e(k2)
    iblock = 10
  Else
    iblock = 12
    If (rrr<(brel+brsgm)) Then
!        nchrg: Net charges of the two incoming particles.
!           Kbar + N -> Sigma + PI
      em1 = asa
      em2 = 0.138

!bz3/8/99 neutralk
      lb1 = 15 + int(3*ranart(nseed))
      lb2 = 3 + int(3*ranart(nseed))
    Else
!           Kbar + N -> Lambda + PI
      em1 = ala
      em2 = 0.138
!     LAmbda
      lb1 = 14
!bz3/8/99 neutralk
      lb2 = 3 + int(3*ranart(nseed))
!           if(nchrg.eq.1)  lb2=5  ! K- + D++ -> Lambda + PI+
!           if(nchrg.eq.0)  lb2=4  ! K- + p(D+,N*+) -> Lambda + PI0
!          if(nchrg.eq.-1) lb2=3 ! K- + n(D,N*) -> Lambda + PI-
!bz3/8/99 neutralk

    End If
  End If
  lb(k1) = lb1
  lb(k2) = lb2

!*******Now, antikaon will be created.
!        call antikaon_fstate(iseed,srt,dm1,dm2,dm3,dm4,px,py,pz,icou1)
!        pkmax: the maximum momentum of anti-kaon
!        write(63,*)'srt,em1,em2',srt,em1,em2
!        write(63,*)'-srt,em1,em2',srt,em1,em2
  pkmax = sqrt((srt**2-(em1+em2)**2)*(srt**2-(em1-em2)**2))/2./srt
  pk = pkmax
!-----------------------------------------------------
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1.-css**2)
  fai = 2*3.1415926*ranart(nseed)
  p1(1) = pk*sss*cos(fai)
  p1(2) = pk*sss*sin(fai)
  p1(3) = pk*css
  Do i = 1, 3
    p2(i) = -1.*p1(i)
  End Do
!        p1,p2: the momenta of the two particles in their cms
!        p1: momentum of kaon
!        p2: momentum of Kbar

! Rotate the momenta of particles in the cms of I1 & I2

!lin-10/28/02 get rid of argument usage mismatch in rotate():
  pxrota = p1(1)
  pyrota = p1(2)
  pzrota = p1(3)
!        call rotate(pcx,pcy,pcz,p1(1),p1(2),p1(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p1(1) = pxrota
  p1(2) = pyrota
  p1(3) = pzrota
!
  pxrota = p2(1)
  pyrota = p2(2)
  pzrota = p2(3)
!        call rotate(pcx,pcy,pcz,p2(1),p2(2),p2(3))
  Call rotate(pcx, pcy, pcz, pxrota, pyrota, pzrota)
  p2(1) = pxrota
  p2(2) = pyrota
  p2(3) = pzrota
!lin-10/28/02-end

! Find the momenta of particles in the final state in the nucleus_nucleus
! cms frame.   Lorentz transformation into lab frame.
  e1cm = sqrt(em1**2+p1(1)**2+p1(2)**2+p1(3)**2)
  p1beta = p1(1)*betax + p1(2)*betay + p1(3)*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + p1(1)
  pt2i1 = betay*transf + p1(2)
  pt3i1 = betaz*transf + p1(3)
  eti1 = em1

! For second kaon, same
  e2cm = sqrt(em2**2+p2(1)**2+p2(2)**2+p2(3)**2)
  p2beta = p2(1)*betax + p2(2)*betay + p2(3)*betaz
  transf = gamma*(gamma*p2beta/(gamma+1)+e2cm)
  pt1i2 = betax*transf + p2(1)
  pt2i2 = betay*transf + p2(2)
  pt3i2 = betaz*transf + p2(3)
  eti2 = em2

!        if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0.)then
!        k1=i1
!        k2=i2
!       k1 stand for bayron, k2 stand for meson.
  p(1, k1) = pt1i1
  p(2, k1) = pt2i1
  p(3, k1) = pt3i1
  e(k1) = eti1
  p(1, k2) = pt1i2
  p(2, k2) = pt2i2
  p(3, k2) = pt3i2
  e(k2) = eti2

!c                iblock = 101  ! K(+)K(-) production
! Get Kaons' momenta and coordinates in nucleus-nucleus cms. frame.
  Return
End Subroutine kaonn

!=======================================================================

!lin Below is the previous artana.f:
!=======================================================================

!.....analysis subroutine before the hadronic space-time evolution

Subroutine artan1
  Parameter (maxstr=150001, maxr=1)
!.....y cut for mt spectrum
!bz3/17/99
!      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
  Parameter (ymt1=-1.0, ymt2=1.0)
!bz3/17/99 end
!.....bin width for mt spectrum and y spectrum
!lin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
!      PARAMETER (BMT = 0.05, BY = 0.2)
  Parameter (bmt=0.05, by=0.4)
  Common /run/num
!c      SAVE /RUN/
  Common /arerc1/multi1(maxr)
!c      SAVE /ARERC1/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
!bz3/17/99
!     &     dm1k0s(50), DMT1LA(50), DMT1LB(50)
!c      SAVE /ARPRC1/
  Common /arana1/dy1ntb(50), dy1ntp(50), dy1hm(50), dy1kp(50), dy1km(50), dy1k0s(50), dy1la(50), dy1lb(50), dy1phi(50), dm1pip(50), dm1pim(50), dmt1pr(50), dmt1pb(50), dmt1kp(50), dm1km(50), dm1k0s(50), dmt1la(50), dmt1lb(50), dy1msn(50), dy1pip(50), dy1pim(50), dy1pi0(50), dy1pr(50), dy1pb(50), dy1neg(50), dy1ch(50), de1neg(50), de1ch(50)
!c      SAVE /ARANA1/
  Save

!bz3/17/99 end
  Do j = 1, num
    Do i = 1, multi1(j)
      ityp = ityp1(i, j)
      px = px1(i, j)
      py = py1(i, j)
      pz = pz1(i, j)
      ee = ee1(i, j)
      xm = xm1(i, j)
!     2/24/03 leptons and photons:
      If (xm<0.01) Goto 200
      ptot = sqrt(px**2+py**2+pz**2)
!lin-9/2012 determine pseudo-rapidity more generally:
!            eta = 0.5*alog((Ptot+pz+1e-5)/(ptot-pz+1e-5))
      If ((px**2+py**2)>0.) Then
        eta = asinh(pz/sqrt(px**2+py**2))
      Else
        eta = 1000000.0*sign(1., pz)
!lin-2/2013 for spectator target nucleons in LAB frame,
!     note that finite precision of HBOOST
!     would give spectator target nucleons a small but non-zero pz:
        If (abs(pz)<=1E-3) eta = 0.
      End If

      xmt = sqrt(px**2+py**2+xm**2)
      dxmt = xmt - xm
!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZ) .GE. EE) THEN
!               PRINT *, 'IN ARTAN1'
!               PRINT *, 'PARTICLE ', I, ' RUN ', J, 'PREC ERR'
!cbzdbg2/16/99
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!cbzdbg2/16/99
!cbzdbg2/15/99
!               PRINT *, ' PZ = ', PZ, ' EE = ', EE
!cbzdbg2/16/99
!               PRINT *, ' XM = ', XM
!cbzdbg2/16/99end
!               GOTO 200
!            else
!c            Y = 0.5 * LOG((EE + PZ) / (EE - PZ))
!               Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ +1e-5))
!c               STOP
!cbzdbg2/15/99end
!            END IF
      If (xmt>0.) Then
        y = asinh(pz/xmt)
      Else
        Print *, ' IN ARTAN1 mt=0'
        y = 1000000.0*sign(1., pz)
      End If

!.....rapidity cut for the rapidity distribution
      If (abs(y)>=10.0) Goto 100
!lin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
!            IY = 1 + int(ABS(Y) / BY)
!            Ieta = 1 + int(ABS(eta) / BY)
      If (abs(eta)>=10.0) Goto 100
      iy = 1 + int((y+10.)/by)
      ieta = 1 + int((eta+10.)/by)

      If (ityp<-1000) Then
        dy1ntb(iy) = dy1ntb(iy) - 1.0
      End If
      If (ityp>1000) Then
        dy1ntb(iy) = dy1ntb(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy1ntp(iy) = dy1ntp(iy) - 1.0
      End If
      If (ityp==2212) Then
        dy1ntp(iy) = dy1ntp(iy) + 1.0
      End If
!            IF (ITYP .EQ. -211 .OR. ITYP .EQ. -321 .OR.
!     &         ITYP .EQ. -2212) THEN
      If (ityp==-2112) Then
        dy1hm(iy) = dy1hm(iy) + 1.0
      End If
!
      If (luchge(ityp)/=0) Then
        dy1ch(iy) = dy1ch(iy) + 1.0
        de1ch(ieta) = de1ch(ieta) + 1.0
        If (luchge(ityp)<0) Then
          dy1neg(iy) = dy1neg(iy) + 1.0
          de1neg(ieta) = de1neg(ieta) + 1.0
        End If
      End If

!bz3/17/99
      If ((ityp>=100 .And. ityp<1000) .Or. (ityp>-1000 .And. ityp<=-100)) Then
        dy1msn(iy) = dy1msn(iy) + 1.0
      End If
      If (ityp==211) Then
        dy1pip(iy) = dy1pip(iy) + 1.0
      End If
      If (ityp==-211) Then
        dy1pim(iy) = dy1pim(iy) + 1.0
      End If
      If (ityp==111) Then
        dy1pi0(iy) = dy1pi0(iy) + 1.0
      End If
      If (ityp==2212) Then
        dy1pr(iy) = dy1pr(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy1pb(iy) = dy1pb(iy) + 1.0
      End If
!bz3/17/99 end
      If (ityp==321) Then
        dy1kp(iy) = dy1kp(iy) + 1.0
      End If
      If (ityp==-321) Then
        dy1km(iy) = dy1km(iy) + 1.0
      End If
!lin-4/24/03 evaluate K0L instead of K0S, since sometimes we may decay K0S:
!            IF (ITYP .EQ. 310) THEN
      If (ityp==130) Then
        dy1k0s(iy) = dy1k0s(iy) + 1.0
      End If
      If (ityp==3122) Then
        dy1la(iy) = dy1la(iy) + 1.0
      End If
      If (ityp==-3122) Then
        dy1lb(iy) = dy1lb(iy) + 1.0
      End If
      If (ityp==333) Then
        dy1phi(iy) = dy1phi(iy) + 1.0
      End If

!.....insert rapidity cut for mt spectrum here
      100 If (y<ymt1 .Or. y>ymt2) Goto 200
      If (dxmt>=50.0*bmt .Or. dxmt==0) Goto 200
      imt = 1 + int(dxmt/bmt)
      If (ityp==211) Then
        dm1pip(imt) = dm1pip(imt) + 1.0/xmt
      End If
      If (ityp==-211) Then
        dm1pim(imt) = dm1pim(imt) + 1.0/xmt
      End If
      If (ityp==2212) Then
        dmt1pr(imt) = dmt1pr(imt) + 1.0/xmt
      End If
      If (ityp==-2212) Then
        dmt1pb(imt) = dmt1pb(imt) + 1.0/xmt
      End If
      If (ityp==321) Then
        dmt1kp(imt) = dmt1kp(imt) + 1.0/xmt
      End If
      If (ityp==-321) Then
        dm1km(imt) = dm1km(imt) + 1.0/xmt
      End If
!lin-4/24/03:
!            IF (ITYP .EQ. 310) THEN
      If (ityp==130) Then
        dm1k0s(imt) = dm1k0s(imt) + 1.0/xmt
      End If
      If (ityp==3122) Then
        dmt1la(imt) = dmt1la(imt) + 1.0/xmt
      End If
      If (ityp==-3122) Then
        dmt1lb(imt) = dmt1lb(imt) + 1.0/xmt
      End If

      200 Continue
    End Do
  End Do

  Return
End Subroutine artan1

!-----------------------------------------------------------------------

!.....analysis subroutine after the hadronic space-time evolution

Subroutine artan2

  Parameter (maxstr=150001, maxr=1)
!.....y cut for mt spectrum
!bz3/17/99
!      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
  Parameter (ymt1=-1.0, ymt2=1.0)
!bz3/17/99 end
!.....bin width for mt spectrum and y spectrum
!      PARAMETER (BMT = 0.05, BY = 0.2)
  Parameter (bmt=0.05, by=0.4)
  Common /run/num
!c      SAVE /RUN/
  Common /arerc1/multi1(maxr)
!c      SAVE /ARERC1/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
!bz3/17/99
!     &     dm2k0s(50), DMT2LA(50), DMT2LB(50)
!c      SAVE /ARPRC1/
  Common /arana2/dy2ntb(50), dy2ntp(50), dy2hm(50), dy2kp(50), dy2km(50), dy2k0s(50), dy2la(50), dy2lb(50), dy2phi(50), dm2pip(50), dm2pim(50), dmt2pr(50), dmt2pb(50), dmt2kp(50), dm2km(50), dm2k0s(50), dmt2la(50), dmt2lb(50), dy2msn(50), dy2pip(50), dy2pim(50), dy2pi0(50), dy2pr(50), dy2pb(50), dy2neg(50), dy2ch(50), de2neg(50), de2ch(50)
!bz3/17/99 end
!c      SAVE /ARANA2/
  Save

  Do j = 1, num
    Do i = 1, multi1(j)
      ityp = ityp1(i, j)
      px = px1(i, j)
      py = py1(i, j)
      pz = pz1(i, j)
      ee = ee1(i, j)
      xm = xm1(i, j)
      xmt = sqrt(px**2+py**2+xm**2)
!     2/24/03 leptons and photons:
      If (xm<0.01) Goto 200
      dxmt = xmt - xm
      ptot = sqrt(px**2+py**2+pz**2)
!lin-9/2012 determine pseudo-rapidity more generally:
!            eta = 0.5*alog((Ptot+pz+1e-5)/(ptot-pz+1e-5))
      If ((px**2+py**2)>0.) Then
        eta = asinh(pz/sqrt(px**2+py**2))
      Else
        eta = 1000000.0*sign(1., pz)
!lin-2/2013 for spectator target nucleons in LAB frame:
        If (abs(pz)<=1E-3) eta = 0.
      End If

!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZ) .GE. EE) THEN
!               PRINT *, 'IN ARTAN2'
!               PRINT *, 'PARTICLE ', I, ' RUN ', J, 'PREC ERR'
!cbzdbg2/16/99
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!cbzdbg2/16/99
!cbzdbg2/15/99
!               PRINT *, ' PZ = ', PZ, ' EE = ', EE
!cbzdbg2/16/99
!               PRINT *, ' XM = ', XM
!cbzdbg2/16/99end
!               GOTO 200
!c               STOP
!cbzdbg2/15/99end
!            END IF
!            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
      If (xmt>0.) Then
        y = asinh(pz/xmt)
      Else
        Print *, ' IN ARTAN2 mt=0'
        y = 1000000.0*sign(1., pz)
      End If

!.....rapidity cut for the rapidity distribution
      If (abs(y)>=10.0) Goto 100
!            IY = 1 + int(ABS(Y) / BY)
!            Ieta = 1 + int(ABS(eta) / BY)
      If (abs(eta)>=10.0) Goto 100
      iy = 1 + int((y+10.)/by)
      ieta = 1 + int((eta+10.)/by)

      If (ityp<-1000) Then
        dy2ntb(iy) = dy2ntb(iy) - 1.0
      End If
      If (ityp>1000) Then
        dy2ntb(iy) = dy2ntb(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy2ntp(iy) = dy2ntp(iy) - 1.0
      End If
      If (ityp==2212) Then
        dy2ntp(iy) = dy2ntp(iy) + 1.0
      End If
      If (ityp==-2112) Then
        dy2hm(iy) = dy2hm(iy) + 1.0
      End If

      If (luchge(ityp)/=0) Then
        dy2ch(iy) = dy2ch(iy) + 1.0
        de2ch(ieta) = de2ch(ieta) + 1.0
        If (luchge(ityp)<0) Then
          dy2neg(iy) = dy2neg(iy) + 1.0
          de2neg(ieta) = de2neg(ieta) + 1.0
        End If
      End If

!bz3/17/99
      If ((ityp>=100 .And. ityp<1000) .Or. (ityp>-1000 .And. ityp<=-100)) Then
        dy2msn(iy) = dy2msn(iy) + 1.0
      End If
      If (ityp==211) Then
        dy2pip(iy) = dy2pip(iy) + 1.0
      End If
      If (ityp==-211) Then
        dy2pim(iy) = dy2pim(iy) + 1.0
      End If
      If (ityp==111) Then
        dy2pi0(iy) = dy2pi0(iy) + 1.0
      End If
      If (ityp==2212) Then
        dy2pr(iy) = dy2pr(iy) + 1.0
      End If
      If (ityp==-2212) Then
        dy2pb(iy) = dy2pb(iy) + 1.0
      End If
!bz3/17/99 end
      If (ityp==321) Then
        dy2kp(iy) = dy2kp(iy) + 1.0
      End If
      If (ityp==-321) Then
        dy2km(iy) = dy2km(iy) + 1.0
      End If
!lin-4/24/03:
!            IF (ITYP .EQ. 310) THEN
      If (ityp==130) Then
        dy2k0s(iy) = dy2k0s(iy) + 1.0
      End If
      If (ityp==3122) Then
        dy2la(iy) = dy2la(iy) + 1.0
      End If
      If (ityp==-3122) Then
        dy2lb(iy) = dy2lb(iy) + 1.0
      End If
      If (ityp==333) Then
        dy2phi(iy) = dy2phi(iy) + 1.0
      End If

!.....insert rapidity cut for mt spectrum here
      100 If (y<ymt1 .Or. y>ymt2) Goto 200
      If (dxmt>=50.0*bmt .Or. dxmt==0) Goto 200
      imt = 1 + int(dxmt/bmt)
      If (ityp==211) Then
        dm2pip(imt) = dm2pip(imt) + 1.0/xmt
      End If
      If (ityp==-211) Then
        dm2pim(imt) = dm2pim(imt) + 1.0/xmt
      End If
      If (ityp==2212) Then
        dmt2pr(imt) = dmt2pr(imt) + 1.0/xmt
      End If
      If (ityp==-2212) Then
        dmt2pb(imt) = dmt2pb(imt) + 1.0/xmt
      End If
      If (ityp==321) Then
        dmt2kp(imt) = dmt2kp(imt) + 1.0/xmt
      End If
      If (ityp==-321) Then
        dm2km(imt) = dm2km(imt) + 1.0/xmt
      End If
!lin-4/24/03:
!            IF (ITYP .EQ. 310) THEN
      If (ityp==130) Then
        dm2k0s(imt) = dm2k0s(imt) + 1.0/xmt
      End If
      If (ityp==3122) Then
        dmt2la(imt) = dmt2la(imt) + 1.0/xmt
      End If
      If (ityp==-3122) Then
        dmt2lb(imt) = dmt2lb(imt) + 1.0/xmt
      End If

      200 Continue
    End Do
  End Do

  Return
End Subroutine artan2

!-----------------------------------------------------------------------

!.....output analysis results at the end of the simulation

Subroutine artout(nevnt)

  Parameter (maxstr=150001, maxr=1)
!.....y cut for mt spectrum
!bz3/17/99
!      PARAMETER (YMT1 = -0.4, YMT2 = 0.4)
  Parameter (ymt1=-1.0, ymt2=1.0)
!bz3/17/99 end
!.....bin width for mt spectrum and y spectrum
!      PARAMETER (BMT = 0.05, BY = 0.2)
  Parameter (bmt=0.05, by=0.4)
  Common /run/num
!c      SAVE /RUN/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
!bz3/17/99
!     &     dm1k0s(50), DMT1LA(50), DMT1LB(50)
!c      SAVE /ARPRC1/
  Common /arana1/dy1ntb(50), dy1ntp(50), dy1hm(50), dy1kp(50), dy1km(50), dy1k0s(50), dy1la(50), dy1lb(50), dy1phi(50), dm1pip(50), dm1pim(50), dmt1pr(50), dmt1pb(50), dmt1kp(50), dm1km(50), dm1k0s(50), dmt1la(50), dmt1lb(50), dy1msn(50), dy1pip(50), dy1pim(50), dy1pi0(50), dy1pr(50), dy1pb(50), dy1neg(50), dy1ch(50), de1neg(50), de1ch(50)
!bz3/17/99 end
!c      SAVE /ARANA1/
!bz3/17/99
!     &     dm2k0s(50), DMT2LA(50), DMT2LB(50)
  Common /arana2/dy2ntb(50), dy2ntp(50), dy2hm(50), dy2kp(50), dy2km(50), dy2k0s(50), dy2la(50), dy2lb(50), dy2phi(50), dm2pip(50), dm2pim(50), dmt2pr(50), dmt2pb(50), dmt2kp(50), dm2km(50), dm2k0s(50), dmt2la(50), dmt2lb(50), dy2msn(50), dy2pip(50), dy2pim(50), dy2pi0(50), dy2pr(50), dy2pb(50), dy2neg(50), dy2ch(50), de2neg(50), de2ch(50)
!c      SAVE /ARANA2/
  Save
!bz3/17/99 end
  Open (30, File='ana/dndy_netb.dat', Status='UNKNOWN')
  Open (31, File='ana/dndy_netp.dat', Status='UNKNOWN')
  Open (32, File='ana/dndy_nb.dat', Status='UNKNOWN')
  Open (33, File='ana/dndy_neg.dat', Status='UNKNOWN')
  Open (34, File='ana/dndy_ch.dat', Status='UNKNOWN')
  Open (35, File='ana/dnde_neg.dat', Status='UNKNOWN')
  Open (36, File='ana/dnde_ch.dat', Status='UNKNOWN')
  Open (37, File='ana/dndy_kp.dat', Status='UNKNOWN')
  Open (38, File='ana/dndy_km.dat', Status='UNKNOWN')
!lin-4/24/03
!      OPEN (39, FILE = 'ana/dndy_k0s.dat', STATUS = 'UNKNOWN')
  Open (39, File='ana/dndy_k0l.dat', Status='UNKNOWN')
  Open (40, File='ana/dndy_la.dat', Status='UNKNOWN')
  Open (41, File='ana/dndy_lb.dat', Status='UNKNOWN')
  Open (42, File='ana/dndy_phi.dat', Status='UNKNOWN')
!bz3/17/99
  Open (43, File='ana/dndy_meson.dat', Status='UNKNOWN')
  Open (44, File='ana/dndy_pip.dat', Status='UNKNOWN')
  Open (45, File='ana/dndy_pim.dat', Status='UNKNOWN')
  Open (46, File='ana/dndy_pi0.dat', Status='UNKNOWN')
  Open (47, File='ana/dndy_pr.dat', Status='UNKNOWN')
  Open (48, File='ana/dndy_pb.dat', Status='UNKNOWN')
!bz3/17/99 end

  Open (50, File='ana/dndmtdy_pip.dat', Status='UNKNOWN')
  Open (51, File='ana/dndmtdy_0_1_pim.dat', Status='UNKNOWN')
  Open (52, File='ana/dndmtdy_pr.dat', Status='UNKNOWN')
  Open (53, File='ana/dndmtdy_pb.dat', Status='UNKNOWN')
  Open (54, File='ana/dndmtdy_kp.dat', Status='UNKNOWN')
  Open (55, File='ana/dndmtdy_0_5_km.dat', Status='UNKNOWN')
  Open (56, File='ana/dndmtdy_k0s.dat', Status='UNKNOWN')
  Open (57, File='ana/dndmtdy_la.dat', Status='UNKNOWN')
  Open (58, File='ana/dndmtdy_lb.dat', Status='UNKNOWN')
!lin-9/26/03 no symmetrization in y (or eta) for ana/*.dat:
!      SCALE1 = 1. / REAL(NEVNT * NUM) / BY / 2.0
  scale1 = 1./real(nevnt*num)/by
  scale2 = 1./real(nevnt*num)/bmt/(ymt2-ymt1)
!
  Do i = 1, 50
    ymid = -10. + by*(i-0.5)
    Write (30, 333) ymid, scale1*dy1ntb(i)
    Write (31, 333) ymid, scale1*dy1ntp(i)
    Write (32, 333) ymid, scale1*dy1hm(i)
    Write (37, 333) ymid, scale1*dy1kp(i)
    Write (38, 333) ymid, scale1*dy1km(i)
    Write (39, 333) ymid, scale1*dy1k0s(i)
    Write (40, 333) ymid, scale1*dy1la(i)
    Write (41, 333) ymid, scale1*dy1lb(i)
    Write (42, 333) ymid, scale1*dy1phi(i)
    Write (33, 333) ymid, scale1*dy1neg(i)
    Write (34, 333) ymid, scale1*dy1ch(i)
    Write (35, 333) ymid, scale1*de1neg(i)
    Write (36, 333) ymid, scale1*de1ch(i)
    Write (43, 333) ymid, scale1*dy1msn(i)
    Write (44, 333) ymid, scale1*dy1pip(i)
    Write (45, 333) ymid, scale1*dy1pim(i)
    Write (46, 333) ymid, scale1*dy1pi0(i)
    Write (47, 333) ymid, scale1*dy1pr(i)
    Write (48, 333) ymid, scale1*dy1pb(i)

    If (dm1pip(i)/=0.0) Then
      Write (50, 333) bmt*(i-0.5), scale2*dm1pip(i)
    End If
    If (dm1pim(i)/=0.0) Then
      Write (51, 333) bmt*(i-0.5), scale2*0.1*dm1pim(i)
    End If
    If (dmt1pr(i)/=0.0) Then
      Write (52, 333) bmt*(i-0.5), scale2*dmt1pr(i)
    End If
    If (dmt1pb(i)/=0.0) Then
      Write (53, 333) bmt*(i-0.5), scale2*dmt1pb(i)
    End If
    If (dmt1kp(i)/=0.0) Then
      Write (54, 333) bmt*(i-0.5), scale2*dmt1kp(i)
    End If
    If (dm1km(i)/=0.0) Then
      Write (55, 333) bmt*(i-0.5), scale2*0.5*dm1km(i)
    End If
    If (dm1k0s(i)/=0.0) Then
      Write (56, 333) bmt*(i-0.5), scale2*dm1k0s(i)
    End If
    If (dmt1la(i)/=0.0) Then
      Write (57, 333) bmt*(i-0.5), scale2*dmt1la(i)
    End If
    If (dmt1lb(i)/=0.0) Then
      Write (58, 333) bmt*(i-0.5), scale2*dmt1lb(i)
    End If
  End Do
!
  Do i = 30, 48
    Write (i, *) 'after hadron evolution'
  End Do
  Do i = 50, 58
    Write (i, *) 'after hadron evolution'
  End Do

  Do i = 1, 50
    ymid = -10. + by*(i-0.5)
    Write (30, 333) ymid, scale1*dy2ntb(i)
    Write (31, 333) ymid, scale1*dy2ntp(i)
    Write (32, 333) ymid, scale1*dy2hm(i)
    Write (37, 333) ymid, scale1*dy2kp(i)
    Write (38, 333) ymid, scale1*dy2km(i)
    Write (39, 333) ymid, scale1*dy2k0s(i)
    Write (40, 333) ymid, scale1*dy2la(i)
    Write (41, 333) ymid, scale1*dy2lb(i)
    Write (42, 333) ymid, scale1*dy2phi(i)
    Write (33, 333) ymid, scale1*dy2neg(i)
    Write (34, 333) ymid, scale1*dy2ch(i)
    Write (35, 333) ymid, scale1*de2neg(i)
    Write (36, 333) ymid, scale1*de2ch(i)
    Write (43, 333) ymid, scale1*dy2msn(i)
    Write (44, 333) ymid, scale1*dy2pip(i)
    Write (45, 333) ymid, scale1*dy2pim(i)
    Write (46, 333) ymid, scale1*dy2pi0(i)
    Write (47, 333) ymid, scale1*dy2pr(i)
    Write (48, 333) ymid, scale1*dy2pb(i)
!
    If (dm2pip(i)/=0.0) Then
      Write (50, 333) bmt*(i-0.5), scale2*dm2pip(i)
    End If
    If (dm2pim(i)/=0.0) Then
      Write (51, 333) bmt*(i-0.5), scale2*0.1*dm2pim(i)
    End If
    If (dmt2pr(i)/=0.0) Then
      Write (52, 333) bmt*(i-0.5), scale2*dmt2pr(i)
    End If
    If (dmt2pb(i)/=0.0) Then
      Write (53, 333) bmt*(i-0.5), scale2*dmt2pb(i)
    End If
    If (dmt2kp(i)/=0.0) Then
      Write (54, 333) bmt*(i-0.5), scale2*dmt2kp(i)
    End If
    If (dm2km(i)/=0.0) Then
      Write (55, 333) bmt*(i-0.5), scale2*0.5*dm2km(i)
    End If
    If (dm2k0s(i)/=0.0) Then
      Write (56, 333) bmt*(i-0.5), scale2*dm2k0s(i)
    End If
    If (dmt2la(i)/=0.0) Then
      Write (57, 333) bmt*(i-0.5), scale2*dmt2la(i)
    End If
    If (dmt2lb(i)/=0.0) Then
      Write (58, 333) bmt*(i-0.5), scale2*dmt2lb(i)
    End If
  End Do

  Return
  333 Format (2(F12.5,1X))
End Subroutine artout

!-----------------------------------------------------------------------

!.....analysis subroutine in HIJING before parton cascade evolution
Subroutine hjana1

  Parameter (ymax=1.0, ymin=-1.0)
  Parameter (dmt=0.05, dy=0.2)
  Parameter (dr=0.2)
  Parameter (maxptn=400001, maxstr=150001)
  Dimension dyp1(50), dmyp1(200), deyp1(50)
  Dimension dyg1(50), dmyg1(200), deyg1(50)
  Dimension snyp1(50), smyp1(200), seyp1(50)
  Dimension snyg1(50), smyg1(200), seyg1(50)
  Dimension dnrpj1(50), dnrtg1(50), dnrin1(50), dnrtt1(50)
  Dimension dyg1c(50), dmyg1c(50), deyg1c(50)
  Dimension snrpj1(50), snrtg1(50), snrin1(50), snrtt1(50)
  Dimension snyg1c(50), smyg1c(50), seyg1c(50)
  Double Precision gx0, gy0, gz0, ft0, px0, py0, pz0, e0, xmass0

  Common /para1/mul
!c      SAVE /PARA1/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
!c      SAVE /prec1/
  Common /arevt/iaevt, iarun, miss
!c      SAVE /AREVT/
  Common /arout/iout
!c      SAVE /AROUT/
  Save
  Data iw/0/

  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 200
      dmyp1(i) = smyp1(i)
      dmyg1(i) = smyg1(i)
    End Do

    Do i = 1, 50
      dyp1(i) = snyp1(i)
      deyp1(i) = seyp1(i)
      dyg1(i) = snyg1(i)
      deyg1(i) = seyg1(i)
      dnrpj1(i) = snrpj1(i)
      dnrtg1(i) = snrtg1(i)
      dnrin1(i) = snrin1(i)
      dnrtt1(i) = snrtt1(i)
      dyg1c(i) = snyg1c(i)
      dmyg1c(i) = smyg1c(i)
      deyg1c(i) = seyg1c(i)
    End Do
    nsubp = nsubps
    nsubg = nsubgs
    nisg = nisgs
  Else
    Do i = 1, 200
      smyp1(i) = dmyp1(i)
      smyg1(i) = dmyg1(i)
    End Do

    Do i = 1, 50
      snyp1(i) = dyp1(i)
      seyp1(i) = deyp1(i)
      snyg1(i) = dyg1(i)
      seyg1(i) = deyg1(i)
      snrpj1(i) = dnrpj1(i)
      snrtg1(i) = dnrtg1(i)
      snrin1(i) = dnrin1(i)
      snrtt1(i) = dnrtt1(i)
      snyg1c(i) = dyg1c(i)
      smyg1c(i) = dmyg1c(i)
      seyg1c(i) = deyg1c(i)
    End Do
    nsubps = nsubp
    nsubgs = nsubg
    nisgs = nisg
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
!.....analysis
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      ityp = kfpj(i, j)
      px = pjpx(i, j)
      py = pjpy(i, j)
      pz = pjpz(i, j)
      pe = pjpe(i, j)
      pm = pjpm(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZ) .GE. PE) THEN
!               PRINT *, ' IN HJANA1, PROJ STR ', I, ' PART ', J
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', PE
!               PRINT *, ' XM = ', PM
!               GOTO 200
!            END IF
!            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        rap = 1000000.0*sign(1., pz)
      End If

      iy = 1 + int(abs(rap)/dy)
!lin-8/2014 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 100
      If (iy<1 .Or. iy>50) Goto 100
      dyp1(iy) = dyp1(iy) + 1.0
      deyp1(iy) = deyp1(iy) + xmt
      If (ityp==21) Then
        dyg1(iy) = dyg1(iy) + 1.0
        deyg1(iy) = deyg1(iy) + xmt
      End If
      100 Continue
      imt = 1 + int(dxmt/dmt)
      If (rap>ymax .Or. rap<=ymin) Goto 200
      If (imt>200) Goto 200
      dmyp1(imt) = dmyp1(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg1(imt) = dmyg1(imt) + 1.0/xmt
      End If
      200 Continue
    End Do
  End Do

  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      ityp = kftj(i, j)
      px = pjtx(i, j)
      py = pjty(i, j)
      pz = pjtz(i, j)
      pe = pjte(i, j)
      pm = pjtm(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZ) .GE. PE) THEN
!               PRINT *, ' IN HJANA1, TARG STR ', I, ' PART ', J
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', PE
!               PRINT *, ' XM = ', PM
!               GOTO 400
!            END IF
!            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA1 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If

      iy = 1 + int(abs(rap)/dy)
!lin-8/2014 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 300
      If (iy<1 .Or. iy>50) Goto 300
      dyp1(iy) = dyp1(iy) + 1.0
      deyp1(iy) = deyp1(iy) + xmt
      If (ityp==21) Then
        dyg1(iy) = dyg1(iy) + 1.0
        deyg1(iy) = deyg1(iy) + xmt
      End If
      300 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 400
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 400
      dmyp1(imt) = dmyp1(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg1(imt) = dmyg1(imt) + 1.0/xmt
      End If
      400 Continue
    End Do
  End Do

  Do i = 1, nsg
    Do j = 1, njsg(i)
      ityp = k2sg(i, j)
      px = pxsg(i, j)
      py = pysg(i, j)
      pz = pzsg(i, j)
      pe = pesg(i, j)
      pm = pmsg(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZ) .GE. PE) THEN
!               PRINT *, ' IN HJANA1, INDP STR ', I, ' PART ', J
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', PE
!               PRINT *, ' XM = ', PM
!               GOTO 600
!            END IF
!            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA1 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If

      iy = 1 + int(abs(rap)/dy)
!lin-8/2014 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 500
      If (iy<1 .Or. iy>50) Goto 500
      dyp1(iy) = dyp1(iy) + 1.0
      deyp1(iy) = deyp1(iy) + xmt
      If (ityp==21) Then
        dyg1(iy) = dyg1(iy) + 1.0
        deyg1(iy) = deyg1(iy) + xmt
      End If
      500 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 600
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 600
      dmyp1(imt) = dmyp1(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg1(imt) = dmyg1(imt) + 1.0/xmt
      End If
      600 Continue
    End Do
  End Do

  Do i = 1, ihnt2(1)
    yr = sqrt(yp(1,i)**2+yp(2,i)**2)
    ir = 1 + int(yr/dr)
!lin-4/2008 protect against out-of-bound errors:
!         IF (IR .GT. 50) GOTO 601
    If (ir>50 .Or. ir<1) Goto 601
    dnrpj1(ir) = dnrpj1(ir) + 1.0
    dnrtt1(ir) = dnrtt1(ir) + 1.0
    601 Continue
  End Do

  Do i = 1, ihnt2(3)
    yr = sqrt(yt(1,i)**2+yt(2,i)**2)
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 602
    dnrtg1(ir) = dnrtg1(ir) + 1.0
    dnrtt1(ir) = dnrtt1(ir) + 1.0
    602 Continue
  End Do

  Do i = 1, nsg
    y1 = 0.5*(yp(1,iasg(i,1))+yt(1,iasg(i,2)))
    y2 = 0.5*(yp(2,iasg(i,1))+yt(2,iasg(i,2)))
    yr = sqrt(y1**2+y2**2)
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 603
    dnrin1(ir) = dnrin1(ir) + 1.0
    dnrtt1(ir) = dnrtt1(ir) + 1.0
    603 Continue
  End Do

  Do i = 1, mul
    ityp = ityp0(i)
    px = sngl(px0(i))
    py = sngl(py0(i))
    pz = sngl(pz0(i))
    pe = sngl(e0(i))
    pm = sngl(xmass0(i))
    xmt = sqrt(px**2+py**2+pm**2)
    dxmt = xmt - pm
!lin-9/2012 determine rapidity more generally:
!         IF (ABS(PZ) .GE. PE) THEN
!            PRINT *, ' IN HJANA1, GLUON ', I
!            PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!            PRINT *, ' PZ = ', PZ, ' EE = ', PE
!            PRINT *, ' XM = ', PM
!            GOTO 800
!         END IF
!         RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
    If (xmt>0.) Then
      rap = asinh(pz/xmt)
    Else
      Print *, ' IN HJANA1 mt=0'
      rap = 1000000.0*sign(1., pz)
    End If

    iy = 1 + int(abs(rap)/dy)
!lin-8/2014 prevent possible segmentation fault (due to IY<=0):
!         IF (IY .GT. 50) GOTO 700
    If (iy<1 .Or. iy>50) Goto 700
    dyg1c(iy) = dyg1c(iy) + 1.0
    deyg1c(iy) = deyg1c(iy) + xmt
    700 Continue
    If (rap>ymax .Or. rap<=ymin) Goto 800
    imt = 1 + int(dxmt/dmt)
    If (imt>50) Goto 800
    dmyg1c(imt) = dmyg1c(imt) + 1.0/xmt
    800 Continue
  End Do
!.....count number of particles
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      nsubp = nsubp + 1
      If (kfpj(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do

  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      nsubp = nsubp + 1
      If (kftj(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do

  Do i = 1, nsg
    Do j = 1, njsg(i)
      nsubp = nsubp + 1
      If (k2sg(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do
  nisg = nisg + nsg
  If (iout==1) Then
!bzdbg2/16/99
!      PRINT *, ' in HJANA1 '
!      PRINT *, ' total number of partons = ', nsubp
!      PRINT *, ' total number of gluons = ', nsubg, MUL
!      PRINT *, ' number of projectile strings = ', IHNT2(1)
!      PRINT *, ' number of target strings = ', IHNT2(3)
!      PRINT *, ' number of independent strings = ', NSG
    Print *, ' in HJANA1 '
    Print *, ' total number of partons = ', nsubp/iw
    Print *, ' total number of gluons = ', nsubg/iw
!      PRINT *, ' number of projectile strings = ', IHNT2(1)
!      PRINT *, ' number of target strings = ', IHNT2(3)
    Print *, ' number of independent strings = ', nisg/iw
!bzdbg2/16/99end
  End If
!
  Return
End Subroutine hjana1

!-----------------------------------------------------------------------

!.....analysis subroutine in ZPC after generation of additional initial
!.....phase space distributions.

Subroutine hjan1a
  Parameter (maxptn=400001)
  Parameter (dgx=0.2, dgy=0.2, dt=0.2)
  Dimension dgxg1a(50), dgyg1a(50), dtg1a(50)
  Dimension sgxg1a(50), sgyg1a(50), stg1a(50)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
!c      SAVE /PARA1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
!c      SAVE /prec2/
  Common /arevt/iaevt, iarun, miss
!c      SAVE /AREVT/
  Common /arout/iout
!c      SAVE /AROUT/
  Save
  Data iw/0/

  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dgxg1a(i) = sgxg1a(i)
      dgyg1a(i) = sgyg1a(i)
      dtg1a(i) = stg1a(i)
    End Do
  Else
    Do i = 1, 50
      sgxg1a(i) = dgxg1a(i)
      sgyg1a(i) = dgyg1a(i)
      stg1a(i) = dtg1a(i)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
!.....analysis
  Do i = 1, mul
    igx = 1 + int(sngl(abs(gx5(i)))/dgx)
!lin-4/2008 protect against out-of-bound errors:
!         IF (IGX .GT. 50) GOTO 100
    If (igx>50 .Or. igx<1) Goto 100
    dgxg1a(igx) = dgxg1a(igx) + 1.0
    100 Continue
    igy = 1 + int(sngl(abs(gy5(i)))/dgy)
    If (igy>50 .Or. igy<1) Goto 200
    dgyg1a(igy) = dgyg1a(igy) + 1.0
    200 Continue
!lin-9/2015 to avoid Floating-Point Exception:
!         IT = 1 + int(sngl(SQRT(FT5(I) ** 2 - GZ5(I) ** 2)) / DT)
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '1:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      it = 1
    Else
      it = 1 + int(sqrt(diff2)/dt)
    End If
!
    If (it>50 .Or. it<1) Goto 300
    dtg1a(it) = dtg1a(it) + 1.0
    300 Continue
  End Do
  Call hjan1b
!
  Return
End Subroutine hjan1a

!-----------------------------------------------------------------------

!.....analysis subroutine in HJAN1A

Subroutine hjan1b
  Parameter (maxptn=400001, maxstr=150001)
  Parameter (dr=0.2, dt=0.2)
  Dimension dnrg1b(50), dtg1b(50)
  Dimension snrg1b(50), stg1b(50)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
!c      SAVE /PARA1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
!c      SAVE /prec2/
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
!c      SAVE /ilist8/
  Common /srec1/nsp, nst, nsi
!c      SAVE /SREC1/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /arevt/iaevt, iarun, miss
!c      SAVE /AREVT/
  Common /arout/iout
!c      SAVE /AROUT/
  Save
  Data iw/0/

  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dnrg1b(i) = snrg1b(i)
      dtg1b(i) = stg1b(i)
    End Do
  Else
    Do i = 1, 50
      snrg1b(i) = dnrg1b(i)
      stg1b(i) = dtg1b(i)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
!.....analysis
  Do i = 1, mul
    j = lstrg1(i)

    If (j<=nsp) Then
      k = j
      gx0 = yp(1, j)
      gy0 = yp(2, j)
    Else If (j<=nsp+nst) Then
      k = j - nsp
      gx0 = yt(1, k)
      gy0 = yt(2, k)
    Else
      k = j - nsp - nst
      gx0 = 0.5*(yp(1,iasg(k,1))+yt(1,iasg(k,2)))
      gy0 = 0.5*(yp(2,iasg(k,1))+yt(2,iasg(k,2)))
    End If
    r0 = sqrt((sngl(gx5(i))-gx0)**2+(sngl(gy5(i))-gy0)**2)
    ir = 1 + int(r0/dr)
    If (ir>50 .Or. ir<1) Goto 100
    dnrg1b(ir) = dnrg1b(ir) + 1.0
    100 Continue
!lin-9/2015 to avoid Floating-Point Exception:
!         TAU7 = SQRT(sngl(FT5(I) ** 2 - GZ5(I) ** 2))
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '5:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      tau7 = 1E-6
    Else
      tau7 = sqrt(diff2)
    End If
!
    it = 1 + int(tau7/dt)
    If (it>50 .Or. it<1) Goto 200
    dtg1b(it) = dtg1b(it) + 1.0
    200 Continue
  End Do
!
  Return
End Subroutine hjan1b



!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!.....analysis subroutine in HIJING after parton cascade evolution
Subroutine hjana2
!
  Parameter (ymax=1.0, ymin=-1.0)
  Parameter (dmt=0.05, dy=0.2)
  Parameter (dr=0.2, dt=0.2)
  Parameter (maxptn=400001)
  Parameter (maxstr=150001)
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Dimension dyp2(50), dmyp2(200), deyp2(50)
  Dimension dyg2(50), dmyg2(200), deyg2(50)
  Dimension snyp2(50), smyp2(200), seyp2(50)
  Dimension snyg2(50), smyg2(200), seyg2(50)
  Dimension dnrpj2(50), dnrtg2(50), dnrin2(50), dnrtt2(50)
  Dimension dtpj2(50), dttg2(50), dtin2(50), dttot2(50)
  Dimension dyg2c(50), dmyg2c(50), deyg2c(50)
  Dimension snrpj2(50), snrtg2(50), snrin2(50), snrtt2(50)
  Dimension stpj2(50), sttg2(50), stin2(50), sttot2(50)
  Dimension snyg2c(50), smyg2c(50), seyg2c(50)
  Double Precision ataui, zt1, zt2, zt3
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
!c      SAVE /PARA1/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
!c      SAVE /SREC2/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
!c      SAVE /prec2/
  Common /arevt/iaevt, iarun, miss
!c      SAVE /AREVT/
  Common /arout/iout
!c      SAVE /AROUT/
  Common /anim/nevent, isoft, isflag, izpc
!c      SAVE /anim/
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
!c      SAVE /SOFT/
  Save
  Data iw/0/

  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 200
      dmyp2(i) = smyp2(i)
      dmyg2(i) = smyg2(i)
    End Do

    Do i = 1, 50
      dyp2(i) = snyp2(i)
      deyp2(i) = seyp2(i)
      dyg2(i) = snyg2(i)
      deyg2(i) = seyg2(i)
      dnrpj2(i) = snrpj2(i)
      dnrtg2(i) = snrtg2(i)
      dnrin2(i) = snrin2(i)
      dnrtt2(i) = snrtt2(i)
      dtpj2(i) = stpj2(i)
      dttg2(i) = sttg2(i)
      dtin2(i) = stin2(i)
      dttot2(i) = sttot2(i)
      dyg2c(i) = snyg2c(i)
      dmyg2c(i) = smyg2c(i)
      deyg2c(i) = seyg2c(i)
    End Do
    nsubp = nsubps
    nsubg = nsubgs
    nisg = nisgs
  Else
    Do i = 1, 200
      smyp2(i) = dmyp2(i)
      smyg2(i) = dmyg2(i)
    End Do

    Do i = 1, 50
      snyp2(i) = dyp2(i)
      seyp2(i) = deyp2(i)
      snyg2(i) = dyg2(i)
      seyg2(i) = deyg2(i)
      snrpj2(i) = dnrpj2(i)
      snrtg2(i) = dnrtg2(i)
      snrin2(i) = dnrin2(i)
      snrtt2(i) = dnrtt2(i)
      stpj2(i) = dtpj2(i)
      sttg2(i) = dttg2(i)
      stin2(i) = dtin2(i)
      sttot2(i) = dttot2(i)
      snyg2c(i) = dyg2c(i)
      smyg2c(i) = dmyg2c(i)
      seyg2c(i) = deyg2c(i)
    End Do
    nsubps = nsubp
    nsubgs = nsubg
    nisgs = nisg
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If

!lin-4/28/01:
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Goto 510

!.....analysis
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      ityp = kfpj(i, j)
      px = pjpx(i, j)
      py = pjpy(i, j)
      pz = pjpz(i, j)
      pe = pjpe(i, j)
      pm = pjpm(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
!lin-9/2012 determine rapidity more generally:
!cbzdbg2/16/99
!c            IF (ABS(PZ) .GE. PE) GOTO 200
!            IF (ABS(PZ) .GE. PE) THEN
!               PRINT *, ' IN HJANA2, PROJ STR ', I, ' PART ', J
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', PE
!               PRINT *, ' XM = ', PM
!               GOTO 200
!            END IF
!cbzdbg2/16/99end
!            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA2 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If

      iy = 1 + int(abs(rap)/dy)
!lin-8/2014 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 100
      If (iy<1 .Or. iy>50) Goto 100
      dyp2(iy) = dyp2(iy) + 1.0
      deyp2(iy) = deyp2(iy) + xmt
      If (ityp==21) Then
        dyg2(iy) = dyg2(iy) + 1.0
        deyg2(iy) = deyg2(iy) + xmt
      End If
      100 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 200
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 200
      dmyp2(imt) = dmyp2(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg2(imt) = dmyg2(imt) + 1.0/xmt
      End If
      200 Continue
    End Do
  End Do

  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      ityp = kftj(i, j)
      px = pjtx(i, j)
      py = pjty(i, j)
      pz = pjtz(i, j)
      pe = pjte(i, j)
      pm = pjtm(i, j)
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
!lin-9/2012 determine rapidity more generally:
!cbzdbg2/16/99
!c            IF (ABS(PZ) .GE. PE) GOTO 400
!            IF (ABS(PZ) .GE. PE) THEN
!               PRINT *, ' IN HJANA2, TARG STR ', I, ' PART ', J
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', PE
!               PRINT *, ' XM = ', PM
!               GOTO 400
!            END IF
!cbzdbg2/16/99end
!            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA2 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If

      iy = 1 + int(abs(rap)/dy)
!lin-8/2014 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 300
      If (iy<1 .Or. iy>50) Goto 300
      dyp2(iy) = dyp2(iy) + 1.0
      deyp2(iy) = deyp2(iy) + xmt
      If (ityp==21) Then
        dyg2(iy) = dyg2(iy) + 1.0
        deyg2(iy) = deyg2(iy) + xmt
      End If
      300 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 400
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 400
      dmyp2(imt) = dmyp2(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg2(imt) = dmyg2(imt) + 1.0/xmt
      End If
      400 Continue
    End Do
  End Do

!lin-4/28/01:
  510 Continue

  Do i = 1, nsg
!lin-4/25/01 soft3:
!         DO J = 1, NJSG(I)
    nj = njsg(i)
    If (isoft==3 .Or. isoft==4 .Or. isoft==5) nj = njsgs(i)
    Do j = 1, nj
!lin-4/25/01-end

      ityp = k2sg(i, j)
      px = pxsg(i, j)
      py = pysg(i, j)
      pz = pzsg(i, j)
      pe = pesg(i, j)
      pm = pmsg(i, j)
!lin-4/25/01 soft3:
      If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
        ityp = k2sgs(i, j)
        px = sngl(pxsgs(i,j))
        py = sngl(pysgs(i,j))
        pz = sngl(pzsgs(i,j))
        pe = sngl(pesgs(i,j))
        pm = sngl(pmsgs(i,j))
      End If
!lin-4/25/01-end

!lin-9/2012 determine rapidity more generally:
      xmt = sqrt(px**2+py**2+pm**2)
      dxmt = xmt - pm
!cbzdbg2/16/99
!c            IF (ABS(PZ) .GE. PE) GOTO 600
!            IF (ABS(PZ) .GE. PE) THEN
!               PRINT *, ' IN HJANA2, INDP STR ', I, ' PART ', J
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', PE
!               PRINT *, ' XM = ', PM
!               GOTO 600
!            END IF
!cbzdbg2/16/99end
!            RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
      If (xmt>0.) Then
        rap = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA2 mt=0'
        rap = 1000000.0*sign(1., pz)
      End If

      iy = 1 + int(abs(rap)/dy)
!lin-8/2014 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 500
      If (iy<1 .Or. iy>50) Goto 500
      dyp2(iy) = dyp2(iy) + 1.0
      deyp2(iy) = deyp2(iy) + xmt
      If (ityp==21) Then
        dyg2(iy) = dyg2(iy) + 1.0
        deyg2(iy) = deyg2(iy) + xmt
      End If
      500 Continue
      If (rap>ymax .Or. rap<=ymin) Goto 600
      imt = 1 + int(dxmt/dmt)
      If (imt>200) Goto 600
      dmyp2(imt) = dmyp2(imt) + 1.0/xmt
      If (ityp==21) Then
        dmyg2(imt) = dmyg2(imt) + 1.0/xmt
      End If
      600 Continue
    End Do
  End Do

!lin-4/28/01:
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Goto 520

  Do i = 1, ihnt2(1)
    j = i
    yr = sqrt(sngl(zt1(j)**2+zt2(j)**2))
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 601
    dnrpj2(ir) = dnrpj2(ir) + 1.0
    dnrtt2(ir) = dnrtt2(ir) + 1.0
    601 Continue
    it = 1 + int(sngl(ataui(j))/dt)
    If (it>50 .Or. it<1) Goto 602
    dtpj2(it) = dtpj2(it) + 1.0
    dttot2(it) = dttot2(it) + 1.0
    602 Continue
  End Do

  Do i = 1, ihnt2(3)
    j = i + ihnt2(1)
    yr = sqrt(sngl(zt1(j)**2+zt2(j)**2))
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 603
    dnrtg2(ir) = dnrtg2(ir) + 1.0
    dnrtt2(ir) = dnrtt2(ir) + 1.0
    603 Continue
    it = 1 + int(sngl(ataui(j))/dt)
    If (it>50 .Or. it<1) Goto 604
    dttg2(it) = dttg2(it) + 1.0
    dttot2(it) = dttot2(it) + 1.0
    604 Continue
  End Do

!lin-4/28/01:
  520 Continue

  Do i = 1, nsg
    j = i + ihnt2(1) + ihnt2(3)
!lin-4/28/01:
    If (isoft==3 .Or. isoft==4 .Or. isoft==5) j = i

    yr = sqrt(sngl(zt1(j)**2+zt2(j)**2))
    ir = 1 + int(yr/dr)
    If (ir>50 .Or. ir<1) Goto 605
    dnrin2(ir) = dnrin2(ir) + 1.0
    dnrtt2(ir) = dnrtt2(ir) + 1.0
    605 Continue
    it = 1 + int(sngl(ataui(j))/dt)
    If (it>50 .Or. it<1) Goto 606
    dtin2(it) = dtin2(it) + 1.0
    dttot2(it) = dttot2(it) + 1.0
    606 Continue
  End Do

  Do i = 1, mul
    ityp = ityp5(i)
    px = sngl(px5(i))
    py = sngl(py5(i))
    pz = sngl(pz5(i))
    pe = sngl(e5(i))
    pm = sngl(xmass5(i))
!lin-9/2012 determine rapidity more generally:
    xmt = sqrt(px**2+py**2+pm**2)
    dxmt = xmt - pm
!cbzdbg2/16/99
!c            IF (ABS(PZ) .GE. PE) GOTO 800
!         IF (ABS(PZ) .GE. PE) THEN
!            PRINT *, ' IN HJANA2, GLUON ', I
!            PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!            PRINT *, ' PZ = ', PZ, ' EE = ', PE
!            PRINT *, ' XM = ', PM
!            GOTO 800
!         END IF
!cbzdbg2/16/99end
!         RAP = 0.5 * LOG((PE + PZ +1e-5) / (PE - PZ + 1e-5))
    If (xmt>0.) Then
      rap = asinh(pz/xmt)
    Else
      Print *, ' IN HJANA2 mt=0'
      rap = 1000000.0*sign(1., pz)
    End If

    iy = 1 + int(abs(rap)/dy)
!lin-9/2012 prevent possible segmentation fault (due to IY<=0):
!         IF (IY .GT. 50) GOTO 700
    If (iy<1 .Or. iy>50) Goto 700
    dyg2c(iy) = dyg2c(iy) + 1.0
    deyg2c(iy) = deyg2c(iy) + xmt
    700 Continue
    If (rap>ymax .Or. rap<=ymin) Goto 800
    imt = 1 + int(dxmt/dmt)
    If (imt>50) Goto 800
    dmyg2c(imt) = dmyg2c(imt) + 1.0/xmt
    800 Continue
  End Do

!lin-4/25/01 soft3:
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Goto 530

!.....count number of particles
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      nsubp = nsubp + 1
      If (kfpj(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do

  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      nsubp = nsubp + 1
      If (kftj(i,j)==21) nsubg = nsubg + 1
    End Do
  End Do

!lin-4/25/01 soft3:
  530 Continue

  Do i = 1, nsg
!lin-4/25/01 soft3:
!         DO J = 1, NJSG(I)
    nj = njsg(i)
    If (isoft==3 .Or. isoft==4 .Or. isoft==5) nj = njsgs(i)
    Do j = 1, nj
!lin-4/25/01-end

      nsubp = nsubp + 1

!lin-4/25/01
!            IF (K2SG(I, J) .EQ. 21) nsubg = nsubg + 1
      If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
        If (k2sgs(i,j)==21) nsubg = nsubg + 1
      Else
        If (k2sg(i,j)==21) nsubg = nsubg + 1
      End If
!lin-4/25/01-end
    End Do
  End Do
!bzdbg2/16/99
  nisg = nisg + nsg

  If (iout==1) Then
!bzdbg2/16/99end
!bzdbg2/16/99
!      PRINT *, ' in HJANA2 '
!      PRINT *, ' total number of partons = ', nsubp
!      PRINT *, ' total number of gluons = ', nsubg, MUL
!      PRINT *, ' number of projectile strings = ', IHNT2(1)
!      PRINT *, ' number of target strings = ', IHNT2(3)
!      PRINT *, ' number of independent strings = ', NSG
    Print *, ' in HJANA2 '
    Print *, ' total number of partons = ', nsubp/iw
    Print *, ' total number of gluons = ', nsubg/iw
!      PRINT *, ' number of projectile strings = ', IHNT2(1)
!      PRINT *, ' number of target strings = ', IHNT2(3)
    Print *, ' number of independent strings = ', nisg/iw
  End If

  Call hjan2a
  Call hjan2b

  Return
End Subroutine hjana2

!-----------------------------------------------------------------------

!.....subroutine called by HJANA2
Subroutine hjan2a

  Parameter (dgx=0.2, dgy=0.2, dt=0.2)
  Parameter (maxptn=400001, maxstr=150001)
  Dimension dgxp2a(50), dgyp2a(50), dtp2a(50)
  Dimension dgxg2a(50), dgyg2a(50), dtg2a(50)
  Dimension sgxp2a(50), sgyp2a(50), stp2a(50)
  Dimension sgxg2a(50), sgyg2a(50), stg2a(50)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Common /para1/mul
!c      SAVE /PARA1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
!c      SAVE /prec2/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /arevt/iaevt, iarun, miss
!c      SAVE /AREVT/
  Common /arout/iout
!c      SAVE /AROUT/
  Save
  Data iw/0/

  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dgxp2a(i) = sgxp2a(i)
      dgyp2a(i) = sgyp2a(i)
      dtp2a(i) = stp2a(i)
      dgxg2a(i) = sgxg2a(i)
      dgyg2a(i) = sgyg2a(i)
      dtg2a(i) = stg2a(i)
    End Do
  Else
    Do i = 1, 50
      sgxp2a(i) = dgxp2a(i)
      sgyp2a(i) = dgyp2a(i)
      stp2a(i) = dtp2a(i)
      sgxg2a(i) = dgxg2a(i)
      sgyg2a(i) = dgyg2a(i)
      stg2a(i) = dtg2a(i)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
!.....analysis
  Do i = 1, ihnt2(1)
    Do j = 1, npj(i)
      If (kfpj(i,j)/=21) Then
        igx = 1 + int(abs(yp(1,i))/dgx)
        If (igx>50 .Or. igx<1) Goto 100
        dgxp2a(igx) = dgxp2a(igx) + 1.0
        100 Continue
        igy = 1 + int(abs(yp(2,i))/dgy)
        If (igy>50 .Or. igy<1) Goto 200
        dgyp2a(igy) = dgyp2a(igy) + 1.0
        200 Continue
        it = 1
        dtp2a(it) = dtp2a(it) + 1.0
      End If
    End Do
  End Do

  Do i = 1, ihnt2(3)
    Do j = 1, ntj(i)
      If (kftj(i,j)/=21) Then
        igx = 1 + int(abs(yt(1,i))/dgx)
        If (igx>50 .Or. igx<1) Goto 300
        dgxp2a(igx) = dgxp2a(igx) + 1.0
        300 Continue
        igy = 1 + int(abs(yt(2,i))/dgy)
        If (igy>50 .Or. igy<1) Goto 400
        dgyp2a(igy) = dgyp2a(igy) + 1.0
        400 Continue
        it = 1
        dtp2a(it) = dtp2a(it) + 1.0
      End If
    End Do
  End Do

  Do i = 1, nsg
    Do j = 1, njsg(i)
      If (k2sg(i,j)/=21) Then
        igx = 1 + int(abs(0.5*(yp(1,iasg(i,1))+yt(1,iasg(i,2))))/dgx)
        If (igx>50 .Or. igx<1) Goto 500
        dgxp2a(igx) = dgxp2a(igx) + 1.0
        500 Continue
        igy = 1 + int(abs(0.5*(yp(2,iasg(i,1))+yt(2,iasg(i,2))))/dgy)
        If (igy>50 .Or. igy<1) Goto 600
        dgyp2a(igy) = dgyp2a(igy) + 1.0
        600 Continue
        it = 1
        dtp2a(it) = dtp2a(it) + 1.0
      End If
    End Do
  End Do

  Do i = 1, mul
    igx = 1 + int(abs(sngl(gx5(i)))/dgx)
    If (igx>50 .Or. igx<1) Goto 700
    dgxg2a(igx) = dgxg2a(igx) + 1.0
    dgxp2a(igx) = dgxp2a(igx) + 1.0
    700 Continue
    igy = 1 + int(abs(sngl(gy5(i)))/dgy)
    If (igy>50 .Or. igy<1) Goto 800
    dgyg2a(igy) = dgyg2a(igy) + 1.0
    dgyp2a(igy) = dgyp2a(igy) + 1.0
    800 Continue
!lin-9/2015 to avoid Floating-Point Exception:
!         IT = 1 + int(SQRT(sngl(FT5(I) ** 2 - GZ5(I) ** 2)) / DT)
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '3:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      it = 1
    Else
      it = 1 + int(sqrt(diff2)/dt)
    End If
!
    If (it>50 .Or. it<1) Goto 900
    dtg2a(it) = dtg2a(it) + 1.0
    dtp2a(it) = dtp2a(it) + 1.0
    900 Continue
  End Do
!
  Return
End Subroutine hjan2a

!-----------------------------------------------------------------------

!.....analysis subroutine in HJANA2

Subroutine hjan2b

  Parameter (maxptn=400001)
  Parameter (maxstr=150001)
  Parameter (dr=0.2, dt=0.2)
  Dimension dnrg2b(50), dtg2b(-24:25)
  Dimension snrg2b(50), stg2b(-24:25)
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Double Precision ataui, zt1, zt2, zt3
  Common /para1/mul
!c      SAVE /PARA1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
!c      SAVE /prec2/
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
!c      SAVE /ilist8/
  Common /srec1/nsp, nst, nsi
!c      SAVE /SREC1/
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
!c      SAVE /SREC2/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /arevt/iaevt, iarun, miss
!c      SAVE /AREVT/
  Common /arout/iout
!c      SAVE /AROUT/
  Save
  Data iw/0/

  If (isevt==iaevt .And. isrun==iarun) Then
    Do i = 1, 50
      dnrg2b(i) = snrg2b(i)
      dtg2b(i-25) = stg2b(i-25)
    End Do
  Else
    Do i = 1, 50
      snrg2b(i) = dnrg2b(i)
      stg2b(i-25) = dtg2b(i-25)
    End Do
    isevt = iaevt
    isrun = iarun
    iw = iw + 1
  End If
!.....analysis
  Do i = 1, mul
    j = lstrg1(i)
    gx0 = sngl(zt1(j))
    gy0 = sngl(zt2(j))
    r0 = sqrt((sngl(gx5(i))-gx0)**2+(sngl(gy5(i))-gy0)**2)
    ir = 1 + int(r0/dr)
    If (ir>50 .Or. ir<1) Goto 100
    dnrg2b(ir) = dnrg2b(ir) + 1.0
    100 Continue
!lin-9/2015 to avoid Floating-Point Exception:
!         TAU7 = SQRT(sngl(FT5(I) ** 2 - GZ5(I) ** 2))
    diff2 = sngl(ft5(i)**2-gz5(i)**2)
    If (diff2<0.) Then
      Write (6, *) '4:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      tau7 = 1E-6
    Else
      tau7 = sqrt(diff2)
    End If
!
    dtau = tau7 - sngl(ataui(j))
    it = 1 + int(dtau/dt)
!bzdbg2/21/99
!         IF (ABS(IT) .GT. 25) GOTO 200
    If (it>25 .Or. it<-24) Goto 200
!bzdbg2/21/99end
    dtg2b(it) = dtg2b(it) + 1.0
    200 Continue
  End Do
!
  Return
End Subroutine hjan2b

!-----------------------------------------------------------------------

!.....analysis subroutine before ARTMN
Subroutine hjana3
!
  Parameter (maxstr=150001, maxr=1)
!.....y cut for mt spectrum
  Parameter (ymin=-1.0, ymax=1.0)
!bz11/7/99 end
!.....bin width for mt spectrum and y spectrum
  Parameter (dmt=0.05, dy=0.2)
  Double Precision v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
  Dimension dndyh3(50), dmyh3(50), deyh3(50)
  Common /run/num
!c      SAVE /RUN/
  Common /arerc1/multi1(maxr)
!c      SAVE /ARERC1/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
!c      SAVE /ARPRC1/
  Common /arout/iout
!c      SAVE /AROUT/
  Common /iflow/v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
!c      SAVE /iflow/
  Save
  Data iw/0/

  iw = iw + 1
  Do j = 1, num
    Do i = 1, multi1(j)
      ityp = ityp1(i, j)
      If (ityp>-100 .And. ityp<100) Goto 200
      px = px1(i, j)
      py = py1(i, j)
      pz = pz1(i, j)
      ee = ee1(i, j)
      xm = xm1(i, j)
      xmt = sqrt(px**2+py**2+xm**2)
      dxmt = xmt - xm
!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZ) .GE. EE) THEN
!               PRINT *, 'IN HJANA3'
!               PRINT *, ' PARTICLE ', I, ' RUN ', J, 'PREC ERR'
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', EE
!               PRINT *, ' XM = ', XM
!               GOTO 200
!            END IF
!            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
      If (xmt>0.) Then
        y = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA3 mt=0'
        y = 1000000.0*sign(1., pz)
      End If

!.....rapidity cut for the rapidity distribution
!            IY = 1 + int(ABS(Y) / DY)
!lin-8/2014 no rapidity shift here:
!            IY = 1 + int((Y+10.) / DY)
      iy = 1 + int(y/dy)
!lin-9/2012 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 100
      If (iy<1 .Or. iy>50) Goto 100
      dndyh3(iy) = dndyh3(iy) + 1.0
      deyh3(iy) = deyh3(iy) + xmt
      100 Continue
!.....insert rapidity cut for mt spectrum here
      If (y<ymin .Or. y>=ymax) Goto 200
      imt = 1 + int(dxmt/dmt)
      If (imt>50) Goto 200
      dmyh3(imt) = dmyh3(imt) + 1.0/xmt
      200 Continue
    End Do
  End Do
!
  Return
End Subroutine hjana3

!-----------------------------------------------------------------------

!.....analysis subroutine after ARTMN
Subroutine hjana4
  Parameter (maxstr=150001, maxr=1)
!.....y cut for mt spectrum
!bz11/7/99
!      PARAMETER (YMIN = -0.5, YMAX = 0.5)
  Parameter (ymin=-1.0, ymax=1.0)
!bz11/7/99 end
!.....bin width for mt spectrum and y spectrum
  Parameter (dmt=0.05, dy=0.2)

  Dimension dndyh4(50), dmyh4(50), deyh4(50)
  Common /run/num
!c      SAVE /RUN/
  Common /arerc1/multi1(maxr)
!c      SAVE /ARERC1/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
!c      SAVE /ARPRC1/
  Common /arout/iout
!c      SAVE /AROUT/
  Common /fflow/v2f, etf, xmultf, v2fpi, xmulpi
!c      SAVE /fflow/
  Save
  Data iw/0/

  iw = iw + 1
  Do j = 1, num
    Do i = 1, multi1(j)
      ityp = ityp1(i, j)
      If (ityp>-100 .And. ityp<100) Goto 200
      px = px1(i, j)
      py = py1(i, j)
      pz = pz1(i, j)
      ee = ee1(i, j)
      xm = xm1(i, j)
      xmt = sqrt(px**2+py**2+xm**2)
      dxmt = xmt - xm
!lin-9/2012 determine rapidity more generally:
!            IF (ABS(PZ) .GE. EE) THEN
!               PRINT *, 'IN HJANA4'
!               PRINT *, ' PARTICLE ', I, ' RUN ', J, 'PREC ERR'
!               PRINT *, ' FLAV = ', ITYP, ' PX = ', PX, ' PY = ', PY
!               PRINT *, ' PZ = ', PZ, ' EE = ', EE
!               PRINT *, ' XM = ', XM
!               GOTO 200
!            END IF
!            Y = 0.5 * LOG((EE + PZ +1e-5) / (EE - PZ + 1e-5))
      If (xmt>0.) Then
        y = asinh(pz/xmt)
      Else
        Print *, ' IN HJANA4 mt=0'
        y = 1000000.0*sign(1., pz)
      End If

!.....rapidity cut for the rapidity distribution
!            IY = 1 + int(ABS(Y) / DY)
!lin-8/2014 no rapidity shift here:
!            IY = 1 + int((Y+10.) / DY)
      iy = 1 + int(y/dy)
!lin-9/2012 prevent possible segmentation fault (due to IY<=0):
!            IF (IY .GT. 50) GOTO 100
      If (iy<1 .Or. iy>50) Goto 100
      dndyh4(iy) = dndyh4(iy) + 1.0
      deyh4(iy) = deyh4(iy) + xmt
      100 Continue
!.....insert rapidity cut for mt spectrum here
      If (y<ymin .Or. y>=ymax) Goto 200
      imt = 1 + int(dxmt/dmt)
      If (imt>50) Goto 200
      dmyh4(imt) = dmyh4(imt) + 1.0/xmt
      200 Continue
    End Do
  End Do
!
  Return
End Subroutine hjana4

!=======================================================================

!.....subroutine to get average values for different strings

Subroutine zpstrg

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Parameter (maxstr=150001)
!      REAL*4 YP, YT, PXSG, PYSG, PZSG, PESG, PMSG, HIPR1, HINT1, BB
  Real yp, yt, pxsg, pysg, pzsg, pesg, pmsg, hipr1, hint1, bb

  Common /para1/mul
!c      SAVE /PARA1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
!c      SAVE /prec2/
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
!c      SAVE /ilist8/
  Common /srec1/nsp, nst, nsi
!c      SAVE /SREC1/
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
!c      SAVE /SREC2/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
!bz6/28/99 flow1
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
!bz6/28/99 flow1 end
  Common /anim/nevent, isoft, isflag, izpc
!c      SAVE /anim/
  Common /strg/np(maxstr)
!c      SAVE /strg/
!lin-6/06/02 test local freezeout:
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
!c      SAVE /frzprc/
  Save

!lin-6/06/02 test local freezeout for string melting,
!     use space-time values at local freezeout saved in /frzprc/:
  If (isoft==5) Then
    Do i = 1, mul
      ityp5(i) = idfrz(i)
      gx5(i) = gxfrz(i)
      gy5(i) = gyfrz(i)
      gz5(i) = gzfrz(i)
      ft5(i) = ftfrz(i)
      px5(i) = pxfrz(i)
      py5(i) = pyfrz(i)
      pz5(i) = pzfrz(i)
      e5(i) = efrz(i)
      xmass5(i) = xmfrz(i)
    End Do
  End If
!lin-6/06/02-end

  Do i = 1, maxstr
    ataui(i) = 0D0
    zt1(i) = 0D0
    zt2(i) = 0D0
!lin-4/25/03 add zt3(I) to track longitudinal positions of partons/strings:
    zt3(i) = 0D0
    np(i) = 0
  End Do
  Do i = 1, mul
    istrg = lstrg1(i)
!lin-9/2015 to avoid Floating-Point Exception:
!         TAU7 = SQRT(FT5(I) ** 2 - GZ5(I) ** 2)
    diff2 = ft5(i)**2 - gz5(i)**2
    If (diff2<0D0) Then
      Write (6, *) '2:I,ft5,gz5,diff2=', i, ft5(i), gz5(i), diff2
      tau7 = 1D-6
    Else
      tau7 = dsqrt(diff2)
    End If
!
    ataui(istrg) = ataui(istrg) + tau7
    zt1(istrg) = zt1(istrg) + gx5(i)
    zt2(istrg) = zt2(istrg) + gy5(i)
    zt3(istrg) = zt3(istrg) + gz5(i)
    np(istrg) = np(istrg) + 1
  End Do

  nstr = nsp + nst + nsi

!lin-7/03/01 correct averaging on transverse coordinates, no shift needed:
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
    Do i = 1, nstr
      If (np(i)/=0) Then
        ataui(i) = ataui(i)/np(i)
        zt1(i) = zt1(i)/np(i)
        zt2(i) = zt2(i)/np(i)
        zt3(i) = zt3(i)/np(i)
      End If
    End Do
    Return
  End If
!lin-7/03/01-end

  Do i = 1, nstr
    If (np(i)/=0) Then
      ataui(i) = ataui(i)/np(i)
      zt1(i) = zt1(i)/np(i)
      zt2(i) = zt2(i)/np(i)
      zt3(i) = zt3(i)/np(i)
    Else
      If (i<=nsp) Then
        j = i
        zt1(i) = dble(yp(1,j))
        zt2(i) = dble(yp(2,j))
        zt3(i) = 0D0
      Else If (i>nsp .And. i<=nsp+nst) Then
        j = i - nsp
        zt1(i) = dble(yt(1,j))
        zt2(i) = dble(yt(2,j))
        zt3(i) = 0D0
      Else
        j = i - nsp - nst
        zt1(i) = 0.5D0*dble((yp(1,iasg(j,1))+yt(1,iasg(j,2))))
        zt2(i) = 0.5D0*dble((yp(2,iasg(j,1))+yt(2,iasg(j,2))))
        zt3(i) = 0D0
      End If
    End If
  End Do

!bz6/28/99 flow1
  bb = hint1(19)
  Do i = 1, nstr
    If (np(i)/=0) Then
      shift = 0D0
    Else
      shift = 0.5D0*dble(bb)
    End If
    If (i<=nsp) Then
      zt1(i) = zt1(i) + shift
    Else If (i>nsp .And. i<=nsp+nst) Then
      zt1(i) = zt1(i) - shift
    End If
  End Do
!bz6/28/99 flow1 end
!
  Return
End Subroutine zpstrg

!lin-10/01/03 random number generator for f77:
Function ranart(nseed)
  Save
!lin-4/2008 ran(nseed) is renamed to avoid conflict with system functions:
!      ran=rand()
  ranart = rand(0)
!     one may also use the following random number generator in PYTHIA/JETSET:
!      ranart=rlu(0)
  Return
End Function ranart

!lin-3/2009
!     Initialize hadron weights;
!     Can add initial hadrons before the hadron cascade starts (but after ZPC).
Subroutine addhad
  Parameter (maxstr=150001, maxr=1, xmd=1.8756)
  Double Precision smearp, smearh
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /smearz/smearp, smearh
  Common /rndf77/nseed
  Common /para8/idpert, npertd, idxsec
  Save

!     All hadrons at the start of hadron cascade have the weight of 1
!     except those inserted by the user in this subroutine:
  np0 = iaint2(1)
  Do i = 1, np0
    dpertp(i) = 1.
  End Do
!     Specify number, species, weight, initial x,p,m for inserted hadrons here:
  nadd = 0
  tau0 = arpar1(1)
  Do i = np0 + 1, np0 + nadd
    itypar(i) = 42
!lin-5/2012 fix type mismatch:
!         dpertp(I)=1d0/dble(nadd)
    dpertp(i) = 1./float(nadd)
    gxar(i) = 5.*(1.-2.*ranart(nseed))
    gyar(i) = 5.*(1.-2.*ranart(nseed))
    gzar(i) = 2.*(1.-2.*ranart(nseed))
    ftar(i) = 0.
    pxar(i) = 1.
    pyar(i) = 0.
    pzar(i) = 1.
    xmar(i) = xmd
!
    pear(i) = sqrt(pxar(i)**2+pyar(i)**2+pzar(i)**2+xmar(i)**2)
!lin-9/2012 determine rapidity more generally:
!         RAP=0.5*alog((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
    rap = asinh(pzar(i)/sqrt(xmar(i)**2+pxar(i)**2+pyar(i)**2))
!
    vx = pxar(i)/pear(i)
    vy = pyar(i)/pear(i)
!.....give initial formation time shift and boost according to rapidity:
    taui = ftar(i) + tau0
    ftar(i) = taui*cosh(rap)
    gxar(i) = gxar(i) + vx*tau0*cosh(rap)
    gyar(i) = gyar(i) + vy*tau0*cosh(rap)
!     Allow the intial z-position to be different from the Bjorken picture:
    gzar(i) = taui*sinh(rap) + gzar(i)
!         GZAR(I)=TAUI*SINH(RAP)
    zsmear = sngl(smearh)*(2.*ranart(nseed)-1.)
    gzar(i) = gzar(i) + zsmear
  End Do
  iaint2(1) = iaint2(1) + nadd
!
  If (nadd>=1 .And. idpert/=1 .And. idpert/=2) Then
    Write (16, *) 'IDPERT must be 1 or 2 to add initial hadrons,      set NPERTD to 0 if you do not need perturbative deuterons'
    Stop
  End If
  If (iaint2(1)>maxstr) Then
    Write (16, *) 'Too many initial hadrons, array size is exceeded!'
    Stop
  End If
!
  Return
End Subroutine addhad

!lin-8/2014 define function asinh():
Function asinh(x)
  Save
  If (x>0) Then
    asinh = alog(x+sqrt(x**2+1.))
  Else
!     a la suggestion de YP Liu:
    asinh = -alog(-x+sqrt(x**2+1.))
  End If
  Return
End Function asinh


