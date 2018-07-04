!.................... zpc.f
!	PROGRAM ZPC
Subroutine zpcmn
  !       Version: 1.0.1
  !       Author: Bin Zhang
  !       (suggestions, problems -> bzhang@nt1.phys.columbia.edu)
  Implicit Double Precision (A-H, O-Z)
  !lin-4/20/01        PARAMETER (NMAXGL = 16000)
  Parameter (maxptn=400001)
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Save
  !
  !       loop over events
  Do i = 1, nevnt
     ievt = i
     !       generation of the initial condition for one event
     Call inievt
     !      loop over many runs of the same event
     Do j = 1, nsbrun
        isbrun = j
        !       initialization for one run of an event
        Call inirun
        !lin-4/2008 not used:
        !             CALL HJAN1A
3000    Continue
        !       do one collision
        Call zpcrun(*4000)
        Call zpca1
        Goto 3000
4000    Continue
        Call zpca2
     End Do
  End Do
  Call zpcou
  !lin-5/2009 ctest off
  !     5/17/01 calculate v2 for parton already frozen out:
  !        call flowp(3)
  !.....to get average values for different strings
  Call zpstrg
  Return
End Subroutine zpcmn

!*****************************************************************************
!*****************************************************************************

Block Data zpcbdt
   !       set initial values in block data

   Implicit Double Precision (A-H, O-Z)
   Parameter (maxptn=400001)
   Parameter (maxstr=150001)
   Common /para1/mul
   !c      SAVE /para1/
   Common /para2/xmp, xmu, alpha, rscut2, cutof2
   !c      SAVE /para2/
   Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
   !c      SAVE /para3/
   Common /para4/iftflg, ireflg, igeflg, ibstfg
   !c      SAVE /para4/
   Common /para5/iconfg, iordsc
   !c      SAVE /para5/
   Common /para6/centy
   !c      SAVE /para6/
   !lin-6/2009 nsmbbbar and nsmmeson respectively give the total number of
   !     baryons/anti-baryons and mesons for each event:
   !        common /para7/ ioscar
   Common /para7/ioscar, nsmbbbar, nsmmeson
   !c      SAVE /para7/
   Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
   !c      SAVE /prec1/
   Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
   !c      SAVE /prec2/
   Common /prec3/gxs(maxptn), gys(maxptn), gzs(maxptn), fts(maxptn), pxs(maxptn), pys(maxptn), pzs(maxptn), es(maxptn), xmasss(maxptn), ityps(maxptn)
   !c      SAVE /prec3/
   Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
   !c      SAVE /prec4/
   Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
   !c      SAVE /prec5/
   Common /prec6/etas(maxptn), raps(maxptn), taus(maxptn)
   !c      SAVE /prec6/
   Common /aurec1/jxa, jya, jza
   !c      SAVE /aurec1/
   Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
   !c      SAVE /aurec2/
   Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
   !c      SAVE /ilist1/
   Common /ilist2/icell, icel(10, 10, 10)
   !c      SAVE /ilist2/
   Common /ilist3/size1, size2, size3, v1, v2, v3, size
   !c      SAVE /ilist3/
   Common /ilist4/ifmpt, ichkpt, indx(maxptn)
   !c      SAVE /ilist4/
   !     6/07/02 initialize in ftime to expedite compiling:
   !        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
   !c      SAVE /ilist5/
   Common /ilist6/t, iopern, icolln
   !c      SAVE /ilist6/
   Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
   !c      SAVE /ilist7/
   Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
   !c      SAVE /ilist8/
   Common /rndm1/number
   !c      SAVE /rndm1/
   Common /rndm2/iff
   !c      SAVE /rndm2/
   Common /rndm3/iseedp
   !c      SAVE /rndm3/
   Common /ana1/ts(12)
   !c      SAVE /ana1/
   Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
   !c      SAVE /ana2/
   Common /ana3/em(4, 4, 12)
   !c      SAVE /ana3/
   Common /ana4/fdetdy(24), fdndy(24), fdndpt(12)
   !c      SAVE /ana4/
   Save
   Data centy/0D0/
   !     6/07/02 initialize in ftime to expedite compiling:
   !        data (ct(i), i = 1, MAXPTN)/MAXPTN*0d0/
   !        data (ot(i), i = 1, MAXPTN)/MAXPTN*0d0/
   !        data tlarge/1000000.d0/
   Data number/0/
   Data ts/0.11D0, 0.12D0, 0.15D0, 0.2D0, 0.3D0, 0.4D0, 0.6D0, 0.8D0, 1D0, 2D0, 4D0, 6D0/
   !
End Block Data

!*****************************************************************************
!*****************************************************************************

Subroutine inizpc

  Implicit Double Precision (A-H, O-Z)
  Save

  Call readpa

  Call inipar

  Call inian1

  Return
End Subroutine inizpc

Subroutine readpa

  Implicit Double Precision (A-H, O-Z)

  External ran1

  Character *50 str

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  !c      SAVE /para4/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /para7/ioscar, nsmbbbar, nsmmeson
  !c      SAVE /para7/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /rndm1/number
  !c      SAVE /rndm1/
  Common /rndm2/iff
  !c      SAVE /rndm2/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  iseed = iseedp
  !       this is the initialization file containing the initial values of
  !          the parameters
  !bz1/31/99
  !        open (5, file = 'zpc.ini', status = 'unknown')
  !bz1/31/99end

  !       this is the final data file containing general info about the cascade
  !bz1/31/99
  !        open (6, file = 'zpc.res', status = 'unknown')
  Open (25, File='ana/zpc.res', Status='unknown')
  !bz1/31/99end

  !       this is the input file containing initial particle records
  !bz1/25/99
  !        open (7, file = 'zpc.inp', status = 'unknown')
  !bz1/25/99end

  !       this gives the optional OSCAR standard output
  !bz1/31/99
  !        open (8, file = 'zpc.oscar', status = 'unknown')
  If (ioscar==1) Then
     Open (26, File='ana/parton.oscar', Status='unknown')
     Open (19, File='ana/hadron.oscar', Status='unknown')
  End If
  !bz1/31/99end

  !     2/11/03 combine zpc initialization into ampt.ini:
  !        open (29, file = 'zpc.ini', status = 'unknown')
  !        read (29, *) str, xmp
  xmp = 0D0
  !        read (29, *) str, xmu
  !        read (29, *) str, alpha
  cutof2 = 4.5D0*(alpha/xmu)**2
  !        read (29, *) str, rscut2
  rscut2 = 0.01D0
  !        read (29, *) str, nsevt
  nsevt = 1
  !        read (29, *) str, nevnt
  nevnt = 1
  !        read (29, *) str, nsbrun
  nsbrun = 1
  !        read (29, *) str, iftflg
  iftflg = 0
  !        read (29, *) str, ireflg
  ireflg = 1
  !bz1/31/99
  If (ireflg==0) Then
     Open (27, File='zpc.inp', Status='UNKNOWN')
  End If
  !bz1/31/99end
  !        read (29, *) str, igeflg
  igeflg = 0
  !        read (29, *) str, ibstfg
  ibstfg = 0
  !        read (29, *) str, iconfg
  iconfg = 1
  !        read (29, *) str, iordsc
  iordsc = 11
  !        read (29, *) str, ioscar
  !        read (29, *) str, v1, v2, v3
  v1 = 0.2D0
  v2 = 0.2D0
  v3 = 0.2D0
  !        read (29, *) str, size1, size2, size3
  size1 = 1.5D0
  size2 = 1.5D0
  size3 = 0.7D0
  If (size1==0D0 .Or. size2==0D0 .Or. size3==0D0) Then
     If (size1/=0D0 .Or. size2/=0D0 .Or. size3/=0D0 .Or. v1/=0D0 .Or. v2/=0D0 .Or. v3/=0D0) Then
        Print *, 'to get rid of space division:'
        Print *, 'set all sizes and vs to 0'
        Stop 'chker'
     End If
  End If
  size = min(size1, size2, size3)
  !        read (29, *) str, iff
  iff = -1
  !        read (29, *) str, iseed

  !     10/24/02 get rid of argument usage mismatch in ran1():
  isedng = -iseed
  !        a = ran1(-iseed)
  a = ran1(isedng)
  !        read (29, *) str, irused
  irused = 2
  Do i = 1, irused - 1
     !           a = ran1(2)
     iseed2 = 2
     a = ran1(iseed2)
  End Do
  !     10/24/02-end

  If (iconfg==2 .Or. iconfg==3) Then
     v1 = 0D0
     v2 = 0D0
  End If

  If (iconfg==4 .Or. iconfg==5) Then
     v1 = 0D0
     v2 = 0D0
     v3 = 0D0
  End If

  Close (5)

  Return
End Subroutine readpa

Subroutine inipar

  Implicit Double Precision (A-H, O-Z)

  Common /para4/iftflg, ireflg, igeflg, ibstfg
  !c      SAVE /para4/
  Common /para6/centy
  !c      SAVE /para6/
  Save

  If (ibstfg/=0) Then
     centy = -6D0
  End If

  Return
End Subroutine inipar

Subroutine inian1

  Implicit Double Precision (A-H, O-Z)

  Common /para4/iftflg, ireflg, igeflg, ibstfg
  !c      SAVE /para4/
  Common /ana1/ts(12)
  !c      SAVE /ana1/
  Save
  If (ibstfg/=0) Then
     a = cosh(6D0)
     Do i = 1, 12
        ts(i) = ts(i)*a
     End Do
  End If

  Return
End Subroutine inian1

!*****************************************************************************

Subroutine inievt

  Implicit Double Precision (A-H, O-Z)

  Common /para1/mul
  !c      SAVE /para1/
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  !c      SAVE /para4/
  Save

  !bz1/25/99
  !        mul = 0
  !bz1/25/99
  If (ireflg==0) Call readi
  If (igeflg/=0) Call genei
  If (ibstfg/=0) Call boosti

  Return
End Subroutine inievt

Subroutine readi

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Double Precision field(9)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  !c      SAVE /prec1/
  Save
  Do i = 1, maxptn
     If (ievt/=1 .And. i==1) Then
        ityp0(i) = ntyp
        gx0(1) = field(1)
        gy0(1) = field(2)
        gz0(1) = field(3)
        ft0(1) = field(4)
        px0(1) = field(5)
        py0(1) = field(6)
        pz0(1) = field(7)
        e0(1) = field(8)
        xmass0(i) = field(9)
        mul = 1
     Else
900     Read (27, *, End=1000) neve, ntyp, field
        If (neve<nsevt) Goto 900
        If (neve>nsevt+ievt-1) Goto 1000
        ityp0(i) = ntyp
        gx0(i) = field(1)
        gy0(i) = field(2)
        gz0(i) = field(3)
        ft0(i) = field(4)
        px0(i) = field(5)
        py0(i) = field(6)
        pz0(i) = field(7)
        e0(i) = field(8)
        xmass0(i) = field(9)
        mul = mul + 1
     End If
  End Do

1000 Continue

  Return
End Subroutine readi

Subroutine genei

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  !c      SAVE /prec1/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /lor/enenew, pxnew, pynew, pznew
  !c      SAVE /lor/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  External ran1

  iseed = iseedp
  incmul = 4000
  temp = 0.5D0
  etamin = -5D0
  etamax = 5D0
  r0 = 5D0
  tau0 = 0.1D0
  deta = etamax - etamin

  Do i = mul + 1, mul + incmul
     ityp0(i) = 21
     xmass0(i) = xmp
     Call energy(e, temp)
     Call momntm(px, py, pz, e)
     !     7/20/01:
     !           e = sqrt(e ** 2 + xmp ** 2)
     e = dsqrt(e**2+xmp**2)
     If (iconfg<=3) Then
        eta(i) = etamin + deta*ran1(iseed)
        bex = 0D0
        bey = 0D0
        bez = -tanh(eta(i))
        Call lorenz(e, px, py, pz, bex, bey, bez)
        px0(i) = pxnew
        py0(i) = pynew
        pz0(i) = pznew
        e0(i) = enenew
     Else
        px0(i) = px
        py0(i) = py
        pz0(i) = pz
        e0(i) = e
     End If
  End Do

  Do i = mul + 1, mul + incmul
     If (iconfg<=3) Then
        gz0(i) = tau0*sinh(eta(i))
        ft0(i) = tau0*cosh(eta(i))
        If (iconfg==1) Then
           Call posit1(x, y, r0)
           gx0(i) = x + px0(i)*ft0(i)/e0(i)
           gy0(i) = y + py0(i)*ft0(i)/e0(i)
        Else If (iconfg==2 .Or. iconfg==3) Then
           Call posit2(x, y)
           gx0(i) = x
           gy0(i) = y
        End If
     Else
        ft0(i) = 0D0
        Call posit3(x, y, z)
        gx0(i) = x
        gy0(i) = y
        gz0(i) = z
     End If
  End Do

  mul = mul + incmul

  !       check if it's necessary to adjust array size 'adarr'
  If (mul>=maxptn .Or. mul==0) Then
     Print *, 'event', ievt, 'has', mul, 'number of gluon', 'adjusting counting is necessary'
     Stop 'adarr'
  End If

  Return
End Subroutine genei

Subroutine posit1(x, y, r0)

  Implicit Double Precision (A-H, O-Z)

  External ran1
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  iseed = iseedp
10 x = 2D0*ran1(iseed) - 1D0
  y = 2D0*ran1(iseed) - 1D0
  If (x**2+y**2>1D0) Goto 10
  x = x*r0
  y = y*r0

  Return
End Subroutine posit1

Subroutine posit2(x, y)

  Implicit Double Precision (A-H, O-Z)

  External ran1

  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save
  iseed = iseedp
  x = 2D0*ran1(iseed) - 1D0
  y = 2D0*ran1(iseed) - 1D0
  x = x*5D0*size1
  y = y*5D0*size2

  Return
End Subroutine posit2

Subroutine posit3(x, y, z)

  Implicit Double Precision (A-H, O-Z)

  External ran1

  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  iseed = iseedp
  x = 2D0*ran1(iseed) - 1D0
  y = 2D0*ran1(iseed) - 1D0
  z = 2D0*ran1(iseed) - 1D0
  x = x*5D0*size1
  y = y*5D0*size2
  z = z*5D0*size3

  Return
End Subroutine posit3

Subroutine energy(e, temp)

  !       to generate the magnitude of the momentum e,
  !       knowing the temperature of the local thermal distribution temp

  Implicit Double Precision (A-H, O-Z)

  External ran1

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  iseed = iseedp
1000 Continue

  e = ran1(iseed)
  e = e*ran1(iseed)
  e = e*ran1(iseed)

  If (e<=0D0) Goto 1000
  e = -temp*log(e)
  If (ran1(iseed)>exp((e-dsqrt(e**2+xmp**2))/temp)) Then
     Goto 1000
  End If

  Return
End Subroutine energy

Subroutine momntm(px, py, pz, e)

  !       to generate the 3 components of the momentum px, py, pz,
  !       from the magnitude of the momentum e

  Implicit Double Precision (A-H, O-Z)

  External ran1

  Parameter (pi=3.14159265358979D0)
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  iseed = iseedp
  cost = 2D0*ran1(iseed) - 1D0
  !     7/20/01:
  !        sint = sqrt(1d0 - cost ** 2)
  sint = dsqrt(1D0-cost**2)
  phi = 2D0*pi*ran1(iseed)

  px = e*sint*cos(phi)
  py = e*sint*sin(phi)
  pz = e*cost

  Return
End Subroutine momntm

Subroutine boosti

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para6/centy
  !c      SAVE /para6/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  !c      SAVE /prec1/
  Common /lor/enenew, pxnew, pynew, pznew
  !c      SAVE /lor/
  Save

  External lorenz

  bex = 0D0
  bey = 0D0
  bez = -tanh(centy)

  !       save data for many runs of the same initial condition
  Do i = 1, mul
     px1 = gx0(i)
     py1 = gy0(i)
     pz1 = gz0(i)
     e1 = ft0(i)
     Call lorenz(e1, px1, py1, pz1, bex, bey, bez)
     gx0(i) = pxnew
     gy0(i) = pynew
     gz0(i) = pznew
     ft0(i) = enenew
     px1 = px0(i)
     py1 = py0(i)
     pz1 = pz0(i)
     e1 = e0(i)
     Call lorenz(e1, px1, py1, pz1, bex, bey, bez)
     px0(i) = pxnew
     py0(i) = pynew
     pz0(i) = pznew
     e0(i) = enenew
  End Do

  Return
End Subroutine boosti

!*****************************************************************************

Subroutine inirun
  Save

  !       sort prec2 according to increasing formation time
  Call ftime
  Call inirec
  Call iilist
  Call inian2

  Return
End Subroutine inirun

Subroutine ftime
  !       this subroutine generates formation time for the particles
  !       indexing ft(i)
  !       input e(i)
  !       output ft(i), indx(i)

  Implicit Double Precision (A-H, O-Z)

  External ftime1
  Parameter (maxptn=400001)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  !c      SAVE /para4/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  !c      SAVE /prec1/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Common /par1/formt
  !c      SAVE /par1/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  iseed = iseedp
  !lin-6/07/02 initialize here to expedite compiling, instead in zpcbdt:
  Do i = 1, maxptn
     ct(i) = 0D0
     ot(i) = 0D0
  End Do
  tlarge = 1000000.D0
  !lin-6/07/02-end

  If (iftflg==0) Then
     !     5/01/01 different prescription for parton initial formation time:
     If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
        Do i = 1, mul
           If (ft0(i)>tlarge) ft0(i) = tlarge
        End Do
        Goto 150
     Else
        !     5/01/01-end

        Do i = 1, maxptn
           ft0(i) = tlarge
        End Do
        Do i = 1, mul
           xmt2 = px0(i)**2 + py0(i)**2 + xmp**2
           formt = xmt2/e0(i)
           ft0(i) = ftime1(iseed)
           If (ft0(i)>tlarge) ft0(i) = tlarge
        End Do
        !     5/01/01:
     End If

  End If

  !     5/01/01:
150 Continue

  !        call index1(MAXPTN, mul, ft0, indx)
  If (mul>1) Then
     Call index1(maxptn, mul, ft0, indx)
  Else
     !lin-7/09/03: need to set value for mul=1:
     indx(1) = 1
  End If
  !
  Return
End Subroutine ftime

Subroutine inirec

  Implicit Double Precision (A-H, O-Z)
  External ran1
  Parameter (maxptn=400001)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para4/iftflg, ireflg, igeflg, ibstfg
  !c      SAVE /para4/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  !c      SAVE /prec1/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec3/gxs(maxptn), gys(maxptn), gzs(maxptn), fts(maxptn), pxs(maxptn), pys(maxptn), pzs(maxptn), es(maxptn), xmasss(maxptn), ityps(maxptn)
  !c      SAVE /prec3/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /prec6/etas(maxptn), raps(maxptn), taus(maxptn)
  !c      SAVE /prec6/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  !bz1/25/99
  Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
  !c      SAVE /ilist7/
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
  !c      SAVE /ilist8/
  !bz1/25/99end
  Common /smearz/smearp, smearh
  !c      SAVE /smearz/
  !lin-8/2015:
  !        dimension vxp(MAXPTN), vyp(MAXPTN), vzp(MAXPTN)
  Common /precpb/vxp(maxptn), vyp(maxptn), vzp(maxptn)
  !lin-8/2015:
  Common /precpa/vxp0(maxptn), vyp0(maxptn), vzp0(maxptn), xstrg0(maxptn), ystrg0(maxptn), xstrg(maxptn), ystrg(maxptn), istrg0(maxptn), istrg(maxptn)
  !        common /precpa/ vxp0(MAXPTN), vyp0(MAXPTN), vzp0(MAXPTN)
  !c      SAVE /precpa/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  !lin-6/06/02 local parton freezeout:
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  !c      SAVE /frzprc/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /arevt/iaevt, iarun, miss
  Save
  iseed = iseedp
  !lin-6/06/02 local freezeout initialization:
  If (isoft==5) Then
     itlast = 0
     Call inifrz
  End If

  Do i = 1, mul
     !lin-7/09/01 define indx(i) to save time:
     !           ityp(i) = ityp0(indx(i))
     !           gx(i) = gx0(indx(i))
     !           gy(i) = gy0(indx(i))
     !           gz(i) = gz0(indx(i))
     !           ft(i) = ft0(indx(i))
     !           px(i) = px0(indx(i))
     !           py(i) = py0(indx(i))
     !           pz(i) = pz0(indx(i))
     !           e(i) = e0(indx(i))
     !           xmass(i) = xmass0(indx(i))
     !cbz1/25/99
     !           LSTRG1(I) = LSTRG0(INDX(I))
     !           LPART1(I) = LPART0(INDX(I))
     !cbz1/25/99end
     indxi = indx(i)
     ityp(i) = ityp0(indxi)
     gx(i) = gx0(indxi)
     gy(i) = gy0(indxi)
     gz(i) = gz0(indxi)
     ft(i) = ft0(indxi)
     px(i) = px0(indxi)
     py(i) = py0(indxi)
     pz(i) = pz0(indxi)
     e(i) = e0(indxi)
     xmass(i) = xmass0(indxi)
     lstrg1(i) = lstrg0(indxi)
     lpart1(i) = lpart0(indxi)
     vxp(i) = vxp0(indxi)
     vyp(i) = vyp0(indxi)
     vzp(i) = vzp0(indxi)
     !lin-8/2015:
     xstrg0(i) = xstrg(indxi)
     ystrg0(i) = ystrg(indxi)
     istrg0(i) = istrg(indxi)
     !lin-7/09/01-end
     !
     !lin-6/06/02 local freezeout initialization:
     If (isoft==5) Then
        idfrz(i) = ityp(i)
        gxfrz(i) = gx(i)
        gyfrz(i) = gy(i)
        gzfrz(i) = gz(i)
        ftfrz(i) = ft(i)
        pxfrz(i) = px(i)
        pyfrz(i) = py(i)
        pzfrz(i) = pz(i)
        efrz(i) = e(i)
        xmfrz(i) = xmass(i)
        ifrz(i) = 0
     End If
     !lin-6/06/02-end
  End Do

  !       save particle info for fixed time analysis
  Do i = 1, mul
     ityps(i) = ityp(i)
     gxs(i) = gx(i)
     gys(i) = gy(i)
     gzs(i) = gz(i)
     fts(i) = ft(i)
     pxs(i) = px(i)
     pys(i) = py(i)
     pzs(i) = pz(i)
     es(i) = e(i)
     xmasss(i) = xmass(i)
  End Do

  !lin-6/2009
  If (isoft==1 .And. (ioscar==2 .Or. ioscar==3)) Write (92, *) iaevt, miss, mul

  Do i = 1, mul
     energy = e(i)
     vx(i) = px(i)/energy
     vy(i) = py(i)/energy
     vz(i) = pz(i)/energy
     If (iftflg==0) Then
        formt = ft(i)
        !     7/09/01 propagate partons with parent velocity till formation
        !     so that partons in same hadron have 0 distance:
        !            gx(i) = gx(i) + vx(i) * formt
        !            gy(i) = gy(i) + vy(i) * formt
        !            gz(i) = gz(i) + vz(i) * formt
        If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
           gx(i) = gx(i) + vxp(i)*formt
           gy(i) = gy(i) + vyp(i)*formt
           gz(i) = gz(i) + vzp(i)*formt
        Else
           gx(i) = gx(i) + vx(i)*formt
           gy(i) = gy(i) + vy(i)*formt
           gz(i) = gz(i) + vz(i)*formt
        End If
        !     7/09/01-end
        !
        !     3/27/00-ctest off no smear z on partons to avoid eta overflow:
        !              gz(i) = gz(i)+smearp*(2d0 * ran1(iseed) - 1d0)
        !     to give eta=y +- smearp*random:
        !              smeary=smearp*(2d0 * ran1(iseed) - 1d0)
        !              smearf=dexp(2*smeary)*(1+vz(i))/(1-vz(i)+1.d-8)
        !              gz(i) = gz(i)+formt*(smearf-1)/(smearf+1)
        !     3/27/00-end
     End If

     !lin-6/2009 write out initial parton information after string melting
     !     and after propagating to its format time:
     If (ioscar==2 .Or. ioscar==3) Then
        If (dmax1(abs(gx(i)),abs(gy(i)),abs(gz(i)),abs(ft(i)))<9999) Then
           !lin-8/2015:
           Write (92, 200) ityp(i), px(i), py(i), pz(i), xmass(i), gx(i), gy(i), gz(i), ft(i), istrg0(i), xstrg0(i), ystrg0(i)
        Else
           !lin-8/2015:
           Write (92, 201) ityp(i), px(i), py(i), pz(i), xmass(i), gx(i), gy(i), gz(i), ft(i), istrg0(i), xstrg0(i), ystrg0(i)
        End If
     End If
     !
  End Do

  If (iconfg<=3) Then
     Do i = 1, mul
        If (ft(i)<=abs(gz(i))) Then
           eta(i) = 1000000.D0
        Else
           eta(i) = 0.5D0*log((ft(i)+gz(i))/(ft(i)-gz(i)))
        End If
        If (e(i)<=abs(pz(i))) Then
           rap(i) = 1000000.D0
        Else
           rap(i) = 0.5D0*log((e(i)+pz(i))/(e(i)-pz(i)))
        End If
        !lin-8/2015 to avoid IEEE_OVERFLOW_FLAG:
        !              tau(i) = ft(i) / cosh(eta(i))
        If (eta(i)<1000000.D0) Then
           tau(i) = ft(i)/cosh(eta(i))
        Else
           tau(i) = 1D-10
        End If
        !
     End Do

     Do i = 1, mul
        etas(i) = eta(i)
        raps(i) = rap(i)
        taus(i) = tau(i)
     End Do
  End If

  Return
  !lin-8/2015:
  ! 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
  ! 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
  !     reduce file size:
  ! 200       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f9.3),
  !     1          1x,I6,2(1x,f8.3))
  ! 201       format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e9.3),
  !     1          1x,I6,2(1x,f8.3))
200 Format (I3, 2(1X,F7.2), 1X, F8.2, 1X, F6.3, 4(1X,F8.2), 1X, I5, 2(1X,F7.2))
201 Format (I3, 2(1X,F7.2), 1X, F8.2, 1X, F6.3, 4(1X,E8.2), 1X, I5, 2(1X,F7.2))
End Subroutine inirec

Subroutine iilist

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  !c      SAVE /para1/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist2/icell, icel(10, 10, 10)
  !c      SAVE /ilist2/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Common /ilist6/t, iopern, icolln
  !c      SAVE /ilist6/
  Save

  iscat = maxptn
  jscat = maxptn

  Do i = 1, mul
     next(i) = 0
     last(i) = 0
     icsta(i) = 0
     nic(i) = 0
     icels(i) = 0
  End Do

  icell = 0
  Do i1 = 1, 10
     Do i2 = 1, 10
        Do i3 = 1, 10
           icel(i1, i2, i3) = 0
        End Do
     End Do
  End Do

  ichkpt = 0
  ifmpt = 1

  Do i = 1, mul
     ct(i) = tlarge
     ot(i) = tlarge
  End Do

  iopern = 0
  icolln = 0
  t = 0.D0

  Return
End Subroutine iilist

Subroutine inian2

  Implicit Double Precision (A-H, O-Z)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  !c      SAVE /ana2/
  Save

  If (iconfg<=3) Then
     Do i = 1, 12
        det(i) = 0D0
        dn(i) = 0D0
        det1(i) = 0D0
        dn1(i) = 0D0
        det2(i) = 0D0
        dn2(i) = 0D0
     End Do
  End If

  Return
End Subroutine inian2

!*****************************************************************************
!*****************************************************************************

Subroutine zpcrun(*)

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Parameter (tend1=250D0)
  Parameter (tend2=6.1D0)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para5/iconfg, iordsc
  Common /para7/ioscar, nsmbbbar, nsmmeson
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Common /ilist6/t, iopern, icolln
  !c      SAVE /ilist6/
  Common /ana1/ts(12)
  !c      SAVE /ana1/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  Common /arevt/iaevt, iarun, miss
  Save

  !       save last collision info
  If (mod(ictype,2)==0) Then
     Call savrec(iscat)
     Call savrec(jscat)
  End If

  !1      get operation type
  Call getict(t1)
  !2      check freezeout condition
  If (iconfg==1 .And. t1>tlarge/2D0) Return 1
  If (iconfg==2 .Or. iconfg==3) Then
     If (t1>tend1) Return 1
     !           if (ichkpt .eq. mul) then
     !              ii = 0
     !              do i = 1, mul
     !                 gztemp = gz(i) + vz(i) * (t1 - ft(i))
     !                 if (sqrt(t1 ** 2 - gztemp ** 2) .lt. tend) then
     !                    ii = 1
     !                    goto 1000
     !                 end if
     !              end do
     ! 1000              continue
     !              if (ii .eq. 0) return 1
     !           end if
  End If
  If (iconfg==4 .Or. iconfg==5) Then
     If (t1>tend2) Return 1
  End If

  !lin-6/06/02 local freezeout for string melting,
  !     decide what partons have frozen out at time t1:
  If (isoft==5) Then
     Call local(t1)
  End If

  !3      update iopern, t

  iopern = iopern + 1
  t = t1
  If (mod(ictype,2)==0) Then
     icolln = icolln + 1

     !     4/18/01-ctest off
     !           write (2006, 1233) 'iscat=', iscat, 'jscat=', jscat,
     !           write (2006, *) 'iscat=', iscat, ' jscat=', jscat,
     !     1 ityp(iscat), ityp(jscat)
     !           write (2006, 1233) 'iscat=', max(indx(iscat), indx(jscat)),
     !     &        'jscat=', min(indx(iscat), indx(jscat))

     !           write (2006, 1234) ' icolln=', icolln, 't=', t

     ! 1233           format (a10, i10, a10, i10)
     ! 1234           format (a15, i10, a5, f23.17, a5, f23.17)
  End If

  !4.1    deal with formation
  If (iconfg==1 .Or. iconfg==2 .Or. iconfg==4) Then
     If (ictype==1 .Or. ictype==2 .Or. ictype==5 .Or. ictype==6) Then
        Call celasn
     End If
  End If

  !4.2    deal with collisions

  If (ictype/=1) Then

     iscat0 = iscat
     jscat0 = jscat

     !        iscat is the larger one so that if it's a wall collision,
     !       it's still ok
     iscat = max0(iscat0, jscat0)
     jscat = min0(iscat0, jscat0)

     !test off check icsta(i): 0 with f77 compiler
     !        write(9,*) 'BB:ictype,t1,iscat,jscat,icsta(i)=',
     !     1 ictype,t1,iscat,jscat,icsta(iscat)

     !       check collision time table error 'tterr'
     !lin-4/2008 to avoid out-of-bound error in next():
     !           if (jscat .ne. 0 .and. next(jscat) .ne. iscat)
     !     &        then
     !              print *, 'iscat=', iscat, 'jscat=', jscat,
     !     &             'next(', jscat, ')=', next(jscat)
     !
     !              if (ct(iscat) .lt. tlarge / 2d0) stop 'tterr'
     !              if (ct(jscat) .lt. tlarge / 2d0) stop 'tterr'
     !           end if
     If (jscat/=0) Then
        If (next(jscat)/=iscat) Then
           Print *, 'iscat=', iscat, 'jscat=', jscat, 'next(', jscat, ')=', next(jscat)
           If (ct(iscat)<tlarge/2D0) Stop 'tterr'
           If (ct(jscat)<tlarge/2D0) Stop 'tterr'
        End If
     End If
     !lin-4/2008-end

     !4.2.1     collisions with wall

     !     8/19/02 avoid actual argument in common blocks of cellre:
     niscat = iscat
     njscat = jscat
     !           if (icsta(iscat) .ne. 0) call cellre(iscat, t)
     !           if (jscat .ne. 0) then
     !              if (icsta(jscat) .ne. 0) call cellre(jscat, t)
     !           end if
     If (icsta(iscat)/=0) Call cellre(niscat, t)
     If (jscat/=0) Then
        If (icsta(jscat)/=0) Call cellre(njscat, t)
     End If

     !4.2.2     collision between particles

     !lin-6/2009 write out info for each collision:
     !           if (mod(ictype, 2) .eq. 0) call scat(t, iscat, jscat)
     If (mod(ictype,2)==0) Then
        If (ioscar==3) Then
           Write (95, *) 'event,miss,iscat,jscat=', iaevt, miss, iscat, jscat
           If (dmax1(abs(gx(iscat)),abs(gy(iscat)),abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))<9999) Then
              Write (95, 200) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 200) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           Else
              Write (95, 201) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 201) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           End If
        End If
        !
        Call scat(t, iscat, jscat)
        !
        If (ioscar==3) Then
           If (dmax1(abs(gx(iscat)),abs(gy(iscat)),abs(gz(iscat)),abs(ft(iscat)),abs(gx(jscat)),abs(gy(jscat)),abs(gz(jscat)),abs(ft(jscat)))<9999) Then
              Write (95, 200) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 200) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           Else
              Write (95, 201) ityp(iscat), px(iscat), py(iscat), pz(iscat), xmass(iscat), gx(iscat), gy(iscat), gz(iscat), ft(iscat)
              Write (95, 201) ityp(jscat), px(jscat), py(jscat), pz(jscat), xmass(jscat), gx(jscat), gy(jscat), gz(jscat), ft(jscat)
           End If
        End If
     End If

  End If

  !5      update the interaction list
  Call ulist(t)

  !6      update ifmpt. ichkpt
  !       old ichkpt and ifmpt are more conveniently used in ulist
  If (ifmpt<=mul) Then
     If (ictype/=0 .And. ictype/=3 .And. ictype/=4) Then
        ichkpt = ichkpt + 1
        ifmpt = ifmpt + 1
     End If
  End If

  Return
200 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,F8.2))
201 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,E8.2))
End Subroutine zpcrun

Subroutine savrec(i)

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec3/gxs(maxptn), gys(maxptn), gzs(maxptn), fts(maxptn), pxs(maxptn), pys(maxptn), pzs(maxptn), es(maxptn), xmasss(maxptn), ityps(maxptn)
  !c      SAVE /prec3/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /prec6/etas(maxptn), raps(maxptn), taus(maxptn)
  !c      SAVE /prec6/
  Save

  ityps(i) = ityp(i)
  gxs(i) = gx(i)
  gys(i) = gy(i)
  gzs(i) = gz(i)
  fts(i) = ft(i)
  pxs(i) = px(i)
  pys(i) = py(i)
  pzs(i) = pz(i)
  es(i) = e(i)
  xmasss(i) = xmass(i)
  etas(i) = eta(i)
  raps(i) = rap(i)
  taus(i) = tau(i)

  Return
End Subroutine savrec

Subroutine getict(t1)
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  !c      SAVE /para1/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  !       neglect possibility of 2 collisions at the same time
  !0       set initial conditions

  t1 = tlarge
  iscat = 0
  jscat = 0

  !1      get next collision between particles
  Do i = 1, ichkpt
     If (ot(i)<t1) Then
        t1 = ot(i)
        iscat = i
     End If
  End Do
  If (iscat/=0) jscat = next(iscat)

  !2      get ictype
  !     10/30/02 ictype=0:collision; 1:parton formation
  If (iscat/=0 .And. jscat/=0) Then
     If (icsta(iscat)==0 .And. icsta(jscat)==0) Then
        ictype = 0
     Else
        ictype = 4
     End If
  Else If (iscat/=0 .Or. jscat/=0) Then
     ictype = 3
  End If
  !
  If (ifmpt<=mul) Then
     If (ft(ifmpt)<t1) Then
        ictype = 1
        t1 = ft(ifmpt)
     Else If (ft(ifmpt)==t1) Then
        If (ictype==0) ictype = 2
        If (ictype==3) ictype = 5
        If (ictype==4) ictype = 6
     End If
  End If

  Return
End Subroutine getict

Subroutine celasn
  !       this subroutine is used to assign a cell for a newly formed particle
  !       output: nic(MAXPTN) icels(MAXPTN) in the common /ilist1/
  !       icell, and icel(10,10,10) in the common /ilist2/

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist2/icell, icel(10, 10, 10)
  !c      SAVE /ilist2/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Save

  External integ

  i = ifmpt
  tt = ft(i)
  td = tt - size
  If (iconfg==1 .And. (size1==0D0 .Or. size2==0D0 .Or. size3==0D0)) Then
     i1 = 11
     i2 = 11
     i3 = 11
  Else If (iconfg==4 .Or. td<=0D0) Then
     i1 = integ(gx(i)/size1) + 6
     i2 = integ(gy(i)/size2) + 6
     i3 = integ(gz(i)/size3) + 6
     If (integ(gx(i)/size1)==gx(i)/size1 .And. vx(i)<0D0) i1 = i1 - 1
     If (integ(gy(i)/size2)==gy(i)/size2 .And. vy(i)<0D0) i2 = i2 - 1
     If (integ(gz(i)/size3)==gz(i)/size3 .And. vz(i)<0D0) i3 = i3 - 1
  Else
     i1 = integ(gx(i)/(size1+v1*td)) + 6
     i2 = integ(gy(i)/(size2+v2*td)) + 6
     i3 = integ(gz(i)/(size3+v3*td)) + 6
     If (integ(gx(i)/(size1+v1*td))==gx(i)/(size1+v1*td) .And. vx(i)<(i1-6)*v1) i1 = i1 - 1
     If (integ(gy(i)/(size2+v2*td))==gy(i)/(size2+v2*td) .And. vy(i)<(i2-6)*v2) i2 = i2 - 1
     If (integ(gz(i)/(size3+v3*td))==gz(i)/(size3+v3*td) .And. vz(i)<(i3-6)*v3) i3 = i3 - 1
  End If

  If (i1<=0 .Or. i1>=11 .Or. i2<=0 .Or. i2>=11 .Or. i3<=0 .Or. i3>=11) Then
     i1 = 11
     i2 = 11
     i3 = 11
  End If

  If (i1==11) Then
     j = icell
     Call newcre(i, j)
     icell = j
     icels(i) = 111111
  Else
     j = icel(i1, i2, i3)
     Call newcre(i, j)
     icel(i1, i2, i3) = j
     icels(i) = i1*10000 + i2*100 + i3
  End If

  Return
End Subroutine celasn

Integer Function integ(x)
  !       this function is used to get the largest integer that is smaller than
  !       x

  Implicit Double Precision (A-H, O-Z)
  Save

  If (x<0D0) Then
     integ = int(x-1D0)
  Else
     integ = int(x)
  End If

  Return
End Function integ

Subroutine cellre(i, t)
  !       this subroutine is used for changing the cell of a particle that
  !       collide with the wall

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist2/icell, icel(10, 10, 10)
  !c      SAVE /ilist2/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical good

  External integ

  !       this happens before update the /prec2/ common; in contrast with
  !       scat which happens after updating the glue common

  t0 = t

1000 Continue

  If (iconfg==3 .Or. iconfg==5) Then
     k = mod(icsta(i), 10)

     If (k==1) Then
        gx(i) = gx(i) - 10D0*size1
        dgxa(i) = dgxa(i) + 10D0*size1
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgxa(ii) = dgxa(ii) - 10D0*size1
           End If
        End Do
     End If
     If (k==2) Then
        gx(i) = gx(i) + 10D0*size1
        dgxa(i) = dgxa(i) - 10D0*size1
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgxa(ii) = dgxa(ii) + 10D0*size1
           End If
        End Do
     End If
     If (k==3) Then
        gy(i) = gy(i) - 10D0*size2
        dgya(i) = dgya(i) + 10D0*size2
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgya(ii) = dgya(ii) - 10D0*size2
           End If
        End Do
     End If
     If (k==4) Then
        gy(i) = gy(i) + 10D0*size2
        dgya(i) = dgya(i) - 10D0*size2
        Do ii = 1, ichkpt
           If (next(ii)==i) Then
              dgya(ii) = dgya(ii) + 10D0*size2
           End If
        End Do
     End If
     If (iconfg==5) Then
        If (k==5) Then
           gz(i) = gz(i) - 10D0*size3
           dgza(i) = dgza(i) + 10D0*size3
           Do ii = 1, ichkpt
              If (next(ii)==i) Then
                 dgza(ii) = dgza(ii) - 10D0*size3
              End If
           End Do
        End If
        If (k==6) Then
           gz(i) = gz(i) + 10D0*size3
           dgza(i) = dgza(i) - 10D0*size3
           Do ii = 1, ichkpt
              If (next(ii)==i) Then
                 dgza(ii) = dgza(ii) + 10D0*size3
              End If
           End Do
        End If
     End If
  Else
     icels0 = icels(i)

     i1 = icels0/10000
     i2 = (icels0-i1*10000)/100
     i3 = icels0 - i1*10000 - i2*100

     !c       for particle inside the cube
     If (i1>=1 .And. i1<=10 .And. i2>=1 .And. i2<=10 .And. i3>=1 .And. i3<=10) Then

        !       this assignment takes care of nic(i)=0 automatically
        If (icel(i1,i2,i3)==i) icel(i1, i2, i3) = nic(i)

        !1      rearrange the old cell

        Call oldcre(i)

        !2      rearrange the new cell

        k = mod(icsta(i), 10)

        !2.1    particle goes out of the cube
        If (iconfg==1) Then
           good = (i1==1 .And. k==2) .Or. (i1==10 .And. k==1) .Or. (i2==1 .And. k==4) .Or. (i2==10 .And. k==3) .Or. (i3==1 .And. k==6) .Or. (i3==10 .And. k==5)
        End If
        If (iconfg==2) Then
           good = (i3==1 .And. k==6) .Or. (i3==10 .And. k==5)
        End If
        If (good) Then

           !                j = icell
           Call newcre(i, icell)
           !                 icell = j

           icels(i) = 111111

           !2.2    particle moves inside the cube
        Else

           If (k==1) i1 = i1 + 1
           If (k==2) i1 = i1 - 1
           If (k==3) i2 = i2 + 1
           If (k==4) i2 = i2 - 1
           If (k==5) i3 = i3 + 1
           If (k==6) i3 = i3 - 1

           If (iconfg==2 .Or. iconfg==4) Then
              If (i1==0) Then
                 i1 = 10
                 gx(i) = gx(i) + 10D0*size1
              End If
              If (i1==11) Then
                 i1 = 1
                 gx(i) = gx(i) - 10D0*size1
              End If
              If (i2==0) Then
                 i2 = 10
                 gy(i) = gy(i) + 10D0*size2
              End If
              If (i2==11) Then
                 i2 = 1
                 gy(i) = gy(i) - 10D0*size2
              End If
              If (iconfg==4) Then
                 If (i3==0) Then
                    i3 = 10
                    gz(i) = gz(i) + 10D0*size3
                 End If
                 If (i3==11) Then
                    i3 = 1
                    gz(i) = gz(i) - 10D0*size3
                 End If
              End If
           End If

           j = icel(i1, i2, i3)

           Call newcre(i, j)
           !       in case icel changes

           icel(i1, i2, i3) = j

           icels(i) = i1*10000 + i2*100 + i3

        End If

        !c       for particles outside the cube
     Else

        If (icell==i) icell = nic(i)

        Call oldcre(i)

        k = mod(icsta(i), 10)

        ddt = t - ft(i)
        dtt = t - size
        If (dtt<=0D0) Then
           i1 = integ((gx(i)+vx(i)*ddt)/size1) + 6
           i2 = integ((gy(i)+vy(i)*ddt)/size2) + 6
           i3 = integ((gz(i)+vz(i)*ddt)/size3) + 6
        Else
           i1 = integ((gx(i)+vx(i)*ddt)/(size1+v1*dtt)) + 6
           i2 = integ((gy(i)+vy(i)*ddt)/(size2+v2*dtt)) + 6
           i3 = integ((gz(i)+vz(i)*ddt)/(size3+v3*dtt)) + 6
        End If


        If (k==1) i1 = 1
        If (k==2) i1 = 10
        If (k==3) i2 = 1
        If (k==4) i2 = 10
        If (k==5) i3 = 1
        If (k==6) i3 = 10

        j = icel(i1, i2, i3)
        Call newcre(i, j)
        icel(i1, i2, i3) = j

        icels(i) = i1*10000 + i2*100 + i3

     End If
  End If

  If (next(i)/=0) Then
     otmp = ot(next(i))
     ctmp = ct(next(i))
  End If

  If (i1==11 .And. i2==11 .And. i3==11) Then
     Call dchout(i, k, t)
  Else
     If (iconfg==1) Then
        Call dchin1(i, k, i1, i2, i3, t)
     Else If (iconfg==2) Then
        Call dchin2(i, k, i1, i2, i3, t)
     Else If (iconfg==4) Then
        Call dchin3(i, k, i1, i2, i3, t)
     End If
  End If

  If (icsta(i)/10==11) Then
     ot(next(i)) = otmp
     ct(next(i)) = ctmp
     next(next(i)) = i
     Call wallc(i, i1, i2, i3, t0, tmin1)
     If (tmin1<ct(i)) Then
        icsta(i) = icsta(i) + 10
        t0 = tmin1
        Goto 1000
     End If
  End If

  Return
End Subroutine cellre

Subroutine oldcre(i)
  !       this subroutine is used to rearrange the old cell nic when a particle
  !       goes out of the cell

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Save

  If (nic(i)==0) Return

  j = nic(i)

  If (nic(j)==i) Then
     nic(j) = 0
     Return
  End If

  Do While (nic(j)/=i)
     j = nic(j)
  End Do

  nic(j) = nic(i)

  Return
End Subroutine oldcre


Subroutine newcre(i, k)
  !       this subroutine is used to mk rearrange of the new cell a particle
  !       enters,
  !       input i
  !       output nic(i)

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Save

  If (k==0) Then
     k = i
     nic(i) = 0
  Else If (nic(k)==0) Then
     nic(k) = i
     nic(i) = k
  Else
     j = k
     Do While (nic(j)/=k)
        j = nic(j)
     End Do

     nic(j) = i
     nic(i) = k

  End If

  Return
End Subroutine newcre

Subroutine scat(t, iscat, jscat)

  !       this subroutine is used to calculate the 2 particle scattering

  Implicit Double Precision (A-H, O-Z)
  Save

  Call newpos(t, iscat)
  Call newpos(t, jscat)
  Call newmom(t)

  Return
End Subroutine scat

Subroutine newpos(t, i)

  !       this subroutine is used to calculate the 2 particle scattering
  !       get new position

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  dt1 = ct(i) - ft(i)

  gx(i) = gx(i) + vx(i)*dt1
  gy(i) = gy(i) + vy(i)*dt1
  gz(i) = gz(i) + vz(i)*dt1
  ft(i) = ct(i)

  If (iconfg<=3) Then
     If (ft(i)<=abs(gz(i))) Then
        eta(i) = 1000000.D0
     Else
        eta(i) = 0.5D0*log((ft(i)+gz(i))/(ft(i)-gz(i)))
     End If
     !lin-8/2015 to avoid IEEE_OVERFLOW_FLAG:
     !           tau(i) = ft(i) / cosh(eta(i))
     If (eta(i)<1000000.D0) Then
        tau(i) = ft(i)/cosh(eta(i))
     Else
        tau(i) = 1D-10
     End If
     !
  End If

  Return
End Subroutine newpos



!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine newmom(t)

  !       this subroutine is used to calculate the 2 particle scattering

  Implicit Double Precision (A-H, O-Z)

  Parameter (hbarc=0.197327054D0)
  Parameter (maxptn=400001)
  Parameter (pi=3.14159265358979D0)
  Common /para1/mul
  !c      SAVE /para1/
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  !trans
  Common /para6/centy
  !c      SAVE /para6/
  !transend
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /lor/enenew, pxnew, pynew, pznew
  !c      SAVE /lor/
  Common /cprod/xn1, xn2, xn3
  !c      SAVE /cprod/
  Common /rndm2/iff
  !c      SAVE /rndm2/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  !c      SAVE /frzprc/
  Save

  !lin-6/06/02 no momentum change for partons already frozen out,
  !     however, spatial upgrade is needed to ensure overall system freezeout:
  If (isoft==5) Then
     If (ifrz(iscat)==1 .Or. ifrz(jscat)==1) Then
        last(iscat) = jscat
        last(jscat) = iscat
        Return
     End If
  End If
  !lin-6/06/02-end

  !       iff is used to randomize the interaction to have both attractive and
  !        repulsive

  iff = -iff

  If (iconfg==2 .Or. iconfg==4) Then
     icels1 = icels(iscat)
     i1 = icels1/10000
     j1 = (icels1-i1*10000)/100
     icels2 = icels(jscat)
     i2 = icels2/10000
     j2 = (icels2-i2*10000)/100
     If (iconfg==4) Then
        k1 = icels1 - i1*10000 - j1*100
        k2 = icels2 - i2*10000 - j2*100
     End If
  End If

  px1 = px(iscat)
  py1 = py(iscat)
  pz1 = pz(iscat)
  e1 = e(iscat)
  x1 = gx(iscat)
  y1 = gy(iscat)
  z1 = gz(iscat)
  t1 = ft(iscat)
  px2 = px(jscat)
  py2 = py(jscat)
  pz2 = pz(jscat)
  e2 = e(jscat)

  If (iconfg==1) Then
     x2 = gx(jscat)
     y2 = gy(jscat)
     z2 = gz(jscat)
  Else If (iconfg==2 .Or. iconfg==4) Then
     If (i1-i2>5) Then
        x2 = gx(jscat) + 10D0*size1
     Else If (i1-i2<-5) Then
        x2 = gx(jscat) - 10D0*size1
     Else
        x2 = gx(jscat)
     End If
     If (j1-j2>5) Then
        y2 = gy(jscat) + 10D0*size2
     Else If (j1-j2<-5) Then
        y2 = gy(jscat) - 10D0*size2
     Else
        y2 = gy(jscat)
     End If
     If (iconfg==4) Then
        If (k1-k2>5) Then
           z2 = gz(jscat) + 10D0*size3
        Else If (k1-k2<-5) Then
           z2 = gz(jscat) - 10D0*size3
        Else
           z2 = gz(jscat)
        End If
     Else
        z2 = gz(jscat)
     End If
  Else If (iconfg==3 .Or. iconfg==5) Then
     x2 = gx(jscat) + dgxa(jscat)
     y2 = gy(jscat) + dgya(jscat)
     If (iconfg==5) Then
        z2 = gz(jscat) + dgza(jscat)
     Else
        z2 = gz(jscat)
     End If
  End If
  t2 = ft(jscat)
  !trans
  rts2 = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
  !transend
  bex = (px1+px2)/(e1+e2)
  bey = (py1+py2)/(e1+e2)
  bez = (pz1+pz2)/(e1+e2)

  !lin-11/2015-ctest off
  !        write(99,*) 'iscat,jscat,etotalA=',iscat,jscat,e1+e2

  Call lorenz(e1, px1, py1, pz1, bex, bey, bez)
  !c      SAVE pxnew, ..., values for later use.
  px1 = pxnew
  py1 = pynew
  pz1 = pznew
  e1 = enenew

  pp2 = pxnew**2 + pynew**2 + pznew**2
  Call getht(iscat, jscat, pp2, that)
  theta = dacos(that/(2D0*pp2)+1D0)
  theta = dble(iff)*theta

  !       we boost to the cm frame, get rotation axis, and rotate 1 particle
  !       momentum

  Call lorenz(t1, x1, y1, z1, bex, bey, bez)

  x1 = pxnew
  y1 = pynew
  z1 = pznew

  Call lorenz(t2, x2, y2, z2, bex, bey, bez)

  x2 = pxnew
  y2 = pynew
  z2 = pznew

  !       notice now pxnew, ..., are new positions
  Call cropro(x1-x2, y1-y2, z1-z2, px1, py1, pz1)

  Call xnormv(xn1, xn2, xn3)

  !bz1/29/99
  !        call rotate(xn1, xn2, xn3, theta, px1, py1, pz1)
  Call zprota(xn1, xn2, xn3, theta, px1, py1, pz1)
  !bz1/29/99end

  !       we invert the momentum to get the other particle's momentum
  px2 = -px1
  py2 = -py1
  pz2 = -pz1
  !lin-4/13/01: modify in case m1, m2 are different:
  !        e2 = e1
  e2 = dsqrt(px2**2+py2**2+pz2**2+xmass(jscat)**2)

  !lin-11/2015-ctest off
  !        write(99,*) 'iscat,jscat,masses= ',iscat,jscat,
  !     1       xmass(iscat),xmass(jscat)

  !       boost the 2 particle 4 momentum back to lab frame
  Call lorenz(e1, px1, py1, pz1, -bex, -bey, -bez)
  px(iscat) = pxnew
  py(iscat) = pynew
  pz(iscat) = pznew
  e(iscat) = enenew
  Call lorenz(e2, px2, py2, pz2, -bex, -bey, -bez)
  px(jscat) = pxnew
  py(jscat) = pynew
  pz(jscat) = pznew
  e(jscat) = enenew

  !lin-11/2015-ctest off
  !        write(99,*) 'iscat,jscat,etotalB= ',iscat,jscat,
  !     1       e(iscat)+e(jscat)

  vx(iscat) = px(iscat)/e(iscat)
  vy(iscat) = py(iscat)/e(iscat)
  vz(iscat) = pz(iscat)/e(iscat)
  vx(jscat) = px(jscat)/e(jscat)
  vy(jscat) = py(jscat)/e(jscat)
  vz(jscat) = pz(jscat)/e(jscat)

  last(iscat) = jscat
  last(jscat) = iscat

  If (iconfg<=3) Then
     If (e(iscat)<=abs(pz(iscat))) Then
        rap(iscat) = 1000000.D0
     Else
        rap(iscat) = 0.5D0*log((e(iscat)+pz(iscat))/(e(iscat)-pz(iscat)))
     End If

     If (e(jscat)<=abs(pz(jscat))) Then
        rap(jscat) = 1000000.D0
     Else
        rap(jscat) = 0.5D0*log((e(jscat)+pz(jscat))/(e(jscat)-pz(jscat)))
     End If

     !trans
     rap1 = rap(iscat)
     rap2 = rap(jscat)

     If ((rap1<centy+0.5D0 .And. rap1>centy-0.5D0)) Then
        !              write (9, *) sqrt(ft(iscat) ** 2 - gz(iscat) ** 2), rts2
     End If
     If ((rap2<centy+0.5D0 .And. rap2>centy-0.5D0)) Then
        !              write (9, *) sqrt(ft(jscat) ** 2 - gz(jscat) ** 2), rts2
     End If
     !transend
  End If

  !lin-11/2015-ctest off
  !        write(99,*) 'iscat,jscat,xmp,xmu,that=',iscat,jscat,xmp,xmu,that

  Return
End Subroutine newmom

Subroutine getht(iscat, jscat, pp2, that)

  !       this subroutine is used to get \hat{t} for a particular processes

  Implicit Double Precision (A-H, O-Z)

  Parameter (hbarc=0.197327054D0)
  Parameter (maxptn=400001)
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  External ran1
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Save

  iseed = iseedp
  xmu2 = (hbarc*xmu)**2
  xmp2 = xmp**2
  xm2 = xmu2 + xmp2
  rx = ran1(iseed)
  that = xm2*(1D0+1D0/((1D0-xm2/(4D0*pp2+xm2))*rx-1D0))
  !test off isotropic scattering:
  !     &     + 1d0/((1d0 - xm2 / (4d0 * pp2 + xm2)) * ran1(2) - 1d0))
  !        if(izpc.eq.100) that=-4d0*pp2*ran1(2)
  If (izpc==100) that = -4D0*pp2*rx

  Return
End Subroutine getht

Subroutine ulist(t)
  !     this subroutine is used to update a new collision time list
  !       notice this t has been updated

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Save

  If (ictype==1 .Or. ictype==2 .Or. ictype==5 .Or. ictype==6) Then
     l = ifmpt
     Call ulist1(l, t)
  End If
  If (ictype/=1) Then
     l = iscat
     Call ulist1(l, t)
     If (jscat/=0) Then
        l = jscat
        Call ulist1(l, t)
     End If
  End If

  Return
End Subroutine ulist

Subroutine ulist1(l, t)
  !       this subroutine is used to update the interaction list when particle
  !       l is disturbed.

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  icels0 = icels(l)
  i1 = icels0/10000
  i2 = (icels0-i1*10000)/100
  i3 = icels0 - i1*10000 - i2*100
  !       save collision info for use when the collision is a collision with wall
  !       otherwise wallc will change icsta
  k = mod(icsta(l), 10)

  Call wallc(l, i1, i2, i3, t, tmin1)
  tmin = tmin1
  nc = 0

  If (i1==11 .And. i2==11 .And. i3==11) Then
     Call chkout(l, t, tmin, nc)
  Else
     If (iconfg==1) Then
        Call chkin1(l, i1, i2, i3, t, tmin, nc)
     Else If (iconfg==2) Then
        Call chkin2(l, i1, i2, i3, t, tmin, nc)
     Else If (iconfg==4) Then
        Call chkin3(l, i1, i2, i3, t, tmin, nc)
     Else If (iconfg==3 .Or. iconfg==5) Then
        Call chkcel(l, i1, i2, i3, t, tmin, nc)
     End If
  End If

  Call fixtim(l, t, tmin1, tmin, nc)

  Return
End Subroutine ulist1

Subroutine wallc(i, i1, i2, i3, t, tmin)
  !       this subroutine calculates the next time for collision with wall
  !       for particle i
  !       input particle label i,t
  !       output tmin collision time with wall, icsta(i) wall collision
  !       information

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  tmin = tlarge

  If (iconfg<=2 .Or. iconfg==4) Then
     !       if particle is inside the cube
     If ((i1>=1 .And. i1<=10) .Or. (i2>=1 .And. i2<=10) .Or. (i3>=1 .And. i3<=10)) Then
        Call wallc1(i, i1, i2, i3, t, tmin)
        !       if particle is outside the cube
     Else
        Call wallcb(i, t, tmin)
     End If
  Else If (iconfg==3 .Or. iconfg==5) Then
     Call wallc2(i, i1, i2, i3, t, tmin)
  End If

  Return
End Subroutine wallc

Subroutine wallc1(i, i1, i2, i3, t, tmin)
  !       this subroutine is used to get wall collision time
  !       when particle is inside the cube, it sets the icsta at the same time
  !       input i,i1,i2,i3,t
  !       output tmin, icsta(i)
  !       note the icsta is not finally set. we need further judgement in
  !       fixtim

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  x1p = gx(i)
  x2p = gy(i)
  x3p = gz(i)
  tf = ft(i)
  v1p = vx(i)
  v2p = vy(i)
  v3p = vz(i)

  If (t<size .And. tf<size) Then

     If (v1p>0D0) Then
        t1 = ((dble(i1)-5D0)*size1-x1p)/v1p + tf
     Else If (v1p<0D0) Then
        t1 = ((dble(i1)-6D0)*size1-x1p)/v1p + tf
     Else
        t1 = tlarge
     End If

     If (v2p>0D0) Then
        t2 = ((dble(i2)-5D0)*size2-x2p)/v2p + tf
     Else If (v2p<0D0) Then
        t2 = ((dble(i2)-6D0)*size2-x2p)/v2p + tf
     Else
        t2 = tlarge
     End If

     If (v3p>0D0) Then
        t3 = ((dble(i3)-5D0)*size3-x3p)/v3p + tf
     Else If (v3p<0D0) Then
        t3 = ((dble(i3)-6D0)*size3-x3p)/v3p + tf
     Else
        t3 = tlarge
     End If

     !       if a particle is on the wall, we don't collide it on the same wall

     !        if (t1 .eq. 0d0) t1 = tlarge
     !        if (t2 .eq. 0d0) t2 = tlarge
     !        if (t3 .eq. 0d0) t3 = tlarge

     tmin = min(t1, t2, t3)

     !       set icsta,
     !       after checking this is not an earlier collision comparing with
     !       a collision with another particle, we need to set icsta=0
     !       after checking whether there is also a particle collision
     !       at the same time, we need to reset the second bit of icsta

     If (tmin==t1) Then
        If (v1p>0D0) Then
           icsta(i) = 101
        Else
           icsta(i) = 102
        End If
     End If

     If (tmin==t2) Then
        If (v2p>0D0) Then
           icsta(i) = 103
        Else
           icsta(i) = 104
        End If
     End If

     If (tmin==t3) Then
        If (v3p>0D0) Then
           icsta(i) = 105
        Else
           icsta(i) = 106
        End If
     End If

     If (tmin<=size) Return

  End If

  If (v1p>(i1-5)*v1) Then
     t1 = ((i1-5)*(size1-v1*size)+v1p*tf-x1p)/(v1p-(i1-5)*v1)
  Else If (v1p<(i1-6)*v1) Then
     t1 = ((i1-6)*(size1-v1*size)+v1p*tf-x1p)/(v1p-(i1-6)*v1)
  Else
     t1 = tlarge
  End If

  If (v2p>(i2-5)*v2) Then
     t2 = ((i2-5)*(size2-v2*size)+v2p*tf-x2p)/(v2p-(i2-5)*v2)
  Else If (v2p<(i2-6)*v2) Then
     t2 = ((i2-6)*(size2-v2*size)+v2p*tf-x2p)/(v2p-(i2-6)*v2)
  Else
     t2 = tlarge
  End If

  If (v3p>(i3-5)*v3) Then
     t3 = ((i3-5)*(size3-v3*size)+v3p*tf-x3p)/(v3p-(i3-5)*v3)
  Else If (v3p<(i3-6)*v3) Then
     t3 = ((i3-6)*(size3-v3*size)+v3p*tf-x3p)/(v3p-(i3-6)*v3)
  Else
     t3 = tlarge
  End If

  !       if a particle is on the wall, we don't collide it on the same wall

  !        if (t1 .eq. 0d0) t1 = tlarge
  !        if (t2 .eq. 0d0) t2 = tlarge
  !        if (t3 .eq. 0d0) t3 = tlarge

  tmin = min(t1, t2, t3)

  !       set icsta,
  !       after checking this is not an earlier collision comparing with
  !       a collision with another particle, we need to set icsta=0
  !       after checking whether there is also a particle collision
  !       at the same time, we need to reset the second bit of icsta

  If (tmin==t1) Then
     If (v1p>(i1-5)*v1) Then
        icsta(i) = 101
     Else
        icsta(i) = 102
     End If
  End If

  If (tmin==t2) Then
     If (v2p>(i2-5)*v2) Then
        icsta(i) = 103
     Else
        icsta(i) = 104
     End If
  End If

  If (tmin==t3) Then
     If (v3p>(i3-5)*v3) Then
        icsta(i) = 105
     Else
        icsta(i) = 106
     End If
  End If

  Return
End Subroutine wallc1

Subroutine wallc2(i, i1, i2, i3, t, tmin)
  !       this subroutine is used to get wall collision time
  !       when particle is inside the cube, it sets the icsta at the same time
  !       input i,i1,i2,i3,t
  !       output tmin, icsta(i)
  !       note the icsta is not finally set. we need further judgement in
  !       fixtim

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  x1p = gx(i)
  x2p = gy(i)
  x3p = gz(i)
  tf = ft(i)
  v1p = vx(i)
  v2p = vy(i)
  v3p = vz(i)

  If (v1p>0D0) Then
     t1 = (5D0*size1-x1p)/v1p + tf
  Else If (v1p<0D0) Then
     t1 = (-5D0*size1-x1p)/v1p + tf
  Else
     t1 = tlarge
  End If

  If (v2p>0D0) Then
     t2 = (5D0*size2-x2p)/v2p + tf
  Else If (v2p<0D0) Then
     t2 = (-5D0*size2-x2p)/v2p + tf
  Else
     t2 = tlarge
  End If

  If (iconfg==5) Then
     If (v3p>0D0) Then
        t3 = (5D0*size3-x3p)/v3p + tf
     Else If (v3p<0D0) Then
        t3 = (-5D0*size3-x3p)/v3p + tf
     Else
        t3 = tlarge
     End If
  Else
     t3 = tlarge
  End If

  tmin = min(t1, t2, t3)

  !       set icsta,
  !       after checking this is not an earlier collision comparing with
  !       a collision with another particle, we need to set icsta=0
  !       after checking whether there is also a particle collision
  !       at the same time, we need to reset the second bit of icsta

  If (tmin==t1) Then
     If (v1p>0D0) Then
        icsta(i) = 101
     Else
        icsta(i) = 102
     End If
  End If

  If (tmin==t2) Then
     If (v2p>0D0) Then
        icsta(i) = 103
     Else
        icsta(i) = 104
     End If
  End If

  If (tmin==t3) Then
     If (v3p>0D0) Then
        icsta(i) = 105
     Else
        icsta(i) = 106
     End If
  End If

  Return
End Subroutine wallc2

Subroutine wallcb(i, t, tmin)
  !       this subroutine is used to calculate the wall collision time
  !       when the particle is outside the cube
  !       input i,t
  !       output tmin,icsta(i)
  !       note the icsta is not finally set. we need further judgement in
  !       fixtim

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  !       check if there is a collision by looking at the closest approach point
  !       and see if it's inside the cube

  If (size1==0D0 .Or. size2==0D0 .Or. size3==0D0) Return

  x1p = gx(i)
  x2p = gy(i)
  x3p = gz(i)
  v1p = vx(i)
  v2p = vy(i)
  v3p = vz(i)
  tf = ft(i)

  If (t<size .And. tf<size) Then
     If (x1p<-5D0*size1 .And. v1p>0D0) Then
        t1 = (-5D0*size1-x1p)/v1p + tf
     Else If (x1p>5D0*size1 .And. v1p<0D0) Then
        t1 = -(x1p-5D0*size1)/v1p + tf
     Else
        t1 = tlarge
     End If

     If (t1/=tlarge) Then
        x2pp = x2p + v2p*(t1-tf)
        x3pp = x3p + v3p*(t1-tf)
        If (x2pp<=-5D0*size2 .Or. x2pp>=5D0*size2 .Or. x3pp<=-5D0*size3 .Or. x3pp>=5D0*size3) t1 = tlarge
     End If

     If (x2p<-5D0*size2 .And. v2p>0D0) Then
        t2 = (-5D0*size2-x2p)/v2p + tf
     Else If (x2p>5D0*size2 .And. v2p<0D0) Then
        t2 = -(x2p-5D0*size2)/v2p + tf
     Else
        t2 = tlarge
     End If

     If (t2/=tlarge) Then
        x1pp = x1p + v1p*(t2-tf)
        x3pp = x3p + v3p*(t2-tf)
        If (x1pp<=-5D0*size1 .Or. x1pp>=5D0*size1 .Or. x3pp<=-5D0*size3 .Or. x3pp>=5D0*size3) t2 = tlarge
     End If

     If (x3p<-5D0*size3 .And. v3p>0D0) Then
        t3 = (-5D0*size3-x3p)/v3p + tf
     Else If (x3p>5D0*size3 .And. v3p<0D0) Then
        t3 = -(x3p-5D0*size3)/v3p + tf
     Else
        t3 = tlarge
     End If

     If (t3/=tlarge) Then
        x1pp = x1p + v1p*(t3-tf)
        x2pp = x2p + v2p*(t3-tf)
        If (x1pp<=-5D0*size1 .Or. x1pp>=5D0*size1 .Or. x2pp<=-5D0*size2 .Or. x2pp>=5D0*size2) t3 = tlarge
     End If

     tmin = min(t1, t2, t3)

     !       set icsta,
     !       after checking this is not an earlier collision comparing with
     !       a collision with another particle, we need to set icsta=0
     !       after checking whether there is also a particle collision
     !       at the same time, we need to reset the second bit of icsta

     If (tmin==t1) Then
        If (v1p>0D0) Then
           icsta(i) = 101
        Else
           icsta(i) = 102
        End If
     End If

     If (tmin==t2) Then
        If (v2p>0D0) Then
           icsta(i) = 103
        Else
           icsta(i) = 104
        End If
     End If

     If (tmin==t3) Then
        If (v3p>0D0) Then
           icsta(i) = 105
        Else
           icsta(i) = 106
        End If
     End If

     If (tmin<=size) Return

  End If

  !       notice now x1q, x2q, x3q are coordinates at time t
  x1q = x1p + v1p*(t-tf)
  x2q = x2p + v2p*(t-tf)
  x3q = x3p + v3p*(t-tf)

  If (x1q<-5D0*(size1+v1*(t-size)) .And. v1p>-5D0*v1) Then
     t1 = (-5D0*(size1-v1*size)+v1p*tf-x1p)/(v1p-(-5D0)*v1)
     icsta1 = 101
  Else If (x1q>5D0*(size1+v1*(t-size)) .And. v1p<5D0*v1) Then
     t1 = (5D0*(size1-v1*size)+v1p*tf-x1p)/(v1p-5D0*v1)
     icsta1 = 102
  Else
     t1 = tlarge
  End If

  If (t1/=tlarge) Then
     x2pp = x2p + v2p*(t1-tf)
     x3pp = x3p + v3p*(t1-tf)
     If (x2pp<=-5D0*(size2+v2*(t1-size)) .Or. x2pp>=5D0*(size2+v2*(t1-size)) .Or. x3pp<=-5D0*(size3+v3*(t1-size)) .Or. x3pp>=5D0*(size3+v3*(t1-size))) t1 = tlarge
  End If

  If (x2q<-5D0*(size2+v2*(t-size)) .And. v2p>-5D0*v2) Then
     t2 = (-5D0*(size2-v2*size)+v2p*tf-x2p)/(v2p-(-5D0)*v2)
     icsta2 = 103
  Else If (x2q>5D0*(size2+v2*(t-size)) .And. v2p<5D0*v2) Then
     t2 = (5D0*(size2-v2*size)+v2p*tf-x2p)/(v2p-5D0*v2)
     icsta2 = 104
  Else
     t2 = tlarge
  End If

  If (t2/=tlarge) Then
     x1pp = x1p + v1p*(t2-tf)
     x3pp = x3p + v3p*(t2-tf)
     If (x1pp<=-5D0*(size1+v1*(t2-size)) .Or. x1pp>=5D0*(size1+v1*(t2-size)) .Or. x3pp<=-5D0*(size3+v3*(t2-size)) .Or. x3pp>=5D0*(size3+v3*(t2-size))) t2 = tlarge
  End If

  If (x3q<-5D0*(size3+v3*(t-size)) .And. v3p>-5D0*v3) Then
     t3 = (-5D0*(size3-v3*size)+v3p*tf-x3p)/(v3p-(-5D0)*v3)
     icsta3 = 105
  Else If (x3q>5D0*(size3+v3*(t-size)) .And. v3p<5D0*v3) Then
     t3 = (5D0*(size3-v3*size)+v3p*tf-x3p)/(v3p-5D0*v3)
     icsta3 = 106
  Else
     t3 = tlarge
  End If

  If (t3/=tlarge) Then
     x2pp = x2p + v2p*(t3-tf)
     x1pp = x1p + v1p*(t3-tf)
     If (x2pp<=-5D0*(size2+v2*(t3-size)) .Or. x2pp>=5D0*(size2+v2*(t3-size)) .Or. x1pp<=-5D0*(size1+v1*(t3-size)) .Or. x1pp>=5D0*(size1+v1*(t3-size))) t3 = tlarge
  End If

  tmin = min(t1, t2, t3)

  !       set icsta,
  !       after checking this is not an earlier collision comparing with
  !       a collision with another particle, we need to set icsta=0
  !       after checking whether there is also a particle collision
  !       at the same time, we need to reset the second bit of icsta

  If (tmin==t1) Then
     icsta(i) = icsta1
  Else If (tmin==t2) Then
     icsta(i) = icsta2
  Else If (tmin==t3) Then
     icsta(i) = icsta3
  End If

  Return
End Subroutine wallcb

Subroutine chkout(l, t, tmin, nc)
  !       this subroutine is used to check the collisions with particles in
  !       surface cells to see if we can get a smaller collision time than tmin
  !       with particle nc, when the colliding particle is outside the cube
  !       input l,t,tmin,nc
  !       output tmin, nc

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Save

  m1 = 11
  m2 = 11
  m3 = 11
  Call chkcel(l, m1, m2, m3, t, tmin, nc)

  Do i = 1, 10
     Do j = 1, 10
        Do k = 1, 10
           If (i==1 .Or. i==10 .Or. j==1 .Or. j==10 .Or. k==1 .Or. k==10) Call chkcel(l, i, j, k, t, tmin, nc)
        End Do
     End Do
  End Do

  Return
End Subroutine chkout

Subroutine chkin1(l, i1, i2, i3, t, tmin, nc)
  !       this subroutine is used to check collisions for particle inside
  !       the cube
  !       and update the afftected particles through chkcel

  Implicit Double Precision (A-H, O-Z)
  Save

  !       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0

  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Then
              Call chkcel(l, i, j, k, t, tmin, nc)
           Else If (itest==0) Then
              m1 = 11
              m2 = 11
              m3 = 11
              Call chkcel(l, m1, m2, m3, t, tmin, nc)
              itest = 1
           End If
        End Do
     End Do
  End Do

  Return
End Subroutine chkin1

Subroutine chkin2(l, i1, i2, i3, t, tmin, nc)
  !       this subroutine is used to check collisions for particle inside
  !       the cube
  !       and update the afftected particles through chkcel

  Implicit Double Precision (A-H, O-Z)
  Save

  !       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0

  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ib = j
           ic = k
           If (k>=1 .And. k<=10) Then
              If (i==0) ia = 10
              If (i==11) ia = 1
              If (j==0) ib = 10
              If (j==11) ib = 1
              Call chkcel(l, ia, ib, ic, t, tmin, nc)
           End If
        End Do
     End Do
  End Do

  Return
End Subroutine chkin2

Subroutine chkin3(l, i1, i2, i3, t, tmin, nc)
  !       this subroutine is used to check collisions for particle inside
  !       the cube
  !       and update the afftected particles through chkcel

  Implicit Double Precision (A-H, O-Z)
  Save

  !       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0

  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (i==0) Then
              ia = 10
           Else If (i==11) Then
              ia = 1
           Else
              ia = i
           End If
           If (j==0) Then
              ib = 10
           Else If (j==11) Then
              ib = 1
           Else
              ib = j
           End If
           If (k==0) Then
              ic = 10
           Else If (k==11) Then
              ic = 1
           Else
              ic = k
           End If
           Call chkcel(l, ia, ib, ic, t, tmin, nc)
        End Do
     End Do
  End Do

  Return
End Subroutine chkin3

Subroutine chkcel(il, i1, i2, i3, t, tmin, nc)
  !       this program is used to check through all the particles
  !       in the cell (i1,i2,i3) and see if we can get a particle collision
  !       with time less than the original input tmin ( the collision time of
  !       il with the wall
  !       and update the affected particles

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist2/icell, icel(10, 10, 10)
  !c      SAVE /ilist2/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Save

  If (iconfg==3 .Or. iconfg==5) Then
     jj = ichkpt
     Do j = 1, jj
        Call ck(j, ick)
        !     10/24/02 get rid of argument usage mismatch in ud2():
        jud2 = j
        !              if (ick .eq. 1) call ud2(j, il, t, tmin, nc)
        If (ick==1) Call ud2(jud2, il, t, tmin, nc)
     End Do
     Return
  End If

  If (i1==11 .And. i2==11 .And. i3==11) Then
     l = icell
  Else
     l = icel(i1, i2, i3)
  End If

  !       if there is no particle
  If (l==0) Then
     Return
  End If
  j = nic(l)
  !       if there is only one particle
  If (j==0) Then
     Call ck(l, ick)
     If (ick==1) Call ud2(l, il, t, tmin, nc)

     !       if there are many particles
  Else

     !       we don't worry about the other colliding particle because it's
     !       set in last(), and will be checked in ud2

     Call ck(l, ick)
     If (ick==1) Call ud2(l, il, t, tmin, nc)

     Do While (j/=l)
        Call ck(j, ick)
        If (ick==1) Call ud2(j, il, t, tmin, nc)
        j = nic(j)
     End Do
  End If

  Return
End Subroutine chkcel

Subroutine ck(l, ick)
  !       this subroutine is used for chcell to check whether l should be
  !       checked or not for updating tmin, nc
  !       input l
  !       output ick
  !       if ick=1, l should be checked, otherwise it should not be.

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Save

  ick = 1
  If (ictype==1) Then
     If (l==ifmpt) ick = 0
  Else If (ictype==0 .Or. ictype==3 .Or. ictype==4) Then
     If (l==iscat .Or. l==jscat) ick = 0
  Else
     If (l==iscat .Or. l==jscat .Or. l==ifmpt) ick = 0
  End If
  !       notice il is either iscat or jscat, or ifmpt, we deal with them
  !       seperately according to ictype

  Return
End Subroutine ck

Subroutine dchout(l, ii, t)
  !       this subroutine is used to check collisions of l with particles when
  !       l is outside the cube and the collision just happened is a collision
  !       including a collision with wall (hence we need to use dcheck to throw
  !       away old collisions that are not in the new neighboring cells.

  !       input l,t

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  External integ

  tt = ft(l)
  td = t - size
  x1 = gx(l) + vx(l)*(t-tt)
  x2 = gy(l) + vy(l)*(t-tt)
  x3 = gz(l) + vz(l)*(t-tt)
  If (td<=0D0) Then
     i1 = integ(x1/size1) + 6
     i2 = integ(x2/size2) + 6
     i3 = integ(x3/size3) + 6
     If (integ(x1/size1)==x1/size1 .And. vx(l)<0D0) i1 = i1 - 1
     If (integ(x2/size2)==x2/size2 .And. vy(l)<0D0) i2 = i2 - 1
     If (integ(x3/size3)==x3/size3 .And. vz(l)<0D0) i3 = i3 - 1
  Else
     i1 = integ(x1/(size1+v1*td)) + 6
     i2 = integ(x2/(size2+v2*td)) + 6
     i3 = integ(x3/(size3+v3*td)) + 6
     !     10/24/02 (i) below should be (l):
     If (integ(x1/(size1+v1*td))==x1/(size1+v1*td) .And. vx(l)<(i1-6)*v1) i1 = i1 - 1
     !     &        vx(i) .lt. (i1 - 6) * v1) i1 = i1 - 1
     If (integ(x2/(size2+v2*td))==x2/(size2+v2*td) .And. vy(l)<(i2-6)*v2) i2 = i2 - 1
     !     &        vy(i) .lt. (i2 - 6) * v2) i2 = i2 - 1
     If (integ(x3/(size3+v3*td))==x3/(size3+v3*td) .And. vz(l)<(i3-6)*v3) i3 = i3 - 1
     !     &        vz(i) .lt. (i3 - 6) * v3) i3 = i3 - 1
  End If

  If (ii==1) Then
     i = 9
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If

  If (ii==2) Then
     i = 2
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If

  If (ii==3) Then
     j = 9
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If

  If (ii==4) Then
     j = 2
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If

  If (ii==5) Then
     k = 9
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If

  If (ii==6) Then
     k = 2
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Then
              Call dchcel(l, i, j, k, t)
           End If
        End Do
     End Do
  End If

  Return
End Subroutine dchout

Subroutine dchin1(l, ii, i1, i2, i3, t)
  !       this subroutine is used to check collisions for particle inside
  !       the cube when the collision just happened is a collision including
  !       collision with wall
  !       and update the afftected particles through chkcel

  !       input l,ii(specifying the direction of the wall collision),
  !          i1,i2,i3, (specifying the position of the cell
  !                    we are going to check)
  !          t

  Implicit Double Precision (A-H, O-Z)
  Save

  !       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0

  If (ii==1) Then
     If (i1==1) Goto 100
     If (i1==2) Then
        If (i2>=2 .And. i2<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     i = i1 - 2
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If

  If (ii==2) Then
     If (i1==10) Goto 100
     If (i1==9) Then
        If (i2>=2 .And. i2<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     i = i1 + 2
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If

  If (ii==3) Then
     If (i2==1) Goto 100
     If (i2==2) Then
        If (i1>=2 .And. i1<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     j = i2 - 2
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If

  If (ii==4) Then
     If (i2==10) Goto 100
     If (i2==9) Then
        If (i1>=2 .And. i1<=9 .And. i3>=2 .And. i3<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     j = i2 + 2
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. k>=1 .And. k<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If

  If (ii==5) Then
     If (i3==1) Goto 100
     If (i3==2) Then
        If (i1>=2 .And. i1<=9 .And. i2>=2 .And. i2<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     k = i3 - 2
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If

  If (ii==6) Then
     If (i3==10) Goto 100
     If (i3==9) Then
        If (i1>=2 .And. i1<=9 .And. i2>=2 .And. i2<=9) Then
           i = 11
           j = 11
           k = 11
           Call dchcel(l, i, j, k, t)
        End If
        Goto 100
     End If
     k = i3 + 2
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10) Call dchcel(l, i, j, k, t)
        End Do
     End Do
  End If

100 Continue

  Return
End Subroutine dchin1

Subroutine dchin2(l, ii, i1, i2, i3, t)
  !       this subroutine is used to check collisions for particle inside
  !       the cube when the collision just happened is a collision including
  !       collision with wall
  !       and update the afftected particles through chkcel

  !       input l,ii(specifying the direction of the wall collision),
  !          i1,i2,i3, (specifying the position of the cell
  !                    we are going to check)
  !          t

  Implicit Double Precision (A-H, O-Z)
  Save

  If (ii==1) Then
     i = i1 - 2
     If (i<=0) i = i + 10
     ia = i
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ib = j
           ic = k
           If (j==0) ib = 10
           If (j==11) ib = 1
           If (k>=1 .And. k<=10) Then
              Call dchcel(l, ia, ib, ic, t)
           End If
        End Do
     End Do
  End If

  If (ii==2) Then
     i = i1 + 2
     If (i>=11) i = i - 10
     ia = i
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ib = j
           ic = k
           If (j==0) ib = 10
           If (j==11) ib = 1
           If (k>=1 .And. k<=10) Then
              Call dchcel(l, ia, ib, ic, t)
           End If
        End Do
     End Do
  End If

  If (ii==3) Then
     j = i2 - 2
     If (j<=0) j = j + 10
     ib = j
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ic = k
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (k>=1 .And. k<=10) Then
              Call dchcel(l, ia, ib, ic, t)
           End If
        End Do
     End Do
  End If

  If (ii==4) Then
     j = i2 + 2
     If (j>=11) j = j - 10
     ib = j
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ic = k
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (k>=1 .And. k<=10) Then
              Call dchcel(l, ia, ib, ic, t)
           End If
        End Do
     End Do
  End If

  If (ii==5) Then
     If (i3==2) Goto 100
     k = i3 - 2
     ic = k
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           ia = i
           ib = j
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (j==0) ib = 10
           If (j==11) ib = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If

  If (ii==6) Then
     If (i3==9) Goto 100
     k = i3 + 2
     ic = k
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           ia = i
           ib = j
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (j==0) ib = 10
           If (j==11) ib = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If

100 Continue

  Return
End Subroutine dchin2

Subroutine dchin3(l, ii, i1, i2, i3, t)
  !       this subroutine is used to check collisions for particle inside
  !       the cube when the collision just happened is a collision including
  !       collision with wall
  !       and update the afftected particles through chkcel

  !       input l,ii(specifying the direction of the wall collision),
  !          i1,i2,i3, (specifying the position of the cell
  !                    we are going to check)
  !          t

  Implicit Double Precision (A-H, O-Z)
  Save

  If (ii==1) Then
     i = i1 - 2
     If (i<=0) i = i + 10
     ia = i
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ib = j
           ic = k
           If (j==0) ib = 10
           If (j==11) ib = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If

  If (ii==2) Then
     i = i1 + 2
     If (i>=11) i = i - 10
     ia = i
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ib = j
           ic = k
           If (j==0) ib = 10
           If (j==11) ib = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If

  If (ii==3) Then
     j = i2 - 2
     If (j<=0) j = j + 10
     ib = j
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ic = k
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If

  If (ii==4) Then
     j = i2 + 2
     If (j>=11) j = j - 10
     ib = j
     Do i = i1 - 1, i1 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ic = k
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (k==0) ic = 10
           If (k==11) ic = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If

  If (ii==5) Then
     k = i3 - 2
     If (k<=0) k = k + 10
     ic = k
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           ia = i
           ib = j
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (j==0) ib = 10
           If (j==11) ib = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If

  If (ii==6) Then
     k = i3 + 2
     If (k>=11) k = k - 10
     ic = k
     Do i = i1 - 1, i1 + 1
        Do j = i2 - 1, i2 + 1
           ia = i
           ib = j
           If (i==0) ia = 10
           If (i==11) ia = 1
           If (j==0) ib = 10
           If (j==11) ib = 1
           Call dchcel(l, ia, ib, ic, t)
        End Do
     End Do
  End If
  !
  Return
End Subroutine dchin3

Subroutine dchcel(l, i, j, k, t)
  !       this subroutine is used to recalculate next collision time for
  !       particles in the cell i,j,k if the next collision partener is l

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist2/icell, icel(10, 10, 10)
  !c      SAVE /ilist2/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  If (i==11 .Or. j==11 .Or. k==11) Then
     If (.Not. (i==11 .And. j==11 .And. k==11)) Stop 'cerr'
     m = icell
  Else
     m = icel(i, j, k)
  End If

  If (m==0) Return
  If (next(m)==l) Then
     tm = tlarge
     last0 = 0
     Call reor(t, tm, m, last0)
  End If
  n = nic(m)
  If (n==0) Return
  Do While (n/=m)
     If (next(n)==l) Then
        tm = tlarge
        last0 = 0
        Call reor(t, tm, n, last0)
     End If
     n = nic(n)
  End Do

  Return
End Subroutine dchcel

Subroutine fixtim(l, t, tmin1, tmin, nc)
  !       this subroutine is used to compare the collision time with wall tmin1
  !       and new collision time with particles for particle l
  !       when used in ulist, input nc may be 0, which indicates no particle
  !       collisions happen before wall collision, of course, then tmin=tmin1

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  k = nc
  If (tmin<tmin1) Then
     ot(l) = tmin
     If (ct(l)<tmin1) Then
        icsta(l) = 0
     Else
        icsta(l) = icsta(l) + 10
     End If
     next(l) = k
  Else If (tmin==tmin1) Then
     ot(l) = tmin
     If (nc==0) Then
        next(l) = 0
     Else
        icsta(l) = icsta(l) + 10
        next(l) = k
     End If
  Else
     ot(l) = tmin1
     next(l) = 0
  End If

  Return
End Subroutine fixtim

Subroutine ud2(i, j, t, tmin, nc)
  !       this subroutine is used to update next(i), ct(i), ot(i),
  !        and get tmin, nc for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  Call isco(i, j, allok, tm, t1, t2)

  If (allok) Then
     !       tm eq tmin, change nc to make sure fixtime get the collision with both
     !       wall and particle

     If (tm<tmin) Then
        tmin = tm
        ct(j) = t2
        nc = i
        If (iconfg==3 .Or. iconfg==5) Then
           dgxa(j) = jxa*10D0*size1
           dgya(j) = jya*10D0*size2
           If (iconfg==5) Then
              dgza(j) = jza*10D0*size3
           End If
        End If
     End If

     If (tm<=ot(i)) Then
        ct(i) = t1
        icels0 = icels(i)
        i1 = icels0/10000
        i2 = (icels0-i1*10000)/100
        i3 = icels0 - i1*10000 - i2*100
        Call wallc(i, i1, i2, i3, t, tmin1)
        Call fixtim(i, t, tmin1, tm, j)
        If (iconfg==3 .Or. iconfg==5) Then
           dgxa(i) = -jxa*10D0*size1
           dgya(i) = -jya*10D0*size2
           If (iconfg==5) Then
              dgza(i) = -jza*10D0*size3
           End If
        End If
     End If

     If (tm>ot(i) .And. next(i)==j) Then
        ct(i) = t1
        Call reor(t, tm, i, j)
     End If

  Else If (next(i)==j) Then

     tm = tlarge

     Call reor(t, tm, i, j)

  End If

  Return
End Subroutine ud2





!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn


Subroutine isco(i, j, allok, tm, t1, t2)

  Implicit Double Precision (A-H, O-Z)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Save

  Logical allok

  iorder = iordsc/10
  If (iconfg==1) Then
     If (iorder==1) Then
        Call isco1(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco2(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco3(i, j, allok, tm, t1, t2)
     End If
  Else If (iconfg==2 .Or. iconfg==4) Then
     If (iorder==1) Then
        Call isco4(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco5(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco6(i, j, allok, tm, t1, t2)
     End If
  Else If (iconfg==3) Then
     If (iorder==1) Then
        Call isco7(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco8(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco9(i, j, allok, tm, t1, t2)
     End If
  Else If (iconfg==5) Then
     If (iorder==1) Then
        Call isco10(i, j, allok, tm, t1, t2)
     Else If (iorder==2) Then
        Call isco11(i, j, allok, tm, t1, t2)
     Else If (iorder==3) Then
        Call isco12(i, j, allok, tm, t1, t2)
     End If
  End If

  Return
End Subroutine isco

Subroutine isco1(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations
  i1 = i
  i2 = j

  p4 = ft(i2) - ft(i1)
  p1 = gx(i2) - gx(i1)
  p2 = gy(i2) - gy(i1)
  p3 = gz(i2) - gz(i1)

  q4 = e(i1)
  q1 = px(i1)
  q2 = py(i1)
  q3 = pz(i1)

  r4 = e(i2)
  r1 = px(i2)
  r2 = py(i2)
  r3 = pz(i2)

  a = p4*q4 - p1*q1 - p2*q2 - p3*q3
  b = p4*r4 - p1*r1 - p2*r2 - p3*r3
  c = q4*q4 - q1*q1 - q2*q2 - q3*q3
  d = r4*r4 - r1*r1 - r2*r2 - r3*r3
  ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
  f = p4*p4 - p1*p1 - p2*p2 - p3*p3

  !       make sure particle 2 formed early
  h = a + b
  If (h>0D0) Then
     g = a
     a = -b
     b = -g

     g = c
     c = d
     d = g

     i1 = j
     i2 = i
  End If

  !       check the approaching criteria
  If (allok) Then

     vp = a*d - b*ee

     allok = allok .And. vp < 0D0

  End If

  !       check the closest approach distance criteria
  If (allok) Then

     dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)

     allok = allok .And. dm2 < cutof2

  End If

  !       check the time criteria
  If (allok) Then

     tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
     tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
     tm = 0.5D0*(tc1+tc2)

     allok = allok .And. tm > ft(i) .And. tm > ft(j)

  End If

  !        check rts cut
  If (allok) Then

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (h>0D0) Then
     t1 = tm
     t2 = tm
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco1

Subroutine isco2(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations
  i1 = i
  i2 = j

  p4 = ft(i2) - ft(i1)
  p1 = gx(i2) - gx(i1)
  p2 = gy(i2) - gy(i1)
  p3 = gz(i2) - gz(i1)

  q4 = e(i1)
  q1 = px(i1)
  q2 = py(i1)
  q3 = pz(i1)

  r4 = e(i2)
  r1 = px(i2)
  r2 = py(i2)
  r3 = pz(i2)

  a = p4*q4 - p1*q1 - p2*q2 - p3*q3
  b = p4*r4 - p1*r1 - p2*r2 - p3*r3
  c = q4*q4 - q1*q1 - q2*q2 - q3*q3
  d = r4*r4 - r1*r1 - r2*r2 - r3*r3
  ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
  f = p4*p4 - p1*p1 - p2*p2 - p3*p3

  !       make sure particle 2 formed early
  h = a + b
  If (h>0D0) Then
     g = a
     a = -b
     b = -g

     g = c
     c = d
     d = g

     i1 = j
     i2 = i
  End If

  !       check the approaching criteria
  If (allok) Then

     vp = a*d - b*ee

     allok = allok .And. vp < 0D0

  End If

  !       check the closest approach distance criteria
  If (allok) Then

     dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)

     allok = allok .And. dm2 < cutof2

  End If

  !       check the time criteria
  If (allok) Then

     tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
     tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
     If (iordsc==20) Then
        tm = min(tc1, tc2)
     Else If (iordsc==21) Then
        tm = 0.5D0*(tc1+tc2)
     Else
        tm = max(tc1, tc2)
     End If

     allok = allok .And. tm > ft(i) .And. tm > ft(j)

  End If

  !        check rts cut
  If (allok) Then

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (h>0D0) Then
     t1 = tc2
     t2 = tc1
  Else
     t1 = tc1
     t2 = tc2
  End If

  Return
End Subroutine isco2

Subroutine isco3(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  If (ft(i)>=ft(j)) Then
     i1 = j
     i2 = i
  Else
     i1 = i
     i2 = j
  End If

  If (allok) Then

     t1 = ft(i1)
     vx1 = vx(i1)
     vy1 = vy(i1)
     vz1 = vz(i1)

     t2 = ft(i2)

     dvx = vx(i2) - vx1
     dvy = vy(i2) - vy1
     dvz = vz(i2) - vz1

     dt = t2 - t1

     dx = gx(i2) - gx(i1) - vx1*dt
     dy = gy(i2) - gy(i1) - vy1*dt
     dz = gz(i2) - gz(i1) - vz1*dt

     vp = dvx*dx + dvy*dy + dvz*dz

     allok = allok .And. vp < 0D0

  End If

  If (allok) Then

     v2 = dvx*dvx + dvy*dvy + dvz*dvz

     If (v2==0D0) Then
        tm = tlarge
     Else
        tm = t2 - vp/v2
     End If

     !       note now tm is the absolute time

     allok = allok .And. tm > t1 .And. tm > t2

  End If

  If (allok) Then

     dgx = dx - dvx*t2
     dgy = dy - dvy*t2
     dgz = dz - dvz*t2

     dm2 = -v2*tm**2 + dgx*dgx + dgy*dgy + dgz*dgz

     allok = allok .And. dm2 < cutof2

  End If

  If (allok) Then

     e1 = e(i1)
     px1 = px(i1)
     py1 = py(i1)
     pz1 = pz(i1)
     e2 = e(i2)
     px2 = px(i2)
     py2 = py(i2)
     pz2 = pz(i2)

     rts2 = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco3

Subroutine isco4(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations

  icels1 = icels(i)
  ii1 = icels1/10000
  jj1 = (icels1-ii1*10000)/100
  kk1 = icels1 - ii1*10000 - jj1*100
  icels2 = icels(j)
  ii2 = icels2/10000
  jj2 = (icels2-ii2*10000)/100
  kk2 = icels2 - ii2*10000 - jj2*100

  i1 = i
  i2 = j

  p4 = ft(i2) - ft(i1)
  p1 = gx(i2) - gx(i1)
  p2 = gy(i2) - gy(i1)
  p3 = gz(i2) - gz(i1)

  If (ii1-ii2>5) Then
     p1 = p1 + 10D0*size1
  Else If (ii1-ii2<-5) Then
     p1 = p1 - 10D0*size1
  End If
  If (jj1-jj2>5) Then
     p2 = p2 + 10D0*size2
  Else If (jj1-jj2<-5) Then
     p2 = p2 - 10D0*size2
  End If
  If (kk1-kk2>5) Then
     p3 = p3 + 10D0*size3
  Else If (kk1-kk2<-5) Then
     p3 = p3 - 10D0*size3
  End If

  q4 = e(i1)
  q1 = px(i1)
  q2 = py(i1)
  q3 = pz(i1)

  r4 = e(i2)
  r1 = px(i2)
  r2 = py(i2)
  r3 = pz(i2)

  a = p4*q4 - p1*q1 - p2*q2 - p3*q3
  b = p4*r4 - p1*r1 - p2*r2 - p3*r3
  c = q4*q4 - q1*q1 - q2*q2 - q3*q3
  d = r4*r4 - r1*r1 - r2*r2 - r3*r3
  ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
  f = p4*p4 - p1*p1 - p2*p2 - p3*p3

  !       make sure particle 2 formed early
  h = a + b
  If (h>0D0) Then
     g = a
     a = -b
     b = -g

     g = c
     c = d
     d = g

     i1 = j
     i2 = i
  End If

  !       check the approaching criteria
  If (allok) Then

     vp = a*d - b*ee

     allok = allok .And. vp < 0D0

  End If

  !       check the closest approach distance criteria
  If (allok) Then

     dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)

     allok = allok .And. dm2 < cutof2

  End If

  !       check the time criteria
  If (allok) Then

     tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
     tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
     tm = 0.5D0*(tc1+tc2)

     allok = allok .And. tm > ft(i) .And. tm > ft(j)

  End If

  !        check rts cut
  If (allok) Then

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (h>0D0) Then
     t1 = tm
     t2 = tm
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco4

Subroutine isco5(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations

  icels1 = icels(i)
  ii1 = icels1/10000
  jj1 = (icels1-ii1*10000)/100
  kk1 = icels1 - ii1*10000 - jj1*100
  icels2 = icels(j)
  ii2 = icels2/10000
  jj2 = (icels2-ii2*10000)/100
  kk2 = icels2 - ii2*10000 - jj2*100

  i1 = i
  i2 = j

  p4 = ft(i2) - ft(i1)
  p1 = gx(i2) - gx(i1)
  p2 = gy(i2) - gy(i1)
  p3 = gz(i2) - gz(i1)

  If (ii1-ii2>5) Then
     p1 = p1 + 10D0*size1
  Else If (ii1-ii2<-5) Then
     p1 = p1 - 10D0*size1
  End If
  If (jj1-jj2>5) Then
     p2 = p2 + 10D0*size2
  Else If (jj1-jj2<-5) Then
     p2 = p2 - 10D0*size2
  End If
  If (kk1-kk2>5) Then
     p3 = p3 + 10D0*size3
  Else If (kk1-kk2<-5) Then
     p3 = p3 - 10D0*size3
  End If

  q4 = e(i1)
  q1 = px(i1)
  q2 = py(i1)
  q3 = pz(i1)

  r4 = e(i2)
  r1 = px(i2)
  r2 = py(i2)
  r3 = pz(i2)

  a = p4*q4 - p1*q1 - p2*q2 - p3*q3
  b = p4*r4 - p1*r1 - p2*r2 - p3*r3
  c = q4*q4 - q1*q1 - q2*q2 - q3*q3
  d = r4*r4 - r1*r1 - r2*r2 - r3*r3
  ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
  f = p4*p4 - p1*p1 - p2*p2 - p3*p3

  !       make sure particle 2 formed early
  h = a + b
  If (h>0D0) Then
     g = a
     a = -b
     b = -g

     g = c
     c = d
     d = g

     i1 = j
     i2 = i
  End If

  !       check the approaching criteria
  If (allok) Then

     vp = a*d - b*ee

     allok = allok .And. vp < 0D0

  End If

  !       check the closest approach distance criteria
  If (allok) Then

     dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)

     allok = allok .And. dm2 < cutof2

  End If

  !       check the time criteria
  If (allok) Then

     tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
     tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
     If (iordsc==20) Then
        tm = min(tc1, tc2)
     Else If (iordsc==21) Then
        tm = 0.5D0*(tc1+tc2)
     Else
        tm = max(tc1, tc2)
     End If

     allok = allok .And. tm > ft(i) .And. tm > ft(j)

  End If

  !        check rts cut
  If (allok) Then

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (h>0D0) Then
     t1 = tc2
     t2 = tc1
  Else
     t1 = tc1
     t2 = tc2
  End If

  Return
End Subroutine isco5

Subroutine isco6(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  If (ft(i)>=ft(j)) Then
     i1 = j
     i2 = i
  Else
     i1 = i
     i2 = j
  End If

  icels1 = icels(i1)
  ii1 = icels1/10000
  jj1 = (icels1-ii1*10000)/100
  kk1 = icels1 - ii1*10000 - jj1*100
  icels2 = icels(i2)
  ii2 = icels2/10000
  jj2 = (icels2-ii2*10000)/100
  kk2 = icels2 - ii2*10000 - jj2*100

  If (allok) Then

     t1 = ft(i1)
     vx1 = vx(i1)
     vy1 = vy(i1)
     vz1 = vz(i1)

     t2 = ft(i2)

     dvx = vx(i2) - vx1
     dvy = vy(i2) - vy1
     dvz = vz(i2) - vz1

     dt = t2 - t1

     dx = gx(i2) - gx(i1) - vx1*dt
     dy = gy(i2) - gy(i1) - vy1*dt
     dz = gz(i2) - gz(i1) - vz1*dt

     If (ii1-ii2>5) Then
        dx = dx + 10D0*size1
     Else If (ii1-ii2<-5) Then
        dx = dx - 10D0*size1
     End If

     If (jj1-jj2>5) Then
        dy = dy + 10D0*size2
     Else If (jj1-jj2<-5) Then
        dy = dy - 10D0*size2
     End If

     If (kk1-kk2>5) Then
        dz = dz + 10D0*size3
     Else If (kk1-kk2<-5) Then
        dz = dz - 10D0*size3
     End If

     vp = dvx*dx + dvy*dy + dvz*dz

     allok = allok .And. vp < 0D0

  End If

  If (allok) Then

     v2p = dvx*dvx + dvy*dvy + dvz*dvz

     If (v2p==0D0) Then
        tm = tlarge
     Else
        tm = t2 - vp/v2p
     End If

     !       note now tm is the absolute time

     allok = allok .And. tm > t1 .And. tm > t2

  End If

  If (allok) Then

     dgx = dx - dvx*t2
     dgy = dy - dvy*t2
     dgz = dz - dvz*t2

     dm2 = -v2p*tm**2 + dgx*dgx + dgy*dgy + dgz*dgz

     allok = allok .And. dm2 < cutof2

  End If

  If (allok) Then

     e1 = e(i1)
     px1 = px(i1)
     py1 = py(i1)
     pz1 = pz(i1)
     e2 = e(i2)
     px2 = px(i2)
     py2 = py(i2)
     pz2 = pz(i2)

     rts2 = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco6

Subroutine isco7(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok, allokp

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations

  tm = tlarge

  If (allok) Then
     Do ii = -1, 1
        Do jj = -1, 1

           allokp = .True.

           i1 = i
           i2 = j

           p4 = ft(j) - ft(i)
           p1 = gx(j) - gx(i)
           p2 = gy(j) - gy(i)
           p3 = gz(j) - gz(i)

           p1 = p1 + ii*10D0*size1
           p2 = p2 + jj*10D0*size2

           q4 = e(i)
           q1 = px(i)
           q2 = py(i)
           q3 = pz(i)

           r4 = e(j)
           r1 = px(j)
           r2 = py(j)
           r3 = pz(j)

           a = p4*q4 - p1*q1 - p2*q2 - p3*q3
           b = p4*r4 - p1*r1 - p2*r2 - p3*r3
           c = q4*q4 - q1*q1 - q2*q2 - q3*q3
           d = r4*r4 - r1*r1 - r2*r2 - r3*r3
           ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
           f = p4*p4 - p1*p1 - p2*p2 - p3*p3

           !       make sure particle 2 formed early
           h = a + b
           If (h>0D0) Then
              g = a
              a = -b
              b = -g
              g = c
              c = d
              d = g
              i1 = j
              i2 = i
           End If

           !       check the approaching criteria
           If (allokp) Then
              vp = a*d - b*ee
              allokp = allokp .And. vp < 0D0
           End If

           !       check the closest approach distance criteria
           If (allokp) Then
              dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)
              allokp = allokp .And. dm2 < cutof2
           End If

           !       check the time criteria
           If (allokp) Then
              tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
              tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
              tmp = 0.5D0*(tc1+tc2)
              allokp = allokp .And. tmp > ft(i) .And. tmp > ft(j)
           End If

           If (allokp .And. tmp<tm) Then
              tm = tmp
              jxa = ii
              jya = jj
              !d                    dgxa(j) = ii * 10d0 * size1
              !d                    dgya(j) = jj * 10d0 * size2
              !d                    dgxa(i) = - dgxa(j)
              !d                    dgya(i) = - dgya(j)
           End If

        End Do
     End Do

     If (tm==tlarge) Then
        allok = .False.
     End If

  End If

  !        check rts cut
  If (allok) Then

     q4 = e(i1)
     q1 = px(i1)
     q2 = py(i1)
     q3 = pz(i1)

     r4 = e(i2)
     r1 = px(i2)
     r2 = py(i2)
     r3 = pz(i2)

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (h>0D0) Then
     t1 = tm
     t2 = tm
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco7

Subroutine isco8(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok, allokp

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations

  tm = tlarge

  If (allok) Then
     Do ii = -1, 1
        Do jj = -1, 1

           allokp = .True.

           i1 = i
           i2 = j

           p4 = ft(j) - ft(i)
           p1 = gx(j) - gx(i)
           p2 = gy(j) - gy(i)
           p3 = gz(j) - gz(i)

           p1 = p1 + ii*10D0*size1
           p2 = p2 + jj*10D0*size2

           q4 = e(i)
           q1 = px(i)
           q2 = py(i)
           q3 = pz(i)

           r4 = e(j)
           r1 = px(j)
           r2 = py(j)
           r3 = pz(j)

           a = p4*q4 - p1*q1 - p2*q2 - p3*q3
           b = p4*r4 - p1*r1 - p2*r2 - p3*r3
           c = q4*q4 - q1*q1 - q2*q2 - q3*q3
           d = r4*r4 - r1*r1 - r2*r2 - r3*r3
           ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
           f = p4*p4 - p1*p1 - p2*p2 - p3*p3

           !       make sure particle 2 formed early
           h = a + b
           If (h>0D0) Then
              g = a
              a = -b
              b = -g
              g = c
              c = d
              d = g
              i1 = j
              i2 = i
           End If

           !       check the approaching criteria
           If (allokp) Then
              vp = a*d - b*ee
              allokp = allokp .And. vp < 0D0
           End If

           !       check the closest approach distance criteria
           If (allokp) Then
              dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)
              allokp = allokp .And. dm2 < cutof2
           End If

           !       check the time criteria
           If (allokp) Then
              tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
              tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
              If (iordsc==20) Then
                 tmp = min(tc1, tc2)
              Else If (iordsc==21) Then
                 tmp = 0.5D0*(tc1+tc2)
              Else
                 tmp = max(tc1, tc2)
              End If
              allokp = allokp .And. tmp > ft(i) .And. tmp > ft(j)
           End If

           If (allokp .And. tmp<tm) Then
              tm = tmp
              jxa = ii
              jya = jj
              ha = h
              tc1a = tc1
              tc2a = tc2
              !d                    dgxa(j) = ii * 10d0 * size1
              !d                    dgya(j) = jj * 10d0 * size2
              !d                    dgxa(i) = - dgxa(j)
              !d                    dgya(i) = - dgya(j)
           End If

        End Do
     End Do

     If (tm==tlarge) Then
        allok = .False.
     End If

  End If

  !        check rts cut
  If (allok) Then

     q4 = e(i1)
     q1 = px(i1)
     q2 = py(i1)
     q3 = pz(i1)

     r4 = e(i2)
     r1 = px(i2)
     r2 = py(i2)
     r3 = pz(i2)

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (ha>0D0) Then
     t1 = tc2a
     t2 = tc1a
  Else
     t1 = tc1a
     t2 = tc2a
  End If

  Return
End Subroutine isco8

Subroutine isco9(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok, allokp

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  If (ft(i)>=ft(j)) Then
     i1 = j
     i2 = i
     isign = -1
  Else
     i1 = i
     i2 = j
     isign = 1
  End If

  If (allok) Then
     tm = tlarge

     t1 = ft(i1)
     vx1 = vx(i1)
     vy1 = vy(i1)
     vz1 = vz(i1)

     t2 = ft(i2)

     dvx = vx(i2) - vx1
     dvy = vy(i2) - vy1
     dvz = vz(i2) - vz1

     dt = t2 - t1

     Do ii = -1, 1
        Do jj = -1, 1

           allokp = .True.

           dx = gx(i2) - gx(i1) - vx1*dt
           dy = gy(i2) - gy(i1) - vy1*dt
           dz = gz(i2) - gz(i1) - vz1*dt

           dx = dx + ii*10D0*size1
           dy = dy + jj*10D0*size2

           vp = dvx*dx + dvy*dy + dvz*dz

           allokp = allokp .And. vp < 0D0

           If (allokp) Then

              v2 = dvx*dvx + dvy*dvy + dvz*dvz

              If (v2==0D0) Then
                 tmp = tlarge
              Else
                 tmp = t2 - vp/v2
              End If

              !       note now tm is the absolute time

              allokp = allokp .And. tmp > t1 .And. tmp > t2

           End If

           If (allokp) Then

              dgx = dx - dvx*t2
              dgy = dy - dvy*t2
              dgz = dz - dvz*t2

              dm2 = -v2*tmp**2 + dgx*dgx + dgy*dgy + dgz*dgz

              allokp = allokp .And. dm2 < cutof2

           End If

           If (allokp .And. tmp<tm) Then
              tm = tmp
              jxa = isign*ii
              jya = isign*jj
           End If

        End Do
     End Do

     If (tm==tlarge) Then
        allok = .False.
     End If
  End If

  If (allok) Then

     e1 = e(i1)
     px1 = px(i1)
     py1 = py(i1)
     pz1 = pz(i1)
     e2 = e(i2)
     px2 = px(i2)
     py2 = py(i2)
     pz2 = pz(i2)

     rts2 = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco9

Subroutine isco10(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok, allokp

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations

  tm = tlarge

  If (allok) Then
     Do ii = -1, 1
        Do jj = -1, 1
           Do kk = -1, 1
              allokp = .True.

              i1 = i
              i2 = j

              p4 = ft(j) - ft(i)
              p1 = gx(j) - gx(i)
              p2 = gy(j) - gy(i)
              p3 = gz(j) - gz(i)

              p1 = p1 + ii*10D0*size1
              p2 = p2 + jj*10D0*size2
              p3 = p3 + kk*10D0*size3

              q4 = e(i)
              q1 = px(i)
              q2 = py(i)
              q3 = pz(i)

              r4 = e(j)
              r1 = px(j)
              r2 = py(j)
              r3 = pz(j)

              a = p4*q4 - p1*q1 - p2*q2 - p3*q3
              b = p4*r4 - p1*r1 - p2*r2 - p3*r3
              c = q4*q4 - q1*q1 - q2*q2 - q3*q3
              d = r4*r4 - r1*r1 - r2*r2 - r3*r3
              ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
              f = p4*p4 - p1*p1 - p2*p2 - p3*p3

              !       make sure particle 2 formed early
              h = a + b
              If (h>0D0) Then
                 g = a
                 a = -b
                 b = -g
                 g = c
                 c = d
                 d = g
                 i1 = j
                 i2 = i
              End If

              !       check the approaching criteria
              If (allokp) Then
                 vp = a*d - b*ee
                 allokp = allokp .And. vp < 0D0
              End If

              !       check the closest approach distance criteria
              If (allokp) Then
                 dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)
                 allokp = allokp .And. dm2 < cutof2
              End If

              !       check the time criteria
              If (allokp) Then
                 tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
                 tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
                 tmp = 0.5D0*(tc1+tc2)
                 allokp = allokp .And. tmp > ft(i) .And. tmp > ft(j)
              End If

              If (allokp .And. tmp<tm) Then
                 tm = tmp
                 jxa = ii
                 jya = jj
                 jza = kk
                 !d                    dgxa(j) = ii * 10d0 * size1
                 !d                    dgya(j) = jj * 10d0 * size2
                 !d                    dgxa(i) = - dgxa(j)
                 !d                    dgya(i) = - dgya(j)
              End If

           End Do
        End Do
     End Do

     If (tm==tlarge) Then
        allok = .False.
     End If

  End If

  !        check rts cut
  If (allok) Then

     q4 = e(i1)
     q1 = px(i1)
     q2 = py(i1)
     q3 = pz(i1)

     r4 = e(i2)
     r1 = px(i2)
     r2 = py(i2)
     r3 = pz(i2)

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (h>0D0) Then
     t1 = tm
     t2 = tm
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco10

Subroutine isco11(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok, allokp

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  !       set up numbers for later calculations

  tm = tlarge

  If (allok) Then
     Do ii = -1, 1
        Do jj = -1, 1
           Do kk = -1, 1

              allokp = .True.

              i1 = i
              i2 = j

              p4 = ft(j) - ft(i)
              p1 = gx(j) - gx(i)
              p2 = gy(j) - gy(i)
              p3 = gz(j) - gz(i)

              p1 = p1 + ii*10D0*size1
              p2 = p2 + jj*10D0*size2
              p3 = p3 + kk*10D0*size3

              q4 = e(i)
              q1 = px(i)
              q2 = py(i)
              q3 = pz(i)

              r4 = e(j)
              r1 = px(j)
              r2 = py(j)
              r3 = pz(j)

              a = p4*q4 - p1*q1 - p2*q2 - p3*q3
              b = p4*r4 - p1*r1 - p2*r2 - p3*r3
              c = q4*q4 - q1*q1 - q2*q2 - q3*q3
              d = r4*r4 - r1*r1 - r2*r2 - r3*r3
              ee = q4*r4 - q1*r1 - q2*r2 - q3*r3
              f = p4*p4 - p1*p1 - p2*p2 - p3*p3

              !       make sure particle 2 formed early
              h = a + b
              If (h>0D0) Then
                 g = a
                 a = -b
                 b = -g
                 g = c
                 c = d
                 d = g
                 i1 = j
                 i2 = i
              End If

              !       check the approaching criteria
              If (allokp) Then
                 vp = a*d - b*ee
                 allokp = allokp .And. vp < 0D0
              End If

              !       check the closest approach distance criteria
              If (allokp) Then
                 dm2 = -f - (a**2*d+b**2*c-2D0*a*b*ee)/(ee**2-c*d)
                 allokp = allokp .And. dm2 < cutof2
              End If

              !       check the time criteria
              If (allokp) Then
                 tc1 = ft(i1) - e(i1)*(a*d-b*ee)/(ee**2-c*d)
                 tc2 = ft(i2) + e(i2)*(b*c-a*ee)/(ee**2-c*d)
                 If (iordsc==20) Then
                    tmp = min(tc1, tc2)
                 Else If (iordsc==21) Then
                    tmp = 0.5D0*(tc1+tc2)
                 Else
                    tmp = max(tc1, tc2)
                 End If
                 allokp = allokp .And. tmp > ft(i) .And. tmp > ft(j)
              End If

              If (allokp .And. tmp<tm) Then
                 tm = tmp
                 jxa = ii
                 jya = jj
                 jza = kk
                 ha = h
                 tc1a = tc1
                 tc2a = tc2
                 !d                    dgxa(j) = ii * 10d0 * size1
                 !d                    dgya(j) = jj * 10d0 * size2
                 !d                    dgxa(i) = - dgxa(j)
                 !d                    dgya(i) = - dgya(j)
              End If

           End Do
        End Do
     End Do

     If (tm==tlarge) Then
        allok = .False.
     End If

  End If

  !        check rts cut
  If (allok) Then

     q4 = e(i1)
     q1 = px(i1)
     q2 = py(i1)
     q3 = pz(i1)

     r4 = e(i2)
     r1 = px(i2)
     r2 = py(i2)
     r3 = pz(i2)

     rts2 = (q4+r4)**2 - (q1+r1)**2 - (q2+r2)**2 - (q3+r3)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else If (ha>0D0) Then
     t1 = tc2a
     t2 = tc1a
  Else
     t1 = tc1a
     t2 = tc2a
  End If

  Return
End Subroutine isco11

Subroutine isco12(i, j, allok, tm, t1, t2)
  !       this subroutine is used to decide whether there is a collision between
  !       particle i and j, if there is one allok=1, and tm gives the
  !       collision time, t1 the collision time for i,
  !       t2 the collision time for j

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok, allokp

  !       preventing consecutive collisions
  allok = last(i) /= j .Or. last(j) /= i

  If (ft(i)>=ft(j)) Then
     i1 = j
     i2 = i
     isign = -1
  Else
     i1 = i
     i2 = j
     isign = 1
  End If

  If (allok) Then
     tm = tlarge

     t1 = ft(i1)
     vx1 = vx(i1)
     vy1 = vy(i1)
     vz1 = vz(i1)

     t2 = ft(i2)

     dvx = vx(i2) - vx1
     dvy = vy(i2) - vy1
     dvz = vz(i2) - vz1

     dt = t2 - t1

     Do ii = -1, 1
        Do jj = -1, 1
           Do kk = -1, 1

              allokp = .True.

              dx = gx(i2) - gx(i1) - vx1*dt
              dy = gy(i2) - gy(i1) - vy1*dt
              dz = gz(i2) - gz(i1) - vz1*dt

              dx = dx + ii*10D0*size1
              dy = dy + jj*10D0*size2
              dz = dz + kk*10D0*size3

              vp = dvx*dx + dvy*dy + dvz*dz

              allokp = allokp .And. vp < 0D0

              If (allokp) Then

                 v2 = dvx*dvx + dvy*dvy + dvz*dvz

                 If (v2==0D0) Then
                    tmp = tlarge
                 Else
                    tmp = t2 - vp/v2
                 End If

                 !       note now tm is the absolute time

                 allokp = allokp .And. tmp > t1 .And. tmp > t2

              End If

              If (allokp) Then

                 dgx = dx - dvx*t2
                 dgy = dy - dvy*t2
                 dgz = dz - dvz*t2

                 dm2 = -v2*tmp**2 + dgx*dgx + dgy*dgy + dgz*dgz

                 allokp = allokp .And. dm2 < cutof2

              End If

              If (allokp .And. tmp<tm) Then
                 tm = tmp
                 jxa = isign*ii
                 jya = isign*jj
                 jza = isign*kk
              End If

           End Do
        End Do
     End Do

     If (tm==tlarge) Then
        allok = .False.
     End If
  End If

  If (allok) Then

     e1 = e(i1)
     px1 = px(i1)
     py1 = py(i1)
     pz1 = pz(i1)
     e2 = e(i2)
     px2 = px(i2)
     py2 = py(i2)
     pz2 = pz(i2)

     rts2 = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2

     allok = allok .And. rts2 > rscut2
  End If

  If (.Not. allok) Then
     tm = tlarge
     t1 = tlarge
     t2 = tlarge
  Else
     t1 = tm
     t2 = tm
  End If

  Return
End Subroutine isco12

Subroutine reor(t, tmin, j, last0)
  !       this subroutine is used to fix ct(i) when tm is greater than ct(i)
  !       next(i) is last1 or last2

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  !d        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
  !c      SAVE /ilist5/
  Save

  icels0 = icels(j)

  i1 = icels0/10000
  i2 = (icels0-i1*10000)/100
  i3 = icels0 - i1*10000 - i2*100

  Call wallc(j, i1, i2, i3, t, tmin1)

  If (tmin<=tmin1) Then
     nc = last0
  Else
     tmin = tmin1
     nc = 0
  End If

  If (iconfg==3 .Or. iconfg==5) Then
     Call chcell(j, i1, i2, i3, last0, t, tmin, nc)
  Else
     If (i1==11 .And. i2==11 .And. i3==11) Then
        Call chout(j, last0, t, tmin, nc)
     Else
        If (iconfg==1) Then
           Call chin1(j, i1, i2, i3, last0, t, tmin, nc)
        Else If (iconfg==2) Then
           Call chin2(j, i1, i2, i3, last0, t, tmin, nc)
        Else If (iconfg==4) Then
           Call chin3(j, i1, i2, i3, last0, t, tmin, nc)
        End If
     End If
  End If

  Call fixtim(j, t, tmin1, tmin, nc)

  Return
End Subroutine reor

Subroutine chout(l, last0, t, tmin, nc)
  !       this subroutine is used to check the surface when the colliding
  !       particle is outside the cube

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Save

  m1 = 11
  m2 = 11
  m3 = 11
  Call chcell(l, m1, m2, m3, last0, t, tmin, nc)

  Do i = 1, 10
     Do j = 1, 10
        Do k = 1, 10
           If (i==1 .Or. i==10 .Or. j==1 .Or. j==10 .Or. k==1 .Or. k==10) Call chcell(l, i, j, k, last0, t, tmin, nc)
        End Do
     End Do
  End Do

  Return
End Subroutine chout

Subroutine chin1(l, i1, i2, i3, last0, t, tmin, nc)
  !       this subroutine is used to check collisions for particle inside
  !       the cube
  !       and update the afftected particles through chcell

  Implicit Double Precision (A-H, O-Z)
  Save

  !       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0

  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (i>=1 .And. i<=10 .And. j>=1 .And. j<=10 .And. k>=1 .And. k<=10) Then
              Call chcell(l, i, j, k, last0, t, tmin, nc)
           Else If (itest==0) Then
              m1 = 11
              m2 = 11
              m3 = 11
              Call chcell(l, m1, m2, m3, last0, t, tmin, nc)
              itest = 1
           End If
        End Do
     End Do
  End Do

  Return
End Subroutine chin1

Subroutine chin2(l, i1, i2, i3, last0, t, tmin, nc)
  !       this subroutine is used to check collisions for particle inside
  !       the cube
  !       and update the afftected particles through chcell

  Implicit Double Precision (A-H, O-Z)
  Save

  !       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0

  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           ia = i
           ib = j
           ic = k
           If (k>=1 .And. k<=10) Then
              If (i==0) ia = 10
              If (i==11) ia = 1
              If (j==0) ib = 10
              If (j==11) ib = 1
              Call chcell(l, ia, ib, ic, last0, t, tmin, nc)
           End If
        End Do
     End Do
  End Do

  Return
End Subroutine chin2

Subroutine chin3(l, i1, i2, i3, last0, t, tmin, nc)
  !       this subroutine is used to check collisions for particle inside
  !       the cube
  !       and update the afftected particles through chcell

  Implicit Double Precision (A-H, O-Z)
  Save

  !       itest is a flag to make sure the 111111 cell is checked only once
  itest = 0

  Do i = i1 - 1, i1 + 1
     Do j = i2 - 1, i2 + 1
        Do k = i3 - 1, i3 + 1
           If (i==0) Then
              ia = 10
           Else If (i==11) Then
              ia = 1
           Else
              ia = i
           End If
           If (j==0) Then
              ib = 10
           Else If (j==11) Then
              ib = 1
           Else
              ib = j
           End If
           If (k==0) Then
              ic = 10
           Else If (k==11) Then
              ic = 1
           Else
              ic = k
           End If
           Call chcell(l, ia, ib, ic, last0, t, tmin, nc)
        End Do
     End Do
  End Do

  Return
End Subroutine chin3

Subroutine chcell(il, i1, i2, i3, last0, t, tmin, nc)
  !       this program is used to check through all the particles, except last0
  !       in the cell (i1,i2,i3) and see if we can get a particle collision
  !       with time less than the original input tmin ( the collision time of
  !       il with the wall
  !       last0 cas be set to 0 if we don't want to exclude last0

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist2/icell, icel(10, 10, 10)
  !c      SAVE /ilist2/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Save

  If (iconfg==3 .Or. iconfg==5) Then
     jj = ichkpt
     Do j = 1, jj
        !     10/24/02 get rid of argument usage mismatch in mintm():
        jmintm = j
        If (j/=il .And. j/=last0) Call mintm(il, jmintm, tmin, nc)
        !     &          call mintm(il, j, tmin, nc)

     End Do
     Return
  End If

  !       set l
  If (i1==11 .And. i2==11 .And. i3==11) Then
     l = icell
  Else
     l = icel(i1, i2, i3)
  End If

  If (l==0) Return

  j = nic(l)

  !       if there is only one particle
  If (j==0) Then

     !       if it's not il or last0,when last is not wall
     If (l==il .Or. l==last0) Return
     Call mintm(il, l, tmin, nc)

     !       if there are many particles
  Else
     If (l/=il .And. l/=last0) Call mintm(il, l, tmin, nc)
     Do While (j/=l)
        If (j/=il .And. j/=last0) Call mintm(il, j, tmin, nc)
        j = nic(j)
     End Do
  End If

  Return
End Subroutine chcell

Subroutine mintm(i, j, tmin, nc)
  !       this subroutine is used to check whether particle j has smaller
  !       collision time with particle i than other particles
  !       or in other words, update next(i)

  !       input i,j,tmin,nc
  !       output tmin,nc

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /aurec1/jxa, jya, jza
  !c      SAVE /aurec1/
  Common /aurec2/dgxa(maxptn), dgya(maxptn), dgza(maxptn)
  !c      SAVE /aurec2/
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Save

  Logical allok

  Call isco(i, j, allok, tm, t1, t2)

  If (allok .And. tm<tmin) Then
     tmin = tm
     ct(i) = t1
     nc = j
     If (iconfg==3 .Or. iconfg==5) Then
        dgxa(i) = -jxa*10D0*size1
        dgya(i) = -jya*10D0*size2
        If (iconfg==5) Then
           dgza(i) = -jza*10D0*size3
        End If
     End If
  End If

  Return
End Subroutine mintm

!*****************************************************************************
!*****************************************************************************

Subroutine zpca1

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Save

  If (mod(ictype,2)==0) Then
     Call zpca1a(iscat)
     Call zpca1a(jscat)
     !lin-5/2009 ctest off v2 for parton:
     !           call flowp(1)
  End If

  Return
End Subroutine zpca1

Subroutine zpca1a(i)

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec3/gxs(maxptn), gys(maxptn), gzs(maxptn), fts(maxptn), pxs(maxptn), pys(maxptn), pzs(maxptn), es(maxptn), xmasss(maxptn), ityps(maxptn)
  !c      SAVE /prec3/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /prec6/etas(maxptn), raps(maxptn), taus(maxptn)
  !c      SAVE /prec6/
  Common /ana1/ts(12)
  !c      SAVE /ana1/
  Save

  If (iconfg==1) Then
     t1 = fts(i)
     t2 = ft(i)
     ipic = 11
  Else If (iconfg==2 .Or. iconfg==3) Then
     !d           t1 = fts(i)
     !d           t2 = ft(i)
     t1 = taus(i)
     t2 = tau(i)
     ipic = 12
  Else If (iconfg==4 .Or. iconfg==5) Then
     t1 = fts(i)
     t2 = ft(i)
     ipic = 12
  End If

  If (iconfg<=3) Then
     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           rapi = raps(i)
           !     7/20/01:
           !                 et = sqrt(pxs(i) ** 2 + pys(i) ** 2 + xmp ** 2)
           et = dsqrt(pxs(i)**2+pys(i)**2+xmp**2)
           Call zpca1b(rapi, et, ian)
        End If
     End Do
  Else
     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           p0 = es(i)
           p1 = pxs(i)
           p2 = pys(i)
           p3 = pzs(i)
           Call zpca1c(p0, p1, p2, p3, ian)
        End If
     End Do
  End If

  Return
End Subroutine zpca1a

Subroutine zpca1b(rapi, et, ian)

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para6/centy
  !c      SAVE /para6/
  Common /ilist6/t, iopern, icolln
  !c      SAVE /ilist6/
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  !c      SAVE /ana2/
  Save

  If (rapi>centy-0.5D0 .And. rapi<centy+0.5D0) Then
     det2(ian) = det2(ian) + et
     dn2(ian) = dn2(ian) + 1D0
     !dtrans
     If (ian==10) Then
        !d              write (10, *) t, det2(ian)
     End If
     If (ian==11) Then
        !d              write (11, *) t, det2(ian)
     End If
     If (ian==12) Then
        !d              write (12, *) t, det2(ian)
     End If
     !dtransend
     If (rapi>centy-0.25D0 .And. rapi<centy+0.25D0) Then
        det1(ian) = det1(ian) + et
        dn1(ian) = dn1(ian) + 1D0
        If (rapi>centy-0.1D0 .And. rapi<centy+0.1D0) Then
           det(ian) = det(ian) + et
           dn(ian) = dn(ian) + 1D0
        End If
     End If
  End If

  Return
End Subroutine zpca1b

Subroutine zpca1c(p0, p1, p2, p3, ian)

  Implicit Double Precision (A-H, O-Z)

  Common /ana3/em(4, 4, 12)
  !c      SAVE /ana3/

  Dimension en(4)
  Save

  en(1) = p0
  en(2) = p1
  en(3) = p2
  en(4) = p3

  Do i = 1, 4
     Do j = 1, 4
        em(i, j, ian) = em(i, j, ian) + en(i)*en(j)/p0
     End Do
  End Do

  Return
End Subroutine zpca1c

!*****************************************************************************

Subroutine zpca2

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /para7/ioscar, nsmbbbar, nsmmeson
  !c      SAVE /para7/
  Common /ilist6/t, iopern, icolln
  !c      SAVE /ilist6/
  Common /rndm1/number
  !c      SAVE /rndm1/
  Common /rndm2/iff
  !c      SAVE /rndm2/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Common /arevt/iaevt, iarun, miss
  !c      SAVE /AREVT/
  Save

  If (iconfg<=3) Then
     Call zpca2a
  Else
     Call zpca2b
  End If

  If (ioscar==1) Then
     Call zpca2c
  End If

  !bzdbg2/17/99
  !        write (25, *) 'Event', nsevt - 1 + ievt,
  !    &         ', run', isbrun,
  !        WRITE (25, *) ' Event ', IAEVT, ', run ', IARUN,
  !     &     ',\n\t number of operations = ', iopern,
  !     &     ',\n\t number of collisions between particles = ',
  !     &         icolln,
  !     &     ',\n\t freezeout time=', t,
  !     &     ',\n\t ending at the ', number, 'th random number',
  !     &     ',\n\t ending collision iff=', iff
  Write (25, *) ' Event ', iaevt, ', run ', iarun
  Write (25, *) '    number of operations = ', iopern
  Write (25, *) '    number of collisions between particles = ', icolln
  Write (25, *) '    freezeout time=', t
  Write (25, *) '    ending at the ', number, 'th random number'
  Write (25, *) '    ending collision iff=', iff

  Return
End Subroutine zpca2

Subroutine zpca2a

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /para1/mul
  !c      SAVE /para1/
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  !c      SAVE /para2/
  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Common /para6/centy
  !c      SAVE /para6/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Common /ilist6/t, iopern, icolln
  !c      SAVE /ilist6/
  Common /rndm1/number
  !c      SAVE /rndm1/
  Common /rndm2/iff
  !c      SAVE /rndm2/
  Common /rndm3/iseedp
  !c      SAVE /rndm3/
  Common /ana1/ts(12)
  !c      SAVE /ana1/
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  !c      SAVE /ana2/
  Common /ana4/fdetdy(24), fdndy(24), fdndpt(12)
  !c      SAVE /ana4/
  Save

  Do i = 1, ichkpt
     rapi = rap(i)
     !     7/20/01:
     !           et = sqrt(px(i) ** 2 + py(i) ** 2 + xmp ** 2)
     et = dsqrt(px(i)**2+py(i)**2+xmp**2)

     Do j = 1, 24
        If (rapi>j+centy-13D0 .And. rapi<j+centy-12D0) Then
           fdetdy(j) = fdetdy(j) + et
           fdndy(j) = fdndy(j) + 1D0
        End If
     End Do

     Do j = 1, 12
        If (et>0.5D0*(j-1) .And. et<0.5D0*j) Then
           fdndpt(j) = fdndpt(j) + 1D0
        End If
     End Do

     If (iconfg==1) Then
        t1 = ft(i)
        t2 = tlarge
        ipic = 11
     Else
        t1 = tau(i)
        t2 = tlarge
        ipic = 12
     End If

     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           Call zpca1b(rapi, et, ian)
        End If
     End Do

     If (iconfg==1) Then
        Call zpca1b(rapi, et, 12)
     End If
  End Do

  Do ian = 1, 12
     If (dn(ian)==0D0 .Or. dn1(ian)==0D0 .Or. dn2(ian)==0D0) Then
        !lin-9/2012 suppress output:
        !              print *, 'event=', ievt
        !              print *, 'dn(', ian, ')=', dn(ian), 'dn1(', ian,
        !     &           ')=', dn1(ian), 'dn2(', ian, ')=', dn2(ian)
     End If
     detdy(ian) = detdy(ian) + det(ian)
     If (dn(ian)/=0) Then
        detdn(ian) = detdn(ian) + det(ian)/dn(ian)
     End If
     dndy(ian) = dndy(ian) + dn(ian)
     detdy1(ian) = detdy1(ian) + det1(ian)
     If (dn1(ian)/=0) Then
        detdn1(ian) = detdn1(ian) + det1(ian)/dn1(ian)
     End If
     dndy1(ian) = dndy1(ian) + dn1(ian)
     detdy2(ian) = detdy2(ian) + det2(ian)
     If (dn2(ian)/=0) Then
        detdn2(ian) = detdn2(ian) + det2(ian)/dn2(ian)
     End If
     dndy2(ian) = dndy2(ian) + dn2(ian)
  End Do

  Return
End Subroutine zpca2a

Subroutine zpca2b

  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)

  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Common /ilist4/ifmpt, ichkpt, indx(maxptn)
  !c      SAVE /ilist4/
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Common /ana1/ts(12)
  !c      SAVE /ana1/
  Save

  Do i = 1, ichkpt
     t1 = ft(i)
     t2 = tlarge
     ipic = 12

     Do ian = 1, ipic
        If (t1<=ts(ian) .And. t2>ts(ian)) Then
           p0 = e(i)
           p1 = px(i)
           p2 = py(i)
           p3 = pz(i)
           Call zpca1c(p0, p1, p2, p3, ian)
        End If
     End Do
  End Do

  Return
End Subroutine zpca2b

Subroutine zpca2c

  Implicit Double Precision (A-H, O-Z)

  Character *8 code, versn
  Character *4 reffra
  Integer aproj, zproj, atarg, ztarg, event
  Parameter (maxptn=400001)

  Common /para1/mul
  !c      SAVE /para1/
  Common /prec2/gx(maxptn), gy(maxptn), gz(maxptn), ft(maxptn), px(maxptn), py(maxptn), pz(maxptn), e(maxptn), xmass(maxptn), ityp(maxptn)
  !c      SAVE /prec2/
  Save
  Data nff/0/

  !       file header
  If (nff==0) Then
     Write (26, 101) 'OSCAR1997A'
     Write (26, 101) 'final_id_p_x'
     code = 'ZPC'
     versn = '1.0.1'
     aproj = -1
     zproj = -1
     atarg = -1
     ztarg = -1
     reffra = 'cm'
     ebeam = 0D0
     ntestp = 1
     Write (26, 102) code, versn, aproj, zproj, atarg, ztarg, reffra, ebeam, ntestp
     nff = 1
     event = 1
     bimp = 0D0
     phi = 0D0
  End If

  !       comment

  !       event header
  Write (26, 103) event, mul, bimp, phi

  !       particles
  Do i = 1, mul
     Write (26, 104) i, ityp(i), px(i), py(i), pz(i), e(i), xmass(i), gx(i), gy(i), gz(i), ft(i)
  End Do

  event = event + 1

  Return

101 Format (A12)
102 Format (2(A8,2X), '(', I3, ',', I6, ')+(', I3, ',', I6, ')', 2X, A4, 2X, E10.4, 2X, I8)
103 Format (I10, 2X, I10, 2X, F8.3, 2X, F8.3)
104 Format (I10, 2X, I10, 2X, 9(E12.6,2X))
End Subroutine zpca2c

!*****************************************************************************

Subroutine zpcou

  Implicit Double Precision (A-H, O-Z)

  Common /para5/iconfg, iordsc
  !c      SAVE /para5/
  Save

  If (iconfg<=3) Then
     Call zpcou1
  Else
     Call zpcou2
  End If

  Return
End Subroutine zpcou

Subroutine zpcou1

  Implicit Double Precision (A-H, O-Z)

  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Common /ana1/ts(12)
  !c      SAVE /ana1/
  Common /ana2/det(12), dn(12), detdy(12), detdn(12), dndy(12), det1(12), dn1(12), detdy1(12), detdn1(12), dndy1(12), det2(12), dn2(12), detdy2(12), detdn2(12), dndy2(12)
  !c      SAVE /ana2/
  Common /ana4/fdetdy(24), fdndy(24), fdndpt(12)
  !c      SAVE /ana4/
  Save
  !
  dpt = 0.5D0
  dy2 = 1D0
  dy1 = 0.5D0
  dy = 0.2D0
  ntotal = nevnt*nsbrun
  !
  Return
End Subroutine zpcou1

Subroutine zpcou2

  Implicit Double Precision (A-H, O-Z)

  Common /para3/nsevt, nevnt, nsbrun, ievt, isbrun
  !c      SAVE /para3/
  Common /ilist3/size1, size2, size3, v1, v2, v3, size
  !c      SAVE /ilist3/
  Common /ana1/ts(12)
  !c      SAVE /ana1/
  Common /ana3/em(4, 4, 12)
  !c      SAVE /ana3/
  Save
  !
  Open (28, File='ana4/em.dat', Status='unknown')
  vol = 1000.D0*size1*size2*size3
  ntotal = nevnt*nsbrun

  Do ian = 1, 12
     Write (28, *) '*** for time ', ts(ian), 'fm(s)'
     Do i = 1, 4
        Write (28, *) em(i, 1, ian)/vol/ntotal, em(i, 2, ian)/vol/ntotal, em(i, 3, ian)/vol/ntotal, em(i, 4, ian)/vol/ntotal
     End Do
  End Do

  Return
End Subroutine zpcou2

!*****************************************************************************

Subroutine lorenz(energy, px, py, pz, bex, bey, bez)

  !     add in a cut for beta2 to prevent gam to be nan (infinity)

  Implicit Double Precision (A-H, O-Z)

  Common /lor/enenew, pxnew, pynew, pznew
  !c      SAVE /lor/
  Save

  beta2 = bex**2 + bey**2 + bez**2
  If (beta2==0D0) Then
     enenew = energy
     pxnew = px
     pynew = py
     pznew = pz
  Else
     If (beta2>0.999999999999999D0) Then
        beta2 = 0.999999999999999D0
        Print *, 'beta2=0.999999999999999'
     End If
     !lin-7/20/01:
     !         gam = 1.d0 / sqrt(1.d0 - beta2)
     gam = 1.D0/dsqrt(1.D0-beta2)
     enenew = gam*(energy-bex*px-bey*py-bez*pz)
     pxnew = -gam*bex*energy + (1.D0+(gam-1.D0)*bex**2/beta2)*px + (gam-1.D0)*bex*bey/beta2*py + (gam-1.D0)*bex*bez/beta2*pz
     pynew = -gam*bey*energy + (gam-1.D0)*bex*bey/beta2*px + (1.D0+(gam-1.D0)*bey**2/beta2)*py + (gam-1.D0)*bey*bez/beta2*pz
     pznew = -gam*bez*energy + (gam-1.D0)*bex*bez/beta2*px + (gam-1.D0)*bey*bez/beta2*py + (1.D0+(gam-1.D0)*bez**2/beta2)*pz
  End If

  Return
End Subroutine lorenz

Subroutine index1(n, m, arrin, indx)
  !     indexes the first m elements of ARRIN of length n, i.e., outputs INDX
  !     such that ARRIN(INDEX(J)) is in ascending order for J=1,...,m

  Implicit Double Precision (A-H, O-Z)

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

End Subroutine index1


Double Precision Function ftime1(iseed)

  !       this program is used to generate formation time
  !       the calling program needs a common /par1/
  !       and declare external ftime1

  !lin-8/19/02
  Implicit Double Precision (A-H, O-Z)

  External ran1

  Parameter (hbarc=0.197327054D0)

  Common /par1/formt
  !c      SAVE /par1/
  Save

  aa = hbarc/formt

  !lin7/20/01:
  !        ftime1 = aa * sqrt(1d0 / ran1(iseed) - 1d0)
  ftime1 = aa*dsqrt(1D0/ran1(iseed)-1D0)
  Return
End Function ftime1


Subroutine cropro(vx1, vy1, vz1, vx2, vy2, vz2)

  !     this subroutine is used to calculate the cross product of
  !     (vx1,vy1,vz1) and (vx2,vy2,vz2) and get the result (vx3,vy3,vz3)
  !     and put the vector into common /cprod/

  Implicit Double Precision (A-H, O-Z)

  Common /cprod/vx3, vy3, vz3
  !c      SAVE /cprod/
  Save

  vx3 = vy1*vz2 - vz1*vy2
  vy3 = vz1*vx2 - vx1*vz2
  vz3 = vx1*vy2 - vy1*vx2

  Return
End Subroutine cropro

Subroutine xnormv(vx, vy, vz)

  !      this subroutine is used to get a normalized vector

  Implicit Double Precision (A-H, O-Z)
  Save

  !lin-7/20/01:
  !      vv = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
  vv = dsqrt(vx**2+vy**2+vz**2)
  vx = vx/vv
  vy = vy/vv
  vz = vz/vv

  Return
End Subroutine xnormv

!bz1/29/99
!      subroutine rotate(xn1, xn2, xn3, theta, v1, v2, v3)
Subroutine zprota(xn1, xn2, xn3, theta, v1, v2, v3)
  !bz1/29/99end

  !     this subroutine is used to rotate the vector (v1,v2,v3) by an angle theta
  !     around the unit vector (xn1, xn2, xn3)

  Implicit Double Precision (A-H, O-Z)
  Save

  vx = v1
  vy = v2
  vz = v3
  c = cos(theta)
  omc = 1D0 - c
  s = sin(theta)
  a11 = xn1**2*omc + c
  a12 = xn1*xn2*omc - s*xn3
  a13 = xn1*xn3*omc + s*xn2
  a21 = xn1*xn2*omc + s*xn3
  a22 = xn2**2*omc + c
  a23 = xn2*xn3*omc - s*xn1
  a31 = xn1*xn3*omc - s*xn2
  a32 = xn3*xn2*omc + s*xn1
  a33 = xn3**2*omc + c
  v1 = vx*a11 + vy*a12 + vz*a13
  v2 = vx*a21 + vy*a22 + vz*a23
  v3 = vx*a31 + vy*a32 + vz*a33

  Return
End Subroutine zprota

Double Precision Function ran1(idum)

  !     return a uniform random deviate between 0.0 and 1.0. set idum to
  !     any negative value to initialize or reinitialize the sequence.

  Implicit Double Precision (A-H, O-Z)

  Dimension r(97)

  Common /rndm1/number
  !c      SAVE /rndm1/
  Parameter (m1=259200, ia1=7141, ic1=54773, rm1=1D0/m1)
  Parameter (m2=134456, ia2=8121, ic2=28411, rm2=1D0/m2)
  Parameter (m3=243000, ia3=4561, ic3=51349)
  !lin-6/23/00 save ix1-3:
  !lin-10/30/02 r unsaved, causing wrong values for ran1 when compiled with f77:
  !c      SAVE ix1,ix2,ix3,r
  Save
  Data iff/0/

  If (idum<0 .Or. iff==0) Then
     iff = 1
     ix1 = mod(ic1-idum, m1)
     ix1 = mod(ia1*ix1+ic1, m1)
     ix2 = mod(ix1, m2)
     ix1 = mod(ia1*ix1+ic1, m1)
     ix3 = mod(ix1, m3)
     Do j = 1, 97
        ix1 = mod(ia1*ix1+ic1, m1)
        ix2 = mod(ia2*ix2+ic2, m2)
        r(j) = (dble(ix1)+dble(ix2)*rm2)*rm1
     End Do
     idum = 1
  End If
  ix1 = mod(ia1*ix1+ic1, m1)
  ix2 = mod(ia2*ix2+ic2, m2)
  ix3 = mod(ia3*ix3+ic3, m3)
  !lin-7/01/02       j = 1 + (97 * i x 3) / m3
  j = 1 + (97*ix3)/m3
  !lin-4/2008:
  !      if (j .gt. 97 .or. j .lt. 1) pause
  If (j>97 .Or. j<1) Print *, 'In zpc ran1, j<1 or j>97', j
  ran1 = r(j)
  r(j) = (dble(ix1)+dble(ix2)*rm2)*rm1

  !lin-6/23/00 check random number generator:
  number = number + 1
  !      if(number.le.100000) write(99,*) 'number, ran1=', number,ran1

  Return
End Function ran1
