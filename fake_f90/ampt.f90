!>==========================================================
!> \file
!> 真个程序
!>==========================================================
Program ampt
  !
  Double Precision xmp, xmu, alpha, rscut2, cutof2, dshadow
  Double Precision smearp, smearh, dpcoal, drcoal, ecritl
  Character frame*8, proj*8, targ*8
  Character *25 amptvn
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  Common /arout/iout
  Common /arevt/iaevt, iarun, miss
  Common /smearz/smearp, smearh
  Common /rndf77/nseed
  Common /anim/nevent, isoft, isflag, izpc
  !     parton coalescence radii in case of string melting:
  Common /coal/dpcoal, drcoal, ecritl
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  !     initialization value for parton cascade:
  Common /para2/xmp, xmu, alpha, rscut2, cutof2
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /para8/idpert, npertd, idxsec
  Common /rndm3/iseedp
  !     initialization value for hadron cascade:
  Common /run/num
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  Common /oscar1/iap, izp, iat, izt
  Common /oscar2/frame, amptvn
  Common /resdcy/nsav, iksdcy
  !lin-4/2012-6/2009:
  !      common/phidcy/iphidcy
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  !lin-7/2009:
  Common /cmsflag/dshadow, ishadow
  !lin-2/2012 allow random orientation of reaction plane:
  Common /phihj/iphirp, phirp

  External hidata, pydata, ludata, ardata, ppbdat, zpcbdt
  Save
  !****************
  Open (24, File='input.ampt', Status='UNKNOWN')
  Open (12, File='ana/version', Status='UNKNOWN')
  Read (24, *) efrm
  !     format-read characters (for ALPHA compilers):
  Read (24, 111) frame
  Read (24, 111) proj
  Read (24, 111) targ
  Read (24, *) iap
  Read (24, *) izp
  Read (24, *) iat
  Read (24, *) izt
  Read (24, *) nevnt
  Read (24, *) bmin
  Read (24, *) bmax
  !     flag to select default AMPT or string melting:
  Read (24, *) isoft
  !     read initialization value for hadron cascade:
  Read (24, *) ntmax
  Read (24, *) dt
  !     parj(41) and (42) are a and b parameters in Lund string fragmentation:
  Read (24, *) parj(41)
  Read (24, *) parj(42)
  !     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and
  !     increase the net-baryon stopping in rapidity (value HIJING is 1):
  Read (24, *) ipop
  If (ipop==1) ihpr2(11) = 3
  !     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn:
  Read (24, *) parj(5)
  !     shadowing flag in HIJING:
  Read (24, *) ihpr2(6)
  !     quenching flag in HIJING:
  Read (24, *) ihpr2(4)
  !     quenching rate when quenching flag is on (=1.0 GeV/fm):
  Read (24, *) hipr1(14)
  !     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV.
  Read (24, *) hipr1(8)
  !     read initialization value for parton cascade:
  Read (24, *) xmu
  Read (24, *) izpc
  Read (24, *) alpha
  !     quark coalescence radii in momentum and space for string melting:
  Read (24, *) dpcoal
  Read (24, *) drcoal
  !     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0):
  Read (24, *) ihjsed
  !     2 seeds for random number generators in HIJING/hadron cascade and ZPC:
  Read (24, *) nseed
  Read (24, *) iseedp
  Read (24, *) iksdcy
  Read (24, *) iphidcy
  Read (24, *) ipi0dcy
  !     flag for OSCAR output for final partons and hadrons:
  Read (24, *) ioscar
  !lin-5/2008     flag for perturbative treatment of deuterons:
  Read (24, *) idpert
  Read (24, *) npertd
  Read (24, *) idxsec
  !lin-6/2009 To select events that have at least 1 high-Pt minijet parton:
  Read (24, *) pttrig
  Read (24, *) maxmiss
  Read (24, *) ihpr2(2)
  Read (24, *) ihpr2(5)
  !lin-6/2009 To embed a back-to-back q/qbar pair into each event:
  Read (24, *) iembed
  Read (24, *) pxqembd, pyqembd
  Read (24, *) xembd, yembd
  Read (24, *) nsembd, psembd, tmaxembd
  !lin-7/2009 Allow modification of nuclear shadowing:
  Read (24, *) ishadow
  Read (24, *) dshadow
  Read (24, *) iphirp
  !
  Close (24)
  !lin-6/2009 ctest off turn on jet triggering:
  !      IHPR2(3)=1
  !     Trigger Pt of high-pt jets in HIJING:
  !      HIPR1(10)=7.
  !
  If (isoft==1) Then
     amptvn = '1.26t7 (Default)'
  Else If (isoft==4) Then
     amptvn = '2.26t7 (StringMelting)'
  Else
     amptvn = 'Test-Only'
  End If
  ! Write (6, 50) amptvn
  ! Write (12, 50) amptvn
  !     when ihjsed=11: use environment variable at run time for HIJING nseed:
  If (ihjsed==11) Then
     Print *, '# Read in NSEED in HIJING at run time (e.g. 20030819):'
  End If
  Read (*, *) nseedr
  If (ihjsed==11) Then
     nseed = nseedr
  End If
  If (ihjsed==11) Then
     Print *, '#   read in: ', nseed
     Write (12, *) '# Read in NSEED in HIJING at run time:', nseed
  End If
  Close (12)
  !lin-5/2015 an odd number is needed for the random number generator:
  !      if(mod(NSEED,2).eq.0) NSEED=NSEED+1
  nseed = 2*nseed + 1
  !     9/26/03 random number generator for f77 compiler:
  Call srand(nseed)
  !
  !.....turn on warning messages in nohup.out when an event is repeated:
  ihpr2(10) = 1
  !     string formation time:
  arpar1(1) = 0.7
  !     smearp is the smearing halfwidth on parton z0,
  !     set to 0 for now to avoid overflow in eta.
  !     smearh is the smearing halfwidth on string production point z0.
  smearp = 0D0
  iamax = max(iap, iat)
  smearh = 1.2D0*iamax**0.3333D0/(dble(efrm)/2/0.938D0)
  nevent = nevnt
  !
  !     AMPT momentum and space info at freezeout:
  Open (16, File='ana/ampt.dat', Status='UNKNOWN')
  Open (14, File='ana/zpc.dat', Status='UNKNOWN')
  !test off for resonance (phi, K*) studies:
  !      OPEN (17, FILE = 'ana/res-gain.dat', STATUS = 'UNKNOWN')
  !      OPEN (18, FILE = 'ana/res-loss.dat', STATUS = 'UNKNOWN')
  Call hijset(efrm, frame, proj, targ, iap, izp, iat, izt)
  Call artset
  Call inizpc
  !lin-5/2009 ctest off:
  !      call flowp(0)
  !      call flowh0(NEVNT,0)
  !      call iniflw(NEVNT,0)
  !      call frztm(NEVNT,0)
  !
  Do j = 1, nevnt
     iaevt = j
     Do k = 1, num
        iarun = k
        If (iaevt==nevnt .And. iarun==num) Then
           iout = 1
        End If
        Print *, ' EVENT ', j, ', RUN ', k
        imiss = 0
100     Call hijing(frame, bmin, bmax)
        iaint2(1) = natt

        !lin-6/2009 ctest off
        If (j==-2) Then
           Write (98, *) hipr1
           Write (98, *) ' '
           Write (98, *) ihpr2
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=1, 20)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=21, 40)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=41, 60)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=61, 80)
           Write (98, *) ' '
           Write (98, *)(hint1(i), i=81, 100)
           Write (98, *) ' '
           Write (98, *) ihnt2
        End If

        !     evaluate Npart (from primary NN collisions) for both proj and targ:
        Call getnp
        !     switch for final parton fragmentation:
        If (ihpr2(20)==0) Goto 2000
        !     In the unlikely case of no interaction (even after loop of 20 in HIJING),
        !     still repeat the event to get an interaction
        !     (this may have an additional "trigger" effect):
        If (natt==0) Then
           imiss = imiss + 1
           If (imiss<=20) Then
              Write (6, *) 'repeated event: natt=0,j,imiss=', j, imiss
              Goto 100
           Else
              Write (6, *) 'missed event: natt=0,j=', j
              Goto 2000
           End If
        End If
        !.....ART initialization and run
        Call arini
        Call arini2(k)
     End Do
     !
     Call artan1
     !lin-9/2012 Analysis is not used:
     !          CALL HJANA3
     Call artmn
     !lin-9/2012 Analysis is not used:
     !          CALL HJANA4
     Call artan2
2000 End Do
  !
  Call artout(nevnt)
  !lin-5/2009 ctest off:
  !       call flowh0(NEVNT,2)
  !       call flowp(2)
  !       call iniflw(NEVNT,2)
  !       call frztm(NEVNT,2)
  !
  Stop
111 Format (A8)
End Program ampt
