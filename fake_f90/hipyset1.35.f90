!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!.................... hipyset1.35.f
!
!
!
!     Modified for HIJING program
!
!    modification July 22, 1997  In pyremnn put an upper limit
!     on the total pt kick the parton can accumulate via multiple
!     scattering. Set the upper limit to be the sqrt(s)/2,
!     this is fix cronin bug for Pb+Pb events at SPS energy.
!
!
! Last modification Oct. 1993 to comply with non-vax
! machines' compiler
!
!*********************************************************************

Subroutine lu2ent(ip, kf1, kf2, pecm)

!...Purpose: to store two partons/particles in their CM frame,
!...with the first along the +z axis.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/

!...Standard checks.
  mstu(28) = 0
  If (mstu(12)>=1) Call lulist(0)
  ipa = max(1, iabs(ip))
  If (ipa>mstu(4)-1) Call luerrm(21, '(LU2ENT:) writing outside LUJETS memory')
  kc1 = lucomp(kf1)
  kc2 = lucomp(kf2)
  If (kc1==0 .Or. kc2==0) Call luerrm(12, '(LU2ENT:) unknown flavour code')

!...Find masses. Reset K, P and V vectors.
  pm1 = 0.
  If (mstu(10)==1) pm1 = p(ipa, 5)
  If (mstu(10)>=2) pm1 = ulmass(kf1)
  pm2 = 0.
  If (mstu(10)==1) pm2 = p(ipa+1, 5)
  If (mstu(10)>=2) pm2 = ulmass(kf2)
  Do i = ipa, ipa + 1
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
  End Do

!...Check flavours.
  kq1 = kchg(kc1, 2)*isign(1, kf1)
  kq2 = kchg(kc2, 2)*isign(1, kf2)
  If (kq1+kq2/=0 .And. kq1+kq2/=4) Call luerrm(2, '(LU2ENT:) unphysical flavour combination')
  k(ipa, 2) = kf1
  k(ipa+1, 2) = kf2

!...Store partons/particles in K vectors for normal case.
  If (ip>=0) Then
    k(ipa, 1) = 1
    If (kq1/=0 .And. kq2/=0) k(ipa, 1) = 2
    k(ipa+1, 1) = 1

!...Store partons in K vectors for parton shower evolution.
  Else
    If (kq1==0 .Or. kq2==0) Call luerrm(2, '(LU2ENT:) requested flavours can not develop parton shower')
    k(ipa, 1) = 3
    k(ipa+1, 1) = 3
    k(ipa, 4) = mstu(5)*(ipa+1)
    k(ipa, 5) = k(ipa, 4)
    k(ipa+1, 4) = mstu(5)*ipa
    k(ipa+1, 5) = k(ipa+1, 4)
  End If

!...Check kinematics and store partons/particles in P vectors.
  If (pecm<=pm1+pm2) Call luerrm(13, '(LU2ENT:) energy smaller than sum of masses')
  pa = sqrt(max(0.,(pecm**2-pm1**2-pm2**2)**2-(2.*pm1*pm2)**2))/(2.*pecm)
  p(ipa, 3) = pa
  p(ipa, 4) = sqrt(pm1**2+pa**2)
  p(ipa, 5) = pm1
  p(ipa+1, 3) = -pa
  p(ipa+1, 4) = sqrt(pm2**2+pa**2)
  p(ipa+1, 5) = pm2

!...Set N. Optionally fragment/decay.
  n = ipa + 1
  If (ip==0) Call luexec

  Return
End Subroutine lu2ent

!*********************************************************************

Subroutine lugive(chin)

!...Purpose: to set values of commonblock variables.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /ludat4/chaf(500)
  Character chaf*8
  Save /ludat4/
  Character chin*(*), chfix*104, chbit*104, chold*8, chnew*8, chnam*4, chvar(17)*4, chalp(2)*26, chind*8, chini*10, chinr*16
  Data chvar/'N', 'K', 'P', 'V', 'MSTU', 'PARU', 'MSTJ', 'PARJ', 'KCHG', 'PMAS', 'PARF', 'VCKM', 'MDCY', 'MDME', 'BRAT', 'KFDP', 'CHAF'/
  Data chalp/'abcdefghijklmnopqrstuvwxyz', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

!...Length of character variable. Subdivide it into instructions.
  If (mstu(12)>=1) Call lulist(0)
  chbit = chin // ' '
  lbit = 101
  100 lbit = lbit - 1
  If (chbit(lbit:lbit)==' ') Goto 100
  ltot = 0
  Do lcom = 1, lbit
    If (chbit(lcom:lcom)==' ') Goto 110
    ltot = ltot + 1
    chfix(ltot:ltot) = chbit(lcom:lcom)
  110 End Do
  llow = 0
  120 lhig = llow + 1
  130 lhig = lhig + 1
  If (lhig<=ltot .And. chfix(lhig:lhig)/=';') Goto 130
  lbit = lhig - llow - 1
  chbit(1:lbit) = chfix(llow+1:lhig-1)

!...Identify commonblock variable.
  lnam = 1
  140 lnam = lnam + 1
  If (chbit(lnam:lnam)/='(' .And. chbit(lnam:lnam)/='=' .And. lnam<=4) Goto 140
  chnam = chbit(1:lnam-1) // ' '
  Do lcom = 1, lnam - 1
    Do lalp = 1, 26
      If (chnam(lcom:lcom)==chalp(1)(lalp:lalp)) chnam(lcom:lcom) = chalp(2)(lalp:lalp)
    End Do
  End Do
  ivar = 0
  Do iv = 1, 17
    If (chnam==chvar(iv)) ivar = iv
  End Do
  If (ivar==0) Then
    Call luerrm(18, '(LUGIVE:) do not recognize variable '//chnam)
    llow = lhig
    If (llow<ltot) Goto 120
    Return
  End If

!...Identify any indices.
  i = 0
  j = 0
  If (chbit(lnam:lnam)=='(') Then
    lind = lnam
    170 lind = lind + 1
    If (chbit(lind:lind)/=')' .And. chbit(lind:lind)/=',') Goto 170
    chind = ' '
    If ((chbit(lnam+1:lnam+1)=='C' .Or. chbit(lnam+1:lnam+1)=='c') .And. (ivar==9 .Or. ivar==10 .Or. ivar==13 .Or. ivar==17)) Then
      chind(lnam-lind+11:8) = chbit(lnam+2:lind-1)
      Read (chind, '(I8)') i1
      i = lucomp(i1)
    Else
      chind(lnam-lind+10:8) = chbit(lnam+1:lind-1)
      Read (chind, '(I8)') i
    End If
    lnam = lind
    If (chbit(lnam:lnam)==')') lnam = lnam + 1
  End If
  If (chbit(lnam:lnam)==',') Then
    lind = lnam
    180 lind = lind + 1
    If (chbit(lind:lind)/=')' .And. chbit(lind:lind)/=',') Goto 180
    chind = ' '
    chind(lnam-lind+10:8) = chbit(lnam+1:lind-1)
    Read (chind, '(I8)') j
    lnam = lind + 1
  End If

!...Check that indices allowed and save old value.
  ierr = 1
  If (chbit(lnam:lnam)/='=') Goto 190
  If (ivar==1) Then
    If (i/=0 .Or. j/=0) Goto 190
    iold = n
  Else If (ivar==2) Then
    If (i<1 .Or. i>mstu(4) .Or. j<1 .Or. j>5) Goto 190
    iold = k(i, j)
  Else If (ivar==3) Then
    If (i<1 .Or. i>mstu(4) .Or. j<1 .Or. j>5) Goto 190
    rold = p(i, j)
  Else If (ivar==4) Then
    If (i<1 .Or. i>mstu(4) .Or. j<1 .Or. j>5) Goto 190
    rold = v(i, j)
  Else If (ivar==5) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    iold = mstu(i)
  Else If (ivar==6) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    rold = paru(i)
  Else If (ivar==7) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    iold = mstj(i)
  Else If (ivar==8) Then
    If (i<1 .Or. i>200 .Or. j/=0) Goto 190
    rold = parj(i)
  Else If (ivar==9) Then
    If (i<1 .Or. i>mstu(6) .Or. j<1 .Or. j>3) Goto 190
    iold = kchg(i, j)
  Else If (ivar==10) Then
    If (i<1 .Or. i>mstu(6) .Or. j<1 .Or. j>4) Goto 190
    rold = pmas(i, j)
  Else If (ivar==11) Then
    If (i<1 .Or. i>2000 .Or. j/=0) Goto 190
    rold = parf(i)
  Else If (ivar==12) Then
    If (i<1 .Or. i>4 .Or. j<1 .Or. j>4) Goto 190
    rold = vckm(i, j)
  Else If (ivar==13) Then
    If (i<1 .Or. i>mstu(6) .Or. j<1 .Or. j>3) Goto 190
    iold = mdcy(i, j)
  Else If (ivar==14) Then
    If (i<1 .Or. i>mstu(7) .Or. j<1 .Or. j>2) Goto 190
    iold = mdme(i, j)
  Else If (ivar==15) Then
    If (i<1 .Or. i>mstu(7) .Or. j/=0) Goto 190
    rold = brat(i)
  Else If (ivar==16) Then
    If (i<1 .Or. i>mstu(7) .Or. j<1 .Or. j>5) Goto 190
    iold = kfdp(i, j)
  Else If (ivar==17) Then
    If (i<1 .Or. i>mstu(6) .Or. j/=0) Goto 190
    chold = chaf(i)
  End If
  ierr = 0
  190 If (ierr==1) Then
    Call luerrm(18, '(LUGIVE:) unallowed indices for '//chbit(1:lnam-1))
    llow = lhig
    If (llow<ltot) Goto 120
    Return
  End If

!...Print current value of variable. Loop back.
  If (lnam>=lbit) Then
    chbit(lnam:14) = ' '
    chbit(15:60) = ' has the value                                '
    If (ivar==1 .Or. ivar==2 .Or. ivar==5 .Or. ivar==7 .Or. ivar==9 .Or. ivar==13 .Or. ivar==14 .Or. ivar==16) Then
      Write (chbit(51:60), '(I10)') iold
    Else If (ivar/=17) Then
      Write (chbit(47:60), '(F14.5)') rold
    Else
      chbit(53:60) = chold
    End If
    If (mstu(13)>=1) Write (mstu(11), 1000) chbit(1:60)
    llow = lhig
    If (llow<ltot) Goto 120
    Return
  End If

!...Read in new variable value.
  If (ivar==1 .Or. ivar==2 .Or. ivar==5 .Or. ivar==7 .Or. ivar==9 .Or. ivar==13 .Or. ivar==14 .Or. ivar==16) Then
    chini = ' '
    chini(lnam-lbit+11:10) = chbit(lnam+1:lbit)
    Read (chini, '(I10)') inew
  Else If (ivar/=17) Then
    chinr = ' '
    chinr(lnam-lbit+17:16) = chbit(lnam+1:lbit)
    Read (chinr, '(F16.2)') rnew
  Else
    chnew = chbit(lnam+1:lbit) // ' '
  End If

!...Store new variable value.
  If (ivar==1) Then
    n = inew
  Else If (ivar==2) Then
    k(i, j) = inew
  Else If (ivar==3) Then
    p(i, j) = rnew
  Else If (ivar==4) Then
    v(i, j) = rnew
  Else If (ivar==5) Then
    mstu(i) = inew
  Else If (ivar==6) Then
    paru(i) = rnew
  Else If (ivar==7) Then
    mstj(i) = inew
  Else If (ivar==8) Then
    parj(i) = rnew
  Else If (ivar==9) Then
    kchg(i, j) = inew
  Else If (ivar==10) Then
    pmas(i, j) = rnew
  Else If (ivar==11) Then
    parf(i) = rnew
  Else If (ivar==12) Then
    vckm(i, j) = rnew
  Else If (ivar==13) Then
    mdcy(i, j) = inew
  Else If (ivar==14) Then
    mdme(i, j) = inew
  Else If (ivar==15) Then
    brat(i) = rnew
  Else If (ivar==16) Then
    kfdp(i, j) = inew
  Else If (ivar==17) Then
    chaf(i) = chnew
  End If

!...Write old and new value. Loop back.
  chbit(lnam:14) = ' '
  chbit(15:60) = ' changed from                to               '
  If (ivar==1 .Or. ivar==2 .Or. ivar==5 .Or. ivar==7 .Or. ivar==9 .Or. ivar==13 .Or. ivar==14 .Or. ivar==16) Then
    Write (chbit(33:42), '(I10)') iold
    Write (chbit(51:60), '(I10)') inew
  Else If (ivar/=17) Then
    Write (chbit(29:42), '(F14.5)') rold
    Write (chbit(47:60), '(F14.5)') rnew
  Else
    chbit(35:42) = chold
    chbit(53:60) = chnew
  End If
  If (mstu(13)>=1) Write (mstu(11), 1000) chbit(1:60)
  llow = lhig
  If (llow<ltot) Goto 120

  Return

!...Format statement for output on unit MSTU(11) (by default 6).
  1000 Format (5X, A60)
End Subroutine lugive

!*********************************************************************

Subroutine luexec

!...Purpose: to administrate the fragmentation and decay chain.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Dimension ps(2, 6)

!...Initialize and reset.
  mstu(24) = 0
  If (mstu(12)>=1) Call lulist(0)
  mstu(31) = mstu(31) + 1
  mstu(1) = 0
  mstu(2) = 0
  mstu(3) = 0
  mcons = 1

!...Sum up momentum, energy and charge for starting entries.
  nsav = n
  Do i = 1, 2
    Do j = 1, 6
      ps(i, j) = 0.
    End Do
  End Do
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 120
    Do j = 1, 4
      ps(1, j) = ps(1, j) + p(i, j)
    End Do
    ps(1, 6) = ps(1, 6) + luchge(k(i,2))
  120 End Do
  paru(21) = ps(1, 4)

!...Prepare system for subsequent fragmentation/decay.
  Call luprep(0)

!...Loop through jet fragmentation and particle decays.
  mbe = 0
  130 mbe = mbe + 1
  ip = 0
  140 ip = ip + 1
  kc = 0
  If (k(ip,1)>0 .And. k(ip,1)<=10) kc = lucomp(k(ip,2))
  If (kc==0) Then

!...Particle decay if unstable and allowed. Save long-lived particle
!...decays until second pass after Bose-Einstein effects.
  Else If (kchg(kc,2)==0) Then
!lin-4/2008 break up compound IF statements:
!        IF(MSTJ(21).GE.1.AND.MDCY(KC,1).GE.1.AND.(MSTJ(51).LE.0.OR.MBE.
!     &  EQ.2.OR.PMAS(KC,2).GE.PARJ(91).OR.IABS(K(IP,2)).EQ.311))
!     &  CALL LUDECY(IP)
    If (mstj(21)>=1 .And. mdcy(kc,1)>=1) Then
      If (mstj(51)<=0 .Or. mbe==2 .Or. pmas(kc,2)>=parj(91) .Or. iabs(k(ip,2))==311) Call ludecy(ip)
    End If
!
!...Decay products may develop a shower.
    If (mstj(92)>0) Then
      ip1 = mstj(92)
      qmax = sqrt(max(0.,(p(ip1,4)+p(ip1+1,4))**2-(p(ip1,1)+p(ip1+1,1))**2-(p(ip1,2)+p(ip1+1,2))**2-(p(ip1,3)+p(ip1+1,3))**2))
      Call lushow(ip1, ip1+1, qmax)
      Call luprep(ip1)
      mstj(92) = 0
    Else If (mstj(92)<0) Then
      ip1 = -mstj(92)
!lin-8/19/02 avoid actual argument in common blocks of LUSHOW:
!          CALL LUSHOW(IP1,-3,P(IP,5))
      pip5 = p(ip, 5)
      Call lushow(ip1, -3, pip5)
      Call luprep(ip1)
      mstj(92) = 0
    End If

!...Jet fragmentation: string or independent fragmentation.
  Else If (k(ip,1)==1 .Or. k(ip,1)==2) Then
    mfrag = mstj(1)
    If (mfrag>=1 .And. k(ip,1)==1) mfrag = 2
    If (mstj(21)>=2 .And. k(ip,1)==2 .And. n>ip) Then
      If (k(ip+1,1)==1 .And. k(ip+1,3)==k(ip,3) .And. k(ip,3)>0 .And. k(ip,3)<ip) Then
        If (kchg(lucomp(k(k(ip,3),2)),2)==0) mfrag = min(1, mfrag)
      End If
    End If
    If (mfrag==1) Then
      Call lustrf(ip)
    End If
    If (mfrag==2) Call luindf(ip)
    If (mfrag==2 .And. k(ip,1)==1) mcons = 0
    If (mfrag==2 .And. (mstj(3)<=0 .Or. mod(mstj(3),5)==0)) mcons = 0
  End If

!...Loop back if enough space left in LUJETS and no error abort.
  If (mstu(24)/=0 .And. mstu(21)>=2) Then
  Else If (ip<n .And. n<mstu(4)-20-mstu(32)) Then
    Goto 140
  Else If (ip<n) Then
    Call luerrm(11, '(LUEXEC:) no more memory left in LUJETS')
  End If

!...Include simple Bose-Einstein effect parametrization if desired.
  If (mbe==1 .And. mstj(51)>=1) Then
    Call luboei(nsav)
    Goto 130
  End If

!...Check that momentum, energy and charge were conserved.
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 160
    Do j = 1, 4
      ps(2, j) = ps(2, j) + p(i, j)
    End Do
    ps(2, 6) = ps(2, 6) + luchge(k(i,2))
  160 End Do
  pdev = (abs(ps(2,1)-ps(1,1))+abs(ps(2,2)-ps(1,2))+abs(ps(2,3)-ps(1,3))+abs(ps(2,4)-ps(1,4)))/(1.+abs(ps(2,4))+abs(ps(1,4)))
  If (mcons==1 .And. pdev>paru(11)) Call luerrm(15, '(LUEXEC:) four-momentum was not conserved')
!      IF(MCONS.EQ.1.AND.PDEV.GT.PARU(11)) then
!         CALL LUERRM(15,
!     &'(LUEXEC:) four-momentum was not conserved')
!         write(6,*) 'PS1,2=',PS(1,1),PS(1,2),PS(1,3),PS(1,4),
!     1        '*',PS(2,1),PS(2,2),PS(2,3),PS(2,4)
!      endif

  If (mcons==1 .And. abs(ps(2,6)-ps(1,6))>0.1) Call luerrm(15, '(LUEXEC:) charge was not conserved')

  Return
End Subroutine luexec

!*********************************************************************

Subroutine luprep(ip)

!...Purpose: to rearrange partons along strings, to allow small systems
!...to collapse into one or two particles and to check flavours.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Dimension dps(5), dpc(5), ue(3)

!...Rearrange parton shower product listing along strings: begin loop.
  i1 = n
  Do mqgst = 1, 2
    Do i = max(1, ip), n
      If (k(i,1)/=3) Goto 120
      kc = lucomp(k(i,2))
      If (kc==0) Goto 120
      kq = kchg(kc, 2)
      If (kq==0 .Or. (mqgst==1 .And. kq==2)) Goto 120

!...Pick up loose string end.
      kcs = 4
      If (kq*isign(1,k(i,2))<0) kcs = 5
      ia = i
      nstp = 0
      100 nstp = nstp + 1
      If (nstp>4*n) Then
        Call luerrm(14, '(LUPREP:) caught in infinite loop')
        Return
      End If

!...Copy undecayed parton.
      If (k(ia,1)==3) Then
        If (i1>=mstu(4)-mstu(32)-5) Then
          Call luerrm(11, '(LUPREP:) no more memory left in LUJETS')
          Return
        End If
        i1 = i1 + 1
        k(i1, 1) = 2
        If (nstp>=2 .And. iabs(k(ia,2))/=21) k(i1, 1) = 1
        k(i1, 2) = k(ia, 2)
        k(i1, 3) = ia
        k(i1, 4) = 0
        k(i1, 5) = 0
        Do j = 1, 5
          p(i1, j) = p(ia, j)
          v(i1, j) = v(ia, j)
        End Do
        k(ia, 1) = k(ia, 1) + 10
        If (k(i1,1)==1) Goto 120
      End If

!...Go to next parton in colour space.
      ib = ia
      If (mod(k(ib,kcs)/mstu(5)**2,2)==0 .And. mod(k(ib,kcs),mstu(5))/=0) Then
        ia = mod(k(ib,kcs), mstu(5))
        k(ib, kcs) = k(ib, kcs) + mstu(5)**2
        mrev = 0
      Else
        If (k(ib,kcs)>=2*mstu(5)**2 .Or. mod(k(ib,kcs)/mstu(5),mstu(5))==0) kcs = 9 - kcs
        ia = mod(k(ib,kcs)/mstu(5), mstu(5))
        k(ib, kcs) = k(ib, kcs) + 2*mstu(5)**2
        mrev = 1
      End If
      If (ia<=0 .Or. ia>n) Then
        Call luerrm(12, '(LUPREP:) colour rearrangement failed')
        Return
      End If
      If (mod(k(ia,4)/mstu(5),mstu(5))==ib .Or. mod(k(ia,5)/mstu(5),mstu(5))==ib) Then
        If (mrev==1) kcs = 9 - kcs
        If (mod(k(ia,kcs)/mstu(5),mstu(5))/=ib) kcs = 9 - kcs
        k(ia, kcs) = k(ia, kcs) + 2*mstu(5)**2
      Else
        If (mrev==0) kcs = 9 - kcs
        If (mod(k(ia,kcs),mstu(5))/=ib) kcs = 9 - kcs
        k(ia, kcs) = k(ia, kcs) + mstu(5)**2
      End If
      If (ia/=i) Goto 100
      k(i1, 1) = 1
    120 End Do
  End Do
  n = i1

!...Find lowest-mass colour singlet jet system, OK if above thresh.
  If (mstj(14)<=0) Goto 320
  ns = n
  140 nsin = n - ns
  pdm = 1. + parj(32)
  ic = 0
  Do i = max(1, ip), ns
    If (k(i,1)/=1 .And. k(i,1)/=2) Then
    Else If (k(i,1)==2 .And. ic==0) Then
      nsin = nsin + 1
      ic = i
      Do j = 1, 4
        dps(j) = dble(p(i,j))
      End Do
      mstj(93) = 1
      dps(5) = dble(ulmass(k(i,2)))
    Else If (k(i,1)==2) Then
      Do j = 1, 4
        dps(j) = dps(j) + dble(p(i,j))
      End Do
    Else If (ic/=0 .And. kchg(lucomp(k(i,2)),2)/=0) Then
      Do j = 1, 4
        dps(j) = dps(j) + dble(p(i,j))
      End Do
      mstj(93) = 1
      dps(5) = dps(5) + dble(ulmass(k(i,2)))
      pd = sngl(sqrt(max(0D0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2))-dps(5))
      If (pd<pdm) Then
        pdm = pd
        Do j = 1, 5
          dpc(j) = dps(j)
        End Do
        ic1 = ic
        ic2 = i
      End If
      ic = 0
    Else
      nsin = nsin + 1
    End If
  End Do
  If (pdm>=parj(32)) Goto 320

!...Fill small-mass system as cluster.
  nsav = n
  pecm = sngl(sqrt(max(0D0,dpc(4)**2-dpc(1)**2-dpc(2)**2-dpc(3)**2)))
  k(n+1, 1) = 11
  k(n+1, 2) = 91
  k(n+1, 3) = ic1
  k(n+1, 4) = n + 2
  k(n+1, 5) = n + 3
  p(n+1, 1) = sngl(dpc(1))
  p(n+1, 2) = sngl(dpc(2))
  p(n+1, 3) = sngl(dpc(3))
  p(n+1, 4) = sngl(dpc(4))
  p(n+1, 5) = pecm

!...Form two particles from flavours of lowest-mass system, if feasible.
  k(n+2, 1) = 1
  k(n+3, 1) = 1
  If (mstu(16)/=2) Then
    k(n+2, 3) = n + 1
    k(n+3, 3) = n + 1
  Else
    k(n+2, 3) = ic1
    k(n+3, 3) = ic2
  End If
  k(n+2, 4) = 0
  k(n+3, 4) = 0
  k(n+2, 5) = 0
  k(n+3, 5) = 0
  If (iabs(k(ic1,2))/=21) Then
    kc1 = lucomp(k(ic1,2))
    kc2 = lucomp(k(ic2,2))
    If (kc1==0 .Or. kc2==0) Goto 320
    kq1 = kchg(kc1, 2)*isign(1, k(ic1,2))
    kq2 = kchg(kc2, 2)*isign(1, k(ic2,2))
    If (kq1+kq2/=0) Goto 320
    200 Call lukfdi(k(ic1,2), 0, kfln, k(n+2,2))
    Call lukfdi(k(ic2,2), -kfln, kfldmp, k(n+3,2))
    If (k(n+2,2)==0 .Or. k(n+3,2)==0) Goto 200
  Else
    If (iabs(k(ic2,2))/=21) Goto 320
    210 Call lukfdi(1+int((2.+parj(2))*rlu(0)), 0, kfln, kfdmp)
    Call lukfdi(kfln, 0, kflm, k(n+2,2))
    Call lukfdi(-kfln, -kflm, kfldmp, k(n+3,2))
    If (k(n+2,2)==0 .Or. k(n+3,2)==0) Goto 210
  End If
  p(n+2, 5) = ulmass(k(n+2,2))
  p(n+3, 5) = ulmass(k(n+3,2))
  If (p(n+2,5)+p(n+3,5)+parj(64)>=pecm .And. nsin==1) Goto 320
  If (p(n+2,5)+p(n+3,5)+parj(64)>=pecm) Goto 260

!...Perform two-particle decay of jet system, if possible.
!lin-5/2012:
!      IF(PECM.GE.0.02d0*DPC(4)) THEN
  If (dble(pecm)>=0.02D0*dpc(4)) Then
    pa = sqrt((pecm**2-(p(n+2,5)+p(n+3,5))**2)*(pecm**2-(p(n+2,5)-p(n+3,5))**2))/(2.*pecm)
    ue(3) = 2.*rlu(0) - 1.
    phi = paru(2)*rlu(0)
    ue(1) = sqrt(1.-ue(3)**2)*cos(phi)
    ue(2) = sqrt(1.-ue(3)**2)*sin(phi)
    Do j = 1, 3
      p(n+2, j) = pa*ue(j)
      p(n+3, j) = -pa*ue(j)
    End Do
    p(n+2, 4) = sqrt(pa**2+p(n+2,5)**2)
    p(n+3, 4) = sqrt(pa**2+p(n+3,5)**2)
    Call ludbrb(n+2, n+3, 0., 0., dpc(1)/dpc(4), dpc(2)/dpc(4), dpc(3)/dpc(4))
  Else
    np = 0
    Do i = ic1, ic2
      If (k(i,1)==1 .Or. k(i,1)==2) np = np + 1
    End Do
    ha = p(ic1, 4)*p(ic2, 4) - p(ic1, 1)*p(ic2, 1) - p(ic1, 2)*p(ic2, 2) - p(ic1, 3)*p(ic2, 3)
    If (np>=3 .Or. ha<=1.25*p(ic1,5)*p(ic2,5)) Goto 260
    hd1 = 0.5*(p(n+2,5)**2-p(ic1,5)**2)
    hd2 = 0.5*(p(n+3,5)**2-p(ic2,5)**2)
    hr = sqrt(max(0.,((ha-hd1-hd2)**2-(p(n+2,5)*p(n+3,5))**2)/(ha**2-(p(ic1,5)*p(ic2,5))**2))) - 1.
    hc = p(ic1, 5)**2 + 2.*ha + p(ic2, 5)**2
    hk1 = ((p(ic2,5)**2+ha)*hr+hd1-hd2)/hc
    hk2 = ((p(ic1,5)**2+ha)*hr+hd2-hd1)/hc
    Do j = 1, 4
      p(n+2, j) = (1.+hk1)*p(ic1, j) - hk2*p(ic2, j)
      p(n+3, j) = (1.+hk2)*p(ic2, j) - hk1*p(ic1, j)
    End Do
  End If
  Do j = 1, 4
    v(n+1, j) = v(ic1, j)
    v(n+2, j) = v(ic1, j)
    v(n+3, j) = v(ic2, j)
  End Do
  v(n+1, 5) = 0.
  v(n+2, 5) = 0.
  v(n+3, 5) = 0.
  n = n + 3
  Goto 300

!...Else form one particle from the flavours available, if possible.
  260 k(n+1, 5) = n + 2
  If (iabs(k(ic1,2))>100 .And. iabs(k(ic2,2))>100) Then
    Goto 320
  Else If (iabs(k(ic1,2))/=21) Then
    Call lukfdi(k(ic1,2), k(ic2,2), kfldmp, k(n+2,2))
  Else
    kfln = 1 + int((2.+parj(2))*rlu(0))
    Call lukfdi(kfln, -kfln, kfldmp, k(n+2,2))
  End If
  If (k(n+2,2)==0) Goto 260
  p(n+2, 5) = ulmass(k(n+2,2))

!...Find parton/particle which combines to largest extra mass.
  ir = 0
  ha = 0.
  Do mcomb = 1, 3
    If (ir/=0) Goto 280
    Do i = max(1, ip), n
      If (k(i,1)<=0 .Or. k(i,1)>10 .Or. (i>=ic1 .And. i<=ic2 .And. k(i,1)>=1 .And. k(i,1)<=2)) Goto 270
      If (mcomb==1) kci = lucomp(k(i,2))
      If (mcomb==1 .And. kci==0) Goto 270
      If (mcomb==1 .And. kchg(kci,2)==0 .And. i<=ns) Goto 270
      If (mcomb==2 .And. iabs(k(i,2))>10 .And. iabs(k(i,2))<=100) Goto 270
      hcr = sngl(dpc(4))*p(i, 4) - sngl(dpc(1))*p(i, 1) - sngl(dpc(2))*p(i, 2) - sngl(dpc(3))*p(i, 3)
      If (hcr>ha) Then
        ir = i
        ha = hcr
      End If
    270 End Do
  280 End Do

!...Shuffle energy and momentum to put new particle on mass shell.
  hb = pecm**2 + ha
  hc = p(n+2, 5)**2 + ha
  hd = p(ir, 5)**2 + ha
!******************CHANGES BY HIJING************
  hk2 = 0.0
  If (ha**2-(pecm*p(ir,5))**2==0.0 .Or. hb+hd==0.0) Goto 285
!******************
  hk2 = 0.5*(hb*sqrt(((hb+hc)**2-4.*(hb+hd)*p(n+2,5)**2)/(ha**2-(pecm*p(ir,5))**2))-(hb+hc))/(hb+hd)
  285 hk1 = (0.5*(p(n+2,5)**2-pecm**2)+hd*hk2)/hb
  Do j = 1, 4
    p(n+2, j) = (1.+hk1)*sngl(dpc(j)) - hk2*p(ir, j)
    p(ir, j) = (1.+hk2)*p(ir, j) - hk1*sngl(dpc(j))
    v(n+1, j) = v(ic1, j)
    v(n+2, j) = v(ic1, j)
  End Do
  v(n+1, 5) = 0.
  v(n+2, 5) = 0.
  n = n + 2

!...Mark collapsed system and store daughter pointers. Iterate.
  300 Do i = ic1, ic2
    If ((k(i,1)==1 .Or. k(i,1)==2) .And. kchg(lucomp(k(i,2)),2)/=0) Then
      k(i, 1) = k(i, 1) + 10
      If (mstu(16)/=2) Then
        k(i, 4) = nsav + 1
        k(i, 5) = nsav + 1
      Else
        k(i, 4) = nsav + 2
        k(i, 5) = n
      End If
    End If
  End Do
  If (n<mstu(4)-mstu(32)-5) Goto 140

!...Check flavours and invariant masses in parton systems.
  320 np = 0
  kfn = 0
  kqs = 0
  Do j = 1, 5
    dps(j) = 0D0
  End Do
  Do i = max(1, ip), n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 360
    kc = lucomp(k(i,2))
    If (kc==0) Goto 360
    kq = kchg(kc, 2)*isign(1, k(i,2))
    If (kq==0) Goto 360
    np = np + 1
    If (kq/=2) Then
      kfn = kfn + 1
      kqs = kqs + kq
      mstj(93) = 1
      dps(5) = dps(5) + dble(ulmass(k(i,2)))
    End If
    Do j = 1, 4
      dps(j) = dps(j) + dble(p(i,j))
    End Do

!lin-4/12/01:
!     np: # of partons, KFN: number of quarks and diquarks,
!     KC=0 for color singlet system, -1 for quarks and anti-diquarks,
!     1 for quarks and anti-diquarks, and 2 for gluons:
    If (k(i,1)==1) Then
!lin-4/12/01     end of color singlet system.
      If (np/=1 .And. (kfn==1 .Or. kfn>=3 .Or. kqs/=0)) Call luerrm(2, '(LUPREP:) unphysical flavour combination')

!lin-4/16/01: 'jet system' should be defined as np.ne.2:
!        IF(NP.NE.1.AND.DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2.LT.
!     &  (0.9*PARJ(32)+DPS(5))**2) CALL LUERRM(3,
!     &  '(LUPREP:) too small mass in jet system')
      If (np/=2 .And. dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2<(0.9D0*dble(parj(32))+dps(5))**2) Then
        Call luerrm(3, '(LUPREP:) too small mass in jet system')
        Write (6, *) 'DPS(1-5),KI1-5=', dps(1), dps(2), dps(3), dps(4), dps(5), '*', k(i, 1), k(i, 2), k(i, 3), k(i, 4), k(i, 5)
      End If

      np = 0
      kfn = 0
      kqs = 0
      Do j = 1, 5
        dps(j) = 0D0
      End Do
    End If
  360 End Do

  Return
End Subroutine luprep

!*********************************************************************

Subroutine lustrf(ip)
!...Purpose: to handle the fragmentation of an arbitrary colour singlet
!...jet system according to the Lund string fragmentation model.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension dps(5), kfl(3), pmq(3), px(3), py(3), gam(3), ie(2), pr(2), in(9), dhm(4), dhg(4), dp(5, 5), irank(2), mju(4), iju(3), pju(5, 5), tju(5), kfjh(2), njs(2), kfjs(2), pjs(4, 5)

!...Function: four-product of two vectors.
  four(i, j) = p(i, 4)*p(j, 4) - p(i, 1)*p(j, 1) - p(i, 2)*p(j, 2) - p(i, 3)*p(j, 3)
  dfour(i, j) = dp(i, 4)*dp(j, 4) - dp(i, 1)*dp(j, 1) - dp(i, 2)*dp(j, 2) - dp(i, 3)*dp(j, 3)

!...Reset counters. Identify parton system.
  mstj(91) = 0
  nsav = n
  np = 0
  kqsum = 0
  Do j = 1, 5
    dps(j) = 0D0
  End Do
  mju(1) = 0
  mju(2) = 0
  i = ip - 1
  110 i = i + 1
  If (i>min(n,mstu(4)-mstu(32))) Then
    Call luerrm(12, '(LUSTRF:) failed to reconstruct jet system')
    If (mstu(21)>=1) Return
  End If
  If (k(i,1)/=1 .And. k(i,1)/=2 .And. k(i,1)/=41) Goto 110
  kc = lucomp(k(i,2))
  If (kc==0) Goto 110
  kq = kchg(kc, 2)*isign(1, k(i,2))
  If (kq==0) Goto 110
  If (n+5*np+11>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSTRF:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If

!...Take copy of partons to be considered. Check flavour sum.
  np = np + 1
  Do j = 1, 5
    k(n+np, j) = k(i, j)
    p(n+np, j) = p(i, j)
    dps(j) = dps(j) + dble(p(i,j))
  End Do
  k(n+np, 3) = i
  If (p(n+np,4)**2<p(n+np,1)**2+p(n+np,2)**2+p(n+np,3)**2) Then
    p(n+np, 4) = sqrt(p(n+np,1)**2+p(n+np,2)**2+p(n+np,3)**2+p(n+np,5)**2)
    dps(4) = dps(4) + dble(max(0.,p(n+np,4)-p(i,4)))
  End If
  If (kq/=2) kqsum = kqsum + kq
  If (k(i,1)==41) Then
    kqsum = kqsum + 2*kq
    If (kqsum==kq) mju(1) = n + np
    If (kqsum/=kq) mju(2) = n + np
  End If
  If (k(i,1)==2 .Or. k(i,1)==41) Goto 110
  If (kqsum/=0) Then
    Call luerrm(12, '(LUSTRF:) unphysical flavour combination')
    If (mstu(21)>=1) Return
  End If

!...Boost copied system to CM frame (for better numerical precision).
  Call ludbrb(n+1, n+np, 0., 0., -dps(1)/dps(4), -dps(2)/dps(4), -dps(3)/dps(4))

!...Search for very nearby partons that may be recombined.
  ntryr = 0
  paru12 = paru(12)
  paru13 = paru(13)
  mju(3) = mju(1)
  mju(4) = mju(2)
  nr = np
  130 If (nr>=3) Then
    pdrmin = 2.*paru12
    Do i = n + 1, n + nr
      If (i==n+nr .And. iabs(k(n+1,2))/=21) Goto 140
      i1 = i + 1
      If (i==n+nr) i1 = n + 1
      If (k(i,1)==41 .Or. k(i1,1)==41) Goto 140
      If (mju(1)/=0 .And. i1<mju(1) .And. iabs(k(i1,2))/=21) Goto 140
      If (mju(2)/=0 .And. i>mju(2) .And. iabs(k(i,2))/=21) Goto 140
      pap = sqrt((p(i,1)**2+p(i,2)**2+p(i,3)**2)*(p(i1,1)**2+p(i1,2)**2+p(i1,3)**2))
      pvp = p(i, 1)*p(i1, 1) + p(i, 2)*p(i1, 2) + p(i, 3)*p(i1, 3)
      pdr = 4.*(pap-pvp)**2/(paru13**2*pap+2.*(pap-pvp))
      If (pdr<pdrmin) Then
        ir = i
        pdrmin = pdr
      End If
    140 End Do

!...Recombine very nearby partons to avoid machine precision problems.
    If (pdrmin<paru12 .And. ir==n+nr) Then
      Do j = 1, 4
        p(n+1, j) = p(n+1, j) + p(n+nr, j)
      End Do
      p(n+1, 5) = sqrt(max(0.,p(n+1,4)**2-p(n+1,1)**2-p(n+1,2)**2-p(n+1,3)**2))
      nr = nr - 1
      Goto 130
    Else If (pdrmin<paru12) Then
      Do j = 1, 4
        p(ir, j) = p(ir, j) + p(ir+1, j)
      End Do
      p(ir, 5) = sqrt(max(0.,p(ir,4)**2-p(ir,1)**2-p(ir,2)**2-p(ir,3)**2))
      Do i = ir + 1, n + nr - 1
        k(i, 2) = k(i+1, 2)
        Do j = 1, 5
          p(i, j) = p(i+1, j)
        End Do
      End Do
      If (ir==n+nr-1) k(ir, 2) = k(n+nr, 2)
      nr = nr - 1
      If (mju(1)>ir) mju(1) = mju(1) - 1
      If (mju(2)>ir) mju(2) = mju(2) - 1
      Goto 130
    End If
  End If
  ntryr = ntryr + 1

!...Reset particle counter. Skip ahead if no junctions are present;
!...this is usually the case!
  nrs = max(5*nr+11, np)
  ntry = 0
  180 ntry = ntry + 1
  If (ntry>100 .And. ntryr<=4) Then
    paru12 = 4.*paru12
    paru13 = 2.*paru13
    Goto 130
  Else If (ntry>100) Then
    Call luerrm(14, '(LUSTRF:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  i = n + nrs
  If (mju(1)==0 .And. mju(2)==0) Goto 500
  Do jt = 1, 2
    njs(jt) = 0
    If (mju(jt)==0) Goto 490
    js = 3 - 2*jt

!...Find and sum up momentum on three sides of junction. Check flavours.
    Do iu = 1, 3
      iju(iu) = 0
      Do j = 1, 5
        pju(iu, j) = 0.
      End Do
    End Do
    iu = 0
    Do i1 = n + 1 + (jt-1)*(nr-1), n + nr + (jt-1)*(1-nr), js
      If (k(i1,2)/=21 .And. iu<=2) Then
        iu = iu + 1
        iju(iu) = i1
      End If
      Do j = 1, 4
        pju(iu, j) = pju(iu, j) + p(i1, j)
      End Do
    End Do
    Do iu = 1, 3
      pju(iu, 5) = sqrt(pju(iu,1)**2+pju(iu,2)**2+pju(iu,3)**2)
    End Do
    If (k(iju(3),2)/100/=10*k(iju(1),2)+k(iju(2),2) .And. k(iju(3),2)/100/=10*k(iju(2),2)+k(iju(1),2)) Then
      Call luerrm(12, '(LUSTRF:) unphysical flavour combination')
      If (mstu(21)>=1) Return
    End If

!...Calculate (approximate) boost to rest frame of junction.
    t12 = (pju(1,1)*pju(2,1)+pju(1,2)*pju(2,2)+pju(1,3)*pju(2,3))/(pju(1,5)*pju(2,5))
    t13 = (pju(1,1)*pju(3,1)+pju(1,2)*pju(3,2)+pju(1,3)*pju(3,3))/(pju(1,5)*pju(3,5))
    t23 = (pju(2,1)*pju(3,1)+pju(2,2)*pju(3,2)+pju(2,3)*pju(3,3))/(pju(2,5)*pju(3,5))
    t11 = sqrt((2./3.)*(1.-t12)*(1.-t13)/(1.-t23))
    t22 = sqrt((2./3.)*(1.-t12)*(1.-t23)/(1.-t13))
    tsq = sqrt((2.*t11*t22+t12-1.)*(1.+t12))
    t1f = (tsq-t22*(1.+t12))/(1.-t12**2)
    t2f = (tsq-t11*(1.+t12))/(1.-t12**2)
    Do j = 1, 3
      tju(j) = -(t1f*pju(1,j)/pju(1,5)+t2f*pju(2,j)/pju(2,5))
    End Do
    tju(4) = sqrt(1.+tju(1)**2+tju(2)**2+tju(3)**2)
    Do iu = 1, 3
      pju(iu, 5) = tju(4)*pju(iu, 4) - tju(1)*pju(iu, 1) - tju(2)*pju(iu, 2) - tju(3)*pju(iu, 3)
    End Do

!...Put junction at rest if motion could give inconsistencies.
    If (pju(1,5)+pju(2,5)>pju(1,4)+pju(2,4)) Then
      Do j = 1, 3
        tju(j) = 0.
      End Do
      tju(4) = 1.
      pju(1, 5) = pju(1, 4)
      pju(2, 5) = pju(2, 4)
      pju(3, 5) = pju(3, 4)
    End If

!...Start preparing for fragmentation of two strings from junction.
    ista = i
    Do iu = 1, 2
      ns = iju(iu+1) - iju(iu)

!...Junction strings: find longitudinal string directions.
      Do is = 1, ns
        is1 = iju(iu) + is - 1
        is2 = iju(iu) + is
        Do j = 1, 5
          dp(1, j) = dble(0.5*p(is1,j))
          If (is==1) dp(1, j) = dble(p(is1,j))
          dp(2, j) = dble(0.5*p(is2,j))
          If (is==ns) dp(2, j) = -dble(pju(iu,j))
        End Do
        If (is==ns) dp(2, 4) = dble(sqrt(pju(iu,1)**2+pju(iu,2)**2+pju(iu,3)**2))
        If (is==ns) dp(2, 5) = 0D0
        dp(3, 5) = dfour(1, 1)
        dp(4, 5) = dfour(2, 2)
        dhkc = dfour(1, 2)
        If (dp(3,5)+2D0*dhkc+dp(4,5)<=0D0) Then
          dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
          dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
          dp(3, 5) = 0D0
          dp(4, 5) = 0D0
          dhkc = dfour(1, 2)
        End If
        dhks = sqrt(dhkc**2-dp(3,5)*dp(4,5))
        dhk1 = 0.5D0*((dp(4,5)+dhkc)/dhks-1D0)
        dhk2 = 0.5D0*((dp(3,5)+dhkc)/dhks-1D0)
        in1 = n + nr + 4*is - 3
        p(in1, 5) = sngl(sqrt(dp(3,5)+2D0*dhkc+dp(4,5)))
        Do j = 1, 4
          p(in1, j) = sngl((1D0+dhk1)*dp(1,j)-dhk2*dp(2,j))
          p(in1+1, j) = sngl((1D0+dhk2)*dp(2,j)-dhk1*dp(1,j))
        End Do
      End Do

!...Junction strings: initialize flavour, momentum and starting pos.
      isav = i
      270 ntry = ntry + 1
      If (ntry>100 .And. ntryr<=4) Then
        paru12 = 4.*paru12
        paru13 = 2.*paru13
        Goto 130
      Else If (ntry>100) Then
        Call luerrm(14, '(LUSTRF:) caught in infinite loop')
        If (mstu(21)>=1) Return
      End If
      i = isav
      irankj = 0
      ie(1) = k(n+1+(jt/2)*(np-1), 3)
      in(4) = n + nr + 1
      in(5) = in(4) + 1
      in(6) = n + nr + 4*ns + 1
      Do jq = 1, 2
        Do in1 = n + nr + 2 + jq, n + nr + 4*ns - 2 + jq, 4
          p(in1, 1) = 2 - jq
          p(in1, 2) = jq - 1
          p(in1, 3) = 1.
        End Do
      End Do
      kfl(1) = k(iju(iu), 2)
      px(1) = 0.
      py(1) = 0.
      gam(1) = 0.
      Do j = 1, 5
        pju(iu+3, j) = 0.
      End Do

!...Junction strings: find initial transverse directions.
      Do j = 1, 4
        dp(1, j) = dble(p(in(4),j))
        dp(2, j) = dble(p(in(4)+1,j))
        dp(3, j) = 0D0
        dp(4, j) = 0D0
      End Do
      dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
      dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
      dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
      dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
      dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
      If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1D0
      If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1D0
      If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1D0
      If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1D0
      dhc12 = dfour(1, 2)
      dhcx1 = dfour(3, 1)/dhc12
      dhcx2 = dfour(3, 2)/dhc12
      dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
      dhcy1 = dfour(4, 1)/dhc12
      dhcy2 = dfour(4, 2)/dhc12
      dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
      dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
      Do j = 1, 4
        dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
        p(in(6), j) = sngl(dp(3,j))
        p(in(6)+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
      End Do

!...Junction strings: produce new particle, origin.
      320 i = i + 1
      If (2*i-nsav>=mstu(4)-mstu(32)-5) Then
        Call luerrm(11, '(LUSTRF:) no more memory left in LUJETS')
        If (mstu(21)>=1) Return
      End If
      irankj = irankj + 1
      k(i, 1) = 1
      k(i, 3) = ie(1)
      k(i, 4) = 0
      k(i, 5) = 0

!...Junction strings: generate flavour, hadron, pT, z and Gamma.
      330 Call lukfdi(kfl(1), 0, kfl(3), k(i,2))
      If (k(i,2)==0) Goto 270
      If (mstj(12)>=3 .And. irankj==1 .And. iabs(kfl(1))<=10 .And. iabs(kfl(3))>10) Then
        If (rlu(0)>parj(19)) Goto 330
      End If
      p(i, 5) = ulmass(k(i,2))
      Call luptdi(kfl(1), px(3), py(3))
      pr(1) = p(i, 5)**2 + (px(1)+px(3))**2 + (py(1)+py(3))**2
      Call luzdis(kfl(1), kfl(3), pr(1), z)
      gam(3) = (1.-z)*(gam(1)+pr(1)/z)
      Do j = 1, 3
        in(j) = in(3+j)
      End Do

!...Junction strings: stepping within or from 'low' string region easy.
      If (in(1)+1==in(2) .And. z*p(in(1)+2,3)*p(in(2)+2,3)*p(in(1),5)**2>=pr(1)) Then
        p(in(1)+2, 4) = z*p(in(1)+2, 3)
        p(in(2)+2, 4) = pr(1)/(p(in(1)+2,4)*p(in(1),5)**2)
        Do j = 1, 4
          p(i, j) = (px(1)+px(3))*p(in(3), j) + (py(1)+py(3))*p(in(3)+1, j)
        End Do
        Goto 420
      Else If (in(1)+1==in(2)) Then
        p(in(2)+2, 4) = p(in(2)+2, 3)
        p(in(2)+2, 1) = 1.
        in(2) = in(2) + 4
        If (in(2)>n+nr+4*ns) Goto 270
        If (four(in(1),in(2))<=1E-2) Then
          p(in(1)+2, 4) = p(in(1)+2, 3)
          p(in(1)+2, 1) = 0.
          in(1) = in(1) + 4
        End If
      End If

!...Junction strings: find new transverse directions.
      360 If (in(1)>n+nr+4*ns .Or. in(2)>n+nr+4*ns .Or. in(1)>in(2)) Goto 270
      If (in(1)/=in(4) .Or. in(2)/=in(5)) Then
        Do j = 1, 4
          dp(1, j) = dble(p(in(1),j))
          dp(2, j) = dble(p(in(2),j))
          dp(3, j) = 0D0
          dp(4, j) = 0D0
        End Do
        dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
        dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
        dhc12 = dfour(1, 2)
!lin-5/2012:
!        IF(DHC12.LE.1E-2) THEN
        If (dhc12<=1D-2) Then
          p(in(1)+2, 4) = p(in(1)+2, 3)
          p(in(1)+2, 1) = 0.
          in(1) = in(1) + 4
          Goto 360
        End If
        in(3) = n + nr + 4*ns + 5
        dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
        dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
        dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
        If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1D0
        If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1D0
        If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1D0
        If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1D0
        dhcx1 = dfour(3, 1)/dhc12
        dhcx2 = dfour(3, 2)/dhc12
        dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
        dhcy1 = dfour(4, 1)/dhc12
        dhcy2 = dfour(4, 2)/dhc12
        dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
        dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
        Do j = 1, 4
          dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
          p(in(3), j) = sngl(dp(3,j))
          p(in(3)+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
        End Do
!...Express pT with respect to new axes, if sensible.
        pxp = -(px(3)*four(in(6),in(3))+py(3)*four(in(6)+1,in(3)))
        pyp = -(px(3)*four(in(6),in(3)+1)+py(3)*four(in(6)+1,in(3)+1))
        If (abs(pxp**2+pyp**2-px(3)**2-py(3)**2)<0.01) Then
          px(3) = pxp
          py(3) = pyp
        End If
      End If

!...Junction strings: sum up known four-momentum, coefficients for m2.
      Do j = 1, 4
        dhg(j) = 0D0
        p(i, j) = px(1)*p(in(6), j) + py(1)*p(in(6)+1, j) + px(3)*p(in(3), j) + py(3)*p(in(3)+1, j)
        Do in1 = in(4), in(1) - 4, 4
          p(i, j) = p(i, j) + p(in1+2, 3)*p(in1, j)
        End Do
        Do in2 = in(5), in(2) - 4, 4
          p(i, j) = p(i, j) + p(in2+2, 3)*p(in2, j)
        End Do
      End Do
      dhm(1) = dble(four(i,i))
      dhm(2) = dble(2.*four(i,in(1)))
      dhm(3) = dble(2.*four(i,in(2)))
      dhm(4) = dble(2.*four(in(1),in(2)))

!...Junction strings: find coefficients for Gamma expression.
      Do in2 = in(1) + 1, in(2), 4
        Do in1 = in(1), in2 - 1, 4
          dhc = dble(2.*four(in1,in2))
          dhg(1) = dhg(1) + dble(p(in1+2,1)*p(in2+2,1))*dhc
          If (in1==in(1)) dhg(2) = dhg(2) - dble(p(in2+2,1))*dhc
          If (in2==in(2)) dhg(3) = dhg(3) + dble(p(in1+2,1))*dhc
          If (in1==in(1) .And. in2==in(2)) dhg(4) = dhg(4) - dhc
        End Do
      End Do

!...Junction strings: solve (m2, Gamma) equation system for energies.
      dhs1 = dhm(3)*dhg(4) - dhm(4)*dhg(3)
!lin-5/2012:
!      IF(ABS(DHS1).LT.1E-4) GOTO 270
      If (dabs(dhs1)<1D-4) Goto 270
      dhs2 = dhm(4)*(dble(gam(3))-dhg(1)) - dhm(2)*dhg(3) - dhg(4)*(dble(p(i,5))**2-dhm(1)) + dhg(2)*dhm(3)
      dhs3 = dhm(2)*(dble(gam(3))-dhg(1)) - dhg(2)*(dble(p(i,5))**2-dhm(1))
      p(in(2)+2, 4) = 0.5*sngl(sqrt(max(0D0,dhs2**2-4D0*dhs1*dhs3))/abs(dhs1)-dhs2/dhs1)
      If (dhm(2)+dhm(4)*dble(p(in(2)+2,4))<=0D0) Goto 270
      p(in(1)+2, 4) = (p(i,5)**2-sngl(dhm(1))-sngl(dhm(3))*p(in(2)+2,4))/(sngl(dhm(2))+sngl(dhm(4))*p(in(2)+2,4))

!...Junction strings: step to new region if necessary.
      If (p(in(2)+2,4)>p(in(2)+2,3)) Then
        p(in(2)+2, 4) = p(in(2)+2, 3)
        p(in(2)+2, 1) = 1.
        in(2) = in(2) + 4
        If (in(2)>n+nr+4*ns) Goto 270
        If (four(in(1),in(2))<=1E-2) Then
          p(in(1)+2, 4) = p(in(1)+2, 3)
          p(in(1)+2, 1) = 0.
          in(1) = in(1) + 4
        End If
        Goto 360
      Else If (p(in(1)+2,4)>p(in(1)+2,3)) Then
        p(in(1)+2, 4) = p(in(1)+2, 3)
        p(in(1)+2, 1) = 0.
        in(1) = in(1) + js
        Goto 710
      End If

!...Junction strings: particle four-momentum, remainder, loop back.
      420 Do j = 1, 4
        p(i, j) = p(i, j) + p(in(1)+2, 4)*p(in(1), j) + p(in(2)+2, 4)*p(in(2), j)
        pju(iu+3, j) = pju(iu+3, j) + p(i, j)
      End Do
      If (p(i,4)<=0.) Goto 270
      pju(iu+3, 5) = tju(4)*pju(iu+3, 4) - tju(1)*pju(iu+3, 1) - tju(2)*pju(iu+3, 2) - tju(3)*pju(iu+3, 3)
      If (pju(iu+3,5)<pju(iu,5)) Then
        kfl(1) = -kfl(3)
        px(1) = -px(3)
        py(1) = -py(3)
        gam(1) = gam(3)
        If (in(3)/=in(6)) Then
          Do j = 1, 4
            p(in(6), j) = p(in(3), j)
            p(in(6)+1, j) = p(in(3)+1, j)
          End Do
        End If
        Do jq = 1, 2
          in(3+jq) = in(jq)
          p(in(jq)+2, 3) = p(in(jq)+2, 3) - p(in(jq)+2, 4)
          p(in(jq)+2, 1) = p(in(jq)+2, 1) - (3-2*jq)*p(in(jq)+2, 4)
        End Do
        Goto 320
      End If

!...Junction strings: save quantities left after each string.
      If (iabs(kfl(1))>10) Goto 270
      i = i - 1
      kfjh(iu) = kfl(1)
      Do j = 1, 4
        pju(iu+3, j) = pju(iu+3, j) - p(i+1, j)
      End Do
    End Do

!...Junction strings: put together to new effective string endpoint.
    njs(jt) = i - ista
    kfjs(jt) = k(k(mju(jt+2),3), 2)
    kfls = 2*int(rlu(0)+3.*parj(4)/(1.+3.*parj(4))) + 1
    If (kfjh(1)==kfjh(2)) kfls = 3
    If (ista/=i) kfjs(jt) = isign(1000*max(iabs(kfjh(1)),iabs(kfjh(2)))+100*min(iabs(kfjh(1)),iabs(kfjh(2)))+kfls, kfjh(1))
    Do j = 1, 4
      pjs(jt, j) = pju(1, j) + pju(2, j) + p(mju(jt), j)
      pjs(jt+2, j) = pju(4, j) + pju(5, j)
    End Do
    pjs(jt, 5) = sqrt(max(0.,pjs(jt,4)**2-pjs(jt,1)**2-pjs(jt,2)**2-pjs(jt,3)**2))
  490 End Do

!...Open versus closed strings. Choose breakup region for latter.
  500 If (mju(1)/=0 .And. mju(2)/=0) Then
    ns = mju(2) - mju(1)
    nb = mju(1) - n
  Else If (mju(1)/=0) Then
    ns = n + nr - mju(1)
    nb = mju(1) - n
  Else If (mju(2)/=0) Then
    ns = mju(2) - n
    nb = 1
  Else If (iabs(k(n+1,2))/=21) Then
    ns = nr - 1
    nb = 1
  Else
    ns = nr + 1
    w2sum = 0.
    Do is = 1, nr
      p(n+nr+is, 1) = 0.5*four(n+is, n+is+1-nr*(is/nr))
      w2sum = w2sum + p(n+nr+is, 1)
    End Do
    w2ran = rlu(0)*w2sum
    nb = 0
    520 nb = nb + 1
    w2sum = w2sum - p(n+nr+nb, 1)
    If (w2sum>w2ran .And. nb<nr) Goto 520
  End If

!...Find longitudinal string directions (i.e. lightlike four-vectors).
  Do is = 1, ns
    is1 = n + is + nb - 1 - nr*((is+nb-2)/nr)
    is2 = n + is + nb - nr*((is+nb-1)/nr)
    Do j = 1, 5
      dp(1, j) = dble(p(is1,j))
      If (iabs(k(is1,2))==21) dp(1, j) = 0.5D0*dp(1, j)
      If (is1==mju(1)) dp(1, j) = dble(pjs(1,j)-pjs(3,j))
      dp(2, j) = dble(p(is2,j))
      If (iabs(k(is2,2))==21) dp(2, j) = 0.5D0*dp(2, j)
      If (is2==mju(2)) dp(2, j) = dble(pjs(2,j)-pjs(4,j))
    End Do
    dp(3, 5) = dfour(1, 1)
    dp(4, 5) = dfour(2, 2)
    dhkc = dfour(1, 2)
    If (dp(3,5)+2.D0*dhkc+dp(4,5)<=0.D0) Then
      dp(3, 5) = dp(1, 5)**2
      dp(4, 5) = dp(2, 5)**2
      dp(1, 4) = sqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2+dp(1,5)**2)
      dp(2, 4) = sqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2+dp(2,5)**2)
      dhkc = dfour(1, 2)
    End If
    dhks = sqrt(dhkc**2-dp(3,5)*dp(4,5))
    dhk1 = 0.5D0*((dp(4,5)+dhkc)/dhks-1.D0)
    dhk2 = 0.5D0*((dp(3,5)+dhkc)/dhks-1.D0)
    in1 = n + nr + 4*is - 3
    p(in1, 5) = sqrt(sngl(dp(3,5)+2.D0*dhkc+dp(4,5)))
    Do j = 1, 4
      p(in1, j) = sngl((1.D0+dhk1)*dp(1,j)-dhk2*dp(2,j))
      p(in1+1, j) = sngl((1.D0+dhk2)*dp(2,j)-dhk1*dp(1,j))
    End Do
  End Do

!...Begin initialization: sum up energy, set starting position.
  isav = i
  550 ntry = ntry + 1
  If (ntry>100 .And. ntryr<=4) Then
    paru12 = 4.*paru12
    paru13 = 2.*paru13
    Goto 130
  Else If (ntry>100) Then
    Call luerrm(14, '(LUSTRF:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  i = isav
  Do j = 1, 4
    p(n+nrs, j) = 0.
    Do is = 1, nr
      p(n+nrs, j) = p(n+nrs, j) + p(n+is, j)
    End Do
  End Do
  Do jt = 1, 2
    irank(jt) = 0
    If (mju(jt)/=0) irank(jt) = njs(jt)
    If (ns>nr) irank(jt) = 1
    ie(jt) = k(n+1+(jt/2)*(np-1), 3)
    in(3*jt+1) = n + nr + 1 + 4*(jt/2)*(ns-1)
    in(3*jt+2) = in(3*jt+1) + 1
    in(3*jt+3) = n + nr + 4*ns + 2*jt - 1
    Do in1 = n + nr + 2 + jt, n + nr + 4*ns - 2 + jt, 4
      p(in1, 1) = 2 - jt
      p(in1, 2) = jt - 1
      p(in1, 3) = 1.
    End Do
  End Do

!...Initialize flavour and pT variables for open string.
  If (ns<nr) Then
    px(1) = 0.
    py(1) = 0.
    If (ns==1 .And. mju(1)+mju(2)==0) Call luptdi(0, px(1), py(1))
    px(2) = -px(1)
    py(2) = -py(1)
    Do jt = 1, 2
      kfl(jt) = k(ie(jt), 2)
      If (mju(jt)/=0) kfl(jt) = kfjs(jt)
      mstj(93) = 1
      pmq(jt) = ulmass(kfl(jt))
      gam(jt) = 0.
    End Do

!...Closed string: random initial breakup flavour, pT and vertex.
  Else
    kfl(3) = int(1.+(2.+parj(2))*rlu(0))*(-1)**int(rlu(0)+0.5)
    Call lukfdi(kfl(3), 0, kfl(1), kdump)
    kfl(2) = -kfl(1)
    If (iabs(kfl(1))>10 .And. rlu(0)>0.5) Then
      kfl(2) = -(kfl(1)+isign(10000,kfl(1)))
    Else If (iabs(kfl(1))>10) Then
      kfl(1) = -(kfl(2)+isign(10000,kfl(2)))
    End If
    Call luptdi(kfl(1), px(1), py(1))
    px(2) = -px(1)
    py(2) = -py(1)
    pr3 = min(25., 0.1*p(n+nr+1,5)**2)
    590 Call luzdis(kfl(1), kfl(2), pr3, z)
    zr = pr3/(z*p(n+nr+1,5)**2)
    If (zr>=1.) Goto 590

    Do jt = 1, 2
      mstj(93) = 1
      pmq(jt) = ulmass(kfl(jt))
      gam(jt) = pr3*(1.-z)/z
      in1 = n + nr + 3 + 4*(jt/2)*(ns-1)
      p(in1, jt) = 1. - z
      p(in1, 3-jt) = jt - 1
      p(in1, 3) = (2-jt)*(1.-z) + (jt-1)*z
      p(in1+1, jt) = zr
      p(in1+1, 3-jt) = 2 - jt
      p(in1+1, 3) = (2-jt)*(1.-zr) + (jt-1)*zr
    End Do
  End If

!...Find initial transverse directions (i.e. spacelike four-vectors).
  Do jt = 1, 2
    If (jt==1 .Or. ns==nr-1) Then
      in1 = in(3*jt+1)
      in3 = in(3*jt+3)
      Do j = 1, 4
        dp(1, j) = dble(p(in1,j))
        dp(2, j) = dble(p(in1+1,j))
        dp(3, j) = 0.D0
        dp(4, j) = 0.D0
      End Do
      dp(1, 4) = dsqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
      dp(2, 4) = dsqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
      dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
      dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
      dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
      If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1.D0
      If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1.D0
      If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1.D0
      If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1.D0
      dhc12 = dfour(1, 2)
      dhcx1 = dfour(3, 1)/dhc12
      dhcx2 = dfour(3, 2)/dhc12
      dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
      dhcy1 = dfour(4, 1)/dhc12
      dhcy2 = dfour(4, 2)/dhc12
      dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
      dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
      Do j = 1, 4
        dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
        p(in3, j) = sngl(dp(3,j))
        p(in3+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
      End Do
    Else
      Do j = 1, 4
        p(in3+2, j) = p(in3, j)
        p(in3+3, j) = p(in3+1, j)
      End Do
    End If
  End Do

!...Remove energy used up in junction string fragmentation.
  If (mju(1)+mju(2)>0) Then
    Do jt = 1, 2
      If (njs(jt)==0) Goto 660
      Do j = 1, 4
        p(n+nrs, j) = p(n+nrs, j) - pjs(jt+2, j)
      End Do
    660 End Do
  End If

!...Produce new particle: side, origin.
  670 i = i + 1
  If (2*i-nsav>=mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSTRF:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If
  jt = int(1.5+rlu(0))
  If (iabs(kfl(3-jt))>10) jt = 3 - jt
  jr = 3 - jt
  js = 3 - 2*jt
  irank(jt) = irank(jt) + 1
  k(i, 1) = 1
  k(i, 3) = ie(jt)
  k(i, 4) = 0
  k(i, 5) = 0

!...Generate flavour, hadron and pT.
  680 Call lukfdi(kfl(jt), 0, kfl(3), k(i,2))
  If (k(i,2)==0) Goto 550
  If (mstj(12)>=3 .And. irank(jt)==1 .And. iabs(kfl(jt))<=10 .And. iabs(kfl(3))>10) Then
    If (rlu(0)>parj(19)) Goto 680
  End If
  p(i, 5) = ulmass(k(i,2))
  Call luptdi(kfl(jt), px(3), py(3))
  pr(jt) = p(i, 5)**2 + (px(jt)+px(3))**2 + (py(jt)+py(3))**2

!...Final hadrons for small invariant mass.
  mstj(93) = 1
  pmq(3) = ulmass(kfl(3))
  wmin = parj(32+mstj(11)) + pmq(1) + pmq(2) + parj(36)*pmq(3)
  If (iabs(kfl(jt))>10 .And. iabs(kfl(3))>10) wmin = wmin - 0.5*parj(36)*pmq(3)
  wrem2 = four(n+nrs, n+nrs)
  If (wrem2<0.10) Goto 550
  If (wrem2<max(wmin*(1.+(2.*rlu(0)-1.)*parj(37)),parj(32)+pmq(1)+pmq(2))**2) Goto 810

!...Choose z, which gives Gamma. Shift z for heavy flavours.
  Call luzdis(kfl(jt), kfl(3), pr(jt), z)

  kfl1a = iabs(kfl(1))
  kfl2a = iabs(kfl(2))
  If (max(mod(kfl1a,10),mod(kfl1a/1000,10),mod(kfl2a,10),mod(kfl2a/1000,10))>=4) Then
    pr(jr) = (pmq(jr)+pmq(3))**2 + (px(jr)-px(3))**2 + (py(jr)-py(3))**2
    pw12 = sqrt(max(0.,(wrem2-pr(1)-pr(2))**2-4.*pr(1)*pr(2)))
    z = (wrem2+pr(jt)-pr(jr)+pw12*(2.*z-1.))/(2.*wrem2)
    pr(jr) = (pmq(jr)+parj(32+mstj(11)))**2 + (px(jr)-px(3))**2 + (py(jr)-py(3))**2
    If ((1.-z)*(wrem2-pr(jt)/z)<pr(jr)) Goto 810
  End If
  gam(3) = (1.-z)*(gam(jt)+pr(jt)/z)
  Do j = 1, 3
    in(j) = in(3*jt+j)
  End Do

!...Stepping within or from 'low' string region easy.
  If (in(1)+1==in(2) .And. z*p(in(1)+2,3)*p(in(2)+2,3)*p(in(1),5)**2>=pr(jt)) Then
    p(in(jt)+2, 4) = z*p(in(jt)+2, 3)
    p(in(jr)+2, 4) = pr(jt)/(p(in(jt)+2,4)*p(in(1),5)**2)
    Do j = 1, 4
      p(i, j) = (px(jt)+px(3))*p(in(3), j) + (py(jt)+py(3))*p(in(3)+1, j)
    End Do
    Goto 770
  Else If (in(1)+1==in(2)) Then
    p(in(jr)+2, 4) = p(in(jr)+2, 3)
    p(in(jr)+2, jt) = 1.
    in(jr) = in(jr) + 4*js
    If (js*in(jr)>js*in(4*jr)) Goto 550
    If (four(in(1),in(2))<=1E-2) Then
      p(in(jt)+2, 4) = p(in(jt)+2, 3)
      p(in(jt)+2, jt) = 0.
      in(jt) = in(jt) + 4*js
    End If
  End If

!...Find new transverse directions (i.e. spacelike string vectors).
  710 If (js*in(1)>js*in(3*jr+1) .Or. js*in(2)>js*in(3*jr+2) .Or. in(1)>in(2)) Goto 550
  If (in(1)/=in(3*jt+1) .Or. in(2)/=in(3*jt+2)) Then
    Do j = 1, 4
      dp(1, j) = dble(p(in(1),j))
      dp(2, j) = dble(p(in(2),j))
      dp(3, j) = 0.D0
      dp(4, j) = 0.D0
    End Do
    dp(1, 4) = dsqrt(dp(1,1)**2+dp(1,2)**2+dp(1,3)**2)
    dp(2, 4) = dsqrt(dp(2,1)**2+dp(2,2)**2+dp(2,3)**2)
    dhc12 = dfour(1, 2)
!lin-5/2012:
!        IF(DHC12.LE.1E-2) THEN
    If (dhc12<=1D-2) Then
      p(in(jt)+2, 4) = p(in(jt)+2, 3)
      p(in(jt)+2, jt) = 0.
      in(jt) = in(jt) + 4*js
      Goto 710
    End If
    in(3) = n + nr + 4*ns + 5
    dp(5, 1) = dp(1, 1)/dp(1, 4) - dp(2, 1)/dp(2, 4)
    dp(5, 2) = dp(1, 2)/dp(1, 4) - dp(2, 2)/dp(2, 4)
    dp(5, 3) = dp(1, 3)/dp(1, 4) - dp(2, 3)/dp(2, 4)
    If (dp(5,1)**2<=dp(5,2)**2+dp(5,3)**2) dp(3, 1) = 1.D0
    If (dp(5,1)**2>dp(5,2)**2+dp(5,3)**2) dp(3, 3) = 1.D0
    If (dp(5,2)**2<=dp(5,1)**2+dp(5,3)**2) dp(4, 2) = 1.D0
    If (dp(5,2)**2>dp(5,1)**2+dp(5,3)**2) dp(4, 3) = 1.D0
    dhcx1 = dfour(3, 1)/dhc12
    dhcx2 = dfour(3, 2)/dhc12
    dhcxx = 1D0/sqrt(1D0+2D0*dhcx1*dhcx2*dhc12)
    dhcy1 = dfour(4, 1)/dhc12
    dhcy2 = dfour(4, 2)/dhc12
    dhcyx = dhcxx*(dhcx1*dhcy2+dhcx2*dhcy1)*dhc12
    dhcyy = 1D0/sqrt(1D0+2D0*dhcy1*dhcy2*dhc12-dhcyx**2)
    Do j = 1, 4
      dp(3, j) = dhcxx*(dp(3,j)-dhcx2*dp(1,j)-dhcx1*dp(2,j))
      p(in(3), j) = sngl(dp(3,j))
      p(in(3)+1, j) = sngl(dhcyy*(dp(4,j)-dhcy2*dp(1,j)-dhcy1*dp(2,j)-dhcyx*dp(3,j)))
    End Do
!...Express pT with respect to new axes, if sensible.
    pxp = -(px(3)*four(in(3*jt+3),in(3))+py(3)*four(in(3*jt+3)+1,in(3)))
    pyp = -(px(3)*four(in(3*jt+3),in(3)+1)+py(3)*four(in(3*jt+3)+1,in(3)+1))
    If (abs(pxp**2+pyp**2-px(3)**2-py(3)**2)<0.01) Then
      px(3) = pxp
      py(3) = pyp
    End If
  End If

!...Sum up known four-momentum. Gives coefficients for m2 expression.
  Do j = 1, 4
    dhg(j) = 0.D0
    p(i, j) = px(jt)*p(in(3*jt+3), j) + py(jt)*p(in(3*jt+3)+1, j) + px(3)*p(in(3), j) + py(3)*p(in(3)+1, j)
    Do in1 = in(3*jt+1), in(1) - 4*js, 4*js
      p(i, j) = p(i, j) + p(in1+2, 3)*p(in1, j)
    End Do
    Do in2 = in(3*jt+2), in(2) - 4*js, 4*js
      p(i, j) = p(i, j) + p(in2+2, 3)*p(in2, j)
    End Do
  End Do
  dhm(1) = dble(four(i,i))
  dhm(2) = dble(2.*four(i,in(1)))
  dhm(3) = dble(2.*four(i,in(2)))
  dhm(4) = dble(2.*four(in(1),in(2)))

!...Find coefficients for Gamma expression.
  Do in2 = in(1) + 1, in(2), 4
    Do in1 = in(1), in2 - 1, 4
      dhc = dble(2.*four(in1,in2))
      dhg(1) = dhg(1) + dble(p(in1+2,jt)*p(in2+2,jt))*dhc
      If (in1==in(1)) dhg(2) = dhg(2) - dble(float(js)*p(in2+2,jt))*dhc
      If (in2==in(2)) dhg(3) = dhg(3) + dble(float(js)*p(in1+2,jt))*dhc
      If (in1==in(1) .And. in2==in(2)) dhg(4) = dhg(4) - dhc
    End Do
  End Do

!...Solve (m2, Gamma) equation system for energies taken.
  dhs1 = dhm(jr+1)*dhg(4) - dhm(4)*dhg(jr+1)
!lin-5/2012:
!      IF(ABS(DHS1).LT.1E-4) GOTO 550
  If (dabs(dhs1)<1D-4) Goto 550
  dhs2 = dhm(4)*(dble(gam(3))-dhg(1)) - dhm(jt+1)*dhg(jr+1) - dhg(4)*(dble(p(i,5))**2-dhm(1)) + dhg(jt+1)*dhm(jr+1)
  dhs3 = dhm(jt+1)*(dble(gam(3))-dhg(1)) - dhg(jt+1)*(dble(p(i,5))**2-dhm(1))
  p(in(jr)+2, 4) = 0.5*sngl((sqrt(max(0D0,dhs2**2-4.D0*dhs1*dhs3)))/abs(dhs1)-dhs2/dhs1)
  If (dhm(jt+1)+dhm(4)*dble(p(in(jr)+2,4))<=0.D0) Goto 550
  p(in(jt)+2, 4) = (p(i,5)**2-sngl(dhm(1))-sngl(dhm(jr+1))*p(in(jr)+2,4))/(sngl(dhm(jt+1))+sngl(dhm(4))*p(in(jr)+2,4))

!...Step to new region if necessary.
  If (p(in(jr)+2,4)>p(in(jr)+2,3)) Then
    p(in(jr)+2, 4) = p(in(jr)+2, 3)
    p(in(jr)+2, jt) = 1.
    in(jr) = in(jr) + 4*js
    If (js*in(jr)>js*in(4*jr)) Goto 550
    If (four(in(1),in(2))<=1E-2) Then
      p(in(jt)+2, 4) = p(in(jt)+2, 3)
      p(in(jt)+2, jt) = 0.
      in(jt) = in(jt) + 4*js
    End If
    Goto 710
  Else If (p(in(jt)+2,4)>p(in(jt)+2,3)) Then
    p(in(jt)+2, 4) = p(in(jt)+2, 3)
    p(in(jt)+2, jt) = 0.
    in(jt) = in(jt) + 4*js
    Goto 710
  End If

!...Four-momentum of particle. Remaining quantities. Loop back.
  770 Do j = 1, 4
    p(i, j) = p(i, j) + p(in(1)+2, 4)*p(in(1), j) + p(in(2)+2, 4)*p(in(2), j)
    p(n+nrs, j) = p(n+nrs, j) - p(i, j)
  End Do
  If (p(i,4)<=0.) Goto 550
  kfl(jt) = -kfl(3)
  pmq(jt) = pmq(3)
  px(jt) = -px(3)
  py(jt) = -py(3)
  gam(jt) = gam(3)
  If (in(3)/=in(3*jt+3)) Then
    Do j = 1, 4
      p(in(3*jt+3), j) = p(in(3), j)
      p(in(3*jt+3)+1, j) = p(in(3)+1, j)
    End Do
  End If
  Do jq = 1, 2
    in(3*jt+jq) = in(jq)
    p(in(jq)+2, 3) = p(in(jq)+2, 3) - p(in(jq)+2, 4)
    p(in(jq)+2, jt) = p(in(jq)+2, jt) - js*(3-2*jq)*p(in(jq)+2, 4)
  End Do
  Goto 670

!...Final hadron: side, flavour, hadron, mass.
  810 i = i + 1
  k(i, 1) = 1
  k(i, 3) = ie(jr)
  k(i, 4) = 0
  k(i, 5) = 0
  Call lukfdi(kfl(jr), -kfl(3), kfldmp, k(i,2))
  If (k(i,2)==0) Goto 550
  p(i, 5) = ulmass(k(i,2))
  pr(jr) = p(i, 5)**2 + (px(jr)-px(3))**2 + (py(jr)-py(3))**2

!...Final two hadrons: find common setup of four-vectors.
  jq = 1
  If (p(in(4)+2,3)*p(in(5)+2,3)*four(in(4),in(5))<p(in(7),3)*p(in(8),3)*four(in(7),in(8))) jq = 2
  dhc12 = dble(four(in(3*jq+1),in(3*jq+2)))
  dhr1 = dble(four(n+nrs,in(3*jq+2)))/dhc12
  dhr2 = dble(four(n+nrs,in(3*jq+1)))/dhc12
  If (in(4)/=in(7) .Or. in(5)/=in(8)) Then
    px(3-jq) = -four(n+nrs, in(3*jq+3)) - px(jq)
    py(3-jq) = -four(n+nrs, in(3*jq+3)+1) - py(jq)
    pr(3-jq) = p(i+(jt+jq-3)**2-1, 5)**2 + (px(3-jq)+(2*jq-3)*js*px(3))**2 + (py(3-jq)+(2*jq-3)*js*py(3))**2
  End If

!...Solve kinematics for final two hadrons, if possible.
  wrem2 = wrem2 + (px(1)+px(2))**2 + (py(1)+py(2))**2
  fd = (sqrt(pr(1))+sqrt(pr(2)))/sqrt(wrem2)
  If (mju(1)+mju(2)/=0 .And. i==isav+2 .And. fd>=1.) Goto 180
  If (fd>=1.) Goto 550
  fa = wrem2 + pr(jt) - pr(jr)
  If (mstj(11)==2) prev = 0.5*fd**parj(37+mstj(11))
  If (mstj(11)/=2) prev = 0.5*exp(max(-100.,log(fd)*parj(37+mstj(11))*(pr(1)+pr(2))**2))
  fb = sign(sqrt(max(0.,fa**2-4.*wrem2*pr(jt))), js*(rlu(0)-prev))
  kfl1a = iabs(kfl(1))
  kfl2a = iabs(kfl(2))
  If (max(mod(kfl1a,10),mod(kfl1a/1000,10),mod(kfl2a,10),mod(kfl2a/1000,10))>=6) fb = sign(sqrt(max(0.,fa**2-4.*wrem2*pr(jt))), float(js))
  Do j = 1, 4
    p(i-1, j) = (px(jt)+px(3))*p(in(3*jq+3), j) + (py(jt)+py(3))*p(in(3*jq+3)+1, j) + 0.5*(sngl(dhr1)*(fa+fb)*p(in(3*jq+1),j)+sngl(dhr2)*(fa-fb)*p(in(3*jq+2),j))/wrem2
    p(i, j) = p(n+nrs, j) - p(i-1, j)
  End Do

!...Mark jets as fragmented and give daughter pointers.
  n = i - nrs + 1
  Do i = nsav + 1, nsav + np
    im = k(i, 3)
    k(im, 1) = k(im, 1) + 10
    If (mstu(16)/=2) Then
      k(im, 4) = nsav + 1
      k(im, 5) = nsav + 1
    Else
      k(im, 4) = nsav + 2
      k(im, 5) = n
    End If
  End Do

!...Document string system. Move up particles.
  nsav = nsav + 1
  k(nsav, 1) = 11
  k(nsav, 2) = 92
  k(nsav, 3) = ip
  k(nsav, 4) = nsav + 1
  k(nsav, 5) = n
  Do j = 1, 4
    p(nsav, j) = sngl(dps(j))
    v(nsav, j) = v(ip, j)
  End Do
  p(nsav, 5) = sqrt(sngl(max(0D0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2)))
  v(nsav, 5) = 0.
  Do i = nsav + 1, n

    Do j = 1, 5
      k(i, j) = k(i+nrs-1, j)
      p(i, j) = p(i+nrs-1, j)
      v(i, j) = 0.
    End Do
  End Do

!...Order particles in rank along the chain. Update mother pointer.
  Do i = nsav + 1, n
    Do j = 1, 5
      k(i-nsav+n, j) = k(i, j)
      p(i-nsav+n, j) = p(i, j)
    End Do
  End Do
  i1 = nsav
  Do i = n + 1, 2*n - nsav
    If (k(i,3)/=ie(1)) Goto 880
    i1 = i1 + 1
    Do j = 1, 5
      k(i1, j) = k(i, j)
      p(i1, j) = p(i, j)
    End Do
    If (mstu(16)/=2) k(i1, 3) = nsav
  880 End Do
  Do i = 2*n - nsav, n + 1, -1
    If (k(i,3)==ie(1)) Goto 900
    i1 = i1 + 1
    Do j = 1, 5
      k(i1, j) = k(i, j)
      p(i1, j) = p(i, j)
    End Do
    If (mstu(16)/=2) k(i1, 3) = nsav
  900 End Do

!...Boost back particle system. Set production vertices.
  Call ludbrb(nsav+1, n, 0., 0., dps(1)/dps(4), dps(2)/dps(4), dps(3)/dps(4))
  Do i = nsav + 1, n

    Do j = 1, 4
      v(i, j) = v(ip, j)
    End Do
  End Do

  Return
End Subroutine lustrf

!*********************************************************************

Subroutine luindf(ip)

!...Purpose: to handle the fragmentation of a jet system (or a single
!...jet) according to independent fragmentation models.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension dps(5), psi(4), nfi(3), nfl(3), ifet(3), kflf(3), kflo(2), pxo(2), pyo(2), wo(2)

!...Reset counters. Identify parton system and take copy. Check flavour.
  nsav = n
  njet = 0
  kqsum = 0
  Do j = 1, 5
    dps(j) = 0.D0
  End Do
  i = ip - 1
  110 i = i + 1
  If (i>min(n,mstu(4)-mstu(32))) Then
    Call luerrm(12, '(LUINDF:) failed to reconstruct jet system')
    If (mstu(21)>=1) Return
  End If
  If (k(i,1)/=1 .And. k(i,1)/=2) Goto 110
  kc = lucomp(k(i,2))
  If (kc==0) Goto 110
  kq = kchg(kc, 2)*isign(1, k(i,2))
  If (kq==0) Goto 110
  njet = njet + 1
  If (kq/=2) kqsum = kqsum + kq
  Do j = 1, 5
    k(nsav+njet, j) = k(i, j)
    p(nsav+njet, j) = p(i, j)
    dps(j) = dps(j) + dble(p(i,j))
  End Do
  k(nsav+njet, 3) = i
  If (k(i,1)==2 .Or. (mstj(3)<=5 .And. n>i .And. k(i+1,1)==2)) Goto 110
  If (njet/=1 .And. kqsum/=0) Then
    Call luerrm(12, '(LUINDF:) unphysical flavour combination')
    If (mstu(21)>=1) Return
  End If

!...Boost copied system to CM frame. Find CM energy and sum flavours.
  If (njet/=1) Call ludbrb(nsav+1, nsav+njet, 0., 0., -dps(1)/dps(4), -dps(2)/dps(4), -dps(3)/dps(4))
  pecm = 0.
  Do j = 1, 3
    nfi(j) = 0
  End Do
  Do i = nsav + 1, nsav + njet
    pecm = pecm + p(i, 4)
    kfa = iabs(k(i,2))
    If (kfa<=3) Then
      nfi(kfa) = nfi(kfa) + isign(1, k(i,2))
    Else If (kfa>1000) Then
      kfla = mod(kfa/1000, 10)
      kflb = mod(kfa/100, 10)
      If (kfla<=3) nfi(kfla) = nfi(kfla) + isign(1, k(i,2))
      If (kflb<=3) nfi(kflb) = nfi(kflb) + isign(1, k(i,2))
    End If
  End Do

!...Loop over attempts made. Reset counters.
  ntry = 0
  150 ntry = ntry + 1
  n = nsav + njet
  If (ntry>200) Then
    Call luerrm(14, '(LUINDF:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  Do j = 1, 3
    nfl(j) = nfi(j)
    ifet(j) = 0
    kflf(j) = 0
  End Do

!...Loop over jets to be fragmented.
  Do ip1 = nsav + 1, nsav + njet
    mstj(91) = 0
    nsav1 = n

!...Initial flavour and momentum values. Jet along +z axis.
    kflh = iabs(k(ip1,2))
    If (kflh>10) kflh = mod(kflh/1000, 10)
    kflo(2) = 0
    wf = p(ip1, 4) + sqrt(p(ip1,1)**2+p(ip1,2)**2+p(ip1,3)**2)

!...Initial values for quark or diquark jet.
    170 If (iabs(k(ip1,2))/=21) Then
      nstr = 1
      kflo(1) = k(ip1, 2)
      Call luptdi(0, pxo(1), pyo(1))
      wo(1) = wf

!...Initial values for gluon treated like random quark jet.
    Else If (mstj(2)<=2) Then
      nstr = 1
      If (mstj(2)==2) mstj(91) = 1
      kflo(1) = int(1.+(2.+parj(2))*rlu(0))*(-1)**int(rlu(0)+0.5)
      Call luptdi(0, pxo(1), pyo(1))
      wo(1) = wf

!...Initial values for gluon treated like quark-antiquark jet pair,
!...sharing energy according to Altarelli-Parisi splitting function.
    Else
      nstr = 2
      If (mstj(2)==4) mstj(91) = 1
      kflo(1) = int(1.+(2.+parj(2))*rlu(0))*(-1)**int(rlu(0)+0.5)
      kflo(2) = -kflo(1)
      Call luptdi(0, pxo(1), pyo(1))
      pxo(2) = -pxo(1)
      pyo(2) = -pyo(1)
      wo(1) = wf*rlu(0)**(1./3.)
      wo(2) = wf - wo(1)
    End If

!...Initial values for rank, flavour, pT and W+.
    Do istr = 1, nstr
      180 i = n
      irank = 0
      kfl1 = kflo(istr)
      px1 = pxo(istr)
      py1 = pyo(istr)
      w = wo(istr)

!...New hadron. Generate flavour and hadron species.
      190 i = i + 1
      If (i>=mstu(4)-mstu(32)-njet-5) Then
        Call luerrm(11, '(LUINDF:) no more memory left in LUJETS')
        If (mstu(21)>=1) Return
      End If
      irank = irank + 1
      k(i, 1) = 1
      k(i, 3) = ip1
      k(i, 4) = 0
      k(i, 5) = 0
      200 Call lukfdi(kfl1, 0, kfl2, k(i,2))
      If (k(i,2)==0) Goto 180
      If (mstj(12)>=3 .And. irank==1 .And. iabs(kfl1)<=10 .And. iabs(kfl2)>10) Then
        If (rlu(0)>parj(19)) Goto 200
      End If

!...Find hadron mass. Generate four-momentum.
      p(i, 5) = ulmass(k(i,2))
      Call luptdi(kfl1, px2, py2)
      p(i, 1) = px1 + px2
      p(i, 2) = py1 + py2
      pr = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
      Call luzdis(kfl1, kfl2, pr, z)
      p(i, 3) = 0.5*(z*w-pr/(z*w))
      p(i, 4) = 0.5*(z*w+pr/(z*w))
      If (mstj(3)>=1 .And. irank==1 .And. kflh>=4 .And. p(i,3)<=0.001) Then
        If (w>=p(i,5)+0.5*parj(32)) Goto 180
        p(i, 3) = 0.0001
        p(i, 4) = sqrt(pr)
        z = p(i, 4)/w
      End If

!...Remaining flavour and momentum.
      kfl1 = -kfl2
      px1 = -px2
      py1 = -py2
      w = (1.-z)*w
      Do j = 1, 5
        v(i, j) = 0.
      End Do

!...Check if pL acceptable. Go back for new hadron if enough energy.
      If (mstj(3)>=0 .And. p(i,3)<0.) i = i - 1
      If (w>parj(31)) Goto 190
      n = i
    End Do
    If (mod(mstj(3),5)==4 .And. n==nsav1) wf = wf + 0.1*parj(32)
    If (mod(mstj(3),5)==4 .And. n==nsav1) Goto 170

!...Rotate jet to new direction.
    the = ulangl(p(ip1,3), sqrt(p(ip1,1)**2+p(ip1,2)**2))
    phi = ulangl(p(ip1,1), p(ip1,2))
    Call ludbrb(nsav1+1, n, the, phi, 0D0, 0D0, 0D0)
    k(k(ip1,3), 4) = nsav1 + 1
    k(k(ip1,3), 5) = n

!...End of jet generation loop. Skip conservation in some cases.
  End Do
  If (njet==1 .Or. mstj(3)<=0) Goto 470
  If (mod(mstj(3),5)/=0 .And. n-nsav-njet<2) Goto 150

!...Subtract off produced hadron flavours, finished if zero.
  Do i = nsav + njet + 1, n
    kfa = iabs(k(i,2))
    kfla = mod(kfa/1000, 10)
    kflb = mod(kfa/100, 10)
    kflc = mod(kfa/10, 10)
    If (kfla==0) Then
      If (kflb<=3) nfl(kflb) = nfl(kflb) - isign(1, k(i,2))*(-1)**kflb
      If (kflc<=3) nfl(kflc) = nfl(kflc) + isign(1, k(i,2))*(-1)**kflb
    Else
      If (kfla<=3) nfl(kfla) = nfl(kfla) - isign(1, k(i,2))
      If (kflb<=3) nfl(kflb) = nfl(kflb) - isign(1, k(i,2))
      If (kflc<=3) nfl(kflc) = nfl(kflc) - isign(1, k(i,2))
    End If
  End Do
  nreq = (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+nfl(2)+nfl(3)))/2 + iabs(nfl(1)+nfl(2)+nfl(3))/3
  If (nreq==0) Goto 320

!...Take away flavour of low-momentum particles until enough freedom.
  nrem = 0
  250 irem = 0
  p2min = pecm**2
  Do i = nsav + njet + 1, n
    p2 = p(i, 1)**2 + p(i, 2)**2 + p(i, 3)**2
    If (k(i,1)==1 .And. p2<p2min) irem = i
    If (k(i,1)==1 .And. p2<p2min) p2min = p2
  End Do
  If (irem==0) Goto 150
  k(irem, 1) = 7
  kfa = iabs(k(irem,2))
  kfla = mod(kfa/1000, 10)
  kflb = mod(kfa/100, 10)
  kflc = mod(kfa/10, 10)
  If (kfla>=4 .Or. kflb>=4) k(irem, 1) = 8
  If (k(irem,1)==8) Goto 250
  If (kfla==0) Then
    isgn = isign(1, k(irem,2))*(-1)**kflb
    If (kflb<=3) nfl(kflb) = nfl(kflb) + isgn
    If (kflc<=3) nfl(kflc) = nfl(kflc) - isgn
  Else
    If (kfla<=3) nfl(kfla) = nfl(kfla) + isign(1, k(irem,2))
    If (kflb<=3) nfl(kflb) = nfl(kflb) + isign(1, k(irem,2))
    If (kflc<=3) nfl(kflc) = nfl(kflc) + isign(1, k(irem,2))
  End If
  nrem = nrem + 1
  nreq = (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+nfl(2)+nfl(3)))/2 + iabs(nfl(1)+nfl(2)+nfl(3))/3
  If (nreq>nrem) Goto 250
  Do i = nsav + njet + 1, n
    If (k(i,1)==8) k(i, 1) = 1
  End Do

!...Find combination of existing and new flavours for hadron.
  280 nfet = 2
  If (nfl(1)+nfl(2)+nfl(3)/=0) nfet = 3
  If (nreq<nrem) nfet = 1
  If (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))==0) nfet = 0
  Do j = 1, nfet
    ifet(j) = 1 + int((iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3)))*rlu(0))
    kflf(j) = isign(1, nfl(1))
    If (ifet(j)>iabs(nfl(1))) kflf(j) = isign(2, nfl(2))
    If (ifet(j)>iabs(nfl(1))+iabs(nfl(2))) kflf(j) = isign(3, nfl(3))
  End Do
  If (nfet==2 .And. (ifet(1)==ifet(2) .Or. kflf(1)*kflf(2)>0)) Goto 280
  If (nfet==3 .And. (ifet(1)==ifet(2) .Or. ifet(1)==ifet(3) .Or. ifet(2)==ifet(3) .Or. kflf(1)*kflf(2)<0 .Or. kflf(1)*kflf(3)<0 .Or. kflf(1)*(nfl(1)+nfl(2)+nfl(3))<0)) Goto 280
  If (nfet==0) kflf(1) = 1 + int((2.+parj(2))*rlu(0))
  If (nfet==0) kflf(2) = -kflf(1)
  If (nfet==1) kflf(2) = isign(1+int((2.+parj(2))*rlu(0)), -kflf(1))
  If (nfet<=2) kflf(3) = 0
  If (kflf(3)/=0) Then
    kflfc = isign(1000*max(iabs(kflf(1)),iabs(kflf(3)))+100*min(iabs(kflf(1)),iabs(kflf(3)))+1, kflf(1))
    If (kflf(1)==kflf(3) .Or. (1.+3.*parj(4))*rlu(0)>1.) kflfc = kflfc + isign(2, kflfc)
  Else
    kflfc = kflf(1)
  End If
  Call lukfdi(kflfc, kflf(2), kfldmp, kf)
  If (kf==0) Goto 280
  Do j = 1, max(2, nfet)
    nfl(iabs(kflf(j))) = nfl(iabs(kflf(j))) - isign(1, kflf(j))
  End Do

!...Store hadron at random among free positions.
  npos = min(1+int(rlu(0)*nrem), nrem)
  Do i = nsav + njet + 1, n
    If (k(i,1)==7) npos = npos - 1
    If (k(i,1)==1 .Or. npos/=0) Goto 310
    k(i, 1) = 1
    k(i, 2) = kf
    p(i, 5) = ulmass(k(i,2))
    p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
  310 End Do
  nrem = nrem - 1
  nreq = (iabs(nfl(1))+iabs(nfl(2))+iabs(nfl(3))-iabs(nfl(1)+nfl(2)+nfl(3)))/2 + iabs(nfl(1)+nfl(2)+nfl(3))/3
  If (nrem>0) Goto 280

!...Compensate for missing momentum in global scheme (3 options).
  320 If (mod(mstj(3),5)/=0 .And. mod(mstj(3),5)/=4) Then
    Do j = 1, 3
      psi(j) = 0.
      Do i = nsav + njet + 1, n
        psi(j) = psi(j) + p(i, j)
      End Do
    End Do
    psi(4) = psi(1)**2 + psi(2)**2 + psi(3)**2
    pws = 0.
    Do i = nsav + njet + 1, n
      If (mod(mstj(3),5)==1) pws = pws + p(i, 4)
      If (mod(mstj(3),5)==2) pws = pws + sqrt(p(i,5)**2+(psi(1)*p(i,1)+psi(2)*p(i,2)+psi(3)*p(i,3))**2/psi(4))
      If (mod(mstj(3),5)==3) pws = pws + 1.
    End Do
    Do i = nsav + njet + 1, n
      If (mod(mstj(3),5)==1) pw = p(i, 4)
      If (mod(mstj(3),5)==2) pw = sqrt(p(i,5)**2+(psi(1)*p(i,1)+psi(2)*p(i,2)+psi(3)*p(i,3))**2/psi(4))
      If (mod(mstj(3),5)==3) pw = 1.
      Do j = 1, 3
        p(i, j) = p(i, j) - psi(j)*pw/pws
      End Do
      p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
    End Do

!...Compensate for missing momentum withing each jet separately.
  Else If (mod(mstj(3),5)==4) Then
    Do i = n + 1, n + njet
      k(i, 1) = 0
      Do j = 1, 5
        p(i, j) = 0.
      End Do
    End Do
    Do i = nsav + njet + 1, n
      ir1 = k(i, 3)
      ir2 = n + ir1 - nsav
      k(ir2, 1) = k(ir2, 1) + 1
      pls = (p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/(p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
      Do j = 1, 3
        p(ir2, j) = p(ir2, j) + p(i, j) - pls*p(ir1, j)
      End Do
      p(ir2, 4) = p(ir2, 4) + p(i, 4)
      p(ir2, 5) = p(ir2, 5) + pls
    End Do
    pss = 0.
    Do i = n + 1, n + njet
      If (k(i,1)/=0) pss = pss + p(i, 4)/(pecm*(0.8*p(i,5)+0.2))
    End Do
    Do i = nsav + njet + 1, n
      ir1 = k(i, 3)
      ir2 = n + ir1 - nsav
      pls = (p(i,1)*p(ir1,1)+p(i,2)*p(ir1,2)+p(i,3)*p(ir1,3))/(p(ir1,1)**2+p(ir1,2)**2+p(ir1,3)**2)
      Do j = 1, 3
        p(i, j) = p(i, j) - p(ir2, j)/k(ir2, 1) + (1./(p(ir2,5)*pss)-1.)*pls*p(ir1, j)
      End Do
      p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
    End Do
  End If

!...Scale momenta for energy conservation.
  If (mod(mstj(3),5)/=0) Then
    pms = 0.
    pes = 0.
    pqs = 0.
    Do i = nsav + njet + 1, n
      pms = pms + p(i, 5)
      pes = pes + p(i, 4)
      pqs = pqs + p(i, 5)**2/p(i, 4)
    End Do
    If (pms>=pecm) Goto 150
    neco = 0
    440 neco = neco + 1
    pfac = (pecm-pqs)/(pes-pqs)
    pes = 0.
    pqs = 0.
    Do i = nsav + njet + 1, n
      Do j = 1, 3
        p(i, j) = pfac*p(i, j)
      End Do
      p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
      pes = pes + p(i, 4)
      pqs = pqs + p(i, 5)**2/p(i, 4)
    End Do
    If (neco<10 .And. abs(pecm-pes)>2E-6*pecm) Goto 440
  End If

!...Origin of produced particles and parton daughter pointers.
  470 Do i = nsav + njet + 1, n
    If (mstu(16)/=2) k(i, 3) = nsav + 1
    If (mstu(16)==2) k(i, 3) = k(k(i,3), 3)
  End Do
  Do i = nsav + 1, nsav + njet
    i1 = k(i, 3)
    k(i1, 1) = k(i1, 1) + 10
    If (mstu(16)/=2) Then
      k(i1, 4) = nsav + 1
      k(i1, 5) = nsav + 1
    Else
      k(i1, 4) = k(i1, 4) - njet + 1
      k(i1, 5) = k(i1, 5) - njet + 1
      If (k(i1,5)<k(i1,4)) Then
        k(i1, 4) = 0
        k(i1, 5) = 0
      End If
    End If
  End Do

!...Document independent fragmentation system. Remove copy of jets.
  nsav = nsav + 1
  k(nsav, 1) = 11
  k(nsav, 2) = 93
  k(nsav, 3) = ip
  k(nsav, 4) = nsav + 1
  k(nsav, 5) = n - njet + 1
  Do j = 1, 4
    p(nsav, j) = sngl(dps(j))
    v(nsav, j) = v(ip, j)
  End Do
  p(nsav, 5) = sqrt(sngl(max(0D0,dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2)))
  v(nsav, 5) = 0.
  Do i = nsav + njet, n
    Do j = 1, 5
      k(i-njet+1, j) = k(i, j)
      p(i-njet+1, j) = p(i, j)
      v(i-njet+1, j) = v(i, j)
    End Do
  End Do
  n = n - njet + 1

!...Boost back particle system. Set production vertices.
  If (njet/=1) Call ludbrb(nsav+1, n, 0., 0., dps(1)/dps(4), dps(2)/dps(4), dps(3)/dps(4))
  Do i = nsav + 1, n
    Do j = 1, 4
      v(i, j) = v(ip, j)
    End Do
  End Do

  Return
End Subroutine luindf

!*********************************************************************


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine ludecy(ip)

!...Purpose: to handle the decay of unstable particles.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Dimension vdcy(4), kflo(4), kfl1(4), pv(10, 5), rord(10), ue(3), be(3), wtcor(10)
!lin-2/18/03 for resonance decay in hadron cascade:
  Common /resdcy/nsav, iksdcy
  Save /resdcy/
  Data wtcor/2., 5., 15., 60., 250., 1500., 1.2E4, 1.2E5, 150., 16./

!...Functions: momentum in two-particle decays, four-product and
!...matrix element times phase space in weak decays.
  pawt(a, b, c) = sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.*a)
  four(i, j) = p(i, 4)*p(j, 4) - p(i, 1)*p(j, 1) - p(i, 2)*p(j, 2) - p(i, 3)*p(j, 3)
  hmeps(ha) = ((1.-hrq-ha)**2+3.*ha*(1.+hrq-ha))*sqrt((1.-hrq-ha)**2-4.*hrq*ha)

!...Initial values.
  ntry = 0
  nsav = n
  kfa = iabs(k(ip,2))
  kfs = isign(1, k(ip,2))
  kc = lucomp(kfa)
  mstj(92) = 0

!...Choose lifetime and determine decay vertex.
  If (k(ip,1)==5) Then
    v(ip, 5) = 0.
  Else If (k(ip,1)/=4) Then
    v(ip, 5) = -pmas(kc, 4)*log(rlu(0))
  End If
  Do j = 1, 4
    vdcy(j) = v(ip, j) + v(ip, 5)*p(ip, j)/p(ip, 5)
  End Do

!...Determine whether decay allowed or not.
  mout = 0
  If (mstj(22)==2) Then
    If (pmas(kc,4)>parj(71)) mout = 1
  Else If (mstj(22)==3) Then
    If (vdcy(1)**2+vdcy(2)**2+vdcy(3)**2>parj(72)**2) mout = 1
  Else If (mstj(22)==4) Then
    If (vdcy(1)**2+vdcy(2)**2>parj(73)**2) mout = 1
    If (abs(vdcy(3))>parj(74)) mout = 1
  End If
  If (mout==1 .And. k(ip,1)/=5) Then
    k(ip, 1) = 4
    Return
  End If

!...Check existence of decay channels. Particle/antiparticle rules.
  kca = kc
  If (mdcy(kc,2)>0) Then
    mdmdcy = mdme(mdcy(kc,2), 2)
    If (mdmdcy>80 .And. mdmdcy<=90) kca = mdmdcy
  End If
  If (mdcy(kca,2)<=0 .Or. mdcy(kca,3)<=0) Then
    Call luerrm(9, '(LUDECY:) no decay channel defined')
    Return
  End If
  If (mod(kfa/1000,10)==0 .And. (kca==85 .Or. kca==87)) kfs = -kfs
  If (kchg(kc,3)==0) Then
    kfsp = 1
    kfsn = 0
    If (rlu(0)>0.5) kfs = -kfs
  Else If (kfs>0) Then
    kfsp = 1
    kfsn = 0
  Else
    kfsp = 0
    kfsn = 1
  End If

!...Sum branching ratios of allowed decay channels.
!lin  110 NOPE=0
  nope = 0
  brsu = 0.
  Do idl = mdcy(kca, 2), mdcy(kca, 2) + mdcy(kca, 3) - 1
    If (mdme(idl,1)/=1 .And. kfsp*mdme(idl,1)/=2 .And. kfsn*mdme(idl,1)/=3) Goto 120
    If (mdme(idl,2)>100) Goto 120
    nope = nope + 1
    brsu = brsu + brat(idl)
  120 End Do
  If (nope==0) Then
    Call luerrm(2, '(LUDECY:) all decay channels closed by user')
    Return
  End If

!...Select decay channel among allowed ones.
  130 rbr = brsu*rlu(0)
  idl = mdcy(kca, 2) - 1
  140 idl = idl + 1
  If (mdme(idl,1)/=1 .And. kfsp*mdme(idl,1)/=2 .And. kfsn*mdme(idl,1)/=3) Then
    If (idl<mdcy(kca,2)+mdcy(kca,3)-1) Goto 140
  Else If (mdme(idl,2)>100) Then
    If (idl<mdcy(kca,2)+mdcy(kca,3)-1) Goto 140
  Else
    idc = idl
    rbr = rbr - brat(idl)
    If (idl<mdcy(kca,2)+mdcy(kca,3)-1 .And. rbr>0.) Goto 140
  End If

!...Start readout of decay channel: matrix element, reset counters.
  mmat = mdme(idc, 2)
  150 ntry = ntry + 1
  If (ntry>1000) Then
    Call luerrm(14, '(LUDECY:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  i = n
  np = 0
  nq = 0
  mbst = 0
  If (mmat>=11 .And. mmat/=46 .And. p(ip,4)>20.*p(ip,5)) mbst = 1
  Do j = 1, 4
    pv(1, j) = 0.
    If (mbst==0) pv(1, j) = p(ip, j)
  End Do
  If (mbst==1) pv(1, 4) = p(ip, 5)
  pv(1, 5) = p(ip, 5)
  ps = 0.
  psq = 0.
  mrem = 0

!...Read out decay products. Convert to standard flavour code.
  jtmax = 5
  If (mdme(idc+1,2)==101) jtmax = 10
  Do jt = 1, jtmax
    If (jt<=5) kp = kfdp(idc, jt)
    If (jt>=6) kp = kfdp(idc+1, jt-5)
    If (kp==0) Goto 170
    kpa = iabs(kp)
    kcp = lucomp(kpa)
    If (kchg(kcp,3)==0 .And. kpa/=81 .And. kpa/=82) Then
      kfp = kp
    Else If (kpa/=81 .And. kpa/=82) Then
      kfp = kfs*kp
    Else If (kpa==81 .And. mod(kfa/1000,10)==0) Then
      kfp = -kfs*mod(kfa/10, 10)
    Else If (kpa==81 .And. mod(kfa/100,10)>=mod(kfa/10,10)) Then
      kfp = kfs*(100*mod(kfa/10,100)+3)
    Else If (kpa==81) Then
      kfp = kfs*(1000*mod(kfa/10,10)+100*mod(kfa/100,10)+1)
    Else If (kp==82) Then
      Call lukfdi(-kfs*int(1.+(2.+parj(2))*rlu(0)), 0, kfp, kdump)
      If (kfp==0) Goto 150
      mstj(93) = 1
      If (pv(1,5)<parj(32)+2.*ulmass(kfp)) Goto 150
    Else If (kp==-82) Then
      kfp = -kfp
      If (iabs(kfp)>10) kfp = kfp + isign(10000, kfp)
    End If
    If (kpa==81 .Or. kpa==82) kcp = lucomp(kfp)

!...Add decay product to event record or to quark flavour list.
    kfpa = iabs(kfp)
    kqp = kchg(kcp, 2)
    If (mmat>=11 .And. mmat<=30 .And. kqp/=0) Then
      nq = nq + 1
      kflo(nq) = kfp
      mstj(93) = 2
      psq = psq + ulmass(kflo(nq))
    Else If (mmat>=42 .And. mmat<=43 .And. np==3 .And. mod(nq,2)==1) Then
      nq = nq - 1
      ps = ps - p(i, 5)
      k(i, 1) = 1
      kfi = k(i, 2)
      Call lukfdi(kfp, kfi, kfldmp, k(i,2))
      If (k(i,2)==0) Goto 150
      mstj(93) = 1
      p(i, 5) = ulmass(k(i,2))
      ps = ps + p(i, 5)
    Else
      i = i + 1
      np = np + 1
      If (mmat/=33 .And. kqp/=0) nq = nq + 1
      If (mmat==33 .And. kqp/=0 .And. kqp/=2) nq = nq + 1
      k(i, 1) = 1 + mod(nq, 2)
      If (mmat==4 .And. jt<=2 .And. kfp==21) k(i, 1) = 2
      If (mmat==4 .And. jt==3) k(i, 1) = 1
      k(i, 2) = kfp
      k(i, 3) = ip
      k(i, 4) = 0
      k(i, 5) = 0
      p(i, 5) = ulmass(kfp)
      If (mmat==45 .And. kfpa==89) p(i, 5) = parj(32)
      ps = ps + p(i, 5)
    End If
  170 End Do

!...Choose decay multiplicity in phase space model.
  180 If (mmat>=11 .And. mmat<=30) Then
    psp = ps
    cnde = parj(61)*log(max((pv(1,5)-ps-psq)/parj(62),1.1))
    If (mmat==12) cnde = cnde + parj(63)
    190 ntry = ntry + 1
    If (ntry>1000) Then
      Call luerrm(14, '(LUDECY:) caught in infinite loop')
      If (mstu(21)>=1) Return
    End If
    If (mmat<=20) Then
      gauss = sqrt(-2.*cnde*log(max(1E-10,rlu(0))))*sin(paru(2)*rlu(0))
      nd = int(0.5+0.5*np+0.25*nq+cnde+gauss)
      If (nd<np+nq/2 .Or. nd<2 .Or. nd>10) Goto 190
      If (mmat==13 .And. nd==2) Goto 190
      If (mmat==14 .And. nd<=3) Goto 190
      If (mmat==15 .And. nd<=4) Goto 190
    Else
      nd = mmat - 20
    End If

!...Form hadrons from flavour content.
    Do jt = 1, 4
      kfl1(jt) = kflo(jt)
    End Do
    If (nd==np+nq/2) Goto 220
    Do i = n + np + 1, n + nd - nq/2
      jt = 1 + int((nq-1)*rlu(0))
      Call lukfdi(kfl1(jt), 0, kfl2, k(i,2))
      If (k(i,2)==0) Goto 190
      kfl1(jt) = -kfl2
    End Do
    220 jt = 2
    jt2 = 3
    jt3 = 4
    If (nq==4 .And. rlu(0)<parj(66)) jt = 4
    If (jt==4 .And. isign(1,kfl1(1)*(10-iabs(kfl1(1))))*isign(1,kfl1(jt)*(10-iabs(kfl1(jt))))>0) jt = 3
    If (jt==3) jt2 = 2
    If (jt==4) jt3 = 2
    Call lukfdi(kfl1(1), kfl1(jt), kfldmp, k(n+nd-nq/2+1,2))
    If (k(n+nd-nq/2+1,2)==0) Goto 190
    If (nq==4) Call lukfdi(kfl1(jt2), kfl1(jt3), kfldmp, k(n+nd,2))
    If (nq==4 .And. k(n+nd,2)==0) Goto 190

!...Check that sum of decay product masses not too large.
    ps = psp
    Do i = n + np + 1, n + nd
      k(i, 1) = 1
      k(i, 3) = ip
      k(i, 4) = 0
      k(i, 5) = 0
      p(i, 5) = ulmass(k(i,2))
      ps = ps + p(i, 5)
    End Do
    If (ps+parj(64)>pv(1,5)) Goto 190

!...Rescale energy to subtract off spectator quark mass.
  Else If ((mmat==31 .Or. mmat==33 .Or. mmat==44 .Or. mmat==45) .And. np>=3) Then
    ps = ps - p(n+np, 5)
    pqt = (p(n+np,5)+parj(65))/pv(1, 5)
    Do j = 1, 5
      p(n+np, j) = pqt*pv(1, j)
      pv(1, j) = (1.-pqt)*pv(1, j)
    End Do
    If (ps+parj(64)>pv(1,5)) Goto 150
    nd = np - 1
    mrem = 1

!...Phase space factors imposed in W decay.
  Else If (mmat==46) Then
    mstj(93) = 1
    psmc = ulmass(k(n+1,2))
    mstj(93) = 1
    psmc = psmc + ulmass(k(n+2,2))
    If (max(ps,psmc)+parj(32)>pv(1,5)) Goto 130
    hr1 = (p(n+1,5)/pv(1,5))**2
    hr2 = (p(n+2,5)/pv(1,5))**2
    If ((1.-hr1-hr2)*(2.+hr1+hr2)*sqrt((1.-hr1-hr2)**2-4.*hr1*hr2)<2.*rlu(0)) Goto 130
    nd = np

!...Fully specified final state: check mass broadening effects.
  Else
    If (np>=2 .And. ps+parj(64)>pv(1,5)) Goto 150
    nd = np
  End If

!...Select W mass in decay Q -> W + q, without W propagator.
  If (mmat==45 .And. mstj(25)<=0) Then
    hlq = (parj(32)/pv(1,5))**2
    huq = (1.-(p(n+2,5)+parj(64))/pv(1,5))**2
    hrq = (p(n+2,5)/pv(1,5))**2
    250 hw = hlq + rlu(0)*(huq-hlq)
    If (hmeps(hw)<rlu(0)) Goto 250
    p(n+1, 5) = pv(1, 5)*sqrt(hw)

!...Ditto, including W propagator. Divide mass range into three regions.
  Else If (mmat==45) Then
    hqw = (pv(1,5)/pmas(24,1))**2
    hlw = (parj(32)/pmas(24,1))**2
    huw = ((pv(1,5)-p(n+2,5)-parj(64))/pmas(24,1))**2
    hrq = (p(n+2,5)/pv(1,5))**2
    hg = pmas(24, 2)/pmas(24, 1)
    hatl = atan((hlw-1.)/hg)
    hm = min(1., huw-0.001)
    hmv1 = hmeps(hm/hqw)/((hm-1.)**2+hg**2)
    260 hm = hm - hg
    hmv2 = hmeps(hm/hqw)/((hm-1.)**2+hg**2)
    hsav1 = hmeps(hm/hqw)
    hsav2 = 1./((hm-1.)**2+hg**2)
    If (hmv2>hmv1 .And. hm-hg>hlw) Then
      hmv1 = hmv2
      Goto 260
    End If
    hmv = min(2.*hmv1, hmeps(hm/hqw)/hg**2)
    hm1 = 1. - sqrt(1./hmv-hg**2)
    If (hm1>hlw .And. hm1<hm) Then
      hm = hm1
    Else If (hmv2<=hmv1) Then
      hm = max(hlw, hm-min(0.1,1.-hm))
    End If
    hatm = atan((hm-1.)/hg)
    hwt1 = (hatm-hatl)/hg
    hwt2 = hmv*(min(1.,huw)-hm)
    hwt3 = 0.
    If (huw>1.) Then
      hatu = atan((huw-1.)/hg)
      hmp1 = hmeps(1./hqw)
      hwt3 = hmp1*hatu/hg
    End If

!...Select mass region and W mass there. Accept according to weight.
    270 hreg = rlu(0)*(hwt1+hwt2+hwt3)
    If (hreg<=hwt1) Then
      hw = 1. + hg*tan(hatl+rlu(0)*(hatm-hatl))
      hacc = hmeps(hw/hqw)
    Else If (hreg<=hwt1+hwt2) Then
      hw = hm + rlu(0)*(min(1.,huw)-hm)
      hacc = hmeps(hw/hqw)/((hw-1.)**2+hg**2)/hmv
    Else
      hw = 1. + hg*tan(rlu(0)*hatu)
      hacc = hmeps(hw/hqw)/hmp1
    End If
    If (hacc<rlu(0)) Goto 270
    p(n+1, 5) = pmas(24, 1)*sqrt(hw)
  End If

!...Determine position of grandmother, number of sisters, Q -> W sign.
  nm = 0
  msgn = 0
  If (mmat==3 .Or. mmat==46) Then
    im = k(ip, 3)
    If (im<0 .Or. im>=ip) im = 0
    If (im/=0) kfam = iabs(k(im,2))
    If (im/=0 .And. mmat==3) Then
      Do il = max(ip-2, im+1), min(ip+2, n)
        If (k(il,3)==im) nm = nm + 1
      End Do
      If (nm/=2 .Or. kfam<=100 .Or. mod(kfam,10)/=1 .Or. mod(kfam/1000,10)/=0) nm = 0
    Else If (im/=0 .And. mmat==46) Then
      msgn = isign(1, k(im,2)*k(ip,2))
      If (kfam>100 .And. mod(kfam/1000,10)==0) msgn = msgn*(-1)**mod(kfam/100, 10)
    End If
  End If

!...Kinematics of one-particle decays.
  If (nd==1) Then
    Do j = 1, 4
      p(n+1, j) = p(ip, j)
    End Do
    Goto 510
  End If

!...Calculate maximum weight ND-particle decay.
  pv(nd, 5) = p(n+nd, 5)
  If (nd>=3) Then
    wtmax = 1./wtcor(nd-2)
    pmax = pv(1, 5) - ps + p(n+nd, 5)
    pmin = 0.
    Do il = nd - 1, 1, -1
      pmax = pmax + p(n+il, 5)
      pmin = pmin + p(n+il+1, 5)
      wtmax = wtmax*pawt(pmax, pmin, p(n+il,5))
    End Do
  End If

!...Find virtual gamma mass in Dalitz decay.
  310 If (nd==2) Then
  Else If (mmat==2) Then
    pmes = 4.*pmas(11, 1)**2
    pmrho2 = pmas(131, 1)**2
    pgrho2 = pmas(131, 2)**2
    320 pmst = pmes*(p(ip,5)**2/pmes)**rlu(0)
    wt = (1+0.5*pmes/pmst)*sqrt(max(0.,1.-pmes/pmst))*(1.-pmst/p(ip,5)**2)**3*(1.+pgrho2/pmrho2)/((1.-pmst/pmrho2)**2+pgrho2/pmrho2)
    If (wt<rlu(0)) Goto 320
    pv(2, 5) = max(2.00001*pmas(11,1), sqrt(pmst))

!...M-generator gives weight. If rejected, try again.
  Else
    330 rord(1) = 1.
    Do il1 = 2, nd - 1
      rsav = rlu(0)
      Do il2 = il1 - 1, 1, -1
        If (rsav<=rord(il2)) Goto 350
        rord(il2+1) = rord(il2)
      End Do
      350 rord(il2+1) = rsav
    End Do
    rord(nd) = 0.
    wt = 1.
    Do il = nd - 1, 1, -1
      pv(il, 5) = pv(il+1, 5) + p(n+il, 5) + (rord(il)-rord(il+1))*(pv(1,5)-ps)
      wt = wt*pawt(pv(il,5), pv(il+1,5), p(n+il,5))
    End Do
    If (wt<rlu(0)*wtmax) Goto 330
  End If

!...Perform two-particle decays in respective CM frame.
  370 Do il = 1, nd - 1
    pa = pawt(pv(il,5), pv(il+1,5), p(n+il,5))
    ue(3) = 2.*rlu(0) - 1.
    phi = paru(2)*rlu(0)
    ue(1) = sqrt(1.-ue(3)**2)*cos(phi)
    ue(2) = sqrt(1.-ue(3)**2)*sin(phi)
    Do j = 1, 3
      p(n+il, j) = pa*ue(j)
      pv(il+1, j) = -pa*ue(j)
    End Do
    p(n+il, 4) = sqrt(pa**2+p(n+il,5)**2)
    pv(il+1, 4) = sqrt(pa**2+pv(il+1,5)**2)
  End Do

!...Lorentz transform decay products to lab frame.
  Do j = 1, 4
    p(n+nd, j) = pv(nd, j)
  End Do
  Do il = nd - 1, 1, -1
    Do j = 1, 3
      be(j) = pv(il, j)/pv(il, 4)
    End Do
    ga = pv(il, 4)/pv(il, 5)
    Do i = n + il, n + nd
      bep = be(1)*p(i, 1) + be(2)*p(i, 2) + be(3)*p(i, 3)
      Do j = 1, 3
        p(i, j) = p(i, j) + ga*(ga*bep/(1.+ga)+p(i,4))*be(j)
      End Do
      p(i, 4) = ga*(p(i,4)+bep)
    End Do
  End Do

!...Matrix elements for omega and phi decays.
  If (mmat==1) Then
    wt = (p(n+1,5)*p(n+2,5)*p(n+3,5))**2 - (p(n+1,5)*four(n+2,n+3))**2 - (p(n+2,5)*four(n+1,n+3))**2 - (p(n+3,5)*four(n+1,n+2))**2 + 2.*four(n+1, n+2)*four(n+1, n+3)*four(n+2, n+3)
    If (max(wt*wtcor(9)/p(ip,5)**6,0.001)<rlu(0)) Goto 310

!...Matrix elements for pi0 or eta Dalitz decay to gamma e+ e-.
  Else If (mmat==2) Then
    four12 = four(n+1, n+2)
    four13 = four(n+1, n+3)
    four23 = 0.5*pmst - 0.25*pmes
    wt = (pmst-0.5*pmes)*(four12**2+four13**2) + pmes*(four12*four13+four12**2+four13**2)
    If (wt<rlu(0)*0.25*pmst*(p(ip,5)**2-pmst)**2) Goto 370

!...Matrix element for S0 -> S1 + V1 -> S1 + S2 + S3 (S scalar,
!...V vector), of form cos**2(theta02) in V1 rest frame.
  Else If (mmat==3 .And. nm==2) Then
    If ((p(ip,5)**2*four(im,n+1)-four(ip,im)*four(ip,n+1))**2<=rlu(0)*(four(ip,im)**2-(p(ip,5)*p(im,5))**2)*(four(ip,n+1)**2-(p(ip,5)*p(n+1,5))**2)) Goto 370

!...Matrix element for "onium" -> g + g + g or gamma + g + g.
  Else If (mmat==4) Then
    hx1 = 2.*four(ip, n+1)/p(ip, 5)**2
    hx2 = 2.*four(ip, n+2)/p(ip, 5)**2
    hx3 = 2.*four(ip, n+3)/p(ip, 5)**2
    wt = ((1.-hx1)/(hx2*hx3))**2 + ((1.-hx2)/(hx1*hx3))**2 + ((1.-hx3)/(hx1*hx2))**2
    If (wt<2.*rlu(0)) Goto 310
    If (k(ip+1,2)==22 .And. (1.-hx1)*p(ip,5)**2<4.*parj(32)**2) Goto 310

!...Effective matrix element for nu spectrum in tau -> nu + hadrons.
  Else If (mmat==41) Then
    hx1 = 2.*four(ip, n+1)/p(ip, 5)**2
    If (8.*hx1*(3.-2.*hx1)/9.<rlu(0)) Goto 310

!...Matrix elements for weak decays (only semileptonic for c and b)
  Else If (mmat>=42 .And. mmat<=44 .And. nd==3) Then
    If (mbst==0) wt = four(ip, n+1)*four(n+2, n+3)
    If (mbst==1) wt = p(ip, 5)*p(n+1, 4)*four(n+2, n+3)
    If (wt<rlu(0)*p(ip,5)*pv(1,5)**3/wtcor(10)) Goto 310
  Else If (mmat>=42 .And. mmat<=44) Then
    Do j = 1, 4
      p(n+np+1, j) = 0.
      Do is = n + 3, n + np
        p(n+np+1, j) = p(n+np+1, j) + p(is, j)
      End Do
    End Do
    If (mbst==0) wt = four(ip, n+1)*four(n+2, n+np+1)
    If (mbst==1) wt = p(ip, 5)*p(n+1, 4)*four(n+2, n+np+1)
    If (wt<rlu(0)*p(ip,5)*pv(1,5)**3/wtcor(10)) Goto 310

!...Angular distribution in W decay.
  Else If (mmat==46 .And. msgn/=0) Then
    If (msgn>0) wt = four(im, n+1)*four(n+2, ip+1)
    If (msgn<0) wt = four(im, n+2)*four(n+1, ip+1)
    If (wt<rlu(0)*p(im,5)**4/wtcor(10)) Goto 370
  End If

!...Scale back energy and reattach spectator.
  If (mrem==1) Then
    Do j = 1, 5
      pv(1, j) = pv(1, j)/(1.-pqt)
    End Do
    nd = nd + 1
    mrem = 0
  End If

!...Low invariant mass for system with spectator quark gives particle,
!...not two jets. Readjust momenta accordingly.
  If ((mmat==31 .Or. mmat==45) .And. nd==3) Then
    mstj(93) = 1
    pm2 = ulmass(k(n+2,2))
    mstj(93) = 1
    pm3 = ulmass(k(n+3,2))
    If (p(n+2,5)**2+p(n+3,5)**2+2.*four(n+2,n+3)>=(parj(32)+pm2+pm3)**2) Goto 510
    k(n+2, 1) = 1
    kftemp = k(n+2, 2)
    Call lukfdi(kftemp, k(n+3,2), kfldmp, k(n+2,2))
    If (k(n+2,2)==0) Goto 150
    p(n+2, 5) = ulmass(k(n+2,2))
    ps = p(n+1, 5) + p(n+2, 5)
    pv(2, 5) = p(n+2, 5)
    mmat = 0
    nd = 2
    Goto 370
  Else If (mmat==44) Then
    mstj(93) = 1
    pm3 = ulmass(k(n+3,2))
    mstj(93) = 1
    pm4 = ulmass(k(n+4,2))
    If (p(n+3,5)**2+p(n+4,5)**2+2.*four(n+3,n+4)>=(parj(32)+pm3+pm4)**2) Goto 480
    k(n+3, 1) = 1
    kftemp = k(n+3, 2)
    Call lukfdi(kftemp, k(n+4,2), kfldmp, k(n+3,2))
    If (k(n+3,2)==0) Goto 150
    p(n+3, 5) = ulmass(k(n+3,2))
    Do j = 1, 3
      p(n+3, j) = p(n+3, j) + p(n+4, j)
    End Do
    p(n+3, 4) = sqrt(p(n+3,1)**2+p(n+3,2)**2+p(n+3,3)**2+p(n+3,5)**2)
    ha = p(n+1, 4)**2 - p(n+2, 4)**2
    hb = ha - (p(n+1,5)**2-p(n+2,5)**2)
    hc = (p(n+1,1)-p(n+2,1))**2 + (p(n+1,2)-p(n+2,2))**2 + (p(n+1,3)-p(n+2,3))**2
    hd = (pv(1,4)-p(n+3,4))**2
    he = ha**2 - 2.*hd*(p(n+1,4)**2+p(n+2,4)**2) + hd**2
    hf = hd*hc - hb**2
    hg = hd*hc - ha*hb
    hh = (sqrt(hg**2+he*hf)-hg)/(2.*hf)
    Do j = 1, 3
      pcor = hh*(p(n+1,j)-p(n+2,j))
      p(n+1, j) = p(n+1, j) + pcor
      p(n+2, j) = p(n+2, j) - pcor
    End Do
    p(n+1, 4) = sqrt(p(n+1,1)**2+p(n+1,2)**2+p(n+1,3)**2+p(n+1,5)**2)
    p(n+2, 4) = sqrt(p(n+2,1)**2+p(n+2,2)**2+p(n+2,3)**2+p(n+2,5)**2)
    nd = nd - 1
  End If

!...Check invariant mass of W jets. May give one particle or start over.
  480 If (mmat>=42 .And. mmat<=44 .And. iabs(k(n+1,2))<10) Then
    pmr = sqrt(max(0.,p(n+1,5)**2+p(n+2,5)**2+2.*four(n+1,n+2)))
    mstj(93) = 1
    pm1 = ulmass(k(n+1,2))
    mstj(93) = 1
    pm2 = ulmass(k(n+2,2))
    If (pmr>parj(32)+pm1+pm2) Goto 490
    kfldum = int(1.5+rlu(0))
    Call lukfdi(k(n+1,2), -isign(kfldum,k(n+1,2)), kfldmp, kf1)
    Call lukfdi(k(n+2,2), -isign(kfldum,k(n+2,2)), kfldmp, kf2)
    If (kf1==0 .Or. kf2==0) Goto 150
    psm = ulmass(kf1) + ulmass(kf2)
    If (mmat==42 .And. pmr>parj(64)+psm) Goto 490
    If (mmat>=43 .And. pmr>0.2*parj(32)+psm) Goto 490
    If (nd==4 .Or. kfa==15) Goto 150
    k(n+1, 1) = 1
    kftemp = k(n+1, 2)
    Call lukfdi(kftemp, k(n+2,2), kfldmp, k(n+1,2))
    If (k(n+1,2)==0) Goto 150
    p(n+1, 5) = ulmass(k(n+1,2))
    k(n+2, 2) = k(n+3, 2)
    p(n+2, 5) = p(n+3, 5)
    ps = p(n+1, 5) + p(n+2, 5)
    pv(2, 5) = p(n+3, 5)
    mmat = 0
    nd = 2
    Goto 370
  End If

!...Phase space decay of partons from W decay.
  490 If (mmat==42 .And. iabs(k(n+1,2))<10) Then
    kflo(1) = k(n+1, 2)
    kflo(2) = k(n+2, 2)
    k(n+1, 1) = k(n+3, 1)
    k(n+1, 2) = k(n+3, 2)
    Do j = 1, 5
      pv(1, j) = p(n+1, j) + p(n+2, j)
      p(n+1, j) = p(n+3, j)
    End Do
    pv(1, 5) = pmr
    n = n + 1
    np = 0
    nq = 2
    ps = 0.
    mstj(93) = 2
    psq = ulmass(kflo(1))
    mstj(93) = 2
    psq = psq + ulmass(kflo(2))
    mmat = 11
    Goto 180
  End If

!...Boost back for rapidly moving particle.
  510 n = n + nd
  If (mbst==1) Then
    Do j = 1, 3
      be(j) = p(ip, j)/p(ip, 4)
    End Do
    ga = p(ip, 4)/p(ip, 5)
    Do i = nsav + 1, n
      bep = be(1)*p(i, 1) + be(2)*p(i, 2) + be(3)*p(i, 3)
      Do j = 1, 3
        p(i, j) = p(i, j) + ga*(ga*bep/(1.+ga)+p(i,4))*be(j)
      End Do
      p(i, 4) = ga*(p(i,4)+bep)
    End Do
  End If

!...Fill in position of decay vertex.
  Do i = nsav + 1, n
    Do j = 1, 4
      v(i, j) = vdcy(j)
    End Do
    v(i, 5) = 0.
  End Do

!...Set up for parton shower evolution from jets.
  If (mstj(23)>=1 .And. mmat==4 .And. k(nsav+1,2)==21) Then
    k(nsav+1, 1) = 3
    k(nsav+2, 1) = 3
    k(nsav+3, 1) = 3
    k(nsav+1, 4) = mstu(5)*(nsav+2)
    k(nsav+1, 5) = mstu(5)*(nsav+3)
    k(nsav+2, 4) = mstu(5)*(nsav+3)
    k(nsav+2, 5) = mstu(5)*(nsav+1)
    k(nsav+3, 4) = mstu(5)*(nsav+1)
    k(nsav+3, 5) = mstu(5)*(nsav+2)
    mstj(92) = -(nsav+1)
  Else If (mstj(23)>=1 .And. mmat==4) Then
    k(nsav+2, 1) = 3
    k(nsav+3, 1) = 3
    k(nsav+2, 4) = mstu(5)*(nsav+3)
    k(nsav+2, 5) = mstu(5)*(nsav+3)
    k(nsav+3, 4) = mstu(5)*(nsav+2)
    k(nsav+3, 5) = mstu(5)*(nsav+2)
    mstj(92) = nsav + 2
  Else If (mstj(23)>=1 .And. (mmat==32 .Or. mmat==44 .Or. mmat==46) .And. iabs(k(nsav+1,2))<=10 .And. iabs(k(nsav+2,2))<=10) Then
    k(nsav+1, 1) = 3
    k(nsav+2, 1) = 3
    k(nsav+1, 4) = mstu(5)*(nsav+2)
    k(nsav+1, 5) = mstu(5)*(nsav+2)
    k(nsav+2, 4) = mstu(5)*(nsav+1)
    k(nsav+2, 5) = mstu(5)*(nsav+1)
    mstj(92) = nsav + 1
  Else If (mstj(23)>=1 .And. mmat==33 .And. iabs(k(nsav+2,2))==21) Then
    k(nsav+1, 1) = 3
    k(nsav+2, 1) = 3
    k(nsav+3, 1) = 3
    kcp = lucomp(k(nsav+1,2))
    kqp = kchg(kcp, 2)*isign(1, k(nsav+1,2))
    jcon = 4
    If (kqp<0) jcon = 5
    k(nsav+1, jcon) = mstu(5)*(nsav+2)
    k(nsav+2, 9-jcon) = mstu(5)*(nsav+1)
    k(nsav+2, jcon) = mstu(5)*(nsav+3)
    k(nsav+3, 9-jcon) = mstu(5)*(nsav+2)
    mstj(92) = nsav + 1
  Else If (mstj(23)>=1 .And. mmat==33) Then
    k(nsav+1, 1) = 3
    k(nsav+3, 1) = 3
    k(nsav+1, 4) = mstu(5)*(nsav+3)
    k(nsav+1, 5) = mstu(5)*(nsav+3)
    k(nsav+3, 4) = mstu(5)*(nsav+1)
    k(nsav+3, 5) = mstu(5)*(nsav+1)
    mstj(92) = nsav + 1
  End If

!...Mark decayed particle.
  If (k(ip,1)==5) k(ip, 1) = 15
  If (k(ip,1)<=10) k(ip, 1) = 11
  k(ip, 4) = nsav + 1
  k(ip, 5) = n

  Return
End Subroutine ludecy

!*********************************************************************

Subroutine lukfdi(kfl1, kfl2, kfl3, kf)

!...Purpose: to generate a new flavour pair and combine off a hadron.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/

!...Default flavour values. Input consistency checks.
  kf1a = iabs(kfl1)
  kf2a = iabs(kfl2)
  kfl3 = 0
  kf = 0
  If (kf1a==0) Return
  If (kf2a/=0) Then
    If (kf1a<=10 .And. kf2a<=10 .And. kfl1*kfl2>0) Return
    If (kf1a>10 .And. kf2a>10) Return
    If ((kf1a>10 .Or. kf2a>10) .And. kfl1*kfl2<0) Return
  End If

!...Check if tabulated flavour probabilities are to be used.
  If (mstj(15)==1) Then
    ktab1 = -1
    If (kf1a>=1 .And. kf1a<=6) ktab1 = kf1a
    kfl1a = mod(kf1a/1000, 10)
    kfl1b = mod(kf1a/100, 10)
    kfl1s = mod(kf1a, 10)
    If (kfl1a>=1 .And. kfl1a<=4 .And. kfl1b>=1 .And. kfl1b<=4) ktab1 = 6 + kfl1a*(kfl1a-2) + 2*kfl1b + (kfl1s-1)/2
    If (kfl1a>=1 .And. kfl1a<=4 .And. kfl1a==kfl1b) ktab1 = ktab1 - 1
    If (kf1a>=1 .And. kf1a<=6) kfl1a = kf1a
    ktab2 = 0
    If (kf2a/=0) Then
      ktab2 = -1
      If (kf2a>=1 .And. kf2a<=6) ktab2 = kf2a
      kfl2a = mod(kf2a/1000, 10)
      kfl2b = mod(kf2a/100, 10)
      kfl2s = mod(kf2a, 10)
      If (kfl2a>=1 .And. kfl2a<=4 .And. kfl2b>=1 .And. kfl2b<=4) ktab2 = 6 + kfl2a*(kfl2a-2) + 2*kfl2b + (kfl2s-1)/2
      If (kfl2a>=1 .And. kfl2a<=4 .And. kfl2a==kfl2b) ktab2 = ktab2 - 1
    End If
    If (ktab1>=0 .And. ktab2>=0) Goto 140
  End If

!...Parameters and breaking diquark parameter combinations.
  100 par2 = parj(2)
  par3 = parj(3)
  par4 = 3.*parj(4)
  If (mstj(12)>=2) Then
    par3m = sqrt(parj(3))
    par4m = 1./(3.*sqrt(parj(4)))
    pardm = parj(7)/(parj(7)+par3m*parj(6))
    pars0 = parj(5)*(2.+(1.+par2*par3m*parj(7))*(1.+par4m))
    pars1 = parj(7)*pars0/(2.*par3m) + parj(5)*(parj(6)*(1.+par4m)+par2*par3m*parj(6)*parj(7))
    pars2 = parj(5)*2.*parj(6)*parj(7)*(par2*parj(7)+(1.+par4m)/par3m)
    parsm = max(pars0, pars1, pars2)
    par4 = par4*(1.+parsm)/(1.+parsm/(3.*par4m))
  End If

!...Choice of whether to generate meson or baryon.
  mbary = 0
  kfda = 0
  If (kf1a<=10) Then
    If (kf2a==0 .And. mstj(12)>=1 .And. (1.+parj(1))*rlu(0)>1.) mbary = 1
    If (kf2a>10) mbary = 2
    If (kf2a>10 .And. kf2a<=10000) kfda = kf2a
  Else
    mbary = 2
    If (kf1a<=10000) kfda = kf1a
  End If

!...Possibility of process diquark -> meson + new diquark.
  If (kfda/=0 .And. mstj(12)>=2) Then
    kflda = mod(kfda/1000, 10)
    kfldb = mod(kfda/100, 10)
    kflds = mod(kfda, 10)
    wtdq = pars0
    If (max(kflda,kfldb)==3) wtdq = pars1
    If (min(kflda,kfldb)==3) wtdq = pars2
    If (kflds==1) wtdq = wtdq/(3.*par4m)
    If ((1.+wtdq)*rlu(0)>1.) mbary = -1
    If (mbary==-1 .And. kf2a/=0) Return
  End If

!...Flavour for meson, possibly with new flavour.
  If (mbary<=0) Then
    kfs = isign(1, kfl1)
    If (mbary==0) Then
      If (kf2a==0) kfl3 = isign(1+int((2.+par2)*rlu(0)), -kfl1)
      kfla = max(kf1a, kf2a+iabs(kfl3))
      kflb = min(kf1a, kf2a+iabs(kfl3))
      If (kfla/=kf1a) kfs = -kfs

!...Splitting of diquark into meson plus new diquark.
    Else
      kfl1a = mod(kf1a/1000, 10)
      kfl1b = mod(kf1a/100, 10)
      110 kfl1d = kfl1a + int(rlu(0)+0.5)*(kfl1b-kfl1a)
      kfl1e = kfl1a + kfl1b - kfl1d
      If ((kfl1d==3 .And. rlu(0)>pardm) .Or. (kfl1e==3 .And. rlu(0)<pardm)) Then
        kfl1d = kfl1a + kfl1b - kfl1d
        kfl1e = kfl1a + kfl1b - kfl1e
      End If
      kfl3a = 1 + int((2.+par2*par3m*parj(7))*rlu(0))
      If ((kfl1e/=kfl3a .And. rlu(0)>(1.+par4m)/max(2.,1.+par4m)) .Or. (kfl1e==kfl3a .And. rlu(0)>2./max(2.,1.+par4m))) Goto 110
      kflds = 3
      If (kfl1e/=kfl3a) kflds = 2*int(rlu(0)+1./(1.+par4m)) + 1
      kfl3 = isign(10000+1000*max(kfl1e,kfl3a)+100*min(kfl1e,kfl3a)+kflds, -kfl1)
      kfla = max(kfl1d, kfl3a)
      kflb = min(kfl1d, kfl3a)
      If (kfla/=kfl1d) kfs = -kfs
    End If

!...Form meson, with spin and flavour mixing for diagonal states.
    If (kfla<=2) kmul = int(parj(11)+rlu(0))
    If (kfla==3) kmul = int(parj(12)+rlu(0))
    If (kfla>=4) kmul = int(parj(13)+rlu(0))
    If (kmul==0 .And. parj(14)>0.) Then
      If (rlu(0)<parj(14)) kmul = 2
    Else If (kmul==1 .And. parj(15)+parj(16)+parj(17)>0.) Then
      rmul = rlu(0)
      If (rmul<parj(15)) kmul = 3
      If (kmul==1 .And. rmul<parj(15)+parj(16)) kmul = 4
      If (kmul==1 .And. rmul<parj(15)+parj(16)+parj(17)) kmul = 5
    End If
    kfls = 3
    If (kmul==0 .Or. kmul==3) kfls = 1
    If (kmul==5) kfls = 5
    If (kfla/=kflb) Then
      kf = (100*kfla+10*kflb+kfls)*kfs*(-1)**kfla
    Else
      rmix = rlu(0)
      imix = 2*kfla + 10*kmul
      If (kfla<=3) kf = 110*(1+int(rmix+parf(imix-1))+int(rmix+parf(imix))) + kfls
      If (kfla>=4) kf = 110*kfla + kfls
    End If
    If (kmul==2 .Or. kmul==3) kf = kf + isign(10000, kf)
    If (kmul==4) kf = kf + isign(20000, kf)

!...Generate diquark flavour.
  Else
    120 If (kf1a<=10 .And. kf2a==0) Then
      kfla = kf1a
      130 kflb = 1 + int((2.+par2*par3)*rlu(0))
      kflc = 1 + int((2.+par2*par3)*rlu(0))
      kflds = 1
      If (kflb>=kflc) kflds = 3
      If (kflds==1 .And. par4*rlu(0)>1.) Goto 130
      If (kflds==3 .And. par4<rlu(0)) Goto 130
      kfl3 = isign(1000*max(kflb,kflc)+100*min(kflb,kflc)+kflds, kfl1)

!...Take diquark flavour from input.
    Else If (kf1a<=10) Then
      kfla = kf1a
      kflb = mod(kf2a/1000, 10)
      kflc = mod(kf2a/100, 10)
      kflds = mod(kf2a, 10)

!...Generate (or take from input) quark to go with diquark.
    Else
      If (kf2a==0) kfl3 = isign(1+int((2.+par2)*rlu(0)), kfl1)
      kfla = kf2a + iabs(kfl3)
      kflb = mod(kf1a/1000, 10)
      kflc = mod(kf1a/100, 10)
      kflds = mod(kf1a, 10)
    End If

!...SU(6) factors for formation of baryon. Try again if fails.
    kbary = kflds
    If (kflds==3 .And. kflb/=kflc) kbary = 5
    If (kfla/=kflb .And. kfla/=kflc) kbary = kbary + 1
    wt = parf(60+kbary) + parj(18)*parf(70+kbary)
    If (mbary==1 .And. mstj(12)>=2) Then
      wtdq = pars0
      If (max(kflb,kflc)==3) wtdq = pars1
      If (min(kflb,kflc)==3) wtdq = pars2
      If (kflds==1) wtdq = wtdq/(3.*par4m)
      If (kflds==1) wt = wt*(1.+wtdq)/(1.+parsm/(3.*par4m))
      If (kflds==3) wt = wt*(1.+wtdq)/(1.+parsm)
    End If
    If (kf2a==0 .And. wt<rlu(0)) Goto 120

!...Form baryon. Distinguish Lambda- and Sigmalike baryons.
    kfld = max(kfla, kflb, kflc)
    kflf = min(kfla, kflb, kflc)
    kfle = kfla + kflb + kflc - kfld - kflf
    kfls = 2
    If ((parf(60+kbary)+parj(18)*parf(70+kbary))*rlu(0)>parf(60+kbary)) kfls = 4
    kfll = 0
    If (kfls==2 .And. kfld>kfle .And. kfle>kflf) Then
      If (kflds==1 .And. kfla==kfld) kfll = 1
      If (kflds==1 .And. kfla/=kfld) kfll = int(0.25+rlu(0))
      If (kflds==3 .And. kfla/=kfld) kfll = int(0.75+rlu(0))
    End If
    If (kfll==0) kf = isign(1000*kfld+100*kfle+10*kflf+kfls, kfl1)
    If (kfll==1) kf = isign(1000*kfld+100*kflf+10*kfle+kfls, kfl1)
  End If
  Return

!...Use tabulated probabilities to select new flavour and hadron.
  140 If (ktab2==0 .And. mstj(12)<=0) Then
    kt3l = 1
    kt3u = 6
  Else If (ktab2==0 .And. ktab1>=7 .And. mstj(12)<=1) Then
    kt3l = 1
    kt3u = 6
  Else If (ktab2==0) Then
    kt3l = 1
    kt3u = 22
  Else
    kt3l = ktab2
    kt3u = ktab2
  End If
  rfl = 0.
  Do kts = 0, 2
    Do kt3 = kt3l, kt3u
      rfl = rfl + parf(120+80*ktab1+25*kts+kt3)
    End Do
  End Do
  rfl = rlu(0)*rfl
  Do kts = 0, 2
    ktabs = kts
    Do kt3 = kt3l, kt3u
      ktab3 = kt3
      rfl = rfl - parf(120+80*ktab1+25*kts+kt3)
      If (rfl<=0.) Goto 170
    End Do
  End Do
  170 Continue

!...Reconstruct flavour of produced quark/diquark.
  If (ktab3<=6) Then
    kfl3a = ktab3
    kfl3b = 0
    kfl3 = isign(kfl3a, kfl1*(2*ktab1-13))
  Else
    kfl3a = 1
    If (ktab3>=8) kfl3a = 2
    If (ktab3>=11) kfl3a = 3
    If (ktab3>=16) kfl3a = 4
    kfl3b = (ktab3-6-kfl3a*(kfl3a-2))/2
    kfl3 = 1000*kfl3a + 100*kfl3b + 1
    If (kfl3a==kfl3b .Or. ktab3/=6+kfl3a*(kfl3a-2)+2*kfl3b) kfl3 = kfl3 + 2
    kfl3 = isign(kfl3, kfl1*(13-2*ktab1))
  End If

!...Reconstruct meson code.
  If (kfl3a==kfl1a .And. kfl3b==kfl1b .And. (kfl3a<=3 .Or. kfl3b/=0)) Then
    rfl = rlu(0)*(parf(143+80*ktab1+25*ktabs)+parf(144+80*ktab1+25*ktabs)+parf(145+80*ktab1+25*ktabs))
    kf = 110 + 2*ktabs + 1
    If (rfl>parf(143+80*ktab1+25*ktabs)) kf = 220 + 2*ktabs + 1
    If (rfl>parf(143+80*ktab1+25*ktabs)+parf(144+80*ktab1+25*ktabs)) kf = 330 + 2*ktabs + 1
  Else If (ktab1<=6 .And. ktab3<=6) Then
    kfla = max(ktab1, ktab3)
    kflb = min(ktab1, ktab3)
    kfs = isign(1, kfl1)
    If (kfla/=kf1a) kfs = -kfs
    kf = (100*kfla+10*kflb+2*ktabs+1)*kfs*(-1)**kfla
  Else If (ktab1>=7 .And. ktab3>=7) Then
    kfs = isign(1, kfl1)
    If (kfl1a==kfl3a) Then
      kfla = max(kfl1b, kfl3b)
      kflb = min(kfl1b, kfl3b)
      If (kfla/=kfl1b) kfs = -kfs
    Else If (kfl1a==kfl3b) Then
      kfla = kfl3a
      kflb = kfl1b
      kfs = -kfs
    Else If (kfl1b==kfl3a) Then
      kfla = kfl1a
      kflb = kfl3b
    Else If (kfl1b==kfl3b) Then
      kfla = max(kfl1a, kfl3a)
      kflb = min(kfl1a, kfl3a)
      If (kfla/=kfl1a) kfs = -kfs
    Else
      Call luerrm(2, '(LUKFDI:) no matching flavours for qq -> qq')
      Goto 100
    End If
    kf = (100*kfla+10*kflb+2*ktabs+1)*kfs*(-1)**kfla

!...Reconstruct baryon code.
  Else
    If (ktab1>=7) Then
      kfla = kfl3a
      kflb = kfl1a
      kflc = kfl1b
    Else
      kfla = kfl1a
      kflb = kfl3a
      kflc = kfl3b
    End If
    kfld = max(kfla, kflb, kflc)
    kflf = min(kfla, kflb, kflc)
    kfle = kfla + kflb + kflc - kfld - kflf
    If (ktabs==0) kf = isign(1000*kfld+100*kflf+10*kfle+2, kfl1)
    If (ktabs>=1) kf = isign(1000*kfld+100*kfle+10*kflf+2*ktabs, kfl1)
  End If

!...Check that constructed flavour code is an allowed one.
  If (kfl2/=0) kfl3 = 0
  kc = lucomp(kf)
  If (kc==0) Then
    Call luerrm(2, '(LUKFDI:) user-defined flavour probabilities '//'failed')
    Goto 100
  End If

  Return
End Subroutine lukfdi

!*********************************************************************

Subroutine luptdi(kfl, px, py)

!...Purpose: to generate transverse momentum according to a Gaussian.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/

!...Generate p_T and azimuthal angle, gives p_x and p_y.
  kfla = iabs(kfl)
  pt = parj(21)*sqrt(-log(max(1E-10,rlu(0))))
  If (mstj(91)==1) pt = parj(22)*pt
  If (kfla==0 .And. mstj(13)<=0) pt = 0.
  phi = paru(2)*rlu(0)
  px = pt*cos(phi)
  py = pt*sin(phi)

  Return
End Subroutine luptdi

!*********************************************************************

Subroutine luzdis(kfl1, kfl2, pr, z)

!...Purpose: to generate the longitudinal splitting variable z.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/

!...Check if heavy flavour fragmentation.
  kfla = iabs(kfl1)
  kflb = iabs(kfl2)
  kflh = kfla
  If (kfla>=10) kflh = mod(kfla/1000, 10)

!...Lund symmetric scaling function: determine parameters of shape.
  If (mstj(11)==1 .Or. (mstj(11)==3 .And. kflh<=3)) Then
    fa = parj(41)
    If (mstj(91)==1) fa = parj(43)
    If (kflb>=10) fa = fa + parj(45)
    fb = parj(42)*pr
    If (mstj(91)==1) fb = parj(44)*pr
    fc = 1.
    If (kfla>=10) fc = fc - parj(45)
    If (kflb>=10) fc = fc + parj(45)
    mc = 1
    If (abs(fc-1.)>0.01) mc = 2

!...Determine position of maximum. Special cases for a = 0 or a = c.
    If (fa<0.02) Then
      ma = 1
      zmax = 1.
      If (fc>fb) zmax = fb/fc
    Else If (abs(fc-fa)<0.01) Then
      ma = 2
      zmax = fb/(fb+fc)
    Else
      ma = 3
      zmax = 0.5*(fb+fc-sqrt((fb-fc)**2+4.*fa*fb))/(fc-fa)
      If (zmax>0.99 .And. fb>100.) zmax = 1. - fa/fb
    End If

!...Subdivide z range if distribution very peaked near endpoint.
    mmax = 2
    If (zmax<0.1) Then
      mmax = 1
      zdiv = 2.75*zmax
      If (mc==1) Then
        fint = 1. - log(zdiv)
      Else
        zdivc = zdiv**(1.-fc)
        fint = 1. + (1.-1./zdivc)/(fc-1.)
      End If
    Else If (zmax>0.85 .And. fb>1.) Then
      mmax = 3
      fscb = sqrt(4.+(fc/fb)**2)
      zdiv = fscb - 1./zmax - (fc/fb)*log(zmax*0.5*(fscb+fc/fb))
      If (ma>=2) zdiv = zdiv + (fa/fb)*log(1.-zmax)
      zdiv = min(zmax, max(0.,zdiv))
      fint = 1. + fb*(1.-zdiv)
    End If

!...Choice of z, preweighted for peaks at low or high z.
    100 z = rlu(0)
    fpre = 1.
    If (mmax==1) Then
      If (fint*rlu(0)<=1.) Then
        z = zdiv*z
      Else If (mc==1) Then
        z = zdiv**z
        fpre = zdiv/z
      Else
        z = 1./(zdivc+z*(1.-zdivc))**(1./(1.-fc))
        fpre = (zdiv/z)**fc
      End If
    Else If (mmax==3) Then
      If (fint*rlu(0)<=1.) Then
        z = zdiv + log(z)/fb
        fpre = exp(fb*(z-zdiv))
      Else
        z = zdiv + z*(1.-zdiv)
      End If
    End If

!...Weighting according to correct formula.
    If (z<=fb/(50.+fb) .Or. z>=1.) Goto 100
    fval = (zmax/z)**fc*exp(fb*(1./zmax-1./z))
    If (ma>=2) fval = ((1.-z)/(1.-zmax))**fa*fval
    If (fval<rlu(0)*fpre) Goto 100

!...Generate z according to Field-Feynman, SLAC, (1-z)**c OR z**c.
  Else
    fc = parj(50+max(1,kflh))
    If (mstj(91)==1) fc = parj(59)
    110 z = rlu(0)
    If (fc>=0. .And. fc<=1.) Then
      If (fc>rlu(0)) z = 1. - z**(1./3.)
    Else If (fc>-1.) Then
      If (-4.*fc*z*(1.-z)**2<rlu(0)*((1.-z)**2-fc*z)**2) Goto 110
    Else
      If (fc>0.) z = 1. - z**(1./fc)
      If (fc<0.) z = z**(-1./fc)
    End If
  End If

  Return
End Subroutine luzdis

!*********************************************************************

Subroutine lushow(ip1, ip2, qmax)

!...Purpose: to generate timelike parton showers from given partons.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension pmth(5, 40), ps(5), pma(4), pmsd(4), iep(4), ipa(4), kfla(4), kfld(4), kfl(4), itry(4), isi(4), isl(4), dp(4), dpt(5, 4)

!...Initialization of cutoff masses etc.
  If (mstj(41)<=0 .Or. (mstj(41)==1 .And. qmax<=parj(82)) .Or. qmax<=min(parj(82),parj(83)) .Or. mstj(41)>=3) Return
  pmth(1, 21) = ulmass(21)
  pmth(2, 21) = sqrt(pmth(1,21)**2+0.25*parj(82)**2)
  pmth(3, 21) = 2.*pmth(2, 21)
  pmth(4, 21) = pmth(3, 21)
  pmth(5, 21) = pmth(3, 21)
  pmth(1, 22) = ulmass(22)
  pmth(2, 22) = sqrt(pmth(1,22)**2+0.25*parj(83)**2)
  pmth(3, 22) = 2.*pmth(2, 22)
  pmth(4, 22) = pmth(3, 22)
  pmth(5, 22) = pmth(3, 22)
  pmqth1 = parj(82)
  If (mstj(41)==2) pmqth1 = min(parj(82), parj(83))
  pmqth2 = pmth(2, 21)
  If (mstj(41)==2) pmqth2 = min(pmth(2,21), pmth(2,22))
  Do if = 1, 8
    pmth(1, if) = ulmass(if)
    pmth(2, if) = sqrt(pmth(1,if)**2+0.25*pmqth1**2)
    pmth(3, if) = pmth(2, if) + pmqth2
    pmth(4, if) = sqrt(pmth(1,if)**2+0.25*parj(82)**2) + pmth(2, 21)
    pmth(5, if) = sqrt(pmth(1,if)**2+0.25*parj(83)**2) + pmth(2, 22)
  End Do
  pt2min = max(0.5*parj(82), 1.1*parj(81))**2
  alams = parj(81)**2
  alfm = log(pt2min/alams)

!...Store positions of shower initiating partons.
  m3jc = 0
  If (ip1>0 .And. ip1<=min(n,mstu(4)-mstu(32)) .And. ip2==0) Then
    npa = 1
    ipa(1) = ip1
  Else If (min(ip1,ip2)>0 .And. max(ip1,ip2)<=min(n,mstu(4)-mstu(32))) Then
    npa = 2
    ipa(1) = ip1
    ipa(2) = ip2
  Else If (ip1>0 .And. ip1<=min(n,mstu(4)-mstu(32)) .And. ip2<0 .And. ip2>=-3) Then
    npa = iabs(ip2)
    Do i = 1, npa
      ipa(i) = ip1 + i - 1
    End Do
  Else
    Call luerrm(12, '(LUSHOW:) failed to reconstruct showering system')
    If (mstu(21)>=1) Return
  End If

!...Check on phase space available for emission.
  irej = 0
  Do j = 1, 5
    ps(j) = 0.
  End Do
  pm = 0.
  Do i = 1, npa
    kfla(i) = iabs(k(ipa(i),2))
    pma(i) = p(ipa(i), 5)
    If (kfla(i)/=0 .And. (kfla(i)<=8 .Or. kfla(i)==21)) pma(i) = pmth(3, kfla(i))
    pm = pm + pma(i)
    If (kfla(i)==0 .Or. (kfla(i)>8 .And. kfla(i)/=21) .Or. pma(i)>qmax) irej = irej + 1
    Do j = 1, 4
      ps(j) = ps(j) + p(ipa(i), j)
    End Do
  End Do
  If (irej==npa) Return
  ps(5) = sqrt(max(0.,ps(4)**2-ps(1)**2-ps(2)**2-ps(3)**2))
  If (npa==1) ps(5) = ps(4)
  If (ps(5)<=pm+pmqth1) Return
  If (npa==2 .And. mstj(47)>=1) Then
    If (kfla(1)>=1 .And. kfla(1)<=8 .And. kfla(2)>=1 .And. kfla(2)<=8) m3jc = 1
    If (mstj(47)>=2) m3jc = 1
  End If

!...Define imagined single initiator of shower for parton system.
  ns = n
  If (n>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSHOW:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If
  If (npa>=2) Then
    k(n+1, 1) = 11
    k(n+1, 2) = 21
    k(n+1, 3) = 0
    k(n+1, 4) = 0
    k(n+1, 5) = 0
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = 0.
    p(n+1, 4) = ps(5)
    p(n+1, 5) = ps(5)
    v(n+1, 5) = ps(5)**2
    n = n + 1
  End If

!...Loop over partons that may branch.
  nep = npa
  im = ns
  If (npa==1) im = ns - 1
  140 im = im + 1
  If (n>ns) Then
    If (im>n) Goto 380
    kflm = iabs(k(im,2))
    If (kflm==0 .Or. (kflm>8 .And. kflm/=21)) Goto 140
    If (p(im,5)<pmth(2,kflm)) Goto 140
    igm = k(im, 3)
  Else
    igm = -1
  End If
  If (n+nep>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSHOW:) no more memory left in LUJETS')
    If (mstu(21)>=1) Return
  End If

!...Position of aunt (sister to branching parton).
!...Origin and flavour of daughters.
  iau = 0
  If (igm>0) Then
    If (k(im-1,3)==igm) iau = im - 1
    If (n>=im+1 .And. k(im+1,3)==igm) iau = im + 1
  End If
  If (igm>=0) Then
    k(im, 4) = n + 1
    Do i = 1, nep
      k(n+i, 3) = im
    End Do
  Else
    k(n+1, 3) = ipa(1)
  End If
  If (igm<=0) Then
    Do i = 1, nep
      k(n+i, 2) = k(ipa(i), 2)
    End Do
  Else If (kflm/=21) Then
    k(n+1, 2) = k(im, 2)
    k(n+2, 2) = k(im, 5)
  Else If (k(im,5)==21) Then
    k(n+1, 2) = 21
    k(n+2, 2) = 21
  Else
    k(n+1, 2) = k(im, 5)
    k(n+2, 2) = -k(im, 5)
  End If

!...Reset flags on daughers and tries made.
  Do ip = 1, nep
    k(n+ip, 1) = 3
    k(n+ip, 4) = 0
    k(n+ip, 5) = 0
    kfld(ip) = iabs(k(n+ip,2))
    itry(ip) = 0
    isl(ip) = 0
    isi(ip) = 0
    If (kfld(ip)>0 .And. (kfld(ip)<=8 .Or. kfld(ip)==21)) isi(ip) = 1
  End Do
  islm = 0

!...Maximum virtuality of daughters.
  If (igm<=0) Then
    Do i = 1, npa
      If (npa>=3) p(n+i, 4) = (ps(4)*p(ipa(i),4)-ps(1)*p(ipa(i),1)-ps(2)*p(ipa(i),2)-ps(3)*p(ipa(i),3))/ps(5)
      p(n+i, 5) = min(qmax, ps(5))
      If (npa>=3) p(n+i, 5) = min(p(n+i,5), p(n+i,4))
      If (isi(i)==0) p(n+i, 5) = p(ipa(i), 5)
    End Do
  Else
    If (mstj(43)<=2) pem = v(im, 2)
    If (mstj(43)>=3) pem = p(im, 4)
    p(n+1, 5) = min(p(im,5), v(im,1)*pem)
    p(n+2, 5) = min(p(im,5), (1.-v(im,1))*pem)
    If (k(n+2,2)==22) p(n+2, 5) = pmth(1, 22)
  End If
  Do i = 1, nep
    pmsd(i) = p(n+i, 5)
    If (isi(i)==1) Then
      If (p(n+i,5)<=pmth(3,kfld(i))) p(n+i, 5) = pmth(1, kfld(i))
    End If
    v(n+i, 5) = p(n+i, 5)**2
  End Do

!...Choose one of the daughters for evolution.
  200 inum = 0
  If (nep==1) inum = 1
  Do i = 1, nep
    If (inum==0 .And. isl(i)==1) inum = i
  End Do
  Do i = 1, nep
    If (inum==0 .And. itry(i)==0 .And. isi(i)==1) Then
      If (p(n+i,5)>=pmth(2,kfld(i))) inum = i
    End If
  End Do
  If (inum==0) Then
    rmax = 0.
    Do i = 1, nep
      If (isi(i)==1 .And. pmsd(i)>=pmqth2) Then
        rpm = p(n+i, 5)/pmsd(i)
        If (rpm>rmax .And. p(n+i,5)>=pmth(2,kfld(i))) Then
          rmax = rpm
          inum = i
        End If
      End If
    End Do
  End If

!...Store information on choice of evolving daughter.
  inum = max(1, inum)
  iep(1) = n + inum
  Do i = 2, nep
    iep(i) = iep(i-1) + 1
    If (iep(i)>n+nep) iep(i) = n + 1
  End Do
  Do i = 1, nep
    kfl(i) = iabs(k(iep(i),2))
  End Do
  itry(inum) = itry(inum) + 1
  If (itry(inum)>200) Then
    Call luerrm(14, '(LUSHOW:) caught in infinite loop')
    If (mstu(21)>=1) Return
  End If
  z = 0.5
  If (kfl(1)==0 .Or. (kfl(1)>8 .And. kfl(1)/=21)) Goto 300
  If (p(iep(1),5)<pmth(2,kfl(1))) Goto 300

!...Calculate allowed z range.
  If (nep==1) Then
    pmed = ps(4)
  Else If (igm==0 .Or. mstj(43)<=2) Then
    pmed = p(im, 5)
  Else
    If (inum==1) pmed = v(im, 1)*pem
    If (inum==2) pmed = (1.-v(im,1))*pem
  End If
  If (mod(mstj(43),2)==1) Then
    zc = pmth(2, 21)/pmed
    zce = pmth(2, 22)/pmed
  Else
    zc = 0.5*(1.-sqrt(max(0.,1.-(2.*pmth(2,21)/pmed)**2)))
    If (zc<1E-4) zc = (pmth(2,21)/pmed)**2
    zce = 0.5*(1.-sqrt(max(0.,1.-(2.*pmth(2,22)/pmed)**2)))
    If (zce<1E-4) zce = (pmth(2,22)/pmed)**2
  End If
  zc = min(zc, 0.491)
  zce = min(zce, 0.491)
  If ((mstj(41)==1 .And. zc>0.49) .Or. (mstj(41)==2 .And. min(zc,zce)>0.49)) Then
    p(iep(1), 5) = pmth(1, kfl(1))
    v(iep(1), 5) = p(iep(1), 5)**2
    Goto 300
  End If

!...Integral of Altarelli-Parisi z kernel for QCD.
  If (mstj(49)==0 .And. kfl(1)==21) Then
    fbr = 6.*log((1.-zc)/zc) + mstj(45)*(0.5-zc)
  Else If (mstj(49)==0) Then
    fbr = (8./3.)*log((1.-zc)/zc)

!...Integral of Altarelli-Parisi z kernel for scalar gluon.
  Else If (mstj(49)==1 .And. kfl(1)==21) Then
    fbr = (parj(87)+mstj(45)*parj(88))*(1.-2.*zc)
  Else If (mstj(49)==1) Then
    fbr = (1.-2.*zc)/3.
    If (igm==0 .And. m3jc==1) fbr = 4.*fbr

!...Integral of Altarelli-Parisi z kernel for Abelian vector gluon.
  Else If (kfl(1)==21) Then
    fbr = 6.*mstj(45)*(0.5-zc)
  Else
    fbr = 2.*log((1.-zc)/zc)
  End If

!...Integral of Altarelli-Parisi kernel for photon emission.
  If (mstj(41)==2 .And. kfl(1)>=1 .And. kfl(1)<=8) fbre = (kchg(kfl(1),1)/3.)**2*2.*log((1.-zce)/zce)

!...Inner veto algorithm starts. Find maximum mass for evolution.
  260 pms = v(iep(1), 5)
  If (igm>=0) Then
    pm2 = 0.
    Do i = 2, nep
      pm = p(iep(i), 5)
      If (kfl(i)>0 .And. (kfl(i)<=8 .Or. kfl(i)==21)) pm = pmth(2, kfl(i))
      pm2 = pm2 + pm
    End Do
    pms = min(pms, (p(im,5)-pm2)**2)
  End If

!...Select mass for daughter in QCD evolution.
  b0 = 27./6.
  Do if = 4, mstj(45)
    If (pms>4.*pmth(2,if)**2) b0 = (33.-2.*if)/6.
  End Do
  If (mstj(44)<=0) Then
    pmsqcd = pms*exp(max(-100.,log(rlu(0))*paru(2)/(paru(111)*fbr)))
  Else If (mstj(44)==1) Then
    pmsqcd = 4.*alams*(0.25*pms/alams)**(rlu(0)**(b0/fbr))
  Else
    pmsqcd = pms*rlu(0)**(alfm*b0/fbr)
  End If
  If (zc>0.49 .Or. pmsqcd<=pmth(4,kfl(1))**2) pmsqcd = pmth(2, kfl(1))**2
  v(iep(1), 5) = pmsqcd
  mce = 1

!...Select mass for daughter in QED evolution.
  If (mstj(41)==2 .And. kfl(1)>=1 .And. kfl(1)<=8) Then
    pmsqed = pms*exp(max(-100.,log(rlu(0))*paru(2)/(paru(101)*fbre)))
    If (zce>0.49 .Or. pmsqed<=pmth(5,kfl(1))**2) pmsqed = pmth(2, kfl(1))**2
    If (pmsqed>pmsqcd) Then
      v(iep(1), 5) = pmsqed
      mce = 2
    End If
  End If

!...Check whether daughter mass below cutoff.
  p(iep(1), 5) = sqrt(v(iep(1),5))
  If (p(iep(1),5)<=pmth(3,kfl(1))) Then
    p(iep(1), 5) = pmth(1, kfl(1))
    v(iep(1), 5) = p(iep(1), 5)**2
    Goto 300
  End If

!...Select z value of branching: q -> qgamma.
  If (mce==2) Then
    z = 1. - (1.-zce)*(zce/(1.-zce))**rlu(0)
    If (1.+z**2<2.*rlu(0)) Goto 260
    k(iep(1), 5) = 22

!...Select z value of branching: q -> qg, g -> gg, g -> qqbar.
  Else If (mstj(49)/=1 .And. kfl(1)/=21) Then
    z = 1. - (1.-zc)*(zc/(1.-zc))**rlu(0)
    If (1.+z**2<2.*rlu(0)) Goto 260
    k(iep(1), 5) = 21
  Else If (mstj(49)==0 .And. mstj(45)*(0.5-zc)<rlu(0)*fbr) Then
    z = (1.-zc)*(zc/(1.-zc))**rlu(0)
    If (rlu(0)>0.5) z = 1. - z
    If ((1.-z*(1.-z))**2<rlu(0)) Goto 260
    k(iep(1), 5) = 21
  Else If (mstj(49)/=1) Then
    z = zc + (1.-2.*zc)*rlu(0)
    If (z**2+(1.-z)**2<rlu(0)) Goto 260
    kflb = 1 + int(mstj(45)*rlu(0))
    pmq = 4.*pmth(2, kflb)**2/v(iep(1), 5)
    If (pmq>=1.) Goto 260
    pmq0 = 4.*pmth(2, 21)**2/v(iep(1), 5)
    If (mod(mstj(43),2)==0 .And. (1.+0.5*pmq)*sqrt(1.-pmq)<rlu(0)*(1.+0.5*pmq0)*sqrt(1.-pmq0)) Goto 260
    k(iep(1), 5) = kflb

!...Ditto for scalar gluon model.
  Else If (kfl(1)/=21) Then
    z = 1. - sqrt(zc**2+rlu(0)*(1.-2.*zc))
    k(iep(1), 5) = 21
  Else If (rlu(0)*(parj(87)+mstj(45)*parj(88))<=parj(87)) Then
    z = zc + (1.-2.*zc)*rlu(0)
    k(iep(1), 5) = 21
  Else
    z = zc + (1.-2.*zc)*rlu(0)
    kflb = 1 + int(mstj(45)*rlu(0))
    pmq = 4.*pmth(2, kflb)**2/v(iep(1), 5)
    If (pmq>=1.) Goto 260
    k(iep(1), 5) = kflb
  End If
  If (mce==1 .And. mstj(44)>=2) Then
    If (z*(1.-z)*v(iep(1),5)<pt2min) Goto 260
    If (alfm/log(v(iep(1),5)*z*(1.-z)/alams)<rlu(0)) Goto 260
  End If

!...Check if z consistent with chosen m.
  If (kfl(1)==21) Then
    kflgd1 = iabs(k(iep(1),5))
    kflgd2 = kflgd1
  Else
    kflgd1 = kfl(1)
    kflgd2 = iabs(k(iep(1),5))
  End If
  If (nep==1) Then
    ped = ps(4)
  Else If (nep>=3) Then
    ped = p(iep(1), 4)
  Else If (igm==0 .Or. mstj(43)<=2) Then
    ped = 0.5*(v(im,5)+v(iep(1),5)-pm2**2)/p(im, 5)
  Else
    If (iep(1)==n+1) ped = v(im, 1)*pem
    If (iep(1)==n+2) ped = (1.-v(im,1))*pem
  End If
  If (mod(mstj(43),2)==1) Then
    pmqth3 = 0.5*parj(82)
    If (kflgd2==22) pmqth3 = 0.5*parj(83)
    pmq1 = (pmth(1,kflgd1)**2+pmqth3**2)/v(iep(1), 5)
    pmq2 = (pmth(1,kflgd2)**2+pmqth3**2)/v(iep(1), 5)
    zd = sqrt(max(0.,(1.-v(iep(1),5)/ped**2)*((1.-pmq1-pmq2)**2-4.*pmq1*pmq2)))
    zh = 1. + pmq1 - pmq2
  Else
    zd = sqrt(max(0.,1.-v(iep(1),5)/ped**2))
    zh = 1.
  End If
  zl = 0.5*(zh-zd)
  zu = 0.5*(zh+zd)
  If (z<zl .Or. z>zu) Goto 260
  If (kfl(1)==21) v(iep(1), 3) = log(zu*(1.-zl)/max(1E-20,zl*(1.-zu)))
  If (kfl(1)/=21) v(iep(1), 3) = log((1.-zl)/max(1E-10,1.-zu))

!...Three-jet matrix element correction.
  If (igm==0 .And. m3jc==1) Then
    x1 = z*(1.+v(iep(1),5)/v(ns+1,5))
    x2 = 1. - v(iep(1), 5)/v(ns+1, 5)
    x3 = (1.-x1) + (1.-x2)
    If (mce==2) Then
      ki1 = k(ipa(inum), 2)
      ki2 = k(ipa(3-inum), 2)
      qf1 = kchg(iabs(ki1), 1)*isign(1, ki1)/3.
      qf2 = kchg(iabs(ki2), 1)*isign(1, ki2)/3.
      wshow = qf1**2*(1.-x1)/x3*(1.+(x1/(2.-x2))**2) + qf2**2*(1.-x2)/x3*(1.+(x2/(2.-x1))**2)
      wme = (qf1*(1.-x1)/x3-qf2*(1.-x2)/x3)**2*(x1**2+x2**2)
    Else If (mstj(49)/=1) Then
      wshow = 1. + (1.-x1)/x3*(x1/(2.-x2))**2 + (1.-x2)/x3*(x2/(2.-x1))**2
      wme = x1**2 + x2**2
    Else
      wshow = 4.*x3*((1.-x1)/(2.-x2)**2+(1.-x2)/(2.-x1)**2)
      wme = x3**2
    End If
    If (wme<rlu(0)*wshow) Goto 260

!...Impose angular ordering by rejection of nonordered emission.
  Else If (mce==1 .And. igm>0 .And. mstj(42)>=2) Then
    maom = 1
    zm = v(im, 1)
    If (iep(1)==n+2) zm = 1. - v(im, 1)
    the2id = z*(1.-z)*(zm*p(im,4))**2/v(iep(1), 5)
    iaom = im
    290 If (k(iaom,5)==22) Then
      iaom = k(iaom, 3)
      If (k(iaom,3)<=ns) maom = 0
      If (maom==1) Goto 290
    End If
    If (maom==1) Then
      the2im = v(iaom, 1)*(1.-v(iaom,1))*p(iaom, 4)**2/v(iaom, 5)
      If (the2id<the2im) Goto 260
    End If
  End If

!...Impose user-defined maximum angle at first branching.
  If (mstj(48)==1) Then
    If (nep==1 .And. im==ns) Then
      the2id = z*(1.-z)*ps(4)**2/v(iep(1), 5)
      If (the2id<1./parj(85)**2) Goto 260
    Else If (nep==2 .And. iep(1)==ns+2) Then
      the2id = z*(1.-z)*(0.5*p(im,4))**2/v(iep(1), 5)
      If (the2id<1./parj(85)**2) Goto 260
    Else If (nep==2 .And. iep(1)==ns+3) Then
      the2id = z*(1.-z)*(0.5*p(im,4))**2/v(iep(1), 5)
      If (the2id<1./parj(86)**2) Goto 260
    End If
  End If

!...End of inner veto algorithm. Check if only one leg evolved so far.
  300 v(iep(1), 1) = z
  isl(1) = 0
  isl(2) = 0
  If (nep==1) Goto 330
  If (nep==2 .And. p(iep(1),5)+p(iep(2),5)>=p(im,5)) Goto 200
  Do i = 1, nep
    If (itry(i)==0 .And. kfld(i)>0 .And. (kfld(i)<=8 .Or. kfld(i)==21)) Then
      If (p(n+i,5)>=pmth(2,kfld(i))) Goto 200
    End If
  End Do

!...Check if chosen multiplet m1,m2,z1,z2 is physical.
  If (nep==3) Then
    pa1s = (p(n+1,4)+p(n+1,5))*(p(n+1,4)-p(n+1,5))
    pa2s = (p(n+2,4)+p(n+2,5))*(p(n+2,4)-p(n+2,5))
    pa3s = (p(n+3,4)+p(n+3,5))*(p(n+3,4)-p(n+3,5))
    pts = 0.25*(2.*pa1s*pa2s+2.*pa1s*pa3s+2.*pa2s*pa3s-pa1s**2-pa2s**2-pa3s**2)/pa1s
    If (pts<=0.) Goto 200
  Else If (igm==0 .Or. mstj(43)<=2 .Or. mod(mstj(43),2)==0) Then
    Do i1 = n + 1, n + 2
      kflda = iabs(k(i1,2))
      If (kflda==0 .Or. (kflda>8 .And. kflda/=21)) Goto 320
      If (p(i1,5)<pmth(2,kflda)) Goto 320
      If (kflda==21) Then
        kflgd1 = iabs(k(i1,5))
        kflgd2 = kflgd1
      Else
        kflgd1 = kflda
        kflgd2 = iabs(k(i1,5))
      End If
      i2 = 2*n + 3 - i1
      If (igm==0 .Or. mstj(43)<=2) Then
        ped = 0.5*(v(im,5)+v(i1,5)-v(i2,5))/p(im, 5)
      Else
        If (i1==n+1) zm = v(im, 1)
        If (i1==n+2) zm = 1. - v(im, 1)
        pml = sqrt((v(im,5)-v(n+1,5)-v(n+2,5))**2-4.*v(n+1,5)*v(n+2,5))
        ped = pem*(0.5*(v(im,5)-pml+v(i1,5)-v(i2,5))+pml*zm)/v(im, 5)
      End If
      If (mod(mstj(43),2)==1) Then
        pmqth3 = 0.5*parj(82)
        If (kflgd2==22) pmqth3 = 0.5*parj(83)
        pmq1 = (pmth(1,kflgd1)**2+pmqth3**2)/v(i1, 5)
        pmq2 = (pmth(1,kflgd2)**2+pmqth3**2)/v(i1, 5)
        zd = sqrt(max(0.,(1.-v(i1,5)/ped**2)*((1.-pmq1-pmq2)**2-4.*pmq1*pmq2)))
        zh = 1. + pmq1 - pmq2
      Else
        zd = sqrt(max(0.,1.-v(i1,5)/ped**2))
        zh = 1.
      End If
      zl = 0.5*(zh-zd)
      zu = 0.5*(zh+zd)
      If (i1==n+1 .And. (v(i1,1)<zl .Or. v(i1,1)>zu)) isl(1) = 1
      If (i1==n+2 .And. (v(i1,1)<zl .Or. v(i1,1)>zu)) isl(2) = 1
      If (kflda==21) v(i1, 4) = log(zu*(1.-zl)/max(1E-20,zl*(1.-zu)))
      If (kflda/=21) v(i1, 4) = log((1.-zl)/max(1E-10,1.-zu))
    320 End Do
    If (isl(1)==1 .And. isl(2)==1 .And. islm/=0) Then
      isl(3-islm) = 0
      islm = 3 - islm
    Else If (isl(1)==1 .And. isl(2)==1) Then
      zdr1 = max(0., v(n+1,3)/v(n+1,4)-1.)
      zdr2 = max(0., v(n+2,3)/v(n+2,4)-1.)
      If (zdr2>rlu(0)*(zdr1+zdr2)) isl(1) = 0
      If (isl(1)==1) isl(2) = 0
      If (isl(1)==0) islm = 1
      If (isl(2)==0) islm = 2
    End If
    If (isl(1)==1 .Or. isl(2)==1) Goto 200
  End If
  If (igm>0 .And. mod(mstj(43),2)==1 .And. (p(n+1,5)>=pmth(2,kfld(1)) .Or. p(n+2,5)>=pmth(2,kfld(2)))) Then
    pmq1 = v(n+1, 5)/v(im, 5)
    pmq2 = v(n+2, 5)/v(im, 5)
    zd = sqrt(max(0.,(1.-v(im,5)/pem**2)*((1.-pmq1-pmq2)**2-4.*pmq1*pmq2)))
    zh = 1. + pmq1 - pmq2
    zl = 0.5*(zh-zd)
    zu = 0.5*(zh+zd)
    If (v(im,1)<zl .Or. v(im,1)>zu) Goto 200
  End If

!...Accepted branch. Construct four-momentum for initial partons.
  330 mazip = 0
  mazic = 0
  If (nep==1) Then
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = sqrt(max(0.,(p(ipa(1),4)+p(n+1,5))*(p(ipa(1),4)-p(n+1,5))))
    p(n+1, 4) = p(ipa(1), 4)
    v(n+1, 2) = p(n+1, 4)
  Else If (igm==0 .And. nep==2) Then
    ped1 = 0.5*(v(im,5)+v(n+1,5)-v(n+2,5))/p(im, 5)
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = sqrt(max(0.,(ped1+p(n+1,5))*(ped1-p(n+1,5))))
    p(n+1, 4) = ped1
    p(n+2, 1) = 0.
    p(n+2, 2) = 0.
    p(n+2, 3) = -p(n+1, 3)
    p(n+2, 4) = p(im, 5) - ped1
    v(n+1, 2) = p(n+1, 4)
    v(n+2, 2) = p(n+2, 4)
  Else If (nep==3) Then
    p(n+1, 1) = 0.
    p(n+1, 2) = 0.
    p(n+1, 3) = sqrt(max(0.,pa1s))
    p(n+2, 1) = sqrt(pts)
    p(n+2, 2) = 0.
    p(n+2, 3) = 0.5*(pa3s-pa2s-pa1s)/p(n+1, 3)
    p(n+3, 1) = -p(n+2, 1)
    p(n+3, 2) = 0.
    p(n+3, 3) = -(p(n+1,3)+p(n+2,3))
    v(n+1, 2) = p(n+1, 4)
    v(n+2, 2) = p(n+2, 4)
    v(n+3, 2) = p(n+3, 4)

!...Construct transverse momentum for ordinary branching in shower.
  Else
    zm = v(im, 1)
    pzm = sqrt(max(0.,(pem+p(im,5))*(pem-p(im,5))))
    pmls = (v(im,5)-v(n+1,5)-v(n+2,5))**2 - 4.*v(n+1, 5)*v(n+2, 5)
    If (pzm<=0.) Then
      pts = 0.
    Else If (mod(mstj(43),2)==1) Then
      pts = (pem**2*(zm*(1.-zm)*v(im,5)-(1.-zm)*v(n+1,5)-zm*v(n+2,5))-0.25*pmls)/pzm**2
    Else
      pts = pmls*(zm*(1.-zm)*pem**2/v(im,5)-0.25)/pzm**2
    End If
    pt = sqrt(max(0.,pts))

!...Find coefficient of azimuthal asymmetry due to gluon polarization.
    hazip = 0.
    If (mstj(49)/=1 .And. mod(mstj(46),2)==1 .And. k(im,2)==21 .And. iau/=0) Then
      If (k(igm,3)/=0) mazip = 1
      zau = v(igm, 1)
      If (iau==im+1) zau = 1. - v(igm, 1)
      If (mazip==0) zau = 0.
      If (k(igm,2)/=21) Then
        hazip = 2.*zau/(1.+zau**2)
      Else
        hazip = (zau/(1.-zau*(1.-zau)))**2
      End If
      If (k(n+1,2)/=21) Then
        hazip = hazip*(-2.*zm*(1.-zm))/(1.-2.*zm*(1.-zm))
      Else
        hazip = hazip*(zm*(1.-zm)/(1.-zm*(1.-zm)))**2
      End If
    End If

!...Find coefficient of azimuthal asymmetry due to soft gluon
!...interference.
    hazic = 0.
    If (mstj(46)>=2 .And. (k(n+1,2)==21 .Or. k(n+2,2)==21) .And. iau/=0) Then
      If (k(igm,3)/=0) mazic = n + 1
      If (k(igm,3)/=0 .And. k(n+1,2)/=21) mazic = n + 2
      If (k(igm,3)/=0 .And. k(n+1,2)==21 .And. k(n+2,2)==21 .And. zm>0.5) mazic = n + 2
      If (k(iau,2)==22) mazic = 0
      zs = zm
      If (mazic==n+2) zs = 1. - zm
      zgm = v(igm, 1)
      If (iau==im-1) zgm = 1. - v(igm, 1)
      If (mazic==0) zgm = 1.
      hazic = (p(im,5)/p(igm,5))*sqrt((1.-zs)*(1.-zgm)/(zs*zgm))
      hazic = min(0.95, hazic)
    End If
  End If

!...Construct kinematics for ordinary branching in shower.
  340 If (nep==2 .And. igm>0) Then
    If (mod(mstj(43),2)==1) Then
      p(n+1, 4) = pem*v(im, 1)
    Else
      p(n+1, 4) = pem*(0.5*(v(im,5)-sqrt(pmls)+v(n+1,5)-v(n+2,5))+sqrt(pmls)*zm)/v(im, 5)
    End If
    phi = paru(2)*rlu(0)
    p(n+1, 1) = pt*cos(phi)
    p(n+1, 2) = pt*sin(phi)
    If (pzm>0.) Then
      p(n+1, 3) = 0.5*(v(n+2,5)-v(n+1,5)-v(im,5)+2.*pem*p(n+1,4))/pzm
    Else
      p(n+1, 3) = 0.
    End If
    p(n+2, 1) = -p(n+1, 1)
    p(n+2, 2) = -p(n+1, 2)
    p(n+2, 3) = pzm - p(n+1, 3)
    p(n+2, 4) = pem - p(n+1, 4)
    If (mstj(43)<=2) Then
      v(n+1, 2) = (pem*p(n+1,4)-pzm*p(n+1,3))/p(im, 5)
      v(n+2, 2) = (pem*p(n+2,4)-pzm*p(n+2,3))/p(im, 5)
    End If
  End If

!...Rotate and boost daughters.
  If (igm>0) Then
    If (mstj(43)<=2) Then
      bex = p(igm, 1)/p(igm, 4)
      bey = p(igm, 2)/p(igm, 4)
      bez = p(igm, 3)/p(igm, 4)
      ga = p(igm, 4)/p(igm, 5)
      gabep = ga*(ga*(bex*p(im,1)+bey*p(im,2)+bez*p(im,3))/(1.+ga)-p(im,4))
    Else
      bex = 0.
      bey = 0.
      bez = 0.
      ga = 1.
      gabep = 0.
    End If
    the = ulangl(p(im,3)+gabep*bez, sqrt((p(im,1)+gabep*bex)**2+(p(im,2)+gabep*bey)**2))
    phi = ulangl(p(im,1)+gabep*bex, p(im,2)+gabep*bey)
    Do i = n + 1, n + 2
      dp(1) = dble(cos(the)*cos(phi)*p(i,1)-sin(phi)*p(i,2)+sin(the)*cos(phi)*p(i,3))
      dp(2) = dble(cos(the)*sin(phi)*p(i,1)+cos(phi)*p(i,2)+sin(the)*sin(phi)*p(i,3))
      dp(3) = dble(-sin(the)*p(i,1)+cos(the)*p(i,3))
      dp(4) = dble(p(i,4))
      dbp = dble(bex)*dp(1) + dble(bey)*dp(2) + dble(bez)*dp(3)
      dgabp = dble(ga)*(dble(ga)*dbp/(1D0+dble(ga))+dp(4))
      p(i, 1) = sngl(dp(1)+dgabp*dble(bex))
      p(i, 2) = sngl(dp(2)+dgabp*dble(bey))
      p(i, 3) = sngl(dp(3)+dgabp*dble(bez))
      p(i, 4) = ga*sngl(dp(4)+dbp)
    End Do
  End If

!...Weight with azimuthal distribution, if required.
  If (mazip/=0 .Or. mazic/=0) Then
    Do j = 1, 3
      dpt(1, j) = dble(p(im,j))
      dpt(2, j) = dble(p(iau,j))
      dpt(3, j) = dble(p(n+1,j))
    End Do
    dpma = dpt(1, 1)*dpt(2, 1) + dpt(1, 2)*dpt(2, 2) + dpt(1, 3)*dpt(2, 3)
    dpmd = dpt(1, 1)*dpt(3, 1) + dpt(1, 2)*dpt(3, 2) + dpt(1, 3)*dpt(3, 3)
    dpmm = dpt(1, 1)**2 + dpt(1, 2)**2 + dpt(1, 3)**2
    Do j = 1, 3
      dpt(4, j) = dpt(2, j) - dpma*dpt(1, j)/dpmm
      dpt(5, j) = dpt(3, j) - dpmd*dpt(1, j)/dpmm
    End Do
    dpt(4, 4) = dsqrt(dpt(4,1)**2+dpt(4,2)**2+dpt(4,3)**2)
    dpt(5, 4) = dsqrt(dpt(5,1)**2+dpt(5,2)**2+dpt(5,3)**2)
!lin-5/2012:
!        IF(MIN(DPT(4,4),DPT(5,4)).GT.0.1*PARJ(82)) THEN
    If (sngl(min(dpt(4,4),dpt(5,4)))>(0.1*parj(82))) Then
      cad = sngl((dpt(4,1)*dpt(5,1)+dpt(4,2)*dpt(5,2)+dpt(4,3)*dpt(5,3))/(dpt(4,4)*dpt(5,4)))
      If (mazip/=0) Then
        If (1.+hazip*(2.*cad**2-1.)<rlu(0)*(1.+abs(hazip))) Goto 340
      End If
      If (mazic/=0) Then
        If (mazic==n+2) cad = -cad
        If ((1.-hazic)*(1.-hazic*cad)/(1.+hazic**2-2.*hazic*cad)<rlu(0)) Goto 340
      End If
    End If
  End If

!...Continue loop over partons that may branch, until none left.
  If (igm>=0) k(im, 1) = 14
  n = n + nep
  nep = 2
  If (n>mstu(4)-mstu(32)-5) Then
    Call luerrm(11, '(LUSHOW:) no more memory left in LUJETS')
    If (mstu(21)>=1) n = ns
    If (mstu(21)>=1) Return
  End If
  Goto 140

!...Set information on imagined shower initiator.
  380 If (npa>=2) Then
    k(ns+1, 1) = 11
    k(ns+1, 2) = 94
    k(ns+1, 3) = ip1
    If (ip2>0 .And. ip2<ip1) k(ns+1, 3) = ip2
    k(ns+1, 4) = ns + 2
    k(ns+1, 5) = ns + 1 + npa
    iim = 1
  Else
    iim = 0
  End If

!...Reconstruct string drawing information.
  Do i = ns + 1 + iim, n
    If (k(i,1)<=10 .And. k(i,2)==22) Then
      k(i, 1) = 1
    Else If (k(i,1)<=10) Then
      k(i, 4) = mstu(5)*(k(i,4)/mstu(5))
      k(i, 5) = mstu(5)*(k(i,5)/mstu(5))
    Else If (k(mod(k(i,4),mstu(5))+1,2)/=22) Then
      id1 = mod(k(i,4), mstu(5))
      If (k(i,2)>=1 .And. k(i,2)<=8) id1 = mod(k(i,4), mstu(5)) + 1
      id2 = 2*mod(k(i,4), mstu(5)) + 1 - id1
      k(i, 4) = mstu(5)*(k(i,4)/mstu(5)) + id1
      k(i, 5) = mstu(5)*(k(i,5)/mstu(5)) + id2
      k(id1, 4) = k(id1, 4) + mstu(5)*i
      k(id1, 5) = k(id1, 5) + mstu(5)*id2
      k(id2, 4) = k(id2, 4) + mstu(5)*id1
      k(id2, 5) = k(id2, 5) + mstu(5)*i
    Else
      id1 = mod(k(i,4), mstu(5))
      id2 = id1 + 1
      k(i, 4) = mstu(5)*(k(i,4)/mstu(5)) + id1
      k(i, 5) = mstu(5)*(k(i,5)/mstu(5)) + id1
      k(id1, 4) = k(id1, 4) + mstu(5)*i
      k(id1, 5) = k(id1, 5) + mstu(5)*i
      k(id2, 4) = 0
      k(id2, 5) = 0
    End If
  End Do

!...Transformation from CM frame.
  If (npa>=2) Then
    bex = ps(1)/ps(4)
    bey = ps(2)/ps(4)
    bez = ps(3)/ps(4)
    ga = ps(4)/ps(5)
    gabep = ga*(ga*(bex*p(ipa(1),1)+bey*p(ipa(1),2)+bez*p(ipa(1),3))/(1.+ga)-p(ipa(1),4))
  Else
    bex = 0.
    bey = 0.
    bez = 0.
    gabep = 0.
  End If
  the = ulangl(p(ipa(1),3)+gabep*bez, sqrt((p(ipa(1),1)+gabep*bex)**2+(p(ipa(1),2)+gabep*bey)**2))
  phi = ulangl(p(ipa(1),1)+gabep*bex, p(ipa(1),2)+gabep*bey)
  If (npa==3) Then
    chi = ulangl(cos(the)*cos(phi)*(p(ipa(2),1)+gabep*bex)+cos(the)*sin(phi)*(p(ipa(2),2)+gabep*bey)-sin(the)*(p(ipa(2),3)+gabep*bez), -sin(phi)*(p(ipa(2),1)+gabep*bex)+cos(phi)*(p(ipa(2),2)+gabep*bey))
    Call ludbrb(ns+1, n, 0., chi, 0D0, 0D0, 0D0)
  End If
  dbex = dble(bex)
  dbey = dble(bey)
  dbez = dble(bez)
  Call ludbrb(ns+1, n, the, phi, dbex, dbey, dbez)

!...Decay vertex of shower.
  Do i = ns + 1, n
    Do j = 1, 5
      v(i, j) = v(ip1, j)
    End Do
  End Do

!...Delete trivial shower, else connect initiators.
  If (n==ns+npa+iim) Then
    n = ns
  Else
    Do ip = 1, npa
      k(ipa(ip), 1) = 14
      k(ipa(ip), 4) = k(ipa(ip), 4) + ns + iim + ip
      k(ipa(ip), 5) = k(ipa(ip), 5) + ns + iim + ip
      k(ns+iim+ip, 3) = ipa(ip)
      If (iim==1 .And. mstu(16)/=2) k(ns+iim+ip, 3) = ns + 1
      k(ns+iim+ip, 4) = mstu(5)*ipa(ip) + k(ns+iim+ip, 4)
      k(ns+iim+ip, 5) = mstu(5)*ipa(ip) + k(ns+iim+ip, 5)
    End Do
  End If

  Return
End Subroutine lushow

!*********************************************************************

Subroutine luboei(nsav)

!...Purpose: to modify event so as to approximately take into account
!...Bose-Einstein effects according to a simple phenomenological
!...parametrization.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Dimension dps(4), kfbe(9), nbe(0:9), bei(100)
  Data kfbe/211, -211, 111, 321, -321, 130, 310, 221, 331/

!...Boost event to overall CM frame. Calculate CM energy.
  If ((mstj(51)/=1 .And. mstj(51)/=2) .Or. n-nsav<=1) Return
  Do j = 1, 4
    dps(j) = 0.D0
  End Do
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 120
    Do j = 1, 4
      dps(j) = dps(j) + dble(p(i,j))
    End Do
  120 End Do
  Call ludbrb(0, 0, 0., 0., -dps(1)/dps(4), -dps(2)/dps(4), -dps(3)/dps(4))
  pecm = 0.
  Do i = 1, n
    If (k(i,1)>=1 .And. k(i,1)<=10) pecm = pecm + p(i, 4)
  End Do

!...Reserve copy of particles by species at end of record.
  nbe(0) = n + mstu(3)
  Do ibe = 1, min(9, mstj(51))
    nbe(ibe) = nbe(ibe-1)
    Do i = nsav + 1, n
      If (k(i,2)/=kfbe(ibe)) Goto 150
      If (k(i,1)<=0 .Or. k(i,1)>10) Goto 150
      If (nbe(ibe)>=mstu(4)-mstu(32)-5) Then
        Call luerrm(11, '(LUBOEI:) no more memory left in LUJETS')
        Return
      End If
      nbe(ibe) = nbe(ibe) + 1
      k(nbe(ibe), 1) = i
      Do j = 1, 3
        p(nbe(ibe), j) = 0.
      End Do
    150 End Do
  End Do

!...Tabulate integral for subsequent momentum shift.
  Do ibe = 1, min(9, mstj(51))
    If (ibe/=1 .And. ibe/=4 .And. ibe<=7) Goto 180
    If (ibe==1 .And. max(nbe(1)-nbe(0),nbe(2)-nbe(1),nbe(3)-nbe(2))<=1) Goto 180
    If (ibe==4 .And. max(nbe(4)-nbe(3),nbe(5)-nbe(4),nbe(6)-nbe(5),nbe(7)-nbe(6))<=1) Goto 180
    If (ibe>=8 .And. nbe(ibe)-nbe(ibe-1)<=1) Goto 180
    If (ibe==1) pmhq = 2.*ulmass(211)
    If (ibe==4) pmhq = 2.*ulmass(321)
    If (ibe==8) pmhq = 2.*ulmass(221)
    If (ibe==9) pmhq = 2.*ulmass(331)
    qdel = 0.1*min(pmhq, parj(93))
    If (mstj(51)==1) Then
      nbin = min(100, nint(9.*parj(93)/qdel))
      beex = exp(0.5*qdel/parj(93))
      bert = exp(-qdel/parj(93))
    Else
      nbin = min(100, nint(3.*parj(93)/qdel))
    End If
    Do ibin = 1, nbin
      qbin = qdel*(ibin-0.5)
      bei(ibin) = qdel*(qbin**2+qdel**2/12.)/sqrt(qbin**2+pmhq**2)
      If (mstj(51)==1) Then
        beex = beex*bert
        bei(ibin) = bei(ibin)*beex
      Else
        bei(ibin) = bei(ibin)*exp(-(qbin/parj(93))**2)
      End If
      If (ibin>=2) bei(ibin) = bei(ibin) + bei(ibin-1)
    End Do

!...Loop through particle pairs and find old relative momentum.
    180 Do i1m = nbe(ibe-1) + 1, nbe(ibe) - 1
      i1 = k(i1m, 1)
      Do i2m = i1m + 1, nbe(ibe)
        i2 = k(i2m, 1)
        q2old = max(0., (p(i1,4)+p(i2,4))**2-(p(i1,1)+p(i2,1))**2-(p(i1,2)+p(i2,2))**2-(p(i1,3)+p(i2,3))**2-(p(i1,5)+p(i2,5))**2)
        qold = sqrt(q2old)

!...Calculate new relative momentum.
        If (qold<0.5*qdel) Then
          qmov = qold/3.
        Else If (qold<(nbin-0.1)*qdel) Then
          rbin = qold/qdel
          ibin = int(rbin)
          rinp = (rbin**3-ibin**3)/(3*ibin*(ibin+1)+1)
          qmov = (bei(ibin)+rinp*(bei(ibin+1)-bei(ibin)))*sqrt(q2old+pmhq**2)/q2old
        Else
          qmov = bei(nbin)*sqrt(q2old+pmhq**2)/q2old
        End If
        q2new = q2old*(qold/(qold+3.*parj(92)*qmov))**(2./3.)

!...Calculate and save shift to be performed on three-momenta.
        hc1 = (p(i1,4)+p(i2,4))**2 - (q2old-q2new)
        hc2 = (q2old-q2new)*(p(i1,4)-p(i2,4))**2
        ha = 0.5*(1.-sqrt(hc1*q2new/(hc1*q2old-hc2)))
        Do j = 1, 3
          pd = ha*(p(i2,j)-p(i1,j))
          p(i1m, j) = p(i1m, j) + pd
          p(i2m, j) = p(i2m, j) - pd
        End Do
      End Do
    End Do
  End Do

!...Shift momenta and recalculate energies.
  Do im = nbe(0) + 1, nbe(min(9,mstj(51)))
    i = k(im, 1)
    Do j = 1, 3
      p(i, j) = p(i, j) + p(im, j)
    End Do
    p(i, 4) = sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  End Do

!...Rescale all momenta for energy conservation.
  pes = 0.
  pqs = 0.
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 240
    pes = pes + p(i, 4)
    pqs = pqs + p(i, 5)**2/p(i, 4)
  240 End Do
  fac = (pecm-pqs)/(pes-pqs)
  Do i = 1, n
    If (k(i,1)<=0 .Or. k(i,1)>10) Goto 260
    Do j = 1, 3
      p(i, j) = fac*p(i, j)
    End Do
    p(i, 4) = sqrt(p(i,5)**2+p(i,1)**2+p(i,2)**2+p(i,3)**2)
  260 End Do

!...Boost back to correct reference frame.
  Call ludbrb(0, 0, 0., 0., dps(1)/dps(4), dps(2)/dps(4), dps(3)/dps(4))

  Return
End Subroutine luboei

!*********************************************************************



!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Function ulmass(kf)

!...Purpose: to give the mass of a particle/parton.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/

!...Reset variables. Compressed code.
  ulmass = 0.
  kfa = iabs(kf)
  kc = lucomp(kf)
  If (kc==0) Return
  parf(106) = pmas(6, 1)
  parf(107) = pmas(7, 1)
  parf(108) = pmas(8, 1)

!...Guarantee use of constituent masses for internal checks.
  If ((mstj(93)==1 .Or. mstj(93)==2) .And. kfa<=10) Then
    ulmass = parf(100+kfa)
    If (mstj(93)==2) ulmass = max(0., ulmass-parf(121))

!...Masses that can be read directly off table.
  Else If (kfa<=100 .Or. kc<=80 .Or. kc>100) Then
    ulmass = pmas(kc, 1)

!...Find constituent partons and their masses.
  Else
    kfla = mod(kfa/1000, 10)
    kflb = mod(kfa/100, 10)
    kflc = mod(kfa/10, 10)
    kfls = mod(kfa, 10)
    kflr = mod(kfa/10000, 10)
    pma = parf(100+kfla)
    pmb = parf(100+kflb)
    pmc = parf(100+kflc)

!...Construct masses for various meson, diquark and baryon cases.
    If (kfla==0 .And. kflr==0 .And. kfls<=3) Then
      If (kfls==1) pmspl = -3./(pmb*pmc)
      If (kfls>=3) pmspl = 1./(pmb*pmc)
      ulmass = parf(111) + pmb + pmc + parf(113)*parf(101)**2*pmspl
    Else If (kfla==0) Then
      kmul = 2
      If (kfls==1) kmul = 3
      If (kflr==2) kmul = 4
      If (kfls==5) kmul = 5
      ulmass = parf(113+kmul) + pmb + pmc
    Else If (kflc==0) Then
      If (kfls==1) pmspl = -3./(pma*pmb)
      If (kfls==3) pmspl = 1./(pma*pmb)
      ulmass = 2.*parf(112)/3. + pma + pmb + parf(114)*parf(101)**2*pmspl
      If (mstj(93)==1) ulmass = pma + pmb
      If (mstj(93)==2) ulmass = max(0., ulmass-parf(122)-2.*parf(112)/3.)
    Else
      If (kfls==2 .And. kfla==kflb) Then
        pmspl = 1./(pma*pmb) - 2./(pma*pmc) - 2./(pmb*pmc)
      Else If (kfls==2 .And. kflb>=kflc) Then
        pmspl = -2./(pma*pmb) - 2./(pma*pmc) + 1./(pmb*pmc)
      Else If (kfls==2) Then
        pmspl = -3./(pmb*pmc)
      Else
        pmspl = 1./(pma*pmb) + 1./(pma*pmc) + 1./(pmb*pmc)
      End If
      ulmass = parf(112) + pma + pmb + pmc + parf(114)*parf(101)**2*pmspl
    End If
  End If

!...Optional mass broadening according to truncated Breit-Wigner
!...(either in m or in m^2).
  If (mstj(24)>=1 .And. pmas(kc,2)>1E-4) Then
    If (mstj(24)==1 .Or. (mstj(24)==2 .And. kfa>100)) Then
      ulmass = ulmass + 0.5*pmas(kc, 2)*tan((2.*rlu(0)-1.)*atan(2.*pmas(kc,3)/pmas(kc,2)))
    Else
      pm0 = ulmass
      pmlow = atan((max(0.,pm0-pmas(kc,3))**2-pm0**2)/(pm0*pmas(kc,2)))
      pmupp = atan((pm0+pmas(kc,3))**2-pm0**2)/(pm0*pmas(kc,2))
      ulmass = sqrt(max(0.,pm0**2+pm0*pmas(kc,2)*tan(pmlow+(pmupp-pmlow)*rlu(0))))
    End If
  End If
  mstj(93) = 0

  Return
End Function ulmass

!*********************************************************************

Subroutine luname(kf, chau)

!...Purpose: to give the particle/parton name as a character string.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat4/chaf(500)
  Character chaf*8
  Save /ludat4/
  Character chau*16

!...Initial values. Charge. Subdivide code.
  chau = ' '
  kfa = iabs(kf)
  kc = lucomp(kf)
  If (kc==0) Return
  kq = luchge(kf)
  kfla = mod(kfa/1000, 10)
  kflb = mod(kfa/100, 10)
  kflc = mod(kfa/10, 10)
  kfls = mod(kfa, 10)
  kflr = mod(kfa/10000, 10)

!...Read out root name and spin for simple particle.
  If (kfa<=100 .Or. (kfa>100 .And. kc>100)) Then
    chau = chaf(kc)
    len = 0
    Do lem = 1, 8
      If (chau(lem:lem)/=' ') len = lem
    End Do

!...Construct root name for diquark. Add on spin.
  Else If (kflc==0) Then
    chau(1:2) = chaf(kfla)(1:1) // chaf(kflb)(1:1)
    If (kfls==1) chau(3:4) = '_0'
    If (kfls==3) chau(3:4) = '_1'
    len = 4

!...Construct root name for heavy meson. Add on spin and heavy flavour.
  Else If (kfla==0) Then
    If (kflb==5) chau(1:1) = 'B'
    If (kflb==6) chau(1:1) = 'T'
    If (kflb==7) chau(1:1) = 'L'
    If (kflb==8) chau(1:1) = 'H'
    len = 1
    If (kflr==0 .And. kfls==1) Then
    Else If (kflr==0 .And. kfls==3) Then
      chau(2:2) = '*'
      len = 2
    Else If (kflr==1 .And. kfls==3) Then
      chau(2:3) = '_1'
      len = 3
    Else If (kflr==1 .And. kfls==1) Then
      chau(2:4) = '*_0'
      len = 4
    Else If (kflr==2) Then
      chau(2:4) = '*_1'
      len = 4
    Else If (kfls==5) Then
      chau(2:4) = '*_2'
      len = 4
    End If
    If (kflc>=3 .And. kflr==0 .And. kfls<=3) Then
      chau(len+1:len+2) = '_' // chaf(kflc)(1:1)
      len = len + 2
    Else If (kflc>=3) Then
      chau(len+1:len+1) = chaf(kflc)(1:1)
      len = len + 1
    End If

!...Construct root name and spin for heavy baryon.
  Else
    If (kflb<=2 .And. kflc<=2) Then
      chau = 'Sigma '
      If (kflc>kflb) chau = 'Lambda'
      If (kfls==4) chau = 'Sigma*'
      len = 5
      If (chau(6:6)/=' ') len = 6
    Else If (kflb<=2 .Or. kflc<=2) Then
      chau = 'Xi '
      If (kfla>kflb .And. kflb>kflc) chau = 'Xi'''
      If (kfls==4) chau = 'Xi*'
      len = 2
      If (chau(3:3)/=' ') len = 3
    Else
      chau = 'Omega '
      If (kfla>kflb .And. kflb>kflc) chau = 'Omega'''
      If (kfls==4) chau = 'Omega*'
      len = 5
      If (chau(6:6)/=' ') len = 6
    End If

!...Add on heavy flavour content for heavy baryon.
    chau(len+1:len+2) = '_' // chaf(kfla)(1:1)
    len = len + 2
    If (kflb>=kflc .And. kflc>=4) Then
      chau(len+1:len+2) = chaf(kflb)(1:1) // chaf(kflc)(1:1)
      len = len + 2
    Else If (kflb>=kflc .And. kflb>=4) Then
      chau(len+1:len+1) = chaf(kflb)(1:1)
      len = len + 1
    Else If (kflc>kflb .And. kflb>=4) Then
      chau(len+1:len+2) = chaf(kflc)(1:1) // chaf(kflb)(1:1)
      len = len + 2
    Else If (kflc>kflb .And. kflc>=4) Then
      chau(len+1:len+1) = chaf(kflc)(1:1)
      len = len + 1
    End If
  End If

!...Add on bar sign for antiparticle (where necessary).
  If (kf>0 .Or. len==0) Then
  Else If (kfa>10 .And. kfa<=40 .And. kq/=0) Then
  Else If (kfa==89 .Or. (kfa>=91 .And. kfa<=99)) Then
  Else If (kfa>100 .And. kfla==0 .And. kq/=0) Then
  Else If (mstu(15)<=1) Then
    chau(len+1:len+1) = '~'
    len = len + 1
  Else
    chau(len+1:len+3) = 'bar'
    len = len + 3
  End If

!...Add on charge where applicable (conventional cases skipped).
  If (kq==6) chau(len+1:len+2) = '++'
  If (kq==-6) chau(len+1:len+2) = '--'
  If (kq==3) chau(len+1:len+1) = '+'
  If (kq==-3) chau(len+1:len+1) = '-'
  If (kq==0 .And. (kfa<=22 .Or. len==0)) Then
  Else If (kq==0 .And. (kfa>=81 .And. kfa<=100)) Then
  Else If (kfa>100 .And. kfla==0 .And. kflb==kflc .And. kflb/=1) Then
  Else If (kq==0) Then
    chau(len+1:len+1) = '0'
  End If

  Return
End Subroutine luname

!*********************************************************************

Function luchge(kf)

!...Purpose: to give three times the charge for a particle/parton.
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/

!...Initial values. Simple case of direct readout.
  luchge = 0
  kfa = iabs(kf)
  kc = lucomp(kfa)
  If (kc==0) Then
  Else If (kfa<=100 .Or. kc<=80 .Or. kc>100) Then
    luchge = kchg(kc, 1)

!...Construction from quark content for heavy meson, diquark, baryon.
  Else If (mod(kfa/1000,10)==0) Then
    luchge = (kchg(mod(kfa/100,10),1)-kchg(mod(kfa/10,10),1))*(-1)**mod(kfa/100, 10)
  Else If (mod(kfa/10,10)==0) Then
    luchge = kchg(mod(kfa/1000,10), 1) + kchg(mod(kfa/100,10), 1)
  Else
    luchge = kchg(mod(kfa/1000,10), 1) + kchg(mod(kfa/100,10), 1) + kchg(mod(kfa/10,10), 1)
  End If

!...Add on correct sign.
  luchge = luchge*isign(1, kf)

  Return
End Function luchge

!*********************************************************************

Function lucomp(kf)

!...Purpose: to compress the standard KF codes for use in mass and decay
!...arrays; also to check whether a given code actually is defined.
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/

!...Subdivide KF code into constituent pieces.
  lucomp = 0
  kfa = iabs(kf)
  kfla = mod(kfa/1000, 10)
  kflb = mod(kfa/100, 10)
  kflc = mod(kfa/10, 10)
  kfls = mod(kfa, 10)
  kflr = mod(kfa/10000, 10)

!...Simple cases: direct translation or special codes.
  If (kfa==0 .Or. kfa>=100000) Then
  Else If (kfa<=100) Then
    lucomp = kfa
    If (kf<0 .And. kchg(kfa,3)==0) lucomp = 0
  Else If (kfls==0) Then
    If (kf==130) lucomp = 221
    If (kf==310) lucomp = 222
    If (kfa==210) lucomp = 281
    If (kfa==2110) lucomp = 282
    If (kfa==2210) lucomp = 283

!...Mesons.
  Else If (kfa-10000*kflr<1000) Then
    If (kflb==0 .Or. kflb==9 .Or. kflc==0 .Or. kflc==9) Then
    Else If (kflb<kflc) Then
    Else If (kf<0 .And. kflb==kflc) Then
    Else If (kflb==kflc) Then
      If (kflr==0 .And. kfls==1) Then
        lucomp = 110 + kflb
      Else If (kflr==0 .And. kfls==3) Then
        lucomp = 130 + kflb
      Else If (kflr==1 .And. kfls==3) Then
        lucomp = 150 + kflb
      Else If (kflr==1 .And. kfls==1) Then
        lucomp = 170 + kflb
      Else If (kflr==2 .And. kfls==3) Then
        lucomp = 190 + kflb
      Else If (kflr==0 .And. kfls==5) Then
        lucomp = 210 + kflb
      End If
    Else If (kflb<=5 .And. kflc<=3) Then
      If (kflr==0 .And. kfls==1) Then
        lucomp = 100 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==0 .And. kfls==3) Then
        lucomp = 120 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==1 .And. kfls==3) Then
        lucomp = 140 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==1 .And. kfls==1) Then
        lucomp = 160 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==2 .And. kfls==3) Then
        lucomp = 180 + ((kflb-1)*(kflb-2))/2 + kflc
      Else If (kflr==0 .And. kfls==5) Then
        lucomp = 200 + ((kflb-1)*(kflb-2))/2 + kflc
      End If
    Else If ((kfls==1 .And. kflr<=1) .Or. (kfls==3 .And. kflr<=2) .Or. (kfls==5 .And. kflr==0)) Then
      lucomp = 80 + kflb
    End If

!...Diquarks.
  Else If ((kflr==0 .Or. kflr==1) .And. kflc==0) Then
    If (kfls/=1 .And. kfls/=3) Then
    Else If (kfla==9 .Or. kflb==0 .Or. kflb==9) Then
    Else If (kfla<kflb) Then
    Else If (kfls==1 .And. kfla==kflb) Then
    Else
      lucomp = 90
    End If

!...Spin 1/2 baryons.
  Else If (kflr==0 .And. kfls==2) Then
    If (kfla==9 .Or. kflb==0 .Or. kflb==9 .Or. kflc==9) Then
    Else If (kfla<=kflc .Or. kfla<kflb) Then
    Else If (kfla>=6 .Or. kflb>=4 .Or. kflc>=4) Then
      lucomp = 80 + kfla
    Else If (kflb<kflc) Then
      lucomp = 300 + ((kfla+1)*kfla*(kfla-1))/6 + (kflc*(kflc-1))/2 + kflb
    Else
      lucomp = 330 + ((kfla+1)*kfla*(kfla-1))/6 + (kflb*(kflb-1))/2 + kflc
    End If

!...Spin 3/2 baryons.
  Else If (kflr==0 .And. kfls==4) Then
    If (kfla==9 .Or. kflb==0 .Or. kflb==9 .Or. kflc==9) Then
    Else If (kfla<kflb .Or. kflb<kflc) Then
    Else If (kfla>=6 .Or. kflb>=4) Then
      lucomp = 80 + kfla
    Else
      lucomp = 360 + ((kfla+1)*kfla*(kfla-1))/6 + (kflb*(kflb-1))/2 + kflc
    End If
  End If

  Return
End Function lucomp

!*********************************************************************

Subroutine luerrm(merr, chmess)

!...Purpose: to inform user of errors in program execution.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Character chmess*(*)

  Write (6, *) 'merr,chmess=', merr, chmess

!...Write first few warnings, then be silent.
  If (merr<=10) Then
    mstu(27) = mstu(27) + 1
    mstu(28) = merr
    If (mstu(25)==1 .And. mstu(27)<=mstu(26)) Write (mstu(11), 1000) merr, mstu(31), chmess

!...Write first few errors, then be silent or stop program.
  Else If (merr<=20) Then
    mstu(23) = mstu(23) + 1
    mstu(24) = merr - 10
    If (mstu(21)>=1 .And. mstu(23)<=mstu(22)) Write (mstu(11), 1100) merr - 10, mstu(31), chmess
    If (mstu(21)>=2 .And. mstu(23)>mstu(22)) Then
      Write (mstu(11), 1100) merr - 10, mstu(31), chmess
      Write (mstu(11), 1200)
      If (merr/=17) Call lulist(2)
      Stop
    End If

!...Stop program in case of irreparable error.
  Else
    Write (mstu(11), 1300) merr - 20, mstu(31), chmess
    Stop
  End If

  Return

!...Formats for output.
  1000 Format (/5X, 'Advisory warning type', I2, ' given after', I6, ' LUEXEC calls:'/5X, A)
  1100 Format (/5X, 'Error type', I2, ' has occured after', I6, ' LUEXEC calls:'/5X, A)
  1200 Format (5X, 'Execution will be stopped after listing of last ', 'event!')
  1300 Format (/5X, 'Fatal error type', I2, ' has occured after', I6, ' LUEXEC calls:'/5X, A/5X, 'Execution will now be stopped!')
End Subroutine luerrm

!*********************************************************************

Function ulalps(q2)

!...Purpose: to give the value of alpha_strong.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/

!...Constant alpha_strong trivial.
  If (mstu(111)<=0) Then
    ulalps = paru(111)
    mstu(118) = mstu(112)
    paru(117) = 0.
    paru(118) = paru(111)
    Return
  End If

!...Find effective Q2, number of flavours and Lambda.
  q2eff = q2
  If (mstu(115)>=2) q2eff = max(q2, paru(114))
  nf = mstu(112)
  alam2 = paru(112)**2
  100 If (nf>max(2,mstu(113))) Then
    q2thr = paru(113)*pmas(nf, 1)**2
    If (q2eff<q2thr) Then
      nf = nf - 1
      alam2 = alam2*(q2thr/alam2)**(2./(33.-2.*nf))
      Goto 100
    End If
  End If
  110 If (nf<min(8,mstu(114))) Then
    q2thr = paru(113)*pmas(nf+1, 1)**2
    If (q2eff>q2thr) Then
      nf = nf + 1
      alam2 = alam2*(alam2/q2thr)**(2./(33.-2.*nf))
      Goto 110
    End If
  End If
  If (mstu(115)==1) q2eff = q2eff + alam2
  paru(117) = sqrt(alam2)

!...Evaluate first or second order alpha_strong.
  b0 = (33.-2.*nf)/6.
  algq = log(q2eff/alam2)
  If (mstu(111)==1) Then
    ulalps = paru(2)/(b0*algq)
  Else
    b1 = (153.-19.*nf)/6.
    ulalps = paru(2)/(b0*algq)*(1.-b1*log(algq)/(b0**2*algq))
  End If
  mstu(118) = nf
  paru(118) = ulalps

  Return
End Function ulalps

!*********************************************************************

Function ulangl(x, y)

!...Purpose: to reconstruct an angle from given x and y coordinates.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/

  ulangl = 0.
  r = sqrt(x**2+y**2)
  If (r<1E-20) Return
  If (abs(x)/r<0.8) Then
    ulangl = sign(acos(x/r), y)
  Else
    ulangl = asin(y/r)
    If (x<0. .And. ulangl>=0.) Then
      ulangl = paru(1) - ulangl
    Else If (x<0.) Then
      ulangl = -paru(1) - ulangl
    End If
  End If

  Return
End Function ulangl

!*********************************************************************

Function rlu(idum)

!...Purpose: to generate random numbers uniformly distributed between
!...0 and 1, excluding the endpoints.
  Common /ludatr/mrlu(6), rrlu(100)
  Save /ludatr/
  Equivalence (mrlu1, mrlu(1)), (mrlu2, mrlu(2)), (mrlu3, mrlu(3)), (mrlu4, mrlu(4)), (mrlu5, mrlu(5)), (mrlu6, mrlu(6)), (rrlu98, rrlu(98)), (rrlu99, rrlu(99)), (rrlu00, rrlu(100))

!...Initialize generation from given seed.
  If (mrlu2==0) Then
    ij = mod(mrlu1/30082, 31329)
    kl = mod(mrlu1, 30082)
    i = mod(ij/177, 177) + 2
    j = mod(ij, 177) + 2
    k = mod(kl/169, 178) + 1
    l = mod(kl, 169)
    Do ii = 1, 97
      s = 0.
      t = 0.5
      Do jj = 1, 24
        m = mod(mod(i*j,179)*k, 179)
        i = j
        j = k
        k = m
        l = mod(53*l+1, 169)
        If (mod(l*m,64)>=32) s = s + t
        t = 0.5*t
      End Do
      rrlu(ii) = s
    End Do
    twom24 = 1.
    Do i24 = 1, 24
      twom24 = 0.5*twom24
    End Do
    rrlu98 = 362436.*twom24
    rrlu99 = 7654321.*twom24
    rrlu00 = 16777213.*twom24
    mrlu2 = 1
    mrlu3 = 0
    mrlu4 = 97
    mrlu5 = 33
  End If

!...Generate next random number.
  130 runi = rrlu(mrlu4) - rrlu(mrlu5)
  If (runi<0.) runi = runi + 1.
  rrlu(mrlu4) = runi
  mrlu4 = mrlu4 - 1
  If (mrlu4==0) mrlu4 = 97
  mrlu5 = mrlu5 - 1
  If (mrlu5==0) mrlu5 = 97
  rrlu98 = rrlu98 - rrlu99
  If (rrlu98<0.) rrlu98 = rrlu98 + rrlu00
  runi = runi - rrlu98
  If (runi<0.) runi = runi + 1.
  If (runi<=0 .Or. runi>=1.) Goto 130

!...Update counters. Random number to output.
  mrlu3 = mrlu3 + 1
  If (mrlu3==1000000000) Then
    mrlu2 = mrlu2 + 1
    mrlu3 = 0
  End If
  rlu = runi

  Return
End Function rlu

!*********************************************************************

Subroutine lurobo(the, phi, bex, bey, bez)

!...Purpose: to perform rotations and boosts.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Dimension rot(3, 3), pr(3), vr(3), dp(4), dv(4)

!...Find range of rotation/boost. Convert boost to double precision.
  imin = 1
  If (mstu(1)>0) imin = mstu(1)
  imax = n
  If (mstu(2)>0) imax = mstu(2)
  dbx = dble(bex)
  dby = dble(bey)
  dbz = dble(bez)
  Goto 100

!...Entry for specific range and double precision boost.
  Entry ludbrb(imi, ima, the, phi, dbex, dbey, dbez)
  imin = imi
  If (imin<=0) imin = 1
  imax = ima
  If (imax<=0) imax = n
  dbx = dbex
  dby = dbey
  dbz = dbez

!...Check range of rotation/boost.
  100 If (imin>mstu(4) .Or. imax>mstu(4)) Then
    Call luerrm(11, '(LUROBO:) range outside LUJETS memory')
    Return
  End If

!...Rotate, typically from z axis to direction (theta,phi).
!lin-5/2012:
!      IF(THE**2+PHI**2.GT.1E-20) THEN
  If ((the**2+phi**2)>1E-20) Then
    rot(1, 1) = cos(the)*cos(phi)
    rot(1, 2) = -sin(phi)
    rot(1, 3) = sin(the)*cos(phi)
    rot(2, 1) = cos(the)*sin(phi)
    rot(2, 2) = cos(phi)
    rot(2, 3) = sin(the)*sin(phi)
    rot(3, 1) = -sin(the)
    rot(3, 2) = 0.
    rot(3, 3) = cos(the)
    Do i = imin, imax
      If (k(i,1)<=0) Goto 130
      Do j = 1, 3
        pr(j) = p(i, j)
        vr(j) = v(i, j)
      End Do
      Do j = 1, 3
        p(i, j) = rot(j, 1)*pr(1) + rot(j, 2)*pr(2) + rot(j, 3)*pr(3)
        v(i, j) = rot(j, 1)*vr(1) + rot(j, 2)*vr(2) + rot(j, 3)*vr(3)
      End Do
    130 End Do
  End If

!...Boost, typically from rest to momentum/energy=beta.
!lin-5/2012:
!      IF(DBX**2+DBY**2+DBZ**2.GT.1E-20) THEN
  If ((dbx**2+dby**2+dbz**2)>1D-20) Then
    db = sqrt(dbx**2+dby**2+dbz**2)
    If (db>0.99999999D0) Then
!...Rescale boost vector if too close to unity.
      Call luerrm(3, '(LUROBO:) boost vector too large')
      dbx = dbx*(0.99999999D0/db)
      dby = dby*(0.99999999D0/db)
      dbz = dbz*(0.99999999D0/db)
      db = 0.99999999D0
    End If
    dga = 1D0/sqrt(1D0-db**2)
    Do i = imin, imax
      If (k(i,1)<=0) Goto 150
      Do j = 1, 4
        dp(j) = dble(p(i,j))
        dv(j) = dble(v(i,j))
      End Do
      dbp = dbx*dp(1) + dby*dp(2) + dbz*dp(3)
      dgabp = dga*(dga*dbp/(1D0+dga)+dp(4))
      p(i, 1) = sngl(dp(1)+dgabp*dbx)
      p(i, 2) = sngl(dp(2)+dgabp*dby)
      p(i, 3) = sngl(dp(3)+dgabp*dbz)
      p(i, 4) = sngl(dga*(dp(4)+dbp))
      dbv = dbx*dv(1) + dby*dv(2) + dbz*dv(3)
      dgabv = dga*(dga*dbv/(1D0+dga)+dv(4))
      v(i, 1) = sngl(dv(1)+dgabv*dbx)
      v(i, 2) = sngl(dv(2)+dgabv*dby)
      v(i, 3) = sngl(dv(3)+dgabv*dbz)
      v(i, 4) = sngl(dga*(dv(4)+dbv))
    150 End Do
  End If

  Return
End Subroutine lurobo

!*********************************************************************
! THIS SUBROUTINE IS ONLY FOR THE USE OF HIJING TO ROTATE OR BOOST
!        THE FOUR MOMENTUM ONLY
!*********************************************************************

Subroutine hirobo(the, phi, bex, bey, bez)

!...Purpose: to perform rotations and boosts.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Dimension rot(3, 3), pr(3), vr(3), dp(4), dv(4)

!...Find range of rotation/boost. Convert boost to double precision.
  imin = 1
  If (mstu(1)>0) imin = mstu(1)
  imax = n
  If (mstu(2)>0) imax = mstu(2)
  dbx = dble(bex)
  dby = dble(bey)
  dbz = dble(bez)

!...Check range of rotation/boost.
  If (imin>mstu(4) .Or. imax>mstu(4)) Then
    Call luerrm(11, '(LUROBO:) range outside LUJETS memory')
    Return
  End If

!...Rotate, typically from z axis to direction (theta,phi).
!lin-5/2012:
!      IF(THE**2+PHI**2.GT.1E-20) THEN
  If ((the**2+phi**2)>1E-20) Then
    rot(1, 1) = cos(the)*cos(phi)
    rot(1, 2) = -sin(phi)
    rot(1, 3) = sin(the)*cos(phi)
    rot(2, 1) = cos(the)*sin(phi)
    rot(2, 2) = cos(phi)
    rot(2, 3) = sin(the)*sin(phi)
    rot(3, 1) = -sin(the)
    rot(3, 2) = 0.
    rot(3, 3) = cos(the)
    Do i = imin, imax
      If (k(i,1)<=0) Goto 130
      Do j = 1, 3
        pr(j) = p(i, j)
      End Do
      Do j = 1, 3
        p(i, j) = rot(j, 1)*pr(1) + rot(j, 2)*pr(2) + rot(j, 3)*pr(3)
      End Do
    130 End Do
  End If

!...Boost, typically from rest to momentum/energy=beta.
!lin-5/2012:
!      IF(DBX**2+DBY**2+DBZ**2.GT.1E-20) THEN
  If ((dbx**2+dby**2+dbz**2)>1D-20) Then
    db = sqrt(dbx**2+dby**2+dbz**2)
    If (db>0.99999999D0) Then
!...Rescale boost vector if too close to unity.
      Call luerrm(3, '(LUROBO:) boost vector too large')
      dbx = dbx*(0.99999999D0/db)
      dby = dby*(0.99999999D0/db)
      dbz = dbz*(0.99999999D0/db)
      db = 0.99999999D0
    End If
    dga = 1D0/sqrt(1D0-db**2)
    Do i = imin, imax
      If (k(i,1)<=0) Goto 150
      Do j = 1, 4
        dp(j) = dble(p(i,j))
      End Do
      dbp = dbx*dp(1) + dby*dp(2) + dbz*dp(3)
      dgabp = dga*(dga*dbp/(1D0+dga)+dp(4))
      p(i, 1) = sngl(dp(1)+dgabp*dbx)
      p(i, 2) = sngl(dp(2)+dgabp*dby)
      p(i, 3) = sngl(dp(3)+dgabp*dbz)
      p(i, 4) = sngl(dga*(dp(4)+dbp))
    150 End Do
  End If

  Return
End Subroutine hirobo

!*********************************************************************

Subroutine luedit(medit)

!...Purpose: to perform global manipulations on the event record,
!...in particular to exclude unstable or undetectable partons/particles.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension ns(2), pts(2), pls(2)

!...Remove unwanted partons/particles.
  If ((medit>=0 .And. medit<=3) .Or. medit==5) Then
    imax = n
    If (mstu(2)>0) imax = mstu(2)
    i1 = max(1, mstu(1)) - 1
    Do i = max(1, mstu(1)), imax
      If (k(i,1)==0 .Or. k(i,1)>20) Goto 110
      If (medit==1) Then
        If (k(i,1)>10) Goto 110
      Else If (medit==2) Then
        If (k(i,1)>10) Goto 110
        kc = lucomp(k(i,2))
        If (kc==0 .Or. kc==12 .Or. kc==14 .Or. kc==16 .Or. kc==18) Goto 110
      Else If (medit==3) Then
        If (k(i,1)>10) Goto 110
        kc = lucomp(k(i,2))
        If (kc==0) Goto 110
        If (kchg(kc,2)==0 .And. luchge(k(i,2))==0) Goto 110
      Else If (medit==5) Then
        If (k(i,1)==13 .Or. k(i,1)==14) Goto 110
        kc = lucomp(k(i,2))
        If (kc==0) Goto 110
        If (k(i,1)>=11 .And. kchg(kc,2)==0) Goto 110
      End If

!...Pack remaining partons/particles. Origin no longer known.
      i1 = i1 + 1
      Do j = 1, 5
        k(i1, j) = k(i, j)
        p(i1, j) = p(i, j)
        v(i1, j) = v(i, j)
      End Do
      k(i1, 3) = 0
    110 End Do
    n = i1

!...Selective removal of class of entries. New position of retained.
  Else If (medit>=11 .And. medit<=15) Then
    i1 = 0
    Do i = 1, n
      k(i, 3) = mod(k(i,3), mstu(5))
      If (medit==11 .And. k(i,1)<0) Goto 120
      If (medit==12 .And. k(i,1)==0) Goto 120
      If (medit==13 .And. (k(i,1)==11 .Or. k(i,1)==12 .Or. k(i,1)==15) .And. k(i,2)/=94) Goto 120
      If (medit==14 .And. (k(i,1)==13 .Or. k(i,1)==14 .Or. k(i,2)==94)) Goto 120
      If (medit==15 .And. k(i,1)>=21) Goto 120
      i1 = i1 + 1
      k(i, 3) = k(i, 3) + mstu(5)*i1
    120 End Do

!...Find new event history information and replace old.
    Do i = 1, n
      If (k(i,1)<=0 .Or. k(i,1)>20 .Or. k(i,3)/mstu(5)==0) Goto 140
      id = i
      130 im = mod(k(id,3), mstu(5))
      If (medit==13 .And. im>0 .And. im<=n) Then
        If ((k(im,1)==11 .Or. k(im,1)==12 .Or. k(im,1)==15) .And. k(im,2)/=94) Then
          id = im
          Goto 130
        End If
      Else If (medit==14 .And. im>0 .And. im<=n) Then
        If (k(im,1)==13 .Or. k(im,1)==14 .Or. k(im,2)==94) Then
          id = im
          Goto 130
        End If
      End If
      k(i, 3) = mstu(5)*(k(i,3)/mstu(5))
      If (im/=0) k(i, 3) = k(i, 3) + k(im, 3)/mstu(5)
      If (k(i,1)/=3 .And. k(i,1)/=13 .And. k(i,1)/=14) Then
        If (k(i,4)>0 .And. k(i,4)<=mstu(4)) k(i, 4) = k(k(i,4), 3)/mstu(5)
        If (k(i,5)>0 .And. k(i,5)<=mstu(4)) k(i, 5) = k(k(i,5), 3)/mstu(5)
      Else
        kcm = mod(k(i,4)/mstu(5), mstu(5))
        If (kcm>0 .And. kcm<=mstu(4)) kcm = k(kcm, 3)/mstu(5)
        kcd = mod(k(i,4), mstu(5))
        If (kcd>0 .And. kcd<=mstu(4)) kcd = k(kcd, 3)/mstu(5)
        k(i, 4) = mstu(5)**2*(k(i,4)/mstu(5)**2) + mstu(5)*kcm + kcd
        kcm = mod(k(i,5)/mstu(5), mstu(5))
        If (kcm>0 .And. kcm<=mstu(4)) kcm = k(kcm, 3)/mstu(5)
        kcd = mod(k(i,5), mstu(5))
        If (kcd>0 .And. kcd<=mstu(4)) kcd = k(kcd, 3)/mstu(5)
        k(i, 5) = mstu(5)**2*(k(i,5)/mstu(5)**2) + mstu(5)*kcm + kcd
      End If
    140 End Do

!...Pack remaining entries.
    i1 = 0
    Do i = 1, n
      If (k(i,3)/mstu(5)==0) Goto 160
      i1 = i1 + 1
      Do j = 1, 5
        k(i1, j) = k(i, j)
        p(i1, j) = p(i, j)
        v(i1, j) = v(i, j)
      End Do
      k(i1, 3) = mod(k(i1,3), mstu(5))
    160 End Do
    n = i1

!...Save top entries at bottom of LUJETS commonblock.
  Else If (medit==21) Then
    If (2*n>=mstu(4)) Then
      Call luerrm(11, '(LUEDIT:) no more memory left in LUJETS')
      Return
    End If
    Do i = 1, n
      Do j = 1, 5
        k(mstu(4)-i, j) = k(i, j)
        p(mstu(4)-i, j) = p(i, j)
        v(mstu(4)-i, j) = v(i, j)
      End Do
    End Do
    mstu(32) = n

!...Restore bottom entries of commonblock LUJETS to top.
  Else If (medit==22) Then
    Do i = 1, mstu(32)
      Do j = 1, 5
        k(i, j) = k(mstu(4)-i, j)
        p(i, j) = p(mstu(4)-i, j)
        v(i, j) = v(mstu(4)-i, j)
      End Do
    End Do
    n = mstu(32)

!...Mark primary entries at top of commonblock LUJETS as untreated.
  Else If (medit==23) Then
    i1 = 0
    Do i = 1, n
      kh = k(i, 3)
      If (kh>=1) Then
        If (k(kh,1)>20) kh = 0
      End If
      If (kh/=0) Goto 200
      i1 = i1 + 1
      If (k(i,1)>10 .And. k(i,1)<=20) k(i, 1) = k(i, 1) - 10
    End Do
    200 n = i1

!...Place largest axis along z axis and second largest in xy plane.
  Else If (medit==31 .Or. medit==32) Then
    Call ludbrb(1, n+mstu(3), 0., -ulangl(p(mstu(61),1),p(mstu(61),2)), 0D0, 0D0, 0D0)
    Call ludbrb(1, n+mstu(3), -ulangl(p(mstu(61),3),p(mstu(61),1)), 0., 0D0, 0D0, 0D0)
    Call ludbrb(1, n+mstu(3), 0., -ulangl(p(mstu(61)+1,1),p(mstu(61)+1,2)), 0D0, 0D0, 0D0)
    If (medit==31) Return

!...Rotate to put slim jet along +z axis.
    Do is = 1, 2
      ns(is) = 0
      pts(is) = 0.
      pls(is) = 0.
    End Do
    Do i = 1, n
      If (k(i,1)<=0 .Or. k(i,1)>10) Goto 220
      If (mstu(41)>=2) Then
        kc = lucomp(k(i,2))
        If (kc==0 .Or. kc==12 .Or. kc==14 .Or. kc==16 .Or. kc==18) Goto 220
        If (mstu(41)>=3 .And. kchg(kc,2)==0 .And. luchge(k(i,2))==0) Goto 220
      End If
      is = int(2.-sign(0.5,p(i,3)))
      ns(is) = ns(is) + 1
      pts(is) = pts(is) + sqrt(p(i,1)**2+p(i,2)**2)
    220 End Do
    If (ns(1)*pts(2)**2<ns(2)*pts(1)**2) Call ludbrb(1, n+mstu(3), paru(1), 0., 0D0, 0D0, 0D0)

!...Rotate to put second largest jet into -z,+x quadrant.
    Do i = 1, n
      If (p(i,3)>=0.) Goto 230
      If (k(i,1)<=0 .Or. k(i,1)>10) Goto 230
      If (mstu(41)>=2) Then
        kc = lucomp(k(i,2))
        If (kc==0 .Or. kc==12 .Or. kc==14 .Or. kc==16 .Or. kc==18) Goto 230
        If (mstu(41)>=3 .And. kchg(kc,2)==0 .And. luchge(k(i,2))==0) Goto 230
      End If
      is = int(2.-sign(0.5,p(i,1)))
      pls(is) = pls(is) - p(i, 3)
    230 End Do
    If (pls(2)>pls(1)) Call ludbrb(1, n+mstu(3), 0., paru(1), 0D0, 0D0, 0D0)
  End If

  Return
End Subroutine luedit

!*********************************************************************

Subroutine lulist(mlist)

!...Purpose: to give program heading, or list an event, or particle
!...data, or current parameter values.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Character chap*16, chac*16, chan*16, chad(5)*16, chmo(12)*3, chdl(7)*4
  Dimension ps(6)
  Data chmo/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/, chdl/'(())', ' ', '()', '!!', '<>', '==', '(==)'/

!...Initialization printout: version number and date of last change.
!      IF(MLIST.EQ.0.OR.MSTU(12).EQ.1) THEN
!        WRITE(MSTU(11),1000) MSTU(181),MSTU(182),MSTU(185),
!     &  CHMO(MSTU(184)),MSTU(183)
!        MSTU(12)=0
!        IF(MLIST.EQ.0) RETURN
!      ENDIF

!...List event data, including additional lines after N.
  If (mlist>=1 .And. mlist<=3) Then
    If (mlist==1) Write (mstu(11), 1100)
    If (mlist==2) Write (mstu(11), 1200)
    If (mlist==3) Write (mstu(11), 1300)
    lmx = 12
    If (mlist>=2) lmx = 16
    istr = 0
    imax = n
    If (mstu(2)>0) imax = mstu(2)
    Do i = max(1, mstu(1)), max(imax, n+max(0,mstu(3)))
      If ((i>imax .And. i<=n) .Or. k(i,1)<0) Goto 120

!...Get particle name, pad it and check it is not too long.
      Call luname(k(i,2), chap)
      len = 0
      Do lem = 1, 16
        If (chap(lem:lem)/=' ') len = lem
      End Do
      mdl = (k(i,1)+19)/10
      ldl = 0
      If (mdl==2 .Or. mdl>=8) Then
        chac = chap
        If (len>lmx) chac(lmx:lmx) = '?'
      Else
        ldl = 1
        If (mdl==1 .Or. mdl==7) ldl = 2
        If (len==0) Then
          chac = chdl(mdl)(1:2*ldl) // ' '
        Else
          chac = chdl(mdl)(1:ldl) // chap(1:min(len,lmx-2*ldl)) // chdl(mdl)(ldl+1:2*ldl) // ' '
          If (len+2*ldl>lmx) chac(lmx:lmx) = '?'
        End If
      End If

!...Add information on string connection.
      If (k(i,1)==1 .Or. k(i,1)==2 .Or. k(i,1)==11 .Or. k(i,1)==12) Then
        kc = lucomp(k(i,2))
        kcc = 0
        If (kc/=0) kcc = kchg(kc, 2)
        If (kcc/=0 .And. istr==0) Then
          istr = 1
          If (len+2*ldl+3<=lmx) chac(lmx-1:lmx-1) = 'A'
        Else If (kcc/=0 .And. (k(i,1)==2 .Or. k(i,1)==12)) Then
          If (len+2*ldl+3<=lmx) chac(lmx-1:lmx-1) = 'I'
        Else If (kcc/=0) Then
          istr = 0
          If (len+2*ldl+3<=lmx) chac(lmx-1:lmx-1) = 'V'
        End If
      End If

!...Write data for particle/jet.
      If (mlist==1 .And. abs(p(i,4))<9999.) Then
        Write (mstu(11), 1400) i, chac(1:12), (k(i,j1), j1=1, 3), (p(i,j2), j2=1, 5)
      Else If (mlist==1 .And. abs(p(i,4))<99999.) Then
        Write (mstu(11), 1500) i, chac(1:12), (k(i,j1), j1=1, 3), (p(i,j2), j2=1, 5)
      Else If (mlist==1) Then
        Write (mstu(11), 1600) i, chac(1:12), (k(i,j1), j1=1, 3), (p(i,j2), j2=1, 5)
      Else If (mstu(5)==10000 .And. (k(i,1)==3 .Or. k(i,1)==13 .Or. k(i,1)==14)) Then
        Write (mstu(11), 1700) i, chac, (k(i,j1), j1=1, 3), k(i, 4)/100000000, mod(k(i,4)/10000, 10000), mod(k(i,4), 10000), k(i, 5)/100000000, mod(k(i,5)/10000, 10000), mod(k(i,5), 10000), (p(i,j2), j2=1, 5)
      Else
        Write (mstu(11), 1800) i, chac, (k(i,j1), j1=1, 5), (p(i,j2), j2=1, 5)
      End If
      If (mlist==3) Write (mstu(11), 1900)(v(i,j), j=1, 5)

!...Insert extra separator lines specified by user.
      If (mstu(70)>=1) Then
        isep = 0
        Do j = 1, min(10, mstu(70))
          If (i==mstu(70+j)) isep = 1
        End Do
        If (isep==1 .And. mlist==1) Write (mstu(11), 2000)
        If (isep==1 .And. mlist>=2) Write (mstu(11), 2100)
      End If
    120 End Do

!...Sum of charges and momenta.
    Do j = 1, 6
      ps(j) = plu(0, j)
    End Do
    If (mlist==1 .And. abs(ps(4))<9999.) Then
      Write (mstu(11), 2200) ps(6), (ps(j), j=1, 5)
    Else If (mlist==1 .And. abs(ps(4))<99999.) Then
      Write (mstu(11), 2300) ps(6), (ps(j), j=1, 5)
    Else If (mlist==1) Then
      Write (mstu(11), 2400) ps(6), (ps(j), j=1, 5)
    Else
      Write (mstu(11), 2500) ps(6), (ps(j), j=1, 5)
    End If

!...Give simple list of KF codes defined in program.
  Else If (mlist==11) Then
    Write (mstu(11), 2600)
    Do kf = 1, 40
      Call luname(kf, chap)
      Call luname(-kf, chan)
      If (chap/=' ' .And. chan==' ') Write (mstu(11), 2700) kf, chap
      If (chan/=' ') Write (mstu(11), 2700) kf, chap, -kf, chan
    End Do
    Do kfls = 1, 3, 2
      Do kfla = 1, 8
        Do kflb = 1, kfla - (3-kfls)/2
          kf = 1000*kfla + 100*kflb + kfls
          Call luname(kf, chap)
          Call luname(-kf, chan)
          Write (mstu(11), 2700) kf, chap, -kf, chan
        End Do
      End Do
    End Do
    Do kmul = 0, 5
      kfls = 3
      If (kmul==0 .Or. kmul==3) kfls = 1
      If (kmul==5) kfls = 5
      kflr = 0
      If (kmul==2 .Or. kmul==3) kflr = 1
      If (kmul==4) kflr = 2
      Do kflb = 1, 8
        Do kflc = 1, kflb - 1
          kf = 10000*kflr + 100*kflb + 10*kflc + kfls
          Call luname(kf, chap)
          Call luname(-kf, chan)
          Write (mstu(11), 2700) kf, chap, -kf, chan
        End Do
        kf = 10000*kflr + 110*kflb + kfls
        Call luname(kf, chap)
        Write (mstu(11), 2700) kf, chap
      End Do
    End Do
    kf = 130
    Call luname(kf, chap)
    Write (mstu(11), 2700) kf, chap
    kf = 310
    Call luname(kf, chap)
    Write (mstu(11), 2700) kf, chap
    Do kflsp = 1, 3
      kfls = 2 + 2*(kflsp/3)
      Do kfla = 1, 8
        Do kflb = 1, kfla
          Do kflc = 1, kflb
            If (kflsp==1 .And. (kfla==kflb .Or. kflb==kflc)) Goto 180
            If (kflsp==2 .And. kfla==kflc) Goto 180
            If (kflsp==1) kf = 1000*kfla + 100*kflc + 10*kflb + kfls
            If (kflsp>=2) kf = 1000*kfla + 100*kflb + 10*kflc + kfls
            Call luname(kf, chap)
            Call luname(-kf, chan)
            Write (mstu(11), 2700) kf, chap, -kf, chan
          180 End Do
        End Do
      End Do
    End Do

!...List parton/particle data table. Check whether to be listed.
  Else If (mlist==12) Then
    Write (mstu(11), 2800)
    mstj24 = mstj(24)
    mstj(24) = 0
    kfmax = 20883
    If (mstu(2)/=0) kfmax = mstu(2)
    Do kf = max(1, mstu(1)), kfmax
      kc = lucomp(kf)
      If (kc==0) Goto 220
      If (mstu(14)==0 .And. kf>100 .And. kc<=100) Goto 220
      If (mstu(14)>0 .And. kf>100 .And. max(mod(kf/1000,10),mod(kf/100,10))>mstu(14)) Goto 220

!...Find particle name and mass. Print information.
      Call luname(kf, chap)
      If (kf<=100 .And. chap==' ' .And. mdcy(kc,2)==0) Goto 220
      Call luname(-kf, chan)
      pm = ulmass(kf)
      Write (mstu(11), 2900) kf, kc, chap, chan, kchg(kc, 1), kchg(kc, 2), kchg(kc, 3), pm, pmas(kc, 2), pmas(kc, 3), pmas(kc, 4), mdcy(kc, 1)

!...Particle decay: channel number, branching ration, matrix element,
!...decay products.
      If (kf>100 .And. kc<=100) Goto 220
      Do idc = mdcy(kc, 2), mdcy(kc, 2) + mdcy(kc, 3) - 1
        Do j = 1, 5
          Call luname(kfdp(idc,j), chad(j))
        End Do
        Write (mstu(11), 3000) idc, mdme(idc, 1), mdme(idc, 2), brat(idc), (chad(j), j=1, 5)
      End Do
    220 End Do
    mstj(24) = mstj24

!...List parameter value table.
  Else If (mlist==13) Then
    Write (mstu(11), 3100)
    Do i = 1, 200
      Write (mstu(11), 3200) i, mstu(i), paru(i), mstj(i), parj(i), parf(i)
    End Do
  End If

  Return

!...Format statements for output on unit MSTU(11) (by default 6).
!lin 1000 FORMAT(///20X,'The Lund Monte Carlo - JETSET version ',I1,'.',I1/
!lin     &20X,'**  Last date of change:  ',I2,1X,A3,1X,I4,'  **'/)
  1100 Format (///28X, 'Event listing (summary)'//4X, 'I  particle/jet KS', 5X, 'KF orig    p_x      p_y      p_z       E        m'/)
  1200 Format (///28X, 'Event listing (standard)'//4X, 'I  particle/jet', '  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', '       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/)
  1300 Format (///28X, 'Event listing (with vertices)'//4X, 'I  particle/j', 'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)', '       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/73X, 'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)'/)
  1400 Format (1X, I4, 2X, A12, 1X, I2, 1X, I6, 1X, I4, 5F9.3)
  1500 Format (1X, I4, 2X, A12, 1X, I2, 1X, I6, 1X, I4, 5F9.2)
  1600 Format (1X, I4, 2X, A12, 1X, I2, 1X, I6, 1X, I4, 5F9.1)
  1700 Format (1X, I4, 2X, A16, 1X, I3, 1X, I8, 2X, I4, 2(3X,I1,2I4), 5F13.5)
  1800 Format (1X, I4, 2X, A16, 1X, I3, 1X, I8, 2X, I4, 2(3X,I9), 5F13.5)
  1900 Format (66X, 5(1X,F12.3))
  2000 Format (1X, 78('='))
  2100 Format (1X, 130('='))
  2200 Format (19X, 'sum:', F6.2, 5X, 5F9.3)
  2300 Format (19X, 'sum:', F6.2, 5X, 5F9.2)
  2400 Format (19X, 'sum:', F6.2, 5X, 5F9.1)
  2500 Format (19X, 'sum charge:', F6.2, 3X, 'sum momentum and inv. mass:', 5F13.5)
  2600 Format (///20X, 'List of KF codes in program'/)
  2700 Format (4X, I6, 4X, A16, 6X, I6, 4X, A16)
  2800 Format (///30X, 'Particle/parton data table'//5X, 'KF', 5X, 'KC', 4X, 'particle', 8X, 'antiparticle', 6X, 'chg  col  anti', 8X, 'mass', 7X, 'width', 7X, 'w-cut', 5X, 'lifetime', 1X, 'decay'/11X, 'IDC', 1X, 'on/off', 1X, 'ME', 3X, 'Br.rat.', 4X, 'decay products')
  2900 Format (/1X, I6, 3X, I4, 4X, A16, A16, 3I5, 1X, F12.5, 2(1X,F11.5), 2X, F12.5, 3X, I2)
  3000 Format (10X, I4, 2X, I3, 2X, I3, 2X, F8.5, 4X, 5A16)
  3100 Format (///20X, 'Parameter value table'//4X, 'I', 3X, 'MSTU(I)', 8X, 'PARU(I)', 3X, 'MSTJ(I)', 8X, 'PARJ(I)', 8X, 'PARF(I)')
  3200 Format (1X, I4, 1X, I9, 1X, F14.5, 1X, I9, 1X, F14.5, 1X, F14.5)
End Subroutine lulist

!*********************************************************************

Function plu(i, j)

!...Purpose: to provide various real-valued event related data.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Dimension psum(4)

!...Set default value. For I = 0 sum of momenta or charges,
!...or invariant mass of system.
  plu = 0.
  If (i<0 .Or. i>mstu(4) .Or. j<=0) Then
  Else If (i==0 .And. j<=4) Then
    Do i1 = 1, n
      If (k(i1,1)>0 .And. k(i1,1)<=10) plu = plu + p(i1, j)
    End Do
  Else If (i==0 .And. j==5) Then
    Do j1 = 1, 4
      psum(j1) = 0.
      Do i1 = 1, n
        If (k(i1,1)>0 .And. k(i1,1)<=10) psum(j1) = psum(j1) + p(i1, j1)
      End Do
    End Do
    plu = sqrt(max(0.,psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2))
  Else If (i==0 .And. j==6) Then
    Do i1 = 1, n
      If (k(i1,1)>0 .And. k(i1,1)<=10) plu = plu + luchge(k(i1,2))/3.
    End Do
  Else If (i==0) Then

!...Direct readout of P matrix.
  Else If (j<=5) Then
    plu = p(i, j)

!...Charge, total momentum, transverse momentum, transverse mass.
  Else If (j<=12) Then
    If (j==6) plu = luchge(k(i,2))/3.
    If (j==7 .Or. j==8) plu = p(i, 1)**2 + p(i, 2)**2 + p(i, 3)**2
    If (j==9 .Or. j==10) plu = p(i, 1)**2 + p(i, 2)**2
    If (j==11 .Or. j==12) plu = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
    If (j==8 .Or. j==10 .Or. j==12) plu = sqrt(plu)

!...Theta and phi angle in radians or degrees.
  Else If (j<=16) Then
    If (j<=14) plu = ulangl(p(i,3), sqrt(p(i,1)**2+p(i,2)**2))
    If (j>=15) plu = ulangl(p(i,1), p(i,2))
    If (j==14 .Or. j==16) plu = plu*180./paru(1)

!...True rapidity, rapidity with pion mass, pseudorapidity.
  Else If (j<=19) Then
    pmr = 0.
    If (j==17) pmr = p(i, 5)
    If (j==18) pmr = ulmass(211)
    pr = max(1E-20, pmr**2+p(i,1)**2+p(i,2)**2)
    plu = sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),1E20)), p(i,3))

!...Energy and momentum fractions (only to be used in CM frame).
  Else If (j<=25) Then
    If (j==20) plu = 2.*sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)/paru(21)
    If (j==21) plu = 2.*p(i, 3)/paru(21)
    If (j==22) plu = 2.*sqrt(p(i,1)**2+p(i,2)**2)/paru(21)
    If (j==23) plu = 2.*p(i, 4)/paru(21)
    If (j==24) plu = (p(i,4)+p(i,3))/paru(21)
    If (j==25) plu = (p(i,4)-p(i,3))/paru(21)
  End If

  Return
End Function plu

!*********************************************************************

Block Data ludata

!...Purpose: to give default values to parameters and particle and
!...decay data.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /ludat4/chaf(500)
  Character chaf*8
  Save /ludat4/
  Common /ludatr/mrlu(6), rrlu(100)
  Save /ludatr/

!...LUDAT1, containing status codes and most parameters.
  Data mstu/0, 0, 0, 9000, 10000, 500, 2000, 0, 0, 2, 6, 1, 1, 0, 1, 1, 0, 0, 0, 0, 2, 10, 0, 0, 1, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 4, 2, 1, 1, 0, 0, 0, 25, 24, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40*0, 1, 5, 3, 5, 0, 0, 0, 0, 0, 0, 60*0, 7, 2, 1989, 11, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  Data paru/3.1415927, 6.2831854, 0.1973, 5.068, 0.3894, 2.568, 4*0., 0.001, 0.09, 0.01, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.0, 1.0, 0.25, 2.5, 0.05, 0., 0., 0.0001, 0., 0., 2.5, 1.5, 7.0, 1.0, 0.5, 2.0, 3.2, 0., 0., 0., 40*0., 0.0072974, 0.230, 0., 0., 0., 0., 0., 0., 0., 0., 0.20, 0.25, 1.0, 4.0, 0., 0., 0., 0., 0., 0., 1.0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 70*0./
  Data mstj/1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 2, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4, 2, 5, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 40*0, 5, 2, 7, 5, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 80*0/
  Data parj/0.10, 0.30, 0.40, 0.05, 0.50, 0.50, 0.50, 0., 0., 0., 0.50, 0.60, 0.75, 0., 0., 0., 0., 1.0, 1.0, 0., 0.35, 1.0, 0., 0., 0., 0., 0., 0., 0., 0., 0.10, 1.0, 0.8, 1.5, 0.8, 2.0, 0.2, 2.5, 0.6, 2.5, 0.5, 0.9, 0.5, 0.9, 0.5, 0., 0., 0., 0., 0., 0.77, 0.77, 0.77, 0., 0., 0., 0., 0., 1.0, 0., 4.5, 0.7, 0., 0.003, 0.5, 0.5, 0., 0., 0., 0., 10., 1000., 100., 1000., 0., 0., 0., 0., 0., 0., 0.4, 1.0, 1.0, 0., 10., 10., 0., 0., 0., 0., 0.02, 1.0, 0.2, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.5, 0.5, 91.2, 2.40, 0.02, 2.0, 1.0, 0.25, 0.002, 0., 0., 0., 0., 0., 0.01, 0.99, 0., 0., 0.2, 0., 60*0./

!...LUDAT2, with particle data and flavour treatment parameters.
  Data (kchg(i,1), i=1, 500)/ -1, 2, -1, 2, -1, 2, -1, 2, 2*0, -3, 0, -3, 0, -3, 0, -3, 6*0, 3, 9*0, 3, 2*0, 3, 46*0, 2, -1, 2, -1, 2, 3, 11*0, 3, 0, 2*3, 0, 3, 0, 3, 12*0, 3, 0, 2*3, 0, 3, 0, 3, 12*0, 3, 0, 2*3, 0, 3, 0, 3, 12*0, 3, 0, 2*3, 0, 3, 0, 3, 12*0, 3, 0, 2*3, 0, 3, 0, 3, 12*0, 3, 0, 2*3, 0, 3, 0, 3, 72*0, 3, 0, 3, 28*0, 3, 2*0, 3, 8*0, -3, 8*0, 3, 0, -3, 0, 3, -3, 3*0, 3, 6, 0, 3, 5*0, -3, 0, 3, -3, 0, -3, 4*0, -3, 0, 3, 6, -3, 0, 3, -3, 0, -3, 0, 3, 6, 0, 3, 5*0, -3, 0, 3, -3, 0, -3, 114*0/
  Data (kchg(i,2), i=1, 500)/8*1, 12*0, 2, 68*0, -1, 410*0/
  Data (kchg(i,3), i=1, 500)/8*1, 2*0, 8*1, 5*0, 1, 9*0, 1, 2*0, 1, 2*0, 1, 41*0, 1, 0, 7*1, 10*0, 9*1, 11*0, 9*1, 11*0, 9*1, 11*0, 9*1, 11*0, 9*1, 11*0, 9*1, 71*0, 3*1, 22*0, 1, 5*0, 1, 0, 2*1, 6*0, 1, 0, 2*1, 6*0, 2*1, 0, 5*1, 0, 6*1, 4*0, 6*1, 4*0, 16*1, 4*0, 6*1, 114*0/
  Data (pmas(i,1), i=1, 500)/.0099, .0056, .199, 1.35, 5., 90., 120., 200., 2*0., .00051, 0., .1057, 0., 1.7841, 0., 60., 5*0., 91.2, 80., 15., 6*0., 300., 900., 600., 300., 900., 300., 2*0., 5000., 60*0., .1396, .4977, .4936, 1.8693, 1.8645, 1.9693, 5.2794, 5.2776, 5.47972, 0., .135, .5488, .9575, 2.9796, 9.4, 117.99, 238., 397., 2*0., .7669, .8962, .8921, 2.0101, 2.0071, 2.1127, 2*5.3354, 5.5068, 0., .77, .782, 1.0194, 3.0969, 9.4603, 118., 238., 397., 2*0., 1.233, 2*1.3, 2*2.322, 2.51, 2*5.73, 5.97, 0., 1.233, 1.17, 1.41, 3.46, 9.875, 118.42, 238.42, 397.42, 2*0., .983, 2*1.429, 2*2.272, 2.46, 2*5.68, 5.92, 0., .983, 1., 1.4, 3.4151, 9.8598, 118.4, 238.4, 397.4, 2*0., 1.26, 2*1.401, 2*2.372, 2.56, 2*5.78, 6.02, 0., 1.26, 1.283, 1.422, 3.5106, 9.8919, 118.5, 238.5, 397.5, 2*0., 1.318, 2*1.426, 2*2.422, 2.61, 2*5.83, 6.07, 0., 1.318, 1.274, 1.525, 3.5563, 9.9132, 118.45, 238.45, 397.45, 2*0., 2*.4977, 83*0., 1.1156, 5*0., 2.2849, 0., 2*2.46, 6*0., 5.62, 0., 2*5.84, 6*0., .9396, .9383, 0., 1.1974, 1.1926, &
    1.1894, 1.3213, 1.3149, 0., 2.454, 2.4529, 2.4522, 2*2.55, 2.73, 4*0., 3*5.8, 2*5.96, 6.12, 4*0., 1.234, 1.233, 1.232, 1.231, 1.3872, 1.3837, 1.3828, 1.535, 1.5318, 1.6724, 3*2.5, 2*2.63, 2.8, 4*0., 3*5.81, 2*5.97, 6.13, 114*0./
  Data (pmas(i,2), i=1, 500)/22*0., 2.4, 2.3, 88*0., .0002, .001, 6*0., .149, .0505, .0513, 7*0., .153, .0085, .0044, 7*0., .15, 2*.09, 2*.06, .04, 3*.1, 0., .15, .335, .08, 2*.01, 5*0., .057, 2*.287, 2*.06, .04, 3*.1, 0., .057, 0., .25, .0135, 6*0., .4, 2*.184, 2*.06, .04, 3*.1, 0., .4, .025, .055, .0135, 6*0., .11, .115, .099, 2*.06, 4*.1, 0., .11, .185, .076, .0026, 146*0., 4*.115, .039, 2*.036, .0099, .0091, 131*0./
  Data (pmas(i,3), i=1, 500)/22*0., 2*20., 88*0., .002, .005, 6*0., .4, 2*.2, 7*0., .4, .1, .015, 7*0., .25, 2*.01, 3*.08, 2*.2, .12, 0., .25, .2, .001, 2*.02, 5*0., .05, 2*.4, 3*.08, 2*.2, .12, 0., .05, 0., .35, .05, 6*0., 3*.3, 2*.08, .06, 2*.2, .12, 0., .3, .05, .025, .001, 6*0., .25, 4*.12, 4*.2, 0., .25, .17, .2, .01, 146*0., 4*.14, .04, 2*.035, 2*.05, 131*0./
  Data (pmas(i,4), i=1, 500)/12*0., 658650., 0., .091, 68*0., .1, .43, 15*0., 7803., 0., 3709., .32, .128, .131, 3*.393, 84*0., .004, 26*0., 15540., 26.75, 83*0., 78.88, 5*0., .054, 0., 2*.13, 6*0., .393, 0., 2*.393, 9*0., 44.3, 0., 24., 49.1, 86.9, 6*0., .13, 9*0., .393, 13*0., 24.6, 130*0./
  Data parf/0.5, 0.25, 0.5, 0.25, 1., 0.5, 0., 0., 0., 0., 0.5, 0., 0.5, 0., 1., 1., 0., 0., 0., 0., 0.5, 0., 0.5, 0., 1., 1., 0., 0., 0., 0., 0.5, 0., 0.5, 0., 1., 1., 0., 0., 0., 0., 0.5, 0., 0.5, 0., 1., 1., 0., 0., 0., 0., 0.5, 0., 0.5, 0., 1., 1., 0., 0., 0., 0., 0.75, 0.5, 0., 0.1667, 0.0833, 0.1667, 0., 0., 0., 0., 0., 0., 1., 0.3333, 0.6667, 0.3333, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.325, 0.325, 0.5, 1.6, 5.0, 0., 0., 0., 0., 0., 0., 0.11, 0.16, 0.048, 0.50, 0.45, 0.55, 0.60, 0., 0., 0.2, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 1870*0./
  Data ((vckm(i,j),j=1,4), i=1, 4)/0.95150, 0.04847, 0.00003, 0.00000, 0.04847, 0.94936, 0.00217, 0.00000, 0.00003, 0.00217, 0.99780, 0.00000, 0.00000, 0.00000, 0.00000, 1.00000/

!...LUDAT3, with particle decay parameters and data.
  Data (mdcy(i,1), i=1, 500)/14*0, 1, 0, 1, 5*0, 3*1, 6*0, 1, 4*0, 1, 2*0, 1, 42*0, 7*1, 12*0, 1, 0, 6*1, 0, 8*1, 2*0, 9*1, 0, 8*1, 2*0, 9*1, 0, 8*1, 2*0, 9*1, 0, 8*1, 2*0, 9*1, 0, 8*1, 2*0, 9*1, 0, 8*1, 3*0, 1, 83*0, 1, 5*0, 1, 0, 2*1, 6*0, 1, 0, 2*1, 9*0, 5*1, 0, 6*1, 4*0, 6*1, 4*0, 16*1, 4*0, 6*1, 114*0/
  Data (mdcy(i,2), i=1, 500)/1, 9, 17, 25, 33, 41, 49, 57, 2*0, 65, 69, 71, 76, 78, 118, 120, 125, 2*0, 127, 136, 149, 166, 186, 6*0, 203, 4*0, 219, 2*0, 227, 42*0, 236, 237, 241, 250, 252, 254, 256, 11*0, 276, 277, 279, 285, 406, 574, 606, 607, 608, 0, 609, 611, 617, 623, 624, 625, 626, 627, 2*0, 628, 629, 632, 635, 638, 640, 641, 642, 643, 0, 644, 645, 650, 658, 661, 670, 685, 686, 2*0, 687, 688, 693, 698, 700, 702, 703, 705, 707, 0, 709, 710, 713, 717, 718, 719, 721, 722, 2*0, 723, 726, 728, 730, 734, 738, 740, 744, 748, 0, 752, 755, 759, 763, 765, 767, 769, 770, 2*0, 771, 773, 775, 777, 779, 781, 784, 786, 788, 0, 791, 793, 806, 810, 812, 814, 816, 817, 2*0, 818, 824, 835, 846, 854, 862, 867, 875, 883, 0, 888, 895, 903, 905, 907, 909, 911, 912, 2*0, 913, 921, 83*0, 923, 5*0, 927, 0, 1001, 1002, 6*0, 1003, 0, 1004, 1005, 9*0, 1006, 1008, 1009, 1012, 1013, 0, 1015, 1016, 1017, 1018, 1019, 1020, 4*0, 1021, 1022, 1023, 1024, 1025, 1026, 4*0, 1027, 1028, 1031, 1034, 1035, 1038, 1041, 1044, 1046, 1048, 1052, &
    1053, 1054, 1055, 1057, 1059, 4*0, 1060, 1061, 1062, 1063, 1064, 1065, 114*0/
  Data (mdcy(i,3), i=1, 500)/8*8, 2*0, 4, 2, 5, 2, 40, 2, 5, 2, 2*0, 9, 13, 17, 20, 17, 6*0, 16, 4*0, 8, 2*0, 9, 42*0, 1, 4, 9, 3*2, 20, 11*0, 1, 2, 6, 121, 168, 32, 3*1, 0, 2, 2*6, 5*1, 2*0, 1, 3*3, 2, 4*1, 0, 1, 5, 8, 3, 9, 15, 2*1, 2*0, 1, 2*5, 2*2, 1, 3*2, 0, 1, 3, 4, 2*1, 2, 2*1, 2*0, 3, 2*2, 2*4, 2, 3*4, 0, 3, 2*4, 3*2, 2*1, 2*0, 5*2, 3, 2*2, 3, 0, 2, 13, 4, 3*2, 2*1, 2*0, 6, 2*11, 2*8, 5, 2*8, 5, 0, 7, 8, 4*2, 2*1, 2*0, 8, 2, 83*0, 4, 5*0, 74, 0, 2*1, 6*0, 1, 0, 2*1, 9*0, 2, 1, 3, 1, 2, 0, 6*1, 4*0, 6*1, 4*0, 1, 2*3, 1, 3*3, 2*2, 4, 3*1, 2*2, 1, 4*0, 6*1, 114*0/
  Data (mdme(i,1), i=1, 2000)/6*1, -1, 7*1, -1, 7*1, -1, 7*1, -1, 7*1, -1, 7*1, -1, 85*1, 2* -1, 7*1, 2* -1, 3*1, 2* -1, 6*1, 2* -1, 6*1, 3* -1, 3*1, -1, 3*1, -1, 3*1, 5* -1, 3*1, -1, 6*1, 2* -1, 3*1, -1, 11*1, 2* -1, 6*1, 2* -1, 3*1, -1, 3*1, -1, 4*1, 2* -1, 2*1, -1, 488*1, 2*0, 1275*1/
  Data (mdme(i,2), i=1, 2000)/70*102, 42, 6*102, 2*42, 2*0, 7*41, 2*0, 23*41, 6*102, 45, 28*102, 8*32, 9*0, 16*32, 4*0, 8*32, 4*0, 32, 4*0, 8*32, 8*0, 4*32, 4*0, 6*32, 3*0, 12, 2*42, 2*11, 9*42, 6*45, 20*46, 7*0, 34*42, 86*0, 2*25, 26, 24*42, 142*0, 25, 26, 0, 10*42, 19*0, 2*13, 3*85, 0, 2, 4*0, 2, 8*0, 2*32, 87, 88, 3*3, 0, 2*3, 0, 2*3, 0, 3, 5*0, 3, 1, 0, 3, 2*0, 2*3, 3*0, 1, 4*0, 12, 3*0, 4*32, 2*4, 6*0, 5*32, 2*4, 2*45, 87, 88, 30*0, 12, 32, 0, 32, 87, 88, 41*0, 12, 0, 32, 0, 32, 87, 88, 40*0, 12, 0, 32, 0, 32, 87, 88, 88*0, 12, 0, 32, 0, 32, 87, 88, 2*0, 4*42, 8*0, 14*42, 50*0, 10*13, 2*84, 3*85, 14*0, 84, 5*0, 85, 974*0/
  Data (brat(i), i=1, 525)/70*0., 1., 6*0., 2*.177, .108, .225, .003, .06, .02, .025, .013, 2*.004, .007, .014, 2*.002, 2*.001, .054, .014, .016, .005, 2*.012, 5*.006, .002, 2*.001, 5*.002, 6*0., 1., 28*0., .143, .111, .143, .111, .143, .085, 2*0., .03, .058, .03, .058, .03, .058, 3*0., .25, .01, 2*0., .01, .25, 4*0., .24, 5*0., 3*.08, 3*0., .01, .08, .82, 5*0., .09, 6*0., .143, .111, .143, .111, .143, .085, 2*0., .03, .058, .03, .058, .03, .058, 4*0., 1., 5*0., 4*.215, 2*0., 2*.07, 0., 1., 2*.08, .76, .08, 2*.112, .05, .476, .08, .14, .01, .015, .005, 1., 0., 1., 0., 1., 0., .25, .01, 2*0., .01, .25, 4*0., .24, 5*0., 3*.08, 0., 1., 2*.5, .635, .212, .056, .017, .048, .032, .035, .03, 2*.015, .044, 2*.022, 9*.001, .035, .03, 2*.015, .044, 2*.022, 9*.001, .028, .017, .066, .02, .008, 2*.006, .003, .001, 2*.002, .003, .001, 2*.002, .005, .002, .005, .006, .004, .012, 2*.005, .008, 2*.005, .037, .004, .067, 2*.01, 2*.001, 3*.002, .003, 8*.002, .005, 4*.004, .015, .005, .027, 2*.005, .007, .014, .007, .01, .008, &
    .012, .015, 11*.002, 3*.004, .002, .004, 6*.002, 2*.004, .005, .011, .005, .015, .02, 2*.01, 3*.004, 5*.002, .015, .02, 2*.01, 3*.004, 5*.002, .038, .048, .082, .06, .028, .021, 2*.005, 2*.002, .005, .018, .005, .01, .008, .005, 3*.004, .001, 3*.003, .001, 2*.002, .003, 2*.002, 2*.001, .002, .001, .002, .001, .005, 4*.003, .001, 2*.002, .003, 2*.001, .013, .03, .058, .055, 3*.003, 2*.01, .007, .019, 4*.005, .015, 3*.005, 8*.002, 3*.001, .002, 2*.001, .003, 16*.001/
  Data (brat(i), i=526, 893)/.019, 2*.003, .002, .005, .004, .008, .003, .006, .003, .01, 5*.002, 2*.001, 2*.002, 11*.001, .002, 14*.001, .018, .005, .01, 2*.015, .017, 4*.015, .017, 3*.015, .025, .08, 2*.025, .04, .001, 2*.005, .02, .04, 2*.06, .04, .01, 4*.005, .25, .115, 3*1., .988, .012, .389, .319, .237, .049, .005, .001, .441, .205, .301, .03, .022, .001, 6*1., .665, .333, .002, .666, .333, .001, .49, .34, .17, .52, .48, 5*1., .893, .08, .017, 2*.005, .495, .343, 3*.043, .019, .013, .001, 2*.069, .862, 3*.027, .015, .045, .015, .045, .77, .029, 6*.02, 5*.05, .115, .015, .5, 0., 3*1., .28, .14, .313, .157, .11, .28, .14, .313, .157, .11, .667, .333, .667, .333, 1., .667, .333, .667, .333, 2*.5, 1., .333, .334, .333, 4*.25, 2*1., .3, .7, 2*1., .8, 2*.1, .667, .333, .667, .333, .6, .3, .067, .033, .6, .3, .067, .033, 2*.5, .6, .3, .067, .033, .6, .3, .067, .033, 2*.4, 2*.1, .8, 2*.1, .52, .26, 2*.11, .62, .31, 2*.035, .007, .993, .02, .98, .3, .7, 2*1., 2*.5, .667, .333, .667, .333, .667, .333, .667, .333, &
    2*.35, .3, .667, .333, .667, .333, 2*.35, .3, 2*.5, 3*.14, .1, .05, 4*.08, .028, .027, .028, .027, 4*.25, .273, .727, .35, .65, .3, .7, 2*1., 2*.35, .144, .105, .048, .003, .332, .166, .168, .084, .086, .043, .059, 2*.029, 2*.002, .332, .166, .168, .084, .086, .043, .059, 2*.029, 2*.002, .3, .15, .16, .08, .13, .06, .08, .04, .3, .15, .16, .08, .13, .06, .08, .04, 2*.4, .1, 2*.05, .3, .15, .16, .08, .13, .06, .08, .04, .3, .15, .16, .08, .13, .06, .08, .04, 2*.4, .1, 2*.05, 2*.35, .144, .105, 2*.024/
  Data (brat(i), i=894, 2000)/.003, .573, .287, .063, .028, 2*.021, .004, .003, 2*.5, .15, .85, .22, .78, .3, .7, 2*1., .217, .124, 2*.193, 2*.135, .002, .001, .686, .314, .641, .357, 2*.001, .018, 2*.005, .003, .002, 2*.006, .018, 2*.005, .003, .002, 2*.006, .005, .025, .015, .006, 2*.005, .004, .005, 5*.004, 2*.002, 2*.004, .003, .002, 2*.003, 3*.002, 2*.001, .002, 2*.001, 2*.002, 5*.001, 4*.003, 2*.005, 2*.002, 2*.001, 2*.002, 2*.001, .255, .057, 2*.035, .15, 2*.075, .03, 2*.015, 5*1., .999, .001, 1., .516, .483, .001, 1., .995, .005, 13*1., .331, .663, .006, .663, .331, .006, 1., .88, 2*.06, .88, 2*.06, .88, 2*.06, .667, 2*.333, .667, .676, .234, .085, .005, 3*1., 4*.5, 7*1., 935*0./
  Data (kfdp(i,1), i=1, 499)/21, 22, 23, 4* -24, 25, 21, 22, 23, 4*24, 25, 21, 22, 23, 4* -24, 25, 21, 22, 23, 4*24, 25, 21, 22, 23, 4* -24, 25, 21, 22, 23, 4*24, 25, 21, 22, 23, 4* -24, 25, 21, 22, 23, 4*24, 25, 22, 23, -24, 25, 23, 24, -12, 22, 23, -24, 25, 23, 24, -12, -14, 34*16, 22, 23, -24, 25, 23, 24, -89, 22, 23, -24, 25, 23, 24, 1, 2, 3, 4, 5, 6, 7, 8, 21, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 15, 17, 37, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 37, 4* -1, 4* -3, 4* -5, 4* -7, -11, -13, -15, -17, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 15, 17, 21, 2*22, 23, 24, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, -1, -3, -5, -7, -11, -13, -15, -17, 1, 2, 3, 4, 5, 6, 11, 13, 15, 82, -11, -13, 2*2, -12, -14, -16, 2* -2, 2* -4, -2, -4, 2*89, 2* -89, 2*89, 4* -1, 4* -3, 4* -5, 4* -7, -11, -13, -15, -17, -13, 130, 310, -13, 3*211, 12, 14, 16* -11, 16* -13, -311, -313, -311, -313, -311, -313, -311, -313, 2*111, 2*221, 2*331, 2*113, 2*223, 2*333, -311, -313, 2* -311, -313, 3* -311, -321, -323, -321, &
    2*211, 2*213, -213, 113, 3*213, 3*211, 2*213, 2* -311, -313, -321, 2* -311, -313, -311, -313, 4* -311, -321, -323, 2* -321, 3*211, 213, 2*211, 213, 5*211, 213, 4*211, 3*213, 211, 213, 321, 311, 3, 2*2, 12* -11, 12* -13, -321, -323, -321, -323, -311, -313, -311, -313, -311, -313, -311, -313, -311, -313, -311, -321, -323, -321, -323, 211, 213, 211, 213, 111, 221, 331, 113, 223, 333, 221, 331, 113, 223, 113, 223, 113, 223, 333, 223, 333, 321, 323, 321, 323, 311, 313, -321, -323, 3* -321, -323, 2* -321, -323, -321, -311, -313, 3* -311, -313, 2* -311, -313, -321, -323, 3* -321/
  Data (kfdp(i,1), i=500, 873)/ -323, 2* -321, -311, 2*333, 211, 213, 2*211, 2*213, 4*211, 10*111, -321, -323, 5* -321, -323, 2* -321, -311, -313, 4* -311, -313, 4* -311, -321, -323, 2* -321, -323, -321, -313, -311, -313, -311, 211, 213, 2*211, 213, 4*211, 111, 221, 113, 223, 113, 223, 2*3, -15, 5* -11, 5* -13, 221, 331, 333, 221, 331, 333, 211, 213, 211, 213, 321, 323, 321, 323, 2212, 221, 331, 333, 221, 2*2, 3*0, 3*22, 111, 211, 2*22, 2*211, 111, 3*22, 111, 3*21, 2*0, 211, 321, 3*311, 2*321, 421, 2*411, 2*421, 431, 511, 521, 531, 2*211, 22, 211, 2*111, 321, 130, -213, 113, 213, 211, 22, 111, 11, 13, 82, 11, 13, 15, 1, 2, 3, 4, 21, 22, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 21, 22, 2*89, 2*0, 223, 321, 311, 323, 313, 2*311, 321, 313, 323, 321, 421, 2*411, 421, 433, 521, 2*511, 521, 523, 513, 223, 213, 113, -213, 313, -313, 323, -323, 82, 21, 663, 21, 2*0, 221, 213, 113, 321, 2*311, 321, 421, 411, 423, 413, 411, 421, 413, 423, 431, 433, 521, 511, 523, 513, 511, 521, 513, 523, 521, 511, 531, 533, 221, 213, &
    -213, 211, 111, 321, 130, 211, 111, 321, 130, 443, 82, 553, 21, 663, 21, 2*0, 113, 213, 323, 2*313, 323, 423, 2*413, 423, 421, 411, 433, 523, 2*513, 523, 521, 511, 533, 213, -213, 10211, 10111, -10211, 2*221, 213, 2*113, -213, 2*321, 2*311, 313, -313, 323, -323, 443, 82, 553, 21, 663, 21, 2*0, 213, 113, 221, 223, 321, 211, 321, 311, 323, 313, 323, 313, 321, 5*311, 321, 313, 323, 313, 323, 311, 4*321, 421, 411, 423, 413, 423, 413, 421, 2*411, 421, 413, 423, 413, 423, 411, 2*421, 411, 433, 2*431, 521, 511, 523, 513, 523, 513, 521/
  Data (kfdp(i,1), i=874, 2000)/2*511, 521, 513, 523, 513, 523, 511, 2*521, 511, 533, 2*531, 213, -213, 221, 223, 321, 130, 111, 211, 111, 2*211, 321, 130, 221, 111, 321, 130, 443, 82, 553, 21, 663, 21, 2*0, 111, 211, -12, 12, -14, 14, 211, 111, 211, 111, 2212, 2*2112, -12, 7* -11, 7* -13, 2*2224, 2*2212, 2*2214, 2*3122, 2*3212, 2*3214, 5*3222, 4*3224, 2*3322, 3324, 2*2224, 5*2212, 5*2214, 2*2112, 2*2114, 2*3122, 2*3212, 2*3214, 2*3222, 2*3224, 4*2, 3, 2*2, 1, 2*2, 5*0, 2112, -12, 3122, 2212, 2112, 2212, 3*3122, 3*4122, 4132, 4232, 0, 3*5122, 5132, 5232, 0, 2112, 2212, 2*2112, 2212, 2112, 2*2212, 3122, 3212, 3112, 3122, 3222, 3112, 3122, 3222, 3212, 3322, 3312, 3322, 3312, 3122, 3322, 3312, -12, 3*4122, 2*4132, 2*4232, 4332, 3*5122, 5132, 5232, 5332, 935*0/
  Data (kfdp(i,2), i=1, 496)/3*1, 2, 4, 6, 8, 1, 3*2, 1, 3, 5, 7, 2, 3*3, 2, 4, 6, 8, 3, 3*4, 1, 3, 5, 7, 4, 3*5, 2, 4, 6, 8, 5, 3*6, 1, 3, 5, 7, 6, 3*7, 2, 4, 6, 8, 7, 3*8, 1, 3, 5, 7, 8, 2*11, 12, 11, 12, 2*11, 2*13, 14, 13, 14, 13, 11, 13, -211, -213, -211, -213, -211, -213, 3* -211, -321, -323, -321, -323, 2* -321, 4* -211, -213, -211, -213, -211, -213, -211, -213, -211, -213, 6* -211, 2*15, 16, 15, 16, 15, 18, 2*17, 18, 17, 18, 17, -1, -2, -3, -4, -5, -6, -7, -8, 21, -1, -2, -3, -4, -5, -6, -7, -8, -11, -13, -15, -17, -37, -1, -2, -3, -4, -5, -6, -7, -8, -11, -12, -13, -14, -15, -16, -17, -18, -37, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 12, 14, 16, 18, -1, -2, -3, -4, -5, -6, -7, -8, -11, -13, -15, -17, 21, 22, 2*23, -24, -1, -2, -3, -4, -5, -6, -7, -8, -11, -12, -13, -14, -15, -16, -17, -18, 2, 4, 6, 8, 12, 14, 16, 18, -3, -4, -5, -6, -7, -8, -13, -15, -17, -82, 12, 14, -1, -3, 11, 13, 15, 1, 4, 3, 4, 1, 3, 5, 3, 6, 4, 7, 5, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 12, 14, 16, 18, 14, &
    2*0, 14, 111, 211, 111, -11, -13, 16*12, 16*14, 2*211, 2*213, 2*321, 2*323, 211, 213, 211, 213, 211, 213, 211, 213, 211, 213, 211, 213, 2*211, 213, 7*211, 213, 211, 111, 211, 111, 2*211, -213, 213, 2*113, 223, 2*113, 221, 321, 2*311, 321, 313, 4*211, 213, 113, 213, -213, 2*211, 213, 113, 111, 221, 331, 111, 113, 223, 4*113, 223, 6*211, 213, 4*211, -321, -311, 3* -1, 12*12, 12*14, 2*211, 2*213, 2*111, 2*221, 2*331, 2*113, 2*223, 333, 2*321, 2*323, 2* -211, 2* -213, 6*111, 4*221, 2*331, 3*113, 2*223, 2* -211, 2* -213, 113, 111, 2*211, 213, 6*211, 321, 2*211, 213, 211, 2*111, 113, 2*223, 2*321/
  Data (kfdp(i,2), i=497, 863)/323, 321, 2*311, 313, 2*311, 111, 211, 2* -211, -213, -211, -213, -211, -213, 3* -211, 5*111, 2*113, 223, 113, 223, 2*211, 213, 5*211, 213, 3*211, 213, 2*211, 2*111, 221, 113, 223, 3*321, 323, 2*321, 323, 311, 313, 311, 313, 3*211, 2* -211, -213, 3* -211, 4*111, 2*113, 2* -1, 16, 5*12, 5*14, 3*211, 3*213, 2*111, 2*113, 2* -311, 2* -313, -2112, 3*321, 323, 2* -1, 3*0, 22, 11, 22, 111, -211, 211, 11, 2* -211, 111, 113, 223, 22, 111, 3*21, 2*0, 111, -211, 111, 22, 211, 111, 22, 211, 111, 22, 111, 5*22, 2* -211, 111, -211, 2*111, -321, 310, 211, 111, 2* -211, 221, 22, -11, -13, -82, -11, -13, -15, -1, -2, -3, -4, 2*21, -11, -12, -13, -14, -15, -16, -1, -2, -3, -4, -5, 2*21, 5, 3, 2*0, 211, -213, 113, -211, 111, 223, 211, 111, 211, 111, 223, 211, 111, -211, 2*111, -211, 111, 211, 111, -321, -311, 111, -211, 111, 211, -311, 311, -321, 321, -82, 21, 22, 21, 2*0, 211, 111, 211, -211, 111, 211, 111, 211, 111, 211, 111, -211, 111, -211, 3*111, -211, 111, -211, 111, 211, 111, 211, 111, &
    -321, -311, 3*111, -211, 211, -211, 111, -321, 310, -211, 111, -321, 310, 22, -82, 22, 21, 22, 21, 2*0, 211, 111, -211, 111, 211, 111, 211, 111, -211, 111, 321, 311, 111, -211, 111, 211, 111, -321, -311, 111, -211, 211, -211, 111, 2*211, 111, -211, 211, 111, 211, -321, 2* -311, -321, -311, 311, -321, 321, 22, -82, 22, 21, 22, 21, 2*0, 111, 3*211, -311, 22, -211, 111, -211, 111, -211, 211, -213, 113, 223, 221, 22, 211, 111, 211, 111, 2*211, 213, 113, 223, 221, 22, 211, 111, 211, 111, 4*211, -211, 111, -211, 111, -211, 211, -211, 211, 321, 311/
  Data (kfdp(i,2), i=864, 2000)/2*111, 211, -211, 111, -211, 111, -211, 211, -211, 2*211, 111, 211, 111, 4*211, -321, -311, 2*111, 211, -211, 211, 111, 211, -321, 310, 22, -211, 111, 2* -211, -321, 310, 221, 111, -321, 310, 22, -82, 22, 21, 22, 21, 2*0, 111, -211, 11, -11, 13, -13, -211, 111, -211, 111, -211, 111, 22, 11, 7*12, 7*14, -321, -323, -311, -313, -311, -313, 211, 213, 211, 213, 211, 213, 111, 221, 331, 113, 223, 111, 221, 113, 223, 321, 323, 321, -211, -213, 111, 221, 331, 113, 223, 111, 221, 331, 113, 223, 211, 213, 211, 213, 321, 323, 321, 323, 321, 323, 311, 313, 311, 313, 2* -1, -3, -1, 2203, 2*3201, 2203, 2101, 2103, 5*0, -211, 11, 22, 111, 211, 22, -211, 111, 22, -211, 111, 211, 2*22, 0, -211, 111, 211, 2*22, 0, 2* -211, 111, 22, 111, 211, 22, 211, 2* -211, 2*111, -211, 2*211, 111, 211, -211, 2*111, 211, -321, -211, 111, 11, -211, 111, 211, 111, 22, 111, 2*22, -211, 111, 211, 3*22, 935*0/
  Data (kfdp(i,3), i=1, 918)/70*0, 14, 6*0, 2*16, 2*0, 5*111, 310, 130, 2*0, 2*111, 310, 130, 113, 211, 223, 221, 2*113, 2*211, 2*223, 2*221, 2*113, 221, 113, 2*213, -213, 123*0, 4*3, 4*4, 1, 4, 3, 2*2, 6*81, 25*0, -211, 3*111, -311, -313, -311, 2* -321, 2* -311, 111, 221, 331, 113, 223, 211, 111, 211, 111, -311, -313, -311, 2* -321, 2* -311, 111, 221, 331, 113, 223, 211, 111, 211, 111, 20*0, 3*111, 2*221, 331, 113, 223, 3*211, -211, 111, -211, 111, 211, 111, 211, -211, 111, 113, 111, 223, 2*111, -311, 4*211, 2*111, 2*211, 111, 7*211, 7*111, 113, 221, 2*223, 2* -211, -213, 4* -211, -213, -211, -213, -211, 2*211, 2, 2*0, -321, -323, -311, -321, -311, 2* -321, -211, -213, 2* -211, 211, -321, -323, -311, -321, -311, 2* -321, -211, -213, 2* -211, 211, 46*0, 3*111, 113, 2*221, 331, 2*223, -311, 3* -211, -213, 8*111, 113, 3*211, 213, 2*111, -211, 3*111, 113, 111, 2*113, 221, 331, 223, 111, 221, 331, 113, 223, 113, 2*223, 2*221, 3*111, 221, 113, 223, 4*211, 3* -211, -213, -211, 5*111, -321, 3*211, 3*111, 2*211, &
    2*111, 2* -211, -213, 3*111, 221, 113, 223, 6*111, 3*0, 221, 331, 333, 321, 311, 221, 331, 333, 321, 311, 19*0, 3, 5*0, -11, 0, 2*111, -211, -11, 11, 2*221, 3*0, 111, 22*0, 111, 2*0, 22, 111, 5*0, 111, 12*0, 2*21, 11*0, 2*21, 2* -6, 111*0, -211, 2*111, -211, 3*111, -211, 111, 211, 15*0, 111, 6*0, 111, -211, 9*0, 111, -211, 9*0, 111, -211, 111, -211, 4*0, 111, -211, 111, -211, 4*0, -211, 4*0, 111, -211, 111, -211, 4*0, 111, -211, 111, -211, 4*0, -211, 3*0, -211, 5*0, 111, 211, 3*0, 111, 10*0, 2*111, 211, -211, 211, -211/
  Data (kfdp(i,3), i=919, 2000)/7*0, 2212, 3122, 3212, 3214, 2112, 2114, 2212, 2112, 3122, 3212, 3214, 2112, 2114, 2212, 2112, 50*0, 3*3, 1, 12*0, 2112, 43*0, 3322, 949*0/
  Data (kfdp(i,4), i=1, 2000)/83*0, 3*111, 9*0, -211, 3*0, 111, 2* -211, 0, 111, 0, 2*111, 113, 221, 111, -213, -211, 211, 123*0, 13*81, 37*0, 111, 3*211, 111, 5*0, -211, 111, -211, 111, 2*0, 111, 3*211, 111, 5*0, -211, 111, -211, 111, 50*0, 2*111, 2* -211, 2*111, -211, 211, 3*111, 211, 14*111, 221, 113, 223, 2*111, 2*113, 223, 2*111, -1, 4*0, -211, 111, -211, 211, 111, 2*0, 2*111, -211, 2*0, -211, 111, -211, 211, 111, 2*0, 2*111, -211, 96*0, 6*111, 3* -211, -213, 4*111, 113, 6*111, 3* -211, 3*111, 2* -211, 2*111, 3* -211, 12*111, 6*0, -321, -311, 3*0, -321, -311, 19*0, -3, 11*0, -11, 280*0, 111, -211, 3*0, 111, 29*0, -211, 111, 5*0, -211, 111, 50*0, 2101, 2103, 2*2101, 1006*0/
  Data (kfdp(i,5), i=1, 2000)/85*0, 111, 15*0, 111, 7*0, 111, 0, 2*111, 175*0, 111, -211, 111, 7*0, 2*111, 4*0, 111, -211, 111, 7*0, 2*111, 93*0, 111, -211, 111, 3*0, 111, -211, 4*0, 111, -211, 111, 3*0, 111, -211, 1571*0/

!...LUDAT4, with character strings.
  Data (chaf(i), i=1, 331)/'d', 'u', 's', 'c', 'b', 't', 'l', 'h', 2*' ', 'e', 'nu_e', 'mu', 'nu_mu', 'tau', 'nu_tau', 'chi', 'nu_chi', 2*' ', 'g', 'gamma', 'Z', 'W', 'H', 6*' ', 'Z''', 'Z"', 'W''', 'H''', 'H"', 'H', 2*' ', 'R', 40*' ', 'specflav', 'rndmflav', 'phasespa', 'c-hadron', 'b-hadron', 't-hadron', 'l-hadron', 'h-hadron', 'Wvirt', 'diquark', 'cluster', 'string', 'indep.', 'CMshower', 'SPHEaxis', 'THRUaxis', 'CLUSjet', 'CELLjet', 'table', ' ', 'pi', 2*'K', 2*'D', 'D_s', 2*'B', 'B_s', ' ', 'pi', 'eta', 'eta''', 'eta_c', 'eta_b', 'eta_t', 'eta_l', 'eta_h', 2*' ', 'rho', 2*'K*', 2*'D*', 'D*_s', 2*'B*', 'B*_s', ' ', 'rho', 'omega', 'phi', 'J/psi', 'Upsilon', 'Theta', 'Theta_l', 'Theta_h', 2*' ', 'b_1', 2*'K_1', 2*'D_1', 'D_1s', 2*'B_1', 'B_1s', ' ', 'b_1', 'h_1', 'h''_1', 'h_1c', 'h_1b', 'h_1t', 'h_1l', 'h_1h', 2*' ', 'a_0', 2*'K*_0', 2*'D*_0', 'D*_0s', 2*'B*_0', 'B*_0s', ' ', 'a_0', 'f_0', 'f''_0', 'chi_0c', 'chi_0b', 'chi_0t', 'chi_0l', 'chi_0h', 2*' ', 'a_1', 2*'K*_1', 2*'D*_1', 'D*_1s', 2*'B*_1', &
    'B*_1s', ' ', 'a_1', 'f_1', 'f''_1', 'chi_1c', 'chi_1b', 'chi_1t', 'chi_1l', 'chi_1h', 2*' ', 'a_2', 2*'K*_2', 2*'D*_2', 'D*_2s', 2*'B*_2', 'B*_2s', ' ', 'a_2', 'f_2', 'f''_2', 'chi_2c', 'chi_2b', 'chi_2t', 'chi_2l', 'chi_2h', 2*' ', 'K_L', 'K_S', 58*' ', 'pi_diffr', 'n_diffr', 'p_diffr', 22*' ', 'Lambda', 5*' ', 'Lambda_c', ' ', 2*'Xi_c', 6*' ', 'Lambda_b', ' ', 2*'Xi_b', 6*' '/
  Data (chaf(i), i=332, 500)/'n', 'p', ' ', 3*'Sigma', 2*'Xi', ' ', 3*'Sigma_c', 2*'Xi''_c', 'Omega_c', 4*' ', 3*'Sigma_b', 2*'Xi''_b', 'Omega_b', 4*' ', 4*'Delta', 3*'Sigma*', 2*'Xi*', 'Omega', 3*'Sigma*_c', 2*'Xi*_c', 'Omega*_c', 4*' ', 3*'Sigma*_b', 2*'Xi*_b', 'Omega*_b', 114*' '/

!...LUDATR, with initial values for the random number generator.
  Data mrlu/19780503, 0, 0, 97, 33, 0/

End Block Data ludata
Subroutine pyinit(frame, beam, target, win)

!...Initializes the generation procedure; finds maxima of the
!...differential cross-sections to be used for weighting.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /ludat4/chaf(500)
  Character chaf*8
  Save /ludat4/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Character *(*) frame, beam, target
  Character chfram*8, chbeam*8, chtarg*8, chmo(12)*3, chlh(2)*6
  Data chmo/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/, chlh/'lepton', 'hadron'/

!lin-12/2012 correct NN differential cross section in HIJING:
  Write (mstu(11), *) 'In PYINIT: BEAM,TARGET= ', beam, target

!...Write headers.
!      IF(MSTP(122).GE.1) WRITE(MSTU(11),1000) MSTP(181),MSTP(182),
!     &MSTP(185),CHMO(MSTP(184)),MSTP(183)
  Call lulist(0)
!      IF(MSTP(122).GE.1) WRITE(MSTU(11),1100)

!...Identify beam and target particles and initialize kinematics.
  chfram = frame // ' '
  chbeam = beam // ' '
  chtarg = target // ' '
  Call pyinki(chfram, chbeam, chtarg, win)

!...Select partonic subprocesses to be included in the simulation.
  If (msel/=0) Then
    Do i = 1, 200
      msub(i) = 0
    End Do
  End If
  If (mint(43)==1 .And. (msel==1 .Or. msel==2)) Then
!...Lepton+lepton -> gamma/Z0 or W.
    If (mint(11)+mint(12)==0) msub(1) = 1
    If (mint(11)+mint(12)/=0) msub(2) = 1
  Else If (msel==1) Then
!...High-pT QCD processes:
    msub(11) = 1
    msub(12) = 1
    msub(13) = 1
    msub(28) = 1
    msub(53) = 1
    msub(68) = 1
    If (mstp(82)<=1 .And. ckin(3)<parp(81)) msub(95) = 1
    If (mstp(82)>=2 .And. ckin(3)<parp(82)) msub(95) = 1
  Else If (msel==2) Then
!...All QCD processes:
    msub(11) = 1
    msub(12) = 1
    msub(13) = 1
    msub(28) = 1
    msub(53) = 1
    msub(68) = 1
    msub(91) = 1
    msub(92) = 1
    msub(93) = 1
    msub(95) = 1
  Else If (msel>=4 .And. msel<=8) Then
!...Heavy quark production.
    msub(81) = 1
    msub(82) = 1
    Do j = 1, min(8, mdcy(21,3))
      mdme(mdcy(21,2)+j-1, 1) = 0
    End Do
    mdme(mdcy(21,2)+msel-1, 1) = 1
  Else If (msel==10) Then
!...Prompt photon production:
    msub(14) = 1
    msub(18) = 1
    msub(29) = 1
  Else If (msel==11) Then
!...Z0/gamma* production:
    msub(1) = 1
  Else If (msel==12) Then
!...W+/- production:
    msub(2) = 1
  Else If (msel==13) Then
!...Z0 + jet:
    msub(15) = 1
    msub(30) = 1
  Else If (msel==14) Then
!...W+/- + jet:
    msub(16) = 1
    msub(31) = 1
  Else If (msel==15) Then
!...Z0 & W+/- pair production:
    msub(19) = 1
    msub(20) = 1
    msub(22) = 1
    msub(23) = 1
    msub(25) = 1
  Else If (msel==16) Then
!...H0 production:
    msub(3) = 1
    msub(5) = 1
    msub(8) = 1
    msub(102) = 1
  Else If (msel==17) Then
!...H0 & Z0 or W+/- pair production:
    msub(24) = 1
    msub(26) = 1
  Else If (msel==21) Then
!...Z'0 production:
    msub(141) = 1
  Else If (msel==22) Then
!...H+/- production:
    msub(142) = 1
  Else If (msel==23) Then
!...R production:
    msub(143) = 1
  End If

!...Count number of subprocesses on.
  mint(44) = 0
  Do isub = 1, 200
    If (mint(43)<4 .And. isub>=91 .And. isub<=96 .And. msub(isub)==1) Then
      Write (mstu(11), 1200) isub, chlh(mint(41)), chlh(mint(42))
      Stop
    Else If (msub(isub)==1 .And. iset(isub)==-1) Then
      Write (mstu(11), 1300) isub
      Stop
    Else If (msub(isub)==1 .And. iset(isub)<=-2) Then
      Write (mstu(11), 1400) isub
      Stop
    Else If (msub(isub)==1) Then
      mint(44) = mint(44) + 1
    End If
  End Do
  If (mint(44)==0) Then
    Write (mstu(11), 1500)
    Stop
  End If
  mint(45) = mint(44) - msub(91) - msub(92) - msub(93) - msub(94)

!...Maximum 4 generations; set maximum number of allowed flavours.
  mstp(1) = min(4, mstp(1))
  mstu(114) = min(mstu(114), 2*mstp(1))
  mstp(54) = min(mstp(54), 2*mstp(1))

!...Sum up Cabibbo-Kobayashi-Maskawa factors for each quark/lepton.
  Do i = -20, 20
    vint(180+i) = 0.
    ia = iabs(i)
    If (ia>=1 .And. ia<=2*mstp(1)) Then
      Do j = 1, mstp(1)
        ib = 2*j - 1 + mod(ia, 2)
        ipm = (5-isign(1,i))/2
        idc = j + mdcy(ia, 2) + 2
        If (mdme(idc,1)==1 .Or. mdme(idc,1)==ipm) vint(180+i) = vint(180+i) + vckm((ia+1)/2, (ib+1)/2)
      End Do
    Else If (ia>=11 .And. ia<=10+2*mstp(1)) Then
      vint(180+i) = 1.
    End If
  End Do

!...Choose Lambda value to use in alpha-strong.
  mstu(111) = mstp(2)
  If (mstp(3)>=1) Then
    alam = parp(1)
    If (mstp(51)==1) alam = 0.2
    If (mstp(51)==2) alam = 0.29
    If (mstp(51)==3) alam = 0.2
    If (mstp(51)==4) alam = 0.4
    If (mstp(51)==11) alam = 0.16
    If (mstp(51)==12) alam = 0.26
    If (mstp(51)==13) alam = 0.36
    parp(1) = alam
    parp(61) = alam
    paru(112) = alam
    parj(81) = alam
  End If

!...Initialize widths and partial widths for resonances.
  Call pyinre

!...Reset variables for cross-section calculation.
  Do i = 0, 200
    Do j = 1, 3
      ngen(i, j) = 0
      xsec(i, j) = 0.
    End Do
  End Do
  vint(108) = 0.

!...Find parametrized total cross-sections.
  If (mint(43)==4) Call pyxtot

!...Maxima of differential cross-sections.
  If (mstp(121)<=0) Call pymaxi

!...Initialize possibility of overlayed events.
  If (mstp(131)/=0) Call pyovly(1)

!...Initialize multiple interactions with variable impact parameter.
  If (mint(43)==4 .And. (mint(45)/=0 .Or. mstp(131)/=0) .And. mstp(82)>=2) Call pymult(1)
!lin 1600 FORMAT(/1X,22('*'),1X,'PYINIT: initialization completed',1X,
!lin     &22('*'))

  Return
!      IF(MSTP(122).GE.1) WRITE(MSTU(11),1600)

!...Formats for initialization information.
!lin 1000 FORMAT(///20X,'The Lund Monte Carlo - PYTHIA version ',I1,'.',I1/
!lin     &20X,'**  Last date of change:  ',I2,1X,A3,1X,I4,'  **'/)
!lin 1100 FORMAT('1',18('*'),1X,'PYINIT: initialization of PYTHIA ',
!lin     &'routines',1X,17('*'))
  1200 Format (1X, 'Error: process number ', I3, ' not meaningful for ', A6, '-', A6, ' interactions.'/1X, 'Execution stopped!')
  1300 Format (1X, 'Error: requested subprocess', I4, ' not implemented.'/1X, 'Execution stopped!')
  1400 Format (1X, 'Error: requested subprocess', I4, ' not existing.'/1X, 'Execution stopped!')
  1500 Format (1X, 'Error: no subprocess switched on.'/1X, 'Execution stopped.')
End Subroutine pyinit

!*********************************************************************

!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine pythia

!...Administers the generation of a high-pt event via calls to a number
!...of subroutines; also computes cross-sections.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/

!...Loop over desired number of overlayed events (normally 1).
  mint(7) = 0
  mint(8) = 0
  novl = 1
  If (mstp(131)/=0) Call pyovly(2)
  If (mstp(131)/=0) novl = mint(81)
  mint(83) = 0
  mint(84) = mstp(126)
  mstu(70) = 0
  Do iovl = 1, novl
    If (mint(84)+100>=mstu(4)) Then
      Call luerrm(11, '(PYTHIA:) no more space in LUJETS for overlayed events')
      If (mstu(21)>=1) Goto 200
    End If
    mint(82) = iovl

!...Generate variables of hard scattering.
    100 Continue
    If (iovl==1) ngen(0, 2) = ngen(0, 2) + 1
    mint(31) = 0
    mint(51) = 0
    Call pyrand
    isub = mint(1)
    If (iovl==1) Then
      ngen(isub, 2) = ngen(isub, 2) + 1

!...Store information on hard interaction.
      Do j = 1, 200
        msti(j) = 0
        pari(j) = 0.
      End Do
      msti(1) = mint(1)
      msti(2) = mint(2)
      msti(11) = mint(11)
      msti(12) = mint(12)
      msti(15) = mint(15)
      msti(16) = mint(16)
      msti(17) = mint(17)
      msti(18) = mint(18)
      pari(11) = vint(1)
      pari(12) = vint(2)
      If (isub/=95) Then
        Do j = 13, 22
          pari(j) = vint(30+j)
        End Do
        pari(33) = vint(41)
        pari(34) = vint(42)
        pari(35) = pari(33) - pari(34)
        pari(36) = vint(21)
        pari(37) = vint(22)
        pari(38) = vint(26)
        pari(41) = vint(23)
      End If
    End If

    If (mstp(111)==-1) Goto 160
    If (isub<=90 .Or. isub>=95) Then
!...Hard scattering (including low-pT):
!...reconstruct kinematics and colour flow of hard scattering.
      Call pyscat
      If (mint(51)==1) Goto 100

!...Showering of initial state partons (optional).
      ipu1 = mint(84) + 1
      ipu2 = mint(84) + 2
      If (mstp(61)>=1 .And. mint(43)/=1 .And. isub/=95) Call pysspa(ipu1, ipu2)
      nsav1 = n

!...Multiple interactions.
      If (mstp(81)>=1 .And. mint(43)==4 .And. isub/=95) Call pymult(6)
      mint(1) = isub
      nsav2 = n

!...Hadron remnants and primordial kT.
      Call pyremn(ipu1, ipu2)
      If (mint(51)==1) Goto 100
      nsav3 = n

!...Showering of final state partons (optional).
      ipu3 = mint(84) + 3
      ipu4 = mint(84) + 4
      If (mstp(71)>=1 .And. isub/=95 .And. k(ipu3,1)>0 .And. k(ipu3,1)<=10 .And. k(ipu4,1)>0 .And. k(ipu4,1)<=10) Then
        qmax = sqrt(parp(71)*vint(52))
        If (isub==5) qmax = sqrt(pmas(23,1)**2)
        If (isub==8) qmax = sqrt(pmas(24,1)**2)
        Call lushow(ipu3, ipu4, qmax)
      End If

!...Sum up transverse and longitudinal momenta.
      If (iovl==1) Then
        pari(65) = 2.*pari(17)
        Do i = mstp(126) + 1, n
          If (k(i,1)<=0 .Or. k(i,1)>10) Goto 130
          pt = sqrt(p(i,1)**2+p(i,2)**2)
          pari(69) = pari(69) + pt
          If (i<=nsav1 .Or. i>nsav3) pari(66) = pari(66) + pt
          If (i>nsav1 .And. i<=nsav2) pari(68) = pari(68) + pt
        130 End Do
        pari(67) = pari(68)
        pari(71) = vint(151)
        pari(72) = vint(152)
        pari(73) = vint(151)
        pari(74) = vint(152)
      End If

!...Decay of final state resonances.
      If (mstp(41)>=1 .And. isub/=95) Call pyresd

    Else
!...Diffractive and elastic scattering.
      Call pydiff
      If (iovl==1) Then
        pari(65) = 2.*pari(17)
        pari(66) = pari(65)
        pari(69) = pari(65)
      End If
    End If

!...Recalculate energies from momenta and masses (if desired).
    If (mstp(113)>=1) Then
      Do i = mint(83) + 1, n
        If (k(i,1)>0 .And. k(i,1)<=10) p(i, 4) = sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2+p(i,5)**2)
      End Do
    End If

!...Rearrange partons along strings, check invariant mass cuts.
    mstu(28) = 0
    Call luprep(mint(84)+1)
    If (mstp(112)==1 .And. mstu(28)==3) Goto 100
    If (mstp(125)==0 .Or. mstp(125)==1) Then
      Do i = mint(84) + 1, n
        If (k(i,2)/=94) Goto 150
        k(i+1, 3) = mod(k(i+1,4)/mstu(5), mstu(5))
        k(i+2, 3) = mod(k(i+2,4)/mstu(5), mstu(5))
      150 End Do
      Call luedit(12)
      Call luedit(14)
      If (mstp(125)==0) Call luedit(15)
      If (mstp(125)==0) mint(4) = 0
    End If

!...Introduce separators between sections in LULIST event listing.
    If (iovl==1 .And. mstp(125)<=0) Then
      mstu(70) = 1
      mstu(71) = n
    Else If (iovl==1) Then
      mstu(70) = 3
      mstu(71) = 2
      mstu(72) = mint(4)
      mstu(73) = n
    End If

!...Perform hadronization (if desired).
    If (mstp(111)>=1) Call luexec
    If (mstp(125)==0 .Or. mstp(125)==1) Call luedit(14)

!...Calculate Monte Carlo estimates of cross-sections.
    160 If (iovl==1) Then
      If (mstp(111)/=-1) ngen(isub, 3) = ngen(isub, 3) + 1
      ngen(0, 3) = ngen(0, 3) + 1
      xsec(0, 3) = 0.
      Do i = 1, 200
        If (i==96) Then
          xsec(i, 3) = 0.
        Else If (msub(95)==1 .And. (i==11 .Or. i==12 .Or. i==13 .Or. i==28 .Or. i==53 .Or. i==68)) Then
          xsec(i, 3) = xsec(96, 2)*ngen(i, 3)/max(1., float(ngen(96,1))*float(ngen(96,2)))
        Else If (ngen(i,1)==0) Then
          xsec(i, 3) = 0.
        Else If (ngen(i,2)==0) Then
          xsec(i, 3) = xsec(i, 2)*ngen(0, 3)/(float(ngen(i,1))*float(ngen(0,2)))
        Else
          xsec(i, 3) = xsec(i, 2)*ngen(i, 3)/(float(ngen(i,1))*float(ngen(i,2)))
        End If
        xsec(0, 3) = xsec(0, 3) + xsec(i, 3)
      End Do
      If (msub(95)==1) Then
        ngens = ngen(91, 3) + ngen(92, 3) + ngen(93, 3) + ngen(94, 3) + ngen(95, 3)
        xsecs = xsec(91, 3) + xsec(92, 3) + xsec(93, 3) + xsec(94, 3) + xsec(95, 3)
        xmaxs = xsec(95, 1)
        If (msub(91)==1) xmaxs = xmaxs + xsec(91, 1)
        If (msub(92)==1) xmaxs = xmaxs + xsec(92, 1)
        If (msub(93)==1) xmaxs = xmaxs + xsec(93, 1)
        If (msub(94)==1) xmaxs = xmaxs + xsec(94, 1)
        fac = 1.
        If (ngens<ngen(0,3)) fac = (xmaxs-xsecs)/(xsec(0,3)-xsecs)
        xsec(11, 3) = fac*xsec(11, 3)
        xsec(12, 3) = fac*xsec(12, 3)
        xsec(13, 3) = fac*xsec(13, 3)
        xsec(28, 3) = fac*xsec(28, 3)
        xsec(53, 3) = fac*xsec(53, 3)
        xsec(68, 3) = fac*xsec(68, 3)
        xsec(0, 3) = xsec(91, 3) + xsec(92, 3) + xsec(93, 3) + xsec(94, 3) + xsec(95, 1)
      End If

!...Store final information.
      mint(5) = mint(5) + 1
      msti(3) = mint(3)
      msti(4) = mint(4)
      msti(5) = mint(5)
      msti(6) = mint(6)
      msti(7) = mint(7)
      msti(8) = mint(8)
      msti(13) = mint(13)
      msti(14) = mint(14)
      msti(21) = mint(21)
      msti(22) = mint(22)
      msti(23) = mint(23)
      msti(24) = mint(24)
      msti(25) = mint(25)
      msti(26) = mint(26)
      msti(31) = mint(31)
      pari(1) = xsec(0, 3)
      pari(2) = xsec(0, 3)/mint(5)
      pari(31) = vint(141)
      pari(32) = vint(142)
      If (isub/=95 .And. mint(7)*mint(8)/=0) Then
        pari(42) = 2.*vint(47)/vint(1)
        Do is = 7, 8
          pari(36+is) = p(mint(is), 3)/vint(1)
          pari(38+is) = p(mint(is), 4)/vint(1)
          i = mint(is)
          pr = max(1E-20, p(i,5)**2+p(i,1)**2+p(i,2)**2)
          pari(40+is) = sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),1E20)), p(i,3))
          pr = max(1E-20, p(i,1)**2+p(i,2)**2)
          pari(42+is) = sign(log(min((sqrt(pr+p(i,3)**2)+abs(p(i,3)))/sqrt(pr),1E20)), p(i,3))
          pari(44+is) = p(i, 3)/sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
          pari(46+is) = ulangl(p(i,3), sqrt(p(i,1)**2+p(i,2)**2))
          pari(48+is) = ulangl(p(i,1), p(i,2))
        End Do
      End If
      pari(61) = vint(148)
      If (iset(isub)==1 .Or. iset(isub)==3) Then
        mstu(161) = mint(21)
        mstu(162) = 0
      Else
        mstu(161) = mint(21)
        mstu(162) = mint(22)
      End If
    End If

!...Prepare to go to next overlayed event.
    msti(41) = iovl
    If (iovl>=2 .And. iovl<=10) msti(40+iovl) = isub
    If (mstu(70)<10) Then
      mstu(70) = mstu(70) + 1
      mstu(70+mstu(70)) = n
    End If
    mint(83) = n
    mint(84) = n + mstp(126)
  End Do

!...Information on overlayed events.
  If (mstp(131)==1 .And. mstp(133)>=1) Then
    pari(91) = vint(132)
    pari(92) = vint(133)
    pari(93) = vint(134)
    If (mstp(133)==2) pari(93) = pari(93)*xsec(0, 3)/vint(131)
  End If

!...Transform to the desired coordinate frame.
  200 Call pyfram(mstp(124))

  Return
End Subroutine pythia

!*********************************************************************

Subroutine pyinki(chfram, chbeam, chtarg, win)

!...Identifies the two incoming particles and sets up kinematics,
!...including rotations and boosts to/from CM frame.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Character chfram*8, chbeam*8, chtarg*8, chcom(3)*8, chalp(2)*26, chidnt(3)*8, chtemp*8, chcde(18)*8, chinit*76
  Dimension len(3), kcde(18)
  Data chalp/'abcdefghijklmnopqrstuvwxyz', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
  Data chcde/'e-      ', 'e+      ', 'nue     ', 'nue~    ', 'mu-     ', 'mu+     ', 'numu    ', 'numu~   ', 'tau-    ', 'tau+    ', 'nutau   ', 'nutau~  ', 'pi+     ', 'pi-     ', 'n       ', 'n~      ', 'p       ', 'p~      '/
  Data kcde/11, -11, 12, -12, 13, -13, 14, -14, 15, -15, 16, -16, 211, -211, 2112, -2112, 2212, -2212/

!...Convert character variables to lowercase and find their length.
  chcom(1) = chfram
  chcom(2) = chbeam
  chcom(3) = chtarg
  Do i = 1, 3
    len(i) = 8
    Do ll = 8, 1, -1
      If (len(i)==ll .And. chcom(i)(ll:ll)==' ') len(i) = ll - 1
      Do la = 1, 26
        If (chcom(i)(ll:ll)==chalp(2)(la:la)) chcom(i)(ll:ll) = chalp(1)(la:la)
      End Do
    End Do
    chidnt(i) = chcom(i)
    Do ll = 1, 6
      If (chidnt(i)(ll:ll+2)=='bar') Then
        chtemp = chidnt(i)
        chidnt(i) = chtemp(1:ll-1) // '~' // chtemp(ll+3:8) // '  '
      End If
    End Do
    Do ll = 1, 8
      If (chidnt(i)(ll:ll)=='_') Then
        chtemp = chidnt(i)
        chidnt(i) = chtemp(1:ll-1) // chtemp(ll+1:8) // ' '
      End If
    End Do
  End Do

!...Set initial state. Error for unknown codes. Reset variables.
  n = 2
  Do i = 1, 2
    k(i, 2) = 0
    Do j = 1, 18
      If (chidnt(i+1)==chcde(j)) k(i, 2) = kcde(j)
    End Do
    p(i, 5) = ulmass(k(i,2))
    mint(40+i) = 1
    If (iabs(k(i,2))>100) mint(40+i) = 2
    Do j = 1, 5
      v(i, j) = 0.
    End Do
  End Do
  If (k(1,2)==0) Write (mstu(11), 1000) chbeam(1:len(2))
  If (k(2,2)==0) Write (mstu(11), 1100) chtarg(1:len(3))
  If (k(1,2)==0 .Or. k(2,2)==0) Stop
  Do j = 6, 10
    vint(j) = 0.
  End Do
  chinit = ' '

!...Set up kinematics for events defined in CM frame.
  If (chcom(1)(1:2)=='cm') Then
    If (chcom(2)(1:1)/='e') Then
      loffs = (34-(len(2)+len(3)))/2
      chinit(loffs+1:76) = 'PYTHIA will be initialized for a ' // chcom(2)(1:len(2)) // '-' // chcom(3)(1:len(3)) // ' collider' // ' '
    Else
      loffs = (33-(len(2)+len(3)))/2
      chinit(loffs+1:76) = 'PYTHIA will be initialized for an ' // chcom(2)(1:len(2)) // '-' // chcom(3)(1:len(3)) // ' collider' // ' '
    End If
!        WRITE(MSTU(11),1200) CHINIT
!        WRITE(MSTU(11),1300) WIN
    s = win**2
    p(1, 1) = 0.
    p(1, 2) = 0.
    p(2, 1) = 0.
    p(2, 2) = 0.
    p(1, 3) = sqrt(((s-p(1,5)**2-p(2,5)**2)**2-(2.*p(1,5)*p(2,5))**2)/(4.*s))
    p(2, 3) = -p(1, 3)
    p(1, 4) = sqrt(p(1,3)**2+p(1,5)**2)
    p(2, 4) = sqrt(p(2,3)**2+p(2,5)**2)

!...Set up kinematics for fixed target events.
  Else If (chcom(1)(1:3)=='fix') Then
    loffs = (29-(len(2)+len(3)))/2
    chinit(loffs+1:76) = 'PYTHIA will be initialized for ' // chcom(2)(1:len(2)) // ' on ' // chcom(3)(1:len(3)) // ' fixed target' // ' '
!        WRITE(MSTU(11),1200) CHINIT
!        WRITE(MSTU(11),1400) WIN
    p(1, 1) = 0.
    p(1, 2) = 0.
    p(2, 1) = 0.
    p(2, 2) = 0.
    p(1, 3) = win
    p(1, 4) = sqrt(p(1,3)**2+p(1,5)**2)
    p(2, 3) = 0.
    p(2, 4) = p(2, 5)
    s = p(1, 5)**2 + p(2, 5)**2 + 2.*p(2, 4)*p(1, 4)
    vint(10) = p(1, 3)/(p(1,4)+p(2,4))
    Call lurobo(0., 0., 0., 0., -vint(10))
!        WRITE(MSTU(11),1500) SQRT(S)

!...Set up kinematics for events in user-defined frame.
  Else If (chcom(1)(1:3)=='use') Then
    loffs = (13-(len(1)+len(2)))/2
    chinit(loffs+1:76) = 'PYTHIA will be initialized for ' // chcom(2)(1:len(2)) // ' on ' // chcom(3)(1:len(3)) // 'user-specified configuration' // ' '
!        WRITE(MSTU(11),1200) CHINIT
!        WRITE(MSTU(11),1600)
!        WRITE(MSTU(11),1700) CHCOM(2),P(1,1),P(1,2),P(1,3)
!        WRITE(MSTU(11),1700) CHCOM(3),P(2,1),P(2,2),P(2,3)
    p(1, 4) = sqrt(p(1,1)**2+p(1,2)**2+p(1,3)**2+p(1,5)**2)
    p(2, 4) = sqrt(p(2,1)**2+p(2,2)**2+p(2,3)**2+p(2,5)**2)
    Do j = 1, 3
      vint(7+j) = sngl((dble(p(1,j))+dble(p(2,j)))/dble(p(1,4)+p(2,4)))
    End Do
    Call lurobo(0., 0., -vint(8), -vint(9), -vint(10))
    vint(7) = ulangl(p(1,1), p(1,2))
    Call lurobo(0., -vint(7), 0., 0., 0.)
    vint(6) = ulangl(p(1,3), p(1,1))
    Call lurobo(-vint(6), 0., 0., 0., 0.)
    s = p(1, 5)**2 + p(2, 5)**2 + 2.*(p(1,4)*p(2,4)-p(1,3)*p(2,3))
!        WRITE(MSTU(11),1500) SQRT(S)

!...Unknown frame. Error for too low CM energy.
  Else
    Write (mstu(11), 1800) chfram(1:len(1))
    Stop
  End If
  If (s<parp(2)**2) Then
    Write (mstu(11), 1900) sqrt(s)
    Stop
  End If

!...Save information on incoming particles.
  mint(11) = k(1, 2)
  mint(12) = k(2, 2)
  mint(43) = 2*mint(41) + mint(42) - 2
  vint(1) = sqrt(s)
  vint(2) = s
  vint(3) = p(1, 5)
  vint(4) = p(2, 5)
  vint(5) = p(1, 3)

!...Store constants to be used in generation.
  If (mstp(82)<=1) vint(149) = 4.*parp(81)**2/s
  If (mstp(82)>=2) vint(149) = 4.*parp(82)**2/s

  Return

!...Formats for initialization and error information.
  1000 Format (1X, 'Error: unrecognized beam particle ''', A, '''.'/1X, 'Execution stopped!')
  1100 Format (1X, 'Error: unrecognized target particle ''', A, '''.'/1X, 'Execution stopped!')
!lin 1200 FORMAT(/1X,78('=')/1X,'I',76X,'I'/1X,'I',A76,'I')
! 1300 FORMAT(1X,'I',18X,'at',1X,F10.3,1X,'GeV center-of-mass energy',
!     &19X,'I'/1X,'I',76X,'I'/1X,78('='))
! 1400 FORMAT(1X,'I',22X,'at',1X,F10.3,1X,'GeV/c lab-momentum',22X,'I')
! 1500 FORMAT(1X,'I',76X,'I'/1X,'I',11X,'corresponding to',1X,F10.3,1X,
!     &'GeV center-of-mass energy',12X,'I'/1X,'I',76X,'I'/1X,78('='))
! 1600 FORMAT(1X,'I',76X,'I'/1X,'I',24X,'px (GeV/c)',3X,'py (GeV/c)',3X,
!     &'pz (GeV/c)',16X,'I')
!lin 1700 FORMAT(1X,'I',15X,A8,3(2X,F10.3,1X),15X,'I')
  1800 Format (1X, 'Error: unrecognized coordinate frame ''', A, '''.'/1X, 'Execution stopped!')
  1900 Format (1X, 'Error: too low CM energy,', F8.3, ' GeV for event ', 'generation.'/1X, 'Execution stopped!')
End Subroutine pyinki

!*********************************************************************

Subroutine pyinre

!...Calculates full and effective widths of guage bosons, stores masses
!...and widths, rescales coefficients to be used for resonance
!...production generation.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint6/proc(0:200)
  Character proc*28
  Save /pyint6/
  Dimension wdtp(0:40), wdte(0:40, 0:5)

!...Calculate full and effective widths of gauge bosons.
  aem = paru(101)
  xw = paru(102)
  Do i = 21, 40
    Do j = 0, 40
      widp(i, j) = 0.
      wide(i, j) = 0.
    End Do
  End Do

!...W+/-:
  wmas = pmas(24, 1)
  wfac = aem/(24.*xw)*wmas
  Call pywidt(24, wmas, wdtp, wdte)
  wids(24, 1) = ((wdte(0,1)+wdte(0,2))*(wdte(0,1)+wdte(0,3))+(wdte(0,1)+wdte(0,2)+wdte(0,1)+wdte(0,3))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(24, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(24, 3) = (wdte(0,1)+wdte(0,3)+wdte(0,4))/wdtp(0)
  Do i = 0, 40
    widp(24, i) = wfac*wdtp(i)
    wide(24, i) = wfac*wdte(i, 0)
  End Do

!...H+/-:
  hcmas = pmas(37, 1)
  hcfac = aem/(8.*xw)*(hcmas/wmas)**2*hcmas
  Call pywidt(37, hcmas, wdtp, wdte)
  wids(37, 1) = ((wdte(0,1)+wdte(0,2))*(wdte(0,1)+wdte(0,3))+(wdte(0,1)+wdte(0,2)+wdte(0,1)+wdte(0,3))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(37, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(37, 3) = (wdte(0,1)+wdte(0,3)+wdte(0,4))/wdtp(0)
  Do i = 0, 40
    widp(37, i) = hcfac*wdtp(i)
    wide(37, i) = hcfac*wdte(i, 0)
  End Do

!...Z0:
  zmas = pmas(23, 1)
  zfac = aem/(48.*xw*(1.-xw))*zmas
  Call pywidt(23, zmas, wdtp, wdte)
  wids(23, 1) = ((wdte(0,1)+wdte(0,2))**2+2.*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(23, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(23, 3) = 0.
  Do i = 0, 40
    widp(23, i) = zfac*wdtp(i)
    wide(23, i) = zfac*wdte(i, 0)
  End Do

!...H0:
  hmas = pmas(25, 1)
  hfac = aem/(8.*xw)*(hmas/wmas)**2*hmas
  Call pywidt(25, hmas, wdtp, wdte)
  wids(25, 1) = ((wdte(0,1)+wdte(0,2))**2+2.*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(25, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(25, 3) = 0.
  Do i = 0, 40
    widp(25, i) = hfac*wdtp(i)
    wide(25, i) = hfac*wdte(i, 0)
  End Do

!...Z'0:
  zpmas = pmas(32, 1)
  zpfac = aem/(48.*xw*(1.-xw))*zpmas
  Call pywidt(32, zpmas, wdtp, wdte)
  wids(32, 1) = ((wdte(0,1)+wdte(0,2)+wdte(0,3))**2+2.*(wdte(0,1)+wdte(0,2))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(32, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(32, 3) = 0.
  Do i = 0, 40
    widp(32, i) = zpfac*wdtp(i)
    wide(32, i) = zpfac*wdte(i, 0)
  End Do

!...R:
  rmas = pmas(40, 1)
  rfac = 0.08*rmas/((mstp(1)-1)*(1.+6.*(1.+ulalps(rmas**2)/paru(1))))
  Call pywidt(40, rmas, wdtp, wdte)
  wids(40, 1) = ((wdte(0,1)+wdte(0,2))*(wdte(0,1)+wdte(0,3))+(wdte(0,1)+wdte(0,2)+wdte(0,1)+wdte(0,3))*(wdte(0,4)+wdte(0,5))+2.*wdte(0,4)*wdte(0,5))/wdtp(0)**2
  wids(40, 2) = (wdte(0,1)+wdte(0,2)+wdte(0,4))/wdtp(0)
  wids(40, 3) = (wdte(0,1)+wdte(0,3)+wdte(0,4))/wdtp(0)
  Do i = 0, 40
    widp(40, i) = wfac*wdtp(i)
    wide(40, i) = wfac*wdte(i, 0)
  End Do

!...Q:
  kflqm = 1
  Do i = 1, min(8, mdcy(21,3))
    idc = i + mdcy(21, 2) - 1
    If (mdme(idc,1)<=0) Goto 170
    kflqm = i
  170 End Do
  mint(46) = kflqm
  kfpr(81, 1) = kflqm
  kfpr(81, 2) = kflqm
  kfpr(82, 1) = kflqm
  kfpr(82, 2) = kflqm

!...Set resonance widths and branching ratios in JETSET.
  Do i = 1, 6
    If (i<=3) kc = i + 22
    If (i==4) kc = 32
    If (i==5) kc = 37
    If (i==6) kc = 40
    pmas(kc, 2) = widp(kc, 0)
    pmas(kc, 3) = min(0.9*pmas(kc,1), 10.*pmas(kc,2))
    Do j = 1, mdcy(kc, 3)
      idc = j + mdcy(kc, 2) - 1
      brat(idc) = wide(kc, j)/wide(kc, 0)
    End Do
  End Do

!...Special cases in treatment of gamma*/Z0: redefine process name.
  If (mstp(43)==1) Then
    proc(1) = 'f + fb -> gamma*'
  Else If (mstp(43)==2) Then
    proc(1) = 'f + fb -> Z0'
  Else If (mstp(43)==3) Then
    proc(1) = 'f + fb -> gamma*/Z0'
  End If

!...Special cases in treatment of gamma*/Z0/Z'0: redefine process name.
  If (mstp(44)==1) Then
    proc(141) = 'f + fb -> gamma*'
  Else If (mstp(44)==2) Then
    proc(141) = 'f + fb -> Z0'
  Else If (mstp(44)==3) Then
    proc(141) = 'f + fb -> Z''0'
  Else If (mstp(44)==4) Then
    proc(141) = 'f + fb -> gamma*/Z0'
  Else If (mstp(44)==5) Then
    proc(141) = 'f + fb -> gamma*/Z''0'
  Else If (mstp(44)==6) Then
    proc(141) = 'f + fb -> Z0/Z''0'
  Else If (mstp(44)==7) Then
    proc(141) = 'f + fb -> gamma*/Z0/Z''0'
  End If

  Return
End Subroutine pyinre

!*********************************************************************

Subroutine pyxtot

!...Parametrizes total, double diffractive, single diffractive and
!...elastic cross-sections for different energies and beams.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension bcs(5, 8), bcb(2, 5), bcc(3)

!...The following data lines are coefficients needed in the
!...Block, Cahn parametrization of total cross-section and nuclear
!...slope parameter; see below.
  Data ((bcs(i,j),j=1,8), i=1, 5)/41.74, 0.66, 0.0000, 337., 0.0, 0.0, -39.3, 0.48, 41.66, 0.60, 0.0000, 306., 0.0, 0.0, -34.6, 0.51, 41.36, 0.63, 0.0000, 299., 7.3, 0.5, -40.4, 0.47, 41.68, 0.63, 0.0083, 330., 0.0, 0.0, -39.0, 0.48, 41.13, 0.59, 0.0074, 278., 10.5, 0.5, -41.2, 0.46/
  Data ((bcb(i,j),j=1,5), i=1, 2)/10.79, -0.049, 0.040, 21.5, 1.23, 9.92, -0.027, 0.013, 18.9, 1.07/
  Data bcc/2.0164346, -0.5590311, 0.0376279/

!...Total cross-section and nuclear slope parameter for pp and p-pbar
  nfit = min(5, max(1,mstp(31)))
  sigp = bcs(nfit, 1) + bcs(nfit, 2)*(-0.25*paru(1)**2*(1.-0.25*bcs(nfit,3)*paru(1)**2)+(1.+0.5*bcs(nfit,3)*paru(1)**2)*(log(vint(2)/bcs(nfit,4)))**2+bcs(nfit,3)*(log(vint(2)/bcs(nfit,4)))**4)/((1.-0.25*bcs(nfit,3)*paru(1)**2)**2+2.*bcs(nfit,3)*(1.+0.25*bcs(nfit,3)*paru(1)**2)*(log(vint(2)/bcs(nfit,4)))**2+bcs(nfit,3)**2*(log(vint(2)/bcs(nfit,4)))**4) + bcs(nfit, 5)*vint(2)**(bcs(nfit,6)-1.)*sin(0.5*paru(1)*bcs(nfit,6))
  sigm = -bcs(nfit, 7)*vint(2)**(bcs(nfit,8)-1.)*cos(0.5*paru(1)*bcs(nfit,8))
  refp = bcs(nfit, 2)*paru(1)*log(vint(2)/bcs(nfit,4))/((1.-0.25*bcs(nfit,3)*paru(1)**2)**2+2.*bcs(nfit,3)*(1.+0.25*bcs(nfit,3)*paru(1)**2)+(log(vint(2)/bcs(nfit,4)))**2+bcs(nfit,3)**2*(log(vint(2)/bcs(nfit,4)))**4) - bcs(nfit, 5)*vint(2)**(bcs(nfit,6)-1.)*cos(0.5*paru(1)*bcs(nfit,6))
  refm = -bcs(nfit, 7)*vint(2)**(bcs(nfit,8)-1.)*sin(0.5*paru(1)*bcs(nfit,8))
  sigma = sigp - isign(1, mint(11)*mint(12))*sigm
  rho = (refp-isign(1,mint(11)*mint(12))*refm)/sigma

!...Nuclear slope parameter B, curvature C:
  nfit = 1
  If (mstp(31)>=4) nfit = 2
  bp = bcb(nfit, 1) + bcb(nfit, 2)*log(vint(2)) + bcb(nfit, 3)*(log(vint(2)))**2
  bm = bcb(nfit, 4) + bcb(nfit, 5)*log(vint(2))
  b = bp - isign(1, mint(11)*mint(12))*sigm/sigp*(bm-bp)
  vint(121) = b
  c = -0.5*bcc(2)/bcc(3)*(1.-sqrt(max(0.,1.+4.*bcc(3)/bcc(2)**2*(1.E-03*vint(1)-bcc(1)))))
  vint(122) = c

!...Elastic scattering cross-section (fixed by sigma-tot, rho and B).
  sigel = sigma**2*(1.+rho**2)/(16.*paru(1)*paru(5)*b)

!...Single diffractive scattering cross-section from Goulianos:
  sigsd = 2.*0.68*(1.+36./vint(2))*log(0.6+0.1*vint(2))

!...Double diffractive scattering cross-section (essentially fixed by
!...sigma-sd and sigma-el).
  sigdd = sigsd**2/(3.*sigel)

!...Total non-elastic, non-diffractive cross-section.
  signd = sigma - sigdd - sigsd - sigel

!...Rescale for pions.
  If (iabs(mint(11))==211 .And. iabs(mint(12))==211) Then
    sigma = 4./9.*sigma
    sigdd = 4./9.*sigdd
    sigsd = 4./9.*sigsd
    sigel = 4./9.*sigel
    signd = 4./9.*signd
  Else If (iabs(mint(11))==211 .Or. iabs(mint(12))==211) Then
    sigma = 2./3.*sigma
    sigdd = 2./3.*sigdd
    sigsd = 2./3.*sigsd
    sigel = 2./3.*sigel
    signd = 2./3.*signd
  End If

!...Save cross-sections in common block PYPARA.
  vint(101) = sigma
  vint(102) = sigel
  vint(103) = sigsd
  vint(104) = sigdd
  vint(106) = signd
  xsec(95, 1) = signd

  Return
End Subroutine pyxtot

!*********************************************************************

Subroutine pymaxi

!...Finds optimal set of coefficients for kinematical variable selection
!...and the maximum of the part of the differential cross-section used
!...in the event weighting.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Common /pyint6/proc(0:200)
  Character proc*28
  Save /pyint6/
  Character cvar(4)*4
  Dimension npts(4), mvarpt(200, 4), vintpt(200, 30), sigspt(200), narel(6), wtrel(6), wtmat(6, 6), coefu(6), iaccmx(4), sigsmx(4), sigssm(3)
  Data cvar/'tau ', 'tau''', 'y*  ', 'cth '/

!...Select subprocess to study: skip cases not applicable.
  vint(143) = 1.
  vint(144) = 1.
  xsec(0, 1) = 0.
  Do isub = 1, 200
    If (isub>=91 .And. isub<=95) Then
      xsec(isub, 1) = vint(isub+11)
      If (msub(isub)/=1) Goto 350
      Goto 340
    Else If (isub==96) Then
      If (mint(43)/=4) Goto 350
      If (msub(95)/=1 .And. mstp(81)<=0 .And. mstp(131)<=0) Goto 350
    Else If (isub==11 .Or. isub==12 .Or. isub==13 .Or. isub==28 .Or. isub==53 .Or. isub==68) Then
      If (msub(isub)/=1 .Or. msub(95)==1) Goto 350
    Else
      If (msub(isub)/=1) Goto 350
    End If
    mint(1) = isub
    istsb = iset(isub)
    If (isub==96) istsb = 2
    If (mstp(122)>=2) Write (mstu(11), 1000) isub

!...Find resonances (explicit or implicit in cross-section).
    mint(72) = 0
    kfr1 = 0
    If (istsb==1 .Or. istsb==3) Then
      kfr1 = kfpr(isub, 1)
    Else If (isub>=71 .And. isub<=77) Then
      kfr1 = 25
    End If
    If (kfr1/=0) Then
      taur1 = pmas(kfr1, 1)**2/vint(2)
      gamr1 = pmas(kfr1, 1)*pmas(kfr1, 2)/vint(2)
      mint(72) = 1
      mint(73) = kfr1
      vint(73) = taur1
      vint(74) = gamr1
    End If
    If (isub==141) Then
      kfr2 = 23
      taur2 = pmas(kfr2, 1)**2/vint(2)
      gamr2 = pmas(kfr2, 1)*pmas(kfr2, 2)/vint(2)
      mint(72) = 2
      mint(74) = kfr2
      vint(75) = taur2
      vint(76) = gamr2
    End If

!...Find product masses and minimum pT of process.
    sqm3 = 0.
    sqm4 = 0.
    mint(71) = 0
    vint(71) = ckin(3)
    If (istsb==2 .Or. istsb==4) Then
      If (kfpr(isub,1)/=0) sqm3 = pmas(kfpr(isub,1), 1)**2
      If (kfpr(isub,2)/=0) sqm4 = pmas(kfpr(isub,2), 1)**2
      If (min(sqm3,sqm4)<ckin(6)**2) mint(71) = 1
      If (mint(71)==1) vint(71) = max(ckin(3), ckin(5))
      If (isub==96 .And. mstp(82)<=1) vint(71) = parp(81)
      If (isub==96 .And. mstp(82)>=2) vint(71) = 0.08*parp(82)
    End If
    vint(63) = sqm3
    vint(64) = sqm4

!...Number of points for each variable: tau, tau', y*, cos(theta-hat).
    npts(1) = 2 + 2*mint(72)
    If (mint(43)==1 .And. (istsb==1 .Or. istsb==2)) npts(1) = 1
    npts(2) = 1
    If (mint(43)>=2 .And. (istsb==3 .Or. istsb==4)) npts(2) = 2
    npts(3) = 1
    If (mint(43)==4) npts(3) = 3
    npts(4) = 1
    If (istsb==2 .Or. istsb==4) npts(4) = 5
    ntry = npts(1)*npts(2)*npts(3)*npts(4)

!...Reset coefficients of cross-section weighting.
    Do j = 1, 20
      coef(isub, j) = 0.
    End Do
    coef(isub, 1) = 1.
    coef(isub, 7) = 0.5
    coef(isub, 8) = 0.5
    coef(isub, 10) = 1.
    coef(isub, 15) = 1.
    mcth = 0
    mtaup = 0
    cth = 0.
    taup = 0.
    sigsam = 0.

!...Find limits and select tau, y*, cos(theta-hat) and tau' values,
!...in grid of phase space points.
    Call pyklim(1)
    nacc = 0
    Do itry = 1, ntry
      If (mod(itry-1,npts(2)*npts(3)*npts(4))==0) Then
        mtau = 1 + (itry-1)/(npts(2)*npts(3)*npts(4))
        Call pykmap(1, mtau, 0.5)
        If (istsb==3 .Or. istsb==4) Call pyklim(4)
      End If
      If ((istsb==3 .Or. istsb==4) .And. mod(itry-1,npts(3)*npts(4))==0) Then
        mtaup = 1 + mod((itry-1)/(npts(3)*npts(4)), npts(2))
        Call pykmap(4, mtaup, 0.5)
      End If
      If (mod(itry-1,npts(3)*npts(4))==0) Call pyklim(2)
      If (mod(itry-1,npts(4))==0) Then
        myst = 1 + mod((itry-1)/npts(4), npts(3))
        Call pykmap(2, myst, 0.5)
        Call pyklim(3)
      End If
      If (istsb==2 .Or. istsb==4) Then
        mcth = 1 + mod(itry-1, npts(4))
        Call pykmap(3, mcth, 0.5)
      End If
      If (isub==96) vint(25) = vint(21)*(1.-vint(23)**2)

!...Calculate and store cross-section.
      mint(51) = 0
      Call pyklim(0)
      If (mint(51)==1) Goto 120
      nacc = nacc + 1
      mvarpt(nacc, 1) = mtau
      mvarpt(nacc, 2) = mtaup
      mvarpt(nacc, 3) = myst
      mvarpt(nacc, 4) = mcth
      Do j = 1, 30
        vintpt(nacc, j) = vint(10+j)
      End Do
      Call pysigh(nchn, sigs)
      sigspt(nacc) = sigs
      If (sigs>sigsam) sigsam = sigs
      If (mstp(122)>=2) Write (mstu(11), 1100) mtau, mtaup, myst, mcth, vint(21), vint(22), vint(23), vint(26), sigs
    120 End Do
    If (sigsam==0.) Then
      Write (mstu(11), 1200) isub
      Stop
    End If

!...Calculate integrals in tau and y* over maximal phase space limits.
    taumin = vint(11)
    taumax = vint(31)
    atau1 = log(taumax/taumin)
    atau2 = (taumax-taumin)/(taumax*taumin)
    If (npts(1)>=3) Then
      atau3 = log(taumax/taumin*(taumin+taur1)/(taumax+taur1))/taur1
      atau4 = (atan((taumax-taur1)/gamr1)-atan((taumin-taur1)/gamr1))/gamr1
    End If
    If (npts(1)>=5) Then
      atau5 = log(taumax/taumin*(taumin+taur2)/(taumax+taur2))/taur2
      atau6 = (atan((taumax-taur2)/gamr2)-atan((taumin-taur2)/gamr2))/gamr2
    End If
    ystmin = 0.5*log(taumin)
    ystmax = -ystmin
    ayst0 = ystmax - ystmin
    ayst1 = 0.5*(ystmax-ystmin)**2
    ayst3 = 2.*(atan(exp(ystmax))-atan(exp(ystmin)))

!...Reset. Sum up cross-sections in points calculated.
    Do ivar = 1, 4
      If (npts(ivar)==1) Goto 230
      If (isub==96 .And. ivar==4) Goto 230
      nbin = npts(ivar)
      Do j1 = 1, nbin
        narel(j1) = 0
        wtrel(j1) = 0.
        coefu(j1) = 0.
        Do j2 = 1, nbin
          wtmat(j1, j2) = 0.
        End Do
      End Do
      Do iacc = 1, nacc
        ibin = mvarpt(iacc, ivar)
        narel(ibin) = narel(ibin) + 1
        wtrel(ibin) = wtrel(ibin) + sigspt(iacc)

!...Sum up tau cross-section pieces in points used.
        If (ivar==1) Then
          tau = vintpt(iacc, 11)
          wtmat(ibin, 1) = wtmat(ibin, 1) + 1.
          wtmat(ibin, 2) = wtmat(ibin, 2) + (atau1/atau2)/tau
          If (nbin>=3) Then
            wtmat(ibin, 3) = wtmat(ibin, 3) + (atau1/atau3)/(tau+taur1)
            wtmat(ibin, 4) = wtmat(ibin, 4) + (atau1/atau4)*tau/((tau-taur1)**2+gamr1**2)
          End If
          If (nbin>=5) Then
            wtmat(ibin, 5) = wtmat(ibin, 5) + (atau1/atau5)/(tau+taur2)
            wtmat(ibin, 6) = wtmat(ibin, 6) + (atau1/atau6)*tau/((tau-taur2)**2+gamr2**2)
          End If

!...Sum up tau' cross-section pieces in points used.
        Else If (ivar==2) Then
          tau = vintpt(iacc, 11)
          taup = vintpt(iacc, 16)
          taupmn = vintpt(iacc, 6)
          taupmx = vintpt(iacc, 26)
          ataup1 = log(taupmx/taupmn)
          ataup2 = ((1.-tau/taupmx)**4-(1.-tau/taupmn)**4)/(4.*tau)
          wtmat(ibin, 1) = wtmat(ibin, 1) + 1.
          wtmat(ibin, 2) = wtmat(ibin, 2) + (ataup1/ataup2)*(1.-tau/taup)**3/taup

!...Sum up y* and cos(theta-hat) cross-section pieces in points used.
        Else If (ivar==3) Then
          yst = vintpt(iacc, 12)
          wtmat(ibin, 1) = wtmat(ibin, 1) + (ayst0/ayst1)*(yst-ystmin)
          wtmat(ibin, 2) = wtmat(ibin, 2) + (ayst0/ayst1)*(ystmax-yst)
          wtmat(ibin, 3) = wtmat(ibin, 3) + (ayst0/ayst3)/cosh(yst)
        Else
          rm34 = 2.*sqm3*sqm4/(vintpt(iacc,11)*vint(2))**2
          rsqm = 1. + rm34
          cthmax = sqrt(1.-4.*vint(71)**2/(taumax*vint(2)))
          cthmin = -cthmax
          If (cthmax>0.9999) rm34 = max(rm34, 2.*vint(71)**2/(taumax*vint(2)))
          acth1 = cthmax - cthmin
          acth2 = log(max(rm34,rsqm-cthmin)/max(rm34,rsqm-cthmax))
          acth3 = log(max(rm34,rsqm+cthmax)/max(rm34,rsqm+cthmin))
          acth4 = 1./max(rm34, rsqm-cthmax) - 1./max(rm34, rsqm-cthmin)
          acth5 = 1./max(rm34, rsqm+cthmin) - 1./max(rm34, rsqm+cthmax)
          cth = vintpt(iacc, 13)
          wtmat(ibin, 1) = wtmat(ibin, 1) + 1.
          wtmat(ibin, 2) = wtmat(ibin, 2) + (acth1/acth2)/max(rm34, rsqm-cth)
          wtmat(ibin, 3) = wtmat(ibin, 3) + (acth1/acth3)/max(rm34, rsqm+cth)
          wtmat(ibin, 4) = wtmat(ibin, 4) + (acth1/acth4)/max(rm34, rsqm-cth)**2
          wtmat(ibin, 5) = wtmat(ibin, 5) + (acth1/acth5)/max(rm34, rsqm+cth)**2
        End If
      End Do

!...Check that equation system solvable; else trivial way out.
      If (mstp(122)>=2) Write (mstu(11), 1300) cvar(ivar)
      msolv = 1
      Do ibin = 1, nbin
        If (mstp(122)>=2) Write (mstu(11), 1400)(wtmat(ibin,ired), ired=1, nbin), wtrel(ibin)
        If (narel(ibin)==0) msolv = 0
      End Do
      If (msolv==0) Then
        Do ibin = 1, nbin
          coefu(ibin) = 1.
        End Do

!...Solve to find relative importance of cross-section pieces.
      Else
        Do ired = 1, nbin - 1
          Do ibin = ired + 1, nbin
            rqt = wtmat(ibin, ired)/wtmat(ired, ired)
            wtrel(ibin) = wtrel(ibin) - rqt*wtrel(ired)
            Do icoe = ired, nbin
              wtmat(ibin, icoe) = wtmat(ibin, icoe) - rqt*wtmat(ired, icoe)
            End Do
          End Do
        End Do
        Do ired = nbin, 1, -1
          Do icoe = ired + 1, nbin
            wtrel(ired) = wtrel(ired) - wtmat(ired, icoe)*coefu(icoe)
          End Do
          coefu(ired) = wtrel(ired)/wtmat(ired, ired)
        End Do
      End If

!...Normalize coefficients, with piece shared democratically.
      coefsu = 0.
      Do ibin = 1, nbin
        coefu(ibin) = max(0., coefu(ibin))
        coefsu = coefsu + coefu(ibin)
      End Do
      If (ivar==1) ioff = 0
      If (ivar==2) ioff = 14
      If (ivar==3) ioff = 6
      If (ivar==4) ioff = 9
      If (coefsu>0.) Then
        Do ibin = 1, nbin
          coef(isub, ioff+ibin) = parp(121)/nbin + (1.-parp(121))*coefu(ibin)/coefsu
        End Do
      Else
        Do ibin = 1, nbin
          coef(isub, ioff+ibin) = 1./nbin
        End Do
      End If
      If (mstp(122)>=2) Write (mstu(11), 1500) cvar(ivar), (coef(isub,ioff+ibin), ibin=1, nbin)
    230 End Do

!...Find two most promising maxima among points previously determined.
    Do j = 1, 4
      iaccmx(j) = 0
      sigsmx(j) = 0.
    End Do
    nmax = 0
    Do iacc = 1, nacc
      Do j = 1, 30
        vint(10+j) = vintpt(iacc, j)
      End Do
      Call pysigh(nchn, sigs)
      ieq = 0
      Do imv = 1, nmax
        If (abs(sigs-sigsmx(imv))<1E-4*(sigs+sigsmx(imv))) ieq = imv
      End Do
      If (ieq==0) Then
        Do imv = nmax, 1, -1
          iin = imv + 1
          If (sigs<=sigsmx(imv)) Goto 280
          iaccmx(imv+1) = iaccmx(imv)
          sigsmx(imv+1) = sigsmx(imv)
        End Do
        iin = 1
        280 iaccmx(iin) = iacc
        sigsmx(iin) = sigs
        If (nmax<=1) nmax = nmax + 1
      End If
    End Do

!...Read out starting position for search.
    If (mstp(122)>=2) Write (mstu(11), 1600)
    sigsam = sigsmx(1)
    Do imax = 1, nmax
      iacc = iaccmx(imax)
      mtau = mvarpt(iacc, 1)
      mtaup = mvarpt(iacc, 2)
      myst = mvarpt(iacc, 3)
      mcth = mvarpt(iacc, 4)
      vtau = 0.5
      vyst = 0.5
      vcth = 0.5
      vtaup = 0.5

!...Starting point and step size in parameter space.
      Do irpt = 1, 2
        Do ivar = 1, 4
          If (npts(ivar)==1) Goto 310
          If (ivar==1) vvar = vtau
          If (ivar==2) vvar = vtaup
          If (ivar==3) vvar = vyst
          If (ivar==4) vvar = vcth
          If (ivar==1) mvar = mtau
          If (ivar==2) mvar = mtaup
          If (ivar==3) mvar = myst
          If (ivar==4) mvar = mcth
          If (irpt==1) vdel = 0.1
          If (irpt==2) vdel = max(0.01, min(0.05,vvar-0.02,0.98-vvar))
          If (irpt==1) vmar = 0.02
          If (irpt==2) vmar = 0.002
          imov0 = 1
          If (irpt==1 .And. ivar==1) imov0 = 0
          Do imov = imov0, 8

!...Define new point in parameter space.
            If (imov==0) Then
              inew = 2
              vnew = vvar
            Else If (imov==1) Then
              inew = 3
              vnew = vvar + vdel
            Else If (imov==2) Then
              inew = 1
              vnew = vvar - vdel
            Else If (sigssm(3)>=max(sigssm(1),sigssm(2)) .And. vvar+2.*vdel<1.-vmar) Then
              vvar = vvar + vdel
              sigssm(1) = sigssm(2)
              sigssm(2) = sigssm(3)
              inew = 3
              vnew = vvar + vdel
            Else If (sigssm(1)>=max(sigssm(2),sigssm(3)) .And. vvar-2.*vdel>vmar) Then
              vvar = vvar - vdel
              sigssm(3) = sigssm(2)
              sigssm(2) = sigssm(1)
              inew = 1
              vnew = vvar - vdel
            Else If (sigssm(3)>=sigssm(1)) Then
              vdel = 0.5*vdel
              vvar = vvar + vdel
              sigssm(1) = sigssm(2)
              inew = 2
              vnew = vvar
            Else
              vdel = 0.5*vdel
              vvar = vvar - vdel
              sigssm(3) = sigssm(2)
              inew = 2
              vnew = vvar
            End If

!...Convert to relevant variables and find derived new limits.
            If (ivar==1) Then
              vtau = vnew
              Call pykmap(1, mtau, vtau)
              If (istsb==3 .Or. istsb==4) Call pyklim(4)
            End If
            If (ivar<=2 .And. (istsb==3 .Or. istsb==4)) Then
              If (ivar==2) vtaup = vnew
              Call pykmap(4, mtaup, vtaup)
            End If
            If (ivar<=2) Call pyklim(2)
            If (ivar<=3) Then
              If (ivar==3) vyst = vnew
              Call pykmap(2, myst, vyst)
              Call pyklim(3)
            End If
            If (istsb==2 .Or. istsb==4) Then
              If (ivar==4) vcth = vnew
              Call pykmap(3, mcth, vcth)
            End If
            If (isub==96) vint(25) = vint(21)*(1.-vint(23)**2)

!...Evaluate cross-section. Save new maximum. Final maximum.
            Call pysigh(nchn, sigs)
            sigssm(inew) = sigs
            If (sigs>sigsam) sigsam = sigs
            If (mstp(122)>=2) Write (mstu(11), 1700) imax, ivar, mvar, imov, vnew, vint(21), vint(22), vint(23), vint(26), sigs
          End Do
        310 End Do
      End Do
      If (imax==1) sigs11 = sigsam
    End Do
    xsec(isub, 1) = 1.05*sigsam
    340 If (isub/=96) xsec(0, 1) = xsec(0, 1) + xsec(isub, 1)
  350 End Do

!...Print summary table.
  If (mstp(122)>=1) Then
    Write (mstu(11), 1800)
    Write (mstu(11), 1900)
    Do isub = 1, 200
      If (msub(isub)/=1 .And. isub/=96) Goto 360
      If (isub==96 .And. mint(43)/=4) Goto 360
      If (isub==96 .And. msub(95)/=1 .And. mstp(81)<=0) Goto 360
      If (msub(95)==1 .And. (isub==11 .Or. isub==12 .Or. isub==13 .Or. isub==28 .Or. isub==53 .Or. isub==68)) Goto 360
      Write (mstu(11), 2000) isub, proc(isub), xsec(isub, 1)
    360 End Do
    Write (mstu(11), 2100)
  End If

  Return

!...Format statements for maximization results.
  1000 Format (/1X, 'Coefficient optimization and maximum search for ', 'subprocess no', I4/1X, 'Coefficient modes     tau', 10X, 'y*', 9X, 'cth', 9X, 'tau''', 7X, 'sigma')
  1100 Format (1X, 4I4, F12.8, F12.6, F12.7, F12.8, 1P, E12.4)
  1200 Format (1X, 'Error: requested subprocess ', I3, ' has vanishing ', 'cross-section.'/1X, 'Execution stopped!')
  1300 Format (1X, 'Coefficients of equation system to be solved for ', A4)
  1400 Format (1X, 1P, 7E11.3)
  1500 Format (1X, 'Result for ', A4, ':', 6F9.4)
  1600 Format (1X, 'Maximum search for given coefficients'/2X, 'MAX VAR ', 'MOD MOV   VNEW', 7X, 'tau', 7X, 'y*', 8X, 'cth', 7X, 'tau''', 7X, 'sigma')
  1700 Format (1X, 4I4, F8.4, F11.7, F9.3, F11.6, F11.7, 1P, E12.4)
  1800 Format (/1X, 8('*'), 1X, 'PYMAXI: summary of differential ', 'cross-section maximum search', 1X, 8('*'))
  1900 Format (/11X, 58('=')/11X, 'I', 38X, 'I', 17X, 'I'/11X, 'I  ISUB  ', 'Subprocess name', 15X, 'I  Maximum value  I'/11X, 'I', 38X, 'I', 17X, 'I'/11X, 58('=')/11X, 'I', 38X, 'I', 17X, 'I')
  2000 Format (11X, 'I', 2X, I3, 3X, A28, 2X, 'I', 2X, 1P, E12.4, 3X, 'I')
  2100 Format (11X, 'I', 38X, 'I', 17X, 'I'/11X, 58('='))
End Subroutine pymaxi

!*********************************************************************

Subroutine pyovly(movly)

!...Initializes multiplicity distribution and selects mutliplicity
!...of overlayed events, i.e. several events occuring at the same
!...beam crossing.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Dimension wti(0:100)
  Save imax, wti, wts

!...Sum of allowed cross-sections for overlayed events.
  If (movly==1) Then
    vint(131) = vint(106)
    If (mstp(132)>=2) vint(131) = vint(131) + vint(104)
    If (mstp(132)>=3) vint(131) = vint(131) + vint(103)
    If (mstp(132)>=4) vint(131) = vint(131) + vint(102)

!...Initialize multiplicity distribution for unbiased events.
    If (mstp(133)==1) Then
      xnave = vint(131)*parp(131)
      If (xnave>40.) Write (mstu(11), 1000) xnave
      wti(0) = exp(-min(50.,xnave))
      wts = 0.
      wtn = 0.
      Do i = 1, 100
        wti(i) = wti(i-1)*xnave/i
        If (i-2.5>xnave .And. wti(i)<1E-6) Goto 110
        wts = wts + wti(i)
        wtn = wtn + wti(i)*i
        imax = i
      End Do
      110 vint(132) = xnave
      vint(133) = wtn/wts
      vint(134) = wts

!...Initialize mutiplicity distribution for biased events.
    Else If (mstp(133)==2) Then
      xnave = vint(131)*parp(131)
      If (xnave>40.) Write (mstu(11), 1000) xnave
      wti(1) = exp(-min(50.,xnave))*xnave
      wts = wti(1)
      wtn = wti(1)
      Do i = 2, 100
        wti(i) = wti(i-1)*xnave/(i-1)
        If (i-2.5>xnave .And. wti(i)<1E-6) Goto 130
        wts = wts + wti(i)
        wtn = wtn + wti(i)*i
        imax = i
      End Do
      130 vint(132) = xnave
      vint(133) = wtn/wts
      vint(134) = wts
    End If

!...Pick multiplicity of overlayed events.
  Else
    If (mstp(133)==0) Then
      mint(81) = max(1, mstp(134))
    Else
      wtr = wts*rlu(0)
      Do i = 1, imax
        mint(81) = i
        wtr = wtr - wti(i)
        If (wtr<=0.) Goto 150
      End Do
      150 Continue
    End If
  End If

  Return

!...Format statement for error message.
  1000 Format (1X, 'Warning: requested average number of events per bunch', 'crossing too large, ', 1P, E12.4)
End Subroutine pyovly

!*********************************************************************

Subroutine pyrand

!...Generates quantities characterizing the high-pT scattering at the
!...parton level according to the matrix elements. Chooses incoming,
!...reacting partons, their momentum fractions and one of the possible
!...subprocesses.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/

!...Initial values, specifically for (first) semihard interaction.
  mint(17) = 0
  mint(18) = 0
  vint(143) = 1.
  vint(144) = 1.
  If (msub(95)==1 .Or. mint(82)>=2) Call pymult(2)
  isub = 0
  100 mint(51) = 0

!...Choice of process type - first event of overlay.
  If (mint(82)==1 .And. (isub<=90 .Or. isub>96)) Then
    rsub = xsec(0, 1)*rlu(0)
    Do i = 1, 200
      If (msub(i)/=1) Goto 110
      isub = i
      rsub = rsub - xsec(i, 1)
      If (rsub<=0.) Goto 120
    110 End Do
    120 If (isub==95) isub = 96

!...Choice of inclusive process type - overlayed events.
  Else If (mint(82)>=2 .And. isub==0) Then
    rsub = vint(131)*rlu(0)
    isub = 96
    If (rsub>vint(106)) isub = 93
    If (rsub>vint(106)+vint(104)) isub = 92
    If (rsub>vint(106)+vint(104)+vint(103)) isub = 91
  End If
  If (mint(82)==1) ngen(0, 1) = ngen(0, 1) + 1
  If (mint(82)==1) ngen(isub, 1) = ngen(isub, 1) + 1
  mint(1) = isub

!...Find resonances (explicit or implicit in cross-section).
  mint(72) = 0
  kfr1 = 0
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    kfr1 = kfpr(isub, 1)
  Else If (isub>=71 .And. isub<=77) Then
    kfr1 = 25
  End If
  If (kfr1/=0) Then
    taur1 = pmas(kfr1, 1)**2/vint(2)
    gamr1 = pmas(kfr1, 1)*pmas(kfr1, 2)/vint(2)
    mint(72) = 1
    mint(73) = kfr1
    vint(73) = taur1
    vint(74) = gamr1
  End If
  If (isub==141) Then
    kfr2 = 23
    taur2 = pmas(kfr2, 1)**2/vint(2)
    gamr2 = pmas(kfr2, 1)*pmas(kfr2, 2)/vint(2)
    mint(72) = 2
    mint(74) = kfr2
    vint(75) = taur2
    vint(76) = gamr2
  End If

!...Find product masses and minimum pT of process,
!...optionally with broadening according to a truncated Breit-Wigner.
  vint(63) = 0.
  vint(64) = 0.
  mint(71) = 0
  vint(71) = ckin(3)
  If (mint(82)>=2) vint(71) = 0.
  If (iset(isub)==2 .Or. iset(isub)==4) Then
    Do i = 1, 2
      If (kfpr(isub,i)==0) Then
      Else If (mstp(42)<=0) Then
        vint(62+i) = pmas(kfpr(isub,i), 1)**2
      Else
        vint(62+i) = ulmass(kfpr(isub,i))**2
      End If
    End Do
    If (min(vint(63),vint(64))<ckin(6)**2) mint(71) = 1
    If (mint(71)==1) vint(71) = max(ckin(3), ckin(5))
  End If

  If (iset(isub)==0) Then
!...Double or single diffractive, or elastic scattering:
!...choose m^2 according to 1/m^2 (diffractive), constant (elastic)
    is = int(1.5+rlu(0))
    vint(63) = vint(3)**2
    vint(64) = vint(4)**2
    If (isub==92 .Or. isub==93) vint(62+is) = parp(111)**2
    If (isub==93) vint(65-is) = parp(111)**2
    sh = vint(2)
    sqm1 = vint(3)**2
    sqm2 = vint(4)**2
    sqm3 = vint(63)
    sqm4 = vint(64)
    sqla12 = (sh-sqm1-sqm2)**2 - 4.*sqm1*sqm2
    sqla34 = (sh-sqm3-sqm4)**2 - 4.*sqm3*sqm4
    thter1 = sqm1 + sqm2 + sqm3 + sqm4 - (sqm1-sqm2)*(sqm3-sqm4)/sh - sh
    thter2 = sqrt(max(0.,sqla12))*sqrt(max(0.,sqla34))/sh
    thl = 0.5*(thter1-thter2)
    thu = 0.5*(thter1+thter2)
    thm = min(max(thl,parp(101)), thu)
    jtmax = 0
    If (isub==92 .Or. isub==93) jtmax = isub - 91
    Do jt = 1, jtmax
      mint(13+3*jt-is*(2*jt-3)) = 1
      sqmmin = vint(59+3*jt-is*(2*jt-3))
      sqmi = vint(8-3*jt+is*(2*jt-3))**2
      sqmj = vint(3*jt-1-is*(2*jt-3))**2
      sqmf = vint(68-3*jt+is*(2*jt-3))
      squa = 0.5*sh/sqmi*((1.+(sqmi-sqmj)/sh)*thm+sqmi-sqmf-sqmj**2/sh+(sqmi+sqmj)*sqmf/sh+(sqmi-sqmj)**2/sh**2*sqmf)
      quar = sh/sqmi*(thm*(thm+sh-sqmi-sqmj-sqmf*(1.-(sqmi-sqmj)/sh))+sqmi*sqmj-sqmj*sqmf*(1.+(sqmi-sqmj-sqmf)/sh))
      sqmmax = squa + sqrt(max(0.,squa**2-quar))
      If (abs(quar/squa**2)<1.E-06) sqmmax = 0.5*quar/squa
      sqmmax = min(sqmmax, (vint(1)-sqrt(sqmf))**2)
      vint(59+3*jt-is*(2*jt-3)) = sqmmin*(sqmmax/sqmmin)**rlu(0)
    End Do
!...Choose t-hat according to exp(B*t-hat+C*t-hat^2).
    sqm3 = vint(63)
    sqm4 = vint(64)
    sqla34 = (sh-sqm3-sqm4)**2 - 4.*sqm3*sqm4
    thter1 = sqm1 + sqm2 + sqm3 + sqm4 - (sqm1-sqm2)*(sqm3-sqm4)/sh - sh
    thter2 = sqrt(max(0.,sqla12))*sqrt(max(0.,sqla34))/sh
    thl = 0.5*(thter1-thter2)
    thu = 0.5*(thter1+thter2)
    b = vint(121)
    c = vint(122)
    If (isub==92 .Or. isub==93) Then
      b = 0.5*b
      c = 0.5*c
    End If
    thm = min(max(thl,parp(101)), thu)
    expth = 0.
    tharg = b*(thm-thu)
    If (tharg>-20.) expth = exp(tharg)
    150 th = thu + log(expth+(1.-expth)*rlu(0))/b
    th = max(thm, min(thu,th))
    ratlog = min((b+c*(th+thm))*(th-thm), (b+c*(th+thu))*(th-thu))
    If (ratlog<log(rlu(0))) Goto 150
    vint(21) = 1.
    vint(22) = 0.
    vint(23) = min(1., max(-1.,(2.*th-thter1)/thter2))

!...Note: in the following, by In is meant the integral over the
!...quantity multiplying coefficient cn.
!...Choose tau according to h1(tau)/tau, where
!...h1(tau) = c0 + I0/I1*c1*1/tau + I0/I2*c2*1/(tau+tau_R) +
!...I0/I3*c3*tau/((s*tau-m^2)^2+(m*Gamma)^2) +
!...I0/I4*c4*1/(tau+tau_R') +
!...I0/I5*c5*tau/((s*tau-m'^2)^2+(m'*Gamma')^2), and
!...c0 + c1 + c2 + c3 + c4 + c5 = 1
  Else If (iset(isub)>=1 .And. iset(isub)<=4) Then
    Call pyklim(1)
    If (mint(51)/=0) Goto 100
    rtau = rlu(0)
    mtau = 1
    If (rtau>coef(isub,1)) mtau = 2
    If (rtau>coef(isub,1)+coef(isub,2)) mtau = 3
    If (rtau>coef(isub,1)+coef(isub,2)+coef(isub,3)) mtau = 4
    If (rtau>coef(isub,1)+coef(isub,2)+coef(isub,3)+coef(isub,4)) mtau = 5
    If (rtau>coef(isub,1)+coef(isub,2)+coef(isub,3)+coef(isub,4)+coef(isub,5)) mtau = 6
    Call pykmap(1, mtau, rlu(0))

!...2 -> 3, 4 processes:
!...Choose tau' according to h4(tau,tau')/tau', where
!...h4(tau,tau') = c0 + I0/I1*c1*(1 - tau/tau')^3/tau', and
!...c0 + c1 = 1.
    If (iset(isub)==3 .Or. iset(isub)==4) Then
      Call pyklim(4)
      If (mint(51)/=0) Goto 100
      rtaup = rlu(0)
      mtaup = 1
      If (rtaup>coef(isub,15)) mtaup = 2
      Call pykmap(4, mtaup, rlu(0))
    End If

!...Choose y* according to h2(y*), where
!...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +
!...I0/I3*c3*1/cosh(y*), I0 = y*max-y*min, and c1 + c2 + c3 = 1.
    Call pyklim(2)
    If (mint(51)/=0) Goto 100
    ryst = rlu(0)
    myst = 1
    If (ryst>coef(isub,7)) myst = 2
    If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
    Call pykmap(2, myst, rlu(0))

!...2 -> 2 processes:
!...Choose cos(theta-hat) (cth) according to h3(cth), where
!...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +
!...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,
!...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products),
!...and c0 + c1 + c2 + c3 + c4 = 1.
    Call pyklim(3)
    If (mint(51)/=0) Goto 100
    If (iset(isub)==2 .Or. iset(isub)==4) Then
      rcth = rlu(0)
      mcth = 1
      If (rcth>coef(isub,10)) mcth = 2
      If (rcth>coef(isub,10)+coef(isub,11)) mcth = 3
      If (rcth>coef(isub,10)+coef(isub,11)+coef(isub,12)) mcth = 4
      If (rcth>coef(isub,10)+coef(isub,11)+coef(isub,12)+coef(isub,13)) mcth = 5
      Call pykmap(3, mcth, rlu(0))
    End If

!...Low-pT or multiple interactions (first semihard interaction).
  Else If (iset(isub)==5) Then
    Call pymult(3)
    isub = mint(1)
  End If

!...Choose azimuthal angle.
  vint(24) = paru(2)*rlu(0)

!...Check against user cuts on kinematics at parton level.
  mint(51) = 0
  If (isub<=90 .Or. isub>100) Call pyklim(0)
  If (mint(51)/=0) Goto 100
  If (mint(82)==1 .And. mstp(141)>=1) Then
    mcut = 0
    If (msub(91)+msub(92)+msub(93)+msub(94)+msub(95)==0) Call pykcut(mcut)
    If (mcut/=0) Goto 100
  End If

!...Calculate differential cross-section for different subprocesses.
  Call pysigh(nchn, sigs)

!...Calculations for Monte Carlo estimate of all cross-sections.
  If (mint(82)==1 .And. isub<=90 .Or. isub>=96) Then
    xsec(isub, 2) = xsec(isub, 2) + sigs
  Else If (mint(82)==1) Then
    xsec(isub, 2) = xsec(isub, 2) + xsec(isub, 1)
  End If

!...Multiple interactions: store results of cross-section calculation.
  If (mint(43)==4 .And. mstp(82)>=3) Then
    vint(153) = sigs
    Call pymult(4)
  End If

!...Weighting using estimate of maximum of differential cross-section.
  viol = sigs/xsec(isub, 1)
  If (viol<rlu(0)) Goto 100

!...Check for possible violation of estimated maximum of differential
!...cross-section used in weighting.
  If (mstp(123)<=0) Then
    If (viol>1.) Then
      Write (mstu(11), 1000) viol, ngen(0, 3) + 1
      Write (mstu(11), 1100) isub, vint(21), vint(22), vint(23), vint(26)
      Stop
    End If
  Else If (mstp(123)==1) Then
    If (viol>vint(108)) Then
      vint(108) = viol
!          IF(VIOL.GT.1.) THEN
!            WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1
!            WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),
!     &      VINT(26)
!          ENDIF
    End If
  Else If (viol>vint(108)) Then
    vint(108) = viol
    If (viol>1.) Then
      xdif = xsec(isub, 1)*(viol-1.)
      xsec(isub, 1) = xsec(isub, 1) + xdif
      If (msub(isub)==1 .And. (isub<=90 .Or. isub>96)) xsec(0, 1) = xsec(0, 1) + xdif
!          WRITE(MSTU(11),1200) VIOL,NGEN(0,3)+1
!          WRITE(MSTU(11),1100) ISUB,VINT(21),VINT(22),VINT(23),VINT(26)
!          IF(ISUB.LE.9) THEN
!            WRITE(MSTU(11),1300) ISUB,XSEC(ISUB,1)
!          ELSEIF(ISUB.LE.99) THEN
!            WRITE(MSTU(11),1400) ISUB,XSEC(ISUB,1)
!          ELSE
!            WRITE(MSTU(11),1500) ISUB,XSEC(ISUB,1)
!          ENDIF
      vint(108) = 1.
    End If
  End If

!...Multiple interactions: choose impact parameter.
  vint(148) = 1.
  If (mint(43)==4 .And. (isub<=90 .Or. isub>=96) .And. mstp(82)>=3) Then
    Call pymult(5)
    If (vint(150)<rlu(0)) Goto 100
  End If
  If (mint(82)==1 .And. msub(95)==1) Then
    If (isub<=90 .Or. isub>=95) ngen(95, 1) = ngen(95, 1) + 1
    If (isub<=90 .Or. isub>=96) ngen(96, 2) = ngen(96, 2) + 1
  End If
  If (isub<=90 .Or. isub>=96) mint(31) = mint(31) + 1

!...Choose flavour of reacting partons (and subprocess).
  rsigs = sigs*rlu(0)
  qt2 = vint(48)
  rqqbar = parp(87)*(1.-(qt2/(qt2+(parp(88)*parp(82))**2))**2)
  If (isub/=95 .And. (isub/=96 .Or. mstp(82)<=1 .Or. rlu(0)>rqqbar)) Then
    Do ichn = 1, nchn
      kfl1 = isig(ichn, 1)
      kfl2 = isig(ichn, 2)
      mint(2) = isig(ichn, 3)
      rsigs = rsigs - sigh(ichn)
      If (rsigs<=0.) Goto 210
    End Do

!...Multiple interactions: choose qqbar preferentially at small pT.
  Else If (isub==96) Then
    Call pyspli(mint(11), 21, kfl1, kfldum)
    Call pyspli(mint(12), 21, kfl2, kfldum)
    mint(1) = 11
    mint(2) = 1
    If (kfl1==kfl2 .And. rlu(0)<0.5) mint(2) = 2

!...Low-pT: choose string drawing configuration.
  Else
    kfl1 = 21
    kfl2 = 21
    rsigs = 6.*rlu(0)
    mint(2) = 1
    If (rsigs>1.) mint(2) = 2
    If (rsigs>2.) mint(2) = 3
  End If

!...Reassign QCD process. Partons before initial state radiation.
  210 If (mint(2)>10) Then
    mint(1) = mint(2)/10
    mint(2) = mod(mint(2), 10)
  End If
  mint(15) = kfl1
  mint(16) = kfl2
  mint(13) = mint(15)
  mint(14) = mint(16)
  vint(141) = vint(41)
  vint(142) = vint(42)
!lin 1200 FORMAT(1X,'Warning: maximum violated by',1P,E11.3,1X,
!     &'in event',1X,I7)
! 1300 FORMAT(1X,'XSEC(',I1,',1) increased to',1P,E11.3)
! 1400 FORMAT(1X,'XSEC(',I2,',1) increased to',1P,E11.3)
!lin 1500 FORMAT(1X,'XSEC(',I3,',1) increased to',1P,E11.3)

  Return

!...Format statements for differential cross-section maximum violations.
  1000 Format (1X, 'Error: maximum violated by', 1P, E11.3, 1X, 'in event', 1X, I7, '.'/1X, 'Execution stopped!')
  1100 Format (1X, 'ISUB = ', I3, '; Point of violation:'/1X, 'tau=', 1P, E11.3, ', y* =', E11.3, ', cthe = ', 0P, F11.7, ', tau'' =', 1P, E11.3)
End Subroutine pyrand

!*********************************************************************

Subroutine pyscat

!...Finds outgoing flavours and event type; sets up the kinematics
!...and colour flow of the hard scattering.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension wdtp(0:40), wdte(0:40, 0:5), pmq(2), z(2), cthe(2), phi(2)

!...Choice of subprocess, number of documentation lines.
  isub = mint(1)
  idoc = 6 + iset(isub)
  If (isub==95) idoc = 8
  mint(3) = idoc - 6
  If (idoc>=9) idoc = idoc + 2
  mint(4) = idoc
  ipu1 = mint(84) + 1
  ipu2 = mint(84) + 2
  ipu3 = mint(84) + 3
  ipu4 = mint(84) + 4
  ipu5 = mint(84) + 5
  ipu6 = mint(84) + 6

!...Reset K, P and V vectors. Store incoming particles.
  Do jt = 1, mstp(126) + 10
    i = mint(83) + jt
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
  End Do
  Do jt = 1, 2
    i = mint(83) + jt
    k(i, 1) = 21
    k(i, 2) = mint(10+jt)
    p(i, 1) = 0.
    p(i, 2) = 0.
    p(i, 5) = vint(2+jt)
    p(i, 3) = vint(5)*(-1)**(jt+1)
    p(i, 4) = sqrt(p(i,3)**2+p(i,5)**2)
  End Do
  mint(6) = 2
  kfres = 0

!...Store incoming partons in their CM-frame.
  sh = vint(44)
  shr = sqrt(sh)
  shp = vint(26)*vint(2)
  shpr = sqrt(shp)
  shuser = shr
  If (iset(isub)>=3) shuser = shpr
  Do jt = 1, 2
    i = mint(84) + jt
    k(i, 1) = 14
    k(i, 2) = mint(14+jt)
    k(i, 3) = mint(83) + 2 + jt
    p(i, 5) = ulmass(k(i,2))
  End Do
  If (p(ipu1,5)+p(ipu2,5)>=shuser) Then
    p(ipu1, 5) = 0.
    p(ipu2, 5) = 0.
  End If
  p(ipu1, 4) = 0.5*(shuser+(p(ipu1,5)**2-p(ipu2,5)**2)/shuser)
  p(ipu1, 3) = sqrt(max(0.,p(ipu1,4)**2-p(ipu1,5)**2))
  p(ipu2, 4) = shuser - p(ipu1, 4)
  p(ipu2, 3) = -p(ipu1, 3)

!...Copy incoming partons to documentation lines.
  Do jt = 1, 2
    i1 = mint(83) + 4 + jt
    i2 = mint(84) + jt
    k(i1, 1) = 21
    k(i1, 2) = k(i2, 2)
    k(i1, 3) = i1 - 2
    Do j = 1, 5
      p(i1, j) = p(i2, j)
    End Do
  End Do

!...Choose new quark flavour for relevant annihilation graphs.
  If (isub==12 .Or. isub==53) Then
    Call pywidt(21, shr, wdtp, wdte)
    rkfl = (wdte(0,1)+wdte(0,2)+wdte(0,4))*rlu(0)
    Do i = 1, 2*mstp(1)
      kflq = i
      rkfl = rkfl - (wdte(i,1)+wdte(i,2)+wdte(i,4))
      If (rkfl<=0.) Goto 150
    End Do
    150 Continue
  End If

!...Final state flavours and colour flow: default values.
  js = 1
  mint(21) = mint(15)
  mint(22) = mint(16)
  mint(23) = 0
  mint(24) = 0
  kcc = 20
  kcs = isign(1, mint(15))

  If (isub<=10) Then
    If (isub==1) Then
!...f + fb -> gamma*/Z0.
      kfres = 23

    Else If (isub==2) Then
!...f + fb' -> W+/- .
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      kfres = isign(24, kch1+kch2)

    Else If (isub==3) Then
!...f + fb -> H0.
      kfres = 25

    Else If (isub==4) Then
!...gamma + W+/- -> W+/-.

    Else If (isub==5) Then
!...Z0 + Z0 -> H0.
      xh = sh/shp
      mint(21) = mint(15)
      mint(22) = mint(16)
      pmq(1) = ulmass(mint(21))
      pmq(2) = ulmass(mint(22))
      240 jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 240
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 240
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 240
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 240
      kcc = 22
      kfres = 25

    Else If (isub==6) Then
!...Z0 + W+/- -> W+/-.

    Else If (isub==7) Then
!...W+ + W- -> Z0.

    Else If (isub==8) Then
!...W+ + W- -> H0.
      xh = sh/shp
      250 Do jt = 1, 2
        i = mint(14+jt)
        ia = iabs(i)
        If (ia<=10) Then
          rvckm = vint(180+i)*rlu(0)
          Do j = 1, mstp(1)
            ib = 2*j - 1 + mod(ia, 2)
            ipm = (5-isign(1,i))/2
            idc = j + mdcy(ia, 2) + 2
            If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 270
            mint(20+jt) = isign(ib, i)
            rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
            If (rvckm<=0.) Goto 280
          270 End Do
        Else
          ib = 2*((ia+1)/2) - 1 + mod(ia, 2)
          mint(20+jt) = isign(ib, i)
        End If
        280 pmq(jt) = ulmass(mint(20+jt))
      End Do
      jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 250
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 250
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 250
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 250
      kcc = 22
      kfres = 25
    End If

  Else If (isub<=20) Then
    If (isub==11) Then
!...f + f' -> f + f'; th = (p(f)-p(f))**2.
      kcc = mint(2)
      If (mint(15)*mint(16)<0) kcc = kcc + 2

    Else If (isub==12) Then
!...f + fb -> f' + fb'; th = (p(f)-p(f'))**2.
      mint(21) = isign(kflq, mint(15))
      mint(22) = -mint(21)
      kcc = 4

    Else If (isub==13) Then
!...f + fb -> g + g; th arbitrary.
      mint(21) = 21
      mint(22) = 21
      kcc = mint(2) + 4

    Else If (isub==14) Then
!...f + fb -> g + gam; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 22
      kcc = 17 + js

    Else If (isub==15) Then
!...f + fb -> g + Z0; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 23
      kcc = 17 + js

    Else If (isub==16) Then
!...f + fb' -> g + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)<0) js = 2
      mint(20+js) = 21
      mint(23-js) = isign(24, kch1+kch2)
      kcc = 17 + js

    Else If (isub==17) Then
!...f + fb -> g + H0; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 25
      kcc = 17 + js

    Else If (isub==18) Then
!...f + fb -> gamma + gamma; th arbitrary.
      mint(21) = 22
      mint(22) = 22

    Else If (isub==19) Then
!...f + fb -> gamma + Z0; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 22
      mint(23-js) = 23

    Else If (isub==20) Then
!...f + fb' -> gamma + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)<0) js = 2
      mint(20+js) = 22
      mint(23-js) = isign(24, kch1+kch2)
    End If

  Else If (isub<=30) Then
    If (isub==21) Then
!...f + fb -> gamma + H0; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 22
      mint(23-js) = 25

    Else If (isub==22) Then
!...f + fb -> Z0 + Z0; th arbitrary.
      mint(21) = 23
      mint(22) = 23

    Else If (isub==23) Then
!...f + fb' -> Z0 + W+/-; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)<0) js = 2
      mint(20+js) = 23
      mint(23-js) = isign(24, kch1+kch2)

    Else If (isub==24) Then
!...f + fb -> Z0 + H0; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 23
      mint(23-js) = 25

    Else If (isub==25) Then
!...f + fb -> W+ + W-; th = (p(f)-p(W-))**2.
      mint(21) = -isign(24, mint(15))
      mint(22) = -mint(21)

    Else If (isub==26) Then
!...f + fb' -> W+/- + H0; th = (p(f)-p(W-))**2 or (p(fb')-p(W+))**2.
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      If (mint(15)*(kch1+kch2)>0) js = 2
      mint(20+js) = isign(24, kch1+kch2)
      mint(23-js) = 25

    Else If (isub==27) Then
!...f + fb -> H0 + H0.

    Else If (isub==28) Then
!...f + g -> f + g; th = (p(f)-p(f))**2.
      kcc = mint(2) + 6
      If (mint(15)==21) kcc = kcc + 2
      If (mint(15)/=21) kcs = isign(1, mint(15))
      If (mint(16)/=21) kcs = isign(1, mint(16))

    Else If (isub==29) Then
!...f + g -> f + gamma; th = (p(f)-p(f))**2.
      If (mint(15)==21) js = 2
      mint(23-js) = 22
      kcc = 15 + js
      kcs = isign(1, mint(14+js))

    Else If (isub==30) Then
!...f + g -> f + Z0; th = (p(f)-p(f))**2.
      If (mint(15)==21) js = 2
      mint(23-js) = 23
      kcc = 15 + js
      kcs = isign(1, mint(14+js))
    End If

  Else If (isub<=40) Then
    If (isub==31) Then
!...f + g -> f' + W+/-; th = (p(f)-p(f'))**2; choose flavour f'.
      If (mint(15)==21) js = 2
      i = mint(14+js)
      ia = iabs(i)
      mint(23-js) = isign(24, kchg(ia,1)*i)
      rvckm = vint(180+i)*rlu(0)
      Do j = 1, mstp(1)
        ib = 2*j - 1 + mod(ia, 2)
        ipm = (5-isign(1,i))/2
        idc = j + mdcy(ia, 2) + 2
        If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 220
        mint(20+js) = isign(ib, i)
        rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
        If (rvckm<=0.) Goto 230
      220 End Do
      230 kcc = 15 + js
      kcs = isign(1, mint(14+js))

    Else If (isub==32) Then
!...f + g -> f + H0; th = (p(f)-p(f))**2.
      If (mint(15)==21) js = 2
      mint(23-js) = 25
      kcc = 15 + js
      kcs = isign(1, mint(14+js))

    Else If (isub==33) Then
!...f + gamma -> f + g.

    Else If (isub==34) Then
!...f + gamma -> f + gamma.

    Else If (isub==35) Then
!...f + gamma -> f + Z0.

    Else If (isub==36) Then
!...f + gamma -> f' + W+/-.

    Else If (isub==37) Then
!...f + gamma -> f + H0.

    Else If (isub==38) Then
!...f + Z0 -> f + g.

    Else If (isub==39) Then
!...f + Z0 -> f + gamma.

    Else If (isub==40) Then
!...f + Z0 -> f + Z0.
    End If

  Else If (isub<=50) Then
    If (isub==41) Then
!...f + Z0 -> f' + W+/-.

    Else If (isub==42) Then
!...f + Z0 -> f + H0.

    Else If (isub==43) Then
!...f + W+/- -> f' + g.

    Else If (isub==44) Then
!...f + W+/- -> f' + gamma.

    Else If (isub==45) Then
!...f + W+/- -> f' + Z0.

    Else If (isub==46) Then
!...f + W+/- -> f' + W+/-.

    Else If (isub==47) Then
!...f + W+/- -> f' + H0.

    Else If (isub==48) Then
!...f + H0 -> f + g.

    Else If (isub==49) Then
!...f + H0 -> f + gamma.

    Else If (isub==50) Then
!...f + H0 -> f + Z0.
    End If

  Else If (isub<=60) Then
    If (isub==51) Then
!...f + H0 -> f' + W+/-.

    Else If (isub==52) Then
!...f + H0 -> f + H0.

    Else If (isub==53) Then
!...g + g -> f + fb; th arbitrary.
      kcs = (-1)**int(1.5+rlu(0))
      mint(21) = isign(kflq, kcs)
      mint(22) = -mint(21)
      kcc = mint(2) + 10

    Else If (isub==54) Then
!...g + gamma -> f + fb.

    Else If (isub==55) Then
!...g + Z0 -> f + fb.

    Else If (isub==56) Then
!...g + W+/- -> f + fb'.

    Else If (isub==57) Then
!...g + H0 -> f + fb.

    Else If (isub==58) Then
!...gamma + gamma -> f + fb.

    Else If (isub==59) Then
!...gamma + Z0 -> f + fb.

    Else If (isub==60) Then
!...gamma + W+/- -> f + fb'.
    End If

  Else If (isub<=70) Then
    If (isub==61) Then
!...gamma + H0 -> f + fb.

    Else If (isub==62) Then
!...Z0 + Z0 -> f + fb.

    Else If (isub==63) Then
!...Z0 + W+/- -> f + fb'.

    Else If (isub==64) Then
!...Z0 + H0 -> f + fb.

    Else If (isub==65) Then
!...W+ + W- -> f + fb.

    Else If (isub==66) Then
!...W+/- + H0 -> f + fb'.

    Else If (isub==67) Then
!...H0 + H0 -> f + fb.

    Else If (isub==68) Then
!...g + g -> g + g; th arbitrary.
      kcc = mint(2) + 12
      kcs = (-1)**int(1.5+rlu(0))

    Else If (isub==69) Then
!...gamma + gamma -> W+ + W-.

    Else If (isub==70) Then
!...gamma + W+/- -> gamma + W+/-
    End If

  Else If (isub<=80) Then
    If (isub==71 .Or. isub==72) Then
!...Z0 + Z0 -> Z0 + Z0; Z0 + Z0 -> W+ + W-.
      xh = sh/shp
      mint(21) = mint(15)
      mint(22) = mint(16)
      pmq(1) = ulmass(mint(21))
      pmq(2) = ulmass(mint(22))
      290 jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 290
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 290
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 290
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 290
      kcc = 22

    Else If (isub==73) Then
!...Z0 + W+/- -> Z0 + W+/-.
      xh = sh/shp
      300 jt = int(1.5+rlu(0))
      i = mint(14+jt)
      ia = iabs(i)
      If (ia<=10) Then
        rvckm = vint(180+i)*rlu(0)
        Do j = 1, mstp(1)
          ib = 2*j - 1 + mod(ia, 2)
          ipm = (5-isign(1,i))/2
          idc = j + mdcy(ia, 2) + 2
          If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 320
          mint(20+jt) = isign(ib, i)
          rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
          If (rvckm<=0.) Goto 330
        320 End Do
      Else
        ib = 2*((ia+1)/2) - 1 + mod(ia, 2)
        mint(20+jt) = isign(ib, i)
      End If
      330 pmq(jt) = ulmass(mint(20+jt))
      mint(23-jt) = mint(17-jt)
      pmq(3-jt) = ulmass(mint(23-jt))
      jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 300
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 300
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 300
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(23,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 300
      kcc = 22

    Else If (isub==74) Then
!...Z0 + H0 -> Z0 + H0.

    Else If (isub==75) Then
!...W+ + W- -> gamma + gamma.

    Else If (isub==76 .Or. isub==77) Then
!...W+ + W- -> Z0 + Z0; W+ + W- -> W+ + W-.
      xh = sh/shp
      340 Do jt = 1, 2
        i = mint(14+jt)
        ia = iabs(i)
        If (ia<=10) Then
          rvckm = vint(180+i)*rlu(0)
          Do j = 1, mstp(1)
            ib = 2*j - 1 + mod(ia, 2)
            ipm = (5-isign(1,i))/2
            idc = j + mdcy(ia, 2) + 2
            If (mdme(idc,1)/=1 .And. mdme(idc,1)/=ipm) Goto 360
            mint(20+jt) = isign(ib, i)
            rvckm = rvckm - vckm((ia+1)/2, (ib+1)/2)
            If (rvckm<=0.) Goto 370
          360 End Do
        Else
          ib = 2*((ia+1)/2) - 1 + mod(ia, 2)
          mint(20+jt) = isign(ib, i)
        End If
        370 pmq(jt) = ulmass(mint(20+jt))
      End Do
      jt = int(1.5+rlu(0))
      zmin = 2.*pmq(jt)/shpr
      zmax = 1. - pmq(3-jt)/shpr - (sh-pmq(jt)**2)/(shpr*(shpr-pmq(3-jt)))
      zmax = min(1.-xh, zmax)
      z(jt) = zmin + (zmax-zmin)*rlu(0)
      If (-1.+(1.+xh)/(1.-z(jt))-xh/(1.-z(jt))**2<(1.-xh)**2/(4.*xh)*rlu(0)) Goto 340
      sqc1 = 1. - 4.*pmq(jt)**2/(z(jt)**2*shp)
      If (sqc1<1.E-8) Goto 340
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(jt)**2)/(z(jt)*shp)
      cthe(jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(jt) = min(1., max(-1.,cthe(jt)))
      z(3-jt) = 1. - xh/(1.-z(jt))
      sqc1 = 1. - 4.*pmq(3-jt)**2/(z(3-jt)**2*shp)
      If (sqc1<1.E-8) Goto 340
      c1 = sqrt(sqc1)
      c2 = 1. + 2.*(pmas(24,1)**2-pmq(3-jt)**2)/(z(3-jt)*shp)
      cthe(3-jt) = (c2-(c2**2-c1**2)/(c2+(2.*rlu(0)-1.)*c1))/c1
      cthe(3-jt) = min(1., max(-1.,cthe(3-jt)))
      phir = paru(2)*rlu(0)
      cphi = cos(phir)
      ang = cthe(1)*cthe(2) - sqrt(1.-cthe(1)**2)*sqrt(1.-cthe(2)**2)*cphi
      z1 = 2. - z(jt)
      z2 = ang*sqrt(z(jt)**2-4.*pmq(jt)**2/shp)
      z3 = 1. - z(jt) - xh + (pmq(1)**2+pmq(2)**2)/shp
      z(3-jt) = 2./(z1**2-z2**2)*(z1*z3+z2*sqrt(z3**2-(z1**2-z2**2)*pmq(3-jt)**2/shp))
      zmin = 2.*pmq(3-jt)/shpr
      zmax = 1. - pmq(jt)/shpr - (sh-pmq(3-jt)**2)/(shpr*(shpr-pmq(jt)))
      zmax = min(1.-xh, zmax)
      If (z(3-jt)<zmin .Or. z(3-jt)>zmax) Goto 340
      kcc = 22

    Else If (isub==78) Then
!...W+/- + H0 -> W+/- + H0.

    Else If (isub==79) Then
!...H0 + H0 -> H0 + H0.
    End If

  Else If (isub<=90) Then
    If (isub==81) Then
!...q + qb -> Q' + Qb'; th = (p(q)-p(q'))**2.
      mint(21) = isign(mint(46), mint(15))
      mint(22) = -mint(21)
      kcc = 4

    Else If (isub==82) Then
!...g + g -> Q + Qb; th arbitrary.
      kcs = (-1)**int(1.5+rlu(0))
      mint(21) = isign(mint(46), kcs)
      mint(22) = -mint(21)
      kcc = mint(2) + 10
    End If

  Else If (isub<=100) Then
    If (isub==95) Then
!...Low-pT ( = energyless g + g -> g + g).
      kcc = mint(2) + 12
      kcs = (-1)**int(1.5+rlu(0))

    Else If (isub==96) Then
!...Multiple interactions (should be reassigned to QCD process).
    End If

  Else If (isub<=110) Then
    If (isub==101) Then
!...g + g -> gamma*/Z0.
      kcc = 21
      kfres = 22

    Else If (isub==102) Then
!...g + g -> H0.
      kcc = 21
      kfres = 25
    End If

  Else If (isub<=120) Then
    If (isub==111) Then
!...f + fb -> g + H0; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(20+js) = 21
      mint(23-js) = 25
      kcc = 17 + js

    Else If (isub==112) Then
!...f + g -> f + H0; th = (p(f) - p(f))**2.
      If (mint(15)==21) js = 2
      mint(23-js) = 25
      kcc = 15 + js
      kcs = isign(1, mint(14+js))

    Else If (isub==113) Then
!...g + g -> g + H0; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(23-js) = 25
      kcc = 22 + js
      kcs = (-1)**int(1.5+rlu(0))

    Else If (isub==114) Then
!...g + g -> gamma + gamma; th arbitrary.
      If (rlu(0)>0.5) js = 2
      mint(21) = 22
      mint(22) = 22
      kcc = 21

    Else If (isub==115) Then
!...g + g -> gamma + Z0.

    Else If (isub==116) Then
!...g + g -> Z0 + Z0.

    Else If (isub==117) Then
!...g + g -> W+ + W-.
    End If

  Else If (isub<=140) Then
    If (isub==121) Then
!...g + g -> f + fb + H0.
    End If

  Else If (isub<=160) Then
    If (isub==141) Then
!...f + fb -> gamma*/Z0/Z'0.
      kfres = 32

    Else If (isub==142) Then
!...f + fb' -> H+/-.
      kch1 = kchg(iabs(mint(15)), 1)*isign(1, mint(15))
      kch2 = kchg(iabs(mint(16)), 1)*isign(1, mint(16))
      kfres = isign(37, kch1+kch2)

    Else If (isub==143) Then
!...f + fb' -> R.
      kfres = isign(40, mint(15)+mint(16))
    End If

  Else
    If (isub==161) Then
!...g + f -> H+/- + f'; th = (p(f)-p(f))**2.
      If (mint(16)==21) js = 2
      ia = iabs(mint(17-js))
      mint(20+js) = isign(37, kchg(ia,1)*mint(17-js))
      ja = ia + mod(ia, 2) - mod(ia+1, 2)
      mint(23-js) = isign(ja, mint(17-js))
      kcc = 18 - js
      If (mint(15)/=21) kcs = isign(1, mint(15))
      If (mint(16)/=21) kcs = isign(1, mint(16))
    End If
  End If

  If (idoc==7) Then
!...Resonance not decaying: store colour connection indices.
    i = mint(83) + 7
    k(ipu3, 1) = 1
    k(ipu3, 2) = kfres
    k(ipu3, 3) = i
    p(ipu3, 4) = shuser
    p(ipu3, 5) = shuser
    k(ipu1, 4) = ipu2
    k(ipu1, 5) = ipu2
    k(ipu2, 4) = ipu1
    k(ipu2, 5) = ipu1
    k(i, 1) = 21
    k(i, 2) = kfres
    p(i, 4) = shuser
    p(i, 5) = shuser
    n = ipu3
    mint(21) = kfres
    mint(22) = 0

  Else If (idoc==8) Then
!...2 -> 2 processes: store outgoing partons in their CM-frame.
    Do jt = 1, 2
      i = mint(84) + 2 + jt
      k(i, 1) = 1
      If (iabs(mint(20+jt))<=10 .Or. mint(20+jt)==21) k(i, 1) = 3
      k(i, 2) = mint(20+jt)
      k(i, 3) = mint(83) + idoc + jt - 2
      If (iabs(k(i,2))<=10 .Or. k(i,2)==21) Then
        p(i, 5) = ulmass(k(i,2))
      Else
        p(i, 5) = sqrt(vint(63+mod(js+jt,2)))
      End If
    End Do
    If (p(ipu3,5)+p(ipu4,5)>=shr) Then
      kfa1 = iabs(mint(21))
      kfa2 = iabs(mint(22))
      If ((kfa1>3 .And. kfa1/=21) .Or. (kfa2>3 .And. kfa2/=21)) Then
        mint(51) = 1
        Return
      End If
      p(ipu3, 5) = 0.
      p(ipu4, 5) = 0.
    End If
    p(ipu3, 4) = 0.5*(shr+(p(ipu3,5)**2-p(ipu4,5)**2)/shr)
    p(ipu3, 3) = sqrt(max(0.,p(ipu3,4)**2-p(ipu3,5)**2))
    p(ipu4, 4) = shr - p(ipu3, 4)
    p(ipu4, 3) = -p(ipu3, 3)
    n = ipu4
    mint(7) = mint(83) + 7
    mint(8) = mint(83) + 8

!...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4).
    Call ludbrb(ipu3, ipu4, acos(vint(23)), vint(24), 0D0, 0D0, 0D0)

  Else If (idoc==9) Then
!'''2 -> 3 processes:

  Else If (idoc==11) Then
!...Z0 + Z0 -> H0, W+ + W- -> H0: store Higgs and outgoing partons.
    phi(1) = paru(2)*rlu(0)
    phi(2) = phi(1) - phir
    Do jt = 1, 2
      i = mint(84) + 2 + jt
      k(i, 1) = 1
      If (iabs(mint(20+jt))<=10 .Or. mint(20+jt)==21) k(i, 1) = 3
      k(i, 2) = mint(20+jt)
      k(i, 3) = mint(83) + idoc + jt - 2
      p(i, 5) = ulmass(k(i,2))
      If (0.5*shpr*z(jt)<=p(i,5)) p(i, 5) = 0.
      pabs = sqrt(max(0.,(0.5*shpr*z(jt))**2-p(i,5)**2))
      ptabs = pabs*sqrt(max(0.,1.-cthe(jt)**2))
      p(i, 1) = ptabs*cos(phi(jt))
      p(i, 2) = ptabs*sin(phi(jt))
      p(i, 3) = pabs*cthe(jt)*(-1)**(jt+1)
      p(i, 4) = 0.5*shpr*z(jt)
      izw = mint(83) + 6 + jt
      k(izw, 1) = 21
      k(izw, 2) = 23
      If (isub==8) k(izw, 2) = isign(24, luchge(mint(14+jt)))
      k(izw, 3) = izw - 2
      p(izw, 1) = -p(i, 1)
      p(izw, 2) = -p(i, 2)
      p(izw, 3) = (0.5*shpr-pabs*cthe(jt))*(-1)**(jt+1)
      p(izw, 4) = 0.5*shpr*(1.-z(jt))
      p(izw, 5) = -sqrt(max(0.,p(izw,3)**2+ptabs**2-p(izw,4)**2))
    End Do
    i = mint(83) + 9
    k(ipu5, 1) = 1
    k(ipu5, 2) = kfres
    k(ipu5, 3) = i
    p(ipu5, 5) = shr
    p(ipu5, 1) = -p(ipu3, 1) - p(ipu4, 1)
    p(ipu5, 2) = -p(ipu3, 2) - p(ipu4, 2)
    p(ipu5, 3) = -p(ipu3, 3) - p(ipu4, 3)
    p(ipu5, 4) = shpr - p(ipu3, 4) - p(ipu4, 4)
    k(i, 1) = 21
    k(i, 2) = kfres
    Do j = 1, 5
      p(i, j) = p(ipu5, j)
    End Do
    n = ipu5
    mint(23) = kfres

  Else If (idoc==12) Then
!...Z0 and W+/- scattering: store bosons and outgoing partons.
    phi(1) = paru(2)*rlu(0)
    phi(2) = phi(1) - phir
    Do jt = 1, 2
      i = mint(84) + 2 + jt
      k(i, 1) = 1
      If (iabs(mint(20+jt))<=10 .Or. mint(20+jt)==21) k(i, 1) = 3
      k(i, 2) = mint(20+jt)
      k(i, 3) = mint(83) + idoc + jt - 2
      p(i, 5) = ulmass(k(i,2))
      If (0.5*shpr*z(jt)<=p(i,5)) p(i, 5) = 0.
      pabs = sqrt(max(0.,(0.5*shpr*z(jt))**2-p(i,5)**2))
      ptabs = pabs*sqrt(max(0.,1.-cthe(jt)**2))
      p(i, 1) = ptabs*cos(phi(jt))
      p(i, 2) = ptabs*sin(phi(jt))
      p(i, 3) = pabs*cthe(jt)*(-1)**(jt+1)
      p(i, 4) = 0.5*shpr*z(jt)
      izw = mint(83) + 6 + jt
      k(izw, 1) = 21
      If (mint(14+jt)==mint(20+jt)) Then
        k(izw, 2) = 23
      Else
        k(izw, 2) = isign(24, luchge(mint(14+jt))-luchge(mint(20+jt)))
      End If
      k(izw, 3) = izw - 2
      p(izw, 1) = -p(i, 1)
      p(izw, 2) = -p(i, 2)
      p(izw, 3) = (0.5*shpr-pabs*cthe(jt))*(-1)**(jt+1)
      p(izw, 4) = 0.5*shpr*(1.-z(jt))
      p(izw, 5) = -sqrt(max(0.,p(izw,3)**2+ptabs**2-p(izw,4)**2))
      ipu = mint(84) + 4 + jt
      k(ipu, 1) = 3
      k(ipu, 2) = kfpr(isub, jt)
      k(ipu, 3) = mint(83) + 8 + jt
      If (iabs(k(ipu,2))<=10 .Or. k(ipu,2)==21) Then
        p(ipu, 5) = ulmass(k(ipu,2))
      Else
        p(ipu, 5) = sqrt(vint(63+mod(js+jt,2)))
      End If
      mint(22+jt) = k(izw, 2)
    End Do
    If (isub==72) k(mint(84)+4+int(1.5+rlu(0)), 2) = -24
!...Find rotation and boost for hard scattering subsystem.
    i1 = mint(83) + 7
    i2 = mint(83) + 8
    bexcm = (p(i1,1)+p(i2,1))/(p(i1,4)+p(i2,4))
    beycm = (p(i1,2)+p(i2,2))/(p(i1,4)+p(i2,4))
    bezcm = (p(i1,3)+p(i2,3))/(p(i1,4)+p(i2,4))
    gamcm = (p(i1,4)+p(i2,4))/shr
    bepcm = bexcm*p(i1, 1) + beycm*p(i1, 2) + bezcm*p(i1, 3)
    px = p(i1, 1) + gamcm*(gamcm/(1.+gamcm)*bepcm-p(i1,4))*bexcm
    py = p(i1, 2) + gamcm*(gamcm/(1.+gamcm)*bepcm-p(i1,4))*beycm
    pz = p(i1, 3) + gamcm*(gamcm/(1.+gamcm)*bepcm-p(i1,4))*bezcm
    thecm = ulangl(pz, sqrt(px**2+py**2))
    phicm = ulangl(px, py)
!...Store hard scattering subsystem. Rotate and boost it.
    sqlam = (sh-p(ipu5,5)**2-p(ipu6,5)**2)**2 - 4.*p(ipu5, 5)**2*p(ipu6, 5)**2
    pabs = sqrt(max(0.,sqlam/(4.*sh)))
    cthwz = vint(23)
    sthwz = sqrt(max(0.,1.-cthwz**2))
    phiwz = vint(24) - phicm
    p(ipu5, 1) = pabs*sthwz*cos(phiwz)
    p(ipu5, 2) = pabs*sthwz*sin(phiwz)
    p(ipu5, 3) = pabs*cthwz
    p(ipu5, 4) = sqrt(pabs**2+p(ipu5,5)**2)
    p(ipu6, 1) = -p(ipu5, 1)
    p(ipu6, 2) = -p(ipu5, 2)
    p(ipu6, 3) = -p(ipu5, 3)
    p(ipu6, 4) = sqrt(pabs**2+p(ipu6,5)**2)
    Call ludbrb(ipu5, ipu6, thecm, phicm, dble(bexcm), dble(beycm), dble(bezcm))
    Do jt = 1, 2
      i1 = mint(83) + 8 + jt
      i2 = mint(84) + 4 + jt
      k(i1, 1) = 21
      k(i1, 2) = k(i2, 2)
      Do j = 1, 5
        p(i1, j) = p(i2, j)
      End Do
    End Do
    n = ipu6
    mint(7) = mint(83) + 9
    mint(8) = mint(83) + 10
  End If

  If (idoc>=8) Then
!...Store colour connection indices.
    Do j = 1, 2
      jc = j
      If (kcs==-1) jc = 3 - j
      If (icol(kcc,1,jc)/=0 .And. k(ipu1,1)==14) k(ipu1, j+3) = k(ipu1, j+3) + mint(84) + icol(kcc, 1, jc)
      If (icol(kcc,2,jc)/=0 .And. k(ipu2,1)==14) k(ipu2, j+3) = k(ipu2, j+3) + mint(84) + icol(kcc, 2, jc)
      If (icol(kcc,3,jc)/=0 .And. k(ipu3,1)==3) k(ipu3, j+3) = mstu(5)*(mint(84)+icol(kcc,3,jc))
      If (icol(kcc,4,jc)/=0 .And. k(ipu4,1)==3) k(ipu4, j+3) = mstu(5)*(mint(84)+icol(kcc,4,jc))
    End Do

!...Copy outgoing partons to documentation lines.
    Do i = 1, 2
      i1 = mint(83) + idoc - 2 + i
      i2 = mint(84) + 2 + i
      k(i1, 1) = 21
      k(i1, 2) = k(i2, 2)
      If (idoc<=9) k(i1, 3) = 0
      If (idoc>=11) k(i1, 3) = mint(83) + 2 + i
      Do j = 1, 5
        p(i1, j) = p(i2, j)
      End Do
    End Do
  End If
  mint(52) = n

!...Low-pT events: remove gluons used for string drawing purposes.
  If (isub==95) Then
    k(ipu3, 1) = k(ipu3, 1) + 10
    k(ipu4, 1) = k(ipu4, 1) + 10
    Do j = 41, 66
      vint(j) = 0.
    End Do
    Do i = mint(83) + 5, mint(83) + 8
      Do j = 1, 5
        p(i, j) = 0.
      End Do
    End Do
  End If

  Return
End Subroutine pyscat

!*********************************************************************

!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine pysspa(ipu1, ipu2)

!...Generates spacelike parton showers.
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Dimension kfls(4), is(2), xs(2), zs(2), q2s(2), tevs(2), robo(5), xfs(2, -6:6), xfa(-6:6), xfb(-6:6), xfn(-6:6), wtap(-6:6), wtsf(-6:6), the2(2), alam(2), dq2(3), dpc(3), dpd(4), dpb(4)

!...Calculate maximum virtuality and check that evolution possible.
  ipus1 = ipu1
  ipus2 = ipu2
  isub = mint(1)
  q2e = vint(52)
  If (iset(isub)==1) Then
    q2e = q2e/parp(67)
  Else If (iset(isub)==3 .Or. iset(isub)==4) Then
    q2e = pmas(23, 1)**2
    If (isub==8 .Or. isub==76 .Or. isub==77) q2e = pmas(24, 1)**2
  End If
  tmax = log(parp(67)*parp(63)*q2e/parp(61)**2)
  If (parp(67)*q2e<max(parp(62)**2,2.*parp(61)**2) .Or. tmax<0.2) Return

!...Common constants and initial values. Save normal Lambda value.
  xe0 = 2.*parp(65)/vint(1)
  alams = paru(111)
  paru(111) = parp(61)
  ns = n
  100 n = ns
  Do jt = 1, 2
    kfls(jt) = mint(14+jt)
    kfls(jt+2) = kfls(jt)
    xs(jt) = vint(40+jt)
    zs(jt) = 1.
    q2s(jt) = parp(67)*q2e
    tevs(jt) = tmax
    alam(jt) = parp(61)
    the2(jt) = 100.
    Do kfl = -6, 6
      xfs(jt, kfl) = xsfx(jt, kfl)
    End Do
  End Do
  dsh = dble(vint(44))
  If (iset(isub)==3 .Or. iset(isub)==4) dsh = dble(vint(26)*vint(2))

!...Pick up leg with highest virtuality.
  120 n = n + 1
  jt = 1
  If (n>ns+1 .And. q2s(2)>q2s(1)) jt = 2
  kflb = kfls(jt)
  xb = xs(jt)
  Do kfl = -6, 6
    xfb(kfl) = xfs(jt, kfl)
  End Do
  dshr = 2D0*sqrt(dsh)
  dshz = dsh/dble(zs(jt))
  xe = max(xe0, xb*(1./(1.-parp(66))-1.))
  If (xb+xe>=0.999) Then
    q2b = 0.
    Goto 220
  End If

!...Maximum Q2 without or with Q2 ordering. Effective Lambda and n_f.
  If (mstp(62)<=1) Then
    q2b = 0.5*(1./zs(jt)+1.)*q2s(jt) + 0.5*(1./zs(jt)-1.)*(q2s(3-jt)-sngl(dsh)+sqrt((sngl(dsh)+q2s(1)+q2s(2))**2+8.*q2s(1)*q2s(2)*zs(jt)/(1.-zs(jt))))
    tevb = log(parp(63)*q2b/alam(jt)**2)
  Else
    q2b = q2s(jt)
    tevb = tevs(jt)
  End If
  alsdum = ulalps(parp(63)*q2b)
  tevb = tevb + 2.*log(alam(jt)/paru(117))
  tevbsv = tevb
  alam(jt) = paru(117)
  b0 = (33.-2.*mstu(118))/6.

!...Calculate Altarelli-Parisi and structure function weights.
  Do kfl = -6, 6
    wtap(kfl) = 0.
    wtsf(kfl) = 0.
  End Do
  If (kflb==21) Then
    wtapq = 16.*(1.-sqrt(xb+xe))/(3.*sqrt(xb))
    Do kfl = -mstp(54), mstp(54)
      If (kfl==0) wtap(kfl) = 6.*log((1.-xb)/xe)
      If (kfl/=0) wtap(kfl) = wtapq
    End Do
  Else
    wtap(0) = 0.5*xb*(1./(xb+xe)-1.)
    wtap(kflb) = 8.*log((1.-xb)*(xb+xe)/xe)/3.
  End If
  160 wtsum = 0.
  If (kflb/=21) xfbo = xfb(kflb)
  If (kflb==21) xfbo = xfb(0)
!***************************************************************
!**********ERROR HAS OCCURED HERE
  If (xfbo==0.0) Then
    Write (mstu(11), 1000)
    Write (mstu(11), 1001) kflb, xfb(kflb)
    xfbo = 0.00001
  End If
!****************************************************************
  Do kfl = -mstp(54), mstp(54)
    wtsf(kfl) = xfb(kfl)/xfbo
    wtsum = wtsum + wtap(kfl)*wtsf(kfl)
  End Do
  wtsum = max(0.0001, wtsum)

!...Choose new t: fix alpha_s, alpha_s(Q2), alpha_s(k_T2).
  180 If (mstp(64)<=0) Then
    tevb = tevb + log(rlu(0))*paru(2)/(paru(111)*wtsum)
  Else If (mstp(64)==1) Then
    tevb = tevb*exp(max(-100.,log(rlu(0))*b0/wtsum))
  Else
    tevb = tevb*exp(max(-100.,log(rlu(0))*b0/(5.*wtsum)))
  End If
  190 q2ref = alam(jt)**2*exp(tevb)
  q2b = q2ref/parp(63)

!...Evolution ended or select flavour for branching parton.
  If (q2b<parp(62)**2) Then
    q2b = 0.
  Else
    wtran = rlu(0)*wtsum
    kfla = -mstp(54) - 1
    200 kfla = kfla + 1
    wtran = wtran - wtap(kfla)*wtsf(kfla)
    If (kfla<mstp(54) .And. wtran>0.) Goto 200
    If (kfla==0) kfla = 21

!...Choose z value and corrective weight.
    If (kflb==21 .And. kfla==21) Then
      z = 1./(1.+((1.-xb)/xb)*(xe/(1.-xb))**rlu(0))
      wtz = (1.-z*(1.-z))**2
    Else If (kflb==21) Then
      z = xb/(1.-rlu(0)*(1.-sqrt(xb+xe)))**2
      wtz = 0.5*(1.+(1.-z)**2)*sqrt(z)
    Else If (kfla==21) Then
      z = xb*(1.+rlu(0)*(1./(xb+xe)-1.))
      wtz = 1. - 2.*z*(1.-z)
    Else
      z = 1. - (1.-xb)*(xe/((xb+xe)*(1.-xb)))**rlu(0)
      wtz = 0.5*(1.+z**2)
    End If

!...Option with resummation of soft gluon emission as effective z shift.
    If (mstp(65)>=1) Then
      rsoft = 6.
      If (kflb/=21) rsoft = 8./3.
      z = z*(tevb/tevs(jt))**(rsoft*xe/((xb+xe)*b0))
      If (z<=xb) Goto 180
    End If

!...Option with alpha_s(k_T2)Q2): demand k_T2 > cutoff, reweight.
    If (mstp(64)>=2) Then
      If ((1.-z)*q2b<parp(62)**2) Goto 180
      alprat = tevb/(tevb+log(1.-z))
      If (alprat<5.*rlu(0)) Goto 180
      If (alprat>5.) wtz = wtz*alprat/5.
    End If

!...Option with angular ordering requirement.
    If (mstp(62)>=3) Then
      the2t = (4.*z**2*q2b)/(vint(2)*(1.-z)*xb**2)
      If (the2t>the2(jt)) Goto 180
    End If

!...Weighting with new structure functions.
    Call pystfu(mint(10+jt), xb, q2ref, xfn, jt)
    If (kflb/=21) xfbn = xfn(kflb)
    If (kflb==21) xfbn = xfn(0)
    If (xfbn<1E-20) Then
      If (kfla==kflb) Then
        tevb = tevbsv
        wtap(kflb) = 0.
        Goto 160
      Else If (tevbsv-tevb>0.2) Then
        tevb = 0.5*(tevbsv+tevb)
        Goto 190
      Else
        xfbn = 1E-10
      End If
    End If
    Do kfl = -mstp(54), mstp(54)
      xfb(kfl) = xfn(kfl)
    End Do
    xa = xb/z
    Call pystfu(mint(10+jt), xa, q2ref, xfa, jt)
    If (kfla/=21) xfan = xfa(kfla)
    If (kfla==21) xfan = xfa(0)
    If (xfan<1E-20) Goto 160
    If (kfla/=21) wtsfa = wtsf(kfla)
    If (kfla==21) wtsfa = wtsf(0)
    If (wtz*xfan/xfbn<rlu(0)*wtsfa) Goto 160
  End If

!...Define two hard scatterers in their CM-frame.
  220 If (n==ns+2) Then
    dq2(jt) = dble(q2b)
    dplcm = dsqrt((dsh+dq2(1)+dq2(2))**2-4D0*dq2(1)*dq2(2))/dshr
    Do jr = 1, 2
      i = ns + jr
      If (jr==1) ipo = ipus1
      If (jr==2) ipo = ipus2
      Do j = 1, 5
        k(i, j) = 0
        p(i, j) = 0.
        v(i, j) = 0.
      End Do
      k(i, 1) = 14
      k(i, 2) = kfls(jr+2)
      k(i, 4) = ipo
      k(i, 5) = ipo
      p(i, 3) = sngl(dplcm)*(-1)**(jr+1)
      p(i, 4) = sngl((dsh+dq2(3-jr)-dq2(jr))/dshr)
      p(i, 5) = -sqrt(sngl(dq2(jr)))
      k(ipo, 1) = 14
      k(ipo, 3) = i
      k(ipo, 4) = mod(k(ipo,4), mstu(5)) + mstu(5)*i
      k(ipo, 5) = mod(k(ipo,5), mstu(5)) + mstu(5)*i
    End Do

!...Find maximum allowed mass of timelike parton.
  Else If (n>ns+2) Then
    jr = 3 - jt
    dq2(3) = dble(q2b)
    dpc(1) = dble(p(is(1),4))
    dpc(2) = dble(p(is(2),4))
    dpc(3) = dble(0.5*(abs(p(is(1),3))+abs(p(is(2),3))))
    dpd(1) = dsh + dq2(jr) + dq2(jt)
    dpd(2) = dshz + dq2(jr) + dq2(3)
    dpd(3) = sqrt(dpd(1)**2-4D0*dq2(jr)*dq2(jt))
    dpd(4) = sqrt(dpd(2)**2-4D0*dq2(jr)*dq2(3))
    ikin = 0
    If (q2s(jr)>=(0.5*parp(62))**2 .And. dpd(1)-dpd(3)>=1D-10*dpd(1)) ikin = 1
    If (ikin==0) dmsma = (dq2(jt)/dble(zs(jt))-dq2(3))*(dsh/(dsh+dq2(jt))-dsh/(dshz+dq2(3)))
    If (ikin==1) dmsma = (dpd(1)*dpd(2)-dpd(3)*dpd(4))/(2.D0*dq2(jr)) - dq2(jt) - dq2(3)

!...Generate timelike parton shower (if required).
    it = n
    Do j = 1, 5
      k(it, j) = 0
      p(it, j) = 0.
      v(it, j) = 0.
    End Do
    k(it, 1) = 3
    k(it, 2) = 21
    If (kflb==21 .And. kfls(jt+2)/=21) k(it, 2) = -kfls(jt+2)
    If (kflb/=21 .And. kfls(jt+2)==21) k(it, 2) = kflb
    p(it, 5) = ulmass(k(it,2))
    If (sngl(dmsma)<=p(it,5)**2) Goto 100
    If (mstp(63)>=1) Then
      p(it, 4) = sngl((dshz-dsh-dble(p(it,5))**2)/dshr)
      p(it, 3) = sqrt(p(it,4)**2-p(it,5)**2)
      If (mstp(63)==1) Then
        q2tim = sngl(dmsma)
      Else If (mstp(63)==2) Then
        q2tim = min(sngl(dmsma), parp(71)*q2s(jt))
      Else
!'''Here remains to introduce angular ordering in first branching.
        q2tim = sngl(dmsma)
      End If
      Call lushow(it, 0, sqrt(q2tim))
      If (n>=it+1) p(it, 5) = p(it+1, 5)
    End If

!...Reconstruct kinematics of branching: timelike parton shower.
    dms = dble(p(it,5)**2)
    If (ikin==0) dpt2 = (dmsma-dms)*(dshz+dq2(3))/(dsh+dq2(jt))
    If (ikin==1) dpt2 = (dmsma-dms)*(0.5D0*dpd(1)*dpd(2)+0.5D0*dpd(3)*dpd(4)-dq2(jr)*(dq2(jt)+dq2(3)+dms))/(4.D0*dsh*dpc(3)**2)
    If (dpt2<0.D0) Goto 100
    dpb(1) = (0.5D0*dpd(2)-dpc(jr)*(dshz+dq2(jr)-dq2(jt)-dms)/dshr)/dpc(3) - dpc(3)
    p(it, 1) = sqrt(sngl(dpt2))
    p(it, 3) = sngl(dpb(1))*(-1)**(jt+1)
    p(it, 4) = sngl((dshz-dsh-dms)/dshr)
    If (n>=it+1) Then
      dpb(1) = sqrt(dpb(1)**2+dpt2)
      dpb(2) = sqrt(dpb(1)**2+dms)
      dpb(3) = dble(p(it+1,3))
      dpb(4) = sqrt(dpb(3)**2+dms)
      dbez = (dpb(4)*dpb(1)-dpb(3)*dpb(2))/(dpb(4)*dpb(2)-dpb(3)*dpb(1))
      Call ludbrb(it+1, n, 0., 0., 0D0, 0D0, dbez)
      the = ulangl(p(it,3), p(it,1))
      Call ludbrb(it+1, n, the, 0., 0D0, 0D0, 0D0)
    End If

!...Reconstruct kinematics of branching: spacelike parton.
    Do j = 1, 5
      k(n+1, j) = 0
      p(n+1, j) = 0.
      v(n+1, j) = 0.
    End Do
    k(n+1, 1) = 14
    k(n+1, 2) = kflb
    p(n+1, 1) = p(it, 1)
    p(n+1, 3) = p(it, 3) + p(is(jt), 3)
    p(n+1, 4) = p(it, 4) + p(is(jt), 4)
    p(n+1, 5) = -sqrt(sngl(dq2(3)))

!...Define colour flow of branching.
    k(is(jt), 3) = n + 1
    k(it, 3) = n + 1
    id1 = it
    If ((k(n+1,2)>0 .And. k(n+1,2)/=21 .And. k(id1,2)>0 .And. k(id1,2)/=21) .Or. (k(n+1,2)<0 .And. k(id1,2)==21) .Or. (k(n+1,2)==21 .And. k(id1,2)==21 .And. rlu(0)>0.5) .Or. (k(n+1,2)==21 .And. k(id1,2)<0)) id1 = is(jt)
    id2 = it + is(jt) - id1
    k(n+1, 4) = k(n+1, 4) + id1
    k(n+1, 5) = k(n+1, 5) + id2
    k(id1, 4) = k(id1, 4) + mstu(5)*(n+1)
    k(id1, 5) = k(id1, 5) + mstu(5)*id2
    k(id2, 4) = k(id2, 4) + mstu(5)*id1
    k(id2, 5) = k(id2, 5) + mstu(5)*(n+1)
    n = n + 1

!...Boost to new CM-frame.
    Call ludbrb(ns+1, n, 0., 0., -dble((p(n,1)+p(is(jr),1))/(p(n,4)+p(is(jr),4))), 0D0, -dble((p(n,3)+p(is(jr),3))/(p(n,4)+p(is(jr),4))))
    ir = n + (jt-1)*(is(1)-n)
    Call ludbrb(ns+1, n, -ulangl(p(ir,3),p(ir,1)), paru(2)*rlu(0), 0D0, 0D0, 0D0)
  End If

!...Save quantities, loop back.
  is(jt) = n
  q2s(jt) = q2b
  dq2(jt) = dble(q2b)
  If (mstp(62)>=3) the2(jt) = the2t
  dsh = dshz
  If (q2b>=(0.5*parp(62))**2) Then
    kfls(jt+2) = kfls(jt)
    kfls(jt) = kfla
    xs(jt) = xa
    zs(jt) = z
    Do kfl = -6, 6
      xfs(jt, kfl) = xfa(kfl)
    End Do
    tevs(jt) = tevb
  Else
    If (jt==1) ipu1 = n
    If (jt==2) ipu2 = n
  End If
  If (n>mstu(4)-mstu(32)-10) Then
    Call luerrm(11, '(PYSSPA:) no more memory left in LUJETS')
    If (mstu(21)>=1) n = ns
    If (mstu(21)>=1) Return
  End If
  If (max(q2s(1),q2s(2))>=(0.5*parp(62))**2 .Or. n<=ns+1) Goto 120

!...Boost hard scattering partons to frame of shower initiators.
  Do j = 1, 3
    robo(j+2) = (p(ns+1,j)+p(ns+2,j))/(p(ns+1,4)+p(ns+2,4))
  End Do
  Do j = 1, 5
    p(n+2, j) = p(ns+1, j)
  End Do
  robot = robo(3)**2 + robo(4)**2 + robo(5)**2
  If (robot>=0.999999) Then
    robot = 1.00001*sqrt(robot)
    robo(3) = robo(3)/robot
    robo(4) = robo(4)/robot
    robo(5) = robo(5)/robot
  End If
  Call ludbrb(n+2, n+2, 0., 0., -dble(robo(3)), -dble(robo(4)), -dble(robo(5)))
  robo(2) = ulangl(p(n+2,1), p(n+2,2))
  robo(1) = ulangl(p(n+2,3), sqrt(p(n+2,1)**2+p(n+2,2)**2))
  Call ludbrb(mint(83)+5, ns, robo(1), robo(2), dble(robo(3)), dble(robo(4)), dble(robo(5)))

!...Store user information. Reset Lambda value.
  k(ipu1, 3) = mint(83) + 3
  k(ipu2, 3) = mint(83) + 4
  Do jt = 1, 2
    mint(12+jt) = kfls(jt)
    vint(140+jt) = xs(jt)
  End Do
  paru(111) = alams

  Return
  1000 Format (5X, 'structure function has a zero point here')
  1001 Format (5X, 'xf(x,i=', I5, ')=', F10.5)
End Subroutine pysspa

!*********************************************************************

Subroutine pymult(mmul)

!...Initializes treatment of multiple interactions, selects kinematics
!...of hardest interaction if low-pT physics included in run, and
!...generates all non-hardest interactions.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension nmul(20), sigm(20), kstr(500, 2)
  Save xt2, xt2fac, xc2, xts, irbin, rbin, nmul, sigm

!...Initialization of multiple interaction treatment.
  If (mmul==1) Then
    If (mstp(122)>=1) Write (mstu(11), 1000) mstp(82)
    isub = 96
    mint(1) = 96
    vint(63) = 0.
    vint(64) = 0.
    vint(143) = 1.
    vint(144) = 1.

!...Loop over phase space points: xT2 choice in 20 bins.
    100 sigsum = 0.
    Do ixt2 = 1, 20
      nmul(ixt2) = mstp(83)
      sigm(ixt2) = 0.
      Do itry = 1, mstp(83)
        rsca = 0.05*((21-ixt2)-rlu(0))
        xt2 = vint(149)*(1.+vint(149))/(vint(149)+rsca) - vint(149)
        xt2 = max(0.01*vint(149), xt2)
        vint(25) = xt2

!...Choose tau and y*. Calculate cos(theta-hat).
        If (rlu(0)<=coef(isub,1)) Then
          taup = (2.*(1.+sqrt(1.-xt2))/xt2-1.)**rlu(0)
          tau = xt2*(1.+taup)**2/(4.*taup)
        Else
          tau = xt2*(1.+tan(rlu(0)*atan(sqrt(1./xt2-1.)))**2)
        End If
        vint(21) = tau
        Call pyklim(2)
        ryst = rlu(0)
        myst = 1
        If (ryst>coef(isub,7)) myst = 2
        If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
        Call pykmap(2, myst, rlu(0))
        vint(23) = sqrt(max(0.,1.-xt2/tau))*(-1)**int(1.5+rlu(0))

!...Calculate differential cross-section.
        vint(71) = 0.5*vint(1)*sqrt(xt2)
        Call pysigh(nchn, sigs)
        sigm(ixt2) = sigm(ixt2) + sigs
      End Do
      sigsum = sigsum + sigm(ixt2)
    End Do
    sigsum = sigsum/(20.*mstp(83))

!...Reject result if sigma(parton-parton) is smaller than hadronic one.
    If (sigsum<1.1*vint(106)) Then
      If (mstp(122)>=1) Write (mstu(11), 1100) parp(82), sigsum
      parp(82) = 0.9*parp(82)
      vint(149) = 4.*parp(82)**2/vint(2)
      Goto 100
    End If
    If (mstp(122)>=1) Write (mstu(11), 1200) parp(82), sigsum

!...Start iteration to find k factor.
    yke = sigsum/vint(106)
    so = 0.5
    xi = 0.
    yi = 0.
    xk = 0.5
    iit = 0
    130 If (iit==0) Then
      xk = 2.*xk
    Else If (iit==1) Then
      xk = 0.5*xk
    Else
      xk = xi + (yke-yi)*(xf-xi)/(yf-yi)
    End If

!...Evaluate overlap integrals.
    If (mstp(82)==2) Then
      sp = 0.5*paru(1)*(1.-exp(-xk))
      sop = sp/paru(1)
    Else
      If (mstp(82)==3) deltab = 0.02
      If (mstp(82)==4) deltab = min(0.01, 0.05*parp(84))
      sp = 0.
      sop = 0.
      b = -0.5*deltab
      140 b = b + deltab
      If (mstp(82)==3) Then
        ov = exp(-b**2)/paru(2)
      Else
        cq2 = parp(84)**2
        ov = ((1.-parp(83))**2*exp(-min(100.,b**2))+2.*parp(83)*(1.-parp(83))*2./(1.+cq2)*exp(-min(100.,b**2*2./(1.+cq2)))+parp(83)**2/cq2*exp(-min(100.,b**2/cq2)))/paru(2)
      End If
      pacc = 1. - exp(-min(100.,paru(1)*xk*ov))
      sp = sp + paru(2)*b*deltab*pacc
      sop = sop + paru(2)*b*deltab*ov*pacc
      If (b<1. .Or. b*pacc>1E-6) Goto 140
    End If
    yk = paru(1)*xk*so/sp

!...Continue iteration until convergence.
    If (yk<yke) Then
      xi = xk
      yi = yk
      If (iit==1) iit = 2
    Else
      xf = xk
      yf = yk
      If (iit==0) iit = 1
    End If
    If (abs(yk-yke)>=1E-5*yke) Goto 130

!...Store some results for subsequent use.
    vint(145) = sigsum
    vint(146) = sop/so
    vint(147) = sop/sp

!...Initialize iteration in xT2 for hardest interaction.
  Else If (mmul==2) Then
    If (mstp(82)<=0) Then
    Else If (mstp(82)==1) Then
      xt2 = 1.
      xt2fac = xsec(96, 1)/vint(106)*vint(149)/(1.-vint(149))
    Else If (mstp(82)==2) Then
      xt2 = 1.
      xt2fac = vint(146)*xsec(96, 1)/vint(106)*vint(149)*(1.+vint(149))
    Else
      xc2 = 4.*ckin(3)**2/vint(2)
      If (ckin(3)<=ckin(5) .Or. mint(82)>=2) xc2 = 0.
    End If

  Else If (mmul==3) Then
!...Low-pT or multiple interactions (first semihard interaction):
!...choose xT2 according to dpT2/pT2**2*exp(-(sigma above pT2)/norm)
!...or (MSTP(82)>=2) dpT2/(pT2+pT0**2)**2*exp(-....).
    isub = mint(1)
    If (mstp(82)<=0) Then
      xt2 = 0.
    Else If (mstp(82)==1) Then
      xt2 = xt2fac*xt2/(xt2fac-xt2*log(rlu(0)))
    Else If (mstp(82)==2) Then
      If (xt2<1. .And. exp(-xt2fac*xt2/(vint(149)*(xt2+vint(149))))>rlu(0)) xt2 = 1.
      If (xt2>=1.) Then
        xt2 = (1.+vint(149))*xt2fac/(xt2fac-(1.+vint(149))*log(1.-rlu(0)*(1.-exp(-xt2fac/(vint(149)*(1.+vint(149))))))) - vint(149)
      Else
        xt2 = -xt2fac/log(exp(-xt2fac/(xt2+vint(149)))+rlu(0)*(exp(-xt2fac/vint(149))-exp(-xt2fac/(xt2+vint(149))))) - vint(149)
      End If
      xt2 = max(0.01*vint(149), xt2)
    Else
      xt2 = (xc2+vint(149))*(1.+vint(149))/(1.+vint(149)-rlu(0)*(1.-xc2)) - vint(149)
      xt2 = max(0.01*vint(149), xt2)
    End If
    vint(25) = xt2

!...Low-pT: choose xT2, tau, y* and cos(theta-hat) fixed.
    If (mstp(82)<=1 .And. xt2<vint(149)) Then
      If (mint(82)==1) ngen(0, 1) = ngen(0, 1) - 1
      If (mint(82)==1) ngen(isub, 1) = ngen(isub, 1) - 1
      isub = 95
      mint(1) = isub
      vint(21) = 0.01*vint(149)
      vint(22) = 0.
      vint(23) = 0.
      vint(25) = 0.01*vint(149)

    Else
!...Multiple interactions (first semihard interaction).
!...Choose tau and y*. Calculate cos(theta-hat).
      If (rlu(0)<=coef(isub,1)) Then
        taup = (2.*(1.+sqrt(1.-xt2))/xt2-1.)**rlu(0)
        tau = xt2*(1.+taup)**2/(4.*taup)
      Else
        tau = xt2*(1.+tan(rlu(0)*atan(sqrt(1./xt2-1.)))**2)
      End If
      vint(21) = tau
      Call pyklim(2)
      ryst = rlu(0)
      myst = 1
      If (ryst>coef(isub,7)) myst = 2
      If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
      Call pykmap(2, myst, rlu(0))
      vint(23) = sqrt(max(0.,1.-xt2/tau))*(-1)**int(1.5+rlu(0))
    End If
    vint(71) = 0.5*vint(1)*sqrt(vint(25))

!...Store results of cross-section calculation.
  Else If (mmul==4) Then
    isub = mint(1)
    xts = vint(25)
    If (iset(isub)==1) xts = vint(21)
    If (iset(isub)==2) xts = (4.*vint(48)+2.*vint(63)+2.*vint(64))/vint(2)
    If (iset(isub)==3 .Or. iset(isub)==4) xts = vint(26)
    rbin = max(0.000001, min(0.999999,xts*(1.+vint(149))/(xts+vint(149))))
    irbin = int(1.+20.*rbin)
    If (isub==96) nmul(irbin) = nmul(irbin) + 1
    If (isub==96) sigm(irbin) = sigm(irbin) + vint(153)

!...Choose impact parameter.
  Else If (mmul==5) Then
    If (mstp(82)==3) Then
      vint(148) = rlu(0)/(paru(2)*vint(147))
    Else
      rtype = rlu(0)
      cq2 = parp(84)**2
      If (rtype<(1.-parp(83))**2) Then
        b2 = -log(rlu(0))
      Else If (rtype<1.-parp(83)**2) Then
        b2 = -0.5*(1.+cq2)*log(rlu(0))
      Else
        b2 = -cq2*log(rlu(0))
      End If
      vint(148) = ((1.-parp(83))**2*exp(-min(100.,b2))+2.*parp(83)*(1.-parp(83))*2./(1.+cq2)*exp(-min(100.,b2*2./(1.+cq2)))+parp(83)**2/cq2*exp(-min(100.,b2/cq2)))/(paru(2)*vint(147))
    End If

!...Multiple interactions (variable impact parameter) : reject with
!...probability exp(-overlap*cross-section above pT/normalization).
    rncor = (irbin-20.*rbin)*nmul(irbin)
    sigcor = (irbin-20.*rbin)*sigm(irbin)
    Do ibin = irbin + 1, 20
      rncor = rncor + nmul(ibin)
      sigcor = sigcor + sigm(ibin)
    End Do
    sigabv = (sigcor/rncor)*vint(149)*(1.-xts)/(xts+vint(149))
    vint(150) = exp(-min(100.,vint(146)*vint(148)*sigabv/vint(106)))

!...Generate additional multiple semihard interactions.
  Else If (mmul==6) Then

!...Reconstruct strings in hard scattering.
    isub = mint(1)
    nmax = mint(84) + 4
    If (iset(isub)==1) nmax = mint(84) + 2
    nstr = 0
    Do i = mint(84) + 1, nmax
      kcs = kchg(lucomp(k(i,2)), 2)*isign(1, k(i,2))
      If (kcs==0) Goto 170
      Do j = 1, 4
        If (kcs==1 .And. (j==2 .Or. j==4)) Goto 160
        If (kcs==-1 .And. (j==1 .Or. j==3)) Goto 160
        If (j<=2) Then
          ist = mod(k(i,j+3)/mstu(5), mstu(5))
        Else
          ist = mod(k(i,j+1), mstu(5))
        End If
        If (ist<mint(84) .Or. ist>i) Goto 160
        If (kchg(lucomp(k(ist,2)),2)==0) Goto 160
        nstr = nstr + 1
        If (j==1 .Or. j==4) Then
          kstr(nstr, 1) = i
          kstr(nstr, 2) = ist
        Else
          kstr(nstr, 1) = ist
          kstr(nstr, 2) = i
        End If
      160 End Do
    170 End Do

!...Set up starting values for iteration in xT2.
    xt2 = vint(25)
    If (iset(isub)==1) xt2 = vint(21)
    If (iset(isub)==2) xt2 = (4.*vint(48)+2.*vint(63)+2.*vint(64))/vint(2)
    If (iset(isub)==3 .Or. iset(isub)==4) xt2 = vint(26)
    isub = 96
    mint(1) = 96
    If (mstp(82)<=1) Then
      xt2fac = xsec(isub, 1)*vint(149)/((1.-vint(149))*vint(106))
    Else
      xt2fac = vint(146)*vint(148)*xsec(isub, 1)/vint(106)*vint(149)*(1.+vint(149))
    End If
    vint(63) = 0.
    vint(64) = 0.
    vint(151) = 0.
    vint(152) = 0.
    vint(143) = 1. - vint(141)
    vint(144) = 1. - vint(142)

!...Iterate downwards in xT2.
    180 If (mstp(82)<=1) Then
      xt2 = xt2fac*xt2/(xt2fac-xt2*log(rlu(0)))
      If (xt2<vint(149)) Goto 220
    Else
      If (xt2<=0.01*vint(149)) Goto 220
      xt2 = xt2fac*(xt2+vint(149))/(xt2fac-(xt2+vint(149))*log(rlu(0))) - vint(149)
      If (xt2<=0.) Goto 220
      xt2 = max(0.01*vint(149), xt2)
    End If
    vint(25) = xt2

!...Choose tau and y*. Calculate cos(theta-hat).
    If (rlu(0)<=coef(isub,1)) Then
      taup = (2.*(1.+sqrt(1.-xt2))/xt2-1.)**rlu(0)
      tau = xt2*(1.+taup)**2/(4.*taup)
    Else
      tau = xt2*(1.+tan(rlu(0)*atan(sqrt(1./xt2-1.)))**2)
    End If
    vint(21) = tau
    Call pyklim(2)
    ryst = rlu(0)
    myst = 1
    If (ryst>coef(isub,7)) myst = 2
    If (ryst>coef(isub,7)+coef(isub,8)) myst = 3
    Call pykmap(2, myst, rlu(0))
    vint(23) = sqrt(max(0.,1.-xt2/tau))*(-1)**int(1.5+rlu(0))

!...Check that x not used up. Accept or reject kinematical variables.
    x1m = sqrt(tau)*exp(vint(22))
    x2m = sqrt(tau)*exp(-vint(22))
    If (vint(143)-x1m<0.01 .Or. vint(144)-x2m<0.01) Goto 180
    vint(71) = 0.5*vint(1)*sqrt(xt2)
    Call pysigh(nchn, sigs)
    If (sigs<xsec(isub,1)*rlu(0)) Goto 180

!...Reset K, P and V vectors. Select some variables.
    Do i = n + 1, n + 2
      Do j = 1, 5
        k(i, j) = 0
        p(i, j) = 0.
        v(i, j) = 0.
      End Do
    End Do
    rflav = rlu(0)
    pt = 0.5*vint(1)*sqrt(xt2)
    phi = paru(2)*rlu(0)
    cth = vint(23)

!...Add first parton to event record.
    k(n+1, 1) = 3
    k(n+1, 2) = 21
    If (rflav>=max(parp(85),parp(86))) k(n+1, 2) = 1 + int((2.+parj(2))*rlu(0))
    p(n+1, 1) = pt*cos(phi)
    p(n+1, 2) = pt*sin(phi)
    p(n+1, 3) = 0.25*vint(1)*(vint(41)*(1.+cth)-vint(42)*(1.-cth))
    p(n+1, 4) = 0.25*vint(1)*(vint(41)*(1.+cth)+vint(42)*(1.-cth))
    p(n+1, 5) = 0.

!...Add second parton to event record.
    k(n+2, 1) = 3
    k(n+2, 2) = 21
    If (k(n+1,2)/=21) k(n+2, 2) = -k(n+1, 2)
    p(n+2, 1) = -p(n+1, 1)
    p(n+2, 2) = -p(n+1, 2)
    p(n+2, 3) = 0.25*vint(1)*(vint(41)*(1.-cth)-vint(42)*(1.+cth))
    p(n+2, 4) = 0.25*vint(1)*(vint(41)*(1.-cth)+vint(42)*(1.+cth))
    p(n+2, 5) = 0.

    If (rflav<parp(85) .And. nstr>=1) Then
!....Choose relevant string pieces to place gluons on.
      Do i = n + 1, n + 2
        dmin = 1E8
        Do istr = 1, nstr
          i1 = kstr(istr, 1)
          i2 = kstr(istr, 2)
          dist = (p(i,4)*p(i1,4)-p(i,1)*p(i1,1)-p(i,2)*p(i1,2)-p(i,3)*p(i1,3))*(p(i,4)*p(i2,4)-p(i,1)*p(i2,1)-p(i,2)*p(i2,2)-p(i,3)*p(i2,3))/max(1., p(i1,4)*p(i2,4)-p(i1,1)*p(i2,1)-p(i1,2)*p(i2,2)-p(i1,3)*p(i2,3))
          If (istr==1 .Or. dist<dmin) Then
            dmin = dist
            ist1 = i1
            ist2 = i2
            istm = istr
          End If
        End Do

!....Colour flow adjustments, new string pieces.
        If (k(ist1,4)/mstu(5)==ist2) k(ist1, 4) = mstu(5)*i + mod(k(ist1,4), mstu(5))
        If (mod(k(ist1,5),mstu(5))==ist2) k(ist1, 5) = mstu(5)*(k(ist1,5)/mstu(5)) + i
        k(i, 5) = mstu(5)*ist1
        k(i, 4) = mstu(5)*ist2
        If (k(ist2,5)/mstu(5)==ist1) k(ist2, 5) = mstu(5)*i + mod(k(ist2,5), mstu(5))
        If (mod(k(ist2,4),mstu(5))==ist1) k(ist2, 4) = mstu(5)*(k(ist2,4)/mstu(5)) + i
        kstr(istm, 2) = i
        kstr(nstr+1, 1) = i
        kstr(nstr+1, 2) = ist2
        nstr = nstr + 1
      End Do

!...String drawing and colour flow for gluon loop.
    Else If (k(n+1,2)==21) Then
      k(n+1, 4) = mstu(5)*(n+2)
      k(n+1, 5) = mstu(5)*(n+2)
      k(n+2, 4) = mstu(5)*(n+1)
      k(n+2, 5) = mstu(5)*(n+1)
      kstr(nstr+1, 1) = n + 1
      kstr(nstr+1, 2) = n + 2
      kstr(nstr+2, 1) = n + 2
      kstr(nstr+2, 2) = n + 1
      nstr = nstr + 2

!...String drawing and colour flow for q-qbar pair.
    Else
      k(n+1, 4) = mstu(5)*(n+2)
      k(n+2, 5) = mstu(5)*(n+1)
      kstr(nstr+1, 1) = n + 1
      kstr(nstr+1, 2) = n + 2
      nstr = nstr + 1
    End If

!...Update remaining energy; iterate.
    n = n + 2
    If (n>mstu(4)-mstu(32)-10) Then
      Call luerrm(11, '(PYMULT:) no more memory left in LUJETS')
      If (mstu(21)>=1) Return
    End If
    mint(31) = mint(31) + 1
    vint(151) = vint(151) + vint(41)
    vint(152) = vint(152) + vint(42)
    vint(143) = vint(143) - vint(41)
    vint(144) = vint(144) - vint(42)
    If (mint(31)<240) Goto 180
    220 Continue
  End If

  Return

!...Format statements for printout.
  1000 Format (/1X, '****** PYMULT: initialization of multiple inter', 'actions for MSTP(82) =', I2, ' ******')
  1100 Format (8X, 'pT0 =', F5.2, ' GeV gives sigma(parton-parton) =', 1P, E9.2, ' mb: rejected')
  1200 Format (8X, 'pT0 =', F5.2, ' GeV gives sigma(parton-parton) =', 1P, E9.2, ' mb: accepted')
End Subroutine pymult

!*********************************************************************

Subroutine pyremn(ipu1, ipu2)

!...Adds on target remnants (one or two from each side) and
!...includes primordial kT.
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save /hparnt/
  Common /hstrng/nfp(300, 15), pphi(300, 15), nft(300, 15), pthi(300, 15)
  Save /hstrng/
!...COMMON BLOCK FROM HIJING
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Dimension kflch(2), kflsp(2), chi(2), pms(6), is(2), robo(5)

!...Special case for lepton-lepton interaction.
  If (mint(43)==1) Then
    Do jt = 1, 2
      i = mint(83) + jt + 2
      k(i, 1) = 21
      k(i, 2) = k(i-2, 2)
      k(i, 3) = i - 2
      Do j = 1, 5
        p(i, j) = p(i-2, j)
      End Do
    End Do
  End If

!...Find event type, set pointers.
  If (ipu1==0 .And. ipu2==0) Return
  isub = mint(1)
  ilep = 0
  If (ipu1==0) ilep = 1
  If (ipu2==0) ilep = 2
  If (isub==95) ilep = -1
  If (ilep==1) iq = mint(84) + 1
  If (ilep==2) iq = mint(84) + 2
  ip = max(ipu1, ipu2)
  ilepr = mint(83) + 5 - ilep
  ns = n

!...Define initial partons, including primordial kT.
  110 Do jt = 1, 2
    i = mint(83) + jt + 2
    If (jt==1) ipu = ipu1
    If (jt==2) ipu = ipu2
    k(i, 1) = 21
    k(i, 3) = i - 2
    If (isub==95) Then
      k(i, 2) = 21
      shs = 0.
    Else If (mint(40+jt)==1 .And. ipu/=0) Then
      k(i, 2) = k(ipu, 2)
      p(i, 5) = p(ipu, 5)
      p(i, 1) = 0.
      p(i, 2) = 0.
      pms(jt) = p(i, 5)**2
    Else If (ipu/=0) Then
      k(i, 2) = k(ipu, 2)
      p(i, 5) = p(ipu, 5)
!...No primordial kT or chosen according to truncated Gaussian or
!...exponential.
!
!     X.N. Wang (7.22.97)
!
      rpt1 = 0.0
      rpt2 = 0.0
      ssw2 = (pphi(ihnt2(11),4)+pthi(ihnt2(12),4))**2 - (pphi(ihnt2(11),1)+pthi(ihnt2(12),1))**2 - (pphi(ihnt2(11),2)+pthi(ihnt2(12),2))**2 - (pphi(ihnt2(11),3)+pthi(ihnt2(12),3))**2
!
!********this is s of the current NN collision
      If (ssw2<=4.0*parp(93)**2) Goto 1211
!
      If (ihpr2(5)<=0) Then
        120 If (mstp(91)<=0) Then
          pt = 0.
        Else If (mstp(91)==1) Then
          pt = parp(91)*sqrt(-log(rlu(0)))
        Else
          rpt1 = rlu(0)
          rpt2 = rlu(0)
          pt = -parp(92)*log(rpt1*rpt2)
        End If
        If (pt>parp(93)) Goto 120
        phi = paru(2)*rlu(0)
        rpt1 = pt*cos(phi)
        rpt2 = pt*sin(phi)
      Else If (ihpr2(5)==1) Then
        If (jt==1) jpt = nfp(ihnt2(11), 11)
        If (jt==2) jpt = nft(ihnt2(12), 11)
        1205 ptgs = parp(91)*sqrt(-log(rlu(0)))
        If (ptgs>parp(93)) Goto 1205
        phi = 2.0*hipr1(40)*rlu(0)
        rpt1 = ptgs*cos(phi)
        rpt2 = ptgs*sin(phi)
        Do iint = 1, jpt - 1
          pkcsq = parp(91)*sqrt(-log(rlu(0)))
          phi = 2.0*hipr1(40)*rlu(0)
          rpt1 = rpt1 + pkcsq*cos(phi)
          rpt2 = rpt2 + pkcsq*sin(phi)
        End Do
        If (rpt1**2+rpt2**2>=ssw2/4.0) Goto 1205
      End If
!     X.N. Wang
!                     ********When initial interaction among soft partons is
!                             assumed the primordial pt comes from the sum of
!                             pt of JPT-1 number of initial interaction, JPT
!                             is the number of interaction including present
!                             one that nucleon hassuffered
      1211 p(i, 1) = rpt1
      p(i, 2) = rpt2
      pms(jt) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
    Else
      k(i, 2) = k(iq, 2)
      q2 = vint(52)
      p(i, 5) = -sqrt(q2)
      pms(jt) = -q2
      shs = (1.-vint(43-jt))*q2/vint(43-jt) + vint(5-jt)**2
    End If
  End Do

!...Kinematics construction for initial partons.
  i1 = mint(83) + 3
  i2 = mint(83) + 4
  If (ilep==0) shs = vint(141)*vint(142)*vint(2) + (p(i1,1)+p(i2,1))**2 + (p(i1,2)+p(i2,2))**2
  shr = sqrt(max(0.,shs))
  If (ilep==0) Then
    If ((shs-pms(1)-pms(2))**2-4.*pms(1)*pms(2)<=0.) Goto 110
    p(i1, 4) = 0.5*(shr+(pms(1)-pms(2))/shr)
    p(i1, 3) = sqrt(max(0.,p(i1,4)**2-pms(1)))
    p(i2, 4) = shr - p(i1, 4)
    p(i2, 3) = -p(i1, 3)
  Else If (ilep==1) Then
    p(i1, 4) = p(iq, 4)
    p(i1, 3) = p(iq, 3)
    p(i2, 4) = p(ip, 4)
    p(i2, 3) = p(ip, 3)
  Else If (ilep==2) Then
    p(i1, 4) = p(ip, 4)
    p(i1, 3) = p(ip, 3)
    p(i2, 4) = p(iq, 4)
    p(i2, 3) = p(iq, 3)
  End If
  If (mint(43)==1) Return

!...Transform partons to overall CM-frame (not for leptoproduction).
  If (ilep==0) Then
    robo(3) = (p(i1,1)+p(i2,1))/shr
    robo(4) = (p(i1,2)+p(i2,2))/shr
    Call ludbrb(i1, i2, 0., 0., -dble(robo(3)), -dble(robo(4)), 0D0)
    robo(2) = ulangl(p(i1,1), p(i1,2))
    Call ludbrb(i1, i2, 0., -robo(2), 0D0, 0D0, 0D0)
    robo(1) = ulangl(p(i1,3), p(i1,1))
    Call ludbrb(i1, i2, -robo(1), 0., 0D0, 0D0, 0D0)
    nmax = max(mint(52), ipu1, ipu2)
    Call ludbrb(i1, nmax, robo(1), robo(2), dble(robo(3)), dble(robo(4)), 0D0)
    robo(5) = max(-0.999999, min(0.999999,(vint(141)-vint(142))/(vint(141)+vint(142))))
    Call ludbrb(i1, nmax, 0., 0., 0D0, 0D0, dble(robo(5)))
  End If

!...Check invariant mass of remnant system:
!...hadronic events or leptoproduction.
  If (ilep<=0) Then
    If (mstp(81)<=0 .Or. mstp(82)<=0 .Or. isub==95) Then
      vint(151) = 0.
      vint(152) = 0.
    End If
    peh = p(i1, 4) + p(i2, 4) + 0.5*vint(1)*(vint(151)+vint(152))
    pzh = p(i1, 3) + p(i2, 3) + 0.5*vint(1)*(vint(151)-vint(152))
    shh = (vint(1)-peh)**2 - (p(i1,1)+p(i2,1))**2 - (p(i1,2)+p(i2,2))**2 - pzh**2
    pmmin = p(mint(83)+1, 5) + p(mint(83)+2, 5) + ulmass(k(i1,2)) + ulmass(k(i2,2))
    If (shr>=vint(1) .Or. shh<=(pmmin+parp(111))**2) Then
      mint(51) = 1
      Return
    End If
    shr = sqrt(shh+(p(i1,1)+p(i2,1))**2+(p(i1,2)+p(i2,2))**2)
  Else
    pei = p(iq, 4) + p(ip, 4)
    pzi = p(iq, 3) + p(ip, 3)
    pms(ilep) = max(0., pei**2-pzi**2)
    pmmin = p(ilepr-2, 5) + ulmass(k(ilepr,2)) + sqrt(pms(ilep))
    If (shr<=pmmin+parp(111)) Then
      mint(51) = 1
      Return
    End If
  End If

!...Subdivide remnant if necessary, store first parton.
  140 i = ns
  Do jt = 1, 2
    If (jt==ilep) Goto 190
    If (jt==1) ipu = ipu1
    If (jt==2) ipu = ipu2
    Call pyspli(mint(10+jt), mint(12+jt), kflch(jt), kflsp(jt))
    i = i + 1
    is(jt) = i
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
    k(i, 1) = 3
    k(i, 2) = kflsp(jt)
    k(i, 3) = mint(83) + jt
    p(i, 5) = ulmass(k(i,2))

!...First parton colour connections and transverse mass.
    kfls = (3-kchg(lucomp(kflsp(jt)),2)*isign(1,kflsp(jt)))/2
    k(i, kfls+3) = ipu
    k(ipu, 6-kfls) = mod(k(ipu,6-kfls), mstu(5)) + mstu(5)*i
    If (kflch(jt)==0) Then
      p(i, 1) = -p(mint(83)+jt+2, 1)
      p(i, 2) = -p(mint(83)+jt+2, 2)
      pms(jt) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2

!...When extra remnant parton or hadron: find relative pT, store.
    Else
      Call luptdi(1, p(i,1), p(i,2))
      pms(jt+2) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
      i = i + 1
      Do j = 1, 5
        k(i, j) = 0
        p(i, j) = 0.
        v(i, j) = 0.
      End Do
      k(i, 1) = 1
      k(i, 2) = kflch(jt)
      k(i, 3) = mint(83) + jt
      p(i, 5) = ulmass(k(i,2))
      p(i, 1) = -p(mint(83)+jt+2, 1) - p(i-1, 1)
      p(i, 2) = -p(mint(83)+jt+2, 2) - p(i-1, 2)
      pms(jt+4) = p(i, 5)**2 + p(i, 1)**2 + p(i, 2)**2
!...Relative distribution of energy for particle into two jets.
      imb = 1
      If (mod(mint(10+jt)/1000,10)/=0) imb = 2
      If (iabs(kflch(jt))<=10 .Or. kflch(jt)==21) Then
        chik = parp(92+2*imb)
        If (mstp(92)<=1) Then
          If (imb==1) chi(jt) = rlu(0)
          If (imb==2) chi(jt) = 1. - sqrt(rlu(0))
        Else If (mstp(92)==2) Then
          chi(jt) = 1. - rlu(0)**(1./(1.+chik))
        Else If (mstp(92)==3) Then
          cut = 2.*0.3/vint(1)
          170 chi(jt) = rlu(0)**2
          If ((chi(jt)**2/(chi(jt)**2+cut**2))**0.25*(1.-chi(jt))**chik<rlu(0)) Goto 170
        Else
          cut = 2.*0.3/vint(1)
          cutr = (1.+sqrt(1.+cut**2))/cut
          180 chir = cut*cutr**rlu(0)
          chi(jt) = (chir**2-cut**2)/(2.*chir)
          If ((1.-chi(jt))**chik<rlu(0)) Goto 180
        End If
!...Relative distribution of energy for particle into jet plus particle.
      Else
        If (mstp(92)<=1) Then
          If (imb==1) chi(jt) = rlu(0)
          If (imb==2) chi(jt) = 1. - sqrt(rlu(0))
        Else
          chi(jt) = 1. - rlu(0)**(1./(1.+parp(93+2*imb)))
        End If
        If (mod(kflch(jt)/1000,10)/=0) chi(jt) = 1. - chi(jt)
      End If
      pms(jt) = pms(jt+4)/chi(jt) + pms(jt+2)/(1.-chi(jt))
      kfls = kchg(lucomp(kflch(jt)), 2)*isign(1, kflch(jt))
      If (kfls/=0) Then
        k(i, 1) = 3
        kfls = (3-kfls)/2
        k(i, kfls+3) = ipu
        k(ipu, 6-kfls) = mod(k(ipu,6-kfls), mstu(5)) + mstu(5)*i
      End If
    End If
  190 End Do
  If (shr<=sqrt(pms(1))+sqrt(pms(2))) Goto 140
  n = i

!...Reconstruct kinematics of remnants.
  Do jt = 1, 2
    If (jt==ilep) Goto 200
    pe = 0.5*(shr+(pms(jt)-pms(3-jt))/shr)
    pz = sqrt(pe**2-pms(jt))
    If (kflch(jt)==0) Then
      p(is(jt), 4) = pe
      p(is(jt), 3) = pz*(-1)**(jt-1)
    Else
      pw1 = chi(jt)*(pe+pz)
      p(is(jt)+1, 4) = 0.5*(pw1+pms(jt+4)/pw1)
      p(is(jt)+1, 3) = 0.5*(pw1-pms(jt+4)/pw1)*(-1)**(jt-1)
      p(is(jt), 4) = pe - p(is(jt)+1, 4)
      p(is(jt), 3) = pz*(-1)**(jt-1) - p(is(jt)+1, 3)
    End If
  200 End Do

!...Hadronic events: boost remnants to correct longitudinal frame.
  If (ilep<=0) Then
    Call ludbrb(ns+1, n, 0., 0., 0D0, 0D0, -dble(pzh/(vint(1)-peh)))
!...Leptoproduction events: boost colliding subsystem.
  Else
    nmax = max(ip, mint(52))
    pef = shr - pe
    pzf = pz*(-1)**(ilep-1)
    pt2 = p(ilepr, 1)**2 + p(ilepr, 2)**2
    phipt = ulangl(p(ilepr,1), p(ilepr,2))
    Call ludbrb(mint(84)+1, nmax, 0., -phipt, 0D0, 0D0, 0D0)
    rqp = p(iq, 3)*(pt2+pei**2) - p(iq, 4)*pei*pzi
    sinth = p(iq, 4)*sqrt(pt2*(pt2+pei**2)/(rqp**2+pt2*p(iq,4)**2*pzi**2))*sign(1., -rqp)
    Call ludbrb(mint(84)+1, nmax, asin(sinth), 0., 0D0, 0D0, 0D0)
    betax = (-pei*pzi*sinth+sqrt(pt2*(pt2+pei**2-(pzi*sinth)**2)))/(pt2+pei**2)
    Call ludbrb(mint(84)+1, nmax, 0., 0., dble(betax), 0D0, 0D0)
    Call ludbrb(mint(84)+1, nmax, 0., phipt, 0D0, 0D0, 0D0)
    pem = p(iq, 4) + p(ip, 4)
    pzm = p(iq, 3) + p(ip, 3)
    betaz = (-pem*pzm+pzf*sqrt(pzf**2+pem**2-pzm**2))/(pzf**2+pem**2)
    Call ludbrb(mint(84)+1, nmax, 0., 0., 0D0, 0D0, dble(betaz))
    Call ludbrb(i1, i2, asin(sinth), 0., dble(betax), 0D0, 0D0)
    Call ludbrb(i1, i2, 0., phipt, 0D0, 0D0, dble(betaz))
  End If

  Return
End Subroutine pyremn

!*********************************************************************

Subroutine pyresd

!...Allows resonances to decay (including parton showers for hadronic
!...channels).
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Dimension iref(10, 6), kdcy(2), kfl1(2), kfl2(2), nsd(2), ilin(6), coup(6, 4), pk(6, 4), pkk(6, 6), cthe(2), phi(2), wdtp(0:40), wdte(0:40, 0:5)
  Complex fgk, ha(6, 6), hc(6, 6)

!...The F, Xi and Xj functions of Gunion and Kunszt
!...(Phys. Rev. D33, 665, plus errata from the authors).
  fgk(i1, i2, i3, i4, i5, i6) = 4.*ha(i1, i3)*hc(i2, i6)*(ha(i1,i5)*hc(i1,i4)+ha(i3,i5)*hc(i3,i4))
  digk(dt, du) = -4.D0*d34*d56 + dt*(3.D0*dt+4.D0*du) + dt**2*(dt*du/(d34*d56)-2.D0*(1.D0/d34+1.D0/d56)*(dt+du)+2.D0*(d34/d56+d56/d34))
  djgk(dt, du) = 8.D0*(d34+d56)**2 - 8.D0*(d34+d56)*(dt+du) - 6.D0*dt*du - 2.D0*dt*du*(dt*du/(d34*d56)-2.D0*(1.D0/d34+1.D0/d56)*(dt+du)+2.D0*(d34/d56+d56/d34))

!...Define initial two objects, initialize loop.
  isub = mint(1)
  sh = vint(44)
  iref(1, 5) = 0
  iref(1, 6) = 0
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    iref(1, 1) = mint(84) + 2 + iset(isub)
    iref(1, 2) = 0
    iref(1, 3) = mint(83) + 6 + iset(isub)
    iref(1, 4) = 0
  Else If (iset(isub)==2 .Or. iset(isub)==4) Then
    iref(1, 1) = mint(84) + 1 + iset(isub)
    iref(1, 2) = mint(84) + 2 + iset(isub)
    iref(1, 3) = mint(83) + 5 + iset(isub)
    iref(1, 4) = mint(83) + 6 + iset(isub)
  End If
  np = 1
  ip = 0
  100 ip = ip + 1
  ninh = 0

!...Loop over one/two resonances; reset decay rates.
  jtmax = 2
  If (ip==1 .And. (iset(isub)==1 .Or. iset(isub)==3)) jtmax = 1
  Do jt = 1, jtmax
    kdcy(jt) = 0
    kfl1(jt) = 0
    kfl2(jt) = 0
    nsd(jt) = iref(ip, jt)
    id = iref(ip, jt)
    If (id==0) Goto 140
    kfa = iabs(k(id,2))
    If (kfa<23 .Or. kfa>40) Goto 140
    If (mdcy(kfa,1)/=0) Then
      If (isub==1 .Or. isub==141) mint(61) = 1
      Call pywidt(kfa, p(id,5), wdtp, wdte)
      If (kchg(kfa,3)==0) Then
        ipm = 2
      Else
        ipm = (5+isign(1,k(id,2)))/2
      End If
      If (jtmax==1 .Or. iabs(k(iref(ip,1),2))/=iabs(k(iref(ip,2),2))) Then
        i12 = 4
      Else
        If (jt==1) i12 = int(4.5+rlu(0))
        i12 = 9 - i12
      End If
      rkfl = (wdte(0,1)+wdte(0,ipm)+wdte(0,i12))*rlu(0)
      Do i = 1, mdcy(kfa, 3)
        idc = i + mdcy(kfa, 2) - 1
        kfl1(jt) = kfdp(idc, 1)*isign(1, k(id,2))
        kfl2(jt) = kfdp(idc, 2)*isign(1, k(id,2))
        rkfl = rkfl - (wdte(i,1)+wdte(i,ipm)+wdte(i,i12))
        If (rkfl<=0.) Goto 130
      End Do
      130 Continue
    End If

!...Summarize result on decay channel chosen.
    If ((kfa==23 .Or. kfa==24) .And. kfl1(jt)==0) ninh = ninh + 1
    If (kfl1(jt)==0) Goto 140
    kdcy(jt) = 2
    If (iabs(kfl1(jt))<=10 .Or. kfl1(jt)==21) kdcy(jt) = 1
    If ((iabs(kfl1(jt))>=23 .And. iabs(kfl1(jt))<=25) .Or. (iabs(kfl1(jt))==37)) kdcy(jt) = 3
    nsd(jt) = n

!...Fill decay products, prepared for parton showers for quarks.
!lin-8/19/02 avoid actual argument in common blocks of LU2ENT:
    pid5 = p(id, 5)
    If (kdcy(jt)==1) Then
!        CALL LU2ENT(-(N+1),KFL1(JT),KFL2(JT),P(ID,5))
      Call lu2ent(-(n+1), kfl1(jt), kfl2(jt), pid5)
    Else
!        CALL LU2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))
      Call lu2ent(n+1, kfl1(jt), kfl2(jt), pid5)
    End If

    If (jtmax==1) Then
      cthe(jt) = vint(13) + (vint(33)-vint(13)+vint(34)-vint(14))*rlu(0)
      If (cthe(jt)>vint(33)) cthe(jt) = cthe(jt) + vint(14) - vint(33)
      phi(jt) = vint(24)
    Else
      cthe(jt) = 2.*rlu(0) - 1.
      phi(jt) = paru(2)*rlu(0)
    End If
  140 End Do
  If (mint(3)==1 .And. ip==1) Then
    mint(25) = kfl1(1)
    mint(26) = kfl2(1)
  End If
  If (jtmax==1 .And. kdcy(1)==0) Goto 530
  If (jtmax==2 .And. kdcy(1)==0 .And. kdcy(2)==0) Goto 530
  If (mstp(45)<=0 .Or. iref(ip,2)==0 .Or. ninh>=1) Goto 500
  If (k(iref(1,1),2)==25 .And. ip==1) Goto 500
  If (k(iref(1,1),2)==25 .And. kdcy(1)*kdcy(2)==0) Goto 500

!...Order incoming partons and outgoing resonances.
  ilin(1) = mint(84) + 1
  If (k(mint(84)+1,2)>0) ilin(1) = mint(84) + 2
  If (k(ilin(1),2)==21) ilin(1) = 2*mint(84) + 3 - ilin(1)
  ilin(2) = 2*mint(84) + 3 - ilin(1)
  imin = 1
  If (iref(ip,5)==25) imin = 3
  imax = 2
  iord = 1
  If (k(iref(ip,1),2)==23) iord = 2
  If (k(iref(ip,1),2)==24 .And. k(iref(ip,2),2)==-24) iord = 2
  If (iabs(k(iref(ip,iord),2))==25) iord = 3 - iord
  If (kdcy(iord)==0) iord = 3 - iord

!...Order decay products of resonances.
  Do jt = iord, 3 - iord, 3 - 2*iord
    If (kdcy(jt)==0) Then
      ilin(imax+1) = nsd(jt)
      imax = imax + 1
    Else If (k(nsd(jt)+1,2)>0) Then
      ilin(imax+1) = n + 2*jt - 1
      ilin(imax+2) = n + 2*jt
      imax = imax + 2
      k(n+2*jt-1, 2) = k(nsd(jt)+1, 2)
      k(n+2*jt, 2) = k(nsd(jt)+2, 2)
    Else
      ilin(imax+1) = n + 2*jt
      ilin(imax+2) = n + 2*jt - 1
      imax = imax + 2
      k(n+2*jt-1, 2) = k(nsd(jt)+1, 2)
      k(n+2*jt, 2) = k(nsd(jt)+2, 2)
    End If
  End Do

!...Find charge, isospin, left- and righthanded couplings.
  xw = paru(102)
  Do i = imin, imax
    Do j = 1, 4
      coup(i, j) = 0.
    End Do
    kfa = iabs(k(ilin(i),2))
    If (kfa>20) Goto 410
    coup(i, 1) = luchge(kfa)/3.
    coup(i, 2) = (-1)**mod(kfa, 2)
    coup(i, 4) = -2.*coup(i, 1)*xw
    coup(i, 3) = coup(i, 2) + coup(i, 4)
  410 End Do
  sqmz = pmas(23, 1)**2
  gzmz = pmas(23, 1)*pmas(23, 2)
  sqmw = pmas(24, 1)**2
  gzmw = pmas(24, 1)*pmas(24, 2)
  sqmzp = pmas(32, 1)**2
  gzmzp = pmas(32, 1)*pmas(32, 2)

!...Select random angles; construct massless four-vectors.
  420 Do i = n + 1, n + 4
    k(i, 1) = 1
    Do j = 1, 5
      p(i, j) = 0.
    End Do
  End Do
  Do jt = 1, jtmax
    If (kdcy(jt)==0) Goto 440
    id = iref(ip, jt)
    p(n+2*jt-1, 3) = 0.5*p(id, 5)
    p(n+2*jt-1, 4) = 0.5*p(id, 5)
    p(n+2*jt, 3) = -0.5*p(id, 5)
    p(n+2*jt, 4) = 0.5*p(id, 5)
    cthe(jt) = 2.*rlu(0) - 1.
    phi(jt) = paru(2)*rlu(0)
    Call ludbrb(n+2*jt-1, n+2*jt, acos(cthe(jt)), phi(jt), dble(p(id,1)/p(id,4)), dble(p(id,2)/p(id,4)), dble(p(id,3)/p(id,4)))
  440 End Do

!...Store incoming and outgoing momenta, with random rotation to
!...avoid accidental zeroes in HA expressions.
  Do i = 1, imax
    k(n+4+i, 1) = 1
    p(n+4+i, 4) = sqrt(p(ilin(i),1)**2+p(ilin(i),2)**2+p(ilin(i),3)**2+p(ilin(i),5)**2)
    p(n+4+i, 5) = p(ilin(i), 5)
    Do j = 1, 3
      p(n+4+i, j) = p(ilin(i), j)
    End Do
  End Do
  therr = acos(2.*rlu(0)-1.)
  phirr = paru(2)*rlu(0)
  Call ludbrb(n+5, n+4+imax, therr, phirr, 0D0, 0D0, 0D0)
  Do i = 1, imax
    Do j = 1, 4
      pk(i, j) = p(n+4+i, j)
    End Do
  End Do

!...Calculate internal products.
  If (isub==22 .Or. isub==23 .Or. isub==25) Then
    Do i1 = imin, imax - 1
      Do i2 = i1 + 1, imax
        ha(i1, i2) = sqrt((pk(i1,4)-pk(i1,3))*(pk(i2,4)+pk(i2,3))/(1E-20+pk(i1,1)**2+pk(i1,2)**2))*cmplx(pk(i1,1), pk(i1,2)) - sqrt((pk(i1,4)+pk(i1,3))*(pk(i2,4)-pk(i2,3))/(1E-20+pk(i2,1)**2+pk(i2,2)**2))*cmplx(pk(i2,1), pk(i2,2))
        hc(i1, i2) = conjg(ha(i1,i2))
        If (i1<=2) ha(i1, i2) = cmplx(0., 1.)*ha(i1, i2)
        If (i1<=2) hc(i1, i2) = cmplx(0., 1.)*hc(i1, i2)
        ha(i2, i1) = -ha(i1, i2)
        hc(i2, i1) = -hc(i1, i2)
      End Do
    End Do
  End If
  Do i = 1, 2
    Do j = 1, 4
      pk(i, j) = -pk(i, j)
    End Do
  End Do
  Do i1 = imin, imax - 1
    Do i2 = i1 + 1, imax
      pkk(i1, i2) = 2.*(pk(i1,4)*pk(i2,4)-pk(i1,1)*pk(i2,1)-pk(i1,2)*pk(i2,2)-pk(i1,3)*pk(i2,3))
      pkk(i2, i1) = pkk(i1, i2)
    End Do
  End Do

  If (iref(ip,5)==25) Then
!...Angular weight for H0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons
    wt = 16.*pkk(3, 5)*pkk(4, 6)
    If (ip==1) wtmax = sh**2
    If (ip>=2) wtmax = p(iref(ip,6), 5)**4

  Else If (isub==1) Then
    If (kfa/=37) Then
!...Angular weight for gamma*/Z0 -> 2 quarks/leptons
      ei = kchg(iabs(mint(15)), 1)/3.
      ai = sign(1., ei+0.1)
      vi = ai - 4.*ei*xw
      ef = kchg(kfa, 1)/3.
      af = sign(1., ef+0.1)
      vf = af - 4.*ef*xw
      gg = 1.
      gz = 1./(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gzmz**2)
      zz = 1./(16.*xw*(1.-xw))**2*sh**2/((sh-sqmz)**2+gzmz**2)
      If (mstp(43)==1) Then
!...Only gamma* production included
        gz = 0.
        zz = 0.
      Else If (mstp(43)==2) Then
!...Only Z0 production included
        gg = 0.
        gz = 0.
      End If
      asym = 2.*(ei*ai*gz*ef*af+4.*vi*ai*zz*vf*af)/(ei**2*gg*ef**2+ei*vi*gz*ef*vf+(vi**2+ai**2)*zz*(vf**2+af**2))
      wt = 1. + asym*cthe(jt) + cthe(jt)**2
      wtmax = 2. + abs(asym)
    Else
!...Angular weight for gamma*/Z0 -> H+ + H-
      wt = 1. - cthe(jt)**2
      wtmax = 1.
    End If

  Else If (isub==2) Then
!...Angular weight for W+/- -> 2 quarks/leptons
    wt = (1.+cthe(jt))**2
    wtmax = 4.

  Else If (isub==15 .Or. isub==19) Then
!...Angular weight for f + fb -> gluon/gamma + Z0 ->
!...-> gluon/gamma + 2 quarks/leptons
    wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*(pkk(1,3)**2+pkk(2,4)**2) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*(pkk(1,4)**2+pkk(2,3)**2)
    wtmax = (coup(1,3)**2+coup(1,4)**2)*(coup(3,3)**2+coup(3,4)**2)*((pkk(1,3)+pkk(1,4))**2+(pkk(2,3)+pkk(2,4))**2)

  Else If (isub==16 .Or. isub==20) Then
!...Angular weight for f + fb' -> gluon/gamma + W+/- ->
!...-> gluon/gamma + 2 quarks/leptons
    wt = pkk(1, 3)**2 + pkk(2, 4)**2
    wtmax = (pkk(1,3)+pkk(1,4))**2 + (pkk(2,3)+pkk(2,4))**2

  Else If (isub==22) Then
!...Angular weight for f + fb -> Z0 + Z0 -> 4 quarks/leptons
    s34 = p(iref(ip,iord), 5)**2
    s56 = p(iref(ip,3-iord), 5)**2
    ti = pkk(1, 3) + pkk(1, 4) + s34
    ui = pkk(1, 5) + pkk(1, 6) + s56
    wt = coup(1, 3)**4*((coup(3,3)*coup(5,3)*abs(fgk(1,2,3,4,5,6)/ti+fgk(1,2,5,6,3,4)/ui))**2+(coup(3,4)*coup(5,3)*abs(fgk(1,2,4,3,5,6)/ti+fgk(1,2,5,6,4,3)/ui))**2+(coup(3,3)*coup(5,4)*abs(fgk(1,2,3,4,6,5)/ti+fgk(1,2,6,5,3,4)/ui))**2+(coup(3,4)*coup(5,4)*abs(fgk(1,2,4,3,6,5)/ti+fgk(1,2,6,5,4,3)/ui))**2) + coup(1, 4)**4*((coup(3,3)*coup(5,3)*abs(fgk(2,1,5,6,3,4)/ti+fgk(2,1,3,4,5,6)/ui))**2+(coup(3,4)*coup(5,3)*abs(fgk(2,1,6,5,3,4)/ti+fgk(2,1,3,4,6,5)/ui))**2+(coup(3,3)*coup(5,4)*abs(fgk(2,1,5,6,4,3)/ti+fgk(2,1,4,3,5,6)/ui))**2+(coup(3,4)*coup(5,4)*abs(fgk(2,1,6,5,4,3)/ti+fgk(2,1,4,3,6,5)/ui))**2)
    wtmax = 4.*s34*s56*(coup(1,3)**4+coup(1,4)**4)*(coup(3,3)**2+coup(3,4)**2)*(coup(5,3)**2+coup(5,4)**2)*4.*(ti/ui+ui/ti+2.*sh*(s34+s56)/(ti*ui)-s34*s56*(1./ti**2+1./ui**2))

  Else If (isub==23) Then
!...Angular weight for f + fb' -> Z0 + W +/- -> 4 quarks/leptons
    d34 = dble(p(iref(ip,iord),5)**2)
    d56 = dble(p(iref(ip,3-iord),5)**2)
    dt = dble(pkk(1,3)+pkk(1,4)) + d34
    du = dble(pkk(1,5)+pkk(1,6)) + d56
    cawz = coup(2, 3)/sngl(dt) - 2.*(1.-xw)*coup(1, 2)/(sh-sqmw)
    cbwz = coup(1, 3)/sngl(du) + 2.*(1.-xw)*coup(1, 2)/(sh-sqmw)
    wt = coup(5, 3)**2*abs(cawz*fgk(1,2,3,4,5,6)+cbwz*fgk(1,2,5,6,3,4))**2 + coup(5, 4)**2*abs(cawz*fgk(1,2,3,4,6,5)+cbwz*fgk(1,2,6,5,3,4))**2
    wtmax = 4.*sngl(d34*d56)*(coup(5,3)**2+coup(5,4)**2)*(cawz**2*sngl(digk(dt,du))+cbwz**2*sngl(digk(du,dt))+cawz*cbwz*sngl(djgk(dt,du)))

  Else If (isub==24) Then
!...Angular weight for f + fb -> Z0 + H0 -> 2 quarks/leptons + H0
    wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*pkk(1, 3)*pkk(2, 4) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*pkk(1, 4)*pkk(2, 3)
    wtmax = (coup(1,3)**2+coup(1,4)**2)*(coup(3,3)**2+coup(3,4)**2)*(pkk(1,3)+pkk(1,4))*(pkk(2,3)+pkk(2,4))

  Else If (isub==25) Then
!...Angular weight for f + fb -> W+ + W- -> 4 quarks/leptons
    d34 = dble(p(iref(ip,iord),5)**2)
    d56 = dble(p(iref(ip,3-iord),5)**2)
    dt = dble(pkk(1,3)+pkk(1,4)) + d34
    du = dble(pkk(1,5)+pkk(1,6)) + d56
    cdww = (coup(1,3)*sqmz/(sh-sqmz)+coup(1,2))/sh
    caww = cdww + 0.5*(coup(1,2)+1.)/sngl(dt)
    cbww = cdww + 0.5*(coup(1,2)-1.)/sngl(du)
    ccww = coup(1, 4)*sqmz/(sh-sqmz)/sh
    wt = abs(caww*fgk(1,2,3,4,5,6)-cbww*fgk(1,2,5,6,3,4))**2 + ccww**2*abs(fgk(2,1,5,6,3,4)-fgk(2,1,3,4,5,6))**2
    wtmax = 4.*sngl(d34*d56)*(caww**2*sngl(digk(dt,du))+cbww**2*sngl(digk(du,dt))-caww*cbww*sngl(djgk(dt,du))+ccww**2*sngl(digk(dt,du)+digk(du,dt)-djgk(dt,du)))

  Else If (isub==26) Then
!...Angular weight for f + fb' -> W+/- + H0 -> 2 quarks/leptons + H0
    wt = pkk(1, 3)*pkk(2, 4)
    wtmax = (pkk(1,3)+pkk(1,4))*(pkk(2,3)+pkk(2,4))

  Else If (isub==30) Then
!...Angular weight for f + g -> f + Z0 -> f + 2 quarks/leptons
    If (k(ilin(1),2)>0) wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*(pkk(1,4)**2+pkk(3,5)**2) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*(pkk(1,3)**2+pkk(4,5)**2)
    If (k(ilin(1),2)<0) wt = ((coup(1,3)*coup(3,3))**2+(coup(1,4)*coup(3,4))**2)*(pkk(1,3)**2+pkk(4,5)**2) + ((coup(1,3)*coup(3,4))**2+(coup(1,4)*coup(3,3))**2)*(pkk(1,4)**2+pkk(3,5)**2)
    wtmax = (coup(1,3)**2+coup(1,4)**2)*(coup(3,3)**2+coup(3,4)**2)*((pkk(1,3)+pkk(1,4))**2+(pkk(3,5)+pkk(4,5))**2)

  Else If (isub==31) Then
!...Angular weight for f + g -> f' + W+/- -> f' + 2 quarks/leptons
    If (k(ilin(1),2)>0) wt = pkk(1, 4)**2 + pkk(3, 5)**2
    If (k(ilin(1),2)<0) wt = pkk(1, 3)**2 + pkk(4, 5)**2
    wtmax = (pkk(1,3)+pkk(1,4))**2 + (pkk(3,5)+pkk(4,5))**2

  Else If (isub==141) Then
!...Angular weight for gamma*/Z0/Z'0 -> 2 quarks/leptons
    ei = kchg(iabs(mint(15)), 1)/3.
    ai = sign(1., ei+0.1)
    vi = ai - 4.*ei*xw
    api = sign(1., ei+0.1)
    vpi = api - 4.*ei*xw
    ef = kchg(kfa, 1)/3.
    af = sign(1., ef+0.1)
    vf = af - 4.*ef*xw
    apf = sign(1., ef+0.1)
    vpf = apf - 4.*ef*xw
    gg = 1.
    gz = 1./(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gzmz**2)
    gzp = 1./(8.*xw*(1.-xw))*sh*(sh-sqmzp)/((sh-sqmzp)**2+gzmzp**2)
    zz = 1./(16.*xw*(1.-xw))**2*sh**2/((sh-sqmz)**2+gzmz**2)
    zzp = 2./(16.*xw*(1.-xw))**2*sh**2*((sh-sqmz)*(sh-sqmzp)+gzmz*gzmzp)/(((sh-sqmz)**2+gzmz**2)*((sh-sqmzp)**2+gzmzp**2))
    zpzp = 1./(16.*xw*(1.-xw))**2*sh**2/((sh-sqmzp)**2+gzmzp**2)
    If (mstp(44)==1) Then
!...Only gamma* production included
      gz = 0.
      gzp = 0.
      zz = 0.
      zzp = 0.
      zpzp = 0.
    Else If (mstp(44)==2) Then
!...Only Z0 production included
      gg = 0.
      gz = 0.
      gzp = 0.
      zzp = 0.
      zpzp = 0.
    Else If (mstp(44)==3) Then
!...Only Z'0 production included
      gg = 0.
      gz = 0.
      gzp = 0.
      zz = 0.
      zzp = 0.
    Else If (mstp(44)==4) Then
!...Only gamma*/Z0 production included
      gzp = 0.
      zzp = 0.
      zpzp = 0.
    Else If (mstp(44)==5) Then
!...Only gamma*/Z'0 production included
      gz = 0.
      zz = 0.
      zzp = 0.
    Else If (mstp(44)==6) Then
!...Only Z0/Z'0 production included
      gg = 0.
      gz = 0.
      gzp = 0.
    End If
    asym = 2.*(ei*ai*gz*ef*af+ei*api*gzp*ef*apf+4.*vi*ai*zz*vf*af+(vi*api+vpi*ai)*zzp*(vf*apf+vpf*af)+4.*vpi*api*zpzp*vpf*apf)/(ei**2*gg*ef**2+ei*vi*gz*ef*vf+ei*vpi*gzp*ef*vpf+(vi**2+ai**2)*zz*(vf**2+af**2)+(vi*vpi+ai*api)*zzp*(vf*vpf+af*apf)+(vpi**2+api**2)*zpzp*(vpf**2+apf**2))
    wt = 1. + asym*cthe(jt) + cthe(jt)**2
    wtmax = 2. + abs(asym)

  Else
    wt = 1.
    wtmax = 1.
  End If
!...Obtain correct angular distribution by rejection techniques.
  If (wt<rlu(0)*wtmax) Goto 420

!...Construct massive four-vectors using angles chosen. Mark decayed
!...resonances, add documentation lines. Shower evolution.
  500 Do jt = 1, jtmax
    If (kdcy(jt)==0) Goto 520
    id = iref(ip, jt)
    Call ludbrb(nsd(jt)+1, nsd(jt)+2, acos(cthe(jt)), phi(jt), dble(p(id,1)/p(id,4)), dble(p(id,2)/p(id,4)), dble(p(id,3)/p(id,4)))
    k(id, 1) = k(id, 1) + 10
    k(id, 4) = nsd(jt) + 1
    k(id, 5) = nsd(jt) + 2
    idoc = mint(83) + mint(4)
    Do i = nsd(jt) + 1, nsd(jt) + 2
      mint(4) = mint(4) + 1
      i1 = mint(83) + mint(4)
      k(i, 3) = i1
      k(i1, 1) = 21
      k(i1, 2) = k(i, 2)
      k(i1, 3) = iref(ip, jt+2)
      Do j = 1, 5
        p(i1, j) = p(i, j)
      End Do
    End Do
    If (jtmax==1) Then
      mint(7) = mint(83) + 6 + 2*iset(isub)
      mint(8) = mint(83) + 7 + 2*iset(isub)
    End If
!lin-8/19/02 avoid actual argument in common blocks of LUSHOW:
!      IF(MSTP(71).GE.1.AND.KDCY(JT).EQ.1) CALL LUSHOW(NSD(JT)+1,
!     &NSD(JT)+2,P(ID,5))
    pid5 = p(id, 5)
    If (mstp(71)>=1 .And. kdcy(jt)==1) Call lushow(nsd(jt)+1, nsd(jt)+2, pid5)

!...Check if new resonances were produced, loop back if needed.
    If (kdcy(jt)/=3) Goto 520
    np = np + 1
    iref(np, 1) = nsd(jt) + 1
    iref(np, 2) = nsd(jt) + 2
    iref(np, 3) = idoc + 1
    iref(np, 4) = idoc + 2
    iref(np, 5) = k(iref(ip,jt), 2)
    iref(np, 6) = iref(ip, jt)
  520 End Do
  530 If (ip<np) Goto 100

  Return
End Subroutine pyresd

!*********************************************************************

Subroutine pydiff

!...Handles diffractive and elastic scattering.
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  Save /lujets/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/

!...Reset K, P and V vectors. Store incoming particles.
  Do jt = 1, mstp(126) + 10
    i = mint(83) + jt
    Do j = 1, 5
      k(i, j) = 0
      p(i, j) = 0.
      v(i, j) = 0.
    End Do
  End Do
  n = mint(84)
  mint(3) = 0
  mint(21) = 0
  mint(22) = 0
  mint(23) = 0
  mint(24) = 0
  mint(4) = 4
  Do jt = 1, 2
    i = mint(83) + jt
    k(i, 1) = 21
    k(i, 2) = mint(10+jt)
    p(i, 5) = vint(2+jt)
    p(i, 3) = vint(5)*(-1)**(jt+1)
    p(i, 4) = sqrt(p(i,3)**2+p(i,5)**2)
  End Do
  mint(6) = 2

!...Subprocess; kinematics.
  isub = mint(1)
  sqlam = (vint(2)-vint(63)-vint(64))**2 - 4.*vint(63)*vint(64)
  pz = sqrt(sqlam)/(2.*vint(1))
  Do jt = 1, 2
    i = mint(83) + jt
    pe = (vint(2)+vint(62+jt)-vint(65-jt))/(2.*vint(1))

!...Elastically scattered particle.
    If (mint(16+jt)<=0) Then
      n = n + 1
      k(n, 1) = 1
      k(n, 2) = k(i, 2)
      k(n, 3) = i + 2
      p(n, 3) = pz*(-1)**(jt+1)
      p(n, 4) = pe
      p(n, 5) = p(i, 5)

!...Diffracted particle: valence quark kicked out.
    Else If (mstp(101)==1) Then
      n = n + 2
      k(n-1, 1) = 2
      k(n, 1) = 1
      k(n-1, 3) = i + 2
      k(n, 3) = i + 2
      Call pyspli(k(i,2), 21, k(n,2), k(n-1,2))
      p(n-1, 5) = ulmass(k(n-1,2))
      p(n, 5) = ulmass(k(n,2))
      sqlam = (vint(62+jt)-p(n-1,5)**2-p(n,5)**2)**2 - 4.*p(n-1, 5)**2*p(n, 5)**2
      p(n-1, 3) = (pe*sqrt(sqlam)+pz*(vint(62+jt)+p(n-1,5)**2-p(n,5)**2))/(2.*vint(62+jt))*(-1)**(jt+1)
      p(n-1, 4) = sqrt(p(n-1,3)**2+p(n-1,5)**2)
      p(n, 3) = pz*(-1)**(jt+1) - p(n-1, 3)
      p(n, 4) = sqrt(p(n,3)**2+p(n,5)**2)

!...Diffracted particle: gluon kicked out.
    Else
      n = n + 3
      k(n-2, 1) = 2
      k(n-1, 1) = 2
      k(n, 1) = 1
      k(n-2, 3) = i + 2
      k(n-1, 3) = i + 2
      k(n, 3) = i + 2
      Call pyspli(k(i,2), 21, k(n,2), k(n-2,2))
      k(n-1, 2) = 21
      p(n-2, 5) = ulmass(k(n-2,2))
      p(n-1, 5) = 0.
      p(n, 5) = ulmass(k(n,2))
!...Energy distribution for particle into two jets.
      120 imb = 1
      If (mod(k(i,2)/1000,10)/=0) imb = 2
      chik = parp(92+2*imb)
      If (mstp(92)<=1) Then
        If (imb==1) chi = rlu(0)
        If (imb==2) chi = 1. - sqrt(rlu(0))
      Else If (mstp(92)==2) Then
        chi = 1. - rlu(0)**(1./(1.+chik))
      Else If (mstp(92)==3) Then
        cut = 2.*0.3/vint(1)
        130 chi = rlu(0)**2
        If ((chi**2/(chi**2+cut**2))**0.25*(1.-chi)**chik<rlu(0)) Goto 130
      Else
        cut = 2.*0.3/vint(1)
        cutr = (1.+sqrt(1.+cut**2))/cut
        140 chir = cut*cutr**rlu(0)
        chi = (chir**2-cut**2)/(2.*chir)
        If ((1.-chi)**chik<rlu(0)) Goto 140
      End If
      If (chi<p(n,5)**2/vint(62+jt) .Or. chi>1.-p(n-2,5)**2/vint(62+jt)) Goto 120
      sqm = p(n-2, 5)**2/(1.-chi) + p(n, 5)**2/chi
      If ((sqrt(sqm)+parj(32))**2>=vint(62+jt)) Goto 120
      pzi = (pe*(vint(62+jt)-sqm)+pz*(vint(62+jt)+sqm))/(2.*vint(62+jt))
      pei = sqrt(pzi**2+sqm)
      pqqp = (1.-chi)*(pei+pzi)
      p(n-2, 3) = 0.5*(pqqp-p(n-2,5)**2/pqqp)*(-1)**(jt+1)
      p(n-2, 4) = sqrt(p(n-2,3)**2+p(n-2,5)**2)
      p(n-1, 3) = (pz-pzi)*(-1)**(jt+1)
      p(n-1, 4) = abs(p(n-1,3))
      p(n, 3) = pzi*(-1)**(jt+1) - p(n-2, 3)
      p(n, 4) = sqrt(p(n,3)**2+p(n,5)**2)
    End If

!...Documentation lines.
    k(i+2, 1) = 21
    If (mint(16+jt)==0) k(i+2, 2) = mint(10+jt)
    If (mint(16+jt)/=0) k(i+2, 2) = 10*(mint(10+jt)/10)
    k(i+2, 3) = i
    p(i+2, 3) = pz*(-1)**(jt+1)
    p(i+2, 4) = pe
    p(i+2, 5) = sqrt(vint(62+jt))
  End Do

!...Rotate outgoing partons/particles using cos(theta).
  Call ludbrb(mint(83)+3, n, acos(vint(23)), vint(24), 0D0, 0D0, 0D0)

  Return
End Subroutine pydiff

!*********************************************************************

Subroutine pyfram(iframe)

!...Performs transformations between different coordinate frames.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/

  If (iframe<1 .Or. iframe>2) Then
    Write (mstu(11), 1000) iframe, mint(6)
    Return
  End If
  If (iframe==mint(6)) Return

  If (mint(6)==1) Then
!...Transform from fixed target or user specified frame to
!...CM-frame of incoming particles.
    Call lurobo(0., 0., -vint(8), -vint(9), -vint(10))
    Call lurobo(0., -vint(7), 0., 0., 0.)
    Call lurobo(-vint(6), 0., 0., 0., 0.)
    mint(6) = 2

  Else
!...Transform from particle CM-frame to fixed target or user specified
!...frame.
    Call lurobo(vint(6), vint(7), vint(8), vint(9), vint(10))
    mint(6) = 1
  End If
  msti(6) = mint(6)

  Return

  1000 Format (1X, 'Error: illegal values in subroutine PYFRAM.', 1X, 'No transformation performed.'/1X, 'IFRAME =', 1X, I5, '; MINT(6) =', 1X, I5)
End Subroutine pyfram

!*********************************************************************

Subroutine pywidt(kflr, rmas, wdtp, wdte)

!...Calculates full and partial widths of resonances.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Dimension wdtp(0:40), wdte(0:40, 0:5)

!...Some common constants.
  kfla = iabs(kflr)
  sqm = rmas**2
  as = ulalps(sqm)
  aem = paru(101)
  xw = paru(102)
  radc = 1. + as/paru(1)

!...Reset width information.
  Do i = 0, 40
    wdtp(i) = 0.
    Do j = 0, 5
      wdte(i, j) = 0.
    End Do
  End Do

  If (kfla==21) Then
!...QCD:
    Do i = 1, mdcy(21, 3)
      idc = i + mdcy(21, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 110
      If (i<=8) Then
!...QCD -> q + qb
        wdtp(i) = (1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    110 End Do

  Else If (kfla==23) Then
!...Z0:
    If (mint(61)==1) Then
      ei = kchg(iabs(mint(15)), 1)/3.
      ai = sign(1., ei)
      vi = ai - 4.*ei*xw
      sqmz = pmas(23, 1)**2
      gzmz = pmas(23, 2)*pmas(23, 1)
      ggi = ei**2
      gzi = ei*vi/(8.*xw*(1.-xw))*sqm*(sqm-sqmz)/((sqm-sqmz)**2+gzmz**2)
      zzi = (vi**2+ai**2)/(16.*xw*(1.-xw))**2*sqm**2/((sqm-sqmz)**2+gzmz**2)
      If (mstp(43)==1) Then
!...Only gamma* production included
        gzi = 0.
        zzi = 0.
      Else If (mstp(43)==2) Then
!...Only Z0 production included
        ggi = 0.
        gzi = 0.
      End If
    Else If (mint(61)==2) Then
      vint(111) = 0.
      vint(112) = 0.
      vint(114) = 0.
    End If
    Do i = 1, mdcy(23, 3)
      idc = i + mdcy(23, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 120
      If (i<=8) Then
!...Z0 -> q + qb
        ef = kchg(i, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
        If (mint(61)==0) Then
          wdtp(i) = 3.*(vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==1) Then
          wdtp(i) = 3.*((ggi*ef**2+gzi*ef*vf+zzi*vf**2)*(1.+2.*rm1)+zzi*af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==2) Then
          ggf = 3.*ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          gzf = 3.*ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          zzf = 3.*(vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        End If
        wid2 = 1.
      Else If (i<=16) Then
!...Z0 -> l+ + l-, nu + nub
        ef = kchg(i+2, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
        wdtp(i) = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        If (mint(61)==0) Then
          wdtp(i) = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==1) Then
          wdtp(i) = ((ggi*ef**2+gzi*ef*vf+zzi*vf**2)*(1.+2.*rm1)+zzi*af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==2) Then
          ggf = ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzf = ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          zzf = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        End If
        wid2 = 1.
      Else
!...Z0 -> H+ + H-
        cf = 2.*(1.-2.*xw)
        If (mint(61)==0) Then
          wdtp(i) = 0.25*cf**2*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==1) Then
          wdtp(i) = 0.25*(ggi+gzi*cf+zzi*cf**2)*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==2) Then
          ggf = 0.25*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzf = 0.25*cf*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
          zzf = 0.25*cf**2*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        End If
        wid2 = wids(37, 1)
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
!lin-4/2008 modified a la pythia6115.f to avoid undefined values (GGF,GZF,ZZF):
!          VINT(111)=VINT(111)+GGF*WID2
!          VINT(112)=VINT(112)+GZF*WID2
!          VINT(114)=VINT(114)+ZZF*WID2
        If (mint(61)==2) Then
          vint(111) = vint(111) + ggf*wid2
          vint(112) = vint(112) + gzf*wid2
          vint(114) = vint(114) + zzf*wid2
        End If
!lin-4/2008-end
      End If
    120 End Do
    If (mstp(43)==1) Then
!...Only gamma* production included
      vint(112) = 0.
      vint(114) = 0.
    Else If (mstp(43)==2) Then
!...Only Z0 production included
      vint(111) = 0.
      vint(112) = 0.
    End If

  Else If (kfla==24) Then
!...W+/-:
    Do i = 1, mdcy(24, 3)
      idc = i + mdcy(24, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 130
      If (i<=16) Then
!...W+/- -> q + qb'
        wdtp(i) = 3.*(2.-rm1-rm2-(rm1-rm2)**2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))*vckm((i-1)/4+1, mod(i-1,4)+1)*radc
        wid2 = 1.
      Else
!...W+/- -> l+/- + nu
        wdtp(i) = (2.-rm1-rm2-(rm1-rm2)**2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    130 End Do

  Else If (kfla==25) Then
!...H0:
    Do i = 1, mdcy(25, 3)
      idc = i + mdcy(25, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 170
      If (i<=8) Then
!...H0 -> q + qb
        wdtp(i) = 3.*rm1*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
        wid2 = 1.
      Else If (i<=12) Then
!...H0 -> l+ + l-
        wdtp(i) = rm1*(1.-4.*rm1)*sqrt(max(0.,1.-4.*rm1))
        wid2 = 1.
      Else If (i==13) Then
!...H0 -> g + g; quark loop contribution only
        etare = 0.
        etaim = 0.
        Do j = 1, 2*mstp(1)
          eps = (2.*pmas(j,1)/rmas)**2
          If (eps<=1.) Then
            If (eps>1.E-4) Then
              root = sqrt(1.-eps)
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./eps-2.)
            End If
            phire = 0.25*(rln**2-paru(1)**2)
            phiim = 0.5*paru(1)*rln
          Else
            phire = -(asin(1./sqrt(eps)))**2
            phiim = 0.
          End If
          etare = etare + 0.5*eps*(1.+(eps-1.)*phire)
          etaim = etaim + 0.5*eps*(eps-1.)*phiim
        End Do
        eta2 = etare**2 + etaim**2
        wdtp(i) = (as/paru(1))**2*eta2
        wid2 = 1.
      Else If (i==14) Then
!...H0 -> gamma + gamma; quark, charged lepton and W loop contributions
        etare = 0.
        etaim = 0.
        Do j = 1, 3*mstp(1) + 1
          If (j<=2*mstp(1)) Then
            ej = kchg(j, 1)/3.
            eps = (2.*pmas(j,1)/rmas)**2
          Else If (j<=3*mstp(1)) Then
            jl = 2*(j-2*mstp(1)) - 1
            ej = kchg(10+jl, 1)/3.
            eps = (2.*pmas(10+jl,1)/rmas)**2
          Else
            eps = (2.*pmas(24,1)/rmas)**2
          End If
          If (eps<=1.) Then
            If (eps>1.E-4) Then
              root = sqrt(1.-eps)
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./eps-2.)
            End If
            phire = 0.25*(rln**2-paru(1)**2)
            phiim = 0.5*paru(1)*rln
          Else
            phire = -(asin(1./sqrt(eps)))**2
            phiim = 0.
          End If
          If (j<=2*mstp(1)) Then
            etare = etare + 0.5*3.*ej**2*eps*(1.+(eps-1.)*phire)
            etaim = etaim + 0.5*3.*ej**2*eps*(eps-1.)*phiim
          Else If (j<=3*mstp(1)) Then
            etare = etare + 0.5*ej**2*eps*(1.+(eps-1.)*phire)
            etaim = etaim + 0.5*ej**2*eps*(eps-1.)*phiim
          Else
            etare = etare - 0.5 - 0.75*eps*(1.+(eps-2.)*phire)
            etaim = etaim + 0.75*eps*(eps-2.)*phiim
          End If
        End Do
        eta2 = etare**2 + etaim**2
        wdtp(i) = (aem/paru(1))**2*0.5*eta2
        wid2 = 1.
      Else If (i==15) Then
!...H0 -> gamma + Z0; quark, charged lepton and W loop contributions
        etare = 0.
        etaim = 0.
        Do j = 1, 3*mstp(1) + 1
          If (j<=2*mstp(1)) Then
            ej = kchg(j, 1)/3.
            aj = sign(1., ej+0.1)
            vj = aj - 4.*ej*xw
            eps = (2.*pmas(j,1)/rmas)**2
            epsp = (2.*pmas(j,1)/pmas(23,1))**2
          Else If (j<=3*mstp(1)) Then
            jl = 2*(j-2*mstp(1)) - 1
            ej = kchg(10+jl, 1)/3.
            aj = sign(1., ej+0.1)
            vj = ai - 4.*ej*xw
            eps = (2.*pmas(10+jl,1)/rmas)**2
            epsp = (2.*pmas(10+jl,1)/pmas(23,1))**2
          Else
            eps = (2.*pmas(24,1)/rmas)**2
            epsp = (2.*pmas(24,1)/pmas(23,1))**2
          End If
          If (eps<=1.) Then
            root = sqrt(1.-eps)
            If (eps>1.E-4) Then
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./eps-2.)
            End If
            phire = 0.25*(rln**2-paru(1)**2)
            phiim = 0.5*paru(1)*rln
            psire = -(1.+0.5*root*rln)
            psiim = 0.5*paru(1)*root
          Else
            phire = -(asin(1./sqrt(eps)))**2
            phiim = 0.
            psire = -(1.+sqrt(eps-1.)*asin(1./sqrt(eps)))
            psiim = 0.
          End If
          If (epsp<=1.) Then
            root = sqrt(1.-epsp)
            If (epsp>1.E-4) Then
              rln = log((1.+root)/(1.-root))
            Else
              rln = log(4./epsp-2.)
            End If
            phirep = 0.25*(rln**2-paru(1)**2)
            phiimp = 0.5*paru(1)*rln
            psirep = -(1.+0.5*root*rln)
            psiimp = 0.5*paru(1)*root
          Else
            phirep = -(asin(1./sqrt(epsp)))**2
            phiimp = 0.
            psirep = -(1.+sqrt(epsp-1.)*asin(1./sqrt(epsp)))
            psiimp = 0.
          End If
          fxyre = eps*epsp/(8.*(eps-epsp))*(1.-eps*epsp/(eps-epsp)*(phire-phirep)+2.*eps/(eps-epsp)*(psire-psirep))
          fxyim = eps*epsp/(8.*(eps-epsp))*(-eps*epsp/(eps-epsp)*(phiim-phiimp)+2.*eps/(eps-epsp)*(psiim-psiimp))
          f1re = eps*epsp/(2.*(eps-epsp))*(phire-phirep)
          f1im = eps*epsp/(2.*(eps-epsp))*(phiim-phiimp)
          If (j<=2*mstp(1)) Then
            etare = etare - 3.*ej*vj*(fxyre-0.25*f1re)
            etaim = etaim - 3.*ej*vj*(fxyim-0.25*f1im)
          Else If (j<=3*mstp(1)) Then
            etare = etare - ej*vj*(fxyre-0.25*f1re)
            etaim = etaim - ej*vj*(fxyim-0.25*f1im)
          Else
            etare = etare - sqrt(1.-xw)*(((1.+2./eps)*xw/sqrt(1.-xw)-(5.+2./eps))*fxyre+(3.-xw/sqrt(1.-xw))*f1re)
            etaim = etaim - sqrt(1.-xw)*(((1.+2./eps)*xw/sqrt(1.-xw)-(5.+2./eps))*fxyim+(3.-xw/sqrt(1.-xw))*f1im)
          End If
        End Do
        eta2 = etare**2 + etaim**2
        wdtp(i) = (aem/paru(1))**2*(1.-(pmas(23,1)/rmas)**2)**3/xw*eta2
        wid2 = wids(23, 2)
      Else
!...H0 -> Z0 + Z0, W+ + W-
        wdtp(i) = (1.-4.*rm1+12.*rm1**2)*sqrt(max(0.,1.-4.*rm1))/(2.*(18-i))
        wid2 = wids(7+i, 1)
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    170 End Do

  Else If (kfla==32) Then
!...Z'0:
    If (mint(61)==1) Then
      ei = kchg(iabs(mint(15)), 1)/3.
      ai = sign(1., ei)
      vi = ai - 4.*ei*xw
      sqmz = pmas(23, 1)**2
      gzmz = pmas(23, 2)*pmas(23, 1)
      api = sign(1., ei)
      vpi = api - 4.*ei*xw
      sqmzp = pmas(32, 1)**2
      gzpmzp = pmas(32, 2)*pmas(32, 1)
      ggi = ei**2
      gzi = ei*vi/(8.*xw*(1.-xw))*sqm*(sqm-sqmz)/((sqm-sqmz)**2+gzmz**2)
      gzpi = ei*vpi/(8.*xw*(1.-xw))*sqm*(sqm-sqmzp)/((sqm-sqmzp)**2+gzpmzp**2)
      zzi = (vi**2+ai**2)/(16.*xw*(1.-xw))**2*sqm**2/((sqm-sqmz)**2+gzmz**2)
      zzpi = 2.*(vi*vpi+ai*api)/(16.*xw*(1.-xw))**2*sqm**2*((sqm-sqmz)*(sqm-sqmzp)+gzmz*gzpmzp)/(((sqm-sqmz)**2+gzmz**2)*((sqm-sqmzp)**2+gzpmzp**2))
      zpzpi = (vpi**2+api**2)/(16.*xw*(1.-xw))**2*sqm**2/((sqm-sqmzp)**2+gzpmzp**2)
      If (mstp(44)==1) Then
!...Only gamma* production included
        gzi = 0.
        gzpi = 0.
        zzi = 0.
        zzpi = 0.
        zpzpi = 0.
      Else If (mstp(44)==2) Then
!...Only Z0 production included
        ggi = 0.
        gzi = 0.
        gzpi = 0.
        zzpi = 0.
        zpzpi = 0.
      Else If (mstp(44)==3) Then
!...Only Z'0 production included
        ggi = 0.
        gzi = 0.
        gzpi = 0.
        zzi = 0.
        zzpi = 0.
      Else If (mstp(44)==4) Then
!...Only gamma*/Z0 production included
        gzpi = 0.
        zzpi = 0.
        zpzpi = 0.
      Else If (mstp(44)==5) Then
!...Only gamma*/Z'0 production included
        gzi = 0.
        zzi = 0.
        zzpi = 0.
      Else If (mstp(44)==6) Then
!...Only Z0/Z'0 production included
        ggi = 0.
        gzi = 0.
        gzpi = 0.
      End If
    Else If (mint(61)==2) Then
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
      vint(114) = 0.
      vint(115) = 0.
      vint(116) = 0.
    End If
    Do i = 1, mdcy(32, 3)
      idc = i + mdcy(32, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 180
      If (i<=8) Then
!...Z'0 -> q + qb
        ef = kchg(i, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
        apf = sign(1., ef+0.1)
        vpf = apf - 4.*ef*xw
        If (mint(61)==0) Then
          wdtp(i) = 3.*(vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==1) Then
          wdtp(i) = 3.*((ggi*ef**2+gzi*ef*vf+gzpi*ef*vpf+zzi*vf**2+zzpi*vf*vpf+zpzpi*vpf**2)*(1.+2.*rm1)+(zzi*af**2+zzpi*af*apf+zpzpi*apf**2)*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        Else If (mint(61)==2) Then
          ggf = 3.*ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          gzf = 3.*ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          gzpf = 3.*ef*vpf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))*radc
          zzf = 3.*(vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
          zzpf = 3.*(vf*vpf*(1.+2.*rm1)+af*apf*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
          zpzpf = 3.*(vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))*radc
        End If
        wid2 = 1.
      Else
!...Z'0 -> l+ + l-, nu + nub
        ef = kchg(i+2, 1)/3.
        af = sign(1., ef+0.1)
        vf = af - 4.*ef*xw
!lin-4/2008 modified above a la pythia6115.f to avoid undefined variable API:
!          APF=SIGN(1.,EF+0.1)
!          VPF=API-4.*EF*XW
        If (i<=10) Then
          vpf = paru(127-2*mod(i,2))
          apf = paru(128-2*mod(i,2))
        Else If (i<=12) Then
          vpf = parj(186-2*mod(i,2))
          apf = parj(187-2*mod(i,2))
        Else
          vpf = parj(194-2*mod(i,2))
          apf = parj(195-2*mod(i,2))
        End If
!lin-4/2008-end
        If (mint(61)==0) Then
          wdtp(i) = (vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==1) Then
          wdtp(i) = ((ggi*ef**2+gzi*ef*vf+gzpi*ef*vpf+zzi*vf**2+zzpi*vf*vpf+zpzpi*vpf**2)*(1.+2.*rm1)+(zzi*af**2+zzpi*af*apf+zpzpi*apf**2)*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        Else If (mint(61)==2) Then
          ggf = ef**2*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzf = ef*vf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          gzpf = ef*vpf*(1.+2.*rm1)*sqrt(max(0.,1.-4.*rm1))
          zzf = (vf**2*(1.+2.*rm1)+af**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
          zzpf = (vf*vpf*(1.+2.*rm1)+af*apf*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
          zpzpf = (vpf**2*(1.+2.*rm1)+apf**2*(1.-4.*rm1))*sqrt(max(0.,1.-4.*rm1))
        End If
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
!lin-4/2008:
!          VINT(111)=VINT(111)+GGF
!          VINT(112)=VINT(112)+GZF
!          VINT(113)=VINT(113)+GZPF
!          VINT(114)=VINT(114)+ZZF
!          VINT(115)=VINT(115)+ZZPF
!          VINT(116)=VINT(116)+ZPZPF
        If (mint(61)==2) Then
          vint(111) = vint(111) + ggf
          vint(112) = vint(112) + gzf
          vint(113) = vint(113) + gzpf
          vint(114) = vint(114) + zzf
          vint(115) = vint(115) + zzpf
          vint(116) = vint(116) + zpzpf
        End If
!lin-4/2008-end
      End If
    180 End Do
    If (mstp(44)==1) Then
!...Only gamma* production included
      vint(112) = 0.
      vint(113) = 0.
      vint(114) = 0.
      vint(115) = 0.
      vint(116) = 0.
    Else If (mstp(44)==2) Then
!...Only Z0 production included
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
      vint(115) = 0.
      vint(116) = 0.
    Else If (mstp(44)==3) Then
!...Only Z'0 production included
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
      vint(114) = 0.
      vint(115) = 0.
    Else If (mstp(44)==4) Then
!...Only gamma*/Z0 production included
      vint(113) = 0.
      vint(115) = 0.
      vint(116) = 0.
    Else If (mstp(44)==5) Then
!...Only gamma*/Z'0 production included
      vint(112) = 0.
      vint(114) = 0.
      vint(115) = 0.
    Else If (mstp(44)==6) Then
!...Only Z0/Z'0 production included
      vint(111) = 0.
      vint(112) = 0.
      vint(113) = 0.
    End If

  Else If (kfla==37) Then
!...H+/-:
    Do i = 1, mdcy(37, 3)
      idc = i + mdcy(37, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 190
      If (i<=4) Then
!...H+/- -> q + qb'
        wdtp(i) = 3.*((rm1*paru(121)+rm2/paru(121))*(1.-rm1-rm2)-4.*rm1*rm2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))*radc
        wid2 = 1.
      Else
!...H+/- -> l+/- + nu
        wdtp(i) = ((rm1*paru(121)+rm2/paru(121))*(1.-rm1-rm2)-4.*rm1*rm2)*sqrt(max(0.,(1.-rm1-rm2)**2-4.*rm1*rm2))
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    190 End Do

  Else If (kfla==40) Then
!...R:
    Do i = 1, mdcy(40, 3)
      idc = i + mdcy(40, 2) - 1
      rm1 = (pmas(iabs(kfdp(idc,1)),1)/rmas)**2
      rm2 = (pmas(iabs(kfdp(idc,2)),1)/rmas)**2
      If (sqrt(rm1)+sqrt(rm2)>1. .Or. mdme(idc,1)<0) Goto 200
      If (i<=4) Then
!...R -> q + qb'
        wdtp(i) = 3.*radc
        wid2 = 1.
      Else
!...R -> l+ + l'-
        wdtp(i) = 1.
        wid2 = 1.
      End If
      wdtp(0) = wdtp(0) + wdtp(i)
      If (mdme(idc,1)>0) Then
        wdte(i, mdme(idc,1)) = wdtp(i)*wid2
        wdte(0, mdme(idc,1)) = wdte(0, mdme(idc,1)) + wdte(i, mdme(idc,1))
        wdte(i, 0) = wdte(i, mdme(idc,1))
        wdte(0, 0) = wdte(0, 0) + wdte(i, 0)
      End If
    200 End Do

  End If
  mint(61) = 0

  Return
End Subroutine pywidt

!***********************************************************************

!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine pyklim(ilim)

!...Checks generated variables against pre-set kinematical limits;
!...also calculates limits on variables used in generation.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/

!...Common kinematical expressions.
  isub = mint(1)
  If (isub==96) Goto 110
  sqm3 = vint(63)
  sqm4 = vint(64)
  If (ilim/=1) Then
    tau = vint(21)
    rm3 = sqm3/(tau*vint(2))
    rm4 = sqm4/(tau*vint(2))
    be34 = sqrt((1.-rm3-rm4)**2-4.*rm3*rm4)
  End If
  pthmin = ckin(3)
  If (min(sqm3,sqm4)<ckin(6)**2) pthmin = max(ckin(3), ckin(5))
  If (ilim==0) Then
!...Check generated values of tau, y*, cos(theta-hat), and tau' against
!...pre-set kinematical limits.
    yst = vint(22)
    cth = vint(23)
    taup = vint(26)
    If (iset(isub)<=2) Then
      x1 = sqrt(tau)*exp(yst)
      x2 = sqrt(tau)*exp(-yst)
    Else
      x1 = sqrt(taup)*exp(yst)
      x2 = sqrt(taup)*exp(-yst)
    End If
    xf = x1 - x2
    If (tau*vint(2)<ckin(1)**2) mint(51) = 1
    If (ckin(2)>=0. .And. tau*vint(2)>ckin(2)**2) mint(51) = 1
    If (x1<ckin(21) .Or. x1>ckin(22)) mint(51) = 1
    If (x2<ckin(23) .Or. x2>ckin(24)) mint(51) = 1
    If (xf<ckin(25) .Or. xf>ckin(26)) mint(51) = 1
    If (yst<ckin(7) .Or. yst>ckin(8)) mint(51) = 1
    If (iset(isub)==2 .Or. iset(isub)==4) Then
      pth = 0.5*be34*sqrt(tau*vint(2)*(1.-cth**2))
      y3 = yst + 0.5*log((1.+rm3-rm4+be34*cth)/(1.+rm3-rm4-be34*cth))
      y4 = yst + 0.5*log((1.+rm4-rm3-be34*cth)/(1.+rm4-rm3+be34*cth))
      ylarge = max(y3, y4)
      ysmall = min(y3, y4)
      etalar = 10.
      etasma = -10.
      sth = sqrt(1.-cth**2)
      If (sth<1.E-6) Goto 100
      expet3 = ((1.+rm3-rm4)*sinh(yst)+be34*cosh(yst)*cth+sqrt(((1.+rm3-rm4)*cosh(yst)+be34*sinh(yst)*cth)**2-4.*rm3))/(be34*sth)
      expet4 = ((1.-rm3+rm4)*sinh(yst)-be34*cosh(yst)*cth+sqrt(((1.-rm3+rm4)*cosh(yst)-be34*sinh(yst)*cth)**2-4.*rm4))/(be34*sth)
      eta3 = log(min(1.E10,max(1.E-10,expet3)))
      eta4 = log(min(1.E10,max(1.E-10,expet4)))
      etalar = max(eta3, eta4)
      etasma = min(eta3, eta4)
      100 cts3 = ((1.+rm3-rm4)*sinh(yst)+be34*cosh(yst)*cth)/sqrt(((1.+rm3-rm4)*cosh(yst)+be34*sinh(yst)*cth)**2-4.*rm3)
      cts4 = ((1.-rm3+rm4)*sinh(yst)-be34*cosh(yst)*cth)/sqrt(((1.-rm3+rm4)*cosh(yst)-be34*sinh(yst)*cth)**2-4.*rm4)
      ctslar = max(cts3, cts4)
      ctssma = min(cts3, cts4)
      If (pth<pthmin) mint(51) = 1
      If (ckin(4)>=0. .And. pth>ckin(4)) mint(51) = 1
      If (ylarge<ckin(9) .Or. ylarge>ckin(10)) mint(51) = 1
      If (ysmall<ckin(11) .Or. ysmall>ckin(12)) mint(51) = 1
      If (etalar<ckin(13) .Or. etalar>ckin(14)) mint(51) = 1
      If (etasma<ckin(15) .Or. etasma>ckin(16)) mint(51) = 1
      If (ctslar<ckin(17) .Or. ctslar>ckin(18)) mint(51) = 1
      If (ctssma<ckin(19) .Or. ctssma>ckin(20)) mint(51) = 1
      If (cth<ckin(27) .Or. cth>ckin(28)) mint(51) = 1
    End If
    If (iset(isub)==3 .Or. iset(isub)==4) Then
      If (taup*vint(2)<ckin(31)**2) mint(51) = 1
      If (ckin(32)>=0. .And. taup*vint(2)>ckin(32)**2) mint(51) = 1
    End If

  Else If (ilim==1) Then
!...Calculate limits on tau
!...0) due to definition
    taumn0 = 0.
    taumx0 = 1.
!...1) due to limits on subsystem mass
    taumn1 = ckin(1)**2/vint(2)
    taumx1 = 1.
    If (ckin(2)>=0.) taumx1 = ckin(2)**2/vint(2)
!...2) due to limits on pT-hat (and non-overlapping rapidity intervals)
    tm3 = sqrt(sqm3+pthmin**2)
    tm4 = sqrt(sqm4+pthmin**2)
    ydcosh = 1.
    If (ckin(9)>ckin(12)) ydcosh = cosh(ckin(9)-ckin(12))
    taumn2 = (tm3**2+2.*tm3*tm4*ydcosh+tm4**2)/vint(2)
    taumx2 = 1.
!...3) due to limits on pT-hat and cos(theta-hat)
    cth2mn = min(ckin(27)**2, ckin(28)**2)
    cth2mx = max(ckin(27)**2, ckin(28)**2)
    taumn3 = 0.
    If (ckin(27)*ckin(28)>0.) taumn3 = (sqrt(sqm3+pthmin**2/(1.-cth2mn))+sqrt(sqm4+pthmin**2/(1.-cth2mn)))**2/vint(2)
    taumx3 = 1.
    If (ckin(4)>=0. .And. cth2mx<1.) taumx3 = (sqrt(sqm3+ckin(4)**2/(1.-cth2mx))+sqrt(sqm4+ckin(4)**2/(1.-cth2mx)))**2/vint(2)
!...4) due to limits on x1 and x2
    taumn4 = ckin(21)*ckin(23)
    taumx4 = ckin(22)*ckin(24)
!...5) due to limits on xF
    taumn5 = 0.
    taumx5 = max(1.-ckin(25), 1.+ckin(26))
    vint(11) = max(taumn0, taumn1, taumn2, taumn3, taumn4, taumn5)
    vint(31) = min(taumx0, taumx1, taumx2, taumx3, taumx4, taumx5)
    If (mint(43)==1 .And. (iset(isub)==1 .Or. iset(isub)==2)) Then
      vint(11) = 0.99999
      vint(31) = 1.00001
    End If
    If (vint(31)<=vint(11)) mint(51) = 1

  Else If (ilim==2) Then
!...Calculate limits on y*
    If (iset(isub)==3 .Or. iset(isub)==4) tau = vint(26)
    taurt = sqrt(tau)
!...0) due to kinematics
    ystmn0 = log(taurt)
    ystmx0 = -ystmn0
!...1) due to explicit limits
    ystmn1 = ckin(7)
    ystmx1 = ckin(8)
!...2) due to limits on x1
    ystmn2 = log(max(tau,ckin(21))/taurt)
    ystmx2 = log(max(tau,ckin(22))/taurt)
!...3) due to limits on x2
    ystmn3 = -log(max(tau,ckin(24))/taurt)
    ystmx3 = -log(max(tau,ckin(23))/taurt)
!...4) due to limits on xF
    yepmn4 = 0.5*abs(ckin(25))/taurt
    ystmn4 = sign(log(sqrt(1.+yepmn4**2)+yepmn4), ckin(25))
    yepmx4 = 0.5*abs(ckin(26))/taurt
    ystmx4 = sign(log(sqrt(1.+yepmx4**2)+yepmx4), ckin(26))
!...5) due to simultaneous limits on y-large and y-small
    yepsmn = (rm3-rm4)*sinh(ckin(9)-ckin(11))
    yepsmx = (rm3-rm4)*sinh(ckin(10)-ckin(12))
    ydifmn = abs(log(sqrt(1.+yepsmn**2)-yepsmn))
    ydifmx = abs(log(sqrt(1.+yepsmx**2)-yepsmx))
    ystmn5 = 0.5*(ckin(9)+ckin(11)-ydifmn)
    ystmx5 = 0.5*(ckin(10)+ckin(12)+ydifmx)
!...6) due to simultaneous limits on cos(theta-hat) and y-large or
!...   y-small
    cthlim = sqrt(1.-4.*pthmin**2/(be34*tau*vint(2)))
    rzmn = be34*max(ckin(27), -cthlim)
    rzmx = be34*min(ckin(28), cthlim)
    yex3mx = (1.+rm3-rm4+rzmx)/max(1E-10, 1.+rm3-rm4-rzmx)
    yex4mx = (1.+rm4-rm3-rzmn)/max(1E-10, 1.+rm4-rm3+rzmn)
    yex3mn = max(1E-10, 1.+rm3-rm4+rzmn)/(1.+rm3-rm4-rzmn)
    yex4mn = max(1E-10, 1.+rm4-rm3-rzmx)/(1.+rm4-rm3+rzmx)
    ystmn6 = ckin(9) - 0.5*log(max(yex3mx,yex4mx))
    ystmx6 = ckin(12) - 0.5*log(min(yex3mn,yex4mn))
    vint(12) = max(ystmn0, ystmn1, ystmn2, ystmn3, ystmn4, ystmn5, ystmn6)
    vint(32) = min(ystmx0, ystmx1, ystmx2, ystmx3, ystmx4, ystmx5, ystmx6)
    If (mint(43)==1) Then
      vint(12) = -0.00001
      vint(32) = 0.00001
    Else If (mint(43)==2) Then
      vint(12) = 0.99999*ystmx0
      vint(32) = 1.00001*ystmx0
    Else If (mint(43)==3) Then
      vint(12) = -1.00001*ystmx0
      vint(32) = -0.99999*ystmx0
    End If
    If (vint(32)<=vint(12)) mint(51) = 1

  Else If (ilim==3) Then
!...Calculate limits on cos(theta-hat)
    yst = vint(22)
!...0) due to definition
    ctnmn0 = -1.
    ctnmx0 = 0.
    ctpmn0 = 0.
    ctpmx0 = 1.
!...1) due to explicit limits
    ctnmn1 = min(0., ckin(27))
    ctnmx1 = min(0., ckin(28))
    ctpmn1 = max(0., ckin(27))
    ctpmx1 = max(0., ckin(28))
!...2) due to limits on pT-hat
    ctnmn2 = -sqrt(1.-4.*pthmin**2/(be34**2*tau*vint(2)))
    ctpmx2 = -ctnmn2
    ctnmx2 = 0.
    ctpmn2 = 0.
    If (ckin(4)>=0.) Then
      ctnmx2 = -sqrt(max(0.,1.-4.*ckin(4)**2/(be34**2*tau*vint(2))))
      ctpmn2 = -ctnmx2
    End If
!...3) due to limits on y-large and y-small
    ctnmn3 = min(0., max((1.+rm3-rm4)/be34*tanh(ckin(11)-yst),-(1.-rm3+rm4)/be34*tanh(ckin(10)-yst)))
    ctnmx3 = min(0., (1.+rm3-rm4)/be34*tanh(ckin(12)-yst), -(1.-rm3+rm4)/be34*tanh(ckin(9)-yst))
    ctpmn3 = max(0., (1.+rm3-rm4)/be34*tanh(ckin(9)-yst), -(1.-rm3+rm4)/be34*tanh(ckin(12)-yst))
    ctpmx3 = max(0., min((1.+rm3-rm4)/be34*tanh(ckin(10)-yst),-(1.-rm3+rm4)/be34*tanh(ckin(11)-yst)))
    vint(13) = max(ctnmn0, ctnmn1, ctnmn2, ctnmn3)
    vint(33) = min(ctnmx0, ctnmx1, ctnmx2, ctnmx3)
    vint(14) = max(ctpmn0, ctpmn1, ctpmn2, ctpmn3)
    vint(34) = min(ctpmx0, ctpmx1, ctpmx2, ctpmx3)
    If (vint(33)<=vint(13) .And. vint(34)<=vint(14)) mint(51) = 1

  Else If (ilim==4) Then
!...Calculate limits on tau'
!...0) due to kinematics
    tapmn0 = tau
    tapmx0 = 1.
!...1) due to explicit limits
    tapmn1 = ckin(31)**2/vint(2)
    tapmx1 = 1.
    If (ckin(32)>=0.) tapmx1 = ckin(32)**2/vint(2)
    vint(16) = max(tapmn0, tapmn1)
    vint(36) = min(tapmx0, tapmx1)
    If (mint(43)==1) Then
      vint(16) = 0.99999
      vint(36) = 1.00001
    End If
    If (vint(36)<=vint(16)) mint(51) = 1

  End If
  Return

!...Special case for low-pT and multiple interactions:
!...effective kinematical limits for tau, y*, cos(theta-hat).
  110 If (ilim==0) Then
  Else If (ilim==1) Then
    If (mstp(82)<=1) vint(11) = 4.*parp(81)**2/vint(2)
    If (mstp(82)>=2) vint(11) = parp(82)**2/vint(2)
    vint(31) = 1.
  Else If (ilim==2) Then
    vint(12) = 0.5*log(vint(21))
    vint(32) = -vint(12)
  Else If (ilim==3) Then
    If (mstp(82)<=1) st2eff = 4.*parp(81)**2/(vint(21)*vint(2))
    If (mstp(82)>=2) st2eff = 0.01*parp(82)**2/(vint(21)*vint(2))
    vint(13) = -sqrt(max(0.,1.-st2eff))
    vint(33) = 0.
    vint(14) = 0.
    vint(34) = -vint(13)
  End If

  Return
End Subroutine pyklim

!*********************************************************************

Subroutine pykmap(ivar, mvar, vvar)

!...Maps a uniform distribution into a distribution of a kinematical
!...variable according to one of the possibilities allowed. It is
!...assumed that kinematical limits have been set by a PYKLIM call.
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/

!...Convert VVAR to tau variable.
  isub = mint(1)
  If (ivar==1) Then
    taumin = vint(11)
    taumax = vint(31)
    If (mvar==3 .Or. mvar==4) Then
      taure = vint(73)
      gamre = vint(74)
    Else If (mvar==5 .Or. mvar==6) Then
      taure = vint(75)
      gamre = vint(76)
    End If
    If (mint(43)==1 .And. (iset(isub)==1 .Or. iset(isub)==2)) Then
      tau = 1.
    Else If (mvar==1) Then
      tau = taumin*(taumax/taumin)**vvar
    Else If (mvar==2) Then
      tau = taumax*taumin/(taumin+(taumax-taumin)*vvar)
    Else If (mvar==3 .Or. mvar==5) Then
      ratgen = (taure+taumax)/(taure+taumin)*taumin/taumax
      tau = taure*taumin/((taure+taumin)*ratgen**vvar-taumin)
    Else
      aupp = atan((taumax-taure)/gamre)
      alow = atan((taumin-taure)/gamre)
      tau = taure + gamre*tan(alow+(aupp-alow)*vvar)
    End If
    vint(21) = min(taumax, max(taumin,tau))

!...Convert VVAR to y* variable.
  Else If (ivar==2) Then
    ystmin = vint(12)
    ystmax = vint(32)
    If (mint(43)==1) Then
      yst = 0.
    Else If (mint(43)==2) Then
      If (iset(isub)<=2) yst = -0.5*log(vint(21))
      If (iset(isub)>=3) yst = -0.5*log(vint(26))
    Else If (mint(43)==3) Then
      If (iset(isub)<=2) yst = 0.5*log(vint(21))
      If (iset(isub)>=3) yst = 0.5*log(vint(26))
    Else If (mvar==1) Then
      yst = ystmin + (ystmax-ystmin)*sqrt(vvar)
    Else If (mvar==2) Then
      yst = ystmax - (ystmax-ystmin)*sqrt(1.-vvar)
    Else
      aupp = atan(exp(ystmax))
      alow = atan(exp(ystmin))
      yst = log(tan(alow+(aupp-alow)*vvar))
    End If
    vint(22) = min(ystmax, max(ystmin,yst))

!...Convert VVAR to cos(theta-hat) variable.
  Else If (ivar==3) Then
    rm34 = 2.*vint(63)*vint(64)/(vint(21)*vint(2))**2
    rsqm = 1. + rm34
    If (2.*vint(71)**2/(vint(21)*vint(2))<0.0001) rm34 = max(rm34, 2.*vint(71)**2/(vint(21)*vint(2)))
    ctnmin = vint(13)
    ctnmax = vint(33)
    ctpmin = vint(14)
    ctpmax = vint(34)
    If (mvar==1) Then
      aneg = ctnmax - ctnmin
      apos = ctpmax - ctpmin
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = ctnmin + (ctnmax-ctnmin)*vctn
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = ctpmin + (ctpmax-ctpmin)*vctp
      End If
    Else If (mvar==2) Then
      rmnmin = max(rm34, rsqm-ctnmin)
      rmnmax = max(rm34, rsqm-ctnmax)
      rmpmin = max(rm34, rsqm-ctpmin)
      rmpmax = max(rm34, rsqm-ctpmax)
      aneg = log(rmnmin/rmnmax)
      apos = log(rmpmin/rmpmax)
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = rsqm - rmnmin*(rmnmax/rmnmin)**vctn
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = rsqm - rmpmin*(rmpmax/rmpmin)**vctp
      End If
    Else If (mvar==3) Then
      rmnmin = max(rm34, rsqm+ctnmin)
      rmnmax = max(rm34, rsqm+ctnmax)
      rmpmin = max(rm34, rsqm+ctpmin)
      rmpmax = max(rm34, rsqm+ctpmax)
      aneg = log(rmnmax/rmnmin)
      apos = log(rmpmax/rmpmin)
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = rmnmin*(rmnmax/rmnmin)**vctn - rsqm
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = rmpmin*(rmpmax/rmpmin)**vctp - rsqm
      End If
    Else If (mvar==4) Then
      rmnmin = max(rm34, rsqm-ctnmin)
      rmnmax = max(rm34, rsqm-ctnmax)
      rmpmin = max(rm34, rsqm-ctpmin)
      rmpmax = max(rm34, rsqm-ctpmax)
      aneg = 1./rmnmax - 1./rmnmin
      apos = 1./rmpmax - 1./rmpmin
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = rsqm - 1./(1./rmnmin+aneg*vctn)
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = rsqm - 1./(1./rmpmin+apos*vctp)
      End If
    Else If (mvar==5) Then
      rmnmin = max(rm34, rsqm+ctnmin)
      rmnmax = max(rm34, rsqm+ctnmax)
      rmpmin = max(rm34, rsqm+ctpmin)
      rmpmax = max(rm34, rsqm+ctpmax)
      aneg = 1./rmnmin - 1./rmnmax
      apos = 1./rmpmin - 1./rmpmax
      If (aneg>0. .And. vvar*(aneg+apos)<=aneg) Then
        vctn = vvar*(aneg+apos)/aneg
        cth = 1./(1./rmnmin-aneg*vctn) - rsqm
      Else
        vctp = (vvar*(aneg+apos)-aneg)/apos
        cth = 1./(1./rmpmin-apos*vctp) - rsqm
      End If
    End If
    If (cth<0.) cth = min(ctnmax, max(ctnmin,cth))
    If (cth>0.) cth = min(ctpmax, max(ctpmin,cth))
    vint(23) = cth

!...Convert VVAR to tau' variable.
  Else If (ivar==4) Then
    tau = vint(11)
    taupmn = vint(16)
    taupmx = vint(36)
    If (mint(43)==1) Then
      taup = 1.
    Else If (mvar==1) Then
      taup = taupmn*(taupmx/taupmn)**vvar
    Else
      aupp = (1.-tau/taupmx)**4
      alow = (1.-tau/taupmn)**4
      taup = tau/(1.-(alow+(aupp-alow)*vvar)**0.25)
    End If
    vint(26) = min(taupmx, max(taupmn,taup))
  End If

  Return
End Subroutine pykmap

!***********************************************************************

Subroutine pysigh(nchn, sigs)

!...Differential matrix elements for all included subprocesses.
!...Note that what is coded is (disregarding the COMFAC factor)
!...1) for 2 -> 1 processes: s-hat/pi*d(sigma-hat), where,
!...when d(sigma-hat) is given in the zero-width limit, the delta
!...function in tau is replaced by a Breit-Wigner:
!...1/pi*(s*m_res*Gamma_res)/((s*tau-m_res^2)^2+(m_res*Gamma_res)^2);
!...2) for 2 -> 2 processes: (s-hat)**2/pi*d(sigma-hat)/d(t-hat);
!...i.e., dimensionless quantities. COMFAC contains the factor
!...pi/s and the conversion factor from GeV^-2 to mb.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  Save /ludat3/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Dimension x(2), xpq(-6:6), kfac(2, -40:40), wdtp(0:40), wdte(0:40, 0:5)

!...Reset number of channels and cross-section.
  nchn = 0
  sigs = 0.

!...Read kinematical variables and limits.
  isub = mint(1)
  taumin = vint(11)
  ystmin = vint(12)
  ctnmin = vint(13)
  ctpmin = vint(14)
  xt2min = vint(15)
  taupmn = vint(16)
  tau = vint(21)
  yst = vint(22)
  cth = vint(23)
  xt2 = vint(25)
  taup = vint(26)
  taumax = vint(31)
  ystmax = vint(32)
  ctnmax = vint(33)
  ctpmax = vint(34)
  xt2max = vint(35)
  taupmx = vint(36)

!...Derive kinematical quantities.
  If (iset(isub)<=2 .Or. iset(isub)==5) Then
    x(1) = sqrt(tau)*exp(yst)
    x(2) = sqrt(tau)*exp(-yst)
  Else
    x(1) = sqrt(taup)*exp(yst)
    x(2) = sqrt(taup)*exp(-yst)
  End If
  If (mint(43)==4 .And. iset(isub)>=1 .And. (x(1)>0.999 .Or. x(2)>0.999)) Return
  sh = tau*vint(2)
  sqm3 = vint(63)
  sqm4 = vint(64)
  rm3 = sqm3/sh
  rm4 = sqm4/sh
  be34 = sqrt((1.-rm3-rm4)**2-4.*rm3*rm4)
  rpts = 4.*vint(71)**2/sh
  be34l = sqrt(max(0.,(1.-rm3-rm4)**2-4.*rm3*rm4-rpts))
  rm34 = 2.*rm3*rm4
  rsqm = 1. + rm34
  rthm = (4.*rm3*rm4+rpts)/(1.-rm3-rm4+be34l)
  th = -0.5*sh*max(rthm, 1.-rm3-rm4-be34*cth)
  uh = -0.5*sh*max(rthm, 1.-rm3-rm4+be34*cth)
  sqpth = 0.25*sh*be34**2*(1.-cth**2)
  sh2 = sh**2
  th2 = th**2
  uh2 = uh**2

!...Choice of Q2 scale.
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    q2 = sh
  Else If (mod(iset(isub),2)==0 .Or. iset(isub)==5) Then
    If (mstp(32)==1) Then
      q2 = 2.*sh*th*uh/(sh**2+th**2+uh**2)
    Else If (mstp(32)==2) Then
      q2 = sqpth + 0.5*(sqm3+sqm4)
    Else If (mstp(32)==3) Then
      q2 = min(-th, -uh)
    Else If (mstp(32)==4) Then
      q2 = sh
    End If
    If (iset(isub)==5 .And. mstp(82)>=2) q2 = q2 + parp(82)**2
  End If

!...Store derived kinematical quantities.
  vint(41) = x(1)
  vint(42) = x(2)
  vint(44) = sh
  vint(43) = sqrt(sh)
  vint(45) = th
  vint(46) = uh
  vint(48) = sqpth
  vint(47) = sqrt(sqpth)
  vint(50) = taup*vint(2)
  vint(49) = sqrt(max(0.,vint(50)))
  vint(52) = q2
  vint(51) = sqrt(q2)

!...Calculate parton structure functions.
  If (iset(isub)<=0) Goto 145
  If (mint(43)>=2) Then
    q2sf = q2
    If (iset(isub)==3 .Or. iset(isub)==4) Then
      q2sf = pmas(23, 1)**2
      If (isub==8 .Or. isub==76 .Or. isub==77) q2sf = pmas(24, 1)**2
    End If
    Do i = 3 - mint(41), mint(42)
      xsf = x(i)
      If (iset(isub)==5) xsf = x(i)/vint(142+i)
      Call pystfu(mint(10+i), xsf, q2sf, xpq, i)
      Do kfl = -6, 6
        xsfx(i, kfl) = xpq(kfl)
      End Do
    End Do
  End If

!...Calculate alpha_strong and K-factor.
  If (mstp(33)/=3) as = ulalps(q2)
  fack = 1.
  faca = 1.
  If (mstp(33)==1) Then
    fack = parp(31)
  Else If (mstp(33)==2) Then
    fack = parp(31)
    faca = parp(32)/parp(31)
  Else If (mstp(33)==3) Then
    q2as = parp(33)*q2
    If (iset(isub)==5 .And. mstp(82)>=2) q2as = q2as + paru(112)*parp(82)
    as = ulalps(q2as)
  End If
  radc = 1. + as/paru(1)

!...Set flags for allowed reacting partons/leptons.
  Do i = 1, 2
    Do j = -40, 40
      kfac(i, j) = 0
    End Do
    If (mint(40+i)==1) Then
      kfac(i, mint(10+i)) = 1
    Else
      Do j = -40, 40
        kfac(i, j) = kfin(i, j)
        If (abs(j)>mstp(54) .And. j/=21) kfac(i, j) = 0
        If (abs(j)<=6) Then
          If (xsfx(i,j)<1.E-10) kfac(i, j) = 0
        Else If (j==21) Then
          If (xsfx(i,0)<1.E-10) kfac(i, 21) = 0
        End If
      End Do
    End If
  End Do

!...Lower and upper limit for flavour loops.
  min1 = 0
  max1 = 0
  min2 = 0
  max2 = 0
  Do j = -20, 20
    If (kfac(1,-j)==1) min1 = -j
    If (kfac(1,j)==1) max1 = j
    If (kfac(2,-j)==1) min2 = -j
    If (kfac(2,j)==1) max2 = j
  End Do
  mina = min(min1, min2)
  maxa = max(max1, max2)

!...Common conversion factors (including Jacobian) for subprocesses.
  sqmz = pmas(23, 1)**2
  gmmz = pmas(23, 1)*pmas(23, 2)
  sqmw = pmas(24, 1)**2
  gmmw = pmas(24, 1)*pmas(24, 2)
  sqmh = pmas(25, 1)**2
  gmmh = pmas(25, 1)*pmas(25, 2)
  sqmzp = pmas(32, 1)**2
  gmmzp = pmas(32, 1)*pmas(32, 2)
  sqmhc = pmas(37, 1)**2
  gmmhc = pmas(37, 1)*pmas(37, 2)
  sqmr = pmas(40, 1)**2
  gmmr = pmas(40, 1)*pmas(40, 2)
  aem = paru(101)
  xw = paru(102)

!...Phase space integral in tau and y*.
  comfac = paru(1)*paru(5)/vint(2)
  If (mint(43)==4) comfac = comfac*fack
  If ((mint(43)>=2 .Or. iset(isub)==3 .Or. iset(isub)==4) .And. iset(isub)/=5) Then
    atau0 = log(taumax/taumin)
    atau1 = (taumax-taumin)/(taumax*taumin)
    h1 = coef(isub, 1) + (atau0/atau1)*coef(isub, 2)/tau
    If (mint(72)>=1) Then
      taur1 = vint(73)
      gamr1 = vint(74)
      atau2 = log(taumax/taumin*(taumin+taur1)/(taumax+taur1))/taur1
      atau3 = (atan((taumax-taur1)/gamr1)-atan((taumin-taur1)/gamr1))/gamr1
      h1 = h1 + (atau0/atau2)*coef(isub, 3)/(tau+taur1) + (atau0/atau3)*coef(isub, 4)*tau/((tau-taur1)**2+gamr1**2)
    End If
    If (mint(72)==2) Then
      taur2 = vint(75)
      gamr2 = vint(76)
      atau4 = log(taumax/taumin*(taumin+taur2)/(taumax+taur2))/taur2
      atau5 = (atan((taumax-taur2)/gamr2)-atan((taumin-taur2)/gamr2))/gamr2
      h1 = h1 + (atau0/atau4)*coef(isub, 5)/(tau+taur2) + (atau0/atau5)*coef(isub, 6)*tau/((tau-taur2)**2+gamr2**2)
    End If
    comfac = comfac*atau0/(tau*h1)
  End If
  If (mint(43)==4 .And. iset(isub)/=5) Then
    ayst0 = ystmax - ystmin
    ayst1 = 0.5*(ystmax-ystmin)**2
    ayst2 = ayst1
    ayst3 = 2.*(atan(exp(ystmax))-atan(exp(ystmin)))
    h2 = (ayst0/ayst1)*coef(isub, 7)*(yst-ystmin) + (ayst0/ayst2)*coef(isub, 8)*(ystmax-yst) + (ayst0/ayst3)*coef(isub, 9)/cosh(yst)
    comfac = comfac*ayst0/h2
  End If

!...2 -> 1 processes: reduction in angular part of phase space integral
!...for case of decaying resonance.
  acth0 = ctnmax - ctnmin + ctpmax - ctpmin
!lin-4/2008 modified a la pythia6115.f to avoid invalid MDCY subcript#1,
!     also break up compound IF statements:
!      IF((ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3).AND.
!     &MDCY(KFPR(ISUB,1),1).EQ.1) THEN
!        IF(KFPR(ISUB,1).EQ.25.OR.KFPR(ISUB,1).EQ.37) THEN
!          COMFAC=COMFAC*0.5*ACTH0
!        ELSE
!          COMFAC=COMFAC*0.125*(3.*ACTH0+CTNMAX**3-CTNMIN**3+
!     &    CTPMAX**3-CTPMIN**3)
!        ENDIF
  If (iset(isub)==1 .Or. iset(isub)==3) Then
    If (mdcy(lucomp(kfpr(isub,1)),1)==1) Then
      If (kfpr(isub,1)==25 .Or. kfpr(isub,1)==37) Then
        comfac = comfac*0.5*acth0
      Else
        comfac = comfac*0.125*(3.*acth0+ctnmax**3-ctnmin**3+ctpmax**3-ctpmin**3)
      End If
    End If
!
!...2 -> 2 processes: angular part of phase space integral.
  Else If (iset(isub)==2 .Or. iset(isub)==4) Then
    acth1 = log((max(rm34,rsqm-ctnmin)*max(rm34,rsqm-ctpmin))/(max(rm34,rsqm-ctnmax)*max(rm34,rsqm-ctpmax)))
    acth2 = log((max(rm34,rsqm+ctnmax)*max(rm34,rsqm+ctpmax))/(max(rm34,rsqm+ctnmin)*max(rm34,rsqm+ctpmin)))
    acth3 = 1./max(rm34, rsqm-ctnmax) - 1./max(rm34, rsqm-ctnmin) + 1./max(rm34, rsqm-ctpmax) - 1./max(rm34, rsqm-ctpmin)
    acth4 = 1./max(rm34, rsqm+ctnmin) - 1./max(rm34, rsqm+ctnmax) + 1./max(rm34, rsqm+ctpmin) - 1./max(rm34, rsqm+ctpmax)
    h3 = coef(isub, 10) + (acth0/acth1)*coef(isub, 11)/max(rm34, rsqm-cth) + (acth0/acth2)*coef(isub, 12)/max(rm34, rsqm+cth) + (acth0/acth3)*coef(isub, 13)/max(rm34, rsqm-cth)**2 + (acth0/acth4)*coef(isub, 14)/max(rm34, rsqm+cth)**2
    comfac = comfac*acth0*0.5*be34/h3
  End If

!...2 -> 3, 4 processes: phace space integral in tau'.
  If (mint(43)>=2 .And. (iset(isub)==3 .Or. iset(isub)==4)) Then
    ataup0 = log(taupmx/taupmn)
    ataup1 = ((1.-tau/taupmx)**4-(1.-tau/taupmn)**4)/(4.*tau)
    h4 = coef(isub, 15) + ataup0/ataup1*coef(isub, 16)/taup*(1.-tau/taup)**3
    If (1.-tau/taup>1.E-4) Then
      fzw = (1.+tau/taup)*log(taup/tau) - 2.*(1.-tau/taup)
    Else
      fzw = 1./6.*(1.-tau/taup)**3*tau/taup
    End If
    comfac = comfac*ataup0*fzw/h4
  End If

!...Phase space integral for low-pT and multiple interactions.
  If (iset(isub)==5) Then
    comfac = paru(1)*paru(5)*fack*0.5*vint(2)/sh2
    atau0 = log(2.*(1.+sqrt(1.-xt2))/xt2-1.)
    atau1 = 2.*atan(1./xt2-1.)/sqrt(xt2)
    h1 = coef(isub, 1) + (atau0/atau1)*coef(isub, 2)/sqrt(tau)
    comfac = comfac*atau0/h1
    ayst0 = ystmax - ystmin
    ayst1 = 0.5*(ystmax-ystmin)**2
    ayst3 = 2.*(atan(exp(ystmax))-atan(exp(ystmin)))
    h2 = (ayst0/ayst1)*coef(isub, 7)*(yst-ystmin) + (ayst0/ayst1)*coef(isub, 8)*(ystmax-yst) + (ayst0/ayst3)*coef(isub, 9)/cosh(yst)
    comfac = comfac*ayst0/h2
    If (mstp(82)<=1) comfac = comfac*xt2**2*(1./vint(149)-1.)
!...For MSTP(82)>=2 an additional factor (xT2/(xT2+VINT(149))**2 is
!...introduced to make cross-section finite for xT2 -> 0.
    If (mstp(82)>=2) comfac = comfac*xt2**2/(vint(149)*(1.+vint(149)))
  End If

!...A: 2 -> 1, tree diagrams.

  145 If (isub<=10) Then
    If (isub==1) Then
!...f + fb -> gamma*/Z0.
      mint(61) = 2
      Call pywidt(23, sqrt(sh), wdtp, wdte)
      facz = comfac*aem**2*4./3.
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 150
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        facf = 1.
        If (iabs(i)<=10) facf = faca/3.
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facf*facz*(ei**2*vint(111)+ei*vi/(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)*vint(112)+(vi**2+ai**2)/(16.*xw*(1.-xw))**2*sh2/((sh-sqmz)**2+gmmz**2)*vint(114))
      150 End Do

    Else If (isub==2) Then
!...f + fb' -> W+/-.
      Call pywidt(24, sqrt(sh), wdtp, wdte)
      facw = comfac*(aem/xw)**2*1./24*sh2/((sh-sqmw)**2+gmmw**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 170
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 160
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 160
          If ((ia<=10 .And. ja>10) .Or. (ia>10 .And. ja<=10)) Goto 160
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          facf = 1.
          If (ia<=10) facf = vckm((ia+1)/2, (ja+1)/2)*faca/3.
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facf*facw*(wdte(0,1)+wdte(0,(5-kchw)/2)+wdte(0,4))
        160 End Do
      170 End Do

    Else If (isub==3) Then
!...f + fb -> H0.
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      fach = comfac*(aem/xw)**2*1./48.*(sh/sqmw)**2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 180
        rmq = pmas(iabs(i), 1)**2/sh
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = fach*rmq*sqrt(max(0.,1.-4.*rmq))
      180 End Do

    Else If (isub==4) Then
!...gamma + W+/- -> W+/-.

    Else If (isub==5) Then
!...Z0 + Z0 -> H0.
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      fach = comfac*1./(128.*paru(1)**2*16.*(1.-xw)**3)*(aem/xw)**4*(sh/sqmw)**2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 200
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 190
          ei = kchg(iabs(i), 1)/3.
          ai = sign(1., ei)
          vi = ai - 4.*ei*xw
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*(vi**2+ai**2)*(vj**2+aj**2)
        190 End Do
      200 End Do

    Else If (isub==6) Then
!...Z0 + W+/- -> W+/-.

    Else If (isub==7) Then
!...W+ + W- -> Z0.

    Else If (isub==8) Then
!...W+ + W- -> H0.
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      fach = comfac*1./(128*paru(1)**2)*(aem/xw)**4*(sh/sqmw)**2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 220
        ei = sign(1., float(i))*kchg(iabs(i), 1)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 210
          ej = sign(1., float(j))*kchg(iabs(j), 1)
          If (ei*ej>0.) Goto 210
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*vint(180+i)*vint(180+j)
        210 End Do
      220 End Do
    End If

!...B: 2 -> 2, tree diagrams.

  Else If (isub<=20) Then
    If (isub==11) Then
!...f + f' -> f + f'.
      facqq1 = comfac*as**2*4./9.*(sh2+uh2)/th2
      facqqb = comfac*as**2*4./9.*((sh2+uh2)/th2*faca-mstp(34)*2./3.*uh2/(sh*th))
      facqq2 = comfac*as**2*4./9.*((sh2+th2)/uh2-mstp(34)*2./3.*sh2/(th*uh))
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 240
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 230
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facqq1
          If (i==-j) sigh(nchn) = facqqb
          If (i==j) Then
            sigh(nchn) = 0.5*sigh(nchn)
            nchn = nchn + 1
            isig(nchn, 1) = i
            isig(nchn, 2) = j
            isig(nchn, 3) = 2
            sigh(nchn) = 0.5*facqq2
          End If
        230 End Do
      240 End Do

    Else If (isub==12) Then
!...f + fb -> f' + fb' (q + qb -> q' + qb' only).
      Call pywidt(21, sqrt(sh), wdtp, wdte)
      facqqb = comfac*as**2*4./9.*(th2+uh2)/sh2*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 250
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facqqb
      250 End Do

    Else If (isub==13) Then
!...f + fb -> g + g (q + qb -> g + g only).
      facgg1 = comfac*as**2*32./27.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)
      facgg2 = comfac*as**2*32./27.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 260
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = 0.5*facgg1
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 2
        sigh(nchn) = 0.5*facgg2
      260 End Do

    Else If (isub==14) Then
!...f + fb -> g + gamma (q + qb -> g + gamma only).
      facgg = comfac*as*aem*8./9.*(th2+uh2)/(th*uh)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 270
        ei = kchg(iabs(i), 1)/3.
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgg*ei**2
      270 End Do

    Else If (isub==15) Then
!...f + fb -> g + Z0 (q + qb -> g + Z0 only).
      faczg = comfac*as*aem/(xw*(1.-xw))*1./18.*(th2+uh2+2.*sqm4*sh)/(th*uh)
      faczg = faczg*wids(23, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 280
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = faczg*(vi**2+ai**2)
      280 End Do

    Else If (isub==16) Then
!...f + fb' -> g + W+/- (q + qb' -> g + W+/- only).
      facwg = comfac*as*aem/xw*2./9.*(th2+uh2+2.*sqm4*sh)/(th*uh)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 300
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 290
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 290
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facwg*fckm*wids(24, (5-kchw)/2)
        290 End Do
      300 End Do

    Else If (isub==17) Then
!...f + fb -> g + H0 (q + qb -> g + H0 only).

    Else If (isub==18) Then
!...f + fb -> gamma + gamma.
      facgg = comfac*faca*aem**2*1./3.*(th2+uh2)/(th*uh)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 310
        ei = kchg(iabs(i), 1)/3.
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgg*ei**4
      310 End Do

    Else If (isub==19) Then
!...f + fb -> gamma + Z0.
      facgz = comfac*faca*aem**2/(xw*(1.-xw))*1./24.*(th2+uh2+2.*sqm4*sh)/(th*uh)
      facgz = facgz*wids(23, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 320
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgz*ei**2*(vi**2+ai**2)
      320 End Do

    Else If (isub==20) Then
!...f + fb' -> gamma + W+/-.
      facgw = comfac*faca*aem**2/xw*1./6.*((2.*uh-th)/(3.*(sh-sqm4)))**2*(th2+uh2+2.*sqm4*sh)/(th*uh)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 340
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 330
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 330
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facgw*fckm*wids(24, (5-kchw)/2)
        330 End Do
      340 End Do
    End If

  Else If (isub<=30) Then
    If (isub==21) Then
!...f + fb -> gamma + H0.

    Else If (isub==22) Then
!...f + fb -> Z0 + Z0.
      faczz = comfac*faca*(aem/(xw*(1.-xw)))**2*1./768.*(uh/th+th/uh+2.*(sqm3+sqm4)*sh/(th*uh)-sqm3*sqm4*(1./th2+1./uh2))
      faczz = faczz*wids(23, 1)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 350
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = faczz*(vi**4+6.*vi**2*ai**2+ai**4)
      350 End Do

    Else If (isub==23) Then
!...f + fb' -> Z0 + W+/-.
      faczw = comfac*faca*(aem/xw)**2*1./6.
      faczw = faczw*wids(23, 2)
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 370
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 360
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 360
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          ei = kchg(ia, 1)/3.
          ai = sign(1., ei)
          vi = ai - 4.*ei*xw
          ej = kchg(ja, 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          If (vi+ai>0) Then
            visav = vi
            aisav = ai
            vi = vj
            ai = aj
            vj = visav
            aj = aisav
          End If
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = faczw*fckm*(1./(sh-sqmw)**2*((9.-8.*xw)/4.*thuh+(8.*xw-6.)/4.*sh*(sqm3+sqm4))+(thuh-sh*(sqm3+sqm4))/(2.*(sh-sqmw))*((vj+aj)/th-(vi+ai)/uh)+thuh/(16.*(1.-xw))*((vj+aj)**2/th2+(vi+ai)**2/uh2)+sh*(sqm3+sqm4)/(8.*(1.-xw))*(vi+ai)*(vj+aj)/(th*uh))*wids(24, (5-kchw)/2)
        360 End Do
      370 End Do

    Else If (isub==24) Then
!...f + fb -> Z0 + H0.
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      fachz = comfac*faca*(aem/(xw*(1.-xw)))**2*1./96.*(thuh+2.*sh*sqmz)/(sh-sqmz)**2
      fachz = fachz*wids(23, 2)*wids(25, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 380
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = fachz*(vi**2+ai**2)
      380 End Do

    Else If (isub==25) Then
!...f + fb -> W+ + W-.
      facww = comfac*faca*(aem/xw)**2*1./12.
      facww = facww*wids(24, 1)
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 390
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        dsigww = thuh/sh2*(3.-(sh-3.*(sqm3+sqm4))/(sh-sqmz)*(vi+ai)/(2.*ai*(1.-xw))+(sh/(sh-sqmz))**2*(1.-2.*(sqm3+sqm4)/sh+12.*sqm3*sqm4/sh2)*(vi**2+ai**2)/(8.*(1.-xw)**2)) - 2.*sqmz/(sh-sqmz)*(vi+ai)/ai + sqmz*sh/(sh-sqmz)**2*(1.-2.*(sqm3+sqm4)/sh)*(vi**2+ai**2)/(2.*(1.-xw))
        If (kchg(iabs(i),1)<0) Then
          dsigww = dsigww + 2.*(1.+sqmz/(sh-sqmz)*(vi+ai)/(2.*ai))*(thuh/(sh*th)-(sqm3+sqm4)/th) + thuh/th2
        Else
          dsigww = dsigww + 2.*(1.+sqmz/(sh-sqmz)*(vi+ai)/(2.*ai))*(thuh/(sh*uh)-(sqm3+sqm4)/uh) + thuh/uh2
        End If
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facww*dsigww
      390 End Do

    Else If (isub==26) Then
!...f + fb' -> W+/- + H0.
      thuh = max(th*uh-sqm3*sqm4, sh*ckin(3)**2)
      fachw = comfac*faca*(aem/xw)**2*1./24.*(thuh+2.*sh*sqmw)/(sh-sqmw)**2
      fachw = fachw*wids(25, 2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 410
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(1,j)==0) Goto 400
          ja = iabs(j)
          If (i*j>0 .Or. mod(ia+ja,2)==0) Goto 400
          kchw = (kchg(ia,1)*isign(1,i)+kchg(ja,1)*isign(1,j))/3
          fckm = 1.
          If (mint(43)==4) fckm = vckm((ia+1)/2, (ja+1)/2)
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fachw*fckm*wids(24, (5-kchw)/2)
        400 End Do
      410 End Do

    Else If (isub==27) Then
!...f + fb -> H0 + H0.

    Else If (isub==28) Then
!...f + g -> f + g (q + g -> q + g only).
      facqg1 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*uh2/th2-uh/sh)*faca
      facqg2 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*sh2/th2-sh/uh)
      Do i = mina, maxa
        If (i==0) Goto 430
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 420
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 420
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facqg1
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 2
          sigh(nchn) = facqg2
        420 End Do
      430 End Do

    Else If (isub==29) Then
!...f + g -> f + gamma (q + g -> q + gamma only).
      fgq = comfac*faca*as*aem*1./3.*(sh2+uh2)/(-sh*uh)
      Do i = mina, maxa
        If (i==0) Goto 450
        ei = kchg(iabs(i), 1)/3.
        facgq = fgq*ei**2
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 440
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 440
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facgq
        440 End Do
      450 End Do

    Else If (isub==30) Then
!...f + g -> f + Z0 (q + g -> q + Z0 only).
      fzq = comfac*faca*as*aem/(xw*(1.-xw))*1./48.*(sh2+uh2+2.*sqm4*th)/(-sh*uh)
      fzq = fzq*wids(23, 2)
      Do i = mina, maxa
        If (i==0) Goto 470
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        faczq = fzq*(vi**2+ai**2)
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 460
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 460
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = faczq
        460 End Do
      470 End Do
    End If

  Else If (isub<=40) Then
    If (isub==31) Then
!...f + g -> f' + W+/- (q + g -> q' + W+/- only).
      facwq = comfac*faca*as*aem/xw*1./12.*(sh2+uh2+2.*sqm4*th)/(-sh*uh)
      Do i = mina, maxa
        If (i==0) Goto 490
        ia = iabs(i)
        kchw = isign(1, kchg(ia,1)*isign(1,i))
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 480
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 480
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facwq*vint(180+i)*wids(24, (5-kchw)/2)
        480 End Do
      490 End Do

    Else If (isub==32) Then
!...f + g -> f + H0 (q + g -> q + H0 only).

    Else If (isub==33) Then
!...f + gamma -> f + g (q + gamma -> q + g only).

    Else If (isub==34) Then
!...f + gamma -> f + gamma.

    Else If (isub==35) Then
!...f + gamma -> f + Z0.

    Else If (isub==36) Then
!...f + gamma -> f' + W+/-.

    Else If (isub==37) Then
!...f + gamma -> f + H0.

    Else If (isub==38) Then
!...f + Z0 -> f + g (q + Z0 -> q + g only).

    Else If (isub==39) Then
!...f + Z0 -> f + gamma.

    Else If (isub==40) Then
!...f + Z0 -> f + Z0.
    End If

  Else If (isub<=50) Then
    If (isub==41) Then
!...f + Z0 -> f' + W+/-.

    Else If (isub==42) Then
!...f + Z0 -> f + H0.

    Else If (isub==43) Then
!...f + W+/- -> f' + g (q + W+/- -> q' + g only).

    Else If (isub==44) Then
!...f + W+/- -> f' + gamma.

    Else If (isub==45) Then
!...f + W+/- -> f' + Z0.

    Else If (isub==46) Then
!...f + W+/- -> f' + W+/-.

    Else If (isub==47) Then
!...f + W+/- -> f' + H0.

    Else If (isub==48) Then
!...f + H0 -> f + g (q + H0 -> q + g only).

    Else If (isub==49) Then
!...f + H0 -> f + gamma.

    Else If (isub==50) Then
!...f + H0 -> f + Z0.
    End If

  Else If (isub<=60) Then
    If (isub==51) Then
!...f + H0 -> f' + W+/-.

    Else If (isub==52) Then
!...f + H0 -> f + H0.

    Else If (isub==53) Then
!...g + g -> f + fb (g + g -> q + qb only).
      Call pywidt(21, sqrt(sh), wdtp, wdte)
      facqq1 = comfac*as**2*1./6.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      facqq2 = comfac*as**2*1./6.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      If (kfac(1,21)*kfac(2,21)==0) Goto 500
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facqq1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 2
      sigh(nchn) = facqq2
      500 Continue

    Else If (isub==54) Then
!...g + gamma -> f + fb (g + gamma -> q + qb only).

    Else If (isub==55) Then
!...g + gamma -> f + fb (g + gamma -> q + qb only).

    Else If (isub==56) Then
!...g + gamma -> f + fb (g + gamma -> q + qb only).

    Else If (isub==57) Then
!...g + gamma -> f + fb (g + gamma -> q + qb only).

    Else If (isub==58) Then
!...gamma + gamma -> f + fb.

    Else If (isub==59) Then
!...gamma + Z0 -> f + fb.

    Else If (isub==60) Then
!...gamma + W+/- -> f + fb'.
    End If

  Else If (isub<=70) Then
    If (isub==61) Then
!...gamma + H0 -> f + fb.

    Else If (isub==62) Then
!...Z0 + Z0 -> f + fb.

    Else If (isub==63) Then
!...Z0 + W+/- -> f + fb'.

    Else If (isub==64) Then
!...Z0 + H0 -> f + fb.

    Else If (isub==65) Then
!...W+ + W- -> f + fb.

    Else If (isub==66) Then
!...W+/- + H0 -> f + fb'.

    Else If (isub==67) Then
!...H0 + H0 -> f + fb.

    Else If (isub==68) Then
!...g + g -> g + g.
      facgg1 = comfac*as**2*9./4.*(sh2/th2+2.*sh/th+3.+2.*th/sh+th2/sh2)*faca
      facgg2 = comfac*as**2*9./4.*(uh2/sh2+2.*uh/sh+3.+2.*sh/uh+sh2/uh2)*faca
      facgg3 = comfac*as**2*9./4.*(th2/uh2+2.*th/uh+3+2.*uh/th+uh2/th2)
      If (kfac(1,21)*kfac(2,21)==0) Goto 510
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = 0.5*facgg1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 2
      sigh(nchn) = 0.5*facgg2
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 3
      sigh(nchn) = 0.5*facgg3
      510 Continue

    Else If (isub==69) Then
!...gamma + gamma -> W+ + W-.

    Else If (isub==70) Then
!...gamma + W+/- -> gamma + W+/-.
    End If

  Else If (isub<=80) Then
    If (isub==71) Then
!...Z0 + Z0 -> Z0 + Z0.
      be2 = 1. - 4.*sqmz/sh
      th = -0.5*sh*be2*(1.-cth)
      uh = -0.5*sh*be2*(1.+cth)
      shang = 1./(1.-xw)*sqmw/sqmz*(1.+be2)**2
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      thang = 1./(1.-xw)*sqmw/sqmz*(be2-cth)**2
      athre = (th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
      athim = -gmmh/((th-sqmh)**2+gmmh**2)*thang
      uhang = 1./(1.-xw)*sqmw/sqmz*(be2+cth)**2
      auhre = (uh-sqmh)/((uh-sqmh)**2+gmmh**2)*uhang
      auhim = -gmmh/((uh-sqmh)**2+gmmh**2)*uhang
      fach = 0.5*comfac*1./(4096.*paru(1)**2*16.*(1.-xw)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+athre+auhre)**2+(ashim+athim+auhim)**2)*sqmz/sqmw
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 530
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        avi = ai**2 + vi**2
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 520
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          avj = aj**2 + vj**2
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*avi*avj
        520 End Do
      530 End Do

    Else If (isub==72) Then
!...Z0 + Z0 -> W+ + W-.
      be2 = sqrt((1.-4.*sqmw/sh)*(1.-4.*sqmz/sh))
      cth2 = cth**2
      th = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh-be2*cth)
      uh = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh+be2*cth)
      shang = 4.*sqrt(sqmw/(sqmz*(1.-xw)))*(1.-2.*sqmw/sh)*(1.-2.*sqmz/sh)
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      atwre = (1.-xw)/sqmz*sh/(th-sqmw)*((cth-be2)**2*(3./2.+be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2+2.*(sqmw+sqmz)/sh*be2*cth))
      atwim = 0.
      auwre = (1.-xw)/sqmz*sh/(uh-sqmw)*((cth+be2)**2*(3./2.-be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2-2.*(sqmw+sqmz)/sh*be2*cth))
      auwim = 0.
      a4re = 2.*(1.-xw)/sqmz*(3.-cth2-4.*(sqmw+sqmz)/sh)
      a4im = 0.
      fach = comfac*1./(4096.*paru(1)**2*16.*(1.-xw)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+atwre+auwre+a4re)**2+(ashim+atwim+auwim+a4im)**2)*sqmz/sqmw
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 550
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        avi = ai**2 + vi**2
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 540
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = aj - 4.*ej*xw
          avj = aj**2 + vj**2
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*avi*avj
        540 End Do
      550 End Do

    Else If (isub==73) Then
!...Z0 + W+/- -> Z0 + W+/-.
      be2 = 1. - 2.*(sqmz+sqmw)/sh + ((sqmz-sqmw)/sh)**2
      ep1 = 1. + (sqmz-sqmw)/sh
      ep2 = 1. - (sqmz-sqmw)/sh
      th = -0.5*sh*be2*(1.-cth)
      uh = (sqmz-sqmw)**2/sh - 0.5*sh*be2*(1.+cth)
      thang = sqrt(sqmw/(sqmz*(1.-xw)))*(be2-ep1*cth)*(be2-ep2*cth)
      athre = (th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
      athim = -gmmh/((th-sqmh)**2+gmmh**2)*thang
      aswre = (1.-xw)/sqmz*sh/(sh-sqmw)*(-be2*(ep1+ep2)**4*cth+1./4.*(be2+ep1*ep2)**2*((ep1-ep2)**2-4.*be2*cth)+2.*be2*(be2+ep1*ep2)*(ep1+ep2)**2*cth-1./16.*sh/sqmw*(ep1**2-ep2**2)**2*(be2+ep1*ep2)**2)
      aswim = 0.
      auwre = (1.-xw)/sqmz*sh/(uh-sqmw)*(-be2*(ep2+ep1*cth)*(ep1+ep2*cth)*(be2+ep1*ep2)+be2*(ep2+ep1*cth)*(be2+ep1*ep2*cth)*(2.*ep2-ep2*cth+ep1)-be2*(ep2+ep1*cth)**2*(be2-ep2**2*cth)-1./8.*(be2+ep1*ep2*cth)**2*((ep1+ep2)**2+2.*be2*(1.-cth))+1./32.*sh/sqmw*(be2+ep1*ep2*cth)**2*(ep1**2-ep2**2)**2-be2*(ep1+ep2*cth)*(ep2+ep1*cth)*(be2+ep1*ep2)+be2*(ep1+ep2*cth)*(be2+ep1*ep2*cth)*(2.*ep1-ep1*cth+ep2)-be2*(ep1+ep2*cth)**2*(be2-ep1**2*cth)-1./8.*(be2+ep1*ep2*cth)**2*((ep1+ep2)**2+2.*be2*(1.-cth))+1./32.*sh/sqmw*(be2+ep1*ep2*cth)**2*(ep1**2-ep2**2)**2)
      auwim = 0.
      a4re = (1.-xw)/sqmz*(ep1**2*ep2**2*(cth**2-1.)-2.*be2*(ep1**2+ep2**2+ep1*ep2)*cth-2.*be2*ep1*ep2)
      a4im = 0.
      fach = comfac*1./(4096.*paru(1)**2*4.*(1.-xw))*(aem/xw)**4*(sh/sqmw)**2*((athre+aswre+auwre+a4re)**2+(athim+aswim+auwim+a4im)**2)*sqrt(sqmz/sqmw)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 570
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        avi = ai**2 + vi**2
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 560
          ej = kchg(iabs(j), 1)/3.
          aj = sign(1., ej)
          vj = ai - 4.*ej*xw
          avj = aj**2 + vj**2
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*(avi*vint(180+j)+vint(180+i)*avj)
        560 End Do
      570 End Do

    Else If (isub==75) Then
!...W+ + W- -> gamma + gamma.

    Else If (isub==76) Then
!...W+ + W- -> Z0 + Z0.
      be2 = sqrt((1.-4.*sqmw/sh)*(1.-4.*sqmz/sh))
      cth2 = cth**2
      th = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh-be2*cth)
      uh = -0.5*sh*(1.-2.*(sqmw+sqmz)/sh+be2*cth)
      shang = 4.*sqrt(sqmw/(sqmz*(1.-xw)))*(1.-2.*sqmw/sh)*(1.-2.*sqmz/sh)
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      atwre = (1.-xw)/sqmz*sh/(th-sqmw)*((cth-be2)**2*(3./2.+be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2+2.*(sqmw+sqmz)/sh*be2*cth))
      atwim = 0.
      auwre = (1.-xw)/sqmz*sh/(uh-sqmw)*((cth+be2)**2*(3./2.-be2/2.*cth-(sqmw+sqmz)/sh+(sqmw-sqmz)**2/(sh*sqmw))+4.*((sqmw+sqmz)/sh*(1.-3.*cth2)+8.*sqmw*sqmz/sh2*(2.*cth2-1.)+4.*(sqmw**2+sqmz**2)/sh2*cth2-2.*(sqmw+sqmz)/sh*be2*cth))
      auwim = 0.
      a4re = 2.*(1.-xw)/sqmz*(3.-cth2-4.*(sqmw+sqmz)/sh)
      a4im = 0.
      fach = 0.5*comfac*1./(4096.*paru(1)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+atwre+auwre+a4re)**2+(ashim+atwim+auwim+a4im)**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 590
        ei = sign(1., float(i))*kchg(iabs(i), 1)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 580
          ej = sign(1., float(j))*kchg(iabs(j), 1)
          If (ei*ej>0.) Goto 580
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*vint(180+i)*vint(180+j)
        580 End Do
      590 End Do

    Else If (isub==77) Then
!...W+/- + W+/- -> W+/- + W+/-.
      be2 = 1. - 4.*sqmw/sh
      be4 = be2**2
      cth2 = cth**2
      cth3 = cth**3
      th = -0.5*sh*be2*(1.-cth)
      uh = -0.5*sh*be2*(1.+cth)
      shang = (1.+be2)**2
      ashre = (sh-sqmh)/((sh-sqmh)**2+gmmh**2)*shang
      ashim = -gmmh/((sh-sqmh)**2+gmmh**2)*shang
      thang = (be2-cth)**2
      athre = (th-sqmh)/((th-sqmh)**2+gmmh**2)*thang
      athim = -gmmh/((th-sqmh)**2+gmmh**2)*thang
      sgzang = 1./sqmw*be2*(3.-be2)**2*cth
      asgre = xw*sgzang
      asgim = 0.
      aszre = (1.-xw)*sh/(sh-sqmz)*sgzang
      aszim = 0.
      tgzang = 1./sqmw*(be2*(4.-2.*be2+be4)+be2*(4.-10.*be2+be4)*cth+(2.-11.*be2+10.*be4)*cth2+be2*cth3)
      atgre = 0.5*xw*sh/th*tgzang
      atgim = 0.
      atzre = 0.5*(1.-xw)*sh/(th-sqmz)*tgzang
      atzim = 0.
      a4re = 1./sqmw*(1.+2.*be2-6.*be2*cth-cth2)
      a4im = 0.
      fach = comfac*1./(4096.*paru(1)**2)*(aem/xw)**4*(sh/sqmw)**2*((ashre+athre+asgre+aszre+atgre+atzre+a4re)**2+(ashim+athim+asgim+aszim+atgim+atzim+a4im)**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 610
        ei = sign(1., float(i))*kchg(iabs(i), 1)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 600
          ej = sign(1., float(j))*kchg(iabs(j), 1)
          If (ei*ej>0.) Goto 600
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = fach*vint(180+i)*vint(180+j)
        600 End Do
      610 End Do

    Else If (isub==78) Then
!...W+/- + H0 -> W+/- + H0.

    Else If (isub==79) Then
!...H0 + H0 -> H0 + H0.

    End If

!...C: 2 -> 2, tree diagrams with masses.

  Else If (isub<=90) Then
    If (isub==81) Then
!...q + qb -> Q + QB.
      facqqb = comfac*as**2*4./9.*(((th-sqm3)**2+(uh-sqm3)**2)/sh2+2.*sqm3/sh)
      If (mstp(35)>=1) Then
        If (mstp(35)==1) Then
          alssg = parp(35)
        Else
          mst115 = mstu(115)
          mstu(115) = mstp(36)
          q2bn = sqrt(sqm3*((sqrt(sh)-2.*sqrt(sqm3))**2+parp(36)**2))
          alssg = ulalps(q2bn)
          mstu(115) = mst115
        End If
        xrepu = paru(1)*alssg/(6.*sqrt(max(1E-20,1.-4.*sqm3/sh)))
        frepu = xrepu/(exp(min(100.,xrepu))-1.)
        pari(81) = frepu
        facqqb = facqqb*frepu
      End If
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 620
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facqqb
      620 End Do

    Else If (isub==82) Then
!...g + g -> Q + QB.
      facqq1 = comfac*faca*as**2*1./6.*((uh-sqm3)/(th-sqm3)-2.*(uh-sqm3)**2/sh2+4.*sqm3/sh*(th*uh-sqm3**2)/(th-sqm3)**2)
      facqq2 = comfac*faca*as**2*1./6.*((th-sqm3)/(uh-sqm3)-2.*(th-sqm3)**2/sh2+4.*sqm3/sh*(th*uh-sqm3**2)/(uh-sqm3)**2)
      If (mstp(35)>=1) Then
        If (mstp(35)==1) Then
          alssg = parp(35)
        Else
          mst115 = mstu(115)
          mstu(115) = mstp(36)
          q2bn = sqrt(sqm3*((sqrt(sh)-2.*sqrt(sqm3))**2+parp(36)**2))
          alssg = ulalps(q2bn)
          mstu(115) = mst115
        End If
        xattr = 4.*paru(1)*alssg/(3.*sqrt(max(1E-20,1.-4.*sqm3/sh)))
        fattr = xattr/(1.-exp(-min(100.,xattr)))
        xrepu = paru(1)*alssg/(6.*sqrt(max(1E-20,1.-4.*sqm3/sh)))
        frepu = xrepu/(exp(min(100.,xrepu))-1.)
        fatre = (2.*fattr+5.*frepu)/7.
        pari(81) = fatre
        facqq1 = facqq1*fatre
        facqq2 = facqq2*fatre
      End If
      If (kfac(1,21)*kfac(2,21)==0) Goto 630
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facqq1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 2
      sigh(nchn) = facqq2
      630 Continue

    End If

!...D: Mimimum bias processes.

  Else If (isub<=100) Then
    If (isub==91) Then
!...Elastic scattering.
      sigs = xsec(isub, 1)

    Else If (isub==92) Then
!...Single diffractive scattering.
      sigs = xsec(isub, 1)

    Else If (isub==93) Then
!...Double diffractive scattering.
      sigs = xsec(isub, 1)

    Else If (isub==94) Then
!...Central diffractive scattering.
      sigs = xsec(isub, 1)

    Else If (isub==95) Then
!...Low-pT scattering.
      sigs = xsec(isub, 1)

    Else If (isub==96) Then
!...Multiple interactions: sum of QCD processes.
      Call pywidt(21, sqrt(sh), wdtp, wdte)

!...q + q' -> q + q'.
      facqq1 = comfac*as**2*4./9.*(sh2+uh2)/th2
      facqqb = comfac*as**2*4./9.*((sh2+uh2)/th2*faca-mstp(34)*2./3.*uh2/(sh*th))
      facqq2 = comfac*as**2*4./9.*((sh2+th2)/uh2-mstp(34)*2./3.*sh2/(th*uh))
      Do i = -3, 3
        If (i==0) Goto 650
        Do j = -3, 3
          If (j==0) Goto 640
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 111
          sigh(nchn) = facqq1
          If (i==-j) sigh(nchn) = facqqb
          If (i==j) Then
            sigh(nchn) = 0.5*sigh(nchn)
            nchn = nchn + 1
            isig(nchn, 1) = i
            isig(nchn, 2) = j
            isig(nchn, 3) = 112
            sigh(nchn) = 0.5*facqq2
          End If
        640 End Do
      650 End Do

!...q + qb -> q' + qb' or g + g.
      facqqb = comfac*as**2*4./9.*(th2+uh2)/sh2*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))
      facgg1 = comfac*as**2*32./27.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)
      facgg2 = comfac*as**2*32./27.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)
      Do i = -3, 3
        If (i==0) Goto 660
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 121
        sigh(nchn) = facqqb
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 131
        sigh(nchn) = 0.5*facgg1
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 132
        sigh(nchn) = 0.5*facgg2
      660 End Do

!...q + g -> q + g.
      facqg1 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*uh2/th2-uh/sh)*faca
      facqg2 = comfac*as**2*4./9.*((2.+mstp(34)*1./4.)*sh2/th2-sh/uh)
      Do i = -3, 3
        If (i==0) Goto 680
        Do isde = 1, 2
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 281
          sigh(nchn) = facqg1
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 282
          sigh(nchn) = facqg2
        End Do
      680 End Do

!...g + g -> q + qb or g + g.
      facqq1 = comfac*as**2*1./6.*(uh/th-(2.+mstp(34)*1./4.)*uh2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      facqq2 = comfac*as**2*1./6.*(th/uh-(2.+mstp(34)*1./4.)*th2/sh2)*(wdte(0,1)+wdte(0,2)+wdte(0,3)+wdte(0,4))*faca
      facgg1 = comfac*as**2*9./4.*(sh2/th2+2.*sh/th+3.+2.*th/sh+th2/sh2)*faca
      facgg2 = comfac*as**2*9./4.*(uh2/sh2+2.*uh/sh+3.+2.*sh/uh+sh2/uh2)*faca
      facgg3 = comfac*as**2*9./4.*(th2/uh2+2.*th/uh+3+2.*uh/th+uh2/th2)
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 531
      sigh(nchn) = facqq1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 532
      sigh(nchn) = facqq2
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 681
      sigh(nchn) = 0.5*facgg1
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 682
      sigh(nchn) = 0.5*facgg2
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 683
      sigh(nchn) = 0.5*facgg3
    End If

!...E: 2 -> 1, loop diagrams.

  Else If (isub<=110) Then
    If (isub==101) Then
!...g + g -> gamma*/Z0.

    Else If (isub==102) Then
!...g + g -> H0.
      Call pywidt(25, sqrt(sh), wdtp, wdte)
      etare = 0.
      etaim = 0.
      Do i = 1, 2*mstp(1)
        eps = 4.*pmas(i, 1)**2/sh
        If (eps<=1.) Then
          If (eps>1.E-4) Then
            root = sqrt(1.-eps)
            rln = log((1.+root)/(1.-root))
          Else
            rln = log(4./eps-2.)
          End If
          phire = 0.25*(rln**2-paru(1)**2)
          phiim = 0.5*paru(1)*rln
        Else
          phire = -(asin(1./sqrt(eps)))**2
          phiim = 0.
        End If
        etare = etare + 0.5*eps*(1.+(eps-1.)*phire)
        etaim = etaim + 0.5*eps*(eps-1.)*phiim
      End Do
      eta2 = etare**2 + etaim**2
      fach = comfac*faca*(as/paru(1)*aem/xw)**2*1./512.*(sh/sqmw)**2*eta2*sh2/((sh-sqmh)**2+gmmh**2)*(wdte(0,1)+wdte(0,2)+wdte(0,4))
      If (kfac(1,21)*kfac(2,21)==0) Goto 700
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = fach
      700 Continue

    End If

!...F: 2 -> 2, box diagrams.

  Else If (isub<=120) Then
    If (isub==111) Then
!...f + fb -> g + H0 (q + qb -> g + H0 only).
      a5stur = 0.
      a5stui = 0.
      Do i = 1, 2*mstp(1)
        sqmq = pmas(i, 1)**2
        epss = 4.*sqmq/sh
        epsh = 4.*sqmq/sqmh
        a5stur = a5stur + sqmq/sqmh*(4.+4.*sh/(th+uh)*(pyw1au(epss,1)-pyw1au(epsh,1))+(1.-4.*sqmq/(th+uh))*(pyw2au(epss,1)-pyw2au(epsh,1)))
        a5stui = a5stui + sqmq/sqmh*(4.*sh/(th+uh)*(pyw1au(epss,2)-pyw1au(epsh,2))+(1.-4.*sqmq/(th+uh))*(pyw2au(epss,2)-pyw2au(epsh,2)))
      End Do
      facgh = comfac*faca/(144.*paru(1)**2)*aem/xw*as**3*sqmh/sqmw*sqmh/sh*(uh**2+th**2)/(uh+th)**2*(a5stur**2+a5stui**2)
      facgh = facgh*wids(25, 2)
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 720
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = facgh
      720 End Do

    Else If (isub==112) Then
!...f + g -> f + H0 (q + g -> q + H0 only).
      a5tsur = 0.
      a5tsui = 0.
      Do i = 1, 2*mstp(1)
        sqmq = pmas(i, 1)**2
        epst = 4.*sqmq/th
        epsh = 4.*sqmq/sqmh
        a5tsur = a5tsur + sqmq/sqmh*(4.+4.*th/(sh+uh)*(pyw1au(epst,1)-pyw1au(epsh,1))+(1.-4.*sqmq/(sh+uh))*(pyw2au(epst,1)-pyw2au(epsh,1)))
        a5tsui = a5tsui + sqmq/sqmh*(4.*th/(sh+uh)*(pyw1au(epst,2)-pyw1au(epsh,2))+(1.-4.*sqmq/(sh+uh))*(pyw2au(epst,2)-pyw2au(epsh,2)))
      End Do
      facqh = comfac*faca/(384.*paru(1)**2)*aem/xw*as**3*sqmh/sqmw*sqmh/(-th)*(uh**2+sh**2)/(uh+sh)**2*(a5tsur**2+a5tsui**2)
      facqh = facqh*wids(25, 2)
      Do i = mina, maxa
        If (i==0) Goto 750
        Do isde = 1, 2
          If (isde==1 .And. kfac(1,i)*kfac(2,21)==0) Goto 740
          If (isde==2 .And. kfac(1,21)*kfac(2,i)==0) Goto 740
          nchn = nchn + 1
          isig(nchn, isde) = i
          isig(nchn, 3-isde) = 21
          isig(nchn, 3) = 1
          sigh(nchn) = facqh
        740 End Do
      750 End Do

    Else If (isub==113) Then
!...g + g -> g + H0.
      a2stur = 0.
      a2stui = 0.
      a2ustr = 0.
      a2usti = 0.
      a2tusr = 0.
      a2tusi = 0.
      a4stur = 0.
      a4stui = 0.
      Do i = 6, 2*mstp(1)
!'''Only t-quarks yet included
        sqmq = pmas(i, 1)**2
        epss = 4.*sqmq/sh
        epst = 4.*sqmq/th
        epsu = 4.*sqmq/uh
        epsh = 4.*sqmq/sqmh
        If (epsh<1.E-6) Goto 760
        bestu = 0.5*(1.+sqrt(1.+epss*th/uh))
        beust = 0.5*(1.+sqrt(1.+epsu*sh/th))
        betus = 0.5*(1.+sqrt(1.+epst*uh/sh))
        beuts = bestu
        betsu = beust
        besut = betus
        w3stur = pyi3au(bestu, epsh, 1) - pyi3au(bestu, epss, 1) - pyi3au(bestu, epsu, 1)
        w3stui = pyi3au(bestu, epsh, 2) - pyi3au(bestu, epss, 2) - pyi3au(bestu, epsu, 2)
        w3sutr = pyi3au(besut, epsh, 1) - pyi3au(besut, epss, 1) - pyi3au(besut, epst, 1)
        w3suti = pyi3au(besut, epsh, 2) - pyi3au(besut, epss, 2) - pyi3au(besut, epst, 2)
        w3tsur = pyi3au(betsu, epsh, 1) - pyi3au(betsu, epst, 1) - pyi3au(betsu, epsu, 1)
        w3tsui = pyi3au(betsu, epsh, 2) - pyi3au(betsu, epst, 2) - pyi3au(betsu, epsu, 2)
        w3tusr = pyi3au(betus, epsh, 1) - pyi3au(betus, epst, 1) - pyi3au(betus, epss, 1)
        w3tusi = pyi3au(betus, epsh, 2) - pyi3au(betus, epst, 2) - pyi3au(betus, epss, 2)
        w3ustr = pyi3au(beust, epsh, 1) - pyi3au(beust, epsu, 1) - pyi3au(beust, epst, 1)
        w3usti = pyi3au(beust, epsh, 2) - pyi3au(beust, epsu, 2) - pyi3au(beust, epst, 2)
        w3utsr = pyi3au(beuts, epsh, 1) - pyi3au(beuts, epsu, 1) - pyi3au(beuts, epss, 1)
        w3utsi = pyi3au(beuts, epsh, 2) - pyi3au(beuts, epsu, 2) - pyi3au(beuts, epss, 2)
        b2stur = sqmq/sqmh**2*(sh*(uh-sh)/(sh+uh)+2.*th*uh*(uh+2.*sh)/(sh+uh)**2*(pyw1au(epst,1)-pyw1au(epsh,1))+(sqmq-sh/4.)*(0.5*pyw2au(epss,1)+0.5*pyw2au(epsh,1)-pyw2au(epst,1)+w3stur)+sh**2*(2.*sqmq/(sh+uh)**2-0.5/(sh+uh))*(pyw2au(epst,1)-pyw2au(epsh,1))+0.5*th*uh/sh*(pyw2au(epsh,1)-2.*pyw2au(epst,1))+0.125*(sh-12.*sqmq-4.*th*uh/sh)*w3tsur)
        b2stui = sqmq/sqmh**2*(2.*th*uh*(uh+2.*sh)/(sh+uh)**2*(pyw1au(epst,2)-pyw1au(epsh,2))+(sqmq-sh/4.)*(0.5*pyw2au(epss,2)+0.5*pyw2au(epsh,2)-pyw2au(epst,2)+w3stui)+sh**2*(2.*sqmq/(sh+uh)**2-0.5/(sh+uh))*(pyw2au(epst,2)-pyw2au(epsh,2))+0.5*th*uh/sh*(pyw2au(epsh,2)-2.*pyw2au(epst,2))+0.125*(sh-12.*sqmq-4.*th*uh/sh)*w3tsui)
        b2sutr = sqmq/sqmh**2*(sh*(th-sh)/(sh+th)+2.*uh*th*(th+2.*sh)/(sh+th)**2*(pyw1au(epsu,1)-pyw1au(epsh,1))+(sqmq-sh/4.)*(0.5*pyw2au(epss,1)+0.5*pyw2au(epsh,1)-pyw2au(epsu,1)+w3sutr)+sh**2*(2.*sqmq/(sh+th)**2-0.5/(sh+th))*(pyw2au(epsu,1)-pyw2au(epsh,1))+0.5*uh*th/sh*(pyw2au(epsh,1)-2.*pyw2au(epsu,1))+0.125*(sh-12.*sqmq-4.*uh*th/sh)*w3ustr)
        b2suti = sqmq/sqmh**2*(2.*uh*th*(th+2.*sh)/(sh+th)**2*(pyw1au(epsu,2)-pyw1au(epsh,2))+(sqmq-sh/4.)*(0.5*pyw2au(epss,2)+0.5*pyw2au(epsh,2)-pyw2au(epsu,2)+w3suti)+sh**2*(2.*sqmq/(sh+th)**2-0.5/(sh+th))*(pyw2au(epsu,2)-pyw2au(epsh,2))+0.5*uh*th/sh*(pyw2au(epsh,2)-2.*pyw2au(epsu,2))+0.125*(sh-12.*sqmq-4.*uh*th/sh)*w3usti)
        b2tsur = sqmq/sqmh**2*(th*(uh-th)/(th+uh)+2.*sh*uh*(uh+2.*th)/(th+uh)**2*(pyw1au(epss,1)-pyw1au(epsh,1))+(sqmq-th/4.)*(0.5*pyw2au(epst,1)+0.5*pyw2au(epsh,1)-pyw2au(epss,1)+w3tsur)+th**2*(2.*sqmq/(th+uh)**2-0.5/(th+uh))*(pyw2au(epss,1)-pyw2au(epsh,1))+0.5*sh*uh/th*(pyw2au(epsh,1)-2.*pyw2au(epss,1))+0.125*(th-12.*sqmq-4.*sh*uh/th)*w3stur)
        b2tsui = sqmq/sqmh**2*(2.*sh*uh*(uh+2.*th)/(th+uh)**2*(pyw1au(epss,2)-pyw1au(epsh,2))+(sqmq-th/4.)*(0.5*pyw2au(epst,2)+0.5*pyw2au(epsh,2)-pyw2au(epss,2)+w3tsui)+th**2*(2.*sqmq/(th+uh)**2-0.5/(th+uh))*(pyw2au(epss,2)-pyw2au(epsh,2))+0.5*sh*uh/th*(pyw2au(epsh,2)-2.*pyw2au(epss,2))+0.125*(th-12.*sqmq-4.*sh*uh/th)*w3stui)
        b2tusr = sqmq/sqmh**2*(th*(sh-th)/(th+sh)+2.*uh*sh*(sh+2.*th)/(th+sh)**2*(pyw1au(epsu,1)-pyw1au(epsh,1))+(sqmq-th/4.)*(0.5*pyw2au(epst,1)+0.5*pyw2au(epsh,1)-pyw2au(epsu,1)+w3tusr)+th**2*(2.*sqmq/(th+sh)**2-0.5/(th+sh))*(pyw2au(epsu,1)-pyw2au(epsh,1))+0.5*uh*sh/th*(pyw2au(epsh,1)-2.*pyw2au(epsu,1))+0.125*(th-12.*sqmq-4.*uh*sh/th)*w3utsr)
        b2tusi = sqmq/sqmh**2*(2.*uh*sh*(sh+2.*th)/(th+sh)**2*(pyw1au(epsu,2)-pyw1au(epsh,2))+(sqmq-th/4.)*(0.5*pyw2au(epst,2)+0.5*pyw2au(epsh,2)-pyw2au(epsu,2)+w3tusi)+th**2*(2.*sqmq/(th+sh)**2-0.5/(th+sh))*(pyw2au(epsu,2)-pyw2au(epsh,2))+0.5*uh*sh/th*(pyw2au(epsh,2)-2.*pyw2au(epsu,2))+0.125*(th-12.*sqmq-4.*uh*sh/th)*w3utsi)
        b2ustr = sqmq/sqmh**2*(uh*(th-uh)/(uh+th)+2.*sh*th*(th+2.*uh)/(uh+th)**2*(pyw1au(epss,1)-pyw1au(epsh,1))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,1)+0.5*pyw2au(epsh,1)-pyw2au(epss,1)+w3ustr)+uh**2*(2.*sqmq/(uh+th)**2-0.5/(uh+th))*(pyw2au(epss,1)-pyw2au(epsh,1))+0.5*sh*th/uh*(pyw2au(epsh,1)-2.*pyw2au(epss,1))+0.125*(uh-12.*sqmq-4.*sh*th/uh)*w3sutr)
        b2usti = sqmq/sqmh**2*(2.*sh*th*(th+2.*uh)/(uh+th)**2*(pyw1au(epss,2)-pyw1au(epsh,2))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,2)+0.5*pyw2au(epsh,2)-pyw2au(epss,2)+w3usti)+uh**2*(2.*sqmq/(uh+th)**2-0.5/(uh+th))*(pyw2au(epss,2)-pyw2au(epsh,2))+0.5*sh*th/uh*(pyw2au(epsh,2)-2.*pyw2au(epss,2))+0.125*(uh-12.*sqmq-4.*sh*th/uh)*w3suti)
        b2utsr = sqmq/sqmh**2*(uh*(sh-uh)/(uh+sh)+2.*th*sh*(sh+2.*uh)/(uh+sh)**2*(pyw1au(epst,1)-pyw1au(epsh,1))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,1)+0.5*pyw2au(epsh,1)-pyw2au(epst,1)+w3utsr)+uh**2*(2.*sqmq/(uh+sh)**2-0.5/(uh+sh))*(pyw2au(epst,1)-pyw2au(epsh,1))+0.5*th*sh/uh*(pyw2au(epsh,1)-2.*pyw2au(epst,1))+0.125*(uh-12.*sqmq-4.*th*sh/uh)*w3tusr)
        b2utsi = sqmq/sqmh**2*(2.*th*sh*(sh+2.*uh)/(uh+sh)**2*(pyw1au(epst,2)-pyw1au(epsh,2))+(sqmq-uh/4.)*(0.5*pyw2au(epsu,2)+0.5*pyw2au(epsh,2)-pyw2au(epst,2)+w3utsi)+uh**2*(2.*sqmq/(uh+sh)**2-0.5/(uh+sh))*(pyw2au(epst,2)-pyw2au(epsh,2))+0.5*th*sh/uh*(pyw2au(epsh,2)-2.*pyw2au(epst,2))+0.125*(uh-12.*sqmq-4.*th*sh/uh)*w3tusi)
        b4stur = sqmq/sqmh*(-2./3.+(sqmq/sqmh-1./4.)*(pyw2au(epss,1)-pyw2au(epsh,1)+w3stur))
        b4stui = sqmq/sqmh*(sqmq/sqmh-1./4.)*(pyw2au(epss,2)-pyw2au(epsh,2)+w3stui)
        b4tusr = sqmq/sqmh*(-2./3.+(sqmq/sqmh-1./4.)*(pyw2au(epst,1)-pyw2au(epsh,1)+w3tusr))
        b4tusi = sqmq/sqmh*(sqmq/sqmh-1./4.)*(pyw2au(epst,2)-pyw2au(epsh,2)+w3tusi)
        b4ustr = sqmq/sqmh*(-2./3.+(sqmq/sqmh-1./4.)*(pyw2au(epsu,1)-pyw2au(epsh,1)+w3ustr))
        b4usti = sqmq/sqmh*(sqmq/sqmh-1./4.)*(pyw2au(epsu,2)-pyw2au(epsh,2)+w3usti)
        a2stur = a2stur + b2stur + b2sutr
        a2stui = a2stui + b2stui + b2suti
        a2ustr = a2ustr + b2ustr + b2utsr
        a2usti = a2usti + b2usti + b2utsi
        a2tusr = a2tusr + b2tusr + b2tsur
        a2tusi = a2tusi + b2tusi + b2tsui
        a4stur = a4stur + b4stur + b4ustr + b4tusr
        a4stui = a4stui + b4stui + b4usti + b4tusi
      760 End Do
      facgh = comfac*faca*3./(128.*paru(1)**2)*aem/xw*as**3*sqmh/sqmw*sqmh**3/(sh*th*uh)*(a2stur**2+a2stui**2+a2ustr**2+a2usti**2+a2tusr**2+a2tusi**2+a4stur**2+a4stui**2)
      facgh = facgh*wids(25, 2)
      If (kfac(1,21)*kfac(2,21)==0) Goto 770
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facgh
      770 Continue

    Else If (isub==114) Then
!...g + g -> gamma + gamma.
      asre = 0.
      asim = 0.
      Do i = 1, 2*mstp(1)
        ei = kchg(iabs(i), 1)/3.
        sqmq = pmas(i, 1)**2
        epss = 4.*sqmq/sh
        epst = 4.*sqmq/th
        epsu = 4.*sqmq/uh
        If (epss+abs(epst)+abs(epsu)<3.E-6) Then
          a0stur = 1. + (th-uh)/sh*log(th/uh) + 0.5*(th2+uh2)/sh2*(log(th/uh)**2+paru(1)**2)
          a0stui = 0.
          a0tsur = 1. + (sh-uh)/th*log(-sh/uh) + 0.5*(sh2+uh2)/th2*log(-sh/uh)**2
          a0tsui = -paru(1)*((sh-uh)/th+(sh2+uh2)/th2*log(-sh/uh))
          a0utsr = 1. + (th-sh)/uh*log(-th/sh) + 0.5*(th2+sh2)/uh2*log(-th/sh)**2
          a0utsi = paru(1)*((th-sh)/uh+(th2+sh2)/uh2*log(-th/sh))
          a1stur = -1.
          a1stui = 0.
          a2stur = -1.
          a2stui = 0.
        Else
          bestu = 0.5*(1.+sqrt(1.+epss*th/uh))
          beust = 0.5*(1.+sqrt(1.+epsu*sh/th))
          betus = 0.5*(1.+sqrt(1.+epst*uh/sh))
          beuts = bestu
          betsu = beust
          besut = betus
          a0stur = 1. + (1.+2.*th/sh)*pyw1au(epst, 1) + (1.+2.*uh/sh)*pyw1au(epsu, 1) + 0.5*((th2+uh2)/sh2-epss)*(pyw2au(epst,1)+pyw2au(epsu,1)) - 0.25*epst*(1.-0.5*epss)*(pyi3au(besut,epss,1)+pyi3au(besut,epst,1)) - 0.25*epsu*(1.-0.5*epss)*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1)) + 0.25*(-2.*(th2+uh2)/sh2+4.*epss+epst+epsu+0.5*epst*epsu)*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1))
          a0stui = (1.+2.*th/sh)*pyw1au(epst, 2) + (1.+2.*uh/sh)*pyw1au(epsu, 2) + 0.5*((th2+uh2)/sh2-epss)*(pyw2au(epst,2)+pyw2au(epsu,2)) - 0.25*epst*(1.-0.5*epss)*(pyi3au(besut,epss,2)+pyi3au(besut,epst,2)) - 0.25*epsu*(1.-0.5*epss)*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2)) + 0.25*(-2.*(th2+uh2)/sh2+4.*epss+epst+epsu+0.5*epst*epsu)*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2))
          a0tsur = 1. + (1.+2.*sh/th)*pyw1au(epss, 1) + (1.+2.*uh/th)*pyw1au(epsu, 1) + 0.5*((sh2+uh2)/th2-epst)*(pyw2au(epss,1)+pyw2au(epsu,1)) - 0.25*epss*(1.-0.5*epst)*(pyi3au(betus,epst,1)+pyi3au(betus,epss,1)) - 0.25*epsu*(1.-0.5*epst)*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1)) + 0.25*(-2.*(sh2+uh2)/th2+4.*epst+epss+epsu+0.5*epss*epsu)*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1))
          a0tsui = (1.+2.*sh/th)*pyw1au(epss, 2) + (1.+2.*uh/th)*pyw1au(epsu, 2) + 0.5*((sh2+uh2)/th2-epst)*(pyw2au(epss,2)+pyw2au(epsu,2)) - 0.25*epss*(1.-0.5*epst)*(pyi3au(betus,epst,2)+pyi3au(betus,epss,2)) - 0.25*epsu*(1.-0.5*epst)*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2)) + 0.25*(-2.*(sh2+uh2)/th2+4.*epst+epss+epsu+0.5*epss*epsu)*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2))
          a0utsr = 1. + (1.+2.*th/uh)*pyw1au(epst, 1) + (1.+2.*sh/uh)*pyw1au(epss, 1) + 0.5*((th2+sh2)/uh2-epsu)*(pyw2au(epst,1)+pyw2au(epss,1)) - 0.25*epst*(1.-0.5*epsu)*(pyi3au(beust,epsu,1)+pyi3au(beust,epst,1)) - 0.25*epss*(1.-0.5*epsu)*(pyi3au(beuts,epsu,1)+pyi3au(beuts,epss,1)) + 0.25*(-2.*(th2+sh2)/uh2+4.*epsu+epst+epss+0.5*epst*epss)*(pyi3au(betus,epst,1)+pyi3au(betus,epss,1))
          a0utsi = (1.+2.*th/uh)*pyw1au(epst, 2) + (1.+2.*sh/uh)*pyw1au(epss, 2) + 0.5*((th2+sh2)/uh2-epsu)*(pyw2au(epst,2)+pyw2au(epss,2)) - 0.25*epst*(1.-0.5*epsu)*(pyi3au(beust,epsu,2)+pyi3au(beust,epst,2)) - 0.25*epss*(1.-0.5*epsu)*(pyi3au(beuts,epsu,2)+pyi3au(beuts,epss,2)) + 0.25*(-2.*(th2+sh2)/uh2+4.*epsu+epst+epss+0.5*epst*epss)*(pyi3au(betus,epst,2)+pyi3au(betus,epss,2))
          a1stur = -1. - 0.25*(epss+epst+epsu)*(pyw2au(epss,1)+pyw2au(epst,1)+pyw2au(epsu,1)) + 0.25*(epsu+0.5*epss*epst)*(pyi3au(besut,epss,1)+pyi3au(besut,epst,1)) + 0.25*(epst+0.5*epss*epsu)*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1)) + 0.25*(epss+0.5*epst*epsu)*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1))
          a1stui = -0.25*(epss+epst+epsu)*(pyw2au(epss,2)+pyw2au(epst,2)+pyw2au(epsu,2)) + 0.25*(epsu+0.5*epss*epst)*(pyi3au(besut,epss,2)+pyi3au(besut,epst,2)) + 0.25*(epst+0.5*epss*epsu)*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2)) + 0.25*(epss+0.5*epst*epsu)*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2))
          a2stur = -1. + 0.125*epss*epst*(pyi3au(besut,epss,1)+pyi3au(besut,epst,1)) + 0.125*epss*epsu*(pyi3au(bestu,epss,1)+pyi3au(bestu,epsu,1)) + 0.125*epst*epsu*(pyi3au(betsu,epst,1)+pyi3au(betsu,epsu,1))
          a2stui = 0.125*epss*epst*(pyi3au(besut,epss,2)+pyi3au(besut,epst,2)) + 0.125*epss*epsu*(pyi3au(bestu,epss,2)+pyi3au(bestu,epsu,2)) + 0.125*epst*epsu*(pyi3au(betsu,epst,2)+pyi3au(betsu,epsu,2))
        End If
        asre = asre + ei**2*(a0stur+a0tsur+a0utsr+4.*a1stur+a2stur)
        asim = asim + ei**2*(a0stui+a0tsui+a0utsi+4.*a1stui+a2stui)
      End Do
      facgg = comfac*faca/(8.*paru(1)**2)*as**2*aem**2*(asre**2+asim**2)
      If (kfac(1,21)*kfac(2,21)==0) Goto 790
      nchn = nchn + 1
      isig(nchn, 1) = 21
      isig(nchn, 2) = 21
      isig(nchn, 3) = 1
      sigh(nchn) = facgg
      790 Continue

    Else If (isub==115) Then
!...g + g -> gamma + Z0.

    Else If (isub==116) Then
!...g + g -> Z0 + Z0.

    Else If (isub==117) Then
!...g + g -> W+ + W-.

    End If

!...G: 2 -> 3, tree diagrams.

  Else If (isub<=140) Then
    If (isub==121) Then
!...g + g -> f + fb + H0.

    End If

!...H: 2 -> 1, tree diagrams, non-standard model processes.

  Else If (isub<=160) Then
    If (isub==141) Then
!...f + fb -> gamma*/Z0/Z'0.
      mint(61) = 2
      Call pywidt(32, sqrt(sh), wdtp, wdte)
      faczp = comfac*aem**2*4./9.
      Do i = mina, maxa
        If (i==0 .Or. kfac(1,i)*kfac(2,-i)==0) Goto 800
        ei = kchg(iabs(i), 1)/3.
        ai = sign(1., ei)
        vi = ai - 4.*ei*xw
        api = sign(1., ei)
        vpi = api - 4.*ei*xw
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = faczp*(ei**2*vint(111)+ei*vi/(8.*xw*(1.-xw))*sh*(sh-sqmz)/((sh-sqmz)**2+gmmz**2)*vint(112)+ei*vpi/(8.*xw*(1.-xw))*sh*(sh-sqmzp)/((sh-sqmzp)**2+gmmzp**2)*vint(113)+(vi**2+ai**2)/(16.*xw*(1.-xw))**2*sh2/((sh-sqmz)**2+gmmz**2)*vint(114)+2.*(vi*vpi+ai*api)/(16.*xw*(1.-xw))**2*sh2*((sh-sqmz)*(sh-sqmzp)+gmmz*gmmzp)/(((sh-sqmz)**2+gmmz**2)*((sh-sqmzp)**2+gmmzp**2))*vint(115)+(vpi**2+api**2)/(16.*xw*(1.-xw))**2*sh2/((sh-sqmzp)**2+gmmzp**2)*vint(116))
      800 End Do

    Else If (isub==142) Then
!...f + fb' -> H+/-.
      Call pywidt(37, sqrt(sh), wdtp, wdte)
      fhc = comfac*(aem/xw)**2*1./48.*(sh/sqmw)**2*sh2/((sh-sqmhc)**2+gmmhc**2)
!'''No construction yet for leptons
      Do i = 1, mstp(54)/2
        il = 2*i - 1
        iu = 2*i
        rmql = pmas(il, 1)**2/sh
        rmqu = pmas(iu, 1)**2/sh
        fachc = fhc*((rmql*paru(121)+rmqu/paru(121))*(1.-rmql-rmqu)-4.*rmql*rmqu)/sqrt(max(0.,(1.-rmql-rmqu)**2-4.*rmql*rmqu))
        If (kfac(1,il)*kfac(2,-iu)==0) Goto 810
        kchhc = (kchg(il,1)-kchg(iu,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = il
        isig(nchn, 2) = -iu
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        810 If (kfac(1,-il)*kfac(2,iu)==0) Goto 820
        kchhc = (-kchg(il,1)+kchg(iu,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = -il
        isig(nchn, 2) = iu
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        820 If (kfac(1,iu)*kfac(2,-il)==0) Goto 830
        kchhc = (kchg(iu,1)-kchg(il,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = iu
        isig(nchn, 2) = -il
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        830 If (kfac(1,-iu)*kfac(2,il)==0) Goto 840
        kchhc = (-kchg(iu,1)+kchg(il,1))/3
        nchn = nchn + 1
        isig(nchn, 1) = -iu
        isig(nchn, 2) = il
        isig(nchn, 3) = 1
        sigh(nchn) = fachc*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
      840 End Do

    Else If (isub==143) Then
!...f + fb -> R.
      Call pywidt(40, sqrt(sh), wdtp, wdte)
      facr = comfac*(aem/xw)**2*1./9.*sh2/((sh-sqmr)**2+gmmr**2)
      Do i = min1, max1
        If (i==0 .Or. kfac(1,i)==0) Goto 860
        ia = iabs(i)
        Do j = min2, max2
          If (j==0 .Or. kfac(2,j)==0) Goto 850
          ja = iabs(j)
          If (i*j>0 .Or. iabs(ia-ja)/=2) Goto 850
          nchn = nchn + 1
          isig(nchn, 1) = i
          isig(nchn, 2) = j
          isig(nchn, 3) = 1
          sigh(nchn) = facr*(wdte(0,1)+wdte(0,(10-(i+j))/4)+wdte(0,4))
        850 End Do
      860 End Do

    End If

!...I: 2 -> 2, tree diagrams, non-standard model processes.

  Else
    If (isub==161) Then
!...f + g -> f' + H+/- (q + g -> q' + H+/- only).
      fhcq = comfac*faca*as*aem/xw*1./24
      Do i = 1, mstp(54)
        iu = i + mod(i, 2)
        sqmq = pmas(iu, 1)**2
        fachcq = fhcq/paru(121)*sqmq/sqmw*(sh/(sqmq-uh)+2.*sqmq*(sqmhc-uh)/(sqmq-uh)**2+(sqmq-uh)/sh+2.*sqmq/(sqmq-uh)+2.*(sqmhc-uh)/(sqmq-uh)*(sqmhc-sqmq-sh)/sh)
        If (kfac(1,-i)*kfac(2,21)==0) Goto 870
        kchhc = isign(1, -kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = -i
        isig(nchn, 2) = 21
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        870 If (kfac(1,i)*kfac(2,21)==0) Goto 880
        kchhc = isign(1, kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = i
        isig(nchn, 2) = 21
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        880 If (kfac(1,21)*kfac(2,-i)==0) Goto 890
        kchhc = isign(1, -kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = 21
        isig(nchn, 2) = -i
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
        890 If (kfac(1,21)*kfac(2,i)==0) Goto 900
        kchhc = isign(1, kchg(i,1))
        nchn = nchn + 1
        isig(nchn, 1) = 21
        isig(nchn, 2) = i
        isig(nchn, 3) = 1
        sigh(nchn) = fachcq*(wdte(0,1)+wdte(0,(5-kchhc)/2)+wdte(0,4))
      900 End Do

    End If
  End If

!...Multiply with structure functions.
  If (isub<=90 .Or. isub>=96) Then
    Do ichn = 1, nchn
      If (mint(41)==2) Then
        kfl1 = isig(ichn, 1)
        If (kfl1==21) kfl1 = 0
        sigh(ichn) = sigh(ichn)*xsfx(1, kfl1)
      End If
      If (mint(42)==2) Then
        kfl2 = isig(ichn, 2)
        If (kfl2==21) kfl2 = 0
        sigh(ichn) = sigh(ichn)*xsfx(2, kfl2)
      End If
      sigs = sigs + sigh(ichn)
    End Do
  End If

  Return
End Subroutine pysigh

!*********************************************************************


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn





Subroutine pystfu(kf, x, q2, xpq, jbt)

!                        *******JBT specifies beam or target of the particle
!...Gives proton and pi+ parton structure functions according to a few
!...different parametrizations. Note that what is coded is x times the
!...probability distribution, i.e. xq(x,Q2) etc.
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Save /hparnt/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Save /hjcrdn/
!                        ********COMMON BLOCK FROM HIJING
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Dimension xpq(-6:6), xq(6), tx(6), tt(6), ts(6), nehlq(8, 2), cehlq(6, 6, 2, 8, 2), cdo(3, 6, 5, 2), cow(3, 5, 4, 2)

!...The following data lines are coefficients needed in the
!...Eichten, Hinchliffe, Lane, Quigg proton structure function
!...parametrizations, see below.
!...Powers of 1-x in different cases.
  Data nehlq/3, 4, 7, 5, 7, 7, 7, 7, 3, 4, 7, 6, 7, 7, 7, 7/
!...Expansion coefficients for up valence quark distribution.
  Data (((cehlq(ix,it,nx,1,1),ix=1,6),it=1,6), nx=1, 2)/7.677E-01, -2.087E-01, -3.303E-01, -2.517E-02, -1.570E-02, -1.000E-04, -5.326E-01, -2.661E-01, 3.201E-01, 1.192E-01, 2.434E-02, 7.620E-03, 2.162E-01, 1.881E-01, -8.375E-02, -6.515E-02, -1.743E-02, -5.040E-03, -9.211E-02, -9.952E-02, 1.373E-02, 2.506E-02, 8.770E-03, 2.550E-03, 3.670E-02, 4.409E-02, 9.600E-04, -7.960E-03, -3.420E-03, -1.050E-03, -1.549E-02, -2.026E-02, -3.060E-03, 2.220E-03, 1.240E-03, 4.100E-04, 2.395E-01, 2.905E-01, 9.778E-02, 2.149E-02, 3.440E-03, 5.000E-04, 1.751E-02, -6.090E-03, -2.687E-02, -1.916E-02, -7.970E-03, -2.750E-03, -5.760E-03, -5.040E-03, 1.080E-03, 2.490E-03, 1.530E-03, 7.500E-04, 1.740E-03, 1.960E-03, 3.000E-04, -3.400E-04, -2.900E-04, -1.800E-04, -5.300E-04, -6.400E-04, -1.700E-04, 4.000E-05, 6.000E-05, 4.000E-05, 1.700E-04, 2.200E-04, 8.000E-05, 1.000E-05, -1.000E-05, -1.000E-05/
  Data (((cehlq(ix,it,nx,1,2),ix=1,6),it=1,6), nx=1, 2)/7.237E-01, -2.189E-01, -2.995E-01, -1.909E-02, -1.477E-02, 2.500E-04, -5.314E-01, -2.425E-01, 3.283E-01, 1.119E-01, 2.223E-02, 7.070E-03, 2.289E-01, 1.890E-01, -9.859E-02, -6.900E-02, -1.747E-02, -5.080E-03, -1.041E-01, -1.084E-01, 2.108E-02, 2.975E-02, 9.830E-03, 2.830E-03, 4.394E-02, 5.116E-02, -1.410E-03, -1.055E-02, -4.230E-03, -1.270E-03, -1.991E-02, -2.539E-02, -2.780E-03, 3.430E-03, 1.720E-03, 5.500E-04, 2.410E-01, 2.884E-01, 9.369E-02, 1.900E-02, 2.530E-03, 2.400E-04, 1.765E-02, -9.220E-03, -3.037E-02, -2.085E-02, -8.440E-03, -2.810E-03, -6.450E-03, -5.260E-03, 1.720E-03, 3.110E-03, 1.830E-03, 8.700E-04, 2.120E-03, 2.320E-03, 2.600E-04, -4.900E-04, -3.900E-04, -2.300E-04, -6.900E-04, -8.200E-04, -2.000E-04, 7.000E-05, 9.000E-05, 6.000E-05, 2.400E-04, 3.100E-04, 1.100E-04, 0.000E+00, -2.000E-05, -2.000E-05/
!...Expansion coefficients for down valence quark distribution.
  Data (((cehlq(ix,it,nx,2,1),ix=1,6),it=1,6), nx=1, 2)/3.813E-01, -8.090E-02, -1.634E-01, -2.185E-02, -8.430E-03, -6.200E-04, -2.948E-01, -1.435E-01, 1.665E-01, 6.638E-02, 1.473E-02, 4.080E-03, 1.252E-01, 1.042E-01, -4.722E-02, -3.683E-02, -1.038E-02, -2.860E-03, -5.478E-02, -5.678E-02, 8.900E-03, 1.484E-02, 5.340E-03, 1.520E-03, 2.220E-02, 2.567E-02, -3.000E-05, -4.970E-03, -2.160E-03, -6.500E-04, -9.530E-03, -1.204E-02, -1.510E-03, 1.510E-03, 8.300E-04, 2.700E-04, 1.261E-01, 1.354E-01, 3.958E-02, 8.240E-03, 1.660E-03, 4.500E-04, 3.890E-03, -1.159E-02, -1.625E-02, -9.610E-03, -3.710E-03, -1.260E-03, -1.910E-03, -5.600E-04, 1.590E-03, 1.590E-03, 8.400E-04, 3.900E-04, 6.400E-04, 4.900E-04, -1.500E-04, -2.900E-04, -1.800E-04, -1.000E-04, -2.000E-04, -1.900E-04, 0.000E+00, 6.000E-05, 4.000E-05, 3.000E-05, 7.000E-05, 8.000E-05, 2.000E-05, -1.000E-05, -1.000E-05, -1.000E-05/
  Data (((cehlq(ix,it,nx,2,2),ix=1,6),it=1,6), nx=1, 2)/3.578E-01, -8.622E-02, -1.480E-01, -1.840E-02, -7.820E-03, -4.500E-04, -2.925E-01, -1.304E-01, 1.696E-01, 6.243E-02, 1.353E-02, 3.750E-03, 1.318E-01, 1.041E-01, -5.486E-02, -3.872E-02, -1.038E-02, -2.850E-03, -6.162E-02, -6.143E-02, 1.303E-02, 1.740E-02, 5.940E-03, 1.670E-03, 2.643E-02, 2.957E-02, -1.490E-03, -6.450E-03, -2.630E-03, -7.700E-04, -1.218E-02, -1.497E-02, -1.260E-03, 2.240E-03, 1.120E-03, 3.500E-04, 1.263E-01, 1.334E-01, 3.732E-02, 7.070E-03, 1.260E-03, 3.400E-04, 3.660E-03, -1.357E-02, -1.795E-02, -1.031E-02, -3.880E-03, -1.280E-03, -2.100E-03, -3.600E-04, 2.050E-03, 1.920E-03, 9.800E-04, 4.400E-04, 7.700E-04, 5.400E-04, -2.400E-04, -3.900E-04, -2.400E-04, -1.300E-04, -2.600E-04, -2.300E-04, 2.000E-05, 9.000E-05, 6.000E-05, 4.000E-05, 9.000E-05, 1.000E-04, 2.000E-05, -2.000E-05, -2.000E-05, -1.000E-05/
!...Expansion coefficients for up and down sea quark distributions.
  Data (((cehlq(ix,it,nx,3,1),ix=1,6),it=1,6), nx=1, 2)/6.870E-02, -6.861E-02, 2.973E-02, -5.400E-03, 3.780E-03, -9.700E-04, -1.802E-02, 1.400E-04, 6.490E-03, -8.540E-03, 1.220E-03, -1.750E-03, -4.650E-03, 1.480E-03, -5.930E-03, 6.000E-04, -1.030E-03, -8.000E-05, 6.440E-03, 2.570E-03, 2.830E-03, 1.150E-03, 7.100E-04, 3.300E-04, -3.930E-03, -2.540E-03, -1.160E-03, -7.700E-04, -3.600E-04, -1.900E-04, 2.340E-03, 1.930E-03, 5.300E-04, 3.700E-04, 1.600E-04, 9.000E-05, 1.014E+00, -1.106E+00, 3.374E-01, -7.444E-02, 8.850E-03, -8.700E-04, 9.233E-01, -1.285E+00, 4.475E-01, -9.786E-02, 1.419E-02, -1.120E-03, 4.888E-02, -1.271E-01, 8.606E-02, -2.608E-02, 4.780E-03, -6.000E-04, -2.691E-02, 4.887E-02, -1.771E-02, 1.620E-03, 2.500E-04, -6.000E-05, 7.040E-03, -1.113E-02, 1.590E-03, 7.000E-04, -2.000E-04, 0.000E+00, -1.710E-03, 2.290E-03, 3.800E-04, -3.500E-04, 4.000E-05, 1.000E-05/
  Data (((cehlq(ix,it,nx,3,2),ix=1,6),it=1,6), nx=1, 2)/1.008E-01, -7.100E-02, 1.973E-02, -5.710E-03, 2.930E-03, -9.900E-04, -5.271E-02, -1.823E-02, 1.792E-02, -6.580E-03, 1.750E-03, -1.550E-03, 1.220E-02, 1.763E-02, -8.690E-03, -8.800E-04, -1.160E-03, -2.100E-04, -1.190E-03, -7.180E-03, 2.360E-03, 1.890E-03, 7.700E-04, 4.100E-04, -9.100E-04, 2.040E-03, -3.100E-04, -1.050E-03, -4.000E-04, -2.400E-04, 1.190E-03, -1.700E-04, -2.000E-04, 4.200E-04, 1.700E-04, 1.000E-04, 1.081E+00, -1.189E+00, 3.868E-01, -8.617E-02, 1.115E-02, -1.180E-03, 9.917E-01, -1.396E+00, 4.998E-01, -1.159E-01, 1.674E-02, -1.720E-03, 5.099E-02, -1.338E-01, 9.173E-02, -2.885E-02, 5.890E-03, -6.500E-04, -3.178E-02, 5.703E-02, -2.070E-02, 2.440E-03, 1.100E-04, -9.000E-05, 8.970E-03, -1.392E-02, 2.050E-03, 6.500E-04, -2.300E-04, 2.000E-05, -2.340E-03, 3.010E-03, 5.000E-04, -3.900E-04, 6.000E-05, 1.000E-05/
!...Expansion coefficients for gluon distribution.
  Data (((cehlq(ix,it,nx,4,1),ix=1,6),it=1,6), nx=1, 2)/9.482E-01, -9.578E-01, 1.009E-01, -1.051E-01, 3.456E-02, -3.054E-02, -9.627E-01, 5.379E-01, 3.368E-01, -9.525E-02, 1.488E-02, -2.051E-02, 4.300E-01, -8.306E-02, -3.372E-01, 4.902E-02, -9.160E-03, 1.041E-02, -1.925E-01, -1.790E-02, 2.183E-01, 7.490E-03, 4.140E-03, -1.860E-03, 8.183E-02, 1.926E-02, -1.072E-01, -1.944E-02, -2.770E-03, -5.200E-04, -3.884E-02, -1.234E-02, 5.410E-02, 1.879E-02, 3.350E-03, 1.040E-03, 2.948E+01, -3.902E+01, 1.464E+01, -3.335E+00, 5.054E-01, -5.915E-02, 2.559E+01, -3.955E+01, 1.661E+01, -4.299E+00, 6.904E-01, -8.243E-02, -1.663E+00, 1.176E+00, 1.118E+00, -7.099E-01, 1.948E-01, -2.404E-02, -2.168E-01, 8.170E-01, -7.169E-01, 1.851E-01, -1.924E-02, -3.250E-03, 2.088E-01, -4.355E-01, 2.239E-01, -2.446E-02, -3.620E-03, 1.910E-03, -9.097E-02, 1.601E-01, -5.681E-02, -2.500E-03, 2.580E-03, -4.700E-04/
  Data (((cehlq(ix,it,nx,4,2),ix=1,6),it=1,6), nx=1, 2)/2.367E+00, 4.453E-01, 3.660E-01, 9.467E-02, 1.341E-01, 1.661E-02, -3.170E+00, -1.795E+00, 3.313E-02, -2.874E-01, -9.827E-02, -7.119E-02, 1.823E+00, 1.457E+00, -2.465E-01, 3.739E-02, 6.090E-03, 1.814E-02, -1.033E+00, -9.827E-01, 2.136E-01, 1.169E-01, 5.001E-02, 1.684E-02, 5.133E-01, 5.259E-01, -1.173E-01, -1.139E-01, -4.988E-02, -2.021E-02, -2.881E-01, -3.145E-01, 5.667E-02, 9.161E-02, 4.568E-02, 1.951E-02, 3.036E+01, -4.062E+01, 1.578E+01, -3.699E+00, 6.020E-01, -7.031E-02, 2.700E+01, -4.167E+01, 1.770E+01, -4.804E+00, 7.862E-01, -1.060E-01, -1.909E+00, 1.357E+00, 1.127E+00, -7.181E-01, 2.232E-01, -2.481E-02, -2.488E-01, 9.781E-01, -8.127E-01, 2.094E-01, -2.997E-02, -4.710E-03, 2.506E-01, -5.427E-01, 2.672E-01, -3.103E-02, -1.800E-03, 2.870E-03, -1.128E-01, 2.087E-01, -6.972E-02, -2.480E-03, 2.630E-03, -8.400E-04/
!...Expansion coefficients for strange sea quark distribution.
  Data (((cehlq(ix,it,nx,5,1),ix=1,6),it=1,6), nx=1, 2)/4.968E-02, -4.173E-02, 2.102E-02, -3.270E-03, 3.240E-03, -6.700E-04, -6.150E-03, -1.294E-02, 6.740E-03, -6.890E-03, 9.000E-04, -1.510E-03, -8.580E-03, 5.050E-03, -4.900E-03, -1.600E-04, -9.400E-04, -1.500E-04, 7.840E-03, 1.510E-03, 2.220E-03, 1.400E-03, 7.000E-04, 3.500E-04, -4.410E-03, -2.220E-03, -8.900E-04, -8.500E-04, -3.600E-04, -2.000E-04, 2.520E-03, 1.840E-03, 4.100E-04, 3.900E-04, 1.600E-04, 9.000E-05, 9.235E-01, -1.085E+00, 3.464E-01, -7.210E-02, 9.140E-03, -9.100E-04, 9.315E-01, -1.274E+00, 4.512E-01, -9.775E-02, 1.380E-02, -1.310E-03, 4.739E-02, -1.296E-01, 8.482E-02, -2.642E-02, 4.760E-03, -5.700E-04, -2.653E-02, 4.953E-02, -1.735E-02, 1.750E-03, 2.800E-04, -6.000E-05, 6.940E-03, -1.132E-02, 1.480E-03, 6.500E-04, -2.100E-04, 0.000E+00, -1.680E-03, 2.340E-03, 4.200E-04, -3.400E-04, 5.000E-05, 1.000E-05/
  Data (((cehlq(ix,it,nx,5,2),ix=1,6),it=1,6), nx=1, 2)/6.478E-02, -4.537E-02, 1.643E-02, -3.490E-03, 2.710E-03, -6.700E-04, -2.223E-02, -2.126E-02, 1.247E-02, -6.290E-03, 1.120E-03, -1.440E-03, -1.340E-03, 1.362E-02, -6.130E-03, -7.900E-04, -9.000E-04, -2.000E-04, 5.080E-03, -3.610E-03, 1.700E-03, 1.830E-03, 6.800E-04, 4.000E-04, -3.580E-03, 6.000E-05, -2.600E-04, -1.050E-03, -3.800E-04, -2.300E-04, 2.420E-03, 9.300E-04, -1.000E-04, 4.500E-04, 1.700E-04, 1.100E-04, 9.868E-01, -1.171E+00, 3.940E-01, -8.459E-02, 1.124E-02, -1.250E-03, 1.001E+00, -1.383E+00, 5.044E-01, -1.152E-01, 1.658E-02, -1.830E-03, 4.928E-02, -1.368E-01, 9.021E-02, -2.935E-02, 5.800E-03, -6.600E-04, -3.133E-02, 5.785E-02, -2.023E-02, 2.630E-03, 1.600E-04, -8.000E-05, 8.840E-03, -1.416E-02, 1.900E-03, 5.800E-04, -2.500E-04, 1.000E-05, -2.300E-03, 3.080E-03, 5.500E-04, -3.700E-04, 7.000E-05, 1.000E-05/
!...Expansion coefficients for charm sea quark distribution.
  Data (((cehlq(ix,it,nx,6,1),ix=1,6),it=1,6), nx=1, 2)/9.270E-03, -1.817E-02, 9.590E-03, -6.390E-03, 1.690E-03, -1.540E-03, 5.710E-03, -1.188E-02, 6.090E-03, -4.650E-03, 1.240E-03, -1.310E-03, -3.960E-03, 7.100E-03, -3.590E-03, 1.840E-03, -3.900E-04, 3.400E-04, 1.120E-03, -1.960E-03, 1.120E-03, -4.800E-04, 1.000E-04, -4.000E-05, 4.000E-05, -3.000E-05, -1.800E-04, 9.000E-05, -5.000E-05, -2.000E-05, -4.200E-04, 7.300E-04, -1.600E-04, 5.000E-05, 5.000E-05, 5.000E-05, 8.098E-01, -1.042E+00, 3.398E-01, -6.824E-02, 8.760E-03, -9.000E-04, 8.961E-01, -1.217E+00, 4.339E-01, -9.287E-02, 1.304E-02, -1.290E-03, 3.058E-02, -1.040E-01, 7.604E-02, -2.415E-02, 4.600E-03, -5.000E-04, -2.451E-02, 4.432E-02, -1.651E-02, 1.430E-03, 1.200E-04, -1.000E-04, 1.122E-02, -1.457E-02, 2.680E-03, 5.800E-04, -1.200E-04, 3.000E-05, -7.730E-03, 7.330E-03, -7.600E-04, -2.400E-04, 1.000E-05, 0.000E+00/
  Data (((cehlq(ix,it,nx,6,2),ix=1,6),it=1,6), nx=1, 2)/9.980E-03, -1.945E-02, 1.055E-02, -6.870E-03, 1.860E-03, -1.560E-03, 5.700E-03, -1.203E-02, 6.250E-03, -4.860E-03, 1.310E-03, -1.370E-03, -4.490E-03, 7.990E-03, -4.170E-03, 2.050E-03, -4.400E-04, 3.300E-04, 1.470E-03, -2.480E-03, 1.460E-03, -5.700E-04, 1.200E-04, -1.000E-05, -9.000E-05, 1.500E-04, -3.200E-04, 1.200E-04, -6.000E-05, -4.000E-05, -4.200E-04, 7.600E-04, -1.400E-04, 4.000E-05, 7.000E-05, 5.000E-05, 8.698E-01, -1.131E+00, 3.836E-01, -8.111E-02, 1.048E-02, -1.300E-03, 9.626E-01, -1.321E+00, 4.854E-01, -1.091E-01, 1.583E-02, -1.700E-03, 3.057E-02, -1.088E-01, 8.022E-02, -2.676E-02, 5.590E-03, -5.600E-04, -2.845E-02, 5.164E-02, -1.918E-02, 2.210E-03, -4.000E-05, -1.500E-04, 1.311E-02, -1.751E-02, 3.310E-03, 5.100E-04, -1.200E-04, 5.000E-05, -8.590E-03, 8.380E-03, -9.200E-04, -2.600E-04, 1.000E-05, -1.000E-05/
!...Expansion coefficients for bottom sea quark distribution.
  Data (((cehlq(ix,it,nx,7,1),ix=1,6),it=1,6), nx=1, 2)/9.010E-03, -1.401E-02, 7.150E-03, -4.130E-03, 1.260E-03, -1.040E-03, 6.280E-03, -9.320E-03, 4.780E-03, -2.890E-03, 9.100E-04, -8.200E-04, -2.930E-03, 4.090E-03, -1.890E-03, 7.600E-04, -2.300E-04, 1.400E-04, 3.900E-04, -1.200E-03, 4.400E-04, -2.500E-04, 2.000E-05, -2.000E-05, 2.600E-04, 1.400E-04, -8.000E-05, 1.000E-04, 1.000E-05, 1.000E-05, -2.600E-04, 3.200E-04, 1.000E-05, -1.000E-05, 1.000E-05, -1.000E-05, 8.029E-01, -1.075E+00, 3.792E-01, -7.843E-02, 1.007E-02, -1.090E-03, 7.903E-01, -1.099E+00, 4.153E-01, -9.301E-02, 1.317E-02, -1.410E-03, -1.704E-02, -1.130E-02, 2.882E-02, -1.341E-02, 3.040E-03, -3.600E-04, -7.200E-04, 7.230E-03, -5.160E-03, 1.080E-03, -5.000E-05, -4.000E-05, 3.050E-03, -4.610E-03, 1.660E-03, -1.300E-04, -1.000E-05, 1.000E-05, -4.360E-03, 5.230E-03, -1.610E-03, 2.000E-04, -2.000E-05, 0.000E+00/
  Data (((cehlq(ix,it,nx,7,2),ix=1,6),it=1,6), nx=1, 2)/8.980E-03, -1.459E-02, 7.510E-03, -4.410E-03, 1.310E-03, -1.070E-03, 5.970E-03, -9.440E-03, 4.800E-03, -3.020E-03, 9.100E-04, -8.500E-04, -3.050E-03, 4.440E-03, -2.100E-03, 8.500E-04, -2.400E-04, 1.400E-04, 5.300E-04, -1.300E-03, 5.600E-04, -2.700E-04, 3.000E-05, -2.000E-05, 2.000E-04, 1.400E-04, -1.100E-04, 1.000E-04, 0.000E+00, 0.000E+00, -2.600E-04, 3.200E-04, 0.000E+00, -3.000E-05, 1.000E-05, -1.000E-05, 8.672E-01, -1.174E+00, 4.265E-01, -9.252E-02, 1.244E-02, -1.460E-03, 8.500E-01, -1.194E+00, 4.630E-01, -1.083E-01, 1.614E-02, -1.830E-03, -2.241E-02, -5.630E-03, 2.815E-02, -1.425E-02, 3.520E-03, -4.300E-04, -7.300E-04, 8.030E-03, -5.780E-03, 1.380E-03, -1.300E-04, -4.000E-05, 3.460E-03, -5.380E-03, 1.960E-03, -2.100E-04, 1.000E-05, 1.000E-05, -4.850E-03, 5.950E-03, -1.890E-03, 2.600E-04, -3.000E-05, 0.000E+00/
!...Expansion coefficients for top sea quark distribution.
  Data (((cehlq(ix,it,nx,8,1),ix=1,6),it=1,6), nx=1, 2)/4.410E-03, -7.480E-03, 3.770E-03, -2.580E-03, 7.300E-04, -7.100E-04, 3.840E-03, -6.050E-03, 3.030E-03, -2.030E-03, 5.800E-04, -5.900E-04, -8.800E-04, 1.660E-03, -7.500E-04, 4.700E-04, -1.000E-04, 1.000E-04, -8.000E-05, -1.500E-04, 1.200E-04, -9.000E-05, 3.000E-05, 0.000E+00, 1.300E-04, -2.200E-04, -2.000E-05, -2.000E-05, -2.000E-05, -2.000E-05, -7.000E-05, 1.900E-04, -4.000E-05, 2.000E-05, 0.000E+00, 0.000E+00, 6.623E-01, -9.248E-01, 3.519E-01, -7.930E-02, 1.110E-02, -1.180E-03, 6.380E-01, -9.062E-01, 3.582E-01, -8.479E-02, 1.265E-02, -1.390E-03, -2.581E-02, 2.125E-02, 4.190E-03, -4.980E-03, 1.490E-03, -2.100E-04, 7.100E-04, 5.300E-04, -1.270E-03, 3.900E-04, -5.000E-05, -1.000E-05, 3.850E-03, -5.060E-03, 1.860E-03, -3.500E-04, 4.000E-05, 0.000E+00, -3.530E-03, 4.460E-03, -1.500E-03, 2.700E-04, -3.000E-05, 0.000E+00/
  Data (((cehlq(ix,it,nx,8,2),ix=1,6),it=1,6), nx=1, 2)/4.260E-03, -7.530E-03, 3.830E-03, -2.680E-03, 7.600E-04, -7.300E-04, 3.640E-03, -6.050E-03, 3.030E-03, -2.090E-03, 5.900E-04, -6.000E-04, -9.200E-04, 1.710E-03, -8.200E-04, 5.000E-04, -1.200E-04, 1.000E-04, -5.000E-05, -1.600E-04, 1.300E-04, -9.000E-05, 3.000E-05, 0.000E+00, 1.300E-04, -2.100E-04, -1.000E-05, -2.000E-05, -2.000E-05, -1.000E-05, -8.000E-05, 1.800E-04, -5.000E-05, 2.000E-05, 0.000E+00, 0.000E+00, 7.146E-01, -1.007E+00, 3.932E-01, -9.246E-02, 1.366E-02, -1.540E-03, 6.856E-01, -9.828E-01, 3.977E-01, -9.795E-02, 1.540E-02, -1.790E-03, -3.053E-02, 2.758E-02, 2.150E-03, -4.880E-03, 1.640E-03, -2.500E-04, 9.200E-04, 4.200E-04, -1.340E-03, 4.600E-04, -8.000E-05, -1.000E-05, 4.230E-03, -5.660E-03, 2.140E-03, -4.300E-04, 6.000E-05, 0.000E+00, -3.890E-03, 5.000E-03, -1.740E-03, 3.300E-04, -4.000E-05, 0.000E+00/

!...The following data lines are coefficients needed in the
!...Duke, Owens proton structure function parametrizations, see below.
!...Expansion coefficients for (up+down) valence quark distribution.
  Data ((cdo(ip,is,1,1),is=1,6), ip=1, 3)/4.190E-01, 3.460E+00, 4.400E+00, 0.000E+00, 0.000E+00, 0.000E+00, 4.000E-03, 7.240E-01, -4.860E+00, 0.000E+00, 0.000E+00, 0.000E+00, -7.000E-03, -6.600E-02, 1.330E+00, 0.000E+00, 0.000E+00, 0.000E+00/
  Data ((cdo(ip,is,1,2),is=1,6), ip=1, 3)/3.740E-01, 3.330E+00, 6.030E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.400E-02, 7.530E-01, -6.220E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, -7.600E-02, 1.560E+00, 0.000E+00, 0.000E+00, 0.000E+00/
!...Expansion coefficients for down valence quark distribution.
  Data ((cdo(ip,is,2,1),is=1,6), ip=1, 3)/7.630E-01, 4.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, -2.370E-01, 6.270E-01, -4.210E-01, 0.000E+00, 0.000E+00, 0.000E+00, 2.600E-02, -1.900E-02, 3.300E-02, 0.000E+00, 0.000E+00, 0.000E+00/
  Data ((cdo(ip,is,2,2),is=1,6), ip=1, 3)/7.610E-01, 3.830E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, -2.320E-01, 6.270E-01, -4.180E-01, 0.000E+00, 0.000E+00, 0.000E+00, 2.300E-02, -1.900E-02, 3.600E-02, 0.000E+00, 0.000E+00, 0.000E+00/
!...Expansion coefficients for (up+down+strange) sea quark distribution.
  Data ((cdo(ip,is,3,1),is=1,6), ip=1, 3)/1.265E+00, 0.000E+00, 8.050E+00, 0.000E+00, 0.000E+00, 0.000E+00, -1.132E+00, -3.720E-01, 1.590E+00, 6.310E+00, -1.050E+01, 1.470E+01, 2.930E-01, -2.900E-02, -1.530E-01, -2.730E-01, -3.170E+00, 9.800E+00/
  Data ((cdo(ip,is,3,2),is=1,6), ip=1, 3)/1.670E+00, 0.000E+00, 9.150E+00, 0.000E+00, 0.000E+00, 0.000E+00, -1.920E+00, -2.730E-01, 5.300E-01, 1.570E+01, -1.010E+02, 2.230E+02, 5.820E-01, -1.640E-01, -7.630E-01, -2.830E+00, 4.470E+01, -1.170E+02/
!...Expansion coefficients for charm sea quark distribution.
  Data ((cdo(ip,is,4,1),is=1,6), ip=1, 3)/0.000E+00, -3.600E-02, 6.350E+00, 0.000E+00, 0.000E+00, 0.000E+00, 1.350E-01, -2.220E-01, 3.260E+00, -3.030E+00, 1.740E+01, -1.790E+01, -7.500E-02, -5.800E-02, -9.090E-01, 1.500E+00, -1.130E+01, 1.560E+01/
  Data ((cdo(ip,is,4,2),is=1,6), ip=1, 3)/0.000E+00, -1.200E-01, 3.510E+00, 0.000E+00, 0.000E+00, 0.000E+00, 6.700E-02, -2.330E-01, 3.660E+00, -4.740E-01, 9.500E+00, -1.660E+01, -3.100E-02, -2.300E-02, -4.530E-01, 3.580E-01, -5.430E+00, 1.550E+01/
!...Expansion coefficients for gluon distribution.
  Data ((cdo(ip,is,5,1),is=1,6), ip=1, 3)/1.560E+00, 0.000E+00, 6.000E+00, 9.000E+00, 0.000E+00, 0.000E+00, -1.710E+00, -9.490E-01, 1.440E+00, -7.190E+00, -1.650E+01, 1.530E+01, 6.380E-01, 3.250E-01, -1.050E+00, 2.550E-01, 1.090E+01, -1.010E+01/
  Data ((cdo(ip,is,5,2),is=1,6), ip=1, 3)/8.790E-01, 0.000E+00, 4.000E+00, 9.000E+00, 0.000E+00, 0.000E+00, -9.710E-01, -1.160E+00, 1.230E+00, -5.640E+00, -7.540E+00, -5.960E-01, 4.340E-01, 4.760E-01, -2.540E-01, -8.170E-01, 5.500E+00, 1.260E-01/

!...The following data lines are coefficients needed in the
!...Owens pion structure function parametrizations, see below.
!...Expansion coefficients for up and down valence quark distributions.
  Data ((cow(ip,is,1,1),is=1,5), ip=1, 3)/4.0000E-01, 7.0000E-01, 0.0000E+00, 0.0000E+00, 0.0000E+00, -6.2120E-02, 6.4780E-01, 0.0000E+00, 0.0000E+00, 0.0000E+00, -7.1090E-03, 1.3350E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00/
  Data ((cow(ip,is,1,2),is=1,5), ip=1, 3)/4.0000E-01, 6.2800E-01, 0.0000E+00, 0.0000E+00, 0.0000E+00, -5.9090E-02, 6.4360E-01, 0.0000E+00, 0.0000E+00, 0.0000E+00, -6.5240E-03, 1.4510E-02, 0.0000E+00, 0.0000E+00, 0.0000E+00/
!...Expansion coefficients for gluon distribution.
  Data ((cow(ip,is,2,1),is=1,5), ip=1, 3)/8.8800E-01, 0.0000E+00, 3.1100E+00, 6.0000E+00, 0.0000E+00, -1.8020E+00, -1.5760E+00, -1.3170E-01, 2.8010E+00, -1.7280E+01, 1.8120E+00, 1.2000E+00, 5.0680E-01, -1.2160E+01, 2.0490E+01/
  Data ((cow(ip,is,2,2),is=1,5), ip=1, 3)/7.9400E-01, 0.0000E+00, 2.8900E+00, 6.0000E+00, 0.0000E+00, -9.1440E-01, -1.2370E+00, 5.9660E-01, -3.6710E+00, -8.1910E+00, 5.9660E-01, 6.5820E-01, -2.5500E-01, -2.3040E+00, 7.7580E+00/
!...Expansion coefficients for (up+down+strange) quark sea distribution.
  Data ((cow(ip,is,3,1),is=1,5), ip=1, 3)/9.0000E-01, 0.0000E+00, 5.0000E+00, 0.0000E+00, 0.0000E+00, -2.4280E-01, -2.1200E-01, 8.6730E-01, 1.2660E+00, 2.3820E+00, 1.3860E-01, 3.6710E-03, 4.7470E-02, -2.2150E+00, 3.4820E-01/
  Data ((cow(ip,is,3,2),is=1,5), ip=1, 3)/9.0000E-01, 0.0000E+00, 5.0000E+00, 0.0000E+00, 0.0000E+00, -1.4170E-01, -1.6970E-01, -2.4740E+00, -2.5340E+00, 5.6210E-01, -1.7400E-01, -9.6230E-02, 1.5750E+00, 1.3780E+00, -2.7010E-01/
!...Expansion coefficients for charm quark sea distribution.
  Data ((cow(ip,is,4,1),is=1,5), ip=1, 3)/0.0000E+00, -2.2120E-02, 2.8940E+00, 0.0000E+00, 0.0000E+00, 7.9280E-02, -3.7850E-01, 9.4330E+00, 5.2480E+00, 8.3880E+00, -6.1340E-02, -1.0880E-01, -1.0852E+01, -7.1870E+00, -1.1610E+01/
  Data ((cow(ip,is,4,2),is=1,5), ip=1, 3)/0.0000E+00, -8.8200E-02, 1.9240E+00, 0.0000E+00, 0.0000E+00, 6.2290E-02, -2.8920E-01, 2.4240E-01, -4.4630E+00, -8.3670E-01, -4.0990E-02, -1.0820E-01, 2.0360E+00, 5.2090E+00, -4.8400E-02/

!...Euler's beta function, requires ordinary Gamma function
!lin-10/25/02 get rid of argument usage mismatch in PYGAMM():
!      EULBT(X,Y)=PYGAMM(X)*PYGAMM(Y)/PYGAMM(X+Y)

!...Reset structure functions, check x and hadron flavour.
  alam = 0.
  Do kfl = -6, 6
    xpq(kfl) = 0.
  End Do
  If (x<0. .Or. x>1.) Then
    Write (mstu(11), 1000) x
    Return
  End If
  kfa = iabs(kf)
  If (kfa/=211 .And. kfa/=2212 .And. kfa/=2112) Then
    Write (mstu(11), 1100) kf
    Return
  End If

!...Call user-supplied structure function. Select proton/neutron/pion.
  If (mstp(51)==0 .Or. mstp(52)>=2) Then
    kfe = kfa
    If (kfa==2112) kfe = 2212
    Call pystfe(kfe, x, q2, xpq)
    Goto 230
  End If
  If (kfa==211) Goto 200

  If (mstp(51)==1 .Or. mstp(51)==2) Then
!...Proton structure functions from Eichten, Hinchliffe, Lane, Quigg.
!...Allowed variable range: 5 GeV2 < Q2 < 1E8 GeV2; 1E-4 < x < 1

!...Determine set, Lamdba and x and t expansion variables.
    nset = mstp(51)
    If (nset==1) alam = 0.2
    If (nset==2) alam = 0.29
    tmin = log(5./alam**2)
    tmax = log(1E8/alam**2)
    If (mstp(52)==0) Then
      t = tmin
    Else
      t = log(q2/alam**2)
    End If
    vt = max(-1., min(1.,(2.*t-tmax-tmin)/(tmax-tmin)))
    nx = 1
    If (x<=0.1) nx = 2
    If (nx==1) vx = (2.*x-1.1)/0.9
    If (nx==2) vx = max(-1., (2.*log(x)+11.51293)/6.90776)
    cxs = 1.
    If (x<1E-4 .And. abs(parp(51)-1.)>0.01) cxs = (1E-4/x)**(parp(51)-1.)

!...Chebyshev polynomials for x and t expansion.
    tx(1) = 1.
    tx(2) = vx
    tx(3) = 2.*vx**2 - 1.
    tx(4) = 4.*vx**3 - 3.*vx
    tx(5) = 8.*vx**4 - 8.*vx**2 + 1.
    tx(6) = 16.*vx**5 - 20.*vx**3 + 5.*vx
    tt(1) = 1.
    tt(2) = vt
    tt(3) = 2.*vt**2 - 1.
    tt(4) = 4.*vt**3 - 3.*vt
    tt(5) = 8.*vt**4 - 8.*vt**2 + 1.
    tt(6) = 16.*vt**5 - 20.*vt**3 + 5.*vt

!...Calculate structure functions.
    Do kfl = 1, 6
      xqsum = 0.
      Do it = 1, 6
        Do ix = 1, 6
          xqsum = xqsum + cehlq(ix, it, nx, kfl, nset)*tx(ix)*tt(it)
        End Do
      End Do
      xq(kfl) = xqsum*(1.-x)**nehlq(kfl, nset)*cxs
    End Do

!...Put into output array.
    xpq(0) = xq(4)
    xpq(1) = xq(2) + xq(3)
    xpq(2) = xq(1) + xq(3)
    xpq(3) = xq(5)
    xpq(4) = xq(6)
    xpq(-1) = xq(3)
    xpq(-2) = xq(3)
    xpq(-3) = xq(5)
    xpq(-4) = xq(6)

!...Special expansion for bottom (thresh effects).
    If (mstp(54)>=5) Then
      If (nset==1) tmin = 8.1905
      If (nset==2) tmin = 7.4474
      If (t<=tmin) Goto 140
      vt = max(-1., min(1.,(2.*t-tmax-tmin)/(tmax-tmin)))
      tt(1) = 1.
      tt(2) = vt
      tt(3) = 2.*vt**2 - 1.
      tt(4) = 4.*vt**3 - 3.*vt
      tt(5) = 8.*vt**4 - 8.*vt**2 + 1.
      tt(6) = 16.*vt**5 - 20.*vt**3 + 5.*vt
      xqsum = 0.
      Do it = 1, 6
        Do ix = 1, 6
          xqsum = xqsum + cehlq(ix, it, nx, 7, nset)*tx(ix)*tt(it)
        End Do
      End Do
      xpq(5) = xqsum*(1.-x)**nehlq(7, nset)
      xpq(-5) = xpq(5)
      140 Continue
    End If

!...Special expansion for top (thresh effects).
    If (mstp(54)>=6) Then
      If (nset==1) tmin = 11.5528
      If (nset==2) tmin = 10.8097
      tmin = tmin + 2.*log(pmas(6,1)/30.)
      tmax = tmax + 2.*log(pmas(6,1)/30.)
      If (t<=tmin) Goto 160
      vt = max(-1., min(1.,(2.*t-tmax-tmin)/(tmax-tmin)))
      tt(1) = 1.
      tt(2) = vt
      tt(3) = 2.*vt**2 - 1.
      tt(4) = 4.*vt**3 - 3.*vt
      tt(5) = 8.*vt**4 - 8.*vt**2 + 1.
      tt(6) = 16.*vt**5 - 20.*vt**3 + 5.*vt
      xqsum = 0.
      Do it = 1, 6
        Do ix = 1, 6
          xqsum = xqsum + cehlq(ix, it, nx, 8, nset)*tx(ix)*tt(it)
        End Do
      End Do
      xpq(6) = xqsum*(1.-x)**nehlq(8, nset)
      xpq(-6) = xpq(6)
      160 Continue
    End If

  Else If (mstp(51)==3 .Or. mstp(51)==4) Then
!...Proton structure functions from Duke, Owens.
!...Allowed variable range: 4 GeV2 < Q2 < approx 1E6 GeV2.

!...Determine set, Lambda and s expansion parameter.
    nset = mstp(51) - 2
    If (nset==1) alam = 0.2
    If (nset==2) alam = 0.4
    If (mstp(52)<=0) Then
      sd = 0.
    Else
      sd = log(log(max(q2,4.)/alam**2)/log(4./alam**2))
    End If

!...Calculate structure functions.
    Do kfl = 1, 5
      Do is = 1, 6
        ts(is) = cdo(1, is, kfl, nset) + cdo(2, is, kfl, nset)*sd + cdo(3, is, kfl, nset)*sd**2
      End Do
      If (kfl<=2) Then

!lin-10/25/02 evaluate EULBT(TS(1),TS(2)+1.):
!          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)*(1.+TS(3)*X)/(EULBT(TS(1),
!     &    TS(2)+1.)*(1.+TS(3)*TS(1)/(TS(1)+TS(2)+1.)))
        eulbt1 = pygamm(ts(1))*pygamm(ts(2)+1.)/pygamm(ts(1)+ts(2)+1.)
        xq(kfl) = x**ts(1)*(1.-x)**ts(2)*(1.+ts(3)*x)/(eulbt1*(1.+ts(3)*ts(1)/(ts(1)+ts(2)+1.)))
      Else
        xq(kfl) = ts(1)*x**ts(2)*(1.-x)**ts(3)*(1.+ts(4)*x+ts(5)*x**2+ts(6)*x**3)
      End If


    End Do

!...Put into output arrays.
    xpq(0) = xq(5)
    xpq(1) = xq(2) + xq(3)/6.
    xpq(2) = 3.*xq(1) - xq(2) + xq(3)/6.
    xpq(3) = xq(3)/6.
    xpq(4) = xq(4)
    xpq(-1) = xq(3)/6.
    xpq(-2) = xq(3)/6.
    xpq(-3) = xq(3)/6.
    xpq(-4) = xq(4)

!...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli.
!...These are accessed via PYSTFE since the files needed may not always
!...available.
  Else If (mstp(51)>=11 .And. mstp(51)<=13) Then
    Call pystfe(2212, x, q2, xpq)

!...Unknown proton parametrization.
  Else
    Write (mstu(11), 1200) mstp(51)
  End If
  Goto 230

  200 If ((mstp(51)>=1 .And. mstp(51)<=4) .Or. (mstp(51)>=11 .And. mstp(51)<=13)) Then
!...Pion structure functions from Owens.
!...Allowed variable range: 4 GeV2 < Q2 < approx 2000 GeV2.

!...Determine set, Lambda and s expansion variable.
    nset = 1
    If (mstp(51)==2 .Or. mstp(51)==4 .Or. mstp(51)==13) nset = 2
    If (nset==1) alam = 0.2
    If (nset==2) alam = 0.4
    If (mstp(52)<=0) Then
      sd = 0.
    Else
      sd = log(log(max(q2,4.)/alam**2)/log(4./alam**2))
    End If

!...Calculate structure functions.
    Do kfl = 1, 4
      Do is = 1, 5
        ts(is) = cow(1, is, kfl, nset) + cow(2, is, kfl, nset)*sd + cow(3, is, kfl, nset)*sd**2
      End Do
      If (kfl==1) Then

!lin-10/25/02 get rid of argument usage mismatch in PYGAMM():
!          XQ(KFL)=X**TS(1)*(1.-X)**TS(2)/EULBT(TS(1),TS(2)+1.)
        eulbt2 = pygamm(ts(1))*pygamm(ts(2)+1.)/pygamm(ts(1)+ts(2)+1.)
        xq(kfl) = x**ts(1)*(1.-x)**ts(2)/eulbt2
      Else
        xq(kfl) = ts(1)*x**ts(2)*(1.-x)**ts(3)*(1.+ts(4)*x+ts(5)*x**2)
      End If
    End Do

!...Put into output arrays.
    xpq(0) = xq(2)
    xpq(1) = xq(3)/6.
    xpq(2) = xq(1) + xq(3)/6.
    xpq(3) = xq(3)/6.
    xpq(4) = xq(4)
    xpq(-1) = xq(1) + xq(3)/6.
    xpq(-2) = xq(3)/6.
    xpq(-3) = xq(3)/6.
    xpq(-4) = xq(4)

!...Unknown pion parametrization.
  Else
    Write (mstu(11), 1200) mstp(51)
  End If

!...Isospin conjugation for neutron, charge conjugation for antipart.
  230 If (kfa==2112) Then
    xps = xpq(1)
    xpq(1) = xpq(2)
    xpq(2) = xps
    xps = xpq(-1)
    xpq(-1) = xpq(-2)
    xpq(-2) = xps
  End If
  If (kf<0) Then
    Do kfl = 1, 4
      xps = xpq(kfl)
      xpq(kfl) = xpq(-kfl)
      xpq(-kfl) = xps
    End Do
  End If

!...Check positivity and reset above maximum allowed flavour.
  Do kfl = -6, 6
    xpq(kfl) = max(0., xpq(kfl))
    If (iabs(kfl)>mstp(54)) xpq(kfl) = 0.
  End Do

!...consider nuclear effect on the structure function
  If ((jbt/=1 .And. jbt/=2) .Or. ihpr2(6)==0 .Or. ihnt2(16)==1) Goto 400
  atnm = ihnt2(2*jbt-1)
  If (atnm<=1.0) Goto 400
  If (jbt==1) Then
    bbr2 = (yp(1,ihnt2(11))**2+yp(2,ihnt2(11))**2)/1.44/atnm**0.66666
  Else If (jbt==2) Then
    bbr2 = (yt(1,ihnt2(12))**2+yt(2,ihnt2(12))**2)/1.44/atnm**0.66666
  End If
  bbr2 = min(1.0, bbr2)
  abx = (atnm**0.33333333-1.0)
  apx = hipr1(6)*4.0/3.0*abx*sqrt(1.0-bbr2)
  aax = 1.192*alog(atnm)**0.1666666
  rrx = aax*(x**3-1.2*x**2+0.21*x) + 1.0 - (apx-1.079*abx*sqrt(x)/alog(atnm+1.0))*exp(-x**2.0/0.01)
  Do kfl = -6, 6
    xpq(kfl) = xpq(kfl)*rrx
  End Do
!                        ********consider the nuclear effect on the structure
!                                function which also depends on the impact
!                                parameter of the nuclear reaction

  400 Continue

  Return
!...Formats for error printouts.
  1000 Format (' Error: x value outside physical range, x =', 1P, E12.3)
  1100 Format (' Error: illegal particle code for structure function,', ' KF =', I5)
  1200 Format (' Error: bad value of parameter MSTP(51) in PYSTFU,', ' MSTP(51) =', I5)
End Subroutine pystfu

!*********************************************************************

Subroutine pyspli(kf, kflin, kflch, kflsp)

!...In case of a hadron remnant which is more complicated than just a
!...quark or a diquark, split it into two (partons or hadron + parton).
  Dimension kfl(3)

!...Preliminaries. Parton composition.
  kfa = iabs(kf)
  kfs = isign(1, kf)
  kfl(1) = mod(kfa/1000, 10)
  kfl(2) = mod(kfa/100, 10)
  kfl(3) = mod(kfa/10, 10)
  kflr = kflin*kfs
  kflch = 0

!...Subdivide meson.
  If (kfl(1)==0) Then
    kfl(2) = kfl(2)*(-1)**kfl(2)
    kfl(3) = -kfl(3)*(-1)**iabs(kfl(2))
    If (kflr==kfl(2)) Then
      kflsp = kfl(3)
    Else If (kflr==kfl(3)) Then
      kflsp = kfl(2)
    Else If (iabs(kflr)==21 .And. rlu(0)>0.5) Then
      kflsp = kfl(2)
      kflch = kfl(3)
    Else If (iabs(kflr)==21) Then
      kflsp = kfl(3)
      kflch = kfl(2)
    Else If (kflr*kfl(2)>0) Then
      Call lukfdi(-kflr, kfl(2), kfdump, kflch)
      kflsp = kfl(3)
    Else
      Call lukfdi(-kflr, kfl(3), kfdump, kflch)
      kflsp = kfl(2)
    End If

!...Subdivide baryon.
  Else
    nagr = 0
    Do j = 1, 3
      If (kflr==kfl(j)) nagr = nagr + 1
    End Do
    If (nagr>=1) Then
      ragr = 0.00001 + (nagr-0.00002)*rlu(0)
      iagr = 0
      Do j = 1, 3
        If (kflr==kfl(j)) ragr = ragr - 1.
        If (iagr==0 .And. ragr<=0.) iagr = j
      End Do
    Else
      iagr = int(1.00001+2.99998*rlu(0))
    End If
    id1 = 1
    If (iagr==1) id1 = 2
    If (iagr==1 .And. kfl(3)>kfl(2)) id1 = 3
    id2 = 6 - iagr - id1
    ksp = 3
    If (mod(kfa,10)==2 .And. kfl(1)==kfl(2)) Then
      If (iagr/=3 .And. rlu(0)>0.25) ksp = 1
    Else If (mod(kfa,10)==2 .And. kfl(2)>=kfl(3)) Then
      If (iagr/=1 .And. rlu(0)>0.25) ksp = 1
    Else If (mod(kfa,10)==2) Then
      If (iagr==1) ksp = 1
      If (iagr/=1 .And. rlu(0)>0.75) ksp = 1
    End If
    kflsp = 1000*kfl(id1) + 100*kfl(id2) + ksp
    If (kflin==21) Then
      kflch = kfl(iagr)
    Else If (nagr==0 .And. kflr>0) Then
      Call lukfdi(-kflr, kfl(iagr), kfdump, kflch)
    Else If (nagr==0) Then
      Call lukfdi(10000+kflsp, -kflr, kfdump, kflch)
      kflsp = kfl(iagr)
    End If
  End If

!...Add on correct sign for result.
  kflch = kflch*kfs
  kflsp = kflsp*kfs

  Return
End Subroutine pyspli

!*********************************************************************

Function pygamm(x)

!...Gives ordinary Gamma function Gamma(x) for positive, real arguments;
!...see M. Abramowitz, I. A. Stegun: Handbook of Mathematical Functions
!...(Dover, 1965) 6.1.36.
  Dimension b(8)
!lin      DATA B/-0.577191652,0.988205891,-0.897056937,0.918206857,
!lin     &-0.756704078,0.482199394,-0.193527818,0.035868343/
  Data b/ -0.57719165, 0.98820589, -0.89705694, 0.91820686, -0.75670408, 0.48219939, -0.19352782, 0.03586834/

  nx = int(x)
  dx = x - nx

  pygamm = 1.
  Do i = 1, 8
    pygamm = pygamm + b(i)*dx**i
  End Do
  If (x<1.) Then
    pygamm = pygamm/x
  Else
    Do ix = 1, nx - 1
      pygamm = (x-ix)*pygamm
    End Do
  End If

  Return
End Function pygamm

!***********************************************************************

Function pyw1au(eps, ireim)

!...Calculates real and imaginary parts of the auxiliary function W1;
!...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
!...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/

!lin-8/2014:
!      ASINH(X)=LOG(X+SQRT(X**2+1.))
  acosh(x) = log(x+sqrt(x**2-1.))

  If (eps<0.) Then
    w1re = 2.*sqrt(1.-eps)*asinh(sqrt(-1./eps))
    w1im = 0.
  Else If (eps<1.) Then
    w1re = 2.*sqrt(1.-eps)*acosh(sqrt(1./eps))
    w1im = -paru(1)*sqrt(1.-eps)
  Else
    w1re = 2.*sqrt(eps-1.)*asin(sqrt(1./eps))
    w1im = 0.
  End If

  If (ireim==1) pyw1au = w1re
  If (ireim==2) pyw1au = w1im

  Return
End Function pyw1au

!***********************************************************************

Function pyw2au(eps, ireim)

!...Calculates real and imaginary parts of the auxiliary function W2;
!...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
!...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/

!lin-8/2014:
!      ASINH(X)=LOG(X+SQRT(X**2+1.))
  acosh(x) = log(x+sqrt(x**2-1.))

  If (eps<0.) Then
    w2re = 4.*(asinh(sqrt(-1./eps)))**2
    w2im = 0.
  Else If (eps<1.) Then
    w2re = 4.*(acosh(sqrt(1./eps)))**2 - paru(1)**2
    w2im = -4.*paru(1)*acosh(sqrt(1./eps))
  Else
    w2re = -4.*(asin(sqrt(1./eps)))**2
    w2im = 0.
  End If

  If (ireim==1) pyw2au = w2re
  If (ireim==2) pyw2au = w2im

  Return
End Function pyw2au

!***********************************************************************

Function pyi3au(be, eps, ireim)

!...Calculates real and imaginary parts of the auxiliary function I3;
!...see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van der Bij,
!...FERMILAB-Pub-87/100-T, LBL-23504, June, 1987
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/

  If (eps<1.) ga = 0.5*(1.+sqrt(1.-eps))

  If (eps<0.) Then
    f3re = pyspen((ga-1.)/(ga+be-1.), 0., 1) - pyspen(ga/(ga+be-1.), 0., 1) + pyspen((be-ga)/be, 0., 1) - pyspen((be-ga)/(be-1.), 0., 1) + (log(be)**2-log(be-1.)**2)/2. + log(ga)*log((ga+be-1.)/be) + log(ga-1.)*log((be-1.)/(ga+be-1.))
    f3im = 0.
  Else If (eps<1.) Then
    f3re = pyspen((ga-1.)/(ga+be-1.), 0., 1) - pyspen(ga/(ga+be-1.), 0., 1) + pyspen(ga/(ga-be), 0., 1) - pyspen((ga-1.)/(ga-be), 0., 1) + log(ga/(1.-ga))*log((ga+be-1.)/(be-ga))
    f3im = -paru(1)*log((ga+be-1.)/(be-ga))
  Else
    rsq = eps/(eps-1.+(2.*be-1.)**2)
    rcthe = rsq*(1.-2.*be/eps)
    rsthe = sqrt(rsq-rcthe**2)
    rcphi = rsq*(1.+2.*(be-1.)/eps)
    rsphi = sqrt(rsq-rcphi**2)
    r = sqrt(rsq)
    the = acos(rcthe/r)
    phi = acos(rcphi/r)
    f3re = pyspen(rcthe, rsthe, 1) + pyspen(rcthe, -rsthe, 1) - pyspen(rcphi, rsphi, 1) - pyspen(rcphi, -rsphi, 1) + (phi-the)*(phi+the-paru(1))
    f3im = pyspen(rcthe, rsthe, 2) + pyspen(rcthe, -rsthe, 2) - pyspen(rcphi, rsphi, 2) - pyspen(rcphi, -rsphi, 2)
  End If

  If (ireim==1) pyi3au = 2./(2.*be-1.)*f3re
  If (ireim==2) pyi3au = 2./(2.*be-1.)*f3im

  Return
End Function pyi3au

!***********************************************************************

Function pyspen(xrein, ximin, ireim)

!...Calculates real and imaginary part of Spence function; see
!...G. 't Hooft and M. Veltman, Nucl. Phys. B153 (1979) 365.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Dimension b(0:14)

  Data b/1.000000E+00, -5.000000E-01, 1.666667E-01, 0.000000E+00, -3.333333E-02, 0.000000E+00, 2.380952E-02, 0.000000E+00, -3.333333E-02, 0.000000E+00, 7.575757E-02, 0.000000E+00, -2.531135E-01, 0.000000E+00, 1.166667E+00/

  xre = xrein
  xim = ximin
  If (abs(1.-xre)<1.E-6 .And. abs(xim)<1.E-6) Then
    If (ireim==1) pyspen = paru(1)**2/6.
    If (ireim==2) pyspen = 0.
    Return
  End If

  xmod = sqrt(xre**2+xim**2)
  If (xmod<1.E-6) Then
    If (ireim==1) pyspen = 0.
    If (ireim==2) pyspen = 0.
    Return
  End If

  xarg = sign(acos(xre/xmod), xim)
  sp0re = 0.
  sp0im = 0.
  sgn = 1.
  If (xmod>1.) Then
    algxre = log(xmod)
    algxim = xarg - sign(paru(1), xarg)
    sp0re = -paru(1)**2/6. - (algxre**2-algxim**2)/2.
    sp0im = -algxre*algxim
    sgn = -1.
    xmod = 1./xmod
    xarg = -xarg
    xre = xmod*cos(xarg)
    xim = xmod*sin(xarg)
  End If
  If (xre>0.5) Then
    algxre = log(xmod)
    algxim = xarg
    xre = 1. - xre
    xim = -xim
    xmod = sqrt(xre**2+xim**2)
    xarg = sign(acos(xre/xmod), xim)
    algyre = log(xmod)
    algyim = xarg
    sp0re = sp0re + sgn*(paru(1)**2/6.-(algxre*algyre-algxim*algyim))
    sp0im = sp0im - sgn*(algxre*algyim+algxim*algyre)
    sgn = -sgn
  End If

  xre = 1. - xre
  xim = -xim
  xmod = sqrt(xre**2+xim**2)
  xarg = sign(acos(xre/xmod), xim)
  zre = -log(xmod)
  zim = -xarg

  spre = 0.
  spim = 0.
  savere = 1.
  saveim = 0.
  Do i = 0, 14
    termre = (savere*zre-saveim*zim)/float(i+1)
    termim = (savere*zim+saveim*zre)/float(i+1)
    savere = termre
    saveim = termim
    spre = spre + b(i)*termre
    spim = spim + b(i)*termim
  End Do

  If (ireim==1) pyspen = sp0re + sgn*spre
  If (ireim==2) pyspen = sp0im + sgn*spim

  Return
End Function pyspen

!*********************************************************************

Block Data pydata

!...Give sensible default values to all status codes and parameters.
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
  Save /pysubs/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Common /pyint1/mint(400), vint(400)
  Save /pyint1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
  Save /pyint2/
  Common /pyint3/xsfx(2, -40:40), isig(1000, 3), sigh(1000)
  Save /pyint3/
  Common /pyint4/widp(21:40, 0:40), wide(21:40, 0:40), wids(21:40, 3)
  Save /pyint4/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
  Save /pyint5/
  Common /pyint6/proc(0:200)
  Character proc*28
  Save /pyint6/

!...Default values for allowed processes and kinematics constraints.
  Data msel/1/
  Data msub/200*0/
  Data ((kfin(i,j),j=-40,40), i=1, 2)/40*1, 0, 80*1, 0, 40*1/
  Data ckin/2.0, -1.0, 0.0, -1.0, 1.0, 1.0, -10., 10., -10., 10., -10., 10., -10., 10., -10., 10., -1.0, 1.0, -1.0, 1.0, 0.0, 1.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0., 0., 2.0, -1.0, 0., 0., 0., 0., 0., 0., 0., 0., 160*0./

!...Default values for main switches and parameters. Reset information.
  Data (mstp(i), i=1, 100)/3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 3, 7, 1, 0, 0, 0, 0, 0, 1, 1, 20, 6, 0, 0, 0, 0, 0, 0, 1, 2, 2, 2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 100, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0/
  Data (mstp(i), i=101, 200)/1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 20, 0, 0, 0, 0, 0, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 3, 1989, 11, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  Data (parp(i), i=1, 100)/0.25, 10., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.5, 2.0, 0.075, 0., 0.2, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.0, 2.26, 1.E4, 1.E-4, 0., 0., 0., 0., 0., 0., 0.25, 1.0, 0.25, 1.0, 2.0, 1.E-3, 4.0, 0., 0., 0., 4.0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.6, 1.85, 0.5, 0.2, 0.33, 0.66, 0.7, 0.5, 0., 0., 0.44, 0.44, 2.0, 1.0, 0., 3.0, 1.0, 0.75, 0., 0./
  Data (parp(i), i=101, 200)/ -0.02, 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.0, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.4, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.01, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./
  Data msti/200*0/
  Data pari/200*0./
  Data mint/400*0/
  Data vint/400*0./

!...Constants for the generation of the various processes.
  Data (iset(i), i=1, 100)/1, 1, 1, -1, 3, -1, -1, 3, -2, -2, 2, 2, 2, 2, 2, 2, -1, 2, 2, 2, -1, 2, 2, 2, 2, 2, -1, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, -1, 4, 4, 4, -1, -1, 4, 4, -1, -1, -2, 2, 2, -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 0, -1, 0, 5, -2, -2, -2, -2/
  Data (iset(i), i=101, 200)/ -1, 1, -2, -2, -2, -2, -2, -2, -2, -2, 2, 2, 2, 2, -1, -1, -1, -2, -2, -2, -1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 1, 1, 1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2/
  Data ((kfpr(i,j),j=1,2), i=1, 50)/23, 0, 24, 0, 25, 0, 24, 0, 25, 0, 24, 0, 23, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 21, 21, 22, 21, 23, 21, 24, 21, 25, 22, 22, 22, 23, 22, 24, 22, 25, 23, 23, 23, 24, 23, 25, 24, 24, 24, 25, 25, 25, 0, 21, 0, 22, 0, 23, 0, 24, 0, 25, 0, 21, 0, 22, 0, 23, 0, 24, 0, 25, 0, 21, 0, 22, 0, 23, 0, 24, 0, 25, 0, 21, 0, 22, 0, 23, 0, 24, 0, 25, 0, 21, 0, 22, 0, 23/
  Data ((kfpr(i,j),j=1,2), i=51, 100)/0, 24, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 21, 24, 24, 22, 24, 23, 23, 24, 24, 23, 24, 23, 25, 22, 22, 23, 23, 24, 24, 24, 25, 25, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  Data ((kfpr(i,j),j=1,2), i=101, 150)/23, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 25, 0, 25, 21, 25, 22, 22, 22, 23, 23, 23, 24, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 37, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  Data ((kfpr(i,j),j=1,2), i=151, 200)/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
  Data coef/4000*0./
  Data (((icol(i,j,k),k=1,2),j=1,4), i=1, 40)/4, 0, 3, 0, 2, 0, 1, 0, 3, 0, 4, 0, 1, 0, 2, 0, 2, 0, 0, 1, 4, 0, 0, 3, 3, 0, 0, 4, 1, 0, 0, 2, 3, 0, 0, 4, 1, 4, 3, 2, 4, 0, 0, 3, 4, 2, 1, 3, 2, 0, 4, 1, 4, 0, 2, 3, 4, 0, 3, 4, 2, 0, 1, 2, 3, 2, 1, 0, 1, 4, 3, 0, 4, 3, 3, 0, 2, 1, 1, 0, 3, 2, 1, 4, 1, 0, 0, 2, 2, 4, 3, 1, 2, 0, 0, 1, 3, 2, 1, 4, 1, 4, 3, 2, 4, 2, 1, 3, 4, 2, 1, 3, 3, 4, 4, 3, 1, 2, 2, 1, 2, 0, 3, 1, 2, 0, 0, 0, 4, 2, 1, 0, 0, 0, 1, 0, 3, 0, 0, 3, 1, 2, 0, 0, 4, 0, 0, 4, 0, 0, 1, 2, 2, 0, 0, 1, 4, 4, 3, 3, 2, 2, 1, 1, 4, 4, 3, 3, 3, 3, 4, 4, 1, 1, 2, 2, 3, 2, 1, 3, 1, 2, 0, 0, 4, 2, 1, 4, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/

!...Character constants: name of processes.
  Data proc(0)/'All included subprocesses   '/
  Data (proc(i), i=1, 20)/'f + fb -> gamma*/Z0         ', 'f + fb'' -> W+/-             ', 'f + fb -> H0                ', 'gamma + W+/- -> W+/-        ', 'Z0 + Z0 -> H0               ', 'Z0 + W+/- -> W+/-           ', '                            ', 'W+ + W- -> H0               ', '                            ', '                            ', 'f + f'' -> f + f''            ', 'f + fb -> f'' + fb''          ', 'f + fb -> g + g             ', 'f + fb -> g + gamma         ', 'f + fb -> g + Z0            ', 'f + fb'' -> g + W+/-         ', 'f + fb -> g + H0            ', 'f + fb -> gamma + gamma     ', 'f + fb -> gamma + Z0        ', 'f + fb'' -> gamma + W+/-     '/
  Data (proc(i), i=21, 40)/'f + fb -> gamma + H0        ', 'f + fb -> Z0 + Z0           ', 'f + fb'' -> Z0 + W+/-        ', 'f + fb -> Z0 + H0           ', 'f + fb -> W+ + W-           ', 'f + fb'' -> W+/- + H0        ', 'f + fb -> H0 + H0           ', 'f + g -> f + g              ', 'f + g -> f + gamma          ', 'f + g -> f + Z0             ', 'f + g -> f'' + W+/-          ', 'f + g -> f + H0             ', 'f + gamma -> f + g          ', 'f + gamma -> f + gamma      ', 'f + gamma -> f + Z0         ', 'f + gamma -> f'' + W+/-      ', 'f + gamma -> f + H0         ', 'f + Z0 -> f + g             ', 'f + Z0 -> f + gamma         ', 'f + Z0 -> f + Z0            '/
  Data (proc(i), i=41, 60)/'f + Z0 -> f'' + W+/-         ', 'f + Z0 -> f + H0            ', 'f + W+/- -> f'' + g          ', 'f + W+/- -> f'' + gamma      ', 'f + W+/- -> f'' + Z0         ', 'f + W+/- -> f'' + W+/-       ', 'f + W+/- -> f'' + H0         ', 'f + H0 -> f + g             ', 'f + H0 -> f + gamma         ', 'f + H0 -> f + Z0            ', 'f + H0 -> f'' + W+/-         ', 'f + H0 -> f + H0            ', 'g + g -> f + fb             ', 'g + gamma -> f + fb         ', 'g + Z0 -> f + fb            ', 'g + W+/- -> f + fb''         ', 'g + H0 -> f + fb            ', 'gamma + gamma -> f + fb     ', 'gamma + Z0 -> f + fb        ', 'gamma + W+/- -> f + fb''     '/
  Data (proc(i), i=61, 80)/'gamma + H0 -> f + fb        ', 'Z0 + Z0 -> f + fb           ', 'Z0 + W+/- -> f + fb''        ', 'Z0 + H0 -> f + fb           ', 'W+ + W- -> f + fb           ', 'W+/- + H0 -> f + fb''        ', 'H0 + H0 -> f + fb           ', 'g + g -> g + g              ', 'gamma + gamma -> W+ + W-    ', 'gamma + W+/- -> gamma + W+/-', 'Z0 + Z0 -> Z0 + Z0          ', 'Z0 + Z0 -> W+ + W-          ', 'Z0 + W+/- -> Z0 + W+/-      ', 'Z0 + Z0 -> Z0 + H0          ', 'W+ + W- -> gamma + gamma    ', 'W+ + W- -> Z0 + Z0          ', 'W+/- + W+/- -> W+/- + W+/-  ', 'W+/- + H0 -> W+/- + H0      ', 'H0 + H0 -> H0 + H0          ', '                            '/
  Data (proc(i), i=81, 100)/'q + qb -> Q + QB, massive   ', 'g + g -> Q + QB, massive    ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', 'Elastic scattering          ', 'Single diffractive          ', 'Double diffractive          ', 'Central diffractive         ', 'Low-pT scattering           ', 'Semihard QCD 2 -> 2         ', '                            ', '                            ', '                            ', '                            '/
  Data (proc(i), i=101, 120)/'g + g -> gamma*/Z0          ', 'g + g -> H0                 ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', 'f + fb -> g + H0            ', 'q + g -> q + H0             ', 'g + g -> g + H0             ', 'g + g -> gamma + gamma      ', 'g + g -> gamma + Z0         ', 'g + g -> Z0 + Z0            ', 'g + g -> W+ + W-            ', '                            ', '                            ', '                            '/
  Data (proc(i), i=121, 140)/'g + g -> f + fb + H0        ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            '/
  Data (proc(i), i=141, 160)/'f + fb -> gamma*/Z0/Z''0     ', 'f + fb'' -> H+/-             ', 'f + fb -> R                 ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            '/
  Data (proc(i), i=161, 180)/'f + g -> f'' + H+/-          ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            ', '                            '/
  Data (proc(i), i=181, 200)/20*'                            '/

End Block Data pydata

!*********************************************************************

Subroutine pykcut(mcut)

!...Dummy routine, which the user can replace in order to make cuts on
!...the kinematics on the parton level before the matrix elements are
!...evaluated and the event is generated. The cross-section estimates
!...will automatically take these cuts into account, so the given
!...values are for the allowed phase space region only. MCUT=0 means
!...that the event has passed the cuts, MCUT=1 that it has failed.
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/

  mcut = 0

  Return
End Subroutine pykcut

!*********************************************************************

Subroutine pystfe(kf, x, q2, xpq)

!...This is a dummy routine, where the user can introduce an interface
!...to his own external structure function parametrization.
!...Arguments in:
!...KF : 2212 for p, 211 for pi+; isospin conjugation for n and charge
!...    conjugation for pbar, nbar or pi- is performed by PYSTFU.
!...X : x value.
!...Q2 : Q^2 value.
!...Arguments out:
!...XPQ(-6:6) : x * f(x,Q2), with index according to KF code,
!...    except that gluon is placed in 0. Thus XPQ(0) = xg,
!...    XPQ(1) = xd, XPQ(-1) = xdbar, XPQ(2) = xu, XPQ(-2) = xubar,
!...    XPQ(3) = xs, XPQ(-3) = xsbar, XPQ(4) = xc, XPQ(-4) = xcbar,
!...    XPQ(5) = xb, XPQ(-5) = xbbar, XPQ(6) = xt, XPQ(-6) = xtbar.
!...
!...One such interface, to the Diemos, Ferroni, Longo, Martinelli
!...proton structure functions, already comes with the package. What
!...the user needs here is external files with the three routines
!...FXG160, FXG260 and FXG360 of the authors above, plus the
!...interpolation routine FINT, which is part of the CERN library
!...KERNLIB package. To avoid problems with unresolved external
!...references, the external calls are commented in the current
!...version. To enable this option, remove the C* at the beginning
!...of the relevant lines.
!...
!...Alternatively, the routine can be used as an interface to the
!...structure function evolution program of Tung. This can be achieved
!...by removing C* at the beginning of some of the lines below.
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  Save /ludat1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  Save /ludat2/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
  Save /pypars/
  Dimension xpq(-6:6), xfdflm(9)
  Character chdflm(9)*5, header*40
  Data chdflm/'UPVAL', 'DOVAL', 'GLUON', 'QBAR ', 'UBAR ', 'SBAR ', 'CBAR ', 'BBAR ', 'TBAR '/
  Data header/'Tung evolution package has been invoked'/
  Data init/0/

!...Proton structure functions from Diemoz, Ferroni, Longo, Martinelli.
!...Allowed variable range 10 GeV2 < Q2 < 1E8 GeV2, 5E-5 < x < .95.
  If (mstp(51)>=11 .And. mstp(51)<=13 .And. mstp(52)<=1) Then
    xdflm = max(0.51E-4, x)
    q2dflm = max(10., min(1E8,q2))
    If (mstp(52)==0) q2dflm = 10.
    Do j = 1, 9
      If (mstp(52)==1 .And. j==9) Then
        q2dflm = q2dflm*(40./pmas(6,1))**2
        q2dflm = max(10., min(1E8,q2))
      End If
      xfdflm(j) = 0.
!...Remove C* on following three lines to enable the DFLM options.
!*      IF(MSTP(51).EQ.11) CALL FXG160(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
!*      IF(MSTP(51).EQ.12) CALL FXG260(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
!*      IF(MSTP(51).EQ.13) CALL FXG360(XDFLM,Q2DFLM,CHDFLM(J),XFDFLM(J))
    End Do
    If (x<0.51E-4 .And. abs(parp(51)-1.)>0.01) Then
      cxs = (0.51E-4/x)**(parp(51)-1.)
      Do j = 1, 7
        xfdflm(j) = xfdflm(j)*cxs
      End Do
    End If
    xpq(0) = xfdflm(3)
    xpq(1) = xfdflm(2) + xfdflm(5)
    xpq(2) = xfdflm(1) + xfdflm(5)
    xpq(3) = xfdflm(6)
    xpq(4) = xfdflm(7)
    xpq(5) = xfdflm(8)
    xpq(6) = xfdflm(9)
    xpq(-1) = xfdflm(5)
    xpq(-2) = xfdflm(5)
    xpq(-3) = xfdflm(6)
    xpq(-4) = xfdflm(7)
    xpq(-5) = xfdflm(8)
    xpq(-6) = xfdflm(9)

!...Proton structure function evolution from Wu-Ki Tung: parton
!...distribution functions incorporating heavy quark mass effects.
!...Allowed variable range: PARP(52) < Q < PARP(53); PARP(54) < x < 1.
  Else
    If (init==0) Then
      i1 = 0
      If (mstp(52)==4) i1 = 1
      ihdrn = 1
      nu = mstp(53)
      i2 = mstp(51)
      If (mstp(51)>=11) i2 = mstp(51) - 3
      i3 = 0
      If (mstp(52)==3) i3 = 1

!...Convert to Lambda in CWZ scheme (approximately linear relation).
      alam = 0.75*parp(1)
      tpms = pmas(6, 1)
      qini = parp(52)
      qmax = parp(53)
      xmin = parp(54)

!...Initialize evolution (perform calculation or read results from
!...file).
!...Remove C* on following two lines to enable Tung initialization.
!*        CALL PDFSET(I1,IHDRN,ALAM,TPMS,QINI,QMAX,XMIN,NU,HEADER,
!*   &    I2,I3,IRET,IRR)
      init = 1
    End If

!...Put into output array.
    q = sqrt(q2)
    Do i = -6, 6
      fixq = 0.
!...Remove C* on following line to enable structure function call.
!*      FIXQ=MAX(0.,PDF(10,1,I,X,Q,IR))
      xpq(i) = x*fixq
    End Do

!...Change order of u and d quarks from Tung to PYTHIA convention.
    xps = xpq(1)
    xpq(1) = xpq(2)
    xpq(2) = xps
    xps = xpq(-1)
    xpq(-1) = xpq(-2)
    xpq(-2) = xps
  End If

  Return
End Subroutine pystfe


