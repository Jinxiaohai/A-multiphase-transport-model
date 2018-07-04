!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine artmn
  Parameter (maxstr=150001, maxr=1, amu=0.9383, aka=0.498, etam=0.5475)
  Parameter (maxx=20, maxz=24)
  Parameter (isum=1001, igam=1100)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  !lin      PARAMETER (MAXP = 14000)
  !----------------------------------------------------------------------*
  Integer outpar, zta, zpr
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /hh/proper(maxstr)
  !c      SAVE /HH/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  !c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  !c      SAVE /gg/
  Common /input/nstar, ndirct, dir
  !c      SAVE /INPUT/
  Common /pp/prho(-20:20, -24:24)
  Common /qq/phrho(-maxz:maxz, -24:24)
  Common /rr/massr(0:maxr)
  !c      SAVE /RR/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Common /zz/zta, zpr
  !c      SAVE /zz/
  Common /run/num
  !c      SAVE /RUN/
  !lin-4/2008:
  !      COMMON  /KKK/     TKAON(7),EKAON(7,0:200)
  Common /kkk/tkaon(7), ekaon(7, 0:2000)
  !c      SAVE /KKK/
  Common /kaon/ak(3, 50, 36), speck(50, 36, 7), mf
  !c      SAVE /KAON/
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /ddpi/pirho(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DDpi/
  Common /tt/pel(-maxx:maxx, -maxx:maxx, -maxz:maxz), rxy(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /tt/
  !lin-4/2008:
  !      DIMENSION TEMP(3,MAXSTR),SKAON(7),SEKAON(7,0:200)
  Dimension temp(3, maxstr), skaon(7), sekaon(7, 0:2000)
  !bz12/2/98
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  !c      SAVE /INPUT2/
  Common /input3/plab, elab, zeropt, b0, bi, bm, dencut, cycbox
  !c      SAVE /INPUT3/
  !bz12/2/98end
  !bz11/16/98
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  !c      SAVE /ARPRNT/

  !.....note in the below, since a common block in ART is called EE,
  !.....the variable EE in /ARPRC/is changed to PEAR.
  !lin-9/29/03 changed name in order to distinguish from /prec2/
  !        COMMON /ARPRC/ ITYPAR(MAXSTR),
  !     &       GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
  !     &       PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
  !     &       XMAR(MAXSTR)
  !c      SAVE /ARPRC/
  !lin-9/29/03-end
  Common /arercp/pro1(maxstr, maxr)
  !c      SAVE /ARERCP/
  Common /arerc1/multi1(maxr)
  !c      SAVE /ARERC1/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  !c      SAVE /ARPRC1/
  !
  Dimension npi(maxr)
  Dimension rt(3, maxstr, maxr), pt(3, maxstr, maxr), et(maxstr, maxr), lt(maxstr, maxr), prot(maxstr, maxr)

  External iarflv, invflv
  !bz11/16/98end
  Common /lastt/itimeh, bimp
  !c      SAVE /lastt/
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  !c      SAVE /snn/
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  !c      SAVE /hbt/
  Common /resdcy/nsav, iksdcy
  !c      SAVE /resdcy/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /ftmax/ftsv(maxstr), ftsvt(maxstr, maxr)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  !lin-4/2008 zet() expanded to avoid out-of-bound errors:
  Real zet(-45:45)
  Save
  Data zet/1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., -1., 0., -2., -1., 0., 1., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., -1., 0., 1., 2., 0., 1., 0., 1., 0., -1., 0., 1., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1./

  nlast = 0
  Do i = 1, maxstr
     ftsv(i) = 0.
     Do irun = 1, maxr
        ftsvt(i, irun) = 0.
     End Do
     lblast(i) = 999
     Do j = 1, 4
        !lin-4/2008 bugs pointed out by Vander Molen & Westfall:
        !            xlast(i,j)=0.
        !            plast(i,j)=0.
        xlast(j, i) = 0.
        plast(j, i) = 0.
     End Do
  End Do

  Call tablem
  ! several control parameters, keep them fixed in this code.
  ikaon = 1
  nstar = 1
  ndirct = 0
  dir = 0.02
  asy = 0.032
  esbin = 0.04
  mf = 36
  !----------------------------------------------------------------------*
  !      CALL FRONT(12,MASSTA,MASSPR,ELAB)
  !----------------------------------------------------------------------*
  radta = 1.124*float(massta)**(1./3.)
  radpr = 1.124*float(masspr)**(1./3.)
  zdist = radta + radpr
  !      if ( cycbox.ne.0 ) zdist=0
  bmax = radta + radpr
  mass = massta + masspr
  ntotal = num*mass
  !
  If (ntotal>maxstr) Then
     Write (12, '(//10X,''**** FATAL ERROR: TOO MANY TEST PART. ****'//' '')')
     Stop
  End If
  !
  !-----------------------------------------------------------------------
  !       RELATIVISTIC KINEMATICS
  !
  !       1) LABSYSTEM
  !
  eta = float(massta)*amu
  pzta = 0.0
  betata = 0.0
  gammta = 1.0
  !
  epr = float(masspr)*(amu+0.001*elab)
  pzpr = sqrt(epr**2-(amu*float(masspr))**2)
  betapr = pzpr/epr
  gammpr = 1.0/sqrt(1.0-betapr**2)
  !
  ! BETAC AND GAMMAC OF THE C.M. OBSERVED IN THE LAB. FRAME
  betac = (pzpr+pzta)/(epr+eta)
  gammc = 1.0/sqrt(1.-betac**2)
  If (insys/=0) Then
     !
     !       2) C.M. SYSTEM
     !
     s = (epr+eta)**2 - pzpr**2
     xx1 = 4.*alog(float(massta))
     xx2 = 4.*alog(float(masspr))
     xx1 = exp(xx1)
     xx2 = exp(xx2)
     psqare = (s**2+(xx1+xx2)*amu**4-2.0*s*amu**2*float(massta**2+masspr**2)-2.0*float(massta**2*masspr**2)*amu**4)/(4.0*s)
     !
     eta = sqrt(psqare+(float(massta)*amu)**2)
     pzta = -sqrt(psqare)
     betata = pzta/eta
     gammta = 1.0/sqrt(1.0-betata**2)
     !
     epr = sqrt(psqare+(float(masspr)*amu)**2)
     pzpr = sqrt(psqare)
     betapr = pzpr/epr
     gammpr = 1.0/sqrt(1.0-betapr**2)
     !
  Else
     !        WRITE(12,'(/10x,''*** CALCULATION DONE IN LAB-FRAME ***''/)')
  End If
  ! MOMENTUM PER PARTICLE
  pzta = pzta/float(massta)
  pzpr = pzpr/float(masspr)
  ! total initial energy in the N-N cms frame
  ecms0 = eta + epr
  Do imany = 1, manyb
     !------------------------------------------------------------------------
     ! Initialize the impact parameter B
     If (manyb>1) Then
111     bx = 1.0 - 2.0*ranart(nseed)
        by = 1.0 - 2.0*ranart(nseed)
        b2 = bx*bx + by*by
        If (b2>1.0) Goto 111
        b = sqrt(b2)*(bm-bi) + bi
     Else
        b = b0
     End If
     !1 INITIALIZATION IN ISOSPIN SPACE FOR BOTH THE PROJECTILE AND TARGET
     Call coulin(masspr, massta, num)
     !2 INITIALIZATION IN PHASE SPACE FOR THE TARGET
     Call init(1, massta, num, radta, b/2., zeropt+zdist/2., pzta, gammta, iseed, mass, imomen)
     !3.1 INITIALIZATION IN PHASE SPACE FOR THE PROJECTILE
     Call init(1+massta, mass, num, radpr, -b/2., zeropt-zdist/2., pzpr, gammpr, iseed, mass, imomen)
     !3.2 OUTPAR IS THE NO. OF ESCAPED PARTICLES
     outpar = 0
     massr(0) = 0
     Do ir = 1, num
        massr(ir) = mass
     End Do
     !3.4 INITIALIZation FOR THE KAON SPECTRUM
     !      CALL KSPEC0(BETAC,GAMMC)
     ! calculate the local baryon density matrix
     Call dens(ipot, mass, num, outpar)
     If (icoll/=-1) Then
        Do i = 1, ntotal
           ix = nint(r(1,i))
           iy = nint(r(2,i))
           iz = nint(r(3,i))
           !lin-4/2008 check bounds:
           If (ix>=maxx .Or. iy>=maxx .Or. iz>=maxz .Or. ix<=-maxx .Or. iy<=-maxx .Or. iz<=-maxz) Goto 700
           Call gradu(ipot, ix, iy, iz, gradx, grady, gradz)
           p(1, i) = p(1, i) - (0.5*dt)*gradx
           p(2, i) = p(2, i) - (0.5*dt)*grady
           p(3, i) = p(3, i) - (0.5*dt)*gradz
700     End Do
     End If
     rcnne = 0
     rdd = 0
     rpp = 0
     rppk = 0
     rpn = 0
     rpd = 0
     rkn = 0
     rnnk = 0
     rddk = 0
     rndk = 0
     rcnnd = 0
     rcndn = 0
     rcoll = 0
     rbloc = 0
     rdirt = 0
     rdecay = 0
     rres = 0
     !4.11 KAON PRODUCTION PROBABILITY COUNTER FOR PERTURBATIVE CALCULATIONS ONLY
     Do kkk = 1, 5
        skaon(kkk) = 0
        Do is = 1, 2000
           sekaon(kkk, is) = 0
        End Do
     End Do
     !4.12 anti-proton and anti-kaon counters
     pr0 = 0.
     pr1 = 0.
     ska0 = 0.
     ska1 = 0.
     !       ============== LOOP OVER ALL TIME STEPS ================       *
     !                             STARTS HERE                              *
     !       ========================================================       *
     !bz11/16/98
     If (iapar2(1)/=1) Then
        Do i = 1, maxstr
           Do j = 1, 3
              r(j, i) = 0.
              p(j, i) = 0.
           End Do
           e(i) = 0.
           lb(i) = 0
           !bz3/25/00
           id(i) = 0
           !     sp 12/19/00
           proper(i) = 1.
        End Do
        mass = 0
        !bz12/22/98
        !         MASSR(1) = 0
        !         NP = 0
        !         NPI = 1
        np = 0
        Do j = 1, num
           massr(j) = 0
           npi(j) = 1
        End Do
        Do i = 1, maxr
           Do j = 1, maxstr
              rt(1, j, i) = 0.
              rt(2, j, i) = 0.
              rt(3, j, i) = 0.
              pt(1, j, i) = 0.
              pt(2, j, i) = 0.
              pt(3, j, i) = 0.
              et(j, i) = 0.
              lt(j, i) = 0
              !     sp 12/19/00
              prot(j, i) = 1.
           End Do
        End Do
        !bz12/22/98end
     End If
     !bz11/16/98end

     Do nt = 1, ntmax
        !TEMPORARY PARTICLE COUNTERS
        !4.2 PION COUNTERS : LP1,LP2 AND LP3 ARE THE NO. OF P+,P0 AND P-
        lp1 = 0
        lp2 = 0
        lp3 = 0
        !4.3 DELTA COUNTERS : LD1,LD2,LD3 AND LD4 ARE THE NO. OF D++,D+,D0 AND D-
        ld1 = 0
        ld2 = 0
        ld3 = 0
        ld4 = 0
        !4.4 N*(1440) COUNTERS : LN1 AND LN2 ARE THE NO. OF N*+ AND N*0
        ln1 = 0
        ln2 = 0
        !4.5 N*(1535) counters
        ln5 = 0
        !4.6 ETA COUNTERS
        le = 0
        !4.7 KAON COUNTERS
        lkaon = 0

        !lin-11/09/00:
        ! KAON* COUNTERS
        lkaons = 0

        !-----------------------------------------------------------------------
        If (icoll/=1) Then
           ! STUDYING BINARY COLLISIONS AMONG PARTICLES DURING THIS TIME INTERVAL *
           !lin-10/25/02 get rid of argument usage mismatch in relcol(.nt.):
           numnt = nt
           Call relcol(lcoll, lbloc, lcnne, ldd, lpp, lppk, lpn, lpd, lrho, lomega, lkn, lnnk, lddk, lndk, lcnnd, lcndn, ldirt, ldecay, lres, ldou, lddrho, lnnrho, lnnom, numnt, ntmax, sp, akaon, sk)
           rcoll = rcoll + float(lcoll)/num
           rbloc = rbloc + float(lbloc)/num
           rcnne = rcnne + float(lcnne)/num
           rdd = rdd + float(ldd)/num
           rpp = rpp + float(lpp)/num
           rppk = rppk + float(lppk)/num
           rpn = rpn + float(lpn)/num
           rpd = rpd + float(lpd)/num
           rkn = rkn + float(lkn)/num
           rnnk = rnnk + float(lnnk)/num
           rddk = rddk + float(lddk)/num
           rndk = rndk + float(lndk)/num
           rcnnd = rcnnd + float(lcnnd)/num
           rcndn = rcndn + float(lcndn)/num
           rdirt = rdirt + float(ldirt)/num
           rdecay = rdecay + float(ldecay)/num
           rres = rres + float(lres)/num
           ! AVERAGE RATES OF VARIOUS COLLISIONS IN THE CURRENT TIME STEP
           adirt = ldirt/dt/num
           acoll = (lcoll-lbloc)/dt/num
           acnnd = lcnnd/dt/num
           acndn = lcndn/dt/num
           adecay = ldecay/dt/num
           ares = lres/dt/num
           adou = ldou/dt/num
           addrho = lddrho/dt/num
           annrho = lnnrho/dt/num
           annom = lnnom/dt/num
           add = ldd/dt/num
           app = lpp/dt/num
           appk = lppk/dt/num
           apn = lpn/dt/num
           apd = lpd/dt/num
           arh = lrho/dt/num
           aom = lomega/dt/num
           akn = lkn/dt/num
           annk = lnnk/dt/num
           addk = lddk/dt/num
           andk = lndk/dt/num

        End If
        !
        !       UPDATE BARYON DENSITY
        !
        Call dens(ipot, mass, num, outpar)
        !
        !       UPDATE POSITIONS FOR ALL THE PARTICLES PRESENT AT THIS TIME
        !
        sumene = 0
        iso = 0
        Do mrun = 1, num
           iso = iso + massr(mrun-1)
           Do i0 = 1, massr(mrun)
              i = i0 + iso
              etotal = sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
              sumene = sumene + etotal
              ! for kaons, if there is a potential
              ! CALCULATE THE ENERGY OF THE KAON ACCORDING TO THE IMPULSE APPROXIMATION
              ! REFERENCE: B.A. LI AND C.M. KO, PHYS. REV. C 54 (1996) 3283.
              If (kpoten/=0 .And. lb(i)==23) Then
                 den = 0.
                 ix = nint(r(1,i))
                 iy = nint(r(2,i))
                 iz = nint(r(3,i))
                 !lin-4/2008:
                 !       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
                 !     & ABS(IZ) .LT. MAXZ) den=rho(ix,iy,iz)
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) den = rho(ix, iy, iz)
                 !         ecor=0.1973**2*0.255*kmul*4*3.14159*(1.+0.4396/0.938)
                 !         etotal=sqrt(etotal**2+ecor*den)
                 !** G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV, m^*=m
                 !     GeV^2 fm^3
                 akg = 0.1727
                 !     GeV fm^3
                 bkg = 0.333
                 rnsg = den
                 ecor = -akg*rnsg + (bkg*den)**2
                 etotal = sqrt(etotal**2+ecor)
              End If
              !
              If (kpoten/=0 .And. lb(i)==21) Then
                 den = 0.
                 ix = nint(r(1,i))
                 iy = nint(r(2,i))
                 iz = nint(r(3,i))
                 !lin-4/2008:
                 !       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
                 !     & ABS(IZ) .LT. MAXZ) den=rho(ix,iy,iz)
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) den = rho(ix, iy, iz)
                 !* for song potential no effect on position
                 !** G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV, m^*=m
                 !     GeV^2 fm^3
                 akg = 0.1727
                 !     GeV fm^3
                 bkg = 0.333
                 rnsg = den
                 ecor = -akg*rnsg + (bkg*den)**2
                 etotal = sqrt(etotal**2+ecor)
              End If
              !
              ! UPDATE POSITIONS
              r(1, i) = r(1, i) + dt*p(1, i)/etotal
              r(2, i) = r(2, i) + dt*p(2, i)/etotal
              r(3, i) = r(3, i) + dt*p(3, i)/etotal
              ! use cyclic boundary conitions
              If (cycbox/=0) Then
                 If (r(1,i)>cycbox/2) r(1, i) = r(1, i) - cycbox
                 If (r(1,i)<=-cycbox/2) r(1, i) = r(1, i) + cycbox
                 If (r(2,i)>cycbox/2) r(2, i) = r(2, i) - cycbox
                 If (r(2,i)<=-cycbox/2) r(2, i) = r(2, i) + cycbox
                 If (r(3,i)>cycbox/2) r(3, i) = r(3, i) - cycbox
                 If (r(3,i)<=-cycbox/2) r(3, i) = r(3, i) + cycbox
              End If
              ! UPDATE THE DELTA, N* AND PION COUNTERS
              lb1 = lb(i)
              ! 1. FOR DELTA++
              If (lb1==9) ld1 = ld1 + 1
              ! 2. FOR DELTA+
              If (lb1==8) ld2 = ld2 + 1
              ! 3. FOR DELTA0
              If (lb1==7) ld3 = ld3 + 1
              ! 4. FOR DELTA-
              If (lb1==6) ld4 = ld4 + 1
              ! 5. FOR N*+(1440)
              If (lb1==11) ln1 = ln1 + 1
              ! 6. FOR N*0(1440)
              If (lb1==10) ln2 = ln2 + 1
              ! 6.1 FOR N*(1535)
              If ((lb1==13) .Or. (lb1==12)) ln5 = ln5 + 1
              ! 6.2 FOR ETA
              If (lb1==0) le = le + 1
              ! 6.3 FOR KAONS
              If (lb1==23) lkaon = lkaon + 1
              !lin-11/09/00: FOR KAON*
              If (lb1==30) lkaons = lkaons + 1

              ! UPDATE PION COUNTER
              ! 7. FOR PION+
              If (lb1==5) lp1 = lp1 + 1
              ! 8. FOR PION0
              If (lb1==4) lp2 = lp2 + 1
              ! 9. FOR PION-
              If (lb1==3) lp3 = lp3 + 1
           End Do
        End Do
        lp = lp1 + lp2 + lp3
        ld = ld1 + ld2 + ld3 + ld4
        ln = ln1 + ln2
        alp = float(lp)/float(num)
        ald = float(ld)/float(num)
        aln = float(ln)/float(num)
        aln5 = float(ln5)/float(num)
        atotal = alp + ald + aln + 0.5*aln5
        ale = float(le)/float(num)
        alkaon = float(lkaon)/float(num)
        ! UPDATE MOMENTUM DUE TO COULOMB INTERACTION
        If (icou==1) Then
           !       with Coulomb interaction
           iso = 0
           Do irun = 1, num
              iso = iso + massr(irun-1)
              Do il = 1, massr(irun)
                 temp(1, il) = 0.
                 temp(2, il) = 0.
                 temp(3, il) = 0.
              End Do
              Do il = 1, massr(irun)
                 i = iso + il
                 If (zet(lb(i))/=0) Then
                    Do jl = 1, il - 1
                       j = iso + jl
                       If (zet(lb(j))/=0) Then
                          ddx = r(1, i) - r(1, j)
                          ddy = r(2, i) - r(2, j)
                          ddz = r(3, i) - r(3, j)
                          rdiff = sqrt(ddx**2+ddy**2+ddz**2)
                          If (rdiff<=1.) rdiff = 1.
                          grp = zet(lb(i))*zet(lb(j))/rdiff**3
                          ddx = ddx*grp
                          ddy = ddy*grp
                          ddz = ddz*grp
                          temp(1, il) = temp(1, il) + ddx
                          temp(2, il) = temp(2, il) + ddy
                          temp(3, il) = temp(3, il) + ddz
                          temp(1, jl) = temp(1, jl) - ddx
                          temp(2, jl) = temp(2, jl) - ddy
                          temp(3, jl) = temp(3, jl) - ddz
                       End If
                    End Do
                 End If
              End Do
              Do il = 1, massr(irun)
                 i = iso + il
                 If (zet(lb(i))/=0) Then
                    Do idir = 1, 3
                       p(idir, i) = p(idir, i) + temp(idir, il)*dt*0.00144
                    End Do
                 End If
              End Do
           End Do
        End If
        !       In the following, we shall:
        !       (1) UPDATE MOMENTA DUE TO THE MEAN FIELD FOR BARYONS AND KAONS,
        !       (2) calculate the thermalization, temperature in a sphere of
        !           radius 2.0 fm AROUND THE CM
        !       (3) AND CALCULATE THE NUMBER OF PARTICLES IN THE HIGH DENSITY REGION
        spt = 0
        spz = 0
        ncen = 0
        ekin = 0
        nlost = 0
        mean = 0
        nquark = 0
        nbaryn = 0
        !sp06/18/01
        rads = 2.
        zras = 0.1
        denst = 0.
        edenst = 0.
        !sp06/18/01 end
        Do irun = 1, num
           mean = mean + massr(irun-1)
           Do j = 1, massr(irun)
              i = j + mean
              !
              !sp06/18/01
              radut = sqrt(r(1,i)**2+r(2,i)**2)
              If (radut<=rads) Then
                 If (abs(r(3,i))<=zras*nt*dt) Then
                    !         vols = 3.14159*radut**2*abs(r(3,i))      ! cylinder pi*r^2*l
                    !     cylinder pi*r^2*l
                    vols = 3.14159*rads**2*zras
                    engs = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
                    gammas = 1.
                    If (e(i)/=0.) gammas = engs/e(i)
                    !     rho
                    denst = denst + 1./gammas/vols
                    !     energy density
                    edenst = edenst + engs/gammas/gammas/vols
                 End If
              End If
              !sp06/18/01 end
              !
              drr = sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
              If (drr<=2.0) Then
                 spt = spt + p(1, i)**2 + p(2, i)**2
                 spz = spz + p(3, i)**2
                 ncen = ncen + 1
                 ekin = ekin + sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2) - e(i)
              End If
              ix = nint(r(1,i))
              iy = nint(r(2,i))
              iz = nint(r(3,i))
              ! calculate the No. of particles in the high density region
              !lin-4/2008:
              !              IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
              !     & ABS(IZ) .LT. MAXZ) THEN
              If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                 If (rho(ix,iy,iz)/0.168>dencut) Goto 5800
                 If ((rho(ix,iy,iz)/0.168>5.) .And. (e(i)>0.9)) nbaryn = nbaryn + 1
                 If (pel(ix,iy,iz)>2.0) nquark = nquark + 1
              End If
              !*
              ! If there is a kaon potential, propogating kaons
              If (kpoten/=0 .And. lb(i)==23) Then
                 den = 0.
                 !lin-4/2008:
                 !       IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
                 !     & ABS(IZ) .LT. MAXZ)then
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                    den = rho(ix, iy, iz)
                    !        ecor=0.1973**2*0.255*kmul*4*3.14159*(1.+0.4396/0.938)
                    !       etotal=sqrt(P(1,i)**2+p(2,I)**2+p(3,i)**2+e(i)**2+ecor*den)
                    !** for G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV
                    !     !! GeV^2 fm^3
                    akg = 0.1727
                    !     !! GeV fm^3
                    bkg = 0.333
                    rnsg = den
                    ecor = -akg*rnsg + (bkg*den)**2
                    etotal = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2+ecor)
                    ecor = -akg + 2.*bkg**2*den + 2.*bkg*etotal
                    !** G.Q. Li potential (END)
                    Call graduk(ix, iy, iz, gradxk, gradyk, gradzk)
                    p(1, i) = p(1, i) - dt*gradxk*ecor/(2.*etotal)
                    p(2, i) = p(2, i) - dt*gradyk*ecor/(2.*etotal)
                    p(3, i) = p(3, i) - dt*gradzk*ecor/(2.*etotal)
                 End If
              End If
              !
              If (kpoten/=0 .And. lb(i)==21) Then
                 den = 0.
                 !lin-4/2008:
                 !           IF (ABS(IX) .LT. MAXX .AND. ABS(IY) .LT. MAXX .AND.
                 !     &        ABS(IZ) .LT. MAXZ)then
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                    den = rho(ix, iy, iz)
                    Call graduk(ix, iy, iz, gradxk, gradyk, gradzk)
                    !        P(1,I) = P(1,I) - DT * GRADXk*(-0.12/0.168)    !! song potential
                    !        P(2,I) = P(2,I) - DT * GRADYk*(-0.12/0.168)
                    !        P(3,I) = P(3,I) - DT * GRADZk*(-0.12/0.168)
                    !** for G.Q Li potential form with n_s = n_b and pot(n_0)=29 MeV
                    !    !! GeV^2 fm^3
                    akg = 0.1727
                    !     !! GeV fm^3
                    bkg = 0.333
                    rnsg = den
                    ecor = -akg*rnsg + (bkg*den)**2
                    etotal = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2+ecor)
                    ecor = -akg + 2.*bkg**2*den - 2.*bkg*etotal
                    p(1, i) = p(1, i) - dt*gradxk*ecor/(2.*etotal)
                    p(2, i) = p(2, i) - dt*gradyk*ecor/(2.*etotal)
                    p(3, i) = p(3, i) - dt*gradzk*ecor/(2.*etotal)
                    !** G.Q. Li potential (END)
                 End If
              End If
              !
              ! for other mesons, there is no potential
              If (j>mass) Goto 5800
              !  with mean field interaction for baryons   (open endif below) !!sp05
              !*      if( (iabs(lb(i)).eq.1.or.iabs(lb(i)).eq.2) .or.
              !*    &     (iabs(lb(i)).ge.6.and.iabs(lb(i)).le.17) .or.
              !*    &      iabs(lb(i)).eq.40.or.iabs(lb(i)).eq.41 )then
              If (icoll/=-1) Then
                 If (ix<maxx .And. iy<maxx .And. iz<maxz .And. ix>-maxx .And. iy>-maxx .And. iz>-maxz) Then
                    Call gradu(ipot, ix, iy, iz, gradx, grady, gradz)
                    tz = 0.
                    gradxn = 0
                    gradyn = 0
                    gradzn = 0
                    gradxp = 0
                    gradyp = 0
                    gradzp = 0
                    If (icou==1) Then
                       Call gradup(ix, iy, iz, gradxp, gradyp, gradzp)
                       Call gradun(ix, iy, iz, gradxn, gradyn, gradzn)
                       If (zet(lb(i))/=0) tz = -1
                       If (zet(lb(i))==0) tz = 1
                    End If
                    If (iabs(lb(i))>=14 .And. iabs(lb(i))<=17) Then
                       facl = 2./3.
                    Else If (iabs(lb(i))==40 .Or. iabs(lb(i))==41) Then
                       facl = 1./3.
                    Else
                       facl = 1.
                    End If
                    p(1, i) = p(1, i) - facl*dt*(gradx+asy*(gradxn-gradxp)*tz)
                    p(2, i) = p(2, i) - facl*dt*(grady+asy*(gradyn-gradyp)*tz)
                    p(3, i) = p(3, i) - facl*dt*(gradz+asy*(gradzn-gradzp)*tz)
                 End If
              End If
              !*          endif          !!sp05
5800       End Do
        End Do
        cden = rho(0, 0, 0)/0.168
        If ((nt/nfreq)*nfreq==nt) Then
           If (icflow==1) Call flow(nt)
        End If
        If (iapar2(1)/=1) Then
           ct = nt*dt
           ia = 0
           Do irun = 1, num
              Do ic = 1, massr(irun)
                 ie = ia + ic
                 rt(1, ic, irun) = r(1, ie)
                 rt(2, ic, irun) = r(2, ie)
                 rt(3, ic, irun) = r(3, ie)
                 pt(1, ic, irun) = p(1, ie)
                 pt(2, ic, irun) = p(2, ie)
                 pt(3, ic, irun) = p(3, ie)
                 et(ic, irun) = e(ie)
                 lt(ic, irun) = lb(ie)
                 !         !! sp 12/19/00
                 prot(ic, irun) = proper(ie)
                 !lin-5/2008:
                 dpertt(ic, irun) = dpertp(ie)
              End Do
              np = massr(irun)
              np1 = npi(irun)
              ctlong = ct
              If (nt==(ntmax-1)) Then
                 ctlong = 1.E30
              Else If (nt==ntmax) Then
                 Goto 1111
              End If
              !
              Do While (np1<=multi1(irun) .And. ft1(np1,irun)>((nt-1)*dt) .And. ft1(np1,irun)<=ctlong)
                 np = np + 1
                 udt = (ct-ft1(np1,irun))/ee1(np1, irun)
                 If (nt==(ntmax-1)) Then
                    ftsvt(np, irun) = ft1(np1, irun)
                    If (ft1(np1,irun)>ct) udt = 0.
                 End If
                 rt(1, np, irun) = gx1(np1, irun) + px1(np1, irun)*udt
                 rt(2, np, irun) = gy1(np1, irun) + py1(np1, irun)*udt
                 rt(3, np, irun) = gz1(np1, irun) + pz1(np1, irun)*udt
                 pt(1, np, irun) = px1(np1, irun)
                 pt(2, np, irun) = py1(np1, irun)
                 pt(3, np, irun) = pz1(np1, irun)
                 et(np, irun) = xm1(np1, irun)
                 lt(np, irun) = iarflv(ityp1(np1,irun))
                 !lin-5/2008:
                 dpertt(np, irun) = dpp1(np1, irun)
                 np1 = np1 + 1
                 !     !! sp 12/19/00
                 prot(np, irun) = 1.
              End Do
              !
1111          Continue
              npi(irun) = np1
              ia = ia + massr(irun)
              massr(irun) = np
           End Do
           ia = 0
           Do irun = 1, num
              ia = ia + massr(irun-1)
              Do ic = 1, massr(irun)
                 ie = ia + ic
                 r(1, ie) = rt(1, ic, irun)
                 r(2, ie) = rt(2, ic, irun)
                 r(3, ie) = rt(3, ic, irun)
                 p(1, ie) = pt(1, ic, irun)
                 p(2, ie) = pt(2, ic, irun)
                 p(3, ie) = pt(3, ic, irun)
                 e(ie) = et(ic, irun)
                 lb(ie) = lt(ic, irun)
                 !     !! sp 12/19/00
                 proper(ie) = prot(ic, irun)
                 If (nt==(ntmax-1)) ftsv(ie) = ftsvt(ic, irun)
                 !lin-5/2008:
                 dpertp(ie) = dpertt(ic, irun)
              End Do
              !lin-3/2009 Moved here to better take care of freezeout spacetime:
              Call hbtout(massr(irun), nt, ntmax)
           End Do
           !bz12/22/98end
        End If
        !bz11/16/98end

        !lin-5/2009 ctest off:
        !      call flowh(ct)

     End Do
     iss = 0
     Do lrun = 1, num
        iss = iss + massr(lrun-1)
        Do l0 = 1, massr(lrun)
           ipart = iss + l0
        End Do
     End Do

     !bz11/16/98
     If (iapar2(1)/=1) Then
        ia = 0
        Do irun = 1, num
           ia = ia + massr(irun-1)
           np1 = npi(irun)
           nsh = massr(irun) - np1 + 1
           multi1(irun) = multi1(irun) + nsh
           !.....to shift the unformed particles to the end of the common block
           If (nsh>0) Then
              ib = multi1(irun)
              ie = massr(irun) + 1
              ii = -1
           Else If (nsh<0) Then
              ib = massr(irun) + 1
              ie = multi1(irun)
              ii = 1
           End If
           If (nsh/=0) Then
              Do i = ib, ie, ii
                 j = i - nsh
                 ityp1(i, irun) = ityp1(j, irun)
                 gx1(i, irun) = gx1(j, irun)
                 gy1(i, irun) = gy1(j, irun)
                 gz1(i, irun) = gz1(j, irun)
                 ft1(i, irun) = ft1(j, irun)
                 px1(i, irun) = px1(j, irun)
                 py1(i, irun) = py1(j, irun)
                 pz1(i, irun) = pz1(j, irun)
                 ee1(i, irun) = ee1(j, irun)
                 xm1(i, irun) = xm1(j, irun)
                 !     !! sp 12/19/00
                 pro1(i, irun) = pro1(j, irun)
                 !lin-5/2008:
                 dpp1(i, irun) = dpp1(j, irun)
              End Do
           End If

           !.....to copy ART particle info to COMMON /ARPRC1/
           Do i = 1, massr(irun)
              ib = ia + i
              ityp1(i, irun) = invflv(lb(ib))
              gx1(i, irun) = r(1, ib)
              gy1(i, irun) = r(2, ib)
              gz1(i, irun) = r(3, ib)
              !lin-10/28/03:
              ! since all unformed hadrons at time ct are read in at nt=ntmax-1,
              ! their formation time ft1 should be kept to determine their freezeout(x,t):
              !              FT1(I, IRUN) = CT
              If (ft1(i,irun)<ct) ft1(i, irun) = ct
              px1(i, irun) = p(1, ib)
              py1(i, irun) = p(2, ib)
              pz1(i, irun) = p(3, ib)
              xm1(i, irun) = e(ib)
              ee1(i, irun) = sqrt(px1(i,irun)**2+py1(i,irun)**2+pz1(i,irun)**2+xm1(i,irun)**2)
              !     !! sp 12/19/00
              pro1(i, irun) = proper(ib)
           End Do
        End Do
        !bz12/22/98end
     End If
  End Do
  Return
  !bz11/16/98end
End Subroutine artmn
!*********************************
Subroutine coulin(masspr, massta, num)
  Integer zta, zpr
  Parameter (maxstr=150001)
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /zz/zta, zpr
  !c      SAVE /zz/
  Save
  mass = massta + masspr
  Do irun = 1, num
     Do i = 1 + (irun-1)*mass, zta + (irun-1)*mass
        lb(i) = 1
     End Do
     Do i = zta + 1 + (irun-1)*mass, massta + (irun-1)*mass
        lb(i) = 2
     End Do
     Do i = massta + 1 + (irun-1)*mass, massta + zpr + (irun-1)*mass
        lb(i) = 1
     End Do
     Do i = massta + zpr + 1 + (irun-1)*mass, massta + masspr + (irun-1)*mass
        lb(i) = 2
     End Do
  End Do
  Return
End Subroutine coulin
!*********************************
!

!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine relcol(lcoll, lbloc, lcnne, ldd, lpp, lppk, lpn, lpd, lrho, lomega, lkn, lnnk, lddk, lndk, lcnnd, lcndn, ldirt, ldecay, lres, ldou, lddrho, lnnrho, lnnom, nt, ntmax, sp, akaon, sk)
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, aks=0.895)
  Parameter (aa1=1.26, aphi=1.02, ap1=0.13496)
  Parameter (maxx=20, maxz=24)
  Parameter (rrkk=0.6, prkk=0.3, srhoks=5., esbin=0.04)
  Dimension massrn(0:maxr), rt(3, maxstr), pt(3, maxstr), et(maxstr)
  Dimension lt(maxstr), prot(maxstr)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /hh/proper(maxstr)
  !c      SAVE /HH/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  !c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  !c      SAVE /gg/
  Common /input/nstar, ndirct, dir
  !c      SAVE /INPUT/
  Common /nn/nnn
  !c      SAVE /NN/
  Common /rr/massr(0:maxr)
  !c      SAVE /RR/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Common /bg/betax, betay, betaz, gamma
  !c      SAVE /BG/
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
  Common /pe/propi(maxstr, maxr)
  !c      SAVE /PE/
  Common /kkk/tkaon(7), ekaon(7, 0:2000)
  !c      SAVE /KKK/
  Common /kaon/ak(3, 50, 36), speck(50, 36, 7), mf
  !c      SAVE /KAON/
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  !c      SAVE /leadng/
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  !c      SAVE /tdecay/
  Common /lastt/itimeh, bimp
  !c      SAVE /lastt/
  !
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
  !c      SAVE /ppbmas/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  !c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  !c      SAVE /ppmm/
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  !c      SAVE /hbt/
  Common /resdcy/nsav, iksdcy
  !c      SAVE /resdcy/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /ftmax/ftsv(maxstr), ftsvt(maxstr, maxr)
  Dimension ftpisv(maxstr, maxr), fttemp(maxstr)
  Common /dpi/em2, lb2
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  !lin-5/2008:
  Dimension dptemp(maxstr)
  Common /para8/idpert, npertd, idxsec
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  !
  Real zet(-45:45)
  Save
  Data zet/1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., -1., 0., -2., -1., 0., 1., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., -1., 0., 1., 2., 0., 1., 0., 1., 0., -1., 0., 1., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1./

  Call inidcy
  resona = 5.
  !-----------------------------------------------------------------------
  !     INITIALIZATION OF COUNTING VARIABLES
  nodelt = 0
  sumsrt = 0.
  lcoll = 0
  lbloc = 0
  lcnne = 0
  ldd = 0
  lpp = 0
  lpd = 0
  lpdr = 0
  lrho = 0
  lrhor = 0
  lomega = 0
  lomgar = 0
  lpn = 0
  lkn = 0
  lnnk = 0
  lddk = 0
  lndk = 0
  lppk = 0
  lcnnd = 0
  lcndn = 0
  ldirt = 0
  ldecay = 0
  lres = 0
  ldou = 0
  lddrho = 0
  lnnrho = 0
  lnnom = 0
  msum = 0
  massrn(0) = 0
  ! COM: MSUM IS USED TO COUNT THE TOTAL NO. OF PARTICLES
  !      IN PREVIOUS IRUN-1 RUNS
  ! KAON COUNTERS
  Do il = 1, 5
     tkaon(il) = 0
     Do is = 1, 2000
        ekaon(il, is) = 0
     End Do
  End Do
  ! sp 12/19/00
  Do i = 1, num
     Do j = 1, maxstr
        propi(j, i) = 1.
     End Do
  End Do

  Do i = 1, maxstr
     fttemp(i) = 0.
     Do irun = 1, maxr
        ftpisv(i, irun) = 0.
     End Do
  End Do

  ! sp 12/19/00 end
  sp = 0
  ! antikaon counters
  akaon = 0
  sk = 0
  !-----------------------------------------------------------------------
  !     LOOP OVER ALL PARALLEL RUNS
  !bz11/17/98
  !      MASS=MASSPR+MASSTA
  mass = 0
  !bz11/17/98end
  Do irun = 1, num
     nnn = 0
     msum = msum + massr(irun-1)
     !     LOOP OVER ALL PSEUDOPARTICLES 1 IN THE SAME RUN
     j10 = 2
     If (nt==ntmax) j10 = 1
     !
     !test off skips the check of energy conservation after each timestep:
     !         enetot=0.
     !         do ip=1,MASSR(IRUN)
     !            if(e(ip).ne.0.or.lb(ip).eq.10022) enetot=enetot
     !     1           +sqrt(p(1,ip)**2+p(2,ip)**2+p(3,ip)**2+e(ip)**2)
     !         enddo
     !         write(91,*) 'A:',nt,enetot,massr(irun),bimp

     Do j1 = j10, massr(irun)
        i1 = j1 + msum
        ! E(I)=0 are for pions having been absorbed or photons which do not enter here:
        !lin-4/2012 option of pi0 decays:
        !            IF(E(I1).EQ.0.)GO TO 800
        If (e(i1)==0.) Goto 798
        !     To include anti-(Delta,N*1440 and N*1535):
        !          IF ((LB(I1) .LT. -13 .OR. LB(I1) .GT. 28)
        !     1         .and.iabs(LB(I1)) .ne. 30 ) GOTO 800
        !lin-4/2012 option of pi0 decays:
        !            IF (LB(I1) .LT. -45 .OR. LB(I1) .GT. 45) GOTO 800
        If (lb(i1)<-45 .Or. lb(i1)>45) Goto 798
        x1 = r(1, i1)
        y1 = r(2, i1)
        z1 = r(3, i1)
        px1 = p(1, i1)
        py1 = p(2, i1)
        pz1 = p(3, i1)
        em1 = e(i1)
        am1 = em1
        e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
        id1 = id(i1)
        lb1 = lb(i1)

        !     generate k0short and k0long from K+ and K- at the last timestep:
        If (nt==ntmax .And. (lb1==21 .Or. lb1==23)) Then
           pk0 = ranart(nseed)
           If (pk0<0.25) Then
              lb(i1) = 22
           Else If (pk0<0.50) Then
              lb(i1) = 24
           End If
           lb1 = lb(i1)
        End If

        !lin-8/07/02 these particles don't decay strongly, so skip decay routines:
        !            IF( (lb1.ge.-2.and.lb1.le.5) .OR. lb1.eq.31 .OR.
        !     &           (iabs(lb1).ge.14.and.iabs(lb1).le.24) .OR.
        !     &           (iabs(lb1).ge.40.and.iabs(lb1).le.45) .or.
        !     &           lb1.eq.31)GO TO 1
        !     only decay K0short when iksdcy=1:
        If (lb1==0 .Or. lb1==25 .Or. lb1==26 .Or. lb1==27 .Or. lb1==28 .Or. lb1==29 .Or. iabs(lb1)==30 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13) .Or. (iksdcy==1 .And. lb1==24) .Or. iabs(lb1)==16 .Or. (ipi0dcy==1 .And. nt==ntmax .And. lb1==4)) Then
           !lin-4/2012-above for option of pi0 decay:
           !     &           .or.iabs(lb1).eq.16) then
           Continue
        Else
           Goto 1
        End If
        ! IF I1 IS A RESONANCE, CHECK WHETHER IT DECAYS DURING THIS TIME STEP
        If (lb1>=25 .And. lb1<=27) Then
           wid = 0.151
        Else If (lb1==28) Then
           wid = 0.00841
        Else If (lb1==29) Then
           wid = 0.00443
        Else If (iabs(lb1)==30) Then
           wid = 0.051
        Else If (lb1==0) Then
           wid = 1.18E-6
           !     to give K0short ct0=2.676cm:
        Else If (iksdcy==1 .And. lb1==24) Then
           wid = 7.36E-15
           !lin-4/29/03 add Sigma0 decay to Lambda, ct0=2.22E-11m:
        Else If (iabs(lb1)==16) Then
           wid = 8.87E-6
           !sp-07/25/01 test a1 resonance:
           !c          ELSEIF(LB1.EQ.32) then
           !c             WID=0.40
        Else If (lb1==32) Then
           Call wida1(em1, rhomp, wid, iseed)
        Else If (iabs(lb1)>=6 .And. iabs(lb1)<=9) Then
           wid = width(em1)
        Else If ((iabs(lb1)==10) .Or. (iabs(lb1)==11)) Then
           wid = w1440(em1)
        Else If ((iabs(lb1)==12) .Or. (iabs(lb1)==13)) Then
           wid = w1535(em1)
           !lin-4/2012 for option of pi0 decay:
        Else If (ipi0dcy==1 .And. nt==ntmax .And. lb1==4) Then
           wid = 7.85E-9
        End If

        ! if it is the last time step, FORCE all resonance to strong-decay
        ! and go out of the loop
        If (nt==ntmax) Then
           pdecay = 1.1
           !lin-5b/2008 forbid phi decay at the end of hadronic cascade:
           If (iphidcy==0 .And. iabs(lb1)==29) pdecay = 0.
           !test off clin-9/2012 forbid long-time decays (eta,omega,K*,Sigma0)
           !     at the end of hadronic cascade to analyze freezeout time:
           !             if(LB1.eq.0.or.LB1.eq.28.or.iabs(LB1).eq.30
           !     1            .or.iabs(LB1).eq.16) pdecay=0.
        Else
           t0 = 0.19733/wid
           gfactr = e1/em1
           t0 = t0*gfactr
           If (t0>0.) Then
              pdecay = 1. - exp(-dt/t0)
           Else
              pdecay = 0.
           End If
        End If
        xdecay = ranart(nseed)

        !c dilepton production from rho0, omega, phi decay
        !c        if(lb1.eq.26 .or. lb1.eq.28 .or. lb1.eq.29)
        !c     &   call dec_ceres(nt,ntmax,irun,i1)
        !c
        If (xdecay<pdecay) Then
           !lin-10/25/02 get rid of argument usage mismatch in rhocay():
           idecay = irun
           tfnl = nt*dt
           !lin-10/28/03 keep formation time of hadrons unformed at nt=ntmax-1:
           If (nt==ntmax .And. ftsv(i1)>((ntmax-1)*dt)) tfnl = ftsv(i1)
           xfnl = x1
           yfnl = y1
           zfnl = z1
           ! use PYTHIA to perform decays of eta,rho,omega,phi,K*,(K0s) and Delta:
           If (lb1==0 .Or. lb1==25 .Or. lb1==26 .Or. lb1==27 .Or. lb1==28 .Or. lb1==29 .Or. iabs(lb1)==30 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=9) .Or. (iksdcy==1 .And. lb1==24) .Or. iabs(lb1)==16 .Or. (ipi0dcy==1 .And. nt==ntmax .And. lb1==4)) Then
              !lin-4/2012 Above for option of pi0 decay:
              !     &           .or.iabs(lb1).eq.16) then
              !     previous rho decay performed in rhodecay():
              !                nnn=nnn+1
              !                call rhodecay(idecay,i1,nnn,iseed)
              !
              !test off record decays of phi,K*,Lambda(1520) resonances:
              !                if(lb1.eq.29.or.iabs(lb1).eq.30)
              !     1               write(18,112) 'decay',lb1,px1,py1,pz1,am1,nt
              !
              !lin-4/2012 option of pi0 decays:
              !                call resdec(i1,nt,nnn,wid,idecay)
              Call resdec(i1, nt, nnn, wid, idecay, 0)
              p(1, i1) = px1n
              p(2, i1) = py1n
              p(3, i1) = pz1n
              !lin-5/2008:
              dpertp(i1) = dp1n
              !     add decay time to freezeout positions & time at the last timestep:
              If (nt==ntmax) Then
                 r(1, i1) = xfnl
                 r(2, i1) = yfnl
                 r(3, i1) = zfnl
                 tfdcy(i1) = tfnl
              End If
              !
              ! decay number for baryon resonance or L/S decay
              If (iabs(lb1)>=6 .And. iabs(lb1)<=9) Then
                 ldecay = ldecay + 1
              End If

              ! for a1 decay
              !             elseif(lb1.eq.32)then
              !                NNN=NNN+1
              !                call a1decay(idecay,i1,nnn,iseed,rhomp)

              ! FOR N*(1440)
           Else If (iabs(lb1)==10 .Or. iabs(lb1)==11) Then
              nnn = nnn + 1
              ldecay = ldecay + 1
              pnstar = 1.
              If (e(i1)>1.22) pnstar = 0.6
              If (ranart(nseed)<=pnstar) Then
                 ! (1) DECAY TO SINGLE PION+NUCLEON
                 Call decay(idecay, i1, nnn, iseed, wid, nt)
              Else
                 ! (2) DECAY TO TWO PIONS + NUCLEON
                 Call decay2(idecay, i1, nnn, iseed, wid, nt)
                 nnn = nnn + 1
              End If
              ! for N*(1535) decay
           Else If (iabs(lb1)==12 .Or. iabs(lb1)==13) Then
              nnn = nnn + 1
              Call decay(idecay, i1, nnn, iseed, wid, nt)
              ldecay = ldecay + 1
           End If
           If (nt==ntmax) Then
              If (lb(i1)==25 .Or. lb(i1)==26 .Or. lb(i1)==27) Then
                 wid = 0.151
              Else If (lb(i1)==0) Then
                 wid = 1.18E-6
              Else If (lb(i1)==24 .And. iksdcy==1) Then
                 !lin-4/2012 corrected K0s decay width:
                 !                   wid=7.36e-17
                 wid = 7.36E-15
                 !lin-4/2012 option of pi0 decays:
              Else If (ipi0dcy==1 .And. lb(i1)==4) Then
                 wid = 7.85E-9
              Else
                 Goto 9000
              End If
              lb1 = lb(i1)
              px1 = p(1, i1)
              py1 = p(2, i1)
              pz1 = p(3, i1)
              em1 = e(i1)
              e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
              !lin-4/2012 option of pi0 decays:
              !                call resdec(i1,nt,nnn,wid,idecay)
              Call resdec(i1, nt, nnn, wid, idecay, 0)
              p(1, i1) = px1n
              p(2, i1) = py1n
              p(3, i1) = pz1n
              r(1, i1) = xfnl
              r(2, i1) = yfnl
              r(3, i1) = zfnl
              tfdcy(i1) = tfnl
              !lin-5/2008:
              dpertp(i1) = dp1n
           End If

           !     Decay daughter of the above decay in lb(i1) may be a pi0:
           If (nt==ntmax .And. ipi0dcy==1 .And. lb(i1)==4) Then
              wid = 7.85E-9
              lb1 = lb(i1)
              px1 = p(1, i1)
              py1 = p(2, i1)
              pz1 = p(3, i1)
              em1 = e(i1)
              e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
              Call resdec(i1, nt, nnn, wid, idecay, 0)
              p(1, i1) = px1n
              p(2, i1) = py1n
              p(3, i1) = pz1n
              r(1, i1) = xfnl
              r(2, i1) = yfnl
              r(3, i1) = zfnl
              tfdcy(i1) = tfnl
              dpertp(i1) = dp1n
           End If

           ! negelecting the Pauli blocking at high energies
           !lin-4/2012 option of pi0 decays:
           ! 9000        go to 800
9000       Goto 798

        End If
        ! LOOP OVER ALL PSEUDOPARTICLES 2 IN THE SAME RUN
        ! SAVE ALL THE COORDINATES FOR POSSIBLE CHANGE IN THE FOLLOWING COLLISION
        !lin-4/2012 option of pi0 decays:
        ! 1        if(nt.eq.ntmax)go to 800
1       If (nt==ntmax) Goto 798

        x1 = r(1, i1)
        y1 = r(2, i1)
        z1 = r(3, i1)
        !
        Do j2 = 1, j1 - 1
           i2 = j2 + msum
           ! IF I2 IS A MESON BEING ABSORBED, THEN GO OUT OF THE LOOP
           If (e(i2)==0.) Goto 600
           !lin-5/2008 in case the first particle is already destroyed:
           If (e(i1)==0.) Goto 800
           !lin-4/2012 option of pi0 decays:
           If (lb(i2)<-45 .Or. lb(i2)>45) Goto 600
           !lin-7/26/03 improve speed
           x2 = r(1, i2)
           y2 = r(2, i2)
           z2 = r(3, i2)
           dr0max = 5.
           !lin-9/2008 deuteron+nucleon elastic cross sections could reach ~2810mb:
           ilb1 = iabs(lb(i1))
           ilb2 = iabs(lb(i2))
           If (ilb1==42 .Or. ilb2==42) Then
              If ((ilb1>=1 .And. ilb1<=2) .Or. (ilb1>=6 .And. ilb1<=13) .Or. (ilb2>=1 .And. ilb2<=2) .Or. (ilb2>=6 .And. ilb2<=13)) Then
                 If ((lb(i1)*lb(i2))>0) dr0max = 10.
              End If
           End If
           !
           If (((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)>dr0max**2) Goto 600
           If (id(i1)*id(i2)==iavoid) Goto 400
           id1 = id(i1)
           id2 = id(i2)
           !
           ix1 = nint(x1/dx)
           iy1 = nint(y1/dy)
           iz1 = nint(z1/dz)
           px1 = p(1, i1)
           py1 = p(2, i1)
           pz1 = p(3, i1)
           em1 = e(i1)
           am1 = em1
           lb1 = lb(i1)
           e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
           ipx1 = nint(px1/dpx)
           ipy1 = nint(py1/dpy)
           ipz1 = nint(pz1/dpz)
           lb2 = lb(i2)
           px2 = p(1, i2)
           py2 = p(2, i2)
           pz2 = p(3, i2)
           em2 = e(i2)
           am2 = em2
           lb1i = lb(i1)
           lb2i = lb(i2)
           px1i = p(1, i1)
           py1i = p(2, i1)
           pz1i = p(3, i1)
           em1i = e(i1)
           px2i = p(1, i2)
           py2i = p(2, i2)
           pz2i = p(3, i2)
           em2i = e(i2)
           !lin-2/26/03 ctest off check energy conservation after each binary search:
           eini = sqrt(e(i1)**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2) + sqrt(e(i2)**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
           pxini = p(1, i1) + p(1, i2)
           pyini = p(2, i1) + p(2, i2)
           pzini = p(3, i1) + p(3, i2)
           nnnini = nnn
           !
           !lin-4/30/03 initialize value:
           iblock = 0
           !
           ! TO SAVE COMPUTING TIME we do the following
           ! (1) make a ROUGH estimate to see whether particle i2 will collide with
           ! particle I1, and (2) skip the particle pairs for which collisions are
           ! not modeled in the code.
           ! FOR MESON-BARYON AND MESON-MESON COLLISIONS, we use a maximum
           ! interaction distance DELTR0=2.6
           ! for ppbar production from meson (pi rho omega) interactions:
           !
           deltr0 = 3.
           If ((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. (iabs(lb1)>=30 .And. iabs(lb1)<=45)) deltr0 = 5.0
           If ((iabs(lb2)>=14 .And. iabs(lb2)<=17) .Or. (iabs(lb2)>=30 .And. iabs(lb2)<=45)) deltr0 = 5.0

           If (lb1==28 .And. lb2==28) deltr0 = 4.84
           !lin-10/08/00 to include pi pi -> rho rho:
           If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
              e2 = sqrt(em2**2+px2**2+py2**2+pz2**2)
              spipi = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
              If (spipi>=(4*0.77**2)) deltr0 = 3.5
           End If

           ! khyperon
           If (lb1==23 .And. (lb2>=14 .And. lb2<=17)) Goto 3699
           If (lb2==23 .And. (lb1>=14 .And. lb1<=17)) Goto 3699

           ! K(K*) + Kbar(K*bar) scattering including
           !     K(K*) + Kbar(K*bar) --> phi + pi(rho,omega) and pi pi(rho,omega)
           If (lb1==21 .And. lb2==23) Goto 3699
           If (lb2==21 .And. lb1==23) Goto 3699
           If (lb1==30 .And. lb2==21) Goto 3699
           If (lb2==30 .And. lb1==21) Goto 3699
           If (lb1==-30 .And. lb2==23) Goto 3699
           If (lb2==-30 .And. lb1==23) Goto 3699
           If (lb1==-30 .And. lb2==30) Goto 3699
           If (lb2==-30 .And. lb1==30) Goto 3699
           !
           !lin-12/15/00
           !     kaon+rho(omega,eta) collisions:
           If (lb1==21 .Or. lb1==23) Then
              If (lb2==0 .Or. (lb2>=25 .And. lb2<=28)) Then
                 Goto 3699
              End If
           Else If (lb2==21 .Or. lb2==23) Then
              If (lb1==0 .Or. (lb1>=25 .And. lb1<=28)) Then
                 Goto 3699
              End If
           End If

           !lin-8/14/02 K* (pi, rho, omega, eta) collisions:
           If (iabs(lb1)==30 .And. (lb2==0 .Or. (lb2>=25 .And. lb2<=28) .Or. (lb2>=3 .And. lb2<=5))) Then
              Goto 3699
           Else If (iabs(lb2)==30 .And. (lb1==0 .Or. (lb1>=25 .And. lb1<=28) .Or. (lb1>=3 .And. lb1<=5))) Then
              Goto 3699
              !lin-8/14/02-end
              ! K*/K*-bar + baryon/antibaryon collisions:
           Else If (iabs(lb1)==30 .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Then
              Goto 3699
           End If
           If (iabs(lb2)==30 .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Then
              Goto 3699
           End If
           ! K^+ baryons and antibaryons:
           !** K+ + B-bar  --> La(Si)-bar + pi
           ! K^- and antibaryons, note K^- and baryons are included in newka():
           ! note that we fail to satisfy charge conjugation for these cross sections:
           If ((lb1==23 .Or. lb1==21) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Then
              Goto 3699
           Else If ((lb2==23 .Or. lb2==21) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Then
              Goto 3699
           End If
           !
           ! For anti-nucleons annihilations:
           ! Assumptions:
           ! (1) for collisions involving a p_bar or n_bar,
           ! we allow only collisions between a p_bar and a baryon or a baryon
           ! resonance (as well as a n_bar and a baryon or a baryon resonance),
           ! we skip all other reactions involving a p_bar or n_bar,
           ! such as collisions between p_bar (n_bar) and mesons,
           ! and collisions between two p_bar's (n_bar's).
           ! (2) we introduce a new parameter rppmax: the maximum interaction
           ! distance to make the quick collision check,rppmax=3.57 fm
           ! corresponding to a cutoff of annihilation xsection= 400mb which is
           ! also used consistently in the actual annihilation xsection to be
           ! used in the following as given in the subroutine xppbar(srt)
           rppmax = 3.57
           ! anti-baryon on baryons
           If ((lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)) .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13))) Then
              deltr0 = rppmax
              Goto 2699
           Else If ((lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6)) .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=13))) Then
              deltr0 = rppmax
              Goto 2699
           End If

           !*  ((anti) lambda, cascade, omega  should not be rejected)
           If ((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. (iabs(lb2)>=14 .And. iabs(lb2)<=17)) Goto 3699
           !
           !lin-9/2008 maximum sigma~2810mb for deuteron+nucleon elastic collisions:
           If (iabs(lb1)==42 .Or. iabs(lb2)==42) Then
              ilb1 = iabs(lb1)
              ilb2 = iabs(lb2)
              If ((ilb1>=1 .And. ilb1<=2) .Or. (ilb1>=6 .And. ilb1<=13) .Or. (ilb2>=1 .And. ilb2<=2) .Or. (ilb2>=6 .And. ilb2<=13)) Then
                 If ((lb1*lb2)>0) deltr0 = 9.5
              End If
           End If
           !
           If ((iabs(lb1)>=40 .And. iabs(lb1)<=45) .Or. (iabs(lb2)>=40 .And. iabs(lb2)<=45)) Goto 3699
           !
           !* phi channel --> elastic + inelastic scatt.
           If ((lb1==29 .And. ((lb2>=1 .And. lb2<=13) .Or. (lb2>=21 .And. lb2<=28) .Or. iabs(lb2)==30)) .Or. (lb2==29 .And. ((lb1>=1 .And. lb1<=13) .Or. (lb1>=21 .And. lb1<=28) .Or. iabs(lb1)==30))) Then
              deltr0 = 3.0
              Goto 3699
           End If
           !
           !  La/Si, Cas, Om (bar)-meson elastic colln
           ! pion vs. La & Ca (bar) coll. are treated in resp. subroutines

           ! SKIP all other K* RESCATTERINGS
           If (iabs(lb1)==30 .Or. iabs(lb2)==30) Goto 400
           ! SKIP KAON(+) RESCATTERINGS WITH particles other than pions and baryons
           If (lb1==23 .And. (lb2<1 .Or. lb2>17)) Goto 400
           If (lb2==23 .And. (lb1<1 .Or. lb1>17)) Goto 400
           !
           ! anti-baryon proccess: B-bar+M, N-bar+R-bar, N-bar+N-bar, R-bar+R-bar
           !  R = (D,N*)
           If (((lb1<=-1 .And. lb1>=-13) .And. (lb2==0 .Or. (lb2>=3 .And. lb2<=5) .Or. (lb2>=25 .And. lb2<=28))) .Or. ((lb2<=-1 .And. lb2>=-13) .And. (lb1==0 .Or. (lb1>=3 .And. lb1<=5) .Or. (lb1>=25 .And. lb1<=28)))) Then
           Else If (((lb1==-1 .Or. lb1==-2) .And. (lb2<-5 .And. lb2>=-13)) .Or. ((lb2==-1 .Or. lb2==-2) .And. (lb1<-5 .And. lb1>=-13))) Then
           Else If ((lb1==-1 .Or. lb1==-2) .And. (lb2==-1 .Or. lb2==-2)) Then
           Else If ((lb1<-5 .And. lb1>=-13) .And. (lb2<-5 .And. lb2>=-13)) Then
              !        elseif((lb1.lt.0).or.(lb2.lt.0)) then
              !         go to 400
           End If

2699       Continue
           ! for baryon-baryon collisions
           If (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=17)) Then
              If (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=17)) Then
                 deltr0 = 2.
              End If
           End If
           !
3699       rsqare = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
           If (rsqare>deltr0**2) Goto 400
           !NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER !
           ! KEEP ALL COORDINATES FOR POSSIBLE PHASE SPACE CHANGE
           ix2 = nint(x2/dx)
           iy2 = nint(y2/dy)
           iz2 = nint(z2/dz)
           ipx2 = nint(px2/dpx)
           ipy2 = nint(py2/dpy)
           ipz2 = nint(pz2/dpz)
           ! FIND MOMENTA OF PARTICLES IN THE CMS OF THE TWO COLLIDING PARTICLES
           ! AND THE CMS ENERGY SRT
           Call cms(i1, i2, pcx, pcy, pcz, srt)

           !lin-7/26/03 improve speed
           drmax = dr0max
           Call distc0(drmax, deltr0, dt, ifirst, pcx, pcy, pcz, x1, y1, z1, px1, py1, pz1, em1, x2, y2, z2, px2, py2, pz2, em2)
           If (ifirst==-1) Goto 400

           iss = nint(srt/esbin)
           !lin-4/2008 use last bin if ISS is out of EKAON's upper bound of 2000:
           If (iss>2000) iss = 2000
           !Sort collisions
           !
           !lin-8/2008 Deuteron+Meson->B+B;
           !     meson=(pi,rho,omega,eta), B=(n,p,Delta,N*1440,N*1535):
           If (iabs(lb1)==42 .Or. iabs(lb2)==42) Then
              ilb1 = iabs(lb1)
              ilb2 = iabs(lb2)
              If (lb1==0 .Or. (lb1>=3 .And. lb1<=5) .Or. (lb1>=25 .And. lb1<=28) .Or. lb2==0 .Or. (lb2>=3 .And. lb2<=5) .Or. (lb2>=25 .And. lb2<=28)) Then
                 Goto 505
                 !lin-9/2008 Deuteron+Baryon or antiDeuteron+antiBaryon elastic collisions:
              Else If (((ilb1>=1 .And. ilb1<=2) .Or. (ilb1>=6 .And. ilb1<=13) .Or. (ilb2>=1 .And. ilb2<=2) .Or. (ilb2>=6 .And. ilb2<=13)) .And. (lb1*lb2)>0) Then
                 Goto 506
              Else
                 Goto 400
              End If
           End If
           !
           ! K+ + (N,N*,D)-bar --> L/S-bar + pi
           If (((lb1==23 .Or. lb1==30) .And. (lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6))) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)))) Then
              bmass = 0.938
              If (srt<=(bmass+aka)) Then
                 pkaon = 0.
              Else
                 pkaon = sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
              End If
              !lin-10/31/02 cross sections are isospin-averaged, same as those in newka
              !     for K- + (N,N*,D) --> L/S + pi:
              sigela = 0.5*(akpel(pkaon)+aknel(pkaon))
              sigsgm = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
              sig = sigela + sigsgm + akplam(pkaon)
              If (sig>1.E-7) Then
                 !     ! K+ + N-bar reactions
                 icase = 3
                 brel = sigela/sig
                 brsgm = sigsgm/sig
                 brsig = sig
                 nchrg = 1
                 Goto 3555
              End If
              Goto 400
           End If
           !
           !
           !  meson + hyperon-bar -> K+ + N-bar
           If (((lb1>=-17 .And. lb1<=-14) .And. (lb2>=3 .And. lb2<=5)) .Or. ((lb2>=-17 .And. lb2<=-14) .And. (lb1>=3 .And. lb1<=5))) Then
              nchrg = -100

              !*       first classify the reactions due to total charge.
              If ((lb1==-15 .And. (lb2==5 .Or. lb2==27)) .Or. (lb2==-15 .And. (lb1==5 .Or. lb1==27))) Then
                 nchrg = -2
                 !     ! D-(bar)
                 bmass = 1.232
                 Goto 110
              End If
              If ((lb1==-15 .And. (lb2==0 .Or. lb2==4 .Or. lb2==26 .Or. lb2==28)) .Or. (lb2==-15 .And. (lb1==0 .Or. lb1==4 .Or. lb1==26 .Or. lb1==28)) .Or. ((lb1==-14 .Or. lb1==-16) .And. (lb2==5 .Or. lb2==27)) .Or. ((lb2==-14 .Or. lb2==-16) .And. (lb1==5 .Or. lb1==27))) Then
                 nchrg = -1
                 !     ! n-bar
                 bmass = 0.938
                 Goto 110
              End If
              If ((lb1==-15 .And. (lb2==3 .Or. lb2==25)) .Or. (lb2==-15 .And. (lb1==3 .Or. lb1==25)) .Or. (lb1==-17 .And. (lb2==5 .Or. lb2==27)) .Or. (lb2==-17 .And. (lb1==5 .Or. lb1==27)) .Or. ((lb1==-14 .Or. lb1==-16) .And. (lb2==0 .Or. lb2==4 .Or. lb2==26 .Or. lb2==28)) .Or. ((lb2==-14 .Or. lb2==-16) .And. (lb1==0 .Or. lb1==4 .Or. lb1==26 .Or. lb1==28))) Then
                 nchrg = 0
                 !     ! p-bar
                 bmass = 0.938
                 Goto 110
              End If
              If ((lb1==-17 .And. (lb2==0 .Or. lb2==4 .Or. lb2==26 .Or. lb2==28)) .Or. (lb2==-17 .And. (lb1==0 .Or. lb1==4 .Or. lb1==26 .Or. lb1==28)) .Or. ((lb1==-14 .Or. lb1==-16) .And. (lb2==3 .Or. lb2==25)) .Or. ((lb2==-14 .Or. lb2==-16) .And. (lb1==3 .Or. lb1==25))) Then
                 nchrg = 1
                 !     ! D++(bar)
                 bmass = 1.232
              End If
              !
              ! 110     if(nchrg.ne.-100.and.srt.ge.(aka+bmass))then !! for elastic
110           sig = 0.
              ! !! for elastic
              If (nchrg/=-100 .And. srt>=(aka+bmass)) Then
                 !c110        if(nchrg.eq.-100.or.srt.lt.(aka+bmass)) go to 400
                 !             ! PI + La(Si)-bar => K+ + N-bar reactions
                 icase = 4
                 !c       pkaon=sqrt(((srt**2-(aka**2+bmass**2))/2./bmass)**2-aka**2)
                 pkaon = sqrt(((srt**2-(aka**2+0.938**2))/2./0.938)**2-aka**2)
                 ! ! lambda-bar + Pi
                 If (lb1==-14 .Or. lb2==-14) Then
                    If (nchrg>=0) sigma0 = akplam(pkaon)
                    If (nchrg<0) sigma0 = aknlam(pkaon)
                    !                ! sigma-bar + pi
                 Else
                    ! !K-p or K-D++
                    If (nchrg>=0) sigma0 = akpsgm(pkaon)
                    ! !K-n or K-D-
                    If (nchrg<0) sigma0 = aknsgm(pkaon)
                    sigma0 = 1.5*akpsgm(pkaon) + aknsgm(pkaon)
                 End If
                 sig = (srt**2-(aka+bmass)**2)*(srt**2-(aka-bmass)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)*sigma0
                 ! ! K0barD++, K-D-
                 If (nchrg==-2 .Or. nchrg==2) sig = 2.*sig
                 !*     the factor 2 comes from spin of delta, which is 3/2
                 !*     detailed balance. copy from Page 423 of N.P. A614 1997
                 If (lb1==-14 .Or. lb2==-14) Then
                    sig = 4.0/3.0*sig
                 Else If (nchrg==-2 .Or. nchrg==2) Then
                    sig = 8.0/9.0*sig
                 Else
                    sig = 4.0/9.0*sig
                 End If
                 !c        brel=0.
                 !c        brsgm=0.
                 !c        brsig = sig
                 !c          if(sig.lt.1.e-7) go to 400
                 !-
              End If
              !                ! PI + La(Si)-bar => elastic included
              icase = 4
              sigela = 10.
              sig = sig + sigela
              brel = sigela/sig
              brsgm = 0.
              brsig = sig
              !-
              Goto 3555
           End If

           !* MULTISTRANGE PARTICLE (Cas,Omega -bar) PRODUCTION - (NON)PERTURBATIVE

           ! K-/K*0bar + La/Si --> cascade + pi/eta
           If (((lb1==21 .Or. lb1==-30) .And. (lb2>=14 .And. lb2<=17)) .Or. ((lb2==21 .Or. lb2==-30) .And. (lb1>=14 .And. lb1<=17))) Then
              kp = 0
              Goto 3455
           End If
           ! K+/K*0 + La/Si(bar) --> cascade-bar + pi/eta
           If (((lb1==23 .Or. lb1==30) .And. (lb2<=-14 .And. lb2>=-17)) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1<=-14 .And. lb1>=-17))) Then
              kp = 1
              Goto 3455
           End If
           ! K-/K*0bar + cascade --> omega + pi
           If (((lb1==21 .Or. lb1==-30) .And. (lb2==40 .Or. lb2==41)) .Or. ((lb2==21 .Or. lb2==-30) .And. (lb1==40 .Or. lb1==41))) Then
              kp = 0
              Goto 3455
           End If
           ! K+/K*0 + cascade-bar --> omega-bar + pi
           If (((lb1==23 .Or. lb1==30) .And. (lb2==-40 .Or. lb2==-41)) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1==-40 .Or. lb1==-41))) Then
              kp = 1
              Goto 3455
           End If
           ! Omega + Omega --> Di-Omega + photon(eta)
           !c        if( lb1.eq.45.and.lb2.eq.45 ) go to 3455

           ! annhilation of cascade(bar), omega(bar)
           kp = 3
           ! K- + L/S <-- cascade(bar) + pi/eta
           If ((((lb1>=3 .And. lb1<=5) .Or. lb1==0) .And. (iabs(lb2)==40 .Or. iabs(lb2)==41)) .Or. (((lb2>=3 .And. lb2<=5) .Or. lb2==0) .And. (iabs(lb1)==40 .Or. iabs(lb1)==41))) Goto 3455
           ! K- + cascade(bar) <-- omega(bar) + pi
           !         if(  (lb1.eq.0.and.iabs(lb2).eq.45)
           !    &       .OR. (lb2.eq.0.and.iabs(lb1).eq.45) )go to 3455
           If (((lb1>=3 .And. lb1<=5) .And. iabs(lb2)==45) .Or. ((lb2>=3 .And. lb2<=5) .And. iabs(lb1)==45)) Goto 3455
           !

           !**  MULTISTRANGE PARTICLE PRODUCTION  (END)

           !* K+ + La(Si) --> Meson + B
           If (lb1==23 .And. (lb2>=14 .And. lb2<=17)) Goto 5699
           If (lb2==23 .And. (lb1>=14 .And. lb1<=17)) Goto 5699
           !* K- + La(Si)-bar --> Meson + B-bar
           If (lb1==21 .And. (lb2>=-17 .And. lb2<=-14)) Goto 5699
           If (lb2==21 .And. (lb1>=-17 .And. lb1<=-14)) Goto 5699

           ! La/Si-bar + B --> pi + K+
           If ((((lb1==1 .Or. lb1==2) .Or. (lb1>=6 .And. lb1<=13)) .And. (lb2>=-17 .And. lb2<=-14)) .Or. (((lb2==1 .Or. lb2==2) .Or. (lb2>=6 .And. lb2<=13)) .And. (lb1>=-17 .And. lb1<=-14))) Goto 5999
           ! La/Si + B-bar --> pi + K-
           If ((((lb1==-1 .Or. lb1==-2) .Or. (lb1<=-6 .And. lb1>=-13)) .And. (lb2>=14 .And. lb2<=17)) .Or. (((lb2==-1 .Or. lb2==-2) .Or. (lb2<=-6 .And. lb2>=-13)) .And. (lb1>=14 .And. lb1<=17))) Goto 5999
           !
           !
           ! K(K*) + Kbar(K*bar) --> phi + pi(rho,omega), M + M (M=pi,rho,omega,eta)
           If (lb1==21 .And. lb2==23) Goto 8699
           If (lb2==21 .And. lb1==23) Goto 8699
           If (lb1==30 .And. lb2==21) Goto 8699
           If (lb2==30 .And. lb1==21) Goto 8699
           If (lb1==-30 .And. lb2==23) Goto 8699
           If (lb2==-30 .And. lb1==23) Goto 8699
           If (lb1==-30 .And. lb2==30) Goto 8699
           If (lb2==-30 .And. lb1==30) Goto 8699
           !* (K,K*)-bar + rho(omega) --> phi +(K,K*)-bar, piK and elastic
           If (((lb1==23 .Or. lb1==21 .Or. iabs(lb1)==30) .And. (lb2>=25 .And. lb2<=28)) .Or. ((lb2==23 .Or. lb2==21 .Or. iabs(lb2)==30) .And. (lb1>=25 .And. lb1<=28))) Goto 8799
           !
           !* K*(-bar) + pi --> phi + (K,K*)-bar
           If ((iabs(lb1)==30 .And. (lb2>=3 .And. lb2<=5)) .Or. (iabs(lb2)==30 .And. (lb1>=3 .And. lb1<=5))) Goto 8799
           !
           !
           !* phi + N --> pi+N(D),  rho+N(D),  K+ +La
           !* phi + D --> pi+N(D),  rho+N(D)
           If ((lb1==29 .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=9))) .Or. (lb2==29 .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=9)))) Goto 7222
           !
           !* phi + (pi,rho,ome,K,K*-bar) --> K+K, K+K*, K*+K*, (pi,rho,omega)+(K,K*-bar)
           If ((lb1==29 .And. ((lb2>=3 .And. lb2<=5) .Or. (lb2>=21 .And. lb2<=28) .Or. iabs(lb2)==30)) .Or. (lb2==29 .And. ((lb1>=3 .And. lb1<=5) .Or. (lb1>=21 .And. lb1<=28) .Or. iabs(lb1)==30))) Then
              Goto 7444
           End If
           !
           !
           ! La/Si, Cas, Om (bar)-(rho,omega,phi) elastic colln
           ! pion vs. La, Ca, Omega-(bar) elastic coll. treated in resp. subroutines
           If (((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. iabs(lb1)>=40) .And. ((lb2>=25 .And. lb2<=29) .Or. lb2==0)) Goto 888
           If (((iabs(lb2)>=14 .And. iabs(lb2)<=17) .Or. iabs(lb2)>=40) .And. ((lb1>=25 .And. lb1<=29) .Or. lb1==0)) Goto 888
           !
           ! K+/K* (N,R)  OR   K-/K*- (N,R)-bar  elastic scatt
           If (((lb1==23 .Or. lb1==30) .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13))) .Or. ((lb2==23 .Or. lb2==30) .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=13)))) Goto 888
           If (((lb1==21 .Or. lb1==-30) .And. (lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6))) .Or. ((lb2==21 .Or. lb2==-30) .And. (lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)))) Goto 888
           !
           ! L/S-baryon elastic collision
           If (((lb1>=14 .And. lb1<=17) .And. (lb2>=6 .And. lb2<=13)) .Or. ((lb2>=14 .And. lb2<=17) .And. (lb1>=6 .And. lb1<=13))) Goto 7799
           If (((lb1<=-14 .And. lb1>=-17) .And. (lb2<=-6 .And. lb2>=-13)) .Or. ((lb2<=-14 .And. lb2>=-17) .And. (lb1<=-6 .And. lb1>=-13))) Goto 7799
           !
           ! skip other collns with perturbative particles or hyperon-bar
           If (iabs(lb1)>=40 .Or. iabs(lb2)>=40 .Or. (lb1<=-14 .And. lb1>=-17) .Or. (lb2<=-14 .And. lb2>=-17)) Goto 400
           !
           !
           ! anti-baryon on baryon resonaces
           If ((lb1==-1 .Or. lb1==-2 .Or. (lb1>=-13 .And. lb1<=-6)) .And. (lb2==1 .Or. lb2==2 .Or. (lb2>=6 .And. lb2<=13))) Then
              Goto 2799
           Else If ((lb2==-1 .Or. lb2==-2 .Or. (lb2>=-13 .And. lb2<=-6)) .And. (lb1==1 .Or. lb1==2 .Or. (lb1>=6 .And. lb1<=13))) Then
              Goto 2799
           End If
           !
           !lin-10/25/02 get rid of argument usage mismatch in newka():
           inewka = irun
           !        call newka(icase,irun,iseed,dt,nt,
           !lin-5/01/03 set iblock value in art1f.f, necessary for resonance studies:
           !        call newka(icase,inewka,iseed,dt,nt,
           !     &                  ictrl,i1,i2,srt,pcx,pcy,pcz)
           Call newka(icase, inewka, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, iblock)

           !lin-10/25/02-end
           If (ictrl==1) Goto 400
           !
           ! SEPARATE NUCLEON+NUCLEON( BARYON RESONANCE+ BARYON RESONANCE ELASTIC
           ! COLLISION), BARYON RESONANCE+NUCLEON AND BARYON-PION
           ! COLLISIONS INTO THREE PARTS TO CHECK IF THEY ARE GOING TO SCATTER,
           ! WE only allow L/S to COLLIDE elastically with a nucleon and meson
           If ((iabs(lb1)>=14 .And. iabs(lb1)<=17) .Or. (iabs(lb2)>=14 .And. iabs(lb2)<=17)) Goto 400
           ! IF PION+PION COLLISIONS GO TO 777
           ! if pion+eta, eta+eta to create kaons go to 777
           If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Goto 777
           If (lb1==0 .And. (lb2>=3 .And. lb2<=5)) Goto 777
           If (lb2==0 .And. (lb1>=3 .And. lb1<=5)) Goto 777
           If (lb1==0 .And. lb2==0) Goto 777
           ! we assume that rho and omega behave the same way as pions in
           ! kaon production
           ! (1) rho(omega)+rho(omega)
           If ((lb1>=25 .And. lb1<=28) .And. (lb2>=25 .And. lb2<=28)) Goto 777
           ! (2) rho(omega)+pion
           If ((lb1>=25 .And. lb1<=28) .And. (lb2>=3 .And. lb2<=5)) Goto 777
           If ((lb2>=25 .And. lb2<=28) .And. (lb1>=3 .And. lb1<=5)) Goto 777
           ! (3) rho(omega)+eta
           If ((lb1>=25 .And. lb1<=28) .And. lb2==0) Goto 777
           If ((lb2>=25 .And. lb2<=28) .And. lb1==0) Goto 777
           !
           ! if kaon+pion collisions go to 889
           If ((lb1==23 .Or. lb1==21) .And. (lb2>=3 .And. lb2<=5)) Goto 889
           If ((lb2==23 .Or. lb2==21) .And. (lb1>=3 .And. lb1<=5)) Goto 889
           !
           !lin-2/06/03 skip all other (K K* Kbar K*bar) channels:
           ! SKIP all other K and K* RESCATTERINGS
           If (iabs(lb1)==30 .Or. iabs(lb2)==30) Goto 400
           If (lb1==21 .Or. lb2==21) Goto 400
           If (lb1==23 .Or. lb2==23) Goto 400
           !
           ! IF PION+baryon COLLISION GO TO 3
           If ((lb1>=3 .And. lb1<=5) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Goto 3
           If ((lb2>=3 .And. lb2<=5) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Goto 3
           !
           ! IF rho(omega)+NUCLEON (baryon resonance) COLLISION GO TO 33
           If ((lb1>=25 .And. lb1<=28) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Goto 33
           If ((lb2>=25 .And. lb2<=28) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Goto 33
           !
           ! IF ETA+NUCLEON (baryon resonance) COLLISIONS GO TO 547
           If (lb1==0 .And. (iabs(lb2)==1 .Or. iabs(lb2)==2 .Or. (iabs(lb2)>=6 .And. iabs(lb2)<=13))) Goto 547
           If (lb2==0 .And. (iabs(lb1)==1 .Or. iabs(lb1)==2 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=13))) Goto 547
           !
           ! IF NUCLEON+BARYON RESONANCE COLLISION GO TO 44
           If ((lb1==1 .Or. lb1==2) .And. (lb2>5 .And. lb2<=13)) Goto 44
           If ((lb2==1 .Or. lb2==2) .And. (lb1>5 .And. lb1<=13)) Goto 44
           If ((lb1==-1 .Or. lb1==-2) .And. (lb2<-5 .And. lb2>=-13)) Goto 44
           If ((lb2==-1 .Or. lb2==-2) .And. (lb1<-5 .And. lb1>=-13)) Goto 44
           !
           ! IF NUCLEON+NUCLEON COLLISION GO TO 4
           If ((lb1==1 .Or. lb1==2) .And. (lb2==1 .Or. lb2==2)) Goto 4
           If ((lb1==-1 .Or. lb1==-2) .And. (lb2==-1 .Or. lb2==-2)) Goto 4
           !
           ! IF BARYON RESONANCE+BARYON RESONANCE COLLISION GO TO 444
           If ((lb1>5 .And. lb1<=13) .And. (lb2>5 .And. lb2<=13)) Goto 444
           If ((lb1<-5 .And. lb1>=-13) .And. (lb2<-5 .And. lb2>=-13)) Goto 444
           !
           ! if L/S+L/S or L/s+nucleon go to 400
           ! otherwise, develop a model for their collisions
           If ((lb1<3) .And. (lb2>=14 .And. lb2<=17)) Goto 400
           If ((lb2<3) .And. (lb1>=14 .And. lb1<=17)) Goto 400
           If ((lb1>=14 .And. lb1<=17) .And. (lb2>=14 .And. lb2<=17)) Goto 400
           !
           ! otherwise, go out of the loop
           Goto 400
           !
           !
547        If (lb1*lb2==0) Then
              ! (1) FOR ETA+NUCLEON SYSTEM, we allow both elastic collision,
              !     i.e. N*(1535) formation and kaon production
              !     the total kaon production cross section is
              !     ASSUMED to be THE SAME AS PION+NUCLEON COLLISIONS
              ! (2) for eta+baryon resonance we only allow kaon production
              ece = (em1+em2+0.02)**2
              xkaon0 = 0.
              If (srt>=1.63 .And. srt<=1.7) xkaon0 = pnlka(srt)
              If (srt>1.7) xkaon0 = pnlka(srt) + pnska(srt)
              !bz3/7/99 neutralk
              xkaon0 = 2.0*xkaon0
              !bz3/7/99 neutralk end

              ! Here we negelect eta+n inelastic collisions other than the
              ! kaon production, therefore the total inelastic cross section
              ! xkaon equals to the xkaon0 (kaon production cross section)
              xkaon = xkaon0
              ! note here the xkaon is in unit of fm**2
              xeta = xn1535(i1, i2, 0)
              If ((iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13) .Or. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) xeta = 0.
              If ((xeta+xkaon)<=1.E-06) Goto 400
              dse = sqrt((xeta+xkaon)/pi)
              deltre = dse + 0.1
              px1cm = pcx
              py1cm = pcy
              pz1cm = pcz
              ! CHECK IF N*(1535) resonance CAN BE FORMED
              Call distce(i1, i2, deltre, dse, dt, ece, srt, ic, pcx, pcy, pcz)
              If (ic==-1) Goto 400
              ekaon(4, iss) = ekaon(4, iss) + 1
              If (xkaon0/(xkaon+xeta)>ranart(nseed)) Then
                 ! kaon production, USE CREN TO CALCULATE THE MOMENTUM OF L/S K+
                 Call cren(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
                 ! kaon production
                 If (iblock==7) Then
                    lpn = lpn + 1
                 Else If (iblock==-7) Then
                 End If
                 !
                 em1 = e(i1)
                 em2 = e(i2)
                 Goto 440
              End If
              ! N*(1535) FORMATION
              resona = 1.
              Goto 98
           End If
           !IF PION+NUCLEON (baryon resonance) COLLISION THEN
3          Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ! the total kaon production cross section for pion+baryon (resonance) is
           ! assumed to be the same as in pion+nucleon
           xkaon0 = 0.
           If (srt>=1.63 .And. srt<=1.7) xkaon0 = pnlka(srt)
           If (srt>1.7) xkaon0 = pnlka(srt) + pnska(srt)
           xkaon0 = 2.0*xkaon0
           !
           ! sp11/21/01  phi production: pi +N(D) -> phi + N(D)
           xphi = 0.
           If ((((lb1>=1 .And. lb1<=2) .Or. (lb1>=6 .And. lb1<=9)) .Or. ((lb2>=1 .And. lb2<=2) .Or. (lb2>=6 .And. lb2<=9))) .And. srt>1.958) Call pibphi(srt, lb1, lb2, em1, em2, xphi, xphin)
           ! !! in fm^2 above

           ! if a pion collide with a baryon resonance,
           ! we only allow kaon production AND the reabsorption
           ! processes: Delta+pion-->N+pion, N*+pion-->N+pion
           ! Later put in pion+baryon resonance elastic
           ! cross through forming higher resonances implicitly.
           !          If(em1.gt.1.or.em2.gt.1.)go to 31
           If ((iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13) .Or. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) Goto 31
           ! For pion+nucleon collisions:
           ! using the experimental pion+nucleon inelastic cross section, we assume it
           ! is exhausted by the Delta+pion, Delta+rho and Delta+omega production
           ! and kaon production. In the following we first check whether
           ! inelastic pion+n collision can happen or not, then determine in
           ! crpn whether it is through pion production or through kaon production
           ! note that the xkaon0 is the kaon production cross section
           ! Note in particular that:
           ! xkaon in the following is the total pion+nucleon inelastic cross section
           ! note here the xkaon is in unit of fm**2, xnpi is also in unit of fm**2
           ! FOR PION+NUCLEON SYSTEM, THE MINIMUM S IS 1.2056 the minimum srt for
           ! elastic scattering, and it is 1.60 for pion production, 1.63 for LAMBDA+kaon
           ! production and 1.7 FOR SIGMA+KAON
           ! (EC = PION MASS+NUCLEON MASS+20MEV)**2
           ec = (em1+em2+0.02)**2
           xkaon = 0.
           If (srt>1.23) xkaon = (pionpp(srt)+pipp1(srt))/2.
           ! pion+nucleon elastic cross section is divided into two parts:
           ! (1) forming D(1232)+N*(1440) +N*(1535)
           ! (2) cross sections forming higher resonances are calculated as
           !     the difference between the total elastic and (1), this part is
           !     treated as direct process since we do not explicitLY include
           !     higher resonances.
           ! the following is the resonance formation cross sections.
           !1. PION(+)+PROTON-->DELTA++,PION(-)+NEUTRON-->DELTA(-)
           If ((lb1*lb2==5 .Or. ((lb1*lb2==6) .And. (lb1==3 .Or. lb2==3))) .Or. (lb1*lb2==-3 .Or. ((lb1*lb2==-10) .And. (lb1==5 .Or. lb2==5)))) Then
              xmax = 190.
              xmaxn = 0
              xmaxn1 = 0
              xdirct = dirct1(srt)
              Goto 678
           End If
           !2. PION(-)+PROTON-->DELTA0,PION(+)+NEUTRON-->DELTA+
           !   or N*(+)(1440) or N*(+)(1535)
           ! note the factor 2/3 is from the isospin consideration and
           ! the factor 0.6 or 0.5 is the branching ratio for the resonance to decay
           ! into pion+nucleon
           If ((lb1*lb2==3 .Or. ((lb1*lb2==10) .And. (lb1==5 .Or. lb2==5))) .Or. (lb1*lb2==-5 .Or. ((lb1*lb2==-6) .And. (lb1==3 .Or. lb2==3)))) Then
              xmax = 27.
              xmaxn = 2./3.*25.*0.6
              xmaxn1 = 2./3.*40.*0.5
              xdirct = dirct2(srt)
              Goto 678
           End If
           !3. PION0+PROTON-->DELTA+,PION0+NEUTRON-->DELTA0, or N*(0)(1440) or N*(0)(1535)
           If ((lb1==4 .Or. lb2==4) .And. (iabs(lb1*lb2)==4 .Or. iabs(lb1*lb2)==8)) Then
              xmax = 50.
              xmaxn = 1./3.*25*0.6
              xmaxn1 = 1/3.*40.*0.5
              xdirct = dirct3(srt)
              Goto 678
           End If
678        xnpin1 = 0
           xnpin = 0
           xnpid = xnpi(i1, i2, 1, xmax)
           If (xmaxn1/=0) xnpin1 = xnpi(i1, i2, 2, xmaxn1)
           If (xmaxn/=0) xnpin = xnpi(i1, i2, 0, xmaxn)
           ! the following
           xres = xnpid + xnpin + xnpin1
           xnelas = xres + xdirct
           icheck = 1
           Goto 34
           ! For pion + baryon resonance the reabsorption
           ! cross section is calculated from the detailed balance
           ! using reab(i1,i2,srt,ictrl), ictrl=1, 2 and 3
           ! for pion, rho and omega + baryon resonance
31         ec = (em1+em2+0.02)**2
           xreab = reab(i1, i2, srt, 1)

           !lin-12/02/00 to satisfy detailed balance, forbid N* absorptions:
           If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) xreab = 0.

           xkaon = xkaon0 + xreab
           ! a constant of 10 mb IS USED FOR PION + N* RESONANCE,
           If ((iabs(lb1)>9 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>9 .And. iabs(lb2)<=13)) Then
              xnelas = 1.0
           Else
              xnelas = dpion(em1, em2, lb1, lb2, srt)
           End If
           icheck = 2
34         If ((xnelas+xkaon+xphi)<=0.000001) Goto 400
           ds = sqrt((xnelas+xkaon+xphi)/pi)
           !sp09/20/01
           !           totcr = xnelas+xkaon
           !           if(srt .gt. 3.5)totcr = max1(totcr,3.)
           !           DS=SQRT(totcr/PI)
           !sp09/20/01 end

           deltar = ds + 0.1
           Call distce(i1, i2, deltar, ds, dt, ec, srt, ic, pcx, pcy, pcz)
           If (ic==-1) Goto 400
           ekaon(4, iss) = ekaon(4, iss) + 1
           !***
           ! check what kind of collision has happened
           ! (1) pion+baryon resonance
           ! if direct elastic process
           If (icheck==2) Then
              !  !!sp11/21/01
              If (xnelas/(xnelas+xkaon+xphi)>=ranart(nseed)) Then
                 !               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2)
                 Call crdir(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
                 Goto 440
              Else
                 ! for inelastic process, go to 96 to check
                 ! kaon production and pion reabsorption : pion+D(N*)-->pion+N
                 Goto 96
              End If
           End If
           !(2) pion+n
           ! CHECK IF inELASTIC COLLISION IS POSSIBLE FOR PION+N COLLISIONS
           !lin-8/17/00 typo corrected, many other occurences:
           !        IF(XKAON/(XKAON+Xnelas).GT.RANART(NSEED))GO TO 95
           If ((xkaon+xphi)/(xkaon+xphi+xnelas)>ranart(nseed)) Goto 95

           ! direct process
           If (xdirct/xnelas>=ranart(nseed)) Then
              !               call Crdir(PX1CM,PY1CM,PZ1CM,SRT,I1,I2)
              Call crdir(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
              Goto 440
           End If
           ! now resonance formation or direct process (higher resonances)
           If ((lb1*lb2==5 .Or. ((lb1*lb2==6) .And. (lb1==3 .Or. lb2==3))) .Or. (lb1*lb2==-3 .Or. ((lb1*lb2==-10) .And. (lb1==5 .Or. lb2==5)))) Then
              !
              ! ONLY DELTA RESONANCE IS POSSIBLE, go to 99
              Goto 99
           Else
              ! NOW BOTH DELTA AND N* RESORANCE ARE POSSIBLE
              ! DETERMINE THE RESORANT STATE BY USING THE MONTRE CARLO METHOD
              xx = (xnpin+xnpin1)/xres
              If (ranart(nseed)<xx) Then
                 ! N* RESONANCE IS SELECTED
                 ! decide N*(1440) or N*(1535) formation
                 xx0 = xnpin/(xnpin+xnpin1)
                 If (ranart(nseed)<xx0) Then
                    resona = 0.
                    ! N*(1440) formation
                    Goto 97
                 Else
                    ! N*(1535) formation
                    resona = 1.
                    Goto 98
                 End If
              Else
                 ! DELTA RESONANCE IS SELECTED
                 Goto 99
              End If
           End If
97         Continue
           If (resona==0.) Then
              !N*(1440) IS PRODUCED,WE DETERMINE THE CHARGE STATE OF THE PRODUCED N*
              i = i1
              If (em1<0.6) i = i2
              ! (0.1) n+pion(+)-->N*(+)
              If ((lb1*lb2==10 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==-6 .And. (lb1==3 .Or. lb2==3))) Then
                 lb(i) = 11
                 Goto 303
              End If
              ! (0.2) p+pion(0)-->N*(+)
              !            IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))THEN
              If (iabs(lb(i1)*lb(i2))==4 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 11
                 Goto 303
              End If
              ! (0.3) n+pion(0)-->N*(0)
              !            IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
              If (iabs(lb(i1)*lb(i2))==8 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 10
                 Goto 303
              End If
              ! (0.4) p+pion(-)-->N*(0)
              !            IF(LB(I1)*LB(I2).EQ.3)THEN
              If ((lb(i1)*lb(i2)==3) .Or. (lb(i1)*lb(i2)==-5)) Then
                 lb(i) = 10
              End If
303           Call dreson(i1, i2)
              If (lb1<0 .Or. lb2<0) lb(i) = -lb(i)
              lres = lres + 1
              Goto 101
              !COM: GO TO 101 TO CHANGE THE PHASE SPACE DENSITY OF THE NUCLEON
           End If
98         If (resona==1.) Then
              !N*(1535) IS PRODUCED, WE DETERMINE THE CHARGE STATE OF THE PRODUCED N*
              i = i1
              If (em1<0.6) i = i2
              ! note: this condition applies to both eta and pion
              ! (0.1) n+pion(+)-->N*(+)
              !            IF(LB1*LB2.EQ.10.AND.(LB1.EQ.2.OR.LB2.EQ.2))THEN
              If ((lb1*lb2==10 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==-6 .And. (lb1==3 .Or. lb2==3))) Then
                 lb(i) = 13
                 Goto 304
              End If
              ! (0.2) p+pion(0)-->N*(+)
              !            IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))THEN
              If (iabs(lb(i1)*lb(i2))==4 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 13
                 Goto 304
              End If
              ! (0.3) n+pion(0)-->N*(0)
              !            IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
              If (iabs(lb(i1)*lb(i2))==8 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
                 lb(i) = 12
                 Goto 304
              End If
              ! (0.4) p+pion(-)-->N*(0)
              !            IF(LB(I1)*LB(I2).EQ.3)THEN
              If ((lb(i1)*lb(i2)==3) .Or. (lb(i1)*lb(i2)==-5)) Then
                 lb(i) = 12
                 Goto 304
              End If
              ! (0.5) p+eta-->N*(+)(1535),n+eta-->N*(0)(1535)
              If (lb(i1)*lb(i2)==0) Then
                 !            if((lb(i1).eq.1).or.(lb(i2).eq.1))then
                 If (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1) Then
                    lb(i) = 13
                    Goto 304
                 Else
                    lb(i) = 12
                 End If
              End If
304           Call dreson(i1, i2)
              If (lb1<0 .Or. lb2<0) lb(i) = -lb(i)
              lres = lres + 1
              Goto 101
              !COM: GO TO 101 TO CHANGE THE PHASE SPACE DENSITY OF THE NUCLEON
           End If
           !DELTA IS PRODUCED,IN THE FOLLOWING WE DETERMINE THE
           !CHARGE STATE OF THE PRODUCED DELTA
99         lres = lres + 1
           i = i1
           If (em1<=0.6) i = i2
           ! (1) p+pion(+)-->DELTA(++)
           !        IF(LB(I1)*LB(I2).EQ.5)THEN
           If ((lb(i1)*lb(i2)==5) .Or. (lb(i1)*lb(i2)==-3)) Then
              lb(i) = 9
              Goto 305
           End If
           ! (2) p+pion(0)-->delta(+)
           !        IF(LB(I1)*LB(I2).EQ.4.AND.(LB(I1).EQ.1.OR.LB(I2).EQ.1))then
           If (iabs(lb(i1)*lb(i2))==4 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
              lb(i) = 8
              Goto 305
           End If
           ! (3) n+pion(+)-->delta(+)
           !        IF(LB(I1)*LB(I2).EQ.10.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
           If ((lb(i1)*lb(i2)==10 .And. (lb(i1)==5 .Or. lb(i2)==5)) .Or. (lb(i1)*lb(i2)==-6 .And. (lb(i1)==3 .Or. lb(i2)==3))) Then
              lb(i) = 8
              Goto 305
           End If
           ! (4) n+pion(0)-->delta(0)
           !        IF(LB(I1)*LB(I2).EQ.8.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
           If (iabs(lb(i1)*lb(i2))==8 .And. (lb(i1)==4 .Or. lb(i2)==4)) Then
              lb(i) = 7
              Goto 305
           End If
           ! (5) p+pion(-)-->delta(0)
           !        IF(LB(I1)*LB(I2).EQ.3)THEN
           If ((lb(i1)*lb(i2)==3) .Or. (lb(i1)*lb(i2)==-5)) Then
              lb(i) = 7
              Goto 305
           End If
           ! (6) n+pion(-)-->delta(-)
           !        IF(LB(I1)*LB(I2).EQ.6.AND.(LB(I1).EQ.2.OR.LB(I2).EQ.2))THEN
           If ((lb(i1)*lb(i2)==6 .And. (lb(i1)==3 .Or. lb(i2)==3)) .Or. (lb(i1)*lb(i2)==-10 .And. (lb(i1)==5 .Or. lb(i2)==5))) Then
              lb(i) = 6
           End If
305        Call dreson(i1, i2)
           If (lb1<0 .Or. lb2<0) lb(i) = -lb(i)
           Goto 101

           !sp-11/08/01 K*
           ! FOR kaON+pion COLLISIONS, form K* (bar) or
           ! La/Si-bar + N <-- pi + K+
           ! La/Si + N-bar <-- pi + K-
           ! phi + K <-- pi + K
           !lin (rho,omega) + K* <-- pi + K
889        Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           ! the cross section is from C.M. Ko, PRC 23, 2760 (1981).
           spika = 60./(1.+4.*(srt-0.895)**2/(0.05)**2)
           !
           !c       if(lb(i1).eq.23.or.lb(i2).eq.23)then   !! block  K- + pi->La + B-bar

           Call crkpla(px1cm, py1cm, pz1cm, ec, srt, spika, emm1, emm2, lbp1, lbp2, i1, i2, icase, srhoks)
           !c
           !* only K* or K*bar formation
           !       else
           !      DSkn=SQRT(spika/PI/10.)
           !      dsknr=dskn+0.1
           !      CALL DISTCE(I1,I2,dsknr,DSkn,DT,EC,SRT,IC,
           !    1     PX1CM,PY1CM,PZ1CM)
           !        IF(IC.EQ.-1) GO TO 400
           !       icase = 1
           !      endif
           !
           If (icase==0) Then
              iblock = 0
              Goto 400
           End If

           If (icase==1) Then
              Call ksreso(i1, i2)
              !lin-4/30/03 give non-zero iblock for resonance selections:
              iblock = 171
              !test off for resonance (phi, K*) studies:
              !             if(iabs(lb(i1)).eq.30) then
              !             write(17,112) 'ks',lb(i1),p(1,i1),p(2,i1),p(3,i1),e(i1),nt
              !             elseif(iabs(lb(i2)).eq.30) then
              !             write(17,112) 'ks',lb(i2),p(1,i2),p(2,i2),p(3,i2),e(i2),nt
              !             endif
              !
              lres = lres + 1
              Goto 101
           Else If (icase==2) Then
              iblock = 71
              !
              ! La/Si (bar) formation

           Else If (iabs(icase)==5) Then
              iblock = 88

           Else
              !
              ! phi formation
              iblock = 222
           End If
           lb(i1) = lbp1
           lb(i2) = lbp2
           e(i1) = emm1
           e(i2) = emm2
           em1 = e(i1)
           em2 = e(i2)
           ntag = 0
           Goto 440
           !
33         Continue
           em1 = e(i1)
           em2 = e(i2)
           ! (1) if rho or omega collide with a nucleon we allow both elastic
           !     scattering and kaon production to happen if collision conditions
           !     are satisfied.
           ! (2) if rho or omega collide with a baryon resonance we allow
           !     kaon production, pion reabsorption: rho(omega)+D(N*)-->pion+N
           !     and NO elastic scattering to happen
           xelstc = 0
           If ((lb1>=25 .And. lb1<=28) .And. (iabs(lb2)==1 .Or. iabs(lb2)==2)) xelstc = erhon(em1, em2, lb1, lb2, srt)
           If ((lb2>=25 .And. lb2<=28) .And. (iabs(lb1)==1 .Or. iabs(lb1)==2)) xelstc = erhon(em1, em2, lb1, lb2, srt)
           ec = (em1+em2+0.02)**2
           ! the kaon production cross section is
           xkaon0 = 0
           If (srt>=1.63 .And. srt<=1.7) xkaon0 = pnlka(srt)
           If (srt>1.7) xkaon0 = pnlka(srt) + pnska(srt)
           If (xkaon0<0) xkaon0 = 0

           !bz3/7/99 neutralk
           xkaon0 = 2.0*xkaon0
           !bz3/7/99 neutralk end

           ! the total inelastic cross section for rho(omega)+N is
           xkaon = xkaon0
           ichann = 0
           ! the total inelastic cross section for rho (omega)+D(N*) is
           ! xkaon=xkaon0+reab(**)

           ! sp11/21/01  phi production: rho + N(D) -> phi + N(D)
           xphi = 0.
           If (((((lb1>=1 .And. lb1<=2) .Or. (lb1>=6 .And. lb1<=9)) .And. (lb2>=25 .And. lb2<=27)) .Or. (((lb2>=1 .And. lb2<=2) .Or. (lb2>=6 .And. lb2<=9)) .And. (lb1>=25 .And. lb1<=27))) .And. srt>1.958) Call pibphi(srt, lb1, lb2, em1, em2, xphi, xphin)
           ! !! in fm^2 above
           !
           If ((iabs(lb1)>=6 .And. lb2>=25) .Or. (lb1>=25 .And. iabs(lb2)>=6)) Then
              ichann = 1
              ictrl = 2
              If (lb1==28 .Or. lb2==28) ictrl = 3
              xreab = reab(i1, i2, srt, ictrl)

              !lin-12/02/00 to satisfy detailed balance, forbid N* absorptions:
              If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) xreab = 0.

              If (xreab<0) xreab = 1.E-06
              xkaon = xkaon0 + xreab
              xelstc = 1.0
           End If
           ds = sqrt((xkaon+xphi+xelstc)/pi)
           !
           !sp09/20/01
           !           totcr = xelstc+xkaon
           !           if(srt .gt. 3.5)totcr = max1(totcr,3.)
           !           DS=SQRT(totcr/PI)
           !sp09/20/01 end
           !
           deltar = ds + 0.1
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ! CHECK IF the collision can happen
           Call distce(i1, i2, deltar, ds, dt, ec, srt, ic, pcx, pcy, pcz)
           If (ic==-1) Goto 400
           ekaon(4, iss) = ekaon(4, iss) + 1
           !*
           ! NOW rho(omega)+N or D(N*) COLLISION IS POSSIBLE
           ! (1) check elastic collision
           If (xelstc/(xelstc+xkaon+xphi)>ranart(nseed)) Then
              !       call crdir(px1CM,py1CM,pz1CM,srt,I1,i2)
              Call crdir(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
              Goto 440
           End If
           ! (2) check pion absorption or kaon production
           Call crrd(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)

           ! kaon production
           !sp05/16/01
           If (iblock==7) Then
              lpn = lpn + 1
           Else If (iblock==-7) Then
           End If
           !sp05/16/01 end
           ! rho obsorption
           If (iblock==81) lrhor = lrhor + 1
           ! omega obsorption
           If (iblock==82) lomgar = lomgar + 1
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           ! for pion+n now using the subroutine crpn to change
           ! the particle label and set the new momentum of L/S+K final state
95         Continue
           ! NOW PION+N INELASTIC COLLISION IS POSSIBLE
           ! check pion production or kaon production
           Call crpn(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)

           ! kaon production
           !sp05/16/01
           If (iblock==7) Then
              lpn = lpn + 1
           Else If (iblock==-7) Then
           End If
           !sp05/16/01 end
           ! pion production
           If (iblock==77) lpd = lpd + 1
           ! rho production
           If (iblock==78) lrho = lrho + 1
           ! omega production
           If (iblock==79) lomega = lomega + 1
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           ! for pion+D(N*) now using the subroutine crpd to
           ! (1) check kaon production or pion reabsorption
           ! (2) change the particle label and set the new
           !     momentum of L/S+K final state
96         Continue
           Call crpd(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)

           ! kaon production
           !sp05/16/01
           If (iblock==7) Then
              lpn = lpn + 1
           Else If (iblock==-7) Then
           End If
           !sp05/16/01 end
           ! pion obserption
           If (iblock==80) lpdr = lpdr + 1
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           ! CALCULATE KAON PRODUCTION PROBABILITY FROM PION + N COLLISIONS
           !        IF(SRT.GT.1.615)THEN
           !        CALL PKAON(SRT,XXp,PK)
           !        TKAON(7)=TKAON(7)+PK
           !        EKAON(7,ISS)=EKAON(7,ISS)+1
           !        CALL KSPEC1(SRT,PK)
           !        call LK(3,srt,iseed,pk)
           !        ENDIF
           ! negelecting the pauli blocking at high energies

101        Continue
           If (e(i2)==0.) Goto 600
           If (e(i1)==0.) Goto 800
           ! IF NUCLEON+BARYON RESONANCE COLLISIONS
44         Continue
           ! CALCULATE THE TOTAL CROSS SECTION OF NUCLEON+ BARYON RESONANCE COLLISION
           ! WE ASSUME THAT THE ELASTIC CROSS SECTION IS THE SAME AS NUCLEON+NUCLEON
           ! COM: WE USE THE PARAMETERISATION BY CUGNON FOR LOW ENERGIES
           !      AND THE PARAMETERIZATIONS FROM CERN DATA BOOK FOR HIGHER
           !      ENERGIES. THE CUTOFF FOR THE TOTAL CROSS SECTION IS 55 MB
           cutoff = em1 + em2 + 0.02
           If (srt<=cutoff) Goto 400
           If (srt>2.245) Then
              signn = pp2(srt)
           Else
              signn = 35.0/(1.+(srt-cutoff)*100.0) + 20.0
           End If
           Call xnd(pcx, pcy, pcz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
           sig = signn + xinel
           ! For nucleon+baryon resonance collision, the minimum cms**2 energy is
           ec = (em1+em2+0.02)**2
           ! CHECK THE DISTENCE BETWEEN THE TWO PARTICLES
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz

           !lin-6/2008 Deuteron production:
           ianti = 0
           If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
           Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
           sig = sig + sdprod
           !lin-6/2008 perturbative treatment of deuterons:
           ipdflag = 0
           If (idpert==1) Then
              ipert1 = 1
              sigr0 = sig
              dspert = sqrt(sigr0/pi/10.)
              dsrpert = dspert + 0.1
              Call distce(i1, i2, dsrpert, dspert, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
              If (ic==-1) Goto 363
              signn0 = 0.
              Call crnd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, signn0, sigr0, sigk, xsk1, xsk2, xsk3, xsk4, xsk5, nt, ipert1)
              !     &  IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
              ipdflag = 1
363           Continue
              ipert1 = 0
           End If
           If (idpert==2) ipert1 = 1
           !
           ds = sqrt(sig/(10.*pi))
           deltar = ds + 0.1
           Call distce(i1, i2, deltar, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           !        IF(IC.EQ.-1)GO TO 400
           If (ic==-1) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If

           ekaon(3, iss) = ekaon(3, iss) + 1
           ! CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON + BARYON RESONANCE
           ! COLLISIONS
           Goto 361

           ! CHECK WHAT KIND OF COLLISION HAS HAPPENED
361        Continue
           Call crnd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, signn, sig, sigk, xsk1, xsk2, xsk3, xsk4, xsk5, nt, ipert1)
           !     &  IBLOCK,SIGNN,SIG,sigk,xsk1,xsk2,xsk3,xsk4,xsk5)
           If (iblock==0 .And. ipdflag==1) iblock = 501
           If (iblock==11) Then
              lndk = lndk + 1
              Goto 400
              !        elseIF(IBLOCK.EQ.-11) then
           Else If (iblock==-11 .Or. iblock==501) Then
              Goto 400
           End If
           If (iblock==222) Then
              !    !! sp12/17/01
              Goto 400
           End If
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           ! IF NUCLEON+NUCLEON OR BARYON RESONANCE+BARYON RESONANCE COLLISIONS
4          Continue
           ! PREPARE THE EALSTIC CROSS SECTION FOR BARYON+BARYON COLLISIONS
           ! COM: WE USE THE PARAMETERISATION BY CUGNON FOR SRT LEQ 2.0 GEV
           !      AND THE PARAMETERIZATIONS FROM CERN DATA BOOK FOR HIGHER
           !      ENERGIES. THE CUTOFF FOR THE TOTAL CROSS SECTION IS 55 MB
           !      WITH LOW-ENERGY-CUTOFF
           cutoff = em1 + em2 + 0.14
           ! AT HIGH ENERGIES THE ISOSPIN DEPENDENCE IS NEGLIGIBLE
           ! THE TOTAL CROSS SECTION IS TAKEN AS THAT OF THE PP
           ! ABOVE E_KIN=800 MEV, WE USE THE ISOSPIN INDEPENDNET XSECTION
           If (srt>2.245) Then
              sig = ppt(srt)
              signn = sig - pp1(srt)
           Else
              ! AT LOW ENERGIES THE ISOSPIN DEPENDENCE FOR NN COLLISION IS STRONG
              sig = xpp(srt)
              If (zet(lb(i1))*zet(lb(i2))<=0) sig = xnp(srt)
              If (zet(lb(i1))*zet(lb(i2))>0) sig = xpp(srt)
              If (zet(lb(i1))==0 .And. zet(lb(i2))==0) sig = xpp(srt)
              If ((lb(i1)==-1 .And. lb(i2)==-2) .Or. (lb(i2)==-1 .And. lb(i1)==-2)) sig = xnp(srt)
              !     WITH LOW-ENERGY-CUTOFF
              If (srt<1.897) Then
                 signn = sig
              Else
                 signn = 35.0/(1.+(srt-1.897)*100.0) + 20.0
              End If
           End If
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           !lin-5/2008 Deuteron production cross sections were not included
           !     in the previous parameterized inelastic cross section of NN collisions
           !     (SIGinel=SIG-SIGNN), so they are added here:
           ianti = 0
           If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
           Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
           sig = sig + sdprod
           !
           !lin-5/2008 perturbative treatment of deuterons:
           ipdflag = 0
           If (idpert==1) Then
              !     For idpert=1: ipert1=1 means we will first treat deuteron perturbatively,
              !     then we set ipert1=0 to treat regular NN or NbarNbar collisions including
              !     the regular deuteron productions.
              !     ipdflag=1 means perturbative deuterons are produced here:
              ipert1 = 1
              ec = 2.012**2
              !     Use the same cross section for NN/NNBAR collisions
              !     to trigger perturbative production
              sigr0 = sig
              !     One can also trigger with X*sbbdm() so the weight will not be too small;
              !     but make sure to limit the maximum trigger Xsec:
              !           sigr0=sdprod*25.
              !           if(sigr0.ge.100.) sigr0=100.
              dspert = sqrt(sigr0/pi/10.)
              dsrpert = dspert + 0.1
              Call distce(i1, i2, dsrpert, dspert, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
              If (ic==-1) Goto 365
              signn0 = 0.
              Call crnn(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn0, sigr0, nt, ipert1)
              ipdflag = 1
365           Continue
              ipert1 = 0
           End If
           If (idpert==2) ipert1 = 1
           !
           !lin-5/2008 in case perturbative deuterons are produced for idpert=1:
           !        IF(SIGNN.LE.0)GO TO 400
           If (signn<=0) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If
           !
           ec = 3.59709
           ds = sqrt(sig/pi/10.)
           dsr = ds + 0.1
           If ((e(i1)>=1.) .And. (e(i2)>=1.)) ec = 4.75
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           !lin-5/2008 in case perturbative deuterons are produced above:
           !        IF(IC.EQ.-1) GO TO 400
           If (ic==-1) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If
           !
           ! CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON+NUCLEON OR
           ! RESONANCE+RESONANCE COLLISIONS
           Goto 362

           ! CHECK WHAT KIND OF COLLISION HAS HAPPENED
362        ekaon(1, iss) = ekaon(1, iss) + 1
           Call crnn(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
           !lin-5/2008 give iblock # in case pert deuterons are produced for idpert=1:
           If (iblock==0 .And. ipdflag==1) iblock = 501
           !lin-5/2008 add iblock # for deuteron formation:
           !        IF(IBLOCK.EQ.4.OR.IBLOCK.Eq.9.or.iblock.ge.44.OR.IBLOCK.EQ.-9
           !     &       .or.iblock.eq.222)THEN
           If (iblock==4 .Or. iblock==9 .Or. iblock>=44 .Or. iblock==-9 .Or. iblock==222 .Or. iblock==501) Then
              !
              !     !! sp12/17/01 above
              ! momentum of the three particles in the final state have been calculated
              ! in the crnn, go out of the loop
              lcoll = lcoll + 1
              If (iblock==4) Then
                 ldirt = ldirt + 1
              Else If (iblock==44) Then
                 lddrho = lddrho + 1
              Else If (iblock==45) Then
                 lnnrho = lnnrho + 1
              Else If (iblock==46) Then
                 lnnom = lnnom + 1
              Else If (iblock==222) Then
              Else If (iblock==9) Then
                 lnnk = lnnk + 1
              Else If (iblock==-9) Then
              End If
              Goto 400
           End If

           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !lin-8/2008 B+B->Deuteron+Meson over
           !
           !lin-8/2008 Deuteron+Meson->B+B collisions:
505        Continue
           ianti = 0
           If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
           Call sdmbb(srt, sdm, ianti)
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           !     minimum srt**2, note a 2.012GeV lower cutoff is used in N+N->Deuteron+pi:
           ec = 2.012**2
           ds = sqrt(sdm/31.4)
           dsr = ds + 0.1
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crdmbb(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, sdm, nt, ianti)
           lcoll = lcoll + 1
           Goto 400
           !lin-8/2008 Deuteron+Meson->B+B collisions over
           !
           !lin-9/2008 Deuteron+Baryon elastic collisions:
506        Continue
           ianti = 0
           If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
           Call sdbelastic(srt, sdb)
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           !     minimum srt**2, note a 2.012GeV lower cutoff is used in N+N->Deuteron+pi:
           ec = 2.012**2
           ds = sqrt(sdb/31.4)
           dsr = ds + 0.1
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crdbel(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, sdb, nt, ianti)
           lcoll = lcoll + 1
           Goto 400
           !lin-9/2008 Deuteron+Baryon elastic collisions over
           !
           ! IF BARYON RESONANCE+BARYON RESONANCE COLLISIONS
444        Continue
           ! PREPARE THE EALSTIC CROSS SECTION FOR BARYON+BARYON COLLISIONS
           cutoff = em1 + em2 + 0.02
           ! AT HIGH ENERGIES THE ISOSPIN DEPENDENCE IS NEGLIGIBLE
           ! THE TOTAL CROSS SECTION IS TAKEN AS THAT OF THE PP
           If (srt<=cutoff) Goto 400
           If (srt>2.245) Then
              signn = pp2(srt)
           Else
              signn = 35.0/(1.+(srt-cutoff)*100.0) + 20.0
           End If
           If (signn<=0) Goto 400
           Call xddin(pcx, pcy, pcz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
           sig = signn + xinel
           ec = (em1+em2+0.02)**2
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz

           !lin-6/2008 Deuteron production:
           ianti = 0
           If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
           Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
           sig = sig + sdprod
           !lin-6/2008 perturbative treatment of deuterons:
           ipdflag = 0
           If (idpert==1) Then
              ipert1 = 1
              sigr0 = sig
              dspert = sqrt(sigr0/pi/10.)
              dsrpert = dspert + 0.1
              Call distce(i1, i2, dsrpert, dspert, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
              If (ic==-1) Goto 367
              signn0 = 0.
              Call crdd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn0, sigr0, nt, ipert1)
              !     1          IBLOCK,NTAG,SIGNN,SIG)
              ipdflag = 1
367           Continue
              ipert1 = 0
           End If
           If (idpert==2) ipert1 = 1
           !
           ds = sqrt(sig/31.4)
           dsr = ds + 0.1
           Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           !        IF(IC.EQ.-1) GO TO 400
           If (ic==-1) Then
              If (ipdflag==1) iblock = 501
              Goto 400
           End If

           ! CALCULATE KAON PRODUCTION PROBABILITY FROM NUCLEON+NUCLEON OR
           ! RESONANCE+RESONANCE COLLISIONS
           Goto 364

           ! CHECK WHAT KIND OF COLLISION HAS HAPPENED
364        ekaon(2, iss) = ekaon(2, iss) + 1
           ! for resonance+resonance
           !lin-6/2008:
           Call crdd(irun, px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
           !     1  IBLOCK,NTAG,SIGNN,SIG)
           If (iblock==0 .And. ipdflag==1) iblock = 501
           !
           If (iabs(iblock)==10) Then
              ! momentum of the three particles in the final state have been calculated
              ! in the crnn, go out of the loop
              lcoll = lcoll + 1
              If (iblock==10) Then
                 lddk = lddk + 1
              Else If (iblock==-10) Then
              End If
              Goto 400
           End If
           !lin-6/2008
           !        if(iblock .eq. 222)then
           If (iblock==222 .Or. iblock==501) Then
              !    !! sp12/17/01
              Goto 400
           End If
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           ! FOR PION+PION,pion+eta, eta+eta and rho(omega)+pion(rho,omega) or eta
777        Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ! energy thresh for collisions
           ec0 = em1 + em2 + 0.02
           If (srt<=ec0) Goto 400
           ec = (em1+em2+0.02)**2
           ! we negelect the elastic collision between mesons except that betwen
           ! two pions because of the lack of information about these collisions
           ! However, we do let them to collide inelastically to produce kaons
           !lin-8/15/02       ppel=1.e-09
           ppel = 20.
           ipp = 1
           If (lb1<3 .Or. lb1>5 .Or. lb2<3 .Or. lb2>5) Goto 778
           Call ppxs(lb1, lb2, srt, ppsig, spprho, ipp)
           ppel = ppsig
778        ppink = pipik(srt)

           ! pi+eta and eta+eta are assumed to be the same as pipik( for pi+pi -> K+K-)
           ! estimated from Ko's paper:
           ppink = 2.0*ppink
           If (lb1>=25 .And. lb2>=25) ppink = rrkk

           !lin-2/13/03 include omega the same as rho, eta the same as pi:
           !        if(((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.25.and.lb2.le.27))
           !     1  .or.((lb2.ge.3.and.lb2.le.5).and.(lb1.ge.25.and.lb1.le.27)))
           If (((lb1==0 .Or. (lb1>=3 .And. lb1<=5)) .And. (lb2>=25 .And. lb2<=28)) .Or. ((lb2==0 .Or. (lb2>=3 .And. lb2<=5)) .And. (lb1>=25 .And. lb1<=28))) Then
              ppink = 0.
              If (srt>=(aka+aks)) ppink = prkk
           End If

           ! pi pi <-> rho rho:
           Call spprr(lb1, lb2, srt)
           !lin-4/03/02 pi pi <-> eta eta:
           Call sppee(lb1, lb2, srt)
           !lin-4/03/02 pi pi <-> pi eta:
           Call spppe(lb1, lb2, srt)
           !lin-4/03/02 rho pi <-> rho eta:
           Call srpre(lb1, lb2, srt)
           !lin-4/03/02 omega pi <-> omega eta:
           Call sopoe(lb1, lb2, srt)
           !lin-4/03/02 rho rho <-> eta eta:
           Call srree(lb1, lb2, srt)

           ppinnb = 0.
           If (srt>thresh(1)) Then
              Call getnst(srt)
              If (lb1>=3 .And. lb1<=5 .And. lb2>=3 .And. lb2<=5) Then
                 ppinnb = ppbbar(srt)
              Else If ((lb1>=3 .And. lb1<=5 .And. lb2>=25 .And. lb2<=27) .Or. (lb2>=3 .And. lb2<=5 .And. lb1>=25 .And. lb1<=27)) Then
                 ppinnb = prbbar(srt)
              Else If (lb1>=25 .And. lb1<=27 .And. lb2>=25 .And. lb2<=27) Then
                 ppinnb = rrbbar(srt)
              Else If ((lb1>=3 .And. lb1<=5 .And. lb2==28) .Or. (lb2>=3 .And. lb2<=5 .And. lb1==28)) Then
                 ppinnb = pobbar(srt)
              Else If ((lb1>=25 .And. lb1<=27 .And. lb2==28) .Or. (lb2>=25 .And. lb2<=27 .And. lb1==28)) Then
                 ppinnb = robbar(srt)
              Else If (lb1==28 .And. lb2==28) Then
                 ppinnb = oobbar(srt)
              Else
                 If (lb1/=0 .And. lb2/=0) Write (6, *) 'missed MM lb1,lb2=', lb1, lb2
              End If
           End If
           ppin = ppink + ppinnb + pprr + ppee + pppe + rpre + xopoe + rree

           ! check if a collision can happen
           If ((ppel+ppin)<=0.01) Goto 400
           dspp = sqrt((ppel+ppin)/31.4)
           dsppr = dspp + 0.1
           Call distce(i1, i2, dsppr, dspp, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           If (ppel==0) Goto 400
           ! the collision can happen
           ! check what kind collision has happened
           ekaon(5, iss) = ekaon(5, iss) + 1
           Call crpp(px1cm, py1cm, pz1cm, srt, i1, i2, iblock, ppel, ppin, spprho, ipp)

           ! rho formation, go to 400
           !       if(iblock.eq.666)go to 600
           If (iblock==666) Goto 555
           If (iblock==6) lpp = lpp + 1
           If (iblock==66) Then
              lppk = lppk + 1
           Else If (iblock==366) Then
              lppk = lppk + 1
           Else If (iblock==367) Then
              lppk = lppk + 1
           End If
           em1 = e(i1)
           em2 = e(i2)
           Goto 440

           ! In this block we treat annihilations of
           !lin-9/28/00* an anti-nucleon and a baryon or baryon resonance
           ! an anti-baryon and a baryon (including resonances)
2799       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           !lin assume the same cross section (as a function of sqrt s) as for PPbar:

           !lin-ctest annih maximum
           !        DSppb=SQRT(amin1(xppbar(srt),30.)/PI/10.)
           dsppb = sqrt(xppbar(srt)/pi/10.)
           dsppbr = dsppb + 0.1
           Call distce(i1, i2, dsppbr, dsppb, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crppba(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !
3555       px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           dskk = sqrt(sig/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crlaba(px1cm, py1cm, pz1cm, srt, brel, brsgm, i1, i2, nt, iblock, nchrg, icase)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !
           ! perturbative production of cascade and omega
3455       px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           Call pertur(px1cm, py1cm, pz1cm, srt, irun, i1, i2, nt, kp, icontp)
           If (icontp==0) Then
              !     inelastic collisions:
              em1 = e(i1)
              em2 = e(i2)
              iblock = 727
              Goto 440
           End If
           !     elastic collisions:
           If (e(i1)==0.) Goto 800
           If (e(i2)==0.) Goto 600
           Goto 400
           !
           !* phi + N --> pi+N(D),  N(D,N*)+N(D,N*),  K+ +La
           !* phi + D --> pi+N(D)
7222       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call xphib(lb1, lb2, em1, em2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, sigp)
           dskk = sqrt(sigp/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crphib(px1cm, py1cm, pz1cm, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, sigp, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !
           !* phi + M --> K+ + K* .....
7444       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call phimes(i1, i2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, sigphi)
           dskk = sqrt(sigphi/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           !*---
           pzrt = p(3, i1) + p(3, i2)
           er1 = sqrt(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+e(i1)**2)
           er2 = sqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+e(i2)**2)
           ert = er1 + er2
           yy = 0.5*log((ert+pzrt)/(ert-pzrt))
           !*------
           Call crphim(px1cm, py1cm, pz1cm, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, sigphi, ikkg, ikkl, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !
           ! lambda-N elastic xsection, Li & Ko, PRC 54(1996)1897.
7799       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call lambar(i1, i2, srt, siglab)
           dshn = sqrt(siglab/pi/10.)
           dshnr = dshn + 0.1
           Call distce(i1, i2, dshnr, dshn, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crhb(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !
           !* K+ + La(Si) --> Meson + B
           !* K- + La(Si)-bar --> Meson + B-bar
5699       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           Call xkhype(i1, i2, srt, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk)
           dskk = sqrt(sigk/pi)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           !
           If (lb(i1)==23 .Or. lb(i2)==23) Then
              ikmp = 1
           Else
              ikmp = -1
           End If
           Call crkhyp(px1cm, py1cm, pz1cm, srt, i1, i2, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk, ikmp, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           ! khyperon end
           !
           !sp11/03/01 La/Si-bar + N --> pi + K+
           !  La/Si + N-bar --> pi + K-
5999       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           sigkp = 15.
           !      if((lb1.ge.14.and.lb1.le.17)
           !     &    .or.(lb2.ge.14.and.lb2.le.17))sigkp=10.
           dskk = sqrt(sigkp/pi/10.)
           dskk0 = dskk + 0.1
           Call distce(i1, i2, dskk0, dskk, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           !
           Call crlan(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !
           !*
           ! K(K*) + K(K*) --> phi + pi(rho,omega)
8699       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           !  CALL CROSSKKPHI(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)  used for KK*->phi+rho

           Call crkphi(px1cm, py1cm, pz1cm, ec, srt, iblock, emm1, emm2, lbp1, lbp2, i1, i2, ikk, icase, rrkk, prkk)
           If (icase==0) Then
              iblock = 0
              Goto 400
           End If

           !*---
           If (lbp1==29 .Or. lbp2==29) Then
              pzrt = p(3, i1) + p(3, i2)
              er1 = sqrt(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+e(i1)**2)
              er2 = sqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+e(i2)**2)
              ert = er1 + er2
              yy = 0.5*log((ert+pzrt)/(ert-pzrt))
              !*------
              iblock = 222
              ntag = 0
           End If

           lb(i1) = lbp1
           lb(i2) = lbp2
           e(i1) = emm1
           e(i2) = emm2
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !*
           ! rho(omega) + K(K*)  --> phi + K(K*)
8799       Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           !  CALL CROSSKKPHI(PX1CM,PY1CM,PZ1CM,SRT,I1,I2,IBLOCK)  used for KK*->phi+rho
           Call crksph(px1cm, py1cm, pz1cm, ec, srt, emm1, emm2, lbp1, lbp2, i1, i2, ikkg, ikkl, iblock, icase, srhoks)
           If (icase==0) Then
              iblock = 0
              Goto 400
           End If
           !
           If (lbp1==29 .Or. lbp2==20) Then
              !*---
              pzrt = p(3, i1) + p(3, i2)
              er1 = sqrt(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2+e(i1)**2)
              er2 = sqrt(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2+e(i2)**2)
              ert = er1 + er2
              yy = 0.5*log((ert+pzrt)/(ert-pzrt))
           End If

           lb(i1) = lbp1
           lb(i2) = lbp2
           e(i1) = emm1
           e(i2) = emm2
           em1 = e(i1)
           em2 = e(i2)
           Goto 440

           ! for kaon+baryon scattering, using a constant xsection of 10 mb.
888        Continue
           px1cm = pcx
           py1cm = pcy
           pz1cm = pcz
           ec = (em1+em2+0.02)**2
           sig = 10.
           If (iabs(lb1)==14 .Or. iabs(lb2)==14 .Or. iabs(lb1)==30 .Or. iabs(lb2)==30) sig = 20.
           If (lb1==29 .Or. lb2==29) sig = 5.0

           dskn = sqrt(sig/pi/10.)
           dsknr = dskn + 0.1
           Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
           If (ic==-1) Goto 400
           Call crkn(px1cm, py1cm, pz1cm, srt, i1, i2, iblock)
           em1 = e(i1)
           em2 = e(i2)
           Goto 440
           !**

440        Continue
           If (iblock==0) Goto 400
           !COM: FOR DIRECT PROCESS WE HAVE TREATED THE PAULI BLOCKING AND FIND
           !     THE MOMENTUM OF PARTICLES IN THE ''LAB'' FRAME. SO GO TO 400
           ! A COLLISION HAS TAKEN PLACE !!
           lcoll = lcoll + 1
           ! WAS COLLISION PAULI-FORBIDEN? IF YES, NTAG = -1
           ntag = 0
           !
           !             LORENTZ-TRANSFORMATION INTO CMS FRAME
           e1cm = sqrt(em1**2+px1cm**2+py1cm**2+pz1cm**2)
           p1beta = px1cm*betax + py1cm*betay + pz1cm*betaz
           transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
           pt1i1 = betax*transf + px1cm
           pt2i1 = betay*transf + py1cm
           pt3i1 = betaz*transf + pz1cm
           ! negelect the pauli blocking at high energies
           Goto 90002

           !lin-10/25/02-comment out following, since there is no path to it:
           !*CHECK IF PARTICLE #1 IS PAULI BLOCKED
           !              CALL PAULat(I1,occup)
           !              if (RANART(NSEED) .lt. occup) then
           !                ntag = -1
           !              else
           !                ntag = 0
           !              end if
           !lin-10/25/02-end

90002      Continue
           !IF PARTICLE #1 IS NOT PAULI BLOCKED
           !              IF (NTAG .NE. -1) THEN
           e2cm = sqrt(em2**2+px1cm**2+py1cm**2+pz1cm**2)
           transf = gamma*(-gamma*p1beta/(gamma+1.)+e2cm)
           pt1i2 = betax*transf - px1cm
           pt2i2 = betay*transf - py1cm
           pt3i2 = betaz*transf - pz1cm
           Goto 90003


90003      If (iblock==1) lcnne = lcnne + 1
           If (iblock==5) ldd = ldd + 1
           If (iblock==2) lcnnd = lcnnd + 1
           If (iblock==8) lkn = lkn + 1
           If (iblock==43) ldou = ldou + 1
           If (iblock==3) lcndn = lcndn + 1
           p(1, i1) = pt1i1
           p(2, i1) = pt2i1
           p(3, i1) = pt3i1
           p(1, i2) = pt1i2
           p(2, i2) = pt2i2
           p(3, i2) = pt3i2
           px1 = p(1, i1)
           py1 = p(2, i1)
           pz1 = p(3, i1)
           em1 = e(i1)
           em2 = e(i2)
           lb1 = lb(i1)
           lb2 = lb(i2)
           id(i1) = 2
           id(i2) = 2
           e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
           id1 = id(i1)
           Goto 90004

90004      Continue
           am1 = em1
           am2 = em2
           !            END IF


400        Continue
555        Continue

600     End Do
798     If (nt==ntmax .And. ipi0dcy==1 .And. i1==(massr(irun)+msum)) Then
           Do ipion = 1, nnn
              If (lpion(ipion,irun)==4) Then
                 wid = 7.85E-9
                 Call resdec(i1, nt, nnn, wid, idecay, ipion)
              End If
           End Do
        End If

800  End Do
     n0 = mass + msum
     Do n = n0 + 1, massr(irun) + msum
        If (e(n)>0. .Or. lb(n)>5000) Then
           !bz11/25/98end
           nnn = nnn + 1
           rpion(1, nnn, irun) = r(1, n)
           rpion(2, nnn, irun) = r(2, n)
           rpion(3, nnn, irun) = r(3, n)
           !lin-10/28/03:
           If (nt==ntmax) Then
              ftpisv(nnn, irun) = ftsv(n)
              tfdpi(nnn, irun) = tfdcy(n)
           End If
           !
           ppion(1, nnn, irun) = p(1, n)
           ppion(2, nnn, irun) = p(2, n)
           ppion(3, nnn, irun) = p(3, n)
           epion(nnn, irun) = e(n)
           lpion(nnn, irun) = lb(n)
           !       !! sp 12/19/00
           propi(nnn, irun) = proper(n)
           !lin-5/2008:
           dppion(nnn, irun) = dpertp(n)
           !        if(lb(n) .eq. 45)
           !    &   write(*,*)'IN-1  NT,NNN,LB,P ',nt,NNN,lb(n),proper(n)
        End If
     End Do
     massrn(irun) = nnn + mass
     !        write(*,*)'F: NNN,massrn ', nnn,massrn(irun)
  End Do
  ia = 0
  ib = 0
  Do irun = 1, num
     ia = ia + massr(irun-1)
     ib = ib + massrn(irun-1)
     Do ic = 1, massrn(irun)
        ie = ia + ic
        ig = ib + ic
        If (ic<=mass) Then
           rt(1, ig) = r(1, ie)
           rt(2, ig) = r(2, ie)
           rt(3, ig) = r(3, ie)
           !lin-10/28/03:
           If (nt==ntmax) Then
              fttemp(ig) = ftsv(ie)
              tft(ig) = tfdcy(ie)
           End If
           !
           pt(1, ig) = p(1, ie)
           pt(2, ig) = p(2, ie)
           pt(3, ig) = p(3, ie)
           et(ig) = e(ie)
           lt(ig) = lb(ie)
           prot(ig) = proper(ie)
           !lin-5/2008:
           dptemp(ig) = dpertp(ie)
        Else
           i0 = ic - mass
           rt(1, ig) = rpion(1, i0, irun)
           rt(2, ig) = rpion(2, i0, irun)
           rt(3, ig) = rpion(3, i0, irun)
           !lin-10/28/03:
           If (nt==ntmax) Then
              fttemp(ig) = ftpisv(i0, irun)
              tft(ig) = tfdpi(i0, irun)
           End If
           !
           pt(1, ig) = ppion(1, i0, irun)
           pt(2, ig) = ppion(2, i0, irun)
           pt(3, ig) = ppion(3, i0, irun)
           et(ig) = epion(i0, irun)
           lt(ig) = lpion(i0, irun)
           prot(ig) = propi(i0, irun)
           !lin-5/2008:
           dptemp(ig) = dppion(i0, irun)
        End If
     End Do
  End Do
  !
  il = 0
  !lin-10/26/01-hbt:
  !        DO 10002 IRUN=1,NUM
  Do irun = 1, num

     massr(irun) = massrn(irun)
     il = il + massr(irun-1)
     Do im = 1, massr(irun)
        in = il + im
        r(1, in) = rt(1, in)
        r(2, in) = rt(2, in)
        r(3, in) = rt(3, in)
        !lin-10/28/03:
        If (nt==ntmax) Then
           ftsv(in) = fttemp(in)
           tfdcy(in) = tft(in)
        End If
        p(1, in) = pt(1, in)
        p(2, in) = pt(2, in)
        p(3, in) = pt(3, in)
        e(in) = et(in)
        lb(in) = lt(in)
        proper(in) = prot(in)
        !lin-5/2008:
        dpertp(in) = dptemp(in)
        If (lb(in)<1 .Or. lb(in)>2) id(in) = 0
     End Do
  End Do
  !
  Return
End Subroutine relcol

!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine cms(i1, i2, px1cm, py1cm, pz1cm, srt)
  ! PURPOSE : FIND THE MOMENTA OF PARTICLES IN THE CMS OF THE
  !          TWO COLLIDING PARTICLES
  ! VARIABLES :
  !****************************************
  Parameter (maxstr=150001)
  Double Precision px1, py1, pz1, px2, py2, pz2, em1, em2, e1, e2, s, etotal, p1beta, transf, dbetax, dbetay, dbetaz, dgamma, scheck
  Common /bb/p(3, maxstr)
  Common /cc/e(maxstr)
  Common /bg/betax, betay, betaz, gamma
  Save
  px1 = dble(p(1,i1))
  py1 = dble(p(2,i1))
  pz1 = dble(p(3,i1))
  px2 = dble(p(1,i2))
  py2 = dble(p(2,i2))
  pz2 = dble(p(3,i2))
  em1 = dble(e(i1))
  em2 = dble(e(i2))
  e1 = dsqrt(em1**2+px1**2+py1**2+pz1**2)
  e2 = dsqrt(em2**2+px2**2+py2**2+pz2**2)
  s = (e1+e2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2
  If (s<=0) s = 0D0
  srt = sngl(dsqrt(s))
  !LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM
  etotal = e1 + e2
  dbetax = (px1+px2)/etotal
  dbetay = (py1+py2)/etotal
  dbetaz = (pz1+pz2)/etotal
  !lin-9/2012: check argument in sqrt():
  scheck = 1.D0 - dbetax**2 - dbetay**2 - dbetaz**2
  If (scheck<=0D0) Then
     Write (99, *) 'scheck1: ', scheck
     Stop
  End If
  dgamma = 1.D0/dsqrt(scheck)
  !TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM)
  p1beta = px1*dbetax + py1*dbetay + pz1*dbetaz
  transf = dgamma*(dgamma*p1beta/(dgamma+1D0)-e1)
  px1cm = sngl(dbetax*transf+px1)
  py1cm = sngl(dbetay*transf+py1)
  pz1cm = sngl(dbetaz*transf+pz1)
  betax = sngl(dbetax)
  betay = sngl(dbetay)
  betaz = sngl(dbetaz)
  gamma = sngl(dgamma)
  Return
End Subroutine cms
!lin-9/2012-end

!**************************************
Subroutine distce(i1, i2, deltar, ds, dt, ec, srt, ic, px1cm, py1cm, pz1cm)
  ! PURPOSE : CHECK IF THE COLLISION BETWEEN TWO PARTICLES CAN HAPPEN
  !           BY CHECKING
  !                      (1) IF THE DISTANCE BETWEEN THEM IS SMALLER
  !           THAN THE MAXIMUM DISTANCE DETERMINED FROM THE CROSS SECTION.
  !                      (2) IF PARTICLE WILL PASS EACH OTHER WITHIN
  !           TWO HARD CORE RADIUS.
  !                      (3) IF PARTICLES WILL GET CLOSER.
  ! VARIABLES :
  !           IC=1 COLLISION HAPPENED
  !           IC=-1 COLLISION CAN NOT HAPPEN
  !****************************************
  Parameter (maxstr=150001)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /bg/betax, betay, betaz, gamma
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /BG/
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Save
  ic = 0
  x1 = r(1, i1)
  y1 = r(2, i1)
  z1 = r(3, i1)
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  x2 = r(1, i2)
  y2 = r(2, i2)
  z2 = r(3, i2)
  px2 = p(1, i2)
  py2 = p(2, i2)
  pz2 = p(3, i2)
  em1 = e(i1)
  em2 = e(i2)
  e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
  !            IF (ABS(X1-X2) .GT. DELTAR) GO TO 400
  !            IF (ABS(Y1-Y2) .GT. DELTAR) GO TO 400
  !            IF (ABS(Z1-Z2) .GT. DELTAR) GO TO 400
  rsqare = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
  If (rsqare>deltar**2) Goto 400
  !NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER !
  e2 = sqrt(em2**2+px2**2+py2**2+pz2**2)
  s = srt*srt
  If (s<ec) Goto 400
  !NOW THERE IS ENOUGH ENERGY AVAILABLE !
  !LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM
  ! BETAX, BETAY, BETAZ AND GAMMA HAVE BEEN GIVEN IN THE SUBROUTINE CMS
  !TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM)
  p1beta = px1*betax + py1*betay + pz1*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)-e1)
  prcm = sqrt(px1cm**2+py1cm**2+pz1cm**2)
  If (prcm<=0.00001) Goto 400
  !TRANSFORMATION OF SPATIAL DISTANCE
  drbeta = betax*(x1-x2) + betay*(y1-y2) + betaz*(z1-z2)
  transf = gamma*gamma*drbeta/(gamma+1)
  dxcm = betax*transf + x1 - x2
  dycm = betay*transf + y1 - y2
  dzcm = betaz*transf + z1 - z2
  !DETERMINING IF THIS IS THE POINT OF CLOSEST APPROACH
  drcm = sqrt(dxcm**2+dycm**2+dzcm**2)
  dzz = (px1cm*dxcm+py1cm*dycm+pz1cm*dzcm)/prcm
  If ((drcm**2-dzz**2)<=0.) Then
     bbb = 0.
  Else
     bbb = sqrt(drcm**2-dzz**2)
  End If
  !WILL PARTICLE PASS EACH OTHER WITHIN 2 * HARD CORE RADIUS ?
  If (bbb>ds) Goto 400
  relvel = prcm*(1.0/e1+1.0/e2)
  ddd = relvel*dt*0.5
  !WILL PARTICLES GET CLOSER ?
  If (abs(ddd)<abs(dzz)) Goto 400
  ic = 1
  Goto 500
400 ic = -1
500 Continue
  Return
End Subroutine distce
!***************************************
!                                                                      *
!                                                                      *
Subroutine crnn(irun, px, py, pz, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383, aphi=1.020)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (xmd=1.8756, npdmax=10000)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  !c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  !c      SAVE /gg/
  Common /input/nstar, ndirct, dir
  !c      SAVE /INPUT/
  Common /nn/nnn
  !c      SAVE /NN/
  Common /bg/betax, betay, betaz, gamma
  !c      SAVE /BG/
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
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  !c      SAVE /leadng/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /dpi/em2, lb2
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /para8/idpert, npertd, idxsec
  Dimension ppd(3, npdmax), lbpd(npdmax)
  Save
  !-----------------------------------------------------------------------
  n12 = 0
  m12 = 0
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
  c2 = pz/pr
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .And. lb(i2)<0) ianti = 1
  Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  !lin-5/2008 Production of perturbative deuterons for idpert=1:
  If (idpert==1 .And. ipert1==1) Then
     If (srt<2.012) Return
     If ((iabs(lb(i1))==1 .Or. iabs(lb(i1))==2) .And. (iabs(lb(i2))==1 .Or. iabs(lb(i2))==2)) Then
        Goto 108
     Else
        Return
     End If
  End If
  !
  !-----------------------------------------------------------------------
  !COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R
  !      N-DELTA OR N*-N* or N*-Delta)
  !      IF (X1 .LE. SIGNN/SIG) THEN
  If (x1<=(signn/sig)) Then
     !COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER
     as = (3.65*(srt-1.8766))**6
     a = 6.0*as/(1.0+as)
     ta = -2.0*pr**2
     x = ranart(nseed)
     !lin-10/24/02        T1  = DLOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A
     t1 = sngl(dlog(dble(1.-x)*dexp(dble(a)*dble(ta))+dble(x)))/a
     c1 = 1.0 - t1/ta
     t1 = 2.0*pi*ranart(nseed)
     iblock = 1
     Goto 107
  Else
     !COM: TEST FOR INELASTIC SCATTERING
     !     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING
     !     CAN HAPPEN ANY MORE ==> RETURN (2.012 = 2*AVMASS + PI-MASS)
     !lin-5/2008: Mdeuteron+Mpi=2.0106 to 2.0152 GeV/c2, so we can still use this:
     If (srt<2.012) Return
     !     calculate the N*(1535) production cross section in N+N collisions
     !     note that the cross sections in this subroutine are in units of mb
     !     as only ratios of the cross sections are used to determine the
     !     reaction channels
     Call n1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
     !COM: HERE WE HAVE A PROCESS N+N ==> N+DELTA,OR N+N==>N+N*(144) or N*(1535)
     !     OR
     ! 3 pi channel : N+N==>d1+d2+PION
     sig3 = 3.*(x3pi(srt)+x33pi(srt))
     ! 2 pi channel : N+N==>d1+d2+d1*n*+n*n*
     sig4 = 4.*x2pi(srt)
     ! 4 pi channel : N+N==>d1+d2+rho
     s4pi = x4pi(srt)
     ! N+N-->NN+rho channel
     srho = xrho(srt)
     ! N+N-->NN+omega
     somega = omega(srt)
     ! CROSS SECTION FOR KAON PRODUCTION from the four channels
     ! for NLK channel
     akp = 0.498
     ak0 = 0.498
     ana = 0.94
     ada = 1.232
     al = 1.1157
     as = 1.1197
     xsk1 = 0
     xsk2 = 0
     xsk3 = 0
     xsk4 = 0
     xsk5 = 0
     t1nlk = ana + al + akp
     If (srt<=t1nlk) Goto 222
     xsk1 = 1.5*pplpk(srt)
     ! for DLK channel
     t1dlk = ada + al + akp
     t2dlk = ada + al - akp
     If (srt<=t1dlk) Goto 222
     es = srt
     pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
     pmdlk = sqrt(pmdlk2)
     xsk3 = 1.5*pplpk(srt)
     ! for NSK channel
     t1nsk = ana + as + akp
     t2nsk = ana + as - akp
     If (srt<=t1nsk) Goto 222
     pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
     pmnsk = sqrt(pmnsk2)
     xsk2 = 1.5*(ppk1(srt)+ppk0(srt))
     ! for DSK channel
     t1dsk = ada + as + akp
     t2dsk = ada + as - akp
     If (srt<=t1dsk) Goto 222
     pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
     pmdsk = sqrt(pmdsk2)
     xsk4 = 1.5*(ppk1(srt)+ppk0(srt))
     !sp11/21/01
     ! phi production
     If (srt<=(2.*amn+aphi)) Goto 222
     !  !! mb put the correct form
     xsk5 = 0.0001
     !sp11/21/01 end
     !
     ! THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN
222  sigk = xsk1 + xsk2 + xsk3 + xsk4

     !bz3/7/99 neutralk
     xsk1 = 2.0*xsk1
     xsk2 = 2.0*xsk2
     xsk3 = 2.0*xsk3
     xsk4 = 2.0*xsk4
     sigk = 2.0*sigk + xsk5
     !bz3/7/99 neutralk end
     !
     !* FOR P+P or L/S+L/S COLLISION:
     !       lb1=lb(i1)
     !       lb2=lb(i2)
     lb1 = iabs(lb(i1))
     lb2 = iabs(lb(i2))
     If ((lb(i1)*lb(i2)==1) .Or. ((lb1<=17 .And. lb1>=14) .And. (lb2<=17 .And. lb2>=14)) .Or. ((lb1<=2) .And. (lb2<=17 .And. lb2>=14)) .Or. ((lb2<=2) .And. (lb1<=17 .And. lb1>=14))) Then
        !lin-8/2008 PP->d+meson here:
        If (x1<=((signn+sdprod)/sig)) Goto 108
        sig1 = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sig2 = 1.5*sigma(srt, 1, 1, 1)
        signd = sig1 + sig2 + sig3 + sig4 + x1535 + sigk + s4pi + srho + somega
        !lin-5/2008:
        !           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
        If (x1>(signn+signd+sdprod)/sig) Return
        dir = sig3/signd
        If (ranart(nseed)<=dir) Goto 106
        If (ranart(nseed)<=sigk/(sigk+x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 306
        If (ranart(nseed)<=s4pi/(x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 307
        If (ranart(nseed)<=srho/(x1535+sig4+sig2+sig1+srho+somega)) Goto 308
        If (ranart(nseed)<=somega/(x1535+sig4+sig2+sig1+somega)) Goto 309
        If (ranart(nseed)<=x1535/(sig1+sig2+sig4+x1535)) Then
           ! N*(1535) production
           n12 = 9
        Else
           If (ranart(nseed)<=sig4/(sig1+sig2+sig4)) Then
              ! DOUBLE DELTA PRODUCTION
              n12 = 66
              Goto 1012
           Else
              !DELTA PRODUCTION
              n12 = 3
              If (ranart(nseed)>sig1/(sig1+sig2)) n12 = 4
           End If
        End If
        Goto 1011
     End If
     !* FOR N+N COLLISION:
     If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
        !lin-8/2008 NN->d+meson here:
        If (x1<=((signn+sdprod)/sig)) Goto 108
        sig1 = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sig2 = 1.5*sigma(srt, 1, 1, 1)
        signd = sig1 + sig2 + x1535 + sig3 + sig4 + sigk + s4pi + srho + somega
        !lin-5/2008:
        !           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
        If (x1>(signn+signd+sdprod)/sig) Return
        dir = sig3/signd
        If (ranart(nseed)<=dir) Goto 106
        If (ranart(nseed)<=sigk/(sigk+x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 306
        If (ranart(nseed)<=s4pi/(x1535+sig4+sig2+sig1+s4pi+srho+somega)) Goto 307
        If (ranart(nseed)<=srho/(x1535+sig4+sig2+sig1+srho+somega)) Goto 308
        If (ranart(nseed)<=somega/(x1535+sig4+sig2+sig1+somega)) Goto 309
        If (ranart(nseed)<=x1535/(x1535+sig1+sig2+sig4)) Then
           ! N*(1535) PRODUCTION
           n12 = 10
        Else
           If (ranart(nseed)<=sig4/(sig1+sig2+sig4)) Then
              ! double delta production
              n12 = 67
              Goto 1013
           Else
              ! DELTA PRODUCTION
              n12 = 6
              If (ranart(nseed)>sig1/(sig1+sig2)) n12 = 5
           End If
        End If
        Goto 1011
     End If
     !* FOR N+P COLLISION
     If (lb(i1)*lb(i2)==2) Then
        !lin-5/2008 NP->d+meson here:
        If (x1<=((signn+sdprod)/sig)) Goto 108
        sig1 = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
        If (nstar==1) Then
           sig2 = (3./4.)*sigma(srt, 2, 0, 1)
        Else
           sig2 = 0.
        End If
        signd = 2.*(sig1+sig2+x1535) + sig3 + sig4 + sigk + s4pi + srho + somega
        !lin-5/2008:
        !           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
        If (x1>(signn+signd+sdprod)/sig) Return
        dir = sig3/signd
        If (ranart(nseed)<=dir) Goto 106
        If (ranart(nseed)<=sigk/(signd-sig3)) Goto 306
        If (ranart(nseed)<=s4pi/(signd-sig3-sigk)) Goto 307
        If (ranart(nseed)<=srho/(signd-sig3-sigk-s4pi)) Goto 308
        If (ranart(nseed)<=somega/(signd-sig3-sigk-s4pi-srho)) Goto 309
        If (ranart(nseed)<x1535/(sig1+sig2+x1535+0.5*sig4)) Then
           ! N*(1535) PRODUCTION
           n12 = 11
           If (ranart(nseed)<=0.5) n12 = 12
        Else
           If (ranart(nseed)<=sig4/(sig4+2.*(sig1+sig2))) Then
              ! double resonance production
              n12 = 68
              Goto 1014
           Else
              If (ranart(nseed)<=sig1/(sig1+sig2)) Then
                 ! DELTA PRODUCTION
                 n12 = 2
                 If (ranart(nseed)>=0.5) n12 = 1
              Else
                 ! N*(1440) PRODUCTION
                 n12 = 8
                 If (ranart(nseed)>=0.5) n12 = 7
              End If
           End If
        End If
     End If
1011 iblock = 2
     Continue
     !PARAMETRIZATION OF THE SHAPE OF THE DELTA RESONANCE ACCORDING
     !     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER
     !     FORMULA FOR N* RESORANCE
     !     DETERMINE DELTA MASS VIA REJECTION METHOD.
     dmax = srt - avmass - 0.005
     dmax = srt - avmass - 0.005
     dmin = 1.078
     If (n12<7) Then
        ! Delta(1232) production
        If (dmax<1.232) Then
           fm = fde(dmax, srt, 0.)
        Else

           !lin-10/25/02 get rid of argument usage mismatch in FDE():
           xdmass = 1.232
           !          FM=FDE(1.232,SRT,1.)
           fm = fde(xdmass, srt, 1.)
           !lin-10/25/02-end

        End If
        If (fm==0.) fm = 1.E-09
        ntry1 = 0
10      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry1 = ntry1 + 1
        If ((ranart(nseed)>fde(dm,srt,1.)/fm) .And. (ntry1<=30)) Goto 10

        !lin-2/26/03 limit the Delta mass below a certain value
        !     (here taken as its central value + 2* B-W fullwidth):
        If (dm>1.47) Goto 10

        Goto 13
     End If
     If ((n12==7) .Or. (n12==8)) Then
        ! N*(1440) production
        If (dmax<1.44) Then
           fm = fns(dmax, srt, 0.)
        Else

           !lin-10/25/02 get rid of argument usage mismatch in FNS():
           xdmass = 1.44
           !          FM=FNS(1.44,SRT,1.)
           fm = fns(xdmass, srt, 1.)
           !lin-10/25/02-end

        End If
        If (fm==0.) fm = 1.E-09
        ntry2 = 0
11      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry2 = ntry2 + 1
        If ((ranart(nseed)>fns(dm,srt,1.)/fm) .And. (ntry2<=10)) Goto 11

        !lin-2/26/03 limit the N* mass below a certain value
        !     (here taken as its central value + 2* B-W fullwidth):
        If (dm>2.14) Goto 11

        Goto 13
     End If
     If (n12>=17) Then
        ! N*(1535) production
        If (dmax<1.535) Then
           fm = fd5(dmax, srt, 0.)
        Else

           !lin-10/25/02 get rid of argument usage mismatch in FNS():
           xdmass = 1.535
           !          FM=FD5(1.535,SRT,1.)
           fm = fd5(xdmass, srt, 1.)
           !lin-10/25/02-end

        End If
        If (fm==0.) fm = 1.E-09
        ntry1 = 0
12      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry1 = ntry1 + 1
        If ((ranart(nseed)>fd5(dm,srt,1.)/fm) .And. (ntry1<=10)) Goto 12

        !lin-2/26/03 limit the N* mass below a certain value
        !     (here taken as its central value + 2* B-W fullwidth):
        If (dm>1.84) Goto 12

        Goto 13
     End If
     ! CALCULATE THE MASSES OF BARYON RESONANCES IN THE DOUBLE RESONANCE
     ! PRODUCTION PROCESS AND RELABLE THE PARTICLES
1012 iblock = 43
     Call rmasdd(srt, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
     Call rmasdd(srt, 1.232, 1.44, 1.08, 1.08, iseed, 3, dm1n, dm2n)
     If (n12==66) Then
        !(1) PP-->DOUBLE RESONANCES
        ! DETERMINE THE FINAL STATE
        xfinal = ranart(nseed)
        If (xfinal<=0.25) Then
           ! (1.1) D+++D0
           lb(i1) = 9
           lb(i2) = 7
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If ((xfinal>0.25) .And. (xfinal<=0.5)) Then
           ! (1.2) D++D+
           lb(i1) = 8
           lb(i2) = 8
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If ((xfinal>0.5) .And. (xfinal<=0.75)) Then
           ! (1.3) D+++N*0
           lb(i1) = 9
           lb(i2) = 10
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If (xfinal>0.75) Then
           ! (1.4) D++N*+
           lb(i1) = 8
           lb(i2) = 11
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
           ! go to 200 to set the new momentum
        End If
     End If
1013 iblock = 43
     Call rmasdd(srt, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
     Call rmasdd(srt, 1.232, 1.44, 1.08, 1.08, iseed, 3, dm1n, dm2n)
     If (n12==67) Then
        !(2) NN-->DOUBLE RESONANCES
        ! DETERMINE THE FINAL STATE
        xfinal = ranart(nseed)
        If (xfinal<=0.25) Then
           ! (2.1) D0+D0
           lb(i1) = 7
           lb(i2) = 7
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If ((xfinal>0.25) .And. (xfinal<=0.5)) Then
           ! (2.2) D++D+
           lb(i1) = 6
           lb(i2) = 8
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If ((xfinal>0.5) .And. (xfinal<=0.75)) Then
           ! (2.3) D0+N*0
           lb(i1) = 7
           lb(i2) = 10
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If (xfinal>0.75) Then
           ! (2.4) D++N*+
           lb(i1) = 8
           lb(i2) = 11
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
           ! go to 200 to set the new momentum
        End If
     End If
1014 iblock = 43
     Call rmasdd(srt, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
     Call rmasdd(srt, 1.232, 1.44, 1.08, 1.08, iseed, 3, dm1n, dm2n)
     If (n12==68) Then
        !(3) NP-->DOUBLE RESONANCES
        ! DETERMINE THE FINAL STATE
        xfinal = ranart(nseed)
        If (xfinal<=0.25) Then
           ! (3.1) D0+D+
           lb(i1) = 7
           lb(i2) = 8
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If ((xfinal>0.25) .And. (xfinal<=0.5)) Then
           ! (3.2) D+++D-
           lb(i1) = 9
           lb(i2) = 6
           e(i1) = dm1
           e(i2) = dm2
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If ((xfinal>0.5) .And. (xfinal<=0.75)) Then
           ! (3.3) D0+N*+
           lb(i1) = 7
           lb(i2) = 11
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
           ! go to 200 to set the new momentum
        End If
        If (xfinal>0.75) Then
           ! (3.4) D++N*0
           lb(i1) = 8
           lb(i2) = 10
           e(i1) = dm1n
           e(i2) = dm2n
           Goto 200
           ! go to 200 to set the new momentum
        End If
     End If
13   Continue
     !-------------------------------------------------------
     ! RELABLE BARYON I1 AND I2
     !1. p+n-->delta(+)+n
     If (n12==1) Then
        If (iabs(lb(i1))==1) Then
           lb(i2) = 2
           lb(i1) = 8
           e(i1) = dm
        Else
           lb(i1) = 2
           lb(i2) = 8
           e(i2) = dm
        End If
        Goto 200
     End If
     !2 p+n-->delta(0)+p
     If (n12==2) Then
        If (iabs(lb(i1))==2) Then
           lb(i2) = 1
           lb(i1) = 7
           e(i1) = dm
        Else
           lb(i1) = 1
           lb(i2) = 7
           e(i2) = dm
        End If
        Goto 200
     End If
     !3 p+p-->delta(++)+n
     If (n12==3) Then
        lb(i1) = 9
        e(i1) = dm
        lb(i2) = 2
        e(i2) = amn
        Goto 200
     End If
     !4 p+p-->delta(+)+p
     If (n12==4) Then
        lb(i2) = 1
        lb(i1) = 8
        e(i1) = dm
        Goto 200
     End If
     !5 n+n--> delta(0)+n
     If (n12==5) Then
        lb(i2) = 2
        lb(i1) = 7
        e(i1) = dm
        Goto 200
     End If
     !6 n+n--> delta(-)+p
     If (n12==6) Then
        lb(i1) = 6
        e(i1) = dm
        lb(i2) = 1
        e(i2) = amp
        Goto 200
     End If
     !7 n+p--> N*(0)+p
     If (n12==7) Then
        If (iabs(lb(i1))==1) Then
           lb(i1) = 1
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 1
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 200
     End If
     !8 n+p--> N*(+)+n
     If (n12==8) Then
        If (iabs(lb(i1))==1) Then
           lb(i2) = 2
           lb(i1) = 11
           e(i1) = dm
        Else
           lb(i1) = 2
           lb(i2) = 11
           e(i2) = dm
        End If
        Goto 200
     End If
     !9 p+p--> N*(+)(1535)+p
     If (n12==9) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           lb(i1) = 13
           e(i1) = dm
        Else
           lb(i1) = 1
           lb(i2) = 13
           e(i2) = dm
        End If
        Goto 200
     End If
     !10 n+n--> N*(0)(1535)+n
     If (n12==10) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 2
           lb(i1) = 12
           e(i1) = dm
        Else
           lb(i1) = 2
           lb(i2) = 12
           e(i2) = dm
        End If
        Goto 200
     End If
     !11 n+p--> N*(+)(1535)+n
     If (n12==11) Then
        If (iabs(lb(i1))==2) Then
           lb(i1) = 2
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 200
     End If
     !12 n+p--> N*(0)(1535)+p
     If (n12==12) Then
        If (iabs(lb(i1))==1) Then
           lb(i1) = 1
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           lb(i1) = 12
           e(i1) = dm
        End If
     End If
  End If
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
200 em1 = e(i1)
  em2 = e(i2)
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  If (srt<=2.14) c1 = 1.0 - 2.0*ranart(nseed)
  If (srt>2.14 .And. srt<=2.4) c1 = ang(srt, iseed)
  If (srt>2.4) Then

     !lin-10/25/02 get rid of argument usage mismatch in PTR():
     xptr = 0.33*pr
     !         cc1=ptr(0.33*pr,iseed)
     cc1 = ptr(xptr, iseed)
     !lin-10/25/02-end

     !lin-9/2012: check argument in sqrt():
     scheck = pr**2 - cc1**2
     If (scheck<0) Then
        Write (99, *) 'scheck2: ', scheck
        scheck = 0.
     End If
     c1 = sqrt(scheck)/pr
     !             c1=sqrt(pr**2-cc1**2)/pr

  End If
  t1 = 2.0*pi*ranart(nseed)
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
  End If
  Goto 107
  !FOR THE NN-->D1+D2+PI PROCESS, FIND MOMENTUM OF THE FINAL TWO
  !DELTAS AND PION IN THE NUCLEUS-NUCLEUS CMS.
106 Continue
  ntry1 = 0
123 Call ddp2(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=40)) Goto 123
  !       if(icou1.lt.0)return
  ! ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  ! DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  ! (1) FOR P+P
  xdir = ranart(nseed)
  If (lb(i1)*lb(i2)==1) Then
     If (xdir<=0.2) Then
        ! (1.1)P+P-->D+++D0+PION(0)
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 9
        lb(i2) = 7
        Goto 205
     End If
     ! (1.2)P+P -->D++D+PION(0)
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 8
        lb(i2) = 8
        Goto 205
     End If
     ! (1.3)P+P-->D+++D+PION(-)
     If ((xdir<=0.6) .And. (xdir>0.4)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 9
        lb(i2) = 8
        Goto 205
     End If
     If ((xdir<=0.8) .And. (xdir>0.6)) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 9
        lb(i2) = 6
        Goto 205
     End If
     If (xdir>0.8) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
  End If
  ! (2)FOR N+N
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     If (xdir<=0.2) Then
        ! (2.1)N+N-->D++D-+PION(0)
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 6
        lb(i2) = 7
        Goto 205
     End If
     ! (2.2)N+N -->D+++D-+PION(-)
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 6
        lb(i2) = 9
        Goto 205
     End If
     ! (2.3)P+P-->D0+D-+PION(+)
     If ((xdir>0.4) .And. (xdir<=0.6)) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 9
        lb(i2) = 8
        Goto 205
     End If
     ! (2.4)P+P-->D0+D0+PION(0)
     If ((xdir>0.6) .And. (xdir<=0.8)) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 7
        lb(i2) = 7
        Goto 205
     End If
     ! (2.5)P+P-->D0+D++PION(-)
     If (xdir>0.8) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
  End If
  ! (3)FOR N+P
  If (lb(i1)*lb(i2)==2) Then
     If (xdir<=0.17) Then
        ! (3.1)N+P-->D+++D-+PION(0)
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lb(i1) = 6
        lb(i2) = 9
        Goto 205
     End If
     ! (3.2)N+P -->D+++D0+PION(-)
     If ((xdir<=0.34) .And. (xdir>0.17)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 9
        Goto 205
     End If
     ! (3.3)N+P-->D++D-+PION(+)
     If ((xdir>0.34) .And. (xdir<=0.51)) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
     ! (3.4)N+P-->D++D++PION(-)
     If ((xdir>0.51) .And. (xdir<=0.68)) Then
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lb(i1) = 8
        lb(i2) = 8
        Goto 205
     End If
     ! (3.5)N+P-->D0+D++PION(0)
     If ((xdir>0.68) .And. (xdir<=0.85)) Then
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 8
        Goto 205
     End If
     ! (3.6)N+P-->D0+D0+PION(+)
     If (xdir>0.85) Then
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lb(i1) = 7
        lb(i2) = 7
     End If
  End If
  ! FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
  ! NUCLEUS CMS. FRAME
  !             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
205 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  !
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==3) Then
        lpion(nnn, irun) = 5
     Else If (lpion(nnn,irun)==5) Then
        lpion(nnn, irun) = 3
     End If
  End If
  !
  lb1 = lb(i1)
  ! FOR DELTA2
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  ! assign delta1 and delta2 to i1 or i2 to keep the leadng particle
  ! behaviour
  !              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
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
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 4
  ! GET PION'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  !lin-5/2008 do not allow smearing in position of produced particles
  !     to avoid immediate reinteraction with the particle I1, I2 or themselves:
  !2002        X01 = 1.0 - 2.0 * RANART(NSEED)
  !            Y01 = 1.0 - 2.0 * RANART(NSEED)
  !            Z01 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2002
  !                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
  !                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
  !                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  !
  Goto 90005
  !lin-5/2008 N+N->Deuteron+pi:
  !     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
108 Continue
  If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
     !     For idpert=1: we produce npertd pert deuterons:
     ndloop = npertd
  Else If (idpert==2 .And. npertd>=1) Then
     !     For idpert=2: we first save information for npertd pert deuterons;
     !     at the last ndloop we create the regular deuteron+pi
     !     and those pert deuterons:
     ndloop = npertd + 1
  Else
     !     Just create the regular deuteron+pi:
     ndloop = 1
  End If
  !
  dprob1 = sdprod/sig/float(npertd)
  Do idloop = 1, ndloop
     Call bbdangle(pxd, pyd, pzd, nt, ipert1, ianti, idloop, pfinal, dprob1, lbm)
     Call rotate(px, py, pz, pxd, pyd, pzd)
     !     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE
     !     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME:
     !     For the Deuteron:
     xmass = xmd
     e1dcm = sqrt(xmass**2+pxd**2+pyd**2+pzd**2)
     p1dbeta = pxd*betax + pyd*betay + pzd*betaz
     transf = gamma*(gamma*p1dbeta/(gamma+1.)+e1dcm)
     pxi1 = betax*transf + pxd
     pyi1 = betay*transf + pyd
     pzi1 = betaz*transf + pzd
     If (ianti==0) Then
        lbd = 42
     Else
        lbd = -42
     End If
     If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
        !ccc  Perturbative production for idpert=1:
        nnn = nnn + 1
        ppion(1, nnn, irun) = pxi1
        ppion(2, nnn, irun) = pyi1
        ppion(3, nnn, irun) = pzi1
        epion(nnn, irun) = xmd
        lpion(nnn, irun) = lbd
        rpion(1, nnn, irun) = r(1, i1)
        rpion(2, nnn, irun) = r(2, i1)
        rpion(3, nnn, irun) = r(3, i1)
        !lin-5/2008 assign the perturbative probability:
        dppion(nnn, irun) = sdprod/sig/float(npertd)
     Else If (idpert==2 .And. idloop<=npertd) Then
        !lin-5/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons
        !     only when a regular (anti)deuteron+pi is produced in NN collisions.
        !     First save the info for the perturbative deuterons:
        ppd(1, idloop) = pxi1
        ppd(2, idloop) = pyi1
        ppd(3, idloop) = pzi1
        lbpd(idloop) = lbd
     Else
        !ccc  Regular production:
        !     For the regular pion: do LORENTZ-TRANSFORMATION:
        e(i1) = xmm
        e2picm = sqrt(xmm**2+pxd**2+pyd**2+pzd**2)
        p2pibeta = -pxd*betax - pyd*betay - pzd*betaz
        transf = gamma*(gamma*p2pibeta/(gamma+1.)+e2picm)
        pxi2 = betax*transf - pxd
        pyi2 = betay*transf - pyd
        pzi2 = betaz*transf - pzd
        p(1, i1) = pxi2
        p(2, i1) = pyi2
        p(3, i1) = pzi2
        !     Remove regular pion to check the equivalence
        !     between the perturbative and regular deuteron results:
        !                 E(i1)=0.
        !
        lb(i1) = lbm
        px1 = p(1, i1)
        py1 = p(2, i1)
        pz1 = p(3, i1)
        em1 = e(i1)
        id(i1) = 2
        id1 = id(i1)
        e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
        lb1 = lb(i1)
        !     For the regular deuteron:
        p(1, i2) = pxi1
        p(2, i2) = pyi1
        p(3, i2) = pzi1
        lb(i2) = lbd
        lb2 = lb(i2)
        e(i2) = xmd
        eti2 = e(i2)
        id(i2) = 2
        !     For idpert=2: create the perturbative deuterons:
        If (idpert==2 .And. idloop==ndloop) Then
           Do ipertd = 1, npertd
              nnn = nnn + 1
              ppion(1, nnn, irun) = ppd(1, ipertd)
              ppion(2, nnn, irun) = ppd(2, ipertd)
              ppion(3, nnn, irun) = ppd(3, ipertd)
              epion(nnn, irun) = xmd
              lpion(nnn, irun) = lbpd(ipertd)
              rpion(1, nnn, irun) = r(1, i1)
              rpion(2, nnn, irun) = r(2, i1)
              rpion(3, nnn, irun) = r(3, i1)
              !lin-5/2008 assign the perturbative probability:
              dppion(nnn, irun) = 1./float(npertd)
           End Do
        End If
     End If
  End Do
  iblock = 501
  Goto 90005
  !lin-5/2008 N+N->Deuteron+pi over
  ! FOR THE NN-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN
  ! THE NUCLEUS-NUCLEUS CMS.
306 Continue
  !sp11/21/01 phi production
  If (xsk5/sigk>ranart(nseed)) Then
     pz1 = p(3, i1)
     pz2 = p(3, i2)
     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 1 + int(2*ranart(nseed))
     nnn = nnn + 1
     lpion(nnn, irun) = 29
     epion(nnn, irun) = aphi
     iblock = 222
     Goto 208
  End If
  !
  iblock = 9
  If (ianti==1) iblock = -9
  !
  pz1 = p(3, i1)
  pz2 = p(3, i2)
  ! DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  nnn = nnn + 1
  lpion(nnn, irun) = 23
  epion(nnn, irun) = aka
  If (srt<=2.63) Then
     ! only lambda production is possible
     ! (1.1)P+P-->p+L+kaon+
     ic = 1
     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 14
     Goto 208
  End If
  If (srt<=2.74 .And. srt>2.63) Then
     ! both Lambda and sigma production are possible
     If (xsk1/(xsk1+xsk2)>ranart(nseed)) Then
        ! lambda production
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
     Else
        ! sigma production
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 15 + int(3*ranart(nseed))
        ic = 2
     End If
     Goto 208
  End If
  If (srt<=2.77 .And. srt>2.74) Then
     ! then pp-->Delta lamda kaon can happen
     If (xsk1/(xsk1+xsk2+xsk3)>ranart(nseed)) Then
        ! * (1.1)P+P-->p+L+kaon+
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk2/(xsk2+xsk3)>ranart(nseed)) Then
           ! pp-->psk
           ic = 2
           lb(i1) = 1 + int(2*ranart(nseed))
           lb(i2) = 15 + int(3*ranart(nseed))
        Else
           ! pp-->D+l+k
           ic = 3
           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
        End If
        Goto 208
     End If
  End If
  If (srt>2.77) Then
     ! all four channels are possible
     If (xsk1/(xsk1+xsk2+xsk3+xsk4)>ranart(nseed)) Then
        ! p lambda k production
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk3/(xsk2+xsk3+xsk4)>ranart(nseed)) Then
           ! delta l K production
           ic = 3
           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
           Goto 208
        Else
           If (xsk2/(xsk2+xsk4)>ranart(nseed)) Then
              ! n sigma k production
              lb(i1) = 1 + int(2*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))
              ic = 2
           Else
              ic = 4
              lb(i1) = 6 + int(4*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))
           End If
           Goto 208
        End If
     End If
  End If
208 Continue
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==23) lpion(nnn, irun) = 21
  End If
  ! KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE
  ntry1 = 0
127 Call bbkaon(ic, srt, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 127
  !       if(icou1.lt.0)return
  ! ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  ! FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
  ! NUCLEUS CMS. FRAME
  ! (1) for the necleon/delta
  !             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
  e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  lbi1 = lb(i1)
  ! (2) for the lambda/sigma
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lbi2 = lb(i2)
  ! GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
  epcm = sqrt(aka**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  !lin-5/2008
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  !lin-5/2008
  !2003        X01 = 1.0 - 2.0 * RANART(NSEED)
  !            Y01 = 1.0 - 2.0 * RANART(NSEED)
  !            Z01 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2003
  !                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
  !                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
  !                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  !
  ! assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the
  ! leadng particle behaviour
  !              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lbi1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lbi2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  Goto 90005
  ! FOR THE NN-->Delta+Delta+rho PROCESS, FIND MOMENTUM OF THE FINAL
  ! PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
307 Continue
  ntry1 = 0
125 Call ddrho(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, amrho, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 125
  !       if(icou1.lt.0)return
  ! ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  arho = amrho
  ! DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  ! (1) FOR P+P
  xdir = ranart(nseed)
  If (lb(i1)*lb(i2)==1) Then
     If (xdir<=0.2) Then
        ! (1.1)P+P-->D+++D0+rho(0)
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 7
        Goto 2051
     End If
     ! (1.2)P+P -->D++D+rho(0)
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 8
        lb(i2) = 8
        Goto 2051
     End If
     ! (1.3)P+P-->D+++D+arho(-)
     If ((xdir<=0.6) .And. (xdir>0.4)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 8
        Goto 2051
     End If
     If ((xdir<=0.8) .And. (xdir>0.6)) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 6
        Goto 2051
     End If
     If (xdir>0.8) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
  End If
  ! (2)FOR N+N
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     If (xdir<=0.2) Then
        ! (2.1)N+N-->D++D-+rho(0)
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 6
        lb(i2) = 7
        Goto 2051
     End If
     ! (2.2)N+N -->D+++D-+rho(-)
     If ((xdir<=0.4) .And. (xdir>0.2)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 6
        lb(i2) = 9
        Goto 2051
     End If
     ! (2.3)P+P-->D0+D-+rho(+)
     If ((xdir>0.4) .And. (xdir<=0.6)) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 9
        lb(i2) = 8
        Goto 2051
     End If
     ! (2.4)P+P-->D0+D0+rho(0)
     If ((xdir>0.6) .And. (xdir<=0.8)) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 7
        Goto 2051
     End If
     ! (2.5)P+P-->D0+D++rho(-)
     If (xdir>0.8) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
  End If
  ! (3)FOR N+P
  If (lb(i1)*lb(i2)==2) Then
     If (xdir<=0.17) Then
        ! (3.1)N+P-->D+++D-+rho(0)
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 6
        lb(i2) = 9
        Goto 2051
     End If
     ! (3.2)N+P -->D+++D0+rho(-)
     If ((xdir<=0.34) .And. (xdir>0.17)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 9
        Goto 2051
     End If
     ! (3.3)N+P-->D++D-+rho(+)
     If ((xdir>0.34) .And. (xdir<=0.51)) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
     ! (3.4)N+P-->D++D++rho(-)
     If ((xdir>0.51) .And. (xdir<=0.68)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 8
        lb(i2) = 8
        Goto 2051
     End If
     ! (3.5)N+P-->D0+D++rho(0)
     If ((xdir>0.68) .And. (xdir<=0.85)) Then
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 8
        Goto 2051
     End If
     ! (3.6)N+P-->D0+D0+rho(+)
     If (xdir>0.85) Then
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 7
        lb(i2) = 7
     End If
  End If
  ! FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
  ! NUCLEUS CMS. FRAME
  !             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
2051 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  !
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==25) Then
        lpion(nnn, irun) = 27
     Else If (lpion(nnn,irun)==27) Then
        lpion(nnn, irun) = 25
     End If
  End If
  !
  lb1 = lb(i1)
  ! FOR DELTA2
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  ! assign delta1 and delta2 to i1 or i2 to keep the leadng particle
  ! behaviour
  !              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
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
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 44
  ! GET rho'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  !lin-5/2008:
  !2004        X01 = 1.0 - 2.0 * RANART(NSEED)
  !            Y01 = 1.0 - 2.0 * RANART(NSEED)
  !            Z01 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2004
  !                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
  !                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
  !                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  !
  Goto 90005
  ! FOR THE NN-->N+N+rho PROCESS, FIND MOMENTUM OF THE FINAL
  ! PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
308 Continue
  ntry1 = 0
126 Call pprho(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, amrho, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 126
  !       if(icou1.lt.0)return
  ! ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  arho = amrho
  ! DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  ! (1) FOR P+P
  xdir = ranart(nseed)
  If (lb(i1)*lb(i2)==1) Then
     If (xdir<=0.5) Then
        ! (1.1)P+P-->P+P+rho(0)
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 1
        Goto 2052
     Else
        ! (1.2)P+P -->p+n+rho(+)
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 2
        Goto 2052
     End If
  End If
  ! (2)FOR N+N
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     If (xdir<=0.5) Then
        ! (2.1)N+N-->N+N+rho(0)
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 2
        lb(i2) = 2
        Goto 2052
     Else
        ! (2.2)N+N -->N+P+rho(-)
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 2
        Goto 2052
     End If
  End If
  ! (3)FOR N+P
  If (lb(i1)*lb(i2)==2) Then
     If (xdir<=0.33) Then
        ! (3.1)N+P-->N+P+rho(0)
        lpion(nnn, irun) = 26
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 2
        Goto 2052
        ! (3.2)N+P -->P+P+rho(-)
     Else If ((xdir<=0.67) .And. (xdir>0.34)) Then
        lpion(nnn, irun) = 25
        epion(nnn, irun) = arho
        lb(i1) = 1
        lb(i2) = 1
        Goto 2052
     Else
        ! (3.3)N+P-->N+N+rho(+)
        lpion(nnn, irun) = 27
        epion(nnn, irun) = arho
        lb(i1) = 2
        lb(i2) = 2
        Goto 2052
     End If
  End If
  ! FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
  ! NUCLEUS CMS. FRAME
  !             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
2052 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  !
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==25) Then
        lpion(nnn, irun) = 27
     Else If (lpion(nnn,irun)==27) Then
        lpion(nnn, irun) = 25
     End If
  End If
  !
  lb1 = lb(i1)
  ! FOR p2
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  ! assign p1 and p2 to i1 or i2 to keep the leadng particle
  ! behaviour
  !              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
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
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 45
  ! GET rho'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  !lin-5/2008:
  !2005        X01 = 1.0 - 2.0 * RANART(NSEED)
  !            Y01 = 1.0 - 2.0 * RANART(NSEED)
  !            Z01 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2005
  !                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
  !                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
  !                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  !
  Goto 90005
  ! FOR THE NN-->p+p+omega PROCESS, FIND MOMENTUM OF THE FINAL
  ! PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
309 Continue
  ntry1 = 0
138 Call ppomga(srt, iseed, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 138
  !       if(icou1.lt.0)return
  ! ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  nnn = nnn + 1
  aomega = 0.782
  ! DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  ! (1) FOR P+P
  If (lb(i1)*lb(i2)==1) Then
     ! (1.1)P+P-->P+P+omega(0)
     lpion(nnn, irun) = 28
     epion(nnn, irun) = aomega
     lb(i1) = 1
     lb(i2) = 1
     Goto 2053
  End If
  ! (2)FOR N+N
  If (iabs(lb(i1))==2 .And. iabs(lb(i2))==2) Then
     ! (2.1)N+N-->N+N+omega(0)
     lpion(nnn, irun) = 28
     epion(nnn, irun) = aomega
     lb(i1) = 2
     lb(i2) = 2
     Goto 2053
  End If
  ! (3)FOR N+P
  If (lb(i1)*lb(i2)==2) Then
     ! (3.1)N+P-->N+P+omega(0)
     lpion(nnn, irun) = 28
     epion(nnn, irun) = aomega
     lb(i1) = 1
     lb(i2) = 2
     Goto 2053
  End If
  ! FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
  ! NUCLEUS CMS. FRAME
  !             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
2053 e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
  End If
  lb1 = lb(i1)
  ! FOR DELTA2
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  lb2 = lb(i2)
  ! assign delta1 and delta2 to i1 or i2 to keep the leadng particle
  ! behaviour
  !              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
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
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  iblock = 46
  ! GET omega'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
  epcm = sqrt(epion(nnn,irun)**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  !lin-5/2008:
  !2006        X01 = 1.0 - 2.0 * RANART(NSEED)
  !            Y01 = 1.0 - 2.0 * RANART(NSEED)
  !            Z01 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2006
  !                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
  !                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
  !                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  !
  Goto 90005
  ! change phase space density FOR NUCLEONS AFTER THE PROCESS

  !lin-10/25/02-comment out following, since there is no path to it:
  !lin-8/16/02 used before set
  !     IX1,IY1,IZ1,IPX1,IPY1,IPZ1, IX2,IY2,IZ2,IPX2,IPY2,IPZ2:
  !                if ((abs(ix1).le.mx) .and. (abs(iy1).le.my) .and.
  !     &              (abs(iz1).le.mz)) then
  !                  ipx1p = nint(p(1,i1)/dpx)
  !                  ipy1p = nint(p(2,i1)/dpy)
  !                  ipz1p = nint(p(3,i1)/dpz)
  !                  if ((ipx1p.ne.ipx1) .or. (ipy1p.ne.ipy1) .or.
  !     &                (ipz1p.ne.ipz1)) then
  !                    if ((abs(ipx1).le.mpx) .and. (abs(ipy1).le.my)
  !     &                .and. (ipz1.ge.-mpz) .and. (ipz1.le.mpzp))
  !     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) =
  !     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) - 1.
  !                    if ((abs(ipx1p).le.mpx) .and. (abs(ipy1p).le.my)
  !     &                .and. (ipz1p.ge.-mpz).and. (ipz1p.le.mpzp))
  !     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) =
  !     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) + 1.
  !                  end if
  !                end if
  !                if ((abs(ix2).le.mx) .and. (abs(iy2).le.my) .and.
  !     &              (abs(iz2).le.mz)) then
  !                  ipx2p = nint(p(1,i2)/dpx)
  !                  ipy2p = nint(p(2,i2)/dpy)
  !                  ipz2p = nint(p(3,i2)/dpz)
  !                  if ((ipx2p.ne.ipx2) .or. (ipy2p.ne.ipy2) .or.
  !     &                (ipz2p.ne.ipz2)) then
  !                    if ((abs(ipx2).le.mpx) .and. (abs(ipy2).le.my)
  !     &                .and. (ipz2.ge.-mpz) .and. (ipz2.le.mpzp))
  !     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) =
  !     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) - 1.
  !                    if ((abs(ipx2p).le.mpx) .and. (abs(ipy2p).le.my)
  !     &                .and. (ipz2p.ge.-mpz) .and. (ipz2p.le.mpzp))
  !     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) =
  !     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) + 1.
  !                  end if
  !                end if
  !lin-10/25/02-end

90005 Continue
  Return
  !-----------------------------------------------------------------------
  !COM: SET THE NEW MOMENTUM COORDINATES
107 If (px==0.0 .And. py==0.0) Then
     t2 = 0.0
  Else
     t2 = atan2(py, px)
  End If
  s1 = 1.0 - c1**2
  If (s1<=0) s1 = 0
  s1 = sqrt(s1)

  !lin-9/2012: check argument in sqrt():
  scheck = 1.0 - c2**2
  If (scheck<0) Then
     Write (99, *) 'scheck3: ', scheck
     scheck = 0.
  End If
  s2 = sqrt(scheck)
  !       S2  =  SQRT( 1.0 - C2**2 )

  ct1 = cos(t1)
  st1 = sin(t1)
  ct2 = cos(t2)
  st2 = sin(t2)
  pz = pr*(c1*c2-s1*s2*ct1)
  ss = c2*s1*ct1 + s2*c1
  px = pr*(ss*ct2-s1*st1*st2)
  py = pr*(ss*st2+s1*st1*ct2)
  Return
End Subroutine crnn
!lin-5/2008 CRNN over


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine crpp(px, py, pz, srt, i1, i2, iblock, ppel, ppin, spprho, ipp)
  !     PURPOSE:                                                         *
  !             DEALING WITH PION-PION COLLISIONS                        *
  !     NOTE   :                                                         *
  !           VALID ONLY FOR PION-PION-DISTANCES LESS THAN 2.5 FM        *
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                     6-> Meson+Meson elastic
  !                     66-> Meson+meson-->K+K-
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
  !c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
  !c      SAVE /ppmm/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save

  lb1i = lb(i1)
  lb2i = lb(i2)

  px0 = px
  py0 = py
  pz0 = pz
  iblock = 1
  !-----------------------------------------------------------------------
  ! check Meson+Meson inelastic collisions
  !lin-9/28/00
  !        if((srt.gt.1.).and.(ppin/(ppin+ppel).gt.RANART(NSEED)))then
  !        iblock=66
  !        e(i1)=0.498
  !        e(i2)=0.498
  !        lb(i1)=21
  !        lb(i2)=23
  !        go to 10
  !lin-11/07/00
  !        if(srt.gt.1.and.(ppin/(ppin+ppel)).gt.RANART(NSEED)) then
  !lin-4/03/02
  If (srt>(2*aka) .And. (ppin/(ppin+ppel))>ranart(nseed)) Then
     !        if(ppin/(ppin+ppel).gt.RANART(NSEED)) then
     !lin-10/08/00

     ranpi = ranart(nseed)
     If ((pprr/ppin)>=ranpi) Then

        !     1) pi pi <-> rho rho:
        Call pi2ro2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)

        !lin-4/03/02 eta equilibration:
     Else If ((pprr+ppee)/ppin>=ranpi) Then
        !     4) pi pi <-> eta eta:
        Call pi2et2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe)/ppin)>=ranpi) Then
        !     5) pi pi <-> pi eta:
        Call pi3eta(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe+rpre)/ppin)>=ranpi) Then
        !     6) rho pi <-> pi eta:
        Call rpiret(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe+rpre+xopoe)/ppin)>=ranpi) Then
        !     7) omega pi <-> omega eta:
        Call opioet(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
     Else If (((pprr+ppee+pppe+rpre+xopoe+rree)/ppin)>=ranpi) Then
        !     8) rho rho <-> eta eta:
        Call ro2et2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
        !lin-4/03/02-end

        !     2) BBbar production:
     Else If (((pprr+ppee+pppe+rpre+xopoe+rree+ppinnb)/ppin)>=ranpi) Then

        Call bbarfs(lbb1, lbb2, ei1, ei2, iblock, iseed)
        !     3) KKbar production:
     Else
        iblock = 66
        ei1 = aka
        ei2 = aka
        lbb1 = 21
        lbb2 = 23
        !lin-11/07/00 pi rho -> K* Kbar and K*bar K productions:
        lb1 = lb(i1)
        lb2 = lb(i2)
        !lin-2/13/03 include omega the same as rho, eta the same as pi:
        !        if(((lb1.ge.3.and.lb1.le.5).and.(lb2.ge.25.and.lb2.le.27))
        !     1  .or.((lb2.ge.3.and.lb2.le.5).and.(lb1.ge.25.and.lb1.le.27)))
        If (((lb1==0 .Or. (lb1>=3 .And. lb1<=5)) .And. (lb2>=25 .And. lb2<=28)) .Or. ((lb2==0 .Or. (lb2>=3 .And. lb2<=5)) .And. (lb1>=25 .And. lb1<=28))) Then
           ei1 = aks
           ei2 = aka
           If (ranart(nseed)>=0.5) Then
              iblock = 366
              lbb1 = 30
              lbb2 = 21
           Else
              iblock = 367
              lbb1 = -30
              lbb2 = 23
           End If
        End If
        !lin-11/07/00-end
     End If
     !lin-ppbar-8/25/00
     e(i1) = ei1
     e(i2) = ei2
     lb(i1) = lbb1
     lb(i2) = lbb2
     !lin-10/08/00-end

  Else
     !bzdbg10/15/99
     !.....for meson+meson elastic srt.le.2Mk, if not pi+pi collision return
     If ((lb(i1)<3 .Or. lb(i1)>5) .And. (lb(i2)<3 .Or. lb(i2)>5)) Return
     !bzdbg10/15/99 end

     ! check Meson+Meson elastic collisions
     iblock = 6
     ! direct process
     If (ipp==1 .Or. ipp==4 .Or. ipp==6) Goto 10
     If (spprho/ppel>ranart(nseed)) Goto 20
  End If
10 ntag = 0
  em1 = e(i1)
  em2 = e(i2)

  !-----------------------------------------------------------------------
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  ! for isotropic distribution no need to ROTATE THE MOMENTUM

  ! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)

  Return
20 Continue
  iblock = 666
  ! treat rho formation in pion+pion collisions
  ! calculate the mass and momentum of rho in the nucleus-nucleus frame
  Call rhores(i1, i2)
  If (ipp==2) lb(i1) = 27
  If (ipp==3) lb(i1) = 26
  If (ipp==5) lb(i1) = 25
  Return
End Subroutine crpp
!*********************************
!*********************************
!                                                                      *
!                                                                      *
Subroutine crnd(irun, px, py, pz, srt, i1, i2, iblock, signn, sig, sigk, xsk1, xsk2, xsk3, xsk4, xsk5, nt, ipert1)
  !     PURPOSE:                                                         *
  !             DEALING WITH NUCLEON-BARYON RESONANCE COLLISIONS         *
  !     NOTE   :                                                         *
  !           VALID ONLY FOR BARYON-BARYON-DISTANCES LESS THAN 1.32 FM   *
  !           (1.32 = 2 * HARD-CORE-RADIUS [HRC] )                       *
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   *
  !           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                      0-> COLLISION CANNOT HAPPEN                     *
  !                      1-> N-N ELASTIC COLLISION                       *
  !                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          *
  !                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          *
  !                      4-> N+N->N+N+PION,DIRTCT PROCESS                *
  !           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      *
  !                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    *
  !                      N12,                                            *
  !                      M12=1 FOR p+n-->delta(+)+ n                     *
  !                          2     p+n-->delta(0)+ p                     *
  !                          3     p+p-->delta(++)+n                     *
  !                          4     p+p-->delta(+)+p                      *
  !                          5     n+n-->delta(0)+n                      *
  !                          6     n+n-->delta(-)+p                      *
  !                          7     n+p-->N*(0)(1440)+p                   *
  !                          8     n+p-->N*(+)(1440)+n                   *
  !                        9     p+p-->N*(+)(1535)+p                     *
  !                        10    n+n-->N*(0)(1535)+n                     *
  !                         11    n+p-->N*(+)(1535)+n                     *
  !                        12    n+p-->N*(0)(1535)+p
  !                        13    D(++)+D(-)-->N*(+)(1440)+n
  !                         14    D(++)+D(-)-->N*(0)(1440)+p
  !                        15    D(+)+D(0)--->N*(+)(1440)+n
  !                        16    D(+)+D(0)--->N*(0)(1440)+p
  !                        17    D(++)+D(0)-->N*(+)(1535)+p
  !                        18    D(++)+D(-)-->N*(0)(1535)+p
  !                        19    D(++)+D(-)-->N*(+)(1535)+n
  !                        20    D(+)+D(+)-->N*(+)(1535)+p
  !                        21    D(+)+D(0)-->N*(+)(1535)+n
  !                        22    D(+)+D(0)-->N*(0)(1535)+p
  !                        23    D(+)+D(-)-->N*(0)(1535)+n
  !                        24    D(0)+D(0)-->N*(0)(1535)+n
  !                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p
  !                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n
  !                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n
  !                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p
  !                        29    N*(+)(14)+D+-->N*(+)(15)+p
  !                        30    N*(+)(14)+D0-->N*(+)(15)+n
  !                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n
  !                        32    N*(0)(14)+D++--->N*(+)(15)+p
  !                        33    N*(0)(14)+D+--->N*(+)(15)+n
  !                        34    N*(0)(14)+D+--->N*(0)(15)+p
  !                        35    N*(0)(14)+D0-->N*(0)(15)+n
  !                        36    N*(+)(14)+D0--->N*(0)(15)+p
  !                        ++    see the note book for more listing
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (xmd=1.8756, npdmax=10000)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  !c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  !c      SAVE /gg/
  Common /input/nstar, ndirct, dir
  !c      SAVE /INPUT/
  Common /nn/nnn
  !c      SAVE /NN/
  Common /bg/betax, betay, betaz, gamma
  !c      SAVE /BG/
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
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  !c      SAVE /leadng/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Dimension ppd(3, npdmax), lbpd(npdmax)
  Save
  !-----------------------------------------------------------------------
  n12 = 0
  m12 = 0
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
  c2 = pz/pr
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .And. lb(i2)<0) ianti = 1

  !lin-6/2008 Production of perturbative deuterons for idpert=1:
  Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  If (idpert==1 .And. ipert1==1) Then
     If (srt<2.012) Return
     If ((iabs(lb(i1))==1 .Or. iabs(lb(i1))==2) .And. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) Then
        Goto 108
     Else If ((iabs(lb(i2))==1 .Or. iabs(lb(i2))==2) .And. (iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13)) Then
        Goto 108
     Else
        Return
     End If
  End If
  !-----------------------------------------------------------------------
  !COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R
  !      N-DELTA OR N*-N* or N*-Delta)
  If (x1<=signn/sig) Then
     !COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER
     as = (3.65*(srt-1.8766))**6
     a = 6.0*as/(1.0+as)
     ta = -2.0*pr**2
     x = ranart(nseed)
     !lin-10/24/02        T1  = ALOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A
     t1 = sngl(dlog(dble(1.-x)*dexp(dble(a)*dble(ta))+dble(x)))/a
     c1 = 1.0 - t1/ta
     t1 = 2.0*pi*ranart(nseed)
     iblock = 1
     Goto 107
  Else
     !COM: TEST FOR INELASTIC SCATTERING
     !     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING
     !     CAN HAPPEN ANY MORE ==> RETURN (2.04 = 2*AVMASS + PI-MASS+0.02)
     If (srt<2.04) Return
     !lin-6/2008 add d+meson production for n*N*(0)(1440) and p*N*(+)(1440) channels
     !     (they did not have any inelastic reactions before):
     If (((iabs(lb(i1))==2 .Or. iabs(lb(i2))==2) .And. (lb(i1)*lb(i2))==20) .Or. (lb(i1)*lb(i2))==13) Then
        If (x1<=((signn+sdprod)/sig)) Goto 108
     End If
     !
     ! Resonance absorption or Delta + N-->N*(1440), N*(1535)
     ! COM: TEST FOR DELTA OR N* ABSORPTION
     !      IN THE PROCESS DELTA+N-->NN, N*+N-->NN
     prf = sqrt(0.25*srt**2-avmass**2)
     If (em1>1.) Then
        deltam = em1
     Else
        deltam = em2
     End If
     renom = deltam*prf**2/denom(srt, 1.)/pr
     renomn = deltam*prf**2/denom(srt, 2.)/pr
     renom1 = deltam*prf**2/denom(srt, -1.)/pr
     ! avoid the inelastic collisions between n+delta- -->N+N
     !       and p+delta++ -->N+N due to charge conservation,
     !       but they can scatter to produce kaons
     If ((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) renom = 0.
     If ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) renom = 0.
     If ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) renom = 0.
     If ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9)) renom = 0.
     Call m1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
     x1440 = (3./4.)*sigma(srt, 2, 0, 1)
     ! CROSS SECTION FOR KAON PRODUCTION from the four channels
     ! for NLK channel
     ! avoid the inelastic collisions between n+delta- -->N+N
     !       and p+delta++ -->N+N due to charge conservation,
     !       but they can scatter to produce kaons
     If (((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) .Or. ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) .Or. ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) .Or. ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9))) Then
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !          IF((SIGK+SIGNN)/SIG.GE.X1)GO TO 306
        If ((sigk+signn+sdprod)/sig>=x1) Goto 306
        !
     End If
     ! WE DETERMINE THE REACTION CHANNELS IN THE FOLLOWING
     ! FOR n+delta(++)-->p+p or n+delta(++)-->n+N*(+)(1440),n+N*(+)(1535)
     ! REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN,
     If (lb(i1)*lb(i2)==18 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        ! REABSORPTION:
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 3
           Goto 206
        Else
           ! N* PRODUCTION
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              ! N*(1440)
              m12 = 37
           Else
              ! N*(1535)       M12=38
              !lin-2/26/03 why is the above commented out? leads to M12=0 but
              !     particle mass is changed after 204 (causes energy violation).
              !     replace by elastic process (return):
              Return

           End If
           Goto 204
        End If
     End If
     ! FOR p+delta(-)-->n+n or p+delta(-)-->n+N*(0)(1440),n+N*(0)(1535)
     ! REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN,
     If (lb(i1)*lb(i2)==6 .And. ((iabs(lb(i1))==1) .Or. (iabs(lb(i2))==1))) Then
        signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF (X1.GT.(SIGNN+SIGDN+X1440+X1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        ! REABSORPTION:
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 6
           Goto 206
        Else
           ! N* PRODUCTION
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              ! N*(1440)
              m12 = 47
           Else
              ! N*(1535)       M12=48
              !lin-2/26/03 causes energy violation, replace by elastic process (return):
              Return

           End If
           Goto 204
        End If
     End If
     ! FOR p+delta(+)-->p+p, N*(+)(144)+p, N*(+)(1535)+p
     If (lb(i1)*lb(i2)==8 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
        signd = 1.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 4
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              ! N*(144)
              m12 = 39
           Else
              m12 = 40
           End If
           Goto 204
        End If
     End If
     ! FOR n+delta(0)-->n+n, N*(0)(144)+n, N*(0)(1535)+n
     If (lb(i1)*lb(i2)==14 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = 1.5*sigma(srt, 1, 1, 1)
        sigdn = 0.25*signd*renom
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+x1440+x1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+x1440+x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+x1440+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1440+x1535)) Then
           m12 = 5
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              ! N*(144)
              m12 = 48
           Else
              m12 = 49
           End If
           Goto 204
        End If
     End If
     ! FOR n+delta(+)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
     !                       N*(+)(1535)+n,N*(0)(1535)+p
     If (lb(i1)*lb(i2)==16 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
        sigdn = 0.5*signd*renom
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+2.*x1440+2.*x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+2*x1440+2*x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+2.*x1440+2.*x1535)) Then
           m12 = 1
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 41
              If (ranart(nseed)<=0.5) m12 = 43
           Else
              m12 = 42
              If (ranart(nseed)<=0.5) m12 = 44
           End If
           Goto 204
        End If
     End If
     ! FOR p+delta(0)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
     !                       N*(+)(1535)+n,N*(0)(1535)+p
     If (lb(i1)*lb(i2)==7) Then
        signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
        sigdn = 0.5*signd*renom
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+2.*x1440+2.*x1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+2.*x1440+2.*x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+2*x1440+2*x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+2.*x1440+2.*x1535)) Then
           m12 = 2
           Goto 206
        Else
           If (ranart(nseed)<x1440/(x1440+x1535)) Then
              m12 = 50
              If (ranart(nseed)<=0.5) m12 = 51
           Else
              m12 = 52
              If (ranart(nseed)<=0.5) m12 = 53
           End If
           Goto 204
        End If
     End If
     ! FOR p+N*(0)(14)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
     ! OR  P+N*(0)(14)-->D(+)+N, D(0)+P,
     If (lb(i1)*lb(i2)==10 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
        signd = (3./4.)*sigma(srt, 2, 0, 1)
        sigdn = signd*renomn
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1535)) Then
           m12 = 7
           Goto 206
        Else
           m12 = 54
           If (ranart(nseed)<=0.5) m12 = 55
        End If
        Goto 204
     End If
     ! FOR n+N*(+)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
     If (lb(i1)*lb(i2)==22 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
        signd = (3./4.)*sigma(srt, 2, 0, 1)
        sigdn = signd*renomn
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+X1535+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+x1535+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn+x1535)>ranart(nseed)) Goto 306
        If (ranart(nseed)<sigdn/(sigdn+x1535)) Then
           m12 = 8
           Goto 206
        Else
           m12 = 56
           If (ranart(nseed)<=0.5) m12 = 57
        End If
        Goto 204
     End If
     ! FOR N*(1535)+N-->N+N COLLISIONS
     If ((iabs(lb(i1))==12) .Or. (iabs(lb(i1))==13) .Or. (iabs(lb(i2))==12) .Or. (iabs(lb(i2))==13)) Then
        signd = x1535
        sigdn = signd*renom1
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGDN+SIGK)/SIG)RETURN
        If (x1>(signn+sigdn+sigk+sdprod)/sig) Return
        !
        If (sigk/(sigk+sigdn)>ranart(nseed)) Goto 306
        If (lb(i1)*lb(i2)==24) m12 = 10
        If (lb(i1)*lb(i2)==12) m12 = 12
        If (lb(i1)*lb(i2)==26) m12 = 11
        If (lb(i1)*lb(i2)==13) m12 = 9
        Goto 206
     End If
204  Continue
     ! (1) GENERATE THE MASS FOR THE N*(1440) AND N*(1535)
     ! (2) CALCULATE THE FINAL MOMENTUM OF THE n+N* SYSTEM
     ! (3) RELABLE THE FINAL STATE PARTICLES
     !PARAMETRIZATION OF THE SHAPE OF THE N* RESONANCE ACCORDING
     !     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER
     !     FORMULA FOR N* RESORANCE
     !     DETERMINE DELTA MASS VIA REJECTION METHOD.
     dmax = srt - avmass - 0.005
     dmin = 1.078
     If ((m12==37) .Or. (m12==39) .Or. (m12==41) .Or. (m12==43) .Or. (m12==46) .Or. (m12==48) .Or. (m12==50) .Or. (m12==51)) Then
        ! N*(1440) production
        If (dmax<1.44) Then
           fm = fns(dmax, srt, 0.)
        Else

           !lin-10/25/02 get rid of argument usage mismatch in FNS():
           xdmass = 1.44
           !          FM=FNS(1.44,SRT,1.)
           fm = fns(xdmass, srt, 1.)
           !lin-10/25/02-end

        End If
        If (fm==0.) fm = 1.E-09
        ntry2 = 0
11      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry2 = ntry2 + 1
        If ((ranart(nseed)>fns(dm,srt,1.)/fm) .And. (ntry2<=10)) Goto 11

        !lin-2/26/03 limit the N* mass below a certain value
        !     (here taken as its central value + 2* B-W fullwidth):
        If (dm>2.14) Goto 11

        Goto 13
     Else
        ! N*(1535) production
        If (dmax<1.535) Then
           fm = fd5(dmax, srt, 0.)
        Else

           !lin-10/25/02 get rid of argument usage mismatch in FNS():
           xdmass = 1.535
           !          FM=FD5(1.535,SRT,1.)
           fm = fd5(xdmass, srt, 1.)
           !lin-10/25/02-end

        End If
        If (fm==0.) fm = 1.E-09
        ntry1 = 0
12      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry1 = ntry1 + 1
        If ((ranart(nseed)>fd5(dm,srt,1.)/fm) .And. (ntry1<=10)) Goto 12

        !lin-2/26/03 limit the N* mass below a certain value
        !     (here taken as its central value + 2* B-W fullwidth):
        If (dm>1.84) Goto 12

     End If
13   Continue
     ! (2) DETERMINE THE FINAL MOMENTUM
     prf = 0.
     pf2 = ((srt**2-dm**2+avmass**2)/(2.*srt))**2 - avmass**2
     If (pf2>0.) prf = sqrt(pf2)
     ! (3) RELABLE FINAL STATE PARTICLES
     ! 37 D(++)+n-->N*(+)(14)+p
     If (m12==37) Then
        If (iabs(lb(i1))==9) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 38 D(++)+n-->N*(+)(15)+p
     If (m12==38) Then
        If (iabs(lb(i1))==9) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 39 D(+)+P-->N*(+)(14)+p
     If (m12==39) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 40 D(+)+P-->N*(+)(15)+p
     If (m12==40) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 41 D(+)+N-->N*(+)(14)+N
     If (m12==41) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 42 D(+)+N-->N*(+)(15)+N
     If (m12==42) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 43 D(+)+N-->N*(0)(14)+P
     If (m12==43) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 44 D(+)+N-->N*(0)(15)+P
     If (m12==44) Then
        If (iabs(lb(i1))==8) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 46 D(-)+P-->N*(0)(14)+N
     If (m12==46) Then
        If (iabs(lb(i1))==6) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 47 D(-)+P-->N*(0)(15)+N
     If (m12==47) Then
        If (iabs(lb(i1))==6) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 48 D(0)+N-->N*(0)(14)+N
     If (m12==48) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 49 D(0)+N-->N*(0)(15)+N
     If (m12==49) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 50 D(0)+P-->N*(0)(14)+P
     If (m12==50) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 10
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 10
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 51 D(0)+P-->N*(+)(14)+N
     If (m12==51) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 11
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 11
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 52 D(0)+P-->N*(0)(15)+P
     If (m12==52) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 53 D(0)+P-->N*(+)(15)+N
     If (m12==53) Then
        If (iabs(lb(i1))==7) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 54 N*(0)(14)+P-->N*(+)(15)+N
     If (m12==54) Then
        If (iabs(lb(i1))==10) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 55 N*(0)(14)+P-->N*(0)(15)+P
     If (m12==55) Then
        If (iabs(lb(i1))==10) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 56 N*(+)(14)+N-->N*(+)(15)+N
     If (m12==56) Then
        If (iabs(lb(i1))==11) Then
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 13
           e(i2) = dm
        Else
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 13
           e(i1) = dm
        End If
        Goto 207
     End If
     ! 57 N*(+)(14)+N-->N*(0)(15)+P
     If (m12==57) Then
        If (iabs(lb(i1))==11) Then
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 12
           e(i2) = dm
        Else
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 12
           e(i1) = dm
        End If
     End If
     Goto 207
     !------------------------------------------------
     ! RELABLE NUCLEONS AFTER DELTA OR N* BEING ABSORBED
     !(1) n+delta(+)-->n+p
206  If (m12==1) Then
        If (iabs(lb(i1))==8) Then
           lb(i2) = 2
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 2
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 207
     End If
     !(2) p+delta(0)-->p+n
     If (m12==2) Then
        If (iabs(lb(i1))==7) Then
           lb(i2) = 1
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 207
     End If
     !(3) n+delta(++)-->p+p
     If (m12==3) Then
        lb(i1) = 1
        lb(i2) = 1
        e(i1) = amp
        e(i2) = amp
        Goto 207
     End If
     !(4) p+delta(+)-->p+p
     If (m12==4) Then
        lb(i1) = 1
        lb(i2) = 1
        e(i1) = amp
        e(i2) = amp
        Goto 207
     End If
     !(5) n+delta(0)-->n+n
     If (m12==5) Then
        lb(i1) = 2
        lb(i2) = 2
        e(i1) = amn
        e(i2) = amn
        Goto 207
     End If
     !(6) p+delta(-)-->n+n
     If (m12==6) Then
        lb(i1) = 2
        lb(i2) = 2
        e(i1) = amn
        e(i2) = amn
        Goto 207
     End If
     !(7) p+N*(0)-->n+p
     If (m12==7) Then
        If (iabs(lb(i1))==1) Then
           lb(i1) = 1
           lb(i2) = 2
           e(i1) = amp
           e(i2) = amn
        Else
           lb(i1) = 2
           lb(i2) = 1
           e(i1) = amn
           e(i2) = amp
        End If
        Goto 207
     End If
     !(8) n+N*(+)-->n+p
     If (m12==8) Then
        If (iabs(lb(i1))==2) Then
           lb(i1) = 2
           lb(i2) = 1
           e(i1) = amn
           e(i2) = amp
        Else
           lb(i1) = 1
           lb(i2) = 2
           e(i1) = amp
           e(i2) = amn
        End If
        Goto 207
     End If
     !lin-6/2008
     !*(9) N*(+)p-->pp
     !(9) N*(+)(1535) p-->pp
     If (m12==9) Then
        lb(i1) = 1
        lb(i2) = 1
        e(i1) = amp
        e(i2) = amp
        Goto 207
     End If
     !(12) N*(0)P-->nP
     If (m12==12) Then
        lb(i1) = 2
        lb(i2) = 1
        e(i1) = amn
        e(i2) = amp
        Goto 207
     End If
     !(11) N*(+)n-->nP
     If (m12==11) Then
        lb(i1) = 2
        lb(i2) = 1
        e(i1) = amn
        e(i2) = amp
        Goto 207
     End If
     !lin-6/2008
     !*(12) N*(0)p-->Np
     !(12) N*(0)(1535) p-->Np
     If (m12==12) Then
        lb(i1) = 1
        lb(i2) = 2
        e(i1) = amp
        e(i2) = amn
     End If
     !----------------------------------------------
207  pr = prf
     c1 = 1.0 - 2.0*ranart(nseed)
     If (srt<=2.14) c1 = 1.0 - 2.0*ranart(nseed)
     If (srt>2.14 .And. srt<=2.4) c1 = ang(srt, iseed)
     If (srt>2.4) Then

        !lin-10/25/02 get rid of argument usage mismatch in PTR():
        xptr = 0.33*pr
        !         cc1=ptr(0.33*pr,iseed)
        cc1 = ptr(xptr, iseed)
        !lin-10/25/02-end

        !lin-9/2012: check argument in sqrt():
        scheck = pr**2 - cc1**2
        If (scheck<0) Then
           Write (99, *) 'scheck4: ', scheck
           scheck = 0.
        End If
        c1 = sqrt(scheck)/pr
        !         c1=sqrt(pr**2-cc1**2)/pr

     End If
     t1 = 2.0*pi*ranart(nseed)
     iblock = 3
  End If
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
  End If

  !-----------------------------------------------------------------------
  !COM: SET THE NEW MOMENTUM COORDINATES
107 If (px==0.0 .And. py==0.0) Then
     t2 = 0.0
  Else
     t2 = atan2(py, px)
  End If

  !lin-9/2012: check argument in sqrt():
  scheck = 1.0 - c1**2
  If (scheck<0) Then
     Write (99, *) 'scheck5: ', scheck
     scheck = 0.
  End If
  s1 = sqrt(scheck)
  !      S1   = SQRT( 1.0 - C1**2 )

  !lin-9/2012: check argument in sqrt():
  scheck = 1.0 - c2**2
  If (scheck<0) Then
     Write (99, *) 'scheck6: ', scheck
     scheck = 0.
  End If
  s2 = sqrt(scheck)
  !      S2  =  SQRT( 1.0 - C2**2 )

  ct1 = cos(t1)
  st1 = sin(t1)
  ct2 = cos(t2)
  st2 = sin(t2)
  pz = pr*(c1*c2-s1*s2*ct1)
  ss = c2*s1*ct1 + s2*c1
  px = pr*(ss*ct2-s1*st1*st2)
  py = pr*(ss*st2+s1*st1*ct2)
  Return
  ! FOR THE NN-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN
  ! THE NUCLEUS-NUCLEUS CMS.
306 Continue
  !sp11/21/01 phi production
  If (xsk5/sigk>ranart(nseed)) Then
     pz1 = p(3, i1)
     pz2 = p(3, i2)
     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 1 + int(2*ranart(nseed))
     nnn = nnn + 1
     lpion(nnn, irun) = 29
     epion(nnn, irun) = aphi
     iblock = 222
     Goto 208
  End If
  !sp11/21/01 end
  iblock = 11
  If (ianti==1) iblock = -11
  !
  pz1 = p(3, i1)
  pz2 = p(3, i2)
  ! DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  nnn = nnn + 1
  lpion(nnn, irun) = 23
  epion(nnn, irun) = aka
  If (srt<=2.63) Then
     ! only lambda production is possible
     ! (1.1)P+P-->p+L+kaon+
     ic = 1

     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 14
     Goto 208
  End If
  If (srt<=2.74 .And. srt>2.63) Then
     ! both Lambda and sigma production are possible
     If (xsk1/(xsk1+xsk2)>ranart(nseed)) Then
        ! lambda production
        ic = 1

        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
     Else
        ! sigma production

        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 15 + int(3*ranart(nseed))
        ic = 2
     End If
     Goto 208
  End If
  If (srt<=2.77 .And. srt>2.74) Then
     ! then pp-->Delta lamda kaon can happen
     If (xsk1/(xsk1+xsk2+xsk3)>ranart(nseed)) Then
        ! * (1.1)P+P-->p+L+kaon+
        ic = 1

        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk2/(xsk2+xsk3)>ranart(nseed)) Then
           ! pp-->psk
           ic = 2

           lb(i1) = 1 + int(2*ranart(nseed))
           lb(i2) = 15 + int(3*ranart(nseed))

        Else
           ! pp-->D+l+k
           ic = 3

           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
        End If
        Goto 208
     End If
  End If
  If (srt>2.77) Then
     ! all four channels are possible
     If (xsk1/(xsk1+xsk2+xsk3+xsk4)>ranart(nseed)) Then
        ! p lambda k production
        ic = 1

        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk3/(xsk2+xsk3+xsk4)>ranart(nseed)) Then
           ! delta l K production
           ic = 3

           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
           Goto 208
        Else
           If (xsk2/(xsk2+xsk4)>ranart(nseed)) Then
              ! n sigma k production

              lb(i1) = 1 + int(2*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))

              ic = 2
           Else
              ic = 4

              lb(i1) = 6 + int(4*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))

           End If
           Goto 208
        End If
     End If
  End If
208 Continue
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==23) lpion(nnn, irun) = 21
  End If
  lbi1 = lb(i1)
  lbi2 = lb(i2)
  ! KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE
  ntry1 = 0
128 Call bbkaon(ic, srt, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 128
  !       if(icou1.lt.0)return
  ! ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  ! FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
  ! NUCLEUS CMS. FRAME
  ! (1) for the necleon/delta
  !             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
  e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  ! (2) for the lambda/sigma
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  ! GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
  epcm = sqrt(aka**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  !lin-5/2008:
  !2008        X01 = 1.0 - 2.0 * RANART(NSEED)
  !            Y01 = 1.0 - 2.0 * RANART(NSEED)
  !            Z01 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2008
  !                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
  !                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
  !                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  !
  ! assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the
  ! leadng particle behaviour
  !              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lbi1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lbi2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  If (lpion(nnn,irun)/=29) iblock = 11
  lb1 = lb(i1)
  lb2 = lb(i2)
  am1 = em1
  am2 = em2
  e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
  Return

  !lin-6/2008 N+D->Deuteron+pi:
  !     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
108 Continue
  If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
     !     For idpert=1: we produce npertd pert deuterons:
     ndloop = npertd
  Else If (idpert==2 .And. npertd>=1) Then
     !     For idpert=2: we first save information for npertd pert deuterons;
     !     at the last ndloop we create the regular deuteron+pi
     !     and those pert deuterons:
     ndloop = npertd + 1
  Else
     !     Just create the regular deuteron+pi:
     ndloop = 1
  End If
  !
  dprob1 = sdprod/sig/float(npertd)
  Do idloop = 1, ndloop
     Call bbdangle(pxd, pyd, pzd, nt, ipert1, ianti, idloop, pfinal, dprob1, lbm)
     Call rotate(px, py, pz, pxd, pyd, pzd)
     !     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE
     !     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME:
     !     For the Deuteron:
     xmass = xmd
     e1dcm = sqrt(xmass**2+pxd**2+pyd**2+pzd**2)
     p1dbeta = pxd*betax + pyd*betay + pzd*betaz
     transf = gamma*(gamma*p1dbeta/(gamma+1.)+e1dcm)
     pxi1 = betax*transf + pxd
     pyi1 = betay*transf + pyd
     pzi1 = betaz*transf + pzd
     If (ianti==0) Then
        lbd = 42
     Else
        lbd = -42
     End If
     If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
        !ccc  Perturbative production for idpert=1:
        nnn = nnn + 1
        ppion(1, nnn, irun) = pxi1
        ppion(2, nnn, irun) = pyi1
        ppion(3, nnn, irun) = pzi1
        epion(nnn, irun) = xmd
        lpion(nnn, irun) = lbd
        rpion(1, nnn, irun) = r(1, i1)
        rpion(2, nnn, irun) = r(2, i1)
        rpion(3, nnn, irun) = r(3, i1)
        !lin-6/2008 assign the perturbative probability:
        dppion(nnn, irun) = sdprod/sig/float(npertd)
     Else If (idpert==2 .And. idloop<=npertd) Then
        !lin-6/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons
        !     only when a regular (anti)deuteron+pi is produced in NN collisions.
        !     First save the info for the perturbative deuterons:
        ppd(1, idloop) = pxi1
        ppd(2, idloop) = pyi1
        ppd(3, idloop) = pzi1
        lbpd(idloop) = lbd
     Else
        !ccc  Regular production:
        !     For the regular pion: do LORENTZ-TRANSFORMATION:
        e(i1) = xmm
        e2picm = sqrt(xmm**2+pxd**2+pyd**2+pzd**2)
        p2pibeta = -pxd*betax - pyd*betay - pzd*betaz
        transf = gamma*(gamma*p2pibeta/(gamma+1.)+e2picm)
        pxi2 = betax*transf - pxd
        pyi2 = betay*transf - pyd
        pzi2 = betaz*transf - pzd
        p(1, i1) = pxi2
        p(2, i1) = pyi2
        p(3, i1) = pzi2
        !     Remove regular pion to check the equivalence
        !     between the perturbative and regular deuteron results:
        !                 E(i1)=0.
        !
        lb(i1) = lbm
        px1 = p(1, i1)
        py1 = p(2, i1)
        pz1 = p(3, i1)
        em1 = e(i1)
        id(i1) = 2
        id1 = id(i1)
        e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
        lb1 = lb(i1)
        !     For the regular deuteron:
        p(1, i2) = pxi1
        p(2, i2) = pyi1
        p(3, i2) = pzi1
        lb(i2) = lbd
        lb2 = lb(i2)
        e(i2) = xmd
        eti2 = e(i2)
        id(i2) = 2
        !     For idpert=2: create the perturbative deuterons:
        If (idpert==2 .And. idloop==ndloop) Then
           Do ipertd = 1, npertd
              nnn = nnn + 1
              ppion(1, nnn, irun) = ppd(1, ipertd)
              ppion(2, nnn, irun) = ppd(2, ipertd)
              ppion(3, nnn, irun) = ppd(3, ipertd)
              epion(nnn, irun) = xmd
              lpion(nnn, irun) = lbpd(ipertd)
              rpion(1, nnn, irun) = r(1, i1)
              rpion(2, nnn, irun) = r(2, i1)
              rpion(3, nnn, irun) = r(3, i1)
              !lin-6/2008 assign the perturbative probability:
              dppion(nnn, irun) = 1./float(npertd)
           End Do
        End If
     End If
  End Do
  iblock = 501
  Return
  !lin-6/2008 N+D->Deuteron+pi over

End Subroutine crnd


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!*********************************
!                                                                      *
!                                                                      *
Subroutine crdd(irun, px, py, pz, srt, i1, i2, iblock, ntag, signn, sig, nt, ipert1)
  !     1NTAG,SIGNN,SIG)
  !     PURPOSE:                                                         *
  !             DEALING WITH BARYON RESONANCE-BARYON RESONANCE COLLISIONS*
  !     NOTE   :                                                         *
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   *
  !           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                      0-> COLLISION CANNOT HAPPEN                     *
  !                      1-> N-N ELASTIC COLLISION                       *
  !                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          *
  !                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          *
  !                      4-> N+N->N+N+PION,DIRTCT PROCESS                *
  !                     5-> DELTA(N*)+DELTA(N*)   TOTAL   COLLISIONS    *
  !           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      *
  !                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    *
  !                      N12,                                            *
  !                      M12=1 FOR p+n-->delta(+)+ n                     *
  !                          2     p+n-->delta(0)+ p                     *
  !                          3     p+p-->delta(++)+n                     *
  !                          4     p+p-->delta(+)+p                      *
  !                          5     n+n-->delta(0)+n                      *
  !                          6     n+n-->delta(-)+p                      *
  !                          7     n+p-->N*(0)(1440)+p                   *
  !                          8     n+p-->N*(+)(1440)+n                   *
  !                        9     p+p-->N*(+)(1535)+p                     *
  !                        10    n+n-->N*(0)(1535)+n                     *
  !                         11    n+p-->N*(+)(1535)+n                     *
  !                        12    n+p-->N*(0)(1535)+p
  !                        13    D(++)+D(-)-->N*(+)(1440)+n
  !                         14    D(++)+D(-)-->N*(0)(1440)+p
  !                        15    D(+)+D(0)--->N*(+)(1440)+n
  !                        16    D(+)+D(0)--->N*(0)(1440)+p
  !                        17    D(++)+D(0)-->N*(+)(1535)+p
  !                        18    D(++)+D(-)-->N*(0)(1535)+p
  !                        19    D(++)+D(-)-->N*(+)(1535)+n
  !                        20    D(+)+D(+)-->N*(+)(1535)+p
  !                        21    D(+)+D(0)-->N*(+)(1535)+n
  !                        22    D(+)+D(0)-->N*(0)(1535)+p
  !                        23    D(+)+D(-)-->N*(0)(1535)+n
  !                        24    D(0)+D(0)-->N*(0)(1535)+n
  !                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p
  !                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n
  !                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n
  !                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p
  !                        29    N*(+)(14)+D+-->N*(+)(15)+p
  !                        30    N*(+)(14)+D0-->N*(+)(15)+n
  !                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n
  !                        32    N*(0)(14)+D++--->N*(+)(15)+p
  !                        33    N*(0)(14)+D+--->N*(+)(15)+n
  !                        34    N*(0)(14)+D+--->N*(0)(15)+p
  !                        35    N*(0)(14)+D0-->N*(0)(15)+n
  !                        36    N*(+)(14)+D0--->N*(0)(15)+p
  !                        +++
  !               AND MORE CHANNELS AS LISTED IN THE NOTE BOOK
  !
  ! NOTE ABOUT N*(1440) RESORANCE:                                       *
  !     As it has been discussed in VerWest's paper,I= 1 (initial isospin)
  !     channel can all be attributed to delta resorance while I= 0      *
  !     channel can all be  attribured to N* resorance.Only in n+p       *
  !     one can have I=0 channel so is the N*(1440) resorance            *
  ! REFERENCES:    J. CUGNON ET AL., NUCL. PHYS. A352, 505 (1981)        *
  !                    Y. KITAZOE ET AL., PHYS. LETT. 166B, 35 (1986)    *
  !                    B. VerWest el al., PHYS. PRV. C25 (1982)1979      *
  !                    Gy. Wolf  et al, Nucl Phys A517 (1990) 615        *
  !                    CUTOFF = 2 * AVMASS + 20 MEV                      *
  !                                                                      *
  !       for N*(1535) we use the parameterization by Gy. Wolf et al     *
  !       Nucl phys A552 (1993) 349, added May 18, 1994                  *
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (xmd=1.8756, npdmax=10000)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
  !c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
  !c      SAVE /gg/
  Common /input/nstar, ndirct, dir
  !c      SAVE /INPUT/
  Common /nn/nnn
  !c      SAVE /NN/
  Common /bg/betax, betay, betaz, gamma
  !c      SAVE /BG/
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
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  !c      SAVE /leadng/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Dimension ppd(3, npdmax), lbpd(npdmax)
  Save
  !-----------------------------------------------------------------------
  n12 = 0
  m12 = 0
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
  c2 = pz/pr
  If (px==0.0 .And. py==0.0) Then
     t2 = 0.0
  Else
     t2 = atan2(py, px)
  End If
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .And. lb(i2)<0) ianti = 1

  !lin-6/2008 Production of perturbative deuterons for idpert=1:
  Call sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  If (idpert==1 .And. ipert1==1) Then
     If (srt<2.012) Return
     If ((iabs(lb(i1))>=6 .And. iabs(lb(i1))<=13) .And. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=13)) Then
        Goto 108
     Else
        Return
     End If
  End If

  !-----------------------------------------------------------------------
  !COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R
  !      N-DELTA OR N*-N* or N*-Delta)
  If (x1<=signn/sig) Then
     !COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER
     as = (3.65*(srt-1.8766))**6
     a = 6.0*as/(1.0+as)
     ta = -2.0*pr**2
     x = ranart(nseed)
     !lin-10/24/02        T1  = DLOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A
     t1 = sngl(dlog(dble(1.-x)*dexp(dble(a)*dble(ta))+dble(x)))/a
     c1 = 1.0 - t1/ta
     t1 = 2.0*pi*ranart(nseed)
     iblock = 20
     Goto 107
  Else
     !COM: TEST FOR INELASTIC SCATTERING
     !     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING
     !     CAN HAPPEN ANY MORE ==> RETURN (2.15 = 2*AVMASS +2*PI-MASS)
     If (srt<2.15) Return
     !     IF THERE WERE 2 N*(1535) AND THEY DIDN'T SCATT. ELAST.,
     !     ALLOW THEM TO PRODUCE KAONS. NO OTHER INELASTIC CHANNELS
     !     ARE KNOWN
     !       if((lb(i1).ge.12).and.(lb(i2).ge.12))return
     !     ALL the inelastic collisions between N*(1535) and Delta as well
     !     as N*(1440) TO PRODUCE KAONS, NO OTHER CHANNELS ARE KNOWN
     !       if((lb(i1).ge.12).and.(lb(i2).ge.3))return
     !       if((lb(i2).ge.12).and.(lb(i1).ge.3))return
     !     calculate the N*(1535) production cross section in I1+I2 collisions
     Call n1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)

     ! for Delta+Delta-->N*(1440 OR 1535)+N AND N*(1440)+N*(1440)-->N*(1535)+X
     !     AND DELTA+N*(1440)-->N*(1535)+X
     ! WE ASSUME THEY HAVE THE SAME CROSS SECTIONS as CORRESPONDING N+N COLLISION):
     ! FOR D++D0, D+D+,D+D-,D0D0,N*+N*+,N*0N*0,N*(+)D+,N*(+)D(-),N*(0)D(0)
     ! N*(1535) production, kaon production and reabsorption through
     ! D(N*)+D(N*)-->NN are ALLOWED.
     ! CROSS SECTION FOR KAON PRODUCTION from the four channels are
     ! for NLK channel
     akp = 0.498
     ak0 = 0.498
     ana = 0.938
     ada = 1.232
     al = 1.1157
     as = 1.1197
     xsk1 = 0
     xsk2 = 0
     xsk3 = 0
     xsk4 = 0
     xsk5 = 0
     t1nlk = ana + al + akp
     If (srt<=t1nlk) Goto 222
     xsk1 = 1.5*pplpk(srt)
     ! for DLK channel
     t1dlk = ada + al + akp
     t2dlk = ada + al - akp
     If (srt<=t1dlk) Goto 222
     es = srt
     pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
     pmdlk = sqrt(pmdlk2)
     xsk3 = 1.5*pplpk(srt)
     ! for NSK channel
     t1nsk = ana + as + akp
     t2nsk = ana + as - akp
     If (srt<=t1nsk) Goto 222
     pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
     pmnsk = sqrt(pmnsk2)
     xsk2 = 1.5*(ppk1(srt)+ppk0(srt))
     ! for DSK channel
     t1dsk = ada + as + akp
     t2dsk = ada + as - akp
     If (srt<=t1dsk) Goto 222
     pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
     pmdsk = sqrt(pmdsk2)
     xsk4 = 1.5*(ppk1(srt)+ppk0(srt))
     !sp11/21/01
     ! phi production
     If (srt<=(2.*amn+aphi)) Goto 222
     !  !! mb put the correct form
     xsk5 = 0.0001
     !sp11/21/01 end
     ! THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN
222  sigk = xsk1 + xsk2 + xsk3 + xsk4

     !bz3/7/99 neutralk
     xsk1 = 2.0*xsk1
     xsk2 = 2.0*xsk2
     xsk3 = 2.0*xsk3
     xsk4 = 2.0*xsk4
     sigk = 2.0*sigk + xsk5
     !bz3/7/99 neutralk end

     ! The reabsorption cross section for the process
     ! D(N*)D(N*)-->NN is
     s2d = reab2d(i1, i2, srt)

     !bz3/16/99 pion
     s2d = 0.
     !bz3/16/99 pion end

     !(1) N*(1535)+D(N*(1440)) reactions
     !    we allow kaon production and reabsorption only
     If (((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=12)) .Or. ((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=6)) .Or. ((iabs(lb(i2))>=12) .And. (iabs(lb(i1))>=6))) Then
        signd = sigk + s2d
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !       if(x1.gt.(signd+signn)/sig)return
        If (x1>(signd+signn+sdprod)/sig) Return
        !
        ! if kaon production
        !lin-6/2008
        !       IF(SIGK/SIG.GE.RANART(NSEED))GO TO 306
        If ((sigk+sdprod)/sig>=ranart(nseed)) Goto 306
        !
        ! if reabsorption
        Goto 1012
     End If
     idd = iabs(lb(i1)*lb(i2))
     ! channels have the same charge as pp
     If ((idd==63) .Or. (idd==64) .Or. (idd==48) .Or. (idd==49) .Or. (idd==11*11) .Or. (idd==10*10) .Or. (idd==88) .Or. (idd==66) .Or. (idd==90) .Or. (idd==70)) Then
        signd = x1535 + sigk + s2d
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
        If (x1>(signn+signd+sdprod)/sig) Return
        !
        ! if kaon production
        If (sigk/signd>ranart(nseed)) Goto 306
        ! if reabsorption
        If (s2d/(x1535+s2d)>ranart(nseed)) Goto 1012
        ! if N*(1535) production
        If (idd==63) n12 = 17
        If (idd==64) n12 = 20
        If (idd==48) n12 = 23
        If (idd==49) n12 = 24
        If (idd==121) n12 = 25
        If (idd==100) n12 = 26
        If (idd==88) n12 = 29
        If (idd==66) n12 = 31
        If (idd==90) n12 = 32
        If (idd==70) n12 = 35
        Goto 1011
     End If
     ! IN DELTA+N*(1440) and N*(1440)+N*(1440) COLLISIONS,
     ! N*(1535), kaon production and reabsorption are ALLOWED
     ! IN N*(1440)+N*(1440) COLLISIONS, ONLY N*(1535) IS ALLOWED
     If ((idd==110) .Or. (idd==77) .Or. (idd==80)) Then
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !       IF(X1.GT.(SIGNN+X1535+SIGK+s2d)/SIG)RETURN
        If (x1>(signn+x1535+sigk+s2d+sdprod)/sig) Return
        !
        If (sigk/(x1535+sigk+s2d)>ranart(nseed)) Goto 306
        If (s2d/(x1535+s2d)>ranart(nseed)) Goto 1012
        If (idd==77) n12 = 30
        If ((idd==77) .And. (ranart(nseed)<=0.5)) n12 = 36
        If (idd==80) n12 = 34
        If ((idd==80) .And. (ranart(nseed)<=0.5)) n12 = 35
        If (idd==110) n12 = 27
        If ((idd==110) .And. (ranart(nseed)<=0.5)) n12 = 28
        Goto 1011
     End If
     If ((idd==54) .Or. (idd==56)) Then
        ! LIKE FOR N+P COLLISION,
        ! IN DELTA+DELTA COLLISIONS BOTH N*(1440) AND N*(1535) CAN BE PRODUCED
        sig2 = (3./4.)*sigma(srt, 2, 0, 1)
        signd = 2.*(sig2+x1535) + sigk + s2d
        !lin-6/2008
        If (x1<=((signn+sdprod)/sig)) Goto 108
        !        IF(X1.GT.(SIGNN+SIGND)/SIG)RETURN
        If (x1>(signn+signd+sdprod)/sig) Return
        !
        If (sigk/signd>ranart(nseed)) Goto 306
        If (s2d/(2.*(sig2+x1535)+s2d)>ranart(nseed)) Goto 1012
        If (ranart(nseed)<x1535/(sig2+x1535)) Then
           ! N*(1535) PRODUCTION
           If (idd==54) n12 = 18
           If ((idd==54) .And. (ranart(nseed)<=0.5)) n12 = 19
           If (idd==56) n12 = 21
           If ((idd==56) .And. (ranart(nseed)<=0.5)) n12 = 22
        Else
           ! N*(144) PRODUCTION
           If (idd==54) n12 = 13
           If ((idd==54) .And. (ranart(nseed)<=0.5)) n12 = 14
           If (idd==56) n12 = 15
           If ((idd==56) .And. (ranart(nseed)<=0.5)) n12 = 16
        End If
     End If
1011 Continue
     iblock = 5
     !PARAMETRIZATION OF THE SHAPE OF THE N*(1440) AND N*(1535)
     ! RESONANCE ACCORDING
     !     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER
     !     FORMULA FOR N* RESORANCE
     !     DETERMINE DELTA MASS VIA REJECTION METHOD.
     dmax = srt - avmass - 0.005
     dmin = 1.078
     If ((n12>=13) .And. (n12<=16)) Then
        ! N*(1440) production
        If (dmax<1.44) Then
           fm = fns(dmax, srt, 0.)
        Else

           !lin-10/25/02 get rid of argument usage mismatch in FNS():
           xdmass = 1.44
           !          FM=FNS(1.44,SRT,1.)
           fm = fns(xdmass, srt, 1.)
           !lin-10/25/02-end

        End If
        If (fm==0.) fm = 1.E-09
        ntry2 = 0
11      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry2 = ntry2 + 1
        If ((ranart(nseed)>fns(dm,srt,1.)/fm) .And. (ntry2<=10)) Goto 11

        !lin-2/26/03 limit the N* mass below a certain value
        !     (here taken as its central value + 2* B-W fullwidth):
        If (dm>2.14) Goto 11

        Goto 13
     End If
     If ((n12>=17) .And. (n12<=36)) Then
        ! N*(1535) production
        If (dmax<1.535) Then
           fm = fd5(dmax, srt, 0.)
        Else

           !lin-10/25/02 get rid of argument usage mismatch in FNS():
           xdmass = 1.535
           !          FM=FD5(1.535,SRT,1.)
           fm = fd5(xdmass, srt, 1.)
           !lin-10/25/02-end

        End If
        If (fm==0.) fm = 1.E-09
        ntry1 = 0
12      dm = ranart(nseed)*(dmax-dmin) + dmin
        ntry1 = ntry1 + 1
        If ((ranart(nseed)>fd5(dm,srt,1.)/fm) .And. (ntry1<=10)) Goto 12

        !lin-2/26/03 limit the N* mass below a certain value
        !     (here taken as its central value + 2* B-W fullwidth):
        If (dm>1.84) Goto 12

     End If
13   Continue
     !-------------------------------------------------------
     ! RELABLE BARYON I1 AND I2
     !13 D(++)+D(-)--> N*(+)(14)+n
     If (n12==13) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 11
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 11
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !14 D(++)+D(-)--> N*(0)(14)+P
     If (n12==14) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 10
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 10
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !15 D(+)+D(0)--> N*(+)(14)+n
     If (n12==15) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 11
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 11
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !16 D(+)+D(0)--> N*(0)(14)+P
     If (n12==16) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 10
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 10
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !17 D(++)+D(0)--> N*(+)(14)+P
     If (n12==17) Then
        lb(i2) = 13
        e(i2) = dm
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     !18 D(++)+D(-)--> N*(0)(15)+P
     If (n12==18) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !19 D(++)+D(-)--> N*(+)(15)+N
     If (n12==19) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !20 D(+)+D(+)--> N*(+)(15)+P
     If (n12==20) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !21 D(+)+D(0)--> N*(+)(15)+N
     If (n12==21) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !22 D(+)+D(0)--> N*(0)(15)+P
     If (n12==22) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !23 D(+)+D(-)--> N*(0)(15)+N
     If (n12==23) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !24 D(0)+D(0)--> N*(0)(15)+N
     If (n12==24) Then
        lb(i2) = 12
        e(i2) = dm
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     !25 N*(+)+N*(+)--> N*(0)(15)+P
     If (n12==25) Then
        lb(i2) = 12
        e(i2) = dm
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     !26 N*(0)+N*(0)--> N*(0)(15)+N
     If (n12==26) Then
        lb(i2) = 12
        e(i2) = dm
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     !27 N*(+)+N*(0)--> N*(+)(15)+N
     If (n12==27) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !28 N*(+)+N*(0)--> N*(0)(15)+P
     If (n12==28) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !27 N*(+)+N*(0)--> N*(+)(15)+N
     If (n12==27) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !29 N*(+)+D(+)--> N*(+)(15)+P
     If (n12==29) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !30 N*(+)+D(0)--> N*(+)(15)+N
     If (n12==30) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !31 N*(+)+D(-)--> N*(0)(15)+N
     If (n12==31) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !32 N*(0)+D(++)--> N*(+)(15)+P
     If (n12==32) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !33 N*(0)+D(+)--> N*(+)(15)+N
     If (n12==33) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 13
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 13
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !34 N*(0)+D(+)--> N*(0)(15)+P
     If (n12==34) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     !35 N*(0)+D(0)--> N*(0)(15)+N
     If (n12==35) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !36 N*(+)+D(0)--> N*(0)(15)+P
     If (n12==36) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 12
           e(i2) = dm
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 12
           e(i1) = dm
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
1012 Continue
     iblock = 55
     lb1 = lb(i1)
     lb2 = lb(i2)
     ich = iabs(lb1*lb2)
     !-------------------------------------------------------
     ! RELABLE BARYON I1 AND I2 in the reabsorption processes
     !37 D(++)+D(-)--> n+p
     If (ich==9*6) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !38 D(+)+D(0)--> n+p
     If (ich==8*7) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !39 D(++)+D(0)--> p+p
     If (ich==9*7) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     !40 D(+)+D(+)--> p+p
     If (ich==8*8) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     !41 D(+)+D(-)--> n+n
     If (ich==8*6) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     !42 D(0)+D(0)--> n+n
     If (ich==6*6) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     !43 N*(+)+N*(+)--> p+p
     If (ich==11*11 .Or. ich==13*13 .Or. ich==11*13) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     !44 N*(0)(1440)+N*(0)--> n+n
     If (ich==10*10 .Or. ich==12*12 .Or. ich==10*12) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     !45 N*(+)+N*(0)--> n+p
     If (ich==10*11 .Or. ich==12*13 .Or. ich==10*13 .Or. ich==11*12) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !46 N*(+)+D(+)--> p+p
     If (ich==11*8 .Or. ich==13*8) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     !47 N*(+)+D(0)--> n+p
     If (ich==11*7 .Or. ich==13*7) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 1
           e(i2) = amp
           lb(i1) = 2
           e(i1) = amn
        Else
           lb(i1) = 1
           e(i1) = amp
           lb(i2) = 2
           e(i2) = amn
        End If
        Goto 200
     End If
     !48 N*(+)+D(-)--> n+n
     If (ich==11*6 .Or. ich==13*6) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     !49 N*(0)+D(++)--> p+p
     If (ich==10*9 .Or. ich==12*9) Then
        lb(i2) = 1
        e(i2) = amp
        lb(i1) = 1
        e(i1) = amp
        Goto 200
     End If
     !50 N*(0)+D(0)--> n+n
     If (ich==10*7 .Or. ich==12*7) Then
        lb(i2) = 2
        e(i2) = amn
        lb(i1) = 2
        e(i1) = amn
        Goto 200
     End If
     !51 N*(0)+D(+)--> n+p
     If (ich==10*8 .Or. ich==12*8) Then
        If (ranart(nseed)<=0.5) Then
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 1
           e(i1) = amp
        Else
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 1
           e(i2) = amp
        End If
        Goto 200
     End If
     lb(i1) = 1
     e(i1) = amp
     lb(i2) = 2
     e(i2) = amn
     ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
     ! ENERGY CONSERVATION
     ! resonance production or absorption in resonance+resonance collisions is
     ! assumed to have the same pt distribution as pp
200  em1 = e(i1)
     em2 = e(i2)
     pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
     If (pr2<=0.) pr2 = 1.E-09
     pr = sqrt(pr2)/(2.*srt)
     If (srt<=2.14) c1 = 1.0 - 2.0*ranart(nseed)
     If (srt>2.14 .And. srt<=2.4) c1 = ang(srt, iseed)
     If (srt>2.4) Then

        !lin-10/25/02 get rid of argument usage mismatch in PTR():
        xptr = 0.33*pr
        !         cc1=ptr(0.33*pr,iseed)
        cc1 = ptr(xptr, iseed)
        !lin-10/25/02-end

        !lin-9/2012: check argument in sqrt():
        scheck = pr**2 - cc1**2
        If (scheck<0) Then
           Write (99, *) 'scheck7: ', scheck
           scheck = 0.
        End If
        c1 = sqrt(scheck)/pr
        !         c1=sqrt(pr**2-cc1**2)/pr

     End If
     t1 = 2.0*pi*ranart(nseed)
     If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
        lb(i1) = -lb(i1)
        lb(i2) = -lb(i2)
     End If
  End If
  !COM: SET THE NEW MOMENTUM COORDINATES

  !lin-9/2012: check argument in sqrt():
107 scheck = 1.0 - c1**2
  If (scheck<0) Then
     Write (99, *) 'scheck8: ', scheck
     scheck = 0.
  End If
  s1 = sqrt(scheck)
  !107   S1   = SQRT( 1.0 - C1**2 )

  !lin-9/2012: check argument in sqrt():
  scheck = 1.0 - c2**2
  If (scheck<0) Then
     Write (99, *) 'scheck9: ', scheck
     scheck = 0.
  End If
  s2 = sqrt(scheck)
  !      S2  =  SQRT( 1.0 - C2**2 )

  ct1 = cos(t1)
  st1 = sin(t1)
  ct2 = cos(t2)
  st2 = sin(t2)
  pz = pr*(c1*c2-s1*s2*ct1)
  ss = c2*s1*ct1 + s2*c1
  px = pr*(ss*ct2-s1*st1*st2)
  py = pr*(ss*st2+s1*st1*ct2)
  Return
  ! FOR THE DD-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN
  ! THE NUCLEUS-NUCLEUS CMS.
306 Continue
  !sp11/21/01 phi production
  If (xsk5/sigk>ranart(nseed)) Then
     pz1 = p(3, i1)
     pz2 = p(3, i2)
     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 1 + int(2*ranart(nseed))
     nnn = nnn + 1
     lpion(nnn, irun) = 29
     epion(nnn, irun) = aphi
     iblock = 222
     Goto 208
  End If
  iblock = 10
  If (ianti==1) iblock = -10
  pz1 = p(3, i1)
  pz2 = p(3, i2)
  ! DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  nnn = nnn + 1
  lpion(nnn, irun) = 23
  epion(nnn, irun) = aka
  If (srt<=2.63) Then
     ! only lambda production is possible
     ! (1.1)P+P-->p+L+kaon+
     ic = 1
     lb(i1) = 1 + int(2*ranart(nseed))
     lb(i2) = 14
     Goto 208
  End If
  If (srt<=2.74 .And. srt>2.63) Then
     ! both Lambda and sigma production are possible
     If (xsk1/(xsk1+xsk2)>ranart(nseed)) Then
        ! lambda production
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
     Else
        ! sigma production
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 15 + int(3*ranart(nseed))
        ic = 2
     End If
     Goto 208
  End If
  If (srt<=2.77 .And. srt>2.74) Then
     ! then pp-->Delta lamda kaon can happen
     If (xsk1/(xsk1+xsk2+xsk3)>ranart(nseed)) Then
        ! * (1.1)P+P-->p+L+kaon+
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk2/(xsk2+xsk3)>ranart(nseed)) Then
           ! pp-->psk
           ic = 2
           lb(i1) = 1 + int(2*ranart(nseed))
           lb(i2) = 15 + int(3*ranart(nseed))
        Else
           ! pp-->D+l+k
           ic = 3
           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
        End If
        Goto 208
     End If
  End If
  If (srt>2.77) Then
     ! all four channels are possible
     If (xsk1/(xsk1+xsk2+xsk3+xsk4)>ranart(nseed)) Then
        ! p lambda k production
        ic = 1
        lb(i1) = 1 + int(2*ranart(nseed))
        lb(i2) = 14
        Goto 208
     Else
        If (xsk3/(xsk2+xsk3+xsk4)>ranart(nseed)) Then
           ! delta l K production
           ic = 3
           lb(i1) = 6 + int(4*ranart(nseed))
           lb(i2) = 14
           Goto 208
        Else
           If (xsk2/(xsk2+xsk4)>ranart(nseed)) Then
              ! n sigma k production
              lb(i1) = 1 + int(2*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))
              ic = 2
           Else
              ! D sigma K
              ic = 4
              lb(i1) = 6 + int(4*ranart(nseed))
              lb(i2) = 15 + int(3*ranart(nseed))
           End If
           Goto 208
        End If
     End If
  End If
208 Continue
  If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
     lb(i1) = -lb(i1)
     lb(i2) = -lb(i2)
     If (lpion(nnn,irun)==23) lpion(nnn, irun) = 21
  End If
  lbi1 = lb(i1)
  lbi2 = lb(i2)
  ! KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE
  ntry1 = 0
129 Call bbkaon(ic, srt, px3, py3, pz3, dm3, px4, py4, pz4, dm4, ppx, ppy, ppz, icou1)
  ntry1 = ntry1 + 1
  If ((icou1<0) .And. (ntry1<=20)) Goto 129
  !       if(icou1.lt.0)return
  ! ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
  Call rotate(px, py, pz, px3, py3, pz3)
  Call rotate(px, py, pz, px4, py4, pz4)
  Call rotate(px, py, pz, ppx, ppy, ppz)
  ! FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
  ! NUCLEUS CMS. FRAME
  ! (1) for the necleon/delta
  !             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
  e1cm = sqrt(dm3**2+px3**2+py3**2+pz3**2)
  p1beta = px3*betax + py3*betay + pz3*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  pt1i1 = betax*transf + px3
  pt2i1 = betay*transf + py3
  pt3i1 = betaz*transf + pz3
  eti1 = dm3
  ! (2) for the lambda/sigma
  e2cm = sqrt(dm4**2+px4**2+py4**2+pz4**2)
  p2beta = px4*betax + py4*betay + pz4*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf + px4
  pt2i2 = betay*transf + py4
  pt3i2 = betaz*transf + pz4
  eti2 = dm4
  ! GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
  epcm = sqrt(aka**2+ppx**2+ppy**2+ppz**2)
  ppbeta = ppx*betax + ppy*betay + ppz*betaz
  transf = gamma*(gamma*ppbeta/(gamma+1.)+epcm)
  ppion(1, nnn, irun) = betax*transf + ppx
  ppion(2, nnn, irun) = betay*transf + ppy
  ppion(3, nnn, irun) = betaz*transf + ppz
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  !lin-5/2008:
  !2007        X01 = 1.0 - 2.0 * RANART(NSEED)
  !            Y01 = 1.0 - 2.0 * RANART(NSEED)
  !            Z01 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2007
  !                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
  !                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
  !                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
  rpion(1, nnn, irun) = r(1, i1)
  rpion(2, nnn, irun) = r(2, i1)
  rpion(3, nnn, irun) = r(3, i1)
  !
  ! assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the
  ! leadng particle behaviour
  !              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
  e(i1) = eti1
  lb(i1) = lbi1
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
  e(i2) = eti2
  lb(i2) = lbi2
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  id(i1) = 2
  id(i2) = 2
  id1 = id(i1)
  lb1 = lb(i1)
  lb2 = lb(i2)
  am1 = em1
  am2 = em2
  e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
  Return

  !lin-6/2008 D+D->Deuteron+pi:
  !     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
108 Continue
  If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
     !     For idpert=1: we produce npertd pert deuterons:
     ndloop = npertd
  Else If (idpert==2 .And. npertd>=1) Then
     !     For idpert=2: we first save information for npertd pert deuterons;
     !     at the last ndloop we create the regular deuteron+pi
     !     and those pert deuterons:
     ndloop = npertd + 1
  Else
     !     Just create the regular deuteron+pi:
     ndloop = 1
  End If
  !
  dprob1 = sdprod/sig/float(npertd)
  Do idloop = 1, ndloop
     Call bbdangle(pxd, pyd, pzd, nt, ipert1, ianti, idloop, pfinal, dprob1, lbm)
     Call rotate(px, py, pz, pxd, pyd, pzd)
     !     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE
     !     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME:
     !     For the Deuteron:
     xmass = xmd
     e1dcm = sqrt(xmass**2+pxd**2+pyd**2+pzd**2)
     p1dbeta = pxd*betax + pyd*betay + pzd*betaz
     transf = gamma*(gamma*p1dbeta/(gamma+1.)+e1dcm)
     pxi1 = betax*transf + pxd
     pyi1 = betay*transf + pyd
     pzi1 = betaz*transf + pzd
     If (ianti==0) Then
        lbd = 42
     Else
        lbd = -42
     End If
     If (idpert==1 .And. ipert1==1 .And. npertd>=1) Then
        !ccc  Perturbative production for idpert=1:
        nnn = nnn + 1
        ppion(1, nnn, irun) = pxi1
        ppion(2, nnn, irun) = pyi1
        ppion(3, nnn, irun) = pzi1
        epion(nnn, irun) = xmd
        lpion(nnn, irun) = lbd
        rpion(1, nnn, irun) = r(1, i1)
        rpion(2, nnn, irun) = r(2, i1)
        rpion(3, nnn, irun) = r(3, i1)
        !lin-6/2008 assign the perturbative probability:
        dppion(nnn, irun) = sdprod/sig/float(npertd)
     Else If (idpert==2 .And. idloop<=npertd) Then
        !lin-6/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons
        !     only when a regular (anti)deuteron+pi is produced in NN collisions.
        !     First save the info for the perturbative deuterons:
        ppd(1, idloop) = pxi1
        ppd(2, idloop) = pyi1
        ppd(3, idloop) = pzi1
        lbpd(idloop) = lbd
     Else
        !ccc  Regular production:
        !     For the regular pion: do LORENTZ-TRANSFORMATION:
        e(i1) = xmm
        e2picm = sqrt(xmm**2+pxd**2+pyd**2+pzd**2)
        p2pibeta = -pxd*betax - pyd*betay - pzd*betaz
        transf = gamma*(gamma*p2pibeta/(gamma+1.)+e2picm)
        pxi2 = betax*transf - pxd
        pyi2 = betay*transf - pyd
        pzi2 = betaz*transf - pzd
        p(1, i1) = pxi2
        p(2, i1) = pyi2
        p(3, i1) = pzi2
        !     Remove regular pion to check the equivalence
        !     between the perturbative and regular deuteron results:
        !                 E(i1)=0.
        !
        lb(i1) = lbm
        px1 = p(1, i1)
        py1 = p(2, i1)
        pz1 = p(3, i1)
        em1 = e(i1)
        id(i1) = 2
        id1 = id(i1)
        e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
        lb1 = lb(i1)
        !     For the regular deuteron:
        p(1, i2) = pxi1
        p(2, i2) = pyi1
        p(3, i2) = pzi1
        lb(i2) = lbd
        lb2 = lb(i2)
        e(i2) = xmd
        eti2 = e(i2)
        id(i2) = 2
        !     For idpert=2: create the perturbative deuterons:
        If (idpert==2 .And. idloop==ndloop) Then
           Do ipertd = 1, npertd
              nnn = nnn + 1
              ppion(1, nnn, irun) = ppd(1, ipertd)
              ppion(2, nnn, irun) = ppd(2, ipertd)
              ppion(3, nnn, irun) = ppd(3, ipertd)
              epion(nnn, irun) = xmd
              lpion(nnn, irun) = lbpd(ipertd)
              rpion(1, nnn, irun) = r(1, i1)
              rpion(2, nnn, irun) = r(2, i1)
              rpion(3, nnn, irun) = r(3, i1)
              !lin-6/2008 assign the perturbative probability:
              dppion(nnn, irun) = 1./float(npertd)
           End Do
        End If
     End If
  End Do
  iblock = 501
  Return
  !lin-6/2008 D+D->Deuteron+pi over

End Subroutine crdd
!*********************************
!*********************************
!                                                                      *
Subroutine init(minnum, maxnum, num, radius, x0, z0, p0, gamma, iseed, mass, iopt)
  !                                                                      *
  !       PURPOSE:     PROVIDING INITIAL CONDITIONS FOR PHASE-SPACE      *
  !                    DISTRIBUTION OF TESTPARTICLES                     *
  !       VARIABLES:   (ALL INPUT)                                       *
  !         MINNUM  - FIRST TESTPARTICLE TREATED IN ONE RUN    (INTEGER) *
  !         MAXNUM  - LAST TESTPARTICLE TREATED IN ONE RUN     (INTEGER) *
  !         NUM     - NUMBER OF TESTPARTICLES PER NUCLEON      (INTEGER) *
  !         RADIUS  - RADIUS OF NUCLEUS "FM"                      (REAL) *
  !         X0,Z0   - DISPLACEMENT OF CENTER OF NUCLEUS IN X,Z-          *
  !                   DIRECTION "FM"                              (REAL) *
  !         P0      - MOMENTUM-BOOST IN C.M. FRAME "GEV/C"        (REAL) *
  !         GAMMA   - RELATIVISTIC GAMMA-FACTOR                   (REAL) *
  !         ISEED   - SEED FOR RANDOM-NUMBER GENERATOR         (INTEGER) *
  !         MASS    - TOTAL MASS OF THE SYSTEM                 (INTEGER) *
  !         IOPT    - OPTION FOR DIFFERENT OCCUPATION OF MOMENTUM        *
  !                   SPACE                                    (INTEGER) *
  !                                                                      *
  !*********************************
  Parameter (maxstr=150001, amu=0.9383)
  Parameter (maxx=20, maxz=24)
  Parameter (pi=3.1415926)
  !
  Real ptot(3)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  !----------------------------------------------------------------------
  !     PREPARATION FOR LORENTZ-TRANSFORMATIONS
  !
  If (p0/=0.) Then
     sign = p0/abs(p0)
  Else
     sign = 0.
  End If

  !lin-9/2012: check argument in sqrt():
  scheck = gamma**2 - 1.
  If (scheck<0) Then
     Write (99, *) 'scheck10: ', scheck
     scheck = 0.
  End If
  beta = sign*sqrt(scheck)/gamma
  !      BETA = SIGN * SQRT(GAMMA**2-1.)/GAMMA

  !-----------------------------------------------------------------------
  !     TARGET-ID = 1 AND PROJECTILE-ID = -1
  !
  If (minnum==1) Then
     idnum = 1
  Else
     idnum = -1
  End If
  !-----------------------------------------------------------------------
  !     IDENTIFICATION OF TESTPARTICLES AND ASSIGMENT OF RESTMASS
  !
  !     LOOP OVER ALL PARALLEL RUNS:
  Do irun = 1, num
     Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
        id(i) = idnum
        e(i) = amu
     End Do
     !-----------------------------------------------------------------------
     !       OCCUPATION OF COORDINATE-SPACE
     !
     Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
200     Continue
        x = 1.0 - 2.0*ranart(nseed)
        y = 1.0 - 2.0*ranart(nseed)
        z = 1.0 - 2.0*ranart(nseed)
        If ((x*x+y*y+z*z)>1.0) Goto 200
        r(1, i) = x*radius
        r(2, i) = y*radius
        r(3, i) = z*radius
     End Do
  End Do
  !=======================================================================
  If (iopt/=3) Then
     !-----
     !     OPTION 1: USE WOODS-SAXON PARAMETRIZATION FOR DENSITY AND
     !-----          CALCULATE LOCAL FERMI-MOMENTUM
     !
     rhow0 = 0.168
     Do irun = 1, num
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
500        Continue
           px = 1.0 - 2.0*ranart(nseed)
           py = 1.0 - 2.0*ranart(nseed)
           pz = 1.0 - 2.0*ranart(nseed)
           If (px*px+py*py+pz*pz>1.0) Goto 500
           rdist = sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
           rhows = rhow0/(1.0+exp((rdist-radius)/0.55))
           pfermi = 0.197*(1.5*pi*pi*rhows)**(1./3.)
           !-----
           !     OPTION 2: NUCLEAR MATTER CASE
           If (iopt==2) pfermi = 0.27
           If (iopt==4) pfermi = 0.
           !-----
           p(1, i) = pfermi*px
           p(2, i) = pfermi*py
           p(3, i) = pfermi*pz
        End Do
        !
        !         SET TOTAL MOMENTUM TO 0 IN REST FRAME AND BOOST
        !
        Do idir = 1, 3
           ptot(idir) = 0.0
        End Do
        npart = 0
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
           npart = npart + 1
           Do idir = 1, 3
              ptot(idir) = ptot(idir) + p(idir, i)
           End Do
        End Do
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
           Do idir = 1, 3
              p(idir, i) = p(idir, i) - ptot(idir)/float(npart)
           End Do
           !           BOOST
           If ((iopt==1) .Or. (iopt==2)) Then
              epart = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+amu**2)
              p(3, i) = gamma*(p(3,i)+beta*epart)
           Else
              p(3, i) = p(3, i) + p0
           End If
        End Do
     End Do
     !-----
  Else
     !-----
     !     OPTION 3: GIVE ALL NUCLEONS JUST A Z-MOMENTUM ACCORDING TO
     !               THE BOOST OF THE NUCLEI
     !
     Do irun = 1, num
        Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
           p(1, i) = 0.0
           p(2, i) = 0.0
           p(3, i) = p0
        End Do
     End Do
     !-----
  End If
  !=======================================================================
  !     PUT PARTICLES IN THEIR POSITION IN COORDINATE-SPACE
  !     (SHIFT AND RELATIVISTIC CONTRACTION)
  !
  Do irun = 1, num
     Do i = minnum + (irun-1)*mass, maxnum + (irun-1)*mass
        r(1, i) = r(1, i) + x0
        ! two nuclei in touch after contraction
        r(3, i) = (r(3,i)+z0)/gamma
        ! two nuclei in touch before contraction
        !          R(3,I) = R(3,I) / GAMMA + Z0
     End Do
  End Do
  !
  Return
End Subroutine init
!*********************************
!                                                                      *
Subroutine dens(ipot, mass, num, nesc)
  !                                                                      *
  !       PURPOSE:     CALCULATION OF LOCAL BARYON, MESON AND ENERGY     *
  !                    DENSITY FROM SPATIAL DISTRIBUTION OF TESTPARTICLES*
  !                                                                      *
  !       VARIABLES (ALL INPUT, ALL INTEGER)                             *
  !         MASS    -  MASS NUMBER OF THE SYSTEM                         *
  !         NUM     -  NUMBER OF TESTPARTICLES PER NUCLEON               *
  !                                                                      *
  !         NESC    -  NUMBER OF ESCAPED PARTICLES      (INTEGER,OUTPUT) *
  !                                                                      *
  !*********************************
  Parameter (maxstr=150001, maxr=1)
  Parameter (maxx=20, maxz=24)
  !
  Dimension pxl(-maxx:maxx, -maxx:maxx, -maxz:maxz), pyl(-maxx:maxx, -maxx:maxx, -maxz:maxz), pzl(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ddpi/pirho(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DDpi/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Common /rr/massr(0:maxr)
  !c      SAVE /RR/
  Common /tt/pel(-maxx:maxx, -maxx:maxx, -maxz:maxz), rxy(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /tt/
  Common /bbb/bxx(-maxx:maxx, -maxx:maxx, -maxz:maxz), byy(-maxx:maxx, -maxx:maxx, -maxz:maxz), bzz(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !
  Real zet(-45:45)
  Save
  Data zet/1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., -1., 0., -2., -1., 0., 1., 0., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., -1., 0., 1., 2., 0., 1., 0., 1., 0., -1., 0., 1., 0., 0., 0., -1., 0., 1., 0., -1., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1./

  Do iz = -maxz, maxz
     Do iy = -maxx, maxx
        Do ix = -maxx, maxx
           rho(ix, iy, iz) = 0.0
           rhon(ix, iy, iz) = 0.0
           rhop(ix, iy, iz) = 0.0
           pirho(ix, iy, iz) = 0.0
           pxl(ix, iy, iz) = 0.0
           pyl(ix, iy, iz) = 0.0
           pzl(ix, iy, iz) = 0.0
           pel(ix, iy, iz) = 0.0
           bxx(ix, iy, iz) = 0.0
           byy(ix, iy, iz) = 0.0
           bzz(ix, iy, iz) = 0.0
        End Do
     End Do
  End Do
  !
  nesc = 0
  big = 1.0/(3.0*float(num))
  small = 1.0/(9.0*float(num))
  !
  msum = 0
  Do irun = 1, num
     msum = msum + massr(irun-1)
     Do j = 1, massr(irun)
        i = j + msum
        ix = nint(r(1,i))
        iy = nint(r(2,i))
        iz = nint(r(3,i))
        If (ix<=-maxx .Or. ix>=maxx .Or. iy<=-maxx .Or. iy>=maxx .Or. iz<=-maxz .Or. iz>=maxz) Then
           nesc = nesc + 1
        Else
           !
           !sp01/04/02 include baryon density
           If (j>mass) Goto 30
           !         if( (lb(i).eq.1.or.lb(i).eq.2) .or.
           !    &    (lb(i).ge.6.and.lb(i).le.17) )then
           ! (1) baryon density
           rho(ix, iy, iz) = rho(ix, iy, iz) + big
           rho(ix+1, iy, iz) = rho(ix+1, iy, iz) + small
           rho(ix-1, iy, iz) = rho(ix-1, iy, iz) + small
           rho(ix, iy+1, iz) = rho(ix, iy+1, iz) + small
           rho(ix, iy-1, iz) = rho(ix, iy-1, iz) + small
           rho(ix, iy, iz+1) = rho(ix, iy, iz+1) + small
           rho(ix, iy, iz-1) = rho(ix, iy, iz-1) + small
           ! (2) CALCULATE THE PROTON DENSITY
           If (zet(lb(i))/=0) Then
              rhop(ix, iy, iz) = rhop(ix, iy, iz) + big
              rhop(ix+1, iy, iz) = rhop(ix+1, iy, iz) + small
              rhop(ix-1, iy, iz) = rhop(ix-1, iy, iz) + small
              rhop(ix, iy+1, iz) = rhop(ix, iy+1, iz) + small
              rhop(ix, iy-1, iz) = rhop(ix, iy-1, iz) + small
              rhop(ix, iy, iz+1) = rhop(ix, iy, iz+1) + small
              rhop(ix, iy, iz-1) = rhop(ix, iy, iz-1) + small
              Goto 40
           End If
           ! (3) CALCULATE THE NEUTRON DENSITY
           If (zet(lb(i))==0) Then
              rhon(ix, iy, iz) = rhon(ix, iy, iz) + big
              rhon(ix+1, iy, iz) = rhon(ix+1, iy, iz) + small
              rhon(ix-1, iy, iz) = rhon(ix-1, iy, iz) + small
              rhon(ix, iy+1, iz) = rhon(ix, iy+1, iz) + small
              rhon(ix, iy-1, iz) = rhon(ix, iy-1, iz) + small
              rhon(ix, iy, iz+1) = rhon(ix, iy, iz+1) + small
              rhon(ix, iy, iz-1) = rhon(ix, iy, iz-1) + small
              Goto 40
           End If
           !           else    !! sp01/04/02
           ! (4) meson density
30         pirho(ix, iy, iz) = pirho(ix, iy, iz) + big
           pirho(ix+1, iy, iz) = pirho(ix+1, iy, iz) + small
           pirho(ix-1, iy, iz) = pirho(ix-1, iy, iz) + small
           pirho(ix, iy+1, iz) = pirho(ix, iy+1, iz) + small
           pirho(ix, iy-1, iz) = pirho(ix, iy-1, iz) + small
           pirho(ix, iy, iz+1) = pirho(ix, iy, iz+1) + small
           pirho(ix, iy, iz-1) = pirho(ix, iy, iz-1) + small
           !           endif    !! sp01/04/02
           ! to calculate the Gamma factor in each cell
           !(1) PX
40         pxl(ix, iy, iz) = pxl(ix, iy, iz) + p(1, i)*big
           pxl(ix+1, iy, iz) = pxl(ix+1, iy, iz) + p(1, i)*small
           pxl(ix-1, iy, iz) = pxl(ix-1, iy, iz) + p(1, i)*small
           pxl(ix, iy+1, iz) = pxl(ix, iy+1, iz) + p(1, i)*small
           pxl(ix, iy-1, iz) = pxl(ix, iy-1, iz) + p(1, i)*small
           pxl(ix, iy, iz+1) = pxl(ix, iy, iz+1) + p(1, i)*small
           pxl(ix, iy, iz-1) = pxl(ix, iy, iz-1) + p(1, i)*small
           !(2) PY
           pyl(ix, iy, iz) = pyl(ix, iy, iz) + p(2, i)*big
           pyl(ix+1, iy, iz) = pyl(ix+1, iy, iz) + p(2, i)*small
           pyl(ix-1, iy, iz) = pyl(ix-1, iy, iz) + p(2, i)*small
           pyl(ix, iy+1, iz) = pyl(ix, iy+1, iz) + p(2, i)*small
           pyl(ix, iy-1, iz) = pyl(ix, iy-1, iz) + p(2, i)*small
           pyl(ix, iy, iz+1) = pyl(ix, iy, iz+1) + p(2, i)*small
           pyl(ix, iy, iz-1) = pyl(ix, iy, iz-1) + p(2, i)*small
           ! (3) PZ
           pzl(ix, iy, iz) = pzl(ix, iy, iz) + p(3, i)*big
           pzl(ix+1, iy, iz) = pzl(ix+1, iy, iz) + p(3, i)*small
           pzl(ix-1, iy, iz) = pzl(ix-1, iy, iz) + p(3, i)*small
           pzl(ix, iy+1, iz) = pzl(ix, iy+1, iz) + p(3, i)*small
           pzl(ix, iy-1, iz) = pzl(ix, iy-1, iz) + p(3, i)*small
           pzl(ix, iy, iz+1) = pzl(ix, iy, iz+1) + p(3, i)*small
           pzl(ix, iy, iz-1) = pzl(ix, iy, iz-1) + p(3, i)*small
           ! (4) ENERGY
           pel(ix, iy, iz) = pel(ix, iy, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*big
           pel(ix+1, iy, iz) = pel(ix+1, iy, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix-1, iy, iz) = pel(ix-1, iy, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy+1, iz) = pel(ix, iy+1, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy-1, iz) = pel(ix, iy-1, iz) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy, iz+1) = pel(ix, iy, iz+1) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
           pel(ix, iy, iz-1) = pel(ix, iy, iz-1) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)*small
        End If
     End Do
  End Do
  !
  Do iz = -maxz, maxz
     Do iy = -maxx, maxx
        Do ix = -maxx, maxx
           If ((rho(ix,iy,iz)==0) .Or. (pel(ix,iy,iz)==0)) Goto 101
           smass2 = pel(ix, iy, iz)**2 - pxl(ix, iy, iz)**2 - pyl(ix, iy, iz)**2 - pzl(ix, iy, iz)**2
           If (smass2<=0) smass2 = 1.E-06
           smass = sqrt(smass2)
           If (smass==0.) smass = 1.E-06
           gamma = pel(ix, iy, iz)/smass
           If (gamma==0) Goto 101
           bxx(ix, iy, iz) = pxl(ix, iy, iz)/pel(ix, iy, iz)
           byy(ix, iy, iz) = pyl(ix, iy, iz)/pel(ix, iy, iz)
           bzz(ix, iy, iz) = pzl(ix, iy, iz)/pel(ix, iy, iz)
           rho(ix, iy, iz) = rho(ix, iy, iz)/gamma
           rhon(ix, iy, iz) = rhon(ix, iy, iz)/gamma
           rhop(ix, iy, iz) = rhop(ix, iy, iz)/gamma
           pirho(ix, iy, iz) = pirho(ix, iy, iz)/gamma
           pel(ix, iy, iz) = pel(ix, iy, iz)/(gamma**2)
           rho0 = 0.163
           If (ipot==0) Then
              u = 0
              Goto 70
           End If
           If (ipot==1 .Or. ipot==6) Then
              a = -0.1236
              b = 0.0704
              s = 2
              Goto 60
           End If
           If (ipot==2 .Or. ipot==7) Then
              a = -0.218
              b = 0.164
              s = 4./3.
           End If
           If (ipot==3) Then
              a = -0.3581
              b = 0.3048
              s = 1.167
              Goto 60
           End If
           If (ipot==4) Then
              denr = rho(ix, iy, iz)/rho0
              b = 0.3048
              s = 1.167
              If (denr<=4 .Or. denr>7) Then
                 a = -0.3581
              Else
                 a = -b*denr**(1./6.) - 2.*0.036/3.*denr**(-0.333)
              End If
              Goto 60
           End If
60         u = 0.5*a*rho(ix, iy, iz)**2/rho0 + b/(1+s)*(rho(ix,iy,iz)/rho0)**s*rho(ix, iy, iz)
70         pel(ix, iy, iz) = pel(ix, iy, iz) + u
101     End Do
     End Do
  End Do
  Return
End Subroutine dens

!*********************************
!                                                                      *
Subroutine gradu(iopt, ix, iy, iz, gradx, grady, gradz)
  !                                                                      *
  !       PURPOSE:     DETERMINE GRAD(U(RHO(X,Y,Z)))                     *
  !       VARIABLES:                                                     *
  !         IOPT                - METHOD FOR EVALUATING THE GRADIENT     *
  !                                                      (INTEGER,INPUT) *
  !         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) *
  !         GRADX, GRADY, GRADZ - GRADIENT OF U            (REAL,OUTPUT) *
  !                                                                      *
  !*********************************
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.167)
  !
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Common /tt/pel(-maxx:maxx, -maxx:maxx, -maxz:maxz), rxy(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /tt/
  Save
  !
  rxplus = rho(ix+1, iy, iz)/rho0
  rxmins = rho(ix-1, iy, iz)/rho0
  ryplus = rho(ix, iy+1, iz)/rho0
  rymins = rho(ix, iy-1, iz)/rho0
  rzplus = rho(ix, iy, iz+1)/rho0
  rzmins = rho(ix, iy, iz-1)/rho0
  den0 = rho(ix, iy, iz)/rho0
  ene0 = pel(ix, iy, iz)
  !-----------------------------------------------------------------------
  Goto (1, 2, 3, 4, 5) iopt
  If (iopt==6) Goto 6
  If (iopt==7) Goto 7
  !
1 Continue
  !       POTENTIAL USED IN 1) (STIFF):
  !       U = -.124 * RHO/RHO0 + .0705 (RHO/RHO0)**2 GEV
  !
  gradx = -0.062*(rxplus-rxmins) + 0.03525*(rxplus**2-rxmins**2)
  grady = -0.062*(ryplus-rymins) + 0.03525*(ryplus**2-rymins**2)
  gradz = -0.062*(rzplus-rzmins) + 0.03525*(rzplus**2-rzmins**2)
  Return
  !
2 Continue
  !       POTENTIAL USED IN 2):
  !       U = -.218 * RHO/RHO0 + .164 (RHO/RHO0)**(4/3) GEV
  !
  expnt = 1.3333333
  gradx = -0.109*(rxplus-rxmins) + 0.082*(rxplus**expnt-rxmins**expnt)
  grady = -0.109*(ryplus-rymins) + 0.082*(ryplus**expnt-rymins**expnt)
  gradz = -0.109*(rzplus-rzmins) + 0.082*(rzplus**expnt-rzmins**expnt)
  Return
  !
3 Continue
  !       POTENTIAL USED IN 3) (SOFT):
  !       U = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV
  !
  expnt = 1.1666667
  acoef = 0.178
  gradx = -acoef*(rxplus-rxmins) + 0.1515*(rxplus**expnt-rxmins**expnt)
  grady = -acoef*(ryplus-rymins) + 0.1515*(ryplus**expnt-rymins**expnt)
  gradz = -acoef*(rzplus-rzmins) + 0.1515*(rzplus**expnt-rzmins**expnt)
  Return
  !
  !
4 Continue
  !       POTENTIAL USED IN 4) (super-soft in the mixed phase of 4 < rho/rho <7):
  !       U1 = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV
  !       normal phase, soft eos of iopt=3
  !       U2 = -.02 * (RHO/RHO0)**(2/3) -0.0253 * (RHO/RHO0)**(7/6)  GEV
  !
  eh = 4.
  eqgp = 7.
  acoef = 0.178
  expnt = 1.1666667
  denr = rho(ix, iy, iz)/rho0
  If (denr<=eh .Or. denr>=eqgp) Then
     gradx = -acoef*(rxplus-rxmins) + 0.1515*(rxplus**expnt-rxmins**expnt)
     grady = -acoef*(ryplus-rymins) + 0.1515*(ryplus**expnt-rymins**expnt)
     gradz = -acoef*(rzplus-rzmins) + 0.1515*(rzplus**expnt-rzmins**expnt)
  Else
     acoef1 = 0.178
     acoef2 = 0.0
     expnt2 = 2./3.
     gradx = -acoef1*(rxplus**expnt-rxmins**expnt) - acoef2*(rxplus**expnt2-rxmins**expnt2)
     grady = -acoef1*(ryplus**expnt-rymins**expnt) - acoef2*(ryplus**expnt2-rymins**expnt2)
     gradz = -acoef1*(rzplus**expnt-rzmins**expnt) - acoef2*(rzplus**expnt2-rzmins**expnt2)
  End If
  Return
  !
5 Continue
  !       POTENTIAL USED IN 5) (SUPER STIFF):
  !       U = -.10322 * RHO/RHO0 + .04956 * (RHO/RHO0)**(2.77)  GEV
  !
  expnt = 2.77
  gradx = -0.0516*(rxplus-rxmins) + 0.02498*(rxplus**expnt-rxmins**expnt)
  grady = -0.0516*(ryplus-rymins) + 0.02498*(ryplus**expnt-rymins**expnt)
  gradz = -0.0516*(rzplus-rzmins) + 0.02498*(rzplus**expnt-rzmins**expnt)
  Return
  !
6 Continue
  !       POTENTIAL USED IN 6) (STIFF-qgp):
  !       U = -.124 * RHO/RHO0 + .0705 (RHO/RHO0)**2 GEV
  !
  If (ene0<=0.5) Then
     gradx = -0.062*(rxplus-rxmins) + 0.03525*(rxplus**2-rxmins**2)
     grady = -0.062*(ryplus-rymins) + 0.03525*(ryplus**2-rymins**2)
     gradz = -0.062*(rzplus-rzmins) + 0.03525*(rzplus**2-rzmins**2)
     Return
  End If
  If (ene0>0.5 .And. ene0<=1.5) Then
     !       U=c1-ef*rho/rho0**2/3
     ef = 36./1000.
     gradx = -0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = -0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = -0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
  If (ene0>1.5) Then
     ! U=800*(rho/rho0)**1/3.-Ef*(rho/rho0)**2/3.-c2
     ef = 36./1000.
     cf0 = 0.8
     gradx = 0.5*cf0*(rxplus**0.333-rxmins**0.333) - 0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = 0.5*cf0*(ryplus**0.333-rymins**0.333) - 0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = 0.5*cf0*(rzplus**0.333-rzmins**0.333) - 0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
  !
7 Continue
  !       POTENTIAL USED IN 7) (Soft-qgp):
  If (den0<=4.5) Then
     !       POTENTIAL USED is the same as IN 3) (SOFT):
     !       U = -.356 * RHO/RHO0 + .303 * (RHO/RHO0)**(7/6)  GEV
     !
     expnt = 1.1666667
     acoef = 0.178
     gradx = -acoef*(rxplus-rxmins) + 0.1515*(rxplus**expnt-rxmins**expnt)
     grady = -acoef*(ryplus-rymins) + 0.1515*(ryplus**expnt-rymins**expnt)
     gradz = -acoef*(rzplus-rzmins) + 0.1515*(rzplus**expnt-rzmins**expnt)
     Return
  End If
  If (den0>4.5 .And. den0<=5.1) Then
     !       U=c1-ef*rho/rho0**2/3
     ef = 36./1000.
     gradx = -0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = -0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = -0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
  If (den0>5.1) Then
     ! U=800*(rho/rho0)**1/3.-Ef*(rho/rho0)**2/3.-c2
     ef = 36./1000.
     cf0 = 0.8
     gradx = 0.5*cf0*(rxplus**0.333-rxmins**0.333) - 0.5*ef*(rxplus**0.67-rxmins**0.67)
     grady = 0.5*cf0*(ryplus**0.333-rymins**0.333) - 0.5*ef*(ryplus**0.67-rymins**0.67)
     gradz = 0.5*cf0*(rzplus**0.333-rzmins**0.333) - 0.5*ef*(rzplus**0.67-rzmins**0.67)
     Return
  End If
End Subroutine gradu
!*********************************
!                                                                      *
Subroutine graduk(ix, iy, iz, gradxk, gradyk, gradzk)
  !                                                                      *
  !       PURPOSE:     DETERMINE the baryon density gradient for         *
  !                    proporgating kaons in a mean field caused by      *
  !                    surrounding baryons                               *
  !       VARIABLES:                                                     *
  !         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) *
  !         GRADXk, GRADYk, GRADZk                       (REAL,OUTPUT)   *
  !                                                                      *
  !*********************************
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.168)
  !
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Save
  !
  rxplus = rho(ix+1, iy, iz)
  rxmins = rho(ix-1, iy, iz)
  ryplus = rho(ix, iy+1, iz)
  rymins = rho(ix, iy-1, iz)
  rzplus = rho(ix, iy, iz+1)
  rzmins = rho(ix, iy, iz-1)
  gradxk = (rxplus-rxmins)/2.
  gradyk = (ryplus-rymins)/2.
  gradzk = (rzplus-rzmins)/2.
  Return
End Subroutine graduk
!-----------------------------------------------------------------------
Subroutine gradup(ix, iy, iz, gradxp, gradyp, gradzp)
  !                                                                      *
  !       PURPOSE:     DETERMINE THE GRADIENT OF THE PROTON DENSITY      *
  !       VARIABLES:                                                     *
  !                                                                           *
  !         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) *
  !         GRADXP, GRADYP, GRADZP - GRADIENT OF THE PROTON              *
  !                                  DENSITY(REAL,OUTPUT)                *
  !                                                                      *
  !*********************************
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.168)
  !
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Save
  !
  rxplus = rhop(ix+1, iy, iz)/rho0
  rxmins = rhop(ix-1, iy, iz)/rho0
  ryplus = rhop(ix, iy+1, iz)/rho0
  rymins = rhop(ix, iy-1, iz)/rho0
  rzplus = rhop(ix, iy, iz+1)/rho0
  rzmins = rhop(ix, iy, iz-1)/rho0
  !-----------------------------------------------------------------------
  !
  gradxp = (rxplus-rxmins)/2.
  gradyp = (ryplus-rymins)/2.
  gradzp = (rzplus-rzmins)/2.
  Return
End Subroutine gradup
!-----------------------------------------------------------------------
Subroutine gradun(ix, iy, iz, gradxn, gradyn, gradzn)
  !                                                                      *
  !       PURPOSE:     DETERMINE THE GRADIENT OF THE NEUTRON DENSITY     *
  !       VARIABLES:                                                     *
  !                                                                           *
  !         IX, IY, IZ          - COORDINATES OF POINT   (INTEGER,INPUT) *
  !         GRADXN, GRADYN, GRADZN - GRADIENT OF THE NEUTRON             *
  !                                  DENSITY(REAL,OUTPUT)                *
  !                                                                      *
  !*********************************
  Parameter (maxx=20, maxz=24)
  Parameter (rho0=0.168)
  !
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
  !c      SAVE /DD/
  Common /ss/inout(20)
  !c      SAVE /ss/
  Save
  !
  rxplus = rhon(ix+1, iy, iz)/rho0
  rxmins = rhon(ix-1, iy, iz)/rho0
  ryplus = rhon(ix, iy+1, iz)/rho0
  rymins = rhon(ix, iy-1, iz)/rho0
  rzplus = rhon(ix, iy, iz+1)/rho0
  rzmins = rhon(ix, iy, iz-1)/rho0
  !-----------------------------------------------------------------------
  !
  gradxn = (rxplus-rxmins)/2.
  gradyn = (ryplus-rymins)/2.
  gradzn = (rzplus-rzmins)/2.
  Return
End Subroutine gradun

!-----------------------------------------------------------------------------
!FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF
!KITAZOE'S FORMULA
Real Function fde(dmass, srt, con)
  Save
  amn = 0.938869
  avpi = 0.13803333
  am0 = 1.232
  fd = 4.*(am0**2)*width(dmass)/((dmass**2-1.232**2)**2+am0**2*width(dmass)**2)
  If (con==1.) Then
     p11 = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (p11<=0) p11 = 1.E-06
     p1 = sqrt(p11)
  Else
     dmass = amn + avpi
     p11 = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (p11<=0) p11 = 1.E-06
     p1 = sqrt(p11)
  End If
  fde = fd*p1*dmass
  Return
End Function fde
!-------------------------------------------------------------
!FUNCTION FDE(DMASS) GIVES N*(1535) MASS DISTRIBUTION BY USING OF
!KITAZOE'S FORMULA
Real Function fd5(dmass, srt, con)
  Save
  amn = 0.938869
  avpi = 0.13803333
  am0 = 1.535
  fd = 4.*(am0**2)*w1535(dmass)/((dmass**2-1.535**2)**2+am0**2*w1535(dmass)**2)
  If (con==1.) Then

     !lin-9/2012: check argument in sqrt():
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck11: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
     !           P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
     !     1          /(4.*SRT**2)-DMASS**2)

  Else
     dmass = amn + avpi

     !lin-9/2012: check argument in sqrt():
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck12: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
     !        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
     !     1  /(4.*SRT**2)-DMASS**2)

  End If
  fd5 = fd*p1*dmass
  Return
End Function fd5


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!*********************************
!*********************************
!                                                                      *
!                                                                      *
!
!--------------------------------------------------------------------------
!FUNCTION FNS(DMASS) GIVES N* MASS DISTRIBUTION
!     BY USING OF BREIT-WIGNER FORMULA
Real Function fns(dmass, srt, con)
  Save
  width = 0.2
  amn = 0.938869
  avpi = 0.13803333
  an0 = 1.43
  fn = 4.*(an0**2)*width/((dmass**2-1.44**2)**2+an0**2*width**2)
  If (con==1.) Then

     !lin-9/2012: check argument in sqrt():
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck13: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
     !        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
     !     1  /(4.*SRT**2)-DMASS**2)

  Else
     dmass = amn + avpi
     !lin-9/2012: check argument in sqrt():
     scheck = (srt**2+dmass**2-amn**2)**2/(4.*srt**2) - dmass**2
     If (scheck<0) Then
        Write (99, *) 'scheck14: ', scheck
        scheck = 0.
     End If
     p1 = sqrt(scheck)
     !        P1=SQRT((SRT**2+DMASS**2-AMN**2)**2
     !     1  /(4.*SRT**2)-DMASS**2)

  End If
  fns = fn*p1*dmass
  Return
End Function fns
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! PURPOSE:1. SORT N*(1440) and N*(1535) 2-body DECAY PRODUCTS
!         2. DETERMINE THE MOMENTUM AND COORDINATES OF NUCLEON AND PION
!            AFTER THE DELTA OR N* DECAYING
! DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA
Subroutine decay(irun, i, nnn, iseed, wid, nt)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, etam=0.5475, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  !c      SAVE /INPUT2/
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  !c      SAVE /RNDF77/
  Save
  lbanti = lb(i)
  !
  dm = e(i)
  !1. FOR N*+(1440) DECAY
  If (iabs(lb(i))==11) Then
     x3 = ranart(nseed)
     If (x3>(1./3.)) Then
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
     Else
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
     End If
     !2. FOR N*0(1440) DECAY
  Else If (iabs(lb(i))==10) Then
     x4 = ranart(nseed)
     If (x4>(1./3.)) Then
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
     Else
        lb(i) = 2
        nalb = 2
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
     End If
     ! N*(1535) CAN DECAY TO A PION OR AN ETA IF DM > 1.49 GeV
     !3 N*(0)(1535) DECAY
  Else If (iabs(lb(i))==12) Then
     ctrl = 0.65
     If (dm<=1.49) ctrl = -1.
     x5 = ranart(nseed)
     If (x5>=ctrl) Then
        ! DECAY TO PION+NUCLEON
        x6 = ranart(nseed)
        If (x6>(1./3.)) Then
           lb(i) = 1
           nlab = 1
           lpion(nnn, irun) = 3
           epion(nnn, irun) = ap2
        Else
           lb(i) = 2
           nalb = 2
           lpion(nnn, irun) = 4
           epion(nnn, irun) = ap1
        End If
     Else
        ! DECAY TO ETA+NEUTRON
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 0
        epion(nnn, irun) = etam
     End If
     !4. FOR N*+(1535) DECAY
  Else If (iabs(lb(i))==13) Then
     ctrl = 0.65
     If (dm<=1.49) ctrl = -1.
     x5 = ranart(nseed)
     If (x5>=ctrl) Then
        ! DECAY TO PION+NUCLEON
        x8 = ranart(nseed)
        If (x8>(1./3.)) Then
           lb(i) = 2
           nlab = 2
           lpion(nnn, irun) = 5
           epion(nnn, irun) = ap2
        Else
           lb(i) = 1
           nlab = 1
           lpion(nnn, irun) = 4
           epion(nnn, irun) = ap1
        End If
     Else
        ! DECAY TO ETA+NUCLEON
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 0
        epion(nnn, irun) = etam
     End If
  End If
  !
  Call dkine(irun, i, nnn, nlab, iseed, wid, nt)
  !
  !     anti-particle ID for anti-N* decays:
  If (lbanti<0) Then
     lbi = lb(i)
     If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     Else If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     End If
     lb(i) = lbi
     !
     lbi = lpion(nnn, irun)
     If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     Else If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     End If
     lpion(nnn, irun) = lbi
  End If
  !
  If (nt==ntmax) Then
     !     at the last timestep, assign rho or eta (decay daughter)
     !     to lb(i1) only (not to lpion) in order to decay them again:
     lbm = lpion(nnn, irun)
     If (lbm==0 .Or. lbm==25 .Or. lbm==26 .Or. lbm==27) Then
        !     switch rho or eta with baryon, positions are the same (no change needed):
        lbsave = lbm
        xmsave = epion(nnn, irun)
        pxsave = ppion(1, nnn, irun)
        pysave = ppion(2, nnn, irun)
        pzsave = ppion(3, nnn, irun)
        !lin-5/2008:
        dpsave = dppion(nnn, irun)
        lpion(nnn, irun) = lb(i)
        epion(nnn, irun) = e(i)
        ppion(1, nnn, irun) = p(1, i)
        ppion(2, nnn, irun) = p(2, i)
        ppion(3, nnn, irun) = p(3, i)
        !lin-5/2008:
        dppion(nnn, irun) = dpertp(i)
        lb(i) = lbsave
        e(i) = xmsave
        p(1, i) = pxsave
        p(2, i) = pysave
        p(3, i) = pzsave
        !lin-5/2008:
        dpertp(i) = dpsave
     End If
  End If

  Return
End Subroutine decay

!-------------------------------------------------------------------
!-------------------------------------------------------------------
! PURPOSE:
!         CALCULATE THE MOMENTUM OF NUCLEON AND PION (OR ETA)
!         IN THE LAB. FRAME AFTER DELTA OR N* DECAY
! DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA PRODUCTION
Subroutine dkine(irun, i, nnn, nlab, iseed, wid, nt)
  Parameter (hbarc=0.19733)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, etam=0.5475, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  !c      SAVE /leadng/
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  !c      SAVE /tdecay/
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  !c      SAVE /INPUT2/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  External iarflv, invflv
  Save
  ! READ IN THE COORDINATES OF DELTA OR N* UNDERGOING DECAY
  px = p(1, i)
  py = p(2, i)
  pz = p(3, i)
  rx = r(1, i)
  ry = r(2, i)
  rz = r(3, i)
  dm = e(i)
  edelta = sqrt(dm**2+px**2+py**2+pz**2)
  pm = epion(nnn, irun)
  am = amp
  If (nlab==2) am = amn
  ! FIND OUT THE MOMENTUM AND ENERGY OF PION AND NUCLEON IN DELTA REST FRAME
  ! THE MAGNITUDE OF MOMENTUM IS DETERMINED BY ENERGY CONSERVATION ,THE FORMULA
  ! CAN BE FOUND ON PAGE 716,W BAUER P.R.C40,1989
  ! THE DIRECTION OF THE MOMENTUM IS ASSUMED ISOTROPIC. NOTE THAT P(PION)=-P(N)
  q2 = ((dm**2-am**2+pm**2)/(2.*dm))**2 - pm**2
  If (q2<=0.) q2 = 1.E-09
  q = sqrt(q2)
11 qx = 1. - 2.*ranart(nseed)
  qy = 1. - 2.*ranart(nseed)
  qz = 1. - 2.*ranart(nseed)
  qs = qx**2 + qy**2 + qz**2
  If (qs>1.) Goto 11
  pxp = q*qx/sqrt(qs)
  pyp = q*qy/sqrt(qs)
  pzp = q*qz/sqrt(qs)
  ep = sqrt(q**2+pm**2)
  pxn = -pxp
  pyn = -pyp
  pzn = -pzp
  en = sqrt(q**2+am**2)
  ! TRANSFORM INTO THE LAB. FRAME. THE GENERAL LORENTZ TRANSFORMATION CAN
  ! BE FOUND ON PAGE 34 OF R. HAGEDORN " RELATIVISTIC KINEMATICS"
  gd = edelta/dm
  fgd = gd/(1.+gd)
  bdx = px/edelta
  bdy = py/edelta
  bdz = pz/edelta
  bpp = bdx*pxp + bdy*pyp + bdz*pzp
  bpn = bdx*pxn + bdy*pyn + bdz*pzn
  p(1, i) = pxn + bdx*gd*(fgd*bpn+en)
  p(2, i) = pyn + bdy*gd*(fgd*bpn+en)
  p(3, i) = pzn + bdz*gd*(fgd*bpn+en)
  e(i) = am
  ! WE ASSUME THAT THE SPACIAL COORDINATE OF THE NUCLEON
  ! IS THAT OF THE DELTA
  ppion(1, nnn, irun) = pxp + bdx*gd*(fgd*bpp+ep)
  ppion(2, nnn, irun) = pyp + bdy*gd*(fgd*bpp+ep)
  ppion(3, nnn, irun) = pzp + bdz*gd*(fgd*bpp+ep)
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i)
  ! WE ASSUME THE PION OR ETA COMING FROM DELTA DECAY IS LOCATED ON THE SPHERE
  ! OF RADIUS 0.5FM AROUND DELTA, THIS POINT NEED TO BE CHECKED
  ! AND OTHER CRIERTION MAY BE TRIED
  !lin-2/20/03 no additional smearing for position of decay daughters:
  !200         X0 = 1.0 - 2.0 * RANART(NSEED)
  !            Y0 = 1.0 - 2.0 * RANART(NSEED)
  !            Z0 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X0*X0+Y0*Y0+Z0*Z0) .GT. 1.0) GOTO 200
  !        RPION(1,NNN,IRUN)=R(1,I)+0.5*x0
  !        RPION(2,NNN,IRUN)=R(2,I)+0.5*y0
  !        RPION(3,NNN,IRUN)=R(3,I)+0.5*z0
  rpion(1, nnn, irun) = r(1, i)
  rpion(2, nnn, irun) = r(2, i)
  rpion(3, nnn, irun) = r(3, i)
  !
  devio = sqrt(epion(nnn,irun)**2+ppion(1,nnn,irun)**2+ppion(2,nnn,irun)**2+ppion(3,nnn,irun)**2) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2) - e1
  !        if(abs(devio).gt.0.02) write(93,*) 'decay(): nt=',nt,devio,lb1

  !     add decay time to daughter's formation time at the last timestep:
  If (nt==ntmax) Then
     tau0 = hbarc/wid
     taudcy = tau0*(-1.)*alog(1.-ranart(nseed))
     !     lorentz boost:
     taudcy = taudcy*e1/em1
     tfnl = tfnl + taudcy
     xfnl = xfnl + px1/e1*taudcy
     yfnl = yfnl + py1/e1*taudcy
     zfnl = zfnl + pz1/e1*taudcy
     r(1, i) = xfnl
     r(2, i) = yfnl
     r(3, i) = zfnl
     tfdcy(i) = tfnl
     rpion(1, nnn, irun) = xfnl
     rpion(2, nnn, irun) = yfnl
     rpion(3, nnn, irun) = zfnl
     tfdpi(nnn, irun) = tfnl
  End If

  Return

200 Format (A30, 2(1X,E10.4))
210 Format (I6, 5(1X,F8.3))
220 Format (A2, I5, 5(1X,F8.3))
End Subroutine dkine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! PURPOSE:1. N*-->N+PION+PION  DECAY PRODUCTS
!         2. DETERMINE THE MOMENTUM AND COORDINATES OF NUCLEON AND PION
!            AFTER THE DELTA OR N* DECAYING
! DATE   : NOV.7,1994
!----------------------------------------------------------------------------
Subroutine decay2(irun, i, nnn, iseed, wid, nt)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, etam=0.5475, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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

  lbanti = lb(i)
  !
  dm = e(i)
  ! DETERMINE THE DECAY PRODUCTS
  ! FOR N*+(1440) DECAY
  If (iabs(lb(i))==11) Then
     x3 = ranart(nseed)
     If (x3<(1./3)) Then
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     Else If (x3<2./3 .And. x3>1./3.) Then
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 3
        epion(nnn+1, irun) = ap2
     Else
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     End If
     ! FOR N*0(1440) DECAY
  Else If (iabs(lb(i))==10) Then
     x3 = ranart(nseed)
     If (x3<(1./3)) Then
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 4
        epion(nnn, irun) = ap1
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     Else If (x3<2./3 .And. x3>1./3.) Then
        lb(i) = 1
        nlab = 1
        lpion(nnn, irun) = 3
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 4
        epion(nnn+1, irun) = ap1
     Else
        lb(i) = 2
        nlab = 2
        lpion(nnn, irun) = 5
        epion(nnn, irun) = ap2
        lpion(nnn+1, irun) = 3
        epion(nnn+1, irun) = ap2
     End If
  End If

  Call dkine2(irun, i, nnn, nlab, iseed, wid, nt)
  !
  !     anti-particle ID for anti-N* decays:
  If (lbanti<0) Then
     lbi = lb(i)
     If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     Else If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     End If
     lb(i) = lbi
     !
     lbi = lpion(nnn, irun)
     If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     Else If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     End If
     lpion(nnn, irun) = lbi
     !
     lbi = lpion(nnn+1, irun)
     If (lbi==3) Then
        lbi = 5
     Else If (lbi==5) Then
        lbi = 3
     Else If (lbi==1 .Or. lbi==2) Then
        lbi = -lbi
     End If
     lpion(nnn+1, irun) = lbi
  End If
  !
  Return
End Subroutine decay2
!-------------------------------------------------------------------
!--------------------------------------------------------------------------
!         CALCULATE THE MOMENTUM OF NUCLEON AND PION (OR ETA)
!         IN THE LAB. FRAME AFTER DELTA OR N* DECAY
! DATE   : JAN. 24,1990, MODIFIED ON MAY 17, 1994 TO INCLUDE ETA PRODUCTION
!--------------------------------------------------------------------------
Subroutine dkine2(irun, i, nnn, nlab, iseed, wid, nt)
  Parameter (hbarc=0.19733)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, etam=0.5475, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  !c      SAVE /leadng/
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  !c      SAVE /tdecay/
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  !c      SAVE /INPUT2/
  External iarflv, invflv
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save

  ! READ IN THE COORDINATES OF THE N*(1440) UNDERGOING DECAY
  px = p(1, i)
  py = p(2, i)
  pz = p(3, i)
  rx = r(1, i)
  ry = r(2, i)
  rz = r(3, i)
  dm = e(i)
  edelta = sqrt(dm**2+px**2+py**2+pz**2)
  pm1 = epion(nnn, irun)
  pm2 = epion(nnn+1, irun)
  am = amn
  If (nlab==1) am = amp
  ! THE MAXIMUM MOMENTUM OF THE NUCLEON FROM THE DECAY OF A N*
  pmax2 = (dm**2-(am+pm1+pm2)**2)*(dm**2-(am-pm1-pm2)**2)/4/dm**2

  !lin-9/2012: check argument in sqrt():
  scheck = pmax2
  If (scheck<0) Then
     Write (99, *) 'scheck15: ', scheck
     scheck = 0.
  End If
  pmax = sqrt(scheck)
  !       PMAX=SQRT(PMAX2)

  ! GENERATE THE MOMENTUM OF THE NUCLEON IN THE N* REST FRAME
  css = 1. - 2.*ranart(nseed)
  sss = sqrt(1-css**2)
  fai = 2*pi*ranart(nseed)
  px0 = pmax*sss*cos(fai)
  py0 = pmax*sss*sin(fai)
  pz0 = pmax*css
  ep0 = sqrt(px0**2+py0**2+pz0**2+am**2)
  !lin-5/23/01 bug: P0 for pion0 is equal to PMAX, leaving pion+ and pion-
  !     without no relative momentum, thus producing them with equal momenta,
  ! BETA AND GAMMA OF THE CMS OF PION+-PION-
  betax = -px0/(dm-ep0)
  betay = -py0/(dm-ep0)
  betaz = -pz0/(dm-ep0)

  !lin-9/2012: check argument in sqrt():
  scheck = 1 - betax**2 - betay**2 - betaz**2
  If (scheck<=0) Then
     Write (99, *) 'scheck16: ', scheck
     Stop
  End If
  gd1 = 1./sqrt(scheck)
  !       GD1=1./SQRT(1-BETAX**2-BETAY**2-BETAZ**2)

  fgd1 = gd1/(1+gd1)
  ! GENERATE THE MOMENTA OF PIONS IN THE CMS OF PION+PION-
  q2 = ((dm-ep0)/(2.*gd1))**2 - pm1**2
  If (q2<=0.) q2 = 1.E-09
  q = sqrt(q2)
11 qx = 1. - 2.*ranart(nseed)
  qy = 1. - 2.*ranart(nseed)
  qz = 1. - 2.*ranart(nseed)
  qs = qx**2 + qy**2 + qz**2
  If (qs>1.) Goto 11
  pxp = q*qx/sqrt(qs)
  pyp = q*qy/sqrt(qs)
  pzp = q*qz/sqrt(qs)
  ep = sqrt(q**2+pm1**2)
  pxn = -pxp
  pyn = -pyp
  pzn = -pzp
  en = sqrt(q**2+pm2**2)
  ! TRANSFORM THE MOMENTA OF PION+PION- INTO THE N* REST FRAME
  bpp1 = betax*pxp + betay*pyp + betaz*pzp
  bpn1 = betax*pxn + betay*pyn + betaz*pzn
  ! FOR PION-
  p1m = pxn + betax*gd1*(fgd1*bpn1+en)
  p2m = pyn + betay*gd1*(fgd1*bpn1+en)
  p3m = pzn + betaz*gd1*(fgd1*bpn1+en)
  epn = sqrt(p1m**2+p2m**2+p3m**2+pm2**2)
  ! FOR PION+
  p1p = pxp + betax*gd1*(fgd1*bpp1+ep)
  p2p = pyp + betay*gd1*(fgd1*bpp1+ep)
  p3p = pzp + betaz*gd1*(fgd1*bpp1+ep)
  epp = sqrt(p1p**2+p2p**2+p3p**2+pm1**2)
  ! TRANSFORM MOMENTA OF THE THREE PIONS INTO THE
  ! THE NUCLEUS-NUCLEUS CENTER OF MASS  FRAME.
  ! THE GENERAL LORENTZ TRANSFORMATION CAN
  ! BE FOUND ON PAGE 34 OF R. HAGEDORN " RELATIVISTIC KINEMATICS"
  gd = edelta/dm
  fgd = gd/(1.+gd)
  bdx = px/edelta
  bdy = py/edelta
  bdz = pz/edelta
  bp0 = bdx*px0 + bdy*py0 + bdz*pz0
  bpp = bdx*p1p + bdy*p2p + bdz*p3p
  bpn = bdx*p1m + bdy*p2m + bdz*p3m
  ! FOR THE NUCLEON
  p(1, i) = px0 + bdx*gd*(fgd*bp0+ep0)
  p(2, i) = py0 + bdy*gd*(fgd*bp0+ep0)
  p(3, i) = pz0 + bdz*gd*(fgd*bp0+ep0)
  e(i) = am
  id(i) = 0
  enucl = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
  ! WE ASSUME THAT THE SPACIAL COORDINATE OF THE PION0
  ! IS in a sphere of radius 0.5 fm around N*
  ! FOR PION+
  ppion(1, nnn, irun) = p1p + bdx*gd*(fgd*bpp+epp)
  ppion(2, nnn, irun) = p2p + bdy*gd*(fgd*bpp+epp)
  ppion(3, nnn, irun) = p3p + bdz*gd*(fgd*bpp+epp)
  epion1 = sqrt(ppion(1,nnn,irun)**2+ppion(2,nnn,irun)**2+ppion(3,nnn,irun)**2+epion(nnn,irun)**2)
  !lin-2/20/03 no additional smearing for position of decay daughters:
  !200         X0 = 1.0 - 2.0 * RANART(NSEED)
  !            Y0 = 1.0 - 2.0 * RANART(NSEED)
  !            Z0 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X0*X0+Y0*Y0+Z0*Z0) .GT. 1.0) GOTO 200
  !        RPION(1,NNN,IRUN)=R(1,I)+0.5*x0
  !        RPION(2,NNN,IRUN)=R(2,I)+0.5*y0
  !        RPION(3,NNN,IRUN)=R(3,I)+0.5*z0
  rpion(1, nnn, irun) = r(1, i)
  rpion(2, nnn, irun) = r(2, i)
  rpion(3, nnn, irun) = r(3, i)
  ! FOR PION-
  ppion(1, nnn+1, irun) = p1m + bdx*gd*(fgd*bpn+epn)
  ppion(2, nnn+1, irun) = p2m + bdy*gd*(fgd*bpn+epn)
  ppion(3, nnn+1, irun) = p3m + bdz*gd*(fgd*bpn+epn)
  !lin-5/2008:
  dppion(nnn, irun) = dpertp(i)
  dppion(nnn+1, irun) = dpertp(i)
  !
  epion2 = sqrt(ppion(1,nnn+1,irun)**2+ppion(2,nnn+1,irun)**2+ppion(3,nnn+1,irun)**2+epion(nnn+1,irun)**2)
  !lin-2/20/03 no additional smearing for position of decay daughters:
  !300         X0 = 1.0 - 2.0 * RANART(NSEED)
  !            Y0 = 1.0 - 2.0 * RANART(NSEED)
  !            Z0 = 1.0 - 2.0 * RANART(NSEED)
  !        IF ((X0*X0+Y0*Y0+Z0*Z0) .GT. 1.0) GOTO 300
  !        RPION(1,NNN+1,IRUN)=R(1,I)+0.5*x0
  !        RPION(2,NNN+1,IRUN)=R(2,I)+0.5*y0
  !        RPION(3,NNN+1,IRUN)=R(3,I)+0.5*z0
  rpion(1, nnn+1, irun) = r(1, i)
  rpion(2, nnn+1, irun) = r(2, i)
  rpion(3, nnn+1, irun) = r(3, i)
  !
  ! check energy conservation in the decay
  !       efinal=enucl+epion1+epion2
  !       DEEE=(EDELTA-EFINAL)/EDELTA
  !       IF(ABS(DEEE).GE.1.E-03)write(6,*)1,edelta,efinal

  devio = sqrt(epion(nnn,irun)**2+ppion(1,nnn,irun)**2+ppion(2,nnn,irun)**2+ppion(3,nnn,irun)**2) + sqrt(e(i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2) + sqrt(epion(nnn+1,irun)**2+ppion(1,nnn+1,irun)**2+ppion(2,nnn+1,irun)**2+ppion(3,nnn+1,irun)**2) - e1
  !        if(abs(devio).gt.0.02) write(93,*) 'decay2(): nt=',nt,devio,lb1

  !     add decay time to daughter's formation time at the last timestep:
  If (nt==ntmax) Then
     tau0 = hbarc/wid
     taudcy = tau0*(-1.)*alog(1.-ranart(nseed))
     !     lorentz boost:
     taudcy = taudcy*e1/em1
     tfnl = tfnl + taudcy
     xfnl = xfnl + px1/e1*taudcy
     yfnl = yfnl + py1/e1*taudcy
     zfnl = zfnl + pz1/e1*taudcy
     r(1, i) = xfnl
     r(2, i) = yfnl
     r(3, i) = zfnl
     tfdcy(i) = tfnl
     rpion(1, nnn, irun) = xfnl
     rpion(2, nnn, irun) = yfnl
     rpion(3, nnn, irun) = zfnl
     tfdpi(nnn, irun) = tfnl
     rpion(1, nnn+1, irun) = xfnl
     rpion(2, nnn+1, irun) = yfnl
     rpion(3, nnn+1, irun) = zfnl
     tfdpi(nnn+1, irun) = tfnl
  End If

  Return

200 Format (A30, 2(1X,E10.4))
210 Format (I6, 5(1X,F8.3))
220 Format (A2, I5, 5(1X,F8.3))
End Subroutine dkine2
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! PURPOSE : CALCULATE THE MASS AND MOMENTUM OF BARYON RESONANCE
!           AFTER PION OR ETA BEING ABSORBED BY A NUCLEON
! NOTE    :
!
! DATE    : JAN.29,1990
Subroutine dreson(i1, i2)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  !lin-9/2012: improve precision for argument in sqrt():
  Double Precision e10, e20, scheck, p1, p2, p3
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Save
  ! 1. DETERMINE THE MOMENTUM COMPONENT OF DELTA/N* IN THE LAB. FRAME
  !lin-9/2012: improve precision for argument in sqrt():
  !        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
  !        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
  e10 = dsqrt(dble(e(i1))**2+dble(p(1,i1))**2+dble(p(2,i1))**2+dble(p(3,i1))**2)
  e20 = dsqrt(dble(e(i2))**2+dble(p(1,i2))**2+dble(p(2,i2))**2+dble(p(3,i2))**2)
  p1 = dble(p(1,i1)) + dble(p(1,i2))
  p2 = dble(p(2,i1)) + dble(p(2,i2))
  p3 = dble(p(3,i1)) + dble(p(3,i2))

  If (iabs(lb(i2))==1 .Or. iabs(lb(i2))==2 .Or. (iabs(lb(i2))>=6 .And. iabs(lb(i2))<=17)) Then
     e(i1) = 0.
     i = i2
  Else
     e(i2) = 0.
     i = i1
  End If
  p(1, i) = p(1, i1) + p(1, i2)
  p(2, i) = p(2, i1) + p(2, i2)
  p(3, i) = p(3, i1) + p(3, i2)
  ! 2. DETERMINE THE MASS OF DELTA/N* BY USING THE REACTION KINEMATICS

  !lin-9/2012: check argument in sqrt():
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
     Write (99, *) 'scheck17: ', scheck
     Write (99, *) 'scheck17', scheck, e10, e20, p(1, i), p(2, i), p(3, i)
     Write (99, *) 'scheck17-1', e(i1), p(1, i1), p(2, i1), p(3, i1)
     Write (99, *) 'scheck17-2', e(i2), p(1, i2), p(2, i2), p(3, i2)
     scheck = 0.D0
  End If
  dm = sqrt(sngl(scheck))
  !        DM=SQRT((E10+E20)**2-P(1,I)**2-P(2,I)**2-P(3,I)**2)

  e(i) = dm
  Return
End Subroutine dreson
!---------------------------------------------------------------------------
! PURPOSE : CALCULATE THE MASS AND MOMENTUM OF RHO RESONANCE
!           AFTER PION + PION COLLISION
! DATE    : NOV. 30,1994
Subroutine rhores(i1, i2)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  !lin-9/2012: improve precision for argument in sqrt():
  Double Precision e10, e20, scheck, p1, p2, p3
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Save
  ! 1. DETERMINE THE MOMENTUM COMPONENT OF THE RHO IN THE CMS OF NN FRAME
  !    WE LET I1 TO BE THE RHO AND ABSORB I2
  !lin-9/2012: improve precision for argument in sqrt():
  !        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
  !        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
  e10 = dsqrt(dble(e(i1))**2+dble(p(1,i1))**2+dble(p(2,i1))**2+dble(p(3,i1))**2)
  e20 = dsqrt(dble(e(i2))**2+dble(p(1,i2))**2+dble(p(2,i2))**2+dble(p(3,i2))**2)
  p1 = dble(p(1,i1)) + dble(p(1,i2))
  p2 = dble(p(2,i1)) + dble(p(2,i2))
  p3 = dble(p(3,i1)) + dble(p(3,i2))

  p(1, i1) = p(1, i1) + p(1, i2)
  p(2, i1) = p(2, i1) + p(2, i2)
  p(3, i1) = p(3, i1) + p(3, i2)
  ! 2. DETERMINE THE MASS OF THE RHO BY USING THE REACTION KINEMATICS

  !lin-9/2012: check argument in sqrt():
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
     Write (99, *) 'scheck18: ', scheck
     scheck = 0.D0
  End If
  dm = sqrt(sngl(scheck))
  !        DM=SQRT((E10+E20)**2-P(1,I1)**2-P(2,I1)**2-P(3,I1)**2)

  e(i1) = dm
  e(i2) = 0
  Return
End Subroutine rhores
!---------------------------------------------------------------------------
! PURPOSE : CALCULATE THE PION+NUCLEON CROSS SECTION ACCORDING TO THE
!           BREIT-WIGNER FORMULA/(p*)**2
! VARIABLE : LA = 1 FOR DELTA RESONANCE
!            LA = 0 FOR N*(1440) RESONANCE
!            LA = 2 FRO N*(1535) RESONANCE
! DATE    : JAN.29,1990
Real Function xnpi(i1, i2, la, xmax)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  !lin-9/2012: improve precision for argument in sqrt():
  Double Precision e10, e20, scheck, p1, p2, p3
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Save
  avmass = 0.5*(amn+amp)
  avpi = (2.*ap2+ap1)/3.
  ! 1. DETERMINE THE MOMENTUM COMPONENT OF DELTA IN THE LAB. FRAME
  !lin-9/2012: improve precision for argument in sqrt():
  !        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
  !        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
  e10 = dsqrt(dble(e(i1))**2+dble(p(1,i1))**2+dble(p(2,i1))**2+dble(p(3,i1))**2)
  e20 = dsqrt(dble(e(i2))**2+dble(p(1,i2))**2+dble(p(2,i2))**2+dble(p(3,i2))**2)
  !        P1=P(1,I1)+P(1,I2)
  !        P2=P(2,I1)+P(2,I2)
  !        P3=P(3,I1)+P(3,I2)
  p1 = dble(p(1,i1)) + dble(p(1,i2))
  p2 = dble(p(2,i1)) + dble(p(2,i2))
  p3 = dble(p(3,i1)) + dble(p(3,i2))

  ! 2. DETERMINE THE MASS OF DELTA BY USING OF THE REACTION KINEMATICS

  !lin-9/2012: check argument in sqrt():
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
     Write (99, *) 'scheck19: ', scheck
     scheck = 0.D0
  End If
  dm = sqrt(sngl(scheck))
  !        DM=SQRT((E10+E20)**2-P1**2-P2**2-P3**2)

  If (dm<=1.1) Then
     xnpi = 1.E-09
     Return
  End If
  ! 3. DETERMINE THE PION+NUCLEON CROSS SECTION ACCORDING TO THE
  !    BREIT-WIGNER FORMULA IN UNIT OF FM**2
  If (la==1) Then
     gam = width(dm)
     f1 = 0.25*gam**2/(0.25*gam**2+(dm-1.232)**2)
     pdelt2 = 0.051622
     Goto 10
  End If
  If (la==0) Then
     gam = w1440(dm)
     f1 = 0.25*gam**2/(0.25*gam**2+(dm-1.440)**2)
     pdelt2 = 0.157897
     Goto 10
  End If
  If (la==2) Then
     gam = w1535(dm)
     f1 = 0.25*gam**2/(0.25*gam**2+(dm-1.535)**2)
     pdelt2 = 0.2181
  End If
10 pstar2 = ((dm**2-avmass**2+avpi**2)/(2.*dm))**2 - avpi**2
  If (pstar2<=0.) Then
     xnpi = 1.E-09
  Else
     ! give the cross section in unit of fm**2
     xnpi = f1*(pdelt2/pstar2)*xmax/10.
  End If
  Return
End Function xnpi
!------------------------------------------------------------------------------
!****************************************
Real Function sigma(srt, id, ioi, iof)
  !PURPOSE : THIS IS THE PROGRAM TO CALCULATE THE ISOSPIN DECOMPOSED CROSS
  !       SECTION BY USING OF B.J.VerWEST AND R.A.ARNDT'S PARAMETERIZATION
  !REFERENCE: PHYS. REV. C25(1982)1979
  !QUANTITIES: IOI -- INITIAL ISOSPIN OF THE TWO NUCLEON SYSTEM
  !            IOF -- FINAL   ISOSPIN -------------------------
  !            ID -- =1 FOR DELTA RESORANCE
  !                  =2 FOR N*    RESORANCE
  !DATE : MAY 15,1990
  !****************************************
  Parameter (amu=0.9383, amp=0.1384, pi=3.1415926, hc=0.19733)
  Save
  If (id==1) Then
     amass0 = 1.22
     t0 = 0.12
  Else
     amass0 = 1.43
     t0 = 0.2
  End If
  If ((ioi==1) .And. (iof==1)) Then
     alfa = 3.772
     beta = 1.262
     am0 = 1.188
     t = 0.09902
  End If
  If ((ioi==1) .And. (iof==0)) Then
     alfa = 15.28
     beta = 0.
     am0 = 1.245
     t = 0.1374
  End If
  If ((ioi==0) .And. (iof==1)) Then
     alfa = 146.3
     beta = 0.
     am0 = 1.472
     t = 0.02649
  End If
  zplus = (srt-amu-amass0)*2./t0
  zminus = (amu+amp-amass0)*2./t0
  deln = atan(zplus) - atan(zminus)
  If (deln==0) deln = 1.E-06
  amass = amass0 + (t0/4.)*alog((1.+zplus**2)/(1.+zminus**2))/deln
  s = srt**2
  p2 = s/4. - amu**2
  s0 = (amu+am0)**2
  p02 = s0/4. - amu**2
  p0 = sqrt(p02)
  pr2 = (s-(amu-amass)**2)*(s-(amu+amass)**2)/(4.*s)
  If (pr2>1.E-06) Then
     pr = sqrt(pr2)
  Else
     pr = 0.
     sigma = 1.E-06
     Return
  End If
  ss = amass**2
  q2 = (ss-(amu-amp)**2)*(ss-(amu+amp)**2)/(4.*ss)
  If (q2>1.E-06) Then
     q = sqrt(q2)
  Else
     q = 0.
     sigma = 1.E-06
     Return
  End If
  ss0 = am0**2
  q02 = (ss0-(amu-amp)**2)*(ss0-(amu+amp)**2)/(4.*ss0)

  !lin-9/2012: check argument in sqrt():
  scheck = q02
  If (scheck<0) Then
     Write (99, *) 'scheck20: ', scheck
     scheck = 0.
  End If
  q0 = sqrt(scheck)
  !        Q0=SQRT(Q02)

  sigma = pi*(hc)**2/(2.*p2)*alfa*(pr/p0)**beta*am0**2*t**2*(q/q0)**3/((ss-am0**2)**2+am0**2*t**2)
  sigma = sigma*10.
  If (sigma==0) sigma = 1.E-06
  Return
End Function sigma

!****************************
Real Function denom(srt, con)
  ! NOTE: CON=1 FOR DELTA RESONANCE, CON=2 FOR N*(1440) RESONANCE
  !       con=-1 for N*(1535)
  ! PURPOSE : CALCULATE THE INTEGRAL IN THE DETAILED BALANCE
  !
  ! DATE : NOV. 15, 1991
  !******************************
  Parameter (ap1=0.13496, ap2=0.13957, pi=3.1415926, avmass=0.9383)
  Save
  avpi = (ap1+2.*ap2)/3.
  am0 = 1.232
  amn = avmass
  amp = avpi
  amax = srt - avmass
  amin = avmass + avpi
  nmax = 200
  dmass = (amax-amin)/float(nmax)
  sum = 0.
  Do i = 1, nmax + 1
     dm = amin + float(i-1)*dmass
     If (con==1.) Then
        q2 = ((dm**2-amn**2+amp**2)/(2.*dm))**2 - amp**2
        If (q2>0.) Then
           q = sqrt(q2)
        Else
           q = 1.E-06
        End If
        tq = 0.47*(q**3)/(amp**2*(1.+0.6*(q/amp)**2))
     Else If (con==2) Then
        tq = 0.2
        am0 = 1.44
     Else If (con==-1.) Then
        tq = 0.1
        am0 = 1.535
     End If
     a1 = 4.*tq*am0**2/(am0**2*tq**2+(dm**2-am0**2)**2)
     s = srt**2
     p0 = (s+dm**2-amn**2)**2/(4.*s) - dm**2
     If (p0<=0.) Then
        p1 = 1.E-06
     Else
        p1 = sqrt(p0)
     End If
     f = dm*a1*p1
     If ((i==1) .Or. (i==(nmax+1))) Then
        sum = sum + f*0.5
     Else
        sum = sum + f
     End If
  End Do
  denom = sum*dmass/(2.*pi)
  Return
End Function denom
!*********************************
! subroutine : ang.FOR
! PURPOSE : Calculate the angular distribution of Delta production process
! DATE    : Nov. 19, 1992
! REFERENCE: G. WOLF ET. AL., NUCL. PHYS. A517 (1990) 615
! Note: this function applies when srt is larger than 2.14 GeV,
! for less energetic reactions, we assume the angular distribution
! is isotropic.
!**********************************
Real Function ang(srt, iseed)
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  !        if(srt.le.2.14)then
  !       b1s=0.5
  !       b2s=0.
  !      endif
  If ((srt>2.14) .And. (srt<=2.4)) Then
     b1s = 29.03 - 23.75*srt + 4.865*srt**2
     b2s = -30.33 + 25.53*srt - 5.301*srt**2
  End If
  If (srt>2.4) Then
     b1s = 0.06
     b2s = 0.4
  End If
  x = ranart(nseed)
  p = b1s/b2s
  q = (2.*x-1.)*(b1s+b2s)/b2s
  If ((-q/2.+sqrt((q/2.)**2+(p/3.)**3))>=0.) Then
     ang1 = (-q/2.+sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  Else
     ang1 = -(q/2.-sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  End If
  If ((-q/2.-sqrt((q/2.)**2+(p/3.)**3)>=0.)) Then
     ang2 = (-q/2.-sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  Else
     ang2 = -(q/2.+sqrt((q/2.)**2+(p/3.)**3))**(1./3.)
  End If
  ang = ang1 + ang2
  Return
End Function ang
!--------------------------------------------------------------------------
!****subprogram * kaon production from pi+B collisions *******************
Real Function pnlka(srt)
  Save
  ! units: fm**2
  !**********************************C
  ala = 1.116
  aka = 0.498
  ana = 0.939
  t1 = ala + aka
  If (srt<=t1) Then
     pnlka = 0
  Else
     If (srt<1.7) sbbk = (0.9/0.091)*(srt-t1)
     If (srt>=1.7) sbbk = 0.09/(srt-1.6)
     pnlka = 0.25*sbbk
     ! give the cross section in units of fm**2
     pnlka = pnlka/10.
  End If
  Return
End Function pnlka
!-------------------------------------------------------------------------
!****subprogram * kaon production from pi+B collisions *******************
Real Function pnska(srt)
  Save
  !**********************************
  If (srt>3.0) Then
     pnska = 0
     Return
  End If
  ala = 1.116
  aka = 0.498
  ana = 0.939
  asa = 1.197
  t1 = asa + aka
  If (srt<=t1) Then
     pnska = 0
     Return
  End If
  If (srt<1.9) sbb1 = (0.7/0.218)*(srt-t1)
  If (srt>=1.9) sbb1 = 0.14/(srt-1.7)
  sbb2 = 0.
  If (srt>1.682) sbb2 = 0.5*(1.-0.75*(srt-1.682))
  pnska = 0.25*(sbb1+sbb2)
  ! give the cross section in fm**2
  pnska = pnska/10.
  Return
End Function pnska

!*******************************
!
!       Kaon momentum distribution in baryon-baryon-->N lamda K process
!
!       NOTE: dsima/dp is prototional to (1-p/p_max)(p/p_max)^2
!              we use rejection method to generate kaon momentum
!
!       Variables: Fkaon = F(p)/F_max
!                 srt   = cms energy of the colliding pair,
!                          used to calculate the P_max
!       Date: Feb. 8, 1994
!
!       Reference: C. M. Ko et al.
!*******************************
Real Function fkaon(p, pmax)
  Save
  fmax = 0.148
  If (pmax==0.) pmax = 0.000001
  fkaon = (1.-p/pmax)*(p/pmax)**2
  If (fkaon>fmax) fkaon = fmax
  fkaon = fkaon/fmax
  Return
End Function fkaon

!************************
! cross section for N*(1535) production in ND OR NN* collisions
! VARIABLES:
! LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES
! SRT IS THE CMS ENERGY
! X1535 IS THE N*(1535) PRODUCTION CROSS SECTION
! NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA
! PRODUCTION CROSS SECTION
! DATE: MAY 18, 1994
! ***********************
Subroutine m1535(lb1, lb2, srt, x1535)
  Save
  s0 = 2.424
  x1535 = 0.
  If (srt<=s0) Return
  sigma = 2.*0.102*(srt-s0)/(0.058+(srt-s0)**2)
  ! I N*(1535) PRODUCTION IN NUCLEON-DELTA COLLISIONS
  !(1) nD(++)->pN*(+)(1535), pD(-)->nN*(0)(1535),pD(+)-->N*(+)p
  !bz11/25/98
  !       IF((LB1*LB2.EQ.18).OR.(LB1*LB2.EQ.6).
  !     1  or.(lb1*lb2).eq.8)then
  If ((lb1*lb2==18 .And. (lb1==2 .Or. lb2==2)) .Or. (lb1*lb2==6 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==8 .And. (lb1==1 .Or. lb2==1))) Then
     !bz11/25/98end
     x1535 = sigma
     Return
  End If
  !(2) pD(0)->pN*(0)(1535),pD(0)->nN*(+)(1535)
  If (lb1*lb2==7) Then
     x1535 = 3.*sigma
     Return
  End If
  ! II N*(1535) PRODUCTION IN N*(1440)+NUCLEON REACTIONS
  !(3) N*(+)(1440)p->N*(0+)(1535)p, N*(0)(1440)n->N*(0)(1535)
  !bz11/25/98
  !       IF((LB1*LB2.EQ.11).OR.(LB1*LB2.EQ.20))THEN
  If ((lb1*lb2==11) .Or. (lb1*lb2==20 .And. (lb1==2 .Or. lb2==2))) Then
     !bz11/25/98end
     x1535 = sigma
     Return
  End If
  !(4) N*(0)(1440)p->N*(0+) or N*(+)(1440)n->N*(0+)(1535)
  !bz11/25/98
  !       IF((LB1*LB2.EQ.10).OR.(LB1*LB2.EQ.22))X1535=3.*SIGMA
  If ((lb1*lb2==10 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==22 .And. (lb1==2 .Or. lb2==2))) x1535 = 3.*sigma
  !bz11/25/98end
  Return
End Subroutine m1535
!************************
! cross section for N*(1535) production in NN collisions
! VARIABLES:
! LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES
! SRT IS THE CMS ENERGY
! X1535 IS THE N*(1535) PRODUCTION CROSS SECTION
! NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA
! PRODUCTION CROSS SECTION
! DATE: MAY 18, 1994
! ***********************
Subroutine n1535(lb1, lb2, srt, x1535)
  Save
  s0 = 2.424
  x1535 = 0.
  If (srt<=s0) Return
  sigma = 2.*0.102*(srt-s0)/(0.058+(srt-s0)**2)
  ! I N*(1535) PRODUCTION IN NUCLEON-NUCLEON COLLISIONS
  !(1) pp->pN*(+)(1535), nn->nN*(0)(1535)
  !bdbg11/25/98
  !       IF((LB1*LB2.EQ.1).OR.(LB1*LB2.EQ.4))then
  If ((lb1*lb2==1) .Or. (lb1==2 .And. lb2==2)) Then
     !bz11/25/98end
     x1535 = sigma
     Return
  End If
  !(2) pn->pN*(0)(1535),pn->nN*(+)(1535)
  If (lb1*lb2==2) Then
     x1535 = 3.*sigma
     Return
  End If
  ! III N*(1535) PRODUCTION IN DELTA+DELTA REACTIONS
  ! (5) D(++)+D(0), D(+)+D(+),D(+)+D(-),D(0)+D(0)
  !bz11/25/98
  !       IF((LB1*LB2.EQ.63).OR.(LB1*LB2.EQ.64).OR.(LB1*LB2.EQ.48).
  !     1  OR.(LB1*LB2.EQ.49))then
  If ((lb1*lb2==63 .And. (lb1==7 .Or. lb2==7)) .Or. (lb1*lb2==64 .And. (lb1==8 .Or. lb2==8)) .Or. (lb1*lb2==48 .And. (lb1==6 .Or. lb2==6)) .Or. (lb1*lb2==49 .And. (lb1==7 .Or. lb2==7))) Then
     !bz11/25/98end
     x1535 = sigma
     Return
  End If
  ! (6) D(++)+D(-),D(+)+D(0)
  !bz11/25/98
  !       IF((LB1*LB2.EQ.54).OR.(LB1*LB2.EQ.56))then
  If ((lb1*lb2==54 .And. (lb1==6 .Or. lb2==6)) .Or. (lb1*lb2==56 .And. (lb1==7 .Or. lb2==7))) Then
     !bz11/25/98end
     x1535 = 3.*sigma
     Return
  End If
  ! IV N*(1535) PRODUCTION IN N*(1440)+N*(1440) REACTIONS
  !bz11/25/98
  !       IF((LB1*LB2.EQ.100).OR.(LB1*LB2.EQ.11*11))X1535=SIGMA
  If ((lb1==10 .And. lb2==10) .Or. (lb1==11 .And. lb2==11)) x1535 = sigma
  !       IF(LB1*LB2.EQ.110)X1535=3.*SIGMA
  If (lb1*lb2==110 .And. (lb1==10 .Or. lb2==10)) x1535 = 3.*sigma
  !bdbg11/25/98end
  Return
End Subroutine n1535
!***********************************
! FUNCTION WA1(DMASS) GIVES THE A1 DECAY WIDTH

Subroutine wida1(dmass, rhomp, wa1, iseed)
  Save
  !
  pimass = 0.137265
  coupa = 14.8
  !
  rhomax = dmass - pimass - 0.02
  If (rhomax<=0) Then
     rhomp = 0.
     !   !! no decay
     wa1 = -10.
  End If
  icount = 0
711 rhomp = rhomas(rhomax, iseed)
  icount = icount + 1
  If (dmass<=(pimass+rhomp)) Then
     If (icount<=100) Then
        Goto 711
     Else
        rhomp = 0.
        !   !! no decay
        wa1 = -10.
        Return
     End If
  End If
  qqp2 = (dmass**2-(rhomp+pimass)**2)*(dmass**2-(rhomp-pimass)**2)
  qqp = sqrt(qqp2)/(2.0*dmass)
  epi = sqrt(pimass**2+qqp**2)
  erho = sqrt(rhomp**2+qqp**2)
  epirho = 2.0*(epi*erho+qqp**2)**2 + rhomp**2*epi**2
  wa1 = coupa**2*qqp*epirho/(24.0*3.1416*dmass**2)
  Return
End Subroutine wida1
!***********************************
! FUNCTION W1535(DMASS) GIVES THE N*(1535) DECAY WIDTH
!     FOR A GIVEN N*(1535) MASS
! HERE THE FORMULA GIVEN BY KITAZOE IS USED
Real Function w1535(dmass)
  Save
  avmass = 0.938868
  pimass = 0.137265
  aux = 0.25*(dmass**2-avmass**2-pimass**2)**2 - (avmass*pimass)**2
  If (aux>0.) Then
     qavail = sqrt(aux/dmass**2)
  Else
     qavail = 1.E-06
  End If
  w1535 = 0.15*qavail/0.467
  !       W1535=0.15
  Return
End Function w1535
!***********************************
! FUNCTION W1440(DMASS) GIVES THE N*(1440) DECAY WIDTH
!     FOR A GIVEN N*(1535) MASS
! HERE THE FORMULA GIVEN BY KITAZOE IS USED
Real Function w1440(dmass)
  Save
  avmass = 0.938868
  pimass = 0.137265
  aux = 0.25*(dmass**2-avmass**2-pimass**2)**2 - (avmass*pimass)**2
  If (aux>0.) Then
     qavail = sqrt(aux)/dmass
  Else
     qavail = 1.E-06
  End If
  !              w1440=0.2
  w1440 = 0.2*(qavail/0.397)**3
  Return
End Function w1440
!***************
! PURPOSE : CALCULATE THE PION(ETA)+NUCLEON CROSS SECTION
!           ACCORDING TO THE BREIT-WIGNER FORMULA,
!           NOTE THAT N*(1535) IS S_11
! VARIABLE : LA = 1 FOR PI+N
!            LA = 0 FOR ETA+N
! DATE    : MAY 16, 1994
!***************
Real Function xn1535(i1, i2, la)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, etam=0.5475, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  !lin-9/2012: improve precision for argument in sqrt():
  Double Precision e10, e20, scheck, p1, p2, p3
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Save
  avmass = 0.5*(amn+amp)
  avpi = (2.*ap2+ap1)/3.
  ! 1. DETERMINE THE MOMENTUM COMPONENT OF N*(1535) IN THE LAB. FRAME
  !lin-9/2012: improve precision for argument in sqrt():
  !        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
  !        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
  e10 = dsqrt(dble(e(i1))**2+dble(p(1,i1))**2+dble(p(2,i1))**2+dble(p(3,i1))**2)
  e20 = dsqrt(dble(e(i2))**2+dble(p(1,i2))**2+dble(p(2,i2))**2+dble(p(3,i2))**2)
  !        P1=P(1,I1)+P(1,I2)
  !        P2=P(2,I1)+P(2,I2)
  !        P3=P(3,I1)+P(3,I2)
  p1 = dble(p(1,i1)) + dble(p(1,i2))
  p2 = dble(p(2,i1)) + dble(p(2,i2))
  p3 = dble(p(3,i1)) + dble(p(3,i2))

  ! 2. DETERMINE THE MASS OF DELTA BY USING OF THE REACTION KINEMATICS

  !lin-9/2012: check argument in sqrt():
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
     Write (99, *) 'scheck21: ', scheck
     scheck = 0.D0
  End If
  dm = sqrt(sngl(scheck))
  !        DM=SQRT((E10+E20)**2-P1**2-P2**2-P3**2)

  If (dm<=1.1) Then
     xn1535 = 1.E-06
     Return
  End If
  ! 3. DETERMINE THE PION(ETA)+NUCLEON->N*(1535) CROSS SECTION ACCORDING TO THE
  !    BREIT-WIGNER FORMULA IN UNIT OF FM**2
  gam = w1535(dm)
  gam0 = 0.15
  f1 = 0.25*gam0**2/(0.25*gam**2+(dm-1.535)**2)
  If (la==1) Then
     xmax = 11.3
  Else
     xmax = 74.
  End If
  xn1535 = f1*xmax/10.
  Return
End Function xn1535
!**************************8
!FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF
!KITAZOE'S FORMULA
Real Function fdelta(dmass)
  Save
  amn = 0.938869
  avpi = 0.13803333
  am0 = 1.232
  fd = 0.25*width(dmass)**2/((dmass-1.232)**2+0.25*width(dmass)**2)
  fdelta = fd
  Return
End Function fdelta
! FUNCTION WIDTH(DMASS) GIVES THE DELTA DECAY WIDTH FOR A GIVEN DELTA MASS
! HERE THE FORMULA GIVEN BY KITAZOE IS USED
Real Function width(dmass)
  Save
  avmass = 0.938868
  pimass = 0.137265
  aux = 0.25*(dmass**2-avmass**2-pimass**2)**2 - (avmass*pimass)**2
  If (aux>0.) Then
     qavail = sqrt(aux/dmass**2)
  Else
     qavail = 1.E-06
  End If
  width = 0.47*qavail**3/(pimass**2*(1.+0.6*(qavail/pimass)**2))
  !       width=0.115
  Return
End Function width
!***********************************
Subroutine ddp2(srt, iseed, px, py, pz, dm1, pnx, pny, pnz, dm2, ppx, ppy, ppz, icou1)
  ! PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM
  ! THE PROCESS N+N--->D1+D2+PION
  !       DATE : July 25, 1994
  ! Generate the masses and momentum for particles in the NN-->DDpi process
  ! for a given center of mass energy srt, the momenta are given in the center
  ! of mass of the NN
  !****************************************
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  icou1 = 0
  pi = 3.1415926
  amn = 938.925/1000.
  amp = 137.265/1000.
  ! (1) GENGRATE THE MASS OF DELTA1 AND DELTA2 USING
  srt1 = srt - amp - 0.02
  ntrym = 0
8 Call rmasdd(srt1, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
  ntrym = ntrym + 1
  ! CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM
  ! FOR ONE OF THE RESONANCES
  v = 0.43
  w = -0.84
  ! (2) Generate the transverse momentum
  !     OF DELTA1
  ! (2.1) estimate the maximum transverse momentum
  ptmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2
  If (ptmax2<=0) Goto 8
  ptmax = sqrt(ptmax2)*1./3.
7 pt = ptr(ptmax, iseed)
  ! (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS
  pzmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2 - pt**2
  If ((pzmax2<0.) .And. ntrym<=100) Then
     Goto 7
  Else
     pzmax2 = 1.E-09
  End If
  pzmax = sqrt(pzmax2)
  xmax = 2.*pzmax/srt
  ! (3.2) THE GENERATED X IS
  ! THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056
  ntryx = 0
  fmax00 = 1.056
  x00 = 0.26
  If (abs(xmax)>0.26) Then
     f00 = fmax00
  Else
     f00 = 1. + v*abs(xmax) + w*xmax**2
  End If
9 x = xmax*(1.-2.*ranart(nseed))
  ntryx = ntryx + 1
  xratio = (1.+v*abs(x)+w*x**2)/f00
  !lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9
  If (xratio<ranart(nseed) .And. ntryx<=50) Goto 9
  ! (3.5) THE PZ IS
  pz = 0.5*srt*x
  ! The x and y components of the deltA1
  fai = 2.*pi*ranart(nseed)
  px = pt*cos(fai)
  py = pt*sin(fai)
  ! find the momentum of delta2 and pion
  ! the energy of the delta1
  ek = sqrt(dm1**2+pt**2+pz**2)
  ! (1) Generate the momentum of the delta2 in the cms of delta2 and pion
  !     the energy of the cms of DP
  eln = srt - ek
  If (eln<=0) Then
     icou1 = -1
     Return
  End If
  ! beta and gamma of the cms of delta2+pion
  bx = -px/eln
  by = -py/eln
  bz = -pz/eln

  !lin-9/2012: check argument in sqrt():
  scheck = 1. - bx**2 - by**2 - bz**2
  If (scheck<=0) Then
     Write (99, *) 'scheck22: ', scheck
     Stop
  End If
  ga = 1./sqrt(scheck)
  !       ga=1./sqrt(1.-bx**2-by**2-bz**2)

  ! the momentum of delta2 and pion in their cms frame
  elnc = eln/ga
  pn2 = ((elnc**2+dm2**2-amp**2)/(2.*elnc))**2 - dm2**2
  If (pn2<=0) Then
     icou1 = -1
     Return
  End If
  pn = sqrt(pn2)

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pn
  !       PNT=PTR(0.33*PN,ISEED)
  pnt = ptr(xptr, iseed)
  !lin-10/25/02-end

  fain = 2.*pi*ranart(nseed)
  pnx = pnt*cos(fain)
  pny = pnt*sin(fain)
  sig = 1
  If (x>0) sig = -1

  !lin-9/2012: check argument in sqrt():
  scheck = pn**2 - pnt**2
  If (scheck<0) Then
     Write (99, *) 'scheck23: ', scheck
     scheck = 0.
  End If
  pnz = sig*sqrt(scheck)
  !       pnz=SIG*SQRT(pn**2-PNT**2)

  en = sqrt(dm2**2+pnx**2+pny**2+pnz**2)
  ! (2) the momentum for the pion
  ppx = -pnx
  ppy = -pny
  ppz = -pnz
  ep = sqrt(amp**2+ppx**2+ppy**2+ppz**2)
  ! (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  pbeta = pnx*bx + pny*by + pnz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+en)
  pnx = bx*trans0 + pnx
  pny = by*trans0 + pny
  pnz = bz*trans0 + pnz
  ! (4) for the pion, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  If (ep==0.) ep = 1.E-09
  pbeta = ppx*bx + ppy*by + ppz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+ep)
  ppx = bx*trans0 + ppx
  ppy = by*trans0 + ppy
  ppz = bz*trans0 + ppz
  Return
End Subroutine ddp2
!***************************************
Subroutine ddrho(srt, iseed, px, py, pz, dm1, pnx, pny, pnz, dm2, ppx, ppy, ppz, amp, icou1)
  ! PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM
  ! THE PROCESS N+N--->D1+D2+rho
  !       DATE : Nov.5, 1994
  ! Generate the masses and momentum for particles in the NN-->DDrho process
  ! for a given center of mass energy srt, the momenta are given in the center
  ! of mass of the NN
  !****************************************
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  icou1 = 0
  pi = 3.1415926
  amn = 938.925/1000.
  amp = 770./1000.
  ! (1) GENGRATE THE MASS OF DELTA1 AND DELTA2 USING
  srt1 = srt - amp - 0.02
  ntrym = 0
8 Call rmasdd(srt1, 1.232, 1.232, 1.08, 1.08, iseed, 1, dm1, dm2)
  ntrym = ntrym + 1
  ! GENERATE THE MASS FOR THE RHO
  rhomax = srt - dm1 - dm2 - 0.02
  If (rhomax<=0 .And. ntrym<=20) Goto 8
  amp = rhomas(rhomax, iseed)
  ! CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM
  ! FOR ONE OF THE RESONANCES
  v = 0.43
  w = -0.84
  ! (2) Generate the transverse momentum
  !     OF DELTA1
  ! (2.1) estimate the maximum transverse momentum
  ptmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2

  !lin-9/2012: check argument in sqrt():
  scheck = ptmax2
  If (scheck<0) Then
     Write (99, *) 'scheck24: ', scheck
     scheck = 0.
  End If
  ptmax = sqrt(scheck)*1./3.
  !       PTMAX=SQRT(PTMAX2)*1./3.

7 pt = ptr(ptmax, iseed)
  ! (3) GENGRATE THE LONGITUDINAL MOMENTUM FOR DM1
  !     USING THE GIVEN DISTRIBUTION
  ! (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS
  pzmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2 - pt**2
  If ((pzmax2<0.) .And. ntrym<=100) Then
     Goto 7
  Else
     pzmax2 = 1.E-06
  End If
  pzmax = sqrt(pzmax2)
  xmax = 2.*pzmax/srt
  ! (3.2) THE GENERATED X IS
  ! THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056
  ntryx = 0
  fmax00 = 1.056
  x00 = 0.26
  If (abs(xmax)>0.26) Then
     f00 = fmax00
  Else
     f00 = 1. + v*abs(xmax) + w*xmax**2
  End If
9 x = xmax*(1.-2.*ranart(nseed))
  ntryx = ntryx + 1
  xratio = (1.+v*abs(x)+w*x**2)/f00
  !lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9
  If (xratio<ranart(nseed) .And. ntryx<=50) Goto 9
  ! (3.5) THE PZ IS
  pz = 0.5*srt*x
  ! The x and y components of the delta1
  fai = 2.*pi*ranart(nseed)
  px = pt*cos(fai)
  py = pt*sin(fai)
  ! find the momentum of delta2 and rho
  ! the energy of the delta1
  ek = sqrt(dm1**2+pt**2+pz**2)
  ! (1) Generate the momentum of the delta2 in the cms of delta2 and rho
  !     the energy of the cms of Drho
  eln = srt - ek
  If (eln<=0) Then
     icou1 = -1
     Return
  End If
  ! beta and gamma of the cms of delta2 and rho
  bx = -px/eln
  by = -py/eln
  bz = -pz/eln

  !lin-9/2012: check argument in sqrt():
  scheck = 1. - bx**2 - by**2 - bz**2
  If (scheck<=0) Then
     Write (99, *) 'scheck25: ', scheck
     Stop
  End If
  ga = 1./sqrt(scheck)
  !       ga=1./sqrt(1.-bx**2-by**2-bz**2)

  elnc = eln/ga
  pn2 = ((elnc**2+dm2**2-amp**2)/(2.*elnc))**2 - dm2**2
  If (pn2<=0) Then
     icou1 = -1
     Return
  End If
  pn = sqrt(pn2)

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pn
  !       PNT=PTR(0.33*PN,ISEED)
  pnt = ptr(xptr, iseed)
  !lin-10/25/02-end

  fain = 2.*pi*ranart(nseed)
  pnx = pnt*cos(fain)
  pny = pnt*sin(fain)
  sig = 1
  If (x>0) sig = -1

  !lin-9/2012: check argument in sqrt():
  scheck = pn**2 - pnt**2
  If (scheck<0) Then
     Write (99, *) 'scheck26: ', scheck
     scheck = 0.
  End If
  pnz = sig*sqrt(scheck)
  !       pnz=SIG*SQRT(pn**2-PNT**2)

  en = sqrt(dm2**2+pnx**2+pny**2+pnz**2)
  ! (2) the momentum for the rho
  ppx = -pnx
  ppy = -pny
  ppz = -pnz
  ep = sqrt(amp**2+ppx**2+ppy**2+ppz**2)
  ! (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  pbeta = pnx*bx + pny*by + pnz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+en)
  pnx = bx*trans0 + pnx
  pny = by*trans0 + pny
  pnz = bz*trans0 + pnz
  ! (4) for the rho, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  If (ep==0.) ep = 1.E-09
  pbeta = ppx*bx + ppy*by + ppz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+ep)
  ppx = bx*trans0 + ppx
  ppy = by*trans0 + ppy
  ppz = bz*trans0 + ppz
  Return
End Subroutine ddrho
!***************************************
Subroutine pprho(srt, iseed, px, py, pz, dm1, pnx, pny, pnz, dm2, ppx, ppy, ppz, amp, icou1)
  ! PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM
  ! THE PROCESS N+N--->N1+N2+rho
  !       DATE : Nov.5, 1994
  ! Generate the masses and momentum for particles in the NN--> process
  ! for a given center of mass energy srt, the momenta are given in the center
  ! of mass of the NN
  !****************************************
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  ntrym = 0
  icou1 = 0
  pi = 3.1415926
  amn = 938.925/1000.
  !        AMP=770./1000.
  dm1 = amn
  dm2 = amn
  ! GENERATE THE MASS FOR THE RHO
  rhomax = srt - dm1 - dm2 - 0.02
  If (rhomax<=0) Then
     icou = -1
     Return
  End If
  amp = rhomas(rhomax, iseed)
  ! CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM
  ! FOR ONE OF THE nucleons
  v = 0.43
  w = -0.84
  ! (2) Generate the transverse momentum
  !     OF p1
  ! (2.1) estimate the maximum transverse momentum
  ptmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2

  !lin-9/2012: check argument in sqrt():
  scheck = ptmax2
  If (scheck<0) Then
     Write (99, *) 'scheck27: ', scheck
     scheck = 0.
  End If
  ptmax = sqrt(scheck)*1./3.
  !       PTMAX=SQRT(PTMAX2)*1./3.

7 pt = ptr(ptmax, iseed)
  ! (3) GENGRATE THE LONGITUDINAL MOMENTUM FOR DM1
  !     USING THE GIVEN DISTRIBUTION
  ! (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS
  pzmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2 - pt**2
  ntrym = ntrym + 1
  If ((pzmax2<0.) .And. ntrym<=100) Then
     Goto 7
  Else
     pzmax2 = 1.E-06
  End If
  pzmax = sqrt(pzmax2)
  xmax = 2.*pzmax/srt
  ! (3.2) THE GENERATED X IS
  ! THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056
  ntryx = 0
  fmax00 = 1.056
  x00 = 0.26
  If (abs(xmax)>0.26) Then
     f00 = fmax00
  Else
     f00 = 1. + v*abs(xmax) + w*xmax**2
  End If
9 x = xmax*(1.-2.*ranart(nseed))
  ntryx = ntryx + 1
  xratio = (1.+v*abs(x)+w*x**2)/f00
  !lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9
  If (xratio<ranart(nseed) .And. ntryx<=50) Goto 9
  ! (3.5) THE PZ IS
  pz = 0.5*srt*x
  ! The x and y components of the delta1
  fai = 2.*pi*ranart(nseed)
  px = pt*cos(fai)
  py = pt*sin(fai)
  ! find the momentum of delta2 and rho
  ! the energy of the delta1
  ek = sqrt(dm1**2+pt**2+pz**2)
  ! (1) Generate the momentum of the delta2 in the cms of delta2 and rho
  !     the energy of the cms of Drho
  eln = srt - ek
  If (eln<=0) Then
     icou1 = -1
     Return
  End If
  ! beta and gamma of the cms of the two partciles
  bx = -px/eln
  by = -py/eln
  bz = -pz/eln

  !lin-9/2012: check argument in sqrt():
  scheck = 1. - bx**2 - by**2 - bz**2
  If (scheck<=0) Then
     Write (99, *) 'scheck28: ', scheck
     Stop
  End If
  ga = 1./sqrt(scheck)
  !       ga=1./sqrt(1.-bx**2-by**2-bz**2)

  elnc = eln/ga
  pn2 = ((elnc**2+dm2**2-amp**2)/(2.*elnc))**2 - dm2**2
  If (pn2<=0) Then
     icou1 = -1
     Return
  End If
  pn = sqrt(pn2)

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pn
  !       PNT=PTR(0.33*PN,ISEED)
  pnt = ptr(xptr, iseed)
  !lin-10/25/02-end

  fain = 2.*pi*ranart(nseed)
  pnx = pnt*cos(fain)
  pny = pnt*sin(fain)
  sig = 1
  If (x>0) sig = -1

  !lin-9/2012: check argument in sqrt():
  scheck = pn**2 - pnt**2
  If (scheck<0) Then
     Write (99, *) 'scheck29: ', scheck
     scheck = 0.
  End If
  pnz = sig*sqrt(scheck)
  !       pnz=SIG*SQRT(pn**2-PNT**2)

  en = sqrt(dm2**2+pnx**2+pny**2+pnz**2)
  ! (2) the momentum for the rho
  ppx = -pnx
  ppy = -pny
  ppz = -pnz
  ep = sqrt(amp**2+ppx**2+ppy**2+ppz**2)
  ! (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  pbeta = pnx*bx + pny*by + pnz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+en)
  pnx = bx*trans0 + pnx
  pny = by*trans0 + pny
  pnz = bz*trans0 + pnz
  ! (4) for the rho, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  If (ep==0.) ep = 1.E-09
  pbeta = ppx*bx + ppy*by + ppz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+ep)
  ppx = bx*trans0 + ppx
  ppy = by*trans0 + ppy
  ppz = bz*trans0 + ppz
  Return
End Subroutine pprho
!**************************8
!***************************************

!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine ppomga(srt, iseed, px, py, pz, dm1, pnx, pny, pnz, dm2, ppx, ppy, ppz, icou1)
  ! PURPOSE : CALCULATE MOMENTUM OF PARTICLES IN THE FINAL SATAT FROM
  ! THE PROCESS N+N--->N1+N2+OMEGA
  !       DATE : Nov.5, 1994
  ! Generate the masses and momentum for particles in the NN--> process
  ! for a given center of mass energy srt, the momenta are given in the center
  ! of mass of the NN
  !****************************************
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  ntrym = 0
  icou1 = 0
  pi = 3.1415926
  amn = 938.925/1000.
  amp = 782./1000.
  dm1 = amn
  dm2 = amn
  ! CONSTANTS FOR GENERATING THE LONGITUDINAL MOMENTUM
  ! FOR ONE OF THE nucleons
  v = 0.43
  w = -0.84
  ! (2) Generate the transverse momentum
  !     OF p1
  ! (2.1) estimate the maximum transverse momentum
  ptmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2

  !lin-9/2012: check argument in sqrt():
  scheck = ptmax2
  If (scheck<0) Then
     Write (99, *) 'scheck30: ', scheck
     scheck = 0.
  End If
  ptmax = sqrt(scheck)*1./3.
  !       PTMAX=SQRT(PTMAX2)*1./3.

7 pt = ptr(ptmax, iseed)
  ! (3) GENGRATE THE LONGITUDINAL MOMENTUM FOR DM1
  !     USING THE GIVEN DISTRIBUTION
  ! (3.1) THE MAXIMUM LONGITUDINAL MOMENTUM IS
  pzmax2 = (srt**2-(dm1+dm2+amp)**2)*(srt**2-(dm1-amp-dm2)**2)/4./srt**2 - pt**2
  ntrym = ntrym + 1
  If ((pzmax2<0.) .And. ntrym<=100) Then
     Goto 7
  Else
     pzmax2 = 1.E-09
  End If
  pzmax = sqrt(pzmax2)
  xmax = 2.*pzmax/srt
  ! (3.2) THE GENERATED X IS
  ! THE DSTRIBUTION HAS A MAXIMUM AT X0=-V/(2*w), f(X0)=1.056
  ntryx = 0
  fmax00 = 1.056
  x00 = 0.26
  If (abs(xmax)>0.26) Then
     f00 = fmax00
  Else
     f00 = 1. + v*abs(xmax) + w*xmax**2
  End If
9 x = xmax*(1.-2.*ranart(nseed))
  ntryx = ntryx + 1
  xratio = (1.+v*abs(x)+w*x**2)/f00
  !lin-8/17/00       IF(xratio.LT.RANART(NSEED).and.ntryx.le.50)GO TO 9
  If (xratio<ranart(nseed) .And. ntryx<=50) Goto 9
  ! (3.5) THE PZ IS
  pz = 0.5*srt*x
  ! The x and y components of the delta1
  fai = 2.*pi*ranart(nseed)
  px = pt*cos(fai)
  py = pt*sin(fai)
  ! find the momentum of delta2 and rho
  ! the energy of the delta1
  ek = sqrt(dm1**2+pt**2+pz**2)
  ! (1) Generate the momentum of the delta2 in the cms of delta2 and rho
  !     the energy of the cms of Drho
  eln = srt - ek
  If (eln<=0) Then
     icou1 = -1
     Return
  End If
  bx = -px/eln
  by = -py/eln
  bz = -pz/eln

  !lin-9/2012: check argument in sqrt():
  scheck = 1. - bx**2 - by**2 - bz**2
  If (scheck<=0) Then
     Write (99, *) 'scheck31: ', scheck
     Stop
  End If
  ga = 1./sqrt(scheck)
  !       ga=1./sqrt(1.-bx**2-by**2-bz**2)

  elnc = eln/ga
  pn2 = ((elnc**2+dm2**2-amp**2)/(2.*elnc))**2 - dm2**2
  If (pn2<=0) Then
     icou1 = -1
     Return
  End If
  pn = sqrt(pn2)

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pn
  !       PNT=PTR(0.33*PN,ISEED)
  pnt = ptr(xptr, iseed)
  !lin-10/25/02-end

  fain = 2.*pi*ranart(nseed)
  pnx = pnt*cos(fain)
  pny = pnt*sin(fain)
  sig = 1
  If (x>0) sig = -1

  !lin-9/2012: check argument in sqrt():
  scheck = pn**2 - pnt**2
  If (scheck<0) Then
     Write (99, *) 'scheck32: ', scheck
     scheck = 0.
  End If
  pnz = sig*sqrt(scheck)
  !       pnz=SIG*SQRT(pn**2-PNT**2)

  en = sqrt(dm2**2+pnx**2+pny**2+pnz**2)
  ! (2) the momentum for the rho
  ppx = -pnx
  ppy = -pny
  ppz = -pnz
  ep = sqrt(amp**2+ppx**2+ppy**2+ppz**2)
  ! (3) for the delta2, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  pbeta = pnx*bx + pny*by + pnz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+en)
  pnx = bx*trans0 + pnx
  pny = by*trans0 + pny
  pnz = bz*trans0 + pnz
  ! (4) for the rho, LORENTZ-TRANSFORMATION INTO nn cms FRAME
  If (ep==0.) ep = 1.E-09
  pbeta = ppx*bx + ppy*by + ppz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+ep)
  ppx = bx*trans0 + ppx
  ppy = by*trans0 + ppy
  ppz = bz*trans0 + ppz
  Return
End Subroutine ppomga
!**************************8
!**************************8
!   DELTA MASS GENERATOR
Real Function rmass(dmax, iseed)
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  ! THE MINIMUM MASS FOR DELTA
  dmin = 1.078
  ! Delta(1232) production
  If (dmax<1.232) Then
     fm = fdelta(dmax)
  Else
     fm = 1.
  End If
  If (fm==0.) fm = 1.E-06
  ntry1 = 0
10 dm = ranart(nseed)*(dmax-dmin) + dmin
  ntry1 = ntry1 + 1
  If ((ranart(nseed)>fdelta(dm)/fm) .And. (ntry1<=10)) Goto 10
  !lin-2/26/03 sometimes Delta mass can reach very high values (e.g. 15.GeV),
  !     thus violating the thresh of the collision which produces it
  !     and leads to large violation of energy conservation.
  !     To limit the above, limit the Delta mass below a certain value
  !     (here taken as its central value + 2* B-W fullwidth):
  If (dm>1.47) Goto 10

  rmass = dm
  Return
End Function rmass

!------------------------------------------------------------------
! THE Breit Wigner FORMULA
Real Function frho(dmass)
  Save
  am0 = 0.77
  wid = 0.153
  fd = 0.25*wid**2/((dmass-am0)**2+0.25*wid**2)
  frho = fd
  Return
End Function frho
!**************************8
!   RHO MASS GENERATOR
Real Function rhomas(dmax, iseed)
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  ! THE MINIMUM MASS FOR DELTA
  dmin = 0.28
  ! RHO(770) production
  If (dmax<0.77) Then
     fm = frho(dmax)
  Else
     fm = 1.
  End If
  If (fm==0.) fm = 1.E-06
  ntry1 = 0
10 dm = ranart(nseed)*(dmax-dmin) + dmin
  ntry1 = ntry1 + 1
  If ((ranart(nseed)>frho(dm)/fm) .And. (ntry1<=10)) Goto 10
  !lin-2/26/03 limit the rho mass below a certain value
  !     (here taken as its central value + 2* B-W fullwidth):
  If (dm>1.07) Goto 10

  rhomas = dm
  Return
End Function rhomas
!*****************************************
! for pp-->pp+2pi
!      real*4 function X2pi(srt)
Real Function x2pi(srt)
  !  This function contains the experimental
  !     total pp-pp+pi(+)pi(-) Xsections    *
  !  srt    = DSQRT(s) in GeV                                                  *
  !  xsec   = production cross section in mb                                   *
  !  earray = EXPerimental table with proton momentum in GeV/c                 *
  !  xarray = EXPerimental table with cross sections in mb (curve to guide eye)*
  !                                                                            *
  !*****************************************
  !      real*4   xarray(15), earray(15)
  Real xarray(15), earray(15)
  Save
  Data earray/2.23, 2.81, 3.67, 4.0, 4.95, 5.52, 5.97, 6.04, 6.6, 6.9, 7.87, 8.11, 10.01, 16.0, 19./
  Data xarray/1.22, 2.51, 2.67, 2.95, 2.96, 2.84, 2.8, 3.2, 2.7, 3.0, 2.54, 2.46, 2.4, 1.66, 1.5/

  pmass = 0.9383
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  x2pi = 0.000001
  If (srt<=2.2) Return
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     x2pi = xarray(1)
     Return
  End If
  !
  ! 2.Interpolate double logarithmically to find sigma(srt)
  !
  Do ie = 1, 15
     If (earray(ie)==plab) Then
        x2pi = xarray(ie)
        Return
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        x2pi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Return
     End If
  End Do
  Return
End Function x2pi
!*****************************************
! for pp-->pn+pi(+)pi(+)pi(-)
!      real*4 function X3pi(srt)
Real Function x3pi(srt)
  !  This function contains the experimental pp->pp+3pi cross sections          *
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !  earray = EXPerimental table with proton energies in MeV                    *
  !  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
  !                                                                             *
  !*****************************************
  !      real*4   xarray(12), earray(12)
  Real xarray(12), earray(12)
  Save
  Data xarray/0.02, 0.4, 1.15, 1.60, 2.19, 2.85, 2.30, 3.10, 2.47, 2.60, 2.40, 1.70/
  Data earray/2.23, 2.81, 3.67, 4.00, 4.95, 5.52, 5.97, 6.04, 6.60, 6.90, 10.01, 19./

  pmass = 0.9383
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  x3pi = 1.E-06
  If (srt<=2.3) Return
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     x3pi = xarray(1)
     Return
  End If
  !
  ! 2.Interpolate double logarithmically to find sigma(srt)
  !
  Do ie = 1, 12
     If (earray(ie)==plab) Then
        x3pi = xarray(ie)
        Return
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        x3pi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Return
     End If
  End Do
  Return
End Function x3pi
!*****************************************
!*****************************************
! for pp-->pp+pi(+)pi(-)pi(0)
!      real*4 function X33pi(srt)
Real Function x33pi(srt)
  !  This function contains the experimental pp->pp+3pi cross sections          *
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !  earray = EXPerimental table with proton energies in MeV                    *
  !  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
  !                                                                             *
  !*****************************************
  !      real*4   xarray(12), earray(12)
  Real xarray(12), earray(12)
  Save
  Data xarray/0.02, 0.22, 0.74, 1.10, 1.76, 1.84, 2.20, 2.40, 2.15, 2.60, 2.30, 1.70/
  Data earray/2.23, 2.81, 3.67, 4.00, 4.95, 5.52, 5.97, 6.04, 6.60, 6.90, 10.01, 19./

  pmass = 0.9383
  x33pi = 1.E-06
  If (srt<=2.3) Return
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     x33pi = xarray(1)
     Return
  End If
  !
  ! 2.Interpolate double logarithmically to find sigma(srt)
  !
  Do ie = 1, 12
     If (earray(ie)==plab) Then
        x33pi = xarray(ie)
        Return
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        x33pi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Return
     End If
  End Do
  Return
End Function x33pi
!*****************************************
!       REAL*4 FUNCTION X4pi(SRT)
Real Function x4pi(srt)
  Save
  !       CROSS SECTION FOR NN-->DD+rho PROCESS
  ! *****************************
  akp = 0.498
  ak0 = 0.498
  ana = 0.94
  ada = 1.232
  al = 1.1157
  as = 1.1197
  pmass = 0.9383
  es = srt
  If (es<=4) Then
     x4pi = 0.
  Else
     ! cross section for two resonance pp-->DD+DN*+N*N*
     xpp2pi = 4.*x2pi(es)
     ! cross section for pp-->pp+spi
     xpp3pi = 3.*(x3pi(es)+x33pi(es))
     ! cross section for pp-->pD+ and nD++
     pps1 = sigma(es, 1, 1, 0) + 0.5*sigma(es, 1, 1, 1)
     pps2 = 1.5*sigma(es, 1, 1, 1)
     ppsngl = pps1 + pps2 + s1535(es)
     ! CROSS SECTION FOR KAON PRODUCTION from the four channels
     ! for NLK channel
     xk1 = 0
     xk2 = 0
     xk3 = 0
     xk4 = 0
     t1nlk = ana + al + akp
     t2nlk = ana + al - akp
     If (es<=t1nlk) Goto 333
     pmnlk2 = (es**2-t1nlk**2)*(es**2-t2nlk**2)/(4.*es**2)
     pmnlk = sqrt(pmnlk2)
     xk1 = pplpk(es)
     ! for DLK channel
     t1dlk = ada + al + akp
     t2dlk = ada + al - akp
     If (es<=t1dlk) Goto 333
     pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
     pmdlk = sqrt(pmdlk2)
     xk3 = pplpk(es)
     ! for NSK channel
     t1nsk = ana + as + akp
     t2nsk = ana + as - akp
     If (es<=t1nsk) Goto 333
     pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
     pmnsk = sqrt(pmnsk2)
     xk2 = ppk1(es) + ppk0(es)
     ! for DSK channel
     t1dsk = ada + as + akp
     t2dsk = ada + as - akp
     If (es<=t1dsk) Goto 333
     pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
     pmdsk = sqrt(pmdsk2)
     xk4 = ppk1(es) + ppk0(es)
     ! THE TOTAL KAON+ AND KAON0 PRODUCTION CROSS SECTION IS THEN
333  xkaon = 3.*(xk1+xk2+xk3+xk4)
     ! cross section for pp-->DD+rho
     x4pi = pp1(es) - ppsngl - xpp2pi - xpp3pi - xkaon
     If (x4pi<=0) x4pi = 1.E-06
  End If
  Return
End Function x4pi
!*****************************************
! for pp-->inelastic
!      real*4 function pp1(srt)
Real Function pp1(srt)
  Save
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !  earray = EXPerimental table with proton energies in MeV                    *
  !  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
  !                                                                             *
  !*****************************************
  pmass = 0.9383
  pp1 = 0.
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  plab2 = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (plab2<=0) Return
  plab = sqrt(plab2)
  pmin = 0.968
  pmax = 2080
  If ((plab<pmin) .Or. (plab>pmax)) Then
     pp1 = 0.
     Return
  End If
  !* fit parameters
  a = 30.9
  b = -28.9
  c = 0.192
  d = -0.835
  an = -2.46
  pp1 = a + b*(plab**an) + c*(alog(plab))**2
  If (pp1<=0) pp1 = 0.0
  Return
End Function pp1
!*****************************************
! for pp-->elastic
!      real*4 function pp2(srt)
Real Function pp2(srt)
  Save
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !  earray = EXPerimental table with proton energies in MeV                    *
  !  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
  !                                                                             *
  !*****************************************
  pmass = 0.9383
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)

  !lin-9/2012: check argument in sqrt():
  scheck = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (scheck<0) Then
     Write (99, *) 'scheck33: ', scheck
     scheck = 0.
  End If
  plab = sqrt(scheck)
  !      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)

  pmin = 2.
  pmax = 2050
  If (plab>pmax) Then
     pp2 = 8.
     Return
  End If
  If (plab<pmin) Then
     pp2 = 25.
     Return
  End If
  !* fit parameters
  a = 11.2
  b = 25.5
  c = 0.151
  d = -1.62
  an = -1.12
  pp2 = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (pp2<=0) pp2 = 0
  Return
End Function pp2

!*****************************************
! for pp-->total
!      real*4 function ppt(srt)
Real Function ppt(srt)
  Save
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !  earray = EXPerimental table with proton energies in MeV                    *
  !  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
  !                                                                             *
  !*****************************************
  pmass = 0.9383
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)

  !lin-9/2012: check argument in sqrt():
  scheck = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (scheck<0) Then
     Write (99, *) 'scheck34: ', scheck
     scheck = 0.
  End If
  plab = sqrt(scheck)
  !      plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)

  pmin = 3.
  pmax = 2100
  If ((plab<pmin) .Or. (plab>pmax)) Then
     ppt = 55.
     Return
  End If
  !* fit parameters
  a = 45.6
  b = 219.0
  c = 0.410
  d = -3.41
  an = -4.23
  ppt = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (ppt<=0) ppt = 0.0
  Return
End Function ppt

!************************
! cross section for N*(1535) production in PP collisions
! VARIABLES:
! LB1,LB2 ARE THE LABLES OF THE TWO COLLIDING PARTICLES
! SRT IS THE CMS ENERGY
! X1535 IS THE N*(1535) PRODUCTION CROSS SECTION
! NOTE THAT THE N*(1535) PRODUCTION CROSS SECTION IS 2 TIMES THE ETA
! PRODUCTION CROSS SECTION
! DATE: Aug. 1 , 1994
! ********************************
Real Function s1535(srt)
  Save
  s0 = 2.424
  s1535 = 0.
  If (srt<=s0) Return
  s1535 = 2.*0.102*(srt-s0)/(0.058+(srt-s0)**2)
  Return
End Function s1535
!***************************************
! generate a table for pt distribution for
Subroutine tablem
  ! THE PROCESS N+N--->N+N+PION
  !       DATE : July 11, 1994
  !****************************************
  Common /table/xarray(0:1000), earray(0:1000)
  !c      SAVE /TABLE/
  Save
  ptmax = 2.01
  anorm = ptdis(ptmax)
  Do l = 0, 200
     x = 0.01*float(l+1)
     rr = ptdis(x)/anorm
     earray(l) = rr
     xarray(l) = x
  End Do
  Return
End Subroutine tablem
!********************************
Real Function ptdis(x)
  Save
  ! NUCLEON TRANSVERSE MOMENTUM DISTRIBUTION AT HIGH ENERGIES
  ! DATE: Aug. 11, 1994
  !********************************
  b = 3.78
  c = 0.47
  d = 3.60
  !       b=b*3
  !       d=d*3
  ptdis = 1./(2.*b)*(1.-exp(-b*x**2)) - c/d*x*exp(-d*x) - c/d**2*(exp(-d*x)-1.)
  Return
End Function ptdis
!****************************
Subroutine ppxs(lb1, lb2, srt, ppsig, spprho, ipp)
  ! purpose: this subroutine gives the cross section for pion+pion
  !          elastic collision
  ! variables:
  !       input: lb1,lb2 and srt are the labels and srt for I1 and I2
  !       output: ppsig: pp xsection
  !               ipp: label for the pion+pion channel
  !               Ipp=0 NOTHING HAPPEND
  !                  1 for Pi(+)+PI(+) DIRECT
  !                   2     PI(+)+PI(0) FORMING RHO(+)
  !                  3     PI(+)+PI(-) FORMING RHO(0)
  !                   4     PI(0)+PI(O) DIRECT
  !                  5     PI(0)+PI(-) FORMING RHO(-)
  !                  6     PI(-)+PI(-) DIRECT
  ! reference: G.F. Bertsch, Phys. Rev. D37 (1988) 1202.
  ! date : Aug 29, 1994
  !****************************
  Parameter (amp=0.14, pi=3.1415926)
  Save
  ppsig = 0.0

  !bzdbg10/15/99
  spprho = 0.0
  !bzdbg10/15/99 end

  ipp = 0
  If (srt<=0.3) Return
  q = sqrt((srt/2)**2-amp**2)
  esigma = 5.8*amp
  tsigma = 2.06*q
  erho = 0.77
  trho = 0.095*q*(q/amp/(1.+(q/erho)**2))**2
  esi = esigma - srt
  If (esi==0) Then
     d00 = pi/2.
     Goto 10
  End If
  d00 = atan(tsigma/2./esi)
10 erh = erho - srt
  If (erh==0.) Then
     d11 = pi/2.
     Goto 20
  End If
  d11 = atan(trho/2./erh)
20 d20 = -0.12*q/amp
  s0 = 8.*pi*sin(d00)**2/q**2
  s1 = 8*pi*3*sin(d11)**2/q**2
  s2 = 8*pi*5*sin(d20)**2/q**2
  !    !! GeV^-2 to mb
  s0 = s0*0.197**2*10.
  s1 = s1*0.197**2*10.
  s2 = s2*0.197**2*10.
  !       ppXS=s0/9.+s1/3.+s2*0.56
  !       if(ppxs.le.0)ppxs=0.00001
  spprho = s1/2.
  ! (1) PI(+)+PI(+)
  If (lb1==5 .And. lb2==5) Then
     ipp = 1
     ppsig = s2
     Return
  End If
  ! (2) PI(+)+PI(0)
  If ((lb1==5 .And. lb2==4) .Or. (lb1==4 .And. lb2==5)) Then
     ipp = 2
     ppsig = s2/2. + s1/2.
     Return
  End If
  ! (3) PI(+)+PI(-)
  If ((lb1==5 .And. lb2==3) .Or. (lb1==3 .And. lb2==5)) Then
     ipp = 3
     ppsig = s2/6. + s1/2. + s0/3.
     Return
  End If
  ! (4) PI(0)+PI(0)
  If (lb1==4 .And. lb2==4) Then
     ipp = 4
     ppsig = 2*s2/3. + s0/3.
     Return
  End If
  ! (5) PI(0)+PI(-)
  If ((lb1==4 .And. lb2==3) .Or. (lb1==3 .And. lb2==4)) Then
     ipp = 5
     ppsig = s2/2. + s1/2.
     Return
  End If
  ! (6) PI(-)+PI(-)
  If (lb1==3 .And. lb2==3) Then
     ipp = 6
     ppsig = s2
  End If
  Return
End Subroutine ppxs
!*********************************
! elementary kaon production cross sections
!  from the CERN data book
!  date: Sept.2, 1994
!  for pp-->pLK+
!      real*4 function pplpk(srt)
Real Function pplpk(srt)
  Save
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !  earray = EXPerimental table with proton energies in MeV                    *
  !  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
  !                                                                             *
  !*****************************************
  pmass = 0.9383
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !   find the center of mass energy corresponding to the given pm as
  !   if Lambda+N+K are produced
  pplpk = 0.

  !lin-9/2012: check argument in sqrt():
  scheck = ((srt**2-2.*pmass**2)/(2.*pmass))**2 - pmass**2
  If (scheck<0) Then
     Write (99, *) 'scheck35: ', scheck
     scheck = 0.
  End If
  plab = sqrt(scheck)
  !        plab=sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)

  pmin = 2.82
  pmax = 25.0
  If (plab>pmax) Then
     pplpk = 0.036
     Return
  End If
  If (plab<pmin) Then
     pplpk = 0.
     Return
  End If
  !* fit parameters
  a = 0.0654
  b = -3.16
  c = -0.0029
  an = -4.14
  pplpk = a + b*(plab**an) + c*(alog(plab))**2
  If (pplpk<=0) pplpk = 0
  Return
End Function pplpk

!*****************************************
! for pp-->pSigma+K0
!      real*4 function ppk0(srt)
Real Function ppk0(srt)
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !                                                                             *
  !*****************************************
  !      real*4   xarray(7), earray(7)
  Real xarray(7), earray(7)
  Save
  Data xarray/0.030, 0.025, 0.025, 0.026, 0.02, 0.014, 0.06/
  Data earray/3.67, 4.95, 5.52, 6.05, 6.92, 7.87, 10./

  pmass = 0.9383
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  ppk0 = 0
  If (srt<=2.63) Return
  If (srt>4.54) Then
     ppk0 = 0.037
     Return
  End If
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     ppk0 = xarray(1)
     Return
  End If
  !
  ! 2.Interpolate double logarithmically to find sigma(srt)
  !
  Do ie = 1, 7
     If (earray(ie)==plab) Then
        ppk0 = xarray(ie)
        Goto 10
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        ppk0 = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Goto 10
     End If
  End Do
10 Continue
  Return
End Function ppk0
!*****************************************
! for pp-->pSigma0K+
!      real*4 function ppk1(srt)
Real Function ppk1(srt)
  !  srt    = DSQRT(s) in GeV                                                   *
  !  xsec   = production cross section in mb                                    *
  !                                                                             *
  !*****************************************
  !      real*4   xarray(7), earray(7)
  Real xarray(7), earray(7)
  Save
  Data xarray/0.013, 0.025, 0.016, 0.012, 0.017, 0.029, 0.025/
  Data earray/3.67, 4.95, 5.52, 5.97, 6.05, 6.92, 7.87/

  pmass = 0.9383
  ! 1.Calculate p(lab)  from srt [GeV]
  !   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  !      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  ppk1 = 0.
  If (srt<=2.63) Return
  If (srt>4.08) Then
     ppk1 = 0.025
     Return
  End If
  plab = sqrt(((srt**2-2.*pmass**2)/(2.*pmass))**2-pmass**2)
  If (plab<earray(1)) Then
     ppk1 = xarray(1)
     Return
  End If
  !
  ! 2.Interpolate double logarithmically to find sigma(srt)
  !
  Do ie = 1, 7
     If (earray(ie)==plab) Then
        ppk1 = xarray(ie)
        Goto 10
     Else If (earray(ie)>plab) Then
        ymin = alog(xarray(ie-1))
        ymax = alog(xarray(ie))
        xmin = alog(earray(ie-1))
        xmax = alog(earray(ie))
        ppk1 = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
        Goto 10
     End If
  End Do
10 Continue
  Return
End Function ppk1
!*********************************
!                                                                      *
!                                                                      *
Subroutine crpn(px, py, pz, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)
  !     PURPOSE:                                                         *
  !           DEALING WITH PION+N-->L/S+KAON PROCESS AND PION PRODUCTION *
  !     NOTE   :                                                         *
  !
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                     7  PION+N-->L/S+KAON
  !           iblock   - 77 pion+N-->Delta+pion
  !           iblock   - 78 pion+N-->Delta+RHO
  !           iblock   - 79 pion+N-->Delta+OMEGA
  !           iblock   - 222 pion+N-->Phi
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  iblock = 1
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
  If (xkaon0/(xkaon+xphi)>=x1) Then
     ! kaon production
     !-----------------------------------------------------------------------
     iblock = 7
     If (ianti==1) iblock = -7
     ntag = 0
     ! RELABLE PARTICLES FOR THE PROCESS PION+n-->LAMBDA K OR SIGMA k
     ! DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW
     ! MOMENTA FOR PARTICLES IN THE FINAL STATE.
     kaonc = 0
     If (pnlka(srt)/(pnlka(srt)+pnska(srt))>ranart(nseed)) kaonc = 1
     If (e(i1)<=0.2) Then
        lb(i1) = 23
        e(i1) = aka
        If (kaonc==1) Then
           lb(i2) = 14
           e(i2) = ala
        Else
           lb(i2) = 15 + int(3*ranart(nseed))
           e(i2) = asa
        End If
        If (ianti==1) Then
           lb(i1) = 21
           lb(i2) = -lb(i2)
        End If
     Else
        lb(i2) = 23
        e(i2) = aka
        If (kaonc==1) Then
           lb(i1) = 14
           e(i1) = ala
        Else
           lb(i1) = 15 + int(3*ranart(nseed))
           e(i1) = asa
        End If
        If (ianti==1) Then
           lb(i2) = 21
           lb(i1) = -lb(i1)
        End If
     End If
     em1 = e(i1)
     em2 = e(i2)
     Goto 50
     ! to gererate the momentum for the kaon and L/S
  Else If (xphi/(xkaon+xphi)>=x1) Then
     iblock = 222
     If (xphin/xphi>=ranart(nseed)) Then
        lb(i1) = 1 + int(2*ranart(nseed))
        e(i1) = amn
     Else
        lb(i1) = 6 + int(4*ranart(nseed))
        e(i1) = am0
     End If
     !  !! at present only baryon
     If (ianti==1) lb(i1) = -lb(i1)
     lb(i2) = 29
     e(i2) = aphi
     em1 = e(i1)
     em2 = e(i2)
     Goto 50
  Else
     ! CHECK WHAT KIND OF PION PRODUCTION PROCESS HAS HAPPENED
     If (ranart(nseed)<=twopi(srt)/(twopi(srt)+threpi(srt)+fourpi(srt))) Then
        iblock = 77
     Else
        If (threpi(srt)/(threpi(srt)+fourpi(srt))>ranart(nseed)) Then
           iblock = 78
        Else
           iblock = 79
        End If
     End If
     ntag = 0
     ! pion production (Delta+pion/rho/omega in the final state)
     ! generate the mass of the delta resonance
     x2 = ranart(nseed)
     ! relable the particles
     If (iblock==77) Then
        ! GENERATE THE DELTA MASS
        dmax = srt - ap1 - 0.02
        dm = rmass(dmax, iseed)
        ! pion+baryon-->pion+delta
        ! Relable particles, I1 is assigned to the Delta and I2 is assigned to the
        ! meson
        !(1) for pi(+)+p-->D(+)+P(+) OR D(++)+p(0)
        If (((lb(i1)==1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              Else
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 4
                 ipi = 4
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              Else
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        !(2) for pi(-)+p-->D(0)+P(0) OR D(+)+p(-),or D(-)+p(+)
        If (((lb(i1)==1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        !(3) for pi(+)+n-->D(+)+Pi(0) OR D(++)+p(-) or D(0)+pi(+)
        If (((lb(i1)==2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        !(4) for pi(0)+p-->D(+)+Pi(0) OR D(++)+p(-) or D(0)+pi(+)
        If ((iabs(lb(i1))==1 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==1)) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        !(5) for pi(-)+n-->D(-)+P(0) OR D(0)+p(-)
        If (((lb(i1)==2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              Else
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              Else
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
        !(6) for pi(0)+n-->D(0)+P(0), D(-)+p(+) or D(+)+p(-)
        If ((iabs(lb(i1))==2 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==2)) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 4
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2<=0.67 .And. x2>0.33) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 5
                 e(i2) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 3
                 e(i2) = ap1
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 4
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2<=0.67 .And. x2>0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 5
                 e(i1) = ap1
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 3
                 e(i1) = ap1
                 Goto 40
              End If
           End If
        End If
     End If
     If (iblock==78) Then
        Call rmasdd(srt, 1.232, 0.77, 1.08, 0.28, iseed, 4, dm, ameson)
        arho = ameson
        ! pion+baryon-->Rho+delta
        !(1) for pi(+)+p-->D(+)+rho(+) OR D(++)+rho(0)
        If (((lb(i1)==1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              Else
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              Else
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        !(2) for pi(-)+p-->D(+)+rho(-) OR D(0)+rho(0) or D(-)+rho(+)
        If (((lb(i1)==1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        !(3) for pi(+)+n-->D(+)+rho(0) OR D(++)+rho(-) or D(0)+rho(+)
        If (((lb(i1)==2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        !(4) for pi(0)+p-->D(+)+rho(0) OR D(++)+rho(-) or D(0)+rho(+)
        If ((iabs(lb(i1))==1 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==1)) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 9
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 9
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        !(5) for pi(-)+n-->D(-)+rho(0) OR D(0)+rho(-)
        If (((lb(i1)==2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.5) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              Else
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
                 Goto 40
              End If
           Else
              ii = i2
              If (x2<=0.5) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              Else
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
                 Goto 40
              End If
           End If
        End If
        !(6) for pi(0)+n-->D(0)+rho(0), D(-)+rho(+) and D(+)+rho(-)
        If ((iabs(lb(i1))==2 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==2)) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              If (x2<=0.33) Then
                 lb(i1) = 7
                 e(i1) = dm
                 lb(i2) = 26
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.33 .And. x2<=0.67) Then
                 lb(i1) = 6
                 e(i1) = dm
                 lb(i2) = 27
                 e(i2) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i1) = 8
                 e(i1) = dm
                 lb(i2) = 25
                 e(i2) = arho
              End If
           Else
              ii = i2
              If (x2<=0.33) Then
                 lb(i2) = 7
                 e(i2) = dm
                 lb(i1) = 26
                 e(i1) = arho
                 Goto 40
              End If
              If (x2<=0.67 .And. x2>0.33) Then
                 lb(i2) = 6
                 e(i2) = dm
                 lb(i1) = 27
                 e(i1) = arho
                 Goto 40
              End If
              If (x2>0.67) Then
                 lb(i2) = 8
                 e(i2) = dm
                 lb(i1) = 25
                 e(i1) = arho
              End If
           End If
        End If
     End If
     If (iblock==79) Then
        aomega = 0.782
        ! GENERATE THE DELTA MASS
        dmax = srt - 0.782 - 0.02
        dm = rmass(dmax, iseed)
        ! pion+baryon-->omega+delta
        !(1) for pi(+)+p-->D(++)+omega(0)
        If (((lb(i1)==1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              lb(i1) = 9
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 9
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        !(2) for pi(-)+p-->D(0)+omega(0)
        If (((lb(i1)==1 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==1)) .Or. ((lb(i1)==-1 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-1))) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              lb(i1) = 7
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 7
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        !(3) for pi(+)+n-->D(+)+omega(0)
        If (((lb(i1)==2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              lb(i1) = 8
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 8
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        !(4) for pi(0)+p-->D(+)+omega(0)
        If ((iabs(lb(i1))==1 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==1)) Then
           If (iabs(lb(i1))==1) Then
              ii = i1
              lb(i1) = 8
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 8
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
              Goto 40
           End If
        End If
        !(5) for pi(-)+n-->D(-)+omega(0)
        If (((lb(i1)==2 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==2)) .Or. ((lb(i1)==-2 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-2))) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              lb(i1) = 6
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 6
              e(i2) = dm
              lb(i1) = 28
              e(i1) = aomega
           End If
        End If
        !(6) for pi(0)+n-->D(0)+omega(0)
        If ((iabs(lb(i1))==2 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==2)) Then
           If (iabs(lb(i1))==2) Then
              ii = i1
              lb(i1) = 7
              e(i1) = dm
              lb(i2) = 28
              e(i2) = aomega
              Goto 40
           Else
              ii = i2
              lb(i2) = 7
              e(i2) = dm
              lb(i1) = 26
              e(i1) = arho
              Goto 40
           End If
        End If
     End If
40   em1 = e(i1)
     em2 = e(i2)
     If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
        lb(ii) = -lb(ii)
        jj = i2
        If (ii==i2) jj = i1
        If (iblock==77) Then
           If (lb(jj)==3) Then
              lb(jj) = 5
           Else If (lb(jj)==5) Then
              lb(jj) = 3
           End If
        Else If (iblock==78) Then
           If (lb(jj)==25) Then
              lb(jj) = 27
           Else If (lb(jj)==27) Then
              lb(jj) = 25
           End If
        End If
     End If
  End If
  !-----------------------------------------------------------------------
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
50 pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 0.00000001
  pr = sqrt(pr2)/(2.*srt)
  ! here we use the same transverse momentum distribution as for
  ! pp collisions, it might be necessary to use a different distribution

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pr
  !         cc1=ptr(0.33*pr,iseed)
  cc1 = ptr(xptr, iseed)
  !lin-10/25/02-end

  !lin-9/2012: check argument in sqrt():
  scheck = pr**2 - cc1**2
  If (scheck<0) Then
     Write (99, *) 'scheck36: ', scheck
     scheck = 0.
  End If
  c1 = sqrt(scheck)/pr
  !         c1=sqrt(pr**2-cc1**2)/pr

  !          C1   = 1.0 - 2.0 * RANART(NSEED)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  ! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  ! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crpn
!*********************************
!                                                                      *
!                                                                      *
Subroutine cren(px, py, pz, srt, i1, i2, iblock)
  !     PURPOSE:                                                         *
  !             DEALING WITH ETA+N-->L/S+KAON PROCESS                   *
  !     NOTE   :                                                         *
  !
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                     7  ETA+N-->L/S+KAON
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  ntag = 0
  iblock = 7
  ianti = 0
  If (lb(i1)<0 .Or. lb(i2)<0) Then
     ianti = 1
     iblock = -7
  End If
  ! RELABLE PARTICLES FOR THE PROCESS eta+n-->LAMBDA K OR SIGMA k
  ! DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW
  ! MOMENTA FOR PARTICLES IN THE FINAL STATE.
  kaonc = 0
  If (pnlka(srt)/(pnlka(srt)+pnska(srt))>ranart(nseed)) kaonc = 1
  If (e(i1)<=0.6) Then
     lb(i1) = 23
     e(i1) = aka
     If (kaonc==1) Then
        lb(i2) = 14
        e(i2) = ala
     Else
        lb(i2) = 15 + int(3*ranart(nseed))
        e(i2) = asa
     End If
     If (ianti==1) Then
        lb(i1) = 21
        lb(i2) = -lb(i2)
     End If
  Else
     lb(i2) = 23
     e(i2) = aka
     If (kaonc==1) Then
        lb(i1) = 14
        e(i1) = ala
     Else
        lb(i1) = 15 + int(3*ranart(nseed))
        e(i1) = asa
     End If
     If (ianti==1) Then
        lb(i2) = 21
        lb(i1) = -lb(i1)
     End If
  End If
  em1 = e(i1)
  em2 = e(i2)
  !-----------------------------------------------------------------------
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  ! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  ! FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE
  Return
End Subroutine cren
!*********************************
!                                                                      *
!                                                                      *
!      SUBROUTINE Crdir(PX,PY,PZ,SRT,I1,I2)
Subroutine crdir(px, py, pz, srt, i1, i2, iblock)
  !     PURPOSE:                                                         *
  !             DEALING WITH pion+N-->pion+N PROCESS                   *
  !     NOTE   :                                                         *
  !
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  iblock = 999
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  !-----------------------------------------------------------------------
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pr
  !         cc1=ptr(0.33*pr,iseed)
  cc1 = ptr(xptr, iseed)
  !lin-10/25/02-end

  !lin-9/2012: check argument in sqrt():
  scheck = pr**2 - cc1**2
  If (scheck<0) Then
     Write (99, *) 'scheck37: ', scheck
     scheck = 0.
  End If
  c1 = sqrt(scheck)/pr
  !         c1=sqrt(pr**2-cc1**2)/pr

  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  ! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  ! ROTATE the momentum
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crdir
!*********************************


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine crpd(px, py, pz, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)
  !     PURPOSE:                                                         *
  !     DEALING WITH PION+D(N*)-->PION +N OR
  !                                             L/S+KAON PROCESS         *
  !     NOTE   :                                                         *
  !
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                     7  PION+D(N*)-->L/S+KAON
  !           iblock   - 80 pion+D(N*)-->pion+N
  !           iblock   - 81 RHO+D(N*)-->PION+N
  !           iblock   - 82 OMEGA+D(N*)-->PION+N
  !                     222  PION+D --> PHI
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  iblock = 1
  x1 = ranart(nseed)
  ianti = 0
  If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
  If (xkaon0/(xkaon+xphi)>=x1) Then
     ! kaon production
     !-----------------------------------------------------------------------
     iblock = 7
     If (ianti==1) iblock = -7
     ntag = 0
     ! RELABLE PARTICLES FOR THE PROCESS PION+n-->LAMBDA K OR SIGMA k
     ! DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW
     ! MOMENTA FOR PARTICLES IN THE FINAL STATE.
     kaonc = 0
     If (pnlka(srt)/(pnlka(srt)+pnska(srt))>ranart(nseed)) kaonc = 1
     !lin-8/17/00     & +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1
     If (e(i1)<=0.2) Then
        lb(i1) = 23
        e(i1) = aka
        If (kaonc==1) Then
           lb(i2) = 14
           e(i2) = ala
        Else
           lb(i2) = 15 + int(3*ranart(nseed))
           e(i2) = asa
        End If
        If (ianti==1) Then
           lb(i1) = 21
           lb(i2) = -lb(i2)
        End If
     Else
        lb(i2) = 23
        e(i2) = aka
        If (kaonc==1) Then
           lb(i1) = 14
           e(i1) = ala
        Else
           lb(i1) = 15 + int(3*ranart(nseed))
           e(i1) = asa
        End If
        If (ianti==1) Then
           lb(i2) = 21
           lb(i1) = -lb(i1)
        End If
     End If
     em1 = e(i1)
     em2 = e(i2)
     Goto 50
     ! to gererate the momentum for the kaon and L/S
     !
     !* Phi production
  Else If (xphi/(xkaon+xphi)>=x1) Then
     iblock = 222
     If (xphin/xphi>=ranart(nseed)) Then
        lb(i1) = 1 + int(2*ranart(nseed))
        e(i1) = amn
     Else
        lb(i1) = 6 + int(4*ranart(nseed))
        e(i1) = am0
     End If
     !   !! at present only baryon
     If (ianti==1) lb(i1) = -lb(i1)
     lb(i2) = 29
     e(i2) = aphi
     em1 = e(i1)
     em2 = e(i2)
     Goto 50
  Else
     ! PION REABSORPTION HAS HAPPENED
     x2 = ranart(nseed)
     iblock = 80
     ntag = 0
     ! Relable particles, I1 is assigned to the nucleon
     ! and I2 is assigned to the pion
     ! for the reverse of the following process
     !(1) for D(+)+P(+)-->p+pion(+)
     If (((lb(i1)==8 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==8)) .Or. ((lb(i1)==-8 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-8))) Then
        If (iabs(lb(i1))==8) Then
           ii = i1
           lb(i1) = 1
           e(i1) = amn
           lb(i2) = 5
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 1
           e(i2) = amn
           lb(i1) = 5
           e(i1) = ap1
           Goto 40
        End If
     End If
     !
     !(2) for D(0)+P(0)-->n+pi(0) or p+pi(-)
     If ((iabs(lb(i1))==7 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==7)) Then
        If (iabs(lb(i1))==7) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(3) for D(+)+Pi(0)-->pi(+)+n or pi(0)+p
     If ((iabs(lb(i1))==8 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==8)) Then
        If (iabs(lb(i1))==8) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(4) for D(-)+Pi(0)-->n+pi(-)
     If ((iabs(lb(i1))==6 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==6)) Then
        If (iabs(lb(i1))==6) Then
           ii = i1
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 3
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 3
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(5) for D(+)+Pi(-)-->pi(0)+n or pi(-)+p
     If (((lb(i1)==8 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==8)) .Or. ((lb(i1)==-8 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-8))) Then
        If (iabs(lb(i1))==8) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(6) D(0)+P(+)-->n+pi(+) or p+pi(0)
     If (((lb(i1)==7 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==7)) .Or. ((lb(i1)==-7 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-7))) Then
        If (iabs(lb(i1))==7) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(7) for D(0)+Pi(-)-->n+pi(-)
     If (((lb(i1)==7 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==7)) .Or. ((lb(i1)==-7 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-7))) Then
        If (iabs(lb(i1))==7) Then
           ii = i1
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 3
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 3
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(8) D(-)+P(+)-->n+pi(0) or p+pi(-)
     If (((lb(i1)==6 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==6)) .Or. ((lb(i1)==-6 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-6))) Then
        If (iabs(lb(i1))==6) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !
     !(9) D(++)+P(-)-->n+pi(+) or p+pi(0)
     If (((lb(i1)==9 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==9)) .Or. ((lb(i1)==-9 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-9))) Then
        If (iabs(lb(i1))==9) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(10) for D(++)+Pi(0)-->p+pi(+)
     If ((iabs(lb(i1))==9 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==9)) Then
        If (iabs(lb(i1))==9) Then
           ii = i1
           lb(i1) = 1
           e(i1) = amn
           lb(i2) = 5
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 1
           e(i2) = amn
           lb(i1) = 5
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(11) for N*(1440)(+)or N*(1535)(+)+P(+)-->p+pion(+)
     If (((lb(i1)==11 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==11) .Or. (lb(i1)==13 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==13)) .Or. ((lb(i1)==-11 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-11) .Or. (lb(i1)==-13 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-13))) Then
        If (iabs(lb(i1))==11 .Or. iabs(lb(i1))==13) Then
           ii = i1
           lb(i1) = 1
           e(i1) = amn
           lb(i2) = 5
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 1
           e(i2) = amn
           lb(i1) = 5
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(12) for N*(1440) or N*(1535)(0)+P(0)-->n+pi(0) or p+pi(-)
     If ((iabs(lb(i1))==10 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==10) .Or. (lb(i1)==4 .And. iabs(lb(i2))==12) .Or. (lb(i2)==4 .And. iabs(lb(i1))==12)) Then
        If (iabs(lb(i1))==10 .Or. iabs(lb(i1))==12) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(13) for N*(1440) or N*(1535)(+)+Pi(0)-->pi(+)+n or pi(0)+p
     If ((iabs(lb(i1))==11 .And. lb(i2)==4) .Or. (lb(i1)==4 .And. iabs(lb(i2))==11) .Or. (lb(i1)==4 .And. iabs(lb(i2))==13) .Or. (lb(i2)==4 .And. iabs(lb(i1))==13)) Then
        If (iabs(lb(i1))==11 .Or. iabs(lb(i1))==13) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(14) for N*(1440) or N*(1535)(+)+Pi(-)-->pi(0)+n or pi(-)+p
     If (((lb(i1)==11 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==11) .Or. (lb(i1)==3 .And. lb(i2)==13) .Or. (lb(i2)==3 .And. lb(i1)==13)) .Or. ((lb(i1)==-11 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-11) .Or. (lb(i1)==5 .And. lb(i2)==-13) .Or. (lb(i2)==5 .And. lb(i1)==-13))) Then
        If (iabs(lb(i1))==11 .Or. iabs(lb(i1))==13) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(15) N*(1440) or N*(1535)(0)+P(+)-->n+pi(+) or p+pi(0)
     If (((lb(i1)==10 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==10) .Or. (lb(i1)==12 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==12)) .Or. ((lb(i1)==-10 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-10) .Or. (lb(i1)==-12 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==-12))) Then
        If (iabs(lb(i1))==10 .Or. iabs(lb(i1))==12) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(16) for N*(1440) or N*(1535) (0)+Pi(-)-->n+pi(-)
     If (((lb(i1)==10 .And. lb(i2)==3) .Or. (lb(i1)==3 .And. lb(i2)==10) .Or. (lb(i1)==3 .And. lb(i2)==12) .Or. (lb(i1)==12 .And. lb(i2)==3)) .Or. ((lb(i1)==-10 .And. lb(i2)==5) .Or. (lb(i1)==5 .And. lb(i2)==-10) .Or. (lb(i1)==5 .And. lb(i2)==-12) .Or. (lb(i1)==-12 .And. lb(i2)==5))) Then
        If (iabs(lb(i1))==10 .Or. iabs(lb(i1))==12) Then
           ii = i1
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 3
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 3
           e(i1) = ap1
           Goto 40
        End If
     End If
40   em1 = e(i1)
     em2 = e(i2)
     If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
        lb(ii) = -lb(ii)
        jj = i2
        If (ii==i2) jj = i1
        If (lb(jj)==3) Then
           lb(jj) = 5
        Else If (lb(jj)==5) Then
           lb(jj) = 3
        End If
     End If
  End If
  !-----------------------------------------------------------------------
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
50 pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pr
  !         cc1=ptr(0.33*pr,iseed)
  cc1 = ptr(xptr, iseed)
  !lin-10/25/02-end

  !lin-9/2012: check argument in sqrt():
  scheck = pr**2 - cc1**2
  If (scheck<0) Then
     Write (99, *) 'scheck38: ', scheck
     scheck = 0.
  End If
  c1 = sqrt(scheck)/pr
  !         c1=sqrt(pr**2-cc1**2)/pr

  !         C1   = 1.0 - 2.0 * RANART(NSEED)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  ! rotate the momentum
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crpd
!*********************************
!                                                                      *
!                                                                      *
Subroutine crrd(px, py, pz, srt, i1, i2, iblock, xkaon0, xkaon, xphi, xphin)
  !     PURPOSE:                                                         *
  !     DEALING WITH rho(omega)+N or D(N*)-->PION +N OR
  !                                             L/S+KAON PROCESS         *
  !     NOTE   :                                                         *
  !
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                     7  rho(omega)+N or D(N*)-->L/S+KAON
  !           iblock   - 80 pion+D(N*)-->pion+N
  !           iblock   - 81 RHO+D(N*)-->PION+N
  !           iblock   - 82 OMEGA+D(N*)-->PION+N
  !           iblock   - 222 pion+N-->Phi
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, aphi=1.02)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  iblock = 1
  ianti = 0
  If (lb(i1)<0 .Or. lb(i2)<0) ianti = 1
  x1 = ranart(nseed)
  If (xkaon0/(xkaon+xphi)>=x1) Then
     ! kaon production
     !-----------------------------------------------------------------------
     iblock = 7
     If (ianti==1) iblock = -7
     ntag = 0
     ! RELABLE PARTICLES FOR THE PROCESS PION+n-->LAMBDA K OR SIGMA k
     ! DECIDE LAMBDA OR SIGMA PRODUCTION, AND TO CALCULATE THE NEW
     ! MOMENTA FOR PARTICLES IN THE FINAL STATE.
     kaonc = 0
     If (pnlka(srt)/(pnlka(srt)+pnska(srt))>ranart(nseed)) kaonc = 1
     !lin-8/17/00     & +PNSKA(SRT)).GT.RANART(NSEED))KAONC=1
     If (e(i1)<=0.92) Then
        lb(i1) = 23
        e(i1) = aka
        If (kaonc==1) Then
           lb(i2) = 14
           e(i2) = ala
        Else
           lb(i2) = 15 + int(3*ranart(nseed))
           e(i2) = asa
        End If
        If (ianti==1) Then
           lb(i1) = 21
           lb(i2) = -lb(i2)
        End If
     Else
        lb(i2) = 23
        e(i2) = aka
        If (kaonc==1) Then
           lb(i1) = 14
           e(i1) = ala
        Else
           lb(i1) = 15 + int(3*ranart(nseed))
           e(i1) = asa
        End If
        If (ianti==1) Then
           lb(i2) = 21
           lb(i1) = -lb(i1)
        End If
     End If
     em1 = e(i1)
     em2 = e(i2)
     Goto 50
     ! to gererate the momentum for the kaon and L/S
     !
     !* Phi production
  Else If (xphi/(xkaon+xphi)>=x1) Then
     iblock = 222
     If (xphin/xphi>=ranart(nseed)) Then
        lb(i1) = 1 + int(2*ranart(nseed))
        e(i1) = amn
     Else
        lb(i1) = 6 + int(4*ranart(nseed))
        e(i1) = am0
     End If
     !   !! at present only baryon
     If (ianti==1) lb(i1) = -lb(i1)
     lb(i2) = 29
     e(i2) = aphi
     em1 = e(i1)
     em2 = e(i2)
     Goto 50
  Else
     ! rho(omega) REABSORPTION HAS HAPPENED
     x2 = ranart(nseed)
     iblock = 81
     ntag = 0
     If (lb(i1)==28 .Or. lb(i2)==28) Goto 60
     ! we treat Rho reabsorption in the following
     ! Relable particles, I1 is assigned to the Delta
     ! and I2 is assigned to the meson
     ! for the reverse of the following process
     !(1) for D(+)+rho(+)-->p+pion(+)
     If (((lb(i1)==8 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==8)) .Or. ((lb(i1)==-8 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==-8))) Then
        If (iabs(lb(i1))==8) Then
           ii = i1
           lb(i1) = 1
           e(i1) = amn
           lb(i2) = 5
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 1
           e(i2) = amn
           lb(i1) = 5
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(2) for D(0)+rho(0)-->n+pi(0) or p+pi(-)
     If ((iabs(lb(i1))==7 .And. lb(i2)==26) .Or. (lb(i1)==26 .And. iabs(lb(i2))==7)) Then
        If (iabs(lb(i1))==7) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(3) for D(+)+rho(0)-->pi(+)+n or pi(0)+p
     If ((iabs(lb(i1))==8 .And. lb(i2)==26) .Or. (lb(i1)==26 .And. iabs(lb(i2))==8)) Then
        If (iabs(lb(i1))==8) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(4) for D(-)+rho(0)-->n+pi(-)
     If ((iabs(lb(i1))==6 .And. lb(i2)==26) .Or. (lb(i1)==26 .And. iabs(lb(i2))==6)) Then
        If (iabs(lb(i1))==6) Then
           ii = i1
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 3
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 3
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(5) for D(+)+rho(-)-->pi(0)+n or pi(-)+p
     If (((lb(i1)==8 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==8)) .Or. ((lb(i1)==-8 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==-8))) Then
        If (iabs(lb(i1))==8) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(6) D(0)+rho(+)-->n+pi(+) or p+pi(0)
     If (((lb(i1)==7 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==7)) .Or. ((lb(i1)==-7 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==-7))) Then
        If (iabs(lb(i1))==7) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(7) for D(0)+rho(-)-->n+pi(-)
     If (((lb(i1)==7 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==7)) .Or. ((lb(i1)==-7 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==-7))) Then
        If (iabs(lb(i1))==7) Then
           ii = i1
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 3
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 3
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(8) D(-)+rho(+)-->n+pi(0) or p+pi(-)
     If (((lb(i1)==6 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==6)) .Or. ((lb(i1)==-6 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==-6))) Then
        If (iabs(lb(i1))==6) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(9) D(++)+rho(-)-->n+pi(+) or p+pi(0)
     If (((lb(i1)==9 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==9)) .Or. ((lb(i1)==-9 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==-9))) Then
        If (iabs(lb(i1))==9) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(10) for D(++)+rho(0)-->p+pi(+)
     If ((iabs(lb(i1))==9 .And. lb(i2)==26) .Or. (lb(i1)==26 .And. iabs(lb(i2))==9)) Then
        If (iabs(lb(i1))==9) Then
           ii = i1
           lb(i1) = 1
           e(i1) = amn
           lb(i2) = 5
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 1
           e(i2) = amn
           lb(i1) = 5
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(11) for N*(1440)(+)or N*(1535)(+)+rho(+)-->p+pion(+)
     If (((lb(i1)==11 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==11) .Or. (lb(i1)==13 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==13)) .Or. ((lb(i1)==-11 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==-11) .Or. (lb(i1)==-13 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==-13))) Then
        If (iabs(lb(i1))==11 .Or. iabs(lb(i1))==13) Then
           ii = i1
           lb(i1) = 1
           e(i1) = amn
           lb(i2) = 5
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 1
           e(i2) = amn
           lb(i1) = 5
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(12) for N*(1440) or N*(1535)(0)+rho(0)-->n+pi(0) or p+pi(-)
     If ((iabs(lb(i1))==10 .And. lb(i2)==26) .Or. (lb(i1)==26 .And. iabs(lb(i2))==10) .Or. (lb(i1)==26 .And. iabs(lb(i2))==12) .Or. (lb(i2)==26 .And. iabs(lb(i1))==12)) Then
        If (iabs(lb(i1))==10 .Or. iabs(lb(i1))==12) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(13) for N*(1440) or N*(1535)(+)+rho(0)-->pi(+)+n or pi(0)+p
     If ((iabs(lb(i1))==11 .And. lb(i2)==26) .Or. (lb(i1)==26 .And. iabs(lb(i2))==11) .Or. (lb(i1)==26 .And. iabs(lb(i2))==13) .Or. (lb(i2)==26 .And. iabs(lb(i1))==13)) Then
        If (iabs(lb(i1))==11 .Or. iabs(lb(i1))==13) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(14) for N*(1440) or N*(1535)(+)+rho(-)-->pi(0)+n or pi(-)+p
     If (((lb(i1)==11 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==11) .Or. (lb(i1)==25 .And. lb(i2)==13) .Or. (lb(i2)==25 .And. lb(i1)==13)) .Or. ((lb(i1)==-11 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==-11) .Or. (lb(i1)==27 .And. lb(i2)==-13) .Or. (lb(i2)==27 .And. lb(i1)==-13))) Then
        If (iabs(lb(i1))==11 .Or. iabs(lb(i1))==13) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(15) N*(1440) or N*(1535)(0)+rho(+)-->n+pi(+) or p+pi(0)
     If (((lb(i1)==10 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==10) .Or. (lb(i1)==12 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==12)) .Or. ((lb(i1)==-10 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==-10) .Or. (lb(i1)==-12 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==-12))) Then
        If (iabs(lb(i1))==10 .Or. iabs(lb(i1))==12) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(16) for N*(1440) or N*(1535) (0)+rho(-)-->n+pi(-)
     If (((lb(i1)==10 .And. lb(i2)==25) .Or. (lb(i1)==25 .And. lb(i2)==10) .Or. (lb(i1)==25 .And. lb(i2)==12) .Or. (lb(i1)==12 .And. lb(i2)==25)) .Or. ((lb(i1)==-10 .And. lb(i2)==27) .Or. (lb(i1)==27 .And. lb(i2)==-10) .Or. (lb(i1)==27 .And. lb(i2)==-12) .Or. (lb(i1)==-12 .And. lb(i2)==27))) Then
        If (iabs(lb(i1))==10 .Or. iabs(lb(i1))==12) Then
           ii = i1
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 3
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 3
           e(i1) = ap1
           Goto 40
        End If
     End If
60   iblock = 82
     ! FOR OMEGA REABSORPTION
     ! Relable particles, I1 is assigned to the Delta
     ! and I2 is assigned to the meson
     ! for the reverse of the following process
     !(1) for D(0)+OMEGA(0)-->n+pi(0) or p+pi(-)
     If ((iabs(lb(i1))==7 .And. lb(i2)==28) .Or. (lb(i1)==28 .And. iabs(lb(i2))==7)) Then
        If (iabs(lb(i1))==7) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(2) for D(+)+OMEGA(0)-->pi(+)+n or pi(0)+p
     If ((iabs(lb(i1))==8 .And. lb(i2)==28) .Or. (lb(i1)==28 .And. iabs(lb(i2))==8)) Then
        If (iabs(lb(i1))==8) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(3) for D(-)+OMEGA(0)-->n+pi(-)
     If ((iabs(lb(i1))==6 .And. lb(i2)==28) .Or. (lb(i1)==28 .And. iabs(lb(i2))==6)) Then
        If (iabs(lb(i1))==6) Then
           ii = i1
           lb(i1) = 2
           e(i1) = amn
           lb(i2) = 3
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 2
           e(i2) = amn
           lb(i1) = 3
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(4) for D(++)+OMEGA(0)-->p+pi(+)
     If ((iabs(lb(i1))==9 .And. lb(i2)==28) .Or. (lb(i1)==28 .And. iabs(lb(i2))==9)) Then
        If (iabs(lb(i1))==9) Then
           ii = i1
           lb(i1) = 1
           e(i1) = amn
           lb(i2) = 5
           e(i2) = ap1
           Goto 40
        Else
           ii = i2
           lb(i2) = 1
           e(i2) = amn
           lb(i1) = 5
           e(i1) = ap1
           Goto 40
        End If
     End If
     !(5) for N*(1440) or N*(1535)(0)+omega(0)-->n+pi(0) or p+pi(-)
     If ((iabs(lb(i1))==10 .And. lb(i2)==28) .Or. (lb(i1)==28 .And. iabs(lb(i2))==10) .Or. (lb(i1)==28 .And. iabs(lb(i2))==12) .Or. (lb(i2)==28 .And. iabs(lb(i1))==12)) Then
        If (iabs(lb(i1))==10 .Or. iabs(lb(i1))==12) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 3
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 3
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
     !(6) for N*(1440) or N*(1535)(+)+omega(0)-->pi(+)+n or pi(0)+p
     If ((iabs(lb(i1))==11 .And. lb(i2)==28) .Or. (lb(i1)==28 .And. iabs(lb(i2))==11) .Or. (lb(i1)==28 .And. iabs(lb(i2))==13) .Or. (lb(i2)==28 .And. iabs(lb(i1))==13)) Then
        If (iabs(lb(i1))==11 .Or. iabs(lb(i1))==13) Then
           ii = i1
           If (x2<=0.5) Then
              lb(i1) = 2
              e(i1) = amn
              lb(i2) = 5
              e(i2) = ap1
              Goto 40
           Else
              lb(i1) = 1
              e(i1) = amn
              lb(i2) = 4
              e(i2) = ap1
              Goto 40
           End If
        Else
           ii = i2
           If (x2<=0.5) Then
              lb(i2) = 2
              e(i2) = amn
              lb(i1) = 5
              e(i1) = ap1
              Goto 40
           Else
              lb(i2) = 1
              e(i2) = amn
              lb(i1) = 4
              e(i1) = ap1
              Goto 40
           End If
        End If
     End If
40   em1 = e(i1)
     em2 = e(i2)
     If (ianti==1 .And. lb(i1)>=1 .And. lb(i2)>=1) Then
        lb(ii) = -lb(ii)
        jj = i2
        If (ii==i2) jj = i1
        If (lb(jj)==3) Then
           lb(jj) = 5
        Else If (lb(jj)==5) Then
           lb(jj) = 3
        End If
     End If
  End If
  !-----------------------------------------------------------------------
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
50 pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  !          C1   = 1.0 - 2.0 * RANART(NSEED)

  !lin-10/25/02 get rid of argument usage mismatch in PTR():
  xptr = 0.33*pr
  !         cc1=ptr(0.33*pr,iseed)
  cc1 = ptr(xptr, iseed)
  !lin-10/25/02-end

  !lin-9/2012: check argument in sqrt():
  scheck = pr**2 - cc1**2
  If (scheck<0) Then
     Write (99, *) 'scheck39: ', scheck
     scheck = 0.
  End If
  c1 = sqrt(scheck)/pr
  !         c1=sqrt(pr**2-cc1**2)/pr

  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  ! ROTATE THE MOMENTUM
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crrd
!*********************************
! sp 03/19/01                                                          *
!                                                                      *
Subroutine crlaba(px, py, pz, srt, brel, brsgm, i1, i2, nt, iblock, nchrg, icase)
  !     PURPOSE:                                                         *
  !            DEALING WITH   K+ + N(D,N*)-bar <-->  La(Si)-bar + pi     *
  !     NOTE   :                                                         *
  !                                                                      *
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                     8-> elastic scatt                               *
  !                     100-> K+ + N-bar -> Sigma-bar + PI
  !                     102-> PI + Sigma(Lambda)-bar -> K+ + N-bar
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (etam=0.5475, aomega=0.782, arho=0.77)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save
  !
  px0 = px
  py0 = py
  pz0 = pz
  !
  If (icase==3) Then
     rrr = ranart(nseed)
     If (rrr<brel) Then
        !            !! elastic scat.  (avoid in reverse process)
        iblock = 8
     Else
        iblock = 100
        If (rrr<(brel+brsgm)) Then
           !*    K+ + N-bar -> Sigma-bar + PI
           lb(i1) = -15 - int(3*ranart(nseed))

           e(i1) = asa
        Else
           !*    K+ + N-bar -> Lambda-bar + PI
           lb(i1) = -14
           e(i1) = ala
        End If
        lb(i2) = 3 + int(3*ranart(nseed))
        e(i2) = 0.138
     End If
  End If
  !
  !
  If (icase==4) Then
     rrr = ranart(nseed)
     If (rrr<brel) Then
        !            !! elastic scat.
        iblock = 8
     Else
        iblock = 102
        !    PI + Sigma(Lambda)-bar -> K+ + N-bar
        !         ! K+
        lb(i1) = 23
        lb(i2) = -1 - int(2*ranart(nseed))
        If (nchrg==-2) lb(i2) = -6
        If (nchrg==1) lb(i2) = -9
        e(i1) = aka
        e(i2) = 0.938
        If (nchrg==-2 .Or. nchrg==1) e(i2) = 1.232
     End If
  End If
  !
  em1 = e(i1)
  em2 = e(i2)
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
  ! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  ! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crlaba
!*********************************
!                                                                      *
!                                                                      *
Subroutine crkn(px, py, pz, srt, i1, i2, iblock)
  !     PURPOSE:                                                         *
  !             DEALING WITH kaON+N/pi-->KAON +N/pi elastic PROCESS      *
  !     NOTE   :                                                         *
  !
  !     QUANTITIES:                                                 *
  !           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
  !           SRT      - SQRT OF S                                       *
  !           IBLOCK   - THE INFORMATION BACK                            *
  !                     8-> PION+N-->L/S+KAON
  !*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  !-----------------------------------------------------------------------
  iblock = 8
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  !-----------------------------------------------------------------------
  ! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  Return
End Subroutine crkn
!*********************************
!                                                                      *
!                                                                      *
Subroutine crppba(px, py, pz, srt, i1, i2, iblock)
!     PURPOSE:                                                         *

!lin-8/29/00*             DEALING WITH anti-nucleon annihilation with
!             DEALING WITH anti-baryon annihilation with

!             nucleons or baryon resonances
!             Determine:                                               *
!             (1) no. of pions in the final state
!             (2) relable particles in the final state
!             (3) new momenta of final state particles                 *
!
!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - INFORMATION about the reaction channel          *
!
!           iblock   - 1902 annihilation-->pion(+)+pion(-)   (2 pion)
!           iblock   - 1903 annihilation-->pion(+)+rho(-)    (3 pion)
!           iblock   - 1904 annihilation-->rho(+)+rho(-)     (4 pion)
!           iblock   - 1905 annihilation-->rho(0)+omega      (5 pion)
!           iblock   - 1906 annihilation-->omega+omega       (6 pion)
!       charge conservation is enforced in relabling particles
!       in the final state (note: at the momentum we don't check the
!       initial charges while dealing with annihilation, since some
!       annihilation channels between antinucleons and nucleons (baryon
!       resonances) might be forbiden by charge conservation, this effect
!       should be small, but keep it in mind.
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
! determine the no. of pions in the final state using a
! statistical model
  Call pbarfs(srt, npion, iseed)
! find the masses of the final state particles before calculate
! their momenta, and relable them. The masses of rho and omega
! will be generated according to the Breit Wigner formula       (NOTE!!!
! NOT DONE YET, AT THE MOMENT LET US USE FIXED RHO AND OMEGA MAEES)
!bali2/22/99
! Here we generate two stes of integer random numbers (3,4,5)
! one or both of them are used directly as the lables of pions
! similarly, 22+nchrg1 and 22+nchrg2 are used directly
! to label rhos
  nchrg1 = 3 + int(3*ranart(nseed))
  nchrg2 = 3 + int(3*ranart(nseed))
! the corresponding masses of pions
  pmass1 = ap1
  pmass2 = ap1
  If (nchrg1==3 .Or. nchrg1==5) pmass1 = ap2
  If (nchrg2==3 .Or. nchrg2==5) pmass2 = ap2
! (1) for 2 pion production
  If (npion==2) Then
    iblock = 1902
! randomly generate the charges of final state particles,
    lb(i1) = nchrg1
    e(i1) = pmass1
    lb(i2) = nchrg2
    e(i2) = pmass2
! TO CALCULATE THE FINAL MOMENTA
    Goto 50
  End If
! (2) FOR 3 PION PRODUCTION
  If (npion==3) Then
    iblock = 1903
    lb(i1) = nchrg1
    e(i1) = pmass1
    lb(i2) = 22 + nchrg2
    e(i2) = amrho
    Goto 50
  End If
! (3) FOR 4 PION PRODUCTION
! we allow both rho+rho and pi+omega with 50-50% probability
  If (npion==4) Then
    iblock = 1904
! determine rho+rho or pi+omega
    If (ranart(nseed)>=0.5) Then
! rho+rho
      lb(i1) = 22 + nchrg1
      e(i1) = amrho
      lb(i2) = 22 + nchrg2
      e(i2) = amrho
    Else
! pion+omega
      lb(i1) = nchrg1
      e(i1) = pmass1
      lb(i2) = 28
      e(i2) = amomga
    End If
    Goto 50
  End If
! (4) FOR 5 PION PRODUCTION
  If (npion==5) Then
    iblock = 1905
! RHO AND OMEGA
    lb(i1) = 22 + nchrg1
    e(i1) = amrho
    lb(i2) = 28
    e(i2) = amomga
    Goto 50
  End If
! (5) FOR 6 PION PRODUCTION
  If (npion==6) Then
    iblock = 1906
! OMEGA AND OMEGA
    lb(i1) = 28
    e(i1) = amomga
    lb(i2) = 28
    e(i2) = amomga
  End If
!bali2/22/99
  50 em1 = e(i1)
  em2 = e(i2)
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-08
  pr = sqrt(pr2)/(2.*srt)
! WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crppba
!bali2/7/99end
!bali3/5/99
!*********************************
!     PURPOSE:                                                         *
!     assign final states for K+K- --> light mesons
!
Subroutine crkkpi(i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigk, iblock, lbp1, lbp2, emm1, emm2)
!
!     QUANTITIES:                                                     *
!           IBLOCK   - INFORMATION about the reaction channel          *
!
!             iblock   - 1907
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, ameta=0.5473, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  iblock = 1907
  x1 = ranart(nseed)*sigk
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  xsk5 = xsk4 + xsk5
  xsk6 = xsk5 + xsk6
  xsk7 = xsk6 + xsk7
  xsk8 = xsk7 + xsk8
  xsk9 = xsk8 + xsk9
  xsk10 = xsk9 + xsk10
  If (x1<=xsk1) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 3 + int(3*ranart(nseed))
    e(i1) = ap2
    e(i2) = ap2
    Goto 100
  Else If (x1<=xsk2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = ap2
    e(i2) = amrho
    Goto 100
  Else If (x1<=xsk3) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 28
    e(i1) = ap2
    e(i2) = amomga
    Goto 100
  Else If (x1<=xsk4) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 0
    e(i1) = ap2
    e(i2) = ameta
    Goto 100
  Else If (x1<=xsk5) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = amrho
    e(i2) = amrho
    Goto 100
  Else If (x1<=xsk6) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 28
    e(i1) = amrho
    e(i2) = amomga
    Goto 100
  Else If (x1<=xsk7) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 0
    e(i1) = amrho
    e(i2) = ameta
    Goto 100
  Else If (x1<=xsk8) Then
    lb(i1) = 28
    lb(i2) = 28
    e(i1) = amomga
    e(i2) = amomga
    Goto 100
  Else If (x1<=xsk9) Then
    lb(i1) = 28
    lb(i2) = 0
    e(i1) = amomga
    e(i2) = ameta
    Goto 100
  Else If (x1<=xsk10) Then
    lb(i1) = 0
    lb(i2) = 0
    e(i1) = ameta
    e(i2) = ameta
  Else
    iblock = 222
    Call rhores(i1, i2)
!     !! phi
    lb(i1) = 29
!          return
    e(i2) = 0.
  End If

  100 Continue
  lbp1 = lb(i1)
  lbp2 = lb(i2)
  emm1 = e(i1)
  emm2 = e(i2)

  Return
End Subroutine crkkpi
!*********************************
!     PURPOSE:                                                         *
!             DEALING WITH K+Y -> piN scattering
!
Subroutine crkhyp(px, py, pz, srt, i1, i2, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk, ikmp, iblock)
!
!             Determine:                                               *
!             (1) relable particles in the final state                 *
!             (2) new momenta of final state particles                 *
!                                                                        *
!     QUANTITIES:                                                    *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - INFORMATION about the reaction channel          *
!                                                                     *
!             iblock   - 1908                                          *
!             iblock   - 222   !! phi                                  *
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, aphi=1.02, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (pimass=0.140, ameta=0.5473, aka=0.498, aml=1.116, ams=1.193, am1440=1.44, am1535=1.535)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  iblock = 1908
!
  x1 = ranart(nseed)*sigk
  xky2 = xky1 + xky2
  xky3 = xky2 + xky3
  xky4 = xky3 + xky4
  xky5 = xky4 + xky5
  xky6 = xky5 + xky6
  xky7 = xky6 + xky7
  xky8 = xky7 + xky8
  xky9 = xky8 + xky9
  xky10 = xky9 + xky10
  xky11 = xky10 + xky11
  xky12 = xky11 + xky12
  xky13 = xky12 + xky13
  xky14 = xky13 + xky14
  xky15 = xky14 + xky15
  xky16 = xky15 + xky16
  If (x1<=xky1) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = pimass
    e(i2) = amp
    Goto 100
  Else If (x1<=xky2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = pimass
    e(i2) = am0
    Goto 100
  Else If (x1<=xky3) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = pimass
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky4) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = pimass
    e(i2) = am1535
    Goto 100
  Else If (x1<=xky5) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = amrho
    e(i2) = amp
    Goto 100
  Else If (x1<=xky6) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = amrho
    e(i2) = am0
    Goto 100
  Else If (x1<=xky7) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = amrho
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky8) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = amrho
    e(i2) = am1535
    Goto 100
  Else If (x1<=xky9) Then
    lb(i1) = 28
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = amomga
    e(i2) = amp
    Goto 100
  Else If (x1<=xky10) Then
    lb(i1) = 28
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = amomga
    e(i2) = am0
    Goto 100
  Else If (x1<=xky11) Then
    lb(i1) = 28
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = amomga
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky12) Then
    lb(i1) = 28
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = amomga
    e(i2) = am1535
    Goto 100
  Else If (x1<=xky13) Then
    lb(i1) = 0
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = ameta
    e(i2) = amp
    Goto 100
  Else If (x1<=xky14) Then
    lb(i1) = 0
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = ameta
    e(i2) = am0
    Goto 100
  Else If (x1<=xky15) Then
    lb(i1) = 0
    lb(i2) = 10 + int(2*ranart(nseed))
    e(i1) = ameta
    e(i2) = am1440
    Goto 100
  Else If (x1<=xky16) Then
    lb(i1) = 0
    lb(i2) = 12 + int(2*ranart(nseed))
    e(i1) = ameta
    e(i2) = am1535
    Goto 100
  Else
    lb(i1) = 29
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = aphi
    e(i2) = amn
    iblock = 222
    Goto 100
  End If

  100 Continue
  If (ikmp==-1) lb(i2) = -lb(i2)

  em1 = e(i1)
  em2 = e(i2)
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-08
  pr = sqrt(pr2)/(2.*srt)
! WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crkhyp
!*********************************
!                                                                      *
!                                                                      *


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Subroutine crlan(px, py, pz, srt, i1, i2, iblock)
!     PURPOSE:                                                         *
!      DEALING WITH La/Si-bar + N --> K+ + pi PROCESS                  *
!                   La/Si + N-bar --> K- + pi                          *
!     NOTE   :                                                         *
!
!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - THE INFORMATION BACK                            *
!                      71
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
  iblock = 71
  ntag = 0
  If ((lb(i1)>=14 .And. lb(i1)<=17) .Or. (lb(i2)>=14 .And. lb(i2)<=17)) Then
    lb(i1) = 21
  Else
    lb(i1) = 23
  End If
  lb(i2) = 3 + int(3*ranart(nseed))
  e(i1) = aka
  e(i2) = 0.138
  em1 = e(i1)
  em2 = e(i2)
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE
  Return
End Subroutine crlan
!sp11/03/01 end
!*********************************
!*********************************
!                                                                      *
!                                                                      *
Subroutine crkpla(px, py, pz, ec, srt, spika, emm1, emm2, lbp1, lbp2, i1, i2, icase, srhoks)

!     PURPOSE:                                                         *
!     DEALING WITH  K+ + Pi ---> La/Si-bar + B, phi+K, phi+K* OR  K* *
!                   K- + Pi ---> La/Si + B-bar  OR   K*-bar          *

!     NOTE   :                                                         *
!
!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - THE INFORMATION BACK                            *
!                      71
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, amrho=0.769, amomga=0.782, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895, ala=1.1157, asa=1.1974, aphi=1.02)
  Parameter (am1440=1.44, am1535=1.535)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  emm1 = 0.
  emm2 = 0.
  lbp1 = 0
  lbp2 = 0
  xkp0 = spika
  xkp1 = 0.
  xkp2 = 0.
  xkp3 = 0.
  xkp4 = 0.
  xkp5 = 0.
  xkp6 = 0.
  xkp7 = 0.
  xkp8 = 0.
  xkp9 = 0.
  xkp10 = 0.
  sigm = 15.
!         if(lb(i1).eq.21.or.lb(i2).eq.21)sigm=10.
  pdd = (srt**2-(aka+ap1)**2)*(srt**2-(aka-ap1)**2)
!
  If (srt<(ala+amn)) Goto 70
  xkp1 = sigm*(4./3.)*(srt**2-(ala+amn)**2)*(srt**2-(ala-amn)**2)/pdd
  If (srt>(ala+am0)) Then
    xkp2 = sigm*(16./3.)*(srt**2-(ala+am0)**2)*(srt**2-(ala-am0)**2)/pdd
  End If
  If (srt>(ala+am1440)) Then
    xkp3 = sigm*(4./3.)*(srt**2-(ala+am1440)**2)*(srt**2-(ala-am1440)**2)/pdd
  End If
  If (srt>(ala+am1535)) Then
    xkp4 = sigm*(4./3.)*(srt**2-(ala+am1535)**2)*(srt**2-(ala-am1535)**2)/pdd
  End If
!
  If (srt>(asa+amn)) Then
    xkp5 = sigm*4.*(srt**2-(asa+amn)**2)*(srt**2-(asa-amn)**2)/pdd
  End If
  If (srt>(asa+am0)) Then
    xkp6 = sigm*16.*(srt**2-(asa+am0)**2)*(srt**2-(asa-am0)**2)/pdd
  End If
  If (srt>(asa+am1440)) Then
    xkp7 = sigm*4.*(srt**2-(asa+am1440)**2)*(srt**2-(asa-am1440)**2)/pdd
  End If
  If (srt>(asa+am1535)) Then
    xkp8 = sigm*4.*(srt**2-(asa+am1535)**2)*(srt**2-(asa-am1535)**2)/pdd
  End If
  70 Continue
  sig1 = 195.639
  sig2 = 372.378
  If (srt>aphi+aka) Then
    pff = sqrt((srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2))

!lin-9/2012: check argument in sqrt():
    scheck = pdd
    If (scheck<=0) Then
      Write (99, *) 'scheck40: ', scheck
      Stop
    End If

    xkp9 = sig1*pff/sqrt(pdd)*1./32./pi/srt**2
    If (srt>aphi+aks) Then
      pff = sqrt((srt**2-(aphi+aks)**2)*(srt**2-(aphi-aks)**2))

!lin-9/2012: check argument in sqrt():
      scheck = pdd
      If (scheck<=0) Then
        Write (99, *) 'scheck41: ', scheck
        Stop
      End If

      xkp10 = sig2*pff/sqrt(pdd)*3./32./pi/srt**2
    End If
  End If

!lin-8/15/02 K pi -> K* (rho omega), from detailed balance,
! neglect rho and omega mass difference for now:
  sigpik = 0.
  If (srt>(amrho+aks)) Then
    sigpik = srhoks*9.*(srt**2-(0.77-aks)**2)*(srt**2-(0.77+aks)**2)/4/srt**2/(px**2+py**2+pz**2)
    If (srt>(amomga+aks)) sigpik = sigpik*12./9.
  End If

!
  sigkp = xkp0 + xkp1 + xkp2 + xkp3 + xkp4 + xkp5 + xkp6 + xkp7 + xkp8 + xkp9 + xkp10 + sigpik
  icase = 0
  dskn = sqrt(sigkp/pi/10.)
  dsknr = dskn + 0.1
  Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
!
  randu = ranart(nseed)*sigkp
  xkp1 = xkp0 + xkp1
  xkp2 = xkp1 + xkp2
  xkp3 = xkp2 + xkp3
  xkp4 = xkp3 + xkp4
  xkp5 = xkp4 + xkp5
  xkp6 = xkp5 + xkp6
  xkp7 = xkp6 + xkp7
  xkp8 = xkp7 + xkp8
  xkp9 = xkp8 + xkp9

  xkp10 = xkp9 + xkp10
!
!   !! K* formation
  If (randu<=xkp0) Then
    icase = 1
    Return
  Else
! La/Si-bar + B formation
    icase = 2
    If (randu<=xkp1) Then
      lbp1 = -14
      lbp2 = 1 + int(2*ranart(nseed))
      emm1 = ala
      emm2 = amn
      Goto 60
    Else If (randu<=xkp2) Then
      lbp1 = -14
      lbp2 = 6 + int(4*ranart(nseed))
      emm1 = ala
      emm2 = am0
      Goto 60
    Else If (randu<=xkp3) Then
      lbp1 = -14
      lbp2 = 10 + int(2*ranart(nseed))
      emm1 = ala
      emm2 = am1440
      Goto 60
    Else If (randu<=xkp4) Then
      lbp1 = -14
      lbp2 = 12 + int(2*ranart(nseed))
      emm1 = ala
      emm2 = am1535
      Goto 60
    Else If (randu<=xkp5) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 1 + int(2*ranart(nseed))
      emm1 = asa
      emm2 = amn
      Goto 60
    Else If (randu<=xkp6) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 6 + int(4*ranart(nseed))
      emm1 = asa
      emm2 = am0
      Goto 60
    Else If (randu<xkp7) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 10 + int(2*ranart(nseed))
      emm1 = asa
      emm2 = am1440
      Goto 60
    Else If (randu<xkp8) Then
      lbp1 = -15 - int(3*ranart(nseed))
      lbp2 = 12 + int(2*ranart(nseed))
      emm1 = asa
      emm2 = am1535
      Goto 60
    Else If (randu<xkp9) Then
!       !! phi +K  formation (iblock=224)
      icase = 3
      lbp1 = 29
      lbp2 = 23
      emm1 = aphi
      emm2 = aka
      If (lb(i1)==21 .Or. lb(i2)==21) Then
!         !! phi +K-bar  formation (iblock=124)
        lbp2 = 21
        icase = -3
      End If
      Goto 60
    Else If (randu<xkp10) Then
!       !! phi +K* formation (iblock=226)
      icase = 4
      lbp1 = 29
      lbp2 = 30
      emm1 = aphi
      emm2 = aks
      If (lb(i1)==21 .Or. lb(i2)==21) Then
        lbp2 = -30
        icase = -4
      End If
      Goto 60

    Else
!       !! (rho,omega) +K* formation (iblock=88)
      icase = 5
      lbp1 = 25 + int(3*ranart(nseed))
      lbp2 = 30
      emm1 = amrho
      emm2 = aks
      If (srt>(amomga+aks) .And. ranart(nseed)<0.25) Then
        lbp1 = 28
        emm1 = amomga
      End If
      If (lb(i1)==21 .Or. lb(i2)==21) Then
        lbp2 = -30
        icase = -5
      End If

    End If
  End If
!
  60 If (icase==2 .And. (lb(i1)==21 .Or. lb(i2)==21)) Then
    lbp1 = -lbp1
    lbp2 = -lbp2
  End If
  px0 = px
  py0 = py
  pz0 = pz
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-emm1**2-emm2**2)**2 - 4.0*(emm1*emm2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE
  Return
End Subroutine crkpla
!*********************************
!                                                                      *
!                                                                      *
Subroutine crkphi(px, py, pz, ec, srt, iblock, emm1, emm2, lbp1, lbp2, i1, i2, ikk, icase, rrkk, prkk)

!     PURPOSE:                                                         *
!     DEALING WITH   KKbar, KK*bar, KbarK*, K*K*bar --> Phi + pi(rho,omega)
!     and KKbar --> (pi eta) (pi eta), (rho omega) (rho omega)
!     and KK*bar or Kbar K* --> (pi eta) (rho omega)
!
!     NOTE   :                                                         *
!
!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - THE INFORMATION BACK                            *
!                      222
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, aphi=1.02, am0=1.232, amns=1.52, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, acas=1.3213)
  Parameter (aks=0.895, aomega=0.7819, arho=0.77)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  lb1 = lb(i1)
  lb2 = lb(i2)
  icase = 0

!        if(srt .lt. aphi+ap1)return
!c        if(srt .lt. aphi+ap1) then
  If (srt<(aphi+ap1)) Then
    sig1 = 0.
    sig2 = 0.
    sig3 = 0.
  Else
!
    If ((lb1==23 .And. lb2==21) .Or. (lb2==23 .And. lb1==21)) Then
      dnr = 4.
      ikk = 2
    Else If ((lb1==21 .And. lb2==30) .Or. (lb2==21 .And. lb1==30) .Or. (lb1==23 .And. lb2==-30) .Or. (lb2==23 .And. lb1==-30)) Then
      dnr = 12.
      ikk = 1
    Else
      dnr = 36.
      ikk = 0
    End If

    sig1 = 0.
    sig2 = 0.
    sig3 = 0.
    srri = e(i1) + e(i2)
    srr1 = aphi + ap1
    srr2 = aphi + aomega
    srr3 = aphi + arho
!
    pii = (srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2)
    srrt = srt - amax1(srri, srr1)
!c   to avoid divergent/negative values at small srrt:
!          if(srrt .lt. 0.3)then
    If (srrt<0.3 .And. srrt>0.01) Then
      sig = 1.69/(srrt**0.141-0.407)
    Else
      sig = 3.74 + 0.008*srrt**1.9
    End If
    sig1 = sig*(9./dnr)*(srt**2-(aphi+ap1)**2)*(srt**2-(aphi-ap1)**2)/pii
    If (srt>aphi+aomega) Then
      srrt = srt - amax1(srri, srr2)
!c         if(srrt .lt. 0.3)then
      If (srrt<0.3 .And. srrt>0.01) Then
        sig = 1.69/(srrt**0.141-0.407)
      Else
        sig = 3.74 + 0.008*srrt**1.9
      End If
      sig2 = sig*(9./dnr)*(srt**2-(aphi+aomega)**2)*(srt**2-(aphi-aomega)**2)/pii
    End If
    If (srt>aphi+arho) Then
      srrt = srt - amax1(srri, srr3)
!c         if(srrt .lt. 0.3)then
      If (srrt<0.3 .And. srrt>0.01) Then
        sig = 1.69/(srrt**0.141-0.407)
      Else
        sig = 3.74 + 0.008*srrt**1.9
      End If
      sig3 = sig*(27./dnr)*(srt**2-(aphi+arho)**2)*(srt**2-(aphi-arho)**2)/pii
    End If
!         sig1 = amin1(20.,sig1)
!         sig2 = amin1(20.,sig2)
!         sig3 = amin1(20.,sig3)
  End If

  rrkk0 = rrkk
  prkk0 = prkk
  sigm = 0.
  If ((lb1==23 .And. lb2==21) .Or. (lb2==23 .And. lb1==21)) Then
    Call xkkann(srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigm, rrkk0)
  Else If ((lb1==21 .And. lb2==30) .Or. (lb2==21 .And. lb1==30) .Or. (lb1==23 .And. lb2==-30) .Or. (lb2==23 .And. lb1==-30)) Then
    Call xkksan(i1, i2, srt, sigks1, sigks2, sigks3, sigks4, sigm, prkk0)
  Else
  End If
!
!         sigks = sig1 + sig2 + sig3
  sigm0 = sigm
  sigks = sig1 + sig2 + sig3 + sigm
  dskn = sqrt(sigks/pi/10.)
  dsknr = dskn + 0.1
  Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  icase = 1
  ranx = ranart(nseed)

  lbp1 = 29
  emm1 = aphi
  If (ranx<=sig1/sigks) Then
    lbp2 = 3 + int(3*ranart(nseed))
    emm2 = ap1
  Else If (ranx<=(sig1+sig2)/sigks) Then
    lbp2 = 28
    emm2 = aomega
  Else If (ranx<=(sig1+sig2+sig3)/sigks) Then
    lbp2 = 25 + int(3*ranart(nseed))
    emm2 = arho
  Else
    If ((lb1==23 .And. lb2==21) .Or. (lb2==23 .And. lb1==21)) Then
      Call crkkpi(i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigm0, iblock, lbp1, lbp2, emm1, emm2)
    Else If ((lb1==21 .And. lb2==30) .Or. (lb2==21 .And. lb1==30) .Or. (lb1==23 .And. lb2==-30) .Or. (lb2==23 .And. lb1==-30)) Then
      Call crkspi(i1, i2, sigks1, sigks2, sigks3, sigks4, sigm0, iblock, lbp1, lbp2, emm1, emm2)
    Else
    End If
  End If
!
  px0 = px
  py0 = py
  pz0 = pz
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-emm1**2-emm2**2)**2 - 4.0*(emm1*emm2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE
  Return
End Subroutine crkphi
!sp11/21/01 end
!*********************************
!                                                                      *
!                                                                      *
Subroutine crksph(px, py, pz, ec, srt, emm1, emm2, lbp1, lbp2, i1, i2, ikkg, ikkl, iblock, icase, srhoks)

!     PURPOSE:                                                         *
!     DEALING WITH   K + rho(omega) or K* + pi(rho,omega)
!                    --> Phi + K(K*), pi + K* or pi + K, and elastic
!     NOTE   :                                                         *
!
!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - THE INFORMATION BACK                            *
!                      222
!                      223 --> phi + pi(rho,omega)
!                      224 --> phi + K <-> K + pi(rho,omega)
!                      225 --> phi + K <-> K* + pi(rho,omega)
!                      226 --> phi + K* <-> K + pi(rho,omega)
!                      227 --> phi + K* <-> K* + pi(rho,omega)
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, aphi=1.02, am0=1.232, amns=1.52, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, acas=1.3213)
  Parameter (aks=0.895, aomega=0.7819, arho=0.77)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  lb1 = lb(i1)
  lb2 = lb(i2)
  icase = 0
  sigela = 10.
  sigkm = 0.
!     K(K*) + rho(omega) -> pi K*(K)
  If ((lb1>=25 .And. lb1<=28) .Or. (lb2>=25 .And. lb2<=28)) Then
    If (iabs(lb1)==30 .Or. iabs(lb2)==30) Then
      sigkm = srhoks
!lin-2/26/03 check whether (rho K) is above the (pi K*) thresh:
    Else If ((lb1==23 .Or. lb1==21 .Or. lb2==23 .Or. lb2==21) .And. srt>(ap2+aks)) Then
      sigkm = srhoks
    End If
  End If

!        if(srt .lt. aphi+aka)return
  If (srt<(aphi+aka)) Then
    sig11 = 0.
    sig22 = 0.
  Else

! K*-bar +pi --> phi + (K,K*)-bar
    If ((iabs(lb1)==30 .And. (lb2>=3 .And. lb2<=5)) .Or. (iabs(lb2)==30 .And. (lb1>=3 .And. lb1<=5))) Then
      dnr = 18.
      ikkl = 0
      iblock = 225
!               sig1 = 15.0
!               sig2 = 30.0
!lin-2/06/03 these large values reduces to ~10 mb for sig11 or sig22
!     due to the factors of ~1/(32*pi*s)~1/200:
      sig1 = 2047.042
      sig2 = 1496.692
! K(-bar)+rho --> phi + (K,K*)-bar
    Else If ((lb1==23 .Or. lb1==21 .And. (lb2>=25 .And. lb2<=27)) .Or. (lb2==23 .Or. lb2==21 .And. (lb1>=25 .And. lb1<=27))) Then
      dnr = 18.
      ikkl = 1
      iblock = 224
!               sig1 = 3.5
!               sig2 = 9.0
      sig1 = 526.702
      sig2 = 1313.960
! K*(-bar) +rho
    Else If ((iabs(lb1)==30 .And. (lb2>=25 .And. lb2<=27)) .Or. (iabs(lb2)==30 .And. (lb1>=25 .And. lb1<=27))) Then
      dnr = 54.
      ikkl = 0
      iblock = 225
!               sig1 = 3.5
!               sig2 = 9.0
      sig1 = 1371.257
      sig2 = 6999.840
! K(-bar) + omega
    Else If (((lb1==23 .Or. lb1==21) .And. lb2==28) .Or. ((lb2==23 .Or. lb2==21) .And. lb1==28)) Then
      dnr = 6.
      ikkl = 1
      iblock = 224
!               sig1 = 3.5
!               sig2 = 6.5
      sig1 = 355.429
      sig2 = 440.558
! K*(-bar) +omega
    Else
      dnr = 18.
      ikkl = 0
      iblock = 225
!               sig1 = 3.5
!               sig2 = 15.0
      sig1 = 482.292
      sig2 = 1698.903
    End If

    sig11 = 0.
    sig22 = 0.
!         sig11=sig1*(6./dnr)*(srt**2-(aphi+aka)**2)*
!    &           (srt**2-(aphi-aka)**2)/(srt**2-(e(i1)+e(i2))**2)/
!    &           (srt**2-(e(i1)-e(i2))**2)

!lin-9/2012: check argument in sqrt():
    scheck = (srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2)
    If (scheck<=0) Then
      Write (99, *) 'scheck42: ', scheck
      Stop
    End If
    pii = sqrt(scheck)
!        pii = sqrt((srt**2-(e(i1)+e(i2))**2)*(srt**2-(e(i1)-e(i2))**2))

!lin-9/2012: check argument in sqrt():
    scheck = (srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2)
    If (scheck<0) Then
      Write (99, *) 'scheck43: ', scheck
      scheck = 0.
    End If
    pff = sqrt(scheck)
!        pff = sqrt((srt**2-(aphi+aka)**2)*(srt**2-(aphi-aka)**2))

    sig11 = sig1*pff/pii*6./dnr/32./pi/srt**2
!
    If (srt>aphi+aks) Then
!         sig22=sig2*(18./dnr)*(srt**2-(aphi+aks)**2)*
!    &           (srt**2-(aphi-aks)**2)/(srt**2-(e(i1)+e(i2))**2)/
!    &           (srt**2-(e(i1)-e(i2))**2)
      pff = sqrt((srt**2-(aphi+aks)**2)*(srt**2-(aphi-aks)**2))
      sig22 = sig2*pff/pii*18./dnr/32./pi/srt**2
    End If
!         sig11 = amin1(20.,sig11)
!         sig22 = amin1(20.,sig22)
!
  End If

!         sigks = sig11 + sig22
  sigks = sig11 + sig22 + sigela + sigkm
!
  dskn = sqrt(sigks/pi/10.)
  dsknr = dskn + 0.1
  Call distce(i1, i2, dsknr, dskn, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  icase = 1
  ranx = ranart(nseed)

  If (ranx<=(sigela/sigks)) Then
    lbp1 = lb1
    emm1 = e(i1)
    lbp2 = lb2
    emm2 = e(i2)
    iblock = 111
  Else If (ranx<=((sigela+sigkm)/sigks)) Then
    lbp1 = 3 + int(3*ranart(nseed))
    emm1 = 0.14
    If (lb1==23 .Or. lb2==23) Then
      lbp2 = 30
      emm2 = aks
    Else If (lb1==21 .Or. lb2==21) Then
      lbp2 = -30
      emm2 = aks
    Else If (lb1==30 .Or. lb2==30) Then
      lbp2 = 23
      emm2 = aka
    Else
      lbp2 = 21
      emm2 = aka
    End If
    iblock = 112
  Else If (ranx<=((sigela+sigkm+sig11)/sigks)) Then
    lbp2 = 23
    emm2 = aka
    ikkg = 1
    If (lb1==21 .Or. lb2==21 .Or. lb1==-30 .Or. lb2==-30) Then
      lbp2 = 21
      iblock = iblock - 100
    End If
    lbp1 = 29
    emm1 = aphi
  Else
    lbp2 = 30
    emm2 = aks
    ikkg = 0
    iblock = iblock + 2
    If (lb1==21 .Or. lb2==21 .Or. lb1==-30 .Or. lb2==-30) Then
      lbp2 = -30
      iblock = iblock - 100
    End If
    lbp1 = 29
    emm1 = aphi
  End If
!
  px0 = px
  py0 = py
  pz0 = pz
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-emm1**2-emm2**2)**2 - 4.0*(emm1*emm2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! FOR THE ISOTROPIC DISTRIBUTION THERE IS NO NEED TO ROTATE
  Return
End Subroutine crksph
!sp11/21/01 end
!*********************************
!*********************************
Subroutine bbkaon(ic, srt, px, py, pz, ana, plx, ply, plz, ala, pkx, pky, pkz, icou1)
! purpose: generate the momenta for kaon,lambda/sigma and nucleon/delta
!          in the BB-->nlk process
! date: Sept. 9, 1994
!
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  pi = 3.1415962
  icou1 = 0
  aka = 0.498
  ala = 1.116
  If (ic==2 .Or. ic==4) ala = 1.197
  ana = 0.939
! generate the mass of the delta
  If (ic>2) Then
    dmax = srt - aka - ala - 0.02
    dm1 = rmass(dmax, iseed)
    ana = dm1
  End If
  t1 = aka + ana + ala
  t2 = ana + ala - aka
  If (srt<=t1) Then
    icou1 = -1
    Return
  End If
  pmax = sqrt((srt**2-t1**2)*(srt**2-t2**2))/(2.*srt)
  If (pmax==0.) pmax = 1.E-09
! (1) Generate the momentum of the kaon according to the distribution Fkaon
!     and assume that the angular distribution is isotropic
!     in the cms of the colliding pair
  ntry = 0
  1 pk = pmax*ranart(nseed)
  ntry = ntry + 1
  prob = fkaon(pk, pmax)
  If ((prob<ranart(nseed)) .And. (ntry<=40)) Goto 1
  cs = 1. - 2.*ranart(nseed)
  ss = sqrt(1.-cs**2)
  fai = 2.*3.14*ranart(nseed)
  pkx = pk*ss*cos(fai)
  pky = pk*ss*sin(fai)
  pkz = pk*cs
! the energy of the kaon
  ek = sqrt(aka**2+pk**2)
! (2) Generate the momentum of the nucleon/delta in the cms of N/delta
!     and lamda/sigma
!  the energy of the cms of NL
  eln = srt - ek
  If (eln<=0) Then
    icou1 = -1
    Return
  End If
! beta and gamma of the cms of L/S+N
  bx = -pkx/eln
  by = -pky/eln
  bz = -pkz/eln

!lin-9/2012: check argument in sqrt():
  scheck = 1. - bx**2 - by**2 - bz**2
  If (scheck<=0) Then
    Write (99, *) 'scheck44: ', scheck
    Stop
  End If
  ga = 1./sqrt(scheck)
!       ga=1./sqrt(1.-bx**2-by**2-bz**2)

  elnc = eln/ga
  pn2 = ((elnc**2+ana**2-ala**2)/(2.*elnc))**2 - ana**2
  If (pn2<=0.) pn2 = 1.E-09
  pn = sqrt(pn2)
  csn = 1. - 2.*ranart(nseed)
  ssn = sqrt(1.-csn**2)
  fain = 2.*3.14*ranart(nseed)
  px = pn*ssn*cos(fain)
  py = pn*ssn*sin(fain)
  pz = pn*csn
  en = sqrt(ana**2+pn2)
! the momentum of the lambda/sigma in the n-l cms frame is
  plx = -px
  ply = -py
  plz = -pz
! (3) LORENTZ-TRANSFORMATION INTO nn cms FRAME for the neutron/delta
  pbeta = px*bx + py*by + pz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+en)
  px = bx*trans0 + px
  py = by*trans0 + py
  pz = bz*trans0 + pz
! (4) Lorentz-transformation for the lambda/sigma
  el = sqrt(ala**2+plx**2+ply**2+plz**2)
  pbeta = plx*bx + ply*by + plz*bz
  trans0 = ga*(ga*pbeta/(ga+1.)+el)
  plx = bx*trans0 + plx
  ply = by*trans0 + ply
  plz = bz*trans0 + plz
  Return
End Subroutine bbkaon
!*****************************************
! for pion+pion-->K+K-
!      real*4 function pipik(srt)
Real Function pipik(srt)
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in mb                                    *
!  NOTE: DEVIDE THE CROSS SECTION TO OBTAIN K+ PRODUCTION                     *
!*****************************************
!      real*4   xarray(5), earray(5)
  Real xarray(5), earray(5)
  Save
  Data xarray/0.001, 0.7, 1.5, 1.7, 2.0/
  Data earray/1., 1.2, 1.6, 2.0, 2.4/

  pmass = 0.9383
! 1.Calculate p(lab)  from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
!      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  pipik = 0.
  If (srt<=1.) Return
  If (srt>2.4) Then
    pipik = 2.0/2.
    Return
  End If
  If (srt<earray(1)) Then
    pipik = xarray(1)/2.
    Return
  End If
!
! 2.Interpolate double logarithmically to find sigma(srt)
!
  Do ie = 1, 5
    If (earray(ie)==srt) Then
      pipik = xarray(ie)
      Goto 10
    Else If (earray(ie)>srt) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      pipik = exp(ymin+(alog(srt)-xmin)*(ymax-ymin)/(xmax-xmin))
      Goto 10
    End If
  End Do
  10 pipik = pipik/2.
  Continue
  Return
End Function pipik
!*********************************
! TOTAL PION-P INELASTIC CROSS SECTION
!  from the CERN data book
!  date: Sept.2, 1994
!  for pion++p-->Delta+pion
!      real*4 function pionpp(srt)
Real Function pionpp(srt)
  Save
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in fm**2                                 *
!  earray = EXPerimental table with proton energies in MeV                    *
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!                                                                             *
!*****************************************
  pmass = 0.14
  pmass1 = 0.938
  pionpp = 0.00001
  If (srt<=1.22) Return
! 1.Calculate p(lab)  from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
!      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  plab = sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
  pmin = 0.3
  pmax = 25.0
  If (plab>pmax) Then
    pionpp = 20./10.
    Return
  End If
  If (plab<pmin) Then
    pionpp = 0.
    Return
  End If
!* fit parameters
  a = 24.3
  b = -12.3
  c = 0.324
  an = -1.91
  d = -2.44
  pionpp = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (pionpp<=0) pionpp = 0
  pionpp = pionpp/10.
  Return
End Function pionpp
!*********************************
! elementary cross sections
!  from the CERN data book
!  date: Sept.2, 1994
!  for pion-+p-->INELASTIC
!      real*4 function pipp1(srt)
Real Function pipp1(srt)
  Save
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in fm**2                                 *
!  earray = EXPerimental table with proton energies in MeV                    *
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!  UNITS: FM**2
!*****************************************
  pmass = 0.14
  pmass1 = 0.938
  pipp1 = 0.0001
  If (srt<=1.22) Return
! 1.Calculate p(lab)  from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
!      ekin = 2.*pmass*((srt/(2.*pmass))**2 - 1.)
  plab = sqrt(((srt**2-pmass**2-pmass1**2)/(2.*pmass1))**2-pmass**2)
  pmin = 0.3
  pmax = 25.0
  If (plab>pmax) Then
    pipp1 = 20./10.
    Return
  End If
  If (plab<pmin) Then
    pipp1 = 0.
    Return
  End If
!* fit parameters
  a = 26.6
  b = -7.18
  c = 0.327
  an = -1.86
  d = -2.81
  pipp1 = a + b*(plab**an) + c*(alog(plab))**2 + d*alog(plab)
  If (pipp1<=0) pipp1 = 0
  pipp1 = pipp1/10.
  Return
End Function pipp1
! *****************************
!       real*4 function xrho(srt)
Real Function xrho(srt)
  Save
!       xsection for pp-->pp+rho
! *****************************
  pmass = 0.9383
  rmass = 0.77
  trho = 0.151
  xrho = 0.000000001
  If (srt<=2.67) Return
  esmin = 2.*0.9383 + rmass - trho/2.
  es = srt
! the cross section for tho0 production is
  xrho0 = 0.24*(es-esmin)/(1.4+(es-esmin)**2)
  xrho = 3.*xrho0
  Return
End Function xrho
! *****************************
!       real*4 function omega(srt)
Real Function omega(srt)
  Save
!       xsection for pp-->pp+omega
! *****************************
  pmass = 0.9383
  omass = 0.782
  tomega = 0.0084
  omega = 0.00000001
  If (srt<=2.68) Return
  esmin = 2.*0.9383 + omass - tomega/2.
  es = srt
  omega = 0.36*(es-esmin)/(1.25+(es-esmin)**2)
  Return
End Function omega
!*****************************************
! for ppi(+)-->DELTA+pi
!      real*4 function TWOPI(srt)
Real Function twopi(srt)
!  This function contains the experimental pi+p-->DELTA+PION cross sections   *
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in mb                                    *
!  earray = EXPerimental table with proton energies in MeV                    *
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!                                                                             *
!*****************************************
!      real*4   xarray(19), earray(19)
  Real xarray(19), earray(19)
  Save
  Data xarray/0.300E-05, 0.187E+01, 0.110E+02, 0.149E+02, 0.935E+01, 0.765E+01, 0.462E+01, 0.345E+01, 0.241E+01, 0.185E+01, 0.165E+01, 0.150E+01, 0.132E+01, 0.117E+01, 0.116E+01, 0.100E+01, 0.856E+00, 0.745E+00, 0.300E-05/
  Data earray/0.122E+01, 0.147E+01, 0.172E+01, 0.197E+01, 0.222E+01, 0.247E+01, 0.272E+01, 0.297E+01, 0.322E+01, 0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01, 0.472E+01, 0.497E+01, 0.522E+01, 0.547E+01, 0.572E+01/

  pmass = 0.14
  pmass1 = 0.938
  twopi = 0.000001
  If (srt<=1.22) Return
! 1.Calculate p(lab)  from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  plab = srt
  If (plab<earray(1)) Then
    twopi = 0.00001
    Return
  End If
!
! 2.Interpolate double logarithmically to find sigma(srt)
!
  Do ie = 1, 19
    If (earray(ie)==plab) Then
      twopi = xarray(ie)
      Return
    Else If (earray(ie)>plab) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      twopi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
      Return
    End If
  End Do
  Return
End Function twopi
!*****************************************
!*****************************************
! for ppi(+)-->DELTA+RHO
!      real*4 function THREPI(srt)
Real Function threpi(srt)
!  This function contains the experimental pi+p-->DELTA + rho cross sections  *
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in mb                                    *
!  earray = EXPerimental table with proton energies in MeV                    *
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!                                                                             *
!*****************************************
!      real*4   xarray(15), earray(15)
  Real xarray(15), earray(15)
  Save
  Data xarray/8.0000000E-06, 6.1999999E-05, 1.881940, 5.025690, 11.80154, 13.92114, 15.07308, 11.79571, 11.53772, 10.01197, 9.792673, 9.465264, 8.970490, 7.944254, 6.886320/
  Data earray/0.122E+01, 0.147E+01, 0.172E+01, 0.197E+01, 0.222E+01, 0.247E+01, 0.272E+01, 0.297E+01, 0.322E+01, 0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01, 0.472E+01/

  pmass = 0.14
  pmass1 = 0.938
  threpi = 0.000001
  If (srt<=1.36) Return
! 1.Calculate p(lab)  from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  plab = srt
  If (plab<earray(1)) Then
    threpi = 0.00001
    Return
  End If
!
! 2.Interpolate double logarithmically to find sigma(srt)
!
  Do ie = 1, 15
    If (earray(ie)==plab) Then
      threpi = xarray(ie)
      Return
    Else If (earray(ie)>plab) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      threpi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
      Return
    End If
  End Do
  Return
End Function threpi
!*****************************************
!*****************************************
! for ppi(+)-->DELTA+omega
!      real*4 function FOURPI(srt)
Real Function fourpi(srt)
!  This function contains the experimental pi+p-->DELTA+PION cross sections   *
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in mb                                    *
!  earray = EXPerimental table with proton energies in MeV                    *
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!                                                                             *
!*****************************************
!      real*4   xarray(10), earray(10)
  Real xarray(10), earray(10)
  Save
  Data xarray/0.0001, 1.986597, 6.411932, 7.636956, 9.598362, 9.889740, 10.24317, 10.80138, 11.86988, 12.83925/
  Data earray/2.468, 2.718, 2.968, 0.322E+01, 0.347E+01, 0.372E+01, 0.397E+01, 0.422E+01, 0.447E+01, 0.472E+01/

  pmass = 0.14
  pmass1 = 0.938
  fourpi = 0.000001
  If (srt<=1.52) Return
! 1.Calculate p(lab)  from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  plab = srt
  If (plab<earray(1)) Then
    fourpi = 0.00001
    Return
  End If
!
! 2.Interpolate double logarithmically to find sigma(srt)
!
  Do ie = 1, 10
    If (earray(ie)==plab) Then
      fourpi = xarray(ie)
      Return
    Else If (earray(ie)>plab) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      fourpi = exp(ymin+(alog(plab)-xmin)*(ymax-ymin)/(xmax-xmin))
      Return
    End If
  End Do
  Return
End Function fourpi
!*****************************************
!*****************************************
! for pion (rho or omega)+baryon resonance collisions
!      real*4 function reab(i1,i2,srt,ictrl)
Real Function reab(i1, i2, srt, ictrl)
!  This function calculates the cross section for
!  pi+Delta(N*)-->N+PION process                                              *
!  srt    = DSQRT(s) in GeV                                                   *
!  reab   = cross section in fm**2                                            *
!  ictrl=1,2,3 for pion, rho and omega+D(N*)
!***************************************
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (amn=0.938, ap1=0.14, arho=0.77, aomega=0.782)
  Parameter (maxx=20, maxz=24)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
!c      SAVE /DD/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Save
  lb1 = lb(i1)
  lb2 = lb(i2)
  reab = 0
  If (ictrl==1 .And. srt<=(amn+2.*ap1+0.02)) Return
  If (ictrl==3 .And. srt<=(amn+ap1+aomega+0.02)) Return
  pin2 = ((srt**2+ap1**2-amn**2)/(2.*srt))**2 - ap1**2
  If (pin2<=0) Return
! for pion+D(N*)-->pion+N
  If (ictrl==1) Then
    If (e(i1)>1) Then
      ed = e(i1)
    Else
      ed = e(i2)
    End If
    pout2 = ((srt**2+ap1**2-ed**2)/(2.*srt))**2 - ap1**2
    If (pout2<=0) Return
    xpro = twopi(srt)/10.
    factor = 1/3.
    If (((lb1==8 .And. lb2==5) .Or. (lb1==5 .And. lb2==8)) .Or. ((lb1==-8 .And. lb2==3) .Or. (lb1==3 .And. lb2==-8))) factor = 1/4.
    If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) factor = 1.
    reab = factor*pin2/pout2*xpro
    Return
  End If
! for rho reabsorption
  If (ictrl==2) Then
    If (lb(i2)>=25) Then
      ed = e(i1)
      arho1 = e(i2)
    Else
      ed = e(i2)
      arho1 = e(i1)
    End If
    If (srt<=(amn+ap1+arho1+0.02)) Return
    pout2 = ((srt**2+arho1**2-ed**2)/(2.*srt))**2 - arho1**2
    If (pout2<=0) Return
    xpro = threpi(srt)/10.
    factor = 1/3.
    If (((lb1==8 .And. lb2==27) .Or. (lb1==27 .And. lb2==8)) .Or. ((lb1==-8 .And. lb2==25) .Or. (lb1==25 .And. lb2==-8))) factor = 1/4.
    If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) factor = 1.
    reab = factor*pin2/pout2*xpro
    Return
  End If
! for omega reabsorption
  If (ictrl==3) Then
    If (e(i1)>1) ed = e(i1)
    If (e(i2)>1) ed = e(i2)
    pout2 = ((srt**2+aomega**2-ed**2)/(2.*srt))**2 - aomega**2
    If (pout2<=0) Return
    xpro = fourpi(srt)/10.
    factor = 1/6.
    If ((iabs(lb1)>=10 .And. iabs(lb1)<=13) .Or. (iabs(lb2)>=10 .And. iabs(lb2)<=13)) factor = 1./3.
    reab = factor*pin2/pout2*xpro
  End If
  Return
End Function reab
!*****************************************
! for the reabsorption of two resonances
! This function calculates the cross section for
! DD-->NN, N*N*-->NN and DN*-->NN
!      real*4 function reab2d(i1,i2,srt)
Real Function reab2d(i1, i2, srt)
!  srt    = DSQRT(s) in GeV                                                   *
!  reab   = cross section in mb
!***************************************
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (amn=0.938, ap1=0.14, arho=0.77, aomega=0.782)
  Parameter (maxx=20, maxz=24)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
!c      SAVE /DD/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Save
  reab2d = 0
  lb1 = iabs(lb(i1))
  lb2 = iabs(lb(i2))
  ed1 = e(i1)
  ed2 = e(i2)
  pin2 = (srt/2.)**2 - amn**2
  pout2 = ((srt**2+ed1**2-ed2**2)/(2.*srt))**2 - ed1**2
  If (pout2<=0) Return
  xpro = x2pi(srt)
  factor = 1/4.
  If ((lb1>=10 .And. lb1<=13) .And. (lb2>=10 .And. lb2<=13)) factor = 1.
  If ((lb1>=6 .And. lb1<=9) .And. (lb2>10 .And. lb2<=13)) factor = 1/2.
  If ((lb2>=6 .And. lb2<=9) .And. (lb1>10 .And. lb1<=13)) factor = 1/2.
  reab2d = factor*pin2/pout2*xpro
  Return
End Function reab2d
!**************************************
Subroutine rotate(px0, py0, pz0, px, py, pz)
  Save
! purpose: rotate the momentum of a particle in the CMS of p1+p2 such that
! the x' y' and z' in the cms of p1+p2 is the same as the fixed x y and z
! quantities:
!            px0,py0 and pz0 are the cms momentum of the incoming colliding
!            particles
!            px, py and pz are the cms momentum of any one of the particles
!            after the collision to be rotated
!**************************************
! the momentum, polar and azimuthal angles of the incoming momentm
  pr0 = sqrt(px0**2+py0**2+pz0**2)
  If (pr0==0) pr0 = 0.00000001
  c2 = pz0/pr0
  If (px0==0.0 .And. py0==0.0) Then
    t2 = 0.0
  Else
    t2 = atan2(py0, px0)
  End If

!lin-9/2012: check argument in sqrt():
  scheck = 1.0 - c2**2
  If (scheck<0) Then
    Write (99, *) 'scheck45: ', scheck
    scheck = 0.
  End If
  s2 = sqrt(scheck)
!      S2  =  SQRT( 1.0 - C2**2 )

  ct2 = cos(t2)
  st2 = sin(t2)
! the momentum, polar and azimuthal angles of the momentum to be rotated
  pr = sqrt(px**2+py**2+pz**2)
  If (pr==0) pr = 0.0000001
  c1 = pz/pr
  If (px==0 .And. py==0) Then
    t1 = 0.
  Else
    t1 = atan2(py, px)
  End If

!lin-9/2012: check argument in sqrt():
  scheck = 1.0 - c1**2
  If (scheck<0) Then
    Write (99, *) 'scheck46: ', scheck
    scheck = 0.
  End If
  s1 = sqrt(scheck)
!      S1   = SQRT( 1.0 - C1**2 )

  ct1 = cos(t1)
  st1 = sin(t1)
  ss = c2*s1*ct1 + s2*c1
! THE MOMENTUM AFTER ROTATION
  px = pr*(ss*ct2-s1*st1*st2)
  py = pr*(ss*st2+s1*st1*ct2)
  pz = pr*(c1*c2-s1*s2*ct1)
  Return
End Subroutine rotate
!*****************************************
!      real*4 function Xpp(srt)
Real Function xpp(srt)
!  This function contains the experimental total n-p cross sections           *
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in mb                                    *
!  earray = EXPerimental table with proton energies in MeV                    *
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!  WITH A CUTOFF AT 55MB                                                      *
!*****************************************
!      real*4   xarray(14), earray(14)
  Real xarray(14), earray(14)
  Save
  Data earray/20., 30., 40., 60., 80., 100., 170., 250., 310., 350., 460., 560., 660., 800./
  Data xarray/150., 90., 80.6, 48.0, 36.6, 31.6, 25.9, 24.0, 23.1, 24.0, 28.3, 33.6, 41.5, 47/

  pmass = 0.9383
! 1.Calculate E_kin(lab) [MeV] from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  ekin = 2000.*pmass*((srt/(2.*pmass))**2-1.)
  If (ekin<earray(1)) Then
    xpp = xarray(1)
    If (xpp>55) xpp = 55
    Return
  End If
  If (ekin>earray(14)) Then
    xpp = xarray(14)
    Return
  End If
!
!
! 2.Interpolate double logarithmically to find sigma(srt)
!
  Do ie = 1, 14
    If (earray(ie)==ekin) Then
      xpp = xarray(ie)
      If (xpp>55) xpp = 55.
      Return
    End If
    If (earray(ie)>ekin) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      xpp = exp(ymin+(alog(ekin)-xmin)*(ymax-ymin)/(xmax-xmin))
      If (xpp>55) xpp = 55.
      Goto 50
    End If
  End Do
  50 Continue
  Return
End Function xpp
!*****************************************
Real Function xnp(srt)
!  This function contains the experimental total n-p cross sections           *
!  srt    = DSQRT(s) in GeV                                                   *
!  xsec   = production cross section in mb                                    *
!  earray = EXPerimental table with proton energies in MeV                    *
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!  WITH  A CUTOFF AT 55MB                                                *
!*****************************************
!      real*4   xarray(11), earray(11)
  Real xarray(11), earray(11)
  Save
  Data earray/20., 30., 40., 60., 90., 135.0, 200., 300., 400., 600., 800./
  Data xarray/410., 270., 214.5, 130., 78., 53.5, 41.6, 35.9, 34.2, 34.3, 34.9/

  pmass = 0.9383
! 1.Calculate E_kin(lab) [MeV] from srt [GeV]
!   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
  ekin = 2000.*pmass*((srt/(2.*pmass))**2-1.)
  If (ekin<earray(1)) Then
    xnp = xarray(1)
    If (xnp>55) xnp = 55
    Return
  End If
  If (ekin>earray(11)) Then
    xnp = xarray(11)
    Return
  End If
!
!Interpolate double logarithmically to find sigma(srt)
!
  Do ie = 1, 11
    If (earray(ie)==ekin) Then
      xnp = xarray(ie)
      If (xnp>55) xnp = 55.
      Return
    End If
    If (earray(ie)>ekin) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      xnp = exp(ymin+(alog(ekin)-xmin)*(ymax-ymin)/(xmax-xmin))
      If (xnp>55) xnp = 55
      Goto 50
    End If
  End Do
  50 Continue
  Return
End Function xnp
!******************************
Function ptr(ptmax, iseed)
! (2) Generate the transverse momentum
!     OF nucleons
!******************************
  Common /table/xarray(0:1000), earray(0:1000)
!c      SAVE /TABLE/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
  ptr = 0.
  If (ptmax<=1.E-02) Then
    ptr = ptmax
    Return
  End If
  If (ptmax>2.01) ptmax = 2.01
  tryial = ptdis(ptmax)/ptdis(2.01)
  xt = ranart(nseed)*tryial
! look up the table and
!Interpolate double logarithmically to find pt
  Do ie = 1, 200
    If (earray(ie)==xt) Then
      ptr = xarray(ie)
      Return
    End If
    If (xarray(ie-1)<=0.00001) Goto 50
    If (xarray(ie)<=0.00001) Goto 50
    If (earray(ie-1)<=0.00001) Goto 50
    If (earray(ie)<=0.00001) Goto 50
    If (earray(ie)>xt) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      ptr = exp(ymin+(alog(xt)-xmin)*(ymax-ymin)/(xmax-xmin))
      If (ptr>ptmax) ptr = ptmax
      Return
    End If
  50 End Do
  Return
End Function ptr

!*********************************
!*********************************
!                                                                      *
!                                                                      *
Subroutine xnd(px, py, pz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
!     PURPOSE:                                                         *
!             calculate NUCLEON-BARYON RESONANCE inelatic Xsection     *
!     NOTE   :                                                         *
!     QUANTITIES:                                                 *
!                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    *
!                      N12,                                            *
!                      M12=1 FOR p+n-->delta(+)+ n                     *
!                          2     p+n-->delta(0)+ p                     *
!                          3     p+p-->delta(++)+n                     *
!                          4     p+p-->delta(+)+p                      *
!                          5     n+n-->delta(0)+n                      *
!                          6     n+n-->delta(-)+p                      *
!                          7     n+p-->N*(0)(1440)+p                   *
!                          8     n+p-->N*(+)(1440)+n                   *
!                        9     p+p-->N*(+)(1535)+p                     *
!                        10    n+n-->N*(0)(1535)+n                     *
!                         11    n+p-->N*(+)(1535)+n                     *
!                        12    n+p-->N*(0)(1535)+p
!                        13    D(++)+D(-)-->N*(+)(1440)+n
!                         14    D(++)+D(-)-->N*(0)(1440)+p
!                        15    D(+)+D(0)--->N*(+)(1440)+n
!                        16    D(+)+D(0)--->N*(0)(1440)+p
!                        17    D(++)+D(0)-->N*(+)(1535)+p
!                        18    D(++)+D(-)-->N*(0)(1535)+p
!                        19    D(++)+D(-)-->N*(+)(1535)+n
!                        20    D(+)+D(+)-->N*(+)(1535)+p
!                        21    D(+)+D(0)-->N*(+)(1535)+n
!                        22    D(+)+D(0)-->N*(0)(1535)+p
!                        23    D(+)+D(-)-->N*(0)(1535)+n
!                        24    D(0)+D(0)-->N*(0)(1535)+n
!                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p
!                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n
!                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n
!                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p
!                        29    N*(+)(14)+D+-->N*(+)(15)+p
!                        30    N*(+)(14)+D0-->N*(+)(15)+n
!                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n
!                        32    N*(0)(14)+D++--->N*(+)(15)+p
!                        33    N*(0)(14)+D+--->N*(+)(15)+n
!                        34    N*(0)(14)+D+--->N*(0)(15)+p
!                        35    N*(0)(14)+D0-->N*(0)(15)+n
!                        36    N*(+)(14)+D0--->N*(0)(15)+p
!                            and more
!**********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
!c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
!c      SAVE /gg/
  Common /input/nstar, ndirct, dir
!c      SAVE /INPUT/
  Common /nn/nnn
!c      SAVE /NN/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
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
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Save

!-----------------------------------------------------------------------
  xinel = 0.
  sigk = 0
  xsk1 = 0
  xsk2 = 0
  xsk3 = 0
  xsk4 = 0
  xsk5 = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
!     CAN HAPPEN ANY MORE ==> RETURN (2.04 = 2*AVMASS + PI-MASS+0.02)
  If (srt<2.04) Return
! Resonance absorption or Delta + N-->N*(1440), N*(1535)
! COM: TEST FOR DELTA OR N* ABSORPTION
!      IN THE PROCESS DELTA+N-->NN, N*+N-->NN
  prf = sqrt(0.25*srt**2-avmass**2)
  If (em1>1.) Then
    deltam = em1
  Else
    deltam = em2
  End If
  renom = deltam*prf**2/denom(srt, 1.)/pr
  renomn = deltam*prf**2/denom(srt, 2.)/pr
  renom1 = deltam*prf**2/denom(srt, -1.)/pr
! avoid the inelastic collisions between n+delta- -->N+N
!       and p+delta++ -->N+N due to charge conservation,
!       but they can scatter to produce kaons
  If ((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) renom = 0.
  If ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) renom = 0.
  If ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) renom = 0.
  If ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9)) renom = 0.
  Call m1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
  x1440 = (3./4.)*sigma(srt, 2, 0, 1)
! CROSS SECTION FOR KAON PRODUCTION from the four channels
! for NLK channel
  akp = 0.498
  ak0 = 0.498
  ana = 0.94
  ada = 1.232
  al = 1.1157
  as = 1.1197
  xsk1 = 0
  xsk2 = 0
  xsk3 = 0
  xsk4 = 0
!      !! phi production
  xsk5 = 0
  t1nlk = ana + al + akp
  If (srt<=t1nlk) Goto 222
  xsk1 = 1.5*pplpk(srt)
! for DLK channel
  t1dlk = ada + al + akp
  t2dlk = ada + al - akp
  If (srt<=t1dlk) Goto 222
  es = srt
  pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
  pmdlk = sqrt(pmdlk2)
  xsk3 = 1.5*pplpk(srt)
! for NSK channel
  t1nsk = ana + as + akp
  t2nsk = ana + as - akp
  If (srt<=t1nsk) Goto 222
  pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
  pmnsk = sqrt(pmnsk2)
  xsk2 = 1.5*(ppk1(srt)+ppk0(srt))
! for DSK channel
  t1dsk = ada + as + akp
  t2dsk = ada + as - akp
  If (srt<=t1dsk) Goto 222
  pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
  pmdsk = sqrt(pmdsk2)
  xsk4 = 1.5*(ppk1(srt)+ppk0(srt))
!sp11/21/01
! phi production
  If (srt<=(2.*amn+aphi)) Goto 222
!  !! mb put the correct form
  xsk5 = 0.0001
!sp11/21/01 end

! THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN
  222 sigk = xsk1 + xsk2 + xsk3 + xsk4

!bz3/7/99 neutralk
  xsk1 = 2.0*xsk1
  xsk2 = 2.0*xsk2
  xsk3 = 2.0*xsk3
  xsk4 = 2.0*xsk4
  sigk = 2.0*sigk + xsk5
!bz3/7/99 neutralk end

! avoid the inelastic collisions between n+delta- -->N+N
!       and p+delta++ -->N+N due to charge conservation,
!       but they can scatter to produce kaons
  If (((iabs(lb(i1))==2) .And. (iabs(lb(i2))==6)) .Or. ((iabs(lb(i2))==2) .And. (iabs(lb(i1))==6)) .Or. ((iabs(lb(i1))==1) .And. (iabs(lb(i2))==9)) .Or. ((iabs(lb(i2))==1) .And. (iabs(lb(i1))==9))) Then
    xinel = sigk
    Return
  End If
! WE DETERMINE THE REACTION CHANNELS IN THE FOLLOWING
! FOR n+delta(++)-->p+p or n+delta(++)-->n+N*(+)(1440),n+N*(+)(1535)
! REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN,
  If (lb(i1)*lb(i2)==18 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
    signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
! FOR p+delta(-)-->n+n or p+delta(-)-->n+N*(0)(1440),n+N*(0)(1535)
! REABSORPTION OR N*(1535) PRODUCTION LIKE IN P+P OR N*(1440) LIKE PN,
  If (lb(i1)*lb(i2)==6 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
    signd = sigma(srt, 1, 1, 0) + 0.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
! FOR p+delta(+)-->p+p, N*(+)(144)+p, N*(+)(1535)+p
!bz11/25/98
  If (lb(i1)*lb(i2)==8 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
    signd = 1.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
! FOR n+delta(0)-->n+n, N*(0)(144)+n, N*(0)(1535)+n
  If (lb(i1)*lb(i2)==14 .And. (iabs(lb(i1))==2 .And. iabs(lb(i2))==2)) Then
    signd = 1.5*sigma(srt, 1, 1, 1)
    sigdn = 0.25*signd*renom
    xinel = sigdn + x1440 + x1535 + sigk
    Return
  End If
! FOR n+delta(+)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
!                       N*(+)(1535)+n,N*(0)(1535)+p
  If (lb(i1)*lb(i2)==16 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
    signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
    sigdn = 0.5*signd*renom
    xinel = sigdn + 2.*x1440 + 2.*x1535 + sigk
    Return
  End If
! FOR p+delta(0)-->n+p, N*(+)(1440)+n,N*(0)(1440)+p,
!                       N*(+)(1535)+n,N*(0)(1535)+p
  If (lb(i1)*lb(i2)==7) Then
    signd = 0.5*sigma(srt, 1, 1, 1) + 0.25*sigma(srt, 1, 1, 0)
    sigdn = 0.5*signd*renom
    xinel = sigdn + 2.*x1440 + 2.*x1535 + sigk
    Return
  End If
! FOR p+N*(0)(14)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
! OR  P+N*(0)(14)-->D(+)+N, D(0)+P,
  If (lb(i1)*lb(i2)==10 .And. (iabs(lb(i1))==1 .Or. iabs(lb(i2))==1)) Then
    signd = (3./4.)*sigma(srt, 2, 0, 1)
    sigdn = signd*renomn
    xinel = sigdn + x1535 + sigk
    Return
  End If
! FOR n+N*(+)-->p+n, N*(+)(1535)+n,N*(0)(1535)+p
  If (lb(i1)*lb(i2)==22 .And. (iabs(lb(i1))==2 .Or. iabs(lb(i2))==2)) Then
    signd = (3./4.)*sigma(srt, 2, 0, 1)
    sigdn = signd*renomn
    xinel = sigdn + x1535 + sigk
    Return
  End If
! FOR N*(1535)+N-->N+N COLLISIONS
  If ((iabs(lb(i1))==12) .Or. (iabs(lb(i1))==13) .Or. (iabs(lb(i2))==12) .Or. (iabs(lb(i2))==13)) Then
    signd = x1535
    sigdn = signd*renom1
    xinel = sigdn + sigk
    Return
  End If
  Return
End Subroutine xnd
!*********************************
!                                                                      *
!                                                                      *
Subroutine xddin(px, py, pz, srt, i1, i2, xinel, sigk, xsk1, xsk2, xsk3, xsk4, xsk5)
!     PURPOSE:                                                         *
!             DEALING WITH BARYON RESONANCE-BARYON RESONANCE COLLISIONS*
!     NOTE   :                                                         *
!           VALID ONLY FOR BARYON-BARYON-DISTANCES LESS THAN 1.32 FM   *
!           (1.32 = 2 * HARD-CORE-RADIUS [HRC] )                       *
!     QUANTITIES:                                                 *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   *
!           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         *
!           IBLOCK   - THE INFORMATION BACK                            *
!                      0-> COLLISION CANNOT HAPPEN                     *
!                      1-> N-N ELASTIC COLLISION                       *
!                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          *
!                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          *
!                      4-> N+N->N+N+PION,DIRTCT PROCESS                *
!                     5-> DELTA(N*)+DELTA(N*)   TOTAL   COLLISIONS    *
!           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      *
!                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    *
!                      N12,                                            *
!                      M12=1 FOR p+n-->delta(+)+ n                     *
!                          2     p+n-->delta(0)+ p                     *
!                          3     p+p-->delta(++)+n                     *
!                          4     p+p-->delta(+)+p                      *
!                          5     n+n-->delta(0)+n                      *
!                          6     n+n-->delta(-)+p                      *
!                          7     n+p-->N*(0)(1440)+p                   *
!                          8     n+p-->N*(+)(1440)+n                   *
!                        9     p+p-->N*(+)(1535)+p                     *
!                        10    n+n-->N*(0)(1535)+n                     *
!                         11    n+p-->N*(+)(1535)+n                     *
!                        12    n+p-->N*(0)(1535)+p
!                        13    D(++)+D(-)-->N*(+)(1440)+n
!                         14    D(++)+D(-)-->N*(0)(1440)+p
!                        15    D(+)+D(0)--->N*(+)(1440)+n
!                        16    D(+)+D(0)--->N*(0)(1440)+p
!                        17    D(++)+D(0)-->N*(+)(1535)+p
!                        18    D(++)+D(-)-->N*(0)(1535)+p
!                        19    D(++)+D(-)-->N*(+)(1535)+n
!                        20    D(+)+D(+)-->N*(+)(1535)+p
!                        21    D(+)+D(0)-->N*(+)(1535)+n
!                        22    D(+)+D(0)-->N*(0)(1535)+p
!                        23    D(+)+D(-)-->N*(0)(1535)+n
!                        24    D(0)+D(0)-->N*(0)(1535)+n
!                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p
!                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n
!                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n
!                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p
!                        29    N*(+)(14)+D+-->N*(+)(15)+p
!                        30    N*(+)(14)+D0-->N*(+)(15)+n
!                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n
!                        32    N*(0)(14)+D++--->N*(+)(15)+p
!                        33    N*(0)(14)+D+--->N*(+)(15)+n
!                        34    N*(0)(14)+D+--->N*(0)(15)+p
!                        35    N*(0)(14)+D0-->N*(0)(15)+n
!                        36    N*(+)(14)+D0--->N*(0)(15)+p
!                        +++
!               AND MORE CHANNELS AS LISTED IN THE NOTE BOOK
!
! NOTE ABOUT N*(1440) RESORANCE:                                       *
!     As it has been discussed in VerWest's paper,I= 1 (initial isospin)
!     channel can all be attributed to delta resorance while I= 0      *
!     channel can all be  attribured to N* resorance.Only in n+p       *
!     one can have I=0 channel so is the N*(1440) resorance            *
! REFERENCES:    J. CUGNON ET AL., NUCL. PHYS. A352, 505 (1981)        *
!                    Y. KITAZOE ET AL., PHYS. LETT. 166B, 35 (1986)    *
!                    B. VerWest el al., PHYS. PRV. C25 (1982)1979      *
!                    Gy. Wolf  et al, Nucl Phys A517 (1990) 615        *
!                    CUTOFF = 2 * AVMASS + 20 MEV                      *
!                                                                      *
!       for N*(1535) we use the parameterization by Gy. Wolf et al     *
!       Nucl phys A552 (1993) 349, added May 18, 1994                  *
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, aka=0.498, aphi=1.020, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
!c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
!c      SAVE /gg/
  Common /input/nstar, ndirct, dir
!c      SAVE /INPUT/
  Common /nn/nnn
!c      SAVE /NN/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
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
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Save
!-----------------------------------------------------------------------
  xinel = 0
  sigk = 0
  xsk1 = 0
  xsk2 = 0
  xsk3 = 0
  xsk4 = 0
  xsk5 = 0
  em1 = e(i1)
  em2 = e(i2)
  pr = sqrt(px**2+py**2+pz**2)
!     IF THERE WERE 2 N*(1535) AND THEY DIDN'T SCATT. ELAST.,
!     ALLOW THEM TO PRODUCE KAONS. NO OTHER INELASTIC CHANNELS
!     ARE KNOWN
!       if((lb(i1).ge.12).and.(lb(i2).ge.12))return
!     ALL the inelastic collisions between N*(1535) and Delta as well
!     as N*(1440) TO PRODUCE KAONS, NO OTHER CHANNELS ARE KNOWN
!       if((lb(i1).ge.12).and.(lb(i2).ge.3))return
!       if((lb(i2).ge.12).and.(lb(i1).ge.3))return
!     calculate the N*(1535) production cross section in I1+I2 collisions
  Call n1535(iabs(lb(i1)), iabs(lb(i2)), srt, x1535)
!
! for Delta+Delta-->N*(1440 OR 1535)+N AND N*(1440)+N*(1440)-->N*(1535)+X
!     AND DELTA+N*(1440)-->N*(1535)+X
! WE ASSUME THEY HAVE THE SAME CROSS SECTIONS as CORRESPONDING N+N COLLISION):
! FOR D++D0, D+D+,D+D-,D0D0,N*+N*+,N*0N*0,N*(+)D+,N*(+)D(-),N*(0)D(0)
! N*(1535) production, kaon production and reabsorption through
! D(N*)+D(N*)-->NN are ALLOWED.
! CROSS SECTION FOR KAON PRODUCTION from the four channels are
! for NLK channel
  akp = 0.498
  ak0 = 0.498
  ana = 0.94
  ada = 1.232
  al = 1.1157
  as = 1.1197
  xsk1 = 0
  xsk2 = 0
  xsk3 = 0
  xsk4 = 0
  t1nlk = ana + al + akp
  If (srt<=t1nlk) Goto 222
  xsk1 = 1.5*pplpk(srt)
! for DLK channel
  t1dlk = ada + al + akp
  t2dlk = ada + al - akp
  If (srt<=t1dlk) Goto 222
  es = srt
  pmdlk2 = (es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
  pmdlk = sqrt(pmdlk2)
  xsk3 = 1.5*pplpk(srt)
! for NSK channel
  t1nsk = ana + as + akp
  t2nsk = ana + as - akp
  If (srt<=t1nsk) Goto 222
  pmnsk2 = (es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
  pmnsk = sqrt(pmnsk2)
  xsk2 = 1.5*(ppk1(srt)+ppk0(srt))
! for DSK channel
  t1dsk = ada + as + akp
  t2dsk = ada + as - akp
  If (srt<=t1dsk) Goto 222
  pmdsk2 = (es**2-t1dsk**2)*(es**2-t2dsk**2)/(4.*es**2)
  pmdsk = sqrt(pmdsk2)
  xsk4 = 1.5*(ppk1(srt)+ppk0(srt))
!sp11/21/01
! phi production
  If (srt<=(2.*amn+aphi)) Goto 222
!  !! mb put the correct form
  xsk5 = 0.0001
!sp11/21/01 end
! THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN
  222 sigk = xsk1 + xsk2 + xsk3 + xsk4

!bz3/7/99 neutralk
  xsk1 = 2.0*xsk1
  xsk2 = 2.0*xsk2
  xsk3 = 2.0*xsk3
  xsk4 = 2.0*xsk4
  sigk = 2.0*sigk + xsk5
!bz3/7/99 neutralk end

  idd = iabs(lb(i1)*lb(i2))
! The reabsorption cross section for the process
! D(N*)D(N*)-->NN is
  s2d = reab2d(i1, i2, srt)

!bz3/16/99 pion
  s2d = 0.
!bz3/16/99 pion end

!(1) N*(1535)+D(N*(1440)) reactions
!    we allow kaon production and reabsorption only
  If (((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=12)) .Or. ((iabs(lb(i1))>=12) .And. (iabs(lb(i2))>=6)) .Or. ((iabs(lb(i2))>=12) .And. (iabs(lb(i1))>=6))) Then
    xinel = sigk + s2d
    Return
  End If
! channels have the same charge as pp
  If ((idd==63) .Or. (idd==64) .Or. (idd==48) .Or. (idd==49) .Or. (idd==11*11) .Or. (idd==10*10) .Or. (idd==88) .Or. (idd==66) .Or. (idd==90) .Or. (idd==70)) Then
    xinel = x1535 + sigk + s2d
    Return
  End If
! IN DELTA+N*(1440) and N*(1440)+N*(1440) COLLISIONS,
! N*(1535), kaon production and reabsorption are ALLOWED
! IN N*(1440)+N*(1440) COLLISIONS, ONLY N*(1535) IS ALLOWED
  If ((idd==110) .Or. (idd==77) .Or. (idd==80)) Then
    xinel = x1535 + sigk + s2d
    Return
  End If
  If ((idd==54) .Or. (idd==56)) Then
! LIKE FOR N+P COLLISION,
! IN DELTA+DELTA COLLISIONS BOTH N*(1440) AND N*(1535) CAN BE PRODUCED
    sig2 = (3./4.)*sigma(srt, 2, 0, 1)
    xinel = 2.*(sig2+x1535) + sigk + s2d
    Return
  End If
  Return
End Subroutine xddin
!*****************************************


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
Real Function dirct1(srt)
!  This function contains the experimental, direct pion(+) + p cross sections *
!  srt    = DSQRT(s) in GeV                                                   *
!  dirct1  = cross section in fm**2                                     *
!  earray = EXPerimental table with the srt
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!*****************************************
!      real*4   xarray(122), earray(122)
  Real xarray(122), earray(122)
  Save
  Data earray/1.568300, 1.578300, 1.588300, 1.598300, 1.608300, 1.618300, 1.628300, 1.638300, 1.648300, 1.658300, 1.668300, 1.678300, 1.688300, 1.698300, 1.708300, 1.718300, 1.728300, 1.738300, 1.748300, 1.758300, 1.768300, 1.778300, 1.788300, 1.798300, 1.808300, 1.818300, 1.828300, 1.838300, 1.848300, 1.858300, 1.868300, 1.878300, 1.888300, 1.898300, 1.908300, 1.918300, 1.928300, 1.938300, 1.948300, 1.958300, 1.968300, 1.978300, 1.988300, 1.998300, 2.008300, 2.018300, 2.028300, 2.038300, 2.048300, 2.058300, 2.068300, 2.078300, 2.088300, 2.098300, 2.108300, 2.118300, 2.128300, 2.138300, 2.148300, 2.158300, 2.168300, 2.178300, 2.188300, 2.198300, 2.208300, 2.218300, 2.228300, 2.238300, 2.248300, 2.258300, 2.268300, 2.278300, 2.288300, 2.298300, 2.308300, 2.318300, 2.328300, 2.338300, 2.348300, 2.358300, 2.368300, 2.378300, 2.388300, 2.398300, 2.408300, 2.418300, 2.428300, 2.438300, 2.448300, 2.458300, 2.468300, 2.478300, 2.488300, 2.498300, 2.508300, 2.518300, 2.528300, 2.538300, 2.548300, 2.558300, 2.568300, &
    2.578300, 2.588300, 2.598300, 2.608300, 2.618300, 2.628300, 2.638300, 2.648300, 2.658300, 2.668300, 2.678300, 2.688300, 2.698300, 2.708300, 2.718300, 2.728300, 2.738300, 2.748300, 2.758300, 2.768300, 2.778300/
  Data xarray/1.7764091E-02, 0.5643668, 0.8150568, 1.045565, 2.133695, 3.327922, 4.206488, 3.471242, 4.486876, 5.542213, 6.800052, 7.192446, 6.829848, 6.580306, 6.868410, 8.527946, 10.15720, 9.716511, 9.298335, 8.901310, 10.31213, 10.52185, 11.17630, 11.61639, 12.05577, 12.71596, 13.46036, 14.22060, 14.65449, 14.94775, 14.93310, 15.32907, 16.56481, 16.29422, 15.18548, 14.12658, 13.72544, 13.24488, 13.31003, 14.42680, 12.84423, 12.49025, 12.14858, 11.81870, 11.18993, 11.35816, 11.09447, 10.83873, 10.61592, 10.53754, 9.425521, 8.195912, 9.661075, 9.696192, 9.200142, 8.953734, 8.715461, 8.484999, 8.320765, 8.255512, 8.190969, 8.127125, 8.079508, 8.073004, 8.010611, 7.948909, 7.887895, 7.761005, 7.626290, 7.494696, 7.366132, 7.530178, 8.392097, 9.046881, 8.962544, 8.879403, 8.797427, 8.716601, 8.636904, 8.558312, 8.404368, 8.328978, 8.254617, 8.181265, 8.108907, 8.037527, 7.967100, 7.897617, 7.829057, 7.761405, 7.694647, 7.628764, 7.563742, 7.499570, 7.387562, 7.273281, 7.161334, 6.973375, 6.529592, 6.280323, &
    6.293136, 6.305725, 6.318097, 6.330258, 6.342214, 6.353968, 6.365528, 6.376895, 6.388079, 6.399081, 6.409906, 6.420560, 6.431045, 6.441367, 6.451529, 6.461533, 6.471386, 6.481091, 6.490650, 6.476413, 6.297259, 6.097826/

  If (srt<earray(1)) Then
    dirct1 = 0.00001
    Return
  End If
  If (srt>earray(122)) Then
    dirct1 = xarray(122)
    dirct1 = dirct1/10.
    Return
  End If
!
!Interpolate double logarithmically to find xdirct2(srt)
!
  Do ie = 1, 122
    If (earray(ie)==srt) Then
      dirct1 = xarray(ie)
      dirct1 = dirct1/10.
      Return
    End If
    If (earray(ie)>srt) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      dirct1 = exp(ymin+(alog(srt)-xmin)*(ymax-ymin)/(xmax-xmin))
      dirct1 = dirct1/10.
      Goto 50
    End If
  End Do
  50 Continue
  Return
End Function dirct1
!******************************
!*****************************************
Real Function dirct2(srt)
!  This function contains the experimental, direct pion(-) + p cross sections *
!  srt    = DSQRT(s) in GeV                                                   *
!  dirct2 = cross section in fm**2
!  earray = EXPerimental table with the srt
!  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
!*****************************************
!      real*4   xarray(122), earray(122)
  Real xarray(122), earray(122)
  Save
  Data earray/1.568300, 1.578300, 1.588300, 1.598300, 1.608300, 1.618300, 1.628300, 1.638300, 1.648300, 1.658300, 1.668300, 1.678300, 1.688300, 1.698300, 1.708300, 1.718300, 1.728300, 1.738300, 1.748300, 1.758300, 1.768300, 1.778300, 1.788300, 1.798300, 1.808300, 1.818300, 1.828300, 1.838300, 1.848300, 1.858300, 1.868300, 1.878300, 1.888300, 1.898300, 1.908300, 1.918300, 1.928300, 1.938300, 1.948300, 1.958300, 1.968300, 1.978300, 1.988300, 1.998300, 2.008300, 2.018300, 2.028300, 2.038300, 2.048300, 2.058300, 2.068300, 2.078300, 2.088300, 2.098300, 2.108300, 2.118300, 2.128300, 2.138300, 2.148300, 2.158300, 2.168300, 2.178300, 2.188300, 2.198300, 2.208300, 2.218300, 2.228300, 2.238300, 2.248300, 2.258300, 2.268300, 2.278300, 2.288300, 2.298300, 2.308300, 2.318300, 2.328300, 2.338300, 2.348300, 2.358300, 2.368300, 2.378300, 2.388300, 2.398300, 2.408300, 2.418300, 2.428300, 2.438300, 2.448300, 2.458300, 2.468300, 2.478300, 2.488300, 2.498300, 2.508300, 2.518300, 2.528300, 2.538300, 2.548300, 2.558300, 2.568300, &
    2.578300, 2.588300, 2.598300, 2.608300, 2.618300, 2.628300, 2.638300, 2.648300, 2.658300, 2.668300, 2.678300, 2.688300, 2.698300, 2.708300, 2.718300, 2.728300, 2.738300, 2.748300, 2.758300, 2.768300, 2.778300/
  Data xarray/0.5773182, 1.404156, 2.578629, 3.832013, 4.906011, 9.076963, 13.10492, 10.65975, 15.31156, 19.77611, 19.92874, 18.68979, 19.80114, 18.39536, 14.34269, 13.35353, 13.58822, 14.57031, 10.24686, 11.23386, 9.764803, 10.35652, 10.53539, 10.07524, 9.582198, 9.596469, 9.818489, 9.012848, 9.378012, 9.529244, 9.529698, 8.835624, 6.671396, 8.797758, 8.133437, 7.866227, 7.823946, 7.808504, 7.791755, 7.502062, 7.417275, 7.592349, 7.752028, 7.910585, 8.068122, 8.224736, 8.075289, 7.895902, 7.721359, 7.551512, 7.386224, 7.225343, 7.068739, 6.916284, 6.767842, 6.623294, 6.482520, 6.345404, 6.211833, 7.339510, 7.531462, 7.724824, 7.919620, 7.848021, 7.639856, 7.571083, 7.508881, 7.447474, 7.386855, 7.327011, 7.164454, 7.001266, 6.842526, 6.688094, 6.537823, 6.391583, 6.249249, 6.110689, 5.975790, 5.894200, 5.959503, 6.024602, 6.089505, 6.154224, 6.218760, 6.283128, 6.347331, 6.297411, 6.120248, 5.948606, 6.494864, 6.357106, 6.222824, 6.091910, 5.964267, 5.839795, 5.718402, 5.599994, 5.499146, 5.451325, 5.404156, &
    5.357625, 5.311721, 5.266435, 5.301964, 5.343963, 5.385833, 5.427577, 5.469200, 5.510702, 5.552088, 5.593359, 5.634520, 5.675570, 5.716515, 5.757356, 5.798093, 5.838732, 5.879272, 5.919717, 5.960068, 5.980941/

  If (srt<earray(1)) Then
    dirct2 = 0.00001
    Return
  End If
  If (srt>earray(122)) Then
    dirct2 = xarray(122)
    dirct2 = dirct2/10.
    Return
  End If
!
!Interpolate double logarithmically to find xdirct2(srt)
!
  Do ie = 1, 122
    If (earray(ie)==srt) Then
      dirct2 = xarray(ie)
      dirct2 = dirct2/10.
      Return
    End If
    If (earray(ie)>srt) Then
      ymin = alog(xarray(ie-1))
      ymax = alog(xarray(ie))
      xmin = alog(earray(ie-1))
      xmax = alog(earray(ie))
      dirct2 = exp(ymin+(alog(srt)-xmin)*(ymax-ymin)/(xmax-xmin))
      dirct2 = dirct2/10.
      Goto 50
    End If
  End Do
  50 Continue
  Return
End Function dirct2
!******************************
!*****************************
! this program calculates the elastic cross section for rho+nucleon
! through higher resonances
!       real*4 function ErhoN(em1,em2,lb1,lb2,srt)
Real Function erhon(em1, em2, lb1, lb2, srt)
! date : Dec. 19, 1994
! ****************************
!       implicit real*4 (a-h,o-z)
  Dimension arrayj(19), arrayl(19), arraym(19), arrayw(19), arrayb(19)
  Save
  Data arrayj/0.5, 1.5, 0.5, 0.5, 2.5, 2.5, 1.5, 0.5, 1.5, 3.5, 1.5, 0.5, 1.5, 0.5, 2.5, 0.5, 1.5, 2.5, 3.5/
  Data arrayl/1, 2, 0, 0, 2, 3, 2, 1, 1, 3, 1, 0, 2, 0, 3, 1, 1, 2, 3/
  Data arraym/1.44, 1.52, 1.535, 1.65, 1.675, 1.68, 1.70, 1.71, 1.72, 1.99, 1.60, 1.62, 1.70, 1.90, 1.905, 1.910, 1.86, 1.93, 1.95/
  Data arrayw/0.2, 0.125, 0.15, 0.15, 0.155, 0.125, 0.1, 0.11, 0.2, 0.29, 0.25, 0.16, 0.28, 0.15, 0.3, 0.22, 0.25, 0.25, 0.24/
  Data arrayb/0.15, 0.20, 0.05, 0.175, 0.025, 0.125, 0.1, 0.20, 0.53, 0.34, 0.05, 0.07, 0.15, 0.45, 0.45, 0.058, 0.08, 0.12, 0.08/

! the minimum energy for pion+delta collision
  pi = 3.1415926
  xs = 0
! include contribution from each resonance
  Do ir = 1, 19
!bz11/25/98
    If (ir<=8) Then
!       if(lb1*lb2.eq.27.OR.LB1*LB2.EQ.25*2)branch=0.
!       if(lb1*lb2.eq.26.OR.LB1*LB2.EQ.26*2)branch=1./3.
!       if(lb1*lb2.eq.27*2.OR.LB1*LB2.EQ.25)branch=2./3.
!       ELSE
!       if(lb1*lb2.eq.27.OR.LB1*LB2.EQ.25*2)branch=1.
!       if(lb1*lb2.eq.26.OR.LB1*LB2.EQ.26*2)branch=2./3.
!       if(lb1*lb2.eq.27*2.OR.LB1*LB2.EQ.25)branch=1./3.
!       ENDIF
      If (((lb1*lb2==27 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==25*2 .And. (lb1==2 .Or. lb2==2))) .Or. ((lb1*lb2==-25 .And. (lb1==-1 .Or. lb2==-1)) .Or. (lb1*lb2==-27*2 .And. (lb1==-2 .Or. lb2==-2)))) branch = 0.
      If ((iabs(lb1*lb2)==26 .And. (iabs(lb1)==1 .Or. iabs(lb2)==1)) .Or. (iabs(lb1*lb2)==26*2 .And. (iabs(lb1)==2 .Or. iabs(lb2)==2))) branch = 1./3.
      If (((lb1*lb2==27*2 .And. (lb1==2 .Or. lb2==2)) .Or. (lb1*lb2==25 .And. (lb1==1 .Or. lb2==1))) .Or. ((lb1*lb2==-25*2 .And. (lb1==-2 .Or. lb2==-2)) .Or. (lb1*lb2==-27 .And. (lb1==-1 .Or. lb2==-1)))) branch = 2./3.
    Else
      If (((lb1*lb2==27 .And. (lb1==1 .Or. lb2==1)) .Or. (lb1*lb2==25*2 .And. (lb1==2 .Or. lb2==2))) .Or. ((lb1*lb2==-25 .And. (lb1==-1 .Or. lb2==-1)) .Or. (lb1*lb2==-27*2 .And. (lb1==-2 .Or. lb2==-2)))) branch = 1.
      If ((iabs(lb1*lb2)==26 .And. (iabs(lb1)==1 .Or. iabs(lb2)==1)) .Or. (iabs(lb1*lb2)==26*2 .And. (iabs(lb1)==2 .Or. iabs(lb2)==2))) branch = 2./3.
      If (((lb1*lb2==27*2 .And. (lb1==2 .Or. lb2==2)) .Or. (lb1*lb2==25 .And. (lb1==1 .Or. lb2==1))) .Or. ((lb1*lb2==-25*2 .And. (lb1==-2 .Or. lb2==-2)) .Or. (lb1*lb2==-27 .And. (lb1==-1 .Or. lb2==-1)))) branch = 1./3.
    End If
!bz11/25/98end
    xs0 = fdr(arraym(ir), arrayj(ir), arrayl(ir), arrayw(ir), arrayb(ir), srt, em1, em2)
    xs = xs + 1.3*pi*branch*xs0*(0.1973)**2
  End Do
  erhon = xs
  Return
End Function erhon
!**************************8
!FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF
!KITAZOE'S FORMULA
!        REAL*4 FUNCTION FDR(DMASS,aj,al,width,widb0,srt,em1,em2)
Real Function fdr(dmass, aj, al, width, widb0, srt, em1, em2)
  Save
  amd = em1
  amp = em2
  ak02 = 0.25*(dmass**2-amd**2-amp**2)**2 - (amp*amd)**2
  If (ak02>0.) Then
    q0 = sqrt(ak02/dmass)
  Else
    q0 = 0.0
    fdr = 0
    Return
  End If
  ak2 = 0.25*(srt**2-amd**2-amp**2)**2 - (amp*amd)**2
  If (ak2>0.) Then
    q = sqrt(ak2/dmass)
  Else
    q = 0.00
    fdr = 0
    Return
  End If
  b = widb0*1.2*dmass/srt*(q/q0)**(2.*al+1)/(1.+0.2*(q/q0)**(2*al))
  fdr = (2.*aj+1)*width**2*b/((srt-dmass)**2+0.25*width**2)/(6.*q**2)
  Return
End Function fdr
!*****************************
! this program calculates the elastic cross section for pion+delta
! through higher resonances
!       REAL*4 FUNCTION DIRCT3(SRT)
Real Function dirct3(srt)
! date : Dec. 19, 1994
! ****************************
!     implicit real*4 (a-h,o-z)
  Dimension arrayj(17), arrayl(17), arraym(17), arrayw(17), arrayb(17)
  Save
  Data arrayj/1.5, 0.5, 2.5, 2.5, 1.5, 0.5, 1.5, 3.5, 1.5, 0.5, 1.5, 0.5, 2.5, 0.5, 1.5, 2.5, 3.5/
  Data arrayl/2, 0, 2, 3, 2, 1, 1, 3, 1, 0, 2, 0, 3, 1, 1, 2, 3/
  Data arraym/1.52, 1.65, 1.675, 1.68, 1.70, 1.71, 1.72, 1.99, 1.60, 1.62, 1.70, 1.90, 1.905, 1.910, 1.86, 1.93, 1.95/
  Data arrayw/0.125, 0.15, 0.155, 0.125, 0.1, 0.11, 0.2, 0.29, 0.25, 0.16, 0.28, 0.15, 0.3, 0.22, 0.25, 0.25, 0.24/
  Data arrayb/0.55, 0.6, 0.375, 0.6, 0.1, 0.15, 0.15, 0.05, 0.35, 0.3, 0.15, 0.1, 0.1, 0.22, 0.2, 0.09, 0.4/

! the minimum energy for pion+delta collision
  pi = 3.1415926
  amn = 0.938
  amp = 0.138
  xs = 0
! include contribution from each resonance
  branch = 1./3.
  Do ir = 1, 17
    If (ir>8) branch = 2./3.
    xs0 = fd1(arraym(ir), arrayj(ir), arrayl(ir), arrayw(ir), arrayb(ir), srt)
    xs = xs + 1.3*pi*branch*xs0*(0.1973)**2
  End Do
  dirct3 = xs
  Return
End Function dirct3
!**************************8
!FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF
!KITAZOE'S FORMULA
!        REAL*4 FUNCTION FD1(DMASS,aj,al,width,widb0,srt)
Real Function fd1(dmass, aj, al, width, widb0, srt)
  Save
  amn = 0.938
  amp = 0.138
  amd = amn
  ak02 = 0.25*(dmass**2-amd**2-amp**2)**2 - (amp*amd)**2
  If (ak02>0.) Then
    q0 = sqrt(ak02/dmass)
  Else
    q0 = 0.0
    fd1 = 0
    Return
  End If
  ak2 = 0.25*(srt**2-amd**2-amp**2)**2 - (amp*amd)**2
  If (ak2>0.) Then
    q = sqrt(ak2/dmass)
  Else
    q = 0.00
    fd1 = 0
    Return
  End If
  b = widb0*1.2*dmass/srt*(q/q0)**(2.*al+1)/(1.+0.2*(q/q0)**(2*al))
  fd1 = (2.*aj+1)*width**2*b/((srt-dmass)**2+0.25*width**2)/(2.*q**2)
  Return
End Function fd1
!*****************************
! this program calculates the elastic cross section for pion+delta
! through higher resonances
!       REAL*4 FUNCTION DPION(EM1,EM2,LB1,LB2,SRT)
Real Function dpion(em1, em2, lb1, lb2, srt)
! date : Dec. 19, 1994
! ****************************
!     implicit real*4 (a-h,o-z)
  Dimension arrayj(19), arrayl(19), arraym(19), arrayw(19), arrayb(19)
  Save
  Data arrayj/0.5, 1.5, 0.5, 0.5, 2.5, 2.5, 1.5, 0.5, 1.5, 3.5, 1.5, 0.5, 1.5, 0.5, 2.5, 0.5, 1.5, 2.5, 3.5/
  Data arrayl/1, 2, 0, 0, 2, 3, 2, 1, 1, 3, 1, 0, 2, 0, 3, 1, 1, 2, 3/
  Data arraym/1.44, 1.52, 1.535, 1.65, 1.675, 1.68, 1.70, 1.71, 1.72, 1.99, 1.60, 1.62, 1.70, 1.90, 1.905, 1.910, 1.86, 1.93, 1.95/
  Data arrayw/0.2, 0.125, 0.15, 0.15, 0.155, 0.125, 0.1, 0.11, 0.2, 0.29, 0.25, 0.16, 0.28, 0.15, 0.3, 0.22, 0.25, 0.25, 0.24/
  Data arrayb/0.15, 0.25, 0., 0.05, 0.575, 0.125, 0.379, 0.10, 0.10, 0.062, 0.45, 0.60, 0.6984, 0.05, 0.25, 0.089, 0.19, 0.2, 0.13/

! the minimum energy for pion+delta collision
  pi = 3.1415926
  amn = 0.94
  amp = 0.14
  xs = 0
! include contribution from each resonance
  Do ir = 1, 19
    branch = 0.
!bz11/25/98
    If (ir<=8) Then
!       IF(LB1*LB2.EQ.5*7.OR.LB1*LB2.EQ.3*8)branch=1./6.
!       IF(LB1*LB2.EQ.4*7.OR.LB1*LB2.EQ.4*8)branch=1./3.
!       IF(LB1*LB2.EQ.5*6.OR.LB1*LB2.EQ.3*9)branch=1./2.
!       ELSE
!       IF(LB1*LB2.EQ.5*8.OR.LB1*LB2.EQ.5*6)branch=2./5.
!       IF(LB1*LB2.EQ.3*9.OR.LB1*LB2.EQ.3*7)branch=2./5.
!       IF(LB1*LB2.EQ.5*7.OR.LB1*LB2.EQ.3*8)branch=8./15.
!       IF(LB1*LB2.EQ.4*7.OR.LB1*LB2.EQ.4*8)branch=1./15.
!       IF(LB1*LB2.EQ.4*9.OR.LB1*LB2.EQ.4*6)branch=3./5.
!       ENDIF
      If (((lb1*lb2==5*7 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==3*8 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-3*7 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-5*8 .And. (lb1==5 .Or. lb2==5)))) branch = 1./6.
      If ((iabs(lb1*lb2)==4*7 .And. (lb1==4 .Or. lb2==4)) .Or. (iabs(lb1*lb2)==4*8 .And. (lb1==4 .Or. lb2==4))) branch = 1./3.
      If (((lb1*lb2==5*6 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==3*9 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-3*6 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-5*9 .And. (lb1==5 .Or. lb2==5)))) branch = 1./2.
    Else
      If (((lb1*lb2==5*8 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==5*6 .And. (lb1==5 .Or. lb2==5))) .Or. ((lb1*lb2==-3*8 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-3*6 .And. (lb1==3 .Or. lb2==3)))) branch = 2./5.
      If (((lb1*lb2==3*9 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==3*7 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-5*9 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==-5*7 .And. (lb1==5 .Or. lb2==5)))) branch = 2./5.
      If (((lb1*lb2==5*7 .And. (lb1==5 .Or. lb2==5)) .Or. (lb1*lb2==3*8 .And. (lb1==3 .Or. lb2==3))) .Or. ((lb1*lb2==-3*7 .And. (lb1==3 .Or. lb2==3)) .Or. (lb1*lb2==-5*8 .And. (lb1==5 .Or. lb2==5)))) branch = 8./15.
      If ((iabs(lb1*lb2)==4*7 .And. (lb1==4 .Or. lb2==4)) .Or. (iabs(lb1*lb2)==4*8 .And. (lb1==4 .Or. lb2==4))) branch = 1./15.
      If ((iabs(lb1*lb2)==4*9 .And. (lb1==4 .Or. lb2==4)) .Or. (iabs(lb1*lb2)==4*6 .And. (lb1==4 .Or. lb2==4))) branch = 3./5.
    End If
!bz11/25/98end
    xs0 = fd2(arraym(ir), arrayj(ir), arrayl(ir), arrayw(ir), arrayb(ir), em1, em2, srt)
    xs = xs + 1.3*pi*branch*xs0*(0.1973)**2
  End Do
  dpion = xs
  Return
End Function dpion
!**************************8
!FUNCTION FDE(DMASS) GIVES DELTA MASS DISTRIBUTION BY USING OF
!KITAZOE'S FORMULA
!        REAL*4 FUNCTION FD2(DMASS,aj,al,width,widb0,EM1,EM2,srt)
Real Function fd2(dmass, aj, al, width, widb0, em1, em2, srt)
  Save
  amp = em1
  amd = em2
  ak02 = 0.25*(dmass**2-amd**2-amp**2)**2 - (amp*amd)**2
  If (ak02>0.) Then
    q0 = sqrt(ak02/dmass)
  Else
    q0 = 0.0
    fd2 = 0
    Return
  End If
  ak2 = 0.25*(srt**2-amd**2-amp**2)**2 - (amp*amd)**2
  If (ak2>0.) Then
    q = sqrt(ak2/dmass)
  Else
    q = 0.00
    fd2 = 0
    Return
  End If
  b = widb0*1.2*dmass/srt*(q/q0)**(2.*al+1)/(1.+0.2*(q/q0)**(2*al))
  fd2 = (2.*aj+1)*width**2*b/((srt-dmass)**2+0.25*width**2)/(4.*q**2)
  Return
End Function fd2
!**************************8
!   MASS GENERATOR for two resonances simultaneously
Subroutine rmasdd(srt, am10, am20, dmin1, dmin2, iseed, ic, dm1, dm2)
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
  amn = 0.94
  amp = 0.14
! the maximum mass for resonance 1
  dmax1 = srt - dmin2
! generate the mass for the first resonance
  5 ntry1 = 0
  ntry2 = 0
  ntry = 0
  ictrl = 0
  10 dm1 = ranart(nseed)*(dmax1-dmin1) + dmin1
  ntry1 = ntry1 + 1
! the maximum mass for resonance 2
  If (ictrl==0) dmax2 = srt - dm1
! generate the mass for the second resonance
  20 dm2 = ranart(nseed)*(dmax2-dmin2) + dmin2
  ntry2 = ntry2 + 1
! check the energy-momentum conservation with two masses
! q2 in the following is q**2*4*srt**2
  q2 = ((srt**2-dm1**2-dm2**2)**2-4.*dm1**2*dm2**2)
  If (q2<=0) Then
    dmax2 = dm2 - 0.01
!         dmax1=dm1-0.01
    ictrl = 1
    Goto 20
  End If
! determine the weight of the mass pair
  If (dmax1<am10) Then
    If (ic==1) fm1 = fmassd(dmax1)
    If (ic==2) fm1 = fmassn(dmax1)
    If (ic==3) fm1 = fmassd(dmax1)
    If (ic==4) fm1 = fmassd(dmax1)
  Else
    If (ic==1) fm1 = fmassd(am10)
    If (ic==2) fm1 = fmassn(am10)
    If (ic==3) fm1 = fmassd(am10)
    If (ic==4) fm1 = fmassd(am10)
  End If
  If (dmax2<am20) Then
    If (ic==1) fm2 = fmassd(dmax2)
    If (ic==2) fm2 = fmassn(dmax2)
    If (ic==3) fm2 = fmassn(dmax2)
    If (ic==4) fm2 = fmassr(dmax2)
  Else
    If (ic==1) fm2 = fmassd(am20)
    If (ic==2) fm2 = fmassn(am20)
    If (ic==3) fm2 = fmassn(am20)
    If (ic==4) fm2 = fmassr(am20)
  End If
  If (fm1==0.) fm1 = 1.E-04
  If (fm2==0.) fm2 = 1.E-04
  prob0 = fm1*fm2
  If (ic==1) prob = fmassd(dm1)*fmassd(dm2)
  If (ic==2) prob = fmassn(dm1)*fmassn(dm2)
  If (ic==3) prob = fmassd(dm1)*fmassn(dm2)
  If (ic==4) prob = fmassd(dm1)*fmassr(dm2)
  If (prob<=1.E-06) prob = 1.E-06
  fff = prob/prob0
  ntry = ntry + 1
  If (ranart(nseed)>fff .And. ntry<=20) Goto 10

!lin-2/26/03 limit the mass of (rho,Delta,N*1440) below a certain value
!     (here taken as its central value + 2* B-W fullwidth):
  If ((abs(am10-0.77)<=0.01 .And. dm1>1.07) .Or. (abs(am10-1.232)<=0.01 .And. dm1>1.47) .Or. (abs(am10-1.44)<=0.01 .And. dm1>2.14)) Goto 5
  If ((abs(am20-0.77)<=0.01 .And. dm2>1.07) .Or. (abs(am20-1.232)<=0.01 .And. dm2>1.47) .Or. (abs(am20-1.44)<=0.01 .And. dm2>2.14)) Goto 5

  Return
End Subroutine rmasdd
!FUNCTION Fmassd(DMASS) GIVES the delta MASS DISTRIBUTION
Real Function fmassd(dmass)
  Save
  am0 = 1.232
  fmassd = am0*width(dmass)/((dmass**2-am0**2)**2+am0**2*width(dmass)**2)
  Return
End Function fmassd
!FUNCTION Fmassn(DMASS) GIVES the N* MASS DISTRIBUTION
Real Function fmassn(dmass)
  Save
  am0 = 1.44
  fmassn = am0*w1440(dmass)/((dmass**2-am0**2)**2+am0**2*w1440(dmass)**2)
  Return
End Function fmassn
!FUNCTION Fmassr(DMASS) GIVES the rho MASS DISTRIBUTION
Real Function fmassr(dmass)
  Save
  am0 = 0.77
  wid = 0.153
  fmassr = am0*wid/((dmass**2-am0**2)**2+am0**2*wid**2)
  Return
End Function fmassr
!*********************************
! PURPOSE : flow analysis
! DATE : Feb. 1, 1995
!**********************************
Subroutine flow(nt)
!       IMPLICIT REAL*4 (A-H,O-Z)
  Parameter (pi=3.1415926, apion=0.13957, aka=0.498)
  Parameter (maxstr=150001, maxr=1, amu=0.9383, etam=0.5475)
  Dimension ypion(-80:80), ypr(-80:80), ykaon(-80:80)
  Dimension pxpion(-80:80), pxpro(-80:80), pxkaon(-80:80)
!----------------------------------------------------------------------*
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /rr/massr(0:maxr)
!c      SAVE /RR/
  Common /run/num
!c      SAVE /RUN/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Save
!----------------------------------------------------------------------*
  ycut1 = -2.6
  ycut2 = 2.6
  dy = 0.2
  ly = nint((ycut2-ycut1)/dy)
!**********************************
! initialize the transverse momentum counters
  Do kk = -80, 80
    pxpion(kk) = 0
    pxpro(kk) = 0
    pxkaon(kk) = 0
  End Do
  Do j = -ly, ly
    ypion(j) = 0
    ykaon(j) = 0
    ypr(j) = 0
  End Do
  nkaon = 0
  npr = 0
  npion = 0
  is = 0
  Do nrun = 1, num
    is = is + massr(nrun-1)
    Do j = 1, massr(nrun)
      i = j + is
! for protons go to 200 to calculate its rapidity and transvese momentum
! distributions
      e00 = sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+e(i)**2)
      y00 = 0.5*alog((e00+p(3,i))/(e00-p(3,i)))
      If (abs(y00)>=ycut2) Goto 20
      iy = nint(y00/dy)
      If (abs(iy)>=80) Goto 20
      If (e(i)==0) Goto 20
      If (lb(i)>=25) Goto 20
      If ((lb(i)<=5) .And. (lb(i)>=3)) Goto 50
      If (lb(i)==1 .Or. lb(i)==2) Goto 200
!bz3/10/99
!       if(lb(i).ge.6.and.lb(i).le.15)go to 200
      If (lb(i)>=6 .And. lb(i)<=17) Goto 200
!bz3/10/99 end
      If (lb(i)==23) Goto 400
      Goto 20
! calculate rapidity and transverse momentum distribution for pions
      50 npion = npion + 1
! (2) rapidity distribution in the cms frame
      ypion(iy) = ypion(iy) + 1
      pxpion(iy) = pxpion(iy) + p(1, i)/e(i)
      Goto 20
! calculate rapidity and transverse energy distribution for baryons
      200 npr = npr + 1
      pxpro(iy) = pxpro(iy) + p(1, i)/e(i)
      ypr(iy) = ypr(iy) + 1.
      Goto 20
      400 nkaon = nkaon + 1
      ykaon(iy) = ykaon(iy) + 1.
      pxkaon(iy) = pxkaon(iy) + p(1, i)/e(i)
    20 End Do
  End Do
! PRINT OUT NUCLEON'S TRANSVERSE MOMENTUM distribution
!       write(1041,*)Nt
!       write(1042,*)Nt
!       write(1043,*)Nt
!       write(1090,*)Nt
!       write(1091,*)Nt
!       write(1092,*)Nt
  Do npt = -10, 10
    If (ypr(npt)==0) Goto 101
    pxpro(npt) = -pxpro(npt)/ypr(npt)
    dnuc = pxpro(npt)/sqrt(ypr(npt))
!       WRITE(1041,*)NPT*DY,Pxpro(NPT),DNUC
! print pion's transverse momentum distribution
    101 If (ypion(npt)==0) Goto 102
    pxpion(npt) = -pxpion(npt)/ypion(npt)
    dnucp = pxpion(npt)/sqrt(ypion(npt))
!       WRITE(1042,*)NPT*DY,Pxpion(NPT),DNUCp
! kaons
    102 If (ykaon(npt)==0) Goto 3
    pxkaon(npt) = -pxkaon(npt)/ykaon(npt)
    dnuck = pxkaon(npt)/sqrt(ykaon(npt))
!       WRITE(1043,*)NPT*DY,Pxkaon(NPT),DNUCk
  3 End Do
!*******************************
! OUTPUT PION AND PROTON RAPIDITY DISTRIBUTIONS
  Do m = -ly, ly
! PROTONS
    dypr = 0
    If (ypr(m)/=0) dypr = sqrt(ypr(m))/float(nrun)/dy
    ypr(m) = ypr(m)/float(nrun)/dy
!       WRITE(1090,'(E11.3,2X,E11.3,2X,E11.3)')m*DY,YPR(M),DYPR
! PIONS
    dypion = 0
    If (ypion(m)/=0) dypion = sqrt(ypion(m))/float(nrun)/dy
    ypion(m) = ypion(m)/float(nrun)/dy
!       WRITE(1091,'(E11.3,2X,E11.3,2X,E11.3)')m*DY,YPION(M),DYPION
! KAONS
    dykaon = 0
    If (ykaon(m)/=0) dykaon = sqrt(ykaon(m))/float(nrun)/dy
    ykaon(m) = ykaon(m)/float(nrun)/dy
!       WRITE(1092,'(E11.3,2X,E11.3,2X,E11.3)')m*DY,YKAON(M),DYKAON
  End Do
  Return
End Subroutine flow
!bali1/16/99
!*******************************************
! Purpose: pp_bar annihilation cross section as a functon of their cms energy
!      real*4 function xppbar(srt)
Real Function xppbar(srt)
!  srt    = DSQRT(s) in GeV                                                   *
!  xppbar = pp_bar annihilation cross section in mb                           *
!
!  Reference: G.J. Wang, R. Bellwied, C. Pruneau and G. Welke
!             Proc. of the 14th Winter Workshop on Nuclear Dynamics,
!             Snowbird, Utah 31, Eds. W. Bauer and H.G. Ritter
!             (Plenum Publishing, 1998)                             *
!
!*****************************************
  Parameter (pmass=0.9383, xmax=400.)
  Save
! Note:
! (1) we introduce a new parameter xmax=400 mb:
!     the maximum annihilation xsection
! there are shadowing effects in pp_bar annihilation, with this parameter
! we can probably look at these effects
! (2) Calculate p(lab) from srt [GeV], since the formular in the
! reference applies only to the case of a p_bar on a proton at rest
! Formula used: srt**2=2.*pmass*(pmass+sqrt(pmass**2+plab**2))
  xppbar = 1.E-06
  plab2 = (srt**2/(2.*pmass)-pmass)**2 - pmass**2
  If (plab2>0) Then
    plab = sqrt(plab2)
    xppbar = 67./(plab**0.7)
    If (xppbar>xmax) xppbar = xmax
  End If
  Return
End Function xppbar
!bali1/16/99 end
!*********************************
!bali2/6/99
!*******************************************
! Purpose: To generate randomly the no. of pions in the final
!          state of pp_bar annihilation according to a statistical
!          model by using of the rejection method.
!bz2/25/99
!      real*4 function pbarfs(srt,npion,iseed)
Subroutine pbarfs(srt, npion, iseed)
!bz2/25/99end
! Quantities:
!  srt: DSQRT(s) in GeV                                                    *
!  npion: No. of pions produced in the annihilation of ppbar at srt        *
!  nmax=6, cutoff of the maximum no. of n the code can handle
!
!  Reference: C.M. Ko and R. Yuan, Phys. Lett. B192 (1987) 31      *
!
!*****************************************
  Parameter (pimass=0.140, pi=3.1415926)
  Dimension factor(6), pnpi(6)
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
! the factorial coefficients in the pion no. distribution
! from n=2 to 6 calculated use the formula in the reference
  factor(2) = 1.
  factor(3) = 1.17E-01
  factor(4) = 3.27E-03
  factor(5) = 3.58E-05
  factor(6) = 1.93E-07
  ene = (srt/pimass)**3/(6.*pi**2)
! the relative probability from n=2 to 6
  Do n = 2, 6
    pnpi(n) = ene**n*factor(n)
  End Do
! find the maximum of the probabilities, I checked a
! Fortan manual: max() returns the maximum value of
! the same type as in the argument list
  pmax = max(pnpi(2), pnpi(3), pnpi(4), pnpi(5), pnpi(6))
! randomly generate n between 2 and 6
  ntry = 0
  10 npion = 2 + int(5*ranart(nseed))
!lin-4/2008 check bounds:
  If (npion>6) Goto 10
  thisp = pnpi(npion)/pmax
  ntry = ntry + 1
! decide whether to take this npion according to the distribution
! using rejection method.
  If ((thisp<ranart(nseed)) .And. (ntry<=20)) Goto 10
! now take the last generated npion and return
  Return
End Subroutine pbarfs
!*********************************
!bali2/6/99 end
!bz3/9/99 kkbar
!bali3/5/99
!*****************************************
! purpose: Xsection for K+ K- to pi+ pi-
!      real*4 function xkkpi(srt)
!  srt    = DSQRT(s) in GeV                                  *
!  xkkpi   = xsection in mb obtained from
!           the detailed balance                             *
! ******************************************
!          parameter (pimass=0.140,aka=0.498)
!       xkkpi=1.e-08
!       ppi2=(srt/2)**2-pimass**2
!       pk2=(srt/2)**2-aka**2
!       if(ppi2.le.0.or.pk2.le.0)return
!bz3/9/99 kkbar
!       xkkpi=ppi2/pk2*pipik(srt)
!       xkkpi=9.0 / 4.0 * ppi2/pk2*pipik(srt)
!        xkkpi = 2.0 * xkkpi
!bz3/9/99 kkbar end

!bz3/9/99 kkbar
!       end
!       return
!        END
!bz3/9/99 kkbar end

!bali3/5/99 end
!bz3/9/99 kkbar end

!bz3/9/99 kkbar
!****************************
! purpose: Xsection for K+ K- to pi+ pi-
Subroutine xkkann(srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8, xsk9, xsk10, xsk11, sigk, rrkk)
!  srt    = DSQRT(s) in GeV                                       *
!  xsk1   = annihilation into pi pi                               *
!  xsk2   = annihilation into pi rho (shifted to XKKSAN)         *
!  xsk3   = annihilation into pi omega (shifted to XKKSAN)       *
!  xsk4   = annihilation into pi eta                              *
!  xsk5   = annihilation into rho rho                             *
!  xsk6   = annihilation into rho omega                           *
!  xsk7   = annihilation into rho eta (shifted to XKKSAN)        *
!  xsk8   = annihilation into omega omega                         *
!  xsk9   = annihilation into omega eta (shifted to XKKSAN)      *
!  xsk10  = annihilation into eta eta                             *
!  sigk   = xsection in mb obtained from                          *
!           the detailed balance                                  *
! ***************************
  Parameter (maxstr=150001, maxx=20, maxz=24)
  Parameter (aka=0.498, pimass=0.140, rhom=0.770, omegam=0.7819, etam=0.5473, aphi=1.02)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
!c      SAVE /DD/
  Save

  s = srt**2
  sigk = 1.E-08
  xsk1 = 0.0
  xsk2 = 0.0
  xsk3 = 0.0
  xsk4 = 0.0
  xsk5 = 0.0
  xsk6 = 0.0
  xsk7 = 0.0
  xsk8 = 0.0
  xsk9 = 0.0
  xsk10 = 0.0
  xsk11 = 0.0

  xpion0 = pipik(srt)
!.....take into account both K+ and K0
  xpion0 = 2.0*xpion0
  pi2 = s*(s-4.0*aka**2)
  If (pi2<=0.0) Return

  xm1 = pimass
  xm2 = pimass
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk1 = 9.0/4.0*pf2/pi2*xpion0
  End If

!lin-8/28/00 (pi eta) eta -> K+K- is assumed the same as pi pi -> K+K-:
  xm1 = pimass
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk4 = 3.0/4.0*pf2/pi2*xpion0
  End If

  xm1 = etam
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk10 = 1.0/4.0*pf2/pi2*xpion0
  End If

  xpion0 = rrkk

!lin-11/07/00: (pi eta) (rho omega) -> K* Kbar (or K*bar K) instead to K Kbar:
!        XM1 = PIMASS
!        XM2 = RHOM
!        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
!        IF (PF2 .GT. 0.0) THEN
!           XSK2 = 27.0 / 4.0 * PF2 / PI2 * XPION0
!        END IF

!        XM1 = PIMASS
!        XM2 = OMEGAM
!        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
!        IF (PF2 .GT. 0.0) THEN
!           XSK3 = 9.0 / 4.0 * PF2 / PI2 * XPION0
!        END IF

  xm1 = rhom
  xm2 = rhom
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk5 = 81.0/4.0*pf2/pi2*xpion0
  End If

  xm1 = rhom
  xm2 = omegam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk6 = 27.0/4.0*pf2/pi2*xpion0
  End If

!        XM1 = RHOM
!        XM2 = ETAM
!        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
!        IF (PF2 .GT. 0.0) THEN
!           XSK7 = 9.0 / 4.0 * PF2 / PI2 * XPION0
!        END IF

  xm1 = omegam
  xm2 = omegam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xsk8 = 9.0/4.0*pf2/pi2*xpion0
  End If

!        XM1 = OMEGAM
!        XM2 = ETAM
!        PF2 = (S - (XM1 + XM2) ** 2) * (S - (XM1 - XM2) ** 2)
!        IF (PF2 .GT. 0.0) THEN
!           XSK9 = 3.0 / 4.0 * PF2 / PI2 * XPION0
!        END IF

!* K+ + K- --> phi
  fwdp = 1.68*(aphi**2-4.*aka**2)**1.5/6./aphi/aphi

!lin-9/2012: check argument in sqrt():
  scheck = srt**2 - 4.0*aka**2
  If (scheck<=0) Then
    Write (99, *) 'scheck47: ', scheck
    Stop
  End If
  pkaon = 0.5*sqrt(scheck)
!          pkaon=0.5*sqrt(srt**2-4.0*aka**2)

  xsk11 = 30.*3.14159*0.1973**2*(aphi*fwdp)**2/((srt**2-aphi**2)**2+(aphi*fwdp)**2)/pkaon**2
!
  sigk = xsk1 + xsk2 + xsk3 + xsk4 + xsk5 + xsk6 + xsk7 + xsk8 + xsk9 + xsk10 + xsk11

  Return
End Subroutine xkkann
!bz3/9/99 kkbar end

!****************************
! purpose: Xsection for Phi + B
Subroutine xphib(lb1, lb2, em1, em2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, sigp)
!
! ***************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Parameter (aka=0.498, ala=1.1157, pimass=0.140, aphi=1.02)
  Parameter (arho=0.77)
  Save

  sigp = 1.E-08
  xsk1 = 0.0
  xsk2 = 0.0
  xsk3 = 0.0
  xsk4 = 0.0
  xsk5 = 0.0
  xsk6 = 0.0
  srrt = srt - (em1+em2)

!* phi + N(D) -> elastic scattering
!            XSK1 = 0.56  !! mb
!  !! mb  (photo-production xsecn used)
  xsk1 = 8.00
!
!* phi + N(D) -> pi + N
  If (srt>(ap1+amn)) Then
    xsk2 = 0.0235*srrt**(-0.519)
  End If
!
!* phi + N(D) -> pi + D
  If (srt>(ap1+am0)) Then
    If (srrt<0.7) Then
      xsk3 = 0.0119*srrt**(-0.534)
    Else
      xsk3 = 0.0130*srrt**(-0.304)
    End If
  End If
!
!* phi + N(D) -> rho + N
  If (srt>(arho+amn)) Then
    If (srrt<0.7) Then
      xsk4 = 0.0166*srrt**(-0.786)
    Else
      xsk4 = 0.0189*srrt**(-0.277)
    End If
  End If
!
!* phi + N(D) -> rho + D   (same as pi + D)
  If (srt>(arho+am0)) Then
    If (srrt<0.7) Then
      xsk5 = 0.0119*srrt**(-0.534)
    Else
      xsk5 = 0.0130*srrt**(-0.304)
    End If
  End If
!
!* phi + N -> K+ + La
  If ((lb1>=1 .And. lb1<=2) .Or. (lb2>=1 .And. lb2<=2)) Then
    If (srt>(aka+ala)) Then
      xsk6 = 1.715/((srrt+3.508)**2-12.138)
    End If
  End If
  sigp = xsk1 + xsk2 + xsk3 + xsk4 + xsk5 + xsk6
  Return
End Subroutine xphib
!
!*********************************
!
Subroutine crphib(px, py, pz, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, sigp, iblock)
!
!     PURPOSE:                                                         *
!             DEALING WITH PHI + N(D) --> pi+N(D), rho+N(D),  K+ + La
!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - INFORMATION about the reaction channel          *
!
!             iblock   - 20  elastic
!             iblock   - 221  K+ formation
!             iblock   - 223  others
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, arho=0.77)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
  px0 = px
  py0 = py
  pz0 = pz
  iblock = 223
!
  x1 = ranart(nseed)*sigp
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  xsk5 = xsk4 + xsk5
!
!  !! elastic scatt.
  If (x1<=xsk1) Then
    iblock = 20
    Goto 100
  Else If (x1<=xsk2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = ap1
    e(i2) = amn
    Goto 100
  Else If (x1<=xsk3) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = ap1
    e(i2) = am0
    Goto 100
  Else If (x1<=xsk4) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 1 + int(2*ranart(nseed))
    e(i1) = arho
    e(i2) = amn
    Goto 100
  Else If (x1<=xsk5) Then
    lb(i1) = 25 + int(3*ranart(nseed))
    lb(i2) = 6 + int(4*ranart(nseed))
    e(i1) = arho
    e(i2) = am0
    Goto 100
  Else
    lb(i1) = 23
    lb(i2) = 14
    e(i1) = aka
    e(i2) = ala
    iblock = 221
  End If
  100 Continue
  em1 = e(i1)
  em2 = e(i2)
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-08
  pr = sqrt(pr2)/(2.*srt)
! WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crphib
!
!****************************
! purpose: Xsection for Phi + B
!!! in fm^2
Subroutine pibphi(srt, lb1, lb2, em1, em2, xphi, xphin)
!
!      phi + N(D) <- pi + N
!      phi + N(D) <- pi + D
!      phi + N(D) <- rho + N
!      phi + N(D) <- rho + D   (same as pi + D)
!
! ***************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
  Parameter (aka=0.498, ala=1.1157, pimass=0.140, aphi=1.02)
  Parameter (arho=0.77)
  Save

  xphi = 0.0
  xphin = 0.0
  xphid = 0.0
!
  If ((lb1>=3 .And. lb1<=5) .Or. (lb2>=3 .And. lb2<=5)) Then
!
    If ((iabs(lb1)>=1 .And. iabs(lb1)<=2) .Or. (iabs(lb2)>=1 .And. iabs(lb2)<=2)) Then
!* phi + N <- pi + N
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        sig = 0.0235*srrt**(-0.519)
        xphin = sig*1.*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
!* phi + D <- pi + N
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        sig = 0.0235*srrt**(-0.519)
        xphid = sig*4.*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    Else
!* phi + N <- pi + D
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphin = sig*(1./4.)*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
!* phi + D <- pi + D
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphid = sig*1.*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    End If
!
!
!** for rho + N(D) colln
!
  Else
!
    If ((iabs(lb1)>=1 .And. iabs(lb1)<=2) .Or. (iabs(lb2)>=1 .And. iabs(lb2)<=2)) Then
!
!* phi + N <- rho + N
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        If (srrt<0.7) Then
          sig = 0.0166*srrt**(-0.786)
        Else
          sig = 0.0189*srrt**(-0.277)
        End If
        xphin = sig*(1./3.)*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
!* phi + D <- rho + N
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        If (srrt<0.7) Then
          sig = 0.0166*srrt**(-0.786)
        Else
          sig = 0.0189*srrt**(-0.277)
        End If
        xphid = sig*(4./3.)*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    Else
!* phi + N <- rho + D  (same as pi+D->phi+N)
      If (srt>(aphi+amn)) Then
        srrt = srt - (aphi+amn)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphin = sig*(1./12.)*(srt**2-(aphi+amn)**2)*(srt**2-(aphi-amn)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
!* phi + D <- rho + D  (same as pi+D->phi+D)
      If (srt>(aphi+am0)) Then
        srrt = srt - (aphi+am0)
        If (srrt<0.7) Then
          sig = 0.0119*srrt**(-0.534)
        Else
          sig = 0.0130*srrt**(-0.304)
        End If
        xphid = sig*(1./3.)*(srt**2-(aphi+am0)**2)*(srt**2-(aphi-am0)**2)/(srt**2-(em1+em2)**2)/(srt**2-(em1-em2)**2)
      End If
    End If
  End If
!   !! in fm^2
  xphin = xphin/10.
!   !! in fm^2
  xphid = xphid/10.
  xphi = xphin + xphid

  Return
End Subroutine pibphi
!
!****************************
! purpose: Xsection for phi +M to K+K etc
Subroutine phimes(i1, i2, srt, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, sigphi)

!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - THE INFORMATION BACK                            *
!                      223 --> phi destruction
!                      20 -->  elastic
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895, aomega=0.7819, arho=0.77, aphi=1.02)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (maxx=20, maxz=24)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /dd/rho(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhop(-maxx:maxx, -maxx:maxx, -maxz:maxz), rhon(-maxx:maxx, -maxx:maxx, -maxz:maxz)
!c      SAVE /DD/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Save

  s = srt**2
  sigphi = 1.E-08
  xsk1 = 0.0
  xsk2 = 0.0
  xsk3 = 0.0
  xsk4 = 0.0
  xsk5 = 0.0
  xsk6 = 0.0
  xsk7 = 0.0
  em1 = e(i1)
  em2 = e(i2)
  lb1 = lb(i1)
  lb2 = lb(i2)
  akap = aka
!******
!
!   !! mb, elastic
  xsk1 = 5.0

!lin-9/2012: check argument in sqrt():
  scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
  If (scheck<=0) Then
    Write (99, *) 'scheck48: ', scheck
    Stop
  End If
  pii = sqrt(scheck)
!           pii = sqrt((S-(em1+em2)**2)*(S-(em1-em2)**2))

! phi + K(-bar) channel
  If (lb1==23 .Or. lb2==23 .Or. lb1==21 .Or. lb2==21) Then
    If (srt>(ap1+akap)) Then
!             XSK2 = 2.5
      pff = sqrt((s-(ap1+akap)**2)*(s-(ap1-akap)**2))
      xsk2 = 195.639*pff/pii/32./pi/s
    End If
    If (srt>(arho+akap)) Then
!              XSK3 = 3.5
      pff = sqrt((s-(arho+akap)**2)*(s-(arho-akap)**2))
      xsk3 = 526.702*pff/pii/32./pi/s
    End If
    If (srt>(aomega+akap)) Then
!               XSK4 = 3.5
      pff = sqrt((s-(aomega+akap)**2)*(s-(aomega-akap)**2))
      xsk4 = 355.429*pff/pii/32./pi/s
    End If
    If (srt>(ap1+aks)) Then
!           XSK5 = 15.0
      pff = sqrt((s-(ap1+aks)**2)*(s-(ap1-aks)**2))
      xsk5 = 2047.042*pff/pii/32./pi/s
    End If
    If (srt>(arho+aks)) Then
!            XSK6 = 3.5
      pff = sqrt((s-(arho+aks)**2)*(s-(arho-aks)**2))
      xsk6 = 1371.257*pff/pii/32./pi/s
    End If
    If (srt>(aomega+aks)) Then
!            XSK7 = 3.5
      pff = sqrt((s-(aomega+aks)**2)*(s-(aomega-aks)**2))
      xsk7 = 482.292*pff/pii/32./pi/s
    End If
!
  Else If (iabs(lb1)==30 .Or. iabs(lb2)==30) Then
! phi + K*(-bar) channel
!
    If (srt>(ap1+akap)) Then
!             XSK2 = 3.5
      pff = sqrt((s-(ap1+akap)**2)*(s-(ap1-akap)**2))
      xsk2 = 372.378*pff/pii/32./pi/s
    End If
    If (srt>(arho+akap)) Then
!              XSK3 = 9.0
      pff = sqrt((s-(arho+akap)**2)*(s-(arho-akap)**2))
      xsk3 = 1313.960*pff/pii/32./pi/s
    End If
    If (srt>(aomega+akap)) Then
!               XSK4 = 6.5
      pff = sqrt((s-(aomega+akap)**2)*(s-(aomega-akap)**2))
      xsk4 = 440.558*pff/pii/32./pi/s
    End If
    If (srt>(ap1+aks)) Then
!           XSK5 = 30.0 !wrong
      pff = sqrt((s-(ap1+aks)**2)*(s-(ap1-aks)**2))
      xsk5 = 1496.692*pff/pii/32./pi/s
    End If
    If (srt>(arho+aks)) Then
!            XSK6 = 9.0
      pff = sqrt((s-(arho+aks)**2)*(s-(arho-aks)**2))
      xsk6 = 6999.840*pff/pii/32./pi/s
    End If
    If (srt>(aomega+aks)) Then
!            XSK7 = 15.0
      pff = sqrt((s-(aomega+aks)**2)*(s-(aomega-aks)**2))
      xsk7 = 1698.903*pff/pii/32./pi/s
    End If
  Else
!
! phi + rho(pi,omega) channel
!
    srr1 = em1 + em2
    If (srt>(akap+akap)) Then
      srrt = srt - srr1
!c          if(srrt .lt. 0.3)then
      If (srrt<0.3 .And. srrt>0.01) Then
        xsk2 = 1.69/(srrt**0.141-0.407)
      Else
        xsk2 = 3.74 + 0.008*srrt**1.9
      End If
    End If
    If (srt>(akap+aks)) Then
      srr2 = akap + aks
      srr = amax1(srr1, srr2)
      srrt = srt - srr
!c          if(srrt .lt. 0.3)then
      If (srrt<0.3 .And. srrt>0.01) Then
        xsk3 = 1.69/(srrt**0.141-0.407)
      Else
        xsk3 = 3.74 + 0.008*srrt**1.9
      End If
    End If
    If (srt>(aks+aks)) Then
      srr2 = aks + aks
      srr = amax1(srr1, srr2)
      srrt = srt - srr
!c          if(srrt .lt. 0.3)then
      If (srrt<0.3 .And. srrt>0.01) Then
        xsk4 = 1.69/(srrt**0.141-0.407)
      Else
        xsk4 = 3.74 + 0.008*srrt**1.9
      End If
    End If
!          xsk2 = amin1(20.,xsk2)
!          xsk3 = amin1(20.,xsk3)
!          xsk4 = amin1(20.,xsk4)
  End If

  sigphi = xsk1 + xsk2 + xsk3 + xsk4 + xsk5 + xsk6 + xsk7

  Return
End Subroutine phimes

!*********************************
!     PURPOSE:                                                         *
!             DEALING WITH phi+M  scatt.
!
Subroutine crphim(px, py, pz, srt, i1, i2, xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, sigphi, ikkg, ikkl, iblock)
!
!     QUANTITIES:                                                      *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - THE INFORMATION BACK                            *
!                      20 -->  elastic
!                      223 --> phi + pi(rho,omega)
!                      224 --> phi + K -> K + pi(rho,omega)
!                      225 --> phi + K -> K* + pi(rho,omega)
!                      226 --> phi + K* -> K + pi(rho,omega)
!                      227 --> phi + K* -> K* + pi(rho,omega)
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, arho=0.77, aomega=0.7819, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, aks=0.895)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
  px0 = px
  py0 = py
  pz0 = pz
  lb1 = lb(i1)
  lb2 = lb(i2)

  x1 = ranart(nseed)*sigphi
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  xsk5 = xsk4 + xsk5
  xsk6 = xsk5 + xsk6
  If (x1<=xsk1) Then
!        !! elastic scatt
    iblock = 20
    Goto 100
  Else
!
!phi + (K,K*)-bar
    If (lb1==23 .Or. lb1==21 .Or. iabs(lb1)==30 .Or. lb2==23 .Or. lb2==21 .Or. iabs(lb2)==30) Then
!
      If (lb1==23 .Or. lb2==23) Then
        ikkl = 1
        iblock = 224
        iad1 = 23
        iad2 = 30
      Else If (lb1==30 .Or. lb2==30) Then
        ikkl = 0
        iblock = 226
        iad1 = 23
        iad2 = 30
      Else If (lb1==21 .Or. lb2==21) Then
        ikkl = 1
        iblock = 124
        iad1 = 21
        iad2 = -30
!         !! -30
      Else
        ikkl = 0
        iblock = 126
        iad1 = 21
        iad2 = -30
      End If
      If (x1<=xsk2) Then
        lb(i1) = 3 + int(3*ranart(nseed))
        lb(i2) = iad1
        e(i1) = ap1
        e(i2) = aka
        ikkg = 1
        Goto 100
      Else If (x1<=xsk3) Then
        lb(i1) = 25 + int(3*ranart(nseed))
        lb(i2) = iad1
        e(i1) = arho
        e(i2) = aka
        ikkg = 1
        Goto 100
      Else If (x1<=xsk4) Then
        lb(i1) = 28
        lb(i2) = iad1
        e(i1) = aomega
        e(i2) = aka
        ikkg = 1
        Goto 100
      Else If (x1<=xsk5) Then
        lb(i1) = 3 + int(3*ranart(nseed))
        lb(i2) = iad2
        e(i1) = ap1
        e(i2) = aks
        ikkg = 0
        iblock = iblock + 1
        Goto 100
      Else If (x1<=xsk6) Then
        lb(i1) = 25 + int(3*ranart(nseed))
        lb(i2) = iad2
        e(i1) = arho
        e(i2) = aks
        ikkg = 0
        iblock = iblock + 1
        Goto 100
      Else
        lb(i1) = 28
        lb(i2) = iad2
        e(i1) = aomega
        e(i2) = aks
        ikkg = 0
        iblock = iblock + 1
        Goto 100
      End If
    Else
!      !! phi destruction via (pi,rho,omega)
      iblock = 223
!phi + pi(rho,omega)
      If (x1<=xsk2) Then
        lb(i1) = 23
        lb(i2) = 21
        e(i1) = aka
        e(i2) = aka
        ikkg = 2
        ikkl = 0
        Goto 100
      Else If (x1<=xsk3) Then
        lb(i1) = 23
!           LB(I2) = 30
        lb(i2) = -30
!lin-2/10/03 currently take XSK3 to be the sum of KK*bar & KbarK*:
        If (ranart(nseed)<=0.5) Then
          lb(i1) = 21
          lb(i2) = 30
        End If

        e(i1) = aka
        e(i2) = aks
        ikkg = 1
        ikkl = 0
        Goto 100
      Else If (x1<=xsk4) Then
        lb(i1) = 30
!           LB(I2) = 30
        lb(i2) = -30
        e(i1) = aks
        e(i2) = aks
        ikkg = 0
        ikkl = 0
        Goto 100
      End If
    End If
  End If
!
  100 Continue
  em1 = e(i1)
  em2 = e(i2)

!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-08
  pr = sqrt(pr2)/(2.*srt)
! WE ASSUME AN ISOTROPIC ANGULAR DISTRIBUTION IN THE CMS
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)
  Return
End Subroutine crphim
!*********************************
!*********************************
!bz3/9/99 khyperon
!************************************
! purpose: Xsection for K+Y ->  piN                                       *
!          Xsection for K+Y-bar ->  piN-bar   !! sp03/29/01               *
!
Subroutine xkhype(i1, i2, srt, xky1, xky2, xky3, xky4, xky5, xky6, xky7, xky8, xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17, sigk)
!      subroutine xkhype(i1, i2, srt, sigk)
!  srt    = DSQRT(s) in GeV                                               *
!  xkkpi   = xsection in mb obtained from                                 *
!           the detailed balance                                          *
! ***********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, amrho=0.769, amomga=0.782, aphi=1.02, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (pimass=0.140, ameta=0.5473, aka=0.498, aml=1.116, ams=1.193, am1440=1.44, am1535=1.535)
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Save

  s = srt**2
  sigk = 1.E-08
  xky1 = 0.0
  xky2 = 0.0
  xky3 = 0.0
  xky4 = 0.0
  xky5 = 0.0
  xky6 = 0.0
  xky7 = 0.0
  xky8 = 0.0
  xky9 = 0.0
  xky10 = 0.0
  xky11 = 0.0
  xky12 = 0.0
  xky13 = 0.0
  xky14 = 0.0
  xky15 = 0.0
  xky16 = 0.0
  xky17 = 0.0

  lb1 = lb(i1)
  lb2 = lb(i2)
  If (iabs(lb1)==14 .Or. iabs(lb2)==14) Then
    xkaon0 = pnlka(srt)
    xkaon0 = 2.0*xkaon0
    pi2 = (s-(aml+aka)**2)*(s-(aml-aka)**2)
  Else
    xkaon0 = pnska(srt)
    xkaon0 = 2.0*xkaon0
    pi2 = (s-(ams+aka)**2)*(s-(ams-aka)**2)
  End If
  If (pi2<=0.0) Return

  xm1 = pimass
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky1 = 3.0*pf2/pi2*xkaon0
  End If

  xm1 = pimass
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky2 = 12.0*pf2/pi2*xkaon0
  End If

  xm1 = pimass
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky3 = 3.0*pf2/pi2*xkaon0
  End If

  xm1 = pimass
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky4 = 3.0*pf2/pi2*xkaon0
  End If

  xm1 = amrho
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky5 = 9.0*pf2/pi2*xkaon0
  End If

  xm1 = amrho
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky6 = 36.0*pf2/pi2*xkaon0
  End If

  xm1 = amrho
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky7 = 9.0*pf2/pi2*xkaon0
  End If

  xm1 = amrho
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky8 = 9.0*pf2/pi2*xkaon0
  End If

  xm1 = amomga
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky9 = 3.0*pf2/pi2*xkaon0
  End If

  xm1 = amomga
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky10 = 12.0*pf2/pi2*xkaon0
  End If

  xm1 = amomga
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky11 = 3.0*pf2/pi2*xkaon0
  End If

  xm1 = amomga
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky12 = 3.0*pf2/pi2*xkaon0
  End If

  xm1 = ameta
  xm2 = amp
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky13 = 1.0*pf2/pi2*xkaon0
  End If

  xm1 = ameta
  xm2 = am0
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky14 = 4.0*pf2/pi2*xkaon0
  End If

  xm1 = ameta
  xm2 = am1440
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky15 = 1.0*pf2/pi2*xkaon0
  End If

  xm1 = ameta
  xm2 = am1535
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    xky16 = 1.0*pf2/pi2*xkaon0
  End If

!sp11/21/01  K+ + La --> phi + N
  If (lb1==14 .Or. lb2==14) Then
    If (srt>(aphi+amn)) Then
      srrt = srt - (aphi+amn)
      sig = 1.715/((srrt+3.508)**2-12.138)
      xm1 = amn
      xm2 = aphi
      pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
!     ! fm^-1
      xky17 = 3.0*pf2/pi2*sig/10.
    End If
  End If
!sp11/21/01  end
!

  If ((iabs(lb1)>=15 .And. iabs(lb1)<=17) .Or. (iabs(lb2)>=15 .And. iabs(lb2)<=17)) Then
    ddf = 3.0
    xky1 = xky1/ddf
    xky2 = xky2/ddf
    xky3 = xky3/ddf
    xky4 = xky4/ddf
    xky5 = xky5/ddf
    xky6 = xky6/ddf
    xky7 = xky7/ddf
    xky8 = xky8/ddf
    xky9 = xky9/ddf
    xky10 = xky10/ddf
    xky11 = xky11/ddf
    xky12 = xky12/ddf
    xky13 = xky13/ddf
    xky14 = xky14/ddf
    xky15 = xky15/ddf
    xky16 = xky16/ddf
  End If

  sigk = xky1 + xky2 + xky3 + xky4 + xky5 + xky6 + xky7 + xky8 + xky9 + xky10 + xky11 + xky12 + xky13 + xky14 + xky15 + xky16 + xky17

  Return
End Subroutine xkhype

!*******************************
Block Data ppbdat

  Parameter (amp=0.93828, amn=0.939457, am0=1.232, am1440=1.44, am1535=1.535)

!     to give default values to parameters for BbarB production from mesons
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
!c      SAVE /ppbmas/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save
!     thresh(i) gives the mass thresh for final channel i:
  Data thresh/1.87656, 1.877737, 1.878914, 2.17028, 2.171457, 2.37828, 2.379457, 2.464, 2.47328, 2.474457, 2.672, 2.767, 2.88, 2.975, 3.07/
!     ppbm(i,j=1,2) gives masses for the two final baryons of channel i,
!     with j=1 for the lighter baryon:
  Data (ppbm(i,1), i=1, 15)/amp, amp, amn, amp, amn, amp, amn, am0, amp, amn, am0, am0, am1440, am1440, am1535/
  Data (ppbm(i,2), i=1, 15)/amp, amn, amn, am0, am0, am1440, am1440, am0, am1535, am1535, am1440, am1535, am1440, am1535, am1535/
!     factr2(i) gives weights for producing i pions from ppbar annihilation:
  Data factr2/0, 1, 1.17E-01, 3.27E-03, 3.58E-05, 1.93E-07/
!     niso(i) gives the degeneracy factor for final channel i:
  Data niso/1, 2, 1, 16, 16, 4, 4, 64, 4, 4, 32, 32, 4, 8, 4/

End Block Data ppbdat


!****************************************
! get the number of BbarB states available for mm collisions of energy srt
Subroutine getnst(srt)
!  srt    = DSQRT(s) in GeV                                                   *
!****************************************
  Parameter (pimass=0.140, pi=3.1415926)
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
!c      SAVE /ppbmas/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  s = srt**2
  nstate = 0
  wtot = 0.
  If (srt<=thresh(1)) Return
  Do i = 1, 15
    weight(i) = 0.
    If (srt>thresh(i)) nstate = i
  End Do
  Do i = 1, nstate
    pf2 = (s-(ppbm(i,1)+ppbm(i,2))**2)*(s-(ppbm(i,1)-ppbm(i,2))**2)/4/s
    weight(i) = pf2*niso(i)
    wtot = wtot + weight(i)
  End Do
  ene = (srt/pimass)**3/(6.*pi**2)
  fsum = factr2(2) + factr2(3)*ene + factr2(4)*ene**2 + factr2(5)*ene**3 + factr2(6)*ene**4

  Return
End Subroutine getnst

!****************************************
! for pion+pion-->Bbar B                                                      *
!      real*4 function ppbbar(srt)
Real Function ppbbar(srt)
!****************************************
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  sppb2p = xppbar(srt)*factr2(2)/fsum
  pi2 = (s-4*pimass**2)/4
  ppbbar = 4./9.*sppb2p/pi2*wtot

  Return
End Function ppbbar

!****************************************
! for pion+rho-->Bbar B                                                      *
!      real*4 function prbbar(srt)
Real Function prbbar(srt)
!****************************************
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  sppb3p = xppbar(srt)*factr2(3)*ene/fsum
  pi2 = (s-(pimass+arho)**2)*(s-(pimass-arho)**2)/4/s
  prbbar = 4./27.*sppb3p/pi2*wtot

  Return
End Function prbbar

!****************************************
! for rho+rho-->Bbar B                                                      *
!      real*4 function rrbbar(srt)
Real Function rrbbar(srt)
!****************************************
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  sppb4p = xppbar(srt)*factr2(4)*ene**2/fsum
  pi2 = (s-4*arho**2)/4
  rrbbar = 4./81.*(sppb4p/2)/pi2*wtot

  Return
End Function rrbbar

!****************************************
! for pi+omega-->Bbar B                                                      *
!      real*4 function pobbar(srt)
Real Function pobbar(srt)
!****************************************
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  sppb4p = xppbar(srt)*factr2(4)*ene**2/fsum
  pi2 = (s-(pimass+aomega)**2)*(s-(pimass-aomega)**2)/4/s
  pobbar = 4./9.*(sppb4p/2)/pi2*wtot

  Return
End Function pobbar

!****************************************
! for rho+omega-->Bbar B                                                      *
!      real*4 function robbar(srt)
Real Function robbar(srt)
!****************************************
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  sppb5p = xppbar(srt)*factr2(5)*ene**3/fsum
  pi2 = (s-(arho+aomega)**2)*(s-(arho-aomega)**2)/4/s
  robbar = 4./27.*sppb5p/pi2*wtot

  Return
End Function robbar

!****************************************
! for omega+omega-->Bbar B                                                    *
!      real*4 function oobbar(srt)
Real Function oobbar(srt)
!****************************************
  Parameter (pimass=0.140, arho=0.77, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  sppb6p = xppbar(srt)*factr2(6)*ene**4/fsum
  pi2 = (s-4*aomega**2)/4
  oobbar = 4./9.*sppb6p/pi2*wtot

  Return
End Function oobbar

!****************************************
! Generate final states for mm-->Bbar B                                       *
Subroutine bbarfs(lbb1, lbb2, ei1, ei2, iblock, iseed)
!****************************************
  Common /ppbmas/niso(15), nstate, ppbm(15, 2), thresh(15), weight(15)
!c      SAVE /ppbmas/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

!     determine which final BbarB channel occurs:
  rd = ranart(nseed)
  wsum = 0.
  Do i = 1, nstate
    wsum = wsum + weight(i)
    If (rd<=(wsum/wtot)) Then
      ifs = i
      ei1 = ppbm(i, 1)
      ei2 = ppbm(i, 2)
      Goto 10
    End If
  End Do
  10 Continue

!1    pbar p
  If (ifs==1) Then
    iblock = 1801
    lbb1 = -1
    lbb2 = 1
  Else If (ifs==2) Then
!2    pbar n
    If (ranart(nseed)<=0.5) Then
      iblock = 18021
      lbb1 = -1
      lbb2 = 2
!2    nbar p
    Else
      iblock = 18022
      lbb1 = 1
      lbb2 = -2
    End If
!3    nbar n
  Else If (ifs==3) Then
    iblock = 1803
    lbb1 = -2
    lbb2 = 2
!4&5  (pbar nbar) Delta, (p n) anti-Delta
  Else If (ifs==4 .Or. ifs==5) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
!     (pbar nbar) Delta
      If (ifs==4) Then
        iblock = 18041
        lbb1 = -1
      Else
        iblock = 18051
        lbb1 = -2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.25) Then
        lbb2 = 6
      Else If (rd2<=0.5) Then
        lbb2 = 7
      Else If (rd2<=0.75) Then
        lbb2 = 8
      Else
        lbb2 = 9
      End If
    Else
!     (p n) anti-Delta
      If (ifs==4) Then
        iblock = 18042
        lbb1 = 1
      Else
        iblock = 18052
        lbb1 = 2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.25) Then
        lbb2 = -6
      Else If (rd2<=0.5) Then
        lbb2 = -7
      Else If (rd2<=0.75) Then
        lbb2 = -8
      Else
        lbb2 = -9
      End If
    End If
!6&7  (pbar nbar) N*(1440), (p n) anti-N*(1440)
  Else If (ifs==6 .Or. ifs==7) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
!     (pbar nbar) N*(1440)
      If (ifs==6) Then
        iblock = 18061
        lbb1 = -1
      Else
        iblock = 18071
        lbb1 = -2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = 10
      Else
        lbb2 = 11
      End If
    Else
!     (p n) anti-N*(1440)
      If (ifs==6) Then
        iblock = 18062
        lbb1 = 1
      Else
        iblock = 18072
        lbb1 = 2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = -10
      Else
        lbb2 = -11
      End If
    End If
!8    Delta anti-Delta
  Else If (ifs==8) Then
    iblock = 1808
    rd1 = ranart(nseed)
    If (rd1<=0.25) Then
      lbb1 = 6
    Else If (rd1<=0.5) Then
      lbb1 = 7
    Else If (rd1<=0.75) Then
      lbb1 = 8
    Else
      lbb1 = 9
    End If
    rd2 = ranart(nseed)
    If (rd2<=0.25) Then
      lbb2 = -6
    Else If (rd2<=0.5) Then
      lbb2 = -7
    Else If (rd2<=0.75) Then
      lbb2 = -8
    Else
      lbb2 = -9
    End If
!9&10 (pbar nbar) N*(1535), (p n) anti-N*(1535)
  Else If (ifs==9 .Or. ifs==10) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
!     (pbar nbar) N*(1440)
      If (ifs==9) Then
        iblock = 18091
        lbb1 = -1
      Else
        iblock = 18101
        lbb1 = -2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = 12
      Else
        lbb2 = 13
      End If
    Else
!     (p n) anti-N*(1535)
      If (ifs==9) Then
        iblock = 18092
        lbb1 = 1
      Else
        iblock = 18102
        lbb1 = 2
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = -12
      Else
        lbb2 = -13
      End If
    End If
!11&12 anti-Delta N*, Delta anti-N*
  Else If (ifs==11 .Or. ifs==12) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
!     anti-Delta N*
      rd1 = ranart(nseed)
      If (rd1<=0.25) Then
        lbb1 = -6
      Else If (rd1<=0.5) Then
        lbb1 = -7
      Else If (rd1<=0.75) Then
        lbb1 = -8
      Else
        lbb1 = -9
      End If
      If (ifs==11) Then
        iblock = 18111
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = 10
        Else
          lbb2 = 11
        End If
      Else
        iblock = 18121
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = 12
        Else
          lbb2 = 13
        End If
      End If
    Else
!     Delta anti-N*
      rd1 = ranart(nseed)
      If (rd1<=0.25) Then
        lbb1 = 6
      Else If (rd1<=0.5) Then
        lbb1 = 7
      Else If (rd1<=0.75) Then
        lbb1 = 8
      Else
        lbb1 = 9
      End If
      If (ifs==11) Then
        iblock = 18112
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = -10
        Else
          lbb2 = -11
        End If
      Else
        iblock = 18122
        rd2 = ranart(nseed)
        If (rd2<=0.5) Then
          lbb2 = -12
        Else
          lbb2 = -13
        End If
      End If
    End If
!13   N*(1440) anti-N*(1440)
  Else If (ifs==13) Then
    iblock = 1813
    rd1 = ranart(nseed)
    If (rd1<=0.5) Then
      lbb1 = 10
    Else
      lbb1 = 11
    End If
    rd2 = ranart(nseed)
    If (rd2<=0.5) Then
      lbb2 = -10
    Else
      lbb2 = -11
    End If
!14   anti-N*(1440) N*(1535), N*(1440) anti-N*(1535)
  Else If (ifs==14) Then
    rd = ranart(nseed)
    If (rd<=0.5) Then
!     anti-N*(1440) N*(1535)
      iblock = 18141
      rd1 = ranart(nseed)
      If (rd1<=0.5) Then
        lbb1 = -10
      Else
        lbb1 = -11
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = 12
      Else
        lbb2 = 13
      End If
    Else
!     N*(1440) anti-N*(1535)
      iblock = 18142
      rd1 = ranart(nseed)
      If (rd1<=0.5) Then
        lbb1 = 10
      Else
        lbb1 = 11
      End If
      rd2 = ranart(nseed)
      If (rd2<=0.5) Then
        lbb2 = -12
      Else
        lbb2 = -13
      End If
    End If
!15   N*(1535) anti-N*(1535)
  Else If (ifs==15) Then
    iblock = 1815
    rd1 = ranart(nseed)
    If (rd1<=0.5) Then
      lbb1 = 12
    Else
      lbb1 = 13
    End If
    rd2 = ranart(nseed)
    If (rd2<=0.5) Then
      lbb2 = -12
    Else
      lbb2 = -13
    End If
  Else
  End If

  Return
End Subroutine bbarfs

!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!****************************************
! for pi pi <-> rho rho cross sections
Subroutine spprr(lb1, lb2, srt)
  Parameter (arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  pprr = 0.
  If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
!     for now, rho mass taken to be the central value in these two processes
    If (srt>(2*arho)) pprr = ptor(srt)
  Else If ((lb1>=25 .And. lb1<=27) .And. (lb2>=25 .And. lb2<=27)) Then
    pprr = rtop(srt)
  End If
!
  Return
End Subroutine spprr

!****************************************
! for pi pi -> rho rho, determined from detailed balance
Real Function ptor(srt)
!****************************************
  Parameter (pimass=0.140, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  s2 = srt**2
  ptor = 9*(s2-4*arho**2)/(s2-4*pimass**2)*rtop(srt)

  Return
End Function ptor

!****************************************
! for rho rho -> pi pi, assumed a constant cross section (in mb)
Real Function rtop(srt)
!****************************************
  rtop = 5.
  Return
End Function rtop

!****************************************
! for pi pi <-> rho rho final states
Subroutine pi2ro2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (ap1=0.13496, ap2=0.13957)
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  If ((lb(i1)>=3 .And. lb(i1)<=5) .And. (lb(i2)>=3 .And. lb(i2)<=5)) Then
    iblock = 1850
    ei1 = 0.77
    ei2 = 0.77
!     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
!     thus the cross sections used are considered as the isospin-averaged ones.
    lbb1 = 25 + int(3*ranart(nseed))
    lbb2 = 25 + int(3*ranart(nseed))
  Else If ((lb(i1)>=25 .And. lb(i1)<=27) .And. (lb(i2)>=25 .And. lb(i2)<=27)) Then
    iblock = 1851
    lbb1 = 3 + int(3*ranart(nseed))
    lbb2 = 3 + int(3*ranart(nseed))
    ei1 = ap2
    ei2 = ap2
    If (lbb1==4) ei1 = ap1
    If (lbb2==4) ei2 = ap1
  End If

  Return
End Subroutine pi2ro2

!****************************************
! for pi pi <-> eta eta cross sections
Subroutine sppee(lb1, lb2, srt)
  Parameter (etam=0.5475)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  ppee = 0.
  If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
    If (srt>(2*etam)) ppee = ptoe(srt)
  Else If (lb1==0 .And. lb2==0) Then
    ppee = etop(srt)
  End If

  Return
End Subroutine sppee

!****************************************
! for pi pi -> eta eta, determined from detailed balance, spin-isospin averaged
Real Function ptoe(srt)
!****************************************
  Parameter (pimass=0.140, etam=0.5475)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  s2 = srt**2
  ptoe = 1./9.*(s2-4*etam**2)/(s2-4*pimass**2)*etop(srt)

  Return
End Function ptoe
!****************************************
! for eta eta -> pi pi, assumed a constant cross section (in mb)
Real Function etop(srt)
!****************************************

!     eta equilibration:
!     most important channel is found to be pi pi <-> pi eta, then
!     rho pi <-> rho eta.
  etop = 5.
  Return
End Function etop

!****************************************
! for pi pi <-> eta eta final states
Subroutine pi2et2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (ap1=0.13496, ap2=0.13957, etam=0.5475)
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  If ((lb(i1)>=3 .And. lb(i1)<=5) .And. (lb(i2)>=3 .And. lb(i2)<=5)) Then
    iblock = 1860
    ei1 = etam
    ei2 = etam
!     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
!     thus the cross sections used are considered as the isospin-averaged ones.
    lbb1 = 0
    lbb2 = 0
  Else If (lb(i1)==0 .And. lb(i2)==0) Then
    iblock = 1861
    lbb1 = 3 + int(3*ranart(nseed))
    lbb2 = 3 + int(3*ranart(nseed))
    ei1 = ap2
    ei2 = ap2
    If (lbb1==4) ei1 = ap1
    If (lbb2==4) ei2 = ap1
  End If

  Return
End Subroutine pi2et2

!****************************************
! for pi pi <-> pi eta cross sections
Subroutine spppe(lb1, lb2, srt)
  Parameter (pimass=0.140, etam=0.5475)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  pppe = 0.
  If ((lb1>=3 .And. lb1<=5) .And. (lb2>=3 .And. lb2<=5)) Then
    If (srt>(etam+pimass)) pppe = pptope(srt)
  Else If ((lb1>=3 .And. lb1<=5) .And. lb2==0) Then
    pppe = petopp(srt)
  Else If ((lb2>=3 .And. lb2<=5) .And. lb1==0) Then
    pppe = petopp(srt)
  End If

  Return
End Subroutine spppe

!****************************************
! for pi pi -> pi eta, determined from detailed balance, spin-isospin averaged
Real Function pptope(srt)
!****************************************
  Parameter (pimass=0.140, etam=0.5475)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  s2 = srt**2
  pf2 = (s2-(pimass+etam)**2)*(s2-(pimass-etam)**2)/2/sqrt(s2)
  pi2 = (s2-4*pimass**2)*s2/2/sqrt(s2)
  pptope = 1./3.*pf2/pi2*petopp(srt)

  Return
End Function pptope
!****************************************
! for pi eta -> pi pi, assumed a constant cross section (in mb)
Real Function petopp(srt)
!****************************************

!     eta equilibration:
  petopp = 5.
  Return
End Function petopp

!****************************************
! for pi pi <-> pi eta final states
Subroutine pi3eta(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (ap1=0.13496, ap2=0.13957, etam=0.5475)
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  If ((lb(i1)>=3 .And. lb(i1)<=5) .And. (lb(i2)>=3 .And. lb(i2)<=5)) Then
    iblock = 1870
    ei1 = ap2
    ei2 = etam
!     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
!     thus the cross sections used are considered as the isospin-averaged ones.
    lbb1 = 3 + int(3*ranart(nseed))
    If (lbb1==4) ei1 = ap1
    lbb2 = 0
  Else If ((lb(i1)>=3 .And. lb(i1)<=5 .And. lb(i2)==0) .Or. (lb(i2)>=3 .And. lb(i2)<=5 .And. lb(i1)==0)) Then
    iblock = 1871
    lbb1 = 3 + int(3*ranart(nseed))
    lbb2 = 3 + int(3*ranart(nseed))
    ei1 = ap2
    ei2 = ap2
    If (lbb1==4) ei1 = ap1
    If (lbb2==4) ei2 = ap1
  End If

  Return
End Subroutine pi3eta

!****************************************
! for rho pi <-> rho eta cross sections
Subroutine srpre(lb1, lb2, srt)
  Parameter (pimass=0.140, etam=0.5475, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  rpre = 0.
  If (lb1>=25 .And. lb1<=27 .And. lb2>=3 .And. lb2<=5) Then
    If (srt>(etam+arho)) rpre = rptore(srt)
  Else If (lb2>=25 .And. lb2<=27 .And. lb1>=3 .And. lb1<=5) Then
    If (srt>(etam+arho)) rpre = rptore(srt)
  Else If (lb1>=25 .And. lb1<=27 .And. lb2==0) Then
    If (srt>(pimass+arho)) rpre = retorp(srt)
  Else If (lb2>=25 .And. lb2<=27 .And. lb1==0) Then
    If (srt>(pimass+arho)) rpre = retorp(srt)
  End If

  Return
End Subroutine srpre

!****************************************
! for rho pi->rho eta, determined from detailed balance, spin-isospin averaged
Real Function rptore(srt)
!****************************************
  Parameter (pimass=0.140, etam=0.5475, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  s2 = srt**2
  pf2 = (s2-(arho+etam)**2)*(s2-(arho-etam)**2)/2/sqrt(s2)
  pi2 = (s2-(arho+pimass)**2)*(s2-(arho-pimass)**2)/2/sqrt(s2)
  rptore = 1./3.*pf2/pi2*retorp(srt)

  Return
End Function rptore
!****************************************
! for rho eta -> rho pi, assumed a constant cross section (in mb)
Real Function retorp(srt)
!****************************************

!     eta equilibration:
  retorp = 5.
  Return
End Function retorp

!****************************************
! for rho pi <-> rho eta final states
Subroutine rpiret(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (ap1=0.13496, ap2=0.13957, etam=0.5475, arho=0.77)
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  If ((lb(i1)>=25 .And. lb(i1)<=27 .And. lb(i2)>=3 .And. lb(i2)<=5) .Or. (lb(i1)>=3 .And. lb(i1)<=5 .And. lb(i2)>=25 .And. lb(i2)<=27)) Then
    iblock = 1880
    ei1 = arho
    ei2 = etam
!     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
!     thus the cross sections used are considered as the isospin-averaged ones.
    lbb1 = 25 + int(3*ranart(nseed))
    lbb2 = 0
  Else If ((lb(i1)>=25 .And. lb(i1)<=27 .And. lb(i2)==0) .Or. (lb(i2)>=25 .And. lb(i2)<=27 .And. lb(i1)==0)) Then
    iblock = 1881
    lbb1 = 25 + int(3*ranart(nseed))
    lbb2 = 3 + int(3*ranart(nseed))
    ei1 = arho
    ei2 = ap2
    If (lbb2==4) ei2 = ap1
  End If

  Return
End Subroutine rpiret

!****************************************
! for omega pi <-> omega eta cross sections
Subroutine sopoe(lb1, lb2, srt)
  Parameter (etam=0.5475, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  xopoe = 0.
  If ((lb1==28 .And. lb2>=3 .And. lb2<=5) .Or. (lb2==28 .And. lb1>=3 .And. lb1<=5)) Then
    If (srt>(aomega+etam)) xopoe = xop2oe(srt)
  Else If ((lb1==28 .And. lb2==0) .Or. (lb1==0 .And. lb2==28)) Then
    If (srt>(aomega+etam)) xopoe = xoe2op(srt)
  End If

  Return
End Subroutine sopoe

!****************************************
! for omega pi -> omega eta,
!     determined from detailed balance, spin-isospin averaged
Real Function xop2oe(srt)
!****************************************
  Parameter (pimass=0.140, etam=0.5475, aomega=0.782)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  s2 = srt**2
  pf2 = (s2-(aomega+etam)**2)*(s2-(aomega-etam)**2)/2/sqrt(s2)
  pi2 = (s2-(aomega+pimass)**2)*(s2-(aomega-pimass)**2)/2/sqrt(s2)
  xop2oe = 1./3.*pf2/pi2*xoe2op(srt)

  Return
End Function xop2oe
!****************************************
! for omega eta -> omega pi, assumed a constant cross section (in mb)
Real Function xoe2op(srt)
!****************************************

!     eta equilibration:
  xoe2op = 5.
  Return
End Function xoe2op

!****************************************
! for omega pi <-> omega eta final states
Subroutine opioet(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (ap1=0.13496, ap2=0.13957, etam=0.5475, aomega=0.782)
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  If ((lb(i1)>=3 .And. lb(i1)<=5 .And. lb(i2)==28) .Or. (lb(i2)>=3 .And. lb(i2)<=5 .And. lb(i1)==28)) Then
    iblock = 1890
    ei1 = aomega
    ei2 = etam
!     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
!     thus the cross sections used are considered as the isospin-averaged ones.
    lbb1 = 28
    lbb2 = 0
  Else If ((lb(i1)==28 .And. lb(i2)==0) .Or. (lb(i1)==0 .And. lb(i2)==28)) Then
    iblock = 1891
    lbb1 = 28
    lbb2 = 3 + int(3*ranart(nseed))
    ei1 = aomega
    ei2 = ap2
    If (lbb2==4) ei2 = ap1
  End If

  Return
End Subroutine opioet

!****************************************
! for rho rho <-> eta eta cross sections
Subroutine srree(lb1, lb2, srt)
  Parameter (etam=0.5475, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  rree = 0.
  If (lb1>=25 .And. lb1<=27 .And. lb2>=25 .And. lb2<=27) Then
    If (srt>(2*etam)) rree = rrtoee(srt)
  Else If (lb1==0 .And. lb2==0) Then
    If (srt>(2*arho)) rree = eetorr(srt)
  End If

  Return
End Subroutine srree

!****************************************
! for eta eta -> rho rho
!     determined from detailed balance, spin-isospin averaged
Real Function eetorr(srt)
!****************************************
  Parameter (etam=0.5475, arho=0.77)
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Save

  s2 = srt**2
  eetorr = 81.*(s2-4*arho**2)/(s2-4*etam**2)*rrtoee(srt)

  Return
End Function eetorr
!****************************************
! for rho rho -> eta eta, assumed a constant cross section (in mb)
Real Function rrtoee(srt)
!****************************************

!     eta equilibration:
  rrtoee = 5.
  Return
End Function rrtoee

!****************************************
! for rho rho <-> eta eta final states
Subroutine ro2et2(i1, i2, lbb1, lbb2, ei1, ei2, iblock, iseed)
  Parameter (maxstr=150001)
  Parameter (etam=0.5475, arho=0.77)
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /ppb1/ene, factr2(6), fsum, ppinnb, s, wtot
!c      SAVE /ppb1/
  Common /ppmm/pprr, ppee, pppe, rpre, xopoe, rree
!c      SAVE /ppmm/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  If (lb(i1)>=25 .And. lb(i1)<=27 .And. lb(i2)>=25 .And. lb(i2)<=27) Then
    iblock = 1895
    ei1 = etam
    ei2 = etam
!     for now, we don't check isospin states(allowing pi+pi+ & pi0pi0 -> 2rho)
!     thus the cross sections used are considered as the isospin-averaged ones.
    lbb1 = 0
    lbb2 = 0
  Else If (lb(i1)==0 .And. lb(i2)==0) Then
    iblock = 1896
    lbb1 = 25 + int(3*ranart(nseed))
    lbb2 = 25 + int(3*ranart(nseed))
    ei1 = arho
    ei2 = arho
  End If

  Return
End Subroutine ro2et2

!****************************
! purpose: Xsection for K* Kbar or K*bar K to pi(eta) rho(omega)
Subroutine xkksan(i1, i2, srt, sigks1, sigks2, sigks3, sigks4, sigk, prkk)
!  srt    = DSQRT(s) in GeV                                       *
!  sigk   = xsection in mb obtained from                          *
!           the detailed balance                                  *
! ***************************
  Parameter (aka=0.498, pimass=0.140, rhom=0.770, aks=0.895, omegam=0.7819, etam=0.5473)
  Parameter (maxstr=150001)
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Save

  s = srt**2
  sigks1 = 1.E-08
  sigks2 = 1.E-08
  sigks3 = 1.E-08
  sigks4 = 1.E-08

  xpion0 = prkk
!lin note that prkk is for pi (rho omega) -> K* Kbar (AND!) K*bar K:
  xpion0 = xpion0/2

!c
!        PI2 = (S - (aks + AKA) ** 2) * (S - (aks - AKA) ** 2)
  pi2 = (s-(e(i1)+e(i2))**2)*(s-(e(i1)-e(i2))**2)
  sigk = 1.E-08
  If (pi2<=0.0) Return

  xm1 = pimass
  xm2 = rhom
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pi2>0.0 .And. pf2>0.0) Then
    sigks1 = 27.0/4.0*pf2/pi2*xpion0
  End If

  xm1 = pimass
  xm2 = omegam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pi2>0.0 .And. pf2>0.0) Then
    sigks2 = 9.0/4.0*pf2/pi2*xpion0
  End If

  xm1 = rhom
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    sigks3 = 9.0/4.0*pf2/pi2*xpion0
  End If

  xm1 = omegam
  xm2 = etam
  pf2 = (s-(xm1+xm2)**2)*(s-(xm1-xm2)**2)
  If (pf2>0.0) Then
    sigks4 = 3.0/4.0*pf2/pi2*xpion0
  End If

  sigk = sigks1 + sigks2 + sigks3 + sigks4

  Return
End Subroutine xkksan

!*********************************
!     PURPOSE:                                                         *
!     assign final states for KK*bar or K*Kbar --> light mesons
!
!      SUBROUTINE Crkspi(PX,PY,PZ,SRT,I1,I2,IBLOCK)
Subroutine crkspi(i1, i2, xsk1, xsk2, xsk3, xsk4, sigk, iblock, lbp1, lbp2, emm1, emm2)
!             iblock   - 466
!*********************************
  Parameter (maxstr=150001, maxr=1)
  Parameter (ap1=0.13496, ap2=0.13957, rhom=0.770, pi=3.1415926)
  Parameter (aeta=0.548, amomga=0.782)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  iblock = 466
! charges of final state mesons:

  x1 = ranart(nseed)*sigk
  xsk2 = xsk1 + xsk2
  xsk3 = xsk2 + xsk3
  xsk4 = xsk3 + xsk4
  If (x1<=xsk1) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = ap2
    e(i2) = rhom
  Else If (x1<=xsk2) Then
    lb(i1) = 3 + int(3*ranart(nseed))
    lb(i2) = 28
    e(i1) = ap2
    e(i2) = amomga
  Else If (x1<=xsk3) Then
    lb(i1) = 0
    lb(i2) = 25 + int(3*ranart(nseed))
    e(i1) = aeta
    e(i2) = rhom
  Else
    lb(i1) = 0
    lb(i2) = 28
    e(i1) = aeta
    e(i2) = amomga
  End If

  If (lb(i1)==4) e(i1) = ap1
  lbp1 = lb(i1)
  lbp2 = lb(i2)
  emm1 = e(i1)
  emm2 = e(i2)

  Return
End Subroutine crkspi

!---------------------------------------------------------------------------
! PURPOSE : CALCULATE THE MASS AND MOMENTUM OF K* RESONANCE
!           AFTER PION + KAON COLLISION
!clin only here the K* mass may be different from aks=0.895
Subroutine ksreso(i1, i2)
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926)
!lin-9/2012: improve precision for argument in sqrt():
  Double Precision e10, e20, scheck, p1, p2, p3
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
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
  Save
! 1. DETERMINE THE MOMENTUM COMPONENT OF THE K* IN THE CMS OF PI-K FRAME
!    WE LET I1 TO BE THE K* AND ABSORB I2

!lin-9/2012: improve precision for argument in sqrt():
!        E10=SQRT(E(I1)**2+P(1,I1)**2+P(2,I1)**2+P(3,I1)**2)
!        E20=SQRT(E(I2)**2+P(1,I2)**2+P(2,I2)**2+P(3,I2)**2)
  e10 = dsqrt(dble(e(i1))**2+dble(p(1,i1))**2+dble(p(2,i1))**2+dble(p(3,i1))**2)
  e20 = dsqrt(dble(e(i2))**2+dble(p(1,i2))**2+dble(p(2,i2))**2+dble(p(3,i2))**2)
  p1 = dble(p(1,i1)) + dble(p(1,i2))
  p2 = dble(p(2,i1)) + dble(p(2,i2))
  p3 = dble(p(3,i1)) + dble(p(3,i2))

  If (lb(i2)==21 .Or. lb(i2)==23) Then
    e(i1) = 0.
    i = i2
  Else
    e(i2) = 0.
    i = i1
  End If
  If (lb(i)==23) Then
    lb(i) = 30
  Else If (lb(i)==21) Then
    lb(i) = -30
  End If
  p(1, i) = p(1, i1) + p(1, i2)
  p(2, i) = p(2, i1) + p(2, i2)
  p(3, i) = p(3, i1) + p(3, i2)
! 2. DETERMINE THE MASS OF K* BY USING THE REACTION KINEMATICS

!lin-9/2012: check argument in sqrt():
  scheck = (e10+e20)**2 - p1**2 - p2**2 - p3**2
  If (scheck<0) Then
    Write (99, *) 'scheck49: ', scheck
    Write (99, *) 'scheck49', scheck, e10, e20, p(1, i), p(2, i), p(3, i)
    Write (99, *) 'scheck49-1', e(i1), p(1, i1), p(2, i1), p(3, i1)
    Write (99, *) 'scheck49-2', e(i2), p(1, i2), p(2, i2), p(3, i2)
  End If
  dm = sqrt(sngl(scheck))
!        DM=SQRT((E10+E20)**2-P(1,I)**2-P(2,I)**2-P(3,I)**2)

  e(i) = dm
  Return
End Subroutine ksreso

!--------------------------------------------------------
!************************************
!                                                                         *
Subroutine pertur(px, py, pz, srt, irun, i1, i2, nt, kp, icont)
!                                                                         *
!       PURPOSE:   TO PRODUCE CASCADE AND OMEGA PERTURBATIVELY            *
! sp 01/03/01
!                   40 cascade-
!                  -40 cascade-(bar)
!                   41 cascade0
!                  -41 cascade0(bar)
!                   45 Omega baryon
!                  -45 Omega baryon(bar)
!                   44 Di-Omega
!*********************************
  Parameter (maxstr=150001, maxr=1, pi=3.1415926)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Parameter (amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974, aks=0.895)
  Parameter (acas=1.3213, aome=1.6724, amrho=0.769, amomga=0.782)
  Parameter (aeta=0.548, adiomg=3.2288)
  Parameter (maxx=20, maxz=24)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /hh/proper(maxstr)
!c      SAVE /HH/
  Common /ff/f(-mx:mx, -my:my, -mz:mz, -mpx:mpx, -mpy:mpy, -mpz:mpzp)
!c      SAVE /ff/
  Common /gg/dx, dy, dz, dpx, dpy, dpz
!c      SAVE /gg/
  Common /input/nstar, ndirct, dir
!c      SAVE /INPUT/
  Common /nn/nnn
!c      SAVE /NN/
  Common /pa/rpion(3, maxstr, maxr)
!c      SAVE /PA/
  Common /pb/ppion(3, maxstr, maxr)
!c      SAVE /PB/
  Common /pc/epion(maxstr, maxr)
!c      SAVE /PC/
  Common /pd/lpion(maxstr, maxr)
!c      SAVE /PD/
  Common /pe/propi(maxstr, maxr)
!c      SAVE /PE/
  Common /rr/massr(0:maxr)
!c      SAVE /RR/
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
!     perturbative method is disabled:
!      common /imulst/ iperts
!
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save

  px0 = px
  py0 = py
  pz0 = pz
  lb1 = lb(i1)
  em1 = e(i1)
  x1 = r(1, i1)
  y1 = r(2, i1)
  z1 = r(3, i1)
  prob1 = proper(i1)
!
  lb2 = lb(i2)
  em2 = e(i2)
  x2 = r(1, i2)
  y2 = r(2, i2)
  z2 = r(3, i2)
  prob2 = proper(i2)
!
!                 !! flag for real 2-body process (1/0=no/yes)
  icont = 1
!                !! flag for elastic scatt only (-1=no)
  icsbel = -1

! K-/K*0bar + La/Si --> cascade + pi
! K+/K*0 + La/Si (bar) --> cascade-bar + pi
  If ((lb1==21 .Or. lb1==23 .Or. iabs(lb1)==30) .And. (iabs(lb2)>=14 .And. iabs(lb2)<=17)) Goto 60
  If ((lb2==21 .Or. lb2==23 .Or. iabs(lb2)==30) .And. (iabs(lb1)>=14 .And. iabs(lb1)<=17)) Goto 60
! K-/K*0bar + cascade --> omega + pi
! K+/K*0 + cascade-bar --> omega-bar + pi
  If ((lb1==21 .Or. lb1==23 .Or. iabs(lb1)==30) .And. (iabs(lb2)==40 .Or. iabs(lb2)==41)) Goto 70
  If ((lb2==21 .Or. lb2==23 .Or. iabs(lb2)==30) .And. (iabs(lb1)==40 .Or. iabs(lb1)==41)) Goto 70
!
! annhilation of cascade,cascade-bar, omega,omega-bar
!
! K- + La/Si <-- cascade + pi(eta,rho,omega)
! K+ + La/Si(bar) <-- cascade-bar + pi(eta,rho,omega)
  If ((((lb1>=3 .And. lb1<=5) .Or. lb1==0) .And. (iabs(lb2)==40 .Or. iabs(lb2)==41)) .Or. (((lb2>=3 .And. lb2<=5) .Or. lb2==0) .And. (iabs(lb1)==40 .Or. iabs(lb1)==41))) Goto 90
! K- + cascade <-- omega + pi
! K+ + cascade-bar <-- omega-bar + pi
!         if( (lb1.eq.0.and.iabs(lb2).eq.45)
!    &    .OR. (lb2.eq.0.and.iabs(lb1).eq.45) ) go to 110
  If (((lb1>=3 .And. lb1<=5) .And. iabs(lb2)==45) .Or. ((lb2>=3 .And. lb2<=5) .And. iabs(lb1)==45)) Goto 110
!

!----------------------------------------------------
!  for process:  K-bar + L(S) --> Ca + pi
!
  60 If (iabs(lb1)>=14 .And. iabs(lb1)<=17) Then
    asap = e(i1)
    akap = e(i2)
    idp = i1
  Else
    asap = e(i2)
    akap = e(i1)
    idp = i2
  End If
  app = 0.138
  If (srt<(acas+app)) Return
  srrt = srt - (acas+app) + (amn+akap)
  pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2-akap**2)
  sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
!lin pii & pff should be each divided by (4*srt**2),
!     but these two factors cancel out in the ratio pii/pff:
  pii = sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))
  pff = sqrt((srt**2-(asap+app)**2)*(srt**2-(asap-app)**2))
  cmat = sigca*pii/pff
  sigpi = cmat*sqrt((srt**2-(acas+app)**2)*(srt**2-(acas-app)**2))/sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
!
  sigeta = 0.
  If (srt>(acas+aeta)) Then
    srrt = srt - (acas+aeta) + (amn+akap)
    pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2-akap**2)
    sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
    cmat = sigca*pii/pff
    sigeta = cmat*sqrt((srt**2-(acas+aeta)**2)*(srt**2-(acas-aeta)**2))/sqrt((srt**2-(asap+akap)**2)*(srt**2-(asap-akap)**2))
  End If
!
  sigca = sigpi + sigeta
  sigpe = 0.
!lin-2/25/03 disable the perturb option:
!        if(iperts .eq. 1) sigpe = 40.   !! perturbative xsecn
  sig = amax1(sigpe, sigca)
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  brpp = sigca/sig
!
! else particle production
  If ((lb1>=14 .And. lb1<=17) .Or. (lb2>=14 .And. lb2<=17)) Then
!   !! cascade- or cascde0
    lbpp1 = 40 + int(2*ranart(nseed))
  Else
! elseif(lb1 .eq. -14 .or. lb2 .eq. -14)
!     !! cascade-bar- or cascde0 -bar
    lbpp1 = -40 - int(2*ranart(nseed))
  End If
  empp1 = acas
  If (ranart(nseed)<sigpi/sigca) Then
!    !! pion
    lbpp2 = 3 + int(3*ranart(nseed))
    empp2 = 0.138
  Else
!    !! eta
    lbpp2 = 0
    empp2 = aeta
  End If
!* check real process of cascade(bar) and pion formation
  If (ranart(nseed)<brpp) Then
!       !! real process flag
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
!  !! cascade formed with prob Gam
    proper(i1) = brpp
    lb(i2) = lbpp2
    e(i2) = empp2
!         !! pion/eta formed with prob 1.
    proper(i2) = 1.
  End If
! else only cascade(bar) formed perturbatively
  Goto 700

!----------------------------------------------------
!  for process:  Cas(bar) + K_bar(K) --> Om(bar) + pi  !! eta
!
  70 If (iabs(lb1)==40 .Or. iabs(lb1)==41) Then
    acap = e(i1)
    akap = e(i2)
    idp = i1
  Else
    acap = e(i2)
    akap = e(i1)
    idp = i2
  End If
  app = 0.138
!         ames = aeta
!  !! only pion
  ames = 0.138
  If (srt<(aome+ames)) Return
  srrt = srt - (aome+ames) + (amn+akap)
  pkaon = sqrt(((srrt**2-(amn**2+akap**2))/2./amn)**2-akap**2)
! use K(bar) + Ca --> Om + eta  xsecn same as  K(bar) + N --> Si + Pi
!  as Omega have no resonances
!** using same matrix elements as K-bar + N -> Si + pi
  sigomm = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
  cmat = sigomm*sqrt((srt**2-(amn+akap)**2)*(srt**2-(amn-akap)**2))/sqrt((srt**2-(asa+app)**2)*(srt**2-(asa-app)**2))
  sigom = cmat*sqrt((srt**2-(aome+ames)**2)*(srt**2-(aome-ames)**2))/sqrt((srt**2-(acap+akap)**2)*(srt**2-(acap-akap)**2))
  sigpe = 0.
!lin-2/25/03 disable the perturb option:
!         if(iperts .eq. 1) sigpe = 40.   !! perturbative xsecn
  sig = amax1(sigpe, sigom)
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
  If (ic==-1) Return
  brpp = sigom/sig
!
! else particle production
  If ((lb1>=40 .And. lb1<=41) .Or. (lb2>=40 .And. lb2<=41)) Then
!    !! omega
    lbpp1 = 45
  Else
! elseif(lb1 .eq. -40 .or. lb2 .eq. -40)
!    !! omega-bar
    lbpp1 = -45
  End If
  empp1 = aome
!           lbpp2 = 0    !! eta
!    !! pion
  lbpp2 = 3 + int(3*ranart(nseed))
  empp2 = ames
!
!* check real process of omega(bar) and pion formation
  xrand = ranart(nseed)
  If (xrand<(proper(idp)*brpp)) Then
!       !! real process flag
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
!  !! P_Om = P_Cas*Gam
    proper(i1) = proper(idp)*brpp
    lb(i2) = lbpp2
    e(i2) = empp2
!   !! pion formed with prob 1.
    proper(i2) = 1.
  Else If (xrand<brpp) Then
! else omega(bar) formed perturbatively and cascade destroyed
    e(idp) = 0.
  End If
  Goto 700

!-----------------------------------------------------------
!  for process:  Ca + pi/eta --> K-bar + L(S)
!
  90 If (iabs(lb1)==40 .Or. iabs(lb1)==41) Then
    acap = e(i1)
    app = e(i2)
    idp = i1
    idn = i2
  Else
    acap = e(i2)
    app = e(i1)
    idp = i2
    idn = i1
  End If
!            akal = (aka+aks)/2.  !! average of K and K* taken
!  !! using K only
  akal = aka
!
  alas = ala
  If (srt<=(alas+aka)) Return
  srrt = srt - (acap+app) + (amn+aka)
  pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2-aka**2)
!** using same matrix elements as K-bar + N -> La/Si + pi
  sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
  cmat = sigca*sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
  sigca = cmat*sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
!    !! pi
  dfr = 1./3.
!       !! eta
  If (lb(idn)==0) dfr = 1.
  sigcal = sigca*dfr*(srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/(srt**2-(acap-app)**2)
!
  alas = asa
  If (srt<=(alas+aka)) Then
    sigcas = 0.
  Else
    srrt = srt - (acap+app) + (amn+aka)
    pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2-aka**2)
! use K(bar) + La/Si --> Ca + Pi  xsecn same as  K(bar) + N --> Si + Pi
!** using same matrix elements as K-bar + N -> La/Si + pi
    sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
    cmat = sigca*sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/sqrt((srt**2-(alas+0.138)**2)*(srt**2-(alas-0.138)**2))
    sigca = cmat*sqrt((srt**2-(acap+app)**2)*(srt**2-(acap-app)**2))/sqrt((srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2))
!    !! pi
    dfr = 1.
!    !! eta
    If (lb(idn)==0) dfr = 3.
    sigcas = sigca*dfr*(srt**2-(alas+aka)**2)*(srt**2-(alas-aka)**2)/(srt**2-(acap+app)**2)/(srt**2-(acap-app)**2)
  End If
!
  sig = sigcal + sigcas
  brpp = 1.
  ds = sqrt(sig/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
!
!lin-2/25/03: checking elastic scatt after failure of inelastic scatt gives
!     conditional probability (in general incorrect), tell Pal to correct:
  If (ic==-1) Then
! check for elastic scatt, no particle annhilation
!  !! elastic cross section of 20 mb
    ds = sqrt(20.0/31.4)
    dsr = ds + 0.1
    Call distce(i1, i2, dsr, ds, dt, ec, srt, icsbel, px, py, pz)
    If (icsbel==-1) Return
    empp1 = em1
    empp2 = em2
    Goto 700
  End If
!
! else pert. produced cascade(bar) is annhilated OR real process
!
! DECIDE LAMBDA OR SIGMA PRODUCTION
!
  If (sigcal/sig>ranart(nseed)) Then
    If (lb1==40 .Or. lb1==41 .Or. lb2==40 .Or. lb2==41) Then
      lbpp1 = 21
      lbpp2 = 14
    Else
      lbpp1 = 23
      lbpp2 = -14
    End If
    alas = ala
  Else
    If (lb1==40 .Or. lb1==41 .Or. lb2==40 .Or. lb2==41) Then
      lbpp1 = 21
      lbpp2 = 15 + int(3*ranart(nseed))
    Else
      lbpp1 = 23
      lbpp2 = -15 - int(3*ranart(nseed))
    End If
    alas = asa
  End If
  empp1 = aka
  empp2 = alas
!
! check for real process for L/S(bar) and K(bar) formation
  If (ranart(nseed)<proper(idp)) Then
! real process
!       !! real process flag
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
!   !! K(bar) formed with prob 1.
    proper(i1) = 1.
    lb(i2) = lbpp2
    e(i2) = empp2
!   !! L/S(bar) formed with prob 1.
    proper(i2) = 1.
    Goto 700
  Else
! else only cascade(bar) annhilation & go out
    e(idp) = 0.
  End If
  Return
!
!----------------------------------------------------
!  for process:  Om(bar) + pi --> Cas(bar) + K_bar(K)
!
  110 If (lb1==45 .Or. lb1==-45) Then
    aomp = e(i1)
    app = e(i2)
    idp = i1
    idn = i2
  Else
    aomp = e(i2)
    app = e(i1)
    idp = i2
    idn = i1
  End If
!            akal = (aka+aks)/2.  !! average of K and K* taken
!  !! using K only
  akal = aka
  If (srt<=(acas+aka)) Return
  srrt = srt - (aome+app) + (amn+aka)
  pkaon = sqrt(((srrt**2-(amn**2+aka**2))/2./amn)**2-aka**2)
! use K(bar) + Ca --> Om + eta  xsecn same as  K(bar) + N --> Si + Pi
!** using same matrix elements as K-bar + N -> La/Si + pi
  sigca = 1.5*(aknpsg(pkaon)+aknpsg(pkaon))
  cmat = sigca*sqrt((srt**2-(amn+aka)**2)*(srt**2-(amn-aka)**2))/sqrt((srt**2-(asa+0.138)**2)*(srt**2-(asa-0.138)**2))
  sigom = cmat*sqrt((srt**2-(aomp+app)**2)*(srt**2-(aomp-app)**2))/sqrt((srt**2-(acas+aka)**2)*(srt**2-(acas-aka)**2))
!            dfr = 2.    !! eta
!    !! pion
  dfr = 2./3.
  sigom = sigom*dfr*(srt**2-(acas+aka)**2)*(srt**2-(acas-aka)**2)/(srt**2-(aomp+app)**2)/(srt**2-(aomp-app)**2)
!
  brpp = 1.
  ds = sqrt(sigom/31.4)
  dsr = ds + 0.1
  ec = (em1+em2+0.02)**2
  Call distce(i1, i2, dsr, ds, dt, ec, srt, ic, px, py, pz)
!
!lin-2/25/03: checking elastic scatt after failure of inelastic scatt gives
!     conditional probability (in general incorrect), tell Pal to correct:
  If (ic==-1) Then
! check for elastic scatt, no particle annhilation
!  !! elastic cross section of 20 mb
    ds = sqrt(20.0/31.4)
    dsr = ds + 0.1
    Call distce(i1, i2, dsr, ds, dt, ec, srt, icsbel, px, py, pz)
    If (icsbel==-1) Return
    empp1 = em1
    empp2 = em2
    Goto 700
  End If
!
! else pert. produced omega(bar) annhilated  OR real process
! annhilate only pert. omega, rest from hijing go out WITHOUT annhil.
  If (lb1==45 .Or. lb2==45) Then
!  !! Ca
    lbpp1 = 40 + int(2*ranart(nseed))
!   !! K-
    lbpp2 = 21
  Else
! elseif(lb1 .eq. -45 .or. lb2 .eq. -45)
!    !! Ca-bar
    lbpp1 = -40 - int(2*ranart(nseed))
!      !! K+
    lbpp2 = 23
  End If
  empp1 = acas
  empp2 = aka
!
! check for real process for Cas(bar) and K(bar) formation
  If (ranart(nseed)<proper(idp)) Then
!       !! real process flag
    icont = 0
    lb(i1) = lbpp1
    e(i1) = empp1
!   !! P_Cas(bar) = P_Om(bar)
    proper(i1) = proper(idp)
    lb(i2) = lbpp2
    e(i2) = empp2
!   !! K(bar) formed with prob 1.
    proper(i2) = 1.
!
  Else
! else Cascade(bar)  produced and Omega(bar) annhilated
    e(idp) = 0.
  End If
!   !! for produced particles
  Goto 700
!
!-----------------------------------------------------------
  700 Continue
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-empp1**2-empp2**2)**2 - 4.0*(empp1*empp2)**2
  If (pr2<=0.) pr2 = 0.00000001
  pr = sqrt(pr2)/(2.*srt)
! using isotropic
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
! ROTATE IT
  Call rotate(px0, py0, pz0, px, py, pz)
  If (icont==0) Return
!
! LORENTZ-TRANSFORMATION INTO CMS FRAME
  e1cm = sqrt(empp1**2+px**2+py**2+pz**2)
  p1beta = px*betax + py*betay + pz*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)+e1cm)
  ppt11 = betax*transf + px
  ppt12 = betay*transf + py
  ppt13 = betaz*transf + pz
!
!c** for elastic scattering update the momentum of pertb particles
  If (icsbel/=-1) Then
!            if(EMpp1 .gt. 0.9)then
    p(1, i1) = ppt11
    p(2, i1) = ppt12
    p(3, i1) = ppt13
!            else
    e2cm = sqrt(empp2**2+px**2+py**2+pz**2)
    transf = gamma*(-gamma*p1beta/(gamma+1)+e2cm)
    ppt21 = betax*transf - px
    ppt22 = betay*transf - py
    ppt23 = betaz*transf - pz
    p(1, i2) = ppt21
    p(2, i2) = ppt22
    p(3, i2) = ppt23
!            endif
    Return
  End If
!lin-5/2008:
!2008        X01 = 1.0 - 2.0 * RANART(NSEED)
!            Y01 = 1.0 - 2.0 * RANART(NSEED)
!            Z01 = 1.0 - 2.0 * RANART(NSEED)
!        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2008
!                Xpt=X1+0.5*x01
!                Ypt=Y1+0.5*y01
!                Zpt=Z1+0.5*z01
  xpt = x1
  ypt = y1
  zpt = z1
!
!
!          if(lbpp1 .eq. 45)then
!           write(*,*)'II lb1,lb2,lbpp1,empp1,proper(idp),brpp'
!           write(*,*)lb1,lb2,lbpp1,empp1,proper(idp),brpp
!          endif
!
  nnn = nnn + 1
  propi(nnn, irun) = proper(idp)*brpp
  lpion(nnn, irun) = lbpp1
  epion(nnn, irun) = empp1
  rpion(1, nnn, irun) = xpt
  rpion(2, nnn, irun) = ypt
  rpion(3, nnn, irun) = zpt
  ppion(1, nnn, irun) = ppt11
  ppion(2, nnn, irun) = ppt12
  ppion(3, nnn, irun) = ppt13
!lin-5/2008:
  dppion(nnn, irun) = dpertp(i1)*dpertp(i2)
  Return
End Subroutine pertur
!*********************************
!  sp 12/08/00                                                         *
Subroutine crhb(px, py, pz, srt, i1, i2, iblock)
!     PURPOSE:                                                         *
!        DEALING WITH hyperon+N(D,N*)->hyp+N(D,N*) elastic PROCESS     *
!     NOTE   :                                                         *
!
!     QUANTITIES:                                                 *
!           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
!           SRT      - SQRT OF S                                       *
!           IBLOCK   - THE INFORMATION BACK                            *
!                     144-> hyp+N(D,N*)->hyp+N(D,N*)
!*********************************
  Parameter (maxstr=150001, maxr=1, amn=0.939457, amp=0.93828, ap1=0.13496, ap2=0.13957, am0=1.232, pi=3.1415926, cutoff=1.8966, avmass=0.9383)
  Parameter (aka=0.498, ala=1.1157, asa=1.1974)
  Parameter (mx=4, my=4, mz=8, mpx=4, mpy=4, mpz=10, mpzp=10)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Common /input1/masspr, massta, iseed, iavoid, dt
!c      SAVE /input1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  px0 = px
  py0 = py
  pz0 = pz
!-----------------------------------------------------------------------
  iblock = 144
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
!-----------------------------------------------------------------------
! CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
! ENERGY CONSERVATION
  pr2 = (srt**2-em1**2-em2**2)**2 - 4.0*(em1*em2)**2
  If (pr2<=0.) pr2 = 1.E-09
  pr = sqrt(pr2)/(2.*srt)
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
  pz = pr*c1
  px = pr*s1*ct1
  py = pr*s1*st1
  Return
End Subroutine crhb
!***************************************
! sp 04/05/01
! Purpose: lambda-baryon elastic xsection as a functon of their cms energy
Subroutine lambar(i1, i2, srt, siglab)
!  srt    = DSQRT(s) in GeV                                               *
!  siglab = lambda-nuclar elastic cross section in mb
!         = 12 + 0.43/p_lab**3.3 (mb)
!
! (2) Calculate p(lab) from srt [GeV], since the formular in the
! reference applies only to the case of a p_bar on a proton at rest
! Formula used: srt**2=2.*pmass*(pmass+sqrt(pmass**2+plab**2))
!****************************
  Parameter (maxstr=150001)
  Common /aa/r(3, maxstr)
!c      SAVE /AA/
  Common /bb/p(3, maxstr)
!c      SAVE /BB/
  Common /cc/e(maxstr)
!c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
!c      SAVE /EE/
  Save

  siglab = 1.E-06
  If (iabs(lb(i1))>=14 .And. iabs(lb(i1))<=17) Then
    eml = e(i1)
    emb = e(i2)
  Else
    eml = e(i2)
    emb = e(i1)
  End If
  pthr = srt**2 - eml**2 - emb**2
  If (pthr>0.) Then
    plab2 = (pthr/2./emb)**2 - eml**2
    If (plab2>0) Then
      plab = sqrt(plab2)
      siglab = 12. + 0.43/(plab**3.3)
      If (siglab>200.) siglab = 200.
    End If
  End If
  Return
End Subroutine lambar
!------------------------------------------------------------------
!lin-7/26/03 improve speed
!**************************************
Subroutine distc0(drmax, deltr0, dt, ifirst, px1cm, py1cm, pz1cm, x1, y1, z1, px1, py1, pz1, em1, x2, y2, z2, px2, py2, pz2, em2)
! PURPOSE : CHECK IF THE COLLISION BETWEEN TWO PARTICLES CAN HAPPEN
!           BY CHECKING
!                      (2) IF PARTICLE WILL PASS EACH OTHER WITHIN
!           TWO HARD CORE RADIUS.
!                      (3) IF PARTICLES WILL GET CLOSER.
! VARIABLES :
!           Ifirst=1 COLLISION may HAPPENED
!           Ifirst=-1 COLLISION CAN NOT HAPPEN
!****************************************
  Common /bg/betax, betay, betaz, gamma
!c      SAVE /BG/
  Save
  ifirst = -1
  e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
!NOW PARTICLES ARE CLOSE ENOUGH TO EACH OTHER !
  e2 = sqrt(em2**2+px2**2+py2**2+pz2**2)
!NOW THERE IS ENOUGH ENERGY AVAILABLE !
!LORENTZ-TRANSFORMATION IN I1-I2-C.M. SYSTEM
! BETAX, BETAY, BETAZ AND GAMMA HAVE BEEN GIVEN IN THE SUBROUTINE CMS
!TRANSFORMATION OF MOMENTA (PX1CM = - PX2CM)
  p1beta = px1*betax + py1*betay + pz1*betaz
  transf = gamma*(gamma*p1beta/(gamma+1)-e1)
  prcm = sqrt(px1cm**2+py1cm**2+pz1cm**2)
  If (prcm<=0.00001) Return
!TRANSFORMATION OF SPATIAL DISTANCE
  drbeta = betax*(x1-x2) + betay*(y1-y2) + betaz*(z1-z2)
  transf = gamma*gamma*drbeta/(gamma+1)
  dxcm = betax*transf + x1 - x2
  dycm = betay*transf + y1 - y2
  dzcm = betaz*transf + z1 - z2
!DETERMINING IF THIS IS THE POINT OF CLOSEST APPROACH
  drcm = sqrt(dxcm**2+dycm**2+dzcm**2)
  dzz = (px1cm*dxcm+py1cm*dycm+pz1cm*dzcm)/prcm
  If ((drcm**2-dzz**2)<=0.) Then
    bbb = 0.
  Else
    bbb = sqrt(drcm**2-dzz**2)
  End If
!WILL PARTICLE PASS EACH OTHER WITHIN 2 * HARD CORE RADIUS ?
  If (bbb>drmax) Return
  relvel = prcm*(1.0/e1+1.0/e2)
  ddd = relvel*dt*0.5
!WILL PARTICLES GET CLOSER ?
  If (abs(ddd)<abs(dzz)) Return
  ifirst = 1
  Return
End Subroutine distc0
!---------------------------------------------------------------------------
!
!lin-8/2008 B+B->Deuteron+Meson cross section in mb:
Subroutine sbbdm(srt, sdprod, ianti, lbm, xmm, pfinal)
  Parameter (xmd=1.8756, ap1=0.13496, ap2=0.13957, xmrho=0.770, xmomega=0.782, xmeta=0.548, srt0=2.012)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Common /rndf77/nseed
  Save
!
  sdprod = 0.
  sbbdpi = 0.
  sbbdrho = 0.
  sbbdomega = 0.
  sbbdeta = 0.
  If (srt<=(em1+em2)) Return
!
  ilb1 = iabs(lb1)
  ilb2 = iabs(lb2)
!test off check Xsec using fixed mass for resonances:
!      if(ilb1.ge.6.and.ilb1.le.9) then
!         em1=1.232
!      elseif(ilb1.ge.10.and.ilb1.le.11) then
!         em1=1.44
!      elseif(ilb1.ge.12.and.ilb1.le.13) then
!         em1=1.535
!      endif
!      if(ilb2.ge.6.and.ilb2.le.9) then
!         em2=1.232
!      elseif(ilb2.ge.10.and.ilb2.le.11) then
!         em2=1.44
!      elseif(ilb2.ge.12.and.ilb2.le.13) then
!         em2=1.535
!      endif
!
  s = srt**2

!lin-9/2012: check argument in sqrt():
  scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
  If (scheck<=0) Then
    Write (99, *) 'scheck50: ', scheck
    Stop
  End If
  pinitial = sqrt(scheck)/2./srt
!      pinitial=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt

  fs = fnndpi(s)
!     Determine isospin and spin factors for the ratio between
!     BB->Deuteron+Meson and Deuteron+Meson->BB cross sections:
  If (idxsec==1 .Or. idxsec==2) Then
!     Assume B+B -> d+Meson has the same cross sections as N+N -> d+pi:
  Else
!     Assume d+Meson -> B+B has the same cross sections as d+pi -> N+N,
!     then determine B+B -> d+Meson cross sections:
    If (ilb1>=1 .And. ilb1<=2 .And. ilb2>=1 .And. ilb2<=2) Then
      pifactor = 9./8.
    Else If ((ilb1>=1 .And. ilb1<=2 .And. ilb2>=6 .And. ilb2<=9) .Or. (ilb2>=1 .And. ilb2<=2 .And. ilb1>=6 .And. ilb1<=9)) Then
      pifactor = 9./64.
    Else If ((ilb1>=1 .And. ilb1<=2 .And. ilb2>=10 .And. ilb2<=13) .Or. (ilb2>=1 .And. ilb2<=2 .And. ilb1>=10 .And. ilb1<=13)) Then
      pifactor = 9./16.
    Else If (ilb1>=6 .And. ilb1<=9 .And. ilb2>=6 .And. ilb2<=9) Then
      pifactor = 9./128.
    Else If ((ilb1>=6 .And. ilb1<=9 .And. ilb2>=10 .And. ilb2<=13) .Or. (ilb2>=6 .And. ilb2<=9 .And. ilb1>=10 .And. ilb1<=13)) Then
      pifactor = 9./64.
    Else If ((ilb1>=10 .And. ilb1<=11 .And. ilb2>=10 .And. ilb2<=11) .Or. (ilb2>=12 .And. ilb2<=13 .And. ilb1>=12 .And. ilb1<=13)) Then
      pifactor = 9./8.
    Else If ((ilb1>=10 .And. ilb1<=11 .And. ilb2>=12 .And. ilb2<=13) .Or. (ilb2>=10 .And. ilb2<=11 .And. ilb1>=12 .And. ilb1<=13)) Then
      pifactor = 9./16.
    End If
  End If
!     d pi: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
!     (1) FOR P+P->Deuteron+pi+:
  If ((ilb1*ilb2)==1) Then
    lbm = 5
    If (ianti==1) lbm = 3
    xmm = ap2
!     (2)FOR N+N->Deuteron+pi-:
  Else If (ilb1==2 .And. ilb2==2) Then
    lbm = 3
    If (ianti==1) lbm = 5
    xmm = ap2
!     (3)FOR N+P->Deuteron+pi0:
  Else If ((ilb1*ilb2)==2) Then
    lbm = 4
    xmm = ap1
  Else
!     For baryon resonances, use isospin-averaged cross sections:
    lbm = 3 + int(3*ranart(nseed))
    If (lbm==4) Then
      xmm = ap1
    Else
      xmm = ap2
    End If
  End If
!
  If (srt>=(xmd+xmm)) Then
    pfinal = sqrt((s-(xmd+xmm)**2)*(s-(xmd-xmm)**2))/2./srt
    If ((ilb1==1 .And. ilb2==1) .Or. (ilb1==2 .And. ilb2==2)) Then
!     for pp or nn initial states:
      sbbdpi = fs*pfinal/pinitial/4.
    Else If ((ilb1==1 .And. ilb2==2) .Or. (ilb1==2 .And. ilb2==1)) Then
!     factor of 1/2 for pn or np initial states:
      sbbdpi = fs*pfinal/pinitial/4./2.
    Else
!     for other BB initial states (spin- and isospin averaged):
      If (idxsec==1) Then
!     1: assume the same |matrix element|**2/s (after averaging over initial
!     spins and isospins) for B+B -> deuteron+meson at the same sqrt(s);
        sbbdpi = fs*pfinal/pinitial*3./16.
      Else If (idxsec==2 .Or. idxsec==4) Then
        threshold = amax1(xmd+xmm, em1+em2)
        snew = (srt-threshold+srt0)**2
        If (idxsec==2) Then
!     2: assume the same |matrix element|**2/s for B+B -> deuteron+meson
!     at the same sqrt(s)-threshold:
          sbbdpi = fnndpi(snew)*pfinal/pinitial*3./16.
        Else If (idxsec==4) Then
!     4: assume the same |matrix element|**2/s for B+B <- deuteron+meson
!     at the same sqrt(s)-threshold:
          sbbdpi = fnndpi(snew)*pfinal/pinitial/6.*pifactor
        End If
      Else If (idxsec==3) Then
!     3: assume the same |matrix element|**2/s for B+B <- deuteron+meson
!     at the same sqrt(s):
        sbbdpi = fs*pfinal/pinitial/6.*pifactor
      End If
!
    End If
  End If
!
!     d rho: DETERMINE THE CROSS SECTION TO THIS FINAL STATE:
  If (srt>(xmd+xmrho)) Then
    pfinal = sqrt((s-(xmd+xmrho)**2)*(s-(xmd-xmrho)**2))/2./srt
    If (idxsec==1) Then
      sbbdrho = fs*pfinal/pinitial*3./16.
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmd+xmrho, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sbbdrho = fnndpi(snew)*pfinal/pinitial*3./16.
      Else If (idxsec==4) Then
!     The spin- and isospin-averaged factor is 3-times larger for rho:
        sbbdrho = fnndpi(snew)*pfinal/pinitial/6.*(pifactor*3.)
      End If
    Else If (idxsec==3) Then
      sbbdrho = fs*pfinal/pinitial/6.*(pifactor*3.)
    End If
  End If
!
!     d omega: DETERMINE THE CROSS SECTION TO THIS FINAL STATE:
  If (srt>(xmd+xmomega)) Then
    pfinal = sqrt((s-(xmd+xmomega)**2)*(s-(xmd-xmomega)**2))/2./srt
    If (idxsec==1) Then
      sbbdomega = fs*pfinal/pinitial*3./16.
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmd+xmomega, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sbbdomega = fnndpi(snew)*pfinal/pinitial*3./16.
      Else If (idxsec==4) Then
        sbbdomega = fnndpi(snew)*pfinal/pinitial/6.*pifactor
      End If
    Else If (idxsec==3) Then
      sbbdomega = fs*pfinal/pinitial/6.*pifactor
    End If
  End If
!
!     d eta: DETERMINE THE CROSS SECTION TO THIS FINAL STATE:
  If (srt>(xmd+xmeta)) Then
    pfinal = sqrt((s-(xmd+xmeta)**2)*(s-(xmd-xmeta)**2))/2./srt
    If (idxsec==1) Then
      sbbdeta = fs*pfinal/pinitial*3./16.
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmd+xmeta, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sbbdeta = fnndpi(snew)*pfinal/pinitial*3./16.
      Else If (idxsec==4) Then
        sbbdeta = fnndpi(snew)*pfinal/pinitial/6.*(pifactor/3.)
      End If
    Else If (idxsec==3) Then
      sbbdeta = fs*pfinal/pinitial/6.*(pifactor/3.)
    End If
  End If
!
  sdprod = sbbdpi + sbbdrho + sbbdomega + sbbdeta
!test off
!      write(99,111) srt,sbbdpi,sbbdrho,sbbdomega,sbbdeta,sdprod
! 111  format(6(f8.2,1x))
!
  If (sdprod<=0) Return
!
!     choose final state and assign masses here:
  x1 = ranart(nseed)
  If (x1<=sbbdpi/sdprod) Then
!     use the above-determined lbm and xmm.
  Else If (x1<=(sbbdpi+sbbdrho)/sdprod) Then
    lbm = 25 + int(3*ranart(nseed))
    xmm = xmrho
  Else If (x1<=(sbbdpi+sbbdrho+sbbdomega)/sdprod) Then
    lbm = 28
    xmm = xmomega
  Else
    lbm = 0
    xmm = xmeta
  End If
!
  Return
End Subroutine sbbdm
!
!     Generate angular distribution of Deuteron in the CMS frame:
Subroutine bbdangle(pxd, pyd, pzd, nt, ipert1, ianti, idloop, pfinal, dprob1, lbm)
  Parameter (pi=3.1415926)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /rndf77/nseed
  Common /para8/idpert, npertd, idxsec
  Common /arevt/iaevt, iarun, miss
  Save
!     take isotropic distribution for now:
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pzd = pfinal*c1
  pxd = pfinal*s1*ct1
  pyd = pfinal*s1*st1
!lin-5/2008 track the number of produced deuterons:
  If (idpert==1 .And. npertd>=1) Then
    dprob = dprob1
  Else If (idpert==2 .And. npertd>=1) Then
    dprob = 1./float(npertd)
  End If
  If (ianti==0) Then
    If (idpert==0 .Or. (idpert==1 .And. ipert1==0) .Or. (idpert==2 .And. idloop==(npertd+1))) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (regular d prodn)    @evt#', iaevt, ' @nt=', nt
    Else If ((idpert==1 .Or. idpert==2) .And. idloop==npertd) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (pert d prodn)       @evt#', iaevt, ' @nt=', nt, ' @prob=', dprob
    End If
  Else
    If (idpert==0 .Or. (idpert==1 .And. ipert1==0) .Or. (idpert==2 .And. idloop==(npertd+1))) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (regular dbar prodn) @evt#', iaevt, ' @nt=', nt
    Else If ((idpert==1 .Or. idpert==2) .And. idloop==npertd) Then
      Write (91, *) lb1, ' *', lb2, ' ->d+', lbm, ' (pert dbar prodn)    @evt#', iaevt, ' @nt=', nt, ' @prob=', dprob
    End If
  End If
!
  Return
End Subroutine bbdangle
!
!     Deuteron+Meson->B+B cross section (in mb)
Subroutine sdmbb(srt, sdm, ianti)
  Parameter (amn=0.939457, amp=0.93828, am0=1.232, am1440=1.44, am1535=1.535, srt0=2.012)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /dpifsl/lbnn1, lbnn2, lbnd1, lbnd2, lbns1, lbns2, lbnp1, lbnp2, lbdd1, lbdd2, lbds1, lbds2, lbdp1, lbdp2, lbss1, lbss2, lbsp1, lbsp2, lbpp1, lbpp2
  Common /dpifsm/xmnn1, xmnn2, xmnd1, xmnd2, xmns1, xmns2, xmnp1, xmnp2, xmdd1, xmdd2, xmds1, xmds2, xmdp1, xmdp2, xmss1, xmss2, xmsp1, xmsp2, xmpp1, xmpp2
  Common /dpisig/sdmel, sdmnn, sdmnd, sdmns, sdmnp, sdmdd, sdmds, sdmdp, sdmss, sdmsp, sdmpp
  Common /para8/idpert, npertd, idxsec
  Common /rndf77/nseed
  Save
!
  sdm = 0.
  sdmel = 0.
  sdmnn = 0.
  sdmnd = 0.
  sdmns = 0.
  sdmnp = 0.
  sdmdd = 0.
  sdmds = 0.
  sdmdp = 0.
  sdmss = 0.
  sdmsp = 0.
  sdmpp = 0.
!test off check Xsec using fixed mass for resonances:
!      if(lb1.ge.25.and.lb1.le.27) then
!         em1=0.776
!      elseif(lb1.eq.28) then
!         em1=0.783
!      elseif(lb1.eq.0) then
!         em1=0.548
!      endif
!      if(lb2.ge.25.and.lb2.le.27) then
!         em2=0.776
!      elseif(lb2.eq.28) then
!         em2=0.783
!      elseif(lb2.eq.0) then
!         em2=0.548
!      endif
!
  If (srt<=(em1+em2)) Return
  s = srt**2
  pinitial = sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt
  fs = fnndpi(s)
!     Determine isospin and spin factors for the ratio between
!     Deuteron+Meson->BB and BB->Deuteron+Meson cross sections:
  If (idxsec==1 .Or. idxsec==2) Then
!     Assume B+B -> d+Meson has the same cross sections as N+N -> d+pi,
!     then determine d+Meson -> B+B cross sections:
    If ((lb1>=3 .And. lb1<=5) .Or. (lb2>=3 .And. lb2<=5)) Then
      xnnfactor = 8./9.
    Else If ((lb1>=25 .And. lb1<=27) .Or. (lb2>=25 .And. lb2<=27)) Then
      xnnfactor = 8./27.
    Else If (lb1==28 .Or. lb2==28) Then
      xnnfactor = 8./9.
    Else If (lb1==0 .Or. lb2==0) Then
      xnnfactor = 8./3.
    End If
  Else
!     Assume d+Meson -> B+B has the same cross sections as d+pi -> N+N:
  End If
!lin-9/2008 For elastic collisions:
  If (idxsec==1 .Or. idxsec==3) Then
!     1/3: assume the same |matrix element|**2/s (after averaging over initial
!     spins and isospins) for d+Meson elastic at the same sqrt(s);
    sdmel = fdpiel(s)
  Else If (idxsec==2 .Or. idxsec==4) Then
!     2/4: assume the same |matrix element|**2/s (after averaging over initial
!     spins and isospins) for d+Meson elastic at the same sqrt(s)-threshold:
    threshold = em1 + em2
    snew = (srt-threshold+srt0)**2
    sdmel = fdpiel(snew)
  End If
!
!     NN: DETERMINE THE CHARGE STATES OF PARTICLESIN THE FINAL STATE
  If (((lb1==5 .Or. lb2==5 .Or. lb1==27 .Or. lb2==27) .And. ianti==0) .Or. ((lb1==3 .Or. lb2==3 .Or. lb1==25 .Or. lb2==25) .And. ianti==1)) Then
!     (1) FOR Deuteron+(pi+,rho+) -> P+P or DeuteronBar+(pi-,rho-)-> PBar+PBar:
    lbnn1 = 1
    lbnn2 = 1
    xmnn1 = amp
    xmnn2 = amp
  Else If (lb1==3 .Or. lb2==3 .Or. lb1==26 .Or. lb2==26 .Or. lb1==28 .Or. lb2==28 .Or. lb1==0 .Or. lb2==0) Then
!     (2) FOR Deuteron+(pi0,rho0,omega,eta) -> N+P
!     or DeuteronBar+(pi0,rho0,omega,eta) ->NBar+PBar:
    lbnn1 = 2
    lbnn2 = 1
    xmnn1 = amn
    xmnn2 = amp
  Else
!     (3) FOR Deuteron+(pi-,rho-) -> N+N or DeuteronBar+(pi+,rho+)-> NBar+NBar:
    lbnn1 = 2
    lbnn2 = 2
    xmnn1 = amn
    xmnn2 = amn
  End If
  If (srt>(xmnn1+xmnn2)) Then
    pfinal = sqrt((s-(xmnn1+xmnn2)**2)*(s-(xmnn1-xmnn2)**2))/2./srt
    If (idxsec==1) Then
!     1: assume the same |matrix element|**2/s (after averaging over initial
!     spins and isospins) for B+B -> deuteron+meson at the same sqrt(s);
      sdmnn = fs*pfinal/pinitial*3./16.*xnnfactor
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmnn1+xmnn2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
!     2: assume the same |matrix element|**2/s for B+B -> deuteron+meson
!     at the same sqrt(s)-threshold:
        sdmnn = fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
      Else If (idxsec==4) Then
!     4: assume the same |matrix element|**2/s for B+B <- deuteron+meson
!     at the same sqrt(s)-threshold:
        sdmnn = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
!     3: assume the same |matrix element|**2/s for B+B <- deuteron+meson
!     at the same sqrt(s):
      sdmnn = fs*pfinal/pinitial/6.
    End If
  End If
!
!     ND: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbnd1 = 1 + int(2*ranart(nseed))
  lbnd2 = 6 + int(4*ranart(nseed))
  If (lbnd1==1) Then
    xmnd1 = amp
  Else If (lbnd1==2) Then
    xmnd1 = amn
  End If
  xmnd2 = am0
  If (srt>(xmnd1+xmnd2)) Then
    pfinal = sqrt((s-(xmnd1+xmnd2)**2)*(s-(xmnd1-xmnd2)**2))/2./srt
    If (idxsec==1) Then
!     The spin- and isospin-averaged factor is 8-times larger for ND:
      sdmnd = fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmnd1+xmnd2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmnd = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
      Else If (idxsec==4) Then
        sdmnd = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmnd = fs*pfinal/pinitial/6.
    End If
  End If
!
!     NS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbns1 = 1 + int(2*ranart(nseed))
  lbns2 = 10 + int(2*ranart(nseed))
  If (lbns1==1) Then
    xmns1 = amp
  Else If (lbns1==2) Then
    xmns1 = amn
  End If
  xmns2 = am1440
  If (srt>(xmns1+xmns2)) Then
    pfinal = sqrt((s-(xmns1+xmns2)**2)*(s-(xmns1-xmns2)**2))/2./srt
    If (idxsec==1) Then
      sdmns = fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmns1+xmns2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmns = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
      Else If (idxsec==4) Then
        sdmns = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmns = fs*pfinal/pinitial/6.
    End If
  End If
!
!     NP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbnp1 = 1 + int(2*ranart(nseed))
  lbnp2 = 12 + int(2*ranart(nseed))
  If (lbnp1==1) Then
    xmnp1 = amp
  Else If (lbnp1==2) Then
    xmnp1 = amn
  End If
  xmnp2 = am1535
  If (srt>(xmnp1+xmnp2)) Then
    pfinal = sqrt((s-(xmnp1+xmnp2)**2)*(s-(xmnp1-xmnp2)**2))/2./srt
    If (idxsec==1) Then
      sdmnp = fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmnp1+xmnp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmnp = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
      Else If (idxsec==4) Then
        sdmnp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmnp = fs*pfinal/pinitial/6.
    End If
  End If
!
!     DD: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbdd1 = 6 + int(4*ranart(nseed))
  lbdd2 = 6 + int(4*ranart(nseed))
  xmdd1 = am0
  xmdd2 = am0
  If (srt>(xmdd1+xmdd2)) Then
    pfinal = sqrt((s-(xmdd1+xmdd2)**2)*(s-(xmdd1-xmdd2)**2))/2./srt
    If (idxsec==1) Then
      sdmdd = fs*pfinal/pinitial*3./16.*(xnnfactor*16.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmdd1+xmdd2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmdd = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*16.)
      Else If (idxsec==4) Then
        sdmdd = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmdd = fs*pfinal/pinitial/6.
    End If
  End If
!
!     DS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbds1 = 6 + int(4*ranart(nseed))
  lbds2 = 10 + int(2*ranart(nseed))
  xmds1 = am0
  xmds2 = am1440
  If (srt>(xmds1+xmds2)) Then
    pfinal = sqrt((s-(xmds1+xmds2)**2)*(s-(xmds1-xmds2)**2))/2./srt
    If (idxsec==1) Then
      sdmds = fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmds1+xmds2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmds = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
      Else If (idxsec==4) Then
        sdmds = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmds = fs*pfinal/pinitial/6.
    End If
  End If
!
!     DP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbdp1 = 6 + int(4*ranart(nseed))
  lbdp2 = 12 + int(2*ranart(nseed))
  xmdp1 = am0
  xmdp2 = am1535
  If (srt>(xmdp1+xmdp2)) Then
    pfinal = sqrt((s-(xmdp1+xmdp2)**2)*(s-(xmdp1-xmdp2)**2))/2./srt
    If (idxsec==1) Then
      sdmdp = fs*pfinal/pinitial*3./16.*(xnnfactor*8.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmdp1+xmdp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmdp = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*8.)
      Else If (idxsec==4) Then
        sdmdp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmdp = fs*pfinal/pinitial/6.
    End If
  End If
!
!     SS: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbss1 = 10 + int(2*ranart(nseed))
  lbss2 = 10 + int(2*ranart(nseed))
  xmss1 = am1440
  xmss2 = am1440
  If (srt>(xmss1+xmss2)) Then
    pfinal = sqrt((s-(xmss1+xmss2)**2)*(s-(xmss1-xmss2)**2))/2./srt
    If (idxsec==1) Then
      sdmss = fs*pfinal/pinitial*3./16.*xnnfactor
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmss1+xmss2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmss = fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
      Else If (idxsec==4) Then
        sdmss = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmns = fs*pfinal/pinitial/6.
    End If
  End If
!
!     SP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbsp1 = 10 + int(2*ranart(nseed))
  lbsp2 = 12 + int(2*ranart(nseed))
  xmsp1 = am1440
  xmsp2 = am1535
  If (srt>(xmsp1+xmsp2)) Then
    pfinal = sqrt((s-(xmsp1+xmsp2)**2)*(s-(xmsp1-xmsp2)**2))/2./srt
    If (idxsec==1) Then
      sdmsp = fs*pfinal/pinitial*3./16.*(xnnfactor*2.)
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmsp1+xmsp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmsp = fnndpi(snew)*pfinal/pinitial*3./16.*(xnnfactor*2.)
      Else If (idxsec==4) Then
        sdmsp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmsp = fs*pfinal/pinitial/6.
    End If
  End If
!
!     PP: DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
  lbpp1 = 12 + int(2*ranart(nseed))
  lbpp2 = 12 + int(2*ranart(nseed))
  xmpp1 = am1535
  xmpp2 = am1535
  If (srt>(xmpp1+xmpp2)) Then
    pfinal = sqrt((s-(xmpp1+xmpp2)**2)*(s-(xmpp1-xmpp2)**2))/2./srt
    If (idxsec==1) Then
      sdmpp = fs*pfinal/pinitial*3./16.*xnnfactor
    Else If (idxsec==2 .Or. idxsec==4) Then
      threshold = amax1(xmpp1+xmpp2, em1+em2)
      snew = (srt-threshold+srt0)**2
      If (idxsec==2) Then
        sdmpp = fnndpi(snew)*pfinal/pinitial*3./16.*xnnfactor
      Else If (idxsec==4) Then
        sdmpp = fnndpi(snew)*pfinal/pinitial/6.
      End If
    Else If (idxsec==3) Then
      sdmpp = fs*pfinal/pinitial/6.
    End If
  End If
!
  sdm = sdmel + sdmnn + sdmnd + sdmns + sdmnp + sdmdd + sdmds + sdmdp + sdmss + sdmsp + sdmpp
  If (ianti==1) Then
    lbnn1 = -lbnn1
    lbnn2 = -lbnn2
    lbnd1 = -lbnd1
    lbnd2 = -lbnd2
    lbns1 = -lbns1
    lbns2 = -lbns2
    lbnp1 = -lbnp1
    lbnp2 = -lbnp2
    lbdd1 = -lbdd1
    lbdd2 = -lbdd2
    lbds1 = -lbds1
    lbds2 = -lbds2
    lbdp1 = -lbdp1
    lbdp2 = -lbdp2
    lbss1 = -lbss1
    lbss2 = -lbss2
    lbsp1 = -lbsp1
    lbsp2 = -lbsp2
    lbpp1 = -lbpp1
    lbpp2 = -lbpp2
  End If
!test off
!      write(98,100) srt,sdmnn,sdmnd,sdmns,sdmnp,sdmdd,sdmds,sdmdp,
!     1     sdmss,sdmsp,sdmpp,sdm
! 100  format(f5.2,11(1x,f5.1))
!
  Return
End Subroutine sdmbb
!
!lin-9/2008 Deuteron+Meson ->B+B and elastic collisions
Subroutine crdmbb(px, py, pz, srt, i1, i2, iblock, ntag, sig, nt, ianti)
  Parameter (maxstr=150001, maxr=1)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /bg/betax, betay, betaz, gamma
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /arevt/iaevt, iarun, miss
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Common /dpifsl/lbnn1, lbnn2, lbnd1, lbnd2, lbns1, lbns2, lbnp1, lbnp2, lbdd1, lbdd2, lbds1, lbds2, lbdp1, lbdp2, lbss1, lbss2, lbsp1, lbsp2, lbpp1, lbpp2
  Common /dpifsm/xmnn1, xmnn2, xmnd1, xmnd2, xmns1, xmns2, xmnp1, xmnp2, xmdd1, xmdd2, xmds1, xmds2, xmdp1, xmdp2, xmss1, xmss2, xmsp1, xmsp2, xmpp1, xmpp2
  Common /dpisig/sdmel, sdmnn, sdmnd, sdmns, sdmnp, sdmdd, sdmds, sdmdp, sdmss, sdmsp, sdmpp
  Common /rndf77/nseed
  Save
!-----------------------------------------------------------------------
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  s = srt**2
  If (sig<=0) Return
!
  If (iabs(lb1)==42) Then
    ideut = i1
    lbm = lb2
    idm = i2
  Else
    ideut = i2
    lbm = lb1
    idm = i1
  End If
!ccc  Elastic collision or destruction of perturbatively-produced deuterons:
  If ((idpert==1 .Or. idpert==2) .And. dpertp(ideut)/=1.) Then
!     choose reaction channels:
    x1 = ranart(nseed)
    If (x1<=sdmel/sig) Then
!     Elastic collisions:
      If (ianti==0) Then
        Write (91, *) '  d+', lbm, ' (pert d M elastic) @nt=', nt, ' @prob=', dpertp(ideut)
      Else
        Write (91, *) '  d+', lbm, ' (pert dbar M elastic) @nt=', nt, ' @prob=', dpertp(ideut)
      End If

!lin-9/2012: check argument in sqrt():
      scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
      If (scheck<0) Then
        Write (99, *) 'scheck51: ', scheck
        scheck = 0.
      End If
      pfinal = sqrt(scheck)/2./srt
!            pfinal=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt

      Call dmelangle(pxn, pyn, pzn, pfinal)
      Call rotate(px, py, pz, pxn, pyn, pzn)
      edcm = sqrt(e(ideut)**2+pxn**2+pyn**2+pzn**2)
      pdbeta = pxn*betax + pyn*betay + pzn*betaz
      transf = gamma*(gamma*pdbeta/(gamma+1.)+edcm)
      pt1d = betax*transf + pxn
      pt2d = betay*transf + pyn
      pt3d = betaz*transf + pzn
      p(1, ideut) = pt1d
      p(2, ideut) = pt2d
      p(3, ideut) = pt3d
      iblock = 504
      px1 = p(1, i1)
      py1 = p(2, i1)
      pz1 = p(3, i1)
      id(i1) = 2
      id(i2) = 2
!     Change the position of the perturbative deuteron to that of
!     the meson to avoid consecutive collisions between them:
      r(1, ideut) = r(1, idm)
      r(2, ideut) = r(2, idm)
      r(3, ideut) = r(3, idm)
    Else
!     Destruction of deuterons:
      If (ianti==0) Then
        Write (91, *) '  d+', lbm, ' ->BB (pert d destrn) @nt=', nt, ' @prob=', dpertp(ideut)
      Else
        Write (91, *) '  d+', lbm, ' ->BB (pert dbar destrn) @nt=', nt, ' @prob=', dpertp(ideut)
      End If
      e(ideut) = 0.
      iblock = 502
    End If
    Return
  End If
!
!ccc  Destruction of regularly-produced deuterons:
  iblock = 502
!     choose final state and assign masses here:
  x1 = ranart(nseed)
  If (x1<=sdmnn/sig) Then
    lbb1 = lbnn1
    lbb2 = lbnn2
    xmb1 = xmnn1
    xmb2 = xmnn2
  Else If (x1<=(sdmnn+sdmnd)/sig) Then
    lbb1 = lbnd1
    lbb2 = lbnd2
    xmb1 = xmnd1
    xmb2 = xmnd2
  Else If (x1<=(sdmnn+sdmnd+sdmns)/sig) Then
    lbb1 = lbns1
    lbb2 = lbns2
    xmb1 = xmns1
    xmb2 = xmns2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp)/sig) Then
    lbb1 = lbnp1
    lbb2 = lbnp2
    xmb1 = xmnp1
    xmb2 = xmnp2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd)/sig) Then
    lbb1 = lbdd1
    lbb2 = lbdd2
    xmb1 = xmdd1
    xmb2 = xmdd2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds)/sig) Then
    lbb1 = lbds1
    lbb2 = lbds2
    xmb1 = xmds1
    xmb2 = xmds2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp)/sig) Then
    lbb1 = lbdp1
    lbb2 = lbdp2
    xmb1 = xmdp1
    xmb2 = xmdp2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp+sdmss)/sig) Then
    lbb1 = lbss1
    lbb2 = lbss2
    xmb1 = xmss1
    xmb2 = xmss2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp+sdmss+sdmsp)/sig) Then
    lbb1 = lbsp1
    lbb2 = lbsp2
    xmb1 = xmsp1
    xmb2 = xmsp2
  Else If (x1<=(sdmnn+sdmnd+sdmns+sdmnp+sdmdd+sdmds+sdmdp+sdmss+sdmsp+sdmpp)/sig) Then
    lbb1 = lbpp1
    lbb2 = lbpp2
    xmb1 = xmpp1
    xmb2 = xmpp2
  Else
!     Elastic collision:
    lbb1 = lb1
    lbb2 = lb2
    xmb1 = em1
    xmb2 = em2
    iblock = 504
  End If
  lb(i1) = lbb1
  e(i1) = xmb1
  lb(i2) = lbb2
  e(i2) = xmb2
  lb1 = lb(i1)
  lb2 = lb(i2)

!lin-9/2012: check argument in sqrt():
  scheck = (s-(xmb1+xmb2)**2)*(s-(xmb1-xmb2)**2)
  If (scheck<0) Then
    Write (99, *) 'scheck52: ', scheck
    scheck = 0.
  End If
  pfinal = sqrt(scheck)/2./srt
!      pfinal=sqrt((s-(xmb1+xmb2)**2)*(s-(xmb1-xmb2)**2))/2./srt

  If (iblock==502) Then
    Call dmangle(pxn, pyn, pzn, nt, ianti, pfinal, lbm)
  Else If (iblock==504) Then
    If (ianti==0) Then
      Write (91, *) ' d+', lbm, ' (regular d M elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
    Else
      Write (91, *) ' d+', lbm, ' (regular dbar M elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
    End If
    Call dmelangle(pxn, pyn, pzn, pfinal)
  Else
    Print *, 'Wrong iblock number in crdmbb()'
    Stop
  End If
!     ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
!     (This is not needed for isotropic distributions)
  Call rotate(px, py, pz, pxn, pyn, pzn)
!     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE
!     FROM THE NUCLEUS-NUCLEUS CMS. FRAME INTO LAB FRAME:
!     For the 1st baryon:
  e1cm = sqrt(e(i1)**2+pxn**2+pyn**2+pzn**2)
  p1beta = pxn*betax + pyn*betay + pzn*betaz
  transf = gamma*(gamma*p1beta/(gamma+1.)+e1cm)
  pt1i1 = betax*transf + pxn
  pt2i1 = betay*transf + pyn
  pt3i1 = betaz*transf + pzn
!
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
!     For the 2nd baryon:
  e2cm = sqrt(e(i2)**2+pxn**2+pyn**2+pzn**2)
  p2beta = -pxn*betax - pyn*betay - pzn*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf - pxn
  pt2i2 = betay*transf - pyn
  pt3i2 = betaz*transf - pzn
!
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
!
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  em2 = e(i2)
  id(i1) = 2
  id(i2) = 2
  Return
End Subroutine crdmbb
!
!     Generate angular distribution of BB from d+meson in the CMS frame:
Subroutine dmangle(pxn, pyn, pzn, nt, ianti, pfinal, lbm)
  Parameter (pi=3.1415926)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /arevt/iaevt, iarun, miss
  Common /rndf77/nseed
  Save
!     take isotropic distribution for now:
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pzn = pfinal*c1
  pxn = pfinal*s1*ct1
  pyn = pfinal*s1*st1
!lin-5/2008 track the number of regularly-destructed deuterons:
  If (ianti==0) Then
    Write (91, *) ' d+', lbm, ' ->BB (regular d destrn) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  Else
    Write (91, *) ' d+', lbm, ' ->BB (regular dbar destrn) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  End If
!
  Return
End Subroutine dmangle
!
!     Angular distribution of d+meson elastic collisions in the CMS frame:
Subroutine dmelangle(pxn, pyn, pzn, pfinal)
  Parameter (pi=3.1415926)
  Common /rndf77/nseed
  Save
!     take isotropic distribution for now:
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pzn = pfinal*c1
  pxn = pfinal*s1*ct1
  pyn = pfinal*s1*st1
  Return
End Subroutine dmelangle
!
!lin-9/2008 Deuteron+Baryon elastic cross section (in mb)
Subroutine sdbelastic(srt, sdb)
  Parameter (srt0=2.012)
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Save
!
  sdb = 0.
  sdbel = 0.
  If (srt<=(em1+em2)) Return
  s = srt**2
!     For elastic collisions:
  If (idxsec==1 .Or. idxsec==3) Then
!     1/3: assume the same |matrix element|**2/s (after averaging over initial
!     spins and isospins) for d+Baryon elastic at the same sqrt(s);
    sdbel = fdbel(s)
  Else If (idxsec==2 .Or. idxsec==4) Then
!     2/4: assume the same |matrix element|**2/s (after averaging over initial
!     spins and isospins) for d+Baryon elastic at the same sqrt(s)-threshold:
    threshold = em1 + em2
    snew = (srt-threshold+srt0)**2
    sdbel = fdbel(snew)
  End If
  sdb = sdbel
  Return
End Subroutine sdbelastic
!lin-9/2008 Deuteron+Baryon elastic collisions
Subroutine crdbel(px, py, pz, srt, i1, i2, iblock, ntag, sig, nt, ianti)
  Parameter (maxstr=150001, maxr=1)
  Common /aa/r(3, maxstr)
  Common /bb/p(3, maxstr)
  Common /bg/betax, betay, betaz, gamma
  Common /cc/e(maxstr)
  Common /ee/id(maxstr), lb(maxstr)
  Common /arevt/iaevt, iarun, miss
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  Common /dpi/em2, lb2
  Common /para8/idpert, npertd, idxsec
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  Save
!-----------------------------------------------------------------------
  iblock = 0
  ntag = 0
  em1 = e(i1)
  em2 = e(i2)
  s = srt**2
  If (sig<=0) Return
  iblock = 503
!
  If (iabs(lb1)==42) Then
    ideut = i1
    lbb = lb2
    idb = i2
  Else
    ideut = i2
    lbb = lb1
    idb = i1
  End If
!ccc  Elastic collision of perturbatively-produced deuterons:
  If ((idpert==1 .Or. idpert==2) .And. dpertp(ideut)/=1.) Then
    If (ianti==0) Then
      Write (91, *) '  d+', lbb, ' (pert d B elastic) @nt=', nt, ' @prob=', dpertp(ideut), p(1, idb), p(2, idb), p(1, ideut), p(2, ideut)
    Else
      Write (91, *) '  d+', lbb, ' (pert dbar Bbar elastic) @nt=', nt, ' @prob=', dpertp(ideut), p(1, idb), p(2, idb), p(1, ideut), p(2, ideut)
    End If

!lin-9/2012: check argument in sqrt():
    scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
    If (scheck<0) Then
      Write (99, *) 'scheck53: ', scheck
      scheck = 0.
    End If
    pfinal = sqrt(scheck)/2./srt
!         pfinal=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt

    Call dbelangle(pxn, pyn, pzn, pfinal)
    Call rotate(px, py, pz, pxn, pyn, pzn)
    edcm = sqrt(e(ideut)**2+pxn**2+pyn**2+pzn**2)
    pdbeta = pxn*betax + pyn*betay + pzn*betaz
    transf = gamma*(gamma*pdbeta/(gamma+1.)+edcm)
    pt1d = betax*transf + pxn
    pt2d = betay*transf + pyn
    pt3d = betaz*transf + pzn
    p(1, ideut) = pt1d
    p(2, ideut) = pt2d
    p(3, ideut) = pt3d
    px1 = p(1, i1)
    py1 = p(2, i1)
    pz1 = p(3, i1)
    id(i1) = 2
    id(i2) = 2
!     Change the position of the perturbative deuteron to that of
!     the baryon to avoid consecutive collisions between them:
    r(1, ideut) = r(1, idb)
    r(2, ideut) = r(2, idb)
    r(3, ideut) = r(3, idb)
    Return
  End If
!
!     Elastic collision of regularly-produced deuterons:
  If (ianti==0) Then
    Write (91, *) ' d+', lbb, ' (regular d B elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  Else
    Write (91, *) ' d+', lbb, ' (regular dbar Bbar elastic) @evt#', iaevt, ' @nt=', nt, ' lb1,2=', lb1, lb2
  End If
!lin-9/2012: check argument in sqrt():
  scheck = (s-(em1+em2)**2)*(s-(em1-em2)**2)
  If (scheck<0) Then
    Write (99, *) 'scheck54: ', scheck
    scheck = 0.
  End If
  pfinal = sqrt(scheck)/2./srt
!      pfinal=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2))/2./srt

  Call dbelangle(pxn, pyn, pzn, pfinal)
!     ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
!     (This is not needed for isotropic distributions)
  Call rotate(px, py, pz, pxn, pyn, pzn)
!     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE
!     FROM THE NUCLEUS-NUCLEUS CMS. FRAME INTO LAB FRAME:
!     For the 1st baryon:
  e1cm = sqrt(e(i1)**2+pxn**2+pyn**2+pzn**2)
  p1beta = pxn*betax + pyn*betay + pzn*betaz
  transf = gamma*(gamma*p1beta/(gamma+1.)+e1cm)
  pt1i1 = betax*transf + pxn
  pt2i1 = betay*transf + pyn
  pt3i1 = betaz*transf + pzn
!
  p(1, i1) = pt1i1
  p(2, i1) = pt2i1
  p(3, i1) = pt3i1
!     For the 2nd baryon:
  e2cm = sqrt(e(i2)**2+pxn**2+pyn**2+pzn**2)
  p2beta = -pxn*betax - pyn*betay - pzn*betaz
  transf = gamma*(gamma*p2beta/(gamma+1.)+e2cm)
  pt1i2 = betax*transf - pxn
  pt2i2 = betay*transf - pyn
  pt3i2 = betaz*transf - pzn
!
  p(1, i2) = pt1i2
  p(2, i2) = pt2i2
  p(3, i2) = pt3i2
!
  px1 = p(1, i1)
  py1 = p(2, i1)
  pz1 = p(3, i1)
  em1 = e(i1)
  em2 = e(i2)
  id(i1) = 2
  id(i2) = 2
  Return
End Subroutine crdbel
!
!     Part of the cross section function of NN->Deuteron+Pi (in mb):
Function fnndpi(s)
  Parameter (srt0=2.012)
  If (s<=srt0**2) Then
    fnndpi = 0.
  Else
    fnndpi = 26.*exp(-(s-4.65)**2/0.1) + 4.*exp(-(s-4.65)**2/2.) + 0.28*exp(-(s-6.)**2/10.)
  End If
  Return
End Function fnndpi
!
!     Angular distribution of d+baryon elastic collisions in the CMS frame:
Subroutine dbelangle(pxn, pyn, pzn, pfinal)
  Parameter (pi=3.1415926)
  Common /rndf77/nseed
  Save
!     take isotropic distribution for now:
  c1 = 1.0 - 2.0*ranart(nseed)
  t1 = 2.0*pi*ranart(nseed)
  s1 = sqrt(1.0-c1**2)
  ct1 = cos(t1)
  st1 = sin(t1)
! THE MOMENTUM IN THE CMS IN THE FINAL STATE
  pzn = pfinal*c1
  pxn = pfinal*s1*ct1
  pyn = pfinal*s1*st1
  Return
End Subroutine dbelangle
!
!     Cross section of Deuteron+Pi elastic (in mb):
Function fdpiel(s)
  Parameter (srt0=2.012)
  If (s<=srt0**2) Then
    fdpiel = 0.
  Else
    fdpiel = 63.*exp(-(s-4.67)**2/0.15) + 15.*exp(-(s-6.25)**2/0.3)
  End If
  Return
End Function fdpiel
!
!     Cross section of Deuteron+N elastic (in mb):
Function fdbel(s)
  Parameter (srt0=2.012)
  If (s<=srt0**2) Then
    fdbel = 0.
  Else
    fdbel = 2500.*exp(-(s-7.93)**2/0.003) + 300.*exp(-(s-7.93)**2/0.1) + 10.
  End If
  Return
End Function fdbel
