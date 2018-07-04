!.................... hijing1.383_ampt.f
!     Version 1.383
!     The variables isng in HIJSFT and JL in ATTRAD were not initialized.
!     The version initialize them. (as found by Fernando Marroquim)
!
!
!
!     Version 1.382
!     Nuclear distribution for deuteron is taken as the Hulthen wave
!     function as provided by Brian Cole (Columbia)
!lin     used my own implementation of impact parameter
!lin     & proton-neutron distance within a deuteron.
!
!
!     Version 1.381
!
!     The parameters for Wood-Saxon distribution for deuteron are
!     constrained to give the right rms ratius 2.116 fm
!     (R=0.0, D=0.5882)
!
!
!     Version 1.38
!
!     The following common block is added to record the number of elastic
!     (NELT, NELP) and inelastic (NINT, NINP) participants
!
!        COMMON/HJGLBR/NELT,NINT,NELP,NINP
!        SAVE /HJGLBR/
!
!     Version 1.37
!
!     A bug in the quenching subroutine is corrected. When calculating the
!     distance between two wounded nucleons, the displacement of the
!     impact parameter was not inculded. This bug was discovered by
!     Dr. V.Uzhinskii JINR, Dubna, Russia
!
!
!     Version 1.36
!
!     Modification Oct. 8, 1998. In hijing, log(ran(nseed)) occasionally
!     causes overfloat. It is modified to log(max(ran(nseed),1.0e-20)).
!
!
!     Nothing important has been changed here. A few 'garbage' has been
!     cleaned up here, like common block HJJET3 for the sea quark strings
!     which were originally created to implement the DPM scheme which
!     later was abadoned in the final version. The lines which operate
!     on these data are also deleted in the program.
!
!
!     Version 1.35
!     There are some changes in the program: subroutine HARDJET is now
!     consolidated with HIJHRD. HARDJET is used to re-initiate PYTHIA
!     for the triggered hard processes. Now that is done  altogether
!     with other normal hard processes in modified JETINI. In the new
!     version one calls JETINI every time one calls HIJHRD. In the new
!     version the effect of the isospin of the nucleon on hard processes,
!     especially direct photons is correctly considered.
!     For A+A collisions, one has to initilize pythia
!     separately for each type of collisions, pp, pn,np and nn,
!     or hp and hn for hA collisions. In JETINI we use the following
!     catalogue for different types of collisions:
!     h+h: h+h (itype=1)
!     h+A: h+p (itype=1), h+n (itype=2)
!     A+h: p+h (itype=1), n+h (itype=2)
!     A+A: p+p (itype=1), p+n (itype=2), n+p (itype=3), n+n (itype=4)
!*****************************************************************
!
!
!     Version 1.34
!     Last modification on January 5, 1998. Two mistakes are corrected in
!     function G. A Mistake in the subroutine Parton is also corrected.
!     (These are pointed out by Ysushi Nara).
!
!
!       Last modifcation on April 10, 1996. To conduct final
!       state radiation, PYTHIA reorganize the two scattered
!       partons and their final momenta will be a little
!       different. The summed total momenta of the partons
!       from the final state radiation are stored in HINT1(26-29)
!       and HINT1(36-39) which are little different from
!       HINT1(21-24) and HINT1(41-44).
!
!       Version 1.33
!
!       Last modfication  on September 11, 1995. When HIJING and
!       PYTHIA are initialized, the shadowing is evaluated at
!       b=0 which is the maximum. This will cause overestimate
!       of shadowing for peripheral interactions. To correct this
!       problem, shadowing is set to zero when initializing. Then
!       use these maximum  cross section without shadowing as a
!       normalization of the Monte Carlo. This however increase
!       the computing time. IHNT2(16) is used to indicate whether
!       the sturcture function is called for (IHNT2(16)=1) initialization
!       or for (IHNT2(16)=0)normal collisions simulation
!
!       Last modification on Aagust 28, 1994. Two bugs associate
!       with the impact parameter dependence of the shadowing is
!       corrected.
!
!
!       Last modification on October 14, 1994. One bug is corrected
!       in the direct photon production option in subroutine
!       HIJHRD.( this problem was reported by Jim Carroll and Mike Beddo).
!       Another bug associated with keeping the decay history
!       in the particle information is also corrected.(this problem
!       was reported by Matt Bloomer)
!
!
!       Last modification on July 15, 1994. The option to trig on
!       heavy quark production (charm IHPR2(18)=0 or beauty IHPR2(18)=1)
!       is added. To do this, set IHPR2(3)=3. For inclusive production,
!       one should reset HIPR1(10)=0.0. One can also trig larger pt
!       QQbar production by giving HIPR1(10) a nonvanishing value.
!       The mass of the heavy quark in the calculation of the cross
!       section (HINT1(59)--HINT1(65)) is given by HIPR1(7) (the
!       default is the charm mass D=1.5). We also include a separate
!       K-factor for heavy quark and direct photon production by
!       HIPR1(23)(D=2.0).
!
!       Last modification on May 24, 1994.  The option to
!       retain the information of all particles including those
!       who have decayed is IHPR(21)=1 (default=0). KATT(I,3) is
!       added to contain the line number of the parent particle
!       of the current line which is produced via a decay.
!       KATT(I,4) is the status number of the particle: 11=particle
!       which has decayed; 1=finally produced particle.
!
!
!       Last modification on May 24, 1994( in HIJSFT when valence quark
!       is quenched, the following error is corrected. 1.2*IHNT2(1) -->
!       1.2*IHNT2(1)**0.333333, 1.2*IHNT2(3) -->1.2*IHNT(3)**0.333333)
!
!
!       Last modification on March 16, 1994 (heavy flavor production
!       processes MSUB(81)=1 MSUB(82)=1 have been switched on,
!       charm production is the default, B-quark option is
!       IHPR2(18), when it is switched on, charm quark is
!       automatically off)
!
!
!       Last modification on March 23, 1994 (an error is corrected
!       in the impact parameter dependence of the jet cross section)
!
!       Last modification Oct. 1993 to comply with non-vax
!       machines' compiler
!
!*********************************************
!	LAST MODIFICATION April 5, 1991
!QUARK DISTRIBUTIOIN (1-X)**A/(X**2+C**2/S)**B
!(A=HIPR1(44),B=HIPR1(46),C=HIPR1(45))
! STRING FLIP, VENUS OPTION IHPR2(15)=1,IN WHICH ONE CAN HAVE ONE AND
! TWO COLOR CHANGES, (1-W)**2,W*(1-W),W*(1-W),AND W*2, W=HIPR1(18),
! AMONG PT DISTRIBUTION OF SEA QUARKS IS CONTROLLED BY HIPR1(42)
!
!	gluon jets can form a single string system
!
!	initial state radiation is included
!
!	all QCD subprocesses are included
!
!	direct particles production is included(currently only direct
!		photon)
!
!	Effect of high P_T trigger bias on multiple jets distribution
!
!******************************************************************
!	                        HIJING.10                         *
!	          Heavy Ion Jet INteraction Generator        	  *
!	                           by                       	  *
!		   X. N. Wang      and   M. Gyulassy           	  *
!	 	      Lawrence Berkeley Laboratory		  *
!								  *
!******************************************************************
!
!******************************************************************
! NFP(K,1),NFP(K,2)=flavor of q and di-q, NFP(K,3)=present ID of  *
! proj, NFP(K,4) original ID of proj.  NFP(K,5)=colli status(0=no,*
! 1=elastic,2=the diffrac one in single-diffrac,3= excited string.*
! |NFP(K,6)| is the total # of jet production, if NFP(K,6)<0 it   *
! can not produce jet anymore. NFP(K,10)=valence quarks scattering*
! (0=has not been,1=is going to be, -1=has already been scattered *
! NFP(k,11) total number of interactions this proj has suffered   *
! PP(K,1)=PX,PP(K,2)=PY,PP(K,3)=PZ,PP(K,4)=E,PP(K,5)=M(invariant  *
! mass), PP(K,6,7),PP(K,8,9)=transverse momentum of quark and     *
! diquark,PP(K,10)=PT of the hard scattering between the valence  *
! quarks; PP(K,14,15)=the mass of quark,diquark.       		  *
!******************************************************************
!
!****************************************************************
!
!	SUBROUTINE HIJING
!
!****************************************************************
Subroutine hijing(frame, bmin0, bmax0)

!bz1/25/99
  Parameter (maxptn=400001)
!lin-4/20/01        PARAMETER (MAXSTR = 1600)
  Parameter (maxstr=150001)
!bz1/25/99end
!lin-4/26/01:
  Parameter (maxidl=4001)

!bz1/31/99
  Double Precision gx0, gy0, gz0, ft0, px0, py0, pz0, e0, xmass0
  Double Precision gx5, gy5, gz5, ft5, px5, py5, pz5, e5, xmass5
  Double Precision ataui, zt1, zt2, zt3
  Double Precision xnprod, etprod, xnfrz, etfrz, dnprod, detpro, dnfrz, detfrz
!lin-8/2015:
  Double Precision vxp0, vyp0, vzp0, xstrg0, ystrg0, xstrg, ystrg

!bz1/31/99end

  Character frame*8
  Dimension scip(300, 300), rnip(300, 300), sjip(300, 300), jtp(3), ipcol(90000), itcol(90000)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
!
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
!lin-7/16/03 NINT is a intrinsic fortran function, rename it to NINTHJ
!        COMMON/HJGLBR/NELT,NINT,NELP,NINP
  Common /hjglbr/nelt, ninthj, nelp, ninp
!c      SAVE /HJGLBR/
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
!c      SAVE /HMAIN1/
!lin-4/26/01
!        COMMON/HMAIN2/KATT(130000,4),PATT(130000,4)
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
!c      SAVE /HMAIN2/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
!lin-4/2008
!        COMMON/HJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
!     &       K2SG(900,100),PXSG(900,100),PYSG(900,100),
!     &       PZSG(900,100),PESG(900,100),PMSG(900,100)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /hjjet4/ndr, iadr(maxstr, 2), kfdr(maxstr), pdr(maxstr, 5)
!lin-4/2008:
!        common/xydr/rtdr(900,2)
  Common /xydr/rtdr(maxstr, 2)
!c      SAVE /HJJET4/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
!
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
!c      SAVE /LUDAT1/

!lin-9/29/03 changed name in order to distinguish from /prec2/
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
!cbz11/11/98
!        COMMON /ARPRC/ ITYP(MAXSTR),
!     &     GX(MAXSTR), GY(MAXSTR), GZ(MAXSTR), FT(MAXSTR),
!     &     PX(MAXSTR), PY(MAXSTR), PZ(MAXSTR), EE(MAXSTR),
!     &     XM(MAXSTR)
!c      SAVE /ARPRC/
!cbz11/11/98end

!bz1/25/99
  Common /para1/mul
!c      SAVE /PARA1/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
!c      SAVE /prec1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
!c      SAVE /prec2/
  Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
!c      SAVE /ilist7/
  Common /ilist8/lstrg1(maxptn), lpart1(maxptn)
!c      SAVE /ilist8/
  Common /srec1/nsp, nst, nsi
!c      SAVE /SREC1/
  Common /srec2/ataui(maxstr), zt1(maxstr), zt2(maxstr), zt3(maxstr)
!c      SAVE /SREC2/
!bz1/25/99end

!lin-2/25/00
  Common /frzout/xnprod(30), etprod(30), xnfrz(30), etfrz(30), dnprod(30), detpro(30), dnfrz(30), detfrz(30)
!c      SAVE /frzout/
!lin-4/11/01 soft:
  Common /anim/nevent, isoft, isflag, izpc
!c      SAVE /anim/
!lin-4/25/01 soft3:
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
!c      SAVE /SOFT/
!lin-4/26/01 lepton and photon info:
  Common /noprec/nnozpc, itypn(maxidl), gxn(maxidl), gyn(maxidl), gzn(maxidl), ftn(maxidl), pxn(maxidl), pyn(maxidl), pzn(maxidl), een(maxidl), xmn(maxidl)
!c      SAVE /NOPREC/
!lin-6/22/01:
  Common /lastt/itimeh, bimp
!c      SAVE /lastt/
  Common /arevt/iaevt, iarun, miss
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
!lin-7/2011 ioscar value is needed:
  Common /para7/ioscar, nsmbbbar, nsmmeson
!lin-2/2012 allow random orientation of reaction plane:
  Common /phihj/iphirp, phirp
!lin-8/2015:
  Common /precpa/vxp0(maxptn), vyp0(maxptn), vzp0(maxptn), xstrg0(maxptn), ystrg0(maxptn), xstrg(maxptn), ystrg(maxptn), istrg0(maxptn), istrg(maxptn)
  Save

  bmax = min(bmax0, hipr1(34)+hipr1(35))
  bmin = min(bmin0, bmax)
  If (ihnt2(1)<=1 .And. ihnt2(3)<=1) Then
    bmin = 0.0
    bmax = 2.5*sqrt(hipr1(31)*0.1/hipr1(40))
  End If
!                        ********HIPR1(31) is in mb =0.1fm**2
!*******THE FOLLOWING IS TO SELECT THE COORDINATIONS OF NUCLEONS
!       BOTH IN PROJECTILE AND TARGET NUCLEAR( in fm)
!
  yp(1, 1) = 0.0
  yp(2, 1) = 0.0
  yp(3, 1) = 0.0
  If (ihnt2(1)<=1) Goto 14
  Do kp = 1, ihnt2(1)
    5 r = hirnd(1)
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
!                ********choose theta from uniform cos(theta) distr
    phi = ranart(nseed)*2.0*hipr1(40)
!                ********choose phi form uniform phi distr 0 to 2*pi
    yp(1, kp) = r*sx*cos(phi)
    yp(2, kp) = r*sx*sin(phi)
    yp(3, kp) = r*cx
    If (hipr1(29)==0.0) Goto 10
    Do kp2 = 1, kp - 1
      dnbp1 = (yp(1,kp)-yp(1,kp2))**2
      dnbp2 = (yp(2,kp)-yp(2,kp2))**2
      dnbp3 = (yp(3,kp)-yp(3,kp2))**2
      dnbp = dnbp1 + dnbp2 + dnbp3
      If (dnbp<hipr1(29)*hipr1(29)) Goto 5
!                        ********two neighbors cannot be closer than
!                                HIPR1(29)
    End Do
  10 End Do

!lin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f,
!     but modified [divide by 2, & x(p)=-x(n)]:
!     (Note: hijing1.383.f has corrected this bug in hijing1.382.f)
  If (ihnt2(1)==2) Then
    rnd1 = max(ranart(nseed), 1.0E-20)
    rnd2 = max(ranart(nseed), 1.0E-20)
    rnd3 = max(ranart(nseed), 1.0E-20)
    r = -(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0+4.38*0.85*log(rnd3)/(4.38+0.85))
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
    phi = ranart(nseed)*2.0*hipr1(40)
!     R above is the relative distance between p & n in a deuteron:
    r = r/2.
    yp(1, 1) = r*sx*cos(phi)
    yp(2, 1) = r*sx*sin(phi)
    yp(3, 1) = r*cx
!     p & n has opposite coordinates in the deuteron frame:
    yp(1, 2) = -yp(1, 1)
    yp(2, 2) = -yp(2, 1)
    yp(3, 2) = -yp(3, 1)
  End If

  Do i = 1, ihnt2(1) - 1
    Do j = i + 1, ihnt2(1)
      If (yp(3,i)>yp(3,j)) Goto 12
      y1 = yp(1, i)
      y2 = yp(2, i)
      y3 = yp(3, i)
      yp(1, i) = yp(1, j)
      yp(2, i) = yp(2, j)
      yp(3, i) = yp(3, j)
      yp(1, j) = y1
      yp(2, j) = y2
      yp(3, j) = y3
    12 End Do
  End Do
!
!******************************
  14 yt(1, 1) = 0.0
  yt(2, 1) = 0.0
  yt(3, 1) = 0.0
  If (ihnt2(3)<=1) Goto 24
  Do kt = 1, ihnt2(3)
    15 r = hirnd(2)
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
!                ********choose theta from uniform cos(theta) distr
    phi = ranart(nseed)*2.0*hipr1(40)
!                ********chose phi form uniform phi distr 0 to 2*pi
    yt(1, kt) = r*sx*cos(phi)
    yt(2, kt) = r*sx*sin(phi)
    yt(3, kt) = r*cx
    If (hipr1(29)==0.0) Goto 20
    Do kt2 = 1, kt - 1
      dnbt1 = (yt(1,kt)-yt(1,kt2))**2
      dnbt2 = (yt(2,kt)-yt(2,kt2))**2
      dnbt3 = (yt(3,kt)-yt(3,kt2))**2
      dnbt = dnbt1 + dnbt2 + dnbt3
      If (dnbt<hipr1(29)*hipr1(29)) Goto 15
!                        ********two neighbors cannot be closer than
!                                HIPR1(29)
    End Do
  20 End Do
!
!lin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f,
!     but modified [divide by 2, & x(p)=-x(n)]:
  If (ihnt2(3)==2) Then
    rnd1 = max(ranart(nseed), 1.0E-20)
    rnd2 = max(ranart(nseed), 1.0E-20)
    rnd3 = max(ranart(nseed), 1.0E-20)
    r = -(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0+4.38*0.85*log(rnd3)/(4.38+0.85))
    x = ranart(nseed)
    cx = 2.0*x - 1.0
    sx = sqrt(1.0-cx*cx)
    phi = ranart(nseed)*2.0*hipr1(40)
    r = r/2.
    yt(1, 1) = r*sx*cos(phi)
    yt(2, 1) = r*sx*sin(phi)
    yt(3, 1) = r*cx
    yt(1, 2) = -yt(1, 1)
    yt(2, 2) = -yt(2, 1)
    yt(3, 2) = -yt(3, 1)
  End If
!
  Do i = 1, ihnt2(3) - 1
    Do j = i + 1, ihnt2(3)
      If (yt(3,i)<yt(3,j)) Goto 22
      y1 = yt(1, i)
      y2 = yt(2, i)
      y3 = yt(3, i)
      yt(1, i) = yt(1, j)
      yt(2, i) = yt(2, j)
      yt(3, i) = yt(3, j)
      yt(1, j) = y1
      yt(2, j) = y2
      yt(3, j) = y3
    22 End Do
  End Do

!********************
  24 miss = -1
  50 miss = miss + 1

!lin-6/2009
!        IF(MISS.GT.50) THEN
  If (miss>maxmiss) Then
    Write (6, *) 'infinite loop happened in  HIJING'
    Stop
  End If

!lin-4/30/01:
  itest = 0

  natt = 0
  jatt = 0
  eatt = 0.0
  Call hijini
  nlop = 0
!                        ********Initialize for a new event
  60 nt = 0
  np = 0
  n0 = 0
  n01 = 0
  n10 = 0
  n11 = 0
  nelt = 0
  ninthj = 0
  nelp = 0
  ninp = 0
  nsg = 0
  ncolt = 0

!****        BB IS THE ABSOLUTE VALUE OF IMPACT PARAMETER,BB**2 IS
!       RANDOMLY GENERATED AND ITS ORIENTATION IS RANDOMLY SET
!       BY THE ANGLE PHI  FOR EACH COLLISION.******************
!
  bb = sqrt(bmin**2+ranart(nseed)*(bmax**2-bmin**2))
!bz6/28/99 flow1
!lin-2/2012:
  phi = 0.
  If (iphirp==1) phi = 2.0*hipr1(40)*ranart(nseed)
  phirp = phi
!bz6/28/99 flow1 end
  bbx = bb*cos(phi)
  bby = bb*sin(phi)
  hint1(19) = bb
  hint1(20) = phi
!
  Do jp = 1, ihnt2(1)
    Do jt = 1, ihnt2(3)
      scip(jp, jt) = -1.0
      b2 = (yp(1,jp)+bbx-yt(1,jt))**2 + (yp(2,jp)+bby-yt(2,jt))**2
      r2 = b2*hipr1(40)/hipr1(31)/0.1
!                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
      rrb1 = min((yp(1,jp)**2+yp(2,jp)**2)/1.2**2/real(ihnt2(1))**0.6666667, 1.0)
      rrb2 = min((yt(1,jt)**2+yt(2,jt)**2)/1.2**2/real(ihnt2(3))**0.6666667, 1.0)
      aphx1 = hipr1(6)*4.0/3.0*(ihnt2(1)**0.3333333-1.0)*sqrt(1.0-rrb1)
      aphx2 = hipr1(6)*4.0/3.0*(ihnt2(3)**0.3333333-1.0)*sqrt(1.0-rrb2)
      hint1(18) = hint1(14) - aphx1*hint1(15) - aphx2*hint1(16) + aphx1*aphx2*hint1(17)
      If (ihpr2(14)==0 .Or. (ihnt2(1)==1 .And. ihnt2(3)==1)) Then
        gs = 1.0 - exp(-(hipr1(30)+hint1(18))*romg(r2)/hipr1(31))
        rantot = ranart(nseed)
        If (rantot>gs) Goto 70
        Goto 65
      End If
      gstot0 = 2.0*(1.0-exp(-(hipr1(30)+hint1(18))/hipr1(31)/2.0*romg(0.0)))
      r2 = r2/gstot0
      gs = 1.0 - exp(-(hipr1(30)+hint1(18))/hipr1(31)*romg(r2))
      gstot = 2.0*(1.0-sqrt(1.0-gs))
      rantot = ranart(nseed)*gstot0
      If (rantot>gstot) Goto 70
      If (rantot>gs) Then
        Call hijcsc(jp, jt)
        Goto 70
!                        ********perform elastic collisions
      End If
      65 scip(jp, jt) = r2
      rnip(jp, jt) = rantot
      sjip(jp, jt) = hint1(18)
      ncolt = ncolt + 1
      ipcol(ncolt) = jp
      itcol(ncolt) = jt
    70 End Do
  End Do
!                ********total number interactions proj and targ has
!                                suffered

!lin-5/22/01 write impact parameter:
  bimp = bb
  Write (6, *) '#impact parameter,nlop,ncolt=', bimp, nlop, ncolt

  If (ncolt==0) Then
    nlop = nlop + 1
    If (nlop<=20 .Or. (ihnt2(1)==1 .And. ihnt2(3)==1)) Goto 60
    Return
  End If
!               ********At large impact parameter, there maybe no
!                       interaction at all. For NN collision
!                       repeat the event until interaction happens
!
  If (ihpr2(3)/=0) Then
    nhard = 1 + int(ranart(nseed)*(ncolt-1)+0.5)
    nhard = min(nhard, ncolt)
    jphard = ipcol(nhard)
    jthard = itcol(nhard)
!lin-6/2009 ctest off:
!           write(99,*) IAEVT,NHARD,NCOLT,JPHARD,JTHARD
  End If
!
  If (ihpr2(9)==1) Then
    nmini = 1 + int(ranart(nseed)*(ncolt-1)+0.5)
    nmini = min(nmini, ncolt)
    jpmini = ipcol(nmini)
    jtmini = itcol(nmini)
  End If
!                ********Specifying the location of the hard and
!                        minijet if they are enforced by user
!
  Do jp = 1, ihnt2(1)
    Do jt = 1, ihnt2(3)
      If (scip(jp,jt)==-1.0) Goto 200
      nfp(jp, 11) = nfp(jp, 11) + 1
      nft(jt, 11) = nft(jt, 11) + 1
      If (nfp(jp,5)<=1 .And. nft(jt,5)>1) Then
        np = np + 1
        n01 = n01 + 1
      Else If (nfp(jp,5)>1 .And. nft(jt,5)<=1) Then
        nt = nt + 1
        n10 = n10 + 1
      Else If (nfp(jp,5)<=1 .And. nft(jt,5)<=1) Then
        np = np + 1
        nt = nt + 1
        n0 = n0 + 1
      Else If (nfp(jp,5)>1 .And. nft(jt,5)>1) Then
        n11 = n11 + 1
      End If
      jout = 0
      nfp(jp, 10) = 0
      nft(jt, 10) = 0
!*****************************************************************
      If (ihpr2(8)==0 .And. ihpr2(3)==0) Goto 160
!                ********When IHPR2(8)=0 no jets are produced
      If (nfp(jp,6)<0 .Or. nft(jt,6)<0) Goto 160
!                ********jets can not be produced for (JP,JT)
!                        because not enough energy avaible for
!                                JP or JT
      r2 = scip(jp, jt)
      hint1(18) = sjip(jp, jt)
      tt = romg(r2)*hint1(18)/hipr1(31)
      tts = hipr1(30)*romg(r2)/hipr1(31)
      njet = 0

      If (ihpr2(3)/=0 .And. jp==jphard .And. jt==jthard) Then
        Call jetini(jp, jt, 1)
        Call hijhrd(jp, jt, 0, jflg, 0)
        hint1(26) = hint1(47)
        hint1(27) = hint1(48)
        hint1(28) = hint1(49)
        hint1(29) = hint1(50)
        hint1(36) = hint1(67)
        hint1(37) = hint1(68)
        hint1(38) = hint1(69)
        hint1(39) = hint1(70)
!
        If (abs(hint1(46))>hipr1(11) .And. jflg==2) nfp(jp, 7) = 1
        If (abs(hint1(56))>hipr1(11) .And. jflg==2) nft(jt, 7) = 1
        If (max(abs(hint1(46)),abs(hint1(56)))>hipr1(11) .And. jflg>=3) iasg(nsg, 3) = 1
        ihnt2(9) = ihnt2(14)
        ihnt2(10) = ihnt2(15)
        Do i05 = 1, 5
          hint1(20+i05) = hint1(40+i05)
          hint1(30+i05) = hint1(50+i05)
        End Do
!lin-6/2009 ctest off:
!           write(99,*) jp,jt,IHPR2(3),HIPR1(10),njet,
!     1          ihnt2(9),hint1(21),hint1(22),hint1(23),
!     2          ihnt2(10),hint1(31),hint1(32),hint1(33)
!           write(99,*) ' '
        jout = 1
        If (ihpr2(8)==0) Goto 160
        rrb1 = min((yp(1,jp)**2+yp(2,jp)**2)/1.2**2/real(ihnt2(1))**0.6666667, 1.0)
        rrb2 = min((yt(1,jt)**2+yt(2,jt)**2)/1.2**2/real(ihnt2(3))**0.6666667, 1.0)
        aphx1 = hipr1(6)*4.0/3.0*(ihnt2(1)**0.3333333-1.0)*sqrt(1.0-rrb1)
        aphx2 = hipr1(6)*4.0/3.0*(ihnt2(3)**0.3333333-1.0)*sqrt(1.0-rrb2)
        hint1(65) = hint1(61) - aphx1*hint1(62) - aphx2*hint1(63) + aphx1*aphx2*hint1(64)
        ttrig = romg(r2)*hint1(65)/hipr1(31)
        njet = -1
!                ********subtract the trigger jet from total number
!                        of jet production  to be done since it has
!                                already been produced here
        xr1 = -alog(exp(-ttrig)+ranart(nseed)*(1.0-exp(-ttrig)))
        106 njet = njet + 1
        xr1 = xr1 - alog(max(ranart(nseed),1.0E-20))
        If (xr1<ttrig) Goto 106
        xr = 0.0
        107 njet = njet + 1
        xr = xr - alog(max(ranart(nseed),1.0E-20))
        If (xr<tt-ttrig) Goto 107
        njet = njet - 1
        Goto 112
      End If
!                ********create a hard interaction with specified P_T
!                                 when IHPR2(3)>0
      If (ihpr2(9)==1 .And. jp==jpmini .And. jt==jtmini) Goto 110
!                ********create at least one pair of mini jets
!                        when IHPR2(9)=1
!
!lin-4/15/2010 changed .LT. to .LE. to avoid problem when two sides are equal;
!     this problem may lead to a jet production when there should be none and
!     crash the run; crashes at low energies were reported by P. Bhaduri.
!        IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LT.EXP(-TT)*
!     &                (1.0-EXP(-TTS))) GO TO 160
      If (ihpr2(8)>0 .And. rnip(jp,jt)<=exp(-tt)*(1.0-exp(-tts))) Goto 160
!
!                ********this is the probability for no jet production
      110 xr = -alog(exp(-tt)+ranart(nseed)*(1.0-exp(-tt)))
      111 njet = njet + 1
      xr = xr - alog(max(ranart(nseed),1.0E-20))
      If (xr<tt) Goto 111
      112 njet = min(njet, ihpr2(8))
      If (ihpr2(8)<0) njet = abs(ihpr2(8))
!                ******** Determine number of mini jet production
!
      Do ijet = 1, njet
        Call jetini(jp, jt, 0)
        Call hijhrd(jp, jt, jout, jflg, 1)
!                ********JFLG=1 jets valence quarks, JFLG=2 with
!                        gluon jet, JFLG=3 with q-qbar prod for
!                        (JP,JT). If JFLG=0 jets can not be produced
!                        this time. If JFLG=-1, error occured abandon
!                        this event. JOUT is the total hard scat for
!                        (JP,JT) up to now.
        If (jflg==0) Goto 160
        If (jflg<0) Then
          If (ihpr2(10)/=0) Write (6, *) 'error occured in HIJHRD'
          Goto 50
        End If
        jout = jout + 1
        If (abs(hint1(46))>hipr1(11) .And. jflg==2) nfp(jp, 7) = 1
        If (abs(hint1(56))>hipr1(11) .And. jflg==2) nft(jt, 7) = 1
        If (max(abs(hint1(46)),abs(hint1(56)))>hipr1(11) .And. jflg>=3) iasg(nsg, 3) = 1
!                ******** jet with PT>HIPR1(11) will be quenched
      End Do
      160 Continue

      Call hijsft(jp, jt, jout, ierror)
      If (ierror/=0) Then
        If (ihpr2(10)/=0) Write (6, *) 'error occured in HIJSFT'
        Goto 50
      End If
!
!                ********conduct soft scattering between JP and JT
      jatt = jatt + jout
    200 End Do
  End Do
!
!**************************
!
!lin-6/2009 write out initial minijet information:
!lin-2/2012:
!           call minijet_out(BB)
  Call minijet_out(bb, phirp)
  If (pttrig>0 .And. ntrig==0) Goto 50
!lin-4/2012
!lin-6/2009 write out initial transverse positions of initial nucleons:
!           write(94,*) IAEVT,MISS,IHNT2(1),IHNT2(3)
  Do jp = 1, ihnt2(1)
!lin-6/2009:
!           write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5)
!lin-2/2012:
!       write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5),yp(3,jp)
!lin-4/2012:
!           write(94,203) YP(1,JP)+0.5*BB*cos(phiRP),
!     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
    If (nfp(jp,5)>2) Then
      ninp = ninp + 1
    Else If (nfp(jp,5)==2 .Or. nfp(jp,5)==1) Then
      nelp = nelp + 1
    End If
  End Do
  Do jt = 1, ihnt2(3)
!lin-6/2009 target nucleon # has a minus sign for distinction from projectile:
!           write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5)
!lin-2/2012:
!       write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5),yt(3,jt)
!lin-4/2012:
!           write(94,203) YT(1,JT)-0.5*BB*cos(phiRP),
!     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
    If (nft(jt,5)>2) Then
      ninthj = ninthj + 1
    Else If (nft(jt,5)==2 .Or. nft(jt,5)==1) Then
      nelt = nelt + 1
    End If
  End Do
! 203    format(f10.3,1x,f10.3,2(1x,I5))
! 203    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
!
!*******************************


!********perform jet quenching for jets with PT>HIPR1(11)**********

  If ((ihpr2(8)/=0 .Or. ihpr2(3)/=0) .And. ihpr2(4)>0 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
    Do i = 1, ihnt2(1)
      If (nfp(i,7)==1) Call quench(i, 1)
    End Do
    Do i = 1, ihnt2(3)
      If (nft(i,7)==1) Call quench(i, 2)
    End Do
    Do isg = 1, nsg
      If (iasg(isg,3)==1) Call quench(isg, 3)
    End Do
  End If

!lin*****4/09/01-soft1, default way of treating strings:
  If (isoft==1) Then
!lin-4/16/01 allow fragmentation:
    isflag = 1

!bz1/25/99
!.....transfer data from HIJING to ZPC
    nsp = ihnt2(1)
    nst = ihnt2(3)
    nsi = nsg
    istr = 0
    npar = 0
    Do i = 1, ihnt2(1)
      istr = istr + 1
      Do j = 1, npj(i)
!bz1/27/99
!.....for now only consider gluon cascade
        If (kfpj(i,j)==21) Then
!bz1/27/99end

          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = j
          ityp0(npar) = kfpj(i, j)
!bz6/28/99 flow1
!lin-7/20/01 add dble or sngl to make precisions consistent
!              GX0(NPAR) = YP(1, I)
!lin-2/2012:
!              GX0(NPAR) = dble(YP(1, I) + 0.5 * BB)
          gx0(npar) = dble(yp(1,i)+0.5*bb*cos(phirp))
!bz6/28/99 flow1 end
!              GY0(NPAR) = dble(YP(2, I))
          gy0(npar) = dble(yp(2,i)+0.5*bb*sin(phirp))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(pjpx(i,j))
          py0(npar) = dble(pjpy(i,j))
          pz0(npar) = dble(pjpz(i,j))
          xmass0(npar) = dble(pjpm(i,j))
!              E0(NPAR) = dble(PJPE(I, J))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
!lin-7/20/01-end

!bz1/27/99
!.....end gluon selection
        End If
!bz1/27/99end
      End Do
    End Do
    Do i = 1, ihnt2(3)
      istr = istr + 1
      Do j = 1, ntj(i)
!bz1/27/99
!.....for now only consider gluon cascade
        If (kftj(i,j)==21) Then
!bz1/27/99end
          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = j
          ityp0(npar) = kftj(i, j)
!bz6/28/99 flow1
!lin-7/20/01 add dble or sngl to make precisions consistent
!              GX0(NPAR) = YT(1, I)
!lin-2/2012:
!              GX0(NPAR) = dble(YT(1, I) - 0.5 * BB)
          gx0(npar) = dble(yt(1,i)-0.5*bb*cos(phirp))
!bz6/28/99 flow1 end
!              GY0(NPAR) = dble(YT(2, I))
          gy0(npar) = dble(yt(2,i)-0.5*bb*sin(phirp))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(pjtx(i,j))
          py0(npar) = dble(pjty(i,j))
          pz0(npar) = dble(pjtz(i,j))
          xmass0(npar) = dble(pjtm(i,j))
!              E0(NPAR) = dble(PJTE(I, J))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)

!bz1/27/99
!.....end gluon selection
        End If
!bz1/27/99end
      End Do
    End Do
    Do i = 1, nsg
      istr = istr + 1
      Do j = 1, njsg(i)
!bz1/27/99
!.....for now only consider gluon cascade
        If (k2sg(i,j)==21) Then
!bz1/27/99end
          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = j
          ityp0(npar) = k2sg(i, j)
!lin-7/20/01 add dble or sngl to make precisions consistent:
          gx0(npar) = 0.5D0*dble(yp(1,iasg(i,1))+yt(1,iasg(i,2)))
          gy0(npar) = 0.5D0*dble(yp(2,iasg(i,1))+yt(2,iasg(i,2)))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(pxsg(i,j))
          py0(npar) = dble(pysg(i,j))
          pz0(npar) = dble(pzsg(i,j))
          xmass0(npar) = dble(pmsg(i,j))
!              E0(NPAR) = dble(PESG(I, J))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
!bz1/27/99
!.....end gluon selection
        End If
!bz1/27/99end
      End Do
    End Do
    mul = npar

!bz2/4/99
    Call hjana1
!bz2/4/99end

!lin-6/2009:
    If (ioscar==3) Write (95, *) iaevt, mul
!.....call ZPC for parton cascade
    Call zpcmn

!     write out parton and wounded nucleon information to ana/zpc1.mom:
!lin-6/2009:
!        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
    Write (14, 395) iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj
    Do i = 1, mul
!c           WRITE (14, 411) PX5(I), PY5(I), PZ5(I), ITYP5(I),
!     &        XMASS5(I), E5(I)
      If (dmax1(abs(gx5(i)),abs(gy5(i)),abs(gz5(i)),abs(ft5(i)))<9999) Then
        Write (14, 210) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      Else
!     change format for large numbers:
        Write (14, 211) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      End If

    End Do

!lin-4/09/01:
    itest = itest + 1
! 411    FORMAT(1X, 3F10.3, I6, 2F10.3)
!bz3/19/99 end

!lin-5/2009 ctest off:
!        call frztm(1,1)

!.....transfer data back from ZPC to HIJING
    Do i = 1, mul
      If (lstrg1(i)<=nsp) Then
        nstrg = lstrg1(i)
        npart = lpart1(i)
        kfpj(nstrg, npart) = ityp5(i)
!lin-7/20/01 add dble or sngl to make precisions consistent
        pjpx(nstrg, npart) = sngl(px5(i))
        pjpy(nstrg, npart) = sngl(py5(i))
        pjpz(nstrg, npart) = sngl(pz5(i))
        pjpe(nstrg, npart) = sngl(e5(i))
        pjpm(nstrg, npart) = sngl(xmass5(i))
      Else If (lstrg1(i)<=nsp+nst) Then
        nstrg = lstrg1(i) - nsp
        npart = lpart1(i)
        kftj(nstrg, npart) = ityp5(i)
        pjtx(nstrg, npart) = sngl(px5(i))
        pjty(nstrg, npart) = sngl(py5(i))
        pjtz(nstrg, npart) = sngl(pz5(i))
        pjte(nstrg, npart) = sngl(e5(i))
        pjtm(nstrg, npart) = sngl(xmass5(i))
      Else
        nstrg = lstrg1(i) - nsp - nst
        npart = lpart1(i)
        k2sg(nstrg, npart) = ityp5(i)
        pxsg(nstrg, npart) = sngl(px5(i))
        pysg(nstrg, npart) = sngl(py5(i))
        pzsg(nstrg, npart) = sngl(pz5(i))
        pesg(nstrg, npart) = sngl(e5(i))
        pmsg(nstrg, npart) = sngl(xmass5(i))
      End If
    End Do
!bz1/25/99end

!bz2/4/99
    Call hjana2
!bz2/4/99end

!lin*****4/09/01-soft2, put q+dq+X in strings into ZPC:
  Else If (isoft==2) Then
    nsp = ihnt2(1)
    nst = ihnt2(3)
!lin-4/27/01:
    nsi = nsg
    npar = 0
    istr = 0
!
!lin  No fragmentation to hadrons, only on parton level,
!     and transfer minijet and string data from HIJING to ZPC:
    mstj(1) = 0
!lin-4/12/01 forbid soft radiation before ZPC to avoid small-mass strings,
!     and forbid jet order reversal before ZPC to avoid unphysical flavors:
    ihpr2(1) = 0
    isflag = 0

    If (ihpr2(20)/=0) Then
      Do ntp = 1, 2
        Do jjtp = 1, ihnt2(2*ntp-1)
          istr = istr + 1
! change: do gluon kink only once: either here or in fragmentation.
          Call hijfrg(jjtp, ntp, ierror)
!                 call lulist(1)
          If (ntp==1) Then
! 354                continue
            npj(jjtp) = max0(n-2, 0)

!lin-4/12/01:                    NPJ(jjtp)=MAX0(ipartn-2,0)
          Else
! 355                continue
            ntj(jjtp) = max0(n-2, 0)
!lin-4/12/01:                    NTJ(jjtp)=MAX0(ipartn-2,0)
          End If

          Do ii = 1, n
            npar = npar + 1
            lstrg0(npar) = istr
            lpart0(npar) = ii
            ityp0(npar) = k(ii, 2)
            gz0(npar) = 0D0
            ft0(npar) = 0D0
!lin-7/20/01 add dble or sngl to make precisions consistent
            px0(npar) = dble(p(ii,1))
            py0(npar) = dble(p(ii,2))
            pz0(npar) = dble(p(ii,3))
            xmass0(npar) = dble(p(ii,5))
!                 E0(NPAR) = dble(P(II,4))
            e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
            If (ntp==1) Then
!lin-7/20/01 add dble or sngl to make precisions consistent
!lin-2/2012:
!                    GX0(NPAR) = dble(YP(1, jjtp)+0.5 * BB)
!                    GY0(NPAR) = dble(YP(2, jjtp))
              gx0(npar) = dble(yp(1,jjtp)+0.5*bb*cos(phirp))
              gy0(npar) = dble(yp(2,jjtp)+0.5*bb*sin(phirp))

              iityp = ityp0(npar)
              nstrg = lstrg0(npar)
              If (iityp==2112 .Or. iityp==2212) Then
              Else If ((iityp==1 .Or. iityp==2) .And. (ii==1 .Or. ii==n)) Then
                pp(nstrg, 6) = sngl(px0(npar))
                pp(nstrg, 7) = sngl(py0(npar))
                pp(nstrg, 14) = sngl(xmass0(npar))
              Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (ii==1 .Or. ii==n)) Then
                pp(nstrg, 8) = sngl(px0(npar))
                pp(nstrg, 9) = sngl(py0(npar))
                pp(nstrg, 15) = sngl(xmass0(npar))
              Else
                npart = lpart0(npar) - 1
                kfpj(nstrg, npart) = ityp0(npar)
                pjpx(nstrg, npart) = sngl(px0(npar))
                pjpy(nstrg, npart) = sngl(py0(npar))
                pjpz(nstrg, npart) = sngl(pz0(npar))
                pjpe(nstrg, npart) = sngl(e0(npar))
                pjpm(nstrg, npart) = sngl(xmass0(npar))
              End If
            Else
!lin-2/2012:
!                    GX0(NPAR) = dble(YT(1, jjtp)-0.5 * BB)
!                    GY0(NPAR) = dble(YT(2, jjtp))
              gx0(npar) = dble(yt(1,jjtp)-0.5*bb*cos(phirp))
              gy0(npar) = dble(yt(2,jjtp)-0.5*bb*sin(phirp))
              iityp = ityp0(npar)
              nstrg = lstrg0(npar) - nsp
              If (iityp==2112 .Or. iityp==2212) Then
              Else If ((iityp==1 .Or. iityp==2) .And. (ii==1 .Or. ii==n)) Then
                pt(nstrg, 6) = sngl(px0(npar))
                pt(nstrg, 7) = sngl(py0(npar))
                pt(nstrg, 14) = sngl(xmass0(npar))
              Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (ii==1 .Or. ii==n)) Then
                pt(nstrg, 8) = sngl(px0(npar))
                pt(nstrg, 9) = sngl(py0(npar))
                pt(nstrg, 15) = sngl(xmass0(npar))
              Else
                npart = lpart0(npar) - 1
                kftj(nstrg, npart) = ityp0(npar)
                pjtx(nstrg, npart) = sngl(px0(npar))
                pjty(nstrg, npart) = sngl(py0(npar))
                pjtz(nstrg, npart) = sngl(pz0(npar))
                pjte(nstrg, npart) = sngl(e0(npar))
                pjtm(nstrg, npart) = sngl(xmass0(npar))
              End If
            End If
          End Do
        End Do
      End Do
      Do isg = 1, nsg
        istr = istr + 1
        Call hijfrg(isg, 3, ierror)
!              call lulist(2)
!
        njsg(isg) = n
!
        Do ii = 1, n
          npar = npar + 1
          lstrg0(npar) = istr
          lpart0(npar) = ii
          ityp0(npar) = k(ii, 2)
          gx0(npar) = 0.5D0*dble(yp(1,iasg(isg,1))+yt(1,iasg(isg,2)))
          gy0(npar) = 0.5D0*dble(yp(2,iasg(isg,1))+yt(2,iasg(isg,2)))
          gz0(npar) = 0D0
          ft0(npar) = 0D0
          px0(npar) = dble(p(ii,1))
          py0(npar) = dble(p(ii,2))
          pz0(npar) = dble(p(ii,3))
          xmass0(npar) = dble(p(ii,5))
!                 E0(NPAR) = dble(P(II,4))
          e0(npar) = dsqrt(px0(npar)**2+py0(npar)**2+pz0(npar)**2+xmass0(npar)**2)
        End Do
      End Do
    End If

    mul = npar
!bz2/4/99
    Call hjana1
!bz2/4/99end
!lin-6/2009:
    If (ioscar==3) Write (95, *) iaevt, mul
!.....call ZPC for parton cascade
    Call zpcmn
!bz3/19/99
!lin-6/2009:
!        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
    Write (14, 395) iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj
    itest = itest + 1

    Do i = 1, mul
!           WRITE (14, 311) PX5(I), PY5(I), PZ5(I), ITYP5(I),
!     &        XMASS5(I), E5(I)
!lin-4/2012 write parton freeze-out position in zpc.dat for this test scenario:
!           WRITE (14, 312) PX5(I), PY5(I), PZ5(I), ITYP5(I),
!     &        XMASS5(I), E5(I),LSTRG1(I), LPART1(I)
      If (dmax1(abs(gx5(i)),abs(gy5(i)),abs(gz5(i)),abs(ft5(i)))<9999) Then
        Write (14, 210) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      Else
        Write (14, 211) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      End If
!
    End Do
! 311    FORMAT(1X, 3F10.4, I6, 2F10.4)
! 312    FORMAT(1X, 3F10.3, I6, 2F10.3,1X,I6,1X,I3)
!bz3/19/99 end

!lin-5/2009 ctest off:
!        call frztm(1,1)

!lin-4/13/01 initialize four momenta and invariant mass of strings after ZPC:
    Do nmom = 1, 5
      Do nstrg = 1, nsp
        pp(nstrg, nmom) = 0.
      End Do
      Do nstrg = 1, nst
        pt(nstrg, nmom) = 0.
      End Do
    End Do
!lin-4/13/01-end

    Do i = 1, mul
      iityp = ityp5(i)
      If (lstrg1(i)<=nsp) Then
        nstrg = lstrg1(i)
!     nucleons without interactions:
        If (iityp==2112 .Or. iityp==2212) Then
!lin-7/20/01 add dble or sngl to make precisions consistent
          pp(nstrg, 1) = sngl(px5(i))
          pp(nstrg, 2) = sngl(py5(i))
          pp(nstrg, 3) = sngl(pz5(i))
          pp(nstrg, 4) = sngl(e5(i))
          pp(nstrg, 5) = sngl(xmass5(i))
!     valence quark:
        Else If ((iityp==1 .Or. iityp==2) .And. (lpart1(i)==1 .Or. lpart1(i)==(npj(nstrg)+2))) Then
          pp(nstrg, 6) = sngl(px5(i))
          pp(nstrg, 7) = sngl(py5(i))
          pp(nstrg, 14) = sngl(xmass5(i))
          pp(nstrg, 1) = pp(nstrg, 1) + sngl(px5(i))
          pp(nstrg, 2) = pp(nstrg, 2) + sngl(py5(i))
          pp(nstrg, 3) = pp(nstrg, 3) + sngl(pz5(i))
          pp(nstrg, 4) = pp(nstrg, 4) + sngl(e5(i))
          pp(nstrg, 5) = sqrt(pp(nstrg,4)**2-pp(nstrg,1)**2-pp(nstrg,2)**2-pp(nstrg,3)**2)
!     diquark:
        Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (lpart1(i)==1 .Or. lpart1(i)==(npj(nstrg)+2))) Then
          pp(nstrg, 8) = sngl(px5(i))
          pp(nstrg, 9) = sngl(py5(i))
          pp(nstrg, 15) = sngl(xmass5(i))
          pp(nstrg, 1) = pp(nstrg, 1) + sngl(px5(i))
          pp(nstrg, 2) = pp(nstrg, 2) + sngl(py5(i))
          pp(nstrg, 3) = pp(nstrg, 3) + sngl(pz5(i))
          pp(nstrg, 4) = pp(nstrg, 4) + sngl(e5(i))
          pp(nstrg, 5) = sqrt(pp(nstrg,4)**2-pp(nstrg,1)**2-pp(nstrg,2)**2-pp(nstrg,3)**2)
!     partons in projectile or target strings:
        Else
          npart = lpart1(i) - 1
          kfpj(nstrg, npart) = ityp5(i)
          pjpx(nstrg, npart) = sngl(px5(i))
          pjpy(nstrg, npart) = sngl(py5(i))
          pjpz(nstrg, npart) = sngl(pz5(i))
          pjpe(nstrg, npart) = sngl(e5(i))
          pjpm(nstrg, npart) = sngl(xmass5(i))
        End If
      Else If (lstrg1(i)<=nsp+nst) Then
        nstrg = lstrg1(i) - nsp
        If (iityp==2112 .Or. iityp==2212) Then
          pt(nstrg, 1) = sngl(px5(i))
          pt(nstrg, 2) = sngl(py5(i))
          pt(nstrg, 3) = sngl(pz5(i))
          pt(nstrg, 4) = sngl(e5(i))
          pt(nstrg, 5) = sngl(xmass5(i))
        Else If ((iityp==1 .Or. iityp==2) .And. (lpart1(i)==1 .Or. lpart1(i)==(ntj(nstrg)+2))) Then
          pt(nstrg, 6) = sngl(px5(i))
          pt(nstrg, 7) = sngl(py5(i))
          pt(nstrg, 14) = sngl(xmass5(i))
          pt(nstrg, 1) = pt(nstrg, 1) + sngl(px5(i))
          pt(nstrg, 2) = pt(nstrg, 2) + sngl(py5(i))
          pt(nstrg, 3) = pt(nstrg, 3) + sngl(pz5(i))
          pt(nstrg, 4) = pt(nstrg, 4) + sngl(e5(i))
          pt(nstrg, 5) = sqrt(pt(nstrg,4)**2-pt(nstrg,1)**2-pt(nstrg,2)**2-pt(nstrg,3)**2)
        Else If ((iityp==1103 .Or. iityp==2101 .Or. iityp==2103 .Or. iityp==2203. .Or. iityp==3101 .Or. iityp==3103. .Or. iityp==3201 .Or. iityp==3203 .Or. iityp==3303) .And. (lpart1(i)==1 .Or. lpart1(i)==(ntj(nstrg)+2))) Then
          pt(nstrg, 8) = sngl(px5(i))
          pt(nstrg, 9) = sngl(py5(i))
          pt(nstrg, 15) = sngl(xmass5(i))
          pt(nstrg, 1) = pt(nstrg, 1) + sngl(px5(i))
          pt(nstrg, 2) = pt(nstrg, 2) + sngl(py5(i))
          pt(nstrg, 3) = pt(nstrg, 3) + sngl(pz5(i))
          pt(nstrg, 4) = pt(nstrg, 4) + sngl(e5(i))
          pt(nstrg, 5) = sqrt(pt(nstrg,4)**2-pt(nstrg,1)**2-pt(nstrg,2)**2-pt(nstrg,3)**2)
        Else
          npart = lpart1(i) - 1
          kftj(nstrg, npart) = ityp5(i)
          pjtx(nstrg, npart) = sngl(px5(i))
          pjty(nstrg, npart) = sngl(py5(i))
          pjtz(nstrg, npart) = sngl(pz5(i))
          pjte(nstrg, npart) = sngl(e5(i))
          pjtm(nstrg, npart) = sngl(xmass5(i))
        End If
      Else
        nstrg = lstrg1(i) - nsp - nst
        npart = lpart1(i)
        k2sg(nstrg, npart) = ityp5(i)
        pxsg(nstrg, npart) = sngl(px5(i))
        pysg(nstrg, npart) = sngl(py5(i))
        pzsg(nstrg, npart) = sngl(pz5(i))
        pesg(nstrg, npart) = sngl(e5(i))
        pmsg(nstrg, npart) = sngl(xmass5(i))
      End If
    End Do
!bz1/25/99end

!lin-4/09/01  turn on fragmentation with soft radiation
!     and jet order reversal to form hadrons after ZPC:
    mstj(1) = 1
    ihpr2(1) = 1
    isflag = 1
!lin-4/13/01 allow small mass strings (D=1.5GeV):
    hipr1(1) = 0.94

!bz2/4/99
    Call hjana2
!bz2/4/99end

!lin-4/19/01-soft3, fragment strings, then convert hadrons to partons
!     and input to ZPC:
  Else If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
!lin-4/24/01 normal fragmentation first:
    isflag = 0
!        write(99,*) 'IAEVT,NSG,NDR=',IAEVT,NSG,NDR

    If (ihpr2(20)/=0) Then
      Do isg = 1, nsg
        Call hijfrg(isg, 3, ierror)
!
        nsbst = 1
        idstr = 92
        If (ihpr2(21)==0) Then
          Call luedit(2)
        Else
          551 nsbst = nsbst + 1
          If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 551
          idstr = k(nsbst, 2)
          nsbst = nsbst + 1
        End If

        If (frame=='LAB') Then
          Call hboost
        End If
!                ******** boost back to lab frame(if it was in)
!
        nsbstr = 0
        Do i = nsbst, n
          If (k(i,2)==idstr) Then
            nsbstr = nsbstr + 1
            Goto 560
          End If
          k(i, 4) = nsbstr
          natt = natt + 1
          katt(natt, 1) = k(i, 2)
          katt(natt, 2) = 20
          katt(natt, 4) = k(i, 1)
!     from Yasushi, to avoid violation of array limits:
!                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
!lin-4/2008 to avoid out-of-bound error in K():
!                   IF(K(I,3).EQ.0 .OR.
!     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
!                      KATT(NATT,3)=0
          If (k(i,3)==0) Then
            katt(natt, 3) = 0
          Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
            katt(natt, 3) = 0
!lin-4/2008-end
          Else
            katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
          End If

!       ****** identify the mother particle
          patt(natt, 1) = p(i, 1)
          patt(natt, 2) = p(i, 2)
          patt(natt, 3) = p(i, 3)
          patt(natt, 4) = p(i, 4)
          eatt = eatt + p(i, 4)
          gxar(natt) = 0.5*(yp(1,iasg(isg,1))+yt(1,iasg(isg,2)))
          gyar(natt) = 0.5*(yp(2,iasg(isg,1))+yt(2,iasg(isg,2)))
          gzar(natt) = 0.
          ftar(natt) = 0.
          itypar(natt) = k(i, 2)
          pxar(natt) = p(i, 1)
          pyar(natt) = p(i, 2)
          pzar(natt) = p(i, 3)
          pear(natt) = p(i, 4)
          xmar(natt) = p(i, 5)
!lin-8/2015: record hadron information, to be used for its constituent partons:
          xstrg0(natt) = dble(gxar(natt))
          ystrg0(natt) = dble(gyar(natt))
          istrg0(natt) = isg
!                   write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT),
!     1                  K(I,2),P(I, 1),P(I, 2),P(I, 3)
!bz11/11/98end

        560 End Do
      End Do
!                ********Fragment the q-qbar jets systems *****
!
      jtp(1) = ihnt2(1)
      jtp(2) = ihnt2(3)
      Do ntp = 1, 2
        Do jjtp = 1, jtp(ntp)
          Call hijfrg(jjtp, ntp, ierror)
!
          nsbst = 1
          idstr = 92
          If (ihpr2(21)==0) Then
            Call luedit(2)
          Else
            581 nsbst = nsbst + 1
            If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 581
            idstr = k(nsbst, 2)
            nsbst = nsbst + 1
          End If
          If (frame=='LAB') Then
            Call hboost
          End If
!                ******** boost back to lab frame(if it was in)
!
          nftp = nfp(jjtp, 5)
          If (ntp==2) nftp = 10 + nft(jjtp, 5)
          nsbstr = 0
          Do i = nsbst, n
            If (k(i,2)==idstr) Then
              nsbstr = nsbstr + 1
              Goto 590
            End If
            k(i, 4) = nsbstr
            natt = natt + 1
            katt(natt, 1) = k(i, 2)
            katt(natt, 2) = nftp
            katt(natt, 4) = k(i, 1)
!                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
!lin-4/2008
!                   IF(K(I,3).EQ.0 .OR.
!     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
!                      KATT(NATT,3)=0
            If (k(i,3)==0) Then
              katt(natt, 3) = 0
            Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
              katt(natt, 3) = 0
!lin-4/2008-end
            Else
              katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
            End If

!       ****** identify the mother particle
            patt(natt, 1) = p(i, 1)
            patt(natt, 2) = p(i, 2)
            patt(natt, 3) = p(i, 3)
            patt(natt, 4) = p(i, 4)
            eatt = eatt + p(i, 4)
            If (ntp==1) Then
!lin-2/2012:
!                      GXAR(NATT) = YP(1, jjtp)+0.5 * BB
!                      GYAR(NATT) = YP(2, jjtp)
              gxar(natt) = yp(1, jjtp) + 0.5*bb*cos(phirp)
              gyar(natt) = yp(2, jjtp) + 0.5*bb*sin(phirp)

            Else
!lin-2/2012:
!                      GXAR(NATT) = YT(1, jjtp)-0.5 * BB
!                      GYAR(NATT) = YT(2, jjtp)
              gxar(natt) = yt(1, jjtp) - 0.5*bb*cos(phirp)
              gyar(natt) = yt(2, jjtp) - 0.5*bb*sin(phirp)
            End If
            gzar(natt) = 0.
            ftar(natt) = 0.
            itypar(natt) = k(i, 2)
            pxar(natt) = p(i, 1)
            pyar(natt) = p(i, 2)
            pzar(natt) = p(i, 3)
            pear(natt) = p(i, 4)
            xmar(natt) = p(i, 5)
!lin-8/2015: record hadron information, to be used for its constituent partons:
            xstrg0(natt) = dble(gxar(natt))
            ystrg0(natt) = dble(gyar(natt))
!     String ID is separated for projectile/target strings:
            istrg0(natt) = ntp*10000 + jjtp
!              if(N.eq.nsbst.and.(K(I,2).eq.2112.or.K(I,2).eq.2212)) then
!                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
!     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3),'spectator'
!                   else
!                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
!     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3)
!                   endif
!bz11/11/98end

          590 End Do
        End Do
      End Do
!     ********Fragment the q-qq related string systems
    End If
!lin-4/2008 check for zero NDR value:
    If (ndr>=1) Then
!
      Do i = 1, ndr
        natt = natt + 1
        katt(natt, 1) = kfdr(i)
        katt(natt, 2) = 40
        katt(natt, 3) = 0
        patt(natt, 1) = pdr(i, 1)
        patt(natt, 2) = pdr(i, 2)
        patt(natt, 3) = pdr(i, 3)
        patt(natt, 4) = pdr(i, 4)
        eatt = eatt + pdr(i, 4)
!lin-11/11/03     set direct photons positions and time at formation:
        gxar(natt) = rtdr(i, 1)
        gyar(natt) = rtdr(i, 2)
        gzar(natt) = 0.
        ftar(natt) = 0.
        itypar(natt) = katt(natt, 1)
        pxar(natt) = patt(natt, 1)
        pyar(natt) = patt(natt, 2)
        pzar(natt) = patt(natt, 3)
        pear(natt) = patt(natt, 4)
        xmar(natt) = pdr(i, 5)
      End Do
!lin-4/2008:
    End If
!lin-6/2009
    Call embedhighpt
!
    Call hjana1

!lin-4/19/01 convert hadrons to partons for ZPC (with GX0 given):
    Call htop

!lin-7/03/01 move up, used in zpstrg (otherwise not set and incorrect):
    nsp = 0
    nst = 0
    nsg = natt
    nsi = nsg
!lin-7/03/01-end

!lin-6/2009:
    If (ioscar==3) Write (95, *) iaevt, mul

!.....call ZPC for parton cascade
    Call zpcmn
!lin-6/2009:
!        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
    Write (14, 395) iaevt, miss, mul, bimp, nelp, ninp, nelt, ninthj
    itest = itest + 1

    Do i = 1, mul
!           WRITE (14, 511) PX5(I), PY5(I), PZ5(I), ITYP5(I),
!     &        XMASS5(I), E5(I)
!lin-4/2012 write parton freeze-out position in zpc.dat
!     for string melting version:
!           WRITE (14, 512) ITYP5(I), PX5(I), PY5(I), PZ5(I),
!     &        XMASS5(I), LSTRG1(I), LPART1(I), FT5(I)
      If (dmax1(abs(gx5(i)),abs(gy5(i)),abs(gz5(i)),abs(ft5(i)))<9999) Then
        Write (14, 210) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      Else
        Write (14, 211) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx5(i), gy5(i), gz5(i), ft5(i)
      End If
!
    End Do
! 511    FORMAT(1X, 3F10.4, I6, 2F10.4)
! 512    FORMAT(I6,4(1X,F10.3),1X,I6,1X,I3,1X,F10.3)
! 513    FORMAT(1X, 4F10.4)

!lin-5/2009 ctest off:
!        call frztm(1,1)

!lin  save data after ZPC for fragmentation purpose:
!.....transfer data back from ZPC to HIJING
    Do i = 1, maxstr
      Do j = 1, 3
        k1sgs(i, j) = 0
        k2sgs(i, j) = 0
        pxsgs(i, j) = 0D0
        pysgs(i, j) = 0D0
        pzsgs(i, j) = 0D0
        pesgs(i, j) = 0D0
        pmsgs(i, j) = 0D0
        gxsgs(i, j) = 0D0
        gysgs(i, j) = 0D0
        gzsgs(i, j) = 0D0
        ftsgs(i, j) = 0D0
      End Do
    End Do
    Do i = 1, mul
      iityp = ityp5(i)
      nstrg = lstrg1(i)
      npart = lpart1(i)
      k2sgs(nstrg, npart) = ityp5(i)
      pxsgs(nstrg, npart) = px5(i)
      pysgs(nstrg, npart) = py5(i)
      pzsgs(nstrg, npart) = pz5(i)
      pmsgs(nstrg, npart) = xmass5(i)
!lin-7/20/01 E5(I) does no include the finite parton mass XMASS5(I),
!     so define it anew:
!           PESGS(NSTRG, NPART) = E5(I)
!           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0)
!     1          write(91,*) 'a',PX5(i),PY5(i),XMASS5(i),PZ5(i),E5(i)
      e5(i) = dsqrt(px5(i)**2+py5(i)**2+pz5(i)**2+xmass5(i)**2)
      pesgs(nstrg, npart) = e5(i)
!           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0)
!     1          write(91,*) 'b: new E5(I)=',E5(i)
!lin-7/20/01-end
      gxsgs(nstrg, npart) = gx5(i)
      gysgs(nstrg, npart) = gy5(i)
      gzsgs(nstrg, npart) = gz5(i)
      ftsgs(nstrg, npart) = ft5(i)
    End Do
    Call hjana2

!lin-4/19/01-end

  End If
!lin-4/09/01-end

!
!**************fragment all the string systems in the following*****
!
!********nsbst is where particle information starts
!********nsbstR+1 is the number of strings in fragmentation
!********the number of strings before a line is stored in K(I,4)
!********IDSTR is id number of the string system (91,92 or 93)
!
!lin-4/30/01 convert partons to hadrons after ZPC:
  If (isoft==3 .Or. isoft==4 .Or. isoft==5) Then
    natt = 0
    eatt = 0.
    Call ptoh
    Do i = 1, nnozpc
      natt = natt + 1
      katt(natt, 1) = itypn(i)
      patt(natt, 1) = pxn(i)
      patt(natt, 2) = pyn(i)
      patt(natt, 3) = pzn(i)
      patt(natt, 4) = een(i)
      eatt = eatt + een(i)
      gxar(natt) = gxn(i)
      gyar(natt) = gyn(i)
      gzar(natt) = gzn(i)
      ftar(natt) = ftn(i)
      itypar(natt) = itypn(i)
      pxar(natt) = pxn(i)
      pyar(natt) = pyn(i)
      pzar(natt) = pzn(i)
      pear(natt) = een(i)
      xmar(natt) = xmn(i)
    End Do
    Goto 565
  End If
!lin-4/30/01-end
  If (ihpr2(20)/=0) Then
    Do isg = 1, nsg
      Call hijfrg(isg, 3, ierror)
      If (mstu(24)/=0 .Or. ierror>0) Then
        mstu(24) = 0
        mstu(28) = 0
        If (ihpr2(10)/=0) Then
!                      call lulist(2)
          Write (6, *) 'error occured ISG, repeat the event'
          Write (6, *) isg

        End If
        Goto 50
      End If
!                        ********Check errors
!
      nsbst = 1
      idstr = 92
      If (ihpr2(21)==0) Then
        Call luedit(2)
      Else
        351 nsbst = nsbst + 1
        If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 351
        idstr = k(nsbst, 2)
        nsbst = nsbst + 1
      End If
!
      If (frame=='LAB') Then
        Call hboost
      End If
!                ******** boost back to lab frame(if it was in)
!
      nsbstr = 0
      Do i = nsbst, n
        If (k(i,2)==idstr) Then
          nsbstr = nsbstr + 1
          Goto 360
        End If
        k(i, 4) = nsbstr
        natt = natt + 1
        katt(natt, 1) = k(i, 2)
        katt(natt, 2) = 20
        katt(natt, 4) = k(i, 1)
!                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
!lin-4/2008:
!                   IF(K(I,3).EQ.0 .OR.
!     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
!                      KATT(NATT,3)=0
        If (k(i,3)==0) Then
          katt(natt, 3) = 0
        Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
          katt(natt, 3) = 0
!lin-4/2008-end
        Else
          katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
        End If

!       ****** identify the mother particle
        patt(natt, 1) = p(i, 1)
        patt(natt, 2) = p(i, 2)
        patt(natt, 3) = p(i, 3)
        patt(natt, 4) = p(i, 4)
        eatt = eatt + p(i, 4)

!bz11/11/98
!bz1/25/99
!                   GXAR(NATT) = 0.5 * (YP(1, IASG(ISG, 1)) +
!     &                YT(1, IASG(ISG, 2)))
!                   GYAR(NATT) = 0.5 * (YP(2, IASG(ISG, 1)) +
!     &                YT(2, IASG(ISG, 2)))
        lsg = nsp + nst + isg
        gxar(natt) = sngl(zt1(lsg))
        gyar(natt) = sngl(zt2(lsg))
        gzar(natt) = sngl(zt3(lsg))
        ftar(natt) = sngl(ataui(lsg))
!bz1/25/99end
        itypar(natt) = k(i, 2)
        pxar(natt) = p(i, 1)
        pyar(natt) = p(i, 2)
        pzar(natt) = p(i, 3)
        pear(natt) = p(i, 4)
        xmar(natt) = p(i, 5)
!bz11/11/98end

      360 End Do
    End Do
!                ********Fragment the q-qbar jets systems *****
!
    jtp(1) = ihnt2(1)
    jtp(2) = ihnt2(3)
    Do ntp = 1, 2
      Do jjtp = 1, jtp(ntp)
        Call hijfrg(jjtp, ntp, ierror)
        If (mstu(24)/=0 .Or. ierror>0) Then
          mstu(24) = 0
          mstu(28) = 0
          If (ihpr2(10)/=0) Then
!                  call lulist(2)
            Write (6, *) 'error occured P&T, repeat the event'
            Write (6, *) ntp, jjtp
!lin-6/2009 when this happens, the event will be repeated,
!     and another record for the same event number will be written into
!     zpc.dat, zpc.res, minijet-initial-beforePropagation.dat,
!     parton-initial-afterPropagation.dat, parton-after-coalescence.dat,
!     and parton-collisionsHistory.dat.
          End If
          Goto 50
        End If
!                        ********check errors
!
        nsbst = 1
        idstr = 92
        If (ihpr2(21)==0) Then
          Call luedit(2)
        Else
          381 nsbst = nsbst + 1
          If (k(nsbst,2)<91 .Or. k(nsbst,2)>93) Goto 381
          idstr = k(nsbst, 2)
          nsbst = nsbst + 1
        End If
        If (frame=='LAB') Then
          Call hboost
        End If
!                ******** boost back to lab frame(if it was in)
!
        nftp = nfp(jjtp, 5)
        If (ntp==2) nftp = 10 + nft(jjtp, 5)
        nsbstr = 0
        Do i = nsbst, n
          If (k(i,2)==idstr) Then
            nsbstr = nsbstr + 1
            Goto 390
          End If
          k(i, 4) = nsbstr
          natt = natt + 1
          katt(natt, 1) = k(i, 2)
          katt(natt, 2) = nftp
          katt(natt, 4) = k(i, 1)
!                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
!lin-4/2008:
!                   IF(K(I,3).EQ.0 .OR.
!     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
!                      KATT(NATT,3)=0
          If (k(i,3)==0) Then
            katt(natt, 3) = 0
          Else If (k(i,3)/=0 .And. k(k(i,3),2)==idstr) Then
            katt(natt, 3) = 0
!lin-4/2008-end
          Else
            katt(natt, 3) = natt - i + k(i, 3) + nsbstr - k(k(i,3), 4)
          End If
!       ****** identify the mother particle
          patt(natt, 1) = p(i, 1)
          patt(natt, 2) = p(i, 2)
          patt(natt, 3) = p(i, 3)
          patt(natt, 4) = p(i, 4)
          eatt = eatt + p(i, 4)
!bz11/11/98
!bz1/25/99
!                   IF (NTP .EQ. 1) THEN
!                      GXAR(NATT) = YP(1, jjtp)
!                   ELSE
!                      GXAR(NATT) = YT(1, jjtp)
!                   END IF
!                   IF (NTP .EQ. 1) THEN
!                      GYAR(NATT) = YP(2, jjtp)
!                   ELSE
!                      GYAR(NATT) = YT(2, jjtp)
!                   END IF
          If (ntp==1) Then
            lsg = jjtp
          Else
            lsg = jjtp + nsp
          End If
          gxar(natt) = sngl(zt1(lsg))
          gyar(natt) = sngl(zt2(lsg))
          gzar(natt) = sngl(zt3(lsg))
          ftar(natt) = sngl(ataui(lsg))
!bz1/25/99end
          itypar(natt) = k(i, 2)
          pxar(natt) = p(i, 1)
          pyar(natt) = p(i, 2)
          pzar(natt) = p(i, 3)
          pear(natt) = p(i, 4)
          xmar(natt) = p(i, 5)
!bz11/11/98end

        390 End Do
      End Do
    End Do
!     ********Fragment the q-qq related string systems
  End If

  Do i = 1, ndr
    natt = natt + 1
    katt(natt, 1) = kfdr(i)
    katt(natt, 2) = 40
    katt(natt, 3) = 0
    patt(natt, 1) = pdr(i, 1)
    patt(natt, 2) = pdr(i, 2)
    patt(natt, 3) = pdr(i, 3)
    patt(natt, 4) = pdr(i, 4)
    eatt = eatt + pdr(i, 4)
!lin-11/11/03     set direct photons positions and time at formation:
    gxar(natt) = rtdr(i, 1)
    gyar(natt) = rtdr(i, 2)
    gzar(natt) = 0.
    ftar(natt) = 0.
    itypar(natt) = katt(natt, 1)
    pxar(natt) = patt(natt, 1)
    pyar(natt) = patt(natt, 2)
    pzar(natt) = patt(natt, 3)
    pear(natt) = patt(natt, 4)
    xmar(natt) = pdr(i, 5)
  End Do

!                        ********store the direct-produced particles
!

!lin-4/19/01 soft3:
  565 Continue

  dengy = eatt/(ihnt2(1)*hint1(6)+ihnt2(3)*hint1(7)) - 1.0
  If (abs(dengy)>hipr1(43) .And. ihpr2(20)/=0 .And. ihpr2(21)==0) Then
    If (ihpr2(10)/=0) Write (6, *) 'Energy not conserved, repeat the event'
!                call lulist(1)
    Write (6, *) 'violated:EATT(GeV),NATT,B(fm)=', eatt, natt, bimp
    Goto 50
  End If
  Write (6, *) 'satisfied:EATT(GeV),NATT,B(fm)=', eatt, natt, bimp
  Write (6, *) ' '
!
!lin-4/2012 write out initial transverse positions of initial nucleons:
  Write (94, *) iaevt, miss, ihnt2(1), ihnt2(3), bimp
  Do jp = 1, ihnt2(1)
!lin-12/2012 write out present and original flavor code of nucleons:
!           write(94,243) YP(1,JP)+0.5*BB*cos(phiRP),
!     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
    Write (94, 243) yp(1, jp) + 0.5*bb*cos(phirp), yp(2, jp) + 0.5*bb*sin(phirp), jp, nfp(jp, 5), yp(3, jp), nfp(jp, 3), nfp(jp, 4)
  End Do
  Do jt = 1, ihnt2(3)
! target nucleon # has a minus sign for distinction from projectile:
!lin-12/2012 write out present and original flavor code of nucleons:
!           write(94,243) YT(1,JT)-0.5*BB*cos(phiRP),
!     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
    Write (94, 243) yt(1, jt) - 0.5*bb*cos(phirp), yt(2, jt) - 0.5*bb*sin(phirp), -jt, nft(jt, 5), yt(3, jt), nft(jt, 3), nft(jt, 4)
  End Do
!lin-4/2012-end

  Return
  210 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,F8.2))
  211 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,E8.2))
  395 Format (3I8, F10.4, 4I5)
!lin-12/2012 write out present and original flavor code of nucleons:
! 243    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
  243 Format (F10.3, 1X, F10.3, 2(1X,I5), 1X, F10.3, 2(1X,I5))
End Subroutine hijing



!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn

Subroutine hijset(efrm, frame, proj, targ, iap, izp, iat, izt)
  Character frame*4, proj*4, targ*4, eframe*4
  Double Precision dd1, dd2, dd3, dd4
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hijdat/hidat0(10, 10), hidat(10)
!c      SAVE /HIJDAT/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
!c      SAVE /LUDAT1/
  External fnkick, fnkc2, fnstru, fnstrm, fnstrs
  Save

  Call title
  ihnt2(1) = iap
  ihnt2(2) = izp
  ihnt2(3) = iat
  ihnt2(4) = izt
  ihnt2(5) = 0
  ihnt2(6) = 0
!
  hint1(8) = max(ulmass(2112), ulmass(2212))
  hint1(9) = hint1(8)
!
  If (proj/='A') Then
    If (proj=='P') Then
      ihnt2(5) = 2212
    Else If (proj=='PBAR') Then
      ihnt2(5) = -2212
    Else If (proj=='PI+') Then
      ihnt2(5) = 211
    Else If (proj=='PI-') Then
      ihnt2(5) = -211
    Else If (proj=='K+') Then
      ihnt2(5) = 321
    Else If (proj=='K-') Then
      ihnt2(5) = -321
    Else If (proj=='N') Then
      ihnt2(5) = 2112
    Else If (proj=='NBAR') Then
      ihnt2(5) = -2112
    Else
      Write (6, *) proj, 'wrong or unavailable proj name'
      Stop
    End If
    hint1(8) = ulmass(ihnt2(5))
  End If
  If (targ/='A') Then
    If (targ=='P') Then
      ihnt2(6) = 2212
    Else If (targ=='PBAR') Then
      ihnt2(6) = -2212
    Else If (targ=='PI+') Then
      ihnt2(6) = 211
    Else If (targ=='PI-') Then
      ihnt2(6) = -211
    Else If (targ=='K+') Then
      ihnt2(6) = 321
    Else If (targ=='K-') Then
      ihnt2(6) = -321
    Else If (targ=='N') Then
      ihnt2(6) = 2112
    Else If (targ=='NBAR') Then
      ihnt2(6) = -2112
    Else
      Write (6, *) targ, 'wrong or unavailable targ name'
      Stop
    End If
    hint1(9) = ulmass(ihnt2(6))
  End If

!...Switch off decay of pi0, K0S, Lambda, Sigma+-, Xi0-, Omega-.
  If (ihpr2(12)>0) Then
    Call lugive('MDCY(C221,1)=0')
!lin-11/07/00 no K* decays:
    Call lugive('MDCY(C313,1)=0')
    Call lugive('MDCY(C-313,1)=0')
    Call lugive('MDCY(C323,1)=0')
    Call lugive('MDCY(C-323,1)=0')
!lin-1/04/01 no K0 and K0bar decays so K0L and K0S do not appear,
!     this way the K/Kbar difference is accounted for exactly:
    Call lugive('MDCY(C311,1)=0')
    Call lugive('MDCY(C-311,1)=0')
!lin-11/08/00 no Delta decays:
    Call lugive('MDCY(C1114,1)=0')
    Call lugive('MDCY(C2114,1)=0')
    Call lugive('MDCY(C2214,1)=0')
    Call lugive('MDCY(C2224,1)=0')
    Call lugive('MDCY(C-1114,1)=0')
    Call lugive('MDCY(C-2114,1)=0')
    Call lugive('MDCY(C-2214,1)=0')
    Call lugive('MDCY(C-2224,1)=0')
!lin-11/07/00-end
!bz12/4/98
    Call lugive('MDCY(C213,1)=0')
    Call lugive('MDCY(C-213,1)=0')
    Call lugive('MDCY(C113,1)=0')
    Call lugive('MDCY(C223,1)=0')
    Call lugive('MDCY(C333,1)=0')
!bz12/4/98end
    Call lugive('MDCY(C111,1)=0')
    Call lugive('MDCY(C310,1)=0')
    Call lugive('MDCY(C411,1)=0;MDCY(C-411,1)=0')
    Call lugive('MDCY(C421,1)=0;MDCY(C-421,1)=0')
    Call lugive('MDCY(C431,1)=0;MDCY(C-431,1)=0')
    Call lugive('MDCY(C511,1)=0;MDCY(C-511,1)=0')
    Call lugive('MDCY(C521,1)=0;MDCY(C-521,1)=0')
    Call lugive('MDCY(C531,1)=0;MDCY(C-531,1)=0')
    Call lugive('MDCY(C3122,1)=0;MDCY(C-3122,1)=0')
    Call lugive('MDCY(C3112,1)=0;MDCY(C-3112,1)=0')
    Call lugive('MDCY(C3212,1)=0;MDCY(C-3212,1)=0')
    Call lugive('MDCY(C3222,1)=0;MDCY(C-3222,1)=0')
    Call lugive('MDCY(C3312,1)=0;MDCY(C-3312,1)=0')
    Call lugive('MDCY(C3322,1)=0;MDCY(C-3322,1)=0')
    Call lugive('MDCY(C3334,1)=0;MDCY(C-3334,1)=0')
!lin-7/2011-no HQ(charm or bottom) decays in order to get net-HQ conservation:
    Call lugive('MDCY(C441,1)=0')
    Call lugive('MDCY(C443,1)=0')
    Call lugive('MDCY(C413,1)=0;MDCY(C-413,1)=0')
    Call lugive('MDCY(C423,1)=0;MDCY(C-423,1)=0')
    Call lugive('MDCY(C433,1)=0;MDCY(C-433,1)=0')
    Call lugive('MDCY(C4112,1)=0;MDCY(C-4112,1)=0')
    Call lugive('MDCY(C4114,1)=0;MDCY(C-4114,1)=0')
    Call lugive('MDCY(C4122,1)=0;MDCY(C-4122,1)=0')
    Call lugive('MDCY(C4212,1)=0;MDCY(C-4212,1)=0')
    Call lugive('MDCY(C4214,1)=0;MDCY(C-4214,1)=0')
    Call lugive('MDCY(C4222,1)=0;MDCY(C-4222,1)=0')
    Call lugive('MDCY(C4224,1)=0;MDCY(C-4224,1)=0')
    Call lugive('MDCY(C4132,1)=0;MDCY(C-4132,1)=0')
    Call lugive('MDCY(C4312,1)=0;MDCY(C-4312,1)=0')
    Call lugive('MDCY(C4314,1)=0;MDCY(C-4314,1)=0')
    Call lugive('MDCY(C4232,1)=0;MDCY(C-4232,1)=0')
    Call lugive('MDCY(C4322,1)=0;MDCY(C-4322,1)=0')
    Call lugive('MDCY(C4324,1)=0;MDCY(C-4324,1)=0')
    Call lugive('MDCY(C4332,1)=0;MDCY(C-4332,1)=0')
    Call lugive('MDCY(C4334,1)=0;MDCY(C-4334,1)=0')
    Call lugive('MDCY(C551,1)=0')
    Call lugive('MDCY(C553,1)=0')
    Call lugive('MDCY(C513,1)=0;MDCY(C-513,1)=0')
    Call lugive('MDCY(C523,1)=0;MDCY(C-523,1)=0')
    Call lugive('MDCY(C533,1)=0;MDCY(C-533,1)=0')
    Call lugive('MDCY(C5112,1)=0;MDCY(C-5112,1)=0')
    Call lugive('MDCY(C5114,1)=0;MDCY(C-5114,1)=0')
    Call lugive('MDCY(C5122,1)=0;MDCY(C-5122,1)=0')
    Call lugive('MDCY(C5212,1)=0;MDCY(C-5212,1)=0')
    Call lugive('MDCY(C5214,1)=0;MDCY(C-5214,1)=0')
    Call lugive('MDCY(C5222,1)=0;MDCY(C-5222,1)=0')
    Call lugive('MDCY(C5224,1)=0;MDCY(C-5224,1)=0')
!lin-7/2011-end
  End If
  mstu(12) = 0
  mstu(21) = 1
  If (ihpr2(10)==0) Then
    mstu(22) = 0
    mstu(25) = 0
    mstu(26) = 0
  End If

!lin    parj(41) and (42) are a, b parameters in Lund, read from input.ampt:
!        PARJ(41)=HIPR1(3)
!        PARJ(42)=HIPR1(4)
!        PARJ(41)=2.2
!        PARJ(42)=0.5

!lin  2 popcorn parameters read from input.ampt:
!        IHPR2(11) = 3
!        PARJ(5) = 0.5
  mstj(12) = ihpr2(11)

!lin  parj(21) gives the mean gaussian width for hadron Pt:
  parj(21) = hipr1(2)
!lin  parj(2) is gamma_s=P(s)/P(u), kappa propto 1/b/(2+a) assumed.
  rkp = hipr1(4)*(2+hipr1(3))/parj(42)/(2+parj(41))
  parj(2) = parj(2)**(1./rkp)
  parj(21) = parj(21)*sqrt(rkp)
!lin-10/31/00 update when string tension is changed:
  hipr1(2) = parj(21)

!lin-4/2015: set upper limit for gamma_s=P(s)/P(u) to 0.4
!     (to limit strangeness enhancement when string tension is strongly
!     increased due to using a very low value of parameter b in Lund
!     symmetric splitting function as done in arXiv:1403.6321):
  parj(2) = min(parj(2), 0.4)

!                        ******** set up for jetset
  If (frame=='LAB') Then
    dd1 = dble(efrm)
    dd2 = dble(hint1(8))
    dd3 = dble(hint1(9))
    hint1(1) = sqrt(hint1(8)**2+2.0*hint1(9)*efrm+hint1(9)**2)
    dd4 = dsqrt(dd1**2-dd2**2)/(dd1+dd3)
    hint1(2) = sngl(dd4)
    hint1(3) = 0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    dd4 = dsqrt(dd1**2-dd2**2)/dd1
    hint1(4) = 0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    hint1(5) = 0.0
    hint1(6) = efrm
    hint1(7) = hint1(9)
  Else If (frame=='CMS') Then
    hint1(1) = efrm
    hint1(2) = 0.0
    hint1(3) = 0.0
    dd1 = dble(hint1(1))
    dd2 = dble(hint1(8))
    dd3 = dble(hint1(9))
    dd4 = dsqrt(1.D0-4.D0*dd2**2/dd1**2)
    hint1(4) = 0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    dd4 = dsqrt(1.D0-4.D0*dd3**2/dd1**2)
    hint1(5) = -0.5*sngl(dlog((1.D0+dd4)/(1.D0-dd4)))
    hint1(6) = hint1(1)/2.0
    hint1(7) = hint1(1)/2.0
  End If
!                ********define Lorentz transform to lab frame
!
!                ********calculate the cross sections involved with
!                        nucleon collisions.
  If (ihnt2(1)>1) Then
    Call hijwds(ihnt2(1), 1, rmax)
    hipr1(34) = rmax
!                        ********set up Wood-Sax distr for proj.
  End If
  If (ihnt2(3)>1) Then
    Call hijwds(ihnt2(3), 2, rmax)
    hipr1(35) = rmax
!                        ********set up Wood-Sax distr for  targ.
  End If
!
!
  i = 0
  20 i = i + 1
  If (i==10) Goto 30
  If (hidat0(10,i)<=hint1(1)) Goto 20
  30 If (i==1) i = 2
  Do j = 1, 9
    hidat(j) = hidat0(j, i-1) + (hidat0(j,i)-hidat0(j,i-1))*(hint1(1)-hidat0(10,i-1))/(hidat0(10,i)-hidat0(10,i-1))
  End Do
  hipr1(31) = hidat(5)
  hipr1(30) = 2.0*hidat(5)
!
!
  Call hijcrs
!
  If (ihpr2(5)/=0) Then
    Call hifun(3, 0.0, 36.0, fnkick)
!                ********booking for generating pt**2 for pt kick
  End If
  Call hifun(7, 0.0, 6.0, fnkc2)
  Call hifun(4, 0.0, 1.0, fnstru)
  Call hifun(5, 0.0, 1.0, fnstrm)
  Call hifun(6, 0.0, 1.0, fnstrs)
!                ********booking for x distribution of valence quarks
  eframe = 'Ecm'
  If (frame=='LAB') eframe = 'Elab'
  Write (6, 100) eframe, efrm, proj, ihnt2(1), ihnt2(2), targ, ihnt2(3), ihnt2(4)
  Return
  100 Format (10X, '**************************************************'/10X, '*', 48X, '*'/10X, '*         HIJING has been initialized at         *'/10X, '*', 13X, A4, '= ', F10.2, ' GeV/n', 13X, '*'/10X, '*', 48X, '*'/10X, '*', 8X, 'for ', A4, '(', I3, ',', I3, ')', ' + ', A4, '(', I3, ',', I3, ')', 7X, '*'/10X, '**************************************************')
End Subroutine hijset
!
!
!
Function fnkick(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  fnkick = 1.0/(x+hipr1(19)**2)/(x+hipr1(20)**2)/(1+exp((sqrt(x)-hipr1(20))/0.4))
  Return
End Function fnkick
!
!
Function fnkc2(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  fnkc2 = x*exp(-2.0*x/hipr1(42))
  Return
End Function fnkc2
!
!
!
Function fnstru(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  fnstru = (1.0-x)**hipr1(44)/(x**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)
  Return
End Function fnstru
!
!
!
Function fnstrm(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  fnstrm = 1.0/((1.0-x)**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)/(x**2+hipr1(45)**2/hint1(1)**2)**hipr1(46)
  Return
End Function fnstrm
!
!
Function fnstrs(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  fnstrs = (1.0-x)**hipr1(47)/(x**2+hipr1(45)**2/hint1(1)**2)**hipr1(48)
  Return
End Function fnstrs
!
!
!
!
Subroutine hboost
  Implicit Double Precision (D)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
!c      SAVE /LUDAT1/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  Do i = 1, n
    dbeta = dble(p(i,3)/p(i,4))
    If (abs(dbeta)>=1.D0) Then
      db = dble(hint1(2))
      If (db>0.99999999D0) Then
!                ********Rescale boost vector if too close to unity.
        Write (6, *) '(HIBOOT:) boost vector too large'
        db = 0.99999999D0
      End If
      dga = 1D0/sqrt(1D0-db**2)
      dp3 = dble(p(i,3))
      dp4 = dble(p(i,4))
      p(i, 3) = sngl((dp3+db*dp4)*dga)
      p(i, 4) = sngl((dp4+db*dp3)*dga)
      Goto 100
    End If
    y = 0.5*sngl(dlog((1.D0+dbeta)/(1.D0-dbeta)))
    amt = sqrt(p(i,1)**2+p(i,2)**2+p(i,5)**2)
    p(i, 3) = amt*sinh(y+hint1(3))
    p(i, 4) = amt*cosh(y+hint1(3))
  100 End Do
  Return
End Subroutine hboost
!
!
!
!
Subroutine quench(jpjt, ntp)
  Parameter (maxstr=150001)
  Dimension rdp(300), lqp(300), rdt(300), lqt(300)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
!
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
!     Uzhi:
  bb = hint1(19)
  phi = hint1(20)
  bbx = bb*cos(phi)
  bby = bb*sin(phi)
!
  If (ntp==2) Goto 400
  If (ntp==3) Goto 2000
!*******************************************************
! Jet interaction for proj jet in the direction PHIP
!******************************************************
!
  If (nfp(jpjt,7)/=1) Return

  jp = jpjt
  Do i = 1, npj(jp)
    ptjet0 = sqrt(pjpx(jp,i)**2+pjpy(jp,i)**2)
    If (ptjet0<=hipr1(11)) Goto 290
    ptot = sqrt(ptjet0*ptjet0+pjpz(jp,i)**2)
    If (ptot<hipr1(8)) Goto 290
    phip = ulangl(pjpx(jp,i), pjpy(jp,i))
!******* find the wounded proj which can interact with jet***
    kp = 0
    Do i2 = 1, ihnt2(1)
      If (nfp(i2,5)/=3 .Or. i2==jp) Goto 100
      dx = yp(1, i2) - yp(1, jp)
      dy = yp(2, i2) - yp(2, jp)
      phi = ulangl(dx, dy)
      dphi = abs(phi-phip)
!     Uzhi:
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>=hipr1(40)/2.0) Goto 100
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 100
      kp = kp + 1
      lqp(kp) = i2
      rdp(kp) = cos(dphi)*rd0
    100 End Do
!*******        rearrange according decending rd************
    Do i2 = 1, kp - 1
      Do j2 = i2 + 1, kp
        If (rdp(i2)<rdp(j2)) Goto 110
        rd = rdp(i2)
        lq = lqp(i2)
        rdp(i2) = rdp(j2)
        lqp(i2) = lqp(j2)
        rdp(j2) = rd
        lqp(j2) = lq
      110 End Do
    End Do
!****** find wounded targ which can interact with jet********
    kt = 0
    Do i2 = 1, ihnt2(3)
      If (nft(i2,5)/=3) Goto 120
      dx = yt(1, i2) - yp(1, jp) - bbx
      dy = yt(2, i2) - yp(2, jp) - bby
      phi = ulangl(dx, dy)
      dphi = abs(phi-phip)
!     Uzhi:
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 120
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 120
      kt = kt + 1
      lqt(kt) = i2
      rdt(kt) = cos(dphi)*rd0
    120 End Do
!*******        rearrange according decending rd************
    Do i2 = 1, kt - 1
      Do j2 = i2 + 1, kt
        If (rdt(i2)<rdt(j2)) Goto 130
        rd = rdt(i2)
        lq = lqt(i2)
        rdt(i2) = rdt(j2)
        lqt(i2) = lqt(j2)
        rdt(j2) = rd
        lqt(j2) = lq
      130 End Do
    End Do

    mp = 0
    mt = 0
    r0 = 0.0
    nq = 0
    dp = 0.0
    ptot = sqrt(pjpx(jp,i)**2+pjpy(jp,i)**2+pjpz(jp,i)**2)
    v1 = pjpx(jp, i)/ptot
    v2 = pjpy(jp, i)/ptot
    v3 = pjpz(jp, i)/ptot

    200 rn = ranart(nseed)
    210 If (mt>=kt .And. mp>=kp) Goto 290
    If (mt>=kt) Goto 220
    If (mp>=kp) Goto 240
    If (rdp(mp+1)>rdt(mt+1)) Goto 240
    220 mp = mp + 1
    drr = rdp(mp) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 210
    dp = drr*hipr1(14)
    If (kfpj(jp,i)/=21) dp = 0.5*dp
!        ********string tension of quark jet is 0.5 of gluon's
    If (dp<=0.2) Goto 210
    If (ptot<=0.4) Goto 290
    If (ptot<=dp) dp = ptot - 0.2
    de = dp

    If (kfpj(jp,i)/=21) Then
      prshu = pp(lqp(mp), 1)**2 + pp(lqp(mp), 2)**2 + pp(lqp(mp), 3)**2
      de = sqrt(pjpm(jp,i)**2+ptot**2) - sqrt(pjpm(jp,i)**2+(ptot-dp)**2)
      ershu = (pp(lqp(mp),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 210
      pp(lqp(mp), 4) = sqrt(ershu)
      pp(lqp(mp), 5) = sqrt(amshu)
    End If
!                ********reshuffle the energy when jet has mass
    r0 = rdp(mp)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
!                ********momentum and energy transfer from jet

    npj(lqp(mp)) = npj(lqp(mp)) + 1
    kfpj(lqp(mp), npj(lqp(mp))) = 21
    pjpx(lqp(mp), npj(lqp(mp))) = dp1
    pjpy(lqp(mp), npj(lqp(mp))) = dp2
    pjpz(lqp(mp), npj(lqp(mp))) = dp3
    pjpe(lqp(mp), npj(lqp(mp))) = dp
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0
    Goto 260

    240 mt = mt + 1
    drr = rdt(mt) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 210
    dp = drr*hipr1(14)
    If (dp<=0.2) Goto 210
    If (ptot<=0.4) Goto 290
    If (ptot<=dp) dp = ptot - 0.2
    de = dp

    If (kfpj(jp,i)/=21) Then
      prshu = pt(lqt(mt), 1)**2 + pt(lqt(mt), 2)**2 + pt(lqt(mt), 3)**2
      de = sqrt(pjpm(jp,i)**2+ptot**2) - sqrt(pjpm(jp,i)**2+(ptot-dp)**2)
      ershu = (pt(lqt(mt),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 210
      pt(lqt(mt), 4) = sqrt(ershu)
      pt(lqt(mt), 5) = sqrt(amshu)
    End If
!                ********reshuffle the energy when jet has mass

    r0 = rdt(mt)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
!                ********momentum and energy transfer from jet
    ntj(lqt(mt)) = ntj(lqt(mt)) + 1
    kftj(lqt(mt), ntj(lqt(mt))) = 21
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1
    pjty(lqt(mt), ntj(lqt(mt))) = dp2
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3
    pjte(lqt(mt), ntj(lqt(mt))) = dp
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0

    260 pjpx(jp, i) = (ptot-dp)*v1
    pjpy(jp, i) = (ptot-dp)*v2
    pjpz(jp, i) = (ptot-dp)*v3
    pjpe(jp, i) = pjpe(jp, i) - de

    ptot = ptot - dp
    nq = nq + 1
    Goto 200
  290 End Do

  Return

!*******************************************************
! Jet interaction for target jet in the direction PHIT
!******************************************************
!
!******* find the wounded proj which can interact with jet***

  400 If (nft(jpjt,7)/=1) Return
  jt = jpjt
  Do i = 1, ntj(jt)
    ptjet0 = sqrt(pjtx(jt,i)**2+pjty(jt,i)**2)
    If (ptjet0<=hipr1(11)) Goto 690
    ptot = sqrt(ptjet0*ptjet0+pjtz(jt,i)**2)
    If (ptot<hipr1(8)) Goto 690
    phit = ulangl(pjtx(jt,i), pjty(jt,i))
    kp = 0
    Do i2 = 1, ihnt2(1)
      If (nfp(i2,5)/=3) Goto 500
      dx = yp(1, i2) + bbx - yt(1, jt)
      dy = yp(2, i2) + bby - yt(2, jt)
      phi = ulangl(dx, dy)
      dphi = abs(phi-phit)
!     Uzhi:
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 500
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 500
      kp = kp + 1
      lqp(kp) = i2
      rdp(kp) = cos(dphi)*rd0
    500 End Do
!*******        rearrange according to decending rd************
    Do i2 = 1, kp - 1
      Do j2 = i2 + 1, kp
        If (rdp(i2)<rdp(j2)) Goto 510
        rd = rdp(i2)
        lq = lqp(i2)
        rdp(i2) = rdp(j2)
        lqp(i2) = lqp(j2)
        rdp(j2) = rd
        lqp(j2) = lq
      510 End Do
    End Do
!****** find wounded targ which can interact with jet********
    kt = 0
    Do i2 = 1, ihnt2(3)
      If (nft(i2,5)/=3 .Or. i2==jt) Goto 520
      dx = yt(1, i2) - yt(1, jt)
      dy = yt(2, i2) - yt(2, jt)
      phi = ulangl(dx, dy)
      dphi = abs(phi-phit)
!     Uzhi:
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 520
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 520
      kt = kt + 1
      lqt(kt) = i2
      rdt(kt) = cos(dphi)*rd0
    520 End Do
!*******        rearrange according to decending rd************
    Do i2 = 1, kt - 1
      Do j2 = i2 + 1, kt
        If (rdt(i2)<rdt(j2)) Goto 530
        rd = rdt(i2)
        lq = lqt(i2)
        rdt(i2) = rdt(j2)
        lqt(i2) = lqt(j2)
        rdt(j2) = rd
        lqt(j2) = lq
      530 End Do
    End Do

    mp = 0
    mt = 0
    nq = 0
    dp = 0.0
    r0 = 0.0
    ptot = sqrt(pjtx(jt,i)**2+pjty(jt,i)**2+pjtz(jt,i)**2)
    v1 = pjtx(jt, i)/ptot
    v2 = pjty(jt, i)/ptot
    v3 = pjtz(jt, i)/ptot

    600 rn = ranart(nseed)
    610 If (mt>=kt .And. mp>=kp) Goto 690
    If (mt>=kt) Goto 620
    If (mp>=kp) Goto 640
    If (rdp(mp+1)>rdt(mt+1)) Goto 640
    620 mp = mp + 1
    drr = rdp(mp) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 610
    dp = drr*hipr1(14)
    If (kftj(jt,i)/=21) dp = 0.5*dp
!        ********string tension of quark jet is 0.5 of gluon's
    If (dp<=0.2) Goto 610
    If (ptot<=0.4) Goto 690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
!
    If (kftj(jt,i)/=21) Then
      prshu = pp(lqp(mp), 1)**2 + pp(lqp(mp), 2)**2 + pp(lqp(mp), 3)**2
      de = sqrt(pjtm(jt,i)**2+ptot**2) - sqrt(pjtm(jt,i)**2+(ptot-dp)**2)
      ershu = (pp(lqp(mp),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 610
      pp(lqp(mp), 4) = sqrt(ershu)
      pp(lqp(mp), 5) = sqrt(amshu)
    End If
!                ********reshuffle the energy when jet has mass
!
    r0 = rdp(mp)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
!                ********momentum and energy transfer from jet
    npj(lqp(mp)) = npj(lqp(mp)) + 1
    kfpj(lqp(mp), npj(lqp(mp))) = 21
    pjpx(lqp(mp), npj(lqp(mp))) = dp1
    pjpy(lqp(mp), npj(lqp(mp))) = dp2
    pjpz(lqp(mp), npj(lqp(mp))) = dp3
    pjpe(lqp(mp), npj(lqp(mp))) = dp
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0

    Goto 660

    640 mt = mt + 1
    drr = rdt(mt) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 610
    dp = drr*hipr1(14)
    If (dp<=0.2) Goto 610
    If (ptot<=0.4) Goto 690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp

    If (kftj(jt,i)/=21) Then
      prshu = pt(lqt(mt), 1)**2 + pt(lqt(mt), 2)**2 + pt(lqt(mt), 3)**2
      de = sqrt(pjtm(jt,i)**2+ptot**2) - sqrt(pjtm(jt,i)**2+(ptot-dp)**2)
      ershu = (pt(lqt(mt),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 610
      pt(lqt(mt), 4) = sqrt(ershu)
      pt(lqt(mt), 5) = sqrt(amshu)
    End If
!                ********reshuffle the energy when jet has mass

    r0 = rdt(mt)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
!                ********momentum and energy transfer from jet
    ntj(lqt(mt)) = ntj(lqt(mt)) + 1
    kftj(lqt(mt), ntj(lqt(mt))) = 21
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1
    pjty(lqt(mt), ntj(lqt(mt))) = dp2
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3
    pjte(lqt(mt), ntj(lqt(mt))) = dp
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0

    660 pjtx(jt, i) = (ptot-dp)*v1
    pjty(jt, i) = (ptot-dp)*v2
    pjtz(jt, i) = (ptot-dp)*v3
    pjte(jt, i) = pjte(jt, i) - de

    ptot = ptot - dp
    nq = nq + 1
    Goto 600
  690 End Do
  Return
!********************************************************
!        Q-QBAR jet interaction
!********************************************************
  2000 isg = jpjt
  If (iasg(isg,3)/=1) Return
!
  jp = iasg(isg, 1)
  jt = iasg(isg, 2)
  xj = (yp(1,jp)+bbx+yt(1,jt))/2.0
  yj = (yp(2,jp)+bby+yt(2,jt))/2.0
  Do i = 1, njsg(isg)
    ptjet0 = sqrt(pxsg(isg,i)**2+pysg(isg,i)**2)
    If (ptjet0<=hipr1(11) .Or. pesg(isg,i)<hipr1(1)) Goto 2690
    ptot = sqrt(ptjet0*ptjet0+pzsg(isg,i)**2)
    If (ptot<max(hipr1(1),hipr1(8))) Goto 2690
    phiq = ulangl(pxsg(isg,i), pysg(isg,i))
    kp = 0
    Do i2 = 1, ihnt2(1)
      If (nfp(i2,5)/=3 .Or. i2==jp) Goto 2500
      dx = yp(1, i2) + bbx - xj
      dy = yp(2, i2) + bby - yj
      phi = ulangl(dx, dy)
      dphi = abs(phi-phiq)
!     Uzhi:
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 2500
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 2500
      kp = kp + 1
      lqp(kp) = i2
      rdp(kp) = cos(dphi)*rd0
    2500 End Do
!*******        rearrange according to decending rd************
    Do i2 = 1, kp - 1
      Do j2 = i2 + 1, kp
        If (rdp(i2)<rdp(j2)) Goto 2510
        rd = rdp(i2)
        lq = lqp(i2)
        rdp(i2) = rdp(j2)
        lqp(i2) = lqp(j2)
        rdp(j2) = rd
        lqp(j2) = lq
      2510 End Do
    End Do
!****** find wounded targ which can interact with jet********
    kt = 0
    Do i2 = 1, ihnt2(3)
      If (nft(i2,5)/=3 .Or. i2==jt) Goto 2520
      dx = yt(1, i2) - xj
      dy = yt(2, i2) - yj
      phi = ulangl(dx, dy)
      dphi = abs(phi-phiq)
!     Uzhi:
      If (dphi>=hipr1(40)) dphi = 2.*hipr1(40) - dphi
      If (dphi>hipr1(40)/2.0) Goto 2520
      rd0 = sqrt(dx*dx+dy*dy)
      If (rd0*sin(dphi)>hipr1(12)) Goto 2520
      kt = kt + 1
      lqt(kt) = i2
      rdt(kt) = cos(dphi)*rd0
    2520 End Do
!*******        rearrange according to decending rd************
    Do i2 = 1, kt - 1
      Do j2 = i2 + 1, kt
        If (rdt(i2)<rdt(j2)) Goto 2530
        rd = rdt(i2)
        lq = lqt(i2)
        rdt(i2) = rdt(j2)
        lqt(i2) = lqt(j2)
        rdt(j2) = rd
        lqt(j2) = lq
      2530 End Do
    End Do

    mp = 0
    mt = 0
    nq = 0
    dp = 0.0
    r0 = 0.0
    ptot = sqrt(pxsg(isg,i)**2+pysg(isg,i)**2+pzsg(isg,i)**2)
    v1 = pxsg(isg, i)/ptot
    v2 = pysg(isg, i)/ptot
    v3 = pzsg(isg, i)/ptot

    2600 rn = ranart(nseed)
    2610 If (mt>=kt .And. mp>=kp) Goto 2690
    If (mt>=kt) Goto 2620
    If (mp>=kp) Goto 2640
    If (rdp(mp+1)>rdt(mt+1)) Goto 2640
    2620 mp = mp + 1
    drr = rdp(mp) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 2610
    dp = drr*hipr1(14)/2.0
    If (dp<=0.2) Goto 2610
    If (ptot<=0.4) Goto 2690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp
!
    If (k2sg(isg,i)/=21) Then
      If (ptot<dp+hipr1(1)) Goto 2690
      prshu = pp(lqp(mp), 1)**2 + pp(lqp(mp), 2)**2 + pp(lqp(mp), 3)**2
      de = sqrt(pmsg(isg,i)**2+ptot**2) - sqrt(pmsg(isg,i)**2+(ptot-dp)**2)
      ershu = (pp(lqp(mp),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 2610
      pp(lqp(mp), 4) = sqrt(ershu)
      pp(lqp(mp), 5) = sqrt(amshu)
    End If
!                ********reshuffle the energy when jet has mass
!
    r0 = rdp(mp)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
!                ********momentum and energy transfer from jet
    npj(lqp(mp)) = npj(lqp(mp)) + 1
    kfpj(lqp(mp), npj(lqp(mp))) = 21
    pjpx(lqp(mp), npj(lqp(mp))) = dp1
    pjpy(lqp(mp), npj(lqp(mp))) = dp2
    pjpz(lqp(mp), npj(lqp(mp))) = dp3
    pjpe(lqp(mp), npj(lqp(mp))) = dp
    pjpm(lqp(mp), npj(lqp(mp))) = 0.0

    Goto 2660

    2640 mt = mt + 1
    drr = rdt(mt) - r0
    If (rn>=1.0-exp(-drr/hipr1(13))) Goto 2610
    dp = drr*hipr1(14)
    If (dp<=0.2) Goto 2610
    If (ptot<=0.4) Goto 2690
    If (ptot<=dp) dp = ptot - 0.2
    de = dp

    If (k2sg(isg,i)/=21) Then
      If (ptot<dp+hipr1(1)) Goto 2690
      prshu = pt(lqt(mt), 1)**2 + pt(lqt(mt), 2)**2 + pt(lqt(mt), 3)**2
      de = sqrt(pmsg(isg,i)**2+ptot**2) - sqrt(pmsg(isg,i)**2+(ptot-dp)**2)
      ershu = (pt(lqt(mt),4)+de-dp)**2
      amshu = ershu - prshu
      If (amshu<hipr1(1)*hipr1(1)) Goto 2610
      pt(lqt(mt), 4) = sqrt(ershu)
      pt(lqt(mt), 5) = sqrt(amshu)
    End If
!               ********reshuffle the energy when jet has mass

    r0 = rdt(mt)
    dp1 = dp*v1
    dp2 = dp*v2
    dp3 = dp*v3
!                ********momentum and energy transfer from jet
    ntj(lqt(mt)) = ntj(lqt(mt)) + 1
    kftj(lqt(mt), ntj(lqt(mt))) = 21
    pjtx(lqt(mt), ntj(lqt(mt))) = dp1
    pjty(lqt(mt), ntj(lqt(mt))) = dp2
    pjtz(lqt(mt), ntj(lqt(mt))) = dp3
    pjte(lqt(mt), ntj(lqt(mt))) = dp
    pjtm(lqt(mt), ntj(lqt(mt))) = 0.0

    2660 pxsg(isg, i) = (ptot-dp)*v1
    pysg(isg, i) = (ptot-dp)*v2
    pzsg(isg, i) = (ptot-dp)*v3
    pesg(isg, i) = pesg(isg, i) - de

    ptot = ptot - dp
    nq = nq + 1
    Goto 2600
  2690 End Do
  Return
End Subroutine quench

!
!
!
!
Subroutine hijfrg(jtp, ntp, ierror)
!        NTP=1, fragment proj string, NTP=2, targ string,
!       NTP=3, independent
!        strings from jets.  JTP is the line number of the string
!*******Fragment all leadng strings of proj and targ**************
!        IHNT2(1)=atomic #, IHNT2(2)=proton #(=-1 if anti-proton)  *
!******************************************************************
  Parameter (maxstr=150001)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hijdat/hidat0(10, 10), hidat(10)
!c      SAVE /HIJDAT/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
!
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
!c      SAVE /LUDAT1/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
!lin-4/11/01 soft:
  Common /anim/nevent, isoft, isflag, izpc
!c      SAVE /anim/
  Save

!bz3/12/99
!.....set up fragmentation function according to the number of collisions
!.....a wounded nucleon has suffered
!        IF (NTP .EQ. 1) THEN
!           NCOLL = NFP(JTP, 11)
!        ELSE IF (NTP .EQ. 2) THEN
!           NCOLL = NFT(JTP, 11)
!        ELSE IF (NTP .EQ. 3) THEN
!           NCOLL = (NFP(IASG(JTP,1), 11) + NFT(IASG(JTP,2), 11)) / 2
!        END IF
!        IF (NCOLL .LE. 1) THEN
!           PARJ(5) = 0.5
!        ELSE IF (NCOLL .EQ. 2) THEN
!           PARJ(5) = 0.75
!        ELSE IF (NCOLL .EQ. 3) THEN
!           PARJ(5) = 1.17
!        ELSE IF (NCOLL .EQ. 4) THEN
!           PARJ(5) = 2.0
!        ELSE IF (NCOLL .EQ. 5) THEN
!           PARJ(5) = 4.5
!        ELSE IF (NCOLL .GE. 6) THEN
!           PARJ(5) = 49.5
!        END IF
!        PARJ(5) = 0.5
!bz3/12/99 end

  ierror = 0
  Call luedit(0)
  n = 0
!                        ********initialize the document lines
  If (ntp==3) Then
    isg = jtp
    n = njsg(isg)
    Do i = 1, njsg(isg)
      k(i, 1) = k1sg(isg, i)
      k(i, 2) = k2sg(isg, i)
      p(i, 1) = pxsg(isg, i)
      p(i, 2) = pysg(isg, i)
      p(i, 3) = pzsg(isg, i)
      p(i, 4) = pesg(isg, i)
      p(i, 5) = pmsg(isg, i)
    End Do

!                IF(IHPR2(1).GT.0) CALL ATTRAD(IERROR)
!                IF(IERROR.NE.0) RETURN
!                CALL LULIST(1)
    If (isoft/=2 .Or. isflag/=0) Call luexec
    Return
  End If
!
  If (ntp==2) Goto 200
  If (jtp>ihnt2(1)) Return
  If (nfp(jtp,5)/=3 .And. nfp(jtp,3)/=0 .And. npj(jtp)==0 .And. nfp(jtp,10)==0) Goto 1000
  If (nfp(jtp,15)==-1) Then
    kf1 = nfp(jtp, 2)
    kf2 = nfp(jtp, 1)
    pq21 = pp(jtp, 6)
    pq22 = pp(jtp, 7)
    pq11 = pp(jtp, 8)
    pq12 = pp(jtp, 9)
    am1 = pp(jtp, 15)
    am2 = pp(jtp, 14)
  Else
    kf1 = nfp(jtp, 1)
    kf2 = nfp(jtp, 2)
    pq21 = pp(jtp, 8)
    pq22 = pp(jtp, 9)
    pq11 = pp(jtp, 6)
    pq12 = pp(jtp, 7)
    am1 = pp(jtp, 14)
    am2 = pp(jtp, 15)
  End If

!        ********for NFP(JTP,15)=-1 NFP(JTP,1) IS IN -Z DIRECTION
  pb1 = pq11 + pq21
  pb2 = pq12 + pq22
  pb3 = pp(jtp, 3)
  pecm = pp(jtp, 5)
  btz = pb3/pp(jtp, 4)
  If ((abs(pb1-pp(jtp,1))>0.01 .Or. abs(pb2-pp(jtp,2))>0.01) .And. ihpr2(10)/=0) Write (6, *) '  Pt of Q and QQ do not sum to the total', jtp, ntp, pq11, pq21, pb1, '*', pq12, pq22, pb2, '*', pp(jtp, 1), pp(jtp, 2)
  Goto 300

  200 If (jtp>ihnt2(3)) Return
  If (nft(jtp,5)/=3 .And. nft(jtp,3)/=0 .And. ntj(jtp)==0 .And. nft(jtp,10)==0) Goto 1200
  If (nft(jtp,15)==1) Then
    kf1 = nft(jtp, 1)
    kf2 = nft(jtp, 2)
    pq11 = pt(jtp, 6)
    pq12 = pt(jtp, 7)
    pq21 = pt(jtp, 8)
    pq22 = pt(jtp, 9)
    am1 = pt(jtp, 14)
    am2 = pt(jtp, 15)
  Else
    kf1 = nft(jtp, 2)
    kf2 = nft(jtp, 1)
    pq11 = pt(jtp, 8)
    pq12 = pt(jtp, 9)
    pq21 = pt(jtp, 6)
    pq22 = pt(jtp, 7)
    am1 = pt(jtp, 15)
    am2 = pt(jtp, 14)
  End If
!        ********for NFT(JTP,15)=1 NFT(JTP,1) IS IN +Z DIRECTION
  pb1 = pq11 + pq21
  pb2 = pq12 + pq22
  pb3 = pt(jtp, 3)
  pecm = pt(jtp, 5)
  btz = pb3/pt(jtp, 4)

  If ((abs(pb1-pt(jtp,1))>0.01 .Or. abs(pb2-pt(jtp,2))>0.01) .And. ihpr2(10)/=0) Write (6, *) '  Pt of Q and QQ do not sum to the total', jtp, ntp, pq11, pq21, pb1, '*', pq12, pq22, pb2, '*', pt(jtp, 1), pt(jtp, 2)
  300 If (pecm<hipr1(1)) Then
    ierror = 1
    If (ihpr2(10)==0) Return
    Write (6, *) ' ECM=', pecm, ' energy of the string is too small'
!lin:
    Write (6, *) 'JTP,NTP,pq=', jtp, ntp, pq11, pq12, pq21, pq22
    Return
  End If
  amt = pecm**2 + pb1**2 + pb2**2
  amt1 = am1**2 + pq11**2 + pq12**2
  amt2 = am2**2 + pq21**2 + pq22**2
  pzcm = sqrt(abs(amt**2+amt1**2+amt2**2-2.0*amt*amt1-2.0*amt*amt2-2.0*amt1*amt2))/2.0/sqrt(amt)
!                *******PZ of end-partons in c.m. frame of the string
  k(1, 1) = 2
  k(1, 2) = kf1
  p(1, 1) = pq11
  p(1, 2) = pq12
  p(1, 3) = pzcm
  p(1, 4) = sqrt(amt1+pzcm**2)
  p(1, 5) = am1
  k(2, 1) = 1
  k(2, 2) = kf2
  p(2, 1) = pq21
  p(2, 2) = pq22
  p(2, 3) = -pzcm
  p(2, 4) = sqrt(amt2+pzcm**2)
  p(2, 5) = am2
  n = 2
!*****
  Call hirobo(0.0, 0.0, 0.0, 0.0, btz)
  jetot = 0
  If ((pq21**2+pq22**2)>(pq11**2+pq12**2)) Then
    pmax1 = p(2, 1)
    pmax2 = p(2, 2)
    pmax3 = p(2, 3)
  Else
    pmax1 = p(1, 1)
    pmax2 = p(1, 2)
    pmax3 = p(1, 3)
  End If
  If (ntp==1) Then
    pp(jtp, 10) = pmax1
    pp(jtp, 11) = pmax2
    pp(jtp, 12) = pmax3
  Else If (ntp==2) Then
    pt(jtp, 10) = pmax1
    pt(jtp, 11) = pmax2
    pt(jtp, 12) = pmax3
  End If
!*******************attach produced jets to the leadng partons****
  If (ntp==1 .And. npj(jtp)/=0) Then
    jetot = npj(jtp)
!                IF(NPJ(JTP).GE.2) CALL HIJSRT(JTP,1)
!                        ********sort jets in order of y
    iex = 0
    If ((abs(kf1)>1000 .And. kf1<0) .Or. (abs(kf1)<1000 .And. kf1>0)) iex = 1
    Do i = n, 2, -1
      Do j = 1, 5
        ii = npj(jtp) + i
        k(ii, j) = k(i, j)
        p(ii, j) = p(i, j)
        v(ii, j) = v(i, j)
      End Do
    End Do

    Do i = 1, npj(jtp)
      Do j = 1, 5
        k(i+1, j) = 0
        v(i+1, j) = 0
      End Do
      i0 = i
!lin-4/12/01:                        IF(IEX.EQ.1) I0=NPJ(JTP)-I+1
      If (iex==1 .And. (isoft/=2 .Or. isflag/=0)) i0 = npj(jtp) - i + 1
!                                ********reverse the order of jets
      kk1 = kfpj(jtp, i0)
      k(i+1, 1) = 2
      k(i+1, 2) = kk1
      If (kk1/=21 .And. kk1/=0) k(i+1, 1) = 1 + (abs(kk1)+(2*iex-1)*kk1)/2/abs(kk1)
      p(i+1, 1) = pjpx(jtp, i0)
      p(i+1, 2) = pjpy(jtp, i0)
      p(i+1, 3) = pjpz(jtp, i0)
      p(i+1, 4) = pjpe(jtp, i0)
      p(i+1, 5) = pjpm(jtp, i0)
    End Do
    n = n + npj(jtp)
  Else If (ntp==2 .And. ntj(jtp)/=0) Then
    jetot = ntj(jtp)
!                IF(NTJ(JTP).GE.2)  CALL HIJSRT(JTP,2)
!                        ********sort jets in order of y
    iex = 1
    If ((abs(kf2)>1000 .And. kf2<0) .Or. (abs(kf2)<1000 .And. kf2>0)) iex = 0
    Do i = n, 2, -1
      Do j = 1, 5
        ii = ntj(jtp) + i
        k(ii, j) = k(i, j)
        p(ii, j) = p(i, j)
        v(ii, j) = v(i, j)
      End Do
    End Do
    Do i = 1, ntj(jtp)
      Do j = 1, 5
        k(i+1, j) = 0
        v(i+1, j) = 0
      End Do
      i0 = i
!lin-4/12/01:                        IF(IEX.EQ.1) I0=NTJ(JTP)-I+1
      If (iex==1 .And. (isoft/=2 .Or. isflag/=0)) i0 = ntj(jtp) - i + 1
!                                ********reverse the order of jets
      kk1 = kftj(jtp, i0)
      k(i+1, 1) = 2
      k(i+1, 2) = kk1
      If (kk1/=21 .And. kk1/=0) k(i+1, 1) = 1 + (abs(kk1)+(2*iex-1)*kk1)/2/abs(kk1)
      p(i+1, 1) = pjtx(jtp, i0)
      p(i+1, 2) = pjty(jtp, i0)
      p(i+1, 3) = pjtz(jtp, i0)
      p(i+1, 4) = pjte(jtp, i0)
      p(i+1, 5) = pjtm(jtp, i0)
    End Do
    n = n + ntj(jtp)
  End If
  If (ihpr2(1)>0 .And. ranart(nseed)<=hidat(3)) Then
    hdat20 = hidat(2)
    hpr150 = hipr1(5)
    If (ihpr2(8)==0 .And. ihpr2(3)==0 .And. ihpr2(9)==0) hidat(2) = 2.0
    If (hint1(1)>=1000.0 .And. jetot==0) Then
      hidat(2) = 3.0
      hipr1(5) = 5.0
    End If
    Call attrad(ierror)
    hidat(2) = hdat20
    hipr1(5) = hpr150
  Else If (jetot==0 .And. ihpr2(1)>0 .And. hint1(1)>=1000.0 .And. ranart(nseed)<=0.8) Then
    hdat20 = hidat(2)
    hpr150 = hipr1(5)
    hidat(2) = 3.0
    hipr1(5) = 5.0
    If (ihpr2(8)==0 .And. ihpr2(3)==0 .And. ihpr2(9)==0) hidat(2) = 2.0
    Call attrad(ierror)
    hidat(2) = hdat20
    hipr1(5) = hpr150
  End If
  If (ierror/=0) Return
!                ******** conduct soft radiations
!****************************
!
!
!lin-4/11/01 soft:
!        CALL LUEXEC
  If (isoft/=2 .Or. isflag/=0) Call luexec

  Return

  1000 n = 1
  k(1, 1) = 1
  k(1, 2) = nfp(jtp, 3)
  Do jj = 1, 5
    p(1, jj) = pp(jtp, jj)
  End Do
!                        ********proj remain as a nucleon or delta
!lin-4/11/01 soft:
!        CALL LUEXEC
  If (isoft/=2 .Or. isflag/=0) Call luexec

!        call lulist(1)
  Return
!
  1200 n = 1
  k(1, 1) = 1
  k(1, 2) = nft(jtp, 3)
  Do jj = 1, 5
    p(1, jj) = pt(jtp, jj)
  End Do
!                        ********targ remain as a nucleon or delta
!lin-4/11/01 soft:
!        CALL LUEXEC
  If (isoft/=2 .Or. isflag/=0) Call luexec

!        call lulist(1)
  Return
End Subroutine hijfrg
!
!
!
!
!****************************************************************
!        conduct soft radiation according to dipole approxiamtion
!****************************************************************
Subroutine attrad(ierror)
!
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hijdat/hidat0(10, 10), hidat(10)
!c      SAVE /HIJDAT/
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
  ierror = 0

!.....S INVARIANT MASS-SQUARED BETWEEN PARTONS I AND I+1......
!.....SM IS THE LARGEST MASS-SQUARED....

  40 sm = 0.
  jl = 1
  Do i = 1, n - 1
    s = 2.*(p(i,4)*p(i+1,4)-p(i,1)*p(i+1,1)-p(i,2)*p(i+1,2)-p(i,3)*p(i+1,3)) + p(i, 5)**2 + p(i+1, 5)**2
    If (s<0.) s = 0.
    wp = sqrt(s) - 1.5*(p(i,5)+p(i+1,5))
    If (wp>sm) Then
      pbt1 = p(i, 1) + p(i+1, 1)
      pbt2 = p(i, 2) + p(i+1, 2)
      pbt3 = p(i, 3) + p(i+1, 3)
      pbt4 = p(i, 4) + p(i+1, 4)
      btt = (pbt1**2+pbt2**2+pbt3**2)/pbt4**2
      If (btt>=1.0-1.0E-10) Goto 30
      If ((i/=1 .Or. i/=n-1) .And. (k(i,2)/=21 .And. k(i+1,2)/=21)) Goto 30
      jl = i
      sm = wp
    End If
  30 End Do
  s = (sm+1.5*(p(jl,5)+p(jl+1,5)))**2
  If (sm<hipr1(5)) Goto 2

!.....MAKE PLACE FOR ONE GLUON.....
  If (jl+1==n) Goto 190
  Do j = n, jl + 2, -1
    k(j+1, 1) = k(j, 1)
    k(j+1, 2) = k(j, 2)
    Do m = 1, 5
      p(j+1, m) = p(j, m)
    End Do
  End Do
  190 n = n + 1

!.....BOOST TO REST SYSTEM FOR PARTICLES JL AND JL+1.....
  p1 = p(jl, 1) + p(jl+1, 1)
  p2 = p(jl, 2) + p(jl+1, 2)
  p3 = p(jl, 3) + p(jl+1, 3)
  p4 = p(jl, 4) + p(jl+1, 4)
  bex = -p1/p4
  bey = -p2/p4
  bez = -p3/p4
  imin = jl
  imax = jl + 1
  Call atrobo(0., 0., bex, bey, bez, imin, imax, ierror)
  If (ierror/=0) Return
!.....ROTATE TO Z-AXIS....
  cth = p(jl, 3)/sqrt(p(jl,4)**2-p(jl,5)**2)
  If (abs(cth)>1.0) cth = max(-1., min(1.,cth))
  theta = acos(cth)
  phi = ulangl(p(jl,1), p(jl,2))
  Call atrobo(0., -phi, 0., 0., 0., imin, imax, ierror)
  Call atrobo(-theta, 0., 0., 0., 0., imin, imax, ierror)

!.....CREATE ONE GLUON AND ORIENTATE.....

  1 Call ar3jet(s, x1, x3, jl)
  Call arorie(s, x1, x3, jl)
  If (hidat(2)>0.0) Then
    ptg1 = sqrt(p(jl,1)**2+p(jl,2)**2)
    ptg2 = sqrt(p(jl+1,1)**2+p(jl+1,2)**2)
    ptg3 = sqrt(p(jl+2,1)**2+p(jl+2,2)**2)
    ptg = max(ptg1, ptg2, ptg3)
    If (ptg>hidat(2)) Then
      fmfact = exp(-(ptg**2-hidat(2)**2)/hipr1(2)**2)
      If (ranart(nseed)>fmfact) Goto 1
    End If
  End If
!.....ROTATE AND BOOST BACK.....
  imin = jl
  imax = jl + 2
  Call atrobo(theta, phi, -bex, -bey, -bez, imin, imax, ierror)
  If (ierror/=0) Return
!.....ENUMERATE THE GLUONS.....
  k(jl+2, 1) = k(jl+1, 1)
  k(jl+2, 2) = k(jl+1, 2)
  k(jl+2, 3) = k(jl+1, 3)
  k(jl+2, 4) = k(jl+1, 4)
  k(jl+2, 5) = k(jl+1, 5)
  p(jl+2, 5) = p(jl+1, 5)
  k(jl+1, 1) = 2
  k(jl+1, 2) = 21
  k(jl+1, 3) = 0
  k(jl+1, 4) = 0
  k(jl+1, 5) = 0
  p(jl+1, 5) = 0.
!----THETA FUNCTION DAMPING OF THE EMITTED GLUONS. FOR HADRON-HADRON.
!----R0=VFR(2)
!              IF(VFR(2).GT.0.) THEN
!              PTG=SQRT(P(JL+1,1)**2+P(JL+1,2)**2)
!              PTGMAX=WSTRI/2.
!              DOPT=SQRT((4.*PAR(71)*VFR(2))/WSTRI)
!              PTOPT=(DOPT*WSTRI)/(2.*VFR(2))
!              IF(PTG.GT.PTOPT) IORDER=IORDER-1
!              IF(PTG.GT.PTOPT) GOTO 1
!              ENDIF
!-----
  If (sm>=hipr1(5)) Goto 40

  2 k(1, 1) = 2
  k(1, 3) = 0
  k(1, 4) = 0
  k(1, 5) = 0
  k(n, 1) = 1
  k(n, 3) = 0
  k(n, 4) = 0
  k(n, 5) = 0

  Return
End Subroutine attrad


Subroutine ar3jet(s, x1, x3, jl)
!
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
  c = 1./3.
  If (k(jl,2)/=21 .And. k(jl+1,2)/=21) c = 8./27.
  exp1 = 3
  exp3 = 3
  If (k(jl,2)/=21) exp1 = 2
  If (k(jl+1,2)/=21) exp3 = 2
  a = 0.24**2/s
  yma = alog(.5/sqrt(a)+sqrt(.25/a-1))
  d = 4.*c*yma
  sm1 = p(jl, 5)**2/s
  sm3 = p(jl+1, 5)**2/s
  xt2m = (1.-2.*sqrt(sm1)+sm1-sm3)*(1.-2.*sqrt(sm3)-sm1+sm3)
  xt2m = min(.25, xt2m)
  ntry = 0
  1 If (ntry==5000) Then
    x1 = .5*(2.*sqrt(sm1)+1.+sm1-sm3)
    x3 = .5*(2.*sqrt(sm3)+1.-sm1+sm3)
    Return
  End If
  ntry = ntry + 1

  xt2 = a*(xt2m/a)**(ranart(nseed)**(1./d))

  ymax = alog(.5/sqrt(xt2)+sqrt(.25/xt2-1.))
  y = (2.*ranart(nseed)-1.)*ymax
  x1 = 1. - sqrt(xt2)*exp(y)
  x3 = 1. - sqrt(xt2)*exp(-y)
  x2 = 2. - x1 - x3
  neg = 0
  If (k(jl,2)/=21 .Or. k(jl+1,2)/=21) Then
    If ((1.-x1)*(1.-x2)*(1.-x3)-x2*sm1*(1.-x1)-x2*sm3*(1.-x3)<=0. .Or. x1<=2.*sqrt(sm1)-sm1+sm3 .Or. x3<=2.*sqrt(sm3)-sm3+sm1) neg = 1
    x1 = x1 + sm1 - sm3
    x3 = x3 - sm1 + sm3
  End If
  If (neg==1) Goto 1

  fg = 2.*ymax*c*(x1**exp1+x3**exp3)/d
  xt2m = xt2
  If (fg<ranart(nseed)) Goto 1

  Return
End Subroutine ar3jet
!*************************************************************


Subroutine arorie(s, x1, x3, jl)
!
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
  w = sqrt(s)
  x2 = 2. - x1 - x3
  e1 = .5*x1*w
  e3 = .5*x3*w
  p1 = sqrt(e1**2-p(jl,5)**2)
  p3 = sqrt(e3**2-p(jl+1,5)**2)
  cbet = 1.
  If (p1>0. .And. p3>0.) cbet = (p(jl,5)**2+p(jl+1,5)**2+2.*e1*e3-s*(1.-x2))/(2.*p1*p3)
  If (abs(cbet)>1.0) cbet = max(-1., min(1.,cbet))
  bet = acos(cbet)

!.....MINIMIZE PT1-SQUARED PLUS PT3-SQUARED.....
  If (p1>=p3) Then
    psi = .5*ulangl(p1**2+p3**2*cos(2.*bet), -p3**2*sin(2.*bet))
    pt1 = p1*sin(psi)
    pz1 = p1*cos(psi)
    pt3 = p3*sin(psi+bet)
    pz3 = p3*cos(psi+bet)
  Else If (p3>p1) Then
    psi = .5*ulangl(p3**2+p1**2*cos(2.*bet), -p1**2*sin(2.*bet))
    pt1 = p1*sin(bet+psi)
    pz1 = -p1*cos(bet+psi)
    pt3 = p3*sin(psi)
    pz3 = -p3*cos(psi)
  End If

  del = 2.0*hipr1(40)*ranart(nseed)
  p(jl, 4) = e1
  p(jl, 1) = pt1*sin(del)
  p(jl, 2) = -pt1*cos(del)
  p(jl, 3) = pz1
  p(jl+2, 4) = e3
  p(jl+2, 1) = pt3*sin(del)
  p(jl+2, 2) = -pt3*cos(del)
  p(jl+2, 3) = pz3
  p(jl+1, 4) = w - e1 - e3
  p(jl+1, 1) = -p(jl, 1) - p(jl+2, 1)
  p(jl+1, 2) = -p(jl, 2) - p(jl+2, 2)
  p(jl+1, 3) = -p(jl, 3) - p(jl+2, 3)
  Return
End Subroutine arorie


!
!*******************************************************************
!        make  boost and rotation to entries from IMIN to IMAX
!*******************************************************************
Subroutine atrobo(the, phi, bex, bey, bez, imin, imax, ierror)
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Dimension rot(3, 3), pv(3)
  Double Precision dp(4), dbex, dbey, dbez, dga, dga2, dbep, dgabep
  Save
  ierror = 0

  If (imin<=0 .Or. imax>n .Or. imin>imax) Return

  If (the**2+phi**2>1E-20) Then
!...ROTATE (TYPICALLY FROM Z AXIS TO DIRECTION THETA,PHI)
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
!**************           IF(MOD(K(I,1)/10000,10).GE.6) GOTO 120
      Do j = 1, 3
        pv(j) = p(i, j)
      End Do
      Do j = 1, 3
        p(i, j) = rot(j, 1)*pv(1) + rot(j, 2)*pv(2) + rot(j, 3)*pv(3)
      End Do
    End Do
  End If

  If (bex**2+bey**2+bez**2>1E-20) Then
!...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
    dbex = dble(bex)
    dbey = dble(bey)
    dbez = dble(bez)
    dga2 = 1D0 - dbex**2 - dbey**2 - dbez**2
    If (dga2<=0D0) Then
      ierror = 1
      Return
    End If
    dga = 1D0/dsqrt(dga2)
    Do i = imin, imax
!*************           IF(MOD(K(I,1)/10000,10).GE.6) GOTO 140
      Do j = 1, 4
        dp(j) = dble(p(i,j))
      End Do
      dbep = dbex*dp(1) + dbey*dp(2) + dbez*dp(3)
      dgabep = dga*(dga*dbep/(1D0+dga)+dp(4))
      p(i, 1) = sngl(dp(1)+dgabep*dbex)
      p(i, 2) = sngl(dp(2)+dgabep*dbey)
      p(i, 3) = sngl(dp(3)+dgabep*dbez)
      p(i, 4) = sngl(dga*(dp(4)+dbep))
    End Do
  End If

  Return
End Subroutine atrobo
!
!
!
Subroutine hijhrd(jp, jt, jout, jflg, iopjt0)
!
!        IOPTJET=1, ALL JET WILL FORM SINGLE STRING SYSTEM
!                0, ONLY Q-QBAR JET FORM SINGLE STRING SYSTEM
!*******Perform jets production and fragmentation when JP JT *******
!     scatter. JOUT-> number of hard scatterings precede this one  *
!     for the the same pair(JP,JT). JFLG->a flag to show whether   *
!     jets can be produced (with valence quark=1,gluon=2, q-qbar=3)*
!     or not(0). Information of jets are in  COMMON/ATTJET and     *
!     /MINJET. ABS(NFP(JP,6)) is the total number jets produced by *
!    JP. If NFP(JP,6)<0 JP can not produce jet anymore.                   *
!*******************************************************************
  Parameter (maxstr=150001)
  Dimension ip(100, 2), ipq(50), ipb(50), it(100, 2), itq(50), itb(50)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hijdat/hidat0(10, 10), hidat(10)
!c      SAVE /HIJDAT/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
!        COMMON/HJJET4/NDR,IADR(900,2),KFDR(900),PDR(900,5)
  Common /hjjet4/ndr, iadr(maxstr, 2), kfdr(maxstr), pdr(maxstr, 5)
  Common /xydr/rtdr(maxstr, 2)
!c      SAVE /HJJET4/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
!************************************ HIJING common block
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
!c      SAVE /LUJETS/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
!c      SAVE /LUDAT1/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
!c      SAVE /PYSUBS/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
!c      SAVE /PYPARS/
  Common /pyint1/mint(400), vint(400)
!c      SAVE /PYINT1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
!c      SAVE /PYINT2/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
!c      SAVE /PYINT5/
  Common /hpint/mint4, mint5, atco(200, 20), atxs(0:200)
!c      SAVE /HPINT/
!lin-2/2012 correction:
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Save
!*********************************** LU common block
  mxjt = 500
!                SIZE OF COMMON BLOCK FOR # OF PARTON PER STRING
  mxsg = 900
!                SIZE OF COMMON BLOCK FOR # OF SINGLE STRINGS
  mxsj = 100
!                SIZE OF COMMON BLOCK FOR # OF PARTON PER SINGLE
!                STRING
  jflg = 0
  ihnt2(11) = jp
  ihnt2(12) = jt
!
  iopjet = iopjt0
  If (iopjet==1 .And. (nfp(jp,6)/=0 .Or. nft(jt,6)/=0)) iopjet = 0
  If (jp>ihnt2(1) .Or. jt>ihnt2(3)) Return
  If (nfp(jp,6)<0 .Or. nft(jt,6)<0) Return
!                ******** JP or JT can not produce jet anymore
!
  If (jout==0) Then
    epp = pp(jp, 4) + pp(jp, 3)
    epm = pp(jp, 4) - pp(jp, 3)
    etp = pt(jt, 4) + pt(jt, 3)
    etm = pt(jt, 4) - pt(jt, 3)
    If (epp<0.0) Goto 1000
    If (epm<0.0) Goto 1000
    If (etp<0.0) Goto 1000
    If (etm<0.0) Goto 1000
    If (epp/(epm+0.01)<=etp/(etm+0.01)) Return
  End If
!                ********for the first hard scattering of (JP,JT)
!                        have collision only when Ycm(JP)>Ycm(JT)

  ecut1 = hipr1(1) + hipr1(8) + pp(jp, 14) + pp(jp, 15)
  ecut2 = hipr1(1) + hipr1(8) + pt(jt, 14) + pt(jt, 15)
  If (pp(jp,4)<=ecut1) Then
    nfp(jp, 6) = -abs(nfp(jp,6))
    Return
  End If
  If (pt(jt,4)<=ecut2) Then
    nft(jt, 6) = -abs(nft(jt,6))
    Return
  End If
!                *********must have enough energy to produce jets

  miss = 0
  misp = 0
  mist = 0
!
  If (nfp(jp,10)==0 .And. nft(jt,10)==0) Then
    mint(44) = mint4
    mint(45) = mint5
    xsec(0, 1) = atxs(0)
    xsec(11, 1) = atxs(11)
    xsec(12, 1) = atxs(12)
    xsec(28, 1) = atxs(28)
    Do i = 1, 20
      coef(11, i) = atco(11, i)
      coef(12, i) = atco(12, i)
      coef(28, i) = atco(28, i)
    End Do
  Else
    isub11 = 0
    isub12 = 0
    isub28 = 0
    If (xsec(11,1)/=0) isub11 = 1
    If (xsec(12,1)/=0) isub12 = 1
    If (xsec(28,1)/=0) isub28 = 1
    mint(44) = mint4 - isub11 - isub12 - isub28
    mint(45) = mint5 - isub11 - isub12 - isub28
    xsec(0, 1) = atxs(0) - atxs(11) - atxs(12) - atxs(28)
    xsec(11, 1) = 0.0
    xsec(12, 1) = 0.0
    xsec(28, 1) = 0.0
    Do i = 1, 20
      coef(11, i) = 0.0
      coef(12, i) = 0.0
      coef(28, i) = 0.0
    End Do
  End If
!        ********Scatter the valence quarks only once per NN
!       collision,
!                afterwards only gluon can have hard scattering.
  155 Call pythia
  jj = mint(31)
  If (jj/=1) Goto 155
!                *********one hard collision at a time
  If (k(7,2)==-k(8,2)) Then
    qmass2 = (p(7,4)+p(8,4))**2 - (p(7,1)+p(8,1))**2 - (p(7,2)+p(8,2))**2 - (p(7,3)+p(8,3))**2
    qm = ulmass(k(7,2))
    If (qmass2<(2.0*qm+hipr1(1))**2) Goto 155
  End If
!                ********q-qbar jets must has minimum mass HIPR1(1)
  pxp = pp(jp, 1) - p(3, 1)
  pyp = pp(jp, 2) - p(3, 2)
  pzp = pp(jp, 3) - p(3, 3)
  pep = pp(jp, 4) - p(3, 4)
  pxt = pt(jt, 1) - p(4, 1)
  pyt = pt(jt, 2) - p(4, 2)
  pzt = pt(jt, 3) - p(4, 3)
  pet = pt(jt, 4) - p(4, 4)

  If (pep<=ecut1) Then
    misp = misp + 1
    If (misp<50) Goto 155
    nfp(jp, 6) = -abs(nfp(jp,6))
    Return
  End If
  If (pet<=ecut2) Then
    mist = mist + 1
    If (mist<50) Goto 155
    nft(jt, 6) = -abs(nft(jt,6))
    Return
  End If
!                ******** if the remain energy<ECUT the proj or targ
!                         can not produce jet anymore

  wp = pep + pzp + pet + pzt
  wm = pep - pzp + pet - pzt
  If (wp<0.0 .Or. wm<0.0) Then
    miss = miss + 1
!lin-6/2009 Let user set the limit when selecting high-Pt events
!     because more attempts may be needed:
!                IF(MISS.LT.50) GO TO 155
    If (pttrig>0) Then
      If (miss<maxmiss) Then
        Write (6, *) 'Failed to generate minijet Pt>', pttrig, 'GeV'
        Goto 155
      End If
    Else
      If (miss<50) Goto 155
    End If

    Return
  End If
!                ********the total W+, W- must be positive
  sw = wp*wm
  ampx = sqrt((ecut1-hipr1(8))**2+pxp**2+pyp**2+0.01)
  amtx = sqrt((ecut2-hipr1(8))**2+pxt**2+pyt**2+0.01)
  sxx = (ampx+amtx)**2
  If (sw<sxx .Or. vint(43)<hipr1(1)) Then
    miss = miss + 1
!lin-6/2009
!                IF(MISS.LT.50) GO TO 155
    If (miss>maxmiss) Goto 155
    Return
  End If
!                ********the proj and targ remnants must have at least
!                        a CM energy that can produce two strings
!                        with minimum mass HIPR1(1)(see HIJSFT HIJFRG)
!
  hint1(41) = p(7, 1)
  hint1(42) = p(7, 2)
  hint1(43) = p(7, 3)
  hint1(44) = p(7, 4)
  hint1(45) = p(7, 5)
  hint1(46) = sqrt(p(7,1)**2+p(7,2)**2)
  hint1(51) = p(8, 1)
  hint1(52) = p(8, 2)
  hint1(53) = p(8, 3)
  hint1(54) = p(8, 4)
  hint1(55) = p(8, 5)
  hint1(56) = sqrt(p(8,1)**2+p(8,2)**2)
  ihnt2(14) = k(7, 2)
  ihnt2(15) = k(8, 2)
!
  pinird = (1.0-exp(-2.0*(vint(47)-hidat(1))))/(1.0+exp(-2.0*(vint(47)-hidat(1))))
  iinird = 0
  If (ranart(nseed)<=pinird) iinird = 1
  If (k(7,2)==-k(8,2)) Goto 190
  If (k(7,2)==21 .And. k(8,2)==21 .And. iopjet==1) Goto 190
!*******************************************************************
!        gluon  jets are going to be connectd with
!        the final leadng string of quark-aintquark
!*******************************************************************
  jflg = 2
  jpp = 0
  lpq = 0
  lpb = 0
  jtt = 0
  ltq = 0
  ltb = 0
  is7 = 0
  is8 = 0
  hint1(47) = 0.0
  hint1(48) = 0.0
  hint1(49) = 0.0
  hint1(50) = 0.0
  hint1(67) = 0.0
  hint1(68) = 0.0
  hint1(69) = 0.0
  hint1(70) = 0.0
  Do i = 9, n
    If (k(i,3)==1 .Or. k(i,3)==2 .Or. abs(k(i,2))>30) Goto 180
!************************************************************
    If (k(i,3)==7) Then
      hint1(47) = hint1(47) + p(i, 1)
      hint1(48) = hint1(48) + p(i, 2)
      hint1(49) = hint1(49) + p(i, 3)
      hint1(50) = hint1(50) + p(i, 4)
    End If
    If (k(i,3)==8) Then
      hint1(67) = hint1(67) + p(i, 1)
      hint1(68) = hint1(68) + p(i, 2)
      hint1(69) = hint1(69) + p(i, 3)
      hint1(70) = hint1(70) + p(i, 4)
    End If
!************************modifcation made on Apr 10. 1996*****
    If (k(i,2)>21 .And. k(i,2)<=30) Then
      ndr = ndr + 1
      iadr(ndr, 1) = jp
      iadr(ndr, 2) = jt
      kfdr(ndr) = k(i, 2)
      pdr(ndr, 1) = p(i, 1)
      pdr(ndr, 2) = p(i, 2)
      pdr(ndr, 3) = p(i, 3)
      pdr(ndr, 4) = p(i, 4)
      pdr(ndr, 5) = p(i, 5)
      rtdr(ndr, 1) = 0.5*(yp(1,jp)+yt(1,jt))
      rtdr(ndr, 2) = 0.5*(yp(2,jp)+yt(2,jt))
!************************************************************
      Goto 180
!************************correction made on Oct. 14,1994*****
    End If
    If (k(i,3)==7 .Or. k(i,3)==3) Then
      If (k(i,3)==7 .And. k(i,2)/=21 .And. k(i,2)==k(7,2) .And. is7==0) Then
        pp(jp, 10) = p(i, 1)
        pp(jp, 11) = p(i, 2)
        pp(jp, 12) = p(i, 3)
        pzp = pzp + p(i, 3)
        pep = pep + p(i, 4)
        nfp(jp, 10) = 1
        is7 = 1
        Goto 180
      End If
      If (k(i,3)==3 .And. (k(i,2)/=21 .Or. iinird==0)) Then
        pxp = pxp + p(i, 1)
        pyp = pyp + p(i, 2)
        pzp = pzp + p(i, 3)
        pep = pep + p(i, 4)
        Goto 180
      End If
      jpp = jpp + 1
      ip(jpp, 1) = i
      ip(jpp, 2) = 0
      If (k(i,2)/=21) Then
        If (k(i,2)>0) Then
          lpq = lpq + 1
          ipq(lpq) = jpp
          ip(jpp, 2) = lpq
        Else If (k(i,2)<0) Then
          lpb = lpb + 1
          ipb(lpb) = jpp
          ip(jpp, 2) = -lpb
        End If
      End If
    Else If (k(i,3)==8 .Or. k(i,3)==4) Then
      If (k(i,3)==8 .And. k(i,2)/=21 .And. k(i,2)==k(8,2) .And. is8==0) Then
        pt(jt, 10) = p(i, 1)
        pt(jt, 11) = p(i, 2)
        pt(jt, 12) = p(i, 3)
        pzt = pzt + p(i, 3)
        pet = pet + p(i, 4)
        nft(jt, 10) = 1
        is8 = 1
        Goto 180
      End If
      If (k(i,3)==4 .And. (k(i,2)/=21 .Or. iinird==0)) Then
        pxt = pxt + p(i, 1)
        pyt = pyt + p(i, 2)
        pzt = pzt + p(i, 3)
        pet = pet + p(i, 4)
        Goto 180
      End If
      jtt = jtt + 1
      it(jtt, 1) = i
      it(jtt, 2) = 0
      If (k(i,2)/=21) Then
        If (k(i,2)>0) Then
          ltq = ltq + 1
          itq(ltq) = jtt
          it(jtt, 2) = ltq
        Else If (k(i,2)<0) Then
          ltb = ltb + 1
          itb(ltb) = jtt
          it(jtt, 2) = -ltb
        End If
      End If
    End If
  180 End Do
!
!
  If (lpq/=lpb .Or. ltq/=ltb) Then
    miss = miss + 1
!lin-6/2009
!                IF(MISS.LE.50) GO TO 155
    If (miss<=maxmiss) Goto 155
    Write (6, *) ' Q -QBAR NOT MATCHED IN HIJHRD'
    jflg = 0
    Return
  End If
!****The following will rearrange the partons so that a quark is***
!****allways followed by an anti-quark ****************************

  j = 0
  181 j = j + 1
  If (j>jpp) Goto 182
  If (ip(j,2)==0) Then
    Goto 181
  Else If (ip(j,2)/=0) Then
    lp = abs(ip(j,2))
    ip1 = ip(j, 1)
    ip2 = ip(j, 2)
    ip(j, 1) = ip(ipq(lp), 1)
    ip(j, 2) = ip(ipq(lp), 2)
    ip(ipq(lp), 1) = ip1
    ip(ipq(lp), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipq(lp)
    Else If (ip2<0) Then
      ipb(-ip2) = ipq(lp)
    End If
!                ********replace J with a quark
    ip1 = ip(j+1, 1)
    ip2 = ip(j+1, 2)
    ip(j+1, 1) = ip(ipb(lp), 1)
    ip(j+1, 2) = ip(ipb(lp), 2)
    ip(ipb(lp), 1) = ip1
    ip(ipb(lp), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipb(lp)
    Else If (ip2<0) Then
      ipb(-ip2) = ipb(lp)
    End If
!                ******** replace J+1 with anti-quark
    j = j + 1
    Goto 181
  End If

  182 j = 0
  183 j = j + 1
  If (j>jtt) Goto 184
  If (it(j,2)==0) Then
    Goto 183
  Else If (it(j,2)/=0) Then
    lt = abs(it(j,2))
    it1 = it(j, 1)
    it2 = it(j, 2)
    it(j, 1) = it(itq(lt), 1)
    it(j, 2) = it(itq(lt), 2)
    it(itq(lt), 1) = it1
    it(itq(lt), 2) = it2
    If (it2>0) Then
      itq(it2) = itq(lt)
    Else If (it2<0) Then
      itb(-it2) = itq(lt)
    End If
!                ********replace J with a quark
    it1 = it(j+1, 1)
    it2 = it(j+1, 2)
    it(j+1, 1) = it(itb(lt), 1)
    it(j+1, 2) = it(itb(lt), 2)
    it(itb(lt), 1) = it1
    it(itb(lt), 2) = it2
    If (it2>0) Then
      itq(it2) = itb(lt)
    Else If (it2<0) Then
      itb(-it2) = itb(lt)
    End If
!                ******** replace J+1 with anti-quark
    j = j + 1
    Goto 183

  End If

  184 Continue
  If (npj(jp)+jpp>mxjt .Or. ntj(jt)+jtt>mxjt) Then
    jflg = 0
    Write (6, *) 'number of partons per string exceeds'
    Write (6, *) 'the common block size'
    Return
  End If
!                        ********check the bounds of common blocks
  Do j = 1, jpp
    kfpj(jp, npj(jp)+j) = k(ip(j,1), 2)
    pjpx(jp, npj(jp)+j) = p(ip(j,1), 1)
    pjpy(jp, npj(jp)+j) = p(ip(j,1), 2)
    pjpz(jp, npj(jp)+j) = p(ip(j,1), 3)
    pjpe(jp, npj(jp)+j) = p(ip(j,1), 4)
    pjpm(jp, npj(jp)+j) = p(ip(j,1), 5)
  End Do
  npj(jp) = npj(jp) + jpp
  Do j = 1, jtt
    kftj(jt, ntj(jt)+j) = k(it(j,1), 2)
    pjtx(jt, ntj(jt)+j) = p(it(j,1), 1)
    pjty(jt, ntj(jt)+j) = p(it(j,1), 2)
    pjtz(jt, ntj(jt)+j) = p(it(j,1), 3)
    pjte(jt, ntj(jt)+j) = p(it(j,1), 4)
    pjtm(jt, ntj(jt)+j) = p(it(j,1), 5)
  End Do
  ntj(jt) = ntj(jt) + jtt
  Goto 900
!*****************************************************************
!This is the case of a quark-antiquark jet it will fragment alone
!****************************************************************
  190 jflg = 3
  If (k(7,2)/=21 .And. k(8,2)/=21 .And. k(7,2)*k(8,2)>0) Goto 155
  jpp = 0
  lpq = 0
  lpb = 0
  Do i = 9, n
    If (k(i,3)==1 .Or. k(i,3)==2 .Or. abs(k(i,2))>30) Goto 200
    If (k(i,2)>21 .And. k(i,2)<=30) Then
      ndr = ndr + 1
      iadr(ndr, 1) = jp
      iadr(ndr, 2) = jt
      kfdr(ndr) = k(i, 2)
      pdr(ndr, 1) = p(i, 1)
      pdr(ndr, 2) = p(i, 2)
      pdr(ndr, 3) = p(i, 3)
      pdr(ndr, 4) = p(i, 4)
      pdr(ndr, 5) = p(i, 5)
      rtdr(ndr, 1) = 0.5*(yp(1,jp)+yt(1,jt))
      rtdr(ndr, 2) = 0.5*(yp(2,jp)+yt(2,jt))
!************************************************************
      Goto 200
!************************correction made on Oct. 14,1994*****
    End If
    If (k(i,3)==3 .And. (k(i,2)/=21 .Or. iinird==0)) Then
      pxp = pxp + p(i, 1)
      pyp = pyp + p(i, 2)
      pzp = pzp + p(i, 3)
      pep = pep + p(i, 4)
      Goto 200
    End If
    If (k(i,3)==4 .And. (k(i,2)/=21 .Or. iinird==0)) Then
      pxt = pxt + p(i, 1)
      pyt = pyt + p(i, 2)
      pzt = pzt + p(i, 3)
      pet = pet + p(i, 4)
      Goto 200
    End If
    jpp = jpp + 1
    ip(jpp, 1) = i
    ip(jpp, 2) = 0
    If (k(i,2)/=21) Then
      If (k(i,2)>0) Then
        lpq = lpq + 1
        ipq(lpq) = jpp
        ip(jpp, 2) = lpq
      Else If (k(i,2)<0) Then
        lpb = lpb + 1
        ipb(lpb) = jpp
        ip(jpp, 2) = -lpb
      End If
    End If
  200 End Do
  If (lpq/=lpb) Then
    miss = miss + 1
!lin-6/2009
!           IF(MISS.LE.50) GO TO 155
    If (miss<=maxmiss) Goto 155
    Write (6, *) lpq, lpb, 'Q-QBAR NOT CONSERVED OR NOT MATCHED'
    jflg = 0
    Return
  End If

!**** The following will rearrange the partons so that a quark is***
!**** allways followed by an anti-quark ****************************
  j = 0
  220 j = j + 1
  If (j>jpp) Goto 222
  If (ip(j,2)==0) Goto 220
  lp = abs(ip(j,2))
  ip1 = ip(j, 1)
  ip2 = ip(j, 2)
  ip(j, 1) = ip(ipq(lp), 1)
  ip(j, 2) = ip(ipq(lp), 2)
  ip(ipq(lp), 1) = ip1
  ip(ipq(lp), 2) = ip2
  If (ip2>0) Then
    ipq(ip2) = ipq(lp)
  Else If (ip2<0) Then
    ipb(-ip2) = ipq(lp)
  End If
  ipq(lp) = j
!                ********replace J with a quark
  ip1 = ip(j+1, 1)
  ip2 = ip(j+1, 2)
  ip(j+1, 1) = ip(ipb(lp), 1)
  ip(j+1, 2) = ip(ipb(lp), 2)
  ip(ipb(lp), 1) = ip1
  ip(ipb(lp), 2) = ip2
  If (ip2>0) Then
    ipq(ip2) = ipb(lp)
  Else If (ip2<0) Then
    ipb(-ip2) = ipb(lp)
  End If
!                ******** replace J+1 with an anti-quark
  ipb(lp) = j + 1
  j = j + 1
  Goto 220

  222 Continue
  If (lpq>=1) Then
    Do l0 = 2, lpq
      ip1 = ip(2*l0-3, 1)
      ip2 = ip(2*l0-3, 2)
      ip(2*l0-3, 1) = ip(ipq(l0), 1)
      ip(2*l0-3, 2) = ip(ipq(l0), 2)
      ip(ipq(l0), 1) = ip1
      ip(ipq(l0), 2) = ip2
      If (ip2>0) Then
        ipq(ip2) = ipq(l0)
      Else If (ip2<0) Then
        ipb(-ip2) = ipq(l0)
      End If
      ipq(l0) = 2*l0 - 3
!
      ip1 = ip(2*l0-2, 1)
      ip2 = ip(2*l0-2, 2)
      ip(2*l0-2, 1) = ip(ipb(l0), 1)
      ip(2*l0-2, 2) = ip(ipb(l0), 2)
      ip(ipb(l0), 1) = ip1
      ip(ipb(l0), 2) = ip2
      If (ip2>0) Then
        ipq(ip2) = ipb(l0)
      Else If (ip2<0) Then
        ipb(-ip2) = ipb(l0)
      End If
      ipb(l0) = 2*l0 - 2
    End Do
!                ********move all the qqbar pair to the front of
!                                the list, except the first pair
    ip1 = ip(2*lpq-1, 1)
    ip2 = ip(2*lpq-1, 2)
    ip(2*lpq-1, 1) = ip(ipq(1), 1)
    ip(2*lpq-1, 2) = ip(ipq(1), 2)
    ip(ipq(1), 1) = ip1
    ip(ipq(1), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipq(1)
    Else If (ip2<0) Then
      ipb(-ip2) = ipq(1)
    End If
    ipq(1) = 2*lpq - 1
!                ********move the first quark to the beginning of
!                                the last string system
    ip1 = ip(jpp, 1)
    ip2 = ip(jpp, 2)
    ip(jpp, 1) = ip(ipb(1), 1)
    ip(jpp, 2) = ip(ipb(1), 2)
    ip(ipb(1), 1) = ip1
    ip(ipb(1), 2) = ip2
    If (ip2>0) Then
      ipq(ip2) = ipb(1)
    Else If (ip2<0) Then
      ipb(-ip2) = ipb(1)
    End If
    ipb(1) = jpp
!                ********move the first anti-quark to the end of the
!                        last string system
  End If
  If (nsg>=mxsg) Then
    jflg = 0
    Write (6, *) 'number of jets forming single strings exceeds'
    Write (6, *) 'the common block size'
    Return
  End If
  If (jpp>mxsj) Then
    jflg = 0
    Write (6, *) 'number of partons per single jet system'
    Write (6, *) 'exceeds the common block size'
    Return
  End If
!                ********check the bounds of common block size
  nsg = nsg + 1
  njsg(nsg) = jpp
  iasg(nsg, 1) = jp
  iasg(nsg, 2) = jt
  iasg(nsg, 3) = 0
  Do i = 1, jpp
    k1sg(nsg, i) = 2
    k2sg(nsg, i) = k(ip(i,1), 2)
    If (k2sg(nsg,i)<0) k1sg(nsg, i) = 1
    pxsg(nsg, i) = p(ip(i,1), 1)
    pysg(nsg, i) = p(ip(i,1), 2)
    pzsg(nsg, i) = p(ip(i,1), 3)
    pesg(nsg, i) = p(ip(i,1), 4)
    pmsg(nsg, i) = p(ip(i,1), 5)
  End Do
  k1sg(nsg, 1) = 2
  k1sg(nsg, jpp) = 1
!******* reset the energy-momentum of incoming particles ********
  900 pp(jp, 1) = pxp
  pp(jp, 2) = pyp
  pp(jp, 3) = pzp
  pp(jp, 4) = pep
  pp(jp, 5) = 0.0
  pt(jt, 1) = pxt
  pt(jt, 2) = pyt
  pt(jt, 3) = pzt
  pt(jt, 4) = pet
  pt(jt, 5) = 0.0

  nfp(jp, 6) = nfp(jp, 6) + 1
  nft(jt, 6) = nft(jt, 6) + 1
  Return
!
  1000 jflg = -1
  If (ihpr2(10)==0) Return
  Write (6, *) 'Fatal HIJHRD error'
  Write (6, *) jp, ' proj E+,E-', epp, epm, ' status', nfp(jp, 5)
  Write (6, *) jt, ' targ E+,E_', etp, etm, ' status', nft(jt, 5)
  Return
End Subroutine hijhrd
!
!


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!
!
!
Subroutine jetini(jp, jt, itrig)
!*******Initialize PYTHIA for jet production**********************
!        itrig=0: for normal processes
!        itrig=1: for triggered processes
!       JP: sequence number of the projectile
!       JT: sequence number of the target
!     For A+A collisions, one has to initilize pythia
!     separately for each type of collisions, pp, pn,np and nn,
!     or hp and hn for hA collisions. In this subroutine we use the following
!     catalogue for different type of collisions:
!     h+h: h+h (itype=1)
!     h+A: h+p (itype=1), h+n (itype=2)
!     A+h: p+h (itype=1), n+h (itype=2)
!     A+A: p+p (itype=1), p+n (itype=2), n+p (itype=3), n+n (itype=4)
!*****************************************************************
  Character beam*16, targ*16
  Dimension xsec0(8, 0:200), coef0(8, 200, 20), ini(8), mint44(8), mint45(8)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /hpint/mint4, mint5, atco(200, 20), atxs(0:200)
!c      SAVE /HPINT/
!
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
!c      SAVE /LUDAT1/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
!c      SAVE /LUDAT3/
  Common /pysubs/msel, msub(200), kfin(2, -40:40), ckin(200)
!c      SAVE /PYSUBS/
  Common /pypars/mstp(200), parp(200), msti(200), pari(200)
!c      SAVE /PYPARS/
  Common /pyint1/mint(400), vint(400)
!c      SAVE /PYINT1/
  Common /pyint2/iset(200), kfpr(200, 2), coef(200, 20), icol(40, 4, 2)
!c      SAVE /PYINT2/
  Common /pyint5/ngen(0:200, 3), xsec(0:200, 3)
!c      SAVE /PYINT5/
  Save
!lin        DATA INI/8*0/ilast/-1/
  Data ini/8*0/, ilast/ -1/
!
  ihnt2(11) = jp
  ihnt2(12) = jt
  If (ihnt2(5)/=0 .And. ihnt2(6)/=0) Then
    itype = 1
  Else If (ihnt2(5)/=0 .And. ihnt2(6)==0) Then
    itype = 1
    If (nft(jt,4)==2112) itype = 2
  Else If (ihnt2(5)==0 .And. ihnt2(6)/=0) Then
    itype = 1
    If (nfp(jp,4)==2112) itype = 2
  Else
    If (nfp(jp,4)==2212 .And. nft(jt,4)==2212) Then
      itype = 1
    Else If (nfp(jp,4)==2212 .And. nft(jt,4)==2112) Then
      itype = 2
    Else If (nfp(jp,4)==2112 .And. nft(jt,4)==2212) Then
      itype = 3
    Else
      itype = 4
    End If
  End If

!lin-12/2012 correct NN differential cross section in HIJING:
!        write(94,*) 'In JETINI: ',jp,jt,NFP(JP,4),NFT(JT,4),itype

!
  If (itrig/=0) Goto 160
  If (itrig==ilast) Goto 150
  mstp(2) = 2
!                        ********second order running alpha_strong
  mstp(33) = 1
  parp(31) = hipr1(17)
!                        ********inclusion of K factor
  mstp(51) = 3
!                        ********Duke-Owens set 1 structure functions
  mstp(61) = 1
!                        ********INITIAL STATE RADIATION
  mstp(71) = 1
!                        ********FINAL STATE RADIATION
  If (ihpr2(2)==0 .Or. ihpr2(2)==2) mstp(61) = 0
  If (ihpr2(2)==0 .Or. ihpr2(2)==1) mstp(71) = 0
!
  mstp(81) = 0
!                        ******** NO MULTIPLE INTERACTION
  mstp(82) = 1
!                        *******STRUCTURE OF MUTLIPLE INTERACTION
  mstp(111) = 0
!                ********frag off(have to be done by local call)
  If (ihpr2(10)==0) mstp(122) = 0
!                ********No printout of initialization information
  parp(81) = hipr1(8)
  ckin(5) = hipr1(8)
  ckin(3) = hipr1(8)
  ckin(4) = hipr1(9)
  If (hipr1(9)<=hipr1(8)) ckin(4) = -1.0
  ckin(9) = -10.0
  ckin(10) = 10.0
  msel = 0
  Do isub = 1, 200
    msub(isub) = 0
  End Do
  msub(11) = 1
  msub(12) = 1
  msub(13) = 1
  msub(28) = 1
  msub(53) = 1
  msub(68) = 1
  msub(81) = 1
  msub(82) = 1
  Do j = 1, min(8, mdcy(21,3))
    mdme(mdcy(21,2)+j-1, 1) = 0
  End Do
  isel = 4
  If (hint1(1)>=20.0 .And. ihpr2(18)==1) isel = 5
  mdme(mdcy(21,2)+isel-1, 1) = 1
!                        ********QCD subprocesses
  msub(14) = 1
  msub(18) = 1
  msub(29) = 1
!                       ******* direct photon production
  150 If (ini(itype)/=0) Goto 800
  Goto 400
!
!        *****triggered subprocesses, jet, photon, heavy quark and DY
!
  160 itype = 4 + itype
  If (itrig==ilast) Goto 260
  parp(81) = abs(hipr1(10)) - 0.25
  ckin(5) = abs(hipr1(10)) - 0.25
  ckin(3) = abs(hipr1(10)) - 0.25
  ckin(4) = abs(hipr1(10)) + 0.25
  If (hipr1(10)<hipr1(8)) ckin(4) = -1.0
!
  msel = 0
  Do isub = 1, 200
    msub(isub) = 0
  End Do
  If (ihpr2(3)==1) Then
    msub(11) = 1
    msub(12) = 1
    msub(13) = 1
    msub(28) = 1
    msub(53) = 1
    msub(68) = 1
    msub(81) = 1
    msub(82) = 1
    msub(14) = 1
    msub(18) = 1
    msub(29) = 1
    Do j = 1, min(8, mdcy(21,3))
      mdme(mdcy(21,2)+j-1, 1) = 0
    End Do
    isel = 4
    If (hint1(1)>=20.0 .And. ihpr2(18)==1) isel = 5
    mdme(mdcy(21,2)+isel-1, 1) = 1
!                        ********QCD subprocesses
  Else If (ihpr2(3)==2) Then
    msub(14) = 1
    msub(18) = 1
    msub(29) = 1
!                ********Direct photon production
!                q+qbar->g+gamma,q+qbar->gamma+gamma, q+g->q+gamma
  Else If (ihpr2(3)==3) Then
    ckin(3) = max(0.0, hipr1(10))
    ckin(5) = hipr1(8)
    parp(81) = hipr1(8)
    msub(81) = 1
    msub(82) = 1
    Do j = 1, min(8, mdcy(21,3))
      mdme(mdcy(21,2)+j-1, 1) = 0
    End Do
    isel = 4
    If (hint1(1)>=20.0 .And. ihpr2(18)==1) isel = 5
    mdme(mdcy(21,2)+isel-1, 1) = 1
!             **********Heavy quark production
  End If
  260 If (ini(itype)/=0) Goto 800
!
!
  400 ini(itype) = 1
  If (ihpr2(10)==0) mstp(122) = 0
  If (nfp(jp,4)==2212) Then
    beam = 'P'
  Else If (nfp(jp,4)==-2212) Then
    beam = 'P~'
  Else If (nfp(jp,4)==2112) Then
    beam = 'N'
  Else If (nfp(jp,4)==-2112) Then
    beam = 'N~'
  Else If (nfp(jp,4)==211) Then
    beam = 'PI+'
  Else If (nfp(jp,4)==-211) Then
    beam = 'PI-'
  Else If (nfp(jp,4)==321) Then
    beam = 'PI+'
  Else If (nfp(jp,4)==-321) Then
    beam = 'PI-'
  Else
    Write (6, *) 'unavailable beam type', nfp(jp, 4)
  End If
  If (nft(jt,4)==2212) Then
    targ = 'P'
  Else If (nft(jt,4)==-2212) Then
    targ = 'P~'
  Else If (nft(jt,4)==2112) Then
    targ = 'N'
  Else If (nft(jt,4)==-2112) Then
    targ = 'N~'
  Else If (nft(jt,4)==211) Then
    targ = 'PI+'
  Else If (nft(jt,4)==-211) Then
    targ = 'PI-'
  Else If (nft(jt,4)==321) Then
    targ = 'PI+'
  Else If (nft(jt,4)==-321) Then
    targ = 'PI-'
  Else
    Write (6, *) 'unavailable target type', nft(jt, 4)
  End If
!
  ihnt2(16) = 1
!       ******************indicate for initialization use when
!                         structure functions are called in PYTHIA
!
  Call pyinit('CMS', beam, targ, hint1(1))
  mint4 = mint(44)
  mint5 = mint(45)
  mint44(itype) = mint(44)
  mint45(itype) = mint(45)
  atxs(0) = xsec(0, 1)
  xsec0(itype, 0) = xsec(0, 1)
  Do i = 1, 200
    atxs(i) = xsec(i, 1)
    xsec0(itype, i) = xsec(i, 1)
    Do j = 1, 20
      atco(i, j) = coef(i, j)
      coef0(itype, i, j) = coef(i, j)
    End Do
  End Do
!
  ihnt2(16) = 0
!
  Return
!                ********Store the initialization information for
!                                late use
!
!
  800 mint(44) = mint44(itype)
  mint(45) = mint45(itype)
  mint4 = mint(44)
  mint5 = mint(45)
  xsec(0, 1) = xsec0(itype, 0)
  atxs(0) = xsec(0, 1)
  Do i = 1, 200
    xsec(i, 1) = xsec0(itype, i)
    atxs(i) = xsec(i, 1)
    Do j = 1, 20
      coef(i, j) = coef0(itype, i, j)
      atco(i, j) = coef(i, j)
    End Do
  End Do
  ilast = itrig
  mint(11) = nfp(jp, 4)
  mint(12) = nft(jt, 4)
  Return
End Subroutine jetini
!
!
!
Subroutine hijini
  Parameter (maxstr=150001)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
!        COMMON/HJJET4/NDR,IADR(900,2),KFDR(900),PDR(900,5)
  Common /hjjet4/ndr, iadr(maxstr, 2), kfdr(maxstr), pdr(maxstr, 5)
!c      SAVE /HJJET4/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!****************Reset the momentum of initial particles************
!             and assign flavors to the proj and targ string       *
!*******************************************************************
  nsg = 0
  ndr = 0
  ipp = 2212
  ipt = 2212
  If (ihnt2(5)/=0) ipp = ihnt2(5)
  If (ihnt2(6)/=0) ipt = ihnt2(6)
!                ********in case the proj or targ is a hadron.
!
  Do i = 1, ihnt2(1)
    pp(i, 1) = 0.0
    pp(i, 2) = 0.0
    pp(i, 3) = sqrt(hint1(1)**2/4.0-hint1(8)**2)
    pp(i, 4) = hint1(1)/2
    pp(i, 5) = hint1(8)
    pp(i, 6) = 0.0
    pp(i, 7) = 0.0
    pp(i, 8) = 0.0
    pp(i, 9) = 0.0
    pp(i, 10) = 0.0
!bzdbg2/22/99
!test OFF
    pp(i, 11) = 0.0
    pp(i, 12) = 0.0
!bzdbg2/22/99end
    nfp(i, 3) = ipp
    nfp(i, 4) = ipp
    nfp(i, 5) = 0
    nfp(i, 6) = 0
    nfp(i, 7) = 0
    nfp(i, 8) = 0
    nfp(i, 9) = 0
    nfp(i, 10) = 0
    nfp(i, 11) = 0
    npj(i) = 0
    If (i>abs(ihnt2(2))) nfp(i, 3) = 2112

!lin-12/2012 correct NN differential cross section in HIJING:
    If (i>abs(ihnt2(2))) nfp(i, 4) = 2112

    Call attflv(nfp(i,3), idq, idqq)
    nfp(i, 1) = idq
    nfp(i, 2) = idqq
    nfp(i, 15) = -1
    If (abs(idq)>1000 .Or. (abs(idq*idqq)<100 .And. ranart(nseed)<0.5)) nfp(i, 15) = 1
    pp(i, 14) = ulmass(idq)
    pp(i, 15) = ulmass(idqq)
  End Do
!
  Do i = 1, ihnt2(3)
    pt(i, 1) = 0.0
    pt(i, 2) = 0.0
    pt(i, 3) = -sqrt(hint1(1)**2/4.0-hint1(9)**2)
    pt(i, 4) = hint1(1)/2.0
    pt(i, 5) = hint1(9)
    pt(i, 6) = 0.0
    pt(i, 7) = 0.0
    pt(i, 8) = 0.0
    pt(i, 9) = 0.0
    pt(i, 10) = 0.0
!test OFF
!bzdbg2/22/99
    pt(i, 11) = 0.0
    pt(i, 12) = 0.0
!bzdbg2/22/99end
    nft(i, 3) = ipt
    nft(i, 4) = ipt
    nft(i, 5) = 0
    nft(i, 6) = 0
    nft(i, 7) = 0
    nft(i, 8) = 0
    nft(i, 9) = 0
    nft(i, 10) = 0
    nft(i, 11) = 0
    ntj(i) = 0
    If (i>abs(ihnt2(4))) nft(i, 3) = 2112

!lin-12/2012 correct NN differential cross section in HIJING:
    If (i>abs(ihnt2(4))) nft(i, 4) = 2112

    Call attflv(nft(i,3), idq, idqq)
    nft(i, 1) = idq
    nft(i, 2) = idqq
    nft(i, 15) = 1
    If (abs(idq)>1000 .Or. (abs(idq*idqq)<100 .And. ranart(nseed)<0.5)) nft(i, 15) = -1
    pt(i, 14) = ulmass(idq)
    pt(i, 15) = ulmass(idqq)
  End Do
  Return
End Subroutine hijini
!
!
!
Subroutine attflv(id, idq, idqq)
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
  If (abs(id)<100) Then
    nsign = 1
    idq = id/100
    idqq = -id/10 + idq*10
    If (abs(idq)==3) nsign = -1
    idq = nsign*idq
    idqq = nsign*idqq
    If (idq<0) Then
      id0 = idq
      idq = idqq
      idqq = id0
    End If
    Return
  End If
!                ********return ID of quark(IDQ) and anti-quark(IDQQ)
!                        for pions and kaons
!
!        Return LU ID for quarks and diquarks for proton(ID=2212)
!        anti-proton(ID=-2212) and nuetron(ID=2112)
!        LU ID for d=1,u=2, (ud)0=2101, (ud)1=2103,
!       (dd)1=1103,(uu)1=2203.
!        Use SU(6)  weight  proton=1/3d(uu)1 + 1/6u(ud)1 + 1/2u(ud)0
!                          nurtron=1/3u(dd)1 + 1/6d(ud)1 + 1/2d(ud)0
!
  idq = 2
  If (abs(id)==2112) idq = 1
  idqq = 2101
  x = ranart(nseed)
  If (x<=0.5) Goto 30
  If (x>0.666667) Goto 10
  idqq = 2103
  Goto 30
  10 idq = 1
  idqq = 2203
  If (abs(id)==2112) Then
    idq = 2
    idqq = 1103
  End If
  30 If (id<0) Then
    id00 = idqq
    idqq = -idq
    idq = -id00
  End If
  Return
End Subroutine attflv
!
!*******************************************************************
!        This subroutine performs elastic scatterings and possible
!        elastic cascading within their own nuclei
!*******************************************************************
Subroutine hijcsc(jp, jt)
  Dimension psc1(5), psc2(5)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Save
  If (jp==0 .Or. jt==0) Goto 25
  Do i = 1, 5
    psc1(i) = pp(jp, i)
    psc2(i) = pt(jt, i)
  End Do
  Call hijels(psc1, psc2)
  dpp1 = psc1(1) - pp(jp, 1)
  dpp2 = psc1(2) - pp(jp, 2)
  dpt1 = psc2(1) - pt(jt, 1)
  dpt2 = psc2(2) - pt(jt, 2)
  pp(jp, 6) = pp(jp, 6) + dpp1/2.0
  pp(jp, 7) = pp(jp, 7) + dpp2/2.0
  pp(jp, 8) = pp(jp, 8) + dpp1/2.0
  pp(jp, 9) = pp(jp, 9) + dpp2/2.0
  pt(jt, 6) = pt(jt, 6) + dpt1/2.0
  pt(jt, 7) = pt(jt, 7) + dpt2/2.0
  pt(jt, 8) = pt(jt, 8) + dpt1/2.0
  pt(jt, 9) = pt(jt, 9) + dpt2/2.0
  Do i = 1, 4
    pp(jp, i) = psc1(i)
    pt(jt, i) = psc2(i)
  End Do
  nfp(jp, 5) = max(1, nfp(jp,5))
  nft(jt, 5) = max(1, nft(jt,5))
!                ********Perform elastic scattering between JP and JT
  Return
!                ********The following is for possible elastic cascade
!
  25 If (jp==0) Goto 45
  pabs = sqrt(pp(jp,1)**2+pp(jp,2)**2+pp(jp,3)**2)
  bx = pp(jp, 1)/pabs
  by = pp(jp, 2)/pabs
  bz = pp(jp, 3)/pabs
  Do i = 1, ihnt2(1)
    If (i==jp) Goto 40
    dx = yp(1, i) - yp(1, jp)
    dy = yp(2, i) - yp(2, jp)
    dz = yp(3, i) - yp(3, jp)
    dis = dx*bx + dy*by + dz*bz
    If (dis<=0) Goto 40
    bb = dx**2 + dy**2 + dz**2 - dis**2
    r2 = bb*hipr1(40)/hipr1(31)/0.1
!                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
    gs = 1.0 - exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(r2))**2
    gs0 = 1.0 - exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(0.0))**2
    If (ranart(nseed)>gs/gs0) Goto 40
    Do k = 1, 5
      psc1(k) = pp(jp, k)
      psc2(k) = pp(i, k)
    End Do
    Call hijels(psc1, psc2)
    dpp1 = psc1(1) - pp(jp, 1)
    dpp2 = psc1(2) - pp(jp, 2)
    dpt1 = psc2(1) - pp(i, 1)
    dpt2 = psc2(2) - pp(i, 2)
    pp(jp, 6) = pp(jp, 6) + dpp1/2.0
    pp(jp, 7) = pp(jp, 7) + dpp2/2.0
    pp(jp, 8) = pp(jp, 8) + dpp1/2.0
    pp(jp, 9) = pp(jp, 9) + dpp2/2.0
    pp(i, 6) = pp(i, 6) + dpt1/2.0
    pp(i, 7) = pp(i, 7) + dpt2/2.0
    pp(i, 8) = pp(i, 8) + dpt1/2.0
    pp(i, 9) = pp(i, 9) + dpt2/2.0
    Do k = 1, 5
      pp(jp, k) = psc1(k)
      pp(i, k) = psc2(k)
    End Do
    nfp(i, 5) = max(1, nfp(i,5))
    Goto 45
  40 End Do
  45 If (jt==0) Goto 80
!lin 50        PABS=SQRT(PT(JT,1)**2+PT(JT,2)**2+PT(JT,3)**2)
  pabs = sqrt(pt(jt,1)**2+pt(jt,2)**2+pt(jt,3)**2)
  bx = pt(jt, 1)/pabs
  by = pt(jt, 2)/pabs
  bz = pt(jt, 3)/pabs
  Do i = 1, ihnt2(3)
    If (i==jt) Goto 70
    dx = yt(1, i) - yt(1, jt)
    dy = yt(2, i) - yt(2, jt)
    dz = yt(3, i) - yt(3, jt)
    dis = dx*bx + dy*by + dz*bz
    If (dis<=0) Goto 70
    bb = dx**2 + dy**2 + dz**2 - dis**2
    r2 = bb*hipr1(40)/hipr1(31)/0.1
!                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
    gs = (1.0-exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(r2)))**2
    gs0 = (1.0-exp(-(hipr1(30)+hint1(11))/hipr1(31)/2.0*romg(0.0)))**2
    If (ranart(nseed)>gs/gs0) Goto 70
    Do k = 1, 5
      psc1(k) = pt(jt, k)
      psc2(k) = pt(i, k)
    End Do
    Call hijels(psc1, psc2)
    dpp1 = psc1(1) - pt(jt, 1)
    dpp2 = psc1(2) - pt(jt, 2)
    dpt1 = psc2(1) - pt(i, 1)
    dpt2 = psc2(2) - pt(i, 2)
    pt(jt, 6) = pt(jt, 6) + dpp1/2.0
    pt(jt, 7) = pt(jt, 7) + dpp2/2.0
    pt(jt, 8) = pt(jt, 8) + dpp1/2.0
    pt(jt, 9) = pt(jt, 9) + dpp2/2.0
    pt(i, 6) = pt(i, 6) + dpt1/2.0
    pt(i, 7) = pt(i, 7) + dpt2/2.0
    pt(i, 8) = pt(i, 8) + dpt1/2.0
    pt(i, 9) = pt(i, 9) + dpt2/2.0
    Do k = 1, 5
      pt(jt, k) = psc1(k)
      pt(i, k) = psc2(k)
    End Do
    nft(i, 5) = max(1, nft(i,5))
    Goto 80
  70 End Do
  80 Return
End Subroutine hijcsc
!
!
!*******************************************************************
!This subroutine performs elastic scattering between two nucleons
!
!*******************************************************************
Subroutine hijels(psc1, psc2)
  Implicit Double Precision (D)
  Dimension psc1(5), psc2(5)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
!
  cc = 1.0 - hint1(12)/hint1(13)
  rr = (1.0-cc)*hint1(13)/hint1(12)/(1.0-hipr1(33)) - 1.0
  bb = 0.5*(3.0+rr+sqrt(9.0+10.0*rr+rr**2))
  ep = sqrt((psc1(1)-psc2(1))**2+(psc1(2)-psc2(2))**2+(psc1(3)-psc2(3))**2)
  If (ep<=0.1) Return
  els0 = 98.0/ep + 52.0*(1.0+rr)**2
  pcm1 = psc1(1) + psc2(1)
  pcm2 = psc1(2) + psc2(2)
  pcm3 = psc1(3) + psc2(3)
  ecm = psc1(4) + psc2(4)
  am1 = psc1(5)**2
  am2 = psc2(5)**2
  amm = ecm**2 - pcm1**2 - pcm2**2 - pcm3**2
  If (amm<=psc1(5)+psc2(5)) Return
!                ********elastic scattering only when approaching
!                                to each other
  pmax = (amm**2+am1**2+am2**2-2.0*amm*am1-2.0*amm*am2-2.0*am1*am2)/4.0/amm
  pmax = abs(pmax)
  20 tt = ranart(nseed)*min(pmax, 1.5)
  els = 98.0*exp(-2.8*tt)/ep + 52.0*exp(-9.2*tt)*(1.0+rr*exp(-4.6*(bb-1.0)*tt))**2
  If (ranart(nseed)>els/els0) Goto 20
  phi = 2.0*hipr1(40)*ranart(nseed)
!
  dbx = dble(pcm1/ecm)
  dby = dble(pcm2/ecm)
  dbz = dble(pcm3/ecm)
  db = dsqrt(dbx**2+dby**2+dbz**2)
  If (db>0.99999999D0) Then
    dbx = dbx*(0.99999999D0/db)
    dby = dby*(0.99999999D0/db)
    dbz = dbz*(0.99999999D0/db)
    db = 0.99999999D0
    Write (6, *) ' (HIJELS) boost vector too large'
!                ********Rescale boost vector if too close to unity.
  End If
  dga = 1D0/sqrt(1D0-db**2)
!
  dp1 = dble(sqrt(tt)*sin(phi))
  dp2 = dble(sqrt(tt)*cos(phi))
  dp3 = dble(sqrt(pmax-tt))
  dp4 = dble(sqrt(pmax+am1))
  dbp = dbx*dp1 + dby*dp2 + dbz*dp3
  dgabp = dga*(dga*dbp/(1D0+dga)+dp4)
  psc1(1) = sngl(dp1+dgabp*dbx)
  psc1(2) = sngl(dp2+dgabp*dby)
  psc1(3) = sngl(dp3+dgabp*dbz)
  psc1(4) = sngl(dga*(dp4+dbp))
!
  dp1 = -dble(sqrt(tt)*sin(phi))
  dp2 = -dble(sqrt(tt)*cos(phi))
  dp3 = -dble(sqrt(pmax-tt))
  dp4 = dble(sqrt(pmax+am2))
  dbp = dbx*dp1 + dby*dp2 + dbz*dp3
  dgabp = dga*(dga*dbp/(1D0+dga)+dp4)
  psc2(1) = sngl(dp1+dgabp*dbx)
  psc2(2) = sngl(dp2+dgabp*dby)
  psc2(3) = sngl(dp3+dgabp*dbz)
  psc2(4) = sngl(dga*(dp4+dbp))
  Return
End Subroutine hijels
!
!
!*******************************************************************
!                                                                      *
!                Subroutine HIJSFT                                   *
!                                                                   *
!  Scatter two excited strings, JP from proj and JT from target    *
!*******************************************************************
Subroutine hijsft(jp, jt, jout, ierror)
  Parameter (maxstr=150001)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hijdat/hidat0(10, 10), hidat(10)
!c      SAVE /HIJDAT/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
!lin-4/25/01
!        COMMON/HJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
!     &                K2SG(900,100),PXSG(900,100),PYSG(900,100),
!     &                PZSG(900,100),PESG(900,100),PMSG(900,100)
!c      SAVE /HJJET2/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /dpmcm1/jjp, jjt, amp, amt, apx0, atx0, ampn, amtn, amp0, amt0, nfdp, nfdt, wp, wm, sw, xremp, xremt, dpkc1, dpkc2, pp11, pp12, pt11, pt12, ptp2, ptt2
!c      SAVE /DPMCM1/
  Common /dpmcm2/ndpm, kdpm(20, 2), pdpm1(20, 5), pdpm2(20, 5)
!c      SAVE /DPMCM2/
  Save
!*******************************************************************
!        JOUT-> the number
!        of hard scatterings preceding this soft collision.
!       IHNT2(13)-> 1=
!        double diffrac 2=single diffrac, 3=non-single diffrac.
!*******************************************************************
  ierror = 0
  jjp = jp
  jjt = jt
  ndpm = 0
!        IOPMAIN=0
  If (jp>ihnt2(1) .Or. jt>ihnt2(3)) Return

  epp = pp(jp, 4) + pp(jp, 3)
  epm = pp(jp, 4) - pp(jp, 3)
  etp = pt(jt, 4) + pt(jt, 3)
  etm = pt(jt, 4) - pt(jt, 3)

  wp = epp + etp
  wm = epm + etm
  sw = wp*wm
!                ********total W+,W- and center-of-mass energy

  If (wp<0.0 .Or. wm<0.0) Goto 1000

  If (jout==0) Then
    If (epp<0.0) Goto 1000
    If (epm<0.0) Goto 1000
    If (etp<0.0) Goto 1000
    If (etm<0.0) Goto 1000
    If (epp/(epm+0.01)<=etp/(etm+0.01)) Return
  End If
!                ********For strings which does not follow a jet-prod,
!                        scatter only if Ycm(JP)>Ycm(JT). When jets
!                        are produced just before this collision
!                        this requirement has already be enforced
!                        (see SUBROUTINE HIJHRD)
  ihnt2(11) = jp
  ihnt2(12) = jt
!
!
!
  miss = 0
  pkc1 = 0.0
  pkc2 = 0.0
  pkc11 = 0.0
  pkc12 = 0.0
  pkc21 = 0.0
  pkc22 = 0.0
  dpkc11 = 0.0
  dpkc12 = 0.0
  dpkc21 = 0.0
  dpkc22 = 0.0
  If (nfp(jp,10)==1 .Or. nft(jt,10)==1) Then
    If (nfp(jp,10)==1) Then
      phi1 = ulangl(pp(jp,10), pp(jp,11))
      ppjet = sqrt(pp(jp,10)**2+pp(jp,11)**2)
      pkc1 = ppjet
      pkc11 = pp(jp, 10)
      pkc12 = pp(jp, 11)
    End If
    If (nft(jt,10)==1) Then
      phi2 = ulangl(pt(jt,10), pt(jt,11))
      ptjet = sqrt(pt(jt,10)**2+pt(jt,11)**2)
      pkc2 = ptjet
      pkc21 = pt(jt, 10)
      pkc22 = pt(jt, 11)
    End If
    If (ihpr2(4)>0 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
      If (nfp(jp,10)==0) Then
        phi = -phi2
      Else If (nft(jt,10)==0) Then
        phi = phi1
      Else
        phi = (phi1+phi2-hipr1(40))/2.0
      End If
      bx = hint1(19)*cos(hint1(20))
      by = hint1(19)*sin(hint1(20))
      xp0 = yp(1, jp)
      yp0 = yp(2, jp)
      xt0 = yt(1, jt) + bx
      yt0 = yt(2, jt) + by
      r1 = max(1.2*ihnt2(1)**0.3333333, sqrt(xp0**2+yp0**2))
      r2 = max(1.2*ihnt2(3)**0.3333333, sqrt((xt0-bx)**2+(yt0-by)**2))
      If (abs(cos(phi))<1.0E-5) Then
        dd1 = r1
        dd2 = r1
        dd3 = abs(by+sqrt(r2**2-(xp0-bx)**2)-yp0)
        dd4 = abs(by-sqrt(r2**2-(xp0-bx)**2)-yp0)
        Goto 5
      End If
      bb = 2.0*sin(phi)*(cos(phi)*yp0-sin(phi)*xp0)
      cc = (yp0**2-r1**2)*cos(phi)**2 + xp0*sin(phi)*(xp0*sin(phi)-2.0*yp0*cos(phi))
      dd = bb**2 - 4.0*cc
      If (dd<0.0) Goto 10
      xx1 = (-bb+sqrt(dd))/2.0
      xx2 = (-bb-sqrt(dd))/2.0
      dd1 = abs((xx1-xp0)/cos(phi))
      dd2 = abs((xx2-xp0)/cos(phi))
!
      bb = 2.0*sin(phi)*(cos(phi)*(yt0-by)-sin(phi)*xt0) - 2.0*bx
      cc = (bx**2+(yt0-by)**2-r2**2)*cos(phi)**2 + xt0*sin(phi)*(xt0*sin(phi)-2.0*cos(phi)*(yt0-by)) - 2.0*bx*sin(phi)*(cos(phi)*(yt0-by)-sin(phi)*xt0)
      dd = bb**2 - 4.0*cc
      If (dd<0.0) Goto 10
      xx1 = (-bb+sqrt(dd))/2.0
      xx2 = (-bb-sqrt(dd))/2.0
      dd3 = abs((xx1-xt0)/cos(phi))
      dd4 = abs((xx2-xt0)/cos(phi))
!
      5 dd1 = min(dd1, dd3)
      dd2 = min(dd2, dd4)
      If (dd1<hipr1(13)) dd1 = 0.0
      If (dd2<hipr1(13)) dd2 = 0.0
      If (nfp(jp,10)==1 .And. ppjet>hipr1(11)) Then
        dp1 = dd1*hipr1(14)/2.0
        dp1 = min(dp1, ppjet-hipr1(11))
        pkc1 = ppjet - dp1
        dpx1 = cos(phi1)*dp1
        dpy1 = sin(phi1)*dp1
        pkc11 = pp(jp, 10) - dpx1
        pkc12 = pp(jp, 11) - dpy1
        If (dp1>0.0) Then
          cthep = pp(jp, 12)/sqrt(pp(jp,12)**2+ppjet**2)
          dpz1 = dp1*cthep/sqrt(1.0-cthep**2)
          dpe1 = sqrt(dpx1**2+dpy1**2+dpz1**2)
          eppprm = pp(jp, 4) + pp(jp, 3) - dpe1 - dpz1
          epmprm = pp(jp, 4) - pp(jp, 3) - dpe1 + dpz1
          If (eppprm<=0.0 .Or. epmprm<=0.0) Goto 15
          epp = eppprm
          epm = epmprm
          pp(jp, 10) = pkc11
          pp(jp, 11) = pkc12
          npj(jp) = npj(jp) + 1
          kfpj(jp, npj(jp)) = 21
          pjpx(jp, npj(jp)) = dpx1
          pjpy(jp, npj(jp)) = dpy1
          pjpz(jp, npj(jp)) = dpz1
          pjpe(jp, npj(jp)) = dpe1
          pjpm(jp, npj(jp)) = 0.0
          pp(jp, 3) = pp(jp, 3) - dpz1
          pp(jp, 4) = pp(jp, 4) - dpe1
        End If
      End If
      15 If (nft(jt,10)==1 .And. ptjet>hipr1(11)) Then
        dp2 = dd2*hipr1(14)/2.0
        dp2 = min(dp2, ptjet-hipr1(11))
        pkc2 = ptjet - dp2
        dpx2 = cos(phi2)*dp2
        dpy2 = sin(phi2)*dp2
        pkc21 = pt(jt, 10) - dpx2
        pkc22 = pt(jt, 11) - dpy2
        If (dp2>0.0) Then
          cthet = pt(jt, 12)/sqrt(pt(jt,12)**2+ptjet**2)
          dpz2 = dp2*cthet/sqrt(1.0-cthet**2)
          dpe2 = sqrt(dpx2**2+dpy2**2+dpz2**2)
          etpprm = pt(jt, 4) + pt(jt, 3) - dpe2 - dpz2
          etmprm = pt(jt, 4) - pt(jt, 3) - dpe2 + dpz2
          If (etpprm<=0.0 .Or. etmprm<=0.0) Goto 16
          etp = etpprm
          etm = etmprm
          pt(jt, 10) = pkc21
          pt(jt, 11) = pkc22
          ntj(jt) = ntj(jt) + 1
          kftj(jt, ntj(jt)) = 21
          pjtx(jt, ntj(jt)) = dpx2
          pjty(jt, ntj(jt)) = dpy2
          pjtz(jt, ntj(jt)) = dpz2
          pjte(jt, ntj(jt)) = dpe2
          pjtm(jt, ntj(jt)) = 0.0
          pt(jt, 3) = pt(jt, 3) - dpz2
          pt(jt, 4) = pt(jt, 4) - dpe2
        End If
      End If
      16 dpkc11 = -(pp(jp,10)-pkc11)/2.0
      dpkc12 = -(pp(jp,11)-pkc12)/2.0
      dpkc21 = -(pt(jt,10)-pkc21)/2.0
      dpkc22 = -(pt(jt,11)-pkc22)/2.0
      wp = epp + etp
      wm = epm + etm
      sw = wp*wm
    End If
  End If
!                ********If jet is quenched the pt from valence quark
!                        hard scattering has to reduced by d*kapa
!
!
  10 ptp02 = pp(jp, 1)**2 + pp(jp, 2)**2
  ptt02 = pt(jt, 1)**2 + pt(jt, 2)**2
!
  amq = max(pp(jp,14)+pp(jp,15), pt(jt,14)+pt(jt,15))
  amx = hipr1(1) + amq
!                ********consider mass cut-off for strings which
!                        must also include quark's mass
  amp0 = amx
  dpm0 = amx
  nfdp = 0
  If (nfp(jp,5)<=2 .And. nfp(jp,3)/=0) Then
    amp0 = ulmass(nfp(jp,3))
    nfdp = nfp(jp, 3) + 2*nfp(jp, 3)/abs(nfp(jp,3))
    dpm0 = ulmass(nfdp)
    If (dpm0<=0.0) Then
      nfdp = nfdp - 2*nfdp/abs(nfdp)
      dpm0 = ulmass(nfdp)
    End If
  End If
  amt0 = amx
  dtm0 = amx
  nfdt = 0
  If (nft(jt,5)<=2 .And. nft(jt,3)/=0) Then
    amt0 = ulmass(nft(jt,3))
    nfdt = nft(jt, 3) + 2*nft(jt, 3)/abs(nft(jt,3))
    dtm0 = ulmass(nfdt)
    If (dtm0<=0.0) Then
      nfdt = nfdt - 2*nfdt/abs(nfdt)
      dtm0 = ulmass(nfdt)
    End If
  End If
!
  ampn = sqrt(amp0**2+ptp02)
  amtn = sqrt(amt0**2+ptt02)
  snn = (ampn+amtn)**2 + 0.001
!
  If (sw<snn+0.001) Goto 4000
!                ********Scatter only if SW>SNN
!*****give some PT kick to the two exited strings******************
!lin 20        SWPTN=4.0*(MAX(AMP0,AMT0)**2+MAX(PTP02,PTT02))
  swptn = 4.0*(max(amp0,amt0)**2+max(ptp02,ptt02))
  swptd = 4.0*(max(dpm0,dtm0)**2+max(ptp02,ptt02))
  swptx = 4.0*(amx**2+max(ptp02,ptt02))
  If (sw<=swptn) Then
    pkcmx = 0.0
  Else If (sw>swptn .And. sw<=swptd .And. npj(jp)==0 .And. ntj(jt)==0) Then
    pkcmx = sqrt(sw/4.0-max(amp0,amt0)**2) - sqrt(max(ptp02,ptt02))
  Else If (sw>swptd .And. sw<=swptx .And. npj(jp)==0 .And. ntj(jt)==0) Then
    pkcmx = sqrt(sw/4.0-max(dpm0,dtm0)**2) - sqrt(max(ptp02,ptt02))
  Else If (sw>swptx) Then
    pkcmx = sqrt(sw/4.0-amx**2) - sqrt(max(ptp02,ptt02))
  End If
!                ********maximun PT kick
!*********************************************************
!
  If (nfp(jp,10)==1 .Or. nft(jt,10)==1) Then
    If (pkc1>pkcmx) Then
      pkc1 = pkcmx
      pkc11 = pkc1*cos(phi1)
      pkc12 = pkc1*sin(phi1)
      dpkc11 = -(pp(jp,10)-pkc11)/2.0
      dpkc12 = -(pp(jp,11)-pkc12)/2.0
    End If
    If (pkc2>pkcmx) Then
      pkc2 = pkcmx
      pkc21 = pkc2*cos(phi2)
      pkc22 = pkc2*sin(phi2)
      dpkc21 = -(pt(jt,10)-pkc21)/2.0
      dpkc22 = -(pt(jt,11)-pkc22)/2.0
    End If
    dpkc1 = dpkc11 + dpkc21
    dpkc2 = dpkc12 + dpkc22
    nfp(jp, 10) = -nfp(jp, 10)
    nft(jt, 10) = -nft(jt, 10)
    Goto 40
  End If
!                ********If the valence quarks had a hard-collision
!                        the pt kick is the pt from hard-collision.
  isng = 0
  If (ihpr2(13)/=0 .And. ranart(nseed)<=hidat(4)) isng = 1
  If ((nfp(jp,5)==3 .Or. nft(jt,5)==3) .Or. (npj(jp)/=0 .Or. nfp(jp,10)/=0) .Or. (ntj(jt)/=0 .Or. nft(jt,10)/=0)) isng = 0
!
!               ********decite whether to have single-diffractive
  If (ihpr2(5)==0) Then
    pkc = hipr1(2)*sqrt(-alog(1.0-ranart(nseed)*(1.0-exp(-pkcmx**2/hipr1(2)**2))))
    Goto 30
  End If

!lin-10/28/02 get rid of argument usage mismatch in HIRND2():
!        PKC=HIRND2(3,0.0,PKCMX**2)
  xminhi = 0.0
  xmaxhi = pkcmx**2
  pkc = hirnd2(3, xminhi, xmaxhi)

  pkc = sqrt(pkc)
  If (pkc>hipr1(20)) pkc = hipr1(2)*sqrt(-alog(exp(-hipr1(20)**2/hipr1(2)**2)-ranart(nseed)*(exp(-hipr1(20)**2/hipr1(2)**2)-exp(-pkcmx**2/hipr1(2)**2))))
!
  If (isng==1) pkc = 0.65*sqrt(-alog(1.0-ranart(nseed)*(1.0-exp(-pkcmx**2/0.65**2))))
!                        ********select PT kick
  30 phi0 = 2.0*hipr1(40)*ranart(nseed)
  pkc11 = pkc*sin(phi0)
  pkc12 = pkc*cos(phi0)
  pkc21 = -pkc11
  pkc22 = -pkc12
  dpkc1 = 0.0
  dpkc2 = 0.0
  40 pp11 = pp(jp, 1) + pkc11 - dpkc1
  pp12 = pp(jp, 2) + pkc12 - dpkc2
  pt11 = pt(jt, 1) + pkc21 - dpkc1
  pt12 = pt(jt, 2) + pkc22 - dpkc2
  ptp2 = pp11**2 + pp12**2
  ptt2 = pt11**2 + pt12**2
!
  ampn = sqrt(amp0**2+ptp2)
  amtn = sqrt(amt0**2+ptt2)
  snn = (ampn+amtn)**2 + 0.001
!***************************************
  wp = epp + etp
  wm = epm + etm
  sw = wp*wm
!****************************************
  If (sw<snn) Then
    miss = miss + 1
    If (miss<=100) Then
      pkc = 0.0
      Goto 30
    End If
    If (ihpr2(10)/=0) Write (6, *) 'Error occured in Pt kick section of HIJSFT'
    Goto 4000
  End If
!******************************************************************
  ampd = sqrt(dpm0**2+ptp2)
  amtd = sqrt(dtm0**2+ptt2)

  ampx = sqrt(amx**2+ptp2)
  amtx = sqrt(amx**2+ptt2)

  dpn = ampn**2/sw
  dtn = amtn**2/sw
  dpd = ampd**2/sw
  dtd = amtd**2/sw
  dpx = ampx**2/sw
  dtx = amtx**2/sw
!
  spntd = (ampn+amtd)**2
  spntx = (ampn+amtx)**2
!                        ********CM energy if proj=N,targ=N*
  spdtn = (ampd+amtn)**2
  spxtn = (ampx+amtn)**2
!                        ********CM energy if proj=N*,targ=N
  spdtx = (ampd+amtx)**2
  spxtd = (ampx+amtd)**2
  sdd = (ampd+amtd)**2
  sxx = (ampx+amtx)**2

!
!
!                ********CM energy if proj=delta, targ=delta
!****************There are many different cases**********
!        IF(IHPR2(15).EQ.1) GO TO 500
!
!                ********to have DPM type soft interactions
!
!lin 45        CONTINUE
  If (sw>sxx+0.001) Then
    If (isng==0) Then
      d1 = dpx
      d2 = dtx
      nfp3 = 0
      nft3 = 0
      Goto 400
    Else
!**** 5/30/1998 this is identical to the above statement. Added to
!**** avoid questional branching to block.
      If ((nfp(jp,5)==3 .And. nft(jt,5)==3) .Or. (npj(jp)/=0 .Or. nfp(jp,10)/=0) .Or. (ntj(jt)/=0 .Or. nft(jt,10)/=0)) Then
        d1 = dpx
        d2 = dtx
        nfp3 = 0
        nft3 = 0
        Goto 400
      End If
!                ********do not allow excited strings to have
!                        single-diffr
      If (ranart(nseed)>0.5 .Or. (nft(jt,5)>2 .Or. ntj(jt)/=0 .Or. nft(jt,10)/=0)) Then
        d1 = dpn
        d2 = dtx
        nfp3 = nfp(jp, 3)
        nft3 = 0
        Goto 220
      Else
        d1 = dpx
        d2 = dtn
        nfp3 = 0
        nft3 = nft(jt, 3)
        Goto 240
      End If
!                ********have single diffractive collision
    End If
  Else If (sw>max(spdtx,spxtd)+0.001 .And. sw<=sxx+0.001) Then
    If (((npj(jp)==0 .And. ntj(jt)==0 .And. ranart(nseed)>0.5) .Or. (npj(jp)==0 .And. ntj(jt)/=0)) .And. nfp(jp,5)<=2) Then
      d1 = dpd
      d2 = dtx
      nfp3 = nfdp
      nft3 = 0
      Goto 220
    Else If (ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtd
      nfp3 = 0
      nft3 = nfdt
      Goto 240
    End If
    Goto 4000
  Else If (sw>min(spdtx,spxtd)+0.001 .And. sw<=max(spdtx,spxtd)+0.001) Then
    If (spdtx<=spxtd .And. npj(jp)==0 .And. nfp(jp,5)<=2) Then
      d1 = dpd
      d2 = dtx
      nfp3 = nfdp
      nft3 = 0
      Goto 220
    Else If (spdtx>spxtd .And. ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtd
      nfp3 = 0
      nft3 = nfdt
      Goto 240
    End If
!*** 5/30/1998 added to avoid questional branching to another block
!*** this is identical to the statement following the next ELSE IF
    If (((npj(jp)==0 .And. ntj(jt)==0 .And. ranart(nseed)>0.5) .Or. (npj(jp)==0 .And. ntj(jt)/=0)) .And. nfp(jp,5)<=2) Then
      d1 = dpn
      d2 = dtx
      nfp3 = nfp(jp, 3)
      nft3 = 0
      Goto 220
    Else If (ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtn
      nfp3 = 0
      nft3 = nft(jt, 3)
      Goto 240
    End If
    Goto 4000
  Else If (sw>max(spntx,spxtn)+0.001 .And. sw<=min(spdtx,spxtd)+0.001) Then
    If (((npj(jp)==0 .And. ntj(jt)==0 .And. ranart(nseed)>0.5) .Or. (npj(jp)==0 .And. ntj(jt)/=0)) .And. nfp(jp,5)<=2) Then
      d1 = dpn
      d2 = dtx
      nfp3 = nfp(jp, 3)
      nft3 = 0
      Goto 220
    Else If (ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtn
      nfp3 = 0
      nft3 = nft(jt, 3)
      Goto 240
    End If
    Goto 4000
  Else If (sw>min(spntx,spxtn)+0.001 .And. sw<=max(spntx,spxtn)+0.001) Then
    If (spntx<=spxtn .And. npj(jp)==0 .And. nfp(jp,5)<=2) Then
      d1 = dpn
      d2 = dtx
      nfp3 = nfp(jp, 3)
      nft3 = 0
      Goto 220
    Else If (spntx>spxtn .And. ntj(jt)==0 .And. nft(jt,5)<=2) Then
      d1 = dpx
      d2 = dtn
      nfp3 = 0
      nft3 = nft(jt, 3)
      Goto 240
    End If
    Goto 4000
  Else If (sw<=min(spntx,spxtn)+0.001 .And. (npj(jp)/=0 .Or. ntj(jt)/=0)) Then
    Goto 4000
  Else If (sw<=min(spntx,spxtn)+0.001 .And. nfp(jp,5)>2 .And. nft(jt,5)>2) Then
    Goto 4000
  Else If (sw>sdd+0.001 .And. sw<=min(spntx,spxtn)+0.001) Then
    d1 = dpd
    d2 = dtd
    nfp3 = nfdp
    nft3 = nfdt
    Goto 100
  Else If (sw>max(spntd,spdtn)+0.001 .And. sw<=sdd+0.001) Then
    If (ranart(nseed)>0.5) Then
      d1 = dpd
      d2 = dtn
      nfp3 = nfdp
      nft3 = nft(jt, 3)
      Goto 100
    Else
      d1 = dpn
      d2 = dtd
      nfp3 = nfp(jp, 3)
      nft3 = nfdt
      Goto 100
    End If
  Else If (sw>min(spntd,spdtn)+0.001 .And. sw<=max(spntd,spdtn)+0.001) Then
    If (spntd>spdtn) Then
      d1 = dpd
      d2 = dtn
      nfp3 = nfdp
      nft3 = nft(jt, 3)
      Goto 100
    Else
      d1 = dpn
      d2 = dtd
      nfp3 = nfp(jp, 3)
      nft3 = nfdt
      Goto 100
    End If
  Else If (sw<=min(spntd,spdtn)+0.001) Then
    d1 = dpn
    d2 = dtn
    nfp3 = nfp(jp, 3)
    nft3 = nft(jt, 3)
    Goto 100
  End If
  Write (6, *) ' Error in HIJSFT: There is no path to here'
  Return
!
!***************  elastic scattering ***************
!        this is like elastic, both proj and targ mass
!        must be fixed
!***************************************************
  100 nfp5 = max(2, nfp(jp,5))
  nft5 = max(2, nft(jt,5))
  bb1 = 1.0 + d1 - d2
  bb2 = 1.0 + d2 - d1
  If (bb1**2<4.0*d1 .Or. bb2**2<4.0*d2) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  If (ranart(nseed)<0.5) Then
    x1 = (bb1-sqrt(bb1**2-4.0*d1))/2.0
    x2 = (bb2-sqrt(bb2**2-4.0*d2))/2.0
  Else
    x1 = (bb1+sqrt(bb1**2-4.0*d1))/2.0
    x2 = (bb2+sqrt(bb2**2-4.0*d2))/2.0
  End If
  ihnt2(13) = 2
  Goto 600
!
!********** Single diffractive ***********************
! either proj or targ's mass is fixed
!*****************************************************
  220 nfp5 = max(2, nfp(jp,5))
  nft5 = 3
  If (nfp3==0) nfp5 = 3
  bb2 = 1.0 + d2 - d1
  If (bb2**2<4.0*d2) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  xmin = (bb2-sqrt(bb2**2-4.0*d2))/2.0
  xmax = (bb2+sqrt(bb2**2-4.0*d2))/2.0
  miss4 = 0
  222 x2 = hirnd2(6, xmin, xmax)
  x1 = d1/(1.0-x2)
  If (x2*(1.0-x1)<(d2+1.E-4/sw)) Then
    miss4 = miss4 + 1
    If (miss4<=1000) Goto 222
    Goto 5000
  End If
  ihnt2(13) = 2
  Goto 600
!                        ********Fix proj mass*********
  240 nfp5 = 3
  nft5 = max(2, nft(jt,5))
  If (nft3==0) nft5 = 3
  bb1 = 1.0 + d1 - d2
  If (bb1**2<4.0*d1) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  xmin = (bb1-sqrt(bb1**2-4.0*d1))/2.0
  xmax = (bb1+sqrt(bb1**2-4.0*d1))/2.0
  miss4 = 0
  242 x1 = hirnd2(6, xmin, xmax)
  x2 = d2/(1.0-x1)
  If (x1*(1.0-x2)<(d1+1.E-4/sw)) Then
    miss4 = miss4 + 1
    If (miss4<=1000) Goto 242
    Goto 5000
  End If
  ihnt2(13) = 2
  Goto 600
!                        ********Fix targ mass*********
!
!*************non-single diffractive**********************
!        both proj and targ may not be fixed in mass
!*********************************************************
!
  400 nfp5 = 3
  nft5 = 3
  bb1 = 1.0 + d1 - d2
  bb2 = 1.0 + d2 - d1
  If (bb1**2<4.0*d1 .Or. bb2**2<4.0*d2) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 3000
    pkc = pkc*0.5
    Goto 30
  End If
  xmin1 = (bb1-sqrt(bb1**2-4.0*d1))/2.0
  xmax1 = (bb1+sqrt(bb1**2-4.0*d1))/2.0
  xmin2 = (bb2-sqrt(bb2**2-4.0*d2))/2.0
  xmax2 = (bb2+sqrt(bb2**2-4.0*d2))/2.0
  miss4 = 0
  410 x1 = hirnd2(4, xmin1, xmax1)
  x2 = hirnd2(4, xmin2, xmax2)
  If (nfp(jp,5)==3 .Or. nft(jt,5)==3) Then
    x1 = hirnd2(6, xmin1, xmax1)
    x2 = hirnd2(6, xmin2, xmax2)
  End If
!                        ********
  If (abs(nfp(jp,1)*nfp(jp,2))>1000000 .Or. abs(nfp(jp,1)*nfp(jp,2))<100) Then
    x1 = hirnd2(5, xmin1, xmax1)
  End If
  If (abs(nft(jt,1)*nft(jt,2))>1000000 .Or. abs(nft(jt,1)*nft(jt,2))<100) Then
    x2 = hirnd2(5, xmin2, xmax2)
  End If
!        IF(IOPMAIN.EQ.3) X1=HIRND2(6,XMIN1,XMAX1)
!        IF(IOPMAIN.EQ.2) X2=HIRND2(6,XMIN2,XMAX2)
!        ********For q-qbar or (qq)-(qq)bar system use symetric
!                distribution, for q-(qq) or qbar-(qq)bar use
!                unsymetrical distribution
!
  If (abs(nfp(jp,1)*nfp(jp,2))>1000000) x1 = 1.0 - x1
  xxp = x1*(1.0-x2)
  xxt = x2*(1.0-x1)
  If (xxp<(d1+1.E-4/sw) .Or. xxt<(d2+1.E-4/sw)) Then
    miss4 = miss4 + 1
    If (miss4<=1000) Goto 410
    Goto 5000
  End If
  ihnt2(13) = 3
!***************************************************
!***************************************************
  600 Continue
  If (x1*(1.0-x2)<(ampn**2-1.E-4)/sw .Or. x2*(1.0-x1)<(amtn**2-1.E-4)/sw) Then
    miss = miss + 1
    If (miss>100 .Or. pkc==0.0) Goto 2000
    pkc = 0.0
    Goto 30
  End If
!
  epp = (1.0-x2)*wp
  epm = x1*wm
  etp = x2*wp
  etm = (1.0-x1)*wm
  pp(jp, 3) = (epp-epm)/2.0
  pp(jp, 4) = (epp+epm)/2.0
  If (epp*epm-ptp2<0.0) Goto 6000
  pp(jp, 5) = sqrt(epp*epm-ptp2)
  nfp(jp, 3) = nfp3
  nfp(jp, 5) = nfp5

  pt(jt, 3) = (etp-etm)/2.0
  pt(jt, 4) = (etp+etm)/2.0
  If (etp*etm-ptt2<0.0) Goto 6000
  pt(jt, 5) = sqrt(etp*etm-ptt2)
  nft(jt, 3) = nft3
  nft(jt, 5) = nft5
!*****recoil PT from hard-inter is shared by two end-partons
!       so that pt=p1+p2
  pp(jp, 1) = pp11 - pkc11
  pp(jp, 2) = pp12 - pkc12

  kcdip = 1
  kcdit = 1
  If (abs(nfp(jp,1)*nfp(jp,2))>1000000 .Or. abs(nfp(jp,1)*nfp(jp,2))<100) Then
    kcdip = 0
  End If
  If (abs(nft(jt,1)*nft(jt,2))>1000000 .Or. abs(nft(jt,1)*nft(jt,2))<100) Then
    kcdit = 0
  End If
  If ((kcdip==0 .And. ranart(nseed)<0.5) .Or. (kcdip/=0 .And. ranart(nseed)<0.5/(1.0+(pkc11**2+pkc12**2)/hipr1(22)**2))) Then
    pp(jp, 6) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 6)
    pp(jp, 7) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 7)
    pp(jp, 8) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 8) + pkc11
    pp(jp, 9) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 9) + pkc12
  Else
    pp(jp, 8) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 8)
    pp(jp, 9) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 9)
    pp(jp, 6) = (pp(jp,1)-pp(jp,6)-pp(jp,8)-dpkc1)/2.0 + pp(jp, 6) + pkc11
    pp(jp, 7) = (pp(jp,2)-pp(jp,7)-pp(jp,9)-dpkc2)/2.0 + pp(jp, 7) + pkc12
  End If
  pp(jp, 1) = pp(jp, 6) + pp(jp, 8)
  pp(jp, 2) = pp(jp, 7) + pp(jp, 9)
!                                ********pt kick for proj
  pt(jt, 1) = pt11 - pkc21
  pt(jt, 2) = pt12 - pkc22
  If ((kcdit==0 .And. ranart(nseed)<0.5) .Or. (kcdit/=0 .And. ranart(nseed)<0.5/(1.0+(pkc21**2+pkc22**2)/hipr1(22)**2))) Then
    pt(jt, 6) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 6)
    pt(jt, 7) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 7)
    pt(jt, 8) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 8) + pkc21
    pt(jt, 9) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 9) + pkc22
  Else
    pt(jt, 8) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 8)
    pt(jt, 9) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 9)
    pt(jt, 6) = (pt(jt,1)-pt(jt,6)-pt(jt,8)-dpkc1)/2.0 + pt(jt, 6) + pkc21
    pt(jt, 7) = (pt(jt,2)-pt(jt,7)-pt(jt,9)-dpkc2)/2.0 + pt(jt, 7) + pkc22
  End If
  pt(jt, 1) = pt(jt, 6) + pt(jt, 8)
  pt(jt, 2) = pt(jt, 7) + pt(jt, 9)
!                        ********pt kick for targ

  If (npj(jp)/=0) nfp(jp, 5) = 3
  If (ntj(jt)/=0) nft(jt, 5) = 3
!                        ********jets must be connected to string
  If (epp/(epm+0.0001)<etp/(etm+0.0001) .And. abs(nfp(jp,1)*nfp(jp,2))<1000000) Then
    Do jsb = 1, 15
      psb = pp(jp, jsb)
      pp(jp, jsb) = pt(jt, jsb)
      pt(jt, jsb) = psb
      nsb = nfp(jp, jsb)
      nfp(jp, jsb) = nft(jt, jsb)
      nft(jt, jsb) = nsb
    End Do
!                ********when Ycm(JP)<Ycm(JT) after the collision
!                        exchange the positions of the two
  End If
!
  Return
!**************************************************
!**************************************************
  1000 ierror = 1
  If (ihpr2(10)==0) Return
  Write (6, *) '     Fatal HIJSFT start error,abandon this event'
  Write (6, *) '     PROJ E+,E-,W+', epp, epm, wp
  Write (6, *) '     TARG E+,E-,W-', etp, etm, wm
  Write (6, *) '     W+*W-, (APN+ATN)^2', sw, snn
  Return
  2000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     (2)energy partition fail,'
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     MP1,MPN', x1*(1.0-x2)*sw, ampn**2
  Write (6, *) '     MT2,MTN', x2*(1.0-x1)*sw, amtn**2
  Return
  3000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     (3)something is wrong with the pt kick, '
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     D1=', d1, ' D2=', d2, ' SW=', sw
  Write (6, *) '     HISTORY NFP5=', nfp(jp, 5), ' NFT5=', nft(jt, 5)
  Write (6, *) '     THIS COLLISON NFP5=', nfp5, ' NFT5=', nft5
  Write (6, *) '     # OF JET IN PROJ', npj(jp), ' IN TARG', ntj(jt)
  Return
  4000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     (4)unable to choose process, but not harmful'
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     PTP=', sqrt(ptp2), ' PTT=', sqrt(ptt2), ' SW=', sw
  Write (6, *) '     AMCUT=', amx, ' JP=', jp, ' JT=', jt
  Write (6, *) '     HISTORY NFP5=', nfp(jp, 5), ' NFT5=', nft(jt, 5)
  Return
  5000 ierror = 0
  If (ihpr2(10)==0) Return
  Write (6, *) '     energy partition failed(5),for limited try'
  Write (6, *) '     HIJSFT not performed, but continue'
  Write (6, *) '     NFP5=', nfp5, ' NFT5=', nft5
  Write (6, *) '     D1', d1, ' X1(1-X2)', x1*(1.0-x2)
  Write (6, *) '     D2', d2, ' X2(1-X1)', x2*(1.0-x1)
  Return
  6000 pkc = 0.0
  miss = miss + 1
  If (miss<100) Goto 30
  ierror = 1
  If (ihpr2(10)==0) Return
  Write (6, *) ' ERROR OCCURED, HIJSFT NOT PERFORMED'
  Write (6, *) ' Abort this event'
  Write (6, *) 'MTP,PTP2', epp*epm, ptp2, '  MTT,PTT2', etp*etm, ptt2
  Return
End Subroutine hijsft
!
!
!
! ********************************************************
! ************************              WOOD-SAX
Subroutine hijwds(ia, idh, xhigh)
!     SETS UP HISTOGRAM IDH WITH RADII FOR
!     NUCLEUS IA DISTRIBUTED ACCORDING TO THREE PARAM WOOD SAXON
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /wood/r, d, fnorm, w
!c      SAVE /WOOD/
!        DIMENSION IAA(20),RR(20),DD(20),WW(20),RMS(20)
  Dimension iaa(20), rr(20), dd(20), ww(20)
  External rwdsax, wdsax
  Save
!
!   PARAMETERS OF SPECIAL NUCLEI FROM ATOMIC DATA AND NUC DATA TABLES
!     VOL 14, 5-6 1974
  Data iaa/2, 4, 12, 16, 27, 32, 40, 56, 63, 93, 184, 197, 208, 7*0./
  Data rr/0.01, .964, 2.355, 2.608, 2.84, 3.458, 3.766, 3.971, 4.214, 4.87, 6.51, 6.38, 6.624, 7*0./
  Data dd/0.5882, .322, .522, .513, .569, .61, .586, .5935, .586, .573, .535, .535, .549, 7*0./
  Data ww/0.0, .517, -0.149, -0.051, 0., -0.208, -0.161, 13*0./
!        DATA RMS/2.11,1.71,2.46,2.73,3.05,3.247,3.482,3.737,3.925,4.31,
!     1        5.42,5.33,5.521,7*0./
!
  a = ia
!
!                 ********SET WOOD-SAX PARAMS FIRST  AS IN DATE ET AL
  d = 0.54
!                        ********D IS WOOD SAX DIFFUSE PARAM IN FM
  r = 1.19*a**(1./3.) - 1.61*a**(-1./3.)
!                         ********R IS RADIUS PARAM
  w = 0.
!                 ********W IS The third of three WOOD-SAX PARAM
!
!                      ********CHECK TABLE FOR SPECIAL CASES
  Do i = 1, 13
    If (ia==iaa(i)) Then
      r = rr(i)
      d = dd(i)
      w = ww(i)
!lin RS not used                              RS=RMS(I)
    End If
  End Do
!                             ********FNORM is the normalize factor
  fnorm = 1.0
  xlow = 0.
  xhigh = r + 12.*d
  If (w<-0.01) Then
    If (xhigh>r/sqrt(abs(w))) xhigh = r/sqrt(abs(w))
  End If
  fgaus = gauss1(rwdsax, xlow, xhigh, 0.001)
  fnorm = 1./fgaus
!
  If (idh==1) Then
    hint1(72) = r
    hint1(73) = d
    hint1(74) = w
    hint1(75) = fnorm/4.0/hipr1(40)
  Else If (idh==2) Then
    hint1(76) = r
    hint1(77) = d
    hint1(78) = w
    hint1(79) = fnorm/4.0/hipr1(40)
  End If
!
!             NOW SET UP HBOOK FUNCTIONS IDH FOR  R**2*RHO(R)
!             THESE HISTOGRAMS ARE USED TO GENERATE RANDOM RADII
  Call hifun(idh, xlow, xhigh, rwdsax)
  Return
End Subroutine hijwds
!
!
Function wdsax(x)
!                             ********THREE PARAMETER WOOD SAXON
  Common /wood/r, d, fnorm, w
!c      SAVE /WOOD/
  Save
  wdsax = fnorm*(1.+w*(x/r)**2)/(1+exp((x-r)/d))
  If (w<0.) Then
    If (x>=r/sqrt(abs(w))) wdsax = 0.
  End If
  Return
End Function wdsax
!
!
Function rwdsax(x)
  Save
  rwdsax = x*x*wdsax(x)
  Return
End Function rwdsax
!
!
!
!
! The next three subroutines are for Monte Carlo generation
! according to a given function FHB. One calls first HIFUN
! with assigned channel number I, low and up limits. Then to
! generate the distribution one can call HIRND(I) which gives
! you a random number generated according to the given function.
!
Subroutine hifun(i, xmin, xmax, fhb)
  Common /hijhb/rr(10, 201), xx(10, 201)
!c      SAVE /HIJHB/
  External fhb
  Save
  fnorm = gauss1(fhb, xmin, xmax, 0.001)
  Do j = 1, 201
    xx(i, j) = xmin + (xmax-xmin)*(j-1)/200.0
    xdd = xx(i, j)
    rr(i, j) = gauss1(fhb, xmin, xdd, 0.001)/fnorm
  End Do
  Return
End Subroutine hifun
!
!
!
Function hirnd(i)
  Common /hijhb/rr(10, 201), xx(10, 201)
!c      SAVE /HIJHB/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
  rx = ranart(nseed)
  jl = 0
  ju = 202
  10 If (ju-jl>1) Then
    jm = (ju+jl)/2
    If ((rr(i,201)>rr(i,1)) .Eqv. (rx>rr(i,jm))) Then
      jl = jm
    Else
      ju = jm
    End If
    Goto 10
  End If
  j = jl
  If (j<1) j = 1
  If (j>=201) j = 200
  hirnd = (xx(i,j)+xx(i,j+1))/2.0
  Return
End Function hirnd
!
!
!
!
!        This generate random number between XMIN and XMAX
Function hirnd2(i, xmin, xmax)
  Common /hijhb/rr(10, 201), xx(10, 201)
!c      SAVE /HIJHB/
  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save
  If (xmin<xx(i,1)) xmin = xx(i, 1)
  If (xmax>xx(i,201)) xmax = xx(i, 201)
  jmin = 1 + int(200*(xmin-xx(i,1))/(xx(i,201)-xx(i,1)))
  jmax = 1 + int(200*(xmax-xx(i,1))/(xx(i,201)-xx(i,1)))
  rx = rr(i, jmin) + (rr(i,jmax)-rr(i,jmin))*ranart(nseed)
  jl = 0
  ju = 202
  10 If (ju-jl>1) Then
    jm = (ju+jl)/2
    If ((rr(i,201)>rr(i,1)) .Eqv. (rx>rr(i,jm))) Then
      jl = jm
    Else
      ju = jm
    End If
    Goto 10
  End If
  j = jl
  If (j<1) j = 1
  If (j>=201) j = 200
  hirnd2 = (xx(i,j)+xx(i,j+1))/2.0
  Return
End Function hirnd2
!
!
!
!
Subroutine hijcrs
!        THIS IS TO CALCULATE THE CROSS SECTIONS OF JET PRODUCTION AND
!        THE TOTAL INELASTIC CROSS SECTIONS.
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /njet/n, ipcrs
!c      SAVE /NJET/
  External fhin, ftot, fnjet, ftotjt, ftotrg
  Save
  If (hint1(1)>=10.0) Call crsjet
!                        ********calculate jet cross section(in mb)
!
!lin-7/2009 these are related to nuclear shadowing:
  aphx1 = hipr1(6)*(ihnt2(1)**0.3333333-1.0)
  aphx2 = hipr1(6)*(ihnt2(3)**0.3333333-1.0)
  hint1(11) = hint1(14) - aphx1*hint1(15) - aphx2*hint1(16) + aphx1*aphx2*hint1(17)
  hint1(10) = gauss1(ftotjt, 0.0, 20.0, 0.01)
  hint1(12) = gauss1(fhin, 0.0, 20.0, 0.01)
  hint1(13) = gauss1(ftot, 0.0, 20.0, 0.01)
  hint1(60) = hint1(61) - aphx1*hint1(62) - aphx2*hint1(63) + aphx1*aphx2*hint1(64)
  hint1(59) = gauss1(ftotrg, 0.0, 20.0, 0.01)
  If (hint1(59)==0.0) hint1(59) = hint1(60)
  If (hint1(1)>=10.0) Then
    Do i = 0, 20
      n = i
      hint1(80+i) = gauss1(fnjet, 0.0, 20.0, 0.01)/hint1(12)
    End Do
  End If
  hint1(10) = hint1(10)*hipr1(31)
  hint1(12) = hint1(12)*hipr1(31)
  hint1(13) = hint1(13)*hipr1(31)
  hint1(59) = hint1(59)*hipr1(31)
!                ********Total and Inel cross section are calculated
!                        by Gaussian integration.
  If (ihpr2(13)/=0) Then
    hipr1(33) = 1.36*(1.0+36.0/hint1(1)**2)*alog(0.6+0.1*hint1(1)**2)
    hipr1(33) = hipr1(33)/hint1(12)
  End If
!                ********Parametrized cross section for single
!                        diffractive reaction(Goulianos)
  Return
End Subroutine hijcrs
!
!
!
!
Function ftot(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  omg = omg0(x)*(hipr1(30)+hint1(11))/hipr1(31)/2.0
  ftot = 2.0*(1.0-exp(-omg))
  Return
End Function ftot
!
!
!
Function fhin(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  omg = omg0(x)*(hipr1(30)+hint1(11))/hipr1(31)/2.0
  fhin = 1.0 - exp(-2.0*omg)
  Return
End Function fhin
!
!
!
Function ftotjt(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  omg = omg0(x)*hint1(11)/hipr1(31)/2.0
  ftotjt = 1.0 - exp(-2.0*omg)
  Return
End Function ftotjt
!
!
!
Function ftotrg(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Save
  omg = omg0(x)*hint1(60)/hipr1(31)/2.0
  ftotrg = 1.0 - exp(-2.0*omg)
  Return
End Function ftotrg
!
!
!
!
Function fnjet(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /njet/n, ipcrs
!c      SAVE /NJET/
  Save
  omg1 = omg0(x)*hint1(11)/hipr1(31)
!lin-8/2015 could cause IEEE_UNDERFLOW, does not seem to affect results:
  c0 = exp(n*alog(omg1)-sgmin(n+1))
  If (n==0) c0 = 1.0 - exp(-2.0*omg0(x)*hipr1(30)/hipr1(31)/2.0)
  fnjet = c0*exp(-omg1)
  Return
End Function fnjet
!
!
!
!
!
Function sgmin(n)
  Save
  ga = 0.
  If (n<=2) Goto 20
  Do i = 1, n - 1
    z = i
    ga = ga + alog(z)
  End Do
  20 sgmin = ga
  Return
End Function sgmin
!
!
!
Function omg0(x)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /besel/x4
!c      SAVE /BESEL/
  External bk
  Save
  x4 = hipr1(32)*sqrt(x)
  omg0 = hipr1(32)**2*gauss2(bk, x4, x4+20.0, 0.01)/96.0
  Return
End Function omg0
!
!
!
Function romg(x)
!                ********This gives the eikonal function from a table
!                        calculated in the first call
  Dimension fr(0:1000)
!lin-10/29/02 unsaved FR causes wrong values for ROMG with f77 compiler:
!c        SAVE FR
  Save
  Data i0/0/

  If (i0/=0) Goto 100
  Do i = 1, 1001
    xr = (i-1)*0.01
    fr(i-1) = omg0(xr)
  End Do
  100 i0 = 1
  If (x>=10.0) Then
    romg = 0.0
    Return
  End If
  ix = int(x*100)
  romg = (fr(ix)*((ix+1)*0.01-x)+fr(ix+1)*(x-ix*0.01))/0.01
  Return
End Function romg
!
!
!
Function bk(x)
  Common /besel/x4
!c      SAVE /BESEL/
  Save
  bk = exp(-x)*(x**2-x4**2)**2.50/15.0
  Return
End Function bk
!
!
!        THIS PROGRAM IS TO CALCULATE THE JET CROSS SECTION
!        THE INTEGRATION IS DONE BY USING VEGAS
!
Subroutine crsjet
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
!c      SAVE /HPARNT/
  Common /njet/n, ipcrs
!c      SAVE /NJET/
  Common /bveg1/xl(10), xu(10), acc, ndim, ncall, itmx, nprn
!c      SAVE /BVEG1/
  Common /bveg2/xi(50, 10), si, si2, swgt, schi, ndo, it
!c      SAVE /BVEG2/
  Common /bveg3/f, ti, tsi
!c      SAVE /BVEG3/
  Common /sedvax/num1
!c      SAVE /SEDVAX/
  External fjet, fjetrg
  Save
!
!************************
!        NCALL give the number of inner-iteration, ITMX
!       gives the limit of out-iteration. Nprn is an option
!       ( 1: print the integration process. 0: do not print)
!
  ndim = 3
  ipcrs = 0
  Call vegas(fjet, avgi, sd, chi2a)
  hint1(14) = sngl(avgi)/2.5682
  If (ihpr2(6)==1 .And. ihnt2(1)>1) Then
    ipcrs = 1
    Call vegas(fjet, avgi, sd, chi2a)
    hint1(15) = sngl(avgi)/2.5682
  End If
  If (ihpr2(6)==1 .And. ihnt2(3)>1) Then
    ipcrs = 2
    Call vegas(fjet, avgi, sd, chi2a)
    hint1(16) = sngl(avgi)/2.5682
  End If
  If (ihpr2(6)==1 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
    ipcrs = 3
    Call vegas(fjet, avgi, sd, chi2a)
    hint1(17) = sngl(avgi)/2.5682
  End If
!                ********Total inclusive jet cross section(Pt>P0)
!
  If (ihpr2(3)/=0) Then
    ipcrs = 0
    Call vegas(fjetrg, avgi, sd, chi2a)
    hint1(61) = sngl(avgi)/2.5682
    If (ihpr2(6)==1 .And. ihnt2(1)>1) Then
      ipcrs = 1
      Call vegas(fjetrg, avgi, sd, chi2a)
      hint1(62) = sngl(avgi)/2.5682
    End If
    If (ihpr2(6)==1 .And. ihnt2(3)>1) Then
      ipcrs = 2
      Call vegas(fjetrg, avgi, sd, chi2a)
      hint1(63) = sngl(avgi)/2.5682
    End If
    If (ihpr2(6)==1 .And. ihnt2(1)>1 .And. ihnt2(3)>1) Then
      ipcrs = 3
      Call vegas(fjetrg, avgi, sd, chi2a)
      hint1(64) = sngl(avgi)/2.5682
    End If
  End If
!                        ********cross section of trigger jet
!
  Return
End Subroutine crsjet
!
!
!
Function fjet(x, wgt)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
!c      SAVE /HPARNT/
  Dimension x(10)
  Save
  pt2 = dble(hint1(1)**2/4.0-hipr1(8)**2)*x(1) + dble(hipr1(8))**2
  xt = 2.0D0*dsqrt(pt2)/dble(hint1(1))
  ymx1 = dlog(1.0D0/xt+dsqrt(1.0D0/xt**2-1.0D0))
  y1 = 2.0D0*ymx1*x(2) - ymx1
  ymx2 = dlog(2.0D0/xt-dexp(y1))
  ymn2 = dlog(2.0D0/xt-dexp(-y1))
  y2 = (ymx2+ymn2)*x(3) - ymn2
  fjet = 2.0D0*ymx1*(ymx2+ymn2)*dble(hint1(1)**2/4.0-hipr1(8)**2)*g(y1, y2, pt2)/2.0D0
  Return
End Function fjet
!
!
!
Function fjetrg(x, wgt)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100), ptmax, ptmin
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
!c      SAVE /HPARNT/
  Dimension x(10)
  Save
  ptmin = abs(hipr1(10)) - 0.25
  ptmin = max(ptmin, hipr1(8))
  am2 = 0.D0
  If (ihpr2(3)==3) Then
    am2 = dble(hipr1(7)**2)
    ptmin = max(0.0, hipr1(10))
  End If
  ptmax = abs(hipr1(10)) + 0.25
  If (hipr1(10)<=0.0) ptmax = hint1(1)/2.0 - sngl(am2)
  If (ptmax<=ptmin) ptmax = ptmin + 0.25
  pt2 = dble(ptmax**2-ptmin**2)*x(1) + dble(ptmin)**2
  amt2 = pt2 + am2
  xt = 2.0D0*dsqrt(amt2)/dble(hint1(1))
  ymx1 = dlog(1.0D0/xt+dsqrt(1.0D0/xt**2-1.0D0))
  y1 = 2.0D0*ymx1*x(2) - ymx1
  ymx2 = dlog(2.0D0/xt-dexp(y1))
  ymn2 = dlog(2.0D0/xt-dexp(-y1))
  y2 = (ymx2+ymn2)*x(3) - ymn2
  If (ihpr2(3)==3) Then
    gtrig = 2.0D0*ghvq(y1, y2, amt2)
  Else If (ihpr2(3)==2) Then
    gtrig = 2.0D0*gphotn(y1, y2, pt2)
  Else
    gtrig = g(y1, y2, pt2)
  End If
  fjetrg = 2.0D0*ymx1*(ymx2+ymn2)*dble(ptmax**2-ptmin**2)*gtrig/2.0D0
  Return
End Function fjetrg
!
!
!
Function ghvq(y1, y2, amt2)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
!c      SAVE /HPARNT/
  Dimension f(2, 7)
  Save
  xt = 2.0D0*dsqrt(amt2)/dble(hint1(1))
  x1 = 0.5D0*xt*(dexp(y1)+dexp(y2))
  x2 = 0.5D0*xt*(dexp(-y1)+dexp(-y2))
  ss = x1*x2*dble(hint1(1))**2
  af = 4.0D0
  If (ihpr2(18)/=0) af = 5.0D0
  dlam = dble(hipr1(15))
  aph = 12.0D0*3.1415926D0/(33.D0-2.D0*af)/dlog(amt2/dlam**2)
!
  Call parton(f, x1, x2, amt2)
!
  gqq = 4.D0*(dcosh(y1-y2)+dble(hipr1(7))**2/amt2)/(1.D0+dcosh(y1-y2))/9.D0*(f(1,1)*f(2,2)+f(1,2)*f(2,1)+f(1,3)*f(2,4)+f(1,4)*f(2,3)+f(1,5)*f(2,6)+f(1,6)*f(2,5))
  ggg = (8.D0*dcosh(y1-y2)-1.D0)*(dcosh(y1-y2)+2.D0*dble(hipr1(7))**2/amt2-2.D0*dble(hipr1(7))**4/amt2**2)/(1.D0+dcosh(y1-y2))/24.D0*f(1, 7)*f(2, 7)
!
  ghvq = (gqq+ggg)*dble(hipr1(23))*3.14159D0*aph**2/ss**2
  Return
End Function ghvq
!
!
!
Function gphotn(y1, y2, pt2)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
!c      SAVE /HPARNT/
  Dimension f(2, 7)
  Save
  xt = 2.D0*dsqrt(pt2)/dble(hint1(1))
  x1 = 0.5D0*xt*(dexp(y1)+dexp(y2))
  x2 = 0.5D0*xt*(dexp(-y1)+dexp(-y2))
  z = dsqrt(1.D0-xt**2/x1/x2)
  ss = x1*x2*dble(hint1(1))**2
  t = -(1.D0-z)/2.D0
  u = -(1.D0+z)/2.D0
  af = 3.D0
  dlam = dble(hipr1(15))
  aph = 12.D0*3.1415926D0/(33.D0-2.D0*af)/dlog(pt2/dlam**2)
  aphem = 1.D0/137.D0
!
  Call parton(f, x1, x2, pt2)
!
  g11 = -(u**2+1.D0)/u/3.D0*f(1, 7)*(4.D0*f(2,1)+4.D0*f(2,2)+f(2,3)+f(2,4)+f(2,5)+f(2,6))/9.D0
  g12 = -(t**2+1.D0)/t/3.D0*f(2, 7)*(4.D0*f(1,1)+4.D0*f(1,2)+f(1,3)+f(1,4)+f(1,5)+f(1,6))/9.D0
  g2 = 8.D0*(u**2+t**2)/u/t/9.D0*(4.D0*f(1,1)*f(2,2)+4.D0*f(1,2)*f(2,1)+f(1,3)*f(2,4)+f(1,4)*f(2,3)+f(1,5)*f(2,6)+f(1,6)*f(2,5))/9.D0
!
  gphotn = (g11+g12+g2)*dble(hipr1(23))*3.14159D0*aph*aphem/ss**2
  Return
End Function gphotn
!
!
!
!
Function g(y1, y2, pt2)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
!c      SAVE /HPARNT/
  Dimension f(2, 7)
  Save
  xt = 2.D0*dsqrt(pt2)/dble(hint1(1))
  x1 = 0.5D0*xt*(dexp(y1)+dexp(y2))
  x2 = 0.5D0*xt*(dexp(-y1)+dexp(-y2))
  z = dsqrt(1.D0-xt**2/x1/x2)
  ss = x1*x2*dble(hint1(1))**2
  t = -(1.D0-z)/2.D0
  u = -(1.D0+z)/2.D0
  af = 3.D0
  dlam = dble(hipr1(15))
  aph = 12.D0*3.1415926D0/(33.D0-2.D0*af)/dlog(pt2/dlam**2)
!
  Call parton(f, x1, x2, pt2)
!
  g11 = ((f(1,1)+f(1,2))*(f(2,3)+f(2,4)+f(2,5)+f(2,6))+(f(1,3)+f(1,4))*(f(2,5)+f(2,6)))*subcr1(t, u)
!
  g12 = ((f(2,1)+f(2,2))*(f(1,3)+f(1,4)+f(1,5)+f(1,6))+(f(2,3)+f(2,4))*(f(1,5)+f(1,6)))*subcr1(u, t)
!
  g13 = (f(1,1)*f(2,1)+f(1,2)*f(2,2)+f(1,3)*f(2,3)+f(1,4)*f(2,4)+f(1,5)*f(2,5)+f(1,6)*f(2,6))*(subcr1(u,t)+subcr1(t,u)-8.D0/t/u/27.D0)
!
  g2 = (af-1)*(f(1,1)*f(2,2)+f(2,1)*f(1,2)+f(1,3)*f(2,4)+f(2,3)*f(1,4)+f(1,5)*f(2,6)+f(2,5)*f(1,6))*subcr2(t, u)
!
  g31 = (f(1,1)*f(2,2)+f(1,3)*f(2,4)+f(1,5)*f(2,6))*subcr3(t, u)
  g32 = (f(2,1)*f(1,2)+f(2,3)*f(1,4)+f(2,5)*f(1,6))*subcr3(u, t)
!
  g4 = (f(1,1)*f(2,2)+f(2,1)*f(1,2)+f(1,3)*f(2,4)+f(2,3)*f(1,4)+f(1,5)*f(2,6)+f(2,5)*f(1,6))*subcr4(t, u)
!
  g5 = af*f(1, 7)*f(2, 7)*subcr5(t, u)
!
  g61 = f(1, 7)*(f(2,1)+f(2,2)+f(2,3)+f(2,4)+f(2,5)+f(2,6))*subcr6(t, u)
  g62 = f(2, 7)*(f(1,1)+f(1,2)+f(1,3)+f(1,4)+f(1,5)+f(1,6))*subcr6(u, t)
!
  g7 = f(1, 7)*f(2, 7)*subcr7(t, u)
!
  g = (g11+g12+g13+g2+g31+g32+g4+g5+g61+g62+g7)*dble(hipr1(17))*3.14159D0*aph**2/ss**2
  Return
End Function g
!
!
!
Function subcr1(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr1 = 4.D0/9.D0*(1.D0+u**2)/t**2
  Return
End Function subcr1
!
!
Function subcr2(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr2 = 4.D0/9.D0*(t**2+u**2)
  Return
End Function subcr2
!
!
Function subcr3(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr3 = 4.D0/9.D0*(t**2+u**2+(1.D0+u**2)/t**2-2.D0*u**2/3.D0/t)
  Return
End Function subcr3
!
!
Function subcr4(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr4 = 8.D0/3.D0*(t**2+u**2)*(4.D0/9.D0/t/u-1.D0)
  Return
End Function subcr4
!
!
!
Function subcr5(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr5 = 3.D0/8.D0*(t**2+u**2)*(4.D0/9.D0/t/u-1.D0)
  Return
End Function subcr5
!
!
Function subcr6(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr6 = (1.D0+u**2)*(1.D0/t**2-4.D0/u/9.D0)
  Return
End Function subcr6
!
!
Function subcr7(t, u)
  Implicit Double Precision (A-H, O-Z)
  subcr7 = 9.D0/2.D0*(3.D0-t*u-u/t**2-t/u**2)
  Return
End Function subcr7
!
!
!
Subroutine parton(f, x1, x2, qq)
  Implicit Double Precision (A-H, O-Z)
  Real hipr1(100), hint1(100)
  Common /hparnt/hipr1, ihpr2(50), hint1, ihnt2(50)
!c      SAVE /HPARNT/
  Common /njet/n, ipcrs
!c      SAVE /NJET/
!lin-7/2009:
  Common /cmsflag/dshadow, ishadow
  Dimension f(2, 7)
  Save
  dlam = dble(hipr1(15))
  q0 = dble(hipr1(16))
  s = dlog(dlog(qq/dlam**2)/dlog(q0**2/dlam**2))
  If (ihpr2(7)==2) Goto 200
!*******************************************************
  at1 = 0.419D0 + 0.004D0*s - 0.007D0*s**2
  at2 = 3.460D0 + 0.724D0*s - 0.066D0*s**2
  gmud = 4.40D0 - 4.86D0*s + 1.33D0*s**2
  at3 = 0.763D0 - 0.237D0*s + 0.026D0*s**2
  at4 = 4.00D0 + 0.627D0*s - 0.019D0*s**2
  gmd = -0.421D0*s + 0.033D0*s**2
!*******************************************************
  cas = 1.265D0 - 1.132D0*s + 0.293D0*s**2
  as = -0.372D0*s - 0.029D0*s**2
  bs = 8.05D0 + 1.59D0*s - 0.153D0*s**2
  aphs = 6.31D0*s - 0.273D0*s**2
  btas = -10.5D0*s - 3.17D0*s**2
  gms = 14.7D0*s + 9.80D0*s**2
!********************************************************
!        CAC=0.135*S-0.075*S**2
!        AC=-0.036-0.222*S-0.058*S**2
!        BC=6.35+3.26*S-0.909*S**2
!        APHC=-3.03*S+1.50*S**2
!        BTAC=17.4*S-11.3*S**2
!        GMC=-17.9*S+15.6*S**2
!***********************************************************
  cag = 1.56D0 - 1.71D0*s + 0.638D0*s**2
  ag = -0.949D0*s + 0.325D0*s**2
  bg = 6.0D0 + 1.44D0*s - 1.05D0*s**2
  aphg = 9.0D0 - 7.19D0*s + 0.255D0*s**2
  btag = -16.5D0*s + 10.9D0*s**2
  gmg = 15.3D0*s - 10.1D0*s**2
  Goto 300
!********************************************************
  200 at1 = 0.374D0 + 0.014D0*s
  at2 = 3.33D0 + 0.753D0*s - 0.076D0*s**2
  gmud = 6.03D0 - 6.22D0*s + 1.56D0*s**2
  at3 = 0.761D0 - 0.232D0*s + 0.023D0*s**2
  at4 = 3.83D0 + 0.627D0*s - 0.019D0*s**2
  gmd = -0.418D0*s + 0.036D0*s**2
!************************************
  cas = 1.67D0 - 1.92D0*s + 0.582D0*s**2
  as = -0.273D0*s - 0.164D0*s**2
  bs = 9.15D0 + 0.530D0*s - 0.763D0*s**2
  aphs = 15.7D0*s - 2.83D0*s**2
  btas = -101.0D0*s + 44.7D0*s**2
  gms = 223.0D0*s - 117.0D0*s**2
!*********************************
!        CAC=0.067*S-0.031*S**2
!        AC=-0.120-0.233*S-0.023*S**2
!        BC=3.51+3.66*S-0.453*S**2
!        APHC=-0.474*S+0.358*S**2
!        BTAC=9.50*S-5.43*S**2
!        GMC=-16.6*S+15.5*S**2
!**********************************
  cag = 0.879D0 - 0.971D0*s + 0.434D0*s**2
  ag = -1.16D0*s + 0.476D0*s**2
  bg = 4.0D0 + 1.23D0*s - 0.254D0*s**2
  aphg = 9.0D0 - 5.64D0*s - 0.817D0*s**2
  btag = -7.54D0*s + 5.50D0*s**2
  gmg = -0.596D0*s + 1.26D0*s**2
!*********************************
  300 b12 = dexp(gmre(at1)+gmre(at2+1.D0)-gmre(at1+at2+1.D0))
  b34 = dexp(gmre(at3)+gmre(at4+1.D0)-gmre(at3+at4+1.D0))
  cnud = 3.D0/b12/(1.D0+gmud*at1/(at1+at2+1.D0))
  cnd = 1.D0/b34/(1.D0+gmd*at3/(at3+at4+1.D0))
!********************************************************
!        FUD=X*(U+D)
!        FS=X*2(UBAR+DBAR+SBAR)  AND UBAR=DBAR=SBAR
!*******************************************************
  fud1 = cnud*x1**at1*(1.D0-x1)**at2*(1.D0+gmud*x1)
  fs1 = cas*x1**as*(1.D0-x1)**bs*(1.D0+aphs*x1+btas*x1**2+gms*x1**3)
  f(1, 3) = cnd*x1**at3*(1.D0-x1)**at4*(1.D0+gmd*x1) + fs1/6.D0
  f(1, 1) = fud1 - f(1, 3) + fs1/3.D0
  f(1, 2) = fs1/6.D0
  f(1, 4) = fs1/6.D0
  f(1, 5) = fs1/6.D0
  f(1, 6) = fs1/6.D0
  f(1, 7) = cag*x1**ag*(1.D0-x1)**bg*(1.D0+aphg*x1+btag*x1**2+gmg*x1**3)
!
  fud2 = cnud*x2**at1*(1.D0-x2)**at2*(1.D0+gmud*x2)
  fs2 = cas*x2**as*(1.D0-x2)**bs*(1.D0+aphs*x2+btas*x2**2+gms*x2**3)
  f(2, 3) = cnd*x2**at3*(1.D0-x2)**at4*(1.D0+gmd*x2) + fs2/6.D0
  f(2, 1) = fud2 - f(2, 3) + fs2/3.D0
  f(2, 2) = fs2/6.D0
  f(2, 4) = fs2/6.D0
  f(2, 5) = fs2/6.D0
  f(2, 6) = fs2/6.D0
  f(2, 7) = cag*x2**ag*(1.D0-x2)**bg*(1.D0+aphg*x2+btag*x2**2+gmg*x2**3)
!***********Nuclear effect on the structure function****************
!
  If (ihpr2(6)==1 .And. ihnt2(1)>1) Then
    aax = 1.193D0*dble(alog(float(ihnt2(1)))**0.16666666)
    rrx = aax*(x1**3-1.2D0*x1**2+0.21D0*x1) + 1.D0 + dble(1.079*(float(ihnt2(1))**0.33333333-1.0))/dble(alog(float(ihnt2(1))+1.0))*dsqrt(x1)*dexp(-x1**2/0.01D0)
!lin-8/2015 DEXP() above may cause IEEE_UNDERFLOW,
!     does not seem to affect results.
!     &          /DLOG(IHNT2(1)+1.0D0)*(DSQRT(X1)*DEXP(-X1**2/0.01)
!lin-7/2009 enable users to modify nuclear shadowing:
    If (ishadow==1) rrx = 1.D0 + dshadow*(rrx-1.D0)
    If (ipcrs==1 .Or. ipcrs==3) rrx = dexp(-x1**2/0.01D0)
!lin-7/2009:
    If ((ipcrs==1 .Or. ipcrs==3) .And. ishadow==1) rrx = dexp(-x1**2/0.01D0)*dshadow
    Do i = 1, 7
      f(1, i) = rrx*f(1, i)
    End Do
  End If
  If (ihpr2(6)==1 .And. ihnt2(3)>1) Then
    aax = 1.193D0*dble(alog(float(ihnt2(3)))**0.16666666)
    rrx = aax*(x2**3-1.2D0*x2**2+0.21D0*x2) + 1.D0 + dble(1.079*(float(ihnt2(3))**0.33333-1.0))/dble(alog(float(ihnt2(3))+1.0))*dsqrt(x2)*dexp(-x2**2/0.01D0)
!     &         /DLOG(IHNT2(3)+1.0D0)*DSQRT(X2)*DEXP(-X2**2/0.01)
!lin-7/2009:
    If (ishadow==1) rrx = 1.D0 + dshadow*(rrx-1.D0)
    If (ipcrs==2 .Or. ipcrs==3) rrx = dexp(-x2**2/0.01D0)
!lin-7/2009:
    If ((ipcrs==2 .Or. ipcrs==3) .And. ishadow==1) rrx = dexp(-x2**2/0.01D0)*dshadow
    Do i = 1, 7
      f(2, i) = rrx*f(2, i)
    End Do
  End If
!
  Return
End Subroutine parton
!
!
!
Function gmre(x)
  Implicit Double Precision (A-H, O-Z)
  Save
  z = x
  If (x>3.0D0) Goto 10
  z = x + 3.D0
  10 gmre = 0.5D0*dlog(2.D0*3.14159265D0/z) + z*dlog(z) - z + dlog(1.D0+1.D0/12.D0/z+1.D0/288.D0/z**2-139.D0/51840.D0/z**3-571.D0/2488320.D0/z**4)
  If (z==x) Goto 20
  gmre = gmre - dlog(z-1.D0) - dlog(z-2.D0) - dlog(z-3.D0)
  20 Continue
  Return
End Function gmre
!
!
!
!***************************************************************

Block Data hidata
  Parameter (maxstr=150001)
  Double Precision xl(10), xu(10), acc
  Common /bveg1/xl, xu, acc, ndim, ncall, itmx, nprn
!c      SAVE /BVEG1/
  Common /sedvax/num1
!c      SAVE /SEDVAX/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
!c      SAVE /HPARNT/
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
!c      SAVE /HMAIN1/
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
!c      SAVE /HMAIN2/
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
!c      SAVE /HSTRNG/
  Common /hjcrdn/yp(3, 300), yt(3, 300)
!c      SAVE /hjcrdn/
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
!c      SAVE /HJJET1/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
!c      SAVE /HJJET2/
  Common /hijdat/hidat0(10, 10), hidat(10)
!c      SAVE /HIJDAT/
  Common /hpint/mint4, mint5, atco(200, 20), atxs(0:200)
!c      SAVE /HPINT/
  Save
  Data num1/30123984/, xl/10*0.D0/, xu/10*1.D0/
  Data ncall/1000/, itmx/100/, acc/0.01/, nprn/0/
!...give all the switchs and parameters the default values
!lin-4/2008 input.ampt provides NSEED for AMPT:
!        DATA NSEED/74769375/
  Data hipr1/1.5, 0.35, 0.5, 0.9, 2.0, 0.1, 1.5, 2.0, -1.0, -2.25, 2.0, 0.5, 1.0, 2.0, 0.2, 2.0, 2.5, 0.3, 0.1, 1.4, 1.6, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 57.0, 28.5, 3.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.14159, 0.0, 0.4, 0.1, 1.5, 0.1, 0.25, 0.0, 0.5, 0.0, 0.0, 50*0.0/

  Data ihpr2/1, 3, 0, 1, 1, 1, 1, 10, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 30*0/

  Data hint1/100*0/
  Data ihnt2/50*0/

!...initialize all the data common blocks
  Data natt/0/, eatt/0.0/, jatt/0/, nt/0/, np/0/, n0/0/, n01/0/, n10/0/, n11/0/
!lin-4/26/01
!        DATA KATT/520000*0/PATT/520000*0.0/
  Data katt/600004*0/, patt/600004*0.0/

  Data nfp/4500*0/, pp/4500*0.0/, nft/4500*0/, pt/4500*0.0/

  Data yp/900*0.0/, yt/900*0.0/

  Data npj/300*0/, kfpj/150000*0/, pjpx/150000*0.0/, pjpy/150000*0.0/, pjpz/150000*0.0/, pjpe/150000*0.0/, pjpm/150000*0.0/
  Data ntj/300*0/, kftj/150000*0/, pjtx/150000*0.0/, pjty/150000*0.0/, pjtz/150000*0.0/, pjte/150000*0.0/, pjtm/150000*0.0/

!lin-4/2008
!        DATA NSG/0/,NJSG/900*0/,IASG/2700*0/,K1SG/90000*0/,K2SG/90000*0/
!     &       ,PXSG/90000*0.0/,PYSG/90000*0.0/,PZSG/90000*0.0/
!     &       ,PESG/90000*0.0/,PMSG/90000*0.0/
  Data nsg/0/, njsg/150001*0/, iasg/450003*0/, k1sg/15000100*0/, k2sg/15000100*0/, pxsg/15000100*0.0/, pysg/15000100*0.0/, pzsg/15000100*0.0/, pesg/15000100*0.0/, pmsg/15000100*0.0/
  Data mint4/0/, mint5/0/, atco/4000*0.0/, atxs/201*0.0/
  Data (hidat0(1,i), i=1, 10)/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.25, 2.5, 4.0, 4.1/
  Data (hidat0(2,i), i=1, 10)/2.0, 3.0, 5.0, 6.0, 7.0, 8.0, 8.0, 10.0, 10.0, 10.0/
  Data (hidat0(3,i), i=1, 10)/1.0, 0.8, 0.8, 0.7, 0.45, 0.215, 0.21, 0.19, 0.19, 0.19/
  Data (hidat0(4,i), i=1, 10)/0.35, 0.35, 0.3, 0.3, 0.3, 0.3, 0.5, 0.6, 0.6, 0.6/
  Data (hidat0(5,i), i=1, 10)/23.8, 24.0, 26.0, 26.2, 27.0, 28.5, 28.5, 28.5, 28.5, 28.5/
  Data ((hidat0(j,i),i=1,10), j=6, 9)/40*0.0/
  Data (hidat0(10,i), i=1, 10)/5.0, 20.0, 53.0, 62.0, 100.0, 200.0, 546.0, 900.0, 1800.0, 4000.0/
  Data hidat/10*0.0/
End Block Data hidata
!*******************************************************************
!
!
!
!
!*******************************************************************
!   SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG'N
!      - BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
!*******************************************************************
!
Subroutine vegas(fxn, avgi, sd, chi2a)
  Implicit Double Precision (A-H, O-Z)
  Common /bveg1/xl(10), xu(10), acc, ndim, ncall, itmx, nprn
!c      SAVE /BVEG1/
  Common /bveg2/xi(50, 10), si, si2, swgt, schi, ndo, it
!c      SAVE /BVEG2/
  Common /bveg3/f, ti, tsi
!c      SAVE /BVEG3/
  External fxn
  Dimension d(50, 10), di(50, 10), xin(50), r(50), dx(10), dt(10), x(10), kg(10), ia(10)
!      REAL*4 QRAN(10)
  Real qran(10)
  Save
  Data ndmx/50/, alph/1.5D0/, one/1.D0/, mds/ -1/
!
  ndo = 1
  Do j = 1, ndim
    xi(1, j) = one
  End Do
!
  Entry vegas1(fxn, avgi, sd, chi2a)
!         - INITIALIZES CUMMULATIVE VARIABLES, BUT NOT GRID
  it = 0
  si = 0.D0
  si2 = si
  swgt = si
  schi = si
!
  Entry vegas2(fxn, avgi, sd, chi2a)
!         - NO INITIALIZATION
  nd = ndmx
  ng = 1
  If (mds==0) Goto 2
  ng = int((real(ncall)/2.)**(1./real(ndim)))
  mds = 1
  If ((2*ng-ndmx)<0) Goto 2
  mds = -1
  npg = ng/ndmx + 1
  nd = ng/npg
  ng = npg*nd
  2 k = ng**ndim
  npg = ncall/k
  If (npg<2) npg = 2
  calls = npg*k
  dxg = one/ng
  dv2g = (calls*dxg**ndim)**2/npg/npg/(npg-one)
  xnd = nd
  ndm = nd - 1
  dxg = dxg*xnd
  xjac = one/calls
  Do j = 1, ndim
!***this is the line 50
    dx(j) = xu(j) - xl(j)
    xjac = xjac*dx(j)
  End Do
!
!   REBIN PRESERVING BIN DENSITY
!
  If (nd==ndo) Goto 8
  rc = ndo/xnd
  Do j = 1, ndim
    k = 0
    xn = 0.D0
    dr = xn
    i = k
    4 k = k + 1
    dr = dr + one
    xo = xn
    xn = xi(k, j)
    5 If (rc>dr) Goto 4
    i = i + 1
    dr = dr - rc
    xin(i) = xn - (xn-xo)*dr
    If (i<ndm) Goto 5
    Do i = 1, ndm
      xi(i, j) = xin(i)
    End Do
    xi(nd, j) = one
  End Do
  ndo = nd
!
  8 Continue
!      IF(NPRN.NE.0) WRITE(16,200) NDIM,CALLS,IT,ITMX,ACC,MDS,ND
!     1                           ,(XL(J),XU(J),J=1,NDIM)
!
  Entry vegas3(fxn, avgi, sd, chi2a)
!         - MAIN INTEGRATION LOOP
  9 it = it + 1
  ti = 0.D0
  tsi = ti
  Do j = 1, ndim
    kg(j) = 1
    Do i = 1, nd
      d(i, j) = ti
      di(i, j) = ti
    End Do
  End Do
!
  11 fb = 0.D0
  f2b = fb
  k = 0
  12 k = k + 1
  Call aran9(qran, ndim)
  wgt = xjac
  Do j = 1, ndim
    xn = dble(float(kg(j))-qran(j))*dxg + one
!*****this is the line 100
    ia(j) = int(xn)
    If (ia(j)>1) Goto 13
    xo = xi(ia(j), j)
    rc = (xn-ia(j))*xo
    Goto 14
    13 xo = xi(ia(j), j) - xi(ia(j)-1, j)
    rc = xi(ia(j)-1, j) + (xn-ia(j))*xo
    14 x(j) = xl(j) + rc*dx(j)
    wgt = wgt*xo*xnd
  End Do
!
  f = wgt
  f = f*fxn(x, wgt)
  f2 = f*f
  fb = fb + f
  f2b = f2b + f2
  Do j = 1, ndim
    di(ia(j), j) = di(ia(j), j) + f
    If (mds>=0) d(ia(j), j) = d(ia(j), j) + f2
  End Do
  If (k<npg) Goto 12
!
  f2b = dsqrt(f2b*npg)
  f2b = (f2b-fb)*(f2b+fb)
  ti = ti + fb
  tsi = tsi + f2b
  If (mds>=0) Goto 18
  Do j = 1, ndim
    d(ia(j), j) = d(ia(j), j) + f2b
  End Do
  18 k = ndim
  19 kg(k) = mod(kg(k), ng) + 1
  If (kg(k)/=1) Goto 11
  k = k - 1
  If (k>0) Goto 19
!
!   FINAL RESULTS FOR THIS ITERATION
!
  tsi = tsi*dv2g
  ti2 = ti*ti
  wgt = ti2/(tsi+1.0D-37)
  si = si + ti*wgt
  si2 = si2 + ti2
  swgt = swgt + wgt
  swgt = swgt + 1.0D-37
  si2 = si2 + 1.0D-37
  schi = schi + ti2*wgt
  avgi = si/swgt
  sd = swgt*it/si2
  chi2a = sd*(schi/swgt-avgi*avgi)/dble(float(it)-.999)
  sd = dsqrt(one/sd)
!****this is the line 150
  If (nprn==0) Goto 21
  tsi = dsqrt(tsi)
!      WRITE(16,201) IT,TI,TSI,AVGI,SD,CHI2A
!      IF(NPRN.GE.0) GO TO 21
!      DO 20 J=1,NDIM
!20    WRITE(16,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
!
!   REFINE GRID
!
  21 Do j = 1, ndim
    xo = d(1, j)
    xn = d(2, j)
    d(1, j) = (xo+xn)/2.D0
    dt(j) = d(1, j)
    Do i = 2, ndm
      d(i, j) = xo + xn
      xo = xn
      xn = d(i+1, j)
      d(i, j) = (d(i,j)+xn)/3.D0
      dt(j) = dt(j) + d(i, j)
    End Do
    d(nd, j) = (xn+xo)/2.D0
    dt(j) = dt(j) + d(nd, j)
  End Do
!
  Do j = 1, ndim
    rc = 0.D0
    Do i = 1, nd
      r(i) = 0.D0
      If (dt(j)>=1.0D18) Then
        Write (6, *) '************** A SINGULARITY >1.0D18'
!      WRITE(5,1111)
!1111  FORMAT(1X,'**************IMPORTANT NOTICE***************')
!      WRITE(5,1112)
!1112  FORMAT(1X,'THE INTEGRAND GIVES RISE A SINGULARITY >1.0D18')
!      WRITE(5,1113)
!1113  FORMAT(1X,'PLEASE CHECK THE INTEGRAND AND THE LIMITS')
!      WRITE(5,1114)
!1114  FORMAT(1X,'**************END NOTICE*************')
      End If
      If (d(i,j)<=1.0D-18) Goto 24
      xo = dt(j)/d(i, j)
      r(i) = ((xo-one)/xo/dlog(xo))**alph
      24 rc = rc + r(i)
    End Do
    rc = rc/xnd
    k = 0
    xn = 0.D0
    dr = xn
    i = k
    25 k = k + 1
    dr = dr + r(k)
    xo = xn
!****this is the line 200
    xn = xi(k, j)
    26 If (rc>dr) Goto 25
    i = i + 1
    dr = dr - rc
    xin(i) = xn - (xn-xo)*dr/(r(k)+1.0D-30)
    If (i<ndm) Goto 26
    Do i = 1, ndm
      xi(i, j) = xin(i)
    End Do
    xi(nd, j) = one
  End Do
!
  If (it<itmx .And. acc*dabs(avgi)<sd) Goto 9
!200   FORMAT('0INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,'  NCALL=',F8.0
!     1    /28X,'  IT=',I5,'  ITMX=',I5/28X,'  ACC=',G9.3
!     2    /28X,'  MDS=',I3,'   ND=',I4/28X,'  (XL,XU)=',
!     3    (T40,'( ',G12.6,' , ',G12.6,' )'))
!201   FORMAT(///' INTEGRATION BY VEGAS' / '0ITERATION NO.',I3,
!     1    ':   INTEGRAL =',G14.8/21X,'STD DEV  =',G10.4 /
!     2    ' ACCUMULATED RESULTS:   INTEGRAL =',G14.8 /
!     3    24X,'STD DEV  =',G10.4 / 24X,'CHI**2 PER IT''N =',G10.4)
!202   FORMAT('0DATA FOR AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
!     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
!     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
!     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
  Return
End Subroutine vegas
!
!
Subroutine aran9(qran, ndim)
  Dimension qran(10)
  Common /sedvax/num1
  Save
  Do i = 1, ndim
    qran(i) = ranart(num1)
  End Do
  Return
End Subroutine aran9

!
!
!*********GAUSSIAN ONE-DIMENSIONAL INTEGRATION PROGRAM*************
!
Function gauss1(f, a, b, eps)
  External f
  Dimension w(12), x(12)
  Save
  Data const/1.0E-12/
  Data w/0.1012285, .2223810, .3137067, .3623838, .0271525, .0622535, 0.0951585, .1246290, .1495960, .1691565, .1826034, .1894506/
  Data x/0.9602899, .7966665, .5255324, .1834346, .9894009, .9445750, 0.8656312, .7554044, .6178762, .4580168, .2816036, .0950125/

  delta = const*abs(a-b)
  gauss1 = 0.0
  aa = a
  5 y = b - aa
  If (abs(y)<=delta) Return
  2 bb = aa + y
  c1 = 0.5*(aa+bb)
  c2 = c1 - aa
  s8 = 0.0
  s16 = 0.0
  Do i = 1, 4
    u = x(i)*c2
    s8 = s8 + w(i)*(f(c1+u)+f(c1-u))
  End Do
  Do i = 5, 12
    u = x(i)*c2
    s16 = s16 + w(i)*(f(c1+u)+f(c1-u))
  End Do
  s8 = s8*c2
  s16 = s16*c2
  If (abs(s16-s8)>eps*(1.+abs(s16))) Goto 4
  gauss1 = gauss1 + s16
  aa = bb
  Goto 5
  4 y = 0.5*y
  If (abs(y)>delta) Goto 2
  Write (6, 7)
  gauss1 = 0.0
  Return
  7 Format (1X, 'GAUSS1....TOO HIGH ACURACY REQUIRED')
End Function gauss1
!
!
!
Function gauss2(f, a, b, eps)
  External f
  Dimension w(12), x(12)
  Save
  Data const/1.0E-12/
  Data w/0.1012285, .2223810, .3137067, .3623838, .0271525, .0622535, 0.0951585, .1246290, .1495960, .1691565, .1826034, .1894506/
  Data x/0.9602899, .7966665, .5255324, .1834346, .9894009, .9445750, 0.8656312, .7554044, .6178762, .4580168, .2816036, .0950125/

  delta = const*abs(a-b)
  gauss2 = 0.0
  aa = a
  5 y = b - aa
  If (abs(y)<=delta) Return
  2 bb = aa + y
  c1 = 0.5*(aa+bb)
  c2 = c1 - aa
  s8 = 0.0
  s16 = 0.0
  Do i = 1, 4
    u = x(i)*c2
    s8 = s8 + w(i)*(f(c1+u)+f(c1-u))
  End Do
  Do i = 5, 12
    u = x(i)*c2
    s16 = s16 + w(i)*(f(c1+u)+f(c1-u))
  End Do
  s8 = s8*c2
  s16 = s16*c2
  If (abs(s16-s8)>eps*(1.+abs(s16))) Goto 4
  gauss2 = gauss2 + s16
  aa = bb
  Goto 5
  4 y = 0.5*y
  If (abs(y)>delta) Goto 2
  Write (6, 7)
  gauss2 = 0.0
  Return
  7 Format (1X, 'GAUSS2....TOO HIGH ACURACY REQUIRED')
End Function gauss2
!
!
!
!
!
Subroutine title

  Common /rndf77/nseed
!c      SAVE /RNDF77/
  Save

  ! Write (6, 200)
  Return
! !lin-8/15/02 f77:
! !200        FORMAT(//10X,
! !     &        '**************************************************'/10X,
! !     &  '*     |      \       _______      /  ------/     *'/10X,
! !     &        '*   ----- ------     |_____|     /_/     /       *'/10X,
! !     &        '*    ||\    /        |_____|      /    / \       *'/10X,
! !     &        '*    /| \  /_/       /_______    /_  /    \_     *'/10X,
! !     &        '*   / |     / /     /  /  / |        -------     *'/10X,
! !     &        '*     |    / /\       /  /  |     /     |        *'/10X,
! !     &        '*     |   / /  \     /  / \_|    /   -------     *'/10X,
!   200 Format (//10X, '**************************************************'/10X, '*     |      |       _______      /  ------/     *'/10X, '*   ----- ------     |_____|     /_/     /       *'/10X, '*    |||    /        |_____|      /    / |       *'/10X, '*    /| |  /_/       /_______    /_  /    |      *'/10X, '*   / |     / /     /  /  / |        -------     *'/10X, '*     |    / /|       /  /  |     /     |        *'/10X, '*     |   / /  |     /  /  _|    /   -------     *'/10X, '*                                                *'/10X, '**************************************************'/10X, '                      HIJING                      '/10X, '       Heavy Ion Jet INteraction Generator        '/10X, '                        by                        '/10X, '            X. N. Wang  and  M. Gyulassy           '/10X, '             Lawrence Berkeley Laboratory           '//)
End Subroutine title
