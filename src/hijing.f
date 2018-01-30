        SUBROUTINE HIJING(FRAME,BMIN0,BMAX0)
cbz1/25/99
        PARAMETER (MAXPTN=400001)
clin-4/20/01        PARAMETER (MAXSTR = 1600)
        PARAMETER (MAXSTR=150001)
cbz1/25/99end
clin-4/26/01:
        PARAMETER (MAXIDL=4001)
cbz1/31/99
        DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
        DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
        DOUBLE PRECISION  ATAUI, ZT1, ZT2, ZT3
        DOUBLE PRECISION  xnprod,etprod,xnfrz,etfrz,
     & dnprod,detpro,dnfrz,detfrz
clin-8/2015:
        DOUBLE PRECISION vxp0,vyp0,vzp0,xstrg0,ystrg0,xstrg,ystrg
cbz1/31/99end
        CHARACTER FRAME*8
        DIMENSION SCIP(300,300),RNIP(300,300),SJIP(300,300),JTP(3),
     &                        IPCOL(90000),ITCOL(90000)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
C
        COMMON/hjcrdn/YP(3,300),YT(3,300)
cc      SAVE /hjcrdn/
clin-7/16/03 NINT is a intrinsic fortran function, rename it to NINTHJ
c        COMMON/HJGLBR/NELT,NINT,NELP,NINP
        COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
cc      SAVE /HJGLBR/
        COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
cc      SAVE /HMAIN1/
clin-4/26/01
c        COMMON/HMAIN2/KATT(130000,4),PATT(130000,4)
        COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
cc      SAVE /HMAIN2/
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
cc      SAVE /HJJET1/
clin-4/2008
c        COMMON/HJJET2/NSG,NJSG(900),IASG(900,3),K1SG(900,100),
c     &       K2SG(900,100),PXSG(900,100),PYSG(900,100),
c     &       PZSG(900,100),PESG(900,100),PMSG(900,100)
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
        COMMON/HJJET4/NDR,IADR(MAXSTR,2),KFDR(MAXSTR),PDR(MAXSTR,5)
clin-4/2008:
c        common/xydr/rtdr(900,2)
        common/xydr/rtdr(MAXSTR,2)
cc      SAVE /HJJET4/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
C
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)
cc      SAVE /LUJETS/
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
cc      SAVE /LUDAT1/
clin-9/29/03 changed name in order to distinguish from /prec2/
        COMMON /ARPRC/ ITYPAR(MAXSTR),
     &       GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &       PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &       XMAR(MAXSTR)
ccbz11/11/98
c        COMMON /ARPRC/ ITYP(MAXSTR),
c     &     GX(MAXSTR), GY(MAXSTR), GZ(MAXSTR), FT(MAXSTR),
c     &     PX(MAXSTR), PY(MAXSTR), PZ(MAXSTR), EE(MAXSTR),
c     &     XM(MAXSTR)
cc      SAVE /ARPRC/
ccbz11/11/98end
cbz1/25/99
        COMMON /PARA1/ MUL
cc      SAVE /PARA1/
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
cc      SAVE /prec1/
        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
cc      SAVE /prec2/
        COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
cc      SAVE /ilist7/
        COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
cc      SAVE /ilist8/
        COMMON /SREC1/ NSP, NST, NSI
cc      SAVE /SREC1/
        COMMON /SREC2/ATAUI(MAXSTR),ZT1(MAXSTR),ZT2(MAXSTR),ZT3(MAXSTR)
cc      SAVE /SREC2/
cbz1/25/99end
clin-2/25/00
        COMMON /frzout/ xnprod(30),etprod(30),xnfrz(30),etfrz(30),
     & dnprod(30),detpro(30),dnfrz(30),detfrz(30)
cc      SAVE /frzout/
clin-4/11/01 soft:
      common/anim/nevent,isoft,isflag,izpc
cc      SAVE /anim/
clin-4/25/01 soft3:
      DOUBLE PRECISION PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
clin-4/26/01 lepton and photon info:
        COMMON /NOPREC/ NNOZPC, ITYPN(MAXIDL),
     &       GXN(MAXIDL), GYN(MAXIDL), GZN(MAXIDL), FTN(MAXIDL),
     &       PXN(MAXIDL), PYN(MAXIDL), PZN(MAXIDL), EEN(MAXIDL),
     &       XMN(MAXIDL)
cc      SAVE /NOPREC/
clin-6/22/01:
        common /lastt/itimeh,bimp
cc      SAVE /lastt/
        COMMON /AREVT/ IAEVT, IARUN, MISS
        common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
clin-7/2011 ioscar value is needed:
        common /para7/ ioscar,nsmbbbar,nsmmeson
clin-2/2012 allow random orientation of reaction plane:
        common /phiHJ/iphirp,phiRP
clin-8/2015:
        common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1       xstrg0(MAXPTN),ystrg0(MAXPTN),
     2       xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
        SAVE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  碰撞参数的最大值和最小值
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        BMAX=MIN(BMAX0,HIPR1(34)+HIPR1(35))
        BMIN=MIN(BMIN0,BMAX)
        IF(IHNT2(1).LE.1 .AND. IHNT2(3).LE.1) THEN
                BMIN=0.0
                BMAX=2.5*SQRT(HIPR1(31)*0.1/HIPR1(40))
        ENDIF
C                        ********HIPR1(31) is in mb =0.1fm**2
C*******THE FOLLOWING IS TO SELECT THE COORDINATIONS OF NUCLEONS
C       BOTH IN PROJECTILE AND TARGET NUCLEAR( in fm)
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  给出弹核的每个核子的三坐标
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        YP(1,1)=0.0
        YP(2,1)=0.0
        YP(3,1)=0.0
        IF(IHNT2(1).LE.1) GO TO 14
        DO 10 KP=1,IHNT2(1)
5        R=HIRND(1)
        X=RANART(NSEED)
        CX=2.0*X-1.0
        SX=SQRT(1.0-CX*CX)
C                ********choose theta from uniform cos(theta) distr
        PHI=RANART(NSEED)*2.0*HIPR1(40)
C                ********choose phi form uniform phi distr 0 to 2*pi
        YP(1,KP)=R*SX*COS(PHI)
        YP(2,KP)=R*SX*SIN(PHI)
        YP(3,KP)=R*CX
        IF(HIPR1(29).EQ.0.0) GO TO 10
        DO 8  KP2=1,KP-1
                DNBP1=(YP(1,KP)-YP(1,KP2))**2
                DNBP2=(YP(2,KP)-YP(2,KP2))**2
                DNBP3=(YP(3,KP)-YP(3,KP2))**2
                DNBP=DNBP1+DNBP2+DNBP3
                IF(DNBP.LT.HIPR1(29)*HIPR1(29)) GO TO 5
C                        ********two neighbors cannot be closer than
C                                HIPR1(29)
8        CONTINUE
10        CONTINUE
clin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f,
c     but modified [divide by 2, & x(p)=-x(n)]:
c     (Note: hijing1.383.f has corrected this bug in hijing1.382.f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  氘核的特殊处理
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(IHNT2(1).EQ.2) then
           rnd1=max(RANART(NSEED),1.0e-20)
           rnd2=max(RANART(NSEED),1.0e-20)
           rnd3=max(RANART(NSEED),1.0e-20)
           R=-(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &          +4.38*0.85*log(rnd3)/(4.38+0.85))
           X=RANART(NSEED)
           CX=2.0*X-1.0
           SX=SQRT(1.0-CX*CX)
           PHI=RANART(NSEED)*2.0*HIPR1(40)
c     R above is the relative distance between p & n in a deuteron:
           R=R/2.
           YP(1,1)=R*SX*COS(PHI)
           YP(2,1)=R*SX*SIN(PHI)
           YP(3,1)=R*CX
c     p & n has opposite coordinates in the deuteron frame:
           YP(1,2)=-YP(1,1)
           YP(2,2)=-YP(2,1)
           YP(3,2)=-YP(3,1)
        endif
        DO 12 I=1,IHNT2(1)-1
        DO 12 J=I+1,IHNT2(1)
        IF(YP(3,I).GT.YP(3,J)) GO TO 12
        Y1=YP(1,I)
        Y2=YP(2,I)
        Y3=YP(3,I)
        YP(1,I)=YP(1,J)
        YP(2,I)=YP(2,J)
        YP(3,I)=YP(3,J)
        YP(1,J)=Y1
        YP(2,J)=Y2
        YP(3,J)=Y3
12        CONTINUE
C
C******************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc 初始化靶核每个核子的三坐标
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
14        YT(1,1)=0.0
        YT(2,1)=0.0
        YT(3,1)=0.0
        IF(IHNT2(3).LE.1) GO TO 24
        DO 20 KT=1,IHNT2(3)
15        R=HIRND(2)
        X=RANART(NSEED)
        CX=2.0*X-1.0
        SX=SQRT(1.0-CX*CX)
C                ********choose theta from uniform cos(theta) distr
        PHI=RANART(NSEED)*2.0*HIPR1(40)
C                ********chose phi form uniform phi distr 0 to 2*pi
        YT(1,KT)=R*SX*COS(PHI)
        YT(2,KT)=R*SX*SIN(PHI)
        YT(3,KT)=R*CX
        IF(HIPR1(29).EQ.0.0) GO TO 20
        DO 18  KT2=1,KT-1
                DNBT1=(YT(1,KT)-YT(1,KT2))**2
                DNBT2=(YT(2,KT)-YT(2,KT2))**2
                DNBT3=(YT(3,KT)-YT(3,KT2))**2
                DNBT=DNBT1+DNBT2+DNBT3
                IF(DNBT.LT.HIPR1(29)*HIPR1(29)) GO TO 15
C                        ********two neighbors cannot be closer than
C                                HIPR1(29)
18        CONTINUE
20        CONTINUE
c
clin-1/27/03 Hulthen wavefn for deuteron borrowed from hijing1.382.f,
c     but modified [divide by 2, & x(p)=-x(n)]:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  氘核的特殊处理
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(IHNT2(3).EQ.2) then
           rnd1=max(RANART(NSEED),1.0e-20)
           rnd2=max(RANART(NSEED),1.0e-20)
           rnd3=max(RANART(NSEED),1.0e-20)
           R=-(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &          +4.38*0.85*log(rnd3)/(4.38+0.85))
           X=RANART(NSEED)
           CX=2.0*X-1.0
           SX=SQRT(1.0-CX*CX)
           PHI=RANART(NSEED)*2.0*HIPR1(40)
           R=R/2.
           YT(1,1)=R*SX*COS(PHI)
           YT(2,1)=R*SX*SIN(PHI)
           YT(3,1)=R*CX
           YT(1,2)=-YT(1,1)
           YT(2,2)=-YT(2,1)
           YT(3,2)=-YT(3,1)
        endif
c
        DO 22 I=1,IHNT2(3)-1
        DO 22 J=I+1,IHNT2(3)
        IF(YT(3,I).LT.YT(3,J)) GO TO 22
        Y1=YT(1,I)
        Y2=YT(2,I)
        Y3=YT(3,I)
        YT(1,I)=YT(1,J)
        YT(2,I)=YT(2,J)
        YT(3,I)=YT(3,J)
        YT(1,J)=Y1
        YT(2,J)=Y2
        YT(3,J)=Y3
22        CONTINUE
C********************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  HIJING初始化失败的次数。(input.ampt中设定的上限是1000)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
24        MISS=-1
50        MISS=MISS+1
clin-6/2009
c        IF(MISS.GT.50) THEN
        IF(MISS.GT.maxmiss) THEN
           WRITE(6,*) 'infinite loop happened in  HIJING'
           STOP
        ENDIF
clin-4/30/01:
        itest=0
        NATT=0
        JATT=0
        EATT=0.0
        CALL HIJINI
        NLOP=0
C                        ********Initialize for a new event
60        NT=0
        NP=0
        N0=0
        N01=0
        N10=0
        N11=0
        NELT=0
        NINTHJ=0
        NELP=0
        NINP=0
        NSG=0
        NCOLT=0
C****        BB IS THE ABSOLUTE VALUE OF IMPACT PARAMETER,BB**2 IS
C       RANDOMLY GENERATED AND ITS ORIENTATION IS RANDOMLY SET
C       BY THE ANGLE PHI  FOR EACH COLLISION.******************
C
        BB=SQRT(BMIN**2+RANART(NSEED)*(BMAX**2-BMIN**2))
cbz6/28/99 flow1
clin-2/2012:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  反应平面的控制。
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        PHI=0.
        if(iphirp.eq.1) PHI=2.0*HIPR1(40)*RANART(NSEED)
        phiRP=phi
cbz6/28/99 flow1 end
        BBX=BB*COS(PHI)
        BBY=BB*SIN(PHI)
        HINT1(19)=BB
        HINT1(20)=PHI
C
        DO 70 JP=1,IHNT2(1)
        DO 70 JT=1,IHNT2(3)
           SCIP(JP,JT)=-1.0
           B2=(YP(1,JP)+BBX-YT(1,JT))**2+(YP(2,JP)+BBY-YT(2,JT))**2
           R2=B2*HIPR1(40)/HIPR1(31)/0.1
C                ********mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
           RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)
     &          /1.2**2/REAL(IHNT2(1))**0.6666667,1.0)
           RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)
     &          /1.2**2/REAL(IHNT2(3))**0.6666667,1.0)
           APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
           APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
           HINT1(18)=HINT1(14)-APHX1*HINT1(15)
     &                        -APHX2*HINT1(16)+APHX1*APHX2*HINT1(17)
           IF(IHPR2(14).EQ.0.OR.
     &          (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) THEN
              GS=1.0-EXP(-(HIPR1(30)+HINT1(18))*ROMG(R2)/HIPR1(31))
              RANTOT=RANART(NSEED)
              IF(RANTOT.GT.GS) GO TO 70
              GO TO 65
           ENDIF
           GSTOT0=2.0*(1.0-EXP(-(HIPR1(30)+HINT1(18))
     &             /HIPR1(31)/2.0*ROMG(0.0)))
           R2=R2/GSTOT0
           GS=1.0-EXP(-(HIPR1(30)+HINT1(18))/HIPR1(31)*ROMG(R2))
           GSTOT=2.0*(1.0-SQRT(1.0-GS))
           RANTOT=RANART(NSEED)*GSTOT0
           IF(RANTOT.GT.GSTOT) GO TO 70
           IF(RANTOT.GT.GS) THEN
              CALL HIJCSC(JP,JT)
              GO TO 70
C                        ********perform elastic collisions
           ENDIF
 65           SCIP(JP,JT)=R2
           RNIP(JP,JT)=RANTOT
           SJIP(JP,JT)=HINT1(18)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  到目前遭遇的碰撞的次數
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           NCOLT=NCOLT+1
           IPCOL(NCOLT)=JP
           ITCOL(NCOLT)=JT
70        CONTINUE
C                ********total number interactions proj and targ has
C                                suffered
clin-5/22/01 write impact parameter:
        bimp=bb
        write(6,*) '#impact parameter,nlop,ncolt=',bimp,nlop,ncolt
        IF(NCOLT.EQ.0) THEN
           NLOP=NLOP+1
           IF(NLOP.LE.20.OR.
     &           (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) GO TO 60
           RETURN
        ENDIF
C               ********At large impact parameter, there maybe no
C                       interaction at all. For NN collision
C                       repeat the event until interaction happens
C
        IF(IHPR2(3).NE.0) THEN
           NHARD=1+INT(RANART(NSEED)*(NCOLT-1)+0.5)
           NHARD=MIN(NHARD,NCOLT)
           JPHARD=IPCOL(NHARD)
           JTHARD=ITCOL(NHARD)
clin-6/2009 ctest off:
c           write(99,*) IAEVT,NHARD,NCOLT,JPHARD,JTHARD
        ENDIF
C
        IF(IHPR2(9).EQ.1) THEN
                NMINI=1+INT(RANART(NSEED)*(NCOLT-1)+0.5)
                NMINI=MIN(NMINI,NCOLT)
                JPMINI=IPCOL(NMINI)
                JTMINI=ITCOL(NMINI)
        ENDIF
C                ********Specifying the location of the hard and
C                        minijet if they are enforced by user
C
        DO 200 JP=1,IHNT2(1)
        DO 200 JT=1,IHNT2(3)
        IF(SCIP(JP,JT).EQ.-1.0) GO TO 200
                NFP(JP,11)=NFP(JP,11)+1
                NFT(JT,11)=NFT(JT,11)+1
        IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).GT.1) THEN
                NP=NP+1
                N01=N01+1
        ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).LE.1) THEN
                NT=NT+1
                N10=N10+1
        ELSE IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).LE.1) THEN
                NP=NP+1
                NT=NT+1
                N0=N0+1
        ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).GT.1) THEN
                N11=N11+1
        ENDIF
        JOUT=0
        NFP(JP,10)=0
        NFT(JT,10)=0
C*****************************************************************
        IF(IHPR2(8).EQ.0 .AND. IHPR2(3).EQ.0) GO TO 160
C                ********When IHPR2(8)=0 no jets are produced
        IF(NFP(JP,6).LT.0 .OR. NFT(JT,6).LT.0) GO TO 160
C                ********jets can not be produced for (JP,JT)
C                        because not enough energy avaible for
C                                JP or JT
        R2=SCIP(JP,JT)
        HINT1(18)=SJIP(JP,JT)
        TT=ROMG(R2)*HINT1(18)/HIPR1(31)
        TTS=HIPR1(30)*ROMG(R2)/HIPR1(31)
        NJET=0
        IF(IHPR2(3).NE.0 .AND. JP.EQ.JPHARD .AND. JT.EQ.JTHARD) THEN
           CALL JETINI(JP,JT,1)
           CALL HIJHRD(JP,JT,0,JFLG,0)
           HINT1(26)=HINT1(47)
           HINT1(27)=HINT1(48)
           HINT1(28)=HINT1(49)
           HINT1(29)=HINT1(50)
           HINT1(36)=HINT1(67)
           HINT1(37)=HINT1(68)
           HINT1(38)=HINT1(69)
           HINT1(39)=HINT1(70)
C
           IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
           IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
           IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &                                JFLG.GE.3) IASG(NSG,3)=1
           IHNT2(9)=IHNT2(14)
           IHNT2(10)=IHNT2(15)
           DO 105 I05=1,5
              HINT1(20+I05)=HINT1(40+I05)
              HINT1(30+I05)=HINT1(50+I05)
 105           CONTINUE
clin-6/2009 ctest off:
c           write(99,*) jp,jt,IHPR2(3),HIPR1(10),njet,
c     1          ihnt2(9),hint1(21),hint1(22),hint1(23),
c     2          ihnt2(10),hint1(31),hint1(32),hint1(33)
c           write(99,*) ' '
           JOUT=1
           IF(IHPR2(8).EQ.0) GO TO 160
           RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)/1.2**2
     &                /REAL(IHNT2(1))**0.6666667,1.0)
           RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)/1.2**2
     &                /REAL(IHNT2(3))**0.6666667,1.0)
           APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
           APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
           HINT1(65)=HINT1(61)-APHX1*HINT1(62)
     &                        -APHX2*HINT1(63)+APHX1*APHX2*HINT1(64)
           TTRIG=ROMG(R2)*HINT1(65)/HIPR1(31)
           NJET=-1
C                ********subtract the trigger jet from total number
C                        of jet production  to be done since it has
C                                already been produced here
           XR1=-ALOG(EXP(-TTRIG)+RANART(NSEED)*(1.0-EXP(-TTRIG)))
 106           NJET=NJET+1
           XR1=XR1-ALOG(max(RANART(NSEED),1.0e-20))
           IF(XR1.LT.TTRIG) GO TO 106
           XR=0.0
 107           NJET=NJET+1
           XR=XR-ALOG(max(RANART(NSEED),1.0e-20))
           IF(XR.LT.TT-TTRIG) GO TO 107
           NJET=NJET-1
           GO TO 112
        ENDIF
C                ********create a hard interaction with specified P_T
c                                 when IHPR2(3)>0
        IF(IHPR2(9).EQ.1.AND.JP.EQ.JPMINI.AND.JT.EQ.JTMINI) GO TO 110
C                ********create at least one pair of mini jets
C                        when IHPR2(9)=1
C
clin-4/15/2010 changed .LT. to .LE. to avoid problem when two sides are equal;
c     this problem may lead to a jet production when there should be none and
c     crash the run; crashes at low energies were reported by P. Bhaduri.
c        IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LT.EXP(-TT)*
c     &                (1.0-EXP(-TTS))) GO TO 160
        IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LE.EXP(-TT)*
     &                 (1.0-EXP(-TTS))) GO TO 160
c
C                ********this is the probability for no jet production
110        XR=-ALOG(EXP(-TT)+RANART(NSEED)*(1.0-EXP(-TT)))
111        NJET=NJET+1
        XR=XR-ALOG(max(RANART(NSEED),1.0e-20))
        IF(XR.LT.TT) GO TO 111
112        NJET=MIN(NJET,IHPR2(8))
        IF(IHPR2(8).LT.0)  NJET=ABS(IHPR2(8))
C                ******** Determine number of mini jet production
C
        DO 150 ijet=1,NJET
           CALL JETINI(JP,JT,0)
           CALL HIJHRD(JP,JT,JOUT,JFLG,1)
C                ********JFLG=1 jets valence quarks, JFLG=2 with
C                        gluon jet, JFLG=3 with q-qbar prod for
C                        (JP,JT). If JFLG=0 jets can not be produced
C                        this time. If JFLG=-1, error occured abandon
C                        this event. JOUT is the total hard scat for
C                        (JP,JT) up to now.
           IF(JFLG.EQ.0) GO TO 160
           IF(JFLG.LT.0) THEN
              IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJHRD'
              GO TO 50
           ENDIF
           JOUT=JOUT+1
           IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
           IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
           IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &                        JFLG.GE.3) IASG(NSG,3)=1
C                ******** jet with PT>HIPR1(11) will be quenched
 150        CONTINUE
 160        CONTINUE
        CALL HIJSFT(JP,JT,JOUT,IERROR)
        IF(IERROR.NE.0) THEN
           IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJSFT'
           GO TO 50
        ENDIF
C
C                ********conduct soft scattering between JP and JT
        JATT=JATT+JOUT
200        CONTINUE
c
c**************************
c
clin-6/2009 write out initial minijet information:
clin-2/2012:
c           call minijet_out(BB)
           call minijet_out(BB,phiRP)
           if(pttrig.gt.0.and.ntrig.eq.0) goto 50
clin-4/2012
clin-6/2009 write out initial transverse positions of initial nucleons:
c           write(94,*) IAEVT,MISS,IHNT2(1),IHNT2(3)
        DO 201 JP=1,IHNT2(1)
clin-6/2009:
c           write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5)
clin-2/2012:
c       write(94,203) YP(1,JP)+0.5*BB, YP(2,JP), JP, NFP(JP,5),yp(3,jp)
clin-4/2012:
c           write(94,203) YP(1,JP)+0.5*BB*cos(phiRP),
c     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
           IF(NFP(JP,5).GT.2) THEN
              NINP=NINP+1
           ELSE IF(NFP(JP,5).EQ.2.OR.NFP(JP,5).EQ.1) THEN
              NELP=NELP+1
           ENDIF
 201    continue
        DO 202 JT=1,IHNT2(3)
clin-6/2009 target nucleon # has a minus sign for distinction from projectile:
c           write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5)
clin-2/2012:
c       write(94,203) YT(1,JT)-0.5*BB, YT(2,JT), -JT, NFT(JT,5),yt(3,jt)
clin-4/2012:
c           write(94,203) YT(1,JT)-0.5*BB*cos(phiRP),
c     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
           IF(NFT(JT,5).GT.2) THEN
              NINTHJ=NINTHJ+1
           ELSE IF(NFT(JT,5).EQ.2.OR.NFT(JT,5).EQ.1) THEN
              NELT=NELT+1
           ENDIF
 202    continue
c 203    format(f10.3,1x,f10.3,2(1x,I5))
c 203    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
c
c*******************************
C********perform jet quenching for jets with PT>HIPR1(11)**********
        IF((IHPR2(8).NE.0.OR.IHPR2(3).NE.0).AND.IHPR2(4).GT.0.AND.
     &                        IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
                DO 271 I=1,IHNT2(1)
                        IF(NFP(I,7).EQ.1) CALL QUENCH(I,1)
271                CONTINUE
                DO 272 I=1,IHNT2(3)
                        IF(NFT(I,7).EQ.1) CALL QUENCH(I,2)
272                CONTINUE
                DO 273 ISG=1,NSG
                        IF(IASG(ISG,3).EQ.1) CALL QUENCH(ISG,3)
273                CONTINUE
        ENDIF
clin*****4/09/01-soft1, default way of treating strings:
        if(isoft.eq.1) then
clin-4/16/01 allow fragmentation:
           isflag=1
cbz1/25/99
c.....transfer data from HIJING to ZPC
        NSP = IHNT2(1)
        NST = IHNT2(3)
        NSI = NSG
        ISTR = 0
        NPAR = 0
        DO 1008 I = 1, IHNT2(1)
           ISTR = ISTR + 1
           DO 1007 J = 1, NPJ(I)
cbz1/27/99
c.....for now only consider gluon cascade
              IF (KFPJ(I, J) .EQ. 21) THEN
cbz1/27/99end
              NPAR = NPAR + 1
              LSTRG0(NPAR) = ISTR
              LPART0(NPAR) = J
              ITYP0(NPAR) = KFPJ(I, J)
cbz6/28/99 flow1
clin-7/20/01 add dble or sngl to make precisions consistent
c              GX0(NPAR) = YP(1, I)
clin-2/2012:
c              GX0(NPAR) = dble(YP(1, I) + 0.5 * BB)
              GX0(NPAR) = dble(YP(1, I)+0.5*BB*cos(phiRP))
cbz6/28/99 flow1 end
c              GY0(NPAR) = dble(YP(2, I))
              GY0(NPAR) = dble(YP(2, I)+0.5*BB*sin(phiRP))
              GZ0(NPAR) = 0d0
              FT0(NPAR) = 0d0
              PX0(NPAR) = dble(PJPX(I, J))
              PY0(NPAR) = dble(PJPY(I, J))
              PZ0(NPAR) = dble(PJPZ(I, J))
              XMASS0(NPAR) = dble(PJPM(I, J))
c              E0(NPAR) = dble(PJPE(I, J))
              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
clin-7/20/01-end
cbz1/27/99
c.....end gluon selection
              END IF
cbz1/27/99end
 1007      CONTINUE
 1008   CONTINUE
        DO 1010 I = 1, IHNT2(3)
           ISTR = ISTR + 1
           DO 1009 J = 1, NTJ(I)
cbz1/27/99
c.....for now only consider gluon cascade
              IF (KFTJ(I, J) .EQ. 21) THEN
cbz1/27/99end
              NPAR = NPAR + 1
              LSTRG0(NPAR) = ISTR
              LPART0(NPAR) = J
              ITYP0(NPAR) = KFTJ(I, J)
cbz6/28/99 flow1
clin-7/20/01 add dble or sngl to make precisions consistent
c              GX0(NPAR) = YT(1, I)
clin-2/2012:
c              GX0(NPAR) = dble(YT(1, I) - 0.5 * BB)
              GX0(NPAR) = dble(YT(1, I)-0.5*BB*cos(phiRP))
cbz6/28/99 flow1 end
c              GY0(NPAR) = dble(YT(2, I))
              GY0(NPAR) = dble(YT(2, I)-0.5*BB*sin(phiRP))
              GZ0(NPAR) = 0d0
              FT0(NPAR) = 0d0
              PX0(NPAR) = dble(PJTX(I, J))
              PY0(NPAR) = dble(PJTY(I, J))
              PZ0(NPAR) = dble(PJTZ(I, J))
              XMASS0(NPAR) = dble(PJTM(I, J))
c              E0(NPAR) = dble(PJTE(I, J))
              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
cbz1/27/99
c.....end gluon selection
              END IF
cbz1/27/99end
 1009      CONTINUE
 1010   CONTINUE
        DO 1012 I = 1, NSG
           ISTR = ISTR + 1
           DO 1011 J = 1, NJSG(I)
cbz1/27/99
c.....for now only consider gluon cascade
              IF (K2SG(I, J) .EQ. 21) THEN
cbz1/27/99end
              NPAR = NPAR + 1
              LSTRG0(NPAR) = ISTR
              LPART0(NPAR) = J
              ITYP0(NPAR) = K2SG(I, J)
clin-7/20/01 add dble or sngl to make precisions consistent:
              GX0(NPAR) = 0.5d0 *
     1             dble(YP(1, IASG(I, 1)) + YT(1, IASG(I, 2)))
              GY0(NPAR) = 0.5d0 *
     2             dble(YP(2, IASG(I, 1)) + YT(2, IASG(I, 2)))
              GZ0(NPAR) = 0d0
              FT0(NPAR) = 0d0
              PX0(NPAR) = dble(PXSG(I, J))
              PY0(NPAR) = dble(PYSG(I, J))
              PZ0(NPAR) = dble(PZSG(I, J))
              XMASS0(NPAR) = dble(PMSG(I, J))
c              E0(NPAR) = dble(PESG(I, J))
              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
cbz1/27/99
c.....end gluon selection
              END IF
cbz1/27/99end
 1011      CONTINUE
 1012   CONTINUE
        MUL = NPAR
cbz2/4/99
        CALL HJANA1
cbz2/4/99end
clin-6/2009:
        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
c.....call ZPC for parton cascade
        CALL ZPCMN
c     write out parton and wounded nucleon information to ana/zpc1.mom:
clin-6/2009:
c        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        DO 1013 I = 1, MUL
cc           WRITE (14, 411) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c     &        XMASS5(I), E5(I)
           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
     1          .lt.9999) then
              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           else
c     change format for large numbers:
              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           endif
 1013   CONTINUE
 210    format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 211    format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
 395    format(3I8,f10.4,4I5)
clin-4/09/01:
        itest=itest+1
c 411    FORMAT(1X, 3F10.3, I6, 2F10.3)
cbz3/19/99 end
clin-5/2009 ctest off:
c        call frztm(1,1)
c.....transfer data back from ZPC to HIJING
        DO 1014 I = 1, MUL
           IF (LSTRG1(I) .LE. NSP) THEN
              NSTRG = LSTRG1(I)
              NPART = LPART1(I)
              KFPJ(NSTRG, NPART) = ITYP5(I)
clin-7/20/01 add dble or sngl to make precisions consistent
              PJPX(NSTRG, NPART) = sngl(PX5(I))
              PJPY(NSTRG, NPART) = sngl(PY5(I))
              PJPZ(NSTRG, NPART) = sngl(PZ5(I))
              PJPE(NSTRG, NPART) = sngl(E5(I))
              PJPM(NSTRG, NPART) = sngl(XMASS5(I))
           ELSE IF (LSTRG1(I) .LE. NSP + NST) THEN
              NSTRG = LSTRG1(I) - NSP
              NPART = LPART1(I)
              KFTJ(NSTRG, NPART) = ITYP5(I)
              PJTX(NSTRG, NPART) = sngl(PX5(I))
              PJTY(NSTRG, NPART) = sngl(PY5(I))
              PJTZ(NSTRG, NPART) = sngl(PZ5(I))
              PJTE(NSTRG, NPART) = sngl(E5(I))
              PJTM(NSTRG, NPART) = sngl(XMASS5(I))
           ELSE
              NSTRG = LSTRG1(I) - NSP - NST
              NPART = LPART1(I)
              K2SG(NSTRG, NPART) = ITYP5(I)
              PXSG(NSTRG, NPART) = sngl(PX5(I))
              PYSG(NSTRG, NPART) = sngl(PY5(I))
              PZSG(NSTRG, NPART) = sngl(PZ5(I))
              PESG(NSTRG, NPART) = sngl(E5(I))
              PMSG(NSTRG, NPART) = sngl(XMASS5(I))
           END IF
 1014   CONTINUE
cbz1/25/99end
cbz2/4/99
        CALL HJANA2
cbz2/4/99end
clin*****4/09/01-soft2, put q+dq+X in strings into ZPC:
        elseif(isoft.eq.2) then
        NSP = IHNT2(1)
        NST = IHNT2(3)
clin-4/27/01:
        NSI = NSG
        NPAR=0
        ISTR=0
C
clin  No fragmentation to hadrons, only on parton level,
c     and transfer minijet and string data from HIJING to ZPC:
        MSTJ(1)=0
clin-4/12/01 forbid soft radiation before ZPC to avoid small-mass strings,
c     and forbid jet order reversal before ZPC to avoid unphysical flavors:
        IHPR2(1)=0
        isflag=0
        IF(IHPR2(20).NE.0) THEN
           DO 320 NTP=1,2
              DO 310 jjtp=1,IHNT2(2*NTP-1)
                 ISTR = ISTR + 1
c change: do gluon kink only once: either here or in fragmentation.
                 CALL HIJFRG(jjtp,NTP,IERROR)
c                 call lulist(1)
                 if(NTP.eq.1) then
c 354                continue
                    NPJ(jjtp)=MAX0(N-2,0)
clin-4/12/01:                    NPJ(jjtp)=MAX0(ipartn-2,0)
                 else
c 355                continue
                    NTJ(jjtp)=MAX0(N-2,0)
clin-4/12/01:                    NTJ(jjtp)=MAX0(ipartn-2,0)
                 endif
                 do 300 ii=1,N
                 NPAR = NPAR + 1
                 LSTRG0(NPAR) = ISTR
                 LPART0(NPAR) = II
                 ITYP0(NPAR) = K(II,2)
                 GZ0(NPAR) = 0d0
                 FT0(NPAR) = 0d0
clin-7/20/01 add dble or sngl to make precisions consistent
                 PX0(NPAR) = dble(P(II,1))
                 PY0(NPAR) = dble(P(II,2))
                 PZ0(NPAR) = dble(P(II,3))
                 XMASS0(NPAR) = dble(P(II,5))
c                 E0(NPAR) = dble(P(II,4))
                 E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1                +PZ0(NPAR)**2+XMASS0(NPAR)**2)
                 IF (NTP .EQ. 1) THEN
clin-7/20/01 add dble or sngl to make precisions consistent
clin-2/2012:
c                    GX0(NPAR) = dble(YP(1, jjtp)+0.5 * BB)
c                    GY0(NPAR) = dble(YP(2, jjtp))
                    GX0(NPAR) = dble(YP(1, jjtp)+0.5*BB*cos(phiRP))
                    GY0(NPAR) = dble(YP(2, jjtp)+0.5*BB*sin(phiRP))
                    IITYP=ITYP0(NPAR)
                    nstrg=LSTRG0(NPAR)
                    if(IITYP.eq.2112.or.IITYP.eq.2212) then
                    elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (II.eq.1.or.II.eq.N)) then
                       PP(nstrg,6)=sngl(PX0(NPAR))
                       PP(nstrg,7)=sngl(PY0(NPAR))
                       PP(nstrg,14)=sngl(XMASS0(NPAR))
                    elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(II.eq.1.or.II.eq.N)) then
                       PP(nstrg,8)=sngl(PX0(NPAR))
                       PP(nstrg,9)=sngl(PY0(NPAR))
                       PP(nstrg,15)=sngl(XMASS0(NPAR))
                    else
                       NPART = LPART0(NPAR)-1
                       KFPJ(NSTRG, NPART) = ITYP0(NPAR)
                       PJPX(NSTRG, NPART) = sngl(PX0(NPAR))
                       PJPY(NSTRG, NPART) = sngl(PY0(NPAR))
                       PJPZ(NSTRG, NPART) = sngl(PZ0(NPAR))
                       PJPE(NSTRG, NPART) = sngl(E0(NPAR))
                       PJPM(NSTRG, NPART) = sngl(XMASS0(NPAR))
                    endif
                 ELSE
clin-2/2012:
c                    GX0(NPAR) = dble(YT(1, jjtp)-0.5 * BB)
c                    GY0(NPAR) = dble(YT(2, jjtp))
                    GX0(NPAR) = dble(YT(1, jjtp)-0.5*BB*cos(phiRP))
                    GY0(NPAR) = dble(YT(2, jjtp)-0.5*BB*sin(phiRP))
                    IITYP=ITYP0(NPAR)
                    nstrg=LSTRG0(NPAR)-NSP
                    if(IITYP.eq.2112.or.IITYP.eq.2212) then
                    elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (II.eq.1.or.II.eq.N)) then
                       PT(nstrg,6)=sngl(PX0(NPAR))
                       PT(nstrg,7)=sngl(PY0(NPAR))
                       PT(nstrg,14)=sngl(XMASS0(NPAR))
                    elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(II.eq.1.or.II.eq.N)) then
                       PT(nstrg,8)=sngl(PX0(NPAR))
                       PT(nstrg,9)=sngl(PY0(NPAR))
                       PT(nstrg,15)=sngl(XMASS0(NPAR))
                    else
                       NPART = LPART0(NPAR)-1
                       KFTJ(NSTRG, NPART) = ITYP0(NPAR)
                       PJTX(NSTRG, NPART) = sngl(PX0(NPAR))
                       PJTY(NSTRG, NPART) = sngl(PY0(NPAR))
                       PJTZ(NSTRG, NPART) = sngl(PZ0(NPAR))
                       PJTE(NSTRG, NPART) = sngl(E0(NPAR))
                       PJTM(NSTRG, NPART) = sngl(XMASS0(NPAR))
                    endif
                 END IF
 300          continue
 310          continue
 320       continue
           DO 330 ISG=1,NSG
              ISTR = ISTR + 1
              CALL HIJFRG(ISG,3,IERROR)
c              call lulist(2)
c
              NJSG(ISG)=N
c
              do 1001 ii=1,N
                 NPAR = NPAR + 1
                 LSTRG0(NPAR) = ISTR
                 LPART0(NPAR) = II
                 ITYP0(NPAR) = K(II,2)
                 GX0(NPAR)=0.5d0*
     1                dble(YP(1,IASG(ISG,1))+YT(1,IASG(ISG,2)))
                 GY0(NPAR)=0.5d0*
     2                dble(YP(2,IASG(ISG,1))+YT(2,IASG(ISG,2)))
                 GZ0(NPAR) = 0d0
                 FT0(NPAR) = 0d0
                 PX0(NPAR) = dble(P(II,1))
                 PY0(NPAR) = dble(P(II,2))
                 PZ0(NPAR) = dble(P(II,3))
                 XMASS0(NPAR) = dble(P(II,5))
c                 E0(NPAR) = dble(P(II,4))
                 E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1                +PZ0(NPAR)**2+XMASS0(NPAR)**2)
 1001         continue
 330       continue
        endif
        MUL = NPAR
cbz2/4/99
        CALL HJANA1
cbz2/4/99end
clin-6/2009:
        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
c.....call ZPC for parton cascade
        CALL ZPCMN
cbz3/19/99
clin-6/2009:
c        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        itest=itest+1
        DO 1015 I = 1, MUL
c           WRITE (14, 311) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c     &        XMASS5(I), E5(I)
clin-4/2012 write parton freeze-out position in zpc.dat for this test scenario:
c           WRITE (14, 312) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c     &        XMASS5(I), E5(I),LSTRG1(I), LPART1(I)
           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
     1          .lt.9999) then
              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           else
              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           endif
c
 1015   CONTINUE
c 311    FORMAT(1X, 3F10.4, I6, 2F10.4)
c 312    FORMAT(1X, 3F10.3, I6, 2F10.3,1X,I6,1X,I3)
cbz3/19/99 end
clin-5/2009 ctest off:
c        call frztm(1,1)
clin-4/13/01 initialize four momenta and invariant mass of strings after ZPC:
        do 1004 nmom=1,5
           do 1002 nstrg=1,nsp
              PP(nstrg,nmom)=0.
 1002      continue
           do 1003 nstrg=1,nst
              PT(nstrg,nmom)=0.
 1003      continue
 1004   continue
clin-4/13/01-end
        DO 1005 I = 1, MUL
           IITYP=ITYP5(I)
           IF (LSTRG1(I) .LE. NSP) THEN
              NSTRG = LSTRG1(I)
c     nucleons without interactions:
              if(IITYP.eq.2112.or.IITYP.eq.2212) then
clin-7/20/01 add dble or sngl to make precisions consistent
                 PP(nstrg,1)=sngl(PX5(I))
                 PP(nstrg,2)=sngl(PY5(I))
                 PP(nstrg,3)=sngl(PZ5(I))
                 PP(nstrg,4)=sngl(E5(I))
                 PP(nstrg,5)=sngl(XMASS5(I))
c     valence quark:
              elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (LPART1(I).eq.1.or.LPART1(I).eq.(NPJ(NSTRG)+2))) then
                 PP(nstrg,6)=sngl(PX5(I))
                 PP(nstrg,7)=sngl(PY5(I))
                 PP(nstrg,14)=sngl(XMASS5(I))
                 PP(nstrg,1)=PP(nstrg,1)+sngl(PX5(I))
                 PP(nstrg,2)=PP(nstrg,2)+sngl(PY5(I))
                 PP(nstrg,3)=PP(nstrg,3)+sngl(PZ5(I))
                 PP(nstrg,4)=PP(nstrg,4)+sngl(E5(I))
                 PP(nstrg,5)=sqrt(PP(nstrg,4)**2-PP(nstrg,1)**2
     1                -PP(nstrg,2)**2-PP(nstrg,3)**2)
c     diquark:
              elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(LPART1(I).eq.1.or.LPART1(I).eq.(NPJ(NSTRG)+2))) then
                 PP(nstrg,8)=sngl(PX5(I))
                 PP(nstrg,9)=sngl(PY5(I))
                 PP(nstrg,15)=sngl(XMASS5(I))
                 PP(nstrg,1)=PP(nstrg,1)+sngl(PX5(I))
                 PP(nstrg,2)=PP(nstrg,2)+sngl(PY5(I))
                 PP(nstrg,3)=PP(nstrg,3)+sngl(PZ5(I))
                 PP(nstrg,4)=PP(nstrg,4)+sngl(E5(I))
                 PP(nstrg,5)=sqrt(PP(nstrg,4)**2-PP(nstrg,1)**2
     1                -PP(nstrg,2)**2-PP(nstrg,3)**2)
c     partons in projectile or target strings:
              else
                 NPART = LPART1(I)-1
                 KFPJ(NSTRG, NPART) = ITYP5(I)
                 PJPX(NSTRG, NPART) = sngl(PX5(I))
                 PJPY(NSTRG, NPART) = sngl(PY5(I))
                 PJPZ(NSTRG, NPART) = sngl(PZ5(I))
                 PJPE(NSTRG, NPART) = sngl(E5(I))
                 PJPM(NSTRG, NPART) = sngl(XMASS5(I))
              endif
           ELSE IF (LSTRG1(I) .LE. NSP + NST) THEN
              NSTRG = LSTRG1(I) - NSP
              if(IITYP.eq.2112.or.IITYP.eq.2212) then
                 PT(nstrg,1)=sngl(PX5(I))
                 PT(nstrg,2)=sngl(PY5(I))
                 PT(nstrg,3)=sngl(PZ5(I))
                 PT(nstrg,4)=sngl(E5(I))
                 PT(nstrg,5)=sngl(XMASS5(I))
              elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (LPART1(I).eq.1.or.LPART1(I).eq.(NTJ(NSTRG)+2))) then
                 PT(nstrg,6)=sngl(PX5(I))
                 PT(nstrg,7)=sngl(PY5(I))
                 PT(nstrg,14)=sngl(XMASS5(I))
                 PT(nstrg,1)=PT(nstrg,1)+sngl(PX5(I))
                 PT(nstrg,2)=PT(nstrg,2)+sngl(PY5(I))
                 PT(nstrg,3)=PT(nstrg,3)+sngl(PZ5(I))
                 PT(nstrg,4)=PT(nstrg,4)+sngl(E5(I))
                 PT(nstrg,5)=sqrt(PT(nstrg,4)**2-PT(nstrg,1)**2
     1                -PT(nstrg,2)**2-PT(nstrg,3)**2)
              elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(LPART1(I).eq.1.or.LPART1(I).eq.(NTJ(NSTRG)+2))) then
                 PT(nstrg,8)=sngl(PX5(I))
                 PT(nstrg,9)=sngl(PY5(I))
                 PT(nstrg,15)=sngl(XMASS5(I))
                 PT(nstrg,1)=PT(nstrg,1)+sngl(PX5(I))
                 PT(nstrg,2)=PT(nstrg,2)+sngl(PY5(I))
                 PT(nstrg,3)=PT(nstrg,3)+sngl(PZ5(I))
                 PT(nstrg,4)=PT(nstrg,4)+sngl(E5(I))
                 PT(nstrg,5)=sqrt(PT(nstrg,4)**2-PT(nstrg,1)**2
     1                -PT(nstrg,2)**2-PT(nstrg,3)**2)
              else
                 NPART = LPART1(I)-1
                 KFTJ(NSTRG, NPART) = ITYP5(I)
                 PJTX(NSTRG, NPART) = sngl(PX5(I))
                 PJTY(NSTRG, NPART) = sngl(PY5(I))
                 PJTZ(NSTRG, NPART) = sngl(PZ5(I))
                 PJTE(NSTRG, NPART) = sngl(E5(I))
                 PJTM(NSTRG, NPART) = sngl(XMASS5(I))
              endif
           ELSE
              NSTRG = LSTRG1(I) - NSP - NST
              NPART = LPART1(I)
              K2SG(NSTRG, NPART) = ITYP5(I)
              PXSG(NSTRG, NPART) = sngl(PX5(I))
              PYSG(NSTRG, NPART) = sngl(PY5(I))
              PZSG(NSTRG, NPART) = sngl(PZ5(I))
              PESG(NSTRG, NPART) = sngl(E5(I))
              PMSG(NSTRG, NPART) = sngl(XMASS5(I))
           END IF
 1005   CONTINUE
cbz1/25/99end
clin-4/09/01  turn on fragmentation with soft radiation
c     and jet order reversal to form hadrons after ZPC:
        MSTJ(1)=1
        IHPR2(1)=1
        isflag=1
clin-4/13/01 allow small mass strings (D=1.5GeV):
        HIPR1(1)=0.94
cbz2/4/99
        CALL HJANA2
cbz2/4/99end
clin-4/19/01-soft3, fragment strings, then convert hadrons to partons
c     and input to ZPC:
        elseif(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
clin-4/24/01 normal fragmentation first:
        isflag=0
c        write(99,*) 'IAEVT,NSG,NDR=',IAEVT,NSG,NDR
        IF(IHPR2(20).NE.0) THEN
           DO 560 ISG=1,NSG
                CALL HIJFRG(ISG,3,IERROR)
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
 551                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  551
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                nsbstR=0
                DO 560 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 560
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=20
                   KATT(NATT,4)=K(I,1)
c     from Yasushi, to avoid violation of array limits:
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008 to avoid out-of-bound error in K():
c                   IF(K(I,3).EQ.0 .OR.
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   GXAR(NATT) = 0.5 * (YP(1, IASG(ISG, 1)) +
     &                YT(1, IASG(ISG, 2)))
                   GYAR(NATT) = 0.5 * (YP(2, IASG(ISG, 1)) +
     &                YT(2, IASG(ISG, 2)))
                   GZAR(NATT) = 0.
                   FTAR(NATT) = 0.
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
clin-8/2015: record hadron information, to be used for its constituent partons:
                   xstrg0(NATT)=dble(GXAR(NATT))
                   ystrg0(NATT)=dble(GYAR(NATT))
                   istrg0(NATT)=ISG
c                   write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT),
c     1                  K(I,2),P(I, 1),P(I, 2),P(I, 3)
cbz11/11/98end
 560            CONTINUE
C                ********Fragment the q-qbar jets systems *****
C
           JTP(1)=IHNT2(1)
           JTP(2)=IHNT2(3)
           DO 600 NTP=1,2
           DO 600 jjtp=1,JTP(NTP)
                CALL HIJFRG(jjtp,NTP,IERROR)
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
 581                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  581
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                NFTP=NFP(jjtp,5)
                IF(NTP.EQ.2) NFTP=10+NFT(jjtp,5)
                nsbstR=0
                DO 590 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 590
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=NFTP
                   KATT(NATT,4)=K(I,1)
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008
c                   IF(K(I,3).EQ.0 .OR.
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   IF (NTP .EQ. 1) THEN
clin-2/2012:
c                      GXAR(NATT) = YP(1, jjtp)+0.5 * BB
c                      GYAR(NATT) = YP(2, jjtp)
                      GXAR(NATT) = YP(1, jjtp)+0.5*BB*cos(phiRP)
                      GYAR(NATT) = YP(2, jjtp)+0.5*BB*sin(phiRP)
                   ELSE
clin-2/2012:
c                      GXAR(NATT) = YT(1, jjtp)-0.5 * BB
c                      GYAR(NATT) = YT(2, jjtp)
                      GXAR(NATT) = YT(1, jjtp)-0.5*BB*cos(phiRP)
                      GYAR(NATT) = YT(2, jjtp)-0.5*BB*sin(phiRP)
                   END IF
                   GZAR(NATT) = 0.
                   FTAR(NATT) = 0.
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
clin-8/2015: record hadron information, to be used for its constituent partons:
                   xstrg0(NATT)=dble(GXAR(NATT))
                   ystrg0(NATT)=dble(GYAR(NATT))
c     String ID is separated for projectile/target strings:
                   istrg0(NATT)=NTP*10000+jjtp
c              if(N.eq.nsbst.and.(K(I,2).eq.2112.or.K(I,2).eq.2212)) then
c                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
c     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3),'spectator'
c                   else
c                      write(99,*) xstrg0(NATT),ystrg0(NATT),istrg0(NATT)
c     1                  ,K(I,2),P(I, 1),P(I, 2),P(I, 3)
c                   endif
cbz11/11/98end
 590                CONTINUE
 600           CONTINUE
C     ********Fragment the q-qq related string systems
        ENDIF
clin-4/2008 check for zero NDR value:
        if(NDR.ge.1) then
c
        DO 650 I=1,NDR
                NATT=NATT+1
                KATT(NATT,1)=KFDR(I)
                KATT(NATT,2)=40
                KATT(NATT,3)=0
                PATT(NATT,1)=PDR(I,1)
                PATT(NATT,2)=PDR(I,2)
                PATT(NATT,3)=PDR(I,3)
                PATT(NATT,4)=PDR(I,4)
                EATT=EATT+PDR(I,4)
clin-11/11/03     set direct photons positions and time at formation:
                GXAR(NATT) = rtdr(I,1)
                GYAR(NATT) = rtdr(I,2)
                GZAR(NATT) = 0.
                FTAR(NATT) = 0.
                ITYPAR(NATT) =KATT(NATT,1)
                PXAR(NATT) = PATT(NATT,1)
                PYAR(NATT) = PATT(NATT,2)
                PZAR(NATT) = PATT(NATT,3)
                PEAR(NATT) = PATT(NATT,4)
                XMAR(NATT) = PDR(I,5)
 650        CONTINUE
clin-4/2008:
         endif
clin-6/2009
         call embedHighPt
c
        CALL HJANA1
clin-4/19/01 convert hadrons to partons for ZPC (with GX0 given):
        call htop
clin-7/03/01 move up, used in zpstrg (otherwise not set and incorrect):
        nsp=0
        nst=0
        nsg=natt
        NSI=NSG
clin-7/03/01-end
clin-6/2009:
        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
c.....call ZPC for parton cascade
        CALL ZPCMN
clin-6/2009:
c        WRITE (14, 395) ITEST, MUL, bimp, NELP,NINP,NELT,NINTHJ
        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        itest=itest+1
        DO 1016 I = 1, MUL
c           WRITE (14, 511) PX5(I), PY5(I), PZ5(I), ITYP5(I),
c     &        XMASS5(I), E5(I)
clin-4/2012 write parton freeze-out position in zpc.dat
c     for string melting version:
c           WRITE (14, 512) ITYP5(I), PX5(I), PY5(I), PZ5(I),
c     &        XMASS5(I), LSTRG1(I), LPART1(I), FT5(I)
           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
     1          .lt.9999) then
              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           else
              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           endif
c
 1016   CONTINUE
c 511    FORMAT(1X, 3F10.4, I6, 2F10.4)
c 512    FORMAT(I6,4(1X,F10.3),1X,I6,1X,I3,1X,F10.3)
c 513    FORMAT(1X, 4F10.4)
clin-5/2009 ctest off:
c        call frztm(1,1)
clin  save data after ZPC for fragmentation purpose:
c.....transfer data back from ZPC to HIJING
        DO 1018 I = 1, MAXSTR
           DO 1017 J = 1, 3
              K1SGS(I, J) = 0
              K2SGS(I, J) = 0
              PXSGS(I, J) = 0d0
              PYSGS(I, J) = 0d0
              PZSGS(I, J) = 0d0
              PESGS(I, J) = 0d0
              PMSGS(I, J) = 0d0
              GXSGS(I, J) = 0d0
              GYSGS(I, J) = 0d0
              GZSGS(I, J) = 0d0
              FTSGS(I, J) = 0d0
 1017      CONTINUE
 1018   CONTINUE
        DO 1019 I = 1, MUL
           IITYP=ITYP5(I)
           NSTRG = LSTRG1(I)
           NPART = LPART1(I)
           K2SGS(NSTRG, NPART) = ITYP5(I)
           PXSGS(NSTRG, NPART) = PX5(I)
           PYSGS(NSTRG, NPART) = PY5(I)
           PZSGS(NSTRG, NPART) = PZ5(I)
           PMSGS(NSTRG, NPART) = XMASS5(I)
clin-7/20/01 E5(I) does no include the finite parton mass XMASS5(I),
c     so define it anew:
c           PESGS(NSTRG, NPART) = E5(I)
c           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0)
c     1          write(91,*) 'a',PX5(i),PY5(i),XMASS5(i),PZ5(i),E5(i)
           E5(I)=dsqrt(PX5(I)**2+PY5(I)**2+PZ5(I)**2+XMASS5(I)**2)
           PESGS(NSTRG, NPART) = E5(I)
c           if(abs(PZ5(i)/E5(i)).gt.0.9999999d0)
c     1          write(91,*) 'b: new E5(I)=',E5(i)
clin-7/20/01-end
           GXSGS(NSTRG, NPART) = GX5(I)
           GYSGS(NSTRG, NPART) = GY5(I)
           GZSGS(NSTRG, NPART) = GZ5(I)
           FTSGS(NSTRG, NPART) = FT5(I)
 1019   CONTINUE
        CALL HJANA2
clin-4/19/01-end
        endif
clin-4/09/01-end
C
C**************fragment all the string systems in the following*****
C
C********nsbst is where particle information starts
C********nsbstR+1 is the number of strings in fragmentation
C********the number of strings before a line is stored in K(I,4)
C********IDSTR is id number of the string system (91,92 or 93)
C
clin-4/30/01 convert partons to hadrons after ZPC:
        if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
           NATT=0
           EATT=0.
           call ptoh
           do 1006 I=1,nnozpc
              NATT=NATT+1
              KATT(NATT,1)=ITYPN(I)
              PATT(NATT,1)=PXN(I)
              PATT(NATT,2)=PYN(I)
              PATT(NATT,3)=PZN(I)
              PATT(NATT,4)=EEN(I)
              EATT=EATT+EEN(I)
              GXAR(NATT)=GXN(I)
              GYAR(NATT)=GYN(I)
              GZAR(NATT)=GZN(I)
              FTAR(NATT)=FTN(I)
              ITYPAR(NATT)=ITYPN(I)
              PXAR(NATT)=PXN(I)
              PYAR(NATT)=PYN(I)
              PZAR(NATT)=PZN(I)
              PEAR(NATT)=EEN(I)
              XMAR(NATT)=XMN(I)
 1006      continue
           goto 565
        endif
clin-4/30/01-end
        IF(IHPR2(20).NE.0) THEN
           DO 360 ISG=1,NSG
                CALL HIJFRG(ISG,3,IERROR)
                IF(MSTU(24).NE.0 .OR.IERROR.GT.0) THEN
                   MSTU(24)=0
                   MSTU(28)=0
                   IF(IHPR2(10).NE.0) THEN
c                      call lulist(2)
                      WRITE(6,*) 'error occured ISG, repeat the event'
                  write(6,*) ISG
                   ENDIF
                   GO TO 50
                ENDIF
C                        ********Check errors
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
351                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  351
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
C
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                nsbstR=0
                DO 360 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 360
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=20
                   KATT(NATT,4)=K(I,1)
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008:
c                   IF(K(I,3).EQ.0 .OR.
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
cbz11/11/98
cbz1/25/99
c                   GXAR(NATT) = 0.5 * (YP(1, IASG(ISG, 1)) +
c     &                YT(1, IASG(ISG, 2)))
c                   GYAR(NATT) = 0.5 * (YP(2, IASG(ISG, 1)) +
c     &                YT(2, IASG(ISG, 2)))
                   LSG = NSP + NST + ISG
                   GXAR(NATT) = sngl(ZT1(LSG))
                   GYAR(NATT) = sngl(ZT2(LSG))
                   GZAR(NATT) = sngl(ZT3(LSG))
                   FTAR(NATT) = sngl(ATAUI(LSG))
cbz1/25/99end
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
cbz11/11/98end
360           CONTINUE
C                ********Fragment the q-qbar jets systems *****
C
           JTP(1)=IHNT2(1)
           JTP(2)=IHNT2(3)
           DO 400 NTP=1,2
           DO 400 jjtp=1,JTP(NTP)
                CALL HIJFRG(jjtp,NTP,IERROR)
                IF(MSTU(24).NE.0 .OR. IERROR.GT.0) THEN
                   MSTU(24)=0
                   MSTU(28)=0
                   IF(IHPR2(10).NE.0) THEN
c                  call lulist(2)
                  WRITE(6,*) 'error occured P&T, repeat the event'
                  WRITE(6,*) NTP,jjtp
clin-6/2009 when this happens, the event will be repeated,
c     and another record for the same event number will be written into
c     zpc.dat, zpc.res, minijet-initial-beforePropagation.dat,
c     parton-initial-afterPropagation.dat, parton-after-coalescence.dat,
c     and parton-collisionsHistory.dat.
                   ENDIF
                   GO TO 50
                ENDIF
C                        ********check errors
C
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
381                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  381
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
C                ******** boost back to lab frame(if it was in)
C
                NFTP=NFP(jjtp,5)
                IF(NTP.EQ.2) NFTP=10+NFT(jjtp,5)
                nsbstR=0
                DO 390 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 390
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=NFTP
                   KATT(NATT,4)=K(I,1)
c                   IF(K(I,3).EQ.0 .OR. K(K(I,3),2).EQ.IDSTR) THEN
clin-4/2008:
c                   IF(K(I,3).EQ.0 .OR.
c     1 (K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR)) THEN
c                      KATT(NATT,3)=0
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
clin-4/2008-end
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
C       ****** identify the mother particle
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
cbz11/11/98
cbz1/25/99
c                   IF (NTP .EQ. 1) THEN
c                      GXAR(NATT) = YP(1, jjtp)
c                   ELSE
c                      GXAR(NATT) = YT(1, jjtp)
c                   END IF
c                   IF (NTP .EQ. 1) THEN
c                      GYAR(NATT) = YP(2, jjtp)
c                   ELSE
c                      GYAR(NATT) = YT(2, jjtp)
c                   END IF
                   IF (NTP .EQ. 1) THEN
                      LSG = jjtp
                   ELSE
                      LSG = jjtp + NSP
                   END IF
                   GXAR(NATT) = sngl(ZT1(LSG))
                   GYAR(NATT) = sngl(ZT2(LSG))
                   GZAR(NATT) = sngl(ZT3(LSG))
                   FTAR(NATT) = sngl(ATAUI(LSG))
cbz1/25/99end
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
cbz11/11/98end
390                CONTINUE
400           CONTINUE
C     ********Fragment the q-qq related string systems
        ENDIF
        DO 450 I=1,NDR
           NATT=NATT+1
           KATT(NATT,1)=KFDR(I)
           KATT(NATT,2)=40
           KATT(NATT,3)=0
           PATT(NATT,1)=PDR(I,1)
           PATT(NATT,2)=PDR(I,2)
           PATT(NATT,3)=PDR(I,3)
           PATT(NATT,4)=PDR(I,4)
           EATT=EATT+PDR(I,4)
clin-11/11/03     set direct photons positions and time at formation:
           GXAR(NATT) = rtdr(I,1)
           GYAR(NATT) = rtdr(I,2)
           GZAR(NATT) = 0.
           FTAR(NATT) = 0.
           ITYPAR(NATT) =KATT(NATT,1)
           PXAR(NATT) = PATT(NATT,1)
           PYAR(NATT) = PATT(NATT,2)
           PZAR(NATT) = PATT(NATT,3)
           PEAR(NATT) = PATT(NATT,4)
           XMAR(NATT) = PDR(I,5)
 450    CONTINUE
C                        ********store the direct-produced particles
C
clin-4/19/01 soft3:
 565    continue
        DENGY=EATT/(IHNT2(1)*HINT1(6)+IHNT2(3)*HINT1(7))-1.0
        IF(ABS(DENGY).GT.HIPR1(43).AND.IHPR2(20).NE.0
     &     .AND.IHPR2(21).EQ.0) THEN
         IF(IHPR2(10).NE.0)
     &        WRITE(6,*) 'Energy not conserved, repeat the event'
c                call lulist(1)
         write(6,*) 'violated:EATT(GeV),NATT,B(fm)=',EATT,NATT,bimp
         GO TO 50
        ENDIF
        write(6,*) 'satisfied:EATT(GeV),NATT,B(fm)=',EATT,NATT,bimp
        write(6,*) ' '
c
clin-4/2012 write out initial transverse positions of initial nucleons:
        write(94,*) IAEVT,MISS,IHNT2(1),IHNT2(3),bimp
        DO JP=1,IHNT2(1)
clin-12/2012 write out present and original flavor code of nucleons:
c           write(94,243) YP(1,JP)+0.5*BB*cos(phiRP),
c     1 YP(2,JP)+0.5*BB*sin(phiRP), JP, NFP(JP,5),yp(3,jp)
           write(94,243) YP(1,JP)+0.5*BB*cos(phiRP),
     1 YP(2,JP)+0.5*BB*sin(phiRP),JP, NFP(JP,5),yp(3,jp),
     2 NFP(JP,3),NFP(JP,4)
        ENDDO
        DO JT=1,IHNT2(3)
c target nucleon # has a minus sign for distinction from projectile:
clin-12/2012 write out present and original flavor code of nucleons:
c           write(94,243) YT(1,JT)-0.5*BB*cos(phiRP),
c     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt)
           write(94,243) YT(1,JT)-0.5*BB*cos(phiRP),
     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt),
     2 NFT(JT,3),NFT(JT,4)
        ENDDO
clin-12/2012 write out present and original flavor code of nucleons:
c 243    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3)
 243    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3,2(1x,I5))
clin-4/2012-end
        RETURN
        END
C
C
C
