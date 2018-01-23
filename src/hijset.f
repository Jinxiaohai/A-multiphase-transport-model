        SUBROUTINE HIJSET(EFRM,FRAME,PROJ,TARG,IAP,IZP,IAT,IZT)
        CHARACTER FRAME*4,PROJ*4,TARG*4,EFRAME*4
        DOUBLE PRECISION  DD1,DD2,DD3,DD4
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
        COMMON/hjcrdn/YP(3,300),YT(3,300)
cc      SAVE /hjcrdn/
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        COMMON/HIJDAT/HIDAT0(10,10),HIDAT(10)
cc      SAVE /HIJDAT/
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
cc      SAVE /LUDAT1/
        EXTERNAL FNKICK,FNKC2,FNSTRU,FNSTRM,FNSTRS
        SAVE
        CALL TITLE
        IHNT2(1)=IAP
        IHNT2(2)=IZP
        IHNT2(3)=IAT
        IHNT2(4)=IZT
        IHNT2(5)=0
        IHNT2(6)=0
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  给出弹核和靶核的静止质量hint1(8), hint1(9)，为啥取最大的质量？？
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        HINT1(8)=MAX(ULMASS(2112),ULMASS(2212))
        HINT1(9)=HINT1(8)
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  弹核和靶核种类的确定
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IF(PROJ.NE.'A') THEN
                IF(PROJ.EQ.'P') THEN
                    IHNT2(5)=2212
                ELSE IF(PROJ.EQ.'PBAR') THEN
                    IHNT2(5)=-2212
                ELSE IF(PROJ.EQ.'PI+') THEN
                    IHNT2(5)=211
                ELSE IF(PROJ.EQ.'PI-') THEN
                    IHNT2(5)=-211
                ELSE IF(PROJ.EQ.'K+') THEN
                    IHNT2(5)=321
                ELSE IF(PROJ.EQ.'K-') THEN
                    IHNT2(5)=-321
                ELSE IF(PROJ.EQ.'N') THEN
                    IHNT2(5)=2112
                ELSE IF(PROJ.EQ.'NBAR') THEN
                    IHNT2(5)=-2112
                ELSE
                    WRITE(6,*) PROJ, 'wrong or unavailable proj name'
                    STOP
                ENDIF
                HINT1(8)=ULMASS(IHNT2(5))
        ENDIF
        IF(TARG.NE.'A') THEN
                IF(TARG.EQ.'P') THEN
                    IHNT2(6)=2212
                ELSE IF(TARG.EQ.'PBAR') THEN
                    IHNT2(6)=-2212
                ELSE IF(TARG.EQ.'PI+') THEN
                    IHNT2(6)=211
                ELSE IF(TARG.EQ.'PI-') THEN
                    IHNT2(6)=-211
                ELSE IF(TARG.EQ.'K+') THEN
                    IHNT2(6)=321
                ELSE IF(TARG.EQ.'K-') THEN
                    IHNT2(6)=-321
                ELSE IF(TARG.EQ.'N') THEN
                    IHNT2(6)=2112
                ELSE IF(TARG.EQ.'NBAR') THEN
                    IHNT2(6)=-2112
                ELSE
                    WRITE(6,*) TARG,'wrong or unavailable targ name'
                    STOP
                ENDIF
                HINT1(9)=ULMASS(IHNT2(6))
        ENDIF
C...Switch off decay of pi0, K0S, Lambda, Sigma+-, Xi0-, Omega-.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  是否关掉下面这些粒子的衰变,缺省值为1.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IF(IHPR2(12).GT.0) THEN
        CALL LUGIVE('MDCY(C221,1)=0')
clin-11/07/00 no K* decays:
        CALL LUGIVE('MDCY(C313,1)=0')
        CALL LUGIVE('MDCY(C-313,1)=0')
        CALL LUGIVE('MDCY(C323,1)=0')
        CALL LUGIVE('MDCY(C-323,1)=0')
clin-1/04/01 no K0 and K0bar decays so K0L and K0S do not appear,
c     this way the K/Kbar difference is accounted for exactly:
        CALL LUGIVE('MDCY(C311,1)=0')
        CALL LUGIVE('MDCY(C-311,1)=0')
clin-11/08/00 no Delta decays:
        CALL LUGIVE('MDCY(C1114,1)=0')
        CALL LUGIVE('MDCY(C2114,1)=0')
        CALL LUGIVE('MDCY(C2214,1)=0')
        CALL LUGIVE('MDCY(C2224,1)=0')
        CALL LUGIVE('MDCY(C-1114,1)=0')
        CALL LUGIVE('MDCY(C-2114,1)=0')
        CALL LUGIVE('MDCY(C-2214,1)=0')
        CALL LUGIVE('MDCY(C-2224,1)=0')
clin-11/07/00-end
cbz12/4/98
        CALL LUGIVE('MDCY(C213,1)=0')
        CALL LUGIVE('MDCY(C-213,1)=0')
        CALL LUGIVE('MDCY(C113,1)=0')
        CALL LUGIVE('MDCY(C223,1)=0')
        CALL LUGIVE('MDCY(C333,1)=0')
cbz12/4/98end
        CALL LUGIVE('MDCY(C111,1)=0')
        CALL LUGIVE('MDCY(C310,1)=0')
        CALL LUGIVE('MDCY(C411,1)=0;MDCY(C-411,1)=0')
        CALL LUGIVE('MDCY(C421,1)=0;MDCY(C-421,1)=0')
        CALL LUGIVE('MDCY(C431,1)=0;MDCY(C-431,1)=0')
        CALL LUGIVE('MDCY(C511,1)=0;MDCY(C-511,1)=0')
        CALL LUGIVE('MDCY(C521,1)=0;MDCY(C-521,1)=0')
        CALL LUGIVE('MDCY(C531,1)=0;MDCY(C-531,1)=0')
        CALL LUGIVE('MDCY(C3122,1)=0;MDCY(C-3122,1)=0')
        CALL LUGIVE('MDCY(C3112,1)=0;MDCY(C-3112,1)=0')
        CALL LUGIVE('MDCY(C3212,1)=0;MDCY(C-3212,1)=0')
        CALL LUGIVE('MDCY(C3222,1)=0;MDCY(C-3222,1)=0')
        CALL LUGIVE('MDCY(C3312,1)=0;MDCY(C-3312,1)=0')
        CALL LUGIVE('MDCY(C3322,1)=0;MDCY(C-3322,1)=0')
        CALL LUGIVE('MDCY(C3334,1)=0;MDCY(C-3334,1)=0')
clin-7/2011-no HQ(charm or bottom) decays in order to get net-HQ conservation:
        CALL LUGIVE('MDCY(C441,1)=0')
        CALL LUGIVE('MDCY(C443,1)=0')
        CALL LUGIVE('MDCY(C413,1)=0;MDCY(C-413,1)=0')
        CALL LUGIVE('MDCY(C423,1)=0;MDCY(C-423,1)=0')
        CALL LUGIVE('MDCY(C433,1)=0;MDCY(C-433,1)=0')
        CALL LUGIVE('MDCY(C4112,1)=0;MDCY(C-4112,1)=0')
        CALL LUGIVE('MDCY(C4114,1)=0;MDCY(C-4114,1)=0')
        CALL LUGIVE('MDCY(C4122,1)=0;MDCY(C-4122,1)=0')
        CALL LUGIVE('MDCY(C4212,1)=0;MDCY(C-4212,1)=0')
        CALL LUGIVE('MDCY(C4214,1)=0;MDCY(C-4214,1)=0')
        CALL LUGIVE('MDCY(C4222,1)=0;MDCY(C-4222,1)=0')
        CALL LUGIVE('MDCY(C4224,1)=0;MDCY(C-4224,1)=0')
        CALL LUGIVE('MDCY(C4132,1)=0;MDCY(C-4132,1)=0')
        CALL LUGIVE('MDCY(C4312,1)=0;MDCY(C-4312,1)=0')
        CALL LUGIVE('MDCY(C4314,1)=0;MDCY(C-4314,1)=0')
        CALL LUGIVE('MDCY(C4232,1)=0;MDCY(C-4232,1)=0')
        CALL LUGIVE('MDCY(C4322,1)=0;MDCY(C-4322,1)=0')
        CALL LUGIVE('MDCY(C4324,1)=0;MDCY(C-4324,1)=0')
        CALL LUGIVE('MDCY(C4332,1)=0;MDCY(C-4332,1)=0')
        CALL LUGIVE('MDCY(C4334,1)=0;MDCY(C-4334,1)=0')
        CALL LUGIVE('MDCY(C551,1)=0')
        CALL LUGIVE('MDCY(C553,1)=0')
        CALL LUGIVE('MDCY(C513,1)=0;MDCY(C-513,1)=0')
        CALL LUGIVE('MDCY(C523,1)=0;MDCY(C-523,1)=0')
        CALL LUGIVE('MDCY(C533,1)=0;MDCY(C-533,1)=0')
        CALL LUGIVE('MDCY(C5112,1)=0;MDCY(C-5112,1)=0')
        CALL LUGIVE('MDCY(C5114,1)=0;MDCY(C-5114,1)=0')
        CALL LUGIVE('MDCY(C5122,1)=0;MDCY(C-5122,1)=0')
        CALL LUGIVE('MDCY(C5212,1)=0;MDCY(C-5212,1)=0')
        CALL LUGIVE('MDCY(C5214,1)=0;MDCY(C-5214,1)=0')
        CALL LUGIVE('MDCY(C5222,1)=0;MDCY(C-5222,1)=0')
        CALL LUGIVE('MDCY(C5224,1)=0;MDCY(C-5224,1)=0')
clin-7/2011-end
        ENDIF
        MSTU(12)=0
        MSTU(21)=1
        IF(IHPR2(10).EQ.0) THEN
                MSTU(22)=0
                MSTU(25)=0
                MSTU(26)=0
        ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  fragmentation parameter a, b.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
clin    parj(41) and (42) are a, b parameters in Lund, read from input.ampt:
c        PARJ(41)=HIPR1(3)
c        PARJ(42)=HIPR1(4)
c        PARJ(41)=2.2
c        PARJ(42)=0.5
clin  2 popcorn parameters read from input.ampt:
c        IHPR2(11) = 3
c        PARJ(5) = 0.5
        MSTJ(12)=IHPR2(11)
clin  parj(21) gives the mean gaussian width for hadron Pt:
        PARJ(21)=HIPR1(2)
clin  parj(2) is gamma_s=P(s)/P(u), kappa propto 1/b/(2+a) assumed.
        rkp=HIPR1(4)*(2+HIPR1(3))/PARJ(42)/(2+PARJ(41))
        PARJ(2)=PARJ(2)**(1./rkp)
        PARJ(21)=PARJ(21)*sqrt(rkp)
clin-10/31/00 update when string tension is changed:
        HIPR1(2)=PARJ(21)
clin-4/2015: set upper limit for gamma_s=P(s)/P(u) to 0.4
c     (to limit strangeness enhancement when string tension is strongly
c     increased due to using a very low value of parameter b in Lund
c     symmetric splitting function as done in arXiv:1403.6321):
        PARJ(2)=min(PARJ(2),0.4)
C                        ******** set up for jetset
        IF(FRAME.EQ.'LAB') THEN
           DD1=dble(EFRM)
           DD2=dble(HINT1(8))
           DD3=dble(HINT1(9))
           HINT1(1)=SQRT(HINT1(8)**2+2.0*HINT1(9)*EFRM+HINT1(9)**2)
           DD4=DSQRT(DD1**2-DD2**2)/(DD1+DD3)
           HINT1(2)=sngl(DD4)
           HINT1(3)=0.5*sngl(DLOG((1.D0+DD4)/(1.D0-DD4)))
           DD4=DSQRT(DD1**2-DD2**2)/DD1
           HINT1(4)=0.5*sngl(DLOG((1.D0+DD4)/(1.D0-DD4)))
           HINT1(5)=0.0
           HINT1(6)=EFRM
           HINT1(7)=HINT1(9)
        ELSE IF(FRAME.EQ.'CMS') THEN
           HINT1(1)=EFRM
           HINT1(2)=0.0
           HINT1(3)=0.0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  质心系中每核子核子碰撞的能量已经静止质量。
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           DD1=dble(HINT1(1))
           DD2=dble(HINT1(8))
           DD3=dble(HINT1(9))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  快度计算
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           DD4=DSQRT(1.D0-4.D0*DD2**2/DD1**2)
           HINT1(4)=0.5*sngl(DLOG((1.D0+DD4)/(1.D0-DD4)))
           DD4=DSQRT(1.D0-4.D0*DD3**2/DD1**2)
           HINT1(5)=-0.5*sngl(DLOG((1.D0+DD4)/(1.D0-DD4)))
           HINT1(6)=HINT1(1)/2.0
           HINT1(7)=HINT1(1)/2.0
        ENDIF
C                ********define Lorentz transform to lab frame
c
C                ********calculate the cross sections involved with
C                        nucleon collisions.
        IF(IHNT2(1).GT.1) THEN
                CALL HIJWDS(IHNT2(1),1,RMAX)
                HIPR1(34)=RMAX
C                        ********set up Wood-Sax distr for proj.
        ENDIF
        IF(IHNT2(3).GT.1) THEN
                CALL HIJWDS(IHNT2(3),2,RMAX)
                HIPR1(35)=RMAX
C                        ********set up Wood-Sax distr for  targ.
        ENDIF
C
C
        I=0
20        I=I+1
        IF(I.EQ.10) GO TO 30
        IF(HIDAT0(10,I).LE.HINT1(1)) GO TO 20
30        IF(I.EQ.1) I=2
        DO 40 J=1,9
           HIDAT(J)=HIDAT0(J,I-1)+(HIDAT0(J,I)-HIDAT0(J,I-1))
     &          *(HINT1(1)-HIDAT0(10,I-1))/(HIDAT0(10,I)-HIDAT0(10,I-1))
40        CONTINUE
        HIPR1(31)=HIDAT(5)
        HIPR1(30)=2.0*HIDAT(5)
C
C
        CALL HIJCRS
C
        IF(IHPR2(5).NE.0) THEN
                CALL HIFUN(3,0.0,36.0,FNKICK)
C                ********booking for generating pt**2 for pt kick
        ENDIF
        CALL HIFUN(7,0.0,6.0,FNKC2)
        CALL HIFUN(4,0.0,1.0,FNSTRU)
        CALL HIFUN(5,0.0,1.0,FNSTRM)
        CALL HIFUN(6,0.0,1.0,FNSTRS)
C                ********booking for x distribution of valence quarks
        EFRAME='Ecm'
        IF(FRAME.EQ.'LAB') EFRAME='Elab'
        WRITE(6,100) EFRAME,EFRM,PROJ,IHNT2(1),IHNT2(2),
     &               TARG,IHNT2(3),IHNT2(4)
100        FORMAT(
     &        10X,'**************************************************'/
     &        10X,'*',48X,'*'/
     &        10X,'*         HIJING has been initialized at         *'/
     &        10X,'*',13X,A4,'= ',F10.2,' GeV/n',13X,'*'/
     &        10X,'*',48X,'*'/
     &        10X,'*',8X,'for ',
     &        A4,'(',I3,',',I3,')',' + ',A4,'(',I3,',',I3,')',7X,'*'/
     &        10X,'**************************************************')
        RETURN
        END
C
C
C
