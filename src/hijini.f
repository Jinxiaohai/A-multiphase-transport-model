        SUBROUTINE HIJINI
        PARAMETER (MAXSTR=150001)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
cc      SAVE /HPARNT/
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
cc      SAVE /HSTRNG/
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
cc      SAVE /HJJET1/
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
c        COMMON/HJJET4/NDR,IADR(900,2),KFDR(900),PDR(900,5)
        COMMON/HJJET4/NDR,IADR(MAXSTR,2),KFDR(MAXSTR),PDR(MAXSTR,5)
cc      SAVE /HJJET4/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
        SAVE
C****************Reset the momentum of initial particles************
C             and assign flavors to the proj and targ string       *
C*******************************************************************
        NSG=0
        NDR=0
        IPP=2212
        IPT=2212
        IF(IHNT2(5).NE.0) IPP=IHNT2(5)
        IF(IHNT2(6).NE.0) IPT=IHNT2(6)
C                ********in case the proj or targ is a hadron.
C


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  遍历弹核的每个核子
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        DO 100 I=1,IHNT2(1)
        PP(I,1)=0.0
        PP(I,2)=0.0
        PP(I,3)=SQRT(HINT1(1)**2/4.0-HINT1(8)**2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  能量和质量
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        PP(I,4)=HINT1(1)/2
        PP(I,5)=HINT1(8)
        PP(I,6)=0.0
        PP(I,7)=0.0
        PP(I,8)=0.0
        PP(I,9)=0.0
        PP(I,10)=0.0
cbzdbg2/22/99
ctest OFF
        PP(I, 11) = 0.0
        PP(I, 12) = 0.0
cbzdbg2/22/99end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  NFP(I,3)是present ID, 而NFP(I,4)是original ID. 先初始质子的。
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        NFP(I,3)=IPP
        NFP(I,4)=IPP
        NFP(I,5)=0
        NFP(I,6)=0
        NFP(I,7)=0
        NFP(I,8)=0
        NFP(I,9)=0
        NFP(I,10)=0
        NFP(I,11)=0
        NPJ(I)=0
        IF(I.GT.ABS(IHNT2(2))) NFP(I,3)=2112
clin-12/2012 correct NN differential cross section in HIJING:
        IF(I.GT.ABS(IHNT2(2))) NFP(I,4)=2112
        CALL ATTFLV(NFP(I,3),IDQ,IDQQ)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  NFP(I,1)是弹核中夸克的味道，NFP(I,2)是弹核中diquark的味道。
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        NFP(I,1)=IDQ
        NFP(I,2)=IDQQ
        NFP(I,15)=-1
        IF(ABS(IDQ).GT.1000.OR.(ABS(IDQ*IDQQ).LT.100.AND.
     &                RANART(NSEED).LT.0.5)) NFP(I,15)=1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  PP(I, 14):the mass of the projectile quark.
cccccc  PP(I, 15):the mass of the projectile diquark.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        PP(I,14)=ULMASS(IDQ)
        PP(I,15)=ULMASS(IDQQ)
100        CONTINUE
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  同理初始化靶核的信息
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        DO 200 I=1,IHNT2(3)
        PT(I,1)=0.0
        PT(I,2)=0.0
        PT(I,3)=-SQRT(HINT1(1)**2/4.0-HINT1(9)**2)
        PT(I,4)=HINT1(1)/2.0
        PT(I,5)=HINT1(9)
        PT(I,6)=0.0
        PT(I,7)=0.0
        PT(I,8)=0.0
        PT(I,9)=0.0
        PT(I,10)=0.0
ctest OFF
cbzdbg2/22/99
        PT(I, 11) = 0.0
        PT(I, 12) = 0.0
cbzdbg2/22/99end
        NFT(I,3)=IPT
        NFT(I,4)=IPT
        NFT(I,5)=0
        NFT(I,6)=0
        NFT(I,7)=0
        NFT(I,8)=0
        NFT(I,9)=0
        NFT(I,10)=0
        NFT(I,11)=0
        NTJ(I)=0
        IF(I.GT.ABS(IHNT2(4))) NFT(I,3)=2112
clin-12/2012 correct NN differential cross section in HIJING:
        IF(I.GT.ABS(IHNT2(4))) NFT(I,4)=2112
        CALL ATTFLV(NFT(I,3),IDQ,IDQQ)
        NFT(I,1)=IDQ
        NFT(I,2)=IDQQ
        NFT(I,15)=1
        IF(ABS(IDQ).GT.1000.OR.(ABS(IDQ*IDQQ).LT.100.AND.
     &       RANART(NSEED).LT.0.5)) NFT(I,15)=-1
        PT(I,14)=ULMASS(IDQ)
        PT(I,15)=ULMASS(IDQQ)
200        CONTINUE
        RETURN
        END
C
C
C
