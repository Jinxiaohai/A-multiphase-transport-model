      SUBROUTINE PYKLIM(ILIM)   
C...Checks generated variables against pre-set kinematical limits;  
C...also calculates limits on variables used in generation. 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
      SAVE /LUDAT3/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200) 
      SAVE /PYSUBS/ 
      COMMON/PYINT1/MINT(400),VINT(400) 
      SAVE /PYINT1/ 
      COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2) 
      SAVE /PYINT2/ 
C...Common kinematical expressions. 
      ISUB=MINT(1)  
      IF(ISUB.EQ.96) GOTO 110   
      SQM3=VINT(63) 
      SQM4=VINT(64) 
      IF(ILIM.NE.1) THEN    
        TAU=VINT(21)    
        RM3=SQM3/(TAU*VINT(2))  
        RM4=SQM4/(TAU*VINT(2))  
        BE34=SQRT((1.-RM3-RM4)**2-4.*RM3*RM4)   
      ENDIF 
      PTHMIN=CKIN(3)    
      IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) PTHMIN=MAX(CKIN(3),CKIN(5))  
      IF(ILIM.EQ.0) THEN    
C...Check generated values of tau, y*, cos(theta-hat), and tau' against 
C...pre-set kinematical limits. 
        YST=VINT(22)    
        CTH=VINT(23)    
        TAUP=VINT(26)   
        IF(ISET(ISUB).LE.2) THEN    
          X1=SQRT(TAU)*EXP(YST) 
          X2=SQRT(TAU)*EXP(-YST)    
        ELSE    
          X1=SQRT(TAUP)*EXP(YST)    
          X2=SQRT(TAUP)*EXP(-YST)   
        ENDIF   
        XF=X1-X2    
        IF(TAU*VINT(2).LT.CKIN(1)**2) MINT(51)=1    
        IF(CKIN(2).GE.0..AND.TAU*VINT(2).GT.CKIN(2)**2) MINT(51)=1  
        IF(X1.LT.CKIN(21).OR.X1.GT.CKIN(22)) MINT(51)=1 
        IF(X2.LT.CKIN(23).OR.X2.GT.CKIN(24)) MINT(51)=1 
        IF(XF.LT.CKIN(25).OR.XF.GT.CKIN(26)) MINT(51)=1 
        IF(YST.LT.CKIN(7).OR.YST.GT.CKIN(8)) MINT(51)=1 
        IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN 
          PTH=0.5*BE34*SQRT(TAU*VINT(2)*(1.-CTH**2))    
          Y3=YST+0.5*LOG((1.+RM3-RM4+BE34*CTH)/(1.+RM3-RM4-BE34*CTH))   
          Y4=YST+0.5*LOG((1.+RM4-RM3-BE34*CTH)/(1.+RM4-RM3+BE34*CTH))   
          YLARGE=MAX(Y3,Y4) 
          YSMALL=MIN(Y3,Y4) 
          ETALAR=10.    
          ETASMA=-10.   
          STH=SQRT(1.-CTH**2)   
          IF(STH.LT.1.E-6) GOTO 100 
          EXPET3=((1.+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH+    
     &    SQRT(((1.+RM3-RM4)*COSH(YST)+BE34*SINH(YST)*CTH)**2-4.*RM3))/ 
     &    (BE34*STH)    
          EXPET4=((1.-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH+    
     &    SQRT(((1.-RM3+RM4)*COSH(YST)-BE34*SINH(YST)*CTH)**2-4.*RM4))/ 
     &    (BE34*STH)    
          ETA3=LOG(MIN(1.E10,MAX(1.E-10,EXPET3)))   
          ETA4=LOG(MIN(1.E10,MAX(1.E-10,EXPET4)))   
          ETALAR=MAX(ETA3,ETA4) 
          ETASMA=MIN(ETA3,ETA4) 
  100     CTS3=((1.+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH)/ 
     &    SQRT(((1.+RM3-RM4)*COSH(YST)+BE34*SINH(YST)*CTH)**2-4.*RM3)   
          CTS4=((1.-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH)/ 
     &    SQRT(((1.-RM3+RM4)*COSH(YST)-BE34*SINH(YST)*CTH)**2-4.*RM4)   
          CTSLAR=MAX(CTS3,CTS4) 
          CTSSMA=MIN(CTS3,CTS4) 
          IF(PTH.LT.PTHMIN) MINT(51)=1  
          IF(CKIN(4).GE.0..AND.PTH.GT.CKIN(4)) MINT(51)=1   
          IF(YLARGE.LT.CKIN(9).OR.YLARGE.GT.CKIN(10)) MINT(51)=1    
          IF(YSMALL.LT.CKIN(11).OR.YSMALL.GT.CKIN(12)) MINT(51)=1   
          IF(ETALAR.LT.CKIN(13).OR.ETALAR.GT.CKIN(14)) MINT(51)=1   
          IF(ETASMA.LT.CKIN(15).OR.ETASMA.GT.CKIN(16)) MINT(51)=1   
          IF(CTSLAR.LT.CKIN(17).OR.CTSLAR.GT.CKIN(18)) MINT(51)=1   
          IF(CTSSMA.LT.CKIN(19).OR.CTSSMA.GT.CKIN(20)) MINT(51)=1   
          IF(CTH.LT.CKIN(27).OR.CTH.GT.CKIN(28)) MINT(51)=1 
        ENDIF   
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) THEN 
          IF(TAUP*VINT(2).LT.CKIN(31)**2) MINT(51)=1    
          IF(CKIN(32).GE.0..AND.TAUP*VINT(2).GT.CKIN(32)**2) MINT(51)=1 
        ENDIF   
      ELSEIF(ILIM.EQ.1) THEN    
C...Calculate limits on tau 
C...0) due to definition    
        TAUMN0=0.   
        TAUMX0=1.   
C...1) due to limits on subsystem mass  
        TAUMN1=CKIN(1)**2/VINT(2)   
        TAUMX1=1.   
        IF(CKIN(2).GE.0.) TAUMX1=CKIN(2)**2/VINT(2) 
C...2) due to limits on pT-hat (and non-overlapping rapidity intervals) 
        TM3=SQRT(SQM3+PTHMIN**2)    
        TM4=SQRT(SQM4+PTHMIN**2)    
        YDCOSH=1.   
        IF(CKIN(9).GT.CKIN(12)) YDCOSH=COSH(CKIN(9)-CKIN(12))   
        TAUMN2=(TM3**2+2.*TM3*TM4*YDCOSH+TM4**2)/VINT(2)    
        TAUMX2=1.   
C...3) due to limits on pT-hat and cos(theta-hat)   
        CTH2MN=MIN(CKIN(27)**2,CKIN(28)**2) 
        CTH2MX=MAX(CKIN(27)**2,CKIN(28)**2) 
        TAUMN3=0.   
        IF(CKIN(27)*CKIN(28).GT.0.) TAUMN3= 
     &  (SQRT(SQM3+PTHMIN**2/(1.-CTH2MN))+  
     &  SQRT(SQM4+PTHMIN**2/(1.-CTH2MN)))**2/VINT(2)    
        TAUMX3=1.   
        IF(CKIN(4).GE.0..AND.CTH2MX.LT.1.) TAUMX3=  
     &  (SQRT(SQM3+CKIN(4)**2/(1.-CTH2MX))+ 
     &  SQRT(SQM4+CKIN(4)**2/(1.-CTH2MX)))**2/VINT(2)   
C...4) due to limits on x1 and x2   
        TAUMN4=CKIN(21)*CKIN(23)    
        TAUMX4=CKIN(22)*CKIN(24)    
C...5) due to limits on xF  
        TAUMN5=0.   
        TAUMX5=MAX(1.-CKIN(25),1.+CKIN(26)) 
        VINT(11)=MAX(TAUMN0,TAUMN1,TAUMN2,TAUMN3,TAUMN4,TAUMN5) 
        VINT(31)=MIN(TAUMX0,TAUMX1,TAUMX2,TAUMX3,TAUMX4,TAUMX5) 
        IF(MINT(43).EQ.1.AND.(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.2)) THEN 
          VINT(11)=0.99999  
          VINT(31)=1.00001  
        ENDIF   
        IF(VINT(31).LE.VINT(11)) MINT(51)=1 
      ELSEIF(ILIM.EQ.2) THEN    
C...Calculate limits on y*  
        IF(ISET(ISUB).EQ.3.OR.ISET(ISUB).EQ.4) TAU=VINT(26) 
        TAURT=SQRT(TAU) 
C...0) due to kinematics    
        YSTMN0=LOG(TAURT)   
        YSTMX0=-YSTMN0  
C...1) due to explicit limits   
        YSTMN1=CKIN(7)  
        YSTMX1=CKIN(8)  
C...2) due to limits on x1  
        YSTMN2=LOG(MAX(TAU,CKIN(21))/TAURT) 
        YSTMX2=LOG(MAX(TAU,CKIN(22))/TAURT) 
C...3) due to limits on x2  
        YSTMN3=-LOG(MAX(TAU,CKIN(24))/TAURT)    
        YSTMX3=-LOG(MAX(TAU,CKIN(23))/TAURT)    
C...4) due to limits on xF  
        YEPMN4=0.5*ABS(CKIN(25))/TAURT  
        YSTMN4=SIGN(LOG(SQRT(1.+YEPMN4**2)+YEPMN4),CKIN(25))    
        YEPMX4=0.5*ABS(CKIN(26))/TAURT  
        YSTMX4=SIGN(LOG(SQRT(1.+YEPMX4**2)+YEPMX4),CKIN(26))    
C...5) due to simultaneous limits on y-large and y-small    
        YEPSMN=(RM3-RM4)*SINH(CKIN(9)-CKIN(11)) 
        YEPSMX=(RM3-RM4)*SINH(CKIN(10)-CKIN(12))    
        YDIFMN=ABS(LOG(SQRT(1.+YEPSMN**2)-YEPSMN))  
        YDIFMX=ABS(LOG(SQRT(1.+YEPSMX**2)-YEPSMX))  
        YSTMN5=0.5*(CKIN(9)+CKIN(11)-YDIFMN)    
        YSTMX5=0.5*(CKIN(10)+CKIN(12)+YDIFMX)   
C...6) due to simultaneous limits on cos(theta-hat) and y-large or  
C...   y-small  
        CTHLIM=SQRT(1.-4.*PTHMIN**2/(BE34*TAU*VINT(2))) 
        RZMN=BE34*MAX(CKIN(27),-CTHLIM) 
        RZMX=BE34*MIN(CKIN(28),CTHLIM)  
        YEX3MX=(1.+RM3-RM4+RZMX)/MAX(1E-10,1.+RM3-RM4-RZMX) 
        YEX4MX=(1.+RM4-RM3-RZMN)/MAX(1E-10,1.+RM4-RM3+RZMN) 
        YEX3MN=MAX(1E-10,1.+RM3-RM4+RZMN)/(1.+RM3-RM4-RZMN) 
        YEX4MN=MAX(1E-10,1.+RM4-RM3-RZMX)/(1.+RM4-RM3+RZMX) 
        YSTMN6=CKIN(9)-0.5*LOG(MAX(YEX3MX,YEX4MX))  
        YSTMX6=CKIN(12)-0.5*LOG(MIN(YEX3MN,YEX4MN)) 
        VINT(12)=MAX(YSTMN0,YSTMN1,YSTMN2,YSTMN3,YSTMN4,YSTMN5,YSTMN6)  
        VINT(32)=MIN(YSTMX0,YSTMX1,YSTMX2,YSTMX3,YSTMX4,YSTMX5,YSTMX6)  
        IF(MINT(43).EQ.1) THEN  
          VINT(12)=-0.00001 
          VINT(32)=0.00001  
        ELSEIF(MINT(43).EQ.2) THEN  
          VINT(12)=0.99999*YSTMX0   
          VINT(32)=1.00001*YSTMX0   
        ELSEIF(MINT(43).EQ.3) THEN  
          VINT(12)=-1.00001*YSTMX0  
          VINT(32)=-0.99999*YSTMX0  
        ENDIF   
        IF(VINT(32).LE.VINT(12)) MINT(51)=1 
      ELSEIF(ILIM.EQ.3) THEN    
C...Calculate limits on cos(theta-hat)  
        YST=VINT(22)    
C...0) due to definition    
        CTNMN0=-1.  
        CTNMX0=0.   
        CTPMN0=0.   
        CTPMX0=1.   
C...1) due to explicit limits   
        CTNMN1=MIN(0.,CKIN(27)) 
        CTNMX1=MIN(0.,CKIN(28)) 
        CTPMN1=MAX(0.,CKIN(27)) 
        CTPMX1=MAX(0.,CKIN(28)) 
C...2) due to limits on pT-hat  
        CTNMN2=-SQRT(1.-4.*PTHMIN**2/(BE34**2*TAU*VINT(2))) 
        CTPMX2=-CTNMN2  
        CTNMX2=0.   
        CTPMN2=0.   
        IF(CKIN(4).GE.0.) THEN  
          CTNMX2=-SQRT(MAX(0.,1.-4.*CKIN(4)**2/(BE34**2*TAU*VINT(2))))  
          CTPMN2=-CTNMX2    
        ENDIF   
C...3) due to limits on y-large and y-small 
        CTNMN3=MIN(0.,MAX((1.+RM3-RM4)/BE34*TANH(CKIN(11)-YST), 
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(10)-YST))) 
        CTNMX3=MIN(0.,(1.+RM3-RM4)/BE34*TANH(CKIN(12)-YST), 
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(9)-YST))   
        CTPMN3=MAX(0.,(1.+RM3-RM4)/BE34*TANH(CKIN(9)-YST),  
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(12)-YST))  
        CTPMX3=MAX(0.,MIN((1.+RM3-RM4)/BE34*TANH(CKIN(10)-YST), 
     &  -(1.-RM3+RM4)/BE34*TANH(CKIN(11)-YST))) 
        VINT(13)=MAX(CTNMN0,CTNMN1,CTNMN2,CTNMN3)   
        VINT(33)=MIN(CTNMX0,CTNMX1,CTNMX2,CTNMX3)   
        VINT(14)=MAX(CTPMN0,CTPMN1,CTPMN2,CTPMN3)   
        VINT(34)=MIN(CTPMX0,CTPMX1,CTPMX2,CTPMX3)   
        IF(VINT(33).LE.VINT(13).AND.VINT(34).LE.VINT(14)) MINT(51)=1    
      ELSEIF(ILIM.EQ.4) THEN    
C...Calculate limits on tau'    
C...0) due to kinematics    
        TAPMN0=TAU  
        TAPMX0=1.   
C...1) due to explicit limits   
        TAPMN1=CKIN(31)**2/VINT(2)  
        TAPMX1=1.   
        IF(CKIN(32).GE.0.) TAPMX1=CKIN(32)**2/VINT(2)   
        VINT(16)=MAX(TAPMN0,TAPMN1) 
        VINT(36)=MIN(TAPMX0,TAPMX1) 
        IF(MINT(43).EQ.1) THEN  
          VINT(16)=0.99999  
          VINT(36)=1.00001  
        ENDIF   
        IF(VINT(36).LE.VINT(16)) MINT(51)=1 
      ENDIF 
      RETURN    
C...Special case for low-pT and multiple interactions:  
C...effective kinematical limits for tau, y*, cos(theta-hat).   
  110 IF(ILIM.EQ.0) THEN    
      ELSEIF(ILIM.EQ.1) THEN    
        IF(MSTP(82).LE.1) VINT(11)=4.*PARP(81)**2/VINT(2)   
        IF(MSTP(82).GE.2) VINT(11)=PARP(82)**2/VINT(2)  
        VINT(31)=1. 
      ELSEIF(ILIM.EQ.2) THEN    
        VINT(12)=0.5*LOG(VINT(21))  
        VINT(32)=-VINT(12)  
      ELSEIF(ILIM.EQ.3) THEN    
        IF(MSTP(82).LE.1) ST2EFF=4.*PARP(81)**2/(VINT(21)*VINT(2))  
        IF(MSTP(82).GE.2) ST2EFF=0.01*PARP(82)**2/(VINT(21)*VINT(2))    
        VINT(13)=-SQRT(MAX(0.,1.-ST2EFF))   
        VINT(33)=0. 
        VINT(14)=0. 
        VINT(34)=-VINT(13)  
      ENDIF 
      RETURN    
      END   
C*********************************************************************  
