      SUBROUTINE CRNN(IRUN,PX,PY,PZ,SRT,I1,I2,IBLOCK,
     1NTAG,SIGNN,SIG,NT,ipert1)
*     PURPOSE:                                                         *
*             DEALING WITH NUCLEON-NUCLEON COLLISIONS                    *
*     NOTE   :                                                         *
*     QUANTITIES:                                                 *
*           PX,PY,PZ - MOMENTUM COORDINATES OF ONE PARTICLE IN CM FRAME*
*           SRT      - SQRT OF S                                       *
*           NSTAR =1 INCLUDING N* RESORANCE,ELSE NOT                   *
*           NDIRCT=1 INCLUDING DIRECT PION PRODUCTION PROCESS         *
*           IBLOCK   - THE INFORMATION BACK                            *
*                      0-> COLLISION CANNOT HAPPEN                     *
*                      1-> N-N ELASTIC COLLISION                       *
*                      2-> N+N->N+DELTA,OR N+N->N+N* REACTION          *
*                      3-> N+DELTA->N+N OR N+N*->N+N REACTION          *
*                      4-> N+N->D+D+pion reaction
*                     43->N+N->D(N*)+D(N*) reaction
*                     44->N+N->D+D+rho reaction
*                     45->N+N->N+N+rho
*                     46->N+N->N+N+omega
*           N12       - IS USED TO SPECIFY BARYON-BARYON REACTION      *
*                      CHANNELS. M12 IS THE REVERSAL CHANNEL OF N12    *
*                      N12,                                            *
*                      M12=1 FOR p+n-->delta(+)+ n                     *
*                          2     p+n-->delta(0)+ p                     *
*                          3     p+p-->delta(++)+n                     *
*                          4     p+p-->delta(+)+p                      *
*                          5     n+n-->delta(0)+n                      *
*                          6     n+n-->delta(-)+p                      *
*                          7     n+p-->N*(0)(1440)+p                   *
*                          8     n+p-->N*(+)(1440)+n                   *
*                        9     p+p-->N*(+)(1535)+p                     *
*                        10    n+n-->N*(0)(1535)+n                     *
*                         11    n+p-->N*(+)(1535)+n                     *
*                        12    n+p-->N*(0)(1535)+p
*                        13    D(++)+D(-)-->N*(+)(1440)+n
*                         14    D(++)+D(-)-->N*(0)(1440)+p
*                        15    D(+)+D(0)--->N*(+)(1440)+n
*                        16    D(+)+D(0)--->N*(0)(1440)+p
*                        17    D(++)+D(0)-->N*(+)(1535)+p
*                        18    D(++)+D(-)-->N*(0)(1535)+p
*                        19    D(++)+D(-)-->N*(+)(1535)+n
*                        20    D(+)+D(+)-->N*(+)(1535)+p
*                        21    D(+)+D(0)-->N*(+)(1535)+n
*                        22    D(+)+D(0)-->N*(0)(1535)+p
*                        23    D(+)+D(-)-->N*(0)(1535)+n
*                        24    D(0)+D(0)-->N*(0)(1535)+n
*                          25    N*(+)(14)+N*(+)(14)-->N*(+)(15)+p
*                          26    N*(0)(14)+N*(0)(14)-->N*(0)(15)+n
*                          27    N*(+)(14)+N*(0)(14)-->N*(+)(15)+n
*                        28    N*(+)(14)+N*(0)(14)-->N*(0)(15)+p
*                        29    N*(+)(14)+D+-->N*(+)(15)+p
*                        30    N*(+)(14)+D0-->N*(+)(15)+n
*                        31    N*(+)(14)+D(-)-->N*(0)(1535)+n
*                        32    N*(0)(14)+D++--->N*(+)(15)+p
*                        33    N*(0)(14)+D+--->N*(+)(15)+n
*                        34    N*(0)(14)+D+--->N*(0)(15)+p
*                        35    N*(0)(14)+D0-->N*(0)(15)+n
*                        36    N*(+)(14)+D0--->N*(0)(15)+p
*                        ++    see the note book for more listing
*                     
*
*     NOTE ABOUT N*(1440) RESORANCE IN Nucleon+NUCLEON COLLISION:      * 
*     As it has been discussed in VerWest's paper,I= 1(initial isospin)*
*     channel can all be attributed to delta resorance while I= 0      *
*     channel can all be  attribured to N* resorance.Only in n+p       *
*     one can have I=0 channel so is the N*(1440) resonance            *
*                                                                      *
*                             REFERENCES:                            *    
*                    J. CUGNON ET AL., NUCL. PHYS. A352, 505 (1981)    *
*                    Y. KITAZOE ET AL., PHYS. LETT. 166B, 35 (1986)    *
*                    B. VerWest el al., PHYS. PRV. C25 (1982)1979      *
*                    Gy. Wolf  et al, Nucl Phys A517 (1990) 615;       *
*                                     Nucl phys A552 (1993) 349.       *
**********************************
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,aka=0.498,AP2=0.13957,AM0=1.232,
     2  PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383,APHI=1.020)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        parameter (xmd=1.8756,npdmax=10000)
        COMMON /AA/ R(3,MAXSTR)
cc      SAVE /AA/
        COMMON /BB/ P(3,MAXSTR)
cc      SAVE /BB/
        COMMON /CC/ E(MAXSTR)
cc      SAVE /CC/
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
cc      SAVE /EE/
        common /ff/f(-mx:mx,-my:my,-mz:mz,-mpx:mpx,-mpy:mpy,-mpz:mpzp)
cc      SAVE /ff/
        common /gg/ dx,dy,dz,dpx,dpy,dpz
cc      SAVE /gg/
        COMMON /INPUT/ NSTAR,NDIRCT,DIR
cc      SAVE /INPUT/
        COMMON /NN/NNN
cc      SAVE /NN/
        COMMON /BG/BETAX,BETAY,BETAZ,GAMMA
cc      SAVE /BG/
        COMMON   /RUN/NUM
cc      SAVE /RUN/
        COMMON   /PA/RPION(3,MAXSTR,MAXR)
cc      SAVE /PA/
        COMMON   /PB/PPION(3,MAXSTR,MAXR)
cc      SAVE /PB/
        COMMON   /PC/EPION(MAXSTR,MAXR)
cc      SAVE /PC/
        COMMON   /PD/LPION(MAXSTR,MAXR)
cc      SAVE /PD/
        COMMON/TABLE/ xarray(0:1000),earray(0:1000)
cc      SAVE /TABLE/
        common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
cc      SAVE /input1/
      common/leadng/lb1,px1,py1,pz1,em1,e1,xfnl,yfnl,zfnl,tfnl,
     1 px1n,py1n,pz1n,dp1n
cc      SAVE /leadng/
      COMMON/RNDF77/NSEED
cc      SAVE /RNDF77/
      common /dpi/em2,lb2
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      common /para8/ idpert,npertd,idxsec
      dimension ppd(3,npdmax),lbpd(npdmax)
      SAVE   
*-----------------------------------------------------------------------
      n12=0
      m12=0
      IBLOCK=0
      NTAG=0
      EM1=E(I1)
      EM2=E(I2)
      PR=SQRT( PX**2 + PY**2 + PZ**2 )
      C2=PZ / PR
      X1=RANART(NSEED)
      ianti=0
      if(lb(i1).lt.0 .and. lb(i2).lt.0) ianti=1
      call sbbdm(srt,sdprod,ianti,lbm,xmm,pfinal)
clin-5/2008 Production of perturbative deuterons for idpert=1:
      if(idpert.eq.1.and.ipert1.eq.1) then
         IF (SRT .LT. 2.012) RETURN
         if((iabs(lb(i1)).eq.1.or.iabs(lb(i1)).eq.2)
     1        .and.(iabs(lb(i2)).eq.1.or.iabs(lb(i2)).eq.2)) then
            goto 108
         else
            return
         endif
      endif
c
*-----------------------------------------------------------------------
*COM: TEST FOR ELASTIC SCATTERING (EITHER N-N OR DELTA-DELTA 0R
*      N-DELTA OR N*-N* or N*-Delta)
c      IF (X1 .LE. SIGNN/SIG) THEN
      IF (X1.LE.(SIGNN/SIG)) THEN
*COM:  PARAMETRISATION IS TAKEN FROM THE CUGNON-PAPER
         AS  = ( 3.65 * (SRT - 1.8766) )**6
         A   = 6.0 * AS / (1.0 + AS)
         TA  = -2.0 * PR**2
         X   = RANART(NSEED)
clin-10/24/02        T1  = DLOG( (1-X) * DEXP(dble(A)*dble(TA)) + X )  /  A
         T1  = sngl(DLOG(dble(1.-X)*DEXP(dble(A)*dble(TA))+dble(X)))/  A
         C1  = 1.0 - T1/TA
         T1  = 2.0 * PI * RANART(NSEED)
         IBLOCK=1
         GO TO 107
      ELSE
*COM: TEST FOR INELASTIC SCATTERING
*     IF THE AVAILABLE ENERGY IS LESS THAN THE PION-MASS, NOTHING
*     CAN HAPPEN ANY MORE ==> RETURN (2.012 = 2*AVMASS + PI-MASS)
clin-5/2008: Mdeuteron+Mpi=2.0106 to 2.0152 GeV/c2, so we can still use this:
         IF (SRT .LT. 2.012) RETURN
*     calculate the N*(1535) production cross section in N+N collisions
*     note that the cross sections in this subroutine are in units of mb
*     as only ratios of the cross sections are used to determine the
*     reaction channels
       call N1535(iabs(lb(i1)),iabs(lb(i2)),srt,x1535)
*COM: HERE WE HAVE A PROCESS N+N ==> N+DELTA,OR N+N==>N+N*(144) or N*(1535)
*     OR 
* 3 pi channel : N+N==>d1+d2+PION
       SIG3=3.*(X3pi(SRT)+x33pi(srt))
* 2 pi channel : N+N==>d1+d2+d1*n*+n*n*
       SIG4=4.*X2pi(srt)
* 4 pi channel : N+N==>d1+d2+rho
       s4pi=x4pi(srt)
* N+N-->NN+rho channel
       srho=xrho(srt)
* N+N-->NN+omega
       somega=omega(srt)
* CROSS SECTION FOR KAON PRODUCTION from the four channels
* for NLK channel
       akp=0.498
       ak0=0.498
       ana=0.94
       ada=1.232
       al=1.1157
       as=1.1197
       xsk1=0
       xsk2=0
       xsk3=0
       xsk4=0
       xsk5=0
       t1nlk=ana+al+akp
       if(srt.le.t1nlk)go to 222
       XSK1=1.5*PPLPK(SRT)
* for DLK channel
       t1dlk=ada+al+akp
       t2dlk=ada+al-akp
       if(srt.le.t1dlk)go to 222
       es=srt
       pmdlk2=(es**2-t1dlk**2)*(es**2-t2dlk**2)/(4.*es**2)
       pmdlk=sqrt(pmdlk2)
       XSK3=1.5*PPLPK(srt)
* for NSK channel
       t1nsk=ana+as+akp
       t2nsk=ana+as-akp
       if(srt.le.t1nsk)go to 222
       pmnsk2=(es**2-t1nsk**2)*(es**2-t2nsk**2)/(4.*es**2)
       pmnsk=sqrt(pmnsk2)
       XSK2=1.5*(PPK1(srt)+PPK0(srt))
* for DSK channel
       t1DSk=aDa+aS+akp
       t2DSk=aDa+aS-akp
       if(srt.le.t1dsk)go to 222
       pmDSk2=(es**2-t1DSk**2)*(es**2-t2DSk**2)/(4.*es**2)
       pmDSk=sqrt(pmDSk2)
       XSK4=1.5*(PPK1(srt)+PPK0(srt))
csp11/21/01
c phi production
       if(srt.le.(2.*amn+aphi))go to 222
c  !! mb put the correct form
       xsk5 = 0.0001
csp11/21/01 end
c
* THE TOTAL KAON+ PRODUCTION CROSS SECTION IS THEN
 222   SIGK=XSK1+XSK2+XSK3+XSK4
cbz3/7/99 neutralk
        XSK1 = 2.0 * XSK1
        XSK2 = 2.0 * XSK2
        XSK3 = 2.0 * XSK3
        XSK4 = 2.0 * XSK4
        SIGK = 2.0 * SIGK + xsk5
cbz3/7/99 neutralk end
c
** FOR P+P or L/S+L/S COLLISION:
c       lb1=lb(i1)
c       lb2=lb(i2)
        lb1=iabs(lb(i1))
        lb2=iabs(lb(i2))
        IF((LB(I1)*LB(I2).EQ.1).or.
     &       ((lb1.le.17.and.lb1.ge.14).and.(lb2.le.17.and.lb2.ge.14)).
     &       or.((lb1.le.2).and.(lb2.le.17.and.lb2.ge.14)).
     &       or.((lb2.le.2).and.(lb1.le.17.and.lb1.ge.14)))THEN
clin-8/2008 PP->d+meson here:
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
           SIG1=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
           SIG2=1.5*SIGMA(SRT,1,1,1)
           SIGND=SIG1+SIG2+SIG3+SIG4+X1535+SIGK+s4pi+srho+somega
clin-5/2008:
c           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
           IF (X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
           DIR=SIG3/SIGND
           IF(RANART(NSEED).LE.DIR)GO TO 106
           IF(RANART(NSEED).LE.SIGK/(SIGK+X1535+SIG4+SIG2+SIG1
     &          +s4pi+srho+somega))GO TO 306
           if(RANART(NSEED).le.s4pi/(x1535+sig4+sig2+sig1
     &          +s4pi+srho+somega))go to 307
           if(RANART(NSEED).le.srho/(x1535+sig4+sig2+sig1
     &          +srho+somega))go to 308
           if(RANART(NSEED).le.somega/(x1535+sig4+sig2+sig1
     &          +somega))go to 309
           if(RANART(NSEED).le.x1535/(sig1+sig2+sig4+x1535))then
* N*(1535) production
              N12=9
           ELSE 
              IF(RANART(NSEED).LE.SIG4/(SIG1+sig2+sig4))THEN
* DOUBLE DELTA PRODUCTION
                 N12=66
                 GO TO 1012
              else
*DELTA PRODUCTION
                 N12=3
                 IF (RANART(NSEED).GT.SIG1/(SIG1+SIG2))N12=4
              ENDIF
           endif
           GO TO 1011
        ENDIF
** FOR N+N COLLISION:
        IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
clin-8/2008 NN->d+meson here:
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
           SIG1=SIGMA(SRT,1,1,0)+0.5*SIGMA(SRT,1,1,1)
           SIG2=1.5*SIGMA(SRT,1,1,1)
           SIGND=SIG1+SIG2+X1535+SIG3+SIG4+SIGK+s4pi+srho+somega
clin-5/2008:
c           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
           IF (X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
           dir=sig3/signd
           IF(RANART(NSEED).LE.DIR)GO TO 106
           IF(RANART(NSEED).LE.SIGK/(SIGK+X1535+SIG4+SIG2+SIG1
     &          +s4pi+srho+somega))GO TO 306
           if(RANART(NSEED).le.s4pi/(x1535+sig4+sig2+sig1
     &          +s4pi+srho+somega))go to 307
           if(RANART(NSEED).le.srho/(x1535+sig4+sig2+sig1
     &          +srho+somega))go to 308
           if(RANART(NSEED).le.somega/(x1535+sig4+sig2+sig1
     &          +somega))go to 309
           IF(RANART(NSEED).LE.X1535/(x1535+sig1+sig2+sig4))THEN
* N*(1535) PRODUCTION
              N12=10
           ELSE 
              if(RANART(NSEED).le.sig4/(sig1+sig2+sig4))then
* double delta production
                 N12=67
                 GO TO 1013
              else
* DELTA PRODUCTION
                 N12=6
                 IF (RANART(NSEED).GT.SIG1/(SIG1+SIG2))N12=5
              ENDIF
           endif
           GO TO 1011
        ENDIF
** FOR N+P COLLISION
        IF(LB(I1)*LB(I2).EQ.2)THEN
clin-5/2008 NP->d+meson here:
           IF(X1.LE.((SIGNN+sdprod)/SIG)) GO TO 108
           SIG1=0.5*SIGMA(SRT,1,1,1)+0.25*SIGMA(SRT,1,1,0)
           IF(NSTAR.EQ.1)THEN
              SIG2=(3./4.)*SIGMA(SRT,2,0,1)
           ELSE
              SIG2=0.
           ENDIF
           SIGND=2.*(SIG1+SIG2+X1535)+sig3+sig4+SIGK+s4pi+srho+somega
clin-5/2008:
c           IF (X1.GT.(SIGNN+SIGND)/SIG)RETURN
           IF (X1.GT.(SIGNN+SIGND+sdprod)/SIG)RETURN
           dir=sig3/signd
           IF(RANART(NSEED).LE.DIR)GO TO 106
           IF(RANART(NSEED).LE.SIGK/(SIGND-SIG3))GO TO 306
           if(RANART(NSEED).le.s4pi/(signd-sig3-sigk))go to 307
           if(RANART(NSEED).le.srho/(signd-sig3-sigk-s4pi))go to 308
           if(RANART(NSEED).le.somega/(signd-sig3-sigk-s4pi-srho))
     1          go to 309
           IF(RANART(NSEED).LT.X1535/(SIG1+SIG2+X1535+0.5*sig4))THEN
* N*(1535) PRODUCTION
              N12=11
              IF(RANART(NSEED).LE.0.5)N12=12
           ELSE 
              if(RANART(NSEED).le.sig4/(sig4+2.*(sig1+sig2)))then
* double resonance production
                 N12=68
                 GO TO 1014
              else
                 IF(RANART(NSEED).LE.SIG1/(SIG1+SIG2))THEN
* DELTA PRODUCTION
                    N12=2
                    IF(RANART(NSEED).GE.0.5)N12=1
                 ELSE
* N*(1440) PRODUCTION
                    N12=8
                    IF(RANART(NSEED).GE.0.5)N12=7
                 ENDIF
              ENDIF
           ENDIF
        endif
 1011   iblock=2
        CONTINUE
*PARAMETRIZATION OF THE SHAPE OF THE DELTA RESONANCE ACCORDING
*     TO kitazoe's or J.D.JACKSON'S MASS FORMULA AND BREIT WIGNER
*     FORMULA FOR N* RESORANCE
*     DETERMINE DELTA MASS VIA REJECTION METHOD.
          DMAX = SRT - AVMASS-0.005
          DMAX = SRT - AVMASS-0.005
          DMIN = 1.078
                   IF(N12.LT.7)THEN
* Delta(1232) production
          IF(DMAX.LT.1.232) THEN
          FM=FDE(DMAX,SRT,0.)
          ELSE
clin-10/25/02 get rid of argument usage mismatch in FDE():
             xdmass=1.232
c          FM=FDE(1.232,SRT,1.)
          FM=FDE(xdmass,SRT,1.)
clin-10/25/02-end
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY1=0
10        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FDE(DM,SRT,1.)/FM).AND.
     1    (NTRY1.LE.30)) GOTO 10
clin-2/26/03 limit the Delta mass below a certain value 
c     (here taken as its central value + 2* B-W fullwidth):
          if(dm.gt.1.47) goto 10
              GO TO 13
              ENDIF
                   IF((n12.eq.7).or.(n12.eq.8))THEN
* N*(1440) production
          IF(DMAX.LT.1.44) THEN
          FM=FNS(DMAX,SRT,0.)
          ELSE
clin-10/25/02 get rid of argument usage mismatch in FNS():
             xdmass=1.44
c          FM=FNS(1.44,SRT,1.)
          FM=FNS(xdmass,SRT,1.)
clin-10/25/02-end
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY2=0
11        DM=RANART(NSEED)*(DMAX-DMIN)+DMIN
          NTRY2=NTRY2+1
          IF((RANART(NSEED).GT.FNS(DM,SRT,1.)/FM).AND.
     1    (NTRY2.LE.10)) GO TO 11
clin-2/26/03 limit the N* mass below a certain value 
c     (here taken as its central value + 2* B-W fullwidth):
          if(dm.gt.2.14) goto 11
              GO TO 13
              ENDIF
                    IF(n12.ge.17)then
* N*(1535) production
          IF(DMAX.LT.1.535) THEN
          FM=FD5(DMAX,SRT,0.)
          ELSE
clin-10/25/02 get rid of argument usage mismatch in FNS():
             xdmass=1.535
c          FM=FD5(1.535,SRT,1.)
          FM=FD5(xdmass,SRT,1.)
clin-10/25/02-end
          ENDIF
          IF(FM.EQ.0.)FM=1.E-09
          NTRY1=0
12        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FD5(DM,SRT,1.)/FM).AND.
     1    (NTRY1.LE.10)) GOTO 12
clin-2/26/03 limit the N* mass below a certain value 
c     (here taken as its central value + 2* B-W fullwidth):
          if(dm.gt.1.84) goto 12
         GO TO 13
             ENDIF
* CALCULATE THE MASSES OF BARYON RESONANCES IN THE DOUBLE RESONANCE
* PRODUCTION PROCESS AND RELABLE THE PARTICLES
1012       iblock=43
       call Rmasdd(srt,1.232,1.232,1.08,
     &  1.08,ISEED,1,dm1,dm2)
       call Rmasdd(srt,1.232,1.44,1.08,
     &  1.08,ISEED,3,dm1n,dm2n)
       IF(N12.EQ.66)THEN
*(1) PP-->DOUBLE RESONANCES
* DETERMINE THE FINAL STATE
       XFINAL=RANART(NSEED)
       IF(XFINAL.LE.0.25)THEN
* (1.1) D+++D0 
       LB(I1)=9
       LB(I2)=7
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF((XFINAL.gt.0.25).and.(xfinal.le.0.5))THEN
* (1.2) D++D+
       LB(I1)=8
       LB(I2)=8
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF((XFINAL.gt.0.5).and.(xfinal.le.0.75))THEN
* (1.3) D+++N*0 
       LB(I1)=9
       LB(I2)=10
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF(XFINAL.gt.0.75)then
* (1.4) D++N*+ 
       LB(I1)=8
       LB(I2)=11
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       ENDIF
1013       iblock=43
       call Rmasdd(srt,1.232,1.232,1.08,
     &  1.08,ISEED,1,dm1,dm2)
       call Rmasdd(srt,1.232,1.44,1.08,
     &  1.08,ISEED,3,dm1n,dm2n)
       IF(N12.EQ.67)THEN
*(2) NN-->DOUBLE RESONANCES
* DETERMINE THE FINAL STATE
       XFINAL=RANART(NSEED)
       IF(XFINAL.LE.0.25)THEN
* (2.1) D0+D0 
       LB(I1)=7
       LB(I2)=7
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
* go to 200 to set the new momentum
        ENDIF
       IF((XFINAL.gt.0.25).and.(xfinal.le.0.5))THEN
* (2.2) D++D+
       LB(I1)=6
       LB(I2)=8
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF((XFINAL.gt.0.5).and.(xfinal.le.0.75))THEN
* (2.3) D0+N*0 
       LB(I1)=7
       LB(I2)=10
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF(XFINAL.gt.0.75)then
* (2.4) D++N*+ 
       LB(I1)=8
       LB(I2)=11
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       ENDIF
1014       iblock=43
       call Rmasdd(srt,1.232,1.232,1.08,
     &  1.08,ISEED,1,dm1,dm2)
       call Rmasdd(srt,1.232,1.44,1.08,
     &  1.08,ISEED,3,dm1n,dm2n)
       IF(N12.EQ.68)THEN
*(3) NP-->DOUBLE RESONANCES
* DETERMINE THE FINAL STATE
       XFINAL=RANART(NSEED)
       IF(XFINAL.LE.0.25)THEN
* (3.1) D0+D+ 
       LB(I1)=7
       LB(I2)=8
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF((XFINAL.gt.0.25).and.(xfinal.le.0.5))THEN
* (3.2) D+++D-
       LB(I1)=9
       LB(I2)=6
       e(i1)=dm1
       e(i2)=dm2
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF((XFINAL.gt.0.5).and.(xfinal.le.0.75))THEN
* (3.3) D0+N*+ 
       LB(I1)=7
       LB(I2)=11
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       IF(XFINAL.gt.0.75)then
* (3.4) D++N*0
       LB(I1)=8
       LB(I2)=10
       e(i1)=dm1n
       e(i2)=dm2n
       GO TO 200
* go to 200 to set the new momentum
       ENDIF
       ENDIF
13       CONTINUE
*-------------------------------------------------------
* RELABLE BARYON I1 AND I2
*1. p+n-->delta(+)+n
          IF(N12.EQ.1)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I2)=2
          LB(I1)=8
          E(I1)=DM
          ELSE
          LB(I1)=2
          LB(I2)=8
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
*2 p+n-->delta(0)+p
          IF(N12.EQ.2)THEN
          IF(iabs(LB(I1)).EQ.2)THEN
          LB(I2)=1
          LB(I1)=7
          E(I1)=DM
          ELSE
          LB(I1)=1
          LB(I2)=7
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
*3 p+p-->delta(++)+n
          IF(N12.EQ.3)THEN
          LB(I1)=9
          E(I1)=DM
          LB(I2)=2
          E(I2)=AMN
         GO TO 200
          ENDIF
*4 p+p-->delta(+)+p
          IF(N12.EQ.4)THEN
          LB(I2)=1
          LB(I1)=8
          E(I1)=DM
         GO TO 200
          ENDIF
*5 n+n--> delta(0)+n
          IF(N12.EQ.5)THEN
          LB(I2)=2
          LB(I1)=7
          E(I1)=DM
         GO TO 200
          ENDIF
*6 n+n--> delta(-)+p
          IF(N12.EQ.6)THEN
          LB(I1)=6
          E(I1)=DM
          LB(I2)=1
          E(I2)=AMP
         GO TO 200
          ENDIF
*7 n+p--> N*(0)+p
          IF(N12.EQ.7)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I1)=1
          LB(I2)=10
          E(I2)=DM
          ELSE
          LB(I2)=1
          LB(I1)=10
          E(I1)=DM
          ENDIF
         GO TO 200
          ENDIF
*8 n+p--> N*(+)+n
          IF(N12.EQ.8)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I2)=2
          LB(I1)=11
          E(I1)=DM
          ELSE
          LB(I1)=2
          LB(I2)=11
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
*9 p+p--> N*(+)(1535)+p
          IF(N12.EQ.9)THEN
          IF(RANART(NSEED).le.0.5)THEN
          LB(I2)=1
          LB(I1)=13
          E(I1)=DM
          ELSE
          LB(I1)=1
          LB(I2)=13
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
*10 n+n--> N*(0)(1535)+n
          IF(N12.EQ.10)THEN
          IF(RANART(NSEED).le.0.5)THEN
          LB(I2)=2
          LB(I1)=12
          E(I1)=DM
          ELSE
          LB(I1)=2
          LB(I2)=12
          E(I2)=DM
          ENDIF
         GO TO 200
          ENDIF
*11 n+p--> N*(+)(1535)+n
          IF(N12.EQ.11)THEN
          IF(iabs(LB(I1)).EQ.2)THEN
          LB(I1)=2
          LB(I2)=13
          E(I2)=DM
          ELSE
          LB(I2)=2
          LB(I1)=13
          E(I1)=DM
          ENDIF
         GO TO 200
          ENDIF
*12 n+p--> N*(0)(1535)+p
          IF(N12.EQ.12)THEN
          IF(iabs(LB(I1)).EQ.1)THEN
          LB(I1)=1
          LB(I2)=12
          E(I2)=DM
          ELSE
          LB(I2)=1
          LB(I1)=12
          E(I1)=DM
          ENDIF
          ENDIF
         endif
* CALCULATE THE MAGNITUDE OF THE FINAL MOMENTUM THROUGH
* ENERGY CONSERVATION
200       EM1=E(I1)
          EM2=E(I2)
          PR2   = (SRT**2 - EM1**2 - EM2**2)**2
     1                - 4.0 * (EM1*EM2)**2
          IF(PR2.LE.0.)PR2=1.e-09
          PR=SQRT(PR2)/(2.*SRT)
              if(srt.le.2.14)C1= 1.0 - 2.0 * RANART(NSEED)
         if(srt.gt.2.14.and.srt.le.2.4)c1=ang(srt,iseed)
         if(srt.gt.2.4)then
clin-10/25/02 get rid of argument usage mismatch in PTR():
             xptr=0.33*pr
c         cc1=ptr(0.33*pr,iseed)
             cc1=ptr(xptr,iseed)
clin-10/25/02-end
clin-9/2012: check argument in sqrt():
             scheck=pr**2-cc1**2
             if(scheck.lt.0) then
                write(99,*) 'scheck2: ', scheck
                scheck=0.
             endif
             c1=sqrt(scheck)/pr
c             c1=sqrt(pr**2-cc1**2)/pr
         endif
          T1   = 2.0 * PI * RANART(NSEED)
       if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
         lb(i1) = -lb(i1)
         lb(i2) = -lb(i2)
       endif
          GO TO 107
*FOR THE NN-->D1+D2+PI PROCESS, FIND MOMENTUM OF THE FINAL TWO
*DELTAS AND PION IN THE NUCLEUS-NUCLEUS CMS.
106     CONTINUE
           NTRY1=0
123        CALL DDP2(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.40))GO TO 123
C       if(icou1.lt.0)return
* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
* (1) FOR P+P
              XDIR=RANART(NSEED)
                IF(LB(I1)*LB(I2).EQ.1)THEN
                IF(XDIR.Le.0.2)then
* (1.1)P+P-->D+++D0+PION(0)
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
              LB(I1)=9
              LB(I2)=7
       GO TO 205
                ENDIF
* (1.2)P+P -->D++D+PION(0)
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
                LB(I1)=8
                LB(I2)=8
       GO TO 205
              ENDIF 
* (1.3)P+P-->D+++D+PION(-)
                IF((XDIR.LE.0.6).AND.(XDIR.GT.0.4))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=9
                LB(I2)=8
       GO TO 205
              ENDIF 
                IF((XDIR.LE.0.8).AND.(XDIR.GT.0.6))THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=9
                LB(I2)=6
       GO TO 205
              ENDIF 
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
               ENDIF
* (2)FOR N+N
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
                IF(XDIR.Le.0.2)then
* (2.1)N+N-->D++D-+PION(0)
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
              LB(I1)=6
              LB(I2)=7
       GO TO 205
                ENDIF
* (2.2)N+N -->D+++D-+PION(-)
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=6
                LB(I2)=9
       GO TO 205
              ENDIF 
* (2.3)P+P-->D0+D-+PION(+)
                IF((XDIR.GT.0.4).AND.(XDIR.LE.0.6))THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=9
                LB(I2)=8
       GO TO 205
              ENDIF 
* (2.4)P+P-->D0+D0+PION(0)
                IF((XDIR.GT.0.6).AND.(XDIR.LE.0.8))THEN
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
                LB(I1)=7
                LB(I2)=7
       GO TO 205
              ENDIF 
* (2.5)P+P-->D0+D++PION(-)
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
              ENDIF
* (3)FOR N+P
                IF(LB(I1)*LB(I2).EQ.2)THEN
                IF(XDIR.Le.0.17)then
* (3.1)N+P-->D+++D-+PION(0)
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP1
              LB(I1)=6
              LB(I2)=9
       GO TO 205
                ENDIF
* (3.2)N+P -->D+++D0+PION(-)
                IF((XDIR.LE.0.34).AND.(XDIR.GT.0.17))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=9
       GO TO 205
              ENDIF 
* (3.3)N+P-->D++D-+PION(+)
                IF((XDIR.GT.0.34).AND.(XDIR.LE.0.51))THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
* (3.4)N+P-->D++D++PION(-)
                IF((XDIR.GT.0.51).AND.(XDIR.LE.0.68))THEN
                LPION(NNN,IRUN)=3
                EPION(NNN,IRUN)=AP2
                LB(I1)=8
                LB(I2)=8
       GO TO 205
              ENDIF 
* (3.5)N+P-->D0+D++PION(0)
                IF((XDIR.GT.0.68).AND.(XDIR.LE.0.85))THEN
                LPION(NNN,IRUN)=4
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=8
       GO TO 205
              ENDIF 
* (3.6)N+P-->D0+D0+PION(+)
                IF(XDIR.GT.0.85)THEN
                LPION(NNN,IRUN)=5
                EPION(NNN,IRUN)=AP2
                LB(I1)=7
                LB(I2)=7
              ENDIF 
                ENDIF
* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
* NUCLEUS CMS. FRAME 
*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
205           E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
c
             if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
                if(LPION(NNN,IRUN) .eq. 3)then
                  LPION(NNN,IRUN)=5
                elseif(LPION(NNN,IRUN) .eq. 5)then
                  LPION(NNN,IRUN)=3
                endif
               endif
c
             lb1=lb(i1)
* FOR DELTA2
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lb2=lb(i2)
* assign delta1 and delta2 to i1 or i2 to keep the leadng particle
* behaviour
C              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=4
* GET PION'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
clin-5/2008:
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
clin-5/2008 do not allow smearing in position of produced particles 
c     to avoid immediate reinteraction with the particle I1, I2 or themselves:
c2002        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2002
c                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
c                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
c                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
                RPION(1,NNN,IRUN)=R(1,I1)
                RPION(2,NNN,IRUN)=R(2,I1)
                RPION(3,NNN,IRUN)=R(3,I1)
c
              go to 90005
clin-5/2008 N+N->Deuteron+pi:
*     FIND MOMENTUM OF THE FINAL PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
 108       CONTINUE
           if(idpert.eq.1.and.ipert1.eq.1.and.npertd.ge.1) then
c     For idpert=1: we produce npertd pert deuterons:
              ndloop=npertd
           elseif(idpert.eq.2.and.npertd.ge.1) then
c     For idpert=2: we first save information for npertd pert deuterons;
c     at the last ndloop we create the regular deuteron+pi 
c     and those pert deuterons:
              ndloop=npertd+1
           else
c     Just create the regular deuteron+pi:
              ndloop=1
           endif
c
           dprob1=sdprod/sig/float(npertd)
           do idloop=1,ndloop
              CALL bbdangle(pxd,pyd,pzd,nt,ipert1,ianti,idloop,pfinal,
     1 dprob1,lbm)
              CALL ROTATE(PX,PY,PZ,PXd,PYd,PZd)
*     LORENTZ-TRANSFORMATION OF THE MOMENTUM OF PARTICLES IN THE FINAL STATE 
*     FROM THE NN CMS FRAME INTO THE GLOBAL CMS FRAME:
*     For the Deuteron:
              xmass=xmd
              E1dCM=SQRT(xmass**2+PXd**2+PYd**2+PZd**2)
              P1dBETA=PXd*BETAX+PYd*BETAY+PZd*BETAZ
              TRANSF=GAMMA*(GAMMA*P1dBETA/(GAMMA+1.)+E1dCM)
              pxi1=BETAX*TRANSF+PXd
              pyi1=BETAY*TRANSF+PYd
              pzi1=BETAZ*TRANSF+PZd
              if(ianti.eq.0)then
                 lbd=42
              else
                 lbd=-42
              endif
              if(idpert.eq.1.and.ipert1.eq.1.and.npertd.ge.1) then
cccc  Perturbative production for idpert=1:
                 nnn=nnn+1
                 PPION(1,NNN,IRUN)=pxi1
                 PPION(2,NNN,IRUN)=pyi1
                 PPION(3,NNN,IRUN)=pzi1
                 EPION(NNN,IRUN)=xmd
                 LPION(NNN,IRUN)=lbd
                 RPION(1,NNN,IRUN)=R(1,I1)
                 RPION(2,NNN,IRUN)=R(2,I1)
                 RPION(3,NNN,IRUN)=R(3,I1)
clin-5/2008 assign the perturbative probability:
                 dppion(NNN,IRUN)=sdprod/sig/float(npertd)
              elseif(idpert.eq.2.and.idloop.le.npertd) then
clin-5/2008 For idpert=2, we produce NPERTD perturbative (anti)deuterons 
c     only when a regular (anti)deuteron+pi is produced in NN collisions.
c     First save the info for the perturbative deuterons:
                 ppd(1,idloop)=pxi1
                 ppd(2,idloop)=pyi1
                 ppd(3,idloop)=pzi1
                 lbpd(idloop)=lbd
              else
cccc  Regular production:
c     For the regular pion: do LORENTZ-TRANSFORMATION:
                 E(i1)=xmm
                 E2piCM=SQRT(xmm**2+PXd**2+PYd**2+PZd**2)
                 P2piBETA=-PXd*BETAX-PYd*BETAY-PZd*BETAZ
                 TRANSF=GAMMA*(GAMMA*P2piBETA/(GAMMA+1.)+E2piCM)
                 pxi2=BETAX*TRANSF-PXd
                 pyi2=BETAY*TRANSF-PYd
                 pzi2=BETAZ*TRANSF-PZd
                 p(1,i1)=pxi2
                 p(2,i1)=pyi2
                 p(3,i1)=pzi2
c     Remove regular pion to check the equivalence 
c     between the perturbative and regular deuteron results:
c                 E(i1)=0.
c
                 LB(I1)=lbm
                 PX1=P(1,I1)
                 PY1=P(2,I1)
                 PZ1=P(3,I1)
                 EM1=E(I1)
                 ID(I1)=2
                 ID1=ID(I1)
                 E1=SQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
                 lb1=lb(i1)
c     For the regular deuteron:
                 p(1,i2)=pxi1
                 p(2,i2)=pyi1
                 p(3,i2)=pzi1
                 lb(i2)=lbd
                 lb2=lb(i2)
                 E(i2)=xmd
                 EtI2=E(I2)
                 ID(I2)=2
c     For idpert=2: create the perturbative deuterons:
                 if(idpert.eq.2.and.idloop.eq.ndloop) then
                    do ipertd=1,npertd
                       nnn=nnn+1
                       PPION(1,NNN,IRUN)=ppd(1,ipertd)
                       PPION(2,NNN,IRUN)=ppd(2,ipertd)
                       PPION(3,NNN,IRUN)=ppd(3,ipertd)
                       EPION(NNN,IRUN)=xmd
                       LPION(NNN,IRUN)=lbpd(ipertd)
                       RPION(1,NNN,IRUN)=R(1,I1)
                       RPION(2,NNN,IRUN)=R(2,I1)
                       RPION(3,NNN,IRUN)=R(3,I1)
clin-5/2008 assign the perturbative probability:
                       dppion(NNN,IRUN)=1./float(npertd)
                    enddo
                 endif
              endif
           enddo
           IBLOCK=501
           go to 90005
clin-5/2008 N+N->Deuteron+pi over
* FOR THE NN-->KAON+X PROCESS, FIND MOMENTUM OF THE FINAL PARTICLES IN 
* THE NUCLEUS-NUCLEUS CMS.
306     CONTINUE
csp11/21/01 phi production
              if(XSK5/sigK.gt.RANART(NSEED))then
              pz1=p(3,i1)
              pz2=p(3,i2)
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 1 + int(2 * RANART(NSEED))
              nnn=nnn+1
                LPION(NNN,IRUN)=29
                EPION(NNN,IRUN)=APHI
                iblock = 222
              GO TO 208
               ENDIF
c
                 IBLOCK=9
                 if(ianti .eq. 1)iblock=-9
c
              pz1=p(3,i1)
              pz2=p(3,i2)
* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
              nnn=nnn+1
                LPION(NNN,IRUN)=23
                EPION(NNN,IRUN)=Aka
              if(srt.le.2.63)then
* only lambda production is possible
* (1.1)P+P-->p+L+kaon+
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              GO TO 208
                ENDIF
       if(srt.le.2.74.and.srt.gt.2.63)then
* both Lambda and sigma production are possible
              if(XSK1/(XSK1+XSK2).gt.RANART(NSEED))then
* lambda production
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              else
* sigma production
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              ic=2
              endif
              GO TO 208
       endif
       if(srt.le.2.77.and.srt.gt.2.74)then
* then pp-->Delta lamda kaon can happen
              if(xsk1/(xsk1+xsk2+xsk3).
     1          gt.RANART(NSEED))then
* * (1.1)P+P-->p+L+kaon+
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              go to 208
              else
              if(xsk2/(xsk2+xsk3).gt.RANART(NSEED))then
* pp-->psk
              ic=2
                LB(I1) = 1 + int(2 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              else
* pp-->D+l+k        
              ic=3
                LB(I1) = 6 + int(4 * RANART(NSEED))
              lb(i2)=14
              endif
              GO TO 208
              endif
       endif
       if(srt.gt.2.77)then
* all four channels are possible
              if(xsk1/(xsk1+xsk2+xsk3+xsk4).gt.RANART(NSEED))then
* p lambda k production
              ic=1
                LB(I1) = 1 + int(2 * RANART(NSEED))
              LB(I2)=14
              go to 208
       else
          if(xsk3/(xsk2+xsk3+xsk4).gt.RANART(NSEED))then
* delta l K production
              ic=3
                LB(I1) = 6 + int(4 * RANART(NSEED))
              lb(i2)=14
              go to 208
          else
              if(xsk2/(xsk2+xsk4).gt.RANART(NSEED))then
* n sigma k production
                   LB(I1) = 1 + int(2 * RANART(NSEED))
                   LB(I2) = 15 + int(3 * RANART(NSEED))
              ic=2
              else
              ic=4
                LB(I1) = 6 + int(4 * RANART(NSEED))
                LB(I2) = 15 + int(3 * RANART(NSEED))
              endif
              go to 208
          endif
       endif
       endif
208             continue
         if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
          lb(i1) = - lb(i1)
          lb(i2) = - lb(i2)
          if(LPION(NNN,IRUN) .eq. 23)LPION(NNN,IRUN)=21
         endif
* KEEP ALL COORDINATES OF PARTICLE 2 FOR POSSIBLE PHASE SPACE CHANGE
           NTRY1=0
127        CALL BBKAON(ic,SRT,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 127
c       if(icou1.lt.0)return
* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
* NUCLEUS CMS. FRAME 
* (1) for the necleon/delta
*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
              E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
             lbi1=lb(i1)
* (2) for the lambda/sigma
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lbi2=lb(i2)
* GET the kaon'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
                EPCM=SQRT(aka**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
clin-5/2008
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
clin-5/2008
c2003        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2003
c                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
c                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
c                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
                RPION(1,NNN,IRUN)=R(1,I1)
                RPION(2,NNN,IRUN)=R(2,I1)
                RPION(3,NNN,IRUN)=R(3,I1)
c
* assign the nucleon/delta and lambda/sigma to i1 or i2 to keep the 
* leadng particle behaviour
C              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lbi1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lbi2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
              go to 90005
* FOR THE NN-->Delta+Delta+rho PROCESS, FIND MOMENTUM OF THE FINAL 
* PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
307     CONTINUE
           NTRY1=0
125        CALL DDrho(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,amrho,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 125
C       if(icou1.lt.0)return
* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
              arho=amrho
* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
* (1) FOR P+P
              XDIR=RANART(NSEED)
                IF(LB(I1)*LB(I2).EQ.1)THEN
                IF(XDIR.Le.0.2)then
* (1.1)P+P-->D+++D0+rho(0)
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=9
              LB(I2)=7
       GO TO 2051
                ENDIF
* (1.2)P+P -->D++D+rho(0)
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
                LB(I1)=8
                LB(I2)=8
       GO TO 2051
              ENDIF 
* (1.3)P+P-->D+++D+arho(-)
                IF((XDIR.LE.0.6).AND.(XDIR.GT.0.4))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=9
                LB(I2)=8
       GO TO 2051
              ENDIF 
                IF((XDIR.LE.0.8).AND.(XDIR.GT.0.6))THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=9
                LB(I2)=6
       GO TO 2051
              ENDIF 
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
               ENDIF
* (2)FOR N+N
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
                IF(XDIR.Le.0.2)then
* (2.1)N+N-->D++D-+rho(0)
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=6
              LB(I2)=7
       GO TO 2051
                ENDIF
* (2.2)N+N -->D+++D-+rho(-)
                IF((XDIR.LE.0.4).AND.(XDIR.GT.0.2))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=6
                LB(I2)=9
       GO TO 2051
              ENDIF 
* (2.3)P+P-->D0+D-+rho(+)
                IF((XDIR.GT.0.4).AND.(XDIR.LE.0.6))THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=9
                LB(I2)=8
       GO TO 2051
              ENDIF 
* (2.4)P+P-->D0+D0+rho(0)
                IF((XDIR.GT.0.6).AND.(XDIR.LE.0.8))THEN
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=7
       GO TO 2051
              ENDIF 
* (2.5)P+P-->D0+D++rho(-)
                IF(XDIR.GT.0.8)THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
              ENDIF
* (3)FOR N+P
                IF(LB(I1)*LB(I2).EQ.2)THEN
                IF(XDIR.Le.0.17)then
* (3.1)N+P-->D+++D-+rho(0)
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
              LB(I1)=6
              LB(I2)=9
       GO TO 2051
                ENDIF
* (3.2)N+P -->D+++D0+rho(-)
                IF((XDIR.LE.0.34).AND.(XDIR.GT.0.17))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=9
       GO TO 2051
              ENDIF 
* (3.3)N+P-->D++D-+rho(+)
                IF((XDIR.GT.0.34).AND.(XDIR.LE.0.51))THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
* (3.4)N+P-->D++D++rho(-)
                IF((XDIR.GT.0.51).AND.(XDIR.LE.0.68))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=8
                LB(I2)=8
       GO TO 2051
              ENDIF 
* (3.5)N+P-->D0+D++rho(0)
                IF((XDIR.GT.0.68).AND.(XDIR.LE.0.85))THEN
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=8
       GO TO 2051
              ENDIF 
* (3.6)N+P-->D0+D0+rho(+)
                IF(XDIR.GT.0.85)THEN
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=7
                LB(I2)=7
              ENDIF 
                ENDIF
* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
* NUCLEUS CMS. FRAME 
*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
2051          E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
c
             if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
                if(LPION(NNN,IRUN) .eq. 25)then
                  LPION(NNN,IRUN)=27
                elseif(LPION(NNN,IRUN) .eq. 27)then
                  LPION(NNN,IRUN)=25
                endif
               endif
c
             lb1=lb(i1)
* FOR DELTA2
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lb2=lb(i2)
* assign delta1 and delta2 to i1 or i2 to keep the leadng particle
* behaviour
C              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=44
* GET rho'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
clin-5/2008:
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
clin-5/2008:
c2004        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2004
c                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
c                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
c                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
                RPION(1,NNN,IRUN)=R(1,I1)
                RPION(2,NNN,IRUN)=R(2,I1)
                RPION(3,NNN,IRUN)=R(3,I1)
c
              go to 90005
* FOR THE NN-->N+N+rho PROCESS, FIND MOMENTUM OF THE FINAL 
* PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
308     CONTINUE
           NTRY1=0
126        CALL pprho(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,amrho,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 126
C       if(icou1.lt.0)return
* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
              arho=amrho
* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
* (1) FOR P+P
              XDIR=RANART(NSEED)
                IF(LB(I1)*LB(I2).EQ.1)THEN
                IF(XDIR.Le.0.5)then
* (1.1)P+P-->P+P+rho(0)
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=1
              LB(I2)=1
       GO TO 2052
                Else
* (1.2)P+P -->p+n+rho(+)
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=1
                LB(I2)=2
       GO TO 2052
              ENDIF 
              endif
* (2)FOR N+N
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
                IF(XDIR.Le.0.5)then
* (2.1)N+N-->N+N+rho(0)
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=2
              LB(I2)=2
       GO TO 2052
                Else
* (2.2)N+N -->N+P+rho(-)
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=1
                LB(I2)=2
       GO TO 2052
              ENDIF 
              endif
* (3)FOR N+P
                IF(LB(I1)*LB(I2).EQ.2)THEN
                IF(XDIR.Le.0.33)then
* (3.1)N+P-->N+P+rho(0)
                LPION(NNN,IRUN)=26
                EPION(NNN,IRUN)=Arho
              LB(I1)=1
              LB(I2)=2
       GO TO 2052
* (3.2)N+P -->P+P+rho(-)
                else IF((XDIR.LE.0.67).AND.(XDIR.GT.0.34))THEN
                LPION(NNN,IRUN)=25
                EPION(NNN,IRUN)=Arho
                LB(I1)=1
                LB(I2)=1
       GO TO 2052
              Else 
* (3.3)N+P-->N+N+rho(+)
                LPION(NNN,IRUN)=27
                EPION(NNN,IRUN)=Arho
                LB(I1)=2
                LB(I2)=2
       GO TO 2052
              ENDIF 
              endif
* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
* NUCLEUS CMS. FRAME 
*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
2052          E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
c
              if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
                if(LPION(NNN,IRUN) .eq. 25)then
                  LPION(NNN,IRUN)=27
                elseif(LPION(NNN,IRUN) .eq. 27)then
                  LPION(NNN,IRUN)=25
                endif
               endif
c
             lb1=lb(i1)
* FOR p2
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
              lb2=lb(i2)
* assign p1 and p2 to i1 or i2 to keep the leadng particle
* behaviour
C              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=45
* GET rho'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
clin-5/2008:
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
clin-5/2008:
c2005        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2005
c                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
c                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
c                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
                RPION(1,NNN,IRUN)=R(1,I1)
                RPION(2,NNN,IRUN)=R(2,I1)
                RPION(3,NNN,IRUN)=R(3,I1)
c
              go to 90005
* FOR THE NN-->p+p+omega PROCESS, FIND MOMENTUM OF THE FINAL 
* PARTICLES IN THE NUCLEUS-NUCLEUS CMS.
309     CONTINUE
           NTRY1=0
138        CALL ppomga(SRT,ISEED,PX3,PY3,PZ3,DM3,PX4,PY4,PZ4,DM4,
     &  PPX,PPY,PPZ,icou1)
       NTRY1=NTRY1+1
       if((icou1.lt.0).AND.(NTRY1.LE.20))GO TO 138
C       if(icou1.lt.0)return
* ROTATE THE MOMENTA OF PARTICLES IN THE CMS OF P1+P2
       CALL ROTATE(PX,PY,PZ,PX3,PY3,PZ3)
       CALL ROTATE(PX,PY,PZ,PX4,PY4,PZ4)
       CALL ROTATE(PX,PY,PZ,PPX,PPY,PPZ)
                NNN=NNN+1
              aomega=0.782
* DETERMINE THE CHARGE STATES OF PARTICLES IN THE FINAL STATE
* (1) FOR P+P
                IF(LB(I1)*LB(I2).EQ.1)THEN
* (1.1)P+P-->P+P+omega(0)
                LPION(NNN,IRUN)=28
                EPION(NNN,IRUN)=Aomega
              LB(I1)=1
              LB(I2)=1
       GO TO 2053
                ENDIF
* (2)FOR N+N
                IF(iabs(LB(I1)).EQ.2.AND.iabs(LB(I2)).EQ.2)THEN
* (2.1)N+N-->N+N+omega(0)
                LPION(NNN,IRUN)=28
                EPION(NNN,IRUN)=Aomega
              LB(I1)=2
              LB(I2)=2
       GO TO 2053
                ENDIF
* (3)FOR N+P
                IF(LB(I1)*LB(I2).EQ.2)THEN
* (3.1)N+P-->N+P+omega(0)
                LPION(NNN,IRUN)=28
                EPION(NNN,IRUN)=Aomega
              LB(I1)=1
              LB(I2)=2
       GO TO 2053
                ENDIF
* FIND THE MOMENTUM OF PARTICLES IN THE FINAL STATE IN THE NUCLEUS-
* NUCLEUS CMS. FRAME 
*             LORENTZ-TRANSFORMATION INTO LAB FRAME FOR DELTA1
2053          E1CM    = SQRT (dm3**2 + PX3**2 + PY3**2 + PZ3**2)
              P1BETA  = PX3*BETAX + PY3*BETAY + PZ3*BETAZ
              TRANSF  = GAMMA * ( GAMMA * P1BETA / (GAMMA + 1) + E1CM )
              Pt1i1 = BETAX * TRANSF + PX3
              Pt2i1 = BETAY * TRANSF + PY3
              Pt3i1 = BETAZ * TRANSF + PZ3
             Eti1   = DM3
              if(ianti.eq.1 .and. lb(i1).ge.1 .and. lb(i2).ge.1)then
               lb(i1) = -lb(i1)
               lb(i2) = -lb(i2)
               endif
             lb1=lb(i1)
* FOR DELTA2
                E2CM    = SQRT (dm4**2 + PX4**2 + PY4**2 + PZ4**2)
                P2BETA  = PX4*BETAX+PY4*BETAY+PZ4*BETAZ
                TRANSF  = GAMMA * (GAMMA*P2BETA / (GAMMA + 1.) + E2CM)
                Pt1I2 = BETAX * TRANSF + PX4
                Pt2I2 = BETAY * TRANSF + PY4
                Pt3I2 = BETAZ * TRANSF + PZ4
              EtI2   = DM4
                lb2=lb(i2)
* assign delta1 and delta2 to i1 or i2 to keep the leadng particle
* behaviour
C              if((pt1i1*px1+pt2i1*py1+pt3i1*pz1).gt.0)then
              p(1,i1)=pt1i1
              p(2,i1)=pt2i1
              p(3,i1)=pt3i1
              e(i1)=eti1
              lb(i1)=lb1
              p(1,i2)=pt1i2
              p(2,i2)=pt2i2
              p(3,i2)=pt3i2
              e(i2)=eti2
              lb(i2)=lb2
                PX1     = P(1,I1)
                PY1     = P(2,I1)
                PZ1     = P(3,I1)
              EM1       = E(I1)
                ID(I1)  = 2
                ID(I2)  = 2
                ID1     = ID(I1)
                IBLOCK=46
* GET omega'S MOMENTUM AND COORDINATES IN NUCLEUS-NUCLEUS CMS. FRAME
                EPCM=SQRT(EPION(NNN,IRUN)**2+PPX**2+PPY**2+PPZ**2)
                PPBETA=PPX*BETAX+PPY*BETAY+PPZ*BETAZ
                TRANSF=GAMMA*(GAMMA*PPBETA/(GAMMA+1.)+EPCM)
                PPION(1,NNN,IRUN)=BETAX*TRANSF+PPX
                PPION(2,NNN,IRUN)=BETAY*TRANSF+PPY
                PPION(3,NNN,IRUN)=BETAZ*TRANSF+PPZ
clin-5/2008:
                dppion(nnn,irun)=dpertp(i1)*dpertp(i2)
clin-5/2008:
c2006        X01 = 1.0 - 2.0 * RANART(NSEED)
c            Y01 = 1.0 - 2.0 * RANART(NSEED)
c            Z01 = 1.0 - 2.0 * RANART(NSEED)
c        IF ((X01*X01+Y01*Y01+Z01*Z01) .GT. 1.0) GOTO 2006
c                RPION(1,NNN,IRUN)=R(1,I1)+0.5*x01
c                RPION(2,NNN,IRUN)=R(2,I1)+0.5*y01
c                RPION(3,NNN,IRUN)=R(3,I1)+0.5*z01
                    RPION(1,NNN,IRUN)=R(1,I1)
                    RPION(2,NNN,IRUN)=R(2,I1)
                    RPION(3,NNN,IRUN)=R(3,I1)
c
              go to 90005
* change phase space density FOR NUCLEONS AFTER THE PROCESS
clin-10/25/02-comment out following, since there is no path to it:
clin-8/16/02 used before set
c     IX1,IY1,IZ1,IPX1,IPY1,IPZ1, IX2,IY2,IZ2,IPX2,IPY2,IPZ2:
c                if ((abs(ix1).le.mx) .and. (abs(iy1).le.my) .and.
c     &              (abs(iz1).le.mz)) then
c                  ipx1p = nint(p(1,i1)/dpx)
c                  ipy1p = nint(p(2,i1)/dpy)
c                  ipz1p = nint(p(3,i1)/dpz)
c                  if ((ipx1p.ne.ipx1) .or. (ipy1p.ne.ipy1) .or.
c     &                (ipz1p.ne.ipz1)) then
c                    if ((abs(ipx1).le.mpx) .and. (abs(ipy1).le.my)
c     &                .and. (ipz1.ge.-mpz) .and. (ipz1.le.mpzp))
c     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) =
c     &                f(ix1,iy1,iz1,ipx1,ipy1,ipz1) - 1.
c                    if ((abs(ipx1p).le.mpx) .and. (abs(ipy1p).le.my)
c     &                .and. (ipz1p.ge.-mpz).and. (ipz1p.le.mpzp))
c     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) =
c     &                f(ix1,iy1,iz1,ipx1p,ipy1p,ipz1p) + 1.
c                  end if
c                end if
c                if ((abs(ix2).le.mx) .and. (abs(iy2).le.my) .and.
c     &              (abs(iz2).le.mz)) then
c                  ipx2p = nint(p(1,i2)/dpx)
c                  ipy2p = nint(p(2,i2)/dpy)
c                  ipz2p = nint(p(3,i2)/dpz)
c                  if ((ipx2p.ne.ipx2) .or. (ipy2p.ne.ipy2) .or.
c     &                (ipz2p.ne.ipz2)) then
c                    if ((abs(ipx2).le.mpx) .and. (abs(ipy2).le.my)
c     &                .and. (ipz2.ge.-mpz) .and. (ipz2.le.mpzp))
c     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) =
c     &                f(ix2,iy2,iz2,ipx2,ipy2,ipz2) - 1.
c                    if ((abs(ipx2p).le.mpx) .and. (abs(ipy2p).le.my)
c     &                .and. (ipz2p.ge.-mpz) .and. (ipz2p.le.mpzp))
c     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) =
c     &                f(ix2,iy2,iz2,ipx2p,ipy2p,ipz2p) + 1.
c                  end if
c                end if
clin-10/25/02-end
90005       continue
       RETURN
*-----------------------------------------------------------------------
*COM: SET THE NEW MOMENTUM COORDINATES
107     IF(PX .EQ. 0.0 .AND. PY .EQ. 0.0) THEN
        T2 = 0.0
      ELSE
        T2=ATAN2(PY,PX)
      END IF
      S1   = 1.0 - C1**2 
       IF(S1.LE.0)S1=0
       S1=SQRT(S1)
clin-9/2012: check argument in sqrt():
       scheck=1.0 - C2**2
       if(scheck.lt.0) then
          write(99,*) 'scheck3: ', scheck
          scheck=0.
       endif
       S2=SQRT(scheck)
c       S2  =  SQRT( 1.0 - C2**2 )
      CT1  = COS(T1)
      ST1  = SIN(T1)
      CT2  = COS(T2)
      ST2  = SIN(T2)
      PZ   = PR * ( C1*C2 - S1*S2*CT1 )
      SS   = C2 * S1 * CT1  +  S2 * C1
      PX   = PR * ( SS*CT2 - S1*ST1*ST2 )
      PY   = PR * ( SS*ST2 + S1*ST1*CT2 )
      RETURN
      END
clin-5/2008 CRNN over
**********************************
**********************************
*                                                                      *
*                                                                      *
c
