      subroutine addhad
      PARAMETER (MAXSTR=150001,MAXR=1,xmd=1.8756)
      double precision  smearp,smearh
      COMMON /ARPRNT/ ARPAR1(100), IAPAR2(50), ARINT1(100), IAINT2(50)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON /dpert/dpertt(MAXSTR,MAXR),dpertp(MAXSTR),dplast(MAXSTR),
     1     dpdcy(MAXSTR),dpdpi(MAXSTR,MAXR),dpt(MAXSTR, MAXR),
     2     dpp1(MAXSTR,MAXR),dppion(MAXSTR,MAXR)
      COMMON /smearz/smearp,smearh
      COMMON/RNDF77/NSEED
      common /para8/ idpert,npertd,idxsec
      SAVE   
c     All hadrons at the start of hadron cascade have the weight of 1
c     except those inserted by the user in this subroutine:
      np0=IAINT2(1)
      DO i=1,np0
         dpertp(I)=1.
      ENDDO
c     Specify number, species, weight, initial x,p,m for inserted hadrons here:
      nadd=0
      tau0=ARPAR1(1)
      DO 100 i=np0+1,np0+nadd
         ITYPAR(I)=42
clin-5/2012 fix type mismatch:
c         dpertp(I)=1d0/dble(nadd)
         dpertp(I)=1./float(nadd)
         GXAR(I)=5.*(1.-2.*RANART(NSEED))
         GYAR(I)=5.*(1.-2.*RANART(NSEED))
         GZAR(I)=2.*(1.-2.*RANART(NSEED))
         FTAR(I)=0.
         PXAR(I)=1.
         PYAR(I)=0.
         PZAR(I)=1.
         XMAR(I)=xmd
c
         PEAR(I)=sqrt(PXAR(I)**2+PYAR(I)**2+PZAR(I)**2+XMAR(I)**2)
clin-9/2012 determine rapidity more generally:
c         RAP=0.5*alog((PEAR(I)+PZAR(I)+1e-5)/(PEAR(I)-PZAR(I)+1e-5))
         RAP=asinh(PZAR(I)/sqrt(XMAR(I)**2+PXAR(I)**2+PYAR(I)**2))
c
         VX=PXAR(I)/PEAR(I)
         VY=PYAR(I)/PEAR(I)
c.....give initial formation time shift and boost according to rapidity:
         TAUI=FTAR(I)+TAU0
         FTAR(I)=TAUI*COSH(RAP)
         GXAR(I)=GXAR(I)+VX*TAU0*COSH(RAP)
         GYAR(I)=GYAR(I)+VY*TAU0*COSH(RAP)
c     Allow the intial z-position to be different from the Bjorken picture:
         GZAR(I)=TAUI*SINH(RAP)+GZAR(I)
c         GZAR(I)=TAUI*SINH(RAP)
         zsmear=sngl(smearh)*(2.*RANART(NSEED)-1.)
         GZAR(I)=GZAR(I)+zsmear
 100  CONTINUE
      IAINT2(1)=IAINT2(1)+nadd
c
      if(nadd.ge.1.and.idpert.ne.1.and.idpert.ne.2) then
         write(16,*) 'IDPERT must be 1 or 2 to add initial hadrons,
     1 set NPERTD to 0 if you do not need perturbative deuterons'
         stop
      endif
      if(IAINT2(1).gt.MAXSTR) then
         write(16,*) 'Too many initial hadrons, array size is exceeded!'
         stop
      endif
c
      return
      end
clin-8/2014 define function asinh():
