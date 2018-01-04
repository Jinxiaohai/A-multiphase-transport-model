      real function Xnp(srt)
*  This function contains the experimental total n-p cross sections           *
*  srt    = DSQRT(s) in GeV                                                   *
*  xsec   = production cross section in mb                                    *
*  earray = EXPerimental table with proton energies in MeV                    *
*  xarray = EXPerimental table with cross sections in mb (curve to guide eye) *
*  WITH  A CUTOFF AT 55MB                                                *
******************************************
c      real*4   xarray(11), earray(11)
      real   xarray(11), earray(11)
      SAVE   
      data   earray /20.,30.,40.,60.,90.,135.0,200.,
     &300.,400.,600.,800./
      data  xarray / 410.,270.,214.5,130.,78.,53.5,
     &41.6,35.9,34.2,34.3,34.9/
       pmass=0.9383
* 1.Calculate E_kin(lab) [MeV] from srt [GeV]
*   Formula used:   DSQRT(s) = 2 m DSQRT(E_kin/(2m) + 1)
      ekin = 2000.*pmass*((srt/(2.*pmass))**2 - 1.)
      if (ekin .lt. earray(1)) then
        xnp = xarray(1)
       IF(XNP.GT.55)XNP=55
        return
      end if
       IF(EKIN.GT.EARRAY(11))THEN
       XNP=XARRAY(11)
       RETURN
       ENDIF
*
*Interpolate double logarithmically to find sigma(srt)
*
      do 1001 ie = 1,11
        if (earray(ie) .eq. ekin) then
          xNP = xarray(ie)
         if(xnp.gt.55)xnp=55.
          return
       endif
        if (earray(ie) .gt. ekin) then
          ymin = alog(xarray(ie-1))
          ymax = alog(xarray(ie))
          xmin = alog(earray(ie-1))
          xmax = alog(earray(ie))
          xNP = exp(ymin + (alog(ekin)-xmin)
     &          *(ymax-ymin)/(xmax-xmin) )
       IF(XNP.GT.55)XNP=55
       go to 50
        end if
 1001 continue
50       continue
        return
        END
*******************************
