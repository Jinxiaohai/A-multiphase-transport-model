!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!.................... linana.f
!=======================================================================
!     10/26/01 update freezeout positions in case of interactions:
!lin-3/2009 Note: freezeout spacetime values cannot be trusted for K0S & K0L
!     as K0S/K0L are converted from K+/K- by hand at the end of hadron cascade.
Subroutine hbtout(nnew, nt, ntmax)
  !
  Parameter (maxstr=150001, maxr=1)
  !lin-5/2008 give tolerance to regular particles (perturbative probability 1):
  Parameter (oneminus=0.99999, oneplus=1.00001)
  Dimension lastkp(maxstr), newkp(maxstr), xnew(3)
  Common /para7/ioscar, nsmbbbar, nsmmeson
  !c      SAVE /para7/
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  !c      SAVE /hbt/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /aa/r(3, maxstr)
  !c      SAVE /AA/
  Common /bb/p(3, maxstr)
  !c      SAVE /BB/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
  Common /lastt/itimeh, bimp
  !c      SAVE /lastt/
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  !c      SAVE /tdecay/
  Common /arevt/iaevt, iarun, miss
  !c      SAVE /AREVT/
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  !c      SAVE /snn/
  Common /hjglbr/nelt, ninthj, nelp, ninp
  !c      SAVE /HJGLBR/
  Common /ftmax/ftsv(maxstr), ftsvt(maxstr, maxr)
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  !lin-12/14/03:
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  External iarflv, invflv
  Common /para8/idpert, npertd, idxsec
  !lin-2/2012:
  Common /phihj/iphirp, phirp
  Save
  !
  Do i = 1, max0(nlast, nnew)
     lastkp(i) = 0
  End Do
  Do i = 1, nnew
     newkp(i) = 0
  End Do
  !     for each of the particles, search the freezeout record (common /hbt/)
  !     to find & keep those which do not have interactions during this timestep:
  Do ip = 1, nnew
     Do iplast = 1, nlast
        If (p(1,ip)==plast(1,iplast) .And. p(2,ip)==plast(2,iplast) .And. p(3,ip)==plast(3,iplast) .And. e(ip)==plast(4,iplast) .And. lb(ip)==lblast(iplast) .And. dpertp(ip)==dplast(iplast) .And. lastkp(iplast)==0) Then
           !lin-5/2008 modified below to the above in case we have perturbative particles:
           !     5           lastkp(iplast).eq.0) then
           deltat = nt*dt - xlast(4, iplast)
           ene = sqrt(plast(1,iplast)**2+plast(2,iplast)**2+plast(3,iplast)**2+plast(4,iplast)**2)
           !     xnew gives the coordinate if a particle free-streams to current time:
           Do ii = 1, 3
              xnew(ii) = xlast(ii, iplast) + plast(ii, iplast)/ene*deltat
           End Do
           dr = sqrt((r(1,ip)-xnew(1))**2+(r(2,ip)-xnew(2))**2+(r(3,ip)-xnew(3))**2)
           !     find particles with dp=0 and dr<0.01, considered to be those
           !     without any interactions during this timestep,
           !     thus keep their last positions and time:
           If (dr<=0.01) Then
              lastkp(iplast) = 1
              newkp(ip) = 1
              !                  if(lb(ip).eq.41) then
              !                write(95,*) 'nt,ip,px,x=',nt,ip,p(1,ip),r(1,ip),ftsv(ip)
              !                write(95,*) 'xnew=',xnew(1),xnew(2),xnew(3),xlast(4,ip)
              !                  endif
              !lin-5/2009 Take care of formation time of particles read in at nt=ntmax-1:
              If (nt==ntmax .And. ftsv(ip)>((ntmax-1)*dt)) xlast(4, iplast) = ftsv(ip)
              Goto 100
           End If
        End If
     End Do
100 End Do
  !     for current particles with interactions, fill their current info in
  !     the freezeout record (if that record entry needs not to be kept):
  Do ip = 1, nnew
     If (newkp(ip)==0) Then
        Do iplast = 1, nnew
           If (lastkp(iplast)==0) Then
              !test off: write collision info
              !                  if(lb(ip).eq.41) then
              !                     write(95,*) 'nt,lb(ip)=',nt,lb(ip)
              !                  write(95,*) '  last p=',plast(1,iplast),
              !     1 plast(2,iplast),plast(3,iplast),plast(4,iplast)
              !                  write(95,*) '  after p=',p(1,ip),p(2,ip),p(3,ip),e(ip)
              !                  write(95,*) 'after x=',r(1,ip),r(2,ip),r(3,ip),ftsv(ip)
              !                  endif
              !
              xlast(1, iplast) = r(1, ip)
              xlast(2, iplast) = r(2, ip)
              xlast(3, iplast) = r(3, ip)
              xlast(4, iplast) = nt*dt
              !
              If (nt==ntmax) Then
                 !     freezeout time for decay daughters at the last timestep
                 !     needs to include the decay time of the parent:
                 If (tfdcy(ip)>(ntmax*dt+0.001)) Then
                    xlast(4, iplast) = tfdcy(ip)
                    !     freezeout time for particles unformed at the next-to-last timestep
                    !     needs to be their formation time instead of (ntmax*dt):
                 Else If (ftsv(ip)>((ntmax-1)*dt)) Then
                    xlast(4, iplast) = ftsv(ip)
                 End If
              End If
              plast(1, iplast) = p(1, ip)
              plast(2, iplast) = p(2, ip)
              plast(3, iplast) = p(3, ip)
              plast(4, iplast) = e(ip)
              lblast(iplast) = lb(ip)
              lastkp(iplast) = 1
              !lin-5/2008:
              dplast(iplast) = dpertp(ip)
              Goto 150
           End If
        End Do
     End If
150 End Do
  !     if the current particle list is shorter than the freezeout record,
  !     condense the last-collision record by filling new record from 1 to nnew,
  !     and label these entries as keep:
  If (nnew<nlast) Then
     Do iplast = 1, nlast
        If (lastkp(iplast)==0) Then
           Do ip2 = iplast + 1, nlast
              If (lastkp(ip2)==1) Then
                 xlast(1, iplast) = xlast(1, ip2)
                 xlast(2, iplast) = xlast(2, ip2)
                 xlast(3, iplast) = xlast(3, ip2)
                 xlast(4, iplast) = xlast(4, ip2)
                 plast(1, iplast) = plast(1, ip2)
                 plast(2, iplast) = plast(2, ip2)
                 plast(3, iplast) = plast(3, ip2)
                 plast(4, iplast) = plast(4, ip2)
                 lblast(iplast) = lblast(ip2)
                 lastkp(iplast) = 1
                 !lin-5/2008:
                 dplast(iplast) = dplast(ip2)
                 Goto 170
              End If
           End Do
        End If
170  End Do
  End If
  nlast = nnew
  !test off look inside each NT timestep (for debugging purpose):
  !      do ip=1,nlast
  !         write(99,*) ' p ',nt,ip,lblast(ip),plast(1,ip),
  !     1        plast(2,ip),plast(3,ip),plast(4,ip),dplast(ip)
  !         write(99,*) '  x ',nt,ip,lblast(ip),xlast(1,ip),
  !     1        xlast(2,ip),xlast(3,ip),xlast(4,ip),dplast(ip)
  !      enddo
  !
  If (nt==ntmax) Then
     !lin-5/2008 find final number of perturbative particles (deuterons only):
     ndpert = 0
     Do ip = 1, nlast
        If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
        Else
           ndpert = ndpert + 1
        End If
     End Do
     !
     !         write(16,190) IAEVT,IARUN,nlast,bimp,npart1,npart2,
     !     1 NELP,NINP,NELT,NINTHJ
     !lin-2/2012:
     !         write(16,190) IAEVT,IARUN,nlast-ndpert,bimp,npart1,npart2,
     !     1 NELP,NINP,NELT,NINTHJ
     Write (16, 191) iaevt, iarun, nlast - ndpert, bimp, npart1, npart2, nelp, ninp, nelt, ninthj, phirp
     !lin-5/2008 write out perturbatively-produced particles (deuterons only):
     If (idpert==1 .Or. idpert==2) Write (90, 190) iaevt, iarun, ndpert, bimp, npart1, npart2, nelp, ninp, nelt, ninthj
     Do ip = 1, nlast
        !lin-12/14/03   No formation time for spectator projectile or target nucleons,
        !     see ARINI1 in 'amptsub.f':
        !lin-3/2009 To be consistent with new particles produced in hadron cascade
        !     that are limited by the time-resolution (DT) of the hadron cascade,
        !     freezeout time of spectator projectile or target nucleons is written as
        !     DT as they are read at the 1st timestep and then propagated to time DT:
        !
        !lin-9/2011 determine spectator nucleons consistently
        !            if(plast(1,ip).eq.0.and.plast(2,ip).eq.0
        !     1           .and.(sqrt(plast(3,ip)**2+plast(4,ip)**2)*2/HINT1(1))
        !     2           .gt.0.99.and.(lblast(ip).eq.1.or.lblast(ip).eq.2)) then
        If (abs(plast(1,ip))<=epsipt .And. abs(plast(2,ip))<=epsipt .And. (plast(3,ip)>amax1(0.,pzproj-epsipz) .Or. plast(3,ip)<(-pztarg+epsipz)) .And. (lblast(ip)==1 .Or. lblast(ip)==2)) Then
           !lin-5/2008 perturbatively-produced particles (currently only deuterons)
           !     are written to ana/ampt_pert.dat (without the column for the mass);
           !     ana/ampt.dat has regularly-produced particles (including deuterons);
           !     these two sets of deuteron data are close to each other(but not the same
           !     because of the bias from triggering the perturbative production);
           !     ONLY use one data set for analysis to avoid double-counting:
           If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
              Write (16, 200) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
              !lin-12/14/03-end
           Else
              If (idpert==1 .Or. idpert==2) Then
                 Write (90, 250) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
              Else
                 Write (99, *) 'Unexpected perturbative particles'
              End If
           End If
        Else If (amax1(abs(xlast(1,ip)),abs(xlast(2,ip)),abs(xlast(3,ip)),abs(xlast(4,ip)))<9999) Then
           If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
              Write (16, 200) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
           Else
              If (idpert==1 .Or. idpert==2) Then
                 Write (90, 250) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip), dplast(ip)
              Else
                 Write (99, *) 'Unexpected perturbative particles'
              End If
           End If
        Else
           !     change format for large numbers:
           If (dplast(ip)>oneminus .And. dplast(ip)<oneplus) Then
              Write (16, 201) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), plast(4, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip)
           Else
              If (idpert==1 .Or. idpert==2) Then
                 Write (90, 251) invflv(lblast(ip)), plast(1, ip), plast(2, ip), plast(3, ip), xlast(1, ip), xlast(2, ip), xlast(3, ip), xlast(4, ip), dplast(ip)
              Else
                 Write (99, *) 'Unexpected perturbative particles'
              End If
           End If
        End If
     End Do
     If (ioscar==1) Call hoscar
  End If
  !
  Return
190 Format (3(I7), F10.4, 5X, 6(I4))
191 Format (3(I7), F10.4, 5X, 6(I4), 5X, F7.4)
  !lin-3/2009 improve the output accuracy of Pz
200 Format (I6, 2(1X,F8.3), 1X, F11.4, 1X, F6.3, 4(1X,F8.2))
201 Format (I6, 2(1X,F8.3), 1X, F11.4, 1X, F6.3, 4(1X,E8.2))
250 Format (I5, 2(1X,F8.3), 1X, F10.3, 2(1X,F7.1), 1X, F8.2, 1X, F7.2, 1X, E10.4)
251 Format (I5, 2(1X,F8.3), 1X, F10.3, 4(1X,E8.2), 1X, E10.4)
End Subroutine hbtout

!=======================================================================
Subroutine decomp(px0, py0, pz0, xm0, i, itq1)
  !
  Implicit Double Precision (D)
  Double Precision enenew, pxnew, pynew, pznew
  !lin-8/2015 changed ptwo(2,5) and related variables to double precision
  !     to avoid IEEE_DIVIDE_BY_ZERO or IEEE_INVALID or IEEE_OVERFLOW_FLAG:
  Double Precision de0, beta2, gam, ptwo, px0, py0, pz0, xm0
  Common /lor/enenew, pxnew, pynew, pznew
  !c      SAVE /lor/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  !c      SAVE /HPARNT/
  Common /decom/ptwo(2, 5)
  !c      SAVE /decom/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  Save
  !
  dcth = dble(ranart(nseed))*2.D0 - 1.D0
  dphi = dble(ranart(nseed)*hipr1(40))*2.D0
  !lin-6/2009 Added if embedding a high-Pt quark pair after string melting:
  If (iembed>=1 .And. iembed<=4) Then
     !     Decompose the parent high-Pt pion to q and qbar with an internal momentum
     !     parallel to the pion direction so that one parton has ~the same hight Pt
     !     and the other parton has a very soft Pt:
     !     Note: htop() decomposes a meson to q as it(1) followed by qbar as it(2):
     If (i==(natt-2*nsembd) .Or. i==(natt-2*nsembd-1)) Then
        dcth = 0.D0
        dphi = dble(phidecomp)
     End If
  End If
  !
  ds = xm0**2
  dpcm = dsqrt((ds-(ptwo(1,5)+ptwo(2,5))**2)*(ds-(ptwo(1,5)-ptwo(2,5))**2)/ds/4D0)
  dpz = dpcm*dcth
  dpx = dpcm*dsqrt(1.D0-dcth**2)*dcos(dphi)
  dpy = dpcm*dsqrt(1.D0-dcth**2)*dsin(dphi)
  de1 = dsqrt(ptwo(1,5)**2+dpcm**2)
  de2 = dsqrt(ptwo(2,5)**2+dpcm**2)
  !
  de0 = dsqrt(px0**2+py0**2+pz0**2+xm0**2)
  dbex = px0/de0
  dbey = py0/de0
  dbez = pz0/de0
  !     boost the reference frame up by beta (pznew=gam(pz+beta e)):
  beta2 = dbex**2 + dbey**2 + dbez**2
  gam = 1.D0/dsqrt(1.D0-beta2)
  If (beta2>=0.9999999999999D0) Then
     Write (6, *) '1', dbex, dbey, dbez, beta2, gam
  End If
  !
  Call lorenz(de1, dpx, dpy, dpz, -dbex, -dbey, -dbez)
  ptwo(1, 1) = pxnew
  ptwo(1, 2) = pynew
  ptwo(1, 3) = pznew
  ptwo(1, 4) = enenew
  Call lorenz(de2, -dpx, -dpy, -dpz, -dbex, -dbey, -dbez)
  ptwo(2, 1) = pxnew
  ptwo(2, 2) = pynew
  ptwo(2, 3) = pznew
  ptwo(2, 4) = enenew
  !
  Return
End Subroutine decomp

!=======================================================================
Subroutine htop
  !
  Parameter (maxstr=150001)
  Parameter (maxptn=400001)
  Parameter (maxidl=4001)
  Double Precision gx0, gy0, gz0, ft0, px0, py0, pz0, e0, xmass0
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs, ptwo, xmdq, ptwox, ptwoy, ptwoz
  Dimension it(4)
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  !c      SAVE /HMAIN2/
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  !c      SAVE /HMAIN1/
  Common /para1/mul
  !c      SAVE /PARA1/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  !c      SAVE /prec1/
  Common /ilist7/lstrg0(maxptn), lpart0(maxptn)
  !c      SAVE /ilist7/
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  !c      SAVE /ARPRC/
  Common /decom/ptwo(2, 5)
  !c      SAVE /decom/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /noprec/nnozpc, itypn(maxidl), gxn(maxidl), gyn(maxidl), gzn(maxidl), ftn(maxidl), pxn(maxidl), pyn(maxidl), pzn(maxidl), een(maxidl), xmn(maxidl)
  !c      SAVE /NOPREC/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  !c      SAVE /HPARNT/
  !     7/20/01: use double precision
  !     otherwise sometimes beta>1 and gamma diverge in lorenz():
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  !c      SAVE /SOFT/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  !lin-8/2015:
  Double Precision vxp0, vyp0, vzp0, xstrg0, ystrg0, xstrg, ystrg
  Common /precpa/vxp0(maxptn), vyp0(maxptn), vzp0(maxptn), xstrg0(maxptn), ystrg0(maxptn), xstrg(maxptn), ystrg(maxptn), istrg0(maxptn), istrg(maxptn)
  !      DOUBLE PRECISION  vxp0,vyp0,vzp0
  !      common /precpa/ vxp0(MAXPTN), vyp0(MAXPTN), vzp0(MAXPTN)
  !c      SAVE /precpa/
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /arevt/iaevt, iarun, miss
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  Save
  !
  npar = 0
  nnozpc = 0
  !lin-5b/2008 calculate the number of hadrons to be converted to q/qbar:
  If ((isoft==4 .Or. isoft==5) .And. (ioscar==2 .Or. ioscar==3)) Then
     nsmbbbar = 0
     nsmmeson = 0
     Do i = 1, natt
        id = itypar(i)
        idabs = iabs(id)
        i2 = mod(idabs/10, 10)
        !lin-9/2011 determine spectator nucleons consistently
        !              if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
        !     1             .ge.(HINT1(1)/2*0.99).and.
        !     2             .and.(id.eq.2112.or.id.eq.2212)) then
        If (abs(pxar(i))<=epsipt .And. abs(pyar(i))<=epsipt .And. (pzar(i)>amax1(0.,pzproj-epsipz) .Or. pzar(i)<(-pztarg+epsipz)) .And. (id==2112 .Or. id==2212)) Then
           !     spectator proj or targ nucleons without interactions, do not enter ZPC:
        Else If (idabs>1000 .And. i2/=0) Then
           !     baryons to be converted to q/qbar:
           nsmbbbar = nsmbbbar + 1
        Else If ((idabs>100 .And. idabs<1000) .Or. idabs>10000) Then
           !     mesons to be converted to q/qbar:
           nsmmeson = nsmmeson + 1
        End If
     End Do

     !lin-6/2009:
     If (ioscar==2 .Or. ioscar==3) Then
        Write (92, *) iaevt, miss, 3*nsmbbbar + 2*nsmmeson, nsmbbbar, nsmmeson, natt, natt - nsmbbbar - nsmmeson
     End If
     !           write(92,*) iaevt, 3*nsmbbbar+2*nsmmeson
     !           write(92,*) ' event#, total # of initial partons after string
     !     1 melting'
     !           write(92,*) 'String melting converts ',nsmbbbar, ' baryons &'
     !     1, nsmmeson, ' mesons'
     !           write(92,*) 'Total # of initial particles= ',natt
     !           write(92,*) 'Total # of initial particles (gamma,e,muon,...)
     !     1 not entering ZPC= ',natt-nsmbbbar-nsmmeson
  End If
  !lin-5b/2008-over
  Do i = 1, natt
     id = itypar(i)
     idabs = iabs(id)
     i4 = mod(idabs/1000, 10)
     i3 = mod(idabs/100, 10)
     i2 = mod(idabs/10, 10)
     i1 = mod(idabs, 10)
     rnum = ranart(nseed)
     ftime = 0.197*pear(i)/(pxar(i)**2+pyar(i)**2+xmar(i)**2)
     inozpc = 0
     it(1) = 0
     it(2) = 0
     it(3) = 0
     it(4) = 0
     !
     !lin-9/2011 determine spectator nucleons consistently
     !           if(PXAR(i).eq.0.and.PYAR(i).eq.0.and.PEAR(i)
     !     1 .ge.(HINT1(1)/2*0.99).and.((id.eq.2112).or.(id.eq.2212))) then
     If (abs(pxar(i))<=epsipt .And. abs(pyar(i))<=epsipt .And. (pzar(i)>amax1(0.,pzproj-epsipz) .Or. pzar(i)<(-pztarg+epsipz)) .And. (id==2112 .Or. id==2212)) Then
        !     spectator proj or targ nucleons without interactions, do not enter ZPC:
        inozpc = 1
     Else If (idabs>1000 .And. i2/=0) Then
        !     baryons:
        If (((i4==1 .Or. i4==2) .And. i4==i3) .Or. (i4==3 .And. i3==3)) Then
           If (i1==2) Then
              If (rnum<=(1./2.)) Then
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 1
              Else If (rnum<=(2./3.)) Then
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              Else
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              End If
           Else If (i1==4) Then
              If (rnum<=(2./3.)) Then
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              Else
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              End If
           End If
        Else If (i4==1 .Or. i4==2) Then
           If (i1==2) Then
              If (rnum<=(1./2.)) Then
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 1
              Else If (rnum<=(2./3.)) Then
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              Else
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              End If
           Else If (i1==4) Then
              If (rnum<=(2./3.)) Then
                 it(1) = i2
                 it(2) = i4*1000 + i3*100 + 3
              Else
                 it(1) = i4
                 it(2) = i3*1000 + i2*100 + 3
              End If
           End If
        Else If (i4>=3) Then
           it(1) = i4
           If (i3<i2) Then
              it(2) = i2*1000 + i3*100 + 1
           Else
              it(2) = i3*1000 + i2*100 + 3
           End If
        End If
        !       antibaryons:
        If (id<0) Then
           it(1) = -it(1)
           it(2) = -it(2)
        End If
        !     isoft=4or5 decompose diquark flavor it(2) to two quarks it(3)&(4):
        If (isoft==4 .Or. isoft==5) Then
           it(3) = mod(it(2)/1000, 10)
           it(4) = mod(it(2)/100, 10)
        End If

     Else If ((idabs>100 .And. idabs<1000) .Or. idabs>10000) Then
        !     mesons:
        If (i3==i2) Then
           If (i3==1 .Or. i3==2) Then
              If (rnum<=0.5) Then
                 it(1) = 1
                 it(2) = -1
              Else
                 it(1) = 2
                 it(2) = -2
              End If
           Else
              it(1) = i3
              it(2) = -i3
           End If
        Else
           If ((isign(1,id)*(-1)**i3)==1) Then
              it(1) = i3
              it(2) = -i2
           Else
              it(1) = i2
              it(2) = -i3
           End If
        End If
     Else
        !     save other particles (leptons and photons) outside of ZPC:
        inozpc = 1
     End If
     !
     If (inozpc==1) Then
        njsgs(i) = 0
        nnozpc = nnozpc + 1
        itypn(nnozpc) = itypar(i)
        pxn(nnozpc) = pxar(i)
        pyn(nnozpc) = pyar(i)
        pzn(nnozpc) = pzar(i)
        een(nnozpc) = pear(i)
        xmn(nnozpc) = xmar(i)
        gxn(nnozpc) = gxar(i)
        gyn(nnozpc) = gyar(i)
        gzn(nnozpc) = gzar(i)
        ftn(nnozpc) = ftar(i)
     Else
        njsgs(i) = 2
        ptwo(1, 5) = dble(ulmass(it(1)))
        ptwo(2, 5) = dble(ulmass(it(2)))
        Call decomp(dble(patt(i,1)), dble(patt(i,2)), dble(patt(i,3)), dble(xmar(i)), i, it(1))
        ipamax = 2
        If ((isoft==4 .Or. isoft==5) .And. iabs(it(2))>1000) ipamax = 1
        Do ipar = 1, ipamax
           npar = npar + 1
           ityp0(npar) = it(ipar)
           px0(npar) = ptwo(ipar, 1)
           py0(npar) = ptwo(ipar, 2)
           pz0(npar) = ptwo(ipar, 3)
           e0(npar) = ptwo(ipar, 4)
           xmass0(npar) = ptwo(ipar, 5)
           gx0(npar) = dble(gxar(i))
           gy0(npar) = dble(gyar(i))
           gz0(npar) = dble(gzar(i))
           ft0(npar) = dble(ftime)
           lstrg0(npar) = i
           lpart0(npar) = ipar
           vxp0(npar) = dble(patt(i,1)/patt(i,4))
           vyp0(npar) = dble(patt(i,2)/patt(i,4))
           vzp0(npar) = dble(patt(i,3)/patt(i,4))
           !lin-8/2015: set parent string information for this parton:
           xstrg(npar) = xstrg0(i)
           ystrg(npar) = ystrg0(i)
           istrg(npar) = istrg0(i)
        End Do
        !
        If ((isoft==4 .Or. isoft==5) .And. iabs(it(2))>1000) Then
           njsgs(i) = 3
           xmdq = ptwo(2, 5)
           ptwo(1, 5) = dble(ulmass(it(3)))
           ptwo(2, 5) = dble(ulmass(it(4)))
           !     8/19/02 avoid actual argument in common blocks of DECOMP:
           !                 call decomp(ptwo(2,1),ptwo(2,2),ptwo(2,3),xmdq)
           ptwox = ptwo(2, 1)
           ptwoy = ptwo(2, 2)
           ptwoz = ptwo(2, 3)
           Call decomp(ptwox, ptwoy, ptwoz, xmdq, i, it(1))
           !
           Do ipar = 1, 2
              npar = npar + 1
              ityp0(npar) = it(ipar+2)
              px0(npar) = ptwo(ipar, 1)
              py0(npar) = ptwo(ipar, 2)
              pz0(npar) = ptwo(ipar, 3)
              e0(npar) = ptwo(ipar, 4)
              xmass0(npar) = ptwo(ipar, 5)
              gx0(npar) = dble(gxar(i))
              gy0(npar) = dble(gyar(i))
              gz0(npar) = dble(gzar(i))
              ft0(npar) = dble(ftime)
              lstrg0(npar) = i
              lpart0(npar) = ipar + 1
              vxp0(npar) = dble(patt(i,1)/patt(i,4))
              vyp0(npar) = dble(patt(i,2)/patt(i,4))
              vzp0(npar) = dble(patt(i,3)/patt(i,4))
              !lin-8/2015: set parent string information for this parton:
              xstrg(npar) = xstrg0(i)
              ystrg(npar) = ystrg0(i)
              istrg(npar) = istrg0(i)
           End Do
        End If
        !
     End If
  End Do
  mul = npar
  !
  !lin-5b/2008:
  If ((isoft==4 .Or. isoft==5) .And. (ioscar==2 .Or. ioscar==3)) Then
     If ((natt-nsmbbbar-nsmmeson)/=nnozpc) Write (92, *) 'Problem with the total # of initial particles (gamma,e,muon,...) not entering ZPC'
     If ((3*nsmbbbar+2*nsmmeson)/=npar) Write (92, *) 'Problem with the total # of initial partons   after string melting'
  End If
  !
  Return
200 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,F8.2))
201 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 4(1X,E8.2))
End Subroutine htop

!=======================================================================
Subroutine ptoh
  !
  Parameter (maxstr=150001)
  Double Precision gxp, gyp, gzp, ftp, pxp, pyp, pzp, pep, pmp
  Double Precision gxp0, gyp0, gzp0, ft0fom, drlocl
  Double Precision enenew, pxnew, pynew, pznew, beta2, gam
  Double Precision ftavg0, gxavg0, gyavg0, gzavg0, bex, bey, bez
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Double Precision xmdiag, px1, py1, pz1, e1, px2, py2, pz2, e2, px3, py3, pz3, e3, xmpair, etot
  !lin-9/2012: improve precision for argument in sqrt():
  Double Precision p1, p2, p3
  Common /loclco/gxp(3), gyp(3), gzp(3), ftp(3), pxp(3), pyp(3), pzp(3), pep(3), pmp(3)
  !c      SAVE /loclco/
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  !c      SAVE /HMAIN1/
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  !c      SAVE /HMAIN2/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  !c      SAVE /HJJET2/
  Common /arprnt/arpar1(100), iapar2(50), arint1(100), iaint2(50)
  !c      SAVE /ARPRNT/
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  !c      SAVE /ARPRC/
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  !c      SAVE /SOFT/
  Common /rndf77/nseed
  !c      SAVE /RNDF77/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  Common /prtn23/gxp0(3), gyp0(3), gzp0(3), ft0fom
  !c      SAVE /prtn23/
  Common /nzpc/nattzp
  !c      SAVE /nzpc/
  Common /lor/enenew, pxnew, pynew, pznew
  !c      SAVE /lor/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  !c      SAVE /LUDAT1/
  !lin 4/19/2006
  Common /lastt/itimeh, bimp
  Common /hjglbr/nelt, ninthj, nelp, ninp
  Common /arevt/iaevt, iarun, miss
  Common /para7/ioscar, nsmbbbar, nsmmeson
  !lin-5/2011
  Common /input1/masspr, massta, iseed, iavoid, dt
  !
  Dimension xmdiag(maxstr), indx(maxstr), ndiag(maxstr)
  Save
  !
  Call coales
  !     obtain particle mass here without broadening by Breit-Wigner width:
  mstj24 = mstj(24)
  mstj(24) = 0
  nuudd = 0
  npich = 0
  nrhoch = 0
  ppi0 = 1.
  prho0 = 0.
  !     determine hadron flavor except (pi0,rho0,eta,omega):
  Do isg = 1, nsg
     If (njsgs(isg)/=0) Then
        natt = natt + 1
        k1 = k2sgs(isg, 1)
        k1abs = iabs(k1)
        px1 = pxsgs(isg, 1)
        py1 = pysgs(isg, 1)
        pz1 = pzsgs(isg, 1)
        k2 = k2sgs(isg, 2)
        k2abs = iabs(k2)
        px2 = pxsgs(isg, 2)
        py2 = pysgs(isg, 2)
        pz2 = pzsgs(isg, 2)
        !     5/02/01 try lowest spin states as first choices,
        !     i.e. octet baryons and pseudoscalar mesons (ibs=2*baryonspin+1):
        e1 = pesgs(isg, 1)
        e2 = pesgs(isg, 2)
        xmpair = dsqrt((e1+e2)**2-(px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2)
        ibs = 2
        imspin = 0
        If (k1==-k2 .And. iabs(k1)<=2 .And. njsgs(isg)==2) Then
           nuudd = nuudd + 1
           xmdiag(nuudd) = xmpair
           ndiag(nuudd) = natt
        End If
        k3 = 0
        If ((isoft==4 .Or. isoft==5) .And. njsgs(isg)==3) Then
           k3 = k2sgs(isg, 3)
           k3abs = iabs(k3)
           px3 = pxsgs(isg, 3)
           py3 = pysgs(isg, 3)
           pz3 = pzsgs(isg, 3)
           e3 = pesgs(isg, 3)
           xmpair = dsqrt((e1+e2+e3)**2-(px1+px2+px3)**2-(py1+py2+py3)**2-(pz1+pz2+pz3)**2)
        End If
        !*****     isoft=3 baryon decomposition is different:
        If (isoft==3 .And. (k1abs>1000 .Or. k2abs>1000)) Then
           If (k1abs>1000) Then
              kdq = k1abs
              kk = k2abs
           Else
              kdq = k2abs
              kk = k1abs
           End If
           ki = mod(kdq/1000, 10)
           kj = mod(kdq/100, 10)
           If (mod(kdq,10)==1) Then
              idqspn = 0
           Else
              idqspn = 1
           End If
           !
           If (kk>ki) Then
              ktemp = kk
              kk = kj
              kj = ki
              ki = ktemp
           Else If (kk>kj) Then
              ktemp = kk
              kk = kj
              kj = ktemp
           End If
           !
           If (ki/=kj .And. ki/=kk .And. kj/=kk) Then
              If (idqspn==0) Then
                 kf = 1000*ki + 100*kk + 10*kj + ibs
              Else
                 kf = 1000*ki + 100*kj + 10*kk + ibs
              End If
           Else If (ki==kj .And. ki==kk) Then
              !     can only be decuplet baryons:
              kf = 1000*ki + 100*kj + 10*kk + 4
           Else
              kf = 1000*ki + 100*kj + 10*kk + ibs
           End If
           !     form a decuplet baryon if the q+diquark mass is closer to its mass
           !     (and if the diquark has spin 1):
           !c     for now only include Delta, which is present in ART:
           !c                 if(idqspn.eq.1.and.MOD(kf,10).eq.2) then
           If (kf==2112 .Or. kf==2212) Then
              If (abs(sngl(xmpair)-ulmass(kf))>abs(sngl(xmpair)-ulmass(kf+2))) kf = kf + 2
           End If
           If (k1<0) kf = -kf
           !lin-6/22/01 isoft=4or5 baryons:
        Else If ((isoft==4 .Or. isoft==5) .And. njsgs(isg)==3) Then
           If (k1abs>k2abs) Then
              ki = k1abs
              kk = k2abs
           Else
              ki = k2abs
              kk = k1abs
           End If
           If (k3abs>ki) Then
              kj = ki
              ki = k3abs
           Else If (k3abs<kk) Then
              kj = kk
              kk = k3abs
           Else
              kj = k3abs
           End If
           !
           If (ki==kj .And. ki==kk) Then
              !     can only be decuplet baryons (Delta-,++, Omega):
              ibs = 4
              kf = 1000*ki + 100*kj + 10*kk + ibs
           Else If (ki/=kj .And. ki/=kk .And. kj/=kk) Then
              !     form Lambda or Sigma according to 3-quark mass,
              !     for now neglect decuplet (Sigma*0 etc) which is absent in ART:
              ibs = 2
              kf1 = 1000*ki + 100*kj + 10*kk + ibs
              kf2 = 1000*ki + 100*kk + 10*kj + ibs
              kf = kf1
              If (abs(sngl(xmpair)-ulmass(kf1))>abs(sngl(xmpair)-ulmass(kf2))) kf = kf2
           Else
              ibs = 2
              kf = 1000*ki + 100*kj + 10*kk + ibs
              !c     for now only include Delta0,+ as decuplets, which are present in ART:
              If (kf==2112 .Or. kf==2212) Then
                 If (abs(sngl(xmpair)-ulmass(kf))>abs(sngl(xmpair)-ulmass(kf+2))) kf = kf + 2
              End If
           End If
           If (k1<0) kf = -kf
           !*****     mesons:
        Else
           If (k1abs==k2abs) Then
              If (k1abs<=2) Then
                 !     treat diagonal mesons later in the subroutine:
                 kf = 0
              Else If (k1abs<=3) Then
                 !     do not form eta', only form phi from s-sbar, since no eta' in ART:
                 kf = 333
              Else
                 kf = 100*k1abs + 10*k1abs + 2*imspin + 1
              End If
           Else
              If (k1abs>k2abs) Then
                 kmax = k1abs
                 kmin = k2abs
              Else If (k1abs<k2abs) Then
                 kmax = k2abs
                 kmin = k1abs
              End If
              kf = (100*kmax+10*kmin+2*imspin+1)*isign(1, k1+k2)*(-1)**kmax
              !     form a vector meson if the q+qbar mass is closer to its mass:
              If (mod(iabs(kf),10)==1) Then
                 If (abs(sngl(xmpair)-ulmass(iabs(kf)))>abs(sngl(xmpair)-ulmass(iabs(kf)+2))) kf = (iabs(kf)+2)*isign(1, kf)
              End If
           End If
        End If
        itypar(natt) = kf
        katt(natt, 1) = kf
        If (iabs(kf)==211) Then
           npich = npich + 1
        Else If (iabs(kf)==213) Then
           nrhoch = nrhoch + 1
        End If
     End If
     !lin-7/2011-check charm hadron flavors:
     !           if(k1abs.eq.4.or.k2abs.eq.4) then
     !              if(k3.eq.0) then
     !                 write(99,*) iaevt,k1,k2,kf,xmpair,
     !     1                ULMASS(iabs(kf)),ULMASS(iabs(kf)+2),isg
     !              else
     !                 write(99,*) iaevt,k1,k2,k3,kf,xmpair,
     !     1                ULMASS(iabs(kf)),ULMASS(iabs(kf)+2),isg
     !              endif
     !           endif
     !lin-7/2011-end
  End Do
  !     assume Npi0=(Npi+ + Npi-)/2, Nrho0=(Nrho+ + Nrho-)/2 on the average:
  If (nuudd/=0) Then
     ppi0 = float(npich/2)/float(nuudd)
     prho0 = float(nrhoch/2)/float(nuudd)
  End If
  !     determine diagonal mesons (pi0,rho0,eta and omega) from uubar/ddbar:
  npi0 = 0
  Do isg = 1, nsg
     If (k2sgs(isg,1)==-k2sgs(isg,2) .And. iabs(k2sgs(isg,1))<=2 .And. njsgs(isg)==2) Then
        If (ranart(nseed)<=ppi0) npi0 = npi0 + 1
     End If
  End Do
  !
  If (nuudd>1) Then
     Call index1(maxstr, nuudd, xmdiag, indx)
  Else
     indx(1) = 1
  End If
  !
  Do ix = 1, nuudd
     iuudd = indx(ix)
     inatt = ndiag(iuudd)
     If (ix<=npi0) Then
        kf = 111
     Else If (ranart(nseed)<=(prho0/(1-ppi0+0.00001))) Then
        kf = 113
     Else
        !     at T=150MeV, thermal weights for eta and omega(spin1) are about the same:
        If (ranart(nseed)<=0.5) Then
           kf = 221
        Else
           kf = 223
        End If
     End If
     itypar(inatt) = kf
     katt(inatt, 1) = kf
  End Do
  !  determine hadron formation time, position and momentum:
  inatt = 0
  !lin-6/2009 write out parton info after coalescence:
  If (ioscar==3) Then
     Write (85, 395) iaevt, 3*nsmbbbar + 2*nsmmeson, nsmbbbar, nsmmeson, bimp, nelp, ninp, nelt, ninthj, miss
  End If
  !
  Do isg = 1, nsg
     If (njsgs(isg)/=0) Then
        inatt = inatt + 1
        k1 = k2sgs(isg, 1)
        k1abs = iabs(k1)
        px1 = pxsgs(isg, 1)
        py1 = pysgs(isg, 1)
        pz1 = pzsgs(isg, 1)
        k2 = k2sgs(isg, 2)
        k2abs = iabs(k2)
        px2 = pxsgs(isg, 2)
        py2 = pysgs(isg, 2)
        pz2 = pzsgs(isg, 2)
        e1 = pesgs(isg, 1)
        e2 = pesgs(isg, 2)
        !
        If (njsgs(isg)==2) Then
           pxar(inatt) = sngl(px1+px2)
           pyar(inatt) = sngl(py1+py2)
           pzar(inatt) = sngl(pz1+pz2)
           patt(inatt, 1) = pxar(inatt)
           patt(inatt, 2) = pyar(inatt)
           patt(inatt, 3) = pzar(inatt)
           etot = e1 + e2
           !lin-9/2012: improve precision for argument in sqrt():
           p1 = px1 + px2
           p2 = py1 + py2
           p3 = pz1 + pz2
           !
        Else If ((isoft==4 .Or. isoft==5) .And. njsgs(isg)==3) Then
           px3 = pxsgs(isg, 3)
           py3 = pysgs(isg, 3)
           pz3 = pzsgs(isg, 3)
           e3 = pesgs(isg, 3)
           pxar(inatt) = sngl(px1+px2+px3)
           pyar(inatt) = sngl(py1+py2+py3)
           pzar(inatt) = sngl(pz1+pz2+pz3)
           patt(inatt, 1) = pxar(inatt)
           patt(inatt, 2) = pyar(inatt)
           patt(inatt, 3) = pzar(inatt)
           etot = e1 + e2 + e3
           !lin-9/2012: improve precision for argument in sqrt():
           p1 = px1 + px2 + px3
           p2 = py1 + py2 + py3
           p3 = pz1 + pz2 + pz3
           !
        End If
        xmar(inatt) = ulmass(itypar(inatt))
        !lin-5/2011-add finite width to resonances (rho,omega,eta,K*,phi,Delta) after formation:
        kf = katt(inatt, 1)
        If (kf==113 .Or. abs(kf)==213 .Or. kf==221 .Or. kf==223 .Or. abs(kf)==313 .Or. abs(kf)==323 .Or. kf==333 .Or. abs(kf)==1114 .Or. abs(kf)==2114 .Or. abs(kf)==2214 .Or. abs(kf)==2224) Then
           xmar(inatt) = resmass(kf)
        End If
        !
        pear(inatt) = sqrt(pxar(inatt)**2+pyar(inatt)**2+pzar(inatt)**2+xmar(inatt)**2)
        patt(inatt, 4) = pear(inatt)
        eatt = eatt + pear(inatt)
        ipartn = njsgs(isg)
        Do i = 1, ipartn
           ftp(i) = ftsgs(isg, i)
           gxp(i) = gxsgs(isg, i)
           gyp(i) = gysgs(isg, i)
           gzp(i) = gzsgs(isg, i)
           pxp(i) = pxsgs(isg, i)
           pyp(i) = pysgs(isg, i)
           pzp(i) = pzsgs(isg, i)
           pmp(i) = pmsgs(isg, i)
           pep(i) = pesgs(isg, i)
        End Do
        Call locldr(ipartn, drlocl)
        !
        tau0 = arpar1(1)
        ftavg0 = ft0fom + dble(tau0)
        gxavg0 = 0D0
        gyavg0 = 0D0
        gzavg0 = 0D0
        Do i = 1, ipartn
           gxavg0 = gxavg0 + gxp0(i)/ipartn
           gyavg0 = gyavg0 + gyp0(i)/ipartn
           gzavg0 = gzavg0 + gzp0(i)/ipartn
        End Do
        !lin-9/2012: improve precision for argument in sqrt():
        !            bex=dble(PXAR(inatt))/etot
        !            bey=dble(PYAR(inatt))/etot
        !            bez=dble(PZAR(inatt))/etot
        bex = p1/etot
        bey = p2/etot
        bez = p3/etot
        !
        beta2 = bex**2 + bey**2 + bez**2
        gam = 1.D0/dsqrt(1.D0-beta2)
        If (beta2>=0.9999999999999D0) Then
           Write (6, *) '2', bex, bey, bez, beta2, gam
        End If
        !
        Call lorenz(ftavg0, gxavg0, gyavg0, gzavg0, -bex, -bey, -bez)
        gxar(inatt) = sngl(pxnew)
        gyar(inatt) = sngl(pynew)
        gzar(inatt) = sngl(pznew)
        ftar(inatt) = sngl(enenew)
        !lin 4/19/2006 write out parton info after coalescence:
        If (ioscar==3) Then
           Write (85, 313) k2sgs(isg, 1), px1, py1, pz1, pmsgs(isg, 1), inatt, katt(inatt, 1), xmar(inatt)
           Write (85, 312) k2sgs(isg, 2), px2, py2, pz2, pmsgs(isg, 2), inatt, katt(inatt, 1)
           If (njsgs(isg)==3) Write (85, 312) k2sgs(isg, 3), px3, py3, pz3, pmsgs(isg, 3), inatt, katt(inatt, 1)
        End If
        !
     End If
  End Do
  !     number of hadrons formed from partons inside ZPC:
  nattzp = natt
  mstj(24) = mstj24
  !
  Return
395 Format (4I8, F10.4, 5I5)
312 Format (I6, 4(1X,F10.3), 1X, I6, 1X, I6)
  !lin-5/02/2011
313 Format (I6, 4(1X,F10.3), 1X, I6, 1X, I6, 1X, F10.3)
End Subroutine ptoh


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!=======================================================================
!lin-5/2011-add finite width to resonances (rho,omega,eta,K*,phi,Delta) after formation:
Function resmass(kf)

  Parameter (arho=0.775, aomega=0.783, aeta=0.548, aks=0.894, aphi=1.019, adelta=1.232)
  Parameter (wrho=0.149, womega=0.00849, weta=1.30E-6, wks=0.0498, wphi=0.00426, wdelta=0.118)
  Common /input1/masspr, massta, iseed, iavoid, dt
  Common /rndf77/nseed
  Save

  If (kf==113 .Or. abs(kf)==213) Then
     amass = arho
     wid = wrho
  Else If (kf==221) Then
     amass = aeta
     wid = weta
  Else If (kf==223) Then
     amass = aomega
     wid = womega
  Else If (abs(kf)==313 .Or. abs(kf)==323) Then
     amass = aks
     wid = wks
  Else If (kf==333) Then
     amass = aphi
     wid = wphi
  Else If (abs(kf)==1114 .Or. abs(kf)==2114 .Or. abs(kf)==2214 .Or. abs(kf)==2224) Then
     amass = adelta
     wid = wdelta
  End If
  dmin = amass - 2*wid
  dmax = amass + 2*wid
  !     Delta mass needs to be big enough to decay to N+pi:
  If (amass==adelta) dmin = 1.078
  !
  fm = 1.
  ntry1 = 0
10 dm = ranart(nseed)*(dmax-dmin) + dmin
  ntry1 = ntry1 + 1
  fmass = (amass*wid)**2/((dm**2-amass**2)**2+(amass*wid)**2)
  !heck      write (99,*) ntry1,kf,amass,wid,fmass,DM
  If ((ranart(nseed)>fmass/fm) .And. (ntry1<=10)) Goto 10
  !
  resmass = dm

  Return
End Function resmass

!=======================================================================
Subroutine coales

  Parameter (maxstr=150001)
  Implicit Double Precision (D)
  Double Precision gxp, gyp, gzp, ftp, pxp, pyp, pzp, pep, pmp
  Dimension iover(maxstr), dp1(2:3), dr1(2:3)
  Double Precision pxsgs, pysgs, pzsgs, pesgs, pmsgs, gxsgs, gysgs, gzsgs, ftsgs
  Double Precision dpcoal, drcoal, ecritl
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  !c      SAVE /SOFT/
  Common /coal/dpcoal, drcoal, ecritl
  !c      SAVE /coal/
  Common /loclco/gxp(3), gyp(3), gzp(3), ftp(3), pxp(3), pyp(3), pzp(3), pep(3), pmp(3)
  !c      SAVE /loclco/
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  !c      SAVE /HJJET2/
  Save
  !
  Do isg = 1, nsg
     iover(isg) = 0
  End Do
  !1     meson q coalesce with all available qbar:
  Do isg = 1, nsg
     If (njsgs(isg)/=2 .Or. iover(isg)==1) Goto 150
     !     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
     If (k2sgs(isg,1)<0) Then
        Write (6, *) 'Antiquark appears in quark loop; stop'
        Stop
     End If
     !
     Do j = 1, 2
        ftp(j) = ftsgs(isg, j)
        gxp(j) = gxsgs(isg, j)
        gyp(j) = gysgs(isg, j)
        gzp(j) = gzsgs(isg, j)
        pxp(j) = pxsgs(isg, j)
        pyp(j) = pysgs(isg, j)
        pzp(j) = pzsgs(isg, j)
        pmp(j) = pmsgs(isg, j)
        pep(j) = pesgs(isg, j)
     End Do
     Call locldr(2, drlocl)
     dr0 = drlocl
     !     dp0^2 defined as (p1+p2)^2-(m1+m2)^2:
     dp0 = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     !
     Do jsg = 1, nsg
        !     skip default or unavailable antiquarks:
        If (jsg==isg .Or. iover(jsg)==1) Goto 120
        If (njsgs(jsg)==2) Then
           ipmin = 2
           ipmax = 2
        Else If (njsgs(jsg)==3 .And. k2sgs(jsg,1)<0) Then
           ipmin = 1
           ipmax = 3
        Else
           Goto 120
        End If
        Do ip = ipmin, ipmax
           dplocl = dsqrt(2*(pep(1)*pesgs(jsg,ip)-pxp(1)*pxsgs(jsg,ip)-pyp(1)*pysgs(jsg,ip)-pzp(1)*pzsgs(jsg,ip)-pmp(1)*pmsgs(jsg,ip)))
           !     skip if outside of momentum radius:
           If (dplocl>dpcoal) Goto 120
           ftp(2) = ftsgs(jsg, ip)
           gxp(2) = gxsgs(jsg, ip)
           gyp(2) = gysgs(jsg, ip)
           gzp(2) = gzsgs(jsg, ip)
           pxp(2) = pxsgs(jsg, ip)
           pyp(2) = pysgs(jsg, ip)
           pzp(2) = pzsgs(jsg, ip)
           pmp(2) = pmsgs(jsg, ip)
           pep(2) = pesgs(jsg, ip)
           Call locldr(2, drlocl)
           !     skip if outside of spatial radius:
           If (drlocl>drcoal) Goto 120
           !     q_isg coalesces with qbar_jsg:
           If ((dp0>dpcoal .Or. dr0>drcoal) .Or. (drlocl<dr0)) Then
              dp0 = dplocl
              dr0 = drlocl
              Call exchge(isg, 2, jsg, ip)
           End If
        End Do
120  End Do
     If (dp0<=dpcoal .And. dr0<=drcoal) iover(isg) = 1
150 End Do
  !
  !2     meson qbar coalesce with all available q:
  Do isg = 1, nsg
     If (njsgs(isg)/=2 .Or. iover(isg)==1) Goto 250
     !     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
     Do j = 1, 2
        ftp(j) = ftsgs(isg, j)
        gxp(j) = gxsgs(isg, j)
        gyp(j) = gysgs(isg, j)
        gzp(j) = gzsgs(isg, j)
        pxp(j) = pxsgs(isg, j)
        pyp(j) = pysgs(isg, j)
        pzp(j) = pzsgs(isg, j)
        pmp(j) = pmsgs(isg, j)
        pep(j) = pesgs(isg, j)
     End Do
     Call locldr(2, drlocl)
     dr0 = drlocl
     dp0 = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     !
     Do jsg = 1, nsg
        If (jsg==isg .Or. iover(jsg)==1) Goto 220
        If (njsgs(jsg)==2) Then
           ipmin = 1
           ipmax = 1
        Else If (njsgs(jsg)==3 .And. k2sgs(jsg,1)>0) Then
           ipmin = 1
           ipmax = 3
        Else
           Goto 220
        End If
        Do ip = ipmin, ipmax
           dplocl = dsqrt(2*(pep(2)*pesgs(jsg,ip)-pxp(2)*pxsgs(jsg,ip)-pyp(2)*pysgs(jsg,ip)-pzp(2)*pzsgs(jsg,ip)-pmp(2)*pmsgs(jsg,ip)))
           !     skip if outside of momentum radius:
           If (dplocl>dpcoal) Goto 220
           ftp(1) = ftsgs(jsg, ip)
           gxp(1) = gxsgs(jsg, ip)
           gyp(1) = gysgs(jsg, ip)
           gzp(1) = gzsgs(jsg, ip)
           pxp(1) = pxsgs(jsg, ip)
           pyp(1) = pysgs(jsg, ip)
           pzp(1) = pzsgs(jsg, ip)
           pmp(1) = pmsgs(jsg, ip)
           pep(1) = pesgs(jsg, ip)
           Call locldr(2, drlocl)
           !     skip if outside of spatial radius:
           If (drlocl>drcoal) Goto 220
           !     qbar_isg coalesces with q_jsg:
           If ((dp0>dpcoal .Or. dr0>drcoal) .Or. (drlocl<dr0)) Then
              dp0 = dplocl
              dr0 = drlocl
              Call exchge(isg, 1, jsg, ip)
           End If
        End Do
220  End Do
     If (dp0<=dpcoal .And. dr0<=drcoal) iover(isg) = 1
250 End Do
  !
  !3     baryon q (antibaryon qbar) coalesce with all available q (qbar):
  Do isg = 1, nsg
     If (njsgs(isg)/=3 .Or. iover(isg)==1) Goto 350
     ibaryn = k2sgs(isg, 1)
     !     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
     Do j = 1, 2
        ftp(j) = ftsgs(isg, j)
        gxp(j) = gxsgs(isg, j)
        gyp(j) = gysgs(isg, j)
        gzp(j) = gzsgs(isg, j)
        pxp(j) = pxsgs(isg, j)
        pyp(j) = pysgs(isg, j)
        pzp(j) = pzsgs(isg, j)
        pmp(j) = pmsgs(isg, j)
        pep(j) = pesgs(isg, j)
     End Do
     Call locldr(2, drlocl)
     dr1(2) = drlocl
     dp1(2) = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     !
     ftp(2) = ftsgs(isg, 3)
     gxp(2) = gxsgs(isg, 3)
     gyp(2) = gysgs(isg, 3)
     gzp(2) = gzsgs(isg, 3)
     pxp(2) = pxsgs(isg, 3)
     pyp(2) = pysgs(isg, 3)
     pzp(2) = pzsgs(isg, 3)
     pmp(2) = pmsgs(isg, 3)
     pep(2) = pesgs(isg, 3)
     Call locldr(2, drlocl)
     dr1(3) = drlocl
     dp1(3) = dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)-pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
     !
     Do jsg = 1, nsg
        If (jsg==isg .Or. iover(jsg)==1) Goto 320
        If (njsgs(jsg)==2) Then
           If (ibaryn>0) Then
              ipmin = 1
           Else
              ipmin = 2
           End If
           ipmax = ipmin
        Else If (njsgs(jsg)==3 .And. (ibaryn*k2sgs(jsg,1))>0) Then
           ipmin = 1
           ipmax = 3
        Else
           Goto 320
        End If
        Do ip = ipmin, ipmax
           dplocl = dsqrt(2*(pep(1)*pesgs(jsg,ip)-pxp(1)*pxsgs(jsg,ip)-pyp(1)*pysgs(jsg,ip)-pzp(1)*pzsgs(jsg,ip)-pmp(1)*pmsgs(jsg,ip)))
           !     skip if outside of momentum radius:
           If (dplocl>dpcoal) Goto 320
           ftp(2) = ftsgs(jsg, ip)
           gxp(2) = gxsgs(jsg, ip)
           gyp(2) = gysgs(jsg, ip)
           gzp(2) = gzsgs(jsg, ip)
           pxp(2) = pxsgs(jsg, ip)
           pyp(2) = pysgs(jsg, ip)
           pzp(2) = pzsgs(jsg, ip)
           pmp(2) = pmsgs(jsg, ip)
           pep(2) = pesgs(jsg, ip)
           Call locldr(2, drlocl)
           !     skip if outside of spatial radius:
           If (drlocl>drcoal) Goto 320
           !     q_isg may coalesce with q_jsg for a baryon:
           ipi = 0
           If (dp1(2)>dpcoal .Or. dr1(2)>drcoal) Then
              ipi = 2
              If ((dp1(3)>dpcoal .Or. dr1(3)>drcoal) .And. dr1(3)>dr1(2)) ipi = 3
           Else If (dp1(3)>dpcoal .Or. dr1(3)>drcoal) Then
              ipi = 3
           Else If (dr1(2)<dr1(3)) Then
              If (drlocl<dr1(3)) ipi = 3
           Else If (dr1(3)<=dr1(2)) Then
              If (drlocl<dr1(2)) ipi = 2
           End If
           If (ipi/=0) Then
              dp1(ipi) = dplocl
              dr1(ipi) = drlocl
              Call exchge(isg, ipi, jsg, ip)
           End If
        End Do
320  End Do
     If (dp1(2)<=dpcoal .And. dr1(2)<=drcoal .And. dp1(3)<=dpcoal .And. dr1(3)<=drcoal) iover(isg) = 1
350 End Do
  !
  Return
End Subroutine coales

!=======================================================================
Subroutine exchge(isg, ipi, jsg, ipj)
  !
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxstr=150001)
  Common /soft/pxsgs(maxstr, 3), pysgs(maxstr, 3), pzsgs(maxstr, 3), pesgs(maxstr, 3), pmsgs(maxstr, 3), gxsgs(maxstr, 3), gysgs(maxstr, 3), gzsgs(maxstr, 3), ftsgs(maxstr, 3), k1sgs(maxstr, 3), k2sgs(maxstr, 3), njsgs(maxstr)
  !c      SAVE /SOFT/
  Save
  !
  k1 = k1sgs(isg, ipi)
  k2 = k2sgs(isg, ipi)
  px = pxsgs(isg, ipi)
  py = pysgs(isg, ipi)
  pz = pzsgs(isg, ipi)
  pe = pesgs(isg, ipi)
  pm = pmsgs(isg, ipi)
  gx = gxsgs(isg, ipi)
  gy = gysgs(isg, ipi)
  gz = gzsgs(isg, ipi)
  ft = ftsgs(isg, ipi)
  k1sgs(isg, ipi) = k1sgs(jsg, ipj)
  k2sgs(isg, ipi) = k2sgs(jsg, ipj)
  pxsgs(isg, ipi) = pxsgs(jsg, ipj)
  pysgs(isg, ipi) = pysgs(jsg, ipj)
  pzsgs(isg, ipi) = pzsgs(jsg, ipj)
  pesgs(isg, ipi) = pesgs(jsg, ipj)
  pmsgs(isg, ipi) = pmsgs(jsg, ipj)
  gxsgs(isg, ipi) = gxsgs(jsg, ipj)
  gysgs(isg, ipi) = gysgs(jsg, ipj)
  gzsgs(isg, ipi) = gzsgs(jsg, ipj)
  ftsgs(isg, ipi) = ftsgs(jsg, ipj)
  k1sgs(jsg, ipj) = k1
  k2sgs(jsg, ipj) = k2
  pxsgs(jsg, ipj) = px
  pysgs(jsg, ipj) = py
  pzsgs(jsg, ipj) = pz
  pesgs(jsg, ipj) = pe
  pmsgs(jsg, ipj) = pm
  gxsgs(jsg, ipj) = gx
  gysgs(jsg, ipj) = gy
  gzsgs(jsg, ipj) = gz
  ftsgs(jsg, ipj) = ft
  !
  Return
End Subroutine exchge

!=======================================================================
Subroutine locldr(icall, drlocl)
  !
  Implicit Double Precision (A-H, O-Z)
  Dimension ftp0(3), pxp0(3), pyp0(3), pzp0(3), pep0(3)
  Common /loclco/gxp(3), gyp(3), gzp(3), ftp(3), pxp(3), pyp(3), pzp(3), pep(3), pmp(3)
  !c      SAVE /loclco/
  Common /prtn23/gxp0(3), gyp0(3), gzp0(3), ft0fom
  !c      SAVE /prtn23/
  Common /lor/enenew, pxnew, pynew, pznew
  !c      SAVE /lor/
  Save
  !     for 2-body kinematics:
  If (icall==2) Then
     etot = pep(1) + pep(2)
     bex = (pxp(1)+pxp(2))/etot
     bey = (pyp(1)+pyp(2))/etot
     bez = (pzp(1)+pzp(2))/etot
     !     boost the reference frame down by beta to get to the pair rest frame:
     Do j = 1, 2
        beta2 = bex**2 + bey**2 + bez**2
        gam = 1.D0/dsqrt(1.D0-beta2)
        If (beta2>=0.9999999999999D0) Then
           Write (6, *) '4', pxp(1), pxp(2), pyp(1), pyp(2), pzp(1), pzp(2), pep(1), pep(2), pmp(1), pmp(2), dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2+pmp(1)**2)/pep(1), dsqrt(pxp(1)**2+pyp(1)**2+pzp(1)**2)/pep(1)
           Write (6, *) '4a', pxp(1) + pxp(2), pyp(1) + pyp(2), pzp(1) + pzp(2), etot
           Write (6, *) '4b', bex, bey, bez, beta2, gam
        End If
        !
        Call lorenz(ftp(j), gxp(j), gyp(j), gzp(j), bex, bey, bez)
        gxp0(j) = pxnew
        gyp0(j) = pynew
        gzp0(j) = pznew
        ftp0(j) = enenew
        Call lorenz(pep(j), pxp(j), pyp(j), pzp(j), bex, bey, bez)
        pxp0(j) = pxnew
        pyp0(j) = pynew
        pzp0(j) = pznew
        pep0(j) = enenew
     End Do
     !
     If (ftp0(1)>=ftp0(2)) Then
        ilate = 1
        iearly = 2
     Else
        ilate = 2
        iearly = 1
     End If
     ft0fom = ftp0(ilate)
     !
     dt0 = ftp0(ilate) - ftp0(iearly)
     gxp0(iearly) = gxp0(iearly) + pxp0(iearly)/pep0(iearly)*dt0
     gyp0(iearly) = gyp0(iearly) + pyp0(iearly)/pep0(iearly)*dt0
     gzp0(iearly) = gzp0(iearly) + pzp0(iearly)/pep0(iearly)*dt0
     drlocl = dsqrt((gxp0(ilate)-gxp0(iearly))**2+(gyp0(ilate)-gyp0(iearly))**2+(gzp0(ilate)-gzp0(iearly))**2)
     !     for 3-body kinematics, used for baryons formation:
  Else If (icall==3) Then
     etot = pep(1) + pep(2) + pep(3)
     bex = (pxp(1)+pxp(2)+pxp(3))/etot
     bey = (pyp(1)+pyp(2)+pyp(3))/etot
     bez = (pzp(1)+pzp(2)+pzp(3))/etot
     beta2 = bex**2 + bey**2 + bez**2
     gam = 1.D0/dsqrt(1.D0-beta2)
     If (beta2>=0.9999999999999D0) Then
        Write (6, *) '5', bex, bey, bez, beta2, gam
     End If
     !     boost the reference frame down by beta to get to the 3-parton rest frame:
     Do j = 1, 3
        Call lorenz(ftp(j), gxp(j), gyp(j), gzp(j), bex, bey, bez)
        gxp0(j) = pxnew
        gyp0(j) = pynew
        gzp0(j) = pznew
        ftp0(j) = enenew
        Call lorenz(pep(j), pxp(j), pyp(j), pzp(j), bex, bey, bez)
        pxp0(j) = pxnew
        pyp0(j) = pynew
        pzp0(j) = pznew
        pep0(j) = enenew
     End Do
     !
     If (ftp0(1)>ftp0(2)) Then
        ilate = 1
        If (ftp0(3)>ftp0(1)) ilate = 3
     Else
        ilate = 2
        If (ftp0(3)>=ftp0(2)) ilate = 3
     End If
     ft0fom = ftp0(ilate)
     !
     If (ilate==1) Then
        imin = 2
        imax = 3
        istep = 1
     Else If (ilate==2) Then
        imin = 1
        imax = 3
        istep = 2
     Else If (ilate==3) Then
        imin = 1
        imax = 2
        istep = 1
     End If
     !
     Do iearly = imin, imax, istep
        dt0 = ftp0(ilate) - ftp0(iearly)
        gxp0(iearly) = gxp0(iearly) + pxp0(iearly)/pep0(iearly)*dt0
        gyp0(iearly) = gyp0(iearly) + pyp0(iearly)/pep0(iearly)*dt0
        gzp0(iearly) = gzp0(iearly) + pzp0(iearly)/pep0(iearly)*dt0
     End Do
  End If
  !
  Return
End Subroutine locldr

!=======================================================================
Subroutine hoscar
  !
  Parameter (maxstr=150001, amn=0.939457, amp=0.93828)
  Character *8 code, reffra, frame
  Character *25 amptvn
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  !c      SAVE /snn/
  Common /lastt/itimeh, bimp
  !c      SAVE /lastt/
  Common /hbt/lblast(maxstr), xlast(4, maxstr), plast(4, maxstr), nlast
  !c      SAVE /hbt/
  Common /oscar1/iap, izp, iat, izt
  !c      SAVE /oscar1/
  Common /oscar2/frame, amptvn
  !c      SAVE /oscar2/
  Save
  Data nff/0/
  !
  !       file header
  If (nff==0) Then
     Write (19, 101) 'OSCAR1997A'
     Write (19, 111) 'final_id_p_x'
     code = 'AMPT'
     If (frame=='CMS') Then
        reffra = 'nncm'
        xmp = (amp*izp+amn*(iap-izp))/iap
        xmt = (amp*izt+amn*(iat-izt))/iat
        ebeam = (efrm**2-xmp**2-xmt**2)/2./xmt
     Else If (frame=='LAB') Then
        reffra = 'lab'
        ebeam = efrm
     Else
        reffra = 'unknown'
        ebeam = 0.
     End If
     ntestp = 1
     Write (19, 102) code, amptvn, iap, izp, iat, izt, reffra, ebeam, ntestp
     nff = 1
     ievent = 1
     phi = 0.
     If (frame=='CMS') Write (19, 112) efrm
  End If
  !       comment
  !       event header
  Write (19, 103) ievent, nlast, bimp, phi
  !       particles
  Do i = 1, nlast
     ene = sqrt(plast(1,i)**2+plast(2,i)**2+plast(3,i)**2+plast(4,i)**2)
     Write (19, 104) i, invflv(lblast(i)), plast(1, i), plast(2, i), plast(3, i), ene, plast(4, i), xlast(1, i), xlast(2, i), xlast(3, i), xlast(4, i)
  End Do
  ievent = ievent + 1
  !
  Return
101 Format (A10)
111 Format (A12)
102 Format (A4, 1X, A20, 1X, '(', I3, ',', I3, ')+(', I3, ',', I3, ')', 2X, A4, 2X, E10.4, 2X, I8)
103 Format (I10, 2X, I10, 2X, F8.3, 2X, F8.3)
104 Format (I10, 2X, I10, 2X, 9(E12.6,2X))
112 Format ('# Center-of-mass energy/nucleon-pair is', F12.3, 'GeV')
End Subroutine hoscar

!=======================================================================
Subroutine getnp

  Parameter (maxstr=150001)
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  !c      SAVE /HMAIN1/
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  !c      SAVE /HMAIN2/
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  !c      SAVE /HPARNT/
  Common /snn/efrm, npart1, npart2, epsipz, epsipt, pzproj, pztarg
  !c      SAVE /snn/
  Save

  If (natt==0) Then
     npart1 = 0
     npart2 = 0
     Return
  End If
  !
  pzproj = sqrt(hint1(6)**2-hint1(8)**2)
  pztarg = sqrt(hint1(7)**2-hint1(9)**2)
  epsipz = 0.01
  !lin-9/2011-add Pt tolerance in determining spectator nucleons
  !     (affect string melting runs when LAB frame is used):
  epsipt = 1E-6
  !
  nspec1 = 0
  nspec2 = 0
  Do i = 1, natt
     !lin-9/2011 determine spectator nucleons consistently
     !           if((KATT(I,1).eq.2112.or.KATT(I,1).eq.2212)
     !     1          .and.PATT(I, 1).eq.0.and.PATT(I, 2).eq.0) then
     If ((katt(i,1)==2112 .Or. katt(i,1)==2212) .And. abs(patt(i,1))<=epsipt .And. abs(patt(i,2))<=epsipt) Then
        If (patt(i,3)>amax1(0.,pzproj-epsipz)) Then
           nspec1 = nspec1 + 1
        Else If (patt(i,3)<(-pztarg+epsipz)) Then
           nspec2 = nspec2 + 1
        End If
     End If
  End Do
  npart1 = ihnt2(1) - nspec1
  npart2 = ihnt2(3) - nspec2

  Return
End Subroutine getnp

!=======================================================================
!     2/18/03 use PYTHIA to decay eta,rho,omega,k*,phi and Delta
!     4/2012 added pi0 decay flag:
!       ipion=0: resonance or pi0 in lb(i1); >0: pi0 in lpion(ipion).
Subroutine resdec(i1, nt, nnn, wid, idecay, ipion)

  Parameter (hbarc=0.19733)
  Parameter (ak0=0.498, apich=0.140, api0=0.135, an=0.940, addm=0.02)
  Parameter (maxstr=150001, maxr=1)
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  !c      SAVE /INPUT2/
  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  !c      SAVE /LUJETS/
  Common /ludat1/mstu(200), paru(200), mstj(200), parj(200)
  !c      SAVE /LUDAT1/
  Common /ludat2/kchg(500, 3), pmas(500, 4), parf(2000), vckm(4, 4)
  !c      SAVE /LUDAT2/
  Common /ludat3/mdcy(500, 3), mdme(2000, 2), brat(2000), kfdp(2000, 5)
  !c      SAVE /LUDAT3/
  Common /cc/e(maxstr)
  !c      SAVE /CC/
  Common /ee/id(maxstr), lb(maxstr)
  !c      SAVE /EE/
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
  Common /resdcy/nsav, iksdcy
  !c      SAVE /resdcy/
  Common /leadng/lb1, px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n, dp1n
  !c      SAVE /leadng/
  External iarflv, invflv
  Common /tdecay/tfdcy(maxstr), tfdpi(maxstr, maxr), tft(maxstr)
  !c      SAVE /tdecay/
  Common /rndf77/nseed
  Common /dpert/dpertt(maxstr, maxr), dpertp(maxstr), dplast(maxstr), dpdcy(maxstr), dpdpi(maxstr, maxr), dpt(maxstr, maxr), dpp1(maxstr, maxr), dppion(maxstr, maxr)
  !c      SAVE /RNDF77/
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Save
  irun = idecay
  !lin-4/2012 for option of pi0 decay:
  If (nt==ntmax .And. ipi0dcy==1 .And. ((lb1==4 .And. ipion==0) .Or. ipion>=1)) Then
     kf = 111
     !        if(lb1.eq.0.or.lb1.eq.25.or.lb1.eq.26.or.lb1.eq.27
  Else If (lb1==0 .Or. lb1==25 .Or. lb1==26 .Or. lb1==27 .Or. lb1==28 .Or. lb1==29 .Or. iabs(lb1)==30 .Or. lb1==24 .Or. (iabs(lb1)>=6 .And. iabs(lb1)<=9) .Or. iabs(lb1)==16) Then
     kf = invflv(lb1)
  Else
     Return
  End If
  !
  ip = 1
  !     label as undecayed and the only particle in the record:
  n = 1
  k(ip, 1) = 1
  k(ip, 3) = 0
  k(ip, 4) = 0
  k(ip, 5) = 0
  !
  k(ip, 2) = kf
  !lin-4/2012 for option of pi0 decay:
  If (ipion==0) Then
     !
     p(ip, 1) = px1
     p(ip, 2) = py1
     p(ip, 3) = pz1
     !        em1a=em1
     !     eta or omega in ART may be below or too close to (pi+pi-pi0) mass,
     !     causing LUDECY error,thus increase their mass ADDM above this thresh,
     !     noting that rho (m=0.281) too close to 2pi thrshold fails to decay:
     If ((lb1==0 .Or. lb1==28) .And. em1<(2*apich+api0+addm)) Then
        em1 = 2*apich + api0 + addm
        !     rho
     Else If (lb1>=25 .And. lb1<=27 .And. em1<(2*apich+addm)) Then
        em1 = 2*apich + addm
        !     K*
     Else If (iabs(lb1)==30 .And. em1<(apich+ak0+addm)) Then
        em1 = apich + ak0 + addm
        !     Delta created in ART may be below (n+pich) mass, causing LUDECY error:
     Else If (iabs(lb1)>=6 .And. iabs(lb1)<=9 .And. em1<(apich+an+addm)) Then
        em1 = apich + an + addm
     End If
     !        if(em1.ge.(em1a+0.01)) write (6,*)
     !     1       'Mass increase in resdec():',nt,em1-em1a,lb1
     e1 = sqrt(em1**2+px1**2+py1**2+pz1**2)
     p(ip, 4) = e1
     p(ip, 5) = em1
     !lin-5/2008:
     dpdecp = dpertp(i1)
     !lin-4/2012 for option of pi0 decay:
  Else If (nt==ntmax .And. ipi0dcy==1 .And. ipion>=1) Then
     p(ip, 1) = ppion(1, ipion, irun)
     p(ip, 2) = ppion(2, ipion, irun)
     p(ip, 3) = ppion(3, ipion, irun)
     p(ip, 5) = epion(ipion, irun)
     p(ip, 4) = sqrt(p(ip,5)**2+p(ip,1)**2+p(ip,2)**2+p(ip,3)**2)
     dpdecp = dppion(ipion, irun)
     !test off
     !           write(99,*) P(IP,4), P(IP,5), dpdecp, ipion, wid
  Else
     Print *, 'stopped in resdec() a'
     Stop
  End If
  !
  Call ludecy(ip)
  !     add decay time to daughter's formation time at the last timestep:
  If (nt==ntmax) Then
     tau0 = hbarc/wid
     taudcy = tau0*(-1.)*alog(1.-ranart(nseed))
     ndaut = n - nsav
     If (ndaut<=1) Then
        Write (10, *) 'note: ndaut(<1)=', ndaut
        Call lulist(2)
        Stop
     End If
     !     lorentz boost:
     !lin-4/2012 for option of pi0 decay:
     If (ipion==0) Then
        taudcy = taudcy*e1/em1
        tfnl = tfnl + taudcy
        xfnl = xfnl + px1/e1*taudcy
        yfnl = yfnl + py1/e1*taudcy
        zfnl = zfnl + pz1/e1*taudcy
     Else If (ipion>=1) Then
        taudcy = taudcy*p(ip, 4)/p(ip, 5)
        tfnl = tfdpi(ipion, irun) + taudcy
        xfnl = rpion(1, ipion, irun) + p(ip, 1)/p(ip, 4)*taudcy
        yfnl = rpion(2, ipion, irun) + p(ip, 2)/p(ip, 4)*taudcy
        zfnl = rpion(3, ipion, irun) + p(ip, 3)/p(ip, 4)*taudcy
     Else
        Print *, 'stopped in resdec() b', ipion, wid, p(ip, 4)
        Stop
     End If
     !     at the last timestep, assign rho, K0S or eta (decay daughter)
     !     to lb(i1) only (not to lpion) in order to decay them again:
     !lin-4/2012 for option of pi0 decay:
     !           if(n.ge.(nsav+2)) then
     If (n>=(nsav+2) .And. ipion==0) Then
        Do idau = nsav + 2, n
           kdaut = k(idau, 2)
           If (kdaut==221 .Or. kdaut==113 .Or. kdaut==213 .Or. kdaut==-213 .Or. kdaut==310) Then
              !     switch idau and i1(nsav+1):
              ksave = kdaut
              pxsave = p(idau, 1)
              pysave = p(idau, 2)
              pzsave = p(idau, 3)
              esave = p(idau, 4)
              xmsave = p(idau, 5)
              k(idau, 2) = k(nsav+1, 2)
              p(idau, 1) = p(nsav+1, 1)
              p(idau, 2) = p(nsav+1, 2)
              p(idau, 3) = p(nsav+1, 3)
              p(idau, 4) = p(nsav+1, 4)
              p(idau, 5) = p(nsav+1, 5)
              k(nsav+1, 2) = ksave
              p(nsav+1, 1) = pxsave
              p(nsav+1, 2) = pysave
              p(nsav+1, 3) = pzsave
              p(nsav+1, 4) = esave
              p(nsav+1, 5) = xmsave
              !     note: phi decay may produce rho, K0s or eta, N*(1535) decay may produce
              !     eta, but only one daughter may be rho, K0s or eta:
              Goto 111
           End If
        End Do
     End If
111  Continue
     !
     enet = 0.
     Do idau = nsav + 1, n
        enet = enet + p(idau, 4)
     End Do
     !           if(abs(enet-e1).gt.0.02)
     !     1          write(93,*) 'resdec(): nt=',nt,enet-e1,lb1
  End If

  Do idau = nsav + 1, n
     kdaut = k(idau, 2)
     lbdaut = iarflv(kdaut)
     !     K0S and K0L are named K+/K- during hadron cascade, and only
     !     at the last timestep they keep their real LB # before output;
     !     K0/K0bar (from K* decay) converted to K0S and K0L at the last timestep:
     If (nt==ntmax .And. (kdaut==130 .Or. kdaut==310 .Or. iabs(kdaut)==311)) Then
        If (kdaut==130) Then
           lbdaut = 22
        Else If (kdaut==310) Then
           lbdaut = 24
        Else If (iabs(kdaut)==311) Then
           If (ranart(nseed)<0.5) Then
              lbdaut = 22
           Else
              lbdaut = 24
           End If
        End If
     End If
     !
     If (idau==(nsav+1)) Then
        !lin-4/2012 for option of pi0 decay:
        If (ipion==0) Then
           lb(i1) = lbdaut
           e(i1) = p(idau, 5)
           px1n = p(idau, 1)
           py1n = p(idau, 2)
           pz1n = p(idau, 3)
           !lin-5/2008:
           dp1n = dpdecp
        Else If (ipion>=1) Then
           lpion(ipion, irun) = lbdaut
           epion(ipion, irun) = p(idau, 5)
           ppion(1, ipion, irun) = p(idau, 1)
           ppion(2, ipion, irun) = p(idau, 2)
           ppion(3, ipion, irun) = p(idau, 3)
           rpion(1, ipion, irun) = xfnl
           rpion(2, ipion, irun) = yfnl
           rpion(3, ipion, irun) = zfnl
           tfdpi(ipion, irun) = tfnl
           dppion(ipion, irun) = dpdecp
        End If
        !
     Else
        nnn = nnn + 1
        lpion(nnn, irun) = lbdaut
        epion(nnn, irun) = p(idau, 5)
        ppion(1, nnn, irun) = p(idau, 1)
        ppion(2, nnn, irun) = p(idau, 2)
        ppion(3, nnn, irun) = p(idau, 3)
        rpion(1, nnn, irun) = xfnl
        rpion(2, nnn, irun) = yfnl
        rpion(3, nnn, irun) = zfnl
        tfdpi(nnn, irun) = tfnl
        !lin-5/2008:
        dppion(nnn, irun) = dpdecp
     End If
  End Do
  Return
End Subroutine resdec

!=======================================================================
Subroutine inidcy

  Common /lujets/n, k(9000, 5), p(9000, 5), v(9000, 5)
  !c      SAVE /LUJETS/
  Common /resdcy/nsav, iksdcy
  !c      SAVE /resdcy/
  Save
  n = 1
  nsav = n
  Return
End Subroutine inidcy

!=======================================================================
!lin-6/06/02 local parton freezeout motivated from critical density:
Subroutine local(t)
  !
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Parameter (r0=1D0)
  Common /para1/mul
  !c      SAVE /para1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  !c      SAVE /prec2/
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  !c      SAVE /frzprc/
  Common /prec4/vx(maxptn), vy(maxptn), vz(maxptn)
  !c      SAVE /prec4/
  Common /prec5/eta(maxptn), rap(maxptn), tau(maxptn)
  !c      SAVE /prec5/
  Common /coal/dpcoal, drcoal, ecritl
  !c      SAVE /coal/
  Save
  !
  Do it = 1, 301
     If (t>=tfrz(it) .And. t<tfrz(it+1)) Then
        If (it==itlast) Then
           Return
        Else
           itlast = it
           Goto 50
        End If
     End If
  End Do
  Write (1, *) 'local time out of range in LOCAL, stop', t, it
  Stop
50 Continue
  !
  Do ip = 1, mul
     !     skip partons which have frozen out:
     If (ifrz(ip)==1) Goto 200
     If (it==301) Then
        !     freezeout all the left partons beyond the time of 3000 fm:
        etcrit = 1D6
        Goto 150
     Else
        !     freezeout when transverse energy density < etcrit:
        etcrit = (ecritl*2D0/3D0)
     End If
     !     skip partons which have not yet formed:
     If (t<ft5(ip)) Goto 200
     rap0 = rap(ip)
     eta0 = eta(ip)
     x0 = gx5(ip) + vx(ip)*(t-ft5(ip))
     y0 = gy5(ip) + vy(ip)*(t-ft5(ip))
     detdy = 0D0
     Do itest = 1, mul
        !     skip self and partons which have not yet formed:
        If (itest==ip .Or. t<ft5(itest)) Goto 100
        ettest = eta(itest)
        xtest = gx5(itest) + vx(itest)*(t-ft5(itest))
        ytest = gy5(itest) + vy(itest)*(t-ft5(itest))
        drt = sqrt((xtest-x0)**2+(ytest-y0)**2)
        !     count partons within drt<1 and -1<(eta-eta0)<1:
        If (dabs(ettest-eta0)<=1D0 .And. drt<=r0) detdy = detdy + dsqrt(px5(itest)**2+py5(itest)**2+xmass5(itest)**2)*0.5D0
100  End Do
     detdy = detdy*(dcosh(eta0)**2)/(t*3.1416D0*r0**2*dcosh(rap0))
     !     when density is below critical density for phase transition, freeze out:
150  If (detdy<=etcrit) Then
        ifrz(ip) = 1
        idfrz(ip) = ityp5(ip)
        pxfrz(ip) = px5(ip)
        pyfrz(ip) = py5(ip)
        pzfrz(ip) = pz5(ip)
        efrz(ip) = e5(ip)
        xmfrz(ip) = xmass5(ip)
        If (t>ft5(ip)) Then
           gxfrz(ip) = x0
           gyfrz(ip) = y0
           gzfrz(ip) = gz5(ip) + vz(ip)*(t-ft5(ip))
           ftfrz(ip) = t
        Else
           !     if this freezeout time < formation time, use formation time & positions.
           !     This ensures the recovery of default hadron when e_crit=infty:
           gxfrz(ip) = gx5(ip)
           gyfrz(ip) = gy5(ip)
           gzfrz(ip) = gz5(ip)
           ftfrz(ip) = ft5(ip)
        End If
     End If
200 End Do
  !
  Return
End Subroutine local

!=======================================================================
!lin-6/06/02 initialization for local parton freezeout
Subroutine inifrz
  !
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Common /ilist5/ct(maxptn), ot(maxptn), tlarge
  !c      SAVE /ilist5/
  Common /frzprc/gxfrz(maxptn), gyfrz(maxptn), gzfrz(maxptn), ftfrz(maxptn), pxfrz(maxptn), pyfrz(maxptn), pzfrz(maxptn), efrz(maxptn), xmfrz(maxptn), tfrz(302), ifrz(maxptn), idfrz(maxptn), itlast
  !c      SAVE /frzprc/
  Save
  !
  !     for freezeout time 0-10fm, use interval of 0.1fm;
  !     for 10-100fm, use interval of 1fm;
  !     for 100-1000fm, use interval of 10fm;
  !     for 1000-3000fm, use interval of 100fm:
  step1 = 0.1D0
  step2 = 1D0
  step3 = 10D0
  step4 = 100D0
  !
  Do it = 1, 101
     tfrz(it) = 0D0 + dble(it-1)*step1
  End Do
  Do it = 102, 191
     tfrz(it) = 10D0 + dble(it-101)*step2
  End Do
  Do it = 192, 281
     tfrz(it) = 100D0 + dble(it-191)*step3
  End Do
  Do it = 282, 301
     tfrz(it) = 1000D0 + dble(it-281)*step4
  End Do
  tfrz(302) = tlarge
  !
  Return
End Subroutine inifrz


!// ForQuill v1.01 Beta www.fcode.cn
!// ForQuill v1.01 Beta www.fcode.cn
!=======================================================================
!     idd=0,1,2,3 specifies different subroutines for partonic flow analysis.
Subroutine flowp(idd)
  !
  Implicit Double Precision (A-H, O-Z)
  Real dt
  Parameter (maxptn=400001)
  !sp
  Parameter (bmt=0.05D0)
  Dimension nlfile(3), nsfile(3), nmfile(3)
  !
  Dimension v2pp(3), xnpp(3), v2psum(3), v2p2sm(3), nfile(3)
  Dimension tsp(31), v2pevt(3), v2pavg(3), varv2p(3)
  Common /ilist1/iscat, jscat, next(maxptn), last(maxptn), ictype, icsta(maxptn), nic(maxptn), icels(maxptn)
  !c      SAVE /ilist1/
  Common /para1/mul
  !c      SAVE /para1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  !c      SAVE /prec2/
  Common /pflow/v2p(30, 3), xnpart(30, 3), etp(30, 3), s2p(30, 3), v2p2(30, 3), nevt(30)
  !c      SAVE /pflow/
  Common /pflowf/v2pf(30, 3), xnpf(30, 3), etpf(30, 3), xncoll(30), s2pf(30, 3), v2pf2(30, 3)
  !c      SAVE /pflowf/
  Common /pfrz/v2pfrz(30, 3), xnpfrz(30, 3), etpfrz(30, 3), s2pfrz(30, 3), v2p2fz(30, 3), tscatt(31), nevtfz(30), iscatt(30)
  !c      SAVE /pfrz/
  Common /hflow/v2h(30, 3), xnhadr(30, 3), eth(30, 3), v2h2(30, 3), s2h(30, 3)
  !c      SAVE /hflow/
  Common /arevt/iaevt, iarun, miss
  !c      SAVE /AREVT/
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  !c      SAVE /INPUT2/
  !c      SAVE itimep,iaevtp,v2pp,xnpp,v2psum,v2p2sm
  !c      SAVE nfile,itanim,nlfile,nsfile,nmfile
  Common /precpb/vxp(maxptn), vyp(maxptn), vzp(maxptn)
  Save
  !sp
  Dimension etpl(30, 3), etps(30, 3), etplf(30, 3), etpsf(30, 3), etlfrz(30, 3), etsfrz(30, 3), xnpl(30, 3), xnps(30, 3), xnplf(30, 3), xnpsf(30, 3), xnlfrz(30, 3), xnsfrz(30, 3), v2pl(30, 3), v2ps(30, 3), v2plf(30, 3), v2psf(30, 3), s2pl(30, 3), s2ps(30, 3), s2plf(30, 3), s2psf(30, 3), dmyil(50, 3), dmyfl(50, 3), dmyis(50, 3), dmyfs(50, 3)
  Data tsp/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30/
  !     idd=0: initialization for flow analysis, called by artdri.f:
  If (idd==0) Then
     nfile(1) = 60
     nfile(2) = 64
     nfile(3) = 20
     Open (nfile(1), File='ana1/v2p.dat', Status='UNKNOWN')
     Open (nfile(1)+1, File='ana1/v2p-formed.dat', Status='UNKNOWN')
     Open (nfile(1)+2, File='ana1/v2p-active.dat', Status='UNKNOWN')
     Open (nfile(1)+3, File='ana1/v2ph.dat', Status='UNKNOWN')
     Open (nfile(2), File='ana1/v2p-y2.dat', Status='UNKNOWN')
     Open (nfile(2)+1, File='ana1/v2p-formed2.dat', Status='UNKNOWN')
     Open (nfile(2)+2, File='ana1/v2p-active2.dat', Status='UNKNOWN')
     Open (nfile(2)+3, File='ana1/v2ph-y2.dat', Status='UNKNOWN')
     Open (nfile(3), File='ana1/v2p-y1.dat', Status='UNKNOWN')
     Open (nfile(3)+1, File='ana1/v2p-formed1.dat', Status='UNKNOWN')
     Open (nfile(3)+2, File='ana1/v2p-active1.dat', Status='UNKNOWN')
     Open (nfile(3)+3, File='ana1/v2ph-y1.dat', Status='UNKNOWN')
     Open (49, File='ana1/v2p-ebe.dat', Status='UNKNOWN')
     Write (49, *) '    ievt,  v2p,  v2p_y2,   v2p_y1'
     !
     Open (59, File='ana1/v2h.dat', Status='UNKNOWN')
     Open (68, File='ana1/v2h-y2.dat', Status='UNKNOWN')
     Open (69, File='ana1/v2h-y1.dat', Status='UNKNOWN')
     Open (88, File='ana1/v2h-ebe.dat', Status='UNKNOWN')
     Write (88, *) '    ievt,  v2h,  v2h_y2,   v2h_y1'
     !sp07/05
     nlfile(1) = 70
     nlfile(2) = 72
     nlfile(3) = 74
     Open (nlfile(1), File='ana1/mtl.dat', Status='UNKNOWN')
     Open (nlfile(1)+1, File='ana1/mtl-formed.dat', Status='UNKNOWN')
     Open (nlfile(2), File='ana1/mtl-y2.dat', Status='UNKNOWN')
     Open (nlfile(2)+1, File='ana1/mtl-formed2.dat', Status='UNKNOWN')
     Open (nlfile(3), File='ana1/mtl-y1.dat', Status='UNKNOWN')
     Open (nlfile(3)+1, File='ana1/mtl-formed1.dat', Status='UNKNOWN')
     nsfile(1) = 76
     nsfile(2) = 78
     nsfile(3) = 80
     Open (nsfile(1), File='ana1/mts.dat', Status='UNKNOWN')
     Open (nsfile(1)+1, File='ana1/mts-formed.dat', Status='UNKNOWN')
     Open (nsfile(2), File='ana1/mts-y2.dat', Status='UNKNOWN')
     Open (nsfile(2)+1, File='ana1/mts-formed2.dat', Status='UNKNOWN')
     Open (nsfile(3), File='ana1/mts-y1.dat', Status='UNKNOWN')
     Open (nsfile(3)+1, File='ana1/mts-formed1.dat', Status='UNKNOWN')
     nmfile(1) = 82
     nmfile(2) = 83
     nmfile(3) = 84
     Open (nmfile(1), File='ana1/Nmt.dat', Status='UNKNOWN')
     Open (nmfile(2), File='ana1/Nmt-y2.dat', Status='UNKNOWN')
     Open (nmfile(3), File='ana1/Nmt-y1.dat', Status='UNKNOWN')
     !lin-8/2015: changed unit number of animation files,
     !test off     turn off animation output (0 to turn off and 1 to turn on):
     ifanim = 0
     !lin-11/27/00 for animation:
     If (ifanim==1) Then
        Open (10, File='ana1/h-animate.dat', Status='UNKNOWN')
        Write (10, *) ntmax, dt
        Open (11, File='ana1/p-animate.dat', Status='UNKNOWN')
        Open (15, File='ana1/p-finalft.dat', Status='UNKNOWN')
     End If
     !lin-10/2014: write out partons at all eta, turn off now:
     !           if(nevent.ge.1)
     If (nevent<1) Open (93, File='ana1/parton-t.dat', Status='UNKNOWN')
     !
     itimep = 0
     itanim = 0
     iaevtp = 0
     !sp
     Do ii = 1, 50
        Do iy = 1, 3
           dmyil(ii, iy) = 0D0
           dmyfl(ii, iy) = 0D0
           dmyis(ii, iy) = 0D0
           dmyfs(ii, iy) = 0D0
        End Do
     End Do
     !
     Do ii = 1, 31
        tscatt(ii) = 0D0
     End Do
     Do ii = 1, 30
        nevt(ii) = 0
        xncoll(ii) = 0D0
        nevtfz(ii) = 0
        iscatt(ii) = 0
        Do iy = 1, 3
           v2p(ii, iy) = 0D0
           v2p2(ii, iy) = 0D0
           s2p(ii, iy) = 0D0
           etp(ii, iy) = 0D0
           xnpart(ii, iy) = 0D0
           v2pf(ii, iy) = 0D0
           v2pf2(ii, iy) = 0D0
           s2pf(ii, iy) = 0D0
           etpf(ii, iy) = 0D0
           xnpf(ii, iy) = 0D0
           v2pfrz(ii, iy) = 0D0
           v2p2fz(ii, iy) = 0D0
           s2pfrz(ii, iy) = 0D0
           etpfrz(ii, iy) = 0D0
           xnpfrz(ii, iy) = 0D0
           !sp07/05
           etpl(ii, iy) = 0D0
           etps(ii, iy) = 0D0
           etplf(ii, iy) = 0D0
           etpsf(ii, iy) = 0D0
           etlfrz(ii, iy) = 0D0
           etsfrz(ii, iy) = 0D0
           xnpl(ii, iy) = 0D0
           xnps(ii, iy) = 0D0
           xnplf(ii, iy) = 0D0
           xnpsf(ii, iy) = 0D0
           xnlfrz(ii, iy) = 0D0
           xnsfrz(ii, iy) = 0D0
           v2pl(ii, iy) = 0D0
           v2ps(ii, iy) = 0D0
           v2plf(ii, iy) = 0D0
           v2psf(ii, iy) = 0D0
           s2pl(ii, iy) = 0D0
           s2ps(ii, iy) = 0D0
           s2plf(ii, iy) = 0D0
           s2psf(ii, iy) = 0D0
        End Do
     End Do
     Do iy = 1, 3
        v2pevt(iy) = 0D0
        v2pavg(iy) = 0D0
        varv2p(iy) = 0D0
        v2pp(iy) = 0.D0
        xnpp(iy) = 0D0
        v2psum(iy) = 0.D0
        v2p2sm(iy) = 0.D0
     End Do
     !     idd=1: calculate parton elliptic flow, called by zpc.f:
  Else If (idd==1) Then
     If (iaevt/=iaevtp .And. ianp==31) itanim = 0
     !
     t2time = ft5(iscat)
     Do ianp = 1, 30
        If (t2time<tsp(ianp+1) .And. t2time>=tsp(ianp)) Then
           !     write flow info only once at each fixed time:
           xncoll(ianp) = xncoll(ianp) + 1D0
           !     to prevent an earlier t2time comes later in the same event
           !     and mess up nevt:
           If (ianp<=itimep .And. iaevt==iaevtp) Goto 101
           nevt(ianp) = nevt(ianp) + 1
           tscatt(ianp+1) = t2time
           iscatt(ianp) = 1
           nevtfz(ianp) = nevtfz(ianp) + 1
           Do i = 1, mul
              !lin-8/2015 to avoid IEEE_DIVIDE_BY_ZERO or IEEE_INVALID:
              !                    ypartn=0.5d0*dlog((E5(i)+PZ5(i))
              !     1                   /(E5(i)-PZ5(i)+1.d-8))
              delta = 1D-8
              If ((e5(i)-dabs(pz5(i))+delta)<=0) Then
                 ypartn = 1000000.D0*sign(1.D0, pz5(i))
                 Write (6, *) 'ypartn error', e5(i) - dabs(pz5(i))
              Else
                 ypartn = 0.5D0*dlog((e5(i)+pz5(i)+delta)/(e5(i)-pz5(i)+delta))
              End If
              pt2 = px5(i)**2 + py5(i)**2
              !test off: initial (pt,y) and (x,y) distribution:
              !                    idtime=1
              !                    if(ianp.eq.idtime) then
              !                       iityp=iabs(ITYP5(I))
              !                       if(iityp.eq.1.or.iityp.eq.2) then
              !                          write(651,*) dsqrt(pt2),ypartn
              !                          write(654,*) GX5(I),GY5(I)
              !                       elseif(iityp.eq.1103.or.iityp.eq.2101
              !     1 .or.iityp.eq.2103.or.iityp.eq.2203.
              !     2 .or.iityp.eq.3101.or.iityp.eq.3103.
              !     3 .or.iityp.eq.3201.or.iityp.eq.3203.or.iityp.eq.3303)
              !     4 then
              !                          write(652,*) dsqrt(pt2),ypartn
              !                          write(655,*) GX5(I),GY5(I)
              !                       elseif(iityp.eq.21) then
              !                          write(653,*) dsqrt(pt2),ypartn
              !                          write(656,*) GX5(I),GY5(I)
              !                       endif
              !                    endif
              !test-end
              !test off density with 2fm radius and z:(-0.1*t,0.1*t):
              !                    gx_now=GX5(i)+(t2time-FT5(i))*PX5(i)/E5(i)
              !                    gy_now=GY5(i)+(t2time-FT5(i))*PY5(i)/E5(i)
              !                    gz_now=GZ5(i)+(t2time-FT5(i))*PZ5(i)/E5(i)
              !                    rt_now=dsqrt(gx_now**2+gy_now**2)
              !                    zmax=0.1d0*t2time
              !                    volume=3.1416d0*(2d0**2)*(2*zmax)
              !                    if(rt_now.gt.2d0.or.dabs(gz_now).gt.zmax)
              !     1                   goto 100
              !test-end
              iloop = 1
              If (dabs(ypartn)<=1D0) Then
                 iloop = 2
                 If (dabs(ypartn)<=0.5D0) Then
                    iloop = 3
                 End If
              End If
              Do iy = 1, iloop
                 !lin-5/2012:
                 !                       if(pt2.gt.0.) then
                 If (pt2>0D0) Then
                    v2prtn = (px5(i)**2-py5(i)**2)/pt2
                    !lin-5/2012:
                    !                          if(abs(v2prtn).gt.1.)
                    If (dabs(v2prtn)>1D0) Write (nfile(iy), *) 'v2prtn>1', v2prtn
                    v2p(ianp, iy) = v2p(ianp, iy) + v2prtn
                    v2p2(ianp, iy) = v2p2(ianp, iy) + v2prtn**2
                 End If
                 xperp2 = gx5(i)**2 + gy5(i)**2
                 !lin-5/2012:
                 !                       if(xperp2.gt.0.)
                 If (xperp2>0D0) s2p(ianp, iy) = s2p(ianp, iy) + (gx5(i)**2-gy5(i)**2)/xperp2
                 xnpart(ianp, iy) = xnpart(ianp, iy) + 1D0
                 etp(ianp, iy) = etp(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                 !test off density:
                 !                       etp(ianp,iy)=etp(ianp,iy)
                 !     1                  +dsqrt(pt2+XMASS5(I)**2+PZ5(i)**2)/volume
                 !lin-2/22/00 to write out parton info only for formed ones:
                 If (ft5(i)<=t2time) Then
                    v2pf(ianp, iy) = v2pf(ianp, iy) + v2prtn
                    v2pf2(ianp, iy) = v2pf2(ianp, iy) + v2prtn**2
                    !lin-5/2012:
                    !                          if(xperp2.gt.0.)
                    If (xperp2>0D0) s2pf(ianp, iy) = s2pf(ianp, iy) + (gx5(i)**2-gy5(i)**2)/xperp2
                    xnpf(ianp, iy) = xnpf(ianp, iy) + 1D0
                    etpf(ianp, iy) = etpf(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                    !test off density:
                    !                  etpf(ianp,iy)=etpf(ianp,iy)
                    !     1                   +dsqrt(pt2+XMASS5(I)**2+PZ5(i)**2)/volume
                 End If
              End Do
           End Do
           itimep = ianp
           iaevtp = iaevt
           !lin-3/30/00 ebe v2 variables:
           If (ianp==30) Then
              Do iy = 1, 3
                 npartn = idint(xnpart(ianp,iy)-xnpp(iy))
                 If (npartn/=0) Then
                    v2pevt(iy) = (v2p(ianp,iy)-v2pp(iy))/npartn
                    v2psum(iy) = v2psum(iy) + v2pevt(iy)
                    v2p2sm(iy) = v2p2sm(iy) + v2pevt(iy)**2
                    v2pp(iy) = v2p(ianp, iy)
                    xnpp(iy) = xnpart(ianp, iy)
                 End If
              End Do
              Write (49, 160) iaevt, v2pevt
           End If
           Goto 101
        End If
     End Do
     !lin-11/28/00 for animation:
     !lin-8/2015 ctest off turn off parton-t.dat:
     ! 101       if(nevent.ge.1) then
101  If (nevent<1) Then
        Do nt = 1, ntmax
           time1 = dble(nt*dt)
           time2 = dble((nt+1)*dt)
           If (t2time<time2 .And. t2time>=time1) Then
              If (nt<=itanim) Return
              If (ifanim==1) Write (11, *) t2time
              iform = 0
              ne1all = 0
              ne1form = 0
              Do i = 1, mul
                 !     Calculate parton coordinates after propagation to current time:
                 gz_now = gz5(i) + (t2time-ft5(i))*pz5(i)/e5(i)
                 If (dabs(gz_now)<t2time) Then
                    etap = 0.5D0*dlog((t2time+gz_now)/(t2time-gz_now))
                 Else
                    etap = 1000000.D0*sign(1.D0, gz_now)
                 End If
                 ne1all = ne1all + 1
                 If (ft5(i)<=t2time) ne1form = ne1form + 1
                 !     write out parton info only for formed ones for animation:
                 If (ft5(i)<=t2time) iform = iform + 1
              End Do
              !lin-8/2015 for animation:
              If (ifanim==1) Write (11, *) iform
              Write (93, 184) 'evt#,t,np,npformed=', iaevt, t2time, ne1all, ne1form
              !
              Do i = 1, mul
                 If (ft5(i)<=t2time) Then
                    !     propagate formed partons to current time t2time using parton v:
                    gz_now = gz5(i) + (t2time-ft5(i))*pz5(i)/e5(i)
                 Else
                    !     back-propagate unformed partons using parent hadron v:
                    gz_now = gz5(i) + (t2time-ft5(i))*vzp(i)
                 End If
                 !
                 If (dabs(gz_now)<t2time) Then
                    etap = 0.5D0*dlog((t2time+gz_now)/(t2time-gz_now))
                 Else
                    etap = 1000000.D0*sign(1.D0, gz_now)
                 End If
                 !     calculate other coordinates of the parton:
                 If (ft5(i)<=t2time) Then
                    gx_now = gx5(i) + (t2time-ft5(i))*px5(i)/e5(i)
                    gy_now = gy5(i) + (t2time-ft5(i))*py5(i)/e5(i)
                 Else
                    gx_now = gx5(i) + (t2time-ft5(i))*vxp(i)
                    gy_now = gy5(i) + (t2time-ft5(i))*vyp(i)
                 End If
                 Write (93, 185) ityp5(i), px5(i), py5(i), pz5(i), xmass5(i), gx_now, gy_now, ft5(i), etap
                 !lin-8/2015 for animation:
                 If (ifanim==1 .And. ft5(i)<=t2time) Then
                    Write (11, 180) ityp5(i), gx5(i), gy5(i), gz5(i), ft5(i), px5(i), py5(i), pz5(i), e5(i)
                 End If
              End Do
              itanim = nt
           End If
        End Do
     End If
     !lin-5/17/01 calculate v2 for active partons (which have not frozen out):
     !     idd=3, called at end of zpc.f:
  Else If (idd==3) Then
     Do ianp = 1, 30
        If (iscatt(ianp)==0) tscatt(ianp+1) = tscatt(ianp)
     End Do
     Do i = 1, mul
        !lin-8/2015 to avoid IEEE_DIVIDE_BY_ZERO or IEEE_INVALID:
        !              ypartn=0.5d0*dlog((E5(i)+PZ5(i))
        !     1             /(E5(i)-PZ5(i)+1.d-8))
        delta = 1D-8
        If ((e5(i)-dabs(pz5(i))+delta)<=0) Then
           Write (6, *) 'ypartn error', e5(i) - dabs(pz5(i))
           ypartn = 1000000.D0*sign(1.D0, pz5(i))
        Else
           ypartn = 0.5D0*dlog((e5(i)+pz5(i)+delta)/(e5(i)-pz5(i)+delta))
        End If
        pt2 = px5(i)**2 + py5(i)**2
        iloop = 1
        If (dabs(ypartn)<=1D0) Then
           iloop = 2
           If (dabs(ypartn)<=0.5D0) Then
              iloop = 3
           End If
        End If
        !
        Do ianp = 1, 30
           If (iscatt(ianp)/=0) Then
              If (ft5(i)<tscatt(ianp+1) .And. ft5(i)>=tscatt(ianp)) Then
                 Do iy = 1, iloop
                    !lin-5/2012:
                    !                          if(pt2.gt.0.) then
                    If (pt2>0D0) Then
                       v2prtn = (px5(i)**2-py5(i)**2)/pt2
                       v2pfrz(ianp, iy) = v2pfrz(ianp, iy) + v2prtn
                       v2p2fz(ianp, iy) = v2p2fz(ianp, iy) + v2prtn**2
                    End If
                    xperp2 = gx5(i)**2 + gy5(i)**2
                    !lin-5/2012:
                    !                          if(xperp2.gt.0.) s2pfrz(ianp,iy)=
                    If (xperp2>0D0) s2pfrz(ianp, iy) = s2pfrz(ianp, iy) + (gx5(i)**2-gy5(i)**2)/xperp2
                    etpfrz(ianp, iy) = etpfrz(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                    xnpfrz(ianp, iy) = xnpfrz(ianp, iy) + 1D0
                    !test off density:
                    !                    etpfrz(ianp,iy)=etpfrz(ianp,iy)
                    !     1                   +dsqrt(pt2+XMASS5(I)**2+PZ5(i)**2)/volume
                    !sp07/05
                    If (ityp5(i)==1 .Or. ityp5(i)==2) Then
                       etlfrz(ianp, iy) = etlfrz(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                       xnlfrz(ianp, iy) = xnlfrz(ianp, iy) + 1D0
                    Else If (ityp5(i)==3) Then
                       etsfrz(ianp, iy) = etsfrz(ianp, iy) + dsqrt(pt2+xmass5(i)**2)
                       xnsfrz(ianp, iy) = xnsfrz(ianp, iy) + 1D0
                    End If
                    !sp07/05 end
                 End Do
                 !     parton freezeout info taken, proceed to next parton:
                 Goto 350
              End If
           End If
        End Do
350  End Do
     !     idd=2: calculate average partonic elliptic flow, called from artdri.f,
  Else If (idd==2) Then
     Do i = 1, 3
        Write (nfile(i), *) '   tsp,   v2p,     v2p2, ' // '   s2p,  etp,   xmult,    nevt,  xnctot'
        Write ((nfile(i)+1), *) '  tsp,   v2pf,   v2pf2, ' // '   s2pf, etpf,  xnform,  xcrate'
        Write ((nfile(i)+2), *) '  tsp,   v2pa,   v2pa2, ' // '   s2pa, etpa,  xmult_ap,  xnform,   nevt'
        Write ((nfile(i)+3), *) '  tsph,  v2ph,   v2ph2, ' // '   s2ph, etph,  xmult_(ap/2+h),xmult_ap/2,nevt'
        !sp
        Write (nlfile(i), *) '   tsp,    v2,     s2,    etp,    xmul'
        Write (nsfile(i), *) '   tsp,    v2,     s2,    etp,    xmul'
        Write (nlfile(i)+1, *) '   tsp,    v2,     s2,    etp,    xmul'
        Write (nsfile(i)+1, *) '   tsp,    v2,     s2,    etp,    xmul'
        !
     End Do
     !lin-3/30/00 ensemble average & variance of v2 (over particles in all events):
     Do ii = 1, 30
        If (nevt(ii)==0) Goto 150
        Do iy = 1, 3
           !lin-5/2012:
           !                 if(xnpart(ii,iy).gt.1) then
           If (xnpart(ii,iy)>1D0) Then
              v2p(ii, iy) = v2p(ii, iy)/xnpart(ii, iy)
              v2p2(ii, iy) = dsqrt((v2p2(ii,iy)/xnpart(ii,iy)-v2p(ii,iy)**2)/(xnpart(ii,iy)-1))
              s2p(ii, iy) = s2p(ii, iy)/xnpart(ii, iy)
              ! xmult and etp are multiplicity and et for an averaged event:
              xmult = dble(xnpart(ii,iy)/dble(nevt(ii)))
              etp(ii, iy) = etp(ii, iy)/dble(nevt(ii))
              !sp
              etpl(ii, iy) = etpl(ii, iy)/dble(nevt(ii))
              etps(ii, iy) = etps(ii, iy)/dble(nevt(ii))
              !
              xnctot = 0D0
              Do inum = 1, ii
                 xnctot = xnctot + xncoll(inum)
              End Do
              If (nevt(1)/=0) xnctot = xnctot/nevt(1)
              Write (nfile(iy), 200) tsp(ii), v2p(ii, iy), v2p2(ii, iy), s2p(ii, iy), etp(ii, iy), xmult, nevt(ii), xnctot
           End If
           If (nevt(ii)/=0) xcrate = xncoll(ii)/(tsp(ii+1)-tsp(ii))/nevt(ii)
           !
           !lin-5/2012:
           !                 if(xnpf(ii,iy).gt.1) then
           If (xnpf(ii,iy)>1D0) Then
              v2pf(ii, iy) = v2pf(ii, iy)/xnpf(ii, iy)
              v2pf2(ii, iy) = dsqrt((v2pf2(ii,iy)/xnpf(ii,iy)-v2pf(ii,iy)**2)/(xnpf(ii,iy)-1))
              s2pf(ii, iy) = s2pf(ii, iy)/xnpf(ii, iy)
              xnform = dble(xnpf(ii,iy)/dble(nevt(ii)))
              etpf(ii, iy) = etpf(ii, iy)/dble(nevt(ii))
              !sp
              etplf(ii, iy) = etplf(ii, iy)/dble(nevt(ii))
              etpsf(ii, iy) = etpsf(ii, iy)/dble(nevt(ii))
              !
              Write (nfile(iy)+1, 210) tsp(ii), v2pf(ii, iy), v2pf2(ii, iy), s2pf(ii, iy), etpf(ii, iy), xnform, xcrate
           End If
           !sp
           !lin-5/2012:
           !                 if(xnpl(ii,iy).gt.1) then
           If (xnpl(ii,iy)>1D0) Then
              v2pl(ii, iy) = v2pl(ii, iy)/xnpl(ii, iy)
              s2pl(ii, iy) = s2pl(ii, iy)/xnpl(ii, iy)
              xmult = dble(xnpl(ii,iy)/dble(nevt(ii)))
              etpl(ii, iy) = etpl(ii, iy)/dble(nevt(ii))
              Write (nlfile(iy), 201) tsp(ii), v2pl(ii, iy), s2pl(ii, iy), etpl(ii, iy), xmult
           End If
           !lin-5/2012:
           !                 if(xnps(ii,iy).gt.1) then
           If (xnps(ii,iy)>1D0) Then
              v2ps(ii, iy) = v2ps(ii, iy)/xnps(ii, iy)
              s2ps(ii, iy) = s2ps(ii, iy)/xnps(ii, iy)
              xmult = dble(xnps(ii,iy)/dble(nevt(ii)))
              etps(ii, iy) = etps(ii, iy)/dble(nevt(ii))
              Write (nsfile(iy), 201) tsp(ii), v2ps(ii, iy), s2ps(ii, iy), etps(ii, iy), xmult
           End If
           !lin-5/2012:
           !                 if(xnplf(ii,iy).gt.1) then
           If (xnplf(ii,iy)>1D0) Then
              v2plf(ii, iy) = v2plf(ii, iy)/xnplf(ii, iy)
              s2plf(ii, iy) = s2plf(ii, iy)/xnplf(ii, iy)
              xmult = dble(xnplf(ii,iy)/dble(nevt(ii)))
              etplf(ii, iy) = etplf(ii, iy)/dble(nevt(ii))
              Write (nlfile(iy)+1, 201) tsp(ii), v2plf(ii, iy), s2plf(ii, iy), etplf(ii, iy), xmult
           End If
           !lin-5/2012:
           !                 if(xnpsf(ii,iy).gt.1) then
           If (xnpsf(ii,iy)>1D0) Then
              v2psf(ii, iy) = v2psf(ii, iy)/xnpsf(ii, iy)
              s2psf(ii, iy) = s2psf(ii, iy)/xnpsf(ii, iy)
              xmult = dble(xnpsf(ii,iy)/dble(nevt(ii)))
              etpsf(ii, iy) = etpsf(ii, iy)/dble(nevt(ii))
              Write (nsfile(iy)+1, 201) tsp(ii), v2psf(ii, iy), s2psf(ii, iy), etpsf(ii, iy), xmult
           End If
           !sp-end
        End Do
150  End Do
     !sp07/05 initial & final mt distrb
     scalei = 0D0
     scalef = 0D0
     If (nevt(1)/=0) scalei = 1D0/dble(nevt(1))/bmt
     If (nevt(30)/=0) scalef = 1D0/dble(nevt(30))/bmt
     Do iy = 2, 3
        yra = 1D0
        If (iy==2) yra = 2D0
        Do i = 1, 50
           Write (nmfile(iy), 251) bmt*dble(i-0.5), scalei*dmyil(i, iy)/yra, scalef*dmyfl(i, iy)/yra, scalei*dmyis(i, iy)/yra, scalef*dmyfs(i, iy)/yra
        End Do
     End Do
     !sp07/05 end
     !lin-3/30/00 event-by-event average & variance of v2:
     If (nevt(30)>=1) Then
        Do iy = 1, 3
           v2pavg(iy) = v2psum(iy)/nevt(30)
           v2var0 = v2p2sm(iy)/nevt(30) - v2pavg(iy)**2
           !lin-5/2012:
           !                 if(v2var0.gt.0) varv2p(iy)=dsqrt(v2var0)
           If (v2var0>0D0) varv2p(iy) = dsqrt(v2var0)
        End Do
        Write (49, 240) 'EBE v2p,v2p(y2),v2p(y1): avg=', v2pavg
        Write (49, 240) 'EBE v2p,v2p(y2),v2p(y1): var=', varv2p
     End If
     !lin-8/2015:
     !lin-11/28/00 for animation:
     If (ifanim==1) Then
        Do i = 1, mul
           If (ft5(i)<=t2time) Then
              Write (15, 140) ityp5(i), gx5(i), gy5(i), gz5(i), ft5(i)
           End If
        End Do
        !lin-11/29/00 signal the end of animation file:
        Write (10, *) - 10.
        Write (10, *) 0
        Write (11, *) - 10.
        Write (11, *) 0
        Close (10)
        Close (11)
        Close (15)
     End If

     !lin-5/18/01 calculate v2 for active partons:
     Do ianp = 1, 30
        Do iy = 1, 3
           v2pact = 0D0
           v2p2ac = 0D0
           s2pact = 0D0
           etpact = 0D0
           xnacti = 0D0
           !lin-5/2012:
           !                 if(xnpf(ianp,iy).gt.1) then
           If (xnpf(ianp,iy)>1D0) Then
              !     reconstruct the sum of v2p, v2p2, s2p, etp, and xnp for formed partons:
              v2pact = v2pf(ianp, iy)*xnpf(ianp, iy)
              v2p2ac = (v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)+v2pf(ianp,iy)**2)*xnpf(ianp, iy)
              s2pact = s2pf(ianp, iy)*xnpf(ianp, iy)
              etpact = etpf(ianp, iy)*dble(nevt(ianp))
              xnpact = xnpf(ianp, iy)
              !
              Do kanp = 1, ianp
                 v2pact = v2pact - v2pfrz(kanp, iy)
                 v2p2ac = v2p2ac - v2p2fz(kanp, iy)
                 s2pact = s2pact - s2pfrz(kanp, iy)
                 etpact = etpact - etpfrz(kanp, iy)
                 xnpact = xnpact - xnpfrz(kanp, iy)
              End Do
              !     save the sum of v2p, v2p2, s2p, etp, and xnp for formed partons:
              v2ph = v2pact
              v2ph2 = v2p2ac
              s2ph = s2pact
              etph = etpact
              xnp2 = xnpact/2D0
              !
              !lin-5/2012:
              !                    if(xnpact.gt.1.and.nevt(ianp).ne.0) then
              If (xnpact>1D0 .And. nevt(ianp)/=0) Then
                 v2pact = v2pact/xnpact
                 v2p2ac = dsqrt((v2p2ac/xnpact-v2pact**2)/(xnpact-1))
                 s2pact = s2pact/xnpact
                 xnacti = dble(xnpact/dble(nevt(ianp)))
                 etpact = etpact/dble(nevt(ianp))
                 Write (nfile(iy)+2, 250) tsp(ianp), v2pact, v2p2ac, s2pact, etpact, xnacti, xnpf(ianp, iy)/dble(nevt(ianp)), nevt(ianp)
              End If
           End If
           !     To calculate combined v2 for active partons plus formed hadrons,
           !     add the sum of v2h, v2h2, s2h, eth, and xnh for formed hadrons:
           !     scale the hadron part in case nevt(ianp) != nevent:
           shadr = dble(nevt(ianp))/dble(nevent)
           ianh = ianp
           v2ph = v2ph + v2h(ianh, iy)*xnhadr(ianh, iy)*shadr
           v2ph2 = v2ph2 + (v2h2(ianh,iy)**2*(xnhadr(ianh,iy)-1)+v2h(ianh,iy)**2)*xnhadr(ianh, iy)*shadr
           s2ph = s2ph + s2h(ianh, iy)*xnhadr(ianh, iy)*shadr
           etph = etph + eth(ianh, iy)*dble(nevent)*shadr
           xnph = xnpact + xnhadr(ianh, iy)*shadr
           xnp2h = xnp2 + xnhadr(ianh, iy)*shadr
           !lin-8/2015 to avoid IEEE_DIVIDE_BY_ZERO:
           !clin-5/2012:
           !c                 if(xnph.gt.1) then
           !                 if(xnph.gt.1d0) then
           If (xnph>1D0 .And. nevt(ianp)/=0) Then
              v2ph = v2ph/xnph
              v2ph2 = dsqrt((v2ph2/xnph-v2ph**2)/(xnph-1))
              s2ph = s2ph/xnph
              etph = etph/dble(nevt(ianp))
              xnp2 = xnp2/dble(nevt(ianp))
              xnp2h = xnp2h/dble(nevent)
              !lin-8/2015
              !                    if(tsp(ianp).le.(ntmax*dt))
              If (tsp(ianp)<=dble(ntmax*dt)) Write (nfile(iy)+3, 250) tsp(ianp), v2ph, v2ph2, s2ph, etph, xnp2h, xnp2, nevt(ianp)
           End If
           !
        End Do
     End Do
     Do ianp = 1, 30
        Do iy = 1, 3
           v2pact = 0D0
           v2p2ac = 0D0
           s2pact = 0D0
           etpact = 0D0
           xnacti = 0D0
           !     reconstruct the sum of v2p, v2p2, s2p, etp, and xnp for formed partons:
           v2pact = v2pf(ianp, iy)*xnpf(ianp, iy)
           v2p2ac = (v2pf2(ianp,iy)**2*(xnpf(ianp,iy)-1)+v2pf(ianp,iy)**2)*xnpf(ianp, iy)
           s2pact = s2pf(ianp, iy)*xnpf(ianp, iy)
           etpact = etpf(ianp, iy)*dble(nevt(ianp))
           xnpact = xnpf(ianp, iy)
        End Do
     End Do
     Close (620)
     Close (630)
     Do nf = 1, 3
        Do ifile = 0, 3
           Close (nfile(nf)+ifile)
        End Do
     End Do
     Do nf = 1, 3
        Close (740+nf)
     End Do
  End If
  !
  Return
184 Format (A20, I7, F8.4, 2(1X,I6))
185 Format (I3, 3(1X,F8.3), 1X, F8.4, 1X, 2(F8.3,1X), F11.4, 1X, F8.3)
  !
160 Format (I10, 3(2X,F9.5))
180 Format (I6, 8(1X,F7.2))
140 Format (I10, 4(2X,F7.2))
200 Format (2X, F5.2, 3(2X,F7.4), 2(2X,F9.2), I6, 2X, F9.2)
210 Format (2X, F5.2, 3(2X,F7.4), 3(2X,F9.2))
240 Format (A30, 3(2X,F9.5))
250 Format (2X, F5.2, 3(2X,F7.4), 3(2X,F9.2), I6)
  !sp
201 Format (2X, F5.2, 4(2X,F9.2))
251 Format (5E15.5)
End Subroutine flowp

!=======================================================================
!     Calculate flow from formed hadrons, called by art1e.f:
!     Note: numbers in art not in double precision!
Subroutine flowh(ct)
  Parameter (maxstr=150001, maxr=1)
  Dimension tsh(31)
  Double Precision v2h, xnhadr, eth, v2h2, s2h
  Double Precision v2hp, xnhadp, v2hsum, v2h2sm, v2hevt(3)
  Double Precision pt2, v2hadr
  Common /hflow/v2h(30, 3), xnhadr(30, 3), eth(30, 3), v2h2(30, 3), s2h(30, 3)
  !c      SAVE /hflow/
  Common /ebe/v2hp(3), xnhadp(3), v2hsum(3), v2h2sm(3)
  !c      SAVE /ebe/
  Common /lastt/itimeh, bimp
  !c      SAVE /lastt/
  Common /run/num
  !c      SAVE /RUN/
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
  Common /anim/nevent, isoft, isflag, izpc
  !c      SAVE /anim/
  Common /arevt/iaevt, iarun, miss
  !c      SAVE /AREVT/
  Save
  !
  Do ii = 1, 31
     tsh(ii) = float(ii-1)
  End Do
  !
  Do ianh = 1, 30
     If ((ct+0.0001)<tsh(ianh+1) .And. (ct+0.0001)>=tsh(ianh)) Then
        If (ianh==itimeh) Goto 101
        ia = 0
        Do j = 1, num
           mult = massr(j)
           ia = ia + massr(j-1)
           Do ic = 1, mult
              i = ia + ic
              !     5/04/01 exclude leptons and photons:
              If (iabs(lb(i)-10000)<100) Goto 100
              px = p(1, i)
              py = p(2, i)
              pt2 = dble(px)**2 + dble(py)**2
              !     2/18/00 Note: e(i) gives the mass in ART:
              ene = sqrt(e(i)**2+sngl(pt2)+p(3,i)**2)
              rap = 0.5*alog((ene+p(3,i))/(ene-p(3,i)))
              !test off density with 2fm radius and z:(-0.1*t,0.1*t):
              !                rt_now=sqrt(r(1,i)**2+r(2,i)**2)
              !                gz_now=r(3,i)
              !                zmax=0.1*ct
              !                volume=3.1416*(2.**2)*2*zmax
              !                if(rt_now.gt.2.or.abs(gz_now).gt.zmax)
              !     1               goto 100
              iloop = 1
              If (abs(rap)<=1) Then
                 iloop = 2
                 If (abs(rap)<=0.5) Then
                    iloop = 3
                 End If
              End If
              Do iy = 1, iloop
                 If (pt2>0D0) Then
                    v2hadr = (dble(px)**2-dble(py)**2)/pt2
                    v2h(ianh, iy) = v2h(ianh, iy) + v2hadr
                    v2h2(ianh, iy) = v2h2(ianh, iy) + v2hadr**2
                    If (dabs(v2hadr)>1D0) Write (1, *) 'v2hadr>1', v2hadr, px, py
                 End If
                 xperp2 = r(1, i)**2 + r(2, i)**2
                 If (xperp2>0.) s2h(ianh, iy) = s2h(ianh, iy) + dble((r(1,i)**2-r(2,i)**2)/xperp2)
                 eth(ianh, iy) = eth(ianh, iy) + dble(sqrt(e(i)**2+sngl(pt2)))
                 !test off density:
                 !               eth(ianh,iy)=eth(ianh,iy)
                 !     1                  +dble(SQRT(e(i)**2+sngl(pt2)+p(3,i)**2))/volume
                 xnhadr(ianh, iy) = xnhadr(ianh, iy) + 1D0
              End Do
100        End Do
        End Do
        itimeh = ianh
        !lin-5/04/01 ebe v2 variables:
        If (ianh==30) Then
           Do iy = 1, 3
              nhadrn = idint(xnhadr(ianh,iy)-xnhadp(iy))
              If (nhadrn/=0) Then
                 v2hevt(iy) = (v2h(ianh,iy)-v2hp(iy))/dble(nhadrn)
                 v2hsum(iy) = v2hsum(iy) + v2hevt(iy)
                 v2h2sm(iy) = v2h2sm(iy) + v2hevt(iy)**2
                 v2hp(iy) = v2h(ianh, iy)
                 xnhadp(iy) = xnhadr(ianh, iy)
              End If
           End Do
           Write (88, 160) iaevt, v2hevt
        End If
        Goto 101
     End If
  End Do
  !lin-8/2015:
  !lin-11/27/00 for animation:
  !test off     turn off animation output (0 to turn off and 1 to turn on):
101 ifanim = 0
  If (ifanim==1) Then
     ia = 0
     Do j = 1, num
        mult = massr(j)
        ia = ia + massr(j-1)
        Write (10, *) ct
        Write (10, *) mult
        Do ic = 1, mult
           i = ia + ic
           !lin-6/2013 for animation:
           If (amax1(abs(r(1,i)),abs(r(2,i)),abs(r(3,i)))<9999) Then
              Write (10, 210) lb(i), r(1, i), r(2, i), r(3, i), p(1, i), p(2, i), p(3, i), e(i)
           Else
              Write (10, 220) lb(i), r(1, i), r(2, i), r(3, i), p(1, i), p(2, i), p(3, i), e(i)
           End If
        End Do
     End Do
     Return
  End If
  Return
160 Format (I10, 3(2X,F9.5))
210 Format (I6, 7(1X,F9.3))
220 Format (I6, 3(1X,E9.3), 4(1X,F9.3))
End Subroutine flowh

!=======================================================================
Subroutine flowh0(nevnt, idd)
  !
  Dimension tsh(31)
  Double Precision v2h, xnhadr, eth, v2h2, s2h
  Double Precision v2hp, xnhadp, v2hsum, v2h2sm, v2havg(3), varv2h(3)
  Common /hflow/v2h(30, 3), xnhadr(30, 3), eth(30, 3), v2h2(30, 3), s2h(30, 3)
  !c      SAVE /hflow/
  Common /ebe/v2hp(3), xnhadp(3), v2hsum(3), v2h2sm(3)
  !c      SAVE /ebe/
  Common /input1/masspr, massta, iseed, iavoid, dt
  !c      SAVE /input1/
  Common /input2/ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, icflow, icrho, icou, kpoten, kmul
  !c      SAVE /INPUT2/
  Common /lastt/itimeh, bimp
  !c      SAVE /lastt/
  Save

  !     idd=0: initialization for flow analysis, called by artdri.f::
  If (idd==0) Then
     itimeh = 0
     !
     Do ii = 1, 31
        tsh(ii) = float(ii-1)
     End Do
     !
     Do ii = 1, 30
        Do iy = 1, 3
           v2h(ii, iy) = 0D0
           xnhadr(ii, iy) = 0D0
           eth(ii, iy) = 0D0
           v2h2(ii, iy) = 0D0
           s2h(ii, iy) = 0D0
        End Do
     End Do
     Do iy = 1, 3
        v2hp(iy) = 0D0
        xnhadp(iy) = 0D0
        v2hsum(iy) = 0D0
        v2h2sm(iy) = 0D0
        If (iy==1) Then
           nunit = 59
        Else If (iy==2) Then
           nunit = 68
        Else
           nunit = 69
        End If
        Write (nunit, *) '   tsh,   v2h,     v2h2,     s2h, ' // ' eth,   xmulth'
     End Do
     !     idd=2: calculate average hadronic elliptic flow, called by artdri.f:
  Else If (idd==2) Then
     Do ii = 1, 30
        Do iy = 1, 3
           If (xnhadr(ii,iy)==0) Then
              xmulth = 0.
           Else If (xnhadr(ii,iy)>1) Then
              v2h(ii, iy) = v2h(ii, iy)/xnhadr(ii, iy)
              eth(ii, iy) = eth(ii, iy)/dble(nevnt)
              v2h2(ii, iy) = dsqrt((v2h2(ii,iy)/xnhadr(ii,iy)-v2h(ii,iy)**2)/(xnhadr(ii,iy)-1))
              s2h(ii, iy) = s2h(ii, iy)/xnhadr(ii, iy)
              xmulth = sngl(xnhadr(ii,iy)/nevnt)
           End If
           If (iy==1) Then
              nunit = 59
           Else If (iy==2) Then
              nunit = 68
           Else
              nunit = 69
           End If
           If (tsh(ii)<=(ntmax*dt)) Write (nunit, 200) tsh(ii), v2h(ii, iy), v2h2(ii, iy), s2h(ii, iy), eth(ii, iy), xmulth
        End Do
     End Do
     !     event-by-event average & variance of v2h:
     Do iy = 1, 3
        v2havg(iy) = v2hsum(iy)/dble(nevnt)
        varv2h(iy) = dsqrt(v2h2sm(iy)/dble(nevnt)-v2havg(iy)**2)
     End Do
     Write (88, 240) 'EBE v2h,v2h(y2),v2h(y1): avg=', v2havg
     Write (88, 240) 'EBE v2h,v2h(y2),v2h(y1): var=', varv2h
  End If
  Return
200 Format (2X, F5.2, 3(2X,F7.4), 2(2X,F9.2))
240 Format (A30, 3(2X,F9.5))
End Subroutine flowh0

!=======================================================================
!     2/23/00 flow from all initial hadrons just before entering ARTMN:
Subroutine iniflw(nevnt, idd)
  Parameter (maxstr=150001, maxr=1)
  Double Precision v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
  Common /run/num
  !c      SAVE /RUN/
  Common /arerc1/multi1(maxr)
  !c      SAVE /ARERC1/
  Common /arprc1/ityp1(maxstr, maxr), gx1(maxstr, maxr), gy1(maxstr, maxr), gz1(maxstr, maxr), ft1(maxstr, maxr), px1(maxstr, maxr), py1(maxstr, maxr), pz1(maxstr, maxr), ee1(maxstr, maxr), xm1(maxstr, maxr)
  !c      SAVE /ARPRC1/
  Common /iflow/v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult
  !c      SAVE /iflow/
  Save
  !
  If (idd==0) Then
     v2i = 0D0
     eti = 0D0
     xmulti = 0D0
     v2mi = 0D0
     s2mi = 0D0
     xmmult = 0D0
     v2bi = 0D0
     s2bi = 0D0
     xbmult = 0D0
  Else If (idd==1) Then
     Do j = 1, num
        Do i = 1, multi1(j)
           ityp = ityp1(i, j)
           !     all hadrons:
           If (ityp>-100 .And. ityp<100) Goto 100
           xmulti = xmulti + 1.D0
           px = px1(i, j)
           py = py1(i, j)
           xm = xm1(i, j)
           pt2 = px**2 + py**2
           xh = gx1(i, j)
           yh = gy1(i, j)
           xt2 = xh**2 + yh**2
           If (pt2>0) v2i = v2i + dble((px**2-py**2)/pt2)
           eti = eti + dble(sqrt(px**2+py**2+xm**2))
           !     baryons only:
           If (ityp<-1000 .Or. ityp>1000) Then
              xbmult = xbmult + 1.D0
              If (pt2>0) v2bi = v2bi + dble((px**2-py**2)/pt2)
              If (xt2>0) s2bi = s2bi + dble((xh**2-yh**2)/xt2)
              !     mesons only:
           Else
              xmmult = xmmult + 1.D0
              If (pt2>0) v2mi = v2mi + dble((px**2-py**2)/pt2)
              If (xt2>0) s2mi = s2mi + dble((xh**2-yh**2)/xt2)
           End If
100        Continue
        End Do
     End Do
  Else If (idd==2) Then
     If (xmulti/=0) v2i = v2i/xmulti
     eti = eti/dble(nevnt)
     xmulti = xmulti/dble(nevnt)
     If (xmmult/=0) Then
        v2mi = v2mi/xmmult
        s2mi = s2mi/xmmult
     End If
     xmmult = xmmult/dble(nevnt)
     If (xbmult/=0) Then
        v2bi = v2bi/xbmult
        s2bi = s2bi/xbmult
     End If
     xbmult = xbmult/dble(nevnt)
  End If
  !
  Return
End Subroutine iniflw

!=======================================================================
!     2/25/00 dN/dt analysis for production (before ZPCMN)
!     and freezeout (right after ZPCMN) for all partons.
Subroutine frztm(nevnt, idd)
  !
  Implicit Double Precision (A-H, O-Z)
  Parameter (maxptn=400001)
  Dimension tsf(31)
  Common /para1/mul
  !c      SAVE /PARA1/
  Common /prec1/gx0(maxptn), gy0(maxptn), gz0(maxptn), ft0(maxptn), px0(maxptn), py0(maxptn), pz0(maxptn), e0(maxptn), xmass0(maxptn), ityp0(maxptn)
  !c      SAVE /prec1/
  Common /prec2/gx5(maxptn), gy5(maxptn), gz5(maxptn), ft5(maxptn), px5(maxptn), py5(maxptn), pz5(maxptn), e5(maxptn), xmass5(maxptn), ityp5(maxptn)
  !c      SAVE /prec2/
  Common /frzout/xnprod(30), etprod(30), xnfrz(30), etfrz(30), dnprod(30), detpro(30), dnfrz(30), detfrz(30)
  !c      SAVE /frzout/
  Save
  Data tsf/0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30/
  !
  If (idd==0) Then
     Do ii = 1, 30
        xnprod(ii) = 0D0
        etprod(ii) = 0D0
        xnfrz(ii) = 0D0
        etfrz(ii) = 0D0
        dnprod(ii) = 0D0
        detpro(ii) = 0D0
        dnfrz(ii) = 0D0
        detfrz(ii) = 0D0
     End Do
     Open (86, File='ana1/production.dat', Status='UNKNOWN')
     Open (87, File='ana1/freezeout.dat', Status='UNKNOWN')
  Else If (idd==1) Then
     Do ip = 1, mul
        Do ii = 1, 30
           eth0 = dsqrt(px0(ip)**2+py0(ip)**2+xmass0(ip)**2)
           eth2 = dsqrt(px5(ip)**2+py5(ip)**2+xmass5(ip)**2)
           !     total number and Et produced by time tsf(ii):
           If (ft0(ip)<tsf(ii+1)) Then
              xnprod(ii) = xnprod(ii) + 1D0
              etprod(ii) = etprod(ii) + eth0
              !     number and Et produced from time tsf(ii) to tsf(ii+1):
              If (ft0(ip)>=tsf(ii)) Then
                 dnprod(ii) = dnprod(ii) + 1D0
                 detpro(ii) = detpro(ii) + eth0
              End If
           End If
           !     total number and Et freezed out by time tsf(ii):
           If (ft5(ip)<tsf(ii+1)) Then
              xnfrz(ii) = xnfrz(ii) + 1D0
              etfrz(ii) = etfrz(ii) + eth2
              !     number and Et freezed out from time tsf(ii) to tsf(ii+1):
              If (ft5(ip)>=tsf(ii)) Then
                 dnfrz(ii) = dnfrz(ii) + 1D0
                 detfrz(ii) = detfrz(ii) + eth2
              End If
           End If
        End Do
     End Do
  Else If (idd==2) Then
     Write (86, *) '       t,       np,       dnp/dt,      etp ' // ' detp/dt'
     Write (87, *) '       t,       nf,       dnf/dt,      etf ' // ' detf/dt'
     Do ii = 1, 30
        xnp = xnprod(ii)/dble(nevnt)
        xnf = xnfrz(ii)/dble(nevnt)
        etp = etprod(ii)/dble(nevnt)
        etf = etfrz(ii)/dble(nevnt)
        dxnp = dnprod(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        dxnf = dnfrz(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        detp = detpro(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        detf = detfrz(ii)/dble(nevnt)/(tsf(ii+1)-tsf(ii))
        Write (86, 200) tsf(ii+1), xnp, dxnp, etp, detp
        Write (87, 200) tsf(ii+1), xnf, dxnf, etf, detf
     End Do
  End If
  !
  Return
200 Format (2X, F9.2, 4(2X,F10.2))
End Subroutine frztm

!=======================================================================
!lin-6/2009 write out initial minijet information
!     before propagating to its formation time:
!lin-2/2012:
!        subroutine minijet_out(BB)
Subroutine minijet_out(bb, phirp)
  Parameter (maxstr=150001)
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /arevt/iaevt, iarun, miss
  Common /para7/ioscar, nsmbbbar, nsmmeson
  Common /phidcy/iphidcy, pttrig, ntrig, maxmiss, ipi0dcy
  Save
  ntrig = 0
  Do i = 1, ihnt2(1)
     Do j = 1, npj(i)
        pt = sqrt(pjpx(i,j)**2+pjpy(i,j)**2)
        If (pt>=pttrig) ntrig = ntrig + 1
     End Do
  End Do
  Do i = 1, ihnt2(3)
     Do j = 1, ntj(i)
        pt = sqrt(pjtx(i,j)**2+pjty(i,j)**2)
        If (pt>=pttrig) ntrig = ntrig + 1
     End Do
  End Do
  Do i = 1, nsg
     Do j = 1, njsg(i)
        pt = sqrt(pxsg(i,j)**2+pysg(i,j)**2)
        If (pt>=pttrig) ntrig = ntrig + 1
     End Do
  End Do
  !     Require at least 1 initial minijet parton above the trigger Pt value:
  If (ntrig==0) Return

  !.....transfer data from HIJING to ZPC
  If (ioscar==3) Write (96, *) iaevt, miss, ihnt2(1), ihnt2(3)
  Do i = 1, ihnt2(1)
     Do j = 1, npj(i)
        ityp = kfpj(i, j)
        !     write out not only gluons:
        !              if(ityp.ne.21) goto 1007
        !lin-2/2012:
        !              gx=YP(1,I)+0.5*BB
        !              gy=YP(2,I)
        gx = yp(1, i) + 0.5*bb*cos(phirp)
        gy = yp(2, i) + 0.5*bb*sin(phirp)
        gz = 0.
        ft = 0.
        px = pjpx(i, j)
        py = pjpy(i, j)
        pz = pjpz(i, j)
        xmass = pjpm(i, j)
        If (ioscar==3) Then
           If (amax1(abs(gx),abs(gy),abs(gz),abs(ft))<9999) Then
              Write (96, 200) ityp, px, py, pz, xmass, gx, gy, gz, ft, 1
           Else
              Write (96, 201) ityp, px, py, pz, xmass, gx, gy, gz, ft, 1
           End If
        End If
     End Do
  End Do
  Do i = 1, ihnt2(3)
     Do j = 1, ntj(i)
        ityp = kftj(i, j)
        !              if(ityp.ne.21) goto 1009
        !lin-2/2012:
        !              gx=YT(1,I)-0.5*BB
        !              gy=YT(2,I)
        gx = yt(1, i) - 0.5*bb*cos(phirp)
        gy = yt(2, i) - 0.5*bb*sin(phirp)
        gz = 0.
        ft = 0.
        px = pjtx(i, j)
        py = pjty(i, j)
        pz = pjtz(i, j)
        xmass = pjtm(i, j)
        If (ioscar==3) Then
           If (amax1(abs(gx),abs(gy),abs(gz),abs(ft))<9999) Then
              Write (96, 200) ityp, px, py, pz, xmass, gx, gy, gz, ft, 2
           Else
              Write (96, 201) ityp, px, py, pz, xmass, gx, gy, gz, ft, 2
           End If
        End If
     End Do
  End Do
  Do i = 1, nsg
     Do j = 1, njsg(i)
        ityp = k2sg(i, j)
        !              if(ityp.ne.21) goto 1011
        gx = 0.5*(yp(1,iasg(i,1))+yt(1,iasg(i,2)))
        gy = 0.5*(yp(2,iasg(i,1))+yt(2,iasg(i,2)))
        gz = 0.
        ft = 0.
        px = pxsg(i, j)
        py = pysg(i, j)
        pz = pzsg(i, j)
        xmass = pmsg(i, j)
        If (ioscar==3) Then
           If (amax1(abs(gx),abs(gy),abs(gz),abs(ft))<9999) Then
              Write (96, 200) ityp, px, py, pz, xmass, gx, gy, gz, ft, 3
           Else
              Write (96, 201) ityp, px, py, pz, xmass, gx, gy, gz, ft, 3
           End If
        End If
     End Do
  End Do
  !
  Return
200 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 2(1X,F8.2), 2(2X,F2.0), 2X, I2)
201 Format (I6, 2(1X,F8.3), 1X, F10.3, 1X, F6.3, 2(1X,E8.2), 2(2X,F2.0), 2X, I2)
End Subroutine minijet_out

!=======================================================================
!lin-6/2009 embed back-to-back high-Pt quark/antiquark pair
!     via embedding back-to-back high-Pt pion pair then melting the pion pair
!     by generating the internal quark and antiquark momentum parallel to
!      the pion momentum (in order to produce a high-Pt and a low Pt parton):
Subroutine embedhighpt
  Parameter (maxstr=150001, maxr=1, pichmass=0.140, pi0mass=0.135, pi=3.1415926, nxymax=10001)
  Common /embed/iembed, nsembd, pxqembd, pyqembd, xembd, yembd, psembd, tmaxembd, phidecomp
  Common /rndf77/nseed
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  Common /arprc/itypar(maxstr), gxar(maxstr), gyar(maxstr), gzar(maxstr), ftar(maxstr), pxar(maxstr), pyar(maxstr), pzar(maxstr), pear(maxstr), xmar(maxstr)
  Common /anim/nevent, isoft, isflag, izpc
  Common /arevt/iaevt, iarun, miss
  Common /xyembed/nxyjet, xyjet(nxymax, 2)
  Save
  !
  If (iembed==1 .Or. iembed==2) Then
     xjet = xembd
     yjet = yembd
  Else If (iembed==3 .Or. iembed==4) Then
     If (nevent<=nxyjet) Then
        Read (97, *) xjet, yjet
     Else
        ixy = mod(iaevt, nxyjet)
        If (ixy==0) ixy = nxyjet
        xjet = xyjet(ixy, 1)
        yjet = xyjet(ixy, 2)
     End If
  Else
     Return
  End If
  !
  ptq = sqrt(pxqembd**2+pyqembd**2)
  If (ptq<(pichmass/2.)) Then
     Print *, 'Embedded quark transverse momentum is too small'
     Stop
  End If
  !     Randomly embed u/ubar or d/dbar at high Pt:
  idqembd = 1 + int(2*ranart(nseed))
  !     Flavor content for the charged pion that contains the leading quark:
  If (idqembd==1) Then
     idqsoft = -2
     idpi1 = -211
  Else If (idqembd==2) Then
     idqsoft = -1
     idpi1 = 211
  Else
     Print *, 'Wrong quark flavor embedded'
     Stop
  End If
  !     Caculate transverse momentum of the parent charged pion:
  xmq = ulmass(idqembd)
  xmqsoft = ulmass(idqsoft)
  ptpi = ((pichmass**2+xmq**2-xmqsoft**2)*ptq-sqrt((xmq**2+ptq**2)*(pichmass**4-2.*pichmass**2*(xmq**2+xmqsoft**2)+(xmq**2-xmqsoft**2)**2)))/(2.*xmq**2)
  If (iembed==1 .Or. iembed==3) Then
     pxpi1 = ptpi*pxqembd/ptq
     pypi1 = ptpi*pyqembd/ptq
     phidecomp = acos(pxqembd/ptq)
     If (pyqembd<0) phidecomp = 2.*pi - phidecomp
  Else
     phidecomp = 2.*pi*ranart(nseed)
     pxpi1 = ptpi*cos(phidecomp)
     pypi1 = ptpi*sin(phidecomp)
  End If
  !     Embedded quark/antiquark are assumed to have pz=0:
  pzpi1 = 0.
  !     Insert the two parent charged pions,
  !     ipion=1 for the pion containing the leading quark,
  !     ipion=2 for the pion containing the leading antiquark of the same flavor:
  Do ipion = 1, 2
     If (ipion==1) Then
        idpi = idpi1
        pxpi = pxpi1
        pypi = pypi1
        pzpi = pzpi1
     Else If (ipion==2) Then
        idpi = -idpi1
        pxpi = -pxpi1
        pypi = -pypi1
        pzpi = -pzpi1
     End If
     natt = natt + 1
     katt(natt, 1) = idpi
     katt(natt, 2) = 40
     katt(natt, 3) = 0
     patt(natt, 1) = pxpi
     patt(natt, 2) = pypi
     patt(natt, 3) = pzpi
     patt(natt, 4) = sqrt(pxpi**2+pypi**2+pzpi**2+pichmass**2)
     eatt = eatt + patt(natt, 4)
     gxar(natt) = xjet
     gyar(natt) = yjet
     gzar(natt) = 0.
     ftar(natt) = 0.
     itypar(natt) = katt(natt, 1)
     pxar(natt) = patt(natt, 1)
     pyar(natt) = patt(natt, 2)
     pzar(natt) = patt(natt, 3)
     pear(natt) = patt(natt, 4)
     xmar(natt) = pichmass
  End Do
  !
  !lin-8/2009
  !     Randomly embed a number of soft pions around each high-Pt quark in pair:
  If (nsembd>0) Then
     Do ipion = 1, 2
        Do ispion = 1, nsembd
           idsart = 3 + int(3*ranart(nseed))
           If (idsart==3) Then
              pimass = pichmass
              idpis = -211
           Else If (idsart==4) Then
              pimass = pi0mass
              idpis = 111
           Else
              pimass = pichmass
              idpis = 211
           End If
           natt = natt + 1
           katt(natt, 1) = idpis
           katt(natt, 2) = 40
           katt(natt, 3) = 0
           !     theta: relative angle between soft pion & associated high-Pt q or qbar,
           !     generate theta and phi uniformly:
           !     Note: it is not generated uniformly in solid angle because that gives
           !     a valley at theta=0, unlike the jet-like correlation (a peak at theta=0).
           theta = tmaxembd*ranart(nseed)
           phi = 2.*pi*ranart(nseed)
           pxspi = psembd*sin(theta)*cos(phi)
           pyspi = psembd*sin(theta)*sin(phi)
           pzspi = psembd*cos(theta)
           If (ipion==1) Then
              Call rotate(pxpi1, pypi1, pzpi1, pxspi, pyspi, pzspi)
           Else
              Call rotate(-pxpi1, -pypi1, -pzpi1, pxspi, pyspi, pzspi)
           End If
           !test off
           !               write(99,*) "2  ",pxspi,pyspi,pzspi
           patt(natt, 1) = pxspi
           patt(natt, 2) = pyspi
           patt(natt, 3) = pzspi
           patt(natt, 4) = sqrt(psembd**2+pimass**2)
           eatt = eatt + patt(natt, 4)
           gxar(natt) = xjet
           gyar(natt) = yjet
           gzar(natt) = 0.
           ftar(natt) = 0.
           itypar(natt) = katt(natt, 1)
           pxar(natt) = patt(natt, 1)
           pyar(natt) = patt(natt, 2)
           pzar(natt) = patt(natt, 3)
           pear(natt) = patt(natt, 4)
           xmar(natt) = pimass
        End Do
     End Do
  End If
  !lin-8/2009-end
  !
  Return
End Subroutine embedhighpt
