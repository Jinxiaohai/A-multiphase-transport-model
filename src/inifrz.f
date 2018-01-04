        subroutine inifrz
c
        implicit double precision  (a-h, o-z)
        PARAMETER (MAXPTN=400001)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
cc      SAVE /ilist5/
        common /frzprc/ 
     &       gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &       pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &       xmfrz(MAXPTN), 
     &       tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
cc      SAVE /frzprc/
        SAVE   
c
c     for freezeout time 0-10fm, use interval of 0.1fm; 
c     for 10-100fm, use interval of 1fm; 
c     for 100-1000fm, use interval of 10fm; 
c     for 1000-3000fm, use interval of 100fm: 
        step1=0.1d0
        step2=1d0
        step3=10d0
        step4=100d0
c     
        do 1001 it=1,101
           tfrz(it)=0d0+dble(it-1)*step1
 1001 continue
        do 1002 it=102,191
           tfrz(it)=10d0+dble(it-101)*step2
 1002   continue
        do 1003 it=192,281
           tfrz(it)=100d0+dble(it-191)*step3
 1003   continue
        do 1004 it=282,301
           tfrz(it)=1000d0+dble(it-281)*step4
 1004   continue
        tfrz(302)=tlarge
c
        return
        end
clin-5/2009 v2 analysis
c=======================================================================
c     idd=0,1,2,3 specifies different subroutines for partonic flow analysis.
