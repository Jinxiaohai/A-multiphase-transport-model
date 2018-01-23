      SUBROUTINE coales
      PARAMETER (MAXSTR=150001)
      IMPLICIT DOUBLE PRECISION(D)
      DOUBLE PRECISION  gxp,gyp,gzp,ftp,pxp,pyp,pzp,pep,pmp
      DIMENSION IOVER(MAXSTR),dp1(2:3),dr1(2:3)
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      double precision  dpcoal,drcoal,ecritl
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
cc      SAVE /SOFT/
      common /coal/dpcoal,drcoal,ecritl
cc      SAVE /coal/
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
cc      SAVE /loclco/
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
cc      SAVE /HJJET2/
      SAVE
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  遍历所有的核子，初始化每个核子的IOVER.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 1001 ISG=1, NSG
         IOVER(ISG)=0
 1001 continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  再AMPT中首先是完成介子的强子化，然后才是重子的强子化。
cccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C1     meson q coalesce with all available qbar:
      do 150 ISG=1,NSG
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc NJSGS(I)(I=1, NSG) 第I个核子中部分子的个数， IOVER是个状态flag.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(NJSGS(ISG).ne.2.or.IOVER(ISG).eq.1) goto 150
C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc K2SGS(ISG, 1) 在第ISG个核子里面的第一个部分子的ID。(即正夸克q)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(K2SGS(ISG,1).lt.0) then
            write(6,*) 'Antiquark appears in quark loop; stop'
            stop
         endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc 得到介子的两个部分子信息。
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do 1002 j=1,2
            ftp(j)=ftsgs(isg,j)
            gxp(j)=gxsgs(isg,j)
            gyp(j)=gysgs(isg,j)
            gzp(j)=gzsgs(isg,j)
            pxp(j)=pxsgs(isg,j)
            pyp(j)=pysgs(isg,j)
            pzp(j)=pzsgs(isg,j)
            pmp(j)=pmsgs(isg,j)
            pep(j)=pesgs(isg,j)
 1002    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  下面的函数只是得到相同时刻两个价夸克的空间距离(二者的质心系)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         call locldr(2,drlocl)
         dr0=drlocl
c     dp0^2 defined as (p1+p2)^2-(m1+m2)^2:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  dp0为动量空间二者的动量距离(二者的实验室系)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         dp0=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         do 120 JSG=1,NSG
c     skip default or unavailable antiquarks:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  跳过自己本身和已经使用过的parton.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if(JSG.eq.ISG.or.IOVER(JSG).eq.1) goto 120
            if(NJSGS(JSG).eq.2) then
               ipmin=2
               ipmax=2
            elseif(NJSGS(JSG).eq.3.and.K2SGS(JSG,1).lt.0) then
               ipmin=1
               ipmax=3
            else
               goto 120
            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  计算qqbar之间的动量距离
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            do 100 ip=ipmin,ipmax
               dplocl=dsqrt(2*(pep(1)*pesgs(jsg,ip)
     1              -pxp(1)*pxsgs(jsg,ip)
     2              -pyp(1)*pysgs(jsg,ip)
     3              -pzp(1)*pzsgs(jsg,ip)
     4              -pmp(1)*pmsgs(jsg,ip)))
c     skip if outside of momentum radius:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  大于组合的动量半径, goto 120
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               if(dplocl.gt.dpcoal) goto 120
               ftp(2)=ftsgs(jsg,ip)
               gxp(2)=gxsgs(jsg,ip)
               gyp(2)=gysgs(jsg,ip)
               gzp(2)=gzsgs(jsg,ip)
               pxp(2)=pxsgs(jsg,ip)
               pyp(2)=pysgs(jsg,ip)
               pzp(2)=pzsgs(jsg,ip)
               pmp(2)=pmsgs(jsg,ip)
               pep(2)=pesgs(jsg,ip)
               call locldr(2,drlocl)
c     skip if outside of spatial radius:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  大于组合的空间半径, goto 120
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               if(drlocl.gt.drcoal) goto 120
c     q_isg coalesces with qbar_jsg:
               if((dp0.gt.dpcoal.or.dr0.gt.drcoal)
     1              .or.(drlocl.lt.dr0)) then
                  dp0=dplocl
                  dr0=drlocl
                  call exchge(isg,2,jsg,ip)
               endif
 100        continue
 120     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  满足条件，进行强子化的。IOVER(ISG) = 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(dp0.le.dpcoal.and.dr0.le.drcoal) IOVER(ISG)=1
 150  continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  下面的过程和上面的正好相反。
cccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C2     meson qbar coalesce with all available q:
      do 250 ISG=1,NSG
         if(NJSGS(ISG).ne.2.or.IOVER(ISG).eq.1) goto 250
C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
         do 1003 j=1,2
            ftp(j)=ftsgs(isg,j)
            gxp(j)=gxsgs(isg,j)
            gyp(j)=gysgs(isg,j)
            gzp(j)=gzsgs(isg,j)
            pxp(j)=pxsgs(isg,j)
            pyp(j)=pysgs(isg,j)
            pzp(j)=pzsgs(isg,j)
            pmp(j)=pmsgs(isg,j)
            pep(j)=pesgs(isg,j)
 1003    continue
         call locldr(2,drlocl)
         dr0=drlocl
         dp0=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         do 220 JSG=1,NSG
            if(JSG.eq.ISG.or.IOVER(JSG).eq.1) goto 220
            if(NJSGS(JSG).eq.2) then
               ipmin=1
               ipmax=1
            elseif(NJSGS(JSG).eq.3.and.K2SGS(JSG,1).gt.0) then
               ipmin=1
               ipmax=3
            else
               goto 220
            endif
            do 200 ip=ipmin,ipmax
               dplocl=dsqrt(2*(pep(2)*pesgs(jsg,ip)
     1              -pxp(2)*pxsgs(jsg,ip)
     2              -pyp(2)*pysgs(jsg,ip)
     3              -pzp(2)*pzsgs(jsg,ip)
     4              -pmp(2)*pmsgs(jsg,ip)))
c     skip if outside of momentum radius:
               if(dplocl.gt.dpcoal) goto 220
               ftp(1)=ftsgs(jsg,ip)
               gxp(1)=gxsgs(jsg,ip)
               gyp(1)=gysgs(jsg,ip)
               gzp(1)=gzsgs(jsg,ip)
               pxp(1)=pxsgs(jsg,ip)
               pyp(1)=pysgs(jsg,ip)
               pzp(1)=pzsgs(jsg,ip)
               pmp(1)=pmsgs(jsg,ip)
               pep(1)=pesgs(jsg,ip)
               call locldr(2,drlocl)
c     skip if outside of spatial radius:
               if(drlocl.gt.drcoal) goto 220
c     qbar_isg coalesces with q_jsg:
               if((dp0.gt.dpcoal.or.dr0.gt.drcoal)
     1              .or.(drlocl.lt.dr0)) then
                  dp0=dplocl
                  dr0=drlocl
                  call exchge(isg,1,jsg,ip)
               endif
 200        continue
 220     continue
         if(dp0.le.dpcoal.and.dr0.le.drcoal) IOVER(ISG)=1
 250  continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  重子的组合
cccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C3     baryon q (antibaryon qbar) coalesce with all available q (qbar):
      do 350 ISG=1,NSG
         if(NJSGS(ISG).ne.3.or.IOVER(ISG).eq.1) goto 350
         ibaryn=K2SGS(ISG,1)
C     DETERMINE CURRENT RELATIVE DISTANCE AND MOMENTUM:
         do 1004 j=1,2
            ftp(j)=ftsgs(isg,j)
            gxp(j)=gxsgs(isg,j)
            gyp(j)=gysgs(isg,j)
            gzp(j)=gzsgs(isg,j)
            pxp(j)=pxsgs(isg,j)
            pyp(j)=pysgs(isg,j)
            pzp(j)=pzsgs(isg,j)
            pmp(j)=pmsgs(isg,j)
            pep(j)=pesgs(isg,j)
 1004    continue
         call locldr(2,drlocl)
         dr1(2)=drlocl
         dp1(2)=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         ftp(2)=ftsgs(isg,3)
         gxp(2)=gxsgs(isg,3)
         gyp(2)=gysgs(isg,3)
         gzp(2)=gzsgs(isg,3)
         pxp(2)=pxsgs(isg,3)
         pyp(2)=pysgs(isg,3)
         pzp(2)=pzsgs(isg,3)
         pmp(2)=pmsgs(isg,3)
         pep(2)=pesgs(isg,3)
         call locldr(2,drlocl)
         dr1(3)=drlocl
         dp1(3)=dsqrt(2*(pep(1)*pep(2)-pxp(1)*pxp(2)
     &        -pyp(1)*pyp(2)-pzp(1)*pzp(2)-pmp(1)*pmp(2)))
c
         do 320 JSG=1,NSG
            if(JSG.eq.ISG.or.IOVER(JSG).eq.1) goto 320
            if(NJSGS(JSG).eq.2) then
               if(ibaryn.gt.0) then
                  ipmin=1
               else
                  ipmin=2
               endif
               ipmax=ipmin
            elseif(NJSGS(JSG).eq.3.and.
     1              (ibaryn*K2SGS(JSG,1)).gt.0) then
               ipmin=1
               ipmax=3
            else
               goto 320
            endif
            do 300 ip=ipmin,ipmax
               dplocl=dsqrt(2*(pep(1)*pesgs(jsg,ip)
     1              -pxp(1)*pxsgs(jsg,ip)
     2              -pyp(1)*pysgs(jsg,ip)
     3              -pzp(1)*pzsgs(jsg,ip)
     4              -pmp(1)*pmsgs(jsg,ip)))
c     skip if outside of momentum radius:
               if(dplocl.gt.dpcoal) goto 320
               ftp(2)=ftsgs(jsg,ip)
               gxp(2)=gxsgs(jsg,ip)
               gyp(2)=gysgs(jsg,ip)
               gzp(2)=gzsgs(jsg,ip)
               pxp(2)=pxsgs(jsg,ip)
               pyp(2)=pysgs(jsg,ip)
               pzp(2)=pzsgs(jsg,ip)
               pmp(2)=pmsgs(jsg,ip)
               pep(2)=pesgs(jsg,ip)
               call locldr(2,drlocl)
c     skip if outside of spatial radius:
               if(drlocl.gt.drcoal) goto 320
c     q_isg may coalesce with q_jsg for a baryon:
               ipi=0
               if(dp1(2).gt.dpcoal.or.dr1(2).gt.drcoal) then
                  ipi=2
                  if((dp1(3).gt.dpcoal.or.dr1(3).gt.drcoal)
     1                 .and.dr1(3).gt.dr1(2)) ipi=3
               elseif(dp1(3).gt.dpcoal.or.dr1(3).gt.drcoal) then
                  ipi=3
               elseif(dr1(2).lt.dr1(3)) then
                  if(drlocl.lt.dr1(3)) ipi=3
               elseif(dr1(3).le.dr1(2)) then
                  if(drlocl.lt.dr1(2)) ipi=2
               endif
               if(ipi.ne.0) then
                  dp1(ipi)=dplocl
                  dr1(ipi)=drlocl
                  call exchge(isg,ipi,jsg,ip)
               endif
 300        continue
 320     continue
         if(dp1(2).le.dpcoal.and.dr1(2).le.drcoal
     1        .and.dp1(3).le.dpcoal.and.dr1(3).le.drcoal)
     2        IOVER(ISG)=1
 350  continue
c
      RETURN
      END
c=======================================================================
