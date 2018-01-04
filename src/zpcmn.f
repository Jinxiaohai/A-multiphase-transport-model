      SUBROUTINE ZPCMN
c       Version: 1.0.1
c       Author: Bin Zhang 
c       (suggestions, problems -> bzhang@nt1.phys.columbia.edu)
        implicit double precision (a-h, o-z)
clin-4/20/01        PARAMETER (NMAXGL = 16000)
        parameter (MAXPTN=400001)
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
cc      SAVE /para3/
        SAVE   
c
c       loop over events
        do 1000 i = 1, nevnt
           ievt = i
c       generation of the initial condition for one event
           call inievt
c      loop over many runs of the same event
           do 2000 j = 1, nsbrun
              isbrun = j
c       initialization for one run of an event
              call inirun
clin-4/2008 not used:
c             CALL HJAN1A
 3000         continue
c       do one collision
              call zpcrun(*4000)
              call zpca1
              goto 3000
 4000         continue
              call zpca2
 2000      continue
 1000   continue
        call zpcou
clin-5/2009 ctest off
c     5/17/01 calculate v2 for parton already frozen out:
c        call flowp(3)
c.....to get average values for different strings
        CALL zpstrg
        RETURN
        end
******************************************************************************
******************************************************************************
