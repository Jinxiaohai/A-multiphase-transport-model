        subroutine oldcre(i) 
c       this subroutine is used to rearrange the old cell nic when a particle 
c       goes out of the cell
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
cc      SAVE /ilist1/
        SAVE   
        if (nic(i) .eq. 0) return
        j = nic(i)
        if (nic(j) .eq. i) then
           nic(j) = 0
           return
        end if
        do 10 while (nic(j) .ne. i)
           j = nic(j)
 10        continue
        nic(j) = nic(i)
        return
        end
