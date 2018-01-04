        subroutine inipar
        implicit double precision (a-h,o-z)
        common /para4/ iftflg, ireflg, igeflg, ibstfg
cc      SAVE /para4/
        common /para6/ centy
cc      SAVE /para6/
        SAVE   
        if (ibstfg .ne. 0) then
           centy = -6d0
        end if
        return
        end
