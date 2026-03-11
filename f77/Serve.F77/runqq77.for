        Integer*2 Function RunQQ(prg,args)
c----------------------------------------------------------
c This is an emulator of Compaq's RunQQ for GNU F77 fortran
c----------------------------------------------------------
        Character       prg*(*),args*(*), nude*24
        Integer         lprg, larg, istat, lnul

        Integer   System        !this function is Int*2 in DOS Fortran & Int*4 in Windows Fortran
c       External  System
c       Intrinsic System

        Character       getOS*8
        External        getOS

        lprg = Len_Trim(prg)
        if (lprg .lt. 1) Then !---+
           RunQQ = -1            !|
           return                !|
        endif !-------------------+

        if (getOS() .eq. 'windows') Then !---+
           nude = ' > NUL 2>&1'             !|
        else !-------------------------------+
           nude = ' > /dev/null 2>&1'       !|
        endif !------------------------------+
        lnul = Len_trim(nude)
        larg = Len_trim(args)
        if (larg .gt. 0) Then !------------------------------------------+
           istat = System(prg(1:lprg)//' '//args(1:larg)//nude(1:lnul)) !|
        else !-----------------------------------------------------------+
           istat = System(prg(1:lprg)//nude(1:lnul))                    !|
        endif  !---------------------------------------------------------+
        RunQQ = INT2(istat)
        return
        end
