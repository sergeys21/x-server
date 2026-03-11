        Integer Function N_prg_running (prg)
c-----------------------------------------------------------------------
c Find the number of running instances of the program
c using the Cygwin "ps":
c c:/cygwin/bin/ps -aW -s | c:/cygwin/bin/grep -i prgname
c
c                 Author: S.Stepanov
c-----------------------------------------------------------------------
        Character prg*(*), str*80
        Integer   lprg, ifail, lun, io_status

        Character       getOS*8
        External        getOS

        Logical*4 FileExist
        External  FileExist

        Integer   FreeLun
        External  FreeLun

        Integer   System
c       External  System
c       Intrinsic System

        N_prg_running = 0
        lprg = Max(Len_Trim(prg),1)

        if (FileExist('Nrunning.tmp')) then !--+
          Call DelFile('Nrunning.tmp',ifail)  !|
          if (ifail.ne.0) goto 1              !|
        endif  !-------------------------------+

        str = prg(1:lprg)
        if(lprg .gt. Len(str)) lprg=Len(str)
        if (getOS() .eq. 'windows') Then !-----------------------------+
           ifail = System('c:\cygwin\bin\ps.exe -aWs | '//            !|
     *                    'c:\cygwin\bin\grep.exe -i '//str(1:lprg)// !|
     *                    ' > Nrunning.tmp 2>NUL')                    !|
        else !---------------------------------------------------------+
           ifail = System('/bin/ps -ef | '//                          !|
     *                    '/bin/grep -i '//str(1:lprg)//              !|
     *                    ' > Nrunning.tmp 2>/dev/null')              !|
        endif !--------------------------------------------------------+
c       if (ifail.ne.0) goto 1

        lun = FreeLun()
        Call OpenFile('Nrunning.tmp',lun,'readwrite','old',io_status,*1)
  10    continue  !<-----------------------------+
        read    (lun,'(a)',err=2,end=2) str     !|
        if (Len_Trim(str).eq.0) goto 10         !|
        N_prg_running = N_prg_running + 1       !|
        goto 10  !-------------------------------+

  2     continue
        close(unit=lun, status='delete',err=1)
c       close(unit=lun)
  1     continue
        return
        end
