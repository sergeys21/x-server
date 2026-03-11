        Subroutine      Message (text,lines,iwait)
c--------------------------------------------------------
c   This subroutine outputs a text block to the screen
c
c text(lines) - text block consisting of lines lines
c iwait       - delay flag (0: message closes after 3 seconds
c                           1: program waits for "Hit any key...",
c                           2: beeps and then waits for "Hit any key..."
c
c                  Author: S.Stepanov
c--------------------------------------------------------
        Integer         lines, iwait, lwr(20),
     *                  lmax, l, i, iscan, iasci
        Character       text(lines)*(*)

c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------

        if (lines.lt.1 .or. lines.gt.20) goto 2  !--Error--+
                                                          !v
c This may cause overflows. Therefore, we siimply replace tabs
c by spaces:
c       Call    Tabexpan (text,lines)

        lmax = 0
        do      l=1,lines   !=================================+
          lwr(l) = Len_Trim(text(l))                         !|
          if (lwr(l).lt.1)      lwr(l)=1                     !|
          if (lwr(l).gt.lmax)   lmax = lwr(l)                !|
          if (lmax.gt.78) Then  !-+                           |
            lmax = 78            !|                           |
            lwr(l)=78            !|                           |
          endif  !----------------+                           |
          do i=1,lwr(l)  !==================================+ |
c           if (text(l)(i:i).lt.' ')      text(l)(i:i)=' ' !| |
c           if (text(l)(i:i).eq.char(11)) text(l)(i:i)=' ' !| |char(11)=vertical tabulation
            if (text(l)(i:i).eq.char(9))  text(l)(i:i)=' ' !| |char(9)=horizontal tabulation
          enddo  !==========================================+ |
        enddo  !==============================================+
c       if (lmax.lt.3 .or. lmax.gt.78)   goto 2  !--Error--+
        if (lmax.lt.3)                   goto 2  !--Error--+
                                                          !v

        if (modebat.eq.0)  Then !------------------------------+
c Write text:                                                  |
          do      l=1,lines  !=============+                   |
            Call  Text0 (lwr(l),text(l))  !|                   |
          enddo  !=========================+                   |
                                                              !|
          if (iwait.gt.0) Then !------------------+            |Hit key..
            Call Pause()                         !|            |
          else   !--------------------------------|            |
            do    i=1,15 !====================+   |            |Wait for 3s
              Call Sleep_tics (20,0)         !|   |            |15*0.20=3
              Call KeyIn  (0,iscan,iasci)    !|   |            |
              if (iscan.ne.0)     goto 7  !---+->-+--+         |
              if (iasci.ne.0)     goto 7  !---+->-+--|         |
            enddo  !==========================+   |  |         |
          endif    !------------------------------+  v         |
                                                    !|         |
  7       continue  !<-------------------------------+         |
          Return                                              !|
                                                              !|
        else  !------------------------------------------------|
                                                              !|
          Call  OpenErrorFile ()                              !|
          do  l=1,lines   !=============================+      |
            write (lunbat,10,err=12) text(l)(1:lwr(l)) !|      |
            cycle                                      !|      |
  12        continue                                   !|      |
            write (0,10) 'Message: error writing [',   !|      |
     *                   text(l)(1:lwr(l)),']'         !|      |
          enddo  !======================================+      |
c         Close (unit=lunbat)                                 !|
                                                              !|
          if (iwait.gt.0) Then !----------------------------+  |Hit key..
            Close (unit=lunbat,err=13)                     !|  |
  13        continue                                       !|  |
            Call exit_quiet()                              !|  |
          else   !------------------------------------------+  |
            Return                                         !|  |
          endif  !------------------------------------------+  |
                                                              !|
        endif  !-----------------------------------------------+
  10    format (10a)

  2     continue
        Stop ' MESSAGE: parameters not in range'
        End

c==============================================================

        Subroutine      OpenErrorFile ()
c -------------------------------------------------------
        Integer         i, j, io_status
        Integer         iarg

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------

        Call GetErrorFileName (j)
        if (j.lt.1) j=1
        if     (ErrorFile(1:j).eq.'stderr') Then !--+
           lunbat = 0                              !|
           goto 5                                  !|
        elseif (ErrorFile(1:j).eq.'stdout') Then !--+
           lunbat = 6                              !|
           goto 5                                  !|
        endif !-------------------------------------+

        if (lunbat.eq.6) Then !---------------------+
           ErrorFile = 'stdout'                    !|
           j = 6                                   !|
           goto 5                                  !|
        endif !-------------------------------------+

        if (lunbat.eq.0) lunbat=252
        Call OpenFile(ErrorFile(1:j),lunbat,'readwrite',
     *                           'unknown',io_status,*1)
c Append:
  3     continue      !<====================+
        read  (lunbat,'(a)',err=2,end=4)   !|
        goto 3  !==============>============+
  4     continue
        Backspace (lunbat,err=5)

  5     continue

c Get program name & path:
        iarg = 0
        Call get_arg_wrapper (iarg,txt(20))
        i = Max(Len_Trim(txt(20)),1)
        write (lunbat,'(2a)',err=6) '==> ', txt(20)(1:i)
  6     continue

        if (istackrezv.gt.0) Then  !-----------------------------+
          do i=1,istackrezv  !================================+  |
            write (lunbat,'(2a)',err=7) '==> ', progrezv(i)  !|  |
          enddo  !============================================+  |
        endif  !-------------------------------------------------+
  7     continue

        write (lunbat,'(2a)',err=8) '==> ', progname
  8     continue

        return
c--------------------------------------------------------------
c ERRORS:
  1     continue
        write (*,*) 'OpenErrorFile: cannot open the ERR file:'
        goto 10

  2     continue
        write (*,*) 'OpenErrorFile: cannot append the ERR file:'
        goto 10

  10    continue
        write (*,*) ErrorFile(1:j)
        stop '*** OpenErrorFile: exiting!'
        end

c==============================================================

        Subroutine      DeleteErrorFile ()
c -------------------------------------------------------
        Logical*4       idv
        Integer         j, lun

        Integer         FreeLun
        External        FreeLun

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------

        Call GetErrorFileName (j)
        if (j.lt.1) j=1

        if (ErrorFile(1:j).eq.'stderr' .OR.
     *      ErrorFile(1:j).eq.'stdout')   return

        Inquire (file=ErrorFile(1:j),exist=idv,err=1)
        if (.NOT.idv)   return

c       if (lunbat.eq.0) lunbat=252             !before introducing stderr and stdout
        lun = FreeLun()
        Open (unit=lun,
     *        file=ErrorFile(1:j),
     *        status='old',
c    *        share='denynone',                 !Compaq/MS Fortran only
c    *        shared,                           !Compaq/MS Fortran only
     *        err=2)
        Close(unit=lun,status='delete',err=3)
        lunbat = 0                              !reset lunbat
        return
c--------------------------------------------------------------
c ERRORS:
  1     continue
        write (*,*) 'DeleteErrorFile: cannot inquire the old ERR file:'
        goto 10

  2     continue
        write (*,*) 'DeleteErrorFile: cannot open the old ERR file:'
        goto 10

  3     continue
        write (*,*) 'DeleteErrorFile: cannot delete the old ERR file:'
        goto 10

  10    continue
        write (0,*) ErrorFile(1:j)
        Stop 'DeleteErrorFile: exiting!'
        end

c==============================================================

        Subroutine      GetErrorFileName (j)
c -------------------------------------------------------
        Integer         j, i, m
        Integer         iarg
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

        Character       getOS*8
        External        getOS
c -------------------------------------------------------
c Define ERR filename as EXE filename with ext=ERR
c (if not specified):
        if (Len_Trim(ErrorFile).eq.0  .Or.
     *      ErrorFile(1:1).eq.char(0)) Then !--------------------------+
c Get program name (incl. path):                                       |
          iarg = 0                                                    !|
          Call get_arg_wrapper (iarg,ErrorFile)                       !|
          j = Len_Trim(ErrorFile)                                     !|
c+-----------------------------------------------------------------+   |
c In principle, this is too much since EXE file                    |   |
c always has some name:                                            |   |
          if (getOS() .eq. 'windows' .OR.                         !|   |
     *        (j.ge.4 .AND. ErrorFile(j-3:j).eq.'.exe')) Then !-+  |   |
            if (j.lt.4) j=4                                    !|  |   |
            ErrorFile(j-3:j)='.err'                            !|  |   |instead of .exe
          else !------------------------------------------------+  |   |
            if (j.lt.2) j=2                                    !|  |   |
            if(j.gt.80-4) j=80-4                               !|  |   |
            ErrorFile(j+1:j+4)='.err'                          !|  |   |
          endif !-----------------------------------------------+  |   |
c+-----------------------------------------------------------------+   |
c         ErrorFile(j-3:j)='.err'                                     !|used before 02.2019
          do i=1,j  !=========================================+        |
            if (ErrorFile(i:i).eq.'\') ErrorFile(i:i) = '/'  !|        |
          enddo  !============================================+        |
c Remove path from the name:                                           |
          if (getOS() .eq. 'windows') Then !---------+                 |
            do i=j,1,-1  !========================+  |                 |
              if (ErrorFile(i:i).eq.'\' .or.     !|  |                 |
     *            ErrorFile(i:i).eq.'/' .or.     !|  |                 |
     *            ErrorFile(i:i).eq.':') goto 1 !-+--+--+              |
            enddo  !==============================+  |  |              |
          else !-------------------------------------+  |              |
            do i=j,1,-1  !========================+  |  v              |
              if (ErrorFile(i:i).eq.'/') goto 1 !-+--+--+              |
            enddo  !==============================+  |  |              |
          endif !------------------------------------+  |              |
c No path information found:                            |              |
          return                                       !|              |
c New length (without path):                            |              |
  1       continue !<-----------------------------------+              |
          j=j-i                                                       !|
c Shift the name to the left:                                          |
          do m=1,j  !===========================+                      |
            ErrorFile(m:m)=ErrorFile(m+i:m+i)  !|                      |
          enddo  !==============================+                      |
          m = Len(ErrorFile)                                          !|
          if (m.gt.j) ErrorFile(j+1:m)=' '                            !|
        endif !--------------------------------------------------------+
        j = Len_Trim (ErrorFile)
        return
        end
