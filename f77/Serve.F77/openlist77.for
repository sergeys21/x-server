        Subroutine      OpenList (name,lun,ifail)
c--------------------------------------------------------
c This subroutine attempts to open new file and if file exists,
c it offers 4 options: Backup, Overwrite, Append, Cancel.
c
c                  Author: S.Stepanov
c -------------------------------------------------------
c This is a version for GNU g77/gfortran
c -------------------------------------------------------
        Logical*4       idv
        Integer*4       ll
        Integer         lun, ifail, i, j, l, lf, lp, lx
        Character       name*(*), ext*4, buf*512

        Integer         Rename
c       External        Rename          !this must be an embedded function,
                                        !do not declare as external!
        Character       getOS*8
        External        getOS
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'OpenList'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        l     = Len_Trim(name)
        ifail = 0

        if (l.eq.0)     goto 11  !--------------------+
                                                     !V

        Inquire (file=name(1:l),exist=idv,err=15)

        if (.NOT.idv)   Then  !--------------------------+New file!
          Open (unit=lun,                               !|
     *          file=name(1:l),                         !|
     *          status='new',                           !|
c    *          share='denynone',                       !|Compaq/MS fortran only
c    *          shared,                                 !|Compaq/MS fortran only
     *          action='readwrite',                     !|
     *          err=12)                                 !|
          goto 100                          !goto return |
        endif  !-----------------------------------------+

c *** File already exists: **

        if (modebat.eq.0)  Then    !---------------------+NOT bat-mode
                                                        !|
          write (*,9) name(1:l)                         !|
  9       format(//' *** File [',a,
     *          '] already exists! ***'/)               !|
  8       continue  !<-------------------------+        !|
          write (*,7,advance='no')            !|        !|
  7       format(' Select action [1=Backup ',
     *    '2=Overwrite 3=Append 4=Cancel]_')  !|        !|
          read (*,'(i6)',err=8) i             !|        !|
          if (i.lt.1 .or. i.gt.4) goto 8 !-----+        !|
                                                        !|
          if (i.eq.4)      Then !---+                    |
            ifail = -1             !|ifail=-1 = Cancel   |
            goto 100               !|                    |
          endif  !------------------+                    |
                                                        !|
        else   !-----------------------------------------+bat-mode
                                                        !|
c BAT mode -- backup file!                              !|
          i = 1                                         !|
                                                        !|
        endif  !-----------------------------------------+

c i=1 Backup
c i=2 Overwrite
c i=3 Append

        if (i.ne.1)  Then  !--------------------------+
c Overwrite or Append:                                |
          Open  (unit=lun,                           !|
     *           file=name(1:l),                     !|
     *           status='unknown',                   !|
c    *           share='denynone',                   !|Compaq/MS fortran only
c    *           shared,                             !|Compaq/MS fortran only
     *           action='readwrite',                 !|
     *           err=12)                             !|
c Appending: position marker at the end of file.      |
          if (i.eq.3)     Then  !-----------------+   |
            ll = 1                               !|   |
  5         continue      !<=================+    |   |
            read  (lun,'(a)',err=13,end=6)  !|    |   |
            ll = ll + 1                     !|    |   |
            goto 5  !==============>=========+    |   |
  6         continue                             !|   |
            if (ll.gt.0)  Backspace (lun,err=14) !|   |
          endif  !--------------------------------+   |
          goto 100                                   !|
                                                     !|
        endif !---------------------------------------+

c Backup:
c ======
c       do j=1,l  !============================+        May not work with GNU Fortran Rename
c         if (name(j:j).eq.'\') name(j:j)='/' !|
c       enddo  !===============================+
c Find path:
        if (getOS() .eq. 'windows') Then !----+
          do j=l,1,-1  !====================+ |
            if (name(j:j).eq.'\' .OR.      !| |
     *          name(j:j).eq.'/' .OR.      !| |
     *          name(j:j).eq.':') goto 3 !--+-+--+
          enddo  !==========================+ |  |
        else !--------------------------------+  |
          do j=l,1,-1  !====================+ |  |
            if (name(j:j).eq.'/') goto 3 !--+-+--+
          enddo  !==========================+ |  |
        endif !-------------------------------+  |
        j = 0           !path not found          |
  3     continue  !<-----------------------------+pos. of 1st symbol after path
        lp = j+1
c Find available filename for backup:
        do j=l,lp,-1  !===================+
          if (name(j:j).eq.'.') goto 22 !-+--+
        enddo  !==========================+  |
        j  = l+1        !"." not found:      |
  22    continue  !<-------------------------+pos. of last symbol before "."
        lx = j-1
        ext = '.bak'
        buf = name(1:lx)//ext
        lf    = Len_Trim (buf)
        Inquire (file=buf(1:lf),exist=idv,err=15)
        if (idv)   Then  !---------------------------------------+
          ext = '.bkN'                                          !|
          do i=1,9  !========================================+   |
            write (ext(4:4),'(i1)')       i                 !|   |
            buf = name(1:lx)//ext                           !|   |
            Inquire (file=buf(1:lf),exist=idv,err=15)       !|   |
            if (.NOT.idv) goto 23   !------>----------->-----+---+---+
          enddo  !===========================================+   |   |
          ext = '.bNN'                                          !|   |
          do i=10,99  !======================================+   |   V
            write (ext(3:4),'(i2)')       i                 !|   |   |
            buf = name(1:lx)//ext                           !|   |   |
            Inquire (file=buf(1:lf),exist=idv,err=15)       !|   |   |
            if (.NOT.idv) goto 23   !------>----------->-----+---+---|
          enddo  !===========================================+   |   |
          goto 16         !all backup filenames are used         |   V
        endif  !-------------------------------------------------+   |
  23    continue  !<----------------<---------------<----------------+

        buf = name(lp:lx)//ext
        lf    = Len_Trim (buf)
        i     = Rename (name(1:l),buf(1:lf))
        if (i.ne.0)     goto 17

        Call Sleep_tics (1,1)   !sleep 1 second (DOS/serve.lib)

        Open (unit=lun,
     *        file=name(1:l),
     *        status='new',
c    *        share='denynone',                         !Compaq/MS fortran only
c    *        shared,                                   !Compaq/MS fortran only
     *        action='readwrite',
     *        err=12)
        goto 100

c ===========================================================
c                             ERRORS:
c ===========================================================
  11    continue
        ifail=1                 !ifail=1 -- no filename!
        buf = 'No filename'
        goto 99

  12    continue
        ifail=2                 !ifail=2 -- Open error
        buf = 'Open error'
        goto 99

  13    continue
        ifail=3                 !ifail=3 -- Read error
        buf = 'Read error'
        Close   (lun,err=99)
        goto 99

  14    continue
        ifail=4                 !ifail=4 -- Rewind error
        Close   (lun,err=99)
        buf = 'Rewind error'
        goto 99

  15    continue
        ifail=5                 !ifail=5 -- Inquire error
        buf = 'Inquire error'
        goto 99

  16    continue
        ifail=6                 !ifail=6 -- All backup names in use
        buf = 'All backup names (BAK, BK1...BK9, B11...B99) in use'
        goto 99

  17    continue
        ifail=7                 !ifail=7 -- Rename error
        buf='Rename error: '//name(1:l)//' to '//buf(1:lf)
        goto 99

  99    continue
        if (modebat.ne.0)       Then  !-----------+
          Call  OpenErrorFile ()                 !|
          write (lunbat,30,err=100) name(1:l)    !|
          j = Len_Trim(buf)                      !|
          write (lunbat,30,err=100) buf(1:j)     !|
  30      format (1x,a)                          !|
          Call exit_quiet()                      !|
        endif  !----------------------------------+

  100   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        Return
        End
