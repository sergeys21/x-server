        Subroutine OpenBinary (filename, lun, filestate,
     *                         file_rw, io_status, ifail)
c--------------------------------------------------------
c Opens binary file for reading or writing depending on the input flags.
c Returns io_status and status
c--------------------------------------------------------
        Integer   lun, io_status, ifail,
     *            lname, lstate, lrw
        Character filename*(*),                 !requsted file name
     *            filestate*(*),                !OLD/NEW/UNKNOWN
     *            file_rw*(*)                   !READ/WRITE/READWRITE (same as "mode")

        ifail = 0
        io_status = 0
        lname  = Len_Trim(filename)
        lstate = Len_Trim(filestate)
        lrw    = Len_Trim(file_rw)
        if (lname .eq. 0) Then !-----+
          ifail = -1                !|          !no name specified
          return                    !|
        endif  !---------------------+
        if (lstate .lt. 3) Then !----+          !OLD/NEW/UNKNOWN
          ifail = -2                !|          !incorrect filestate
          return                    !|
        endif  !---------------------+
        if (lrw .lt. 4) Then !-------+          !READ/WRITE/READWRITE
          ifail = -3                !|          !incorrect filestate
          return                    !|
        endif  !---------------------+

        Open (unit=lun,
     *        file=filename(1:lname),
     *        status=filestate(1:lstate),
c    *        form='binary',                    !Compaq/MS fortran only
     *        form='unformatted',               !'unformatted+stream' in GnuFortran
     *        access='stream',                  !replace 'binary' in Compaq fortran
     *        action=file_rw(1:lrw),            !Compaq/GNU fortran only
c    *        mode=file_rw(1:lrw),              !Compaq/MS fortran only
c    *        share='denynone',                 !Compaq/MS fortran only
     *        iostat=io_status,
     *        err=1)
        return
  1     continue
        ifail = 1
        return
        end

c =====================================================================
        Subroutine      OpenBin (filename,lun,ifail)
c--------------------------------------------------------
c This subroutine attempts to open new BINARY file and if the file
c exists, it offers 4 options: Backup, Overwrite, Cancel. [NO APPEND!]
c
c                  Author: S.Stepanov
c -------------------------------------------------------
c This is a version for GNU g77/gfortran
c -------------------------------------------------------
        Logical*4       idv
        Integer         lun, ifail, i, j, l, lf, lp, lx
        Character       filename*(*), ext*4, buf*512

        Integer         Rename
c       External        Rename          !this must be an embedded function, do not declare as external!

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
          progname = 'OpenBin'            !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        l     = Len_Trim(filename)
        ifail = 0

        if (l.eq.0)     goto 11  !--------------------+
                                                     !V

        Inquire (file=filename(1:l),exist=idv,err=15)

        if (.NOT.idv)   Then  !-----+New file!
          Open (unit=lun,           !|
     *          file=filename(1:l), !|
c    *          form='binary',      !|Compaq/MS fortran only
     *          form='unformatted', !|'unformatted+stream' in gfortran
     *          access='stream',    !|replace 'binary' of Compaq fortran
     *          status='new',       !|
c    *          share='denynone',   !|Compaq/MS fortran only
c    *          shared,             !|Compaq/MS fortran only
c    *          mode='write',       !|Compaq/MS fortran only
     *          action='write',     !|Compaq/Gnu fortran only
     *          err=12)             !|
          goto 100                  !|goto return
        endif  !---------------------+

c *** File already exists: **

        if (modebat.eq.0)  Then    !------------------+NOT bat-mode
            write (*,9) filename(1:l)                !|
  9       format(//' *** File [',a,
     *          '] already exists! ***'/)            !|
  8       continue  !--------------------------+      |
          write (*,7,advance='no')            !|      |no carriage return at the end
c         write (*,7)                         !|      |need to suppress CR in format
  7       format(' Select action [1=Backup ',
     *    '2=Overwrite 3=Cancel]_')           !|      |
c    *    '2=Overwrite 3=Cancel]_',\)         !|      |suppress CR: MS/Compaq
c    *    '2=Overwrite 3=Cancel]_',$)         !|      |suppress CR: GNU/Compaq
c 8       Call Keyin (1,iscan,iasci)          !|      |
c         if (iasci.eq.27) iasci=52           !|      |
c         i = iasci-48                        !|      |
c                                             !|      |
          read (*,'(i6)',err=8) i             !|      |
          if (i.lt.1 .or. i.gt.3) goto 8 !-----+      |
                                                     !|
          if (i.eq.3)  Then !---+                     |
            ifail = -1         !|ifail=-1 = Cancel    |
            goto 100           !|                     |
          endif  !--------------+                     |
                                                     !|
        else   !--------------------------------------+
                                                     !|
c BAT mode -- backup file!                           !|
          i = 1                                      !|
                                                     !|
        endif  !--------------------------------------+

c i=1 Backup
c i=2 Overwrite

        if (i.eq.2)  Then  !----------+
c Overwrite:                          |
          Open (unit=lun,            !|
     *          file=filename(1:l),  !|
c    *          form='binary',       !|Compaq/MS fortran only
     *          form='unformatted',  !|'unformatted+stream' in gfortran
     *          access='stream',     !|replace 'binary' of Compaq fortran
     *          status='unknown',    !|
c    *          share='denynone',    !|Compaq/MS fortran only
c    *          shared,              !|Compaq/MS fortran only
c    *          mode='readwrite',    !|Compaq/MS fortran only
c    *          action='readwrite',  !|Compaq/GNU fortran only
     *          err=12)              !|
          goto 100                   !|
        endif !-----------------------+

c Backup:
c ======
c Find path:
        if (getOS() .eq. 'windows') Then !--------+
          do j=l,1,-1  !========================+ |
            if (filename(j:j).eq.'\' .OR.      !| |
     *          filename(j:j).eq.'/' .OR.      !| |
     *          filename(j:j).eq.':') goto 3 !--+-+--+
          enddo  !==============================+ |  |
        else !------------------------------------+  |
          do j=l,1,-1  !========================+ |  |
            if (filename(j:j).eq.'/') goto 3 !--+-+--+
          enddo  !==============================+ |  |
        endif !-----------------------------------+  |
        j = 0           !path not found              |
  3     continue  !<---------------------------------+pos. of 1st symbol after path
        lp = j+1
c Find available filename for backup:
        do j=l,lp,-1  !=======================+
          if (filename(j:j).eq.'.') goto 22 !-+--+
        enddo  !==============================+  |
        j  = l+1        !"." not found:          |
  22    continue  !<-----------------------------+pos. of last symbol before "."
        lx = j-1
        ext = '.bak'
        buf = filename(1:lx)//ext
        lf    = Len_Trim (buf)
        Inquire (file=buf(1:lf),exist=idv,err=15)
        if (idv)   Then  !---------------------------------------+
          ext = '.bkN'                                          !|
          do i=1,9  !========================================+   |
            write (ext(4:4),'(i1)')       i                 !|   |
            buf = filename(1:lx)//ext                       !|   |
            Inquire (file=buf(1:lf),exist=idv,err=15)       !|   |
            if (.NOT.idv) goto 23   !------>----------->-----+---+---+
          enddo  !===========================================+   |   |
          ext = '.bNN'                                          !|   |
          do i=10,99  !======================================+   |   V
            write (ext(3:4),'(i2)')       i                 !|   |   |
            buf = filename(1:lx)//ext                       !|   |   |
            Inquire (file=buf(1:lf),exist=idv,err=15)       !|   |   |
            if (.NOT.idv) goto 23   !------>----------->-----+---+---|
          enddo  !===========================================+   |   |
          goto 16         !all backup filenames are used         |   V
        endif  !-------------------------------------------------+   |
  23    continue  !<----------------<---------------<----------------+

        buf = filename(lp:lx)//ext
        lf    = Len_Trim (buf)
        i     = Rename (filename(1:l),buf(1:lf))
        if (i.ne.0)     goto 17
        Open (unit=lun,
     *        file=filename(1:l),
c    *        form='binary',        !Compaq/MS fortran only
     *        form='unformatted',   !'unformatted+stream' in gfortran
     *        access='stream',      !replace 'binary' of Compaq fortran
     *        status='new',
c    *        share='denynone',     !Compaq/MS fortran only
c    *        shared,               !Compaq/MS fortran only
c    *        mode='readwrite',     !Compaq/MS fortran only
c    *        action='readwrite',   !Compaq/GNU fortran only
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
c       buf='Rename error: no disk space or cannot start command.com'
        buf='Rename error: '//filename(1:l)//' to '//buf(1:lf)
        goto 99

  99    continue
        if (modebat.ne.0)       Then  !------------+
          Call  OpenErrorFile ()                  !|
          write (lunbat,30,err=100) filename(1:l) !|
          j = Len_Trim(buf)                       !|
          write (lunbat,30,err=100) buf(1:j)      !|
  30      format (1x,a)                           !|
          Call exit_quiet()                       !|
        endif  !-----------------------------------+

  100   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        Return
        End
