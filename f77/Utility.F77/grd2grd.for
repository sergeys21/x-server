        Program grd2grd
c -------------------------------------------------------
c    Converting GRD files from Surf to WinSurf format
c            (adding or correcting z-header)
c
c                Author: Sergey Stepanov
c -------------------------------------------------------
c Command line arguments:
c
c  grd2gr2 file.grd /B
c where:
c       /B - batch (silent) operation
c -------------------------------------------------------
        Integer         ndim, nkeys
        Parameter       (ndim=2048)
        Parameter       (nkeys=1)

        Real*8          yr(ndim),
     *                  ymin, ymax

        Real            xn1, xk1,
     *                  xn2, xk2

        Integer         ivalues(nkeys), io_status,
     *                  l, i, j, lest, lines,
     *                  ln, ls, nx1, nx2,
     *                  ifail, ncc,
     *                  lunin/1/, lunout/2/

        Character       inpname*80,
     *                  TMPfile*80,
     *                  textkeys(nkeys)*1,
     *                  mask*76 /'*.grd'/,
     *                  line(2)*80,
     *                  bak*4 /'.gr_'/,
     *                  data*32766

        Integer         Rename
c       External        Rename  !this must be an embedded function, do not declare as external!

        Integer         n_Column
        External        n_Column

        Character       getOS*8
        External        getOS

        Logical*4       FileExist
        External        FileExist

        Character       txt(20)*80
        Common  /msg/   txt

c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------
        modebat    = 1         !assume batch-mode by default (for FOR32)
        istackrezv = 0
        progname   = 'grd2grd'
        ErrorFile  = 'grd2grd.err'

        ln = 0                                  !keep GNU Fortran happy

c Syntax: grd2grd grdfile [/B]
        inpname     = ' '
        textkeys(1) = "B"               !batch (silent) opeation
c Get input parameters
c from the command line:
        Call GetInput (inpname,ivalues,textkeys,nkeys,ifail)
        if (ifail.ne.0) goto 404
c If /Key is specified, and <ivalue> is given,   ivalue=value
c If /Key is specified, but <ivalue> is missing, ivalue=-1
c If /Key is missing,                            ivalue=-2
        if (ivalues(1).ne.-2) then  !----+if key specified
          modebat = 1                   !|
        else  !--------------------------|
          modebat = 0                   !|
        endif  !-------------------------+

        if (modebat.eq.0) Write (*,123)
  123   Format(1x,'This program adds Z-range',
     *           ' to GRD headers'/1x,40('-'))
c------------------------Read input file:
c Number of command line arguments:
        if (Len_Trim(inpname).eq.0) Then !--------------------------+
          if (modebat.eq.0)  Then !-----------------------------+   |
            Call SelecFil (Mask,       !Dir_Mask                |   |
     *                     inpname,    !Selected_File           |   |
     *                     l,          !Selected_Indexes        |   |
     *                     1,          !Max_Selected            |   |
     *                     i,          !Num_Selected            |   |
     *                     0,          !Multi_Flag              |   |
     *                     5,3,        !Nx,Ny                   |   |
     *                     5,15)       !Ncolumns,Nlines         |   |
            if (i.eq.0)  stop 'Syntax: grd2grd [/B] file.grd'  !|   |
            else   !--------------------------------------------+   |
              Call exit_quiet()                                !|   |
            endif  !--------------------------------------------+   |
          endif  !--------------------------------------------------+
        ln  = Len_Trim(inpname)
        if (ln.eq.0) Stop 'grd2grd: Input file is not specified'
        if (ln.lt.4) Stop 'grd2grd: Input file extention is not .grd'
        if (inpname(1:1).eq.'"' .AND.
     *     inpname(ln:ln).eq.'"') Then  !----------------+
          ln = ln-1                                     !|
          do i=1,ln  !========================+          |
            inpname(i:i) = inpname(i+1:i+1)  !|          |
          enddo  !============================+          |
          inpname(ln:ln)=' '                            !|
          ln = ln-1                                     !|
          if (ln.lt.4)                                  !|
     *      Stop 'Input filename extention is not .grd' !|
        endif  !-----------------------------------------+
        do j=1,ln  !=================================+
          if (inpname(j:j).eq.'\') inpname(j:j)='/' !|
        enddo !======================================+
        data(1:ln) = inpname(1:ln)
        Call Case (inpname,1)                   !to lower case
        if (inpname(ln-3:ln).ne.'.grd')
     *               Stop 'Input filename extention is not .grd'
        inpname(1:ln) = data(1:ln)
c--------------------------
c Distill filename from the full specification which also includes path:
        if (getOS() .eq. 'windows') Then !---------+
          do j=ln,1,-1 !========================+  |
            if (inpname(j:j).eq.'\' .OR.       !|  |
     *          inpname(j:j).eq.'/' .OR.       !|  |
     *          inpname(j:j).eq.':') goto 35 !--+--+--+
          enddo  !==============================+  |  |
        else !-------------------------------===---+  v
          do j=ln,1,-1 !========================+  |  |
            if (inpname(j:j).eq.'/') goto 35 !--+--+--+
          enddo  !==============================+  |  |
        endif !------------------------------------+  v
        j  = 0                                       !|
  35    continue  !<----------------------------------+
        ls = j+1                                !not used anymore (2019.01)
c--------------------------
        if (modebat.eq.0) Write (*,21) inpname(1:ln),'reading'
  21    format  (' Opening [',a,'] for ',a,' ...')
        Call OpenFile(inpname(1:ln),lunin,'read','old',io_status,*401)
        read    (lunin,'(a)',err=403)           line(1)
        read    (lunin,  *,  err=403)           nx1,nx2
        read    (lunin,  *,  err=403)           xn1,xk1
        read    (lunin,  *,  err=403)           xn2,xk2
        read    (lunin,'(a)',err=403)           line(2)

        if (modebat.eq.0) Write (*,22) nx1,nx2
  22    format  (' Number of rows nx1=',i4,',  number of lines nx2=',i4)
        if (nx1.gt.ndim) goto 30
        if (modebat.eq.0) Write (*,*)

        ymax =-1.E+32
        ymin =+1.E+32
        lest = Len(data)
        do      j=1,nx2  !====================================+
                                                             !|
          read  (lunin,'(a)',err=403) data(1:lest)           !|
                                                             !|
          if (j.eq.1)   Then    !---------------------+       |
c Determine the number of columns:                    |       |
            lest = Min0 (Len(data),Len_Trim(data)+50)!|       |
            ncc  = n_Column (data(1:lest))           !|       |
            if (ncc.ne.nx1)  goto 407                !|       |
          endif  !------------------------------------+       |
                                                             !|
          Call  RDreal8 (yr,nx1,data(1:lest),ifail)          !|
          if (ifail.ne.0)       goto 403                     !|
          do    i=1,nx1  !======================+             |
            if   (yr(i).gt.ymax) ymax = yr(i)  !|             |
            if   (yr(i).lt.ymin                !|             |
c For positive GRDs only:                      !|             |
c    *      .AND. yr(i).gt.0.                  !|             |
     *                         ) ymin = yr(i)  !|             |
          enddo  !==============================+             |
        enddo  !==============================================+
        Rewind (unit=lunin,err=403)
        do      j=1,5    !=======================+
          read  (lunin,'(a)',err=403) line(2)   !|
        enddo  !=================================+

        if (modebat.eq.0) Then !--------------------------+
c         Write (*,*)                                    !|
          Write (*,*)   'Min = ',ymin,'  Max = ',ymax    !|
          Write (*,*)                                    !|
c         Write (*,*)                                    !|
        endif !-------------------------------------------+

c 1) Produce .^^^ file with Z-scale determined
        TMPfile = inpname(1:ln-4)//'.^^^'
c See serve/OpenList.FOR (open a new ASCII file
c with the test for "file exists"). Whether the
c program is interactive or not, depends on the
c modebat flag transferred via common /batmode/:
c In the batchmode program backs up existing file as .bak/.bkN/.bNN:
        if (modebat.eq.0) Write (*,21) TMPfile(1:ln),'writing'
        Call    Openlist (TMPfile(1:ln),lunout,ifail)
        if (ifail.ne.0)     goto 402

        if (Len_Trim(line(1)) .eq. 0) line(1) = 'DSAA'
        Write   (lunout, '(a)', err=406)  line(1)(1:Len_Trim(line(1)))
        Write   (lunout,'(2i5)',err=406)  nx1,   nx2
        Write   (lunout,  44,   err=406)  xn1,   xk1
        Write   (lunout,  44,   err=406)  xn2,   xk2
        Write   (lunout,  45,   err=406)  ymin,  ymax
  44    format  (2g16.7)
  45    format  (2g16.7E3)

        do      j=1,nx2  !====================================+
                                                             !|
          read  (lunin,'(a)',err=403) data(1:lest)           !|
          l = Len_Trim (data(1:lest))                        !|
          if (l.lt.1) l = 1                                  !|
          Write (lunout,'(a)',err=406) data(1:l)             !|
        enddo  !==============================================+
        Close   (unit=lunout,err=406)
        Close   (unit=lunin, err=403)


c Now, we want to:
c  2) delete previous .gr_ if exists
c  3) then rename .grd -> .gr_,
c  4) then rename .^^^ -> .grd,
c  5) then delete .gr_

        if (FileExist(inpname(1:ln-4)//bak)) Then  !--------------+
          if (modebat.eq.0) Then  !---------------------+         |
            Write (*,*)                                !|         |
            Write (*,26) inpname(1:ln-4)//bak          !|         |
  26        Format(' Deleting if exists [',a,'] ...')  !|         |
          endif !---------------------------------------+         |
c 2a) Delete old BAK file (if present):                           |
          Call DelFile (inpname(1:ln-4)//bak,ifail)              !|
c 2b) Verify successful deleting of old BAK file:                 |
          if (FileExist(inpname(1:ln-4)//bak)) goto 551          !|
        endif !---------------------------------------------------+

c Rename function status codes:
c i=0          ! Success
c i=EACCES=13  ! Permission Denied
c i=ENOENT=2   ! File or path specified by from could not be found
c i=EXDEV=18   ! Attempt to move a file to a different device.
        j = 0
        i = 1
c 3a) Rename previous GRD file to BAK:
        do while (j .lt. 3 .AND. i .ne. 0)  !==========================+ up to 3 attempts
          if (modebat.eq.0) Then  !-----------------------------+      |
            Write (*,27) inpname(1:ln),inpname(1:ln-4)//bak    !|      |
  27        Format(' Renaming [',a,'] -> ['a,'] ...')          !|      |
          endif !-----------------------------------------------+      |
          i = Rename (inpname(1:ln), inpname(1:ln-4)//bak)            !|
          j = j+1                                                     !|
          Call Sleep_tics (1,1) !sleep 1 second (DOS/serve.lib)       !|
ccc       CAll Sleep (1)        !sleep 1 second (DOS/portlib.lib)     !|
        enddo  !=======================================================+
        if (i.ne.0)                                   goto 552

c 3b) Check renaming:
        if (modebat.eq.0) Write (*,28) inpname(1:ln-4)//bak,'is created'
  28    Format(' Checking if [',a,'] ',a,' ...')
        if (FileExist(inpname(1:ln-4)//bak)) Then !--------------------+
           if (modebat.eq.0) write (*,*) ' ...OK, file is created'    !|
        else !---------------------------------------------------------+
           if (modebat.eq.0) write (*,*) ' ...ERR, no such file'      !|
           goto 553                                                   !|
        endif !--------------------------------------------------------+
        if (modebat.eq.0) Write (*,28) inpname(1:ln),'is renamed'
        if (.NOT. FileExist(inpname(1:ln))) Then !---------------------+
           if (modebat.eq.0) write (*,*) ' ...OK, file is removed'    !|
        else !---------------------------------------------------------+
           if (modebat.eq.0) write (*,*) ' ...ERR, file still exists' !|
           Call DelFile (inpname(1:ln),ifail)                         !|try to delete
           if (FileExist(inpname(1:ln))) goto 554                     !|deleting failed
        endif !--------------------------------------------------------+

c Rename function status codes:
c i=0          ! Success
c i=EACCES=13  ! Permission Denied
c i=ENOENT=2   ! File or path specified by from could not be found
c i=EXDEV=18   ! Attempt to move a file to a different device.

c 4a) Rename temp file .^^^ as new GRD:
        if (modebat.eq.0) Write (*,27) TMPfile(1:ln), inpname(1:ln)
        i = Rename (TMPfile(1:ln), inpname(1:ln))
        if (i.ne.0)                                   goto 555
        Call Sleep_tics (1,1)           !sleep 1 second (DOS/serve.lib)
c ccc   Call Sleep (1)                  !sleep 1 second (Win/portlib.lib)

c 4b) Check renaming:
        if (modebat.eq.0) Write (*,28) inpname(1:ln),'is created'
        if (FileExist(inpname(1:ln)) ) Then !--------------------------+
           if (modebat.eq.0) write (*,*) ' ...OK, file is created'    !|
        else !---------------------------------------------------------+
           if (modebat.eq.0) write (*,*) ' ...ERR, no such file'      !|
           goto 556                                                   !|
        endif !--------------------------------------------------------+
        if (modebat.eq.0) Write (*,28) TMPfile(1:ln),'is renamed'
        if (.NOT. FileExist(TMPfile(1:ln))) Then !---------------------+
           if (modebat.eq.0) write (*,*) ' ...OK, file is removed'    !|
        else !---------------------------------------------------------+
           if (modebat.eq.0) write (*,*) ' ...ERR, file still exists' !|
           Call DelFile (TMPfile(1:ln),ifail)                         !|try to delete
           if (FileExist(TMPfile(1:ln))) goto 557                     !|deleting failed
        endif !--------------------------------------------------------+

c 5) Delete BAK file (.gr_):
        if (modebat.eq.0) Write (*,26) inpname(1:ln-4)//bak
        Call DelFile (inpname(1:ln-4)//bak,ifail)

c 6) Delete temp .^^^ file:
c ccccc Call DelFile (TMPfile(1:ln),ifail)

  100   continue
        if (modebat.eq.0)  Write (*,*)
        Call exit_quiet()
c =======================================================
c       ERROR PROCESSING:

  30    continue
        txt(3) = 'File = '//inpname(1:ln)
        Write   (txt(4),314) nx1,nx2
  314   format  ('Too many points in input file:', 2i8)
        lines  = 4
        goto    998

  401   continue
        txt(3) = 'File = '//inpname(1:ln)
        txt(4) = 'Input file opening error!'
        lines  = 4
        goto    998

  402   continue
        txt(3) = 'File = '//TMPfile(1:ln)
        txt(4) = 'Output file opening error (file already exists?)'
        lines  = 4
        goto    998

  403   continue
        txt(3) = 'File = '//inpname(1:ln)
        txt(4) = 'Input file read error'
        lines  = 4
        goto    998

  404   continue
        txt(3) = 'Error parsing command line arguments'
        lines  = 3
        goto    998

  406   continue
        txt(3) = 'File = '//TMPfile(1:ln)
        txt(4) = 'Output file write error'
        lines  = 4
        goto    998

  407   continue
        write (txt,171) inpname(1:ln), nx1, ncc
  171   format('GRD2GRD  E R R O R:'//
     *         'File = ',a/
     *         'Columns mismatch in input file:'/
     *         'Claimed  = ',i7/
     *         'Found    = ',i7)
        lines  = 6
        goto    999

  551   continue
        txt(3) = 'Deleting previous BAK file failed!'
        txt(4) = 'File = '//inpname(1:ln-4)//bak
        txt(5) = 'See result in '//TMPfile(1:ln)
        lines  = 5
        goto    998

  552   continue
        write (txt(3),662) i
  662   format('Renaming GRD to BAK failed with status = ',i3,'!')
        txt(4) = 'GRD='//inpname(1:ln)
        txt(5) = 'BAK='//inpname(1:ln-4)//bak
        txt(6) = 'See result in '//TMPfile(1:ln)
        lines  = 6
        goto    998

  553   continue
        txt(3) = 'Check renaming GRD to BAK failed (new file not found)'
        txt(4) = 'GRD='//inpname(1:ln)
        txt(5) = 'BAK='//inpname(1:ln-4)//bak
        txt(6) = 'See result in '//TMPfile(1:ln)
        lines  = 6
        goto    998

  554   continue
        txt(3) = 'Check renaming GRD to BAK failed (old file exists)'
        txt(4) = 'GRD='//inpname(1:ln)
        txt(5) = 'BAK='//inpname(1:ln-4)//bak
        txt(6) = 'See result in '//TMPfile(1:ln)
        lines  = 6
        goto    998

  555   continue
        write (txt(3),665) i
  665   format('Renaming TMP to GRD failed with status = ',i3,'!')
        txt(4) = 'TMP='//TMPfile(1:ln)
        txt(5) = 'GRD='//inpname(1:ln)
        txt(6) = 'See result in '//TMPfile(1:ln)
        lines  = 6
        goto    998

  556   continue
        txt(3) = 'Check renaming TMP to GRD failed (new file not found)'
        txt(4) = 'TMP='//TMPfile(1:ln)
        txt(5) = 'GRD='//inpname(1:ln)
        txt(6) = 'See result in '//TMPfile(1:ln)
        lines  = 6
        goto    998

  557   continue
        txt(3) = 'Check renaming TMP to GRD failed (old file exists)'
        txt(4) = 'TMP='//TMPfile(1:ln)
        txt(5) = 'GRD='//inpname(1:ln)
        txt(6) = 'See result in '//TMPfile(1:ln)
        lines  = 6
        goto    998

  998   continue
        txt(1) = 'GRD2GRD  E R R O R:'
        txt(2) = ' '
  999   continue
        Call Message (txt,lines,2)
        goto 100
        end
