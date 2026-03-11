        Program grd2dat
c -------------------------------------------------------
c             ".grd"  --->  ".dat"
c
c  The .dat files are in the ASCII format as follows:
c                <x1>, <y1> <z1>
c                  :     :    :
c                <xN>, <yN> <zN>
c
c
c                Author: Sergey Stepanov
c -------------------------------------------------------
c Commmand line arguments:
c
c  grd2dat file.grd /B /L /R<number> /N /G
c where:
c       /B - batch (silent) opeation
c       /V - verbose opeation forced
c       /L - record GNU-file in logarithmic format (apply Log10)
c       /R - restrict data range to <n> orders of maximum (Log only)
c       /N - normalize logarithm to [0:10] range
c       /G - make speces between different columns for usage with GnuPlot
c -------------------------------------------------------
        Integer         ndim, nkeys
        Parameter       (ndim=2048)
        Parameter       (nkeys=6)

        Real*8          yr8(ndim), y8, Threshold,
     *                  ymin, ymip, ymax, ynorm

        Real            xn1, xk1,
     *                  xn2, xk2,
     *                  dx1, dx2,
     *                  x1,  x2

        Integer         ivalues(nkeys), io_status,
     *                  l, i, j, lest, lines,
     *                  ln, ls, nx1, nx2,
     *                  ifail, ncc,
     *                  lunin/1/, lunout/2/

        Character       inpname*80,outname*80,
     *                  textkeys(nkeys)*1,
     *                  mask*76 /'*.grd'/,
     *                  line(2)*80,
     *                  databuf*32766

        Logical         LogScale, Normalize, GnuSpaces

        Integer         n_Column
        External        n_Column

        Character       getOS*8
        External        getOS

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
        progname   = 'grd2dat'
        ErrorFile  = 'grd2dat.err'

        ln = 0                                  !keep GNU Fortran happy
        y8 = 0.                                 !keep GNU Fortran happy

        inpname     = ' '
        textkeys(1) = 'B'               !batch (silent) opeation
        textkeys(2) = 'V'               !verbose opeation forced
        textkeys(3) = 'L'               !apply logarithm
        textkeys(4) = 'R'               !restrict log range to <n>
        textkeys(5) = 'N'               !apply normalization to logarithm
        textkeys(6) = 'G'               !insert spaces into DAT for using with GnuPlot
c Get input parameters
c from the command line:
        Call GetInput (inpname,ivalues,textkeys,nkeys,ifail)
        if (ifail.ne.0) goto 504
c If /Key is specified, and <ivalue> is given,   ivalue=value
c If /Key is specified, but <ivalue> is missing, ivalue=-1
c If /Key is missing,                            ivalue=-2
        if (ivalues(1).ne.-2) then  !----+if key specified
          modebat = 1                   !|
        else  !--------------------------|
          modebat = 0                   !|
        endif  !-------------------------+
        if (ivalues(2).ne.-2) then  !----+if /V is specified,
          modebat = 0                   !|it overwrites /B
        endif  !-------------------------+
        if (ivalues(3).ne.-2) then  !----+if /L key specified
          LogScale = .true.             !|
        else  !--------------------------|
          LogScale = .false.            !|
        endif  !-------------------------+
        if (ivalues(4).gt.0) then  !-----+if /R key specified with pos.value
        if (ivalues(4).gt.32)           !|
     *                    ivalues(4)=32 !|
          Threshold = 10.0D0**ivalues(4)!|
        else  !--------------------------|
          Threshold = 0.0D0             !|
        endif  !-------------------------+
        if (ivalues(5).ne.-2) then  !----+if /N key specified
          Normalize = .true.            !|
        else  !--------------------------|
          Normalize = .false.           !|
        endif  !-------------------------+
        if (ivalues(6).ne.-2) then  !----+if /G key specified
          GnuSpaces = .true.            !|
        else  !--------------------------|
          GnuSpaces = .false.           !|
        endif  !-------------------------+
        if (.not. Logscale) Then !-------+
          Normalize = .false.           !|
          Threshold = 0.0D0             !|
        endif !--------------------------+

        if (modebat.eq.0) Write (*,123)
  123   Format  (1x,'This program converts GRDs to DAT files'/
     *           1x,39('-'))
c------------------------Read input file:
c Number of passed command line arguments:
        if (Len_Trim(inpname).eq.0)  Then !----------------------------+
          Call SelecFil (Mask,       !Dir_Mask                         |
     *                   inpname,    !Selected_File                    |
     *                   l,          !Selected_Indexes                 |
     *                   1,          !Max_Selected                     |
     *                   i,          !Num_Selected                     |
     *                   0,          !Multi_Flag                       |
     *                   5,3,        !Nx,Ny                            |
     *                   5,15)       !Ncolumns,Nlines                  |
          if (i.eq.0) Then !----------------------------------------+  |
            write (*,101)                                          !|  |
  101       format (
     *      'Syntax: grd2dat [/B /V /L /R<number> /N /G] file.grd'/
     *      '  where:'/
     *      '  /B - batch (silent) opeation'/
     *      '  /V - verbose opeation (overwtites batch)'/
     *      '  /L - record GNU-file in log format (apply Log10)'/
     *      '  /R - restrict data range to max/10^n (Log only)'/
     *      '  /N - normalize logarithm to [0:10] range'/
     *      '  /G - insert space lines for GnuPlot-3D'/)           !|  |
            Call exit_quiet()                                      !|  |
          endif !---------------------------------------------------+  |
        endif  !-------------------------------------------------------+
        ln  = Len_Trim(inpname)
        if (ln.eq.0) Stop 'grd2dat: Input file is not specified'
        if (ln.lt.4) Stop 'grd2dat: Input file extention is not .grd'
c Treat names with spaces presented in "...":
        if (inpname(1:1).eq.'"' .AND.
     *     inpname(ln:ln).eq.'"') Then  !------------------------+
          ln = ln-1                                             !|
          do i=1,ln  !========================+                  |
            inpname(i:i) = inpname(i+1:i+1)  !|                  |
          enddo  !============================+                  |
          inpname(ln:ln)=' '                                    !|
          ln = ln-1                                             !|
          if (ln.lt.4)                                          !|
     *      Stop 'grd2dat: Input filename extention is not .grd'!|
        endif  !-------------------------------------------------+
        do j=1,ln  !=================================+
          if (inpname(j:j).eq.'\') inpname(j:j)='/' !|
        enddo !======================================+
        databuf(1:ln) = inpname(1:ln)
        Call Case (inpname,1)                   !to lower case
        if (inpname(ln-3:ln).ne.'.grd')
     *              Stop 'grd2dat: Input filename extention is not .grd'
        inpname(1:ln) = databuf(1:ln)
c--------------------------
c Determine filename from full file specificatation
c that may also include path:
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
        ls = j+1
c--------------------------
        if (modebat.eq.0) write (*,21)  inpname(1:ln)
  21    format  (' Opening ',a)
        Call OpenFile(inpname(1:ln),lunin,'read','old',io_status,*501)
        read    (lunin,'(a)',err=503)   line(1)
        read    (lunin,  *,  err=503)   nx1, nx2
        read    (lunin,  *,  err=503)   xn1, xk1
        read    (lunin,  *,  err=503)   xn2, xk2
        read    (lunin,'(a)',err=503)   line(2)
        if (modebat.eq.0) write (*,27)  nx1,nx2
  27    format  (' Number of points:  nx1=',i4,'  nx2=',i4)
c Array overflow:
        if (nx1.gt.ndim)        goto 30

        if (Logscale) Then !-----------------------------------------+
          ymax =-1.E+37                                             !|
          ymin = 1.E+37                                             !|
          ymip = 1.E+37                                             !|
          lest = Len(databuf)                                       !|
          if (modebat.eq.0) write (*,*) 'Searching for [Min/Max]...'!|
          do j=1,nx2  !=========================================+    |
                                                               !|    |
            read (lunin,'(a)',err=503) databuf(1:lest)         !|    |
                                                               !|    |
            if (j.eq.1)   Then    !---------------------------+ |    |
c Determine the number of columns:                            | |    |
              lest = Min0(Len(databuf),Len_Trim(databuf)+50) !| |    |
c Column number mismatch:                                    !| |    |
              ncc  = n_Column (databuf(1:lest))              !| |    |
              if (ncc.ne.nx1)  goto 506  !--------------------+-+----+-+
            endif  !------------------------------------------+ |    | v
                                                               !|    |
            Call  RDreal8 (yr8,nx1,databuf(1:lest),ifail)      !|    |
            if (ifail.ne.0)       goto 503                     !|    |
            do    i=1,nx1  !======================+             |    |
              if (yr8(i).gt.ymax) ymax = yr8(i)  !|max number   |    |
              if (yr8(i).lt.ymin) ymin = yr8(i)  !|mim number   |    |
              if (yr8(i).lt.ymip .and.           !|             |    |
     *            yr8(i).gt.0.)   ymip = yr8(i)  !|min positive |    |
            enddo  !==============================+    number   |    |
          enddo  !==============================================+    |
          Rewind (unit=lunin)                                       !|
c Skip the header (i.e. position file cursor at the                  |
c place where it would be if no max/min pass needed:                 |
          do      j=1,5    !===========================+             |
            read  (lunin,'(a)',err=503) databuf(1:1)  !|             |
          enddo  !=====================================+             |
                                                                    !|
          if (modebat.eq.0) Then  !---------------------+            |
            write (*,*)                                !|            |
            write (*,*) 'Min = ',ymin,'  Max = ',ymax  !|            |
          endif  !--------------------------------------+            |
                                                                    !|
          if (abs(ymax-ymin).lt.1.E-32) Then !---------------+       |
            if (modebat.eq.0) Then  !---------------------+  |       |
              write (*,*)                                !|  |       |
              write (*,*) '!!! Plain GRD file: Min=Max!' !|  |       |
            endif  !--------------------------------------+  |       |
            if (Logscale) Then !------------------------+    |       |
              Logscale = .false.                       !|    |       |
              Normalize = .false.                      !|    |       |
              Threshold = 0.0D0                        !|    |       |
              if (modebat.eq.0)                        !|    |       |
     *          write (*,*) ' ...Disabling LOG!'       !|    |       |
            endif  !------------------------------------+    |       |
            if (Normalize) Then !-----------------------+    |       |
              Normalize = .false.                      !|    |       |
              if (modebat.eq.0)                        !|    |       |
     *          write (*,*) ' ...Disabling Normalize!' !|    |       |
            endif  !------------------------------------+    |       |
          endif  !-------------------------------------------+       |
                                                                    !|
          if (ymax .le. 0.) Then !----------------------------+      |
            if (Logscale) Then !----------------------------+ |      |
              Logscale = .false.                           !| |      |
              Normalize = .false.                          !| |      |
              Threshold = 0.0D0                            !| |      |
              ymip = ymin                                  !| |      |
              if (modebat.eq.0) Then  !-------------------+ | |      |
                write (*,*)                              !| | |      |
                write (*,*) '!!! No z>0: Setting NOLOG!' !| | |      |
              endif  !------------------------------------+ | |      |
            endif  !----------------------------------------+ |      |
          endif  !--------------------------------------------+      |
                                                                    !|
          if (LogScale) Then  !-------------------------------+      |
c Threshold correction                                        |      |
            if (ymip.gt.0. .and. Threshold.gt.0.) Then !--+   |      |
              if (modebat.eq.0) Then !----------------+   |   |      |
                write (*,*) 'Range = ',ymax/ymip,    !|   |   |      |
     *                    '  Threshold range = ',    !|   |   |      |
     *                                    Threshold  !|   |   |      |
              endif ! --------------------------------+   |   |      |
              if ((ymax/ymip) .gt. Threshold) Then !--+   |   |      |
                ymip = ymax / Threshold              !|   |   |      |
              endif  !--------------------------------+   |   |      |
              if (modebat.eq.0) Then !----------------+   |   |      |
                write (*,*) 'Min_cutoff = ',ymip     !|   |   |      |
              endif ! --------------------------------+   |   |      |
            endif  !--------------------------------------+   |      |
                                                             !|      |
            if (modebat.eq.0) Then !---+                      |      |
              write (*,*)             !|                      |      |
              write (*,*)             !|                      |      |
            endif  !-------------------+                      |      |
                                                             !|      |
c Cannot normalize plain files:                               |      |
            if (abs(ymax-ymip).lt.1.E-32) Normalize=.false.  !|      |
                                                             !|      |
c Normalization of LOG file:                                  |      |
            if (Normalize) Then !---------------------------+ |      |
c We map the whole LOG range into [0:10] interval:          | |      |
              ynorm = 10.                                  !| |      |
c This is 10^ynorm-1 :                                      | |      |
              y8 = dexp(dble(ynorm)*dlog(10.0D0)) - 1.0D0  !| |      |
              if (modebat.eq.0) Then !--------------------+ | |      |
                write (*,*)                              !| | |      |
                write (*,*) '[DAT] normalized to ',ynorm !| | |      |
                write (*,*)                              !| | |      |
              endif  !------------------------------------+ | |      |
            else  !-----------------------------------------| |      |
              ynorm = 0.                                   !| |      |
              y8 = 0.0D0                                   !| |      |
            endif !-----------------------------------------+ |      |
                                                             !|      |
          endif !---------------------------------------------+      |
                                                                    !|
          endif !----------------------------------------------------+

        if (modebat.eq.0)  write (*,*)

        outname = inpname(1:ln-4)//'.dat'

        if (modebat.eq.0) Then !---+
          write (*,*)             !|
          write (*,*)             !|
          write (*,21)  outname   !|
        endif  !-------------------+

c See serve/OpenList.FOR (open a new ASCII file
c with the test for "file exists"). Whether the
c program is interactive or not, depends on the
c modebat flag transferred via common /batmode/:
        Call    Openlist (outname,lunout,ifail)
        if (ifail.ne.0)     goto 502

        if (nx1.gt.1)   then  !-----------+
          dx1 = (xk1 - xn1)/(nx1-1)      !|
        else  !---------------------------|
          dx1 = 0.                       !|
        endif  !--------------------------+

        if (nx2.gt.1)   then  !-----------+
          dx2 = (xk2 - xn2)/(nx2-1)      !|
        else  !---------------------------|
          dx2 = 0.                       !|
        endif  !--------------------------+

        lest = Len(databuf)
        do j=1,nx2  !============================================+
          x2  =  xn2 + dx2*(j-1)                                !|
                                                                !|
          read  (lunin,'(a)',err=503) databuf(1:lest)           !|
                                                                !|
          if (j.eq.1)   Then    !----------------------------+   |
c Determine the number of columns:                           |   |
            lest = Min0 (Len(databuf),Len_Trim(databuf)+50) !|   |
            ncc  = n_Column (databuf(1:lest))               !|   |
            if (ncc.ne.nx1)  goto 506                       !|   |
          endif  !-------------------------------------------+   |
                                                                !|
          Call  RDreal8 (yr8,nx1,databuf(1:lest),ifail)         !|
          if (ifail.ne.0)       goto 503                        !|
                                                                !|
          if (LogScale) Then !--------------------------+        |
            do    i=1,nx1  !==========================+ |        |
              if (yr8(i).lt.ymip)  yr8(i)=ymip       !| |        |
              if (ynorm.gt.0.) then !---------------+ | |        |
                yr8(i) = y8*(yr8(i)-ymip)          !| | |        |
     /                /    dble(ymax-ymip) + 1.0D0 !| | |        |
              endif !-------------------------------+ | |        |
              yr8(i) = dlog10(yr8(i))                !| |        |
            enddo  !==================================+ |        |
          endif  !--------------------------------------+        |
                                                                !|
          do    i=1,nx1  !============================+          |
            x1 = xn1 + dx1*(i-1)                     !|          |
            write(lunout,444,err=505) x1, x2, yr8(i) !|          |
c 444       format (1x, 2g15.7, g16.7)               !|          |
  444       format (1x, 2g15.7, g16.7E3)             !|          |
          enddo  !====================================+          |
          if (GnuSpaces) Then !------------+                     |
            write (lunout,'(a)',err=505)  !|                     |
          endif !--------------------------+                     |
        enddo  !=================================================+
        close (lunin,err=99)
  99    continue
        close (lunout,err=100)
        if (modebat.eq.0) Then !-----------------------------+
          write (*,*)                                       !|
          write (*,*) 'File ',outname(1:ln),' is created.'  !|
        endif !----------------------------------------------+

  100   continue
        if (modebat.eq.0)  write (*,*)
        if (ifail .eq. 0) Then !--------+
          Call exit_quiet()            !|
        else !--------------------------+
          stop 1                       !|
        endif  !------------------------+
c =======================================================
c       ERRORS PROCESSING:

  30    continue
        txt(3) = 'File = '//inpname(1:ln)
        write   (txt(4),314) nx1,nx2
  314   format  ('Too many points in file:', 2i8)
        lines = 4
        goto    998

  501   continue
        txt(3) = 'File = '//inpname(1:ln)
        txt(4) = 'Input file open error'
        lines  = 4
        goto    998

  502   continue
        txt(3) = 'File = '//outname(1:ln)
        txt(4) = 'Output file open error (files already exist?)'
        lines  = 4
        goto    998

  503   continue
        txt(3) = 'File = '//inpname(1:ln)
        txt(4) = 'Input file read error'
        lines  = 4
        goto    998

  504   continue
        txt(3) = 'Error parsing command line arguments'
        lines  = 3
        goto    998

  505   continue
        txt(3) = 'File = '//outname(1:ln)
        txt(4) = 'Output file write error'
        lines  = 4
        goto    998

  506   continue
        write (txt,171) inpname(1:ln), nx1, ncc
  171   format('GRD2DAT  E R R O R:'//
     *         'File = ',a/
     *         'Columns mismatch in input file:'/
     *         'Claimed  = ',i7/
     *         'Found    = ',i7)
        lines  = 6
        goto    999

  998   continue
        txt(1) = 'GRD2DAT  E R R O R'
        txt(2) = ' '
  999   continue
        Call    Message (txt,lines,2)
        ifail = -1
        goto 100

        end
