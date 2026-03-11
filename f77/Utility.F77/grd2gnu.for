        Program grd2gnu
c -------------------------------------------------------
c             ".grd"  --->  ".GNU" (GnuPlot binary file)
c
c In the GnuPlot binary file single precision floats
c are stored as follows:
c
c     <N+1>  <y0>   <y1>   <y2>  ...  <yN>
c      <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
c      <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
c       :      :      :      :   ...    :
c
c                Author: Sergey Stepanov
c -------------------------------------------------------
c Command line arguments:
c
c  grd2gnu file.grd /B /L /R<number> /N
c where:
c       /B - batch (silent) operation
c       /V - verbose opeation forced
c       /L - record GNU-file in logarithmic format (apply Log10)
c       /R - restrict data range to <n> orders of maximum (Log only)
c       /N - normalize logarithm to [0:10] range
c -------------------------------------------------------
        Integer         ndim, nkeys
        Parameter       (ndim=2048)
        Parameter       (nkeys=5)

        Real*8          yr8(ndim), y8, Threshold,
     *                  ymin, ymip, ymax, ynorm

        Real            yr4(ndim),
     *                  xn1, xk1, xn2, xk2,
     *                  dx1, dx2, x1,  x2,
     *                  nx1_r

        Integer         ivalues(nkeys), io_status,
     *                  l, i, j, lest, ifail,
     *                  ln, ls, nx1, nx2,
     *                  ncc, lines,
     *                  lunin/1/, lunout/2/

        Character       inpname*80,outname*80,
     *                  textkeys(nkeys)*1,
     *                  mask*76/'*.grd'/,
     *                  databuf*32766

        Logical         LogScale, Normalize

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
        progname   = 'grd2gnu'
        ErrorFile  = 'grd2gnu.err'

        ln = 0                                  !keep GNU Fortran happy
        y8 = 0.                                 !keep GNU Fortran happy

        inpname     = ' '
        textkeys(1) = 'B'               !batch (silent) opeation
        textkeys(2) = 'V'               !verbose opeation forced
        textkeys(3) = 'L'               !apply logarithm
        textkeys(4) = 'R'               !restrict log range to <n>
        textkeys(5) = 'N'               !apply normalization to logarithm
c Get input parameters
c from the command line:
        Call GetInput (inpname,ivalues,textkeys,nkeys,ifail)
        if (ifail.ne.0) goto 504
c If /Key is specified, and <ivalue> is given,   ivalue=value
c If /Key is specified, but <ivalue> is missing, ivalue=-1
c If /Key is missing,                            ivalue=-2
        if (ivalues(1).ne.-2) then !-----+if key specified
          modebat = 1                   !|
        else  !--------------------------|
          modebat = 0                   !|
        endif  !-------------------------+
        if (ivalues(2).ne.-2) then !-----+if key specified, it
          modebat = 0                   !|overwrites /B
        endif  !-------------------------+
        if (ivalues(3).ne.-2) then !-----+if key specified
          LogScale = .true.             !|
        else  !--------------------------|
          LogScale = .false.            !|
        endif  !-------------------------+
        if (ivalues(4).gt.0) then !------+if key specified with pos.value
        if (ivalues(4).gt.32)           !|
     *                    ivalues(4)=32 !|
          Threshold = 10.0D0**ivalues(4)!|
        else  !--------------------------|
          Threshold = 0.0D0             !|
        endif  !-------------------------+
        if (ivalues(5).ne.-2) then !-----+if key specified
          Normalize = .true.            !|
        else  !--------------------------|
          Normalize = .false.           !|
        endif  !-------------------------+

        if (modebat.eq.0)  Write (*,123)
  123   Format  (1x,'This program converts GRD to GNUplot binary file'/
     *           1x,48('-'))
c------------------------Read input file:
c Number of passed command line arguments:
        if (Len_Trim(inpname).eq.0)  Then  !---------------------------+
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
     *      'Syntax: grd2gnu [/B /V /L /S /R<number>] file.grd'/
     *      '  where:'/
     *      '  /B - batch (silent) opeation'/
     *      '  /V - verbose opeation (overwtites batch)'/
     *      '  /L - record GNU-file in log format (apply Log10)'/
     *      '  /R - restrict data range to max/10^n (Log only)'/
     *      '  /N - normalize logarithm to [0:10] range'/)         !|  |
            Call exit_quiet()                                      !|  |
          endif !---------------------------------------------------+  |
        endif  !-------------------------------------------------------+
        ln  = Len_Trim(inpname)
        if (ln.eq.0) Stop 'grd2gnu: Input file is not specified'
        if (ln.lt.4) Stop 'grd2gnu: Input file extention is not .grd'
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
     *      Stop 'grd2gnu: Input filename extention is not .grd'!|
        endif  !-------------------------------------------------+
        do j=1,ln  !=================================+
          if (inpname(j:j).eq.'\') inpname(j:j)='/' !|
        enddo !======================================+
        databuf(1:ln) = inpname(1:ln)
        Call Case (inpname,1)                   !to lower case
        if (inpname(ln-3:ln).ne.'.grd')
     *                  Stop 'grd2gnu: Input file extention is not .grd'
        inpname(1:ln) = databuf(1:ln)
c--------------------------
c Determine filename from full file specifications
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
        if (modebat.eq.0) write (*,21) inpname(1:ln)
  21    format  (' Opening ',a)
        Call OpenFile(inpname(1:ln),lunin,'read','old',io_status,*501)
        read    (lunin,'(a)',err=503)           databuf(1:1)
        read    (lunin,  *,  err=503)           nx1, nx2
        read    (lunin,  *,  err=503)           xn1, xk1
        read    (lunin,  *,  err=503)           xn2, xk2
        read    (lunin,'(a)',err=503)           databuf(1:1)
        if (modebat.eq.0) write (*,27) nx1,nx2
  27    format  (' Number of points:  nx1=',i4,'  nx2=',i4)
c Array overflow:
        if (nx1.gt.ndim)        goto 30

        ymax =-1.E+37
        ymin = 1.E+37
        ymip = 1.E+37
        lest = Len(databuf)
        if (modebat.eq.0) write (*,*) 'Searching for [Min/Max]...'
        do j=1,nx2  !==========================================+
                                                              !|
          read (lunin,'(a)',err=503) databuf(1:lest)          !|
                                                              !|
          if (j.eq.1)   Then    !---------------------------+  |
c Determine the number of columns:                          |  |
            lest = Min0(Len(databuf),Len_Trim(databuf)+50) !|  |
c Column number mismatch:                                  !|  |
            ncc  = n_Column (databuf(1:lest))              !|  |
            if (ncc.ne.nx1)  goto 506  !--------------------+--+------+
          endif  !------------------------------------------+  |      v
                                                              !|
          Call  RDreal8 (yr8,nx1,databuf(1:lest),ifail)       !|
          if (ifail.ne.0)       goto 503                      !|
          do    i=1,nx1  !======================+              |
            if (yr8(i).gt.ymax) ymax = yr8(i)  !|max number    |
            if (yr8(i).lt.ymin) ymin = yr8(i)  !|mim number    |
            if (yr8(i).lt.ymip .and.           !|              |
     *          yr8(i).gt.0.)   ymip = yr8(i)  !|min positive  |
          enddo  !==============================+    number    |
        enddo  !===============================================+
        Rewind (unit=lunin)
c Skip the header:
        do      j=1,5    !===========================+
          read  (lunin,'(a)',err=503) databuf(1:1)  !|
        enddo  !=====================================+

        if (modebat.eq.0) Then  !---------------------+
          write (*,*)                                !|
          write (*,*) 'Min = ',ymin,'  Max = ',ymax  !|
        endif  !--------------------------------------+

        if (abs(ymax-ymin).lt.1.E-32) Then !---------------------+
          if (modebat.eq.0) Then  !---------------------------+  |
            write (*,*)                                      !|  |
            write (*,*) '*** Nothing to plot: Min=Max!'      !|  |
            write (*,*) '*** GNU file will not be produced!' !|  |
          endif  !--------------------------------------------+  |
          ifail = -1                                            !|
          close(lunin,err=100)                                  !|
          goto 100                                              !|
        endif  !-------------------------------------------------+

        if (ymax .le. 0.) Then !---------------------------------------------+
          if (Logscale) Then !---------------------------------------------+ |
            if (modebat.eq.0) Then  !------------------------------------+ | |
              write (*,*)                                               !| | |
              write (*,*) '!!! No positive data: turning off logscale!' !| | |
            endif  !-----------------------------------------------------+ | |
            Logscale = .false.                                            !| |
            ymip = ymin                                                   !| |
          endif  !---------------------------------------------------------+ |
        endif  !-------------------------------------------------------------+

        if (LogScale) Then  !------------------------------------------+
c Threshold correction                                                 |
          if (ymip.gt.0. .and. Threshold.gt.0.) Then !--+              |
            if (modebat.eq.0) Then !----------------+   |              |
              write (*,*) 'Range = ',ymax/ymip,    !|   |              |
     *                  '  Threshold range = ',    !|   |              |
     *                                  Threshold  !|   |              |
            endif ! --------------------------------+   |              |
            if ((ymax/ymip) .gt. Threshold) Then !--+   |              |
              ymip = ymax / Threshold              !|   |              |
            endif  !--------------------------------+   |              |
            if (modebat.eq.0) Then !----------------+   |              |
              write (*,*) 'Min_cutoff = ',ymip     !|   |              |
            endif ! --------------------------------+   |              |
          endif  !--------------------------------------+              |
                                                                      !|
          if (modebat.eq.0) Then !---+                                 |
            write (*,*)             !|                                 |
            write (*,*)             !|                                 |
          endif  !-------------------+                                 |
                                                                      !|
c Cannot normalize plain files:                                        |
        if (abs(ymax-ymip).lt.1.E-32) Normalize=.false.               !|
                                                                      !|
c Normalization of LOG file:                                           |
          if (Normalize) Then !----------------------------------+     |
c We map the whole LOG range into [0:10] interval:               |     |
            ynorm = 10.                                         !|     |
c This is 10^ynorm-1 :                                           |     |
            y8 = dexp(dble(ynorm)*dlog(10.0D0)) - 1.0D0         !|     |
            if (modebat.eq.0) Then !---------------------------+ |     |
              write (*,*)                                     !| |     |
              write (*,*) ' Output file normalized to ',ynorm !| |     |
              write (*,*)                                     !| |     |
            endif  !-------------------------------------------+ |     |
          else  !------------------------------------------------|     |
            ynorm = 0.                                          !|     |
            y8 = 0.0D0                                          !|     |
          endif !------------------------------------------------+     |
                                                                      !|
        endif  !-------------------------------------------------------+

        if (modebat.eq.0) write (*,*)

        outname = inpname(1:ln-4)//'.gnu'

        if (modebat.eq.0) Then  !-----+
          write (*,*)                !|
          write (*,*)                !|
          write (*,21)  outname      !|
        endif !-----------------------+

c See serve/openbin.for (open a new binary file
c with the test for "file exists"). Whether the
c program is interactive or not, depends on the
c modebat flag transferred via common /batmode/:
        Call  OpenBIN  (outname(1:ln),lunout,ifail)
        if (ifail.ne.0) goto 502

        if (nx1.gt.1)   then  !-------+
          dx1 = (xk1 - xn1)/(nx1-1)  !|
        else  !-----------------------|
          dx1 = 0.                   !|
        endif  !----------------------+

        if (nx2.gt.1)   then  !-------+
          dx2 = (xk2 - xn2)/(nx2-1)  !|
        else  !-----------------------|
          dx2 = 0.                   !|
        endif  !----------------------+

        do    i=1,nx1  !==============+
          x1    = xn1 + dx1*(i-1)    !|
          yr4(i) = x1                !|
        enddo  !======================+
        nx1_r = nx1
        write(lunout,err=505) nx1_r, (yr4(i),i=1,nx1)

        lest = Len(databuf)
        do j=1,nx2  !============================================+
          x2  =  xn2 + dx2*(j-1)                                !|
                                                                !|
          read  (lunin,'(a)',err=503) databuf(1:lest)           !|
                                                                !|
          if (j.eq.1)   Then    !---------------------------+    |
c Determine the number of columns:                          |    |
            lest = Min0(Len(databuf),Len_Trim(databuf)+50) !|    |
            ncc  = n_Column (databuf(1:lest))              !|    |
c Column number mismatch:                                  !|    |
            if (ncc.ne.nx1)  goto 506 !---------------------+----+----+
          endif  !------------------------------------------+    |    v
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
              yr4(i) = sngl(dlog10(yr8(i)))          !| |        |
            enddo  !==================================+ |        |
          endif  !--------------------------------------+        |
                                                                !|
          write(lunout,err=505) x2, (yr4(i),i=1,nx1)            !|
        enddo  !=================================================+
        close (lunin,err=99)
  99    continue
        close (lunout,err=100)
        if (modebat.eq.0) Then !-----------------------------+
          write (*,*)                                       !|
          write (*,*) "File ",outname(1:ln)," is created."  !|
        endif !----------------------------------------------+

  100   continue
        if (modebat.eq.0) write (*,*)
        if (ifail .eq. 0) Then !-----+
          Call exit_quiet()         !|
        else !-----------------------|
          stop 1                    !|
        endif  !---------------------+

c=======================================================
c                     ERRORS PROCESSING:
c-------------------------------------------------------

  30    continue
        txt(3) = 'File = '//inpname(1:ln)
        write   (txt(4),314) nx1,nx2
  314   format  ('Too many points in file:', 2i8)
        lines  = 4
        goto    998

  501   continue
        txt(3) = 'File = '//inpname(1:ln)
        txt(4) = 'Input file opening error!'
        lines  = 4
        goto    998

  502   continue
        txt(3) = 'File = '//outname(1:ln)
        txt(4) = 'Output file opening error (file already exists?)'
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
  171   format('GRD2GNU  E R R O R:'//
     *         'File = ',a/
     *         'Columns mismatch in input file:'/
     *         'Claimed  = ',i7/
     *         'Found    = ',i7)
        lines  = 6
        goto    999

  998   continue
        txt(1) = 'GRD2GNU  E R R O R'
        txt(2) = ' '
  999   continue
        Call    Message (txt,lines,2)
        ifail = -1
        goto 100

        end
