c ----------------------------------------------------------
c This is set of functions for calculating and then
c interpolating correlation integral. Contents:
c 1. Cos_Integral (p)
c 2. Read_CoFu_TBL (h,lun)
c 3. Calc_n12h (h)
c ----------------------------------------------------------

c =======================================================

        Real*4  Function  Cos_Integral (p)
c+=========================================================+
c|This function is called from trds_xx and gids_xx programs|
c+=========================================================+
        Real*4     p, q, p1, p2, A1, A2, CI
        Integer    n1, n2

        Real*4     Finli
        External   Finli

        Real*4     ArrCos(4001),   !Table of G_n (see Sinha)
     *             RangeCos(2),    !Range of argument (0.:10.)
     *             StepCos         !(Range(2)-Range(1))/(npoints-1)
        Integer*2  nptCos
        Common /ArrcoDat/ ArrCos,RangeCos,StepCos,nptCos

c The upper estimation for p is: p=qx*L. For qx=0.004
c and L=10,000~A  we get: p=40.

        q = Abs(p)

        if (q.le.RangeCos(2))      Then  !--------------+
          Cos_Integral = Finli (q,ArrCos,RangeCos(1),  !|
     *                                StepCos,nptCos)  !|
        else  !-----------------------------------------|
c cccccc  Cos_Integral = 0.                            !|
c cccccc  n1 = 0.8*nptCos                              !|
          n1 = nptCos-1                                !|
          n2 = nptCos                                  !|
          p1 = RangeCos(1) + StepCos*(n1-1)            !|
          p2 = RangeCos(1) + StepCos*(n2-1)            !|
          A1 = ArrCos(n1)                              !|
          A2 = ArrCos(n2)                              !|
          CI = A1 + (A2-A1)*(q-p1)/(p2-p1)             !|
c Avoid possible noise:                                !|
          if (CI.gt.Min(A1,A2)) CI=Min(A1,A2)          !|
          if (CI.lt.0.)         CI=0.                  !|
          Cos_Integral = CI                            !|
        endif  !----------------------------------------+

        return
        end

c =======================================================

        Subroutine Read_CoFu_TBL (h,lun)
c+===========================================================+
c|This subroutine is called from trds_xx and gids_xx programs|
c+===========================================================+
c It reads the following integral previously tabulated
c by CoFu_Enc executable:
c+-----------------------------------------------------------+
c|     +inf                   | -- for p={0.,50.}            |
c|      ?    -x^{2h}          | The details see in S.K.Sinha |
c|G(p)= |dx e      cos(px)    |  -- J.Physique III (France)  |
c|      ?                     | v.4, p.1543-1557 (1994)      |
c|      0                     |  --- Eq.(25) on p.1552.      |
c+-----------------------------------------------------------+
c
c ATTENTION: for alternative (old or unsuccessful) ideas about
c            calculating this integral see UTILITY/CORRFUN0.FOR
c            and  UTILITY/COFFINT8.FOR.
c -------------------------------------------------------
        Integer    ndim
        Parameter  (ndim= 4001)
        Real*4     h,
     *             RangeH(2),           !hmin,hmax
     *             dh,                  !(hmax-hmin)/(nh-1)
     *             x, q, Z
        Integer    i, j, ih,
     *             lun, lug, lwrk,
     *             istat, ifail
        Integer*2  nh
        Character  CoFu_TBL*12 /'corrfun.cft'/

        Character       getOS*8
        External        getOS

        Logical*4  FileExist
        External   FileExist

        Real*4     ArrCos(4001),   !Table of G_n (see Sinha)
     *             RangeCos(2),    !pmin,pmax
     *             StepCos         !(pmax-pmin)/(nptCos-1)
        Integer*2  nptCos
        Common  /ArrcoDat/ ArrCos,RangeCos,StepCos,nptCos

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c x0hpath*80  -  the path to the X)H directory where all the
c                DB file including corrfun.cft are located
c lx0hpat     -  the number of characters in the path
        Integer         lx0hpat
        Character       x0hpath*80
        common  /x0hpa/ x0hpath,lx0hpat
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'read_cofu_tbl'      !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
c -------------------------------------------------------

c Calculate n^{1/2h}  for n=[1-20]:
        Call Calc_n12h (h)

c -------------------------------------------------------
c The CFT file is supposed to be either in X0H directory,
c or in the directory of EXE file. Find this path in the
c case if it has not been determined by X0h software:
        if (x0hpath(1:1).eq.char(0)) Then !------------------------+
                                                                  !|
          Call Defpat (x0hpath)                                   !|
          i = Len_Trim (x0hpath)                                  !|
          lx0hpat = i                                             !|
                                                                  !|
          x0hpath(i+1:80) = CoFu_TBL                              !|
          if (.NOT.FileExist(x0hpath)) Then  !---------------+     |
                                                            !|     |
            Call  GetEnv32 ('X0H',txt(20))                  !|     |
            i = Len_Trim (txt(20))                          !|     |
            if (i.eq.0)   goto 5   !-------------------------+--+  |
                                                            !|  v  |
            x0hpath = txt(20)                               !|     |
            do j=1,i  !==================================+   |     |
              if (x0hpath(j:j).eq.'\') x0hpath(j:j)='/' !|   |     |
            enddo !======================================+   |     |
            if (x0hpath(i:i).ne.'\' .AND.                   !|     |
     *          x0hpath(i:i).ne.'/') Then !---------+        |     |
              i = i+1                              !|        |     |
              if (getOS() .eq. 'windows') Then !-+  |        |     |
c               x0hpath(i:i)='\'                !|  |        |     |
                x0hpath(i:i)='/'                !|  |        |     |
              else !-----------------------------+  |        |     |
                x0hpath(i:i)='/'                !|  |        |     |
              endif !----------------------------+  |        |     |
            endif  !--------------------------------+        |     |
            lx0hpat = i                                     !|     |
                                                            !|     |
            x0hpath(i+1:80) = CoFu_TBL                      !|     |
            if (.NOT.FileExist(x0hpath)) goto 5 !------------+--+  |
                                                            !|  v  |
          endif  !-------------------------------------------+     |
                                                                  !|
        endif  !---------------------------------------------------+

        x0hpath(lx0hpat+1:80) = CoFu_TBL
        lwrk = Len_Trim(x0hpath)

        if (lun.gt.0) Then !-+
          lug = lun         !|
        else !---------------|
          lug = 88          !|
        endif  !-------------+
c Open correlation function table (CFT):
        Call OpenBinary (x0hpath(1:lwrk),lug,'old','read',istat,ifail)
        if (ifail.ne.0) goto 5

        read    (lug,err=6)     nptCos,nh,
     *                          RangeCos(1),RangeCos(2),
     *                          RangeH(1),  RangeH(2)

        if (nptCos.ne.ndim) goto 7

        dh      = (RangeH(2)-RangeH(1)) / (nh-1)
        stepCos = (RangeCos(2)-RangeCos(1)) / (nptCos-1)

        x = (h-RangeH(1))/dh + 1.
        ih = INT(x+0.1)
        if (ih.lt.1)    ih=1
        if (ih.gt.nh-1) ih=nh-1
        q  = x-ih

        do      i=1,ih  !============================+
          read (lug,err=6) (ArrCos(j),j=1,nptCos)   !|
        enddo  !=====================================+
        do      j=1,nptCos  !========================+
          read (lug,err=6,end=9) Z                  !|
          ArrCos(j) = ArrCos(j) + q*(Z-ArrCos(j))   !|
                                                    !|
c cccc    p  = RangeCos(1) + stepCos*(j-1)          !|
c cccc    write (3,*) p, ArrCos(j)                  !|
        enddo  !=====================================+
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return

c#######################################################
c                  Error messages
c########################################################
  5     continue
        txt(1) = 'Cannot open'
        goto 28

  6     continue
        txt(1) = 'Error while reading'
        goto 28

  9     continue
        txt(1) = 'End-of-file while reading'
        goto 28

  7     continue
        write (txt(1),17) nptCos, ndim
  17    format('The table length [',i4,'] does '
     *          'not match buffer [',i4,']')
        goto 28

  28    continue
        txt(2) = ' '
        txt(3) = x0hpath
        Call    Message (txt,3,2)
        Call exit_quiet()
        End

c =======================================================

        Subroutine Calc_n12h (h)
c Calculate n**(1/2h)  for n=[1-20]:
        Real*4     h, x
        Integer    n

        Real*4     n12h(20)        ! n^{1/2h}
        Common  /n12/     n12h

        n12h(1) = 1.
        x       = 1./(2.*h)
        do      n=2,20   !===============+
          n12h(n) = Exp(x*Log(Real(n))) !|n**(1/2h)
        enddo  !=========================+
        return
        end
