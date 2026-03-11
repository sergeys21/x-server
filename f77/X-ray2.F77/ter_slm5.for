        Program Ter_slm5
c       USE MSFLIB                                      !for QuickWin
c+-------------------------------------------------------------------+
c|Nearly same as ter_sl97/ter_sl99 (same algorithm). The difference  |
c|from ter_sl99 is that ter_slm5 uses new engine SFMm5 returning the |
c|the x-ray wavefields in the target.                                |
c+-------------------------------------------------------------------+
c -------------------------------------------------------
c   X-ray total external reflection from multilayers
c
c                     By S.A.Stepanov
c
c      *** Version'm5 optimized for WWW access ***
c
c 1.0a Transition layers PROFILES included,
c      Small-Scale roughness PROFILES included.
c 1.2  Simplified: SL excluded
c      (One general profile for:
c      thickness, x0*exp(-w), transition, & roughness)
c '97  Combined advantages of 1.0 and 1.2: script language for
c      top-layer profile, allowing a specification of periodic
c      structures; improved batch ("silent" or "server") mode
c      for using with CGI-interfaces.
c '99  Same as '97 but:
c      1. added x-ray polarization dependence
c      2. added possibility to specify rho and chemical formulas
c      3. added "Identical" flag to accelerate the calculation of x0.
c '99.2 Added Henke/Cowan option
c 'm5   Added option to calclulate maps of wavefields.
c -------------------------------------------------------
c+===============================+
        Include 'ter_slm5.inc'  !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter       (N_Top_Max     = 1001,
c cc *                   N_Total_Max   = N_Top_Max+1)
c cc    Parameter       (N_Tran_Max    = 101)           !50
c cc    Parameter       (N_Stand_Max   = 201)           !50

        Integer          N_Frac_Max
        Parameter       (N_Frac_Max    = 4)

        Complex*8       x0(N_Total_Max+1), x0_, x0_f

        Real*8          k0,                     ! = 2*pi/wave
     *                  R0s,                    !Reflection coeff.
     *                  SFMm5                   !calc. engine
        External        SFMm5

        Real*4          Energy, rpar(1),
     *                  Thickness_Top,
     *                  f0, f0_, f0_range(2), df0,
     *                  TER_Angle, rho(N_Total_Max),
     *                  W0(N_Total_Max),
     *                  Frac(N_Frac_Max,N_Total_Max),
     *                  pi, AngDegree, AngMinute,
     *                  AngSec, AngMrad, AngUrad,
     *                  xr0, xi0, rho_, rho_f,
     *                  qz_max, f0_max, Sigma_max, SQZ/10./,
     *                  Trans_mean, Sigma_mean, W0_mean,
     *                  z, DDD, xxx, Start, progress

        Integer*4       hh,mm,ss,bias

        Integer         Law_Index, nf0, iUnits, iBeepStep,
     *                  N_Frac(N_Total_Max), Identical(N_Total_Max),
     *                  ipar(1), i, j, l, linp, lout, ltop,
     *                  ipv, line, io_status,
     *                  iBatch_mode, ixway, iPol, mode_x0h,
     *                  iSyngony, nf, nerror,
     *                  iscan, iasci, iii, luninp/2/, narg

        Integer         iarg

        Character       Code(N_Frac_Max,N_Total_Max)*20,
     *                  C00(N_Total_Max)*2,
     *                  Radiat*6, wrk*80,
     *                  InpFile*80, OutFile*80,
     *                  File_Top*80, DoneFile*80,
     *                  ProgramVer*8, Version(2)*80,
     *                  timestamp*9, Comment(3)*80,
     *                  YesNo(2)*3 /
     *                           'No',
     *                           'Yes'/,
     *                  Spm(2)*5 /
     *                           'Sigma',
     *                           ' Pi  '/,
     *                  LawName(2)*6 /
     *                           '[cos]',
     *                           'Linear'/,
     *                  Uni(5)*5 /
     *                           'degr.',
     *                           'min.',
     *                           'mrad.',
     *                           'sec.',
     *                           'urad.'/

c X0h (International Tables) range:
c             (5,000 eV --  25,000 eV):   0.50  --    2.47 Angstrem
c Henke    range (10 eV --  30,000 eV):   0.41  -- 1239.81 Angstrem
c Cowan    range (30 eV -- 694,500 eV):   0.02  --  413.27 Angstrem
c Windt    range (10 eV -- 100,000 eV):   0.12  -- 1239.81 Angstrem
c Chantler range (10 eV -- 450,000 eV):   0.28  -- 1239.81 Angstrem
        Character DBtext(10)*45 /
     *                  'Automatic DB selection',                        !-1 (1)
     *                  'X0h (International Tables), 5-25 KeV',          !0  (2)
     *                  'Henke Tables, 0.01-30 KeV (df1 only)',          !1  (3)
     *                  'Henke Tables, 0.01-30 KeV',                     !2  (4)
     *                  'Brennan-Cowan Tables, 0.03-700 KeV (df1 only)', !3  (5)
     *                  'Brennan-Cowan Tables, 0.03-700 KeV',            !4  (6)
     *                  'Windt Tables, 0.01-100 KeV (df1 only)',         !5  (7)
     *                  'Windt Tables, 0.01-100 KeV',                    !6  (8)
     *                  'Chantler/NIST Tables, 0.01-450 KeV (df1 only)', !7  (9)
     *                  'Chantler/NIST Tables, 0.01-450 KeV'/            !8  (10)

        Real            wave2energy
        External        wave2energy

        Logical*4       FileExist
        External        FileExist

        Logical         YesTr

        Integer*4       nargs
        External        nargs

        Real*8          SW(N_Stand_Max),
     *                  Positions(N_Total_Max)
        Real*4          Wave, UnitCoef,
     *                  Thickness(N_Total_Max),
     *                  Thickness_Tr(N_Total_Max),
     *                  TranLaw(N_Tran_Max),
     *                  Sigma(N_Total_Max),
     *                  xabs,                            ! abs(x0_max)
     *                  TERang,                          ! critical TER angle (radian)
     *                  standing_range(2), standing_step
        Integer         N_Top, N_Total, N_Used, N_Tran,
     *                  MaxOrder, iDebug, ifail,
     *                  i_standing, i_reference, n_standing
        Common /TER_m5/ Wave, UnitCoef, x0,
     *                  Positions,
     *                  Thickness,
     *                  Thickness_Tr,
     *                  TranLaw,
     *                  Sigma,
     *                  xabs, TERang,
     *                  standing_range, standing_step,
     *                  SW,
     *                  N_Top, N_Total, N_Used, N_Tran,
     *                  MaxOrder, iDebug, ifail,
     *                  i_standing, i_reference, n_standing

c       Parameter       (kcompMax = 10)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Real            prcn(kcompMax)
        Common  /x0pa3/ prcn

c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer              iHenkeCowan
        Common  /HenkeCowan/ iHenkeCowan

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

        istackrezv  = 0
        progname    = 'ter_slm5'
c       i           = SetExitQQ(2)        !QWIN$EXITNOPERSIST=2   !for QuickWin
c -------------------------------
        i = N_Top_Max                   !to suppress G77 warnings
        i = N_Total_Max                 !to suppress G77 warnings
        i = N_Tran_Max                  !to suppress G77 warnings
        i = N_Pts_Max                   !to suppress G77 warnings
        i = N_Stand_Max                 !to suppress G77 warnings

        i = kcompMax                    !to suppress G77 warnings
        i = kolMax                      !to suppress G77 warnings
        i = ncodMax                     !to suppress G77 warnings
        i = natmMax                     !to suppress G77 warnings
        i = nwavMax                     !to suppress G77 warnings
        i = nedgMax                     !to suppress G77 warnings
        i = nabcMax                     !to suppress G77 warnings
        i = nedgTb                      !to suppress G77 warnings
        i = maxX0hWarnLines             !to suppress G77 warnings
c -------------------------------
        ifail       = 0
        pi          = 4.*atan(1.)
        AngDegree   = 2.*pi/360.
        AngMinute   = AngDegree/60.
        AngSec      = AngDegree/3600.
        AngMrad     = 1.E-03
        AngUrad     = 1.E-06
        i_standing  = 0
        i_reference = 0
        n_standing  = 0
        do i=1,N_Stand_Max !==+
          SW(i) = 0.         !|
        enddo !===============+

        do      l=1,N_Total_Max  !====+
          N_Frac(l) = 1              !|
          do    j=1,N_Frac_Max  !==+  |
            Frac(j,l) = 0.        !|  |
            Code(j,l) = ' '       !|  |
          enddo  !=================+  |
          Frac(1,l) = 1.             !|
        enddo  !======================+

        ProgramVer = 'ter_slm5'
        ipv = Len_Trim (ProgramVer)
c This is the startup ERR-file (it is used until we've
c verified the output filename):
        ErrorFile = ProgramVer(1:ipv)//'.err'
        Call    DeleteErrorFile ()

        write   (Version,111)
  111   format  ('TER_sl: X-Ray Total External Reflection ',
     *  'from Multilayers with Rough Interfaces.'/
     *  'By S.Stepanov <sstepanov@anl.gov>',24x,
     *  'Version-M5: May, 2005')

        narg = Nargs()-1
        if (narg .eq. 0) Then !---------------------+
                                                   !|
          InpFile = ProgramVer(1:ipv)//'.inp'      !|
                                                   !|
        else  !-------------------------------------+
                                                   !|
          iarg = 1                                 !|
          Call get_arg_wrapper (iarg,InpFile)      !|
                                                   !|
        endif  !------------------------------------+

        modebat = 1                     !assume batch-mode by default
        iBatch_mode = modebat

c #######################################################
c              Read input file
c#######################################################
        linp = Len_Trim(InpFile)
        Call OpenFile(InpFile(1:linp),luninp,
     *                'read','old',io_status,*99)

        Call Make_Err_Filename (InpFile(1:linp))

        line = 0
        txt(7) = 'batch mode flag [0/1/2]'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        iBatch_mode = ipar(1)
        if (iBatch_mode .lt. 0)   goto 101
        if (iBatch_mode .gt. 2)   goto 101
        if (iBatch_mode .eq. 0) iBatch_mode=1
        modebat = iBatch_mode
c-----------------------
        txt(7) = 'output file name'
        Call    LineLoop (luninp,line,wrk,*100)
        OutFile = wrk
        lout = Len_Trim(OutFile)
        if (lout .eq. 0) Then !-----+
          OutFile = '~'            !|
          lout    = 1              !|
        endif  !--------------------+
c       Call    Case    (OutFile(1:lout),1)     !to lower case
        Call    OpenList(OutFile(1:lout)//'.tbl',3,ifail)
        if (ifail .ne. 0)       goto 28
c ccc   lunbat = 3                !this is for ERR file (default is 252)
        Call    OpenList (OutFile(1:lout)//'.dat',1,ifail)
        if (ifail .ne. 0)       goto 28
c-----------------------
        txt(7) = 'surface profile name'
        N_Top   = 0
        N_Total = 0
        Call    LineLoop (luninp,line,wrk,*100)
        File_Top = wrk
        ltop = Max (Len_Trim(File_Top),1)
c       Call    Case (File_Top(1:ltop),1)       !to lower case
c-----------------------
        do i=1,3  !=====================================+
          write (txt(7),'(a,i1)') 'comment line #',i   !|
          Call  LineLoop (luninp,line,Comment(i),*100) !|
        enddo  !========================================+
c-----------------------
        txt(7) = 'substrate data'
        Thickness(N_Total_Max) = 0.
        W0(N_Total_Max)        = 1.
c Set transtions to "yes" while we do not know if we need standing waves:
        YesTr = .True.
        Call Read_Sub99_0 (Code(1,N_Total_Max), rho(N_Total_Max),
     *                     C00(N_Total_Max),
     *                     x0(N_Total_Max+1), W0(N_Total_Max),
     *                     Sigma(N_Total_Max),
     *                     Thickness_Tr(N_Total_Max),YesTr,
     *                     wrk, line, *100, *101, *28)
c-----------------------
        txt(7) = 'x-ray specification mode (1=wave 2=energy 3=line)'
        Call    LineLoop (luninp,line,wrk,*100)
c 1=wave 2=energy 3=x-ray line:
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        ixway = ipar(1)
        if (ixway .lt. 1)       goto 101
        if (ixway .gt. 3)       goto 101
c-----------------------
        txt(7) = 'x-ray wavelength/energy (range=[0.01-1000.0]Angstrom)'
        Call    LineLoop (luninp,line,wrk,*100)
        Wave = 0.
        if (ixway.lt.3) Then  !-----------------------+
          Call  rdReal  (rpar,1,wrk,i)               !|
          if (i .ne. 0)            goto 100          !|
          Wave = rpar(1)                             !|
          if (Wave .le. 0.)        goto 101          !|
            if (ixway .eq. 2) Wave=wave2energy(Wave) !|
            if (Wave .lt. 0.01)    goto 101          !|
            if (Wave .gt. 1000.)   goto 101          !|
        endif  !--------------------------------------+
c-----------------------
        txt(7) = 'x-ray line name'
        Call    LineLoop (luninp,line,wrk,*100)
        if (ixway.eq.3) Then  !-------------------+
          Radiat = wrk(1:6)                      !|
          if (Len_Trim(Radiat).eq.0) goto 131 !---+---+
c Get wavelength from the X-ray line name:        |   v
          Call  Get_Wave (Wave,Radiat,ifail)     !|
          if (ifail.ne.0)          goto 130      !|
          if (Wave .lt. 0.01)      goto 101      !|
          if (Wave .gt. 1000.)     goto 101      !|
        else !------------------------------------|
          Radiat = ' '                           !|
        endif  !----------------------------------+
        k0 = 2.* pi/wave
c-----------------------
        txt(7) = 'x-ray polarization (0,1=sigma 2=pi)'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        iPol = ipar(1)
c Polarization: 0-mixed, 1-sigma, 2-pi:
        if (iPol .eq. 0)        iPol=1
        if (iPol .lt. 1)        goto 101
        if (iPol .gt. 2)        goto 101
c-----------------------
        write (txt(7),'(a,i3,a)') 'number of transition sublayers (1-',
     *                          N_Tran_Max,')'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)               goto 100
        N_Tran = ipar(1)
        if (N_Tran .lt. 1) N_Tran = 1
        if (N_Tran .gt. N_Tran_Max) goto 101
c-----------------------
        txt(7) = 'transition layers type (0=cosine, 1=linear)'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        Law_Index = ipar(1)
        if (Law_Index .lt. 0)   goto 101
        if (Law_Index .gt. 1)   goto 101
c-----------------------
        txt(7) = 'angular units (0=degr,1=min,2=mrad,3=sec,4=urad)'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        iUnits = ipar(1)
        if (iUnits*(4-iUnits) .lt. 0) goto 101
        goto (40,41,42,43,44) iUnits+1 !->-+
  40    continue                !<---------|
        UnitCoef=AngDegree                !|
        goto 49   !->-------------------->-+--->+
  41    continue                !<---------|    |
        UnitCoef = AngMinute              !|    |
        goto 49   !->-------------------->-+--->|
  42    continue                !<---------|    |
        UnitCoef = AngMrad                !|    |
        goto 49   !->-------------------->-+--->|
  43    continue                !<---------|    |
        UnitCoef = AngSec                 !|    |
        goto 49   !->-------------------->-+--->|
  44    continue                !<---------+    |
        UnitCoef = AngUrad                     !|
        goto 49   !->-------------------->----->|
  49    continue   !<--------------<-----------<+
c-----------------------
        txt(7) = 'incidence angle limits (range=[0:90]degr.)'
        Call    LineLoop (luninp,line,wrk,*100)
c Read incidence angle range:
        Call    Rdreal  (f0_range,2,wrk,i)
        if (i .ne. 0)           goto 100
        if (f0_range(1) .lt. 0.)                        goto 101
        if (f0_range(2) .lt. 0.)                        goto 101
        if (UnitCoef*f0_range(1) .gt. 90.001*AngDegree) goto 101
        if (UnitCoef*f0_range(2) .gt. 90.001*AngDegree) goto 101
c-----------------------
        write (txt(7),'(2a,i5,a)') 'number of incidence angle points ',
     *                             '(min=1, max=',N_Pts_Max,')'
        Call    LineLoop (luninp,line,wrk,*100)
c Read incidence angle points:
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        nf0 = ipar(1)
        if (nf0 .lt. 1)         goto 101
        if (nf0 .gt. N_Pts_Max) goto 101
        if (nf0 .gt. 1) Then !-------------------------+
          xxx = f0_range(2) - f0_range(1)             !|
          df0 = xxx / (nf0-1.)                        !|
          if (abs(xxx).lt.1.E-20) goto 109            !|
        else  !----------------------------------------|
          df0 = 0.                                    !|
        endif  !---------------------------------------+
c-----------------------
        txt(7) = 'standing waves printout option (0=no, 1=yes)'
        Call    LineLoop (luninp,line,wrk,*100)
c Read standing waves option:
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        i_standing = ipar(1)
        if (i_standing .lt. 0)  goto 101
        if (i_standing .gt. 1)  goto 101
c-----------------------
        txt(7) = 'standing waves reference interface (0=surface)'
        Call    LineLoop (luninp,line,wrk,*100)
c Read standing waves option:
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        i_reference = ipar(1)
        if (i_reference .lt. 0) goto 101
c-----------------------
        txt(7) = 'standing waves offsets from given interface [from,to]'
        Call    LineLoop (luninp,line,wrk,*100)
        Call  rdReal  (standing_range,2,wrk,i)
        if (i .ne. 0)           goto 100
c-----------------------
        write (txt(7),'(a,i3,a)') 'number of standing wave points (1-',
     *                          N_Stand_Max,')'
        Call    LineLoop (luninp,line,wrk,*100)
c Read standing waves option:
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        n_standing = ipar(1)
        if (n_standing .eq. 0)  n_standing = 1
        if (n_standing .lt. 1)           goto 101
        if (n_standing .gt. N_Stand_Max) goto 101
        if (abs(standing_range(2)-standing_range(1)).lt.1.E-20)
     *                                   n_standing = 1
        if (n_standing .gt. 1) Then !-------------------------------+
          standing_step = (standing_range(2)-standing_range(1))    !|
     /                  / (n_standing-1.)                          !|
        else !------------------------------------------------------|
          standing_step = 0.                                       !|
        endif !-----------------------------------------------------+
c-----------------------
        MaxOrder  = 14
        iDebug    = 0
        iBeepStep = 10
        close   (unit=luninp)
c==============================================================
c Read top layer profile:
        YesTr = .True.
        txt(7) = 'surface layer profile'
c See: Rd99prf0.for:
        Call    Read_Top99_0 (File_Top(1:ltop),
     *                        N_Top, N_Top_Max,
     *                        Thickness, Thickness_Top,
     *                        Code(1,1), Frac(1,1),
     *                        N_Frac(1), N_Frac_Max,
     *                        rho(1), C00(1), x0(2), W0(1),
     *                        Identical(1),
     *                        Sigma(1),
     *                        Thickness_Tr(1), YesTr,
     *                        ifail)
        if (ifail .ne. 0)       goto 28
        N_Total = N_Top+1
c-----------------------
c Copy the substrate data that was previously placed into N_Total_Max:
        if (N_Total .lt. N_Total_Max) Then !---------------------+
          Thickness(N_Total)    = Thickness(N_Total_Max)        !|
          Thickness_Tr(N_Total) = Thickness_Tr(N_Total_Max)     !|
          Sigma(N_Total)        = Sigma(N_Total_Max)            !|
          x0(N_Total+1)         = x0(N_Total_Max+1)             !|
          W0(N_Total)           = W0(N_Total_Max)               !|
          C00(N_Total)          = C00(N_Total_Max)              !|
          rho(N_Total)          = rho(N_Total_Max)              !|
          N_Frac(N_Total)       = N_Frac(N_Total_Max)           !|
          do i=1,N_Frac_Max !==============================+     |
            Code(i,N_Total)     = Code(i,N_Total_Max)     !|     |
            Frac(i,N_Total)     = Frac(i,N_Total_Max)     !|     |
          enddo  !=========================================+     |
c Now for any case wipe the old substrate place:                 |
          Thickness(N_Total_Max)    = 0.                        !|
          Thickness_Tr(N_Total_Max) = 0.                        !|
          Sigma(N_Total_Max)        = 0.                        !|
          x0(N_Total_Max+1)         = (0.,0.)                   !|
          W0(N_Total_Max)           = 0.                        !|
          C00(N_Total_Max)          = ' '                       !|
          rho(N_Total_Max)          = 0.                        !|
          N_Frac(N_Total_Max)       = 0                         !|
          do i=1,N_Frac_Max !=================+                  |
            Code(i,N_Total_Max)     = ' '    !|                  |
            Frac(i,N_Total_Max)     = 0      !|                  |
          enddo  !============================+                  |
        endif !--------------------------------------------------+
c-----------------------
c Transition layers are only allowed if
c no standing waves requested:
        if (i_standing .ne. 0) Then !--------------------------+
          do i=1,N_Total  !=================================+  |
            if (abs(Thickness_Tr(i)).gt.1.E-20) goto 104 !--+--+---+
          enddo !===========================================+  |   v
        endif !------------------------------------------------+
c-----------------------
c Reference interface may not be larger than the number of interfaces:
        if (i_standing .ne. 0 .and.
     *      i_reference .gt. N_Total-1) goto 105
c-----------------------
c Get substrate x0:
        if (C00(N_Total) .ne. 'x')  Then !-------------------+
          mode_x0h  = 0                       !new crystal   |
          iSyngony  = 0                                     !|
          do i=1,kcompMax !==+                               |
            prcn(i) = 0.    !|                               |
          enddo  !===========+                               |
          Call Get_X0 (Wave, Code(1,N_Total), rho(N_Total), !|
     *                 iSyngony, mode_x0h, xr0, xi0,        !|
     *                 x0(N_Total+1), 1, 1, ifail)          !|
          if (ifail .ne. 0)     goto 130                    !|
        endif !----------------------------------------------+
        Identical(N_Total) = 0
c A control for the substrate x0:
        if (abs(Real(x0(N_Total+1))).lt.1.E-32 .OR.
     *              abs(W0(N_Total)).lt.1.E-32)   goto 123
c #######################################################
        do      l=1,N_Top  !=================================+
          if (Identical(l) .eq. 0) Then  !-----------------+ |
                                                          !| |
            if ((Len_Trim(Code(1,l)) .eq. 0) .OR.         !| |
     *          (Code(1,l) .eq. 'Unknown'))  Then !---+    | |
              Code(1,l) = Code(1,N_Total)            !|    | |
            endif  !----------------------------------+    | |
            x0_  = (0., 0.)                               !| |
            rho_ = 0.                                     !| |
            do j=1,N_Frac(l)  !==========================+ | |
              if (j .eq. 1)  Then !-+                   !| | |
                rho_f = rho(l)     !|                   !| | |
              else  !---------------|                   !| | |
                rho_f = 0.         !|                   !| | |
              endif !---------------+                   !| | |
              if (Code(j,l).eq.Code(1,N_Total) .AND.    !| | |
     *              abs(rho_f).lt.1.E-20) Then !-------+ | | |
                x0_f  = x0(N_Total+1)                 !| | | |
c This means, even when the substrate density was     !| | | |
c specified externally, it is assigned to the layer   !| | | |
c with unspecified roughness:                         !| | | |
                rho_f = rho(N_Total)                  !| | | |
              else !-----------------------------------+ | | |
c New material, but same wavelength:                   | | | |
                mode_x0h = 1                          !| | | |
                iSyngony = 0                          !| | | |
                do i=1,kcompMax !==+                   | | | |
                  prcn(i) = 0.    !|                   | | | |
                enddo  !===========+                   | | | |
                Call Get_X0 (Wave, Code(j,l), rho_f,  !| | | |
     *                       iSyngony, mode_x0h,      !| | | |
     *                       xr0, xi0,                !| | | |
     *                       x0_f, 1, 1, ifail)       !| | | |
                if (ifail .ne. 0)    goto 130         !| | | |
              endif !----------------------------------+ | | |
              x0_  = x0_  + x0_f *Frac(j,l)             !| | |
              rho_ = rho_ + rho_f*Frac(j,l)             !| | |
            enddo  !=====================================+ | |
                                                          !| |
            if (C00(l) .eq. ' ') Then  !----------+        | |
c x0 was not specified in profile:                |        | |
              x0(l+1) = x0_                      !|        | |
              rho(l)  = rho_                     !|        | |
            else  !-------------------------------|        | |
c x0 was specified, than Code is not used:        |        | |
              Code(1,l) = 'Unknown'              !|        | |
              Frac(1,l) = 1.                     !|        | |
              N_Frac(l) = 1                      !|        | |
            endif !-------------------------------+        | |
                                                          !| |
            x0(l+1) = x0(l+1)*W0(l)                       !| |
                                                          !| |
          else  !------------------------------------------| |
                                                          !| |
            i = Identical(l)                              !| |
            do j=1,N_Frac(i)  !========+                   | |
              Code(j,l) = Code(j,i)   !|                   | |
              Frac(j,l) = Frac(j,i)   !|                   | |
            enddo  !===================+                   | |
            N_Frac(l) = N_Frac(i)                         !| |
            rho(l)    = rho(i)                            !| |
            x0(l+1)   = x0(i+1)                           !| |
                                                          !| |
          endif !------------------------------------------+ |
                                                            !|
        enddo  !=============================================+
c Corrections for substrate:
        x0(N_Total+1) = x0(N_Total+1)*W0(N_Total)

c 2018.10: changed xabs from |x0_substr| to |x0_max|
c       xabs      = abs(x0(N_Total+1))
c       TERang    = sqrt(xabs)
c       TER_Angle = TERang / UnitCoef
c       if (Abs(xabs) .lt. (1.E-24)) xabs = 1.
        Call xabsmax (x0, N_Total_Max, N_Total, UnitCoef,
     *                xabs, TERang, TER_angle)

c Control for big sigma (huge exponents):
        f0_max = Max(Abs(f0_range(1)),Abs(f0_range(2)))
        if (f0_max .gt. 0.) Then  !----------------------+
c This should be changed after we add scanning over qz:  |
c (See mag_sl99 for an example)                          |
          qz_max    = SNGL(2.*k0*sin(f0_max*UnitCoef))  !|
          Sigma_max = SQZ/qz_max                        !|
          do l=1,N_Total  !=========================+    |
            if (Sigma(l) .gt. Sigma_Max) goto 701 !-+----+-+
          enddo  !==================================+    | v
        endif  !-----------------------------------------+

c Calculate the mean profile parameters:
        Trans_mean = 0.
        Sigma_mean = 0.
        W0_mean    = 0.
        if (N_Top .gt. 0) Then  !--------------------------+
          do    l = 1, N_Top  !========================+   |
            Trans_mean = Trans_mean + Thickness_Tr(l) !|   |
            Sigma_mean = Sigma_mean + Sigma(l)        !|   |
            W0_mean    = W0_mean    + W0(l)           !|   |
          enddo  !=====================================+   |
          Trans_mean = Trans_mean / N_Top                 !|
          Sigma_mean = Sigma_mean / N_Top                 !|
          W0_mean    = W0_mean    / N_Top                 !|
        endif  !-------------------------------------------+

c Calculate positions of interfaces (for standing waves):
        Positions(1) = 0.
        do l=2, N_Total  !================================+
          Positions(l) = Positions(l-1) + Thickness(l-1) !|
        enddo !===========================================+

        Energy = wave2energy(Wave)
        i = Min(lout,40)
        j = Min(ltop,40)
        write (txt,33) OutFile(1:i),
     *                 File_Top(1:j),
     *                 Thickness_Top, N_Top,
     *                 Trans_mean, Sigma_mean,
     *                 Code(1,N_Total), Thickness_Tr(N_Total),
     *                                         Sigma(N_Total),
     *                 Wave, Energy, Radiat,
     *                 Spm(iPol),
     *                 Uni(iUnits+1), f0_range, nf0,
     *                 N_Tran, LawName(1+Law_Index),
     *                 YesNo(1+i_standing),
     *                 i_reference, Positions(i_reference+1),
     *                 (standing_range(i),i=1,2),
     *                 n_standing
  33    format  (
     *  'Output files [.dat & .tbl] will be...|',a/
     *  'Top: Profile name....................|',a/
     *  'Top: Thickness(A), N-sublayers.......|',f10.1,i10/
     *  'Top: <TransLayer>, <roughness>.......|',2f10.2/
     *  'Substrate: Code, TransLayer,Roughness|',a,    2f7.2/
     *  'X-ray wavelength (A), energy (keV)...|',g12.7,1x,g12.7,1x,a/
     *  'X-ray polarization...................|',a/
     *  'Incidence angle range [',a,'], points|',2g12.3,i8/
     *  'Transition layers: N_sublayers, Shape|',i5,10x,a/
     *  'Print standing waves?................|',a/
     *  'Standing waves: reference interface..|',i5,10x,'(',f13.2,'A)'/
     *  'Standing waves: offsets range (A)....|',g12.5,',  ',g12.5/
     *  'Standing waves: number of depth pts..|',i5)

        line = 13
c #######################################################
c                 Open output files:
c#######################################################
        write   (3,55)  Version
        write   (3,11)
  11    format  (78('='))
        write   (3,55)  (txt(i)(1:Len_Trim(txt(i))),i=1,line)
  55    format  (20(a:/))
        Call    Algorithm_Type  (Version(1))
        l = Len_Trim (Version(1))
        i = Len_Trim (DBtext(iHenkeCowan+2))
        write   (3,62)     MaxOrder,
     *                     Version(1)(1:l),
     *                     DBtext(iHenkeCowan+2)(1:i)
  62    format  (
     *  'Maximum order of T1*..Tn:  10**k.....|',i3/
     *  'Algorithm type (Matrix / Recursive)..|',a/
     *  'Database used for dispersion correct.|',a)
        write   (3,11)
        write   (3,55)  ' *** User Comment: '
        do i=1,3  !==================================+
          l = Len_Trim (Comment(i))                 !|
          if (l.gt.0) write (3,55) Comment(i)(1:l)  !|
        enddo  !=====================================+
        write   (3,11)

        if ((Thickness_Tr(N_Total) .gt. 0.)  .OR.
     *                 (Trans_mean .gt. 0.)) Then !-----+
                                                       !|
            do  i=1,N_Tran  !=====================+     |
              z = Real(i) / Real(N_Tran+1)       !|     |
              if (Law_Index .eq. 0) Then !-----+  |     |
c Cosine-law: f(z)=0.5{1+cos(pi*z)}, z:[0:1].  |  |     |
c -> This is equivalent to: cos<(pi*z/2):      |  |     |
c f(z=0)=1, f(z=pi)=0.                         |  |     |
                TranLaw(i) = cos(pi*z/2.)**2  !|  |     |
              else  !--------------------------|  |     |
c Linear-law: f(z)=1-z, z:[0:1].               |  |     |
c f(z=0)=1, f(z=pi)=0.                         |  |     |
                TranLaw(i) = 1.-z             !|  |     |
              endif  !-------------------------+  |     |
            enddo  !==============================+     |
                                                       !|
          write (3,66)  (TranLaw(i),i=1,N_Tran)        !v
  66      format (/
     *    'Transition law (from upper to lower):'/
     *    101(10f7.3/:))                               !^
        endif  !----------------------------------------+

        if (iDebug .ne. 0)
     *  Open    (unit=33,file=OutFile(1:lout)//'.deb',
     *                                  status='unknown')
c #######################################################
c              Print surface layer profile:
c #######################################################
        write   (3,61)
  61    format  (/'The Structure Profile:'/
     *           '====================='/
     *  '   Nr    Position Thickness Roughness Trans.layer ',
     *  '________Code________ Fraction    rho      Exp(-w)',
     *  '   x0=(xr0,xi0)   Identical_to_Nr')

        do      l=1,N_Total   !=========================+
          write (3,64)  l, Positions(l), Thickness(l), !|
     *                  Sigma(l), Thickness_Tr(l),     !|
     *                  Code(1,l), Frac(1,l), rho(l),  !|
     *                  W0(l), x0(l+1), Identical(l)   !|
          do j=2,N_Frac(l)  !===================+       |
            write (3,65) Code(j,l), Frac(j,l)  !|       |
          enddo  !==============================+       |
        enddo  !========================================+
  64    format  (i5,f13.2,f10.2,f8.2,f12.2,3x,a,f7.3,f10.4,f11.3,
     *                                       ' (',g8.2,',',g8.2,')',i6)
  65    format  (37x,                   a,f7.3)

        do      i=1,N_Total-1   !==========+
          DDD = Thickness(i)              !|
     -        - Thickness_Tr(i)           !|
     -        - Thickness_Tr(i+1)         !|
                                          !|
          if (DDD .lt. 0.   .AND.         !|
     *        DDD .gt. (-1.E-3)) DDD = 0. !|
                                          !|
          if (DDD .lt. 0.) goto 103 !------+---+
          Thickness(i) = DDD              !|   |
        enddo  !===========================+   v

        if (i_standing .ne. 0) Then !--------------------------------+
          if (n_standing .gt. 1 .and. nf0 .gt. 1) Then !----------+  |
            Call OpenList (OutFile(1:lout)//'_sw.grd',2,ifail)   !|  |
            if (ifail .ne. 0) goto 28                            !|  |
            write (2,56) n_standing,        nf0,                 !|  |
     *                   standing_range(1), standing_range(2),   !|  |
     *                   f0_range(1),       f0_range(2),         !|  |
     *                       0.,                0.               !|  |
  56        Format ('DSAA'/i8,1x,i8,3(/2g15.7))                  !|  |
          else !--------------------------------------------------|  |
            Call OpenList (OutFile(1:lout)//'_sw.dat',2,ifail)   !|  |
            if (ifail .ne. 0) goto 28                            !|  |
          endif !-------------------------------------------------+  |
        endif !------------------------------------------------------+

c #######################################################
c                  Start of processing:
c #######################################################
        Start  = 0.
        Call Duration (start,bias,hh,mm,ss)
        hh = 0
        mm = 0
        ss = 0
        Call Convert_Time (hh,mm,ss,timestamp)
        nerror = 0
        nf = 0
        f0 = f0_range(1)
        progress = 0.
        if (iBatch_mode .lt. 2)
     *    write (*,88)  nf, progress,
     *                  f0, Uni(iUnits+1),
     *                  timestamp
  88    format (1x,i5,' points done (',f4.0,'%). Incidence=',
     *  g11.4,a,' Elapsed time=',a)
        if (modebat .ne. 0)     Then  !---------------+
          DoneFile = OutFile(1:lout)//'.~~~'         !|
          Call OpenFile(DoneFile(1:lout+4),50,       !|
     *             'write','unknown',io_status,*89)  !|
          write (50,90,err=91)  nf, nf0,             !|
     *                          f0,Uni(iUnits+1),    !|
     *                          progress, timestamp  !|
  91      continue                                   !|
          close(unit=50)                             !|
        endif  !--------------------------------------+
  89    continue
  90    format ('Points done = ',i5,' of ',i5,' <br> ',
     *  'Current angle = ',g12.5,1x,a,' <br> ',
     *  'Progress = ',f10.2,'% <br> ',
     *  'Elapsed time (hh:mm:ss) =',a,' <br> ')
c #######################################################
c                   Loop over f0
c #######################################################

c ccc   Call    DivZ()
c YES message from DivZ!
        do      nf=1,nf0  !==============================+
c NO message from DivZ!                                  |
          f0     = f0_range(1) + df0*(nf-1)             !|
          f0_    = sin (f0 * UnitCoef) / UnitCoef       !|
                                                        !|
c ccccc   iPol = 1  !currently sigma-polarization only  !|
          R0s  = SFMm5(f0_,iPol)                        !|
                                                        !|
          if (ifail .ne. 0) Then !------------+          |
            write   (3,77)       ifail,nf    !v          v
  77        format  ('SFMm5: Failure code ',
     *               i4,' -- at step: ',i5)  !^          ^
          endif  !----------------------------+          |
          if (iDebug .ne. 0)  write (33,*) f0,N_Used    !|
                                                        !|
          if (R0s.lt.1.D-99) R0s=0.                     !|
          write (1,'(2g15.6)',err=161) f0, R0s          !|
                                                        !|
          if (i_standing .ne. 0) Then !---------------+  |
            if (n_standing .gt. 1) Then !-----------+ |  |
              if (nf0 .gt. 1) Then !--------------+ | |  |
                write (2,'(1001e15.6)',err=161)  !| | |  |
     *                  (SW(i),i=1,n_standing)   !| | |  |
              else !------------------------------| | |  |
                do i=1,n_standing  !============+ | | |  |
                  z = sngl(                    !| | | |  |
     (                Positions(i_reference+1) !| | | |  |
     +              + standing_range(1)        !| | | |  |
     +              + (i-1) * standing_step    !| | | |  |
     )                )                        !| | | |  |
                  write (2,'(2e15.6)',err=161) !| | | |  |
     *                  z, SW(i)               !| | | |  |
                enddo !=========================+ | | |  |
              endif !-----------------------------+ | |  |
            else !----------------------------------| |  |
              write (2,'(2e15.6)',err=161)         !| |  |
     *                  f0, SW(1)                  !| |  |
            endif !---------------------------------+ |  |
          endif !-------------------------------------+  |
                                                        !|
          if (nf .eq. (nf/iBeepStep)*iBeepStep) Then !-+ |
            progress = 100.*Real(nf)/Real(nf0)        !| |
            if (progress .gt. 99.99) progress = 99.99 !| |
            Call Duration (start,bias,hh,mm,ss)       !| |
            Call Convert_Time (hh,mm,ss,timestamp)    !| |
            if (iBatch_mode .lt. 2)                   !| |
     *        write (*,88)  nf, progress,             !| |
     *                      f0, Uni(iUnits+1),        !| |
     *                      timestamp                 !| |
            if (modebat .ne. 0) Then  !--------------+ | |
              Call OpenFile(DoneFile(1:lout+4),50,  !| | |
     *             'write','unknown',io_status,*93) !| | |
              write(50,90,err=92)nf, nf0,           !| | |
     *                           f0, Uni(iUnits+1), !| | |
     *                           progress, timestamp!| | |
  92          continue                              !| | |
              close(unit=50)                        !| | |
            endif  !---------------------------------+ | |
  93        continue                                  !| |
          endif  !-------------------------------------+ |
                                                        !|
          Call  KeyIn   (0,iscan,iasci)                 !|
c If interrupted by <Esc> and confirmed:                !|
          if  (iasci .eq. 27) Then !-------------------+ |
            write (3,444) nf, nf0                     !| |
            goto 301  !->------------------------------+-+-+
                                                      !| | |
          elseif (FileExist(OutFile(1:lout)//'.esc')) !| | |
     *                                Then !-----------| | v
            Call DelFile (OutFile(1:lout)//'.esc',iii)!| | |
            write (3,444) nf,nf0                      !| | |
            goto 301  !->------------------------------+-+-|
          endif  !-------------------------------------+ | |
                                                        !| |
        enddo  !=========================================+ v
  444   format(/'User interrupted at point',i6,' of',i6)  !|
                                                          !|
  301   continue    !<-------------------------------------+
        if (nf .gt. nf0) nf = nf0
c #######################################################
        close   (unit=1)
        Call    Duration (start,bias,hh,mm,ss)
        Call    Convert_Time (hh,mm,ss,timestamp)
        write   (3,63)  timestamp
  63    format  (/'   Elapsed time (hh:mm:ss)=',a)
        close   (unit=3)
        if (modebat .ne. 0) Then !--------------------+
          Call OpenFile(DoneFile(1:lout+4),50,       !|
     *             'write','unknown',io_status,*95)  !|
          f0       = f0_range(1) + df0*(nf-1)        !|
          progress = 100.*Real(nf)/Real(nf0)         !|
          write (50,90,err=94) nf, nf0,              !|
     *                         f0,Uni(iUnits+1),     !|
     *                         progress, timestamp   !|
  94      continue                                   !|
          close(unit=50)                             !|
        endif  !--------------------------------------+
  95    continue
        if (nf .lt. nf0)      goto 28

c cc    txt(1)=' *** Work done. Press any key to continue..'
        txt(1)=' *** Work done...'
        if (iBatch_mode .lt. 2) Then !----+
          l = Len_Trim (txt(1))          !|
          write (*,'(1x,a)') txt(1)(1:l) !|
        endif  !--------------------------+
        do      i=1,1   !==================+
          Call  KeyIn   (0,iscan,iasci)   !|
          if (iscan .ne. 0)   goto 28     !|
        enddo  !===========================+
c       if (iBatch_Mode .eq. 0)
c    *  Call  KeyIn   (1,iscan,iasci)

  28    continue
        Call exit_quiet()

c #######################################################
c                        Errors
c#######################################################
  99    continue
        write   (txt,117)  ProgramVer(1:ipv), InpFile(1:linp)
  117   format  (a,': cannot open input file:'//a//
     *  'File not found in current directory.')
        ifail = 5
        goto 130

  100   continue
        write   (txt,118)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line
c 118   format  (a,': Incorrect input parameters in file:'//a//
  118   format  (a,': Incorrect syntax in file:'//a//
     *  'at line:',i4,' -- while reading:'/)
        ifail = 7                                               !+txt(7)
        goto 120

  101   continue
        write   (txt,119)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line
  119   format  (a,': Parameter(s) out of range in file:'//a//
     *  'at line:',i4,' -- while reading:'/)
        ifail = 7                                               !+txt(7)
        goto 120

  104   continue
        write   (txt,114)  ProgramVer(1:ipv), i
  114   format  (a,':'/
     *  'Transition sublayers are requested for layer [',i5,' ].'/
     *  'This cannot be combined with request for standing waves.')
        ifail = 3
        goto 120

  105   continue
        write   (txt,115)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     i_reference, N_Total-1
  115   format  (a,': contradictory input in:'//
     *  a//
     *  'Requested reference interface [',i6,'] for standing waves'/
     *  'exceeds the maximum interface index in the structure [',i6,']')
        ifail = 6
        goto 120

  123   continue
        write   (txt,124)  ProgramVer(1:ipv), InpFile(1:linp)
  124   format  (a,': Zero substrate reflectivity in file:'//a//
     *  'Either a substrate code must be entered or x0#0.'/
     *  'Also make sure that W0>0.'/)
        ifail = 7
        goto 120

  131   continue
        write   (txt,132)  ProgramVer(1:ipv)
  132   format  (a,':'//
     *  'X-ray wavelength was chosen to be specified via'/
     *  'characteristic X-ray line, but no line given.')
        ifail = 4
        goto 130

  109   continue
        write   (txt,209)  ProgramVer(1:ipv), InpFile(1:linp)
  209   format  (a,': Inconsistent parameters in file:'//a//
     *          'Multiple scan points with zero scan range.'/
     *          '-- while reading:')
        ifail = 7                                               !+txt(7)
        goto 120

  120   continue
        if (ifail .lt. 0) ifail = 0
        if (ifail .gt. 18) ifail=18
        txt(ifail+1) = ' '
        txt(ifail+2) = 'Check your input file!'
        ifail  = ifail+2
        goto 130

  701   continue
        qz_max = SQZ/Sigma(l)
c This should be changed after we add scanning over qz:
c (See mag_sl99 for an example)
        f0_max = sngl(qz_max/(2.*k0))
        if (f0_max .lt. 1.) Then  !-------+
          f0_max = Asin(f0_max)/UnitCoef !|
        else !----------------------------|
          f0_max = (pi/2.)/UnitCoef      !|
        endif  !--------------------------+

        write   (txt,702)  ProgramVer(1:ipv),
     *                     Sigma(l), l,
     *                     Sigma_max,
     *                     f0_max, Uni(iUnits+1)
  702   format  (a,':'//
     *  'Requested rms roughness sigma =',g9.3,' Angstrom for layer =',
     *                                               i5,' is too high.'/
     *  'The maximum sigma with the Nevot-Croce model and given',1x,
     *                                          'incidence angles is:'//
     *  'sigma_max =',g9.3,' Angstrom'//
     *  'Please use transition layers to model roughness with',1x,
     *                                                     'high sigma'/
     *  'or restrict the incidence angles by '//
     *  'f0_max =',g9.3,1x,a//
     *  'Note that transitions layers are incompatible with wavefields',
     *                                                    1x,'request.')
        ifail = 13
        goto 130

  103   continue
        write   (txt,203)  ProgramVer(1:ipv),i,
     *                     Thickness(i),
     -                     Thickness_Tr(i),
     -                     Thickness_Tr(i+1),
     *                     DDD
  203   format  (a,':'//
     *  'Negative thickness of layer [',i5,' ] due to the'/
     *  'subtraction of transition domain:'//
     *  't - tr_upper - tr_lower = ',f12.2,2(' - ',f9.2),' = ',f10.2/)
        ifail = 7
        goto 130

  161   continue
        write   (txt,121)  ProgramVer(1:ipv),
     *                     OutFile(1:lout)//'.dat'
  121   format  (a,': Cannot write output file:'//a//
     *  'May be, device full or drive not ready.')
        ifail = 5
        close   (unit=1,status='delete',err=130)
        if (i_standing .ne. 0)
     *  close   (unit=2,status='delete',err=130)
        goto 130

  130   continue
        if (modebat .eq. 1) Then  !--------------+
          do    i=1,Iabs(ifail)  !=============+ |
            l = Max0(Len_Trim(txt(i)),1)      !| |
            write (*,'(1x,a)') txt(i)(1:l)    !| |
          enddo  !=============================+ |
        endif  !---------------------------------+
        Call    Message (txt,Iabs(ifail),2)
        goto 28

        End
