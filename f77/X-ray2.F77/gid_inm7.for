c       Reading input file for Gid_sl, ver.M7
c
c  Contents:
c       Subroutine      Read_gid_InputM7
c Calls:
c    1. algebra.for + algebra3.for
c          -- Subroutines BackLat,VecSca,VecMod,Angle,...
c    2. build_k0.for
c          -- Subroutine Build_k0
c    3. get_xh97.for
c          -- Subroutine Get_Xh97
c    4. bragan.for + chiparse.for + chipars2.for + getwave.for +
c          + IndCheck.for + SigmaPi.for + X0h1.for +
c          + X0hMake.for + X0hRead.for + X0hf1f2.for
c    5. mkerrfile.for
c          -- Subroutine Make_Err_Filename
c    6. rdm7prfh.for
c          -- Subroutine Read_TopM7_0h
c          -- Subroutine Read_SubM7_0h
c    7. surfnorm.for
c          -- Subroutine SurfNormVect
c          -- Subroutine SurfNormCoplanar
c          -- Subroutine FindDisorient
c    8. xinput97.for
c          -- Subroutine LineLoop,
c          -- Subroutine Compress_Line,
c          -- Subroutine ArgSearch,
c          -- Subroutine ChiSearch
c --------------------------------------------------------
        Subroutine Read_GID_InputM7 (indices,
     *                               DisorPlanes8, Pfi, igie,
     *                               QB, fcentre, iPol,
     *                               Code, Frac, N_Frac,
     *                               N_Frac_Max, C0h, Aauto)
c--------------------------------------------------------
c+===============================+
        Include 'gid_slm7.inc'  !|
c+===============================+
c The above INC file contains smth like:
c ccc   Parameter       (N_Top_Max     = 1001,
c ccc*                   N_Total_Max   = N_Top_Max+1)

        Integer         N_Frac_Max

        Complex*8       x0_,  xh_,
     *                  x0_f, xh_f

        Real*8          DisorPlanes8

        Real*4          Frac(N_Frac_Max,N_Total_Max),
     *                  UnitData(6), QB, rho, TER_Angle,
     *                  a_substrate, p_substrate, Pfi,
     *                  xr0_Sub, xi0_Sub,
     *                  xrh_Sub, xih_Sub,
     *                  xr0, xi0, xrh, xih,
     *                  x0_abs, xh_abs_sg, xh_abs_pi,
     *                  xh_abs_sg_max, xh_abs_pi_max,
     *                  fcentre, QB_, cos2QB,
     *                  xf_, a_, p_, xf_f,
     *                  AngMinute, AngSec,
     *                  AngMrad, AngMurad, eV2keV,
     *                  rpar(1), xxx,
     *                  d_spacing, DisorPlanExt4,
     *                  asymmetry_mono,         !M7: should go into /RGI_M7/
     *                  Scan_Max /45./,         !Max (+-) scan angles
     *                  Scan_eV_Max /1000./     !Max (+-) scan energy
        Integer         indices(3), ipar(1), iPol,
     *                  N_Frac(N_Total_Max), igie,
     *                  iSyngony, iSyngony_substrate,
     *                  nonCubic, nonCubic_substrate,
     *                  nreflect, line, linp, lout, ltop, ixway,
     *                  mode_x0h, ifail, i, l, j, m,
     *                  mm1, mm2, luninp /2/,
     *                  mono_convolution,       !M7: flag (0=no,1=yes)
     *                  nreflections_mono,      !M7: range=[1-5]
     *                  iu2warn, lurl, io_status

        Character       Code(N_Frac_Max,N_Total_Max)*20,
     *                  C0h(N_Total_Max)*2,
     *                  Aauto(N_Total_Max+1)*1,
     *                  Radiat*6, wrk*80, url*64,
     *                  Symmetry(9)*24 /
     *                   'Unknown',                                     !0
     *                   'Cubic',                                       !1
     *                   'Hexagonal/Trigonal',                          !2
     *                   'Tetragonal',                                  !3
     *                   'Trigonal(Rhombohedral)',                      !4
     *                   'Orthorhombic',                                !5
     *                   'Monoclinic',                                  !6
     *                   'Triclinic',                                   !7
     *                   '** Amorphous ** '/                            !8

        Real*8          DVecMod, Dangle2
        External        DVecMod, Dangle2

        Real            wave2energy
        External        wave2energy

        Character       Get_Server_URL*64
        External        Get_Server_URL

c----------------------
        Complex*8  x0(N_Total_Max+1),
     *             xh(N_Total_Max+1)
        Real*4     xf(N_Total_Max+1),           !phase shift between xrh & xih (non-cubic)
     *             xs(N_Total_Max+1),           !M7: xh phase shift in layer (Kaganer, monolayers)
     *             Thickness(N_Total_Max),      !thicknesses of layers
     *             Sigma(N_Total_Max),          !rms rougness of upper interface
     *             Daa(N_Total_Max+1),          !normal lattice strain
     *             Wave,                        !X-wavelength
     *             xabs,                        !abs(x0_max)
     *             TERang,                      !critical TER angle
     *             standing_range(2),           !M7: offsets in A from reference interface
     *             standing_step,               !M7: SW depth step in A
     *             standing_phase,              !M7: SW phase in units of pi (-1 = no fixed phase)
     *             mono_asymmetry,              !M7: g0/gh for monochromator
     *             xh_weak_percent,             !M5: the weakest possible xh
     *             xh_weak_limit,               !M5: the weakest possible xh
     *             alpha_max                    !M5: max(alpha/xh) to treat crystal as "crystal"
        Integer    N_Top,                       !top layer sublayers
     *             N_Total,                     !total of layers
     *             i_standing,                  !M7: 1=yes/0=no SW
     *             i_reference,                 !M7: reference interface for SW (0=surface)
     *             n_standing,                  !M7: Number of SW depth points
     *             m_reduction,                 !M7: matrix size reduction flag (0/1/2)
     *             m_type,                      !M7: matrix type 4x4,3x3a,3x3b,2x2a,2x2b.
     *             iDebug                       !Debug flag
c This is common between all 3 programs:
c gid_sl, gid_in and RFM:
        Common /GEN_M7/ x0, xh, xf, xs,
     *                  Thickness, Wave,
     *                  Sigma, Daa,
     *                  xabs, TERang,
     *                  xh_weak_limit,          !M5
     *                  alpha_max,              !M5
     *                  standing_range,         !M7: new
     *                  standing_step,          !M7: new
     *                  standing_phase,         !M7: new
     *                  i_standing,             !M7: new
     *                  i_reference,            !M7: new
     *                  n_standing,             !M7: new
     *                  m_reduction,            !M7: new
     *                  m_type,                 !M7: new
     *                  N_Top, N_Total,
     *                  iDebug

c----------------------
        Real*8          DisorNorm, Surface_Normal(3),
     *                  h_vector(3)
        Real*4          W0(N_Total_Max), Wh(N_Total_Max),
     *                  Scan_limits(2), Thickness_Top,
     *                  dScan_angle, AngDegree, pi,
     *                  UnitCoef(3), Asymmetry_par
        Integer         ipv, Mode_Scan, nScan, iBatch_mode,
     *                  nBaseNormal(3), nReperDirection(3),
     *                  iOther_Scan_Vector(3), Invert_Axis,
     *                  iUnits(3), iBeepStep
        Character       ProgramVer*8,
     *                  InpFile*80,
     *                  OutFile*80,
     *                  File_Top*80,
     *                  place*12, Comment(3)*80
c This is common between
c gid_sl and gid_in only:
        Common /RGI_M7/ ProgramVer, ipv,
     *                  InpFile, OutFile,
     *                  W0, Wh, Thickness_Top,
     *                  nBaseNormal,
     *                  nReperDirection,
     *                  DisorNorm,
     *                  Surface_Normal,
     *                  h_vector, Mode_Scan,
     *                  Scan_limits, nScan,
     *                  iBatch_mode, dScan_angle,
     *                  AngDegree, pi, iUnits,
     *                  UnitCoef, iBeepStep,
     *                  File_Top, place,
     *                  iOther_Scan_Vector,
     *                  Invert_Axis,
     *                  Asymmetry_par,
     *                  Comment

c----------------------
        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        Real            Poisson, edge_nearest
        Character       name_nearest*4
        Common  /x0pa9/ Poisson, edge_nearest, name_nearest

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv .lt. 10) Then  !------+
           iirezv = 1                       !|
           istackrezv = istackrezv+1        !|
           progrezv(istackrezv) = progname  !|
           progname = 'Read_GID_InputM7'    !|
        else !-------------------------------+
           iirezv = 0                       !|
        endif  !-----------------------------+
c -------------------------------------------------------
        ifail     = 0
        pi        = 4.*atan(1.)
        AngDegree = 2.*pi/360.
        AngMinute = AngDegree/60.
        AngSec    = AngDegree/3600.
        AngMrad   = 1.E-03
        AngMurad  = 1.E-06
        eV2keV    = 1.E-03

        UnitData(1) = AngDegree
        UnitData(2) = AngMinute
        UnitData(3) = AngMrad
        UnitData(4) = AngSec
        UnitData(5) = AngMurad
        UnitData(6) = eV2keV

        Daa(1)    = 0.                                          ! in vacuum
        xs(1)     = 0.                                          ! in vacuum; added in M7
        Aauto(1)  = ' '
        do l=1,N_Total_Max  !=========+
           N_Frac(l) = 1             !|
           do j=1,N_Frac_Max  !====+  |
              Frac(j,l) = 0.      !|  |
              Code(j,l) = ' '     !|  |
           enddo  !================+  |
           Frac(1,l)  = 1.           !|
           Daa(l+1)   = 0.           !|
           Aauto(l+1) = ' '          !|
        enddo  !======================+

        i_standing  = 0
        i_reference = 0
        n_standing  = 0

        modebat = 1     !assume batch-mode by default
c#######################################################
c              Read input file
c#######################################################
        linp = Max (Len_Trim(InpFile),1)
        Call OpenFile(InpFile(1:linp),luninp,'read','old',
     *                                     io_status,*99)

        Call Make_Err_Filename (InpFile(1:linp))

        line = 0
c-----------------------
        txt(7) = 'Program version: '//ProgramVer(1:ipv)
        Call LineLoop (luninp,line,wrk,*90)
        Call CaseLat (wrk,1)                    !convert to lower case
        if (wrk(1:ipv) .ne. ProgramVer(1:ipv)) goto 123
c-----------------------
        txt(7) = 'batch mode flag [0/1/2]'
        Call LineLoop (luninp,line,wrk,*90)
        Call rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        iBatch_mode = ipar(1)
        if (iBatch_mode .lt. 0)   goto 101
        if (iBatch_mode .gt. 2 .AND.
     *      iBatch_mode .ne. 10)  goto 101
c (iBatch_mode=10 is used for debugging -- it raises iDebug=1)
        if (iBatch_mode .eq. 0) iBatch_mode=1
        modebat = iBatch_mode
c-----------------------
        txt(7) = 'output file name'
        Call LineLoop (luninp,line,wrk,*90)
        OutFile = wrk
        lout = Len_Trim(OutFile)
        if (lout .eq. 0) Then !-----+
           OutFile = '~'           !|
           lout    = 1             !|
        endif  !--------------------+
c       Call Case    (OutFile(1:lout),1)                        !to lower case
        Call OpenList (OutFile(1:lout)//'.tbl',3,ifail)
        if (ifail .ne. 0) goto 28
c-----------------------
        txt(7) = 'surface profile name'
        N_Top   = 0
        Call LineLoop (luninp,line,wrk,*90)
        File_Top = wrk
        ltop = Max(Len_Trim(File_Top),1)
c       Call Case (File_Top(1:ltop),1)                          !to lower case
c-----------------------
        txt(7) = 'surface layer profile'
        Call Read_TopM7_0h (File_Top(1:ltop),                   !this is in rdm7prfh.for
     *                      N_Top, N_Top_Max,
     *                      Thickness, Thickness_Top,
     *                      Code(1,1), Frac(1,1),
     *                      N_Frac(1), N_Frac_Max,
     *                      x0(2), xh(2), xf(2),
     *                      W0(1), Wh(1),
     *                      Daa(2), Sigma(1),
     *                      xs(2),                              !added im M7
     *                      C0h(1), ifail)
        if (ifail .ne. 0) goto 28
        N_Total = N_Top+1
c-----------------------
        do i=1,3  !======================================+
           write (txt(7),'(a,i1)') 'comment line #',i   !|
           Call LineLoop (luninp,line,Comment(i),*90)   !|
        enddo  !=========================================+
c-----------------------
        txt(7) = 'substrate data'
        Call Read_SubM7_0h (Code(1,N_Total),
c cc *  >>>>>--------->>>>> x0(N_Total+1), xh(N_Total+1), !<--+
     *                      W0(N_Total), Wh(N_Total),        !|
     *                      Daa(N_Total+1), Sigma(N_Total),  !|
     *                      xs(N_Total+1),                   !|added im M7
     *                      wrk, line, *90, *101)            !|
c The substrate can only be specified via its code, because  !|
c this is the only way to get the lattice parameters!!!      !|
        if (Len_Trim(Code(1,N_Total)) .eq. 0) Then !--+   !<--+
           goto 765                                  !|
        endif  !--------------------------------------+
c-----------------------
        txt(7)='x-ray specification mode (1=wave 2=energy 3=line 4=QB)'
        Call LineLoop (luninp,line,wrk,*90)
c 1=wave 2=energy 3=x-ray line 4=QB:
        if (Len_Trim(wrk).eq.0) goto 80
        Call rdInt (ipar,1,wrk,i)
        if (i .ne. 0)     goto 100
        ixway = ipar(1)
        if (ixway .lt. 1) goto 101
        if (ixway .gt. 4) goto 101
c-----------------------
        if (ixway.eq.1) txt(7)='x-ray wavelength (range=[0.1-10]A)'
        if (ixway.eq.2) txt(7)='x-ray energy (range=[1.24-124]keV)'
        if (ixway.eq.3) txt(7)='x-ray wavelength (ignored)'
        if (ixway.eq.4) txt(7)='Bragg angle (range=[0-90]degr)'
        Call    LineLoop (luninp,line,wrk,*90)
        Wave = 0.
        QB   = 0.
        if (ixway.lt.3 .OR. ixway.eq.4) Then !-----------+
           if (Len_Trim(wrk).eq.0) goto 80              !|
           Call  rdReal  (rpar,1,wrk,i)                 !|
           if (i.ne.0)             goto 100             !|
           if (ixway.lt.3) Then !----------------------+ |
              Wave = rpar(1)                          !| |
              if (Wave .le. 0.)     goto 101          !| |
              if (ixway.eq.2) Wave=wave2energy(Wave)  !| |
              if (Wave .lt. 0.1)    goto 101          !| |
              if (Wave .gt. 10.)    goto 101          !| |
           elseif (ixway.eq.4) Then !------------------+ |
              QB = rpar(1)                            !| |
              if (QB .le.  0.)      goto 101          !| |
              if (QB .gt. 90.)      goto 101          !| |
           endif  !------------------------------------+ |
        endif  !-----------------------------------------+
c-----------------------
        txt(7) = 'x-ray line name'
        Call LineLoop (luninp,line,wrk,*90)
        if (ixway.eq.3) Then  !----------------------+
           Radiat = wrk(1:6)                        !|
           if (Len_Trim(Radiat).eq.0) goto 131 !-----+---v
c Get wavelength from the X-ray line name:           |
           Call  Get_Wave (Wave,Radiat,ifail)       !|
           if (ifail.ne.0)         goto 130         !|
           if (Wave .lt. 0.1)      goto 101         !|
           if (Wave .gt. 10.)      goto 101         !|
        else !---------------------------------------+
           Radiat = ' '                             !|
        endif  !-------------------------------------+
c-----------------------
        txt(7) = 'x-ray polarization (1=sigma, 2=pi, 0/3=mixed)'
        Call LineLoop (luninp,line,wrk,*90)
        Call rdInt    (ipar,1,wrk,i)
        if (i .ne. 0)       goto 100
        iPol = ipar(1)
c Polarization: 0/3-mixed, 1-sigma, 2-pi:
        if (iPol .eq. 0)  iPol=3
        if (iPol .lt. 1)    goto 101
        if (iPol .gt. 3)    goto 101
c-----------------------
        txt(7) = 'indices of the Bragg reflection'//
     *           ' (three integers, not 0 0 0)'
        Call LineLoop (luninp,line,wrk,*90)
        if (Len_Trim(wrk).eq.0) goto 80
        Call rdInt    (indices,3,wrk,i)
        if (i .ne. 0)             goto 100
        if ((indices(1) .eq. 0) .AND.
     *      (indices(2) .eq. 0) .AND.
     *      (indices(3) .eq. 0))  goto 101
        if (ixway .ne.4) QB = 0.
        x0(1) = (0.,0.)
        xh(1) = (0.,0.)
        xf(1) = 0.
c Get x0,xh for the substrate:
        mode_x0h = 0            ! 1st calculation (new crystal,...)
        iSyngony_substrate = 0  ! any syngony!
        nonCubic_substrate = 0  ! yet we think this is a cubic substrate
        nreflect = 1
        rho      = 0.
        mm1 = 1
        mm2 = N_Total
        Call Get_Xh97  (Wave, Code(mm1,mm2), rho,
     *                  iSyngony_substrate, nreflect,
     *                  indices, QB, mode_x0h,
     *                  xr0_Sub, xi0_Sub,
     *                  xrh_Sub, xih_Sub, xf(N_Total+1),
     *                  x0(N_Total+1), xh(N_Total+1),
     *                  ifail)
c The structure was not found in the database;
c then it was interpreted as a chemical formula,
c but we do not have the option ot supply density.
c Report this as "no structure in DB":
        if (iSyngony_substrate .eq. 8 .AND.
     *                   ifail .lt. 0 .AND.
     *                Abs(rho) .lt. (1.E-20)) goto 302  !-no-record-in-db-v
        if (iSyngony_substrate .eq. 8)        goto 304  !-not-a-crystal---v
        if (ifail.ne.0) Then !-------------------------------------+
           i = Index(txt(3),'does not exist for this wavelength') !|
           if (QB.le.0. .AND. i.gt.0) Then !---------+             |
              goto 306 !-no-reflection-------v       |             |
           else !------------------------------------+             |
              goto 130 !-some-other-failure--v       |             |
           endif  !----------------------------------+             |
        endif  !---------------------------------------------------+
c Forbidden reflection limit:
        xh_weak_percent = 0.5
c       xh_weak_limit = Amax1( (1.E-8)*sin(QB*AngDegree),       !the order of gamma_h
c    *                  (1.E-2)*xh_weak_percent*abs(xr0_Sub) )  !1.0% reflection - change of 2003/06/30
        xh_weak_limit = (1.E-2)*xh_weak_percent*abs(xr0_Sub)    !1.0% reflection - change of 2020/05/20
c Critical angle for total external reflection:
        TER_Angle = Sqrt (Abs(xr0_Sub)) / AngDegree
c Polarization factor for pi:
        cos2QB = abs( cos(2.*QB*AngDegree) )
        a_substrate = a(1)
        p_substrate = Poisson
        if (iSyngony_substrate .gt. 1) nonCubic_substrate = 1
c-----------------------
        txt(7) = 'units for miscut, incidence, and scan angles [0--4]'
        Call LineLoop (luninp,line,wrk,*90)
        if (Len_Trim(wrk).eq.0) goto 80
        Call rdInt    (iUnits,3,wrk,i)
        if (i.ne.0)             goto 100

c Units for:
c iUnits(1) -- surface miscut angle  (unim)
c iUnits(2) -- incidence/exit angles (unic) - can also be "-1" when does not matter
c iUnits(3) -- scan angles           (unis)
        txt(7) = 'units for miscut, incidence, and scan angles [0--4]'
        iu2warn = 0
        if (iUnits(2) .eq. -1) Then !--+
           iUnits(2) = 0              !| 0=degr (see tests if it is allowed below)
           iu2warn = 1                !| raise flag
        endif !------------------------+
        if (iUnits(3) .eq. -1) Then !--+
           iUnits(3) = 5              !| 5=eV  (see tests if it is allowed below)
        endif !------------------------+
        do i=1,3  !==============================+
           if (iUnits(i).lt.0)   goto 101       !|0=degr
           if (i .ne. 3) Then !--------------+   |1=min
              if (iUnits(i).gt.4) goto 101  !|   |2=mrad
           else !----------------------------+   |3=sec
              if (iUnits(i).gt.5) goto 101  !|   |4=urad
           endif !---------------------------+   |5=eV
           UnitCoef(i) = UnitData(iUnits(i)+1)  !|
        enddo  !=================================+

c-----------------------
        txt(7) = 'the type of diffraction geometry specification [1--9]'
        Call LineLoop (luninp,line,wrk,*90)
c Geometry specification:
c 1 - non-coplanar case via the incidence angle of k0
c 2 - non-coplanar case via the exit angle of kh
c 3 - coplanar case; grazing incidence
c 4 - coplanar case; grazing exit
c 5 - symmetric Bragg case: coplanar/non-coplanar (not GID)
c     ------ The following does not use surface normal ------
c 6 - coplanar case specified via the angle of Bragg planes to the surface
c 7 - coplanar case specified via the incidence angle of k0
c 8 - coplanar case specified via the exit angle of kh
c 9 - coplanar case specified via beta=g0/|gh|
        if (Len_Trim(wrk).eq.0) goto 80
        Call rdInt (ipar,1,wrk,i)
        if (i .ne. 0)    goto 100
        igie = ipar(1)
        if (igie .lt. 1) goto 101
        if (igie .gt. 9) goto 101

c----------------------- *** for igie <= 5 only ***
c Read the 3 parameters required to determine the actual INTERNAL
c surface normal "Surface_Normal":
c 1. nBaseNormal - base surface normal vector indices:
c 2. nReperDirection - indices of vector, pointing out
c    the direction of maximum misorientation of real
c    surface normal with respect to base normal vector:
c 3. DisorNorm - Maximum misorientation angle of real
c    surface normal:
c -- then, calculate the actual surface normal.
c NOTE: for igie>5 these are returned all zero (no normal required)
        Call InputSufNorm (luninp, line, igie, wrk,
     *                     nBaseNormal, nReperDirection,
     *                     DisorNorm, UnitCoef(1)/AngDegree,
     *                     Surface_Normal, *90, *101, *28)

c----------------------- *** for geometries [1-2,6-9] only ***
        if     (igie .eq. 1) Then !--------------------------------------+
           txt(7) = 'incidence angle of K0 at the Bragg position'       !|
        elseif (igie .eq. 2) Then !--------------------------------------+
           txt(7) = 'exit angle of Kh at the Bragg position'            !|
        elseif (igie .eq. 3) Then !--------------------------------------+
           txt(7) = 'dummy asymmetry parameter (not used in this mode)' !|
        elseif (igie .eq. 4) Then !--------------------------------------+
           txt(7) = 'dummy asymmetry parameter (not used in this mode)' !|
        elseif (igie .eq. 5) Then !--------------------------------------+
           txt(7) = 'dummy asymmetry parameter (not used in this mode)' !|
        elseif (igie .eq. 6) Then !--------------------------------------+
           txt(7) = 'angle of the Bragg planes to the surface'          !|
        elseif (igie .eq. 7) Then !--------------------------------------+
           txt(7) = 'incidence angle of K0 at the Bragg position'       !|
        elseif (igie .eq. 8) Then !--------------------------------------+
           txt(7) = 'exit angle of Kh at the Bragg position'            !|
        elseif (igie .eq. 9) Then !--------------------------------------+
           txt(7) = 'asymmetry of the Bragg reflection beta=g0/|gh|'    !|
        endif !----------------------------------------------------------+
        Call LineLoop (luninp,line,wrk,*90)
        Call rdReal   (rpar,1,wrk,i)
        if (i .ne. 0) goto 100
        Asymmetry_par = rpar(1)
        if (igie.ne.3 .AND.             !dummy
     *      igie.ne.4 .AND.             !dummy
     *      igie.ne.5 .AND.             !dummy
     *      igie.ne.9 .AND.             !no units
     *      Abs(Asymmetry_par).gt.(1.E-20) .AND.
     *      iu2warn .ne. 0) Then !------------------------------------+
           txt(7) = 'Units for geometry parameter are not specified' !|
           goto 100                                                  !|
        endif !-------------------------------------------------------+
        if     (igie .eq. 1 .OR. igie .eq. 2) Then !--------------------+
c 1 - non-coplanar case via the incidence angle of k0                   |
c 2 - non-coplanar case via the exit angle of kh                        |
                                                                       !|
c Read central point angle (incidence or exit) for igie=1,2:            |
           fcentre = Asymmetry_par                                     !|
c This principally might happen                                         |
c (therefore, the following                                             |
c checking is commented out):                                           |
c cccc     if (fcentre .le. 0.) goto 101                                |
           if (fcentre*UnitCoef(2)/AngDegree .gt. 90.)       goto 101  !|
           if (fcentre*UnitCoef(2)/AngDegree .lt. -10.*                !|
     *                                            TER_Angle) goto 101  !|
                                                                       !|
        elseif (igie .ge. 3 .AND. igie .le. 5) Then !-------------------+
c 3 - coplanar case; grazing incidence                                  |
c 4 - coplanar case; grazing exit                                       |
c 5 - symmetric Bragg case: coplanar/non-coplanar (not GID)             |
                                                                       !|
           Asymmetry_par = 0.                    !not used in this mode |
           fcentre = 0.                          !not used in this mode |
                                                                       !|
        elseif (igie .ge. 6 .AND. igie .le. 9) Then !-------------------+
c 6 - coplanar case specified via the angle of Bragg planes to surface  |
c 7 - coplanar case specified via the incidence angle of k0             |
c 8 - coplanar case specified via the exit angle of kh                  |
c 9 - coplanar case specified via beta=g0/|gh|                          |
                                                                       !|
c Find angle between Bragg planes and the surface                       |
c (coplanar case and no surface indices entered):                       |
c (DisorPlanes8 is calculated in degrees -- it is not the final value!) |
           Call FindDisorient (QB,                         ! in degr.  !|
     *                         TER_Angle,                  ! in degr.  !|
     *                         Asymmetry_par, UnitCoef(2), ! in Unit2  !|
     *                         igie,                       !  6--9     !|
     *                         DisorPlanes8,               ! in degr.  !|(output)
     *                         ifail)                                  !|
           if (ifail .ne. 0) goto 28                                   !|
c Incidence angle at the Bragg position (in Unit1):                     |
           fcentre = sngl(QB+DisorPlanes8)                             !|
c Added June 9, 2002 (see error in GID08377.zip)                        |
           if (fcentre .gt. 90.) fcentre = 180. - fcentre              !|
c Changed April 12, 2005 (see error in GID26454.zip)                    |
           fcentre = fcentre * AngDegree/UnitCoef(2)                   !|
c Build coplanar surface normal:                                        |
           Call SurfNormCoplanar (indices,                 ! hkl       !|
     *                            DisorPlanes8,            ! in degr.  !|
     *                            Surface_Normal,          ! hkl       !|(output)
     *                            ifail)                               !|
           if (ifail .ne. 0) goto 28                                   !|
                                                                       !|
        else !----------------------------------------------------------+
                                                                       !|
           stop 'gid_inm7: Unexpected igie'                            !|
                                                                       !|
        endif !---------------------------------------------------------+

c-----------------------
c Calculate the misorientation of reciprocal lattice vector
c with respect to the internal surface normal:
        Call DVecCopy_id(indices, h_vector, 3)
        DisorPlanes8 = Dangle2 (Surface_Normal,h_vector,3)      !in degr.

c If h_vector points too steep inward - invert it:
c +--------------------------------------+
c |ATTENTION!!! The inverting procedure  | 21.09.95
c |             inhibits GID reflections |
c |             with large inward        |
c |             misorientation !!!       |
c +--------------------------------------+
        if (DisorPlanes8 .lt. (90.D0-5.0001*TER_Angle)) Then !--+
                                                               !v
           write (txt,141)  ProgramVer(1:ipv)
  141      format (a/'W A R N I N G'/
     *     'At the geometry you specified'/
     *     'the Bragg vector points too steep'/
     *     'inside the crystal, and the reflected'/
     *     'diffracted intensity is negligible.'/
     *     'Inverting the Bragg vector for you...')            !^
           if (modebat .ne. 0) Then !---------------+           |
c Write to .tbl:                                    |           |
              do i=1,7  !========================+  |           |
                 j = Len_Trim(txt(i))           !|  |           |
                 if (j .lt. 1) j=1              !|  |           |
                 write (3,'(1x,a)') txt(i)(1:j) !|  |           |
              enddo  !===========================+  |           |
           else  !----------------------------------+           |
c cccccc      Call Message (txt,7,1)               !|wait key   |
              Call Message (txt,7,0)               !|not wait   |
           endif  !---------------------------------+           |
                                                               !|
           do i=1,3  !======================+                   |
              indices(i)  = -indices(i)    !|                   |
              h_vector(i) = -h_vector(i)   !|                   |
           enddo  !=========================+                   |
           DisorPlanes8 = 180.D0 - DisorPlanes8   ! in degr.    |
        endif  !------------------------------------------------+
c This is the angle between h_vector and the surface:
c (it is negative if the reciprocal lattice vector
c points outside):
        Pfi = sngl(90.D0-DisorPlanes8)                          !in degr.
        if (Pfi .gt. 5.0001*TER_Angle)
     *          Stop 'Bragg Planes Determination Logic Error'
c Save angle of reciprocal lattice vector
c to external normal in degrees:
        DisorPlanExt4 = sngl(180.D0 - DisorPlanes8)             !in degr.
c Recalculate asymmetry angles into
c the incidence/exit angles units:
        DisorPlanes8  = DisorPlanes8 * AngDegree/UnitCoef(2)    !to Units-2
        Pfi           = Pfi          * AngDegree/UnitCoef(2)    !to Units-2
c-----------------------
c Scan mode for non-coplanar cases:
c (use for [1,2] & non-coplanar [5];
c in other cases [k0*h] is taken)
        txt(7) = 'mode of scan for non-coplanar cases (1-6) '
     *         //'or energy (7-8)'
        Call LineLoop (luninp,line,wrk,*90)
c Scan axis indicator:
c  1. Round N_surface
c  2. Round [k0*N_surface]
c  3. Round h
c  4. Round [k0*h]
c  5. Round other axis
c  6. Takeoff spectrum (PSD)
c  7. Energy (eV) with recalculation of x0,xh
c  8. Energy (eV) fast (no recalculation of x0,xh) -- cannot be used for grazing cases
        if (Len_Trim(wrk).eq.0) goto 80
        Call    rdInt   (ipar,1,wrk,i)
        if (i.ne.0)             goto 100
        Mode_Scan = ipar(1)
        if (Mode_Scan .le. 6) Then !----------------------+
           if (igie.eq.1 .or.                            !|non-coplanar case via the incidence angle of k0
     *         igie.eq.2 .or.                            !|non-coplanar case via the exit angle of kh
     *        (igie.eq.5 .and.                           !|symmetric Bragg case: coplanar/non-coplanar
     *         abs(DisorPlanExt4).gt.0.1)) Then !-+       | > 0.1degr.
              if (Mode_Scan.lt.1) goto 101       !|       |
           else !---------------------------------+       |
              Mode_Scan = 4                      !|       |coplanar case: take [k0*h]
           endif  !-------------------------------+       |
           if (iUnits(3).gt.4) Then !------------------+  |scan units cannot be in eV
              txt(7) = 'units for scan angle [0--4],' !|  |
     *              //' 5 for energy scans only'      !|  |
              goto 101                                !|  |scan units cannot be in eV
           endif !-------------------------------------+  |
        elseif (Mode_Scan.eq.7 .OR.                      !|
     *          Mode_Scan.eq.8) Then !--------------------+
           iUnits(3) = 5                                 !|scan units must be in eV
           UnitCoef(3) = UnitData(iUnits(3)+1)           !|
        else   !------------------------------------------+
           goto 101                                      !|
        endif  !------------------------------------------+
c-----------------------
        txt(7) = 'indices defining scan axis'
        Call LineLoop (luninp,line,wrk,*90)
c Indices of scan axis
c for Mode_Scan=5:
        if (Mode_Scan .eq. 5) Then  !-------------------+
           if (Len_Trim(wrk).eq.0) goto 80             !|
           Call rdInt (iOther_Scan_Vector,3,wrk,i)     !|
           if (i .ne. 0)           goto 100            !!
           if (iOther_Scan_Vector(1) .eq. 0 .AND.      !|
     *         iOther_Scan_Vector(2) .eq. 0 .AND.      !|
     *         iOther_Scan_Vector(3) .eq. 0 ) goto 110 !|
        else !------------------------------------------+
           iOther_Scan_Vector(1) = 0                   !|
           iOther_Scan_Vector(2) = 0                   !|
           iOther_Scan_Vector(3) = 0                   !|
        endif  !----------------------------------------+
c-----------------------
        txt(7) = 'invert scan axis flag (0=no, 1=yes)'
        Call LineLoop (luninp,line,wrk,*90)
c Flag for inverting the scan axis:
        Call rdInt (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        Invert_Axis = ipar(1)
        if (Invert_Axis .lt. 0)   goto 101
        if (Invert_Axis .gt. 1)   goto 101
c-----------------------
        if (Mode_Scan .le. 6) Then !--------------------------+
           write (txt(7),11) -Scan_Max, Scan_Max, 'degr.'    !|
        else !------------------------------------------------+
           write (txt(7),11) -Scan_eV_Max, Scan_eV_Max, 'eV' !|
        endif ! ----------------------------------------------+
  11    format ('scan limits (range=[',f6.0,' :',f6.0,'] ',a,').')
        Call LineLoop (luninp,line,wrk,*90)
        if (Len_Trim(wrk).eq.0) goto 80
c Read scan angle limits:
        Call rdReal (Scan_limits,2,wrk,i)
        if (i .ne. 0)             goto 100
        if (Mode_Scan .le. 6) Then !--------------------------------+
           xxx = (Scan_Max + 0.0001) * AngDegree                   !|
           if (Abs(UnitCoef(3)*Scan_limits(1)) .gt. xxx) goto 111  !|
           if (Abs(UnitCoef(3)*Scan_limits(2)) .gt. xxx) goto 111  !|
        else !------------------------------------------------------+
           if (Abs(Scan_limits(1)) .gt. Scan_eV_Max)     goto 111  !|
           if (Abs(Scan_limits(2)) .gt. Scan_eV_Max)     goto 111  !|
        endif !-----------------------------------------------------+
c-----------------------
        write (txt(7),'(2a,i5,a)') 'number of scan points ',
     *                             '(min=1, max=',N_Pts_Max,')'
        Call LineLoop (luninp,line,wrk,*90)
c Read scan points:
        if (Len_Trim(wrk).eq.0) goto 80
        Call rdInt (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        nScan = ipar(1)
        if (nScan .lt. 1)         goto 101
        if (nScan .gt. N_Pts_Max) goto 101
        if (nScan .gt. 1) Then !----------------------------+
           dScan_angle = (Scan_limits(2) - Scan_limits(1)) !|
     /                 / (nScan-1)                         !|
           if (Abs(dScan_angle) .lt. (1.E-20)) goto 109    !|
        else !----------------------------------------------+
           dScan_angle = 0.0                               !|
        endif !---------------------------------------------+
c-----------------------
        txt(7) = 'output data columns [ADSBIEM2 -- up to 10 chars]'
        Call LineLoop (luninp,line,wrk,*90)
c Read the angles to include in the output file:
c 1. A -- scan Angle
c 2. D -- Diffracted intensity  [always included at place-2]
c 3. S -- Specular (transmitted) intensity
c 4. B -- Bragg deviation (alpha-parameter)
c 5. I -- Incidence_angle
c 6. E -- Exit angle
c 7. M -- iMaginary part of exit angle
c 8  2 -- 2*theta (the sum of incident & exit angles)
        Call txtShift (wrk,j,i)         !i=Len_Trim(wrk)
        i = Len(place)
        place = wrk(1:i)
        i = Len_Trim(place)
        if (i .gt. Len(place)-1)  goto 101
c Default (when nothing specified):
        if (i .eq. 0)         place='ADSIEB'
cccc    if (Mode_Scan .eq. 6) place='A'
        if (Mode_Scan .eq. 7 .OR.
     *      Mode_Scan .eq. 8) place='A'
        Call Case (place,0)             !to upper case
        i = Len_Trim(place)
        do j=1,i  !===============================+
           if (place(j:j) .ne. 'A' .AND.         !|
     *         place(j:j) .ne. 'D' .AND.         !|
     *         place(j:j) .ne. 'S' .AND.         !|
     *         place(j:j) .ne. 'B' .AND.         !|
     *         place(j:j) .ne. 'I' .AND.         !|
     *         place(j:j) .ne. 'E' .AND.         !|
     *         place(j:j) .ne. 'M' .AND.         !|
     *         place(j:j) .ne. '2')  goto 101    !|
        enddo  !==================================+
c-----------------------
        txt(7) = 'Max alpha/x0 when crystal is treated as crystal'//
     *                                              '[1e+1 -- 1e+8]'
        Call LineLoop (luninp,line,wrk,*90)
        Call rdReal  (rpar,1,wrk,i)
        if (i .ne. 0)        goto 100
        alpha_max = rpar(1)
        if (Abs(alpha_max) .lt. (1.E-20)) alpha_max=1.E+8
        if (Abs(alpha_max) .lt. (1.E+1))  goto 101
c This control is enforced in Rfms_Urd.for -> Dispersion_Roots:
c       if (Abs(alpha_max) .gt. (1.E+8))        goto 101
c-----------------------
        txt(7) ='Scattering matrix reduction flag: 0=no,1=prescan,2=yes'
        Call LineLoop (luninp,line,wrk,*90)
        Call rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)        goto 100
        m_reduction = ipar(1)
        if (m_reduction .lt. 0 .OR.
     *      m_reduction .gt. 2)     goto 101
c-----------------------
        txt(7) = 'standing waves printout option (0=no, 1=yes)'
        Call LineLoop (luninp,line,wrk,*90)
c Read standing waves option:
        Call rdInt (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        i_standing = ipar(1)
        if (i_standing .lt. 0)  goto 101
        if (i_standing .gt. 1)  goto 101
c-----------------------
        txt(7) = 'standing waves reference interface (0=surface)'
        Call LineLoop (luninp,line,wrk,*90)
        if (i_standing .gt. 0) Then !------------+
c Read standing waves option:                    |
c          if (Len_Trim(wrk).eq.0) goto 80      !|
           Call rdInt (ipar,1,wrk,i)            !|
           if (i .ne. 0)           goto 100     !|
           i_reference = ipar(1)                !|
        else !-----------------------------------+
           i_reference = 0                      !|
        endif !----------------------------------+
        if (i_reference .lt. 0)    goto 101
c-----------------------
        txt(7) = 'standing waves offsets from given interface [from,to]'
        Call LineLoop (luninp,line,wrk,*90)
        if (i_standing .gt. 0) Then !---------------+
c          if (Len_Trim(wrk).eq.0) goto 80         !|
           Call rdReal (standing_range,2,wrk,i)    !|
           if (i .ne. 0)           goto 100        !|
        else !--------------------------------------+
           standing_range(1) = 0.                  !|
           standing_range(2) = 0.                  !|
        endif !-------------------------------------+
c-----------------------
        write (txt(7),'(a,i3,a)') 'number of standing wave points (1-',
     *                          N_Stand_Max,')'
        Call LineLoop (luninp,line,wrk,*90)
        if (i_standing .gt. 0) Then !--------------+
c Read standing waves option:                      |
c          if (Len_Trim(wrk).eq.0) goto 80        !|
           Call rdInt (ipar,1,wrk,i)              !|
           if (i .ne. 0)           goto 100       !|
           n_standing = ipar(1)                   !|
           if (n_standing .eq. 0) n_standing = 1  !|
        else !-------------------------------------+
           n_standing = 1                         !|
        endif !------------------------------------+
        if (n_standing .lt. 1)           goto 101
        if (n_standing .gt. N_Stand_Max) goto 101
        if (Abs(standing_range(2)-standing_range(1)) .lt. (1.E-20))
     *                                                  n_standing = 1
        if (n_standing .gt. 1) Then !-------------------------------+
           standing_step = (standing_range(2)-standing_range(1))   !|
     /                   / (n_standing-1.)                         !|
        else !------------------------------------------------------+
           standing_step = 0.                                      !|
        endif !-----------------------------------------------------+
c-----------------------
c Standing waves phase:
c 0.=at the Bragg planes
c 1.=between the Bragg planes
        txt(7) = 'standing waves phase/pi [0.--2.]'
        Call LineLoop (luninp,line,wrk,*90)
        Call CaseLat (wrk,1)                    !convert to lower case
        if (len_trim(wrk).eq.0 .OR. wrk.eq.'none') Then !--+
c No fixed XSW phase                                       |
           standing_phase = -1.                           !|
        else !---------------------------------------------+
c Fixed XSW phase                                          |
           Call rdReal (rpar,1,wrk,i)                     !|
           if (i .ne. 0)                goto 100          !|
           standing_phase = rpar(1)                       !|
           if (standing_phase .lt. 0.)  goto 101          !|
           if (standing_phase .gt. 2.)  goto 101          !|
        endif  !-------------------------------------------+
c-----------------------
c MONO: curve convolution (1:Yes/0:No)..| 0
c MONO: asymmetry beta=g0/|gh|..........| 1.
c MONO: # of reflections [1--5].........| 1
c Read mono convolution option:
        txt(7) = 'Mono convolution (0=no, 1=yes)'
        Call LineLoop (luninp,line,wrk,*90)
        Call rdInt (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        mono_convolution = ipar(1)
c Mono convolution option is not implemented yet:
        if (mono_convolution .ne. 0) goto 113 !------------------------+
c-----------------------
c MONO: asymmetry beta=g0/|gh|..........| 1.
        txt(7) = 'Mono asymmetry [0.01 - 100.]'
        Call LineLoop (luninp,line,wrk,*90)
        if (mono_convolution .gt. 0) Then !-------------+
c          if (Len_Trim(wrk).eq.0) goto 80             !|
           Call rdReal (rpar,1,wrk,i)                  !|
           if (i .ne. 0)           goto 100            !|
           mono_asymmetry = rpar(1)                    !|
        else !------------------------------------------+
           mono_asymmetry = 1.                         !|
        endif !-----------------------------------------+
        if (mono_asymmetry .lt. 0.01)   goto 101
        if (mono_asymmetry .gt. 100.)   goto 101
c-----------------------
c MONO: # of reflections [1--5].........| 1
        txt(7) = 'Mono number of reflections [1-5]'
        Call LineLoop (luninp,line,wrk,*90)
        if (mono_convolution .gt. 0) Then !-------------+
c          if (Len_Trim(wrk).eq.0) goto 80             !|
           Call rdInt (ipar,1,wrk,i)                   !|
           if (i .ne. 0)           goto 100            !|
           nreflections_mono = ipar(1)                 !|
        else !------------------------------------------+
           nreflections_mono = 1                       !|
        endif !-----------------------------------------+
        if (nreflections_mono .lt. 1) goto 101
        if (nreflections_mono .gt. 5) goto 101
c-----------------------
c Mono convolution option is not implemented yet:
        asymmetry_mono    = 1.
        nreflections_mono = 1
c-----------------------
        iBeepStep = 10
        iDebug    = 0
        if (iBatch_mode.eq.10) then !---+
           iDebug = 1                  !|
        endif !-------------------------+
        Close   (unit=luninp)
c#######################################################
c Added 2011/03/31: backup substrate lattice for Algrebra routines
c so that we can restore it after calling X0h on the profile stack:
        Call    LatticeBackup ()
        do l=1,N_Top  !==============================================+
           if ((Len_Trim(Code(1,l)).eq.0) .OR.                      !|
     *         (Code(1,l).eq.'Unknown'))                            !|
     *                         Code(1,l) = Code(1,N_Total)          !|
           x0_ = (0., 0.)                                           !|
           xh_ = (0., 0.)                                           !|
           xf_ = 0.                                                 !|
           a_  = 0.                                                 !|
           p_  = 0.                                                 !|
c Controls for automatic da/a calculation:                           |
           if (Abs(Daa(l+1)) .gt. (1.E+30)) Then !-------+           |
c Non-cubic substrate:                                   |           |
              if (nonCubic_substrate .eq. 1) goto 115 !--+-----------+---v
           endif  !--------------------------------------+           |
           do j=1,N_Frac(l)  !=============================-======+  |
              mm1 = j                                            !|  |
              mm2 = l                                            !|  |
              if (Code(j,l) .eq. Code(1,N_Total)) Then !-------+  |  |
                 x0_f     = x0(N_Total+1)                     !|  |  |
                 xh_f     = xh(N_Total+1)                     !|  |  |
                 xf_f     = xf(N_Total+1)                     !|  |  |
                 a(1)     = a_substrate                       !|  |  |
                 Poisson  = p_substrate                       !|  |  |
                 iSyngony = iSyngony_substrate                !|  |  |
                 nonCubic = nonCubic_substrate                !|  |  |
              else  !------------------------------------------+  |  |
c No control on the coincidence of Bragg angles:               |  |  |
                 QB_ = 0.                                     !|  |  |
c New crystal, but same reflection & wavelength:               |  |  |
                 mode_x0h = 1                                 !|  |  |
                 iSyngony = 0   !any syngony!                  |  |  |
                 nreflect = 1                                 !|  |  |
                 rho      = 0.                                !|  |  |
                 Call Get_Xh97 (Wave, Code(j,l), rho,         !|  |  |
     *                       iSyngony, nreflect,              !|  |  |
     *                       indices, QB_, mode_x0h,          !|  |  |
     *                       xr0, xi0,                        !|  |  |
     *                       xrh, xih, xf_f,                  !|  |  |
     *                       x0_f, xh_f, ifail)               !|  |  |
                 if (iSyngony .le. 1) Then !---+               |  |  |
                    nonCubic = 0              !|               |  |  |
                 else  !-----------------------+               |  |  |
                    nonCubic = 1              !|               |  |  |
                 endif !-----------------------+               |  |  |
c The structure was not found in the database;                 |  |  |
c then it was interpreted as a chemical formula,               |  |  |
c but we do not have the option ot supply density.             |  |  |
c Report this as "no structure in DB":                         |  |  |
                 if (ifail.lt.0    .AND.                      !|  |  |
     *               iSyngony.eq.8 .AND.                      !|  |  |
     *               Abs(rho).lt.(1.E-20)) goto 302 !----------+--+--+---v
                 if (ifail.ne.0) Then !---------------------+  |  |  |
                    i = Index(txt(3),                      !|  |  |  |
     *                'does not exist for this wavelength')!|  |  |  |
                    if (QB_.le.0. .AND. i.gt.0) Then !---+  |  |  |  |
                       goto 306 !-no-reflection-------v  |  |  |  |  |
                    else !-------------------------------+  |  |  |  |
                       goto 130 !-some-other-failure--v  |  |  |  |  |
                    endif  !-----------------------------+  |  |  |  |
                 endif  !-----------------------------------+  |  |  |
c Daa cannot be specified for amorphous structures:            |  |  |
                 if (Abs(Daa(l+1)) .gt. (1.E-20) .AND.        !|  |  |
     *                    iSyngony .eq. 8) goto 112 !----------+--+--+---v
              endif !------------------------------------------+  |  |
              x0_ = x0_ + x0_f*Frac(j,l)                         !|  |
              xh_ = xh_ + xh_f*Frac(j,l)                         !|  |
              xf_ = xf_ + xf_f*Frac(j,l)                         !|  |
              a_  = a_  + a(1)*Frac(j,l)                         !|  |
              p_  = p_  + Poisson*Frac(j,l)                      !|  |
c Controls for automatic da/a calculation:                        |  |
              if (Abs(Daa(l+1)) .ge. (1.E+30)) Then !------+      |  |
c Non-cubic structure:                                     |      |  |
                 if (nonCubic .eq. 1)  goto 105 !--+-------+------+--+---v
c Cubic structure, but                                     |      |  |
c no Poisson data:                                         |      |  |
                 if (Abs(Poisson).lt.(1.E-20)) goto 102   !+------+--+---v
              endif  !-------------------------------------+      |  |
           enddo  !===============================================+  |
                                                                    !|
           if (C0h(l)(1:1) .eq. ' ') Then !---+                      |
              x0(l+1) = x0_                  !|                      |
           endif  !---------------------------+                      |
           if (C0h(l)(2:2) .eq. ' ') Then !---+                      |
              xh(l+1) = xh_                  !|                      |
              xf(l+1) = xf_                  !|                      |
           endif  !---------------------------+                      |
c If both x0 and xh were specified, -- it means that                 |
c the Code was not actually used:                                    |
           if (C0h(l) .eq. 'xx') Then  !------------------+          |
              Code(1,l) = ' '                            !|          |
              Frac(1,l) = 1.                             !|          |
              N_Frac(l) = 1                              !|          |
c No automatic Daa if                                     |          |
c code is unknown:                                        |          |
              if (Abs(Daa(l+1)).ge.(1.E+30)) goto 103 !---+----------+---v
           endif  !---------------------------------------+          |
                                                                    !|
           x0(l+1) = x0(l+1)*W0(l)                                  !|
           xh(l+1) = xh(l+1)*Wh(l)                                  !|
           if (Abs(xh(l+1)) .ge. Abs(x0(l+1))) goto 400 !------------+---v
                                                                    !|
c Automatic caculation of da/a:                                      |
           if (Abs(Daa(l+1)) .ge. (1.E+30)) Then !-----+             |
c Strained substrate (no reference!):                  |             |
c (no longer needed after 2020/05 when we began        |             |
c shifting all da/a by the substrate da/a)             |             |
c             if (Abs(Daa(N_Total+1)) .gt. (1.E-20))  !|             |
c    *                                    goto 104 !---+-------------+---v
              if (a_substrate .le. 0.)    goto 106 !---+-------------+---v
              if (a_          .le. 0.)    goto 106 !---+-------------+---v
              if (Abs(1.-p_) .lt. 0.01)   goto 107 !---+-------------+---v
c Correction of July 9, 1998                           |             |
c (due to Dr. Marius Grundmann from TU Berlin):        |             |
              Daa(l+1)   = (a_-a_substrate)           !|             |
     /                   /   a_substrate              !|             |
c The version below is less accurate because           |             |
c everywhere in the program we multiply da/a           |             |
c by a_substrate:                                      |             |
c ccc         Daa(l+1)   = (a_-a_substrate)           !|             |
c ccc/                   /((a_+a_substrate)/2.)       !|             |
              Daa(l+1)   = Daa(l+1)*(1.+p_)/(1.-p_)   !|             |
              if (Abs(Daa(l+1)) .gt. 0.5) goto 108 !---+-------------+---v
              Aauto(l+1) = 'A'                        !|            !|
           endif  !------------------------------------+            !|
                                                                    !|
        enddo  !=====================================================+
c Added 2011/03/31: restore substrate lattice for Algrebra routines
        Call    LatticeRestore ()
c Corrections for substrate:
        x0(N_Total+1) = x0(N_Total+1)*W0(N_Total)
c We must have non-zero absorption in the substrate!
        if (Abs(Aimag(x0(N_Total+1))) .lt. (1.E-20)) goto 114  !------v
        xh(N_Total+1) = xh(N_Total+1)*Wh(N_Total)
        C0h(N_Total)  = 'ss'
        if (iirezv .eq. 1)  Then  !----------+
           progname = progrezv(istackrezv)  !|
           istackrezv = istackrezv-1        !|
        endif  !-----------------------------+
c Additional test for x0/xh data consistence:
        xh_abs_sg_max = 0.
        xh_abs_pi_max = 0.
        do      l=1,N_Total  !========================================+
           x0_abs    = Abs(x0(l+1))                                  !|
           xh_abs_sg = Abs(xh(l+1))                                  !|
           xh_abs_pi = xh_abs_sg * cos2QB                            !|
           if (xh_abs_sg .gt. xh_abs_sg_max) xh_abs_sg_max=xh_abs_sg !|
           if (xh_abs_pi .gt. xh_abs_pi_max) xh_abs_pi_max=xh_abs_pi !|
c ??? Should we also control for xh = x0 ???                         !|
           if (xh_abs_sg .ge. x0_abs) goto 400 !----------------------+--v
           if (Abs(Aimag(xh(l+1))) .ge.                              !|
     *         Abs(Aimag(x0(l+1)))) goto 401 !------------------------+--v
        enddo  !======================================================+
        if (xh_abs_sg_max .lt. xh_weak_limit) goto 300 !-----------------v
        if (xh_abs_pi_max .lt. xh_weak_limit) Then !--------+
c If pi-polarization is involved, also check if we have     |
c polarization-forbidden reflection:                        |
c (polarization states: 0/3-mixed, 1-sigma, 2-pi)           |
           if     (iPol .eq. 3) Then !-----------------+    |
              iPol = 1  !change from mixed to sigma    |    |
           elseif (iPol .eq. 2) Then !-----------------+    |
              goto 308 !--nothing-to-do-for-pi-only----+----+------------v
           endif  !------------------------------------+    |
        endif  !--------------------------------------------+

c-----------------------
c Reference interface may not be larger than the number of interfaces:
        if (i_standing .ne. 0 .and.
     *      i_reference .gt. N_Total-1) goto 125

        return
c=================================================================
c                           E R R O R S
c-----------------------------------------------------------------
  99    continue
        write   (txt,117)  ProgramVer(1:ipv), InpFile(1:linp)
  117   format  (a,': cannot open input file:'//
     *           a//
     *          'File not found in current directory.')
        ifail = 5
        goto 130

  90    continue
        i = Max(Len_Trim(wrk),1)
        write   (txt,190)  ProgramVer(1:ipv), InpFile(1:linp), line
  190   format  (a,': error reading input file:'//
     *           a//
     *          'At the file line No: ',i4,', while reading:'/)
        ifail = 7                                               !+txt(7)
        goto 120

  80    continue
        i = Max(Len_Trim(wrk),1)
        write   (txt,180)  ProgramVer(1:ipv), InpFile(1:linp), line
  180   format  (a,': no data in input file:'//
     *           a//
     *          'At the file line No: ',i4,', while reading:'/)
        ifail = 7                                               !+txt(7)
        goto 120

  100   continue
        i = Max(Len_Trim(wrk),1)
        if (i .gt. 72) i=72
        write   (txt,118)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line, wrk(1:i)
  118   format  (a,': Incorrect syntax in file:'/
     *           a/
     *          'At the file line No: ',i4,', line content:'/
     *          '[',a,'],'/
     *          'while reading:'/)
        ifail = 7                                               !+txt(7)
        goto 120

  101   continue
        i = Max(Len_Trim(wrk),1)
        if (i .gt. 72) i=72
        write   (txt,119)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line, wrk(1:i)
  119   format  (a,': Parameter(s) out of range in file:'/
     *           a/
     *           'line No:',i4/
     *           'Line=[',a,']'/
     *           ' -- while reading:'/)
        ifail = 7                                               !+txt(7)
        goto 120

  123   continue
        i = Max(Len_Trim(wrk),1)
        write   (txt,223)  ProgramVer(1:ipv), wrk(1:i), InpFile(1:linp)
  223   format  (a,': Input file version: ',a/
     *          'does not match program version'/
     *          'while reading file:'/
     *          a)
        ifail = 4
        goto 130

  125   continue
        write   (txt,135)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     i_reference, N_Total-1
  135   format  (a,': contradictory input in:'//
     *  a//
     *  'Requested reference interface [',i6,'] for standing waves'/
     *  'exceeds the maximum interface index in the structure [',i6,']')
        ifail = 6
        goto 130

  109   continue
        write   (txt,209)  ProgramVer(1:ipv), InpFile(1:linp)
  209   format  (a,': Inconsistent parameters in file:'//a//
     *          'Multiple scan points with zero scan range'/
     *          'while reading:')
        ifail = 7                                               !+txt(7)
        goto 120

  120   continue
        if (ifail .lt. 0)  ifail = 0
        if (ifail .gt. 18) ifail = 18
        txt(ifail+1) = ' '
        txt(ifail+2) = 'Check your input file!'
        ifail  = ifail+2
        goto 130

  111   continue
        i = Max(Len_Trim(wrk),1)
        if (i .gt. 72) i=72
        write   (txt,119)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line, wrk(1:i)
        txt(8) = ' '
        if (Mode_Scan .le. 6) Then !-------------------------------------+
           write (txt(9),219)  Scan_limits(1)*UnitCoef(3)/AngDegree,    !|
     *                         Scan_limits(2)*UnitCoef(3)/AngDegree,    !|
     *                         'Degr.'                                  !|
  219      format ('-- your input reads: [',f7.0,' :',f7.0,'] ',a)      !|
           txt(10)=' '                                                  !|
           txt(11)='Please note that even when you request the Bragg'   !|
           txt(12)='curve to be plotted as a function of the incidence' !|
           txt(13)='or exit angle, the scan range must be specified'    !|
           txt(14)='as the angular deviations of the incident wave from'!|
           txt(15)='the Bragg angle along the direction of scan axis.'  !|
           ifail  = 15                                                  !|
        else !-----------------------------------------------------------+
           write (txt(9),219)  Scan_limits(1),Scan_limits(2),'eV'       !|
           ifail  = 9                                                   !|
        endif !----------------------------------------------------------+
        goto 130

  765   continue
        write   (txt,753)  ProgramVer(1:ipv)
  753   format  (a,':'//'Substrate code is not specified')
        ifail = 3
        goto 130

  131   continue
        write   (txt,132)  ProgramVer(1:ipv)
  132   format  (a,':'//
     *  'X-ray wavelength was chosen to be specified via'/
     *  'characteristic X-ray line, but no line given.')
        ifail = 4
        goto 130

  300   continue
        i = Max(Len_Trim(Code(1,N_Total)),1)
        url = Get_Server_URL()
        lurl = Max(Len_Trim(url),1)
        write   (txt,754)  ProgramVer(1:ipv), indices,
     *                     Code(1,N_Total)(1:i),
     *                     xh_abs_sg_max,
     *                     100.*xh_abs_sg_max/abs(xr0_Sub),
     *                     xh_weak_limit, xh_weak_percent,
     *                     url(1:lurl),
     *                     xh_weak_percent
  754   format  (a,':'//
     *  'Reflection (',3i4,') is forbidden or very weak for crystal'/
     *  '[',a,']: |xh_max| =',g10.3,' =',f5.2,'%*|x0_max|'/
     *  'which is below the lower limit =',g10.3,' =',f5.2,'%*|x0_max|'/
c    *  'The above limit is max[(1E-8)sin(QB),(1E-3)|x0_max|]')         !discarded 2003/06/30
c    *  'The above limit is max[(1E-8)*sin(QB) and ',f5.2,'*|x0_max|]'/ !discarded 2020/05/20
     *  '_____________'//
     *  'Forbidden reflections are those where the reflected'/
     *  'intensity is zero due to crystal structure symmetry.'/
     *  '_____________'//
     *  'You can search for strong reflections for given crystal &',
     *                                    1x,'x-ray wavelength at:'/
     *  a,'cgi/www_form.pl?template=x0p_form.htm'/
     *  'by specifying Minimum |xh/x0| =',f5.2,'% or higher.')
c This automation would be nice, but it needs a buffer much longer than 80 chars:
c    *  a,'cgi/x0p_form.pl?','xway=1&wave=',a,'&code=',a,'&prcmin=',a,
c    *  '&hkl11=-5&hkl12=-5&hkl13=-5&hkl21=5&hkl22=5&hkl23=5',
c    *  '&qb1=0&qb2=90&q1=0&q2=180&base1=1&base2=0&base3=0',
c    *  '&df1df2=-1&modesearch=3')
        ifail = 14
        goto 130

  308   continue
        i = Max(Len_Trim(Code(1,N_Total)),1)
        write   (txt,309)  ProgramVer(1:ipv), indices,
     *                     Code(1,N_Total)(1:i),
     *                     xh_abs_pi_max,
     *                     xh_weak_limit, xh_weak_percent,
     *                     QB
  309   format  (a,':'//
     *  'Reflection (',3i4,') is forbidden or very weak for crystal'/
     *  '[',a,'] and pi-polarized X-rays: |xh_max(pi)|=',g10.3,' ,'/
c    *  'The above limit is max[(1E-8)sin(QB),(1E-3)|x0_max|].'/        !change of 2003/06/30
c    *  'The above limit is max[(1E-8)*sin(QB) and 1%*|x0_max|].'/      !change of 2020/07/13
     *  'which is below the lower limit =',g10.3,' =',f5.2,'%*|x0_max|'/
     *  'Perhaps the Bragg angle QB=',f5.2,' is too close to 45 degr.'/
     *  'Try the same reflection for sigma-polarized x-rays.')
        ifail = 7
        goto 130

  400   continue
        write   (txt,755)  ProgramVer(1:ipv), l,
     *                     Abs(xh(l+1)), Abs(x0(l+1))
  755   format  (a,':'//
     *  'Non-physical input for layer',i5,':'/
     *  '|xh|=',g10.4,' exceeds or equals |x0|=',g10.4)
        ifail = 4
        if (l .eq. N_Total) Then !----------------------------+
           ifail = ifail+1                                   !|
           txt(ifail) = 'NOTE: this is the substrate layer.' !|
        endif !-----------------------------------------------+
        goto 130

  401   continue
        write   (txt,756)  ProgramVer(1:ipv), l,
     *                     Abs(Aimag(xh(l+1))), Abs(Aimag(x0(l+1)))
  756   format  (a,':'//
     *  'Non-physical input for layer',i5,':'/
     *  '|Im(xh)|=',g10.4,' exceeds or equals |Im(x0)|=',g10.4)
        ifail = 4
        if (l .eq. N_Total) Then !----------------------------+
           ifail = ifail+1                                   !|
           txt(ifail) = 'NOTE: this is the substrate layer.' !|
        endif !-----------------------------------------------+
        goto 130

  114   continue
        write   (txt,214)  ProgramVer(1:ipv)
  214   format  (a,':'//'Non-physical input: zero absorption in the',1x,
     *          'substrate!')
        ifail = 3
        goto 130

  302   continue
        i = Max(Len_Trim(Code(mm1,mm2)),1)
        write   (txt,303)  ProgramVer(1:ipv), Code(mm1,mm2)(1:i)
  303   format  (a,':'//
     *          'Crystal [',a,'] not found in the database.')
        ifail = 3
        goto 130

  304   continue
        i = Max(Len_Trim(Code(mm1,mm2)),1)
        write   (txt,305)  ProgramVer(1:ipv), Code(mm1,mm2)(1:i)
  305   format  (a,':'//
     *  'Substrate code [',a,'] corresponds to amorphous material.'//
     *  'Cannot build diffraction geometry for amorphous substrates!')
        ifail = 5
        goto 130

  306   continue
c We get here after Get_Xh97 reported QB=0.
        d_spacing = sngl(1.0D0/DVecMod(indices,3))
        Call  Pakint  (indices,3,wrk,m)
        i = Max(Len_Trim(Code(mm1,mm2)),1)
        write   (txt,307)  ProgramVer(1:ipv),
     *                     wave,
     *                     Code(mm1,mm2)(1:i),wrk(1:m),
     *                     wave, 2*d_spacing, wave/(2*d_spacing)
  307   format  (a,':'//
     *  'The Bragg condition for given X-ray wavelength=',g11.5,'A and'/
     *  a,' crystal with (',a,') Bragg planes cannot be met:'/
     *  'L=',g11.5,' & 2d=',g11.5,' result in sin(QB)=L/(2d)=',g11.5,
     *                                                         ' > 1.'//
     *  'Please use a shorter X-ray wavelength or different Bragg',1x,
     *                                                        'planes!')
        ifail = 7
        goto 130

  112   continue
        i = Max(Len_Trim(Code(j,l)),1)
        write   (txt,212)  ProgramVer(1:ipv), l, Code(j,l)(1:i)
  212   format  (a,':'//
     *          'Cannot use da/a for layer ',i4/
     *          'since the layer code ',a/
     *          'corresponds to amorphous material'//
     *          'Please remove da/a from the profile.')
        ifail = 7
        goto 130

  102   continue
        i = Max(Len_Trim(Code(j,l)),1)
        write   (txt,202)  ProgramVer(1:ipv), l, Code(j,l)(1:i)
  202   format  (a,':'//
     *          'Cannot automatically calculate da/a for layer ',i4/
     *          'because of zero (unknown) Poisson factor in the DB'/
     *          'file COORD.X0H for crystal ', a, ' .'//
     *          'Please specify da/a manually. Note that this version'/
     *          'of the program is limited to |da/a| less than 0.5.')
        ifail = 8
        goto 130

  103   continue
        write   (txt,203)  ProgramVer(1:ipv), l
  203   format  (a,':'//
     *          'Cannot automatically calculate da/a for layer ',i4/
     *          'since the layer is specified by x0 and xh.'//
     *          'Please specify da/a manually. Note that this version'/
     *          'of the program is limited to |da/a| less than 0.5.')
        ifail = 7
        goto 130

c 104   continue
c       write   (txt,204)  ProgramVer(1:ipv)
c 204   format  (a,':'//
c    *          'Cannot automatically calculate da/a'/
c    *          'since the substrate is specified as strained.'//
c    *          'Please specify da/a manually. Note that this version'/
c    *          'of the program is limited to |da/a| less than 0.5.')
c       ifail = 7
c       goto 130

  105   continue
        i = Max(Len_Trim(Code(j,l)),1)
        m = Max(Len_Trim(Symmetry(iSyngony+1)),1)
        write   (txt,205)  ProgramVer(1:ipv), l, Code(j,l)(1:i),
     *                     Symmetry(iSyngony+1)(1:m)
  205   format  (a,':'//
     *          'Cannot automatically calculate da/a for layer ',i4/
     *          'since the layer code ',a,' is not a cubic crystal'/
     *          '(the code is described as ',a,').'//
     *          'Please specify da/a manually. Note that this version'/
     *          'of the program is limited to |da/a| less than 0.5.')
        ifail = 8
        goto 130

  115   continue
        i = Max(Len_Trim(Code(1,N_Total)),1)
        m = Max(Len_Trim(Symmetry(iSyngony_substrate+1)),1)
        write   (txt,215)  ProgramVer(1:ipv), l, Code(1,N_Total)(1:i),
     *                     Symmetry(iSyngony_substrate+1)(1:m)
  215   format  (a,':'//
     *          'Cannot automatically calculate da/a for layer ',i4/
     *          'since the substrate code ',a,' is not cubic crystal.'/
     *          '(the code is described as ',a,')'//
     *          'Please specify da/a manually. Note that this version'/
     *          'of the program is limited to |da/a| less than 0.5.')
        ifail = 8
        goto 130

  106   continue
        write   (txt,206)  ProgramVer(1:ipv), l
  206   format  (a,':'//
     *          'Cannot automatically calculate da/a for layer ',i4/
     *          'since one of lattice constants is zero.'//
     *          'Please specify da/a manually. Note that this version'/
     *          'of the program is limited to |da/a| less than 0.5.')
        ifail = 7
        goto 130

  107   continue
        write   (txt,207)  ProgramVer(1:ipv), l
  207   format  (a,':'//
     *          'Cannot automatically calculate da/a for layer ',i4/
     *          'since the Poisson ratio is close to 1.'//
     *          'Please verify the Poisson ratio in COORD.X0H.')
        ifail = 6
        goto 130

  108   continue
        i = Max(Len_Trim(Code(1,N_Total)),1)
        write   (txt,208)  ProgramVer(1:ipv), l, Daa(l+1),
     *                     Code(1,N_Total)(1:i), a_substrate
  208   format  (a,':'//
     *          'Automatic calculation of da/a for layer ',i4/
     *          'provides suspiciously big value da/a=',g10.3/
     *          'with respect to a_substrate(',a,')=',g10.3,'A'//
     *          'Please specify da/a manually. Note that this version'/
     *          'of the program is limited to |da/a| less than 0.5.')
        ifail = 8
        goto 130

  110   continue
        write   (txt,210)  ProgramVer(1:ipv)
  210   format  (a,':'//
     *          'Please enter non-zero integer Miller indices',
     *          1x,'to specify OTHER scan axis!')
        ifail = 3
        goto 130

  113   continue
        write   (txt,213)  ProgramVer(1:ipv)
  213   format  (a,':'//
     *          'Mono convolution option is not implemented yet!')
        ifail = 3
        goto 130

  130   continue
        if (modebat .eq. 1) Then  !---------------+
           do  i=1,Iabs(ifail) !===============+  |
              l = Max0(Len_Trim(txt(i)),1)    !|  |
              write (*,'(1x,a)') txt(i)(1:l)  !|  |
           enddo  !============================+  |
        endif  !----------------------------------+
        Call Message (txt,Iabs(ifail),2)
        goto 28

  28    continue
        Call exit_quiet()
        End
