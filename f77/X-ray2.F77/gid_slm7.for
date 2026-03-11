        Program Gid_slm7
c       USE MSFLIB                                      !for QuickWin
c ---------------------------------------------------------
c // If we compile gid_slm7 as a QuickWin application,   //
c // then it can be called from a console application    //
c // using the following comment:                        //
c //                                                     //
c //    str='gid_slm7 xxx.inp'                           //
c //    i = System (str)                                 //
c // Then, the console immediately gets the control back //
c // and it can watch the progress of gid_slm7.          //
c // See also: i = SETEXITQQ(2) (forcing Qwin to close   //
c // window after the application is finished)           //
c ---------------------------------------------------------
c Subroutines: GIDin_m7.for + RFMm7_rx.for + RFMs_urd.for + RFMs_sss.for
c Input file:  gid_slm7.inp
c ---------------------------------------------------------
c   X-ray grazing incidence diffraction in multilayers
c
c                   By S.A.Stepanov
c
c *** Version M7 (i.e. Millennium-7, 2007) optimized for WWW access ***
c
c 1.  Perfect SLs with overlayer.
c 2.  Transition layers included.
c 2.1 Large-scale surface roughness included.
c 3.  Small-Scale interface roughness included.
c 3.1 Small-Scale surface roughness included,
c     Large-scale surface roughness modified.
c 3.2 Small-Scale surface roughness corrected by (1/2),
c     possibility to change exponent index excluded.
c 4.  Special version for EXTREMELY ASYMMETRIC cases (XEAD).
c 4.1 Normal lattice strain included,
c     Exit angle scan replaced by incident beam scan
c 4.2 Large angular deviations from Bragg condition
c     allowed.
c 4.3 Reading of input files transferred to a separate
c     subroutine gid_4inp.for,
c     Control for Laue cases added
c 4.4 Normal strain now generalized for XEAD and GID
c 5.0 Simplified in order to calculate the x-ray wavefields
c     needed for diffuse scattering and secondary emission:
c     large-scale roughness and transition layers excluded
c 5.5 Simplified: SL excluded (for Petrashen)
c 97. Combined advantages of 5.0 and 5.5: script language for
c     top-layer profile, allowing a specification of periodic
c     structures; generalized scanning; improved batch
c     ("silent") mode for using with CGI-interfaces.
c 99.1 Added interface to Henke and Brennan-Cowan interfaces.
c     (this version is also known as gid_sl98; active till 2000/08/28).
c 99.2 Added support for high Bragg angles QB=90.
c M1. (or Millennium-1, 2001) Added comment lines; added support for
c     four more ways to specify coplanar diffraction geometry; added
c     hooks in input file for standing waves and monochromator
c     convolution (not implemented yet!)
c M5. (or Millennium-5, 2005) Added Auto-DB selection; added Alpha-Max
c     for "crystal" layer threshold.
c M7. (or Millennium-7, 2007) Added dynamic selection of scattering matrix
c     size (like extremely asymmetric in v.4); added phases of structural
c     amplitudes (Kaganer's suggestion), implemented standing waves.
c -------------------------------------------------------
c+===============================+
        Include 'gid_slm7.inc'  !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter       (N_Top_Max   = 1001,
c cc *                   N_Total_Max = N_Top_Max+1)

        Integer         iarg
        Integer         N_Frac_Max
        Parameter      (N_Frac_Max  = 4)

        Complex*16      fh_cmplx

        Real*8          Rhs, R0s, Rh,
     *                  SWSUM(N_Stand_Max),     !M7: summ of SW intensity
     *                  DisorPlanes8, h,
     *                  k0_vector(3), k0, k00,
     *                  Gam0, SiQB,
     *                  Scan_vector(3),
     *                  dk_vector(3), dk,
     *                  fh2, fh, fhi, Alpha,
     *                  b_vector(3), b_len,
     *                  c_vector(3), c_len,
     *                  z8,                     !M7:
     *                  cp(2), ttt(12), ooo(12),
     *                  eps /1.D-5/

        Real*8          DVecSca2, DVecMod2, RFM8_m7
        External        DVecSca2, DVecMod2, RFM8_m7

        Real*4          Ratio, Energy, dE, Energy0, Wave0,
     *                  Frac(N_Frac_Max,N_Total_Max),
     *                  QB, QBr, Pfi, fcentre,
     *                  TER_Angle, PsiEff, d_spacing,
     *                  f0_centre, fh_centre,
     *                  Daa_Top_mean, Sigma_Top_mean,
     *                  Asym, Exti, c1, c2, Start,
     *                  x0_max, threshold, TER_threshold,  !M7:
     *                  progress

        Integer*4       hh, mm, ss, bias

        Integer         indices(3), iii(12),
     *                  N_Frac(N_Total_Max),
     *                  iPol, MinPol, MaxPol,
     *                  narg, igie, kkk,
     *                  iscan, iasci,
     *                  i, j, l, io, lout, ltop,
     *                  line, nSc,
     *                  nerror, ierrcount, iskipped,
     *                  isico, iwarn, lelab,
     *                  m_0, m_h, neg_angles,   !M7:
     *                  m_type_previous,        !M7: matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s.
     *                  iostts

        Real            wave2energy, Bragan
        External        wave2energy, Bragan

        Character       Code(N_Frac_Max,N_Total_Max)*20,
     *                  C0h(N_Total_Max)*2,
     *                  Aauto(N_Total_Max+1)*1,
     *                  DoneFile*80, timestamp*9,
     *                  Version(2)*80, scan_label*6,
     *                  wrk*80

        Character       ScanName(8)*22 /
     *                    'Round N_surface',
     *                    'Round [k0*N_surface]',
     *                    'Round h',
     *                    'Round [k0*h]',
     *                    'Round other axis',
     *                    'Takeoff spectrum (PSD)',
     *                    'Energy(eV), x0h recalc',
     *                    'Energy(eV), x0h fixed'/
        Character       Geometry(9)*38 /
     *                    'surface & incidence angle [noncopl.]',       !1
     *                    'surface & exit angle of kh [noncopl.]',      !2
     *                    'surface & coplanar grazing incidence',       !3
     *                    'surface & coplanar grazing exit of kh',      !4
     *                    'surface & symmetric Bragg case',             !5
     *                    'coplanar & Bragg planes disorientation',     !6
     *                    'coplanar & incidence angle',                 !7
     *                    'coplanar & exit angle of kh',                !8
     *                    'coplanar & beta=g0/|gh|'/                    !9
        Character       Spm(3)*5 /'Sigma',' Pi  ','Mixed'/
        Character       Uni(6)*5  /'degr.','min.',
     *                             'mrad.','sec.',
     *                             'urad.', 'eV'/
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
        Character       YesNo(2)*3 /'No','Yes'/
        Character       MatrixMode(3)*20/'Always 4x4',
     *                                   'Prescan reduction',
     *                                   'Auto-reduce @ scan'/
        Character       MatrixType(5)*24/'Full 4x4',                    !0
     *                                   'Grazing incidence 3x3',       !1
     *                                   'Grazing exit 3x3',            !2
     *                                   'Normal diffraction 2x2',      !3
     *                                   'Specular reflection 2x2'/     !4
        Logical*4       FileExist
        External        FileExist

        Logical         Failure

        Integer*4       nargs
        External        nargs

        Character       radiat*6
        Common  /x0pa1/ radiat

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
     *             standing_phase,              !M7: SW phase in units of pi
     *             xh_weak_limit,               !M5: the weakest possible xh
     *             alpha_max                    !M5: max(alpha/xh) to treat crystal as "crystal"
        Integer    N_Top,                       !top layer sublayers
     *             N_Total,                     !total of layers
     *             i_standing,                  !M7: 1=yes/0=no SW
     *             i_reference,                 !M7: reference interface for SW (0=surface)
     *             n_standing,                  !M7: Number of SW depth points
     *             m_reduction,                 !M7: matrix size reduction flag (0/1/2)
     *             m_type,                      !M7: matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s.
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
        Real*8          Positions(N_Total_Max), !depth of interfaces
     *                  Phase(N_Total_Max),     !int(d(h_z*dz))
     *                  SW(N_Stand_Max),        !M7: Standing waves intensity
     *                  f0,                     !incident angle
     *                  Scan_angle,             !scan angle
     *                  R0                      !reflection coefficient
        Real*4          Psi,                    !misorientation angle
     *                  Pol(2)                  !polarization factors
        Logical         Account_Sigma
        Integer         N_Used,                 !number of used layers
     *                  jPol,                   !polariz.index
     *                  ifail                   !failure code
c This is common between
c gid_sl and RFM only:
        Common /RFM7dat/ R0, Positions, Phase,
     *                   SW,                    !M7: new
     *                   f0, Scan_angle,
     *                   Psi, Pol,
     *                   Account_Sigma,
     *                   N_Used, jPol, ifail
c----------------------
c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer             iHenkeCowan
        Common /HenkeCowan/ iHenkeCowan

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        istackrezv = 0
        progname = 'gid_slm7'
c       i = SetExitQQ(2)        !QWIN$EXITNOPERSIST=2   !for QuickWin
c --------------------------------------------------------
        i     = N_Pts_Max                       !keep GNU Fortran happy
        b_len = 0.                              !keep GNU Fortran happy
        c_len = 0.                              !keep GNU Fortran happy
        fh_cmplx = (0.,0.)                      !keep GNU Fortran happy

        ipv = Len (ProgramVer)
        ProgramVer = progname(1:ipv)
        ipv = Len_Trim (ProgramVer)
        Call CaseLat (ProgramVer(1:ipv),1)      !convert to lower case
c This is preliminary (until we've figured out the output filename):
        ErrorFile = ProgramVer(1:ipv)//'.err'
        Call DeleteErrorFile()

        write   (Version,111)
  111   format  ('GID_sl: X-Ray Diffraction from Strained ',
     *  'Crystals at Usual and Grazing Incidence.'/
     *  'By S.Stepanov <sstepanov@anl.gov>',26x,
c    *  'Version M7: Jan, 2007')
     *  'Version M7: Jun, 2020')

        narg = Nargs()-1
        if (narg .eq. 0) Then !---------------------+
                                                   !|
           InpFile = ProgramVer(1:ipv)//'.inp'     !|
                                                   !|
        else  !-------------------------------------+
                                                   !|
           iarg = 1                                !|
           Call get_arg_wrapper (iarg,InpFile)     !|
                                                   !|
        endif  !------------------------------------+

        do i=1,N_Stand_Max !==+
           SW(i) = 0.        !|
        enddo !===============+
        modebat = 1                     !assume batch-mode by default
        iBatch_mode = modebat

c #######################################################
c              Read input file
c #######################################################

        Call Read_GID_InputM7 (indices,
     *                         DisorPlanes8, Pfi, igie,
     *                         QB, fcentre, iPol,
     *                         Code, Frac, N_Frac,
     *                         N_Frac_Max, C0h, Aauto)
c#######################################################
        Pol(1)    = 1.
        Pol(2)    = Abs(cos(2.*QB*AngDegree))
c 2018.10: changed xabs from |x0_substr| to |x0_max|
c       xabs      = abs(x0(N_Total+1))
c       TERang    = sqrt(xabs)
c       TER_Angle = TERang / UnitCoef(2)
c       if (Abs(xabs) .lt. (1.E-24)) xabs = 1.
        Call xabsmax (x0, N_Total_Max, N_Total, UnitCoef(2),
     *                xabs, TERang, TER_angle)

        if (Mode_Scan .le. 6) Then !------+
           scan_label = 'angle'          !|
        else  !---------------------------|
           scan_label = 'energy'         !|
        endif !---------------------------+
        lelab = Len_Trim(scan_label)

c #######################################################
c            Build the incident wave vector:
c #######################################################
        k0   =  2.* pi/wave
        Call Build_k0 (Surface_Normal, DisorPlanes8,
     *                 h_vector, h, Pfi, Psi, PsiEff,
     *                 k0_vector, k0, Gam0, QB, SiQB,
     *                 fcentre, f0_centre, fh_centre,
     *                 TER_Angle, UnitCoef(2), AngDegree,
     *                 nBaseNormal, nReperDirection, igie,
     *                 ProgramVer, ipv, ifail)
        if (ifail .ne. 0) goto 28

c #######################################################
c       Determine the direction of angular scans:
c #######################################################
        goto (601,602,603,604,605,606,607,607)  Mode_Scan
c------
c Round internal normal n:
  601   continue
        Call DUnitVec2 (Surface_Normal,
     *                  Scan_vector,3)
        goto 610
c------
c Round [k0*n]:
  602   continue
        Call DVecVec (k0_vector,
     *                Surface_Normal,
     *                Scan_vector,3)
        goto 610
c------
c Round h:
  603   continue
        Call DUnitVec2 (h_vector,
     *                  Scan_vector,3)
        goto 610
c------
c Round [k0*h]:
  604   continue
        Call DVecVec (k0_vector,
     *                h_vector,
     *                Scan_vector,3)
        goto 610
c------
c Round other -- specified -- axis:
  605   continue
        do i=1,3 !==================================+
           Scan_vector(i) = iOther_Scan_Vector(i)  !|
        enddo  !====================================+
        goto 610
c------
c PSD analysis of diffracted waves over their takeoff angles
c (uncollimated beam):
  606   continue
        Invert_Axis = 0
        goto 616        !-----------------------------------+
c------                                                     v
c Energy scan in eV: 7=x0h recalc, 8=x0h fixed
  607   continue
        Call    DUnitVec2 (k0_vector,
     *                     Scan_vector,3)
        if (Mode_Scan .eq. 8 .AND.              !cannot use same x0h
     *  (f0_centre .lt. 10*TER_Angle .OR.       !in grazing cases
     *   fh_centre .lt. 10*TER_Angle)) Then !---+
           Mode_Scan = 7                       !|
        endif  !--------------------------------+
        Invert_Axis=0
        goto 616        !-----------------------------------+
c------                                                     V
  610   continue
        if (Invert_Axis .eq. 1) Then !-------------+
           do i=1,3 !===========================+  |
              Scan_vector(i) = -Scan_vector(i) !|  |
           enddo  !=============================+  |
        endif  !-----------------------------------+

        Call DUnitVec2 (Scan_vector,Scan_vector,3)

c------
c      ->
c Let  a = unit vector along the scan axis (Scan_vector).
c                    ->    ->    ->
c        We expand:  dk = x*c + y*b,
c        ->  -> ->       ->  -> -> ->    ->  -> -> ->
c where  c = [k*a]  and  b = [a*[k*a]] = k - (k*a)*a
c                                                         ->
c -- are the two vectors, perpendicular to each other and a .
c
c It is easy to find that:
c               x = (b/c)*sin(dQ)
c               y = - 2*sin^2(dQ/2)
c                   ->    ->
c Now, we calculate b and c:

        Call DProject2 (Scan_vector,k0_vector,b_vector,3)
        Call DVecVec   (k0_vector,Scan_vector,c_vector,3)
        b_len = DVecMod2 (b_vector,3)
        c_len = DVecMod2 (c_vector,3)

c------                                                     |
c Depending on Exit angle:                                  V
  616   continue  !-----------------------------------------+

c #######################################################
c  Contributions of Polarizations in the incident beam:
c #######################################################
        goto (21,22,23) iPol !-->---+
  21    continue   !<--------------<|
        Ratio  = 1.e+10            !|
        cp(1)  = 1.                !|
        cp(2)  = 0.                !|
        MinPol = 1                 !V
        MaxPol = 1                 !|
        goto 24    !---------->-----+---------+
  22    continue   !<--------------<|         |
        Ratio  = 1.e+10            !|         |
        cp(1)  = 0.                !|         |
        cp(2)  = 1.                !V         V
        MinPol = 2                 !|         |
        MaxPol = 2                 !|         |
        goto 24    !---------->-----+---------|
  23    continue   !<--------------<+         |
        Ratio = abs(cos(2.*QB*AngDegree))    !|
        if (Mode_Scan.ne.6) Then !---------+  |
           cp(1) =   1.  / (1.+Ratio)     !|  V
           cp(2) = Ratio / (1.+Ratio)     !|  |
        else  !----------------------------|  |
c At fh-scan (Position Sensitive Detector, |  |
c PSD), then mono collimates incidence     |  |
c angle (i.e. it is rotated by 90 degr.):  |  |
           cp(1) = Ratio / (1.+Ratio)     !|  |
           cp(2) =   1.  / (1.+Ratio)     !|  V
        endif  !---------------------------+  |
                                             !|
        MinPol = 1                           !|
        MaxPol = 2                           !|
  24    continue   !<---------------------<---+

c#######################################################
c Reflection asymmetry:
        c1 = Amax1(f0_centre,TER_Angle)
        c1 = sin (c1*UnitCoef(2))
        c2 = Amax1(fh_centre,TER_Angle)
        c2 = sin (c2*UnitCoef(2))
        if (abs(c1) .lt. 1.E-32)  c1 = 1.E-32
        if (abs(c2) .lt. 1.E-32)  c2 = 1.E-32
        Asym = c1/c2
c Extinction length (actually we need |X_{h}*X_{-h}|):
        if (Abs(xh(N_Total+1)) .gt. 0.) Then !-----------------+
           Exti = Wave*Sqrt(c1*c2) / (pi*Abs(xh(N_Total+1)))  !|
           do l=1,N_Total !======================+             |
              if (Sigma(l) .gt. Exti) goto 162 !-+-------------+----+
           enddo  !==============================+             |    v
        else !-------------------------------------------------+
           Exti = 0.                                          !|
        endif !------------------------------------------------+

c #######################################################
c M7: calculate threshold g0,gh for matrix reduction:
        x0_max = 0.
        do  l = 1, N_Total  !================+
           c1 = Abs(x0(l+1))                !|
           if (c1 .gt. x0_max) x0_max = c1  !|
        enddo !==============================+
c Should we make this constant
c defined by a user?
c       threshold = 5.
        threshold = 10.                     !changed from 5 on 2023.04.26
        TER_threshold = threshold * sqrt(x0_max)

c #######################################################
c Mean paramters of surface profile:
        Daa_Top_mean   = 0.
        Sigma_Top_mean = 0.
        if (N_Top .gt. 0) Then  !----------------------------+
           do    l = 1, N_Top  !==========================+  |
              Daa_Top_mean   = Daa_Top_mean + Daa(l+1)   !|  |
              Sigma_Top_mean = Sigma_Top_mean + Sigma(l) !|  |
           enddo  !=======================================+  |
           Daa_Top_mean   = Daa_Top_mean   / N_Top          !|
           Sigma_Top_mean = Sigma_Top_mean / N_Top          !|
        endif  !---------------------------------------------+

        Positions(1) = 0.0D0                !the surface
        Phase(1)     = 0.0D0
        do i=2,N_Total !===========================================+
                                                                  !|
           Positions(i) = Positions(i-1) + Thickness(i-1)         !|
                                                                  !|
c ATTENTION: phase is not used in this program at                  |
c all. Moreover, this definition must be verified,                 |
c because other programs use a different one                       |
c d(Phase) = Thickness*h_z                                         |
c          = Thickness*k0*Psi*(1-Daa)                              |
                                                                  !|
           Phase(i) = Phase(i-1) + Thickness(i-1) * Dble(Daa(i))  !|
        enddo  !===================================================+

        Energy = wave2energy(Wave)

c #######################################################
c                Print summary of parameters:
c#######################################################
        lout   = Max(Len_Trim(OutFile),1)
        i      = Min(lout,40)
        ltop   = Max(Len_Trim(File_Top),1)
        j      = Min(ltop,40)
        write (txt,33) OutFile(1:i),
     *                 File_Top(1:j),
     *                 Thickness_Top, N_Top,
     *                 Daa_Top_mean, Sigma_Top_mean,
     *                 Code(1,N_Total), Daa(N_Total+1), Sigma(N_Total),
     *                 Wave, Energy, Radiat,
     *                 Spm(iPol), cp,
     *                 indices,
     *                 Geometry(igie),
     *                 Uni(iUnits(1)+1), nBaseNormal,
     *                           DisorNorm, nReperDirection,
     *                 (180.-DisorPlanes8), -Pfi, Uni(iUnits(2)+1),
     *                 Uni(iUnits(2)+1), QB, TER_Angle,
     *                 Uni(iUnits(2)+1), f0_centre, fh_centre,
     *                 Asym,
     *                 (1.E-4)*Exti,
     *                 ScanName(Mode_Scan),
     *                 Uni(iUnits(3)+1), Scan_limits, nScan
  33    format  (
     *  'Output files [.dat & .tbl] will be...|',a/
     *  'Top: Profile name....................|',a/
     *  'Top: Thickness(A), N_sublayers.......|',f17.2,i12/
     *  'Top: <d(a)/a>, <roughness(A)>........|',g12.5,f6.2/
     *  'Substrate: Code, d(a)/a, roughness(A)|',a,g12.3,f6.2/
     *  'X-ray wavelength (A), energy (keV)...|',g12.7,1x,g12.7,1x,a/
     *  'Polarization, contributions sigma:pi.|',a,6x,g9.3,' : ',g9.3/
     *  'Bragg reflection ....................|(',i3,',',i3,',',i3,')'/
     *  'Parameters defining diffr. geometry..|',a/
     *  'Surf.plane, miscut (',a,'), direction|(',3i3,')  ',g10.3,'  (',
     *                                                          3i3,')'/
     *  'Bragg planes & vector angles to surf.|',g14.7', ',g14.7,1x,a/
     *  'Bragg angle (degr.), TER angle(',a,')|',2g15.7/
     *  'X-ray incident & exit angles  (',a,')|',2g15.7/
     *  'Asymmetry factor b=gamma_0/|gamma_h|.|',g15.7/
     *  'Extinction length (um)...............|',g15.7/
     *  'Scan mode............................|',a/
     *  'Scan angular limits (',a,'), points..|',2g13.5,i6)

        if (igie .gt. 5) Then  !----------------------------------------+
           line = 10                                                   !|
           if (igie .eq. 6) Then !------------------------------------+ |
              write (txt(line),34) 'Bragg planes angle to surface ', !| |
     *                             Uni(iUnits(2)+1), Asymmetry_par   !| |
           elseif (igie .eq. 7) Then !--------------------------------+ |
              write (txt(line),34) 'Incidence angle of K0.........', !| |
     *                             Uni(iUnits(2)+1), Asymmetry_par   !| |
           elseif (igie .eq. 8) Then !--------------------------------+ |
              write (txt(line),34) 'Exit angle of Kh..............', !| |
     *                             Uni(iUnits(2)+1), Asymmetry_par   !| |
           elseif (igie .eq. 9) Then !--------------------------------+ |
              write (txt(line),35) 'parameter beta=|g0/gh|........', !| |
     *                             '.....', Asymmetry_par            !| |
           endif  !---------------------------------------------------+ |
  34       format (a30,'(',a,')| ',g12.5)                              !|
  35       format (a30,'.',a,'.| ',g12.5)                              !|
        endif  !--------------------------------------------------------+

        if (Mode_Scan .eq. 5)  Then  !-----------------------------+
           line = 16                                              !|
           txt(line)(56:56) = '('                                 !|
           Call Pakint (iOther_Scan_Vector,3,txt(line)(57:76),i)  !|
           txt(line)(57+i:57+i) = ')'                             !|
        endif  !---------------------------------------------------+
        if (Invert_Axis .eq. 1) Then !---------------+
           line = 16                                !|
           i = Len_Trim(txt(line))                  !|
           if (i .lt. 72) Then !----------------+    |
              txt(line)(i+2:76) = '/inverted/' !|    |
           else  !------------------------------+    |
              txt(line)(27:36) = '/inverted/'  !|    |
           endif  !-----------------------------+    |
        endif  !-------------------------------------+

        line = 17
c#######################################################

        write   (3,55)  Version
        write   (3,11)
  11    format  (1x,78('='))
        write   (3,55)  (txt(i)(1:Len_Trim(txt(i))),i=1,line)
  55    format  (20(a:/))                       !1x,a

        Call Algorithm_Type (Version(1))
        l = Len_Trim (Version(1))
        i = Len_Trim (DBtext(iHenkeCowan+2))
        write   (3,62)  place(1:1),place(2:12),
     *                  Uni(iUnits(2)+1), PsiEff,
     *                  Version(1)(1:l),
     *                  DBtext(iHenkeCowan+2)(1:i),
     *                  xh_weak_limit,
     *                  MatrixMode(m_reduction+1),
     *                  TER_threshold/UnitCoef(2), Uni(iUnits(2)+1),
     *                                             threshold,
     *                  YesNo(1+i_standing),
     *                  alpha_max
  62    format  (
     *  'Columns in output (DAT) file (*).....|',a,'D',a/
     *  ' (*) A=scan_Angle                    |'/
     *  ' (*) D=Diffracted_intensity          |'/
     *  ' (*) E=Exit_angle of diffracted wave |'/
     *  ' (*) I=Incidence_angle               |'/
     *  ' (*) 2=2*Theta: incidence+exit angles|'/
     *  ' (*) S=Specularly_reflected_intensity|'/
     *  'Effec.inward planes misorient.(',a,')|',g15.8/
     *  'Algorithm type (Matrix / Recursive)..|',a/
     *  'Database used for dispersion correct.|',a/
     *  'Weak reflection limit (minimum |xh|).|',g15.8/
     *  'Scattering matrix reduction mode.....|',a/
     *  'Threshold f0,fh for matrix reduction.|',g15.8,a,
     *                                        ' (',f4.1,' x TER_angle)'/
     *  'Standing waves printout..............|',a/
     *  '|alpha_max| for amorph.approximation.|',g15.8)

        if (i_standing .ne. 0) Then !---------------------------+
           write   (3,67)                                      !|
     *                  i_reference, Positions(i_reference+1), !|
     *                  (standing_range(i),i=1,2),             !|
     *                  n_standing                             !|
           if (standing_phase .ge. 0.) Then !---------+         |
              Call PakReal(standing_phase,1,wrk,l)   !|         |
           else  !------------------------------------+         |
              wrk = 'Any, i.e. varying with Z'       !|         |
           endif !------------------------------------+         |
           l = Len_Trim (wrk)                                  !|
           write (3,73) wrk(1:l)                               !|
        endif !-------------------------------------------------+
  67    format  (
     *  'Standing waves: reference interface..|',i5,10x,'(',f13.2,'A)'/
     *  'Standing waves: offsets range (A)....|',g12.5,',  ',g12.5/
     *  'Standing waves: number of depth pts..|',i5)
  73    format  (
     *  'Standing waves: phase/pi at location.|',a)

        write   (3,11)
        write   (3,55)  ' *** User Comment: '
        do i=1,3  !====================================+
           l = Len_Trim (Comment(i))                  !|
           if (l .gt. 0) write (3,55) Comment(i)(1:l) !|
        enddo  !=======================================+
        write   (3,11)

        if (Index(place,'D') .eq. 0) Then !-----+
c Insert diffracted intensity into place-2:     |
           do i=12,3,-1 !===================+   |
              place(i:i) = place(i-1:i-1)  !|   |
           enddo  !=========================+   |
           place(2:2) = 'D'                    !|
        endif  !--------------------------------+
c 1. A -- scan Angle
c 2. D -- Diffracted intensity  [always included at place-2]
c 3. S -- Specular (transmitted) intensity
c 4. B -- Bragg deviation (alpha-parameter)
c 5. I -- Incidence_angle
c 6. E -- Exit angle
c 7. M -- iMaginary part of exit angle
c 8  2 -- 2*theta (the sum of incident & exit angles)
        do i=1,12  !===============================+
           iii(i) = Index("ADSBIEM2",place(i:i))  !|
        enddo  !===================================+

c#######################################################

        write   (3,61)
  61    format  (' The Structure Profile:'/
     *           ' ====================='/
     *  '  Nr   Thickness    Da/a   # Roughness _______Code________',1x,
     *  'Fraction #0h#  W0    Wh',11x,'X0',18x,'Xh',7x,'df(Xrh-Xih)/pi',
     *  3x,'Xh_phase_shift/pi')
        do l=1,N_Total !=====================================+
           write (3,64)  l, Thickness(l),                   !|
     *                   Daa(l+1), Aauto(l+1),              !|
     *                   Sigma(l),                          !|
     *                   Code(1,l), Frac(1,l),              !|
     *                   C0h(l), W0(l), Wh(l),              !|
     *                   x0(l+1), xh(l+1), xf(l+1), xs(l+1) !|
           do j=2,N_Frac(l)  !====================+          |
              write (3,65) Code(j,l), Frac(j,l)  !|          |
           enddo  !===============================+          |
        enddo  !=============================================+
  64    format  (i4,1x,f10.1,e11.3,1x,a1,f7.1,4x,a,f7.3,3x,a,2f6.2,
     *          2(' (',g8.2,',',g8.2,')'),f8.3,f18.3)
  65    format  (36x,                     a,f7.3)
        write   (3,66)
  66    format  (1x,40('-')/
     *  ' xx = X0h data are replaced by x0/xh from profile'/
     *  ' ss = the substrate'/
     *  '  A = automatic calculation of da/a via Poisson ratio')

c Change of 2020/05:
        if (N_Top.gt.0 .AND. Abs(Daa(N_Total+1)).gt.(1.E-20)) Then !---+
           write (3,71) Daa(N_Total+1)                                !|
  71       format (/' Since substrate is specified with da/a=',e11.3,
     *     ', adding it to all top layers'/
     *     ' because their strains are treated as relative substrate.'/
     *     1x,40('-'))                                                !|
           do l=1,N_Top !===============================+              |
              if (Abs(xh(l+1)) .gt. 0.) Then !-------+  |              |
                 c1 = Daa(l+1) + Daa(N_Total+1)     !|  |              |
                 write (3,72) l, Daa(l+1), c1       !|  |              |
  72             format(i4,1x,e12.3, ' => ',e12.3)  !|  |              |
                 Daa(l+1) = c1                      !|  |              |
              endif  !-------------------------------+  |              |
           enddo !======================================+              |
        endif !--------------------------------------------------------+

c#######################################################
c                 Open output files:
c#######################################################
        Call OpenList (OutFile(1:lout)//'.dat',1,ifail)
        if (ifail .ne. 0) goto 28

        if (i_standing .ne. 0) Then !----------------------------------+
           if (n_standing .gt. 1 .and. nScan .gt. 1) Then !---------+  |
              Call OpenList (OutFile(1:lout)//'_sw.grd',2,ifail)   !|  |
              if (ifail .ne. 0) goto 28                            !|  |
              write (2,56) n_standing,        nScan,               !|  |
     *                     standing_range(1), standing_range(2),   !|  |
     *                     Scan_limits(1),    Scan_limits(2),      !|  |
     *                       0.,                0.                 !|  |
  56          Format ('DSAA'/2I6,3(/2G15.7))                       !|  |
           else !---------------------------------------------------+  |
              Call OpenList (OutFile(1:lout)//'_sw.dat',2,ifail)   !|  |
              if (ifail .ne. 0) goto 28                            !|  |
           endif !--------------------------------------------------+  |
        endif !--------------------------------------------------------+

        if (iDebug .ne. 0)  Then  !------------------------+
           Open  (unit=33, file=OutFile(1:lout)//'.deb',  !|
     *                                 status='unknown')  !|
        endif  !-------------------------------------------+

c#######################################################
c                          M7: Prescan:
c#######################################################
        if (m_reduction .eq. 1) Then !-------------------------------------+
           m_0 = 0                       !no grazing incidence             |
           m_h = 0                       !no grazing exit                  |
           Energy0 = Energy                                               !|
           Wave0   = Wave                                                 !|
           k00     = k0                                                   !|
           do nSc=1,nScan !=============================================+  |
              Scan_angle = Scan_limits(1) + dScan_angle*(nSc-1)        !|  |
              if (Mode_Scan .lt. 6) Then !---------------------------+  |  |
                 Call CalcScanAngles(Scan_angle, scan_label,  !input |  |  |
     *                               xabs,                    !input |  |  |
     *                               iUnits, UnitCoef,        !input |  |  |
     *                               k0_vector, k0, h_vector, !input |  |  |
     *                               b_vector, b_len,         !input |  |  |
     *                               c_vector, c_len,         !input |  |  |
     *                               Surface_Normal, Psi,     !input |  |  |
     *                               Alpha, f0, fh,           !-outpt|  |  |
     *                               fh_cmplx,                !-outpt|  |  |
     *                               neg_angles)              !-outpt|  |  |
                 if (neg_angles .ne. 0) goto 225 !----->-----v       |  |  |
              elseif (Mode_Scan.eq.6) Then  !------------------------+  |  |
c Scan over fh    --  Mode_Scan=6:                                   |  |  |
                 f0 = Gam0                                          !|  |  |
                 fh = Dsin (Scan_angle*UnitCoef(3))                 !|  |  |
                 if (fh .le. 0.)        goto 225 !----->-----v       |  |  |
              else  !------------------------------------------------+  |  |
c Energy scan:                                                       |  |  |
                 dE = sngl(Scan_angle*UnitCoef(3))                  !|  |  |
                 Energy = Energy0 + dE                              !|  |  |
                 if (Energy .le. 0)     goto 225 !----->-----v       |  |  |
                 Wave   = wave2energy(Energy)                       !|  |  |
                 if (Mode_Scan.eq.7) Then !-----------------------+  |  |  |
                    Call recalcXh (Wave, indices,                !|  |  |  v
     *                             N_Total, N_Frac_Max,          !|  |  |  |
     *                             N_Frac,                       !|  |  |  |
     *                             Code, Frac,                   !|  |  |  |
     *                             x0, xh, xf,                   !|  |  |  |
     *                             W0, Wh,                       !|  |  |  |
     *                             C0h, QB, ifail)               !|  |  |  |
                    if (ifail .ne. 0)    goto 225 !----->-----v   |  |  |  |
                    Call xabsmax (x0, N_Total_Max, N_Total,      !|  |  |  |
     *                            UnitCoef(2),                   !|  |  |  |
     *                            xabs, TERang, TER_angle)       !|  |  |  |
                 elseif (Mode_Scan.eq.8) Then !-------------------+  |  |  |
                    isico = 0            !skip sin/cos calculat.  |  |  |  |
                    iwarn = 2            !take sin(QB)=1, QB=90   |  |  |  |
                    QB = Bragan (Wave,indices,3,iwarn,d_spacing, !|  |  |  |
     *                           isico,c1,c1,c1,c1,c1,c1)        !|  |  |  |
                 endif !------------------------------------------+  |  |  v
                 Psi = 2.*sin(QB*AngDegree) * sin(Pfi*UnitCoef(2))  !|  |  |
c k0 deviation: k0*(dE/E)                                           !|  |  |
                 dk = k00 * dE/Energy0                              !|  |  |
                 Call  DVecCon(Scan_vector, dk, dk_vector, 3)       !|  |  |
c "Deviated" k0:                                                    !|  |  |
                 k0 = k00 + dk                                      !|  |  |
c Deviations from the Bragg condition:                              !|  |  |
                 Alpha = 2.*DVecSca2(dk_vector,h_vector,3)/k0**2    !|  |  |
                 f0 = Gam0                                          !|  |  |
                 fh2 = (f0 + Psi)**2 - Alpha                        !|  |  |
                 if (fh2.gt.0.D0)  Then !------------------------+   |  |  v
                    fh = Dsqrt (fh2)                            !|   |  |  |
                    if (fh.gt.1.D0) fh=1.0D0  !Ignore fh>1 error!|   |  |  |
                 else  !-----------------------------------------+   |  |  |
                    fh = 0.0D0                                  !|   |  |  |
                 endif  !----------------------------------------+   |  |  v
              endif  !-----------------------------------------------+  |  |
              if (f0.lt.TER_threshold) m_0 = 1   !yes grazing incidence |  |
              if (fh.lt.TER_threshold) m_h = 1   !yes grazing exit      |  |
  225         continue  !---------<------------<------------v          !|  |
           enddo !======================================================+  |
c M7: matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s.                       |
           if (m_0.eq.1 .AND.  m_h.eq.1) Then !---+                        |
              m_type = 0                         !| 4x4                    |
           elseif (m_0.eq.1) Then !---------------+                        |
              m_type = 1                         !| 3x3 graz.incidence     |
           elseif (m_h.eq.1) Then !---------------+                        |
              m_type = 2                         !| 3x3 graz.exit          |
           else !---------------------------------+                        |
              m_type = 3                         !| 2x2 normal diffr.      |
           endif !--------------------------------+                        |
           write (3,68)  MatrixType(m_type+1)                             !|
           write (3,680)                                                  !|
     *                   DAsin(f0)/UnitCoef(2),                           !|
     *                   Uni(iUnits(2)+1),                                !|
     *                   DAsin(fh)/UnitCoef(2),                           !|
     *                   Uni(iUnits(2)+1),                                !|
     *                   TER_threshold/UnitCoef(2),                       !|
     *                   Uni(iUnits(2)+1),                                !|
     *                   TER_Angle,                                       !|
     *                   Uni(iUnits(2)+1),                                !|
     *                   threshold                                        !|
        elseif (m_reduction .eq. 2) Then !---------------------------------+
           m_type = 0                                                     !|
           write (3,68)  'Auto-adjusting at each point'                   !|
        else  !------------------------------------------------------------+
c M7: no Prescan:                                                          |
           m_type = 0                                                     !|
           write (3,68)  'Forced 4x4'                                     !|
        endif !------------------------------------------------------------+
  68    format (' Scattering matrix type: ',a)
  680   format (' Matrix choice criterion:',1x,
     *          ' f0_min=',g14.7,a,', fh_min=',g14.7,a,
     *          ', threshold_angle=',g14.7,a,
     *          ', TER_angle=',g14.7,a,
     *          ', threshold_factor=',f4.1)

c#######################################################
c                  Start of processing:
c#######################################################
        Start  = 0.
        Call Duration (start,bias,hh,mm,ss)
        hh = 0
        mm = 0
        ss = 0
        Call Convert_Time (hh,mm,ss,timestamp)
        nerror     = 0
        nSc        = 0
        progress   = 0.
        Scan_angle = Scan_limits(1)
        if (iBatch_mode .lt. 2)
     *     write (*,49)  nSc, progress,
     *                   Scan_angle,
     *                   Uni(iUnits(3)+1),
     *                   timestamp
  49    format (1x,i4,' points done (',f4.0,'%). Scan angle=',
     *  g11.4,1x,a,' Elapsed time=',a)
        if (modebat .ne. 0) Then !------------------------------+
           DoneFile = OutFile(1:lout)//'.~~~'                  !|
           Call OpenFile(DoneFile(1:lout+4),50,                !|
     *                  'write','unknown',iostts,*89)          !|
           write (50,90,err=91)  nSc, nScan,                   !|
     *                           Scan_angle, Uni(iUnits(3)+1), !|
     *                           progress, timestamp           !|
  91       continue                                            !|
           close(unit=50)                                      !|
        endif  !------------------------------------------------+
  89    continue
  90    format ('Points done = ',i4,' of ',i4,' <br> ',
     *          'Current angle = ',g12.5,1x,a,' <br> ',
     *          'Progress = ',f10.2,'% <br> ',
     *          'Elapsed time (hh:mm:ss) =',a,' <br> ')

c#######################################################
c             Loop over scan_angle or energy
c#######################################################
        Account_Sigma   = .True.
        ierrcount       = 0
        iskipped        = 0
        m_type_previous = -1
        Energy0         = Energy
        Wave0           = Wave
        k00             = k0
        do nSc=1,nScan !===============================================+
           Scan_angle = Scan_limits(1) + dScan_angle*(nSc-1)          !|
                                                                      !|
           Rhs = 0.                                                   !|
           R0s = 0.                                                   !|
           if (i_standing .ne. 0) Then !--+                            |
              do i=1,n_standing !===+     |                            |
                 SWSUM(i) = 0.0    !|     |                            |
              enddo !===============+     |                            |
           endif !------------------------+                            |
                                                                      !|
           if (Mode_Scan.lt.6)   Then  !---------------------------+   |
              Call CalcScanAngles(Scan_angle, scan_label,  !input  |   |
     *                           xabs,                     !input  |   |
     *                           iUnits, UnitCoef,         !input  |   |
     *                           k0_vector, k0, h_vector,  !input  |   |
     *                           b_vector, b_len,          !input  |   |
     *                           c_vector, c_len,          !input  |   |
     *                           Surface_Normal, Psi,      !input  |   |
     *                           Alpha, f0, fh,            !output |   |
     *                           fh_cmplx,                 !output |   |
     *                           neg_angles)               !output |   |
              if (neg_angles .ne. 0) goto 25  !------>-------------+---+---+
           elseif (Mode_Scan.eq.6) Then  !-------------------------+   |   |
c Scan over fh    --  Mode_Scan=6:                                 |   |   |
              f0 = Gam0                                           !|   |   v
              fh = Dsin (Scan_angle*UnitCoef(3))                  !|   |   |
c Nothing to calculate at negative fh:                             |   |   |
              if (fh .le. 0.) goto 25 !-------->-------------------+---+---+
              Alpha = (f0 + Psi)**2 - fh**2                       !|   |   |
           else  !-------------------------------------------------+   |   |
c Energy scan:                                                     |   |   |
              dE = sngl(Scan_angle*UnitCoef(3))                   !|   |   |
              Energy = Energy0 + dE                               !|   |   v
              if (Energy .le. 0) Then !------------------------+   |   |   |
                 write(txt,124) ProgramVer(1:ipv),            !|   |   |   |
     *                          nSc, Energy                   !|   |   |   |
  124            format (a,': at point =',i4,' scan resulted'/
     *           'in negative energy E =',f7.3,'KeV.'//
     *           'Please, correct scan limits!')
                 Call Message (txt,4,2)                       !|   |   |   |
                 goto 28                                      !|   |   |   |
              endif !------------------------------------------+   |   |   |
              Wave   = wave2energy(Energy)                        !|   |   |
              if (Mode_Scan.eq.7) Then !-----------------------+   |   |   |
                 Call recalcXh (Wave, indices,                !|   |   |   v
     *                          N_Total, N_Frac_Max,          !|   |   |   |
     *                          N_Frac,                       !|   |   |   |
     *                          Code, Frac,                   !|   |   |   |
     *                          x0, xh, xf,                   !|   |   |   |
     *                          W0, Wh,                       !|   |   |   |
     *                          C0h, QB, ifail)               !|   |   |   |
                 if (ifail .ne. 0) Then !------------------+   |   |   |   v
                    ifail = abs(ifail)                    !|   |   |   |   |
c                   if (Index(txt(1),                     !|   |   |   |   |
c    *               'does not exist for this wavelength')!|   |   |   |   |
c    *                          .ne.0) Then !---+          |   |   |   |   |
c                      iskipped = iskipped+1   !|          |   |   |   |   |
c                      goto ...                !|          |   |   |   |   |
c                   endif !---------------------+          |   |   |   |   |
                    txt(ifail+1) =                        !|   |   |   |   |
     *                  '-- when recalculating x0,xh at'  !|   |   |   |   |
                    write(txt(ifail+2),125) nSc, Energy   !|   |   |   |   |
  125               format ('scan point =',i4,
     *                      ' energy E =',f7.3,'KeV.')    !|   |   |   |   |
                    ifail = ifail+2                       !|   |   |   |   |
                    Call Message (txt,ifail,2)            !|   |   |   |   |
                    goto 28                               !|   |   |   |   |
                 endif !-----------------------------------+   |   |   |   |
c 2018.10: changed xabs from |x0_substr| to |x0_max|           |   |   |   |
c                xabs   = abs(x0(N_Total+1))                  !|   |   |   |
c                TERang = sqrt(xabs)                          !|   |   |   |
c                if (Abs(xabs) .lt. (1.E-24)) xabs = 1.       !|   |   |   |
                 Call xabsmax (x0, N_Total_Max, N_Total,      !|   |   |   |
     *                         UnitCoef(2),                   !|   |   |   |
     *                         xabs, TERang, TER_angle)       !|   |   |   |
              elseif (Mode_Scan.eq.8) Then !-------------------+   |   |   |
                 isico = 0              !skip sin/cos calculat.|   |   |   |
                 iwarn = 2              !take sin(QB)=1, QB=90 |   |   |   |
                 QB = Bragan (Wave,indices,3,iwarn,d_spacing, !|   |   |   |
     *                        isico,c1,c1,c1,c1,c1,c1)        !|   |   |   |
              endif !------------------------------------------+   |   |   v
              Psi = 2.*sin(QB*AngDegree) * sin(Pfi*UnitCoef(2))   !|   |   |
c k0 deviation: k0*(dE/E)                                         !|   |   |
              dk = k00 * dE/Energy0                               !|   |   |
              Call  DVecCon(Scan_vector, dk,                      !|   |   |
     *                      dk_vector, 3)                         !|   |   |
c "Deviated" k0:                                                  !|   |   |
              k0 = k00 + dk                                       !|   |   |
c Deviations from the Bragg condition:                             |   |   |
              Alpha = 2.*DVecSca2(dk_vector,h_vector,3)/k0**2     !|   |   |
              f0 = Gam0                                           !|   |   |
              fh2 = (f0 + Psi)**2 - Alpha                         !|   |   |
              if (fh2.gt.0.D0)  Then !-----------------------+     |   |   v
                 fh = Dsqrt (fh2)                           !|     |   |   |
                 if (fh.gt.1.D0) Then  !------------------+  |     |   |   |
                    if ((fh-1.0D0).gt.(1.0D-6)) Then !-+  |  |     |   |   |
                       write(3,75) fh,                !|  |  |     |   |   |
     *                        scan_label(1:lelab),    !|  |  |     |   |   |
     *                        Scan_angle,             !|  |  |     |   |   |
     *                        Uni(iUnits(3)+1),       !|  |  |     |   |   |
     *                        Alpha/xabs,             !|  |  |     |   |   |
     *                        DAsin(f0)/UnitCoef(2),  !|  |  |     |   |   |
     *                        Uni(iUnits(2)+1)        !|  |  |     |   |   |
  75                   format(/' WARNING: sin(fh)=',
     *                    g15.8,' > 1',
     *                   ' -- at scan_',a,'=',g14.7,a,
     *                   '    Alpha=',g14.7,'* |x0|',
     *                   '    f0=',g14.7,a)           !|  |  |     |   |   |
                    endif  !---------------------------+  |  |     |   |   |
                    fh = 1.0D0     !Ignore fh>1 error!!!  |  |     |   |   |
                 endif  !---------------------------------+  |     |   |   |
                 fh_cmplx = Dcmplx(fh,0.0d0)                !|     |   |   |
              else  !----------------------------------------+     |   |   |
                 fh = 0.0D0                                 !|     |   |   |
                 fh_cmplx = Dcmplx(fh,Dsqrt(-fh2))          !|     |   |   |
              endif  !---------------------------------------+     |   |   v
           endif  !------------------------------------------------+   |   |
           if (m_reduction .eq. 2) Then !----------------+             |   |
c M7: matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s.     |             |   |
              if     (f0.lt.TER_threshold .AND.         !|             |   |
     *                fh.lt.TER_threshold) Then !--+     |             |   |
                 m_type = 0                       !|     |4x4          |   |
              elseif (f0.lt.TER_threshold) Then !--+     |             |   |
                 m_type = 1                       !|     |3x3 graz.inci|   |
              elseif (fh.lt.TER_threshold) Then !--+     |             |   |
                 m_type = 2                       !|     |3x3 graz.exit|   |
              else !-------------------------------+     |             |   |
                 m_type = 3                       !|     |2x2          |   |
              endif !------------------------------+     |             |   |
              if (m_type .ne. m_type_previous) Then !-+  |             |   |
                 write (3,69) nSc,                   !|  |             |   |
     *                 MatrixType(m_type+1),         !|  |             |   |
     *                 DAsin(f0)/UnitCoef(2),        !|  |             |   |
     *                 Uni(iUnits(2)+1),             !|  |             |   |
     *                 DAsin(fh)/UnitCoef(2),        !|  |             |   |
     *                 Uni(iUnits(2)+1),             !|  |             |   |
     *                 TER_threshold/UnitCoef(2),    !|  |             |   |
     *                 Uni(iUnits(2)+1),             !|  |             |   |
     *                 TER_Angle,                    !|  |             |   |
     *                 Uni(iUnits(2)+1),             !|  |             |   |
     *                 threshold                     !|  |             |   |
  69             format (' point=',i5,5x,
     *              ' Scatt.matrix type set to: ',a,
     *              ', at f0=',g14.7,a,
     *              ', fh=',g14.7,a,
     *              ', threshold_angle=',g14.7,a,
     *              ', TER_angle=',g14.7,a,
     *              ', threshold_factor=',f4.1)      !|  |             |   |
                 m_type_previous = m_type            !|  |             |   |
              endif !---------------------------------+  |             |   |
           endif !---------------------------------------+             |   |
                                                                      !|   |
c##################################                                    |   |
c     Loop over polarizations                                          |   |
c##################################                                    |   |
           Failure = .False.                                          !|   |
           do    jPol=MinPol,MaxPol  !=============================+   |   |
c If polarization-1 did not fail:                                 !|   |   |
              if (.NOT.Failure)   Then   !----------------------+  |   |   |
                                                               !|  |   |   |
                 QBr = QB*AngDegree                            !|  |   |   v
c Modification of 22.12.95                                      |  |   |   |
c (see /X-ray2/RFMm7_RX.FOR )                                   |  |   |   |
                 Rh = RFM8_m7(Alpha)   !R0 is returned in common|  |   |   |
                 if (Mode_Scan .eq. 6) Rh = 2.*fh*Rh           !|  |   |   |
                                                               !|  |   |   |
                 if (ifail .eq. 0) Then !-------------------+   |  |   |   |
                    Rhs = Rhs + Rh*cp(jPol)                !|   |  |   |   |
                    R0s = R0s + R0*cp(jPol)                !|   |  |   |   |
                    if (i_standing .ne. 0) Then !-------+   |   |  |   |   |
                       do i=1,n_standing !============+ |   |   |  |   |   |
                          SWSUM(i) = SWSUM(i)        !| |   |   |  |   |   |
     +                             + SW(i)*cp(jPol)  !| |   |   |  |   |   |
                       enddo !========================+ |   |   |  |   |   |
                    endif !-----------------------------+   |   |  |   |   |
                 else  !------------------------------------+   |  |   |   |
c ifail#0:                                                 !|   |  |   |   |
                    ierrcount = ierrcount + 1              !|   |  |   |   |
                    write(3,77) ierrcount,                 !|   |  |   |   |
     *                          ifail, nSc,                !|   |  |   |   |
     *                          scan_label(1:lelab),       !|   |  |   |   |
     *                          Scan_angle,                !|   |  |   |   |
     *                          Uni(iUnits(3)+1),          !|   |  |   |   v
     *                          Alpha/xabs,                !|   |  |   |   |
     *                          DAsin(f0)/UnitCoef(2),     !|   |  |   |   |
     *                          Uni(iUnits(2)+1),          !|   |  |   |   |
     *                          DAsin(fh)/UnitCoef(2),     !|   |  |   |   |
     *                          Uni(iUnits(2)+1),          !|   |  |   |   |
     *                          jpol                       !|   |  |   |   |
  77                format(i4,'. RFM8_m7: Failure code=',i4,
     *                     '   at step=',i4,:
     *                  '   scan_',a,'=',g14.7,a,
     *                       '   Alpha=',g14.7,'* |x0|',
     *                          '   f0=',g14.7,a,
     *                          '   fh=',g14.7,a,
     *                        '   ipol=',i1)               !|   |  |   |   |
                    if (iDebug .ne. 0) Then !----+          |   |  |   |   |
c Keep writing data:                             |          |   |  |   |   |
                       Rhs = Rhs + Rh*cp(jPol)  !|          |   |  |   |   |
                       R0s = R0s + R0*cp(jPol)  !|          |   |  |   |   |
                    else !-----------------------+          |   |  |   |   |
c Clean everything:                             !|          |   |  |   |   |
                       Rhs = 0.                 !|          |   |  |   |   |
                       R0s = 0.                 !|          |   |  |   |   |
                    endif  !---------------------+          |   |  |   |   |
c Clean SW in either case:                                  |   |  |   |   |
                    if (i_standing.ne.0) Then !-+           |   |  |   |   |
                       do i=1,n_standing !==+   |           |   |  |   |   |
                          SWSUM(i) = 0.    !|   |           |   |  |   |   |
                       enddo !==============+   |           |   |  |   |   |
                    endif !---------------------+           |   |  |   |   |
                    Failure = .True.                       !|   |  |   |   |
                 endif  !-----------------------------------+   |  |   |   |
                 if (iDebug .ne. 0) Then !--------+             |  |   |   |
                    write (33,*) Scan_angle,     !|             |  |   |   |
     *                           Alpha/xabs,     !|             |  |   |   |
     *                           Rh,             !|             |  |   |   |
     *                           R0,             !|             |  |   |   |
     *                           N_Used, jPol    !|             |  |   |   |
                 endif !--------------------------+             |  |   |   |
              endif  !------------------------------------------+  |   |   |
c End of loop over polarizations                                   |   |   |
           enddo   !===============================================+   |   |
                                                                      !|   v
  25       continue                    !<------------<-----------------+---+
                                                                      !|
c#######################                                               |
c    Write the curve                                                   |
c#######################                                               |
           if ((f0 .ge. 0.D0) .AND.                                   !|
     *         (fh .ge. 0.D0 .OR. Mode_Scan .ne. 6)) Then !---------+  |
c 1. A -- scan Angle                                                |  |
c 2. D -- Diffracted intensity  [always included at place-2]        |  |
c 3. S -- Specular (transmitted) intensity                          |  |
c 4. B -- Bragg deviation (alpha-parameter)                         |  |
c 5. I -- Incidence_angle                                           |  |
c 6. E -- Exit angle                                                |  |
c 7. M -- iMaginary part of exit angle                              |  |
c 8  2 -- 2*theta (the sum of incident & exit angles)               |  |
              ttt(1) = Scan_angle                                  !|  |
              ttt(2) = Rhs                                         !|  |
              if (Mode_Scan .eq. 6) Then !---+                      |  |
                 ttt(2) = 2.*fh*ttt(2)      !|                      |  |
              endif  !-----------------------+                      |  |
              ttt(3) = R0s                                         !|  |
              if (ttt(2).lt.1.D-99) ttt(2)=0.                      !|  |
              if (ttt(3).lt.1.D-99) ttt(3)=0.                      !|  |
              ttt(4) = Alpha / xabs                                !|  |
              ttt(5) = DAsin(f0) / UnitCoef(2)                     !|  |
              ttt(6) = DAsin(fh) / UnitCoef(2)                     !|  |
              fhi = Dimag(fh_cmplx)                                !|  |
              if (abs(fhi) .lt. 1.0D0) Then !--------+              |  |
                 ttt(7) = DAsin(fhi) / UnitCoef(2)  !|              |  |
              else !---------------------------------+              |  |
c Workaround for cases when fhi>1 (have happened!)   | 2005/11/29   |  |
                 ttt(7) = (90.*fhi) / UnitCoef(2)   !|              |  |
              endif !--------------------------------+              |  |
              ttt(8) = ttt(5)+ttt(6)                               !|  |
              ttt(9) = 0.    !reserved                             !|  |
              ttt(10)= 0.    !reserved                             !|  |
              ttt(11)= 0.    !reserved                             !|  |
              ttt(12)= 0.    !reserved                             !|  |
              io     = 0                                           !|  |
              do j=1,12 !=====================+                     |  |
                 if (iii(j) .gt. 0) Then !--+ |                     |  |
                    io      = io+1         !| |                     |  |
                    ooo(io) = ttt(iii(j))  !| |                     |  |
                 endif  !-------------------+ |                     |  |
              enddo  !========================+                     |  |
              write (1,333,err=161) (ooo(j),j=1,io)                !|  |
                                                                   !|  |
              if (i_standing .ne. 0) Then !--------------------+    |  |
                 if (n_standing .gt. 1) Then !---------------+ |    |  |
                    if (nScan .gt. 1) Then !---------------+ | |    |  |
                       write (2,'(1001e15.6)',err=161)    !| | |    |  |
     *                       (SWSUM(i),i=1,n_standing)    !| | |    |  |
                    else !---------------------------------+ | |    |  |
                       do i=1,n_standing  !==============+ | | |    |  |
                          z8 = Positions(i_reference+1) !| | | |    |  |
     +                       + standing_range(1)        !| | | |    |  |
     +                       + (i-1)*standing_step      !| | | |    |  |
                          write (2,'(2e15.6)',err=161)  !| | | |    |  |
     *                                     z8, SWSUM(i) !| | | |    |  |
                       enddo !===========================+ | | |    |  |
                    endif !--------------------------------+ | |    |  |
                 else !--------------------------------------+ |    |  |
                    write (2,'(2e15.6)',err=161)            !| |    |  |
     *                             Scan_angle, SWSUM(1)     !| |    |  |
                 endif !-------------------------------------+ |    |  |
              endif !------------------------------------------+    |  |
                                                                   !|  |
c Rounding (neglecting tiny precision-loss errors):                 |  |
              if (R0s.gt.1. .AND. R0s.lt.(1.+eps)) R0s=1.          !|  |
              if (Rhs.gt.1. .AND. Rhs.lt.(1.+eps)) Rhs=1.          !|  |
                                                                   !|  |
c This is a reflection coefficient error:                           |  |
c Use mode-10 for debugging (iDebug=1)!                             |  |
              if (iBatch_mode .ne. 0   .AND.                       !|  |
     *            iBatch_mode .ne. 10)  Then !-------------------+  |  |
                 if ((ifail .eq. 880) .OR. (R0s. lt. 0.) .OR.   !|  |  |
     *               (ifail .eq. 881) .OR. (R0s. gt. 1.) .OR.   !|  |  |
     *               (ifail .eq. 990) .OR. (Rhs. lt. 0.) .OR.   !|  |  |
     *               (ifail .eq. 991) .OR. (Rhs. gt. 1.) .OR.   !|  |  |
     *               (ierrcount .gt. 3 )) Then !--------------+  |  |  |
                    write(txt,122) ProgramVer(1:ipv),        !|  |  |  |
     *                             ifail, R0s, Rhs,
     *                             OutFile(1:lout)//'.tbl'   !|  |  |  |
  122               format(
     *              a,': the job was aborted due to numerical'/
     *              'errors detected while calculating reflection'/
     *              'coefficient(s).'/
     *              'ifail=',i3,'  R0s=',g10.4,'  Rhs=',g10.4/
     *              'This may be caused by precision loss because'/
     *              'of weak reflection (too small xh/x0 in a layer)'/
     *              'or large deviations from the Bragg condition.'/
     *              'Check ',a,' for more details.'//)       !|  |  |  |
                    if (m_reduction.gt.0) Then !-----------+  |  |  |  |
                       txt(10) = 'Please report the '//   !|  |  |  |  |
     *                          'problem to the author!'  !|  |  |  |  |
                       txt(11)= ' '                       !|  |  |  |  |
                    else !---------------------------------+  |  |  |  |
                       txt(10) = 'Please note: setting '//!|  |  |  |  |
     *                          'matrix reduction to'     !|  |  |  |  |
                       txt(11)= '"fly" or "prescan" '//   !|  |  |  |  |
     *                          'may resolve the issue.'  !|  |  |  |  |
                    endif !--------------------------------+  |  |  |  |
                    Call Message (txt,11,2)                  !|  |  |  |
                    goto 28                                  !|  |  |  |
                 endif  !-------------------------------------+  |  |  |
              endif !--------------------------------------------+  |  |
           else  !--------------------------------------------------+  |
              iskipped = iskipped+1                                !|  |
              if (iBatch_mode.lt.2 .AND. Mode_Scan.ne.6) Then !--+  |  |
                 write (3,70) nSc,Scan_angle,f0,fh              !|  |  |
  70             format(' pt=',i5,' Scan angle=',g11.4,
     *                  ' f0=',g11.4,' fh=',g11.4,' - skipping')!|  |  |
              endif  !-------------------------------------------+  |  |
           endif !--------------------------------------------------+  |
  333      format (1x,10g17.8)                                        !|
                                                                      !|
           if (nSc .eq. (nSc/iBeepStep)*iBeepStep) Then !---+          |
              Call  Duration (start,bias,hh,mm,ss)         !|          |
              Call Convert_Time (hh,mm,ss,timestamp)       !|          |
              progress = 100.*Real(nSc)/Real(nScan)        !|          |
              if (progress .gt. 99.99) progress = 99.99    !|          |
              if (iBatch_mode .lt. 2)                      !|          |
     *           write (*,49)  nSc, progress,              !|          |
     *                         Scan_angle,                 !|          |
     *                         Uni(iUnits(3)+1),           !|          |
     *                         timestamp                   !|          |
              if (modebat .ne. 0) Then !------------------+ |          |
                 Call OpenFile(DoneFile(1:lout+4),50,    !| |          |
     *                     'write','unknown',iostts,*93) !| |          |
                 write(50,90,err=92)nSc, nScan,          !| |          |
     *                              Scan_angle,          !| |          |
     *                              Uni(iUnits(3)+1),    !| |          |
     *                              progress, timestamp  !| |          |
  92             continue                                !| |          |
                 close(unit=50)                          !| |          |
              endif  !------------------------------------+ |          |
  93          continue                                     !|          |
           endif  !-----------------------------------------+          |
                                                                      !|
           Call  KeyIn   (0,iscan,iasci)                              !|
c If interrupted by <Esc> and confirmed:                              !|
           if (iasci.eq.27) Then !-----------------------------+       |
                                                              !|       |
              write (3,444) nSc,nScan                         !|       |
              goto 29 !->--------------------------------------+-------+--+
                                                              !|       |  |
           elseif (FileExist(OutFile(1:lout)//'.esc')) Then !--+       |  v
              Call DelFile (OutFile(1:lout)//'.esc',kkk)      !|       |  |
              write (3,444) nSc,nScan                         !|       |  |
              goto 29 !->--------------------------------------+-------+--+
           endif  !--------------------------------------------+       |  |
c End of loop over scan angles                                         |  |
        enddo  !=======================================================+  v
  444   format(/'User interrupted at point',i6,' of',i6)                 !|
                                                                         !|
  29    continue    !<-------------------<------------------------<-------+
        if (nSc .gt. nScan) nSc = nScan
c #######################################################
c                     E X I T
c#######################################################
        close(unit=1,err=131)
  131   continue
        if (i_standing .ne. 0) close(unit=2,err=132)
  132   continue
        Call Duration (start,bias,hh,mm,ss)
        Call Convert_Time (hh,mm,ss,timestamp)
        write   (3,63)  timestamp, ierrcount, iskipped
  63    format  (/' Elapsed time (hh:mm:ss)=',a/
     *            ' Total of errors encountered=',i4/
     *            ' Points skipped because of negative incidence=',i4)
        close(unit=3,err=133)
  133   continue
        if (modebat .ne. 0) Then !---------------------------+
           Call OpenFile(DoneFile(1:lout+4),50,             !|
     *                   'write','unknown',iostts,*95)      !|
           Scan_angle = Scan_limits(1)                      !|
     +                + dScan_angle*(nSc-1)                 !|
           progress = 100.*Real(nSc)/Real(nScan)            !|
           write (50,90,err=94) nSc, nScan,                 !|
     *                         Scan_angle,Uni(iUnits(3)+1), !|
     *                         progress, timestamp          !|
  94       continue                                         !|
           close(unit=50)                                   !|
        endif  !---------------------------------------------+
  95    continue
        if (iskipped .eq. nSc) Then  !--------------------------+
           write (txt,121) ProgramVer(1:ipv)                   !|
  121      format(a,':'//
     *    'No output data written because all the points in'/
     *    'the scan correspond to negative incidence angles')  !|
           Call Message (txt,4,2)                              !|
        endif  !------------------------------------------------+
        if (nSc .lt. nScan)     goto 28

c       txt(1)=' *** Work done. Press any key to continue..'
        txt(1)=' *** Work done...'
        if (iBatch_mode .lt. 2) Then !-----+
           l = Len_Trim (txt(1))          !|
           write (*,'(1x,a)') txt(1)(1:l) !|
        endif  !---------------------------+
        do      i=1,1   !===================+
           Call  KeyIn   (0,iscan,iasci)   !|
           if (iscan .ne. 0)     goto 28   !|
        enddo  !============================+
c       if (iBatch_Mode .eq. 0)  Call KeyIn (1,iscan,iasci)
        goto 28

c #######################################################
c                  Error messages
c#######################################################
  161   continue
        write   (txt,120)  ProgramVer(1:ipv),
     *                     OutFile(1:lout)//'.dat'
  120   format  (a,': Cannot write output file.'//a//
     *  'May be, disk is full or drive is not ready.')
        Call Message (txt,5,2)
        close (unit=1,status='delete',err=28)
        if (i_standing .ne. 0)
     *  close (unit=2,status='delete',err=28)
        goto 28

  162   continue
        write   (txt,123)  ProgramVer(1:ipv),
     *                     l, Sigma(l), Exti
  123   format  (a,': in layer=',i5/
     *  'rms roughness Sigma =',g12.5,' A'/
     *  'exceeds'/
     *  'Extinction legnth =',g12.5,' A'/
     *  'This is beyond the assumptions of roughness model.'/
     *  'Please use transition layers to simulate roughness.')
        Call Message (txt,6,2)
        close (unit=1,status='delete',err=28)
        if (i_standing .ne. 0)
     *  close (unit=2,status='delete',err=28)
        goto 28

  28    continue
        Call exit_quiet()
        End

c===================================================================

        Subroutine CalcScanAngles (Scan_angle, scan_label,      !input
     *                             xabs,                        !input
     *                             iUnits, UnitCoef,            !input
     *                             k0_vector, k0, h_vector,     !input
     *                             b_vector, b_len,             !input
     *                             c_vector, c_len,             !input
     *                             Surface_Normal, Psi,         !input
     *                             Alpha, f0, fh,               !-output
     *                             fh_cmplx,                    !-output
     *                             neg_angles)                  !-output

        Complex*16      fh_cmplx
        Real*8          f0,                     !sin of incident angle
     *                  fh,                     !sin of exit angle
     *                  Alpha,                  !deviation from Bragg
     *                  Scan_angle,             !scan angle
     *                  h_vector(3),
     *                  Surface_Normal(3),
     *                  k0_vector(3), k0,
     *                  dk_vector(3),
     *                  wrk_vector(3), k0_,
     *                  b_vector(3), b_len,
     *                  c_vector(3), c_len,
     *                  x8, y8, fh2
        Real*4          Psi,                    !misorientation angle
     *                  xabs,                   !abs(x0_max)
     *                  UnitCoef(3)
        Integer         iUnits(3),
     *                  neg_angles,
     *                  lelab
        Character       scan_label*(*),
     *                  Uni(5)*5  /'degr.','min.',
     *                             'mrad.','sec.', 'urad.'/

        Real*8          DVecSca2, DVecMod2
        External        DVecSca2, DVecMod2

        Alpha      = 0.
        f0         = 0.
        fh         = 0.
        fh_cmplx   = (0.,0.)
        neg_angles = 0
        lelab      = Max(Len_Trim(scan_label),1)

c+---------------------------------+
c|Vector deviation from the Bragg  |
c|condition on the INCIDENCE side: |
c+---------------------------------+
c ->    ->    ->
c dk = x*c + y*b,
c  x = (b/c)*sin(dQ)
c  y = -2*sin^2(dQ/2)
        x8 = (b_len/c_len)*Dsin(Scan_angle*UnitCoef(3))
        y8 =-2.*(Dsin(Scan_angle*UnitCoef(3)/2.))**2
        Call DVecSum2 (c_vector, x8,
     +                 b_vector, y8,
     =                 dk_vector, 3)

c Deviations from the Bragg condition:
        Alpha = 2.*DVecSca2(dk_vector,h_vector,3)/k0**2

c "Deviated" k0_vector:
        Call  DVecSum2 (k0_vector, 1.0D0,
     +                  dk_vector, 1.0D0,
     =                  wrk_vector, 3)

c This must be k0, but for any error case we re-normalize it:
        k0_ = DVecMod2 (wrk_vector,3)
        if (Abs(1.-k0_/k0) .gt. (1.0D-6)) Then !-----+
c Report Warning:                                    |
           write(3,1) k0, k0_,                      !|
     *                scan_label(1:lelab),          !|
     *                Scan_angle,                   !|
     *                Uni(iUnits(3)+1),             !|
     *                Alpha/xabs                    !|
  1        format(/' WARNING:  K0 length mismatch:'/
     *     ' Expected=',g15.8,' Obtained=',g15.8/
     *     ' -- at scan_',a,'=',g14.7,a,
     *             '    Alpha=',g14.7,'* |x0|')     !|
        endif  !-------------------------------------+

c "deviated" Gam0:
        f0 = DVecSca2 (wrk_vector,Surface_Normal,3)/k0_
c Workaround for f0>1 (small round errors, have happened!):
        if (f0. gt. 1.0D0) Then !-----------------------------+
           if (f0 .gt. 1.01D0) Then !----------------------+  |
              write(3,2) f0,                              !|  |
     *                   scan_label(1:lelab),             !|  |
     *                   Scan_angle,                      !|  |
     *                   Uni(iUnits(3)+1)                 !|  |
  2           format(/' WARNING: sin(f0)=',g15.8,' > 1'/
     *                ' -- at scan_',a,'=',g14.7,a)       !|  |
           endif !-----------------------------------------+  |
           f0 = 1.0D0                                        !|
        endif !-----------------------------------------------+

        if (f0 .le. 0.D0) Then !--+
           neg_angles = 1        !|
           return                !|
        endif !-------------------+

        fh2 = (f0 + Psi)**2 - Alpha
        if (fh2 .gt. 0.D0) Then !----------------------------------+
           fh = Dsqrt (fh2)                                       !|
           if (fh .gt. 1.D0) Then  !----------------------------+  |
              if ((fh-1.0D0) .gt. (1.0D-6)) Then !-----------+  |  |
                 write(3,3) fh,                             !|  |  |
     *                      scan_label(1:lelab),            !|  |  |
     *                      Scan_angle,                     !|  |  |
     *                      Uni(iUnits(3)+1),               !|  |  |
     *                      Alpha/xabs,                     !|  |  |
     *                      DAsin(f0)/UnitCoef(2),          !|  |  |
     *                      Uni(iUnits(2)+1)                !|  |  |
  3              format(/' WARNING: sin(fh)=',g15.8,' > 1'/
     *                   ' -- at scan_',a,'=',g14.7,a,
     *                   '    Alpha=',g14.7,'* |x0|',
     *                   '    f0=',g14.7,a)                 !|  |  |
              endif  !---------------------------------------+  |  |
              fh = 1.0D0  !Ignore fh>1 error!!!                 |  |
           endif  !---------------------------------------------+  |
           fh_cmplx = Dcmplx(fh,0.0d0)                            !|
        else  !----------------------------------------------------+
           fh = 0.0D0                                             !|
           fh_cmplx = Dcmplx(fh,Dsqrt(-fh2))                      !|
c If fh=0 (or becomes imaginary), the                              |
c diffracted intensity becomes zero,                               |
c but the reflected is not. Therefore,                             |
c the following flag is commented out                              |
c and we continue:                                                 |
c cccc     neg_angles = 1                                         !|
        endif  !---------------------------------------------------+
        return
        end
