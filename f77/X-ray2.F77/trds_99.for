        Program trds_99
c       USE MSFLIB                                      !for QuickWin
c+=================================================================+
c|       X-ray diffuse scattering from interface roughness         |
c|   under grazing-incidence specular reflection in multilayers    |
c|           *** Version'99 optimized for WWW access ***           |
c|-----------------------------------------------------------------|
c|                       By S.A.Stepanov                           |
c|-----------------------------------------------------------------|
c|'97  Script language for top-layer profile, allowing specifying  |
c|     periodic structures; improved batch ("silent" or "server")  |
c|     mode for using with CGI-interfaces.                         |
c|'99  Same as '97 but:                                            |
c|     1. added x-ray polarization dependence                      |
c|     2. added possibility to specify rho and chemical formulas   |
c|     3. added "Identical" flag to accelerate calculation of x0.  |
c|'99.2 Added Henke/Cowan option                                   |
c|'99.3 Added user comment and startup check for zero roughness    |
c+=================================================================+

c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter      (N_Top_Max     = 150,
c cc *                  N_Total_Max   = N_Top_Max+1)

        Integer         N_Points_Max
        Parameter      (N_Points_Max  = 1001)   ! 1001

        Integer         N_Frac_Max
        Parameter      (N_Frac_Max    = 4)

        Real*8          DSSL(N_Points_Max),     ! array for intensity
     *                  xi2_n, p2_n, S

        Real*8          DSFM                    ! (function)
        External        DSFM

        Real*4          Frac(N_Frac_Max,N_Total_Max),
     *                  rho(N_Total_Max),
     *                  W0(N_Total_Max),
     *                  Rc(2),                  ! refl.coeff.-s
     *                  Theta(2),               ! incid./exit angles
     *                  k4,                     ! (real*4)k0
c cc *                  UnitCoef(3),            ! radians/(angular unit)
     *                  Energy, Full_Points,
     *                  Thickness_Top,
     *                  TER_Angle, Sigma_mean,
     *                  Sigma_max,
     *                  Scan_Range(2,2),
     *                  Scan_Step(2),
     *                  Theta1, Theta2, Theta_Cnt,
     *                  Theta_Smp, Theta_Off,
     *                  x_, y_, z_, r, t, start,
     *                  qx_eps, qz_eps, roo, soo

        Real            wave2energy
        External        wave2energy

        Integer*4       hh, mm, ss, bias,
     *                  nfa,                            ! number of fails
     *                  iii, jjj,
     *                  idat, igrd, jgrd

        Integer         N_Frac(N_Total_Max),
     *                  Identical(N_Total_Max),
     *                  N_Used(2), iUnits(3),
     *                  nScan(2), mode_scan,
     *                  iBeepStep,
     *                  i, j, l, m, n,  ipol,           ! 1=sigma, 2=pi
     *                  ipv,
     *                  iskip, icont, jcont,
     *                  lout, ltop, line, kkk,
     *                  iBatch_mode,
     *                  iskipped, looped, lug,
     *                  nerror, iscan, iasci,
     *                  io_status

        Logical         Ok_angles, Ok_symmetry,
     *                  GRD_Output, Intrp

        Character       Code(N_Frac_Max,N_Total_Max)*20,
     *                  C00(N_Total_Max)*2,
     *                  InpFile*80, OutFile*80,
     *                  File_Top*80, DoneFile*80,
     *                  Radiat*6, timestamp*9,
     *                  engday*12, txdate*10, txtime*5,
     *                  ProgramVer*8, Version(2)*80, buff*5120,
     *                  Uni(5)*5 /
     *                           'degr.',                               !1
     *                           'min.',                                !2
     *                           'mrad.',                               !3
     *                           'sec.',                                !4
     *                           'urad'/,                               !5
     *                  Spm(2)*5 /
     *                           'Sigma',                               !1
     *                           ' Pi  '/,                              !2
     *                  SpecMode(2)*18 /
     *                           'Diffuse scattering',                  !1
     *                           'Reflection coeff.'/,                  !2
     *                  Scan_Type(4)*30 /
     *                           'Q-scans at 2Q=fixed',                 !1
     *                           'Q-2Q scans at fixed Q-offsets',       !2
     *                           '2Q (detector) scans at Q=fixed',      !3
     *                           'qx scans at qz=fixed'/,               !4
     *                  Corfu(20)*31 /
     *                           'NO vertical correlation        ',     !01
     *                           'Complete vertical replication  ',     !02
     *                           'Ming''s model                  ',     !03
     *                           'Lagally''s model               ',     !04
     *                           'Holy''s model    (accumulated) ',     !05
     *                           'Spiller''s model (accumulated) ',     !06
     *                           'Classic Pukite''s model (steps)',     !07
     *                           'Smoothed Pukite''s model(steps)',     !08
     *                           'Pershan''s model        (steps)',     !09
     *                           'Ming + lateral periodicity     ',     !10
     *                           'Not implemented yet..          ',     !11
     *                           'Not implemented yet...         ',     !12
     *                           'Not implemented yet...         ',     !13
     *                           'Not implemented yet...         ',     !14
     *                           'Not implemented yet...         ',     !15
     *                           'Not implemented yet...         ',     !16
     *                           'Not implemented yet...         ',     !17
     *                           'Not implemented yet...         ',     !18
     *                           'Not implemented yet...         ',     !19
     *                           'Not implemented yet...         '/,    !20
     *                  Text_Born(3)*22 /
     *                           'DWBA: distorted wave..',              !0
     *                           'semi-Born approx.: R=0',              !1
     *                           'Born approx.: R=0, T=1'/,             !2
     *                  Text_Skip(2)*34 /
     *                           'Write all points, including Q<0',     !0
     *                           'Skip Q<0 -- Continuation disabled.'/, !1
     *                  Text_Fourier(2)*11 /
     *                           'exp[K(x)]-1',                         !0
     *                           'Fourier(K)'/,                         !1
     *                  Text_Summar(2)*5 /
     *                           '     ',                               !0
     *                           '+Ming'/                               !1

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

        Integer         n_Column, LefTxt
        External        n_Column, LefTxt

        Logical*4       FileExist
        External        FileExist

c--------------------------------------------------------
        Complex*8       x0(N_Total_Max+1)
        Real*8          Z(N_Total_Max),         ! coordinates of interfaces
     *                  k0                      ! k0 = 2*pi/wave
        Real*4          Wave, UnitCoef(3), qz, qx,
     *                  Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Sigma_Accu(N_Total_Max),
     *                  CorreLength,CorreVert,h_jagged,
     *                  xabs, TERang, pi
        Integer         i_Spec, N_Top, N_Total, Index_corfu,
     *                  MaxOrder(2), iDebug, ifail
        Common  /TRDSdat/ Wave, UnitCoef, qz, qx, x0, Z,
     *                    Thickness, Sigma, Sigma_Accu,
     *                    CorreLength, CorreVert, h_jagged,
     *                    xabs, TERang, pi, k0, i_Spec,
     *                    N_Top, N_Total, Index_corfu,
     *                    MaxOrder, iDebug, ifail

        Real*4          CorreCross, RelaxLength,
     *                  EpsRel, RouMiscut,
     *                  SurMiscut, TerraceSpre,
     *                  StepHeight, Reserved1,
     *                  Reserved2, Reserved3,
     *                  S_factor, A_Factor,
     *                  Relief_Period
        Integer         i_Born, i_Fourier,
     *                  isummariz, max_Periods
        Common  /ADDsigma/CorreCross, ! Lateral cor.leng.i#j (Lagally)
     *                    RelaxLength,! Diffusion-like relax.length (Spiller)
     *                    EpsRel,     ! integration precision(Spiller)
     *                    RouMiscut,  ! relief transfer angle
     *                    SurMiscut,  ! Steps (Pukite): surface miscut angle
     *                    TerraceSpre,! Steps (Pukite): terrace spread
     *                    StepHeight, ! Steps (Pukite): step height
     *                    Reserved1,  ! reserved roughness parameter-1
     *                    Reserved2,  ! reserved roughness parameter-2
     *                    Reserved3,  ! reserved roughness parameter-3
     *                    i_Born,     ! approximation for wavefields (0=DWBA,1=semi-Born,2=Born)
     *                    i_Fourier,  ! corr.integral: exp[K(x)]-1 or K(x)
     *                    isummariz,  ! add roughness-to-steps flag
     *                    Relief_Period,!surface relief period
     *                    max_Periods,! number of correlated periods
     *                    S_factor,   !==(Relief_Period/CorreLength)^2
     *                    A_Factor    ! normalization for periodic func.

        Character       Comment(3)*80
        Common /Cmnt/   Comment

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
        progname = 'trds_99'
c       i = SetExitQQ(2)        !QWIN$EXITNOPERSIST=2   !for QuickWin
c #######################################################
        pi     = 4.* atan(1.)

        ProgramVer = 'trds_99'
        ipv = Len_Trim (ProgramVer)
c This is the startup ERR-file (it is used until we've
c verified the output filename):
        ErrorFile = ProgramVer(1:ipv)//'.err'
        Call    DeleteErrorFile ()

        write   (Version,111)
  111   format  (
     *  'TRDS: total-reflection x-ray diffuse scattering from ',
     *  'surface/interface roughness'/
     *  'By S.Stepanov <sstepanov@anl.gov>',24x,
     *  'Version-99.3: Aug, 2012')

        modebat = 1                     !assume batch-mode by default
        iBatch_mode = modebat

c #######################################################
c              Read input file
c#######################################################

        Call    Read_TRDS_Input (InpFile, OutFile, File_Top,
     *                          ProgramVer, ipv, N_Points_Max,
     *                          iUnits, iBeepStep, mode_scan,
     *                          iBatch_mode, GRD_output, iskip,
     *                          icont, jcont, Thickness_Top,
     *                          Scan_Range, Scan_Step, nScan,
     *                          Code, Frac, N_Frac, N_Frac_Max,
     *                          rho, W0, C00, Identical,
     *                          ipol, Radiat)

        Sigma_mean = 0.
        Sigma_max  = 0.
        if (N_Top.gt.0) Then  !-----------------------------------+
          do    l = 1, N_Top  !=================================+ |
            Sigma_mean = Sigma_mean + Sigma(l)                 !| |
            if (Sigma(l) .gt. Sigma_max) Sigma_max = Sigma(l)  !| |
          enddo  !==============================================+ |
          Sigma_mean = Sigma_mean / N_Top                        !|
        endif  !--------------------------------------------------+
        if (Sigma(N_Total) .gt. Sigma_max) Sigma_max = Sigma(N_Total)
        if (Sigma_max .le. 0.) Then !-------------------------------+
          write (txt,122)  ProgramVer(1:ipv)                       !|
  122     format(a,':'//
     *    'Nothing to do because all layers have zero roughness.') !|
          Call  Message (txt,3,2)                                  !|
        endif  !----------------------------------------------------+

        Energy = wave2energy(Wave)
c #######################################################
        lout = Max(Len_Trim(OutFile),1)
        i = Min(lout,40)
        ltop = Max(Len_Trim(File_Top),1)
        j = Min(ltop,40)

        write (txt,33) OutFile(1:i),
     *                 File_Top(1:j),
     *                 Thickness_Top, N_Top, Sigma_mean,
     *                 Code(1,N_Total), Sigma(N_Total),
     *                 Wave,  Energy,   Radiat,
     *                 Spm(iPol),
     *                 Scan_Type(mode_scan),
     *                 (Uni(iUnits(1)+1), (Scan_Range(i,j),i=1,2),
     *                 nscan(j),j=1,2),
     *                 Text_Born(i_Born+1), Text_Fourier(i_Fourier+1),
     *                 Index_corfu, Corfu(Index_corfu),
     *                                      Text_Summar(isummariz+1),
     *                 CorreLength, CorreVert, h_jagged,
     *                 Uni(iUnits(2)+1), RouMiscut,
     *                 Uni(iUnits(3)+1), SurMiscut,
     *                 StepHeight,
     *                 TerraceSpre
  33    format  (
     *  'Output files [.dat & .tbl] will be...|',a/
     *  'Top: Profile name....................|',a/
     *  'Top: Thickness(A), N_sublayers, <rms>|',f14.2,i6,f7.2/
     *  'Substrate:  Code, rms_roughness.(A)..|',a,    f7.2/
     *  'X-ray wavelength (A), energy (keV)...|',g12.7,1x,g12.7,1x,a/
     *  'X-ray polarization...................|',a/
     *  'Scan type............................|',a/
     *  '[1]: Theta(',a,') or q(1/A), N-points| ',2g15.4,i7/
     *  '[2]: Theta(',a,') or q(1/A), N-points| ',2g15.4,i7/
     *  'Method: DWBA or Born; [exp(K)-1] or K|',a,',',3x,a/
     *  'Roughness: correlation function type.|',i2,'.',2a/
     *  ' - correl.lengths(L_h,L_v),jaggedness|',2(1x,f10.1),9x,f5.3/
     *  ' - skew transfer angle (',a,').......|',g11.3/
     *  ' - surface miscut (',a,')  [Pukite]..|',g11.3/
     *  ' - atomic steps height (A) [Pukite]..|',g11.3/
     *  ' - terrace width spread (A) [Pershan]|',g11.3)

        if (CorreVert.ge.1.E+10)  txt(12)(50:60) = '   Infinity'

        line = 16
        if (Index_corfu.lt.7 .Or.
     *      Index_corfu.gt.9) line=line-3

c #######################################################
c                 Open output files:
c#######################################################
        DoneFile    = OutFile(1:lout)//'.~~~'
        Full_Points = (nscan(1)+0.)*(nscan(2)+0.)
        if (Full_Points.le.0.)  Full_Points=1.
        iskipped    = 0
        looped      = 0
        Intrp       = .False.

c  +===========+   +===========+
c  |DAT output:| / |GRD output:|
c  +===========+   +===========+

        if (icont.eq.0)       Then  !-----------------------+
          Call  Openlist (OutFile(1:lout)//'.dat',1,ifail) !|New start!
          if (ifail.ne.0)     goto 160                     !|
          n = 0                                            !|
        else  !---------------------------------------------+
c Append:                                                  !|
          Call OpenFile(OutFile(1:lout)//'.dat',1,         !|Continuation!
     *              'readwrite','unknown',io_status,*160)  !|
          icont = 0                                        !|
          jcont = 0                                        !|
  5       continue      !<========================+         |
          read (1,888,err=162,end=6) x_,y_,z_ !+  |         |
          icont = icont + 1                   !V  |         |
          if (icont.ge.nscan(1)) Then  !--+    |  |         |
           icont = icont-nscan(1)        !|    |  |         |
           jcont = jcont+1               !|    |  |         |
          endif  !------------------------+    |  |         |
          goto 5  !==============>=============+==+         |
  6       continue       !<-----------End------+            |
          if (jcont.gt.nscan(2))        goto 201           !|Continuation mode error
          if (jcont.eq.nscan(2) .AND.                      !|
     *        icont.gt.0)               goto 201           !|Continuation mode error
          if (jcont.eq.nscan(2) .AND.                      !|
     *        icont.eq.0)               Then  !-+           |
            i  = nscan(1)                      !|           |
            j  = nscan(2)                      !|           |
            r  = 100.                          !|           |
            hh = 0                             !|           |
            mm = 0                             !|           |
            ss = 0                             !|           |
            goto 505                           !|           |nothing to do!
          endif  !------------------------------+           |
          idat = icont + jcont*nscan(1)                    !|
          if (idat.gt.0)  Then  !-----+                     |
            Backspace (1,err=163)    !|                     |
          else  !---------------------|                     |
            Rewind (1,err=164)       !|                     |
  164       continue                 !|                     |
            Rewind (3,err=163)       !|                     |no continuation!
          endif  !--------------------+                     |
  163     continue                                         !|
        endif  !--------------------------------------------+

        if (GRD_output) Then  !--------------------------+
c  +===========+                                        !|
c  |GRD output:|                                        !|
c  +===========+                                        !|
          if (n.eq.0)       Then  !-------------------+  |
            Call  Openlist (OutFile(1:lout)//'.grd', !|  |New start!
     *                                     11,ifail) !|  |
            if (ifail.ne.0)     goto 190             !|  |
          else  !-------------------------------------|  |
c Append:                                            !|  |
            Call OpenFile(OutFile(1:lout)//'.grd',11,!|  |Continuation!
     *          'readwrite','unknown',io_status,*190)!|  |
            igrd = 0                                 !|  |
            jgrd = 0                                 !|  |
            do  j=1,5+jcont  !==================+     |  |
              read (11,22,err=192,end=193) !----+-----+--+-+Data mismatch
  22          format(a)                        !|     |  | V in continue!
              if (j.gt.5) jgrd = jgrd+1        !|     |  |
            enddo  !============================+     |  |
c Uncomplete line:                                    |  |
            if (icont.gt.0)     Then  !------------+  |  |
              read (11,22,err=192,end=193) buff   !|  |  |
              i    = n_Column (buff)              !|  |  |
              igrd = jgrd*nscan(1)+i              !|  |  |
              if (i.eq.nscan(1)) Then !--+         |  |  |
                i    = 0                !|         |  |  |
                jgrd = jgrd+1           !|         |  |  |
                goto 193  !--------------+---------+--+--+-+
              endif  !-------------------+         |  |  | |
              if (i.ne.icont)   goto 193  !--------+--+--+-|Data mismatch
              if (igrd.ne.idat) goto 193  !--------+--+--+-| in continue!
              Call Rdreal8 (DSSL,icont,buff,j)    !|  |  | V
              Backspace (11,err=192)              !|  |  |
            endif !--------------------------------+  |  |
          endif  !------------------------------------+  |
        endif  !-----------------------------------------+

        if (icont.eq.0 .AND. jcont.eq.0) Then  !--------------------+
          write (3,55)  (Version(i)(1:Len_Trim(Version(i))),i=1,2) !|
          Call    DateStam (txdate,engday)                         !|
          Call    TimeStam (txtime)                                !|
          write (3,88) engday(1:3), txdate, txtime                 !|
  88      format(a,', ',a,', ',a)                                  !|
          write (3,11)                                             !|
          write (3,55)  (txt(i)(1:Len_Trim(txt(i))),i=1,line)      !|
          write (3,11)                                             !|
        endif  !----------------------------------------------------+
  11    format  (1x,78('='))
  55    format  (20(a:/))            !1x,a:

        if (iDebug.ne.0)
     *  Open    (unit=33,
     *           file=OutFile(1:lout)//'.deb',
     *           status='unknown')
c #######################################################
c               Print structure profiles:
c#######################################################
        k4 = 2.*pi /wave
        k0 = Dble(k4)
c 2018.10: changed xabs from |x0_substr| to |x0_max|
c       xabs      = Abs(x0(N_Top+2))
c       TERang    = Sqrt(xabs)
c       TER_Angle = TERang / UnitCoef(1)
        Call xabsmax (x0, N_Total_Max, N_Total, UnitCoef(1),
     *                xabs, TERang, TER_angle)

        if (icont.eq.0 .and. jcont.eq.0) Then  !----------+
          if (Index_corfu.eq.1 .or.                      !|non-correlated
     *        Index_corfu.eq.5 .or.                      !|Holy
     *        Index_corfu.eq.6) Then  !------------+      |Spiller
            write (3,62)  RelaxLength, EpsRel     !|      |
          endif  !---------------------------------+      |
          if (Index_corfu.eq.4) Then  !------------+      |Lagally
            write (3,67)  CorreCross              !|      |
          endif  !---------------------------------+      |
          if (Index_corfu.eq.10) Then  !-----------+      |laterally
            write (3,167)  Relief_Period,         !|      |periodic
     *                     max_Periods            !|      |
          endif  !---------------------------------+      |
          Call  Algorithm_Type  (Version(1))             !|
          l = Len_Trim (Version(1))                      !|
          i = Len_Trim (DBtext(iHenkeCowan+2))           !|
          write (3,68)  Text_Skip(iskip+1),              !|
     *                  SpecMode(i_Spec+1),              !|
     *                  Uni(iUnits(1)+1), TER_Angle,     !|
     *                  Version(1)(1:l),                 !|
     *                  DBtext(iHenkeCowan+2)(1:i)       !|
                                                         !|
c Print user comment:                                    !|
          write (3,11)                                   !|
          write (3,55)  ' *** User Comment: '            !|
          do i=1,3  !==================================+  |
            l = Len_Trim (Comment(i))                 !|  |
            if (l.gt.0) write (3,55) Comment(i)(1:l)  !|  |
          enddo  !=====================================+  |
          write (3,11)                                   !|
                                                         !|
c Print surface layer profile:                           !|
          write (3,61)                                   !|
          do      l=1,N_Total   !=======================+ |
            write (3,64) l, Thickness(l),              !| |
     *                   Sigma(l),                     !| |
     *                   Code(1,l), Frac(1,l), rho(l), !| |
     *                   W0(l), x0(l+1), Identical(l)  !| |
            do j=2,N_Frac(l)  !===================+     | |
              write (3,65) Code(j,l), Frac(j,l)  !|     | |
            enddo  !==============================+     | |
          enddo  !======================================+ |
        endif  !------------------------------------------+

  61    format  (/'The Structure Profile:'/
     *            '====================='/
     *  '  Nr Thickness Roughness ________Code________ ',
     *  'Fraction    rho      Exp(-w)',7x,'x0       Identical')

  64    format  (i4,f10.2,f8.2,3x,a,f7.3,f10.4,f11.3,
     *                                       ' (',g8.2,',',g8.2,')',i6)
  65    format  (25x,             a,f7.3)
  62    format  (
     *  'Roughness relaxation length (Spiller).[A]|',g10.3/
     *  'Roughness integral precision (Spiller)...|',g10.3)
  67    format  (
     *  'Cross-correlation length (Lagally)....[A]|',g10.3)
  167   format  (
     *  'Period of surface relief..............[A]|',g10.3/
     *  'Number of correlated periods.............|',i5)
  68    format  (
     *  'Processing mode for points with Theta<0).|',a/
     *  'Intensity computed in specular points....|',a/
     *  'Critical TER angle in ',a,'..............|',g10.3/
     *  'Algorithm type (Matrix / Recursive)......|',a/
     *  'Database used for dispersion correct.....|',a)

c #######################################################
c                Roughness parameters
c#######################################################
        Z(1) = 0.0D0
        do      i=2,N_Total  !==============+
          Z(i) = Z(i-1) + Thickness(i-1)   !|
        enddo  !============================+

c If accumulated roughness:

        xi2_n    = CorreLength**2
        do    i=1,N_Total   !==========================+
          if (                                        !|
c ccc*        Index_corfu.eq.1  .OR.    !Not correlated|
     *        Index_corfu.eq.5  .OR.    !Holy          |
     *        Index_corfu.eq.6)         !Spiller       |
     *                          Then  !--------------+ |
            S = 0.                                  !| |
            do  n=i,N_Total  !====================+  | |
              p2_n = 8.*RelaxLength*(Z(n)-Z(i))  !|  | |
              S = S + Sigma(n)**2                !|  | |
     /              / (1.+p2_n/xi2_n)            !|  | |
            enddo  !==============================+  | |
            Sigma_Accu(i) = SNGL(Sqrt(S))           !| |
          else   !-----------------------------------| |
c ccc         Index_corfu.eq.1  .OR.    !Not correlated|
c ccc         Index_corfu.eq.2  .OR.    !Full        | |
c ccc         Index_corfu.eq.3  .OR.    !Ming        | |
c ccc         Index_corfu.eq.4  .OR.    !Lagally     | |
c ccc         Index_corfu.ge.7          !Pukite & >..| |
            Sigma_Accu(i) = Sigma(i)                !| |
          endif  !-----------------------------------+ |
        enddo  !=======================================+
        if (icont.eq.0 .and. jcont.eq.0)  Then  !----------+
          if (                                            !|
c ccc*        Index_corfu.eq.1  .OR.    !Not correlated   !|
     *        Index_corfu.eq.5  .OR.    !Holy             !|
     *        Index_corfu.eq.6)         !Spiller          !|
     *                          Then  !-------------------+|
            write (3,15)                                 !||
            do  i=1,N_Total   !==========================+||
              write (3,16) i,Thickness(i),Sigma_Accu(i) !|||
            enddo  !=====================================+||
          endif  !----------------------------------------+|
        endif  !-------------------------------------------+
  15    format  (/' Accumulated rms roughness profile:'/
     *            ' (from the top to the substrate)'/
     *            '    Nr      Thickness',5x,
     *            'Rms roughness height (Angstr.)')
  16    format  (i5,4x,f9.1,15x,f8.2)
c#######################################################
c Tabulate the terms of expansion of the correlation function.
c The data for 4001 points in the argument range [0.:100.] are
c stored in the array ArrCos in the common block /ArrcoDat/:
        lug = 88
c+============================================+
        Call    Read_CoFu_TBL (h_jagged,lug) !|
c+============================================+
c #######################################################
c                  Start of processing:
c#######################################################
        start  = 0.
        Call Duration (start,bias,hh,mm,ss)
        nerror = 0
        if (icont.ne.0 .or. jcont.ne.0) Then  !----+
           write (3,444)    icont+1, jcont+1      !|
           if (iBatch_mode .lt. 2)                !|
     *        write (*,444) icont+1, jcont+1      !|
        endif  !-----------------------------------+
  444   format  (1x,'...Continuation started at point ('
     *           i4,',',i4,')...')

        i = icont
        j = max0(jcont,1)
        r = 100.*(i+nscan(1)*(j-1.))/Full_Points
        if (r .gt. 99.99) r = 99.99
        hh = 0
        mm = 0
        ss = 0
        Call Convert_Time (hh,mm,ss,timestamp)

        if (iBatch_mode .lt. 2) write (*,333) i, j, r, timestamp
  333   format  (1x,'(',i4,',',i4,') points done (',
     *  f6.2,'%)',6x,'Elapsed time (hh:mm:ss)=',a)

        if (modebat.ne.0)       Then  !----------------+
          Call OpenFile(DoneFile(1:lout+4),50,        !|
     *              'write','unknown',io_status,*89)  !|
          write (50,90,err=91)   i, j, nscan,         !|
     *                           r, timestamp         !|
  91      continue                                    !|
          close(unit=50)                              !|
        endif  !---------------------------------------+
  89    continue
  90    format (
     *  'Points done = (',i4,',',i4,') of (',i4,',',i4,') <br> ',
     *  'Progress = ',f6.2,'% <br> ',
     *  'Elapsed time (hh:mm:ss) =',a,' <br> ')

c Calculation of dead bands for
c the symmetry accelerator:
        qx_eps = 0.
        qz_eps = 0.
        goto (601,610,610,602) mode_scan  !----------+
c              1   2   3   4                        !|
c+----- Mode1: over Q  at 2Q=fix -----|              V
  601   continue  !<---------------------------------|
        t = Amax1(Abs(Scan_Range(1,1)),             !|
     *            Abs(Scan_Range(1,2)),             !|
     *            Abs(Scan_Range(2,1)),             !|
     *            Abs(Scan_Range(2,2)))             !V
                                                    !|
        qx_eps = 0.5*k4*UnitCoef(1)*UnitCoef(1)*t*t !|
        qz_eps = 2.*k4*UnitCoef(1)*t                !|
        goto 610    !---------------------->---------+--+
                                                    !V  |
c+----- Mode4: over qx at fixed qz ---|              |  |
  602   continue  !<---------------------------------+  |
        qx_eps = Amax1(Abs(Scan_Range(1,1)),           !|
     *                 Abs(Scan_Range(2,1)))           !V
        qz_eps = Amax1(Abs(Scan_Range(1,2)),           !|
     *                 Abs(Scan_Range(2,2)))           !|
        goto 610    !---------------------->------------|
                                                       !V
  610   continue   !<-----------------------------------+
        qx_eps = (1.e-4)*qx_eps
        qz_eps = (1.e-4)*qz_eps

c
c#######################################################
c        Loops over (Theta(1),Theta(2)) or (qx,qz)
c#######################################################
c cccc  Call    DivZ()
c YES message from DivZ!
        if (GRD_Output) Then !----------------------------+
        if (icont.eq.0 .and. jcont.eq.0) Then  !--------+ |
          Write (11,56,err=191) nscan(1),nscan(2),     !| |
     *                          Scan_Range(1,1),       !| |
     *                          Scan_Range(2,1),       !| |
     *                          Scan_Range(1,2),       !| |
     *                          Scan_Range(2,2),       !| |
     *                             0.,0.               !| |
  56      Format ('DSAA'/2I6,3(/2G15.7))               !| |
        endif  !----------------------------------------+ |
        endif  !------------------------------------------+

        do      j=1,nscan(2)  !============================+
                                                          !|
          if (j.le.jcont) goto 502  !-------continuation---+-+
          Theta2 = Scan_Range(1,2)+Scan_Step(2)*(j-1)     !| V
                                                          !|
        do      i=1,nscan(1)  !==========================+ |
                                                        !| |
          if (i.le.icont) goto 501  !-------continuation-+-++
          icont     = 0                                 !| ||
          jcont     = 0                                 !| |V
          Theta1 = Scan_Range(1,1)+Scan_Step(1)*(i-1)   !| |
                                                        !| |
          Ok_angles = .True.                            !| |
          ifail     = 0                                 !| |
          nfa       = 0                                 !| |
          goto (701,702,703,704) mode_scan  !----------+ | |
c                1   2   3   4                        !| | |
c+------- Mode1: over Q  at 2Q=fix -----|              | | |
  701     continue                                    !| | |
          Theta(1) = Theta1                           !| | |
          Theta(2) = Theta2 - Theta(1)                !| | |
          goto 711                                    !| | |
                                                      !| | |
c+------- Mode2: Q-2Q serie with offsets|              | | |
  702     continue                                    !| | |
          Theta_Smp = Theta1                          !| | |
          Theta_Off = Theta2                          !| | |
c Here the offset is assumed positive if the angle of  | | |
c incidence is increased:                              | | |
          Theta(1) = Theta_Smp + Theta_Off            !| | |
          Theta(2) = Theta_Smp - Theta_Off            !| | |
          goto 711                                    !| | |
                                                      !| | |
c+------- Mode3: over 2Q  at Q=fix -----|              | | |
  703     continue                                    !| | |
          Theta_Smp = Theta2                          !| | |
          Theta_Cnt = Theta1                          !| | |
          Theta(1) = Theta_Smp                        !| | |
          Theta(2) = Theta_Cnt - Theta_Smp            !| | |
          goto 711                                    !| | |
                                                      !| | |
  711     continue                                    !| | |
          if (Theta(1).lt.0.)     Ok_angles = .False. !| | |
          if (Theta(2).lt.0.)     Ok_angles = .False. !| | |
c ccc     qx     = 0.5*k4*(Theta(1)**2 - Theta(2)**2) !| | |
c ccc*           * UnitCoef(1)*UnitCoef(1)            !| | |
c ccc     qz     = k4*(Theta(1) + Theta(2))           !| | |
c ccc*           * UnitCoef(1)                        !| | |
          qx     =(cos(Theta(2)*UnitCoef(1))          !| | |
     -            -cos(Theta(1)*UnitCoef(1))) * k4    !| | |
          qz     =(sin(Theta(2)*UnitCoef(1))          !| | |
     +            +sin(Theta(1)*UnitCoef(1))) * k4    !| | |
          goto 710                                    !| | |
                                                      !| | |
                                                      !| | |
c+------- Mode4: 3D over (qx,qz) -------|              | | |
  704     continue                                    !| | |
          qx     = Theta1                             !| | |
          qz     = Theta2                             !| | |
          if (       (qz.lt.0.)                       !| | |
     *                  .OR.                          !| | |
     *        (abs(qz).lt.1.E-32 .AND.                !| | |
     *         abs(qx).gt.1.E-32)) Then !------+       | | |
            Ok_angles = .False.               !|       | | |
            Theta(1)    = 0.                  !|       | | |
            Theta(2)    = 0.                  !|       | | |
            goto 710                          !|       | | |
          endif  !-----------------------------+       | | |
          if (qz.gt.0.) Then  !--------------------+   | | |
                                                  !|   | | |
c Approximate equations (small angles!):           |   | | |
c                                                 !|   | | |
c That is in fact sin(Theta(1)):                   |   | | |
c           Theta(1) = (0.5*qz/k4 + qx/qz)        !|   | | |
c    /               / UnitCoef(1)                !|   | | |
c That is in fact sin(Theta(2)):                   |   | | |
c           Theta(2) = (0.5*qz/k4 - qx/qz)        !|   | | |Here is a derivation
c    /               / UnitCoef(1)                !|   | | |(I express q=q/k4):
                                                  !|   | | |------------------
c Exact equations (small and large angles!):       |   | | |qx=cos(Q2)-cos(Q1)
            soo = qz*qz+qx*qx                     !|   | | |qz=sin(Q2)+sin(Q1)
            roo = 1./soo - 1./(4.*k4*k4)          !|   | | |
            if (roo.gt.0.)  Then  !-------------+  |   | | |cos(Q2)=qx+cos(Q1)
c Approximately:                                |  |   | | |sin(Q2)=qz-sin(Q1)
c soo = qz*qz ,                                 |  |   | | |
c roo = sqrt(1/qz*qz)=1/qz ,                    |  |   | | |Sum of squares:
c because k4 >> qz >> qx .                      |  |   | | |qx^2+2qxcos(Q1)+qz^2-2qzsin(Q1)=0
              roo = Sqrt(roo)                  !|  |   | | |
c That is in fact sin(Theta(1)):                |  |   | | |(qz^2+qx^2)-2qzsin(Q1)=-2qxcos(Q1)
              Theta(1) = (0.5*qz/k4 + qx*roo)  !|  |   | | |
     /               / UnitCoef(1)             !|  |   | | |Take the square again:
c That is in fact sin(Theta(2)):                |  |   | | |(qz^2+qx^2)^2-
              Theta(2) = (0.5*qz/k4 - qx*roo)  !|  |   | | |   -4qz(qz^2+qx^2)sin(Q1)+
     /               / UnitCoef(1)             !|  |   | | |   +4qz^2sin^2(Q1)=
                                               !|  |   | | |   =4qx^2(1-sin^2(Q1))
              if (Theta(1).lt.0.)              !|  |   | | |
     *                    Ok_angles = .False.  !|  |   | | |After transformations:
              if (Theta(2).lt.0.)              !|  |   | | |(sin(Q1)-qz/2)^2=
     *                    Ok_angles = .False.  !|  |   | | |  =qx^2*(1/(qz^2+qx^2)-1/4)
            else   !----------------------------|  |   | | |
              Ok_angles = .False.              !|  |   | | |
            endif  !----------------------------+  |   | | |
                                                  !|   | | |
            if (.NOT.Ok_angles)  Then  !--+       !|   | | |
              Theta(1) = 0.              !|        |   | | |
              Theta(2) = 0.              !|        |   | | |
            endif  !----------------------+        |   | | |
          else  !----------------------------------|   | | |
            Theta(1) = 0.                         !|   | | |
            Theta(2) = 0.                         !|   | | |
          endif !----------------------------------+   | | |
                                                      !| | |
  710     continue                               !-----+ | |
                                                        !| |
          if (Ok_angles) Then  !---------------------+   | |
            if (abs(RouMiscut).lt.1.E-32 .AND.      !|   | |
     *          abs(SurMiscut).lt.1.E-32 .AND.      !|   | |
     *            nscan(1).gt.1 .AND.               !|   | |
     *          (mode_scan.eq.1 .OR.                !|   | |
     *           mode_scan.eq.4)) Then !----------+  |   | |
              Call Use_Symmetry8(qz,qx,i,j,      !|  |   | |
     *                           N_Points_Max,   !|  |   | |
     *                           mode_scan,      !|  |   | |
     *                           Scan_Range,     !|  |   | |
     *                           Scan_Step,      !|  |   | |
     *                           DSSL,           !|  |   | |
     *                           k4,UnitCoef(1), !|  |   | |
     *                           qz_eps,qx_eps,  !|  |   | |
     *                           Ok_symmetry)    !|  |   | |
            else   !------------------------------|  |   | |
              Ok_symmetry=.False.                !|  |   | |
            endif  !------------------------------+  |   | |
            if (.NOT.Ok_symmetry)  Then  !---------+ |   | |
              DSSL(i) = DSFM(Theta, N_Used, Rc,   !| |   | |
     *                       ipol, nfa,           !| |   | |
     *                    OutFile(1:lout)//'.esc')!| |   | |
              if (DSSL(i).lt.1.E-32) DSSL(i)=0.   !| |   | |
              if (DSSL(i).gt.1.E+32) goto 165 !----+-+---+-+---+
            endif  !-------------------------------+ |   | |   v
          else  !------------------------------------|   | |
            DSSL(i) = 0.0D0                         !|   | |
            Rc(1)   = 0.0E0                         !|   | |
            Rc(2)   = 0.0E0                         !|   | |
          endif  !-----------------------------------+   | |
                                                        !| |
c If interrupted inside the point calculation:           | |
          if (ifail.eq.-999) Then  !-+                   | |
            Intrp = .True.          !|                   | |
            m     = i-1             !|                   | |
            goto 503    !----->------+------------->-----+-+---+
          endif  !-------------------+                   | |   |
                                                        !| |   |
c If failure in the point calculation:                   | |   V
          if (ifail.ne.0)    Then  !------------------+  | |   |
            write   (3,77)       ifail,i,j,nfa       !V  V V   V
  77        format  (1x,'DSFM last failure code: ',i4,
     *               ' -- at point (',i4,',',i4,').',
     *               '  Non-critical failures: ',i6) !^  ^ ^   V
            if (iBatch_mode .lt. 2) Then !----+       |  | |   |
               write (*,77) nfa,ifail,i,j    !|       |  | |   |
            else !----------------------------+       |  | |   |
               if (nfa.eq.0) goto 166 !-------+-------+--+-+---+--+ff serious error
            endif  !--------------------------+       |  | |   |  V
            ifail = 0                                !|  | |   |
c If serious error, then abort:                       |  | |   V
            if (nfa.eq.0)  goto 28   !------>---------+->+-+-+ |
          endif  !------------------------------------+  | | V |
                                                        !| |   |
c Write the curve:                                       | |   |
c (this writes even if the angles are not in the range,  | |   |
c but otherwice we have a conflict with the continuation | |   |
c mode):                                                 | |   |
          looped = looped+1                             !| |   |
          if ((Ok_angles) .OR. (iskip.eq.0)) Then  !--+  | |   |
                                                     !|  | |   |
            if (nscan(2).eq.1)     Then  !----------+ |  | |   |
              write (1,888,err=161) Theta1,DSSL(i) !| |  | |   |
            elseif (nscan(1).eq.1) Then  !----------| |  | |   |
              write (1,888,err=161) Theta2,DSSL(i) !| |  | |   |
            else  !---------------------------------| |  | |   |
              write (1,888,err=161) Theta1,Theta2, !| |  | |   |
     *                              DSSL(i)        !| |  | |   |
            endif  !--------------------------------+ |  | |   |
                                                     !|  | |   |
          else   !------------------------------------|  | |   |
                                                     !|  | |   |
            iskipped = iskipped+1                    !|  | |   |
                                                     !|  | |   |
          endif  !------------------------------------+  | |   |
                                                        !| |   V
          if (iDebug.ne.0) write (33,*) Theta(1),       !| |   |
     *                                  Theta(2),       !| |   |
     *                                  qz,qx,          !| |   |
     *                                  N_Used          !| |   |
                                                        !| |   V
c Print progress report:                                !| |   |
          jjj = j                                       !| |   |
          iii = i + nscan(1)*(jjj-1)                    !| |   |
          if (iii.eq.(iii/iBeepStep)*iBeepStep) Then !+  | |   |
            r = 100.*iii / Full_Points               !|  | |   |
            if (r .gt. 99.99) r = 99.99              !|  | |   |
            Call  Duration (start,bias,hh,mm,ss)     !|  | |   |
            Call Convert_Time (hh,mm,ss,timestamp)   !|  | |   |
            if (iBatch_mode .lt. 2)                  !|  | |   |
     *         write (*,333) i, j, r, timestamp      !|  | |   |
            if (modebat.ne.0)   Then  !-------------+ |  | |   |
              Call OpenFile(DoneFile(1:lout+4),50, !| |  | |   |
     *            'write','unknown',io_status,*93) !| |  | |   |
              write(50,90,err=92) i, j, nscan,     !| |  | |   |
     *                            r, timestamp     !| |  | |   |
  92          continue                             !| |  | |   |
              close(unit=50)                       !| |  | |   |
            endif  !--------------------------------+ |  | |   |
  93        continue                                 !|  | |   |
          endif  !------------------------------------+  | |   V
                                                        !| |   |
c Verify for user interruption:                         !| |   |
          Call  KeyIn (0,iscan,iasci)                   !| |   |
          if (iasci.eq.27)   Then  !-------------------+ | |   |
            Intrp = .True.                            !| | |   |
            m     = i       !this data pt is written   | | |   |
            goto 503    !----->------------------>-----+-+-+---|
          endif  !-------------------------------------+ | |   |
          if (FileExist(OutFile(1:lout)//'.esc'))       !| |   |
     *                                          Then !-+ | |   |
            Call DelFile (OutFile(1:lout)//'.esc',kkk)!| | |   |
            Intrp = .True.                            !| | |   |
            m     = i         !this data pt is written | | |   |
            goto 503    !------->------------------>---+-+-+---|
          endif  !-------------------------------------+ | |   |
                                                        !| |V  |
  501     continue  !<------------continuation-----------+-++  V
        enddo  !=========================================+ |   |
        m = nscan(1)                                      !|   |
  503   continue  !<--------------interrupted--<-----------+---+
        if (GRD_Output) Then !---------------------+       |
          if (m.gt.0)                             !|       |
     *    Write (11,888,err=191) (DSSL(i),i=1,m)  !|       |
  888     format (1001G12.5)                      !|       |
        endif  !-----------------------------------+       | |
        if (Intrp)  goto 504  !------------>---------------+-+--+
  502   continue  !<--------------continuation-------------+-+  |
        enddo  !===========================================+    |
                                                               !V
c ===========================================================   |
                                                               !|
  504   continue  !<--------------------------------------------+
        if (i.gt.nscan(1)) i=nscan(1)
        if (j.gt.nscan(2)) j=nscan(2)
        r = 100.*(i+nscan(1)*(j-1.))/Full_Points
        Call  Duration (start,bias,hh,mm,ss)
  505   continue
        Close   (unit=1)
        Call Convert_Time (hh,mm,ss,timestamp)
        write   (3,333) i, j, r, timestamp
        write   (3,334) iskipped, looped
  334   format  (' Points skipped because of negative',
     *           ' incidence/exit angles =',i5,
     *           ' (out of',i5,' processed)')
        Close   (unit=3)
        if (modebat.ne.0)       Then  !---------------+
          Call OpenFile(DoneFile(1:lout+4),50,       !|
     *              'write','unknown',io_status,*95) !|
          write (50,90,err=94) i, j, nscan,          !|
     *                         r, timestamp          !|
  94      continue                                   !|
          Close(unit=50)                             !|
        endif  !--------------------------------------+
  95    continue

        if (looped.gt.0 .AND. iskipped.eq.looped) Then  !-----------+
          write (txt,121)  ProgramVer(1:ipv)                       !|
  121     format(a,':'//
     *    'No output data written because all the points in'/
     *    'the scan correspond to negative incidence/exit angles') !|
          Call  Message (txt,4,2)                                  !|
        endif  !----------------------------------------------------+

        if (.NOT.Intrp) Then  !----------------+
          do      i=1,1   !==================+ |
            Call  KeyIn   (0,iscan,iasci)   !| |
            if (iscan.ne.0)       goto 28   !| |
          enddo  !===========================+ |
        endif  !-------------------------------+

  28    continue
        Call exit_quiet()
c #######################################################
c                  Error messages
c########################################################
  161   continue
        write (txt,120)  ProgramVer(1:ipv),
     *                   OutFile(1:lout)//'.dat'
  120   format  (a,': Cannot write to output file.'//a//
     *  'May be, device full or drive not ready.')
        Call    Message (txt,5,2)
        Close   (unit=1,status='delete',err=28)
        goto 28

  191   continue
        write (txt,120)  ProgramVer(1:ipv),
     *                   OutFile(1:lout)//'.grd'
        Call    Message (txt,5,2)
        Close   (unit=11,status='delete',err=28)
        goto 28

  160   continue
        write (txt,130)  ProgramVer(1:ipv),
     *                   OutFile(1:lout)//'.dat'
  130   format  (a,': Cannot open output file.'//a)
        Call    Message (txt,3,2)
        goto 28

  190   continue
        write (txt,130)  ProgramVer(1:ipv),
     *                   OutFile(1:lout)//'.grd'
        Call    Message (txt,3,2)
        goto 28

  162   continue
        write (txt,132)  ProgramVer(1:ipv),
     *                   OutFile(1:lout)//'.dat'
  132   format  (a,': Error reading output file.'//a)
        Call    Message (txt,3,2)
        goto 28

  192   continue
        write (txt,132)  ProgramVer(1:ipv),
     *                   OutFile(1:lout)//'.grd'
        Call    Message (txt,3,2)
        goto 28

  193   continue
        write (txt,133)  ProgramVer(1:ipv), idat, jcont, icont,
     *                                      igrd, jgrd,  i
  133   format  (a,': Data mismatch in the continuation mode:'//
     *  i6,' data points in the previously created DAT file'/
     *  ' (',i4,' full grid lines + 1 line with ',i4,' points)'/
     *  ' -- while'/
     *  i6,' data points in the previously created GRD file'/
     *  ' (',i4,' full grid lines + 1 line with ',i4,' points)')
        Call    Message (txt,7,2)
        goto 28

  201   continue
        write (txt,202)  ProgramVer(1:ipv)
  202   format  (a,': Continuation mode error.'//
     *           'The number of lines in the existing'/
     *           'DAT file is more than needed.'/
     *           'Possibly the input file has been altered'/
     *           'before the start of continuation')
        Call    Message (txt,6,2)
        goto 28

  165   continue
        write (txt,205)  ProgramVer(1:ipv), DSSL(i),
     *                   i, j, Theta1, Theta2
  205   format  (a,': Unreasonably big scatterring data = ',g12.5/
     *           'calculated at scan pts = (',i4,',',i4,'), '/
     *           'scan angles = (',g12.5,',',g12.5,').'/
     *           'Perhaps too high roughness sigma for one of layers'/
     *           'breaks the model assumption that the X-ray'/
     *           'wavefields do not change within roughness height.'/
     *           'Please, reduce the rms height and try again.'//
     *           'Contact the author if you believe it is an error!')
        Call    Message (txt,9,2)
        goto 28

  166   continue
        write (txt,206)  ProgramVer(1:ipv), ifail,
     *                   i, j, Theta1, Theta2
  206   format  (a,': Critical error = ',i4,' from a subroutine'/
     *           'encounted at scan pts = (',i4,',',i4,'), '/
     *           'scan angles = (',g12.5,',',g12.5,').'//
     *           'Please report this error to the author!')
        Call    Message (txt,5,2)
        goto 28

        End
