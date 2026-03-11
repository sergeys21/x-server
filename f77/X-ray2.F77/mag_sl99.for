        Program mag_sl99
c       USE MSFLIB                                      !for QuickWin
c ---------------------------------------------------------------
c   X-ray total external reflection from MAGNETIC multilayers
c      +--------------------------------------------------+
c      | This version works for both soft and hard x-rays |
c      | (no approximations) -- see the paper with Sinha. |
c      +--------------------------------------------------+
c
c                         By S.A.Stepanov
c
c *** Version'99 optimized for WWW access and batch processing ***
c
c '99.1 Derived from ter_sl97 and gid_sl97: script language for
c       top-layer profile allowing to specify periodic structures;
c       improved batch ("silent" or "server") mode for using with
c       CGI-interfaces.
c       The chi_ij matrix is 2x2 (hard x-rays with grazing incidence,
c       no X-component of electric fields)
c '99.2 Extension of Version '99.1 to soft X-rays: chi_ij matrix is
c       3x3 and the X-component of electric fields is taken into
c       account.
c '99.3 Added Henke/Cowan option for X0h database
c '99.4 Added option to scan over qz
c '99.5 Added initializing variables (zeroing them) for GNU Fortran
c ---------------------------------------------------------------
c+===============================+
        Include 'mag_sl99.inc'  !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter       (N_Top_Max     = 1001,
c cc *                   N_Total_Max   = N_Top_Max+1,
c cc *                   N_TT_Max      = N_Top_Max+2)
c cc    Parameter       (N_Tran_Max    = 101)           !50

        Integer          N_Frac_Max
        Parameter       (N_Frac_Max    = 4)

        Complex*16      E0(2), Es(2),
     *                  Dy(4, N_TT_Max),        !sigma-Fields in layers
     *                  Dz(4, N_TT_Max),        !pi-Fields in layers
     *                  Dx(4, N_TT_Max),        !longitudinal Fields in layers
     *                  u (4, N_TT_Max)         !dispersion eqn.roots
     *

        Complex*8       x0_, x0_f,
     *                  F10(N_TT_Max),          !(3*lambda/8*pi)*F_10/e_radius
     *                  F11(N_TT_Max),          !(3*lambda/8*pi)*F_11/e_radius
     *                  F1T(N_TT_Max),          !(3*lambda/8*pi)*F_1T/e_radius
     *                  A8, B8, C8,
     *                  Im1 /(0.,1.)/

        Real*8          k0                      ! = 2*pi/wave

        Real*4          R0s,                    !integral reflectivity
     *                  Phase_shift,            !between Es(sigma) & Es(pi)
     *                  Circ_pm(2),             !reflectivity for +-circ.incident polar.
     *                  Circ_Diff,              !difference in +-circ.reflectivity
     *                  Circ_Ratio,             !ratio (I+ - I-)/(I+ + I-)
     *                  Mvect(3,N_TT_Max),      !magnetic vector coordinates
     *                  MvLen,                  !magnetic vector length
     *                  xm_max(N_TT_Max),x4,    !max magnetic part in chi_ij (%)
     *                  Mdensity(N_TT_Max),     !magnetic atoms density (1/cm^3)
     *                  Mshare(N_TT_Max),       !magnetic atoms share of all atoms
     *                  Thickness_Top,
     *                  f0, f0_range(2), df0,   !incidence angle data
     *                  f0_,                    !sin(f0)
     *                  UnitCoef, xxx,
     *                  ter_Angle,
     *                  rho(N_Total_Max),       !material density (g/cm^3)
     *                  W0(N_Total_Max),
     *                  Frac(N_Frac_Max,N_Total_Max),
     *                  SkewPolAngle,
     *                  e_radius/2.81794e-05/,  !Classical electron radius(A)
     *                  pi, AngDegree, AngMinute,
     *                  AngSec, AngMrad, AngUrad,
     *                  xr0, xi0, rho_, rho_f,
     *                  adn_, adn_f, adn_substrate,
     *                  coeff, Trans_mean,
     *                  Sigma_mean, xm_mean,
     *                  Sigma_max, qz_max, f0_max, SQZ/10./,
     *                  Energy, z, DDD, Start,
     *                  progress,
     *                  Energy_max /10./        !threshold (keV) to use this program.
                                                !For higher energies one must use the
                                                !complimentary algorithm designed for
                                                !hard X-rays.

        Real            rpar(2)

        Logical         YesTr

        Integer*4       hh,mm,ss,bias

        Integer         Law_Index,
     *                  nf0, iUnits,
     *                  N_Frac(N_Total_Max),
     *                  Identical(N_Total_Max),
     *                  iBeepStep, lyr, i, j, m,
     *                  ipv, linp, lout, ltop,
     *                  narg, io_status,
     *                  line, iBatch_mode, ixway,
     *                  iPol_mode, iSyngony, mode_x0h,
     *                  iReqFields, nf,
     *                  nerror, iscan, iasci, iarg,
     *                  i_disableABC(3),
     *                  ipar(1)

        Character       Code(N_Frac_Max,N_Total_Max)*20,
     *                  C00(N_Total_Max)*2,
     *                  Radiat*6, wrk*256,
     *                  InpFile*80, OutFile*80,
     *                  File_Top*80, DoneFile*80,
     *                  ProgramVer*8, Version(2)*80,
     *                  timestamp*9, Comment(3)*80,
     *                  LawName(2)*6 /'[cos]', 'Linear'/,
     *                  YesNo(2)*3   /'No',    'Yes'/,
     *                  Spm(5)*9  /'Sigma', 'Pi', 'Skew',
     *                             'Circular-', 'Circular+'/,
     *                  Uni(6)*5  /'degr.','min.',
     *                             'mrad.','sec.',
     *                             'urad.','1/A'/,
     *                  Mtype(4)*7  /' None', '  Yes',
     *                               'Along_Y', 'Along_Z'/
c    *                  Validation*16, Validcode*16/'abracadabra'/

        Real            wave2energy
        External        wave2energy

        Logical*4       FileExist
        External        FileExist

        Integer*4       Nargs
        External        Nargs

        Integer         Levi
        External        Levi

        Character       txt(20)*80
        Common  /msg/   txt

        Complex*16      x_ij(3,3,N_TT_Max)
        Complex*8       x0(N_TT_Max)
        Real*4          Wave,
     *                  xabs,                           ! abs(x0_max)
     *                  TERang,                         ! critical TER angle (radians)
     *                  Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Thickness_Tr(N_Total_Max),
     *                  TranLaw(N_Tran_Max)
        Integer         Magnetic(N_TT_Max),
     *                  N_Top, N_Total,
     *                  N_TT, N_Used, N_Tran,
     *                  MaxOrder, iDebug, ifail
        Common  /MAGdat/  x_ij,                         ! (3,3,N_TT_Max)
     *                    x0,                           ! (N_TT_Max)
     *                    Wave,
     *                    xabs,                         ! abs(x0_max)
     *                    TERang,                       ! critical TER angle (radians)
     *                    Thickness,                    ! (N_Total_Max),
     *                    Sigma,                        ! (N_Total_Max),
     *                    Thickness_Tr,                 ! (N_Total_Max),
     *                    TranLaw,                      ! (N_Tran_Max),
     *                    Magnetic,                     ! (N_TT_Max),
     *                    N_Top, N_Total,
     *                    N_TT, N_Used, N_Tran,
     *                    MaxOrder, iDebug, ifail

c       Parameter       (kcompMax = 10)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+

c We must clear prcn before each new call of X0h, otherwice
c the of previous material are used:
        Real*4          prcn(kcompMax)
        Common  /x0pa3/ prcn

c This is an interface to internal X0h data.
c We use atom density (1/cm^3) from this common:
        Real*4          ucell_mass_gram, atom_density_cm3
        Integer         n_atoms_ucell
        Common  /x0paA/ ucell_mass_gram, atom_density_cm3, n_atoms_ucell

c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer             iHenkeCowan
        Common /HenkeCowan/ iHenkeCowan

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
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Integer         modebat, lunbat, ierrbat, istackrezv
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        istackrezv = 0
        progname = 'mag_sl99'
c       i = SetExitQQ(2)        !QWIN$EXITNOPERSIST=2   !for QuickWin

c-------------------------------
        i = kcompMax                    !to suppress unused paramters warnings
        i = kolMax                      !to suppress unused paramters warnings
        i = ncodMax                     !to suppress unused paramters warnings
        i = natmMax                     !to suppress unused paramters warnings
        i = nwavMax                     !to suppress unused paramters warnings
        i = nedgMax                     !to suppress unused paramters warnings
        i = nabcMax                     !to suppress unused paramters warnings
        i = nedgTb                      !to suppress unused paramters warnings
        i = maxX0hWarnLines             !to suppress unused paramters warnings
c #######################################################
        ifail      = 0
        pi         = 4.*atan(1.)
        AngDegree  = 2.*pi/360.
        AngMinute  = AngDegree/60.
        AngSec     = AngDegree/3600.
        AngMrad    = 1.E-03
        AngUrad    = 1.E-06

        do      lyr=1,N_Total_Max  !=====+
          rho(lyr)    = 0.              !| material density (g/cm^3)
          W0(lyr)     = 0.              !|
          C00(lyr)    = ' '             !|
          N_Frac(lyr) = 1               !|
          do    j=1,N_Frac_Max  !=====+  |
            Frac(j,lyr) = 0.         !|  |
            Code(j,lyr) = ' '        !|  |
          enddo  !====================+  |
          Frac(1,lyr) = 1.              !|
        enddo  !=========================+

        do lyr=1,N_TT_Max !================+
          Magnetic(lyr) = 0               !|
          xm_max(lyr)   = 0.              !| max magnetic part in chi_ij (%)
          Mdensity(lyr) = 0.              !| magnetic atoms density (1/cm^3)
          Mshare(lyr)   = 0.              !| magnetic atoms share of all atoms
          F10(lyr)      = (0.,0.)         !| (3*lambda/8*pi)*F_10/e_radius
          F11(lyr)      = (0.,0.)         !| (3*lambda/8*pi)*F_11/e_radius
          F1T(lyr)      = (0.,0.)         !| (3*lambda/8*pi)*F_1T/e_radius
          x0(lyr)       = (0.,0)          !|
          do i=1,3 !====================+  |
            Mvect(i,lyr) = 0.          !|  | magnetic vector coordinates
            do j=1,3 !================+ |  |
              x_ij(i,j,lyr) = (0.,0) !| |  |
            enddo !===================+ |  |
          enddo !=======================+  |
        enddo !============================+

        ProgramVer = 'mag_sl99'
        ipv = Len_Trim (ProgramVer)
c This is the startup ERR-file (it is used until we've
c verified the output filename):
        ErrorFile = ProgramVer(1:ipv)//'.err'
        Call    DeleteErrorFile ()

        write   (Version,111)
  111   format  ('MAG_sl: X-Ray Total Reflection from Magnetic ',
     *  'Multilayers with Rough Interfaces.'/
     *  'By Sergey Stepanov <sstepanov@anl.gov>>',14x,
c    *  'Version: 99_4.3,  Nov-2006')
c    *  'Version: 99_4.4,  Sep-2012')
     *  'Version: 99_4.5,  Nov-2013')

        narg = Nargs()-1
        if (narg .eq. 0)  Then  !-------------------+
                                                   !|
          InpFile = ProgramVer(1:ipv)//'.inp'      !|
                                                   !|
        else  !-------------------------------------+
                                                   !|
          iarg = 1                                 !|
          Call get_arg_wrapper (iarg,InpFile)      !|
                                                   !|
        endif  !------------------------------------+
        linp = Len_Trim(InpFile)

        modebat = 1     !assume batch-mode by default
        iBatch_mode = modebat

c #######################################################
c              Read input file
c#######################################################
        Call OpenFile(InpFile(1:linp),2,'read','old',io_status,*99)

        Call Make_Err_Filename (InpFile(1:linp))

        line   = 0
        txt(7) = 'batch mode flag'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt    (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        iBatch_mode = ipar(1)
        if (iBatch_mode .lt. 0)   goto 101
        if (iBatch_mode .gt. 2)   goto 101
        if (iBatch_mode .eq. 0) iBatch_mode=1
c       if (iBatch_mode .eq .2) Then  !-----------------+
c         read(*,'(a)') Validation                     !|
c         if (Validation .ne. Validcode)  Call Check() !|
c         iBatch_mode = 1                              !|
c       endif  !----------------------------------------+
        modebat = iBatch_mode
c-----------------------
        txt(7)  = 'output file name'
        Call    LineLoop (2,line,wrk,*100)
        lout    = Min(Len(OutFile),Len(wrk))
        OutFile = wrk(1:lout)
        lout    = Len_Trim(OutFile)
        if (lout .eq .0)  Then  !---+
          OutFile = '~'            !|
          lout    = 1              !|
        endif  !--------------------+
c       Call    Case    (OutFile(1:lout),1)     !to lower case
        Call    OpenList(OutFile(1:lout)//'.tbl',3,ifail)
        if (ifail .ne. 0)         goto 28
        lunbat = 3
        Call    OpenList (OutFile(1:lout)//'.dat',1,ifail)
        if (ifail .ne. 0)         goto 28
c-----------------------
        txt(7) = 'surface profile name'
        N_Top   = 0
        N_Total = 0
        Call    LineLoop (2,line,wrk,*100)
        ltop = Min(Len(File_Top),Len(wrk))
        File_Top = wrk(1:ltop)
        ltop = Max (Len_Trim(File_Top),1)
c       Call    Case (File_Top(1:ltop),1)       !to lower case
c Allow input of transition layers in the
c Read_Top99_M0 (this might be disabled in
c the diffuse scattering program):
        YesTr = .True.
        txt(7) = 'surface layer profile'
c For vacuum:
        x0(1)       = (0., 0.)
        F10(1)      = (0., 0.)
        F11(1)      = (0., 0.)
        F1T(1)      = (0., 0.)
        Mdensity(1) = 0.
        Mshare(1)   = 0.
        do i=1,3  !========+
          Mvect(i,1) = 0. !|
        enddo  !===========+
c-----------------------
c Now, read the top profile (.prf) -- see: Rd99prfM.for.
c We need to read the profile before the substrate data because the later
c is added at the end of profile. If we want to read the profile after we
c are done with the .inp file, we can place the substrate data into
c N_Total_Max (into the last element of the array) and then  move it to
c the N_Total position after the profile is read. That is how it is
c implemented in ter_sl.
        Call    Read_Top99_M0 (File_Top(1:ltop),
     *                         N_Top, N_Top_Max,
     *                         Thickness, Thickness_Top,
     *                         Code(1,1), rho(1),
     *                         Frac(1,1), N_Frac(1), N_Frac_Max,
     *                         C00(1), W0(1),
     *                         x0(2), F10(2), F11(2), F1T(2),
     *                         Mdensity(2), Mshare(2), Mvect(1,2),
     *                         Sigma(1), Thickness_Tr(1), YesTr,
     *                         Identical(1), ifail)
        if (ifail .ne. 0)         goto 28
        N_Total = N_Top+1               !layers+substrate
        N_TT    = N_Top+2               !layers+substrate+vacuum
c-----------------------
c Back to reading .inp file:
        do i=1,3  !=====================================+
          write (txt(7),'(a,i1)') 'comment line #',i   !|
          Call  LineLoop (2,line,Comment(i),*100)      !|
        enddo  !========================================+
c-----------------------
        txt(7) = 'substrate data'
        Thickness(N_Total) = 0.
        W0(N_Total)        = 1.
        YesTr = .True.
c See: Rd99prfM.for:
        Call Read_Sub99_M0 (Code(1,N_Total), rho(N_Total),
     *                      C00(N_Total), W0(N_Total),
     *                      x0(N_TT), F10(N_TT), F11(N_TT), F1T(N_TT),
     *                      Mdensity(N_TT), Mshare(N_TT), Mvect(1,N_TT),
     *                      Sigma(N_Total), Thickness_Tr(N_Total),
     *                      YesTr, wrk, line, *100, *101, *28)
c-----------------------
        txt(7) = 'x-ray specification mode'
        Call    LineLoop (2,line,wrk,*100)
c 1=wave 2=energy 3=x-ray line:
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        ixway = ipar(1)
        if (ixway .lt. 1)         goto 101
        if (ixway .gt. 3)         goto 101
c-----------------------
        txt(7) = 'x-ray wavelength/energy (range=[0.01-1000.0]Angstrom)'
        Call    LineLoop (2,line,wrk,*100)
        Wave = 0.
        if (ixway .lt. 3) Then  !---------------------+
          Call  rdReal  (rpar,1,wrk,i)               !|
          if (i .ne. 0)           goto 100           !|
          Wave = rpar(1)                             !|
          if (Wave .le. 0.)       goto 101           !|
          if (ixway .eq. 2) Wave=wave2energy(Wave)   !|
          if (Wave .lt. 0.01)     goto 101           !|
          if (Wave .gt. 1000.)    goto 101           !|
        endif  !--------------------------------------+
c-----------------------
        txt(7) = 'x-ray line name'
        Call    LineLoop (2,line,wrk,*100)
        if (ixway.eq.3) Then  !-------------------+
          Radiat = wrk(1:6)                      !|
          if (Len_Trim(Radiat).eq.0) goto 131 !---+---+
c Get wavelength from the X-ray line name:        |   v
          Call  Get_Wave (Wave,Radiat,ifail)     !|
          if (ifail.ne.0)         goto 130       !|
          if (Wave .lt. 0.01)     goto 101       !|
          if (Wave .gt. 1000.)    goto 101       !|
        else !------------------------------------|
          Radiat = ' '                           !|
        endif  !----------------------------------+
        k0 = 2.* pi/wave
c-----------------------
c 1. Sigma    4. Circular-
c 2. Pi       5. Circular+
c 3. Skew
        txt(7) = 'incident x-ray polarization mode [1-5]'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        iPol_mode = ipar(1)
        if ((iPol_mode .lt. 1) .OR. (iPol_mode .gt. 5)) goto 101
c-----------------------
        txt(7) = 'Angle between x-ray polarization & sigma [+-360]'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (rpar,1,wrk,i)
        if (i .ne. 0)             goto 100
        SkewPolAngle = rpar(1)
        if (Abs(SkewPolAngle) .gt. 360.) goto 101
        if (iPol_mode .eq. 1) SkewPolAngle = 0.
        if (iPol_mode .eq. 2) SkewPolAngle = 90.
        if (iPol_mode .lt. 4) Then  !------------+
c 1.Sigma  2.Pi  3.Skew                          |
          E0(1) = cos(SkewPolAngle*AngDegree)   !|
          E0(2) = sin(SkewPolAngle*AngDegree)   !|
        elseif  (iPol_mode .eq. 4) Then  !-------|pi=-i*sigma +----------------+
c 4.Circular-                                    |            |What is the sca-|
          E0(1) = (1.+Im1)/2.                   !|            |ling factor?    |
          E0(2) = (1.-Im1)/2.                   !|            |I use the norma-|
          SkewPolAngle = 0.                     !|            |lization where  |
        elseif  (iPol_mode .eq. 5) Then  !-------|pi=+i*sigma ||Es|^2+|Ep|^2=1 |
c 5. Circular+                                   |            |----------------|
          E0(1) = (1.+Im1)/2.                   !|            |Sunil suggested |
          E0(2) =-(1.-Im1)/2.                   !|            |/sqrt(2) instead|
          SkewPolAngle = 0.                     !|            |of /2 ??????????|
        else  !----------------------------------|            +----------------+
          E0(1) =(0.,0.)                        !|
          E0(2) =(0.,0.)                        !|
          SkewPolAngle = 0.                     !|
        endif  !---------------------------------+
c-----------------------
c Get substrate x0:
        if (C00(N_Total) .ne. 'x')  Then !-------------------+
          mode_x0h = 0                          !new crystal |
          iSyngony = 0                                      !|
          do i=1,kcompMax !==+                               |
            prcn(i) = 0.    !|                               |
          enddo  !===========+                               |
          Call Get_X0 (Wave, Code(1,N_Total), rho(N_Total), !|
     *                 iSyngony, mode_x0h, xr0, xi0,        !|
     *                 x0(N_TT), 1, 1, ifail)               !|
          if (ifail .ne. 0)   goto 130                      !|
c adn = means "Atomic Density" (in 1/cm^3)                   |
          adn_substrate = atom_density_cm3                  !|
        else  !----------------------------------------------|
c adn = means "Atomic Density" (in 1/cm^3)                   |
          adn_substrate = 0.                                !|
        endif !----------------------------------------------+
        Identical(N_Total) = 0
c N_TT=N_Top+2 -- layers+substrate+vacuum:
        if (Mshare(N_TT) .gt. 0.) Then  !--------------------+
          if (abs(adn_substrate).gt.1.E-20) Then !--------+  |
            Mdensity(N_TT) = Mshare(N_TT)*adn_substrate  !|  |
          else  !-----------------------------------------|  |
            lyr = N_Total                                !|  |
            j   = 1                                      !|  |
            goto 125  !------error: no density data-------+--+---+
          endif  !----------------------------------------+  |   v
        endif  !---------------------------------------------+
c-----------------------
        txt(7) = 'number of transition sublayers'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)               goto 100
        N_Tran = ipar(1)
        if (N_Tran .lt. 1) N_Tran = 1
        if (N_Tran .gt. N_Tran_Max) goto 101
c-----------------------
        txt(7) = 'transition layers type'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)               goto 100
        Law_Index = ipar(1)
        if (Law_Index .lt. 0)       goto 101
        if (Law_Index .gt. 1)       goto 101
c-----------------------
        txt(7) = 'angular units'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)               goto 100
        iUnits = ipar(1)
        if (iUnits*(5-iUnits) .lt. 0) goto 101
        goto (40,41,42,43,44,45) iUnits+1 !->-+
                                             !|
  40    continue                   !<---------|
        UnitCoef = AngDegree                 !|
        goto 49   !->----------------------->-+--->+
  41    continue                   !<---------|    |
        UnitCoef = AngMinute                 !|    |
        goto 49   !->----------------------->-+--->|
  42    continue                   !<---------|    |
        UnitCoef = AngMrad                   !|    |
        goto 49   !->----------------------->-+--->|
  43    continue                   !<---------|    |
        UnitCoef = AngSec                    !|    |
        goto 49   !->----------------------->-+--->|
  44    continue                   !<---------|    |
        UnitCoef = AngUrad                   !|    |
        goto 49   !->----------------------->-+--->|
  45    continue                   !<---------+    |
        UnitCoef = sngl(1./(2.*k0))               !| for qz-scan: sin(f0)=qz/(2*k0)
        goto 49   !->----------------------->----->|
                                                  !|
  49    continue   !<-----------------<-----------<+
c-----------------------
        if (iUnits .ne. 5) Then  !-------------------------------+angular scan
          txt(7) = 'incidence angle limits (range=[0:90]degr.)' !|
        else  !--------------------------------------------------|qz-scan
c qz = 2k*sin(Theta) => qz_max = 2k:                             |
          write (txt(7),112) 2.*k0                              !|
  112     format('qz limits ( range = [ 0 : ',g8.3,'] (1/A) )') !|
        endif  !-------------------------------------------------+
        Call    LineLoop (2,line,wrk,*100)
c Read incidence angle range:
        Call    Rdreal  (f0_range,2,wrk,i)
        if (i .ne. 0)                                  goto 100
        if (f0_range(1) .lt. 0.)                       goto 101
        if (f0_range(2) .lt. 0.)                       goto 101
        if (iUnits .ne. 5) Then  !-------------------------------+angular scan
          if (UnitCoef*f0_range(1) .gt. 90.*AngDegree) goto 101 !|
          if (UnitCoef*f0_range(2) .gt. 90.*AngDegree) goto 101 !|
        else  !--------------------------------------------------|qz-scan
c qz = 2k*sin(Theta) => qz_max = 2k:                             |
          if (f0_range(1) .gt. 2.*k0)                  goto 101 !|
          if (f0_range(2) .gt. 2.*k0)                  goto 101 !|
        endif  !-------------------------------------------------+
c-----------------------
        txt(7) = 'number of incidence angle points (>0)'
        Call    LineLoop (2,line,wrk,*100)
c Read incidence angle points:
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        nf0 = ipar(1)
        if (nf0 .lt. 1)           goto 101
        if (nf0 .gt. 1) Then !-------------------------+
          xxx = f0_range(2) - f0_range(1)             !|
          df0 = xxx / (nf0-1)                         !|
          if (abs(xxx).lt.1.E-20) goto 109            !|
        else  !----------------------------------------|
          df0 = 0.                                    !|
        endif  !---------------------------------------+
c-----------------------
        txt(7) = 'request wavefields flag'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)             goto 100
        iReqFields = ipar(1)
        if (iReqFields .lt. 0)    goto 101
        if (iReqFields .gt. 1)    goto 101
c Printout of wavefields is incompatible
c transition layers:
        if (iReqFields .eq. 1)  Then  !-----------------------+
          do    lyr = 1, N_Total !==========================+ |
            if (abs(Thickness_Tr(lyr)).gt.1.E-20) goto 113 !| |
          enddo  !==========================================+ |
        endif  !----------------------------------------------+
c-----------------------
        txt(7) ='flags to disable A,B,C contributions to resonance x_ij'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (i_disableABC,3,wrk,i)
        if (i .ne. 0)               goto 100
        do i=1,3  !===================================+
          if (i_disableABC(i) .lt. 0)    goto 101    !|
          if (i_disableABC(i) .gt. 1)    goto 101    !|
        enddo  !======================================+

c-----------------------
        MaxOrder  = 14
        iDebug    = 0
        iBeepStep = 10
c-----------------------
        close (unit=2)

c Now, copy INP and PRF files into the TBL:
        Call List_INP(3, InpFile, File_Top, wrk)

c A control for the substrate x0:
        if (abs(Real(x0(N_TT))).lt.1.E-32 .OR.
     *         abs(W0(N_Total)).lt.1.E-32)      goto 123
c #######################################################
c Calculate x0 for the top layer
c (including the case of fractions)
        do      lyr=1,N_Top  !===============================+
          if (Identical(lyr) .eq. 0) Then  !---------------+ |
                                                          !| |
            if ((Len_Trim(Code(1,lyr)) .eq. 0) .OR.       !| |
     *          (Code(1,lyr) .eq. 'Unknown'))  Then !--+   | |
              Code(1,lyr) = Code(1,N_Total)           !|   | |
            endif  !-----------------------------------+   | |
            x0_  = (0., 0.)                               !| |
            rho_ = 0.                                     !| |
            adn_ = 0.                                     !| |
            do j=1,N_Frac(lyr)  !========================+ | |
              if (j .eq. 1)  Then !-+                   !| | |
                rho_f = rho(lyr)   !|                   !| | |
              else  !---------------|                   !| | |
                rho_f = 0.         !|                   !| | |
              endif !---------------+                   !| | |
              if (Code(j,lyr).eq.Code(1,N_Total) .AND.  !| | |
     *             abs(rho_f).lt.1.E-20) Then !--------+ | | |
                x0_f  = x0(N_Total+1)                 !| | | |
c This means, even when the substrate density was     !| | | |
c specified externally, it is assigned to the layer   !| | | |
c with unspecified density:                           !| | | |
                rho_f = rho(N_Total)                  !| | | |
                adn_f = adn_substrate                 !| | | |
              else !-----------------------------------| | | |
c New material, but same wavelength:                   | | | |
                mode_x0h = 1                          !| | | |
                iSyngony = 0                          !| | | |
                do i=1,kcompMax !==+                   | | | |
                  prcn(i) = 0.    !|                   | | | |
                enddo  !===========+                   | | | |
                Call Get_X0 (Wave, Code(j,lyr), rho_f,!| | | |
     *                       iSyngony, mode_x0h,      !| | | |
     *                       xr0, xi0,                !| | | |
     *                       x0_f, 1, 1, ifail)       !| | | |
                if (ifail .ne. 0)      goto 130       !| | | |
                adn_f = atom_density_cm3              !| | | |
              endif !----------------------------------+ | | |
              x0_  = x0_  + x0_f *Frac(j,lyr)           !| | |
              rho_ = rho_ + rho_f*Frac(j,lyr)           !| | |
              if (Mshare(lyr+1).gt.0. .AND.             !| | |
     *               abs(adn_f).lt.1.E-20) goto 125 !----+-+-+-+error
              adn_ = adn_ + adn_f*Frac(j,lyr)           !| | | v
            enddo  !=====================================+ | |
                                                          !| |
            if (C00(lyr) .eq. ' ') Then  !----------+      | |
c x0 was not specified in profile:                  |      | |
              x0(lyr+1) = x0_                      !|      | |
              rho(lyr)  = rho_                     !|      | |
            else  !---------------------------------+      | |
c x0 was specified, than Code is not used:          |      | |
              Code(1,lyr) = 'Unknown'              !|      | |
              Frac(1,lyr) = 1.                     !|      | |
              N_Frac(lyr) = 1                      !|      | |
              adn_      = 1                        !|      | |
              j         = 1                        !|      | |
              if (Mshare(lyr+1) .gt. 0.) goto 125 !-+------+-+-+error
            endif !---------------------------------+      | | v
                                                          !| |
            x0(lyr+1) = x0(lyr+1)*W0(lyr)                 !| |
                                                          !| |
            if (Mshare(lyr+1) .gt. 0.) Then  !-------+     | |
              Mdensity(lyr+1) = Mshare(lyr+1)*adn_  !|     | |
            endif  !---------------------------------+     | |
                                                          !| |
          else  !------------------------------------------| |
                                                          !| |
            i = Identical(lyr)                            !| |
            do j=1,N_Frac(i)  !==========+                 | |
              Code(j,lyr)  = Code(j,i)  !|                 | |
              Frac(j,lyr)  = Frac(j,i)  !|                 | |
            enddo  !=====================+                 | |
            N_Frac(lyr)    = N_Frac(i)                    !| |
            rho(lyr)       = rho(i)                       !| |
            x0(lyr+1)      = x0(i+1)                      !| |
            Mdensity(lyr+1)= Mdensity(i+1)                !| |
c This must be already copied in Read_Top99_M0 :           | |
c           Mshare(lyr+1)  = Mshare(i+1)                  !| |
c           F10(lyr+1)     = F10(i+1)                     !| |
c           F11(lyr+1)     = F11(i+1)                     !| |
c           F1T(lyr+1)     = F1T(i+1)                     !| |
c           do j=1,3  !==================+                !| |
c             Mvect(j,lyr) = Mvect(j,i) !|                !| |
c           enddo  !=====================+                !| |
                                                          !| |
          endif !------------------------------------------+ |
                                                            !|
        enddo  !=============================================+
c Corrections for substrate:
        x0(N_TT) = x0(N_TT)*W0(N_Total)


c Control for big sigma (huge exponents):
        f0_max = Max(Abs(f0_range(1)),Abs(f0_range(2)))
        if (f0_max .gt. 0.) Then  !------------------------+
          if (iUnits .ne. 5) Then  !------------------+    |angular scan
            qz_max=sngl(2.*k0*sin(f0_max*UnitCoef))  !|    |
          else !--------------------------------------|    |qz-scan
            qz_max=f0_max                            !|    |
          endif !-------------------------------------+    |
c The condition like sigma*qz < 1:                         |
          Sigma_max = SQZ/qz_max                          !|
          do lyr=1,N_Total  !=========================+    |
            if (Sigma(lyr) .gt. Sigma_Max) goto 701 !-+----+-+
          enddo  !====================================+    | v
        endif  !-------------------------------------------+

c------------------------------------------------------------------
c Fill x_ij matrix (the polarizabilities with account for resonance
c magnetic scattering) for each layer (*):
c------------------------------------------------------------------
c The factor between x0 and Z*r0 in X0h is (Z=f0):
c       4*Pi   Kol        Wave^2  Kol
c  x0 = ----  ---- r0*Z = ----   ---- r0*Z
c       k0^2   V           Pi     V
c
c where Kol/V is the density (1/cm^3 -> 1/A^3).
c Then, the magnetic scattering amplitude contains its own 3/4*(wave/2*pi)  !!!!!!!!!!!????????????
c------------------------------------------------------------------
c (*) The magnetic data are specified for the structure as a full --
c     not for fractions!
c------------------------------------------------------------------
c This coefficient is for ABSOLUTE values of F11,F1T,F10.
c       coeff = (1.E-24)*3.*(wave**3)/(8.*pi**2)

c In the literature it is common to express (wave/(8*pi))*F11,F1T,F10
c in classical electron radius. Then the coefficient is:
        coeff = (1.E-24)*(wave**2)*e_radius/pi
c               +------+
c                  +this comes from the fact that Mdensity is expressed
c                   in 1/cm**3 while we need it in 1/A**3.

        do      lyr=1,N_TT  !========================================+
c Our coordinate system:                                             |
c 1=k=k0=x, 2=s=sigma=y, 3=p=pi=z_internal:                          |
                                                                    !|
          A8 = coeff*Mdensity(lyr)*(F11(lyr)+F1T(lyr))              !|
          B8 = coeff*Mdensity(lyr)*(F11(lyr)-F1T(lyr))              !|
          C8 = coeff*Mdensity(lyr)*(2.*F10(lyr)-F11(lyr)-F1T(lyr))  !|
                                                                    !|
          if (i_disableABC(1) .eq. 1)  A8 = 0.                      !|
          if (i_disableABC(2) .eq. 1)  B8 = 0.                      !|
          if (i_disableABC(3) .eq. 1)  C8 = 0.                      !|
                                                                    !|
          MvLen = Mvect(1,lyr)**2                                   !|
     +          + Mvect(2,lyr)**2                                   !|
     +          + Mvect(3,lyr)**2                                   !|
          if (MvLen .lt. (1.E-20)) MvLen = 1.                       !|
          MvLen = Sqrt(MvLen)                                       !|
                                                                    !|
          do i=1,3  !=========================================+      |
            do j=1,3  !=====================================+ |      |
c Over columns or lines???                                  | |      |This only
              x_ij(i,j,lyr) = C8*Mvect(i,lyr)*Mvect(j,lyr) !| |      |inverts the
c             x_ij(j,i,lyr) = C8*Mvect(i,lyr)*Mvect(j,lyr) !| |      |difference
     /                      / (MvLen*MvLen)                !| |      |between E+
              do m=1,3  !============================+      | |      |and E- circ.
c Over columns or lines???                           |      | |      |polarizations!
                x_ij(i,j,lyr) = x_ij(i,j,lyr)       !|      | |      |
c               x_ij(j,i,lyr) = x_ij(j,i,lyr)       !|      | |      |
     -                        - im1*B8*Levi(i,j,m)  !|      | |      |
     *                        * Mvect(m,lyr)/MvLen  !|      | |      |
              enddo  !===============================+      | |      |
            enddo  !========================================+ |      |
            x_ij(i,i,lyr) = x_ij(i,i,lyr) + A8               !|      |
          enddo  !============================================+      |
                                                                    !|
c Now make a projection of x_ij onto the plane perpendicular to the  |
c incident/reflected waves direction (onto the sigma_&_pi plane):    |
c                                                                    |
c Our coordinate system:                                             |
c 1=k=k0=x, 2=s=sigma=y, 3=p=pi=z_internal:                          |
          xm_max(lyr) = 0.                                          !|
          do i=1,3  !=====================================+          |
            do j=1,3  !================================+  |          |
              x4 = sngl(Abs(x_ij(i,j,lyr)))           !|  |          |
c Find the maximum magnetic contribution:              |  |          |
              if (x4 .gt. xm_max(lyr)) xm_max(lyr)=x4 !|  |          |
            enddo  !===================================+  |          |
            x_ij(i,i,lyr) = x_ij(i,i,lyr) + x0(lyr)      !|          |
          enddo  !========================================+          |
          x4 = Abs(x0(lyr)) + xm_max(lyr)                           !|
          if (x4 .gt. 0.) Then  !--------------+                     |
            xm_max(lyr) = 100.*xm_max(lyr)/x4 !|Express the          |
          else  !------------------------------| ratio in %.         |
            xm_max(lyr) = 0.                  !|                     |
          endif !------------------------------+                     |
                                                                    !|
          if (xm_max(lyr) .gt. (1.E-5)) Then !-------------------+   |
            if   (Abs(Mvect(2,lyr))/MvLen .gt. 0.9999) Then !--+ |   |
              Magnetic(lyr) = 2                 !special M||Y  | |   |
            elseif (Abs(Mvect(3,lyr))/MvLen .gt. 0.9999) Then !| |   |
              Magnetic(lyr) = 3                 !special M||Z  | |   |
            else  !--------------------------------------------| |   |
              Magnetic(lyr) = 1                 !magnetic      | |   |
            endif !--------------------------------------------+ |   |
          else !-------------------------------------------------|   |
            Magnetic(lyr) = 0                    !non-magnetic   |   |
          endif  !-----------------------------------------------+   |
                                                                    !|
c         if (lyr .eq. 1) Then !--------------------------------+    |
c           write (3,127) lyr, "Vacuum          ", xm_max(lyr) !|    |
c         else !------------------------------------------------|    |
c           write (3,127) lyr, Code(1,lyr-1)     , xm_max(lyr) !|    |
c         endif  !----------------------------------------------+    |
c         write (3,128) ((x_ij(i,j,lyr),i=1,3),j=1,3)               !|
c 127     format(' ----------'/1x,i6,2x,a,g12.3)                    !|
c 128     format(3(1x,3('  (',g11.5,',',g11.5,')  ')/))             !|
                                                                    !|
        enddo  !=====================================================+

        Trans_mean = 0.
        Sigma_mean = 0.
        xm_mean    = 0.
        if (N_Top .gt. 0) Then  !----------------------------+
          do    lyr = 1, N_Top  !==========================+ |
            Trans_mean = Trans_mean + Thickness_Tr(lyr)   !| |
            Sigma_mean = Sigma_mean + Sigma(lyr)          !| |
            xm_mean    = xm_mean    + xm_max(lyr)         !| |
          enddo  !=========================================+ |
          Trans_mean = Trans_mean / N_Top                   !|
          Sigma_mean = Sigma_mean / N_Top                   !|
          xm_mean    = xm_mean    / N_Top                   !|
        endif  !---------------------------------------------+

        Energy = wave2energy(Wave)
        i = Min(lout,40)
        j = Min(ltop,40)
        write (txt,33) OutFile(1:i),
     *                 File_Top(1:j),
     *                 Thickness_Top, N_Top,
     *                 Trans_mean, Sigma_mean,
     *                 xm_mean,
     *                 Code(1,N_Total), Thickness_Tr(N_Total),
     *                                         Sigma(N_Total),
     *                 xm_max(N_TT),
     *                 Wave, Energy, Radiat,
     *                 Spm(iPol_mode), SkewPolAngle,
     *                 Uni(iUnits+1), f0_range,nf0,
     *                 N_Tran, LawName(1+Law_Index)
  33    format  (
     *  'Output files [.dat & .tbl] will be...|',a/
     *  'Top: Profile name....................|',a/
     *  'Top: Thickness(A), N-sublayers.......|',f10.1,i10/
     *  'Top: <TransLayer>, <roughness>.......|',2f10.2/
     *  'Top: <Magnetic_contribution>(%)......|',g10.2/
     *  'Substrate: Code, TransLayer,Roughness|',a,    2f7.2/
     *  'Substrate: Magnetic_contribution (%).|',g10.2/
     *  'X-ray: Wavelength (A), Energy (keV)..|',g12.7,1x,g12.7,1x,a/
     *  'X-ray: polarization, (angle to sigma)|',a,'  (',f6.2,' degr.)'/
     *  'Incidence angle range [',a,'], points|',2g12.3,i8/
     *  'Transition layers: N_sublayers, Shape|',i5,10x,a)

        line = 11

c 2018.10: changed xabs from |x0_substr| to |x0_max|
c       xabs      = abs(x0(N_Total+1))
c       TERang    = sqrt(xabs)
c       TER_Angle = TERang / UnitCoef
c       if (Abs(xabs) .lt. (1.E-24)) xabs = 1.
        Call xabsmax (x0, N_Total_Max, N_Total, UnitCoef,
     *                xabs, TERang, TER_angle)

c This is more-or-less voluntary limit (it could be 6keV, 10keV, and
c etc.), but the bottom line is that one must use the high-energy
c approximation for hard X-rays because present algorithm fails
c due to precision loss errors:
        if (Energy .gt. Energy_max) goto 137
c #######################################################
c                 Open output files:
c#######################################################
        write   (3,'(/33x,a)')  'R E S U L T S'
        write   (3,11)
        do i=1,2  !=======================+
          j = Len_Trim(Version(i))       !|
          write (3,55)  Version(i)(1:j)  !|
        enddo  !==========================+
        write   (3,11)
  11    format  (79('-'))
        write   (3,55)  (txt(i)(1:Len_Trim(txt(i))),i=1,line)
  55    format  (20(a:/))
        Call    Algorithm_Type  (Version(1))
        j = Len_Trim (Version(1))
        i = Len_Trim (DBtext(iHenkeCowan+2))
        write   (3,62)   MaxOrder,
     *                   Version(1)(1:j),
     *                   DBtext(iHenkeCowan+2)(1:i),
     *                   YesNo(iReqFields+1),
     *                   (YesNo(i_disableABC(i)+1),i=1,3)
  62    format  (
     *  'Maximum order of T1*..Tn:  10**k.....|',i3/
     *  'Algorithm type (Matrix / Recursive)..|',a/
     *  'Database used for dispersion correct.|',a/
     *  'Print fields.(*).....................|',a/
     *  '  (*) You cannot select printing and |'/
     *  '      transition layers at the same. |'/
     *  '      Also, if printing is ON,the rms|'/
     *  '      roughness is set zero.         |'/
     *  'Disabled A,B,C contributions to reso-|'/
     *  '     nance magnetic amplitudes.......|',3(a,1x))

        write   (3,11)
        write   (3,55)  ' *** User Comment: '
        do i=1,3  !==================================+
          j = Len_Trim (Comment(i))                 !|
          if (j.gt.0) write (3,55) Comment(i)(1:j)  !|
        enddo  !=====================================+
        write   (3,11)

        if (ireqFields .eq. 0) Then  !--------------------+
          if (Thickness_Tr(N_Total) .gt. 0.              !|
     *           .OR.    Trans_mean .gt. 0.) Then !-----+ |
                                                       !| |
            do  i=1,N_Tran  !=====================+     | |
              z = Real(i) / Real(N_Tran+1)       !|     | |
              if (Law_Index .eq. 0)     Then !-+  |     | |
c Cosine-law: f(z)=[1+cos(pi*z)]/2,  z:[0:1].  |  |     | |
c -> This is equivalent to: cos^2(pi*z/2):     |  |     | |
c f(z=0)=1, f(z=pi)=0.                         |  |     | |
                TranLaw(i) = cos(pi*z/2.)**2  !|  |     | |
              else  !--------------------------|  |     | |
c Linear-law: f(z)=1-z, z:[0:1].               |  |     | |
c f(z=0)=1, f(z=pi)=0.                         |  |     | |
                TranLaw(i) = 1.-z             !|  |     | |
              endif  !-------------------------+  |     | |
            enddo  !==============================+     | |
                                                       !| |
            write (3,66)  (TranLaw(i),i=1,N_Tran)      !v v
  66        format (/
     *      'Transition law (from upper to lower):'/
     *      101(10f7.3/:))                             !^ ^
          endif  !--------------------------------------+ |
        endif  !------------------------------------------+

        if (iDebug .ne. 0)
     *  Open    (unit=33,file=OutFile(1:lout)//'.deb',
     *                                  status='unknown')
c #######################################################
c              Print surface layer profile:
c#######################################################
        write   (3,61)
  61    format  (/'The Structure Profile:'/
     *           '====================='/
     *  '  Nr Thickness Roughness Transition ________Code________ ',
     *  'Fraction    rho      Exp(-w)',8x,'x0',10x,
     *  '|Resonance Magn.part    Magn.Density  Share of ',7x,
     *  'F10',16x,'F11',16x,'F1T',13x,'Magnetization         Identical'/
     *  16x,'rms (A)   layer (A)',32x,'(g/cm^3)   for x0',20x,
     *  '|  type    in x_ij(%)     (1/cm^3)   magn.atoms',59x,
     *  '    direction',11x,'to layer..'/)
c    *  1x,240('-'))
        do      lyr=1,N_Total   !============================+
          i = Magnetic(lyr+1)                               !|
          write (3,64)  lyr, Thickness(lyr),                !|
     *                  Sigma(lyr), Thickness_Tr(lyr),      !|
     *                  Code(1,lyr), Frac(1,lyr), rho(lyr), !|
     *                  W0(lyr), x0(lyr+1),                 !|
     *                  Mtype(i+1), xm_max(lyr+1),          !|
     *                  Mdensity(lyr+1), Mshare(lyr+1),     !|
     *                  F10(lyr+1), F11(lyr+1), F1T(lyr+1), !|
     *                  (Mvect(j,lyr+1),j=1,3),             !|
     *                  Identical(lyr)                      !|
          do j=2,N_Frac(lyr)  !=====================+        |
            write (3,65) Code(j,lyr), Frac(j,lyr)  !|        |
          enddo  !==================================+        |
        enddo  !=============================================+
        write (3,67)
  64    format  (i4, f10.2,
     *           f8.2, f11.2, 3x,
     *           a, f7.3, f10.4,
     *           f11.3, ' (',g8.2,',',g8.2,') | ',
     *           a, g12.3,
     *           g13.3, g13.3,
     *           3(' (',g8.2,',',g8.2,')'), 1x,
     *           '(',2(f6.2,','),f6.2,')',
     *           i6)
  65    format  (37x,a,f7.3)
  67    format  (1x,10('_')/
     *           ' Here Magnetic part in x_ij(%) is determined as:'/
     *           ' 100% * x_ij_magnetic_max/(x_0 + x_ij_magnetic_max)'/)

        write   (3,71) ((i,j,i=1,3),j=1,3)
  71    format  (/'Magnetic Structure Profile:'/
     *           '============================'/
c    *           '  Nr',7x,'x(F11)',14x,'x(F1T)',14x,'x(F10)',14x,
     *           '  Nr',5x,'x(F11+F1T)',10x,'x(F11-F1T)',7x,
     *                                      'x(2*F10-F11-F1T)',9x,
     *           9('x(',i1,','i1,')',14x))
        do      lyr=2,N_Total+1 !====================================+
c         A8 = coeff*Mdensity(lyr)*F11(lyr)                         !|
c         B8 = coeff*Mdensity(lyr)*F1T(lyr)                         !|
c         C8 = coeff*Mdensity(lyr)*F10(lyr)                         !|
          A8 = coeff*Mdensity(lyr)*(F11(lyr)+F1T(lyr))              !|
          B8 = coeff*Mdensity(lyr)*(F11(lyr)-F1T(lyr))              !|
          C8 = coeff*Mdensity(lyr)*(2.*F10(lyr)-F11(lyr)-F1T(lyr))  !|
          write(3,72) lyr-1, A8, B8, C8,                            !|
     *                ((x_ij(i,j,lyr),i=1,3),j=1,3)                 !|
        enddo  !=====================================================+
  72    format  (i4, 12(' (',g8.2,',',g8.2,')'))
        write (3,'(1x,a/)') '___________________________'

        do      lyr=1,N_Total-1   !==========+
          DDD = Thickness(lyr)              !|
     -        - Thickness_Tr(lyr)           !|
     -        - Thickness_Tr(lyr+1)         !|
                                            !|
          if (DDD .lt.   0.     .AND.       !|
     *        DDD .gt. -1.E-3)  DDD = 0.    !|
                                            !|
          if (DDD .lt. 0.) goto 103  !-------+---+
c If the wavefields are requested, then      |   |
c the roughness and the transition layers    |   |
c are ignored (no need to subtract!).        |   |
          if (iReqFields .eq. 0) Then !-+    |   |
            Thickness(lyr) = DDD       !|    |   |
          endif  !----------------------+    |   |
        enddo  !=============================+   v

c #######################################################
c                  Start of processing:
c#######################################################
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
  88    format (1x,i4,' points done (',f4.0,'%).  Incidence=',
     *  g11.4,1x,a,' Elapsed time=',a)
        if (modebat .ne. 0) Then  !-------------------+
          DoneFile = OutFile(1:lout)//'.~~~'         !|
          Call OpenFile(DoneFile(1:lout+4),50,       !|
     *              'write','unknown',io_status,*89) !|
          write (50,90,err=91)  nf, nf0,             !|
     *                          f0,Uni(iUnits+1),    !|
     *                          progress, timestamp  !|
  91      continue                                   !|
          close(unit=50)                             !|
        endif  !--------------------------------------+
  89    continue
  90    format ('Points done = ',i4,' of ',i4,' <br> ',
     *  'Current angle = ',g12.5,1x,a,' <br> ',
     *  'Progress = ',f10.2,'% <br> ',
     *  'Elapsed time (hh:mm:ss) =',a,' <br> ')
c#######################################################
c                   Loop over f0
c#######################################################

c NO message from DivZ!
c cc    Call    DivZ()
c YES message from DivZ!
        do      nf=1,nf0  !==============================+
          f0     = f0_range(1) + df0*(nf-1)             !|
c Calculate sin(f0):                                     |
          if (iUnits .ne. 5) Then  !-----+               |
c Scan over incidence angle:            !|               |
            f0_  = sin(f0*UnitCoef)     !|               |
          else  !------------------------|               |
c Scan over qz (in 1/Angstrom):         !|               |
            f0_  = f0*UnitCoef          !|sin(f0)=qz/2k  |
          endif  !-----------------------+               |
c         write (3,96)      char(12), nf,               !|
c    *                      f0, Uni(iUnits+1)           !|
c 96      format (1x,a,' point=',i4,' f0=',g11.4,1x,a)  !|
                                                        !|
c Call magnetic multilayer reflection engine:            |
c (some parameters are sent via Common/MAGdat/)          |
          Call  MMRE99 (f0, f0_, E0, Es, R0s,           !|
     *                  Phase_shift,                    !|
     *                  Circ_pm, Circ_Diff, Circ_Ratio, !|
     *                  iReqFields, Dx, Dy, Dz, u)      !|
                                                        !|
          if (ifail .ne. 0) Then  !-----------+          |
            write   (3,77)       ifail,nf    !v          v
  77        format  ('MMRE: Failure code ',
     *               i4,' -- at step: ',i4)  !^          ^
          endif  !----------------------------+          |
          if (iDebug .ne. 0) write (33,*) f0, N_Used    !|
                                                        !|
c Test for non-physical result:                         !|
          if (  R0s      .lt. 0.  .or.                  !|
     *          R0s      .gt. 1.  .or.                  !|
     *        Circ_pm(1) .lt. 0.  .or.                  !|
     *        Circ_pm(1) .gt. 1.  .or.                  !|
     *        Circ_pm(2) .lt. 0.  .or.                  !|
     *        Circ_pm(2) .gt. 1.) goto 53               !|
                                                        !|
          if (nf .eq. 1) write (3,'(1x,a)',err=161)     !|
     *                   'DAT file columns:'            !|
     *                 , 'f0 - incident angle or qz'    !|
c    *                 , 'f0*f0         '               !|
     *                 , 'R0s - reflectivity'           !|
c    *                 , 'Phase_shift   '               !|
c    *                 , 'Circ_Diff     '               !|
     *                 , 'Circ_pm(1) - e+ reflectivity' !|
     *                 , 'Circ_pm(2) - e- reflectivity' !|
     *                 , 'Circ_Ratio - e+/e- ratio'     !|
c    *                 , 'Abs(Es(1))**2 '               !|
c    *                 , 'Abs(Es(2))**2 '               !|
          write (1,'(1001g14.6)',err=161)               !|
     *                   f0                             !|
c    *                 , f0*f0                          !|
     *                 , R0s                            !|
c    *                 , Phase_shift                    !|
c    *                 , Circ_Diff                      !|
     *                 , Circ_pm(1)                     !|
     *                 , Circ_pm(2)                     !|
     *                 , Circ_Ratio                     !|
c    *                 , Abs(Es(1))**2                  !|
c    *                 , Abs(Es(2))**2                  !|
                                                        !|
          if (nf .eq. (nf/iBeepStep)*iBeepStep) Then !-+ |
            progress = 100.*Real(nf)/Real(nf0)        !| |
            if (progress .gt. 99.99) progress = 99.99 !| |
            Call  Duration (start,bias,hh,mm,ss)      !| |
            Call Convert_Time (hh,mm,ss,timestamp)    !| |
            if (iBatch_mode .lt. 2)                   !| |
     *        write (*,88)  nf, progress,             !| |
     *                      f0, Uni(iUnits+1),        !| |
     *                      timestamp                 !| |
              if (modebat .ne. 0) Then  !------------+ | |
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
                                                        !|
c If interrupted by <Esc> and confirmed:                !|
          if  (iasci .eq. 27) Then !-------------------+ |
                                                      !| |
            write (3,444) nf,nf0                      !| |
            goto 301  !->------------------------------+-+-+
                                                      !| | |
          elseif (FileExist(OutFile(1:lout)//'.esc')) !| | |
     *                                          Then !-| | v
            Call DelFile (OutFile(1:lout)//'.esc',i)  !| | |
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
        Close   (unit=1)
        Call    Duration (start,bias,hh,mm,ss)
        Call    Convert_Time (hh,mm,ss,timestamp)
        write   (3,63)  timestamp
  63    format  (/'   Elapsed time (hh:mm:ss)=',a)
        Close   (unit=3)
        if (modebat .ne. 0) Then !--------------------+
          Call OpenFile(DoneFile(1:lout+4),50,       !|
     *              'write','unknown',io_status,*95) !|
          f0       = f0_range(1) + df0*(nf-1)        !|
          progress = 100.*Real(nf)/Real(nf0)         !|
          write (50,90,err=94) nf, nf0,              !|
     *                         f0,Uni(iUnits+1),     !|
     *                         progress, timestamp   !|
  94      continue                                   !|
          close(unit=50)                             !|
        endif  !--------------------------------------+
  95    continue
        if (nf .lt. nf0)  goto 28

c cc    txt(1)=' *** Work done. Press any key to continue..'
        txt(1)=' *** Work done...'
        if (iBatch_mode .lt. 2) Then  !---+
          j = Len_Trim (txt(1))          !|
          write (*,'(1x,a)') txt(1)(1:j) !|
        endif  !--------------------------+
        do      i=1,1   !==================+
          Call  KeyIn   (0,iscan,iasci)   !|
          if (iscan .ne. 0)     goto 28   !|
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
     *  'at line:',i4,' -- while reading:'/)                    !+txt(7)
        ifail = 7
        goto 120

  101   continue
        write   (txt,119)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line
  119   format  (a,': Parameter(s) out of range in file:'//a//
     *  'at line:',i4,' -- while reading:'/)                    !+txt(7)
        ifail = 7
        goto 120

  113   continue
        write   (txt,213)  ProgramVer(1:ipv), InpFile(1:linp), line
  213   format  (a,': Request for the printout of wavefields in file:'//
     *  a//'at line:',i4,' is incompatible with transition layers')
        ifail = 5
        goto 120

  123   continue
        write   (txt,124)  ProgramVer(1:ipv), InpFile(1:linp)
  124   format  (a,': Zero substrate reflectivity in file:'//a//
     *  'Either a substrate code must be entered or x0#0.'/
     *  'Also make sure that W0>0.')
        ifail = 6
        goto 120

  131   continue
        write   (txt,132)  ProgramVer(1:ipv)
  132   format  (a,':'//
     *  'X-ray wavelength was chosen to be specified via'/
     *  'characteristic X-ray line, but no line given.')
        ifail = 4
        goto 130

  125   continue
        write   (txt,126)  ProgramVer(1:ipv), lyr, j
  126   format  (a,': No atomic density data for'//
     *  'layer =',i5,',  fraction =',i2,'.'//
     *  'Cannot calculate the magnetic atoms density'/
     *  'for this layer as a share of magnetic atoms'/
     *  'because of missing atomic density data.')
        ifail = 7
        goto 120

  53    continue
        write   ( 3, 136)  ProgramVer(1:ipv), f0, Uni(iUnits+1)
        write   (txt,136)  ProgramVer(1:ipv), f0, Uni(iUnits+1)
  136   format  (a,': algorithm internal consistency check failed'/
     *  'at incidence angle =',g12.5,1x,a/
     *  '*** Please, report the problem to the author!')
        ifail = 3
        goto 120

  109   continue
        write   (txt,209)  ProgramVer(1:ipv), InpFile(1:linp)
  209   format  (a,': Inconsistent parameters in file:'//a//
     *          'Multiple scan points with zero scan range'/
     *          'while reading:')                               !+txt(7)
        ifail = 7
        goto 120

  120   continue
        if (ifail .lt. 0) ifail = 0
        if (ifail .gt. 18) ifail=18
        txt(ifail+1) = ' '
        txt(ifail+2) = 'Check your input file!'
        ifail  = ifail+2
        goto 130

  701   continue
        qz_max = SQZ/Sigma(lyr)
        if (iUnits .ne. 5) Then !-------------+angular scan
          f0_max = sngl(qz_max/(2.*k0))      !|
          if (f0_max .lt. 1.) Then  !-------+ |
            f0_max = Asin(f0_max)/UnitCoef !| |
          else !----------------------------| |
            f0_max = (pi/2.)/UnitCoef      !| |
          endif  !--------------------------+ |
        else !--------------------------------|qz-scan
          f0_max = qz_max                    !|
        endif  !------------------------------+

        write   (txt,702)  ProgramVer(1:ipv),
     *                     Sigma(lyr), lyr,
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
     *  'f0_max =',g9.3,1x,a)
        ifail = 11
        goto 130

  103   continue
        write   (txt,203)  ProgramVer(1:ipv),lyr,
     *                     Thickness(lyr),
     -                     Thickness_Tr(lyr),
     -                     Thickness_Tr(lyr+1),
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
        Close   (unit=1,status='delete',err=130)
        goto 130

  137   continue
        write   (txt,138)  ProgramVer(1:ipv), Energy, Energy_max
  138   format  (a,': Requested X-ray Energy = ',g12.5,' keV'/
     *  'exceeds current threshold E_max =',f6.2,' keV for this',1x,
     *  'algorithm.'//
     *  'This algorithm may fail at high energies due to precision'/
     *  'losses when trying to account for small refraction or X-rays.'/
     *  /
     *  'Please use a complimentary algorithm specifically designed'/
     *  'for mid- and high-energy x-rays (6keV and up).')
        ifail = 8
        goto 130

  130   continue
        if (modebat .eq. 1) Then  !--------------+
          do    i=1,Iabs(ifail)  !=============+ |
            j = Max0(Len_Trim(txt(i)),1)      !| |
            write (*,'(1x,a)') txt(i)(1:j)    !| |
          enddo  !=============================+ |
        endif  !---------------------------------+
        Call    Message (txt,Iabs(ifail),2)
        goto 28

        End

c=====================================================

        Integer Function P01AAF (ifail, error, srname)

c Returns the value of NAG-2 error or terminates the program.
c ---------------------------------------------------------
        Integer error, ifail, nout
        Double precision srname

c Test if no error detected:

        if (error .eq. 0) Goto 20

c Determine output unit for message

        Call X04AAF (0,nout)

c Test for soft failure:

        if (mod(ifail,10) .eq. 1) Goto 10

c Hard failure:
c -------------
        Write (nout,99999) srname, error

c Stopping mechanism may also differ:
c ccc   Call exit_quiet()
        Goto 20

c Soft failure:
c -------------
c Test if error messages suppressed:
10      Continue
        if (mod(ifail/10,10) .eq. 0) Goto 20
        Write (nout,99999) srname, error

20      continue
        P01AAF = error
        Return

!99999  Format (/' Error detected by NAG library routine ',
!     *           a8, ' - ifail = ', i5/)
99999   Format (/' Error from ',a8, ' - ifail = ', i5/)
        End

c=====================================================

      Subroutine X04AAF(i,nerr)
      Integer i, nerr
      Integer NERR1
      Data nerr1 /3/
      if (i .eq. 0) nerr = nerr1
      if (i .eq. 1) nerr1 = nerr
      Return
      End


