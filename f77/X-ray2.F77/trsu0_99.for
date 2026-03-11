c
c *** Version 99a ***  for trds_99
c     -- Interface roughness profile
c        (in dependence on the interface number)
c     -- 10 models of roughness correlations
c
c +=====================================================+
c |     These routines implement input for trds_10c     |
c |-----------------------------------------------------|
c |Contents:                                            |
c | 1. Subroutine:  Read_TRDS_Input                     |
c | 2. Subroutine:  Read_Top99_0  => red99prf0.for      |
c | 3. Subroutine:  Read_Sub99_0  => red99prf0.for      |
c | 4. Subroutine:  Use_Symmetry                        |
c |-----------------------------------------------------|
c |Subroutine calls: see  --  TRSU1_97.for,             |
c |                       --  TRSU2_97.for              |
c +=====================================================+

c ####################################################### 1

        Subroutine Read_TRDS_Input (InpFile, OutFile, File_Top,
     *                          ProgramVer, ipv, N_Points_Max,
     *                          iUnits, iBeepStep, mode_scan,
     *                          iBatch_mode, GRD_output, iskip,
     *                          icont, jcont, Thickness_Top,
     *                          Scan_Range, Scan_Step, nScan,
     *                          Code, Frac, N_Frac, N_Frac_Max,
     *                          rho, W0, C00, Identical,
     *                          ipol, Radiat)

c--------------------------------------------------------
c          Reading input file for TRDS, ver.99a
c--------------------------------------------------------
c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c ccc   Parameter      (N_Top_Max     = 150,
c ccc*                  N_Total_Max   = N_Top_Max+1)
        Integer         N_Frac_Max,
     *                  N_Points_Max

        Complex*8       x0_, x0_f

        Real*4          Frac(N_Frac_Max, N_Total_Max),
     *                  rho(N_Total_Max),
     *                  W0(N_Total_Max),
     *                  Scan_Range(2,2), Scan_Rng,
     *                  Scan_Step(2),
     *                  UnitData(5),
     *                  S,                            !==S_factor
     *                  mr,                           !==m/m0
     *                  AngDegree, AngMinute,
     *                  AngSec, AngMrad, AngMurad,
     *                  Thickness_Top, x,
     *                  xr0, xi0, rho_, rho_f

        Real            wave2energy
        External        wave2energy

        Integer         nScan(2),
     *                  mode_scan, iUnits(3),
     *                  N_Frac(N_Total_Max),
     *                  Identical(N_Total_Max),
     *                  i, l, j, m, m0, icont,
     *                  jcont, linp, ipv, line,
     *                  iBatch_mode, lout, ltop,
     *                  ixway, iPol, iskip,
     *                  i_VR, io_status,
     *                  iBeepStep, mode_x0,
     *                  iSyngony, mode_x0h,
     *                  iTmp(1)

        Character       Code(N_Frac_Max,N_Total_Max)*20,
     *                  C00(N_Total_Max)*2,
     *                  Radiat*6, wrk*80,
     *                  InpFile*80, OutFile*80,
     *                  File_Top*80,
     *                  ProgramVer*8

        Logical         YesTr, GRD_output

c--------------------------------------------------------
        Complex*8       x0(N_Total_Max+1)
        Real*8          Z(N_Total_Max),         ! coordinates of interfaces
     *                  k0                      ! k0 = 2*pi/wave
        Real*4          Wave, UnitCoef(3), qz, qx,
     *                  Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Sigma_Accu(N_Total_Max),
     *                  CorreLength, CorreVert, h_jagged,
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

c New of August'99:
c       Parameter       (kcompMax = 10)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Real            prcn(kcompMax)
        Common  /x0pa3/ prcn

        Character       txt(20)*80
        Common  /msg/   txt
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
          progname = 'TRSU0_99'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+

c -------------------------------
        i = N_Top_Max                   !to suppress G77 warnings
        i = N_Total_Max                 !to suppress G77 warnings

        i = kcompMax                    !to suppress G77 warnings
        i = kolMax                      !to suppress G77 warnings
        i = ncodMax                     !to suppress G77 warnings
        i = natmMax                     !to suppress G77 warnings
        i = nwavMax                     !to suppress G77 warnings
        i = nedgMax                     !to suppress G77 warnings
        i = nabcMax                     !to suppress G77 warnings
        i = nedgTb                      !to suppress G77 warnings
        i = maxX0hWarnLines             !to suppress G77 warnings
c ===========================================================
        AngDegree   = 2.*pi/360.
        AngMinute   = AngDegree/60.
        AngSec      = AngDegree/3600.
        AngMrad     = 1.E-03
        AngMurad    = 1.E-06

        UnitData(1) = AngDegree
        UnitData(2) = AngMinute
        UnitData(3) = AngMrad
        UnitData(4) = AngSec
        UnitData(5) = AngMurad

        do      l=1,N_Total_Max  !====+
          x0(l)        = (0.,0.)     !|
          Thickness(l) = 0.          !|
          Z(l)         = 0.          !|
          Sigma(l)     = 0.          !|
          Sigma_Accu(l)= 0.          !|
          W0(i)        = 0.          !|
          rho(i)       = 0.          !|
          N_Frac(l)    = 1           !|
          Identical(l) = 0           !|
          C00(l)       = ' '         !|
          do    j=1,N_Frac_Max  !==+  |
            Frac(j,l)  = 0.       !|  |
            Code(j,l)  = ' '      !|  |
          enddo  !=================+  |
          Frac(1,l) = 1.             !|
        enddo  !======================+
        x0(N_Total_Max+1) = (0.,0.)

        modebat = 1     !assume batch-mode by default
c #######################################################
c #                    Read input file                  #
c #######################################################
        InpFile = ' '
c Analyze command-line arguments (See Serve.lib):
c [Here icont(1) is an integer array of length=1
c returning the number specified in the command
c line with the key '/C'. If '/C' is missing,
c then GetInput returns icont=-2.  If '/C' is
c present, the  GetInput returns icont=-1. The
c other situations (/C<icont>) are not expected.]
        Call    GetInput (InpFile,icont,'C',1,ifail)
        if (ifail.ne.0) goto 28
        if (icont.ne.-1 .and. icont.ne.-2) Then  !------------------+
          ifail = 5                                                !|
          txt(1) = "GetInput: /C flag must not have any argument"  !|
          Call  Message (txt,1,2)                                  !|
          goto 28                                                  !|
        endif  !----------------------------------------------------+

        if (icont.eq.-2) Then !-+
          icont = 0            !|Key: '/C' is not specified
        else  !-----------------|
          icont = 1            !|Key: '/C' is specified
        endif  !----------------+
        jcont = 0

        linp = Len_Trim(InpFile)
        if (linp.eq.0)  Then  !-----------------+
          InpFile = ProgramVer(1:ipv)//'.inp'  !|
          linp = Len_Trim(InpFile)             !|
        endif  !--------------------------------+

        Call OpenFile(InpFile(1:linp),2,'read','old',io_status,*99)

        Call Make_Err_Filename (InpFile(1:linp))

        line = 0
        txt(7) = 'batch mode flag (range=[0:2])'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        iBatch_mode = iTmp(1)
        if (iBatch_mode.lt.0)   goto 101
        if (iBatch_mode.gt.2)   goto 101
        if (iBatch_mode.eq.0) iBatch_mode=1
        modebat = iBatch_mode
c-----------------------
        txt(7) = 'output file name'
        Call    LineLoop (2,line,wrk,*100)
        OutFile = wrk
        lout = Len_Trim(OutFile)
        if (lout.eq.0)  Then  !-----+
          OutFile = '~'            !|
          lout    = 1              !|
        endif  !--------------------+
c       Call    Case    (OutFile(1:lout),1)     !to lower case
  116   continue   !<---------------------------------------------+
        if (icont.eq.0) Then  !-----------------------------+     |
          Call  OpenList(OutFile(1:lout)//'.tbl',3,ifail)  !|     |
          if (ifail.ne.0) goto 28                          !|     |
        else  !---------------------------------------------|     |
          Open (unit=3,file=OutFile(1:lout)//'.tbl',       !|     |
     *          status='old',err=121) !--+no file!          |     |
          goto 122  !-file present!------+---------+        |     |
  121     continue  !<-------------------+        !|        |     |
c No TBL-file present in the continuation mode.    |        |     |
c Switching to the normal mode!                    |        |     |
          icont = 0                               !|        |     |
          goto 116  !------------------------------+--------+-----+
  122     continue  !<-----------------------------+        |
c Append:                                                   |
          l = 1                                            !|
  25      continue      !<===========+                      |
          read  (3,'(a)',end=9) !--+ |                      |
          l = l + 1               !V |                      |
          goto 25  !=======>=======+=+                      |
  9       continue  !<-------------+                        |
          if (l.gt.1)  Backspace (3,err=262)               !|
  262     continue                                         !|
        endif  !--------------------------------------------+
c-----------------------
        txt(7)   = 'surface profile name'
        N_Top    = 0
        N_Total  = 0
        Call    LineLoop (2,line,wrk,*100)
        File_Top = wrk
        ltop     = Max (Len_Trim(File_Top),1)
c       Call    Case    (File_Top(1:ltop),1)    !to lower case
        YesTr    = .False.
        txt(7)   = 'surface layer profile'
c See: Rd99prf0.for:
        Call    Read_Top99_0 (File_Top(1:ltop),
     *                        N_Top, N_Top_Max,
     *                        Thickness, Thickness_Top,
     *                        Code(1,1), Frac(1,1),
     *                        N_Frac(1), N_Frac_Max,
     *                        rho(1), C00(1), x0(2), W0(1),
     *                        Identical(1),
     *                        Sigma(1),
     *                        Sigma(1), YesTr,
     *                        ifail)
        if (ifail.ne.0) goto 28
        N_Total = N_Top+1
c-----------------------
        do i=1,3  !=====================================+
          write (txt(7),'(a,i1)') 'comment line #',i   !|
          Call  LineLoop (2,line,Comment(i),*100)      !|
        enddo  !========================================+
c-----------------------
        txt(7) = 'substrate data'
        Thickness(N_Total) = 0.
        W0(N_Total)        = 1.
        YesTr = .False.
c See: Rd99prf0.for:
        Call Read_Sub99_0 (Code(1,N_Total), rho(N_Total),
     *                     C00(N_Total),
     *                     x0(N_Total+1), W0(N_Total),
     *                     Sigma(N_Total),
     *                     Sigma(N_Total), YesTr,
     *                     wrk, line, *100, *101, *28)
c-----------------------
        txt(7) = 'x-ray specification mode (range=[1:3])'
        Call    LineLoop (2,line,wrk,*100)
c 1=wave 2=energy 3=x-ray line:
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        ixway = iTmp(1)
        if (ixway.lt.1)         goto 101
        if (ixway.gt.3)         goto 101
c-----------------------
        txt(7) = 'x-ray wavelength/energy (range=[0.01-1000.0]Angstrom)'
        Call    LineLoop (2,line,wrk,*100)
        Wave = 0.
        if (ixway.lt.3) Then  !-------------------+
          Call  rdReal  (Wave,1,wrk,i)           !|
          if (i.ne.0)           goto 100         !|
          if (Wave .le. 0.)     goto 101         !|
          if (ixway.eq.2) Wave=wave2energy(Wave) !|
          if (Wave .lt. 0.01)   goto 101         !|
          if (Wave .gt. 1000.)  goto 101         !|
        endif  !----------------------------------+
c-----------------------
        txt(7) = 'x-ray line name'
        Call    LineLoop (2,line,wrk,*100)
        if (ixway.eq.3) Then  !-------------------+
          Radiat = wrk(1:6)                      !|
          if (Len_Trim(Radiat).eq.0) goto 131 !---+---+
c Get wavelength from the X-ray line name:        |   v
          Call  Get_Wave (Wave,Radiat,ifail)     !|
          if (ifail.ne.0)       goto 130         !|
          if (Wave .lt. 0.01)   goto 101         !|
          if (Wave .gt. 1000.)  goto 101         !|
        else !------------------------------------|
          Radiat = ' '                           !|
        endif  !----------------------------------+
c-----------------------
        txt(7) = 'x-ray polarization [1=sigma 2=pi]'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 130
        iPol = iTmp(1)
c Polarization: 0-mixed, 1-sigma, 2-pi:
        if (iPol.eq.0)  iPol=1
        if (iPol.lt.1)  goto 101
        if (iPol.gt.2)  goto 101
c-----------------------
        txt(7) = 'angular units (range=[0:4])'
        Call    LineLoop (2,line,wrk,*100)
c Read angular units (degrees, minutes, seconds, etc):
        Call    rdInt   (iUnits,3,wrk,i)
        if (i.ne.0)             goto 100
c iUnit(1) is for scan angles
c iUnit(2) is for skew roughness transfer
c iUnit(3) is for surface miscut angle
        do i=1,3  !=============================+
          if (iUnits(i).lt.0)   goto 101       !|
          if (iUnits(i).gt.4)   goto 101       !|
          UnitCoef(i) = UnitData(iUnits(i)+1)  !|
        enddo  !================================+
c-----------------------
        txt(7) = 'scan mode (range=[1:4])'
        Call    LineLoop (2,line,wrk,*100)
c Read mode of scan:
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        mode_scan = iTmp(1)
c mode_scan=1: over Q  at 2Q=fix
c mode_scan=2: Q-2Q at fixed offsets
c mode_scan=3: over 2Q  at Q=fix
c mode_scan=4: qx scans at fixed qz (1/A)
        if (mode_scan.lt.1)     goto 101
        if (mode_scan.gt.4)     goto 101
c-----------------------
        txt(7) = 'scan limits'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (Scan_Range(1,1),2,wrk,i)
        if (i.ne.0)             goto 100
        if     (mode_scan.eq.1) Then  !-----------------------------------+
c 1- Q-scans at 2Q=fixed              - Read incidence angle range:      !|
c Theta_In  = Theta_Scan                                                  |
c Theta_Out = Theta_Off - Theta_Scan                                      |
          txt(7) = 'scan limits (range=[0:20]degr.)'                     !|
          if (Scan_Range(1,1).lt.0.)  goto 101                           !|
          if (Scan_Range(2,1).lt.0.)  goto 101                           !|
          if (UnitCoef(1)*Scan_Range(1,1).gt.20.*AngDegree) goto 101     !|
          if (UnitCoef(1)*Scan_Range(2,1).gt.20.*AngDegree) goto 101     !|
        elseif (mode_scan.eq.2) Then  !-----------------------------------|
c 2- Q-2Q scans at fixed Q-offsets    - Read incidence angle range:      !|
c Theta_In  = Theta_Scan + Theta_Off                                      |
c Theta_Out = Theta_Scan - Theta_Off                                      |
          txt(7) = 'scan limits (range=[0:10]degr.)'                     !|
          if (Scan_Range(1,1).lt.0.)  goto 101                           !|
          if (Scan_Range(2,1).lt.0.)  goto 101                           !|
          if (UnitCoef(1)*Scan_Range(1,1).gt.10.*AngDegree) goto 101     !|
          if (UnitCoef(1)*Scan_Range(2,1).gt.10.*AngDegree) goto 101     !|
        elseif (mode_scan.eq.3) Then  !-----------------------------------|
c 3- 2Q (detector) scans at Q=fixed   - Read exit angle range:           !|
c Theta_In  = Theta_Off                                                   |
c Theta_Out = Theta_Scan - Theta_Off                                      |
          txt(7) = 'scan limits (range=[0:20]degr.)'                     !|
          if (Scan_Range(1,1).lt.0.)  goto 101                           !|
          if (Scan_Range(2,1).lt.0.)  goto 101                           !|
          if (UnitCoef(1)*Scan_Range(1,1).gt.20.*AngDegree) goto 101     !|
          if (UnitCoef(1)*Scan_Range(2,1).gt.20.*AngDegree) goto 101     !|
        elseif (mode_scan.eq.4) Then  !-----------------------------------|
c 4- qx scans at qz=fixed             - Read qx range:                   !|
          txt(7) = 'scan limits (range=[-1:1]A^-1)'                      !|
          if (Abs(Scan_Range(1,1)).gt.1.) goto 101                       !|
          if (Abs(Scan_Range(2,1)).gt.1.) goto 101                       !|
        endif  !----------------------------------------------------------+
c-----------------------
        write (txt(7),456) 'scan points', N_Points_Max
  456   format('number of ',a,' (min=1, max=',i4,')')
        Call    LineLoop (2,line,wrk,*100)
c Read the number of points on the scan:
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        nscan(1) = iTmp(1)
        if (nscan(1).lt.1)      nscan(1)=1
        if (nscan(1).gt.N_Points_Max)  goto 101
        if (nscan(1).gt.1)      Then  !-----------------------+
          Scan_Rng     = Scan_Range(2,1) - Scan_Range(1,1)   !|
          Scan_Step(1) = Scan_Rng / (nscan(1)-1)             !|
          if (abs(Scan_Rng).lt.1.E-32) goto 102              !|
        else  !-----------------------------------------------+
          Scan_Step(1)    = 0.                               !|
          Scan_Range(2,1) = Scan_Range(1,1)                  !| 2006/08/01
        endif  !----------------------------------------------+
c-----------------------
        txt(7) = 'Offset Limits'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (Scan_Range(1,2),2,wrk,i)
        if (i.ne.0)             goto 100
        if     (mode_scan.eq.1) Then  !-----------------------------------+
c 1- Q-scans at 2Q=fixed              -Read 2Q-angle range:               |
c Theta_In  = Theta_Scan                                                  |
c Theta_Out = Theta_Off - Theta_Scan                                      |
          txt(7) = '2Q limits (range=[0:20]degr.)'                       !|
          if(Scan_Range(1,2).lt.0.)  goto 101                            !|
          if(Scan_Range(2,2).lt.0.)  goto 101                            !|
          if(UnitCoef(1)*Scan_Range(1,2).gt.20.*AngDegree) goto 101      !|
          if(UnitCoef(1)*Scan_Range(2,2).gt.20.*AngDegree) goto 101      !|
        elseif (mode_scan.eq.2) Then  !-----------------------------------|
c 2- Q-2Q scans at fixed Q-offsets    -Read Q-offset range:               |
c Theta_In  = Theta_Scan + Theta_Off                                      |
c Theta_Out = Theta_Scan - Theta_Off                                      |
          txt(7) = 'Q-offset limits (range=[-10:10]degr.)'               !|
          if(Abs(UnitCoef(1)*Scan_Range(1,2)).gt.10.*AngDegree) goto 101 !|
          if(Abs(UnitCoef(1)*Scan_Range(2,2)).gt.10.*AngDegree) goto 101 !|
        elseif (mode_scan.eq.3) Then  !-----------------------------------|
c 3- 2Q (detector) scans at Q=fixed   -Read Q-angle range:                |
c Theta_In  = Theta_Off                                                   |
c Theta_Out = Theta_Scan - Theta_Off                                      |
          txt(7) = 'Q-fixed limits (range=[0:20]degr.)'                  !|
          if(Scan_Range(1,2).lt.0.)  goto 101                            !|
          if(Scan_Range(2,2).lt.0.)  goto 101                            !|
          if(UnitCoef(1)*Scan_Range(1,2).gt.20.*AngDegree) goto 101      !|
          if(UnitCoef(1)*Scan_Range(2,2).gt.20.*AngDegree) goto 101      !|
        elseif (mode_scan.eq.4) Then  !-----------------------------------|
c 4- qx scans at qz=fixed             -Read qz range:                     |
          txt(7) = 'qz limits (range=[0:1]A^-1)'                         !|
          if(Scan_Range(1,2).lt.0.)  goto 101                            !|
          if(Scan_Range(2,2).lt.0.)  goto 101                            !|
          if (Scan_Range(1,2).gt.1.) goto 101                            !|
          if (Scan_Range(2,2).gt.1.) goto 101                            !|
        endif  ! ---------------------------------------------------------+
c-----------------------
        write (txt(7),456) 'offsets', N_Points_Max
        Call    LineLoop (2,line,wrk,*100)
c Read the number of offsets:
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        nscan(2) = iTmp(1)
        if (nscan(2).lt.1)      nscan(2)=1
        if (nscan(2).gt.N_Points_Max)  goto 101
        if (nscan(2).gt.1)      Then  !-----------------------+
          Scan_Rng     = Scan_Range(2,2) - Scan_Range(1,2)   !|
          Scan_Step(2) = Scan_Rng / (nscan(2)-1)             !|
          if (abs(Scan_Rng).lt.1.E-32) goto 102              !|
        else  !-----------------------------------------------|
          Scan_Step(2)    = 0.                               !|
          Scan_Range(2,2) = Scan_Range(1,2)                  !| 2006/08/01
        endif  !----------------------------------------------+
c-----------------------
c GRD Output:
        if (nscan(1).gt.1 .AND. nscan(2).gt.1) Then  !-+
          GRD_Output = .True.                         !|
        else   !---------------------------------------|
          GRD_Output = .False.                        !|
        endif  !---------------------------------------+
c-----------------------
        txt(7) = '"skip negative angles" flag (1/0)'
        Call    LineLoop (2,line,wrk,*100)
c Skip points corresponding to negative incidence/exit angles?
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        iskip = iTmp(1)
c iskip=0  - No, don't skip!
c iskip=1  - Yes, skip!
        if (iskip.lt.0)         goto 101
        if (iskip.gt.1)         goto 101
c ERROR: skipping is not compatible with continuation mode:
        if (iskip.eq.1 .and.icont.eq.1) Then  !-----------+
          ifail = 101                                    !|
          txt(10)='Skipping negative angles is incom'//  !|
     *            'patible with continuation mode'       !|
          goto 100                                       !|
        endif  !------------------------------------------+
c-----------------------
        txt(7) = 'specular/diffuse flag (1/0)'
        Call    LineLoop (2,line,wrk,*100)
c What compute at specular?
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        i_Spec = iTmp(1)
c i_Spec=1  - Reflection coefficient
c i_Spec=0  - Diffuse scattering
        if (i_Spec.lt.0)        goto 101
        if (i_Spec.gt.1)        goto 101

c-----------------------ACCELERATORS:
        txt(7) = 'Born accelerator flag (0/1/2)'
        Call    LineLoop (2,line,wrk,*100)
c Read calculation method:
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        i_Born = iTmp(1)
c 0 = Distorted wave Born approximation
c 1 = semi-Born approximation, R=0 (no reflected waves)
c 2 = Born approximation, R=0, T=1 (no change in transmitted waves)
        if (i_Born.lt.0)        goto 101
        if (i_Born.gt.2)        goto 101
c-----------------------
        txt(7) = 'Fourier accelerator flag (1/0)'
        Call    LineLoop (2,line,wrk,*100)
c Read the calculation method for the correlation
c integral at large s=qz*sigma:
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        i_Fourier = iTmp(1)

c                              s*K(x)
c i_Fourier=0:    Integral{dx(e      - 1)}
c
c i_Fourier=1:    s*Integral{dx K(x)}

        if (i_Fourier.lt.0)     goto 101
        if (i_Fourier.gt.1)     goto 101

c-----------------------ROUGHNESS:
        txt(7) = 'roughness correlation model [1-20] '
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        Index_corfu = iTmp(1)
c Read index of correlation function type:
c  1 - No vertical correlation
c  2 - Complete vertical correlation
c  3 - Holy & Baumbach`s model
c  4 - Ming,Krol,Soo,Kao,Park,Wang model
c  5 - Spiller and Sterns model
c  6 - Phang,Savage,Kariotis,Lagally model
c  7 - Pukite: monoatomic steps [exact]
c  8 - Pukite: steps+spread [exact]
c  9 - Pukite/Pershan: steps+spread
c 10 - Ming + periodic relief period
        if (Index_corfu.lt.1)   goto 101
        if (Index_corfu.gt.20)  goto 101
        if (Index_corfu.ge.7 .AND.
     *      Index_corfu.le.9)   Then  !---+
          i_Born    = 2                  !|     !for Pukite use Born approximation
c         i_Born    = 1                  !|     !for Pukite use semi-Born approximation
c         i_Born    = 0                  !|     !for Pukite use DBWA approximation
          i_Fourier = 1                  !|     !for Pukite
        endif  !--------------------------+
        if (Index_corfu.eq.10)  Then  !---+
          i_Fourier = 1                  !|     !small roughness approx.
        endif  !--------------------------+
c-----------------------
        txt(7) = 'lateral correlation length (L_h > 0)'
        Call    LineLoop (2,line,wrk,*100)
c Read roughness lateral correllation length (Angstrem)
        Call    rdReal  (CorreLength,1,wrk,i)
        if (i.ne.0)               goto 100
        if (CorreLength.le.0.)    goto 101
c       if (CorreLength.gt.1.E+5) goto 101
c-----------------------
        txt(7) = 'vertical correlation length (L_v > 0)'
        Call    LineLoop (2,line,wrk,*100)
c Read roughness vertical correllation length (Angstrem)
        Call    rdReal  (CorreVert,1,wrk,i)
        if (i.ne.0)             goto 100
        if (CorreVert.lt.0.)    goto 101
        if (CorreVert.gt.1.E+10) CorreVert=1.E+10       ! 1m

        if (Index_corfu.eq.1)   CorreVert=0.         ! No
        if (Index_corfu.eq.2)   CorreVert=1.E+10     ! 1m (inf.) - for complete
        if (Index_corfu.eq.5)   CorreVert=1.E+10     ! 1m (inf.) - for Holy
c For Ming or Lagally:
        if (Index_corfu.eq.3 .OR. Index_corfu.eq.4) Then !-+
          if (abs(CorreVert).lt.1.E-32) Index_corfu=1     !|->(no vertical)
          if (CorreVert.ge.1.E+10)      Index_corfu=2     !|->(complete)
        endif !--------------------------------------------+
c-----------------------
        txt(7) = 'jaggedness (range=[0.1:1])'
        Call    LineLoop (2,line,wrk,*100)
c Read roughness jaggedness (0.1-1.0)
        if (Len_Trim(wrk).eq.0) Then  !-------+
          h_jagged = 1.0                     !| not entered
        else  !-------------------------------|
          Call  rdReal  (h_jagged,1,wrk,i)   !|
          if (i.ne.0)           goto 100     !|
        endif  !------------------------------+
c cccc  if (Index_corfu.eq.4)   h_jagged=0.5         ! Phang,..Lagally
        if (Index_corfu.ge.7 .AND.                   ! Pukite ...
     *      Index_corfu.le.9)   h_jagged=0.5         ! ... Pukite
c cccc  if (Index_corfu.eq.10)  h_jagged=1.0         ! Gaussian for periodic
        if (h_jagged .lt. 0.1)  goto 101
        if (h_jagged .gt. 1.0)  goto 101
c-----------------------
        txt(7) = 'inheritage inclination angle (range=[-89:89]'
        Call    LineLoop (2,line,wrk,*100)
c Read roughness interface-interface inheritage inclination angle:
c (can be used for all models, except for Spiller's one)
        Call    rdReal  (RouMiscut,1,wrk,i)
        if (i.ne.0)              goto 100
        if (RouMiscut .lt.-89.)  goto 101
        if (RouMiscut .gt. 89.)  goto 101
c-----------------------
        txt(7) = 'Lagally''s cross correlation length (L_cross >= L_h)'
        Call    LineLoop (2,line,wrk,*100)
c Read the lateral correlation length for inter-interface
c roughness correlation (Lagally model):
        Call    rdReal   (CorreCross,1,wrk,i)
        if (i.ne.0)             goto 100
        if (Index_corfu.eq.4)   Then  !--------------+
          if (CorreCross.lt.CorreLength)  goto 101  !|
        endif  !-------------------------------------+
c-----------------------
        txt(7) = 'Spiller''s relaxation length (L_relax > 0)'
        Call    LineLoop (2,line,wrk,*100)
c Read diffusion-like relaxation length for the Spillers model:
        Call    rdReal   (RelaxLength,1,wrk,i)
        if (i.ne.0)             goto 100
c This parameter controls the interface-interface roughness
c correlation. It's connection with the vertical correlation
c length is approximately as follows:
c +----------------------------------------+
c |CorreVert = CorreLength**2 / RelaxLength|
c +----------------------------------------+
        if (Index_corfu.eq.6)   Then  !-----------+
            if (RelaxLength.lt.0.)    goto 101   !|
        endif  !----------------------------------+
c-----------------------
        txt(7) = 'Spiller''s relaxation/correlation flag [0-1]'
        Call    LineLoop (2,line,wrk,*100)
c What take for Spiller`s model: CorreVert or RelaxLength?
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        i_VR = iTmp(1)
c i_VR=0  - Vertical Correlation Length
c i_VR=1  - Relaxation Length
        if (Index_corfu.eq.6)  Then  !--------------------+Spiller's model
          if (i_VR.gt.1)        goto 101                 !|
          if (i_VR.lt.0)        goto 101                 !|
        if (i_VR.eq.0)  Then  !--------------------------+|<-CorreVert input
          if (abs(CorreVert).lt.1.E-32) Then  !--------+ ||
            Index_corfu=1                             !| ||no correlation
          elseif (CorreVert.ge.1.E+10)  Then  !--------+ ||
            Index_corfu=5                             !| ||Holy
          else !---------------------------------------+ ||
            RelaxLength = CorreLength**2 / CorreVert  !| ||
          endif  !-------------------------------------+ ||
        else  !------------------------------------------||<-RelaxLength input
          if (abs(RelaxLength).lt.1.E-32) Then !-------+ ||
             Index_corfu=5                            !| ||Holy
          else  !--------------------------------------+ ||
            CorreVert = CorreLength**2 / RelaxLength  !| ||
            if (CorreVert.le.1.E-3)     Index_corfu=1 !| ||no corr.(<0.001A)
            if (CorreVert.ge.1.E+10)    Index_corfu=5 !| ||Holy
          endif  !-------------------------------------+ ||
        endif  !-----------------------------------------+|
        endif  !------------------------------------------+

        if (Index_corfu.eq.1) Then  !---+
          CorreVert   = 0.             !|               ! 1m (if no   correl.)
          RelaxLength = 1.E+10         !|               ! 1m (if no   correl.)
        endif  !------------------------+
        if (Index_corfu.eq.2) Then  !---+
          CorreVert   = 1.E+10         !|               ! 1m (if no   correl.)
          RelaxLength = 0.             !|               ! (if Holy model)
        endif  !------------------------+
        if (Index_corfu.eq.6) Then  !---+
          h_jagged=1.0                 !|               ! Spiller
          RouMiscut=0.                 !|               ! not for Spiller!
        endif  !------------------------+
c-----------------------
        txt(7) = 'integral precision'
        EpsRel = 1.E-3
c cccc  Call    LineLoop (2,line,wrk,*100)
c Read the relative precision of correlation function integral
c (Spiller model):
c cccc  Call    rdReal  (EpsRel,1,wrk,i)
c cccc  if (i.ne.0)             goto 100
        if (EpsRel.lt.(1.E-4))  EpsRel = 1.E-4
        if (EpsRel.gt.(1.E-1))  EpsRel = 1.E-1
c-----------------------
        txt(7) = 'Pukite''s surface miscut angle (range=[-20:20])'
        Call    LineLoop (2,line,wrk,*100)
c Read surface miscut angle (Pukite model):
        Call    rdReal  (SurMiscut,1,wrk,i)
        if (i.ne.0)             goto 100
        if (RouMiscut .lt.-20.)  goto 101
        if (RouMiscut .gt. 20.)  goto 101
        if (Index_corfu.lt.7  .OR.
     *      Index_corfu.gt.9) Then  !--+
          SurMiscut=0.                !|
        endif  !-----------------------+
c-----------------------
        txt(7) = 'Pukite''s rms steps height (h >= 0)'
        Call    LineLoop (2,line,wrk,*100)
c Read steps rms height (Angstrem, Pukite model):
        Call    rdReal  (StepHeight,1,wrk,i)
        if (i.ne.0)             goto 100
        if (Index_corfu.lt.7  .OR.
     *      Index_corfu.gt.9) Then  !--+
          StepHeight=0.               !|
        endif  !-----------------------+
        if (StepHeight.lt.0.)   goto 101
c-----------------------
        txt(7) = 'Pershan''s terrase spread (L_h > Spread > 0)'
        Call    LineLoop (2,line,wrk,*100)
c Read terrace width spread (Angstrem, Pukite/Pershan model):
        Call    rdReal  (TerraceSpre,1,wrk,i)
        if (i.ne.0)             goto 100
        if (Index_corfu.eq.9) Then !-----------------+
          if (TerraceSpre.le.0.)          goto 101  !|2006/08/01
          if (TerraceSpre.ge.CorreLength) goto 101  !|
        else  !--------------------------------------|
          TerraceSpre = 0.                          !|
        endif  !-------------------------------------+
c-----------------------
        txt(7) = 'Pukite+Sinha summation flag (1/0)'
        Call    LineLoop (2,line,wrk,*100)
c Read the flag requesting the summation of scattering from
c steps and self-affine roughness:
        Call    rdInt   (iTmp,1,wrk,i)
        if (i.ne.0)             goto 100
        isummariz = iTmp(1)
        if (isummariz.lt.0)     goto 101
        if (isummariz.gt.1)     goto 101
        if (Index_corfu.lt.7  .OR.
     *      Index_corfu.gt.9) Then  !--+
          isummariz=0                 !|
        endif  !-----------------------+
c-----------------------
c This is a temporary solution (to ensure the compatibility
c on input files with Radicon TRDS program shell):
c (they do do NOT have lateral periodic reliefs
c (Index_corfu=10).
        if (Index_corfu.ne.10)  Then  !--+             ????????????
          Relief_Period = 0.            !|
          max_Periods = 0               !|
          goto 777  !--------------------+-------------------------+
        endif  !-------------------------+                         V
c-----------------------
        txt(7) = 'Number of correlated lateral periods (N_periods > 0)'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt  (iTmp,1,wrk,i)
        if (i.ne.0)              goto 100
        max_Periods = iTmp(1)
        if (max_Periods.lt.0)    goto 101
        if (Index_corfu.eq.10)  Then  !-------+         ! Periodic
c         if (max_Periods.eq.0.) goto 101    !|         ! relief
        else   !------------------------------|         ! model
          max_Periods = 0                    !|
        endif  !------------------------------+
c-----------------------
        txt(7) = 'Lateral period of interface relief (L_relief > 0)'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (Relief_Period,1,wrk,i)
        if (i.ne.0)              goto 100
        if (abs(Relief_Period).lt.1.E-32) Then !--+
          Relief_Period=2.*CorreLength           !|
        endif  !----------------------------------+
        if (Relief_Period.lt.0.) goto 101
        if (Index_corfu.eq.10 .AND.
     *      max_Periods.gt.0)   Then  !---------------+ ! Periodic relief
          if (Relief_Period.lt.CorreLength) goto 101 !|
        else   !--------------------------------------|
          Relief_Period = 0.                         !|
        endif  !--------------------------------------+
c-----------------------                                           V
  777   continue                !<---------------------------------+
        m0  = max_Periods
        if (m0.eq.0) Then  !-------------------+
          S = 0.                              !|
        else  !--------------------------------|
          S = (Relief_Period/CorreLength)**2  !|        ! S = (P/L)**2
        endif  !-------------------------------+
        S_factor = S
c---------------------------------------------------------------------
c           +inf.
c 1/A = 1+ 2*SUM exp[-(m*P/L_m)^{2h}-m/m0],
c            m=1
c
c L_m^2 = L^2*(1+|m/m0|*S),     S = (P/L)^2
c
c                   +inf.
c Then:   1/A = 1+ 2*SUM  exp[-(m^2*S/(1+S*m/m0))^h - m/m0]
c                    m=1
c---------------------------------------------------------------------
        A_factor = 1.
        do m=1,5*m0  !=================+
          mr = Real(m)/Real(m0)       !|
c         x = m*m*S/(1.+S*mr)         !| growing L_m
          x = m*m*S                   !| NON growing L_m
          x = exp(h_jagged*Alog(x))   !| x^h
          x = 2.*exp(-x-mr)           !|
          A_factor = A_factor + x     !|
          if (x.lt.1E-6)  goto 789    !|
        enddo  !=======================+
  789   A_factor = 1./A_factor
c-----------------------
        MaxOrder(1) = 14              !For coherent reflection calcs
        MaxOrder(2) = 7               !12 is bad for MS / 10 is bad for NDP
        iDebug      = 0
c ccc   iDebug      = 1               ! !!!! DEBUG !!!
        iBeepStep   = 1
        txt(7)      = ' '
c-----------------------
        Close   (unit=2)
c-----------------------
c Get substrate x0:
        ifail   = 0
        if (C00(N_Total).ne.'x')  Then !---------------------+
          mode_x0  = 0                        !new crystal   |
          iSyngony = 0                                      !|
          do i=1,kcompMax !==+                               |
            prcn(i) = 0.    !|                               |
          enddo  !===========+                               |
          xr0 = 0.                                          !|
          xi0 = 0.                                          !|
          Call Get_X0 (Wave, Code(1,N_Total), rho(N_Total), !|
     *                 iSyngony, mode_x0, xr0, xi0,         !|
     *                 x0(N_Total+1), 1, 1, ifail)          !|
          if (ifail.ne.0)            goto 130               !|
        endif !----------------------------------------------+
        Identical(N_Total) = 0
c A control for the substrate x0:
        if (abs(Real(x0(N_Total+1))).lt. 1.E-32 .OR.
     *           abs(W0(N_Total))   .lt. 1.E-32) goto 123
c#######################################################
        do      l=1,N_Top  !=================================+
          if (Identical(l).eq.0) Then  !-------------------+ |
                                                          !| |
            if ((Len_Trim(Code(1,l)).eq.0) .OR.           !| |
     *          (Code(1,l).eq.'Unknown'))  Then  !----+    | |
              Code(1,l) = Code(1,N_Total)            !|    | |
            endif  !----------------------------------+    | |
            x0_  = (0., 0.)                               !| |
            rho_ = 0.                                     !| |
            do j=1,N_Frac(l)  !==========================+ | |
              if (j.eq.1)  Then !-+                     !| | |
                rho_f = rho(l)   !|                     !| | |
              else  !-------------|                     !| | |
                rho_f = 0.       !|                     !| | |
              endif !-------------+                     !| | |
              if (Code(j,l).eq.Code(1,N_Total) .AND.    !| | |
     *           abs(rho_f).lt.1.E-32)         Then !--+ | | |
                x0_f  = x0(N_Total+1)                 !| | | |
c This means, even when the substrate density was     !| | | |
c specified externally, it is assigned to the layer   !| | | |
c with unspecified roughness:                         !| | | |
                rho_f = rho(N_Total)                  !| | | |
              else !-----------------------------------| | | |
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
                if (ifail.ne.0)      goto 130         !| | | |
              endif !----------------------------------+ | | |
              x0_  = x0_  + x0_f *Frac(j,l)             !| | |
              rho_ = rho_ + rho_f*Frac(j,l)             !| | |
            enddo  !=====================================+ | |
                                                          !| |
            if (C00(l).eq.' ') Then  !------------+        | |
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
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c #######################################################
c                  Error messages
c########################################################
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
        goto 120

  101   continue
        write   (txt,119)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line
  119   format  (a,': Parameter(s) out of range in file:'//a//
     *  'at line:',i4,' -- while reading:'/)                    !+txt(7)
        goto 120

  123   continue
        write   (txt,124)  ProgramVer(1:ipv), InpFile(1:linp)
  124   format  (a,': Zero substrate reflectivity in file:'//a//
     *  'Either a substrate code must be entered or x0#0.'/
     *  'Also make sure that W0#0.'/)
        goto 120

  131   continue
        write   (txt,132)  ProgramVer(1:ipv)
  132   format  (a,':'//
     *  'X-ray wavelength was chosen to be specified via'/
     *  'characteristic X-ray line, but no line given.')
        ifail = 4
        goto 130

  102   continue
        write   (txt,125)  ProgramVer(1:ipv), InpFile(1:linp)
  125   format  (a,': Inconsistent parameters in file:'//a//
     *  'Multiple scan points with zero scan range'/
     *  'while reading:')                                       !+txt(7)
        goto 120

  120   continue
        txt(8) = ' '
        if (ifail.ne.101)  Then  !-----------------+
          txt(9) = 'Check your input file!'       !|
        else  !------------------------------------|
          txt(9) = txt(10)          !other message!|
        endif  !-----------------------------------+
        ifail = 9
        goto 130

  130   continue
        if (modebat.eq.1) Then  !----------------+
          do    i=1,Iabs(ifail)  !=============+ |
            l = Max0(Len_Trim(txt(i)),1)      !| |
            write (*,'(1x,a)') txt(i)(1:l)    !| |
          enddo  !=============================+ |
        endif  !---------------------------------+
        Call    Message (txt,Iabs(ifail),2)
        goto 28

  28    continue
        Call exit_quiet()
        End
c
c =======================================================
c
        Subroutine Use_Symmetry8 (qz0, qx0, i, j,
     *                            N_Points_Max, mode_scan,
     *                            Scan_Range, Scan_Step,
     *                            DSSL, k4, UnitCoef,
     *                            qz_eps, qx_eps, Ok)

        Integer         N_Points_Max,
     *                  mode_scan,
     *                  i, j, l
        Real*8          DSSL(N_Points_Max)      ! array for intensity
        Real*4          Theta(2),               ! incid./exit angles
     *                  k4,                     ! real*4(k0)
     *                  qx, qx0, qxx, qx_eps,
     *                  qz, qz0, qzz, qz_eps,
     *                  UnitCoef,
     *                  Scan_Range(2,2),
     *                  Scan_Step(2),
     *                  Theta_Off, Theta_Smp
        Logical         Ok

        Ok = .False.

        if (i.lt.2)             goto 100        !return
        if (mode_scan.ne.1  .AND.
     *      mode_scan.ne.3  .AND.
     *      mode_scan.ne.4  .AND.
     *      mode_scan.ne.6)     goto 100        !return

        do      l=1,i-1    !============================+
                                                       !|
          goto (701,800,703,704,800,706) mode_scan  !--+|
c                1   2   3   4   5   6                !||
c+-----   Mode1: over Q  at 2Q=fix -----|              ||
  701     continue                                    !||
          Theta(1) = Scan_Range(1,1)                  !||
     +             + Scan_Step(1) * (l-1)             !||
          Theta(2) = Scan_Range(1,2) - Theta(1)       !||
          goto 711                                    !||
                                                      !||
c+------- Mode3: Q-2Q serie with offsets|              ||
  703     continue                                    !||
          Theta_Off = Scan_Range(1,1)                 !||
     +              + Scan_Step(1) * (l-1)            !||
          Theta_Smp = Scan_Range(1,2)                 !||
     +              + Scan_Step(2) * (j-1)            !||
          Theta(1) = Theta_Smp + Theta_Off            !||
          Theta(2) = Theta_Smp - Theta_Off            !||
                                                      !||
  711     continue                                    !||
          if (Theta(1).lt.0.)     goto 800            !||
          if (Theta(2).lt.0.)     goto 800            !||
          qx     = 0.5*k4*(Theta(1)**2 - Theta(2)**2) !||
     *           * UnitCoef*UnitCoef                  !||
          qz     = k4*(Theta(1) + Theta(2))           !||
     *           * UnitCoef                           !||
          goto 710                                    !||
                                                      !||
c+-----   Mode4: over qx at fixed qz ---|              ||
  704     continue                                    !||
          qx     = Scan_Range(1,1)                    !||
     +           + Scan_Step(1) * (l-1)               !||
          qz     = Scan_Range(1,2)                    !||
          goto 710                                    !||
                                                      !||
c+-----   Mode6: 3D over (qx,qz) -------|              ||
  706     continue                                    !||
          qx     = Scan_Range(1,1)                    !||
     +           + Scan_Step(1) * (l-1)               !||
          qz     = Scan_Range(1,2)                    !||
     +           + Scan_Step(2) * (j-1)               !||
          if (qz.le.0.) goto 800                      !||
                                                      !||
  710     continue                               !-----+|
          qzz = qz-qz0                                 !|
          qxx = abs(qx)-abs(qx0)                       !|
          if ((abs(qzz).lt.qz_eps)  .AND.              !|
     *        (abs(qxx).lt.qx_eps)) Then !---+          |
            DSSL(i) = DSSL(l)               !|          |
            Ok = .True.                     !|          |
            goto 100    !return             !|          |
          endif  !---------------------------+          |
  800     continue                                     !|
        enddo  !========================================+
100     continue
        return
        end
