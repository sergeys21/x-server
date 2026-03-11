c
c
c              +==================================+
c  Apr 2, 1996 |New method of evaluating Spiller's|
c              | correlation function, applicable |
c              |    at great q_x (like in GIDS)   |
c              +==================================+
c
c *** Version 1.0b ***  for TRDS
c                    -- Interface roughness profile
c                       (in dependence on the interface number)
c                    -- 6 models of roughness correlations
c                    -- Roughness asymmetry (b)
c
c +=====================================================+
c |   These routines compute x-ray diffuse scattering   |
c |  under grazing incidence specular reflection due to |
c |        interface roughness in multilayers           |
c | The calculations are carried out in 1 angular point |
c |      corresponding to incident angle = Theta(1)     |
c |              and exit angle Theta(2).               |
c |                                                     |
c |       Input data are transferred to DSFM in         |
c |                 common /TRDSdat/                    |
c |-----------------------------------------------------|
c |Contents:                                            |
c | 1. Real*8     function: DSFM                        |
c | 2. Complex*16 function: CorFun                      |
c | 3. Complex*16 function: Corr1                       |
c | 4. Complex*16 function: Corr2                       |
c | 5. Real*4     function: Spiller                     |
c | 6. Subroutine:          X_Sort                      |
c |-----------------------------------------------------|
c |Subroutine calls: see  --  TRSU2_1r(x).for,          |
c |                       --  RFMs_UNI.for              |
c +=====================================================+

c ####################################################### 1

        Real*8  Function  DSFM (Theta, N_Used, Rc, ipol, nfa, EscFile)
c  ** Diffuse_Scattering_from_Multi **
c--------------------------------------------------------
c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c ccc   Parameter  (N_Top_Max     = 150,
c ccc*              N_Total_Max   = N_Top_Max+1)

        Integer         N_Roots
        Parameter      (N_Roots  = 4)

        Complex*16      E(N_Roots,N_Total_Max+1,2),
     *                  E_i, E_j, E2,
     *                  u(N_Roots,N_Total_Max+1,2),
     *                  Summ, Sum1, Sum2,
     *                  CorFun, CF,
     *                  dx0(2), dx2,
     *                  ex2, q_i, q_j,
     *                  Im /(0.,1.)/

        Real*4          Theta(2),               ! incidence/exit angles
     *                  Rc(2),                  ! reflection coefficient
     *                  Coeff, Err,             ! work cells
     *                  q1max, q2max, q,        ! work cells
     *                  width /1.E-6/           !width of Bragg peak over qx

        Integer*4       nfa, ifa,               ! number of fails, last failure code
     *                  Icount4,                ! loops counter
     *                  Ifull4                  ! loops total

        Integer         N_Used(2),              ! number of used layers
     *                  N_Use,                  ! min number of used layers
     *                  N_Fields(2)/2,2/,       ! # of wavefields for incident and scattered waves
     *                  Mode_Sigma,             ! flag to account for sigma
     *                  ipol, i, j,
     *                  l1_max, l2_max,
     *                  m1_max, m2_max,
     *                  l1, l2, m1, m2,
     *                  iscan, iasci,
     *                  ipcn, iii

        Logical         Compute_Fields          ! for roughness

        Character       EscFile*(*)

        Logical*4       FileExist
        External        FileExist

c------------------------------------------------
        Complex*8       x0(N_Total_Max+1)
        Real*8          Z(N_Total_Max),         ! coordinates of interfaces
     *                  k0                      ! k0 = 2*pi/wave
        Real*4          Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Sigma_Accu(N_Total_Max),
     *                  Wave, UnitCoef(3),
     *                  qz, qx,                 ! (k2+k1)_z & (k2+k1)_x
     *                  CorreLength,            ! in-plane correlation length(A)
     *                  CorreVert,              ! vertical correlation length(A)
     *                  h_jagged,               ! jaggendness of the surface(0.:1.)
     *                  xabs,                   ! abs(x0_max)
     *                  TERang,                 ! critical TER angle
     *                  pi                      ! 3.1416
        Integer         i_Spec,
     *                  N_Top,                  ! top layer sublayers
     *                  N_Total,                ! total of layers
     *                  Index_corfu,            ! index of correlation function
     *                  MaxOrder(2),            ! max allowed S-matrix element (depricated)
     *                  iDebug,                 ! Debug flag
     *                  ifail                   ! failure code
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

        Character  txt(20)*80
        Common  /msg/   txt
c--------------------------------------------------------
        DSFM  = 0.0D0
        ifail = 0
        ifa   = 0
        nfa   = 0                               !errors counter: increment for non-critical errors only!

c Coherent reflection:
        if (i_Spec.eq.1 .AND. Abs(qx).lt.width) Then  !--+
          Mode_Sigma     = 1                            !|
          Compute_Fields = .False.                      !|
          Call  SFM    (Theta(1),                       !|
     *                  Rc(1),                          !|
     *                  N_Fields(1),                    !|
     *                  N_Used(1),                      !|
     *                  E(1,1,1),                       !|
     *                  u(1,1,1),                       !|
     *                  ipol, Mode_Sigma,               !|
     *                  Compute_Fields)                 !|
          Rc(2) = Rc(1)                                 !|
          DSFM  = Dble(Rc(1))                           !|
          return                                        !|
        endif  !-----------------------------------------+

        Mode_Sigma     = 0
        Compute_Fields = .True.
        if (i_Born.lt.2)  Then  !--------------+0=DWBA,1=semi-Born
          do    i=1,2  !=====================+ |
            Call  SFM    (Theta(i),         !| |
     *                    Rc(i),            !| |
     *                    N_Fields(i),      !| |
     *                    N_Used(i),        !| |
     *                    E(1,1,i),         !| |
     *                    u(1,1,i),         !| |
     *                    ipol, Mode_Sigma, !| |
     *                    Compute_Fields)   !| |
            if (ifail.ne.0) Then !--+        | |
c              nfa = nfa + 1       !|        | |combination ifail#0 & nfa=0 is
               return              !|        | |considered as a critical error
            endif !-----------------+        | |
          enddo  !===========================+ |
        else   !-------------------------------+2=Born
          do    i=1,2  !============+          |
            N_Used(i) = N_Total    !|          |
          enddo  !==================+          |
        endif  !-------------------------------+

        N_Use = Min0(N_Used(1),N_Used(2))
c+=========================================================================+
c|The transition from the cross section to the reflection coefficient is   |
c|according to the formulae from the paper by S.K.Sinha, E.B.Sirota,       |
c|S.Garoff, H.Stanley -- Phys.Rev.B v.38(4), p.2297-2311):                 |
c|                                                                         |
c|              1       ?dSigma                                            |
c|   P_r = ------------ |------ dOmega                   -- Sinha,.. (2.13)|
c|         S*sin(Phi_0) ?dOmega                                            |
c|         +----------+                                                    |
c|              +<-incident beam cross section                             |
c|                                                                         |
c|where dOmega is the dimensionless solid angle which can be expressed as: |
c|                                                                         |
c|                                   dq_x         dq_y                     |
c|dOmega = d(Phi_h) * d(Theta) = -------------- * ----   -- Sinha,.. (2.14)|
c|                               k * sin(Phi_h)    k                       |
c|                                                                         |
c|-- the latter equation follows in reflectometry from:                    |
c|q_x = k*(cos(Phi_0)-cos(Phi_h)),  |=>  d(q_x)/k = sin(Phi_h) d(Phi_h)    |
c|-------------------------------------------------------------------------|
c|Then, at specular we can write:            +--------Delta-functions!     |
c|                                      +----------+                       |
c|   dSigma/dOmega = (dSigma/dOmega)'*delta(q_x)*delta(q_y)*k^2            |
c|                   +--------------+                                      |
c|                          +-this is a "specular part" of cross section   |
c|Carrying out the integration in (2.13), we obtain:                       |
c|                                                                         |
c|   P_r = (dSigma/dOmega)' / (S*sin^2(Phi_0)),          -- Sinha,.. (4.32)|
c|                                                                         |
c|since Phi_0 = Phi_h at specular.                                         |
c|                                                                         |
c|To use Eq.(2.13), we must have dimensionless (dSigma/dOmega)/S.          |
c|-------------------------------------------------------------------------|
c|To re-calculate the data of TRDS into reflection coefficient one must    |
c|apply the following relation:                                            |
c|                                                                         |
c| P_r = TRDS * d(Phi_h) / sin(Phi_0)                                      |
c|                                                                         |
c|where d(Phi_h) is the detector angular resolution.                       |
c+=========================================================================+
c Our (dSigma/dOmega) is integrated over q_y. It gives additional 1/k0.
c As CorFun is ? CorreLength * Sigma * Sigma ? (Angstrem)^3, the cross
c section is dimensionless:
        Coeff   = SNGL(k0*k0*k0/(16.*pi*pi))
        Summ    = (0.,0.)
        Icount4 = 0
        Ifull4  = ((N_Use+0)*(N_Use+1)/2)
     *          * (N_Fields(1)*N_Fields(2))**2
        if (i_Born.eq.0)  Then !--------+0=DWBA
          l1_max = N_Fields(1)         !|
          m1_max = N_Fields(2)         !|
          l2_max = N_Fields(1)         !|
          m2_max = N_Fields(2)         !|
        else  !-------------------------+1=semi-Born,2=Born
          l1_max = 1  !! N_Fields(1)/2 !|
          m1_max = 1  !! N_Fields(2)/2 !|
          l2_max = 1  !! N_Fields(1)/2 !|
          m2_max = 1  !! N_Fields(2)/2 !|
        endif  !------------------------+

        do     i = 1,N_Use  !=====================================+
          if (abs(Sigma(i)).lt.1.E-32) Then !-------+             |
            Icount4 = Icount4 + (N_Use-i+1)        !|             |
     *              * (N_Fields(1)*N_Fields(2))**2 !|             |
            goto  3  !---------------->-------------+-->----------+--+
          endif  !----------------------------------+             |  |
          dx0(1) = x0(i)-x0(i+1)                                 !|  V
          dx0(1) = Dconjg(dx0(1))                                !|
          Sum1   = (0.,0.)                                       !|
        do     j = i,N_Use  !===================================+ | From "i"
          if (abs(Sigma(j)).lt.1.E-32) Then !-------+             |
            Icount4 = Icount4                      !|           | |
     *              + (N_Fields(1)*N_Fields(2))**2 !|           | |
            goto  2  !---------------->-------------+-->--------+-+-+
          endif  !----------------------------------+           | | |
c For Spiller model:                                            | | V
          if (Index_corfu.eq.6) Then  !---------------------+   | |
            if (i_Fourier.eq.0) Then  !-------------------+ |   | |
c Estimate Max(q_i*q_j):                                  | |   | |
              q1max = 0.                                 !| |   | |
              q2max = 0.                                 !| |   | |
              do  l1 = 1,N_Fields(1)   !==============+   | |   | |
              do  m1 = 1,N_Fields(2)   !============+ |   | |   | |
                q = SNGL(                          !| |   | |   | |
     (              Abs(u(l1,i+1,1) + u(m1,i+1,2)) !| |   | |   | |
     )              )                              !| |   | |   | |
                if (q.gt.q1max)   q1max = q        !| |   | |   | |
                                                   !| |   | |   | |
                q = SNGL(                          !| |   | |   | |
     (              Abs(u(l1,j+1,1) + u(m1,j+1,2)) !| |   | |   | |
     )              )                              !| |   | |   | |
                if (q.gt.q2max)   q2max = q        !| |   | |   | |
              enddo  !==============================+ |   | |   | |
              enddo  !================================+   | |   | |
              q2max = SNGL(q1max*q2max*((k0*TERang)**2)) !| |   | |
            endif  !--------------------------------------+ |   | |
c Tabulate Spiller's function:                             !|   | |
            Call Spiller_Fill (i,j,q2max)                  !|   | |
          endif  !------------------------------------------+   | |
                                                               !| |
          dx0(2) = x0(j)-x0(j+1)                               !| |
          dx2    = Coeff*dx0(1)*dx0(2)                         !| |
          Sum2   = (0.,0.)                                     !| |
        do    l1 = 1,l1_max  !================================+ | |
        do    m1 = 1,m1_max  !==============================+ | | |
        do    l2 = 1,l2_max  !============================+ | | | |
        do    m2 = 1,m2_max  !==========================+ | | | | |
                                                       !| | | | | |
          Icount4 = Icount4 + 1                        !| | | | | |
          Call  KeyIn (0,iscan,iasci)                  !| | | | | |
          if (iasci.eq.27)      Then  !---------------+ | | | | | |
            ipcn = 100 * Icount4 / Ifull4            !| | | | | | |
            write (txt(9),10)   Icount4,Ifull4,      !| | | | | | |
     *                          ipcn,ifa,nfa         !V V V V V V V
  10        format('Loops=',i6,'/',i6,' (',i3,'%), ',
     *      'ifail=',i3,' (*',i5,');   Interrupt?')  !^ ^ ^ ^ ^ ^ ^
            ifail = -999                             !| | | | | | |
            return                                   !| | | | | | |
          endif  !------------------------------------+ | | | | | |
          if (FileExist(EscFile)) Then !--+             | | | | | |
            Call DelFile (EscFile,iii)   !|             | | | | | |
            ifail = -999                 !|             | | | | | |
            return                       !|             | | | | | |
          endif  !------------------------+             | | | | | |
                                                       !| | | | | |
          if (i_Born.lt.2) Then  !-----------------+    | | | | | |0=DWBA,1=semi-Born
                                                  !|    | | | | | |
            q_i = (u(l1,i+1,1) + u(m1,i+1,2))     !|    | | | | | |
     *          * k0 * TERang                     !|    | | | | | |
            q_i = Dconjg(q_i)                     !|    | | | | | |
                                                  !|    | | | | | |
            q_j = (u(l2,j+1,1) + u(m2,j+1,2))     !|    | | | | | |
     *          * k0 * TERang                     !|    | | | | | |
                                                  !|    | | | | | |
            ex2 = (Sigma_Accu(i)*q_i)**2          !|    | | | | | |
     *          + (Sigma_Accu(j)*q_j)**2          !|    | | | | | |
            ex2 = Exp(-0.5*ex2)                   !|    | | | | | |
                                                  !|    | | | | | |
            E_i = E(l1,i+1,1) * E(m1,i+1,2)       !|    | | | | | |
            E_j = E(l2,j+1,1) * E(m2,j+1,2)       !|    | | | | | |
            E_i = Dconjg(E_i)                     !|    | | | | | |
            E2  = E_i * E_j                       !|    | | | | | |
                                                  !|    | | | | | |
          else  !----------------------------------+    | | | | | |2=Born
                                                  !|    | | | | | |
            q_i = qz                              !|    | | | | | |
            q_j = qz                              !|    | | | | | |
            ex2 = (Sigma_Accu(i)*qz)**2           !|    | | | | | |
     *          + (Sigma_Accu(j)*qz)**2           !|    | | | | | |
            ex2 = Exp(-0.5*ex2)                   !|    | | | | | |
c R=0, T=1 (multiplied by respective phase):      !|    | | | | | |
            E2  = exp(Im*qz*(Z(i)-Z(j)))          !|    | | | | | |
c -- here exp(-i*qz...) provides nearly same result|??? | | | | | |
          endif !----------------------------------+    | | | | | |
                                                       !| | | | | |
          if (abs(E2).lt.1.E-32) goto 1  !----->-+      | | | | | |
          CF    = CorFun(i, j, q_i*q_j)         !|      | | | | | |
c+--------------------------------------+        |      | | | | | |
c|DEBUG: calculations without CorFun:   |        |      | | | | | |
c+-                                     |        |      | | | | | |
c         CF    = 2.* CorreLength      !|        |      | | | | | |
c    *          * Sigma(i) * Sigma(j)  !|        |      | | | | | |
c+-                                     |        |      | | | | | |
c|cc *          * ArrCos(1)?1           |        |      | | | | | |
c+--------------------------------------+        |      | | | | | |
          if (ifail.ne.0)       Then  !--+       |      | | | | | |
            nfa = nfa + 1               !|       |      | | | | | |non-critical error
            ifa = ifail                 !|       |      | | | | | |
          endif  !-----------------------+       |      | | | | | |
          Sum2  = Sum2 + dx2 * ex2 * E2 * CF    !|      | | | | | |
  1       continue  !----------------------------+      | | | | | |
        enddo  !========================================+ | | | | |
        enddo  !==========================================+ | | | |
        enddo  !============================================+ | | |
        enddo  !==============================================+ | |
          if (j.ne.i)   Sum2 = 2.*Real(Sum2)           !        | |
          Sum1  = Sum1 + Sum2                          !        | | V
  2       continue                                     !<-------+-+-+
        enddo  !================================================+ |
          Summ  = Summ + Sum1                          !          |  V
  3       continue                                     !<---------+--+
        enddo  !==================================================+
        ifail = ifa                             !this copies last error

c cc    DSFM = Abs (Summ)

        DSFM = Real (Summ)

        if (DSFM.lt.0.) Then  !----------------------------------+
          write   (3,57) Summ                                   !|
  57      format (' Error 557: Cross section=(',2g10.3,') < 0') !|
          ifail = 557                                           !|
          nfa   = nfa + 1                                       !|non-critical error
          return                                                !|
        endif  !-------------------------------------------------+

        if (abs(DSFM).gt.1.E-32 .AND. Index_corfu.lt.7) Then !---+
          Err = SNGL(100.*Abs(Imag(Summ)/DSFM))                 !|
          if (Err.gt.5.)  Then  !-----------------------+        |
            write   (3,55)      Err                    !|        |
  55        format (' Error 555: Im(Sigma)=',g10.3,'%')!|        |
            ifail = 555                                !|        |
            nfa   = nfa + 1                            !|        |non-critical error
            return                                     !|        |
          endif  !--------------------------------------+        |
        endif  !-------------------------------------------------+

        return
        end

c ####################################################### 2

        Complex*16      Function  CorFun (i,j,q_ij)
c  ** Correlation function of interface roughness **
c--------------------------------------------------------
c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter  (N_Top_Max     = 150,
c cc *              N_Total_Max   = N_Top_Max+1)

        Complex*16      q_ij,                   ! q_i * q_j
     *                  IM/(0.,1.)/,            ! i
     *                  qw                      ! work cell

        Real*8          dZ,                     ! z1 - z2
     *                  Corr1m                  ! corr.integral

        Real*4          Sigma_ij                ! Sigma_i * Sigma_j

        Integer         i, j, m, n
c------------------------------------------------
        Complex*8       x0(N_Total_Max+1)
        Real*8          Z(N_Total_Max),         ! coordinates of interfaces
     *                  k0                      ! k0 = 2*pi/wave
        Real*4          Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Sigma_Accu(N_Total_Max),
     *                  Wave, UnitCoef(3),
     *                  qz, qx,                 ! (k2+k1)_z & (k2+k1)_x
     *                  CorreLength,            ! in-plane correlation length(A)
     *                  CorreVert,              ! vertical correlation length(A)
     *                  h_jagged,               ! jaggendness of the surface(0.:1.)
     *                  xabs,                   ! abs(x0_max)
     *                  TERang,                 ! critical TER angle
     *                  pi                      ! 3.1416
        Integer         i_Spec,
     *                  N_Top,                  ! top layer sublayers
     *                  N_Total,                ! total of layers
     *                  Index_corfu,            ! index of correlation function
     *                  MaxOrder(2),            ! max allowed S-matrix element (depricated)
     *                  iDebug,                 ! Debug flag
     *                  ifail                   ! failure code
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
     *                    i_Born,     ! approximation for wavefields
     *                    i_Fourier,  ! corr.integral: exp[K(x)]-1 or K(x)
     *                    isummariz,  ! add roughness-to-steps flag
     *                    Relief_Period,!surface relief period
     *                    max_Periods,! number of correlated periods
     *                    S_factor,   !==(Relief_Period/CorreLength)^2
     *                    A_Factor    ! normalization for periodic func.


        Complex*16      Corr1,                  ! corr.integral
     *                  Corr2,                  ! corr.integral (Spiller)
     *                  Corr3                   ! corr.integral (Pikute)
        External        Corr1,
     *                  Corr2,
     *                  Corr3
c--------------------------------------------------------
        ifail  = 0
        CorFun = (0.0D0, 0.0D0)
        goto(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
     *                                                     Index_corfu !---------+
c                                                                       VVVVVVVVVV
        ifail  = 700
        write (3,*)     'CorFun: failure code=',ifail,
     *                  '-- at i,j: ',i,j
        return

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  1     continue
c Uncorrelated roughness:

        if (i.ne.j)             return
        Sigma_ij = Sigma(i) * Sigma(j)
        CorFun = Corr1 (Sigma_ij,q_ij,qx,CorreLength,i_Fourier)
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  2     continue
c Complete vertical correlation of interface roughness:

        Sigma_ij = Sigma(i) * Sigma(j)
        CorFun = Corr1 (Sigma_ij,q_ij,qx,CorreLength,i_Fourier)
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  3     continue
c Z.H.Ming, A.Krol, Y.L.Soo, Y.H.Kao, J.S.Park, K.L.Wang --
c Phys.Rev. B v.47(24) p.16373-16381 (1993)

        Sigma_ij = Sigma(i) * Sigma(j)
        if (i.ne.j)     Then  !--------------------------+i#j
          if (abs(CorreVert).lt.1.E-32) Then !---------+ |
            return                                    !| |
          else  !--------------------------------------+ |
            dZ = Abs(Z(i)-Z(j))/CorreVert             !| |
            if (dZ.lt.1.0D+10)  Then !---------------+ | |
              Sigma_ij = SNGL(Sigma_ij * Exp(-dZ))  !| | |
            else !-----------------------------------+ | |
              return                                !| | |
            endif  !---------------------------------+ | |
          endif  !-------------------------------------+ |
        endif  !-----------------------------------------+
        CorFun = Corr1 (Sigma_ij,q_ij,qx,CorreLength,i_Fourier)
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  4     continue
c Y.H.Phang, D.E.Savage, R.Kariotis, M.G.Lagally --
c J.Appl.Phys. v.74(5) p.3181-3188 (1993).
c
c We use here the simplified form of
c Phang`s correlation function
c (cp. p.3185 of their article)

        Sigma_ij = Sigma(i) * Sigma(j)
        if (i.eq.j)     Then  !----------------------------------+i=j
          CorFun = Corr1(Sigma_ij,q_ij,qx,CorreLength,i_Fourier)!|
c                                         ===========            |
          goto 99                                               !|
        else  !--------------------------------------------------|i#j
          if (abs(CorreVert).lt.1.E-32) Then !--------+          |
            return                                   !|          |
          else  !-------------------------------------|          |
            dZ = Abs(Z(i)-Z(j))/CorreVert            !|          |
            if (dZ.lt.1.0D+10)  Then !--------------+ |          |
              Sigma_ij = SNGL(Sigma_ij * Exp(-dZ)) !| |          |
            else   !--------------------------------| |          |
              return                               !| |          |
            endif  !--------------------------------+ |          |
          endif  !------------------------------------+          |Is this model
          CorFun = Corr1 (Sigma_ij,q_ij,qx,CorreCross,i_Fourier)!|compatible
c                                          ==========            |with
          goto 99                                               !|inclination?
        endif  !-------------------------------------------------+?????????

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  5     continue
c V.Holy and T.Baumbach -- Phys.Rev. B v.49(15)
c p.10668-10676 (1994)

        Sigma_ij = 0.
        m = max0(i,j)
        do      n=m,N_Total  !==============+
          Sigma_ij = Sigma_ij              !|
     +             + Sigma(n) * Sigma(n)   !|
        enddo  !============================+
        CorFun = Corr1 (Sigma_ij,q_ij,qx,CorreLength,i_Fourier)
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  6     continue
c E.Spiller, D.Stearns, M.Krumley -- J.Appl.Phys. v.74(1)
c p.107-118 (1993)
c
c We restrict ourselves with the
c case h=1 (Gaussian distribution)
c For this case, the formulae are
c given by Kaganer, Stepanov & Kohler.

        CorFun = Corr2 (i,j,q_ij)
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  7     continue
  8     continue
  9     continue
c Similar to P.R.Pukite, C.S.Lent and P.I.Cohen --
c Surf.Sci. v.161 p.39-68 (1985)

        Sigma_ij = Sigma(i) * Sigma(j)
        if (i.ne.j)     Then  !-------------------+i#j
          if (abs(CorreVert).lt.1.E-32) Then !--+ |
            return                             !| |
          else  !-------------------------------+ |
            dZ = Abs(Z(i)-Z(j))/CorreVert      !| |
            if (dZ.lt.1.0D+10) Then !--+        | |
              dZ = Exp(-dZ)           !|        | |
            else  !--------------------|        | |
              return                  !|        | |
            endif !--------------------+        | |
          endif  !------------------------------+ |
        else  !-----------------------------------+
          dZ = 1.                                !|
        endif !-----------------------------------+
        CorFun = dZ * Corr3 (Sigma_ij,q_ij,qz,qx,CorreLength,
     *                       SurMiscut*UnitCoef(3),TerraceSpre,
     *                       StepHeight,Index_corfu,isummariz)
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  10    continue
c Periodical surface relief with vertical correlation
c like in Ming, et al.

        Sigma_ij = Sigma(i) * Sigma(j)
        if (i.ne.j)     Then  !-------------------------+i#j
          if (abs(CorreVert).lt.1.E-32) Then !--------+ |
            goto 99                                  !| |
          else  !-------------------------------------| |
            dZ = Abs(Z(i)-Z(j))/CorreVert            !| |
            if (dZ.lt.1.0D+10)  Then !--------------+ | |
              Sigma_ij = SNGL(Sigma_ij * Exp(-dZ)) !| | |
            else   !--------------------------------| | |
              goto 99                              !| | |
            endif  !--------------------------------+ | |
          endif  !------------------------------------+ |
        endif  !----------------------------------------+
c q_ij is not needed for Corr1m, because it uses the small
c roughness approximation:
        CorFun   = Corr1m (Sigma_ij,
c ccc*                     q_ij,
     *                     qx, CorreLength,
     *                     Relief_Period, max_Periods,
     *                     A_factor, S_factor,
     *                     i_Fourier)
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  11    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  12    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  13    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  14    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  15    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  16    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  17    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  18    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  19    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++
  20    continue
        goto 99

c ++++++++++++++++++++++++++++++++++++++++++++++++++++

  99    continue
        if (Index_corfu.eq.6)   goto 100
        if (abs(RouMiscut).gt.1.E-32) Then !-----------------+
          qw = IM*qx*Tan(RouMiscut*UnitCoef(2))*(Z(i)-Z(j)) !|
          CorFun = CorFun * exp (qw)                        !|
        endif  !---------------------------------------------+
  100   continue
        return
        end

c ####################################################### 3

        Real*8     Function  Corr1m (Sigma_ij,
c ccc*                               q_ij,      !not used in the small roughn.approx.
     *                               qx,
     *                               CorreLength,
     *                               Relief_Period,
     *                               max_Periods,
     *                               A_factor,
     *                               S_factor,
     *                               i_Fourier)
c  ** The correlation integral of interface roughness **
c--------------------------------------------------------
c       Complex*16      q_ij                 ! q_i*q_j

        Real*8          Corr_m               ! work cells

        Real*4          qx,                  ! (k2+k1)_x
     *                  CorreLength, L, L_m, ! in-plane correl.length(A)
     *                  Relief_Period, P,    ! in-plane period (A)
     *                  Sigma_ij,            ! Sigma_i * Sigma_j
     *                  G_n,                 ! n-th integral
     *                  B2,                  ! work cell
     *                  mr,                  ! m/m0
     *                  S, S_factor,
     *                  A, A_factor

        Integer         max_Periods,
     *                  i_Fourier,
     *                  m0, m

        Real*4          Cos_Integral
        External        Cos_Integral
c+====================================================================+
c|The correlation function is postulated in the form:                 |
c|                                                                    |
c|                  +inf.                        2h                   |
c| K(x) = Sigma_ij * SUM  A_m * exp[-(|x-m*P|/L_m) ]                  |
c|                 m=-inf.                                            |
c|                                                                    |
c| where A_m = A*exp(-|m|/m0) ,  L_m = L*sqrt(1+|m/m0|*S)             |
c|                                            ---- "Random Walk" model|
c| Here we introduced a speed S of smearing. If we want the steps to  |
c| reach overlapping at m0, then: S=(P/L)^2                           |
c| /i.e. at S>>1 and m=m0 we have L_m=P/                              |
c|                                                                    |
c| m0=max_Periods, P=Relief_Period, L=CorreLength                     |
c|                                                                    |
c| The constant A is given by K(0)=Sigma_ij:                          |
c|                                                                    |
c| At m0=0  A=1                                                       |
c|                                                                    |
c| At m0=1  A=1/{1+2*exp[-(P/L*sqrt(2))^2h-1]+                        |
c|                +2*exp[-(2P/L*sqrt(3))^2h-2]}                       |
c|                                                                    |
c| At m0=2  A=1/{1+2*exp[-(P/L*sqrt(1.5))^2h-.5]+                     |
c|                +2*exp[-(2P/L*sqrt(2))^2h-1]}                       |
c|                                                                    |
c| At m0>>1 A=1/{1+2*exp[-(P/L)^2h]}                                  |
c|                                                                    |
c| -- since we presume that P >> L.                                   |
c|                                                                    |
c|The correlation integral is expanded using the                      |
c|SMALL ROUGHNESS APPROXIMATION (we denote: qz^2=q_i * q^{*}_j):      |!!!!!!
c|                                                                    |
c|     +inf   +              +            +inf.                       |
c|  1    |    | qz^2 * K(x)  |  -iq_x x     |                         |
c| ----  | dx |e          - 1| e        =   | dx K(x) cos(q_x*x)  =   |
c| qz^2  |    +              +              |                         |
c|     -inf                               -inf.                       |
c|                                                                    |
c|                    +inf.                                           |
c|     inf.             |                                             |
c| = B*SUM exp(-|m|/m0) |dx exp[-(|x-m*P|/L_m)^{2h}] cos(q_x*x) = ... |
c|    -inf.             |                                             |
c|                    -inf.                                           |
c|                                                                    |
c|where B = A * Sigma_ij                                              |
c|                                                                    |
c|First, we are making the substitution:                              |
c|          y=(x-m*P)/L_m   =>  x=y*L_m+m*P   =>  dx=dy               |
c|                                                                    |
c|                  +inf.                                             |
c|   inf.    -|m|     |            2h                                 |
c|=B*SUM exp(---) L_m |dy exp(-{|y|  ) cos(L_m*q_x*y+m*P*q_x) = ...   |
c|  -inf.    m0       |                                               |
c|                  -inf.                                             |
c|-- now we are expanding the cosine(a+b)=cos(a)*cos(b)-sin(a)*Sin(b),|
c|and the second term vanishes since it is the integral from an odd   |
c|function in symmetric limits. So, we have:                          |
c|                                                                    |
c|                                 +inf.                              |
c|     inf.    -|m|                  |           2h                   |
c|=2*B*SUM exp(---)*cos(m*P*q_x)*L_m |dy exp(-|y|  ) cos(L_m*q_x*y) = |
c|    -inf.    m0                    |                                |
c|                                   0                                |
c|     inf.                                                           |
c|=2*B*SUM exp(-|m|/m0)*cos(m*P*q_x)*L_m*G(L_m*q_x)                   |
c|    -inf.                                                           |
c|                                                                    |
c|where we denoted:                                                   |
c|                                                                    |
c|       +inf   +               +                                     |
c|         |    | -y^{2h}       |                                     |
c| G(p) =  | dy |e       cos(px)|.                                    |
c|         |    +               +                                     |
c|         0                                                          |
c|                                                                    |
c|Integral G(P) is tabulated in subroutine CorrFun2 for p={0.,100.}   |
c|and tranferred as the array "ArrCos" in common /ArrcoDat/           |
c|(see CORRFUN2.FOR).                                                 |
c|--------------------------------------------------------------------|
c|          +inf.                                                     |
c|1/A = 1+ 2*SUM exp[-(m*P/L_m)^{2h}-m/m0],                           |
c|           m=1                                                      |
c|                                                                    |
c|L_m^2 = L^2*(1+|m/m0|*S),     S = (P/L)^2                           |
c|                                                                    |
c|                  +inf.                                             |
c|Then:   1/A = 1+ 2*SUM  exp[-(m^2*S/(1+S*m/m0))^h - m/m0]           |
c|                   m=1                                              |
c| -- this is calculated in GIDS0_xx.for (Read_gids_Input)            |
c+====================================================================+
        Corr1m = 0.0D0

        if (i_Fourier.ne.1) return

c+---------------------------------------------------------------+ |
c|               inf.                                            | |
c|Integral = 2*B*SUM exp(-|m|/m0)*cos(m*P*q_x)*L_m G(L_m*q_x) =  | |
c|              -inf.                                            | |
c|                  inf.                                         | |
c|=2B*[G(L*q_x)*L+2*SUM exp(-|m|/m0)*cos(m*P*q_x)*L_m*G(L_m*q_x)]| |
c|                 -inf.                                         | |
c+---------------------------------------------------------------+ |
        m0 = max_Periods
        P  = Relief_Period
        L  = CorreLength
        S  = S_factor
        A  = A_factor
        B2 = 2.*A * Sigma_ij
        G_n = Cos_Integral(L*qx)
        Corr1m = B2*L*G_n
        if (Abs(Corr1m).lt.1.E-32) return
        do m=1,5*m0  !==================================+
          mr  = Real(m)/Real(m0)                       !|
c         L_m = L*sqrt(1+S*mr)                         !| growing L_m
          L_m = L                                      !|NON growing L_m
          G_n = Cos_Integral(L_m*qx)                   !|
          Corr_m = 2*B2*exp(-mr)*cos(m*P*qx)*L_m*G_n   !|
          Corr1m = Corr1m + Corr_m                     !|
          if (Abs(Corr_m/Corr1m).lt.1.E-5) goto 1  !----+-----+
        enddo  !========================================+     |
                                                             !V
  1     continue  !<------------------------------------------+
        if (Corr1m.le.0.) Then
          write (3,*) '----------'
          write (3,*) Corr1m, Corr_m
          write (3,*) m, m0, mr
          write (3,*) L, L_m
          write (3,*) S, A, B2
          write (3,*) P, m*P*qx
        endif
        return
        end

c ####################################################### 4

        Complex*16      Function  Corr1 (Sigma_ij,q_ij,qx,
     *                                  CorreLength,i_Fourier)
c  ** The correlation integral of interface roughness **
c--------------------------------------------------------
        Complex*16 q,              ! work cells
     *             q_ij,           ! q_i*q_j
     *             q1,             ! (q_ij)*(Sigma_ij)
     *             qn              ! q1**n

        Real*4     qx,             ! (k2+k1)_x
     *             CorreLength,    ! in-plane correlation length(A)
     *             Sigma_ij,       ! Sigma_i * Sigma_j
     *             G_n,            ! n-th integral
     *             p1,             ! CorreLength*qx
     *             pn,             ! CorreLength*qx/n^{1/2h}
     *             n_fact          ! n!

        Integer*4  n

        Integer    i_Fourier

        Real*4     n12h(20)        ! n^{1/2h}
        Common     /n12/ n12h

        Real*4     Cos_Integral
        External   Cos_Integral
c+====================================================================+
c|The correlation integral is expanded into a series of               |
c|up to 10 integrals following to S.K.Sinha -- J.Physique III         |
c|(France) v.4, p.1543-1557 (1994)  --- Eq.(25)                       |
c|                                                                    |
c|   +infinity +   2     2 -(|x|/L)^{2h}  +                           |
c|  1    |     | qz sigma e               |  -iq_x x                  |
c| ----  |  dx |e                      - 1| e        =                |
c| qz^2  |     |                          |                           |
c|   -infinity +                          +                           |
c|                                    n-1                             |
c|            +infinity (qz^2*sigma^2)         + q_x * L  +           |
c| = 2*sigma^2  Summ    --------------- * L * G| -------- |           |
c|              n=1     n! * n^{1/2h}          + n^{1/2h} +           |
c|                          +-------+                                 |
c|  +---------------------------+                                     |
c| Not included in the paper by Sinha (propably, misprint).           |
c|                                                                    |
c|where (qz*sigma)^2 ={q_i * q^{*}_j * sigma_i * sigma_j},  and:      |
c|                                                                    |
c|        +infinity +               +                                 |
c|            |     | -x^{2h}       |                                 |
c| G(p)  =    |  dx |e       cos(px)|                                 |
c|            |     |               |                                 |
c|            0     +               +                                 |
c|                                                                    |
c|Integral G(P) is tabulated in subroutine CorrFun2 for p={0.,10}     |
c|and tranferred as array "ArrCos" in common /ArrcoDat/               |
c|(see the function Cos_Integral in the module CORRFUN2.FOR).         |
c+====================================================================+
        q1    = q_ij * Sigma_ij
        p1    = CorreLength * qx
        Corr1 = (0.0D0, 0.0D0)

        if (Abs(q1).lt.0.01 .OR. i_Fourier.eq.1) Then !-+
c The accuracy is 0.01/2 = 5.E-3                        |
          G_n = Cos_Integral (p1)                      !|
                                                       !|
          Corr1 = 2.* CorreLength * G_n * Sigma_ij     !|
                                                       !|
        else   !----------------------------------------|
                                                       !|
        n_fact = 1                                     !|
        qn     = 1.0D+0                                !|
c The number of steps: < 20 (dimension of n12h)         |
        do   n = 1,10   !=============================+ |
          n_fact = n_fact * n         ! n_fact =n!    | |
          qn  = qn * q1               ! qn=q1**n      | |
          pn  = p1 / n12h(n)          ! pn=p1/n^{1/2h}| |
          G_n = Cos_Integral (pn)                    !| |
          q   = 2. * CorreLength * G_n * qn          !| |
     /        / (n_fact * n12h(n) * q_ij)            !| |
          Corr1 = Corr1 + q                          !| |
c Accuracy:                                           | |
          if (Abs(q).lt.(5.D-3)*Abs(Corr1))  Return  !| |
        enddo  !======================================+ |
        endif  !----------------------------------------+
        return
        end

c ####################################################### 5

        Complex*16      Function  Corr2 (i,j,q_ij)
c  ** The integral from Spiller-type correlation function **
c                **  of interface roughness **
c--------------------------------------------------------
c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter  (N_Top_Max     = 150,
c cc *              N_Total_Max   = N_Top_Max+1)

        Complex*16      q_ij,                   ! q_i * q_j
     *                  Sum16,                  ! work cell
     *                  Sum16_2,                ! work cell
     *                  Sum16_3                 ! work cell

        Real*4          q2abs,                  ! abs(q_ij)
     *                  e4, Up_max,             ! work cells
     *                  Bound_Lo,               ! lower integration bound
     *                  Bound_Up,               ! upper integration bound
     *                  criter

        Integer         N_Sub,
     *                  i, j, l, m

        Complex*16      Spiller_Fourie, Trapez16
        External        Spiller_Fourie, Trapez16

c------------------------------------------------
        Complex*8       x0(N_Total_Max+1)
        Real*8          Z(N_Total_Max),         ! coordinates of interfaces
     *                  k0                      ! k0 = 2*pi/wave
        Real*4          Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Sigma_Accu(N_Total_Max),
     *                  Wave, UnitCoef(3),
     *                  qz, qx,                 ! (k2+k1)_z & (k2+k1)_x
     *                  CorreLength,            ! in-plane correlation length(A)
     *                  CorreVert,              ! vertical correlation length(A)
     *                  h_jagged,               ! jaggendness of the surface(0.:1.)
     *                  xabs,                   ! abs(x0_max)
     *                  TERang,                 ! critical TER angle
     *                  pi                      ! 3.1416
        Integer         i_Spec,
     *                  N_Top,                  ! top layer sublayers
     *                  N_Total,                ! total of layers
     *                  Index_corfu,            ! index of correlation function
     *                  MaxOrder(2),            ! max allowed S-matrix element (depricated)
     *                  iDebug,                 ! Debug flag
     *                  ifail                   ! failure code
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
     *                    i_Born,     ! approximation for wavefields
     *                    i_Fourier,  ! corr.integral: exp[K(x)]-1 or K(x)
     *                    isummariz,  ! add roughness-to-steps flag
     *                    Relief_Period,!surface relief period
     *                    max_Periods,! number of correlated periods
     *                    S_factor,   !==(Relief_Period/CorreLength)^2
     *                    A_Factor    ! normalization for periodic func.

        Integer    nSpil_max
        Parameter (nSpil_max = 2001)
        Real*8     SpilSum8_1,          ! "small"-Spiller integral
     *             SpilSum8_2,          ! -"-, 2nd correction
     *             SpilSum8_3           ! -"-, 3rd correction
        Real*4     x0_Sp, dx_Sp, xk_Sp,
     *             Sigma2_prime(N_Total_Max),
     *             CorLe2_prime(N_Total_Max),
     *             CorLen_prime(N_Total_Max),
     *             CorLen_min, CorLen_max,
     *             ximax, Sigma2sum, Trapint,
     *             Spiller_Array(nSpil_max)
        Integer    nCor, nSpil, Npts
        Common  /SpilDat/ SpilSum8_1, SpilSum8_2, SpilSum8_3,
     *                    Spiller_Array, Sigma2_prime,
     *                    CorLe2_prime, CorLen_prime,
     *                    CorLen_min, CorLen_max,
     *                    x0_Sp, dx_Sp, xk_Sp,
     *                    ximax, Sigma2sum, Trapint,
     *                    nCor, nSpil, Npts

        Complex*16 QQ                   ! q_i * q_j
        Real*4     Omega                ! = qx
        Common  /SpilInt/ QQ, Omega
c+==========================================================================+
c|The integral calculated by this program is:                               |
c|                                                                          |
c|        inf.                                                              |
c|      2  |               +                     +                          |
c|I = ---- | dx cos(q_x*x) |exp[q_ij*K_ij(x)] - 1|                       (1)|
c|    q_ij |               +                     +                          |
c|         0                                                                |
c|                                                                          |
c|where K_ij(x) is the Spiller-type correlation function:                   |
c|                                                                          |
c|             N                 +                 +                        |
c|K_ij(x) =   SUM  Sigma'_n^2 exp| - x^2 / xi'_n^2 |                     (2)|
c|         max(i,j)              +                 +                        |
c|                                                                          |
c|and:                                                                      |
c|                  Sigma_n^2                                               |
c|Sigma'_n^2 = ------------------ ;    xi'_n^2 = xi_n^2 + p_n^2          (3)|
c|             1 + p_n^2 / xi_n^2                                           |
c|                                                                          |
c|--------------------------------------------------------------------------|
c|At small                                                                  |
c|                    |q_ij*Sigma'_n^2| << 1,                            (4)|
c|                                                                          |
c|the exponent in (1) can be expanded and we arrive at:                     |
c|                                                                          |
c|                  inf.                                                    |
c|                   |                                                      |
c|I ~ K_ij(q_x) =  2 | dx cos(q_x*x) K_ij(x) =                              |
c|                   |                                                      |
c|                   0                                                      |
c|               N                         -xi'_n^2 * q_x^2 / 4             |
c|  = sqrt(pi)  SUM  Sigma'_n^2 * xi'_n * e                     =           |
c|            max(i,j)                                                      |
c|                                 +                              +         |
c|      N                          |sqrt(pi)  -(xi'_n * q_x)^2 / 4|         |
c|  =  SUM  Sigma'_n^2 * (2*xi'_n) |-------  e                    | =       |
c|   max(i,j)                      |    2                         |         |
c|                                 +                              +         |
c|     N                                                                    |
c|  = SUM  Sigma'_n^2 * (2*xi'_n) * ArrCos(xi'_n * q_x)                  (5)|
c|  max(i,j)                                                                |
c|                                                                          |
c| - because the array ArrCos just contains a half of a dimensionless       |
c| version of this integral, tabulated in CORRFUN2 (see subr. Spiller_Fill).|
c|--------------------------------------------------------------------------|
c|At large  |q_ij*Sigma'_n^2| ? 1, the integration in (1) must be carried   |
c|out numerically.  However, the integration is simplified due to the       |
c|expansion of the exponent at great                                        |
c|                   +                                                +     |
c|                   |          +                                    +|     |
c|x^2 < x_max^2 = max|xi'_n^2 ln|100*(N-max(i,j)+1)*|q_ij|*Sigma'_n^2||  (6)|
c|                   |          +                                    +|     |
c|                   +                                                +     |
c|Then:                                                                     |
c|                 x_max           +                                   +    |
c|                  ?              | 1  +                     +        |    |
c|I = K_ij(q_x) + 2 | dx cos(q_x*x)|----|exp[q_ij*K_ij(x)] - 1|-K_ij(x)| (7)|
c|                  ?              |q_ij+                     +        |    |
c|                  0              +                                   +    |
c|--------------------------------------------------------------------------|
c|A further development of this idea is the expansion of the exponent to    |
c|higher orders:                                                            |
c|                                                                          |
c| 1  + q_ij*K_ij(x) +             q_ij            q_ij^2                   |
c|----|e          - 1| = K_ij(x) + ----K_ij^2(x) + ------K_ij^3(x) + ...    |
c|q_ij+              +               2                6                     |
c|                                                                          |
c|Then:                                                                     |
c|              +     2       -q_x^2*xi_n^2/4+                              |
c|I=sqrt(pi)*SUM|Sigma * xi  e               | +                            |
c|            n +     n    n                 +                              |
c|                                                                          |
c|              +     2        2         -q_x^2*xi_nm^2/4+ q_ij             |
c| +sqrt(pi)*SUM|Sigma  * Sigma  * xi   e                | ---- +           |
c|           n,m+     n        m     nm                  +   2              |
c|                                                                          |
c|               +     2       2       2         -q_x^2*xi_nml^2/4+ q_ij^2  |
c| +sqrt(pi)*SUM*|Sigma * Sigma * Sigma * xi    e                 | ----- + |
c|          n,m,l+     n       m       l    nml                   +   6     |
c|                                                                          |
c|  x_max           +                                                       |
c|   |              | 1  + q_ij*K_ij(x) +                                   |
c| +2|dx cos(q_x*x)*|----|e          - 1| -                                 |
c|   |              |q_ij+              +                                   |
c|   0              +                                                     + |
c|                                         q_ij            q_ij^2         | |
c|                             - K_ij(x) - ----K_ij^2(x) - ------K_ij^3(x)| |
c|                                           2               6            | |
c|                                                                        + |
c|where:                                                                    |
c|                1          1        1                                     |
c|             -------  = ------ + ------                                   |
c|             xi_nm^2    xi_n^2   xi_m^2                                   |
c|                                                                          |
c|               1           1        1        1                            |
c|            --------  = ------ + ------ + ------                          |
c|            xi_nml^2    xi_n^2   xi_m^2   xi_l^2                          |
c|                                                                          |
c|The sums over n,m,l do not depend on (i,j) and can be tabulated for all   |
c|the roots (i,j) at one pair of interfaces (see subr.<Spiller_Fill>). The  |
c|numerical part of the integral is evaluated by subr.<Trapez> and the      |
c|integrand (the difference) is given by <Spiller_Fourie>.                  |
c+==========================================================================+
        Corr2    = (0.0D0, 0.0D0)

c The estimation for the integral at small
c |q_ij|*Summ(Sigma^2) (if (4) is satisfied):

        Corr2 = SpilSum8_1
        if (i_Fourier.eq.1)     return

c+---------------------------+
c| Estimation of the maximum |
c|  of correlation function  |
c| (checking condition (4)): |
c+---------------------------+

        q2abs  = SNGL(Abs(q_ij))
        criter = Sigma2sum * q2abs
        if (criter.lt.0.01)     return

        Sum16_2 = SpilSum8_2 * q_ij/2.
        Corr2   = Corr2 + Sum16_2
        if (Abs(Sum16_2).lt.0.05*Abs(Corr2))    return

        if (iDebug.eq.5)  Then  !------------------------+
        Sum16_3 = SpilSum8_3 * (q_ij**2)/6.             !|
        Corr2   = Corr2 + Sum16_3                       !|
        if (Abs(Sum16_3).lt.0.05*Abs(Corr2))    return  !|
        endif  !-----------------------------------------+

c##########################################################
c The calculations at great |q_ij|*Summ(Sigma^2).
c The integration interval is divided onto subintervals
c in order to avoid the oscillations of the integrand.
c##########################################################
c Upper bound of integration:
        Up_max = Min( CorLen_max * Sqrt(Log((1.0E6)*criter)),
     *                CorLen_min * 20. )

c The maximum number of sub-integrals:
        N_Sub = INT(Up_Max / Trapint + 1.5)
        QQ    = q_ij
        Omega = qx              ! that is from cos(omega*x)
        e4    = 0.01 / N_Sub

        if (iDebug.eq.5)  Then  !-------------------------------+
          write (3,*) ' Up_max/Corlen_max=',Up_max/CorLen_max, !|
     *                ' Points=',Int((Npts-1)*Up_max/Trapint), !|
     *                ' Inpol=',nSpil                          !|
          Open  (unit=12,file='bound.dat',status='unknown')    !|
          Open  (unit=14,file='func.dat',status='unknown')     !|
          m = 0                                                !|
        endif  !------------------------------------------------+

c       l = 0
c 1     l = l + 1   !<--------------------------------+
        do      l = 1,N_Sub   !=======================+
          Bound_Lo = Trapint * (l-1)                 !|
          Bound_Up = Trapint *   l                   !|
                                                     !|
          Sum16 = Trapez16 (Bound_Lo, Bound_Up,      !|
     *                      Spiller_Fourie, Npts)    !|
          Corr2  = Corr2 + Sum16                     !|
          if (Abs(Sum16).lt.e4*Abs(Corr2)) goto 20 !--+->+
c cccc    if (Bound_Up.gt. 20.*CorLen_min) goto 20 !--+->|
                                                     !|  |
          if (iDebug.eq.5)  Then  !-------------+     |  |
            m = m + 1                          !|     |  |
            write (12,52) i,j,m,               !|     |  |
     *                    Bound_Lo,Bound_Up,   !|     |  |
     *                    Abs(Real(Sum16)),    !|     |  |
     *                    Abs(Imag(Sum16)),    !|     |  |
     *                    Abs(Real(Corr2)),    !|     |  |
     *                    Abs(Imag(Corr2))     !|     |  |
  52        format(3i5,10g15.7)                !|     |  |
        endif  !--------------------------------+     |  |
                                                     !|  |
        enddo  !======================================+  |
c ccc   goto 1 !==================>===================+  |
                                                        !|
  20    continue  !<------------------<------------------+

        if (iDebug.eq.5)  Then  !--+
          Close (unit=12)         !|
          Close (unit=14)         !|
        endif  !-------------------+
        return
        end

c ####################################################### 6

        Complex*16 Function Spiller_Fourie (x)
c---------------------------------------------------------------
c          ** The Fourie-integrand corresponding **
c           ** to Spiller's correlation function **
c---------------------------------------------------------------
c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter       (N_Top_Max     = 150,
c cc *                   N_Total_Max   = N_Top_Max+1)

        Complex*16 w16                  ! work cell

        Real*4     x,                   ! the argument
     *             e4                   ! work cell

        Integer    iDebug

        Real*4     Finli, Spiller_Func
        External   Finli, Spiller_Func

        Integer    nSpil_max
        Parameter (nSpil_max = 2001)
        Real*8     SpilSum8_1,          ! "small"-Spiller integral
     *             SpilSum8_2,          ! -"-, 2nd correction
     *             SpilSum8_3           ! -"-, 3rd correction
        Real*4     x0_Sp, dx_Sp, xk_Sp,
     *             Sigma2_prime(N_Total_Max),
     *             CorLe2_prime(N_Total_Max),
     *             CorLen_prime(N_Total_Max),
     *             CorLen_min, CorLen_max,
     *             ximax, Sigma2sum, Trapint,
     *             Spiller_Array(nSpil_max)
        Integer    nCor, nSpil, Npts
        Common  /SpilDat/ SpilSum8_1, SpilSum8_2, SpilSum8_3,
     *                    Spiller_Array, Sigma2_prime,
     *                    CorLe2_prime, CorLen_prime,
     *                    CorLen_min, CorLen_max,
     *                    x0_Sp, dx_Sp, xk_Sp,
     *                    ximax, Sigma2sum, Trapint,
     *                    nCor, nSpil, Npts

        Complex*16 QQ                   ! q_i * q_j
        Real*4     Omega                ! = qx
        Common  /SpilInt/ QQ, Omega
c----------------------------------------------------------------

        if (x.lt.xk_Sp) Then  !---------------------------+
          e4 = Finli (x,Spiller_Array,x0_Sp,dx_Sp,nSpil) !|
        else  !-------------------------------------------|
          e4 = Spiller_Func(x)                           !|
        endif  !------------------------------------------+

        w16 = QQ*e4
        iDebug = 0
ccccc   iDebug = 5
        if (iDebug.eq.5)        Then  !-----------------------------+
c DEBUG case (improved precision, longer time) [see Spiller_Fill]   |
          if (Abs(w16).gt.1.D-4)  Then  !----------------+          |
            w16 = exp(w16)-1.-w16-w16*w16/2.-(w16**3)/6.!| CDexp(x) |
            w16 = w16/QQ                                !|          |
          else  !----------------------------------------| inf x^n  |
            w16 = e4*(w16**3)/24.                       !| Sum ---  |
          endif  !---------------------------------------+ n=0  n!  |
                                                                   !|
        else  !-----------------------------------------------------|
c NORMAL case (normal precision, shorter time) [see Spiller_Fill]   |
          if (Abs(w16).gt.1.D-4)  Then  !----------------+          |
            w16 = exp(w16)-1.-w16-w16*w16/2.            !| CDexp(x) |
            w16 = w16/QQ                                !|          |
          else  !----------------------------------------| inf x^n  |
            w16 = e4*w16*(w16/6.+ w16*w16/24.)          !| Sum ---  |
          endif  !---------------------------------------+ n=0  n!  |
                                                                   !|
        endif  !----------------------------------------------------+
c +=========================+
c |The factor "2" is due    |
c |to integrating over the  |
c |range (0,+infinity) in-  |
c |stead of integrating over|
c |(-infinity,+infinity):   |
c +=========================+
        w16 = 2. * w16 * cos(Omega*x)

        Spiller_Fourie = w16

cc      write   (14,54)  x, Abs(Real(w16)), Abs(Imag(w16))
cc54    format  (10g15.7)

        return
        end

c ####################################################### 7

        Real*4  Function Spiller_Func (x)
c--------------------------------------------------------
c ** Calculates Spiller's correlation function in one point **
c--------------------------------------------------------

c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter       (N_Top_Max     = 150,
c cc *                   N_Total_Max   = N_Top_Max+1)

        Real*8     Summ8           ! work cell

        Real*4     e4, x
        Integer    n

        Integer    nSpil_max
        Parameter (nSpil_max = 2001)
        Real*8     SpilSum8_1,          ! "small"-Spiller integral
     *             SpilSum8_2,          ! -"-, 2nd correction
     *             SpilSum8_3           ! -"-, 3rd correction
        Real*4     x0_Sp, dx_Sp, xk_Sp,
     *             Sigma2_prime(N_Total_Max),
     *             CorLe2_prime(N_Total_Max),
     *             CorLen_prime(N_Total_Max),
     *             CorLen_min, CorLen_max,
     *             ximax, Sigma2sum, Trapint,
     *             Spiller_Array(nSpil_max)
        Integer    nCor, nSpil, Npts
        Common  /SpilDat/ SpilSum8_1, SpilSum8_2, SpilSum8_3,
     *                    Spiller_Array, Sigma2_prime,
     *                    CorLe2_prime, CorLen_prime,
     *                    CorLen_min, CorLen_max,
     *                    x0_Sp, dx_Sp, xk_Sp,
     *                    ximax, Sigma2sum, Trapint,
     *                    nCor, nSpil, Npts

        Summ8 = 0.0D0
        do      n=1,nCor  !================+
          e4 = x*x / CorLe2_prime(n)      !|
                                          !|
c cc      if (e4.lt.707.) Then  !-----+    |for real*8
          if (e4.lt.87.)  Then  !-----+    |for real*4
            e4 = exp(-e4)            !|    |
          else  !---------------------|    |
            e4 = 0.                  !|    |
          endif  !--------------------+    |
                                          !|
          e4 = Sigma2_prime(n) * e4       !|
                                          !|
          Summ8 = Summ8 + e4              !|
        enddo  !===========================+
        Spiller_Func = SNGL(Summ8)

        return
        end

c ####################################################### 8

        Subroutine      Spiller_Fill (i,j,q2max)
c--------------------------------------------------------
c                 ** Fills Spiller_Array **
c ** (the array of Spiller's correlation function values) **
c--------------------------------------------------------

c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
c cc    Parameter       (N_Top_Max     = 150,
c cc *                   N_Total_Max   = N_Top_Max+1)

        Real*4     xi2_n, p2_n,
     *             e4, criter, x,
     *             q2max, Corlen
        Integer    i, j, l, m, n

        Real*4     Spiller_Func, Cos_Integral
        External   Spiller_Func, Cos_Integral

c------------------------------------------------
        Complex*8       x0(N_Total_Max+1)
        Real*8          Z(N_Total_Max),         ! coordinates of interfaces
     *                  k0                      ! k0 = 2*pi/wave
        Real*4          Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Sigma_Accu(N_Total_Max),
     *                  Wave, UnitCoef(3),
     *                  qz, qx,                 ! (k2+k1)_z & (k2+k1)_x
     *                  CorreLength,            ! in-plane correlation length(A)
     *                  CorreVert,              ! vertical correlation length(A)
     *                  h_jagged,               ! jaggendness of the surface(0.:1.)
     *                  xabs,                   ! abs(x0_max)
     *                  TERang,                 ! critical TER angle
     *                  pi                      ! 3.1416
        Integer         i_Spec,
     *                  N_Top,                  ! top layer sublayers
     *                  N_Total,                ! total of layers
     *                  Index_corfu,            ! index of correlation function
     *                  MaxOrder(2),            ! max allowed S-matrix element (depricated)
     *                  iDebug,                 ! Debug flag
     *                  ifail                   ! failure code
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
     *                    i_Born,     ! approximation for wavefields
     *                    i_Fourier,  ! corr.integral: exp[K(x)]-1 or K(x)
     *                    isummariz,  ! add roughness-to-steps flag
     *                    Relief_Period,!surface relief period
     *                    max_Periods,! number of correlated periods
     *                    S_factor,   !==(Relief_Period/CorreLength)^2
     *                    A_Factor    ! normalization for periodic func.

        Integer    nSpil_max
        Parameter (nSpil_max = 2001)
        Real*8     SpilSum8_1,          ! "small"-Spiller integral
     *             SpilSum8_2,          ! -"-, 2nd correction
     *             SpilSum8_3           ! -"-, 3rd correction
        Real*4     x0_Sp, dx_Sp, xk_Sp,
     *             Sigma2_prime(N_Total_Max),
     *             CorLe2_prime(N_Total_Max),
     *             CorLen_prime(N_Total_Max),
     *             CorLen_min, CorLen_max,
     *             ximax, Sigma2sum, Trapint,
     *             Spiller_Array(nSpil_max)
        Integer    nCor, nSpil, Npts
        Common  /SpilDat/ SpilSum8_1, SpilSum8_2, SpilSum8_3,
     *                    Spiller_Array, Sigma2_prime,
     *                    CorLe2_prime, CorLen_prime,
     *                    CorLen_min, CorLen_max,
     *                    x0_Sp, dx_Sp, xk_Sp,
     *                    ximax, Sigma2sum, Trapint,
     *                    nCor, nSpil, Npts
c--------------------------------------------------------

        xi2_n      = CorreLength**2
        SpilSum8_1 = 0.0D00
        SpilSum8_2 = 0.0D00
        SpilSum8_3 = 0.0D00
        Sigma2sum  = 0.0E00
        CorLen_min = 1.0E38
        CorLen_max = 0.0E00
        xk_Sp      = 0.
        m          = Max(i,j)
        nCor       = 0

        do     n = m, N_Total  !==================================+
          nCor = nCor + 1                                        !|
          p2_n = SNGL(4.*RelaxLength*((Z(n)-Z(i))+(Z(n)-Z(j))))  !|
          CorLe2_prime(nCor) = xi2_n + p2_n                      !|
          CorLen_prime(nCor) = Sqrt(CorLe2_prime(nCor))          !|
          Sigma2_prime(nCor) = Sigma(n)**2 / (1.+p2_n/xi2_n)     !|
                                                                 !|
          Sigma2sum = Sigma2sum + Sigma2_prime(nCor)             !|
                                                                 !|
          e4 = CorLen_prime(nCor)                                !|
          if (e4.lt.CorLen_min) CorLen_min = e4                  !|
          if (e4.gt.CorLen_max) CorLen_max = e4                  !|
                                                                 !|
c+-------------------------+                                      |
c|Calculation of integral  |                                      |
c|from Spiller at small qz:|                                      |
c+-------------------------+                                      |
          e4 = CorLen_prime(nCor) * Abs(qx)                      !|
                                                                 !|
          e4 = 2. * CorLen_prime(nCor)                           !|
     *       * Sigma2_prime(nCor)                                !|
     *       * Cos_Integral (e4)                                 !|
                                                                 !|
          SpilSum8_1 = SpilSum8_1 + e4                           !|
        enddo  !==================================================+
        if (i_Fourier.eq.1)     return

c Criterium of small sigma*q:
        criter = Sigma2sum * q2max
ccccc   if (criter.lt.0.01)     return          ! 1%: like in GIDS
        if (criter.lt.0.05)     return          ! 5%: corr. 04-10-96

c 2nd order corrections:
        do     n = 1, nCor  !============================+
        do     m = n, nCor  !==========================+ |
                                                      !| |
          Corlen = 1./CorLe2_prime(n)                 !| |
     +           + 1./CorLe2_prime(m)                 !| |
          Corlen = 1./Sqrt(CorLen)                    !| |
          e4     = CorLen * Abs(qx)                   !| |
                                                      !| |
          e4 = 2.*Corlen                              !| |
     *       * Sigma2_prime(n)                        !| |
     *       * Sigma2_prime(m)                        !| |
     *       * Cos_Integral (e4)                      !| |
                                                      !| |
          if (n.ne.m)   e4 = 2*e4                     !| |
          SpilSum8_2 = SpilSum8_2 + e4                !| |
        enddo  !=======================================+ |
        enddo  !=========================================+

        if (iDebug.eq.5) Then  !-----------------------------+
c 3rd order corrections:                                     |
        do     n = 1, nCor  !==============================+ |
        do     m = n, nCor  !============================+ | |
        do     l = m, nCor  !==========================+ | | |
                                                      !| | | |
          Corlen = 1./CorLe2_prime(n)                 !| | | |
     +           + 1./CorLe2_prime(m)                 !| | | |
     +           + 1./CorLe2_prime(l)                 !| | | |
          Corlen = 1./Sqrt(CorLen)                    !| | | |
          e4     = CorLen * Abs(qx)                   !| | | |
                                                      !| | | |
          e4 = 2.*Corlen                              !| | | |
     *       * Sigma2_prime(n)                        !| | | |
     *       * Sigma2_prime(m)                        !| | | |
     *       * Sigma2_prime(l)                        !| | | |
     *       * Cos_Integral (e4)                      !| | | |
                                                      !| | | |
          if (n.ne.m .OR.  n.ne.l)  e4 = 3*e4         !| | | |
          if (n.ne.m .AND. n.ne.l)  e4 = 2*e4         !| | | |
          SpilSum8_3 = SpilSum8_3 + e4                !| | | |
        enddo  !=======================================+ | | |
        enddo  !=========================================+ | |
        enddo  !===========================================+ |
        endif  !---------------------------------------------+

c+---------------------------------------------------+
c|The size of 1 integration subinterval;             |
c|no greater than 1/2 oscillations are               |
c|allowed (q*x < pi) in one subinterval:             |
                                                    !|
        Trapint = Min (CorLen_min,                  !|
     *                 pi/Max(Abs(qx),1.E-10))      !|
                                                    !|
c The number of function calls in one integration    |
c (it is optimum to tabulate Spiller in these points)|
                                                    !|
        Npts  = 31                                  !|
c+---------------------------------------------------+

c+------------------------+
c|Tabulating the Spiller's|
c|correlation function    |
c+------------------------+
        dx_Sp = Trapint / (Npts-1)
        x0_Sp = 0.

        ximax = CorLen_max * Sqrt(Log((1.0E6)*criter))
        nSpil = INT(ximax / dx_Sp + 1.5)
        if (nSpil.gt.nSpil_max)
     *  nSpil = nSpil_max
        xk_Sp = dx_Sp*(nSpil-1)

        do      n=1,nSpil  !=======================+
          x = x0_Sp + dx_Sp*(n-1)                 !|
          Spiller_Array(n) = Spiller_Func(x)      !|
        enddo  !===================================+
        return
        end

c ####################################################### 9

        Complex*16 Function Trapez16 (a,b,f,n)
c-----------------------------------------------------
c This subroutine integrated function f in interval [a,b]
c using the trapezium method.
c n - number of interval [a,b] subdivisions during integration.
c-----------------------------------------------------
        Complex*16      f,s
        Real*4          a,b,x,dx
        Integer         n,i
        External        f

        if (a.ge.b) Then  !-----------+
          Trapez16 = (0.0D0, 0.0D0)  !|
          return                     !|
        endif  !----------------------+

        dx = (b-a)/(n-1)
        s  = 0.5*(f(a)+f(b))
        do      i=2,n-1  !====+
          x = a + dx*(i-1)   !|
          s = s + f(x)       !|
        enddo  !==============+
        Trapez16 = s * dx
        return
        end

c ####################################################### 10

        Complex*16 Function  Corr3 (Sigma_ij, q_ij, qz0, qx,
     *                              CorreLength, SurMiscut,
     *                              TerraceSpre, StepHeight,
     *                              Index_corfu,isummariz)
c+===================================================+
c| ** The correlation integral of interface STEPS ** |
c+===================================================+
        Complex*16 q_ij,           ! q_i*q_j
     *             qs2,            ! (q_ij)*(Sigma_ij)
     *             qz,             ! SQRT(q_ij)
     *             dvdr            ! work cell
        Complex*16 Im/(0.,1.)/,    ! i
     *             P, H            ! Fouries of probabilities
        Real*8     px              ! CorreLength*qx
        Real*4     qz0, qx,        ! (k2+k1)_z & (k2+k1)_x
     *             CorreLength,    ! in-plane correlation length(A)
     *             SurMiscut,      ! surface miscut
     *             TerraceSpre,    ! terrace size spread
     *             StepHeight,     ! step height
     *             Sigma_ij,       ! Sigma_i * Sigma_j
     *             L,              ! == CorreLength
     *             S,              ! == TerraceSpre
     *             d,              ! CorreLength*Theta - height of steps
     *             Anm/2./,        ! coefficient,
     *             Corr_r,
     *             s2
        Integer    Index_corfu,
     *             isummariz

        Real*4     Cos_Integral
        External   Cos_Integral
c+====================================================================+
c|The correlation integral is:                                        |
c|                                                                    |
c|    +infinity (            2        )                               |
c|  1     |     |  (qz*sigma) C(x)    |  -iqx*x        2              |
c| ----   | dx  |e                 - 1| e       ~ sigma *C(qx)     (1)|
c| qz^2   |     |                     |         ~                     |
c|    -infinity (                     )                               |
c|                                                                    |
c|where qz^2 = q_i * q^{*}_j ;  and  sigma^2 = sigma_i * sigma_j .    |
c|                                                                    |
c|If C(x)=exp(-|x|/<L>), then C(qx)=2<L>/(1+qx^2<L>^2) ,              |
c|If C(x)=exp(-x^2/<L>^2), then C(qx)=sqrt(pi)*<L>*exp(-qx^2<L>^2/4) .|
c|--------------------------------------------------------------------|
c|We take C(qx) in the form by P.R.Pukite, C.S.Lent and P.I.Cohen --  |
c|Surf.Sci. v.161 p.39-68 (1985):                                     |
c|                   (                     )                          |
c|         2*Anm     | [1-P(qx)]*[1-H(qz)] |                          |
c| C(qx) = ------- Re| ------------------- |                       (2)|
c|         qx^2<L>   |    1-P(qx)*H(qz)    |                          |
c|                   (                     )                          |
c|where <L>=CorreLength is the mean terrace length, and:              |
c|                                                                    |
c|                                                 infty              |
c|      (         -iqx*L        infty       -iqz*hd  (         -iqz*s |
c| P(x)=|dL P(L) e     ,  H(qz)= SUM  H(h) e       * |ds H(s) e     , |
c|      )                      h=-infty              )                |
c|                                                -infty              |
c|                                                                    |
c|-- are the Fourier transforms of a probability per unit length of a |
c|terrace of length L, and a probability distribution for a step of   |
c|height s=hd, respectively (d is an elementary atomic step).         |
c|--------------------------------------------------------------------|
c|We take P in the form corresponding to a geometric staircase:       |
c|                                                                    |
c|         1   -x/<L>                              1                  |
c| P(x) = --- e               ==>      P(qx) = ----------          (3)|
c|        <L>                                  1+i*qx*<L>             |
c|--------------------------------------------------------------------|
c|For H two forms are suggested:                                      |
c|-------------+                                                      |
c|Index_corfu=7|                                                      |
c|-------------+                      -iqz*hd                  2      |
c| H(h) = delta_{h,1}  ,     H(qz) = e       ~ 1-iqz*hd-(qz*hd) /2    |
c|                                           ~                        |
c|The expansion (1) assumes qz*sigma << 1, therefore the cosine in (2)|
c|can be expanded. We take (qz*sigma)=SQRT(q_i * q^{*}_j * sigma_ij^2)|
c|                                                                    |
c|              Anm*<L>*(qz*d)^2                                      |
c| C(qx) = ---------------------------   ,                        (4a)|
c|         (qx*<L>+qz*d)^2 + (qz*d)^4/4                               |
c|                                                                    |
c|where d=sigma=Theta*<L>, and Theta is the miscut angle.             |
c|-------------+ - - - - - - - - - - - - - - - - - - - - - - - - - - -|
c|Index_corfu=8|                           2     2                    |
c|-------------+          -iqz*Theta*<L>-qz sigma /4                  |
c|               H(qz) = e                                            |
c|                                                                    |
c|Under the same approximations as above we obtain:                   |
c|                                       2                            |
c|                 Anm*<L>*(qz*sigma_eff) /2                          |
c| C(qx) = ---------------------------------------------   ,      (4b)|
c|         (qx*<L>+qz*Theta*<L>)^2 + (qz*sigma_eff)^4/16              |
c|-------------+ - - - - - - - - - - - - - - - - - - - - - - - - - - -|
c|Index_corfu=9|                                                      |
c|-------------+      2 2                                  2     2    |
c|        -i*qx*<L>-qx S /4              -i*qz*Theta*<L>-qz sigma /4  |
c| P(qx)=e                 ,      H(qz)=e                             |
c|                                                                    |
c|where S=TerraceSpre - spread of terrace sizes (see [Rabedeau, Tids- |
c|well, Pershan, Bevk, and Freer -- Appl.Phys.Lett. 59 (1991) 3422]). |
c|--------------------------------------------------------------------|
c|According to Sinha et.al.[Physica B, 198 (1994) 72-77], the correla-|
c|tion function for steps at  qx*<L> << 1  must transform  to  a self-|
c|affine function. At qx*<L> << 1 and small sigma^2 << 2*(Theta*<L>)^2|
c|the denominator is reduced to (qz*Theta*sigma)^2 and C(qx)=Anm*<L>. |
c|At the same conditions the correlation function                     |
c|                    C(qx)=2<L>/(1+(qx<L>)^2),                       |
c|and the correlation function of steps transforms to 1/2 of that for |
c|self-affine roughness. Thus, we take Anm=2.                         |
c+====================================================================+
ccc???  qz    = SQRT(q_ij)      !Provides negative intensity!!!
        qz    = qz0
        L     = CorreLength
        d     = SurMiscut*L     !Surface shift at <L> scale
        S     = TerraceSpre
        Corr3 = (0.0D0, 0.0D0)

        if     (Index_corfu.eq.7)  Then !----------+i=7 (exact)
                                                  !|
           P = 1./(1.+Im*qx*L)                    !|
           H = exp(-Im*qz*d)                      !|Sharp peaks: no smoothing!
           if (abs(qx).gt.1.E-32)                 !|
     *           Corr3 = (2.*Anm/(qx*qx*L))       !|
     *                 * Real( (1.-P)*(1.-H)      !|
     /                 /          (1.-P*H) )      !|
                                                  !|
        elseif (Index_corfu.eq.8)  Then !----------|i=8 (exact)
                                                  !|
           P = 1./(1.+Im*qx*L)                    !|
           H = exp(-Im*qz*d-qz*qz*Sigma_ij/4.)    !|Broad peaks: smoothing!
           if (abs(qx).gt.1.E-32)                 !|
     *           Corr3 = (2.*Anm/(qx*qx*L))       !|
     *                 * Real( (1.-P)*(1.-H)      !|
     /                 /          (1.-P*H) )      !|
                                                  !|
        elseif (Index_corfu.eq.9)  Then !----------|i=9 (Pershan)
                                                  !|
           P = exp(-Im*qx*L-qx*qx*S*S/4.)         !|
           H = exp(-Im*qz*d-qz*qz*Sigma_ij/4.)    !|Broad peaks: smoothing!
           if (abs(qx).gt.1.E-32)                 !|
     *           Corr3 = (2.*Anm/(qx*qx*L))       !|
     *                 * Real( (1.-P)*(1.-H)      !|
     /                 /          (1.-P*H) )      !|
                                                  !|
        elseif (Index_corfu.eq.57) Then !----------|i=57 (=7a, approx)
                                                  !|
ccc???+-----------------------+                   !|
ccc???|    d    = StepHeight  |                   !|
ccc???|    s2   = Sigma_ij    |                   !|
ccc???+-----------------------+                   !|
           s2   = 2.*d*d                          !|Sharp peaks: no smoothing!
           px   = Cdabs(qx*L + qz*d)              !|
           qs2  = qz*qz*s2                        !|
           dvdr = px*px + qs2*qs2/16.             !|
           if (Abs(dvdr).gt.0.)                   !|
     *           Corr3 = Anm*L*(qs2/2) / dvdr     !|
                                                  !|
        elseif (Index_corfu.eq.58) Then !----------|i=58 (=8a, approx)
                                                  !|
           s2   = 2.*d*d + Sigma_ij               !|Broad peaks: smoothing!
           px   = Cdabs(qx*L + qz*d)              !|
           qs2  = qz*qz*s2                        !|
           dvdr = px*px + qs2*qs2/16.             !|
           if (Abs(dvdr).gt.0.)                   !|
     *           Corr3 = Anm*L*(qs2/2) / dvdr     !|
                                                  !|
        else   !-----------------------------------|
                                                  !|
           qz    = SQRT(q_ij)                     !|to suppress warnings!!!
                                                  !|
        endif  !-----------------------------------+

cccc    Corr3 = L
        Corr3 = (StepHeight**2) * Corr3

        if (isummariz.ne.0)     Then  !--------------+
c+======================================+            |
c|The correlation function is a summ of |            |
c|staircase and roughness contributions:|            |
c+======================================+            |
cccc For h=0.5                                       |
ccc       Corr_r = 2.*L*Sigma_ij/(1.+(qx*L)**2)     !|
          Corr_r = 2.*L*Sigma_ij*Cos_Integral(qx*L) !|
                                                    !|
          Corr3  = Corr3 + Corr_r                   !|
        endif  !-------------------------------------+

        return
        end
