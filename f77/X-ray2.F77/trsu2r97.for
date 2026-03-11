c
c *** Version 1.0 ***  for TRDS
c                   -- Interface roughness profile
c                      (in dependence on the interface number)
c                   -- 6 models of roughness correlations
c
c--------------------------------------------------------
c                         +------------------+
c ATTENTION (12-10-95) -- |RECURSION VERSION |
c                         +==================+
c ATTENTION (27-09-95) -- based on common unified
c                         subroutines for gid_sl32
c                                         gid_sl50
c                                         gid_sl55
c                                         gid_sl97
c                                         gids_10
c                                         ter_sl10
c                                         ter_sl12
c                                         ter_sl97
c                                         trds_10
c                                         trds_97
c +=====================================================+
c |   These routines compute x-ray wavefields under     |
c |        grazing incidence specular reflection        |
c | The calculations are carried out in 1 angular point |
c |-----------------------------------------------------|
c |Contents:                                            |
c | 1. Subroutine:            SFM                       |
c |-----------------------------------------------------|
c |Subroutine calls: RFMs_UNI.for -> NAG2(C02ADF,F04ADF)|
c |                  FLDs_UNI.for                       |
c +=====================================================+

c ####################################################### 1

        Subroutine      SFM (f0, Rc, N_Fields, N_Used, E, u,
     *                       ipol, Mode_Sigma, Compute_Fields)
c  ** Total_External_Reflex_from_Multi **
c--------------------------------------------------------
c  Mode_Sigma = 0   do NOT account for the roughness when
c                   computing the coherent wavefields
c  Mode_Sigma = 1   DO account for the roughness when
c                   computing the coherent wavefields
c
c ** The calculations have shown that there is no any
c    noticeable difference for diffuse scattering.
c--------------------------------------------------------
c+===============================+
        Include 'trds_97.inc'   !|
c+===============================+
c The above INC file contains smth like:
cc      Parameter  (N_Top_Max     = 150,
cc   *              N_Total_Max   = N_Top_Max+1)

        Integer    N_Roots, N_Half
        Parameter (N_Roots = 4,
     *             N_Half  = N_Roots/2)
        Real       Daa
        Parameter (Daa     = 0.)                    !da/a (not used)

        Complex*16 E(N_Roots,N_Total_Max+1),        ! Fields in layers
     *             u(N_Roots,N_Total_Max+1),        ! dispersion eqn.roots
     *             F(N_Roots,N_Total_Max+1),        ! F-matrix
     *             X(N_Roots,N_Roots),              ! S_{k}/S_{k-1}
     *             M_rr(N_Half,N_Half,N_Total_Max), ! product of operators
     *             M_rt(N_Half,N_Half,N_Total_Max), ! product of operators
     *             M_tr(N_Half,N_Half),             ! product of operators
     *             M_tt(N_Half,N_Half),             ! product of operators
     *             W_rr(N_Half,N_Half),             ! product of operators
     *             W_rt(N_Half,N_Half),             ! product of operators
     *             W_tr(N_Half,N_Half,N_Total_Max), ! product of operators
     *             W_tt(N_Half,N_Half,N_Total_Max), ! product of operators
     *             x0n,                             ! normalized x0
     *             xhn(2)/(0.0D0,0.0D0),
     *                    (0.0D0,0.0D0)/,           ! normalized xh (not used)
     *             Sqrt_Eps,                        ! sigma=1 pi=Sqrt(1+x0)
     *             S(N_Roots,N_Roots),              ! S_matrix 2x2
     *             B(N_Roots,N_Roots),              ! 1/S_matrix
     *             B1(N_Roots,N_Roots),             ! 1/S_matrix (previous)
     *             TTT(N_Roots,N_Roots),            ! T1*...Tn*Ss*exp(ikuz)
     *             Es                               ! vacuum wave amplitude

        Real*8     g0,                              ! normalized angle
     *             gi/0.0D0/,                       !    -"- miscut angle (not used)
     *             Upper_Limit                      ! max matr.element (lim)

        Real*4     f0,                              ! incidence angle
     *             Sgm_i,                           ! Sigma of interface
     *             Rc

        Integer    N_Used,                          ! number of used layers
     *             N_Plus,                          ! # roots with Im(u)>0
     *             N_Fields,                        ! # of wave fields
     *             i, j, k, n, Kod,                 ! work sell
     *             mode_Sigma,
     *             ipol

        Logical    Compute_Fields                   ! for roughness

c------------------------------------------------
        Complex*8       x0(N_Total_Max+1)
        Real*8          Z(N_Total_Max),             ! coordinates of interfaces
     *                  k0                          ! k0 = 2*pi/wave
        Real*4          Thickness(N_Total_Max),
     *                  Sigma(N_Total_Max),
     *                  Sigma_Accu(N_Total_Max),
     *                  Wave, UnitCoef(3),
     *                  qz, qx,                     ! (k2+k1)_z & (k2+k1)_x
     *                  CorreLength,                ! in-plane correlation length(A)
     *                  CorreVert,                  ! vertical correlation length(A)
     *                  h_jagged,                   ! jaggendness of the surface(0.:1.)
     *                  xabs,                       ! abs(x0_max)
     *                  TERang,                     ! critical TER angle
     *                  pi                          ! 3.1416
        Integer         i_Spec,
     *                  N_Top,                      ! top layer sublayers
     *                  N_Total,                    ! total of layers
     *                  Index_corfu,                ! index of correlation function
     *                  MaxOrder(2),                ! max allowed S-matrix element (depricated)
     *                  iDebug,                     ! Debug flag
     *                  ifail                       ! failure code
        Common  /TRDSdat/ Wave, UnitCoef, qz, qx, x0, Z,
     *                    Thickness, Sigma, Sigma_Accu,
     *                    CorreLength, CorreVert, h_jagged,
     *                    xabs, TERang, pi, k0, i_Spec,
     *                    N_Top, N_Total, Index_corfu,
     *                    MaxOrder, iDebug, ifail
c--------------------------------------------------------
        g0     = f0 * UnitCoef(1)/TERang

        Rc     = 0.
        ifail  = 0
        N_Used = 0
        Sgm_i  = 0.
        N_Plus = N_Fields / 2

        if (.NOT.Compute_Fields)  Then  !----+
          Upper_Limit = 10.**MaxOrder(1)    !| 14
        else  !------------------------------|
          Upper_Limit = 10.**MaxOrder(2)    !|  7
        endif  !-----------------------------+

c At f0=0 vacuum matrix becomes singular and cannot be inverted:
c                      |1 1|
c                      |0 0|
c Then, reflection coefficient cannot be determined and we set it to 1:
        do      i=1,N_Fields  !=======+
        do      k=1,N_Total+1 !=====+ |
          E(i,k) = (0.D0, 0.D0)    !| |
        enddo  !====================+ |
        enddo  !======================+

        if (Abs(g0).lt.1.0D-7)  g0 = 1.0D-7

        if (g0.lt.0.)      Then  !-----+
          Rc = 1.                     !|
          Return                      !|
        endif  !-----------------------+

        Kod = 1                         !Don't compute TTT
        Call    M1 (N_Roots,TTT,N_Fields)

c
c                     Process Vacuum
c
c+-------------------------------------------------------+
        x0n        = x0(N_Used+1) / xabs                !|
        if (ipol.eq.1) Then  !-----------+sigma          |
          Sqrt_Eps = (1.,0.)            !|               |
        else !---------------------------|pi             |
          Sqrt_Eps = Sqrt(1.+x0n*xabs)  !|               |
        endif  !-------------------------+               |
        Call    Process_Layer (g0, Sqrt_Eps,            !|see "rfms_uni.for"
     *                         gi, Daa,                 !|not used
     *                         x0n, xhn, k0,            !|
     *                         0.,                      !|Thc(0)
     *                         TERang, S, B,            !|
     *                         u(1,N_Used+1),           !|
     *                         F(1,N_Used+1),           !|
     *                         N_Fields, ifail)         !|
        if (ifail.ne.0) goto 100                        !|
c+-------------------------------------------------------+

c
c                    Process Top Layer
c
        do    i=1,N_Total !==============================+
                                                        !|
          Call  MtoM (N_Roots, N_Roots, B, B1, N_Fields)!|
          N_Used       = N_Used+1                       !|
                                                        !|
          x0n          = x0(N_Used+1) / xabs            !|
          if (ipol.eq.1) Then  !-----------+sigma        |
            Sqrt_Eps = (1.,0.)            !|             |
          else !---------------------------|pi           |
            Sqrt_Eps = Sqrt(1.+x0n*xabs)  !|             |
          endif  !-------------------------+             |
          Call  Process_Layer (g0, Sqrt_Eps,            !|see "rfms_uni.for"
     *                         gi, Daa,                 !|not used
     *                         x0n, xhn, k0,            !|
     *                         Thickness(N_Used),       !|
     *                         TERang, S, B,            !|
     *                         u(1,N_Used+1),           !|
     *                         F(1,N_Used+1),           !|
     *                         N_Fields,ifail)          !|
          if (ifail.ne.0)     goto 100                  !|
          if (mode_Sigma.eq.1) Sgm_i = Sigma(N_Used)    !|
          Call Process_Interface (k0, TERang, Sgm_i,    !|see "rfms_uni.for"
     *                            B1, S,                !|
     *                            u(1,N_Used),          !|
     *                            u(1,N_Used+1),        !|
     *                            F(1,N_Used+1),        !|
     *                            X, TTT,               !|
     *                            N_Fields,             !|
     *                            Upper_Limit,          !|
     *                            Kod, ifail)           !|
          if (ifail.ne.0)     goto 100  !----------------+--+
          Call Process_Interface2(X, F(1,N_Used+1),     !|  v
     *                            M_rr(1,1,N_Used),     !|
     *                            M_rt(1,1,N_Used),     !|
     *                            M_tr,                 !|
     *                            M_tt,                 !|
     *                            W_rr,                 !|
     *                            W_rt,                 !|
     *                            W_tr(1,1,N_Used),     !|
     *                            W_tt(1,1,N_Used),     !|
     *                            N_Fields, N_Used,     !|
     *                            ifail)                !|
          if (ifail.ne.0)     goto 100  !----------------+--+
          if (i.lt.N_Total) Then  !-----------+          |  v
            Call Backup2(W_tr(1,1,N_Used),   !|          |
     *                   W_tt(1,1,N_Used),   !|          |
     *                   W_tr(1,1,N_Used+1), !|          |
     *                   W_tt(1,1,N_Used+1), !|          |
     *                   N_Half,N_Plus)      !|          |
          endif  !----------------------------+          |
        enddo  !=========================================+

c
c                   Find solution Es
c

        Es = W_rt(1,1)
        Rc = SNGL((Abs(Es))**2)

        if (.NOT.Compute_Fields)        Return

c
c                  Compute the wavefields
c

        E(1,1)        = (1.0D0, 0.0D0)
        E(2,1)        = Es
        E(1,N_Used+1) = W_tt(1,1,N_Used)
cc   *                * E(1,1)
cc   /                / F(1,N_Used+1)
        E(2,N_Used+1) = (0.0D0, 0.0D0)

        N = N_Plus + 1
        do      k = N_Used-1,1,-1   !========+
          Call  Find_Fields (E(1,1),        !|T_0
     *                       E(N,k+2),      !|R_{k+1}
     *                       E(1,k+1),      !|T_k
     *                       E(N,k+1),      !|R_k
     *                       W_tt(1,1,k),   !|
     *                       W_tr(1,1,k),   !|
     *                       M_rt(1,1,k+1), !|
     *                       M_rr(1,1,k+1), !|
     *                       N_Plus,ifail)  !|
          if (ifail.ne.0)       goto 100  !--+-+
        enddo  !=============================+ v

        do      k = 2, N_Used   !============+
          Call  Norm_Fields (E(1,k),        !|T_k
     *                       E(N,k),        !|R_k
     *                       F(1,k),        !|Ft_k
     *                       F(N,k),        !|Fr_k
     *                       N_Plus)        !|
        enddo  !=============================+

        if (iDebug.ne.0)        Then !------------------+
          write (3,50)  g0                             !|
  50      format('The wavefields for g0/TERang=',g9.3) !|
                                                       !|
          do    k = 1, N_Used  !=============+          |
            if (k.lt.N_Used)  Then !---+     |          |
              j = N_Fields            !|     |          |
            else  !--------------------|     |          |
              j = N_Plus              !|     |          |
            endif  !-------------------+     |          |
            write (3,55) k,(Abs(E(i,k+1)),  !|          | !!!!!!!!!!! DEBUG
     *                   i=1,j)             !|          | !!!!!!!!!!! DEBUG
  55        format(i5,4g12.5)               !|          | !!!!!!!!!!! DEBUG
          enddo  !===========================+          |
        endif  !----------------------------------------+

        Return

  100   continue
        Return
        end
