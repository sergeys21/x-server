c  Transition layer processing for:  gid_sl32.FOR
c                                    ter_sl10.FOR
c                                    ter_sl12.FOR
c        --  Recursive VERSION --

        Subroutine Transi_r   (Thickness_Tr, Sigma,
     *                         TranLaw, N_Tran,
     *                         x0_1, x0_2,
     *                         xh_1, xh_2,
     *                         Pol, xabs, TERang,
     *                         g0, gh, gi, Daa, k0,
     *                         N_Roots, N_Fields,
     *                         N_Half,  j_Used,
     *                         S, B, B1, u, u1,
     *                         F, X, TTT,
     *                         M_rr, M_rt,
     *                         M_tr, M_tt,
     *                         W_rr, W_rt,
     *                         W_tr, W_tt,
     *                         Upper_Limit, Kod,
     *                         ifail )

        Integer    N_Roots,
     *             N_Half,
     *             N_Fields,
     *             N_Tran,
     *             j_Used,
     *             Kod,
     *             ifail,
     *             j

        Complex*16 x0n,xhn(2),           ! normalized x0,xh
     *             S(N_Roots,N_Roots),   ! S_matrix 4x4
     *             B(N_Roots,N_Roots),   ! 1/S_matrix
     *             B1(N_Roots,N_Roots),  ! 1/S_matrix (previous)
     *             u(N_Roots),           ! dispersion eqn.roots
     *             u1(N_Roots),          ! dispersion eqn.roots (previous)
     *             F(N_Roots),           ! F-matrix
     *             TTT(N_Roots,N_Roots), ! T1*...Tn*Ss*exp(ikuz)
     *             X(N_Roots,N_Roots),   ! X-matrix of interface
     *             M_rr(N_Half,N_Half),  ! product of operators
     *             M_rt(N_Half,N_Half),  ! product of operators
     *             M_tr(N_Half,N_Half),  ! product of operators
     *             M_tt(N_Half,N_Half),  ! product of operators
     *             W_rr(N_Half,N_Half),  ! product of operators
     *             W_rt(N_Half,N_Half),  ! product of operators
     *             W_tr(N_Half,N_Half),  ! product of operators
     *             W_tt(N_Half,N_Half),  ! product of operators
     *             gh,                   ! normalized exit angle
     *             Sqrt_Eps              ! Spec:        1   for Pol=1(Pol=0)
     *                                   !       Sqrt(1+x0) for Pol=2

        Complex*8  x0_1, x0_2,
     *             xh_1, xh_2

        Real*8     g0, gi,         ! normalized angles
     *             k0,             ! k0 = 2*pi/wave
     *             Upper_Limit     ! max matr.element (lim)

        Real*4     Thickness_Tr,   ! Thickness/2 of transition layer
     *             T_tr,           ! thickness of 1 trans. sublayer
     *             TranLaw(N_Tran),! law of transition
     *             Sigma,          ! r.m.s. roughness
     *             xabs,           ! abs(x0_max)
     *             TERang,         ! critical TER angle
     *             Pol,            ! 1. or 2.
     *             Daa,            ! Da/a
     *             c1, c2

        T_Tr   = 2.*Thickness_Tr / N_Tran

        do    j=1,N_Tran    !========================+
          j_Used = j_Used+1                         !|
          Call  Copy_Layer (B,B1,u,u1,N_Fields)     !|
          c1     = TranLaw(j)                       !|
          c2     = 1.-c1                            !|
          x0n    = (c1*x0_1 + c2*x0_2) / xabs       !|
          if (N_Fields .eq. 2) Then  !-----------+   |SPEC
            if (abs(Pol-2.).gt.1.E-20) Then !-+  |   |
              Sqrt_Eps = (1.,0.)             !|  |   |
            else !----------------------------+  |   |
              Sqrt_Eps = Sqrt(1.+x0n*xabs)   !|  |   |
            endif  !--------------------------+  |   |
            Call Process_Layer(g0, Sqrt_Eps,    !|   |
     *                         gi, Daa,         !|   |
     *                         x0n, xhn, k0,    !|   |
     *                         T_Tr,            !|   |
     *                         TERang, S, B,    !|   |
     *                         u,               !|   |
     *                         F,               !|   |
     *                         N_Fields, ifail) !|   |
          else !---------------------------------|   |DIFR
            xhn(1) = (c1*xh_1 + c2*xh_2)        !|   |
     *             * Pol / xabs                 !|   |
            xhn(2) = xhn(1)                     !|   |OLD!!! NO ACCOUNT for
          Call  Process_Layer (g0, gh,          !|   |non-centrosymmetric
     *                         gi, Daa,         !|   |crystals like in gid_sl97!!!
     *                         x0n, xhn, k0,    !|   |
     *                         T_Tr,            !|   |
     *                         TERang, S, B,    !|   |
     *                         u,               !|   |
     *                         F,               !|   |
     *                         N_Fields, ifail) !|   |
          endif  !-------------------------------+   |
          if (ifail.ne.0)     goto 100  !------------+-+
          Call  Process_Interface (k0,TERang,       !| v
     *                         Sigma,               !|
     *                         B1, S,               !|
     *                         u1,                  !|
     *                         u,                   !|
     *                         F,                   !|
     *                         X, TTT,              !|
     *                         N_Fields,            !|
     *                         Upper_Limit,         !|
     *                         Kod, ifail)          !|
            if (ifail.ne.0)     goto 100  !----------+-+
            Call Process_Interface2(X, F,           !| v
     *                         M_rr,                !|
     *                         M_rt,                !|
     *                         M_tr,                !|
     *                         M_tt,                !|
     *                         W_rr,                !|
     *                         W_rt,                !|
     *                         W_tr,                !|
     *                         W_tt,                !|
     *                         N_Fields, j_Used,    !|
     *                         ifail)               !|
          if (ifail.ne.0)     goto 100  !------------+-+
        enddo  !=====================================+ v

  100   return
        end
