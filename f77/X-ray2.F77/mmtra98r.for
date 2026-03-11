c ----------------------------------------------------------------------
c  Transition layer processing for:  mag_sl98.FOR
c                                    mag_sl99.FOR v.1.1b
c                        -- i.e. for magnetic resonance of HARD x-rays.
c
c        --  Recursive VERSION --
c ----------------------------------------------------------------------

        Subroutine Transi_rM98  (Thickness_Tr, Sigma,
     *                           TranLaw, N_Tran,
     *                           x_ij, x1_ij, xabs,
     *                           TERang, g0, k0,
     *                           Magnetic,
     *                           N_Roots, N_Fields,
     *                           N_Half,  j_Used,
     *                           S, B, B1, u, u1,
     *                           F, X, TTT,
     *                           M_rr, M_rt,
     *                           M_tr, M_tt,
     *                           W_rr, W_rt,
     *                           W_tr, W_tt,
     *                           Upper_Limit, Kod,
     *                           ifail )

        Integer    N_Roots,
     *             N_Half,
     *             N_Fields,
     *             N_Tran,
     *             Magnetic,
     *             j_Used,
     *             Kod, ifail

        Complex*16 x_ij(2,2),x1_ij(2,2),   ! normalized x_ij
     *             xt_ij(2,2),             ! x_ij in transition layer
     *             S(N_Roots,N_Roots),     ! S_matrix 4x4
     *             B(N_Roots,N_Roots),     ! 1/S_matrix
     *             B1(N_Roots,N_Roots),    ! 1/S_matrix (previous)
     *             u(N_Roots),             ! dispersion eqn.roots
     *             u1(N_Roots),            ! dispersion eqn.roots (previous)
     *             F(N_Roots),             ! F-matrix
     *             TTT(N_Roots,N_Roots),   ! T1*...Tn*Ss*exp(ikuz)
     *             X(N_Roots,N_Roots),     ! X-matrix of interface
     *             M_rr(N_Half,N_Half),    ! product of operators
     *             M_rt(N_Half,N_Half),    ! product of operators
     *             M_tr(N_Half,N_Half),    ! product of operators
     *             M_tt(N_Half,N_Half),    ! product of operators
     *             W_rr(N_Half,N_Half),    ! product of operators
     *             W_rt(N_Half,N_Half),    ! product of operators
     *             W_tr(N_Half,N_Half),    ! product of operators
     *             W_tt(N_Half,N_Half)     ! product of operators

        Real*8     g0,             ! normalized incidence
     *             k0,             ! k0 = 2*pi/wave
     *             Upper_Limit     ! max matr.element (lim)

        Real*4     xabs,           ! normalization factor (x0 in substrate)
     *             Thickness_Tr,   ! Thickness/2 of transition layer
     *             T_tr,           ! thickness of 1 trans. sublayer
     *             TranLaw(N_Tran),! law of transition
     *             Sigma,          ! r.m.s. roughness
     *             TERang,         ! critical TER angle
     *             c1, c2          ! work cells

        Integer    i, j, l

c Transition layers and roughness
c cannot be used simultaneously:
        if (abs(Sigma).gt.1.E-20) Then !--+
          ifail = 111                    !|
          goto 100                       !|
        endif  !--------------------------+

        T_Tr = 2.*Thickness_Tr / N_Tran

        do    l=1,N_Tran    !================================+
          j_Used = j_Used+1                                 !|
          c1     = TranLaw(l)                               !|
          c2     = 1.-c1                                    !|
          do i=1,2  !====================================+   |
          do j=1,2  !==================================+ |   |
            xt_ij(i,j) = c1*x_ij(i,j) + c2*x1_ij(i,j) !| |   |
          enddo  !=====================================+ |   |
          enddo  !=======================================+   |
          Call  Copy_Layer        (B, B1, u, u1, N_Fields)  !|
          Call  Process_Layer_M   (g0, TERang,              !|
     *                             xt_ij, xabs, k0,         !|
     *                             T_Tr,                    !|
     *                             Magnetic,                !|
     *                             S, B,                    !|
     *                             u, F, N_Fields,          !|
     *                             ifail)                   !|
          if (ifail.ne.0)     goto 100  !--------------------+--+
          Call  Process_Interface (k0, TERang,              !|  v
     *                             Sigma,                   !|this must be 0!
     *                             B1, S, u1, u,            !|
     *                             F, X, TTT,               !|
     *                             N_Fields,                !|
     *                             Upper_Limit,             !|
     *                             Kod, ifail)              !|
          if (ifail.ne.0)     goto 100  !--------------------+--+
          Call Process_Interface2 (X, F,                    !|  v
     *                             M_rr, M_rt,              !|
     *                             M_tr, M_tt,              !|
     *                             W_rr, W_rt,              !|
     *                             W_tr, W_tt,              !|
     *                             N_Fields, j_Used,        !|
     *                             ifail)                   !|
          if (ifail.ne.0)     goto 100  !--------------------+--+
        enddo  !=============================================+  v

  100   return
        end
