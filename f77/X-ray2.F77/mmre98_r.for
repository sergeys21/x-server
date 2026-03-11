c
c
c
c                         +------------------+
c ATTENTION (99/04/26) -- |RECURSION VERSION |
c                         +==================+
c                         based on common unified
c                         subroutines for gid_sl32,
c                                         gid_sl50,
c                                         gid_sl55,
c                                         gid_sl97,
c                                         gids_10
c                                         gids_97
c                                         gids_98
c                                         ter_sl10
c                                         ter_sl12 == ter_sl97
c                                         trds_10
c                                         trds_97
c
c *** Version 98.0    *** -- Same as 99.1.1b (just renamed) - to be
c                            called from mag_sl98 (hard x-rays
c                            approximation).
c *** Version 99.1.1b *** -- Modified function parameters for v.1.1b.
c                            Also the circular polarization is added.
c *** Version 99.1.1 ***  -- Printout for wavefields with simultaneous
c                            switching the roughness transition layers
c                            off.
c *** Version 99.1.0 ***  -- Interface roughness /transition layers
c                            included.
c
c +=====================================================+
c | These routines compute TER reflection coefficients  |
c |   R0 from MAGNETIC multilayer at 1 angular point    |
c |               Incident angle = f0                   |
c |     Input data are transferred to MMRE thru the     |
c |                 common /MAGdat/                     |
c |-----------------------------------------------------|
c |Contents:                                            |
c | 1a. Subroutine            MMRE                      |
c | 2a. Subroutine:           Algorithm_Type            |
c | 3a. Subroutine:           Normalize_x_ij            |
c |-----------------------------------------------------|
c | 1b. Subroutine            Process_Layer_M           |
c | 2b. Subroutine            Smatrix_M                 |
c | 3b. Subroutine            Dispersion_Roots_M        |
c | 4b. Subroutine            Sort_Roots_C              |
c |-----------------------------------------------------|
c |Subroutine calls: see  --  RFMs_sss.for              |
c |                           MMtransr.for              |
c |                           Flds_Uni.for              |
c +=====================================================+

c####################################################### 1
        Subroutine  MMRE98 (f0, f0_, E0, Es, R0s, Phase_shift,
     *                      Circ_pm, Circ_Diff, Circ_Ratio,
     *                      iReqFields, Ds, Dp, u)
c  ** Total_External_Reflex_from_Multi **
c--------------------------------------------------------
c+===============================+
        Include 'mag_sl99.inc'  !|
c+===============================+
c The above INC file contains smth like:
cc      Parameter       (N_Top_Max     = 1001,
cc   *                   N_Total_Max   = N_Top_Max+1)
cc      Parameter       (N_Tran_Max    = 101)           !50

        Integer    N_Roots, N_Half
        Parameter (N_Roots = 4,
     *             N_Half  = N_Roots/2)
        Integer    N_Fields, N_Plus
        Parameter (N_Fields = 4,
     *             N_Plus = N_Fields/2)

        Complex*16 xn_ij(2,2),xn1_ij(2,2),          ! normalized p-z-s matrix
     *             E0(2), Es(2),                    ! incident,reflected waves (sigma & pi)
     *             Ds(N_Roots, N_TT_Max),           ! sigma-Fields in layers
     *             Dp(N_Roots, N_TT_Max),           ! pi-Fields in layers
     *             u (N_Roots, N_TT_Max),           ! dispersion eqn.roots
     *             F (N_Roots, N_TT_Max),           ! F-matrices
     *             M_rr(N_Half,N_Half,N_Total_Max), ! product of operators
     *             M_rt(N_Half,N_Half,N_Total_Max), ! product of operators
     *             M_tr(N_Half,N_Half),             ! product of operators
     *             M_tt(N_Half,N_Half),             ! product of operators
     *             V,                               ! work cell (Dp=V*Ds)
     *             W_rr(N_Half,N_Half),             ! product of operators
     *             W_rt(N_Half,N_Half),             ! product of operators
     *             W_tr(N_Half,N_Half,N_Total_Max), ! product of operators
     *             W_tt(N_Half,N_Half,N_Total_Max), ! product of operators
     *             S(N_Roots,N_Roots),              ! S_matrix 4x4
     *             B(N_Roots,N_Roots),              ! 1/S_matrix
     *             B1(N_Roots,N_Roots),             ! 1/S_matrix (previous)
     *             TTT(N_Roots,N_Roots),            ! T1*...Tn*Ss*exp(ikuz)
     *             X (N_Roots,N_Roots)              ! X-matrix: B1*S

        Real*8     g0,                              ! normalized angle: g0=f0/sqrt(x0)
     *             k0,                              ! k0 = 2*pi/wave
     *             Phs(2),                          ! phase of Es(2)
     *             Upper_Limit                      ! max matr.element (lim)

        Real*4     R0s,                             ! reflection coefficient
     *             f0,                              ! incidence angle in user's units
     *             f0_,                             ! sinus of incidence angle
     *             Phase_shift,                     ! between Es(sigma) & Es(pi)
     *             Circ_pm(2),                      ! reflectivity for +-circ.incident polar.
     *             Circ_Diff,                       ! difference in +-circ.reflectivity
     *             Circ_Summ,                       ! Summ of +-circ.reflectivity
     *             Circ_Ratio,                      ! Ratio (I+ - I-)/(I+ + I-)
     *             pi,                              ! 3.1416
     *             Sgm_i                            ! rms roughness (i-th layer)

        Integer    j_Used,                          ! number of used layers incl.transition layers
     *             iReqFields,                      ! the flag to request fields
     *             Kod,                             ! flag to compute TTT
     *             Magnetic_,                       ! i-th layer magnetic type
     *             i, j, k                          ! work sell

        Complex*16 x_ij(2,2,N_TT_Max)               ! polarizabilities matrix
        Complex*8  x0(N_TT_Max)                     ! Top layer x0
        Real*4     Wave,                            ! X-wavelength
     *             xabs,                            ! abs(x0_max)
     *             TERang,                          ! critical TER angle (radians)
     *             Thickness(N_Total_Max),
     *             Sigma(N_Total_Max),
     *             Thickness_Tr(N_Total_Max),
     *             TranLaw(N_Tran_Max)
        Integer    Magnetic(N_TT_Max),
     *             N_Top,                           ! top layer sublayers
     *             N_Total,                         ! total of layers
     *             N_TT,                            ! total of layers+vacuum
     *             N_Used,                          ! number of used layers
     *             N_Tran,                          ! number of trans.sublr.
     *             MaxOrder,                        ! max order of matrix elem.
     *             iDebug,                          ! Debug flag
     *             ifail                            ! failure code
        Common  /MAGdat/  x_ij,                     ! (2,2,N_TT_Max)
     *                    x0,                       ! (N_TT_Max)
     *                    Wave,
     *                    xabs,                     ! abs(x0_max)
     *                    TERang,                   ! critical TER angle (radians)
     *                    Thickness,                ! (N_Total_Max),
     *                    Sigma,                    ! (N_Total_Max),
     *                    Thickness_Tr,             ! (N_Total_Max),
     *                    TranLaw,                  ! (N_Tran_Max),
     *                    Magnetic,                 ! (N_TT_Max),
     *                    N_Top, N_Total,
     *                    N_TT, N_Used, N_Tran,
     *                    MaxOrder, iDebug, ifail
c--------------------------------------------------------
        R0s    = 0.
        Es(1)  = (0., 0.)
        Es(2)  = (0., 0.)
        ifail  = 0

        do      i=1,N_Fields  !=======+
        do      k=1,N_TT      !=====+ |
          Ds(i,k) = (0.D0, 0.D0)   !| |D_sigma
          Dp(i,k) = (0.D0, 0.D0)   !| |D_pi
          u (i,k) = (0.D0, 0.D0)   !| |
          F (i,k) = (1.D0, 0.D0)   !| |
        enddo  !====================+ |
        enddo  !======================+

c       xabs   = Abs(x0(N_TT))                  !this is calculated in main
c       TERang = Sqrt(xabs)                     !this is calculated in main
        if (abs(xabs).lt.1.E-32) Then !--+
          ifail = 555                   !|
          return                        !|
        endif  !-------------------------+

        pi     = 4.* Atan(1.)
        k0     = 2.* pi/wave
        g0     = f0_ / TERang           !normalized sinus of incidence angle
        Sgm_i  = 0.
c       N_Plus = N_Fields / 2

c+-------------------------+
c ccc   MaxOrder = 7      !| DEBUG !!!
c+-------------------------+
        Upper_Limit = 10.**MaxOrder

c At f0=0 vacuum matrix becomes singular and cannot be inverted:
c                      |1 1|
c                      |0 0|
c Then, reflection coefficient cannot be determined and we set it to 1:

        if (Abs(g0).lt.1.0D-10) g0 = 1.0D-10

        if (g0.lt.0.)      Then  !-----+
          R0s = 1.                    !|
          Return                      !|
        endif  !-----------------------+

        Kod = 1                         !Don't compute TTT
        Call    M1 (N_Roots,TTT,N_Fields)

c
c                     Process Vacuum
c
c+-------------------------------------------------------+
        N_used = 0                                      !|
        j_Used = 0                                      !|
        Call    Normalize_x_ij (x_ij(1,1,N_Used+1),     !|
     *                          xn_ij, xabs, 2, 2)      !|
        Call    Process_Layer_M (g0, TERang,            !|
     *                           xn_ij, xabs, k0,       !|
     *                           0.,                    !|Thickness(vacuum)
     *                           Magnetic(N_Used+1),    !|
     *                           S, B,                  !|
     *                           u(1,N_Used+1),         !|
     *                           F(1,N_Used+1),         !|
     *                           N_Fields,ifail)        !|
        if (ifail.ne.0) goto 100 !-----------------------+-+
c+-------------------------------------------------------+ v

c
c                     Process Top Layer
c
        do    N_Used=1,N_Total !============================+
c+---------------------------------------------+            |
c| Process sublayer-sublayer transition in TOP |            |
c+=============================================+            |
          if (Thickness_Tr(N_Used).gt.0. .AND.             !|
     *        iReqFields.eq.0)  Then  !-------------------+ |
            Call  Normalize_x_ij (x_ij(1,1,N_Used),      !| |
     *                            xn_ij, xabs, 2, 2)     !| |
            Call  Normalize_x_ij (x_ij(1,1,N_Used+1),    !| |
     *                            xn1_ij, xabs, 2, 2)    !| |
c This might be a potential problem if e.g. magnetic_1=1  | |
c and magnetic_2=2. Hopefully, it is not physical!        | |
            Magnetic_ = Max(Magnetic(N_Used),            !| |
     *                      Magnetic(N_Used+1))          !| |
            Call   Transi_rM98(Thickness_Tr(N_Used),     !| |
     *                         Sigma(N_Used),            !| |this must be 0!
     *                         TranLaw, N_Tran,          !| |
     *                         xn_ij, xn1_ij, xabs,      !| |
     *                         TERang, g0, k0,           !| |
     *                         Magnetic_,                !| |
     *                         N_Roots, N_Fields,        !| |
     *                         N_Half,  j_Used,          !| |
     *                         S, B, B1,                 !| |
     *                         u(1,N_Used+1),            !| |
     *                         u(1,N_Used),              !| |
     *                         F(1,N_Used+1), X, TTT,    !| |
     *                         M_rr(1,1,N_Used),         !| |
     *                         M_rt(1,1,N_Used),         !| |
     *                         M_tr, M_tt,               !| |
     *                         W_rr, W_rt,               !| |
     *                         W_tr(1,1,N_Used),         !| |
     *                         W_tt(1,1,N_Used),         !| |
     *                         Upper_Limit, Kod,         !| |
     *                         ifail )                   !| |
              if (ifail.ne.0)     goto 100  !-------------+-+-+
          endif  !----------------------------------------+ | v
                                                           !|
          Call  MtoM (N_Roots,N_Roots,B,B1,N_Fields)       !|
                                                           !|
          Call  Normalize_x_ij (x_ij(1,1,N_Used+1),        !|
     *                          xn_ij, xabs, 2, 2)         !|
          j_Used = j_Used+1                                !|
          Call  Process_Layer_M (g0, TERang,               !|
     *                           xn_ij, xabs, k0,          !|
     *                           Thickness(N_Used),        !|
     *                           Magnetic(N_Used+1),       !|
     *                           S, B,                     !|
     *                           u(1,N_Used+1),            !|
     *                           F(1,N_Used+1),            !|
     *                           N_Fields, ifail)          !|
          if (ifail.ne.0)     goto 100  !-------------------+-+
          if (iReqFields.eq.0) Then !----+                 !| V
            Sgm_i = Sigma(N_Used)       !|                  |
          else  !------------------------|                  |
            Sgm_i = 0.                  !|                  |
          endif  !-----------------------+                  |
          Call  Process_Interface (k0, TERang, Sgm_i,      !|
     *                             B1, S,                  !|
     *                             u(1,N_Used),            !|
     *                             u(1,N_Used+1),          !|
     *                             F(1,N_Used+1),          !|
     *                             X, TTT, N_Fields,       !|
     *                             Upper_Limit,            !|
     *                             Kod, ifail)             !|
          if (ifail.ne.0)       goto 100  !-----------------+-+
          Call Process_Interface2(X, F(1,N_Used+1),        !| v
     *                            M_rr(1,1,N_Used),        !|
     *                            M_rt(1,1,N_Used),        !|
     *                            M_tr, M_tt,              !|
     *                            W_rr, W_rt,              !|
     *                            W_tr(1,1,N_Used),        !|
     *                            W_tt(1,1,N_Used),        !|
     *                            N_Fields, j_Used,        !|j_Used is here to detect whether
     *                            ifail)                   !|this is the first layer or not?
          if (ifail.ne.0)     goto 100  !-------------------+-+
          if (N_Used.lt.N_Total) Then  !---------+          | v
c W_tr(N_Used+1)=W_tr(N_Used)                   !|          |
c W_tt(N_Used+1)=W_tt(N_Used)                   !|          |
              Call Backup2(W_tr(1,1,N_Used),    !|          |
     *                     W_tt(1,1,N_Used),    !|          |
     *                     W_tr(1,1,N_Used+1),  !|          |
     *                     W_tt(1,1,N_Used+1),  !|          |
     *                     N_Half,N_Plus)       !|          |
          endif  !-------------------------------+          |
        enddo  !============================================+

c
c              Find the solutions Es(1),Es(2)
c
c+---------------------------------------------------------------------+
        Es(1) = W_rt(1,1)*E0(1)+W_rt(1,2)*E0(2)                       !|
        Es(2) = W_rt(2,1)*E0(1)+W_rt(2,2)*E0(2)                       !|
                                                                      !|
c Compute the reflection coefficient:                                  |
        R0s  = sngl((CDabs(Es(1)))**2+(CDabs(Es(2)))**2)              !|
                                                                      !|
c Compute the phase shift between Es(1) and Es(2):                    !|
        do i=1,2  !========================================+           |
          if (Cdabs(Es(i)).gt.0.) Then !-----------------+ |           |
            Phs(i) = Datan2 (Dimag(Es(i)),Dreal(Es(i))) !| |           |
          else  !----------------------------------------| |           |
            Phs(i) = 0.0D0                              !| |           |
          endif  !---------------------------------------+ |           |
        enddo  !===========================================+           |
        Phase_shift = sngl((Phs(2)-Phs(1))/pi)                        !|
                                                                      !|
c The difference in reflectivity for E+ and E- circular                |
c polarizations, I^{+} - I^{-} :                                       |
        Circ_Diff = sngl(2.*Dimag(W_rt(1,1)*Dconjg(W_rt(1,2))         !|
     +                          + W_rt(2,1)*Dconjg(W_rt(2,2))))       !|
c This inserted for compatibility with generic mag_sl v.4.0            |
c (perhaps different order of roots sorting) -- 2003/07/02             |
        Circ_Diff =-Circ_Diff                                         !|
c I^{+} + I^{-}                                                        |
        Circ_Summ = sngl(CdAbs(W_rt(1,1))**2 + CdAbs(W_rt(1,2))**2    !|
     +                 + CdAbs(W_rt(2,1))**2 + CdAbs(W_rt(2,2))**2)   !|
c I^{+} and I^{-}                                                      |
        Circ_pm(1) = 0.5*(Circ_Summ+Circ_Diff)                        !|
        Circ_pm(2) = 0.5*(Circ_Summ-Circ_Diff)                        !|
c Ratio [I^{+} + I^{-}] / [I^{+} + I^{-}]                              |
        if (Circ_Summ.gt.0.) Then  !-------------+                     |
          Circ_Ratio = Circ_Diff / Circ_Summ    !|                     |
        else  !----------------------------------|                     |
          Circ_Ratio = 0.                       !|                     |
        endif !----------------------------------+                     |
                                                                      !|
c+---------------------------------------------------------------------+

        if (iReqFields.eq.0) Return

c
c                  Compute the wavefields
c

c Vacuum fields:
        Ds(1,1) = E0(1)
        Ds(2,1) = E0(2)
        Ds(3,1) = Es(1)
        Ds(4,1) = Es(2)
c Substrate fields (T_n = W_tt^n * T_0):
        Ds(1,N_Used+1) = W_tt(1,1,N_Used)*Ds(1,1)/F(1,N_Used+1)  !(F=1)
     +                 + W_tt(1,2,N_Used)*Ds(2,1)/F(2,N_Used+1)  !(F=1)
        Ds(2,N_Used+1) = W_tt(2,1,N_Used)*Ds(1,1)/F(1,N_Used+1)  !(F=1)
     +                 + W_tt(2,2,N_Used)*Ds(2,1)/F(2,N_Used+1)  !(F=1)
        Ds(3,N_Used+1) = (0.0D0, 0.0D0)
        Ds(4,N_Used+1) = (0.0D0, 0.0D0)
c The rest fields:
c R_k = (1-M_rt^{k+1}*W_tr^k)^{-1} * (M_rr^{k+1}*R_{k+1) + M_rt^{k+1}*W_tt^k*T_0)
c T_k = W_tt^k*T_0 + W_tr^k*R_k
        do      k = N_Used-1,1,-1   !===========+
          Call  Find_Fields (Ds(1,1),          !|T_0
     *                       Ds(N_Plus+1,k+2), !|R_{k+1}
     *                       Ds(1,k+1),        !|T_k
     *                       Ds(N_Plus+1,k+1), !|R_k
     *                       W_tt(1,1,k),      !|
     *                       W_tr(1,1,k),      !|
     *                       M_rt(1,1,k+1),    !|
     *                       M_rr(1,1,k+1),    !|
     *                       N_Plus,ifail)     !|
          if (ifail.ne.0)       goto 100  !-----+-+
        enddo  !================================+ V

c The program Norm_Fields calculates [D_k]*[F_k^(U)], so that if
c we look at the paper PRB 54 (1996) 8150-8162, Eq.(24) and below,
c than a part of Q_{nij} is already included into E_{nij} -- namely
c exp(iku^{IN}_{ni}*Z_n) and exp(iku^{OUT}_{nj}*Z_n.
        do      k = 2, N_Used   !=============+
          Call  Norm_Fields (Ds(1,k),        !|T_k
     *                       Ds(N_Plus+1,k), !|R_k
     *                       F(1,k),         !|Ft_k
     *                       F(N_Plus+1,k),  !|Fr_k
     *                       N_Plus)         !|
        enddo  !==============================+

c       iDebug=0
c       iDebug=1
        if (iDebug.ne.0)        Then !---------------------+
          write (3,50)  g0                                !|
  50      format('The Ds-wavefields for g0/TERang=',g9.3) !|
                                                          !|
          do    k = 1, N_Used  !====================+      |
            if (k.lt.N_Used)  Then !---+            |      |
              j = N_Fields            !|            |      |
            else  !--------------------|            |      |
              j = N_Plus              !|            |      |
            endif  !-------------------+            |      |
            write (3,55) k,(Abs(Ds(i,k+1)),i=1,j)  !|      | !!!!!!! DEBUG
  55        format(i5,4g12.5)                      !|      | !!!!!!! DEBUG
          enddo  !==================================+      |
        endif  !-------------------------------------------+

c Calculate Dp:
        do      k=1,N_Used+1   !=======================+
c This was a bug in ver.1b-4 partly fixed 2003/07/02:  |
c The layer may be magnetic but if the polarization    |
c is along Y or Z, we do not see it!                   | TO BE FIXED!!!!!!
          if   ((Magnetic(k) .ne. 0)        .AND.     !|
     *     (abs(x_ij(1,2,k)) .gt. (1.E-20)) .AND.     !|
     *     (abs(x_ij(2,1,k)) .gt. (1.E-20))) Then !--+ |
                                                    !| |
            Call Normalize_x_ij (x_ij(1,1,k),       !| |
     *                           xn_ij, xabs, 2, 2) !| |
                                                    !| |
            do      i=1,N_Fields !==============+    | |
              V      = u(i,k)                  !|    | |
              V      = (V**2-g0**2-xn_ij(1,1)) !|    | |
     /               /  xn_ij(1,2)             !|    | | TO CLARIFY !!!!!!
              Dp(i,k)= V * Ds(i,k)             !|    | |
            enddo  !============================+    | |
            if (iDebug.ne.0)  Then !-----------+     | |
              write (3,51)                    !|     | |
  51          format('The Dp-wavefields:')    !|     | |
              write (3,55) k-1,(Abs(Dp(i,k)), !|     | | !!!!!!!!!!! DEBUG
     *                     i=1,N_Fields)      !|     | | !!!!!!!!!!! DEBUG
            endif  !---------------------------+     | |
          else  !------------------------------------| |
            Dp(1,k) =(0.,0.)                        !| |
            Dp(2,k) = Ds(2,k)                       !| |
            Dp(3,k) =(0.,0.)                        !| |
            Dp(4,k) = Ds(4,k)                       !| |
            Ds(2,k) =(0.,0.)                        !| |
            Ds(4,k) =(0.,0.)                        !| |
          endif  !-----------------------------------+ |
c Sort the wavefields in decreasing order and neglect  |
c the small ones (get Sort_H_Amp from GISU2_98.for)    |
c         Call Sort_H_Amp (D0(1,k), Dh(1,k), u(1,k),  !|
c    *                    N_Fields, N_Fact(k), i_Born)!|
        enddo  !=======================================+

        return

  100   continue
c ifail=100  Dispersion_Roots: Count of Im(u) # 2
c ifail=120  Smatrix_M: x_ij(1,2)=0
c ifail=200  Inverse_Matrix_simple  size < 1
c ifail=201  Inverse_Matrix_simple  Det=0. (for 1x1 matrix)
c ifail=202  Inverse_Matrix_simple  Det=0. (for 2x2 matrix)
        write (3,*)     'MMRE: failure code=',ifail,
     *                  '-- at f0: ',f0
        Return
        end

c =======================================================

        Subroutine  Normalize_x_ij (x,xn,xdel,ndim,nused)
        Integer         ndim, nused, i, j
        Complex*16      x(ndim,ndim), xn(ndim,ndim)
        Real*4          xdel

        do i=1,nused  !==================+
        do j=1,nused  !================+ |
          xn(i,j) = x(i,j) / xdel     !| |
        enddo  !=======================+ |
        enddo  !=========================+
        return
        end

c ####################################################### 1a

        Subroutine  Process_Layer_M (g0, TERang,
     *                               x_ij, xabs, k0,
     *                               Depth, Magnetic,
     *                               S, B, u, F,
     *                               N_Fields, ifail)
c +-------------------------------+
c |For MAGNETIC x-ray reflection. |
c +===============================+
        Integer         N_Msize,
     *                  N_Fields
        Parameter      (N_Msize = 4)
        Complex*16      x_ij(2,2),                 ! normalized polariz.
     *                  S(N_Msize,N_Fields),       ! layer S-matrix
     *                  B(N_Msize,N_Fields),       ! layer (1/S)-matrix
     *                  F(N_Fields),               ! layer F-matrix
     *                  u(N_Fields),               ! layer roots
     *                  q, qw,                     ! work cells
     *                  im/(0.0D0,1.0D0)/


        Real*8          g0,                        ! normalized incid.angle g0=f0/sqrt(x0)
     *                  k0,                        ! k0 = 2*pi/wave
     *                  qr,qi,                     ! Re(qw),Im(qw)
     *                  qmax /100./                ! The upper limit for
                                                   ! exponent value. The
                                                   ! actual limit for PC
                                                   ! is about 707; for
                                                   ! other platforms it
                                                   ! may differ!
        Real*4          Depth,
     *                  TERang,
     *                  xabs                       ! abs(x0_max)

        Integer         Magnetic, N_Plus,
     *                  ifail, i

        ifail = 0

c+----------------------------------------------------------+
                                                           !|
c Compute u-roots and S-matrix of the layer:                |
        Call    Smatrix_M (g0,                             !|
     *                     x_ij, xabs,                     !|
     *                     Magnetic,                       !|
     *                     N_Fields, S, u, ifail)          !|
        if (ifail.ne.0) Return                             !|
c+----------------------------------------------------------+

c+---------------------------------------------------+
c|Find inverse matrix B=1/S of matrix S:             |
        if (N_Fields.le.2)      Then  !-----------+  |
c For N_Msize*N_Msize arrays containing           |  |
c 2x2 or 1x1 matrices                             |  |
c (just faster):                                  |  |
          Call  Inverse_Matrix_simple (S, B,     !|  |
     *                                 N_Msize,  !|  |
     *                                 N_Fields, !|  |
     *                                 ifail)    !|  |
        else  !-----------------------------------+  |
c For 4x4 arrays containng N_Fields*N_Fields      |  |
c matrices:                                       |  |
          Call  Inverse_Matrix        (S, B,     !|  |
     *                                 N_Fields, !|  |
     *                                 ifail)    !|  |
        endif  !----------------------------------+  |
        if (ifail.ne.0) Return                      !|
c+---------------------------------------------------+

c+------------------------------------------------------+
c|Compute F(t) matrix:                                  |
        if (abs(Depth).gt.1.E-32) Then !-------------+  |
                                                    !|  |
          N_Plus = N_Fields / 2                     !|  |
          q      = -im*k0*Depth*TERang              !|  |
c The growing exponents are inverted                 |  |
          do      i=1,N_Plus    !=================+  |  |
            qw =-q*u(i)                          !|  |  |
            qr = dReal(qw)                       !|  |  |
            qi = dImag(qw)                       !|  |  |
            if (qr.gt.+qmax) qw=dCmplx(+qmax,qi) !|  |  |
            if (qr.lt.-qmax) qw=dCmplx(-qmax,qi) !|  |  |
            F(i) = exp(qw)                       !|  |  |
          enddo   !===============================+  |  |
c The decreacing exponents aren't inverted           |  |
          do      i=N_Plus+1,N_Fields !===========+  |  |
            qw =+q*u(i)                          !|  |  |
            qr = dReal(qw)                       !|  |  |
            qi = dImag(qw)                       !|  |  |
            if (qr.gt.+qmax) qw=dCmplx(+qmax,qi) !|  |  |
            if (qr.lt.-qmax) qw=dCmplx(-qmax,qi) !|  |  |
            F(i) = exp(qw)                       !|  |  |
          enddo   !===============================+  |  |
                                                    !|  |
        else  !--------------------------------------|  |
                                                    !|  |
          do      i=1,N_Fields  !===+                |  |
            F(i) = (1.0D0, 0.0D0)  !|                |  |
          enddo  !==================+                |  |
                                                    !|  |
        endif  !-------------------------------------+  |
c+------------------------------------------------------+

        Return
        End

c ####################################################### 2a

        Subroutine  Smatrix_M (g0,
     *                         x_ij, xabs,
     *                         Magnetic,
     *                         N_Fields, S, u, ifail)
c +-------------------------------+
c |For MAGNETIC x-ray reflection. |
c +===============================+
        Integer         N_Msize,
     *                  N_Fields
        Parameter      (N_Msize = 4)
        Complex*16      x_ij(2,2),
     *                  eps,                       !sqrt(1+x0)
     *                  S(N_Msize,N_Fields),       ! S-matrix
     *                  u(N_Fields),               ! dispersion roots
     *                  V,
     *                  W,
     *                  ss

        Real*8          g0                         ! normalized incid.angle g0=f0/sqrt(x0)

        Real*4          xabs

        Integer         Magnetic,
     *                  ifail,
     *                  j, iA, iB

c Find dispersion equation roots:
        Call    Dispersion_Roots_M (g0, u, x_ij,
     *                              Magnetic, N_Fields, ifail)
        if (ifail.ne.0) return

c This was a bug in ver.1b corrected 2003/07/02:
c The layer may be magnetic but if the polarization
c is along Y or Z, we do not see it!
        if     (Magnetic .eq. 1) Then !-----------------------+
                                                             !|
c Magnetic layer:                                             |
c---------------                                              |
          if (abs(x_ij(1,2)) .lt. (1.E-20)) Then !-+          |
            write (3,*) ' Smatrix_M: ',           !|          |
     *                  'x_ij(1,2)=0'             !|          |
            ifail = 120                           !|          |
            return                                !|          |
          endif  !---------------------------------+          |
          ss = g0*g0 + x_ij(1,1)                             !|
c cccc    ss = g0*g0 + x_ij(2,2)                             !|
          do    j=1,N_Fields  !=============+                 |
            V = (u(j)*u(j)-ss)/x_ij(1,2)   !|                 |
c cccc      V = x_ij(2,1)/(u(j)*u(j)-ss)   !|                 |
            W = V * u(j)                   !|                 |
                                           !|                 |     Magnetic
            S(1,j) = 1.                    !|                 |  | 1, 1, 1, 1 |
            S(2,j) = V                     !|                 |  | V, V, V, V |
            S(3,j) = u(j)                  !|                 |  | u, u, u, u |
            S(4,j) = W                     !|                 |  | w, w, w, w |
                                           !|                 |
c           write (3,1) "S-matrix:",       !|                 |
c    *                  (S(i,j),i=1,4)     !|                 |
c 1         format(1x,a,3x,10(' (',g11.5,  !|                 |
c    *                 ',',g11.5,') ',:))  !|                 |
          enddo  !==========================+                 |
                                                             !|
        elseif ((Magnetic .eq. 2) .OR.                       !|
     *          (Magnetic .eq. 3)) Then  !--------------------|
                                                             !|
c Magnetic layer: case M || Y  or M || Z                      |
c---------------------------------------                      |
          do      j=1,N_Fields  !===========+                 |
            iA     = j-2*(j/2)             !|   1, 0, 1, 0    |
            iB     = iA-1                  !|   0,-1, 0,-1    |    Magnetic_Y
                                           !|                 |    Magnetic_Z
            S(1,j) = iA                    !|                 |  | 1, 0, 1, 0 |
            S(2,j) = iB                    !|                 |  | 0,-1, 0,-1 |
            S(3,j) = iA*u(j)               !|                 |  | u, 0, u, 0 |
            S(4,j) = iB*u(j)               !|                 |  | 0,-u, 0,-u |
c           write (3,1) "S-matrix:",       !|                 |
c    *                (S(i,j),i=1,4)       !|                 |
          enddo  !==========================+                 |
                                                             !|
        else !------------------------------------------------|
                                                             !|
c Non-Magnetic layer:                                         |
c-------------------                                          |
          eps = Sqrt(1.+x_ij(1,1)*xabs)                      !|
          do      j=1,N_Fields  !===========+                 |
            iA     = j-2*(j/2)             !|   1, 0, 1, 0    |
            iB     = 1-iA                  !|   0, 1, 0, 1    |
                                           !|                 |   NON_Magnetic
            S(1,j) = iA                    !|                 |  | 1, 0, 1, 0 |
            S(2,j) = iB*eps                !|                 |  | 0, e, 0, e |
            S(3,j) = iA*u(j)               !|                 |  | u, 0, u, 0 |
            S(4,j) = iB*u(j)/eps           !|                 |  | 0,u/e,0,u/e|
c           write (3,1) "S-matrix:",       !|                 |
c    *                (S(i,j),i=1,4)       !|                 |
          enddo  !==========================+                 |
                                                             !|
        endif  !----------------------------------------------+
        return
        end

c ####################################################### 4a

        Subroutine  Dispersion_Roots_M (g0, u, x_ij,
     *                                  Magnetic, N_Fields, ifail)
c +-------------------------------+
c |For MAGNETIC x-ray reflection. |
c +===============================+
c------------------------------------------------------
c      Forms and solves the polynomial dispersion
c    equation of MAGNETIC x-ray reflection.
c
c  On exit the roots are sorted in descending of Im(u).
c------------------------------------------------------
        Integer     N_Fields

        Integer    Magnetic, ifail,
     *             i, nPosImag

        Complex*16 u(N_Fields),
     *             x_ij(2,2),
     *             s0, s1, s2, s3,
     *             ssg, spi
        Real*8     g0                              ! normalized incid.angle g0=f0/sqrt(x0)

        s1 = g0*g0 + x_ij(1,1)
        s2 = g0*g0 + x_ij(2,2)

        if (Magnetic .eq. 1) Then !----------------------------------+
                                                                    !|
c Magnetic layer:                                                    |
c---------------                                                     |
          s0 = (s1+s2)/2.                                           !|
          s3 = (x_ij(1,1)-x_ij(2,2))/2.                             !|
          s3 = s3*s3 + x_ij(1,2)*x_ij(2,1)                          !|
          s3 = Cdsqrt(s3)                                           !|
                                                                    !|
          u(1) = Cdsqrt(s0+s3)                                      !|
          u(2) = Cdsqrt(s0-s3)                                      !|
          u(3) =-u(1)                                               !|
          u(4) =-u(2)                                               !|
c+--------------------------+                                        |
c|Sort the roots in the     |                                        |
c|descending order of Im(u):|                                        |
c+--------------------------+                                        |
          Call  Sort_Roots_C (u,N_Fields)                           !|
                                                                    !|
c+--------------------+                                              |
c|Check for the right |                                              |
c|Number of Im's:     |                                              |
c+--------------------+                                              |
          nPosImag = 0                                              !|
          do      i=1,N_Fields  !===============+                    |
            if (Dimag(u(i)).gt.0.)             !|                    |
     *                   nPosImag = nPosImag+1 !|                    |
          enddo  !==============================+                    |
          if (nPosImag.ne.2)  Then  !-----------+                    |
            write (3,*) ' Dispersion_Roots: ', !|                    |
     *                  'Count of Im(u) # 2'   !|                    |
            ifail = 100                        !|                    |
            return                             !|                    |
          endif  !------------------------------+                    |
                                                                    !|
        elseif ((Magnetic .eq. 2) .OR.                              !|
     *          (Magnetic .eq. 3))  Then  !--------------------------|
                                                                    !|
c Special magnetic case: M||Y or M||Z; de-coupled sigma & pi         |
c-----------------------------------------------------------         |
          ssg   = g0*g0 + x_ij(1,1)                                 !|
          spi   = g0*g0 + x_ij(2,2)                                 !|
c         write (3,'(1001g14.7)') (g0*g0+x_ij(1,1)), ssg, spi       !|
          u(1)  = Sqrt(ssg)                                         !|
          u(2)  = Sqrt(spi)                                         !|
          if (Dimag(u(1)).lt.0.) u(1)=-u(1)                         !|
          if (Dimag(u(2)).lt.0.) u(2)=-u(2)                         !|
          u(3)  = -u(1)                                             !|
          u(4)  = -u(2)                                             !|
                                                                    !|
c+=========================================+                         |
c|ATTENTION !!! These roots are not sorted.|                         |
c|              That might cause some dis- |                         |
c|              continuities if one plots  |                         |
c|              these roots or the x-ray   |                         |
c|              wavefields as a function   |                         |
c|              of layer number            |                         |
c+=========================================+                         |
                                                                    !|
        else  !------------------------------------------------------|
                                                                    !|
c Non-Magnetic layer:                                                |
c-------------------                                                 |
          u(1) = Cdsqrt(s1)                                         !|
          u(2) = Cdsqrt(s2)                                         !|
          u(2) = u(1)                                               !|
          u(3) =-u(1)                                               !|
          u(4) =-u(2)                                               !|
                                                                    !|
c+=========================================+                         |
c|ATTENTION !!! These roots are not sorted.|                         |
c|              That might cause some dis- |                         |
c|              continuities if one plots  |                         |
c|              these roots or the x-ray   |                         |
c|              wavefields as a function   |                         |
c|              of layer number            |                         |
c+=========================================+                         |
                                                                    !|
        endif  !-----------------------------------------------------+

        return
        end

c ####################################################### 3c

        Subroutine  Sort_Roots_C (u,n)
        Integer         n, i, j
        Complex*16      u(n), v
c--------------------------------------------------------
c  Sort roots in descending of Im(u).
c--------------------------------------------------------
        do      j = 2,n  !==============================+
          v = u(j)                                     !|
                                                       !|
          do    i = j-1,1,-1  !====================+    |
            if (Dimag(u(i)).ge.Dimag(v))  goto 3 !-+->+ |
            u(i+1) = u(i)                         !|  | |
          enddo  !=================================+  v |
                                                     !| |
          i = 0                                      !| |
  3       continue      !<----------------------------+ |
          u(i+1) = v                                   !|
        enddo  !========================================+

        return
        end

