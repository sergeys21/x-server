c
c
c
c                         +------------------+
c ATTENTION (12-10-95) -- |RECURSION VERSION |
c                         +==================+
c                         based on common unified
c                         subroutines for gid_sl32,
c                                         gid_sl50,
c                                         gid_sl55,
c                                         gids_10
c                                         TER_sl10
c                                         ter_sl12 == ter_sl97
c                                         trds_10
c
c *** Version 1.2 ***  -- Interface roughness included
c                         Exponent index = 2
c
c +=====================================================+
c | These routines compute TER reflection coefficients  |
c | R0 from multilayered structure at 1 angular point   |
c |               Incident angle = f0                   |
c |       Input data are transferred to SFM in          |
c |                 common /TERdat/                     |
c |-----------------------------------------------------|
c |Contents:                                            |
c | 1. Real*8 function:       SFMm5                     |
c | 2. Subroutine:            Algorithm_Type            |
c |-----------------------------------------------------|
c |Subroutine calls: see  --  RFMs_URD.for              |
c |                           RFMs_sss.for              |
c |                           TRANSI_r.for              |
c +=====================================================+

c ####################################################### 1
        Real*8  Function        SFMm5 (f0,ipol)
c  ** Total_External_Reflex_from_Multi **
c--------------------------------------------------------
c+===============================+
        Include 'ter_slm5.inc'  !|
c+===============================+
c The above INC file contains smth like:
cc      Parameter       (N_Top_Max     = 1001,
cc   *                   N_Total_Max   = N_Top_Max+1)
cc      Parameter       (N_Tran_Max    = 101)           !50

        Real*4          Daa
        Integer         N_Roots,
     *                  N_Half,
     *                  N_Fields,
     *                  N_Plus
        Parameter       (N_Roots = 4,
     *                   N_Half  = N_Roots/2)
        Parameter       (N_Fields = 2,
     *                   N_Plus = N_Fields/2)
        Parameter       (Daa      = 0.)             ! da/a (not used)

        Complex*16 x0n,                             ! normalized x0
     *             xhn(2)/(0.0D0,0.0D0),
     *                    (0.0D0,0.0D0)/,           ! normalized xh (not used)
     *             gh/(0.0D0,0.0D0)/,               ! normalized exit angle (not used)
     *             Sqrt_Eps,                        ! sigma=1 pi=Sqrt(1+x0)
     *             E(N_Roots,N_Total_Max+1),        ! Fields in layers
     *             u(N_Roots,N_Total_Max+1),        ! dispersion eqn.roots
     *             F(N_Roots,N_Total_Max+1),        ! F-matrix
     *             X(N_Roots,N_Roots),              ! X=S_{k}/S_{k-1} (not used)
     *             M_rr(N_Half,N_Half,N_Total_Max), ! product of operators
     *             M_rt(N_Half,N_Half,N_Total_Max), ! product of operators
     *             M_tr(N_Half,N_Half),             ! product of operators
     *             M_tt(N_Half,N_Half),             ! product of operators
     *             W_rr(N_Half,N_Half),             ! product of operators
     *             W_rt(N_Half,N_Half),             ! product of operators
     *             W_tr(N_Half,N_Half,N_Total_Max), ! product of operators
     *             W_tt(N_Half,N_Half,N_Total_Max), ! product of operators
     *             S(N_Roots,N_Roots),              ! S_matrix 2x2
     *             B(N_Roots,N_Roots),              ! 1/S_matrix
     *             B1(N_Roots,N_Roots),             ! 1/S_matrix (previous)
     *             TTT(N_Roots,N_Roots),            ! T1*...Tn*Ss*exp(ikuz)
     *             Es,                              ! vacuum wave amplitude
     *             ikuz,                            ! i*k*u*z
     *             SUMM,                            ! summ of waves
     *             Im /(0.,1.)/                     ! imaginary 1.


        Complex*8  x0(N_Total_Max+1),               ! Top layer x0
     *             xh/(0.0E0,0.0E0)/                ! xh (not used)

        Real*8     g0,                              ! normalized angle
     *             gi/0.0D0/,                       ! normalized miscut angle (not used)
     *             k0,                              ! k0 = 2*pi/wave
     *             Upper_Limit,                     ! max matr.element (lim)
     *             SW(N_Stand_Max),                 ! standing waves
     *             Positions(N_Total_Max),          ! positions of iterfaces
     *             rkz, ikz,                        ! Real(i*k*u*z), Imag(i*k*u*z)
     *             z,                               ! depth inside layers
     *             zi                               ! same, relative to a closest interface

        Real*4     Pol,                             ! polarization index=ipol (sigma=1,pi=2)
     *             Wave,                            ! X-wavelength
     *             f0,                              ! incidence angle
     *             UnitCoef,                        ! radians/(angular unit)
     *             pi,                              ! 3.1416
     *             xabs,                            ! abs(x0_max)
     *             TERang,                          ! critical TER angle (radian)
     *             Thickness(N_Total_Max),
     *             Thickness_Tr(N_Total_Max),
     *             TranLaw(N_Tran_Max),
     *             Sigma(N_Total_Max),
     *             standing_range(2),
     *             standing_step

        Integer    N_Total,                         ! total of layers
     *             N_Used,                          ! number of used layers
     *             j_Used,                          ! index of used layer
     *             N_Top,                           ! top layer sublayers
     *             N_Tran,                          ! number of trans.sublr.
     *             ifail,                           ! failure code
     *             MaxOrder,                        ! max order of mat.elem.
     *             iDebug,                          ! Debug flag
     *             i_standing,                      ! calculate or not standing waves
     *             i_reference,                     ! reference interface for stand.waves
     *             n_standing,                      ! number of standing wave pts.
     *             ipol,                            ! polarization state (1=sg,2=pi)
     *             Kod,                             ! flag (recursive/matrix) for process_interface
     *             N, k, kk, i, j, j2               ! work cells

c--------
        Common /TER_m5/ Wave, UnitCoef, x0,
     *                  Positions,
     *                  Thickness,
     *                  Thickness_Tr,
     *                  TranLaw,
     *                  Sigma,
     *                  xabs, TERang,
     *                  standing_range, standing_step,
     *                  SW,
     *                  N_Top, N_Total, N_Used, N_Tran,
     *                  MaxOrder, iDebug, ifail,
     *                  i_standing, i_reference, n_standing

c--------------------------------------------------------
        i = N_Top_Max                   !to suppress G77 warnings
        i = N_Total_Max                 !to suppress G77 warnings
        i = N_Tran_Max                  !to suppress G77 warnings
        i = N_Pts_Max                   !to suppress G77 warnings
        i = N_Stand_Max                 !to suppress G77 warnings
c--------------------------------------------------------
c       xabs   = Abs(x0(N_Total+1))     !calculated in main
c       TERang = Sqrt(xabs)             !calculated in main
        pi     = 4.* Atan(1.)
        k0     = 2.* pi/wave
        g0     = f0 * UnitCoef/TERang

c+-------------------------+
cccc    MaxOrder = 7      !| DEBUG !!!
c+-------------------------+
        Upper_Limit = 10.**MaxOrder

        SFMm5  = 0.
        ifail  = 0
        Pol    = ipol   !polarization index (sigma=1, pi=2)

c At f0=0 vacuum matrix becomes singular and cannot be inverted:
c                      |1 1|
c                      |0 0|
c Then, reflection coefficient cannot be determined and we set it to 1:

        if (Abs(g0) .lt. 1.0D-7) g0 = 1.0D-7

        if (g0 .lt. 0.) Then !---+
          SFMm5 = 1.            !|
          Return                !|
        endif  !-----------------+

        Kod = 1                         !Don't compute TTT
        Call    M1 (N_Roots,TTT,N_Fields)

c
c                     Process Vacuum
c
c+-------------------------------------------------------+
        N_used = 0                                      !|
        j_Used = 0                                      !|
        x0n    = x0(N_Used+1) / xabs                    !|
        if (ipol .eq. 1) Then  !---------+sigma          |
          Sqrt_Eps = (1.,0.)            !|               |
        else !---------------------------|pi             |
          Sqrt_Eps = Sqrt(1.+x0n*xabs)  !|               |
        endif  !-------------------------+               |
        Call    Process_Layer (g0, Sqrt_Eps,            !|
     *                         gi, Daa,                 !|not used
     *                         x0n, xhn, k0,            !|
     *                         0.,                      !|Thc(0)
     *                         TERang, S, B,            !|
     *                         u(1,N_Used+1),           !|
     *                         F(1,N_Used+1),           !|
     *                         N_Fields, ifail)         !|
        if (ifail .ne. 0) goto 100 !---------------------+-+
c+-------------------------------------------------------+ v

c
c                     Process Top Layer
c
        do    i=1,N_Total !==============================+
          N_Used = N_Used+1                             !|
c+---------------------------------------------+         |
c| Process sublayer-sublayer transition in TOP |         |
c+=============================================+         |
          if (Thickness_Tr(N_Used) .gt. 0. .AND.        !|
     *                  i_standing .eq. 0)  Then  !----+ |
c See Transi_R.FOR:                                    | |
          Call   Transi_r     (Thickness_Tr(N_Used),  !| |
     *                         Sigma(N_Used),         !| |
     *                         TranLaw, N_Tran,       !| |
     *                         x0(N_Used),            !| |
     *                         x0(N_Used+1),          !| |
     *                         xh, xh, Pol,           !| |=ipol
     *                         xabs, TERang,          !| |
     *                         g0,                    !| |
     *                         gh, gi, Daa,           !| |not used
     *                         k0,                    !| |
     *                         N_Roots, N_Fields,     !| |
     *                         N_Half,  j_Used,       !| |
     *                         S, B, B1,              !| |
     *                         u(1,N_Used+1),         !| |
     *                         u(1,N_Used),           !| |
     *                         F(1,N_Used+1),         !| |
     *                         X, TTT,                !| |
     *                         M_rr(1,1,N_Used),      !| |
     *                         M_rt(1,1,N_Used),      !| |
     *                         M_tr,                  !| |
     *                         M_tt,                  !| |
     *                         W_rr,                  !| |
     *                         W_rt,                  !| |
     *                         W_tr(1,1,N_Used),      !| |
     *                         W_tt(1,1,N_Used),      !| |
     *                         Upper_Limit, Kod,      !| |
     *                         ifail )                !| |
            if (ifail .ne. 0) goto 100  !--------------+-+-+
          endif  !-------------------------------------+ | v
                                                        !|
          j_Used = j_Used+1                             !|
          Call  MtoM (N_Roots, N_Roots, B, B1, N_Fields)!|
          x0n = x0(N_Used+1) / xabs                     !|
          if (ipol .eq. 1) Then  !---------+sigma        |
            Sqrt_Eps = (1.,0.)            !|             |
          else !---------------------------|pi           |
            Sqrt_Eps = Sqrt(1.+x0n*xabs)  !|             |
          endif  !-------------------------+             |
          Call  Process_Layer     (g0, Sqrt_Eps,        !|
     *                             gi, Daa,             !|not used
     *                             x0n, xhn, k0,        !|
     *                             Thickness(N_Used),   !|
     *                             TERang, S, B,        !|
     *                             u(1,N_Used+1),       !|
     *                             F(1,N_Used+1),       !|
     *                             N_Fields, ifail)     !|
          if (ifail .ne. 0) goto 100  !------------------+-+
          Call  Process_Interface (k0, TERang,          !| v
     *                             Sigma(N_Used),       !|
     *                             B1, S,               !|
     *                             u(1,N_Used),         !|
     *                             u(1,N_Used+1),       !|
     *                             F(1,N_Used+1),       !|
     *                             X, TTT,              !|
     *                             N_Fields,            !|
     *                             Upper_Limit,         !|
     *                             Kod, ifail)          !|
          if (ifail .ne. 0) goto 100  !------------------+-+
          Call Process_Interface2 (X, F(1,N_Used+1),    !| v
     *                             M_rr(1,1,N_Used),    !|
     *                             M_rt(1,1,N_Used),    !|
     *                             M_tr,                !|
     *                             M_tt,                !|
     *                             W_rr,                !|
     *                             W_rt,                !|
     *                             W_tr(1,1,N_Used),    !|
     *                             W_tt(1,1,N_Used),    !|
     *                             N_Fields, j_Used,    !|
     *                             ifail)               !|
          if (ifail .ne. 0) goto 100  !------------------+-+
          if (i .lt. N_Total) Then  !----------+         | v
c See Flds_Uni.for:                            |         |
            Call Backup2 (W_tr(1,1,N_Used),   !|         |
     *                    W_tt(1,1,N_Used),   !|         |
     *                    W_tr(1,1,N_Used+1), !|         |
     *                    W_tt(1,1,N_Used+1), !|         |
     *                    N_Half, N_Plus)     !|         |
          endif  !-----------------------------+         |
        enddo  !=========================================+

c
c                     Find the solution Es
c
c+-----------------------------------------------------+
        Es = W_rt(1,1)                                !|
                                                      !|
c Compute the reflection coefficient:                  |
                                                      !|
        SFMm5 = (CDabs(Es))**2                        !|
c+-----------------------------------------------------+
        if (i_standing .eq. 0) Return

c
c                  Compute the wavefields
c

        E(1,1)        = (1.0D0, 0.0D0)      !T_0
        E(2,1)        = Es                  !R_0
        E(1,N_Used+1) = W_tt(1,1,N_Used)    !T_n = W_n^tt * T_0 + W_n^tr * R_n =  W_n^tt
cc   *                * E(1,1)              ![see Eq.(24) in PRB, 57 (1998) 4829]
cc   /                / F(1,N_Used+1)
        E(2,N_Used+1) = (0.0D0, 0.0D0)      !R_n

        N = N_Plus + 1
        do      k = N_Used-1,1,-1   !========+
          Call  Find_Fields (E(1,1),        !|T_0
     *                       E(N,k+2),      !|R_{k+1}
     *                       E(1,k+1),      !|T_k
     *                       E(N,k+1),      !|R_k
     *                       W_tt(1,1,k),   !|W_k^tt
     *                       W_tr(1,1,k),   !|W_k^tr
     *                       M_rt(1,1,k+1), !|M_{k+1}^rt
     *                       M_rr(1,1,k+1), !|M_{k+1}^rr
     *                       N_Plus,ifail)  !|
          if (ifail .ne. 0) goto 100  !------+-+
        enddo  !=============================+ v

c This is the renormalization of wavefields from the
c lower interfaces, as they are given by the recursive
c equations, to the upper interface, as they are required
c for the calculations -- see Eq.(24) in PRB, 57 (1998) 4829.
c In other words, this program calculates [D_k]*[F_k^(U)].
        do      k = 2, N_Used   !============+
          Call  Norm_Fields (E(1,k),        !|T_k
     *                       E(N,k),        !|R_k
     *                       F(1,k),        !|Ft_k(i)=exp(+im*k0*TERang*u(i)*t_k), Im(u(i)>0
     *                       F(N,k),        !|Fr_k(i)=exp(-im*k0*TERang*u(i)*t_k), Im(u(i)<0
     *                       N_Plus)        !|where t_k=Thickness(k-1)
        enddo  !=============================+

        if (iDebug .ne. 0) Then !-----------------------+
          write (3,50)  g0                             !|
  50      format('The wavefields for g0/TERang=',g9.3) !|
                                                       !|
          do    k = 1, N_Used  !=============+          |
            if (k .lt. N_Used) Then !--+     |          |
              j = N_Fields            !|     |          |
            else  !--------------------|     |          |
              j = N_Plus              !|     |          |
            endif  !-------------------+     |          |
            write (3,55) k,(Abs(E(i,k+1)),  !|          | !!!!!!!!!!! DEBUG
     *                   i=1,j)             !|          | !!!!!!!!!!! DEBUG
  55        format(i5,4g12.5)               !|          | !!!!!!!!!!! DEBUG
          enddo  !===========================+          |
        endif  !----------------------------------------+

c
c                  Compute standing waves
c

        do i=1,n_standing  !=========================================+
          z = Positions(i_reference+1)                              !|
     +      + standing_range(1) + (i-1) * standing_step             !|
c Now, find in what layer we are:                                   !|
          do    k=N_Used,1,-1  !==============+                      |
            if (Positions(k) .le. z) goto 5 !-+-+                    |
          enddo   !===========================+ |                    |
          k  = 0                               !|                    |
  5       continue    !<------------------------+                    |
                                                                    !|
          kk = k                                     !Vacuum fields !|--+exchange these lines
          if (k .lt. 1) k=1                                         !|  |to use vacuum or subs-
ccc       kk = k                                     !Substrate flds!|<-+trate fields in vacuum
c (NOTE: even when using substrate fields, the exponents must belong |
c to the vacuum fields because otherwise we may get huge numbers in  |
c the vacuum! -- SS 2005-07-29)                                      |
                                                                    !|
c The offset from the closest upper interface:                      !|
          zi = z - Positions(k)                                     !|
          SUMM = (0.0D0, 0.0D0)                                     !|
          do j=1,N_Fields !========================================+ |
            if (abs(E(j,kk+1)).gt.1.E-32) Then !-----------------+ | |
              ikuz = Im * k0*u(j,kk+1)*TERang * zi              !| | |
              rkz = dReal(ikuz)                                 !| | |
              ikz = dImag(ikuz)                                 !| | |
              if     (ikz .gt. (2.*pi)) Then !---+               | | |
                j2  = int(ikz/(2.*pi))          !|               | | |
                ikz = ikz - (2.*pi)*j2          !|               | | |
              elseif (ikz .lt. -(2.*pi)) Then !--|               | | |
                j2  =-int(ikz/(2.*pi))          !|               | | |
                ikz = ikz + (2.*pi)*j2          !|               | | |
              endif  !---------------------------+               | | |
              if     (rkz .gt.  1.D+10) Then !-----------------+ | | |
                write(3,*) 'Warning: too big positive i*kz*z' !| | | |
                rkz = 1.E+10                                  !| | | |
              elseif (rkz .lt. -1.D+10) Then !-----------------| | | |
                write(3,*) 'Warning: too big negative i*kz*z' !| | | |
                rkz = -1.E+10                                 !| | | |
              endif !------------------------------------------+ | | |
              ikuz = dCmplx(rkz,ikz)                            !| | |
              SUMM = SUMM + E(j,kk+1)*exp(ikuz)                 !| | |
            endif  !---------------------------------------------+ | |
          enddo  !=================================================+ |
          SW(i) = abs(SUMM)**2                                      !|
        enddo  !=====================================================+

        Return

c========
c ERRORS:
c========
  100   continue
        write (3,*)     'SFMm5: failure code=',ifail,
     *                  '-- at f0: ',f0
        Return
        end

c =======================================================
c In the case when we need
c wavefields, this is defined
c in Flds_Uni.for:

c       Subroutine Algorithm_Type  (Type)
c       Character       Type*(*)
c       Type = 'Recursive'             !10
c cc    Type = 'Matrix'
c       return
c       end
