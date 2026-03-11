c
c rfms_uuu.for [NAG version, replaces rfms_sss.for for gid_slm7]
c ==============================================================
c
c +========+     +========+         +=======================+
c |gid_slm7|-->--|rfmm7_rx|-->------|rfms_urd7 (GID_slm7)   |-->-+
c +--------+     +--------+         +-----------------------+    |
c +========+     +========+                                      |
c |gid_slm5|-->--|rfmm5_rx|->-+                                  |
c +--------+     +--------+   |                                  |  WE ARE HERE
c +========+     +========+   |     +=======================+    |   +========+
c |ter_slm5|-->--|sfmm5_r |->-+-->--|rfms_urd (GID,TER,TRDS)|->--+->-|RFMs_uuu|
c +--------+     +--------+   |     +-----------------------+    |   +--------+
c +========+     +========+   |                                  |
c |trds_99 |-->--|trsu2r97|->-+                                  |
c +--------+     +--------+                                      |
c +========+     +========+                                      |
c |mag_sl98|-->--|mmre98_r|->-+     +=======================+    |
c +--------+     +--------+   |-->--|rfms_umd (MAG)         |->--+
c |mag_sl99|-->--|mmre98_r|->-+     +-----------------------+
c +--------+     +--------+
c
c     These unified programs implement some assistance
c      operations for RFM's called from < gid_slm7 >.
c                                         ========
c      They are supposed to be more generic and thus replace
c      the pograms contained in rfms_sss.for which are used
c      in the following REFLECTION, MAGNETIC REFLECTION and 
c                   DIFFRACTION routines:
c
c         < ter_sl10, ter_sl12, ter_sl97, trds_97 >
c                       < mag_sl99 >
c     < gid_sl32, gid_sl50, gid_sl55, gids_10, gids_97 >
c
c The module is called by the above programs via either
c 1. RMFs_urd.for (reflection & diffraction)  or
c 2. RFMs_umg.for (resonant magnetic reflection)
c    (-- in the project stage!!! Meanwhile see MMRE99_R.for)
c
c          The x-ray diffaction coefficient under GID
c          or  x-ray reflection coefficient under TER
c          or  x-ray reflection coefficient under mag_TER
c                         is computed
c           with the account for interface roughness.
c
c +======================================================+(old RFMs_uni.for)
c | 1a. Subroutine:     Process_Interface                |
c | 2a. Subroutine:     Process_Interface_no_x           |
c | 3a. Subroutine:     Inverse_Matrix                   |
c | 4a. Subroutine:     Gauss_Solution ----> F04ADF (NAG)|or ZGESV (LAPACK)
c | 5a. Subroutine:     Copy_Layer                       |
c |------------------------------------------------------|(old RFMs_rec.for)
c | 1b. Subroutine:     Process_Interface2               |
c | 2b. Subroutine:     Inverse_Matrix_simple            |
c | 3b. Subroutine:     Backup4                          |
c |------------------------------------------------------|(old Mc16_sub.for)
c | 1c. Subroutine:     MxM                              |
c | 2c. Subroutine:     M0                               |
c | 2c. Subroutine:     M1                               |
c | 4c. Subroutine:     MtoM                             |
c | 5c. Subroutine:     MpmM                             |
c | 6c. Real*8 function Element_Max                      |
c | 7c. Subroutine:     Sort_Roots                       |
c |------------------------------------------------------|(rfms_uuu.for)
c | 1d. Subroutine:     Process_Interface2_rect          |
c | 2d. Subroutine:     Backup4_rect                     |
c | 3d. Subroutine:     MtoM_rect                        |
c | 4d. Subroutine:     MpmM_rect                        |
c +======================================================+
c
c
c
c ---------------------------------------------------------
c                       RFMs_uni.for
c ####################################################### 1a

        Subroutine  Process_Interface (k0, TERang, Sigma,
     *                                 B1, S, u1, u,
     *                                 F, X, TTT, N_Fields,
     *                                 Upper_Limit,
     *                                 Kod, ifail)
        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields

        Complex*16      B1(N_Dim,N_Fields),        ! Previous layer 1/S
     *                  S(N_Dim,N_Fields),         ! layer S-matrix
     *                  u1(N_Fields),              ! Previous layer roots
     *                  u(N_Fields),               ! layer roots
     *                  F(N_Fields),               ! layer F-matrix
     *                  TTT(N_Dim,N_Fields),       ! product of operators
     *                  X(N_Dim,N_Fields),         ! X-matrix of interface
     *                  T(N_Dim,N_Dim),            ! work matrix
     *                  q                          ! work cell

        Real*8          k0,                        ! k0 = 2*pi/wave
     *                  S_Max,                     ! max element of T1*..Tn
     *                  Upper_Limit                ! max matr.element (lim)

        Real*4          TERang,                    ! critical TER angle
     *                  Sigma                      ! rms roughness heigth

        Integer         N_Plus, ifail, Kod,
     *                  i, j

        Real*8          Element_Max                ! max matr.element (fun)
        External        Element_Max

        Call Process_Interface_no_x (k0, TERang, Sigma,
     *                               B1, S, u1, u,
     *                               X, N_Fields)

        ifail = 0
c+----------------------------------+
c|No multiplication of matrices:    |e.g. we return here when we use
        if (Kod.ne.0)   Return     !|the recursive matrix technique
c+----------------------------------+

c Compute T=TTT*X:
        Call    MxM (N_Dim,   N_Dim,   N_Dim,
     *               TTT,     X,       T,
     *               N_Fields,N_Fields,N_Fields)

c Compute (TTT*X)*F(t):
        N_Plus = N_Fields / 2
        do    j=1,N_Fields  !=============+
c The growing exponents were inverted:    |
           if (j.le.N_Plus) Then  !---+   |
              q = 1.0D0 / F(j)       !|   |
           else  !--------------------+   |
              q = F(j)               !|   |
           endif  !-------------------+   |
           do  i=1,N_Fields  !=========+  |
              TTT(i,j) = T(i,j) * q   !|  |
           enddo !=====================+  |
        enddo   !=========================+

        S_Max = Element_Max (N_Dim,    TTT,
     *                       N_Fields, N_Fields)

        if (S_Max .gt. Upper_Limit)  Then  !---+
           Call MtoM (N_Dim, N_Dim,           !|
     *                T, TTT,  N_Fields)      !| TTT=T
           Kod   = 1                          !|
           ifail =-1                          !| Overflow! We must proceed
        endif  !-------------------------------+ to the recursive method!
                                                !(see Process_Interface2).
        Return
        End

c ####################################################### 2a

        Subroutine  Process_Interface_no_x (k0, TERang, Sigma,
     *                                      B1, S, u1, u,
     *                                      X, N_Fields)
        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields

        Complex*16      B1(N_Dim,N_Fields),        ! Previous layer 1/S
     *                  S(N_Dim,N_Fields),         ! layer S-matrix
     *                  u1(N_Fields),              ! Previous layer roots
     *                  u(N_Fields),               ! layer roots
     *                  X(N_Dim,N_Fields),         ! X-matrix of interface
     *                  q                          ! work cell

        Real*8          k0,                        ! k0 = 2*pi/wave
     *                  p, p2, q1                  ! work cells

        Real*4          TERang,                    ! critical TER angle
     *                  Sigma                      ! rms roughness heigth

        Integer         i, j

c+--------------------------------------------------+
c|Compute X-matrix of layer upper interface:        |
c|                     |    |-1                     |
c|            X      = |S   |  * S                  |
c|             k-1,k   | k-1|     k                 |
c+----------------------------------                |
                                                   !|
        Call    MxM (N_Dim,   N_Dim,   N_Dim,      !|
     *               B1,      S,       X,          !|
     *               N_Fields,N_Fields,N_Fields)   !|
c---------------------------------------------------+

c Average X-matrix over interface roughness:

        if (Sigma.gt.0.) Then  !-------------------------+
           p = k0*Sigma*TERang                          !|
           p2 = p * p                                   !|
                                                        !|
           do    i=1,N_Fields  !=======================+ |
              do j=1,N_Fields  !=====================+ | |
                                                    !| | |
                 q = u(j)-u1(i)                     !| | |initial: high peak at b<1
c                q = u(j)-Dconjg(u1(i))             !| | |2007.10.06 no effecr
c                q = u(i)-u1(j)                     !| | |2007.10.06 mess
c                q = u(j)-Dconjg(u(i))              !| | |2007.10.06 zero outside peak
c                q = u(j)-u(i)                      !| | |2007.10.06 zero outside peak
c                q = u1(j)-u1(i)                    !| | |2007.10.06 low peak at b<1
                 q1 = CDabs(q)                      !| | |
                                                    !| | |
                 if (q1 .gt. 1.D-04) Then !--------+ | | |
                                                  !| | | |
                    q  = p2*(q*q)/2.              !| | | |
                    q  = CDexp(-q)                !| | | |exp(-0.5*k0*Sigma*TERang*(u(j)-u1(i))^2)
                    X(i,j) = X(i,j) * q           !| | | |
                                                  !| | | |
c                   p = p2*(q1*q1)/2.             !| | | |2007.10.06 no effect
c                   if (p .gt. pq_max) p=pq_max   !| | | |2007.10.06
c ccccccc           p = Dexp(+p)                  !| | | |2007.10.06
c                   p = Dexp(-p)                  !| | | |2007.10.06
c                   X(i,j) = X(i,j) * p           !| | | |2007.10.06
                                                  !| | | |
                 endif  !--------------------------+ | | |
                                                    !| | |
              enddo !================================+ | |
           enddo   !===================================+ |
        endif  !-----------------------------------------+

        return
        End

c ####################################################### 3a

        Subroutine      Inverse_Matrix  (A,B,N_Fields,ifail)
c B=1/A (square)
        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields

        Complex*16      A(N_Dim,N_Fields),
     *                  B(N_Dim,N_Fields),
     *                  C(N_Dim,N_Dim),
     *                  D(N_Dim,N_Dim)
        Integer         ifail

        ifail = 0
c Backup C=A
        Call    MtoM (N_Dim,N_Dim,A,C,N_Fields)
c Fill D=I
        Call    M1   (N_Dim,D,N_Fields)

        Call    Gauss_Solution  (C, D, B,
     *                           N_Fields,
     *                           N_Fields,
     *                           ifail)

        return
        end

c ####################################################### 4a

        Subroutine      Gauss_Solution  (A, B, X, 
     *                                   N_Fields,
     *                                   Nright,
     *                                   ifail)

        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields, Nright

        Complex*16      A(N_Dim,N_Fields),
     *                  B(N_Dim,Nright),
     *                  X(N_Dim,Nright)

c       Real*8          wkspce(N_Dim)                   !for NAG F04ADF
        
        Integer         IPIV(N_Dim)                     !for LAPACK ZGESV

        Integer         ifail
        
        Logical         use_NAG/.False./


        ifail = 0
        if (use_NAG) Then !--------------------------------------------+
c F04ADF calculates the approximate solution of a set of complex       |
c linear equations with multiple right-hand sides, using an LU         |
c factorization with partial pivoting.                                 |
c          Call F04ADF (A, N_Dim, B, N_Dim, N_Fields,                 !|
c    *                  Nright, X, N_Dim, wkspce, ifail)              !|
        else  !--------------------------------------------------------+
c Use the LAPACK routine ZGESV which computes the solution to          |
c system of linear equations A*X=B for GE matrices (simple driver)     |
c https://www.netlib.org/lapack/explore-html/d1/ddc/zgesv_8f.html      |
           Call ZGESV (N_Fields, Nright, A, N_Dim, IPIV, B, N_Dim,    !|
     *                 ifail)                                         !|
           if (ifail .eq. 0) Then !---------------------------------+  |
              Call MtoM_rect(N_Dim, N_Dim, B, X, N_Fields, Nright) !|  |
           endif !--------------------------------------------------+  |
        endif !--------------------------------------------------------+
        return
        end

c ####################################################### 5a

        Subroutine      Copy_Layer (B,B1,u,u1,N_Fields)

        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields

        Complex*16      B(N_Dim,N_Fields),
     *                  B1(N_Dim,N_Fields),
     *                  u(N_Fields),
     *                  u1(N_Fields)
        Integer         i

        Call    MtoM    (N_Dim, N_Dim, B, B1, N_Fields)

        do      i=1,N_Fields  !==+
          u1(i) = u(i)          !|
        enddo   !================+
        Return
        End

c
c ---------------------------------------------------------
c                       RFMs_rec.for (used in case N_t=N_r)
c ####################################################### 1b

        Subroutine Process_Interface2 (X, F,
     *                                 M_rr, M_rt,
     *                                 M_tr, M_tt,
     *                                 W_rr, W_rt,
     *                                 W_tr, W_tt,
     *                                 N_Fields,
     *                                 LayerNo,
     *                                 ifail)
        Integer         N_Dim, N_Half
        Parameter      (N_Dim  = 4,
     *                  N_Half = 2)                ! N_Half = N_Dim/2

        Integer         LayerNo,                   ! Layer number in the stack
     *                  N_Fields,
     *                  N_t, N_r                   ! Number of transmitted and reflected waves

        Complex*16      F(N_Fields),               ! layer F-matrix
     *                  X(N_Dim,N_Dim),            ! X-matrix of interface
     *                  M_tt(N_Half,N_Half),       ! operator of interface
     *                  M_rt(N_Half,N_Half),       ! operator of interface
     *                  M_tr(N_Half,N_Half),       ! operator of interface
     *                  M_rr(N_Half,N_Half),       ! operator of interface
     *                  W_tt(N_Half,N_Half),       ! product of operators
     *                  W_rt(N_Half,N_Half),       ! product of operators
     *                  W_tr(N_Half,N_Half),       ! product of operators
     *                  W_rr(N_Half,N_Half)        ! product of operators

        Integer         ifail

        N_t = N_Fields / 2                      ! this is for N_t=N_r
        N_r = N_Fields / 2                      ! this is for N_t=N_r
        Call Process_Interface2_rect (X, F,
     *                                M_rr, M_rt, M_tr, M_tt,
     *                                W_rr, W_rt, W_tr, W_tt,
     *                                N_t, N_r, LayerNo, ifail)
        Return
        End

c ####################################################### 2b

        Subroutine      Inverse_Matrix_simple (A, B,
     *                                         N_Dim,
     *                                         N_Fields,
     *                                         ifail)
c B=1/A for N_Fields=1,2

        Integer         N_Dim, N_Fields

        Complex*16      A(N_Dim,N_Fields),
     *                  B(N_Dim,N_Fields),
     *                  C1, C2, Det
        Integer         ifail

        ifail = 0

        if (N_Fields.eq.1) Then !-----------------+
                                                 !|
           Det = A(1,1)                          !|
           if (abs(Det).lt.1.E-32) Then !---+     |
              ifail = 201                  !|     |
           else   !-------------------------+     |
              B(1,1) = 1.0D0 / Det         !|     |
           endif  !-------------------------+     |
                                                 !|
        elseif (N_Fields.eq.2)  Then  !-----------+
                                                 !|
           Det = A(1,1)*A(2,2)-A(1,2)*A(2,1)     !|
           if (abs(Det).lt.1.E-32) Then !---+     |
              ifail = 202                  !|     |
           else   !-------------------------+     |
              Det = 1.0D0 / Det            !|     |
              C1     = A(1,1)              !|     |
              C2     = A(2,2)              !|     |
              B(1,1) =   C2   * Det        !|     |
              B(1,2) =-A(1,2) * Det        !|     |
              B(2,1) =-A(2,1) * Det        !|     |
              B(2,2) =   C1   * Det        !|     |
           endif  !-------------------------+     |
                                                 !|
        else   !----------------------------------+
                                                 !|
           ifail = 200                           !|
                                                 !|
        endif  !----------------------------------+

        return
        end

c ####################################################### 3b

        Subroutine Backup4 (W_rr,W_rt,W_tr,W_tt,
     *                      X_rr,X_rt,X_tr,X_tt,
     *                      N_Half,N_Plus)
c X_rr=W_rr
c X_rt=W_rt
c X_tr=W_tr
c X_tt=W_tt

        Integer         N_Half, N_Plus

        Complex*16      W_rr(N_Half,N_Plus),
     *                  W_rt(N_Half,N_Plus),
     *                  W_tr(N_Half,N_Plus),
     *                  W_tt(N_Half,N_Plus),
     *                  X_rr(N_Half,N_Plus),
     *                  X_rt(N_Half,N_Plus),
     *                  X_tr(N_Half,N_Plus),
     *                  X_tt(N_Half,N_Plus)

        Integer         i, j

        do    i=1,N_Plus    !=======+
        do    j=1,N_Plus    !=====+ |
          X_rr(i,j) = W_rr(i,j)  !| |
          X_rt(i,j) = W_rt(i,j)  !| |
          X_tr(i,j) = W_tr(i,j)  !| |
          X_tt(i,j) = W_tt(i,j)  !| |
        enddo   !=================+ |
        enddo   !===================+

        return
        end

c ---------------------------------------------------------
c                Mc16_sub.for /modified/
c
c ATTENTION: The number of parameters in routines is
c            different from that in mc16_sub.for !!!
c ####################################################### 1c

        Subroutine MxM (Na, Nb, Nc, 
     *                  A,  B,  C, 
     *                  mx, my, mm)
c *** C(mx,my) = A(mx,mm)*B(mm,my)

        Integer         Na, Nb, Nc, mx, my, mm,
     *                  i, j, k
        Complex*16      A(Na,mm), B(Nb,my), C(Nc,my)

        if (mx.gt.Nc .OR. mx.gt.Na .OR. mm.gt.Nb) Then !-------+
           write (0,*) 'MxM: Na=',Na,'  Nb=',Nb,'  Nc=',Nc    !|
           write (0,*) 'MxM: mm=',mm,'  mx=',mx,'  my=',my    !|
           stop 'MxM: oversize'                               !|
        endif !------------------------------------------------+

        Call    M0 (Nc,C,mx,my)

        do      i=1,mx  !========================+
        do      j=1,my  !======================+ |
        do      k=1,mm  !====================+ | |
           C(i,j) = C(i,j) + A(i,k)*B(k,j)  !| | |
        enddo  !=============================+ | |
        enddo  !===============================+ |
        enddo  !=================================+

        return
        end

c ======================================================= 2c

        Subroutine      M0 (Na, A, mx, my)
c M=0 (size: mx x my)
        Integer         Na, mx, my,
     *                  i, j
        Complex*16      A(Na,my)

        if (mx.gt.Na) Then !---------------------------------+
           write (0,*) 'M0: Na=',Na,'  mx=',mx,'  my=',my   !|
           stop 'M0: oversize'                              !|
        endif !----------------------------------------------+

        do      i=1,mx  !===============+
        do      j=1,my  !============+  |
           A(i,j) = (0.0D0, 0.0D0)  !|  |
        enddo  !=====================+  |
        enddo  !========================+

        return
        end

c ======================================================= 3c

        Subroutine      M1 (Na, A, m)
c M=1 (size: m x m)
        Integer         Na, m,
     *                  i
        Complex*16      A(Na,m)

        if (m.gt.Na) Then !-----------------------+
           write (0,*) 'M1: Na=',Na,'  m=',m     !|
           stop 'M1: oversize'                   !|
        endif !-----------------------------------+

        Call    M0 (Na,A,m,m)

        do      i=1,m  !===============+
           A(i,i) = (1.0D0, 0.0D0)    !|
        enddo  !=======================+

        return
        end

c ======================================================= 4c

        Subroutine      MtoM (Na, Nb, A, B, m)
c B = A (size: m x m)
        Integer         Na, Nb, m,
     *                  i, j
        Complex*16      A(Na,m), B(Nb,m)

        if (m.gt.Na .OR. m.gt.Nb) Then !-----------------------+
           write (0,*) 'MtoM: Na=',Na,'  Nb=',Nb,'  m=',m     !|
           stop 'MtoM: oversize'                              !|
        endif !------------------------------------------------+

        do      i=1,m  !===========+
        do      j=1,m  !========+  |
           B(i,j) = A(i,j)     !|  |
        enddo  !================+  |
        enddo  !===================+

        return
        end

c ======================================================= 5c

        Subroutine      MpmM (Na, Nb, Nc, A, B, C, m, mp)
c C = A +- B (size: m x m)
c mp=1 or mp=-1
        Integer         Na, Nb, Nc, m, mp,
     *                  i, j
        Complex*16      A(Na,m), B(Nb,m), C(Nc,m)

        if (m.gt.Na .OR. m.gt.Nb) Then !-----------------------+
           write (0,*) 'MpmM: Na=',Na,'  Nb=',Nb,'  m=',m     !|
           stop 'MpmM: oversize'                              !|
        endif !------------------------------------------------+

        do      i=1,m  !====================+
        do      j=1,m  !=================+  |
           C(i,j) = A(i,j) + mp*B(i,j)  !|  |
        enddo  !=========================+  |
        enddo  !============================+

        return
        end

c ======================================================= 6c

        Real*8  Function  Element_Max (Na, A, mx, my)
c Element_Max = Max(abs(A)), (size: mx x my)

        Integer         Na, mx, my,
     *                  i, j
        Complex*16      A(Na,my)
        Real*8          s

        Element_Max = 0.0D0
        do i=1,mx  !=====================================+
        do j=1,my  !==================================+  |
           s = Abs (A(i,j))                          !|  |
           if (s .gt. Element_Max)  Element_Max = s  !|  |
        enddo  !======================================+  |
        enddo  !=========================================+
        return
        end

c ======================================================= 7c

        Subroutine Sort_Roots (ur,ui,n)

        Integer         n,
     *                  i, j
        Real*8          ur(n), ui(n),
     *                  vr,    vi
c--------------------------------------------------------
c  Sort roots in descending of Im(u).
c--------------------------------------------------------
        do      j = 2,n  !====================+
           vi = ui(j)                        !|
           vr = ur(j)                        !|
                                             !|
           do    i = j-1,1,-1  !=========+    |
              if (ui(i).ge.vi)  goto 3 !-+->+ |
              ui(i+1) = ui(i)           !|  | |
              ur(i+1) = ur(i)           !|  | |
           enddo  !======================+  V |
                                           !| |
           i = 0                           !| |
  3        continue      !<-----------------+ |
           ui(i+1) = vi                      !|
           ur(i+1) = vr                      !|
        enddo  !==============================+

        return
        end

c ---------------------------------------------------------
c                       RFMs_uuu.for
c ####################################################### 1d

        Subroutine Process_Interface2_rect (X, F,
     *                                      M_rr, M_rt,
     *                                      M_tr, M_tt,
     *                                      W_rr, W_rt,
     *                                      W_tr, W_tt,
     *                                      N_t, N_r,
     *                                      LayerNo,
     *                                      ifail)
        Integer         N_Dim, N_Half
        Parameter      (N_Dim  = 4,
     *                  N_Half = 2)                ! N_Half = N_Dim/2

        Integer         LayerNo,                   ! Layer number in the stack
     *                  N_t, N_r                   ! Number of transmitted and reflected waves

        Complex*16      F(N_t+N_r),                ! layer F-matrix
     *                  X(N_Dim,N_Dim),            ! X-matrix of interface
     *                  M_tt(N_Half,N_Half),       ! operator of interface
     *                  M_rt(N_Half,N_Half),       ! operator of interface
     *                  M_tr(N_Half,N_Half),       ! operator of interface
     *                  M_rr(N_Half,N_Half),       ! operator of interface
     *                  W_tt(N_Half,N_Half),       ! product of interface operators
     *                  W_rt(N_Half,N_Half),       ! product of interface operators
     *                  W_tr(N_Half,N_Half),       ! product of interface operators
     *                  W_rr(N_Half,N_Half),       ! product of interface operators
     *                  X_tt(N_Half,N_Half),       ! 1/4 blocks of X-matrix (work matrix)
     *                  X_rt(N_Half,N_Half),       ! 1/4 blocks of X-matrix (work matrix)
     *                  X_tr(N_Half,N_Half),       ! 1/4 blocks of X-matrix (work matrix)
     *                  X_rr(N_Half,N_Half),       ! 1/4 blocks of X-matrix (work matrix)
     *                  T_1(N_Half,N_Half),        ! work matrix
     *                  T_2(N_Half,N_Half)         ! work matrix

        Integer         ifail,
     *                  i, j

        ifail  = 0
        if (N_t .gt. N_Half .OR. N_r .gt. N_Half) Then !--------+
           write (0,*) 'Process_Interface2_rect: N_t=',N_t,    !|
     *                                         ' N_r=',N_r,    !|
     *                                      ' N_half=',N_Half  !|
           stop 'Process_Interface2_rect: oversize'            !|
        endif !-------------------------------------------------+
c+----------------------------------------------------+
c|Find 1/4 matrices of X:                             |
c|       | X_tt | X_rt |                              |
c|   X = | -----+----- |                              |
c|       | X_tr | X_rr |                              |
c+--------------------------                          |
c See PRB v.57, 4829-4841 (1998), Eq.(16):            |
        do      i=1,N_t  !==================+         |
        do      j=1,N_t  !================+ |         |
           X_tt(i,j) = X(i,j)            !| |         |t x t
        enddo  !==========================+ |         |
        do      j=1,N_r  !================+ |         |
           X_tr(i,j) = X(i,j+N_t)        !| |         |t x r
        enddo  !==========================+ |         |
        enddo  !============================+         |
c                                                     |
        do      i=1,N_r  !==================+         |
        do      j=1,N_t  !================+ |         |
           X_rt(i,j) = X(i+N_t,j)        !| |         |r x t 
        enddo  !==========================+ |         |
        do      j=1,N_r  !================+ |         |
           X_rr(i,j) = X(i+N_t,j+N_t)    !| |         |r x r
        enddo  !==========================+ |         |
        enddo  !============================+         |
c-----------------------------------------------------+

c+----------------------------------------------------+
c|Find M-matrices of the interface:                   |
c+--------------------------------                    |
c X_tt = 1/(X_tt)                                     |
        Call    Inverse_Matrix_simple (X_tt,   X_tt, !|
     *                                 N_Half, N_t,  !|
     *                                 ifail)        !|t x t
        if (ifail.ne.0) return                       !|
                                                     !|
                                                     !|
c M_rt = X_rt * (1/X_tt)                              | Eq.(18)
        Call    MxM (N_Half, N_Half, N_Half,         !|
     *               X_rt,   X_tt,   M_rt,           !|
     *               N_r,    N_t,    N_t)            !|r x t
                                                     !|
c T_1 = M_rt * X_tr                                   |
        Call    MxM (N_Half, N_Half, N_Half,         !|
     *               M_rt,   X_tr,   T_1,            !|
     *               N_r,    N_r,    N_t)            !|r x r
                                                     !|
c T_2 = X_rr - [M_rt * X_tr] = X_rr - T_1             | Eq.(18)
        Call    MpmM (N_Half, N_Half, N_Half,        !|
     *                X_rr,   T_1,    T_2,           !|
     *                N_r,    -1)                    !|r x r
                                                     !|
c M_tt = (1/F_+)*(1/X_tt)                             | Eq.(18)
        do    i=1,N_t  !========================+     |
        do    j=1,N_t  !=====================+  |     |
           M_tt(i,j)= F(i)      * X_tt(i,j) !|  |     |
        enddo   !============================+  |     |t x t
        enddo   !===============================+     |
c T_1  = - X_tr * (F_-)                               |
        do    i=1,N_t  !========================+     |
        do    j=1,N_r  !=====================+  |     |
           T_1(i,j) =-X_tr(i,j) * F(j+N_t)  !|  |     |t x r
        enddo   !============================+  |     |
        enddo   !===============================+     |
c M_rr = {X_rr - [M_rt * X_tr]} * (F_-) = T_2 * (F_-) | Eq.(18)
        do    i=1,N_r  !========================+     |
        do    j=1,N_r  !=====================+  |     |
           M_rr(i,j)= T_2(i,j)   * F(j+N_t) !|  |     |r x r
        enddo   !============================+  |     |
        enddo   !===============================+     |
                                                     !|
c M_tr = - M_tt * [X_tr * (F_-)] = M_tt * T_1         |
        Call    MxM (N_Half, N_Half, N_Half,         !|
     *               M_tt,   T_1,    M_tr,           !|
     *               N_t,    N_r,    N_t)            !|t x r
c+----------------------------------------------------+

        if (LayerNo.lt.1) Then !---+
          ifail = -100            !|
          return                  !|
        endif  !-------------------+

        if (LayerNo.eq.1) Then !----------------------+
c+-------------------------------+                    |
c|If the 1st layer: M_ij -> W_ij |                    |
c+-------------------------------+                    |
           Call  Backup4_rect (M_rr,M_rt,M_tr,M_tt,  !|
     *                         W_rr,W_rt,W_tr,W_tt,  !|
     *                         N_Half,N_t,N_r)       !|
           return                                    !|
        endif  !--------------------------------------+

c+----------------------------------------------------+
c|Iterative calculation of W-matrices of the interface|
c| (the product of M-matrices)                        |
c+--------------------------------                    |
c Backup previous W-matrices: W_ij --> X_ij           |
        Call    Backup4_rect (W_rr,W_rt,W_tr,W_tt,   !|
     *                        X_rr,X_rt,X_tr,X_tt,   !|
     *                        N_Half,N_t,N_r)        !|
                                                     !|
c Eq.(23): calculating B_n                            |
c T_1 = 1                                             |
        Call    M1  (N_Half,T_1,N_r)                 !|r x r
                                                     !|
c T_2 = M_rt * X_tr   (X_tr is W_tr(n))               |(M_rt * W_tr(n)) for B_n
        Call    MxM (N_Half, N_Half, N_Half,         !|
     *               M_rt,   X_tr,   T_2,            !|
     *               N_r,    N_r,    N_t)            !|r x r
                                                     !|
c T_1 = 1 - M_rt * X_tr = 1 -T_2 (X_tr is W_tr(n))    |(1 - M_rt * W_tr(n)) for B_n
        Call    MpmM (N_Half, N_Half, N_Half,        !|
     *                T_1,    T_2,    T_1,           !|
     *                N_r,    -1)                    !|r x r
                                                     !|
c T_1 = 1/(1 - M_rt * X_tr) = 1/T_1 (X_tr is W_tr(n)) |1/(1 - M_rt * W_tr(n)) for B_n
        Call    Inverse_Matrix_simple (T_1,    T_1,  !|
     *                                 N_Half, N_r,  !|
     *                                 ifail)        !|r x r
        if (ifail.ne.0) return                       !|
                                                     !|
c T_2 = X_rr * [1/(1 - M_rt * X_tr)] = X_rr * T_1     |B_n = W_rr(n)/(1 - M_rt * W_tr(n))
        Call    MxM (N_Half, N_Half, N_Half,         !|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     *               X_rr,   T_1,    T_2,            !|
     *               N_r,    N_r,    N_r)            !|r x r
                                                     !|
c Eq.(22): calculating W_rr for layer (n+1)           |
c W_rr = {X_rr*[1/(1-M_rt*X_tr)]} * M_rr = T_2 * M_rr |W_rr(n+1) = B_n * M_rr
        Call    MxM (N_Half, N_Half, N_Half,         !|
     *               T_2,    M_rr,   W_rr,           !|
     *               N_r,    N_r,    N_r)            !|r x r
                                                     !|
c T_1 = M_rt * X_tt  (X_tt is W_tt)                   |
        Call    MxM (N_Half, N_Half, N_Half,         !|M_rt * W_tt(n) for Eq.(22)
     *               M_rt,   X_tt,   T_1,            !|
     *               N_r,    N_t,    N_t)            !|r x t
                                                     !|
c W_rt = B_n * (M_rt * X_tt) = T_2 x T_1              |
        Call    MxM (N_Half, N_Half, N_Half,         !|B_n * M_rt * W_tt(n)
     *               T_2,    T_1,    W_rt,           !|
     *               N_r,    N_t,    N_r)            !|r x t
                                                     !|
c Eq.(22): calculating W_rr for layer (n+1)           |
c W_rt = X_rt + B_n * (M_rt * X_tt)                  !|
        Call    MpmM_rect (N_Half, N_Half, N_Half,   !|W_rt(n+1) = W_rt(n) + B_n * M_rt * W_tt(n)
     *                     X_rt,   W_rt,   W_rt,     !|
     *                     N_r,    N_t,  +1)         !|r x t
                                                     !|
c--------------------                                 |
c Eq.(23): calculating A_n                            |
c T_1 = 1                                             |
        Call    M1  (N_Half, T_1, N_t)               !|t x t
                                                     !|
c T_2 = X_tr * M_rt   (X_tr is W_tr(n))               |
        Call    MxM (N_Half, N_Half, N_Half,         !|W_tr(n) * M_rt
     *               X_tr,   M_rt,   T_2,            !|
     *               N_t,    N_t,    N_r)            !|t x t
                                                     !|
c T_1 = 1 - X_tr * M_rt = 1 - T_2  (X_tr is W_tr(n))  |
        Call    MpmM (N_Half, N_Half, N_Half,        !|1 - W_tr(n) * M_rt
     *                T_1,    T_2,    T_1,           !|
     *                N_t,    -1)                    !|t x t
                                                     !|
c T_1 = 1/(1 - X_tr * M_rt)   (X_tr is W_tr(n))       |
        Call    Inverse_Matrix_simple (T_1,    T_1,  !|1/(1 - W_tr(n) * M_rt)
     *                                 N_Half, N_t,  !|
     *                                 ifail)        !|t x t
        if (ifail.ne.0) return                       !|
                                                     !|
c T_2 = M_tt * [1/(1 - X_tr * M_rt)] = M_tt * T_1     |
        Call    MxM (N_Half, N_Half, N_Half,         !|A_n = M_tt/(1 - W_tr(n) * M_rt)
     *               M_tt,   T_1,    T_2,            !|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     *               N_t,    N_t,    N_t)            !|t x t
                                                     !|
c Eq.(22): calculating W_tt for layer (n+1)           |
c W_tt = A_n * X_tt   (X_tt is W_tt(n))               |
        Call    MxM (N_Half, N_Half, N_Half,         !|W_tt(n+1) = A_n * W_tt(n)
     *               T_2,    X_tt,   W_tt,           !|
     *               N_t,    N_t,    N_t)            !|t x t
                                                     !|
c T_1 = X_tr * M_rr   (X_tr is W_tr(n)                |
        Call    MxM (N_Half, N_Half, N_Half,         !|W_tr(n) * M_rr
     *               X_tr,   M_rr,   T_1,            !|
     *               N_t,    N_r,    N_r)            !|t x r
                                                     !|
c W_tr = A_n * (X_tr * M_rr) = T_2 * T_1              |
        Call    MxM (N_Half, N_Half, N_Half,         !|A_n * (W_tr(n) * M_rr)
     *               T_2,    T_1,    W_tr,           !|
     *               N_t,    N_r,    N_t)            !|t x r
                                                     !|
c Eq.(22): calculating W_tr for layer (n+1)           |
c W_tr = M_tr + A_n * (X_tr * M_rr)                   |
        Call    MpmM_rect (N_Half, N_Half, N_Half,   !|W_tr(n+1) = M_tr + A_n * (W_tr(n) * M_rr)
     *                     M_tr,   W_tr,   W_tr,     !|
     *                     N_t,    N_r,   +1)        !|t x r
c+----------------------------------------------------+
        Return
        End

c ======================================================= 2d

        Subroutine Backup4_rect (W_rr,W_rt,W_tr,W_tt,
     *                           X_rr,X_rt,X_tr,X_tt,
     *                           N_Half,N_t,N_r)
c X_rr=W_rr
c X_rt=W_rt
c X_tr=W_tr
c X_tt=W_tt

        Integer         N_Half, N_t, N_r

        Complex*16      W_tt(N_Half,N_t), X_tt(N_Half,N_t),
     *                  W_rt(N_Half,N_t), X_rt(N_Half,N_t),
     *                  W_tr(N_Half,N_r), X_tr(N_Half,N_r),
     *                  W_rr(N_Half,N_r), X_rr(N_Half,N_r)

        Call MtoM_rect (N_Half, N_Half, W_tt, X_tt, N_t, N_t)
        Call MtoM_rect (N_Half, N_Half, W_rt, X_rt, N_r, N_t)
        Call MtoM_rect (N_Half, N_Half, W_tr, X_tr, N_t, N_r)
        Call MtoM_rect (N_Half, N_Half, W_rr, X_rr, N_r, N_r)

        return
        end

c ======================================================= 3d

        Subroutine      MtoM_rect (Na, Nb, A, B, mx, my)
c B(mx,my) = A(mx,my) (size: mx x my)
        Integer         Na, Nb, mx, my,
     *                  i, j
        Complex*16      A(Na,my), B(Nb,my)

        if (mx.gt.Na .OR. mx.gt.Nb) Then !-------------------------+
           write (0,*) 'MtoM_rect: Na=',Na,'  Nb=',Nb,'  mx=',mx  !|
           stop 'MtoM_rect: oversize'                             !|
        endif !----------------------------------------------------+

        do      i=1,mx  !==========+
        do      j=1,my  !=======+  |
           B(i,j) = A(i,j)     !|  |
        enddo  !================+  |
        enddo  !===================+

        return
        end

c ======================================================= 4d

        Subroutine      MpmM_rect (Na, Nb, Nc, A, B, C, mx, my, mp)
c C = A +- B (size: mx x my)
        Integer         Na, Nb, Nc, mx, my, mp,
     *                  i, j
        Complex*16      A(Na,my), B(Nb,my), C(Nc,my)

        if (mx.gt.Na .OR. mx.gt.Nb .OR. mx.gt.Nc) Then !-----+
           write (0,*) 'MpmM_rect: Na=',Na,'  Nb=',Nb,      !|
     *                          '  Nc=',Nc,'  mx=',mx       !|
           stop 'MpmM_rect: oversize'                       !|
        endif !----------------------------------------------+

        do      i=1,mx !====================+
        do      j=1,my !=================+  |
           C(i,j) = A(i,j) + mp*B(i,j)  !|  |
        enddo  !=========================+  |
        enddo  !============================+

        return
        end
