c
c
c +==========================================+
c | These routines find x-ray wavefields for |
c | 1. gid_slm7.for (rectangular matrices)   |
c | (for the case where the reflection from a|
c |  multilayer is reduced to the recursive  |
c |  equations)                              |
c |------------------------------------------|
c |Contents:                                 |
c | 1. Subroutine:      Find_Fields_rect     |
c | 2. Subroutine:      Norm_Fields_rect     |
c | 3. Subroutine       Backup2_rect         |
c +==========================================+

c ####################################################### 1

          Subroutine Find_Fields_rect (T_0,  R_k1,      ! T_0,        R_{n+1}    - known
     *                                 T_k,  R_k,       ! T_n,        R_n        - to be calculated
     *                                 W_tt, W_tr,      ! W_n^tt,     W_n^tr     - known
     *                                 M_rt, M_rr,      ! M_{n+1}^rt, M_{n+1}^rr - known
     *                                 N_t, N_r, ifail)
c--------------------------------------------------------
c This routine implements Eq.(24) from S.A.Stepanov,
c E.A.Kondrashkina, R.Koehler, D.V.Novikov, G.Materlik,
c and S.M.Durbin, Phys. Rev. B, v.57, No 8, p.4829-4841, (1998).
c
c              rt    tr -1    rr          rt    tt
c   R = ( 1 - M   * W  )  * (M   * R   + M   * W  * T )
c    n         n+1   n        n+1   n+1   n+1   n    0
c
c        tt        tr
c   T = W  * T  + W  R
c    n   n    0    n  n
c
c Equations (24) must be progressively applied to all the layers
c starting at the substrate where R_N = 0.
c--------------------------------------------------------
        Integer         N_Roots, N_Half
        Parameter      (N_Roots = 4,
     *                  N_Half  = N_Roots/2)
        Complex*16      T_0(N_Half), R_k1(N_Half),
     *                  T_k(N_Half), R_k(N_Half),
     *                  W_tt(N_Half,N_Half),
     *                  W_tr(N_Half,N_Half),
     *                  M_rt(N_Half,N_Half),
     *                  M_rr(N_Half,N_Half),
     *                  T_1(N_Half,N_Half),        ! work matrix
     *                  T_2(N_Half,N_Half),        ! work matrix
     *                  S(N_Half)                  ! work vector
        Integer         N_t, N_r, ifail, i, j

        if (N_t.gt.N_Half .OR. N_r.gt.N_Half) Then !--------+
           write (0,*) 'Find_Fields_rect: N_t=',N_t        !|
           write (0,*) 'Find_Fields_rect: N_r=',N_r        !|
           write (0,*) 'Find_Fields_rect: N_half=',N_Half  !|
           write (0,*) 'Find_Fields_rect: size overflow'   !|
           ifail = -1                                      !|
           return                                          !|
        endif !---------------------------------------------+

        do    i=1,N_Half  !=========+
           S(i) = (0.0D0, 0.0D0)   !| make compiler happy
        enddo !=====================+

                                              !---------+
c T_1 = 1                                               |
        Call    M1 (N_Half,T_1,N_r)                    !|r x r
                                                       !|
c T_2 = M_rt*W_tr                                       |
        Call    MxM (N_Half, N_Half, N_Half,           !|r x r
     *               M_rt,   W_tr,   T_2,              !|
     *               N_r,    N_r,    N_t)              !|
                                                       !|
c T_1 = 1-M_rt*W_tr                                     |
        Call    MpmM (N_Half, N_Half, N_Half,          !|r x r
     *                T_1,    T_2,    T_1,             !|
     *                N_r,    -1)                      !|
                                                       !|
c T_1 = 1/(1-M_rt*W_tr)                                 |
        Call    Inverse_Matrix_simple (T_1,    T_1,    !|r x r
     *                                 N_Half, N_r,    !|
     *                                 ifail)          !|
        if (ifail .ne. 0) return                       !|
                                              !---------+
                                                      !----+
c T_k = W_tt*T_0 (temporary)                               |
        do    i=1,N_t  !==============================+    |
           T_k(i) = (0.0D0, 0.0D0)                   !|    |t-vector
           do  j=1,N_t  !===========================+ |    |
              T_k(i) = T_k(i) + W_tt(i,j) * T_0(j) !| |    |<-T_k(temporary)
           enddo  !=================================+ |    |
        enddo  !======================================+    |
                                                          !|
c S = M_rt*[W_tt*T_0] + M_rr*R_{k+1}                       |
        do    i=1,N_r  !==============================+    |
           S(i) = (0.0D0, 0.0D0)                     !|    |r-vector
           do  j=1,N_t  !===========================+ |    |
              S(i)   = S(i) + M_rt(i,j) * T_k(j)   !| |    |<-T_k(temporary)
           enddo  !=================================+ |    |
           do  j=1,N_r  !===========================+ |    |
              S(i)   = S(i) + M_rr(i,j) * R_k1(j)  !| |    |R_k1 is r-vector
           enddo  !=================================+ |    |
        enddo  !======================================+    |
                                                          !|
c R_k = [1/(1-M_rt*W_tr)] * {M_rt*[W_tt*T_0]+M_rr*R_{k+1}} |=T_1*S
        do    i=1,N_r  !==============================+    |
           R_k(i)   = (0.0D0, 0.0D0)                 !|    |
           do  j=1,N_r  !===========-===============+ |    |
              R_k(i) = R_k(i) + T_1(i,j) * S(j)    !| |    |r-vector
           enddo  !=================================+ |    |
        enddo  !======================================+    |
                                                          !|
c T_k = [W_tt * T_0] + W_tr * R_k                          |
        do    i=1,N_t  !==============================+    |
           do  j=1,N_r  !===========================+ |    |
              T_k(i) = T_k(i) + W_tr(i,j) * R_k(j) !| |    |t-vector
           enddo  !=================================+ |    |
        enddo  !======================================+    |
                                              !------------+
        Return
        End

c ####################################################### 2

          Subroutine  Norm_Fields_rect (T_k,  R_k,
     *                                  Ft_k, Fr_k,
     *                                  N_t, N_r)
        Integer         N_t, N_r, i
        Complex*16      T_k(N_t), Ft_k(N_t), 
     *                  R_k(N_r), Fr_k(N_r)

c+----------------------+
c|Renormalize wavefields|
c+======================+
c This is the renormalization of wavefields from the
c lower interfaces, as they are given by the recursive
c equations, to the upper interface, as they are required
c for the diffuse scattering calculations. In other words,
c this program calculates [D_k]*[F_k^(U)].
c If one needs the absolute values, he has to divide the
c result by [F_k^(U)] /see RFM52_R.FOR for an example/.

        do    i=1,N_t  !===============+
c The growing T-exponents were         |
c inverted in Process_Layer:           |
           T_k(i) = T_k(i) / Ft_k(i)  !|Ft_k(i)=exp(+im*k0*TERang*u(i)*t_k), Im(u(i)>0
        enddo  !=======================+

        do    i=1,N_t  !===============+
c The damping R-exponents were         |
c NOT inverted in Process_Layer:       |
           R_k(i) = R_k(i) * Fr_k(i)  !|Fr_k(i)=exp(-im*k0*TERang*u(i)*t_k), Im(u(i)<0
        enddo  !=======================+
        Return
        End

c ####################################################### 3

        Subroutine Backup2_rect (W_tr, W_tt,
     *                           X_tr, X_tt,
     *                           N_Half, N_t, N_r)
        Integer         N_Half, N_t, N_r,
     *                  i, j
        Complex*16      W_tr(N_Half,N_r),
     *                  W_tt(N_Half,N_t),
     *                  X_tr(N_Half,N_r),
     *                  X_tt(N_Half,N_t)

        if (N_t.gt.N_Half .OR. N_r.gt.N_Half) Then !------+
           write (0,*) 'Backup2_rect: N_t=',N_t          !|
           write (0,*) 'Backup2_rect: N_r=',N_r          !|
           write (0,*) 'Backup2_rect: N_half=',N_Half    !|
           stop 'Backup2_rect: size overflow'            !|
        endif !-------------------------------------------+

        do    i=1,N_t !===============+
        do    j=1,N_r !============+  |
           X_tr(i,j) = W_tr(i,j)  !|  |
        enddo   !==================+  |
        do    j=1,N_t !============+  |
           X_tt(i,j) = W_tt(i,j)  !|  |
        enddo   !==================+  |
        enddo   !=====================+

        return
        end
