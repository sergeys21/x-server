c
c                  rfms_urd7.for [NAG version]
c                  (replaces obsolete rfms_uni)
c
c                                           WE ARE HERE
c +========+     +========+         +=======================+
c |gid_slm7|-->--|rfmm7_rx|-->------|rfms_urd7 (GID_slm7)   |-->-+
c +--------+     +--------+         +-----------------------+    |
c +========+     +========+                                      |
c |gid_slm5|-->--|rfmm5_rx|->-+                                  |
c +--------+     +--------+   |                                  |
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
c =======================================================
c     These unified programs implement some assistance
c      operations for RFM's called from the following
c            REFLECTION and DIFFRACTION routines:
c
c         < ter_sl10, ter_sl12, ter_sl97, trds_97 >
c     < gid_sl32, gid_sl50, gid_sl55, gids_10, gids_97 >
c
c For the reflection from magnetic multilayers (mag_slxx)
c          see the module called RFMs_UMG.for
c (-- in the project stage!!! Meanwhile see MMRE99_R.for)
c
c       Both RMFs_URD and RFMs_UMG call RFMs_sss.for.
c
c In the pervious versions (before developing mag_slxx and
c  RFMs_UMG in May, 99), the RFMs_URD + RFMs_sss composed
c                one file called RFMs_Uni.for
c
c        The x-ray diffaction coefficient under GID
c        or  x-ray reflection coefficient under TER
c                       is computed
c         with the account for interface roughness.
c
c +==============================================================+ (old RFMs_Uni.for)
c | 1a. Subroutine: Process_Layer                                |
c | 2a. Subroutine: Smatrix                                      |
c | 3a. Subroutine: Dispersion_Roots                             |
c |--------------------------------------------------------------+ (RFMs_URD7 for GID_slm7)
c | 1b. Subroutine: Process_Layer_rect                           |
c | 2b. Subroutine: Smatrix_rect                                 |
c | 3b. Subroutine: Dispersion_Roots_rect                        |
c | 4b. Subroutine: Find_polynomial_roots -> Cpoly or C02ADF(NAG)|
c | 5b. Subroutine: Check_Im_u_count                             |
c +==============================================================+
c
c####################################################### 1a

        Subroutine      Process_Layer   (g0, gh, gi, Daa,
     *                                   x0, xh, k0, Depth,
     *                                   TERang, S, B, u, F,
     *                                   N_Fields,ifail)

        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields, N_Plus

        Complex*16      x0, xh(2),                 ! normalized polariz.
     *                  S(N_Dim,N_Fields),         ! layer S-matrix
     *                  B(N_Dim,N_Fields),         ! layer (1/S)-matrix
     *                  F(N_Fields),               ! layer F-matrix
     *                  u(N_Fields),               ! layer roots
     *                  q, qw,                     ! work cells
     *                  gh,                        ! normalized exit angle
     *                  im/(0.0D0,1.0D0)/

        Real*8          g0,                        ! normalized incid.angle
     *                  gi,                        ! normalized Psi
     *                  k0,                        ! k0 = 2*pi/wave
     *                  qr,qi,                     ! Re(qw),Im(qw)
     *                  qmax /100./                ! the upper limit for
                                                   ! exponent value. The
                                                   ! actual limit for PC
                                                   ! is about 707; for
                                                   ! other platforms it
                                                   ! may differ!

        Real*4          Depth, TERang, Daa

        Integer         ifail, i, j

        ifail = 0

c+--------------------------------------------+
                                             !|
c Compute u-roots and S-matrix of the layer:  |
        Call    Smatrix (g0, gh, gi,         !|
     *                   Daa, x0, xh,        !|
     *                   N_Fields,           !|
     *                   S, u, ifail)        !|
        if (ifail .ne. 0) Return             !|
c+--------------------------------------------+

c+-------------------------------------------------------+
c|Find inverse matrix B=1/S of matrix S:                 |
        if (N_Fields .le. 2) Then !-----------------+    |
           Call  Inverse_Matrix_simple (S, B,      !|    |
     *                                  N_Dim,     !|    |
     *                                  N_Fields,  !|    |
     *                                  ifail)     !|    |
        else  !-------------------------------------+    |
           Call  Inverse_Matrix        (S, B,      !|    |
     *                                  N_Fields,  !|    |
     *                                  ifail)     !|    |
        endif  !------------------------------------+    |
        if (ifail .ne. 0) Return                        !|
c+-------------------------------------------------------+

c+---------------------------------------------------------------------+
c|Compute F(t) matrix:                                                 |
        if (abs(Depth).gt.1.E-32) Then  !---------------------------+  |
           N_Plus = N_Fields / 2                                   !|  |
           q      = -im*k0*Depth*TERang                            !|  |
c The growing exponents are inverted                                |  |
           do      i=1,N_Plus    !===============================+  |  |
              qw =-q*u(i)                                       !|  |  |
              qr = dReal(qw)                                    !|  |  |
              qi = dImag(qw)                                    !|  |  |
c             if (qr .gt. +qmax) qw = dCmplx(+qmax,qi)          !|  |  |
              if (qr. gt. 1.D-64) Then !----------------------+  |  |  |must be negative
                 write(0,*) 'Process_Layer: wrong q-sign'    !|  |  |  |
                 write(0,*) 'i=',i,'  N_plus=',N_plus        !|  |  |  |
                 write(0,*) (u(j),j=1,N_Fields)              !|  |  |  |
              endif !-----------------------------------------+  |  |  |
              if (qr .lt. -qmax) qw = dCmplx(-qmax,qi)          !|  |  |
              F(i) = exp(qw)                                    !|  |  |
           enddo   !=============================================+  |  |
c The decreacing exponents aren't inverted                          |  |
           do      i=N_Plus+1,N_Fields !=========================+  |  |
              qw =+q*u(i)                                       !|  |  |
              qr = dReal(qw)                                    !|  |  |
              qi = dImag(qw)                                    !|  |  |
c             if (qr .gt. +qmax) qw = dCmplx(+qmax,qi)          !|  |  |
              if (qr. gt. 1.D-64) Then !----------------------+  |  |  |must be negative
                  write(0,*) 'Process_Layer: wrong q-sign'   !|  |  |  |
                  write(0,*) 'i=',i,'  N_plus=',N_plus       !|  |  |  |
                  write(0,*) (u(j),j=1,N_Fields)             !|  |  |  |
              endif !-----------------------------------------+  |  |  |
              if (qr .lt. -qmax) qw = dCmplx(-qmax,qi)          !|  |  |
              F(i) = exp(qw)                                    !|  |  |
           enddo   !=============================================+  |  |
        else  !-----------------------------------------------------+  |
           do      i=1,N_Fields  !====+                             |  |
              F(i) = (1.0D0, 0.0D0)  !|                             |  |
            enddo  !==================+                             |  |
        endif  !----------------------------------------------------+  |
c+---------------------------------------------------------------------+

        Return
        End

c####################################################### 2a

        Subroutine      Smatrix (g0, gh, gi,
     *                           Daa, x0, xh,
     *                           N_Fields,
     *                           S, u, ifail)
        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields,
     *                  ifail

        Complex*16      x0, xh(2), V(N_Dim),
     *                  S(N_Dim,N_Fields),         !S-matrix
     *                  u(N_Fields),               !dispersion roots
     *                  gh,                        !normalized exit angle (diffr.)
     *                  Sqrt_Eps                   !1 or sqrt(1+x0) (spec.)

        Real*8          g0,                        !normalized incid.angle
     *                  gi,                        !normalized Psi (diffraction!)
     *                  gi_k,                      !gi/(1+Daa)
     *                  xh2,                       !abs(xh*xh)
     *                  W(N_Dim)

        Real*4          Daa

        Integer         j

        Sqrt_Eps = (0.,0.)                      !make GNU compiler happy

        if (N_Fields .eq. 4) Then !------------------------+ DIFFR
                                                          !|
c+------------------------------------------------------+  |
c|February 22, 1995: that is from dh/h = -da/a          |  |
c|OPPOSITE TO (19) IN THE ARTICLE  !!!!!!!!             |  |
c|J.Phys.D: Appl.Phys. v.{27}, p.1923-1928, (1994)      |  |
c+----------------------------------------------------- |  |
c          gi_k = gi * (1.0D0 - Dble(Daa))             !|  |
c+------------------------------------------------------+  |
c|July 9, 1998: this modification of the above formula  |  |
c|follows the suggestion by Marius Grundmann (TU Berlin)|  |
c|Although the difference seems to be small, it provides|  |
c|quite considerable effect on strongly strained struc- |  |
c|tures like InAs/GaAs. This probably happens because   |  |
c|when gi_k >> xh, small changes in it provide a strong |  |
c|effect on the pattern (the small variations of gi_k   |  |
c|must be rather compared with xh,  than with gi itself!|  |
c+----------------------------------------------------- |  |
           xh2  = Abs(xh(1)*xh(2))                     !|  |
c Change of 2005/04/09 -- version M5. The Daa parameter |  |
c should only apply to crystalline layers. If we have   |  |
c approximated the layer as amorphous (because of small |  |
c xh or largre da/a), then it should be ignored:        |  |
           if (abs(xh2).gt.1.E-32) Then !--------+      |  |
              gi_k = gi / (1.0D0 + Dble(Daa))   !|      |  |
           else !--------------------------------+      |  |
              gi_k = gi                         !|      |  |
           endif !-------------------------------+      |  |
c+------------------------------------------------------+  |
                                                          !|
        else  !--------------------------------------------+ SPEC
                                                          !|
c For specular problem we use:                             |
c gh = 1          for sigma-polarization                   |
c gh = sqrt(1+x0)  for pi-polarization                     |
c This is a modification of 05/99 that allows to take      |
c into account the difference between sigma- and pi-       |
c polarizations for soft x-rays:                           |
           gi_k = 0.                                      !|
           Sqrt_Eps = gh                                  !|
c For compatibility with old programs that transfer gi=0:  |
           if (abs(Sqrt_Eps).lt.1.E-32) Sqrt_Eps=(1.,0.)  !|
                                                          !|
        endif !--------------------------------------------+

c Find dispersion equation roots:
        Call    Dispersion_Roots (g0, gh, gi_k, x0, xh,
     *                            u, V, W, N_Fields,
     *                            ifail)
        if (ifail .ne. 0) return

        if (N_Fields .eq. 2) Then !-------------+ SPEC
                                               !|
           do    j=1,N_Fields  !===========+    |
c cccccccc    S(1,j) = (1.0D0, 0.0D0)    !|     |1,1  -- before 05/99
c cccccccc    S(2,j) = u(j)              !|     |u,u  -- before 05/99
              S(1,j) = Sqrt_Eps          !|     | e , e       e=1          for sigma-
              S(2,j) = u(j)/Sqrt_Eps     !|     |u/e,u/e      e=sqrt(1+x0) for pi-
           enddo  !=======================+     |
                                               !|
        else  !---------------------------------+ GID
                                               !|
           do    j=1,N_Fields  !=============+  |  cr.      am.
              S(1,j) = W(j)                 !|  |1,1,1,1  1,0,1,0
              S(2,j) = V(j)                 !|  |V,V,V,V  0,1,0,1
              S(3,j) = W(j) * u(j)          !|  |u,u,u,u  u,0,u,0
c cc          S(4,j) = V(j) * (u(j) + gi_k) !|  |w,w,w,w  0,w,0,w
c cc Before 09-11-95 !!!!!!!! ------+        |  |
c cc After  09-11-95 !!!!!!!! ------+        |  |
c cc          S(4,j) = V(j) * (u(j) + gi)   !|  |w,w,w,w  0,w,0,w
c cc After  21-01-97 !!!!!!!! ------+        |  |
c cc (the previous version produced |        |  |
c cc   R>1 at large strains!!!)     |        |  |
              S(4,j) = V(j) * (u(j) + gi_k) !|  |w,w,w,w  0,w,0,w
           enddo  !==========================+  |
                                               !|
        endif  !--------------------------------+
        return
        end

c####################################################### 3a

        Subroutine      Dispersion_Roots (g0,gh,gi,x0,xh,
     *                                    u,V,W,N_Fields,
     *                                    ifail)
c------------------------------------------------------
c      Forms and solves the polynomial dispersion
c    equation of grazing-incidence x-ray diffraction.
c
c  On exit the roots are sorted in descending of Im(u).
c------------------------------------------------------
        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         N_Fields

        Complex*16 a(N_Dim+1),
     *             u(N_Fields),
     *             V(N_Fields),
     *             x0, xh(2),
     *             s0, sh,
     *             gh                              ! normalized exit angle

        Real*8     g0,                             ! normalized incid.angle
     *             gi,                             ! norm-zed Psi
     *             alpha,                          !a = (g0+gi)^2 - gh^2
     *             W(N_Fields),
     *             cfnr(4),
     *             ur(N_Dim+1),ui(N_Dim+1)

        Integer    ifail, i, j

        if (N_Fields .eq. 2) Then !----+ SPEC
           s0   = cdsqrt(g0*g0+x0)    !|
           u(1) = s0                  !|
           u(2) =-s0                  !|
           return                     !|
        endif !------------------------+

        alpha = (g0+gi)*(g0+gi) - Real(gh*gh)
        sh    = Sqrt(xh(1)*xh(2))

c Below we presume that |x0|~1:

c       if (Abs(alpha) .gt. (1.E+6)) Then !-----------+before 2005.04.24
        if (Abs(alpha) .gt. (1.E+8)) Then !-----------+after  2005.04.24
c Added 2003/08/06 -- check for big alpha             |
c (for this particular layer) and switch              |
c to "amorphous" model if |alpha| is very big:        |
c                                                     |
c 2005/04/09: this control is obsolete for gid_slM5.  |
c In gid_sl99 or gid_slm1 it was used to avoid nume-  |
c rical errors at big alpha. gid_slm5 dynamically     |
c changes layer Xh to 0 at the RFM level (see file    |
c RFMm5_RX.for).                                      |
           write (3,2) Abs(alpha),abs(sh),abs(x0)    !|
  2        format (
     *         ' --- Dispersion_Roots WARNING:',1x,
     *         'layer treated as amorphous at',1x,
     *         '|alpha/x0_max|=',g12.5,3x,
     *         '|xh/x0_max|=',g12.5,3x,
     *         '|x0/x0_max|=',g12.5)                 !|
           sh = (0.0D0, 0.0D0)                  !xh=0 |
        endif !---------------------------------  ----+
        if (Abs(sh) .gt. (1.E-6)) Then !--------------+GID (diffraction)
c Crystal layer:                                      |
c---------------                                      |
           s0   = g0*g0 + x0                         !|
           sh   = gh*gh + x0 - gi*gi                 !|
c Multipliers at ...                                  |
           a(1) = (1.0D0, 0.0D0)         !u**4        |
           a(2) = 2.*gi                  !u**3        |
           a(3) =-(s0 + sh)              !u**2        |
           a(4) =-2.*gi*s0               !u**1        |
           a(5) = s0*sh - xh(1)*xh(2)    !u**0        |x(+h)*x(-h)
                                                     !|
           cfnr(4) = abs(a(5))                       !|
c Additional normalization of roots. This may be      |
c required for high-order reflections with            |
c QB~90 degr and very small x0.                       |
           if (cfnr(4) .gt. 1.0D+8) Then !---+        |
              cfnr(2) = sqrt (cfnr(4))      !|        |
              cfnr(1) = sqrt (cfnr(2))      !|        |
              cfnr(3) = cfnr(2)*cfnr(1)     !|        |
              do i=1,4 !==================+  |        |
                 a(i+1) = a(i+1)/cfnr(i) !|  |        |
              enddo  !====================+  |        |
           else  !---------------------------+        |
              do i=1,4 !==================+  |        |
                 cfnr(i) = 1.0D0         !|  |        |
              enddo  !====================+  |        |
           endif  !--------------------------+        |
                                                     !|
c+--------------------------+                        !|
c|Find the dispersion roots:|                         |
c+--------------------------+                        !|
           Call Find_polynomial_roots(a,ur,ui,       !|
     *                               N_Fields,ifail) !|
           if (ifail.ne.0) return                    !|
                                                     !|
c+--------------------------+                         |
c|Sort the roots in the     |                         |
c|descending order of Im(u):|                         |
c+--------------------------+                         |
           Call  Sort_Roots (ur,ui,N_Fields)         !|
                                                     !|
c+--------------------+                               |
c|Check for the right |                               |
c|Number of Im's:     |                               |
c+--------------------+                               |
           Call Check_Im_u_count(ui,N_Fields,        !|
     *                             N_Fields/2,ifail) !|
           if (ifail.ne.0) return                    !|
                                                     !|
           do      i=1,N_Fields  !=============+      |
              u(i) = Dcmplx (ur(i),ui(i))     !|      |
     *             * cfnr(1)                  !|      |
              W(i) = 1.0D0                    !|      |
              V(i) = (u(i)*u(i) - s0) / xh(2) !|      |
           enddo  !============================+      |
                                                     !|
        else  !---------------------------------------+ AMORPH
                                                     !|
c Amorph. layer:                                      |
c---------------                                      |
           s0   = cdsqrt(g0*g0+x0)                   !|
           sh   = cdsqrt(gh*gh+x0)                   !|
           u(1) = s0                                 !|
           u(2) = sh - gi                            !|
           u(3) =-s0                                 !|
           u(4) =-sh - gi                            !|
                                                     !|
c+=======================================+            |
c|ATTENTION !!! "Amorphous" roots are not|            |
c|              sorted. That might cause |            |
c|              some discontinuities if  |            |
c|              one plots these roots or |            |
c|              the x-ray wavefields as a|            |
c|              function of layer number |            |
c+=======================================+            |
           do i=1,N_Fields  !=============+           |
              j    = i-2*(i/2)           !|           |
              W(i) = 1.0D0 * j           !|1,0,1,0    |
              V(i) = 1.0D0 - W(i)        !|0,1,0,1    |
           enddo  !=======================+           |
                                                     !|
        endif  !--------------------------------------+

        return
        end

c --------------------------------------------------------
c                      RFMs_urd7 (for gid_slm7)
c####################################################### 1b

        Subroutine      Process_Layer_rect (g0, gh, gi, Daa,
     *                                      x0, xh, k0, Depth,
     *                                      TERang, S, B, u, F,
     *                                      m_type,
     *                                      N_t, N_r,
     *                                      ifail)

        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         m_type,                    ! matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s
     *                  N_t, N_r,                  ! number of transmitted and reflected waves
     *                  N_Fields

        Complex*16      x0, xh(2),                 ! normalized polariz.
     *                  S(N_Dim,N_t+N_r),          ! layer S-matrix
     *                  B(N_Dim,N_t+N_r),          ! layer (1/S)-matrix
     *                  F(N_t+N_r),                ! layer F-matrix
     *                  u(N_t+N_r),                ! layer roots
     *                  q, qw,                     ! work cells
     *                  gh,                        ! normalized exit angle
     *                  im/(0.0D0,1.0D0)/

        Real*8          g0,                        ! normalized incid.angle
     *                  gi,                        ! norm-zed Psi
     *                  k0,                        ! k0 = 2*pi/wave
     *                  qr,qi,                     ! Re(qw),Im(qw)
     *                  qmax /100./                ! the upper limit for
                                                   ! exponent value. The
                                                   ! actual limit for PC
                                                   ! is about 707; for
                                                   ! other platforms it
                                                   ! may differ!

        Real*4          Depth, TERang, Daa

        Integer         ifail, i, j

        ifail = 0
        N_Fields = N_t + N_r

c+--------------------------------------------+
                                             !|
c Compute u-roots and S-matrix of the layer:  |
        Call Smatrix_rect (g0, gh, gi,       !|
     *                     Daa, x0, xh,      !|
     *                     m_type, N_t, N_r, !|
     *                     S, u, ifail)      !|
        if (ifail .ne. 0) Return             !|
c+--------------------------------------------+

c+-------------------------------------------------------+
c|Find inverse matrix B=1/S of matrix S:                 |
        if (N_Fields .le. 2) Then !-----------------+    |
           Call  Inverse_Matrix_simple (S, B,      !|    |simple inversion using determinant
     *                                  N_Dim,     !|    |
     *                                  N_Fields,  !|    |
     *                                  ifail)     !|    |
        else  !-------------------------------------+    |
           Call  Inverse_Matrix        (S, B,      !|    |inversion using Gassian solution
     *                                  N_Fields,  !|    |
     *                                  ifail)     !|    |
        endif  !------------------------------------+    |
        if (ifail .ne. 0) Return                        !|
c+-------------------------------------------------------+

c+----------------------------------------------------------------------+
c|Compute F(t) matrix:                                                  |
        if (abs(Depth).gt.1.E-32) Then  !----------------------------+  |
           q = -im*k0*Depth*TERang                                  !|  |
c The growing exponents are inverted                                 |  |
           do  i=1,N_t  !=========================================+  |  |
              qw =-q*u(i)                                        !|  |  |
              qr = dReal(qw)                                     !|  |  |
              qi = dImag(qw)                                     !|  |  |
c             if (qr .gt. +qmax) qw = dCmplx(+qmax,qi)           !|  |  |
              if (qr. gt. 1.D-64) Then !-----------------------+  |  |  |must be negative
                 write(0,*) 'Process_Layer_rect: wrong q-sign'!|  |  |  |
                 write(0,*) 'i=',i,'  N_t=',N_t,'  N_r=',N_r  !|  |  |  |
                 write(0,*) (u(j),j=1,N_Fields)               !|  |  |  |
              endif !------------------------------------------+  |  |  |
              if (qr .lt. -qmax) qw = dCmplx(-qmax,qi)           !|  |  |
              F(i) = exp(qw)                                     !|  |  |
           enddo   !==============================================+  |  |
c The decreacing exponents aren't inverted                           |  |
           do      i=N_t+1,N_Fields !=============================+  |  |
              qw =+q*u(i)                                        !|  |  |
              qr = dReal(qw)                                     !|  |  |
              qi = dImag(qw)                                     !|  |  |
c             if (qr .gt. +qmax) qw = dCmplx(+qmax,qi)           !|  |  |
              if (qr. gt. 1.D-64) Then !-----------------------+  |  |  |must be negative
                 write(0,*) 'Process_Layer_rect: wrong q-sign'!|  |  |  |
                 write(0,*) 'i=',i,'  N_t=',N_t,'  N_r=',N_r  !|  |  |  |
                 write(0,*) (u(j),j=1,N_Fields)               !|  |  |  |
              endif !------------------------------------------+  |  |  |
              if (qr .lt. -qmax) qw = dCmplx(-qmax,qi)           !|  |  |
              F(i) = exp(qw)                                     !|  |  |
           enddo   !==============================================+  |  |
        else  !------------------------------------------------------+  |
           do      i=1,N_Fields  !===+                               |  |
              F(i) = (1.0D0, 0.0D0) !|                               |  |
           enddo  !==================+                               |  |
        endif  !--------------------------------------------- -------+  |
c+----------------------------------------------------------------------+

        Return
        End

c####################################################### 2b

        Subroutine      Smatrix_rect (g0, gh, gi,
     *                                Daa, x0, xh,
     *                                m_type, N_t, N_r,
     *                                S, u, ifail)
        Integer         N_Dim
        Parameter      (N_Dim = 4)

        Integer         m_type,                    ! matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s
     *                  N_t, N_r,                  ! number of transmitted and reflected waves
     *                  N_Fields,
     *                  ifail

        Complex*16      x0, xh(2),
     *                  V(N_Dim),
     *                  S(N_Dim,N_t+N_r),          ! S-matrix
     *                  u(N_t+N_r),                ! dispersion roots
     *                  gh,                        ! normalized exit angle (diffr.)
     *                  Sqrt_Eps                   ! 1 or sqrt(1+x0) (spec.)

        Real*8          g0,                        ! normalized incid.angle
     *                  gi,                        ! normalized Psi (diffraction!)
     *                  gi_k,                      ! gi/(1+Daa)
     *                  xh2,                       ! abs(xh*xh)
     *                  W(N_Dim)

        Real*4          Daa

        Integer         j

        N_Fields = N_t + N_r
        Sqrt_Eps = (0.,0.)                         ! make GNU compiler happy

        if (m_type .lt. 4) Then !--------------------------+ DIFFR
                                                          !|
c+------------------------------------------------------+  |
c|February 22, 1995: that is from dh/h = -da/a          |  |
c|OPPOSITE TO (19) IN THE ARTICLE  !!!!!!!!             |  |
c|J.Phys.D: Appl.Phys. v.{27}, p.1923-1928, (1994)      |  |
c+----------------------------------------------------- |  |
c          gi_k = gi * (1.0D0 - Dble(Daa))             !|  |
c+------------------------------------------------------+  |
c|July 9, 1998: this modification of the above formula  |  |
c|follows the suggestion by Marius Grundmann (TU Berlin)|  |
c|Although the difference seems to be small, it provides|  |
c|quite considerable effect on strongly strained struc- |  |
c|tures like InAs/GaAs. This probably happens because   |  |
c|when gi_k >> xh, small changes in it provide a strong |  |
c|effect on the pattern (the small variations of gi_k   |  |
c|must be rather compared with xh,  than with gi itself!|  |
c+----------------------------------------------------- |  |
           xh2  = Abs(xh(1)*xh(2))                     !|  |
c Change of 2005/04/09 -- version M5. The Daa parameter |  |
c should only apply to crystalline layers. If we have   |  |
c approximated the layer as amorphous (because of small |  |
c xh or largre da/a), then it should be ignored:        |  |
           if (abs(xh2).gt.1.E-32) Then !--------+      |  |
              gi_k = gi / (1.0D0 + Dble(Daa))   !|      |  |
           else !--------------------------------+      |  |
              gi_k = gi                         !|      |  |
           endif !-------------------------------+      |  |
c+------------------------------------------------------+  |
                                                          !|
        else  !--------------------------------------------+ SPEC
                                                          !|
c For specular problem we use:                             |
c gh = 1          for sigma-polarization                   |
c gh = sqrt(1+x0)  for pi-polarization                     |
c This is a modification of 05/99 that allows to take      |
c into account the difference between sigma- and pi-       |
c polarizations for soft x-rays:                           |
          gi_k = 0.                                       !|
          Sqrt_Eps = gh                                   !|
c For compatibility with old programs that transfer gi=0:  |
          if (abs(Sqrt_Eps).lt.1.E-32) Sqrt_Eps=(1.,0.)   !|
                                                          !|
        endif !--------------------------------------------+

c Find dispersion equation roots:
        Call    Dispersion_Roots_rect(g0, gh, gi_k,
     *                                x0, xh,
     *                                u,                    !u(i)
     *                                V,                    !V(i)=(u(i)^2-g0^2-x0)/x_(-h)
     *                                W,                    !W(i)=1
     *                                m_type,
     *                                N_t, N_r,
     *                                ifail)
        if (ifail .ne. 0) return

        if     (m_type .eq. 0) Then !-----------------+ Diffraction GID, matrix type 0=4x4
                                                     !|
c Note: for amorphous layers:                         |
c u = (s0, sh-gi, -s0, -sh-gi)                        |
           do    j=1,N_Fields  !============+         |  cr.           am.
              S(1,j) = W(j)                !|D0       |1,1,1,1    1,  0, 1,  0
              S(2,j) = V(j)                !|Dh       |V,V,V,V    0,  1, 0,  1
              S(3,j) = W(j) * u(j)         !|d(D0)/dz |u,u,u,u    s0, 0,-s0, 0
              S(4,j) = V(j) *(u(j)+gi_k)   !|d(Dh)/dz |w,w,w,w    0, sh, 0,-sh
           enddo  !=========================+         |
                                                     !|
        elseif (m_type .eq. 1) Then !-----------------+ Asymmetric diffraction, matrix type 1=3x3i
                                                     !| (grazing incidence)
           do    j=1,N_Fields  !============+         | cr.      am.
              S(1,j) = W(j)                !|D0       |1,1,1    1, 1, 0
              S(2,j) = V(j)                !|Dh       |V,V,V    0, 0, 1
              S(3,j) = W(j) * u(j)         !|d(D0)/dz |u,u,u   s0,-s0,0
           enddo  !=========================+         |
                                                     !|
        elseif (m_type .eq. 2) Then !-----------------+ Asymmetric diffraction, matrix type 2=3x3h
                                                     !| (grazing exit
           do    j=1,N_Fields  !============+         |  cr.      am.
              S(1,j) = W(j)                !|D0       |1,1,1    1, 0, 0
              S(2,j) = V(j)                !|Dh       |V,V,V    0, 1, 1
              S(3,j) = V(j) *(u(j)+gi_k)   !|d(Dh)/dz |w,w,w    0,sh,-sh
           enddo  !=========================+         |
                                                     !|
        elseif (m_type .eq. 3) Then !-----------------+ Normal diffraction, matrix type 3=2x2d
                                                     !|
           do    j=1,N_Fields  !============+         |cr.    am.
              S(1,j) = W(j)                !|D0       |1,1    1,0
              S(2,j) = V(j)                !|Dh       |V,V    0,1
           enddo  !=========================+         |
                                                     !|
        elseif (m_type .eq. 4) Then !-----------------+ Specular reflection, matrix type 4=2x2s
                                                     !|
           do    j=1,N_Fields  !============+         |
              S(1,j) = Sqrt_Eps            !|         | e , e       e=1          for sigma-
              S(2,j) = u(j)/Sqrt_Eps       !|         |u/e,u/e      e=sqrt(1+x0) for pi-
           enddo  !=========================+         |
                                                     !|
        else !----------------------------------------+
                                                     !|
           write(0,*) 'Smatrix_rect: m_type=',m_type !|
           stop 'Smatrix_rect: wrong m_type'         !|
                                                     !|
        endif  !--------------------------------------+

        return
        end

c####################################################### 3b

        Subroutine      Dispersion_Roots_rect (g0, gh, gi,
     *                                         x0, xh,
     *                                         u, V, W,
     *                                         m_type,
     *                                         N_t, N_r,
     *                                         ifail)
c------------------------------------------------------
c      Forms and solves the polynomial dispersion
c    equation of grazing-incidence x-ray diffraction.
c
c  On exit the roots are sorted in descending of Im(u).
c------------------------------------------------------
        Integer         N_Dim
        Parameter       (N_Dim = 4)

        Integer         m_type,                    ! matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s
     *                  N_t, N_r,                  ! number of transmitted and reflected waves
     *                  N_Fields

        Complex*16      a(N_Dim+1),
     *                  u(N_t+N_r),
     *                  V(N_t+N_r),
     *                  x0, xh(2), xhh,
     *                  s0, sh, s02, sh2, sg2,     ! s02=s0^2, sh2=sh^2, sg2=sh^2-gi^2
     *                  gh,                        ! normalized exit angle
     *                  p, q

        Real*8          g0,                        ! normalized incid.angle
     *                  gi,                        ! gi=normalized Psi
     *                  alpha,                     ! a = (g0+gi)^2 - gh^2
     *                  W(N_t+N_r),
     *                  cfnr(4),
     *                  ur(N_Dim+1),ui(N_Dim+1)

        Integer         ifail, i, j

        N_Fields = N_t + N_r

        if (N_Fields .gt. N_Dim) Then !--------------------------------+
           write (0,*) 'Dispersion_Roots_rect: N_t=',N_t,' N_r=',N_r  !|
           write (0,*) 'Dispersion_Roots_rect: N_Fields=',N_Fields    !|
           stop 'Dispersion_Roots_rect: oversize'                     !|
        endif !--------------------------------------------------------+


        if (m_type .lt. 4) Then !------------------+ if diffraction
           alpha = (g0+gi)*(g0+gi) - Real(gh*gh)  !|
           xhh   = Sqrt(xh(1)*xh(2))              !|sqrt(x(+h)*x(-h))
        else !-------------------------------------+
           alpha = 0.0D0                          !|
           xhh   = (0.0D0, 0.0D0)                 !|
        endif !------------------------------------+

c Below we presume that |x0|~1:

        if (Abs(alpha) .gt. (1.E+8)) Then !-------------+
c Check for big alpha (for this particular layer) and   |
c switch to "amorphous" model if |alpha| is very big.   |
c This control is obsolete for gid_slM5 and newer.      |
c In gid_sl99 or gid_slm1 it was used to avoid numerical|
c errors at big alpha. gid_slm5 dynamically changes the |
c layer Xh to 0 at the RFM level (see RFMm5_RX.for).    |
           write (3,2) Abs(alpha),abs(xhh),abs(x0)     !|
  2        format (
     *           ' --- Dispersion_Roots WARNING:',1x,
     *           'layer treated as amorphous at',1x,
     *           '|alpha/x0_max|=',g12.5,3x,
     *           '|xh/x0_max|=',g12.5,3x,
     *           '|x0/x0_max|=',g12.5)                 !|
           xhh = (0.0D0, 0.0D0)                 !xh=0   |sqrt(x(+h)*x(-h))=0
        endif !-----------------------------------------+

        s02 = g0*g0 + x0
        s0  = cdsqrt(s02)
        if (m_type.lt.4) Then !----+
           sh2 = gh*gh + x0       !|
           sh  = cdsqrt(sh2)      !|
        else !---------------------+
           sh2 = (0.0, 0.0)       !|to keep compiler happy
           sh  = (0.0, 0.0)       !|to keep compiler happy
        endif !--------------------+

        if (Abs(xhh) .gt. (1.E-6)) Then !------------------------------+ Crystal layer
                                                                      !| =============
           if (m_type .eq. 0) Then !---------------------+GID matrx 4x4|
                                                        !|             |
c (u^2-g0^2-x0)((u+gi)^2-gh^2-x0) = xh(1)*xh(2)          |PRB, Eq.(6)  |
              sg2 = sh2 - gi*gi                         !|             |
c Multipliers at ...                                     |             |
c (u^2 - s02)*(u^2 + 2*gi*u -sg2) = xh(1)*xh(2)          |             |
c u^4 + 2*gi*u^3 - (s02+sg2)*u^2 - 2*gi*s02*u +          |             |
c                            + (s02*sg2-xh(1)*xh(2)) = 0 |             |
              a(1) = (1.0D0, 0.0D0)         !u**4        |             |
              a(2) = 2.*gi                  !u**3        |             |
              a(3) =-(s02 + sg2)            !u**2        |             |
              a(4) =-2.*gi*s02              !u**1        |             |
              a(5) = s02*sg2 - xh(1)*xh(2)  !u**0        |x(+h)*x(-h)  |
                                                        !|             |
              cfnr(4) = abs(a(5))                       !|             |
c Additional normalization of roots. This may be         |             |
c required for high-order reflections with               |             |
c QB~90 degr and very small x0.                          |             |
              if (cfnr(4) .gt. 1.0D+8) Then !---+        |             |
                 cfnr(2) = sqrt (cfnr(4))      !|        |             |
                 cfnr(1) = sqrt (cfnr(2))      !|        |             |
                 cfnr(3) = cfnr(2)*cfnr(1)     !|        |             |
                 do i=1,4 !==================+  |        |             |
                    a(i+1) = a(i+1)/cfnr(i) !|  |        |             |
                 enddo  !====================+  |        |             |
              else  !---------------------------+        |             |
                 do i=1,4 !==================+  |        |             |
                    cfnr(i) = 1.0D0         !|  |        |             |
                 enddo  !====================+  |        |             |
              endif  !--------------------------+        |             |
                                                        !|             |
c+--------------------------+                            |             |
c|Find the dispersion roots:|                            |             |
c+--------------------------+                            |             |
              Call Find_polynomial_roots(a,ur,ui,       !|             |
     *                                  N_Fields,ifail) !|             |
              if (ifail.ne.0) return                    !|             |
                                                        !|             |
              do      i=1,N_Fields  !======+             |             |
                 ur(i) = cfnr(1) * ur(i)  !|             |             |
                 ui(i) = cfnr(1) * ui(i)  !|             |             |
              enddo  !=====================+             |             |
                                                        !|             |
           elseif (m_type .eq. 1) Then !-----------------+extrim.assym.|
                                                        !|grazing incid|
c Multipliers at ...                                     |matrix 3x3i  |
c (u-s0)*(u+s0)*(u+gi-sh)*(u+gi+sh) = xh(1)*xh(2)        |             |
c From the Bragg condition:                              |             |
c gh^2=(g0+gi)^2-alpha, then gh ~ -(g0+gi) ~ -gi,        |             |
c Substitute y=u-s0, y must be small and then large root |             |
c (y+s0+gi-sh) ~ (y+g0+gi-gh) ~ (y+2*(g0+gi)) ~ 0 can be |             |
c discarded and for the rest we can divide by (s0+gi-sh):|             |
c (u^2-s0^2)*(u+gi+sh) = xh(1)*xh(2)/(s0+gi-sh)          |             |
c u^3 + (gi+sh)*u^2 - s0^2*u - s0^2*(sh+gi) - p = 0      |             |
              p    = (xh(1)*xh(2))/(s0+gi-sh)           !|             |
c             p    =-(xh(1)*xh(2))/(2.*sh)              !|             |
              a(1) = (1.0D0, 0.0D0)         !u**3        |             |
              a(2) = sh+gi                  !u**2        |             |
              a(3) =-s02                    !u**1        |             |
              a(4) =-s02*(sh+gi)-p          !u**0        |             |
              a(5) = 0                                  !|             |
                                                        !|             |
c+--------------------------+                            |             |
c|Find the dispersion roots:|                            |             |
c+--------------------------+                            |             |
              Call Find_polynomial_roots(a,ur,ui,       !|             |
     *                                  N_Fields,ifail) !|             |
              if (ifail.ne.0) return                    !|             |
                                                        !|             |
           elseif (m_type .eq. 2) Then !-----------------+extrem.assym.|
                                                        !|grazing exit |
c Multipliers at ...                                     |matrix 3x3h  |
c (u-s0)*(u+s0)*(u+gi-sh)*(u+gi+sh) = xh(1)*xh(2)        |             |
c From Bragg condition:                                  |             |
c gh^2=(g0+gi)^2-alpha, then gh~-(g0+gi)~-gi,            |             |
c Substitute y=u-s0, y must be small and then large root |             |
c (y+2*s0)~0 can be discarded and for the rest we can    |             |
c divide by (2*s0):                                      |             |
c (u-s0)*((u+gi)^2-sh^2) = xh(1)*xh(2)/(2*s0)            |             |
c Substitute: y=u+gi, then u=y-gi                        |             |
c (y-(s0+gi))*(y^2 - sh^2) = p                           |             |
c y^3 - (s0+gi)*y^2 - sh^2*y + (sh^2*(s0+gi)-p) = 0      |             |
              p    = (xh(1)*xh(2))/(2.*s0)              !|             |
              a(1) = (1.0D0, 0.0D0)         !y**3        |             |
              a(2) =-(s0+gi)                !y**2        |             |
              a(3) =-sh2                    !y**1        |             |
              a(4) = sh2*(s0+gi)-p          !y**0        |             |
              a(5) = 0                                  !|             |
                                                        !|             |
c+--------------------------+                            |             |
c|Find the dispersion roots:|                            |             |
c+--------------------------+                            |             |
              Call Find_polynomial_roots(a,ur,ui,       !|             |
     *                                  N_Fields,ifail) !|             |
              if (ifail.ne.0) return                    !|             |
                                                        !|             |
              do      i=1,N_Fields  !====+               |             |
                 ur(i) = ur(i) - gi     !|               |             |
              enddo  !===================+               |             |
                                                        !|             |
           elseif (m_type .eq. 3) Then !-----------------+normal diffr.|
                                                        !|matrix 2x2   |
c Multipliers at ...                                     |             |
c (u-s0)*(u+s0)*(u+gi-sh)*(u+gi+sh) = xh(1)*xh(2)        |             |
c u ~ s0 or (u+gi) ~ -sh and then                        |             |
c (u+s0) ~ 2*s0 and (u+gi-sh) ~ -2*sh                    |             |
c (u-s0)*(u+gi+sh) = -xh(1)*xh(2)/(4*s0*sh)              |             |
c Substitute u=s0+y :                                    |             |
c y*(y+s0+gi+sh) = p                                     |             |
c y^2 + (s0+gi+sh)*y = p                                 |             |
c y^2 + q*y + q^2/4 = p + q^2/4                          |             |
c (y+q/2)^2 = p + q^2/4                                  |             |
c y = -q/2 +- sqrt(p+q^2/4)                              |             |
              p    =-(xh(1)*xh(2))/(4*s0*sh)            !|             |
              q    = s0+gi+sh                           !|             |
              p = cdsqrt(p + 0.25*q*q)                  !|             |
              u(1) = s0 - 0.5*q + p                     !|             |
              u(2) = s0 - 0.5*q - p                     !|             |
              do      i=1,N_Fields  !====+               |             |
                 ur(i) = Dreal(u(i))    !|               |             |
                 ui(i) = dImag(u(i))    !|               |             |
              enddo  !===================+               |             |
                                                        !|             |
           else  !---------------------------------------+             |
                                                        !|             |
              write(0,*) ' Diffr. wrong m_type=',m_type !|             |
              stop 'Dispersion_Roots_rect wrong m_type' !|             |
                                                        !|             |
           endif !---------------------------------------+             |
                                                                      !|
c+--------------------------+                                          |
c|Sort the roots in the     |                                          |
c|descending order of Im(u):|                                          |
c+--------------------------+                                          |
           Call  Sort_Roots (ur,ui,N_Fields)                          !|
                                                                      !|
c+--------------------+                                                |
c|Check for the right |                                                |
c|Number of Im's:     |                                                |
c+--------------------+                                                |
           Call Check_Im_u_count (ui,N_Fields,N_t,ifail)              !|
           if (ifail.ne.0) return                                     !|
                                                                      !|
           do      i=1,N_Fields  !================+                    |
              u(i) = Dcmplx(ur(i),ui(i))         !|                    |
              W(i) = 1.0D0                       !|D0                  |
              V(i) = (u(i)*u(i) - s02) / xh(2)   !|Dh=V*D0 (PRB Eq.(7))|
           enddo  !===============================+                    |
                                                                      !|
        else !---------------------------------------------------------+ Amorphous layer
                                                                      !| ===============
c+==================================================+                  |
c|ATTENTION: "Amorphous" roots are not sorted, which|                  |
c|           might cause some discontinuities if one|                  |
c|           plots these roots or X-ray wavefields  |                  |
c|           as a function of layer number          |                  |
c+--------------------------------------------------+                  |
c|The order of roots u(i) determines what elements  |                  |
c|of W_rt in RFM8 will contain Es and Eh: reflect D0|                  |
c|(if exists) should go first and reflect Dh should |                  |
c|go second.                                        |                  |
c+==================================================+                  |
                                                                      !|
           if (m_type .eq. 0) Then !---------------------+GID Matrx 4x4|
                                                        !|             !
              u(1) = s0                 !transm  D0      |             |
              u(2) = sh - gi            !transm  Dh      |             |
              u(3) =-s0                 !reflect D0      |             |
              u(4) =-sh - gi            !reflect Dh      |             |
                                                        !|             |
              do i=1,N_Fields  !=============+           |             |
                 j    = i-2*(i/2)           !|           |             |
                 W(i) = 1.0D0 * j           !|1,0,1,0    |D0           |
                 V(i) = 1.0D0 - W(i)        !|0,1,0,1    |Dh=V*D0      |
              enddo  !=======================+           |             |
                                                        !|             |
           elseif (m_type .eq. 1) Then !-----------------+extrem.assym.|
                                                        !|grazing incid|
              u(1) = s0                 !transmit D0     |matrix 3x3i  |
              u(2) =-s0                 !reflect  D0     |             |
              u(3) =-sh - gi            !reflect  Dh     |             |
                                                        !|             |
              do i=1,N_Fields  !=============+           |             |
                 if (i.lt.3) then !--+       |           |             |
                    j = 1           !|       |           |             |
                 else !--------------+       |           |             |
                    j = 0           !|       |           |             |
                 endif !-------------+       |           |             |
                 W(i) = 1.0D0 * j           !|1,1,0      |D0           |
                 V(i) = 1.0D0 - W(i)        !|0,0,1      |Dh=V*D0      |
              enddo  !=======================+           |             |
                                                        !|             |
           elseif (m_type .eq. 2) Then !-----------------+extrem.assym.|
                                                        !|grazing exit |
              u(1) = s0                 !transmit D0     |matrix 3x3h  |
              u(2) = sh - gi            !transmit Dh     |             |
              u(3) =-sh - gi            !reflect  Dh     |             |
                                                        !|             |
              do i=1,N_Fields  !=============+           |             |
                 if (i.eq.1) then !--+       |           |             |
                    j = 1           !|       |           |             |
                 else !--------------+       |           |             |
                    j = 0           !|       |           |             |
                 endif !-------------+       |           |             |
                 W(i) = 1.0D0 * j           !|1,0,0      |D0           |
                 V(i) = 1.0D0 - W(i)        !|0,1,1      |Dh=V*D0      |
              enddo  !=======================+           |             |
                                                        !|             |
           elseif (m_type .eq. 3) Then !-----------------+normal diffr.|
                                                        !|matrix 2x2   |
              u(1) = s0                 !transm  D0      |             |
              u(2) = -sh - gi           !reflect Dh      |             |
                                                        !|             |
              do i=1,N_Fields  !==============+          |             |
                 j    = i-2*(i/2)            !|          |             |
                 W(i) = 1.0D0 * j            !|1,0       |D0           |
                 V(i) = 1.0D0 - W(i)         !|0,1       |Dh=V*D0      |
              enddo  !========================+          |             |
                                                        !|             |
           elseif (m_type .eq. 4) Then !-----------------+specular refl|
                                                        !|matrix 2x2   |
              u(1) = s0                 !transmit D0     |             |
              u(2) =-s0                 !reflect  D0    !|             !
                                                        !|             |
           else !----------------------------------------+             |
                                                        !|             |
              write(0,*) ' Amorph.wrong m_type=',m_type !|             |
              stop 'Dispersion_Roots_rect wrong m_type' !|             |
                                                        !|             |
           endif !---------------------------------------+             |
                                                                      !|
        endif  !-------------------------------------------------------+

        return
        end

c####################################################### 4b

        Subroutine Find_polynomial_roots(a,ur,ui,N_Fields,ifail)

        Integer         N_Dim
        Parameter       (N_Dim = 4)

        Integer         N_Fields,
     *                  ifail,
     *                  n, i

        Logical         fail,
     *                  use_NAG/.False./

        Complex*16      a(N_Fields+1)

        Real*8          ar(N_Dim+1), ur(N_Fields+1),
     *                  ai(N_Dim+1), ui(N_Fields+1),
     *                  tol

c Minimal value 1.+tol > 1.  :
        Data    tol     /2.220446049250313D-16/

        if (N_Fields .gt. N_Dim) Then !-----------------------------+
           write (0,*) 'Find_polynomial_roots: N_Fields=',N_Fields !|
           stop 'Find_polynomial_roots: oversize'                  !|
        endif !-----------------------------------------------------+

c+--------------------------+
c|Find the dispersion roots:|
c+--------------------------+
        n     = N_Fields + 1
        ifail = 0
        fail  = .False.
        do i=1,n  !==============+
           ar(i) = Real (a(i))  !|
           ai(i) = Imag (a(i))  !|
        enddo  !=================+

        if (use_NAG) Then  !-----------------------------------------+
c          Call C02ADF (ar,ai,n,ur,ui,tol,ifail)                    !|
        else !-------------------------------------------------------+
c Cpoly is not a part of Lapack. It is a general solver for the roots|
c of a polynomial, highly recommended by Netlib (Lapack distributor).|
c See algoritm 419 at https://netlib.org/toms/                       |
c n is the number of coefficients of the reduced polynomial after    |
c C02ADF. It is not used in Cpoly and therefore needs to be set to 1.|
           n = 1                                                    !|
           Call Cpoly (ar,ai,N_Fields,ur,ui,fail)                   !|
        endif !------------------------------------------------------+

        if (ifail.ne.0 .OR. n.gt.1 .OR. fail) Then !---+
c Possible errors from NAG program C02ADF:             |
c IFAIL=1:                                             |
c EITHER AR(1)=0 AND AC(1)=0, OR N.LT.2, OR N.GT.100.  |
c IFAIL=2:                                             |
c A POSSIBLE SADDLE POINT HAS BEEN DETECTED BY THE     |
c ROUTINE. (AT SUCH A POINT, BOTH THE FIRST AND        |
c SECOND DERIVATIVES OF THE POLYNOMIAL ARE ZERO).      |
                                                      !|
c Write DEBUG info:                                    |
c (these statements may not be reached --              |
c -- see error handling settings in P01ABZ)            |
           write (3,1) (N_Fields-n+1), N_Fields,      !|
     *                 ifail,fail                     !|
           do i=n,N_Fields !==========+                |
              write (3,2) 'U: ',     !|                |
     *                   ur(i),      !|                |
     *                   ui(i)       !|                |
           enddo !====================+                |
           do i=1,N_Fields+1  !=======+                |
              write (3,2) 'A: ',     !|                |
     *                   Real(a(i)), !|                |
     *                   Imag(a(i))  !|                |
           enddo !====================+                |
        endif !----------------------------------------+
        if (n .gt. 1) ifail = ifail+10

  1     format(' *** Find_polynomial_roots ERROR: failed',
     *         ' to find polynomial equation roots.'/
     *         ' Total roots found=',i1,' of max=',i1,
     *         ' ifail=',i3,' fail=',l)
  2     format(1x,a,g14.7,'+ i * ',g14.7)
        return
        end

c####################################################### 5b

        Subroutine Check_Im_u_count(ui,nmax,npos,ifail)

        Integer         nmax,npos,np,ifail,i

        Real*8          ui(nmax)

c+--------------------+
c|Check for the right |
c|Number of Im's:     |
c+--------------------+
        ifail = 0
        np = 0
        do i=1,nmax  !====================+
           if (ui(i) .gt. 0.) np = np+1  !|
        enddo  !==========================+
        if (np .ne. npos) Then !--------+
           write (3,1) np, npos, nmax  !|
           ifail = 100                 !|
        endif  !------------------------+
  1     format (' *** Dispersion_Roots ERROR: Count of Im(u)=',i1,
     *          ' expected=',i1,' among total=',i1)
        return
        end
