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
c                                         gid_sl97,
c                                         gids_10
c                                         ter_sl10
c                                         ter_sl12
c                                         trds_10
c
c *** Version 4.1::4.3 ***  Interface roughness included,
c                           Exponent index = 2
c                           Exponent corrected
c                           Specially for asymmetric cases
c                           Normal lattice strain added
c ***   Version 4.4   ***   Extended for general case (GID, XEAD,...)
c ***   Version 5.0   ***   Simplified in order to calculate the
c                           x-ray wavefields needed for diffuse
c                           scattering and secondary emission:
c                           large-scale roughness and transition
c                           layers excluded
c ***   Version 5.5   ***   Simplified: SL excluded (for Petrashen)
c ***   Version M5    ***   Added dynamic switching to amorphous
c                           layer and back at big alpha.
c ***   Version M7    ***   Added manually specified phases of
c                           structural amplitudes and dynamic choice
c                           of the size of scattering matrix
c                           depending on the case: GID/EAD/normal.
c
c                WE ARE HERE
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
c +====================================================================+
c |This  routine computes GID reflection coefficients                  |
c |Rh, R0 from multilayered structure at 1 angular point               |
c |   ( Incident angle = f0,     Exit angle = fh ).                    |
c |Input data are transferred to RFM in common/RFM7dat/                |
c |--------------------------------------------------------------------|
c |Contents:                                                           |
c | 1. Real*4 function RFM_m7                                          |
c | 2. Subroutine      Algorithm_Type                                  |
c | 3. Subroutine      DetectLayerTypes                                |
c |--------------------------------------------------------------------|
c |Subroutine calls:                                                   |
c | - rfms_urd7.for (Process_Layer_rect, Smatrix_rect,                 |
c |                  Dispersion_Roots_rect)                            |
c | - rfms_uuu.for  (Process_Interface,                                |
c |                 Inverse_Matrix, Gauss_Solution, Copy_Layer,        |
c |                 Process_Interface2, Inverse_Matrix_simple, Backup4,|
c |                 MxM, M0, M1, MtoM, MpmM, Element_Max, Sort_Roots,  |
c |                 Process_Interface2_rect, Backup4_rect, MtoM_rect,  |
c |                 MpmM_rect                                          |
c +====================================================================+

c ======================================================= 1
        Real*8  Function        RFM8_m7 (alpha)
c  ** Reflex_from_Multi **
c--------------------------------------------------------
c+===============================+
        Include 'gid_slm7.inc'  !|
c+===============================+
c The above INC file contains smth like:
c ccc   Parameter       (N_Top_Max     = 1001,
c ccc*                   N_Total_Max   = N_Top_Max+1)

        Integer          N_Dim,
     *                   N_Half,
     *                   N_t, N_r, N_Fields
        Parameter       (N_Dim  = 4,
     *                   N_Half   = N_Dim/2)
c       Parameter       (N_Fields = 4)

        Complex*16 D0(N_Dim,N_Total_Max+1),             !D0-fields in layers
     *             Dh(N_Dim,N_Total_Max+1),             !Dh-fields in layers
     *             u (N_Dim,N_Total_Max+1),             !dispersion eqn.roots
     *             F (N_Dim,N_Total_Max+1),             !F-matrix in layers
     *             X (N_Dim,N_Dim),                     !X-matrix of interface
     *             S(N_Dim,N_Dim),                      !S_matrix 4x4
     *             B(N_Dim,N_Dim),
     *             B1(N_Dim,N_Dim),                     !1/S_matrix
     *             M_rr(N_Half,N_Half,N_Total_Max),     !product of operators
     *             M_rt(N_Half,N_Half,N_Total_Max),     !product of operators
     *             M_tr(N_Half,N_Half),                 !product of operators
     *             M_tt(N_Half,N_Half),                 !product of operators
     *             W_rr(N_Half,N_Half),                 !product of operators
     *             W_rt(N_Half,N_Half),                 !product of operators
     *             W_tr(N_Half,N_Half,N_Total_Max),     !product of operators
     *             W_tt(N_Half,N_Half,N_Total_Max),     !product of operators
     *             x0n,xhn(2),                          !normalized x0,xh,x(-h),
     *             Es, Eh,                              !vacuum wave amplitudes
     *             V,                                   !work cell (Dh=V*D0)
     *             s02,                                 !work cell (g0^2+x0)
     *             gh,                                  !normalized exit angle
     *             ikuz,                                !i*k*u*z
     *             i2pp,                                !i*2*pi*standing_phase
     *             SUMM,                                !summ of waves
     *             Im /(0.,1.)/                         !imaginary 1.

        Complex*8  xe(N_Total_Max+1),                   !xh_effective
     *             xs_phase,                            !i*xs*pi
     *             ii /(0.,1.)/                         !imaginary "1"

        Real*8     g0, gh_r, gh2,                       !normalized angles
     *             alpha,                               !Bragg deviation
     *             gi,                                  !normalized Psi (in TERang units)
     *             gi_k,                                !gi/(1+Daa) gi in layer
     *             k0,                                  !k0 = 2*pi/wave
     *             rkz, ikz,                            !Real(i*k*u*z), Imag(i*k*u*z)
     *             hzz,                                 ! h_z * z
     *             z,                                   !depth inside layers
     *             zi,                                  !same, relative to a closest interface
     *             eps /0.001/,                         !how much Rh & R0 may exceed 1 due to rounding
c    *             small_R8 /2.23D-308/                 !smallest possible Real*8 number in Fortran
     *             small_R8 /2.23D-208/                 !small number to prevent NaN in calculating R0s,Rhs

        Real*4     pi,                                  !3.1416
     *             Xrh,Xih,Xir,Xii,                     !normalized parts of xh
     *             Sgm_i                                !Sigma of interface

        Integer    Kod,                                 !Matrix(0)/Recursive(1) flag
     *             n_cryst,                             !Number of crystalline layers after approximation
     *             i, j, j2, k, kk

c----------------------
        Complex*8  x0(N_Total_Max+1),
     *             xh(N_Total_Max+1)

        Real*4     xf(N_Total_Max+1),           !phase shift between xrh & xih (non-cubic)
     *             xs(N_Total_Max+1),           !M7: xh phase shift in layer (Kaganer, monolayers)
     *             Thickness(N_Total_Max),      !thicknesses of layers
     *             Sigma(N_Total_Max),          !rms rougness of upper interface
     *             Daa(N_Total_Max+1),          !normal lattice strain
     *             Wave,                        !X-ray wavelength
     *             xabs,                        !abs(x0_max)
     *             TERang,                      !critical TER angle
     *             xh_weak_limit,               !M5: the weakest possible xh
     *             alpha_max,                   !M5: max(alpha/xh) to treat crystal as "crystal"
     *             standing_range(2),           !M7: offsets in A from reference interface
     *             standing_step,               !M7: SW depth step in A
     *             standing_phase               !M7: SW phase in units of pi

        Integer    N_Top,                       !top layer sublayers
     *             N_Total,                     !total of layers
     *             i_standing,                  !M7: 1=yes/0=no SW
     *             i_reference,                 !M7: reference interface for SW (0=surface)
     *             n_standing,                  !M7: Number of SW depth points
     *             m_reduction,                 !M7: matrix size reduction flag (0/1/2)
     *             m_type,                      !M7: matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s
     *             iDebug                       !Debug flag
c This is common between all 3 programs:
c gid_sl, gid_in and RFM:
        Common /GEN_M7/ x0, xh, xf, xs,
     *                  Thickness, Wave,
     *                  Sigma, Daa,
     *                  xabs, TERang,
     *                  xh_weak_limit,          !M5
     *                  alpha_max,              !M5
     *                  standing_range,         !M7: new
     *                  standing_step,          !M7: new
     *                  standing_phase,         !M7: new
     *                  i_standing,             !M7: new
     *                  i_reference,            !M7: new
     *                  n_standing,             !M7: new
     *                  m_reduction,            !M7: new
     *                  m_type,                 !M7: new
     *                  N_Top, N_Total,
     *                  iDebug

c----------------------
        Real*8     Positions(N_Total_Max),      !depth of interfaces
     *             Phase(N_Total_Max),          !int(d(h_z*dz))
     *             SW(N_Stand_Max),             !M7: Standing waves intensity
     *             f0,                          !incident angle
     *             Scan_angle,                  !scan angle
     *             R0                           !reflection coefficient

        Real*4     Psi,                         !misorientation angle
     *             Pol(2)                       !polarization factors

        Logical    Account_Sigma

        Integer    N_Used,                      !number of used layers
     *             jPol,                        !polariz.index
     *             ifail                        !failure code
c This is common between
c gid_sl and RFM only:
        Common /RFM7dat/ R0, Positions, Phase,
     *                   SW,                    !M7: new
     *                   f0, Scan_angle,
     *                   Psi, Pol,
     *                   Account_Sigma,
     *                   N_Used, jPol, ifail

c----------------------
        i = N_Top_Max                           !to suppress unused paramters warnings
        i = N_Total_Max                         !to suppress unused paramters warnings
        i = N_Pts_Max                           !to suppress unused paramters warnings
        i = N_Stand_Max                         !to suppress unused paramters warnings
c----------------------
        if     (m_type .eq. 0) Then !----+      !matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s
           N_t = 2                      !|
           N_r = 2                      !|
        elseif (m_type .eq. 1) Then !----+      !3x3i grazing incidence
           N_t = 1                      !|      !  (incident)
           N_r = 2                      !|      !  (reflected+diffracted)
        elseif (m_type .eq. 2) Then !----+      !3x3h grazing exit
           N_t = 2                      !|      !  (incident+diffracted)
           N_r = 1                      !|      !  (diffracted&reflected)
        elseif (m_type .eq. 3) Then !----+      !2x2d normal diffraction
           N_t = 1                      !|      !  (incident)
           N_r = 1                      !|      !  (diffracted)
        elseif (m_type .eq. 4) Then !----+      !2x2s specular reflection
           N_t = 1                      !|      !  (incident)
           N_r = 1                      !|      !  (reflected)
        else !---------------------------+
           write (0,*) 'm_type=',m_type !|
           stop 'RFM8_m7: wrong m_type' !|
        endif !--------------------------+
        N_fields = N_t + N_r
        Kod = 1                                 !this program uses recursive method only
        Sgm_i = 0
c--------------------------------------------------------
c       xabs   = Abs(x0_max)                    !defined in main
c       TERang = Sqrt(xabs)                     !defined in main
        pi     = 4.* Atan(1.)
        k0     = 2.* pi/wave
        g0     = f0 / TERang
        gi     = Psi/ TERang
        gh2    = (g0+gi)**2 - alpha/xabs
        if (gh2 .ge. 0.) Then !---------------+
           gh  = Dsqrt(gh2)                  !|
        else   !------------------------------+
           gh  = Dcmplx(0.0D0,Dsqrt(-gh2))   !|
        endif  !------------------------------+

        R0      = 0.0D0
        RFM8_m7 = 0.0D0
        ifail   = 0

        if (i_standing .ne. 0) Then !-------+
           do i=1,N_Fields  !============+  |
           do k=1,N_Total+1 !==========+ |  |
              D0(i,k) = (0.D0, 0.D0)  !| |  |
              Dh(i,k) = (0.D0, 0.D0)  !| |  |
              u (i,k) = (0.D0, 0.D0)  !| |  |
              F (i,k) = (1.D0, 0.D0)  !| |  |
           enddo  !====================+ |  |
           enddo  !======================+  |
        endif  !----------------------------+

c 2005/04/09: added for version gid_slm5.
c Detect what layers must be approximated
c as amorphous because of large alpha.
c The most of paramters are transferred
c via common /GEN_M7/. This subroutine
c fills the array xe:
        Call DetectLayerTypes(Scan_angle,Pol(jPol),g0,gi,gh,xe,jPol,
     *                                                  n_cryst,ifail)
        if (ifail .ne. 0)   goto 100
        if (n_cryst .eq. 0) goto 100

c For zero angles g0=0 or gh=0 "vacuum" scattering matrix
c cannot be inverted. But, the wave amplitudes in vacuum
c are known to be zero in this case and, therefore, we can
c return from rfm with zero amplitudes.
c-----------
c In case where the fields inside crystal must be computed
c this drawback can be overcome, if we take into account
c that wave amplitude E_0h of incoming diffracted wave
c is equal to zero and, therefore, the respective coeffi-
c cient Sm(..,..)=gh=0 in "vacuum" matrix does not play any
c role and can be changed from zero value to make the
c matrix invertable.

        if (       (g0 .lt. 1.0D-10) .OR.
     *      (CDabs(gh) .lt. 1.0D-10) )  Return

c
c                     Process Vacuum
c
c+-----------------------------------------------------------+
        N_Used = 0                                          !|
        x0n    = x0(N_Used+1) / xabs                        !|
        Xrh    = Real (xe(N_Used+1)) / xabs                 !|
        Xih    = Aimag(xe(N_Used+1)) / xabs                 !|
        Xir    = Xih * sin(xf(N_Used+1)*pi)                 !|
        Xii    = Xih * cos(xf(N_Used+1)*pi)                 !|
        xhn(1) = Cmplx(Xrh+Xir,Xii)                !x(+h)   !|
        xhn(2) = Cmplx(Xrh-Xir,Xii)                !x(-h)   !|
c Extension introduced in M7 by suggestion of V.Kaganer:     |
c lattice of a layer may be sfifted because it may be a      |
c stack of ABAB... or BABA. This is accounted by the phase   |
c factor xs*pi that can be ranged between -2*pi and 2*pi.    |
        if (abs(xs(N_Used+1)) .gt. 1.E-32) Then !--+non-zero |
           xs_phase = ii*xs(N_Used+1)*pi          !|         |
           xhn(1) = xhn(1) * exp(+xs_phase)       !|x(+h)    |
           xhn(2) = xhn(2) * exp(-xs_phase)       !|x(-h)    |
        endif !------------------------------------+         |
                                                            !|
        Call Process_Layer_rect(g0, gh, gi,                 !|see rfms_urd7.for
     *                          Daa(N_Used+1),              !|
     *                          x0n, xhn, k0,               !|
     *                          0.,             !Thickness(0)|
     *                          TERang, S, B,               !|
     *                          u(1,N_Used+1),              !|
     *                          F(1,N_Used+1),              !|
     *                          m_type, N_t, N_r,           !|
     *                          ifail)                      !|
        if (ifail .ne. 0) goto 100                          !|
c+-----------------------------------------------------------+

c
c                     Process Top Layer
c
        Sgm_i = 0.                              !initialize when unused
        do    N_Used=1,N_Total  !=============================+
                                                             !|
           Call MtoM (N_Dim, N_Dim, B, B1, N_Fields)         !|
                                                             !|
           x0n    = x0(N_Used+1) / xabs                      !|
           Xrh    = Real (xe(N_Used+1)) / xabs               !|
           Xih    = Aimag(xe(N_Used+1)) / xabs               !|
           Xir    = Xih * sin(xf(N_Used+1)*pi)               !|
           Xii    = Xih * cos(xf(N_Used+1)*pi)               !|
           xhn(1) = Cmplx(Xrh+Xir,Xii)              !x(+h)   !|
           xhn(2) = Cmplx(Xrh-Xir,Xii)              !x(-h)   !|
c Extension introduced in M7 by suggestion of V.Kaganer:      |
c lattice of a layer may be sfifted because it may be a       |
c stack of ABAB... or BABA. This is accounted by the phase    |
c factor xs*pi that can be ranged between -2*pi and 2*pi.     |
           if (abs(xs(N_Used+1)) .gt. 1.E-32) Then !---+      |
              xs_phase = ii*xs(N_Used+1)*pi           !|      |
              xhn(1) = xhn(1) * exp(+xs_phase)        !|x(+h) |
              xhn(2) = xhn(2) * exp(-xs_phase)        !|x(-h) |
           endif !-------------------------------------+      |
                                                             !|
           Call Process_Layer_rect(g0, gh, gi,               !|see rfms_urd7.for
     *                             Daa(N_Used+1),            !|
     *                             x0n, xhn, k0,             !|
     *                             Thickness(N_Used),        !|
     *                             TERang, S, B,             !|
     *                             u(1,N_Used+1),            !|
     *                             F(1,N_Used+1),            !|
     *                             m_type, N_t, N_r,         !|
     *                             ifail)                    !|
           if (ifail .ne. 0)   goto 100  !--------------------+----+
c Compute X-matrix of layer upper interface and optionally    |    v
c average it over interface roughness:                        |
           if (Account_Sigma) Sgm_i = Sigma(N_Used)          !|
           Call Process_Interface_no_x (k0, TERang, Sgm_i,   !|
     *                                  B1, S,               !|
     *                                  u(1,N_Used),         !|
     *                                  u(1,N_Used+1),       !|
     *                                  X, N_Fields)         !|
c                                                             |
c Split X-matrix into four M-blocks and then calculate        |
c four W-matrices                                             |
           Call Process_Interface2_rect (X, F(1,N_Used+1),   !|
     *                                   M_rr(1,1,N_Used),   !|
     *                                   M_rt(1,1,N_Used),   !|
     *                                   M_tr,               !|
     *                                   M_tt,               !|
     *                                   W_rr,               !|
     *                                   W_rt,               !|
     *                                   W_tr(1,1,N_Used),   !|
     *                                   W_tt(1,1,N_Used),   !|
     *                                   N_t, N_r,           !|
     *                                   N_Used,             !|
     *                                   ifail)              !|
           if (ifail .ne. 0) goto 100  !----------------------+----+
           if (N_Used.lt.N_Total) Then  !--------------+      |    v
              Call Backup2_rect (W_tr(1,1,N_Used),    !|      |
     *                           W_tt(1,1,N_Used),    !|      |
     *                           W_tr(1,1,N_Used+1),  !|      |
     *                           W_tt(1,1,N_Used+1),  !|      |
     *                           N_Half, N_t, N_r)    !|      |
           endif  !------------------------------------+      |
        enddo  !==============================================+
c We may exit from the above loop with N_Used=N_Total+1, so correct:
        N_Used = N_Total
c
c                   Find solutions Es, Eh
c
c What elements of W_rt contain Es and Eh depends on the order of
c reflect D0 and reflect Dh in Dispersion_Roots_rect. When both
c reflect D0 and reflect Dh exits, then if D0 goes first and Dh goes
c second, Es and Eh are given by W_rt(1,1) and W_rt(1,2). Otherwise,
c they are inverted.
        if     (m_type .eq. 0) Then !----+      !matrix type 0=4x4,1=3x3i,2=3x3h,3=2x2d,4=2x2s
           Es = W_rt(1,1)               !|      |  (spec. reflected)
           Eh = W_rt(2,1)               !|      |  (reflected+diffracted)
        elseif (m_type .eq. 1) Then !----+      !3x3 grazing incidence
           Es = W_rt(1,1)               !|      !  (spec. reflected)
           Eh = W_rt(2,1)               !|      !  (diffracted)
        elseif (m_type .eq. 2) Then !----+      !3x3 grazing exit
           Es = 0.                      !|      !  (spec. reflected )
           Eh = W_rt(1,1)               !|      !  (diffracted)
        elseif (m_type .eq. 3) Then !----+      !2x2 normal diffraction
           Es = 0.                      !|      !  (spec.reflected)
           Eh = W_rt(1,1)               !|      !  (diffracted)
        elseif (m_type .eq. 4) Then !----+      !2x2 specular reflection
           Es = W_rt(1,1)               !|
           Eh = 0.                      !|
        else !---------------------------+
           write (0,*) 'm_type=',m_type !|
           stop 'RFM8_m7: wrong m_type' !|
        endif !--------------------------+

c Compute the reflection coefficients:

        R0  = Abs(Es)**2
        if (g0 .gt. 0.) then !------------------------+
           gh_r = Dreal(gh)                          !|
           if (gh_r .gt. 0.0) Then !---------------+  |
              RFM8_m7 = (gh_r/g0) * (Abs(Eh)**2)  !|  |
           else !----------------------------------+  |
              RFM8_m7 = 0.0                       !|  |
           endif !---------------------------------+  |
        endif  !--------------------------------------+
        if (RFM8_m7 .lt. 0.0) then  !-------------+
           write (3,*) 'RFM8_m7: Rh=', RFM8_m7,  !|
     *                 ' at scan=',Scan_angle    !|
           ifail = 990                           !|
           goto 100  !----------------------------+-------------+
        elseif (RFM8_m7 .lt. small_R8) Then !-----+            !v exclude NaN
          RFM8_m7 = 0.0D0                        !|
        endif  !----------------------------------+
        if (RFM8_m7 .gt. 1.0) then  !-------------+
           write (3,*) 'RFM8_m7: Rh=', RFM8_m7,  !|
     *                 ' at scan=',Scan_angle    !|
           if (RFM8_m7 .gt. 1.0+eps) then !---+   |
              ifail = 991                    !|   |
              goto 100  !---------------------+---+-------------+
           else !-----------------------------+   |             v
              RFM8_m7 = 1.0                  !|   |
           endif !----------------------------+   |
        endif  !----------------------------------+
        if (R0 .lt. 0.0) then !-------------------+
           write (3,*) 'RFM8_m7: R0=', R0,       !|
     *                 ' at scan=',Scan_angle    !|
           ifail = 880                           !|
           goto 100  !----------------------------+-------------+
        elseif (R0 .lt. small_R8) Then !----------+            !v exclude NaN
          R0 = 0.0D0                             !|
        endif  !----------------------------------+
        if (R0 .gt. 1.0) then !-------------------+
           write (3,*) 'RFM8_m7: R0=', R0,       !|
     *                 ' at scan=',Scan_angle    !|
           if (R0 .gt. 1.0+eps) then !--------+   |
              ifail = 881                    !|   |
              goto 100  !---------------------+---+-------------+
           else !-----------------------------+   |             v
              R0 = 1.0                       !|   |
           endif !----------------------------+   |
        endif  !----------------------------------+
        if (i_standing .eq. 0) Return

c
c Compute the wavefields
c
        if     (m_type .eq. 0) Then !--------------+ Diffraction GID, matrix type 0=4x4
                                                  !|
           D0(  1,  1) = (1.0D0, 0.0D0)           !| D0
           D0(  2,  1) = (0.0D0, 0.0D0)           !|
           D0(N_t+1,1) = Es                       !|
           D0(N_t+2,1) = Eh                       !|
           D0(  1,  N_Used+1) = W_tt(1,1,N_Used)  !|
c cc *                        * D0(1,1)           !|1
c cc /                        / F(1,N_Used+1)     !|1
c cc +                        + W_tt(1,2,N_Used)  !|
c cc *                        * D0(2,1)           !|0
c cc /                        / F(2,N_Used+1)     !|1
           D0(  2,  N_Used+1) = W_tt(2,1,N_Used)  !|
c cc *                        * D0(1,1)           !|1
c cc /                        / F(1,N_Used+1)     !|1
c cc +                        + W_tt(2,2,N_Used)  !|
c cc *                        * D0(2,1)           !|0
c cc /                        / F(2,N_Used+1)     !|1
           D0(N_t+1,N_Used+1) = (0.0D0, 0.0D0)    !|
           D0(N_t+2,N_Used+1) = (0.0D0, 0.0D0)    !|
                                                  !|
        elseif (m_type .eq. 1) Then !--------------+ Asymmetric diffraction, matrix type 1=3x3i
                                                  !| (grazing incidence)
           D0(  1,  1) = (1.0D0, 0.0D0)           !| D0
           D0(N_t+1,1) = Es                       !|
           D0(N_t+2,1) = Eh                       !|
           D0(  1,  N_Used+1) = W_tt(1,1,N_Used)  !|
c cc *                        * D0(1,1)           !|1
c cc /                        / F(1,N_Used+1)     !|1
           D0(N_t+1,N_Used+1) = (0.0D0, 0.0D0)    !|
           D0(N_t+2,N_Used+1) = (0.0D0, 0.0D0)    !|
                                                  !|
        elseif (m_type .eq. 2) Then !--------------+ Asymmetric diffraction, matrix type 2=3x3h
                                                  !| (grazing exit
           D0(  1,  1) = (1.0D0, 0.0D0)           !| D0
           D0(  2,  1) = (0.0D0, 0.0D0)           !|
           D0(N_t+1,1) = Eh                       !|
           D0(  1,  N_Used+1) = W_tt(1,1,N_Used)  !|
c cc *                        * D0(1,1)           !|1
c cc /                        / F(1,N_Used+1)     !|1
c cc +                        + W_tt(1,2,N_Used)  !|
c cc *                        * D0(2,1)           !|0
c cc /                        / F(2,N_Used+1)     !|1
           D0(  2,  N_Used+1) = W_tt(2,1,N_Used)  !|
c cc *                        * D0(1,1)           !|1
c cc /                        / F(1,N_Used+1)     !|1
c cc +                        + W_tt(2,2,N_Used)  !|
c cc *                        * D0(2,1)           !|0
c cc /                        / F(2,N_Used+1)     !|1
           D0(N_t+1,N_Used+1) = (0.0D0, 0.0D0)    !|
                                                  !|
        elseif (m_type .eq. 3) Then !--------------+ Normal diffraction, matrix type 3=2x2d
                                                  !|
           D0(  1,  1) = (1.0D0, 0.0D0)           !| D0
           D0(N_t+1,1) = Eh                       !|
           D0(  1,  N_Used+1) = W_tt(1,1,N_Used)  !|
c cc *                        * D0(1,1)           !|1
c cc /                        / F(1,N_Used+1)     !|1
           D0(N_t+1,N_Used+1) = (0.0D0, 0.0D0)    !|
                                                  !|
        elseif (m_type .eq. 4) Then !--------------+ Specular reflection, matrix type 4=2x2s
                                                  !|
           D0(  1,  1) = (1.0D0, 0.0D0)           !|
           D0(N_t+1,1) = Es                       !|
           D0(  1,  N_Used+1) = W_tt(1,1,N_Used)  !|
c cc *                        * D0(1,1)           !|1
c cc /                        / F(1,N_Used+1)     !|1
           D0(N_t+1,N_Used+1) = (0.0D0, 0.0D0)    !|
                                                  !|
        endif  !-----------------------------------+

        do      k = N_Used-1,1,-1   !================+
           Call Find_Fields_rect (D0(1,1),          !| T_0
     *                            D0(N_t+1,k+2),    !| R_{k+1}
     *                            D0(  1,  k+1),    !| T_k
     *                            D0(N_t+1,k+1),    !| R_k
     *                            W_tt(1,1,k),      !|
     *                            W_tr(1,1,k),      !|
     *                            M_rt(1,1,k+1),    !|
     *                            M_rr(1,1,k+1),    !|
     *                            N_t, N_r, ifail)  !|
           if (ifail.ne.0) goto 100  !---------------+--+
        enddo  !=====================================+  v

c The program Norm_Fields calculates [D_k]*[F_k^(U)], so that if
c we look at the paper PRB 54 (1996) 8150-8162, Eq.(24) and below,
c than a part of Q_{nij} is already included into E_{nij} -- namely
c exp(iku^{IN}_{ni}*Z_n) and exp(iku^{OUT}_{nj}*Z_n. However,
c the term exp(ik*Psi*Z_n) remains and therefore it is included in
c Phase(i) -- see the definition of phase iun the main program!

        do      k = 2, N_Used   !================+
           Call Norm_Fields_rect (D0(  1,  k),  !| T_k
     *                            D0(N_t+1,k),  !| R_k
     *                            F(  1,  k),   !| Ft_k
     *                            F(N_t+1,k),   !| Fr_k
     *                            N_t, N_r)     !|
        enddo  !=================================+

c
c Calculate Dh for diffraction:
c
        if (m_type .lt. 4) Then !-------------------------------+ Diffraction
                                                               !|
           do      k=1,N_Used+1   !==========================+  |
              if (Abs(xe(k)/xabs) .gt. (1.E-6)) Then !----+  |  |Crystal layer
                                                         !|  |  |
                 x0n    = x0(k) / xabs                   !|  |  |
                 Xrh    = Real (xe(k)) / xabs            !|  |  |
                 Xih    = Aimag(xe(k)) / xabs            !|  |  |
                 Xir    = Xih * sin(xf(k)*pi)            !|  |  |
                 Xii    = Xih * cos(xf(k)*pi)            !|  |  |
                 xhn(1) = Cmplx(Xrh+Xir,Xii)       !x(+h) |  |  |
                 xhn(2) = Cmplx(Xrh-Xir,Xii)       !x(-h) |  |  |
                 s02    = g0**2+x0n                      !|  |  |
                                                         !|  |  |
                 do i=1,N_Fields !===========+            |  |  |
                    V = u(i,k)              !|            |  |  |
                    V = (V**2-s02)/xhn(2)   !|            |  |  |
                    Dh(i,k)= V * D0(i,k)    !|            |  |  |
                 enddo  !====================+            |  |  |
                                                         !|  |  |
              else  !-------------------------------------+  |  |Amorphous layer
                                                         !|  |  |
                 if     (m_type .eq. 0) Then !----+       |  |  | GID 4x4
                    Dh(2,k) = D0(2,k)            !|       |  |  |D0: 1,0,1,0
                    D0(2,k) = (0.0D0, 0.0D0)     !|       |  |  |Dh: 0,1,0,1
                    Dh(4,k) = D0(4,k)            !|       |  |  |
                    D0(4,k) = (0.0D0, 0.0D0)     !|       |  |  |
                 elseif (m_type .eq. 1) Then !----+       |  |  | Assym.graz.inci 3x3i
                    Dh(3,k) = D0(3,k)            !|       |  |  |D0: 1,1,0
                    D0(3,k) = (0.0D0, 0.0D0)     !|       |  |  |Dh: 0,0,1
                 elseif (m_type .eq. 2) Then !----+       |  |  | Assym.graz.exit 3x3h
                    Dh(2,k) = D0(2,k)            !|       |  |  |D0: 1,0,0
                    D0(2,k) = (0.0D0, 0.0D0)     !|       |  |  |Dh: 0,1,1
                    Dh(3,k) = D0(3,k)            !|       |  |  |
                    D0(3,k) = (0.0D0, 0.0D0)     !|       |  |  |
                 elseif (m_type .eq. 3) Then !----+       |  |  | Normal diffrac. 2x2d
                    Dh(2,k) = D0(2,k)            !|       |  |  |D0: 1,0
                    D0(2,k) = (0.0D0, 0.0D0)     !|       |  |  |Dh: 0,1
                 endif !--------------------------+       |  |  | Diffraction
                                                         !|  |  |
              endif  !------------------------------------+  |  |
           enddo  !==========================================+  |
                                                               !|
        endif  !------------------------------------------------+

c
c Compute standing waves
c

        do i=1,n_standing  !===========================================+
           z = Positions(i_reference+1)                               !|
     +       + standing_range(1) + (i-1) * standing_step              !|
c Now, find in what layer we are:                                     !|
           do k=N_Used,1,-1  !=================+                       |
              if (Positions(k).le.z) goto 5 !--+--+                    |
           enddo   !===========================+  |                    |
           k  = 0                                !|                    |
  5        continue !<----------------------------+                    |
                                                                      !|
           kk = k                                       !Vacuum fields |--+ exchange these lines
           if (k.lt.1) k=1                                            !|  | to use vacuum or subs-
c ccc      kk = k                                       !Substrate flds|<-+ trate fields in vacuum
c (NOTE: even when using substrate fields, the exponents must belong   |
c to the vacuum fields because otherwise we may get huge numbers in    |
c the vacuum! -- SS 2005-07-29)                                        |
                                                                      !|
c The offset from the closest upper interface:                        !|
           zi   = z - Positions(k)                                    !|
           gi_k = gi / (1.0D0 + Dble(Daa(k)))                         !|
           SUMM = (0.0D0, 0.0D0)                                      !|
           do j=1,N_Fields !========================================+  |
              if (abs(D0(j,kk+1)).gt.1.D-32 .OR.                   !|  |
     *            abs(Dh(j,kk+1)).gt.1.D-32) Then !--------------+  |  |
                 ikuz = Im * k0*u(j,kk+1)*TERang * zi           !|  |  |
                 rkz = dReal(ikuz)                              !|  |  |
                 ikz = dImag(ikuz)                              !|  |  |
                 if     (ikz .gt. 2.*pi) Then !------+           |  |  |
                    j2  = int(ikz/(2.*pi))          !|           |  |  |
                    ikz = ikz - (2.*pi)*j2          !|           |  |  |
                 elseif (ikz .lt. -2.*pi) Then !-----+           |  |  |
                    j2  =-int(ikz/(2.*pi))          !|           |  |  |
                    ikz = ikz + (2.*pi)*j2          !|           |  |  |
                 endif  !----------------------------+           |  |  |
                 if     (rkz.gt.1.D+10) Then !----------------+  |  |  |
                    write(3,*) 'WARNING: too big positive',  !|  |  |  |
     *                         ' i*kz*z=',rkz,               !|  |  |  |
     *                         ' at scan=',Scan_angle,       !|  |  |  |
     *                         ' z=',z                       !|  |  |  |
                    rkz = 1.E+10                             !|  |  |  |
                 elseif (rkz .lt. -1.D+10) Then !-------------+  |  |  |
                    write(3,*) 'WARNING: too big negative',  !|  |  |  |
     *                         ' i*kz*z=',rkz,               !|  |  |  |
     *                         ' at scan=',Scan_angle,       !|  |  |  |
     *                         ' z=',z                       !|  |  |  |
                    rkz = -1.E+10                            !|  |  |  |
                 endif !--------------------------------------+  |  |  |
                 ikuz = dCmplx(rkz,ikz)                         !|  |  |
                 if (standing_phase.ge.0.) Then !---------+      |  |  |
c If standing wave phase is fixed:                        |      |  |  |
                    i2pp = Im * pi * standing_phase      !|      |  |  |
                 else !-----------------------------------+      |  |  |
c If standing wave phase is changing with z:              |      |  |  |
                    hzz = k0*gi_k*TERang * zi            !|      |  |  |
                    if     (hzz .gt. 2.*pi) Then !----+   |      |  |  |
                       j2  = int(hzz/(2.*pi))        !|   |      |  |  |
                       hzz = hzz - (2.*pi)*j2        !|   |      |  |  |
                    elseif (hzz .gt. -2.*pi) Then !---+   |      |  |  |
                       j2  =-int(hzz/(2.*pi))        !|   |      |  |  |
                       hzz = hzz + (2.*pi)*j2        !|   |      |  |  |
                    endif  !--------------------------+   |      |  |  |
                    i2pp = Im * hzz                      !|      |  |  |
                 endif !----------------------------------+      |  |  |
                 SUMM = SUMM + (D0(j,kk+1)                      !|  |  |
     +                       +  Dh(j,kk+1)*exp(i2pp))           !|  |  |
     *                       *  exp(ikuz)                       !|  |  |
              endif  !-------------------------------------------+  |  |
           enddo  !=================================================+  |
           SW(i) = abs(SUMM)**2                                       !|
        enddo  !=======================================================+

        Return

  100   continue
        Return
        end

c ======================================================= 2
        Subroutine Algorithm_Type  (Type)
        Character       Type*(*)
        Type = 'Recursive'             !10
c cccc  Type = 'Matrix'
c cccc  Type = 'Combined'
        return
        end

c ======================================================= 3
        Subroutine DetectLayerTypes (Scan_angle,polfactor,
     *                                  g0,gi,gh,xe,jPol,
     *                                  n_cryst,ifail)
c--------------------------------------------------------
c+===============================+
        Include 'gid_slm7.inc'  !|
c+===============================+
c The above INC file contains smth like:
c ccc   Parameter       (N_Top_Max     = 1001,
c ccc*                   N_Total_Max   = N_Top_Max+1)

c+------------------------------------------------------+
c| ATTENTION: two big local arrays -- alpha and LayType |
c+------------------------------------------------------+

        Complex*16 gh                           !normalized exit angle

        Complex*8  xe(N_Total_Max+1)            !xh_effective

        Real*8     alpha(N_Total_Max+1),        !Bragg deviations per layer
     *             aa_,                         !|alpha|
     *             alpha_min,                   !min alpha among big alpha values
     *             Scan_angle,                  !scan angle
     *             g0,                          !normalized angles
     *             gi,                          !normalized Psi
     *             gi_k                         !gi/(1+Daa)

        Real*4     polfactor

        Integer    LayType(N_Total_Max+1),      !layer types (ver.M5)
     *             jPol,                        !Polarization state
     *             ifail,                       !failure flag
     *             n_cryst,                     !number of cryst.layers
     *             n_large,                     !number of layers with large alpha
     *             n_weak,                      !number of layers with weak refleactions
     *             i

c----------------------
        Complex*8  x0(N_Total_Max+1),
     *             xh(N_Total_Max+1)

        Real*4     xf(N_Total_Max+1),           !phase shift between xrh & xih (non-cubic)
     *             xs(N_Total_Max+1),           !M7: xh phase shift in layer (Kaganer, monolayers)
     *             Thickness(N_Total_Max),      !thicknesses of layers
     *             Sigma(N_Total_Max),          !rms rougness of upper interface
     *             Daa(N_Total_Max+1),          !normal lattice strain
     *             Wave,                        !X-wavelength
     *             xabs,                        !abs(x0_max)
     *             TERang,                      !critical TER angle
     *             standing_range(2),           !M7: offsets in A from reference interface
     *             standing_step,               !M7: SW depth step in A
     *             standing_phase,              !M7: SW phase in units of pi
     *             xh_weak_limit,               !M5: the weakest possible xh
     *             alpha_max                    !M5: max(alpha/xh) to treat crystal as "crystal"

        Integer    N_Top,                       !top layer sublayers
     *             N_Total,                     !total of layers
     *             i_standing,                  !M7: 1=yes/0=no SW
     *             i_reference,                 !M7: reference interface for SW (0=surface)
     *             n_standing,                  !M7: Number of SW depth points
     *             m_reduction,                 !M7: matrix size reduction flag (0/1/2)
     *             m_type,                      !M7: matrix type 4x4,3x3a,3x3b,2x2a,2x2b.
     *             iDebug                       !Debug flag

c This is common between all 3 programs:
c gid_sl, gid_in and RFM:
        Common /GEN_M7/ x0, xh, xf, xs,
     *                  Thickness, Wave,
     *                  Sigma, Daa,
     *                  xabs, TERang,
     *                  xh_weak_limit,          !M5
     *                  alpha_max,              !M5
     *                  standing_range,         !M7: new
     *                  standing_step,          !M7: new
     *                  standing_phase,         !M7: new
     *                  i_standing,             !M7: new
     *                  i_reference,            !M7: new
     *                  n_standing,             !M7: new
     *                  m_reduction,            !M7: new
     *                  m_type,                 !M7: new
     *                  N_Top, N_Total,
     *                  iDebug

c----------------------
        i = N_Top_Max                           !to suppress unused paramters warnings
        i = N_Total_Max                         !to suppress unused paramters warnings
        i = N_Pts_Max                           !to suppress unused paramters warnings
        i = N_Stand_Max                         !to suppress unused paramters warnings
c----------------------
        ifail     = 0
        n_cryst   = 0
        n_large   = 0
        n_weak    = 0
        alpha_min = 1.0D+37
        do i=1,N_Total+1  !=======================================+
          xe(i) = xh(i) * polfactor                              !|
          if (abs(xe(i)) .lt. xh_weak_limit) Then !------------+  |
c Amorphous layer because of weak reflection:                 !|  |
            xe(i)      = (0.0, 0.0)                           !|  |
            alpha(i)   = 0.0                                  !|  |
            if (abs(xh(i)).gt. 1.E-20) Then !-------------+    |  |
              LayType(i) = -2   !weak x0,treat as amorph. |    |  |
              n_weak     = n_weak + 1                    !|    |  |
            else !----------------------------------------+    |  |
              LayType(i) = 0    !amorphous layer          |    |  |
            endif !---------------------------------------+    |  |
          else !-----------------------------------------------+  |
            gi_k  = gi / (1.0D0 + Dble(Daa(i)))               !|  |
            alpha(i) = (g0+gi_k)*(g0+gi_k) - Real(gh*gh)      !|  |
c Express alpha in times (xh) of particular layer:             |  |
            alpha(i) = alpha(i) * (xabs/Abs(xe(i)))           !|  |
            aa_      = Abs(alpha(i))                          !|  |
            if (aa_ .gt. alpha_max) Then !----------------+    |  |
c Amorphous layer because of large alpha:                !|    |  |
              xe(i)      = (0.0, 0.0)                    !|    |  |
              LayType(i) = -1   !large alpha,treat as amor|    |  |
              n_large    = n_large+1                     !|    |  |
              if (aa_ .lt. alpha_min) alpha_min = aa_    !|    |  |
            else !----------------------------------------+    |  |
              LayType(i) = 1    !crystalline layer        |    |  |
              n_cryst    = n_cryst+1                     !|    |  |
            endif  !--------------------------------------+    |  |
          endif  !---------------------------------------------+  |
        enddo !===================================================+

c Warn about weak diffracting layers:
        if (n_weak .gt. 0) Then !-------------------------------------+
          write (3,5) n_weak,n_cryst+n_weak+n_large,Scan_angle,jPol  !|
  5       format (' --- DetectLayerTypes WARNING: ',i5,' of ',i5,
     *    1x,'crystalline layers approximated amorphous',
     *    1x,'because of small |xh| at scan_angle = ',g12.5,
     *    1x,'polarization=',i1/
     *    '     The layers numbers are: ',$)                         !|Suppress CR (MS/Compaq/Gnu)
          do i=1,N_Total+1 !===========================+              |
            if (LayType(i) .eq. -2) write (3,6) i-1   !|              |
  6         format(1x,i5,$)                           !|              |Suppress CR (MS/Compaq/Gnu)
          enddo !======================================+              |
          write(3,4)                                                 !|
        endif !-------------------------------------------------------+

c If we eliminated all crystalline layers, then we have to restore
c at least those that have minimum alpha in the list:
        if (n_cryst .eq. 0 .and. n_large .gt. 0)  Then !--------------+
          write (3,7) alpha_max, Scan_angle, jPol                    !|
  7       format (' --- DetectLayerTypes WARNING: some crystalline',
     *    1x,'layer(s) have large |alpha|>|alpha_max|=',e10.4,', but',
     *    1x,'still treated as crystalline at scan_angle = ',g12.5,
     *    2x,'polarization=',i1/
     *    '     The layers numbers are: ',$)                         !|Suppress CR (MS/Compaq/Gnu)
          do i=1,N_Total+1 !=======================================+  |
            if (LayType(i) .eq. -1) Then !-----------------------+ |  |
              aa_ = Abs(alpha(i))                               !| |  |
              if (Abs(aa_/alpha_min-1.) .lt. 1.E-5) Then !-----+ | |  |
                xe(i) = xh(i) * polfactor                     !| | |  |
                LayType(i) = 1          !revert to crystal.    | | |  |
                n_cryst    = n_cryst+1                        !| | |  |
                n_large    = n_large-1                        !| | |  |
                write (3,3) i-1, alpha(i)                     !| | |  |
              endif  !-----------------------------------------+ | |  |
            endif  !---------------------------------------------+ |  |
          enddo  !=================================================+  |
          write(3,4)                                                 !|
        endif !-------------------------------------------------------+

c Check if that helped (should!)
        if (n_cryst .eq. 0)  Then !-----------------------------------+
          write (3,1) Scan_angle                                     !|
  1       format (
     *    ' *** DetectLayerTypes ERROR:'/
     *    '  converting crystalline layers to amorphous led to'/
     *    '  no more diffracting layers left at scan_angle=',g12.5/
     *    '  Please reduce scan range or try increasing Alpha_max!') !|
          ifail = 777                                                !|
          return                                                     !|
        endif !-------------------------------------------------------+

c Warn about approximated layers:
        if (n_large .gt. 0) Then !------------------------------------+
          write (3,2) n_large,n_cryst+n_weak+n_large,Scan_angle,jPol !|
  2       format (' --- DetectLayerTypes WARNING: ',i5,' of ',i5,
     *    1x,'crystalline layers approximated amorphous because of',
     *    1x,'large Bragg deviation alpha at scan_angle = ',g12.5,
     *    1x,'polarization=',i1/
     *    '     The layers numbers are: ',$)                         !|Suppress CR (MS/Compaq/Gnu)
          do i=1,N_Total+1 !====================================+     |
            if (LayType(i) .eq. -1) write (3,3) i-1, alpha(i)  !|     |
  3         format(1x,i5,' (alpha=',e10.4,')',$)               !|     |Suppress CR (MS/Compaq/Gnu)
          enddo !===============================================+     |
          write(3,4)                                                 !|
  4       format(//)                                                 !|
        endif !-------------------------------------------------------+
        return
        end
