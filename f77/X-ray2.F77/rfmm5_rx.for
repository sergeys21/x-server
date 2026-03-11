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
c
c +=====================================================+
c |This  routine computes GID reflection coefficients   |
c |Rh, R0 from multilayered structure at 1 angular point|
c |   ( Incident angle = f0,     Exit angle = fh ).     |
c |Input data are transferred to RFM in common/RFM5dat/ |
c |-----------------------------------------------------|
c |Contents:                                            |
c | 1. Real*4 function:       RFM                       |
c | 2. Subroutine:            Algorithm_Type            |
c |-----------------------------------------------------|
c |Subroutine calls: see  --  RFMs_Uni.for              |
c +=====================================================+

c ======================================================= 1
        Real*8  Function        RFM8_m5 (alpha)
c  ** Reflex_from_Multi **
c--------------------------------------------------------
c+===============================+
        Include 'gid_slm7.inc'  !|
c+===============================+
c The above INC file contains smth like:
c ccc   Parameter       (N_Top_Max     = 1001,
c ccc*                   N_Total_Max   = N_Top_Max+1)

        Integer          N_Roots,
     *                   N_Half,
     *                   N_Fields
        Parameter       (N_Roots  = 4,
     *                   N_Half   = N_Roots/2)
        Parameter       (N_Fields = 4)

        Complex*16 X(N_Roots,N_Roots,N_Total_Max), !X-matrix of interface
     *             F(N_Roots,N_Total_Max+1),       !F-matrix
     *             u(N_Roots),                     !dispersion eqn.roots
     *             u1(N_Roots),                    !     --"--  (previous)
     *             S(N_Roots,N_Roots),             !S_matrix 4x4
     *             B(N_Roots,N_Roots),
     *             B1(N_Roots,N_Roots),            !1/S_matrix
     *             TTT(N_Roots,N_Roots),           !T1*...Tn*Ss*exp(ikuz)
     *             M_rr(N_Half,N_Half),            !product of operators
     *             M_rt(N_Half,N_Half),            !product of operators
     *             M_tr(N_Half,N_Half),            !product of operators
     *             M_tt(N_Half,N_Half),            !product of operators
     *             W_rr(N_Half,N_Half),            !product of operators
     *             W_rt(N_Half,N_Half),            !product of operators
     *             W_tr(N_Half,N_Half),            !product of operators
     *             W_tt(N_Half,N_Half),            !product of operators
     *             Determinant,                    !work sell
     *             x0n,xhn(2),                     !normalized x0,xh,x(-h),
     *             Es, Eh,                         !vacuum wave amplitudes
     *             gh                              !normalized exit angle
        Complex*8  xe(N_Total_Max+1)               !xh_effective
        Real*8     g0, gh_r, gh2,                  !normalized angles
     *             alpha,                          !Bragg deviation
     *             gi,                             !norm-zed Psi
     *             k0,                             !k0 = 2*pi/wave
     *             Upper_Limit                     !max matr.elmnt (lim)
        Real*4     pi,                             !3.1416
     *             Xrh,Xih,Xir,Xii,                !normalized parts of xh
     *             Sgm_i                           !Sigma of interface
        Integer    Kod,                            !Matrix(0)/Recursive(1) flag
     *             n_cryst,                        !Number of crystalline layers after approximation
     *             i

c----------------------
        Complex*8  x0(N_Total_Max+1),
     *             xh(N_Total_Max+1)
        Real*4     xf(N_Total_Max+1),           !phase shift between xrh & xih (non-cubic)
     *             Thickness(N_Total_Max),      !thicknesses of layers
     *             Sigma(N_Total_Max),          !rms rougness of upper interface
     *             Daa(N_Total_Max+1),          !normal lattice strain
     *             Wave,                        !X-wavelength
     *             xabs,                        !abs(x0_max)
     *             TERang,                      !critical TER angle
     *             xh_weak_limit,               !M5: the weakest possible xh
     *             alpha_max                    !M5: max(alpha/xh) to treat crystal as "crystal"
        Integer    N_Top,                       !top layer sublayers
     *             N_Total,                     !total of layers
     *             MaxOrder,                    !max order of mat.elem.
     *             iDebug                       !Debug flag
c This is common between all 3 programs:
c gid_sl, gid_in and RFM:
        Common /GEN_M5/ x0, xh, xf,
     *                  Thickness, Wave,
     *                  Sigma, Daa,
     *                  xabs, TERang,
     *                  xh_weak_limit,          !M5
     *                  alpha_max,              !M5
     *                  N_Top, N_Total,
     *                  MaxOrder, iDebug

c----------------------
        Real*8     Positions(N_Total_Max),      !depth of interfaces
     *             Phase(N_Total_Max),          !int(d(h_z*dz))
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
        Common /RFM5dat/ R0, Positions, Phase,
     *                   f0, Scan_angle,
     *                   Psi, Pol,
     *                   Account_Sigma,
     *                   N_Used, jPol, ifail

c----------------------
        i = N_Top_Max                           !to suppress unused paramters warnings
        i = N_Total_Max                         !to suppress unused paramters warnings
        i = N_Pts_Max                           !to suppress unused paramters warnings
        i = N_Stand_Max                         !to suppress unused paramters warnings
c--------------------------------------------------------
c       xabs   = Abs(x0_max)                    !defined in main
c       TERang = Sqrt(xabs)                     !defined in main
        pi     = 4.* Atan(1.)
        k0     = 2.* pi/wave
        g0     = f0 / TERang
        gi     = Psi/ TERang
        gh2    = (g0+gi)**2 - alpha/xabs
        if (gh2 .ge. 0.) Then !---------------+
          gh   = Dsqrt(gh2)                  !|
        else   !------------------------------|
          gh   = Dcmplx(0.0D0,Dsqrt(-gh2))   !|
        endif  !------------------------------+

c+------------------------------------+
        MaxOrder    = 7              !|
        Upper_Limit = 10.**MaxOrder  !|
c+------------------------------------+

        R0      = 0.0D0
        RFM8_m5 = 0.0D0
        ifail   = 0
        Sgm_i   = 0.

c 2005/04/09: added for version gid_slm5.
c Detect what layers must be approximated
c as amorphous because of large alpha.
c The most of paramters are transferred
c via common /GEN_M5/. This subroutine
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

c AT START, use matrix method (before 2020.04.21):
c       Kod = 0
c AT START, use recursive matrix method (after 2020.04.21):
        Kod = 0
c Set TTT as diagonal unit matrix (see rfms_uni.for):
        Call    M1 (N_Roots,TTT,N_Fields)

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
                                                            !|
        Call Process_Layer (g0, gh, gi,                     !|see rfms_urd.for
     *                      Daa(N_Used+1),                  !|
     *                      x0n, xhn, k0,                   !|
     *                      0.,                 !Thickness(0)|
     *                      TERang, S, B, u,                !|
     *                      F(1,N_Used+1),                  !|
     *                      N_Fields,ifail)                 !|
        if (ifail .ne. 0) goto 100                          !|
c+-----------------------------------------------------------+

c
c                     Process Top Layer
c
        do    N_Used=1,N_Total  !============================+
                                                            !|
          Call Copy_Layer (B,B1,u,u1,N_Fields)              !|
                                                            !|
          x0n    = x0(N_Used+1) / xabs                      !|
          Xrh    = Real (xe(N_Used+1)) / xabs               !|
          Xih    = Aimag(xe(N_Used+1)) / xabs               !|
          Xir    = Xih * sin(xf(N_Used+1)*pi)               !|
          Xii    = Xih * cos(xf(N_Used+1)*pi)               !|
          xhn(1) = Cmplx(Xrh+Xir,Xii)              !x(+h)   !|
          xhn(2) = Cmplx(Xrh-Xir,Xii)              !x(-h)   !|
                                                            !|
          Call Process_Layer (g0, gh, gi,                   !|see rfms_urd.for
     *                        Daa(N_Used+1),                !|
     *                        x0n, xhn, k0,                 !|
     *                        Thickness(N_Used),            !|
     *                        TERang, S, B, u,              !|
     *                        F(1,N_Used+1),                !|
     *                        N_Fields,ifail)               !|
          if (ifail .ne .0)   goto 100                      !|
          if (Account_Sigma) Sgm_i = Sigma(N_Used)          !|
          Call Process_Interface(k0, TERang, Sgm_i,         !|
     *                           B1, S, u1, u,              !|
     *                           F(1,N_Used+1),             !|
     *                           X(1,1,N_Used),             !|
     *                           TTT,N_Fields,              !|
     *                           Upper_Limit,               !|
     *                           Kod,ifail)                 !|
c Here ifail=-1 is overflow. At this Kod is automatically    |
c switched to Kod=1 (== do not multiply X-matrices), and     |
c we continue to calculate x-matrices of all the layers.     |
c Then, the interfaces will be processed by the recursive    |
c method (see Process_Interface2 below):                     |
          if (ifail .eq. -1)  ifail = 0                     !|
          if (ifail .ne. 0)   goto 100  !--------------------+----+
        enddo  !=============================================+    V

c
c                   Find solutions Es, Eh
c

        if (Kod .eq. 0) Then !-------------------------------+
                                                            !|
c+---------------------------------+                         |
c|No overflow -- use matrix method:|                         |
c+=================================+                         |
          Determinant = TTT(1,1)*TTT(2,2)                   !|
     -                - TTT(1,2)*TTT(2,1)                   !|
          if (Abs(Determinant) .lt. 1.E-32) Then !--+        |
            ifail = 123                            !|        |
            goto 100                               !|        |
          endif  !----------------------------------+        |
          Es = (TTT(3,1)*TTT(2,2)-TTT(3,2)*TTT(2,1))        !|
     /       / Determinant                                  !|
          Eh = (TTT(4,1)*TTT(2,2)-TTT(4,2)*TTT(2,1))        !|
     /       / Determinant                                  !|
                                                            !|
        else  !----------------------------------------------|
                                                            !|
c+-------------------------------------------+               |
c|Overflow -- switch to recursive equations: |               |
c+===========================================+               |
          do  N_Used=1,N_Total  !======================+     |
            Call Process_Interface2 (X(1,1,N_Used),   !|     |
     *                               F(1,N_Used+1),   !|     |
     *                               M_rr, M_rt,      !|     |
     *                               M_tr, M_tt,      !|     |
     *                               W_rr, W_rt,      !|     |
     *                               W_tr, W_tt,      !|     |
     *                               N_Fields, N_Used,!|     |
     *                               ifail)           !|     |
            if (ifail .ne. 0) goto 100  !--------------+-----+--+
          enddo  !=====================================+     |  V
                                                            !|
          Es = W_rt(1,1)                                    !|
          Eh = W_rt(2,1)                                    !|
        endif  !---------------------------------------------+

c Compute the reflection coefficients:

        R0  = Abs(Es)**2
        if (g0 .gt. 0.) then !------------------+
          gh_r = Dreal(gh)                     !|
          RFM8_m5 = (gh_r/g0) * (Abs(Eh)**2)   !|
        endif  !--------------------------------+
        if (RFM8_m5 .lt. 0.0) then  !-----------+
          write (3,*) 'RFM8_m5: Rh=', RFM8_m5  !|
          ifail = 990                          !|
          goto 100  !---------------------------+---------------+
        endif  !--------------------------------+               v
        if (RFM8_m5 .gt. 1.0) then  !-----------+
          write (3,*) 'RFM8_m5: Rh=', RFM8_m5  !|
          ifail = 991                          !|
          goto 100  !---------------------------+---------------+
        endif  !--------------------------------+               v
        if (R0 .lt. 0.0) then !-----------------+
          write (3,*) 'RFM8_m5: R0=', R0       !|
          ifail = 880                          !|
          goto 100  !---------------------------+---------------+
        endif  !--------------------------------+               v
        if (R0 .gt. 1.0) then !-----------------+
          write (3,*) 'RFM8_m5: R0=', R0       !|
          ifail = 881                          !|
          goto 100  !---------------------------+---------------+
        endif  !--------------------------------+               v

        Return
  100   continue
        Return
        end

c ======================================================= 2
        Subroutine Algorithm_Type  (Type)
        Character       Type*(*)
c cccc  Type = 'Recursive'             !10
c cccc  Type = 'Matrix'
        Type = 'Combined'
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
     *             Thickness(N_Total_Max),      !thicknesses of layers
     *             Sigma(N_Total_Max),          !rms rougness of upper interface
     *             Daa(N_Total_Max+1),          !normal lattice strain
     *             Wave,                        !X-wavelength
     *             xabs,                        !abs(x0_max)
     *             TERang,                      !critical TER angle
     *             xh_weak_limit,               !M5: the weakest possible xh
     *             alpha_max                    !M5: max(alpha/xh) to treat crystal as "crystal"
        Integer    N_Top,                       !top layer sublayers
     *             N_Total,                     !total of layers
     *             MaxOrder,                    !max order of mat.elem.
     *             iDebug                       !Debug flag
c This is common between all 3 programs:
c gid_sl, gid_in and RFM:
        Common /GEN_M5/ x0, xh, xf,
     *                  Thickness, Wave,
     *                  Sigma, Daa,
     *                  xabs, TERang,
     *                  xh_weak_limit,          !M5
     *                  alpha_max,              !M5
     *                  N_Top, N_Total,
     *                  MaxOrder, iDebug

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
            if (abs(xh(i)).gt. 1.E-20) Then !---+              |  |
              LayType(i) = -2      !approx.weak |              |  |
              n_weak     = n_weak + 1          !|              |  |
            else !------------------------------+              |  |
              LayType(i) = 0         !amorphous |              |  |
            endif !-----------------------------+              |  |
          else !-----------------------------------------------|  |
            gi_k  = gi / (1.0D0 + Dble(Daa(i)))               !|  |
            alpha(i) = (g0+gi_k)*(g0+gi_k) - Real(gh*gh)      !|  |
c Express alpha in times (xh) of particular layer:             |  |
            alpha(i) = alpha(i) * (xabs/Abs(xe(i)))           !|  |
            aa_      = Abs(alpha(i))                          !|  |
            if (aa_ .gt. alpha_max) Then !-------------------+ |  |
c Amorphous layer because of large alpha:                   !| |  |
              xe(i)      = (0.0, 0.0)                       !| |  |
              LayType(i) = -1                               !| |  |
              n_large    = n_large+1                        !| |  |
              if (aa_ .lt. alpha_min) alpha_min = aa_       !| |  |
            else !-------------------------------------------+ |  |
              LayType(i) = 1                                !| |  |
              n_cryst    = n_cryst+1                        !| |  |
            endif  !-----------------------------------------+ |  |
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
        if (n_cryst .eq. 0)  Then !-----------------------------------+
          write (3,7) Scan_angle, jPol                               !|
  7       format (' --- DetectLayerTypes WARNING: some crystalline',
     *    1x,'layers have large |alpha|>|alpha_max|, but still',
     *    1x,'treated as crystalline at scan_angle = ',g12.5,
     *    2x,'polarization=',i1/
     *    '     The layers numbers are: ',$)                         !|Suppress CR (MS/Compaq/Gnu)
          do i=1,N_Total+1 !====================================+     |
            if (LayType(i) .eq. -1) Then !--------------------+ |     |
              if (Abs(alpha(i)-alpha_min).lt.1.E-20) Then !-+ | |     |
                xe(i) = xh(i) * polfactor                  !| | |     |
                LayType(i) = 1                             !| | |     |
                n_cryst    = n_cryst+1                     !| | |     |
                n_large    = n_large-1                     !| | |     |
                write (3,3) i-1, alpha(i)                  !| | |     |
              endif  !--------------------------------------+ | |     |
            endif  !------------------------------------------+ |     |
          enddo  !==============================================+     |
          write(3,4)                                                 !|
        endif !-------------------------------------------------------+

c Check if that helped (should!)
        if (n_cryst .eq. 0)  Then !-----------------------------------+
          write (3,1) Scan_angle                                     !|
  1       format (' *** DetectLayerTypes ERROR: converting',1x,
     *    'crystalline layers to amorphous failed at scan_angle=',
     *                                                       g12.5)  !|
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
  3         format(1x,i5,'(a=',e10.4,')',$)                    !|     |Suppress CR (MS/Compaq/Gnu)
          enddo !===============================================+     |
          write(3,4)                                                 !|
  4       format(//)                                                 !|
        endif !-------------------------------------------------------+
        return
        end
