c This module includes the following subroutines:
c   1. recalcXh
c -- which recalculates X0,Xh for the whole stack of layers,
c including substrate, for gid_slM5 and up. It is meant to
c be used with energy scanning where Xh/X0 changes during
c the scan are essential.
c
c =======================================================

        Subroutine recalcXh (Wave, indices,
     *                       N_Total, N_Frac_Max,
     *                       N_Frac,
     *                       Code, Frac,
     *                       x0, xh, xf,
     *                       W0, Wh,
     *                       C0h, QB_new, ifail)
c--------------------------------------------------------
c       Recalculating  x0,xh of top layer profile -- smart version!
c
c                    Author S.Stepanov
c -------------------------------------------------------
        Integer         N_Total, N_Frac_Max,
     *                  N_Frac(N_Total),
     *                  indices(3),
     *                  ifail, nreflect,
     *                  mode_x0h, iSyngony,
     *                  l, i, j

        Complex         x0(N_Total+1),
     *                  xh(N_Total+1),
     *                  x0_f, xh_f

        Real            Wave,
     *                  xf(N_Total+1),
     *                  W0(N_Total),
     *                  Wh(N_Total),
     *                  Frac(N_Frac_Max,N_Total),
     *                  xr0, xi0, xrh, xih, xf_f,
     *                  rho, QB_new, QB

        Character       Code(N_Frac_Max,N_Total)*(*),
     *                  C0h(N_Total)*(*)

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv .lt. 10) Then  !----+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'recalcXh'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
c -------------------------------------------------------
c Get x0,xh for the substrate:
        mode_x0h = 0            ! 1st calculation (new crystal,...)
        QB       = 0.           ! No control on the coincidence of Bragg angles
        iSyngony = 0            ! any syngony!
        nreflect = 1
        rho      = 0.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Call Get_Xh97  (Wave,                           !input, but may be output: if wave=0
     *                  Code(1,N_Total),                !input
     *                  rho, iSyngony,                  !output (for amorphous only)
     *                  nreflect, indices,              !input
     *                  QB,                             !input/output: if QB=0, calc. it
     *                  mode_x0h,                       !input: 0-5 (normally 1st call mode=0
     *                  xr0, xi0,                       !output
     *                  xrh, xih, xf(N_Total+1),        !output
     *                  x0(N_Total+1), xh(N_Total+1),   !output
     *                  ifail)                          !output
        if (ifail.ne.0)   goto 130 !----------------------+
                                                         !v
        QB_new = QB
c Added 2011/03/31: backup substrate lattice for Algrebra routines
c so that we can restore it after calling X0h on the profile stack:
        Call    LatticeBackup ()
        do      l=1,N_Total-1  !=======================+
c This program cannot work if x0 or xh are specified   |
c implicitly (i.e. not via Code and X0h database):     |
          if (Len_Trim(C0h(l)).ne.0) goto 403 !--------+---+
          x0(l+1) = (0., 0.)                          !|   v
          xh(l+1) = (0., 0.)                          !|
          xf(l+1) = 0.                                !|
          do j=1,N_Frac(l)  !========================+ |
            if (Code(j,l).eq.Code(1,N_Total)) Then !+| |
              x0_f = x0(N_Total+1)                 !|| |
              xh_f = xh(N_Total+1)                 !|| |
              xf_f = xf(N_Total+1)                 !|| |
            else !----------------------------------|| |
c No control on the coincidence of Bragg angles:    || |
              QB = 0.                              !|| |
c New crystal, but same reflection & wavelength:    || |
              mode_x0h = 1                         !|| |
              iSyngony = 0   !any syngony!          || |
              nreflect = 1                         !|| |
              rho      = 0.                        !|| |
              Call Get_Xh97 (Wave, Code(j,l), rho, !|| |
     *                    iSyngony, nreflect,      !|| |
     *                    indices, QB, mode_x0h,   !|| |
     *                    xr0, xi0,                !|| |
     *                    xrh, xih, xf_f,          !|| |
     *                    x0_f, xh_f, ifail)       !|| |
              if (ifail.ne.0)   goto 130 !----------++-+---+
            endif !---------------------------------+| |   v
            x0(l+1) = x0(l+1) + x0_f*Frac(j,l)      !| |
            xh(l+1) = xh(l+1) + xh_f*Frac(j,l)      !| |
            xf(l+1) = xf(l+1) + xf_f*Frac(j,l)      !| |
          enddo  !===================================+ |
                                                      !|
          x0(l+1) = x0(l+1)*W0(l)                     !|
          xh(l+1) = xh(l+1)*Wh(l)                     !|
        enddo  !=======================================+
c Added 2011/03/31: restore substrate lattice for Algrebra routines
        Call    LatticeRestore ()
c Corrections for substrate:
        x0(N_Total+1) = x0(N_Total+1)*W0(N_Total)
c We must have non-zero absorption in the substrate!
        if (abs(Aimag(x0(N_Total+1))).lt.1.E-32) goto 402  !----------+
        xh(N_Total+1) = xh(N_Total+1)*Wh(N_Total)                    !v

c Test for x0,xh data consistence:
        do      l=1,N_Total  !========================+
          if (Abs(xh(l+1))        .ge.               !|
     *        Abs(x0(l+1)))        goto 400 !---------+--+
          if (Abs(Aimag(xh(l+1))) .ge.               !|  v
     *        Abs(Aimag(x0(l+1)))) goto 401 !---------+--+
        enddo  !======================================+  v

        if (iirezv .eq. 1)  Then  !--------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return

c=================================================================
c                           E R R O R S
c-----------------------------------------------------------------
  400   continue
        i = Len_Trim(progrezv(1))
        j = Len_Trim(progname)
        write   (txt,500)  progrezv(1)(1:i), progname(1:j), l
  500   format  (a,' -> ',a,':'//
     *  'Non-physical input for layer',i5,
     *                                  ': |xh| exceeds or equals |x0|')
        ifail = 3
        goto 130

  401   continue
        i = Len_Trim(progrezv(1))
        j = Len_Trim(progname)
        write   (txt,501)  progrezv(1)(1:i), progname(1:j), l
  501   format  (a,' -> ',a,':'//
     *  'Non-physical input for layer',i5,
     *                              ': Im(xh) exceeds or equals Im(x0)')
        ifail = 3
        goto 130

  402   continue
        i = Len_Trim(progrezv(1))
        j = Len_Trim(progname)
        write   (txt,502)  progrezv(1)(1:i), progname(1:j),
     *                     x0(N_Total+1),W0(N_Total)
  502   format  (a,' -> ',a,':'//
     *  'Non-physical input: zero absorption in the substrate!'//
     *  'x0=(',g12.5,', ',g12.5,'),   w0=',g12.5)
        ifail = 5
        goto 130

  403   continue
        i = Len_Trim(progrezv(1))
        j = Len_Trim(progname)
        write   (txt,503)  progrezv(1)(1:i), progname(1:j), l
  503   format  (a,' -> ',a,':'//
     *  'Explicit data for x0,xh in layer',i5,
     *               ': Cannot calculate wavelength dependence')
        ifail = 3
        goto 130

  130   continue
        return
        end
