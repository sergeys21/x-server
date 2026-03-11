        Subroutine Get_X0h (
     *    Wave,
     *    Code,
     *    rho,      !density (g/cm^3 -- for amorphous structutes only)
     *    LayType,  !0/1=number of reflections to calculate
     *    indices,
     *    QB,       !define, if QB=0; compare with prev., if QB#0
     *    mode_x0h, !0-5, actually 0 or 1 (may be 2)
     *    Xr0,Xi0,
     *    Xrh,Xih,
     *    X0,XH,
     *    N_Layers,
     *    istep,
     *    ifail)
c--------------------------------------------------------
c This subroutine is called from gids_xx (may be obsolete! - see GetXh97.for)
c--------------------------------------------------------
        Integer         N_Layers, istep, ifail,
     *                  LayType, indices(3), mode_x0h,
     *                  isyng, iwait, i, j,
     *                  n_atoms_ucell

        Complex         X0(1+(N_Layers-1)*istep),
     *                  XH(1+(N_Layers-1)*istep)

        Real            a_lattice(6), QB, QB_, Wave, rho,
     *                  Xr0, Xi0, Xrh, Xih,
     *                  xqh, xrf, xif, xqf,
     *                  ucell_mass_gram, atom_density_cm3

        Character       Code*20

        Character       txt(20)*80
        Common  /msg/   txt

c This is an interface to internal X0h data.
        Common  /x0paA/ ucell_mass_gram, atom_density_cm3, n_atoms_ucell
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Get_X0h'            !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        ifail    = 0
c -------------------------------------------------------
c Possible values of mode_x0h (see x0h1.for)
c       md_NewBase     =-1
c       md_NewAll      = 0
c       md_NewCrystal  = 1      !use new name,prcn,a(1-6) from X0h DB
c       md_NewCrystal_a= 2      !use new name,prcn from X0h DB, use ad(1-6) if non-zero, otherwise from X0h DB
c       md_NewWave     = 3
c       md_NewLine     = 4
c       md_NewHKL      = 5
c Subr. x0h1 will copy this to common /back1/ and thus effectively
c set a new crystal:
        do i=1,6  !===========+
          a_lattice(i) = 0.  !|
        enddo !===============+
        if (mode_x0h.le.2) Then  !-----+
          ucell_mass_gram  = 0.       !|
          atom_density_cm3 = 0.       !|
          n_atoms_ucell    = 0        !|
        endif  !-----------------------+
        Xr0 = 0.
        Xi0 = 0.
        Xrh = 0.
        Xih = 0.
        xqh = 0.
        xrf = 0.
        xif = 0.
        xqf = 0.
        if (Code.ne.'Unknown')  Then  !------------------+
          isyng   = 1                   !cubic crystal   |
          Call X0h1  (a_lattice,Wave,Code,rho,LayType,  !|
     *                indices(1),indices(2),            !|
     *                indices(3),isyng,                 !|
     *                mode_x0h,0,QB_,                   !|
     *                Xr0,Xi0,                          !|
     *                Xrh,Xih,                          !|
     *                xqh,xrf,xif,xqf,ifail)            !|
          if (ifail.ne.0)  goto 100                     !|
c At next call the crystal might be not the same:        |
          mode_x0h = 1                                  !|
        endif  !-----------------------------------------+

c What to do at absorption edges
c where xr0 may change the sign?
c This however should not be a
c problem as long as we use x0h
c since x0h does not process the
c edges!
c       X0(1) = Cmplx(-abs(Xr0),abs(Xi0))
c Changed 2003/07/25:
c We now have control over manual input in Xinput97 -> ChiSearch
c and X0h should calculate correctly!
c Hope this is safe!
        X0(1) = Cmplx(Xr0,abs(Xi0))


        if (LayType.eq.1)  Then  !----------------------+
c This is an old definition corresponding to phase=pi:  |
c         XH(1) = Cmplx(-abs(Xrh),abs(Xih))            !|
c Starting from June-98 in Get_Xh97 we choose:          |Changed in this program 07/2003!
          XH(1) = Cmplx(abs(Xrh),abs(Xih))             !|
c -- then, the real Xh and Xh_ are:                     |
c         Xh = Xrh+Xih*sin(Xdf*pi) + i*Xih*cos(Xdf*pi)  |
c         Xh_= Xrh-Xih*sin(Xdf*pi) + i*Xih*cos(Xdf*pi)  |
c -- where the sign before sin corresponds to our       |
c "inverted" Xdf=Xrf-Xif  !!!!!!!!!!!!!!!!!!!!!         |
          if (abs(QB).lt.1.E-32) Then !---------------+ |
            QB=QB_                                   !| |
          else  !-------------------------------------| |
            if (abs(QB-QB_).gt.1.)  Then  !---------+ | |
              txt(1)=progname//' WARNING: Bragg'// !| | |
     *               ' angle mismatch for '//Code  !| | |
              if (modebat.ne.0) Then !--+           | | |
                iwait = 0              !|           | | |
              else !--------------------|           | | |
                iwait = 2              !|           | | |
              endif  !------------------+           | | |
              Call Message(txt,1,iwait)            !| | |
            endif  !--------------------------------+ | |
          endif  !------------------------------------+ |
        else  !-----------------------------------------|
          XH(1) = (0.,0.)                              !|
        endif  !----------------------------------------+

        j = 1
        do      i=2,N_Layers  !=====+
          j = j + istep            !|
          X0(j) = X0(1)            !|
          XH(j) = XH(1)            !|
        enddo   !===================+

  100   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        Return
        End
