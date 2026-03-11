        Subroutine      Get_X0  (Wave, Code, rho,
     *                           isyng, mode_x0h,
     *                           Xr0, Xi0,
     *                           X0, N_Layers, N_Step,
     *                           ifail)
c--------------------------------------------------------
c This subroutine is called from ter_slxx and trds_xx
c--------------------------------------------------------
        Integer         N_Layers, N_Step, ifail,
     *                  isyng, mode_x0h, LayType,
     *                  jpr, i, j,
     *                  n_atoms_ucell

        Complex         X0(N_Step*N_Layers)

        Real            a_lattice(6), Wave, rho, QB,
     *                  Xr0, Xi0, Xrh, Xih,
     *                  xqh, xrf, xif, xqf,
     *                  ucell_mass_gram, atom_density_cm3

        Character       Code*20

        Character       txt(20)*80
        common  /msg/   txt

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
          progname = 'Get_X0'             !|
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
        if (Code.ne.'Unknown')  Then  !--------------------+
c cccccc  isyng   = 1                  !cubic crystal      |
          isyng   = 0                  !no sym.control     |
          LayType = 0                  !amorph.crystal     |<-actually this is nrefl
          QB      = 0.                                    !|
          jpr     = 0                  !debug print        |
c         jpr = 2                                         !|
          Call X0h1 (a_lattice, Wave, Code, rho, LayType, !|
     *               0, 0, 0, isyng,                      !|
     *               mode_x0h, jpr, QB,                   !|
     *               Xr0, Xi0,                            !|
     *               Xrh, Xih,                            !|
     *               xqh, xrf, xif, xqf, ifail)           !|
          if (ifail.ne.0)  goto 100                       !|
c In the following, crystal is not the same:               |
          mode_x0h = 1                                    !|
        endif  !-------------------------------------------+

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

        j = 1
        do      i=2,N_Layers  !=====+
          j     = j + N_Step       !|
          X0(j) = X0(1)            !|
        enddo   !===================+

  100   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        Return
        End
