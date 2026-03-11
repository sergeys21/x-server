        Subroutine Get_Xh97 (
     *    Wave,         !input, but may be output: if wave=0, then expect QB#0 and get wave from QB
     *    Code,         !input
     *    rho,          !output: density (g/cm^3 -- for amorphous structures only)
     *    iSyngony,     !output
     *    LayType,      !input: 0/1=number of reflections to calculate
     *    indices,      !input, array of 3: only used if LayType=1
     *    QB,           !input/output: if QB=0, calc. it; if QB#0, compare with prev.
     *    mode_x0h,     !input: 0-5 (normally 1st call with 0=md_NewAll, then 1=md_newCrystal)
     *    Xr0,Xi0,      !output
     *    Xrh,Xih,Xdf,  !output
     *    X0,XH,        !output
     *    ifail)        !output: number of lines in error/warning message (ifail<0 means warning)
c--------------------------------------------------------
c This subroutine is called from gid_sl97 and gids_97
c--------------------------------------------------------
        Complex         X0, XH
        Real            Wave, rho, QB, QB_,
     *                  Xr0, Xi0, Xrh, Xih, Xdf,
     *                  xqh, xrf, xif, xqf, a_lattice(6)
        Integer         iSyngony, LayType,
     *                  indices(3), mode_x0h,
     *                  ifail, jpr, iwait, i,
     *                  lcod
        Character       Code*20
        Integer                 iHenkeCowan
        Common  /HenkeCowan/    iHenkeCowan

        Character       txt(20)*80
        Common  /msg/   txt

c New of August'99:
c       Parameter       (kcompMax = 10)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Real            prcn(kcompMax)
        Common  /x0pa3/ prcn
c This is an interface to internal X0h data.
        Real            ucell_mass_gram, atom_density_cm3
        Integer         n_atoms_ucell
        Common  /x0paA/ ucell_mass_gram, atom_density_cm3, n_atoms_ucell

        Real*8          DVecMod, d8, gra8
        External        DVecMod
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
          progname = 'Get_Xh97'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
c Prevents GNU fortran -Wextra warnings:
        i = kcompMax
        i = kolMax
        i = ncodMax
        i = natmMax
        i = nwavMax
        i = nedgMax
        i = nabcMax
        i = nedgTb
        i = maxX0hWarnLines
c -------------------------------------------------------
        ifail    = 0
        Xr0 = 0.
        Xi0 = 0.
        Xrh = 0.
        Xih = 0.
        Xdf = 0.
        xqh = 0.
        xrf = 0.
        xif = 0.
        xqf = 0.
c -------------------------------------------------------
c Possible values of mode_x0h (see x0h1.for)
c       md_NewBase     =-1
c       md_NewAll      = 0
c       md_NewCrystal  = 1      !use new name,prcn,a(1-6) from X0h DB
c       md_NewCrystal_a= 2      !use new name,prcn from X0h DB, use ad(1-6) if non-zero, otherwise from X0h DB
c       md_NewWave     = 3
c       md_NewLine     = 4
c       md_NewHKL      = 5
c -------------------------------------------------------
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

        if (Wave.le.0.) Then !-----------------------------------------+
          if (Code.eq.'Unknown' .OR.                                  !|
     *        LayType.eq.0      .OR.                                  !|
     *        abs(QB).lt.1.E-32 .OR.                                  !|
     *        (indices(1).eq.0 .AND.                                  !|
     *         indices(2).eq.0 .AND.                                  !|
     *         indices(3).eq.0)) Then !------------------------------+ |
            lcod = Max(Len_Trim(Code),1)                            !| |
            write (txt,1) code(1:lcod),LayType,QB,indices           !| |
  1         format('X0H -- E R R O R:'//
     *      'Cannot determine X-ray wavelength under the conditions:'/
     *      'structure=',a,' type=',i1,' Bragg angle=',f5.2,
     *                                    ' reflection=(',3i3,')')  !| |
            ifail = 4                                               !| |
            goto 100                                                !| |
          else !-----------------------------------------------------+ |
            do i=1,kcompMax !==+                                    !| |
              prcn(i) = 0.    !|                                    !| |
            enddo  !===========+                                    !| |
            jpr = 0                                                 !| |
            Wave = 0.1                  !the shortest we can use    !| |
            Call X0h1  (a_lattice, Wave, Code, rho, LayType,        !| |
     *                  indices(1), indices(2),                     !| |
     *                  indices(3), iSyngony,                       !| |
     *                  mode_x0h, jpr, QB_,                         !| |
     *                  Xr0, Xi0,                                   !| |
     *                  Xrh, Xih,                                   !| |
     *                  xqh, xrf, xif, xqf, ifail)                  !| |
            if (ifail.ne.0)  goto 100                               !| |
            gra8 = Datan(1.0D0)/45.                                 !| |
            d8 = 1./DVecMod(indices,3)                              !| |
            Wave = sngl (2.*d8*Dsin(QB*gra8)) !wave is given by QB   | |this is when Wave=0
            QB = 0.                                                 !| |
            QB_= 0.                                                 !| |
          endif !----------------------------------------------------+ |
        endif !--------------------------------------------------------+

        if (Code.ne.'Unknown')  Then  !---------------------+
          do i=1,kcompMax !==+                             !|
            prcn(i) = 0.    !|                             !|
          enddo  !===========+                             !|
          jpr = 0                                          !|
c         jpr = 2                                          !|
          QB_ = 0.0                                        !|
          Call X0h1  (a_lattice, Wave, Code, rho, LayType, !|
     *                indices(1), indices(2), indices(3),  !|
     *                iSyngony, mode_x0h, jpr, QB_,        !|
     *                Xr0, Xi0,                            !|
     *                Xrh, Xih,                            !|
     *                xqh, xrf, xif, xqf, ifail)           !|
          if (ifail.ne.0)  goto 100                        !|

c This is a "normal" sign of the phase difference,          |
c but it is opposite to what is shown by X0h:               |
          Xdf = xif - xrf                                  !|
c Therefore, we invert it!!!                                |
          Xdf = -Xdf                                       !|
c   +=============================+                         |
c   |  ATTENTION: the quadrupole  |                         |
c   | contributions are ignored!  |                         |
c   +=============================+                         |
c At next call the crystal might be not the same:           |
          mode_x0h = 1                                     !|
        endif  !--------------------------------------------+

c What to do at absorption edges
c where xr0 may change the sign?
c This however should not be a
c problem as long as we use x0h
c since x0h does not process the
c edges!
c       X0 = Cmplx(-abs(Xr0),abs(Xi0))
c Changed 2003/07/25:
c We now have control over manual input in Xinput97 -> ChiSearch
c and X0h should calculate correctly!
c Hope it is safe!
        X0 = Cmplx(Xr0,abs(Xi0))

        if (LayType.eq.1)  Then  !----------------------+
c This is an old definition corresponding to phase=pi:  |
c         XH = Cmplx(-abs(Xrh),abs(Xih))               !|
c Starting from June-98 we choose:                      |
          XH = Cmplx(abs(Xrh),abs(Xih))                !|
c -- then, the real Xh and Xh_ are:                     |
c         Xh = Xrh+Xih*sin(Xdf*pi) + i*Xih*cos(Xdf*pi)  |
c         Xh_= Xrh-Xih*sin(Xdf*pi) + i*Xih*cos(Xdf*pi)  |
c -- where the sign before sin corresponds to our       |
c "inverted" Xdf=Xrf-Xif  !!!!!!!!!!!!!!!!!!!!!         |
          if (abs(QB).lt.1.E-32) Then  !--------------+ |
            QB = QB_                                 !| |
          else  !-------------------------------------| |
            if (abs(QB-QB_).gt.1.)  Then  !---------+ | |
              txt(1)=progname//' WARNING: Bragg'// !| | |
     *               ' angle mismatch for '//Code  !| | |
              if (modebat.ne.0) Then !--+           | | |
                iwait=0                !|           | | |
              else !--------------------|           | | |
                iwait=2                !|           | | |
              endif  !------------------+           | | |
              Call Message(txt,1,iwait)            !| | |
            endif  !--------------------------------+ | |
          endif  !------------------------------------+ |
        else  !-----------------------------------------|
          XH = (0.,0.)                                 !|
        endif  !----------------------------------------+

  100   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        Return
        End
