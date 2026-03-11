        Subroutine SigmaPi (isipi, QB, Wave, x0_abs,    !inputs
     *                      xrh, xihD, xihQ,            !inputs
     *                      xrf, xifD, xifQ,            !inputs
     *                      xrhPol, xihPol, xdfPol,     !OUTputs
     *                      xhPol,                      !OUTputs, added 2001/10
     *                      extiPol, halfPol, prcPol)   !OUTputs
c--------------------------------------------------------
c Calculate crystal susceptibilities for sigma- and pi-
c polarized X-ray using data provided by subroutine x0h1.
c
c isipi - polarization flag (isipi=1 - sigma,  isipi=2 - pi)    /Input/
c QB    - Bragg angle in degrees.                               /Input from x0h/
c Wave  - X-ray wavelength in Angstroms.                        /Input/
c x0_abs - 0-th crystal polarizability.                         /Input from x0h/
c xrh   - the module of the real part of the H-th crystal       /Input from x0h/
c         polarizability (xrh_module).
c xihD,xihQ - the modula of the dipole and quadrupole
c         components of the imaginary part of the H-th crystal  /Input from x0h/
c         polarizability (xihD_module, xihQ_module).
c xrf,xifD,xifQ - the phases of the above parameters            /Input from x0h/
c         (xrh_phase,xihD_phase,  xihQ_phase)
c
c xrhPol - the calculated module of the real part of the        /Output/
c         H-th crystal polarizability for sigma- or pi-
c         polarization (depending on isipi=1 or isipi=2)
c xrhPol - the calculated module of the imaginary part of       /Output/
c         the H-th crystal polarizability for sigma- or pi-
c         polarization (depending on isipi=1 or isipi=2)
c xdfPol - the phase difference between the two parameters      /Output/
c         above.
c xhPol  - sqrt[x(h)*x(-h)]  -- added 2001/10                   /Output/
c extiPol - the calculated extinction length of x-rays in       /Output/
c         crystal for sigma- or pi-  polarization (in
c         micrones)
c halfPol - the calculated halfwidth of the symmetric Bragg     /Output/
c         reflection for sigma- or pi-  polarization (in arc
c         seconds)
c prcPol - the calculated relative strength of Bragg            /Output/
c         reflection 'H' with respect to the reflection '0'
c         for sigma- or pi- polarization (in percents).
c--------------------------------------------------------
        Complex xic, xdc, xqc
        Real*8  xhxh, xhxh_re, xhxh_im
        Real    QB, QB_, Wave, x0_abs,
     *          xrh, xihD, xihQ,
     *          xrf, xifD, xifQ,
     *          xrhPol,  xihPol,  xdfPol,
     *          extiPol, halfPol, prcPol,
     *          xifPol, halfPol_max, xhPol,
     *          pi, gra, sec, c2, c4, s2
        Integer isipi

        pi  = 4.*atan(1.)               !the "pi" number
        gra = 2.*pi/360.                !number of radians in 1 degree
        sec = gra/3600.                 !number of radians in 1 arc sec
        QB_ = QB
        if (QB_ .ge. 90.)  QB_ = 90.-1.e-05
        xdc = xihD*cmplx(cos(xifD*pi),sin(xifD*pi))
        xqc = xihQ*cmplx(cos(xifQ*pi),sin(xifQ*pi))
        c2  = abs(cos(2.*QB_*gra))
c - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (isipi .eq. 1) Then  !--------+Sigma-polarization
          xrhPol = xrh                  !|
          xic    = xdc+c2*xqc           !|
        else  !--------------------------|Pi-polarization
          c4     = abs(cos(4.*QB_*gra)) !|
          xrhPol = c2*xrh               !|
          xic    = c2*xdc+c4*xqc        !|
        endif  !-------------------------+
c - - - - - - - - - - - - - - - - - - - - - - - - - - -
        xihPol = abs(xic)
        if (abs(xihPol).lt.1.E-32) Then  !---------+
          xifPol = 0.                             !|
        else  !------------------------------------|
          xifPol = atan2(aimag(xic),real(xic))/pi !|
        endif  !-----------------------------------+
c Phase difference: 
        xdfPol = xrf-xifPol
c Round errors:
        if     (abs(xdfPol)    .lt. 0.01)  Then  !---+
          xdfPol = 0.                               !|
        elseif (abs(xdfPol-1.) .lt. 0.01)  Then  !---|
          xdfPol = 1.                               !|
        elseif (abs(xdfPol+1.) .lt. 0.01)  Then  !---|
          xdfPol =-1.                               !|
        endif  !-------------------------------------+
c Find real and imaginary parts of x(h)*x(-h)
c (see the Pinsker book, p.73):
        xhxh_re = xrhPol*xrhPol - xihPol*xihPol
        xhxh_im = 2.*xrhPol*xihPol*cos(xdfPol*pi)
c |x(h)*x(-h)|^2:
        xhxh    = xhxh_re*xhxh_re + xhxh_im*xhxh_im
        if (xhxh .gt. 0.0D0) Then  !----+
c |x(h)*x(-h)|:                        !|
          xhxh  = dsqrt (xhxh)         !|
c sqrt(|x(h)*x(-h)|):                  !|
          xhPol = Sngl(dsqrt(xhxh))    !|
        else !--------------------------|
          xhPol = 0.0                  !|
        endif  !------------------------+
c Bragg peak FWHM:
        if (xhPol .gt. 0.) Then  !-------------------------+
          s2 = sin(2.*QB_*gra)                            !|
          if (abs(s2).lt.(1.E-20)) s2 = 1.E-20            !|
c Usual halfwidth:                                         |
          halfPol = 2.*xhPol/(s2*sec)                     !|
c The estimate for exact back diffraction:                 |
          halfPol_max = sqrt(2.*xhPol)/sec                !|
          if (halfPol .gt. halfPol_max) Then !--+          |
            halfPol = halfPol_max              !|          |
          endif  !------------------------------+          |
        else  !--------------------------------------------|
          halfPol = 0.                                    !|
        endif  !-------------------------------------------+
c The extinction length (Bragg case!):
        if (abs(xhPol).gt.1.E-32) Then !-------------------+
          extiPol = 0.0001*Wave*sin(QB_*gra)/(pi*xhPol)   !|
        else !---------------------------------------------|
          extiPol = 0.                                    !|
        endif  !-------------------------------------------+
c Relative intensity (%):
        prcPol = 100.*xhPol/abs(x0_abs)
        return
        end

