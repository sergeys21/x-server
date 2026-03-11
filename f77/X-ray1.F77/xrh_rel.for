        Real function Xrh_Rel (ad,wv,code,h,k,l,
     *                         ising,mode,ifail)
c -------------------------------------------------------
c                  Xrh_Rel = xrh /xr0
c -------------------------------------------------------
c       Input and output parameters:
c
c ad(6)  -      unit cell dimensions a,b,c in Angstrom and
c               angles alpha,beta,gamma in degrees (see X0H).
c wv     -      X-ray wavelength in Angstrom.
c code*20  -    codename of crystal structure or amorphous material.
c h,k,l -       Bragg reflections indices.
c ising -       crystal syngony index. On entry it provides control to
c               what syngony the crystal is expected to be (1-7). If
c               zero value is passed, control is off. On exit is shows
c               the structure syngony from "coord.x0h".
c mode     -    operation mode flag:
c             0 - calls X0h1.
c             1 - calls X0h1.
c             2 - calls X0h1.
c             3 - not first call, same crystal as before,
c                 but different X-ray wavelength passed with
c                 either "wave" or "radiat". Take all previously
c                 found crystal structure parameters and only
c                 read "atom.x0h" for scattering properties of
c                 atoms.
c             4 - same as mode=3, but ignore previous X-ray
c                 line name if non-zero wavelength is passed.
c             5 - not the first call; all parameters are the
c                 same, but different Bragg reflection. Do not
c                 read any data files.
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Complex         im, xrhc, xihc, sumexp
        Real*8          xhxh, xhxh_re, xhxh_im
        Real            ad(6), wv, x0_abs, xh_abs,
     *                  xr0, xi0, qb, rho,
     *                  xrh, xih, xqh,
     *                  xrf, xif, xqf,
     *                  wh(kcompMax),
     *                  fh(kcompMax),
     *                  cr, ci, sum,
     *                  s, d_spacing
        Integer         h, k, l, ising, mode, ifail,
     *                  inmax, ind(4), i, j, ipm, itrace,
     *                  iwarn, isico
        Character       code*20
c -------------------------------------------------------
c Common blocks for exchanging data with X0H:
        Character       radiat*6
        common  /x0pa1/ radiat

        Real            prcn(kcompMax)
        common  /x0pa3/ prcn

        Real            wave, df1(kcompMax), df2(kcompMax),
     *                  ab(kcompMax,nabcMax), s1, s21, s31, s32, s41,
     *                  e_radius, Compton, Rydberg, pi,
     *                  Planck_, Boltzman, Avogadro,
     *                  wedg(nedgMax), gedg(nedgMax)
        Integer         indx, ipr, nedg
        common  /x0pa4/ wave, df1, df2,
     *                  ab, s1, s21, s31, s32, s41,
     *                  e_radius, Compton, Rydberg, pi, indx, ipr,
     *                  Planck_, Boltzman, Avogadro,
     *                  wedg, gedg, nedg

        Integer         kcomp, kol(kcompMax)
        common  /x0pa6/ kcomp, kol

        Real            wx(kcompMax,kolMax),
     *                  wy(kcompMax,kolMax),
     *                  wz(kcompMax,kolMax)
        Integer         isym
        common  /x0pa7/ wx, wy, wz, isym

        Real            bdw(kcompMax),f0(kcompMax),
     *                  td(kcompMax),tq(kcompMax)
        Integer         z(kcompMax)
        common  /x0pa8/ bdw, f0, td, tq, z

        Real            a(6), vol
        Integer         lbev
        common  /back1/ a, vol, lbev

        Character       txt(20)*80
        common  /msg/   txt

        Real            Bragan, Forfac
        External        Bragan, Forfac
        Integer         IndCheck
        External        IndCheck
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------
c ipm=1 corresponds to convention by Pinsker (p.18,71). At this
c convention the Bloch waves are expanded as exp(-i*k*r).
c ipm=-1 corresponds to convention by Afanasiev and Baryshevsky.
c At this convention the Bloch waves are expanded as exp(i*k*r).
        Data    ipm     /-1/
        Data    im      /(0.,1.)/

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'X0h_Rel'            !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+

        Xrh_Rel = 0.                            !make GNU compiler happy
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
c--------------------------------------------------------
        pi       = 4.*atan(1.)          !the "pi" number
        e_radius = 2.81794e-05          !classical electron radius (Angstrom)
        wave  = wv
        ifail = 0
c -------------------------------------------------------
        if (radiat(1:1).eq.char(0))     radiat = ' '
        if (mode.eq.4 .and. abs(wave).lt.1.E-32)  radiat = ' '
c Is the structure code given?
        if (code(1:1).eq.char(0))       code   = ' '
c Disable results printout:
        ipr = 0
        rho = 0.

        if (mode.lt.3)  Then  !-----------------------+
          Call  X0h1 (ad, wv, code, rho, 1,          !|
     *                h, k, l,                       !|
     *                ising, mode, ipr,              !|
     *                qb,  xr0, xi0,                 !|
     *                xrh, xih, xqh,                 !|
     *                xrf, xif, xqf, ifail)          !|
          if (ifail.ne.0) goto 100 !------------------+--return--+
          x0_abs = xr0*xr0+xi0*xi0                   !|         !v
          if (x0_abs.gt.0.)  Then  !----------------+ |
            x0_abs = sqrt(x0_abs)                  !| |
c Find real and imaginary parts of x(h)*x(-h)       | |
c (see the Pinsker book, p.73):                     | |
            xhxh_re = xrh*xrh - xih*xih            !| |
            xhxh_im = 2.*xrh*xih*cos((xrf-xif)*pi) !| |
c |x(h)*x(-h)|^2:                                   | |
            xhxh    = xhxh_re*xhxh_re              !| |
     +              + xhxh_im*xhxh_im              !| |
            if (xhxh .gt. 0.0D0) Then  !----+       | |
c |x(h)*x(-h)|:                            !|       | |
              xhxh   = dsqrt (xhxh)        !|       | |
c sqrt(|x(h)*x(-h)|):                      !|       | |
              xh_abs = sngl (dsqrt(xhxh))  !|       | |
            else !--------------------------+       | |
              xh_abs = 0.0                 !|       | |
            endif  !------------------------+       | |
            Xrh_Rel = abs(xh_abs/x0_abs)           !| |
          else !------------------------------------+ |
            Xrh_Rel = 0.                           !| |
          endif  !----------------------------------+ |
          goto 100  !---------------------------------+--return--+
        endif  !--------------------------------------+         !v

c Crystal syngony index:
        if (ising.lt.0 .or. ising.gt.8) ising=0
c .......................................................
c Determine X-ray wavelength:

        itrace = 0
        if (mode.lt.5) Then  !+--------------+
          Call Redwav (wave,itrace,ifail)   !|
          if (ifail.ne.0)  goto 100         !|
        endif  !-----------------------------+
c .......................................................
        if (isym.ne.2) Then !--+
          inmax = 3           !|
        else  !----------------+
          inmax = 4           !|
        endif  !---------------+

        if (ising.ne.0 .and. ising.ne.isym) Then !---+
          Stop 'Xrh_Rel - Syngony mismatch error'   !|
        endif  !-------------------------------------+

        if (ising.eq.8) Then !--+
          Xrh_Rel = 0.         !|
          goto 100 !------------+-----------------------return---+
        endif  !----------------+                                v

c ------------------------------------------------------
c                    Calculate x0
c ......................................................
        cr  = -wave*wave*e_radius/(pi*vol)
        ci  = -ipm*wave/(2.*pi*vol)
        xr0 = 0.
        xi0 = 0.
        do      indx=1,kcomp  !=========================+
          xr0 = xr0 + cr*kol(indx)*f0(indx)*prcn(indx) !|
          xi0 = xi0 + ci*kol(indx)*td(indx)*prcn(indx) !|
        enddo  !========================================+
        x0_abs = xr0*xr0+xi0*xi0
        if (x0_abs.gt.0.)  Then !--+
          x0_abs = sqrt(x0_abs)   !|
        else !---------------------+
          Xrh_Rel = 0.            !|
          goto 100 !---------------+--------------------return---+
        endif  !-------------------+                             v

c ------------------------------------------------------
c            Data input for calculating xh:
c ......................................................
        ind(1)     = h
        ind(2)     = k
        ind(inmax) = l
        if (isym.eq.2)  ind(3)=-(ind(1)+ind(2))
        if (ind(1)    .eq.0  .and.
     *      ind(2)    .eq.0  .and.
     *      ind(inmax).eq.0) Then  !-----------------------+
            txt(1) = 'Reflection indices not specified!'  !|
            ifail  = 1                                    !|
            goto 100  !-----------------------return-------+---+
                                                          !|   |
        endif  !-------------------------------------------+   v

c Calculate Bragg angle. If QB>90, take QB=90 and warn:
        if (modebat.eq.0)  Then  !----------------+
c Take sin(QB)=1., QB=90. and warn:               |
          iwarn = 3                              !|
        else  !-----------------------------------+
c Take sin(QB)=0., QB=0. and warn:                |
          iwarn = 1                              !|
        endif !-----------------------------------+
c Skip calculating sin,cos,... (faster):
        isico = 0
        qb = Bragan (wave,ind,inmax,iwarn,d_spacing,isico,
     *                                s,s,s,s,s,s)
        if (abs(qb).lt.1.E-32) Then !--+
c The length of message in Bragan:     |
          ifail = 5                   !|
          goto 100 !-------------------+----------------return---+
        endif !------------------------+                         v
c -------------------------------------------------------
c Calculate Xh by summation over all crystal components:

c This is equivalent sin(qb)/lambda:
        s = 1./(2.*d_spacing)
        xrhc = (0.,0.)
        xihc = (0.,0.)
        do      indx=1,kcomp  !==========================+
          wh(indx) = exp(-bdw(indx)*s*s)                !|
          fh(indx) = Forfac (s,ifail)                   !|
          if (ifail.ne.0)       goto 100                !|
          sumexp   = (0.,0.)                            !|
          do    j=1,kol(indx)  !======================+  |
            sum = h*wx(indx,j)                       !|  |
     +          + k*wy(indx,j)                       !|  |
     +          + l*wz(indx,j)                       !|  |
            sumexp = sumexp + exp(2.*pi*im*sum*ipm)  !|  |
          enddo  !==============================^=====+  |
c                                               |        |
c  +--------------------------------------------+        |
c ipm=1 corresponds to convention by Pinsker (p.18,71).  |
c At this convention the Bloch waves are expanded as     |
c exp(-i*k*r).                                           |
c ipm=-1 corresponds to convention by Afanasiev and      |
c Baryshevsky. At this convention the Bloch waves are    |
c expanded as exp(i*k*r).                                |
c ..................................................     |
          sumexp = wh(indx) * prcn(indx) * sumexp       !|
          xrhc   = xrhc  + cr * fh(indx) * sumexp       !|
          xihc   = xihc  + ci * td(indx) * sumexp       !|
        enddo  !=========================================+

c Modules of xrh,xih,xqh:
        xrh = abs(xrhc)
        xih = abs(xihc)
c Phases of xrh,xih:
        if (xrh.gt.0.) Then !---------------------+
          xrf = Atan2(aimag(xrhc),real(xrhc))/pi !|
        else  !-----------------------------------+
          xrf = 0.                               !|
        endif !-----------------------------------+
        if (xih.gt.0.) Then !---------------------+
          xif = Atan2(aimag(xihc),real(xihc))/pi !|
        else  !-----------------------------------+
          xif = 0.                               !|
        endif  !----------------------------------+

c Determine real and imaginary parts of x(h)*x(-h)
c (see the Pinsker book, p.73):
        xhxh_re = xrh*xrh - xih*xih
        xhxh_im = 2.*xrh*xih*cos((xrf-xif)*pi)
c |x(h)*x(-h)|^2:
        xhxh    = xhxh_re*xhxh_re
     +          + xhxh_im*xhxh_im
        if (xhxh .gt. 0.0D0) Then !-------+
c |x(h)*x(-h)|:                          !|
          xhxh   = dsqrt (xhxh)          !|
c sqrt(|x(h)*x(-h)|):                    !|
          xh_abs  = sngl (dsqrt(xhxh))   !|
          Xrh_Rel = abs(xh_abs/x0_abs)   !|
        else !----------------------------+
          Xrh_Rel = 0.                   !|
        endif  !--------------------------+

                                                                !v
  100   continue  !----------------------------------------------+
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
