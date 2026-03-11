        Subroutine  X0h2 (ad,wv,code,rho,nref,index,inciden,
     *                    ising,mode_x0h,jpr,QB,xr0,xi0,ndim,
     *                    xrh,xih,xqh,xrf,xif,xqf,
     *                    ifail)
c -------------------------------------------------------
c ATTENTION: additional data is transferred via these common blocks:
c       common  /x0pa1/ radiat
c       common  /x0pa2/ name
c       common  /x0pa3/ prcn(kcompMax)
c       common  /lunds/ luwa,luat,lusp,luco
c       common  /x0hpa/ x0hpath,lx0hpat
c       Common  /HenkeCowan/ iHenkeCowan
c       common  /batmode/ modebat,lunbat,ierrbat,
c    *                    progname, progrezv, istackrezv, ErrorFile
c       common  /x0pa9/ Poisson, edge_nearest, name_nearest
c -------------------------------------------------------
c       Input and output parameters:
c
c ad(6)  -      unit cell dimensions a,b,c in Angstrom and
c               angles alpha,beta,gamma in degrees.
c              -for cubic structures one needs to pass only a(1)=a,
c              -for hexagonal/trigonal and tetragonal: a(1)=a, a(3)=c,
c              -for trigonal(rhombohedral): a(1)=a, a(4)=alpha,
c              -for orthorhombic: a(1)=a, a(2)=b, a(3)=c,
c              -for monoclinic:   a(1)=a, a(2)=b, a(3)=c, a(5)=beta,
c              -for triclinic:    a(1)=a, a(2)=b, a(3)=c, a(4)=alpha,
c               a(5)=beta, a(6)=gamma.
c               If a(i) not given (zeroes are passed), then X0H tries
c               to read them from the "coord.x0h" file and if not found,
c               then requests from terminal (if in interactive mode).
c               At exit all a(i) determined.
c wv     -      X-ray wavelength in Angstrom. If zero wavelength is
c               passed and also if the name of characteristic X-ray
c               line "radiat" is empty, then X0H requests from terminal
c               first the name "radiat" and if not given, then the
c               wavelength. At exit the wavelength is determined.
c code*20  -    codename of crystal structure or amorphous material.
c               If empty code is passed, it is requested from the
c               terminal (in the interactive mode). At exit the codename
c               is determined. Codename must be either a structure name
c               in "coord.x0h" or atom name in "atom.x0h" (then only X0
c               is calculated), or a chemical formula consisting of
c               atoms names from "atom.x0h with their integer repetition
c               factors. In the last case only X0 is calculated and X0H
c               needs specification of material density rho.
c rho    -      material density in gm/cm^3. Used for amorphous
c               materials only. If zero value is passed, X0H tries
c               to find it in "coord.x0h" or "atom.x0h" and if
c               unsuccessful, then X0H asks from terminal (in the
c               interactive mode). If non-zero rho is passed to X0H,
c               it overrides the rho data from "coord.x0h".
c nref   -      number of reflections for which xrh,xih are requested.
c               It must be >= 1 and does not change on exit.
c index(3,nref)-Bragg reflections indices. If zero hkl are passed for a
c               reflection, then X0H requests them from terminal.
c               One exception are index(1:3,inciden) which must be
c               (0,0,0) because they are for zero reflection (incident
c               or specularly reflected wave). For hexagonal and
c               trigonal structures the indices need to be in orthogonal
c               settings. At exit all indices are defined.
c inciden  -    position of 0th reflection (incident wave) in the Bragg
c               indices array (1:nref). If inciden=0, it is taken to be
c               the first array element (inciden=1).
c ising -       crystal syngony index. On entry it provides control to
c               what syngony the crystal is expected to be (1-8). If
c               zero value is passed, control is off. On exit is shows
c               the structure syngony from "coord.x0h". The index may
c               have the following values:
c               1 - cubic syngony,
c               2 - hexagonal/trigonal syngony,
c               3 - tetragonal syngony,
c               4 - trigonal (rhombohedral) syngony,
c               5 - orthorhombic syngony,
c               6 - monoclinic syngony,
c               7 - triclinic syngony,
c               8 - **amorphous**
c mode_x0h -    operation mode flag:
c            -1 - same as 0, but re-read the database index             !md_NewBase
c                 (for X0hWin / Sitnikov)
c                 ATTENTION: on exit X0H overwrites this parameter
c                 value by mode=0. --- DON'T USE CONSTANTS!!!
c             0 - first call; take, if present, ALL                     !md_NewAll
c                 input parameters: ad,wv,code,nref,h,k,l,
c                 and common area parameters:
c                 common /x0pa1/ radiat
c                 common /x0pa2/ name
c                 common /x0pa3/ prcn(kcompMax)
c             1 - not first call, but different structure.              !md_NewCrystal
c                 Take, if present, input parameters
c                 wv,code,nref,h,k,l, as well as these
c                 common values:
c                 common /x0pa1/ radiat
c                 and clear (set to zero) all lattice parameters
c                 <ad> and the parameters from these commons:
c                 common /x0pa2/ name
c                 common /x0pa3/ prcn(kcompMax)
c                 common  /x0pa9/ Poisson, edge_nearest,
c                                          name_nearest
c             2 - same as mode=1, but do not clear lattice              !md_NewCrystal_a
c                 parameters <ad>. If there are non-zero <ad>
c                 in the X0H1 calling command, then take them;
c                 otherwise use the previously found.
c             3 - not first call, same crystal as before,               !md_NewWave
c                 but different X-ray wavelength passed with
c                 either "wave" or "radiat". Take all previously
c                 found crystal structure parameters and only
c                 read "atom.x0h" for scattering properties of
c                 atoms.
c             4 - same as mode=3, but ignore previous X-ray             !md_NewLine
c                 line name if non-zero wavelength is passed.
c             5 - not the first call; all parameters are the            !md_NewHKL
c                 same, but different Bragg reflection. Do not
c                 read any data files.
c jpr      -    results printing level.
c             0 - printing is disabled,
c             1 - print results only,
c             2 - print results, data from X0H files, and
c                 used calculations schemes.
c QB(nref) -    Bragg angles of reflections in degrees (output only).
c               These are calculated in X0H and returned to the caller.
c xr0,xi0 -     real and imaginary parts of the zeroth term of the
c               crystal dielectric susceptibility Fourier expansion
c               (output only). These are calculated in X0H returned
c               to the caller.
c ndim  -       dimension of square matrices xrh,xih,xqh,xrf,xif,xqf.
c               Only nref*nref block is filled in the matrices.
c               Therefore ndim >= nref must be satisfied.
c xrh(nref*2),
c xih(nref*2),
c xqh(nref*2),
c xrf(nref*2),
c xif(nref*2),
c xqf(nref*2) - modula of the real and imaginary parts of dipole
c               component and the imaginary part of the quadruple
c               component of the crystal dielectric susceptibility
c               Fourier expansion for reflection (hkl), as well as
c               their phases in pi-inits (output only). These are
c               calculated in X0H returned to the caller.
c .......................................................
c       Optional parameters passed via common blocks
c           /x0pa1/, /x0pa2/, /x0pa3/, /x0pa9/:
c
c       common  /x0pa1/ radiat
c       common  /x0pa2/ name
c       common  /x0pa3/ prcn()
c       common  /x0pa9/ Poisson, edge_nearest, name_nearest
c
c radiat*6   -  name of characteristic X-ray line. It is used only
c               when the X-ray wavelength "wv" passed zero value.
c               If wv=0 and radiat is empty, X0H requests from
c               terminal first radiat, and if empty, then wv.
c               On exit radiat is filled if was specified on entry
c               or entered from terminal (in the interactive mode).
c name(kcompMax)*4 - names of chemical element composing the crystal
c               as contained in the "atom.x0h" file. If the names are
c               empty, X0H tries to read them from the crystal
c               description in the "coord.x0h" file or if not present,
c               then asks from terminal (in the interactive mode). On
c               exit all names are defined.
c prcn(kcompMax) - percentage (share) or probability of atoms filling
c               their nodes in case of solid solutions. if non-zero
c               values are passed, they override the values from
c               "coord.x0h". If values are beyond the [0-1] range,
c               then X0H requests them from terminal (in the
c               interactive mode).  On exit all shares are defined.
c Poisson    -  Poisson's ratio (for cubic structures only). When
c               this is present for given crystal in "coord.x0h", it
c               is returned to caller. Otherwise, X0H returns zero.
c edge_nearest- energy in KeV of the absorption edge in the crystal
c               structure which closest to the energy of the incident
c               x-rays (output parameter).
c name_nearest- name of the atom in the crystal which has this
c               absorption edge (output parameter).
c .......................................................
c       Optional parameter passed via common /lunds/:
c
c       common  /lunds/ luwa,luat,lusp,luco
c
c luwa,luat,
c lusp,luco  -  Logical Unit Numbers (LUN) of the X0H database files.
c               They must differ from each other (X0H does not control
c               it). If zeroes are passed, X0H uses LUNs 11-14.
c .......................................................
c       Optional parameters passed via common /X0hpa/
c
c       common  /x0hpa/ x0hpath,lx0hpat
c
c x0hpath*80  -  path to directoty with X0H database files.
c lx0hpat     -  number of characters in the path.
c .......................................................
c       Optional parameters returned by X0H via common /x0hwarn/:
c
c       common  /x0hwarn/       txtX0hWarn, linesX0hWarn
c
c txtX0hWarn(maxX0hWarnLines)*80  -  warning messages array.
c linesX0hWarn                    -  number of lines recorded into
c                                    it by X0h
c -------------------------------------------------------
c                       Subroutines:
c   Redcor  --x0hread--         Lengh   --serve.lib--
c   Redlat                      Defpat
c   Redwav                      Pause
c   redpar
c   Forfac  --x0hmake--
c   Disint
c   Sort
c   Gaus
c   Disper
c   Parratt
c   Crosec
c   Fifunc
c               Algebra.for + Algebra3.for
c               Bragan8.for
c               ChiParse.for (or ChiPars2.for)
c               indCheck.for
c -------------------------------------------------------
c This program calculates complex Fourier components x0,xh
c in the X-ray energy range for crystals which belong to any
c syngony. The atomic coordinates are taken by the structure
c codename from the "coord.x0h" file.
c
c          Copyright: S.Stepanov, 1985-1998 (*)
c (*) Some early code (the temperature dependencies of xh and
c the interpolation of dispersion corrections) was written
c in collaboration with Olga Lugovskaya.
c -------------------------------------------------------
c               General calculations formulae:
c .......................................................
c              2
c         -wave *e_radius    n
c xr0  =  --------------- * sum(fh (0))
c              pi*vol       j=1   j
c .......................................................
c          -wave       n
c xi0  =  --------- * sum(tfull )
c         2*pi*vol    j=1      j
c .......................................................
c              2
c         -wave *e_radius    n          sin(QB)
c xrh  =  --------------- * sum(wh *fh (-------)*exp(i*h*r ))
c              pi*vol       j=1   j   j   wave            j
c .......................................................
c    d,q   -wave      n       d,q
c xih   = -------- * sum(wh *t    *exp(i*h*r ))
c         2*pi*vol   j=1   j  j             j
c -------------------------------------------------------
c    sg                        pi
c xrh  = xrh ,              xrh  = xrh * /cos(2*QB)/ ,
c             sg     d                  q
c          xih  = xih             +  xih * /cos(2*QB)/ ,
c             pi     d                  q
c          xih  = xih * /cos(2*QB)/ +  xih * /cos(4*QB)/
c
c n        - number of atoms in unit cell,
c vol      - unit cell volume ,
c h=2*pi/d - reciprocal lattice vector corresponding to hkl reflection,
c d        - interplanar spacing corresponding to hkl reflection.
c
c wh       - Debye-Waller temperature factor for hkl.
c -------------------------------------------------------
c As follows from the above equations, to calculate x0,xh one needs
c to know the following data for atoms composing the crystal:
c - structure amplitudes fh(s) with the account for dispersion
c   corrections
c - dipole, quadrupole, and dipole-octopole absorption cross sections
c   td,tq,to,
c - Debye-Waller temperatures of atoms.
c
c These data are contained in the "atom.x0h" file and can be found by
c the names of atoms composing crystal. To compensate for incompleteness
c of existing tables, each parameter can be calculated by several
c methods:
c
c 1. Calculation of atomic scattering factors fh which lead to
c determining xrh is carried out in subroutine forfac and based on
c the interpolation formula and tables from [1-3]. The fh value depends
c on the atomic number (number of electrons in the atom) and parameter
c s=sin(QB)/wave (wave is the X-ray wavelength). The interpolation
c coefficients are read from "atom.x0h" by given atom name. If these
c coefficients are null (not known) then the second method for
c calculating fh is applied: we take fh(0)=z, (z is the number of
c electrons in the atom), and for calculating fh(s) X0H requests two
c closest s,fh(s) from the tables and applies linear interpolation.
c
c 2. Dispersion corrections are calculated in subroutine disper and
c based on the method described in [4]. To carry out the calculations,
c X0H read from "atom.x0h" the energies of absorption edges in KeV and
c respective oscillator strengths. The energies in "atom.x0h" were
c mostly taken from [5] and the oscillator strengths from [4]. If the
c oscillator strengths are unknown for an atom, then X0H applies the
c second method. Instead of oscillator strengths the atom record in
c "atom.x0h" is filled with tabulated dispersion corrections for 5
c fixed wavelengths: Cr-Ka1, Fe-Ka1, Cu-Ka1, Mo-Ka1, Ag-Ka1, which
c can be taken e.g. from [4,6-7]. The the oscillator strengths are
c determined by applying equations from [4] to calculate known
c tabulated dispersion corrections. This is implemented in subroutines
c disint-sort-gaus.
c
c 3. The calculation of atomic absorption cross sections td and tq
c needed for xi0 and xih is carried out in the subroutine crosec and
c based on the method described in [8]. The screening constants for
c electronic shells are taken from [9]. Two terms are  calculated: td
c (the sum of dipole and octopole contributions) and tq (the quadrupole
c contribution). The screening constants are read from "atom.x0h".
c This method is recommended for atoms with Z in range 6-54. If the
c screening constants are not available of f Z is beyond the 6-54
c range, then there 4 more methods of calculating the dipole cross
c section contribution listed below (in this case tq is taken zero).
c - Method-2: calculate td using formula
c   mu/rho = c * l**3 - d * l**4,
c   td = (mu/rho)*(a/Avogadro),
c   where
c   l is x-ray wavelength, a is the atomic weight, Avogadro is the
c   Avogadro constant, and the coefficients c,d can be taken from [2].
c - Method-3: same as Method-2, but c,d are not known. For determining
c   c,d file "atom.x0h" provides mu/rho for 5 characteristic X-ray
c   lines.
c - Method-4: same as Method-3, but mu/rho in "atom.x0h" are given
c   for 5 specified X-ray wavelengths.
c   The values of mu/rho for methods 3 and 4 are taken from [2,5,6].
c - Method-5: td is calculated using formula
c   td = 2 * l * e_radius * df2,
c   where e_radius is classical radius of electron and df2 is the
c   imaginary part of dispersion correction calculated as described
c   above (see 2.).
c
c 4. The Debye-Waller temperature factors (the Debye coefficients in
c the formula exp(-b*s*s)) are either directly read by atom name from
c "atom.x0h", in which they are taken from [2,9], or they are calculated
c in subroutine fifunc using formula from [10] and the Debye temperature
c of atom from "atom.x0h" and the crystal temperature. The Debye
c temperatures in "atom.x0h" are mostly taken from [10].
c -------------------------------------------------------
c                      REFERENCES
c  [1] P.A.Doyle,P.S.Turner - Acta Crystallogr.(1968),
c      v.A24, p.390-397.
c  [2] International Tables for X-ray Crystallography, v.4,
c      - England, Kynoch press, 1974, (p.99-101,232-246)
c  [3] D.T.Cromer, J.T.Waber - Acta Crystallogr.(1965),
c      v.18, p.104-109.
c  [4] D.T.Cromer - Acta Crystallogr.(1965), v.18, p.17-23.
c  [5] M.A.Blokhin, I.G.Schweitzer - Rentgenospektralnyj spravochnik,
c      - Moscow, "Nauka", 1982 (tables 6-7).
c  [6] D.T.Cromer, D.Liberman - J.Chem.Phys.(1970), v.53,
c      p.1891-1898.
c  [7] H.Wagenfeld - Physical Review,(1966), v.144,
c      p.216-224.
c  [8] G.Hildebrandt, J.D.Stephenson, H.Wagenfeld -
c      Z.Naturforschung,(1975), v.30a, p.697-707.
c  [9] Z.G.Pinsker - Rentgenovskaja kristallooptika, -
c      Moscow, "Nauka", 1982.
c [10] A.A.Rusakov - Rentgenografia metallov, - Moscow,
c      "Atomizdat", 1977, (p.134-137).
c-------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c       Parameter       (maxX0hWarnLines = 10);
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Integer         ndim
        Complex         im,xrhc,xihd,xihq,sumexp
        Real*8          QBij, d8
        Real            ad(6),energy,wv,QB(ndim),
     *                  xr0,xi0,xi0d,xi0q,
     *                  xrh(ndim,ndim),xrf(ndim,ndim),
     *                  xih(ndim,ndim),xif(ndim,ndim),
     *                  xqh(ndim,ndim),xqf(ndim,ndim),
     *                  wh(kcompMax),fh(kcompMax),
     *                  cr,ci,gra,sec,mu0,a0,sum,
     *                  rho, rho_i, r_atoms_ucell,
     *                  TER, dl_TER, delta, eta, d, s
        Integer         index(3,ndim),inciden,incid,jndx,
     *                  nref,ind(4),need_a(6),ir,
     *                  i,j,jpr,ipm,mode_x0h,mode,ising,ifail,itrace,
     *                  iimode,inmax,lcod,lbev,md_NewHKL,
     *                  md_NewBase,md_NewAll,md_NewCrystal,
     *                  md_NewCrystal_a,md_NewWave,md_NewLine,
     *                  n,m,ifa,ll,iref,iwarn,isico,nb
        Character       code*20, smtry(8)*24

        Real*8          Bragan8
        External        Bragan8
        Real            Forfac, wave2energy
        External        Forfac, wave2energy
        Integer         indCheck, CodeType
        External        indCheck, CodeType
c -------------------------------------------------------
c Common blocks for passing parameters from calling programs:
        Character       radiat*6
        Common  /x0pa1/ radiat
        Character       name(kcompMax)*4
        Common  /x0pa2/ name
        Real            prcn(kcompMax)
        Common  /x0pa3/ prcn
        Real            Poisson, edge_nearest
        Character       name_nearest*4
        Common  /x0pa9/ Poisson, edge_nearest, name_nearest
c This is not passed anywhere, just to store between x0h2 calls:
        Real            atom_mass_aem(kcompMax)
        Common  /x0paS/ atom_mass_aem
        Integer         luwa,luat,lusp,luco
        Common  /lunds/ luwa,luat,lusp,luco
        Integer         linesX0hWarn
        Character       txtX0hWarn(maxX0hWarnLines)*80
        Common /x0hwarn/ txtX0hWarn, linesX0hWarn
c
c txtX0hWarn(maxX0hWarnLines)*80  -  warning messages array.
c linesX0hWarn                    -  number of lines recorded into it by X0h
c -------------------------------------------------------
c Common blocks for passing parameters between X0H subroutines:
        Real            wave,df1(kcompMax),df2(kcompMax),
     *                  ab(kcompMax,nabcMax),s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,
     *                  Planck_,Boltzman,Avogadro,
     *                  wedg(nedgMax),gedg(nedgMax)
        Integer         indx,ipr,nedg
        Common  /x0pa4/ wave,df1,df2,ab,s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,indx,ipr,
     *                  Planck_,Boltzman,Avogadro,wedg,gedg,nedg
        Real            ddeb(kcompMax,2)
        Integer         mdeb(kcompMax)
        Common  /x0pa5/ ddeb,mdeb
        Integer         kcomp,kol(kcompMax)
        Common  /x0pa6/ kcomp,kol
        Real            wx(kcompMax,kolMax),
     *                  wy(kcompMax,kolMax),
     *                  wz(kcompMax,kolMax)
        Integer         isym
        Common  /x0pa7/ wx,wy,wz,isym
        Real            bdw(kcompMax),f0(kcompMax),
     *                  td(kcompMax),tq(kcompMax)
        Integer         z(kcompMax)
        Common  /x0pa8/ bdw,f0,td,tq,z
        Real            ucell_mass_gram, atom_density_cm3
        Integer         n_atoms_ucell
        Common  /x0paA/ ucell_mass_gram, atom_density_cm3, n_atoms_ucell
c X0H is can handle crystals with the number of components not
c exceeding kcompMax and the number of atomic position for
c each component not exceeding kolMax.
        Integer         natm,katm,ismatm(natmMax),idatm,
     *                  ncod,kcod,kelem,ismcod(ncodMax),idcod,
     *                  nwav,kwav,ismwav(nwavMax),idwav
        Character       lstatm(natmMax)*4,
     *                  lstcod(ncodMax)*20,
     *                  lstwav(nwavMax)*6
        Common  /atmx0/ lstatm,natm,katm,ismatm,idatm
        Common  /codx0/ lstcod,ncod,kcod,kelem,ismcod,idcod
        Common  /wavx0/ lstwav,nwav,kwav,ismwav,idwav
        Real            a(6),vol
        Common  /back1/ a,vol,lbev
        Character       txt(20)*80
        Common  /msg/   txt
c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer                 iHenkeCowan
        Common  /HenkeCowan/    iHenkeCowan
c Acceptable values of iHenkeCowan:
c -1 - automatic choice of internal or external dispersion correction DB
c  0 - do not use external dispersion correction databases
c  1 - use henke.dat    for df1
c  2 - use henke.dat    for df1, df2 (Crosec)
c  3 - use cowan.dat    for df1
c  4 - use cowan.dat    for df1, df2 (Crosec)
c  5 - use windt.dat    for df1
c  6 - use windt.dat    for df1, df2 (Crosec)
c  7 - use chantler.dat for df1
c  8 - use chantler.dat for df1, df2 (Crosec)
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------
c ATTENTION - 2:
c X0H calculates crystal susceptibilities for expanding the wavefields
c over Bloch waves like exp(ikr). Then crystal susceptibilities are
c calculated by the summation over reciprocal unit cell with exp(-ikr).
c If you use expansion over exp(ikr), then full susceptibilities are
c presented as:
c
c xh =   xrh*exp(i*pi*xrf) +
c    + i*xih*exp(i*pi*xif) +
c    + i*xqh*exp(i*pi*xqf)
c
c If you use expansion over exp(-ikr), (for example, like in the Pinsker
c book on dynamical X-ray diffraction), then full susceptibilities are
c presented as:
c
c xh =   xrh*exp(-i*pi*xrf) -
c    - i*xih*exp(-i*pi*xif) -
c    - i*xqh*exp(-i*pi*xqf)
c
c Another solution for switching to the Pinsker notation is to change
c the factor ipm from ipm=-1 to ipm=1.
c
c                             December 10, 1992 -  Stepanov
c -------------------------------------------------------
        data    ipm     /-1/
        data    im      /(0.,1.)/
        data    smtry   /'Cubic',
     *                   'Hexagonal/Trigonal',
     *                   'Tetragonal',
     *                   'Trigonal(Rhombohedral)',
     *                   'Orthorhombic',
     *                   'Monoclinic',
     *                   'Triclinic',
     *                   '** Amorphous ** '/

c These are constants defining dimension of arrays:
c +=======================+
        natm = natmMax   !|
        ncod = ncodMax   !|
        nwav = nwavMax   !|
c +=======================+
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

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'X0h2'               !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+

        pi    = 4.*atan(1.)             !set the "pi" constant
        gra   = 2.*pi/360.              !number of radians in 1 degree
        sec   = gra/3600.               !number of radians in 1 arc sec
        e_radius  = 2.81794e-05         !classical electron radius (A)
        Compton   = 2.42631e-02         !Compton length (A)
        Planck_   = 1.054573e-27        !Planck constant (erg*c)/2*pi
        Rydberg   = 1.097373e-03        !Ridberg constant (1/A)
        Boltzman  = 1.380658e-16        !Boltzmann constant (erg/k)
        Avogadro  = 6.0221367e+23       !Avogadro constant
c       aem       = 1.6605402e-24       !atomic mass constant (grams)
c Avogadro = the number of atoms or molecules contained in a mole.
c Mole is defined as a mass in grams equal to the atomic or molecular
c weight of a substance.
c The atomic mass constant  (atomnaja edinitca massi) is the 1/12
c of the mass of Si atom expressed in grams [==1/Avogadro]:

        linesX0hWarn   = 0              !no warnings yet
        md_NewBase     =-1
        md_NewAll      = 0
        md_NewCrystal  = 1              !use new name,prcn,a(1-6) from X0h DB
        md_NewCrystal_a= 2              !use new name,prcn from X0h DB, use ad(1-6) if non-zero, otherwise from X0h DB
        md_NewWave     = 3
        md_NewLine     = 4
        md_NewHKL      = 5
        mode           = mode_x0h
        if (mode.eq.md_NewBase) Then  !----------------+
c Re-read COORD.X0h index (special mode for X0hWin):   |
          idcod = 0                                   !|
          mode = md_NewAll                            !|
        endif  !---------------------------------------+
c -------------------------------------------------------
c Set X0H database files LUNs if not passed:
        if (luwa.eq.0)  luwa=11         !"wave.x0h"
        if (luco.eq.0)  luco=12         !"coord.x0h"
        if (luat.eq.0)  luat=13         !"atom.x0h"
        if (lusp.eq.0)  lusp=14         !"space.x0h"
c .......................................................
c Save a copy of X-ray wavelength:
        wave = wv
        if (radiat(1:1).eq.char(0))                     radiat=' '
        if (mode.eq.md_NewLine .OR.
     *      mode.eq.md_NewAll) Then !----------+
           if (wv.gt.0.) radiat=' '           !|
        endif !--------------------------------+

c Is the structure code known?
        if (code(1:1).eq.char(0))       code=' '
        lcod = Max(Len_Trim(code),1)

c Crystal syngony:
        if (ising.lt.0 .or. ising.gt.8) ising=0
        do      indx=1,kcompMax  !========================+
          if (name(indx)(1:1).eq.char(0)) name(indx)=' ' !|
        enddo  !==========================================+
c Place of the incident beam in the indices array:
        incid = inciden
        if (incid.lt.1 .or. incid.gt.nref) incid=1

c Save parameters of the previous X0H runs?
          if (mode.eq.md_NewAll     .OR.
     *        mode.eq.md_NewCrystal .OR.
     *        mode.eq.md_NewCrystal_a) Then !+
          do  indx=1,kcompMax !+             |
            name(indx) = ' '  !|             |
            prcn(indx) = 0.   !|             |
          enddo  !=============+             |
        endif  !-----------------------------+
        if (mode.eq.md_NewCrystal .OR.
     *      mode.eq.md_NewAll) Then !--------------------------------+
          do  i=1,6  !===+                                           |
            ad(i) = 0.  !|                                           |
          enddo  !=======+                                           |
          lbev = 0                                                  !|
          Call  DlatReset ()                                        !|
          Poisson = 0.                                              !|
        else   !-----------------------------------------------------+
c Filter out unreasonable predefines:                                |
          do  i=1,3  !============================================+  |
            if     (ad(i).lt.0.) Then  !-----------------------+  |  |
               write (0,100) 'X0h1',mode,code(1:lcod),i, a(i) !|  |  |
               a(i) = 0.                                      !|  |  |
            elseif (ad(i).gt.1.E+3) Then !---------------------+  |  |
               write (0,100) 'X0h1',mode,code(1:lcod),i, a(i) !|  |  |
               a(i) = 0.                                      !|  |  |
            elseif (ad(i).lt.1.E-37) Then !--------------------+  |  |
               a(i) = 0.                            !no warning|  |  |
            elseif (ad(i).lt.1.E-3) Then !---------------------+  |  |
               write (0,100) 'X0h1',mode,code(1:lcod),i, a(i) !|  |  |
               a(i) = 0.                                      !|  |  |
            endif !--------------------------------------------+  |  |
          enddo  !================================================+  |
          do  i=4,6  !============================================+  |
            if     (ad(i).lt.0.) Then  !-----------------------+  |  |
               write (0,100) 'X0h1',mode,code(1:lcod),i, a(i) !|  |  |
               a(i) = 0.                                      !|  |  |
            elseif (ad(i).gt.179.) Then !----------------------+  |  |
               write (0,100) 'X0h1',mode,code(1:lcod),i, a(i) !|  |  |
               a(i) = 0.                            !no warning|  |  |
            elseif (ad(i).lt.1.E-37) Then !--------------------+  |  |
               a(i) = 0.                                      !|  |  |
            elseif (ad(i).lt.1.) Then !------------------------+  |  |
               write (0,100) 'X0h1',mode,code(1:lcod),i, a(i) !|  |  |
               a(i) = 0.                                      !|  |  |
            endif  !-------------------------------------------+  |  |
          enddo  !================================================+  |
        endif  !-----------------------------------------------------+
  100   format(1x,a,' incorrect predefinition for mode=',i2,
     *         ', code=',a,': a=(',i1,') = ',g12.5)
        if (mode.ne.md_NewHKL) Then !-+
          edge_nearest = -100.       !|
          name_nearest = ' '         !|
        endif  !----------------------+

c Save a copy of passed print option:
        ipr   = jpr
        ifail = 0
c .......................................................
c Make catalogs of the database files:
        if (idcod.eq.0) Then  !------------------------+
          Call Basini ('coord.x0h',lstcod,ncod,kcod,  !|
     *                  ismcod,idcod,1,kelem,ifail)   !|
          if (ifail.ne.0) goto 998                    !|
          kcod = kelem                                !|
        endif  !---------------------------------------+
        if (idatm.eq.0) Then  !------------------------+
          Call Basini ('atom.x0h',lstatm,natm,katm,   !|
     *                  ismatm,idatm,0,i,ifail)       !|
          if (ifail.ne.0) goto 998                    !|
        endif  !---------------------------------------+
        if (idwav.eq.0 .AND. abs(wave).lt.1.E-37) Then !---+
          Call Basini ('wave.x0h',lstwav,nwav,kwav,       !|
     *                  ismwav,idwav,0,i,ifail)           !|
          if (ifail.ne.0) goto 998                        !|
        endif  !-------------------------------------------+
c .......................................................
c Determine X-ray wavelength:

        itrace = 0
        if (mode.ne.md_NewHKL)  Then  !------+
          Call  Redwav (wave,itrace,ifail)  !|
          if (ifail.ne.0)  goto 998         !|
        endif  !-----------------------------+

        energy = wave2energy(wave)
c .......................................................
c Make a copy of crystal lattice parameters:
        do      i=1,6  !=====+
          a(i) = ad(i)      !|
        enddo  !=============+
c ....................................................
c Read the structure description (atomic coordinates):
c The parameters are returned in common /x0pa?/

        if (mode.eq.md_NewAll     .OR.
     *      mode.eq.md_NewCrystal .OR.
     *      mode.eq.md_NewCrystal_a) Then !--------+
          iimode = 0                              !|
          Call  RedCor  (code,iimode,rho,ifail)   !|
          if (ifail.gt.0)  goto 998               !|
        endif  !-----------------------------------+
c (with negative ifail we still may continue...)
        if (isym.ne.2)  Then  !---+
          inmax = 3              !|
        else  !-------------------+
          inmax = 4              !|
        endif  !------------------+

        lcod = Max(Len_Trim(Code),1)
        if (ising.ne.0 .AND. ising.ne.isym)  Then !------+
          write (txt,813) code(1:lcod),smtry(ising)     !|
  813     format('X0H -- E R R O R:'//
     *    'Symmetry of structure [',a,'] is not ',a)    !|
          ifail = 3                                     !|
          goto 998                                      !|
        endif  !-----------------------------------------+

c Poisson is used with cubic structures only:
        if (isym.gt.1)  Poisson = 0.
c ......................................................
c Request crystal lattice parameters is some of them are
c not present in the "coord.x0h" file:

        if (isym.ne.8)  Then  !--------------------------+for non-amorphous
                                                        !|
          if (mode.eq.md_NewAll     .OR.                !|
     *        mode.eq.md_NewCrystal .OR.                !|
     *        mode.eq.md_NewCrystal_a) Then !-+          |
c This may return negative ifail              |          |
c (incomplete data):                          |          |
            Call  RedLat (code,need_a,ifail) !|          |
            if (ifail.gt.0) goto 998         !|          |
c (For negative ifail we continue...)         |          |
          endif  !----------------------------+          |
                                                        !|
c .................................................      |
c Determine the crystal reciprocal lattice parameters:   |
c (returned via common/back1/ and common/back2/)         |
                                                        !|
          if (mode.eq.md_NewAll     .OR.                !|
     *        mode.eq.md_NewCrystal .OR.                !|
     *        mode.eq.md_NewCrystal_a) Then !--+         |
            lbev = 0                          !|         |
            if (ifail.eq.0)  Call DBacklat()  !|         |this calculates vol
          endif  !-----------------------------+         |
                                                        !|
          if (lbev .ne.1) stop 'X0h2: no back lattice'  !|
                                                        !|
        endif  !-----------------------------------------+

c Exit if no structure components are present:
        if (kcomp.eq.0) goto 613

c -------------------------------------------------------
c             Read data for calculating X0:
c ......................................................
        if (mode.ne.md_NewHKL) Then !-----------------------------------------------+
          do      indx=1,kcomp  !================================================+  |
c .......................................................                        |  |
c Read from "atom.x0h" under atom name the atom scattering parameters.           |  |
c The parameters are returned via common/param/. Also calculate f0 using         |  |
c subroutine forfac and check the condition z=f0. Calculate the absorption       |  |
c cross sections td and tq with the help of subroutine crosec and the Debye      |  |
c coefficient (read it from "atom.x0h" or calculate from the Debye temperature   |  |
c in subroutine fifunc):                                                         |  |
                                                                                !|  |
            rho_i = 0.                       !this is used for x0                |  |
c This may return negative ifail                                                 |  |
c (incomplete data):                                                             |  |
            Call  RedPar (code,rho_i,atom_mass_aem(indx),ifail)                 !|  |
            if (ifail.gt.0) goto 998                                            !|  |
c (For negative ifail we continue...)                                            |  |
                                                                                !|  |
            if (abs(prcn(indx)).lt.1.E-37) prcn(indx)=1.                        !|  |
c Commented in August-99. Now we only use range [0.-1.]                         !|  |
c ccc       if (prcn(indx).gt.1.) prcn(indx)=prcn(indx)/100.                    !|  |
                                                                                !|  |
            if (prcn(indx).lt.0. .OR.                                           !|  |
     *          (isym.ne.8 .AND. prcn(indx).gt.1.)) Then !--------------------+  |  |
                                                                             !|  |  |
              if ((isym.ne.8 .AND. prcn(indx).gt.1.)                         !|  |  |
     *                       .AND. modebat.eq.0) Then !--+                    |  |  |
                write (txt,44)  code(1:lcod),           !|                    |  |  |
     *                          name(indx),             !|                    |  |  |
     *                          prcn(indx)              !v                    v  v  v
  44            format('Structure: ',a,'  Element: ',a/
     *          '  Occupation=',f8.3/
     *          ' -- Occupation not in range [0.-1.]')  !^                    ^  ^  ^
                Call  Message (txt,3,2)                 !|                    |  |  |
              endif  !-----------------------------------+                    |  |  |
                                                                             !|  |  |
              txt(1) = 'Structure: '//code                                   !|  |  |
              txt(2) = 'Enter occupation for ['//name(indx)//']:'            !|  |  |
              txt(3)=' '                                                     !|  |  |
              txt(4)='This request is not supported'//                       !|  |  |
     *               'in the current version of x0h!'                        !|  |  |
              ifail = -4                                                     !|  |  |
              goto  995                                                      !|  |  |
                                                                             !|  |  |
            endif  !----------------------------------------------------------+  |  |
                                                                                !|  |
            if (ipr.eq.2)   Then   !-------------------------+                   |  |
              write (txt,115)  name(indx),                  !|                   |  |
     *                         kol(indx),                   !|                   |  |
     *                         prcn(indx)                   !|                   |  |
                Call  Priblo  (txt,1,1,0,0,ifa)             !|                   |  |
              if (isym.ne.8)  Then  !---------------------+  |                   |  |for crystals only
                write (txt,116) (i,                      !|  |                   |  |
     *                         wx(indx,i),               !|  |                   |  |
     *                         wy(indx,i),               !|  |                   |  |
     *                         wz(indx,i),               !|  |                   |  |
     *                         i=1,kol(indx))            !|  |                   |  |
                Call  Priblo  (txt,1+kol(indx),1,0,0,ifa)!|  |                   |  |
              endif  !------------------------------------+  |                   |  |
            endif  !-----------------------------------------+                   |  |
                                                                                !V  V
  115       format('Element: ',a,'  Count of atoms=',i2,
     *      '  Occupation=',f6.3)
  116       format('Coordinates:',100(/3x,i2,
     *      '.  x=',g12.5,'  y=',g12.5,'  z=',g12.5,:))                         !^  ^
                                                                                !|  |
  995       continue                                                            !|  |
          enddo  !===============================================================+  |
                                                                                   !|
c For mono-atomic amorphous materials take density from ATOM.X0H:                   |
          if (abs(rho).lt.1.E-10 .AND.                                             !|
     *      kcomp.eq.1           .AND.                                             !|
     *     kol(1).eq.1           .AND.                                             !|
     *       isym.eq.8) Then  !---------+                                           |
            rho = rho_i                !|                                           |
          endif  !----------------------+                                           |
        endif  !--------------------------------------------------------------------+

        if     (ifail.gt.0) Then  !---------------------------------+
          goto 998   !----------------report-error------------------+--+
        elseif (ifail.lt.0) Then  !---------------------------------+  v
          txt(1)='The following structure parameter(s) '//         !|
     *           'required for ['//Code(1:lcod)//']:'              !|
          do i=2,20  !=====+                                        |
            txt(i) = ' '  !|                                        |
          enddo  !=========+                                        |
          j = 1         !starting position in txt line             !|
          ll= 2         !starting txt line number                  !|
          ir= 0         !there is no error message yet             !|
          if (isym.ne.8) Then  !--------------------------------+   |
            do i=1,6  !=====================================+   |   |
              if (need_a(i).ne.0 .and. a(i).le.0.) Then !-+ |   |   |
                ir= 1            !got error message      !| |   |   |
c cccccccc      m = 2                                    !| |   |   |
c cccccccc      txt(ll)(j:j+m-1) = 'a?'                  !| |   |   |
                m = 17                                   !| |   |   |
                txt(ll)(j:j+m-1) = 'Lattice_constant?'   !| |   |   |
                write (txt(ll)(j+m-1:j+m-1),'(i1)') i    !| |   |   |
                j = j+m-1+2                              !| |   |   |
                if (j.gt.64)  Then  !-+                   | |   |   |
                  j  = 1             !|                   | |   |   |
                  ll = ll+1          !|                   | |   |   |
                endif  !--------------+                   | |   |   |
              endif  !------------------------------------+ |   |   |
            enddo  !========================================+   |   |
          else  !-----------------------------------------------+   |
            if (rho.le.0.) Then  !----------------+             |   |
              ir= 1         !got error message    |             |   |
c cccccccc    m = 3                              !|             |   |
c cccccccc    txt(ll)(j:j+m-1) = 'rho'           !|             |   |
              m = 13                             !|             |   |
              txt(ll)(j:j+m-1) = 'Density (rho)' !|             |   |
              j = j+m-1+2                        !|             |   |
            endif  !------------------------------+             |   |
          endif  !----------------------------------------------+   |
          if (Len_Trim(txt(ll)).gt.0) Then !---+start from          |
            j  = 1                            !|a new line!         |
            ll = ll+1                         !|                    |
          endif  !-----------------------------+                    |
          do i=1,kcomp  !====================================+      |
            if (Len_Trim(name(i)).eq.0)  Then  !-----------+ |      |
              ir= 1         !got error message             | |      |
              m = 6                                       !| |      |
              txt(ll)(j:j+m-1) = 'Atom??'                 !| |      |
              write (txt(ll)(j+m-2:j+m-1),'(i2)') i       !| |      |
              if (txt(ll)(j+m-2:j+m-2).eq.' ')            !| |      |
     *            txt(ll)(j+m-2:j+m-2)='0'                !| |      |
              j = j+m-1+2                                 !| |      |
              if (j.gt.64)  Then  !-+                      | |      |
                j  = 1             !|                      | |      |
                ll = ll+1          !|                      | |      |
              endif  !--------------+                      | |      |
            endif  !---------------------------------------+ |      |
          enddo  !===========================================+      |
          if (Len_Trim(txt(ll)).gt.0) Then !+start from             |
            j  = 1                         !|a new line!            |
            ll = ll+1                      !|                       |
          endif  !--------------------------+                       |
          do i=1,kcomp  !====================================+      |
            if (prcn(i).lt.0. .OR. prcn(i).gt.1.) Then  !--+ |      |
              ir= 1         !got error message             | |      |
              m = 12                                      !| |      |
              txt(ll)(j:j+m-1) = 'Occupation??'           !| |      |
              write (txt(ll)(j+m-2:j+m-1),'(i2)') i       !| |      |
              if (txt(ll)(j+m-2:j+m-2).eq.' ')            !| |      |
     *            txt(ll)(j+m-2:j+m-2)='0'                !| |      |
              j = j+m-1+2                                 !| |      |
              if (j.gt.64)  Then  !-+                      | |      |
                j  = 1             !|                      | |      |
                ll = ll+1          !|                      | |      |
              endif  !--------------+                      | |      |
            endif  !---------------------------------------+ |      |
          enddo  !===========================================+      |
          if (ir.ne.0) Then  !----------------------+we have err msg|
            if (Len_Trim(txt(ll)).eq.0) ll=ll-1    !|               |
            ifail = -ll                            !|               |
            goto 613  !-----to-exit-----------------+---------------+--+
          endif  !----------------------------------+               |  v
          ifail = 0                  !the problem has been resloved |(e.g. rho was found)
        endif  !----------------------------------------------------+

        lcod = Max(Len_Trim(Code),1)
c Acceptable values of iHenkeCowan:
c -1 - automatic choice of internal or external dispersion correction DB
c  0 - do not use external dispersion correction databases
c  1 - use henke.dat    for df1
c  2 - use henke.dat    for df1, df2 (Crosec)
c  3 - use cowan.dat    for df1
c  4 - use cowan.dat    for df1, df2 (Crosec)
c  5 - use windt.dat    for df1
c  6 - use windt.dat    for df1, df2 (Crosec)
c  7 - use chantler.dat for df1
c  8 - use chantler.dat for df1, df2 (Crosec)
        if (iHenkeCowan.eq.0 .AND.
     *      (wave.lt.0.1 .OR. wave.gt.10.)) Then  !--------------------+
          write (txt,111) Code(1:lcod), wave, energy                  !|
          ifail = 6                                                   !|
          if (energy.ge.0.01 .AND. energy.le.700.) Then !------+       |
            ifail = ifail+2                                   !|       |
            txt(ifail) = 'Please choose Henke or '//          !|       |
     *                   'Cowan-Brennan tables.'              !|       |
          endif !----------------------------------------------+       |
          return                                                      !|
        elseif ((iHenkeCowan.eq.1 .OR. iHenkeCowan.eq.2) .AND.        !|
     *      (energy.lt.0.01 .OR. energy.gt.30.)) Then  !---------------+henke.dat
          write (txt,211) Code(1:lcod), wave, energy                  !|
          ifail = 6                                                   !|
          if (energy.ge.0.03 .AND. energy.le.700.) Then !------+       |
            ifail = ifail+2                                   !|       |
            txt(ifail) = 'Please use Cowan-Brennan data '//   !|       |
     *                   'that are valid for [0.003-700] KeV' !|       |
          endif !----------------------------------------------+       |
          return                                                      !|
        elseif ((iHenkeCowan.eq.3 .OR. iHenkeCowan.eq.4) .AND.        !|
     *      (energy.lt.0.03 .OR. energy.gt.700.)) Then !---------------+cowan.dat
          write (txt,311) Code(1:lcod), wave, energy                  !|
          ifail = 6                                                   !|
          if (energy.ge.0.01 .AND. energy.le.30.) Then  !------+       |
            ifail = ifail+2                                   !|       |
            txt(ifail) = 'Please use Henke data that '//      !|       |
     *                   'are valid for [0.001-30] KeV'       !|       |
          endif !----------------------------------------------+       |
          return                                                      !|
        elseif ((iHenkeCowan.eq.5 .OR. iHenkeCowan.eq.6) .AND.        !|
     *      (energy.lt.0.01 .OR. energy.gt.100.)) Then !---------------+windt.dat
          write (txt,411) Code(1:lcod), wave, energy                  !|
          ifail = 6                                                   !|
          if (energy.gt.100. .AND. energy.le.700.) Then !------+       |
            ifail = ifail+2                                   !|       |
            txt(ifail) = 'Please use Cowan-Brennan data '//   !|       |
     *                   'that are valid up to 700 KeV'       !|       |
          endif !----------------------------------------------+       |
          return                                                      !|
        elseif ((iHenkeCowan.eq.7 .OR. iHenkeCowan.eq.8) .AND.        !|
     *      (energy.lt.0.01 .OR. energy.gt.450.)) Then !---------------+chantler.dat
          write (txt,511) Code(1:lcod), wave, energy                  !|
          ifail = 6                                                   !|
          if (energy.gt.450. .AND. energy.le.700.) Then !------+       |
            ifail = ifail+2                                   !|       |
            txt(ifail) = 'Please use Cowan-Brennan data '//   !|       |
     *                   'that are valid up to 700 KeV'       !|       |
          endif !----------------------------------------------+       |
          return                                                      !|
        elseif (iHenkeCowan.eq.-1 .AND.                               !|
     *      (energy.lt.0.01 .OR. energy.gt.700.)) Then !---------------+automatic
          write (txt,911) Code(1:lcod), wave, energy                  !|
          ifail = 6                                                   !|
          return                                                      !|
        endif  !-------------------------------------------------------+
  111   format(
     *  '*** Wavelength not in the [0.1-10.0] Angstroms range ***'/
     *  'X0h cannot calculate scattering factors for material'/
     *  a/
     *  'at wavelength=',g10.4,' A  (energy=',g10.4,' keV)'/
     *  'because there is no data in the database on the'/
     *  'dispersion corrections df1, df2 for this wavelength.'//)
  211   format(
     *  '*** X-ray energy not in the [0.01-30.0] KeV range ***'/
     *  'X0h cannot calculate scattering factors for material'/
     *  a/
     *  'at wavelength=',g10.4,' A  (energy=',g10.4,' keV)'/
     *  'because there is no data in the Henke database on the'/
     *  'dispersion corrections df1, df2 for this wavelength.'//)
  311   format(
     *  '*** X-ray energy not in the [0.03-700.0] KeV range ***'/
     *  'X0h cannot calculate scattering factors for material'/
     *  a/
     *  'at wavelength=',g10.4,' A  (energy=',g10.4,' keV)'/
     *  'because there is no data in the Brennan-Cowan DB on the'/
     *  'dispersion corrections df1, df2 for this wavelength.'//)
  411   format(
     *  '*** X-ray energy not in the [0.01-100.0] KeV range ***'/
     *  'X0h cannot calculate scattering factors for material'/
     *  a/
     *  'at wavelength=',g10.4,' A  (energy=',g10.4,' keV)'/
     *  'because there is no data in the Windt DB on the'/
     *  'dispersion corrections df1, df2 for this wavelength.'//)
  511   format(
     *  '*** X-ray energy not in the [0.01-450.0] KeV range ***'/
     *  'X0h cannot calculate scattering factors for material'/
     *  a/
     *  'at wavelength=',g10.4,' A  (energy=',g10.4,' keV)'/
     *  'because there is no data in the Chantler/NIST DB on the'/
     *  'dispersion corrections df1, df2 for this wavelength.'//)
  911   format(
     *  '*** X-ray energy not in the [0.01-700.0] KeV range ***'/
     *  'X0h cannot calculate scattering factors for material'/
     *  a/
     *  'at wavelength=',g10.4,' A  (energy=',g10.4,' keV)'/
     *  'because there is no data in any available databases on the'/
     *  'dispersion corrections df1, df2 for this wavelength.'//)

c ------------------------------------------------------
c                    Calculation of x0
c ......................................................
        ucell_mass_gram = 0.
        r_atoms_ucell   = 0.0
        do      i=1,kcomp  !=====================================+
          ucell_mass_gram = ucell_mass_gram                     !|
     +                    + kol(i) * prcn(i) * atom_mass_aem(i) !|
          r_atoms_ucell   = r_atoms_ucell                       !|
     +                    + kol(i) * prcn(i)                    !|
        enddo  !=================================================+
        n_atoms_ucell = int(r_atoms_ucell + 0.1)        !changed 2013.08
        ucell_mass_gram = ucell_mass_gram / Avogadro

c (1.e+24) is the recalculation factor from cm^3 to A^3:
        if (isym.ne.8)  Then  !------------------------------------------+
          if (abs(vol).lt.1.E-37) stop 'X0h2: crystal cell vol=0'       !|
          rho = ((1.e+24)*ucell_mass_gram)/vol                          !|
        else  !----------------------------------------------------------+
          if (rho.le.0.) Then  !----------------------------------+      |
            lcod = Max(Len_Trim(Code),1)                         !|      |
            if (CodeType(code).lt.2) Then !----------+            |      |
              write (txt,112)  code(1:lcod), rho    !|            |      |
              ifail = 3                             !|            |      |
            else !-----------------------------------+            |      |
              write (txt,113)  rho, code(1:lcod)    !|            |      |
              ifail = 8                             !|            |      |
            endif !----------------------------------+            |      |
  112       format(' X0h1: Error in db_control for [',a,']:'/
     *      ' Calculated rho=',g10.3,' is not positive number'/
     *      ' Please, report to the author.')                    !|      |
  113       format(' X0h1: Incorrect or unspecified rho=',g10.3/
     *      ' for chemical formula [',a,'].'/
     *      ' This material is not present in the X0H database.'/
     *      ' It is treated as an amorphous compound and requires'/
     *      ' specifying material density rho (g/cm^3).'/
     *      ' If you need treating this material as a crystal,'/
     *      ' please submit the structure via the X0h web form.'/
     *      ' Otherwise, please input correct rho parameter.')   !|      |
            return                                               !|      |
          endif  !------------------------------------------------+      |
          if (abs(rho).lt.1.E-37) stop 'X0h2: amorphous material rho=0' !|
          vol = ((1.e+24)*ucell_mass_gram)/rho                          !|
        endif  !---------------------------------------------------------+
c This is the density per (cm**3):
        atom_density_cm3 = (1.E+24)*n_atoms_ucell / vol

        cr = -wave*wave*e_radius/(pi*vol)
        ci = -ipm*wave/(2.*pi*vol)
c .............|..........................................
c  +-----------+
c ipm=1 corresponds to convention by Pinsker (p.18,71). At this
c convention the Bloch waves are expanded as exp(-i*k*r).
c ipm=-1 corresponds to convention by Afanasiev and Baryshevsky.
c At this convention the Bloch waves are expanded as exp(i*k*r).
c ......................................................
        xr0 = 0.
        xi0d= 0.
        xi0q= 0.
        do      i=1,kcomp  !=========================+
          xr0  = xr0  + cr*kol(i)*prcn(i)*f0(i)     !|
          xi0d = xi0d + ci*kol(i)*prcn(i)*td(i)     !|
          xi0q = xi0q + ci*kol(i)*prcn(i)*tq(i)     !|
        enddo  !=====================================+

        xi0 = xi0d + xi0q

c Absorption coefficient (1/cm):
        if (abs(wave).gt.1.E-37) Then  !--------+
          mu0 = (2.E+08)*pi*abs(xi0)/wave      !|
        else  !---------------------------------+
          mu0 = 0.                             !|
        endif  !--------------------------------+
c Absorption length (micron):
        if (abs(xi0).gt.1.E-37) Then  !---------+
          a0 = (1.E-04)*wave/(2.*pi*abs(xi0))  !|
        else  !---------------------------------+
          a0 = 0.                              !|
        endif  !--------------------------------+
c Critical angle of total reflection:
        TER = Sqrt(Abs(xr0))
c Penetration depth in TER region:
        if (abs(TER).gt.1.E-37) Then  !---------+
          dl_TER = wave / (2.*pi*TER)          !|
        else  !---------------------------------+
          TER = 0.                             !|
        endif  !--------------------------------+
c delta and eta from n=1-delta-ieta:
        delta = -0.5*xr0
        eta   = -0.5*xi0
c ....................................................
        if (ipr.ne.0 .and. mode.ne.md_NewHKL) Then !------+
          write   (txt,101)     code,                    !|
     *                          smtry(isym),             !|
     *                          rho                      !|
          Call  Insert_Poisson (Poisson,isym,txt(3))     !|
          n = 4                                          !|
          if (isym.ne.8)  Then  !---------------+         |
            n = n+1                            !|         |
            write (txt(n),102)  (a(i),i=1,3)   !|         |
            n = n+1                            !|         |
            write (txt(n),103)  (a(i),i=4,6)   !|         |
          endif  !------------------------------+         |
          n = n+1                                        !|
          write   (txt(n),104)  name(1),kol(1),          !|
     *                          prcn(1)                  !|
          do i=2,kcomp  !=======================+         |
            n = n+1                            !|         |
            write (txt(n),105)  name(i),kol(i),!|         |
     *                          prcn(i)        !|         |
          enddo  !==============================+         |
          Call  Priblo  (txt,n,0,0,0,ifa)                !|
          write (txt,110)       wave, radiat, energy,    !|
     *                          xr0, xi0,                !|
     *                          delta, eta,              !|
     *                          TER/gra, TER*1000.,      !|
     *                          dl_TER, mu0, a0          !|
          Call  Insert_Edge    (edge_nearest,            !|
     *                          name_nearest,txt(3))     !|
          Call  Priblo  (txt,10,0,0,0,ifa)               !|
        endif  !------------------------------------------+
  101   format  (
     *  '+',      46('-'),                             '-',
     *  24('-'),         '+'/
     *  '|Material code',33x,                          '|',
     *  a,4x,            '|'/
     *  '|Material symmetry',29x,                      '|',
     *  a,               '|'/
     *  '|Material density (gm/cm^3)',20x,             '|',
     *  g12.5,12x,       '|')
  102   format  (
     *  '|Lattice spacings (Angstrom)',19x,            '|',
     *  3f8.2,           '|')
  103   format  (
     *  '|Lattice angles (Degrees of arc)',15x,        '|',
     *  3f8.2,           '|')
  104   format  (
     *  '|Composition: atom -- N_sites (site occupation)|',
     *  a,' --',i3,4x,'(',f7.3,') |')
  105   format  (
     *  '|',46x,                                       '|',
     *  a,' --',i3,4x,'(',f7.3,') |')
  110   format  (
     *  '|',      46('-'),                             '+',
     *  24('-'),         '|'/
     *  '|X-ray wavelength (Angstrom),    [x-ray line]  |',
     *  g11.5,',  [',a,']  |'/
     *  '|X-ray energy (kev)',     28x,                '|',
     *  g11.5,   13x,    '|'/
     *  '|',      46('-'),                             '+',
     *  24('-'),         '|'/
     *  '|xr0, xi0: (xr0=Re(x0), xi0=Im(x0), n = 1+x0/2)|',
     *  e11.4,', ',e11.4,'|'/
     *  '|delta, eta:  (n = 1-delta-i*eta)              |',
     *  e11.4,', ',e11.4,'|'/
     *  '|Total reflection: critical angle (degr.,mrad) |',
     *  g11.4,', ',g11.4,'|'/
     *  '|Total reflection: penetration depth (Angstrom)|',
     *  g11.4,13x,       '|'/
     *  '|Linear absorption factor (1/cm) & length (um) |',
     *  g11.4,', ',g11.4,'|'/
     *  '+',      46('-'),                             '-',
     *  24('-'),         '+')

c If there are no reflections, then go to exit:
        if (nref.le.1)  goto 613  !------>---------------+
c                                                        v
c ------------------------------------------------------
c             +=========================+
c             |   Calculation of xh:    |
c             +=========================+
c ......................................................
c 0. Verify the the 0th reflection indices are all zero:
        if (index(1,incid).ne.0   .OR.
     *      index(2,incid).ne.0   .OR.
     *      index(3,incid).ne.0)        Then  !-------------------+
          write (txt,913) incid,code(1:lcod)                     !v
  913     format('X0H -- E R R O R:'//
     *    'Zero reflection ',i2,' from [',a,'] is not (0,0,0)')  !^
          ifail = 3                                              !|
          goto 998                                               !|
        endif  !--------------------------------------------------+

c Amorphous material:
        if (isym.eq.8)  Then  !---------------------+
          do    i=1,nref  !======================+  |
            QB(i) = 0.                          !|  |
            do  j=1,nref  !====================+ |  |
              ind(1) = index(1,i)-index(1,j)  !| |  |
              ind(2) = index(2,i)-index(2,j)  !| |  |
              ind(3) = index(3,i)-index(3,j)  !| |  |
              if (ind(1).eq.0 .AND.           !| |  |
     *            ind(2).eq.0 .AND.           !| |  |
     *            ind(3).eq.0) Then !---+      | |  |
                xrh(i,j) = xr0         !|      | |  |
                xih(i,j) = xi0d        !|      | |  |
                xqh(i,j) = xi0q        !|      | |  |
              else  !-------------------+      | |  |
                xrh(i,j) = 0.          !|     !| |  |
                xih(i,j) = 0.          !|     !| |  |
                xqh(i,j) = 0.          !|     !| |  |
                xrf(i,j) = 0.          !|     !| |  |
                xif(i,j) = 0.          !|     !| |  |
                xqf(i,j) = 0.          !|     !| |  |
              endif  !------------------+      | |  |
            enddo  !===========================+ |  |
          enddo  !===============================+  |
          goto 613  !-------------------------------+----+
        endif  !------------------------------------+    v

c 1. Loop over all Bragg reflections and verify that all indices are
c passed. Request missing indices from terminal (in orthogonal form):

        do      iref=1,nref    !=============================+
          ind(1)     = index(1,iref)                        !|
          ind(2)     = index(2,iref)                        !|
          ind(inmax) = index(3,iref)                        !|
          if (isym.eq.2)        Then  !---+                  |
            ind(3)=-(ind(1)+ind(2))      !|                  |
          endif  !------------------------+                  |
c If indices are not passed, then request:                   |
          if (ind(1)     .eq. 0     .AND.                   !|
     *        ind(2)     .eq. 0     .AND.                   !|
     *        ind(inmax) .eq. 0     .AND.                   !|
     *        iref       .ne. incid)    Then   !----------+  |
            txt(1) = 'Reflection    indices are '//      !|  |
     *                               'not specified!'    !|  |
            if (nref.gt.1)                               !|  |
     *      write (txt(1)(12:13),'(i2)') iref            !|  |
            ifail = 1                                    !|  |
            goto 998                                     !|  |
          endif  !----------------------------------------+  |
        enddo  !=============================================+

c=====================================================
c 2. Dual loop to fill matrices xrh(nref,nref)
c                               xih(nref,nref)
c                               xqh(nref,nref)
c - in such a way that xh(i,i) contain x0 and the other elements
c   contain xh(i,j)=xh(i-j):

        do      i=1,nref  !===================================+
        do      j=1,nref  !================================+  |
          ind(1)     = index(1,i)-index(1,j)              !|  |
          ind(2)     = index(2,i)-index(2,j)              !|  |
          ind(inmax) = index(3,i)-index(3,j)              !|  |
          if (isym.eq.2)      Then  !--+                   |  |
            ind(3)=-(ind(1)+ind(2))   !|                   |  |
          endif  !---------------------+                   |  |
                                                          !|  |
          if (ind(1)     .eq. 0     .AND.                 !|  |
     *        ind(2)     .eq. 0     .AND.                 !|  |
     *        ind(inmax) .eq. 0)        Then  !---------+  |  |
c Diagonal matrix elements:                             |  |  |
            xrh(i,j) = xr0                             !|  |  |
            xih(i,j) = xi0d                            !|  |  |
            xqh(i,j) = xi0q                            !|  |  |
            xrf(i,j) = 0.                              !|  |  |added 2014/02/15
            xif(i,j) = 0.                              !|  |  |added 2014/02/15
            xqf(i,j) = 0.                              !|  |  |added 2014/02/15
                                                       !|  |  |
          else  !---------------------------------------+  |  |
                                                       !|  |  |
c Non-diagonal matrix elements:                         |  |  |
c +-----------------------------------------------+     |  |  |
c |Calculate the Bragg angle; if sin(QB)>1, then: |     |  |  |
c |- in interactive mode take QB=90 and warn,     |     |  |  |
c |- in batch mode tak QB=0 and warn.             |     |  |  |
c +-----------------------------------------------+     |  |  |
            if (modebat.eq.0) Then !--+                 |  |  |
c Take sin(QB) = 1., QB=90. and warn: |                 |  |  |
              iwarn = 3              !|                 |  |  |
            else  !-------------------+                 |  |  |
c Take sin(QB) = 0., QB=0. and warn:  |                 |  |  |
              iwarn = 1              !|                 |  |  |
            endif !-------------------+                 |  |  |
c Do not calculate sin,cos,... (works faster):          |  |  |
            isico = 0                                  !|  |  |
            QBij = Bragan8 (wave,ind,inmax,iwarn,d8,   !|  |  |
     *                          isico,s,s,s,s,s,s)     !|  |  |
            d = Sngl(d8)                               !|  |  |
            if (abs(QBij).lt.1.e-30) Then !--------+    |  |  |
              ifail=5 !the length of Bragan message|    |  |  |
              goto 998                            !|    |  |  |
            endif !--------------------------------+    |  |  |
                                                       !|  |  |
c This is equivalent sin(QB)/lambda:                    |  |  |
            s = 1./(2.*d)                              !|  |  |
                                                       !|  |  |
            if (i.ne.incid .and. j.eq.incid)           !|  |  |
     *                              QB(i)=Sngl(QBij)   !|  |  |
                                                       !|  |  |
            xrhc = (0.,0.)                             !|  |  |
            xihd = (0.,0.)                             !|  |  |
            xihq = (0.,0.)                             !|  |  |
                                                       !|  |  |
c +----------------------------------+                  |  |  |
c | Internal loop for calculating    |                  |  |  |
c | xh by summation of contributions +---------+        |  |  |
c | of crystal components:           |         |        |  |  |
c +----------------------------------+         |        |  |  |
            do  indx=1,kcomp  !======================+  |  |  |
c +------------------------------------+             |  |  |  |
c | Calculate the Debye-Waller factor: |             |  |  |  |
c +------------------------------------+             |  |  |  |
              wh(indx) = exp(-bdw(indx)*s*s)        !|  |  |  |
                                                    !|  |  |  |
              if (ipr.eq.2)     Then  !--------+     |  |  |  |
                write (txt,410)  name(indx),  !|     |  |  |  |
     *                             wh(indx)   !|     v  v  v  v
  410           format('Element: ',a,
     *                 '    exp(-w)=',f7.4)   !|     ^  ^  ^  ^
                Call PriBlo (txt,1,1,0,0,ifa) !|     |  |  |  |
              endif  !-------------------------+     |  |  |  |
                                                    !|  |  |  |
c +-------------------------------------+            |  |  |  |
c | Calculate the scattering amplitude  |            |  |  |  |
c | of X-ray at the Bragg angle QB      |            |  |  |  |
c | for the i-th crystal element with   |            |  |  |  |
c | account for dispersion corrections: |            |  |  |  |
c +-------------------------------------+            |  |  |  |
              fh(indx) = Forfac(s,ifail)            !|  |  |  |
            if (ifail.ne.0)     goto 998            !|  |  |  |
                                                    !|  |  |  |
c +-------------------------------+                  |  |  |  |
c | Sum the exponents multiplying |                  |  |  |  |
c | the elements contributions:   |                  |  |  |  |
c +-------------------------------+                  |  |  |  |
              sumexp = (0.,0.)                      !|  |  |  |
              do        jndx=1,kol(indx)   !======+  |  |  |  |
                sum = ind(1)     * wx(indx,jndx) !|  |  |  |  |
     +              + ind(2)     * wy(indx,jndx) !|  |  |  |  |
     +              + ind(inmax) * wz(indx,jndx) !|  |  |  |  |
                sumexp = sumexp +                !|  |  |  |  |
     +                    exp(2.*pi*im*sum*ipm)  !|  |  |  |  |
              enddo  !=======================^====+  |  |  |  |
c ...................                        |       |  |  |  |
c  +-----------------------------------------+       |  |  |  |
c ipm=1 in the exponents correspond to the Pinsker   |  |  |  |
c notation (p.18,71). In this case the Bloch waves   |  |  |  |
c should contain exp(-i*k*r).                        |  |  |  |
c ipm=-1 in the exponents correspond to notation by  |  |  |  |
c Afanasiev and Baryshevsky. In this case the Bloch  |  |  |  |
c wave should  contain exp(i*k*r).                   |  |  |  |
c .................................................. |  |  |  |
              sumexp = wh(indx)*prcn(indx) * sumexp !|  |  |  |
              xrhc = xrhc + cr*fh(indx) * sumexp    !|  |  |  |
              xihd = xihd + ci*td(indx) * sumexp    !|  |  |  |
              xihq = xihq + ci*tq(indx) * sumexp    !|  |  |  |
            enddo   !================================+  |  |  |
c +-----------------------+  |                          |  |  |
c | End of summation over +--+                          |  |  |
c | crystal elements      |                             |  |  |
c +-----------------------+                             |  |  |
                                                       !|  |  |
c Modula of xrh,xih,xqh:                                |  |  |
            xrh(i,j) = Abs(xrhc)                       !|  |  |
            xih(i,j) = Abs(xihd)                       !|  |  |
            xqh(i,j) = Abs(xihq)                       !|  |  |
                                                       !|  |  |
c Phases of xrh,xih,xqh:                                |  |  |
            if (abs(xrhc).lt.1.e-30) Then !--------+    |  |  |
              xrf(i,j)=0.                         !|    |  |  |
            else  !--------------------------------+    |  |  |
              xrf(i,j) = atan2(aimag(xrhc),       !|    |  |  |
     *                          real(xrhc)) / pi  !|    |  |  |
            endif  !-------------------------------+    |  |  |
                                                       !|  |  |
            if (abs(xihd).lt.1.e-30) Then !--------+    |  |  |
              xif(i,j)=0.                         !|    |  |  |
            else   !-------------------------------+    |  |  |
              xif(i,j) = atan2(aimag(xihd),       !|    |  |  |
     *                          real(xihd)) / pi  !|    |  |  |
            endif  !-------------------------------+    |  |  |
                                                       !|  |  |
            if (abs(xihq).lt.1.e-30) Then !--------+    |  |  |
              xqf(i,j)=0.                         !|    |  |  |
            else   !-------------------------------+    |  |  |
              xqf(i,j) = atan2(aimag(xihq),       !|    |  |  |
     *                          real(xihq)) / pi  !|    |  |  |
            endif  !-------------------------------+    |  |  |
                                                       !|  |  |
            if (ipr.ne.0)   Then  !----------------+    |  |  |
              write (txt,18) QBij,i,j,            !|    |  |  |
     *               xrh(i,j),xrf(i,j),           !|    |  |  |
     *               xih(i,j),xif(i,j),           !|    |  |  |
     *               xqh(i,j),xqf(i,j)            !|    |  |  |
              Call PakInt (ind,inmax,txt(20),nb)  !|    |  |  |
              txt(2)(14:35)=txt(20)(1:nb)         !|    |  |  |
              Call PriBlo (txt,8,0,0,0,ifa)       !|    |  |  |
            endif  !-------------------------------+    |  |  |
          endif   !-------------------------------------+  |  |
        enddo  !===========================================+  |
        enddo  !==============================================+
                                               !|
c +------------------------------------------+  |
c | End of dual loop over Bragg reflections. +--+
c +------------------------------------------+

  18    format  (
     *  '+',34('-'),                       '-',24('-'),             '+'/
     *  '|Reflection:',23x,                '|Bragg angle:',g12.4,   '|'/
     *  '|',17('-'),      '-',16('-'),     '-----',20('-'),         '|'/
     *  '| {',i3,',',i3,'}',7x,
     *                    '|',7x,'Module',7x,  '|',6x,'Phase/pi',6x,'|'/
     *  '|Xrh',14x,       '|',g14.5,6x,        '|',g14.5,6x,        '|'/
     *  '|Xih-dipole       |',g14.5,6x,        '|',g14.5,6x,        '|'/
     *  '|Xih-quadrupole   |',g14.5,6x,        '|',g14.5,6x,        '|'/
     *  '+',17('-'),      '-',20('-'),         '-',20('-'),         '+')


  613   continue
        wv = wave
        do      i=1,6  !=====+
          ad(i) = a(i)      !|
        enddo  !=============+
        ising = isym

  998   continue
        if (ipr.eq.2 .AND. linesX0hWarn.gt.0)  Then  !-----------------+
          Call Priblo("X0h2: warnings summary for "//code,1,1,0,0,ifa)!|
          Call Priblo(txtX0hWarn,linesX0hWarn,0,1,0,ifa)              !|
        endif !--------------------------------------------------------+
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
