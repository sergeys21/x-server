c This module contains subroutines for BrLs_Surf:
c 1. Read_Input_File
c 2. Reflex_Sorting
c 3. Sort_by_Epsilon_In_Part (Imaginary Part)

c ====================================================
        Subroutine      Read_Input_File (luninp,luntbl,ifail)
c-------------------------------------------------------
c  This sub reads the BrLs_12 input file and performs
c  preparationally operations for BrLs_12
c
c
c               Author: S. Stepanov (1993-2003)
c-------------------------------------------------------
c The following include contain something like this:
c       Integer         nDimen
c       Parameter       (nDimen = 12)
c       Integer         nDimenTheta
c       Parameter       (nDimenTheta = 501)
c+===============================+
        Include 'brls_12.inc'   !|
c+===============================+

        Complex*16      x0

        Real*8          WaveLength, Wavelength_coplanar,
     *                  Incident_Polarization(3),
     *                  h_vectors(3,nDimen+1), h_module(nDimen+1),
     *                  k_vectors(3,nDimen+1), k_coplanar(3),
     *                  h1xh2(3), h1xh2_inv(3), h1xh2_module,
     *                  k_candidate(3), a1(3), a2(3),
     *                  Uncoplanar, DisorNorm, DisorNorm_gra,
     *                  k_quest(3), Wave_quest,
     *                  Interplane_Angles,
     *                  k_Angles(nDimen), h_Angles(nDimen),
     *                  Vertical(3), Beta(nDimen+1),
     *                  Extinction(nDimen+1), Extinction0, Ext_min,
     *                  Dpi, Dgra, angle8a, angle8b,
     *                  Coef, s8, s8_, x8, Demp8, exp8,
     *                  devmax/500./, b8, d8

        Real*8          DAngle2,DVecMod2,DVecSca2,wave2energy8,Bragan8
        External        DAngle2,DVecMod2,DVecSca2,wave2energy8,Bragan8

c       Real*8          X02AEF          !NAGlib routine
c       External        X02AEF          !Largest negative argument for EXP
c       Real*8          X02AFF          !NAGlib routine
c       External        X02AFF          !Largest positive argument for EXP

        Real            LatticeParameters(6), Wave, w,
     *                  x0r, x0i, mu0, rho,
     *                  xhr(nDimen,nDimen),
     *                  xhi(nDimen,nDimen),
     *                  xfr(nDimen,nDimen),
     *                  xfi(nDimen,nDimen),
     *                  Bragg_Angles(nDimen),
     *                  polang, rpar(2), s4,
     *                  Tetapi, Tetasg, Tetain

        Integer         luninp, luntbl, ifail, lpn,
     *                  ind(3,nDimen+1), iii(3), ipar(2),
     *                  nBaseNormal(3), nReperDirection(3),
     *                  h1xh2_indices(3), k0_indices(3),
     *                  i, j, l, m, n, jpr, nforbidden, nforbid_pi,
     *                  iTet, iSort, iSing, Mode, line,
     *                  i1, i2, i3, k0_choice, io_status

        Character       Code*20, Radiat*6,
     *                  cmTet(3)*32, wrk*80,
     *                  BraggWave(nDimen)*28,
     *                  eqormore*1, attn*17, weak*16

        Character       choice_lbl(4)*26 /'First k0 found',
     *                                    'Second k0 found',
     *                                    'Smaller angle with surface',
     *                                    'Larger angle with surface'/
        Character       Uni(5)*5        /'degr.','min.',
     *                                   'mrad.','sec.', 'urad.'/
c--------------------------------------------------------
        Complex*16      xh(nDimen,nDimen), Imaginary1
        Real*8          k_module,
     *                  ca1(nDimen),ca2(nDimen),
     *                  ca3(nDimen),ca4(nDimen),
     *                  k_unit(3,nDimen+1),
     *                  Sg_unit(3,nDimen),
     *                  Pi_unit(3,nDimen),
     *                  csg, cpi, cSurf1, cSurf2,
     *                  Gamma_Br(nDimen+1),
     *                  Total_Reflex_Angle
        Real            pi, xabs, Thickness,
     *                  xhq(nDimen,nDimen),
     *                  xfq(nDimen,nDimen),
     *                  Tet0(2), Tetk(2), dTet(2),
     *                  UnitData(5), UnitCoef(2)
        Logical         Surface_Case(nDimen)
        Integer         nReflexes, nReflexes2,
     *                  ipol, iout, ii, ipv,
     *                  nMaxPoints(2),
     *                  Surf_Ref_Num,
     *                  Laue_Ref_Num,
     *                  Bragg_Ref_Num,
     *                  Incident_Place,
     *                  iWarning, Task_Size,
     *                  New_ord(nDimen),
     *                  New_pos(nDimen),
     *                  linp, lout, iUnits(2),
     *                  iBatch_mode, ixway
        Character       Refl_Type(nDimen+1)*8,
     *                  ProgramVer*8,
     *                  InpFile*80, OutFile*80,
     *                  Comment(3)*80
        Common  /BrSurf/xh, Imaginary1,                 !complex*16
     *                  k_module,                       !real*8
     *                  ca1, ca2, ca3, ca4,             !real*8
     *                  k_unit, Sg_unit, Pi_unit,       !real*8
     *                  csg, cpi,                       !real*8
     *                  cSurf1, cSurf2,                 !real*8
     *                  Gamma_Br, Total_Reflex_Angle,   !real*8
     *                  Dpi, Dgra,                      !real*8
     *                  pi, xabs,                       !real
     *                  xhq, xfq,                       !real
     *                  Thickness,                      !real
     *                  Tet0, dTet, Tetk,               !real
     *                  UnitData, UnitCoef,             !real
     *                  Surface_Case,                   !Log
     *                  nReflexes, nReflexes2,          !int
     *                  ipol, iout,                     !int
     *                  ii, nMaxPoints,                 !int
     *                  Surf_Ref_Num,                   !int
     *                  Laue_Ref_Num,                   !int
     *                  Bragg_Ref_Num,                  !int
     *                  Incident_Place, iWarning,       !int
     *                  New_ord, New_pos,               !int
     *                  Task_Size, ipv,                 !int
     *                  linp, lout, ixway,              !int
     *                  iBatch_mode, iUnits,            !int
     *                  Refl_Type, ProgramVer,          !char*8
     *                  InpFile, OutFile, Comment       !char*80
c -------------------------------------------------------
        Real*8              Surface_Normal(3)
        Common /SurfNorm/   Surface_Normal
c -------------------------------------------------------
c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer             iHenkeCowan
        Common /HenkeCowan/ iHenkeCowan
c--------------------------------------------------------
        Character       txt(20)*80
        Common  /msg/   txt
c--------------------------------------------------------
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
          progname = 'Read_Input_File'    !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim (progname)
c-------------------------------------------------------
c Minimum (largest negative) argument for EXP:
        x8    = 0.0D0
c X02AEF Largest negative permissible argument for EXP
c       exp8  = X02AEF(x8)              !NAG2 function
c X02AFF Largest positive permissible argument for EXP
c       exp8  = X02AFF(x8)              !NAG2 function
        exp8  = 709.78271289338397      !value for GNU Fortran
c       exp8  = maxexponent(x8)         !GNU Fortran function gives "1024"
        ifail = 0
        cmTet(1) = 'Theta2 along sg(0) vector'
        cmTet(2) = 'Theta2 .ort. Surface Normal'
        cmTet(3) = 'Theta2 .ort. ( ??? ??? ???)'
c Number of Bragg planes indices (3 or 4):
        ii = 3
c The PI number:
        Dpi = (4.D0)*Atan2(1.D0,1.D0)
        pi  = Sngl(Dpi)
c Number of radians in 1 degree, 1 minute, and 1 second:
        Dgra = (2.D0)*Dpi/(360.D0)
        UnitData(1) = Sngl(Dgra)                        !degr
        UnitData(2) = Sngl(Dgra/(60.D0))                !arc min
        UnitData(3) = 1.E-03                            !mrad
        UnitData(4) = Sngl(Dgra/(3600.D0))              !arc sec
        UnitData(5) = 1.E-06                            !murad

c Imaginary 1 (i):
        Imaginary1 = (0.,1.)

c+-------------------------------------------------------+
c Open input file:                                       |
        linp = Max (Len_Trim(InpFile),1)                !|
        Call OpenFile(InpFile(1:linp),luninp,           !|
     *                'read','old',io_status,*99)       !v

        Call Make_Err_Filename (InpFile(1:linp))

        line = 0
        txt(7) = 'batch mode flag [0/1/2]'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i.ne.0)             goto 100
        iBatch_mode = ipar(1)
        if (iBatch_mode.lt.0)   goto 101
        if (iBatch_mode.gt.2)   goto 101
        if (iBatch_mode.eq.0) iBatch_mode=1
        modebat = iBatch_mode
c-----------------------
        txt(7) = 'output file name'
        Call    LineLoop (luninp,line,wrk,*100)
        OutFile = wrk
        lout = Len_Trim(OutFile)
        if (lout.eq.0)  Then  !-----+
          OutFile = '~'            !|
          lout    = 1              !|
        endif  !--------------------+
c Open output log file:
        Call    OpenList (OutFile(1:lout)//'.tbl',luntbl,ifail)
        if (ifail.ne.0) goto 28
c-----------------------
        do i=1,3  !=====================================+
          write (txt(7),'(a,i1)') 'comment line #',i   !|
          Call  LineLoop (luninp,line,Comment(i),*100) !|
        enddo  !========================================+
c-----------------------
        txt(7) = 'units for surface miscut and scan angles [0--4]'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (iUnits,2,wrk,i)
        if (i.ne.0)             goto 100
c Units for:
c iUnits(1) -- surface miscut angle
c iUnits(2) -- scan angles
        do i=1,2  !=============================+0=degr
          if (iUnits(i).lt.0)   goto 101       !|1=min
          if (iUnits(i).gt.4)   goto 101       !|2=mrad
          UnitCoef(i) = UnitData(iUnits(i)+1)  !|3=sec
        enddo  !================================+4=urad
c-----------------------
c Code - crystal code in the x0h DB:
        txt(7) = 'crystal code'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    txtShift (wrk,j,i)
        i = Min(Len(Code),Len(wrk))
        Code = wrk(1:i)
c-----------------------
c Thickness - crystal thickness in micrones:
        txt(7) = 'crystal thickness'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdReal  (rpar,1,wrk,ifail)
        if (ifail.ne.0)         goto 100
        Thickness = rpar(1)
        if (Thickness.le.0.)    goto 101
c-----------------------
c Read the 3 parameters required to determine the actual INTERNAL
c surface normal 'Surface_Normal':

c 1. nBaseNormal - base surface normal vector indices:
        txt(7) = 'base surface normal indices'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (nBaseNormal,3,wrk,ifail)
        if (ifail.ne.0)                 goto 100
        txt(7) = 'base surface normal indices are all zero'
        if ((nBaseNormal(1).eq.0) .AND.
     *      (nBaseNormal(2).eq.0) .AND.
     *      (nBaseNormal(3).eq.0))      goto 101

c 2. nReperDirection - indices of vector, pointing out
c    the direction of maximum misorientation of real
c    surface normal with respect to base normal vector:
        txt(7) = 'surface miscut direction indices'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (nReperDirection,3,wrk,ifail)
        if (ifail.ne.0)                 goto 100

c 3. DisorNorm - Maximum misorientation angle of real
c    surface normal:
        txt(7) = 'surface miscut angle [-45:45] degr.'
        Call    LineLoop (luninp,line,wrk,*100)
c ATTENTION: this is a double-precision READING!
        Call    rdReal8 (DisorNorm,1,wrk,ifail)
        if (ifail.ne.0)                 goto 100
        DisorNorm_gra = DisorNorm * UnitCoef(1) / Dgra
        if ((DisorNorm_gra .lt. -45.) .OR.
     *      (DisorNorm_gra .gt. +45.))  goto 101

        txt(7) = 'surface miscut direction indices are all zero'
        if ((abs(DisorNorm) .gt. 1.0D-10) .AND.
     *      (nReperDirection(1).eq.0)     .AND.
     *      (nReperDirection(2).eq.0)     .AND.
     *      (nReperDirection(3).eq.0))  goto 101
c-----------------------
        txt(7) = 'x-ray wave data (1=wave 2=energy 3=line 4=copl. 5=h3)'
        Call    LineLoop (luninp,line,wrk,*100)
c 1=wave 2=energy 3=x-ray line 4=coplanar case 5=reflex-3
        Call    rdInt   (ipar,1,wrk,i)
        if (i.ne.0)             goto 100
        ixway = ipar(1)
        if (ixway.lt.1)         goto 101
        if (ixway.gt.5)         goto 101
c-----------------------
        txt(7) = 'x-ray wavelength/energy (range=[0.1-10.0]Angstrom)'
        Call    LineLoop (luninp,line,wrk,*100)
        WaveLength = 0.
        if (ixway.le.2) Then  !--------------------------------+
          Call  rdReal8 (WaveLength,1,wrk,i)                  !|
          if (i.ne.0)           goto 100                      !|
          if (WaveLength.le.0.) goto 101                      !|
          if (ixway.eq.2) Then  !------------------+           |
            WaveLength = wave2energy8(WaveLength) !|           |
          endif  !---------------------------------+           |
          if (WaveLength.lt.0.1)    goto 101                  !|
          if (WaveLength.gt.10.)    goto 101                  !|
        endif  !-----------------------------------------------+
        Wave = Sngl(WaveLength)
c-----------------------
        txt(7) = 'x-ray line name'
        Call    LineLoop (luninp,line,wrk,*100)
        if (ixway.eq.3) Then  !-----------+
          i = Min(Len(Radiat),Len(wrk))  !|
          Radiat = wrk(1:i)              !|
        else !----------------------------+
          Radiat = ' '                   !|
        endif  !--------------------------+
        if     (ixway.le.3) Then  !---------------------------------+
          Call Get_Wave (Wave,Radiat,ifail)            !Get Wave!   |
          if (ifail.ne.0) goto 130                                 !|
          WaveLength = Dble(Wave)                                  !|
        elseif (ixway.eq.4) Then  !---------------------------------+
c Wave to be determined in InciVec from the coplanar condition:     |
          Wave = -1.                                               !|
          WaveLength = Dble(Wave)                                  !|
        elseif (ixway.eq.5) Then  !---------------------------------+
c Wave to be determined in InciVec from the Bragg condition for h3: |
          Wave = -2.                                               !|
          WaveLength = Dble(Wave)                                  !|
        endif  !----------------------------------------------------+
c-----------------------
c+=======================================================+
c|ipol - incident wave polarization state                |
c|0 - incident wave in unpolarized,                      |
c|1 - incident wave is polarized, angle will be asked    |
c|2 - incident wave is polarized, E0.ort.k0 & E0.ort.k1  |
c|    (emulation of the Victor Kohn paper in             |
c|     Sov.Phys.Crystallography, 1975).                  |
c+=======================================================+
        txt(7) = 'x-ray polarization(0-both, 1-angle to pi0, 2~[k0*k1])'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i.ne.0)             goto 100
        iPol = ipar(1)
c Polarization: 0-both, 1-angle to pi0, 2~[k0*k1]:
        if (iPol.lt.0)  goto 101
        if (iPol.gt.2)  goto 101
c-----------------------
        txt(7) = 'angle of incident polarization plane to Pi(0) plane'
        Call    LineLoop (luninp,line,wrk,*100)
        polang = 0.
        if (ipol.eq.1) Then  !---------------------------------+
          Call  rdReal  (rpar,1,wrk,i)                        !|
          if (i.ne.0)           goto 100                      !|
          polang = rpar(1)                                    !|
          if (polang.lt.-360.)      goto 101                  !|
          if (polang.gt. 360.)      goto 101                  !|
        endif  !-----------------------------------------------+
c-----------------------
        txt(7) = 'external X0h database flag [-1:4]'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i .ne. 0)           goto 100
        iHenkeCowan = ipar(1)
        if (iHenkeCowan .lt. -1) goto 101
        if (iHenkeCowan .gt.  8) goto 101
c-----------------------
c ind - Bragg reflection indices:
        do      j=1,ii  !==+ The first element in the array of reflections
          ind(j,1) = 0    !| is always (0,0,0)
        enddo  !===========+
        nReflexes = 1
        j = 0                   ! stop flag
        do  i=2,nDimen  !==========================+
          txt(7) = 'Indices of Bragg Reflections' !|
          Call LineLoop (luninp,line,wrk,*100)    !|
          Call rdInt (ind(1,i),ii,wrk,ifail)      !|
          if (ifail .ne. 0)     goto 100          !|
          if (j.eq.0) Then !-------------------+   |
            if (ind(1,i) .ne. 0  .or.         !|   |
     *          ind(2,i) .ne. 0  .or.         !|   |
     *          ind(ii,i).ne. 0) Then  !---+   |   |
              nReflexes = nReflexes+1     !|   |   |
            else !-------------------------+   |   |
              j = 1   !end of reflections !|   |   |
            endif  !-----------------------+   |   |
          else !-------------------------------|   |
            ind(1,i)  = 0                     !|   |
            ind(2,i)  = 0                     !|   |
            ind(ii,i) = 0                     !|   |
          endif  !-----------------------------+   |
        enddo  !===================================+
c nReflexes - number of Bragg reflections:
        if (ixway.ne.5) then !-----------------------------+
          txt(7) = 'number of Reflections (more than 2)'  !|
          if (nReflexes.lt.3) goto 101                    !|
        else  !--------------------------------------------|
          txt(7) = 'number of Reflections (more than 3)'  !|
          if (nReflexes.lt.4) goto 101                    !|
        endif !--------------------------------------------+
        nReflexes2 = 2*nReflexes
c-----------------------
c nMaxPoints - number of scan points over Theta1, Theta2:
        write(txt(7),'(a,i5,a)')
     *                       'number of scan points [1-',nDimenTheta,']'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (nMaxPoints,2,wrk,i)
        if (i .ne. 0)           goto 100
        do i=1,2 !====================================+
          if (nMaxPoints(i).lt.     1     ) goto 101 !|
          if (nMaxPoints(i).gt.nDimenTheta) goto 101 !|
        enddo !=======================================+
c-----------------------
c Tet0 - scan starts for Theta1 and Theta2:
        txt(7) = 'start scan angles [-5:5] degr.'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdReal  (Tet0,2,wrk,i)
        if (i .ne. 0)           goto 100
        do i=1,2 !====================================+
          x8 = Tet0(i) * UnitCoef(2) / Dgra          !|
          if (x8 .lt. -5.) goto 101                  !|
          if (x8 .gt.  5.) goto 101                  !|
        enddo !=======================================+
c-----------------------
c Tetk - scan ends for Theta1, Theta2:
        txt(7) = 'end scan angles [-5:5] degr.'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdReal  (Tetk,2,wrk,i)
        if (i .ne. 0)           goto 100
        do i=1,2 !=======================================+
          x8 = Tetk(i) * UnitCoef(2) / Dgra             !|
          if (x8 .lt. -5.) goto 101                     !|
          if (x8 .gt.  5.) goto 101                     !|
c Added 2005/05/23:                                      |
c Changed to error on 06/06/16:                          |
c (the calling program like brl_frm cannot               |
c know that we are no longer producing GRD):             |
          if (abs(Tetk(i)-Tet0(i)).lt.1.E-10) Then !--+  |if unequal
            if (nMaxPoints(i) .ne. 1) goto 109 !------+--+-----+
          endif !-------------------------------------+  |     v
        enddo !==========================================+

c+======================================================+
c|iTet - flag defining the choice of vectors a1 and a2, |
c        along which Theta1 and Theta2 are measured     |
c|0- Normal choice: a1=pi0=pp(1),  a2=sg0=sg(1).        |
c|1- (a2.ort.norm & a2.ort.k0), (a1.ort.a2 & a1.ort.k0).|
c|2- (a2.ort.k1 & a2.ort.k0), (a1.ort.a2 & a1.ort.k0)   |
c|   (emulation of the V.Kohn paper in Sov. Phys.       |
c|   Crystallography, 1975).                            |
c|   If ipol=2, then iTet=2 is chosen automatically.    |
c+======================================================+
        txt(7) = 'flag choosing scan axes direction [0-2]'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i.ne.0)             goto 100
        iTet = ipar(1)
        if (iTet.lt.0)          goto 101
        if (iTet.gt.2)          goto 101
c-----------------------
c iWarning - warning flag (0,1,2,3):
        txt(7) = 'warning level flag [0-3]'
        Call    LineLoop (luninp,line,wrk,*100)
        Call    rdInt   (ipar,1,wrk,i)
        if (i.ne.0)             goto 100
        iWarning = ipar(1)
        if (iWarning.lt.0)      goto 101
        if (iWarning.gt.3)      goto 101
c-----------------------
c Optional "k0_choice" flag (can be 1,2,3,4):
c This flag directs BRL what it should do if two k0 are found,
c both satisfying the multiple diffraction condition and both
c being able to enter the crystal (i.e. making obtuse angles
c with external surface normal)
c 0 = take the first found (the only behavior before 2/2014)
c 1 = take the first found (the only behavior before 2/2014)
c 2 = take the second found
c 3 = take the one with the smallest |gamma0| (smallest angle to the surface)
c 4 = take the one with the largest |gamma0| (largest angle to the surface)
        txt(7) = 'Multiple k0 choice flag [0-2]'
        k0_choice = 0                                   !old-style default
        Call    LineLoop (luninp,line,wrk,*10)
        Call    rdInt   (ipar,1,wrk,i)
        if (i.ne.0)             goto 100
        k0_choice = ipar(1)
        if (k0_choice.eq.0) k0_choice = 1               !default
        if (k0_choice.lt.1)      goto 101
        if (k0_choice.gt.4)      goto 101
  10    continue
                                                        !^
        Close   (unit=luninp)                           !|
c+-------------------------------------------------------+
c-------------------------------------------------------
c Bragg reflections order (may change later if the reflections get sorted):
        do      i=1,nReflexes !==========+
          New_ord(i) = i                !|
          New_pos(i) = New_ord(i)       !|
        enddo !==========================+
        Incident_Place = New_pos(1)
c-------------------------------------------------------
c Call X0h to calculate LatticeParameters(1-6) and x0 = cmplx(x0r,x0i):
        do i=1,6 !=====================+
           LatticeParameters(i) = 0.  !|
        enddo !========================+
        do i=1,nReflexes !=========+
          Bragg_Angles(i) = 0.    !|
          do j=1,nReflexes !=====+ |
            xhr(i,j) = 0.       !| |
            xhi(i,j) = 0.       !| |
            xfr(i,j) = 0.       !| |
            xfi(i,j) = 0.       !| |
          enddo !================+ |
        enddo !====================+
        iSing = 0
        Mode  = 0       !calling X0h for the first time
        jpr   = 0
        if (ixway.le.3) Then !--------------+
          w = Wave                         !|
        else !------------------------------|
          w = 1.   !just put something!    !|
        endif !-----------------------------+
        ifail = 0
        rho   = 0.
        Call X0h1 (LatticeParameters, w,
     *             Code, rho, 0,
     *             0, 0, 0,
     *             iSing, Mode,
     *             jpr, Bragg_Angles, x0r, x0i,
     *             xhr, xhi, xhq,
     *             xfr, xfi, xfq, ifail)
        if (ifail.ne.0) goto 130
        if (ixway.eq.3 .AND.
     *      abs(Wave).lt.1.E-10) Then !----+
          Wave       = w      !x-ray line  |
          Wavelength = w      !x-ray line  |
        endif !----------------------------+

c-------------------------------------------------------
c Determine INTERNAL normal to crystal surface:
        Call    SurfNormVect    (nBaseNormal,
     *                           nReperDirection,
     *                           DisorNorm,
     *                           Surface_Normal,
     *                           ifail)
        if (ifail.ne.0) goto 28
c-------------------------------------------------------
c Find k0, the wavevector of the incident wave:
c (see incivect.for)
        Call IncidentVector (WaveLength, ii,
     *                       nReflexes+1, ind,
     *                       h_vectors, h_module,
     *                       h1xh2, h1xh2_module,
     *                       k_coplanar,
     *                       Wavelength_coplanar,
     *                       Uncoplanar,
     *                       k_vectors(1,Incident_Place),
     *                       k_module, ifail)
        if (ifail.ne.0) goto 130

        do      i=1,ii  !===============================+
          if (h1xh2(i) .ge. 0.) Then !---------------+  |
             h1xh2_indices(i) = INT(h1xh2(i)+0.5)   !|  |
          else !-------------------------------------+  |
             h1xh2_indices(i) = INT(h1xh2(i)-0.5)   !|  |
          endif !------------------------------------+  |
        enddo !=========================================+

        if (ixway .eq. 4 .OR. ixway .eq. 5) Wave = Sngl(Wavelength)
c Now we must know the wavelength for sure!
        Call X0h1 (LatticeParameters, Wave,
     *             Code, rho, 0,
     *             0, 0, 0,
     *             iSing, Mode,
     *             jpr, Bragg_Angles, x0r, x0i,
     *             xhr, xhi, xhq,
     *             xfr, xfi, xfq, ifail)
        if (ifail.ne.0) goto 130
        Mode  = 5       !further X0h calls will not be reading X0h DB

        xabs  = Sqrt(Abs(x0r**2+x0i**2))
c Critical angle of Total External Reflection (TER) in radians:
        Total_Reflex_Angle = Sqrt(xabs)
        Extinction0 = (1.E-4)*Wavelength/(2.*pi*Total_Reflex_Angle)
c-------------------------------------------------------

c+=====================================================+
c|Flag for output files: '.grd'  or '.dat':            |
c|iout=0 - 2D scanning, write .grd files,              |
c|iout=1 - 1D scanning over Theta1, write .dat files,  |
c|iout=2 - 1D scanning over Theta2, write .dat files.  |
c+=====================================================+
        iout=0
        if (nMaxPoints(2).le.1) iout=1
        if (nMaxPoints(1).le.1) iout=2

        do      i=1,2  !=============================+
          if (nMaxPoints(i).gt.1)       Then  !---+  |
            dTet(i) = (Tetk(i) - Tet0(i)) /      !|  |
     /                (nMaxPoints(i) - 1)        !|  |
          else  !---------------------------------|  |
            dTet(i) = 0.                         !|  |
          endif  !--------------------------------+  |
        enddo  !=====================================+
c-------------------------------------------------------
        isort = 0
c-------------------------------------------------------
c Print input parameters from 'BrL_12.inp':
        write (luntbl,401,err=301) Code, Thickness, Wavelength
        write (luntbl,433,err=301) (Uni(iUnits(2)+1),nMaxPoints(i),
     *                             dTet(i), Tet0(i), Tetk(i),i=1,2)
        write (luntbl,431,err=301) (nBaseNormal(i),i=1,ii),
     *                             DisorNorm,
     *                             (nReperDirection(i),i=1,ii)
        l = Len_Trim(choice_lbl(k0_choice))
        write (luntbl,435,err=301) choice_lbl(k0_choice)(1:l)
c Convolution of INTERNAL surface normal with k0:
        s8 = DVecSca2 (k_vectors(1,Incident_Place),
     *                 Surface_Normal, ii)
        if (s8.lt.0.)   Then   !--------------------------------+
          write (luntbl,441,err=301)                           !|
  441     format(' WARNING: Inverting [H1*H2]')                !|
          Call DVecCon (h1xh2, -1.D0, h1xh2, ii)               !|
c Recaclulate k0 and its convolution                            |
c with INTERNAL surface normal:                                 |
          Call DVecSum2 (k_coplanar, 1.D0,                     !|
     *                   h1xh2,      Uncoplanar,               !|
     *                   k_vectors(1,Incident_Place),          !|
     *                   ii)                                   !|
          s8 = DVecSca2 (k_vectors(1,Incident_Place),          !|
     *                   Surface_Normal, ii)                   !|
c If did not help, restore [h1xh2] and revert surface normal:   |
          if (s8.lt.0.)  Then  !----------------------------+   |
            write (luntbl,442,err=301)                     !|   |
  442       format(' WARNING: Inverting Surface Normal')   !|   |
            Call DVecCon (h1xh2, -1.D0, h1xh2, ii)         !|   |
            s8 = -s8                                       !|   |
            do  i=1,ii   !==============================+   |   |
              nBaseNormal(i)     = -nBaseNormal(i)     !|   |   |
              nReperDirection(i) = -nReperDirection(i) !|   |   |
              Surface_Normal(i)  =- Surface_Normal(i)  !|   |   |
            enddo   !===================================+   |   |
            Call DVecSum2(k_coplanar, 1.D0,                !|   |
     *                    h1xh2,      Uncoplanar,          !|   |
     *                    k_vectors(1,Incident_Place),     !|   |
     *                    ii)                              !|   |
          endif !-------------------------------------------+   |
        endif !-------------------------------------------------+

c If non-coplanar case:
        if (s8.gt.0.) Then !-----------------------------------------+
c Check if we possibly have two solutions for k0:                    |
          Call DVecCon (h1xh2, -1.D0, h1xh2_inv, ii)                !|
c Another k0 candidate:                                              |
          Call DVecSum2 (k_coplanar, 1.D0,                          !|
     *                   h1xh2_inv,  Uncoplanar,                    !|
     *                   k_candidate, ii)                           !|
          s8_ = DVecSca2 (k_candidate, Surface_Normal, ii)          !|
          if (s8_.gt.0) Then !------------------------------------+  |
c Angles of the two k0 candidates with the surface:               |  |
c           write (0,*) ' Read_Input_File: Surf.normal = (',     !|  |
c    *                                     Surface_Normal,')'    !|  |
c           write (0,*) ' Read_Input_File: k_vector = (',        !|  |
c    *                  (k_vectors(i,Incident_Place),i=1,ii),')' !|  |
c           write (0,*) ' Read_Input_File: k_candidate = (',     !|  |
c    *                           k_candidate,')'                 !|  |
            angle8a = 90. - DAngle2(Surface_Normal,              !|  |
     *                        k_vectors(1,Incident_Place), ii)   !|  |
            angle8b = 90. - DAngle2(Surface_Normal,              !|  |
     *                        k_candidate, ii)                   !|  |
            write (luntbl,443,err=301)                           !|  |
     *                   (k_vectors(i,Incident_Place),i=1,ii),   !|  |
     *                    angle8a,                               !|  |
     *                   (k_candidate(i),i=1,ii),                !|  |
     *                    angle8b                                !|  |
  443       format(
     *      '; WARNING: there are two incident wave vectors,',
     *                           1x,'both directed in such a'/
     *      '; way that they can enter the crystal and both',
     *                           1x,'satisfying the multiple'/
     *      '; diffraction condition:'/
     *      ';   k0_1=(',3g10.3,'), angle_with_surface:',
     *                                           f7.2, ' degr'/
     *      ';   k0_2=(',3g10.3,'), angle_with_surface:',
     *                                           f7.2, ' degr')  !|  |
            if     (k0_choice.eq.1) Then !---------------------+  |  |
c 1 = take the first one found -- i.e. do nothing              |  |  |
              write (luntbl,444,err=301) 'k0_1',              !|  |  |
     *        'requsted in the BRL input file'                !|  |  |
            elseif (k0_choice.eq.2) Then !---------------------+  |  |
c 2 = take the second one found                                |  |  |
              write (luntbl,444,err=301) 'k0_2',              !|  |  |
     *        'requsted in the BRL input file'                !|  |  |
              do  i=1,ii   !=====================+             |  |  |
                k_vectors(i,Incident_Place) =   !|             |  |  |
     =                          k_candidate(i)  !|             |  |  |
                h1xh2(i) = h1xh2_inv(i)         !|             |  |  |
              enddo   !==========================+             |  |  |
            elseif (k0_choice.eq.3) Then !---------------------+  |  |
c 3 = take the one with the smallest angle to the surface      |  |  |
              if (abs(angle8b).lt.abs(angle8a)) Then !------+  |  |  |
                write (luntbl,444,err=301) 'k0_2',         !|  |  |  |
     *          'it makes smaller angle with the surface'  !|  |  |  |
                do  i=1,ii   !=====================+        |  |  |  |
                  k_vectors(i,Incident_Place) =   !|        |  |  |  |
     =                            k_candidate(i)  !|        |  |  |  |
                  h1xh2(i) = h1xh2_inv(i)         !|        |  |  |  |
                enddo   !==========================+        |  |  |  |
            else !------------------------------------------+  |  |  |
                write (luntbl,444,err=301) 'k0_1',         !|  |  |  |
     *          'it makes smaller angle with the surface'  !|  |  |  |
            endif !-----------------------------------------+  |  |  |
            elseif (k0_choice.eq.4) Then !---------------------+  |  |
c 4 = take the one with the largest angle to the surface       |  |  |
              if (abs(angle8b).lt.abs(angle8a)) Then !------+  |  |  |
                write (luntbl,444,err=301) 'k0_2',         !|  |  |  |
     *          'it makes larger angle with the surface'   !|  |  |  |
                do  i=1,ii   !=====================+        |  |  |  |
                  k_vectors(i,Incident_Place) =   !|        |  |  |  |
     =                            k_candidate(i)  !|        |  |  |  |
                  h1xh2(i) = h1xh2_inv(i)         !|        |  |  |  |
                enddo   !==========================+        |  |  |  |
            else !------------------------------------------+  |  |  |
                write (luntbl,444,err=301) 'k0_1',         !|  |  |  |
     *          'it makes larger angle with the surface'   !|  |  |  |
            endif !-----------------------------------------+  |  |  |
            else !---------------------------------------------+  |  |
c 0 = take the first found (the only behavior before 2/2014)   |  |  |
              write (luntbl,444,err=301) 'k0_1',              !|  |  |
     *        'the choice is not given in the BRL input file' !|  |  |
            endif !--------------------------------------------+  |  |
          endif !-------------------------------------------------+  |
        endif !------------------------------------------------------+
  444   format(';'/'; Choosing ',a,' as ',a/';')

        s8 = DVecMod2(k_vectors(1,Incident_Place),ii)
        do      i=1,ii  !============================+
          x8 = 100.*k_vectors(i,Incident_Place)/s8  !|
          if (x8 .ge. 0.) Then !-----------+         |
            k0_indices(i) = INT(x8+0.5)   !|         |
          else !---------------------------+         |
            k0_indices(i) = INT(x8-0.5)   !|         |
          endif !--------------------------+         |
        enddo  !=====================================+

        write  (luntbl,58,err=301) k0_indices
  58    format ('; Incident wave vector k0 indices: (',3i5,')'/
     *  ';')

c Commented out by Sergey 2003/07/22:
c       if (ipol.eq.2 .and. iTet.ne.2)  Then  !--------+
c         iTet = 2                                    !|
c         write (luntbl,217,err=301)                  !v
c 217     format('; WARNING: switch to iTet=2 --',
c    *           ' a2.ort.k1 & a2.ort.k0')            !^
c       endif  !---------------------------------------+

c-------------------------------------------------------
c Here starts the second pass after sorting Bragg reflections:

  505   continue     !<------------------------------------------------+
                                                                      !|
        Laue_Ref_Num  = 0                                             !|
        Surf_Ref_Num  = 0                                             !|
        Bragg_Ref_Num = 0                                             !|
        Ext_min       = 1.E10                                         !|
                                                                      !|
        do i=1,nReflexes !==========================================+  |
c Exact wave  vectors of diffracted waves:                          |  ^
          Call DVecSum2 (k_vectors(1,Incident_Place),              !|  |
     *                   1.D0,                                     !|  |
     *                   h_vectors(1,i),                           !|  |
     *                   1.D0,                                     !|  |
     *                   k_vectors(1,i), ii)                       !|  |
                                                                   !|  |
c Deviation of incident vector from the given Bragg condition:      |  |
c ((k0+h)^2-k0^2)/k0^2 in the inits of |x0]:                        |  |
          ca4(i) = ((DVecMod2(k_vectors(1,i),ii)                   !|  |
     /           /   k_module)**2 - 1.D0) / xabs                   !|  |
                                                                   !|  |
c If the reflection is not in the exact Bragg condition:            |  |
          if (Dabs(ca4(i)) .gt. 0.01 .AND.                         !|  |
     *                  i  .ne. Incident_Place) Then  !---------+   |  |
                                                               !|   |  |
c           write (0,*) ' Read_Input_File:  ',                 !|   |  |
c    *                             'h_vectors(',i,') = (',     !|   |  |
c    *                             (h_vectors(j,i),j=1,ii),')' !|   |  |
c           write (0,*) ' Read_Input_File:  h1xh2 = (',        !|   |  |
c    *                                      h1xh2,')'          !|   |  |
            if (h_module(i) .lt. 1E-32) goto 104 !--------------+-+ |  |
            s8 = DAngle2 (h_vectors(1,i),h1xh2, ii)            !| v |  |
                                                               !|   |  |
            if ( (s8.gt.1.) .and. (s8.lt.179.)) Then !-------+  |   |  |
                                                            !|  |   |  |
c The reciprocal lattice vector is NOT coplanar to h1,h2:    |  |   |  |
                                                            !|  |   |  |
c+==========================================+                |  |   |  |
c|Calculate x-ray wavelength at which this  |                |  |   |  |
c|reflection would be in exact Bragg cond.  |                |  |   |  |
c|  From the conditions:                    |                |  |   |  |
c| ->    ->              ->                 |                |  |   |  |
c| k = k_copl + Uncopl*h1xh2 ,              |                |  |   |  |
c|                                          |                |  |   |  |
c|  -> ->    2                              |                |  |   |  |
c| 2(k*h) + h  = 0 ,                        |                |  |   |  |
c|                                          |                |  |   |  |
c|-- we have:                               |                |  |   |  |
c|            2        ->   ->       ->  -> |                |  |   |  |
c| Uncop = -{h  + 2(k_copl*h)} / 2(h1xh2*h) |                |  |   |  |
c+==========================================+                |  |   |  ^
              Uncoplanar=-0.5*(h_module(i)**2               !|  |   |  |
     +           + 2.*DVecSca2(k_coplanar,                  !|  |   |  |
     *                         h_vectors(1,i),              !|  |   |  |
     *                         ii)  )                       !|  |   |  |
     /              / DVecSca2(h1xh2,                       !|  |   |  |
     *                         h_vectors(1,i),              !|  |   |  |
     *                         ii)                          !|  |   |  |
              Call  DVecSum2 (k_coplanar,1.D0,              !|  |   |  |
     *                        h1xh2,Uncoplanar,             !|  |   |  |
     *                        k_quest,ii)                   !|  |   |  |
              Wave_quest = 2.*Dpi                           !|  |   |  |
     /                   / DVecMod2(k_quest,ii)             !|  |   |  |
c----------                                                  |  |   |  |
              if (Wave_quest                                !|  |   |  |
     *                .le.                                  !|  |   |  |
     *        Wavelength_coplanar)  Then  !--------------+   |  |   |  |
                                                        !|   |  |   |  |
                write (BraggWave(i),'(g18.10,1x,a)')    !|   |  |   |  |
     *                           Wave_quest, 'Angstrom' !|   |  |   |  |
                                                        !|   |  |   |  |
              else  !------------------------------------+   |  |   |  ^
                                                        !|   |  |   |  |
                BraggWave(i) = 'NOT EXIST, '//          !|   |  |   |  |
     *                              'Non-coplanar case' !|   |  |   |  |
                                                        !|   |  |   |  |
              endif  !-----------------------------------+   |  |   |  |
                                                            !|  |   |  |
            else  !------------------------------------------+  |   |  |
                                                            !|  |   |  |
c The reciprocal lattice vector IS coplanar to h1,h2:        |  |   |  |
                                                            !|  |   |  |
              BraggWave(i) = 'NOT EXIST; Coplanar case'     !|  |   |  |
                                                            !|  |   |  |
            endif  !-----------------------------------------+  |   |  |
                                                               !|   |  |
          else  !-----------------------------------------------+   |  |
                                                               !|   |  |
c The reflection IS in the exact Bragg condition:               |   |  ^
                                                               !|   |  |
            BraggWave(i) = 'THIS wavelength'                   !|   |  |
                                                               !|   |  |
          endif  !----------------------------------------------+   |  |
                                                                   !|  |
c Unit vectors along incident and diffracted waves:                !|  |
          Call DUnitVec2 (k_vectors(1,i),                          !|  |
     *                    k_unit(1,i),ii)                          !|  |
                                                                   !|  |
c Calculate cosines of the wavevectors with surface normal:        !|  |
          Gamma_Br(i) = DVecSca2 (k_unit(1,i),                     !|  |
     *                            Surface_Normal,ii)               !|  |
                                                                   !|  |
c Asymmetry factors:                                                |  |
          if (abs(Gamma_Br(i)) .gt. 1.D-5) Then !---------------+   |  |
            Beta(i) = Gamma_Br(Incident_Place)                 !|   |  |
     /              / Gamma_Br(i)                              !|   |  |
          else   !----------------------------------------------+   |  ^
            x8 = Abs(Gamma_Br(Incident_Place))                 !|   |  |
     -        - Abs(Gamma_Br(i))                               !|   |  |
            if     (x8.lt.0.D0) Then  !----------------------+  |   |  |
              Beta(i) = 0.D0                                !|  |   |  |
            elseif (abs(x8).lt.1.0D-10) Then  !--------------+  |   |  |
              Beta(i) = 1.D0                                !|  |   |  |
     *                * sign(1.D0,Gamma_Br(Incident_Place)) !|  |   |  |
     *                * sign(1.D0,Gamma_Br(i))              !|  |   |  |
            elseif (x8.gt.0.D0) Then  !----------------------+  |   |  |
              Beta(i) = (1.D+32)                            !|  |   |  |
     *                * sign(1.D0,Gamma_Br(Incident_Place)) !|  |   |  |
     *                * sign(1.D0,Gamma_Br(i))              !|  |   |  |
            endif  !-----------------------------------------+  |   |  |
          endif   !---------------------------------------------+   |  |
                                                                   !|  |
c Calculate extinction length in microns:                           |  |
          if (i.ne.Incident_Place) Then !-----------------------+   |  |
            xhr(1,1) = 0.                                      !|   |  |
            xhi(1,1) = 0.                                      !|   |  |
            xhq(1,1) = 0.                                      !|   |  |
            xfr(1,1) = 0.                                      !|   |  |
            xfi(1,1) = 0.                                      !|   |  |
            xfq(1,1) = 0.                                      !|   |  |
            Bragg_Angles(1) = 0.                               !|   |  |
            Call X0h1 (LatticeParameters, Wave,                !|   |  |
     *             Code, rho, 1,                               !|   |  |
     *             ind(1,i), ind(2,i), ind(ii,i),              !|   |  |
     *             iSing, Mode,                                !|   |  |
     *             jpr, Bragg_Angles, x0r, x0i,                !|   |  |
     *             xhr, xhi, xhq,                              !|   |  |
     *             xfr, xfi, xfq, ifail)                       !|   |  |
            Extinction(i) = (1.E-4)*WaveLength                 !|   |  |
     *       * sqrt(Gamma_Br(Incident_Place)*abs(Gamma_Br(i))) !|   |  |
     /       / abs(xhr(1,1))                                   !|   |  |
            if (Gamma_Br(i).lt.0.) Then !-------------+         |   |  |
c For the Bragg case definition differs: divide by pi |         |   |  |
               Extinction(i) = Extinction(i) / pi    !|         |   |  |
            endif !-----------------------------------+         |   |  |
          else !------------------------------------------------+   |  |
c For reflection-0 it is a projection of the absorption length: |   |  |
            Extinction(i) = (1.E-4)*WaveLength                 !|   |  |
     *                    * Gamma_Br(Incident_Place)           !|   |  |
     /                    / (2.*pi*abs(x0i))                   !|   |  |
          endif !-----------------------------------------------+   |  |
c This is the minimum extinction in the                             |  |
c case of zero-angle incidence or exit:                             |  |
          if (Extinction(i).lt.Extinction0)                        !|  |
     *                             Extinction(i) = Extinction0     !|  |
          if (Ext_min.gt.Extinction(i)) Ext_min = Extinction(i)    !|  !
                                                                   !|  |
c Determine reflection types:                                      !|  |
          if    (Gamma_Br(i) .gt. 10.*Total_Reflex_Angle) Then !+   |  |
c Laue:                                                         |   |  |
            Refl_Type(i)    = 'Laue'                           !|   |  |
            Laue_Ref_Num    = Laue_Ref_Num + 1                 !|   |  |
            Surface_Case(i) = .False.                          !|   |  |
                                                               !|   |  |
          elseif(Gamma_Br(i) .lt.-10.*Total_Reflex_Angle) Then !|   |  |
c Bragg:                                                        |   |  |
            Refl_Type(i)    = 'Bragg'                          !|   |  |
            Bragg_Ref_Num   = Bragg_Ref_Num + 1                !|   |  |
            Surface_Case(i) = .False.                          !|   |  |
                                                               !|   |  |
          else   !----------------------------------------------|   |  |
c Surface:                                                      |   |  |
            Refl_Type (i)   = 'Surface'                        !|   |  |
            Surf_Ref_Num    = Surf_Ref_Num + 1                 !|   |  |
            Surface_Case(i) = .True.                           !|   |  |
          endif  !----------------------------------------------+   |  ^
                                                                   !|  |
          if (isort .eq. 0)  Then  !----------------------------+   |  |
c+=================================+                            |   |  |
c|Print indices of all reflections |                            |   |  |
c|with respective type (Bragg/Laue)|                            |   |  |
c+=================================+                            |   |  |
            l = Max0(Len_Trim(Refl_Type(i)),1)                 !|   |  |
            if (i.le.3 .or. Abs(ca4(i)).lt.1E-3) Then !-----+   |   |  |
              write (luntbl,113,err=301) i-1,              !|   |   |  |
     *                                  (ind(j,i),j=1,ii), !|   |   |  |
     *                                   Refl_Type(i)(1:l) !|   |   |  |
            else !------------------------------------------+   |   |  |
              write (luntbl,113,err=301) i-1,              !|   |   |  |
     *                                  (ind(j,i),j=1,ii), !|   |   |  |
     *                                   Refl_Type(i)(1:l),!|   |   |  |
     *                                   i-1, ca4(i)       !|   |   |  |
            endif !-----------------------------------------+   |   |  |
            write (luntbl,114,err=301) i-1, Gamma_Br(i),       !|   |  |
     *                                 i-1, Beta(i),           !|   |  |
     *                                 i-1, Extinction(i)      !|   |  |
          endif  !----------------------------------------------+   |  |
        enddo  !====================================================+  |
                                                                      !^
  113   format(';',i2,'. Reflection (',3i4,') - [',a,']',:/
     *         ';',26x,'Alpha(',i2,')/x0 = ',g9.3)                    !|
  114   format(';',31x,'Gam(',i2,') = ',g9.3,' Beta=Gam(0)/Gam(',
     *                                               i2,') = ',g9.3:/
     *         ';',31x,'ExtinctionLength(',i2,') = ',g9.3,' um')      !|
                                                                      !|
c+===============================================+                     |
c|Process the linking reflection (h1-h2). For it,|                     |
c|k2= incident wace, k1= diffracted wave, because|                     |
c|k2=k0+h2, k1=k0+h1=k0+h2+(h1-h2)=k2+(h1-h2).   |                     |
c+===============================================+                     |
        if (isort.eq.0) Then   !------------------------------------+  |
c Calculate the cosine of the wavevector with surface normal (=Gam1)|  |
          i = nReflexes + 1                                        !|  |
          Gamma_Br(i) = Gamma_Br(2)                                !|  |
                                                                   !|  |
c The asymmetry factor (Gam2/Gam1):                                 |  |
          if (abs(Gamma_Br(3)) .gt. 1.D-5) Then  !----+             |  |
            Beta(i) = Gamma_Br(3)/Gamma_Br(2)        !|             |  |
          else   !------------------------------------+             |  ^
            x8 = Abs(Gamma_Br(3))-Abs(Gamma_Br(2))   !|             |  |
            if     (x8.lt.0.0D0) Then  !---------+    |             |  |
              Beta(i) = 0.D0                    !|    |             |  |
            elseif (abs(x8).lt.1.0D-10) Then !---+    |             |  |
              Beta(i) = 1.D0                    !|    |             |  |
     *                * sign(1.D0,Gamma_Br(2))  !|    |             |  |
     *                * sign(1.D0,Gamma_Br(3))  !|    |             |  |
            elseif (x8.gt.0.0D0) Then  !---------+    |             |  |
              Beta(i) = (1.D+32)                !|    |             |  |
     *                * sign(1.D0,Gamma_Br(2))  !|    |             |  |
     *                * sign(1.D0,Gamma_Br(3))  !|    |             |  |
            endif  !-----------------------------+    |             |  |
          endif   !-----------------------------------+             |  |
                                                                   !|  |
c Calculate extinction length in microns:                           |  |
          xhr(1,1) = 0.                                            !|  |
          xhi(1,1) = 0.                                            !|  |
          xhq(1,1) = 0.                                            !|  |
          xfr(1,1) = 0.                                            !|  |
          xfi(1,1) = 0.                                            !|  |
          xfq(1,1) = 0.                                            !|  |
          Bragg_Angles(1) = 0.                                     !|  |
          Call X0h1 (LatticeParameters, Wave,                      !|  |
     *           Code, rho, 1,                                     !|  |
     *           ind(1,i), ind(2,i), ind(ii,i),                    !|  |
     *           iSing, Mode,                                      !|  |
     *           jpr, Bragg_Angles, x0r, x0i,                      !|  |
     *           xhr, xhi, xhq,                                    !|  |
     *           xfr, xfi, xfq, ifail)                             !|  |
          Extinction(i) = (1.E-4)*WaveLength                       !|  |
     *     * sqrt(Gamma_Br(Incident_Place)*abs(Gamma_Br(i)))       !|  |
     /     / abs(xhr(1,1))                                         !|  |
          if (Gamma_Br(i).lt.0.) Then !---------------+             |  |
c For the Bragg case definition differs: divide by pi |             |  |
               Extinction(i) = Extinction(i) / pi    !|             |  |
            endif !-----------------------------------+             |  |
c This is the minimum extinction in the                             |  |
c case of zero-angle incidence or exit:                             |  |
          if (Extinction(i).lt.Extinction0)                        !|  |
     *                             Extinction(i) = Extinction0     !|  |
          if (Ext_min.gt.Extinction(i)) Ext_min = Extinction(i)    !|  !
                                                                   !|  |
          xhr(1,1) = 0.         !restore "0" before calling x0h2    |  |
          xhi(1,1) = 0.         !restore "0" before calling x0h2    |  |
          xhq(1,1) = 0.         !restore "0" before calling x0h2    |  |
          xfr(1,1) = 0.         !restore "0" before calling x0h2    |  |
          xfi(1,1) = 0.         !restore "0" before calling x0h2    |  |
          xfq(1,1) = 0.         !restore "0" before calling x0h2    |  |
          Bragg_Angles(1) = 0.                                     !|  |
                                                                   !|  |
c Determine reflection types:                                      !|  |
          if    (Gamma_Br(i) .gt. 10.*Total_Reflex_Angle) Then !-+  |  |
c Laue:                                                          |  |  |
            Refl_Type (i) = 'Laue'                              !|  |  |
                                                                !|  |  |
          elseif(Gamma_Br(i) .lt.-10.*Total_Reflex_Angle) Then !-|  |  |
c Bragg:                                                         |  |  |
            Refl_Type (i) = 'Bragg'                             !|  |  |
                                                                !|  |  |
          else  !------------------------------------------------|  |  |
c Surface:                                                       |  |  |
            Refl_Type (i) = 'Surface'                           !|  |  |
                                                                !|  |  |
          endif !------------------------------------------------+  |  |
                                                                   !|  |
c Print the indices with the type (Bragg/Laue)                      |  |
          l = Max0(Len_Trim(Refl_Type(i)),1)                       !|  |
          write (luntbl,112,err=301) (Ind(j,i),j=1,ii),            !|  |
     *                               Refl_Type(i)(1:l),            !|  |
     *                               i-1,Gamma_Br(i), Beta(i),     !|  |
     *                               i-1, Extinction(i)            !|  |
c         Call  Pause ()                                           !|  |
        endif   !---------------------------------------------------+  ^
  112     format('; Coupling: (',3i4,') - [',a,']'/
     *    ';',31x,'Gam(',i2,') = ',g9.3,' Beta=Gam(2)/Gam( 1) =',g9.3:/
     *    ';',31x,'ExtinctionLength(',i2,') = ',g9.3,' um')           !^
                                                                      !|
        if (iSort.eq.0) Then  !---------------------------------+      |
c Sort Gamma in decreasing order (from Laue to Bragg):          |      |
          isort = 1                                            !|      |
c         write (*,     775)                                   !|      |
          write (luntbl,775,err=301)                           !v      ^
  775     format(';'/
     *    '; Sorting reflections in the Laue-case to',1x,
     *                             'the Bragg-case order...')  !^      ^
                                                               !|      |
          Call Reflex_Sorting (                                !|      |
     *                     Gamma_Br,                           !|      |
     *                     nReflexes ,  Ind,                   !|      |
     *                     h_vectors, h_module,                !|      |
     *                     New_Ord,  New_Pos    )              !|      |
                                                               !|      |
          Incident_Place = New_pos(1)                          !|      |
c Move the incident vector to new place:                        |      ^
          do    i=1,ii  !==================================+    |      |
            k_vectors(i,Incident_Place) = k_vectors(i,1)  !|    |      |
          enddo  !=========================================+    |      |
                                                               !|      |
c Re-calculate the diffracted waves after sorting:              |      |
          goto 505   !---->---------------------->--------------+------+
        endif    !----------------------------------------------+


        if (Ext_min.gt.0.) Then !-------+
          s8 = Thickness / Ext_min     !|
        else !--------------------------+
          s8 = 100.                    !|
        endif !-------------------------+
        s8_ = exp(s8)
        write (luntbl,520,err=301) Thickness, Ext_min, s8, s8_
        if (s8_.gt.(1.E+15)) Then !---------------------------+
          write (luntbl,521,err=301) Thickness, s8, Ext_min  !|
        endif !-----------------------------------------------+
  520   format (';'/
     *  ';                             Crystal Thickness =',g9.3,' um'/
     *  '; Min_Extinction_Length among Bragg reflections =',g9.3,' um'/
     *  ';         Thickness/Min_Extinction_Length ratio =',g8.2/
     *  ';          exp(Thickness/Min_Extinction_Length) =',g8.2)
  521   format (';'/
     *  '; WARNING: crystal thickness =',g9.3,' um is more than ',f6.0/
     *  '; times the extinction length=',g9.3,' um for one or more',1x,
     *                                            'Bragg reflections.'/
     *  '; This may potentially lead to the algorithm failure due',1x,
     *  'to a precision loss because of large exponents.'/';')

        if (isort.eq.1) Then  !-------------------------+
c Print sorting results:                                |
          write (luntbl,517,err=301)                   !|
          write (luntbl,519,err=301)                   !v
  517     format(';'/
     *    '; Sorting Results:'/
     *    '; ',16('=')/
     *    '; Nr',8x,'(indices)',10x,'Gamma',6x,
     *                           'Type   New-Pos.')
  519     format('; ',55('-'))
  518     format('; ',i2,4x,'(',3i5,')',
     *                            4x,g9.3,4x,a,i4)     !^
          do    i=1,nReflexes  !==================+     |
            j = New_Pos(i)                       !|     |
            write (luntbl,518,err=301)           !|     |
     *                   i-1,                    !|     |
     *                  (ind(l,j),l=1,ii),       !|     |
     *                   Gamma_Br(j),            !|     |
     *                   Refl_Type(j),           !|     |
     *                   j                       !|     |
          enddo  !================================+     |
          write (luntbl,519,err=301)                   !|
        endif  !----------------------------------------+

        Task_size = 2 * (nReflexes + Surf_Ref_Num)

c-------------------------------------------------------
        do      i=1,nReflexes  !==============================+
c Angles of wavevectors with the internal normal:            !|
c         write (0,*) ' Read_Input_File:  ',                 !|
c    *                          'k_unit(',i,') = (',         !|
c    *                          (k_unit(j,i),j=1,ii),')'     !|
c         write (0,*) ' Read_Input_File:  ',                 !|
c    *                          'Surface_normal = (',        !|
c    *                           Surface_Normal,')'          !|
          k_Angles(i) = Dangle2 (k_unit(1,i),                !|
     *                           Surface_Normal, ii)         !|
          k_Angles(i) = 90.-k_Angles(i)                      !|
                                                             !|
c Angles for Bragg normals with the internal normal:          |
          if (i.ne.Incident_Place)  Then  !--------------+    |
c           write (0,*) ' Read_Input_File:  ',          !|    |
c    *                      'h_vectors(',i,') = (',     !|    |
c    *                      (h_vectors(j,i),j=1,ii),')' !|    |
c            write (0,*) ' Read_Input_File:  ',         !|    |
c    *                          'Surface_normal = (',   !|    |
c    *                           Surface_Normal,')'     !|    |
            if (h_module(i) .lt. 1E-32) goto 104 !-------+-+  |
            h_Angles(i) = Dangle2 (h_vectors(1,i),      !| v  |
     *                             Surface_Normal,      !|    |
     *                             ii)                  !|    |
                                                        !|    |
            if (h_Angles(i).gt.90.)                     !|    |
     *                h_Angles(i) = 180.-h_Angles(i)    !|    |
          endif  !---------------------------------------+    |
        enddo  !==============================================+

c Calculate the susceptibilities Xh of all used reflections:
        ifail = 0
        Call    x0h2    (LatticeParameters, Wave,
     *                   Code,rho, nReflexes, Ind,
     *                   Incident_Place, iSing,
     *                   Mode, jpr, Bragg_Angles,
     *                   x0r, x0i, nDimen,
     *                   xhr, xhi, xhq, xfr, xfi, xfq, ifail)
        if (ifail.ne.0) goto 130

c +-----------------------------------------+
c | ATTENTION:                              |
c | This program uses the presentation of Xh|
c | in the notation by Afanasiev-Baryshevsky|
c | ( Bloch wave = exp(ikr) )               |
c |     x0r=-abs(x0r)                       |
c |     x0i=+abs(x0i)                       |
c +-----------------------------------------+

        do      i=1,nReflexes  !========================+
c Unit vectors along SIGMA-polarization direction:     !|
          if ((90.-k_Angles(i)) .gt. 5.                !|k_Angles are the angles of
     *                  .AND.                          !|wavevectors with internal normal
     *        (90.-k_Angles(i)) .lt. 175.)   Then !--+  |
c Normally we use: sg=[k*surf_norm]                  |  |
            Call DVecVec (k_unit(1,i),              !|  |
     *                   Surface_Normal,            !|  |
     *                   Sg_unit(1,i),ii)           !|  |
          else   !-----------------------------------+  |
c When k || surf_norm, we cannot use                 |  |
c sg=[k*surf_norm] and then we choose:               |  |
c sg=[k*[h1*h2]]                                     |  |
            Call DVecVec (k_unit(1,i),              !|  |
     *                    h1xh2,                    !|  |
     *                    Sg_unit(1,i),ii)          !|  |
          endif      !-------------------------------+  |
                                                       !|
          Call DUnitVec2 (Sg_unit(1,i),                !|
     *                    Sg_unit(1,i),ii)             !|
                                                       !|
c Unit vectors along PI-polarization dir-on: pi=[sg*k]: |
          Call  DVecVec (Sg_unit(1,i),                 !|
     *                   k_unit(1,i),                  !|
     *                   Pi_unit(1,i), ii)             !|
        enddo  !========================================+

c Define the polarization plane:
        if     (ipol.eq.1) Then !-----------------------+ipol=1: given by angle to pi0
                                                       !|
c Orientation was prompted:                            !|
          csg = Dsin(polang*Dgra)                      !|
          cpi = Dcos(polang*Dgra)                      !|
                                                       !|
        elseif (ipol.eq.2) Then !-----------------------+ipol=2: along [k0*k1]
                                                       !|
c Emulation of V.Kohn paper in Sov.Phys.Cristallography |
c                 (E0.ort.k0 & E0.ort.k1)              !|
          Call DvecVec  (k_vectors(1,Incident_Place),  !|
     *                   k_vectors(1,New_pos(2)),      !|
     *                   Incident_Polarization,        !|
     *                   ii)                           !|
          Call DUnitVec2(Incident_Polarization,        !|
     *                   Incident_Polarization,        !|
     *                   ii)                           !|
          cpi = DVecSca2 (Incident_Polarization,       !|
     *                    Pi_unit(1,Incident_Place),   !|
     *                    ii)                          !|
          csg = DVecSca2 (Incident_Polarization,       !|
     *                    Sg_unit(1,Incident_Place),   !|
     *                    ii)                          !|
          polang = SNGL (Datan2(csg,cpi) / Dgra)       !|
                                                       !|
        endif   !<--------------------------------------+

c+======================================================+
c|Find the coefficients in formula:                     |
c|                                                      |
c|alf(m)= (2/k_module) * [(a1*h(m))*dq1 + (a2*h(m))*dq2]|
c|                                                      |
c|where a1 and a2 are unit vectors perpendicular to k0; |
c|(Pinsker book, page 346).                             |
c+======================================================+
        goto (660,661,662)      iTet+1    !-------------+
                                                       !|
c 0- normal choice: a1=pi0=pp(1),  a2=sg0=sg(1).        |
  660   continue  !<------------------------------------+
        do  i=1,ii !=========================+          |
          a1(i) = Pi_unit(i,Incident_Place) !|          |
          a2(i) = Sg_unit(i,Incident_Place) !|          |
        enddo   !============================+          |
        goto 666  !---------------------------------+   |
c                                                   v   |
c 1- (a2.ort.norm & a2.ort.k0), (a1.ort.a2 & a1.ort.k0).|
  661   continue  !<------------------------------------+
        Call  DVecVec   (k_unit(1,Incident_Place),     !|
     *                   Surface_Normal,               !|
     *                   a2, ii)                       !|
        Call  DUnitVec2 (a2,a2, ii)                    !|
        Call  DVecVec   (k_unit(1,Incident_Place),     !|
     *                   a2,                           !|
     *                   a1, ii)                       !|
        goto 666  !----------------------------------+  |
c                                                    v  |
c 2- (a2.ort.k1 & a2.ort.k0), (a1.ort.a2 & a1.ort.k0)   |
c    i.e. a2 is along E0 (see above); emulation of the  |
c    V.Kohn pager (Sov.Phys.Crystallography, 1975);     |
  662   continue  !<------------------------------------+
        Call  DVecVec   (k_unit(1,Incident_Place),
     *                   k_unit(1,New_pos(2)),
     *                   a2, ii)                    !|
                                                    !v
c ccc   cmTet(3) = 'Theta2 .ort. ( ??? ??? ???)'    !|
        write (cmTet(3)(15:26),'(3i4)')             !|
     *                  (Ind(j,New_pos(2)),j=1,ii)  !|
                                                    !|
        Call  DUnitVec2 (a2, a2, ii)                !|
        Call  DVecVec   (k_unit(1,Incident_Place),  !|
     *                   a2,                        !|
     *                   a1, ii)                    !v
c                                                   !|
c The divider "/xabs" is due to our normalization of |
c the matrix (before it was "(2.*Wavelength)"):      |
c                                                    |
  666   continue  !<---------------------------------+
        Coef = (2./k_module)*(UnitCoef(2)/xabs)

        do      i=1,nReflexes    !=======================+
c Coefficient at Tet1:                                   |
          ca1(i) = Coef*DVecSca2 (a1,                   !|
     *                            h_vectors(1,i),       !|
     *                            ii )                  !|
c Coefficient at Tet2:                                   |
          ca2(i) = Coef*DVecSca2 (a2,                   !|
     *                            h_vectors(1,i),       !|
     *                            ii )                  !|
c Coefficient at (Tet1**2+Tet2**2):                      |
c (that is we use extra UnitCoef factor                  |
c to convert the angles into radians)                    |
          ca3(i) = Coef * UnitCoef(2) *                 !|
     *            (1.-DVecSca2(k_unit(1,Incident_Place),!|
     *                         k_unit(1,i),             !|
     *                         ii) )                    !|
        enddo    !=======================================+

c Coefficient in the equation for changing Gamma0:
        Coef = UnitCoef(2)
        cSurf1 = Coef*DVecSca2(a1,Surface_Normal,ii)
        cSurf2 = Coef*DVecSca2(a2,Surface_Normal,ii)

c ##################################################### 2
c Absorption coefficient in 1/cm:
        mu0 = SNGL ((2.e+08)*pi*abs(x0i) / Wavelength)
c Beam attenuation by the crystal plate:
        if (Gamma_Br(Incident_Place) .gt. Total_Reflex_Angle) Then  !---+
          Demp8 = (1.D-04)*Thickness * mu0 / Gamma_Br(Incident_Place)  !|
        else  !---------------------------------------------------------|
          Demp8 = (1.D-04)*Thickness * mu0 / Total_Reflex_Angle        !|
        endif  !--------------------------------------------------------+
        eqormore = '='
        attn     = ' '
        if (Demp8.gt.Exp8) Then  !--------+
          Demp8 = Exp8                   !|
          eqormore = '>'                 !|
          attn     = '<- thick crystal!' !|
        endif  !--------------------------+
        Demp8 = Exp(Demp8)
c-------------------------------------------------------
        if (iSort.eq.0) Then  !--------------------------------------+
          write (luntbl,401,err=301) Code, Thickness,               !|
     *                               Wavelength,                    !|
     *                               Total_Reflex_Angle/UnitData(2) !|
          write (luntbl,431,err=301) nBaseNormal, DisorNorm,        !|
     *                               nReperDirection                !|
        endif !------------------------------------------------------+
  401   format  (';'/
     *  '; ',63('-')/
     *  ';',40x,' Crystal:  ',a/
     *  '; Crystal thickness = ',g12.5, '   micron'/
     *  '; X-ray  wavelength = ',g18.10,' Angstrem':/
     *  '; Total Refl. Angle = ',g12.5, '  arc min.')
  431   format  (
     *  '; External surface normal - base indices: (',3i4,')'/
     *  '; Misorientation by ',f8.2,' degr. along: (',3i4,')'/';')
  435   format  ('; Choice of k0 if two options found: ',a/';')

        do      i=1,nReflexes  !==============================+
          j = New_pos(i)                                     !|
          l = Max0(Len_Trim(Refl_Type(j)),1)                 !|
          if (j.ne.Incident_Place) Then !-----------------+   |
            write (luntbl,113,err=301) i-1,              !|   |
     *                                (ind(m,j),m=1,ii), !|   |
     *                                 Refl_Type(j)(1:l),!|   |
     *                                 i-1, ca4(j)       !|   |
          else !------------------------------------------+   |
            write (luntbl,113,err=301) i-1,              !|   |
     *                                (ind(m,j),m=1,ii), !|   |
     *                                 Refl_Type(j)(1:l) !|   |
          endif !-----------------------------------------+   |
          write (luntbl,114,err=301) i-1,Gamma_Br(j),        !|
     *                               i-1,Beta(j),            !|
     *                               i-1, Extinction(j)      !|
          if (j.ne.Incident_Place)  Then  !---------------+  !|
            write (luntbl,115,err=301) BraggWave(j)      !|   |
  115       format(';',20x,'"Bragg" wavelength = ',a)    !|   |
            write (luntbl,403,err=301) Bragg_Angles(j),  !|   |
     *                                 h_Angles(j)       !|   |
          endif  !----------------------------------------+   |
          write (luntbl,432,err=301) k_Angles(j)             !|
          write (luntbl,116,err=301) (k_unit(m,i),m=1,ii),   !|
     *                               (Sg_unit(m,i),m=1,ii),  !|
     *                               (Pi_unit(m,i),m=1,ii)   !|
        enddo   !=============================================+

        write (luntbl,112,err=301) (
     *             Ind(j,nReflexes+1),j=1,ii),
     *             Refl_Type(nReflexes+1),
     *             nReflexes, Gamma_Br(nReflexes+1), Beta(nReflexes+1),
     *             nReflexes, Extinction(nReflexes+1)

  403   format  (';',27x,'Bragg angle = ',f7.3,' Degr.'/
     *  ';',9x,'Angle of plane to the surface = ',f7.3,
     *  ' Degr.')
  432   format  (';',6x,'Angle of the wave to the surface = ',
     *  f7.3,' Degr.')
  116   format(';',11x,'Unit vector along  k-vector : (',3g9.3,')'/
     *         ';',11x,'Unit vector along sigma-pol : (',3g9.3,')'/
     *         ';',11x,'Unit vector along    pi-pol : (',3g9.3,')')

c Calculate angles between all planes:
        do      m=1,nReflexes-1  !=====================================+
          do    i=m+1,nReflexes  !===================================+ |
            if (m.ne.Incident_Place                                 !| |
     *                .AND.                                         !| |
     *          i.ne.Incident_Place)  Then  !--------------------+   | |
c             write (0,*) ' Read_Input_File:  ',                !|   | |
c    *                      'h_vectors(',m,') = (',             !|   | |
c    *                      (h_vectors(j,m),j=1,ii),')'         !|   | |
c             write (0,*) ' Read_Input_File:  ',                !|   | |
c    *                      'h_vectors(',i,') = (',             !|   | |
c    *                      (h_vectors(j,i),j=1,ii),')'         !|   | |
              if (h_module(i) .lt. 1E-32) goto 104 !-------------+-+ | |
              Interplane_Angles = 180. -                        !| v | |
     -               Dangle2(h_vectors(1,m),h_vectors(1,i),ii)  !|   | |
              write (luntbl,448,err=301) New_ord(m)-1,          !|   | |
     *                                   New_ord(i)-1,          !|   | |
     *                                   Interplane_Angles      !|   | |
            endif  !---------------------------------------------+   | |
          enddo  !===================================================+ |
        enddo  !=======================================================+
  448   format (';'/
     *  '; Angle between planes (',i2,' -',i2,') = ',g14.6,' Degr.')

        write   (luntbl,55,err=301)     h1xh2_indices
  55    format  ('; Zone vector [ h1 x h2 ] indices: (',3i5,')')
        write   (luntbl,58,err=301)     k0_indices

        write   (luntbl,421,err=301)    x0r,x0i
        nforbidden = 0
        nforbid_pi = 0
        do n=1,nReflexes !=======================================+
        do m=1,nReflexes !====================================+  |
          i1 = Ind(1,m)-Ind(1,n)                             !|  |
          i2 = Ind(2,m)-Ind(2,n)                             !|  |
          i3 = Ind(3,m)-Ind(3,n)                             !|  |
          s8 = 0.                                            !|  |
          x8 = DBLE(xhr(m,n))                                !|  |
          s8 = s8 + x8*x8                                    !|  |
          x8 = DBLE(xhi(m,n))                                !|  |
          s8 = s8 + x8*x8                                    !|  |
          x8 = DBLE(xhq(m,n))                                !|  |
          s8 = s8 + x8*x8                                    !|  |
          s8 = SQRT(s8)  / DBLE(xabs)                        !|  |
          iii(1)=i1                                          !|  |
          iii(2)=i2                                          !|  |
          iii(3)=i3                                          !|  |
          b8 = Bragan8 (Wave,iii,3,2,d8,0,s4,s4,s4,s4,s4,s4) !|  |
          s8_ = s8 * abs(cos(2.*b8*Dgra))                    !|  |
          if (s8. gt. (1.D-3)) Then !--------------------+    |  |
c Reflection is okay (stronger than 0.1%):               |    |  |
            if (s8_ .gt.(1.D-3)) Then !----+             |    |  |
              weak = ' '                  !|             |    |  |
            else !-------------------------+             |    |  |
              weak = '*** pi-weak ***'    !|             |    |  |
              nforbid_pi = nforbid_pi + 1 !|             |    |  |
            endif !------------------------+             |    |  |
            write (luntbl,422,err=301)                  !|    |  |
     *          m-1, n-1, xhr(m,n), xhi(m,n), xhq(m,n), !|    |  |
     *          i1, i2, i3, 100.*s8, 100.*s8_, weak     !|    |  |
           else !----------------------------------------+    |  |
c Reflection is forbidden (weaker than 0.1%):            |    |  |
             write (luntbl,423,err=301)                 !|    |  |
     *          m-1, n-1, xhr(m,n), xhi(m,n), xhq(m,n), !|    |  |
     *          i1, i2, i3                              !|    |  |
            xhr(m,n) = 0.                               !|    |  |
            xhi(m,n) = 0.                               !|    |  |
            xhq(m,n) = 0.                               !|    |  |
            nforbidden = nforbidden + 1                 !|    |  |
c           nforbid_pi = nforbid_pi + 1                 !|    |  |
          endif !----------------------------------------+    |  |
        enddo !===============================================+  |
        enddo !==================================================+
        write   (luntbl,434,err=301)    mu0, eqormore, Demp8, attn
c Error: more than one
c forbidden reflex:
        if (nforbidden .gt. 2) goto 103
        if (nforbid_pi .gt. 0) write (luntbl,424,err=301) nforbid_pi

  421   format  ('; Polarizabilities:'/
     *           ';     X0    = ',g12.5,'  + i *',g12.5)
  422   format  ('; Xh(',i2,',',i2,') = ',g12.5,2('  + i *',g12.5),
     *                        '  hkl = (',3i4,')  strength = ',f5.1,'%',
     *                        2x,'strength(pi) = ',f5.1,'%', 2x,a)
  423   format  ('; Xh(',i2,',',i2,') = ',g12.5,2('  + i *',g12.5),
     *                        '  hkl = (',3i4,') *** forbidden ***')
  424   format  (';'/
     *           '; WARNING: there are ',i2,' reflections weakened for',
     *        1x,'pi-polarization due to Bragg angle close to 45 degr.'/
     *           '; This may lead to the algoritm failure.'/';')
  433   format  (
     *  '; +--------------',       4('+',11('-')),                  '+'/
     *  '; |     Angle    |  Points   |   Step    |   Start   |',
     *                                                   '  Finish   |'/
     *  '; +--------------',       4('+',11('-')),                  '+'/
     *  '; |Theta1 (',a5,')|',  i10, ' |',                3(g10.3,' |')/
     *  '; |Theta2 (',a5,')|',  i10, ' |',                3(g10.3,' |')/
     *  '; +--------------',       4('+',11('-')),                  '+')
  434   format  (';'/
     *  '; Linear  absorption  factor  mu = ',g12.5,' 1/cm'/
     *  '; Attenuation: exp(mu*t/Gamma_Br(0)) ',a,1x, g14.5, 5x, a/
     *  ';')

        if (ipol.eq.0)  Then   !------------------------------------+
c         write (   *,  219)                                       !|
          write (luntbl,219,err=301)                               !|
        else   !----------------------------------------------------+
c         write (   *, 319)          polang                        !|
          write (luntbl,319,err=301) polang                        !|
          do    i=1,nReflexes  !=================================+  |
            j = New_pos(i)                                      !|  |
            if (j.ne.Incident_Place)    Then  !---------------+  |  |
                                                             !|  |  |
c A vector in vertical direction if k0 and k(i)               |  |  |
c lie in horizontal plane:                                    |  |  |
c             write (0,*) ' Read_Input_File:  ',             !|  |  |
c    *        'k_unit(',Incident_Place,') = (',              !|  |  |
c    *        (k_unit(m,Incident_Place),m=1,ii),')'          !|  |  |
c             write (0,*) ' Read_Input_File:  ',             !|  |  |
c    *        'k_unit(',j,') = (',                           !|  |  |
c    *        (k_unit(m,j),m=1,ii),')'                       !|  |  |
              Call DVecVec(k_unit(1,Incident_Place),         !|  |  |
     *                     k_unit(1,j),                      !|  |  |
     *                     Vertical, ii)                     !|  |  |
              if (Dvecmod2(Vertical, ii) > 1.E-20) Then !-+   |  |  |
c Angles of Sg0 and Pi0 with the vertical:                !|  |  |  |
c               write (0,*) ' Read_Input_File:  ',        !|  |  |  |
c    *          'Pi_unit(',Incident_Place,') = (',        !|  |  |  |
c    *          (Pi_unit(m,Incident_Place),m=1,ii),')'    !|  |  |  |
c               write (0,*) ' Read_Input_File:  ',        !|  |  |  |
c    *                     'Vertical = (',Vertical,')'    !|  |  |  |
                TetaPi = SNGL( Dangle2(                   !|  |  |  |
     *                    Pi_unit(1,Incident_Place),      !|  |  |  |
     *                    Vertical, ii ) )                !|  |  |  |
                if (TetaPi.gt.90.) Then !-----+            |  |  |  |
                  TetaPi = TetaPi - 90.      !|            |  |  |  |
                endif !-----------------------+            |  |  |  |
                TetaSg = TetaPi + 90.                     !|  |  |  |
                TetaIn = TetaPi + polang                  !|  |  |  |
                write  (luntbl,320,err=301)               !|  |  |  |
     *                    (Ind(l,j),l=1,ii),              !|  |  |  |
     *                    TetaSg,                         !|  |  |  |
     *                    TetaPi,                         !|  |  |  |
     *                    TetaIn                          !|  |  |  |
              else   !-------------------------------------+  |  |  |
                write  (luntbl,321,err=301)               !|  |  |  |
     *                    (Ind(l,j),l=1,ii),              !|  |  |  |
     *                 (k_unit(l,j),l=1,ii),              !|  |  |  |
     *                 (k_unit(l,Incident_Place),l=1,ii)  !|  |  |  |
              endif  !-------------------------------------+  |  |  |
            endif  !------------------------------------------+  |  |
          enddo  !===============================================+  |
        endif  !----------------------------------------------------+

  219   format('; Averaged over incident polarization')
  319   format('; Angle between Polariz.plane & Pi(0): ',f10.3,' Degr.')
  320   format ('; -- If diffracted vector for (',3i4,
     *      ') propagates horizontally, then:'/
     *      '; Sg(0) makes with Vertical ',f5.1,' Degr.'/
     *      '; Pi(0) makes with Vertical ',f5.1,' Degr.'/
     *      '; Incident polarization -"- ',f5.1,' Degr.')
  321   format ('; !! ** diffracted vector k(',3i4,') = (',
     *      3f6.3,') is parallel to k0 = ('
     *      3f6.3,').'/
     *      '; Cannot choose polarization vector for this reflection.')

c       write   (  *,   419)            cmTet(iTet+1)
        write   (luntbl,419,err=301)    cmTet(iTet+1)
  419   format  ('; Direction of scan axes: ',a)

        x0  = cmplx (x0r,x0i)/xabs
c Complex matrices of suspeptibilities (nReflexes*nReflexes):
        do      i=1,nReflexes  !=========================+
          do    j=1,nReflexes    !=====================+ |
            xhq(i,j) = xhq(i,j) / xabs                !| |
            xh(i,j)  = (xhr(i,j)                      !| |
     *               * exp(+Imaginary1*pi*xfr(i,j))   !| |
     +               + Imaginary1                     !| |
     *               * xhi(i,j)                       !| |
     *               * exp(+Imaginary1*pi*xfi(i,j)))  !| |
     /               / xabs                           !| |
          enddo   !====================================+ |
        enddo   !========================================+

c Filter on large Bragg deviations of incident vectors from the given Bragg conditions:
c ((k0+h)^2-k0^2)/k0^2 in the inits of |x0] should not exceed devmax
        do      i=1,nReflexes  !==================+
          if (Dabs(ca4(i)).gt.devmax)  goto 102  !|
        enddo  !==================================+

        ifail = 0
  28    continue
        if (iirezv .eq. 1)  Then  !--------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c#######################################################
c                       E R R O R S :
c#######################################################
  100   continue
        write   (txt,118)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line
  118   format  (a,': Incorrect input parameters in file:'//a//
     *  'at line:',i4,' -- while reading:'/)
        goto 120

c-------
  101   continue
        write   (txt,119)  ProgramVer(1:ipv), InpFile(1:linp),
     *                     line
  119   format  (a,': Parameter(s) out of range in file:'//a//
     *  'at line:',i4,' -- while reading:'/)
        goto 120

c-------
  104   continue
        write   (txt,121)  ProgramVer(1:ipv), InpFile(1:linp), i
  121   format  (a,': Incorrect data in input file:'//a//
     *  ' -- reflection-',i2,' has zero length.'/)
        ifail = 6
        goto 130

c-------
  120   txt(8) = ' '
        txt(9) = 'Check your input file!'
        ifail  = 9
        goto 130

c-------
  99    continue
        write   (txt,117)  ProgramVer(1:ipv), InpFile(1:linp)
  117   format  (a,': cannot open input file:'//a//
     *  'File not found in current directory.')
        ifail = 5
        goto 130

c-------
  102   continue
        write   (txt,807) ProgramVer(1:ipv), (ind(j,i),j=1,ii),
     *                    Dabs(ca4(i)), devmax
  807   format  (a,': Too large initial Bragg deviation for',
     *                                      ' reflection (',3i4,'):'//
     *  '|alfa|/|x0|=',e8.3,'  >  ',f4.0,',     where',
     *                                 ' alfa=((k0+h)^2-k0^2)/k0^2.'//
     *  'Brl cannot handle reflections which are very far from their',
     *                                            ' Bragg conditions'/
     *  'because their effect is vanishingly small, thus leading to',
     *                                         ' numerical problems.')
        ifail = 6
        goto 130

c-------
  103   continue
        write   (txt,808) ProgramVer(1:ipv)
  808   format  (a,': multiple forbidden reflections found.'/
     *           'This means that either more than one forbidden'/
     *           'reflection is specified or a single forbidden'/
     *           'reflection cannot be excited because its coupling'/
     *           'reflection is also forbidden.'/
     *           /
     *           'See the TBL file for more details.'/
     *           /
     *           'Please remove multiple or non-excitable forbidden'/
     *           'reflections from the reflections list because'/
     *           'they are resulted in singular scattering matrix.')
        ifail = 11
        goto 130

c-------
  109   continue
        write   (txt,209)  ProgramVer(1:ipv), InpFile(1:linp), i
  209   format  (a,': Inconsistent parameters in file:'//a//
     *          'Multiple scan points with zero scan range'/
     *          'while reading Scan Angle-',i1,'.')
        ifail = 6
        goto 130


c-------
  301   continue
        write   (txt,809) ProgramVer(1:ipv), OutFile(1:lout)//'.tbl'
  809   format  (a,': error writing ',a)
        ifail = 1
        goto 130

c-------
  130   continue
        if (modebat.eq.1) Then  !----------------+
          do    i=1,Iabs(ifail)  !=============+ |
            l = Max0(Len_Trim(txt(i)),1)      !| |
            write (*,'(1x,a)') txt(i)(1:l)    !| |
          enddo  !=============================+ |
        endif  !---------------------------------+
        Call    Message (txt,Iabs(ifail),2)
        goto 28

        end

c ====================================================

        Subroutine Reflex_Sorting (gam,ndim,ind,h,hl,no,np)

        Integer ndim, ind(3,ndim), m(3),
     *          no(ndim), np(ndim),
     *          m4, i, j, k
        Real*8  gam(ndim), g, h(3,ndim), r(3),
     *          hl(ndim), r4
c+============================================================+
c|gam(ndim) - initial real-type array of ndim dimension, which|
c|            is re-sorted in this sub in decreasing order    |
c|ind(3,ndim) - initial integer-type array transformed in the |
c|              same order as gam during sorting              |
c|no:    new order of gam elements (e.g. no(1)=4 - former 4th |
c|       element of gam is now on the 1st place)              |
c|np:    new positions of game elements (e.g. np(1)=5 - former|
c|       1th element of gam is now on 5th place)              |
c|On entry: no=1,2,3,4,...                                    |
c+============================================================+
c The main double loop of sorting over decreasing gam
c (the code is adapted from the 'Numerical Recipes' book):
        do j=2,ndim  !===========================+
          g = gam(j)                            !|
          do      k=1,3  !=====+                 |
            m(k) = ind(k,j)   !|                 |
            r(k) = h(k,j)     !|                 |
          enddo  !=============+                 |
          m4 = no(j)                            !|
          r4 = hl(j)                            !|
                                                !|
          do i=j-1,1,-1 !=================+      |
            if (gam(i).ge.g)    goto 3  !>+--+   |
            gam(i+1) = gam(i)            !|  |   |
            do  k=1,3  !==============+   |  |   |
              ind(k,i+1) = ind(k,i)  !|   |  |   |
              h(k,i+1)   = h(k,i)    !|   |  |   |
            enddo  !==================+   |  |   |
            no(i+1) = no(i)              !|  |   |
            hl(i+1) = hl(i)              !|  |   |
          enddo  !========================+  |   |
                                            !|   |
          i=0                               !|   |
  3       continue         !<----------------+   |
          gam(i+1) = g                          !|
          do k=1,3  !============+               |
            ind(k,i+1) = m(k)   !|               |
            h(k,i+1)   = r(k)   !|               |
          enddo  !===============+               |
          no(i+1) = m4                          !|
          hl(i+1) = r4                          !|
        enddo  !=================================+

        do  i=1,ndim !=====+
          np(no(i)) = i   !|
        enddo  !===========+

        return
        end

c ====================================================

        Subroutine Sort_by_Epsilon_Im_Part
     *                           (mDim,nUse,deltr,delti,
     *                            dr,di,wr,wi)

        Integer mDim, nUse, i, j, m
        Real*8  deltr(mDim), delti(mDim),
     *          dr(mDim,mDim), di(mDim,mDim),
     *          wr(mDim), wi(mDim), Vr, Vi
c--------------------------------------------------------
c Sort (deltr,delti) in the decreasing order of the imaginary part "delti":
c--------------------------------------------------------
c The main double loop of sorting over decreasing "delti"
c (the code is adapted from the 'Numerical Recipes' book):
        do j=2,nUse !===========================+
          Vi = delti(j)                        !|
          Vr = deltr(j)                        !|
          do m=1,nUse !=======+                 |
            wr(m) = dr(m,j)  !|                 |
            wi(m) = di(m,j)  !|                 |
          enddo !=============+                 |
                                               !|
          do    i=j-1,1,-1 !===============+    |
            if (delti(i).ge.vi) goto 3 !>+ |    |
            delti(i+1) = delti(i)       !| |    |
            deltr(i+1) = deltr(i)       !| |    |
            do  m=1,nUse   !========+    v |    |
              dr(m,i+1) = dr(m,i)  !|    | |    |
              di(m,i+1) = di(m,i)  !|    | |    |
            enddo !=================+    | |    |
          enddo !========================+=+    |
                                        !|      |
          i = 0                         !|      |
  3       continue !<--------------------+      |
          delti(i+1) = Vi                      !|
          deltr(i+1) = Vr                      !|
                                               !|
          do m=1,nUse !=========+               |
            dr(m,i+1) = Wr(m)  !|               |
            di(m,i+1) = Wi(m)  !|               |
          enddo !===============+               |
        enddo !=================================+

        return
        end
