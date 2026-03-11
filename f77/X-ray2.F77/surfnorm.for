c-----------------------------------------------------------------------
c This file contains the two groups of programs to determine internal
c surface normal to crystal surface.
c
c 1a. InputSufNorm reads from INP file the base surface normal, the
c                  miscut angle, and the direction of miscut (the
c                  reference direction).
c 1b. SurfNormVect uses base surface normal and miscut angle and the
c                  direction of miscut (the reference direction) to
c                  build the actual crystal normal.
c
c 2a. SurfCoplanar  uses the Bragg plane and the asymmetry of the Bragg
c                   reflection (the angle of Bragg planes to the surface)
c                   to build a generally accidental normal from an infinite
c                   ensemble satisfying the condition. This program is
c                   for coplanar case where the real normal is not specified.
c 2b. FindDisorient is the auxiliary program to SurfCoplanar which derives
c                   the Bragg planes angle to the surface using several
c                   different inputs.
c
c These subroutines are called from gid_slXX and gids_XX.
c
c                Copyright Sergey Stepanov 1994--2001
c-----------------------------------------------------------------------

        Subroutine InputSufNorm (lun,            !input file LUN (input)
     *                           line,           !current input file line (input/output)
     *                           igie,           !flag how the geometry is build (input)
     *                           wrk,            !text string work space (input)
     *                           nBaseNormal,    !hkl of base normal (output)
     *                           nReperDirection,!hkl of miscut reper (output)
     *                           DisorNorm8,     !miscut angle in Units1 (output)
     *                           Unit2degree,    !conversion factor Units1->degrees (input)
     *                           Surface_Normal, !hkl of resulted surface normal (output)
     *                           *, *, *)        !labels to proceed on errors
c -------------------------------------------------------
c This subroutine reads the base EXTERNAL surface normal,
c the miscut direction, and the miscut angle for gid_slXX
c programs and then calculates the real INTERNAL normal
c for disoriented surface.
c "igie" is the flag selecting how the diffraction geometry
c is specified:
c  [1]=non-coplanar via surface normal and the angle of k0
c  [2]=non-coplanar via surface normal and the angle of kh
c  [3]=coplanar via surface normal and the condition of grazing incidence
c  [4]=coplanar via surface normal and the condition of grazing exit
c  [5]=symmetric Bragg case (any) via surface normal
c  [6]=coplanar via Bragg planes angle to the surface (no surface normal needed)
c  [7]=coplanar via incidence angle of k0 (no surface normal needed)
c  [8]=coplanar via exit angle of kh (no surface normal needed)
c  [9]=coplanar via beta=g0/|gh| (no surface normal needed)
c -- thus, for [igie > 5] surface normal is not required and
c therefore the consistency of input is not analyzed.
c -------------------------------------------------------
        Real*8          DisorNorm8, DisorNorm8_gra,
     *                  Surface_Normal(3)
        Real*4          Unit2degree
        Integer         lun, line, igie, ifail,
     *                  nBaseNormal(3), nReperDirection(3)
        Character       wrk*(*)

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c External surface normal orientation parameters:
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c 1. nBaseNormal - base surface normal vector indices:
        txt(7) = 'base surface normal indices'
        Call    LineLoop (lun,line,wrk,*100)
        Call    rdInt   (nBaseNormal,3,wrk,ifail)
        if (ifail.ne.0)                 goto 100
        if (igie.le.5) Then  !--------------------------------+
          txt(7) = 'base surface normal indices are all zero'!|
          if ((nBaseNormal(1).eq.0) .AND.                    !|
     *        (nBaseNormal(2).eq.0) .AND.                    !|
     *        (nBaseNormal(3).eq.0))    goto 101             !|
        else !------------------------------------------------|
          nBaseNormal(1) = 0                                 !|
          nBaseNormal(2) = 0                                 !|
          nBaseNormal(3) = 0                                 !|
        endif  !----------------------------------------------+

c 2. nReperDirection - indices of vector, pointing out
c    the direction of maximum misorientation of real
c    surface normal with respect to base normal vector:
        txt(7) = 'surface miscut direction indices'
        Call    LineLoop (lun,line,wrk,*100)
        Call    rdInt   (nReperDirection,3,wrk,ifail)
        if (ifail.ne.0)                 goto 100

c 3. DisorNorm8 - Maximum misorientation angle of real
c    surface normal:
        txt(7) = 'surface miscut angle (range=[-90.:90.]degr)'
        Call    LineLoop (lun,line,wrk,*100)
c ATTENTION: this is a double-precision READING!
        Call    rdReal8 (DisorNorm8,1,wrk,ifail)
        if (ifail.ne.0)                 goto 100
c       Unit2degree = UnitCoef(1) / AngDegree
        DisorNorm8_gra  =  DisorNorm8 * Unit2degree
        if (DisorNorm8_gra .le. -90.0D0 .OR.
     *      DisorNorm8_gra .ge. 90.0D0) goto 101
        if (igie.le.5) Then  !-------------------------------------+
          txt(7) = 'surface miscut direction indices are all zero'!|
          if ((abs(DisorNorm8).gt.1.E-32) .AND.                   !|
     *        (nReperDirection(1).eq.0)   .AND.                   !|
     *        (nReperDirection(2).eq.0)   .AND.                   !|
     *        (nReperDirection(3).eq.0)) goto 101                 !|
        else  !----------------------------------------------------|
          DisorNorm8 = 0.0                                        !|
          nReperDirection(1) = 0                                  !|
          nReperDirection(2) = 0                                  !|
          nReperDirection(3) = 0                                  !|
        endif  !---------------------------------------------------+
c-----------------------
c Calculate internal surface normal vector:
        if (igie.le.5) Then  !--------------------------+
          Call  SurfNormVect (nBaseNormal,             !|
     *                        nReperDirection,         !|
     *                        DisorNorm8_gra,          !|
     *                        Surface_Normal,ifail)    !|
          if (ifail.ne.0) goto 28                      !|
        else !------------------------------------------|
          Surface_Normal(1) = 0.0                      !|
          Surface_Normal(2) = 0.0                      !|
          Surface_Normal(3) = 0.0                      !|
        endif  !----------------------------------------+
        return
c -------------------------------------------------------
  100   continue                !read error
        return 1
  101   continue                !parameters are not in range
        return 2
  28    continue
        return 3                !other error (failure from SurfNormVect)
        End

c ======================================================================

        Subroutine SurfNormVect (nBaseNormal,     !hkl of base EXTERNAL normal (input)
     *                           nReperDirection, !hkl of miscut reper (input)
     *                           Disorient8,      !miscut angle inside the slab in the reper direction in Units1 (input)
     *                           Surface_Normal,  !hkl of resulted surface normal (output)
     *                           ifail)           !failure flag (output)
c +-----------------------------------------------------------+
c |  ATTENTION!!! This subroutine does not want to work with  |
c |               NDP FORTRAN if included into a library!     |
c +-----------------------------------------------------------+

c +====================================================+
c ||||||||||| Determination of internal unit |||||||||||
c |||||||||||    normal to crystal surface   |||||||||||
c +====================================================+
        Real*8  Base_Normal(3),
     *          Projection(3),
     *          Surface_Normal(3),
     *          Disorient8,
     *          Dpi8, Dgra8, s8

        Real*8  DAngle, DAngle2, DVecSca2
        External DAngle, DAngle2, DVecSca2

        Integer nBaseNormal(3),
     *          nReperDirection(3),
     *          ifail, lines, lpn, i

        Character       txt(20)*80
        Common  /msg/   txt

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
          progname = 'SurfNormVect'       !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim(progname)
c-------------------------------------------------------
        ifail = 0
        if (nBaseNormal(1).eq.0 .and.
     *      nBaseNormal(2).eq.0 .and.
     *      nBaseNormal(3).eq.0)        goto 100
        do      i=1,3   !=========================+
          Base_Normal(i) = Real(nBaseNormal(i))  !|
        enddo   !=================================+
        Call    DUnitVec2 (Base_Normal,Base_Normal,3)

        if (abs(Disorient8).lt.1.E-32) Then  !------------------------+
c Unit vector along internal surface normal:                          |
          Call  DVecCon (Base_Normal,-1.D0,                          !|
     *                  Surface_Normal,3)                            !|
        else  !-------------------------------------------------------|
          if (nReperDirection(1).eq.0 .and.                          !|
     *        nReperDirection(2).eq.0 .and.                          !|
     *        nReperDirection(3).eq.0)                     goto 101  !|
          s8 = Dabs(Dangle(nReperDirection,nBaseNormal,3))           !|
          if ((s8.lt.1.0D0) .or. (s8.gt.179.0D0))          goto 102  !|   external
c+=============================================================+      |  Base_Normal
c|First, calculate the EXTERNAL normal. We seek it in the form |      |     ^-------^
c|(not normalized):                                            |      |     |      /|External surface normal
c|             ->            ->               ->               |      |     |Disor/
c|      Surface_Normal = Base_Normal + c * Projection,         |      |     |angl/  |
c|                                                             |      |     |   /
c|          ->                                                 |      |     |  /    |
c|where Projection is a unit vector along the projection of    |      |     | /           Projection
c|         ->                               ->                 |      |    ================>======= Base plane
c|Disorient_Reper onto the plane normal to Base_Normal.        |      |      ~ _ Disor > 0 for inside miscut
c|                                                             |      |          ~ _
c|Therefore:                                                   |      |              ~ _
c|          ->           ->                                    |      |                  ~ _
c|   (Surface_Normal*Projection) = c                           |      |                      ~ surface
c|                                                             |      |
c|       ->           ->                                       |      |
c|   (Surface_Normal*Base_Normal) = 1                          |      |
c|                                                             |      |
c|And (see picture at the right);                              |      |
c|                                                             |      |
c|   c = tan(Disorient8)                                       |      |
c+=============================================================+      |
          Dpi8  = (4.D0)*Datan(1.D0)                                 !|
          Dgra8 = (2.D0)*Dpi8/(360.D0)                               !|
c Calculate the projection of chosen reper on the plane               |
c perpendicular to the base normal:                                   |
          Call DProject  (nBaseNormal,nReperDirection,Projection,3)  !|
          Call DUnitVec2 (Projection,Projection,3)                   !|
          s8 = Dtan (Disorient8*Dgra8)                               !|
          Call DVecSum2  (Base_Normal, 1.D0,                         !|
     *                    Projection,  s8,                           !|
     *                    Surface_Normal,3)                          !|
c Normalize to unit vector                                            |
          Call DUnitVec2 (Surface_Normal,Surface_Normal,3)           !|
c Verification:                                                       |
          s8 = DVecSca2 (Surface_Normal,Projection,3)                !|must be cos(90-Disor)=sin(Disor)
          if (Dabs(s8) .gt. 1D-8 ) Then !------------------+          |
             s8 = Dangle2 (Surface_Normal,Base_Normal,3)  !|          |must be disor wuthiut sign
     *          * s8 / Dabs(s8)                           !|          |correct the sign to what Disor should have
c If verification is not passed:                           |          |
            if (Dabs(s8-Disorient8) .gt. 1.0D0) goto 103 !-+----------+---+
          endif  !-----------------------------------------+          |   v
c Unit vector along internal surface normal (i.e. invert):            |
          Call   DVecCon (Surface_Normal,-1.D0,                      !|
     *                    Surface_Normal, 3)                         !|
        endif  !------------------------------------------------------+

  29    continue  !<---------------------------------------------------+
        if (iirezv.eq.1)  Then  !----------+                           |
          progname = progrezv(istackrezv) !|                           |
          istackrezv = istackrezv-1       !|                           |
        endif  !---------------------------+                           |
        return                                                        !|
                                                                      !|
c-----------------------------------------------------------           |
c                          ERRORS:                                     |
c-----------------------------------------------------------           |
  100   continue                                                      !|
        ifail = 100                                                   !|
        lines = 5                                                     !|
        write   (txt,200)  progname(1:lpn)                            !^
  200   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine surface normal:'//
     *          '** zero normal vector **')                           !^
        goto 300                                                      !|
                                                                      !|
  101   continue                                                      !|
        ifail = 101                                                   !|
        lines = 5                                                     !|
        write   (txt,201)  progname(1:lpn)                            !^
  201   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine surface normal:'//
     *          '** zero miscut vector **')                           !^
        goto 300                                                      !|
                                                                      !|
  102   continue                                                      !|
        ifail = 102                                                   !|
        lines = 5                                                     !|
        write   (txt,202)  progname(1:lpn)                            !^
  202   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine surface normal:'//
     *          '** miscut direction is along the base normal **')    !^
        goto 300                                                      !|
                                                                      !|
  103   continue                                                      !|
        ifail = 103                                                   !|
        lines = 8                                                     !|
        write (txt,203)  progname(1:lpn),                             !|
     *                   nBaseNormal, s8,                             !|
     *                   nBaseNormal, Disorient8                      !^
  203   format  (1x,a,'  E R R O R !!!'//
     *  ' Inconsistent result of building crystal surface normal.'//
     *  ' Computed normal angle to (',3i5,'): ',g10.3,' degr.'/
     *  ' Expected normal angle to (',3i5,'): ',g10.3,' degr.'//
     *  ' ** Please report the problem to the author! ** ')           !^
        goto 300                                                      !|
                                                                      !|
  300   continue                                                      !|
        Call    Message (txt,lines,2)                                 !|
        goto 29  !-----------------------------------------------------+
        end

c ======================================================================

        Subroutine SurfNormCoplanar (nhkl,       !hkl of Bragg planes (input)
     *                            Disorient8,    !Bragg planes miscut in degr. (input)
     *                            Surface_Normal,!hkl of surface normal (output)
     *                            ifail)         !failure flag (output)
        Real*8  rhkl(3),
     *          rany(3),
     *          orth(3),
     *          Surface_Normal(3),
     *          Disorient8,
     *          gra8, x1, x2
        Real*8   DVecMod, DVecMod2, DAngle2
        External DVecMod, DVecMod2, DAngle2

        Integer nhkl(3),
     *          ifail, lines, lpn, i

        Character       txt(20)*80
        Common  /msg/   txt

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
          progname = 'SurfNormCoplanar'   !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim(progname)
c-------------------------------------------------------
        ifail = 0
        gra8 = Datan(1.0D0)/45.         !number of radians in 1 degree

        x1 =  DVecMod(nhkl,3)
c If zero length of hkl:
        if (x1. lt. 0.01) goto 100

c Copy nhkl to real*8:
        Call DVecCopy_id (nhkl,rhkl,3)
c Normalize the length:
        Call DUnitVec2(rhkl,rhkl,3)

c Find any other vector which is non-parallel to rhkl:
c We probe: (h+1,k,l), (h,k+1,l)  and (h,k,l+1)
        do i=1,3  !====================================================+
          Call DVecCopy_dd (rhkl,rany,3)                              !|
          rany(i) = rany(i) + 1.0                                     !|
          if (abs(rany(i)) .lt. 0.01) rany(i) = rany(i) + 1.          !|
          Call DUnitVec2(rany,rany,3)                                 !|
          x1 = DAngle2(rhkl,rany,3)                                   !|
c If such a vector is found, build the orth to hkl:                    |
          if (x1.gt.3. .and. x1.lt.177.) Then !---------------+        |
c Calculate orth=[rhkl*rany]:                                !|        |
            Call DVecVec (rhkl,rany,orth,3)                  !|        |
            Call DUnitVec2(orth,orth,3)                      !|        |
            x1 = dcos(Disorient8*gra8)                       !|        |
            x2 = dsin(Disorient8*gra8)                       !|        |
c Calculate EXTERNAL surface normal:                         !|        |
c Surface_Normal = cos(pfi)*rhkl + sin(pfi)*orth:            !|        |
            Call DVecSum2 (rhkl,x1,orth,x2,Surface_Normal,3) !|        |
c Calculate INTERNAL surface normal:                         !|        |
            x1 = -1.                                         !|        |
            Call DVecCon (Surface_Normal,x1,Surface_Normal,3)!|        |
c Verify the surface normal length:                          !|        |
            x1 =  DVecMod2(Surface_Normal,3)                 !|        |
            if (x1 .lt. 0.99 .or. x1 .gt. 1.01)              !|        |
     *               stop 'SurfNormCoplanar: internal error' !|        |
            Call DUnitVec2(Surface_Normal,Surface_Normal,3)  !|        |
            goto 29 !----------------+                       !|        |
          endif !--------------------+------------------------+        |
        enddo  !=====================+=================================+
c Orth is not found:                 |
        goto 200                    !|
                                    !v
  29    continue  !<--------------<------------<---------------<-------+
        if (iirezv.eq.1)  Then  !----------+                           |
          progname = progrezv(istackrezv) !|                           |
          istackrezv = istackrezv-1       !|                           |
        endif  !---------------------------+                           |
        return                                                        !|
c-----------------------------------------------------------           |
c                          ERRORS:                                     |
c-----------------------------------------------------------           |
  100   continue                                                      !|
        ifail = 101                                                   !|
        lines = 5                                                     !|
        write   (txt,101)  progname(1:lpn)                            !^
  101   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine coplanar surface normal:'//
     *          '** zero Bragg reflection indices **')
        goto 900

  200   continue                                                      !|
        ifail = 102                                                   !|
        lines = 5                                                     !|
        write   (txt,201)  progname(1:lpn)                            !^
  201   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine coplanar surface normal:'//
     *          '** error building orthogonal coordinate system **')
        goto 900

  900   continue                                                      !^
        Call    Message (txt,lines,2)                                 !|
        goto 29  !-----------------------------------------------------+
        end

c ======================================================================

        Subroutine FindDisorient (QB,         !Bragg angle in degrees (input)
     *                            TER_Angle,  !TER critical angle in degrees (input)
     *                            Asymmetry,  !Disorient/incident/exit/beta (input)
     *                            UnitCoef,   !Units2 factor for Asymmetry, if angle (input)
     *                            igie,       !flag to interprete Asymmetry (input)
     *                            Disorient8, !Bragg planes misorientation (output)
     *                            ifail)      !failure flag  (output)
c ----------
c This program determines the angle between the surface and the Bragg planes,
c i.e. the disorientation of the Bragg planes with respect to the surface.
c We consider this angle is positive if g0>gh, i.e. if the incidence angle
c increases.
c ----------
        Real*8  Disorient8, gra8, Qi8, Qh8, beta8, tan8
        Real*4  QB, TER_Angle, Asymmetry, UnitCoef, Unit2degree
        Integer igie, ifail, lines, lpn

        Character       txt(20)*80
        Common  /msg/   txt

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
          progname = 'FindDisorient'      !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim(progname)
c-------------------------------------------------------
        ifail       = 0
        gra8        = Datan(1.0D0)/45.          !number of radians in 1 Degree
        Unit2degree = sngl(UnitCoef/gra8)       !recalc.factor from Unit to Degr.

        if (QB .le. 0. .OR. QB. gt. 90.) goto 100  !-------------+
                                                                !v
        if     (igie .eq. 6) Then !--------------------------+
c  [6]=coplanar via Bragg planes angle to surface            |
          Disorient8 = dble(Asymmetry) * Unit2degree        !|Unit2 -> Degr.
          if (Disorient8 .lt. -90. .OR.                     !|
     *        Disorient8 .gt. +90.) goto 700 !---------------+---+
                                                            !|   v
        elseif (igie .eq. 7) Then !--------------------------|
c  [7]=coplanar via incidence angle                          |
          Qi8 = dble(Asymmetry) * Unit2degree               !|Unit2 -> Degr.
          if (Qi8 .lt. -10.*TER_Angle .OR.                  !|
     *        Qi8 .gt. +90.)        goto 800 !---------------+---+
          Disorient8 = Qi8 - QB                             !|   v
                                                            !|
        elseif (igie .eq. 8) Then !--------------------------|
c  [8]=coplanar via exit angle                               |
          Qh8 = dble(Asymmetry) * Unit2degree               !|Unit2 -> Degr.
          if (Qh8 .lt. -10.*TER_Angle .OR.                  !|
     *        Qh8 .gt. +90.)        goto 850 !---------------+---+
          Disorient8 = QB - Qh8                             !|   v
                                                            !|
        elseif (igie .eq. 9) Then !--------------------------|
c  [9]=coplanar via beta=g0/|gh|                             |
          beta8 = dble(Asymmetry)                           !|
          if (beta8 .lt. 0. .OR. beta8 .gt. 1.D5) goto 200 !-+---+
          if (abs(beta8-1.).lt.1.E-20) Then !-----------+    |   v
            Disorient8 = 0.                            !|    |
          else !----------------------------------------+    |
            if (QB .ge. 90.00) goto 300 !---------------+----+---+
            tan8 = tan(QB*gra8)*(beta8-1.)/(beta8+1.)  !|    |   v
            Disorient8 = datan(tan8) / gra8            !|    |Degr.
          endif  !--------------------------------------+    |
                                                            !|
        else !-----------------------------------------------|
c  Unknown option:                                           |
          goto 400 !-----------------------------------------+---+
        endif  !---------------------------------------------+   v
          Qi8 = QB + Disorient8         !in degr.
          Qh8 = QB - Disorient8         !in degr.
c This principally can be commented out
c since it is checked again in Build_K0:
          if (Qi8 .le. -10.0001*TER_Angle) goto 500 !------------+
          if (Qh8 .le. -10.0001*TER_Angle) goto 600 !------------|
                                                                !v
  29    continue  !<--------------<------------<---------------<-------+
        if (iirezv.eq.1)  Then  !----------+                           |
          progname = progrezv(istackrezv) !|                           |
          istackrezv = istackrezv-1       !|                           |
        endif  !---------------------------+                           |
        return                                                        !|
c-----------------------------------------------------------           |
c                          ERRORS:                                     |
c-----------------------------------------------------------           |
  100   continue                                                      !|
        ifail = 100                                                   !|
        lines = 5                                                     !|
        write   (txt,101)  progname(1:lpn), QB                        !^
  101   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine Bragg planes disorientation:'//
     *          'incorrect Bragg angle QB=',f6.2,' degr.')            !|
        goto 900                                                      !|
                                                                      !|
  200   continue                                                      !|
        ifail = 200                                                   !|
        lines = 5                                                     !|
        write   (txt,201)  progname(1:lpn), beta8                     !^
  201   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine Bragg planes disorientation:'//
     *          'beta=g0/gh=',e8.2,' not in the range [0.-- 1.E5]')   !|
        goto 900                                                      !|
                                                                      !|
  300   continue                                                      !|
        ifail = 300                                                   !|
        lines = 5                                                     !|
        write   (txt,301)  progname(1:lpn)                            !^
  301   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine Bragg planes disorientation:'//
     *          'must specify beta=1 with QB=90 degr.')               !|
        goto 900                                                      !|
                                                                      !|
  400   continue                                                      !|
        ifail = 400                                                   !|
        lines = 5                                                     !|
        write   (txt,401)  progname(1:lpn), igie                      !^
  401   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine Bragg planes disorientation:'//
     *          'input_way =',i3,' not in the range [6-9]')           !|
        goto 900                                                      !|
                                                                      !|
  500   continue                                                      !|
        ifail = 500                                                   !|
        lines = 9                                                     !|
        write   (txt,501)  progname(1:lpn),                           !|
     *                     Qi8, -10.*TER_Angle                        !^
  501   format  (1x,a,':'//
     *          'Unpractical input found while calculating the Bragg',
     *                                     1x,'planes disorientation:'//
     *          'The incidence angle of x-rays corresponding to the',
     *                                            1x,'Bragg condition'/
     *          'appears to have too large negative value of f0=',f8.3,
     *                                                        ' degr.'/
     *          'being higher than  -10*f_critical=',f8.3,' degr.'//
     *          'This will result in negligible reflected intensity') !|
        goto 900                                                      !|
                                                                      !|
  600   continue                                                      !|
        ifail = 600                                                   !|
        lines = 9                                                     !|
        write   (txt,601)  progname(1:lpn),                           !|
     *                     Qh8, -10.*TER_Angle                        !^
  601   format  (1x,a,':'//
     *          'Unpractical input found while calculating the Bragg',
     *                                     1x,'planes disorientation:'//
     *          'The exit angle of x-rays corresponding to the',
     *                                            1x,'Bragg condition'/
     *          'appears to have too large negative value of fh=',f8.3,
     *                                                        ' degr.'/
     *          'being higher than  -10*f_critical=',f8.3,' degr.'//
     *          'This will result in negligible reflected intensity') !^
        goto 900                                                      !|
                                                                      !|
  700   continue                                                      !|
        ifail = 700                                                   !|
        lines = 7                                                     !|
        write   (txt,701)  progname(1:lpn), Disorient8,               !|
     *                     Asymmetry                                  !^
  701   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine Bragg planes disorientation:'//
     *          'Specified angle of Bragg planes to the surface = ',
     *                                               g10.3,' degr.'/
     *          '(your input of ',g10.3,' converted to degr.)'/
     *          'is not in the range [-90 : +90] degr.')              !|
        goto 900                                                      !|
                                                                      !|
  800   continue                                                      !|
        ifail = 800                                                   !|
        lines = 7                                                     !|
        write   (txt,801)  progname(1:lpn), Qi8, Asymmetry,           !|
     *                    -10.*TER_Angle                              !^
  801   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine Bragg planes disorientation:'//
     *          'Specified incidence angle  = ',g10.3,' degr.'/
     *          '(your input of ',g10.3,' converted to degr.)'/
     *          'is not in the range [',g9.3,' : +90] degr.')         !|
        goto 900                                                      !|
                                                                      !|
  850   continue                                                      !|
        ifail = 850                                                   !|
        lines = 6                                                     !|
        write   (txt,851)  progname(1:lpn), Qh8, Asymmetry,           !|
     *                     -10.*TER_Angle                             !^
  851   format  (1x,a,'  E R R O R !!!'//
     *          'Unable to determine Bragg planes disorientation:'//
     *          'Specified exit angle  = ',g10.3,' degr.'/
     *          '(your input of ',g10.3,' converted to degr.)'/
     *          'is not in the range [',g9.3,' : +90] degr.')         !|
        goto 900                                                      !|
                                                                      !|
  900   continue                                                      !|
        Call    Message (txt,lines,2)                                 !|
        goto 29  !-----------------------------------------------------+
        end

