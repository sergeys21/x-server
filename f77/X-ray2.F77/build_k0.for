        Subroutine Build_k0(Surface_Normal8,! input vector(3)
     *                      DisorPlanes8,   ! input scalar [angle of RecLatVec to IntSurfNorm in UnitCoef_IE; > 90degr.]
     *                      h_vector8,      ! input vector(3)
     *                      h8,             ! input scalar [length of h_vector]
     *                      Pfi4,           ! input scalar [angle of h to surface in UnitCoef_IE]
     *                      Psi4,           !OUTput scalar [Psi= 2*sin(QB)*sin(Pfi)]
     *                      PsiEff,         !OUTput scalar [PsiEff= 2*sin(QB)*Pfi] - for reference only
     *                      k0_vector8,     !OUTput vector(3)
     *                      k08,            ! input scalar [length of k0_vector]
     *                      Gam08,          !OUTput scalar [Gam0=sin(f0_centre)]
     *                      QB4,            ! input scalar [Bragg angle in degr]
     *                      SiQB8,          !OUTput scalar [sinus of Bragg angle]
     *                      fcentre4,       ! input scalar [f0_centre or fh_centre in UnitCoef_IE]
     *                      f0_centre4,     !OUTput scalar [in UnitCoef_IE]
     *                      fh_centre4,     !OUTput scalar [in UnitCoef_IE]
     *                      TER_Angle4,     ! input scalar [in UnitCoef_IE]
     *                      UnitCoef_IE,    ! input scalar [for incident/exit angles, i.e. UnitCoef(2)]
     *                      AngDegree,      ! input scalar
     *                      nBaseNormal,    ! input vector(3) [external base surface normal]
     *                      nReperDirection,! input vector(3)
     *                      igie,           ! input scalar [geometry specification, range 1-7]
     *                      ProgramVer,     ! input char*8
     *                      ipv,            ! input scaler [length of ProgramVer]
     *                      ifail)          !OUTput scalar
c######################################################################
c Build the incident wave vector for:    gid_sl50,
c                                        gid_sl52,
c                                        gid_sl55,
c                                        gid_sl97,
c                                        gid_sl98,
c                                        gid_sl99,
c                                        gid_slm1,
c                                        gids_10,
c                                        gids_97,
c                                        gids_98.
c######################################################################
c igie = Geometry specification:
c 1 - noncoplanar case via the incidence angle of k0 --> fcentre=incidence angle
c 2 - noncoplanar case via the exit angle of kh      --> fcentre=exit angle
c 3 - coplanar case; grazing incidence               --> fcentre=0 (not used)
c 4 - coplanar case; grazing exit                    --> fcentre=0 (not used)
c 5 - symmetric Bragg case (not GID)                 --> fcentre=0 (not used)
c 6 - coplanar & Bragg planes disorientation --> fcentre=incidence angle (see gid_IN -> FindDisorient)
c 7 - coplanar & incidence angle             --> fcentre=incidence angle (see gid_IN -> FindDisorient)
c 8 - coplanar & exit angle of kh            --> fcentre=incidence angle (see gid_IN -> FindDisorient)
c 9 - coplanar & beta=g0/|gh|                --> fcentre=incidence angle (see gid_IN -> FindDisorient)
c----------------------------------------------------------------------
        Real*8          Surface_Normal8(3), DisorPlanes8,
     *                  DisorPlanesExt8,
     *                  h_vector8(3), h8, hn8, h_x_n8(3),
     *                  k0_vector8(3), k08, Gam08, SiQB8, CoQB8,
     *                  c0, c1, c2, c3, b0, b1, b2, b3, DDD,
     *                  wrk_vector8(3), wr2_vector8(3)

        Real*8          DVecSca2, DAngle2, DvecMod2
        External        DVecSca2, DAngle2, DvecMod2

        Real*8          Range_IE(2)     ! allowed range of incidence/exit
        Common  /IERAN/ Range_IE        ! angles in degrees (determined here)

        Real*8          f0_centre8, fh_centre8, Pfi8, Psi8, QB8,
     *                  f0_centr_B, f0_centr_L, Allow_IE, Allow_Degr,
     *                  fh_centr_B, fh_centr_L

        Real*4          QB4, TER_Angle4, fcentre4, Pfi4, Psi4,
     *                  f0_centre4, fh_centre4, PsiEff,
     *                  UnitCoef_IE, AngDegree

        Integer         nBaseNormal(3), nReperDirection(3),
     *                  igie, ipv, ifail, i, j, leunit

        Character       ProgramVer*8, UnitName*5

        Character       txt(20)*80
        Common  /msg/   txt
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
          progname = 'Build_k0'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
c -------------------------------------------------------
        ifail      =  0
c Get the name of angular units:
        Call GetUname (UnitCoef_IE, UnitName)
        leunit     = Max (Len_Trim(UnitName),1)

        Allow_IE   =  10.*Dble(TER_Angle4)              ! in UnitCoef_IE
        Allow_Degr =  Allow_IE*UnitCoef_IE/AngDegree    ! convert to degrees
c cccc  k08   = 2.*pi/wave
        Pfi8 =  Dble(Pfi4)                              ! in UnitCoef_IE
        QB8  =  Dble(QB4)                                ! in degr.
        SiQB8 = Dsin(QB8*AngDegree)
        CoQB8 =  Dcos(QB8*AngDegree)
c Angle of reciprocal lattice vector
c to external normal in degrees:
c (== angle of Bragg planes to surface)
        DisorPlanesExt8 = 180. - DisorPlanes8*UnitCoef_IE/AngDegree

c Normalize h_vector:
        Call    DUnitVec2 (h_vector8,h_vector8,3)
        h8 = 2.*k08*SiQB8
        Call DVecCon (h_vector8,h8,h_vector8,3)

c Compute scalar product (h_vector * Surface_Normal):
        hn8 = DVecSca2 (h_vector8,Surface_Normal8,3)

c Compute vector product [h_vector * Surface_Normal]:
        Call DVecVec (h_vector8,Surface_Normal8,h_x_n8,3)

c Compute [h_vector * Surface_Normal]^2:
        c0 = h8*h8 - hn8*hn8

c Compute [h_vector * Surface_Normal]^2 (testing another way):
c ccccc c1 = DVecSca2 (h_x_n,h_x_n,3)

        if (c0 .gt. (1.0D-8)*(h8*h8))  Then !------------------------+
c                                       +--------------------------+ |
c                                       | Reciprocal vector is NOT | |
c                                       | perpendicular to surface | |
c                                       +==========================+ |
c Effective miscut of Bragg planes to internal surface NORMAL        |
c                     -- OR of h_vector to the surface PLANE:        |
c ('+'=IN, '-'=OUT):                                                 |
          Psi8 = 2.*Dsin(QB8*AngDegree)                             !|
     *             *Dsin(Pfi8*UnitCoef_IE)                          !|
c+--------------------------+                                        |
c|Determine a possible range|                                        |
c|of incidence/exit angles: |                                        |
c+--------------------------+                                        |
c The angle between <h_vector> and the surface:                      |
          DDD = DisorPlanes8*UnitCoef_IE/AngDegree - 90.0D0  !in Degr|
c Minimum incidence:                                                !|
          Range_IE(1) = DDD-(90.D0-QB8)                      !in Degr|
c Maximum incidence:                                                 |
          Range_IE(2) = DDD+(90.D0-QB8)                      !in Degr|
          do    i=1,2   !==========================================+ |
ccc         if (Range_IE(i) .lt. 0.D0)                            !| |
ccc  *                             Range_IE(i)=0.D0               !| |
            if (Range_IE(i) .lt.-Allow_Degr)                      !| |
     *                             Range_IE(i)=-Allow_Degr        !| |
            if (Range_IE(i) .gt. 90.D0)                           !| |
     *                             Range_IE(i)=180.D0-Range_IE(i) !| |
          enddo  !=================================================+ |
          c1 = Min(Range_IE(1),Range_IE(2))                         !|
          c2 = Max(Range_IE(1),Range_IE(2))                         !|
          Range_IE(1) = c1                                          !|
          Range_IE(2) = c2                                          !|
          if (igie .lt. 3) Then  !-----------------+                 |
c Specified incidence/exit angle in Degr.:         |                 |
            c1 = fcentre4*UnitCoef_IE/AngDegree   !|                 |
            if (c1 .lt. Range_IE(1) .Or.          !|                 |
     *          c1 .gt. Range_IE(2)) goto 182 !----+-----------------+--+
          endif  !---------------------------------+                 |  v
c+-------------------+                                               |
c|Find the incidence |                                               |
c|and exit angles:   |                                               |
c+-------------------+                                               |
          if (igie .eq. 1) Then !-Via incidence angle of k0--+       |
                                                            !|       |
            f0_centre8 = Dble(fcentre4)                     !|       |
c Here Psi8= 2*sin(QB)*sin(Pfi):                            !|       |
            fh_centr_B =-Psi8 - Dsin(f0_centre8*UnitCoef_IE)!|Bragg  |
            fh_centr_L = Psi8 + Dsin(f0_centre8*UnitCoef_IE)!|Laue   |
            if (fh_centr_B .gt. 0.) Then !---------------+   |       |
              fh_centre8 = Dasin(fh_centr_B)/UnitCoef_IE!|   |Bragg  |in Units
            else  !--------------------------------------|   |       |
              fh_centre8 = Dasin(fh_centr_L)/UnitCoef_IE!|   |Laue   |in Units
              if (fh_centre8 .gt. Allow_IE) goto 162 !---|->-+-+     |
            endif  !-------------------------------------+   | v     |
                                                            !|       |
          elseif (igie.eq.2) Then !----Via exit angle of kh--|       |
                                                            !|       |
            fh_centre8 = Dble(fcentre4)                     !|       |
c Here Psi8= 2*sin(QB)*sin(Pfi):                            !|       |
            f0_centr_B =-Psi8 - Dsin(fh_centre8*UnitCoef_IE)!|Bragg  |
            f0_centr_L =-Psi8 + Dsin(fh_centre8*UnitCoef_IE)!|Laue   |
            if (f0_centr_B .gt. 0.) Then !---------------+   |       |
              f0_centre8 = Dasin(f0_centr_B)/UnitCoef_IE!|   |Bragg  |in Units
            else  !--------------------------------------|   |       |
              f0_centre8 = Dasin(f0_centr_L)/UnitCoef_IE!|   |Laue   |in Units
c cccccccccc  if (f0_centre8 .lt.  0.D0)    goto 172 !---+->-+-+     |
              if (f0_centre8 .lt.-Allow_IE) goto 172 !---+->-+-|     |
              if (f0_centre8 .gt. Allow_IE) goto 172 !---+->-+-|     |
            endif  !-------------------------------------+   | v     |
                                                            !|       |
          elseif (igie .eq. 3) Then !-Coplanar-graz.incidence|       |
                                                            !|       |
            c1 = DDD-(90.D0-QB8)                            !|       |
            c2 = DDD+(90.D0-QB8)                            !|       |
            f0_centre8 = Min(c1,c2)                         !|       |
            fh_centre8 = Max(c1,c2)                         !|       |
            if (fh_centre8 .gt. 90.D0)                      !|       |
     *                        fh_centre8=180.D0-fh_centre8  !|       |
c Convert from Degrees to input Units:                      !|       |
            f0_centre8 = f0_centre8*AngDegree/UnitCoef_IE   !|       |
            fh_centre8 = fh_centre8*AngDegree/UnitCoef_IE   !|       |
            if (f0_centre8 .lt.-Allow_IE) goto 193    !------+-+     |
                                                            !| v     |
          elseif (igie .eq. 4) Then !---Coplanar graz.exit---|       |
                                                            !|       |
            c1 = DDD-(90.D0-QB8)                            !|       |
            c2 = DDD+(90.D0-QB8)                            !|       |
            fh_centre8 = Min(c1,c2)                         !|       |
            f0_centre8 = Max(c1,c2)                         !|       |
            if (f0_centre8 .gt. 90.D0)                      !|       |
     *                        f0_centre8=180.D0-f0_centre8  !|       |
c Convert from Degrees to input Units:                      !|       |
            f0_centre8 = f0_centre8*AngDegree/UnitCoef_IE   !|       |
            fh_centre8 = fh_centre8*AngDegree/UnitCoef_IE   !|       |
            if (fh_centre8 .lt.-Allow_IE) goto 193    !------+-+     |
                                                            !| v     |
          elseif (igie .eq. 5) Then !--------Symmetric case--|       |
                                                            !|       |
            if (Psi8.ge.0.) goto 194 !-----------------------+-+Laue!|
            f0_centre8 = Dasin(-Psi8/2.) / UnitCoef_IE      !| v     |
            fh_centre8 = f0_centre8                         !|       |
                                                            !|       |
          elseif (igie .ge. 6 .and.                         !|       |
     *            igie .le. 9) Then !----Coplanar-case-------|       |
                                                            !|       |
c max & min incidence angles in coplanar case:               |       |
            c1 = DDD-(90.D0-QB8)                            !|       |
            c2 = DDD+(90.D0-QB8)                            !|       |
            if (c1 .gt. 90.D0) c1 = 180.-c1                 !|       |
            if (c2 .gt. 90.D0) c2 = 180.-c2                 !|       |
c Convert from Degrees to input Units:                      !|       |
            c1 = c1*AngDegree/UnitCoef_IE                   !|       |
            c2 = c2*AngDegree/UnitCoef_IE                   !|       |
c 'fcentre' must be the incidence angle (input parameter).   |       |
c Select, what is closer to it:                              |       |
            if (abs(fcentre4-c1) .lt.                       !|       |
     *          abs(fcentre4-c2)) Then !-+                   |       |
              f0_centre8 = c1           !|                   |       |
              fh_centre8 = c2           !|                   |       |
            else !-----------------------|                   |       |
              f0_centre8 = c2           !|                   |       |
              fh_centre8 = c1           !|                   |       |
            endif  !---------------------+                   |       |
c Check for error (we expect that the                        |       |
c difference must not exceed 1 degr.):                       |       |
            c3 = abs(f0_centre8-fcentre4)                   !|       |
     *         * UnitCoef_IE / AngDegree                    !|       |
            if (c3 .gt. 1.) goto 195 !-----------------------+-+     |
                                                            !| v     |
c Check that we are not below -10*TERngle:                   |       |
            if (f0_centre8 .lt. -Allow_IE) goto 193 !--------+-+     |
                                                            !| v     |
          else !---------------------------------------------|       |
                                                            !|       |
            goto 192 !-------------unknown mode--------------+-+     |
                                                            !| v     |
          endif  !-------------------------------------------+       |
                                                                    !|
          Gam08 =  Dsin(f0_centre8*UnitCoef_IE)                     !|
                                                                    !|
c+==============================================================+    |
c|We present the incident wave vector in the form:              |    |
c|(this only works when they are not parallel - checked above)  |    |
c|                                                              |    |
c|->               ->             ->             ->             |    |
c|k0 = c1*k0*Surface_Normal + c2*h_vector + c3*h_x_n            |    |
c|                                                              |    |
c|Then, using:                                                  |    |
c| (k0*n)  = k0*sin(f0)   = c1*k0 + c2*(h*n)                    |    |
c| (k0*h)  =-k0*h*sin(QB) = c1*k0*(h*n) + c2*h^2                |    |
c| (k0*k0) = k0^2         = c1^2*k0^2 + c2^2*h^2 + c3^2*[hxn]^2 |    |
c|                        + 2*c1*c2*k0*(h*n)                    |    |
c|-- We determine the coefficients c1,c2,c3.                    |    |
c|                                                              |    |
c|ATTENTION: Surface_Normal is the INTERNAL normal,             |    |
c|            contrary to nBaseNormal !                         |    |
c|==============================================================|    |
c+==============================================================+    |
c ccccccc c0 = h*h - hn*hn = [h*n]^2                                !|
          c1 = h8*(h8*Gam08 + hn8*SiQB8) / c0                       !|
          c2 =-k08*(h8*SiQB8 + hn8*Gam08) / c0                      !|
          b1 = k08*k08*(1.-c1*c1)                                   !|
          b2 = h8*h8*c2*c2                                          !|
          b3 = 2.*c1*c2*k08*hn8                                     !|
c Calculate c3:                                                      |
          c3 = ( b1 - b2 - b3) / c0                                 !|
          if (c3 .le. 0.0D0) Then !---------------+                  |
c Estimate possible error:                        |                  |
          b0 = Max(abs(b1),abs(b2),abs(b3)) / c0 !|                  |
c This is allowance for numerical errors:         |                  |
            if (c3 .lt. (-1.0D-8)*b0) goto 163 !--+-->---Error-------+-+
            c3 = 0.0D0                           !|                  | v
          else  !---------------------------------|                  |
            c3 =  Dsqrt (c3)                     !|                 !|
          endif  !--------------------------------+                  |
          Call  DVecSum2 (Surface_Normal8,c1*k08,                   !|
     +                    h_vector8,      c2,                       !|
     =                    k0_vector8,     3)                        !|
          Call  DVecSum2 (k0_vector8,     1.D0,                     !|
     +                    h_x_n8,         c3,                       !|
     =                    k0_vector8,     3)                        !|
c Testing:                                                           |
          c1   = DVecSca2 (k0_vector8,k0_vector8,3)                 !|
          c2   = k08*k08                                            !|
          c3   = 100. * Abs(c2-c1)/c2                               !|
          if (c3 .gt. 0.1) goto 164   !---------->------Error--------+-+
c                                                                    | v
        else  !-----------(c0 .le. (1.0D-8)*(h8*h8))-----------------|
c                                                                    |
c                      +--------------------------+ +------------+   |
c                      | Reciprocal vector **is** | |  Symmetric |   |
c                      | perpendicular to surface | | Bragg case |   |
c                      +==========================+ +============+   |
          Pfi8      = -90.D0 * AngDegree / UnitCoef_IE              !|
c Effective misorientation of Bragg planes to surface NORMAL         |
c                    -- OR of h_vector to the surface PLANE:         |
c         Psi8      = 2.*Dsin(QB8*AngDegree)                        !|
c    *                  *Dsin(Pfi8*UnitCoef_IE)                     !|
c Modification of July 9, 1998:                                      |
          Psi8      =-2.*Dsin(QB8*AngDegree)                        !|
c -- since:              Dsin(Pfi8*UnitCoef_IE) = -1                !|
                                                                    !|
c Find the incidence                                                 |
c and exit angles:                                                   |
          if (igie .lt. 5) Then !----------------------------+       |
            write (txt,141)  ProgramVer(1:ipv),igie         !v       v
  141       format (a//
     *      'W A R N I N G :'//
     *      'The Bragg planes are parallel to the surface'/
     *      'Proceeding from mode=',i1,
     *                 ' to symmetric Bragg case (mode=5).'//
     *      ' Hit any key to continue now...')              !^       ^
            if (modebat .ne. 0) Then !-----------+           |       |
              do i=1,6  !=====================+  |           |       |
                j = Len_Trim(txt(i))         !|  |           |       |
                if (j .lt. 1) j=1            !|  |           |       |
                write (3,'(a)') txt(i)(1:j)  !|  |           |       |
              enddo  !========================+  |           |       |
            else  !------------------------------|           |       |
              Call Message (txt,8,0)            !|           |       |
            endif  !-----------------------------+           |       |
            igie = 5                                        !|       |
          endif  !-------------------------------------------+       |
          f0_centre8 = QB8 * AngDegree / UnitCoef_IE                !|
          Gam08      = Dsin(f0_centre8*UnitCoef_IE)                 !|
c ccccccc fh_centre8 = Abs(Psi8) - Gam0                             !|
          fh_centre8 = Dasin(-Psi8-Gam08) / UnitCoef_IE             !|
c+=====================================================+             |
c|Vector h is parallel to the <surface normal>         |             |
c|Use another expansion (decomposition) of k_0:        |             |
c|->               ->             ->   ->              |             |
c|k0 = c1*k0*Surface_Normal + c2*k0*Projection  ,      |             |
c|                                                     |             |
c|where Projection is a projection of ReperDirection to|             |
c|the surface (to the plane normal to Surface_Normal)  |             |
c|                                                     |             |
c|Then, using (k0*Surface_Normal) = k0*sin(QB) = k0*c1 |             |
c|            (k0*Projection)     = k0*cos(QB) = k0*c2 |             |
c|                                                     |             |
c|-- We determine the coefficients c1,c2:              |             |
c+=====================================================+             |
          do i=1,3  !==========================================+     |
            wr2_vector8(i) = Dfloat(nReperDirection(i))       !|     |
     +                    + Dfloat(nBaseNormal(i))            !|     |
c (the second line avoids warning about not using nBaseNormal) |     |
          enddo  !=============================================+     |
          c3 = DvecMod2 (wr2_vector8,3)                             !|
c If the chosen vector is zero, then select ( 1 2 3):                |
          if (c3 .lt. 1.0D-5) Then  !--------+                       |
            do i=1,3   !==================+  |                       |
              wr2_vector8(i) = Dfloat(i) !|  |                       |
            enddo  !======================+  |                       |
          endif  !---------------------------+                       |
          Call DUnitVec2  (wr2_vector8,                             !|
     *                     wr2_vector8,3)                           !|
          Call  DProject2 (Surface_Normal8,                         !|
     *                     wr2_vector8,                             !|
     *                     wrk_vector8,3)                           !|
          c3 = DvecMod2 (wrk_vector8,3)                             !|
c If the chosen vector | Surface_Normal (the length of the           |
c projection wrk_vector is zero), then just select ANY other         |
c vector which is not zero and not parallel to Surface_Normal:       |
c (here 'ABS' guaranties that wr2_vector(3)#0)                       |
          if (c3 .lt. 1.0D-5) Then  !------------------+             |
            wr2_vector8(1) = Abs(wr2_vector8(1))+10.  !|Arbitrary    |
            wr2_vector8(2) = Abs(wr2_vector8(2))+3.14 !|numbers!!    |
            wr2_vector8(3) = Abs(wr2_vector8(3))+17.  !|             |
          endif  !-------------------------------------+             |
          Call  DProject2 (Surface_Normal8,                         !|
     *                     wr2_vector8,                             !|
     *                     wrk_vector8,3)                           !|
          Call  DUnitVec2 (wrk_vector8,                             !|
     *                     wrk_vector8,3)                           !|
          Call  DVecSum2  (Surface_Normal8, SiQB8*k08,              !|
     +                     wrk_vector8,     CoQB8*k08,              !|
     =                     k0_vector8,      3)                      !|
        endif  !-----------------------------------------------------+
c cc The 1st version does not work because it can be: Psi > 1.
c cc    Psi8    = 2.*Dsin(QB8*AngDegree) * Dsin(Pfi8*UnitCoef_IE)
c cc    PsiEff  = Dasin(Psi8) / UnitCoef_IE
c cc That is used solely for printout:
        PsiEff  = sngl (2.*Dsin(QB8*AngDegree) * Pfi8)

  29    continue  !<---------------------<-----------------------+
        f0_centre4 = Sngl(f0_centre8)            !in UnitCoef_IE |
        fh_centre4 = Sngl(fh_centre8)            !in UnitCoef_IE |
        Pfi4      = Sngl(Pfi8)                   !in UnitCoef_IE |
        Psi4      = Sngl(Psi8)                   !no units!      |
        if (iirezv .eq. 1) Then !----------+                    !^
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c #######################################################
c                  Error messages
c#######################################################
  162   continue
        write   (txt,121)  ProgramVer(1:ipv),
     *                     fh_centre8, UnitName(1:leunit),
     *                     Allow_IE, UnitName(1:leunit)
  121   format  (/
     *  a,'  E R R O R :'/
     *  /
     *  'At the specified incidence angle, the diffracted'/
     *  'wave points too steep inside the crystal:'/
     *  'fh=',g10.4,a,' > 10*f_critical=',g10.4,a/
     *  'The reflected diffracted intensity at the'/
     *  'entrance surface of crystal will be negligible.'/
     *  /
     *  'This program does not handle the Laue (transmisssion) cases.'/
     *  'Please, decrease the incidence angle..')
        Call    Message (txt,11,2)
        goto 28

  172   continue
        write   (txt,131)  ProgramVer(1:ipv),
     *                     f0_centre8, UnitName(1:leunit),
     *                     -Allow_IE, UnitName(1:leunit)
  131   format  (/
     *  a,'  E R R O R :'/
     *  /
     *  'The specified exit angle corresponds'/
     *  'to a large negative incidence angle:'/
     *  'f0=',g10.4,a,' < -10*f_critical=',g10.4,a/
     *  /
     *  'Please, decrease the exit angle..')
        Call    Message (txt,8,2)
        goto 28

  182   continue
        write   (txt,142)  ProgramVer(1:ipv),
     *                     fcentre4, UnitName(1:leunit), c1,
     *                     Range_IE, DDD, (90.-QB4)
  142   format  (/
     *  a,' E R R O R :'/
     *  /
     *  'The specified incidence/exit angle:'/
     *   g15.8,a,' =',g15.8,'degr.'/
     *  'is outside the range of possible values:'/
     *  '[',g15.8,'::',g15.8,'] degr.'/
     *  '-- when the reciprocal lattice vector corresponding to the'/
     *  'Bragg planes makes the angle ',g15.8,'degr.'/
     *  'to the crystal surface and the Bragg cone opening is:'/
     *  '(90-QB)=',g15.8,'degr.'/
     *  /
     *  'Please, correct the incidence/exit angle.')
        if (UnitName(1:leunit) .eq. 'degr.') txt(5)(15+leunit+1:80)=' '
        Call    Message (txt,13,2)
        goto 28

  163   continue
        write   (txt,119)  ProgramVer(1:ipv), c3,
     *                     Surface_Normal8, h_vector8,
     *                     f0_centre8, UnitName(1:leunit),
     *                     fh_centre8, UnitName(1:leunit),
     *                     QB4
  119   format  (/
     *  a,': Incident wavevector determination error.'/
     *  /
     *  'The calculations resulted in negative value: [H*N]**2 =',g9.2/
     *  'when trying to build: k0=c1*N+c2*H+c3*[H*N]'/
     *  '  (N=surface normal = (',3g11.3,'),'/
     *  '   H=reciprocal lattice vector = (',3g11.3,')).'/
     *  /
     *  'Conditions: f0=',g10.4,a,'  fh=',g10.4,a,'  QB=',f5.2,'degr.'/
     *  /
     *  'Possible causes may be: unreachable incidence/exit angle,'/
     *  'or other contradictory input, or precision loss in calc.'/
     *  /
     *  'Please, verify your input, try another way to describe the'/
     *  'geometry, or communicate the problem to the author.')
        Call    Message (txt,15,2)
        goto 28

  164   continue
        write   (txt,889)  ProgramVer(1:ipv), c3
  889   format  (/
     *  a,': k0 determination error.'/
     *  /
     *  'Test results: |1-(k0*k0)/k0^2| =',g10.3,'% exceeds 0.1%'/
     *  /
     *  'Please, report the problem to the author.')
        Call    Message (txt,6,2)
        goto 28

  192   continue
        write   (txt,890)  ProgramVer(1:ipv), igie
  890   format  (/
     *  a,'  E R R O R :'/
     *  /
     *  'Unsupported mode=',i3,' of geometry specification'/
     *  /
     *  'Please, report the problem to the author.')
        Call    Message (txt,6,2)
        goto 28

  193   continue
        write   (txt,891)  ProgramVer(1:ipv),
     *                     DisorPlanesExt8, QB8,
     *                     f0_centre8, UnitName(1:leunit),
     *                     fh_centre8, UnitName(1:leunit),
     *                     -Allow_IE, UnitName(1:leunit)
  891   format  (/
     *  a,'  E R R O R :'/
     *  /
     *  'The geometry resulted from the input is:'/
     *  'The angle  of Bragg planes to  the surface =',f6.2,' degr.,'/
     *  'The Bragg angle for given x-ray reflection =',f6.2,' degr.'/
     *  /
     *  'Then, the coplanar mode corresponds to a Laue case:'/
     *  'either the incidence angle (f0 =',g10.4,a,') or the exit'/
     *  'angle (fh =',g10.4,a,') has a negative value that is more'/
     *  'than 10 times bigger than the critical angle for total'/
     *  'external reflection: -10*f_critical=',g10.4,a/
     *  /
     *  'The Laue cases are not handled by this program since it'/
     *  'assumes an infinitely thick crystal substrate while the'/
     *  'diffracted intensity yielding the entrance surface due'/
     *  'to the specular reflection effect will be negligible.')
        Call    Message (txt,17,2)
        goto 28

  195   continue
        write   (txt,893)  ProgramVer(1:ipv),
     *                     f0_centre8, UnitName(1:leunit),
     *                     fcentre4, UnitName(1:leunit)
  893   format  (/
     *  a,'  E R R O R :'/
     *  /
     *  'Cross-checking shows that two ways to build the K0 vector'/
     *  'result in two different incidence angles:'/
     *  /
     *  'f0_1=',g10.4,a,',  f0_2=',g10.4,a/
     *  /
     *  'Please, report the problem to the author.')
        Call    Message (txt,9,2)
        goto 28

  194   continue
        write   (txt,892)  ProgramVer(1:ipv)
  892   format  (/
     *  a,'  E R R O R :'/
     *  /
     *  'The reciprocal lattice vector corresponding to the Bragg'/
     *  'planes directed inwards or along the surface.'/
     *  /
     *  'Symmetric Bragg geometry cannot be built for this case.'/
     *  'This program does not handle Laue cases except for GID.'/
     *  'To build the GID, please, choose non-coplanar diffraction'/
     *  'geometry and specify a small incidence or exit angle.')
        Call    Message (txt,10,2)
        goto 28

  28    continue
        ifail = 1
        goto 29
        end
