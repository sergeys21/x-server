        Subroutine IncidentVector (
     *          Wave,      !x-ray wavelength (Angstrem)          - input/output  real*8
     *          nIndexes,  !number of Miller indices (3 or 4)    - input  integer
     *          nReflex1,  !dimension of ind_h,h_vec,h_mod       - input  integer (nReflect+1)
     *          ind_h,     !ind_h(nIndexes,nReflex1) - h-indeces - input  integer ---->--+
     *          h_vec,     !h_vec(nIndexes,nReflex1) - h-vectors - output real*8  ---->--|
     *          h_mod,     !h_mod(nReflex1) = |h_vec(i)|         - output real*8  ---->--|
     *          Zone,      !Zone(nIndexes)=[h_vec(1)*h_vec(2)]   - output real*8         |
     *          Zone_mod,  !=|Zone|                              - output real*8         |
     *          k0_copl,   !k0_copl(nIndexes) = coplanar k0      - output real*8         |
     *          Wave_copl, !wavelength corresponding to k0_copl  - output real*8         |
     *          Uncopl,    !non-coplanar component of k0         - output real*8         |
     *          k0,        !k0(nIndexes) incident wavevector     - output real*8         |
     *          k0_mod,    !=|k0|                                - output real*8         |
     *          ifail)     !completion flag                      - output integer        |
                                                                                        !|
c -- Order of vector in the array ind_h:  <----------------------------------------------+
c ind_h(1)   = (0,0,0)
c ind_h(2)   = h_1
c ind_h(3)   = h_2
c ind_h(4)   = h_3
c ................
c ind_h(n)   = h_(n-1)
c ind_h(n+1) = h_1 - h_2
c+==========================================================+
c| Determine incident wavevector for crystal with arbitrary |
c| symmetry when one of the three conditions are given:     |
c| 1) X-ray wavelength and indices of two reciprocal lattice|
c|    vectors (diffaction vectors)                          |
c| 2) Condition that the diffraction geometry is complanar  |
c|    indices of two reciprocal lattice vectors (diffaction |
c|    vectors)                                              |
c| 3) Indices of three  reciprocal lattice vectors          |
c|    (diffaction vectors)                                  |
c|                                                          |
c|                               S.Stepanov 1993-2003       |
c|----------------------------------------------------------|
c| This subroutine should be called from the multiple       |
c| diffraction programs like brls.
c|----------------------------------------------------------|
c|ATTENTION: The unit cell paramteters are expected to be   |
c|           already determined & stored in respective      |
c|           Common (e.g. by calling X0h) so that the       |
c|           reciprocal lattice is either defined or can    |
c|           by determined automatically by one of the      |
c|           linear algerbra routines (see algebraNN.for).  |
c|----------------------------------------------------------|
c|ATTENTION: there are several possible options with input  |
c|           wave:                                          |
c|   wave> 0. - normal input                                |
c|   wave=-1. - wave to be calculated as wave_coplanar      |
c|   wave=-2. - wave to be calculated from 3 reflection     |
c+==========================================================+
        Integer   nIndexes,
     *            nReflex1,
     *            ind_h(nIndexes,nReflex1),
     *            ifail, i, j
        Real*8    Wave, pi, c1, c2, ss,
     *            h_vec(nIndexes,nReflex1),
     *            h_mod(nReflex1),
     *            Zone(nIndexes),
     *            Zone_mod,
     *            k0_copl(nIndexes), k0_copl_mod,
     *            Wave_copl, Uncopl,
     *            k0(nIndexes),
     *            k0_mod

        Real*8    DVecMod2, DVecSca2, Dangle2
        External  DVecMod2, DVecSca2, Dangle2

        Integer   lpn, involved, iii(4)

        Character wrk*24

c Added 07/2003: this is internal normal. When the incident vector is
c specified via the wavelength, we have 2 choices (above and below the
c coplanar plane). The normal is used to select one of the two:
        Real*8    Surface_Normal(3),
     *            Surface_Norm(4)
        Common /SurfNorm/ Surface_Normal

        Character txt(20)*80
        Common /msg/ txt

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
          progname = 'IncidentVector'     !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim (progname)

c------------------------------------------------------------------
        pi = 4.*Datan(1.D0)
        ifail = 0

c Surface normal is given by 3 indices
c -- we may need to construct four:
        if (nIndexes .eq. 3) Then  !--------------+
          do i=1,3 !=============================+|
            Surface_Norm(i) = Surface_Normal(i) !||
          enddo  !===============================+|
        else  !-----------------------------------|
         Surface_Norm(1) = Surface_Normal(1)     !|
         Surface_Norm(2) = Surface_Normal(2)     !|
         Surface_Norm(3) =-(Surface_Normal(1)+   !|
     +                      Surface_Normal(2))   !|
         Surface_Norm(4) = Surface_Normal(3)     !|
        endif  !----------------------------------+

c Vector hn=(h1-h2) is needed for the calculation of
c non-coplanar component of the wavevector:
        do      i=1,nIndexes  !============================+
          ind_h (i,nReflex1) = ind_h (i,2) - ind_h (i,3)  !|
        enddo  !===========================================+

c Length of vectors h0=0, h1, h2,..., hn=(h1-h2) in 1/Angstrem:

        do i=1,nIndexes !====+          !0th reflection
          h_vec(i,1) = 0.   !|
        enddo !==============+
        h_mod(1) = 0.

        do      j=2,nReflex1  !===========================+
                                                         !|
          Call Build_h (ind_h(1,j),                      !|
     *                  h_vec(1,j),h_mod(j),             !|
     *                  nIndexes)                        !|
          if (abs(h_mod(j)).lt.1.E-20) Then !--+          |
            if (j .lt. nReflex1) Then !-+      |          |
              goto 51 !-----------------+------+-zero-----+----+
            else !----------------------+      |          |    v
              goto 59 !-----------------+------+-parallel-+----+
            endif !---------------------+      |          |    v
          endif !------------------------------+          |
                                                         !|
        enddo  !==========================================+

c Zone=[h1*h2]:
        Call DVecVec (h_vec(1,2),h_vec(1,3),Zone,nIndexes)
c Zone_mod =|[h1*h2]|:
        Zone_mod = DVecMod2 (Zone,nIndexes)
c If h1 and h2 are parallel or anti-parallel:
        if (Zone_mod/(h_mod(2)*h_mod(3)) .lt. 1.E-5) goto 59 !-+
c+========================================================+    v
c|Incident wavevector in coplanar case:                   |
c|                                                        |
c|                ->          ->      ->                  |
c|                k(copl)= c1*h1 + c2*h2                  |
c|                                                        |
c| Then, the Bragg conditions:                            |
c|                                                        |
c|                  ->      ->     2                      |
c|                2*k(copl)*h1 + h1  = 0                  |
c|                                                        |
c|                  ->      ->     2                      |
c|                2*k(copl)*h2 + h2  = 0                  |
c|                                                        |
c| give:                                                  |
c|           2            2  ->     2            2  ->    |
c|         h2 *((h1*h2)-h1 )*h1 + h1 *((h1*h2)-h2 )*h2    |
c|k(copl)= -------------------------------------------    |
c|                     2*([h1*h2]*[h1*h2])                |
c+========================================================+
        ss =  DVecSca2 (h_vec(1,2),h_vec(1,3),nIndexes)
        c1 = (h_mod(3)**2) * (ss - h_mod(2)**2)
        c2 = (h_mod(2)**2) * (ss - h_mod(3)**2)
        ss = 2.*Zone_mod**2
        c1 = c1/ss
        c2 = c2/ss
c k0_copl - shortest wavevector corresponding to coplanar case:
        Call    DVecSum2 (h_vec(1,2),c1,
     *                    h_vec(1,3),c2,
     *                    k0_copl,
     *                    nIndexes)
c-------------------------------------------------------
c Wave_copl - maximum x-ray wavelength in Angstrom:
        k0_copl_mod = DVecMod2 (k0_copl,nIndexes)
        Wave_copl   = 2.*pi / k0_copl_mod

        if (Wave .gt. 0.)  Then  !-----------------(wave>0)-------+
                                                                 !|
          if (Wave .gt. Wave_copl) goto 52 !-------------------+  |
c Length of incident wavevector in 1/Angstrom:                 v  |
          k0_mod = 2.*pi/Wave                                    !|
c+---------------+                                                |
c|Before 07/2003:|                                                |
c+=====================================================+          |
c| Non-coplanar addition to incident wavevector        |          |
c| according to the Pinsker book, p.333:               |          |
c|                    2         2    2    2         2  |          |
c| Uncopl = sqrt{ 4*k0 * [h1*h2] - h1 * h2 * (h1-h2) } |          |
c+=====================================================+          |
c                                 +-----length of h1*h2           |
c         Uncopl = (2.*k0_mod*Zone_mod)**2                       !|
c    -           - (h_mod(2)*h_mod(3)*                           !|
c    *              h_mod(nReflex1))**2                          !|
c                                                                !|
c         if (Uncopl .le. 0.) goto 53 !------------------------+  |
c         Uncopl = Dsqrt(Uncopl)                              !v  |
c         Uncopl =-Uncopl/(2.*Zone_mod**2)                       !|
c                 +--this must be either "+" or "-" !!!          !|
c+--------------+                                                 |
c|Since 07/2003:|                                                 |
c+=====================================================+          |
c| Non-coplanar addition to the incident wavevector:   |          |
c|                                                     |          |
c|               ->   ->          ->                   |          |
c|               k0 = k_copl + x*Zone                  |          |
c|                                                     |          |
c|                 2        2   2      2               |          |
c|               k0 = k_copl + x * Zone                |          |
c|                                                     |          |
c|                        2        2                   |          |
c|            x = sqrt{ k0 - k_copl } / |Zone|         |          |
c+=====================================================+          |
          Uncopl = k0_mod**2 - k0_copl_mod**2                    !|
                                                                 !|
          if (Uncopl .le. 0.) goto 53 !------------------------+  |
          Uncopl = Dsqrt(Uncopl) / Zone_mod                   !v  |
c                 +--this must be either "+" or "-" !!!          !|
                                                                 !|
c Incident wavevector at the exact Bragg condition:               |
          Call  DVecSum2 (k0_copl,1.D0,                          !|
     *                    Zone,   Uncopl,                        !|
     *                    k0,                                    !|
     *                    nIndexes)                              !|
                                                                 !|
c Added 07/2003: When the incident vector is specified via the    |
c wavelength, we have 2 choices -- above and below the coplanar   |
c plane). The normal is used to select one of the two:            |
          if (DVecMod2(Surface_Norm,nIndexes) .gt. 0.) Then !-+   |
c The angle between k0 and the                                |   |
c internal normal must be < 90 degr.:                         |   |
            ss = Dangle2(k0,Surface_Norm,nIndexes)           !|   |
            if ((ss .gt.  90.) .OR.                          !|   |
     *          (ss .lt. -90.)) Then !--------------+         |   |
              Uncopl = -Uncopl                     !|         |   |
c Incident wavevector at the exact Bragg condition: |         |   |
              Call DVecSum2 (k0_copl,1.D0,         !|         |   |
     *                       Zone,   Uncopl,       !|         |   |
     *                       k0,                   !|         |   |
     *                       nIndexes)             !|         |   |
c ccc         write (4,*) '*Inverting wavevector*' !|         |   |
             endif !--------------------------------+         |   |
           endif !--------------------------------------------+   |
                                                                 !|
          involved = 2   !number of used reflections             !|
                                                                 !|
        elseif (Abs(Wave+1.) .lt. 1.E-5) Then !----(wave=-1)------|
                                                                 !|
c Coplanar case (requested by wave=-1);                           |
          Call  DVecCopy_dd (k0_copl,                            !|
     *                       k0,                                 !|
     *                       nIndexes)                           !|
          Wave     = Wave_copl                                   !|
          k0_mod   = 2.*pi / Wave                                !|
          involved = 2   !number of used reflections             !|
                                                                 !|
        elseif (Abs(Wave+2.) .lt. 1.E-5) Then !----(wave=-2)------|
                                                                 !|
c Wave is given by the Bragg conditions for reflex-3              |
c (requested by wave=-2)                                          |
                                                                 !|
c Test that reflection-3 is present:                             !|
          if ((nReflex1-1) .lt. 4) goto 54 !-------------------+  |
                                                              !v  |
c Calculate angke between h3 and [h1xh2] -- it should not be at   |
c +-90degr; i.e. all three vectors should not lie in one plance.  |
          c1 = DVecSca2(Zone,h_vec(1,4),nIndexes)                !|
          c2 = c1 / (Zone_mod*h_mod(4))                          !|
                                                                 !|
          if (abs(c2) .lt. 1.E-4) goto 55 !--------------------+  |
                                                              !v  |
c+=====================================================+          |
c| Non-coplanar addition to the incident wavevector:   |          |
c|                                                     |          |
c|         -> -> 2  -> 2        -> ->     2            |          |
c|        (k0+h3) = k0    =>  2*k0*h3 + h3 = 0         |          |
c|                                                     |          |
c|               ->   ->          ->                   |          |
c|               k0 = k_copl + x*Zone                  |          |
c|                                                     |          |
c|         ->     ->         ->  ->      2             |          |
c|      2*(k_copl*h3) + 2x*(Zone*h3) + h3 = 0          |          |
c|                                                     |          |
c|                2  ->     ->    ->  ->               |          |
c|     x =-(0.5*h3 + k_copl*h3)/(Zone*h3)              |          |
c|                                                     |          |
c+=====================================================+          |
          c2 = DVecSca2(k0_copl,h_vec(1,4),nIndexes)             !|
          ss =-(0.5*h_mod(4)**2 + c2) / c1                       !|
                                                                 !|
c Incident wavevector at the exact Bragg condition:               |
          Call  DVecSum2 (k0_copl,1.D0,                          !|
     *                    Zone,   ss,                            !|
     *                    k0,                                    !|
     *                    nIndexes)                              !|
          k0_mod   = DVecMod2 (k0,nIndexes)                      !|
          Wave     = 2.*pi / k0_mod                              !|
          involved = 3   !number of used reflections             !|
                                                                 !|
        else !----------------------------------------------------|
                                                                 !|
c Incorrect wavelength:                                           |
          goto 56 !--------------------------------------------+  |
                                                              !v  |
        endif  !--------------------------------------------------+

c Test: is the Bragg condition indeed satisfied?
        do i=2,involved+1 !==============================+
          c1 = 2.*DVecSca2 (k0,h_vec(1,i),nIndexes)     !|
          c2 = h_mod(i)**2                              !|
          ss = abs(c1+c2)/k0_mod**2                     !|
          if (ss .gt. (1.E-5)) goto 57 !-----------------+-----+
        enddo  !=========================================+     v

c Test: is the angle with internal surface normal
c less than 90 degr.?
        if (DVecMod2(Surface_Norm,nIndexes) .gt. 0.) Then !-+
c The angle between k0 and the                              |
c internal normal must be < 90 degr.:                       |
          ss = Dangle2(k0,Surface_Norm,nIndexes)           !|
          if ((ss .gt.  90.) .OR.                          !|
     *        (ss .lt. -90.)) goto 58 !---------------------+--+
        endif !---------------------------------------------+  v

        ifail = 0
  28    continue
        if (iirezv .eq. 1)  Then  !--------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c================================================================
c                          E R R O R S
c----------------------------------------------------------------
  51    continue
        write (txt,100) progname(1:lpn), j-1
  100   format(1x,a,' ERROR:'//
     *  ' Reflex <',i2,'> - zero indexes')
        ifail = 3
        goto 999

  52    continue
        write (txt,1) progname(1:lpn), Wave, Wave_copl
  1     format(1x,a,' ERROR:'//
     *  ' The incident wavelength=',g14.7,' Angstrem is greater'/
     *  ' than the maximum  wavelength=',g14.7,' Angstrem'/
     *  ' for the specified set of diffracting planes.'//
     *  ' The Bragg condition(s) cannot be satisfied.')
        ifail = 7
        goto 999

  53    continue
        write (txt,10) progname(1:lpn)
  10    format(1x,a,' ERROR:'//
     *  'Noncoplanar part of k0 is negative')
        ifail = 3
        goto 999

  54    continue
        write (txt,2) progname(1:lpn), nReflex1-1
  2     format(1x,a,' ERROR:'//
     *  ' Insufficient number of Bragg reflections'/
     *  ' "',i1,'" to determine the x-ray wavelength.'/
     *  ' The minimum number should be "4"')
        ifail = 5
        goto 999

  55    continue
        write (txt,3) progname(1:lpn)
  3     format(1x,a,' ERROR:'//
     *  ' Cannot build the incident x-ray wavelength'/
     *  ' since all 3 Bragg vectors lie in one plane.'/
     *  ' Please specify the wavelength or select'/
     *  ' "coplanar case" diffraction geometry')
        ifail = 6
        goto 999

  56    continue
        write (txt,6) progname(1:lpn), wave
  6     format(1x,a,' ERROR:'//
     *  ' Incorrect x-ray wavelength specification'//
     *  '      wavelength = ', g14.7,' Angstrem')
        ifail = 5
        goto 999

  57    continue
        write (txt,7) progname(1:lpn), ss, i,
     *                ind_h(1,i),ind_h(2,i), ind_h(nIndexes,i)
  7     format(1x,a,' ERROR:'//
     *  ' Cross-checking shows that the built x-ray wavevector'/
     *  ' deviates by |alpha|=',g10.3,' > 1.E-5 from reflection #',i1/
     *  ' with Miller indices (',3i4,').'//
     *  ' Please report to: sstepanov@anl.gov')
        ifail = 7
        goto 999

  58    continue
        c1 = 0.
        do i=1,nIndexes !==========================+
          if (abs(k0(i)) .gt. c1) c1 = abs(k0(i)) !|
        enddo !====================================+
        do i=1,nIndexes  !===============+
          iii(i) = int(100.*k0(i)/c1)   !|
        enddo !==========================+
        Call PakInt (iii,nIndexes,wrk,i)
        if (i.lt.1) i=1
        write (txt,8) progname(1:lpn), ss, wrk(1:i)
  8     format(1x,a,' ERROR:'//
     *  ' The x-ray wavevector built with given conditions makes an'/
     *  ' obtuse angle of ',f7.3,' degr. with internal surface normal.'/
     *  ' The x-rays with such a wavevector cannot enter the crystal.'//
     *  ' Info: the wavevector is approximately directed along:'/
     *  ' (',a,')')
        ifail = 8
        goto 999

  59    continue
        write (txt,9) progname(1:lpn),
     *                ind_h(1,2),ind_h(2,2), ind_h(nIndexes,2),
     *                ind_h(1,3),ind_h(2,3), ind_h(nIndexes,3)
  9     format(1x,a,' ERROR:'//
     *  ' The reciprocal lattice vectors h1=(',3i4,') and h2=(',3i4,')'/
     *  ' seem to be parallel'//
     *  ' The diffractrion geometry cannot be built')
        ifail = 6
        goto 999

  999   continue
c       Call  Message (txt,ifail,2)
        goto 28
        end

c ==================================================================

        Subroutine Build_h (indices,h_vector,h_module,nIndexes)

        Integer nIndexes,
     *          indices(nIndexes),
     *          i
        Real*8  h_vector(nIndexes),
     *          h_module, pi

        Real*8   DVecMod
        External DVecMod

        pi = 4.*Datan(1.D0)

c Since |indices|=1./d (d is Bragg planes spacing, see Bragan),
c then h = 2*pi*indices:

        h_module = 2.*pi*DVecMod(indices,nIndexes)

        do    i=1,nIndexes  !==============+
          h_vector(i) = 2.*pi*indices(i)  !|
        enddo  !===========================+

        return
        end
