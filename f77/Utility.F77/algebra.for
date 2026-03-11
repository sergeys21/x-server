c
c
c          **  Vector operations in crystals  **
c          **    of an arbitrary symmetry.    **
c          **    Integer indices version      **
c                =======================
c    Content:
c    --------
c   0. Subroutine    LatReset   ()
c   1. Subroutine    BackLat    ()
c   2. Real Function VecSca     (iv1,iv2,k)
c   3. Real Function VecMod     (iv,k)
c   4. Real Function Angle      (iv1,iv2,k)
c   5. Subroutine    VecCopy_ii (iv,iu,k)
c            +---------------------------------------+
c            |  ATTENTION:!                          |
c            |  When compiling for Microsoft Fortran |
c            |  library, there is a problem with     |
c            |  combining the /FPc and /Zl keys.     |
c            |  The /FPc key should be dropped.      |
c            +---------------------------------------+
c
c=====================================================
c
        subroutine      LatReset   ()
c-----------------------------------------------------
c Reset crystal lattice
c-----------------------------------------------------
        Integer         i

        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        Real            at(6), g(3,3), gb(3,3)
        Common  /back2/ at, g, gb

        do      i=1,6  !===+
           a(i)  = 0.     !|
           at(i) = 0.     !|
        enddo  !===========+
        vol  = 0.
        lbev = 0               !parameters are reset
        return
        end
c
c=====================================================
c
        subroutine      BackLat ()
c-----------------------------------------------------
c Determine crystal reciprocal lattice and the components 
c of metric tensor
c-----------------------------------------------------
        Real            gra, wol,
     *                  coalf,  cobet,  cogam,
     *                  sialf,  sibet,  sigam,
     *                  coalft, cobett, cogamt

        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        Real            at(6), g(3,3), gb(3,3)
        Common  /back2/ at, g, gb

        gra = Atan(1.)/45.              !number of radians in 1 degree
c Cosines and sinuses of direct lattice unit cell angles:
        coalf = cos(a(4)*gra)
        cobet = cos(a(5)*gra)
        cogam = cos(a(6)*gra)
        sialf = sin(a(4)*gra)
        sibet = sin(a(5)*gra)
        sigam = sin(a(6)*gra)
        if (sialf.lt.1.E-20 .OR. 
     *      sibet.lt.1.E-20 .OR. 
     *      sigam.lt.1E-20) then !----------------------+
           write (0,*) ' BackLat:  a = (',a,')'        !|
           write (0,*) ' BackLat:  angle(s) <= 0'      !|
           stop  'BackLat: angle(s) <= 0'              !|
        endif  !----------------------------------------+
c Elementary cell volume (Sirotin, Shaskolskaya, "Basics of
c Crystal Physics") - "Osnovi kristallophisiki", p.96:
        vol = a(1)*a(2)*a(3)
        wol = 1. - coalf*coalf - cobet*cobet - cogam*cogam
     +      + 2.*coalf*cobet*cogam
        if     (vol.lt.1.E-37)  then   !-----------------+
           write (0,*) ' BackLat:  a = (',a,')'         !|
           write (0,*) ' BackLat:  vol = ',vol,' <= 0'  !|
           stop  'Backlat: vol <= 0'                    !|
        elseif (wol.lt.1.E-37)  then   !-----------------+
           write (0,*) ' BackLat:  a = (',a,')'         !|
           write (0,*) ' BackLat:  wol = ',wol,' <= 0'  !|
           stop  'Backlat: wol <= 0'                    !|
        endif  !-----------------------------------------+
        wol = sqrt(wol)
        vol = vol*wol
c Reciprocal lattice unit cell parameters:
        at(1) = sialf/(a(1)*wol)
        at(2) = sibet/(a(2)*wol)
        at(3) = sigam/(a(3)*wol)
c Cosines of reciprocal lattice angles:
        coalft = (cobet*cogam-coalf) / (sibet*sigam)
        cobett = (cogam*coalf-cobet) / (sigam*sialf)
        cogamt = (coalf*cobet-cogam) / (sialf*sibet)
c Reciprocal lattice angles unit cell angles:
        at(4) = atan2(sqrt(1.-coalft*coalft),coalft)/gra
        at(5) = atan2(sqrt(1.-cobett*cobett),cobett)/gra
        at(6) = atan2(sqrt(1.-cogamt*cogamt),cogamt)/gra
c Components of direct lattice metric tensor (Sirotin, Shaskolskaya, 
c "Basics of Crystal Physics") - "Osnovi kristallophisiki", p.625:
        g(1,1) = a(1)*a(1)
        g(2,1) = a(2)*a(1)*cogam
        g(3,1) = a(3)*a(1)*cobet
        g(1,2) = g(2,1)
        g(2,2) = a(2)*a(2)
        g(3,2) = a(3)*a(2)*coalf
        g(1,3) = g(3,1)
        g(2,3) = g(3,2)
        g(3,3) = a(3)*a(3)
c Components of reciprocal lattice metric tensor (Sirotin, Shaskolskaya,
c "Basics of Crystal Physics") - "Osnovi kristallophisiki", p.625:
        gb(1,1) = at(1)*at(1)
        gb(2,1) = at(2)*at(1)*cogamt
        gb(3,1) = at(3)*at(1)*cobett
        gb(1,2) = gb(2,1)
        gb(2,2) = at(2)*at(2)
        gb(3,2) = at(3)*at(2)*coalft
        gb(1,3) = gb(3,1)
        gb(2,3) = gb(3,2)
        gb(3,3) = at(3)*at(3)
        lbev = 1                !flag indicating that the parameters are calculated
        return
        end
c
c=====================================================
c
        Real Function VecSca (iv1,iv2,k)
c-----------------------------------------------------
c Calculate scalar product of vector V1 (indices iv1)
c and V2 (indices iv2) in a crystal with arbitrary 
c symmetry.

c (V1*V2) = Sum ( iv1(i) * iv2(j) * gb(i,j) )
c           i,j
c-----------------------------------------------------
        Integer         k, iv1(k), iv2(k)

        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        Real            at(6), g(3,3), gb(3,3)
        Common  /back2/ at, g, gb

c Determine reciprocal lattice:
        if (lbev.eq.0)  then  !--------+
           Call Backlat ()            !|
        endif  !-----------------------+
c Calculate scalar product:
        VecSca = iv1(1)*iv2(1)*gb(1,1)+
     +           iv1(2)*iv2(2)*gb(2,2)+
     +           iv1(k)*iv2(k)*gb(3,3)+
     +          (iv1(1)*iv2(2)+iv1(2)*iv2(1))*gb(1,2)+
     +          (iv1(1)*iv2(k)+iv1(k)*iv2(1))*gb(1,3)+
     +          (iv1(2)*iv2(k)+iv1(k)*iv2(2))*gb(2,3)
        return
        end
c
c=====================================================
c
        Real Function VecMod (iv,k)
c-----------------------------------------------------
c Calculate the length of vector V (indices iv) in a
c crystal with arbitrary symmetry. 
c-----------------------------------------------------
        Integer         k, iv(k)
        Real            v2

        Real            VecSca
        External        VecSca

        v2     = VecSca (iv,iv,k)       !square of |V|
        VecMod = sqrt(v2)               ! |V|
        return
        end
c
c=====================================================
c
        Real Function Angle (iv1,iv2,k)
c-----------------------------------------------------
c Determine the angle between two planes given by their
c normals or between two vectors in a crystal with
c arbitrary symmetry.
c-----------------------------------------------------
        Integer         k, iv1(k), iv2(k)
        Real            d1, d2, dd, s, gra

        Real            VecMod, VecSca
        External        VecMod, VecSca

        d1 = VecMod(iv1,k)              !Length of vector V1
        d2 = VecMod(iv2,k)              !Length of vector V2
        dd = d1*d2

        if(dd.lt.1.e-30)  then  !------------------------------+
           write (0,*) ' Angle:  V1 = (',iv1,'), |V1| = ',d1  !|
           write (0,*) ' Angle:  V2 = (',iv2,'), |V1| = ',d2  !|
           write (0,*) ' Angle:  |V1|*|V2| = ',dd             !|
           stop 'Angle: |V1|*|V2|=0 or < 1E-30'               !|
        endif  !-----------------------------------------------+

        s   = VecSca(iv1,iv2,k) / dd    !Cosine of angle between V1 & V2
        gra = Atan(1.) / 45.            !number of radians in 1 degree

        if (s.gt. 0.999) Then  !------------------------+
                                                       !|
           if (s.ge.1.) Then !---------------+          |
              Angle = 0.                    !|          |
           else  !---------------------------+          |
              Angle = Sqrt(2.*(1.-s)) / gra !|          |
           endif  !--------------------------+          |
                                                       !|
        elseif (s.lt.-0.999) Then  !--------------------+
                                                       !|
           if (s.le.-1.) Then !---------------------+   |
              Angle = 180.                         !|   |
           else  !----------------------------------+   |
              Angle = 180. - Sqrt(2.*(1.+s)) / gra !|   |
           endif  !---------------------------------+   |
                                                       !|
        else !------------------------------------------+
                                                       !|
           Angle = Atan2(sqrt(1.-s*s),s) / gra         !|
                                                       !|
        endif  !----------------------------------------+
        return
        end

c=====================================================

        Subroutine      VecCopy_ii (iv,iu,k)
c-----------------------------------------------------
c Copy vector from integer indices to integer
c-----------------------------------------------------
        Integer         k, iv(k), iu(k), i

        do      i=1,k  !=====+
           iu(i) = iv(i)    !|
        enddo  !=============+
        return
        end

