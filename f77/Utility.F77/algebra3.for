c
c
c          **   Vector operations in crystals **
c          **     of an arbitrary symmetry.   **
c          **      Real*8 indices version     **
c                  ======================
c    Contents:
c    --------
c    -1. Subroutine      DlatReset   ()
c     0. Subroutine      DBackLat    ()
c    1a. Real*8 Function DAngle      (iv1,iv2,k)
c    1b. Real*8 Function DAngle2     (rv1,rv2,k)
c    2a. Subroutine      DProject    (iz,iv,pro,k)
c    2b. Subroutine      DProject2   (rz,rv,pro,k)
c    3a.  ------
c    3b. Subroutine      DUnitVec2   (v,u,k)
c    4a.  ------
c    4b. Subroutine      DVecCon     (v,c,u,k)
c    5a. Real*8 Function DVecMod     (iv,k)
c    5b. Real*8 Function DVecMod2    (rv,k)
c    6a. Real*8 Function DVecSca     (iv1,iv2,k)
c    6b. Real*8 Function DVecSca2    (rv1,rv2,k)
c    7a.  ------
c    7b. Subroutine      DVecSum2    (v1,c1,v2,c2,v3,k)
c    8a.  ------
c    8b. Subroutine      DVecVec     (rv1,rv2,rvv,k)
c    9a. Subroutine      DVecCopy_id (iv,u,k)
c    9b. Subroutine      DVecCopy_dd (v,u,k)
c
c ("a": integer indices, "b": real*8 indices)
c
c            +---------------------------------------+
c            |  ATTENTION:!                          |
c            |  When compiling for Microsoft Fortran |
c            |  library, there is a problem with     |
c            |  combining the /FPc and /Zl keys.     |
c            |  The /FPc key should be dropped.      |
c            +---------------------------------------+

c=====================================================
        Subroutine      DlatReset ()
c-----------------------------------------------------
c Reset crystal lattice (double-precision version)
c-----------------------------------------------------
        Integer         i

        Real            a1(6), vol1
        Integer         lbev1
        Common /back1/  a1, vol1, lbev1

        Real*8          at(6), g(3,3), gb(3,3), vol
        Integer         lbev
        Common /Dback2/ at, g, gb, vol, lbev

        do      i=1,6  !===+
           a1(i) = 0.     !|
           at(i) = 0.     !|
        enddo  !===========+
        vol   = 0.
        vol1  = 0.
        lbev  = 0               !parameters are reset
        lbev1 = 0               !parameters are reset
        return
        end

c=====================================================
        Subroutine      DBackLat ()
c-----------------------------------------------------
c Determine crystal reciprocal lattice and the components 
c of metric tensor
c-----------------------------------------------------
        Integer         i
        Real*8          a(6), gra, wol,
     *                  sialf, sibet, sigam,
     *                  coalf, cobet, cogam,
     *                  coalft, cobett, cogamt

        Real            a1(6), vol1
        Integer         lbev1
        Common /back1/  a1, vol1, lbev1

        Real*8          at(6), g(3,3), gb(3,3), vol
        Integer         lbev
        Common /Dback2/ at, g, gb, vol, lbev

        gra = Datan(1.D+0)/45.          !number of radians in 1 degree
        if (lbev1.eq.0) then !--+
           Call Backlat ()     !|see algebra.for
        endif  !----------------+
        do      i=1,6  !===+
           a(i) = a1(i)   !|
        enddo  !===========+
c Cosines and sinuses of direct lattice unit cell angles:
        coalf = Dcos(a(4)*gra)
        cobet = Dcos(a(5)*gra)
        cogam = Dcos(a(6)*gra)
        sialf = Dsin(a(4)*gra)
        sibet = Dsin(a(5)*gra)
        sigam = Dsin(a(6)*gra)
        if (sialf.lt.1.D-20 .OR. 
     *      sibet.lt.1.D-20 .OR. 
     *      sigam.lt.1D-20) then !-----------------------+
           write (0,*) ' DBackLat:  a = (',a,')'        !|
           write (0,*) ' DBackLat:  angle(s) <= 0'      !|
           stop  'DBackLat: angle(s) <= 0'              !|
        endif  !-----------------------------------------+
c Elementary cell volume (Sirotin, Shaskolskaya, "Basics of
c Crystal Physics") - "Osnovi kristallophisiki", p.96:
        vol = a(1)*a(2)*a(3)
        wol = 1. - coalf*coalf - cobet*cobet - cogam*cogam
     +      + 2.*coalf*cobet*cogam
        if     (vol.lt.1.D-37)  then  !------------------+
           write (0,*) ' DBackLat:  a = (',a,')'        !|
           write (0,*) ' DBackLat:  vol = ',vol,' <= 0' !|
           stop  'DBackLat: vol <= 0'                   !|
        elseif (wol.lt.1.D-37)  then  !------------------+
           write (0,*) ' DBackLat:  a = (',a,')'        !|
           write (0,*) ' DBackLat:  wol = ',wol,' <= 0' !|
           stop  'DBackLat: wol <= 0'                   !|
        endif  !-----------------------------------------+
        wol = Dsqrt(wol)
        vol = vol*wol
c Reciprocal lattice unit cell parameters:
        at(1) = sialf/(a(1)*wol)
        at(2) = sibet/(a(2)*wol)
        at(3) = sigam/(a(3)*wol)
c Cosines of reciprocal lattice angles:
        coalft = (cobet*cogam-coalf) / (sibet*sigam)
        cobett = (cogam*coalf-cobet) / (sigam*sialf)
        cogamt = (coalf*cobet-cogam) / (sialf*sibet)
c Reciprocal lattice unit cell angles:
        at(4) = Datan2(Dsqrt(1.-coalft*coalft),coalft)/gra
        at(5) = Datan2(Dsqrt(1.-cobett*cobett),cobett)/gra
        at(6) = Datan2(Dsqrt(1.-cogamt*cogamt),cogamt)/gra
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

c=====================================================
        Real*8 Function  DAngle (iv1,iv2,k)
c-----------------------------------------------------
c Determine the angle between two planes given by their
c normals or between two vectors in a crystal with
c arbitrary symmetry.
c The integer indices version.
c-----------------------------------------------------
        Integer         k, iv1(k), iv2(k)
        Real*8          s, d1, d2, dd, gra

        Real*8          DVecMod, DVecSca
        External        DVecMod, DVecSca

        d1 = DVecMod(iv1,k)             !Length of vector V1
        d2 = DVecMod(iv2,k)             !Length of vector V2
        dd = d1*d2

        if (dd.lt.1.D-30)  then  !------------------------------+
           write (0,*) ' DAngle:  V1 = (',iv1,'), |V1| = ',d1  !|
           write (0,*) ' DAngle:  V2 = (',iv2,'), |V1| = ',d2  !|
           write (0,*) ' DAngle:  |V1|*|V2| = ',dd             !|
           stop 'DAngle: |V1|*|V2|=0 or < 1E-30'               !|
        endif  !------------------------------------------------+

        s   = DVecSca(iv1,iv2,k) / dd   !Cosine of angle between V1 & V2
        gra = Datan(1.D+0)/45.          !number of radians in 1 degree

        if (s.gt. 0.999) Then !-------------------------+
                                                       !|
           if (s.ge.1.D0) Then !--------------+         |
              DAngle = 0.                     !|        |
           else  !----------------------------|         |
              DAngle = Dsqrt(2.*(1.-s)) / gra !|        |
           endif  !---------------------------+         |
                                                       !|
        elseif (s.lt.-0.999) Then !---------------------+
                                                       !|
           if (s.le.-1.D0) Then !---------------------+ |
              DAngle = 180.                          !| |
           else  !------------------------------------+ |
              DAngle = 180. - Dsqrt(2.*(1.+s)) / gra !| |
           endif  !-----------------------------------+ |
                                                       !|
        else  !-----------------------------------------+
                                                       !|
           DAngle = DAtan2(Dsqrt(1.-s*s),s) / gra      !|
                                                       !|
        endif  !----------------------------------------+
        return
        end

c=====================================================
        Real*8 Function  DAngle2 (rv1,rv2,k)
c-----------------------------------------------------
c Determine the angle between two planes given by their
c normals or between two vectors in a crystal with
c arbitrary symmetry. 
c The real*8 indices version.
c-----------------------------------------------------
        Integer         k
        Real*8          rv1(k), rv2(k), s, d1, d2, dd, gra

        Real*8          DVecMod2, DVecSca2
        External        DVecMod2, DVecSca2

        d1 = DVecMod2(rv1,k)            !Length of vector V1
        d2 = DVecMod2(rv2,k)            !Length of vector V2
        dd = d1*d2

        if (dd.lt.1.D-30)  then  !-------------------------------+
           write (0,*) ' DAngle2:  V1 = (',rv1,'), |V1| = ',d1  !|
           write (0,*) ' DAngle2:  V2 = (',rv2,'), |V1| = ',d2  !|
           write (0,*) ' DAngle2:  |V1|*|V2| = ',dd             !|
           stop 'DAngle2: |V1|*|V2|=0 or < 1E-30'               !|
        endif  !-------------------------------------------------+

        s   = DVecSca2(rv1,rv2,k) / dd  !Cosine of angle between V1 & V2
        gra = Datan(1.D+0)/45.          !number of radians in 1 degree

        if (s.gt. 0.999) Then !--------------------------+
                                                        !|
           if (s.ge.1.D0) Then !----------------+        |
              DAngle2 = 0.                     !|        |
           else  !------------------------------+        |
              DAngle2 = Dsqrt(2.*(1.-s)) / gra !|        |
           endif  !-----------------------------+        |
                                                        !|
        elseif (s.lt.-0.999) Then !----------------------+
                                                        !|
           if (s.le.-1.D0) Then !----------------------+ |
              DAngle2 = 180.                          !| |
           else  !-------------------------------------+ |
              DAngle2 = 180. - Dsqrt(2.*(1.+s)) / gra !| |
           endif  !------------------------------------+ |
                                                        !|
        else !-------------------------------------------+
                                                        !|
           DAngle2 = DAtan2(Dsqrt(1.-s*s),s) / gra      !|
                                                        !|
        endif  !-----------------------------------------+
        return
        end

c=====================================================
        Subroutine      DProject (iz,iv,pro,k)
c-----------------------------------------------------
c Determine projection of vector V (indices iv)
c on a plane normal to vector Z (indices iz) 
c in a crystal with arbitrary symmetry.
c
c  -->   ->   ->  ->   ->   2
c  PRO = V - (V * Z) * Z / Z
c
c The integer indices version.
c-----------------------------------------------------
        Integer         k, iz(k), iv(k)
        Real*8          pro(k), z, v

        Real*8          DVecMod, DVecSca
        External        DVecMod, DVecSca

        z = DVecMod (iz,k)              !Length of vector Z
        if (z.lt.1.D-20)  then  !------------------------------+
           write (0,*) ' DProject:  Z = (',iz,'),  |Z| = ',z  !|
           stop 'DProject: Z=0 or < 1E-20'                    !|
        endif  !-----------------------------------------------+
        v = DVecSca (iv,iz,k) / (z*z)   !(V*Z)/Z*Z
        pro(1) = iv(1)-iz(1)*v          !projection-1
        pro(2) = iv(2)-iz(2)*v          !projection-2
        pro(k) = iv(k)-iz(k)*v          !projection-3
        if (k.eq.4) pro(3)=-(pro(1)+pro(2))
        return
        end

c=====================================================
        Subroutine      DProject2 (rz,rv,pro,k)
c-----------------------------------------------------
c Determine projection of vector V (indices rv)
c on a plane normal to vector Z (indices rz) 
c in a crystal with arbitrary symmetry.
c
c  -->   ->   ->  ->   ->   2
c  PRO = V - (V * Z) * Z / Z
c
c The real*8 indices version.
c-----------------------------------------------------
        Integer         k
        Real*8          rz(k),rv(k), pro(k), z, v

        Real*8          DVecMod2, DVecSca2
        External        DVecMod2, DVecSca2

        z = DVecMod2 (rz,k)             !Length of vector Z
        if (z.lt.1.D-20)  then  !-------------------------------+
           write (0,*) ' DProject2:  Z = (',rz,'),  |Z| = ',z  !|
           stop 'DProject2: Z=0 or < 1E-20'                    !|
        endif  !------------------------------------------------+
        v = DVecSca2 (rv,rz,k) / (z*z)  !(V*Z)/Z*Z
        pro(1) = rv(1)-rz(1)*v          !projection-1
        pro(2) = rv(2)-rz(2)*v          !projection-2
        pro(k) = rv(k)-rz(k)*v          !projection-3
        if (k.eq.4) pro(3)=-(pro(1)+pro(2))
        return
        end

c=====================================================
        Subroutine      DUnitVec2 (v,u,k)
c-----------------------------------------------------
c Determine unit-length vector along given vector
c-----------------------------------------------------
        Integer         k, i
        Real*8          v(k), u(k), vm

        Real*8          DVecMod2
        External        DVecMod2

        vm = DVecMod2 (v,k)
        if (Abs(vm).lt.1.D-37)   return
        do      i=1,k  !=======+
           u(i) = v(i) / vm   !|
        enddo  !===============+
        return
        end

c=====================================================
        Subroutine      DVecCon (v,c,u,k)
c-----------------------------------------------------
c Multiply a vector by a constant (change vector length)
c-----------------------------------------------------
        Integer         k, i
        Real*8          v(k),u(k),c

        do      i=1,k  !=====+
           u(i) = c*v(i)    !|
        enddo  !=============+
        return
        end

c=====================================================
        Real*8 Function DVecMod (iv,k)
c-----------------------------------------------------
c Calculate the length of vector V (indices iv) in a
c crystal with arbitrary symmetry. 
c The integer indices version.
c-----------------------------------------------------
        Integer         k, iv(k)
        Real*8          v2

        Real*8          DvecSca
        External        DvecSca

        v2     = DVecSca (iv,iv,k)      !square of |V|
        DVecMod= Dsqrt(v2)              ! |V|
        return
        end

c=====================================================
        Real*8 Function DVecMod2 (rv,k)
c-----------------------------------------------------
c Calculate the length of vector V (indices rv) in a
c crystal with arbitrary symmetry. 
c The real*8 indices version.
c-----------------------------------------------------
        Integer         k
        Real*8          rv(k), v2

        Real*8          DVecSca2
        External        DVecSca2

        v2      = DVecSca2 (rv,rv,k)    !square of |V|
        DVecMod2 = Dsqrt(v2)            ! |V|
        return
        end

c=====================================================
        Real*8 Function DVecSca (iv1,iv2,k)
c-----------------------------------------------------
c Calculate scalar product of vector V1 (indices iv1)
c and V2 (indices iv2) in a crystal with arbitrary 
c symmetry.

c (V1*V2) = Sum ( iv1(i) * iv2(j) * gb(i,j) )
c           i,j
c
c The integer indices version.
c-----------------------------------------------------
        Integer         k, iv1(k), iv2(k)

        Real            a1(6), vol1
        Integer         lbev1
        Common /back1/  a1, vol1, lbev1

        Real*8          at(6), g(3,3), gb(3,3), vol
        Integer         lbev
        Common /Dback2/ at,g,gb,vol,lbev

c Determine reciprocal lattice:
        if (lbev.eq.0 .or. lbev1.eq.0) then !--+
          Call DBacklat ()                    !|
        endif  !-------------------------------+
c Calculate scalar product:
        DVecSca = iv1(1)*iv2(1)*gb(1,1)+
     +            iv1(2)*iv2(2)*gb(2,2)+
     +            iv1(k)*iv2(k)*gb(3,3)+
     +           (iv1(1)*iv2(2)+iv1(2)*iv2(1))*gb(1,2)+
     +           (iv1(1)*iv2(k)+iv1(k)*iv2(1))*gb(1,3)+
     +           (iv1(2)*iv2(k)+iv1(k)*iv2(2))*gb(2,3)
        return
        end

c=====================================================
        Real*8 Function DVecSca2 (rv1,rv2,k)
c-----------------------------------------------------
c Calculate scalar product of vector V1 (indices rv1)
c and V2 (indices rv2) in a crystal with arbitrary 
c symmetry.

c (V1*V2) = Sum ( rv1(i) * rv2(j) * gb(i,j) )
c           i,j
c
c The real*8 indices version.
c-----------------------------------------------------
        Integer         k
        Real*8          rv1(k), rv2(k)

        Real            a1(6), vol1
        Integer         lbev1
        Common /back1/  a1, vol1, lbev1

        Real*8          at(6), g(3,3), gb(3,3), vol
        Integer         lbev
        Common /Dback2/ at, g, gb, vol, lbev

c Determine reciprocal lattice:
        if (lbev.eq.0 .or. lbev1.eq.0) then !--+
           Call DBacklat ()                   !|
        endif  !-------------------------------+
c Calculate scalar product:
        DVecSca2 = rv1(1)*rv2(1)*gb(1,1)+
     +             rv1(2)*rv2(2)*gb(2,2)+
     +             rv1(k)*rv2(k)*gb(3,3)+
     +            (rv1(1)*rv2(2)+rv1(2)*rv2(1))*gb(1,2)+
     +            (rv1(1)*rv2(k)+rv1(k)*rv2(1))*gb(1,3)+
     +            (rv1(2)*rv2(k)+rv1(k)*rv2(2))*gb(2,3)
        return
        end

c=====================================================
        Subroutine      DVecSum2 (v1,c1,v2,c2,v3,k)
c-----------------------------------------------------
c Determine a sum of teo vectors v1 and v2 in a 
c crystal with arbitrary symmetry
c
c     ->      ->      ->
c     v3 = c1*v1 + c2*v2
c-----------------------------------------------------
        Integer         k
        Real*8          v1(k), v2(k), v3(k), c1, c2

        v3(1) = c1*v1(1) + c2*v2(1)
        v3(2) = c1*v1(2) + c2*v2(2)
        v3(k) = c1*v1(k) + c2*v2(k)
        if (k.eq.4) v3(3) = -(v3(1)+v3(2))
        return
        end

c=====================================================
        Subroutine DVecVec (rv1,rv2,rvv,k)
c-----------------------------------------------------
c Calculate vector product of vectors V1 (indices rv1)
c and V2 (indices rv2) in a crystal with arbitrary 
c symmetry.
c
c             1
c [V1*V2]  = --- Summ (rv1(i)*rv2(j)*delt(i,j,k)*g(k,l)),
c       (l)  Vol i,j,k
c
c here delt(i,j,k) is a Levi-Civita symbol.
c-----------------------------------------------------
        Integer         k, i
        Real*8          rv1(k), rv2(k), rvv(k), q(3)

        Real            a1(6), vol1
        Integer         lbev1
        Common /back1/  a1, vol1, lbev1

        Real*8          at(6), g(3,3), gb(3,3), vol
        Integer         lbev
        Common /Dback2/ at, g, gb, vol, lbev

c Determine reciprocal lattice:
        if (lbev.eq.0 .or. lbev1.eq.0) then !--+
           Call DBacklat ()                   !|
        endif  !-------------------------------+
        if (vol.lt.1.D-37) then !-----------------------+
           write (0,*) ' DVecVec:  a = (',at,')'       !|
           write (0,*) ' DVecVec:  vol = ',vol,' <= 0' !|
           stop  'DVecVec: vol <= 0'                   !|
        endif !-----------------------------------------+
c Vector product:
        q(1) = (rv1(2)*rv2(k)-rv1(k)*rv2(2)) / vol
        q(2) = (rv1(k)*rv2(1)-rv1(1)*rv2(k)) / vol
        q(3) = (rv1(1)*rv2(2)-rv1(2)*rv2(1)) / vol
        do      i=1,3  !============+
           rvv(i) = q(1)*g(1,i)    !|
     +            + q(2)*g(2,i)    !|
     +            + q(3)*g(3,i)    !|
        enddo   !===================+
        if (k.eq.4)     then   !--------+ if 4 indices like in hexagonal crystals
           rvv(4) = rvv(3)             !|
           rvv(3) =-(rvv(1)+rvv(2))    !|
        endif  !------------------------+
        return
        end

c=====================================================
        Subroutine      DVecCopy_id (iv,u,k)
c-----------------------------------------------------
c Copy vector from integer indices to real*8
c-----------------------------------------------------
        Integer         k, iv(k), i
        Real*8          u(k)

        do      i=1,k  !=====+
           u(i) = iv(i)     !|
        enddo  !=============+
        return
        end

c=====================================================
        Subroutine      DVecCopy_dd (v,u,k)
c-----------------------------------------------------
c Copy vector from real*8 to real*8 indices
c-----------------------------------------------------
        Integer         k, i
        Real*8          v(k), u(k)

        do      i=1,k  !=====+
           u(i) = v(i)      !|
        enddo  !=============+
        return
        end
