c
c
c          **   Vector operations in crystals **
c          **     of an arbitrary symmetry.   **
c          **      Real*4 indices version     **
c                  ======================
c    Contents:
c    --------
c    1b. Real Function VecSca2    (rv1,rv2,k)
c    2b. Real Function VecMod2    (rv,k)
c    3b. Real Function Angle2     (rv1,rv2,k)
c    4b. Subroutine    VecVec     (rv1,rv2,rvv,k)
c    5a. Subroutine    Project    (iz,iv,pro,k)
c    5b. Subroutine    Project2   (rz,rv,pro,k)
c    6b. Subroutine    VecSum2    (v1,c1,v2,c2,v3,k)
c    7b. Subroutine    UnitVec2   (v,u,k)
c    8b. Subroutine    VecCon     (v,c,u,k)
c    9a. Subroutine    VecCopy_ir (iv,u,k)
c    9b. Subroutine    VecCopy_rr (v,u,k)
c
c ("a": integer indices, "b": real*4 indices)
c
c            +---------------------------------------+
c            |  ATTENTION:!                          |
c            |  When compiling for Microsoft Fortran |
c            |  library, there is a problem with     |
c            |  combining the /FPc and /Zl keys.     |
c            |  The /FPc key should be dropped.      |
c            +---------------------------------------+
 
c=====================================================
        Real Function VecSca2 (rv1,rv2,k)
c-----------------------------------------------------
c Calculate scalar product of vector V1 (indices rv1)
c and V2 (indices rv2) in a crystal with arbitrary 
c symmetry.

c (V1*V2) = Sum ( rv1(i) * rv2(j) * gb(i,j) )
c           i,j
c
c The real*4 indices version.
c-----------------------------------------------------
        Integer         k
        Real            rv1(k), rv2(k)

        Real    a(6), vol
        Integer lbev
        Common  /back1/ a, vol, lbev

        Real    at(6), g(3,3), gb(3,3)
        Common  /back2/ at, g, gb

c Determine reciprocal lattice:
        if (lbev.eq.0)  then  !-----+
           Call Backlat ()         !|see algebra.for
        endif  !--------------------+
c Calculate scalar product:
        VecSca2 = rv1(1)*rv2(1)*gb(1,1)+
     +            rv1(2)*rv2(2)*gb(2,2)+
     +            rv1(k)*rv2(k)*gb(3,3)+
     +           (rv1(1)*rv2(2)+rv1(2)*rv2(1))*gb(1,2)+
     +           (rv1(1)*rv2(k)+rv1(k)*rv2(1))*gb(1,3)+
     +           (rv1(2)*rv2(k)+rv1(k)*rv2(2))*gb(2,3)
        return
        end

c=====================================================
        Real Function VecMod2 (rv,k)
c-----------------------------------------------------
c Calculate the length of vector V (indices rv) in a
c crystal with arbitrary symmetry. 
c The real*4 indices version.
c-----------------------------------------------------
        Integer         k
        Real            rv(k), v2

        Real            VecSca2
        External        VecSca2

        v2      = VecSca2 (rv,rv,k)     ! square of |V|
        VecMod2 = sqrt(v2)              ! |V|
        return
        end

c=====================================================
        Real Function  Angle2 (rv1,rv2,k)
c-----------------------------------------------------
c Determine the angle between two planes given by their
c normals or between two vectors in a crystal with
c arbitrary symmetry. 
c The real*4 indices version.
c-----------------------------------------------------
        Integer         k
        Real            rv1(k), rv2(k), d1, d2, dd, s, gra

        Real            VecMod2, VecSca2
        External        VecMod2, VecSca2

        d1 = VecMod2(rv1,k)             !Length of vector V1
        d2 = VecMod2(rv2,k)             !Length of vector V2
        dd = d1*d2

        if (dd.lt.1.e-30) then  !-------------------------------+
           write (0,*) ' Angle2:  V1 = (',rv1,'), |V1| = ',d1  !|
           write (0,*) ' Angle2:  V2 = (',rv2,'), |V1| = ',d2  !|
           write (0,*) ' Angle2:  |V1|*|V2| = ',dd             !|
           stop 'Angle2: |V1|*|V2|=0 or < 1E-30'               !|
        endif  !------------------------------------------------+

        s   = VecSca2(rv1,rv2,k) / dd   !Cosine of angle between V1 & V2
        gra = atan(1.) / 45.            !number of radians in 1 degree

        if (s.gt. 0.999) Then !-------------------------+
                                                       !|
           if (s.ge.1.)     then  !-----------+         |
              Angle2 = 0.                    !|         |
           else  !----------------------------+         |
              Angle2 = Sqrt(2.*(1.-s)) / gra !|         |
           endif  !---------------------------+         |
                                                       !|
        elseif (s.lt.-0.999) Then !---------------------+
                                                       !|
           if (s.le.-1.)    then  !------------------+  |
              Angle2 = 180.                         !|  |
           else  !-----------------------------------+  |
              Angle2 = 180. - Sqrt(2.*(1.+s)) / gra !|  |
           endif  !----------------------------------+  |
                                                       !|
        else  !-----------------------------------------+
                                                       !|
           Angle2 = Atan2(sqrt(1.-s*s),s) / gra        !|
                                                       !|
        endif  !------------------------------------==--+
        return
        end

c=====================================================
        Subroutine VecVec (rv1,rv2,rvv,k)
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
        Integer k, i
        Real    rv1(k), rv2(k), rvv(k), q(3)

        Real    a(6), vol
        Integer lbev
        Common  /back1/ a, vol, lbev

        Real    at(6), g(3,3), gb(3,3)
        Common  /back2/ at, g, gb

c Determine reciprocal lattice:
        if (lbev.eq.0)  then  !-----+
           Call Backlat ()         !|see algebra.for
        endif  !--------------------+
        if (vol.lt.1.E-37) then !----------------------+
           write (0,*) ' VecVec:  a = (',at,')'       !|
           write (0,*) ' VecVec:  vol = ',vol,' <= 0' !|
           stop  'VecVec: vol <= 0'                   !|
        endif !----------------------------------------+
c Vector product:
        q(1) = (rv1(2)*rv2(k)-rv1(k)*rv2(2)) / vol
        q(2) = (rv1(k)*rv2(1)-rv1(1)*rv2(k)) / vol
        q(3) = (rv1(1)*rv2(2)-rv1(2)*rv2(1)) / vol
        do      i=1,3  !============+
           rvv(i) = q(1)*g(1,i)    !|
     +            + q(2)*g(2,i)    !|
     +            + q(3)*g(3,i)    !|
        enddo   !===================+
        if (k.eq.4)     then   !--------+if 4 indices like in hexagonal crystals
           rvv(4) = rvv(3)             !|
           rvv(3) =-(rvv(1)+rvv(2))    !|
        endif  !------------------------+
        return
        end

c=====================================================
        Subroutine      Project (iz,iv,pro,k)
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
        Real            pro(k), z, v

        Real            VecMod, VecSca
        External        VecMod, VecSca

        z = VecMod (iz,k)                       !Length of vector Z
        if (z.lt.1.e-20) then  !------------------------------+
           write (0,*) ' Project:  Z = (',iz,'),  |Z| = ',z  !|
           stop 'Project: Z=0 or < 1E-20'                    !|
        endif  !----------------------------------------------+
        v = VecSca (iv,iz,k) / (z*z)            !(V*Z)/Z*Z
        pro(1) = iv(1)-iz(1)*v                  !projection-1
        pro(2) = iv(2)-iz(2)*v                  !projection-2
        pro(k) = iv(k)-iz(k)*v                  !projection-3
        if (k.eq.4) pro(3)=-(pro(1)+pro(2))
        return
        end

c=====================================================
        Subroutine      Project2 (rz,rv,pro,k)
c-----------------------------------------------------
c Determine projection of vector V (indices rv)
c on a plane normal to vector Z (indices rz) 
c in a crystal with arbitrary symmetry.
c
c  -->   ->   ->  ->   ->   2
c  PRO = V - (V * Z) * Z / Z
c
c The real*4 indices version.
c-----------------------------------------------------
        Integer         k
        Real            rz(k), rv(k), pro(k), z, v

        Real            VecMod2, VecSca2
        External        VecMod2, VecSca2

        z = VecMod2 (rz,k)                      !Length of vector Z
        if (z.lt.1.e-20)  then  !------------------------------+
           write (0,*) ' Project2:  Z = (',rz,'),  |Z| = ',z  !|
           stop 'Project2: Z=0 or < 1E-20'                    !|
        endif  !-----------------------------------------------+
        v = VecSca2(rv,rz,k) / (z*z)            ! (V*Z)/Z*Z
        Call VecSum2(rv,1.,rz,-v,pro,k)         ! projection
        return
        end

c=====================================================
        Subroutine      VecSum2 (v1,c1,v2,c2,v3,k)
c-----------------------------------------------------
c Determine a sum of teo vectors v1 and v2 in a 
c crystal with arbitrary symmetry
c
c     ->      ->      ->
c     v3 = c1*v1 + c2*v2
c-----------------------------------------------------
        Integer         k
        Real            v1(k), v2(k), v3(k), c1, c2

        v3(1) = c1*v1(1) + c2*v2(1)
        v3(2) = c1*v1(2) + c2*v2(2)
        v3(k) = c1*v1(k) + c2*v2(k)
        if (k.eq.4) v3(3) = -(v3(1)+v3(2))
        return
        end

c=====================================================
        Subroutine      UnitVec2 (v,u,k)
c-----------------------------------------------------
c Determine unit-length vector along given vector
c-----------------------------------------------------
        Integer         k, i
        Real            v(k), u(k), vm

        Real            VecMod2
        External        VecMod2

        vm = VecMod2 (v,k)
        if (abs(vm).lt.1.E-32)   return
        do      i=1,k  !=======+
           u(i) = v(i) / vm   !|
        enddo  !===============+
        return
        end

c=====================================================
        Subroutine      VecCon (v,c,u,k)
c-----------------------------------------------------
c Multiply a vector by a constant (change vector length)
c-----------------------------------------------------
        Integer         k, i
        Real            v(k), u(k), c

        do      i=1,k  !=====+
           u(i) = c*v(i)    !|
        enddo  !=============+
        return
        end

c=====================================================
        Subroutine      VecCopy_ir (iv,u,k)
c-----------------------------------------------------
c Copy vector from integer indices to real*4
c-----------------------------------------------------
        Integer         k, iv(k), i
        Real            u(k)

        do      i=1,k  !=====+
           u(i) = iv(i)     !|
        enddo  !=============+
        return
        end

c=====================================================
        Subroutine      VecCopy_rr (v,u,k)
c-----------------------------------------------------
c Copy vector from real*4 to real*4 indices
c-----------------------------------------------------
        Integer         k, i
        Real            v(k), u(k)

        do      i=1,k  !=====+
           u(i) = v(i)      !|
        enddo  !=============+
        return
        end
