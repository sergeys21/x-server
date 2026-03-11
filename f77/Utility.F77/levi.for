        Integer Function  Levi (i,j,k)
        Integer i,j,k
c------------------------------------------
c Antisymmetric Levi-Civita symbol
c
c e_123 = e_231 = e_312 = 1,
c e_132 = e_213 = e_321 = -1,
c All other e_ijk = 0
c
c 123: (1-2)*(2-3)*(3-1)=(-1)*(-1)*(+2)=2
c 132: (1-3)*(3-2)*(2-1)=(-2)*(+1)*(+1)=-2
c------------------------------------------
        Levi = (i-j)*(j-k)*(k-i)
        Levi = Levi/2
        return
        end
