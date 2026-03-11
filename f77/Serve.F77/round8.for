        real*8  function  round8 (x,n,lowup)
c-----------------------------------------
c Round a number to "n" (1:6)  significant digit
c lowup=-1  - round towards a smaller number,
c lowup=+1  - round towards a larger number.
c-----------------------------------------
        Integer n, lowup
        Real*8  x

        Integer i
        Real*8  q, y, z

        integer  npower8
        External npower8

        if ((n-1)*(6-n).lt.0 .or. lowup*lowup.ne.1) then !-+
          round8 = x                                       !|
          return                                          !|
        endif  !-------------------------------------------+
c The order of the number:
        i = npower8(x)
        q = 1.
        if (n .ne. i)   q = 10.0D0**(n-i)
        y = x*q
c Round to "n" significant digits:
        z = dnint(y)

        if (lowup.eq.+1 .and. z.lt.y)   z=z+1
        if (lowup.eq.-1 .and. z.gt.y)   z=z-1

        round8 = z/q

        return
        end
