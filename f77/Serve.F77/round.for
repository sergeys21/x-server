        real    function  round (x,n,lowup)
c-----------------------------------------
c Round a number to "n" (1:6)  significant digit
c lowup=-1  - round towards a smaller number,
c lowup=+1  - round towards a larger number.
c-----------------------------------------
        Integer n, lowup
        Real    x

        Integer i
        Real    q, y, z

        integer  npower
        External npower

        if ((n-1)*(6-n).lt.0 .or. lowup*lowup.ne.1) then !-+
          round = x                                       !|
          return                                          !|
        endif  !-------------------------------------------+
c The order of the number:
        i = npower(x)
        q=1.
        if (n.ne.i)     q=10.**(n-i)
        y = x*q
c Round to "n" significant digits:
        z = anint(y)

        if (lowup.eq.+1 .and. z.lt.y)   z=z+1
        if (lowup.eq.-1 .and. z.gt.y)   z=z-1

        round=z/q

        return
        end
