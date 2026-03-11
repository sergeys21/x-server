        subroutine round2 (xmin,xmax)
c-----------------------------------------
c       Round interval range
c-----------------------------------------
        Real    xmin, xmax

        Integer n
        Real    xm, dx, x

        Integer  npower
        External npower
        Real     round
        External round

        xm = 0.5*(abs(xmin)+abs(xmax))
        dx = abs(xmax-xmin)
c If the limits differ less than in the 3rd digit:
        if (dx.lt.0.001*xm)     return

c Round lower limit:
        n = npower(abs(xmin)/dx)
        if (n.lt.1)     n=1
  1     continue    !----<------------------------+
        x = round(xmin,n,-1)                     !|
        if (abs(x-xmin).gt.0.2*dx) then  !--+     |
          n = n+1                          !|     ^
          goto 1  !->-----------------------+-----+
        endif  !----------------------------+
        xmin = x

c Round upper limit:
        n = npower(abs(xmax)/dx)
        if (n.lt.1)     n=1
  2     continue    !----<------------------------+
        x = round(xmax,n,+1)                     !|
        if (abs(x-xmax).gt.0.2*dx) then  !--+     |
          n = n+1                          !|     ^
          goto 2   !->----------------------+-----+
        endif  !----------------------------+
        xmax = x

        return
        end
