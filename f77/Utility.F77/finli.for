        Real Function Finli (x,f,xn,dx,n)
        Integer n, j
        Real    x, f(n), xn, dx, xk, z, fz, f1, q

        Finli = 0.
        xk = xn + dx*(n-1)
c       if (x.lt.Amin1(xn,xk))          return
c       if (x.gt.Amax1(xn,xk))          return
        z = Abs(dx)
        if (x.lt.Amin1(xn,xk)-z)        return
        if (x.gt.Amax1(xn,xk)+z)        return
        z  = (x-xn)/dx + 1.
        j  = int(z)
        if (j.lt.1)     j = 1
        if (j.gt.n-1)   j = n-1
        q  = z-j
        fz = f(j)
        f1 = f(j+1)-fz
        Finli = fz + q*f1
        return
        end
