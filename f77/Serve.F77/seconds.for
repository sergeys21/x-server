        Real    function seconds  (secnd0)
c -----------------------------------------------------
c Wrapper for the SECNDS function (PDP, Microsoft, NDP)
c
c                    A T T E N S I O N :
c
c fOR NDP Fortran compile with the -VMS switch because
c otherwise it interprets SECNDS as INTEGER*4
c -----------------------------------------------------
        Real    secnd0,
     *          secnds,
     *          s,
     *          sutki/86400./           !3600*24

        s = secnds (secnd0)

        if (abs(secnd0).gt.1.E-32) Then !---+
          if (s.lt.-1.)  s = s + sutki     !|
        endif  !----------------------------+

        seconds = s

        return
        end
