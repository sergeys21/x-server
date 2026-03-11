        integer function  npower (x)
c-----------------------------------------
c This subroutine determines exponent (order) in the
c floating-point representation of real*4 number, which reads:
c
c  mantissa x 10^exponent
c-----------------------------------------
        real            x
        character       str*14
        write   (str,'(e14.8)') x
        read    (str(12:14),'(i3)',err=1) npower
        return
  1     continue
        write (0,*) 'npower error: x=',x,'  str=',str
        npower = 0
        stop 'npower error'
        end
