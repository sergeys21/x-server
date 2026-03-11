        integer function  npower8 (x8)
c-----------------------------------------
c This subroutine determines exponent (order) in the
c floating-point representation of real*8 number, which reads:
c
c  mantissa x 10^exponent
c-----------------------------------------
        real*8          x8
        character       str*14
        write   (str,'(e14.8)') x8
        read    (str(12:14),'(i3)',err=1) npower8
        return
  1     continue
        write (0,*) 'npower8 error: x8=',x8,'  str=',str
        npower8 = 0
        stop 'npower8 error'
        end
