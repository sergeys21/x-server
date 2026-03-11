c       program testrand
c       integer i
c       real*8 rand8
c       external rand8
c       do i=1,8
c         write (*,*) rand8(-10.,10.)
c       enddo
c       stop
c       end

        real*8 function rand8(rmin,rmax)
        real   rmin,rmax
        real*8 r
        if (rmin >= rmax) stop 'rand: rmin >= rmax'
c GNU fortran also has a RAND function, but random_number
c is considered superior to it
        call random_number(r)
        rand8 = rmin + (rmax-rmin)*r
        return
        end
