        subroutine sleep_tics (itime,iunits)
c -----------------------------------------------------
c       Function for waiting specified itime
c
c iunits=0 - itime in tics (1/100 of second),
c iunits=1 - itime in seconds
c
c                 by S.Stepanov
c -----------------------------------------------------
        integer         itime,iunits
        integer*4       start,finish,interv,
     *                  h1,h2,m1,s1,t1,it,k1,k6,k36
        integer         hour,minute,second,tic
c
        data    k1,k6,k36       /100,60,360000/
c
        it=itime
        if (iunits.eq.1)        it=it*100
c
        call gettim (hour,minute,second,tic)
        h1      =       hour
        m1      =       minute
        s1      =       second
        t1      =       tic
        start   =       t1+k1*(s1+k6*m1)
c
  1     continue  !<--------------------------------+
          call gettim (hour,minute,second,tic)     !|
          h2      =       hour                     !|
          m1      =       minute                   !|
          s1      =       second                   !|
          t1      =       tic                      !|
          finish  =       t1+k1*(s1+k6*m1)         !|
          interv  =       finish-start             !|
          if (h1.ne.h2)   interv=interv+k36        !|
        if (interv.lt.it) goto 1 !------------------+
c
        return
        end
