        Subroutine Duration (start,bias,hh,mm,ss)
        Real            start, Seconds
        Integer*4       hh, mm, ss, time, bias
        External        Seconds

        if (abs(start).lt.1.E-32) Then !---------+
                                                !|
          Start  = Seconds (0.)                 !|
          Bias   = 0                            !|
          hh = 0                                !|
          mm = 0                                !|
          ss = 0                                !|
                                                !|
        else   !---------------------------------|
                                                !|
          time = Int(Seconds(start))            !|
                                                !|
          if (time.gt.72000)  Then  !--+ 20*3600 |
            bias = bias + time        !|         |
            time = 0                  !|         |
            start = Seconds(0.)       !|         |
           endif  !--------------------+         |
                                                !|
          ss = time + bias                      !|
          mm = ss/60                            !|
          ss = ss-60*mm                         !|
          hh = mm/60                            !|
          mm = mm-60*hh                         !|
                                                !|
        endif  !---------------------------------+

        return
        end
