        Subroutine      GetDat (iyear, month, iday)
c-------------------------------------------------------
c    This routine for GNU fortran simulates GetDat of Compaq fortran
c ------------------------------------------------------
        Integer         iyear, month, iday
        Integer         values(8)            !year,month,day_of_month,UTC_offset_minutes,hour,minutes,seconds,milliseconds
        Character       datestamp*8,         !  1    2        3              4             5     6       7         8
     *                  timestamp*10,
     *                  zone*5

        Call date_and_time(datestamp,timestamp,zone,values)     !GNU G77

        iyear = values(1)
        month = values(2)
        iday  = values(3)
        return
        end

c ====================================================================

        Subroutine      GetTim (ihour, minute, isecond, i100th)
c-------------------------------------------------------
c    This routine for GNU fortran simulates GetTim of Compaq fortran
c ------------------------------------------------------
        Integer         ihour, minute, isecond, i100th
        Integer         values(8)            !year,month,day_of_month,UTC_offset_minutes,hour,minutes,seconds,milliseconds
        Character       datestamp*8,         !  1    2        3              4             5     6       7         8
     *                  timestamp*10,
     *                  zone*5

        Call date_and_time(datestamp,timestamp,zone,values)     !GNU G77

        ihour   = values(5)
        minute  = values(6)
        isecond = values(7)
        i100th  = INT (values(8)/10.)
        return
        end
