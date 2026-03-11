        Subroutine      CurrentTime (txdate,txtime)
c-------------------------------------------------------
c    This routine fills up txdate and txtime with
c        current date and time, respectively
c-------------------------------------------------------
c       Character       txdate*8,txtime*5
        Character       txdate*(*),txtime*(*)
        Integer         date, month, year,
     *                  hour, minute, second, i100th

        txdate = ' '
        txtime = ' '
        if (Len(txdate) .lt. 8) Then !----------------------------------+
           write (*,*) 'CurrentTime: too short datestamp field passed' !|
           return                                                      !|
        endif  !--------------------------------------------------------+
        if (Len(txtime) .lt. 5) Then !----------------------------------+
           write (*,*) 'CurrentTime: too short timestamp field passed' !|
           return                                                      !|
        endif  !--------------------------------------------------------+

        Call GetDat (year, month, date)
        Call GetTim (hour, minute, second, i100th)

c Current date in the text form:
        write (txdate(5:8), '(i4)')  year
        txdate(5:6) = ' .'
        write (txdate(4:5), '(i2)')  month
        if (txdate(4:4).eq.' ') txdate(4:4) = '0'
        txdate(3:3) = '.'
        write (txdate(1:2), '(i2)')  date
c       if (txdate(1:1).eq.' ') txdate(1:1) = '0'

c Current time in the text form:
        write (txtime(4:5), '(i2)')  minute
        if (txtime(4:4).eq.' ') txtime(4:4) = '0'
        txtime(3:3) = ':'
        write (txtime(1:2), '(i2)')  hour
c       if (txtime(1:1).eq.' ') txtime(1:1) = '0'

        return
        end
