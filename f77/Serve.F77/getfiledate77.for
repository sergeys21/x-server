        Subroutine GeFilDat (file,datestamp,timestamp)

c This is version for GNU Fortran

        Character       file*(*), datestamp*10, timestamp*5
        Integer         buff(13), timevals(9), status, ll, i

        datestamp = ' '
        timestamp = ' '

        ll = Len_Trim(file)
        if (ll .eq. 0) return

        Call stat(file(1:ll), buff, status)

        if (status .ne. 0) return

c The array "buff" structure:
c  1. Device ID
c  2. Inode number
c  3. File mode (octal)
c  4. Number of links
c  5. Owner''s uid
c  6. Owner''s gid
c  7. Device where located
c  8. File size
c  9. Last access time          [need to use CTIME for stamp or LTIME for values]
c 10. Last modification time    [need to use CTIME for stamp or LTIME for values]
c 11. Last status change time   [need to use CTIME for stamp or LTIME for values]
c 12. Preferred block size
c 13. No. of blocks allocated

c Parse  modification time:
        Call ltime(buff(10),timevals)

c The array "timevals" structure:
c  1. Seconds after the minute, range 0–59 or 0–61 to allow for leap seconds
c  2. Minutes after the hour, range 0–59
c  3. Hours past midnight, range 0–23
c  4. Day of month, range 0–31
c  5. Number of months since January, range 0–12
c  6. Years since 1900
c  7. Number of days since Sunday, range 0–6
c  8. Days since January 1
c  9. Daylight savings indicator: positive if daylight savings is in effect, zero if not, and negative if the information is not available.

        write   (datestamp,1) timevals(4),timevals(5)+1,timevals(6)+1900
  1     format  (i2,'.',i2,'.',i4)
        write   (timestamp,2) timevals(3),timevals(2)
  2     format  (i2,':',i2)
        do i=1,len(datestamp)  !===========================+
           if (datestamp(i:i).eq.' ') datestamp(i:i)='0'  !|
        enddo !============================================+
        do i=1,len(timestamp)  !===========================+
           if (timestamp(i:i).eq.' ') timestamp(4:4)='0'  !|
        enddo !============================================+
        return
        end
