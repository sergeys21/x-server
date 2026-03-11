        Subroutine GetFileDate_NN (file,idate)

c This is version for GNU Fortran

        Character       file*(*)
        Integer         buff(13), timevals(9), status, ll, i
        Integer*2       idate(12)

c    1     2      3      4     5      6
c  iyr_w imon_w iday_w ihr_w imin_w isec_w  -- last write array
c    7     8      9     10    11     12
c  iyr_a imon_a iday_a ihr_a imin_a isec_a  -- last access array
        do      i=1,12 !===+
          idate(i) = 0    !|
        enddo  !===========+
        ll = Len_Trim(file)
        if (ll .eq. 0) return

        Call stat(file(1:ll), buff, status)             !GNU fortran function

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

c Parse  modification time:
        Call ltime(buff(10),timevals)

        idate(1)  = int2(timevals(6)+1900)
        idate(2)  = int2(timevals(5)+1)
        idate(3)  = int2(timevals(4))
        idate(4)  = int2(timevals(3))
        idate(5)  = int2(timevals(2))
        idate(6)  = int2(timevals(1))

c Parse access time:
        Call ltime(buff(9),timevals)

        idate(7)  = int2(timevals(6)+1900)
        idate(8)  = int2(timevals(5)+1)
        idate(9)  = int2(timevals(4))
        idate(10) = int2(timevals(3))
        idate(11) = int2(timevals(2))
        idate(12) = int2(timevals(1))

        return
        end
