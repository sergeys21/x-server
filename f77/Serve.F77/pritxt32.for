        subroutine      pritxt  (txt,icr,ilf,itimo,ifail)
c--------------------------------------------------------
c           This subroutime prints text string  
c txt - text string,     
c icr - carriage return flag,
c ilf - number carriage returns,
c itimo - timeout flag (0 - no timeout, 1 - use current timeout)
c ifail - error status on return
c--------------------------------------------------------
        Integer   icr, ilf, itimo, ifail
        Character txt*(*)

        Integer   lun, l, i
        Logical*4 Opened

        Integer numprn
        Common /prnnum/ numprn

        Integer  nowrite
        External nowrite

c Suppress "unused arguments" warnings:
        i = icr
        i = itimo
        i = ifail

        if (numprn*(6-numprn).lt.0)  stop 'Illegal PRN No'

c 1-3 - lpt1:lpt3, 4-5 - com1:com2  6 - file (lun=66)

        l   = Len_Trim(txt)
        lun = 66
c------------------------------------------------
c Is the file open?
        Inquire (lun,opened=Opened)
        if (.NOT.Opened) Then !-----------------+
          Open (lun, file='Dump.out', err=8)   !|
        endif !---------------------------------+

c------------------------------------------------
c Print into file:
        i = ilf

        if (l.gt.0)   Then  !--------------------+
  4       continue                              !|
          write (lun,'(a)',err=2)  txt(1:l)     !|
          i=i-1                                 !|
          goto 3                                !|
  2       continue                              !|
          if (nowrite().eq.0) goto 4            !|
          goto 8                                !|
  3       continue                              !|
        endif  !---------------------------------+

        if (i.gt.0)   Then  !--------------------+
          do  l=1,i  !=======================+   |
  7         continue                        !|   |
            write (lun,*,err=5)             !|   |
            goto 6                          !|   |
  5         continue                        !|   |
            if (nowrite().eq.0) goto 7      !|   |
            goto 8                          !|   |
  6         continue                        !|   |
          enddo  !===========================+   |
        endif  !---------------------------------+
  8     continue
        return
        end
