c
c This file contains the following subroutines:
c            1. LineLoop,
c            2. Compress_Line,
c            3. ArgSearch,
c            4. ChiSearch
c
c They help to implement data input for:
c
c gid_sl, ter_sl, trds_sl, mag_sl, brls, and other programs
c
c                    Author: S.Stepanov
c ============================================================


        Subroutine      LineLoop (lun,line,wrk,*)
        Integer         lun, line, i, j
        Character       wrk*(*)

  1     line = line+1  !<------------------------------+
        read    (lun,'(40x,a)',err=100,end=101) wrk   !^
c       Call txtShift (wrk,j,i)                       !|i=Lengh(wrk)
c Skip commented lines:                                |
        if ((wrk(1:1).eq.';') .or.                    !|
     *      (wrk(1:1).eq.'!')) goto 1 !----->----------+

        i = Len_Trim (wrk)
        if (i.eq.0) return
c Remove comments at the end of line:
        j = Index(wrk(1:i),';')
        if (j.gt.0) wrk(j:i) = ' '
        j = Index(wrk(1:i),'!')
        if (j.gt.0) wrk(j:i) = ' '

        i = Len_Trim (wrk)
        if (i.gt.0) Call txtShift (wrk,j,i)        !i=Lengh(wrk)
        return
  100   continue
  101   continue
        return 1
        end

c =======================================================

        Subroutine Compress_Line (buffer, lenbuf, string, ioffset)
        Character buffer*(*), string*(*)
        Integer   lenbuf, ioffset, i, j
c Example:
c  Call Compress_Line (buffer(3),lenbuf,'  ',1)
c  Call Compress_Line (buffer(3),lenbuf,' =',0)
c  Call Compress_Line (buffer(3),lenbuf,'= ',1)
  1     continue  !<---------------------------------+
        if (lenbuf.lt.2) return                     !|
c Find the string to be                              |
c compressed in the buffer:                          |
        i = Index(buffer(1:lenbuf),string)          !|
        if (i.gt.0)     Then  !------------------+   |
c Compress the buffer starting at                |   |
c the string position + ioffset:                 |   |
          do j=i+ioffset,lenbuf-1  !=========+   |   ^
           buffer(j:j) = buffer(j+1:j+1)    !|   |   |
          enddo  !===========================+   |   |
         buffer(lenbuf:lenbuf) = ' '            !|   |
         lenbuf = lenbuf-1                      !|   |
         goto 1  !-------------------------------+---+
        endif  !---------------------------------+
        return
        end

c =======================================================

        Subroutine ArgSearch (buffer, keyword,
     *                        stopbyte1, stopbyte2,
     *                        keybeg, iargbeg, iargend)
c------------------------------------------------------
c This program parameters:
c   buffer    -(Input)-  buffer where to search for arguments
c   keyword   -(Input)-  argument keyword
c   stopbyte1 -(Input)-  delimiter of type-1 to search for argument end
c   stopbyte2 -(Input)-  delimiter of type-2 to search for argument end
c   keybeg    -(OUTput)- position of key begin (0=not found, -1=duplicate)
c   iargbeg   -(OUTput)- position of argument begin
c   iargend   -(OUTput)- position of argument end
c------------------------------------------------------
        Character buffer*(*), keyword*(*),
     *            stopbyte1*1, stopbyte2*1
        Integer   keybeg, iargbeg, iargend,
     *            lenbuf, lenkey, keybeg2, kstart

        lenbuf  = Len(buffer)
        lenkey  = Len(keyword)
        keybeg  = 0
        keybeg2 = 0
        iargbeg = 0
        iargend = 0
        if (lenbuf.lt.1) return
        keybeg = Index(buffer(1:lenbuf),keyword)
c Avoid the case that the key is a part of bigger key
c like t= and F1t= may both be triggered by t= search
c (added 2003/07/01):
        if (keybeg.gt.1) Then !-----------------------------+
          if (buffer(keybeg-1:keybeg-1).ne.' ' .AND.       !|
     *        buffer(keybeg-1:keybeg-1).ne.',') keybeg=0   !|
        endif !---------------------------------------------+
c If key found:
        if (keybeg .gt. 0) Then  !--------------------------------------+
c Search for duplicate key:                                             |
          kstart = keybeg                                              !|
  2       continue   !<---------------------------------------------+   |
          if (kstart .lt. lenbuf) Then !-----------------------+    |   |
            keybeg2 = Index(buffer(kstart+1:lenbuf),keyword)  !|    |   |
          endif  !---------------------------------------------+    |   |
c Avoid the case that the key is a part of bigger key               |   |
c like t= and F1t= may both be triggered by t= search               |   |
c (added 2003/07/01):                                               |   |
          if (keybeg2.gt.0) Then !-----------------------------+    |   |
            keybeg2 = keybeg2+kstart                          !|    |   |
            if (buffer(keybeg2-1:keybeg2-1).ne.' ' .AND.      !|    |   |
     *          buffer(keybeg2-1:keybeg2-1).ne.',') Then !-+   |    |   |
c False alarm, keep searching:                             |   |    |   |
              kstart  = kstart + keybeg2                  !|   |    |   |
              keybeg2 = 0                                 !|   |    |   |
              goto 2 !-------------------------------------+---+----+   |
            endif !----------------------------------------+   |        |
c Duplicate key found!                                         |        |
          endif !----------------------------------------------+        |
          if (keybeg2 .eq. 0) Then  !--------------------------------+  |
c No duplicates:                                                     |  |
            iargbeg = keybeg+lenkey                                 !|  |
            do iargend=iargbeg,lenbuf   !========================+   |  |
              if (buffer(iargend:iargend).eq.stopbyte1 .or.     !|   |  |
     *            buffer(iargend:iargend).eq.stopbyte2) goto 1 !-+-+ |  |
            enddo !==============================================+ | |  |
            iargend = lenbuf                                      !v |  |
  1         continue  !<-------------------------------------------+ |  |
          else  !----------------------------------------------------|  |
c Duplicate key found:                                               |  |
            keybeg = -1                                             !|  |
          endif  !---------------------------------------------------+  |
        endif  !--------------------------------------------------------+
        return
        end

c =======================================================

        Subroutine ChiSearch (buffer, keyword, chi, C0h, wrk, oh, *)
        Complex*8       chi
        Real            xx(2)
        Integer         i, m, n, l, ifail, iedge,
     *                  keybeg, iargbeg, iargend
        Character       buffer*(*), keyword*(*), C0h*(*), wrk*(*), oh*1

        Logical*4       FileOpened
        External        FileOpened

        chi = (0.0, 0.0)
        c0h = ' '
        l   = Len(wrk)

c ADDED 2003/07/25: handling resonance
c cases for soft x-rays when x0r might
c be positive!
c Check if there is an "edge" flag
c (if not, x0r must be strictly < 0)
        iedge = 0
        Call ArgSearch (buffer,'edge=',' ',',',
     *                             keybeg, iargbeg, iargend)
        if     (keybeg .gt. 0)  Then  !-----------------------+
          wrk = 'edge'                                       !|
          Call rdInt (iedge,1,buffer(iargbeg:iargend),ifail) !|
          if (ifail .ne. 0)        goto 2 !-------------------+---ERROR+
          buffer(keybeg:iargend) = ' '  !erase the key!       |        v
        elseif (keybeg .eq. -1) Then  !-----------------------|
          wrk = 'duplicate keyword "edge="'                  !|
          goto 2 !------------------------------->------------+---ERROR+
        endif  !----------------------------------------------+        v

        Call ArgSearch (buffer, keyword, ')', ')',
     *                             keybeg, iargbeg, iargend)

        if (keybeg .gt. 0)  Then  !---------------------------+
c Format the beginning of error message string, e.g.: "x0: "  |
          wrk      = keyword                                 !|
          m        = Len_Trim(wrk)                           !|
          wrk(m:l) = ': '                                    !|
          m        = m+2                                     !|
c Check for () and remove () if needed:                      !|
c Error message:  "x0: mising (). Use: x0=(xr,xi)"            |
          wrk(m:l) = 'missing (). Use:'                      !|
          n        = Len_Trim(wrk)+2                         !|
          wrk(n:l) = keyword//'(xr,xi)'                      !|
                                                             !|
          if (buffer(iargbeg:iargbeg).ne.'(') goto 2 !--------+---ERROR+
          if (buffer(iargend:iargend).ne.')') goto 2 !--------+---ERROR+
          iargbeg  = iargbeg + 1                             !|
          iargend  = iargend - 1                             !|
                                                             !|
c Check for the string consistency:                           |
c Error message:  "x0: parse error"                           |
          wrk(m:l) = 'parse error'                           !|
          if (iargend .lt. iargbeg)  goto 2 !-----------------+---ERROR+
                                                             !|
c Process the inputs like: (a+i*b) or (a-i*b):               !|This requires
c         n = Index(buffer(iargbeg:iargend),'+i*')           !|more accurate
c         if (n.gt.0) buffer(iargbeg+n-1:iargbeg+n+1)=', +'  !|treating of
c         n = Index(buffer(iargbeg:iargend),'-i*')           !|spaces!!!
c         if (n.gt.0) buffer(iargbeg+n-1:iargbeg+n+1)=', -'  !|
                                                             !|
c Read xr,xi:                                                !|
c Error message:  "x0: syntax error. Use: x0=(xr,xi)"         |
          wrk(m:l) = 'syntax error. Use:'                    !|
          n        = Len_Trim(wrk)+2                         !|
          wrk(n:l) = keyword//'(xr,xi)'                      !|
          Call rdReal (xx,2,buffer(iargbeg:iargend),ifail)   !|
          if (ifail .ne. 0)  goto 2 !-------------------------+---ERROR+
                                                             !|
          iargend  = iargend + 1        !restore iargend      |
          buffer(keybeg:iargend) = ' '  !erase the key!       |
c Erase comma after the argument:                             |
          n = Len(buffer)                                    !|
          if (iargend .lt. n) Then  !-------------------+     |
            i = iargend+1                              !|     |
            if (buffer(i:i).eq.',') buffer(i:i) = ' '  !|     |
          endif  !--------------------------------------+     |
                                                             !|
c Check for big (xr,xi):                                     !|
          wrk(m:l) = 'too big value (>0.1). '//              !|
     *               'Not in the medium x-ray range?'        !|
          do n=1,2  !================================+        |
            if (Abs(xx(n)) .gt. 0.1) goto 2 !--------+--------+---ERROR+
          enddo  !===================================+        |        v
                                                             !|
          if (oh.eq.'0') then !---------------------------+   |
c x0 (reflection): normally the x0r is negative           |   |
c ADDED 2003/07/25: handling resonance                    |   |
c cases for soft x-rays when x0r might                    |   |
c be positive!                                            |   |
            if (iedge.ne.1) Then !------------------+     |   |
              chi = Cmplx(-Abs(xx(1)),Abs(xx(2)))  !|     |   |
            else !----------------------------------|     |   |
              chi = Cmplx(xx(1),Abs(xx(2)))        !|     |   |
            endif  !--------------------------------+     |   |
          else !------------------------------------------|   |
c xh (diffraction): the sign of xhr is handled by xdf!    |   |
            chi = Cmplx(+Abs(xx(1)),Abs(xx(2)))          !|   |
          endif !-----------------------------------------+   |
c Report in the TBL file that x0r is positive                 |
c (if the TBL file is opened):                                |
          if (oh.eq.'0' .AND. xx(1).gt.0.) then !-----------+ |
            if (FileOpened(3,wrk)) then !-----------------+ | |
              l = Max(4,Len_Trim(wrk))                   !| | |
              Call Case (wrk(1:l),1)       !to lower case | | |
              if (wrk(l-3:l).eq.'.tbl') then !----------+ | | |
                if (iedge.ne.1) Then !----------------+ | | | |
                  write(3,11,err=3) xx(1)            !| | | | |here oh is either "0" or "h"
  11              format(' ** WARNING: x0r = ',g12.5,
     *            ' is specified as positive --',1x,
     *            'inverted')                        !| | | | |
                else !--------------------------------| | | | |
                  write(3,12,err=3) xx(1)            !| | | | |
  12              format(' ** WARNING: x0r = ',g12.5,
     *            ' is specified as positive --',1x,
     *            'preserved due to edge=1 flag')    !| | | | |
                endif  !------------------------------+ | | | |
              endif  !----------------------------------+ | | |
            endif  !--------------------------------------+ | |
          endif  !------------------------------------------+ |
  3       continue                                           !|
          C0h      = 'x'                                     !|
                                                             !|
        elseif (keybeg .eq. -1)  Then  !----------------------|
                                                             !|
          wrk      = keyword                                 !|
          m        = Len_Trim(wrk)                           !|
          wrk(m:l) = ': duplicate keyword'                   !|
          goto 2 !--------------------------------------------+---ERROR+
                                                             !|        v
        endif  !----------------------------------------------+
        return
c Return on error:
  2     continue
        return 1
        end
