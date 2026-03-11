        Subroutine Rdreal (val,nval,wrk,ifail)
c----------------------------------------------------------
c        This program reads the real*4 array "val"
c                from the text string "wrk"
c
c                   By Sergey Stepanov
c----------------------------------------------------------
        Integer         nval, ifail
        Real            val(nval)
        Character       wrk*(*)

        Integer         i, j, l, nbytes, nb, i_begin,
     *                  i_beg, i_end, lnum, iz, ie
        Real*8          val8

        Integer         lnum_max
        Parameter       (lnum_max=28)
        Character       num*28          !=lnum_max

        do      j=1,nval  !=====+
          val(j) = 0.          !|
        enddo  !================+

        nbytes = Len (wrk)
        ifail  = 0
        nb     = Len_Trim (wrk)
        if (nb.eq.0) return

        Do      i=1,nb  !===================================+
          if (wrk(i:i).eq.char(9))      wrk(i:i)=' '       !|<tab>
          if (wrk(i:i).eq.',')          wrk(i:i)=' '       !|<,>
          if (wrk(i:i).eq.'e')          wrk(i:i)='E'       !|
          if (i.lt.nb)  Then  !--------------------------+  |
            if (wrk(i:i+1).eq.'E ')     wrk(i:i+1)='E+' !|  |
          endif  !---------------------------------------+  |
        Enddo  !============================================+

c-------------------------
c Analyze the input string:
        i_begin = 1
        Do      i=1,nval  !======================================+
c If the starting coordinate of the number is outside the        |
c field, then all the rest of numbers are NULLs:                 |
          if (i_begin.gt.nb)            return                  !|
                                                                !|
c Find the limits of the number in the string:                   |
          i_beg = 0                                             !|
          i_end = 0                                             !|
          Do    j=i_begin,nb  !================================+ |
            If (wrk(j:j).eq.' ') Then  !--------------------+  | |
              If (i_beg.ne.0)    Then  !-+ if number started|  | |
                i_end = j-1             !|                  |  | |
                goto 1          !--------+------------------+--+-+-+
              Endif  !-------------------+                  |  | | |
            Else  !-----------------------------------------+  | | |
              If (i_beg.eq.0)    Then  !-+ if at the begin- |  | | |
                i_beg = j               !| ning of number   |  | | |
              Endif  !-------------------+                  |  | | |
            Endif  !----------------------------------------+  | | |
          Enddo  !=============================================+ | |
          i_end = nb                                            !| |
                                                                !| |
c The number position is found:                                 !| |
  1       continue       !<--------------------------------------+-+
          i_begin = i_end+2                                     !|
          lnum    = i_end-i_beg+1                               !|
          If (lnum.gt.lnum_max) goto 11   !----------------------+--+
          num     = wrk(i_beg:i_end)                            !|  |
          Call TxtShift (num,j,l)                               !|  |
                                                                !|  v
          iz = 0                                                !|
          ie = 0                                                !|
          Do    j=1,lnum  !======================+               |
            If (num(j:j).eq.'.')        iz = 1  !|               |
            If (num(j:j).eq.'E')        ie = j  !|               |
          Enddo  !===============================+               |
                                                                !|
c If the number does not contain ".":                            |
          If (iz.eq.0)  Then  !----------------------------+     |
c Add "." to the number:                                  !|     |
            lnum = lnum + 1                               !|     |
            If (lnum.gt.lnum_max) goto 12   !--------------+-----+--+
                                                          !|     |  |
            If (ie.eq.0) Then   !----------+"E" not found  |     |  v
              num(lnum:lnum) = '.'        !|               |     |
            Else    !----------------------|"E" found      |     |
              ie = ie+1                   !|               |     |
              Do  l=lnum,ie,-1  !========+ |               |     |
                num(l:l) = num(l-1:l-1) !| |               |     |
              Enddo  !===================+ |               |     |
              num(ie-1:ie-1) = '.'        !|               |     |
            Endif  !-----------------------+               |     |
          Endif  !-----------------------------------------+     |
                                                                !|
c For MS DOS Fortran 5.1 (does not work with Linux):             |
c         Read  (num(1:lnum),'(g20.8)',err=13) val(i) !----------+--+
          Read  (num(1:lnum),'(g28.16)',err=13) val8  !----------+--+
          val(i) = sngl(val8)                                   !|  |
c For Linux Fortran g77 (does not work with g77):                |  v
c         Read  (num(1:lnum),*,err=13)  val(i)          !--------+--+
        Enddo  !=================================================+  v
        Return
c-------------------------
c ERRORS:                                                        !v
  11    continue   !<---------------------------------------------+
        ifail = 1
        goto 15  !---------------+                               !v
  12    continue   !<------------+--------------------------------+
        ifail = 2               !|
        goto 15  !---------------|                               !v
  13    continue   !<------------+--------------------------------+
        ifail = 3               !|
        goto 15  !---------------|
                                !|
c ERRORS:                       !|
  15    continue   !<------------+
ccc     Call Beep ()
        return
        end
