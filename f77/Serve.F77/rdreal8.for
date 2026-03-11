        Subroutine      Rdreal8 (val,nval,wrk,ifail)
c-----------------------------------------------------
c        This program reads the real*8 array "val"
c                from the text string "wrk"
c
c                   By Sergey Stepanov
c-----------------------------------------------------
        Integer         nval, ifail
        Real*8          val(nval)
        Character       wrk*(*)

        Integer         i, j, l, nbytes, nb, i_begin,
     *                  i_beg, i_end, lnum, iz, ie

        Integer         lnum_max
        Parameter       (lnum_max=28)
        Character       num*28          !=lnum_max

        do      j=1,nval  !=====+
          val(j) = 0.          !|
        enddo  !================+

        nbytes = Len (wrk)
        ifail  = 0
        nb     = Len_Trim (wrk)
        if (nb .eq. 0)  return

        Do      i=1,nb  !===================================+
          if (wrk(i:i).eq.char(9))      wrk(i:i)=' '       !|<tab>
          if (wrk(i:i).eq.',')          wrk(i:i)=' '       !|<,>
          if (wrk(i:i).eq.'e')          wrk(i:i)='D'       !|
          if (wrk(i:i).eq.'E')          wrk(i:i)='D'       !|
          if (i.lt.nb)  Then  !--------------------------+  |
            if (wrk(i:i+1).eq.'D ')     wrk(i:i+1)='D+' !|  |
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
            If (wrk(j:j).eq.' ') Then  !---------------------+ | |
              If (i_beg.ne.0)    Then  !-+if number started  | | |
                i_end = j-1             !|                   | | |
                goto 1          !--------+-------------------+-+-+-+
              Endif  !-------------------+                   | | | |
            Else  !------------------------------------------+ | | |
              If (i_beg.eq.0)    Then  !-+if at the beginning| | | |
                i_beg = j               !|                   | | | |
              Endif  !-------------------+                   | | | |
            Endif  !-----------------------------------------+ | | |
          Enddo  !=============================================+ | |
          i_end = nb                                            !| |
                                                                !| |
c The number position is found:                                 !| |
  1       continue       !<--------------------------------------+-+
          i_begin = i_end+2                                     !|
          lnum    = i_end-i_beg+1                               !|
          If (lnum.gt.lnum_max)         goto 11   !--------------+--+
          num     = wrk(i_beg:i_end)                            !|  |
          Call  TxtShift (num,j,l)                              !|  |
                                                                !|  v
          iz = 0                                                !|
          ie = 0                                                !|
          Do    j=1,lnum  !======================+               |
            If (num(j:j).eq.'.')        iz = 1  !|               |
            If (num(j:j).eq.'D')        ie = j  !|               |
          Enddo  !===============================+               |
                                                                !|
c If the number does not contain ".":                            |
          If (iz.eq.0)  Then  !----------------------------+     |
c Add "." to the number:                                  !|     |
            lnum = lnum + 1                               !|     |
            If (lnum.gt.lnum_max)       goto 11   !--------+-----+--+
                                                          !|     |  |
            If (ie.eq.0) Then   !----------+"D" not found  |     |  v
              num(lnum:lnum) = '.'        !|               |     |
            Else    !----------------------+"D" found      |     |
              ie = ie+1                   !|               |     |
              Do  l=lnum,ie,-1  !========+ |               |     |
                num(l:l) = num(l-1:l-1) !| |               |     |
              Enddo  !===================+ |               |     |
              num(ie-1:ie-1) = '.'        !|               |     |
            Endif  !-----------------------+               |     |
          Endif  !-----------------------------------------+     |
                                                                !|
          Read  (num(1:lnum),'(g28.16)',err=11) val(i)  !--------+--+
        Enddo  !=================================================+  v
        Return
c-------------------------
c Beep on error and return to string editing:                       v
  11    continue   !<-----------------------------------------------+
c       Call    Beep ()
        ifail = 1
        return
        end
