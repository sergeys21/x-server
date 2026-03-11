        Subroutine      RdInt   (ival,nval,wrk,ifail)
c-----------------------------------------------------
c        This program reads the integer array "ival"
c                from the text string "wrk"
c
c                   By Sergey Stepanov
c-----------------------------------------------------
        Integer         nval, ival(nval), ifail,
     *                  nb, ns, i, i1, j
        Character       wrk*(*)

        ifail = 0
        nb = Len_Trim (wrk)
        if (nb.gt.0) Then  !----------------------------------+
                                                             !|
c Preliminary transformations:                                |
          do      i=1,nb  !======================+            |
          if (wrk(i:i).eq.char(9)) wrk(i:i)=' ' !|            |<tab>
          if (wrk(i:i).eq.',')     wrk(i:i)=' ' !|            |<,>
          enddo  !===============================+            |
                                                             !|
c Analyze the text string:                                    |
          i1 = 1                                             !|
          do      i=1,nval   !===========================+    |
c If the starting coordinate of the number is outside the|    |
c field, then all the rest of numbers are NULLs:         |    |
            if (i1.gt.nb)   goto 9  !---------->---------+-+  |
c Shift the number all the way to the left:              | |  |
            Call  TxtShift (wrk(i1:nb),ns,j)            !| |  |
c If the remaining bytes are spaces, then                | |  |
c the remaining numbers are NULL:                        | |  |
            if (j.eq.0)  goto 9  !------------->---------+-|  |
c Find the number end in the string and read the number: | |  |
            do    j=i1,nb  !=================+           | |  |
              if (wrk(j:j).eq.' ') goto 4 !--+-+         | |  |
            enddo  !=========================+ |         | |  |
            j = nb+1                          !|         | v  |
  4         continue   !<----------------------+         | |  |
            read (wrk(i1:j-1),'(i10)',err=11) ival(i) !--+-+--+-+
            i1 = j                                      !| |  | |
          enddo  !=======================================+ |  | |
                                                          !|  | |
          return                                          !|  | |
                                                          !|  | |
        endif  !-------------------------------------------+--+ |
                                                          !|    |
c Values that are not present considered to be all NULLs   v    |
        i = 1                                             !|    |
  9     continue   !<--------------------------------------+    v
        do      j=i,nval !==+                                   |
          ival(j) = 0      !|                                   |
        enddo  !============+                                   |
        return                                                 !|
                                                               !|
c Beep on error and return to string editing:                   |
  11    continue   !<-------------------------------------------+
c       Call    Beep ()
        ifail = 1
        return
        end
