        Subroutine mergeQuery(nkeys,keys,vals,query)
c Merges arrays (keys,vals) into single html-style query
        Integer         nkeys, i, j, l, m, lq

        Character       keys(nkeys)*(*),
     *                  vals(nkeys)*(*),
     *                  query*(*)

        m = Len(query)
        query = ' '
        lq = 0
        do i=1,nkeys  !============================================+
          l = Len_Trim(keys(i))                                   !|
          j = Len_Trim(vals(i))                                   !|
          if (l .eq. 0 .OR. j .eq. 0) goto 1                      !|
          if (vals(i)(1:j) .eq. "toolong") goto 1                 !|
          if (vals(i)(1:j) .eq. "none")    goto 1                 !|
          if (vals(i)(1:j) .eq. "missing") goto 1                 !|
          if (i .eq. 1) then !----------------------------------+  |
            query = keys(i)(1:l)//'='//vals(i)(1:j)            !|  |
          else !------------------------------------------------+  |
            lq = Len_Trim(query)                               !|  |
            if (lq+l+j+2 .le. m) query = query(1:lq)//         !|  |
     *                                   '&'//keys(i)(1:l)//   !|  |
     *                                   '='//vals(i)(1:j)     !|  |
          endif !-----------------------------------------------+  |
  1       continue                                                !|
        enddo !====================================================+
        return
        end
