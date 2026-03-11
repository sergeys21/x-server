        Subroutine      tabExpan (text,lines)
c ------------------------------------------------------
c         Expands tabs in text array into spaces
c ------------------------------------------------------
        Integer   lines, line, ifail
        Character text(lines)*(*)

c Loop over lines:
        do      line=1,lines  !================+
                                              !|
          Call tabExpand_ (text(line),ifail)  !|
                                              !|
        enddo  !===============================+
        return
        end
c
c======================================================================
c
        Subroutine      tabExpand_ (text,ifail)
c ------------------------------------------------------
c              Expands tabs into spaces
c ------------------------------------------------------
        Integer   ifail, klen, kt,
     *            ir, it, i1,
     *            m1, m2, i, j
        Character text*(*)

        ifail = 0
        klen  = Len (text)
        if (klen.lt.1) return
        kt    = Len_Trim (text)
        if (kt.lt.1) return

c ******** Loop over characters in the line: *********
c      Do the analysis while the tabs are found and
c  the line expansion does not cause buffer overflow:
        i = 1
        Do while (i.gt.0 .AND. kt.le.klen)  !=====================+
                                                                 !|
          j = i                                                  !|
          i = Index (text(j:kt),char(9))                         !|
          if (i.gt.0) Then  !----------------------------------+  |
                                                              !|  |
          i = i + (j-1)                                       !|  |
c Replacing tabs by spaces:                                    |  |
            ir = ((i-1)/8)*8+8          !"end" of the tab     !|  |
            if (ir.gt.klen) Then !-+                           |  |
              ir    = klen        !|                           |  |
              ifail = 1           !|                           |  |
            endif  !---------------+                           |  |
            it = ir-i                   !text shift value      |  |
            i1 = i+1                    !next symbol after cur.|  |
                                                              !|  |
c If the next symbol is not the last one:                      |  |
            if (i1.le.kt) Then  !-------------------+          |  |
c Shift the text:                                   |          |  |
              do  m1=kt,i1,-1  !=================+  |          |  |
                m2 = m1+it                      !|  |          |  |
c if in limits:                                 !|  |          |  |
                if (m2.le.klen) Then  !------+   |  |          |  |
                  text(m2:m2)=text(m1:m1)   !|   |  |          |  |
                else  !----------------------|   |  |          |  |
                  ifail=2                   !|   |  |          |  |
                endif  !---------------------+   |  |          |  |
              enddo  !===========================+  |          |  |
            endif  !--------------------------------+          |  |
                                                              !|  |
            text(i:ir) = ' '                                  !|  |
            i  = i +it                                        !|  |
            kt = kt+it                                        !|  |
          endif  !---------------------------------------------+  |
                                                                 !|
        enddo  !==================================================+

c Replace non-printable symbols except for page break:
        do i=1,kt  !=====================================+
          if (text(i:i) .lt. ' '       .and.            !|
     *        text(i:i) .ne. char(12)) text(i:i) = ' '  !|
        enddo  !=========================================+

        return
        end
c
c======================================================================
c
        Subroutine      tabReplace (text,lines)
c ------------------------------------------------------
c         Replace tabs in text array by spaces
c ------------------------------------------------------
        Integer   lines, line, kt, i
        Character text(lines)*(*)

c Loop over lines:
        do line=1,lines  !================================+
          kt = Len_Trim (text(line))                     !|
          if (kt.gt.0) Then !--------------------------+  |
  1         continue !<-----------<------------<----+  |  |
              i = Index (text(line)(1:kt),char(9)) !|  |  |
              if (i.gt.0) Then !---------+          ^  |  |
                text(line)(i:i) = ' '   !|          |  |  |
                goto 1  !-------->-------+----->----+  |  |
              endif !--------------------+             |  |
          endif !--------------------------------------+  |
        enddo  !==========================================+
        return
        end
