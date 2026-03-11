        Integer function nonASCII (txt,ipos)
c--------------------------------------------------
c This function detects if there are any non-ASCII
c characters in the text string. Non-ASCII are the
c charactes with the ASCII code > 127 (which also
c includes cyrillic letters and pseudo-graphics).
c If a character is detected, the function returns
c "1" and the character position in the string.
c
c--------------------------------------------------
        Integer    ipos, nchar, i, n

        Character  txt*(*)

        nonASCII = 0
        ipos     = 0
        
        nchar = Len_Trim(txt)
        if (nchar.lt.1) return

        do i=1,nchar !=================+
           n = iChar(txt(i:i))        !|
           if (n .gt. 127) Then !---+  |
              nonASCII = 1         !|  |
              ipos = i             !|  |
              return               !|  |
          endif  !------------------+  |
        enddo  !=======================+

        return
        end
