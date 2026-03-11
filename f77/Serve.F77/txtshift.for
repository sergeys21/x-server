        Subroutine txtShift (txt, nshift, nchar)
c-----------------------------------------------------
c     Removes spaces at the beginning of string
c         and shifts the string to the left
c nshift - the value for which the text is shifted.
c nchar  - the length of txt string after the shift.
c-----------------------------------------------------
        Integer         nshift, nchar,
     *                  nc, i, j
        Character       txt*(*)

        nshift = 0
        nc = Len_Trim(txt)
        nchar  = nc
        if (nc.le.1)            return
        if (txt(1:1).ne.' ')    return
        do      i=2,nc  !-------------------+
          if (txt(i:i).ne.' ')  goto 2  !---+--+
        enddo  !----------------------------+  |
  2     continue   !<--------------------------+
        nshift = i-1
        nchar  = nc-nshift
        do      j=i,nc  !--------------------+
          txt(j-nshift:j-nshift)=txt(j:j)   !|
        enddo  !-----------------------------+
        txt(nc-nshift+1:nc) = ' '
        return
        end
