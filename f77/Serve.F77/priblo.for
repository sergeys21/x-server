        subroutine      priblo  (txt,lines,lf1,lf2,itimo,ifail)
c--------------------------------------------------------
c       This subroutine prints a block of text (multi-line)
c txt(lines) - text block (array),
c lines      - number of lines,
c lf1        - number of empty lines to print before block,
c lf2        - number of empty lines to print after block,
c itimo      - timeout flag,
c ifail      - error status on exit,
c--------------------------------------------------------
        integer         lines,lf1,lf2,itimo,ifail,i
        character       txt(lines)*(*)

        ifail=0

        if (lf1.gt.0)   then  !---------------------+
          Call  Pritxt  (' ',1,lf1,itimo,ifail)    !|
          if (ifail.ne.0)       return             !|
        endif  !------------------------------------+

        do      i=1,lines   !=======================+
          Call  Pritxt  (txt(i),1,1,itimo,ifail)   !|
          if (ifail.ne.0)       return             !|
        enddo  !====================================+

        if (lf2.lt.1)   return
        Call    Pritxt  (' ',1,lf2,itimo,ifail)
        return
        end
