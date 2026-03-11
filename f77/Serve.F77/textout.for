        Subroutine      TextOut (txt,lines,iborder)
c--------------------------------------------------------
c   This subroutine outputs text block to the screen
c
c iborder: 0 - output with no frame,
c          1 - output with single-line frame,
c          2 - output with dual-line frame.
c
c                  Author: S.Stepanov
c--------------------------------------------------------
        Integer         lines, iborder,
     *                  lwr, l, l1, l2, r, r1, r2,
     *                  irmk, ltxt, lmax
        Character       txt(lines)*(*),wrk*80

        if (iborder.gt.0)  Then  !--+
          irmk = 2                 !|
        else  !---------------------|
          irmk = 0                 !|
        endif  !--------------------+
        if (lines.lt.1 .or.
     *      lines.gt.25-irmk)        goto 2
        ltxt = Len(txt(1))
        lmax = 0

        do      l=1,lines  !=================+
          lwr = Len_Trim(txt(l)(1:ltxt))    !|
          if (lwr.gt.lmax)      lmax=lwr    !|
        enddo  !=============================+

        if (lmax.lt.1+irmk .or.
     *      lmax.gt.80-irmk)         goto 2
        l1 = 0
        l2 = l1+lines+irmk-1
        r1 = 1
        r2 = r1+lmax+irmk-1

        if (iborder .eq. 0) then  !-----------------------------+
c Output text with no frame:                                    |
          do      l=1,lines   !==================+              |
            wrk(r1:r2) = ' '                    !|              |
            wrk(r1:r2) = txt(l)                 !|              |
            Call Text0 (lmax,wrk(r1:r2))        !|              |
          enddo  !===============================+              |
                                                               !|
        elseif (iborder .eq. 1) then  !-------------------------+
c Single-line frame, upper part:                                |
          do      r=r1+1,r2-1  !===+                            |
            wrk(r:r) = '-'        !|      !char(196)            |
          enddo  !=================+                            |
          wrk(r1:r1) = '+'                !char(218)            |
          wrk(r2:r2) = '+'                !char(191)            |
          Call Text0 (lmax+2,wrk(r1:r2))                       !|
          wrk(r1:r1) = '|'                !char(179)            |
          wrk(r2:r2) = '|'                !char(179)            |
c Output the text:                                              |
          do  l=1,lines !=======================+               |
            wrk(r1+1:r2-1)=' '                 !|               |
            wrk(r1+1:r2-1)=txt(l)              !|               |
            Call  Text0   (lmax+2,wrk(r1:r2))  !|               |
          enddo  !==============================+               |
c Single-line frame, bottom part:                               |
          do      r=r1+1,r2-1  !===+                            |
            wrk(r:r) = '-'        !|      !char(205)            |
          enddo  !=================+                            |
          wrk(r1:r1) = '+'                !char(192)            |
          wrk(r2:r2) = '+'                !char(217)            |
          Call  Text0 (lmax+2,wrk(r1:r2))                      !|
                                                               !|
        elseif (iborder .eq. 2) then  !-------------------------+
c Dual-line frame, upper part:                                  |
          do      r=r1+1,r2-1  !===+                            |
            wrk(r:r) = '='        !|      !char(205)            |
          enddo  !=================+                            |
          wrk(r1:r1) = '+'                !char(201)            |
          wrk(r2:r2) = '+'                !char(187)            |
          Call Text0 (lmax+2,wrk(r1:r2))                       !|
          wrk(r1:r1) = '|'                !char(186)            |
          wrk(r2:r2) = '|'                !char(186)            |
c Output the text:                                              |
          do      l=1,lines  !==================+               |
            wrk(r1+1:r2-1)=' '                 !|               |
            wrk(r1+1:r2-1)=txt(l)              !|               |
            Call  Text0   (lmax+2,wrk(r1:r2))  !|               |
          enddo  !==============================+               |
c Dual-line frame, the bottom part:                             |
          do      r=r1+1,r2-1  !===+                            |
            wrk(r:r) = '='        !|      !char(205)            |
          enddo  !=================+                            |
                                                               !|
          wrk(r1:r1) = '+'                !char(200)            |
          wrk(r2:r2) = '+'                !char(188)            |
          Call  Text0 (lmax+2,wrk(r1:r2))                      !|
        endif  !------------------------------------------------+

c----------------------------------------------------
        Return
c----------------------------------------------------
  2     Stop ' TEXTOUT: Params not in range'
        End
