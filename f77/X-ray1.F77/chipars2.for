        Subroutine ChiParse2 (Code,ncompMax,ncomp,Atoms,rkol,ifail)
c--------------------------------------------------------
c   Parses a substance chemical formula into elements
c   and the numbers of their occurances in the compound.
c
c    *** Version with NON-integer number of atoms ***
c
c                     By S.Stepanov
c--------------------------------------------------------
c An example: Cd2H5(Pt3Cl4)5
c  +-----------+
c  |ATTENTION: +-------------------------------------+
c  |This version allows real indices in the formula, |
c  |like:  Al0.3Ga0.7As +----------------------------+
c  +--------------------+
c--------------------------------------------------------
        Integer         ncompMax, ncomp, ifail
        Real            rkol(ncompMax), rloopfactor, rpt
        Character       Code*(*), Atoms(ncompMax)*(*),
     *                  Formula*40

        Character       txt(20)*80
        Common  /msg/   txt

        Integer         lc, lf, is, il, loop, i, j, k, m

        ifail = 0
        lc    = Len_Trim (code)
        lf    = Len (Formula)
        if (lc.eq.0)    goto 11
        if (lc.gt.lf)   goto 16
        Formula = Code(1:lc)
        do i=1,lc  !==================================================+
c Syntax control (the only symbols allowed):                          |
          if (Formula(i:i).eq.'^')     Formula(i:i)=' '              !|
          if (Formula(i:i).eq.'_')     Formula(i:i)=' '              !|
          if (Formula(i:i).eq.char(9)) Formula(i:i)=' '              !|
          if (Formula(i:i).ne.'+' .AND.                              !|
     *        Formula(i:i).ne.'-' .AND.                              !|
     *        Formula(i:i).ne.' ' .AND.                              !|
     *        Formula(i:i).ne.'.' .AND.                              !|
     *        Formula(i:i).ne.'(' .AND.                              !|
     *        Formula(i:i).ne.')' .AND.                              !|
     *       (Formula(i:i).lt.'0' .OR. Formula(i:i).gt.'9') .AND.    !|
     *       (Formula(i:i).lt.'A' .OR. Formula(i:i).gt.'Z') .AND.    !|
     *       (Formula(i:i).lt.'a' .OR. Formula(i:i).gt.'z')) goto 12 !+-ERROR+
        enddo  !======================================================+      v

c Remove all internal spaces:
  1     continue  !<------------------------------+
        i = Index(Formula(1:lc),' ')             !|
        if (i.gt.0) Then  !--------------------+  ^
          Call TxtShift(Formula(i:lc),is,il)  !|  |
          lc = lc-is                          !|  |
          goto 1   ! --------->----------------+--+
        endif  !-------------------------------+

        if (lc.lt.lf) Formula(lc+1:lf) = ' '

c Begin parsing:
c =============
        ncomp = 0
        i     = 1
        loop  = 0               !() loop is closed
  2     continue   !<---------------------<-------------------<------------+
        if (i.gt.lc)    goto 10  !--------->-------+--Normal Exit!----+    |N
                                                  !|                  |    |E
c "(" -- Begin of loop in the formula:             |                  v    |X
        if (Formula(i:i).eq.'(') Then  !----+      |                  |    |T
          if (loop.eq.0)     Then  !--+     |      ^                  |E   ^
            loop = 1       !open loop!|     |      |                  |X   |
          else !----------------------|     |      |                  |I   |
            goto 13    !already open!-+-----+------+-ERROR-+          |T   |
          endif  !--------------------+     |      |       v          |    |
          i = i+1                          !|      |                  |    |
          goto 2  !-------------------------+------+                  |    |
        endif  !----------------------------+                         |    |N
                                                                     !|    |E
c ")" -- End of loop in the formula:                                  |    |X
        if (Formula(i:i).eq.')') Then  !-----------------------+      |    |T
          if (loop.eq.1)  Then  !---------------------------+  |      v    |
            loop = 0                 !close the loop       !|  |      |    ^
            i    = i+1                                     !|  |      |    |
c What is the loop repetition?                             !|  |      |    |
            if (i.gt.lc) Then  !-----------------------+    |  |      |    |
c End of formula text. No factor specified:            |    |  |      |E   |
              rloopfactor = 1                         !|    |  |      |X   |
              i = i-1                                 !|    |  |      |I   |
              goto 3  !--------->------------>---------+-+  |  |      |T   |
            else  !------------------------------------| |  |  |      |    |
              if ((Formula(i:i).ge.'A'  .AND.         !| v  |  |      |    |
     *             Formula(i:i).le.'Z') .OR.          !| |  |  |      |    |
     *            (Formula(i:i).eq.'(')) Then !------+ | |  |  |      v    |N
c Start of another element or another loop.          | | |  |  |      |    |E
c No factor specified:                               | | |  |  |      |    |X
                rloopfactor = 1                     !| | |  |  |      |    |T
                goto 3  !------->------------>-------+-+-|  |  |      |    |
              elseif ((Formula(i:i).ge.'0' .AND.    !| | |  |  |      |    ^
     *                 Formula(i:i).le.'9') .OR.    !| | |  |  |      |    |
     *                (Formula(i:i).eq.'.')) Then !--| | |  |  |      |E   |
c Read the group repetition:                         | | |  |  |      |X   |
                Call GetReal(Formula,i,rloopfactor, !| | |  |  |      |I   |
     *                                       ifail) !| | v  |  |      |T   |
                if (ifail.ne.0) goto 14 !------------+-+-+--+--+-ERROR+--+ |
                goto 3  !------->------------>-------+-+-|  |  |      |  v |
              else  !--------------------------------| | v  |  |      |    |
                goto 14 !unexpected symb.after ")" !-+-+-+--+--+-ERROR+--+ |
              endif  !-------------------------------+ | |  |  |      |  v |
            endif  !-----------------------------------+ v  |  |      v    |
  3         continue  !<-----------<---------------<-----+  |  |      |    |
            k = 0                  !Nr of elements in loop  |  |      |    |N
            do j=1,ncomp !==============================+   |  |      |    |E
              if (rkol(j).lt.0.)  Then  !-------------+ |   |  |      |    |X
                rkol(j) = Abs(rkol(j)) * rloopfactor !| |   |  |      |    |T
                k       = k+1          !element found | |   |  |      v    |
              endif   !-------------------------------+ |   |  |      |    ^
            enddo  !====================================+   |  |      |    |
            if (k.eq.0) goto 14  !No elements inside loop !-+--+-ERROR+--+ |
          else !--------------------------------------------|  |      |  v |
            goto 14              !loop was not open before!-+--+-ERROR+--+ |
          endif  !------------------------------------------+  |      |  v |
          i = i+1                                             !|      |    |
          goto 2  !--------------->------------------->--------+------+->--|
        endif  !-----------------------------------------------+      |    |N
                                                                     !|    |E
c Begin of element name:                                              |    |X
        if (Formula(i:i).ge.'A' .AND.                                !|    |T
     *      Formula(i:i).le.'Z') Then  !-----------------+            |    |
          ncomp = ncomp+1                               !|            |    |
          if (ncomp.gt.ncompMax) goto 15 !buffer overflow+-------ERROR+--+ |
          Atoms(ncomp) = Formula(i:i)                   !|            |  v |
          m            = 1                              !|            |    |
          if (loop.eq.0) Then !--+                       |            |    |
            rkol(ncomp) = 1.    !|outside the loop       |            |    ^
          else  !----------------|                       |            |    |
            rkol(ncomp) = -1.   !|inside the loop        |            |    |
          endif  !---------------+                       |            v    |N
          i = i+1                                       !|            |    |E
          if (i.gt.lc) goto 10  !---------->-------------+->------->--|    |X
c For 2-symbol elements (like Pt, Li...):                |            |E   |T
          if (Formula(i:i).ge.'a' .AND.                 !|            |X   |
     *        Formula(i:i).le.'z') Then  !--------+      |            |I   |
            m = m + 1                            !|      |            |T   |
            Atoms(ncomp)(m:m) = Formula(i:i)     !|      |            |    |
            i = i+1                              !|      |            |    |
            if (i.gt.lc) goto 10 !--------->------+------+->------->--|    |
          endif   !-------------------------------+      |            |    |
c For ions (Li+, H-, F-...):                             |            |    |
          if (Formula(i:i).eq.'+'  .OR.                 !|            |    |
     *        Formula(i:i).eq.'-') Then  !--------+      |            |    |
            m = m + 1                            !|      |            |    |
            Atoms(ncomp)(m:m) = Formula(i:i)     !|      |            |    |
            i = i+1                              !|      |            v    ^
            if (i.gt.lc) goto 10 !--------->------+------+->------->--|    |
          endif   !-------------------------------+      |            |    |
          goto 2   !-------------->------------------->--+------------+->--|
        endif  !-----------------------------------------+            |E   |N
                                                                     !|X   |E
c Begin of element repetition:                                        |I   |X
        if ((Formula(i:i).ge.'0' .AND.                               !|T   |T
     *       Formula(i:i).le.'9') .OR.                               !|    |
     *      (Formula(i:i).eq.'.')) Then !--------------+              |    |
          if (ncomp.eq.0) Then  !-----------------+    |              |    |
            goto 14     !no component to apply----+----+--------ERROR-+--+ |
          else  !---------------------------------+    |              |  v |
c Read the element reperition:                    |    |              |    |
            Call GetReal (Formula,i,rpt,ifail)   !|    |              v    |
            if (ifail.ne.0) goto 14  !-read error-+----+--------ERROR-+--+ |
            if (loop.eq.0) Then !--+              |    |              |  v |
              rkol(ncomp) = rpt   !|outside a loop|    |              |    |
            else  !----------------|              |    |              |    |
              rkol(ncomp) =-rpt   !|inside a loop |    |              |    |
            endif  !---------------+              |    |              |    |N
          endif  !--------------------------------+    |              |    |E
          i = i+1                                     !|              |    |X
          goto 2  !---------->--------------->------->-+--------------+->--+T
        endif  !---------------------------------------+              |
                                                                     !|
c Nothing else can be expected:                                       |
        goto 14  !------------>----------------->---------------ERROR-+--+
                                                                     !|  v
c Exit:                                                               |
  10    continue     !<------------------Normal exit!-----------------+
c Some more controls:
c Check if loop is closed:
        if (loop.ne.0)          goto 14  !----------------------ERROR----+
c Check if all rkol's are defined:                                       v
        do i=1,ncomp  !=====================+
          if (rkol(i).le.0.)    goto 14  !--+-------------------ERROR----+
        enddo  !============================+                            v
  999   continue
        return
c========================================================
c ERRORS:
  11    continue
        txt(1)='No material formula specified!'
        ifail = -1
        goto 999

  12    continue
        txt(1)='Illegal character(s) in material formula'
        ifail = 1
        goto 999

  13    continue
        txt(1)='Enclosed loops not allowed in material formula'
        ifail = 1
        goto 999

  14    continue
        txt(1)='Syntax error in material formula'
        ifail = 1
        goto 999

  15    continue
        write (txt(1),25)       ncompMax
  25    format (
     *  'Too many components in material formula [max=',i3,']')
        ifail = 1
        ncomp = ncompmax
        goto 999

  16    continue
        write (txt(1),26)       Code(1:lc), lf
  26    format (
     *  'Too long structure code=[',a,'] exceeding ',i3,' chars buffer')
        ifail = 1
        goto 999

        end

c ==================================================================

        subroutine GetReal(stroke,icurr,rnumber,ifail)
        Integer         icurr, ifail, lc, i, j
        Real            rnumber
        Character       stroke*(*), tmp*8

        ifail = 1
        lc    = Len_Trim(stroke)
        if (icurr.gt.lc)           goto 2
c Find the end of the number:
        do i=icurr+1,lc !=======================+
          if ((stroke(i:i).lt.'0'  .OR.        !|
     *         stroke(i:i).gt.'9') .AND.       !|
     *         stroke(i:i).ne.'.' )   goto 1 !--+-+
        enddo  !================================+ |
        i = lc+1                                 !v
  1     continue  !<--------------<---------------+
        i = i-1
        tmp = stroke(icurr:i)
c Do we have "." in the number?
        j = Index(tmp,'.')
        if (j.eq.0) Then  !----------+
          j = Len_Trim(tmp)+1       !|
          if (j.gt.Len(tmp)) goto 2 !|
          tmp(j:j) = '.'            !|
        endif  !---------------------+
        read (tmp,'(f8.4)',err=2)  rnumber
        if (rnumber.le.0.) goto 2
        ifail = 0
        icurr = i
  2     continue
        return
        end
