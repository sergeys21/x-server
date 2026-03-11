        Subroutine ChiParse (Code,ncompMax,ncomp,Atoms,kol,ifail)
c--------------------------------------------------------
c   Parses a substance chemical formula into elements
c   and the numbers of their occurances in the compound.
c
c                     By S.Stepanov
c--------------------------------------------------------
c An example: Cd2H5(Pt3Cl4)5
c--------------------------------------------------------
        Integer         ncompMax, kol(ncompMax),
     *                  ncomp, ifail
        Character       code*(*), atoms(ncompMax)*(*)

        Character       txt(20)*80
        Common  /msg/   txt

        Integer lc, is, il, loop, loopfactor, i, j, k, m

        ifail = 0
        lc    = Len_Trim (code)
        if (lc.eq.0)    goto 11
        do i=1,lc  !=============================================+
c Syntax control (the only symbols allowed):                     |
          if (code(i:i).eq.'^')     code(i:i)=' '               !|
          if (code(i:i).eq.'_')     code(i:i)=' '               !|
          if (code(i:i).eq.char(9)) code(i:i)=' '               !|
          if (code(i:i).ne.'+' .AND.                            !|
     *        code(i:i).ne.'-' .AND.                            !|
     *        code(i:i).ne.' ' .AND.                            !|
     *        code(i:i).ne.'(' .AND.                            !|
     *        code(i:i).ne.')' .AND.                            !|
     *       (code(i:i).lt.'0' .OR. code(i:i).gt.'9') .AND.     !|
     *       (code(i:i).lt.'A' .OR. code(i:i).gt.'Z') .AND.     !|
     *       (code(i:i).lt.'a' .OR. code(i:i).gt.'z')) goto 12  !|
        enddo  !=================================================+
c Remove all internal spaces:
  1     continue  !<---------------------------+
        i = Index(code(1:lc),' ')             !|
        if (i.gt.0) Then  !-----------------+  ^
          Call TxtShift(code(i:lc),is,il)  !|  |
          lc = lc-is                       !|  |
          goto 1   !---------->-------------+--+
        endif  !----------------------------+

        i = Len(code)
        if (lc.lt.i) code(lc+1:i) = ' '

c Begin parsing:
c =============
        ncomp = 0
        i     = 1
        loop  = 0               !() loop is closed
  2     continue
        if (i.gt.lc)    goto 10

c "(" -- Begin of loop in the formula:
        if (code(i:i).eq.'(') Then  !----+
          if (loop.eq.0)     Then  !--+  |
            loop = 1       !open loop!|  |
          else !----------------------|  |
            goto 13     !already open!|  |
          endif  !--------------------+  |
          i = i+1                       !|
          goto 2                        !|
        endif  !-------------------------+

c ")" -- End of loop in the formula:
        if (code(i:i).eq.')') Then  !------------------------+
          if (loop.eq.1)  Then  !-------------------------+  |
            loop = 0                 !close the loop     !|  |
            i    = i+1                                   !|  |
c What is the loop repetition?                           !|  |
            if (i.gt.lc) Then  !-----------------------+  |  |
              loopfactor = 1        !end of space, no !|  |  |
              i  = i-1              !factor specified !|  |  |
              goto 3                                  !|  |  |
            else  !------------------------------------|  |  |
              if ((code(i:i).ge.'A'  .AND.            !|  |  |
     *             code(i:i).le.'Z') .OR.             !|  |  |
     *            (code(i:i).eq.'(')) Then  !--------+ |  |  |
                loopfactor = 1      !Another loop or | |  |  |
                goto 3              !element start.  | |  |  |
                                    !No factor given.| |  |  |
              elseif (code(i:i).ge.'0' .AND.        !| |  |  |
     *                code(i:i).le.'9') Then  !------| |  |  |
                Call GetInteger (code,i,   !read fac-| |  |  |
     *                 loopfactor,ifail)   !tor value| |  |  |
                if (ifail.ne.0) goto 14             !| |  |  |
                goto 3                              !| |  |  |
              else  !--------------------------------| |  |  |
                goto 14   !unexpected symbol after ) | |  |  |
              endif  !-------------------------------+ |  |  |
            endif  !-----------------------------------+  |  |
  3         continue                                     !|  |
            k = 0                  !Nr of elements in loop|  |
            do j=1,ncomp !=============================+  |  |
              if (kol(j).lt.0)  Then  !-------------+  |  |  |
                kol(j) = Abs(kol(j)) * loopfactor  !|  |  |  |
                k      = k+1          !element found|  |  |  |
              endif   !-----------------------------+  |  |  |
            enddo  !===================================+  |  |
            if (k.eq.0) goto 14  !No elements inside loop |  |
          else !------------------------------------------|  |
            goto 14              !loop was not open before|  |
          endif  !----------------------------------------+  |
          i = i+1                                           !|
          goto 2                                            !|
        endif  !---------------------------------------------+

c Begin of element name:
        if (code(i:i).ge.'A' .AND. code(i:i).le.'Z') Then  !-+
          ncomp = ncomp+1                                   !|
          if (ncomp.gt.ncompMax) goto 15     !buffer overflow|
          Atoms(ncomp) = code(i:i)                          !|
          m            = 1                                  !|
          if (loop.eq.0) Then !--+                           |
            kol(ncomp) = 1      !|outside the loop           |
          else  !----------------|                           |
            kol(ncomp) = -1     !|inside the loop            |
          endif  !---------------+                           |
          i = i+1                                           !|
          if (i.gt.lc) goto 10                              !|
c For 2-symbol elements (like Pt, Li...):                    |
          if (code(i:i).ge.'a' .AND.                        !|
     *        code(i:i).le.'z') Then  !-----------+          |
            m = m + 1                            !|          |
            Atoms(ncomp)(m:m) = code(i:i)        !|          |
            i = i+1                              !|          |
            if (i.gt.lc) goto 10                 !|          |
          endif   !-------------------------------+          |
c For ions (Li+, H-, F-...):                                 |
          if (code(i:i).eq.'+'  .OR.                        !|
     *        code(i:i).eq.'-') Then  !-----------+          |
            m = m + 1                            !|          |
            Atoms(ncomp)(m:m) = code(i:i)        !|          |
            i = i+1                              !|          |
            if (i.gt.lc) goto 10                 !|          |
          endif   !-------------------------------+          |
          goto 2                                            !|
        endif  !---------------------------------------------+

c Begin of element repetition:
        if (code(i:i).ge.'0' .AND.
     *      code(i:i).le.'9') Then !-------------------------+
          if (ncomp.eq.0) Then  !-----------------+          |
            goto 14     !no component to apply   !|          |
          else  !---------------------------------|          |
            Call GetInteger (code,i,j,ifail)     !|          |
            if (ifail.ne.0)     goto 14  !read err|          |
            if (loop.eq.0) Then !--+              |          |
              kol(ncomp) = j      !|outside a loop|          |
            else  !----------------|              |          |
              kol(ncomp) = -j     !|inside a loop |          |
            endif  !---------------+              |          |
          endif  !--------------------------------+          |
          i = i+1                                           !|
          goto 2                                            !|
        endif  !---------------------------------------------+

c Nothing else can be expected:
        goto 14

c Exit:
  10    continue                        !Normal exit!!!
c Some more controls:
        if (loop.eq.1)  goto 14         !loop not closed!
        do i=1,ncomp  !=====================+
          if (kol(i).lt.1)      goto 14    !|kol's not defined!
        enddo  !============================+
  999   return
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
  13    txt(1)='Enclosed loops not allowed in material formula'
        ifail = 1
        goto 999
  14    txt(1)='Syntax error in material formula'
        ifail = 1
        goto 999
  15    write (txt(1),25)       ncompMax
  25    format (
     *  'Too many components in material formula [max=',i3,']')
        ifail = 1
        ncomp = ncompmax
        goto 999
        end

c ==================================================================

        subroutine GetInteger(code,icurr,number,ifail)
        Integer         icurr, number, ifail,
     *                  lc, i
        Character       code*(*)

        ifail = 1
        lc    = Len_Trim(code)
        if (icurr.gt.lc)           goto 2
        do i=icurr+1,lc !===============+
          if (code(i:i).lt.'0' .OR.    !|
     *        code(i:i).gt.'9') goto 1 !+-+
        enddo  !========================+ |
        i = lc+1                         !v
  1     continue  !<-------------<--------+
        i = i-1
c The repetition number cannot exceed 99:
        if (i-icurr+1.gt.2)        goto 2
        read (code(icurr:i),'(i3)',err=2)       number
        ifail = 0
        icurr = i
  2     return
        end
