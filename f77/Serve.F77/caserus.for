        subroutine      CaseRus (txt,icase)
c-----------------------------------------------
c This routine transforms the RUSSIAN text in the
c string "txt" into the upper or lower case.
c
c icase=0 - transform "txt" into the UPPER case
c icase=1 - transform "txt" into the LOWER case
c
c                 By Sergey Stepanov
c-----------------------------------------------
        integer         icase,
     *                  m1, m2, m3, m4,
     *                  ltxt, i, j
        character       txt*(*)

        ltxt = Len_Trim(txt)
        if (icase.eq.1) goto 2

c Lower -> Upper
        m1 = 160                        !ichar('à')
        m2 = 175                        !ichar('ï')
        m3 = 224                        !ichar('ð')
        m4 = 239                        !ichar('ÿ')
        do      1       i=1,ltxt
          j = ichar(txt(i:i))
          if (j.lt.m1)  goto 1
          if (j.gt.m4)  goto 1
          if (j.le.m2)  then
            j = j-32
            txt(i:i) = char(j)
            goto 1
          endif
          if (j.ge.m3)  then
            j = j-80
            txt(i:i) = char(j)
            goto 1
          endif
  1     continue
        return

c Upper -> Lower
  2     m1 = 128                        !ichar('A')
        m2 = 143                        !ichar('Ï')
        m3 = 144                        !ichar('Ð')
        m4 = 159                        !ichar('ß')
        do      3       i=1,ltxt
          j = ichar(txt(i:i))
          if (j.lt.m1)  goto 3
          if (j.gt.m4)  goto 3
          if (j.le.m2)  then
            j = j+32
            txt(i:i) = char(j)
            goto 3
          endif
          if (j.ge.m3)  then
            j = j+80
            txt(i:i) = char(j)
            goto 3
          endif
  3     continue
        return
        end
