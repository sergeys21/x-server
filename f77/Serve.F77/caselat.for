        subroutine      CaseLat (txt,icase)
c-----------------------------------------------
c This routine transforms ENGLISH text in the
c string "txt" into the upper or lower case.
c
c icase=0 - transform "txt" into the UPPER case
c icase=1 - transform "txt" into the LOWER case
c
c                 By Sergey Stepanov
c-----------------------------------------------
        integer         icase,
     *                  m1, m2,
     *                  ltxt, i, j
        character       txt*(*)

        ltxt = Len_Trim(txt)

        if (icase.eq.1) goto 2

c Lower -> Upper
        m1 = 97                 !ichar('a')
        m2 = 122                !ichar('z')
        do      1       i=1,ltxt
          j = ichar(txt(i:i))
          if (j.lt.m1)  goto 1
          if (j.gt.m2)  goto 1
          j = j-32
          txt(i:i) = char(j)
  1     continue
        return

c Upper -> Lower
  2     m1 = 65                 !ichar('A')
        m2 = 90                 !ichar('Z')
        do      3       i=1,ltxt
          j = ichar(txt(i:i))
          if (j.lt.m1)  goto 3
          if (j.gt.m2)  goto 3
          j = j+32
          txt(i:i) = char(j)
  3     continue
        return
        end
