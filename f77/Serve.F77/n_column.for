        Integer Function n_Column (wrk)
c-----------------------------------------------------
c         This routine calculates the number
c       of numeric columns in given text string
c
c               By Sergey Stepanov
c-----------------------------------------------------
        Integer         nb, i
        Character       wrk*(*),byt*1

        nb = Len_Trim (wrk)
        if (nb.eq.0)    Then  !----+
          n_Column = 0            !|
          return                  !|
        endif  !-------------------+

        do      i=1,nb  !==============================+
          if (wrk(i:i).eq.char(9))  wrk(i:i)=' '      !|
          if (wrk(i:i).eq.',')      wrk(i:i)=' '      !|
          if (wrk(i:i).eq.'e')      wrk(i:i)='E'      !|
          if (i.lt.nb)  Then  !----------------------+ |
            if (wrk(i:i+1).eq.'E ') wrk(i:i+1)='E+' !| |
          endif  !-----------------------------------+ |
        enddo  !=======================================+

c-------------------------
c Analyze the input string:
        if (wrk(1:1).ne.' ') Then !---+
          n_Column = 1               !|
          byt = 'N'                  !|
        else  !-----------------------|
          n_Column = 0               !|
          byt = ' '                  !|
        endif  !----------------------+

        do      i=2,nb   !=====================+
          if (wrk(i:i).ne.' ')  Then  !-----+  |
            if (byt.eq.' ')     Then  !--+  |  |
c Beginning of column:                   |  |  |
              n_Column = n_Column+1     !|  |  |
              byt = 'N'                 !|  |  |
            endif  !---------------------+  |  |
          else  !---------------------------|  |
            byt = ' '                      !|  |
          endif  !--------------------------+  |
        enddo  !===============================+
        return
        end
