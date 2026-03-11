        Subroutine List_INP (luntbl, InpFile, File_Top, wrk)
c Copy INP and PRF files into the TBL:
        Integer   luntbl, luninp,
     *            linp, ltop, i,
     *            io_status
        Character InpFile*(*),
     *            File_Top*(*),
     *            wrk*(*)

        Integer   FreeLun
        External  FreeLun

        luninp = Freelun()
        linp   = Len_Trim(InpFile)
        ltop   = Len_Trim(File_Top)

        Call OpenFile(InpFile(1:linp),luninp,'read','old',
     *                                     io_status,*301)

        write (luntbl,396,err=398)
        write (luntbl,395,err=398) InpFile(1:linp)
        write (luntbl,396,err=398)

  397   continue  !<------------------------------+
        read (luninp,'(a)',err=398,end=398) wrk  !|
        i = Max(Len_Trim(wrk),1)                 !|
        write(luntbl,'(a)',err=398) wrk(1:i)     !|
        goto 397  !-------------------------------+
  398   continue
        Close (unit=luninp,err=301)
  301   continue

        if (Len_Trim(File_Top) .gt. 0) Then  !--------------------------+
          wrk = File_Top(1:ltop)//'.prf'                               !|
          Call OpenFile(wrk(1:ltop+4),luninp,'read','old',             !|
     *                                     io_status,*401)             !|
          write (luntbl,396,err=400)                                   !|
          write (luntbl,395,err=400) wrk(1:ltop+4)                     !|
          write (luntbl,396,err=400)                                   !|
  399     continue  !<--------------------------------+                 |
          read (luninp,'(a)',err=400,end=400) wrk    !|                 |
          i = Max(Len_Trim(wrk),1)                   !|                 |
          write(luntbl,'(a)',err=400) wrk(1:i)       !|                 |
          goto 399  !---------------------------------+                 |
  400     continue                                                     !|
          Close (unit=luninp,err=401)                                  !|
        endif  !--------------------------------------------------------+
  401   continue
        write (luntbl,396,err=402)
        write (luntbl,*,err=402)
  402   continue
        return

  395   format(20x,'Contents of file: ',a)
  396   format(79('-'))
        end
