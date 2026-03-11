        Subroutine      PakInt4 (int,nint,txt,nbytes)
c Packs integer*4 variable "int" or array "int(nint)"
c into a text string "txt":
        Integer         nint, nbytes, nb, i, i1, ishift, length
        Integer*4       int(nint)
        Character       txt*(*), wrk*10

        txt= ' '
        nb = Len(txt)
        i1 = 1
        Do      i=1,nint  !====================+
          write (wrk,'(i10)')   int(i)        !|
          Call  TxtShift (wrk,ishift,length)  !|
          nbytes = i1+length-1                !|
          if (nbytes.gt.nb)     goto 99       !|
          txt(i1:nbytes) = wrk(1:length)      !|
          i1 = nbytes+2                       !|(leave one space)
        enddo  !===============================+
        return

  99    continue
        nbytes = nb
        return
c       stop 'PakInt4: short field'
        end
