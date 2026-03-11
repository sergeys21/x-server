        Subroutine N4_Pkr0 (number,tx)
        Integer*4 number
        Integer   nt, n, i, j
        Character tx*(*)
c--------------------------------------------------
c This routine packs int4 numbers at the right side
c of text string and adds nulls at the left:
c                   0000123
c--------------------------------------------------
        nt = Len(tx)
        Call    PakInt4 (number,1,tx,n)
        if (n.lt.nt)        Then  !--+
          do    i=n,1,-1  !======+   |
            j = i+(nt-n)        !|   |
            tx(j:j) = tx(i:i)   !|   |
          enddo  !===============+   |
          j = nt-n                  !|
          do i=1,j !==========+      |
            tx(i:i) = '0'    !|      |
          enddo  !============+      |
        endif  !---------------------+
        return
        end
