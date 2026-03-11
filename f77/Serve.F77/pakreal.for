        Subroutine      PakReal (val,nval,txt,nbytes)
c----------------------------------------------------
c Packing real values into the char string.
c Max field for 1 value is 10 bytes.
c----------------------------------------------------
        Integer         nval, nbytes
        Real            val(nval)
        Character       txt*(*)

        Integer         i, i1, j, l, nb,
     *                  ishift, length
        Character       wrk*12

        txt = ' '
        nb  = Len(txt)
        i1  = 1
        Do      i=1,nval  !======================================+
          if (abs(val(i)).lt.1.E-32) Then  !-------------------+ |
c Zero number:                                                !|!|
            wrk    ='0.'                                      !|!|
            length = 2                                        !|!|
          else  !----------------------------------------------+ |
c NON-zero number:                                            !|!|
            write (wrk,'(g12.6)') val(i)                      !| |
c Shift the string to the left:                                | |
            Call  TxtShift (wrk,ishift,length)                !| |
c Search for "E" or "e" letters:                               | |
            Do  j=1,length  !================================+ | |
              if (wrk(j:j).eq.char(69) .OR.                 !| | |
     *            wrk(j:j).eq.char(101)) Then  !-----------+ | | |
c+============+                                            | | | |
c|"E" found!  |                                            | | | |
c+============+                                            | | | |
c 1. Shorten the precision by 2 digits:                    | | | |
                wrk(j-2:j-1) = '00'                       !| | | |
c 2. Insert the sign of the "E"-order:                     | | | |
                if (wrk(j+1:j+1).eq.' ') wrk(j+1:j+1)='+' !| | | |
c 3. Analyze the order:                                    | | | |
                if (wrk(j+2:j+2).ne.'0' .AND.             !| | | |
     *              wrk(j+2:j+2).ne.' ') Then  !--------+  | | | |
c The order contains 2 digits (i.e. it is > 10 or <-10) |  | | | |
c Reduce the precision by one more digit:               |  | | | |
                  wrk(j-3:j-3) = '0'                   !|  | | | |
                else  !---------------------------------|  | | | |
c The order contains 1 digit (i.e. it is < 10 and >-10) |  | | | |
c Compress the order:                                   |  | | | |
                  wrk(j+2:j+2) = wrk(j+3:j+3)          !|  | | | |
                  wrk(j+3:j+3) = ' '                   !|  | | | |
                  length = length-1                    !|  | | | |
                endif  !--------------------------------+  | | | |
c 4. Remove zeros in mantissa after comma (if present):    | | | |
                Do    l=j-1,2,-1  !===============+        | | | |
                  if (wrk(l:l).ne.'0') goto 6 !-+ |        | | | |
                  wrk(l:l) = ' '               !| |        | | | |
                enddo  !========================+=+        | | | |
  6             continue  !<--------------------+          | | | |
                Call TxtShift (wrk(l+1:length),ishift,j)  !| | | |
                goto 40  !--------->--------->---------->--+-+-+-+----+
              endif  !-------------------------------------+ | | |    |
            enddo  !=========================================+ | |    |
c+===============+                                             | |    |
c|No "E" found!  |                                             | |    |
c+===============+                                             | |    |
c No "E" found; remove all zeros after comma:                  | |    v
            Do    j=length,1,-1  !===========+                 | |    |
              if (wrk(j:j).ne.'0') goto 40 !-+-----------------+-+----|
              wrk(j:j) = ' '                !|                 | |    |
            enddo  !=========================+                 | |    |
          endif  !---------------------------------------------+ |    v
                                                                !|    |
  40      continue  !<---------------<-------------<--------<----+----+
c Put the number in the pack string:                             |
          length = Len_Trim (wrk)                               !|
          nbytes = i1+length-1                                  !|
          if (nbytes.gt.nb)             goto 99                 !|
          txt(i1:nbytes) = wrk(1:length)                        !|
          i1     = nbytes + 2  !(leave one space)                |
        enddo  !=================================================+
        return

  99    continue
        nbytes = nb
        return
c       stop 'PakReal: short field'
        end
