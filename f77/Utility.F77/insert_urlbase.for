        Subroutine Insert_UrlBase(txt,luno)
c Inserts URL base after the <head> string
        Integer   luno, l, i, m
        Character txt*(*), hdr*64

        Character       Get_Server_Root*64
        External        Get_Server_Root

        m = Len_Trim(txt)
        if (m. gt. 64) m=64
        if (m.ge.6) Then !--------------------------------------+
          hdr = txt(1:m)                                       !|
          Call TxtShift(hdr(1:m),i,l)   !remove left spaces     |
          if (l.       eq. 6  .and.                            !|
     *        hdr(1:1).eq.'<' .and.                            !|
     *        hdr(6:6).eq.'>') Then !-----------------------+   |
            Call Case (hdr(1:l),1)      !to lower case      |   |
            if (hdr(1:l).eq.'<head>') Then !------------+   |   |
              hdr = Get_Server_Root()                  !|   |   |
              m = Len_Trim(hdr)                        !|   |   |
              write (luno,1,err=2) hdr(1:m)            !|   |   |
            endif !-------------------------------------+   |   |
          endif !-------------------------------------------+   |
        endif !-------------------------------------------------+
  1     format('<base href="',a,'" />')
  2     continue
        return
        end
