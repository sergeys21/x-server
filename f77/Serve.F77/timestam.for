        Subroutine      TimeStam (txtime)
c       Character       txtime*5
        Character       txtime*(*)
        Integer         ihr,imin,isec,i100th
c-------------------------------------------------------
        txtime = ' '
        if (Len(txtime) .lt. 5) Then !-------------------------------+
           write (*,*) 'TimeStam: too short timestamp field passed' !|
           return                                                   !|
        endif  !-----------------------------------------------------+

        Call    GetTim (ihr,imin,isec,i100th)

        write (txtime(1:2),'(i2)')      ihr
        txtime(3:3)=':'
        write (txtime(4:5),'(i2)')      imin
        if (imin.lt.10) txtime(4:4)='0'

        return
        end
