        subroutine      text0   (nsy,txt)

        integer         nsy, l
        character       txt*(*)

        if (nsy.eq.0)   return

        l = min(len(txt),nsy)
        if (l.gt.0)   write (*,'(1x,a)') txt(1:l)

        return
        end
