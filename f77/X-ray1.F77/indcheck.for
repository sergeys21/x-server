        integer function indcheck (ind,inmax)
        integer         inmax, ind(inmax)
        character       txt(20)*80
        common  /msg/   txt
        indcheck=0
        if (inmax.eq.3) return
        if (ind(3).eq.-(ind(1)+ind(2))) return
        txt(1)=' ** ERROR **'
        txt(2)=' '
        txt(3)=' ** Illegal indexes: i#-(h+k)! **'
        call    message (txt,3,2)
        indcheck=1
        return
        end
