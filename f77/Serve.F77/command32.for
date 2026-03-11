        Subroutine      Command (str,ifail)
        Integer         System          !this function is Int*2 in DOS Fortran
        External        System          !and Int*4 in Windows Fortran
        Integer         ifail, l
        Character       str*(*)

        ifail = -1

        l = Len_Trim(str)
        if (l.lt.1) return

        ifail = System(str(1:l))

        return
        end
