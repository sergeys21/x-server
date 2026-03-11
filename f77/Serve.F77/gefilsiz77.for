        Subroutine GeFilSiz (filename,filesize)

c This is version for GNU Fortran

        Character       filename*(*)
        Integer*4       filesize
        Integer         fsize
        Logical         ex

        filesize = 0

c Attention: we can also use the 8th element of the array returned
c by the "stat" function. See getfiledate77.for for additional detals.

        inquire(file=filename, size=fsize, exist=ex)

        if (.NOT.ex)    return                !file not found

        filesize = fsize
        return
        end
