        Subroutine DelFile (file,ifail)
c Deletes specified file
        Integer   ifail, lun, io_status
        Character file*(*)
        Logical*4 FileExist
        External  FileExist
        Integer   FreeLun
        External  FreeLun

        ifail = 1
        if (FileExist(file)) Then  !-------------------------------+
c Get the free (not used) lun for the file:                        |
          lun = FreeLun()                                         !|
          Call OpenFile(file,lun,'readwrite','old',io_status,*1)  !|
          Close(unit=lun,                                         !|
     *          status='delete',                                  !|
     *          err=1)                                            !|
        endif  !---------------------------------------------------+
        ifail = 0
  1     continue
        return
        end
