        Logical*4 Function FileExist (filename)
c ---------------------------------------------------
c .TRUE.  If the specified file exists and can be opened
c .FALSE. If the specified file does not exist or if the file exists but cannot be opened
c ---------------------------------------------------
        Character filename*(*)
        Integer   lun, io_status

        Integer   FreeLun
        External  FreeLun

        lun = FreeLun()
        FileExist = .False.

        Call OpenFile(filename,lun,'read','old',io_status,*99)
        FileExist = .True.
        close (unit=lun,err=99)
  99    continue
        return
        end

