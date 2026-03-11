        Logical*4 Function FileReady (file)

        Character file*(*), bit*1
        Integer   lun, io_status

        Logical*4 FileExist
        External  FileExist

        Integer   FreeLun
        External  FreeLun

        FileReady = .False.

        if (.NOT.FileExist(file)) return

        lun = FreeLun()

        Call OpenFile(file,lun,'read','old',io_status,*99)
c Check that there is at least one byte in the file
c (i.e. that the file is not empty):
        read (lun,'(a)',err=2,end=2) bit

        FileReady = .True.

  2     continue
        close (unit=lun,err=99)
  99    continue
        return
        end

