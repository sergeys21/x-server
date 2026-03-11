        Integer Function N_prg_max (filename)
        Integer   lf, lun, i, io_status
        Character filename*(*), buff*16

c Read maximum number of allowed instances of the program
c from file like gid_usr.max, ter_usr.max, trds_usr.max

        Logical*4 FileExist
        External  FileExist

        Integer   FreeLun
        External  FreeLun

        N_prg_max = 0

        lf = Max(Len_Trim(filename),1)

        if (.NOT.FileExist(filename(1:lf))) return

        lun = FreeLun()
        Call OpenFile(filename(1:lf),lun,'read','old',io_status,*1)
  3     continue  !<-----------------------------+
        read(lun,'(a)',err=2,end=2) buff        !|
        i = Index(buff,'#')                     !|strip comment
        if (i.gt.0) buff(i:16)=' '              !|
        i = Index(buff,'!')                     !| strip comment
        if (i.gt.0) buff(i:16)=' '              !|
        if (Len_Trim(buff).eq.0) goto 3 !--------+
        read(buff,'(i5)',err=2,end=2) N_prg_max
  2     continue
        close(unit=lun,err=1)
  1     continue
        if (N_prg_max.lt.0) N_prg_max=0
        return
        end
