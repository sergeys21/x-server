        Integer Function FreeLun ()
        Integer   lun
        Logical*4 Opened

c Get the free (not used) lun:
        lun = 255
  1     continue  !<--------------------+
        Inquire (lun,opened=Opened)    !|
        if (Opened) Then  !----------+  ^
          lun = lun-1               !|  |
          if (lun.gt.0)   goto 1 !---+--+
        endif !----------------------+
        FreeLun = lun
        return
        end

