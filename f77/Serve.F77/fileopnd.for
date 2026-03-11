        Logical*4 Function FileOpened (lun,name)
        Logical*4   fop
        Integer     lun
        Character   name*(*)
        name = ' '
        Inquire (unit=lun, name=name,opened=fop)
        FileOpened = fop
        return
        end

