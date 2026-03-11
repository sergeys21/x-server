        subroutine flushoutput (lun)
        integer lun
c This is currently supported in GNU fortran only
        call flush(lun)
        return
        end
