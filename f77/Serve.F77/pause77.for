        Subroutine Pause ()
        Character byt*1
        write (*,1,advance='no')
  1     format(' Hit ENTER to continue...')
        read (*,'(a)') byt
        return
        end
