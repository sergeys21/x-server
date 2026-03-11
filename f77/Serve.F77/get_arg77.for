        subroutine      get_arg_wrapper  (n,arg)
c ----------------------------------------------------------------------
c This is a wrapper to FORTRAN getarg subroutine which is different
c between FORTRAN compilers:
c -- MS DOS Fortran-5 always needs 3-argument form:
c       call    getarg  (n,arg,istat)
c -- Compaq Fortran takes both 3-argument and 2-argument forms:
c       call    getarg  (n,arg,istat)
c       call    getarg  (n,arg)
C GNU Fortran takes only 2-arguments and also offers a Fortran-2003 standard sub:
c       call    getarg  (n,arg)
c       call    get_command_argument(n,arg)
c ----------------------------------------------------------------------
        character       arg*(*)
        integer         n

        arg = ' '
        call      getarg  (n,arg)
        return
        end
