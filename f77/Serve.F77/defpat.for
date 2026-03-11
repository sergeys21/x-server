        subroutine      defpat  (path)
c    +----------------------------------+
c    |  Determines the directory from   |
c    |   which the program is called    |
c    +----------------------------------+
        character       path*(*)
        integer         i, l
        character       wrk*80

        Character       getOS*8
        External        getOS
c -------------------------------------------------------

        path = ' '
        wrk  = ' '

c Request full specification of EXE-file:
c ATTN: with GNU Fortran under Windows arg(0) only shows the executable
c path when it was given in the command line -- see this discussion:
c https://stackoverflow.com/questions/4021763/how-do-i-get-my-executables-location-using-fortran
        call get_arg_wrapper (0,wrk)
        l = Len_Trim (wrk)
        if (l.eq.0)     return
c       write (*,*) 'defpath: get_arg_wrapper returned arg(0)=',wrk(1:l)

c Separate path from the full EXE-file specification:
        if (getOS() .eq. 'windows') Then !--------+
          do i=l,1,-1  !========================+ |
            if (wrk(i:i).eq.'\' .OR.           !| |char(92)='\'
     *          wrk(i:i).eq.'/' .OR.           !| |char(47)='/'
     *          wrk(i:i).eq.':') goto 2        !| |char(58)=':'
          enddo !===============================+ |
        else !------------------------------------+
          do i=l,1,-1  !========================+ |
            if (wrk(i:i).eq.'/') goto 2        !| |char(47)='/'
          enddo !===============================+ |
        endif !-----------------------------------+
c       write (*,*) 'defpath: no slashed found in the path=',wrk(1:l)
        return

  2     continue
        l = len (path)
c If path fits the buffer:
        if (i.le.l) Then !-----------------------------+
          path=wrk(1:i)                               !|
          l = i                                       !|
          do i=1,l !================================+  |
            if (path(i:i).eq.'\') path(i:i) = '/'  !|  |
          enddo !===================================+  |
        endif !----------------------------------------+

        return
        end
