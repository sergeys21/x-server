        Subroutine      OpenFile (filename,     !input: filename, cannot be blank
     *                            lun,          !input: file lun
     *                            rorw,         !input: read, write, readwrite
     *                            oldnew,       !input: old, new, unknown, append (append is g77 extension)
     *                            io_status,    !output:
     *                            *)            !output: label in calling program where to go on opening error
c--------------------------------------------------------
c Opens file: a wrapper to mitigate difference between Fortran compilers
c--------------------------------------------------------
        Integer   lun, io_status, lname,
     *            lrw, lonu, lshare,
     *            append
        Character filename*(*),                 !requsted file name
     *            rorw*(*),                     !READ/WRITE/READWRITE (same as "mode")
     *            oldnew*(*),                   !OLD/NEW/UNKNOWN
     *            rorw_local*12,                !local copy of rorw
     *            onu*8,                        !old, new, unknown
     *            share*8                       !denywr or denynone

        io_status = 0
        lname  = Len_Trim(filename)
        if (lname .eq. 0) goto 1                !no name specified

        rorw_local = rorw                       !make a local copy because rorw may be a constant
        lrw = Len_Trim(rorw_local)
        if (lrw .eq. 0) Then !---------+
           rorw_local = 'READWRITE'   !|        !default when not specified
           lrw = Len_Trim(rorw_local) !|
        endif !------------------------+
        Call CaseLat(rorw_local(1:lrw),0)       !convert to upper case

        onu = oldnew                            !make a local copy because we shall be changing case

        if     (rorw_local(1:lrw) .eq. 'READ') Then !------+
           share = 'DENYNONE'                             !|
           onu   = 'OLD'                                  !|
        elseif (rorw_local(1:lrw) .eq. 'WRITE') Then !-----+
           share = 'DENYWR'                               !|
        elseif (rorw_local(1:lrw) .eq. 'READWRITE') Then !-+
           share = 'DENYWR'                               !|
        else  !--------------------------------------------+
           share = 'DENYWR'                               !|
           rorw_local = 'READWRITE'                       !|
           lrw = Len_Trim(rorw_local)                     !|
        endif !--------------------------------------------+
        lshare = Len_Trim(share)
        lonu   = Len_Trim(onu)
        if (lonu .eq. 0) Then !-----+
           onu  = 'UNKNOWN'        !|           !default when not specified
           lonu = Len_Trim(onu)    !|
        endif  !--------------------+
        Call CaseLat(onu(1:lonu),0)             !convert to upper case
        append = 0
        if (onu(1:lonu) .eq. 'APPEND' .AND.
     *      rorw_local(1:lrw) .ne. 'READ') Then !-+
           append = 1                            !|
        else !------------------------------------+
           append = 0                            !|
        endif !-----------------------------------+
           
        if (onu(1:lonu) .ne. 'OLD' .AND.
     *      onu(1:lonu) .ne. 'NEW' .AND.
     *      onu(1:lonu) .ne. 'UNKNOWN') Then !------+
           onu  = 'UNKNOWN'                        !|
           lonu = Len_Trim(onu)                    !|
        endif  !------------------------------------+

        if (append .eq.1 ) Then !-----------+
           Open (unit=lun,                 !|
     *           file=filename(1:lname),   !|
     *           status=onu(1:lonu),       !|
c    *           mode=rorw_local(1:lrw),   !|   !Compaq/MS fortran only
     *           action=rorw_local(1:lrw), !|   !Compaq/GNU fortran only
c    *           share=share(1:lshare),    !|   !Compaq/MS fortran only
     *           position='APPEND',        !|   !this is g77 extension, better use scrolling to EOF
     *           iostat=io_status,         !|
     *           err=1)                    !|
        else !------------------------------+
           Open (unit=lun,                 !|
     *           file=filename(1:lname),   !|
     *           status=onu(1:lonu),       !|
c    *           mode=rorw_local(1:lrw),   !|   !Compaq/MS fortran only
     *           action=rorw_local(1:lrw), !|   !Compaq/GNU fortran only
c    *           share=share(1:lshare),    !|   !Compaq/MS fortran only
     *           iostat=io_status,         !|
     *           err=1)                    !|
        endif !-----------------------------+
        return
  1     continue
        return 1
        end

