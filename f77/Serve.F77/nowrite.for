        integer function nowrite ()
c--------------------------------------------------------
c This function request Retry/Cancel when a device or file
c is not ready for writing data. Possible returned values:
c   nowrite=0 - Retry,
c   nowrite=1 - Cancel
c--------------------------------------------------------
        integer         iscan, iasci
c
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------
        if (modebat.ne.0) then    !--+ bat-mode
          nowrite = 1               !| Cancel
          return                    !|
        endif !----------------------+

        write   (*,1)
  1     format  (
     *  'ERROR: Cannot write file: <ENTER>=Rewrite <ESC>=Cancel')
  2     continue
        call    keyin   (1,iscan,iasci)
c
        if (iasci.eq.13)     then !--+ <ENTER>
           nowrite = 0              !|
           return                   !|
        elseif (iasci.eq.27) then !--+ <ESC>
           nowrite = 1              !|
           return                   !|
        else !-----------------------+
           goto 2                   !|
        endif !----------------------+
        end
