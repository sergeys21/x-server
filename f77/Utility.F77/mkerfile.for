        Subroutine Make_Err_Filename (InpFile)
c -------------------------------------------------------
c This subroutine forms the ERR-file path from the
c name of INP-file for gid_sl97, ter_SK97, trds_97 and gids_97
c -------------------------------------------------------
        Integer         i, j, l
        Character       InpFile*(*)

        Character       getOS*8
        External        getOS

c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

        l = Len_Trim(InpFile)
        if (l.eq.0) Return
        do j=1,l  !===================================+
          if (InpFile(j:j).eq.'\') InpFile(j:j)='/'  !|
        enddo !=======================================+

c Delete the the startup ERR-file, if present
c (the startup ERR-file has the program name)
        Call DeleteErrorFile ()

c Remove path from the INP-file specification::
        if (getOS() .eq. 'windows') Then !--------+
          do j=l,1,-1  !=======================+  |
            if (InpFile(j:j).eq.'\'  .or.     !|  |
     *          InpFile(j:j).eq.'/'  .or.     !|  |
     *          InpFile(j:j).eq.':') goto 1 !--+--+--+
          enddo  !=============================+  |  |
        else !------------------------------------+  v
          do j=l,1,-1  !=======================+  |  |
            if (InpFile(j:j).eq.'/') goto 1 !--+--+--+
          enddo  !=============================+  |  |
        endif !-----------------------------------+  v
        j = 0                                       !|
  1     continue  !<---------------------------------+
        j = j+1
        ErrorFile = InpFile(j:l)
        l = l-j+1
c Replace extension by '.err':
        i = Max0(l-3,1)
        do j=l,i,-1  !========================+
          if (ErrorFile(j:j).eq.'.') goto 2 !-+--+
        enddo  !==============================+  |
        j = l+1                                 !|
  2     continue  !<-----------------------------+
        l = Len(ErrorFile)
        ErrorFile(j:l) = '.err'
c Delete the previous version of the ERR-file, if present:
        Call DeleteErrorFile ()
        return
        end
