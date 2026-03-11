c  Here are the subroutines to produce *EXPANDED* X0h_List:
c =======================================================
c
        Subroutine BasList_Ex (filename,
     *                         list, desc, nlist, klist, ledb,
     *                         icoord, iCryst, ifail)
c -------------------------------------------------------
c This subroutine makes catalogs of the X0H database files
c
c PARAMETERS:
c ==========
c filename  - a name of one of the X0h database files           /input/
c list(nlist) - array of found codes from
c             'coord.x0h', 'atom.x0h', or 'wave.x0h'            /output/
c
c desc(nlist) - array of descriptions for found codes           /output/
c              [this is the main difference from subr. BasList]
c
c klist     - the number of found codes.                        /output/
c
c ledb      - the length of codes (20/4/6)                      /output/
c
c icoord=1  - file=coord.x0h                                    /input/
c icoord=2  - file=atom.x0h (!-different from Subroutine BasList)
c icoord=3  - file=wave.x0h (!-different from Subroutine BasList)
c
c iCryst=2  - only CUBIC CRYSTAL structures     /coord.x0h/     /input/
c iCryst=1  - only CRYSTAL structures           /coord.x0h/     /input/
c iCryst=-1 - only AMORPHOUS structures         /coord.x0h/     /input/
c iCryst=0  - no control                        /coord.x0h/     /input/
c
c ifail     - failure code ('0' corresponds to normal return)   /output/
c
c -------------------------------------------------------
c Attention, unlike Subroutine BasList, this program works with:
c
c iX0h==1    - only X0h-type structures        /coord.x0h/
c iFull==1   - only structures with full data  /coord.x0h/
c -------------------------------------------------------
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
c       Integer         kcompMax, kolMax
c       Parameter      (kcompMax = 10,
c    *                  kolMax   = 32)
        Real            prcn, rho,
     *                  atom_mass_aem,
     *                  rpar(7)
        Integer         nlist, klist, z,
     *                  ledb, icoord,
     *                  iCryst, ifail, line,
     *                  n_latt(7) /1,2,2,2,3,4,6/,
     *                  lwrk, jwrk, lcod,
     *                  kcomp, isym, kol, lun,
     *                  nFileData, iwait,
     *                  lfor, mfor, ld, na,
c    *                  ll,
     *                  ipar(4), io_status,
     *                  i, j, k, l, m, mc
        Character       filename*(*),
     *                  list(nlist)*(*),
     *                  desc(nlist)*(*),
     *                  formula*32, bbb*4,
     *                  code*20, name*4, str*80,
     *                  smtry(8)*15 /'Cubic',
     *                               'Hexa/Trigonal',
     *                               'Tetragonal',
     *                               'Rhombohedral',
     *                               'Orthorhombic',
     *                               'Monoclinic',
     *                               'Triclinic',
     *                               '*Amorphous*'/

        Integer         n_Column, nonASCII
        External        n_Column, nonASCII

        Character       getOS*8
        External        getOS

        Logical*4       FileExist
        External        FileExist

        Integer         lx0hpat
        Character       x0hpath*80
        Common  /x0hpa/ x0hpath, lx0hpat

        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

        if (istackrezv.lt.10) Then !-------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'x0h/BasList_Ex'     !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif !----------------------------+
c -------------------------------------------------------
        lfor = 0                        !keep GNU Fortran happy
        mfor = 0                        !keep GNU Fortran happy
        i = kcompMax                    !to suppress unused paramters warnings
        i = kolMax                      !to suppress unused paramters warnings
        i = ncodMax                     !to suppress unused paramters warnings
        i = natmMax                     !to suppress unused paramters warnings
        i = nwavMax                     !to suppress unused paramters warnings
        i = nedgMax                     !to suppress unused paramters warnings
        i = nabcMax                     !to suppress unused paramters warnings
        i = nedgTb                      !to suppress unused paramters warnings
        i = maxX0hWarnLines             !to suppress unused paramters warnings
c -------------------------------------------------------
        if (icoord.lt.1 .or.icoord.gt.3) stop 'BasList_Ex: wrong icoord'
        if (iCryst.lt.-1.or.iCryst.gt.2) stop 'BasList_Ex: wrong iCryst'

        ld = Len(desc(1))
        do i=1,nlist  !=====+
          list(i)  = ' '   !|
          desc(i)  = ' '   !|
        enddo  !============+
        ifail = 0

c Path to X0h database:
        if (x0hpath(1:1).eq.char(0)) Then !------------------------+
          Call Defpat (x0hpath)                                   !|
          i = Len_Trim (x0hpath)                                  !|
          x0hpath(i+1:80) = filename                              !|
          if (.NOT.FileExist(x0hpath)) Then !----------------+     |
            Call  GetEnv32 ('X0H',txt(20))                  !|     |
            i = Len_Trim (txt(20))                          !|     |
            if (i.eq.0)   goto 99 !--------------------------+--+  |
            x0hpath = txt(20)                               !|  V  |
            do j=1,i  !==================================+   |     |
              if (x0hpath(j:j).eq.'\') x0hpath(j:j)='/' !|   |     |
            enddo !======================================+   |     |
            if (x0hpath(i:i).ne.char(92) .AND.              !|     |char(92)='\'
     *          x0hpath(i:i).ne.char(47)) Then !-----+       |     |char(47)='/'
              i = i+1                               !|       |     |
              if (getOS() .eq. 'windows') Then !--+  |       |     |
c               x0hpath(i:i) = char(92)          !|  |       |     |char(92)='\'
                x0hpath(i:i) = char(47)          !|  |       |     |char(92)='\'
              else !------------------------------+  |       |     |
                x0hpath(i:i) = char(47)          !|  |       |     |char(47)='/'
              endif !-----------------------------+  |       |     |
            endif !----------------------------------+       |     |
            x0hpath(i+1:80) = filename                      !|     |
            if (.NOT.FileExist(x0hpath)) goto 99 !-----------+--+  |
          endif !--------------------------------------------+  V  |
          lx0hpat = i                                             !|
        endif !----------------------------------------------------+
        x0hpath(lx0hpat+1:80) = filename
        lwrk = Len_Trim(x0hpath)
        lun = 10
c -------------------------------------------------------
        Call OpenFile(x0hpath(1:lwrk),lun,
     *                'read','old',io_status,*99)

        jwrk  = ledb+1
        klist = 0
        line  = 0

  1     continue !<------------------------------------------+
        line = line+1                                       !|
        read (lun,2,err=20,end=3) txt(20)                   !|
  2     format  (a)                                         !^
c Look for record name:                                      |
        if (txt(20)(1:1).ne.'#')        goto 1 !---->--------+
        code = txt(20)(2:jwrk)                              !|
        lcod = Len_Trim (code)                              !|
                                                            !|
c Controls for coord.x0h:                                    ^
        if (icoord.eq.1) Then !----------------------------+ |
                                                          !| |
c 1. Get the number of elements in the structure:          | |
  4       continue       !<--------------------------+     | |
          line = line+1                             !|     | |
          str = 'Reading number of elements string' !|     | |
          read (lun,2,err=40,end=50) txt(20) !--->---+-----+-+---ERROR-----+
c Skip empty and comment lines:                      |     | |
          if (Len_Trim(txt(20)).eq.0) goto 4 !-->----+     | |
          if (txt(20)(1:1).eq.';')    goto 4 !-->----+     | |
c Strip end-of-line comments:                              | |
          mc = Index(txt(20),';')                         !| |
          if (mc.gt.0) txt(20)(mc:80)=' '                 !| |
          if (nonASCII(txt(20),na).ne.0) goto 90 !---------+-+---ERROR-----+
c kcomp,isym,[rho] -- the first two parameters are         | |
c integers, but we read them as reals:                     | |
          Call TabReplace(txt(20),1)                      !| |
          Call RdReal (rpar,3,txt(20),ifail)              !| ^
          if (ifail.ne.0) Then !-----------+               | |
            write(str,'(a,i4,a)')         !|               | |
     *        'Rdreal fail =',ifail,      !|               | |
     *        ' on kcomp,isym,[rho]'      !|               | |
            goto 40 !----------------------+---------------+-+---ERROR-----+
          endif !--------------------------+               | |
          kcomp = int(rpar(1)+0.1)                        !| |
          isym  = int(rpar(2)+0.1)                        !| |
          rho   = rpar(3)                                 !| |
          if (isym.lt.0 .or. isym.gt.8) goto 141 !---------+-+---ERROR-----+
          if (isym.eq.0)     isym=1                       !| |
c When filter set for cubic cystals only                   | ^
c (added 2003/07/07):                                      | |
          if (isym.ne.1 .AND. iCryst.eq.2) goto 1 !-skip---+-+
          if (kcomp.gt.kcompMax)        goto 61 !----------+-+---ERROR-----+
c Is this structure a 'X0h' structure?                     | |
          if (kcomp.eq.0 .AND. isym.ne.8) goto 1 !--skip---+-+
c We are here if it is a X0h structure:                    | |                 |
c Is this structure a crystal?                             | ^                 v
          if (isym.eq.8) Then !--------------------+       | |                 |
            if (rho.le.0.)    goto 81 !------------+-------+-+---ERROR-----+   |
            if (iCryst.ge.1)  goto 1  !-----skip---+-------+-+                 |
            if (kcomp.eq.0) Then !------+          |       | |                 |
              formula = code           !|          |       | ^                 |
              mfor = Len_Trim(formula) !|          |       | |                 |
              goto 9  !-----------------+---add----+-------+-+-----------------+
            endif  !--------------------+          |       | |                 |
          else !-----------------------------------|       | |                 |
            if (iCryst.eq.-1) goto 1  !-----skip---+-------+-+                 |
          endif !----------------------------------+       | |                 |
c Now check the record for 'completeness':                 | |                 |
          if (isym.ne.8) Then !--------------------------+ | |                 |
c Read lattice constants:                                | | ^                 |
  5         continue      !<--------------------------+  | | |                 |
            line = line+1                            !|  | | |                 |
            str = 'Reading lattice constants string' !|  | | |                 |
            read (lun,2,err=40,end=50) txt(20) !-->---+--+-+-+-ERROR-----+     |
            if (txt(20)(1:1).eq.';')    goto 5 !-->---|  | | |                 |
c cccccc    if (Len_Trim(txt(20)).eq.0) goto 5 !-->---+  | | ^                 |
c Strip end-of-line comments:                            | | |                 |
            mc = Index(txt(20),';')                     !| | |                 |
            if (mc.gt.0) txt(20)(mc:80)=' '             !| | |                 |
            if (nonASCII(txt(20),na).ne.0) goto 90 !-----+-+-+-ERROR-----+     |
            Call TabReplace(txt(20),1)                  !| | |                 |
            nFileData = n_Column(txt(20))               !| | ^                 v
            if (isym.eq.1)  Then !---------------------+ | | |                 |
              if (nFileData.gt.n_latt(1)+1)  goto 143 !| | | |may be Poisson   |
              if (nFileData.lt.1)            goto 143 !| | | |                 |
c             nFileData = 1                           !| | | |no Poisson needed|
            else !-------------------------------------| | | |                 |
              if (nFileData.ne.n_latt(isym)) goto 143 !| | | ^never Poisson!   |
            endif !------------------------------------+ | | |                 v
            Call RdReal(rpar,n_latt(isym),txt(20),ifail)!| | |                 |
            if (ifail.ne.0) Then !---------+             | | |                 |
              write(str,'(a,i4,a,i1,a)')  !|             | | |                 |
     *            'Rdreal fail =',ifail,  !|             | | |                 |
     *            ' on ',n_latt(isym),    !|             | | ^                 |
     *            ' lattice constants'    !|             | | |                 |
              goto 40 !--------------------+-------------+-+-+-ERROR-----+     |
            endif !------------------------+             | | |                 |
            do j=1,n_latt(isym) !==================+     | | |                 |
              if (rpar(j).lt.0.) Then !----------+ |     | | |                 |
                write(str,'(a,i1)')             !| |     | | ^                 |
     *           'Negative lattice constant-',j !| |     | | |                 |
                goto 40 !------------------------+-+-----+-+-+-ERROR-----+     |
              endif !----------------------------+ |     | | |                 |
c If the structure is not full:                    |skip | | |                 |
              if (abs(rpar(j)).lt.1.E-10) goto 1 !-+->---+-+-+                 |
            enddo  !===============================+     | | |                 v
            Call Fill_Cell (rpar,n_latt(isym),isym)     !| | |                 |
          endif !----------------------------------------+ | |                 |
                                                          !| |                 |
c Read structure components information:                   | ^                 |
          formula = ' '                                   !| |                 |
          lfor = Len(formula)                             !| |                 |
          mfor = 0                                        !| |                 |
          do      k=1,kcomp  !===========================+ | ^                 |
c Read names of atoms:                                   | | |                 |
  6         continue      !<----------------------+      | | |                 |
            line = line+1                        !|      | | |                 |
            str = 'Reading atom names string'    !|      | | |                 |
            read (lun,2,err=40,end=50) txt(20)   !^  ----+-+-+-ERROR-----+     |
            if (txt(20)(1:1).eq.';')    goto 6 !->|      | | |                 v
c cccccc    if (Len_Trim(txt(20)).eq.0) goto 6 !->+      | | |                 |
c Strip end-of-line comments:                            | | ^                 |
            mc = Index(txt(20),';')                     !| | |                 |
            if (mc.gt.0) txt(20)(mc:80)=' '             !| | |                 |
            if (nonASCII(txt(20),na).ne.0) goto 90 !-----+-+-+-ERROR-----+     |
            Call TabExpan(txt(20),1)                    !| | |                 |
            name = txt(20)(1:4)                         !| | ^                 |
            l = Len_Trim(name)                          !| | |                 |
c If the structure is not full:                         !| | |                 |
            if (l.eq.0) goto 1 !----------skip---->------+-+-+                 |
            if (mfor+l.gt.lfor) goto 71 !----------------+-+-+-ERROR-----+     |
            formula(mfor+1:mfor+l) = name               !| | |                 |
            mfor = mfor+l                               !| | |                 |
c Read kol, prcn:                                        | | |                 |
  7         continue      !<----------------------+      | | ^                 |
            line = line+1                        !|      | | |                 |
            continue                             !|      | | |                 |
            str = 'Reading kol, [prcn] string'   !|      | | |                 |
            read (lun,2,err=40,end=50) txt(20)   !^      |-+-+-ERROR-----+     |
            if (txt(20)(1:1).eq.';')    goto 7 !->|      | | |                 |
c cccccc    if (Len_Trim(txt(20)).eq.0) goto 7 !->+      | | ^                 |
c Strip end-of-line comments:                            | | |                 |
            mc = Index(txt(20),';')                     !| | |                 |
            if (mc.gt.0) txt(20)(mc:80)=' '             !| | |                 |
c kol,[prcn] -- the first parameter is integer,          | | ^                 |
c but we read it as real:                                | | |                 |
            if (nonASCII(txt(20),na).ne.0) goto 90 !-----+-+-+-ERROR-----+     |
            Call TabReplace(txt(20),1)                  !| | |                 |
            Call  RdReal (rpar,2,txt(20),ifail)         !| | ^                 |
            if (ifail.ne.0) Then !--------------+        | | |                 |
              write(str,'(a,i4,a)')            !|        | | |                 |
     *             'Rdreal fail =',ifail,      !|        | | |                 |
     *             ' on kol, [prcn]'           !|        | | ^                 |
              goto 40 !-------------------------+--------+-+-+---ERROR-----+   |
            endif !-----------------------------+        | | |                 |
            kol  = int(rpar(1)+0.1)                     !| | |                 |
            prcn = rpar(2)                              !| | |                 |
            if (kol.lt.1)                      goto 51 !-+-+-+-ERROR-----+     v
            if (kol.gt.kolMax .AND. isym.ne.8) goto 51 !-+-+-+-ERROR-----+     |
            if (prcn.lt.0. .OR.                         !| | |                 |
     *          prcn.gt.100.)   goto 1 !-----skip--------+-+-+                 |
            if (kol.gt.1) Then !--------------+          | | |                 |
              Call Pakint (kol,1,bbb,l)      !|          | | |                 |
              if (mfor+l.gt.lfor) goto 71 !---+----------+-+-+-ERROR-----+     |
              formula(mfor+1:mfor+l) = bbb   !|          | | |                 |
              mfor = mfor+l                  !|          | | |                 |
            endif !---------------------------+          | | ^                 |
c Read atomic coordinates:                               | | |                 |
            if (isym.ne.8) Then !----------------------+ | | |                 |
              do m=1,kol !============================+| | | |                 |
  8             line = line+1 !<--------------------+ || | | |                 |
                str = 'Reading coords string'      !| || | | ^                 |
                read (lun,2,err=40,end=50) txt(20) !^ || +-+-+-ERROR-----+     v
                if (txt(20)(1:1).eq.';')    goto 8 !| || | | |                 |
c Strip end-of-line comments:                       ^ || | | ^                 |
                j = Index(txt(20),';')             !| || | | |                 |
                if (j.gt.0) txt(20)(j:80) = ' '    !| || | | |                 |
                if (nonASCII(txt(20),na).ne.0)     !| || | | |                 |
     *                                   goto 90 !--+-++-+-+-+-ERROR-----+     |
                Call TabReplace(txt(20),1)         !| || | | |                 |
c We may have x,y,z as empty string:                ^ || | | |                 |
c cccccc        if (Len_Trim(txt(20)).eq.0) goto 8 !+ || | | |                 |
c x,y,z:                                              || | | ^                 |
                Call RdReal (rpar,3,txt(20),ifail)   !|| | | |                 |
                if (ifail.ne.0) Then !-------+        || | | |                 |
                  write(str,'(a,i4,a)')     !|        || | | ^                 |
     *               'Rdreal fail =',ifail, !|        || | | |                 |
     *                ' of coords'          !|        || | | |                 |
                  goto 40 !------------------+--------++-+-+-+-ERROR-----+     |
                endif !----------------------+        || | | |                 |
              enddo !=================================+| | | |                 |
            endif !------------------------------------+ | | |                 |
          enddo !========================================+ | |                 |
                                                          !| |                 |
c Controls for atom.x0h:                                   | ^                 |
        elseif (icoord.eq.2) Then !------------------------+ |                 |
                                                          !| |                 |
c Read Ifactor, Idisper, Isigma, Idebye                    | |                 |
  13      continue      !<-------------------------+       | |                 |
          line = line+1                           !|       | |                 |
          str = 'Reading Ifactor,Idisper,'//      !|       | |                 |
     *          'Isigma,Idebye string'            !|       | |                 |
          read (lun,2,err=40,end=50) txt(20)      !^  -----+-+-ERROR-----+     v
          if (txt(20)(1:1).eq.';')    goto 13 !----+       | |                 |
c Strip end-of-line comments:                      ^       | |                 |
          j = Index(txt(20),';')                  !|       | |                 |
          if (j.gt.0) txt(20)(j:80) = ' '         !|       | |                 |
          if (nonASCII(txt(20),na).ne.0) goto 90 !-+-------+-+-ERROR-----+     |
          Call TabReplace(txt(20),1)              !|       | |                 |
c cccccc  if (Len_Trim(txt(20)).eq.0) goto 13 !----+       | |                 |
          Call RdInt (ipar,4,txt(20),ifail)               !| |                 |
          if (ifail.ne.0) Then !------------------------+  | |                 |
            write(str,'(a,i4,a)')                      !|  | |                 |
     *            'RdInt fail =',ifail,                !|  | |                 |
     *            ' on Ifactor,Idisper,Isigma,Idebye'  !|  | |                 |
            goto 40 !-----------------------------------+--+-+-ERROR-----+     |
          endif !---------------------------------------+  | |                 |
                                                          !| |                 |
c Read Z, atomic_weight, Density(gm/cm^3)                  | |                 |
  14      continue      !<-------------------------+       | |                 |
          line = line+1                           !|       | |                 |
          str = 'Reading Z,a_weight,density'      !|       | |                 |
          read (lun,2,err=40,end=50) txt(20)      !^       +-+-ERROR-----+     v
          if (txt(20)(1:1).eq.';')    goto 14 !----+       | |                 |
c Strip end-of-line comments:                      ^       | |                 |
          j = Index(txt(20),';')                  !|       | |                 |
          if (j.gt.0) txt(20)(j:80) = ' '         !|       | |                 |
          if (nonASCII(txt(20),na).ne.0) goto 90 !-+-------+-+-ERROR-----+     |
          Call TabReplace(txt(20),1)              !|       | |                 |
c cccccc  if (Len_Trim(txt(20)).eq.0) goto 14 !----+       | |                 |
          Call RdReal (rpar,3,txt(20),ifail)              !| |                 |
          if (ifail.ne.0) Then !------------+              | |                 |
            write(str,'(a,i4,a)')          !|              | |                 |
     *        'Rdreal fail = ',ifail,      !|              | |                 |
     *         ' on Z,a_weight,density'    !|              | |                 |
            goto 40 !-----------------------+--------------+-+-ERROR-----+     |
          endif !---------------------------+              | |                 |
          z             = int(rpar(1)+0.1)   !rounding     | |                 |
          atom_mass_aem = rpar(2)                         !| |                 |
          rho           = rpar(3)                         !| |                 |
          isym          = 8                               !| |                 |
          formula       = code                            !| |                 |
          mfor          = Len_Trim(formula)               !| |                 |
                                                          !| |                 |
c Controls for wave.x0h:                                   | ^                 |
        elseif (icoord.eq.3) Then !------------------------+ |                 |
                                                          !| |                 |
c NONE!                                                   !| |                 |
                                                          !| |                 |
        endif !--------------------------------------------+ |                 |
                                                            !^                 |
c Insert the code into catalog:                              |                 |
  9     continue         !<------------<----------------<----+----------<------+
        klist = klist+1                                     !|
        if (klist.gt.nlist) goto 10 !------------------------+--+error
        list(klist) = code(1:lcod)                          !|  v
                                                            !|
c Coord.x0h:                                                 ^
        if (icoord.eq.1) Then !----------------------------+ |
                                                          !| |
          if (isym.lt.8) Then  !-------------------------+ | |
            write (desc(klist),31) smtry(isym),         !| | |
     *                            (a(i),i=1,6)          !| | |
  31        format(1x,a,' cell=[',3(f6.3,1x),2(f6.2,1x),
     *                                        f6.2,']') !| | |
            l = Len_Trim(desc(klist))+1                 !| | |
          else !-----------------------------------------| | |
            write (desc(klist),32) smtry(isym)          !| | |
  32        format(1x,a,' rho=')                        !| | ^
            l = Len_Trim(desc(klist))                   !| | |
            Call PakReal(rho,1,desc(klist)(l+1:ld),i)   !| | |
            l = l + 9                                   !| | |
          endif !----------------------------------------+ | |
          desc(klist)(l+1:ld) = ' /'//formula(1:mfor)//'/'!| |
                                                          !| |
c Controls for atom.x0h:                                   | ^
        elseif (icoord.eq.2) Then !------------------------+ |
                                                          !| |
          write (desc(klist),32) smtry(isym)              !| |
          l = Len_Trim(desc(klist))                       !| ^
          Call PakReal(rho,1,desc(klist)(l+1:ld),i)       !| |
          l = l + 9                                       !| |
          desc(klist)(l+1:ld) = ' /'//formula(1:mfor)//'/'!| |
                                                          !| |
c Controls for wave.x0h:                                   | ^
        elseif (icoord.eq.3) Then !------------------------| |
c Wave length:                                            !| |
        desc(klist) = ' wavelength='//txt(20)(jwrk+2:80)  !| |
                                                          !| |
        endif !--------------------------------------------+ |
                                                            !|
        goto 1 !-------------->--------------------->--------+
  3     continue !<-----[EOF]
        Close   (unit=lun)

c Sort the database index:
c       ll = Len(list(1))
        Call  Sort_codes_Ex (list, desc, klist)         !, ll, ld)

  999   continue
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1) Then !------------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif !----------------------------+
        return
c ----------------------------------------------------
c    ERRORS PROCESSING:
c ........................
  99    continue
        ifail = 99
        lwrk = Len_Trim(x0hpath)
        write   (txt,12)      x0hpath(1:lwrk)
  12    format  (' "',a,'"'//' - open error')
        goto 998
c ........................
  20    continue
        ifail = 20
        write   (txt,53)      x0hpath(1:lwrk),line
  53    format  (' "',a,'"'//' - read error at line=',i4)
        goto 997
c ........................
  40    continue
        ifail = 40
        i = max(Len_Trim(str),1)
        write   (txt,45)        x0hpath(1:lwrk), code(1:lcod),
     *                          line, str(1:i)
  45    format  ('"',a,'"'/
     *           ' - read error on structure <',a,'>, line=',i4/
     *           ' - from: ',a)
        Close   (unit=lun)
        goto 997
c ........................
  50    continue
        ifail = 50
        write   (txt,55)        x0hpath(1:lwrk), code(1:lcod)
  55    format  ('"',a,'"'//' - uncomplete data for structure <',a,'>')
        goto 997
c ........................
  10    continue
        ifail = 10
        write   (txt,11)      x0hpath(1:lwrk)
  11    format  (' "',a,'"'//' - too many elements.',1x,
     *  'Buffer size exceeded')
        goto 997
c ........................
  61    continue
        ifail = 61
        write   (txt,62)        x0hpath(1:lwrk), code(1:lcod), kcompMax
  62    format  ('"',a,'"'//a,' - number of components >',i3)
        goto 997
c ........................
  51    continue
        ifail = 51
        write   (txt,52)        x0hpath(1:lwrk), code(1:lcod),
     *                          k, kolMax
  52    format  ('"',a,'"'//a,' - number of coords for',1x,
     *  'component ',i2,' not in the range [1-',i2,']')
        goto 997
c ........................
  71    continue
        ifail = 71
        write   (txt,72)        x0hpath(1:lwrk), code(1:lcod)
  72    format  ('"',a,'"'//a,' - too short formula buffer')
        goto 997
c ........................
  81    continue
        ifail = 81
        write   (txt,82)        x0hpath(1:lwrk), code(1:lcod), rho
  82    format  ('"',a,'"'//a,' - negative rho=',g10.3)
        goto 997
c ........................
  90    continue
        ifail = 91
        write   (txt,91)      x0hpath(1:lwrk),line,na
  91    format  (' "',a,'"'//
     *           ' - non-ASCII character at line=',i4,', pos=',i3)
        goto 997
c ........................
  141   continue
        ifail = 141
        write   (txt,142)       x0hpath(1:lwrk), code(1:lcod), isym
  142   format  ('"',a,'"'//a,' - invalid symmetry code, (',i3,')')
        goto 997
c ........................
  143   continue
        ifail = 143
        write   (txt,144)       x0hpath(1:lwrk), code(1:lcod),
     *                          nFileData
  144   format  ('"',a,'"'//a,' - wrong number of lattice',1x,
     *  'constants (',i1,')')
        goto 997
c ........................
  997   continue
        Close (unit=lun)
  998   continue
        if (modebat.eq.0) Then !---+
          iwait = 2               !|
        else !---------------------|
          iwait = 0               !|
        endif !--------------------+
        Call    Message (txt,3,iwait)
        goto 999
c ........................
        end

c
c=======================================================
c
        Subroutine Sort_codes_Ex (name, desc, n)        !, ln, ld)

        Integer         n, i, j,
     *                  kn, kd,
     *                  ln, ld
        Character       name(n)*(*),
     *                  desc(n)*(*)
c--------------------------------------------------------
c These are dynamically allocatable strings:
c DO NOT WORK WITH FORTRAN 5.1 or INIX FORTRAN:
c       Character(ln)   ntmp, nam1, nam2
c       Character(ld)   dtmp
c--------------------------------------------------------
        Character       ntmp*32, nam1*32, nam2*32
        Character       dtmp*256
        ln = Len(ntmp)
        ld = Len(dtmp)
c--------------------------------------------------------
c  Sort Names in ascending order with ignoring case
c--------------------------------------------------------
        kn = Len(name(1))
        kd = Len(desc(1))
        if (kn.gt.ln)   Stop 'Sort_Codes: long name'
        if (kd.gt.ld)   Stop 'Sort_Codes: long desciption'
        do      j = 2,n  !================================+
          ntmp(1:kn) = name(j)                           !|
          dtmp(1:kd) = desc(j)                           !|
                                                         !|
c For case-insensitive sorting:                          !|
          nam1(1:kn) = name(j)                           !|
          Call  Case (nam1(1:kn),0)                      !|
                                                         !|
          do    i = j-1,1,-1  !====================+      |
c For case-insensitive sorting:                   !|      |
            nam2(1:kn) = name(i)                  !|      |
            Call Case (nam2(1:kn),0)              !|      |
                                                  !|      |
            if (nam2(1:kn).le.nam1(1:kn)) goto 3 !-+->+   |
                                                  !|  |   |
            name(i+1) = name(i)                   !|  |   |
            desc(i+1) = desc(i)                   !|  |   |
          enddo  !=================================+  V   |
          i = 0                                      !|   |
  3       continue   !<-------------------------------+   |
          name(i+1) = ntmp(1:kn)                         !|
          desc(i+1) = dtmp(1:kd)                         !|
        enddo  !==========================================+
        return
        end

c
c=======================================================
c
        Subroutine Fill_Cell (rpar,n,isym)

        Integer n, isym, j
        Real    rpar(n)

        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        goto (91,92,93,94,95,96,97) isym !-+
c Cubic structure:                         |
c  0,1 - a(1) [, Poisson]                  |
  91    continue !<------------------------|
        a(1) = rpar(1)                    !|
        a(2) = rpar(1)                    !|
        a(3) = rpar(1)                    !|
        a(4) = 90.                        !v
        a(5) = 90.                        !|
        a(6) = 90.                        !|
        goto 98 !--------------------------+---+
c Hexagonal/trigonal structure:            |   |
c  2   - a(1),a(3)                         |   |
  92    continue !<------------------------|   |
        a(1) = rpar(1)                    !|   |
        a(2) = rpar(1)                    !|   |
        a(3) = rpar(2)                    !|   |
        a(4) = 90.                        !|   |
        a(5) = 90.                        !v   v
        a(6) = 120.                       !|   |
        goto 98 !--------------------------+---|
c Tetragonal structure:                    |   |
c  3   - a(1),a(3)                         |   |
  93    continue !<------------------------|   |
        a(1) = rpar(1)                    !|   |
        a(2) = rpar(1)                    !|   |
        a(3) = rpar(2)                    !|   |
        a(4) = 90.                        !|   |
        a(5) = 90.                        !v   v
        a(6) = 90.                        !|   |
        goto 98 !--------------------------+---|
c Trigonal (rhombohedral) structure:       |   |
c  4   - a(1),a(4)                         |   |
  94    continue !<------------------------|   |
        a(1) = rpar(1)                    !|   |
        a(2) = rpar(1)                    !|   |
        a(3) = rpar(1)                    !|   |
        a(4) = rpar(2)                    !|   |
        a(5) = rpar(2)                    !v   v
        a(6) = rpar(2)                    !|   |
        goto 98 !--------------------------+---|
c Orthorhombic structure:                  |   |
c  5   - a(1),a(2),a(3)                    |   |
  95    continue !<------------------------|   |
        a(1) = rpar(1)                    !|   |
        a(2) = rpar(2)                    !|   |
        a(3) = rpar(3)                    !|   |
        a(4) = 90.                        !|   |
        a(5) = 90.                        !v   v
        a(6) = 90.                        !|   |
        goto 98 !--------------------------+---|
c Monoclinic structure:                    |   |
c  6   - a(1),a(2),a(3),a(5)               |   |
  96    continue !<------------------------|   |
        a(1) = rpar(1)                    !|   |
        a(2) = rpar(2)                    !|   |
        a(3) = rpar(3)                    !|   |
        a(4) = 90.                        !|   |
        a(5) = rpar(4)                    !v   v
        a(6) = 90.                        !|   |
        goto 98 !--------------------------+---|
c Triclinic structure:                     |   |
c  7   - a(1),a(2),a(3),a(4),a(5),a(6)     |   |
  97    continue !<------------------------+   |
        do j=1,6  !=========+                  |
          a(j) = rpar(j)   !|                  |
        enddo  !============+                  |
        goto 98 !--------------->--------------|
  98    continue !<----------------------------+
c Determine the volume of unit cell:
c cc    lbev = 0
c cc    Call Backlat()
        return
        end
