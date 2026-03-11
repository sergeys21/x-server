c           Here are the subroutines for X0h_List:
c =======================================================
c
        Subroutine BasList (filename,
     *                      list, nlist, klist, ledb,
     *                      icoord, iX0h, iFull, iCryst, ifail)
c -------------------------------------------------------
c This subroutine compiles the X0H database catalogs
c
c PARAMETERS:
c ==========
c filename  - a name of one of the X0h database files           /input/
c list(nlist) - array of found codes from
c             'coord.x0h', 'atom.x0h', 'wave.x0h', 'group.x0h'  /output/
c
c klist     - the number of found codes.                        /output/
c
c ledb      - the length of codes (20/4/6)                      /output/
c
c icoord=1  - file=coord.x0h                                    /input/
c icoord=0  - file=atom.x0h or file=wave.x0h
c icoord=2  - file=atom.x0h (for compatibility with  subroutine BasList, same as icoord=0)
c icoord=3  - file=wave.x0h (for compatibility with  subroutine BasList, same as icoord=0)
c icoord=4  - file=group.x0h (for compatibility with subroutine BasList, same as icoord=0)
c
c
c iX0h=1    - only X0h-type structures          /coord.x0h/     /input/
c iX0h=-1   - only non-X0h structures           /coord.x0h/     /input/
c iX0h=0    - no control                        /coord.x0h/     /input/
c
c iFull=1   - only structures with full data    /coord.x0h/     /input/
c iFull=-1  - only incomplete structures        /coord.x0h/     /input/
c iFull=0   - no control                        /coord.x0h/     /input/
c
c iCryst=2  - only CUBIC CRYSTALs               /coord.x0h/     /input/ -- new 2003/06/26
c iCryst=1  - only CRYSTALs                     /coord.x0h/     /input/
c iCryst=-1 - only AMORPHOUS structures         /coord.x0h/     /input/
c iCryst=0  - no control                        /coord.x0h/     /input/
c
c ifail     - failure code ('0' corresponds to normal return)   /output/
c
c -------------------------------------------------------
c The following combinations make sense:
c   iX0h   iFull iCryst/Amorph.
c    1       1       2,1        Ok       /full crystals/
c    1       1       -1        test!    /all amorphous must be full/
c    1      -1       2,1       Ok       /uncomplete crystals/
c    1      -1       -1        test!    /uncomplete amorphous must not exist/
c   -1       1       2,1       no sence /not X0h record/
c   -1       1       -1        no sence /not X0h record/
c   -1      -1       2,1       no sence /not X0h record/
c   -1      -1       -1        no sence /not X0h record/
c -------------------------------------------------------
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
c       Integer         kcompMax, kolMax
c       Parameter      (kcompMax = 10,
c    *                  kolMax   = 32)
        Real            prcn, rpar(7)
        Integer         nlist, klist, ledb, icoord,
     *                  iFull, iX0h, iCryst,
     *                  jFull, jX0h, jCryst,
     *                  ifail, line,
     *                  n_latt(7) /1,2,2,2,3,4,6/,
     *                  lwrk, jwrk, lcod,
     *                  kcomp, isym, kol,
     *                  iControl, lun, io_status,
     *                  nFileData, iwait,
     *                  i, j, k, l, m, mc
        Character       filename*(*),
     *                  list(nlist)*(*),
     *                  code*20, name*4

        Integer         n_Column
        External        n_Column

        Character       getOS*8
        External        getOS

        Logical*4       FileExist
        External        FileExist

        Integer         lx0hpat
        Character       x0hpath*80
        Common  /x0hpa/ x0hpath, lx0hpat

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
          progname = 'x0h/BasList'        !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
c Prevents GNU fortran -Wextra warnings:
        i = kcompMax
        i = kolMax
        i = ncodMax
        i = natmMax
        i = nwavMax
        i = nedgMax
        i = nabcMax
        i = nedgTb
        i = maxX0hWarnLines
c -------------------------------------------------------
        ifail = 0

        if (iFull.lt.-1  .or. iFull.gt.1)  stop 'BasList: wrong iFull'
        if (iX0h.lt.-1   .or. iX0h.gt.1)   stop 'BasList: wrong iX0h'
        if (iCryst.lt.-1 .or. iCryst.gt.2) stop 'BasList: wrong iCryst'
        jFull  = iFull
        jX0h   = iX0h
        jCryst = iCryst

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
        if (icoord.ne.1)        Then !------------------------------+
c These controls are relevant to coord.x0h only:                    |
          jX0h = 0                                                 !|
          jFull = 0                                                !|
          jCryst = 0                                               !|
c         if (iX0h.ne.0 .OR. iFull.ne.0 .OR. iCryst.ne.0) goto 200 !|
        endif !-----------------------------------------------------+

c If control is set on Crystals or Full structures, then we must
c select X0h-type records only:
        if (jFull.ne.0 .OR. jCryst.ne.0) Then !---------------------+
          if (jX0h.eq.-1)                                 goto 200 !|
          jX0h = 1                                                 !|
        endif !-----------------------------------------------------+

        if (jX0h.ne.0 .OR. jFull.ne.0 .OR. jCryst.ne.0) Then !-+
          iControl = 1                                        !|
        else !-------------------------------------------------|
          iControl = 0                                        !|
        endif !------------------------------------------------+

        Call OpenFile(x0hpath(1:lwrk),lun,
     *                'read','old',io_status,*99)

        jwrk  = ledb+1
        klist = 0
        line  = 0

  1     continue !<---------------------------------------------+
        line = line+1                                          !|
        read (lun,2,err=20,end=3) txt(20)(1:jwrk)              !|
  2     format  (a)                                            !^
c Look for record name:                                         |
        if (txt(20)(1:1).ne.'#') goto 1 !--------->-------------+
        mc = Index(txt(20),';')                                !|
        if (mc.gt.0) txt(20)(mc:80)=' '                        !|
        Call TabExpan(txt(20),1)                               !|
        code = txt(20)(2:jwrk)                                 !|
        lcod = Max(Len_Trim(code),1)                           !|
                                                               !|
c Controls for coord.x0h:                                       ^
        if (icoord.eq.1 .AND. icontrol.eq.1) Then !---------+   |
                                                           !|   |
c 1. Get the number of elements in the structure:           |   |
  4       continue       !<------------------------+        |   |
          line = line+1                           !|        |   |
          read (lun,2,err=40,end=50) txt(20) !--->-+--------+---+ERROR-----+
c Skip empty and comment lines:                    |        |   |
          if (Len_Trim(txt(20)).eq.0)   goto 4 !->-|        |   |
          if (txt(20)(1:1).eq.';')      goto 4 !->-+        |   |
c Strip end-of-line comments:                               |   |
          mc = Index(txt(20),';')                          !|   |
          if (mc.gt.0) txt(20)(mc:80)=' '                  !|   |
c kcomp,isym,[rho] -- the first two parameters are          |   |
c integers, but we read them as reals:                      |   |
          Call TabReplace(txt(20),1)                       !|   |
          Call RdReal  (rpar,3,txt(20),ifail)              !|   |
          if (ifail.ne.0)               goto 40 !-----------+---+ERROR-----+
          kcomp = int(rpar(1)+0.1)                         !|   |
          isym  = int(rpar(2)+0.1)                         !|   |
c         rho   = rpar(3)    !not needed here              !|   |
          if (isym.lt.0 .or. isym.gt.8) goto 141 !----------+---+ERROR-----+
          if (isym.eq.0)     isym=1                        !|   |
c When filter set for cubic cystals only                    |   |
c (added 2003/06/26):                                       |   |
          if (isym.ne.1 .AND. jCryst.eq.2) goto 1 !-skip----+   |
          if (kcomp.gt.kcompMax)        goto 61 !-----------+---+ERROR-----+
c Is this structure a 'X0h' structure?                      |   |
          if (kcomp.eq.0 .AND. isym.ne.8) Then !---+        |   |
c We are here if it is non-X0h structure:          |        |   |
            if (jX0h.eq.1)  goto 1  !-------skip---+--------+   |
            if (jX0h.eq.-1) goto 9  !-------add----+--------+---+--------------+
          else  !----------------------------------|        |   |              |
            if (jX0h.eq.-1) goto 1  !-------skip---+--------+   |              |
          endif !----------------------------------+        |   |              |
c We are here if it is a X0h structure:                     |   |              |
c Is this structure a crystal?                              |   |              v
          if (isym.eq.8) Then  !-------------------+        |   |              |
            if (jCryst.ge.1)  goto 1  !-----skip---+--------+   |              |
          else !-----------------------------------|        |   |              |
            if (jCryst.eq.-1) goto 1 !------skip---+--------+   |              |
          endif  !---------------------------------+        |   |              |
c Now check coord.x0h records for 'completeness':           |   |              |
         if (isym.ne.8) Then !---------------------------+  |   |              |
c Read lattice constants:                                |  |   |              |
  5         continue      !<--------------------------+  |  |   |              |
            line = line+1                            !|  |  |   |              |
            read (lun,2,err=40,end=50) txt(20) !-->---+--+--+---+ERROR---+     |
            if (txt(20)(1:1).eq.';')    goto 5 !-->---|  |  |   |              |
c cccccc    if (Len_Trim(txt(20)).eq.0) goto 5 !-->---+  |  |   ^              |
c Strip end-of-line comments:                            |  |   |              |
            mc = Index(txt(20),';')                     !|  |   |              |
            if (mc.gt.0) txt(20)(mc:80)=' '             !|  |   |              |
            Call TabReplace(txt(20),1)                  !|  |   |              |
            nFileData = n_Column(txt(20))               !|  |   |              v
            if (isym.eq.1)  Then !---------------------+ |  |   |              |
              if (nFileData.gt.n_latt(1)+1)  goto 143 !| |  |   |              |may be Poisson
              if (nFileData.lt.1)            goto 143 !| |  |   |              |
c             nFileData = 1                           !| |  |   |              |no Poisson needed
            else !-------------------------------------| |  |   |              |
              if (nFileData.ne.n_latt(isym)) goto 143 !| |  |   ^              |never Poisson!
            endif !------------------------------------+ |  |   |              v
            Call RdReal(rpar,n_latt(isym),txt(20),ifail)!|  |   |              |
            if (ifail.ne.0) goto 40 !--------------------+--+---+ERROR---+     |
            do j=1,n_latt(isym) !===================+    |  |   |              |
              if (rpar(j).lt.0.) goto 40 !----------+----+--+---+ERROR---+     |
              if (abs(rpar(j)).lt.1.E-30) Then !-+  |    |  |   ^              |
                if (jFull.eq.1)  goto 1 !--------+--+skip+ -+  -|              |
                if (jFull.eq.-1) goto 9 !--------+--+add-+--+---+--------------|
              endif !----------------------------+  |    |  |   |              |
            enddo  !================================+    |  |   |              v
          endif  !---------------------------------------+  |   |              v
                                                           !|   |              |
c Read structure components information:                    |   ^              |
          do      k=1,kcomp  !===========================+  |   |              |
c Read names of atoms:                                   |  |   |              |
  6         line = line+1   !<--------------------+      |  |   |              |
            read (lun,2,err=40,end=50) txt(20)   !^  ->--+--+---+ERROR---+     |
            if (txt(20)(1:1).eq.';')    goto 6 !->|      |  |   |              v
c cccccc    if (Len_Trim(txt(20)).eq.0) goto 6 !->+      |  |   |              |
c Strip end-of-line comments:                            |  |   |              |
            mc = Index(txt(20),';')                     !|  |   |              |
            if (mc.gt.0) txt(20)(mc:80)=' '             !|  |   |              |
            Call TabExpan(txt(20),1)                    !|  |   |              |
            name = txt(20)(1:4)                         !|  |   ^              |
            l = Len_Trim(name)                          !|  |   |              |
            if (l.eq.0) Then !--------------+            |  |   |              |
              if (jFull.eq.1)  goto 1 !-----+------skip--+ -+  -|              |
              if (jFull.eq.-1) goto 9 !-----+------add---+--+---+--------------|
            endif !-------------------------+            |  |   |              |
c Read kol, prcn:                                        |  |   |              |
  7         continue      !<----------------------+      |  |   |              |
            line = line+1                        !|      |  |   |              |
            read (lun,2,err=40,end=50) txt(20)   !^  ->--+--+---+ERROR---+     |
            if (txt(20)(1:1).eq.';')    goto 7 !->|      |  |   |              |
c cccccc    if (Len_Trim(txt(20)).eq.0) goto 7 !->+      |  |   ^              |
c Strip end-of-line comments:                            |  |   |              |
            mc = Index(txt(20),';')                     !|  |   |              |
            if (mc.gt.0) txt(20)(mc:80)=' '             !|  |   |              |
c kol,[prcn] -- the first parameter is  integer,         |  |   |              |
c but we read it as real:                                |  |   |              |
            Call TabReplace(txt(20),1)                  !|  |   |              |
            Call  RdReal (rpar,2,txt(20),ifail)         !|  |   ^              |
            if (ifail.ne.0)             goto 40 !--------+--+---+ERROR-----+   |
            kol  = int(rpar(1)+0.1)                     !|  |   |              |
            prcn = rpar(2)                              !|  |   |              |
            if (kol.lt.1)      goto 51 !-----------------+--+---+ERROR---+     v
            if (kol.gt.kolMax .AND. isym.ne.8) goto 51 !-+--+---+ERROR---+     |
            if (prcn.lt.0.  .OR.                        !|  |   |              |
     *          prcn.gt.100.) Then !------+              |  |   |              |
              if (jFull.eq.1)  goto 1 !---+--------skip--+ -+  -|              |
              if (jFull.eq.-1) goto 9 !---+--------add---+--+---+--------------|
            endif  !----------------------+              |  |   |              |
c Read atomic coordinates:                               |  |   |              |
            if (isym.ne.8) Then !----------------------+ |  |   |              |
              do m=1,kol !============================+| |  |   |              |
  8             continue        !<------------------+ || |  |   |              |
                line = line+1                      !| || |  |   |              |
                read (lun,2,err=40,end=50) txt(20) !^ || |--+---+ERROR---+     |
                if (txt(20)(1:1).eq.';')    goto 8 !| || |  |   |              |
c Strip end-of-line comments:                       ^ || |  |   |              |
                mc = Index(txt(20),';')            !| || |  |   |              |
                if (mc.gt.0) txt(20)(mc:80)=' '    !| || |  |   |              |
                Call TabReplace(txt(20),1)         !| || |  |   |              |
c We may have x,y,z as empty string:                ^ || |  |   |              |
c cccccc        if (Len_Trim(txt(20)).eq.0) goto 8 !+ || |  |   |              |
c x,y,z:                                              || |  |   |              |
                Call RdReal (rpar,3,txt(20),ifail)   !|| |  |   |              |
                if (ifail.ne.0) goto 40 !-------------++-+--+---+ERROR---+     |
              enddo !=================================+| |  |   |              |
            endif !------------------------------------+ |  |   |              |
          enddo !========================================+  |   |              |
        endif !------------------------------------------- -+   ^              |
                                                               !|              |
        if (jFull.eq.-1) goto 1 !------------------skip---------+------->------|
c Insert into the catalog:                                      |              v
  9     continue         !<------------<----------------<-------+-------<------+
        klist = klist+1                                        !|
        if (klist.gt.nlist) goto 10 !-------------------------+-+ERROR---+
        list(klist) = code(1:lcod)                             !|        v
        goto 1 !-------------->--------------------->-----------+
  3     continue !<-----[EOF]
        Close   (unit=lun)

c Sort the database index:
        if (Index(filename,'group.x0h').eq.0) Then !-----+
           Call  Sort_codes (list,klist)                 !+
        endif !-------------------------------------------+

  999   continue
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1) Then !------------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif !----------------------------+
        return
c ----------------------------------------------------
c    ERRORS processing:
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
  53    format  (' "',a,'"'//' - read error',
     *  ' at line=',i4)
        goto 997
c ........................
  40    continue
        ifail = 40
        write   (txt,45)        x0hpath(1:lwrk), code(1:lcod), line
  45    format  ('"',a,'"'//' - read error',
     *  ' of structure <',a,'>, line=',i4)
        Close   (unit=lun)
        goto 997
c ........................
  50    continue
        ifail = 50
        write   (txt,55)        x0hpath(1:lwrk), code(1:lcod)
  55    format  ('"',a,'"'//' - uncomplete data',
     *  ' for structure <',a,'>')
        goto 997
c ........................
  10    continue
        ifail = 10
        write   (txt,11)      x0hpath(1:lwrk), nlist
  11    format  (' "',a,'"'//' - too many elements.',
     *  ' Buffer size=',i5,' exceeded')
c       goto 997
        klist = nlist
        list(nlist) = '***'
        goto 998
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
  52    format  ('"',a,'"'//a,' - number of coords for',
     *  ' component ',i2,' not in the range [1-',i2,']')
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
  144   format  ('"',a,'"'//a,' - wrong number of lattice',
     *  ' constants (',i1,')')
        goto 997
c ........................
  200   continue
        ifail = 200
        write   (txt,201)       iX0h, iFull, iCryst
  201   format  (
     *  'Conflicting material selection filters specified,'/
     *  'or the filters applied to DB other than coord.x0h.'/
     *  'X0h_filter=',i2,' Full_filter=',i2,' Crystal_filter=',i2)
        goto 998
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
        if (ifail .eq. 10) Then !--+
           ifail = 0              !|
           goto 3                 !|
        else !---------------------+
           goto 999               !|
        endif !--------------------+
c ........................
        end
c
c=======================================================
c
        Subroutine Sort_codes (name,n)

        Integer         n, ln, i, j, kn
        Character       name(n)*(*),
     *                  ntmp*32, nam1*32, nam2*32
c--------------------------------------------------------
c  Sort Names in ascending order with ignoring case
c--------------------------------------------------------
        ln = Len(ntmp)
        kn = Len(name(1))
        if (kn.gt.ln)   Stop 'Sort_Codes: long name'
        do      j = 2,n  !================================+
          ntmp(1:kn) = name(j)                           !|
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
          enddo  !=================================+  V   |
          i = 0                                      !|   |
  3       continue   !<-------------------------------+   |
          name(i+1) = ntmp(1:kn)                         !|
        enddo  !==========================================+
        return
        end

