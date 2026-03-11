c Contents:  subroutine Basini
c            subroutine Redwav
c            subroutine Redcor
c            subroutine Redlat
c            subroutine Redpar
c            subroutine Insert_Poisson
c            subroutine Insert_Edge
c            subroutine CheckChiNames
c            integer function CodeType
c =======================================================
c
        Subroutine      BasIni  (filename,list,nlist,
     *                          klist,ismlist,indx,
     *                          imode,kelem,ifail)
c -------------------------------------------------------
c Make the catalogs of the X0H database files.
c imode=0 - wave.x0h, atom.x0h, tube.x0h
c imode=1 - coord.x0h
c -------------------------------------------------------
        Integer         nlist, klist, ismlist(nlist),
     *                  indx, imode, kelem, ifail,
     *                  jfail, lun, kcomp, isym, kempt,
     *                  lwrk, jwrk, ipar(2), i, j, mc,
     *                  io_status
        Character       filename*(*), list(nlist)*(*)

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

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'x0h/BasIni'         !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+
c -------------------------------------------------------
        lun = 10

c Is the path to the files defined in the calling program?
c If not, then find out:
        if (x0hpath(1:1).eq.char(0)) Then !------------------------+
                                                                  !|
          Call Defpat (x0hpath)                                   !|
          i = Len_Trim (x0hpath)                                  !|
          lx0hpat = i                                             !|
                                                                  !|
          x0hpath(i+1:80) = filename                              !|
          if (.NOT.FileExist(x0hpath)) Then  !---------------+     |
                                                            !|     |
            Call  GetEnv32 ('X0H',txt(20))                  !|     |
            i = Len_Trim (txt(20))                          !|     |
            if (i.eq.0)   goto 99  !-------------------------+--+  |
                                                            !|  V  |
            x0hpath = txt(20)                               !|     |
            do j=1,i  !==================================+   |     |
              if (x0hpath(j:j).eq.'\') x0hpath(j:j)='/' !|   |     |
            enddo !======================================+   |     |
            if (x0hpath(i:i).ne.'\' .AND.                   !|     |
     *          x0hpath(i:i).ne.'/') Then !----------+       |     |
              i = i+1                               !|       |     |
              if (getOS() .eq. 'windows') Then !--+  |       |     |
c               x0hpath(i:i)='\'                 !|  |       |     |
                x0hpath(i:i)='/'                 !|  |       |     |
              else !------------------------------+  |       |     |
                x0hpath(i:i)='/'                 !|  |       |     |
              endif  !----------------------------+  |       |     |
            endif  !---------------------------------+       |     |
            lx0hpat = i                                     !|     |
                                                            !|     |
            x0hpath(i+1:80) = filename                      !|     |
            if (.NOT.FileExist(x0hpath)) goto 99 !-----------+--+  |
                                                            !|  V  |
          endif  !-------------------------------------------+     |
                                                                  !|
        endif  !---------------------------------------------------+

        x0hpath(lx0hpat+1:80) = filename
        lwrk = Len_Trim(x0hpath)
        jwrk = Len(list(1))+1

        Call OpenFile(x0hpath(1:lwrk),lun,'read','old',io_status,*99)
        klist = 0
        i     = 0
  1     continue   !-<------------------------------------+
        i = i+1                                          !|
        read (lun,2,err=20,end=3) txt(20)(1:jwrk)        !|
  2     format  (a)                                      !|
c Look for record name:                                   |
        if (txt(20)(1:1).ne.'#')        goto 1   !-->-----+
c Insert the code into the catalog:                       ^
        klist = klist+1                                  !|
        if (klist.gt.nlist)             goto 10  !--------+--+error
        Call TabExpan (txt(20),1)                        !|
        list(klist) = txt(20)(2:jwrk)                    !|  V
        ismlist(klist) = i                               !|
        if (imode.ne.0) Then  !----------------------+    |coord.x0h
c Get the number of elements in the structure:       |    |
  5       i = i+1  !<-----------------------------+  |    |
          read  (lun,2,err=20,end=20) txt(20)    !|  |    |
c Skip empty and comment lines:                   |  |    |
          if (txt(20)(1:1).eq.';')    goto 5 !-->-+  |    |
          mc = Index(txt(20),';')                !|  |    |
          if (mc.gt.0) txt(20)(mc:80) = ' '      !|  |    |
          Call TabReplace (txt(20),1)            !|  |    |
          if (Len_Trim(txt(20)).eq.0) goto 5 !-->-+  |    |
          Call  RdInt(ipar,2,txt(20),jfail)         !|    ^
          if (jfail.ne.0)           goto 20 !--------+----+--+error
          kcomp = ipar(1)                           !|    |
          isym  = ipar(2)                           !|    |
          if (kcomp.eq.0 .AND. isym.ne.8)           !|    |This is NOT a X0h structure
     *               ismlist(klist)=-ismlist(klist) !|    |(used for Bragg angles only)
        endif  !-------------------------------------+    |
        goto 1  !------------->--------------------->-----+
  3     continue  !<----[EOF]
        indx = 1
c Sort the database index:
        if (imode.eq.0) Then  !----------------------+
          Call  Sort_ABC (list,ismlist,klist)       !|
        else  !--------------------------------------+
c Shift all empty structures to the end:             |
          Call  Sort_Index (list,ismlist,klist)     !|
          do    kelem=klist,1,-1  !==============+   |
            if (ismlist(kelem).gt.0)  goto 6 !->-+-+ |
            ismlist(kelem)=-ismlist(kelem)      !| | |
          enddo  !===============================+ V |
          kelem = 0                               !| |
  6       continue  !<-----------------------------+ |
          kempt = klist-kelem                       !|
          if (kelem.gt.1)                           !|
     *       Call Sort_ABC (list,ismlist,kelem)     !|
          if (kempt.gt.1)                           !|
     *       Call Sort_ABC (list(kelem+1),          !|
     *                   ismlist(kelem+1),kempt)    !|
        endif  !-------------------------------------+
  998   continue
        Close   (unit=lun,err=999)
  999   continue
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c ----------------------------------------------------
c    ERRORS PROCESSING:
c ........................
  99    continue
        lwrk = Len_Trim(x0hpath)
        write   (txt,52)      x0hpath(1:lwrk), io_status
  52    format  (' "',a,'"'//' - open error, io_status=',i3)
        ifail = 3                                !number of lines
        goto 999
c ........................
  10    continue
        write   (txt,51)      x0hpath(1:lwrk)
  51    format  (' "',a,'"'//' - too many elements.',
     *  ' Buffer size exceeded')
        ifail = 3                                !number of lines
        goto 998
c ........................
  20    continue
        write   (txt,53)      x0hpath(1:lwrk),
     *                          txt(20)(1:jwrk)
  53    format  (' "',a,'"'//' - read error',
     *  ' at element [',a,']')
        ifail = 3                                !number of lines
        goto 998
c ........................
        end
c
c =======================================================
c
        Subroutine      RedWav  (wave,iinp,ifail)
c -------------------------------------------------------
c Determine X-ray wavelength either by reading "wave.x0h" by
c given characteristic line name or by requesting from terminal.
c iinp=0 - first call (on entry)
c iinp=1 - new wavelength is selected (on exit)
c iinp=2 - nothing was changed on first call (on exit). Will be called again.
c iinp=3 - not first call (on entry); will remain on exit and no changes.
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Real            wave
        Integer         iinp, ifail,
     *                  lwrk, i, j, mc,
     *                  io_status,
     *                  lrad, lwav
        Character       byt

        Character       radiat*6
        Common  /x0pa1/ radiat

        Integer         nwav, kwav, ismwav(nwavMax), idwav
        Character       lstwav(nwavMax)*6
        Common  /wavx0/ lstwav, nwav, kwav, ismwav, idwav

        Integer         luwa, luat, lusp, luco
        Common  /lunds/ luwa, luat, lusp, luco

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

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'x0h/RedWav'         !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+

        i = kcompMax            ! Prevents GNU fortran -Wextra warnings
        i = kolMax              ! Prevents GNU fortran -Wextra warnings
        i = ncodMax             ! Prevents GNU fortran -Wextra warnings
        i = natmMax             ! Prevents GNU fortran -Wextra warnings
        i = nwavMax             ! Prevents GNU fortran -Wextra warnings
        i = nedgMax             ! Prevents GNU fortran -Wextra warnings
        i = nabcMax             ! Prevents GNU fortran -Wextra warnings
        i = nedgTb              ! Prevents GNU fortran -Wextra warnings
        i = maxX0hWarnLines     ! Prevents GNU fortran -Wextra warnings
c -------------------------------------------------------
        x0hpath(lx0hpat+1:80) = 'wave.x0h'
        lwrk = lx0hpat+8
        iinp = 0
c If wavelength is passed, then nothing to do:              
        if (wave.gt.0.) Then !---+
          iinp = 1              !|
          goto 998        !return|    
        endif !------------------+
c--------------------------------------------------------   
c If line name is passed, then try to read:                 
        if (Len_Trim(radiat).gt.0) Then  !-----------------+
c Compare passed line name with the                        |
c list of codes in 'wave.x0h':                             |
          lrad = Len_Trim(radiat)                         !|
          do    i=1,kwav  !====================+           |
            lwav = max(Len_Trim(lstwav(i)),1) !|           |
            if (radiat(1:lrad) .eq.           !|           |
     *          lstwav(i)(1:lwav)) goto 15 !---+-FOUND--+  |
          enddo  !=============================+        v  |
c Code is not found. Print error and exit                  |
          goto 191   !-------------->------------ERROR-----+--+
c Code passed from calling program matches the DB. Open    |  v
c the file and scroll to the beginning of line data:    v  |
  15      continue   !<-------------------------<-------+  |
          Call OpenFile(x0hpath(1:lwrk),luwa,'read',      !|
     *                  'old',io_status,*186)             !|
          j = ismwav(i)-1                                 !|
          if (j.gt.0)     Then  !-------------------+      |
            do    i=1,j  !========================+ |      |
              read (luwa,7,err=190,end=193) byt  !| |      |
            enddo   !=============================+ |      |
  7         format(a1,a,g16.8)                     !|      |
          endif  !----------------------------------+      |
          read (luwa,'(a)',err=190,end=193) txt(20)       !|
          Close (unit=luwa,err=190)                       !|
          if (txt(20)(1:1).ne.'#') goto 190  !-------------+--+
          mc = Index(txt(20),';')                         !|  v
          if (mc.gt.0) txt(20)(mc:80) = ' '               !|
          Call TabExpan(txt(20),1)                        !|
          read (txt(20),7,err=190,end=193) byt,radiat,wave!|
          iinp = 1                                        !|
          goto 998                                !return  |
        else  !--------------------------------------------+
          lrad = max(Len_Trim(radiat),1)                  !|
          write   (txt,1) wave, radiat(1:lrad)            !|
  1       format  (
     *    ' Select X-ray wavelength by one of:'//
     *    '- Name of characteristic line,'/
     *    '- Wave length in Angstroms,'/
     *    '- X-ray energy in Kev.'//
     *    'No wavelength information provided!'/
     *    'Currently: wavelength=',g12.5,'  line=',a)     !|
          ifail = 8                                       !|
          goto 998                                        !|
        endif  !-------------------------------------------+
                                                             !v
  998   continue   !<--------------------------------RETURN---+
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c---------------------------------------------------------
c ERROR processing:

  186   continue
        write   (txt,500)       x0hpath(1:lwrk), io_status
  500   format  (' "',a,'"'//' - open error, io_status=',i3)
        ifail = 3                                !number of lines
        goto 998
c..........
  190   continue
        write   (txt,600)       x0hpath(1:lwrk)
  600   format  (' "',a,'"'//' - read error')
        ifail = 3                                !number of lines
        Close   (unit=luwa,err=998)
        goto 998
c..........
  193   continue
        Close   (unit=luwa,err=191)
  191   continue
        write   (txt,700)       x0hpath(1:lwrk),radiat
  700   format  (' "',a,'"'//' - x-ray line [',a,'] - not found')                                      !^
        ifail = 3                                !number of lines
        goto 998

        end
c
c =======================================================
c
        Subroutine      RedCor  (code,iimode,rho,ifail)
c -------------------------------------------------------
c Read crystal lattice parameters for crystal which belong
c to any syngony (triclinic, monoclinic, trigonal, tetragonal,
c hexagonal, or cubic). The lattice parameters and atomic
c coordinates are read from 'coord.x0h' under the crystal
c name. The coordinates must be present in the orthogonal
c settings.
c ........................................
c The structure of crystal record in 'coord.x0h':
c
c code  -     crystal name (the database key).
c kcomp,isym,[rho] - number of components, syngony index, and density.
c                    Syngony indices:
c                     0,1 - cubic
c                     2   - hexagonal/trigonal
c                     3   - tetragonal
c                     4   - trigonal (rhombohedral)
c                     5   - orthorhombic
c                     6   - monoclinic
c                     7   - triclinic
c a1...,a6,Poisson - lattice parameters and angles depending on isym:
c                    0,1 - a(1) [, Poisson]              - Cubic
c                    2   - a(1),a(3)                     - Hexagonal/Trigonal
c                    3   - a(1),a(3)                     - Tetragonal
c                    4   - a(1),a(4)                     - Trigonal(Rhombohedral)
c                    5   - a(1),a(2),a(3)                - Orthorhombic
c                    6   - a(1),a(2),a(3),a(5)           - Monoclinic
c                    7   - a(1),a(2),a(3),a(4),a(5),a(6) - Triclinic
c ----------- Then in the loop over kcomp follow the parameter
c ----------- for each crystal crystal component:
c name(i), [mdeb(i),ddeb(i,1),ddeb(i,2)]
c          - name of the component per "atom.x0h"
c          - [Debye mode and the Debye parameters] (optinal).
c            If present, these Debye parameters override
c            respective parametes in "atom.x0h".
c kol(i)   - number of atomic positions of this component in the unit
c             cell.
c wx(i,1:kol(i)),wy(i,1:kol(i)),wz(i,1:kol(i))
c          - atomic coordinates.
c
c Comment lines can be present in between the description lines.
c It must be started with the symbol ';'.
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Real            rho, rkol(kcompMax),
     *                  rpar(7), prc
        Integer         iimode, ifail, jfail,
     *                  lwrk, lcod, i, j, k, l,
     *                  kodin, mc, nFileData, io_status,
     *                  n_latt(7) /2,2,2,2,3,4,6/
        Character       code*(*)
c Number of lattice constants for different syngonies
c (+Poisson for cubic):
        Integer         n_Column
        External        n_Column

        Character       name(kcompMax)*4
        Common  /x0pa2/ name

        Real            prcn(kcompMax)
        Common  /x0pa3/ prcn

        Real            ddeb(kcompMax,2)
        Integer         mdeb(kcompMax)
        Common  /x0pa5/ ddeb, mdeb

        Integer         kcomp, kol(kcompMax)
        Common  /x0pa6/ kcomp, kol

        Real            wx(kcompMax,kolMax),
     *                  wy(kcompMax,kolMax),
     *                  wz(kcompMax,kolMax)
        Integer         isym
        Common  /x0pa7/ wx, wy, wz, isym

        Real            Poisson, edge_nearest
        Character       name_nearest*4
        common  /x0pa9/ Poisson, edge_nearest, name_nearest

        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        Integer         natm, katm, ismatm(natmMax), idatm
        Character       lstatm(natmMax)*4
        Common  /atmx0/ lstatm, natm, katm, ismatm, idatm

        Integer         ncod, kcod, kelem, ismcod(ncodMax), idcod
        Character       lstcod(ncodMax)*20
        Common  /codx0/ lstcod, ncod, kcod, kelem, ismcod, idcod

        Integer         luwa, luat, lusp, luco
        Common  /lunds/ luwa, luat, lusp, luco

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

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'x0h/RedCor'         !|
        else !-----------------------------+
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
c       write (*,*) 'Redcor: code=',code,' iimode=',iimode              !DEBUG
        iimode = 0
        x0hpath(lx0hpat+1:80) = 'coord.x0h'
        lwrk = lx0hpat+9

        lcod = Len_Trim(code)
        if (lcod.eq.0) Then  !-------------------+no code specified!
          txt(1) = 'No crystal code specified!' !|
          ifail  = 1                            !|
          goto 998  !---------->-----------------+--+Return-+
        endif !----------------------------------+          v
c ........................................
c Crystal code is passed:
c (1) Compare the code with the list of keys in 'coord.x0h':
        do      i=1,kcod  !================+
          if (code.eq.lstcod(i)) Then  !-+ |
            kodin  = i                  !| |
            goto 149   !-----------------+-+----------->------+code
          endif  !-----------------------+ |                  |found!
        enddo  !===========================+                  |
c Code not found in 'coord.x0h'.                              v
c (2) Compare the code with the list of keys in 'atom.x0h':   |
        if (lcod.le.4) Then  !------------------------+       |
          do      i=1,katm  !=====================+   |       |
            if (code(1:4).eq.lstatm(i)) Then  !-+ |   |       |
c cc          kodtyp = 2                       !| |   |       |
              isym   = 8           !amorphous! !| |   |       |
              kcomp  = 1           !amorphous! !| |   |       |
              name(1)= code        !amorphous! !| |   |       |
              kol(1) = 1           !amorphous! !| |   |       |
              goto 999   !----------------------+-+---+-->----+Return-+
            endif  !----------------------------+ |   |       |       v
          enddo  !================================+   |       |
        endif  !--------------------------------------+       |
c Code not found again:                                       |
c (3) Trying to parse the code as a chemical formula!         |
        Call ChiParse2 (Code,kcompMax,kcomp,name,rkol,ifail) !|
c Error parsing:                                              |
        if (ifail.ne.0) Then  !-------------------------+     |
          write (txt,75)                               !|     |
     *          x0hpath(1:lx0hpat)//'Coord.X0h',       !|    149
     *          x0hpath(1:lx0hpat)//'Atom.X0h',        !|     |
     *          code(1:lcod)                           !v     v
  75      format (2('"',a,'"'/)/
     *    ' - structure code [',a,'] not found in these databases.'/
     *    'Also it cannot not be parsed as a chemical formula.'/
     *    'Be aware that X0h DB is case-sensitive.')   !v     v
          ifail = 6             !number of lines        |     |
          goto 999  !------------>----------ERROR-------+-+   |
        endif  !----------------------------------------+ v   |
c Verifying if all the atomic names in the                    |
c structure are valid:                                        |
        Call CheckChiNames (kcomp,name,ifail)                !|
        if (ifail.ne.0) Then  !-------------------------+     |
          j = Len_Trim(name(ifail))                    !|     |
          write (txt,76)                               !|     |
     *           code(1:lcod),                         !|    149
     *           name(ifail)(1:j),                     !|     |
     *           x0hpath(1:lx0hpat)//'Atom.X0h'        !v     v
  76      format (
     *    'Structure code [',a,'] parsed as a chemical formula,'/
     *    'but element [',a,'] not found in the database file: '/
     *    /
     *    '"',a,'"')                                   !v     v
          ifail = 4             !number of lines!       |     |
          goto 999  !------------>----------ERROR-------+-+   |
        endif  !----------------------------------------+ v   |
        isym   = 8           !amorphous!                   ! 149
c There are no records in database for this structure.        |
c Therefore, at least density is required:                    |
        if (rho.le.0.)  Then !--------------------------+     |
          write (txt,71) code(1:lcod)                  !v    !v
  71      format (
     *    'Structure code [',a,'] is not found in the database.'/
     *    'It is parsed as a chemical formula, but the material'/
     *    'density is required to complete the calculations.')!
          ifail =-3             !number of lines!       v     |
          goto 999  !------------>----------ERROR-------+-+   |
        endif  !----------------------------------------+ v   |
c Finally all the controls passed!                           !|
                                                           ! 149
c We use 'prcn' instead of kol because we need real numbers.  |
c+-----------------------------------------+                  |
c|This prevents using the calling-program  |                  |
c|specified 'prcn' with chemical formulas! |                  |
c+=========================================+                  |
        do k=1,kcomp !===============================+        |
          kol(k) = 1                                !|        |
          if (abs(prcn(k)-0.).gt.1.E-10 .and.       !|        |
     *        abs(prcn(k)-1.).gt.1.E-10) goto 555 !--+---error+----+
          prcn(k) = rkol(k)                         !|        |    v
        enddo  !=====================================+        |
c Structure is valid:                                         |
        goto 999   !-----------------------------return-->----+----+
c ........................                                    |    v
  149   continue    !<---------<---------------------<----<---+
c ........................................
c Code is entered: open the file and scroll to the
c beginning of crystal description:
        Call OpenFile(x0hpath(1:lwrk),luco,'read','old',
     *                io_status,*30)    !-------------------ERROR---------+
c ........................................
        j = ismcod(kodin)
        do      i=1,j  !================================+
          read  (luco,3,err=40,end=50)  txt(20)(1:1) !--+---ERROR---------+
        enddo  !========================================+
  3     format  (a)
c --------------
c Read crystal parameters:
  4     continue  !<--------------------------------+
        read    (luco,3,err=40,end=50)  txt(20)  !--+-------ERROR---------+
        if (txt(20)(1:1).eq.';')    goto 4 !--------+
c Strip end-of-line comments:                       ^
        mc = Index(txt(20),';')                    !|
        if (mc.gt.0) txt(20)(mc:80) = ' '          !^
        Call TabReplace (txt(20),1)                !|
c cccc  if (Len_Trim(txt(20)).eq.0) goto 4 !--------+
c kcomp,isym,[rho] -- the first two parameters are integers,
c but we read them as reals:
        Call RdReal (rpar,3,txt(20),jfail)
        if (jfail.ne.0)                goto 40   !----------ERROR---------+
        kcomp = int (rpar(1)+0.1)
        isym  = int (rpar(2)+0.1)
        if (rho.le.0.) rho=rpar(3)
        if (rho.le.0. .AND. isym.eq.8) goto 150  !----------ERROR-missing rho--+
        if (isym.lt.0 .or. isym.gt.8)  goto 141  !----------ERROR---------+
        if (isym.eq.0)  isym=1
        if (kcomp.gt.kcompMax)         goto 61   !----------ERROR---------+
c ++++++++++++++++++++++++++++++++++++++++++++++
c Read lattice parameters:
        if (isym.ne.8)  Then !---------------------------------------+
  5       continue  !<--------------------------------+              |
          read (luco,3,err=40,end=50) txt(20)  !------+-------ERROR--+----+
          if (txt(20)(1:1).eq.';')    goto 5  !-------+              |
c Strip end-of-line comments:                         ^              |
          mc = Index(txt(20),';')                    !|              |
          if (mc.gt.0) txt(20)(mc:80) = ' '          !^              |
          Call TabReplace (txt(20),1)                !|              |
c cccc    if (Len_Trim(txt(20)).eq.0) goto 5  !-------+              |
          nFileData = n_Column(txt(20))                             !|
          if (isym.eq.1)  Then !----------------------+              |
            if (nFileData.gt.n_latt(isym))  goto 143 !|may be Poisson|
            if (nFileData.lt.1)             goto 143 !|              |
          else  !-------------------------------------+              |
            if (nFileData.ne.n_latt(isym))  goto 143 !|no Poisson    |
          endif  !------------------------------------+              |
          Call    RdReal  (rpar,7,txt(20),jfail)                    !|
          if (jfail.ne.0) goto 40   !-------------------------ERROR--+----+
          if (abs(a(1)).lt.1.E-10) a(1)=rpar(1)                     !|
                                                                    !|
c 0,1 - cubic                                                        |
c 2   - hexagonal/trigonal                                           |
c 3   - tetragonal                                                   |
c 4   - trigonal (rhombohedral)                                      |
c 5   - orthorhombic                                                 |
c 6   - monoclinic                                                   |
c 7   - triclinic                                                    |
c                1  2  3  4  5  6  7                                 |
          goto (80,81,81,82,83,84,85) isym  !=====================+  |
c 1. Cubic structure (may contain the Poisson ratio)          VVVVV  |
  80      continue                               !================+  |
          if (abs(Poisson).lt.1.E-10) Poisson=rpar(2)               !|
          goto 88     !-------------->----------------+              |
c 2. Hexagonal/Trigonal structure:                    |              |
c 3. Tetragonal structure                             |       VVVVV  |
  81      continue                               !====+===========+  |
          if (abs(a(3)).lt.1.E-10) a(3)=rpar(2)      !|              |
          goto 88     !-------------->----------------+              |
c Trigonal (rhombohedral) structure:                  |       VVVVV  |
  82      continue                               !====+===========+  |
          if (abs(a(4)).lt.1.E-10) a(4)=rpar(2)      !|              |
          goto 88     !-------------->----------------+              |
c 5. Orthorhombic structure:                          |       VVVVV  |
  83      continue                               !====+===========+  |
          if (abs(a(2)).lt.1.E-10) a(2)=rpar(2)      !|              |
          if (abs(a(3)).lt.1.E-10) a(3)=rpar(3)      !v              |
          goto 88     !-------------->----------------+              |
c 6. Monoclinic structure:                            |       VVVVV  |
  84      continue                               !====+===========+  |
          if (abs(a(2)).lt.1.E-10) a(2)=rpar(2)      !|              |
          if (abs(a(3)).lt.1.E-10) a(3)=rpar(3)      !v              |
          if (abs(a(5)).lt.1.E-10) a(5)=rpar(4)      !|              |
          goto 88     !-------------->----------------+              |
c 7. Triclinic structure:                             |       VVVVV  |
  85      continue                               !====+===========+  |
          do      i=2,6 !==========================+  |              |
            if (abs(a(i)).lt.1.E-10) a(i)=rpar(i) !|  v              |
          enddo  !=================================+  |              |
                                                     !|              |
  88      continue   !<-------------------------------+              |
        endif  !-----------------------------------------------------+

        if (kcomp.gt.0) Then  !-------------------------------+
c Read atomic coordinates:   :                                |
          do      i=1,kcomp  !==============================+ |
  6         continue  !<-----------------------------+      | |
            read  (luco,3,err=40,end=50)  txt(20)   !|      |-|ERROR------+
            if (txt(20)(1:1).eq.';')    goto 6  !----+      | |
c Strip end-of-line comments:                        ^      | |
            mc = Index(txt(20),';')                 !|      | |
            if (mc.gt.0) txt(20)(mc:80) = ' '       !|      | |
            Call TabExpan (txt(20),1)               !|      | |
c cccc      if (Len_Trim(txt(20)).eq.0) goto 6  !----+      | |
c                                                           | |
c Name_of_Element [, Debye_mode, par1, par2 ]               | |
c       name      [, mdeb, ddeb(1), ddeb(2) ]               | |
c For Debye_mode=1                                          | |
c           par1=Debye coefficient                          | |
c           par2=Crystal temperature for which it is given  | |
c                (by default room temperature is used).     | |
c For Debye_mode=2,3                                        | |
c           par1=Debye temperature                          | |
c           par2 - not used                                 | |
c                                                           | |
            if (name(i).eq.' ') name(i)=txt(20)(1:4)       !| |
            l = Len_Trim(txt(20)(5:80))                    !| |
            if (l.gt.0)   Then  !----------------------+    | |
c name, [ mdeb, ddeb(1), ddeb(2) ] -- mdeb is integer, |    | |
c but we read it as real:                              |    | |
              Call RdReal (rpar,3,txt(20)(5:80),jfail)!|    | |
              if (jfail.ne.0)               goto 77   !|    |-|ERROR------+
              mdeb(i) = int (rpar(1)+0.1)             !|    | |
              if (mdeb(i)*(3-mdeb(i)).lt.0) goto 77   !|    |-|ERROR------+
              if (mdeb(i).gt.0) Then  !-+              |    |-|ERROR------+
                ddeb(i,1) = rpar(2)    !|              |    | |
                ddeb(i,2) = rpar(3)    !|              |    | |
              endif  !------------------+              |    | |
            endif  !-----------------------------------+    | |
c                                                           | |
  7         continue  !<-----------------------------+      | |
            read (luco,3,err=40,end=50)  txt(20)    !|      |-|ERROR------+
            if (txt(20)(1:1).eq.';')    goto 7  !----+      | |
c Strip end-of-line comments:                        ^      | |
            mc = Index(txt(20),';')                 !|      | |
            if (mc.gt.0) txt(20)(mc:80) = ' '       !^      | |
            Call TabReplace (txt(20),1)             !|      | |
c cccc      if (Len_Trim(txt(20)).eq.0) goto 7  !----+      | |
c                                                           | |
c Nr_of_positions_of_element_1 [, Occupation]               | |
c    kol, [prc] -- kol is integer, but we read it as real:  | |
c                                                           | |
            Call RdReal (rpar,2,txt(20),jfail)             !| |
            if (jfail.ne.0)       goto 40                  !|-|ERROR------+
            kol(i) = int (rpar(1)+0.1)                     !| |
            prc    = rpar(2)                               !|-|ERROR------+
            if (kol(i).lt.1)      goto 50  !----------------+-|ERROR------+
            if (kol(i).gt.kolMax .AND.                     !| |
     *            isym.ne.8)      goto 51  !----------------+-|ERROR------+
c The range is: [0.-1.]:                                   !| |
            if (prc.lt.0. .OR.                             !| |
     *          prc.gt.1.)        goto 91  !----------------+-|ERROR------+
c Take the database value of occupation unless it is        | |
c specified by the calling program:                         | |
            if (abs(prcn(i)).lt.1.E-10) prcn(i)=prc        !| |
            if (isym.ne.8) Then !-------------------------+ | |for crystals only
              do k=1,kol(i)  !==========================+ | | |
  8             continue  !<--------------------------+ | | | |
                read  (luco,3,err=40,end=50) txt(20) !| | | |-|ERROR------+
                if (txt(20)(1:1).eq.';')    goto 8 !--+ | | | |
c Strip end-of-line comments:                         ^ | | | |
                j = Index(txt(20),';')               !| | | | |
                if (j.gt.0) txt(20)(j:80) = ' '      !| | | | |
                Call TabReplace (txt(20),1)          !| | | | |
c We may have x,y,z as empty string:                  ^ | | | |
c cccc          if (Len_Trim(txt(20)).eq.0) goto 8 !--+ | | | |
c x,y,z:                                                | | | |
                Call RdReal  (rpar,3,txt(20),jfail)    !| | | |
                if (jfail.ne.0)     goto 40  !----------+-+-+-|ERROR------+
                wx(i,k) = rpar(1)                      !| | | |
                wy(i,k) = rpar(2)                      !| | | |
                wz(i,k) = rpar(3)                      !| | | |
              enddo  !==================================+ | | |
            endif  !--------------------------------------+ | |
          enddo  !==========================================+ |
                                                             !|
        else  !-----------------------------------------------+
                                                             !|
c If no elements in the structure description:                |
                                                             !|
c If this is a 'no-elements description' structure (template):|
          if (kodin.gt.kelem) goto 998 !------------return-----------+
          if (isym.ne.8)      goto 63  !------------error------------+---+
                                                             !|     !v   v
c Try to parse the amorphous structure name with zero number  |      |
c of components in order to get elements names & their number:|      |
          Call ChiParse2 (Code,kcompMax,kcomp,               !|      |
     *                    name,rkol,ifail)                   !|      |
          if (ifail.ne.0) goto  145  !--------------error-----+------+---+
c Do these names exist in Atom.x0h?                           |      |   v
          Call CheckChiNames (kcomp,name,ifail)              !|      |
          if (ifail.ne.0) goto  147  !--------------error-----+------+---+
c We use 'prcn' instead of kol                                |      |
c because we need real numbers.                               |      |
c+----------------------------------------+                   |      |
c|This prevents using the calling-program |                   |      |
c|specified 'prcn' with chemical formulas!|                   |      |
c+========================================+                   |      |
          do k=1,kcomp !================================+     |      |
            kol(k) = 1                                 !|     |      |
            if (abs(prcn(k)-0.).gt.1.E-10 .and.        !|     |      |
     *          abs(prcn(k)-1.).gt.1.E-10) goto 555 !---+-----+------+---+error
            prcn(k) = rkol(k)                          !|     |      |   v
          enddo  !======================================+     |      |
        endif  !----------------------------------------------+      |
                                                                    !V
  998   continue !<---------------------------------------RETURN-----+
        Close   (unit=luco,err=999)
  999   continue
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
c       write (*,*) 'Redcor: rho=',rho                                  !DEBUG
        return
c ----------------------------------------------------
c    ERRORS processing:
c ........................
  30    continue
        write   (txt,35)        x0hpath(1:lwrk), io_status
  35    format  ('"',a,'"'//' - open error, io_status=',i3)
        ifail = 3                                !number of lines
        goto 999
c ........................
  40    continue
        write   (txt,45)        x0hpath(1:lwrk), code(1:lcod)
  45    format  ('"',a,'"'//' - read error',
     *  ' at structure [',a,']')
        ifail = 3                                !number of lines
        goto 998
c ........................
  50    continue
        write   (txt,55)        x0hpath(1:lwrk), code(1:lcod)
  55    format  ('"',a'"'//a,' - uncomplete data',
     *  ' for the structure')
        ifail = 3                                !number of lines
        goto 998
c ........................
  63    continue
        write   (txt,64)        x0hpath(1:lwrk), code(1:lcod)
  64    format  ('"',a,'"'//a,' - zero number of components')
        ifail = 3                                !number of lines
        goto 998
c ........................
  61    continue
        write   (txt,62)        x0hpath(1:lwrk), code(1:lcod),
     *                          kcompMax
  62    format  ('"',a,'"'//a,' - number of components > ',i3)
        ifail = 3                                !number of lines
        goto 998
c ........................
  51    continue
        write   (txt,52)        x0hpath(1:lwrk), code(1:lcod),
     *                          i, kolMax
  52    format  ('"',a,'"'//a,' - number of coords for',
     *  ' component ',i2,' exceeds ',i3)
        ifail = 3                                !number of lines
        goto 998
c ........................
  91    continue
        write   (txt,92)        x0hpath(1:lwrk), code(1:lcod),
     *                          i,name(i)
  92    format  ('"',a,'"'//a,' - occupation for component ',i2,
     *  ' ( ',a,')'/'not in range [0%-100%]')
        ifail = 4                                !number of lines
        goto 998
c ........................
  141   continue
        write   (txt,142)       x0hpath(1:lwrk), code(1:lcod)
  142   format  ('"',a,'"'//a,' - invalid symmetry code')
        ifail = 3                                !number of lines
        goto 998
c ........................
  143   continue
        write   (txt,144)       x0hpath(1:lwrk), code(1:lcod),
     *                          nFileData
  144   format  ('"',a,'"'//a,' - wrong number of lattice',
     *  ' constants (',i1,')')
        ifail = 3                                !number of lines
        goto 998
c ........................
  77    continue
        write   (txt,87)        x0hpath(1:lwrk), code(1:lcod),
     *                          name(i)
  87    format  ('"',a,'"'//a,' - invalid Debye-Waller',
     *  ' factor for "',a,'"')
        ifail = 3                                !number of lines
        goto 998
c ........................
  145   continue
        txt(5) = txt(1)
        write   (txt,146)       x0hpath(1:lwrk), code(1:lcod)
  146   format  ('"',a,'"'//'Syntax error in chemical formula: [',
     *                                                     a,']'/)
        ifail = 5                                !number of lines
        goto 998
c ........................
  147   continue
        i = Len_Trim(name(ifail))
        write   (txt,247)       x0hpath(1:lwrk), code(1:lcod),
     *                          name(ifail)(1:i)
  247   format  ('"',a,'"'//'Error parsing ',a/
     *           'Atom [',a,'] not found in X0H database')
        ifail = 5                                !number of lines
        goto 998
c ........................
  150   continue
        write   (txt,151)       x0hpath(1:lwrk), code(1:lcod)
  151   format  ('"',a,'"'//'Missing, zero, or negative density data '/
     *           'in the database for the amorphous structure:'/a)
        ifail = 5                                !number of lines
        goto 998
c ........................
  555   continue
        write   (txt,152)       x0hpath(1:lwrk), code(1:lcod), name(k),
     *                          prcn(k)
  152   format  ('"',a,'"'//
     *  'X0H cannot override occupations of chemical formulas!'/
     *  'Chemical formula=',a/
     *  'Element=',a,',  received occupation=',f7.3)
        ifail = 5                                !number of lines
        goto 998
c ............................
        end
c
c =======================================================
c
        Subroutine      RedLat  (code,need_a,ifail)
c -------------------------------------------------------
c Request lattice parameters from terminal with providing
c help data (lattice parameters for various structures)
c from 'space.x0h'.
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Integer         need_a(6), ifail,
     *                  lwrk, jm, i, j
        Character       code*(*)

        Real            wx(kcompMax,kolMax),
     *                  wy(kcompMax,kolMax),
     *                  wz(kcompMax,kolMax)
        Integer         isym
        Common  /x0pa7/ wx, wy, wz, isym

        Real            a(6), vol
        Integer         lbev
        Common  /back1/ a, vol, lbev

        Integer         luwa, luat, lusp, luco
        Common  /lunds/ luwa, luat, lusp, luco

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

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'x0h/RedLat'         !|
        else !-----------------------------+
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
        x0hpath(lx0hpat+1:80) = 'space.x0h'
        lwrk = lx0hpat+9

        need_a(1) = 1
        do      i=2,6  !=====+
          need_a(i)=0       !|
        enddo   !============+
        goto    (88,81,81,82,83,84,85)  isym
c Hexagonal/trigonal or tetragonal structure:
  81    continue
        need_a(3) = 1
        goto 88
c Trigonal (rhombohedral) structure:
  82    continue
        need_a(4) = 1
        goto 88
c Orthorhombic structure:
  83    continue
        need_a(2) = 1
        need_a(3) = 1
        goto 88
c Monoclinic structure:
  84    continue
        need_a(2) = 1
        need_a(3) = 1
        need_a(5) = 1
        goto 88
c Triclinic structure:
  85    continue
        need_a(2) = 1
        need_a(3) = 1
        need_a(4) = 1
        need_a(5) = 1
        need_a(6) = 1
c Loop over needed lattice parameters:
  88    continue
        jm = 0
        do      i=1,6  !================================+
          if (need_a(i).ne.0 .and. a(i).le.0.) Then !-+ |
            jm    = jm+1                             !| |
            if (i.le.3) write   (txt(jm+2),2)  i     !| |
            if (i.gt.3) write   (txt(jm+2),22) i-3   !| |
  2         format('Spacing [',i1,']:  ',10('_'))    !| |
  22        format('  Angle [',i1,']:  ',10('_'))    !| |
          endif  !------------------------------------+ |
        enddo  !========================================+

        if (jm.ne.0)    Then  !-------------------------------+
          j = Len_Trim(code)                                 !|
          if (j.lt.1) j=1                                    !|
          txt(1) = code(1:j)//': Enter lattice constant(s):' !|
          txt(2) = ' '                                       !|
          txt(jm+3) = ' '                                    !|
          txt(jm+4)='Incomplete structures cannot be'//      !|
     *              ' used with this version of X0h!'        !|
          ifail = -(jm+4)                                    !|
          goto 998                                           !|
        endif  !----------------------------------------------+
c --------------
c 1 - cubic
c 2 - hexagonal/trigonal
c 3 - tetragonal
c 4 - trigonal (rhombohedral)
c 5 - orthorhombic
c 6 - monoclinic
c 7 - triclinic
c                 1  2  3  4  5  6  7
        goto    (91,92,93,94,95,96,998)  isym
c 1. Cubic structure:
  91    continue
        a(2) = a(1)
        a(3) = a(1)
        a(4) = 90.
        a(5) = 90.
        a(6) = 90.
        goto 998
c 2. Hexagonal/trigonal structure:
  92    continue
        a(2) = a(1)
        a(4) = 90.
        a(5) = 90.
        a(6) = 120.
        goto 998
c 3. Tetragonal structure:
  93    continue
        a(2) = a(1)
        a(4) = 90.
        a(5) = 90.
        a(6) = 90.
        goto 998
c 4. Trigonal (rhombohedral) structure:
  94    continue
        a(2) = a(1)
        a(3) = a(1)
        a(5) = a(4)
        a(6) = a(4)
        goto 998
c 5. Orthorhombic structure:
  95    continue
        a(4) = 90.
        a(5) = 90.
        a(6) = 90.
        goto 998
c 6. Monoclinic structure:
  96    continue
        a(4) = 90.
        a(6) = 90.
c 7. Triclinic structure:
  998   continue
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
c
c =======================================================
c
        Subroutine      RedPar  (code,rho,atom_mass_aem,ifail)
c -------------------------------------------------------
c Read atoms scattering parameters from 'atom.x0h'
c ........................................
c The structure of records in 'atom.x0h':
c
c namfil - name of atom (database key).
c koza,ifac,idis,isig,ideb - number of record lines for the atom and
c           indices of methods for calculating structure amplitudes,
c           dispersion corrections, absorption correction and the
c           Debye-Waller factors.
c --------- After that the parameters follow depending on chosen methods
c --------- (see the comments below in the process of reading)
c
c Comment lines between data lines are allowed. They must start
c with the ';' symbol.
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c       Parameter       (nedgTb   = 5)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Real            rho, atom_mass_aem, atom_mass_gram,
     *                  energy, gap, Planck_cm, tt,
     *                  lfix(nedgTb), ldat(nedgTb), dmur(nedgTb),
     *                  dfv1(nedgTb), dfv2(nedgTb), del(nedgTb),
     *                  av(2,2), bv(2), det, zz, zerr, c, d, s,
     *                  T_DWF, T_sample, T_Debye,
     *                  rpar(11)
        Integer         ifail, jfail, i_old(nedgTb),
     *                  iab, idb, jHenkeCowan,
     *                  ifac, idis, isig, ideb,
     *                  ipar(4), lwrk, nval,
     *                  mc, i, j, k, m, io_status
        Character       code*(*), byt*1

        Real            FiFunc, wave2energy
        External        FiFunc, wave2energy
        Integer         N_Column
        External        N_Column

        Character       name(kcompMax)*4
        Common  /x0pa2/ name

        Real            wave, df1(kcompMax), df2(kcompMax),
     *                  ab(kcompMax,nabcMax), s1, s21, s31, s32, s41,
     *                  e_radius, Compton, Rydberg, pi,
     *                  Planck_, Boltzman, Avogadro,
     *                  wedg(nedgMax), gedg(nedgMax)
        Integer         indx, ipr, nedg
        Common  /x0pa4/ wave, df1, df2, ab, s1, s21, s31, s32, s41,
     *                  e_radius, Compton, Rydberg, pi, indx, ipr,
     *                  Planck_, Boltzman, Avogadro, wedg, gedg, nedg

        Real            ddeb(kcompMax,2)
        Integer         mdeb(kcompMax)
        Common  /x0pa5/ ddeb, mdeb

        Real            bdw(kcompMax), f0(kcompMax),
     *                  td(kcompMax), tq(kcompMax)
        Integer         z(kcompMax)
        Common  /x0pa8/ bdw, f0, td, tq, z

        Real            Poisson, edge_nearest
        Character       name_nearest*4
        common  /x0pa9/ Poisson, edge_nearest, name_nearest

        Integer         natm, katm, ismatm(natmMax), idatm
        Character       lstatm(natmMax)*4
        Common  /atmx0/ lstatm, natm, katm, ismatm, idatm

        Integer         luwa, luat, lusp, luco
        Common  /lunds/ luwa, luat, lusp, luco

        Integer         lx0hpat
        Character       x0hpath*80
        Common  /x0hpa/ x0hpath, lx0hpat

        Character       txt(20)*80
        Common  /msg/   txt
c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer             iHenkeCowan
        Common /HenkeCowan/ iHenkeCowan
c Allowed iHenkeCowan values:
c -1 - automatic choice of internal or external dispersion correction DB
c  0 - do not use external dispersion correction databases
c  1 - use henke.dat    for df1
c  2 - use henke.dat    for df1, df2 (Crosec)
c  3 - use cowan.dat    for df1
c  4 - use cowan.dat    for df1, df2 (Crosec)
c  5 - use windt.dat    for df1
c  6 - use windt.dat    for df1, df2 (Crosec)
c  7 - use chantler.dat for df1
c  8 - use chantler.dat for df1, df2 (Crosec)
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------
c X-ray wavelengths in Angstrom for the lines
c Cr-Ka1, Fe-Ka1, Cu-Ka1, Mo-Ka1, Ag-Ka1:
        Data    ldat    /2.28962,               ! E= 5.415
     *                   1.93597,               ! E= 6.404
     *                   1.54052,               ! E= 8.048
     *                   0.70926,               ! E=17.480
     *                   0.55936/               ! E=22.165
c -- So:
c X0h (International Tables) range:       0.56  --    2.29 Angstrem
c          (5,000ev -- 25,000ev):         0.50  --    2.47 Angstrem
c Henke range (10ev -- 30,000ev):         0.41  -- 1239.81 Angstrem
c Cowan range (30ev -- 694,500ev):        0.02  --  413.27

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'x0h/RedPar'         !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+

        iab = 0                                 !make GNU compiler happy
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
c------------------------------------------------------------
        T_sample = 0.
c------------------------------------------------------------
        x0hpath(lx0hpat+1:80) = 'atom.x0h'
        lwrk = lx0hpat+8
c Copy the global external database
c flag to the one for a given element:
        jHenkeCowan = iHenkeCowan
        if (jHenkeCowan .eq. -1) Then !-------------+automatic choice
          if     ( wave .gt. 10.) Then !---------+  |
             jHenkeCowan = 2                    !|  |use Henke data
          elseif ( wave .lt. 0.1) Then !---------+  |
             jHenkeCowan = 4                    !|  |use Cowan data
          else  !--------------------------------+  |
             jHenkeCowan = 0                    !|  |use X0h data
          endif !--------------------------------+  |
c Auto-choice does not include Windt or Chantler~   |
        endif !-------------------------------------+
c If the atom code is not passed:
        if (Len_Trim(name(indx)).eq.0) goto 1 !--ERROR---->--+
c If the atom code is passed, the compare                    |
c it with the DB keys in 'atom.x0h':                         |
          do      i=1,katm  !=========================+      v
            if (name(indx).eq.lstatm(i))  goto 149 !--+->-+  |
          enddo  !====================================+   |  |
c Code is not found, report error:                        v  |
          goto 1  !------------------------------ERROR----+--+
                                                         !|  v
  149   continue  !<------------------<-------------------+
        j = ismatm(i)
c ........................................
c Code is present: open 'atom.x0h' and scroll to the beginning
c of atom description (to the scattering parameters):
        Call OpenFile(x0hpath(1:lwrk),luat,'read','old',
     *                io_status,*500)
c ........................................
        do      i=1,j  !======================+
          read  (luat,43,err=11,end=12) byt  !|
        enddo  !==============================+
  43    format  (a)
c --------------
c Read the parameters:
  2     continue  !<--------------------------------+
        read (luat,43,err=11,end=12) txt(20)       !|
        if (txt(20)(1:1).eq.';')    goto 2  !-------+
c Strip end-of-line comments:                       ^
        mc = Index(txt(20),';')                    !|
        if (mc.gt.0) txt(20)(mc:80) = ' '          !^
        Call TabReplace (txt(20),1)                !|
c cccc  if (Len_Trim(txt(20)).eq.0) goto 2  !-------+
        Call RdInt (ipar,4,txt(20),jfail)
        if (jfail.ne.0)             goto 11 !--error--+
        ifac = ipar(1)                               !V
        idis = ipar(2)
        isig = ipar(3)
        ideb = ipar(4)
        if (ifac.lt.1 .or. ifac.gt.2) goto 130 !-error+
        if (idis.lt.1 .or. idis.gt.3) goto 130 !-error|
        if (isig.lt.1 .or. isig.gt.5) goto 130 !-error|
        if (ideb.lt.1 .or. ideb.gt.3) goto 130 !-error|
c ++++++++++++++++++++++++++++++++++++++++++++++      V
c Read the parameters for calculating f0,fh:
c z,             - element order number in the Periodic Table
c                  (the number of electron it has)
c atom_mass_aem, - atom weight in atomic mass units.
c rho            - mass density in grams/cm^3

  131   continue  !<--------------------------------+
        read (luat,43,err=11,end=12) txt(20)       !|
        if (txt(20)(1:1).eq.';')    goto 131  !-----+
c Strip end-of-line comments:                       ^
        mc = Index(txt(20),';')                    !|
        if (mc.gt.0) txt(20)(mc:80) = ' '          !^
        Call TabReplace (txt(20),1)                !|
c cccc  if (Len_Trim(txt(20)).eq.0) goto 131  !-----+
        Call Rdreal (rpar,3,txt(20),jfail)
        if (jfail.ne.0)           goto 11  !---error--+
                                                     !V
        z(indx)       = int(rpar(1)+0.1)   !rounding
        atom_mass_aem = rpar(2)
        rho           = rpar(3)
        if (atom_mass_aem.le.0.) goto 940  !---error--+
        if (rho.le.0.)           goto 942  !---error--+
                                                     !V

        if     (ifac.eq.1) Then  !-----------------------------+
c Method-1. Read coefficients for the 4-term interpolation     |
c formula.                                                     |
          iab = 9                                             !|
          ab(indx,10) = 0.                                    !|
          ab(indx,11) = 0.                                    !|
                                                              !|
        elseif (ifac.eq.2) Then  !-----------------------------+
c Method-2. Read coefficients for the 5-term interpolation     |
c formula.                                                     |
c /New of August 99/                                           |
          iab = 11                                            !|
                                                              !|
        endif  !-----------------------------------------------+

        i = 0
  3     continue  !<-------------------------------+
        read (luat,43,err=11,end=12) txt(20)      !|
        if (txt(20)(1:1).eq.';')    goto 3  !------+
c Strip end-of-line comments:                      ^
        mc = Index(txt(20),';')                   !|
        if (mc.gt.0) txt(20)(mc:80) = ' '         !^
        Call TabReplace (txt(20),1)               !|
        j = N_Column(txt(20))                     !|
        if (j.lt.1)              goto 3  !---------+
        if (j.gt.11)             goto 11 !---error-+-+
c Strict control:                                  | v
        if ((i+j).gt.iab)        goto 11 !---error-+-+
        Call RdReal (rpar,j,txt(20),jfail)        !| v
        if (jfail.ne.0)          goto 11 !---error-+-+
        do      m=1,j !===========+                | v
          ab(indx,i+m) = rpar(m) !|                |
        enddo  !==================+                |
        i = i+j                                   !|
        if (i.lt.iab)            goto 3  !---------+
c Calculate f0 using equation for fh at s=0 (see Forfac):
        f0(indx) = ab(indx,iab)
        do      i=1,iab-2,2  !================+
          f0(indx) = f0(indx) + ab(indx,i)   !|
        enddo  !==============================+
        zz   = z(indx)
        zerr = 0.03
        if (zz.gt.39.)  zerr=0.05
        if (zz.gt.72.)  zerr=0.06
        if (zz.gt.81.)  zerr=0.07
        if (abs(zz-f0(indx)).gt.zerr) goto 53 !-------error-----+
                                                               !V

c ++++++++++++++++++++++++++++++++++++++++++++++
c Read data for calculating dispersion corrections:
c nedg - number of absorption edges for the edges (=<nedgMax=16),
c wedg - energies of absorption edges in KeV.

  132   continue  !<-------------<------------------+
        read (luat,43,err=11,end=12) txt(20)       !|
        if (txt(20)(1:1).eq.';')    goto 132  !-----+
c Strip end-of-line comments:
        mc = Index(txt(20),';')
        if (mc.gt.0) txt(20)(mc:80) = ' '
        Call TabReplace (txt(20),1)
        Call RdInt (ipar,1,txt(20),jfail)
        nedg = ipar(1)
        if (jfail.ne.0)                      goto 11  !----error----+
        if (nedg.lt.0 .or.  nedg.gt.nedgMax) goto 935 !----error----+
        if (nedg.eq.0 .and. z(indx).ge.6)    goto 925 !----error----+
                                                                   !v
        if (nedg.gt.0)      Then  !----------------------+
c Read absorption edges (kev):                           |
          i = 0                                         !|
  4       continue  !<-------------------------------+   |
          read (luat,43,err=11,end=12) txt(20)      !|   |
          if (txt(20)(1:1).eq.';')    goto 4  !------+   v
c Strip end-of-line comments:                        ^   |
          mc = Index(txt(20),';')                   !|   |
          if (mc.gt.0) txt(20)(mc:80) = ' '         !^   |
          Call TabReplace (txt(20),1)               !|   |
          if (Len_Trim(txt(20)).eq.0) goto 4  !------+   |
          j = N_Column(txt(20))                     !|   |
          if (j.lt.1)              goto 4 !----------+   |
          if (j.gt.11)             goto 11 !---error-+---+-+
c Strict control:                                    |   | v
          if ((i+j).gt.nedg)       goto 11 !---error-+---+-+
          Call RdReal (rpar,j,txt(20),jfail)        !|   | v
          if (jfail.ne.0)          goto 11 !---error-+---+-+
          do m=1,j  !==========+                     |   | v
            wedg(i+m)=rpar(m) !|                     |   |
          enddo  !=============+                     |   |
          i = i+j                                   !|   |
          if (i.lt.nedg)           goto 4  !---------+   |
c Find the closest absorption edge                       |
c (this is used for reference and to estimate            |
c the reliability of interpolation)                      |
          if (wave.gt.0.) Then  !-----------------+      |
            energy  = wave2energy(wave)          !|      |
            gap     = Abs(energy-edge_nearest)   !|      |
            do i=1,nedg  !=====================+  |      |
              s = Abs(energy-wedg(i))         !|  |      |
              if (s.lt.gap) Then !----------+  |  |      V
                gap = s                    !|  |  |      |
                edge_nearest = wedg(i)     !|  |  |      |
                name_nearest = name(indx)  !|  |  |      |
              endif  !----------------------+  |  |      |
            enddo  !===========================+  |      |
          endif  !--------------------------------+      |
        endif  !-----------------------------------------+

c 1-read the oscillators strengths,
c 2-read df1,df2 for 5 x-ray lines and then calculate the oscillators
c   strengths from these data,
c 3-read 5 wavelengths, then df1,df2 for them and then calculate the
c   oscillators strengths from these data.
c Set default zero values (2008/05/27):
        do  k=1,nedgTb !============+
c Take tabulated X-ray wavelengths: |
          lfix(k) = ldat(k)        !|
          dfv1(k) = 0.             !|
          dfv2(k) = 0.             !|
        enddo  !====================+
        if (idis.eq.1)  Then  !-----------------------------+
c Method-1. Read the oscillators strengths at the absorption|
c edges. In this method df' and df" are not read.           |
          if (nedg.gt.0) Then  !-------------------------+  |
            i = 0                                       !|  |
  5         continue  !<----------------------------+    |  |
            read (luat,43,err=11,end=12) txt(20)   !|    |  |
            if (txt(20)(1:1).eq.';')    goto 5  !->-+    |  |
c Strip end-of-line comments:                       ^    |  |
            mc = Index(txt(20),';')                !|    |  |
            if (mc.gt.0) txt(20)(mc:80) = ' '      !^    |  |
            Call TabReplace (txt(20),1)            !|    |  |
            if (Len_Trim(txt(20)).eq.0) goto 5  !->-+    |  |
            j = N_Column(txt(20))                  !|    |  |
            if (j.lt.1)          goto 5  !----------+    |  |
            if (j.gt.11)         goto 11 !----error-+-+  |  |
c Strict control:                                   | v  |  |
            if ((i+j).gt.nedg)   goto 11 !----error-+-+  |  |
            Call RdReal (rpar,j,txt(20),jfail)     !| v  |  |
            if (jfail.ne.0)      goto 11 !----error-+-+  |  |
            do m=1,j  !==========+                  | v  |  |
              gedg(i+m)=rpar(m) !|                  |    |  |
            enddo  !=============+                  |    |  |
            i = i+j                                !|    |  |
            if (i.lt.nedg)       goto 5  !---->-----+    |  |
          endif  !---------------------------------------+  |
                                                           !|
        elseif (idis.eq.2) Then  !--------------------------+
c Method-2. Read df',df" for 5 characteristic X-ray lines   |
c Cr-Ka1, Fe-Ka1, Cu-Ka1, Mo-Ka1, Ag-Ka1 and then calculate |
c the oscillators strengths with the help of Disint.        |
c Take the tabulated X-ray wavelengths:                     |
c         do  k=1,nedgTb !=======+                          |
c           lfix(k) = ldat(k)   !|                          |
c         enddo  !===============+                          |
c Values of df' and df" are read further in this routine    |
                                                           !|
        elseif (idis.eq.3) Then  !-------------<------------+
c Method-3. Read df',df" for 5 wavelengths specified in     |
c 'atom.x0h' and then calculate the oscillators strengths   |
c with the help of Disint.                                  |
c Read the wavelengths:                                     |
          i = 0                                            !|
  17      continue  !<---------------<------------+         |
          read (luat,43,err=11,end=12) txt(20)   !|         |
          if (txt(20)(1:1).eq.';')    goto 17 !---+         |
c Strip end-of-line comments:                     ^         |
          mc = Index(txt(20),';')                !|         |
          if (mc.gt.0) txt(20)(mc:80) = ' '      !^         |
          Call TabReplace (txt(20),1)            !|         |
          if (Len_Trim(txt(20)).eq.0) goto 17 !---+         |
          j = N_Column(txt(20))                  !|         |
          if (j.lt.1)            goto 17 !--------+         |
          if (j.gt.11)           goto 11 !--error-+-+       |
c Strict control:                                 | v       |
          if ((i+j).gt.nedgTb)   goto 11 !--error-+-+       |
          Call RdReal (rpar,j,txt(20),jfail)     !| v       |
          if (jfail.ne.0)        goto 11 !--error-+-+       |
          do m=1,j  !==========+                  | v       |
            lfix(i+m)=rpar(m) !|                  |         |
          enddo  !=============+                  |         |
          i = i+j                                !|         |
          if (i.lt.nedgTb)      goto 17 !---------+         |
c Values of df' and df" are read further in this routine    |
                                                           !|
        endif  !--------------------------------------------+

        if     (idis.eq.2 .OR. idis.eq.3) Then  !-----------+
c Read df', df" for 5 wavelengths (methods 2 and 3):        |
          i = 0                                            !|
  6       continue  !<---------------------------+          |
          read (luat,43,err=11,end=12) txt(20)  !|          |
          if (txt(20)(1:1).eq.';')    goto 6  !--+          |
c Strip end-of-line comments:                    ^          |
          mc = Index(txt(20),';')               !|          |
          if (mc.gt.0) txt(20)(mc:80)  =' '     !^          |
          Call TabReplace (txt(20),1)           !|          |
          if (Len_Trim(txt(20)).eq.0) goto 6  !--+          |
          j = N_Column(txt(20))                 !|          |
          if (j.lt.1)          goto 6  !---------+          |
          if (j.gt.11)         goto 11 !---error-+-+        |
c Strict control:                                | v        |
          if ((i+j).gt.nedgTb) goto 11 !---error-+-+        |
          Call RdReal (rpar,j,txt(20),jfail)    !| v        |
          if (jfail.ne.0)      goto 11 !---error-+-+        |
          do  m=1,j !=============+              | v        |
            dfv1(i+m) = rpar(m)  !|              |          |
          enddo  !================+              |          |
          i = i+j                               !|          |
          if (i.lt.nedgTb)     goto 6  !--->-----+          |
          i = 0                                            !|
  16      continue  !<---------------------------+          |
          read (luat,43,err=11,end=12) txt(20)  !|          |
          if (txt(20)(1:1).eq.';')    goto 16 !--+          |
c Strip end-of-line comments:                    ^          |
          mc = Index(txt(20),';')               !|          |
          if (mc.gt.0) txt(20)(mc:80) = ' '     !^          |
          Call TabReplace (txt(20),1)           !|          |
          if (Len_Trim(txt(20)).eq.0) goto 16 !--+          |
          j = N_Column(txt(20))                 !|          |
          if (j.lt.1)          goto 16 !---------+          |
          if (j.gt.11)         goto 11 !---error-+-+        |
c Strict control:                                | v        |
          if ((i+j).gt.nedgTb) goto 11 !---error-+-+        |
          Call RdReal (rpar,j,txt(20),jfail)    !| v        |
          if (jfail.ne.0)      goto 11 !---error-+-+        |
          do  m=1,j !=============+              | v        |
            dfv2(i+m) = rpar(m)  !|              ^          |
          enddo  !================+              |          |
          i = i+j                               !|          |
          if (i.lt.nedgTb)     goto 16 !--->-----+          |
        endif  !--------------------------------------------+

c If external database for df',df" to be used:
c Acceptable values of iHenkeCowan:
c -1 - automatic choice of internal or external dispersion correction DB
c  0 - do not use external dispersion correction databases
c  1 - use henke.dat    for df1
c  2 - use henke.dat    for df1, df2 (Crosec)
c  3 - use cowan.dat    for df1
c  4 - use cowan.dat    for df1, df2 (Crosec)
c  5 - use windt.dat    for df1
c  6 - use windt.dat    for df1, df2 (Crosec)
c  7 - use chantler.dat for df1
c  8 - use chantler.dat for df1, df2 (Crosec)

        idb = 0
        if     (jHenkeCowan.ne.0) Then !----------------------------+
          if     (jHenkeCowan.eq.1 .OR. jHenkeCowan.eq.2) Then !-+  |
            idb = 1        !'henke.dat'                          |  |
          elseif (jHenkeCowan.eq.3 .OR. jHenkeCowan.eq.4) Then !-+  |
            idb = 2        !'cowan.dat' or 'cowanlng.dat'        |  |
          elseif (jHenkeCowan.eq.5 .OR. jHenkeCowan.eq.6) Then !-+  |
            idb = 3        !'windt.dat'                          |  |
          elseif (jHenkeCowan.eq.7 .OR. jHenkeCowan.eq.8) Then !-+  |
            idb = 4        !'chantler.dat'                       |  |
          else  !------------------------------------------------+  |
            goto 943 !-----------------------------------error---+--+--+
          endif  !-----------------------------------------------+  |  v
                                                                   !|
          nval = nedgTb                                            !|
          Call f1f2 (idb, wave, name(indx), Z(indx), nval,         !|
     *               lfix, dfv1, dfv2, ipr, jfail)                 !|
          if     (jfail.gt.0) Then !-------------------+            |
c Error:                                               |            |
            ifail = jfail                             !|            |
            goto 998  !--------------------------------+---return---+--+
          elseif (jfail.lt.0) Then !-------------------+            |  v
c No Henke/Cowan/Windt/Chantler data for this element: |            |
            jHenkeCowan = 0                           !|            |
            jfail = 0                                 !|            |
            idb = 0                                   !|            |
          endif  !-------------------------------------+            |
        endif  !----------------------------------------------------+

c Acceptable values of idis (dispersion mode, 1-3)
c 1: Cromer formula for f" with known oscillator strengths.
c    No reading Cowan or Henke tables is required in this case.
c 2: Cromer formula for f" with interpolation of oscillator
c    strengths using tabulated df',df" for 5 fixed x-ray lines
c    (Cr-Ka1, Fe-Ka1, Cu-Ka1, Mo-Ka1, Ag-Ka1).
c 3: Cromer formula for f" with interpolation of oscillator
c    strengths using tabulated df',df" for 5 x-ray wavelengths
c    listed for this element in atom.x0h

        if (idis.ne.1 .or. jHenkeCowan.ne.0) Then  !----+
c Find (interpolate) the oscillators strengths 'gedg'   |
c by applying the oscillators formula to known          |
c dispersion corrections at 5 wavelengths:              |
c (Methods 2 and 3 or if external DB for                |
c df',df" is used):                                     |
          Call Disint (dfv1,dfv2,lfix,jfail)           !|
          if     (jfail.gt.0) Then !--+                 |
            ifail = jfail            !|                 |
            goto 998  !---------------+------return-----+----+
          elseif (jfail.lt.0) Then !--+                 |    v
            jfail = 0                !|                 |
            do i=1,nedgMax !==+       |                 |
              gedg(i) = 0.   !|       |                 |
            enddo  !==========+       |                 |
          endif !---------------------+                 |
        endif  !----------------------------------------+

c Calculate dispersion correction according to
c D.T.Cromer - Acta Cryst.(1965),18,p.17-23.
c/Change of August 99: the Cramer method does not
c work when atoms does not have absorption edges in
c the X-ray range (applies to light atoms). In this
c case another simple interpolation method is applied --
c see respective changes in subroutine Disper/:
        nval = nedgTb
        Call  Disper (atom_mass_aem, idis, jHenkeCowan,
     *                dfv1, dfv2, lfix, nval, ifail)
        if (ifail.ne.0) goto 998  !--------------return------+
c                                                            v

c ++++++++++++++++++++++++++++++++++++++++++++++
c Read parameters for calculating dipole and quadrupole
c terms of the absorption cross section.

        if     (isig.eq.1) Then  !--------------------------+
c Method-1. Read the screening constants.                   |
c+--------------------------+                               |
c|This is for medium energy |                               |
c|range (5-25kev) and Z=6-54|                               |
c|as declared in the paper! |                               |
c+--------------------------+                               |
  651     continue  !<-----------------<-------------+      |
          read (luat,43,err=11,end=12)  txt(20)     !|      |
          if (txt(20)(1:1).eq.';')    goto 651  !----+      |
c Strip end-of-line comments:                        ^      |
          mc = Index(txt(20),';')                   !|      |
          if (mc.gt.0) txt(20)(mc:80) = ' '         !^      |
          Call TabReplace (txt(20),1)               !|      |
c cccc    if (Len_Trim(txt(20)).eq.0) goto 651  !----+      |
          Call RdReal (rpar,5,txt(20),jfail)               !|
          if (jfail.ne.0)          goto 11 !--error---------+--+
          s1  = rpar(1)                                    !|  v
          s21 = rpar(2)                                    !|
          s31 = rpar(3)                                    !|
          s32 = rpar(4)                                    !|
          s41 = rpar(5)                                    !|
                                                           !|
        elseif (isig.eq.2) Then  !--------------------------+
c Method-2. Read the c,d coefficients.                      |
  652     continue  !<-----------------<-------------+      |
          read (luat,43,err=11,end=12)  txt(20)     !|      |
          if (txt(20)(1:1).eq.';')    goto 652  !----+      |
c Strip end-of-line comments:                        ^      |
          mc = Index(txt(20),';')                   !|      |
          if (mc.gt.0) txt(20)(mc:80) = ' '         !|      |
          Call TabReplace (txt(20),1)               !|      |
c cccc    if (Len_Trim(txt(20)).eq.0) goto 652  !----+      |
          Call RdReal (rpar,2,txt(20),jfail)               !|
          if (jfail.ne.0)          goto 11 !--error---------+--+
          c = rpar(1)                                      !|  v
          d = rpar(2)                                      !|
                                                           !|
        elseif (isig.eq.3) Then  !--------------------------+
c Method-3. Read mu/rho for 5 characteristic X-ray lines    |
c Cr-Ka1, Fe-Ka1, Cu-Ka1, Mo-Ka1, Ag-Ka1                    |
c and then calculate c,d using interpolation.               |
          do  k=1,nedgTb !=======+                          |
            lfix(k) = ldat(k)   !|                          |
          enddo  !===============+                          |
                                                           !|
        elseif (isig.eq.4) Then  !--------------------------+
c Methos-4. Read mu/rho for 5 wavelengths from 'atom.x0h'   |
c and then calculate c,d using interpolation.               |
                                                           !|
c Read the wavelengths:                                     |
          i = 0                                            !|
  7       continue  !<---------------<---------------+      |
          read (luat,43,err=11,end=12) txt(20)      !|      |
          if (txt(20)(1:1).eq.';')    goto 7  !------+      |
c Strip end-of-line comments:                        ^      |
          mc = Index(txt(20),';')                   !|      |
          if (mc.gt.0) txt(20)(mc:80) = ' '         !^      |
          Call TabReplace (txt(20),1)               !|      |
          if (Len_Trim(txt(20)).eq.0) goto 7  !------+      |
          j = N_Column(txt(20))                     !|      |
          if (j.lt.1)              goto 7  !---------+      |
          if (j.gt.11)             goto 11 !---error-+------+--+
c Strict control:                                    |      |  v
          if ((i+j).gt.nedgTb)     goto 11 !-error---+------+--+
          Call RdReal (rpar,j,txt(20),jfail)        !|      |  v
          if (jfail.ne.0)          goto 11 !-error---+------+--+
          do m=1,j  !==========+                     |      |  v
            lfix(i+m)=rpar(m) !|                     |      |
          enddo  !=============+                     |      |
          i = i+j                                   !|      |
          if (i.lt.nedgTb)         goto 7  !---------+      |
                                                           !|
c       elseif (isig.eq.5) Then  !--------------------------+
c Method-5. Do not read anything and determine td from the  |
c imaginary part of dispersion correction - implemented in  |
c subroutine Crosec.                                        |
                                                           !|
        endif  !--------------------------------------------+

        if     (isig.eq.3 .or. isig.eq.4) Then  !-----------+
c Read mu/rho (cm^2/g) for Methods 3,4:                     |
          i = 0                                            !|
  8       continue  !<---------------<---------------+      |
          read (luat,43,err=11,end=12) txt(20)      !^      |
          if (txt(20)(1:1).eq.';')    goto 8  !------+      |
c Strip end-of-line comments:                        ^      |
          mc = Index(txt(20),';')                   !|      |
          if (mc.gt.0) txt(20)(mc:80) = ' '         !|      |
          Call TabReplace (txt(20),1)               !|      |
          if (Len_Trim(txt(20)).eq.0) goto 8  !------+      |
          j = N_Column(txt(20))                     !|      |
          if (j.lt.1)              goto 8  !---------+      |
          if (j.gt.11)             goto 11 !---------+-error+--+
c Strict control:                                    |      |  v
          if ((i+j).gt.nedgTb)     goto 11 !---------+-error+--+
          Call RdReal (rpar,j,txt(20),jfail)        !|      |  v
          if (jfail.ne.0)          goto 11 !---------+-error+--+
          do m=1,j  !==========+                     |      |  v
            dmur(i+m)=rpar(m) !|                     |      |
          enddo  !=============+                     |      |
          i = i+j                                   !|      |
          if (i.lt.nedgTb)         goto 8  !---------+      |
c Interpolate c,d for Mathods 3,4:                          |
          do  i=1,nedgTb !===================+              |
            del(i)   = abs(wave-lfix(i))    !|              |
            i_old(i) = i                    !|              |
          enddo  !===========================+              |
c Sort the wavelength differences in the increasing order:  |
          Call  Sort_differences (del, i_old, nedgTb)      !|
c del(1) is the smalest difference (the closest wavelength),|
c del(2) is the next one, .. and so on.                     |
c i_old(i) contain the numbers del(i) had before sorting.   |
          do      m=1,2  !==============+                   |
            j = i_old(m)               !|                   |
            bv(m)   = dmur(j)          !|                   |
            av(m,1) = lfix(j)**3       !|                   |
            av(m,2) = lfix(j)**4       !|                   |
          enddo  !======================+                   |
c We use:                                                  !|
c mu1 = c*L1^3 - d*L1^4                                    !|
c mu2 = c*L2^3 - d*L2^4                                    !|
c where L1, L2 are the two closest waves.                   |
c This gives c,d:                                           |
          det = av(2,1)*av(1,2)-av(2,2)*av(1,1)            !|
          if (Abs(det).lt.1.E-20)       goto 141 !--error---+--+
          c = (bv(2)*av(1,2)-bv(1)*av(2,2))/det            !|  v
          d = (bv(2)*av(1,1)-bv(1)*av(2,1))/det            !|
                                                           !|
        endif  !--------------------------------------------+

c /New of August-99: ignore isig in
c the X0h database and set isig=5 if
c jHenkeCowan=2 or jHenkeCowan=4/:
        if (jHenkeCowan.eq.2 .OR.
     *      jHenkeCowan.eq.4 .OR.
     *      jHenkeCowan.eq.6 .OR.
     *      jHenkeCowan.eq.8) Then !-----+
          isig = 5                      !|
        endif  !-------------------------+

c Calculate dipole and quadrupole parts of the absorption cross section.
c The screening parameters s1,s21,s31,s32,s41 are passed to Crosec via
c common/x0pa4/.
c The calculation method depends on isig !
        Call    CroSec  (atom_mass_aem,
     *                   c,d,z(indx),
     *                   td(indx),tq(indx),
     *                   isig, idb, atom_mass_aem, jfail)
        if (jfail.gt.0) Then !-------------------+
c Error:                                         |
          ifail = jfail                         !|
          goto 998  !----------------------------+---return---+
        endif !----------------------------------+            v

c ++++++++++++++++++++++++++++++++++++++++++++++
c Read data for calculating the Debye-Waller factor.

c The Debye mode can be specified in structure description in
c 'coord.x0h' (the same atom may be specified with different Debye
c modes in different crystals). If such mode is specified, it overrides
c the mode for the atom in 'atom.x0h'.
        i = mdeb(indx)
        if ((i-1)*(3-i).ge.0)   ideb=i
        if     (ideb.eq.1) Then  !--------------------------+
c Method-1. Read the Debye coefficient and the temperature  |
c for which it is specified. By default room temperature is |
c assumed.                                                  |
  654     continue  !<--------------<----------------+      |
          read (luat,43,err=11,end=12) txt(20)      !|      |
          if (txt(20)(1:1).eq.';')    goto 654  !----+      |
c Strip end-of-line comments:                              !|
          mc = Index(txt(20),';')                          !|
          if (mc.gt.0) txt(20)(mc:80) = ' '                 !|
          Call TabReplace (txt(20),1)                      !|
          Call RdReal (rpar,2,txt(20),jfail)               !|
          if (jfail.ne.0)             goto 11 !-------------+--+ error
          bdw(indx) = rpar(1)                              !|  v
          T_DWF     = rpar(2)                              !|
          if (mdeb(indx).ne.0)    Then  !----+              |
            bdw(indx) = ddeb(indx,1)        !|              |
            T_DWF     = ddeb(indx,2)        !|              |
          endif  !---------------------------+              |
          if (abs(T_DWF)   .lt.1.E-10) T_DWF=293.          !|
          if (abs(T_sample).lt.1.E-10) T_sample=T_DWF      !|
          if (abs(T_sample-T_DWF).gt.1.E-10) goto 938 !-----+--+ error
                                                           !|  v
        elseif (ideb.eq.2) Then  !--------------------------+
c Method-2. Crystal is assumed to be at room temperature and|
c the Debye temperature is read further is this routine.    |
          T_DWF = 293.                         !just for now|
          if (abs(T_sample).lt.1.E-10) T_sample=293.       !|
          if (abs(T_sample-293.).gt.1.E-10) goto 938 !------+--+ error
                                                           !|  v
        elseif (ideb.eq.3) Then  !--------------------+
c Method-3. Request sample temperature (atom's Debye  |
c temperature is read further is this routine).       |
          T_DWF = 293.                  !just for now |
          if (abs(T_sample).lt.1.E-10) Then !------+  |
                                                  !|  |
            txt(1)='Enter sample temperature'//   !|  |
     *                       ' (Kelvin degr.):'   !|  |
            txt(2) = ' '                          !|  |
            txt(3) = 'This request is not '//     !|  |
     *               'supported in the '//        !|  |
     *               'current version of X0h.'    !|  |
              ifail = 3                           !|  |
              goto 998  !------------ERROR---------+--+--+
          endif  !---------------------------------+  |  v
        endif !---------------------------------------+

c+---------------------------------------------------+
c|Methods 2-3. Calculate the Debye coefficient using |
c|the Debye temperature and crystal temperature.     |
c+---------------------------------------------------+
        if     (ideb.eq.2 .or. ideb.eq.3) Then  !-----------+
c Methods 2-3. Read Debue temperature (in Kelvin degrees).  |
  655     continue  !<--------------<----------------+      |
          read (luat,43,err=11,end=12) txt(20)      !|      |
          if (txt(20)(1:1).eq.';')    goto 655  !----+      |
c Strip end-of-line comments:                               |
          mc = Index(txt(20),';')                          !|
          if (mc.gt.0) txt(20)(mc:80) = ' '                !|
          Call TabReplace (txt(20),1)                      !|
          Call RdReal (rpar,1,txt(20),jfail)               !|
          if (jfail.ne.0)             goto 11 !------error--+--+
          T_Debye = rpar(1)                                !|  v
          if (mdeb(indx).ne.0)  T_Debye=ddeb(indx,1)       !|
          if (T_Debye.lt.0.001) T_Debye=0.001              !|
          atom_mass_gram = atom_mass_aem/Avogadro          !|
          if (abs(atom_mass_gram).lt.1.E-37) goto 928 !-----+--+ error
c Calculate the Debye coefficient:                          |  v
          tt = T_Debye/T_sample                            !|
c Planck_cm - the Planck constant in cm without line bar.   |
          Planck_cm = 2.*pi*(1.e+08)*Planck_               !|
          bdw(indx) = 12.*(Planck_cm/atom_mass_gram)       !|
     *              * (Planck_cm/Boltzman)/T_Debye         !|
          bdw(indx) = bdw(indx)*(FiFunc(tt)/tt+0.25)       !|
                                                           !|
        endif !---------------------------------------------+

c Output the result of calculating bdw:
        if (ipr.eq.2)   Then  !------------------------+
          write (txt,577)       name(indx),           !|
     *                          bdw(indx),            !|
     *                          ideb                  !v
  577     format('Element: ',a,' -- RedPar data'/
     *    4x,'B =',g14.6,'   Debye-Waller mode =',i2) !^
          Call  Priblo  (txt,2,1,0,0,i)               !|
        endif  !---------------------------------------+


c ++++++++++++++++++++++++++++++++++++++++++++++
c Account for dispersion correction for f0 (amplitude
c of scattering at angle=0) and print the result:

        f0(indx) = f0(indx) + df1(indx)
        if (ipr.eq.2)   Then  !---------------------------+
          write (txt,57)        name(indx),              !|
     *                          f0(indx),                !|
     *                          ifac                     !v
  57      format('Element: ',a,' -- RedPar data'/
     *    4x,'f0=fh(0)=',f8.4,'   Formfactor mode =',i2) !^
          Call  Priblo  (txt,2,1,0,0,i)                  !|    |
        endif  !------------------------------------------+    |
                                                              !v
  998   continue  !<-------------------------------------------+
        Close   (unit=luat,err=999)                           !v|
  999   continue  !<-------------------------------------------+
        x0hpath(lx0hpat+1:80) = ' '
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c ==========================================================
c                     ERRORS processing:
c ..........................................................
  53    continue
        write   (txt,52)        x0hpath(1:lwrk),name(indx),
     *                          f0(indx),z(indx),
     *                          (ab(indx,i),i=1,9)
  52    format  (' "',a,'"'//' - inconsistent A,B,C for',
     *  ' element [',a,']:'/
     *  15x,'f0=',g15.7,'  #  Z=',i5/
     *  'a1,b1,a2,b2: ',4g12.4/
     *  'a3,b3,a4,b4: ',4g12.4/
     *  'c: ', g12.4)
        ifail = 7
        goto 998
c ........................
  500   continue
        write   (txt,800)       x0hpath(1:lwrk), io_status
  800   format  (' "',a,'"'//' - open error, io_status=',i3)
        ifail = 3
        goto 999
c ........................
  11    continue
        write   (txt,900)       x0hpath(1:lwrk),name(indx)
  900   format  (' "',a,'"'//' - read error',
     *  ' at element [',a,']')
        goto 998
c ........................
  130   continue
        write   (txt,902)       x0hpath(1:lwrk),name(indx),
     *                          ifac,idis,isig,ideb
  902   format(' "',a,'"'//
     *  ' Invalid mode(s) in X0h database for atom [',a,']'/
     *  ' Scattering amplitude calculation mode=',i2,' [1--2]'/
     *  ' Dispersion corrections calculat. mode=',i2,' [1--3]'/
     *  ' Absorbtion cross section calcat. mode=',i2,' [1--5]'/
     *  ' Debye-Waller factor calculation  mode=',i2,' [1--3]')
        ifail = 7
        goto 998
c ........................
  12    continue
        write   (txt,13)        x0hpath(1:lwrk),name(indx)
  13    format  (' "',a,'"'//' - uncomplete data',
     *  ' for element [',a,']')
        ifail = 3
        goto 998
c ........................
  925   continue
        write   (txt,927)       x0hpath(1:lwrk),name(indx)
  927   format  (' "',a,'"'//' - zero absorption edges',
     *  ' for [',a,']')
        ifail = 3
        goto 998
c ........................
  935   continue
        write   (txt,936)       x0hpath(1:lwrk),name(indx),
     *                          nedg, nedgMax
  936   format  (' "',a,'"'//' - icorrect number of absorption',
     *  'edges for [',a,']:'/
     *  ' n_edges=',i2,'  n_edges_max=',i2)
        ifail = 4
        goto 998
c ........................
  928   continue
        write   (txt,929)       x0hpath(1:lwrk),name(indx)
  929   format  (' "',a,'"'//' - zero atomic weight',
     *  ' for [',a,']')
        ifail = 3
        goto 998
c ........................
  938   continue
        write   (txt,939)       x0hpath(1:lwrk),name(indx),
     *                          T_sample, T_DWF, ideb
  939   format  (' "',a,'"'//' - Debye-Waller data',
     *  ' temperature mismatch for [',a,']:'/
     *  ' T_sample = ',g9.3,'  T_DWF=',g9.3,'  DWF_mode=',i1)
        ifail = 4
        goto 998
c ........................
  141   continue
        write   (txt,241)       x0hpath(1:lwrk),name(indx)
  241   format  (' "',a,'"'//
     *  ' - Incorrect dispersion data in the file resulted in'/
     *  'zero determinant while calculating c,d factors for [',a,']')
        ifail = 4
        goto 998
c ........................
  940   continue
        write   (txt,929)       x0hpath(1:lwrk),name(indx)
c 929   format  (' "',a,'"'//' - zero atomic weight',
c    *  ' for [',a,']')
        ifail = -3
        goto 998
c ........................
  942   continue
        write   (txt,923)       x0hpath(1:lwrk),name(indx)
  923   format  (' "',a,'"'//' - zero material density',
     *  ' for [',a,']')
        ifail = -3                              !number of lines
        goto 998
c ........................
  943   continue
        write   (txt,924)       jHenkeCowan
  924   format  ('Flag =',i2,' of external f'',f" database not ',
     *  'in the range 0--8')
        ifail = 1                               !number of lines
        goto 998
c ........................
  1     continue
        i = max(Len_Trim(code),1) 
        write   (txt,33)        x0hpath(1:lwrk),name(indx),code(1:i)
  33    format  (' "',a,'"'//' - element name [',a,
     *  '] - not found for crystal=[',a,']')
        ifail = 3                               !number of lines
        goto 999
c ............................
        end
c
c ==============================================================
c
        Subroutine Insert_Poisson (Poisson,isym,tt)
        Real    Poisson
        Integer isym
        Character tt*(*)
c          1         2         3         4         5         6         7
c 1234567890123456789012345678901234567890123456789012345678901234567890123
c |Material symmetry_____________________________|Cubic___________________|
c |                 , Poisson ratio              |       12345678901234567
        if (isym.le.1) Then !---------------------------+
          tt(19:33)=', Poisson ratio'                  !|
          tt(54:72)=', '                               !|
          if (abs(Poisson).lt.1.E-10) Then  !--------+  |
            tt(56:72)='      **Unknown**'           !|  |
          else  !------------------------------------+  |
            write (tt(56:72),'(f17.4)') Poisson     !|  |
          endif  !-----------------------------------+  |
        endif  !----------------------------------------+
        return
        end
c
c ==============================================================
c
        Subroutine Insert_Edge (edge_nearest, name_nearest,tt)
        Real      edge_nearest
        Character name_nearest*4, tt*(*)
c          1         2         3         4         5         6         7
c 1234567890123456789012345678901234567890123456789012345678901234567890123
c |Energy (Kev)__________________________________|12345678901_____________|
c |            , nearest absorption edge [1234]  |           , 12345678901|
        if (edge_nearest.gt.0.          .AND.
     *      Len_Trim(name_nearest).gt.0) Then  !----------------------+
          tt(14:47)=', nearest absorption edge ['//name_nearest//']' !|
          tt(60:61)=','                                              !|
          write(tt(62:72),'(g11.4)') edge_nearest                    !|
        endif  !------------------------------------------------------+
        return
        end
c
c ==============================================================
c
        Subroutine CheckChiNames (kcomp,name,ifail)
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Integer         kcomp, ifail,
     *                  i, j, k, l
        Character       name(kcomp)*(*)

        Integer         natm, katm, ismatm(natmMax), idatm
        Character       lstatm(natmMax)*4
        Common  /atmx0/ lstatm, natm, katm, ismatm, idatm

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

c Checks if specified atomic name(s) are available in Atom.X0h:
        do      k=1,kcomp  !===========================+
          j = Len_Trim(name(k))                       !|
          if (j.gt.4)   goto 1                        !|
          l = Min (4,Len(name(k)))                    !|
          do      i=1,katm  !======================+   |
            if (name(k)(1:l).eq.lstatm(i)) goto 2 !|   |
          enddo  !=================================+   |
          goto 1                                      !|
  2       continue                                    !|
        enddo  !=======================================+
        ifail = 0
        return
  1     ifail = k
        return
        end
c
c ==============================================================
c
        Integer Function CodeType (code)
c Determines and returns code type:
c 0=crystal, 1=other/atom, 2=formula, -1=error
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Character       code*(*)
        Integer         lcod, lwrk, isym, io_status, ifail,
     *                  mc,  i, l
        Real            rpar(3)

        Integer         natm, katm, ismatm(natmMax), idatm
        Character       lstatm(natmMax)*4
        Common  /atmx0/ lstatm, natm, katm, ismatm, idatm

        Integer         ncod, kcod, kelem, ismcod(ncodMax), idcod
        Character       lstcod(ncodMax)*20
        Common  /codx0/ lstcod, ncod, kcod, kelem, ismcod, idcod

        Integer         luwa, luat, lusp, luco
        Common  /lunds/ luwa, luat, lusp, luco

        Integer         lx0hpat
        Character       x0hpath*80
        Common  /x0hpa/ x0hpath, lx0hpat

        Character       txt(20)*80
        Common  /msg/   txt

        lcod = Len_Trim(code)
c If crystal code is not passed:
        if (lcod.eq.0) Then !---+
          CodeType = -1        !|
          return               !|
        endif !-----------------+
c (1) Compare the code with the list of keys in 'coord.x0h':
        do      i=1,kcod  !==========================================+
          l = Len_Trim(lstcod(i))                                   !|
          if (l.gt.0 .AND. code(1:lcod).eq.lstcod(i)(1:l)) Then !-+  |
            ! Open 'coord.x0h' and scroll to the                  |  |
            ! beginning of crystal description:                   |  |
            x0hpath(lx0hpat+1:80) = 'coord.x0h'                  !|  |
            lwrk = lx0hpat+9                                     !|  |
            Call OpenFile(x0hpath(1:lwrk),luco,'read','old',     !|  |
     *                io_status,*30)    !-------------------ERROR-+--+---+
            ! Read crystal parameters:                            |  |   V
  4         continue  !<------------------------------+           |  |
            read (luco,'(a)',err=40,end=40) txt(20) !-+-----ERROR-+--+---+
            if (txt(20)(1:1).eq.';') goto 4 !---------+           |  |   V
            ! Strip end-of-line comments:                         |  |
            mc = Index(txt(20),';')                              !|  |
            if (mc.gt.0) txt(20)(mc:80) = ' '                    !|  |
            Call TabReplace (txt(20),1)                          !|  |
            ! kcomp,isym,[rho] -- the first two parameters        |  |
            ! are integer, but we read them as real:              |  |
            Call RdReal (rpar,3,txt(20),ifail)                   !|  |
            if (ifail.ne.0)          goto 40   !------------ERROR-+--+---+
            isym  = int (rpar(2)+0.1)                            !|  |   V
            if (isym.lt.8) Then !----+                            |  |
              CodeType = 0          !| crystal                    |  |
            else !-------------------+                            |  |
              CodeType = 1          !| other structure            |  |
            endif !------------------+                            |  |
            return                                               !|  |
          endif !-------------------------------------------------+  |
        enddo  !=====================================================+
c (2) Compare the code with the list of keys in 'atom.x0h':
        if (lcod.le.4) Then  !-----------------------------------------+
          do      i=1,katm  !=========================================+|
            l = Len_Trim(lstatm(i))                                  !||
            if (l.gt.0 .AND. code(1:lcod).eq.lstatm(i)(1:l)) Then !-+ ||
              CodeType = 1      ! atom (i.e. other structure)       | ||
              return                                               !| ||
            endif  !------------------------------------------------+ ||
          enddo  !====================================================+|
        endif  !-------------------------------------------------------+
c Code not found again:
c (3) Presuming the code as a chemical formula!
        CodeType = 2
        return
c ----------------------------------------------------
c    ERRORS processing:
c ........................
  40    continue
        Close   (unit=luco,err=30)
  30    continue
        x0hpath(lx0hpat+1:80) = ' '
        CodeType = -1
        return
        end
