        Subroutine  Get_Wave  (Wave,LineName,ifail)
c -------------------------------------------------------
c Find X-ray wavelength using given name of X-ray characteristic line
c
c                   Author: S.Stepanov
c -------------------------------------------------------
c ATTENTION:       
c Data from calling program is also passed via several common blocks:
c       common  /lunds/ luwa,luat,lusp,luco
c       common  /x0hpa/ x0hpath,lx0hpat
c -------------------------------------------------------
c       Input and output parameters:
c
c Wave       -  X-ray wavelength in Angstrom. If zero wavelength is 
c               passed and also if the name of characteristic X-ray
c               line "radiat" is empty, then X0H requests from terminal
c               first the name "radiat" and if not given, then the 
c               wavelength. At exit the wavelength is determined.
c LineName*6 -  Name of X-ray characteristic line 
c .......................................................
c       Optional parameter passed via common /lunds/:
c
c       common  /lunds/ luwa,luat,lusp,luco
c
c luwa,luat,
c lusp,luco  -  Logical Unit Numbers (LUN) of the X0H database files.
c               They must differ from each other (X0H does not control
c               it). If zeroes are passed, X0H uses LUNs 11-14.
c .......................................................
c       Optional parameters passed via common /X0hpa/
c
c       common  /x0hpa/ x0hpath,lx0hpat
c
c x0hpath*80  -  path to directoty with X0H database files.
c lx0hpat     -  number of characters in the path.
c -------------------------------------------------------
c                       Subroutines:
c   redwav  --x0hread--       &         --serve.lib--
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
        Real            Wave
        Integer         ifail
        Character       LineName*6
     *
c -------------------------------------------------------
c Common blocks for passing parameters from calling programs:
        Character       Radiat*6
        Common  /x0pa1/ Radiat

        Integer         luwa,luat,lusp,luco
        Common  /lunds/ luwa,luat,lusp,luco

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c Common block for passing parameters between X0H subroutines:
        Integer         nwav, kwav, ismwav(nwavMax), idwav
        Character       lstwav(nwavMax)*6
        Common  /wavx0/ lstwav,nwav,kwav,ismwav,idwav
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

        Integer         iirezv, itrace, i

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Get_Wave'           !|
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
        if (abs(Wave).gt.1.E-20) goto 100 !-------------+ nothing to do
                                                       !v
c Set X0H database files LUNs if not passed:
        if (luwa.eq.0)  luwa=11         !"wave.x0h"

c Make a copy of X-ray line name:
        Radiat = LineName
        if (Radiat(1:1).eq.Char(0))     Radiat=' '

c This is the wavelengths array dimension:
        nwav = nwavMax

c Make catalog of the database file 'wave.x0h':
        if (idwav.eq.0) Then  !--------------+
          Call Basini ('wave.x0h',lstwav,   !|
     *                 nwav,kwav,ismwav,    !|
     *                 idwav,0,i,ifail)     !|
          if (ifail.ne.0) goto 100          !|
        endif  !-----------------------------+

c Determine the X-ray wavelength:
        itrace = 0
        Call    Redwav (Wave,itrace,ifail)
c We need to add some error processing here
c (currently handled in the x0hread.for function "Redwav"):
        if (ifail.ne.0) goto 100 !----------------------+
                                                       !v
  100   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
