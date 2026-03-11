        Subroutine f1f2 (idb,         !INPUT:  1=Henke, 2=Cowan, 3=Windt, 4=Chantler
     +                   wave,        !INPUT:  wavelength in Angstrem
     +                   Element,     !INPUT:  element code (name)
     +                   Z,           !INPUT:  element's Z
     +                   nval,        !INPUT:  number of requested f',f"
     +                   wave_n,      !OUTPUT: array of waves corr. f',f"
     +                   f1_n, f2_n,  !OUTPUT: arrays of f',f"
     +                   ipr,         !INPUT:  0,1,2 - debug level
     +                   ifail)       !OUTPUT: failure flag
c-----------------------------------------------------------------------
c This subroutine searches for f1,f2 dispersion corrections
c for element 'Element' with Mendeleev number 'Z' in the file
c henke.dat    (idb=1), contains E_eV, (f+df1), df2
c cowanlng.dat (idb=2), contains E_ev,   df1,   df2
c windt.dat    (idb=3), contains E_eV, (f+df1), df2
c chantler.dat (idb=4), contains E_eV, (f+df1), df2
c
c nodata=-.999900 (in old versions was -9999.00)
c
c These files are supposed to be in the DABAX format.
c The program returns 'nval' tabular values of f1_n, f2_n
c for 'nval' closest x-ray wavelengths 'wave_n' expressed
c in Angstrom.
c
c                    Author: S.A.Stepanov
c-----------------------------------------------------------------------
        Integer         nval

        Real*8          E2W     /12398.1/,
     *                  Eps     /1.D-8/,
     *                  EnergyEV, EKeV,
     *                  En, En_saved,
     *                  Wv, Wv_saved,
     *                  Diff, Diff0

        Real            wave,
     *                  wave_n(nval),
     *                  f1_n(nval),
     *                  f2_n(nval),
     *                  rpar(3),
     *                  nodata /-.999900/

        Integer*4       line, lineStart, LineCenter

        Integer         idb, lun, Z, iZ, ipar(1),
     *                  ipr, ifail, ltx, io_status

        Character       File*12, Element*(*), El*2,
     *                  x0hpath*80, txt(20)*80,
     *                  DBtxt(4)*8 /'Henke',
     *                              'Cowan',
     *                              'Windt',
     *                              'Chantler'/

        Integer         Zmin, Zmax,
     *                  lx0hpat,
     *                  modebat, lunbat, ierrbat, istackrezv,
     *                  iirezv, ishift, lwrk, ii, ii_center,
     *                  i_repeat, nmore, i, j

        Real            wave2energy
        External        wave2energy
        Integer         FreeLun
        External        FreeLun

        Common  /x0hpa/ x0hpath,lx0hpat
        Common  /msg/   txt
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'X0h/f1f2'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+
        lwrk = lx0hpat

        if (idb.lt.1 .OR. idb.gt.4) goto 501 !---------------------------+
                                                                         !v
        if (ipr .ge. 2)  Then !-------------------------------+++--+
          ltx = Len_Trim(DBtxt(idb))                              !|
          write (txt,21) Element, Z, wave, idb, DBtxt(idb)(1:ltx) !|
  21      format('F1F2: Element=',a,' [Z=',i3,'] wave=',g12.6,
     *                              ' idb=',i1,' [',a,' DB]')     !|
          Call  Priblo  (txt,1,1,0,0,i)                           !|
c         txt(1)='Initial (internal) f1 and f2:'                  !|
c         Call  Priblo  (txt,1,0,0,0,i)                           !|
c         write (txt,16) (wave_n(i), wave2energy(wave_n(i)),      !|
c    *                    f1_n(i),f2_n(i),i=1,nval)               !|
c         Call  Priblo  (txt,1+nval,0,0,0,i)                      !|
        endif  !---------------------------------------------------+
c ------
        ifail = 0
        if (wave .le. 0.) Then !-----------------------------------+
          txt(1)='F1F2: exit on zero or negative wavelength'      !|
          if (ipr.ge.2) Call Priblo (txt,1,0,0,0,i)               !|
          ifail = -1                                              !|
          goto 999 !-----------------------------------------------+------+
        endif !----------------------------------------------------+      v
        EnergyEV = E2W / wave
        EKeV = EnergyEV / 1000.

c Check if the wavelength is in the tables range:
c All of these tables are for Z=1:92 or Z=3:92:
        if     (idb.eq.1) Then !---------------------------------------+ Henke
          if (EKeV.lt.0.01 .OR. EKeV.gt.30.) Then !-----------------+  |
            txt(1) = 'F1F2: X-rays not in Henke tables range'//    !|  |
     *               ' of 0.01-30 KeV (0.41-1240 A)'               !|  |
            if (ipr.ge.2) Call Priblo (txt,1,0,0,0,i)              !|  |
            ifail = -1                                             !|  |
            goto 999 !----------------------------------------------+--+---+
          endif !---------------------------------------------------+  |   v
          Zmin = 1                                                    !|
          Zmax = 92                                                   !|
          File = 'henke.dat'                                          !|
        elseif (idb.eq.2) Then !---------------------------------------+ Cowan
          if (EKeV.lt.0.03 .OR. EKeV.gt.700.) Then !----------------+  |
            txt(1) = 'F1F2: X-rays not in the Brennan-Cowan'//     !|  |
     *               ' tables range of 0.03-700 KeV (0.02-413 A)'  !|  |
            if (ipr.ge.2) Call Priblo (txt,1,0,0,0,i)              !|  |
            ifail = -1                                             !|  |
            goto 999 !----------------------------------------------+--+---+
          endif !---------------------------------------------------+  |   v
          Zmin = 3                                                    !|
          Zmax = 92                                                   !|
c         File = 'cowan.dat'                                          !|
          File = 'cowanlng.dat'                                       !|
        elseif (idb.eq.3) Then !---------------------------------------+ Windt
          if (EKeV.lt.0.01 .OR. EKeV.gt.100.) Then !----------------+  |
            txt(1) = 'F1F2: X-rays not in the Windt'//             !|  |
     *               ' tables range of 0.01-100 KeV (0.12-1240 A)' !|  |
            if (ipr.ge.2) Call Priblo (txt,1,0,0,0,i)              !|  |
            ifail = -1                                             !|  |
            goto 999 !----------------------------------------------+--+---+
          endif !---------------------------------------------------+  |   v
          Zmin = 1                                                    !|
          Zmax = 92                                                   !|
          File = 'windt.dat'                                          !|
        elseif (idb.eq.4) Then !---------------------------------------+ Chantler
          if (EKeV.lt.0.01 .OR. EKeV.gt.450.) Then !----------------+  |
            txt(1) = 'F1F2: X-rays not in the Chantler/NIST'//     !|  |
     *               ' tables range of 0.01-450 KeV (0.28-1240 A)' !|  |
            if (ipr.ge.2) Call Priblo (txt,1,0,0,0,i)              !|  |
            ifail = -1                                             !|  |
            goto 999 !----------------------------------------------+--+---+
          endif !---------------------------------------------------+  |   v
          Zmin = 1                                                    !|
          Zmax = 92                                                   !|
          File = 'chantler.dat'                                       !|
        else !---------------------------------------------------------+ ** incorrect DB **
          Zmin = 0                                                    !|
          Zmax = 999                                                  !|
          File = ' '                                                  !|
          goto 501 !---------------------------------------------error-+---+
        endif !--------------------------------------------------------+   v

c This is to detect ions:
        if     (Index(Element,'+').gt.0) Then !--+
          ishift = +1                           !|
        elseif (Index(Element,'-').gt.0) Then !--|
          ishift = -1                           !|
        else  !----------------------------------|
          ishift = 0                            !|
        endif  !---------------------------------+

        if (Z+ishift.lt.Zmin .or. Z+ishift.gt.Zmax) Then  !-------+
c No data for this element (do not overwrite internal f',f")      |
          write (txt(1),22) Z+ishift,                            !|
     *                      ' not in range of external DB:',     !|
     *                      Zmin, Zmax                           !|
  22      format (' F1F2: exit on Z=',i3,a,' [',i1,' - ',i3,']') !|
          if (ipr.ge.2) Call Priblo (txt,1,0,0,0,i)              !|
          ifail = -1                                             !|
          goto 999 !----------------------------------------------+-------+
        endif  !--------------------------------------------------+       v

        do j=1,nval !=====+
          wave_n(j) = 0. !|
          f1_n(j)   = 0. !|
          f2_n(j)   = 0. !|
        enddo  !==========+

c We are interested in the first 2 symbols only (thus, Li+ -> Li, F- -> F)
        El = Element(1:2)
        if (El(2:2).lt.'a' .OR. El(2:2).gt.'z') El(2:2)=' '
        j = Len_Trim (El)
        if (j .eq. 0)  goto 498  !--------------------error-+
                                                           !v
        lun = FreeLun ()    !get free available LUN
        j = Len_Trim(File)
        if (j.eq.0 .or. lun.le.0) goto 499  !---------error-+
        x0hpath(lx0hpat+1:80) = File                       !v
        lwrk = lx0hpat+j
        Call OpenFile(x0hpath(1:lwrk),lun,'read','old',
     *                io_status,*500)
c Now search for the element records:
c Examples of element header structure:   #S 1  H     #S  3 Li     #S  1  H
c                                         #S 92  U    #S  92 U     #S 92  U
        line = 0
  2     continue  !<--------------------------------------+
        line = line + 1                                  !|
        read  (lun,43,err=11,end=12)  txt(20)(1:16)      !^
  43    format(a)                                        !|
        if (txt(20)(1:3).ne.'#S ') goto 2  !--------->----|
        txt(20)(1:2) = ' '                               !|
        j = Index(txt(20)(1:16),El)                      !|
        if (j.eq.0)                goto 2  !--------->----|
        txt(20)(j:j+1) = ' '                             !|
        Call    rdInt (ipar,1,txt(20)(1:16),ifail)       !|
        if (ifail.ne.0) goto 11  !-----------------error--+----+
        iZ = ipar(1)                                     !|    v
        if (Z.ne.iZ .AND. Z+ishift.ne.iZ) goto 2  !-->----+
c We found the right element!
c Now skip the rest of comments:
        txt(20) = '#'
        do while (txt(20)(1:1).eq.'#')  !===============+
          line = line + 1                              !|
          read  (lun,43,err=11,end=12)  txt(20)(1:16)  !|
        enddo  !========================================+
c 'line' always contain the line number we have read,
c not the one we are going to read!
        lineStart = line

c Now find the closest energy in the datafile:
        Call rdReal (rpar,1,txt(20),ifail)
        if (ifail.ne.0)    goto 11  !--------------error-------+
        if (rpar(1).le.0.) goto 11  !--------------error-------+
        En = rpar(1)                                          !v
        Diff0 = Abs(En-EnergyEV)
        LineCenter = line
        do while (txt(20)(1:3)      .ne.'#S ' .and.   !These conditions
     *            Len_Trim(txt(20)) .gt. 0)  !--------+are verified not
          line = line + 1                            !|here, but below
          read  (lun,43,err=11,end=130) txt(20) !-->-+|     |
          if (txt(20)(1:3).eq.'#S ')  goto 13  !-->-+||   <=|
c The space between elements is present             |||     |
c in cowan.dat only!                                |||     |
          if (Len_Trim(txt(20)).eq.0) goto 13  !-->-|||   <=+
          Call rdReal (rpar,1,txt(20),ifail)       !|||
          if (ifail.ne.0)    goto 11  !-------error-+++------+
          if (rpar(1).le.0.) goto 11  !-------error-+++------+
          En = rpar(1)                             !|||      v
          Diff = Abs(En-EnergyEV)                  !|||
          if (Diff.gt.Diff0)   goto 13 !------>-----|||
          Diff0 = Diff                             !|||
          LineCenter = line                        !|||
        enddo  !------------------------------------+++
  130   continue  !<-------------------<------------++
c       line = line - 1                            !|
        Backspace(lun,err=11)                      !|
  13    continue  !<-------------------<------------+

c Here we are one line AFTER the center.
c Backspace by 'nval' DIFFERENT energy lines from the center:
        ii = 1
        En_saved = En
        do while (ii.lt.nval .AND. line.gt.lineStart) !===+
          line = line - 1                                !|
          Backspace(lun,err=11)                          !|
          line = line - 1                                !|
          Backspace(lun,err=11)                          !|
          line = line + 1                                !|
          read  (lun,43,err=11,end=11) txt(20) !-error----+------+
          if (txt(20)(1:3).eq.'#S ')  goto 17  !-->--+    |      v
          Call rdReal (rpar,1,txt(20),ifail)        !|    |
          if (ifail.ne.0)    goto 11  !---------error+----+------+
          if (rpar(1).le.0.) goto 11  !---------error+----+------+
          En = rpar(1)                              !|    |      v
          Diff = 2.*Abs((En-En_saved)/(En+En_saved))!|    |
          if (Diff.gt.Eps)  Then  !----+             |    |
            ii       = ii + 1         !|             |    |
            En_saved = En             !|             v    |
          endif  !---------------------+             |    |
        enddo  !=====================================+====+
        line = line - 1                             !|
        Backspace(lun,err=11)                       !|
  17    continue  !<---------------------<-----------+

c Start reading the data:
c ----------------------
c (first, we try to read 'nval' BEFORE the center. Then,
c  we try to replace a half of them by the values AFTER
c  the center):
        ii        = 0
        Wv_saved  = 0.
        i_repeat  = 1
        ii_center = 1
        do  while (ii.le.nval) !=========================+
          line = line + 1                               !|
          read  (lun,43,err=11,end=11) txt(20)          !|
          Call rdReal (rpar,3,txt(20),ifail)            !|
          if (ifail.ne.0)     goto 11  !----error--------+--+
          if (rpar(1).le.0.)  goto 11  !----error--------+--|
c Special values (see Database file documentation):      |  v
          if (abs(rpar(2)-nodata).lt.1.E-20) Then !--+   |
            rpar(2) = 0.                            !|   |
          else  !------------------------------------+   |
c All tables except for Cowan                        |   |
c contain f+df1 instead of df1:                      |   |
            if (idb.ne.2) rpar(2)=rpar(2)-iZ        !|   |
          endif  !-----------------------------------+   |
c Special values (see Database file documentation):      |
          if (abs(rpar(3)-nodata).lt.1.E-20) Then !--+   |
            rpar(3) = 0.                            !|   |
          endif  !-----------------------------------+   |
          Wv = E2W / rpar(1)                            !|
          Diff = 2.*Abs((Wv-Wv_saved)/(Wv+Wv_saved))    !|
          if (Diff .gt. Eps)  Then  !------------+       |
c Normalize the previous result for the number   |       |
c of repetitions for the previous wave:          |       |
            if (ii.gt.0)  Then  !--------------+ |       |
              f1_n(ii) = f1_n(ii) / i_repeat  !| |       |
              f2_n(ii) = f2_n(ii) / i_repeat  !| |       |
            endif  !---------------------------+ |       |
c Start new count of repetitions:                |       |
            i_repeat = 1                        !|       |
            Wv_saved = Wv                       !|       |
            ii = ii + 1                         !|       |
            wave_n(ii) = Sngl(Wv)               !|       |
            f1_n(ii)   = rpar(2)                !|       |
            f2_n(ii)   = rpar(3)                !|       |
            if (ii.ge.nval)  goto 18  !----->----+--->---+---+
          else  !--------------------------------|       |   |
            i_repeat   = i_repeat + 1           !|       |   |
            f1_n(ii)   = f1_n(ii) + rpar(2)     !|       |   |
            f2_n(ii)   = f2_n(ii) + rpar(3)     !|       |   |
          endif !--------------------------------+       |   |
          if (line.eq.lineCenter) ii_center=ii          !|   |
        enddo  !=========================================+   |
  18    continue  !<---------------------------------<-------+

        nmore = nval/2 + (ii_center - nval)
        ii = 0
        do while (ii.le.nmore)   !=======================+
          line = line + 1                               !|
          read  (lun,43,err=11,end=14) txt(20) !----->---+--->-+
          if (txt(20)(1:3).eq.'#S ') goto 14   !----->---+--->-|
          Call rdReal (rpar,3,txt(20),ifail)            !|     |
          if (ifail.ne.0)     goto 11  !-----------error-+-+   |
          if (rpar(1).le.0.)  goto 11  !-----------error-+-|   |
c Special values (see Database file documentation):      | v   |
          if (abs(rpar(2)-nodata).lt.1.E-20) Then !--+   |     |
            rpar(2) = 0.                            !|   |     |
          else  !------------------------------------+   |     |
c All tables except for Cowan                        |   |     |
c contain f+df1 instead of df1:                      |   |     |
            if (idb.ne.2) rpar(2)=rpar(2)-iZ        !|   |     |
          endif  !-----------------------------------+   |     |
          if (abs(rpar(3)-nodata).lt.1.E-20) Then !--+   |     |
            rpar(3) = 0.                            !|   |     |
          endif  !-----------------------------------+   |     |
          Wv = E2W / rpar(1)                            !|     |
          Diff = 2.*Abs((Wv-Wv_saved)/(Wv+Wv_saved))    !|     |
          if (Diff .gt. Eps)  Then  !--------------+     |     |
            f1_n(nval) = f1_n(nval) / i_repeat    !|     |     |
            f2_n(nval) = f2_n(nval) / i_repeat    !|     |     |
            i_repeat   = 1                        !|     |     |
            Wv_saved   = Wv                       !|     |     |
            ii         = ii + 1                   !|     |     |
            if (ii.gt.nmore)  goto 14  !----->-----+-----+-----|
            do j=2,nval  !==============+          |     |     |
              wave_n(j-1) = wave_n(j)  !|          |     |     |
              f1_n(j-1)   = f1_n(j)    !|          |     |     |
              f2_n(j-1)   = f2_n(j)    !|          |     |     |
            enddo  !====================+          |     |     |
            wave_n(nval) = Sngl(Wv)               !|     |     |
            f1_n(nval)   = rpar(2)                !|     |     |
            f2_n(nval)   = rpar(3)                !|     |     |
          else !-----------------------------------|     |     |
            i_repeat     = i_repeat + 1           !|     |     |
            f1_n(nval)   = f1_n(nval) + rpar(2)   !|     |     |
            f2_n(nval)   = f2_n(nval) + rpar(3)   !|     |     |
          endif !----------------------------------+     |     |
        enddo  !=========================================+     |
  14    continue  !<---------------<-------------------------<-+
        f1_n(nval) = f1_n(nval) / i_repeat
        f2_n(nval) = f2_n(nval) / i_repeat

c Check again for coinciding waves:
        do i=2,nval  !================================================+
          diff = Abs((wave_n(i)-wave_n(i-1))/(wave_n(i)+wave_n(i-1)))!|
          if (diff .lt. Eps) goto 19  !-------------------------------+--+
        enddo  !======================================================+  V

        if (ipr .ge. 2)  Then !---------------------------------------+
          j = Len_Trim(File)                                         !|
          write (txt,15) Element, iZ, File(1:j)                      !|
  15      format('Element: ',a,' [',i3,'],  df1,df2 data from: ',a)  !|
          Call  Priblo  (txt,1,1,0,0,i)                              !|
          write (txt,16) (wave_n(i), wave2energy(wave_n(i)),         !|
     *                    f1_n(i),f2_n(i),i=1,nval)                  !|
  16      format(' - Wavelength:     Energy:       df1:',10x,'df2:'/
     *           21(4g14.6,:/))                                      !|
          Call  Priblo  (txt,1+nval,0,0,0,i)                         !|
        endif  !------------------------------------------------------+
        goto 998  !-------------------------------------exit--+
c ........................                                    V
c     E R R O R S :
c ........................
  498   continue
        write   (txt,798)       x0hpath(1:lwrk)
  798   format  (' "',a,'"'//' - no element specified')
        ifail = 3
        goto 999
c ........................
  499   continue
        write   (txt,799)       x0hpath(1:lwrk)
  799   format  (' "',a,'"'//' - no filename given or lun=0')
        ifail = 3
        goto 999
c ........................
  500   continue
        write   (txt,800)       x0hpath(1:lwrk)
  800   format  (' "',a,'"'//' - open error')
        ifail = 3
        goto 999
c ........................
  501   continue
        write   (txt,801)       idb
  801   format  (' Incorrect Henke/Cowan/Windt/Chantler flag idb =',i2)
        ifail = 1
        goto 999
c ........................
  11    continue
        write   (txt,900)       x0hpath(1:lwrk), line, El
  900   format  (' "',a,'"'//' - read error',
     *  ' at line=',i5,',  element [',a,']')
        ifail = 3
        goto 998
c ........................
  12    continue
        write   (txt,901)       x0hpath(1:lwrk), El
  901   format  (' "',a,'"'//' - element [',a,'] not found')
        ifail = 3
        goto 998
c ........................
  19    continue
        write   (txt,902)       x0hpath(1:lwrk), El
  902   format  (' "',a,'"'//' - coinciding waves for element [',a,']'/
     *           'in the table of dispersion corrections')
        ifail = 4
        goto 998

  998   Close   (unit=lun,err=999)  !<---------------------+
  999   x0hpath(lx0hpat+1:80)=' '
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
