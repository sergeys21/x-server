        Program x0h_web
c **********************************************************************
c * This program reads data from command line in the key=data format   *
c * with the keys from the X0h request web page at X-ray Server and    *
c * calculates X-ray structure factor for given wavelength, crystal    *
c * and Bragg indices. The results are printed to STDOUT as a web page.*
c *                                                                    *
c *                  Author: Sergey Stepanov                           *
c * 2026/01 Version: 1.0:   initial implementation                     *
c **********************************************************************
c The include file contains:
c       Parameter       (kcompMax =  10)
c       Parameter       (kolMax   =  32)
c       Parameter       (ncodMax  = 200)
c       Parameter       (natmMax  = 100)
c       Parameter       (nwavMax  =  50)
c       Parameter       (nedgMax  =  16)
c       Parameter       (nabcMax  =  11)
c       Parameter       (maxX0hWarnLines = 10);
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Real*8    w2e
c https://www.physics.nist.gov/MajResFac/SURF/SURF/schwingersample.html
        Parameter (w2e=12.3984186)      !E=w2e/wave, NIST web site

        Complex*8 xh
        Real*8    qb8, d8, ESinTheta
        Real*4    Lattice_constants(6), wave, energy, x0_abs,
     *            xr0, xi0, xrh, xih, xqh, xrf, xif, xqf, xdf,
     *            mu0, a0, rho, pi, gra, sec, s6, qb,
     *            sin1, cos1, tg1, ctg1, sin2, cos2,
     *            xrsg, xisg, xdfsg, extisg, extisg_l, halfsg, prcsg,
     *            xrpi, xipi, xdfpi, extipi, extipi_l, halfpi, prcpi,
     *            xhsg, xhpi, TER,dl_TER, w, rangmin, rangmax,
     *            ah, ah1sg, ah2sg, ah1pi, ah2pi,
     *            xh_weak_limit,                !the weakest possible xh
     *            shift, sqr2
        Integer   indices(3), iii(4),
     *            isyngony, i, j, l,
     *            icoway, ixway,
     *            ifail, mfail, modeout, mdetail,
     *            nHenkeCowan, npass, ipass, lout,
     *            iref, jpr, inmax, icall, iwarn,
     *            n0, n1, n2, n3, n9, indx, io_status,
     *            iunis, npts, isubway, ipol,
     *            lunout, lundup, lhost, laddr, lurl,
     *            lsp, lq0, lqh, ldb
        Character code*20, txcount*8, pl*2,
     *            quest0*128, questh*128, url*64, head*6,
     *            outfile*32, query*256, parser*8,
     *            hostname*256, address*60, referer*160,
     *            Errors(13)*56 /
     *        ' ',                                                       ! ifail=0
     *        'No material name specified',                              ! ifail=1
     *        'Unexpected type of x-ray input',                          ! ifail=2
     *        'Wrong or missing x-ray wavelength/energy',                ! ifail=3
     *        'Missing name of characteristic x-ray line',               ! ifail=4
     *        'Wrong or missing Bragg reflection indices (hkl)',         ! ifail=5
     *        'No x-ray data specified',                                 ! ifail=6
     *        'Negative X-ray wavelength/energy',                        ! ifail=7
     *        'Requested reflection does not exist for this wavelength', ! ifail=8
     *        'Incorrect or not specified material density',             ! ifail=9
     *        'Unexpected type DB option',                               ! ifail=10
     *        'Ambiguous material specification',                        ! ifail=11
     *        'Incorrect X0h input parameter'/,                          ! ifail=12
     *            Symmetry(9)*24 /
     *        'Unknown',                                                 !0
     *        'Cubic',                                                   !1
     *        'Hexagonal/Trigonal',                                      !2
     *        'Tetragonal',                                              !3
     *        'Trigonal(Rhombohedral)',                                  !4
     *        'Orthorhombic',                                            !5
     *        'Monoclinic',                                              !6
     *        'Triclinic',                                               !7
     *        '** Amorphous ** '/                                        !8
c X0h (International Tables) range:
c             (5,000 eV --  25,000 eV):   0.50  --    2.47 Angstrem
c Henke    range (10 eV --  30,000 eV):   0.41  -- 1239.81 Angstrem
c Cowan    range (30 eV -- 694,500 eV):   0.02  --  413.27 Angstrem
c Windt    range (10 eV -- 100,000 eV):   0.12  -- 1239.81 Angstrem
c Chantler range (10 eV -- 450,000 eV):   0.28  -- 1239.81 Angstrem
        Character DBtext(10)*45 /
     *        'Automatic DB selection',                                  !-1 (1)
     *        'X0h (International Tables), 5-25 KeV',                    !0  (2)
     *        'Henke Tables, 0.01-30 KeV (df1 only)',                    !1  (3)
     *        'Henke Tables, 0.01-30 KeV',                               !2  (4)
     *        'Brennan-Cowan Tables, 0.03-700 KeV (df1 only)',           !3  (5)
     *        'Brennan-Cowan Tables, 0.03-700 KeV',                      !4  (6)
     *        'Windt Tables, 0.01-100 KeV (df1 only)',                   !5  (7)
     *        'Windt Tables, 0.01-100 KeV',                              !6  (8)
     *        'Chantler/NIST Tables, 0.01-450 KeV (df1 only)',           !7  (9)
     *        'Chantler/NIST Tables, 0.01-450 KeV'/                      !8  (10)

        Integer   Nkeys
        Parameter (Nkeys=18)
        Character keys(Nkeys)*8 /'xway', 'wave', 'line',
     *                           'code', 'coway', 'amor', 'chem', 'rho',
     *                           'i1', 'i2', 'i3', 'df1df2',
     *                           'detail', 'modeout',
     *                           'job', 'ip', 'host', 'referrer'/,
     *            vals(Nkeys)*256

        Character tnb*14 /' class="novrb"'/,
     *            ttb*14 /' class="nobtm"'/,
c    *            tbb*14 /' class="notop"'/,
     *            br*4, spacer*24

        Real*8    Bragan8
        External  Bragan8

        Real      round, wave2energy
        External  round, wave2energy

        Character Get_Server_URL*64
        External  Get_Server_URL

        Character radiat*6
        Common /x0pa1/ radiat

        Character name(kcompMax)*4
        Common /x0pa2/ name

        Real*4 prcn(kcompMax)
        Common /x0pa3/ prcn

        Integer kcomp,kol(kcompMax)
        Common /x0pa6/ kcomp,kol

        Real            wx(kcompMax,kolMax),
     *                  wy(kcompMax,kolMax),
     *                  wz(kcompMax,kolMax)

        Integer         isym
        Common  /x0pa7/ wx, wy, wz, isym

        Real*4 Poisson, edge_nearest
        Character name_nearest*4
        Common /x0pa9/ Poisson, edge_nearest, name_nearest

        Real            ucell_mass_gram, atom_density_cm3
        Integer         n_atoms_ucell
        Common  /x0paA/ ucell_mass_gram, atom_density_cm3, n_atoms_ucell

        Integer         linesX0hWarn, linesX0hWarn_
        Character       txtX0hWarn(maxX0hWarnLines)*80,
     *                  txtX0hWarn_(maxX0hWarnLines)*80
        Common /x0hwarn/ txtX0hWarn, linesX0hWarn

        Character txt(20)*80
        Common /msg/ txt
c-----------------------------------------------------------------------
c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer             iHenkeCowan
        Common /HenkeCowan/ iHenkeCowan
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

        istackrezv = 0
        progname   = 'x0h_web'
        ErrorFile  = 'x0h_web.err'
        modebat    = 1
        lunbat     = 0
        Call    DeleteErrorFile ()
c -------------------------------------------------------
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
        pi       = 4.*atan(1.)          !the "pi" constant
        gra      = 2.*pi/360.           !number of radians in a degree
        sec      = gra/3600.            !number of radians in arc sec
        s6       = sec*(1.E+6)
        sqr2     = sqrt(2.)
        TER      = 0.
        dl_TER   = 0.
        mu0      = 0.
        a0       = 0.
        ah       = 0.
        ah1sg    = 0.
        ah2sg    = 0.
        ah1pi    = 0.
        ah2pi    = 0.
        d8       = 0.0

        quest0 = '<a href="javascript:void(0)" onClick="Wfloat'//
     *           '(''images/x0h_help_0.gif'',''x0h_0'',740,357);">'//
     *           '<b>?</b></a> &nbsp;'
        questh = '<a href="javascript:void(0)" onClick="Wfloat'//
     *           '(''images/x0h_help_h.gif'',''x0h_h'',705,853);">'//
     *           '<b>?</b></a> &nbsp;'
        lq0 = Len_Trim(quest0)
        lqh = Len_Trim(questh)

        outfile  = '/dev/stdout'
        lout = Len_Trim(outfile)
        lunout   = 6
        lundup   = lunout + 1
        ifail = 0

        Call Get_CmdLine (nkeys,keys,vals,ifail)

        Call getbykey(nkeys,keys,vals,'host',hostname)
        Call getbykey(nkeys,keys,vals,'ip',address)
        Call getbykey(nkeys,keys,vals,'referrer',referer)
        Call getbykey(nkeys,keys,vals,'job',txcount)

        laddr = Len_Trim(address)
        if (laddr.eq.0) Then !--------+
          address = '127.0.0.1'      !|
          hostname = 'localhost'     !|
          laddr = Len_Trim(address)  !|
        endif  !----------------------+
        lhost = Len_Trim (hostname)
        if (lhost.eq.0) Then !--+
          hostname = address   !|
          lhost    = laddr     !|
        endif  !----------------+
        if (Len_Trim(txcount).eq.0) txcount = 'x0h00000'

c If call for Get_CmdLine above failed:
        if (ifail .ne. 0) Then !----------------------------+
          do i=1,nkeys !=================================+  |
            l = Max(Len_Trim(vals(i)),1)                !|  |
            write(0,'(1x,3a)') keys(i),'=',vals(i)(1:l) !|  |
          enddo !========================================+  |
          txt(1) = 'Error parsing cmdline'                 !|
          Call HTMLstop (lunout, 1)                        !|
          stop 'x0h_web: emergency exit'                   !|
        endif !---------------------------------------------+

        Call OpenFile(txcount//'.htm',lundup,
     *                'readwrite','unknown',
     *                io_status,*4)

        write (lunout,511,err=14)
c 511   format(' Content-type: text/html'/)     !before 2013.11.01
  511   format('Content-type: text/html'/)
        head = '<head>'
        write (lunout,'(a//a/a)') '<!DOCTYPE html>','<html>',head
        Call Insert_UrlBase(head,lunout)                                !see
        write (lunout,512,err=14)
  512   format( ' <title>X0h Results</title>'/
     *  ' <link rel="stylesheet" href="style.css" type="text/css">'/
     *  ' <script language="javascript" type="text/javascript">'/
     *  ' <!-- hide'/
     *  ' function Wfloat(webpage,windowname,pagewidth,pageheight) {'/
     *  ' launch = window.open (webpage,windowname,',
     *    '"width="+pagewidth+",height="+pageheight+",',
     *    'toolbar=0,location=0,directories=0,status=0,',
     *    'menubar=0,scrollbars=1,resizable=1");'/
     *  ' }'/
     *  ' // -->'/
     *  ' </script>'/
     *  ' </head><body>')

        Call mergeQuery(nkeys-4,keys,vals,query)
        j = Max(Len_Trim(query),1)
        url = Get_Server_URL()
        lurl = Max(Len_Trim(url),1)
        write (lundup,513,err=14) url(1:lurl),'cgi/x0h_form.pl',
     *                            query(1:j)
        Call flushoutput(lundup)
  513   format('<!DOCTYPE html>'//'<html>'/
     *  ' <head><title>X0h Results</title></head>'/
     *  ' <body text="#000000" bgcolor="silver">'/
     *  ' <!-- ',2a,'?',a,' -->'/
     *  ' <table border="0"><tr valign="top">')
c-------
        write (lunout,50,err=14)
  50    format(' <table border="0"><tr valign="top"><td>'/
     *         ' &nbsp; <img src="images/x0h_web_th_th.jpg"',1x,
     *                  'align="absmiddle" border="0">&nbsp;</td>')
c-------
        if (hostname .ne. address)  Then  !---------------------+
          do i=lunout,lundup !================================+ |
            write (i,51,err=14) txcount,                     !| |
     *                          ' for ',address(1:laddr),    !| |
     *                          ' [',hostname(1:lhost),']:'  !| |
            Call flushoutput(i)                              !| |
          enddo !=============================================+ |
        else  !-------------------------------------------------+
          do i=lunout,lundup !================================+ |
            write (i,51,err=14) txcount,                     !| |
     *                          ' for ',address(1:laddr),    !| |
     *                          ' :',' ',' '                 !| |
            Call flushoutput(i)                              !| |
          enddo !=============================================+ |
        endif  !------------------------------------------------+
  51    format(' <td style="white-space:nowrap;"><font size="+1">'/
     *         ' &nbsp; Job ID: &nbsp;<b>',a,'</b><br>'/
     *         ' &nbsp; <b><i>X0h</i> Results',5a,'</b>'/
     *         ' </font></td></tr></table><br>')
        if (lhost.eq.0)  hostname=address

        do i=lunout,lundup !=======+
          write (i,52,err=14)     !|
          Call flushoutput(i)     !|
        enddo !====================+
c Allowed colors:
c --------------
c Black, Maroon, Green, Olive, Navy, Purple, Teal, Gray,
c Silver,Red, Lime, Yellow, Blue, Fuchsia, Aqua, White
c
c <table border bgcolor="#d0f0f0"> -- light blue color!
c <table border bgcolor="#f0f0f0"> -- not bright white color!

  52    format(' <center><table border bgcolor="#ffffff">')

        parser = 'keys2x0h'
        write (lundup,*,err=14) '<!-- Calling ',parser,' -->'
        Call flushoutput(lundup)
        Call keys2x0h (nkeys-4, keys, vals,
     *                 icoway, code, rho,
     *                 ixway, radiat, wave, energy,
     *                 indices, iref,
     *                 nHenkeCowan, modeout, mdetail, ifail)
        if (ifail.ne.0)  Then  !--------------------------------+
          mfail = 3                                            !|
            txt(1) = 'E R R O R:'                              !|
            txt(2) = ' '                                       !|
          if (ifail.le.12) Then !-----------------------------+ |
            txt(3) = Errors(ifail+1)                         !| |
          else !----------------------------------------------+ |
            write(txt(3),'(3a,i3)') 'Unexpected error from ',!| |
     *                              parser,' err=',ifail     !| |
          endif !---------------------------------------------+ |
          Call HTMLstop (lunout,mfail)                         !|
          Call HTMLstop_noJPG (lundup,mfail)                   !|
          goto 999                                             !|
        endif  !------------------------------------------------+
c       write (lundup,*,err=14) '<!- ixway=',ixway,' wave=',wave,
c    *                                    ' radiat=',radiat,' -->'
        write (lundup,*,err=14) '<!-- Done with ',parser,' -->'
        Call flushoutput(lundup)

        if (nHenkeCowan.ne.10) Then  !---+
          iHenkeCowan = nHenkeCowan     !|
          npass = 1                     !|
        else  !--------------------------+
          iHenkeCowan = 0               !|
          npass = 5                     !|
        endif  !-------------------------+

        do      i=1,6  !===============+
          Lattice_Constants(i) = 0.   !|
        enddo  !=======================+

        if (modeout.eq.0) Then !---+
          spacer = '<br>'         !|
        else  !--------------------!
          spacer = '</pre>'       !|
        endif !--------------------+
        lsp = Max(Len_Trim(spacer),1)

        do ipass =1,npass  !===============================================+
          if (icoway.eq.0 .or. icoway.eq.-1) Then !---+                    |
            isyngony  = 0       ! no syngony control  |                    |
          else  !-------------------------------------+                    |
            isyngony  = 8       ! must be amorphous  !|                    |
          endif  !------------------------------------+                    |
          icall        = 0                                                !|
          jpr          = 0                                                !|
c         jpr          = 2                                                !|
          qb           = 0.                                               !|
          xr0          = 0.                                               !|
          xi0          = 0.                                               !|
          x0_abs       = 0.                                               !|
          xrh          = 0.                                               !|
          xih          = 0.                                               !|
          xqh          = 0.                                               !|
          xrf          = 0.                                               !|
          xif          = 0.                                               !|
          xqf          = 0.                                               !|
          xdf          = 0.                                               !|
          ucell_mass_gram  = 0.                                           !|
          atom_density_cm3 = 0.                                           !|
          n_atoms_ucell    = 0                                            !|
          do i=1,3 !=============+                                         |
            iii(i) = indices(i) !|                                         |
          enddo  !===============+                                         |
          inmax        = 3                                                !|
          edge_nearest = 0.                                               !|
          name_nearest =' '                                               !|
          do      i=1,kcompMax !==+                                        |
            prcn(i) = 0.         !|                                        |
          enddo  !================+                                        |
          linesX0hWarn_ = 0                                               !|
c ifail=5 Wrong or missing Bragg reflection indices (hkl),                 |
c ifail=8 Requested reflection does not exist for this wavelength,         |
c ifail=9 No density for chemical formula (but it can be atom/crystal!)    |
          if (ifail.eq.0  .OR.                                            !|
     *        ifail.eq.5  .OR.                                            !|
     *        ifail.eq.8  .OR.                                            !|
     *        ifail.eq.9)     Then  !-------------------------------+      |
c Calculate X0:                                                     |      |
            write (lundup,*,err=14) '<!-- Calling X0h1 (000) -->'  !|      |
            Call flushoutput(lundup)                               !|      |
            Call X0h1 (Lattice_constants, wave,                    !|      |
     *                 code, rho, 0,                               !|      |
     *                 0, 0, 0,                                    !|      |
     *                 isyngony, icall,                            !|      |
     *                 jpr, qb, xr0, xi0,                          !|      |
     *                 xrh, xih, xqh,                              !|      |
     *                 xrf, xif, xqf, mfail)                       !|      |
            write (lundup,*,err=14) '<!-- Done X0h1 (000) -->'     !|      |
            Call flushoutput(lundup)                               !|      |
            x0_abs = xr0*xr0+xi0*xi0                               !|      |
            if (x0_abs.gt.0.) Then  !---+                           |      |
              x0_abs = sqrt(x0_abs)    !|                           |      |
            endif  !--------------------+                           |      |
c Critical angle of total reflection:                               |      |
            if (x0_abs.gt.0.) Then !---------+                      |      |
              TER = Sqrt(x0_abs)            !|                      |      |
c Penetration depth in TER region:           |                      |      |
              dl_TER = wave / (2.*pi*TER)   !|                      |      |
            else  !--------------------------+                      |      |
              TER    = 0.                   !|                      |      |
              dl_TER = 0.                   !|                      |      |
            endif  !-------------------------+                      |      |
            if (abs(xi0).gt.1.E-33) Then !---------+                |      |
c Absorption factor(1/cm):                         |                |      |
c mu0 = k0*|xi0|                                   |                |      |
              mu0 = (2.E+08)*pi*abs(xi0)/wave     !|                |      |
c Absorption length(micron):                       |                |      |
c mu0 = 1/(k0*|xi0|)                               |                |      |
              a0 = (1.E+04)/mu0                   !|                |      |
c             a0 = (1.E-04)*wave/(2.*pi*abs(xi0)) !|                |      |
            else  !--------------------------------+                |      |
              mu0    = 0.                         !|                |      |
              a0     = 0.                         !|                |      |
            endif  !-------------------------------+                |      |
            if (mfail.ne.0)  Then  !-----------------+X0h1 failed!  |      |
              mfail = Iabs(mfail)                   !|              |      |
              Call HTMLstop       (lunout,mfail)    !|              |      |
              Call HTMLstop_noJPG (lundup,mfail)    !|              |      |
              goto 999                              !|              |      |
            endif  !---------------------------------+              |      |
            if (isyngony.eq.2)  Then  !----------+                  |      |
              iii(3) = -(indices(1)+indices(2)) !|                  |      |
              iii(4) = indices(3)               !|                  |      |
              inmax  = 4                        !|                  |      |
            endif  !-----------------------------+                  |      |
            if (isyngony.eq.8)  Then  !----+                        |      |
              iref = 0                    !|                        |      |
            endif  !-----------------------+                        |      |
            if (linesX0hWarn .gt. 0) Then !---------+               |      |
              linesX0hWarn_ = linesX0hWarn         !|               |      |
              do i=1,linesX0hWarn !==============+  |               |      |
                txtX0hWarn_(i) = txtX0hWarn(i)  !|  |               |      |
              enddo !============================+  |               |      |
              linesX0hWarn = 0       !reset counter |               |      |
            endif !---------------------------------+               |      |
                                                                   !|      |
            if (iref .ne. 0  .AND.                                 !|      |
     *          ifail .eq. 0)   Then  !-----------------------+     |      |
              iwarn = 0         !take QB=0 if l>2d           !|     |      |
              write (lundup,*,err=14) '<!-- Calling',        !|     |      |
     *                                   ' Bragan8 -->'      !|     |      |
              Call flushoutput(lundup)                       !|     |      |
              qb8 = Bragan8 (wave,iii,inmax,iwarn,d8,1,      !|     |      |
     *                     sin1,cos1,tg1,ctg1,               !|     |      |
     *                     sin2,cos2)                        !|     |      |
              write (lundup,*,err=14) '<!-- Done',           !|     |      |
     *                                   ' Bragan8 -->'      !|     |      |
              Call flushoutput(lundup)                       !|     |      |
              qb = SNGL(qb8)                                 !|     |      |
              if (abs(qb8).lt.1.E-10) Then !---------------+  |     |      |
                ifail = 8     !no Bragg angle for hkl      |  |     |      |
                ESinTheta = 0.0D0                         !|  |     |      |
              else  !--------------------------------------+  |     |      |
c Calculate Xh:                                            |  |     |      |
                ESinTheta = w2e/(2.*d8)                   !|  |     |      |
                icall = 5                                 !|  |     |      |
                write (lundup,*,err=14) '<!-- Calling',   !|  |     |      |
     *                            ' X0h1 (hkl) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
                Call X0h1(Lattice_constants, wave,        !|  |     |      |
     *                   code, rho, iref,                 !|  |     |      |
     *                   iii(1), iii(2), iii(inmax),      !|  |     |      |
     *                   isyngony, icall,                 !|  |     |      |
     *                   jpr, qb, xr0, xi0,               !|  |     |      |
     *                   xrh, xih, xqh,                   !|  |     |      |
     *                   xrf, xif, xqf, mfail)            !|  |     |      |
                write (lundup,*,err=14) '<!-- Done',      !|  |     |      |
     *                            ' X0h1 (hkl) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
                if (mfail.ne.0)  Then  !--------------+    |  |     |      |
                  mfail = Iabs(mfail)                !|    |  |     |      |
                  Call HTMLstop (lunout,mfail)       !|    |  |     |      |
                  Call HTMLstop_noJPG (lundup,mfail) !|    |  |     |      |
                  goto 999                           !|    |  |     |      |
                endif  !------------------------------+    |  |     |      |
                ah = a0*cos1        !this is 0 for QB=90  !|  |     |      |
                xdf = xrf - xif                           !|  |     |      |
                write (lundup,*,err=14) '<!-- Calling',   !|  |     |      |
     *                           ' SigmaPi(sg) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
                Call Sigmapi (1, qb, wave, x0_abs,        !|  |     |      |
     *                   xrh, xih, xqh,                   !|  |     |      |
     *                   xrf, xif, xqf,                   !|  |     |      |
     *                   xrsg, xisg, xdfsg,               !|  |     |      |
     *                   xhsg,         !added 2001/10      |  |     |      |
     *                   extisg, halfsg, prcsg)           !|  |     |      |
                write (lundup,*,err=14) '<!-- Done',      !|  |     |      |
     *                           ' SigmaPi(sg) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
                                                          !|  |     |      |
                write (lundup,*,err=14) '<!-- Calling',   !|  |     |      |
     *                          ' Borrmann(sg) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
c See the code below in this file:                         |  |     |      |
                Call Borrmann (ah, xi0, xisg, xdfsg,      !|  |     |      |
     *                         ah1sg, ah2sg)              !|  |     |      |
                write (lundup,*,err=14) '<!-- Done',      !|  |     |      |
     *                          ' Borrmann(sg) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
                                                          !|  |     |      |
                write (lundup,*,err=14) '<!-- Calling',   !|  |     |      |
     *                           ' SigmaPi(pi) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
                Call Sigmapi (2, qb, wave, x0_abs,        !|  |     |      |
     *                   xrh, xih, xqh,                   !|  |     |      |
     *                   xrf, xif, xqf,                   !|  |     |      |
     *                   xrpi, xipi, xdfpi,               !|  |     |      |
     *                   xhpi,         !added 2001/10      |  |     |      |
     *                   extipi, halfpi, prcpi)           !|  |     |      |
                write (lundup,*,err=14) '<!-- Done',      !|  |     |      |
     *                           ' SigmaPi(pi) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
                                                          !|  |     |      |
                write (lundup,*,err=14) '<!-- Calling',   !|  |     |      |
     *                          ' Borrmann(pi) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
c See the code below in this file:                         |  |     |      |
                Call Borrmann (ah, xi0, xipi, xdfpi,      !|  |     |      |
     *                         ah1pi, ah2pi)              !|  |     |      |
                write (lundup,*,err=14) '<!-- Done',      !|  |     |      |
     *                          ' Borrmann(pi) -->'       !|  |     |      |
                Call flushoutput(lundup)                  !|  |     |      |
c Laue-case extinction lengths                            !|  |     |      |
c (differ by the factor of pi from the Bragg               |  |     |      |
c  case definition of the extinction length):              |  |     |      |
                if (xhsg .gt. 0.) Then !-------------+     |  |     |      |
                  extisg_l = 0.0001*Wave*cos1/xhsg  !|     |  |     |      |
                else !-------------------------------+     |  |     |      |
                  extisg_l = 0.                     !|     |  |     |      |
                endif !------------------------------+     |  |     |      |
                if (xhpi .gt. 0.) Then !-------------+     |  |     |      |
                  extipi_l = 0.0001*Wave*cos1/xhpi  !|     |  |     |      |
                else !-------------------------------+     |  |     |      |
                  extipi_l = 0.                     !|     |  |     |      |
                endif !------------------------------+     |  |     |      |
              endif  !-------------------------------------+  |     |      |
            endif  !------------------------------------------+     |      |
                                                                   !|      |
            energy = wave2energy(wave)                             !|      |
                                                                   !|      |
          elseif (ifail.eq.3      !wrong/missing wave/energy        |      |
     *       .or. ifail.eq.4      !missing x-ray line               |      |
     *       .or. ifail.eq.6      !no x-ray data specified          |      |
c    *       .or. ifail.eq.7      !Wave not in [0.1-10A]            |      |
     *                      )  Then !-------------------------------+      |
c No X-ray data: call X0h to get lattice constants only:            |      |
            w = 1.      !(wave=1)                                  !|      |
            write (lundup,*,err=14) '<!-- Calling X0h1 (abc) -->'  !|      |
            Call flushoutput(lundup)                               !|      |
            Call X0h1 (Lattice_constants, w,                       !|      |
     *                 code, rho, 0,                               !|      |
     *                 0, 0, 0,                                    !|      |
     *                 isyngony, icall,                            !|      |
     *                 jpr, qb, xr0, xi0,                          !|      |
     *                 xrh, xih, xqh,                              !|      |
     *                 xrf, xif, xqf, mfail)                       !|      |
            write (lundup,*,err=14) '<!-- Done X0h1 (abc) -->'     !|      |
            Call flushoutput(lundup)                               !|      |
            if (mfail.ne.0)  Then  !--------------+                 |      |
              mfail = Iabs(mfail)                !|                 |      |
              Call HTMLstop (lunout,mfail)       !|                 |      |
              Call HTMLstop_noJPG (lundup,mfail) !|                 |      |
              goto 999                           !|                 |      |
            endif  !------------------------------+                 |      |
            xr0    = 0.                                            !|      |
            xi0    = 0.                                            !|      |
            x0_abs = 0.                                            !|      |
            if (isyngony.eq.2)  Then  !----------+                  |      |
              iii(3) = -(indices(1)+indices(2)) !|                  |      |
              iii(4) = indices(3)               !|                  |      |
              inmax  = 4                        !|                  |      |
            endif  !-----------------------------+                  |      |
          endif  !--------------------------------------------------+      |
                                                                          !|
          n0 = Len_Trim(code)                                             !|
          if (n0.eq.0)  Then  !-------------------------------------+      |
            code = 'none'                                          !|      |
            n0 = Len_Trim(code)                                    !|      |
            if     (icoway.eq.0) Then !--------------------------+  |      |
              Errors(1+1)='No crystal name specified'           !|  |      |
            elseif (icoway.eq.1) Then !--------------------------+  |      |
              Errors(1+1)='No amorphous material name specified'!|  |      |
            elseif (icoway.eq.2) Then !--------------------------+  |      |
              Errors(1+1)='No chemical formula specified'       !|  |      |
            endif  !---------------------------------------------+  |      |
          endif  !--------------------------------------------------+      |
                                                                          !|
          if (abs(wave).lt.1.E-33) Then  !-----+                           |
            txt(1) = 'not specified'          !|                           |
            n1 = Len_Trim(txt(1))             !|                           |
          else  !------------------------------+                           |
            Call PakReal (wave,1,txt(1),n1)   !|                           |
          endif  !-----------------------------+                           |
                                                                          !|
          if (abs(energy).lt.1.E-33) Then !----+                           |
            txt(2) = 'not specified'          !|                           |
            n2 = Len_Trim(txt(2))             !|                           |
          else  !------------------------------+                           |
            Call PakReal (energy,1,txt(2),n2) !|                           |
          endif  !-----------------------------+                           |
                                                                          !|
          Call PakInt (iii,inmax,txt(3),n3)                               !|
                                                                          !|
          if (ipass.eq.1) Then  !---------------------------------------+  |
            if (modeout.eq.0) Then !----------------+                   |  |
              write(lunout,6,err=14)  code(1:n0)   !|                   |  |
              write(lundup,6,err=14)  code(1:n0)   !|                   |  |
            else !----------------------------------+                   |  |
              write(lunout,506,err=14)  code(1:n0) !|                   |  |
              write(lundup,506,err=14)  code(1:n0) !|                   |  |
            endif !---------------------------------+                   |  |
            Call flushoutput(lundup)                                   !|  |
  6         format(' <tr><td>Structure :   </td><td>',a,'</td></tr>')  !|  |
  506       format(' <tr><td colspan="2"><pre>'/
     *             ' Structure:',25x,'code=',a)                        !|  |
                                                                       !|  |
            if (ifail.eq.0 .or. ifail.gt.1) Then  !------------------+  |  |
c This is printed in any case except for material problem:          !|  |  |
              if (modeout.eq.0) Then !----------------------------+  |  |  |
                do i=lunout,lundup !============================+ |  |  |  |
                  write(i,361,err=14) tnb,tnb,                 !| |  |  |  |
     *                                Symmetry(isyngony+1),    !| |  |  |  |
     *                                tnb,tnb,                 !| |  |  |  |
     *                                rho                      !| |  |  |  |
                  if (isyngony.ne.8)                           !| |  |  |  |
     *              write(i,365,err=14) tnb,tnb,               !| |  |  |  |
     *                            (Lattice_constants(j),j=1,3),!| |  |  |  |
     *                                  tnb,tnb,               !| |  |  |  |
     *                            (Lattice_constants(j),j=4,6) !| |  |  |  |
                  if (isyngony.eq.1.AND.abs(Poisson).gt.1E-33) !| |  |  |  |
     *              write(i,367,err=14) tnb,tnb, Poisson       !| |  |  |  |
                  Call flushoutput(i)                          !| |  |  |  |
                  if (kcomp.gt.0) Then ! --------------------+ !| |  |  |  |
                    write(i,362,err=14) tnb                 !|  | |  |  |  |
                    if (mdetail.eq.1) write(i,462,err=14)   !|  | |  |  |  |
                    write(i,101,err=14) tnb                 !|  | |  |  |  |
                    do indx=1,kcomp !======================+ |  | |  |  |  |
                      if (mdetail.gt.0 .OR.               !| |  | |  |  |  |
     *                       indx.lt.kcomp) Then !-+       | |  | |  |  |  |
                        br = '<br>'               !|       | |  | |  |  |  |
                      else !-----------------------+       | |  | |  |  |  |
                        br = ' '                  !|       | |  | |  |  |  |
                      endif !----------------------+       | |  | |  |  |  |
                      write(i,363,err=14) name(indx),     !| |  | |  |  |  |
     *                                    kol(indx),      !| |  | |  |  |  |
     *                                    prcn(indx),br   !| |  | |  |  |  |
                      if (mdetail.gt.0) Then !----------+  | |  | |  |  |  |
                        do j=1,kol(indx) !===========+  |  | |  | |  |  |  |
                          if (j.lt.kol(indx) .OR.   !|  |  | |  | |  |  |  |
     *                     indx.lt.kcomp) Then !-+   |  |  | |  | |  |  |  |
                            br = '<br>'         !|   |  |  | |  | |  |  |  |
                          else !-----------------+   |  |  | |  | |  |  |  |
                            br = ' '            !|   |  |  | |  | |  |  |  |
                          endif !----------------+   |  |  | |  | |  |  |  |
                          write(i,364,err=14)       !|  |  | |  | |  |  |  |
     *                               wx(indx,j),    !|  |  | |  | |  |  |  |
     *                               wy(indx,j),    !|  |  | |  | |  |  |  |
     *                               wz(indx,j),br  !|  |  | |  | |  |  |  |
                        enddo !====================== = |  | |  | |  |  |  |
                      endif !---------------------------+  | |  | |  |  |  |
                    enddo !================================+ |  | |  |  |  |
                  endif !------------------------------------+  | |  |  |  |
                  write(i,201,err=14)                          !| |  |  |  |
                  Call flushoutput(i)                          !| |  |  |  |
                enddo !=========================================+ |  |  |  |
                                                                 !|  |  |  |
              else !----------------------------------------------+  |  |  |
                                                                 !|  |  |  |
                do i=lunout,lundup !============================+ |  |  |  |
                  write(i,761,err=14) Symmetry(isyngony+1),rho !| |  |  |  |
                  if (isyngony.ne.8)                           !| |  |  |  |
     *            write(i,765,err=14) Lattice_constants        !| |  |  |  |
                  if (isyngony.eq.1 .AND.                      !| |  |  |  |
     *                abs(Poisson).gt.1.E-33)                  !| |  |  |  |
     *                             write(i,767,err=14) Poisson !| |  |  |  |
                  Call flushoutput(i)                          !| |  |  |  |
                  if (kcomp.gt.0) Then ! ---------------------+ | |  |  |  |
                    write(i,762,err=14)   name(1),           !| | |  |  |  |
     *                                    kol(1),            !| | |  |  |  |
     *                                    prcn(1)            !| | |  |  |  |
                    if (mdetail.eq.1) write(i,764,err=14)    !| | |  |  |  |
     *                                   (wx(1,j),           !| | |  |  |  |
     *                                    wy(1,j),           !| | |  |  |  |
     *                                    wz(1,j),           !| | |  |  |  |
     *                                    j=1,kol(1))        !| | |  |  |  |
                    do indx=2,kcomp !========================+| | |  |  |  |
                      write(i,763,err=14) name(indx),       !|| | |  |  |  |
     *                                    kol(indx),        !|| | |  |  |  |
     *                                    prcn(indx)        !|| | |  |  |  |
                      if (mdetail.eq.1) write(i,764,err=14) !|| | |  |  |  |
     *                                   (wx(indx,j),       !|| | |  |  |  |
     *                                    wy(indx,j),       !|| | |  |  |  |
     *                                    wz(indx,j),       !|| | |  |  |  |
     *                                    j=1,kol(indx))    !|| | |  |  |  |
                    enddo !==================================+| | |  |  |  |
                  endif !-------------------------------------+ | |  |  |  |
                  write(i,701,err=14)                          !| |  |  |  |
                  Call flushoutput(i)                          !| |  |  |  |
                enddo !=========================================+ |  |  |  |
              endif !---------------------------------------------+  |  |  |
            endif  !-------------------------------------------------+  |  |
                                                                       !|  |
  361   format(' <tr><td',a,'>Symmetry : </td>',
     *              '<td',a,'>',a,'</td></tr>'/
     *         ' <tr><td',a,'>Density (gm/cm<sup>3</sup>) : </td>',
     *              '<td',a,'>',g12.5,'<sup>&nbsp;</sup></td></tr>')
  365   format(' <tr><td',a,'>Unit cell constants (A) : </td>',
     *              '<td',a,'>',2(g12.5,' , '),g12.5,'</td></tr>'/
     *         ' <tr><td',a,'>Unit cell angles (degr.) : </td>',
     *              '<td',a,'>',2(g12.5,' , '),g12.5,'</td></tr>')
  367   format(' <tr><td',a,'>Poisson Ratio : </td>',
     *              '<td',a,'>',f8.4,'</td></tr>')
  362   format(' <tr><td',a,'>Composition : Atom -- N_sites',1x
     *                                       ,'(site occupation)',$)
  462   format(' <br><span style="padding-left:80px;">... ',
     *                                       'Atom x,y,z :</span>')
  101   format(' </td><td',a,'>',$)
  363   format(a,' --',i3,' (',f7.3,')',a)
  364   format(f8.5,',',f8.5,',',f8.5,a)
  201   format(' </td></tr>')

  761   format(' Symmetry:                       syngony=',a/
     *         ' Density (g/cm^3):                   rho=',g12.5)
  765   format(' Unit cell constants(A):              a1=',g12.5/
     *         '                                      a2=',g12.5/
     *         '                                      a3=',g12.5/
     *         ' Unit cell angles(degr):              a4=',g12.5/
     *         '                                      a5=',g12.5/
     *         '                                      a6=',g12.5)
  767   format(' Poisson Ratio:                  poisson=',f8.4)
  762   format(' Composition:',
     *     23x,'atom=',a,'  nsites=',i3.3,'  occupation=',f5.3)
  763   format(36x,'atom=',a,'  nsites=',i3.3,'  occupation=',f5.3)
  764   format(99(39x,'x=',f7.5,'  y=',f7.5,'  z=',f7.5:/))
  701   format(' </pre>')
                                                                       !|  |
            if (Len_Trim(radiat).eq.0) radiat = 'none'                 !|  |
            l = Len_Trim(name_nearest)                                 !|  |
            Call PakReal(edge_nearest,1,txt(9),n9)                     !|  |
            if (modeout.eq.0) Then !--------------------------------+   |  |
              do i=lunout,lundup !================================+ |   |  |
                if (edge_nearest.gt.0. .AND. l.gt.0) Then !-----+ | |   |  |
                  write(i,7,err=14) ttb,ttb,radiat,            !| | |   |  |
     *                              tnb,tnb,txt(1)(1:n1),      !| | |   |  |
     *                              tnb,tnb,txt(2)(1:n2),      !| | |   |  |
     *                              tnb,tnb,txt(9)(1:n9),      !| | |   |  |
     *                                      name_nearest(1:l)  !| | |   |  |
                else  !-----------------------------------------+ | |   |  |
                  write(i,7,err=14) ttb,ttb,radiat,            !| | |   |  |
     *                              tnb,tnb,txt(1)(1:n1),      !| | |   |  |
     *                              tnb,tnb,txt(2)(1:n2)       !| | |   |  |
                endif  !----------------------------------------+ | |   |  |
                Call flushoutput(i)                              !| |   |  |
              enddo !=============================================+ |   |  |
            else !--------------------------------------------------+   |  |
              do i=lunout,lundup !============================+     |   |  |
                if (edge_nearest.gt.0. .AND. l.gt.0) Then !-+ |     |   |  |
                  write(i,571,err=14) radiat,              !| |     |   |  |
     *                                txt(1)(1:n1),        !| |     |   |  |
     *                                txt(2)(1:n2),        !| |     |   |  |
     *                                txt(9)(1:n9),        !| |     |   |  |
     *                                name_nearest(1:l)    !| |     |   |  |
                else  !-------------------------------------+ |     |   |  |
                  write(i,571,err=14) radiat,              !| |     |   |  |
     *                                txt(1)(1:n1),        !| |     |   |  |
     *                                txt(2)(1:n2)         !| |     |   |  |
                endif  !------------------------------------+ |     |   |  |
                write (i,701,err=14)                         !|     |   |  |
                Call flushoutput(i)                          !|     |   |  |
              enddo !=========================================+     |   |  |
            endif !-------------------------------------------------+   |  |
  7         format(
     *      ' <tr><td',a,'>X-ray line : </td>',
     *           '<td',a,'>',a,'</td></tr>'/
c
     *      ' <tr><td',a,'>Wavelength (A) : </td>',
     *           '<td',a,'>',a,'</td></tr>'/
c
     *      ' <tr><td',a,'>Energy (keV) : </td>',
     *           '<td',a,'>',a,'</td></tr>':/
c
     *      ' <tr><td',a,'>Closest absorption edge (keV) : </td>',
     *           '<td',a,'>',a,' (for atom <b>',a,'</b>)</td></tr>')
  571       format(' <pre>'/
     *      ' X-ray line:                        line=',a/
     *      ' Wavelength (A):                    wave=',a/
     *      ' Energy (keV):                    energy=',a,:/
     *      ' Closest absorption edge:'/
     *      '      -- energy(keV):             edgeEn=',a/
     *      '      -- element:                 edgeEl=',a)             !|  |
c 701       format(' </pre>')                                          !|  |
                                                                       !|  |
          endif  !------------------------------------------------------+  |
                                                                          !|
          if (ifail.ne.10) Then  !--------------------------------+        |
c This is printed if the database option is defined:              |        |
c may be no Bragg reflection:                                     |        |
            ldb = Len_Trim(DBtext(iHenkeCowan+2))                !|        |
            do i=lunout,lundup !==============================+   |        |
              if (modeout.eq.0) Then !---------------------+  |   |        |
                write (i,220,err=14)                      !|  |   |        |
     *                       DBtext(iHenkeCowan+2)(1:ldb) !|  |   |        |
              else  !--------------------------------------+  |   |        |
                write (i,720,err=14)                      !|  |   |        |
     *                       DBtext(iHenkeCowan+2)(1:ldb) !|  |   |        |
              endif  !----------------------------------------+   |        |
            Call flushoutput(i)                              !|   |        |
            enddo  !==========================================+   |        |
          endif  !------------------------------------------------+        |
  220     format(
     *    ' <tr><td>Database for df<sub>1</sub>, df<sub>2</sub> :</td>',
     *         '<td><font color="red"> *** ',a,' *** </font></td></tr>')  !|
  720     format(' <pre>'/' DB for df1,df2:',18x,'dbtype=',a/' </pre>')   !|
                                                                          !|
          if (ifail.eq.0 .or. ifail.eq.8)  Then  !------------------+      |
c This is printed if there is material and x-ray data;              |      |
c may be no Bragg reflection:                                       |      |
            do i=lunout,lundup !===============================+    |      |
              if (modeout.eq.0) Then !--------------------+    |    |      |
                write (i,72,err=14)                      !|    |    |      |
     *                tnb, tnb,  xr0,     xi0,           !|    |    |      |
     *                tnb, tnb, -xr0/2., -xi0/2.,        !|    |    |      |
     *                tnb, quest0(1:lq0), tnb, mu0, a0,  !|    |    |      |
     *                tnb, quest0(1:lq0), tnb, dl_TER,   !|    |    |      |
     *                tnb, quest0(1:lq0), tnb, TER/gra,  !|    |    |      |
     *                                         TER*1000. !|    |    |      |
              else  !-------------------------------------+    |    |      |
                write(i,572,err=14)                      !|    |    |      |
     *                           xr0,     xi0,           !|    |    |      |
     *                          -xr0/2., -xi0/2.,        !|    |    |      |
     *                           mu0, a0,                !|    |    |      |
     *                           dl_TER,                 !|    |    |      |
     *                           TER/gra, TER*1000.,     !|    |    |      |
     *                           quest0(1:lq0)           !|    |    |      |
              endif  !------------------------------------+    |    |      |
c Limitations of ter_sl:                                       |    |      |
              if (TER  .gt. 0.   .AND.                        !|    |      |
     *            wave .ge. 0.01 .AND.                        !|    |      |
     *            wave .lt. 1000.) Then !---------------+      |    |      |
c Round the interval to the upper side and 2 digits:    |      |    |      |
                rangmin = 0.                           !|      |    |      |
                rangmax = round (3.*TER/gra, 2, +1)    !|      |    |      |
                iunis   = 0          !in degrees       !|      |    |      |
                npts    = 401                          !|      |    |      |
c icoway (x0h form):                                    |      |    |      |
c                 0=crystal, 1=amorphous, 2=chem_formula|      |    |      |
c isubway (TER form):                                   |      |    |      |
c                 1=database_code, 2=chem_formula, 3=x0 |      |    |      |
                if (icoway.lt.2) Then  !---+            |      |    |      |
                  isubway = 1             !|database    |      |    |      |
                else !---------------------+            |      |    |      |
                  isubway = 2             !|formula     |      |    |      |
                endif  !-------------------+            |      |    |      |
c Write HTML form input to launch ter_sl:               |      |    |      |
                write (i,73,err=14)                    !|      |    |      |
     *                       spacer(1:lsp),            !|      |    |      |
     *                       txcount,                  !|      |    |      |
     *                       wave, radiat,             !|      |    |      |
     *                       isubway,                  !|      |    |      |
     *                       code(1:n0),               !|      |    |      |
     *                       code(1:n0), rho,          !|      |    |      |
     *                       iHenkeCowan,              !|      |    |      |
     *                       rangmin, rangmax,         !|      |    |      |
     *                       iunis, npts               !|      |    |      |
c End of HTML form for launching ter_sl                 |      |    |      |
              elseif (modeout.ne.0) Then !--------------+      |    |      |
                write (i,701,err=14)                   !|      |    |      |
              endif  !----------------------------------+      |    |      |
              if (modeout.eq.0) write (i,75,err=14)           !|    |      |
            enddo  !===========================================+    |      |
          endif  !--------------------------------------------------+      |
          Call flushoutput(lundup)                                        !|
  73      format(1x,a/
     *    ' <form action="/cgi/ter_form.pl" method="get">'/
     *    ' <input type="hidden" name="comment1"',1x,
     *                            'value="template: x0h_form, ',a,'">'/
     *    ' <input type="hidden" name="xway"    value="1">'/
     *    ' <input type="hidden" name="wave"    value="',g15.8,'">'/
     *    ' <input type="hidden" name="line"    value="',a,'">'/
     *    ' <input type="hidden" name="ipol"    value="1">'/
     *    ' <input type="hidden" name="subway"  value="',i1,'">'/
     *    ' <input type="hidden" name="code"    value="',a,'">'/
     *    ' <input type="hidden" name="chem"    value="',a,'">'/
     *    ' <input type="hidden" name="rho"     value="',g15.8,'">'/
     *    ' <input type="hidden" name="df1df2"  value="',i2,'">'/
     *    ' <input type="hidden" name="x0"      value="(0., 0.)">'/
     *    ' <input type="hidden" name="w0"      value="1.">'/
     *    ' <input type="hidden" name="sigma"   value="0.">'/
     *    ' <input type="hidden" name="tr"      value="0.">'/
     *    ' <input type="hidden" name="scanmin" value="',g12.5,'">'/
     *    ' <input type="hidden" name="scanmax" value="',g12.5,'">'/
     *    ' <input type="hidden" name="unis"    value="',i1,'">'/
     *    ' <input type="hidden" name="nscan"   value="',i4,'">'/
     *    ' <input type="hidden" name="swflag"  value="0">'/
     *    ' <input type="image" src="images/get_the_curve.gif"',
     *                    1x,'border="0" width="102" height="12"',
     *                    1x,'alt="Get the reflectivity curve">'/
     *    ' </form>')
  72      format(
     *    ' <tr><td',a,'><i>x<sub>r0</sub></i> ,',1x,
     *         '<i>x<sub>i0</sub></i> &nbsp; &nbsp; &nbsp; ',
     *         '(<i>n = 1 +   x<sub>r0&nbsp;</sub>/2',1x,
     *                   '+ i*x<sub>i0&nbsp;</sub>/2</i>) : ',
     *                                         '&nbsp; </td>',
     *         '<td',a,'>',g12.5,' , &nbsp; &nbsp; ',g12.5,'</td></tr>'/
c
     *    ' <tr><td',a,'><i>delta</i> , <i>eta</i> &nbsp; &nbsp;',
     *         ' (<i>n = 1 - delta - i*eta</i>) : </td>',
     *         '<td',a,'>',g12.5,' , &nbsp; &nbsp; ',g12.5,'</td></tr>'/
c
     *    ' <tr><td',a,'>',a,'Absorption factor (1/cm) and length ',
     *                                              '(um) : </td>',
     *         '<td',a,'>',g12.5,' , &nbsp; &nbsp; ',g12.5,'</td></tr>'/
c
     *    ' <tr><td',a,'>',a,'Extinction length at TER (A) : </td>'
     *         '<td',a,'>',g12.5,'</td></tr>'/
c
     *    ' <tr><td',a,'>',a,'Critical angle for TER ',
     *                                     '(degr., mrad) : </td>',
     *         '<td',a,'>',g12.5,' , &nbsp; &nbsp; ',g12.5)
c
  75      format(' </td></tr>')
c <tt> forces monospace font; also we set "pre {font-family:Consolas,monospace}" in style.css
  572     format(' <pre>'/
     *    ' <tt><i>n = 1 + xr0/2 + i*xi0/2</i></tt>:',19x,'xr0=',g12.5/
     *    '                                     xi0=',g12.5/
     *    ' <tt><i>n = 1 - delta - i*eta</i></tt>:',19x,'delta=',g12.5/
     *    '                                     eta=',g12.5/
     *    ' Absorption factor(1/cm):            mu0=',g12.5/
     *    ' Absorption length(um):               a0=',g12.5/
     *    ' TER extinction length(A):           ex0=',g12.5/
     *    ' TER critical angle(degr.):      terDegr=',g12.5/
     *    ' TER critical angle(mrad):       terMrad=',g12.5/
     *    ' TER parameters help: ',a)
c 701       format(' </pre>')                                             !|
                                                                          !|
          if (iref.gt.0)      Then  !-----------------------------------+  |
            if (modeout.eq.0) Then !------------------+                 |  |
              write(lunout,8,err=14)  txt(3)(1:n3)   !|                 |  |
              write(lundup,8,err=14)  txt(3)(1:n3)   !|                 |  |
            else !------------------------------------+                 |  |
              write(lunout,508,err=14)  txt(3)(1:n3) !|                 |  |
              write(lundup,508,err=14)  txt(3)(1:n3) !|                 |  |
            endif !-----------------------------------+                 |  |
  8         format(' <tr><td>Reflection : </td><td>(',a,')</td></tr>') !|  |
  508       format(' <pre>'/' Reflection:',25x,'hkl=',a/' </pre>')     !|  |
                                                                       !|  |
            if (x0_abs.gt.0.) Then !-------------------+                |  |
              xh_weak_limit = Amax1( (1.E-8)*sin1,    !|order of gh     |  |
     *                               (1.E-3)*x0_abs ) !|0.1% reflec.    |  |
            else  !------------------------------------+                |  |
              xh_weak_limit = 0.0                     !|                |  |
            endif  !-----------------------------------+                |  |
                                                                       !|  |
            if (ifail.eq.0)   Then  !---------------------------------+ |  |
              if (ipass.eq.1) Then  !------------------------+        | |  |
                if (modeout.eq.0) Then !------------------+  |        | |  |
                  do i=lunout,lundup !=================+  |  |        | |  |
                    write(i,81,err=14)                !|  |  |        | |  |
     *                      tnb, tnb, qb8,            !|  |  |        | |  |
     *                      tnb, tnb, d8,             !|  |  |        | |  |
     *                      tnb, tnb, sin1, cos1,     !|  |  |        | |  |
     *                      tnb, tnb, tg1, ctg1,      !|  |  |        | |  |
     *                      tnb, tnb, sin2, cos2,     !|  |  |        | |  |
     *                      tnb, tnb, ESinTheta       !|  |  |        | |  |
c ccc*                     ,tnb, tnb, xrh, xih, xdf   !|  |  |        | |  |
                  enddo !==============================+  |  |        | |  |
                else  !-----------------------------------+  |        | |  |
                  do i=lunout,lundup !=================+  |  |        | |  |
                    write(i,581,err=14) qb8, d8,      !|  |  |        | |  |
     *                                  sin1, cos1,   !|  |  |        | |  |
     *                                  tg1, ctg1,    !|  |  |        | |  |
     *                                  sin2, cos2,   !|  |  |        | |  |
     *                                  ESinTheta     !|  |  |        | |  |
c ccc*                                 ,xrh, xih, xdf !|  |  |        | |  |
                  enddo !==============================+  |  |        | |  |
                endif  !----------------------------------+  |        | |  |
              endif  !---------------------------------------+        | |  |
  81    format(
     *  ' <tr><td',a,'>Bragg angle (degr.) : </td>',
     *       '<td',a,'> ',g12.5,                  '</td></tr>'/
c
     *  ' <tr><td',a,'>Interplanar spacing (A) : </td>',
     *       '<td',a,'>',g12.5,                   '</td></tr>'/
c
     *  ' <tr><td',a,'>sin(QB) , cos(QB) : </td>',
     *       '<td',a,'>',g12.5,' , &nbsp; ',g12.5,'</td></tr>'/
c
     *  ' <tr><td',a,'>tan(QB) , cotan(QB) : </td>',
     *       '<td',a,'>',g12.5,' , &nbsp; ',g12.5,'</td></tr>'/
c
     *  ' <tr><td',a,'>sin(2*QB) , cos(2*QB) : </td>',
     *       '<td',a,'>',g12.5,' , &nbsp; ',g12.5,'</td></tr>'/
c
     *  ' <tr><td',a,'>ESinTheta=12.398/(2d) : </td>',
     *       '<td',a,'>',g14.7,                   '</td></tr>')

  581   format(' <pre>'/
     *  ' Bragg angle (degr.):                 QB=',g12.5/
     *  ' Interplanar spacing (A):              d=',g12.5/
     *  '                                 sin(QB)=',g12.5/
     *  '                                 cos(QB)=',g12.5/
     *  '                                 tan(QB)=',g12.5/
     *  '                               cotan(QB)=',g12.5/
     *  '                               sin(2*QB)=',g12.5/
     *  '                               cos(2*QB)=',g12.5/
     *  '                               ESinTheta=',g14.7/
     *  ' </pre>')                                                   !| |  |
                                                                     !| |  |
c Show the x0h data only if |xh|/|x0| > 0.001% = 1.e-4               !| |  |
              if (prcsg .gt. 0.01)  Then  !------------------------+  | |  |
                                                                  !|  | |  |
                iunis = 3               !in seconds               !|  | |  |
                npts  = 401                                       !|  | |  |
                if (sin2 .gt. 1E-5) Then !----------+              |  | |  |
                  shift = Abs(xr0) / (sec*sin2)    !|              |  | |  |
                else !------------------------------+              |  | |  |
                  shift = 0.            !QB=90     !|              |  | |  |
                endif !-----------------------------+              |  | |  |
c Round to 2 digits towards a larger number:                      !|  | |  |
                rangmax = round ( 3.*halfsg+shift, 2, +1)         !|  | |  |
c Round to 2 digits towards a smaller number:                     !|  | |  |
                rangmin = round (-3.*halfsg+shift, 2, -1)         !|  | |  |
c               rangmin =-rangmax                                 !|  | |  |
                                                                  !|  | |  |
c SIGMA POLARIZATION:                                              |  | |  |
                ipol = 1                                          !|  | |  |
                pl   = 'sg'                                       !|  | |  |
                xh   = Cmplx(xrsg+xisg*sin(xdfsg*pi),             !|  | |  |
     *                            xisg*cos(xdfsg*pi))             !|  | |  |
                do i=lunout,lundup !============================+  |  | |  |
                  if (modeout.eq.0) Then !-------------------+  |  |  | |  |
                    write(i,184,err=14)                     !|  |  |  | |  |
     *                ttb, ttb, 'Sigma',                    !|  |  |  | |  |
     *                tnb, tnb, xrsg, xisg,                 !|  |  |  | |  |
     *                tnb, tnb, xdfsg,                      !|  |  |  | |  |
     *                tnb, tnb, prcsg                       !|  |  |  | |  |
                    if (ah .gt. 0.)                         !|  |  |  | |  |
     *                write(i,185,err=14)                   !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, ah,        !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, ah1sg,     !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, ah2sg,     !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, extisg_l   !|  |  |  | |  |
                    write(i,186,err=14)                     !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, extisg,    !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb,            !|  |  |  | |  |
     *                                      sqr2*halfsg,    !|  |  |  | |  |
     *                                      sqr2*halfsg*s6, !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb,            !|  |  |  | |  |
     *                                           halfsg,    !|  |  |  | |  |
     *                                           halfsg*s6  !|  |  |  | |  |
                  else !-------------------------------------+  |  |  | |  |
                    write(i,684,err=14)                     !|  |  |  | |  |
     *                       'Sigma',                       !|  |  |  | |  |
     *                       pl, xrsg, pl, xisg,            !|  |  |  | |  |
     *                       pl, xdfsg, pl, prcsg           !|  |  |  | |  |
                    if (ah .gt. 0.)                         !|  |  |  | |  |
     *              write(i,685,err=14)                     !|  |  |  | |  |
     *                       ah, pl, ah1sg, pl, ah2sg,      !|  |  |  | |  |
     *                       pl, extisg_l                   !|  |  |  | |  |
                    write(i,686,err=14)                     !|  |  |  | |  |
     *                       pl, extisg,                    !|  |  |  | |  |
     *                       pl, sqr2*halfsg,               !|  |  |  | |  |
     *                       pl, sqr2*halfsg*s6,            !|  |  |  | |  |
     *                       pl, halfsg,                    !|  |  |  | |  |
     *                       pl, halfsg*s6,              !|  |  |  | |  |
     *                       questh(1:lqh)                  !|  |  |  | |  |
                  endif !------------------------------------+  |  |  | |  |
c Limitations of gid_sl:                                        |  |  | |  |
                  if (Abs(xh) .gt. xh_weak_limit .AND.         !|  |  | |  |
     *              wave.ge.0.1 .AND. wave.lt.10.) Then !+      |  |  | |  |
c Write HTML form input to launch gid_sl:                |      |  |  | |  |
                      write (i,84,err=14)               !|      |  |  | |  |
     *                   spacer(1:lsp),                 !|      |  |  | |  |
     *                   txcount,                       !|      |  |  | |  |
     *                   wave, radiat, ipol,            !|      |  |  | |  |
     *                   code(1:n0),                    !|      |  |  | |  |
     *                   iHenkeCowan,                   !|      |  |  | |  |
     *                   iii(1),iii(2),iii(inmax),      !|      |  |  | |  |
     *                   iii(1),iii(2),iii(inmax),      !|      |  |  | |  |
     *                   rangmin, rangmax,              !|      |  |  | |  |
     *                   iunis, npts,                   !|      |  |  | |  |
     *                   'sigma'                        !|      |  |  | |  |
c End of HTML form for launching gid_s                   |      |  |  | |  |
                  elseif (modeout.ne.0) Then !-----------+      |  |  | |  |
                    write(i,701,err=14)                 !|      |  |  | |  |
                  endif  !-------------------------------+      |  |  | |  |
                  if (modeout.eq.0) write(i,75,err=14)         !|  |  | |  |
                enddo  !========================================+  |  | |  |
                                                                  !|  | |  |
c PI POLARIZATION:                                                 |  | |  |
                ipol = 2                                          !|  | |  |
                pl   = 'pi'                                       !|  | |  |
                xh   = Cmplx(xrpi+xipi*sin(xdfpi*pi),             !|  | |  |
     *                            xipi*cos(xdfsg*pi))             !|  | |  |
                do i=lunout,lundup !============================+  |  | |  |
c Write HTML form input to launch gid_sl:                       |  |  | |  |
c Limitations of gid_sl:                                        |  |  | |  |
                  if (modeout.eq.0) Then !-------------------+  |  |  | |  |
                    write(i,184,err=14)                     !|  |  |  | |  |
     *                ttb, ttb, 'Pi',                       !|  |  |  | |  |
     *                tnb, tnb, xrpi, xipi,                 !|  |  |  | |  |
     *                tnb, tnb, xdfpi,                      !|  |  |  | |  |
     *                tnb, tnb, prcpi                       !|  |  |  | |  |
                    if (ah .gt. 0.)                         !|  |  |  | |  |
     *                write(i,185,err=14)                   !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, ah,        !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, ah1pi,     !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, ah2pi,     !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, extipi_l   !|  |  |  | |  |
                    write(i,186,err=14)                     !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, extipi,    !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb,            !|  |  |  | |  |
     *                                      sqr2*halfpi,    !|  |  |  | |  |
     *                                      sqr2*halfpi*s6, !|  |  |  | |  |
     *                  tnb, questh(1:lqh), tnb, halfpi,    !|  |  |  | |  |
     *                                           halfpi*s6  !|  |  |  | |  |
                  else !-------------------------------------+  |  |  | |  |
                    write(i,684,err=14)                     !|  |  |  | |  |
     *                       'Pi',                          !|  |  |  | |  |
     *                       pl, xrpi, pl, xipi,            !|  |  |  | |  |
     *                       pl, xdfpi, pl, prcpi           !|  |  |  | |  |
                    if (ah .gt. 0.)                         !|  |  |  | |  |
     *              write(i,685,err=14)                     !|  |  |  | |  |
     *                       ah, pl, ah1pi, pl, ah2pi,      !|  |  |  | |  |
     *                       pl, extipi_l                   !|  |  |  | |  |
                    write(i,686,err=14)                     !|  |  |  | |  |
     *                       pl, extipi,                    !|  |  |  | |  |
     *                       pl, sqr2*halfpi,               !|  |  |  | |  |
     *                       pl, sqr2*halfpi*s6,            !|  |  |  | |  |
     *                       pl, halfpi,                    !|  |  |  | |  |
     *                       pl, halfpi*s6,                 !|  |  |  | |  |
     *                       questh(1:lqh)                  !|  |  |  | |  |
                  endif !------------------------------------+  |  |  | |  |
c Limitations of gid_sl:                                        |  |  | |  |
                  if (Abs(xh) .gt. xh_weak_limit .AND.         !|  |  | |  |
     *              wave.ge.0.1 .AND. wave.lt.10.) Then !+      |  |  | |  |
c Write HTML form input to launch gid_sl:                |      |  |  | |  |
                      write (i,84,err=14)               !|      |  |  | |  |
     *                   spacer(1:lsp),                 !|      |  |  | |  |
     *                   txcount,                       !|      |  |  | |  |
     *                   wave, radiat, ipol,            !|      |  |  | |  |
     *                   code(1:n0),                    !|      |  |  | |  |
     *                   iHenkeCowan,                   !|      |  |  | |  |
     *                   iii(1),iii(2),iii(inmax),      !|      |  |  | |  |
     *                   iii(1),iii(2),iii(inmax),      !|      |  |  | |  |
     *                   rangmin, rangmax,              !|      |  |  | |  |
     *                   iunis, npts,                   !|      |  |  | |  |
     *                   'pi'                           !|      |  |  | |  |
c End of HTML form for launching gid_s                   |      |  |  | |  |
                  elseif (modeout.ne.0) Then !-----------+      |  |  | |  |
                    write(i,701,err=14)                 !|      |  |  | |  |
                  endif  !-------------------------------+      |  |  | |  |
                  if (modeout.eq.0) write(i,75,err=14)         !|  |  | |  |
                enddo  !========================================+  |  | |  |
                                                                  !|  | |  |
              else  !----------------------------------------------+  | |  |
                                                                  !|  | |  |
                if (modeout.eq.0) Then !--------+                  |  | |  |
                  write(lunout,385,err=14)     !|                  |  | |  |
                  write(lundup,385,err=14)     !|                  |  | |  |
                else !--------------------------+                  |  | |  |
                  write(lunout,585,err=14)     !|                  |  | |  |
                  write(lundup,585,err=14)     !|                  |  | |  |
                endif !-------------------------+                  |  | |  |
  385           format(' <tr><td><b><i>Attention: </i></b></td>'/
     *          ' <td><b>Forbidden reflection </b></td></tr>')    !|  | |  |
  585           format(' <pre>'/' <b><i>Attention: </i>',1x,
     *          ' Forbidden reflection </b>'/' </pre>')           !|  | |  |
                                                                  !|  | |  |
              endif !----------------------------------------------+  | |  |
            endif  !--------------------------------------------------+ |  |
          endif  !------------------------------------------------------+  |
                                                                          !|
          if (linesX0hWarn_ .gt. 0 .OR.                                   !|
     *        linesX0hWarn  .gt. 0)  Then  !--------------------+          |
            if (modeout.eq.0) Then !--------------------------+ |          |
              do  i=lunout,lundup !=========================+ | |          |
                write(i,91,err=14)                         !| | |          |
                if (linesX0hWarn_.gt.linesX0hWarn) Then !-+ | | |          |
c Print warnings from the 1st X0h call:                   | | | |          |
                  do j=1,linesX0hWarn_  !===============+ | | | |          |
                    n0=Max(Len_Trim(txtX0hWarn_(j)),1) !| | | | |          |
                    write(i,92,err=14)                 !| | | | |          |
     *                            txtX0hWarn_(j)(1:n0) !| | | | |          |
                  enddo !===============================+ | | | |          |
                endif !-----------------------------------+ | | |          |
c Print warnings from the 2nd X0h call:                     | | |          |
                do j=1,linesX0hWarn  !================+     | | |          |
                  n0=Max(Len_Trim(txtX0hWarn(j)),1)  !|     | | |          |
                  write(i,92,err=14)                 !|     | | |          |
     *                          txtX0hWarn(j)(1:n0)  !|     | | |          |
                enddo !===============================+     | | |          |
                write(i,75,err=14)                         !| | |          |
              enddo !=======================================+ | |          |
            else   !------------------------------------------+ |          |
              do  i=lunout,lundup !=========================+ | |          |
                write(i,93,err=14)                         !| | |          |
                if (linesX0hWarn_.gt.linesX0hWarn) Then !-+ | | |          |
c Print warnings from the 1st X0h call:                   | | | |          |
                  do j=1,linesX0hWarn_  !===============+ | | | |          |
                    n0=Max(Len_Trim(txtX0hWarn_(j)),1) !| | | | |          |
                    write(i,92,err=14)                 !| | | | |          |
     *                            txtX0hWarn_(j)(1:n0) !| | | | |          |
                  enddo !===============================+ | | | |          |
                endif !-----------------------------------+ | | |          |
c Print warnings from the 2nd X0h call:                     | | |          |
                do j=1,linesX0hWarn  !================+     | | |          |
                  n0=Max(Len_Trim(txtX0hWarn(j)),1)  !|     | | |          |
                  write(i,92,err=14)                 !|     | | |          |
     *                          txtX0hWarn(j)(1:n0)  !|     | | |          |
                enddo !===============================+     | | |          |
              enddo !=======================================+ | |          |
            endif  !------------------------------------------+ |          |
          endif !-----------------------------------------------+          |
  91      format(' <tr><td colspan="2"><b><i>X0h Warnings: </i></b>')     !|
  92      format(' <br>',a)                                               !|
  93      format(' <b><i>X0h Warnings: </i></b>')                         !|
                                                                          !|
          iHenkeCowan = iHenkeCowan + 2                                   !|
        enddo  !===========================================================+
  84    format(1x,a/
     *  ' <form action="/cgi/gid_form.pl" method="get">'/
     *  ' <input type="hidden" name="prg" value="gid_slm7">'/
     *  ' <input type="hidden" name="reduction" value="2">'/
     *  ' <input type="hidden" name="comment1"',1x,
     *                            'value="template: x0h_form, ',a,'">'/
     *  ' <input type="hidden" name="xway"    value="1">'/
     *  ' <input type="hidden" name="wave"    value="',g15.8,'">'/
     *  ' <input type="hidden" name="line"    value="',a,'">'/
     *  ' <input type="hidden" name="ipol"    value="',i1,'">'/
     *  ' <input type="hidden" name="code"    value="',a,'">'/
     *  ' <input type="hidden" name="df1df2"  value="',i2,'">'/
     *  ' <input type="hidden" name="w0"      value="1.">'/
     *  ' <input type="hidden" name="wh"      value="1.">'/
     *  ' <input type="hidden" name="sigma"   value="0.">'/
     *  ' <input type="hidden" name="daa"     value="0.">'/
     *  ' <input type="hidden" name="i1"      value="',i4,'">'/
     *  ' <input type="hidden" name="i2"      value="',i4,'">'/
     *  ' <input type="hidden" name="i3"      value="',i4,'">'/
     *  ' <input type="hidden" name="n1"      value="',i4,'">'/
     *  ' <input type="hidden" name="n2"      value="',i4,'">'/
     *  ' <input type="hidden" name="n3"      value="',i4,'">'/
     *  ' <input type="hidden" name="m1"      value="0">'/
     *  ' <input type="hidden" name="m2"      value="0">'/
     *  ' <input type="hidden" name="m3"      value="0">'/
     *  ' <input type="hidden" name="miscut"  value="0.">'/
     *  ' <input type="hidden" name="unim"    value="0">'/
     *  ' <input type="hidden" name="igie"    value="5">'/
     *  ' <input type="hidden" name="fcentre" value="0.">'/
     *  ' <input type="hidden" name="unic"    value="0">'/
     *  ' <input type="hidden" name="axis"    value="4">'/
     *  ' <input type="hidden" name="a1"      value="0">'/
     *  ' <input type="hidden" name="a2"      value="0">'/
     *  ' <input type="hidden" name="a3"      value="0">'/
     *  ' <input type="hidden" name="invert"  value="0">'/
     *  ' <input type="hidden" name="scanmin" value="',g12.5,'">'/
     *  ' <input type="hidden" name="scanmax" value="',g12.5,'">'/
     *  ' <input type="hidden" name="unis"    value="',i1,'">'/
     *  ' <input type="hidden" name="nscan"   value="',i4,'">'/
     *  ' <input type="hidden" name="column"  value="A">'/
     *  ' <input type="image" src="images/get_the_curve.gif"',
     *       ' border="0" width="102" height="12"',
     *       ' alt="Get the Bragg curve (',a,')"></form>')
  184   format(
     *  ' <tr><td',a,'>Polarization : <sub>&nbsp;</sub></td>',
     *       '<td',a,'><b><i>',a,'</i></b></td></tr>'/
c
     *  ' <tr><td',a,'>|<i>x<sub>rh</sub></i>| , ',
     *                '|<i>x<sub>ih</sub></i>| : &nbsp;</td>',
     *       '<td',a,'>',g12.5,' , &nbsp; &nbsp; ',g12.5,' </td></tr>'/
c
     *  ' <tr><td',a,'>Phase difference',1x,
     *               '(<i>x<sub>rh</sub></i> - ',
     *                '<i>x<sub>ih</sub></i>) : &nbsp;</td>',
     *       '<td',a,'>',g12.5,' * pi</td></tr>'/
c
     *  ' <tr><td',a,'>Relative intensity ',
     *               '(<i>x<sub>h</sub></i> / ',
     *                '<i>x<sub>0</sub></i>) : &nbsp;</td>',
     *       '<td',a,'>',g12.5,' %</td></tr>')
c
  185   format(
     *  ' <tr><td',a,'>',a,'Laue-case absorption depth',
     *                              ' (off Bragg, um) : </td>',
     *       '<td',a,'>',g12.5,'</td></tr>'/
c
     *  ' <tr><td',a,'>',a,'Laue-case absorption depth',
     *                             ' (Borrmann field) : </td>',
     *       '<td',a,'>',g12.5,'</td>'/
c
     *  ' <tr><td',a,'>',a,'Laue-case absorption depth',
     *                              ' (anti-Borrmann) : </td>',
     *       '<td',a,'>',g12.5,'</td></tr>'/
c
     *  ' <tr><td',a,'>',a,'Symmetric Laue-case extinction',
     *                                  ' length (um) : </td>',
     *       '<td',a,'>',g12.5,'</td></tr>')

  186   format(
     *  ' <tr><td',a,'>',a,'Symmetric Bragg-case extinction',
     *                                  ' length (um) : </td>',
     *       '<td',a,'>',g12.5,'</td></tr>'/
c
     *  ' <tr><td',a,'>',a,'<b>.....</b> Double-crystal curve',
     *                  ' FWHM (arcsec., urad) : &nbsp; </td>',
     *       '<td',a,'>',g12.5,' , &nbsp; &nbsp; ',g12.5,'</td></tr>'/
c
     *  ' <tr><td',a,'>',a,'<b>.....</b> Darwin curve FWHM',
     *                       ' (arcsec., urad) : &nbsp; </td>',
     *       '<td',a,'>',g12.5,' , &nbsp; &nbsp; ',g12.5)

c 75    format(' </td></tr>')

  684   format(' <pre>'/
     *  ' Polarization:                       pol=',a/
     *  ' |xrh|:                            xrh',a,'=',g12.5/
     *  ' |xih|:                            xih',a,'=',g12.5/
     *  ' Phase difference xrh-xih:         xdf',a,'=',g12.5/
     *  ' Relative intensity xh/x0:        prcn',a,'=',g12.5)
  685   format(
     *  ' Symmetric Laue case parameters:'/
     *  '  Absorption depth, off Bragg (um):   ah=',g12.5/
     *  '  Absorp.depth, Borrman field: ah_brmn',a,'=',g12.5/
     *  '  Absorp.depth, anti-Borrman:  ah_anti',a,'=',g12.5/
     *  '  Laue extinction length (um):   ext_L',a,'=',g12.5)
  686   format(
     *  ' Symmetric Bragg case parameters:'/
     *  '  Bragg extinction length (um):  ext_B',a,'=',g12.5/
     *  '  Bragg-case double-crystal curve:'/
     *  '         -- FWHM in arcsec: FWHM2asec_',a,'=',g12.5/
     *  '         -- FWHM in urad:   FWHM2urad_',a,'=',g12.5/
     *  '  Bragg case Darwin curve:'/
     *  '          -- FWHM in arcsec: FWHMasec_',a,'=',g12.5/
     *  '          -- FWHM in urad:   FWHMurad_',a,'=',g12.5/
     *  ' Bragg reflection parameters help: ',a)
c 701   format(' </pre>')
        if (modeout.ne.0) Then !-------+
          write(lunout,75,err=14)     !|
          write(lundup,75,err=14)     !|
        endif !------------------------+

        if     (ifail.eq.8 .and. d8.gt.0.)     Then !-----------------+
          i = Len_Trim(Errors(ifail+1))                              !|
          write(lunout,19,err=14) Errors(ifail+1)(1:i),wave/(2*d8)   !|lambda/2d>1
          write(lundup,19,err=14) Errors(ifail+1)(1:i),wave/(2*d8)   !|
        elseif (ifail.gt.0 .AND. ifail.le.11)  Then !-----------------+only 11 error messages
          i = Len_Trim(Errors(ifail+1))                              !|are defined!
          write(lunout,9,err=14)  Errors(ifail+1)(1:i)               !|
          write(lundup,9,err=14)  Errors(ifail+1)(1:i)               !|
        endif  !------------------------------------------------------+
  19    format(
     *  ' <tr><td><b><i>Error: </i></b></td><td><b>',a,', lambda/2d=',
     *                                           g11.5,'</b></td></tr>')
  9     format(
     *  ' <tr><td><b><i>Error: </i></b></td><td><b>',a,'</b></td></tr>')

  999   continue
        write (lunout,98,err=14)
        write (lundup,98,err=14)
  98    format(' </table></center>')
c ccc*    ' <hr>')
c       if (isyngony.eq.8 .AND. icoway.ne.2) Then !---------------------+
        n0 = Max(Len_Trim(code),1)
        if (isyngony  .eq. 8       .AND.
     *      code(1:n0).ne.'none'   .AND.
     *      ifail     .eq. 0 )     Then !-------------------------------+
          write (lunout,152,err=14) code(1:n0)                         !|
          write (lundup,152,err=14) code(1:n0)                         !|
  152     format(
     *    ' <p><center><table border bgcolor="#ffffff"><tr><td>'/
     *    ' <font color="red"><b>A T T E N T I O N:</b></font> <br>'/
     *    ' <font color="black">'/
     *    ' The structures marked as <b><i>**Amorphous**</i></b> are'/
     *    ' not necessary really amorphous. It only means that the'/
     *    ' <b><i>X0h</i></b> database does not have the data'/
     *    ' on the lattice constants and the unit cell coordinates'/
     *    ' of atoms for this material. Please, submit me these data'/
     *    ' if you want to use <b><i>X0h</i></b> for calculating the '/
     *    ' Bragg reflections from <b><i>',a,'</i></b>.'/
     *    ' </font>'/
     *    ' </td></tr></table></center>')                              !|
        endif  !--------------------------------------------------------+

        if (Len_Trim(referer).eq.0)
     *              referer = '/cgi/www_form.pl?template=x0h_form.htm'

        Call Navigate_Site(lunout,referer,i)
        if (i.ne.0) goto 14
        write (lunout,99,err=14)
        write (lundup,99,err=14)
  99    format(' </body></html>'/)

        if (radiat.eq.'none') radiat=' '
  11    continue
        Call exit_quiet()
c==================================================================
  4     continue
        l = Max (Len_Trim(outfile),1)
        write (*,12) outfile(1:l)
        write (0,12) outfile(1:l)
  12    format (' *** Error opening output: ',a)
        goto 11

  14    continue
        l = Max (Len_Trim(outfile),1)
        write (*,15) outfile(1:l)
        write (0,15) outfile(1:l)
  15    format (' *** Error writing output: ',a)
        goto 11

        end

c=================================================================

        Subroutine keys2x0h (nkeys, keys, vals,                 !inputs
     *                       icoway, code, rho,                 !outputs
     *                       ixway, radiat, wave, energy,       !outputs
     *                       indices, iref,                     !outputs
     *                       nHenkeCowan,                       !outputs
     *                       modeout, mdetail, ifail)           !outputs
        Real*4          wave, energy, rho
        Integer         nkeys, icoway, ixway, nHenkeCowan,
     *                  indices(3), iref, modeout, mdetail,
     *                  ifail, i, m,
     *                  lcrst, lamor, lchem

        Character       keys(nkeys)*(*), vals(nkeys)*(*),
     *                  code*(*), radiat*(*), valbuf*32,
     *                  crst*20, amor*20, chem*20, ii*2
        Logical         needRho

        Real      wave2energy
        External  wave2energy

c ' '                                                        ! ifail=0
c 'No material name specified'                               ! ifail=1
c 'Unexpected type of x-ray input'                           ! ifail=2
c 'Wrong or missing x-ray wavelength/energy'                 ! ifail=3
c 'Missing name of characteristic x-ray line'                ! ifail=4
c 'Wrong or missing Bragg reflection indices (hkl)'          ! ifail=5
c 'No x-ray data specified'                                  ! ifail=6
c 'X-ray wavelength not in the range 0.1-10A'                ! ifail=7
c 'Requested reflection does not exist for this wavelength'  ! ifail=8
c 'Material density not specified'                           ! ifail=9
c 'Unexpected type DB option'                                ! ifail=10
c 'Ambiguous material specification'                         ! ifail=11
c 'Incorrect X0h input parameter'                            ! ifail=12
        code        = ' '
        crst        = ' '
        amor        = ' '
        chem        = ' '
        radiat      = ' '
        rho         = 0.
        wave        = 0.
        energy      = 0.
        ixway       = 0         !1=wave,2=energy,3=line
        iref        = 0         !0/1 reflections
        modeout     = 0         !0=html, 1=text
        mdetail     = 0         !0=no, 1=yes
        icoway      = -1        !0-crystal,1-amorphous,2-chem.formula, -1=auto
        nHenkeCowan = -1        !-1=Auto,0=X0h_DB,1=Henke(f1),2=Henke(f1,f2),
                                !3=Cowan(f1),4=Cowan(f1,f2),10=compare all
        needRho     = .False.
        do      i=1,3 !===+
          indices(i) = 0 !|
        enddo  !==========+
        ifail      = 0
c Query String:
c code=Silicon xway=1 wave=1.54 line=Cu-Ka1 i1=7 i2=8 i3=9 chem=Al2O3 coway=0 rho= modeout=0 detail=0

        Call getbykey(nkeys,keys,vals,'modeout',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !------+
          Call Rdint(modeout,1,valbuf,i)        !|
          if (i.ne.0       .OR.                 !|
     *        modeout.lt.0 .OR.                 !|0=html
     *        modeout.gt.1) Then !-----+         |1=txt
            modeout = -1              !|         |
            ifail  = 12               !|         |Incorrect X0h input parameter
          endif  !---------------------+         |
        endif  !---------------------------------+

        Call getbykey(nkeys,keys,vals,'detail',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !------+
          Call Rdint(mdetail,1,valbuf,i)        !|
          if (i.ne.0       .OR.                 !|
     *        mdetail.lt.0 .OR.                 !|0=no
     *        mdetail.gt.1) Then !-----+         |1=yes
            mdetail = -1              !|         |
            ifail  = 12               !|         |Incorrect X0h input parameter
          endif  !---------------------+         |
        endif  !---------------------------------+

        Call getbykey(nkeys,keys,vals,'coway',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !------+
          Call Rdint(icoway,1,valbuf,i)         !|
          if (i.ne.0       .OR.                 !|
     *        icoway.lt.-1 .OR.                 !|-1=auto, 0=crystal,
     *        icoway.gt.2) Then !------+         |1=amorphous, 2-formula
            icoway = -1               !|         |
            ifail  = 2                !|         |Unexpected type of x-ray input
          endif  !---------------------+         |
        endif  !---------------------------------+

c Take crystal code selected from a list:
        Call getbykey(nkeys,keys,vals,'code',crst)
c Take amorphous code selected from a list:
        Call getbykey(nkeys,keys,vals,'amor',amor)
c Interpret code entered as a chemical formula:
        Call getbykey(nkeys,keys,vals,'chem',chem)

        if     (icoway.eq.0) Then  !--------------------------+Crystal!
          code = crst                                        !|
        elseif (icoway.eq.1) Then  !--------------------------+Amorphous material!
          code = amor                                        !|
        elseif (icoway.eq.2) Then  !--------------------------+Chemical formula!
          code = chem                                        !|
          needRho = .True.                                   !|
        elseif (icoway.eq.-1) Then !--------------------------+auto
          lcrst = Len_trim(crst)                             !|
          lamor = Len_trim(amor)                             !|
          lchem = Len_trim(chem)                             !|
          if     (lcrst.gt.0) Then !---------------------+    |
            code = crst                                 !|    |
            if (lamor.gt.0 .and. amor.ne.code) ifail=11 !|    |Ambigous material specification
            if (lchem.gt.0 .and. chem.ne.code) ifail=11 !|    |Ambigous material specification
          elseif (lamor.gt.0) Then !---------------------+    |
            code = amor                                 !|    |
            if (lchem.gt.0 .and. chem.ne.code) ifail=11 !|    |Ambigous material specification
          elseif (lchem.gt.0) Then !---------------------+    |
            code = chem                                 !|    |
            needRho = .True.                            !|    |
          else !-----------------------------------------+    |missing keyword for material
            ifail = 1                                   !|    |No material name specified
          endif !----------------------------------------+    |
        else  !-----------------------------------------------+keyword not in range for material
          ifail = 1                                          !|No material name specified
        endif  !----------------------------------------------+

        if (Len_Trim(code).eq.0)        Then  !--+missing keyword for material
          ifail = 1                             !|No material name specified
        endif  !---------------------------------+

        if (needRho) Then  !---------------------------+Chemical formula!
          Call getbykey(nkeys,keys,vals,'rho',valbuf) !|
          if (Len_Trim(valbuf) .eq. 0) Then  !------+  |
c This control can be removed: first,               |  |
c density is checked and properly reported          |  |
c by X0h1; second, user may specify a valid         |  |
c crystal or atom here and then density is          |  |
c not required /SS 2006.01.05/:                     |  |No density for
c           if (ifail.eq.0)  ifail = 9             !|  |chemical formula
          else  !-----------------------------------+  |
            Call Rdreal (rho,1,valbuf,i)           !|  |
            if (i.ne.0 .OR. rho.le.0.) Then !---+   |  |
c This control can be removed: first,           |   |  |
c density is checked and properly reported      |   |  |
c by X0h1; second, user may specify a valid     |   |  |
c crystal or atom here and then density is      |   |  |
c not required /SS 2006.01.05/:                 |   |  |No density for
c             if (ifail.eq.0)  ifail = 9       !|   |  |chemical formula
            endif  !----------------------------+   |  |
          endif  !----------------------------------+  |
        endif  !---------------------------------------+

        Call getbykey(nkeys,keys,vals,'xway',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !------+
          Call Rdint(ixway,1,valbuf,i)          !|
          if (i.ne.0       .OR.                 !|
     *        ixway.lt.1   .OR.                 !|1=wave,2=energy,3=line
     *        ixway.gt.3) Then  !-------+        |
            ixway = 0                  !|        |
            if (ifail.eq.0)  ifail = 2 !|        |Unexpected type of x-ray input
          endif  !----------------------+        |
        else !-----------------------------------+xway keyword problem
          ixway = 0                             !|
          if (ifail.eq.0)  ifail = 2            !|No x-ray data specified
        endif  !---------------------------------+

        if (ixway.eq.1 .OR. ixway.eq.2) Then  !---------------------------+wave or energy
                                                                         !|
          radiat = 'none'                                                !|
                                                                         !|
          Call getbykey(nkeys,keys,vals,'wave',valbuf)                   !|
          if (Len_Trim(valbuf) .gt. 0) Then !-------+                     |
            Call Rdreal (wave,1,valbuf,i)      !|   |                     |
            if (i.ne.0 .OR. wave.lt.0.) Then  !-+   |                     |
              if (ifail.eq.0)  ifail = 3       !|   |                     |Wrong or missing x-ray wavelength/energy
            endif  !----------------------------+   |                     |
          else  !-----------------------------------+wave keyword problem |
            if (ifail.eq.0)  ifail = 3             !|                     |Wrong or missing x-ray wavelength/energy
          endif  !----------------------------------+                     |
                                                                         !|
          if (ixway.eq.1)  Then  !----------------------------------+     |
            if (abs(wave).gt.1.E-33) Then !---------------+         |     |
              energy = wave2energy(wave)                 !|         |     |
c             if (wave.lt.0.1 .or. wave.gt.10.) Then  !-+ |         |     |X-ray wavelength not in the range 0.1-10A
              if (wave.lt.0.) Then  !-------------------+ |         |     |
                if (ifail.eq.0) ifail = 3              !| |         |     |Wrong or missing x-ray wavelength/energy
              endif  !----------------------------------+ |         |     |
            else  !---------------------------------------+         |     |
              if (ifail.eq.0)   ifail = 6                !|wave=0   |     |No x-ray data specified
            endif  !--------------------------------------+         |     |
          else  !--(ixway.eq.2)-------------------------------------+     |
            energy = wave                                          !|     |
            wave   = 0.                                            !|     |
            if (abs(energy).gt.1.E-33) Then  !------------+         |     |
              wave = wave2energy(energy)                 !|         |     |
c             if (wave.lt.0.1 .or. wave.gt.10.) Then  !-+ |         |     |X-ray wavelength not in the range 0.1-10A
              if (wave.lt.0.) Then  !-------------------+ |         |     |
                if (ifail.eq.0) ifail = 3              !| |wave not |     |Wrong or missing x-ray wavelength/energy
              endif  !----------------------------------+ |in range |     |
            else  !---------------------------------------+         |     |
              if (ifail.eq.0)   ifail = 6                !|wave=0   |     |No x-ray data specified
            endif  !--------------------------------------+         |     |
          endif  !--------------------------------------------------+     |
                                                                         !|
        elseif (ixway.eq.3)  Then  !--------------------------------------+x-ray line
                                                                         !|
          Call getbykey(nkeys,keys,vals,'line',radiat)                   !|
          if (Len_Trim(radiat) .eq. 0) Then !---+line keyword problem     |
            if (ifail.eq.0)  ifail = 4         !|                         |Missing name of characteristic x-ray line
            radiat = 'none'                    !|                         |
          endif  !------------------------------+                         |
                                                                         !|
          energy = 0.                                                    !|
          wave   = 0.                                                    !|
          if ((radiat.eq.'none' .OR. Len_Trim(radiat).eq.0)              !|
     *       .AND. ifail.eq.0) ifail = 6                         !wave=0  |No x-ray data specified
                                                                         !|
        else  !------------- (ixway.eq.0) --------------------------------+incorrect xway
                                                                         !|
          wave   = 0.                                                    !|
          radiat = 'none'                                                !|
          if (ifail.eq.0)  ifail = 2                                     !|Unexpected type of x-ray input
                                                                         !|
        endif  !----------------------------------------------------------+

        ii= 'i_'
        do      m=1,3  !==============================+
          write(ii(2:2),'(i1)') m                    !|
          Call getbykey(nkeys,keys,vals,ii,valbuf)   !|
          if (Len_Trim(valbuf) .gt. 0) Then !----+    |
            Call Rdint(indices(m),1,valbuf,i)   !|    |
            if (i.ne.0) Then  !------------+     |    |index keyword problem
              indices(m) = 0              !|     |    |
              if (ifail.eq.0)  ifail = 5  !|     |    |Wrong or missing Bragg reflection indices (hkl)
            endif  !-----------------------+     |    |
          endif  !-------------------------------+    |
          if (indices(m).ne.0) iref = 1              !|
        enddo  !======================================+
c Do not calculate reflections for amorphous materials:
        if (icoway.ne.0) iref = 0

        Call getbykey(nkeys,keys,vals,'df1df2',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !-----+
          Call Rdint(nHenkeCowan,1,valbuf,i)   !|
          if (i.ne.0             .OR.          !|
     *        (nHenkeCowan.ne.-1 .AND.         !|-1=Auto            DB
     *         nHenkeCowan.ne.0  .AND.         !| 0=X0h             DB
c    *         nHenkeCowan.ne.1  .AND.         !| 1=Henke   (f1)    DB
     *         nHenkeCowan.ne.2  .AND.         !| 2=Henke   (f1,f2) DB
c    *         nHenkeCowan.ne.3  .AND.         !| 3=Cowan   (f1)    DB
     *         nHenkeCowan.ne.4  .AND.         !| 4=Cowan   (f1,f2) DB
c    *         nHenkeCowan.ne.5  .AND.         !| 5=Windt   (f1)    DB
     *         nHenkeCowan.ne.6  .AND.         !| 6=Windt   (f1,f2) DB
c    *         nHenkeCowan.ne.7  .AND.         !| 7=Chantler(f1)    DB
     *         nHenkeCowan.ne.8  .AND.         !| 8=Chantler(f1,f2) DB
     *         nHenkeCowan.ne.10)) Then !--+    |10=compare all: 0,2,4,6,8
            nHenkeCowan = 0               !|    |
            if (ifail.eq.0)  ifail = 10   !|    |Unexpected type DB option
          endif  !-------------------------+    |
        endif  !--------------------------------+

        return
        end

c=================================================================

        Subroutine Borrmann (ah, xi0, xih, xdf,
     *                         ah_brmn, ah_anti)
        Real*4  pi, ah, xi0, xih, xdf, ah_brmn, ah_anti, tmp

        pi = 4.*atan(1.)                !the "pi" constant

        ah_brmn = 0.
        ah_anti = 0.
        if (abs(xi0) .lt. 1.E-33) return
        tmp = abs(xih/xi0) * abs(cos(pi*xdf))
        if (tmp .lt. 1.0)  Then !----+
           ah_brmn = ah / (1.-tmp)  !|
        else !-----------------------+
          ah_brmn  = 9.9999E+36     !|
        endif !----------------------+
        ah_anti    = ah / (1.+tmp)
        return
        end
