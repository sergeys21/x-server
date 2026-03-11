        Program x0p_web
c **********************************************************************
c * This program reads data from command line in the key=data format   *
c * with the keys from the X0p (X0h+) request web page at X-ray Server *
c * and makes search for the Bragg planes under specified search       *
c * conditions. The results are printed to STDOUT as a web page.       *
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
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Integer   nfoundMax
        Parameter (nfoundMax = 1331)            ! 11*11*11

        Real*4    Lattice_constants(6), rho, wave, energy,
     *            Q_limits(2), QB_limits(2), PRC_min,
     *            Q_found(nfoundMax), QB_found(nfoundMax),
     *            PRC_found(nfoundMax), xr0, xi0, xxx, qb
     *
        Integer   nfound, hkl_base(4),
     *            hkl_limits(4,2),
     *            hkl_found(4,nfoundMax),
     *            isyngony,
     *            i, j, m,
     *            ixway,
     *            mode_search,
     *            ifail, jpr, inmax, icall,
     *            n16, n17, n18, n19,
     *            lunout, lundup,
     *            lhost, laddr, lout, lurl,
     *            io_status

        Character txcount*8, outfile*32, url*64, query*512,
     *            hostname*256, address*60, referer*160, head*6,
c    *            quest*128,
     *            Symmetry(9)*24 /
     *                'Unknown',                                        !0
     *                'Cubic',                                          !1
     *                'Hexagonal/Trigonal',                             !2
     *                'Tetragonal',                                     !3
     *                'Trigonal(Rhombohedral)',                         !4
     *                'Orthorhombic',                                   !5
     *                'Monoclinic',                                     !6
     *                'Triclinic',                                      !7
     *                '** Amorphous ** '/,                              !8
     *            FINDtxt(3)*48 /
     *                'From (BraggAngle-Theta1) to (BraggAngle-Theta2)',!1
     *                'From Theta1 to (BraggAngle-Theta2)',             !2
     *                'From Theta1 to Theta2'/                          !3

c X0h (International Tables) range:       0.56  --    2.29 Angstrem
c             (5,000 eV --  25,000 eV):   0.50  --    2.47 Angstrem
c Henke    range (10 eV --  30,000 eV):   0.41  -- 1239.81 Angstrem
c Cowan    range (30 eV -- 694,500 eV):   0.02  --  413.27 Angstrem
c Windt    range (10 eV -- 100,000 eV):   0.12  -- 1239.81 Angstrem
c Chantler range (10 eV -- 450,000 eV):   0.28  -- 1239.81 Angstrem
        Character DBtext(10)*45 /
     *                'Automatic DB selection',                         !-1 (1)
     *                'X0h (International Tables), 5-25 KeV',           !0  (2)
     *                'Henke Tables, 0.01-30 KeV (df1 only)',           !1  (3)
     *                'Henke Tables, 0.01-30 KeV',                      !2  (4)
     *                'Brennan-Cowan Tables, 0.03-700 KeV (df1 only)',  !3  (5)
     *                'Brennan-Cowan Tables, 0.03-700 KeV',             !4  (6)
     *                'Windt Tables, 0.01-100 KeV (df1 only)',          !5  (7)
     *                'Windt Tables, 0.01-100 KeV',                     !6  (8)
     *                'Chantler/NIST Tables, 0.01-450 KeV (df1 only)',  !7  (9)
     *                'Chantler/NIST Tables, 0.01-450 KeV'/             !8  (10)

        Integer   Nkeys
        Parameter (Nkeys=26)

        Character keys(Nkeys)*10 /'xway', 'wave', 'line', 'code',
     *                            'hkl11',  'hkl12', 'hkl13',
     *                            'hkl21',  'hkl22', 'hkl23',
     *                            'base1',  'base2', 'base3',
     *                            'qb1', 'qb2', 'q1', 'q2',
     *                            'modesearch', 'prcmin', 'df1df2',
     *                            'detail', 'modeout',
     *                            'job', 'ip', 'host', 'referrer'/,
     *            vals(Nkeys)*256

        Character tnb*14 /' class="novrb"'/,
     *            ttb*14 /' class="nobtm"'/
c    *            tbb*14 /' class="notop"'/

        Real      wave2energy
        External  wave2energy

        Character Get_Server_URL*64
        External  Get_Server_URL

        Character code*20, radiat*6
        Common  /x0pa1/ radiat

        Character txt(20)*80
        Common  /msg/   txt

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
        progname = 'x0p_form'
        ErrorFile= 'x0p_form.err'
        modebat  = 1
        lunbat   = 0
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
c       quest = '<a  href="javascript:void(0)" onClick="Wfloat'//
c    *          '(''/images/x0p_help.gif'',''x0p_help'',740,357);">'//
c    *          '<b>?</b></a> &nbsp;'

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
        if (ifail .ne. 0) Then !--------------+
          txt(1) = 'Error parsing cmdline'   !|
          Call HTMLstop (lunout, 1)          !|
          stop 'x0p_web: emergency exit'     !|
        endif !-------------------------------+

        Call OpenFile(txcount//'.htm',lundup,
     *                'readwrite','unknown',
     *                io_status,*11)

        write (lunout,511,err=11)
c 511   format(' Content-type: text/html'/)             !before 2013.11.01
  511   format('Content-type: text/html'/)
        head = '<head>'
        write (lunout,'(a//a/a)') '<!DOCTYPE html>','<html>',head
        Call Insert_UrlBase(head,lunout)
        write (lunout,512,err=11)
  512   format(
     *  ' <title>X0p: Bragg Planes Search Results</title>'/
     *  ' <link rel="stylesheet" href="style.css" type="text/css">'/
     *  ' </head><body>')

        Call mergeQuery(nkeys-4,keys,vals,query)
        j = Max(Len_Trim(query),1)
        url = Get_Server_URL()
        lurl = Max(Len_Trim(url),1)
        write (lundup,513,err=11) url(1:lurl),'cgi/x0p_form.pl',
     *                            query(1:j)
  513   format(' <html>'/
     *  ' <head><title>X0h+ Results</title></head>'/
     *  ' <body text="#000000" bgcolor="silver">'/
     *  ' <!-- ',2a,'?',a,' -->'/
     *  ' <table border=0><tr valign="top">')
c-------
        write (lunout,50,err=11)
  50    format(' <table border="0"><tr valign="top"><td>'/
     *         '&nbsp; <img src="images/x0h_search_th_th.jpg"',1x,
     *                  'align="absmiddle" border="0">&nbsp;</td>')
c-------
        if (hostname .ne. address) Then  !---------------------+
          do i=lunout,lundup !================================+ |
            write (i,51,err=11) txcount,                     !| |
     *                          ' for ',address(1:laddr),    !| |
     *                          ' [',hostname(1:lhost),']:'  !| |
            Call flushoutput(i)                              !| |
          enddo !=============================================+ |
        else  !-------------------------------------------------+
          do i=lunout,lundup !================================+ |
            write (i,51,err=11) txcount,                     !| |
     *                          ' for ',address(1:laddr),    !| |
     *                          ' :',' ',' '                 !| |
            Call flushoutput(i)                              !| |
          enddo !=============================================+ |
        endif  !------------------------------------------------+
  51    format(' <td><font size="+1">'/
     *         ' &nbsp; Job ID: &nbsp;<b>',a,'</b><br>'/
     *         ' &nbsp; <b><i>X0h+</i> results',5a,'</b>'/
     *         ' </font></td></tr></table><br>')

        do i=lunout,lundup !=======+
          write (i,52,err=11)     !|
          Call flushoutput(i)     !|
        enddo !====================+
c Allowed colors:
c --------------
c Black, Maroon, Green, Olive, Navy, Purple, Teal, Gray,
c Silver,Red, Lime, Yellow, Blue, Fuchsia, Aqua, White
c
c <table border bgcolor="#d0f0f0"> -- light blue color!
c <table border bgcolor="#f0f0f0"> -- not bright white color!

  52    format(
     *  ' <center><table border bgcolor="#ffffff" class="vtop">')

        Call keys2x0p (nkeys-4, keys, vals,
     *                 code, iHenkeCowan,
     *                 ixway, wave, energy, radiat,
     *                 mode_search, hkl_base, hkl_limits,
     *                 Q_limits, QB_limits, PRC_min, ifail)
        if (ifail.ne.0) goto 998 !------------------+
                                                   !v
        do      i=1,6  !===============+
          Lattice_Constants(i) = 0.   !|
        enddo  !=======================+
        rho      = 0.
        isyngony = 0
        icall    = 0
        jpr      = 0
c       jpr      = 2
c Calculate X0:
        Call X0h1 (Lattice_constants, wave,
     *             code, rho, 0,
     *             0, 0, 0,
     *             isyngony, icall,
     *             jpr, qb, xr0, xi0,
     *             xxx, xxx, xxx,
     *             xxx, xxx, xxx, ifail)
        if (isyngony.eq.8) Then  !--------------------------------+
          i = Len_Trim(code)                                     !|
          write(txt,1) code(1:i)                                 !|
  1       format('The structure <b>',a,'</b> is not a crystal.'/
     *    'No Bragg planes can be found.')                       !|
           ifail = 2                                             !|
        endif !---------------------------------------------------+
        if (ifail.ne.0) goto 998 !------------------+
                                                   !v

        if (Len_Trim(radiat).eq.0) radiat = 'none'
        energy = wave2energy(wave)


        if (isyngony.eq.2)  Then  !-----------------------------+
          hkl_base(4)     = hkl_base(3)                        !|
          hkl_base(3)     = -(hkl_base(1)+hkl_base(2))         !|
          hkl_limits(4,1) = hkl_limits(3,1)                    !|
          hkl_limits(3,1) = -(hkl_limits(1,1)+hkl_limits(2,1)) !|
          hkl_limits(4,2) = hkl_limits(3,2)                    !|
          hkl_limits(3,2) = -(hkl_limits(1,2)+hkl_limits(2,2)) !|
          inmax  = 4                                           !|
        else  !-------------------------------------------------+
          inmax  = 3                                           !|
        endif  !------------------------------------------------+

        Call X0h_find (Code, Lattice_constants, isyngony, wave,
     *                 mode_search, hkl_base, hkl_limits,
     *                 Q_limits, QB_limits, PRC_min,
     *                 nfoundMax, nfound, hkl_found,
     *                 Q_found, QB_found, PRC_found, ifail)
        if (ifail.ne.0) goto 998 !------------------+
                                                   !v

        Call PakInt (hkl_base,       inmax,txt(17),n17)
        Call PakInt (hkl_limits(1,1),inmax,txt(18),n18)
        Call PakInt (hkl_limits(1,2),inmax,txt(19),n19)
        do j=lunout,lundup  !==========================================+
          write (j,720,err=11)                                        !|
     *                     tnb, tnb, Code,                            !|
     *                     tnb, tnb, Symmetry(isyngony+1),            !|
     *                     tnb, tnb, rho,                             !|
     *                     tnb, tnb, (Lattice_constants(i), i=1,3),   !|
     *                     tnb, tnb, (Lattice_constants(i), i=4,6),   !|
     *                     tnb, tnb, wave,                            !|
     *                     tnb, tnb, energy,                          !|
     *                     tnb, tnb, radiat,                          !|
     *                     tnb, tnb, txt(18)(1:n18), txt(19)(1:n19),  !| hkl_limits(i,1),i=1,inmax, hkl_limits(i,2),i=1,inmax
     *                     tnb, tnb, QB_limits,                       !|
     *                     tnb, tnb, PRC_min,                         !|
     *                     tnb, tnb, txt(17)(1:n17),                  !| hkl_base(i),i=1,inmax
     *                     tnb, tnb, FINDtxt(mode_search),            !|
     *                     tnb, tnb, Q_limits                         !|
          if (PRC_min.gt.0.) Then !---------------------------+        |
            m = Len_Trim(DBtext(iHenkeCowan+2))              !|        |
            write (j,721,err=11) tnb, tnb,                   !|        |
     *                           DBtext(iHenkeCowan+2)(1:m)  !|        |
          endif  !--------------------------------------------+        |
          write (j,722,err=11) ttb                                    !|
          if (nfound.gt.0)      Then !------------------------------+  |
            write (j,723,err=11) ttb,nfound,min0(nfound,nfoundMax) !|  |
            if (PRC_min.gt.0.) Then !-----------------------------+ |  |
              write (j,724,err=11)                               !| |  |
              do m=1,min0(nfound,nfoundMax) !===================+ | |  |
                Call PakInt (hkl_found(1,m),inmax,txt(16),n16) !| | |  |
                write(j,725) txt(16)(1:n16),                   !| | |  |
     *                       Q_found(m), QB_found(m),          !| | |  |
     *                       PRC_found(m)                      !| | |  |
              enddo  !==========================================+ | |  |
            else !------------------------------------------------+ |  |
              write (j,726,err=11)                               !| |  |
              do m=1,min0(nfound,nfoundMax) !===================+ | |  |
                Call PakInt (hkl_found(1,m),inmax,txt(16),n16) !| | |  |
                write(j,727) txt(16)(1:n16),                   !| | |  |
     *                       Q_found(m), QB_found(m)           !| | |  |
              enddo  !==========================================+ | |  |
            endif !-----------------------------------------------+ |  |
          else !----------------------------------------------------+  |
            write (j,728,err=11)                                   !|  |
          endif !---------------------------------------------------+  |
          write (j,729,err=11)                                        !|
        enddo  !=======================================================+

  720   format(
     *  ' <tr><th colspan="2"> search conditions : </th></tr>'/
     *  ' <tr><td',a,'> Crystal:',1x,                 ' </td>',
     *       '<td',a,'> ',          a,                ' </td></tr>'/
     *  ' <tr><td',a,'> Symmetry group:',1x,          ' </td>',
     *       '<td',a,'> ',          a,                ' </td></tr>'/
     *  ' <tr><td',a,'> Density (gm/cm<sup>3</sup>) :   </td>',
     *       '<td',a,'> ',g12.5,    ' <sup>&nbsp;</sup> </td></tr>'/
     *  ' <tr><td',a,'> Unit cell constants (A)  :      </td>',
     *       '<td',a,'> ',2(g12.5,' , '),g12.5,       ' </td></tr>'/
     *  ' <tr><td',a,'> Unit cell angles (degr.) :      </td>',
     *       '<td',a,'> ',2(g12.5,' , '),g12.5,       ' </td></tr>'/
     *  ' <tr><td',a,'> X-ray wavelength (Angstrom) :   </td>',
     *       '<td',a,'> ',      g14.7                 ' </td></tr>'/
     *  ' <tr><td',a,'> X-ray energy (keV) :',1x,     ' </td>',
     *       '<td',a,'> ',      g14.7                 ' </td></tr>'/
     *  ' <tr><td',a,'> X-ray characteristic line :     </td>',
     *       '<td',a,'> ',        a,                  ' </td></tr>'/
     *  ' <tr><td',a,'> Bragg planes range :',1x,     ' </td>',
     *       '<td',a,'> (',    a,   ') -- (',   a,   ') </td></tr>'/
     *  ' <tr><td',a,'> Bragg angles range :',1x,     ' </td>',
     *       '<td',a,'> ',   f9.4,   ' -- ',   f9.4,  ' </td></tr>'/
     *  ' <tr><td',a,'> Minimum intensity (xh/x0) :     </td>',
     *       '<td',a,'> ',         f9.4,'%',1x,       ' </td></tr>'/
     *  ' <tr><td',a,'> Surface :',1x,                ' </td>',
     *       '<td',a,'> (',            a,            ') </td></tr>'/
     *  ' <tr><td',a,'> Planes angles to surface :      </td>',
     *       '<td',a,'> ',             a,             ' </td></tr>'/
     *  ' <tr><td',a,'> Theta1, Theta2 :',1x,         ' </td>',
     *       '<td',a,'> ',   f9.4,   ' -- ',   f9.4,  ' </td></tr>')
  721   format(
     *  ' <tr><td',a,'> Dispersion corrections DB:      </td>',
     *       '<td',a,'> ',              a,        '      </td></tr>')
  722   format(' <tr><th',a,' colspan="2"> SEARCH RESULTS : </th></tr>')
  723   format(' <tr><th',a,' colspan="2"> Planes found: ',i4,
     *              '. &nbsp; Planes being displayed: ',i4,'</th></tr>'/
     *  ' <tr><td colspan="2">',
     *   '<table width="100%" border="0" class="vtop">')
  724   format(' <tr><th>hkl</th>'/
     *             ' <th>Angle to surface</th>'/
     *             ' <th>Bragg angle</th>'/
     *             ' <th>Relative Intensity <br>xh/x0(%)</th></tr>')
  725   format(' <tr><td align="center"> (', a,  ') </td>'/
     *             ' <td align="center"> ', f9.4, ' </td>'/
     *             ' <td align="center"> ', f9.4, ' </td>'/
     *             ' <td align="center"> ', f9.4, ' </td></tr>')
  726   format(' <tr><th>hkl</th>'/
     *         ' <th>Angle to surface</th>'/
     *         ' <th>Bragg angle</th></tr>')
  727   format(' <tr><td align="center"> (', a,  ') </td>'/
     *             ' <td align="center"> ', f9.4, ' </td>'/
     *             ' <td align="center"> ', f9.4, ' </td></tr>')
  728   format(' <tr><th colspan="2">No Bragg planes found to satisfy',
     *                              ' the above conditions!</th></tr>')
  729   format(' </table></td></tr>')

  999   continue
        write (lunout,98,err=11)
        write (lundup,98,err=11)
  98    format(' </table></center>')

        if (Len_Trim(referer).eq.0)
c    *              referer = '/x0p_form.shtml'
     *              referer = '/cgi/www_form.pl?template=x0p_form.htm'

        Call Navigate_Site(lunout,referer,i)
        if (i.ne.0) goto 11

        write (lunout,99,err=11)
        write (lundup,99,err=11)
  99    format(' </body></html>'/)

  11    continue
        Call exit_quiet()
c==================================================================
  998   continue
        ifail = Iabs(ifail)
        Call HTMLstop       (lunout,ifail)
        Call HTMLstop_noJPG (lundup,ifail)                               !^
        goto 999  !--------------------------------------------------+
        end

c=================================================================

        Subroutine keys2x0p (nkeys, keys, vals,                 !inputs
     *                       code, iHenkeCowan,                 !outputs
     *                       ixway, wave, energy, radiat,       !outputs
     *                       mode_search, hkl_base, hkl_limits, !outputs
     *                       Q_limits, QB_limits, PRC_min,      !outputs
     *                       ifail)                             !outputs

        Integer         hkl_MAX
        Parameter       (hkl_MAX = 10)

        Real            wave, energy, Q_limits(2), QB_limits(2),
     *                  PRC_min
        Integer         nkeys, iHenkeCowan, ixway, mode_search,
     *                  hkl_base(4), hkl_limits(4,2), ifail,
     *                  i, l, m, n
        Character       keys(nkeys)*(*), vals(nkeys)*(*),
     *                  code*(*), radiat*(*), valbuf*32, key*8

        Real      wave2energy
        External  wave2energy

        Character       txt(20)*80
        Common  /msg/   txt

        code        = ' '
        radiat      = ' '
        wave        = 0.
        energy      = 0.
        ixway       = 0
        iHenkeCowan = -1
        do i=1,3 !=============+
          hkl_base(i)     = 0 !|
          hkl_limits(i,1) = 0 !|
          hkl_limits(i,2) = 0 !|
        enddo  !===============+
        do i=1,2 !===========+
          Q_limits(i)  = 0. !|
          QB_limits(i) = 0. !|
        enddo  !=============+
        PRC_min     = 0.
        ifail       = 0

c Query String:
c code=Silicon xway=1 wave=1.54 line=Cu-Ka1 i1=7 i2=8 i3=9

c Take crystal code selected from a list:
        txt(3) = 'Crystal code'
        Call getbykey(nkeys,keys,vals,'code',code)
        if (Len_Trim(code) .eq. 0) goto 40

        txt(3) = 'Dispersion corrections database option'
        Call getbykey(nkeys,keys,vals,'df1df2',valbuf)
        if (Len_Trim(valbuf) .eq. 0) goto 40
        txt(2) = 'Parameter is not in range [-1::8]:'
        Call Rdint(iHenkeCowan,1,valbuf,i)
        if (i.ne.0             .OR.
     *      (iHenkeCowan.ne.-1 .AND.              !-1=Auto            DB
     *       iHenkeCowan.ne.0  .AND.              ! 0=X0h             DB
c    *       iHenkeCowan.ne.1  .AND.              ! 1=Henke   (f1)    DB
     *       iHenkeCowan.ne.2  .AND.              ! 2=Henke   (f1,f2) DB
c    *       iHenkeCowan.ne.3  .AND.              ! 3=Cowan   (f1)    DB
     *       iHenkeCowan.ne.4  .AND.              ! 4=Cowan   (f1,f2) DB
c    *       iHenkeCowan.ne.5  .AND.              ! 5=Windt   (f1)    DB
     *       iHenkeCowan.ne.6  .AND.              ! 6=Windt   (f1,f2) DB
c    *       iHenkeCowan.ne.7  .AND.              ! 7=Chantler(f1)    DB
     *       iHenkeCowan.ne.8)) goto 41           ! 8=Chantler(f1,f2) DB

        txt(3) = 'X-ray wavelength specification mode'
        Call getbykey(nkeys,keys,vals,'xway',valbuf)
        if (Len_Trim(valbuf) .eq. 0) goto 40
        txt(2) = 'Parameter is not in range [1::3]:'
        Call Rdint(ixway,1,valbuf,i)
c 1=wave 2=energy 3=line
        if (i.ne.0      .OR.
     *      ixway.lt.1  .OR.
     *      ixway.gt.3) goto 41

        if (ixway.eq.1 .OR. ixway.eq.2) Then !-------------+
c wave or energy:                                         !|
          radiat = 'none'                                 !|
          txt(3) = 'X-ray wavelength/energy value'        !|
          Call getbykey(nkeys,keys,vals,'wave',valbuf)    !|
          if (Len_Trim(valbuf) .eq. 0) goto 40            !|
          txt(2) = 'Parameter is zero or negative:'       !|
          Call Rdreal (wave,1,valbuf,i)                   !|
          if (i.ne.0)     goto 40                         !|
          if (wave.le.0.) goto 41                         !|
                                                          !|
          if (ixway.eq.1)  Then  !--------+                |
            energy = wave2energy(wave)   !|                |
          else  !--(ixway.eq.2)-----------+                |
            energy = wave                !|                |
            wave   = wave2energy(energy) !|                |
          endif  !------------------------+                |
                                                          !|
        else  !-----(ixway.eq.3)---------------------------+
c x-ray line:                                             !|
          txt(3) = 'Characteristic X-ray line'            !|
          Call getbykey(nkeys,keys,vals,'line',radiat)    !|
          if (Len_Trim(radiat) .eq. 0 .OR.                !|
     *        radiat.eq.'none')        goto 40            !|
                                                          !|
        endif  !-------------------------------------------+

        txt(3) = 'Bragg planes search mode'
        Call getbykey(nkeys,keys,vals,'modesearch',valbuf)
        if (Len_Trim(valbuf) .eq. 0) goto 40
        Call Rdint(mode_search,1,valbuf,i)
        if (i.ne.0)     goto 40
c 1. From QB-Theta1 to QB-Theta2:
c 2. From Theta1 to QB-Theta2:
c 3. From Theta1 to Theta2:
        txt(2) = 'Parameter is not in range [1::3]:'
        if (mode_search.lt.1  .or.
     *      mode_search.gt.3) goto 41

        txt(3) = 'Base plane indices hkl'
        key = 'base_'
        n   = Len_Trim(key)
        do  m=1,3  !=======================================+
          write(key(n:n),'(i1)') m                        !|
          Call getbykey(nkeys,keys,vals,key(1:n),valbuf)  !|
          if (Len_Trim(valbuf) .eq. 0) goto 40            !|
          Call Rdint (hkl_base(m),1,valbuf,i)             !|
          if (i .ne. 0) goto 40                           !|
        enddo  !===========================================+
        txt(2) = 'Zero Miller indices specified:'
        if (hkl_base(1) .eq. 0 .AND.
     *      hkl_base(2) .eq. 0 .AND.
     *      hkl_base(3) .eq. 0 ) goto 41

        txt(3) = 'hkl-limits for searching the Bragg planes'
        key = 'hkl__'
        n   = Len_Trim(key)
        do  m=1,3  !=====================================================+
          do  l=1,2  !=======================================+           |
            write(key(n-1:n),'(2i1)') l, m                  !|           |
            Call getbykey(nkeys,keys,vals,key(1:n),valbuf)  !|           |
            if (Len_Trim(valbuf) .eq. 0) goto 40            !|           |
            Call Rdint (hkl_limits(m,l),1,valbuf,i)         !|           |
            if (i .ne. 0) goto 40                           !|           |
          enddo  !===========================================+           |
          if (Iabs(hkl_limits(m,2)-hkl_limits(m,1)).gt.hkl_MAX) goto 50 !|
        enddo  !=========================================================+

        txt(3) = 'Theta-limits for searching the Bragg planes'
        key = 'q_'
        n   = Len_Trim(key)
        do  l=1,2  !=======================================+
          write(key(n:n),'(i1)') l                        !|
          Call getbykey(nkeys,keys,vals,key(1:n),valbuf)  !|
          if (Len_Trim(valbuf) .eq. 0) goto 40            !|
          Call Rdreal (Q_limits(l),1,valbuf,i)            !|
          if (i .ne. 0) goto 40                           !|
          txt(2) = 'Parameter is not in range [0::180]:  '!|
          if (Q_limits(l).lt.0. .or.                      !|
     *        Q_limits(l).gt.180.) goto 41                !|
        enddo  !===========================================+

        txt(3) = 'Bragg-angle limits for searching the Bragg planes'
        key = 'qb_'
        n   = Len_Trim(key)
        do  l=1,2  !=======================================+
          write(key(n:n),'(i1)') l                        !|
          Call getbykey(nkeys,keys,vals,key(1:n),valbuf)  !|
          if (Len_Trim(valbuf) .eq. 0) goto 40            !|
          Call Rdreal(QB_limits(l),1,valbuf,i)            !|
          if (i .ne. 0) goto 40                           !|
          txt(2) = 'Parameter is not in range [0::90]:'   !|
          if (QB_limits(l) .lt. 0. .or.                   !|
     *        QB_limits(l) .gt. 90.) goto 41              !|
        enddo  !===========================================+

        txt(3) = 'Minimum reflectivity for searching the Bragg planes'
        Call getbykey(nkeys,keys,vals,'prcmin',valbuf)
        if (Len_Trim(valbuf) .eq. 0) goto 40
        Call Rdreal(PRC_min,1,valbuf,i)
        if (i .ne. 0)      goto 40
        txt(2) = 'Parameter is negative:'
        if (PRC_min .lt. 0.) goto 41
        return
c------------
c ERRORS:
  40    continue
        txt(1) = 'keys2x0p:'
        txt(2) = 'Missing or incorrect data for'
        ifail  = 3
        return

  41    continue
        txt(1) = 'keys2x0p:'
        ifail  = 3
        return

  50    continue
        write (txt,1)   m,hkl_MAX
  1     format('keys2x0p: requested search range for Miller index-',i1/
     *         'exceeds the maximum allowed value of ',i3)
        ifail  = 2
        return

        end
