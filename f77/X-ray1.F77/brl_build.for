        Program brl_frm1
c **********************************************************************
c * This program reads data from command line in the key=data format   *
c * and builds multiple Bragg/Laue diffraction geometry. The command   *
c * line input is identical to filling BRL step-1 HTML form at X-ray   *
c * server. The output is printed to STDOUT in HTML format.            *
c *                                                                    *
c *                  Author: Sergey Stepanov                           *
c * 2026/01 Version: 1.0:   initial implementation                     *
c **********************************************************************
        Integer   nReflex1
        Parameter (nReflex1 = 5)                !max number of base reflections including (000)
                                                !nReflex1=5 corresponds to 4-wave diffraction

        Real*8    wave8, k0(3,2), k0_mod,
     *            h_vec(3,nReflex1),            !base reflections vectors
     *            h_mod(nReflex1),              !base reflections modules
     *            h1xh2(3), h1xh2_inv(3),
     *            h1xh2_mod, DisorNorm8_gra,
     *            k0_copl(3), wave_copl,
     *            Uncopl, c8, s8(2), angle8(2),
     *            Sg_unit(3,2), Pi_Unit(3,2),
     *            h_v(3), h_m,                  !h_v is h-vector, h_m is h-module
     *            alpha_range /1000./           !alpha/|x0| to be considered in Bragg condition

        Real*8    DVecSca2, DAngle2
        External  DVecSca2, DAngle2

        Real*4    Lattice_constants(6), wave, energy,
     *            xr0, xi0, xrh, xih, xqh, xrf, xif, xqf,
     *            DisorNorm, filter, rho, qb, w

        Real      wave2energy
        External  wave2energy

        Integer   indices(3,nReflex1), iii(4), indrange,
     *            nBaseNormal(3), nMiscutDirection(3),
     *            isyngony, i, j, k, l, m, n,
     *            ixway, ixw2, iunim,
     *            ifail, nref, iref, ifound, iwatch,
     *            jpr, inmax, icall, iwrk(1),
     *            lhost, laddr, ldisor, lrefl(3),
     *            lcod, lrad, lbase, lmisc, io_status,
     *            normal_inverted, h1xh2_inverted,
     *            found_two_k0, lunout, lundup, lout, lurl,
     *            ii/3/                 !Number of Bragg planes indices (3 or 4)

        Character code*20, txcount*8, url*64,
     *            query*256, inputstring*256, outfile*128,
     *            hostname*256, address*60, referer*160,
     *            wrk*24,
     *            base_txt*24, misc_txt*24, disor_txt*24,
     *            refl_txt(3)*24,
     *            Uni(5)*5  /
     *                  'degr.','min.',
     *                  'mrad.','sec.',
     *                  'urad.'/,
     *            Errors(12)*52 /
     *                  'No crystal name specified',                     !ifail=1
     *                  'Incorrect type of x-ray input',                 !ifail=2
     *                  'Wrong or missing x-ray wavelength/energy',      !ifail=3
     *                  'Missing name of characteristic x-ray line',     !ifail=4
     *                  'Wrong or missing reflection indices (hkl)',     !ifail=5
     *                  'Wrong index search range',                      !ifail=6
     *                  'Wrong X0h DB option',                           !ifail=7
     *                  'Wrong or missing surface normal indices',       !ifail=8
     *                  'Wrong surface miscut angle',                    !ifail=9
     *                  'Wrong surface miscut angle units',              !ifail=10
     *                  'Wrong surface miscut axis indices',             !ifail=11
     *                  'Incorrect intensity filter |xh/x0|'/,           !ifail=12
     *            Symmetry(9)*24 /
     *                  'Unknown',                                       !0
     *                  'Cubic',                                         !1
     *                  'Hexagonal/Trigonal',                            !2
     *                  'Tetragonal',                                    !3
     *                  'Trigonal(Rhombohedral)',                        !4
     *                  'Orthorhombic',                                  !5
     *                  'Monoclinic',                                    !6
     *                  'Triclinic',                                     !7
     *                  '** Amorphous ** '/,                             !8
     *            nobr*19 /'<span class="nobr">'/

        Integer   Nkeys
        Parameter (Nkeys=28)

        Character keys(Nkeys)*10 /'xway', 'wave', 'line',
     *                            'code', 'df1df2',
     *                            'n1',  'n2',  'n3',
     *                            'm1',  'm2',  'm3', 'miscut', 'unim',
     *                            'i11', 'i12', 'i13',
     *                            'i21', 'i22', 'i23', 'indrange',
     *                            'i31', 'i32', 'i33', 'filter',
     *                            'job', 'ip', 'host', 'referrer'/,
     *            vals(Nkeys)*256

c X0h (International Tables) range:
c             (5,000 eV --  25,000 eV):   0.50  --    2.47 Angstrem
c Henke    range (10 eV --  30,000 eV):   0.41  -- 1239.81 Angstrem
c Cowan    range (30 eV -- 694,500 eV):   0.02  --  413.27 Angstrem
c Windt    range (10 eV -- 100,000 eV):   0.12  -- 1239.81 Angstrem
c Chantler range (10 eV -- 450,000 eV):   0.28  -- 1239.81 Angstrem
        Character DBtext(10)*45 /
     *                  'Automatic DB selection',                        !-1 (1)
     *                  'X0h (International Tables), 5-25 KeV',          !0  (2)
     *                  'Henke Tables, 0.01-30 KeV (df1 only)',          !1  (3)
     *                  'Henke Tables, 0.01-30 KeV',                     !2  (4)
     *                  'Brennan-Cowan Tables, 0.03-700 KeV (df1 only)', !3  (5)
     *                  'Brennan-Cowan Tables, 0.03-700 KeV',            !4  (6)
     *                  'Windt Tables, 0.01-100 KeV (df1 only)',         !5  (7)
     *                  'Windt Tables, 0.01-100 KeV',                    !6  (8)
     *                  'Chantler/NIST Tables, 0.01-450 KeV (df1 only)', !7  (9)
     *                  'Chantler/NIST Tables, 0.01-450 KeV'/            !8  (10)

        Character Get_Server_URL*64
        External  Get_Server_URL
c -------------------------------------------------------
        Character radiat*6
        Common  /x0pa1/ radiat
c -------------------------------------------------------
        Real*8              Surface_Normal(3)
        Common /SurfNorm/   Surface_Normal
c -------------------------------------------------------
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
        progname = 'brl_build'
        ErrorFile= 'brl_build.err'
        modebat  = 1
        lunbat   = 0
        Call    DeleteErrorFile ()
c -------------------------------------------------------
        outfile  = '/dev/stdout'
        lout = Len_Trim(outfile)
        lunout   = 6
        lundup   = lunout + 1
        ifail    = 0

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
        if (Len_Trim(txcount).eq.0) txcount = 'brl00000'

c If call for Get_CmdLine above failed:
        if (ifail .ne. 0) Then !--------------+
          txt(1) = 'Error parsing cmdline'   !|
          Call HTMLstop (lunout, 1)          !|
          stop 'brl_build: emergency exit'   !|
        endif !-------------------------------+

        ErrorFile = txcount//'.err'
        Call    DeleteErrorFile ()

c Open duplicating HTML page (for LOG purposes):
        Call OpenFile(txcount//'.htm',lundup,
     *                'readwrite','unknown',
     *                io_status,*11)

c       write (lunout,'(1x,a/)',err=11) 'Content-type: text/html'       !before 2013.11.01
        write (lunout,'(a/)',   err=11) 'Content-type: text/html'
        write (lunout,'(a)',err=11) ' <html><head>'
        write (lundup,'(a)',err=11) ' <html><head>'
        Call Insert_UrlBase('<head>',lunout)
        Call Insert_UrlBase('<head>',lundup)
        write (lunout,511,err=11) '>'
        write (lundup,511,err=11) ' text="#000000" bgcolor="silver">'
  511   format(
     *  ' <title>BRL Results</title>'/
     *  ' <link type="text/css" rel="stylesheet" href="style.css">'/
     *  ' <script type="text/javascript" language="JavaScript"',
     *                                   ' src="brl2.js"></script>'/
     *  ' </head>'/
     *  ' <body',a)

        Call mergeQuery(nkeys-4,keys,vals,query)
        j = Max(Len_Trim(query),1)
        url = Get_Server_URL()
        lurl = Max(Len_Trim(url),1)
        write (lundup,513,err=11) url(1:lurl),'cgi/brl_frm1.pl',
     *                            query(1:j)
  513   format(' <!-- ',2a,'?',a,' -->')
c-------
        write (lunout,54,err=11) txcount, hostname(1:lhost)
        write (lundup,54,err=11) txcount, hostname(1:lhost)
  54    format(
     *  ' <b><i>Job ID: &nbsp; ',a,'</i></b><br>'/
     *  ' <b><i>User: &nbsp; &nbsp; &nbsp;',a,'</i></b><br>')

        write (lunout,50,err=11)
  50    format(' <h3>&nbsp;',1x,
     *    '<img src="images/brl_th_th.jpg"',1x,
     *    'border="0" width="152" height="32"',1x,
     *    'align="absmiddle" alt="BRL on the Web"> &nbsp;',1x,
     *    '<u><b><i>BRL</i></b> Request Form: Step-2/2</u></h3>')
        write (lundup,51,err=11)
  51    format(' <h3>&nbsp; <u><b><i>BRL</i></b> Request Form</u></h3>')
c-------

c+============================================================+
        Call keys2brl  (nkeys-4, keys, vals,                 !|
     *                  code,                                !|
     *                  nBaseNormal, nMiscutDirection,       !|
     *                  DisorNorm, DisorNorm8_gra, iunim,    !|
     *                  ixway, wave, energy,                 !|
     *                  nReflex1, indices, indrange,         !|
     *                  filter, ifail)                       !|
c+============================================================+
        if (ifail.gt.0)  Then  !-------------------+
          if (ifail.ge.1   .AND.                  !|
     *        ifail.le.12)  Then  !-------------+  |
            txt(1) = Errors(ifail)             !|  |
          else   !------------------------------+  |
            txt(1) = 'Unknown error'           !|  |
          endif  !------------------------------+  |
          ifail = 1                               !|
          goto 999 !-------------------------------+-------------+
        elseif (ifail.lt.0)  Then  !---------------+             v
          goto 999 !-------------------------------+-------------+
        endif  !-----------------------------------+             v
        rho          = 0.
        isyngony     = 0
        icall        = 0                                !first x0h call
        jpr          = 0
c       jpr          = 2
        qb           = 0.
        xr0          = 0.
        xi0          = 0.
        xrh          = 0.
        xih          = 0.
        xqh          = 0.
        xrf          = 0.
        xif          = 0.
        xqf          = 0.
        inmax        = 3
c Calculate crystal constants:
        if (ixway.gt.3) Then !--+
          w = 1.               !|
        else !------------------+
c This will help to determine   |
c wave via x-ray line:          |
          w = wave             !|
        endif !-----------------+
        do i=1,6  !==================+
          Lattice_constants(i) = 0. !|
        enddo !======================+
        Call X0h1 (Lattice_constants, w,
     *             code, rho, 0,
     *             0, 0, 0,
     *             isyngony, icall,
     *             jpr, qb, xr0, xi0,
     *             xrh, xih, xqh,
     *             xrf, xif, xqf, ifail)
        if (ifail.ne.0)  goto 999 !------------------------------+
c                                                                v
c wave was calculated in x0h by x-ray line:
        if (ixway.eq.3) wave = w

c Build INTERNAL surface normal
c (see the sub below):
        Call    SurfNormVect (nBaseNormal,
     *                        nMiscutDirection,
     *                        DisorNorm8_gra,
     *                        Surface_Normal,ifail)
        if (ifail.ne.0)  goto 999 !------------------------------+
c                                                                v
c Build the incident vector:
c (see incivect.for)
        wave8 = dble(wave)
        Call IncidentVector (wave8, inmax, nReflex1,
     *                       indices, h_vec, h_mod,
     *                       h1xh2, h1xh2_mod,
     *                       k0_copl, wave_copl, Uncopl,
     *                       k0(1,1), k0_mod,
     *                       ifail)
        if (ifail.ne.0) goto 999 !-------------------------------+
c                                                                v
c ================================================================ Begin added 2014.02.25
c Build the second k0:
        Call DVecCon (h1xh2, -1.D0, h1xh2_inv, ii)
        Call DVecSum2 (k0_copl, 1.D0,
     *                 h1xh2_inv,  Uncopl,
     *                 k0(1,2), ii)

c Flags for problematic situations:
        normal_inverted = 0
        h1xh2_inverted  = 0
        found_two_k0    = 0

c Convolution of INTERNAL surface normal with k0(1) and k0(2):
        s8(1) = DVecSca2 (k0(1,1), Surface_Normal, ii)
        s8(2) = DVecSca2 (k0(1,2), Surface_Normal, ii)

        if (s8(1).lt.0. .AND. s8(2).lt.0.) Then !-----------+
c None of the two k0 can enter the crystal. Swap normal     |
          normal_inverted = 1                              !|
          do  i=1,ii   !================================+   |
            nBaseNormal(i)      = -nBaseNormal(i)      !|   |
            nMiscutDirection(i) = -nMiscutDirection(i) !|   |
            Surface_Normal(i)   =- Surface_Normal(i)   !|   |
          enddo   !=====================================+   |
          s8(1) = -s8(1)                                   !|
          s8(2) = -s8(2)                                   !|
        endif !---------------------------------------------+

        if     (s8(1).gt.0. .AND. s8(2).gt.0.) Then !------------------+
c We have two solutions for k0. User will have to choose:              |
          found_two_k0 = 1                                            !|
          do n=1,2 !================================================+  |
c Angles of the two k0 candidates with the surface:                 |  |
            angle8(n) = 90. - DAngle2(Surface_Normal, k0(1,n), ii) !|  |
c Unit vectors along SIGMA-polarization direction:                  |  |
            if (angle8(n).gt.5. .AND. angle8(n).lt.175.) Then !-+   |  |
c Normally we use: sg=[k*surf_norm]                             |   |  |
              Call DVecVec (k0(1,n), Surface_Normal,           !|   |  |
     *                                    Sg_unit(1,n), ii)    !|   |  |
            else !----------------------------------------------+   |  |
c When k || surf_norm, we cannot use  sg=[k*surf_norm]          |   |  |
c and then we choose                  sg=[k*[h1*h2]]            |   |  |
              Call DVecVec (k0(1,n), h1xh2, Sg_unit(1,n), ii)  !|   |  |
            endif !---------------------------------------------+   |  |
c Unit vectors along PI-polarization dir-on: pi=[sg*k]:             |  |
            Call DVecVec (Sg_unit(1,n), k0(1,n), Pi_unit(1,n), ii) !|  |
c Normalize the vectors:                                            |  |
            Call DUnitVec2 (Sg_unit(1,n), Sg_unit(1,n),ii)         !|  |
            Call DUnitVec2 (Pi_unit(1,n), Pi_unit(1,n), ii)        !|  |
          enddo !===================================================+  |
                                                                      !|
        elseif (s8(1).lt.0. .AND. s8(2).gt.0.) Then !------------------+
c Only second k0 can enter the crystal.                                |
c Swap h1 and h2 to make it the main k0                                |
          h1xh2_inverted  = 1                                         !|
          Call DVecCopy_dd(h1xh2_inv, h1xh2, ii)                      !|
          Call DVecCopy_dd(k0(1,2), k0(1,1), ii)                      !|
          do  i=1,ii   !===================+                           |
            j            = indices(i,1)   !|                           |
            indices(i,1) = indices(i,2)   !|                           |
            indices(i,2) = j              !|                           |
          enddo   !========================+                           |
        endif !--------------------------------------------------------+
c ================================================================ Finish added 2014.02.25

        if (ixway.gt.3) Then !-----+
c ixway=4 & ixway=5 the wavelength |
c has just been determined in      |
c IncidentVector:                  |
          wave = Sngl(wave8)      !|
        endif !--------------------+
        energy = wave2energy(wave)

        lcod = Max(Len_Trim(code),1)
        lrad = Max(Len_Trim(radiat),1)
        Call Pakint(nBaseNormal,ii,base_txt,lbase)
        Call Pakint(nMiscutDirection,ii,misc_txt,lmisc)
        Call PakReal(DisorNorm,1,disor_txt,ldisor)

        if (ixway.ne.5) then !--+
          nref = 3             !|
        else !------------------+
          nref = 4             !|
        endif !-----------------+
        do i=1,nref-1  !=======================================+
          Call Pakint(indices(1,i+1),ii,refl_txt(i),lrefl(i)) !|
        enddo !================================================+
        write(lunout,91,err=11)
        write(lundup,92,err=11)
  91    format(' <blockquote>'/
     *  ' <b>Questions? Please, read</b>',1x,
     *      '<a href="/BRL.html"><b><i>BRL</i> documentation</b></a>.')
  92    format(' <blockquote>')

        if (normal_inverted.ne.0 .OR. h1xh2_inverted.ne.0) Then !---+
          do j=lunout,lundup  !=============================+       |
            write(j,16,err=11)                             !|       |
            if (normal_inverted.ne.0) write(j,44,err=11)   !|       |
            if (h1xh2_inverted.ne.0)  write(j,45,err=11)   !|       |
            write(j,42,err=11)                             !|       |
          enddo !===========================================+       |
        endif !-----------------------------------------------------+
  44    format(
     *  ' <p><b><font color="red">WARNING:</font></b> Crystal surface',
     *  1x,'normal was inverted for the resulted wave to be able',
     *  1x,'entering the crystal.</p>')
  45    format(
     *  ' <p><b><font color="red">WARNING:</font></b> The order of h1',
     *  1x,'and h2 was inverted for the resulted wave to be able',
     *  1x,'entering the crystal.</p>')

        ixw2 = ixway
        if (ixw2 .eq. 2) ixw2 = 1       !we always specify wave instead of energy
c The following two lines are commented out
c because we want the information about how
c the wavelength was calculated to propagate
c into brls12:
c       if (ixw2 .eq. 4) ixw2 = 1       !wave was determined from coplanar case
c       if (ixw2 .eq. 5) ixw2 = 1       !wave was determined from reflection-3

c This is needed to prevent IP as hostname in the beamline
c logs when the hostname is not given by the web server:
        do j=lunout,lundup  !=========================================+
          write(j,9,err=11)     hostname(1:lhost),                   !|
     *                          address(1:laddr),                    !|
     *                          txcount,                             !|
     *                          code(1:lcod),                        !|
     *                          iHenkeCowan,                         !|
     *                          base_txt(1:lbase),                   !|
     *                          misc_txt(1:lmisc),                   !|
     *                          disor_txt(1:ldisor), iunim,          !|
     *                          Wave8, radiat(1:lrad), ixw2,         !|
c    *                          indrange, filter,                    !|
     *                          (refl_txt(i)(1:lrefl(i)),i=1,nref-1) !|
          write(j,16,err=11)                                         !|
        enddo  !======================================================+
  9     format(' <p>'/
     *  ' <form',1x,
c    *    'action="https://x-server.gmca.aps.anl.gov/cgi/brl_frm2.pl"',
     *    'action="/cgi/brl_frm2.pl"',
     *   1x,'method="get" name="brlfrm2"',
     *   1x,'onsubmit="return brlfrm2_validate()">'/
     *  ' <input type="hidden" name="host"     value="',a,'">'/
     *  ' <input type="hidden" name="address"  value="',a,'">'/
c    *  ' <input type="hidden" name="referer"  value="',a,'">'/
     *  ' <input type="hidden" name="jobid"    value="',a,'">'/
     *  ' <input type="hidden" name="code"     value="',a,'">'/
     *  ' <input type="hidden" name="df1df2"   value="',i2,'">'/
     *  ' <input type="hidden" name="basenorm" value="',a,'">'/
     *  ' <input type="hidden" name="miscaxis" value="',a,'">'/
     *  ' <input type="hidden" name="miscut"   value="',a,'">'/
     *  ' <input type="hidden" name="unim"     value="',i1,'">'/
     *  ' <input type="hidden" name="wave"     value="',g15.8,'">'/
     *  ' <input type="hidden" name="line"     value="',a,'">'/
     *  ' <input type="hidden" name="xway"     value="',i1,'">'/
c    *  ' <input type="hidden" name="indrange" value="',i2,'">'/
c    *  ' <input type="hidden" name="filter"   value="',g12.5,'">'/
     *  ' <input type="hidden" name="rfx1"     value="',a,'">'/
     *  ' <input type="hidden" name="rfx2"     value="',a,'">'/:
     *  ' <input type="hidden" name="rfx3"     value="',a,'">')
16      format(
     *  ' <table cellspacing="0" cellpadding="0" border="0"',1x,
     *                                    'bgcolor="#c1c1c1"><tr><td>'/
     *  ' <table border="1" bgcolor="#c1c1c1" bordercolor="#800000"',1x,
     *                                      'cellpadding="5"><tr><td>')

        l = Len_Trim(Symmetry(isyngony+1))
        m = Len_Trim(DBtext(iHenkeCowan+2))
        do j=lunout,lundup  !==========================================+
          write(j,6,err=11)  nobr,                                    !|
     *                       code(1:lcod), Symmetry(isyngony+1)(1:l), !|
     *                       DBtext(iHenkeCowan+2)(1:m),              !|
     *                       nobr,                                    !|
     *                       base_txt(1:lbase), misc_txt(1:lmisc),    !|
     *                       disor_txt(1:ldisor), Uni(iunim+1),       !|
     *                       nobr,nobr,                               !|
     *                       wave, energy, radiat(1:lrad),            !|
     *                       nobr                                     !|
          write(j,66,err=11) nobr,nobr,nobr,nobr                      !|
        enddo  !=======================================================+
        radiat = ' '
  6     format(
     *  ' <p>',a,'<b>Crystal: </b><font color="red">',a,'</font>',1x,
     *    '&nbsp; Symmetry: <font color="red">',a,'</font>',1x,
     *    '&nbsp; X0h data: <font color="red">',a,'</font>',1x,
     *    '</span><br>'/
     *  ' ',a,'<b>Surface: </b>',1x,
     *    '&nbsp; base plane=(<font color="red">',a,'</font>)',1x,
     *    '&nbsp; miscut direction=(<font color="red">',a,'</font>)',1x,
     *    '&nbsp; miscut angle=<font color="red">',a,'</font> ',a,1x,
     *    '</span><br>'/
     *  ' ',a,'<b>Thickness </b>(microns): ',1x,
     *    '<input type="text" name="thickness" size="5" value="100.">',
     *    1x,'</span></p>'/
     /  /
     *  ' <p>',a,'<b>X-rays: </b>',1x,
     *    '&nbsp; wavelength=<font color="red">',g12.5,' A</font>',1x,
     *    '&nbsp; energy=<font color="red">',g12.5,' keV</font>',1x,
     *    '&nbsp; line=<font color="red">',a,'</font></span><br>'/
     *  ' ',a,'<b>X-ray polarization:</b> &nbsp; '/
     *    '    <select name="ipol">'/
     *    '      <option value="0" selected>[1] Mixed (unpolarized)'/
     *    '      <option value="1">[2] Linear, given by angle to pi0'/
     *    '      <option value="2">[3] Linear, along [k0*k1]'/
     *    '    </select>'/
     *    ' &nbsp; Angle to pi0 for mode [2]:',1x,
     *       '<input type="text" name="polang" size="3" value="0.">',1x,
     *    ' </span></p>'/)
  66    format(
     *  ' <p>',a,'<b>Scan limits (Theta1):</b> &nbsp; from',1x,
     *    '<input type="text" name="q1min" size="5" value="-50.">',
     *    1x,'&nbsp; to',1x,
     *    '<input type="text" name="q1max" size="5" value="50.">',
     *    1x,'&nbsp; points =',1x,
     *    '<input type="text" name="nq1" size="3" value="101">',
     *    1x,'</span><br>'/
     *  ' ',a,'<b>Scan limits (Theta2):</b> &nbsp; from',1x,
     *    '<input type="text" name="q2min" size="5" value="0.">',
     *    1x,'&nbsp; to',1x,
     *    '<input type="text" name="q2max" size="5" value="0.">',
     *    1x,'&nbsp; points =',1x,
     *    '<input type="text" name="nq2" size="3" value="1">',
     *    1x,'</span><br>'/
     *  ' ',a,'<b>Scan axes:</b> &nbsp; '/
     *    '   <select name="scan">'/
     *    '     <option value="0" selected>[1] Theta2 along sigma0'/
     *    '     <option value="1">[2] Theta2 along [k0*surface_normal]'/
     *    '     <option value="2">[3] Theta2 along [k0*k1]'/
     *    '   </select>'/
     *    ' &nbsp; <b>Scan angle units:</b> &nbsp;'/
     *    '   <select name="units">'/
     *    '     <option value="0">degr.'/
     *    '     <option value="1">min.'/
     *    '     <option value="2">mrad.'/
     *    '     <option value="3" selected>sec.'/
     *    '     <option value="4">urad.'/
     *    '   </select>'/
     *    ' </span></p>'/
     /  /
     *  ' <p>',a,'<b>Specified reflections:</b></span>')

        icall = 3       !(md_NewWave): not first call, same crystal as before, but different X-ray wavelength
        do iref=2,nref  !===============================================+
          Call X0h1 (Lattice_constants, wave,                          !|
     *               code, rho, 1,                                     !|
     *               indices(1,iref),                                  !|
     *               indices(2,iref),                                  !|
     *               indices(3,iref),                                  !|
     *               isyngony, icall,                                  !|
     *               jpr, qb, xr0, xi0,                                !|
     *               xrh, xih, xqh,                                    !|
     *               xrf, xif, xqf, ifail)                             !|
          icall = 5  !(md_NewHKL): not the first call; different hkl    |
c Calculate deviation from the Bragg condition                          |
c expresees in |x0|:                                                    |
          c8 = 2.*DVecSca2(k0(1,1),h_vec(1,iref),ii) + h_mod(iref)**2  !|
          c8 = c8 / (xr0 * k0_mod**2)                                  !|
          c8 = Abs(c8)                                                 !|
                                                                       !|
c         do i=1,ii  !==================+                               |
c           iii(i) = indices(i,iref)   !|                               |
c         enddo  !======================+                               |
c         if (isyngony.eq.2)  Then  !---------------------+             |
c           iii(3) = -(indices(1,iref)+indices(2,iref))  !|             |
c           iii(4) = indices(3,iref)                     !|             |
c           inmax  = 4                                   !|             |
c         endif  !----------------------------------------+             |
c         Call PakInt (iii,inmax,str,l)                                !|
                                                                       !|
          w = 100.*abs(xrh/xr0)                                        !|
          do n=lunout,lundup !=====================================+    |
            write(n,7,err=11) nobr,                               !|    |
     *                        iref-1,                             !|    |
     *                        refl_txt(iref-1)(1:lrefl(iref-1)),  !|    |
     *                        qb, w, c8                           !|    |
          enddo  !=================================================+    |
        enddo  !========================================================+
  7     format(' <br>',a/
     *  ' <b>Reflection',i1,'</b> = (<font color="red">',a,'</font>)',
     *     1x,' &nbsp; &nbsp; QB = ',f6.3,' degr.',
     *     1x,' &nbsp; &nbsp; |xh/x0| = ',f6.3,'%',
     *     1x,' &nbsp; &nbsp; |alpha/x0| = ',g10.3,'</span>')

        do n=lunout,lundup !===============+
          write(n,8,err=11) nobr,         !|
     *                      12-nref,      !|
     *                      nobr,         !|
     *                      -indrange,    !|
     *                      -indrange,    !|
     *                      -indrange,    !|
     *                       indrange,    !|
     *                       indrange,    !|
     *                       indrange,    !|
     *                       filter       !|
        enddo  !===========================+
  8     format(/' </p><p>',a/
     *  ' <b>Additional reflections search results',1x,
     *     ' (you can select up to ',i2,' planes if available):</b>',1x,
     *     '</span></p>'/
     *  ' <p>',a,'&nbsp; &nbsp; Searching',1x,
     *     'from (<font color="red">',3i4,'</font>)',1x,
     *     'to (<font color="red">',3i4,'</font>)',1x,
     *     '<b>&nbsp; &nbsp; </b>Intensity filter |xh/x0| &gt; ',1x,
     *     '<font color="red">',f6.3,'%</font> ...</span></p>'/)

c Search for additional reflections:
        ifound = 0
        do i=-indrange,indrange !==================================+
        do j=-indrange,indrange !===============================+  |
        do k=-indrange,indrange !============================+  |  |
                                                            !|  |  |
          do n=1,nref !========================+             |  |  |
c If the reflection already specified          |             |  |  |
c (including 0,0,0):                           |             |  |  |
            if (i .eq. indices(1,n) .AND.     !|             |  |  |
     *          j .eq. indices(2,n) .AND.     !|             |  |  |
     *          k .eq. indices(3,n)) goto 10 !-+----------+  |  |  |
          enddo !==============================+          |  |  |  |
c Calculate reflection strength:                          |  |  |  |
          Call X0h1 (Lattice_constants, wave,            !|  |  |  |
     *               code, rho, 1,                       !|  |  |  |
     *               i, j, k,                            !|  |  |  |
     *               isyngony, icall,                    !|  |  |  |
     *               jpr, qb, xr0, xi0,                  !v  |  |  |
     *               xrh, xih, xqh,                      !|  |  |  |
     *               xrf, xif, xqf, ifail)               !|  |  |  |
c This can only be something like "no Bragg angle":       |  |  |  |
          if (ifail .ne. 0) goto 10  !--------->----------+  |  |  |
                                                         !|  |  |  |
          w = 100.*abs(xrh/xr0)                          !|  |  |  |
          if (w .lt. filter) goto 10  !-------->----------+  |  |  |
                                                         !|  |  |  |
          iii(1) = i                                     !|  |  |  |
          iii(2) = j                                     !|  |  |  |
          iii(3) = k                                     !|  |  |  |
c Build h-vector:                                        !|  |  |  |
          Call Build_h(iii,h_v,h_m,ii)                   !|  |  |  |
          c8 = 2.*DVecSca2(k0(1,1),h_v,ii) + h_m**2      !|  |  |  |
          c8 = c8 / (xr0 * k0_mod**2)                    !|  |  |  |
          c8 = Abs(c8)                                   !|  |  |  |
c if alpha less than alpha_range*|x0|:                   !v  |  |  |
          if (c8 .lt. alpha_range) Then  !-------------+  |  |  |  |
c Output additional reflection on the HTML form:       |  |  |  |  |
            ifound = ifound+1                         !|  |  |  |  |
            if (ifound.eq.1) Then !-------+            |  |  |  |  |
              do n=lunout,lundup !====+   |            |  |  |  |  |
                write(n,30,err=11)   !|   |            |  |  |  |  |
              enddo !=================+   |            |  |  |  |  |
            endif  !----------------------+            |  |  |  |  |
c           if (isyngony.eq.2) Then !--+               |  |  |  |  |
c             iii(3) = -(i+j)         !|               |  |  |  |  |
c             iii(4) = k              !|               |  |  |  |  |
c           endif  !-------------------+               |  |  |  |  |
c           Call PakInt (iii,inmax,str,l)             !|  |  |  |  |
            iwrk(1) = nref+ifound-1                   !|  |  |  |  |
            Call PakInt (iwrk,1,wrk,m)                !|  |  |  |  |
            Call PakInt (iii,ii,inputstring,l)        !|  |  |  |  |
            do n=lunout,lundup !====================+  |  |  |  |  |
              write(n,29,err=11) wrk(1:m),         !|  |  |  |  |  |
     *                           inputstring(1:l), !|  |  |  |  |  |
     *                           wrk(1:m),         !|  |  |  |  |  |
     *                           inputstring(1:l), !|  |  |  |  |  |
     *                           qb, w, c8         !|  |  |  |  |  |
            enddo !=================================+  |  |  |  |  |
          endif  !-------------------------------------+  |  |  |  |
                                                         !|  |  |  |
  10      continue  !<------------------------------------+  |  |  |
        enddo  !=============================================+  |  |
        enddo  !================================================+  |
        enddo  !===================================================+
        if (ifound.gt.0) Then !-------+
          do n=lunout,lundup !====+   |
            write(n,31,err=11)   !|   |
          enddo !=================+   |
        endif  !----------------------+
  30    format(' <table border="0" cellpadding="1" cellspacing="0">')
  29    format(
     *  ' <tr valign="middle" style="white-space:nowrap;">'/
     *  '   <td><input type="checkbox" name="rfx',a,
     *                                      '" value="',a,'"></td>'/
     *  '   <td><b>Reflection',a,'</b> = (',a,') </td>'/
     *  '   <td> QB = ',f6.3,' degr. </td>'/
     *  '   <td> |xh/x0| = ',f6.3,'% </td>'/
     *  '   <td> |alpha/x0| = ',g10.3,' </td>'/
     *  ' </tr>')
  31    format(' </table>')

        if (ifound. eq. 0) Then  !--------+
          do n=lunout,lundup !=========+  |
            write(n,39,err=11) nobr   !|  |
          enddo  !=====================+  |
        endif  !--------------------------+
  39    format(' <p>',a/
     *  ' <b>*** No additional Bragg planes found. *** </b></span></p>')

        if (found_two_k0.ne.0) Then !---------------------------+
          do n=lunout,lundup !===============================+  |
            write(n,30,err=11)                              !|  |
            write(n,47,err=11)                              !|  |
            write(n,41,err=11) 1, 'checked', 1,             !|  |
     *                         (k0(i,1),i=1,ii), angle8(1), !|  |
     *                         (Sg_unit(i,1),i=1,ii),       !|  |
     *                         (Pi_unit(i,1),i=1,ii)        !|  |
            write(n,41,err=11) 2, ' ', 2,                   !|  |
     *                         (k0(i,2),i=1,ii), angle8(2), !|  |
     *                         (Sg_unit(i,2),i=1,ii),       !|  |
     *                         (Pi_unit(i,2),i=1,ii)        !|  |
            write(n,31,err=11)                              !|  |
          enddo  !===========================================+  |
        else !--------------------------------------------------+
          do n=lunout,lundup !=======+                          |
            write(n,46,err=11)      !|                          |
          enddo  !===================+                          |
        endif !-------------------------------------------------+
  47    format(' <tr valign="middle" style="white-space:nowrap;">'/
     *  ' <td colspan="4">&nbsp;<br>'/
     *  ' <b>Found two K0 satisfying the conditions. Choose one:',1x,
     *                                                 '</b></td></tr>')
  41    format(' <tr valign="middle" style="white-space:nowrap;">'/
     *  ' <td>&nbsp;<input type="radio" name="kchoice" value="',i1,'" ',
     *                        a,'> k0_',i1,'=(',3g10.3,') &nbsp; </td>'/
     *  '  <td>angle_with_surface:',f7.2,' degr, &nbsp; </td>',1x,
     *  '  <td>sg=(',3g10.3,'), &nbsp; </td>',1x,
     *  '  <td>pi=(',3g10.3,') </td></tr>')
  46    format(' <input type="hidden" name="kchoice" value="0">')

        iwatch = 0
        iwrk(1) = ifound + nref - 1     ! how many planes we can expect
        Call Pakint(iwrk,1,wrk,i)
        do n=lunout,lundup !===========================+
            write(n,40,err=11) nobr, wrk(1:i), iwatch !|
            write(n,42,err=11)                        !|
            write(n,43,err=11)                        !|
        enddo  !=======================================+
  40    format(' <p>',a/
     *    ' <input type="hidden" name="maxrfx" value="',a,'">'/
     *    ' <input type="hidden" name="watch"  value="',i1,'">'/
     *    ' <input type="submit" value="Submit">',1x,
     *    ' <input type="reset"> </span></p>')
  42    format(' </td></tr></table>'/
     *         ' </td></tr></table>')
  43    format(' </form></p>'/
     *         ' </blockquote>'/)

  555   continue
        Call Navigate_Site(lunout,referer,i)
        if (i.ne.0) goto 11

        write (lunout,99,err=11)
        write (lundup,99,err=11)
  99    format(' </body></html>'/)

  11    continue
        Call exit_quiet()
c==================================================================
c                            E R R O R S
c==================================================================
  999   continue
        write(lunout,19,err=11)
        write(lundup,19,err=11)
        ifail = iabs(ifail)
        Call HTMLstop       (lunout,ifail)
        Call HTMLstop_noJPG (lundup,ifail)
        write (lunout,98,err=11)
        write (lundup,98,err=11)
  19    format(
     *  ' <center>'/
     *  ' <table border bgcolor="#ffffff"><tr><td>'/
     *  ' <table border="0" bgcolor="#ffffff">')
  98    format(
     *  ' </table>'/
     *  ' </td></tr></table>'/
     *  ' </center>')
        goto 555

        end

c=================================================================

        Subroutine keys2brl (nkeys, keys, vals,
     *                       code,
     *                       nBaseNormal, nMiscutDirection,
     *                       DisorNorm, DisorNorm8_gra, iunim,
     *                       ixway, wave, energy,
     *                       nReflex1, indices, indrange,
     *                       filter, ifail)
        Real*8          DisorNorm8_gra
        Real*4          DisorNorm, wave, energy,
     *                  UnitData(5), pi, filter,
     *                  AngDegree, AngMinute,
     *                  AngSec, AngMrad,
     *                  AngMurad
        Integer         nkeys, ixway, iunim,
     *                  nReflex1, indices(3,nReflex1), indrange,
     *                  nBaseNormal(3), nMiscutDirection(3),
     *                  i, m, n, o, omax, ifail,
     *                  ii/3/                           !Number of Bragg planes indices (3 or 4)
        Integer         iend
        External        iend
        Character       keys(nkeys)*(*), vals(nkeys)*(*),
     *                  code*(*), valbuf*32, key*4
c -------------------------------------------------------
        Character radiat*6
        Common  /x0pa1/ radiat
c -------------------------------------------------------
c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer             iHenkeCowan
        Common /HenkeCowan/ iHenkeCowan
c -------------------------------------------------------
        code        = ' '
        radiat      = ' '
        wave        = 0.
        energy      = 0.
        ixway       = 0
        indrange    = 0
        iHenkeCowan = -1
        do m=1,ii !=================+
          do  o=1,4 !============+  |
            indices(m,o) = 0    !|  |
          enddo  !===============+  |
          nBaseNormal(m)      = 0  !|
          nMiscutDirection(m) = 0  !|
        enddo  !====================+
        DisorNorm   = 0.
        iunim       = 0

        pi          = 4.*atan(1.)               !the "pi" number
        AngDegree   = 2.*pi/360.
        AngMinute   = AngDegree/60.
        AngSec      = AngDegree/3600.
        AngMrad     = 1.E-03
        AngMurad    = 1.E-06

        UnitData(1) = AngDegree
        UnitData(2) = AngMinute
        UnitData(3) = AngMrad
        UnitData(4) = AngSec
        UnitData(5) = AngMurad

        ifail       = 0

c Query String:
c code=Silicon xway=1 wave=1.54 line=Cu-Ka1 i11=7 i12=8 i13=9 i21=7& i22=8 i23=9 i31=7 i32=8 i33=9 indrange=5

c Take crystal code selected from a list:
        Call getbykey(nkeys,keys,vals,'code',code)
        if (Len_Trim(code).eq.0) Then !+
c No crystal name specified            |
          ifail = 1                   !|
          goto 28                     !|
        endif  !-----------------------+

c Base surface normal:
        key = 'n_'
        n   = Len_Trim(key)
        do  m=1,ii !======================================+
          write(key(n:n),'(i1)') m                       !|
          Call getbykey(nkeys,keys,vals,key(1:n),valbuf) !|
          if (Len_Trim(valbuf) .gt. 0) Then  !-------+    |
            Call Rdint(nBaseNormal(m),1,valbuf,i)   !|    |
            if (i.ne.0) Then  !---+                  |    |
                ifail = 8        !|                  |    |
                goto 28          !|                  |    |
            endif  !--------------+                  |    |
          else  !------------------------------------+    |
            ifail = 8                               !|    |
            goto 28                                 !|    |
          endif  !-----------------------------------+    |
        enddo  !==========================================+
        if (nBaseNormal(1).eq.0 .AND.
     *      nBaseNormal(2).eq.0 .AND.
     *      nBaseNormal(3).eq.0) Then !--+
c Wrong or missing normal indices (hkl)  |
            ifail = 8                   !|
            goto 28                     !|
        endif !--------------------------+

c Miscut angle:
        Call getbykey(nkeys,keys,vals,'miscut',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !---+
          Call Rdreal(DisorNorm,1,valbuf,i)  !|
          if (i.ne.0) Then !-+                |
            ifail = 9       !|                |
            goto 28         !|                |
          endif  !-----------+                |
        endif  !------------------------------+

c Miscut angle units:
        Call getbykey(nkeys,keys,vals,'unim',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !---+
          Call Rdint(iunim,1,valbuf,i)       !|
          if (i.ne.0) Then !-+                |
            ifail = 10      !|                |
            goto 28         !|                |
          endif  !-----------+                |
        endif  !------------------------------+

        if (iunim .lt. 0 .OR.
     *      iunim .gt. 4 ) Then !+
c Unexpected type of units       |
          ifail = 10            !|
          goto 28               !|
        endif  !-----------------+

        DisorNorm8_gra = Dble(DisorNorm) * UnitData(iunim+1) / AngDegree
        if (DisorNorm8_gra .lt. -90. .OR.                       !was +-360
     *      DisorNorm8_gra .gt. 90. ) Then !--+
c Wrong miscut angle:                         |
          ifail = 9                          !|
          goto 28                            !|
        endif  !------------------------------+

c Miscut direction:
        key = 'm_'
        n   = Len_Trim(key)
        do  m=1,ii !======================================+
          write(key(n:n),'(i1)') m                       !|
          Call getbykey(nkeys,keys,vals,key(1:n),valbuf) !|
          if (Len_Trim(valbuf) .gt. 0) Then  !----------+ |
            Call Rdint(nMiscutDirection(m),1,valbuf,i) !| |
            if (i.ne.0) Then !--+                       | |
              ifail = 11       !|                       | |
              goto 28          !|                       | |
            endif  !------------+                       | |
          endif  !--------------------------------------+ |
        enddo  !==========================================+
        if (Abs(DisorNorm).gt.1.E-10  .AND.             !non-zero
     *      nMiscutDirection(1).eq.0  .AND.
     *      nMiscutDirection(2).eq.0  .AND.
     *      nMiscutDirection(3).eq.0)  Then !-------+
c Wrong or missing miscut direction indices (hkl)   |
            ifail = 11                             !|
            goto 28                                !|
        endif !-------------------------------------+

        Call getbykey(nkeys,keys,vals,'xway',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !--+
          Call Rdint(ixway,1,valbuf,i)      !|
          if (i.ne.0) Then !-+               |
            ifail = 2       !|               |
            goto 28         !|               |
          endif  !-----------+               |
        endif  !-----------------------------+

        if (ixway .lt. 1 .OR.
     *      ixway .gt. 5 ) Then !+
c Unexpected type of x-ray input |
          ifail = 2             !|
          goto 28               !|
        endif  !-----------------+

        if (ixway.eq.1 .OR. ixway.eq.2) Then  !----------+wave or energy
                                                        !|
          radiat = 'none'                               !|
          Call getbykey(nkeys,keys,vals,'wave',valbuf)  !|
          if (Len_Trim(valbuf) .gt. 0) Then !--+         |
            Call Rdreal (wave,1,valbuf,i)     !|         |
            if (i.ne.0) Then !---+             |         |
              ifail = 3         !|             |         |
              goto 28           !|             |         |
            endif  !-------------+             |         |
          endif  !-----------------------------+         |
          if (wave.le.0.) Then !-+                       |
c Wrong or missing x-ray         |                       |
c wavelength/energy:             |                       |
            ifail = 3           !|                       |
            goto 28             !|                       |
          endif  !---------------+                       |
                                                        !|
          if (ixway.eq.1) Then !-----+                   |
            energy = 12.3981/wave   !|                   |
          else  !--(ixway.eq.2)------+                   |
            energy = wave           !|                   |
            wave   = 12.3981/energy !|                   |
          endif  !-------------------+                   |
                                                        !|
        elseif (ixway.eq.3) Then !-----------------------+x-ray line
                                                        !|
          Call getbykey(nkeys,keys,vals,'line',radiat)  !|
          if (Len_Trim(radiat).eq.0 .OR.                !|
     *                  radiat.eq.'none') Then !--+      |
c missing keyword for radiation line              |      |
            ifail = 4                            !|      |
            goto 28                              !|      |
          endif  !--------------------------------+      |
                                                        !|
        elseif (ixway.eq.4) Then !-----------------------+coplanar case
                                                        !|
          wave   = -1.                                  !|
          radiat = 'none'                               !|
                                                        !|
        elseif (ixway.eq.5) Then !-----------------------+reflection-3
                                                        !|
          wave   = -2.                                  !|
          radiat = 'none'                               !|
                                                        !|
        endif  !-----------------------------------------+

        if (ixway.ne.5) Then !----+
          omax = 2               !|
c Reflection-3 is not involved,   |
c but we need a workaround since  |
c "IncidentVector" checks that    |
c all h-vectors must be non-zero !|
          indices(1,4) = 1       !|
          indices(2,4) = 1       !|
          indices(3,4) = 1       !|
        else !--------------------+
          omax = 3               !|
        endif  !------------------+
        key = 'i__'
        n   = Len_Trim(key)
        do  o=1,omax !==========================================+
          do  m=1,ii !=======================================+  |
            write(key(n-1:n),'(2i1)') o,m                   !|  |
            Call getbykey(nkeys,keys,vals,key(1:n),valbuf)  !|  |
            if (Len_Trim(valbuf) .gt. 0) Then !-------+      |  |
              Call Rdint(indices(m,o+1),1,valbuf,i)  !|      |  |
              if (i.ne.0) Then  !-+                   |      |  |
                ifail = 5        !|                   |      |  |
                goto 28          !|                   |      |  |
              endif  !------------+                   |      |  |
            else  !-----------------------------------+      |  |
              ifail = 5         !(field not filled)  !|      |  |
              goto 28                                !|      |  |
            endif  !----------------------------------+      |  |
          enddo  !===========================================+  |
          if (indices(1,o+1).eq.0 .AND.                        !|
     *        indices(2,o+1).eq.0 .AND.                        !|
     *        indices(3,o+1).eq.0) Then !--+                    |
c Wrong or missing reflex indices (hkl)    |                    |
            ifail = 5                     !|                    |
            goto 28                       !|                    |
          endif !--------------------------+                    |
        enddo  !================================================+

        Call getbykey(nkeys,keys,vals,'indrange',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !--+
          Call Rdint(indrange,1,valbuf,i)   !|
          if (i.ne.0) Then !-+               |
            ifail = 6       !|               |
            goto 28         !|               |
          endif  !-----------+               |
        endif  !-----------------------------+

        if (indrange .lt. 3 .OR.
     *      indrange .gt. 15 ) Then !+
c Wrong index search range           |
          ifail = 6                 !|
          goto 28                   !|
        endif  !---------------------+

        Call getbykey(nkeys,keys,vals,'filter',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !--+
          Call Rdreal(filter,1,valbuf,i)    !|
          if (i.ne.0) Then !-+               |
            ifail = 12      !|               |
            goto 28         !|               |
          endif  !-----------+               |
        endif  !-----------------------------+

        if (filter .lt. 0. .OR.
     *      filter .gt. 100. ) Then !+
c Wrong index search range           |
          ifail = 12                !|
          goto 28                   !|
        endif  !---------------------+

        Call getbykey(nkeys,keys,vals,'df1df2',valbuf)
        if (Len_Trim(valbuf) .gt. 0) Then !---+
          Call Rdint(iHenkeCowan,1,valbuf,i) !|
          if (i.ne.0) Then  !-+               |
            ifail = 7        !|               |
            goto 28          !|               |
          endif  !------------+               |
        else  !-------------------------------+
          ifail = 7        !(field not filled)|
          goto 28                            !|
        endif  !------------------------------+
c -1=Auto,0=X0h,
c  1=Henke(f1),   2=Henke(f1,f2)
c  3=Cowan(f1),   4=Cowan(f1,f2)
c  5=Windt(f1),   6=Cowan(f1,f2)
c  7=Chantler(f1),8=Chantler(f1,f2)
        if (iHenkeCowan.lt.-1 .OR.
     *      iHenkeCowan.gt.8) Then !-+
          ifail = 7                 !|
          goto 28                   !|
        endif  !---------------------+

  28    continue
        return
        end

c======================================================================

        Subroutine      SurfNormVect (nBaseNormal,      !hkl of base normal (input)
     *                                nMiscutDirection, !hkl of miscut reper (input)
     *                                Disorient8,       !miscut angle in Units1 (input)
     *                                Surface_Normal,   !hkl of resulted surface normal (output)
     *                                ifail)            !failure flag (output)
c +-----------------------------------------------------------+
c |  ATTENTION!!! This subroutine does not want to work with  |
c |               NDP FORTRAN if included into a library!     |
c +-----------------------------------------------------------+

c +-----------------------------------------------------------+
c |  ATTENTION!!! This subroutine differs from the one in     |
c |               SurfNorm.FOR by the way it reports errors.  |
c +-----------------------------------------------------------+

c +====================================================+
c ||||||||||| Determination of internal unit |||||||||||
c |||||||||||    normal to crystal surface   |||||||||||
c +====================================================+
        Real*8  Base_Normal(3),
     *          Projection(3),
     *          Surface_Normal(3),
     *          Disorient8,
     *          Dpi8, Dgra8, s8

        Real*8   DAngle, DAngle2, DVecSca2
        External DAngle, DAngle2, DVecSca2

        Integer nBaseNormal(3),
     *          nMiscutDirection(3),
     *          ifail, i,
     *          ii/3/                                   !Number of Bragg planes indices (3 or 4)


        Character       txt(20)*80
        Common  /msg/   txt

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
          progname = 'SurfNormVect'       !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+
c-------------------------------------------------------
        ifail = 0
        if (nBaseNormal(1).eq.0 .and.
     *      nBaseNormal(2).eq.0 .and.
     *      nBaseNormal(3).eq.0)        goto 100
        do      i=1,ii  !=========================+
          Base_Normal(i) = Real(nBaseNormal(i))  !|
        enddo   !=================================+
        Call    DUnitVec2 (Base_Normal,Base_Normal,ii)

        if (Abs(Disorient8).lt.1E-10) Then  !-------------------------+ if zero
c Unit vector along internal surface normal:                          |
          Call  DVecCon (Base_Normal,-1.D0,                          !|
     *                  Surface_Normal,ii)                           !|
        else  !-------------------------------------------------------+
          if (nMiscutDirection(1).eq.0 .and.                         !|
     *        nMiscutDirection(2).eq.0 .and.                         !|
     *        nMiscutDirection(3).eq.0)                    goto 101  !|
          s8 = Dabs(Dangle(nMiscutDirection,nBaseNormal,ii))         !|
          if ((s8.lt.1.0D0) .or. (s8.gt.179.0D0))          goto 102  !|
c+====================================================+               |
c|Calculate the normal; We seek it in the form:       |               |
c|             ->            ->               ->      |               |
c|      Surface_Normal = Base_Normal + c * Projection,|               |
c|                                                    |               |
c|          ->                            ->          |               |
c|where Projection - projection of Disorient_Reper on |               |
c|                            ->                      |               |
c| the plane, normal to   Base_Normal.                |               |
c|                                                    |               |
c|Therefore:                                          |               |
c|                                                    |               |
c|(Surface_Normal*Projection) =                       |               |
c|          = (Base_Normal*Projection) + c            |               |
c|          = c - Surface_Normal*sin(Disorient8) ,    |               |
c|                                                    |               |
c|(Surface_Normal*Base_Normal) =                      |               |
c|          = 1 + c*(Projection*Base_Normal)          |               |
c|          = 1 - Surface_Normal*cos(Disorient8) ,    |               |
c|                                                    |               |
c|            c = tg(Disorient)                       |               |
c+====================================================+               |
          Dpi8  = (4.D0)*Datan(1.D0)                                 !|
          Dgra8 = (2.D0)*Dpi8/(360.D0)                               !|
c Calculate the projection of chosen reper on the plane               |
c perpendicular to the base normal:                                   |
          Call DProject(nBaseNormal,nMiscutDirection,Projection,ii)  !|
          Call DUnitVec2(Projection,Projection,ii)                   !|
          s8 = Dtan(Disorient8*Dgra8)                                !|
          Call DVecSum2(Base_Normal, 1.D0,                           !|
     *                  Projection,  s8,                             !|
     *                  Surface_Normal, ii)                          !|
          Call DUnitVec2(Surface_Normal,Surface_Normal,ii)           !|
c Verification:                                                       |
          s8 = DVecSca2(Surface_Normal,Projection,ii)                !|
          s8 = Dangle2(Surface_Normal,Base_Normal,ii) * s8/Dabs(s8)  !|
c If verification is not passed:                                      |
          if (Dabs(s8-Disorient8).gt.1.) goto 103 !-------------------+---+
c Unit vector along internal surface normal:                          |   v
        Call   DVecCon (Surface_Normal,-1.D0,                        !|
     *                  Surface_Normal, ii)                          !|
        endif  !------------------------------------------------------+

        ifail = 0
  29    continue  !<---------------------------------------------------+
        if (iirezv.eq.1)  Then  !----------+                           |
          progname = progrezv(istackrezv) !|                           |
          istackrezv = istackrezv-1       !|                           |
        endif  !---------------------------+                           |
        return                                                        !|
                                                                      !|
c-----------------------------------------------------------           |
c                          ERRORS:                                     |
c-----------------------------------------------------------           |
  100   continue                                                      !|
        ifail = -5                                                    !|
        write   (txt,200)                                             !^
  200   format  ('SurfNormVect  E R R O R:'//
     *          'Unable to determine surface normal:'//
     *          '** zero normal vector **')                           !^
        goto 29  !-----------------------------------------------------+

  101   continue
        ifail = -5
        write   (txt,201)
  201   format  ('SurfNormVect  E R R O R:'//
     *          'Unable to determine surface normal:'//
     *          '** zero miscut vector **')                           !^
        goto 29  !-----------------------------------------------------+

  102   continue
        ifail = -5
        write   (txt,202)
  202   format  ('SurfNormVect  E R R O R:'//
     *          'Unable to determine surface normal:'//
     *          '** miscut direction is along the base normal **')    !^
        goto 29  !-----------------------------------------------------+

  103   continue
        ifail = -5
        write (txt,203) s8, Disorient8
  203   format  ('SurfNormVect  E R R O R:'//
     *  ' Error computing disorientiented normal'//
     *  ' Computed miscut=',g10.3,'  Expected miscut=',g10.3)         !^
        goto 29  !-----------------------------------------------------+
        end
