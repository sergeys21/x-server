        Subroutine Navigate_Site(lunout,referer,ifail)
        Integer         lunout, ifail, lr, ls, lu, i
        Character       referer*(*), url*64

        Character       Get_Server_URL*64
        External        Get_Server_URL

c Here referer is used only when this sub is called from wwwwatch.pl and
c in this case "referer" is the web form to return for "Another calculation".
c In all other cases pass empty string!!!
        lr = Len_Trim(referer)
        ls = lr
        if (lr.gt.0) Then  !----------------------------+
          if (Index(referer(1:lr),'!').gt.0 .OR.       !|
     *        Index(referer(1:lr),'@').gt.0 .OR.       !|
     *        Index(referer(1:lr),'#').gt.0 .OR.       !|
     *        Index(referer(1:lr),'$').gt.0 .OR.       !|
c    *        Index(referer(1:lr),'?').gt.0 .OR.       !|may be in URL, e.g. /cgi/www_form.pl?template=x0h_form.htm
c    *        Index(referer(1:lr),'&').gt.0 .OR.       !|may be in URL, e.g. /cgi/www_form.pl?template=x0h_form.htm
c    *        Index(referer(1:lr),'%').gt.0 .OR.       !|may be in URL, e.g. /cgi/www_form.pl?template=x0h_form.htm
c    *        Index(referer(1:lr),'=').gt.0 .OR.       !|
     *        Index(referer(1:lr),'^').gt.0 .OR.       !|
     *        Index(referer(1:lr),'*').gt.0 .OR.       !|
     *        Index(referer(1:lr),'(').gt.0 .OR.       !|
     *        Index(referer(1:lr),')').gt.0 .OR.       !|
     *        Index(referer(1:lr),'-').gt.0 .OR.       !|
     *        Index(referer(1:lr),'+').gt.0 .OR.       !|
     *        Index(referer(1:lr),'[').gt.0 .OR.       !|
     *        Index(referer(1:lr),']').gt.0 .OR.       !|
     *        Index(referer(1:lr),'<').gt.0 .OR.       !|
     *        Index(referer(1:lr),'>').gt.0 .OR.       !|
     *        Index(referer(1:lr),'{').gt.0 .OR.       !|
     *        Index(referer(1:lr),'}').gt.0 .OR.       !|
     *        Index(referer(1:lr),'|').gt.0 .OR.       !|
     *        Index(referer(1:lr),';').gt.0 .OR.       !|
     *        Index(referer(1:lr),'"').gt.0 .OR.       !|Double quote
     *        Index(referer(1:lr),char(39)).gt.0 .OR.  !|Single quote
c    *        Index(referer(1:lr),char(9)).gt.0  .OR.  !|Horizontal Tab
c    *        Index(referer(1:lr),char(10)).gt.0 .OR.  !|Line Feed
c    *        Index(referer(1:lr),char(13)).gt.0 .OR.  !|Carriage Return
     *        Index(referer(1:lr),char(92)).gt.0 .OR.  !|backslash "\"
     *        Index(referer(1:lr),' ').gt.0) Then !--+  |
            referer = ' '                           !|  |
            ls = 0                                  !|  !
            goto 5  !--------------------------------+--+-----+
          endif  !-----------------------------------+  |     |
          do i=1,lr !---------------------------------+ |     v
            if (ichar(referer(i:i)).lt.32 .OR.       !| |
     *          ichar(referer(i:i)).gt.122) Then !-+  | |
              referer = ' '                       !|  | |
              ls = 0                              !|  ! |
              goto 5  !----------------------------+--+-+-----+
            endif  !-------------------------------+  | |     |
          enddo !-------------------------------------+ |     v
        endif  !----------------------------------------+
  5     continue
c We need the table tag only if we passed some non-empty referrer
        if (lr .gt. 0) Then !------------------------------------------+
          write (lunout,1,err=10)                                     !|
  1       format(/' <hr>'/' <table>')                                 !|
        endif !--------------------------------------------------------+
c Add "back to referrer" only if it is good/legit
        if (ls .gt. 0) Then !------------------------------------------+
          write (lunout,2,err=10) referer(1:lr), referer(1:lr)        !|
  2       format(' <tr valign="middle">'/
     *    ' <td><a href="',a,'"><img src="images/3darrow.gif"',
     *                    1x,'width=32 height=25 border=0></a></td>'/
     *    ' <td><a href="',a,'"><b>Another calculation</b></a></td>'/
     *    ' </tr>')                                                   !|
        endif  !-------------------------------------------------------+

        url = Get_Server_URL()
        lu = Len_Trim(url)
        if (lu .eq. 0) Then !----+
          url = '/'             !|
          lu = 1                !|
        endif !------------------+
        write (lunout,3,err=10) url(1:lu), url(1:lu),
     *                          url(1:lu), url(1:lu),
     *                          url(1:lu), url(1:lu),
     *                          url(1:lu), url(1:lu)
  3     format(' <tr valign="middle">'/
     *  ' <td><img src="images/x.gif" width="31" height="32"',
     *                                            1x,'border="0"></td>'/
     *  ' <td><form><b>Site navigation:</b>'/
     *  ' <select Name="list">'/
     *  ' <option Value="',a,'index.shtml">Home page</option>'/
     *  ' <option Value="',a,'x0h.html">X-ray scattering',
     *                                      1x,'factors (X0h)</option>'/
     *  ' <option Value="',a,'x0h_search.html">Search for X-ray Bragg',
     *                                1x,'planes (X0h-search)</option>'/
     *  ' <option Value="',a,'GID_sl.html">X-ray Bragg',
     *                               1x,'diffraction (GID_sl)</option>'/
     *  ' <option Value="',a,'TER_sl.html">X-ray specular',
     *                                1x,'reflection (TER_sl)</option>'/
     *  ' <option Value="',a,'MAG_sl.html">X-ray resonant',
     *                                1x,'reflection (MAG_sl)</option>'/
     *  ' <option Value="',a,'TRDS_sl.html">X-ray scattering from',
     *                                1x,'roughness (TRDS_sl)</option>'/

     *  ' <option Value="',a,'BRL.html">Multiple Bragg/Laue',
     *                                  1x,'diffraction (BRL)</option>'/

     *  ' </select>'/
     *  ' <input Type="button" Value="Go!"',1x,
     *           'onClick="window.location.href=this.form.list.options',
     *                         '[this.form.list.selectedIndex].value">'/
     *  ' </form></td>'/
     *  ' </tr>')
        if (lr .gt. 0) Then !----------+
          write(lunout,4,err=10)      !|
  4       format('</table>')          !|
        endif !------------------------+

        ifail = 0
        return

  10    continue
        ifail = 1
        return
        end
