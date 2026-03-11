        Subroutine HTMLstop (lun,lines)
c This sub is called from x0h_form and x0p_form and it prints
c fatal error messages on HTML page. The error text is taken
c from the txt(20)*80.
        Integer   lun, i, l, lines
        Character txt(20)*80
        Common  /msg/   txt

c This is table-in table (the upper-level table header is
c specified in x0h_form):

c 52    format(' <center><font size=+1>'/
c    *        ' <table border bgcolor="#ffffff">')
        write (lun,1)
  1     format(
     *  ' <tr><td><table border=0 cellspacing=16>'/
     *  ' <tr><td><img src="images/stop1.gif" width=64 height=64',
     *                          1x,'border=0 valign="center"></td>'/
     *      ' <td><center><font size=+1>&nbsp;<br>')
        do i=1,lines  !====================+
          l = Max (Len_Trim(txt(i)),1)    !|
          write (lun,2) txt(i)(1:l)       !|
  2       format(1x,a,'&nbsp;<br>')       !|
        enddo  !===========================+
        write (lun,3)
  3     format(' </font></center></td></tr>'/
     *         ' </table></td></tr>')
c       write (lun,98)
c 98    format(' </table></font></center>')
        return
        end

c =================================================================

        Subroutine HTMLstop_noJPG (lun,lines)
c This sub does absolutely the same as HTMLstop, but the error
c message does not contain the stop-sign image. This sub is
c used to produce a copy of output HTML for log purposes.
        Integer   i, l, lun, lines
        Character txt(20)*80
        Common  /msg/   txt

c This is table-in table (the upper-level table header is
c specified in x0h_form):

c 52    format(' <center><font size=+1>'/
c    *         ' <table border bgcolor="#ffffff">')
        write (lun,1)
  1     format(
     *  ' <tr><td><table border=0 cellspacing=16>'/
     *  ' <tr><td><center><font size=+1>&nbsp;<br>')
        do i=1,lines  !====================+
          l = Max (Len_Trim(txt(i)),1)    !|
          write (lun,2) txt(i)(1:l)       !|
  2       format(1x,a,'&nbsp;<br>')       !|
        enddo  !===========================+
        write (lun,3)
  3     format(
     *  ' </font></center></td></tr>'/
     *  ' </table></td></tr>')
c       write (lun,98)
c 98    format(' </table></font></center>')
        return
        end
