        Subroutine ComprSp (txt,nchar)
c-----------------------------------------------------
c     Replaces all repeated spaces by singles
c nchar  - the length of txt string after the compression.
c-----------------------------------------------------
        Integer         nchar, i, is, js, nshift
        Character       txt*(*)

        nchar = Len_Trim(txt)
c Remove tabs and unprintable characters:
        do i=1,nchar  !========================+
          if (txt(i:i).lt.' ') txt(i:i)=' '   !|
        enddo  !===============================+
        is = 1
  2     continue   !<-------------<------------<-----+
        js = is                                     !|
        if (nchar.le.1) return                      !|
        is = Index(txt(is:nchar),' ')               !|
        if (is.eq.0)    return                      !^
        nshift =0                                   !|
        is = is+js                                  !|
        do i=is,nchar  !=====================+       |
          if (txt(i:i).ne.' ')  goto 1  ! ---+-->-+  |
          nshift=nshift+1                   !|    |  |
        enddo  !=============================+    |  |
  1     continue  !<------------------------------+  |
        if (nshift.gt.0) Then !-------------------+  ^
          do i=is,nchar-nshift  !=============+   |  |
            txt(i:i)=txt(i+nshift:i+nshift)  !|   |  |
          enddo  !============================+   |  |
          txt(nchar-nshift+1:nchar)=' '          !|  |
          nchar=nchar-nshift                     !|  |
        endif  !----------------------------------+  |
        goto 2  !----->----------->----------->------+
        end
