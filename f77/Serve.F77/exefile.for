        Logical Function ExeFile (program)
c       use portlib
c       use msflib
c ---------------------------------------------------
c Checks whether the specified program exists in path
c and can be executed
c ---------------------------------------------------
        Integer    i, j, k, l, m, ii, i95
        Character  program*(*), prgcopy*64, sp*1
        Character  exebuf*4096
        Common /which/  exebuf

        Character       getOS*8
        External        getOS
c ---------------------------------------------------
c This is portability lib to check the access to a file:
c result = ACCESS (name, mode), where mode:
c r=readtest, w=writetest, x=exetest, blank=existtest
        Integer    Access
        External   Access

        ExeFile = .False.
        m = Len_Trim(program)
        if (m.eq.0)     goto 1

        Call GetEnv32 ('PATH',exebuf)
        l = Len_Trim(exebuf)
        if (l.eq.0)     goto 1
        do j=1,l  !================================+
          if (exebuf(j:j).eq.'\') exebuf(j:j)='/' !|
        enddo  !===================================+

        if (getOS() .eq. 'windows') Then !---+
c          sp = char(92)                    !| '\'
           sp = char(47)                    !| '/'
        else !-------------------------------+
           sp = char(47)                    !| '/'
        endif !------------------------------+

        i95 = Index (exebuf(1:l),sp)         !find directory separator
        if (getOS() .eq. 'windows' .AND. i95 .eq. 0) Then !-+
        i95 = Index (exebuf(1:l),char(47))                 !|find ALTERNATIVE directory separator
        endif !---------------------------------------------+

c For Win95/NT check the current directory before checking the path:
        if (i95.gt.0) Then  !--------------------------+
          if (Access(program(1:m),' ').eq.0) Then !--+ |if program(1:m) exists
            ExeFile = .True.                        !| |
            exebuf = program(1:m)                   !| |
            return                                  !| |
          endif  !-----------------------------------+ |
        endif  !---------------------------------------+

        i = 1
  3     continue  !<----------------------<----------------------------+
                                                                      !|
        if (i95.gt.0) Then  !-----------------+                        ^
          j = Index (exebuf(i:l),char(59))   !|WIN: ';' path separator |
        else  !-------------------------------+                        |
          j = Index (exebuf(i:l),char(58))   !|UNIX:':' path separator |
        endif  !------------------------------+                        |
        if (j.eq.0) Then !-+                                           |
          j = l+1         !|                                           |
        else  !------------+                                           |
          j = j+i-1       !|                                           |
        endif  !-----------+                                           |
        exebuf(j:j) = sp                        !Add path separator    |
                                                                      !^
c "prgcopy" is used instead of "program" because concatenation with    |
c character*(*) is not allowed in DOS Fortran-5.1:                     |
        prgcopy = program                                             !|
c       write (*,2) exebuf(i:j)//prgcopy(1:m)                         !|
c 2     format(" Checking: ",a/)                                      !|
c       Call Keyin(1,iscan,iascii                                     !|
                                                                      !|
c If program exists:                                                   ^
        if   (Access(exebuf(i:j)//prgcopy(1:m),' ').eq.0) Then !----+  |
c If program executable:                                            |  |
          if (Access(exebuf(i:j)//prgcopy(1:m),'x').eq.0) Then !-+  |  |
            ExeFile = .True.                                    !|  |  |
c Shift to the left:                                             |  |  |
            k = i-1                                             !|  |  |
            do ii=i,j  !===========================+             |  |  |
              exebuf(ii-k:ii-k) = exebuf(ii:ii)   !|             |  |  |
            enddo  !===============================+             |  |  |
c Add the program name:                                          |  |  ^
            k = j-i+1                                           !|  |  |
            do ii=1,m  !===========================+             |  |  |
              exebuf(ii+k:ii+k) = program(ii:ii)  !|             |  |  |
            enddo  !===============================+             |  |  |
            k = k+m+1                                           !|  |  |
            l = Len(exebuf)                                     !|  |  |
            exebuf(k:l) = ' '                                   !|  |  |
            return                                              !|  |  |
          endif  !-----------------------------------------------+  |  |
        endif  !----------------------------------------------------+  |
        i = j+1                                                       !|
        if (i.le.l)     goto 3 !---------------------------------------+

  1     continue
        exebuf = ' '
        return
        end
