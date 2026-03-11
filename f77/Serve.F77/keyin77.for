        Subroutine KeyIn (m,is,ia)

c This is GNU fortran version
c-------------------------------------------------------
c       This returns scan (is) and ASCII (ia) codes
c                of a pressed keyboard key.
c
c  m=0 -- do no wait for the pressed key
c  m=1 -- wait for the pressed key
c
c Unfortunately, the "m" option is ignored: with GNU Fortran
c the function waits for ENTER to get the pressed key.
c A way around it may be implementing it in C:
c
c void initscr_() { initscr(); }
c void getch_() { getch(); }
c
c Then in Fortran:
c call initscr()
c call getch()
c
c Then compiling:
c gcc cstuff.c -fno-leading-underscore -c
c g77 keyIn.for -c
c g77 keyIn.o cstuff.o - o prgname
c-------------------------------------------------------
        character*1 c
        Integer   m,is,ia
        Integer*1 a2s(127) /
     *  30,48,46,32,18,33,34,14,15,36,37,38,28,49,24,25,
     *  16,19,31,20,22,47,17,45,21,44, 1,43,27, 7,12,57,
     *   2,40, 4, 5, 6, 8,40,10,11, 9,13,51,12,52,53,11,
     *   2, 3, 4, 5, 6, 7, 8, 9,10,39,39,51,13,52,53, 3,
     *  30,48,46,32,18,33,34,35,23,36,37,38,50,49,24,25,
     *  16,19,31,20,22,47,17,45,21,44,26,43,27, 7,12,41,
     *  30,48,46,32,18,33,34,35,23,36,37,38,50,49,24,25,
     *  16,19,31,20,22,47,17,45,21,44,26,43,27,41,14/

        is = 0
        ia = 0
        if (m .eq. 0) return

        call fget(c)                            !GNU fortran function
        ia = ichar (c)
        if   (ia.ne.0 .AND. ia.ne.224) Then !----+ #E0
          if (ia.le.127)  is = a2s(ia)          !|
        else  !----------------------------------+
          ia = 0                                !|
          is = Ichar (c)                        !|
        endif  !---------------------------------+
        return
        end
