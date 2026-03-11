c=============================================================
c This is a new-style BRAGAN with Single precision (Bragan,d).
c
c It can be used with Back Diffraction Programs because it
c takes QB=90 even if lambda/2d ~ 1.+OVER
c=============================================================
        Real Function  Bragan   (wave,                          !input,  real
     *                           ind,inmax,iwarn,               !input,  integer
     *                           d,                             !output, real
     *                           isico,                         !input,  integer
     *                           sin1,cos1,tg1,ctg1,sin2,cos2)  !output, real
c-----------------------------------------------------
c  Calculate Bragg angle in a crystal with arbitrary symmetry
c
c iwarn:  What to do when Bragg condition cannot be met (sin(QB)>1)
c         0 - take QB = 0.
c         1 - take QB = 0. and warn.
c         2 - take sin(QB) = 1., QB=90.
c         3 - take sin(QB) = 1., QB=90. and warn.
c         4 - ask what to take.
c
c isico:  0 - do not calculate sin,cos,... (faster).
c         1 - calculate sin,cos,... .
c-----------------------------------------------------
        Real*8          BIG, SMALL, OVER, Bragan8
c       Parameter      (BIG   = 1.0D+37,                !we assign them below
c    *                  SMALL = 1.0D-37,
c    *                  OVER  = 1.0D-3)                 !~100*|x0| or 10eV at 10KeV
        Real*8          d8, s8, f8, x8, gra8, w8
        Real            wave, d,
     *                  sin1, cos1, tg1, ctg1,
     *                  sin2, cos2
        Integer         inmax, ind(inmax), isico,
     *                  iwarn, ifail
c -------------------------------------------------------
        Real*8          DVecMod
        External        DVecMod
c -------------------------------------------------------
        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
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
          progname = 'Bragan'             !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+
c -------------------------------------------------------
        BIG   = 1.0D+37
        SMALL = 1.0D-37
        OVER  = 1.0D-3                  !~100*|x0| or 10eV at 10KeV
c -------------------------------------------------------
        Bragan8 = 0.0D0
        ifail = 0
        gra8 = Datan(1.0D0)/45.         !number of radians in 1 degr.
        d8 = 1./DVecMod(ind,inmax)      !Bragg planes spacing (Angstrom)
        d  = Sngl(d8)
        s8 = 1./(2.*d8)                 !s = sin(QB)/L = 1/2d
        f8 = wave*s8                    !f  = sin(QB) = L/2d
        x8 = 1.0D0 - f8                 !must be positive: 1-sin(QB)=(1-L/2d)

        if (x8 .lt. (-OVER))  Then  !-------------------------+ 1-sin(QB)=(1-L/2d) < -1e-5
c+=========================================+                  |
c|             Case: sin(QB) > 1.          |                  |
c+=========================================+                  |
c Act depending on iwarn:                                     |
          goto (1,2,3,4,5)  iwarn+1 !--------->-----------+   |
c ___iwarn:     0 1 2 3 4                                 |   |
          stop    'Bragan: illegal IWARN'                !|   |
                                                         !|   |
c iwarn=0 -- take QB = 0:                                 |   |
  1       continue       !<-------------------------------+   |
          Bragan8 = 0.0D0                                !|   |
          goto 10                                        !|   |
                                                         !|   |
c iwarn=1 --  take QB = 0 and warn:                       |   |
  2       continue       !<-------------------------------+   |
          Bragan8 = 0.0D0                                !|   |
          Call Bragan_Warn (f8,ind,inmax)                !|   |
          if (modebat.eq.0) Then  !-----------------+     |   |
            Call Message (txt,5,2)                 !|     |   |
          else   !----------------------------------+     |   |
            ifail = 5                              !|     |   |this may be used if we want to
          endif  !----------------------------------+     |   |add ifail to Bragan call. Now
          goto 10                                        !|   |this warning is ignored in BAT-mode
                                                         !|   |
c iwarn=2 --  take QB=90:                                 |   |
  3       continue       !<-------------------------------+   |
          Bragan8 = 90.0D0                               !|   |
          goto 10                                        !|   |
                                                         !|   |
c iwarn=3 --  take QB=90 and warn:                        |   |
  4       continue       !<-------------------------------+   |
          Bragan8 = 90.0D0                               !|   |
          Call Bragan_Warn (f8,ind,inmax)                !|   |
          txt(5)=' ** Accepting: <wave/2d=1,QB=90.> **'  !|   |
          if (modebat.eq.0) Then  !-----------------+     |   |
            Call Message (txt,5,2)                 !|     |   |
          else   !----------------------------------+     |   |
            ifail = 5                              !|     |   |this may be used if we want to
          endif  !----------------------------------+     |   |add ifail to Bragan call. Now
          goto 10                                        !|   |this warning is ignored in BAT-mode
                                                         !|   |
c iwarn=4 -- ask what to take:                            |   |
  5       continue    !<----------------------------------+   |
          Call Bragan_Warn (f8,ind,inmax)                    !|
c Take same action as iwarn=0 -- take QB = 0:                 |
          Bragan8 = 0.0D0                                    !|
  10      continue                                           !|
                                                             !|
        elseif (x8.ge.(-OVER) .AND. x8.le.(0.D0))  Then !-----+ -1E-5 <= 1-sin(QB)=(1-L/2d) <= 0
                                                             !|
c+==========================================+                 |
c|         Case: Bragg angle = 90 degr:     |                 |
c+==========================================+                 |
          Bragan8 = 90.0D0                                   !|
                                                             !|
        elseif (x8.gt.(0.D0) .AND. x8.lt.(OVER))  Then !------+ 0 < 1-sin(QB)=(1-L/2d) < 1E-5
                                                             !|
c+==========================================+                 |
c|         Case: Bragg angle ~ 90 degr:     |                 |
c+==========================================+                 |
c sin(QB) = cos(90-QB) = 1-(90-QB)^2/2 = L/2d                 |
c (90-QB)^2/2 = 1- L/2d = x8                                  |
c (90-QB) = sqrt(2*x8)                                        |
          Bragan8 = 90.0D0 - Dsqrt(2.*x8)/gra8               !|
                                                             !|
        else  ! ----------------------------------------------+ 1-sin(QB)=(1-L/2d) >= 1E-5
                                                             !|
c+===============================================+            |
c|   Normal case: Bragg angle < 90 degr:         |            |
c+===============================================+            |
c Bragan is the calculated Bragg angle:                       |
          Bragan8 = Dasin(f8)/gra8                           !|
                                                             !|
        endif ! ----------------------------------------------+

c If need to calculate sin, cos,...:
        if (isico.gt.0) Then  !-----------------+
          if     (Bragan8.le.0.0D0) Then !----+ |
            sin1 = 0.                        !| |
            cos1 = 0.                        !| |
            tg1  = 0.                        !| |
            ctg1 = 0.                        !| |
            sin2 = 0.                        !| |
            cos2 = 0.                        !| |
          elseif (Bragan8.ge.90.0D0) Then !---+ |
            sin1 = 1.                        !| |
            cos1 = 0.                        !| |
            tg1  = Sngl(BIG)                 !| |
            ctg1 = 0.                        !| |
            sin2 = 0.                        !| |
            cos2 =-1.                        !| |
          else !------------------------------+ |
            sin1 = Sngl(f8)                  !| |
            w8  = 1.-f8*f8                   !| |
            if (w8.gt.0.) Then !-------+      | |
              cos1 = Sngl(Dsqrt(w8))  !|      | |
              tg1  = sin1/cos1        !|      | |
            else  !--------------------+      | |
               cos1 = 0.              !|      | |
               tg1  = Sngl(BIG)       !|      | |
            endif  !-------------------+      | |
            if (sin1.gt.SMALL) Then !--+      | |
              ctg1 = cos1/sin1        !|      | |
            else  !--------------------+      | |
              ctg1 = Sngl(BIG)        !|      | |
            endif  !-------------------+      | |
            sin2 = 2.*sin1*cos1              !| |
            cos2 = cos1*cos1-sin1*sin1       !| |
          endif  !----------------------------+ |
        endif  !--------------------------------+

        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        Bragan = Sngl(Bragan8)
        return
        end

c  ==================================================================

        Subroutine      Bragan_Warn   (f8,ind,inmax)
        Real*8          f8
        Integer         inmax, ind(inmax), nb, i
        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
        txt(1) = 'W A R N I N G :'
        txt(2) = ' '
        txt(3) = 'Reflection ('
        i = Len_Trim(txt(3))
        Call  PakInt  (ind,inmax,txt(3)(i+1:80),nb)
        txt(3)(i+1+nb:80)=') does not exist for this wavelength:'
        write   (txt(4),1)   f8
  1     format  ('wave/2d =',f10.6,' > 1.')
        txt(5) = ' '
        return
        end
