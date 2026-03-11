c ATTENTION::: the full optimization does not work here with the
c -----------  MS Fortran Power Station-4.0. Therefore, the
c              following metacommand is used to suppress it:
c
!MS$OPTIMIZE:'ON:xp'
c
c
c Available optimization options in MS Fortran Power Station:
c /Od   Disables optimization (default).
c /Op   Enables improved floating-point consistency.
c /Ox   Enables speed optimization and denoted inlining,
c       without run-time math error checking.
c/Oxp   Enables speed optimization and denoted inlining,
c       with run-time math error checking.
c =====================================================================

c Contents:  real function ForFac
c            subroutine    DisInt
c            subroutine    Sort_differences
c            subroutine    Gaus
c            subroutine    Disper
c            real*8 function Parratt
c            subroutine    CroSec
c            real function FiFunc
c            subroutine    Sort_ABC
c            subroutine    Sort_Index
c
c =======================================================
c
        Real Function ForFac (s,ifail)
c -------------------------------------------------------
c Calculation of x-ray scattering amplitudes for crystal atoms with
c the help of tables contained in [1-3]. The caclulations use the
c following formulae:
c
c .......................................................
c Method-1 (4-term formula chosen by ab(i,11)=0):
c                4           -b(j)*s*s
c       fh(s) = sum {a(j) * e         } + c               (1)
c               j=1
c Parameters a,b,c are contained in "atom.x0h" in the array ab(i,j):
c ab(i,1),ab(i,3),ab(i,5),ab(i,7) - corresponds to a;
c ab(i,2),ab(i,4),ab(i,6),ab(i,8) - corresponds to b;
c ab(i,9) - corresponds to c -- all from the tables in [1-3].
c .......................................................
c Method-2 (5-term formula chosen by ab(i,11)#0):
c                5           -b(j)*s*s
c       fh(s) = sum {a(j) * e         } + c               (2)
c               j=1
c Parameters a,b,c are contained in "atom.x0h" in the array ab(i,j):
c ab(i,1),ab(i,3),ab(i,5),ab(i,7),ab(i,9)  - corresponds to a;
c ab(i,2),ab(i,4),ab(i,6),ab(i,8),ab(i,10) - corresponds to b;
c ab(i,11) - corresponds to c -- all from more modern tables
c (see e.g. the DABAX server at ESRF).
c .......................................................
c s=sin(Qb)/wavelength - the parameter used in formulae (1-2) and
c expressed in (1/Angstrom).
c .......................................................
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
c +------------------------------------------+
        Real*8          exponent_max /100./ !|
c +------------------------------------------+
        Real*8          ph8, f8
        Real            s
        Integer         ifail, iab, i, j

        Character       name(kcompMax)*4
        Common  /x0pa2/ name

        Real            wave,df1(kcompMax),df2(kcompMax),
     *                  ab(kcompMax,nabcMax),s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,
     *                  Planck_,Boltzman,Avogadro,
     *                  wedg(nedgMax),gedg(nedgMax)
        Integer         indx,ipr,nedg
        Common  /x0pa4/ wave,df1,df2,ab,s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,indx,ipr,
     *                  Planck_,Boltzman,Avogadro,wedg,gedg,nedg

        Integer         luwa,luat,lusp,luco
        Common  /lunds/ luwa,luat,lusp,luco

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
          progname = 'x0h/Forfac'         !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+

c Prevents GNU fortran -Wextra warnings:
        i = kcompMax
        i = kolMax
        i = ncodMax
        i = natmMax
        i = nwavMax
        i = nedgMax
        i = nabcMax
        i = nedgTb
        i = maxX0hWarnLines
c -------------------------------------------------------
        ifail = 0       !this is reserved for future extensions

        if (abs(ab(indx,11)).gt.1.E-37) Then  !--+
          iab = 11                              !|
        else  !----------------------------------+
          iab = 9                               !|
        endif  !---------------------------------+

        f8 = ab(indx,iab)
        do      j=1,iab-2,2  !=======================+
          ph8 = ab(indx,j+1)*s*s                    !|
          if (Dabs(ph8).lt.exponent_max)  Then !--+  |
            f8 = f8 + ab(indx,j)*dexp(-ph8)      !|  |
          endif !---------------------------------+  |
        enddo  !=====================================+

c Account for dispersion corrections calculated in subroutine disper:
        ForFac = sngl(f8 + df1(indx))

        if (ipr.eq.2) Then !------------------------------------+
          write (txt,120)       name(indx), s, ForFac          !v
  120     format('Element: ',a,' -- FORFAC results'/
     *           5x,'s=',g13.6,'  fh(s)=',g13.6)               !^
          Call  Priblo  (txt,2,1,0,0,i)                        !|
        endif  !------------------------------------------------+

        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
c
c =======================================================
c
        Subroutine DisInt (dfv1,dfv2,wvdata,ifail)
c--------------------------------------------------------
c Calculate oscillator strengths 'gedg' for subroutine Disper
c which calculates the dispersion correction.
c The calculated array 'gedg(1-nedg)' is returned in common/x0pa4/.
c
c dfv1(5),dfv2(5) [INPUT]- tabular values of the 1st and 2nd dispersion
c                          corrections for 5 fixed wavelengths from
c                          D.T.Cromer, D.Liberman - J.Chem.Phys. (1970),
c                          v.53, No.5, p.1891-1898, or equivalent.
c wvdata(5)       [INPUT]- corresponding fixed wavelengths from the
c                          same tables.
c
c Disint interpolated dispersion corrections for an arbitrary wavelength
c in the X-ray range. It assumes that on a limited interval the
c following formulaer from Disper remain valid:
c
c       nedg
c df1 = sum ( gedg(k) * Parratt(X_k,n(k)) ),
c       k=1                                                   +------------------------+
c                                                             |ATTENTION: my paper with|
c       nedg                                                  |Olga Lugovskaya there is|
c df2 = sum ( 0.5*pi * (n(k)-1) * gedg(k) / X_k**(n(k)-1) ) . |a typo: instead of      |
c       k=1,                                                  |X_k=w_i/wedg it reads:  |
c    (X_k > 1)                                                |X_k=lambda_i/lambda_k   |
c                                                             +------------------------+
c where n(k) -       11/4 for (1s 1/2) electrons,
c                     7/3 for (2s 1/2) electrons,
c                     5/2 for all other electrons;
c Parratt(X_k,n(k)) - is the intergral calculatted by Parratt in:
c                     L.G.Parratt, C.F.Hempstead - Phys.Rev.
c                     (1954), 94, p.1593-1600.
c gedg(k)  -          the oscillator strengths for respective
c                     absorption edges.
c
c However, the oscillator strengths in these formulae are not the
c tabulated values, but they are calculated in Disint. For this
c purpose with the help of Sort Disint finds the closest fixed
c wavelengths wvdata(i) with known dispersion corrections and applies
c them the above formulae. This make a system of linear equations for
c sought oscillator strengths, which is solved by the Gaussian method.
c in subroutine Gaus. Once the oscillator strengths are found, they
c are applied for given X-ray wavelength instead of tabular values
c from "atom.x0h".
c--------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c       Parameter       (nedgTb   = 5)
c+==============================================+
        Include         'x0h_size.inc'         !|
c In particular, the above include              |
c contains this:                                |
c       Parameter       (nedgTb   = 5)          |
c       Parameter       (maxX0hWarnLines = 10); |
c+==============================================+
        Integer         nedg2x
        Parameter       (nedg2x = 2*nedgTb) !we may use not all 16 edges!

c In reality the maximum expected oscillator strength value should
c be about 15 (see the values of oscillator strengths in Atom.X0h):
        Real*8          oscillator_max  /1000./

        Real*8          a(nedg2x,nedg2x),  b(nedg2x),  x(nedg2x),
     *                  a_(nedg2x,nedg2x), b_(nedg2x), x_(nedg2x),
     *                  xmax
        Real            dfv1(nedgTb), dfv2(nedgTb),
     *                  wvdata(nedgTb), del(nedgTb),
     *                  n1
        Integer         i_old(nedgTb), i, j, k, k_, m,
     *                  ifail, jfail, iname

        Real            wave2energy
        External        wave2energy

        Character       name(kcompMax)*4
        Common  /x0pa2/ name

        Real            wave,df1(kcompMax),df2(kcompMax),
     *                  ab(kcompMax,nabcMax),s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,
     *                  Planck_,Boltzman,Avogadro,
     *                  wedg(nedgMax),gedg(nedgMax)
        Integer         indx,ipr,nedg
        Common  /x0pa4/ wave,df1,df2,ab,s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,indx,ipr,
     *                  Planck_,Boltzman,Avogadro,wedg,gedg,nedg

        Integer         linesX0hWarn
        Character       txtX0hWarn(maxX0hWarnLines)*80
        Common /x0hwarn/ txtX0hWarn, linesX0hWarn

        Integer         luwa,luat,lusp,luco
        Common  /lunds/ luwa,luat,lusp,luco

        Character       txt(20)*80
        Common  /msg/   txt

        Real*8          Parratt, X_k
        External        Parratt

c Prevents GNU fortran -Wextra warnings:
        i = kcompMax
        i = kolMax
        i = ncodMax
        i = natmMax
        i = nwavMax
        i = nedgMax
        i = nabcMax
        i = nedgTb
        i = maxX0hWarnLines

        iname = max (Len_Trim(name(indx)),1)

        do i=1,nedgMax !==+
          gedg(i) = 0.   !|
        enddo  !==========+

        if (nedg .gt. nedg2x) Then  !-------------------------+
c                                                             |
c+==================================================+         |
c| THIS IDEA IS WRONG:                              |         |
c|                                                  |         |
c Cutoff some edges:                                |         |
c -- first exclude those with X_k < 1               |         |
c    (no contribution to df")                       |         |
c         do k=1,nedg-nedg2x  !=================+   |         |
c           X_k = dble( wave2energy(wave)      !|   |         |
c    /                   / wedg(k) )           !|   |         |
c Here we also include the edges where X_k<1,   |   |         |
c but the edge is very close:                   |   |         |
c           if (X_k.gt.0.95) goto 1 !----+      |   |         |
c         enddo  !=======================+======+   |         |
c         k = nedg - nedg2x + 1         !v          |         |
c 1       continue  !<-------------------+          |         |
c         k = k-1                                   |        !|
c         do i=1,nedg-k  !========+                 |         |
c           wedg(i) = wedg(i+k)  !|                 |         |
c         enddo  !================+                 |         |
c |                                                 |        !|
c +=================================================+        !|
                                                             !|
c Truncate 'nedg':                                           !|
          nedg = nedg2x                                      !|
                                                             !|
          if (ipr .eq. 2) Then  !---------------------------+ |
            write (txt,100)  name(indx)(1:iname),          !| |
     *                       nedg2x, wave                  !| |
  100       format('DISINT: [',a,'] count of edges',1x,
     *             'truncated to ',i2,' @ wave=',g12.6)    !| |
            Call Priblo (txt,1,1,0,0,i)                    !| |
          endif  !------------------------------------------+ |
                                                             !|
        elseif (nedg .le. 0)  Then  !-------------------------+
                                                             !|
          goto 999  !------------------nothing-to-do----------+--+
                                                             !|  v
        endif  !----------------------------------------------+

c Calculate the differences
c between x-ray energy and
c the energies of absorption
c edges:
        do i=1,nedgTb  !==============================================+nedgTb<=5
          del(i)   = Abs(wave2energy(wave) - wave2energy(wvdata(i))) !|
          i_old(i) = i                                               !|
        enddo  !======================================================+

c Sort the differences in ascending order:
c (the i_old are sorted synchronously and
c contain the former numbers of del)
        Call Sort_differences (del, i_old, nedgTb)

c Initialize (set zero) all
c parameters used below:
        do m=1,nedg2x  !=========+
          b_(m) = 0.            !|
          x_(m) = 0.            !|
          do k=1,nedg2x  !====+  |
            a_(m,k) = 0.     !|  |
          enddo  !============+  |
        enddo  !=================+

        do      m=1,nedg  !=====================================+nedg<=10
          j = i_old((m-1)/2+1)                                 !|j<=5
          if (m .ne. 2*(m/2)) Then !-+                          |
            b_(m) = dfv1(j)         !|odd (1,3,5,...)           |
          else  !--------------------+                          |
            b_(m) = dfv2(j)         !|even (2,4,6,...)          |
          endif  !-------------------+                          |
          do    k=1,nedg  !===================================+ |nedg<=10
c X=wvdata(j)/wedg(k),                                        | |
c wi     - energy of incident X-rays in KeV:                  | |
c          wi[kev] = wave2energy(wave[a])                     | |
c wedg(k)- energies of absorption edges in KeV.               | |
            X_k = dble( wave2energy(wvdata(j))/wedg(k) )     !| |
            if (m .ne. 2*(m/2))  Then  !------odd(1,3,5..)-+  | |
             k_ = k                                       !|  | |
c+-------------------------------------------+             |  | |
c|        nedg                               |             |  | |
c|df1(j)= sum (gedg(k)*Parratt(X_k(j),n(k))) |             |  | |
c|        k=1                                |             |  | |
c+-------------------------------------------+             |  | |
             a_(m,k) = Parratt(X_k,k_)                    !|  | |
            else  !--------------------------even(2,4,6..)-+  | |
c The condition X_k>1 may cause the interpolation          |  | |
c to fail for soft x-rays where several a(m,k)=0           |  | |
              if (X_k. gt. 1.0D0) Then !-----------------+ |  | |
c n(k):  11/4 for (1s 1/2) electrons,                    | |  | |
c         7/3 for (2s 1/2) electrons,                    | |  | |
c         5/2 for all other electrons;                   | |  | |
                if     (k.eq.1) Then !--+                | |  | |
                   n1 = 11./4.-1.      !|                | |  | |
                elseif (k.eq.2) Then !--+                | |  | |
                   n1 = 7./3.-1.       !|                | |  | |
                else !------------------+                | |  | |
                   n1 = 5./2.-1.       !|                | |  | |
                endif  !----------------+                | |  | |
c+-----------------------------------------------------+ | |  | |
c|       nedg                                          | | |  | |
c|df2(j)=sum( 0.5*pi*(n(k)-1)*gedg(k)/X_k(j)^(n(k)-1) )| | |  | |
c|       k=1,                                          | | |  | |
c|   (X_k(j) > 1)                                      | | |  | |
c+-----------------------------------------------------+ | |  | |
                a_(m,k) = 0.5*pi*n1*dexp(-n1*dlog(X_k)) !| |  | |
              else !-------------------------------------+ |  | |
                a_(m,k) = 0.                            !| |  | |
              endif  !-----------------------------------+ |  | |
            endif  !---------------------------------------+  | |
          enddo  !============================================+ |
        enddo  !================================================+

c Make a copy of matrices since Gaus_8 destroys them:
  1     continue  !<-------------------------------------------------+
        do      m=1,nedg  !=========+nedg<=10                        |
          b(m) = b_(m)             !|                                |
          x(m) = x_(m)             !|                                |
          do    k=1,nedg  !=======+ |nedg<=10                        |
            a(m,k) = a_(m,k)     !| |                                |
          enddo  !================+ |                                ^
        enddo  !====================+                                |
                                                                    !|
c       if (ipr.eq.2) Then  !-------------------------------------+  |
c         write (txt,401) name(indx)(1:iname, nedg,              !|  |
c    *                    (b(i),i=1,nedg2x)                      !|  |
c 401     format('Element: ',a,' -- Interpolation eqn.:'/        !|  |
c    *           10x,'-- Right side (first  ',i2,' elements):'/  !|  |
c    *           2(5g12.3,:/))                                   !|  ^
c         Call Priblo (txt,4,1,0,0,i)                            !|  |
c         write (txt,402) nedg, nedg, ' elements):'              !|  |
c 402     format(10x,'-- Left side (first  ',i2,' x ',i2, a)     !|  |
c         Call Priblo (txt,1,0,0,0,i)                            !|  |
c         do m=1, nedg2x  !============================+          |  |
c           write (txt,403) m,(a(i,m),i=1,nedg2x)     !|          |  |
c           Call Priblo (txt,3,0,0,0,i)               !|          |  ^
c         enddo  !=====================================+          |  |
c 403     format(14x,'line=',i2/2(5g12.3,:/))                    !|  |
c       endif  !--------------------------------------------------+  |
                                                                    !|
c Solve linear equations: A*x=B                                      |
        Call    Gaus_8  (a,b,x,nedg2x,nedg,jfail)                   !|
                                                                    !|
        if (jfail .eq. 0) Then !----------------------------+        |
c Verify for wrong values because                           |        |
c of possible precision losses:                             |        ^
          xmax = 0.                                        !|        |
          do i=1,nedg  !=================================+  |        |
            if (dabs(x(i)) .gt. xmax) xmax = dabs(x(i)) !|  |        |
          enddo  !=======================================+  |        |
c In reality the maximum expected value should be about 15  |        |
c(see the values of oscillator strengths in Atom.X0h):      |        |
          if (xmax .gt. oscillator_max) jfail = -1         !|        |
        endif !---------------------------------------------+        ^
                                                                    !|
        if (jfail.ne.0) Then  !-------------------------------+      |
          write (txt,200)  name(indx)(1:iname), nedg, wave   !v      ^
  200     format('DISINT: [',a,'] oscillator strengths',1x,
     *    'failed with nedg=',i2,' @ wave=',g12.6)           !^      ^
          if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+    |      |
            linesX0hWarn = linesX0hWarn+1               !|    |      |
            txtX0hWarn(linesX0hWarn) = txt(1)           !|    |      ^
          endif !----------------------------------------+    |      |
          if (ipr.eq.2) Call Priblo (txt,1,1,0,0,i)          !|      |
          if (nedg .gt. 1) Then !-+                           |      |
            nedg = nedg - 1      !|                           |      |
c Retry interpolation with a      |                           |      |
c smaller number of edges:        |                           |      |
            goto 1    !-----------+---------------------------+------+
          endif  !----------------+                           |
          ifail =-1                                          !|
          write (txt,201)  name(indx)(1:iname), wave         !v
  201     format('DISINT: [',a,'] oscillator strengths',1x,
     *    'failed with any nedg @ wave=',g12.6)              !^
          if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+    |
            linesX0hWarn = linesX0hWarn+1               !|    |
            txtX0hWarn(linesX0hWarn) = txt(1)           !|    |
          endif !----------------------------------------+    |
          goto 999  !-----------------------------------------+error+
        endif  !----------------------------------------------+     v

        do      i=1,nedg !=========+
          gedg(i) = sngl(x(i))    !|
        enddo  !===================+

        if (ipr.eq.2) Then  !----------------------------------------+
          write (txt,300) name(indx)(1:iname),                      !|
     *                    (wedg(i),gedg(i),i=1,nedg)                !|
300       format('Element: ',a,' -- Oscillator strengths by DISINT'/
     *           2x,'Edge_Energy     Oscillator_Strength'/
     *           10(2g16.6,:/))                                     !|
          Call Priblo (txt,2+nedg,1,0,0,i)                          !|
        endif  !-----------------------------------------------------+

  999   continue
        return
        end
c
c =======================================================
c
        Subroutine Sort_differences (del,ind,N)
c--------------------------------------------------------
c         Sort the wavelengths differences for Disint
c                  in the increasing order
c--------------------------------------------------------
        Integer N, k, m,
     *          ind(N), ind_k
        Real    del(N), del_k
        do      k = 2,N  !=======================+
          del_k = del(k)                        !|
          ind_k = ind(k)                        !|
                                                !|
          do    m = k-1,1,-1  !============+     |
            if (del(m).le.del_k) goto 3 !--+->+  |
            del(m+1) = del(m)             !|  |  |
            ind(m+1) = ind(m)             !|  |  |
          enddo  !=========================+  V  |
          m = 0                              !|  |
  3       continue   !<-----------------------+  |
          del(m+1) = del_k                      !|
          ind(m+1) = ind_k                      !|
        enddo  !=================================+
        return
        end
c
c =======================================================
c
        Subroutine Gaus_8 (a,b,x,m,n,jfail)
c--------------------------------------------------------
c Solve system of linear equations A*x=B for Disint
c (find oscillator strengths from the system of linear
c equations for calculating dispersion corrections)
c--------------------------------------------------------
        Integer m,n,n1,i,i1,j,k,h,jfail
        Real*8  a(m,m),b(m),x(m),v,r,
     *          eps /1.D-16/

        jfail = 0

        if (n.le.1)     Then  !-----------------+
          if (dabs(a(1,1)).lt.eps)  goto 135 !--+-+
          x(1) = b(1)/a(1,1)                   !| |
          return                               !| V
        endif  !--------------------------------+

        n1 = n-1

        do      i=1,n1  !===============================+
          h = i+1                                      !|
          if (dabs(a(i,i)).lt.eps)      Then !-------+  |
            do    k=h,n  !=========================+ |  |
              if (dabs(a(k,i)).gt.eps)  goto 30 !-+| |  |
            enddo  !==============================++ |  |
            goto 135  !-->------------------------+--+--+--+
  30        continue  !--<------------------------+  |  |  V
            do    j=i,n  !========+                  |  |
              v      = a(i,j)    !|                  |  |
              a(i,j) = a(k,j)    !|                  |  |
              a(k,j) = v         !|                  |  |
            enddo  !==============+                  |  |
            v    = b(i)                             !|  |
            b(i) = b(k)                             !|  |
            b(k) = v                                !|  |
          endif  !-----------------------------------+  |
  50      continue  !-<-----------------------------+   |
          if (dabs(a(h,i)).gt.eps) Then  !------+   |   |
          if (dabs(a(i,i)).lt.eps) goto 135  !--+---+---+-+
            r = -a(h,i)/a(i,i)                 !|   |   | V
            do      j=i,n  !=============+      |   |   |
              a(h,j) = a(h,j)+r*a(i,j)  !|      |   |   |
            enddo  !=====================+      |   ^   |
            b(h) = b(h)+r*b(i)                 !|   |   |
          endif  !------------------------------+   |   |
          h = h+1                                  !|   |
          if (h.le.n)   goto 50  !->----------------+   |
        enddo  !========================================+

        if (dabs(a(n,n)).lt.eps) goto 135 !---------------+
        x(n) = b(n)/a(n,n)                               !V
        i = n1

  110   continue   !-<------------------------+
        v  = 0.                              !|
        i1 = i+1                             !|
        do h=i1,n  !==========+               |
          v = v+a(i,h)*x(h)  !|               |
        enddo  !==============+               ^
        if (dabs(a(i,i)).lt.eps) goto 135  !--+-+
        x(i) = (b(i)-v) / a(i,i)             !| V
        i=i-1                                !|
        if (i.ge.1)     goto 110  !-----------+

        do i=1,n !==========================+
          if (abs(x(i)).lt.eps) x(i) = 0.  !|
        enddo  !============================+
        return

  135   continue
        jfail = 1
        do i=1,n !====+
          x(i) = 0.  !|
        enddo  !======+
        return
        end
c
c =======================================================
c
        Subroutine Disper (atom_mass_aem, idis, jHenkeCowan,
     *                     dfv1, dfv2, lfix, nval, ifail)
c -------------------------------------------------------
c atom_mass_aem - input (atomic mass in AEM)
c idis          - input (dispersion mode, 1-3)
c                       1: Cromer formula for f" with known oscillator strengths.
c                          No reading Cowan or Henke tables in this case.
c                       2: Cromer formula for f" with interpolation of oscillator
c                          strengths using tabulated df',df" for 5 fixed x-ray lines
c                          (Cr-Ka1, Fe-Ka1, Cu-Ka1, Mo-Ka1, Ag-Ka1).
c                       3: Cromer formula for f" with interpolation of oscillator
c                          strengths using tabulated df',df" for 5 x-ray wavelengths
c                          listed for this element in atom.x0h.
c jHenkeCowan   - input (0=X0h,1=Henke(f1),   2=Henke(f1,f2),
c                              3=Cowan(f1),   4=Cowan(f1,f2)
c                              5=Windt(f1),   6=Windt(f1,f2)
c                              7=Chantler(f1),8=Chantler(f1,f2))
c dfv1          - input (f' for 5 or less waves)
c dfv2          - input (f" for 5 or less waves)
c lfix          - input (wavelength for 5 or less waves)
c nval          - input (how many waves available -- 5 or less)
c ifail         - output
c
c Results df1(indx) and df2(index) are returned via common /x0pa4/
c -------------------------------------------------------
c Calculate dispersion corrections using the following formula from
c D.T.Cromer - Acta Cryst.(1965),18,p.17-23 :
c
c       nedg
c df1 = sum ( gedg(k) * Parratt(X_k,n(k)) ) ,
c       k=1                                                   +------------------------+
c                                                             |ATTENTION: my paper with|
c       nedg                                                  |Olga Lugovskaya contains|
c df2 = sum ( 0.5*pi * (n(k)-1) * gedg(k) / X_k**(n(k)-1) ) . |a typo: instead of      |
c       k=1,                                                  |X_k=w_i/wedg it reads:  |
c    (X_k > 1)                                                |X_k=lambda_i/lambda_k   |
c                                                             +------------------------+
c where n(k) -       11/4 for (1s 1/2) electrons,
c                     7/3 for (2s 1/2) electrons,
c                     5/2 for all other electrons;
c Parratt(X_k,n(k)) - integral calculated by Parratt in:
c                     L.G.Parratt, C.F.Hempstead - Phys.Rev.
c                     (1954),94, p.1593-1600.
c gedg(k)  -          oscillator strengths for respective
c                     absorption edges.
c X_k = w_i/wedg(k),
c w_i      -          energy of incident x-ray wave in KeV:
c                     w_i[kev] = wave2energy(wave[a])
c wedg(nedg) -        energies of absorption edges in KeV.
c
c The absorption edges can be taken from Table-7 of X-ray Spectral
c Tables by M.A.Blokhin and I.G.Schweitzer (Rentgenospektralnyj
c spravochnik, Moscow, "Nauka", 1982) or elsewhere.
c Data for some hard-to-find absorption edges can be found in:
c D.T.Cromer - Acta Cryst.(1965),18,p.17-23.
c+---------------------------------------------------+
c|ATTENTION: if using Cromer's data, it should be    |
c|converted from the Rydberg units into KeV:         |
c|w[kev] = 1.09737312e-03 * 12.3981 * w[Rydberg] =   |
c|       = 13.60534e-03 * w[Rydberg]                 |
c+---------------------------------------------------+
c The oscillator strengths for many atoms and absorption edges
c can be found in the D.T.Cromer paper. For the other cases
c the oscillator strengths are interpolated in Disint using
c tabulated df1,df2 for several fixed wavelengths.
c -------------------------------------------------------
c New of August-99: for atoms with none or unknown absorption edges
c the above formula naturally do not work. Usually this for light
c elements like Li, H,.. The two other interpolation methods are
c used:
c 1. First we try to determine the coefficients in the formula
c    df = A * wave**B  with the help of substitution:
c               df(wave_1) = A ** (wave_1)**B
c               df(wave_2) = A ** (wave_2)**B
c 2. If the above does not work, for example when df(wave_1)
c    and df(wave_2) have opposite signs or if some of them is zero,
c    then a linear interpolation is applied.
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c+==============================================+
        Include         'x0h_size.inc'         !|
c In particular, the above include              |
c contains this:                                |
c       Parameter       (maxX0hWarnLines = 10); |
c+==============================================+
c Tabulated dispersion corrections for several wavelengths:
        Integer         idis, jHenkeCowan, nval, ifail,
     *                  k, k_, i, j, j1, j2, linpo, iname
c+-----------------------------------------+
c|This is how much the methods may differ: |
        Real            df_err /0.5/      !|
c+-----------------------------------------+

        Real            dfv1(nval), dfv2(nval), lfix(nval),
     *                  atom_mass_aem, w_i, gedg_sum, B,
     *                  td_disp, dmur_disp, n1,
     *                  df1_simple, df2_simple,
     *                  df1_diff, df2_diff

        Character DBtext(10)*12 /
     *        'AutoDB',         !-1 (1)
     *        'X0h',            !0  (2)
     *        'Henke/df1',      !1  (3)
     *        'Henke',          !2  (4)
     *        'Cowan/df1',      !3  (5)
     *        'Cowan',          !4  (6)
     *        'Windt/df1',      !5  (7)
     *        'Windt',          !6  (8)
     *        'Chantler/df1',   !7  (9)
     *        'Chantler'/       !8  (10)

        Character       name(kcompMax)*4
        Common  /x0pa2/ name

        Real            wave,df1(kcompMax),df2(kcompMax),
     *                  ab(kcompMax,nabcMax),s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,
     *                  Planck_,Boltzman,Avogadro,
     *                  wedg(nedgMax),gedg(nedgMax)
        Integer         indx,ipr,nedg
        Common  /x0pa4/ wave,df1,df2,ab,s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,indx,ipr,
     *                  Planck_,Boltzman,Avogadro,wedg,gedg,nedg

        Integer         linesX0hWarn
        Character       txtX0hWarn(maxX0hWarnLines)*80
        Common /x0hwarn/ txtX0hWarn, linesX0hWarn

        Integer         luwa,luat,lusp,luco
        Common  /lunds/ luwa,luat,lusp,luco

        Character       txt(20)*80
        Common  /msg/   txt

        Real*8          Parratt, X_k
        External        Parratt

        Real            wave2energy
        External        wave2energy

        ifail = 0
c Prevents GNU fortran -Wextra warnings:
        i = kcompMax
        i = kolMax
        i = ncodMax
        i = natmMax
        i = nwavMax
        i = nedgMax
        i = nabcMax
        i = nedgTb
        i = maxX0hWarnLines

        iname = max(Len_Trim(name(indx)),1)

c Incident X-ray wave energy in KeV:
        w_i = wave2energy(wave)

        df1(indx) = 0.
        df2(indx) = 0.

        df1_simple = 0.                 !df' data from Henke/Cowan/...
        df2_simple = 0.                 !df" data from Henke/Cowan/...

c In the idis=1 mode we do not have tabulated df',df" values:
        if (idis.ne.1 .or. jHenkeCowan.ne.0) Then !--------------------+
c +-------------------------------------------------------------+      |
c |####### First do a simple non-X0h style interpolation #######|      |2008/05/27
c +-------------------------------------------------------------+      |
c Find 2 closest wavelengths to wave:                                  |
          if (nval .lt. 2) stop 'Disper: nval < 2'                    !|
c The fixed wavelengths are normally sorted in the decreasing order    |
c (i.e. in the order of increasing energy). We want a fork, i.e. to    |
c find two fixed wavelengths above and below the current one.          |
          if (lfix(1).gt.lfix(nval)) Then !--+                         |
c lfix is in decreasing order                |                         |
            do j1=nval-1,1,-1  !==========+  |                         |
              if (lfix(j1).gt.wave) exit !|  |                         |
            enddo !=======================+  |                         |
            if (j1.lt.1) j1=1               !|                         |
            j2 = j1+1                       !|                         |
          else  !----------------------------+                         |
c lfix is in increasing order                |                         |
            do j2=2,nval  !===============+  |                         |
              if (lfix(j2).gt.wave) exit !|  |                         |
            enddo !=======================+  |                         |
            if (j2.gt.nval) j2=nval         !|                         |
            j1 = j2-1                       !|                         |
          endif !----------------------------+                         |
c Check if the wavelength is outside the interpolation range:          |
          if ((wave-lfix(j1))*(wave-lfix(j2)).gt.0) Then !------+      |
            linpo = 1  !yes, outside: use linear interpolation! |      |
          else  !-----------------------------------------------+      |
            linpo = 0  !no, in range: use x^n interpolation!    |      |
          endif  !----------------------------------------------+      |
                                                                      !|
c Find two closest waves, may be from the same side                    |
c (used before 2023.12, caused interpolation problems @ edges):        |
c         if (Abs(wave-lfix(1)) .lt. Abs(wave-lfix(2))) Then !-+       |
c           j1 = 1                                            !|       |
c           j2 = 2                                            !|       |
c         else !-----------------------------------------------+       |
c           j1 = 2                                            !|       |
c           j2 = 1                                            !|       |
c         endif  !---------------------------------------------+       |
c         do j=3,nval  !==========================+                    |
c           if (Abs(wave-lfix(j)).lt.            !|                    |
c    *          Abs(wave-lfix(j1))) Then !------+ |                    |
c             j2 = j1                          !| |                    |
c             j1 = j                           !| |                    |
c           else  !-----------------------------+ |                    |
c             if (Abs(wave-lfix(j)).lt.        !| |                    |
c    *            Abs(wave-lfix(j2))) Then !-+  | |                    |
c               j2 = j                      !|  | |                    |
c             endif  !-----------------------+  | |                    |
c           endif  !----------------------------+ |                    |
c         enddo  !================================+                    |
                                                                      !|
c         if (j1 .eq. j2) stop 'Disper: j1=j2'                        !|
c Check if the wavelength is outside the interpolation range:          |
c         if (Abs(wave-lfix(j1)).gt.Abs(lfix(j2)-lfix(j1))) Then !-+   |
c           linpo = 1  !yes, outside: use linear interpolation!    |   |
c         else  !--------------------------------------------------+   |
c           linpo = 0  !no, in range: use x^n interpolation!       |   |
c         endif  !-------------------------------------------------+   |
                                                                      !|
          if (ipr.eq.2) Then  !-------------------------------------+  |
            write (txt,35)    name(indx)(1:iname),                 !|  |
     *                        jHenkeCowan, wave,                   !|  |
     *                        nval, lfix(1), lfix(nval),           !|  |
     *                        j1, lfix(j1), dfv1(j1), dfv2(j1),    !|  |
     *                        j2, lfix(j2), dfv1(j2), dfv2(j2)     !v  v
  35        format('Element: ',a,' -- DISPER interpolation:'/
     *      'Henke/Cowan flag=',i2,'  wave=',g12.6/
     *      'nval=',i3,' wave(beg)=',g12.6,' wave(end)=',g12.6/
     *      'i1=',i3,' wave1=',g12.6,' f1''=',g10.4,' f1"=',g10.4/
     *      'i2=',i3,' wave2=',g12.6,' f1''=',g10.4,' f1"=',g10.4) !^  ^
            Call  Priblo  (txt,5,1,0,0,i)                          !|  |
          endif  !--------------------------------------------------+  |
                                                                      !|
c df1:                                                                 |
          if (dfv1(j1)*dfv1(j2).gt.0. .AND. linpo.eq.0) Then  !-----+  |
c Interpolation mode: df(wave) = A * (wave)**B                      |  |
c f1 = A * w1^B                                                     |  |
c f2 = A * w2^B                                                     |  |
c f1/f2 = (w1/w2)^B                                                 |  |
c ln(f1/f2) = B * ln(w1/w2)                                         |  |
c B = ln(f1/f2)/ln(w1/ws2)                                          |  |
c f = A * w^B                                                       |  |
c f/f1 = (w/w1)^B                                                   |  |
c f = f1 * (w/w1)^B                                                 |  |
c f = f1 * exp(ln((w/w1)^B))                                        |  |
c f = f1 * exp(B*ln(w/w1))                                          |  |
            B = alog(dfv1(j1)/dfv1(j2))                            !|  |
     /        / alog(lfix(j1)/lfix(j2))                            !|  |
c The B-formula does not apply to absorption edges where            |  |
c dispersion  corrections may not grow with wavelength              |  |
            if (B. gt. 0) Then !----------------------------------+ |  |
              df1_simple = dfv1(j1) * exp(B*alog(wave/lfix(j1))) !| |  |
              if (ipr.eq.2) Then  !--------------------------+    | |  |
                write (txt,36)    "'",B,"'",df1_simple      !|    | |  |
  36            format('B(f',a,')=',g12.6,' f',a,'=',g10.4) !|    | |  |
                Call  Priblo  (txt,1,1,0,0,i)               !|    | |  |
              endif  !---------------------------------------+    | |  |
            endif !-----------------------------------------------+ |  |
          endif  !--------------------------------------------------+  |
c If the above method did not apply (did not assign value):            |
          if (abs(df1_simple).lt.1E-37) Then !-------------------+     |
c Interpolation mode: df(wave) = A + B*wave                      |     |
            df1_simple = dfv1(j1) + (dfv1(j2)-dfv1(j1))         !|     |
     *                 * (wave-lfix(j1)) / (lfix(j2)-lfix(j1))  !|     |
          endif  !-----------------------------------------------+     |
c df2:                                                                 |
          if (dfv2(j1)*dfv2(j2).gt.0. .AND. linpo.eq.0) Then  !-----+  |
c Interpolation mode: df(wave) = A * (wave)**B                      |  |
c B = ln(f1/f2)/ln(w1/ws2)                                          |  |
c f = f1 * exp(B*ln(w/w1))                                          |  |
            B = alog(dfv2(j1)/dfv2(j2))                            !|  |
     /        / alog(lfix(j1)/lfix(j2))                            !|  |
c The B-formula does not apply to absorption edges where            |  |
c dispersion  corrections may not grow with wavelength              |  |
            if (B. gt. 0) Then !----------------------------------+ |  |
              df2_simple = dfv2(j1) * exp(B*alog(wave/lfix(j1))) !| |  |
              if (ipr.eq.2) Then  !--------------------------+    | |  |
                write (txt,36)    '"',B,'"',df2_simple      !|    | |  |
                Call  Priblo  (txt,1,1,0,0,i)               !|    | |  |
              endif  !---------------------------------------+    | |  |
            endif !-----------------------------------------------+ |  |
          endif  !--------------------------------------------------+  |
c If the above method did not apply or led to negative df2:            |
          if (df2_simple.le.1E-37) Then !------------------------+     |
c Interpolation mode: df(wave) = A + B*wave                      |     |
            df2_simple = dfv2(j1) + (dfv2(j2)-dfv2(j1))         !|     |
     *                 * (wave-lfix(j1)) / (lfix(j2)-lfix(j1))  !|     |
          endif  !-----------------------------------------------+     |
        endif !--------------------------------------------------------+

c+-----------------------------------------------------------------+
c|######## Now do a X0h style (Cromer-style) interpolation ########|
c+-----------------------------------------------------------------+
c Here we check if there are ANY oscillator strengths calculated:
        gedg_sum  = 0.
        do k=1,nedg  !=========================+
          gedg_sum = gedg_sum + abs(gedg(k))  !|
        enddo  !===============================+

        if (nedg.gt.0 .AND. abs(gedg_sum).gt.1.E-37) Then  !-----------+2008/05/27: changed from 0 to 1
                                                                      !|2008/06/03: changed back from 1 to 0
          do      k=1,nedg  !=======================+                  |
            X_k = dble(w_i/wedg(k))                !|                  |
            k_  = k                                !|                  |
            df1(indx) = sngl (df1(indx)            !|                  |
     +                + gedg(k) * Parratt(X_k,k_)) !|                  |
            if (X_k.gt.1.) Then  !--------------+   |The condition     |
              if     (k.eq.1) Then !---+        |   |X_k > 1 follows   |
                n1 =11./4.-1.         !|       !|   |from the Cromer   |
              elseif (k.eq.2) Then !---+        |   |paper Eq.(6):     |
                n1 = 7./3.-1.         !|       !|   |  w_i > wedg(k)   |
              else !-------------------+        |   |                  |
                n1 = 5./2.-1.         !|       !|   |                  |
              endif !------------------+        |   |                  |
              df2(indx) = sngl (df2(indx)      !|   |                  |
     +                  + 0.5*pi*n1*gedg(k)    !|   |                  |
     *                  * dexp(-n1*dlog(X_k))) !|   | == /X_k^{n1}     |
            endif  !----------------------------+   |                  |
          enddo  !==================================+                  |
                                                                      !|
          if ((idis.ne.1 .OR. jHenkeCowan.ne.0) .AND.                 !|
     *        abs(df1_simple*df2_simple).gt.1.E-37) Then !-----------+ |both df1 & df2 are non-zero
            df1_diff = abs((df1(indx)-df1_simple)/df1_simple)       !| |
            df2_diff = abs((df2(indx)-df2_simple)/df2_simple)       !| |
            if (ipr.eq.2) Then  !------------------------------+     | |
              write (txt,33)    name(indx)(1:iname),          !|     | |
     *                          df1(indx),df1_simple,         !|     | |
     *                          100.*df1_diff,                !|     | |
     *                          df2(indx),df2_simple,         !|     | |
     *                          100.*df2_diff,                !|     | |
     *                          idis, jHenkeCowan             !v     v v
  33          format('Element: ',a,' -- DISPER data:'/
     *        ' f''(Cramer)=',g10.4,'  f''(interpol)=',g10.4,
     *                                   '  diff.=',g10.4,'%'/
     *        ' f"(Cramer)=',g10.4,'  f"(interpol)=',g10.4,
     *                                   '  diff.=',g10.4,'%'/
     *        ' Dispersion mode =',i2,' Henke/Cowan flag=',i2)!^     ^ ^
              Call  Priblo  (txt,4,1,0,0,i)                   !|     | |
            endif  !-------------------------------------------+     | |
            if (df1_diff .gt. df_err .OR.                           !| |
     *          df2_diff .gt. df_err .OR.                           !| |
     *          df2(indx) .le. 0.) Then  !----------------------+    | |
              write (txt,200) name(indx)(1:iname)              !v    v v
  200         format('DISPER: oscillator strengths method',
     *        1x,'failed for [',a,']. Used interpolation.')    !^    ^ ^
              if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+  |    | |
                linesX0hWarn = linesX0hWarn+1               !|  |    | |
                txtX0hWarn(linesX0hWarn) = txt(1)           !|  |    | |
              endif !----------------------------------------+  |    | |
              if (ipr.eq.2) Call Priblo (txt,1,1,0,0,i)        !|    | |
              df1(indx) = df1_simple                           !|    | |
              df2(indx) = df2_simple                           !|    | |
            endif !---------------------------------------------+    | |
          endif !----------------------------------------------------+ |
                                                                      !|
        else   !-----(nedg.eq.0 .OR. gedg_sum.eq.0.)-------------------+
                                                                      !|
          write (txt,201) name(indx)(1:iname)                         !v
  201     format('DISPER: oscillator strengths method',
     *    1x,'not used for [',a,']: too few edges.')                  !^
          if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+             |
            linesX0hWarn = linesX0hWarn+1               !|             |
            txtX0hWarn(linesX0hWarn) = txt(1)           !|             |
          endif !----------------------------------------+             |
          if (ipr.eq.2) Call Priblo (txt,1,1,0,0,i)                   !|
          df1(indx) = df1_simple                                      !|
          df2(indx) = df2_simple                                      !|
                                                                      !|
        endif  !-------------------------------------------------------+

c Extra check: absorption cannot be negative:
        if (df2(indx) .le. 0.) Then  !----------------------------+
          write (txt,202) name(indx)(1:iname), wave, df2(indx)   !v
  202     format('DISPER: failed for atom [',a,'] @ wave=',g12.6,
     *                                ' with f"=',g12.5,' <= 0.'/
     *                                  'Setting: f'' = f" = 0.'/
     *    'Please try different Henke/Cowan/.. DB option.')      !^
          ifail = 3                                              !|
          if (jHenkeCowan.ge.-1 .AND.                            !|
     *        jHenkeCowan.le.8) Then !-------------------------+  |
            ifail = ifail + 1                                 !|  |
            i = max (Len_Trim(DBtext(jHenkeCowan+2)),1)       !|  |
            write (txt(ifail),203) DBtext(jHenkeCowan+2)(1:i),!|  |
     *                             DBtext(jHenkeCowan+2)(1:i),!|  |
     *                             df2_simple                 !|  |
  203       format ('Current DB selection: ',a,
     *              ', f"(',a,') = ',g12.5)                   !|  |
          endif !----------------------------------------------+  |
          if (df2_simple .lt. 0.) Then !-------------------+      |
            ifail = ifail + 1                             !|      |
            write (txt(ifail),204)                        !|      |
  204       format ('The DB indicates absorption',1x,
     *              'edge at or close to given energy.')  !|      |
          endif !------------------------------------------+      |
          do j=1,ifail !====================================+     |
            if (linesX0hWarn .lt. maxX0hWarnLines) Then !-+ |     |
              linesX0hWarn = linesX0hWarn+1              !| |     |
              txtX0hWarn(linesX0hWarn) = txt(j)          !| |     |
            endif !---------------------------------------+ |     |
          enddo !===========================================+     |
          if (ipr.eq.2) Call Priblo (txt,ifail,1,0,0,i)          !|
          ifail = 0     !clear the error to continue X0h          |
          df1(indx) = 0.                                         !|
          df2(indx) = 0.                                         !|
        endif !---------------------------------------------------+

c Print debug info:
        if (ipr.eq.2) Then  !--------------------------------+
                                                            !|
          td_disp = 2.*wave*e_radius*df2(indx)              !|
                                                            !|
          if (atom_mass_aem .gt. 0.) Then !----------+       |
            dmur_disp = td_disp * Avogadro          !|       |
     /                / ((1.e+16)*atom_mass_aem)    !|       |
          else  !------------------------------------+       |
            dmur_disp = 0.                          !|       |
          endif  !-----------------------------------+       |
                                                            !|
          write (txt,34)        name(indx)(1:iname),        !|
     *                          df1(indx),df2(indx),        !|
     *                          td_disp,                    !|
     *                          dmur_disp,                  !|
     *                          idis, jHenkeCowan           !v
  34      format('Element: ',a,' -- DISPER final results:'/
     *     5x,'Dispersion corrections:  f''=',g13.6/
     *     30x,                        'f"=',g13.6/
     *     4x,'Respective absorption section=',g13.6,' A**2'/
     *     8x,    'Respective mu/rho [cm2/g]=',g13.6/
     *    'Dispersion mode =',i2,' Henke/Cowan flag=',i2)   !^
          Call  Priblo  (txt,6,1,0,0,i)                     !|
        endif  !---------------------------------------------+
        return
        end
c
c =======================================================
c
        Real*8 Function Parratt (X_k,k)
c -------------------------------------------------------
c Calculate the Parratt integral for Disper according to:
c L.G.Parratt, C.F.Hempstead - Phys.Rev.(1954), 94, p.1593-1600.
c -------------------------------------------------------
        Real*8  X_k, z, a, b, c1, c2, c3, c4, PP
        Real    pi, s2, s3
        Integer k

        a  = 0.                                 !make GNU compiler happy
        b  = 0.                                 !make GNU compiler happy
        c1 = 0.                                 !make GNU compiler happy
        c2 = 0.                                 !make GNU compiler happy
        c3 = 0.                                 !make GNU compiler happy
        c4 = 0.                                 !make GNU compiler happy
        PP = 0.                                 !make GNU compiler happy
c
c This function does not work exactly at the edge:
        if (dabs(X_k-1.).lt.1.D-6) Then !----+
          Parratt = 0.                      !|
          return                            !|
        endif  !-----------------------------+

        pi = 4.*atan(1.)
        s2 = sqrt(2.)
        s3 = sqrt(3.)
        z  = 1./(X_k*X_k)

c All the Parratt (and df_r) dependencies on X_k occur
c to be same as for df_i.

        if     (k.eq.1) Then  !----------------------------+
c .....................................                    |
c k=1, 1s 1/2 electrons, n(k)=11/4, m=8                    |
c .....................................                    |
          if     (z.lt.1.) Then  !----------------------+  |
            a  = dexp(dlog(z)/8.)          !a=X_k^(-1/4)|  |
            c1 = dlog((1.+2.*a+a*a)/(1.-2.*a+a*a))/2.  !|  |
            c2 = dlog((1.+s2*a+a*a)/(1.-s2*a+a*a))/s2  !|  |
            c3 = 2.*datan(a)+s2*datan(s2*a/(1.-a*a))   !|  |
            c4 = pi*cos(pi/8.)/sin(pi/8.)              !|  |
            PP = -0.875*(c1+c2+c3-c4)*(a**7)           !|  |Parratt ~
          elseif (z.gt.1.) Then  !----------------------+  |  X_k^(-7/4)
            b  = dexp(-dlog(z)/8.)         !a=X_k^(1/4) |  |
            c1 = dlog((1.+2.*b+b*b)/(1.-2.*b+b*b))/2.  !|  |
            c2 = dlog((1.+s2*b+b*b)/(1.-s2*b+b*b))/s2  !|  |
            c3 = 2.*datan(b)+s2*datan(s2*b/(1.-b*b))   !|  |
            PP = -0.875*(c1+c2-c3)/(b**7)              !|  |
c         else  !---------------------------------------+  |
c           PP    = 0.                                 !|  |
          endif  !--------------------------------------+  |
                                                          !|
        elseif (k.eq.2) Then  !----------------------------+
c .....................................                    |
c k=2, 2s 1/2 electrons, n(k)=7/3, m=3                     |
c .....................................                    |
          if     (z.le.1.) Then  !----------------------+  |
            a  = dexp(dlog(z)/3.)          !a=X_k^(-2/3)|  |
            c1 = dlog((1.+a+a*a)/(1.-2.*a+a*a))/2.     !|  |
            c2 = s3*datan(s3*a/(2.+a))-pi/s3           !|  |
            PP = -(2./3.)*(c1+c2)*(a**2)               !|  |Parratt ~
          elseif (z.gt.1.) Then  !----------------------+  |  X_k^(-4/3)
            b  = dexp(-dlog(z)/3.)         !a=X_k^(2/3) |  |
            c1 = dlog((1.+b+b*b)/(1.-2.*b+b*b))/2.     !|  |
            c2 = s3*datan(s3*b/(2.+b))                 !|  |
            PP = -(2./3.)*(c1-c2)/(b**2)               !|  |
c         else  !---------------------------------------+  |
c           PP    = 0.                                 !|  |
          endif  !--------------------------------------+  |
                                                          !|
        else !---------------------------------------------+
c .....................................                    |
c k>2, all other electrons, n(k)=5/2, m=4                  |
c .....................................                    |
          if     (z.le.1.) Then  !----------------------+  |
            a  = dexp(dlog(z)/4.)          !a=X_k^(-1/2)|  |
            c1 = dlog((1.+2.*a+a*a)/(1.-2.*a+a*a))/2.  !|  |
            c2 = 2.*datan(a)-pi                        !|  |
            PP = -0.75*(c1+c2)*(a**3)                  !|  |Parratt ~
          elseif (z.gt.1.) Then  !----------------------+  |  X_k^(-3/2)
            b  = dexp(-dlog(z)/4.)         !a=X_k^(1/2) |  |
            c1 = dlog((1.+2.*b+b*b)/(1.-2.*b+b*b))/2.  !|  |
            c2 = 2.*datan(b)                           !|  |
            PP = -0.75*(c1-c2)/(b**3)                  !|  |
c         else  !---------------------------------------+  |
c           PP    = 0.                                 !|  |
          endif  !--------------------------------------+  |
                                                          !|
        endif  !-------------------------------------------+

        Parratt = PP
c       write (3,1) X_k, k, PP
c 1     format (' Parratt (X_k, k, PP): ', g14.6, i3, g14.6)
        return
        end
c
c =======================================================
c
        Subroutine CroSec (a,c,d,z,td,tq,isig,idb,atom_mass_aem,ifail)
c -------------------------------------------------------
c Calculate dipole and quadrupole cross sections components of X-ray
c absorption by crystal atoms with the help of one of these methods:
c 1.   Formula from [7] and electron shell screening constants from [8]
c 2-3. Interpolation formula  mu/ro=c*l**3-d*l**4
c 4.   Imaginary part of dispersion correction df2 from Henke/Cowan @ Redpar.
c .......................................................
c This is for medium energy range (5-25kev) and Z=6-54,
c as declared in the paper:
c  [7] H.Wagenfeld - Physical Review,(1966), v.144,
c      p.216-224.
c  [8] G.Hildebrandt, J.D.Stephenson, H.Wagenfeld -
c      Z.Naturforschung,(1975), v.30a, p.697-707.
c .......................................................
c    Subroutine input parameters:
c a     -    atom weight in atomic mass units,
c c,d   -    factors of the interpolation formula,
c z     -    number of electrons in the atom,
c td    -    dipole component of absorption cross section for atomic
c            element with number z and x-rays with wavelengths
c            l in Angstrom.
c tq    -    respective quadrupole component of the cross section.
c isig  -    index of calculation method (1-5).
c    =1 -    use td from Wagenfeld method
c    =2 -    use (mu/ro) = c*l^3 - d*l^4, c,d from Int.Tables
c    =3 -    use (mu/ro) = c*l^3 - d*l^4, c,d from mu/ro at 5 fixed waves
c    =4 -    use (mu/ro) = c*l^3 - d*l^4, c,d from mu/ro at 5 read waves
c    =5 -    use td = 2.*l*e_radius*df2(indx)
c idb   -    index of external database
c    =1 -    henke.dat
c    =2 -    cowan.dat or cowanlng.dat
c    =3 -    windt.dat
c    =4 -    chantler.dat
c atom_mass_aem - input
c ifail -    output
c -------------------------------------------------------
c The include file contains:
c       Parameter       (kcompMax = 10)
c       Parameter       (kolMax   = 32)
c       Parameter       (ncodMax  =200)
c       Parameter       (natmMax  =100)
c       Parameter       (nwavMax  = 50)
c       Parameter       (nedgMax  = 16)
c       Parameter       (nabcMax  = 11)
c       Parameter       (maxX0hWarnLines = 20);
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Integer         z, isig, idb, ifail,
     *                  jsig, jdb,
     *                  n1s, nful, n2p,
     *                  n3s, n3p, n3d, n4s, n4p,
     *                  i, j,
     *                  ndb, nval, jfail,
     *                  nedg_save, inddiff,
     *                  id_DB(5), iname

        Real            a, c, d, td, tq, f, q,
     *                  l, l1, l21, n2s, l31, l32, l41,
     *                  n1, n21, n31, n32, n41, lk,
     *                  tdn, tdm, tdl, tdk, tqk, tql,
     *                  t1sd, t1sq, t2sd, t2sq, t2pd, t2pq,
     *                  t3sd, t3dd, t3pd, t4sd, t4pd,
     *                  df2_cross, dmur_cross, td_f2

        Real            atom_mass_aem, lfix(nedgTb),
     *                  dfv1(nedgTb), dfv2(nedgTb),
     *                  df1_DB(5), df2_DB(5), td_DB(5),
     *                  wedg_save(nedgMax),
     *                  gedg_save(nedgMax),
     *                  diffmin, diffnow,
     *                  error_max /0.25/        !max allowed error for Wagenfeld

        Real            wave2energy
        External        wave2energy

        Character       DBopt(5)*8      /'None',
     *                                   'Henke DB',
     *                                   'Cowan DB',
     *                                   'Windt DB',
     *                                   'Chantler'/

        Character       name(kcompMax)*4
        Common  /x0pa2/ name

        Real            wave,df1(kcompMax),df2(kcompMax),
     *                  ab(kcompMax,nabcMax),s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,
     *                  Planck_,Boltzman,Avogadro,
     *                  wedg(nedgMax),gedg(nedgMax)
        Integer         indx,ipr,nedg
        Common  /x0pa4/ wave,df1,df2,ab,s1,s21,s31,s32,s41,
     *                  e_radius,Compton,Rydberg,pi,indx,ipr,
     *                  Planck_,Boltzman,Avogadro,wedg,gedg,nedg

        Integer         luwa,luat,lusp,luco
        Common  /lunds/ luwa,luat,lusp,luco

        Integer         linesX0hWarn
        Character       txtX0hWarn(maxX0hWarnLines)*80
        Common /x0hwarn/ txtX0hWarn, linesX0hWarn

        Integer                 iHenkeCowan
        Common  /HenkeCowan/    iHenkeCowan

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
          progname = 'x0h/CroSec'         !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+

        ifail      = 0                          !added 2023/12/01
        f          = 0.                         !make GNU compiler happy
        q          = 0.                         !make GNU compiler happy
        l          = 0.                         !make GNU compiler happy
        l1         = 0.                         !make GNU compiler happy
        lk         = 0.                         !make GNU compiler happy
        l21        = 0.                         !make GNU compiler happy
        l31        = 0.                         !make GNU compiler happy
        l32        = 0.                         !make GNU compiler happy
        l41        = 0.                         !make GNU compiler happy
        n1         = 0.                         !make GNU compiler happy
        n21        = 0.                         !make GNU compiler happy
        n31        = 0.                         !make GNU compiler happy
        n32        = 0.                         !make GNU compiler happy
        n41        = 0.                         !make GNU compiler happy
        n2s        = 0.                         !make GNU compiler happy
        tdn        = 0.                         !make GNU compiler happy
        tdm        = 0.                         !make GNU compiler happy
        tdl        = 0.                         !make GNU compiler happy
        tdk        = 0.                         !make GNU compiler happy
        tqk        = 0.                         !make GNU compiler happy
        tql        = 0.                         !make GNU compiler happy
        t1sd       = 0.                         !make GNU compiler happy
        t1sq       = 0.                         !make GNU compiler happy
        t2sd       = 0.                         !make GNU compiler happy
        t2sq       = 0.                         !make GNU compiler happy
        t2pd       = 0.                         !make GNU compiler happy
        t2pq       = 0.                         !make GNU compiler happy
        t3sd       = 0.                         !make GNU compiler happy
        t3dd       = 0.                         !make GNU compiler happy
        t3pd       = 0.                         !make GNU compiler happy
        t4sd       = 0.                         !make GNU compiler happy
        t4pd       = 0.                         !make GNU compiler happy
        td_f2      = 0.                         !make GNU compiler happy
        df2_cross  = 0.                         !make GNU compiler happy
        dmur_cross = 0.                         !make GNU compiler happy
c Prevents GNU fortran -Wextra warnings:
        i = kcompMax
        i = kolMax
        i = ncodMax
        i = natmMax
        i = nwavMax
        i = nedgMax
        i = nabcMax
        i = nedgTb
        i = maxX0hWarnLines
c -------------------------------------------------------
        jsig = isig                             !make local copy
        jdb = idb                               !make local copy
c -------------------------------------------------------
        inddiff = 0                             !initializing just in case
        if (jdb.lt.1 .OR. jdb.gt.4) jdb=0       !unknown
        iname = max (Len_Trim(name(indx)),1)
        l = wave

c  isig=1 - use td from Wagenfeld method
c  isig=2 - use (mu/ro) = c*l^3 - d*l^4, c,d from Int.Tables
c  isig=3 - use (mu/ro) = c*l^3 - d*l^4, c,d from mu/ro at 5 fixed waves
c  isig=4 - use (mu/ro) = c*l^3 - d*l^4, c,d from mu/ro at 5 read waves
c  isig=5 - use td = 2.*l*e_radius*df2(indx)
        goto (302,202,202,202,402)      jsig
c              1   2   3   4   5
c ================================
c METHOD 1: use td from Wagenfeld method
c Check if the parameters are available (passed):
  302   continue
        if (abs(s1).lt.1.E-37) Then  !-----------------------+
          txt(1) = 'CROSEC ['//name(indx)(1:iname)//        !|
     *             ']: unknown screening constants'         !|
          if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+   |
            linesX0hWarn = linesX0hWarn+1               !|   |
            txtX0hWarn(linesX0hWarn) = txt(1)           !|   |
          endif !----------------------------------------+   |
          goto 401 !------------------------------use-df2----+----+
        endif  !---------------------------------------------+    v
c Parameters are available:
        tdk = 0.
        tqk = 0.
        tdl = 0.
        tql = 0.
        tdm = 0.
        tdn = 0.
        if (z.lt.6 .or. z.gt.54) Then  !---------------------------+
          txt(1) = 'CROSEC ['//name(indx)(1:iname)//']: element ' !|
     *           //'is not in [6-54] range of Periodic Table'     !|
          if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+         |
            linesX0hWarn = linesX0hWarn+1               !|         |
            txtX0hWarn(linesX0hWarn) = txt(1)           !|         |
          endif !----------------------------------------+         |
          goto 401 !------------------------------use-df2----------+--+
        endif  !---------------------------------------------------+  v

c+============+
c| k - shell: |
c+============+
        n1s = 2
        if (z.eq.1) n1s=1
        nful = n1s
        l1   = 1./(Rydberg*(z-s1)**2)
        lk   = l1
c lk - tabulated wavelength of element's absorption k-edge. It differs
c from the wavelength l1 calculated according to [7] using hydrogen-like
c approximation: it should be lk>l1.
        if (wedg(1).gt.0.)  lk = wave2energy(wedg(1))
c If l<l1    the shell contributes,
c If lk>l>l1 the shell contributes, but the theory extension is applied,
c If l>lk    the shell does not contribute.
        if (l.le.lk)    Then  !--------------------------------+
          if (l.lt.l1)    Then  !--------------------------+   |
            n1   = sqrt(l/(l1-l))                         !|   |
            t1sd = ((2.**7)/3.)*pi*e_radius*l*((l/l1)**3) !|   |
     *           * exp(-4.*n1*atan2(1.,n1))               !|   |
     *           / (1.-exp(-2.*pi*n1))                    !|   |
          else  !------------------------------------------+   |
            n1   = sqrt(l/(l-l1))                         !|   |
            t1sd = ((2.**7)/3.)*pi*e_radius*l*((l/l1)**3) !|   |
     *           * exp(-2.*n1*alog((n1+1.)/(n1-1.)))      !|   |
          endif  !-----------------------------------------+   |
          t1sq = t1sd*(2./5.)*(Compton/l)*(4.-3.*l/l1)        !|
          t1sd = t1sd*(1.+(1./350.)*((Compton/l)**2)          !|
     *         * (59.-366.*l/l1+296.*(l/l1)**2                !|
     *         +  16.*(4.-3.*l/l1)*(9.-8.*l/l1))              !|
     *         +  (2./5.)*(Compton/l)*(1.-2.*l/l1))           !|
          tdk  = n1s*t1sd                                     !|
          tqk  = n1s*t1sq                                     !|
          if (z.le.nful) goto 100                             !|
        endif  !-----------------------------------------------+

c+============+
c| l - shell: |
c+============+
        n2s  = z-nful
        if (n2s.gt.2) n2s=2
        nful = nful + int (n2s+0.5)
        l21  = 4./(Rydberg*(z-s21)**2)
c Does the shell contribute?
        if (l21.ge.l) Then  !----------------------------------+
          n21  = 2.*sqrt(l/(l21-l))                           !|
          t2sd =((2.**10)/3.)*pi*e_radius*l*((l/l21)**3)      !|
     *         * (1.+3.*l/l21)                                !|
     *         * exp(-4.*n21*atan2(2.,n21))                   !|
     *         / (1.-exp(-2.*pi*n21))                         !|
          t2sq = t2sd*(8./5.)*(Compton/l)*(l/l21-1.)**2       !|
          t2sd = t2sd*(1.+(2./5.)*(Compton/l)*(1.-2.*l/l21))  !|
        endif  !-----------------------------------------------+
        if (z.le.nful) goto 100
        n2p = z-nful
        if (n2p.gt.6) n2p=6
c Does the shell contribute?
        if (l21.ge.l) Then  !----------------------------------+
          t2pd = ((2.**10)/3.)*pi*e_radius*l*((l/l21)**4)     !|
     *         * (1.+(8./3.)*l/l21)                           !|
     *         * exp(-4.*n21*atan2(2.,n21))                   !|
     *         / (1.-exp(-2.*pi*n21))                         !|
          t2pq = t2pd*(4./5.)*(Compton/l)*(11.-6.*l/l21)      !|
     *         * (1.+3.*l/l21)/(3.+8.*l/l21)                  !|
          t2pd = t2pd*(1.+(2./5.)*(Compton/l)*(1.-2.*l/l21))  !|
          tdl  = n2p*t2pd+n2s*t2sd                            !|
          tql  = n2p*t2pq+n2s*t2sq                            !|
        endif  !-----------------------------------------------+
        nful = nful+n2p
        if (z.le.nful) goto 100

c+============+
c| m - shell: |
c+============+
        n3s = z-nful
        if (n3s.gt.2) n3s=2
        l31 = 9./(Rydberg*(z-s31)**2)
c Does the shell contribute?
        if (l31.ge.l) Then  !------------------------------+
          n31  = 3.*sqrt(l/(l31-l))                       !|
          t3sd = (2.**7)*pi*e_radius*l*((l/l31)**3)       !|
     *         * (9.+96.*l/l31+208.*(l/l31)**2            !|
     *         + 128.*(l/l31)**3)                         !|
     *         * exp(-4.*n31*atan2(3.,n31))               !|
     *         / (1.-exp(-2.*pi*n31))                     !|
          tdm  = n3s*t3sd                                 !|
        endif  !-------------------------------------------+

        nful = nful+n3s
        if (z.le.nful) goto 100
        n3p  = z-nful
        if (n3p.gt.6) n3p=6
c Does the shell contribute?
        if (l21.ge.l) Then  !------------------------------+
          t3pd = (2.**10)*pi*e_radius*l*((l/l31)**4)      !|
     *         * (3.+26.*l/l31+28.*(l/l31)**2)            !|
     *         * exp(-4.*n31*atan2(3.,n31))               !|
     *         / (1.-exp(-2.*pi*n31))                     !|
          tdm  = n3p*t3pd+n3s*t3sd                        !|
        endif  !-------------------------------------------+
        nful = nful+n3p
        if (z.le.nful) goto 100
        n3d = z-nful
        if (n3d.gt.10) n3d=10

        l32 = 9./(Rydberg*(z-s32)**2)
c Does the shell contribute?
        if (l32.ge.l) Then  !------------------------------+
          n32  = 9.*sqrt(l/(l32-l))                       !|
          t3dd = ((2.**11)/5.)*pi*e_radius*l*((l/l32)**5) !|
     *         * (5.+46.*l/l32+48.*(l/l32)**2)            !|
     *         * exp(-4.*n32*atan2(3.,n32))               !|
     *         / (1.-exp(-2.*pi*n32))                     !|
          tdm  = n3d*t3dd+n3p*t3pd+n3s*t3sd               !|
        endif  !-------------------------------------------+
        nful = nful+n3d
        if (z.le.nful) goto 100

c+============+
c| n - shell: |
c+============+
        n4s = z-nful
        if (n4s.gt.2) n4s=2
        l41 = 16./(Rydberg*(z-s41)**2)
c Does the shell contribute?
        if (l41.ge.l) Then  !------------------------------+
          n41  = 4.*sqrt(l/(l41-l))                       !|
          t4sd = ((2.**13)/9.)*pi*e_radius*l*((l/l41)**3) !|
     *         * (3.+69.*l/l41+424.*(l/l41)**2            !|
     *         + 1024.*(l/l41)**3+(2944./3.)*(l/l41)**4   !|
     *         + 320.*(l/l41)**5)                         !|
     *         * exp(-4.*n41*atan2(4.,n41))               !|
     *         / (1.-exp(-2.*pi*n41))                     !|
          tdn  = n4s*t4sd                                 !|
        endif  !-------------------------------------------+
        nful = nful+n4s
        if (z.le.nful) goto 100
        n4p  = z-nful
        if (n4p.gt.6) n4p=6
c Does the shell contribute?
        if (l41.ge.l) Then  !------------------------------+
          t4pd = ((2.**13)/15.)*pi*e_radius*l*((l/l41)**4)!|
     *         * (75.+1400.*l/l41+5304.*(l/l41)**2        !|
     *         + 6016.*(l/l41)**3+2112.*(l/l41)**4)       !|
     *         * exp(-4.*n41*atan2(4.,n41))               !|
     *         / (1.-exp(-2.*pi*n41))                     !|
          tdn  = n4p*t4pd+n4s*t4sd                        !|
        else  !--------------------------------------------+
          txt(1)='CROSEC ['//name(indx)(1:iname)//        !|
     *           ']: too big wavelength for this method'  !|
          if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+ |
            linesX0hWarn = linesX0hWarn+1               !| |
            txtX0hWarn(linesX0hWarn) = txt(1)           !| |
          endif !----------------------------------------+ |
          goto 401  !------------------------use-df2-------+---+
        endif  !-------------------------------------------+   v

  100   continue
        td = tdn+tdm+tdl+tdk
        tq =         tql+tqk

c Testing for reasonable limits:
        td_f2 = 2.*l*e_radius*df2(indx)
        if (abs(td_f2) .gt. 1.E-37) Then  !--------------------------+
c If Wagenfeld differs too much (>25%) from dispersion correction,   |
c then use the dispersion correction for td:                         |
          diffmin = abs((td-td_f2)/td_f2)                           !|
          ndb = 0                                                   !|
          if (diffmin .gt. error_max) Then !-----------------------+ |
            ndb = ndb + 1                                         !| |
            df1_DB(ndb) = df1(indx)                               !| |
            df2_DB(ndb) = df2(indx)                               !| |
            td_DB(ndb)  = td_f2                                   !| |
            id_DB(ndb)  = jdb                                     !| |
            inddiff     = jdb                                     !| |
            if (ipr .eq. 2) Then  !----------------------------+   | |
              i = max(Len_Trim(DBopt(inddiff+1)),1)           !|   | |
              write (txt,201)  name(indx)(1:iname),           !|   | |
     *                         td/(2.*l*e_radius),            !|   | |
     *                         DBopt(inddiff+1)(1:i),         !|   | |
     *                         df2(indx),                     !|   | |
     *                         100.*diffmin                   !|   | |
  201         format('CROSEC [',a,']: Wagenfeld f"=',f7.4,
     *        ' differs from ',a,' f"=',f7.4,' by ',f6.0,'%') !|   | |
              Call Priblo (txt,1,1,0,0,i)                     !|   | |
            endif  !-------------------------------------------+   | |
                                                                  !| |
            if (jdb .eq. 0) Then !------------------------------+  | |
c Check again (do not trust dispersion                          |  | |
c corrections from atom.x0h):                                   |  | |
              do i=1,nedgMax !============+     !backup         |  | |
                wedg_save(i) = wedg(i)   !|     !backup         |  | |
                gedg_save(i) = gedg(i)   !|     !backup         |  | |
              enddo !=====================+     !backup         |  | |
              nedg_save = nedg                  !backup        !|  | |
              nval = nedgTb                                    !|  | |
c Calculate df2 according to Henke/Cowan/Windt/Chantler         |  | |
c (1=Henke, 2=Cowan, 3=Windt, 4=Chantler):                      |  | |
              do i=1,4 !======================================+ |  | |
                nval = nedgTb                                !| |  | |
                Call f1f2 (i, wave, name(indx), z, nval,     !| |  | |
     *                     lfix, dfv1, dfv2, ipr, jfail)     !| |  | |
                if (jfail .eq. 0) Then !--------------------+ | |  | |
                  Call Disint (dfv1,dfv2,lfix,jfail)       !| | |  | |
                  if (jfail. eq. 0) Then !-----------------+| | |  | |
                    nval = nedgTb                         !|| | |  | |
                    Call  Disper (atom_mass_aem, 3, i*2,  !|| | |  | |
     *                            dfv1, dfv2, lfix, nval, !|| | |  | |
     *                            jfail)                  !|| | |  | |
                    if (jfail.eq.0 .AND.                  !|| | |  | |
     *                  abs(df2(indx)).gt.1.E-37) Then !-+ || | |  | |
                      ndb = ndb+1                       !| || | |  | |
                      df1_DB(ndb) = df1(indx)           !| || | |  | |
                      df2_DB(ndb) = df2(indx)           !| || | |  | |
                      id_DB(ndb) = i     !DB index      !| || | |  | |
                    endif !------------------------------+ || | |  | |
                  endif !----------------------------------+| | |  | |
                endif !-------------------------------------+ | |  | |
              enddo !=========================================+ |  | |
              df1(indx) = df1_DB(1)             !restore        |  | |
              df2(indx) = df2_DB(1)             !restore        |  | |
              do i=1,nedgMax !============+     !restore        |  | |
                wedg(i) = wedg_save(i)   !|     !restore        |  | |
                gedg(i) = gedg_save(i)   !|     !restore        |  | |
              enddo !=====================+     !restore        |  | |
              nedg = nedg_save                  !restore        |  | |
              do i=1,ndb !============================+         |  | |
                td_DB(i) = 2.*l*e_radius*df2_DB(i)   !|         |  | |
              enddo !=================================+         |  | |
            endif !---------------------------------------------+  | |
            do i=1,ndb  !=============================+            | |
              diffnow = abs((td-td_DB(i))/td_DB(i))  !|            | |
              if (diffnow .lt. diffmin) Then !--+     |            | |
                diffmin = diffnow              !|     |            | |
                inddiff = id_DB(i)             !|     |            | |
                td_f2   = td_DB(i)             !|     |            | |
              endif ! --------------------------+     |            | |
            enddo !===================================+            | |
          endif  !-------------------------------------------------+ |
                                                                    !|
          if (diffmin .gt. error_max)  Then !----------------------+ |
c Still bad -- with other f" too:                                  | |
            i = max(Len_Trim(DBopt(inddiff+1)),1)                 !| |
            if (ipr .eq. 2) Then  !----------------------------+   | |
              write (txt,201)  name(indx)(1:iname),           !|   | |
     *                         td/(2.*l*e_radius),            !|   | |
     *                         DBopt(inddiff+1)(1:i),         !|   | |
     *                         td_f2/(2.*l*e_radius),         !|   | |
     *                         100.*diffmin                   !|   | |
              Call Priblo (txt,1,1,0,0,i)                     !|   | |
            endif  !-------------------------------------------+   | |
            write (txt(1),111,err=101) name(indx)(1:iname),       !| |
     *                                 DBopt(inddiff+1)(1:i),     !| |
     *                                 100.*diffmin               !| |
  111       format('CROSEC [',a,']: Wagenfeld f" differs from',
     *             1x,'f"(',a,') by',f6.0,'%; using DB f"')       !| |
            if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+       | |
              linesX0hWarn = linesX0hWarn+1               !|       | |
              txtX0hWarn(linesX0hWarn) = txt(1)           !|       | |
            endif !----------------------------------------+       | |
            if (ipr.eq.2) Then !-----------+                       | |
              Call Priblo (txt,3,1,0,0,i) !|                       | |
            endif  !-----------------------+                       | |
c Should we scale tq or set it zero?                              !| |
c Since Wagenfeld method failed, we should probably zero it.      !| |
cc          tq = tq * (td_df2/td)                                 !| |
            tq = 0.                                               !| |
            td = td_f2                                            !| |
          else   !-------------------------------------------------+ |
            if (ipr .eq. 2) Then  !----------------------------+   | |
              i = max(Len_Trim(DBopt(inddiff+1)),1)           !|   | |
              write (txt,301)  name(indx)(1:iname),           !|   | |
     *                         td/(2.*l*e_radius),            !|   | |
     *                         DBopt(inddiff+1)(1:i),         !|   | |
     *                         td_f2/(2.*l*e_radius),         !|   | |
     *                         100.*diffmin                   !|   | |
  301         format('CROSEC [',a,']: Wagenfeld f"=',f7.4,
     *        ' & ',a,' f"=',f7.4,' differ OK (',f6.0,'%)')   !|   | |
              Call Priblo (txt,1,1,0,0,i)                     !|   | |
            endif  !-------------------------------------------+   | |
          endif  !-------------------------------------------------+ |
        endif  !-----------------------------------------------------+
        goto 101

c+=====================================================+
c|Methods 2-4.                                         |
c|Calculate td using the formula from the International|
c|Tables for X-ray Crystallography:                    |
c|                                                     |
c|      (mu/ro) = c * l**3 - d * l**4 , ([cm**2/g])    |
c|          td  = (mu/ro)*(a/Avogadro) ,  ([cm**2])    |
c|                                                     |
c|Here mu/ro is the mass-absorption coefficient        |
c|Avogadro=6.025e+23 is the Avogadro constant,         |
c|a - the atom weight in atomic mass units,            |
c|c,d - taken for the Int.Tables (Method-2), or        |
c|    - interpolated from mu/ro for 5 known x-ray      |
c|      wavelengths (Method-3), or                     |
c|    - interpolated from mu/ro 5 specified x-ray      |
c|      wavelengths (Method-4).                        |
c|The interpolation is carried out in subroutine Redpar|
c|(see x0hread.for).                                   |
c+=====================================================+
c Check if the atomic weight is known:
  202   continue
        if (abs(a).lt.1.E-37) Then !-----------------------+
          txt(1) = 'CROSEC ['//name(indx)(1:iname)//      !|
     *             ']: unknown atomic weight'             !|
          if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+ |
            linesX0hWarn = linesX0hWarn+1               !| |
            txtX0hWarn(linesX0hWarn) = txt(1)           !| |
          endif !----------------------------------------+ |
          goto 401  !------------------------use-df2-------+---+
        endif  !-------------------------------------------+   v
        td = (c*l**3-d*l**4) * ((1.e+16)*a) / Avogadro
        tq = 0.
        goto 101

  401   continue
        txt(2) = ' '
        i = max(Len_Trim(DBopt(inddiff+1)),1)
        txt(3) = '-- Using f" of '//DBopt(inddiff+1)(1:i)//' for Xi0!'
        if (linesX0hWarn .lt. maxX0hWarnLines) Then !--+
          linesX0hWarn = linesX0hWarn+1               !|
          txtX0hWarn(linesX0hWarn) = txt(3)           !|
        endif !----------------------------------------+
        if (ipr.eq.2)   Then  !----------+
          Call  Priblo  (txt,3,1,0,0,i) !|
        endif  !-------------------------+
c ================================
c Method 5.
c Calculate td from the imaginary part of dispersion correction:
  402   continue
        jsig = 5
        td = 2.*l*e_radius*df2(indx)
        tq = 0.
c ================================
c Output data:
  101   continue
        f = td+tq
        if (f.gt.0.) Then !------+
          q = tq/f              !|
        else  !------------------+
          q = 0.                !|
        endif  !-----------------+
        if (ipr.eq.2)   Then  !------------------------+
                                                      !|
          df2_cross = td / (2.*l*e_radius)            !|
                                                      !|
          if (a.gt.0.) Then !--------------------+     |!a=atom_mass_aem!
            dmur_cross = td * Avogadro          !|     |
     /                 / ((1.e+16)*a)           !|     |
          else  !--------------------------------+     |
            dmur_cross = 0.                     !|     |
          endif  !-------------------------------+     |
                                                      !|
          write (txt,120)       name(indx)(1:iname),  !|
     *                          f,                    !|
     *                          q,                    !|
     *                          df2_cross,            !|
     *                          dmur_cross,           !|
     *                          isig, jsig            !V
  120     format('Element: ',a,' -- CROSEC data'/
     *      11x,       'Absorption section=',g13.6,' A**2'/
     *      13x,         'Quadrupole quote=',g13.6/
     *       4x,'Respective correction df2=',g13.6/
     *       4x,'Respective mu/rho [cm2/g]=',g13.6/
     *    'Absorption mode =',i1,' (',i1,')')         !^
          Call  Priblo  (txt,6,1,0,0,i)               !|
        endif  !---------------------------------------+
        if (abs(td).lt.1.E-32) Then !-------------------------------+
          write (txt,121)       name(indx)(1:iname), wave,         !|
     *                          isig, jsig, iHenkeCowan            !|
  121     format('CROSEC [',a,']: calculation failed (td=0) at:'/
     *    'wave=',g12.6,', isig=',i1,' (',i1,'),',1x,
     *                           'index(Henke/Cowan/..)=',i1,'.'/
     *    'Please try different Henke/Cowan/.. option.')           !|
          ifail = 3                                                !|
          do j=1,ifail !====================================+       |
            if (linesX0hWarn .lt. maxX0hWarnLines) Then !-+ |       |
              linesX0hWarn = linesX0hWarn+1              !| |       |
              txtX0hWarn(linesX0hWarn) = txt(j)          !| |       |
            endif !---------------------------------------+ |       |
          enddo !===========================================+       |
          if (ipr.eq.2) Call Priblo (txt,ifail,1,0,0,i)            !|
          ifail = 0     !clear the error to continue X0h            |
        endif !-----------------------------------------------------+
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
c ================================
        end
c
c =======================================================
c
        Real Function FiFunc (x)
c -------------------------------------------------------
c Function fifunc calculates integral over t in range (0,x) for the
c expression t/(exp(t)-1). It is used for the calculation of the
c Debye-Waller factor at an arbitrary crystal structure.
c The integral is evaluated the Simpson method.
c f(t)  - the integrand,
c a,b   - the integration limits,
c epsil - absolute precision.
c -------------------------------------------------------
        Integer i, n
        Real    x, f, t, a, b, epsil, s, snext,
     *          h, odd, even

        f(t)  = t/(exp(t)-1.0)
        a     = 0.0
        b     = x
        epsil = 0.001

c ***  Initialization;
c s,snext - integral sums at steps h and h/2.
        h     = (b-a)/2.0
        snext = 1.0e+38

c *** Calculation of integral sum snext:
  11    continue  !----<-------<----------------------------+
        s = snext                                          !|
        h = h/2.0                                          !|
        n = int(0.5*(b-a)/h + 0.1)                         !|changed 2013.08
                                                           !|
c odd: sum of the integrand values at points a+h,a+3*h,...  |
        odd = 0.0                                          !|
        do      i=1,n  !==============+                     |
          odd = odd + f(a+(2*i-1)*h) !|                     ^
        enddo  !======================+                     |
                                                           !|
c even: sum of the integrand values at points a+2*h,a+4*h,..|
        even = 0.0                                         !|
        n    = n-1                                         !|
        do  i=1,n  !=================+                      |
          even = even + f(a+2*i*h)  !|                      ^
        enddo  !=====================+                      |
                                                           !|
c To exclude the critical point we take f(a)=f(0.0)=1.0     |
        snext = (h/3.0*(1.0+f(b)+4.0*odd+2.0*even))/x      !|
                                                           !|
c *** The Runge rule:                                       |
        if      (abs(s-snext) .gt. epsil) goto 11  !--->----+

        Fifunc = snext

        return
        end
c
c =======================================================
c
        Subroutine Sort_ABC (name,list,n)

        Integer         n, list(n), i, j, k, lwrk
        Character       name(n)*(*), wa*32, wb*32, wc*32
c--------------------------------------------------------
c  Sort Names in ascending order with ignoring case
c--------------------------------------------------------
        k = Len(name(1))
        if (k .gt. 32)  Stop 'Sort_ABC: long name'
        do      j = 2,n  !============================+
          lwrk    = list(j)                          !|
          wa(1:k) = name(j)                          !|
          wb(1:k) = name(j)                          !|
          Call  Case (wb(1:k),0)                     !|
                                                     !|
          do    i = j-1,1,-1  !================+      |
            wc(1:k) = name(i)                 !|      |
            Call Case (wc(1:k),0)             !|      |
            if (wc(1:k) .le. wb(1:k)) goto 3 !-+->+   |
            list(i+1) = list(i)               !|  |   |
            name(i+1) = name(i)               !|  |   |
          enddo  !=============================+  V   |
          i=0                                    !|   |
  3       continue   !<---------------------------+   |
          list(i+1) = lwrk                           !|
          name(i+1) = wa(1:k)                        !|
        enddo  !======================================+
        return
        end
c
c =======================================================
c
        Subroutine Sort_Index (name,list,n)

        Integer         n, list(n), lwrk, i, j, k
        Character       name(n)*(*), w*32
c--------------------------------------------------------
c  Sort Names in descending order over Index
c--------------------------------------------------------
        k = Len(name(1))
        if (k. gt. 32)  Stop 'Sort_Index: long name'
        do      j = 2,n  !============================+
          lwrk   = list(j)                           !|
          w(1:k) = name(j)                           !|
                                                     !|
          do    i = j-1,1,-1  !================+      |
            if (list(i) .ge. lwrk) goto 3 !----+->+   |
            list(i+1) = list(i)               !|  |   |
            name(i+1) = name(i)               !|  |   |
          enddo  !=============================+  V   |
          i = 0                                  !|   |
  3       continue   !<---------------------------+   |
          list(i+1) = lwrk                           !|
          name(i+1) = w(1:k)                         !|
        enddo  !======================================+
        return
        end
