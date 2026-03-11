c This module includes the following subroutines:
c   1. Read_Top99_M0
c   2. Read_Sub99_M0
c -- which read the surface layer profile and the
c substrate data, respectively, for ter_sl97 and
c trds_sl97 programs.
c
c =======================================================
c
        Subroutine Read_Top99_M0 (File_Top,
     *                            N_Top, N_Top_Max,
     *                            Thickness, Thickness_Top,
     *                            Code, rho,
     *                            Frac, N_Frac, N_Frac_Max,
     *                            C00, W0,
     *                            x0, F10, F11, F1T,
     *                            Mdensity, Mshare, Mvect,
     *                            Sigma, Transit, YesTr,
     *                            Identical, ifail)
c--------------------------------------------------------
c       Reading of top layer profile -- smart version!
c
c                  Author: S.Stepanov
c -------------------------------------------------------
        Integer         N_Top_Max, N_Frac_Max

        Complex*8       x0(N_Top_Max),
     *                  F10(N_Top_Max),
     *                  F11(N_Top_Max),
     *                  F1T(N_Top_Max)

        Real*8          TTT

        Real*4          Mdensity(N_Top_Max),
     *                  Mshare(N_Top_Max),
     *                  Mvect(3,N_Top_Max),
     *                  Thickness(N_Top_Max),
     *                  Transit(N_Top_Max),     ! == Thickness_Tr
     *                  Sigma(N_Top_Max),
     *                  rho(N_Top_Max),
     *                  W0(N_Top_Max),
     *                  Frac(N_Frac_Max,N_Top_Max),
     *                  Thickness_Top, Frac_Sum

        Logical         YesTr,Mdensity_Specified

        Integer         N_Frac(N_Top_Max), N_Top,
     *                  Identical(N_Top_Max),
     *                  nbegin(2), nperiod(2),
     *                  lpn, iprf, lun, line, lena,
     *                  ncycle, ll, i, j, k, l, m,
     *                  l1, l2, ifail, lfile, ltx,
     *                  izero_frac, mzero_frac,
     *                  merr, t_found, nlines,
     *                  io_status

        Character       File_Top*(*), file*80,
     *                  Code(N_Frac_Max,N_Top_Max)*(*),
     *                  C00(N_Top_Max)*(*),
     *                  Fij*3

        Character       txt(20)*80
        Common  /msg/   txt
c -------------------------------------------------------
c Buffer space:
        Character       buff*4000
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv .lt. 10) Then  !----+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Read_Top99_M0'      !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim (progname)
c -------------------------------------------------------
        i  = 0                                  !keep GNU Fortran happy
        j  = 0                                  !keep GNU Fortran happy
        k  = 0                                  !keep GNU Fortran happy
        l  = 0                                  !keep GNU Fortran happy
        l1 = 0                                  !keep GNU Fortran happy
        l2 = 0                                  !keep GNU Fortran happy
        ll = 0                                  !keep GNU Fortran happy
        m  = 0                                  !keep GNU Fortran happy
        ifail = 0
        Thickness_Top = 0.
        TTT           = 0.
        iprf = Len_Trim(File_Top)
        if (iprf .eq. 0) goto 99  !--------return-------+
        file = File_Top                                !v
        file = file(1:iprf)//'.prf'
        lun = 4
        lfile = Len_Trim(file)

        Call OpenFile(file(1:lfile),lun,'read','old',io_status,*1)
        l       = 0
        line    = 0
        ncycle  = 0  !number of open loops
  5     continue    !<--------------------------------<-----+        <-----+
        t_found = 1    !no problems with layer thickness "t"|              |
        txt(19) = 'file corrupt'                           !|              |
        buff   = ' '                                       !|              |
        line = line + 1                                    !|              |
        Read (lun,'(a)',err=2,end=4) buff   !-------------->+--+           |
        ll = Len_Trim(buff)                                !^  |E          |
        if (ll .eq. 0)              goto 5  !->-------------+  |N  -empty--+
        if ((buff(1:1) .eq. ';') .or.                      !|  |D  comment-+
     *      (buff(1:1) .eq. '!'))  goto 5  !->--------------+  |   comment-+
                                                           !|  v
c Remove hidden/unprintable characters (including <TAB>):   |
        do j=1,ll !================================+        |
          if (buff(j:j) .lt. ' ') buff(j:j) = ' ' !|        |
        enddo  !===================================+        |
c Remove comments at the end of the line:                   |
        j = Index(buff(1:ll),';')                          !|
        if (j .gt. 0) Then !-+                              |
          buff(j:ll) = ' '  !|                              |
          ll = j-1          !|                              |
        endif  !-------------+                              |
        j = Index(buff(1:ll),'!')                          !|
        if (j .gt. 0) Then !-+                              |
          buff(j:ll) = ' '  !|                              |
          ll = j-1          !|                              |
        endif  !-------------+                              |
        if (ll .eq. 0)             goto 5   !->-------------+
c Compress the line:                                        |
c       Call Compress_Line (buff,ll,'  ',1)                !|
        Call ComprSp (buff,ll)                             !|
        Call Compress_Line (buff,ll,' =',0)                !|
        Call Compress_Line (buff,ll,'= ',1)                !|
                                                           !|
c I presume: ll<2000!                                       |
        buff(2001:4000) = buff(1:2000)                     !^
c Switch to lower case:                                     |
        Call Case (buff(1:ll),1)                           !|
c-----------------------------------------------------------+
c Search for the keywords:                                  |
c period= ...  end period                                   |
c t= code= rho= x0= w0= sigma= tr=  [x=                     ^
c                                    code2= x2=             |
c                                    code3= x3=             |
c                                    code4=  ]              |
c              mdensity= mshare= mvector= F10= F11= F1T=    |
c Here:                                                     |
c mshare - the share of magnetic atoms in the material [0:1]|
c mdensity - density of magnetic atoms (1/cm^3).            |
c      ---  Only one of mshare and mdensity can be specified|
c           at a time.                                      |
c mvector - x,y,z components of the magnetization vector in |
c           the material (arbitrary normalization accepted).|
c      ---  Here x=along beam on surface; z=along internal  |
c      ---  normal, y=in the direction perperndicular to the|
c      ---  scattering plane.                              !|
c F10, F11, F1T - magnetic scattering amplitudes (complex!) |
c      ---  { these are actually (3*wave/8*pi)*Fkl/ro }     |
c-----------------------------------------------------------+
                                                           !|
        Call ArgSearch (buff(1:ll),'period=',' ',',',k,i,j)!|
        if (k .gt. 0)  Then  !-------------------------+    |
          ncycle = ncycle+1                           !|    |
          if (ncycle .gt. 2)    goto 6  !--------------+----+-ERROR+
          Call rdInt (nperiod(ncycle),1,buff(i:j),m)  !|    |      v
          txt(19) = 'period'                          !|    |
          if (m .ne. 0)   goto 2  !--------------------+----+-ERROR+
          txt(19) = 'zero/negative period'            !|    |      v
          if (nperiod(ncycle) .lt. 1) goto 2 !---------+----+-ERROR+
          nbegin(ncycle) = l+1                        !|    |      v
          buff(k:j) = ' '                             !|    |
          Call txtShift (buff(1:ll),j,i)              !|    ^
c -- the call returns: i=Len_Trim(vm)                  |    |
          if (i.gt.0)  goto 19  !----------------------+----+-ERROR+
          goto 5   !->-------------------------->------+----+      v
        elseif (k .eq. -1)  Then  !--------------------+    |
          txt(19) = 'duplicate keyword "period="'     !|    |
          goto 2 !------------------------------->-----+----+-ERROR+
        endif  !---------------------------------------+    |      v
                                                           !|
                                                           !|
        i = Index(buff(1:ll),'end period')                 !|
        j = Index(buff(1:ll),'endperiod')                  !|
        if (i .gt. 0 .or. j .gt. 0) Then !-------------+    |
          if (i .gt. 0) Then !---+                     |    |
            buff(i:i+9) = ' '   !|                     |    |
          else  !----------------+                     |    |
            buff(j:j+8) = ' '   !|                     |    |
          endif  !---------------+                     |    |
          ncycle = ncycle-1                           !|    ^
          if (ncycle .lt. 0)    goto 7  !--------------+----+-ERROR+
          lena = l-nbegin(ncycle+1)+1                 !|    |      v
          l  = l + lena*(nperiod(ncycle+1)-1)         !|    |
          if (l .gt. N_Top_Max) goto 3   !-->----------+----+-ERROR+
          do    i=2,nperiod(ncycle+1)  !======+        |    ^      v
          do    j=1,lena   !================+ |        |    |
            l1 = nbegin(ncycle+1)+j-1      !| |        |    |
            l2 = l1 + lena*(i-1)           !| |        |    |
            Thickness(l2) = Thickness(l1)  !| |        |    ^
            TTT = TTT + Thickness(l1)      !| |        |    |
            N_Frac(l2)    = N_Frac(l1)     !| |        |    |
            do m=1,N_Frac_Max !=========+   | |        |    |
              Code(m,l2)  = Code(m,l1) !|   | |        |    |
              Frac(m,l2)  = Frac(m,l1) !|   | |        |    |
            enddo  !====================+   | |        |    |
            Identical(l2) = l1             !| |        |    |
            x0(l2)        = x0(l1)         !| |        |    |
            W0(l2)        = W0(l1)         !| |        |    |
            Sigma(l2)     = Sigma(l1)      !| |        |    |
            if (YesTr)                     !| |        |    |
     *      Transit(l2)   = Transit(l1)    !| |        |    |
            C00(l2)       = C00(l1)        !| |        |    |
            rho(l2)       = rho(l1)        !| |        |    |
            Mdensity(l2)  = Mdensity(l1)   !| |        |    |
            Mshare(l2)    = Mshare(l1)     !| |        |    |
            F10(l2)       = F10(l1)        !| |        |    |
            F11(l2)       = F11(l1)        !| |        |    |
            F1T(l2)       = F1T(l1)        !| |        |    |
            do m=1,3  !==================+  | |        |    |
              Mvect(m,l2) = Mvect(m,l1) !|  | |        |    |
            enddo  !=====================+  | |        |    |
          enddo  !==========================+ |        |    |
          enddo  !============================+        |    |
          Call txtShift (buff(1:ll),j,i)              !|    ^
c -- the call returns: i=Len_Trim(vm)                  |    |
          if (i.gt.0)  goto 19  !----------------------+----+-ERROR+
          goto 5   !->------------------------>--------+----+      v
        elseif (i .eq. -1  .or.  j .eq. -1) Then !-----+    |
          txt(19) = 'duplicate keyword "end period"'  !|    |
          goto 2 !------------------------------>------+----+-ERROR+
        endif  !---------------------------------------+    |      v
                                                           !|
                                                           !^
c Now, since the line does not describe a period,           |
c it should describe a layer:                               |
        l = l + 1                                          !|
        if (l .gt. N_Top_Max)   goto 3   !-->---------------+-ERROR+
        Thickness(l) = 0.                                  !|      v
        N_Frac(l)    = 1                                   !|
        do m=1,N_Frac_Max  !===+                            |
          Code(m,l)  = ' '    !|                            |
          Frac(m,l)  = 0.     !|                            |
        enddo  !===============+                            |
        Frac(1,l)    = 1.                                  !|
        Identical(l) = 0                                   !|
        x0(l)        = (0.,0.)                             !^
        W0(l)        = 1.                                  !|
        Sigma(l)     = 0.                                  !|
        if (YesTr)                                         !|
     *  Transit(l)   = 0.                                  !|
        C00(l)       = ' '                                 !|
        rho(l)       = 0.                                  !|
        Mdensity(l)  = 0.                                  !|
        Mshare(l)    = 0.                                  !|
        F10(l)       = (0., 0.)                            !|
        F11(l)       = (0., 0.)                            !|
        F1T(l)       = (0., 0.)                            !|
        do m=1,3  !========+                                |
          Mvect(m,l) = 0. !|                                |
        enddo  !===========+                                |
                                                           !|
        Call ArgSearch (buff(1:ll),'t=',' ',',',k,i,j)     !|
        if (k .gt. 0)  Then  !--------------------------+   |
          Call rdReal (Thickness(l),1,buff(i:j),m)     !|   |
          txt(19) = 'thickness'                        !|   |
          if (m .ne. 0)  goto 2 !-----------------------+---+-ERROR+
          buff(k:j) = ' '                              !|   |      v
          txt(19) = 'zero layer thickness "t" - '//    !|   |
     *              'please comment out the line!'     !|   |
          if (Abs(Thickness(l)) .lt. 1.E-32) goto 2 !---+->-+-ERROR+
          txt(19) = 'negative layer thickness "t"'     !|   |      v
          if (Thickness(l) .lt. 0.) goto 2 !----->------+---+-ERROR+
          TTT = TTT + Thickness(l)                     !|   |      v
          t_found = 1                                  !|   |
        elseif (k .eq. -1)  Then  !---------------------+   |
          txt(19) = 'duplicate keyword "t="'           !|   |
          goto 2 !------------------------------->------+---+-ERROR+
        elseif (k .eq. 0)   Then  !---------------------+   |      v
          t_found = 0  !unspecified layer thickness "t" |   |
        endif  !----------------------------------------+   |
                                                           !|
                                                           !|
        Call ArgSearch (buff(1:ll),'w0=',' ',',',k,i,j)    !|
        if (k .gt. 0)  Then  !-------------------------+    |
          Call rdReal (W0(l),1,buff(i:j),m)           !|    |
          txt(19) = 'w0'                              !|    |
          if (m .ne. 0)       goto 2 !-----------------+----+-ERROR+
          txt(19) = 'negative "w0"'                   !|    |      v
          if (W0(l) .lt. 0.)  goto 2 !-----------------+----+-ERROR+
          txt(19) = 'w0 > 99'                         !|    |      v
          if (W0(l) .gt. 99.) goto 2 !-----------------+----+-ERROR+
          buff(k:j) = ' '                             !|    |      v
        elseif (k .eq. -1)  Then  !--------------------+    |
          txt(19) = 'duplicate keyword "w0="'         !|    |
          goto 2 !------------------------------->-----+----+-ERROR+
        endif  !---------------------------------------+    |      v
                                                           !|
                                                           !|
        Call ArgSearch (buff(1:ll),'sigma=',' ',',',k,i,j) !|
        if (k .gt. 0)  Then  !-------------------------+    |
          Call rdReal (Sigma(l),1,buff(i:j),m)        !|    |
          txt(19) = 'sigma'                           !|    |
          if (m .ne. 0)                   goto 2 !-----+----+-ERROR+
          txt(19) = 'negative rms roughness "sigma"'  !|    |      v
          if (Sigma(l) .lt. 0.)           goto 2 !-----+----+-ERROR+
          txt(19) = 'rms roughness "sigma" exceeds'// !|    |      v
     *        ' layer thickness by more than 2 times' !|    |
          if (Sigma(l) .gt. 2.*Thickness(l) .AND.     !|    |
     *        t_found .ne. 0)             goto 2 !-----+----+-ERROR+
          buff(k:j) = ' '                             !|    |      v
        elseif (k .eq. -1)  Then  !--------------------+    |
          txt(19) = 'duplicate keyword "sigma="'      !|    |
          goto 2 !------------------------------->-----+----+-ERROR+
        endif  !---------------------------------------+    |      v
                                                           !|
                                                           !|
        if (YesTr) Then  !-------------------------------+  |
         Call ArgSearch(buff(1:ll),'tr=',' ',',',k,i,j) !|  |
          if (k .gt. 0)  Then  !----------------------+  |  |
            Call rdReal (Transit(l),1,buff(i:j),m)   !|  |  |
            txt(19) = 'tr'                           !|  |  |
            if (m .ne. 0)                   goto 2 !--+--+--+-ERROR+
            txt(19) = 'negative transition '//       !|  |  |      v
     *                       'layer thickness "tr"'  !|  |  |
            if (Transit(l) .lt. 0.)         goto 2 !--+--+--+-ERROR+
            txt(19) = 'transition domain "tr" '//    !|  |  |      v
     *                'exceeds layer thickness "t"'  !|  |  |
            if (Transit(l) .gt. Thickness(l) .AND.   !|  |  |
     *          t_found .ne. 0)             goto 2 !--+--+--+-ERROR+
            buff(k:j) = ' '                          !|  |  |      v
            if ((Transit(l) .gt. 0.) .AND.           !|  |  |
     *          (Sigma(l) .gt. 0.))         goto 61 !-+--+--+-ERROR+
          elseif (k .eq. -1)  Then  !-----------------+  |  |      v
            txt(19) = 'duplicate keyword "tr="'      !|  |  |
            goto 2 !----------------------------->----+--+--+-ERROR+
          endif  !------------------------------------+  |  |      v
        endif  !-----------------------------------------+  |
                                                           !^
        Call ChiSearch (buff(1:ll), 'x0=', x0(l),          !|
     *                  C00(l)(1:1), txt(19), '0,', *2)    !|
                                                           !|
        Frac_Sum  = 0.                                     !|
        izero_frac= 0                                      !|
        mzero_frac= 0                                      !|
        do m=1,N_Frac_Max  !===========================+    |
                                                      !|    |
          txt(18) = 'codeN='                          !|    |
          write(txt(18)(5:5),'(i1)') m                !|    |
          ltx = 6                                     !|    |
c Search for "codeN=":                                 |    |
          Call ArgSearch (buff(1:ll),txt(18)(1:ltx),  !|    |
     *                               ' ',',',k,i,j)   !|    |
          if (k .eq. 0 .AND. m .eq. 1) Then  !-----+   |    |
c Retry with "code=" instead of "code1="           |   |    |
            Call ArgSearch (buff(1:ll),'code=',   !|   |    |
     *                             ' ',',',k,i,j) !|   |    |
          endif  !---------------------------------+   |    |
          if (k .gt. 0) Then  !---------------------+  |    |
                                                   !|  |    ^
c One more fraction code found:                     |  |    |
            if (m .gt. 1) N_Frac(l) = N_Frac(l)+1  !|  |    |
                                                   !|  |    |
            if (buff(j:j) .eq. ',') Then !--+       |  |    |
              buff(j:j)           = ' '    !|       |  |    |
              buff(j+2000:j+2000) = ' '    !|       |  |    |
              if (j .gt. i) j = j-1        !|       |  |    |
            endif  !------------------------+       |  |    |
            Code(m,l) = buff(i+2000:j+2000)        !|  |    |
            i = Len_Trim(Code(m,l))                !|  |    |
            if ((Code(m,l)(1:1).eq.char(34) .and.  !|  |    |
     *           Code(m,l)(i:i).eq.char(34)) .or.  !|  |    ^
     *          (Code(m,l)(1:1).eq.char(39) .and.  !|  |    |
     *           Code(m,l)(i:i).eq.char(39)))      !|  |    |
     *      Code(m,l) = Code(m,l)(2:i-1)           !|  |    |
            buff(k:j) = ' '                        !|  |    |
                                                   !|  |    |
          elseif (k .eq. -1)  Then  !---------------+  |    |
                                                   !|  |    |
            write (txt(19),119) txt(18)(1:ltx)     !|  |    |
  119       format('duplicate keyword "',a,'"')    !|  |    |
            goto 2 !----------------------------->--+--+----+-ERROR+
                                                   !|  |    |
          else  !-----------------------------------+  |    ^
                                                   !|  |    |
c Goto exit from the loop:                         !|  |    |
            if (m .gt. 1)  goto 71 !----------------+--+--+ |
c (if Code-1 is missing, the substrate code is used.|  |  | |
c  Here we simply continue)                         |  |  | |
                                                   !|  |  | |
          endif  !----------------------------------+  |  | |
                                                      !|  | |
          txt(18) = 'xN='                             !|  | |
          write(txt(18)(2:2),'(i1)') m                !|  | |
          ltx = 3                                     !|  | |
c Search for "xN=":                                   !|  | |
          Call ArgSearch (buff(1:ll),txt(18)(1:ltx),  !|  | |
     *                               ' ',',',k,i,j)   !|  | |
          if (k .eq. 0 .AND. m .eq. 1) Then  !-----+   |  | |
c Retry with "x=" instead of "x1="                 |   |  | |
            Call ArgSearch (buff(1:ll),'x=',      !|   |  | |
     *                             ' ',',',k,i,j) !|   |  | |
          endif  !---------------------------------+   |  | |
          if (k .gt. 0)  Then !---------------------+  |  | ^
                                                   !|  |  | |
            Call rdReal (Frac(m,l),1,buff(i:j),    !|  |  | |
     *                                     merr)   !|  |  | |
            txt(19) = 'fraction "'//               !|  |  | |
     *                txt(18)(1:ltx-1)//'"'        !|  |  | |
            if (merr .ne. 0)              goto 2 !--+--+--+-+-ERROR+
            txt(19) = 'negative fraction "'//      !|  |  | |      v
     *                txt(18)(1:ltx-1)//'"'        !|  |  | |
            if (Frac(m,l) .lt. 0.)        goto 2 !--+--+--+-+-ERROR+
            txt(19) = 'fraction "'//               !|  |  | |      v
     *                txt(18)(1:ltx-1)//           !|  |  | |
     *                '" exceeds 1'                !|  |  | |
            if (Frac(m,l) .gt. 1.)        goto 2 !--+--+--+-+-ERROR+
            buff(k:j) = ' '                        !|  |  | |      v
                                                   !|  |  v ^
            Frac_Sum = Frac_Sum + Frac(m,l)        !|  |  | |
            txt(19) = 'sum of fractions exceeds'// !|  |  | |
     *                ' 1 or negative fraction xN' !|  |  | |
            if (Frac_Sum .gt. 1.001)      goto 2 !--+--+--+-+-ERROR+
                                                   !|  |  | |      v
          elseif (k .eq. -1)  Then  !---------------+  |  | |
                                                   !|  |  | |
            write (txt(19),119) txt(18)(1:ltx)     !|  |  | |
            goto 2 !----------------------------->--+--+--+-+-ERROR+
                                                   !|  |  | |      v
          else  !-----------------------------------+  |  v ^
                                                   !|  |  | |
c Fraction xN is not specified. This is allowed    !|  |  | |
c only once, otherwice we don't know the sum:      !|  |  | |
            izero_frac = izero_frac+1              !|  |  | |
            mzero_frac = m                         !|  |  | |
            txt(19) = 'more than one unspeci'//    !|  |  | |
     *                'fied fraction xN'           !|  |  | |
            if (izero_frac .gt. 1)        goto 2 !--+--+--+-+-ERROR+
                                                   !|  |  | |      v
          endif  !----------------------------------+  |  | |
        enddo  !=======================================+  | |
  71    continue   !<-------------------------------------+ |
        if (izero_frac .eq. 1) Then !------------------+    |
          Frac(mzero_frac,l) = 1. - Frac_Sum          !|    |
          txt(19) = 'fraction x? is zero or negative' !|    |
          write (txt(19)(11:11),'(i1)') mzero_frac    !|    |
          if (Frac(mzero_frac,l) .le. 0.)   goto 2 !---+----+-ERROR+
        else  !----------------------------------------+    |      v
          txt(19) = 'sum of fractions is not 1'       !|   !|
          if (Abs(Frac_Sum-1.) .gt. 0.001)  goto 2 !---+----+-ERROR+
        endif  !---------------------------------------+   !|      v
                                                           !|
                                                           !|
        Call ArgSearch (buff(1:ll),'rho=',' ',',',k,i,j)   !|
        if (k .gt. 0)  Then  !----------------------------+ |
          Call rdReal (rho(l),1,buff(i:j),m)             !| |
          txt(19) = 'material density rho (g/cm^3)'      !| |
          if (m .ne. 0)                   goto 2 !--------+-+-ERROR+
          txt(19) = 'negative material density (rho)'    !| |      v
          if (rho(l) .lt. 0.)             goto 2 !--------+-+-ERROR+
          txt(19) = 'rho is incompatible with '//        !| |      v
     *              'multiple fractions'                 !| |
          if (rho(l) .gt. 0. .And.                       !| |
     *     N_Frac(l) .gt. 1) goto 2 !---------------------+-+-ERROR+
          buff(k:j) = ' '                                !| |      v
        elseif (k .eq. -1)  Then  !-----------------------+ |
          txt(19) = 'duplicate keyword "rho="'           !| |
          goto 2 !----------------->----------------------+-+-ERROR+
        endif  !------------------------------------------+ |      v
                                                           !|
                                                           !|
        Call ArgSearch (buff(1:ll),'mdensity=',' ',',',    !|
     *                                               k,i,j)!|
        if (k .gt. 0) Then !------------------------------+ |
          Call rdReal (Mdensity(l),1,buff(i:j),m)        !| |
          txt(19) = 'magnetic atoms density (1/cm^3)'    !| |
          if (m .ne. 0)                   goto 2 !--------+-+-ERROR+
          txt(19) = 'negative magnetic atoms density'    !| |      v
          if (Mdensity(l) .lt. 0.)        goto 2 !--------+-+-ERROR+
          buff(k:j) = ' '                                !| |      v
          Mdensity_Specified=.True.                      !| |
        elseif (k .eq. -1)  Then  !-----------------------+ |
          txt(19) = 'duplicate keyword "mdensity="'      !| |
          goto 2 !----------------->----------------------+-+-ERROR+
        else   !------------------------------------------+ |      v
          Mdensity_Specified=.False.                     !| |
        endif  !------------------------------------------+ |
                                                           !|
                                                           !|
        Call ArgSearch (buff(1:ll),'mshare=',' ',',',k,i,j)!|
        if (k .gt. 0)  Then  !----------------------------+ |
          txt(19) = 'magnetic density and share'//       !| |
     *              ' specified simultaneously'          !| |
          if (Mdensity_Specified)        goto 2 !---------+-+-ERROR+
          Call rdReal (Mshare(l),1,buff(i:j),m)          !| |      v
          txt(19) = 'share of magnetic atoms in'//       !| |
     *              ' the material'                      !| |
          if (m .ne. 0)                    goto 2 !-------+-+-ERROR+
          txt(19) = 'share of magnetic atoms not'//      !| |      v
     *              ' in [0.-1.] range'                  !| |
          if (Mshare(l) .lt. 0.)           goto 2 !-------+-+-ERROR+
          if (Mshare(l) .gt. 1.)           goto 2 !-------+-+-ERROR+
          txt(19) = 'share of magnetic atoms'//          !| |      v
     *              ' requires material code'            !| |
          if (Len_Trim(Code(1,l)) .eq. 0 .AND.           !| |
     *                  Mshare(l) .gt. 0.) goto 2 !-------+-+-ERROR+
          buff(k:j) = ' '                                !| |      v
        elseif (k .eq. -1)  Then  !-----------------------+ |
          txt(19) = 'duplicate keyword "mshare="'        !| |
          goto 2 !----------------->----------------------+-+-ERROR+
        endif  !------------------------------------------+ |      v
                                                           !|
                                                           !^
        Call VectorSearch (buff(1:ll), 'mvector=',         !|
     *                             Mvect(1,l), txt(19), *2)!|
                                                           !|
c Below are:                                               !|
c  (3*lambda/8*pi)*F_10/e_radius                           !|
c  (3*lambda/8*pi)*F_11/e_radius                           !|
c  (3*lambda/8*pi)*F_1T/e_radius                           !|
                                                           !|
        Fij = 'F10'                                        !|
        Call ComplexSearch (buff(1:ll), 'f10=', f10(l),    !|
     *                                         txt(19), *2)!|
        if (Aimag(f10(l)).lt.0.)     goto 31               !|
                                                           !|
                                                           !^
        Fij = 'F11'                                        !|
        Call ComplexSearch (buff(1:ll), 'f11=', f11(l),    !|
     *                                         txt(19), *2)!|
        if (Aimag(f11(l)).lt.0.)     goto 31               !|
                                                           !|
                                                           !^
        Fij = 'F1T'                                        !|
        Call ComplexSearch (buff(1:ll), 'f1t=', f1t(l),    !|
     *                                         txt(19), *2)!|
        if (Aimag(f1t(l)).lt.0.)     goto 31               !|
                                                           !|
        if (Abs(f10(l)).gt.1.E-32 .and.                    !|
     *      Abs(f11(l)).lt.1.E-32 .and.                    !|
     *      Abs(f1t(l)).lt.1.E-32)   goto 34               !|
                                                           !|
        if (Abs(f11(l)).gt.1.E-32 .and.                    !|
     *      Abs(f1t(l)).gt.1.E-32) Then  !---------------+  |
          if (Abs((f11(l)-f1t(l))/f11(l)) .lt. 1.E-3)   !|  |
     *                                            goto 9!|  |
        endif !------------------------------------------+  |
                                                           !|
        if (Abs(f10(l)).gt.1.E-32 .or.                     !|
     *      Abs(f11(l)).gt.1.E-32 .or.                     !|
     *      Abs(f1t(l)).gt.1.E-32) Then  !---------------+  |
          if (Abs(Mvect(1,l)) .lt. 1.E-32 .and.         !|  |
     *        Abs(Mvect(2,l)) .lt. 1.E-32 .and.         !|  |
     *        Abs(Mvect(3,l)) .lt. 1.E-32)       goto 10!|  |
        endif !------------------------------------------+  |  |
                                                           !|  vE
c Test for unknown keywords (all the found keywords         |  |N
c and respecive values are erased:                          |  |D
        Call txtShift (buff(1:ll),j,i)      !i=Len_Trim(vm) |  |
        if (i .gt. 0)  goto 11 !--(unknown keyword)---------+--+-ERROR+
        if (t_found .eq. 0)   Then  !------------------+    |  |      v
          txt(19) = 'unspecified layer thickness "t"' !|    |  |
          goto 2 !------------------------------->-----+----+--+-ERROR+
        endif !----------------------------------------+    |  |      v
c Read next profile line:                                   |  |      v
        goto 5   !->----------------------------------------+  |
  4     continue    !<------------------<----------------------+
        if (ncycle .ne. 0)  goto 8 !-(loops not closed)------ERROR+
        Close   (unit=lun)                                       !v
        goto 30 !-------------------------------------------------------+
                                                                       !|
c============================================================           |
c                         E R R O R S                                   |
c------------------------------------------------------------           |
  1     continue                                                       !|
        write   (txt,21) progname(1:lpn), file(1:lfile)                !|
  21    format  (a,' cannot open top layer profile:'//a///)            !|
        goto 29  !----------------------------------------------------+ |
                                                                     !| |
  2     continue                                                     !| |
        i = Len_Trim(txt(19))                                        !| |
        if (i .lt. 1)  i=1                                           !| |
        if (i .gt. 67) i=67           !80-2-Len_Trim('Problem: []')  !| |
        buff(1:i) = txt(19)(1:i)                                     !| |
        j = Len_Trim(buff(2001:4000))                                !| |
        if (j .lt. 1)  j=1                                           !| |
        if (j .gt. 70) j=70           !80-2-Len_Trim('Line: []')     !| |
        write   (txt,22) progname(1:lpn), file(1:lfile),             !| |
     *                   buff(1:i),line, buff(2000+1:2000+j)         !| |
  22    format  (a,' has problem with top layer profile:'/a/
     *           'Problem: [',a,']'/
     *           'Line Nr: [',i4,']'/
     *           'Line: [',a,']')                                    !| |
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  11    continue                                                     !| |
        if (i .lt. 1)  i=1                                           !| |
        j = Index(buff(1:i),' ')   !find 1st word in unrecognized str | |
        if (j .gt. 1)  i=j                                           !| |
        if (i .gt. 62) i=62         !80-2-Len_Trim('Cannot parse: []')| |
        write   (txt,49) progname(1:lpn),file(1:lfile),              !| |
     *                   buff(1:i), line                             !| |
  49    format  (a,' has problem with top layer profile:'/a/
     *           'Unknown word: [',a,']'/
     *           'Line Nr: [',i4,']')                                !| |
        nlines = 4                                                   !| |
        goto 69 !------------------------->---------------------------+ |
                                                                     !| |
  19    continue                                                     !| |
        if (i .lt. 1)  i=1                                           !| |
        j = Index(buff(1:i),' ')   !find 1st word in unrecognized str | |
        if (j .gt. 1)  i=j                                           !| |
        if (i .gt. 64) i=64         !80-2-Len_Trim('Extra chars []')  | |
        write   (txt,59) progname(1:lpn),file(1:lfile),              !| |
     *                   buff(1:i), line                             !| |
  59    format  (a,' has problem with top layer profile:'/a/
     *           'Extra chars [',a,']'/
     *           'found while reading "period=" or "end period" line.'/
     *           'Line Nr: [',i4,']'//
     *           'Period specifications cannot be on the',1x,
     *                         'same line with layers data')         !| |
        nlines = 8                                                   !| |
        goto 69  !--------------------------------------------->------++|
                                                                     !|||
  3     continue                                                     !|||
        write   (txt,23) progname(1:lpn), file(1:lfile), N_Top_Max   !|v|
  23    format  (a,' has problem with top layer profile:'//a//
     *           'Buffer size (',i4,' layers) exceeded')             !| |
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  6     continue                                                     !| |
        write   (txt,26) progname(1:lpn), file(1:lfile), line        !| |
  26    format  (a,' has problem with top layer profile:'//a//
     *           'Attempt to start level-3 ',
     *           'enclosed loop at line=',i4)                        !| |
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  7     continue                                                     !| |
        write   (txt,27) progname(1:lpn), file(1:lfile), line        !| |
  27    format  (a,' has problem with top layer profile:'//a//
     *          'Attempt to close non-existing ',
     *          'loop at line=',i4)                                  !| |
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  8     continue                                                     !| |
        write   (txt,28) progname(1:lpn), file(1:lfile)              !| |
  28    format  (a,' has problem with top layer profile:'//a//
     *          'Loop(s) not closed at EOF')                         !| |
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  9     continue                                                     !| |
c The F11=F1T input is causing singular scattering matrix and as a   !| |
c result the IFAIL=1 error from F04ADF:                              !| |
        write   (txt,24) progname(1:lpn), file(1:lfile), line        !| |
  24    format  (a,' has problem with top layer profile:'//a//
     *          'Non-physical input F11=F1T at line=',i4)            !| |
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  31     continue                                                    !| |
c The Aimag(Fij)<0 input is causing diverging exponents              !| |
c and as a result reflectivity exceeding 1:                          !| |
        write   (txt,41) progname(1:lpn), file(1:lfile), Fij, line   !v v
  41    format  (a,' found non-physical input in top layer profile:'/
     *           a,'.'/
     *          'Negative data for Im(',a,') at line=',i4,'.'/
     *          'This may result in diverging exponents. The sign'/
     *          'of Im(Fij) must be the same as for x0, Im(x0)>0.')  !v v
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  34     continue                                                    !| |
c The F10#0 and F11=F1T=0 input is causing singular scattering matrix!| |
c (because the layer is detected as magnetic!) and as a result       !| |
c the IFAIL=1 error from F04ADF:                                     !| |
        write   (txt,44) progname(1:lpn), file(1:lfile), line        !| |
  44    format  (a,' has problem with top layer profile:'//a//
     *        'Non-physical input F10#0 while F11=F1T=0 at line=',i4)!| |
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  10    continue                                                     !| |
c When one specifies non-zero F10,F11, or F1T and does not specify   !| |
c the magnetizatrion direction, such input is resulted in a singular !| |
c scattering matrix and the IFAIL=1 error from F04ADF:               !| |
        write   (txt,25) progname(1:lpn), file(1:lfile), line        !v v
  25    format  (a,' has problem with top layer profile:'/
     *           a/
     *          'Non-zero F10,F11,F1T are specified, but'/
     *          'no magnetization direction is given'/
     *          'at  line=',i4)                                      !v v
        goto 29 !-----------------------------------------------------+ |
                                                                     !| |
  61    continue                                                     !| |
        j = Len_Trim(buff(2001:4000))                                !| |
        if (j .lt. 1)  j=1                                           !| |
        if (j .gt. 70) j=70           !80-2-Len_Trim('Line: []')     !| |
        write   (txt,66) progname(1:lpn), file(1:lfile),             !| |
     *                   line, buff(2000+1:2000+j)                   !v v
  66    format  (a,' has problem with top layer profile:'/a/
     *  'Problem: [both transition layer and roughness specified]'/
     *  'Line Nr: [',i4,']'/
     *  'Line: [',a,']')                                             !v v
        goto 29  !----------------------------------------------------+ |
                                                                     !| |
  29    continue !<---------------------------------------------------+ |
        nlines = 5                                                     !|
  69    continue                                                       !|
        Call    Message (txt,nlines,2)                                 !|
        ifail = 1                                                      !|
                                                                       !|
  30    continue  !<----------------------------------------------------+
        Thickness_Top = Sngl(TTT)
        N_Top         = l

  99    continue
        if (iirezv .eq. 1) Then  !---------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
c
c =======================================================
c
        Subroutine Read_Sub99_M0 (Code, rho, C00, W0,
     *                            x0, F10, F11, F1T,
     *                            Mdensity, Mshare, Mvect,
     *                            Sigma, Transit, YesTr,
     *                            wrk, line, *, *, *)

        Complex*8       x0,
     *                  F10, F11, F1T

        Real*4          Mdensity, Mshare, Mvect(3),
     *                  Sigma, Transit, rpar(2),
     *                  W0, rho

        Integer         line, lpn, i, j, l, n, iexit

        Integer         N_Column
        External        N_Column

        Logical         YesTr
        Character       Code*(*),
     *                  wrk*(*),
     *                  C00*(*)

        Character       txt(20)*80
        Common  /msg/   txt


c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer             iHenkeCowan
        Common /HenkeCowan/ iHenkeCowan
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv .lt. 10) Then  !----+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Read_Sub99_M0'      !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim (progname)
c-----------------------
        txt(7) = 'substrate code'
        Call    LineLoop (2,line,wrk,*100)
        Call    txtShift (wrk,j,i)
        Code = wrk

c-----------------------
        txt(7) = 'External X0h database flag [-1:4]'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (iHenkeCowan,1,wrk,i)
        if (i .ne. 0)               goto 100
        if (iHenkeCowan .lt. -1)    goto 101
        if (iHenkeCowan .gt. 8)     goto 101

c-----------------------
        txt(7) = 'substrate density rho (g/cm^3)'
        rho = 0.
        Call    LineLoop (2,line,wrk,*100)
        if (Len_Trim(wrk) .gt. 0)  Then  !-----+
          if (n_Column(wrk) .gt. 0) Then !---+ |
             Call  rdReal  (rpar,1,wrk,i)   !| |
             if (i .ne. 0)        goto 100  !| |
             rho = rpar(1)                  !| |
             if (rho .lt. 0.)     goto 101  !| |
          endif  !---------------------------+ |
        endif  !-------------------------------+

c-----------------------
        txt(7) = 'substrate x0 (range: |x0| < 0.1)'
        Call    LineLoop (2,line,wrk,*100)
        n = Len_Trim (wrk)
        do i=1,n  !===============================+
          if (wrk(i:i) .eq. '('  .or.            !|
     *        wrk(i:i) .eq. ')') wrk(i:i) = ' '  !|
        enddo  !==================================+
c Process the inputs like: (a+i*b) or (a-i*b):
c       i = Index(wrk(1:n),'+i*')
c       if (i .gt. 0) wrk(i:i+2) = ','
c       i = Index(wrk(1:n),'-i*')
c       if (i .gt. 0) wrk(i:i+2) = ','
        n = N_Column (wrk)
        if (n .gt. 0)  Then  !-----------------------+
          Call  RdReal (rpar,2,wrk,i)               !|
          if (i .ne. 0)           goto 100          !|
c What to do at absorption edges                     |
c where xr0 may change the sign?                     |
          x0  = Cmplx (-abs(rpar(1)), abs(rpar(2))) !|
          C00 = 'x'                                 !|
          Code= 'Unknown'                           !|
          if (abs(x0).gt.0.1)     goto 101          !|
        else   !-------------------------------------+
          C00 = ' '                                 !|
        endif  !-------------------------------------+
        if (Code .eq. ' '  .AND.
     *       C00 .eq. ' ') Then  !-----------------+
          write   (txt,55) progname(1:lpn)        !|
  55      format  (a,' error:'/
     *           'neither code nor x0'/
     *           'specified for the substrate')   !^
          Call  Message (txt,3,2)                 !|
          goto 28                                 !|
        endif  !-----------------------------------+

c-----------------------
        txt(7) = 'substrate Debye-Waller factor "w0"'//
     *                                      '(range: [0.-99.], w0 > 0)'
        W0 = 1.
        Call    LineLoop (2,line,wrk,*100)
        if (Len_Trim(wrk) .gt. 0)  Then  !-----+
          if (n_Column(wrk) .gt. 0) Then !---+ |
             Call  rdReal  (rpar,1,wrk,i)   !| |
             if (i .ne. 0)        goto 100  !| |
             W0 = rpar(1)                   !| |
             if (W0 .le. 0.)      goto 101  !| |
             if (W0 .gt. 99.)     goto 101  !| |
          endif  !---------------------------+ |
        endif  !-------------------------------+
c-----------------------
        txt(7) = 'substrate rms roughness "sigma" (sigma >= 0)'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (rpar,1,wrk,i)
        if (i .ne. 0)           goto 100
        Sigma = rpar(1)
        if (Sigma .lt. 0.)      goto 101
c-----------------------
        if (YesTr)  Then  !----------------------------------------+
          txt(7) = 'substrate transition layer "tr" (tr >= 0)'    !|
          Call  LineLoop (2,line,wrk,*100)                        !|
          Call  rdReal  (rpar,1,wrk,i)                            !|
          if (i .ne. 0)           goto 100                        !|
          Transit = rpar(1)                                       !|
          if (Transit .lt. 0.)    goto 101                        !|
                                                                  !|
          if (Transit .gt. 0.  .AND.  Sigma .gt. 0.) Then !-----+  |
            write   (txt,66) progname(1:lpn)                   !|  |
  66        format  (a,' error:'/
     *              'both transition layer and roughness'/
     *              'specified for the substrate')             !^  ^
            Call  Message (txt,3,2)                            !|  |
            goto 28                                            !|  |
          endif  !----------------------------------------------+  |
        endif  !---------------------------------------------------+
c-----------------------
        txt(7) = 'substrate magnetization direction (Mk,Ms,Mp)'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (Mvect,3,wrk,i)
        if (i .ne. 0)           goto 100
c-----------------------
        txt(7) = 'normalized substrate magnetic amplitude '//
     *           '(3*l/8*pi)*F10/r0'
        Call    LineLoop (2,line,wrk,*100)
        F10 = (0., 0.)
        l = Len_Trim(wrk)
        do i=1,l  !===============================+
          if (wrk(i:i) .eq. '('  .or.            !|
     *        wrk(i:i) .eq. ')') wrk(i:i) = ' '  !|
        enddo  !==================================+
        Call    rdReal  (rpar,2,wrk,i)
        if (i .ne. 0)           goto 100
        F10 = Cmplx(rpar(1),rpar(2))
        if (Aimag(F10).lt.0.) Then !----------------------------+
c The Aimag(Fij)<0 input is causing diverging exponents         |
c and as a result reflectivity exceeding 1:                     |
          write (txt,41) progname(1:lpn), 'F10'                !|
  41      format  (a,' error:'/
     *    'Non-physical input for the substrate:'/
     *    'Negative data for Im(',a,') may result in diverging'/
     *    'exponents. The sign of Im(Fij) must be the same as'/
     *    'for x0, Im(x0)>0.')                                 !|
          Call  Message (txt,5,2)                              !|
          goto 28                                              !|
        endif !-------------------------------------------------+
c-----------------------
        txt(7) = 'normalized substrate magnetic amplitude '//
     *           '(3*l/8*pi)*F11/r0'
        Call    LineLoop (2,line,wrk,*100)
        F11 = (0., 0.)
        l = Len_Trim(wrk)
        do i=1,l  !===============================+
          if (wrk(i:i) .eq. '('  .or.            !|
     *        wrk(i:i) .eq. ')') wrk(i:i) = ' '  !|
        enddo  !==================================+
        Call    rdReal  (rpar,2,wrk,i)
        if (i .ne. 0)           goto 100
        F11 = Cmplx(rpar(1),rpar(2))
        if (Aimag(F11).lt.0.) Then !----------------------------+
c The Aimag(Fij)<0 input is causing diverging exponents         |
c and as a result reflectivity exceeding 1:                     |
          write (txt,41) progname(1:lpn), 'F11'                !|
          Call  Message (txt,5,2)                              !|
          goto 28                                              !|
        endif !-------------------------------------------------+
c-----------------------
        txt(7) = 'normalized substrate magnetic amplitude '//
     *           '(3*l/8*pi)*F1T/r0'
        Call    LineLoop (2,line,wrk,*100)
        F1T = (0., 0.)
        l = Len_Trim(wrk)
        do i=1,l  !===============================+
          if (wrk(i:i) .eq. '('  .or.            !|
     *        wrk(i:i) .eq. ')') wrk(i:i) = ' '  !|
        enddo  !==================================+
        Call    rdReal  (rpar,2,wrk,i)
        if (i .ne. 0)           goto 100
        F1T = Cmplx(rpar(1),rpar(2))
        if (Aimag(F1T).lt.0.) Then !----------------------------+
c The Aimag(Fij)<0 input is causing diverging exponents         |
c and as a result reflectivity exceeding 1:                     |
          write (txt,41) progname(1:lpn), 'F1T'                !|
          Call  Message (txt,5,2)                              !|
          goto 28                                              !|
        endif !-------------------------------------------------+
c-----------------------
        if (Abs(F10).gt.1.E-32 .and.
     *      Abs(F11).lt.1.E-32 .and.
     *      Abs(F1t).lt.1.E-32) Then  !-------------------------+
c When one specifies non-zero F10,F11, or F1T and does          |
c not specify the magnetizatrion direction, such input          |
c is resulted in a singular scattering matrix and the           |
c IFAIL=1 error from F04ADF:                                    |
          write   (txt,26) progname(1:lpn)                     !|
  26      format  (a,' error:'/
     *    'Non-physical input for the substrate:'/
     *    'F10#0 while F11=F1T=0')                             !|
          Call  Message (txt,3,2)                              !|
          goto 28                                              !|
        endif !-------------------------------------------------+
        if (Abs(F11).gt.1.E-32 .and.
     *      Abs(F1T).gt.1.E-32) Then  !-------------------------+
          if (Abs((F11-F1T)/F11) .lt. 1.E-3) Then !-----------+ |
c The F11=F1T input is causing singular scattering            | |
c matrix and as a result the IFAIL=1 error from F04ADF:       | |
            write   (txt,24) progname(1:lpn)                 !| |
  24        format  (a,' error:'/
     *      'Non-physical input F11=F1T for the substrate')  !| |
            Call  Message (txt,2,2)                          !| |
            goto 28                                          !| |
          endif !---------------------------------------------+ |
        endif !-------------------------------------------------+

        if (Abs(F10).gt.1.E-32 .or.
     *      Abs(F11).gt.1.E-32 .or.
     *      Abs(F1t).gt.1.E-32) Then  !-------------------------+
          if (Abs(Mvect(1)) .lt. 1.E-32 .and.                  !|
     *        Abs(Mvect(2)) .lt. 1.E-32 .and.                  !|
     *        Abs(Mvect(3)) .lt. 1.E-32) Then !---------------+ |
c When one specifies non-zero F10,F11, or F1T and does        | |
c not specify the magnetizatrion direction, such input        | |
c is resulted in a singular scattering matrix and the         | |
c IFAIL=1 error from F04ADF:                                  | |
            write   (txt,25) progname(1:lpn)                 !| |
  25        format  (a,' error:'/
     *      'Non-zero F10,F11,F1T are specified for substrate'/
     *      'but no magnetization direction is given')       !| |
            Call  Message (txt,3,2)                          !| |
            goto 28                                          !| |
          endif !---------------------------------------------+ |
        endif !-------------------------------------------------+
c-----------------------
        txt(7) = 'substrate magnetic atoms density (1/cm^3)'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (rpar,1,wrk,i)
        if (i .ne. 0)           goto 100
        Mdensity = rpar(1)
        if (Mdensity .lt. 0.)   goto 101
c-----------------------
        txt(7) = 'share of magnetic atoms in substrate [0.-1.]'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (rpar,1,wrk,i)
        if (i .ne. 0)           goto 100
        Mshare = rpar(1)
        if (Mshare .lt. 0.)     goto 101
        if (Mshare .gt. 1.)     goto 101
        if (Mdensity .gt. 0.  .AND.  Mshare .gt. 0.) Then !-----+
          write   (txt,77) progname(1:lpn)                     !|
  77      format  (a,' error:'/
     *            'both magnetic density and share simultaneously'/
     *            'specified for the substrate')               !^
          Call  Message (txt,3,2)                              !|
          goto 28                                              !|
        endif  !------------------------------------------------+
        if (Len_Trim(Code) .eq. 0 .AND.  Mshare .gt. 0.) Then !-+
          write   (txt,88) progname(1:lpn)                     !|
  88      format  (a,' error:'/
     *            'specification of magneric share for'/
     *            'the substrate requires material code')      !^
          Call  Message (txt,3,2)                              !|
          goto 28                                              !|
        endif  !------------------------------------------------+
c-----------------------
        iexit = 0
  999   continue
        if (iirezv .eq. 1) Then  !---------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return iexit
  100   iexit = 1
        goto 999
  101   iexit = 2
        goto 999
  28    iexit = 3
        goto 999
        end

c =======================================================

        Subroutine VectorSearch (vm, tx, param, txt, *)
        Real            param(3)
        Integer         i, j, k, l, m, n
        Character       vm*(*), tx*(*), txt*(*)

        Call ArgSearch (vm,tx,')',')',k,i,j)

        if (k .gt. 0) Then  !------------------------------+
                                                          !|
          txt = tx                                        !|
          m   = Len_Trim(txt)                             !|
          l   = Len(txt)                                  !|
          txt(m:m) = ':'                                  !|
          m = m+2                                         !|
c Check for () and remove () if needed:                   !|
          txt(m:l) = 'missing (). Use:'                   !|
          n = Len_Trim(txt)+2                             !|
          txt(n:l) = tx                                   !|
          n = n+(m-2)                                     !|
          txt(n:l) = '(component1,component2,component3)' !|
          if (vm(i:i) .ne. '(')    goto 2 !----------------+---ERROR+
          if (vm(j:j) .ne. ')')    goto 2 !----------------+---ERROR+
          i = i + 1                                       !|        v
          j = j - 1                                       !|
          txt(m:l) = 'parse error'                        !|
          if (j .lt. i)  goto 2 !--------------------------+---ERROR+
                                                          !|        v
c Read param:                                             !|
          txt(m:l) = 'syntax error. Use:'                 !|
          n = Len_Trim(txt)+2                             !|
          txt(n:l) = tx                                   !|
          n = n+(m-2)                                     !|
          txt(n:l) = '(component1,component2,component3)' !|
          Call rdReal (param,3,vm(i:j),n)                 !|
          if (n .ne. 0)  goto 2 !--------------------------+---ERROR+
                                                          !|        v
          j   = j + 1   !restore j                        !|
          vm(k:j) = ' ' !erase the key!                   !|
                                                          !|
        elseif (k .eq. -1) Then  !-------------------------+
                                                          !|
          txt = tx                                        !|
          m   = Len_Trim(txt)+1                           !|
          l   = Len(txt)                                  !|
          txt(m:l) = ': duplicate keyword'                !|
          goto 2 !-----------------------------------------+---ERROR+
                                                          !|        v
        endif  !-------------------------------------------+
        return
c Return on error:
  2     continue
        return 1
        end

c =======================================================

        Subroutine ComplexSearch (vm, tx, param, txt, *)
        Complex         param
        Real            reim(2)
        Integer         i, j, k, l, m, n
        Character       vm*(*), tx*(*), txt*(*)

c First we search for a "normal" complex number (a,b):
        Call ArgSearch (vm,tx,')',')',k,i,j)

        if (k .gt. 0)  Then  !--------------------------+
                                                       !|
          txt = tx                                     !|
          m   = Len_Trim(txt)                          !|
          l   = Len(txt)                               !|
          txt(m:m) = ':'                               !|
          m = m+2                                      !|
c Check for () and remove () if needed:                !|
          txt(m:l) = 'missing (). Use:'                !|
          n = Len_Trim(txt)+2                          !|
          txt(n:l) = tx                                !|
          n = n+(m-2)                                  !|
          txt(n:l) = '(component1,component2)         '!|
          if (vm(i:i) .ne. '(')    goto 2 !-------------+---ERROR+
          if (vm(j:j) .ne. ')')    goto 2 !-------------+---ERROR+
          i = i + 1                                    !|        v
          j = j - 1                                    !|
          txt(m:l) = 'parse error'                     !|
          if (j .lt. i)  goto 2 !-----------------------+---ERROR+
                                                       !|        v
c Read param:                                          !|
          txt(m:l) = 'syntax error. Use:'              !|
          n = Len_Trim(txt)+2                          !|
          txt(n:l) = tx                                !|
          n = n+(m-2)                                  !|
          txt(n:l) = '(Re,Im)'                         !|
          Call rdReal (reim,2,vm(i:j),n)               !|
          if (n .ne. 0)  goto 2 !-----------------------+---ERROR+
          param = Cmplx(reim(1),reim(2))               !|        v
                                                       !|
          j   = j + 1   !restore j                     !|
          vm(k:j) = ' ' !erase the key!                !|
          return                                       !|
                                                       !|
        elseif (k .eq. -1) Then  !----------------------+
                                                       !|
          txt = tx                                     !|
          m   = Len_Trim(txt)+1                        !|
          l   = Len(txt)                               !|
          txt(m:l) = ': duplicate keyword'             !|
          goto 2 !--------------------------------------+---ERROR+
                                                       !|        v
        endif  !----------------------------------------+

c Alternatively search for a REAL number only:
        Call ArgSearch (vm,tx,' ',',',k,i,j)

        if (k .gt. 0)  Then  !--------------------------+
                                                       !|
          txt = tx                                     !|
          m   = Len_Trim(txt)                          !|
          l   = Len(txt)                               !|
          txt(m:m) = ':'                               !|
          m = m+2                                      !|
c Read param:                                          !|
          txt(m:l) = 'syntax error.'                   !|
          Call rdReal (reim,1,vm(i:j),n)               !|
          if (n .ne. 0)  goto 2 !-----------------------+---ERROR+
          param = Cmplx(reim(1),0.)                    !|        v
                                                       !|
          vm(k:j) = ' ' !erase the key!                !|
                                                       !|
        elseif (k .eq. -1) Then  !----------------------+
                                                       !|
          txt = tx                                     !|
          m   = Len_Trim(txt)+1                        !|
          l   = Len(txt)                               !|
          txt(m:l) = ': duplicate keyword'             !|
          goto 2 !--------------------------------------+---ERROR+
                                                       !|        v
        endif  !----------------------------------------+
        return
c Return on error:
  2     continue
        return 1
        end


