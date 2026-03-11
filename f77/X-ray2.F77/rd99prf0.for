c This module includes the following subroutines:
c   1. Read_Top99_0
c   2. Read_Sub99_0
c -- which read the surface layer profile and the
c substrate data, respectively, for ter_sl99 and
c trds_sl99 programs.
c
c=======================================================
c
        Subroutine Read_Top99_0 (File_Top,
     *                           N_Top, N_Top_Max,
     *                           Thickness, Thickness_Top,
     *                           Code, Frac,
     *                           N_Frac, N_Frac_Max,
     *                           rho, C00, x0, W0,
     *                           Identical,
     *                           Sigma,
     *                           Transit, YesTr,
     *                           ifail)
c--------------------------------------------------------
c       Reading of top layer profile -- smart version!
c
c                  Author: S.Stepanov
c -------------------------------------------------------
        Integer         N_Top,  N_Top_Max, N_Frac_Max
        Complex*8       x0(N_Top_Max)

        Real*8          TTT

        Real*4          Thickness(N_Top_Max),
     *                  Transit(N_Top_Max),     ! == Thickness_Tr
     *                  Sigma(N_Top_Max),
     *                  W0(N_Top_Max),
     *                  rho(N_Top_Max),
     *                  Frac(N_Frac_Max,N_Top_Max),
     *                  Thickness_Top, Frac_Sum,
     *                  Tmp(1)

        Logical         YesTr

        Integer         N_Frac(N_Top_Max),
     *                  Identical(N_Top_Max),
     *                  nbegin(2), nperiod(2), io_status,
     *                  lpn, iprf, lun, line, lena,
     *                  ncycle, ll, i, j, k, l, m,
     *                  l1, l2, ifail, iirezv,
     *                  lfile, izero_frac, mzero_frac,
     *                  ltx, merr, t_found, nlines

        Character       File_Top*(*), file*80,
     *                  Code(N_Frac_Max,N_Top_Max)*(*),
     *                  C00(N_Top_Max)*(*),
     *                  txt(20)*80

        Common  /msg/   txt
c -------------------------------------------------------
c Buffer space:
        Character       buff*4000
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv .lt. 10) Then  !----+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Read_Top99_0'       !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim (progname)
c -------------------------------------------------------
        l = 0                                   !make GNU compiler happy
        ifail = 0
        Thickness_Top = 0.
        TTT           = 0.
        iprf = Len_Trim(File_Top)
        if (iprf .eq. 0) goto 99  !-------return-------+
        file = File_Top                               !v
        file = file(1:iprf)//'.prf'
        lun = 4
        lfile = Len_Trim(file)

        Call OpenFile(file(1:lfile),lun,'read','old',io_status,*1)
        l       = 0
        line    = 0
        ncycle  = 0  !number of open loops
  5     continue    !<-----------------------<--------------+
        t_found = 1    !no problems with layer thickness "t"|
        txt(19) =  'file corrupt'                          !|
        buff  = ' '                                        !^
        line = line + 1                                    !|
        Read (lun,'(a)',err=2,end=4) buff !------>----------+--+
        ll = Len_Trim(buff)                                !|  |E
        if (ll .eq. 0) goto 5 !----->---------empty--->-----+  |N
                                                           !|  |D
c Remove hidden/unprintable characters (including <TAB>):   |  v
        do j=1,ll !================================+        |
          if (buff(j:j) .lt. ' ') buff(j:j) = ' ' !|        |
        enddo  !===================================+        |
        ll = Len_Trim(buff)                                !|
        if (ll .eq. 0) goto 5 !----->---------empty--->-----+
                                                           !|
c Remove comments at the end of the line:                   |
        j = Index(buff(1:ll),';')                          !|
        if (j .gt. 0) Then !-------+                        |
          buff(j:ll) = ' '        !|                        ^
          ll = Len_Trim(buff)     !|                       !|
          if (ll .eq. 0) goto 5 !--+----------empty--->-----+
        endif  !-------------------+                        |
        j = Index(buff(1:ll),'!')                          !|
        if (j .gt. 0) Then !-------+                        |
          buff(j:ll) = ' '        !|                        ^
          ll = Len_Trim(buff)     !|                       !|
          if (ll .eq. 0) goto 5 !--+----------empty--->-----+
        endif  !-------------------+                        |
c Compress the line:                                        |
c       Call Compress_Line (buff,ll,'  ',1)                !|
c See Serve.lib:                                            |
        Call ComprSp (buff,ll)                             !|
c See Xinput97.for:                                         |
        Call Compress_Line (buff,ll,' =',0)                !|
        Call Compress_Line (buff,ll,'= ',1)                !|
c Replace comma's by spaces:                                |2000/07/10
        do j=1,ll !================================+        |
          if (buff(j:j) .eq. ',') buff(j:j) = ' ' !|        |2003/04/30
        enddo  !===================================+        |(moved up)
        Call ComprSp (buff,ll)                             !|
                                                           !|
c I presume: ll<2000!                                       |
        if (ll .le. 2000) Then !---------------------------+|2003/04/30
          buff(2001:4000) = buff(1:2000)                  !||
        else !---------------------------------------------+|
          txt(19) = 'profile line exceeds 2000 characters'!||
          goto 2 !------------------------------->---------++-ERROR+
        endif !--------------------------------------------+|      v
c Switch to lower case (see Serve.lib):                     |
        Call Case (buff(1:ll),1)                           !|
c Search for the keywords:                                  |
c period= ...  end period                                   |
c t= code= x0= w0= rho= sigma= tr=  [x=                     ^
c                                    code2= x2=             |
c                                    code3= x3=             |
c                                    code4=  ]              |
                                                           !|
c See Xinput97.for:                                         |
        Call ArgSearch (buff(1:ll),'period=',' ',',',k,i,j)!|
        if (k .gt. 0)  Then  !-------------------------+    |
          ncycle = ncycle+1                           !|    |
          if (ncycle .gt. 2) goto 6  !-----------------+----+-ERROR+
          Call rdInt (nperiod(ncycle),1,buff(i:j),m)  !|    |      v
          txt(19) = 'period'                          !|    |
          if (m .ne. 0) goto 2  !----------------------+----+-ERROR+
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
        if (i .gt. 0  .or.  j .gt. 0)  Then  !---------+    |
          if (i .gt. 0) Then !---+                     |    |
            buff(i:i+9) = ' '   !|                     |    |
          else  !----------------+                     |    |
            buff(j:j+8) = ' '   !|                     |    |
          endif  !---------------+                     |    |
          ncycle = ncycle-1                           !|    ^
          if (ncycle .lt. 0)    goto 7  !--------------+----+-ERROR+
c The size of group to be multiplied:                  |    |      v
          lena = l-nbegin(ncycle+1)+1                 !|    |
c The total of layers after multiplication:            |    |
          l  = l + lena*(nperiod(ncycle+1)-1)         !|    |
          if (l .gt. N_Top_Max) goto 3   !-->----------+----+-ERROR+
          do    i=2,nperiod(ncycle+1)  !=========+     |    ^      v
          do    j=1,lena   !===================+ |     |    |
            l1 = nbegin(ncycle+1)+j-1         !| |     |    |
            l2 = l1 + lena*(i-1)              !| |     |    |
            Thickness(l2) = Thickness(l1)     !| |     |    ^
            TTT = TTT + Thickness(l1)         !| |     |    |
            N_Frac(l2)    = N_Frac(l1)        !| |     |    |
            do m=1,N_Frac_Max !=========+      | |     |    |
              Code(m,l2)  = Code(m,l1) !|      | |     |    |
              Frac(m,l2)  = Frac(m,l1) !|      | |     |    |
            enddo  !====================+      | |     |    |
            Identical(l2) = l1                !| |     |    |
            x0(l2)        = x0(l1)            !| |     |    |
            W0(l2)        = W0(l1)            !| |     |    |
            Sigma(l2)     = Sigma(l1)         !| |     |    |
            if (YesTr)                        !| |     |    |
     *      Transit(l2)   = Transit(l1)       !| |     |    |
            C00(l2)       = C00(l1)           !| |     |    |
            rho(l2)       = rho(l1)           !| |     |    |
          enddo  !=============================+ |     |    |
          enddo  !===============================+     |    |
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
        Thickness(l)   = 0.                                !|
        N_Frac(l)      = 1                                 !|
        do m=1,N_Frac_Max  !===+                            |
          Code(m,l)    = ' '  !|                            |
          Frac(m,l)    = 0.   !|                            |
        enddo  !===============+                            |
        Frac(1,l)      = 1.                                !|
        x0(l)          = (0.,0.)                           !^
        W0(l)          = 1.                                !|
        Sigma(l)       = 0.                                !|
        if (YesTr)                                         !|
     *      Transit(l) = 0.                                !|
        C00(l)         = ' '                               !|
        rho(l)         = 0.                                !|
        Identical(l)   = 0                                 !|
                                                           !|
        Call ArgSearch (buff(1:ll),'t=',' ',',',k,i,j)     !|
        if (k .gt. 0)  Then  !-------------------------+    |
          txt(19) = 'thickness'                       !|    |
          Call rdReal (Tmp,1,buff(i:j),m)             !|    |
          if (m .ne. 0)  goto 2 !----------------------+->--+-ERROR+
          Thickness(l) = Tmp(1)                       !|    |      v
          buff(k:j) = ' '                             !|    |
          txt(19) = 'zero layer thickness "t" - '//   !|    |
     *              'please comment out the line!'    !|    |
          if (Abs(Thickness(l)) .lt. 1.E-32) goto 2 !--+->--+-ERROR+
          txt(19) = 'negative layer thickness "t"'    !|    |      v
          if (Thickness(l) .lt. 0.) goto 2 !-----------+->--+-ERROR+
          TTT = TTT + Thickness(l)                    !|    |      v
          t_found = 1                                 !|    |
        elseif (k .eq. -1)  Then  !--------------------+    |
          txt(19) = 'duplicate keyword "t="'          !|    |
          goto 2 !-------------------------------------+->--+-ERROR+
        elseif (k .eq. 0)   Then  !--------------------+    |      v
          t_found = 0  !unspecified layer thickness "t"|    |
        endif  !---------------------------------------+    |
                                                           !|
                                                           !|
        Call ArgSearch (buff(1:ll),'w0=',' ',',',k,i,j)    !|
        if (k .gt. 0)  Then  !-------------------------+    |
          txt(19) = 'w0'                              !|    |
          Call rdReal (Tmp,1,buff(i:j),m)             !|    |
          if (m .ne. 0)       goto 2 !-----------------+----+-ERROR+
          W0(l)   = Tmp(1)                            !|    |      v
          txt(19) = 'negative "w0"'                   !|    |
          if (W0(l) .lt. 0.)  goto 2 !-----------------+----+-ERROR+
          txt(19) = 'too big "w0" (max=99)'           !|    |      v
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
          txt(19) = 'sigma'                           !|    |
          Call rdReal (Tmp,1,buff(i:j),m)             !|    |
          if (m .ne. 0)                 goto 2 !-------+----+-ERROR+
          Sigma(l) =  Tmp(1)                          !|    |      v
          txt(19) = 'negative rms roughness "sigma"'  !|    |
          if (Sigma(l) .lt. 0.)         goto 2 !-------+----+-ERROR+
          txt(19) = 'rms roughness "sigma" '//        !|    |      v
     *                 'exceeds layer thickness "t"'  !|    |
          if (Sigma(l) .gt. Thickness(l) .AND.        !|    |
     *        t_found .ne. 0)           goto 2 !-------+----+-ERROR+
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
            txt(19) = 'tr'                           !|  |  |
            Call rdReal (Tmp,1,buff(i:j),m)          !|  |  |
            if (m .ne. 0)                   goto 2 !--+--+--+-ERROR+
            Transit(l) = Tmp(1)                      !|  |  |      v
            txt(19) = 'negative transition '//       !|  |  |
     *                       'layer thickness "tr"'  !|  |  |
            if (Transit(l) .lt. 0.)         goto 2 !--+--+--+-ERROR+
            txt(19) = 'transition domain "tr" '//    !|  |  |      v
     *                'exceeds layer thickness "t"'  !|  |  |
            if (Transit(l) .gt. Thickness(l) .AND.   !|  |  |
     *          t_found .ne. 0)             goto 2 !--+--+--+-ERROR+
            buff(k:j) = ' '                          !|  |  |      v
            if ((Transit(l) .gt. 0.) .AND.           !|  |  |
     *          (Sigma(l).gt.0.))           goto 61 !-+--+--+-ERROR+
          elseif (k .eq. -1)  Then  !-----------------+  |  |
            txt(19) = 'duplicate keyword "tr="'      !|  |  |
            goto 2 !----------------------------->----+--+--+-ERROR+
          endif  !------------------------------------+  |  |      v
        endif  !-----------------------------------------+  |
                                                           !^
        Call ChiSearch (buff(1:ll), 'x0=', x0(l),          !|
     *                  C00(l)(1:1), txt(19), '0', *2)     !|
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
            Call ArgSearch (buff(1:ll),'code=',  !|    |    |
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
                                                   !|  |    |      v
          else  !-----------------------------------+  |    ^
                                                   !|  |    |
c Goto exit from the loop:                         !|  |    |
            if (m .gt. 1)  goto 71 !----------------+--+--+ |
c (if Code-1 is missing, the substrate code         |  |  | |
c  code will be used. Here we simply continue)      |  |  | |
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
            txt(19) = 'fraction "'//               !|  |  | |
     *                txt(18)(1:ltx-1)//'"'        !|  |  | |
            Call rdReal (Tmp,1,buff(i:j),merr)     !|  |  | |
            if (merr .ne. 0)            goto 2 !----+--+--+-+-ERROR+
            Frac(m,l) = Tmp(1)                     !|  |  | |      v
            txt(19) = 'negative fraction "'//      !|  |  | |
     *                txt(18)(1:ltx-1)//'"'        !|  |  | |
            if (Frac(m,l) .lt. 0.)      goto 2 !----+--+--+-+-ERROR+
            txt(19) = 'fraction "'//               !|  |  | |      v
     *                txt(18)(1:ltx-1)//           !|  |  | |
     *                '" exceeds 1'                !|  |  | |
            if (Frac(m,l) .gt. 1.)      goto 2 !----+--+--+-+-ERROR+
            buff(k:j) = ' '                        !|  |  | |      v
                                                   !|  |  v ^
            Frac_Sum = Frac_Sum + Frac(m,l)        !|  |  | |
            txt(19) = 'sum of fractions exceeds'// !|  |  | |
     *                ' 1 or negative fraction xN' !|  |  | |
            if (Frac_Sum .gt. 1.001)    goto 2 !----+--+--+-+-ERROR+
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
            if (izero_frac .gt. 1)      goto 2 !----+--+--+-+-ERROR+
                                                   !|  |  | |      v
          endif  !----------------------------------+  |  | |
        enddo  !=======================================+  | |
  71    continue   !<-------------------------------------+ |
        if (izero_frac .eq. 1) Then !------------------+    |
          Frac(mzero_frac,l) = 1. - Frac_Sum          !|    |
          txt(19) = 'fraction x? is zero or negative' !|    |
          write (txt(19)(11:11),'(i1)') mzero_frac    !|    |
          if (Frac(mzero_frac,l) .le. 0.) goto 2  !----+----+-ERROR+
        else  !----------------------------------------+    |      v
          txt(19) = 'sum of fractions is not 1'       !|   !|
          if (Abs(Frac_Sum-1.) .gt. 0.001) goto 2 !----+----+-ERROR+
        endif  !---------------------------------------+   !|      v
                                                           !|
                                                           !|
        Call ArgSearch (buff(1:ll),'rho=',' ',',',k,i,j)   !|
        if (k .gt. 0)  Then  !----------------------------+ |
          txt(19) = 'material density rho (g/cm^3)'      !| |
          Call rdReal (Tmp,1,buff(i:j),m)                !| |
          if (m .ne. 0)                 goto 2 !----------+-+-ERROR+
          rho(l) = Tmp(1)                                !| |      v
          txt(19) = 'negative material density "rho"'    !| |
          if (rho(l) .lt. 0.)           goto 2 !----------+-+-ERROR+
          txt(19) = 'using "rho" is incompatible '//     !| |      v
     *              'with multiple fractions'            !| |
          if (rho(l) .gt. 0. .and.                       !| |
     *     N_Frac(l) .gt. 1) goto 2 !---------------------+-+-ERROR+
          buff(k:j) = ' '                                !| |      v
        elseif (k .eq. -1)  Then  !-----------------------+ |
          txt(19) = 'duplicate keyword "rho="'           !| |
          goto 2 !----------------->----------------------+-+-ERROR+
        endif  !------------------------------------------+ |      v
                                                           !|  vE
c Test for unknown keywords (all the found keywords         |  |N
c and respecive values are erased:                          |  |D
        Call txtShift (buff(1:ll),j,i)      !i=Len_Trim(vm) |  |
        if (i .gt. 0)  goto 9  !--(unknown keyword)---------+--+-ERROR+
        if (t_found .eq. 0)   Then  !------------------+    |  |      v
          txt(19) = 'unspecified layer thickness "t"' !|    |  |
          goto 2 !------------------------------->-----+----+--+-ERROR+
        endif !----------------------------------------+    |  |      v
c Read next profile line:                                   |  |
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
        goto 29 !------------------------->---------------------------+ |
                                                                     !| |
  2     continue                                                     !| |
        i = Len_Trim(txt(19))                                        !| |
        if (i .lt. 1)  i=1                                           !| |
        if (i .gt. 67) i=67                                          !| |
        buff(1:i) = txt(19)(1:i)                                     !| |
        j = Len_Trim(buff(2001:4000))                                !| |
        if (j .lt. 1)  j=1                                           !| |
        if (j .gt. 70) j=70                                          !| |
        write   (txt,22) progname(1:lpn),file(1:lfile),              !| |
     *                   buff(1:i),line,buff(2000+1:2000+j)          !| |
  22    format  (a,' has problem with top layer profile:'/a/
     *           'Problem: [',a,']'/
     *           'Line Nr: [',i4,']'/
     *           'Line: [',a,']')                                    !| |
        goto 29 !------------------------->---------------------------+ |
                                                                     !| |
  9     continue                                                     !| |
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
        nlines = 7                                                   !| |
        goto 69  !--------------------------------------------->------++|
                                                                     !|||
  3     continue                                                     !|||
        write   (txt,23) progname(1:lpn), file(1:lfile),             !|v|
     *                   N_Top_Max                                   !| |
  23    format  (a,' has problem with top layer profile:'//a//
     *           'Buffer size (',i4,' layers) exceeded')             !| |
        goto 29                                                      !| |
                                                                     !| |
  6     continue                                                     !| |
        write   (txt,26) progname(1:lpn), file(1:lfile),             !| |
     *                   line                                        !| |
  26    format  (a,' has problem with top layer profile:'//a//
     *           'Attempt to start level-3 ',
     *           'enclosed loop at line=',i4)                        !| |
        goto 29 !------------------------->---------------------------+ |
                                                                     !| |
  7     continue                                                     !| |
        write   (txt,27) progname(1:lpn), file(1:lfile),             !| |
     *                   line                                        !| |
  27    format  (a,' has problem with top layer profile:'//a//
     *          'Attempt to close non-existing ',
     *          'loop at line=',i4)                                  !| |
        goto 29 !------------------------->---------------------------+ |
                                                                     !| |
  8     continue                                                     !| |
        write   (txt,28) progname(1:lpn),file(1:lfile)               !| |
  28    format  (a,' has problem with top layer profile:'//a//
     *          'Loop(s) not closed at EOF')                         !| |
        goto 29 !------------------------->---------------------------+ |
                                                                     !| |
  61    continue                                                     !| |
        j = Len_Trim(buff(2001:4000))                                !| |
        if (j .lt. 1)  j=1                                           !| |
        if (j .gt. 70) j=70                                          !| |
        write   (txt,66) progname(1:lpn), file(1:lfile),             !| |
     *                   line, buff(2000+1:2000+j)                   !| |
  66    format  (a,' has problem with top layer profile:'/a/
     *  'Problem: [both transition layer and roughness are',1x,
     *  'simultaneosly specified]'/'Line Nr: [',i4,']'/
     *  'Line: [',a,']')                                             !| |
        goto 29 !------------------------->---------------------------+ |
                                                                     !| |
  29    continue !<---------------------------------------------------+ |
        nlines = 5                                                     !|
  69    continue                                                       !|
        Call    Message (txt,nlines,2)                                 !|
        ifail = 1                                                      !|
                                                                       !|
  30    continue !<-----------------------------------------------------+
        Thickness_Top = Sngl(TTT)
        N_Top         = l

  99    continue
        if (iirezv .eq. 1) Then !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
c
c=======================================================
c
        Subroutine Read_Sub99_0 (Code, rho, C00,
     *                           x0, W0,
     *                           Sigma,
     *                           Transit, YesTr,        !Transit may be the same address as Sigma
     *                           wrk, line, *, *, *)

        Complex*8       x0
        Real*4          Sigma, Transit,
     *                  Tmp(2), rho, W0
        Integer         lpn, i, j, n,
     *                  iirezv, line,
     *                  iexit
        Integer         N_Column
        External        N_Column
        Logical         YesTr
        Character       Code*(*),
     *                  wrk*(*),
     *                  C00*(*),
     *                  txt(20)*80

        Common  /msg/   txt

c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
        Integer                 iHenkeCowan
        Common  /HenkeCowan/    iHenkeCowan
c -------------------------------------------------------
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
        if (istackrezv .lt. 10) Then  !----+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Read_Sub99_0'       !|
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
        iHenkeCowan = -1
        Call    LineLoop (2,line,wrk,*100)
        Call    rdInt   (iHenkeCowan,1,wrk,i)
        if (i .ne. 0)             goto 100
        if (iHenkeCowan .lt. -1)  goto 101
        if (iHenkeCowan .gt. 8)   goto 101

c-----------------------
        txt(7) = 'substrate density rho (g/cm^3)'
        rho = 0.
        Call    LineLoop (2,line,wrk,*100)
        if (Len_Trim(wrk) .gt. 0)  Then  !-----+
          if (n_Column(wrk) .gt. 0) Then !---+ |
             Call  rdReal  (Tmp,1,wrk,i)    !| |
             if (i .ne. 0)      goto 100    !| |
             rho = Tmp(1)                   !| |
             if (rho .lt. 0.)   goto 101    !| |
          endif  !---------------------------+ |
        endif  !-------------------------------+

c-----------------------
        txt(7) = 'substrate x0 (range: |x0| < 0.1)'
        x0 = (0.,0.)
        Call    LineLoop (2,line,wrk,*100)
        n = Len_Trim (wrk)
        do i=1,n  !==============================+
          if (wrk(i:i) .eq. '('  .or.           !|
     *        wrk(i:i) .eq. ')') wrk(i:i) = ' ' !|
        enddo  !=================================+
c Process the inputs like: (a+i*b) or (a-i*b):
c       i = Index(wrk(1:n),'+i*')
c       if (i .gt. 0) wrk(i:i+2) = ','
c       i = Index(wrk(1:n),'-i*')
c       if (i .gt. 0) wrk(i:i+2) = ','
        n = N_Column (wrk)
        if (n .gt. 0)  Then  !---------------------+
          Call  RdReal (Tmp,2,wrk,i)              !|
          if (i .ne. 0)         goto 100          !|
c What to do at absorption edges                   |
c where xr0 may change the sign?                   |
          x0  = Cmplx (-abs(Tmp(1)), abs(Tmp(2))) !|
          C00 = 'x'                               !|
          Code= 'Unknown'                         !|
          if (abs(x0).gt.0.1)   goto 101          !|
        else   !-----------------------------------+
          C00 = ' '                               !|
        endif  !-----------------------------------+
        if (Code.eq.' ' .AND. C00.eq.' ') Then  !--+
          write (txt,55)                          !v
  55      format(/'Error: neither code nor x0'//
     *           'specified for the substrate'/)  !^
          Call  Message (txt,5,2)                 !|
          goto 28                                 !|
        endif  !-----------------------------------+

c-----------------------
        txt(7) = 'substrate Debye-Waller factor "w0"'//
     *                                      '(range: [0.-99.], w0 > 0)'
        W0 = 1.
        Call    LineLoop (2,line,wrk,*100)
        if (Len_Trim(wrk) .gt. 0) Then  !------+
          if (n_Column(wrk) .gt. 0) Then !---+ |
             Call  rdReal  (Tmp,1,wrk,i)    !| |
             if (i .ne. 0)      goto 100    !| |
             W0 = Tmp(1)                    !| |
             if (W0 .le. 0.)    goto 101    !| |
             if (W0 .gt. 99.)   goto 101    !| |
          endif  !---------------------------+ |
        endif  !-------------------------------+
c-----------------------
        txt(7) = 'substrate rms roughness "sigma" (sigma >= 0)'
        Sigma   = 0.
        Transit = 0.                                    !this may be the same address as Sigma
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (Tmp,1,wrk,i)
        if (i .ne. 0)           goto 100
        Sigma = Tmp(1)
        if (Sigma .lt. 0.)      goto 101

c-----------------------
        if (YesTr)  Then  !----------------------------------------+
          txt(7) = 'substrate transition layer "tr" (tr >= 0)'    !|
          Call  LineLoop (2,line,wrk,*100)                        !|
          Call  rdReal  (Tmp,1,wrk,i)                             !|
          if (i .ne. 0)         goto 100                          !|
          Transit = Tmp(1)                                        !|
          if (Transit .lt. 0.)  goto 101                          !|
                                                                  !|
          if (Transit.gt.0. .AND. Sigma.gt.0.) Then  !----------+  |
            write (txt,66)                                     !v  v
  66        format(/'Error: both transition layer and roughness'//
     *             'specified for the substrate'/)             !^  ^
            Call  Message (txt,5,2)                            !|  |
            goto 28                                            !|  |
          endif  !----------------------------------------------+  |
        endif  !---------------------------------------------------+

c-----------------------
        iexit = 0
  999   continue
        if (iirezv .eq. 1) Then !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return iexit
  100   continue
        iexit = 1               !Read error
        goto 999
  101   continue
        iexit = 2               !Parameter not in range
        goto 999
  28    continue                !Parameters conflict
        iexit = 3
        goto 999
        end
