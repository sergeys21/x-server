c This file rdm7prfh.for is called from gid_slm7.
c As compared to the previous version rd99prfh.for,
c it adds option for the phase of structure amplitude.
c
c The following 2 subroutines are contained:
c     -- Read_TopM7_0h
c     -- Read_SubM7_0h
c
c=======================================================

        Subroutine Read_TopM7_0h (
     *       File_Top,
     *       N_Top, N_Top_Max,
     *       Thickness,                         ! passed from index=1 (index=1 is top layer); thicknesses of layers
     *       Thickness_Top,                     ! accumulated surface layer thinkness
     *       Code, Frac,                        ! passed from index=1 (index=1 is top layer)
     *       N_Frac,                            ! passed from index=1 (index=1 is top layer)
     *       N_Frac_Max,                        ! N_Frac_Max=4
     *       x0, xh, xf,                        ! passed from index=2 (index=1 is vacuum)    <==
     *       W0, Wh,                            ! passed from index=1 (index=1 is top layer)
     *       Daa,                               ! passed from index=2 (index=1 is vacuum)    <==
     *       Sigma,                             ! passed from index=1 (index=1 is top layer)
     *       xs,                                ! passed from index=2 (index=1 is vacuum)    <==; added in M7
     *       C0h,                               ! passed from index=1 (index=1 is top layer); flags if x0,xh are specified.
     *       ifail)
c--------------------------------------------------------
c       Reading of top layer profile -- smart version!
c
c                  Author: S.Stepanov
c -------------------------------------------------------
        Integer         N_Top,  N_Top_Max, N_Frac_Max
        Complex*8       x0(N_Top_Max),
     *                  xh(N_Top_Max)

        Real*8          TTT

        Real*4          xf(N_Top_Max),          ! phase shift between xrh & xih (non-cubic)
     *                  xs(N_Top_Max),          ! xh phase shift in layer (Kaganer, monolayers, M7)
     *                  Thickness(N_Top_Max),   ! thicknesses of layers
     *                  Sigma(N_Top_Max),       ! rms rougness of upper interface
     *                  Daa(N_Top_Max),         ! normal lattice strain
     *                  W0(N_Top_Max),
     *                  Wh(N_Top_Max),
     *                  Frac(N_Frac_Max,N_Top_Max),
     *                  Thickness_Top, Frac_Sum,
     *                  Tmp(1)

        Integer         N_Frac(N_Top_Max), io_status,
     *                  nbegin(2), nperiod(2),
     *                  lpn, iprf, lun, line, lena,
     *                  ncycle, ll, i, j, k, lyr, m,
     *                  l1, l2, ifail, iirezv,
     *                  lfile, izero_frac, mzero_frac,
     *                  ltx, merr, t_found, nlines

        Character       File_Top*(*),file*80,
     *                  Code(N_Frac_Max,N_Top_Max)*(*),
     *                  C0h(N_Top_Max)*(*),
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
          progname = 'Read_TopM7_0h'      !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+
        lpn = Len_Trim (progname)
c -------------------------------------------------------
        ifail         = 0
        lyr           = 0
        line          = 0
        ncycle        = 0                               !number of open loops
        Thickness_Top = 0.
        TTT           = 0.

        iprf          = Len_Trim(File_Top)
        if (iprf .eq. 0) goto 99                        !return
        file          = File_Top(1:iprf)//'.prf'
        lun           = 4
        lfile         = Len_Trim(file)

        Call OpenFile(file(1:lfile),lun,'read','old',io_status,*1)
  5     continue    !<-----------------------<----------------+
        t_found = 1      !no problems with layer thickness "t"|
        txt(19) =  'file corrupt'                            !|
        buff  = ' '                                          !^
        line = line + 1                                      !|
        Read (lun,'(a)',err=2,end=4) buff !------>------------+--+
        ll = Len_Trim(buff)                                  !|  |E
        if (ll .eq. 0) goto 5 !----->---------empty--->-------+  |N
                                                             !|  |D
c Remove hidden/unprintable characters (including <TAB>):     |  v
        do j=1,ll !================================+          |
          if (buff(j:j) .lt. ' ') buff(j:j) = ' ' !|          |
        enddo  !===================================+          |
        ll = Len_Trim(buff)                                  !|
        if (ll .eq. 0) goto 5 !----->---------empty--->-------+
                                                             !|
c Remove comments at the end of the line:                     |
        j = Index(buff(1:ll),';')                            !|
        if (j .gt. 0) Then !-------+                          |
          buff(j:ll) = ' '        !|                          ^
          ll = Len_Trim(buff)     !|                         !|
          if (ll .eq. 0) goto 5 !--+----------empty--->-------+
        endif  !-------------------+                          |
        j = Index(buff(1:ll),'!')                            !|
        if (j .gt. 0) Then !-------+                          |
          buff(j:ll) = ' '        !|                          ^
          ll = Len_Trim(buff)     !|                         !|
          if (ll .eq. 0) goto 5 !--+----------empty--->-------+
        endif  !-------------------+                          |
c Compress the line:                                          |
c       Call Compress_Line (buff,ll,'  ',1)                  !|
c See Serve.lib:                                              |
        Call ComprSp (buff,ll)                               !|
c See xinput97.for:                                           |
        Call Compress_Line (buff,ll,' =',0)                  !|
        Call Compress_Line (buff,ll,'= ',1)                  !|
c Replace comma's by spaces:                                  |2000/07/10
        do j=1,ll !================================+          |
          if (buff(j:j) .eq. ',') buff(j:j) = ' ' !|          |2003/04/30
        enddo  !===================================+          |(moved up)
        Call ComprSp (buff,ll)                               !|
                                                             !|
c I presume: ll<2000!                                         |
        if (ll .le. 2000) Then !---------------------------+  |2003/04/30
          buff(2001:4000) = buff(1:2000)                  !|  |
        else !---------------------------------------------+  |
          txt(19) = 'profile line exceeds 2000 characters'!|  |
          goto 2 !------------------------------->---------+--+-ERROR+
        endif !--------------------------------------------+  |      v
c Switch to lower case (see Serve.lib):                       |
        Call Case (buff(1:ll),1)                             !|
c Search for the keywords:                                    |
c period= ...  end period                                     |
c t= code= [x= code2= x2= code3= x3= code4= x4=]              |
c x0= xh= xhdf= w0= wh= da/a= sigma= xhphase=                 ^
                                                             !|
c See xinput97.for:                                           |
        Call ArgSearch (buff(1:ll),'period=',' ',',',k,i,j)  !|
        if (k .gt. 0)  Then  !-------------------------+      |
          ncycle = ncycle+1                           !|      |
          if (ncycle .gt. 2)    goto 6  !--------------+------+-ERROR+
          Call rdInt (nperiod(ncycle),1,buff(i:j),m)  !|      |      v
          txt(19) = 'period'                          !|      |
          if (m .ne. 0) goto 2  !----------------------+------+-ERROR+
          txt(19) = 'zero/negative period'            !|      |      v
          if (nperiod(ncycle) .lt. 1) goto 2  !--------+------+-ERROR+
          nbegin(ncycle) = lyr+1                      !|      |      v
          buff(k:j) = ' '                             !|      |
          Call txtShift (buff(1:ll),j,i)              !|      ^
c -- the call returns: i=Len_Trim(vm)                  |      |
          if (i.gt.0)  goto 19  !----------------------+------+-ERROR+
          goto 5   !->-------------------------->------+------+      v
        elseif (k .eq. -1)  Then  !--------------------+      |
          txt(19) = 'duplicate keyword "period="'     !|      |
          goto 2 !------------------------------->-----+------+-ERROR+
        endif  !---------------------------------------+      |      v
                                                             !|
                                                             !|
        i = Index(buff(1:ll),'end period')                   !|
        j = Index(buff(1:ll),'endperiod')                    !|
        if (i .gt. 0  .or.  j .gt. 0)  Then  !---------+      |
          if (i .gt. 0) Then !--+                      |      |
            buff(i:i+9) = ' '  !|                      |      |
          else  !---------------+                      |      |
            buff(j:j+8) = ' '  !|                      |      |
          endif  !--------------+                      |      |
          ncycle = ncycle-1                           !|      ^
          if (ncycle .lt. 0) goto 7 !------->----------+------+-ERROR+
          lena = lyr - nbegin(ncycle+1) + 1           !|      |      v
          lyr  = lyr + lena*(nperiod(ncycle+1)-1)     !|      |
          if (lyr .gt. N_Top_Max) goto 3 !-->----------+------+-ERROR+
          do    i=2,nperiod(ncycle+1)  !======+        |      ^      v
          do    j=1,lena   !================+ |        |      |
            l1 = nbegin(ncycle+1)+j-1      !| |        |      |
            l2 = l1 + lena*(i-1)           !| |        |      |
            Thickness(l2) = Thickness(l1)  !| |        |      ^
            TTT = TTT + Thickness(l1)      !| |        |      |
            N_Frac(l2)    = N_Frac(l1)     !| |        |      |
            do m=1,N_Frac_Max !=========+   | |        |      |
              Code(m,l2)  = Code(m,l1) !|   | |        |      |
              Frac(m,l2)  = Frac(m,l1) !|   | |        |      |
            enddo  !====================+   | |        |      |
            x0(l2)        = x0(l1)         !| |        |      |
            xh(l2)        = xh(l1)         !| |        |      ^
            xf(l2)        = xf(l1)         !| |        |      |
            W0(l2)        = W0(l1)         !| |        |      |
            Wh(l2)        = Wh(l1)         !| |        |      |
            xs(l2)        = xs(l1)         !| |        |      |
            Sigma(l2)     = Sigma(l1)      !| |        |      |
            Daa(l2)       = Daa(l1)        !| |        |      ^
            C0h(l2)       = C0h(l1)        !| |        |      |
          enddo  !==========================+ |        |      |
          enddo  !============================+        |      |
          Call txtShift (buff(1:ll),j,i)              !|      ^
c -- the call returns: i=Len_Trim(vm)                  |      |
          if (i.gt.0)  goto 19  !----------------------+------+-ERROR+
          goto 5   !->------------------------>--------+------+      v
        elseif (i .eq. -1  .or.  j .eq. -1) Then !-----+      |
          txt(19) = 'duplicate keyword "end period"'  !|      |
          goto 2 !------------------------------>------+------+-ERROR+
        endif  !---------------------------------------+      |      v
                                                             !|
                                                             !^
c Now, since the line does not describe a period,             |
c it should describe a layer:                                 |
        lyr = lyr + 1                                        !|
        if (lyr .gt. N_Top_Max) goto 3   !-->-----------------+-ERROR+
        Thickness(lyr) = 0.                                  !|      v
        N_Frac(lyr)    = 1                                   !|
        do m=1,N_Frac_Max  !===+                              |
          Code(m,lyr)  = ' '  !|                              |
          Frac(m,lyr)  = 0.   !|                              |
        enddo  !===============+                              |
        Frac(1,lyr)    = 1.                                  !|
        x0(lyr)        = (0.,0.)                             !^
        xh(lyr)        = (0.,0.)                             !|
        xf(lyr)        = 0.                                  !|
        W0(lyr)        = 1.                                  !|
        Wh(lyr)        = 1.                                  !|
        xs(lyr)        = 0.                                  !|
        Sigma(lyr)     = 0.                                  !|
        Daa(lyr)       = 0.                                  !^
        C0h(lyr)       = ' '                                 !|
                                                             !|
                                                             !|
        Call ArgSearch (buff(1:ll),'t=',' ',',',k,i,j)       !|
        if (k .gt. 0)  Then  !----------------------------+   |
          txt(19) = 'thickness'                          !|   |
          Call rdReal (Tmp,1,buff(i:j),m)                !|   |
          if (m .ne. 0)  goto 2 !-------------------------+->-+-ERROR+
          Thickness(lyr) = Tmp(1)                        !|   |      v
          buff(k:j) = ' '                                !|   |
          txt(19) = 'zero layer thickness "t" - '//      !|   |
     *              'please comment out the line!'       !|   |
          if (Abs(Thickness(lyr)) .lt. 1.E-32) goto 2 !---+->-+-ERROR+
          txt(19) = 'negative layer thickness "t"'       !|   |      v
          if (Thickness(lyr) .lt. 0.) goto 2 !------------+->-+-ERROR+
          TTT = TTT + Thickness(lyr)                     !|   |      v
          t_found = 1                                    !|   |
        elseif (k .eq. -1)  Then  !-----------------------+   |
          txt(19) = 'duplicate keyword "t="'             !|   |
          goto 2 !----------------------------------------+->-+-ERROR+
        elseif (k .eq. 0)   Then  !-----------------------+   |
          t_found = 0  !unspecified layer thickness "t"   |   |
        endif  !------------------------------------------+   |
                                                             !|
                                                             !|
        Call ArgSearch (buff(1:ll),'w0=',' ',',',k,i,j)      !|
        if (k .gt. 0)  Then  !-------------------------+      |
          txt(19) = 'w0'                              !|      |
          Call rdReal (Tmp,1,buff(i:j),m)             !|      |
          if (m .ne. 0)         goto 2 !---------------+------+-ERROR+
          W0(lyr) = Tmp(1)                            !|      |      v
          txt(19) = 'negative "w0"'                   !|      |
          if (W0(lyr) .lt. 0.)  goto 2 !---------------+------+-ERROR+
          txt(19) = 'too big "w0" (max=99)'           !|      |      v
          if (W0(lyr) .gt. 99.) goto 2 !---------------+------+-ERROR+
          buff(k:j) = ' '                             !|      |      v
        elseif (k .eq. -1)  Then  !--------------------+      |
          txt(19) = 'duplicate keyword "w0="'         !|      |
          goto 2 !------------------------------->-----+------+-ERROR+
        endif  !---------------------------------------+      |      v
                                                             !|
                                                             !|
        Call ArgSearch (buff(1:ll),'wh=',' ',',',k,i,j)      !|
        if (k .gt. 0)  Then  !-------------------------+      |
          txt(19) = 'wh'                              !|      |
          Call rdReal (Tmp,1,buff(i:j),m)             !|      |
          if (m .ne. 0)         goto 2 !---------------+------+-ERROR+
          Wh(lyr) = Tmp(1)                            !|      |      v
          txt(19) = 'negative "wh"'                   !|      |
          if (Wh(lyr) .lt. 0.)  goto 2 !---------------+------+-ERROR+
          txt(19) = 'too big "wh" (max=99)'           !|      |      v
          if (Wh(lyr) .gt. 99.) goto 2 !---------------+------+-ERROR+
          buff(k:j) = ' '                             !|      |      v
        elseif (k .eq. -1)  Then  !--------------------+      |
          txt(19) = 'duplicate keyword "wh="'         !|      |
          goto 2 !------------------------------->-----+------+-ERROR+
        endif  !---------------------------------------+      |      v
                                                             !|
c Added in M7:                                               !|
        Call ArgSearch (buff(1:ll),'xhphase=',' ',',',k,i,j) !|
        if (k .gt. 0)  Then  !-------------------------+      |
          txt(19) = 'xhphase'                         !|      |
          Call rdReal (Tmp,1,buff(i:j),m)             !|      |
          if (m .ne. 0)         goto 2 !---------------+------+-ERROR+
          xs(lyr) = Tmp(1)                            !|      |      v
          txt(19) = '"xhphase" not in range [-2:2]'   !|      |
          if (xs(lyr) .lt. -2. .OR.                   !|      |
     *        xs(lyr) .gt. +2.) goto 2 !---------------+------+-ERROR+
          buff(k:j) = ' '                             !|      |      v
        elseif (k .eq. -1)  Then  !--------------------+      |
          txt(19) = 'duplicate keyword "xhphase="'    !|      |
          goto 2 !------------------------------->-----+------+-ERROR+
        endif  !---------------------------------------+      |      v
                                                             !|
                                                             !^
        Call ArgSearch (buff(1:ll),'sigma=',' ',',',k,i,j)   !|
        if (k .gt. 0)  Then  !---------------------------+    |
          txt(19) = 'sigma'                             !|    |
          Call rdReal (Tmp,1,buff(i:j),m)               !|    |
          if (m .ne. 0)                   goto 2 !-------+----+-ERROR+
          Sigma(lyr) = Tmp(1)                           !|    |      v
          txt(19) = 'negative rms roughness "sigma"'    !|    |
          if (Sigma(lyr) .lt. 0.)         goto 2 !-------+----+-ERROR+
          txt(19) = 'rms roughness "sigma" '//          !|    |      v
     *                 'exceeds layer thickness "t"'    !|    |
          if (Sigma(lyr) .gt. Thickness(lyr) .AND.      !|    |
     *        t_found .ne. 0)             goto 2 !-------+----+-ERROR+
          buff(k:j) = ' '                               !|    |      v
        elseif (k .eq. -1)  Then  !----------------------+    |
          txt(19) = 'duplicate keyword "sigma="'        !|    |
          goto 2 !------------------------------->-------+----+-ERROR+
        endif  !-----------------------------------------+    |      v
                                                             !|
                                                             !|
        Call ArgSearch (buff(1:ll),'da/a=',' ',',',k,i,j)    !|
        if (k .gt. 0)  Then  !---------------------------+    |
          if (buff(j:j) .eq. ',') Then !--+              |    |
            buff(j:j) = ' '              !|              |    |
            if (j .gt. i) j=j-1           !|             |    |
          endif  !-------------------------+             |    |
          if (buff(i:j) .eq. 'a' .or.                   !|    |
     *        buff(i:j) .eq. 'au' .or.                  !|    |
     *        buff(i:j) .eq. 'aut' .or.                 !|    |
     *        buff(i:j) .eq. 'auto') Then  !-------+     |    |
            Daa(lyr) = 1.E+30                     !|     |    |
          else  !----------------------------------+     |    |
            txt(19) = 'da/a'                      !|     |    |
            Call rdReal (Tmp,1,buff(i:j),m)       !|     |    |
            if (m .ne. 0)               goto 2 !---+-----+----+-ERROR+
            Daa(lyr) = Tmp(1)                     !|     |    |      v
            txt(19) = 'too large |da/a|>0.5'      !|     |    |
            if (Abs(Daa(lyr)) .gt. 0.5) goto 2 !---+-----+----+-ERROR+
          endif  !---------------------------------+     |    |      v
          buff(k:j) = ' '                               !|    |
        elseif (k .eq. -1)  Then  !----------------------+    |
          txt(19) = 'duplicate keyword "da/a="'         !|    |
          goto 2 !------------------------------->-------+----+-ERROR+
        endif  !-----------------------------------------+    |      v
                                                             !|
                                                             !^
        Call ChiSearch (buff(1:ll), 'x0=', x0(lyr),          !|
     *                  C0h(lyr)(1:1), txt(19), '0', *2)     !|
                                                             !|
        Call ChiSearch (buff(1:ll), 'xh=', xh(lyr),          !|
     *                  C0h(lyr)(2:2), txt(19), 'h', *2)     !|
                                                             !|
        Call ArgSearch (buff(1:ll),'xhdf=',' ',',',k,i,j)    !|
        if (k .gt. 0)  Then  !---------------------------+    |
          txt(19) = 'Xrh/Xih phase (xhdf=)'             !|    |
          Call rdReal (Tmp,1,buff(i:j),m)               !|    |
          if (m .ne. 0)             goto 2 !-------------+----+-ERROR+
          xf(lyr) = Tmp(1)                              !|    |      v
          txt(19) = 'xhdf outside [-1,1] range'         !|    |
          if (Abs(xf(lyr)) .gt. 1.) goto 2 !-------------+----+-ERROR+
          buff(k:j) = ' '                               !|    |      v
          txt(19) = 'xhdf data is only valid with '//   !|    |
     *                'implicit specification of xh#0'  !|    |
          if (Abs(xh(lyr)) .lt. 1.E-32) goto 2 !---------+----+-ERROR+
        elseif (k .eq. -1)  Then  !----------------------+    |      v
          txt(19) = 'duplicate keyword "xhdf="'         !|    |
          goto 2 !----------------->---------------------+----+-ERROR+
        endif  !-----------------------------------------+    |      v
                                                             !^
        Frac_Sum  = 0.                                       !|
        izero_frac= 0                                        !|
        mzero_frac= 0                                        !|
        do m=1,N_Frac_Max  !=============================+    |
                                                        !|    |
          txt(18) = 'codeN='                            !|    |
          write(txt(18)(5:5),'(i1)') m                  !|    |
          ltx = 6                                       !|    |
c Search for "codeN=":                                   |    |
          Call ArgSearch (buff(1:ll),txt(18)(1:ltx),    !|    |
     *                               ' ',',',k,i,j)     !|    |
          if (k .eq. 0 .AND. m .eq. 1) Then  !-----+     |    |
c Retry with "code=" instead of "code1="           |     |    |
            Call ArgSearch (buff(1:ll),'code=',   !|     |    |
     *                             ' ',',',k,i,j) !|     |    |
          endif  !---------------------------------+     |    |
          if (k .gt. 0) Then  !-----------------------+  |    |
                                                     !|  |    ^
c One more fraction code found:                       |  |    |
            if (m .gt. 1) N_Frac(lyr)=N_Frac(lyr)+1  !|  |    |
                                                     !|  |    |
            if (buff(j:j) .eq. ',') Then !--+         |  |    |
              buff(j:j)           = ' '    !|         |  |    |
              buff(j+2000:j+2000) = ' '    !|         |  |    |
              if (j .gt. i) j = j-1        !|         |  |    |
            endif  !------------------------+         |  |    |
            Code(m,lyr) = buff(i+2000:j+2000)        !|  |    |
            i = Len_Trim(Code(m,lyr))                !|  |    |
            if ((Code(m,lyr)(1:1).eq.char(34) .and.  !|  |    |
     *           Code(m,lyr)(i:i).eq.char(34)) .or.  !|  |    ^
     *          (Code(m,lyr)(1:1).eq.char(39) .and.  !|  |    |
     *           Code(m,lyr)(i:i).eq.char(39)))      !|  |    |
     *      Code(m,lyr) = Code(m,lyr)(2:i-1)         !|  |    |
            buff(k:j) =  ' '                         !|  |    |
                                                     !|  |    |
          elseif (k .eq. -1)  Then  !-----------------+  |    |
                                                     !|  |    |
            write (txt(19),119) txt(18)(1:ltx)       !|  |    |
  119       format('duplicate keyword "',a,'"')      !|  |    |
            goto 2 !------------------------------->--+--+----+-ERROR+
                                                     !|  |    |      v
          else  !-------------------------------------+  |    ^
                                                     !|  |    |
c Goto exit from the loop:                           !|  |    |
            if (m .gt. 1)  goto 71 !------------------+--+--+ |
c (if Code-1 is missing, the substrate code           |  |  | |
c  code will be used. Here we simply continue)        |  |  | |
                                                     !|  |  | |
          endif  !------------------------------------+  |  | |
                                                        !|  | |
          txt(18) = 'xN='                               !|  | |
          write(txt(18)(2:2),'(i1)') m                  !|  | |
          ltx = 3                                       !|  | |
c Search for "xN=":                                     !|  | |
          Call ArgSearch (buff(1:ll),txt(18)(1:ltx),    !|  | |
     *                               ' ',',',k,i,j)     !|  | |
          if (k .eq. 0 .AND. m .eq. 1) Then  !-----+     |  | |
c Retry with "x=" instead of "x1="                 |     |  | |
            Call ArgSearch (buff(1:ll),'x=',      !|     |  | |
     *                             ' ',',',k,i,j) !|     |  | |
          endif  !---------------------------------+     |  | |
          if (k .gt. 0)  Then !---------------------+    |  | ^
                                                   !|    |  | |
            txt(19) = 'fraction "'//               !|    |  | |
     *                txt(18)(1:ltx-1)//'"'        !|    |  | |
            Call rdReal (Tmp,1,buff(i:j),merr)     !|    |  | |
            if (merr .ne. 0)            goto 2 !----+----+--+-+-ERROR+
            Frac(m,lyr) = Tmp(1)                   !|    |  | |      v
            txt(19) = 'negative fraction "'//      !|    |  | |
     *                txt(18)(1:ltx-1)//'"'        !|    |  | |
            if (Frac(m,lyr) .lt. 0.)    goto 2 !----+----+--+-+-ERROR+
            txt(19) = 'fraction "'//               !|    |  | |      v
     *                txt(18)(1:ltx-1)//           !|    |  | |
     *                '" exceeds 1'                !|    |  | |
            if (Frac(m,lyr) .gt. 1.)    goto 2 !----+----+--+-+-ERROR+
            buff(k:j) = ' '                        !|    |  | |      v
                                                   !|    |  v ^
            Frac_Sum = Frac_Sum + Frac(m,lyr)      !|    |  | |
            txt(19) = 'sum of fractions exceeds'// !|    |  | |
     *                ' 1 or negative fraction xN' !|    |  | |
            if (Frac_Sum .gt. 1.001)    goto 2 !----+----+--+-+-ERROR+
                                                   !|    |  | |      v
          elseif (k .eq. -1)  Then  !---------------+    |  | |
                                                   !|    |  | |
            write (txt(19),119) txt(18)(1:ltx)     !|    |  | |
            goto 2 !----------------------------->--+----+--+-+-ERROR+
                                                   !|    |  | |      v
          else  !-----------------------------------+    |  v ^
                                                   !|    |  | |
c Fraction xN is not specified. This is allowed    !|    |  | |
c only once, otherwice we don't know the sum:      !|    |  | |
            izero_frac = izero_frac+1              !|    |  | |
            mzero_frac = m                         !|    |  | |
            txt(19) = 'more than one unspeci'//    !|    |  | |
     *                'fied fraction xN'           !|    |  | |
            if (izero_frac .gt. 1)      goto 2 !----+----+--+-+-ERROR+
                                                   !|    |  | |      v
          endif  !----------------------------------+    |  | |
        enddo  !=========================================+  | |
  71    continue   !<---------------------------------------+ |
        if (izero_frac .eq. 1) Then !--------------------+    |
          Frac(mzero_frac,lyr) = 1. - Frac_Sum          !|    |
          txt(19) = 'fraction x? is zero or negative'   !|    |
          write (txt(19)(11:11),'(i1)') mzero_frac      !|    |
          if (Frac(mzero_frac,lyr) .le. 0.) goto 2 !-----+----+-ERROR+
        else  !------------------------------------------+    |      v
          txt(19) = 'sum of fractions is not 1'         !|   !|
          if (Abs(Frac_Sum-1.) .gt. 0.001)  goto 2 !-----+----+-ERROR+
        endif  !-----------------------------------------+   !|      v
                                                             !|  |E
c Test for unknown keywords (all the found keywords           |  vN
c and respecive values are erased):                           |  |D
        Call txtShift (buff(1:ll),j,i)      !i=Len_Trim(vm)   |  |
        if (i .gt. 0)  goto 9  !--(unknown keyword)-----------+--+-ERROR+
c Test if thickness was specified:                            |  |      v
        if (t_found .eq. 0)   Then  !------------------+      |  |
          txt(19) = 'unspecified layer thickness "t"' !|      |  |
          goto 2 !------------------------------->-----+------+--+-ERROR+
        endif !----------------------------------------+      |  |      v
c Read next profile line:                                     |  |
        goto 5   !->------------------------------------------+  |
  4     continue    !<------------------<------------------------+
        if (ncycle .ne. 0)  goto 8 !-(loops not closed)--------ERROR+
        Close   (unit=lun)                                         !v
        goto 30 !-------------------------------------------------------+
                                                                       !|
c----------------------------                                           |
  1     continue                                                       !|
        write   (txt,21) progname(1:lpn),file(1:lfile)                 !|
  21    format  (a,' cannot open top layer profile:'//a///)            !|
        goto 29  !----------------------------------------------------+ |
                                                                     !| |
  2     continue                                                     !| |
        i = Len_Trim(txt(19))                                        !| |
        if (i .lt. 1)  i=1                                           !| |
        if (i .gt. 67) i=67                                          !| |
        buff(1:i) = txt(19)(1:i)                                     !| |
        j = Len_Trim(buff(2001:4000))                                !| |
        if (j .lt. 1)  j=1                                           !v |
        if (j .gt. 70) j=70                                          !| |
        write   (txt,22) progname(1:lpn),file(1:lfile),              !| |
     *                   buff(1:i),line,buff(2000+1:2000+j)        !v v
  22    format  (a,' has problem with top layer profile:'/a/
     *           'Problem: [',a,']'/
     *           'Line Nr: [',i4,']'/
     *           'Line: [',a,']')                                    !| |
        goto 29   !--------------------->-----------------------------+ |
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
        write   (txt,23) progname(1:lpn),file(1:lfile),              !|v|
     *                   N_Top_Max                                   !v v
  23    format  (a,' has problem with top layer profile:'//a//
     *           'Buffer size (',i4,' layers) exceeded')             !| |
        goto 29   !--------------------->-----------------------------+ |
                                                                     !| |
  6     continue                                                     !| |
        write   (txt,26) progname(1:lpn),file(1:lfile),line          !v v
  26    format  (a,' has problem with top layer profile:'//a//
     *          'Attempt to start level-3 enclosed loop at line=',i4)!| |
        goto 29   !--------------------->-----------------------------+ |
                                                                     !| |
  7     continue                                                     !| |
        write   (txt,27) progname(1:lpn),file(1:lfile),line          !v v
  27    format  (a,' has problem with top layer profile:'//a//
     *          'Attempt to close non-existing loop at line=',i4)    !| |
        goto 29   !--------------------->-----------------------------+ |
                                                                     !| |
  8     continue                                                     !| |
        write   (txt,28) progname(1:lpn),file(1:lfile)               !v v
  28    format  (a,' has problem with top layer profile:'//a//
     *          'Loop(s) not closed at EOF')                         !| |
        goto 29   !--------------------->-----------------------------+ |
                                                                     !| |
  29    continue  !<--------------------------------------------------+ |
        nlines = 5                                                     !|
  69    continue                                                       !|
        Call    Message (txt,nlines,2)                                 !|
        ifail = 1                                                      !|
                                                                       !|
  30    continue !<-----------------------------------------------------+
        Thickness_Top = Sngl(TTT)
        N_Top         = lyr

  99    continue
        if (iirezv .eq. 1)  Then  !--------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end

c=======================================================

        Subroutine Read_SubM7_0h (Code, W0, Wh,
     *                            Daa, Sigma, xs,
     *                            wrk, line, *, *)

        Real*4          W0, Wh, Daa, Sigma, xs,
     *                  Tmp(2)
        Integer         lpn, i, j, n,
     *                  iirezv, line,
     *                  iexit
        Integer         N_Column
        External        N_Column
        Character       Code*(*), wrk*(*)

        Character       txt(20)*80
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
        if (istackrezv .lt. 10) Then !-----+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Read_SubM7_0h'      !|
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
        txt(7) = 'substrate Debye-Waller factors "w0","wh" '//
     *                                           '(range: [0-99], w0>0)'
        W0 = 1.
        Wh = 1.
        Call    LineLoop (2,line,wrk,*100)
        n = n_Column (wrk)
        if (n .gt. 0)  Then  !-----------------+
          Call  rdReal  (Tmp,2,wrk,i)         !|
          if (i .ne. 0)     goto 100          !|
          W0 = Tmp(1)                         !|
          if (W0 .le. 0.)   goto 101          !|
          if (W0 .gt. 99.)  goto 101          !|
          if (n .gt. 1) Then  !--------------+ |
            Wh = Tmp(2)                     !| |
            if (Wh .lt. 0.)  goto 101       !| |
            if (Wh .gt. 99.) goto 101       !| |
          endif  !---------------------------+ |
        endif  !-------------------------------+

c-----------------------
        txt(7) = 'substrate strain "da/a" (|da/a|<0.5)'
        Call    LineLoop (2,line,wrk,*100)
        Call    Case    (wrk,1) !to lower case
        if (wrk .eq. 'a'    .or.
     *      wrk .eq. 'au'   .or.
     *      wrk .eq. 'aut'  .or.
     *      wrk .eq. 'auto') Then  !------------------+
          Daa = 0.                                   !|
        else  !---------------------------------------+
          Call  rdReal  (Tmp,1,wrk,i)                !|
          if (i .ne. 0)          goto 100            !|
          Daa = Tmp(1)                               !|
          if (Abs(Daa) .gt. 0.5) goto 101            !|
        endif  !--------------------------------------+

c-----------------------
        txt(7) = 'substrate rms roughness "sigma" (sigma >=0)'
        Call    LineLoop (2,line,wrk,*100)
        Call    rdReal  (Tmp,1,wrk,i)
        if (i .ne. 0)           goto 100
        Sigma  = Tmp(1)
        if (sigma .lt. 0.)      goto 101

c-----------------------
        txt(7) = 'substrate phase shift in units of pi [-2:2]'
c       Call    LineLoop (2,line,wrk,*100)
c       Call    rdReal  (Tmp,1,wrk,i)
c       if (i .ne. 0)           goto 100
c       xs  = Tmp(1)
        xs  = 0.
        if (xs .lt. -2. .OR.
     *      xs .gt. 2.)         goto 101

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
        end
