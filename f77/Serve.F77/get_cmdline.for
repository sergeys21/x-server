        Subroutine Get_CmdLine (nkeys,          !input, number of acceptable keys
     *                          keys,           !input, acceptable keys
     *                          vals,           !output, specified keys values (character
     *                          ifail)          !output, status (at input serves as "strict" flag)
c--------------------------------------------------------
c             Analyzing command line arguments
c
c The command line is assumed to have the structure:
c
c Exefile key1=value1 ... keyN=valueN
c
c If key is specified, and value is given,       val=value
c If key is specified, and value is too long,    val="toolong"
c If key is specified, but value is missing,     value="none"
c If key is missing,                             value="missing"
c
c ATTENTION: here the keys are NOT case sensitive
c
c                    Author: S.Stepanov
c
c See also: Txt_Args (for command lines with single-character arguments)
c--------------------------------------------------------
        Integer*4       Nargs
        External        Nargs

        Integer         nkeys, ifail, narg, nfound, i, j, l, k, n, m,
     *                  lenkey, lenval, strict, lunerr

        Character       keys(nkeys)*(*),
     *                  vals(nkeys)*(*),
     *                  prg*16, key1*32, key2*32

        Character       txt(20)*80
        Common  /msg/   txt

c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c--------------------------------------------------------
        prg = 'Get_CmdLine'
        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = prg                  !|
        else !-----------------------------+
          iirezv = 0                      !|
        endif  !---------------------------+

        lunerr = 0
        strict = ifail
        ifail  = 0
        nfound = 0
        do      i=1,nkeys  !===+
         vals(i) = 'missing'  !|
        enddo  !===============+
        narg = Nargs()-1
        if (narg.lt.1) goto 4  !-------------return-------------+
                                                               !v
        lenkey = min(len(keys(1)),32)
        lenval = len(vals(1))
        do      i=1,narg  !=================================+
c This is the GNU Fortran form. Should also work            |
c for Compaq/MS fortran                                     |
          Call get_arg_wrapper (i,txt(20))                 !|
          l = Index(txt(20),'=')                           !|
          if (l .lt. 2)        goto 1 !---------ERROR-------+---+
          if (l-1 .gt. lenkey) goto 1 !---------ERROR-------+---+
          l = l-1                                          !|   v
          key1 = txt(20)(1:l)                              !|
          Call Case (key1(1:l),1)    !convert to lower case |
          do    j=1,nkeys  !==============================+ |   v
            k = Max(Len_Trim(keys(j)),1)                 !| |
            key2 = keys(j)(1:k)                          !| |
            Call Case (key2(1:k),1) !convert to lower case| |
            if (key1(1:l) .eq. key2(1:k)) Then !-------+  | |
              m = Len_Trim(txt(20))                   !|  | |
              if (m-l-1 .gt. lenval) Then !---------+  |  | |
                vals(j)='toolong'                  !|  |  | |
              elseif (m-l-1 .le. 0) Then !----------+  |  | |
                vals(j)='none'                     !|  |  | |
              else !--------------------------------+  |  | |
                vals(j) = txt(20)(l+2:m)           !|  |  | |
                m = m-l-1                          !|  |  | |
c Remove quotes, if present                         |  |  | |
                if (vals(j)(1:1).eq.'"' .AND.      !|  |  | |
     *              vals(j)(m:m).eq.'"') Then !--+  |  |  | |
                  vals(1:1) = ' '               !|  |  |  | |
                  vals(m:m) = ' '               !|  |  |  | |
                  Call TxtShift(vals(1:m),n,m)  !|  |  |  | |
                endif !--------------------------+  |  |  | |
              endif !-------------------------------+  |  | |
              nfound = nfound + 1                     !|  | |
              goto 7 !--------+                        |  | |
            endif  !----------+------------------------+  | |
          enddo  !============+===========================+ |
c Unknown key in cmdline:     |                             |
          goto 2  !-----------+-----------------ERROR-------+--+
  7       continue  !<--------+                             |  v
        enddo  !============================================+

        if (nfound .eq. nkeys .OR. strict .ne. 1) goto 555  !--return--+
c Not all keys are specified. Find at least one not-specified key      |
        do    j=1,nkeys  !=====================================+       v
          k = Max(Len_Trim(keys(j)),1)                        !|
          key2 = keys(j)(1:k)                                 !|
          Call Case (key2(1:k),1)    !convert to lower case    |
          do i=1,narg  !===================================+   |
            Call get_arg_wrapper (i,txt(20))              !|   |
            l = Index(txt(20),'=')                        !|   |
            l = l-1                                       !|   |
            key1 = txt(20)(1:l)                           !|   |
            Call Case (key1(1:l),1)  !convert to lower case|   |
            if (key1(1:l) .eq. key2(1:k)) goto 10 !--------+-+ |
          enddo !==========================================+ | |
c keys(j) is not present in cmdline (this is strict          | |
c checking only; the value is already set to "missing")      | |
          goto 3 !------------------------------ERROR--------+-+--+
  10      continue !<----------------------------------------+ |  v
        enddo !================================================+
c We should not be here:
        goto 555  !--------------------------------------------return--+
c=========================================================             v
c                       E R R O R S :
c---------------------------------------------------------
  1     continue
        m = Len_Trim(txt(20))
        write   (txt,101)  prg(1:len_Trim(prg)), i, txt(20)(1:m)
  101   format  (a,': incorrect cmdline arg No.',i2/'val=[',a,']')
        n = 2
        ifail = 1
        goto 999

  2     continue
        write   (txt,102)  prg(1:len_Trim(prg)), txt(20)(1:l), i,
     *                     txt(20)(1:m)
  102   format  (a,': unknown cmdline key=[',a,'] in arg No.',i2/
     *           'val=[',a,']')
        n = 2
        ifail = 2
        goto 999

  3     continue
        write   (txt,103)  prg(1:len_Trim(prg)), keys(j)(1:k)
  103   format  (a,': cmdline key=[',a,'] is not specified.')
        n = 1
        ifail = 3
        goto 999

  4     continue
        l = max(len_Trim(progrezv(istackrezv)),1)
        write (lunerr,104) progrezv(istackrezv)(1:l),
     *                     (keys(i)(1:len_Trim(keys(i))),i=1,nkeys)
  104   format(1x,a,': no cmdline specified. Available keys:',100(1x,a))
        write   (txt,105)  progrezv(istackrezv)(1:l)
  105   format  (a,': no cmdline specified.')
        n = 1
        ifail = 4
        goto 999

  999   continue
        Call    Message (txt,n,2)
        goto 555

  555   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end

c=================================================================

        Subroutine getbykey(nkeys,keys,vals,key,kval)
c Finds kval by key in the arrays (keys,vals)
        Integer         nkeys, i, l, m

        Character       keys(nkeys)*(*),
     *                  vals(nkeys)*(*),
     *                  key*(*), kval*(*)

        m = Len_Trim(key)
        if (m .eq. 0) return
        do i=1,nkeys  !=============================================+
          l = Len_Trim(keys(i))                                    !|
          if (l .gt. 0 .AND. key(1:m) .eq. keys(i)(1:l)) Then !--+  |
             kval = vals(i)                                     !|  |
             if (kval .eq. "toolong") kval = ' '                !|  |
             if (kval .eq. "none")    kval = ' '                !|  |
             if (kval .eq. "missing") kval = ' '                !|  |
             return                                             !|  |
          endif !------------------------------------------------+  |
        enddo !=====================================================+
        kval = ' '
        return
        end
