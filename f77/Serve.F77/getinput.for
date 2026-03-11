        Subroutine GetInput (File,              !output, specified file
     *                       ivalues,           !output, specified keys values
     *                       textkeys,          !input, acceptable keys
     *                       nkeys,             !input, number of acceptable keys
     *                       ifail)             !output, status
c--------------------------------------------------------
c             Analyzing command line arguments
c
c The command line is assumed to have the structure:
c
c Exefile <Inpfilename> </Key1<ivalue1>> ... </KeyN<ivalueN>>
c
c where /Key is a text constant consisting of one letter.
c
c If /Key is specified, and <ivalue> is given,   ivalue=value
c If /Key is specified, but <ivalue> is missing, ivalue=-1
c If /Key is missing,                            ivalue=-2
c
c                    Author: S.Stepanov
c--------------------------------------------------------
c Example:
c           Exefile wwww.dat /a13 /b25 -c888
c
c See also: Txt_Args (for command lines without FILE arguments)
c--------------------------------------------------------
        Integer*4       Nargs
        External        Nargs

        Integer         nkeys, ivalues(nkeys), lines,
     *                  narg, ifail, i, j, k, l, ifa

        Character       File*(*), textkeys(nkeys)*1, buf*1, kk*2

        Character       txt(20)*80
        Common  /msg/   txt

c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c -------------------------------------------------------
        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'GetInput'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+

        ifail    = 0
        narg     = Nargs()-1
        if (narg.lt.1)  Then  !----+
c All /Key are missing:            |
          do      i=1,nkeys  !==+  |
           ivalues(i) = -2     !|  |
          enddo  !==============+  |
          goto 555                !|return
        endif  !-------------------+

        do      i=1,nkeys  !====+
         ivalues(i) = -32766   !|
        enddo  !================+
        ifa = 0
        do      i=1,narg  !====================================+
          Call get_arg_wrapper (i,txt(20))                    !|
c Ignore everything after redirection (e.g. > /dev/null 2>&1)  |
          if (txt(20)(1:1).eq.'>') goto 20      !break the loop|
          kk = txt(20)(1:2)                                   !|
          Call  Case (kk,0)                                   !|
          if (kk(1:1).eq.'/' .OR. kk(1:1).eq.'-') Then !----+  |
c This is a key?                                            |  |
            do    j=1,nkeys  !===========================+  |  |
              buf = textkeys(j)                         !|  |  |
              Call Case (buf,0)                         !|  |  |
              if (kk(2:2).eq.buf) Then  !--------------+ |  |  |
c If option is already specified:                      | |  |  |
                if (ivalues(j).ne.-32766) goto 2  !----+-+--+--+-----+2
                l = Len_Trim (txt(20))                !| |  |  |     v
                if (l.lt.3)       Then  !-----------+  | |  |  |
c If <option> is missing:                           |  | |  |  |
                  ivalues(j) = -1                  !|  | |  |  |
                else   !----------------------------+  | |  |  |
                  Call RdInt (ivalues(j),1,        !|  | |  |  |
     *                        txt(20)(3:l),k)      !|  | |  |  |
                  if (k.ne.0)       goto 3   !------+--+-+--+--+----+3
                endif  !----------------------------+  | |  |  |    v
                goto 10  !---------------------+       | |  |  |
              endif  !-------------------------+-------+ |  |  |
            enddo  !===========================+=========+  |  |
c The key is not recognized:                   |           !|  |
            goto 4 !---------------------------+------------+--+---+4
  10        continue    !<---------------------+            |  |   v
                                                           !|  |
          else !--------------------------------------------|  |
                                                           !|  |
c We assume it is a file:                                   |  |
            if (ifa.ne.0)         goto 5  !-----------------+--+-+5
            File = txt(20)                                 !|  | v
            ifa = 1                                        !|  |
          endif !-------------------------------------------+  |
        enddo  !===============================================+
  20    continue
c If /Key is missing:
        do      i=1,nkeys  !==========================+
         if (ivalues(i).eq.-32766) ivalues(i) = -2   !|
        enddo  !======================================+

        goto 555  !--------------------------return-------------+
c=========================================================      v
c                       E R R O R S :
c---------------------------------------------------------
  2     continue
        write   (txt,102)  kk,i
  102   format  (/'GetInput: parameter "',a,'"'/
     *  'is redefined in command line argument Nr ',i2/'CMD args:')
        lines = 4
        ifail = 2
        goto 999

  3     continue
        write   (txt,103)  kk,i
  103   format  (/'GetInput: error reading parameter "',a,'"'/
     *  'in command line argument Nr ',i2/'CMD args:')
        lines = 4
        ifail = 3
        goto 999

  4     continue
        write   (txt,104)  kk,i
  104   format  (/'GetInput: unknown key "',a,'"'/
     *  'in command line argument Nr ',i2/'CMD args:')
        lines = 4
        ifail = 3
        goto 999

  5     continue
        write   (txt,105)  i
  105   format  (/'GetInput: input filename is redefined'/
     *  'in command line argument Nr ',i2/'CMD args:')
        lines = 4
        ifail = 5
        goto 999

  999   continue
        l = min(narg,20-lines-2)
        do      i=1,l  !=======================+
           Call get_arg_wrapper (i,txt(20))   !|
           lines = lines + 1                  !|
           txt(lines) = txt(20)               !|
        enddo !================================+
        lines = lines + 1
        txt(lines) = ' '

        Call    Message (txt,lines,2)
        goto 555

  555   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
