        Subroutine Txt_Args (TextData,Flags,ninps,ifail)
c--------------------------------------------------------
c             Analyzing command line arguments
c
c The command line is assumed to have the structure:
c
c Exefile <Key1=<option1>> ... <KeyN=<optionN>>
c
c where /Key is a text constant consisting of one letter.
c If /option is missed it is assumed to have blank value.
c
c                    Author: S.Stepanov
c--------------------------------------------------------
c Example:
c           Exefile A=zzz B=rrr C:qqq
c
c See also: GetInput (for command lines with FILE arguments)
c--------------------------------------------------------
        Integer*4       Nargs
        External        Nargs

        Integer         ninps, ifail, narg, ifa, i, j, l

        Character       TextData(ninps)*(*),
     *                  Flags(ninps)*1, buf*1

        Character       txt(20)*80
        Common  /msg/   txt

c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv, iirezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile
c--------------------------------------------------------
        if (istackrezv.lt.10) Then  !------+
          iirezv = 1                      !|
          istackrezv = istackrezv+1       !|
          progrezv(istackrezv) = progname !|
          progname = 'Txt_Args'           !|
        else !-----------------------------|
          iirezv = 0                      !|
        endif  !---------------------------+

        ifail = 0
        do      i=1,ninps  !===+
         TextData(i) = ' '    !|
        enddo  !===============+
        narg = Nargs()-1
        if (narg.lt.1)  goto 555  !----------return-------------+
                                                               !v
        ifa = 0
        do      i=1,narg  !==============================+
c This is the GNU Fortran form. Should also work         |
c for Compaq/MS fortran                                  |
          Call get_arg_wrapper (i,txt(20))              !|
c Convert to UPPER case:                                !|
          Call  Case (txt(20)(1:1),0)                   !|
          do    j=1,ninps  !===========================+ |
            buf = Flags(j)                            !| |
c Convert to UPPER case:                              !| |
            Call Case (buf,0)                         !| |
            if ((txt(20)(1:2).eq.buf//'=') .OR.       !| |
     *          (txt(20)(1:2).eq.buf//':')) Then  !-+  | |
              l = Max(Len_Trim(txt(20)),1)         !|  | |
              TextData(j) = txt(20)(3:l)           !|  | |
              goto 3  !-------+                    !|  | |
            endif  !----------+---------------------+  | |
          enddo  !============+========================+ |
c Unknown parameter:          |                          |
          goto 2  !-----------+--------------------------+--+
  3       continue  !<--------+                          |  v
        enddo  !=========================================+

        goto 555  !--------------------------return-------------+
c=========================================================      v
c                       E R R O R S :
c---------------------------------------------------------
c ERRORS:
  2     continue
        l = Len_Trim(txt(20))
        write   (txt,102)  txt(20)(1:1),i,txt(20)(1:l)
  102   format  (/'Txt_Args: unknown parameter "',a,'"'//
     *  'in arg',i3,' "',a,'"'/)
        ifail = 2
        goto 999

  999   continue
        Call    Message (txt,5,2)
        goto 555

  555   continue
        if (iirezv.eq.1)  Then  !----------+
          progname = progrezv(istackrezv) !|
          istackrezv = istackrezv-1       !|
        endif  !---------------------------+
        return
        end
