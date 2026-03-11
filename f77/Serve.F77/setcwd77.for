        Integer*4 Function SetCWD (path)
c------------------------------------------------
c Set CWD to path.
c The 32-bit version (GNU g77/gfortran)
c
c                 By S.Stepanov
c------------------------------------------------
        Integer         ifail, l
        Character       path*(*)

        ifail = 1
        l = Len_trim (path)
        if (l.eq.0)     goto 1

        call chdir(path(1:l),ifail)             !GNU fortran function

  1     continue
        SetCWD = ifail
        return
        end

c======================================================================
c
c GetCWD does not require any additional function
c       i4result = GetCWD(dirname)
