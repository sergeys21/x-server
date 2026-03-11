        Subroutine GetEnv32 (envname,envvalue)
c This is a wrapper to GNU fortran GetEnv
c
c ATTENTION: Microsoft Fortran Power Station has its own GetEnv
c subroutine, but it returns garbage in the case when "envname" is
c not found. That is why for Compaq/MS fortran getenv is emulated by
c GetEnvQQ which works properly.
        Character envname*(*),envvalue*(*)

        Call GetEnv (envname,envvalue)
        return
        end
