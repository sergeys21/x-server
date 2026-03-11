        Character*8 Function getOS ()
        Character envvalue*32

        Logical*4 FileExist
        External  FileExist

c       Call GetEnv32 ('SHELL',envvalue)                !SHELL is not defined in CGI
        Call GetEnv32 ('PATH',envvalue)

        if (Len_Trim(envvalue) .gt. 0) Then !----------+if PATH exists/defined
           if (envvalue(1:1) .eq. '/') Then !-------+  |
              getOS = 'linux'                      !|  |
           else !-----------------------------------+  |
              getOS = 'windows'                    !|  |
           endif  !---------------------------------+  |
        else !-----------------------------------------+
           if (FileExist('/etc/passwd')) Then !-----+  |
              getOS = 'linux'                      !|  |
           else !-----------------------------------+  |
              getOS = 'windows'                    !|  |
           endif  !---------------------------------+  |
        endif !----------------------------------------+

        return
        end
