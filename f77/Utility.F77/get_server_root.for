        Character*64 Function Get_Server_Root ()
c Returns root of X-ray server:
        Integer         nch
        Character       root*64

c Specifying server and especially a protocol may lead to an insecure
c mixed content (e.g. browser will ignore CSS sheets as "insecure" if
c the page itself is loaded as https). Therefore it is better to use
c a folder-only specification.
        Call GetEnv32 ('XSERVER',root)
        if (Len_Trim(root).eq.0) Then !-----------------+ no environment
          root = '/'                                   !|
c         root = '/x-server/'                          !| when piggy-back on some server
        endif !-----------------------------------------+
        nch = Len_Trim(root)
        if (root(nch:nch).ne.'/' .and. nch.lt.Len(root)) Then !--+
          nch = nch+1                                           !|
          root(nch:nch) = '/'                                   !|must end with '/'
        endif  !-------------------------------------------------+
        Get_Server_Root = root
        return
        end
