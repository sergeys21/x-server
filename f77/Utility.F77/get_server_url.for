        Character*64 Function Get_Server_URL ()
c Returns root URL of X-ray server:
        Integer         nch
        Character       xserver*64

c Specifying server and especially a protocol may lead to an insecure
c mixed content (e.g. browser will ignore CSS sheets as "insecure" if
c the page itself is loaded as https). Therefore it is better to use
c a folder-only specification.
        Call GetEnv32 ('XSERVER_URL',xserver)
        if (Len_Trim(xserver).eq.0) Then !-----------------+ no environment
          xserver = 'https://x-server.gmca.aps.anl.gov/'  !|
        endif !--------------------------------------------+
        nch = Len_Trim(xserver)
        if (xserver(nch:nch).ne.'/' .and. nch.lt.Len(xserver)) Then !--+
          nch = nch+1                                                 !|
          xserver(nch:nch) = '/'                                      !|must end with '/'
        endif  !-------------------------------------------------------+
        Get_Server_URL = xserver
        return
        end
