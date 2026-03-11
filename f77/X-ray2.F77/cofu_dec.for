        Program Cofu_DECODE
c This is a stand-alone utility which asks for h_jaggedness,
c and decodes for this h_jaggedness the correlation function
c integral encoded in the file <corrfun.cft>. The encoded
c intepgral G(p,h) is stored in <corrfun.dat>. The procedure
c Cofu_DECODE may be used for debugging purposes. The integral
c is used in trds_sl and gids_sl.
c+-----------------------------------------------------------+
c|       +inf                 | -- for p={0.,50.}            |
c|        ?    -x^{2h}        | The details see in S.K.Sinha |
c|G(p,h)= |dx e      cos(px)  |  -- J.Physique III (France)  |
c|        ?                   | v.4, p.1543-1557 (1994)      |
c|        0                   |  --- Eq.(25) on p.1552.      |
c+-----------------------------------------------------------+
        Integer*2  np, nh
        Integer    i, j, ih, io_status,
     *             lug, ndim, luo,
     *             istat, ifail
        Parameter  (ndim = 4001)
        Real*4     Array(ndim),    !Table of G_n(p) according to Sinha
     *             pmin, pmax,
     *             hmin, hmax,
     *             dp,             !(pmax-pmin)/(np-1)
     *             dh,             !(hmax-hmin)/(nh-1)
     *             p,h,x,q,Z

        lug = 1
        luo = 3
        Call OpenBinary ('corrfun.cft',lug,'old','read',istat,ifail)
        if (ifail.ne.0) goto 5

        read    (lug)   np,nh,
     *                  pmin,pmax,
     *                  hmin,hmax
        if (np.gt.ndim) Then !--------------------------------------+
          write (*,*) ' Number of points  in  file   np = ',np     !|
          write (*,*) ' Number of points in buffer ndim = ',ndim   !|
          stop 'np > ndim'                                         !|
        endif !-----------------------------------------------------+
        if (np .le. 1) stop 'np <=1'
        dp = (pmax-pmin)/(np-1)
        dh = (hmax-hmin)/(nh-1)

        Call OpenFile('corrfun.dat',luo,'write','unknown',io_status,*7)
  2     continue  !<---------------------------------+
        write (*,1) hmin,hmax                       !|
c 1     format (' Enter h [',f4.2,':',f4.2,']_'\)   !|  Compaq/MS fortran suppress carriage return
  1     format (' Enter h [',f4.2,':',f4.2,']_',$)  !|  Compaq/GNU fortran suppress carriage return
        read (*,*) h                                !|
        if (h.lt.hmin  .or. h.gt.hmax) goto 2   !----+

        x = (h-hmin)/dh + 1.
        ih = INT(x+0.1)
        if (ih.lt.1)    ih=1
        if (ih.gt.nh-1) ih=nh-1
        q  = x-ih

        write (luo,4) h
 4      format('# Data from corrfun.cft exported at h=',f5.3)
        do      i=1,ih  !========================+
          read (lug,err=6) (Array(j),j=1,np)    !|
        enddo  !=================================+
        do      j=1,np  !========================+
          read(lug,err=6,end=9) Z               !|
          Array(j) = Array(j) + q*(Z-Array(j))  !|
                                                !|
          p = pmin + dp*(j-1)                   !|
          write (luo,*,err=8) p, Array(j)       !|
        enddo  !=================================+

        stop 'Done!'
  5     continue
        stop 'Error opening corrfun.cft'
  6     continue
        stop 'Error reading corrfun.cft'
  7     continue
        stop 'Error opening corrfun.dat'
  8     continue
        stop 'Error writing corrfun.dat'
  9     continue
        stop 'End-of-file while reading corrfun.cft'
        end

