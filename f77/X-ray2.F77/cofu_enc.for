        Program Cofu_ENCODE
c This is a stand-alone program which calculates the integral
c below for trds_sl and gids_sl.
c+-----------------------------------------------------------+
c|Tabulating the integral:    |                              |
c|                            |                              |
c|     +inf                   | -- for p={0.,50.}            |
c|      |    -x^{2h}          | The details see in S.K.Sinha |
c|G(p)= |dx e      cos(px)    |  -- J.Physique III (France)  |
c|      |                     | v.4, p.1543-1557 (1994)      |
c|      0                     |  --- Eq.(25) on p.1552.      |
c|                            |                              |
c+-----------------------------------------------------------+
c ATTENTION: The upper estimation for p is: p=qx*L.
c            For qx=0.004 and L=10,000~A  we get: p=40. -->--+
c                                                            |
        Integer*2  np, nh                                   !|
        Real*4     pmin, pmax,                              !|
     *             hmin, hmax,                              !V
     *             Eps                                      !|
        Parameter  (  np = 4001,                            !|
     *                nh = 91,                              !|
     *              pmin = 0.,                              !|
     *              pmax = 100.,   !<------------------------+
     *              hmin = 0.1,
     *              hmax = 1.0,
     *               Eps = 1.E-10)

        Integer    iscan, iascii,
     *             lun, istat,
     *             ifail, i, j
        Real*4     zmin, zmax,
     *             Array(np),      !Table of G_n(p) according to Sinha
     *             dp,             !(pmax-pmin)/(np-1)
     *             dh,             !(hmax-hmin)/(nh-1)
     *             p, h, Z4

        Real*8     CorrInt8
        External   CorrInt8

c ---------------
c       lun = 1
c ---------------
c Open correlation function table (CFT):
        Call OpenBinary ('corrfun.cft',lun,'new','write',istat,ifail)
        if (ifail.ne.0) goto 11

        write   (lun,err=13)    np,nh,
     *                          pmin,pmax,
     *                          hmin,hmax

        dp = (pmax-pmin)/(np-1)
        dh = (hmax-hmin)/(nh-1)

        do      i=1,nh  !=======================================+
          h = hmin + dh*(i-1)                                  !|
          do    j=1,np  !=+                                     |
            Array(j)=0.  !|                                     |
          enddo  !========+                                     |
          Z4 = 0.                                              !|
          zmin = 1.E+32                                        !|
          zmax =-1.E+32                                        !|
          do    j=1,np  !===============================+       |
            p = pmin+dp*(j-1)                          !|       |
            Array(j)=SNGL(CorrInt8(h,p,Z4,Eps,ifail))  !|       |
            if (ifail.ne.0) Then  !------------+        |       |
              write (*,*) ' ifail=',ifail,    !|        |       |
     *                    ' h=',h,' p=',p     !|        |       |
              goto 28                         !|        |       |
            endif  !---------------------------+        |       |
            if (Array(j).gt.zmax) zmax=Array(j)        !|       |
            if (Array(j).lt.zmin) zmin=Array(j)        !|       |
          enddo  !======================================+       |
          write (lun,err=13) Array                             !|
          write (*,1) h,zmax,zmin                              !|
  1       format(' Done h=',f4.2,' Max=',g10.3,' Min=',g10.3)  !|
          Call KeyIn (0,iscan,iascii)                          !|
          if (iascii.eq.27) stop 'Interrupted!'                !|
        enddo  !================================================+
        stop 'Done!'
  11    continue
        stop 'Error opening <corrfun.cft>'
  13    continue
        stop 'Error writing <corrfun.cft>'
  28    continue
        end

c =======================================================

        Real*8 Function CorrInt8 (h4, p4, Z4, Eps4, ifail)
c
c ATTENTION: for alternative (old or unsuccessful) ideas about
c            calculating this integral see UTILITY/CORRFUN0.FOR
c            and  UTILITY/COFFINT8.FOR.
c
c+-------------------------------------------------------------+
c|Caculating the integral:                                     |
c|                                                             |
c|     +inf                     +inf                           |
c|      |    -x^{2h}          1  |    -y  (1/2h-1)      1/2h   |
c|G(p)= |dx e      cos(px) = --- |dy e   y        cos(py    )  |
c|      |                    2*h |                             |
c|      0                        0                             |
c|                         +---------------------------------+ |
c|                           (this is currently not used!)     |
c|                                                             |
c|-- for p={0.,100.}                                           |
c|The details see in S.K.Sinha -- J.Physique III (France)      |
c|v.4, p.1543-1557 (1994)  --- Eq.(25) on p.1552.              |
c|-------------------------------------------------------------|
c|For h=0.5 (Gradstein, p.512, [Eq.3.893] -- English Edition): |
c|                                                             |
c|    +infinity                                                |
c|        |      -x             1                              |
c|G(p) =  |  dx e  cos(px) = --------                          |
c|        |                  1 + p*p                           |
c|        0                                                    |
c|-------------------------------------------------------------|
c|For h=1. (Gradstein, p.515, [Eq.3.896] -- English Edition):  |
c|                                                             |
c|    +infinity                                                |
c|        |     -x^2          sqrt(pi)   -(p^2)/4              |
c|G(p) =  | dx e    cos(px) = --------  e                      |
c|        |                      2                             |
c|        0                                                    |
c|-------------------------------------------------------------|
c|For p=0 (should be specified) to CorrInt8 since it is used   |
c|as an upper estimate for G(p):                               |
c|                                                             |
c|    +infinity            +infinity                           |
c|        |     -x^{2h}   1    |     {1/2h-1}  -y   1   ( 1  ) |
c|G(0) =  | dx e       = -- *  | dy y         e  = -- GF| -- | |
c|        |              2h    |                   2h   ( 2h ) |
c|        0                    0                               |
c|                                                             |
c| where GF is a well tabulated Gamma-function                 |
c+-------------------------------------------------------------+
c Possible numerical integrators:
c
c D01AMF 1-D quadrature, adaptive, infinite    |Used in old
c        or semi-infinite interval             |versions
c
c D01ASF 1-D quadrature, adaptive,             |Used in this
c        semi-infinite interval, weight        |version starting
c        function cos(wx) or sin(wx)           |by March, 1998.
c
c This one can be used if we substute y=x^{2h} and limit the
c integration range to, e.g. [0.0-10.0]:
c
c D01AJF 1-D quadrature, adaptive,             |This I did
c        finite interval, strategy due         |not check!
c        to Piessens & de Doncker,             |
c        allowing for badly-behaved            |
c        integrands                            |
c---------------------------------------------------------------

        Real*4     h4, p4, Z4, Eps4
        Real*8     h8, p8, v8
        Common     /CorrFFF/  h8

c       Logical         use_NAG /.True./
c       Logical         use_NAG /.False./

        Integer    ifail

        Real*4     sqpi2           ! sqrt(pi)/2

c       Real*8    S14AAF           ! Gamma-function (GF)
c       External  S14AAF
c--------------------------------------------------
c Definitions for both D01ASF (NAG) and DQAWF (QUADPACK) integral
c calculators (they use the same algorithm as D01ASF is based on the 
c extended version of DQAWF called DQAWFE, but QUADPACK is opensource):

        Integer    LstLimit
        Parameter (LstLimit = 400)              !suggested: LstLimit=50
        Integer    key_Cos,                     !1=cos, 2=sin
     *             LstActual
        Real*8     Bound_Lo8,                   !lower integration boundary
     *             Result8,                     !integration result
     *             EpsRel8,                     !requested relative error
     *             EpsAbs8,                     !requested absolute error
     *             ErrAbs8                      !actual error

        Real*8     Fun2_h
        External   Fun2_h
c--------------------------------------------------
c Definitions for NAG program D01ASF:
c       Integer    LWORK_N, LIWORK_N
c       Parameter (LWORK_N  = 8000,             !suggested: LWORK_N=2000
c    *             LIWORK_N = LWORK_N/2+2)      !suggested: at least LWORK_N/2
c       Integer    IWORK_N(LIWORK_N),
c    *             IerLst(LstLimit)
c       Real*8     WORK_N(LWORK_N),
c    *             ErLst(LstLimit),
c    *             RsLst(LstLimit)
c--------------------------------------------------
c Definitions for QUADPACK program DQAWFE
        Integer    MAXP1, LWORK_Q, LIWORK_Q
        Parameter (MAXP1 = 21,                  !no suggestion, taken from D01ASF
     *             LIWORK_Q = 4000,             !suggested: LIWORK_Q=1002
     *             LWORK_Q = 2*LIWORK_Q+25*MAXP1)  !as suggested
        Integer    IWORK_Q(LIWORK_Q),
     *             NEVAL
        Real*8     WORK_Q(LWORK_Q)
c--------------------------------------------------
        h8       = Dble(h4)
        p8       = Dble(p4)
        EpsRel8  = Eps4
        if (EpsRel8.lt.1.D-14) EpsRel8=1.D-14
        sqpi2    = Sqrt(atan(1.))       ! = Sqrt(pi/4) = Sqrt(pi)/2
        CorrInt8 = 0.0D0
        ifail    = 0

c The value of function at p=0:
        if (abs(Z4).lt.1.E-20) Then !-------------+
          v8 = 1./(2.*h8)                        !|
c         if (use_NAG) Then !------------------+  |
c           Z4 = SNGL(v8 * S14AAF(v8,ifail))  !|  | S14AAF=Gamma Function
c         else !-------------------------------+  |
            Z4 = SNGL(v8 * DGAMMA(v8))        !|  | GNU Fortran function  
c         endif  !-----------------------------+  |
        endif  !----------------------------------+

        if (abs(p8).lt.1.E-20) Then !-------------------------------------------+
                                                                               !|
          CorrInt8 = Z4                                                        !|
                                                                               !|
        else  !-----------------------------------------------------------------+
                                                                               !|
          if (abs(h8-0.5).lt.1.E-20) Then !----------------------------------+  |
                                                                            !|  |
c Gradstein, p.491:                                                         !|  |
            CorrInt8 = 1./(1.+p8*p8)                                        !|  |
                                                                            !|  |
          elseif (abs(h8-1.0).lt.1.E-20) Then !------------------------------+  |
                                                                            !|  |
c Gradstein, p.494:                                                         !|  |
            if (p8.lt.16.) CorrInt8 = Sqpi2*exp(-p8*p8/4.)      !min=exp(-64)|  |
                                                                            !|  |
          else   !-----------------------------------------------------------+  |
                                                                            !|  |
c Numerical evaluation of the integral:                                      |  |
c (NAG D01ASF is based on the QUADPACK routine DQAWFE                        |  |
c (Piessens et al. (1983))                                                   |  |
            Bound_Lo8 = 0.0D+0                                              !|  |
            key_Cos   = 1                   !1=cos, 2=sin in D01ASF          |  |
            EpsAbs8   = EpsRel8 * Abs(Z4)                                   !|  |
            ifail     =-1                   !0,-1,1. -1=request explanation s|  |
                                                                            !|  |
c           if (use_NAG) Then !----------------------------------------+     |  |
c             Call D01asf(Fun2_h,            !function to integrate    |input|  |
c    *                    Bound_Lo8,         !lower integration limit  |input|  |
c    *                    p8,                !in sin(p*x),cos(p*x)     |input|  |
c    *                    key_Cos,           !1=cos or 2=sin           |input|  |
c    *                    EpsAbs8,           !accuracy required        |input|  |
c    *                    Result8,           !integral value           |outpu|  |
c    *                    ErrAbs8,           !result error             |outpu|  |
c    *                    LstLimit,          !NAG suggested LstLimit=50|input|  |
c    *                    LstActual,         !actually used LstLimit   |outpu|  |
c    *                    ErLst,             !array of errs per intrval|outpu|  |
c    *                    RsLst,             !array of contributions   |outpu|  |
c    *                    IerLst,            !array of error flags     |outpu|  |
c    *                    WORK_N, LWORK_N,   !NAG suggested LWORK=2000 |work |  |
c    *                    IWORK_N, LIWORK_N, !NAG sugstd LIWORK=LWORK/2|work |  |
c    *                    ifail)             !exit flag                |I/O  |  |
c           else !-----------------------------------------------------+     |  |
              Call DQAWF (Fun2_h,            !function to integrate    |input|  |
     *                    Bound_Lo8,         !lower integration limit  |input|  |
     *                    p8,                !in sin(p*x),cos(p*x)     |input|  |
     *                    key_Cos,           !1=cos or 2=sin           |input|  |
     *                    EpsAbs8,           !accuracy required        |input|  |
     *                    Result8,           !integral value           |outpu|  |
     *                    ErrAbs8,           !result error             |outpu|  |
     *                    NEVAL,             !No.of integrand evaluat-s|outpu|  |
     *                    ifail,             !exit flag                |I/O  |  |
     *                    LstLimit,          !NAG suggested LstLimit=50|input|  |
     *                    LstActual,         !actually used LstLimit   |outpu|  |
     *                    LIWORK_Q,          !length of IWORK_Q        |work |  |
     *                    MAXP1,             !No.of Chebyshev moments  |work |  |
     *                    LWORK_Q,           !length of WORK_Q         |work |  |
     *                    IWORK_Q,           !int array of L=LIWORK_Q  |work |  |
     *                    WORK_Q)            !real array of L=LWORK_Q  |work |  |
c           endif !----------------------------------------------------+     |  |
c ifail=1 max divisions reached                                              |  |
c ifail=2 round-off error prevents achieving requested tolerance             |  |
c ifail=3 extremely bad local integrand behaviour                            |  |
c ifail=4 cannot achieve requested tolerance due to extrapolation            |  |
c ifail=5 integral is probably divergent, or slowly convergent               |  |
c ifail=6 key # 1 or 2, or LstLimit < 3                                      |  |
c ifail=7 bad integration behaviour within one or more intervals             |  |
c ifail=8 max number of intervals achierved                                  |  |
c ifail=9 extrapolation table does not converge to requested accuracy        |  |
c ifail=10 LWORK<5 or LIWORK<2                                               |  |
c ifail=-99 unexpected error                                                 |  |
c ifail=-999 dynamic RAM allocation faild                                    |  |
            if (ifail.ne.0) ifail = ifail + 900                             !|  |
            CorrInt8 = Result8                                              !|  |
          endif  !-----------------------------------------------------------+  |
                                                                               !|
        endif  !----------------------------------------------------------------+

        return
        end

c =======================================================

        Real*8  Function  Fun2_h(x8)
        Real*8     h8, x8, y8
        Common  /CorrFFF/  h8
c+--------------------------+
c|Calculating the integral: |
c|                          |
c|     +inf                 |
c|      |    -x^{2h}        |
c|G(p)= |dx e      cos(px)  | -- where cos is NOT specified!
c|      |                   |(it is a weight function in D01ASF)
c|      0                   |
c+--------------------------+
c y = x^{2h}
c log(y) = 2h*log(x)
c y = exp(log(y) ) = exp(2h*log(x)
c
        if (x8.gt.0.D+0)  Then   !-------+
          y8     = Dexp(2.*h8*Dlog(x8)) !| x^{2h}
c         if (y8.gt.708.) Then !------+  | max exponent (leads to 3.D-308)
          if (y8.gt.478.) Then !------+  | max exponent (leads to 3.D-208)
            Fun2_h = 0.0D0           !|  |
          else !----------------------+  |
            Fun2_h = Dexp(-y8)       !|  | exp(-x^{2h})
          endif !---------------------+  |
        else  !--------------------------+
          Fun2_h = 1.0D+0               !|
        endif  !-------------------------+
        return
        end

