        Program BrL_Surf
c+=====================================================+
c| This program calculates multiple diffraction of     |
c|                                                     |
c| X-rays in 3 to 12 Bragg and Laue cases including    |
c|                                                     |
c| SURFACE waves, i.e. those parallel to the surface   |
c| =============                                       |
c|                                                     |
c| The code is based on the paper:                     |
c|                                                     |
c| S.A.Stepanov and A.P.Ulyanenkov                     |
c| "A new algorithm for computation of x-ray multiple  |
c| Bragg diffraction"                                  |
c| Acta Crystallogr., v.A50, Part 5, p.579-585, (1994).|
c|                                                     |
c| Initial implementation: Stepanov & Ulyanenkov;      |
c| All versions since 1996: Stepanov                   |
c|                                                     |
c| Since 1997 the program is operational online at:    |
c|  https://x-server.gmca.aps.anl.gov/BRL.html         |
c|                                                     |
c|Subroutines: x0h2 + brl_12s + serve32.lib + nag2.lib |
c|                                                     |
c| A T T E N T I O N:                                  |
c| -----------------                                   |
c| With Microsoft Fortran Compiler, use /Gt key, as    |
c| otherwise data does not fit 64KB.                   |
c|                                                     |
c| With NDP Fortran Compiler, the program GJYF02 called|
c| from F02GJF should be compiled with the -cg1 key    |
c| (likely a compiler bug).                            |
c|---------------------------------------------------- |
c|2002/11: the Win32 version: the NAG library calls    |
c|         replaced by the NAG2 calls (everything went |
c|         to double precision).                       |
c|2013/04: added printing separate polarizations       |
c|         (DAT files only)                            |
c|                                         Sergey      |
c+=====================================================+
c The following include contain something like this:
c       Integer         nDimen
c       Parameter       (nDimen = 12)
c       Integer         nDimenTheta
c       Parameter       (nDimenTheta = 501)
c+===============================+
        Include 'brls_12.inc'   !|defines nDimen,nDimenTheta
c+===============================+
        Integer         nDimen4
        Parameter       (nDimen4 = 4*Ndimen)
        Integer         nWrongMax
        Parameter       (nWrongMax=3)

        Complex*16      Scatt_Matrix(nDimen4,nDimen4),    !<-+The
     *                  Gauss_Matrix(nDimen4,nDimen4),    !<-+Same!
     *                  Incident_Field(nDimen4,2),
     *                  Incident_Field_Saved(nDimen4,2),
     *                  Lambda(nDimen4,2), cf, qf,
     *                  Epsilon(nDimen4),
     *                  Fields(nDimen4,nDimen4),
     *                  Expon(nDimen4), Power,
     *                  Phi(nDimen), Surf_Coeff,
     *                  dlt
        Equivalence (Scatt_Matrix(1,1), Gauss_Matrix(1,1))

        Real*8          dGamma0, Gamma0, Gamma_local,
     *                  ReflexCoeff(nDimenTheta,nDimen), Rc,
     *                  ReflexCoeffPol(nDimenTheta,nDimen,2),           !for individual polarizations of diffracted waves
     *                  ReflexPhaseDiff(nDimenTheta,nDimen),            !phase difference between individual polarizations of diffracted waves
     *                  Left_Re(nDimen4,nDimen4),
     *                  Left_Im(nDimen4,nDimen4),
     *                  Right_Re(nDimen4,nDimen4),
     *                  Right_Im(nDimen4,nDimen4),
     *                  D_Re(nDimen4,nDimen4),
     *                  D_Im(nDimen4,nDimen4),
     *                  Scatt_Re(nDimen4,nDimen4),
     *                  Scatt_Im(nDimen4,nDimen4),
     *                  Gamma_Matr(nDimen4,nDimen4),
     *                  Epsilon_R(nDimen4),
     *                  Epsilon_I(nDimen4),
     *                  Epsilon_I_Min, Epsilon_I_Max,
     *                  Tolerance,
     *                  u, Phase

        Real            alf(nDimen), Tet(2), Renorm,
     *                  Seconds, Start, Finish,
     *                  Power_1, Power_2, c,
     *                  ReflexCoeff_minimum

        Logical         print_polarizations /.True./
c       Logical         print_polarizations /.False./

        Integer         Placement, Reduced_Size, Reduced_Size2,
     *                  Num_Positive, Num_Negative, i_Ref,
     *                  Num_Dropped_Bragg, i_Bragg_Begin,
     *                  Ext_Polar, Int_Polar,
     *                  N_Error, N_Correct, Times_BounSolv,
     *                  i, j, l, m, n, ifail, n2, n1, ni,
     *                  nScan, nAscii, Index_I_Min,
     *                  luninp, luntbl, narg, lines, iDebug,
     *                  nWrong, LineFile(nDimen), lgrd, lnul

        Integer         iarg

        Character       wrk*128, grdprg*24, nude*24

        Real*8          DVecSca2
        External        DVecSca2

        Integer  Nargs
        External Nargs

        Integer   System                !this function is Int*2 in DOS Fortran
c       External  System                !and Int*4 in Windows Fortran
c       Intrinsic System

        Character       getOS*8
        External        getOS
c-------------------------------------------------------
        Complex*16      xh(nDimen,nDimen), Imaginary1
        Real*8          k_module,
     *                  ca1(nDimen),ca2(nDimen),
     *                  ca3(nDimen),ca4(nDimen),
     *                  k_unit(3,nDimen+1),
     *                  Sg_unit(3,nDimen),
     *                  Pi_unit(3,nDimen),
     *                  csg, cpi, cSurf1, cSurf2,
     *                  Gamma_Br(nDimen+1),
     *                  Total_Reflex_Angle,
     *                  Dpi, Dgra
        Real            pi, xabs, Thickness,
     *                  xhq(nDimen,nDimen),
     *                  xfq(nDimen,nDimen),
     *                  Tet0(2), Tetk(2), dTet(2),
     *                  UnitData(5), UnitCoef(2)
        Logical         Surface_Case(nDimen)
        Integer         nReflexes, nReflexes2,
     *                  ipol, iout, ii, ipv,
     *                  nMaxPoints(2),
     *                  Surf_Ref_Num,
     *                  Laue_Ref_Num,
     *                  Bragg_Ref_Num,
     *                  Incident_Place,
     *                  iWarning, Task_Size,
     *                  New_ord(nDimen),
     *                  New_pos(nDimen),
     *                  linp, lout, iUnits(2),
     *                  iBatch_mode, ixway
        Character       Refl_Type(nDimen+1)*8,
     *                  ProgramVer*8,
     *                  InpFile*80, OutFile*80,
     *                  Comment(3)*80
        Common /BrSurf/ xh, Imaginary1,                 !complex*16
     *                  k_module,                       !real*8
     *                  ca1, ca2, ca3, ca4,             !real*8
     *                  k_unit, Sg_unit, Pi_unit,       !real*8
     *                  csg, cpi,                       !real*8
     *                  cSurf1, cSurf2,                 !real*8
     *                  Gamma_Br, Total_Reflex_Angle,   !real*8
     *                  Dpi, Dgra,                      !real*8
     *                  pi, xabs,                       !real
     *                  xhq, xfq,                       !real
     *                  Thickness,                      !real
     *                  Tet0, dTet, Tetk,               !real
     *                  UnitData, UnitCoef,             !real
     *                  Surface_Case,                   !Log
     *                  nReflexes, nReflexes2,          !int
     *                  ipol, iout,                     !int
     *                  ii, nMaxPoints,                 !int
     *                  Surf_Ref_Num,                   !int
     *                  Laue_Ref_Num,                   !int
     *                  Bragg_Ref_Num,                  !int
     *                  Incident_Place, iWarning,       !int
     *                  New_ord, New_pos,               !int
     *                  Task_Size, ipv,                 !int
     *                  linp, lout, ixway,              !int
     *                  iBatch_mode, iUnits,            !int
     *                  Refl_Type, ProgramVer,          !char*8
     *                  InpFile, OutFile, Comment       !char*80


c /New of August'99: interface to Henke and Cowan/Brennan
c databases/:
c       Integer             iHenkeCowan
c       Common /HenkeCowan/ iHenkeCowan

        Character       txt(20)*80
        Common  /msg/   txt
c The batch-mode processing (no dialogs and screen output)
c  -- added October-1997:
        Integer         modebat, lunbat, ierrbat, istackrezv
        Character       progname*16, progrezv(10)*16, ErrorFile*80
        Common  /batmode/ modebat, lunbat, ierrbat,
     *                    progname, progrezv, istackrezv, ErrorFile

c----------------------------------------------------------------------
        if (getOS() .eq. 'windows') Then !---+
           grdprg = 'grd2dat.exe'           !|
           nude   = ' > NUL 2>&1'           !|
        else !-------------------------------+
           grdprg = 'grd2dat'               !|
           nude   = ' > /dev/null 2>&1'     !|
        endif !------------------------------+
        lgrd = Len_Trim(grdprg)
        lnul = Len_Trim(nude)
        iDebug     = 0
        modebat    = 1         !assume batch-mode by default (for FOR32)
        istackrezv = 0
        progname   = 'BrL_Surf'
        ipv = Min(Len(ProgramVer),Len(progname))
        ProgramVer = progname(1:ipv)
        ipv = Len_Trim (ProgramVer)
        nWrong     = 0

        c = 0.
        ReflexCoeff_minimum = 0.0
c       ReflexCoeff_minimum = X02ABE(c)  !minimum (real*4): 1.1754944E-38
        c = 0.

c This is preliminary (until we've figured out the output filename):
        ErrorFile = ProgramVer(1:ipv)//'.err'
        Call    DeleteErrorFile ()
c-------------------------------------------------------
        narg    = Nargs()-1
        if (narg.eq.0)  Then  !---------------------+
                                                   !|
          InpFile = ProgramVer(1:ipv)//'.inp'      !|
                                                   !|
        else  !-------------------------------------+
                                                   !|
          iarg = 1                                 !|
          Call get_arg_wrapper (iarg,InpFile)      !|
                                                   !|
        endif  !------------------------------------+

        ifail  = 0
        luninp = 1
        luntbl = 3
c
c +--------------------------------------------------+
        Call Read_Input_File (luninp,luntbl,ifail)  !|
c +--------------------------------------------------+
c
        if (ifail.ne.0) goto 500  !------------stop--------------------+
                                                                      !v
        Start = Seconds(0.)

c A normalization factor for the Epsilon roots of scattering matrix:
        if (Surf_Ref_Num.eq.0)  Then !---------+
          Renorm = xabs                       !|no surface waves
        else  !--------------------------------+
          Renorm = SNGL(Total_Reflex_Angle)   !|surface waves present
        endif  !-------------------------------+

c+============ Preparation of Common Matrices ===========+
c|  First Step - Scanning Angle Independ Part            |
c|  -------------------------------------------          |
c|    +-Depend on Scanning Angles-+                      |
c|    |                           |                      |
c|    |    Left-Matrix :          |    Right-Matrix :    |
c| +  v                     + +   v                    + |
c| |+==========+|           | |+=========+|  Indentity | |
c| ||  Scatt.  |   0-matrix | ||2*(GamBr+| with respect| |
c| ||  Matrix  || (2k x 2n) | || +D(Gam0)||  to reflex | |
c| ||(2n x 2n) |            | ||(2n x 2n)|    numbers  | |
c| |+==========+|           | |+=========+|  (2k x 2n) | |
c| | - -- -- --- -- -- -- --| | - -- -- -- -- -- -- -- | |
c| |  0-matrix  |    Unit   | |  T-matrix |   0-matrix | |
c| |                matrix  | |   to R-T               | |
c| | (2n x 2k)  | (2k x 2k) | | (2n x 2k) |  (2k x 2k) | |
c| +                        + +                        + |
c| n - order of diffraction problem ( x 2 polarisations) |
c| k - number of surface reflexes   ( x 2 polarisations) |
c+=======================================================+

c                    +================+
c                    |   Left Matrix  |
c                    +================+

c Calculate the constant elements at the left side of the scattering
c matrix according to the book by Pinsker, page :
c Scatt_Matrix(n,m,s,ss) = [ alf(n)*delta(n,m,s,ss) -
c                xh(n,m)*(e(n,s)*e(m,ss) ]
c.......................................................
c Scattering matrix elements with odd (1,3,5,7,..) indices
c correspond SIGMA-polarization, while the elements with
c even (2,4,6,8,..) indices correspond to PI-polarization:

        do      i=1,nReflexes    !=========================+
          do    j=1,nReflexes      !======================+|
            m  = 2*(i-1)                                 !||
            n  = 2*(j-1)                                 !||
            cf = xh (i,j)                                !||
            qf = exp (Imaginary1*pi*xfq(i,j))            !||
            qf = Imaginary1 * xhq(i,j) * qf              !||
            u  = DVecSca2 (k_unit(1,i),k_unit(1,j),ii)   !||
            cf = cf + qf*u                               !||
                                                         !||
            Scatt_Matrix(m+1,n+1) =                      !||
     =       cf* DVecSca2(Sg_unit(1,i),Sg_unit(1,j),ii)+ !||
     +       qf*(DVecSca2(Sg_unit(1,i),k_unit(1,j),ii) * !||
     *           DVecSca2(k_unit (1,i),Sg_unit(1,j),ii)) !||
                                                         !||
            Scatt_Matrix(m+1,n+2) =                      !||
     =       cf* DVecSca2(Sg_unit(1,i),Pi_unit(1,j),ii)+ !||
     +       qf*(DVecSca2(Sg_unit(1,i),k_unit(1,j),ii) * !||
     *           DVecSca2(k_unit (1,i),Pi_unit(1,j),ii)) !||
                                                         !||
            Scatt_Matrix(m+2,n+1) =                      !||
     =       cf* DVecSca2(Pi_unit(1,i),Sg_unit(1,j),ii)+ !||
     +       qf*(DVecSca2(Pi_unit(1,i),k_unit(1,j),ii) * !||
     *           DVecSca2(k_unit (1,i),Sg_unit(1,j),ii)) !||
                                                         !||
            Scatt_Matrix(m+2,n+2) =                      !||
     =       cf* DVecSca2(Pi_unit(1,i),Pi_unit(1,j),ii)+ !||
     +       qf*(DVecSca2(Pi_unit(1,i),k_unit(1,j),ii) * !||
     *           DVecSca2(k_unit (1,i),Pi_unit(1,j),ii)) !||
                                                         !||
          enddo   !=======================================+|
        enddo !============================================+

        do      i=1,nReflexes2   !=======================+
          do    j=1,nReflexes2      !===================+|
            if (Cdabs(Scatt_Matrix(i,j)).lt.1.D-5)     !||
     *                       Scatt_Matrix(i,j)=(0.,0.) !||
            Scatt_Re(i,j) = Real (Scatt_Matrix(i,j))   !||
            Scatt_Im(i,j) = Imag (Scatt_Matrix(i,j))   !||
          enddo    !====================================+|
        enddo    !=======================================+

c Unit matrix at the left lower side with the respect to the diagonal:
          do i = nReflexes2+1, Task_Size  !===+
            Scatt_Re (i,i) = 1.              !|
            Scatt_Im (i,i) = 0.              !|
          enddo    !==========================+

c Zero rectangular matrices at the left along the anti-diagonal:
          do i = nReflexes2+1, Task_Size !====+
            do j = 1, nReflexes2 !====+       |
              Scatt_Re (i,j) = 0.    !|       |
              Scatt_Im (i,j) = 0.    !|       |
              Scatt_Re (j,i) = 0.    !|       |
              Scatt_Im (j,i) = 0.    !|       |
            enddo    !================+       |
          enddo    !==========================+

        write (luntbl,776,err=301)
  776   format  (';'/
     *  '; ========= Scattering matrix G:')
        write (luntbl,'(a)',err=301) '; ---------- Real(G):'
        do      i=1,Task_Size   !========================+
          write (luntbl,777,err=301) (Scatt_Re(i,j),    !|
     *                                j=1,Task_Size)    !|
        enddo   !========================================+
        write (luntbl,'(a)',err=301) '; ---------- Imag(G):'
        do      i=1,Task_Size   !========================+
          write (luntbl,777,err=301) (Scatt_Im(i,j),    !|
     *                                j=1,Task_Size)    !|
        enddo   !========================================+
  777   format  (';',48(1x,g12.4))

c                    +================+
c                    |  Right Matrix  |
c                    +================+

          Placement = 1
          do i_Ref = 1, nReflexes    !================+
            j = 2 * i_Ref                            !|
            Gamma_Matr(j,j)    = 2.* Gamma_Br(i_Ref) !|
     *                         * Renorm / xabs       !|
            Gamma_Matr(j-1,j-1)= Gamma_Matr(j,j)     !|
                                                     !|
            if (Surface_Case(i_Ref)) Then  !---+      |
              m = nReflexes2 + 2*Placement    !|      |
              Gamma_Matr (j  ,m  ) = 1        !|      |
              Gamma_Matr (j-1,m-1) = 1        !|      |
              Gamma_Matr (m  ,j  ) = 1        !|      |
              Gamma_Matr (m-1,j-1) = 1        !|      |
              Placement = Placement + 1       !|      |
            endif       !----------------------+      |
          enddo    !==================================+

c Zero matrix at the right along the diagonal:
          do    i = nReflexes2+1, Task_Size  !=====+
            do  j = nReflexes2+1, Task_Size  !===+ |
              Gamma_Matr (i,j) = 0.             !| |
            enddo  !=============================+ |
          enddo    !===============================+


c+======================================================+
c| The right-size part of the equations set for finding |
c|     the excitation coefficients Lambda of different  |
c|                 solution  branches                   |
c|        (Incident_Field_Saved E Incident_Field)       |
c|  --------------------------------------------------- |
        do      i=1,2       !=======================+   |
          do    j=1,nDimen4     !=================+ |   |
            Incident_Field_Saved(j,i) = (0.,0.)  !| |   |
          enddo     !=============================+ |   |
        enddo     !=================================+   |
c                                                       |
c   The Placement of Incident Wave Fields is            |
c   ----------------------------------------            |
        j = 2*(Incident_Place-1)+1                     !|
                                                       !|
        if (ipol .eq. 0)  Then  !--------------------+  |
c  +------------------------------------+            |  |
c  |Incident radiation is not polarized |            |  |
c  | (the boundary conditions equations |            |  |
c  | need to be solved 2 times)         |            |  |
c  +------------------------------------+            |  |
          Incident_Field_Saved(j  ,1) =             !|  |
     *                       Cmplx(1./sqrt(2.),0.)  !|  |
          Incident_Field_Saved(j+1,2) =             !|  |
     *                       Cmplx(1./sqrt(2.),0.)  !|  |
          Times_BounSolv = 2                        !|  |
        else  !--------------------------------------+  |
c  +------------------------------------+            |  |
c  |Incident radiation is polarized     |            |  |
c  | (the boundary conditions equations |            |  |
c  | need to be solved 1 time)          |            |  |
c  +------------------------------------+            |  |
          Incident_Field_Saved (j  ,1) =            !|  |
     *                              Dcmplx(csg,0.)  !|  |
          Incident_Field_Saved (j+1,1) =            !|  |
     *                              Dcmplx(cpi,0.)  !|  |
          Times_BounSolv = 1                        !|  |
        endif  !-------------------------------------+  |
c+======================================================+

c#######################################################
c Open output files of these types:
c '.grd' - for writing grid data for Surfer
c '.dat' - for writing lines data for Origin, Gnuplot, Grapher...

        lout = Len_Trim(OutFile)
        if (iout .eq. 0)  Then  !-------------+
          wrk = OutFile(1:lout)//'_00.grd'   !|
        else  !-------------------------------|
          wrk = OutFile(1:lout)//'_00.dat'   !|
        endif  !------------------------------+

        do  m=1,nReflexes  !======================================+
          write (wrk(lout+2:lout+3),'(i2)')  m-1                 !|
          if (wrk(lout+2:lout+2).eq.' ') wrk(lout+2:lout+2)='0'  !|
                                                                 !|
c Open output files:                                             !|
          Call OpenList (wrk(1:lout+7),m+6,ifail)                !|
                                                                 !|
          if (ifail.ne.0)       goto 500   !----------------------+-+
                                                                 !| v
          if (iout.eq.0)        Then  !-------------------------+ |
c Write header for GRD file:                                   !| |
            write (m+6,56,err=302) nMaxPoints(1),              !| |
     *                             nMaxPoints(2),              !| |
     *                             Tet0(1),Tetk(1),            !| |
     *                             Tet0(2),Tetk(2),            !| |
     *                                  0.,0.                  !| |
  56        format ('DSAA'/2i4,3(/2g11.3))                     !| |
          endif     ! <-----------------------------------------+ |
        enddo     !===============================================+

c
c
c
c           Main loops over Theta-1 and Theta-2:
c
c

c-------------------------------------------------------
        do      n2=1,nMaxPoints(2) !===================================+
          Tet(2) = Tet0(2)+dTet(2)*(n2-1)                             !|
                                                                      !|
          do n1=1,nMaxPoints(1) !=====================================+|
            Tet(1) = Tet0(1)+dTet(1)*(n1-1)                          !||
c Current matrix index:                                              !||
                            ni=n1                                    !||
            if (iout.eq.2)  ni=n2                                    !||
                                                                     !||
            do i_Ref = 1, nReflexes   !=============+                 ||
              ReflexCoeff   (ni, i_Ref)    = 0.    !|                 ||
              ReflexCoeffPol(ni, i_Ref, 1) = 0.    !|                 ||
              ReflexCoeffPol(ni, i_Ref, 2) = 0.    !|                 ||
              ReflexPhaseDiff(ni, i_Ref)   = 0.    !|                 ||
            enddo !=================================+                 ||
c                                                                     ||
c+-----------------------------------------+                          ||
c| Second Step - Restoring of All Matrices |                          ||
c+-----------------------------------------+                          ||
c                                                                     ||
            do i = 1, Task_Size   !==================+                ||
              do j = 1, Task_Size   !=============+  |                ||
                Left_Re (i,j) = Scatt_Re(i,j)    !|  |                ||
                Left_Im (i,j) = Scatt_Im(i,j)    !|  |                ||
                Right_Re(i,j) = Gamma_Matr(i,j)  !|  |                ||
                Right_Im(i,j) = 0.               !|  |                ||
              enddo  !============================+  |                ||
            enddo  !=================================+                ||
                                                                     !||
c Calculate change in Gamma0:                                         ||
            Dgamma0 = cSurf1*Tet(1) + cSurf2*Tet(2)                  !||
            Gamma0  = Gamma_Br(Incident_Place) + Dgamma0             !||
                                                                     !||
            if (Gamma0.lt.0.)  Then !-----------------------------+   ||
              if (iWarning.ge.1) Then !-------------------------+ |   ||
c               write (   *,  444,err=300) Gamma0/UnitData(2), !| |   ||
c    *                                     n1,n2,              !| |   ||
c    *                                     Tet(1),Tet(2)       !| |   ||
                write (luntbl,444,err=301) Gamma0/UnitData(2), !| |   ||
     *                                     n1,n2,              !| |   ||
     *                                     Tet(1),Tet(2)       !v v   vv
  444           format(' ! Negative Gamma0=',g9.3,
     *          'min  -- at point=(',i3,',',i3,
     *          ')  Q=(',g9.3,',',g9.3,')')                    !^ ^   ^^
              endif  !------------------------------------------+ |   ||
              goto 555  !-------------------------------------+   |   ||
            endif  !------------------------------------------+---+   ||
c                                                             |       ||
c Exclude very small incident angles ( < 0.0001"):            v       ||
            if (Gamma0 .lt. 1.D-8) Gamma0 = 1.D-8                    !||
c                                                                     ||
c Update the matrices with the account for changed Alfa and Gamma0    ||
c                                                                     ||
            do  i_Ref=1,nReflexes  !================+                 ||
c Calculate alf(1:nReflexes):                       |                 ||
              alf(i_Ref) = SNGL(                   !|                 ||
     (                       ca1(i_Ref)*Tet(1)     !|                 ||
     +                     + ca2(i_Ref)*Tet(2)     !|                 ||
     +                     + ca4(i_Ref)            !|                 ||
     -                     + ca3(i_Ref)            !|                 ||
     -                     * (Tet(1)**2+Tet(2)**2) !|                 ||
     )                     )                       !|                 ||
              m = 2*(i_Ref-1)+1                    !|                 ||
c--------------                                    !|                 ||
              Left_Re (m,  m  ) = Left_Re (m,m)    !|                 ||
     +                          - alf(i_Ref)       !|                 ||
              Left_Re (m+1,m+1) = Left_Re (m,m)    !|                 ||
c--------------                                    !|                 ||
              Right_Re(m,  m  ) = Right_Re(m,m)    !|                 ||
     +                          + 2.* Dgamma0      !|                 ||
     *                          * Renorm / xabs    !|                 ||
              Right_Re(m+1,m+1) = Right_Re(m,m)    !|                 ||
            enddo  !================================+                 ||
c                                                                     ||
c+------------------------------------------------+                   ||
c| Solve the generalized eigenvalue problem:      |                   ||
c|                                                |                   ||
c|      Left * D  = Epsilon * Right * D           |                   ||
c|..and check the number of eigenvalues computed: |                   ||
c|                                                |                   ||
c|                |         |                     |                   ||
c|                | D(1:2n) |                     |                   ||
c|   D(2n + 2k) = |         | - Eigenvectors      |                   ||
c|                | D(add)  |   structure         |                   ||
c|                |         |                     |                   ||
c|                                                |                   ||
c|   n = nReflexes,    k = Surf_Ref_Num           |                   ||
c+------------------------------------------------+                   ||
                                                                     !||
            Ifail = 0                                                !||
c This is only used by NAG F02GJF and not used by LAPACK ZGGEV        ||
c Element of matrix is considered negligible if it is                 ||
c less, than Eps times the norm of its matrix.                        ||
c If Eps<=0., the value supplied by NAG X02AAE is used.               ||
c           Tolerance = (1.e-4)*xabs/Renorm     !this is inaccurate   ||
            Tolerance = 0.                                           !||
                                                                     !||
c+--------------------------------------------------------+           ||
c| Solve eigenvalue problem (NAG F02GJF or LAPACK ZGGEV): |           ||
c+--------------------------------------------------------+           ||
            Call Generalized_Eigenvalues (                           !||
     *                   Task_Size,                                  !||
     *                   Left_Re,   Left_Im,    !left matrix (destryd)||
     *                   Right_Re,  Right_Im,   !rght matrix (destryd)||
     *                   Epsilon_R, Epsilon_I,  !eigenvalues (output) ||
     *                   D_Re,      D_Im,       !eigenvectors (output)||
     *                   Tolerance,             !neglection level     ||
     *                   ifail,                 !out: failure flag    ||
     *                   wrk)                   !out: error message   ||
                                                                     !||
            if (ifail.eq.-999)  goto 507  !------------------------+  || insufficient wrk array
            if (ifail.eq.-1000) goto 505  !---------------------+  v  || zero eigenvalues divider
            if (ifail.ne. 0)    goto 503  !------------------+  v     || eigenvalue error
c                                                            v        ||
c+-------------------------------------------------------+            ||
c|Sort the solutions in decreasing order of Im(Epsilon): |            ||
c+-------------------------------------------------------+            ||
            Call Sort_by_Epsilon_Im_Part (nDimen4,                   !||
     *                   Task_Size,                                  !||
     *                   Epsilon_R, Epsilon_I,                       !||
     *                   D_Re,      D_Im,                            !||
     *                   Left_Re,   Left_Im)  !<=used as work arrays  ||
                                                                     !||
            if (iDebug. gt. 0) then !------------------------------+  ||
              write (luntbl,33) 1,n1,1,Tet(1),2,n2,2,Tet(2)       !|  ||
              do  i=1,Task_Size  !=============================+   |  ||
                write (luntbl,778,err=301) i, Epsilon_R(i),   !|   |  ||
     *                                        Epsilon_I(i)    !|   |  ||
                do j=1,Task_Size !=========================+   |   |  ||
                  write (luntbl,779,err=301) D_re(i,j),   !|   |   |  ||
     *                                       D_Im(i,j)    !|  !|   |  ||
                enddo !====================================+   |   |  ||
              enddo !==========================================+   |  ||
  778         format(i4,'. lambda=(',g14.4,',',g12.4,'), Vector:')!|  ||
  779         format(' (',g14.4,',',g12.4,')')                    !|  ||
            endif !------------------------------------------------+  ||
                                                                     !||
c+--------------------------------------------+                       ||
c|  Test on Correlation between the Sign of   |                       ||
c|  Im(Epsilon) and the Bragg Reflection Type |                       ||
c+--------------------------------------------+                       ||
            Num_Positive = 0                                         !||
            Num_Negative = 0                                         !||
            do i = 1, Task_Size  !===================+                ||
              if (Epsilon_I (i) .lt. 0.) Then !--+   |                ||
                 Num_Negative = Num_Negative+1  !|   |                ||
              else  !----------------------------+   |                ||
                 Num_Positive = Num_Positive+1  !|   |                ||
              endif   !--------------------------+   |                ||
            enddo  !=================================+                ||
                                                                     !||
            if (iWarning.ge.1)  Then  !----------------------------+  ||
              if ((2*(Surf_Ref_Num + Laue_Ref_Num)                !|  ||
     *                           .ne.                             !|  ||
     *                      Num_Positive )                        !|  ||
     *                           .OR.                             !|  ||
     *            (2*(Surf_Ref_Num + Bragg_Ref_Num)               !|  ||
     *                           .ne.                             !|  ||
     *                      Num_Negative ) )                      !|  ||
     *                                      Then  !-------------+  |  ||
                N_Error = Num_Positive - 2*                    !|  |  ||
     *                    (Surf_Ref_Num + Laue_Ref_Num)        !|  |  ||
c               write (  *,   341,err=300) N_Error, n1, n2     !|  |  ||
                write (luntbl,341,err=301) N_Error, n1, n2     !|  |  ||
                if (iWarning.ge.2)                             !|  |  ||
     *          write (luntbl,351,err=301)                     !|  |  ||
     *                         Laue_Ref_Num,                   !|  |  ||
     *                         Bragg_Ref_Num,                  !|  |  ||
     *                         Surf_Ref_Num,                   !|  |  ||
     *                         Num_Positive,                   !|  |  ||
     *                         2*(Surf_Ref_Num+Laue_Ref_Num),  !|  |  ||
     *                         Num_Negative,                   !|  |  ||
     *                         2*(Surf_Ref_Num+Bragg_Ref_Num), !|  |  ||
     *                         (Epsilon_I(i),i=1,Task_Size)    !v  v  vv
  341           format(
     *          ' ! Error [',i3,'] in number of Bragg/Laue',
     *          ' roots -- at point (',i3,',',i3,')')
  351           format(
     *          '      Laue-case reflections=',i2/
     *          '     Bragg-case reflections=',i2/
     *          '        Surface reflections=',i2/
     *          '   Positive (Laue)    roots=',i2,
     *          '   (Expected=',i2,')'/
     *          '   Negative (Bragg)   roots=',i2,
     *          '   (Expected=',i2,')'://
     *          ' Roots dump:'//16(4g12.3/))                   !^  ^  ^^
                                                               !|  |  ||
c Correction of possible precision-loss errors:                 |  |  ||
                if (iWarning.ge.3)  Then  !------------------+  |  |  ||
                                                            !|  |  |  ||
                  Epsilon_I_Max = Max (                     !|  |  |  ||
     *                        Abs(Epsilon_I(1)),            !|  |  |  ||
     *                        Abs(Epsilon_I(Task_Size)))    !|  |  |  ||
                                                            !|  |  |  ||
                  N_Correct = 0                             !|  |  |  ||
                                                            !|  |  |  ||
                  if (N_Error.gt.0) Then  !---------------+  |  |  |  ||
c If too many positive roots:                             |  |  |  |  ||
                    do  j=1,N_Error  !=================+  |  |  |  |  ||
                      Epsilon_I_Min = Epsilon_I_Max   !|  |  |  |  |  ||
                      Index_I_Min = 0                 !|  |  |  |  |  ||
c Find Epsilon with minimum positive Im:               |  |  |  |  |  ||
                      do  i=1,Task_Size   !===========+|  |  |  |  |  ||
                        if (Epsilon_I(i).gt.0.       !||  |  |  |  |  ||
     *                          .and.                !||  |  |  |  |  ||
     *                      Epsilon_I(i).lt.         !||  |  |  |  |  ||
     *                          Epsilon_I_Min)       !||  |  |  |  |  ||
     *                                      Then  !-+ ||  |  |  |  |  ||
                          Index_I_Min = i          !| ||  |  |  |  |  ||
                          Epsilon_I_Min =          !| ||  |  |  |  |  ||
     *                                Epsilon_I(i) !| ||  |  |  |  |  ||
                        endif  !--------------------+ ||  |  |  |  |  ||
                      enddo  !========================+|  |  |  |  |  ||
                      if (Index_I_Min.gt.0            !|  |  |  |  |  ||
     *                        .and.                   !|  |  |  |  |  ||
     *                    Epsilon_I_Min.lt.           !|  |  |  |  |  ||
     *                        Epsilon_I_Max/1.e+6)    !|  |  |  |  |  ||
     *                                      Then  !-+  |  |  |  |  |  ||
                        Epsilon_I(Index_I_Min)=    !|  |  |  |  |  |  ||
     -                     -Epsilon_I(Index_I_Min) !|  |  |  |  |  |  ||
                        N_Correct=N_Correct+1      !|  |  |  |  |  |  ||
                      endif  !----------------------+  |  |  |  |  |  ||
                    enddo  !===========================+  |  |  |  |  ||
                  else  !---------------------------------|  |  |  |  ||
                    N_Error = -N_Error                   !|  |  |  |  ||
                    do  j=1,N_Error  !=================+  |  |  |  |  ||
                      Epsilon_I_Min = Epsilon_I_Max   !|  |  |  |  |  ||
                      Index_I_Min = 0                 !|  |  |  |  |  ||
c Find Epsilon with minimum negative Im:               |  |  |  |  |  ||
                      do  i=1,Task_Size   !===========+|  |  |  |  |  ||
                        if (Epsilon_I(i).lt.0.       !||  |  |  |  |  ||
     *                          .and.                !||  |  |  |  |  ||
     *                      Epsilon_I(i).gt.         !||  |  |  |  |  ||
     *                         -Epsilon_I_Min)       !||  |  |  |  |  ||
     *                                      Then  !-+ ||  |  |  |  |  ||
                          Index_I_Min = i          !| ||  |  |  |  |  ||
                          Epsilon_I_Min =          !| ||  |  |  |  |  ||
     *                                Epsilon_I(i) !| ||  |  |  |  |  ||
                        endif  !--------------------+ ||  |  |  |  |  ||
                      enddo  !========================+|  |  |  |  |  ||
                      if (Index_I_Min.gt.0            !|  |  |  |  |  ||
     *                        .and.                   !|  |  |  |  |  ||
     *                    Epsilon_I_Min.gt.           !|  |  |  |  |  ||
     *                       -Epsilon_I_Max/1.e+6)    !|  |  |  |  |  ||
     *                                      Then  !-+  |  |  |  |  |  ||
                        Epsilon_I(Index_I_Min)=    !|  |  |  |  |  |  ||
     -                     -Epsilon_I(Index_I_Min) !|  |  |  |  |  |  ||
                        N_Correct=N_Correct+1      !|  |  |  |  |  |  ||
                      endif  !----------------------+  |  |  |  |  |  ||
                    enddo  !===========================+  |  |  |  |  ||
                  endif  !--------------------------------+  |  |  |  ||
                  if (N_Correct.eq.N_Error) Then  !-------+  |  |  |  ||
                    write (luntbl,*,err=301)             !|  |  |  |  ||
     *                              '! Corrected...'     !|  |  |  |  ||
                  else   !--------------------------------+  |  |  |  ||
                    write (luntbl,*,err=301)             !|  |  |  |  ||
     *                              '! Not corrected...' !|  |  |  |  ||
                    goto 500  !-------------exit----------+--+--+--+--++----+
                  endif    !------------------------------+  |  |  |  ||    v
                endif  !-------------------------------------+  |  |  ||
              endif    !----------------------------------------+  |  ||
            endif    !---------------------------------------------+  ||
c                                                                     ||
c +----------------------------------------+                          ||
c | Merge Epsilon parts into complex form  |                          ||
c | and Re-Initialize Incident Wave Fields |                          ||
c +----------------------------------------+                          ||
c                                                                     ||
            do  i=1,Task_Size   !=====================+               ||
              Epsilon(i) = Dcmplx( Epsilon_R(i),     !|               ||
     *                            Epsilon_I(i) )     !|               ||
              do    j=1,Task_Size     !=============+ |               ||
                Fields(j,i) = Dcmplx( D_Re(j,i),   !| |               ||
     *                               D_Im(j,i) )   !| |               ||
              enddo   !=============================+ |               ||
            enddo   !=================================+               ||
c                                                                     ||
c+====================================================+               ||
c| Structure of Wavefields Matrix:                    |               ||
c| Full matrix has dimensions: Task_Size * Task_Size, |               ||
c| where Task_size = 2*(nReflexes + Surf_Ref_Num)     |               ||
c| -------------------------------------------------  |               ||
c| Matrix structure:                                  |               ||
c|               +-------< Components D(h) (Reflexes) |               ||
c|               |                                    |               ||
c|               |                                    |               ||
c|               v Components D(j) (Epsilon)          |               ||
c|                                                    |               ||
c| The matrix dimension can be reduced as follows:    |               ||
c|                                                    |               ||
c|                   +=  Reduced_Size2 = 2*nReflexes  |               ||
c|               +=======+                            |               ||
c|              +-         -+                         |               ||
c|           += |          X|                         |               ||
c|           |  |           |   +=================+   |               ||
c|           |  |          X|   |lines, that can  |   |               ||
c|           |  |           |   |be dropped for   |   |               ||
c|           |  |          X|   |some Bragg-case  |   |               ||
c|           |  |           |   |reflexes in thick|   |               ||
c|     +=====|  |          X|   |crystal:   2 *   |   |               ||
c|     |     |  |           |   |Num_Dropped_Bragg|   |               ||
c|     |     += |          X| <===================+   |               ||
c|     |        |           |   +===================+ |               ||
c|     |        |X X X X X X| <=|lines, dropped for | |               ||
c|     v        +-         -+   | (2*Surf_Ref_Num)  | |               ||
c|Reduced_Size2=       ^   ^    |solutions with huge| |               ||
c|=2*nReflexes         |   |    |Im(Epsilon) < 0    | |               ||
c|                     |   |    +===================+ |               ||
c|+====================|   |=======================+  |               ||
c||columns, that can be|   | columns, dropped for  |  |               ||
c||dropped for some Br.|   |  (2 x Surf_Ref_Num)   |  |               ||
c||reflexes & thick cr.|   |nonphysical amplitudes |  |               ||
c|+====================+   +=======================+  |               ||
c+====================================================+               ||
c                                                                    !||
c+================================================+                   ||
c| Find the excitation coefficients of Lambda     |                   ||
c| branches according to the Pinsker book p.372   |                   ||
c| by using the Gaussian solution for a ser of    |                   ||
c| linear equations (the boundary conditions):    |                   ||
c|                                                |                   ||
c|<Column_Dim>                                    |                   ||
c|   Sum   Gauss_Matrix(m,s,j)*Lambda(j) =        |                   ||
c|   j=1                    = Incident_Field(m,s) |                   ||
c|                                                |                   ||
c|     -- where:                                  |                   ||
c|                                                |                   ||
c| Gauss_Matrix(m,s,j) = Fields(m,s,j)            |                   ||
c|     -- for Laue-case diffracted waves,         |                   ||
c|                                                |                   ||
c| Gauss_Matrix(m,s,j) = Fields(m,s,j)            |                   ||
c|                     * Exp(-i*k*l*Epsilon)      |                   ||
c|     -- for Bragg-case diffracted waves,        |                   ||
c|                                                |                   ||
c| Gauss_Matrix(m,s,j) = Fields(m,s,j) *          |                   ||
c|                     ( Delta_m0 +               |                   ||
c|                       Epsilon_h(j)/2.*Gamma0 ) |                   ||
c|     -- for surface waves.                      |                   ||
c|                                                |                   ||
c|   Here:                                        |                   ||
c| Epsilon_h(j) = Epsilon(j) + Dgamma0   +        |                   ||
c|               ( Gamma_Br(m) + Phi(m) ) ,       |                   ||
c| Delta_m0 = 1 for m=Incident_Place ,            |                   ||
c| Delta_m0 = 0 for m#Incident_Place ,            |                   ||
c|                                                |                   ||
c| Gamma must be expresssed in the units of       |                   ||
c|   Renorm = Total_Reflex_Angle                  |                   ||
c+================================================+                   ||
c  +--------------------------------+                                 ||
c  | Incident Field Matrix Restoring|                                 ||
c  +--------------------------------+                                 ||
            do      j=1,Times_BounSolv   !=========================+  ||
              do    i=1,Task_Size     !==========================+ |  ||
                Incident_Field(i,j) = Incident_Field_Saved(i,j) !| |  ||
              enddo   !==========================================+ |  ||
            enddo   !==============================================+  ||
c                                                                     ||
c  +-----------------------+                                          ||
c  | Gauss Matrix Preparing|                                          ||
c  +-----------------------+                                          ||
            do      i=1,Task_Size   !==============+                  ||
              do    j=1,Task_Size     !==========+ |                  ||
                Gauss_Matrix(i,j) = Fields(i,j) !| |                  ||
              enddo   !==========================+ |                  ||
            enddo   !==============================+                  ||
c                                                                     ||
c    Then : 1. Handling Bragg Reflexes                                ||
c           2. Handling Surface Reflexes                              ||
c                                                                     ||
c                                                                     ||
c     For Bragg reflexes (if there are) we test                       ||
c     the exponent value for the successive exc-                      ||
c     luding of nonphysical boundary conditions                       ||
c                                                                     ||
c                                                                     ||
c+----------------------------------+                                 ||
c| Field Amplitude Matrix Reducing  |                                 ||
c+----------------------------------+                                 ||
            Reduced_Size  = nReflexes                                !||
            Reduced_Size2 = nReflexes2                               !||
            Num_Dropped_Bragg = 0                                    !||
c                                                                     ||
            if ( Bragg_Ref_Num .ne. 0 )   Then  !---------------+     ||
c +-----------------------------+                               |     ||
c |Test for number of roots with|                               |     ||
c |    Im( Epsilon ) < -10^5    |                               |     ||
c +-----------------------------+                               |     ||
c                                                               |     ||
              do      i = Reduced_Size, 1, -1  !==============+ |     ||
                j = 2 * i                                    !| |     ||
                Power_1 = -(1.e+04) * SNGL( k_module         !| |     ||
     *                                    * Thickness        !| |     ||
     *                                    * Renorm           !| |     ||
     *                                    * Epsilon_I(j) )   !| |     ||
                                                             !| |     ||
                Power_2 = -(1.e+04) * SNGL( k_module         !| |     ||
     *                                    * Thickness        !| |     ||
     *                                    * Renorm           !| |     ||
     *                                    * Epsilon_I(j-1) ) !| |     ||
                                                             !| |     ||
                if ((Power_1 .lt. 15.)                       !| |     ||
     *                      .AND.                            !| |     ||
     *              (Power_2 .lt. 15.)) goto 11 !--->----+    | |     ||
                                                        !|    | |     ||
                  Num_Dropped_Bragg =                   !|    | |     ||
     *                           Num_Dropped_Bragg + 1  !|    | |     ||
              enddo   !==================================+====+ |     ||
c                                                        |      |     ||
c +-----------------------------+                        |      |     ||
c | Exclude Bragg Reflexes for  |                        |      |     ||
c | which Crystal is 'Infinite' |                        |      |     ||
c +-----------------------------+                        |      |     ||
  11          continue !<--------------------------------+      |     ||
              Reduced_Size2 = nReflexes2                       !|     ||
     -                      - 2*Num_Dropped_Bragg              !|     ||
              Reduced_Size  = nReflexes                        !|     ||
     -                      -   Num_Dropped_Bragg              !|     ||
c                                                              !|     ||
            endif !---------------------------------------------+     ||
c                                                                     ||
c +-----------------------------+                                     ||
c | Calculate The Exponents     |                                     ||
c |    and Multiplying          |                                     ||
c |  Bragg Reflexes Remained    |                                     ||
c +-----------------------------+                                     ||
            do    j = 1, Reduced_Size2    !=========+                 ||
              Power = + (1.e+04) * Imaginary1 *    !|                 ||
     *                   k_module  *               !|                 ||
     *                   Thickness *               !|                 ||
     *                   Renorm *                  !|                 ||
     *                   Epsilon(j)                !|                 ||
              Expon(j) = Exp(Power)                !|                 ||
            enddo    !==============================+                 ||
c                                                                     ||
c+---------------------------------------------+                      ||
c|If there are any non-dropped Bragg reflexes, |                      ||
c|for the the matrix of boundary conditions    |                      ||
c|must be multiplied by an exponent (since the |                      ||
c|boundary conditions are at the exit, i.e.    |                      ||
c|(lower surface of the crystal)               |                      ||
c+---------------------------------------------+                      ||
            i_Bragg_Begin = 2*(Laue_Ref_Num+Surf_Ref_Num) + 1        !||
c                                                                     ||
            do i = i_Bragg_Begin, Reduced_Size2  !============+       ||
c+-----------------------------+                              |       ||
c|Shift Bragg-case reflexes, if|                              |       ||
c|the boundary conditions for  |                              |       ||
c|some of them are dropped:    |                              |       ||
c+-----------------------------+                              |       ||
              i_Ref = i + Num_Dropped_Bragg                  !|       ||
              do    j=1,Reduced_Size2  !==================+   |       ||
                Gauss_Matrix(i,j) = Fields(i_Ref,j)      !|   |       ||
     *                            * Expon(j)             !|   |       ||
              enddo  !====================================+   |       ||
            enddo  !==========================================+       ||
c                                                                     ||
c +-------------------------------------+                             ||
c | For Surface Reflexes (if There Are) |                             ||
c +-------------------------------------+                             ||
            do i_Ref = 1, Reduced_Size !==========================+   ||
                                                                 !|   ||
              if (Surface_Case(i_Ref))  Then  !---------------+   |   ||
                c = SNGL( (Gamma_Br(i_Ref)+Dgamma0)**2       !|   |   ||
     -                                  - Xabs*Alf(i_Ref) )  !|   |   ||
                if (c.gt.0.)   Then   !--------+              |   |   ||
                  Phi(i_Ref) = Sqrt(c)        !|              |   |   ||
                else   !-----------------------+              |   |   ||
                  Phi(i_Ref) = Cmplx(0.,      !|              |   |   ||
     *                               Sqrt(-c))!|              |   |   ||
                endif  !-----------------------+              |   |   ||
                                                             !|   |   ||
                l = 2 * i_Ref                                !|   |   ||
                do j = 1, Reduced_Size2  !================+   |   |   ||
                  if (i_Ref.ne.Incident_Place)           !|   |   |   ||
     *                                     Then !-----+   |   |   |   ||
                                                     !|   |   |   |   ||
                   Surf_Coeff = Epsilon(j)*Renorm    !|   |   |   |   ||
     +                        + Dgamma0              !|   |   |   |   ||
     +                        + Gamma_Br(i_Ref)      !|   |   |   |   ||
     +                        + Phi(i_Ref)           !|   |   |   |   ||
                  else  !-----------------------------+   |   |   |   ||
                    Surf_Coeff = 1.                  !|   |   |   |   ||
     *                         + Epsilon(j)*Renorm   !|   |   |   |   ||
     /                         / (2.*Gamma0)         !|   |   |   |   ||
                  endif !-----------------------------+   |   |   |   ||
                                                         !|   |   |   ||
                  Gauss_Matrix(l-1,j)=Fields(l-1,j)      !|   |   |   ||
     *                               *Surf_Coeff         !|   |   |   ||
                  Gauss_Matrix(l,  j)=Fields(l,  j)      !|   |   |   ||
     *                               *Surf_Coeff         !|   |   |   ||
                enddo   !=================================+   |   |   ||
              endif   !---------------------------------------+   |   ||
            enddo   !=============================================+   ||
                                                                     !||
c+---------------------+                                              ||
c| Gauss System Solving|                                              ||
c+---------------------+                                             !||
            ifail = 0                                                !||
            if (iDebug. gt. 0) then !------------------------------+  ||
              write (luntbl,33) 1,n1,1,Tet(1),2,n2,2,Tet(2)       !|  ||
  33          format(2(' n',i1,'=',i4,'  Theta(',i1,')=',g12.5/)) !|  ||
              write (luntbl,'(a)',err=301)                        !|  ||
     *                                 '; ---Real(Gauss_Matrix):' !|  ||
              do  i=1,Task_Size  !=============================+   |  ||
                write (luntbl,777,err=301)                    !|   |  ||
     *                (Real(Gauss_Matrix(i,j)),j=1,Task_Size) !|   |  ||
              enddo !==========================================+   |  ||
              write (luntbl,'(a)',err=301)                        !|  ||
     *                                 '; ---Imag(Gauss_Matrix):' !|  ||
              do  i=1,Task_Size  !=============================+   |  ||
                write (luntbl,777,err=301)                    !|   |  ||
     *                (Imag(Gauss_Matrix(i,j)),j=1,Task_Size) !|   |  ||
              enddo !==========================================+   |  ||
            endif !------------------------------------------------+  ||
                                                                     !||
            Call Gauss_Eigenvalues  (                                !||
     *                    Gauss_Matrix,   !left mtrx (input/destroyed)||
     *                    Incident_Field, !rght mtrx (input/destroyed)||
     *                    Lambda,         !solutions matrix(output)   ||
     *                    Reduced_Size2,  !number of equations(input) ||
     *                    Times_BounSolv, !nmbr of right sides (1,2)  ||
     *                    ifail)          !exit status (output)       ||
                                                                     !||
            if (ifail.ne.0)  goto 504  !-----------------+            ||
c                                                        v            ||
c+================================================+                   ||
c|  Find the reflection coefficients:             |                   ||
c|  ----------------------------------            |                   ||
c|       |Re{g(m)}|  2 |<Red_Dim>                 |                   ||
c| rh(m)=---------- sum|  Sum  D(m,s,j)*Lambda(j)*|                   ||
c|          g(0)    s=1| j = 1                    |                   ||
c|                                            |<  |                   ||
c|               *exp(m,j) - delta_m0*delta_z0|   |                   ||
c|                                            |   |                   ||
c| Where:                                         |                   ||
c| exp(m,j) = 1, for reflected beams,             |                   ||
c| exp(m,j) = Expon(j) for transmitted beams.     |                   ||
c| delta_m0 = 1, for m=Incident_Place,            |                   ||
c| delta_m0 = 0, for any other m.                 |                   ||
c| delta_z0 = 1, for reflected beams,             |                   ||
c| delta_z0 = 0, for transmitted beams.           |                   ||
c+================================================+                   ||
c                                                                     ||
c+---------------------------------+                                  ||
c| Loop over the Bragg reflections:|                                  ||
c+---------------------------------+                                  ||
            do    i_Ref  = 1, nReflexes   !========================+  ||
              ReflexCoeff   (ni,i_Ref)   = 0.                     !|  ||
              ReflexCoeffPol(ni,i_Ref,1) = 0.                     !|  ||
              ReflexCoeffPol(ni,i_Ref,2) = 0.                     !|  ||
              ReflexPhaseDiff(ni,i_Ref)  = 0.                     !|  ||
              if (Surface_Case(i_Ref)) Then !-------+              |  ||
                Gamma_Local = Real( Phi(i_Ref))    !|              |  ||
              else  !-------------------------------+              |  ||
                Gamma_Local = Abs(Gamma_Br(i_Ref)) !|              |  ||
              endif !-------------------------------+              |  ||
                                                                  !|  ||
c+------------------------------------------------+                |  ||
c| Averaging over the incident wave polarizations |                |  ||
c| if the incident radtiation is unpolarized      |                |  ||
c+------------------------------------------------+                |  ||
              do Ext_Polar = 1, Times_BounSolv !==============+    |  ||
c+--------------------------------------------+               |    |  ||
c|Summation of the reflection coefficient over|               |    |  ||
c|the polarizations of the diffracted wave:   |               |    |  ||
c+--------------------------------------------+               |    |  ||
              do   Int_Polar = 1, 2  !=====================+  |    |  ||
                i = 2*(i_Ref - 1) + Int_Polar             !|  |    |  ||
                                                          !|  |    |  ||
                dlt = (0.,0.)                             !|  |    |  ||
c+---------------------------------------------------+     |  |    |  ||
c|Summation over the eigenvalues (over the branches):|     |  |    |  ||
c+---------------------------------------------------+     |  |    |  ||
                if (Refl_Type(i_Ref).ne.'Laue')           !|  |    |  ||
     *                                   Then !------+     |  |    |  ||
                                                    !|     |  |    |  ||
                  if ((Surface_Case(i_Ref))         !|     |  |    |  ||
     *                       .AND.                  !|     |  |    |  ||
     *              i_Ref .eq. Incident_Place)      !|     |  |    |  ||
     *                                 Then !----+   |     |  |    |  ||
c+----------------------------+                  |   |     |  |    |  ||
c|For Incident Surface Reflex:| (specular)       |   |     |  |    |  ||
c+----------------------------+                  |   |     |  |    |  ||
                    do j = 1, Reduced_Size2 !=+  |   |     |  |    |  ||
                    dlt = dlt                !|  |   |     |  |    |  ||
     +                  + Fields(i,j)        !|  |   |     |  |    |  ||
     *                  * Lambda(j,Ext_Polar)!|  |   |     |  |    |  ||
     *                  * Epsilon(j)         !|  |   |     |  |    |  ||
     *                  * Renorm             !|  |   |     |  |    |  ||
     /                  / (2.*Gamma0)        !|  |   |     |  |    |  ||
                    enddo !===================+  |   |     |  |    |  ||
                                                !|   |     |  |    |  ||
                  else   !-----------------------+   |     |  |    |  ||
                                                !|   |     |  |    |  ||
c+---------------------------------------+       |   |     |  |    |  ||
c|For the other Bragg & Surface Reflexes:|       |   |     |  |    |  ||
c+---------------------------------------+       |   |     |  |    |  ||
                    do j = 1, Reduced_Size2 !=+  |   |     |  |    |  ||
                    dlt = dlt                !|  |   |     |  |    |  ||
     +                  + Fields(i,j)        !|  |   |     |  |    |  ||
     *                  * Lambda(j,Ext_Polar)!|  |   |     |  |    |  ||
                    enddo !===================+  |   |     |  |    |  ||
                  endif  !-----------------------+   |     |  |    |  ||
                else !-------------------------------+     |  |    |  ||
c+-------------------+                               |     |  |    |  ||
c|For Laue Reflexes: |                               |     |  |    |  ||
c+-------------------+                               |     |  |    |  ||
                  do j = 1, Reduced_Size2 !====+     |     |  |    |  ||
                    dlt = dlt                 !|     |     |  |    |  ||
     +                  + Fields(i, j)        !|     |     |  |    |  ||
     *                  * Lambda(j,Ext_Polar) !|     |     |  |    |  ||
     *                  * Expon(j)            !|     |     |  |    |  ||
                  enddo !======================+     |     |  |    |  ||
                endif  !-----------------------------+     |  |    |  ||
c+-----------------------------+                          !|  |    |  ||
c|   ATTENTION !!!             |                          !|  |    |  ||
c|Transmission coefficients for|                          !|  |    |  ||
c|surface beams are not compu- |                          !|  |    |  ||
c|ted in this version of BRL.  |                          !|  |    |  ||
c+-----------------------------+                          !|  |    |  ||
                                                          !|  |    |  ||
                Rc = (Gamma_Local/Gamma0) * Cdabs(dlt)**2 !|  |    |  ||
                                                          !|  |    |  ||
                ReflexCoeffPol(ni,i_Ref,Int_Polar) =      !|  |    |  ||
     =          ReflexCoeffPol(ni,i_Ref,Int_Polar) + Rc   !|  |    |  ||
                                                          !|  |    |  ||
                Phase = Atan2(Imag(dlt),Real(dlt)) / Dgra !|  |    |  ||phase in degrees
     -                        - ReflexPhaseDiff(ni,i_Ref) !|  |    |  ||
                if     (Phase.lt.-180.) Then !--+          |  |    |  ||
                  Phase = Phase+360.           !|          |  |    |  ||
                elseif (Phase.gt. 180.) Then !--+          |  |    |  ||
                  Phase = Phase-360.           !|          |  |    |  ||
                endif !-------------------------+          |  |    |  ||
                ReflexPhaseDiff(ni,i_Ref) = Phase         !|  |    |  ||
                                                          !|  |    |  ||
                ReflexCoeff(ni,i_Ref) =                   !|  |    |  ||
     =          ReflexCoeff(ni,i_Ref) + Rc                !|  |    |  ||
                                                          !|  |    |  ||
              enddo !======================================+  |    |  ||
c - End of loop over contributions of the polarizations of    |    |  ||
c   the diffracted wave                                       |    |  ||
                                                             !|    |  ||
c Restrict reflection coefficient by the                      |    |  ||
c minimum real*4 number!!!                                    |    |  ||
c This is a pity since we could do down to 1e-240,            |    |  ||
c but the other software (like grd2xxx) is not aimed          |    |  ||
c working with doubles (real*8):                              |    |  ||
              if (ReflexCoeff(ni,i_Ref) .lt.                 !|    |  ||
     *             ReflexCoeff_minimum)                      !|    |  ||
     *                ReflexCoeff(ni,i_Ref) = 0.             !|    |  ||
                                                             !|    |  ||
              enddo  !========================================+    |  ||
c - End of loop over the polarizations of the incident wave        |  ||
c                                                                  |  ||
              if (ReflexCoeff(ni,i_Ref).lt.0.  .or.               !|  ||
     *            ReflexCoeff(ni,i_Ref).gt.1.) Then !-----------+  |  ||
                if (ReflexCoeff(ni,i_Ref).lt.-1.D-3  .or.      !|  |  ||
     *              ReflexCoeff(ni,i_Ref).gt.1.+1.D-3) Then !-+ |  |  ||
                  nWrong = nWrong + 1                        !| |  |  ||
                  if (nWrong .gt. nWrongMax) goto 303        !| |  |  ||
c                 write (   *,  340,err=300)                 !| |  |  ||
c    *                          New_ord(i_Ref)-1,            !| |  |  ||
c    *                          n1,n2,                       !| |  |  ||
c    *                          ReflexCoeff(ni,i_Ref)        !| |  |  ||
                  write (luntbl,340,err=301)                 !| |  |  ||
     *                          New_ord(i_Ref)-1,            !| |  |  ||
     *                          n1,n2,                       !| |  |  ||
     *                          ReflexCoeff(ni,i_Ref)        !| |  |  ||
                endif !---------------------------------------+ v  v  vv
  340           format(
     *          ' ! Incorrect result for Reflection-',i2,
     *          '   -- at point=(',i3,',',i3,'):',
     *          '  Refl.Coeff.=',g13.6/
     *          ' ! Possible precission loss error',1x,
     *          '- Corrected!')                                !^  ^  ^^
                if (ReflexCoeff(ni,i_Ref).lt.0.)               !|  |  ||
     *              ReflexCoeff(ni,i_Ref) = 0.                 !|  |  ||
                if (ReflexCoeff(ni,i_Ref).gt.1.)               !|  |  ||
     *              ReflexCoeff(ni,i_Ref) = 1.                 !|  |  ||
                do  i = 1, 2  !===========================+     |  |  ||
                  if (ReflexCoeffPol(ni,i_Ref,i).lt.0.)  !|     |  |  ||
     *                ReflexCoeffPol(ni,i_Ref,i) = 0.    !|     |  |  ||
                  if (ReflexCoeffPol(ni,i_Ref,i).gt.1.)  !|     |  |  ||
     *                ReflexCoeffPol(ni,i_Ref,i) = 1.    !|     |  |  ||
                enddo  !==================================+     |  |  ||
              endif  !------------------------------------------+  |  ||
            enddo   !==============================================+  ||
c - End of loop over Bragg reflections                        v       ||
c                                                             |       ||
  555       Continue  !<-----if Gamma0<0------<---------------+       ||
c +----------------------------------+                                ||
c | Check for user hit "interrupt":  |                                ||
c +----------------------------------+                                ||
            Call KeyIn (0,nScan,nAscii)                       !Esc?   ||
            if (nscan .eq. 1)  goto 499 !-----------------------------++--+          ||
                                                                     !||  v
          enddo  !====================================================+|
c -- End of loop over Theta-1                                          |
c                                                                     !|
c+----------------------------------------------+                      |
c|Write reflection coefficients into .grd files:|                      |
c+----------------------------------------------+                      |
          if (iout.eq.0)  Then   !------------------------------+      |
            Finish = Seconds(Start) / 60.                      !|      |
            write (luntbl,339,err=300) n2,Tet(2),Finish        !v      v
  339       format(i5,'. Data written for all Tet1 &',1x,
     *           'Tet2=',g10.3,'   Elapsed time=',f6.2,' min') !^      ^
                                                               !|      |
            do    m = 1, nReflexes   !=====================+    |      |
              write (New_ord(m)+6,57,err=302)             !|    |      |
     *                               (ReflexCoeff(n1,m),  !|    |      |
     *                               n1=1,nMaxPoints(1))  !|    |      |
            enddo   !======================================+    |      |
          endif   !---------------------------------------------+      |
                                                                      !|
        enddo       !==================================================+
c -- End of loop over Theta-2

c A common representation of G format is: Gw.d[Ee], by default e=2.
c If e is specified, then w should be greater than or equal to d+e+5.
c 57    format(5001g12.3)
c 57    format(5001g12.3E3)                     !3+3+5=11
  57    format(5001g14.5E3)                     !5+3+5=13

c +------------------------------------+
c |Write .dat files (1D case; no GRD): |
c +------------------------------------+
        do i_Ref=1,nReflexes  !====+
          LineFile(i_Ref) = 0     !|
        enddo ! ===================+
        if (iout.gt.0)  Then   !---------------------------------------+
          do    n2=1,nMaxPoints(2) !=================================+ |
            Tet(2)   = Tet0(2) + dTet(2)*(n2-1)                     !| |
            do  n1=1,nMaxPoints(1) !==============================+  | |
              Tet(1) = Tet0(1) + dTet(1)*(n1-1)                  !|  | |
c Current matrix index:                                          !|  | |
                              ni=n1                              !|  | |
              if (iout.eq.2)  ni=n2                              !|  | |
              Dgamma0 = cSurf1*Tet(1)                            !|  | |
     +                + cSurf2*Tet(2)                            !|  | |
              Gamma0  = ( Gamma_Br(Incident_Place)               !|  | |
     +                +   Dgamma0 ) / UnitData(2)                !|  | |
              if (Gamma0.ge.0.)     Then  !--------------------+  |  | |
                do  i_Ref=1,nReflexes  !=====================+ |  |  | |
c Calculate alf(1:nReflexes):                                | |  |  | |
                  alf(i_Ref) = SNGL (                       !| |  |  | |
     (                             ca1(i_Ref)*Tet(1)        !| |  |  | |
     +                           + ca2(i_Ref)*Tet(2)        !| |  |  | |
     +                           + ca4(i_Ref)               !| |  |  | |
     -                           + ca3(i_Ref)               !| |  |  | |
     -                           * (Tet(1)**2+Tet(2)**2)    !| |  |  | |
     )                         )                            !| |  |  | |
                                                            !| |  |  | |
                  j = New_ord(i_Ref)+6                      !| |  |  | |
                                                            !| |  |  | |
                  if (LineFile(i_Ref).eq.0) Then !--------+  | |  |  | |
                    write (wrk,58,err=302) iout,i_Ref-1  !|  | |  |  | |
  58                format('#Theta(',i1,')',7x,
     *                     'ReflCoeff(',i2,')')          !|  | |  |  | |
                    Call PakInt (i_Ref-1,1,wrk(27:29),l) !|  | |  |  | |pack ReflCoeff(n) or ReflCoeff(nn)
                    wrk(27+l:27+l+1) = ') '              !|  | |  |  | |
                    l = Len_Trim(wrk)                    !|  | |  |  | |
                    if (print_polarizations) Then !---+   |  | |  |  | |
                       write(wrk(l+1:len(wrk)),55)   !|   |  | |  |  | |
  55                   format(2x,'RC(polariz1)',2x,
     *                           'RC(polariz2)',2x,
     *                 'PhaseDiff(pol2-pol1),degr.') !|   |  | |  |  | |
                       l = Len_Trim(wrk)             !|   |  | |  |  | |
                    endif !---------------------------+   |  | |  |  | |
                    write (j,'(a)',err=302) wrk(1:l)     !|  | |  |  | |
                  endif !---------------------------------+  | |  |  | |
                  LineFile(i_Ref) = LineFile(i_Ref) + 1     !| |  |  | |
                  if (print_polarizations) Then !------+     | |  |  | |
                    write (j,59,err=302)              !|     | |  |  | |
     *                     Tet(iout)                  !|     | |  |  | |
     *                    ,ReflexCoeff(ni,i_Ref)      !|     | |  |  | |
     *                    ,ReflexCoeffPol(ni,i_Ref,1) !|     | |  |  | |
     *                    ,ReflexCoeffPol(ni,i_Ref,2) !|     | |  |  | |
     *                    ,ReflexPhaseDiff(ni,i_Ref)  !|     | |  |  | |
c    *                    ,Gamma0,                    !|     | |  |  | |
c    *                    ,alf(i_Ref)                 !|     | |  |  | |
                  else !-------------------------------+     | |  |  | |
                    write (j,59,err=302)              !|     | |  |  | |
     *                     Tet(iout)                  !|     | |  |  | |
     *                    ,ReflexCoeff(ni,i_Ref)      !|     | |  |  | |
                  endif !------------------------------+     | |  |  | |
                enddo  !=====================================+ |  |  | |
              endif  !-----------------------------------------+  |  | |
            enddo  !==============================================+  | |
          enddo  !===================================================+ |
        endif  !-------------------------------------------------------+
c The common representation of G format is: Gw.d[Ee], by default e=2.
c If e is specified, then w should be greater than or equal to d+e+5.
c 59    format(2g12.3E3,10g14.5E3)                      !3+3+5=11
  59    format(100g14.5E3)                              !5+3+5=13
c#####################################################
        do    m = 1, nReflexes   !===================================+
          close(unit=m+6,err=101)                                   !|
  101     continue                                                  !|
          if (iout.eq.0)  Then  !----------------------------------+ |
c Produce DAT files from GRD (2D-case):                            | |
            wrk = OutFile(1:lout)//'_00.grd'                      !| |
            l = Len_Trim(wrk)                                     !| |
            write (wrk(lout+2:lout+3),'(i2)')  m-1                !| |
            if (wrk(lout+2:lout+2).eq.' ') wrk(lout+2:lout+2)='0' !| |
            i = System(grdprg(1:lgrd)//' /B '//wrk(1:l)           !| |
     *                                            //nude(1:lnul)) !| |
          endif !--------------------------------------------------+ |
        enddo !======================================================+
c#####################################################
        Finish = Seconds(Start) / 60.
c       write   (  *,   269)            Finish
        write   (luntbl,269,err=301)    Finish
  269   format  (';'/'; Working time: ',f8.2,' minutes')
        goto 500

c ######################################################
c #                     E R R O R S :                  #
c ######################################################
  503   continue
        l = Max(Len_Trim(wrk),1)
        write   (txt,995)       ifail, wrk(1:l), n1, n2
c    ,                          , Task_Size, ifail
  995   format(
c    *  1x,'Error ',i3,' from F02GJF'/
     *  1x,'Error ',i3,' from generalized eigenproblem.'/
     *  1x,a/
     *  1x,' - At point: ',i3,',',i3)
c    *  ' -- 30 *',i3,' iterations failed while trying to'/
c    *  ' -- isolate the ',i4,'th eigenvalue')
        lines = 3
        goto    991

c--------
  504   continue
        write   (txt,996)       ifail,n1,n2
  996   format  (
     *  1x,'Error ',i3,' from linear equation solution'/
     *  1x,'at point: (',i3,',',i3,')')
        lines = 2
        if (ifail .eq. 1) then !-----------------------------+
          txt(3) = 'The scattering matrix is singular.'     !|
          lines = 3                                         !|
        endif  !---------------------------------------------+
        goto    991

c--------
  507   continue
        l = Max(Len_Trim(wrk),1)
        write (txt,342) wrk(1:l), n1, n2
  342   format(
     *  1x,'! Insufficient work array for solving generalized',
     *  1x,'eigenproblem.'/
     *  1x,a/
     *  1x,' - At point (',i3,',',i3,')')
        lines = 3
        goto    991

c--------
  505   continue
        l = Max(Len_Trim(wrk),1)
        write (txt,343) wrk(1:l), n1, n2
  343   format(
     *  1x,'! Zero eigenvalues divider after solving generalized',
     *  1x,'eigenproblem.'/
     *  1x,a/
     *  1x,' - At point (',i3,',',i3,')')
        lines = 3
        goto    991

c--------
  300   continue
        txt(1) = 'Error writing output to console'
        lines = 1
        goto    991

c--------
  301   continue
        txt(1) = 'Error writing output to TBL-file'
        lines = 1
        goto    991

c--------
  302   continue
        txt(1) = 'Error writing output to GRD-file or DAT-file'
        lines = 1
        goto    991

c--------
  303   continue
        txt(1) = 'Sorry, algorithm aborted because of too many errors.'
        txt(2) = 'Possible reasons: too big deviations from the Bragg'
        txt(3) = 'condition or too thick crystal. The later may easily'
        txt(4) = 'happen especially when grazing (surface) waves are'
        txt(5) = 'present because the BRLS program does not incorporate'
        txt(6) = 'the "thick crystal" approximation -- it will be added'
        txt(7) = 'in the next release.'
        txt(8) = ' '
        txt(9) = 'Please, report the problem to the author.'
        lines = 9
        goto    991

c--------
  499   continue
        txt(1) = 'Work interrupted'
        lines = 1
        goto    991

c--------
  991   continue
        Call    Message (txt,lines,2)

c--------
  500   continue
        Call exit_quiet()
        End

c####################################################### 4a

        Subroutine Generalized_Eigenvalues (
     *                   Task_Size,             !number of equations
     *                   Left_Re,   Left_Im,    !left matrix A of A*X=lambda*B*X
     *                   Right_Re,  Right_Im,   !right matrix B of A*X=lambda*B*X
     *                   Epsilon_R, Epsilon_I,  !eigenvalues lambda of A*X=lambda*B*X (output)
     *                   D_Re,      D_Im,       !eigenvectors X of A*X=lambda*B*X (output)
     *                   Tolerance,             !tolerance to determine negligible elements
     *                   ifail,                 !exit status
     *                   msg)                   !exit message
c+===============================+
        Include 'brls_12.inc'   !|defines nDimen,nDimenTheta
c+===============================+
        Integer         nDimen4
        Parameter       (nDimen4 = 4*nDimen)

        Integer         Task_Size,
     *                  ifail,
     *                  i, j

        Character       msg*(*)

c       Logical         use_NAG /.True./
        Logical         use_NAG /.False./

        Real*8          Left_Re(nDimen4,Task_Size),
     *                  Left_Im(nDimen4,Task_Size),
     *                  Right_Re(nDimen4,Task_Size),
     *                  Right_Im(nDimen4,Task_Size),
     *                  D_Re(nDimen4,Task_Size),
     *                  D_Im(nDimen4,Task_Size),
     *                  Epsilon_R(Task_Size),
     *                  Epsilon_I(Task_Size),
     *                  Tolerance

c Work arrays for F02GJF:
c       Integer         Iter(nDimen4)                   !wrk for F02GJF
c       Real*8          beta_divider(nDimen4)           !wrk for F02GJF

c Work arrays for ZGGEV:
        Integer         lwork                           !should be at least 2*Task_Size,
c       Parameter       (lwork = 8*nDimen4)             !we use for ZGGEV  (see test at end) - too small!
        Parameter       (lwork = 8*nDimen4*(nDimen4+1)) !we use for ZGGEV3 (see test at end)

        Real*8          rwork(8*nDimen4)                !should be 8*Task_Size

        Complex*16      Left (nDimen4,nDimen4),         !wrk for ZGGEV (left matrix A)
     *                  Right(nDimen4,nDimen4),         !wrk for ZGGEV (right matrix B)
     *                  Alpha(nDimen4),                 !wrk for ZGGEV (lambda=alpha/beta)
     *                  Beta (nDimen4),                 !wrk for ZGGEV (lambda=alpha/beta)
     *                  dummy(1,1),                     !wrk for ZGGEV (left eigenvector)
c    *                  VL   (nDimen4,nDimen4),         !wrk for ZGGEV (right eigenvevtor)
     *                  VR   (nDimen4,nDimen4),         !wrk for ZGGEV (right eigenvevtor)
     *                  cwork(lwork),                   !wrk for ZGGEV (wrk array)
     *                  Lambda                          !wrk for ZGGEV (alpha/beta)

        ifail   = 0
        msg     = ' '

        i = nDimenTheta                                 !make compiler happy

        do i=1,Task_Size  !============+
           Epsilon_R(i) = 0.0D0       !|
           Epsilon_I(i) = 0.0D0       !|
           do j=1,Task_Size  !======+  |
              D_Re(i,j) = 0.0D0    !|  |
              D_Im(i,j) = 0.0D0    !|  |
           enddo  !=================+  |
        enddo  !=======================+

        if (use_NAG) Then !----------------------------------------------+
                                                                        !|
c //F02GJF//  calculates all the eigenvalues  and, if required,          |
c             all the eigenvectors of the complex generalized            |
c             eigenproblem Ax = lamda * Bx where A and B are complex,    |
c             square matrices, using the QZ algorithm.                   |
                                                                        !|
c F02GJF solves matrix equation A*X=lambda*B*X                           |
c F02GJF does not actually produce the eigenvalues lambda(j), but instead|
c returns alpha(j) and beta(j), such that lambda(j)=alpha(j)/beta(j).    |
c The division by beta(j) is the responsibility of your program, since   |
c beta(j) may be zero, indicating an infinite eigenvalue.                |
                                                                        !|
c Errors detected by the routine:                                        |
c IFAIL = I -- more than 30*N  iterations  have been performed           |
c              altogether in the second step of the QZ algorithm;        |
c              this limit has been reached while the  routine is         |
c              trying to isolate  the I(th)  eigenvalue.  On soft        |
c              failure, alpha(J) and beta(J) are correct for             |
c              J=I+1,...,N,  but the arrays VR and VI do not contain     |
c              any correct eigenvectors.                                 |
                                                                        !|
c          do i=1, Task_Size  !===========+                              |
c             beta_divider(i) = 0.0D0    !|                              |
c          enddo  !=======================+                              |
                                                                        !|
c          Call F02GJF (Task_Size,            !number of equations       |input
c    *                  Left_Re,  nDimen4,    !left matrix, Re(A)        |input, destroyed
c    *                  Left_Im,  nDimen4,    !left matrix, Im(A)        |input, destroyed
c    *                  Right_Re, nDimen4,    !rght matrix, Re(B)        |input, destroyed
c    *                  Right_Im, nDimen4,    !rght matrix, Im(B)        |input, destroyed
c    *                  Tolerance,            !tolerance neglect elements|input, if 0, will use machine precision
c    *                  Epsilon_R, Epsilon_I, !Re(alpha), Im(alpha)      |output
c    *                  beta_divider,         !beta                      |output
c    *                  .True.,               !Compute eigenvecors flag  |input
c    *                  D_Re, nDimen4,        !Re(X)                     |output
c    *                  D_Im, nDimen4,        !Im(X)                     |output
c    *                  Iter,                 !iterations needed for each|output
c    *                  ifail)                !exit status               |input/output
                                                                        !|
c          if (ifail.ne.0) return                                       !|
                                                                        !|
c          do i = 1,Task_Size  !=================================+       |
c             if (Abs(beta_divider(i)) .lt. 1.D-38) Then !----+  |       |
c                ifail = -1                                  !|  |       |
c                return                                      !|  |       |
c             endif  !----------------------------------------+  |       |
c             Epsilon_R(i) = Epsilon_R(i) / beta_divider(i)     !|       |
c             Epsilon_I(i) = Epsilon_I(i) / beta_divider(i)     !|       |
c          enddo  !==============================================+       |
                                                                        !|
        else !-----------------------------------------------------------+
                                                                        !|
c LAPACK ZGGEV computes the eigenvalues and, optionally, the left        |
c and/or right eigenvectors for GE matrices                              |
c https://www.netlib.org/lapack/complex16/zggev.f                        |
c https://www.netlib.org/lapack/explore-html/d3/d47/zggev_8f.html        |
                                                                        !|
c ZGGEV solves matrix equation A*X=lambda*B*X                            |
c ZGGEV does not actually produce the eigenvalues lambda(j), but instead |
c returns alpha(j) and beta(j), such that lambda(j)=alpha(j)/beta(j).    |
c The division by beta(j) is the responsibility of your program, since   |
c beta(j) may be zero, indicating an infinite eigenvalue.                |
                                                                        !|
           Lambda = Tolerance   !stop complains about "unused Tolerance" |
           do i=1, Task_Size  !=====================================+    |
           do j=1, Task_Size  !==================================+  |    |
              Left (i,j) = dcmplx(Left_Re (i,j), Left_Im (i,j)) !|  |    |
              Right(i,j) = dcmplx(Right_Re(i,j), Right_Im(i,j)) !|  |    |
           enddo  !==============================================+  |    |
           enddo  !=================================================+    |
                                                                        !|
c See also:                                                              |
c ZGGEV3 - works the same as ZGGEV, requires much larger lwork           |
c ZGGEVX - has more tuning parameters, have not tried                    |
           Call ZGGEV ('N',             !"N" not compute left eigenvectrs|input
     *                 'V',             !"V" do compute right eigenvectrs|input
     *                 Task_Size,       !number of equations             |input
     *                 Left,  nDimen4,  !left matrix A                   |input, destroyed
     *                 Right, nDimen4,  !right matrix B                  |input, destroyed
     *                 Alpha,           !lambda=alpha/beta (output)      |output
     *                 Beta,            !lambda=alpha/beta (output)      |output
     *                 dummy, 1,        !left eigenvectors X*A=lambda*X*B|output (unused)
c    *                 VL, nDimen4,     !left eigenvectors X*A=lambda*X*B|output (unused)
     *                 VR, nDimen4,     !rght eigenvectors A*X=lambda*B*X|output
     *                 cwork, lwork,    !complex wrk array, lwork >= 2*N |output/wrk
     *                 rwork,           !real wrk erray, 8*N             |output/wrk
     *                 ifail)           !exit status                     |input/output
                                                                        !|
c [ ifail < 0 ]       if ifail=-i, the i-th argument has an illegal value|
c [ ifail 1..Task_Size ] The QZ iteration failed. No eigenvectors have   |
c                        been calculated, but ALPHA(j) and BETA(j) should|
c                        be correct for j = ifail+1,...Task_Size         |
c [ ifail=Task_Size+1 ] Unexpected error returned from ZHGEQZ            !
c [ ifail=Task_Size+2 ] Error returned from ZTGEVC                       |
           if (ifail.ne.0) Then !-------------------------------------+  |
              if (ifail.lt.0) Then !--------------------------------+ |  |
                 write(msg,'("ZGGEV: illegal",i3,"-th argument.")')!| |  |
     *                               -ifail                        !| |  |
              elseif (ifail.le.Task_Size) Then !--------------------+ |  |
                 msg = 'QZ iteration failed in ZGGEV.'             !| |  |
              elseif (ifail.eq.Task_Size+1) Then !------------------+ |  |
                 msg = 'ZGGEV: Unexpected error from ZHGEQZ.'      !| |  |
              elseif (ifail.eq.Task_Size+2) Then !------------------+ |  |
                 msg = 'ZGGEV: Error returned from ZTGEVC.'        !| |  |
              endif  !----------------------------------------------+ |  |
              return                                                 !|  |
           endif  !---------------------------------------------------+  |
                                                                        !|
c On exit, if ifail=0, the real part of cwork(1) contains the            |
c minimum value of lwork required for optimal performance:               |
           i = Int(Dreal(cwork(1))+0.5)                                 !|
           if (i.gt.lwork) Then !-------------------------------------+  |
              write (msg,*) 'Needed=',i,' allocated=',lwork          !|  |
              j = Len_Trim(msg)                                      !|  |
              write (0,*) ' Generalized_Eigenvalues: insufficient',  !|  |
     *                   ' wrk array: ',msg(1:j)                     !|  |
              ifail = -999                                           !|  |
              return                                                 !|  |
           endif !----------------------------------------------------+  |
                                                                        !|
           do i=1,Task_Size  !========================================+  |
              if (Abs(Beta(i)) .lt. 1.D-38) Then !-----------------+  |  |
                 write(msg,'("Beta(",i2,") of ",i3," is zero.")') !|
     *                                i,     Task_Size            !|  |  |
                 ifail = -1000                                    !|  |  |
                 return                                           !|  |  |
              endif  !---------------------------------------------+  |  |
              Lambda = Alpha(i)/Beta(i)                              !|  |
              Epsilon_R(i) = Dreal(Lambda)                           !|  |
              Epsilon_I(i) = Dimag(Lambda)                           !|  |
              do j=1,Task_Size  !===============+                     |  |
                 D_Re(i,j) = Dreal(VR(i,j))    !|                     |  |
                 D_Im(i,j) = Dimag(VR(i,j))    !|                     |  |
              enddo  !==========================+                     |  |
           enddo  !===================================================+  |
                                                                        !|
        endif !----------------------------------------------------------+

        return
        end

c ####################################################### 4a

        Subroutine      Gauss_Eigenvalues  (Gauss_Matrix,       !left matrix         (input, destroyed)
     *                                      Incident_Field,     !right sides matrix  (input, destroyed)
     *                                      Lambda,             !solutions matrix    (output)
     *                                      Reduced_Size2,      !number of equations (input)
     *                                      Times_BounSolv,     !nmbr of right sides (input, 1 or 2)
     *                                      ifail)              !exit status         (output)
c+===============================+
        Include 'brls_12.inc'   !|defines nDimen,nDimenTheta
c+===============================+
        Integer         nDimen4
        Parameter       (nDimen4 = 4*Ndimen)

        Integer         Reduced_Size2,
     *                  Times_BounSolv,
     *                  ifail,
     *                  i, j

        Complex*16      Gauss_Matrix(nDimen4,Reduced_Size2),
     *                  Incident_Field(nDimen4,Times_BounSolv),
     *                  Lambda(nDimen4,Times_BounSolv)

c       Real*8          wkspce(nDimen4 )                !for NAG F04ADF

        Integer         IPIV(nDimen4)                   !for LAPACK ZGESV

c       Logical         use_NAG/.True./
        Logical         use_NAG/.False./

        ifail   = 0

        i = nDimenTheta                                 !make compiler happy

        if (use_NAG) Then !---------------------------------------------+
c //F04ADF// approximate solution of complex simultaneous linear        |
c            equations with multiple right-hand sides, using an LU      |
c            factorization with partial pivoting.                       |
                                                                       !|
c Possible errors detected by the routine:                              |
c IFAIL = 1 -- failure in  //F03AHF//,  the matrix A is singular,       |
c possibly  due to rounding errors.                                     |
                                                                       !|
c          Call F04ADF (Gauss_Matrix,   !left matrix                    |
c    *                  nDimen4,        !left matrix 1st dimension      |
c    *                  Incident_Field, !right sides matrix             |
c    *                  nDimen4,        !right sides matrx 1st dimension|
c    *                  Reduced_Size2,  !number of equations            |
c    *                  Times_BounSolv, !number of right sides          |
c    *                  Lambda,         !solutions matrix               |
c    *                  nDimen4,        !solutions matrix 1st dimension |
c    *                  wkspce,         !wrk array (number of equations)|
c    *                  ifail)          !exit status                   !|
                                                                       !|
        else  !---------------------------------------------------------+
                                                                       !|
c Use the LAPACK routine ZGESV which computes the solution to           |
c system of linear equations A*X=B for GE matrices (simple driver)      |
c https://www.netlib.org/lapack/explore-html/d1/ddc/zgesv_8f.html       |
                                                                       !|
           Call ZGESV (Reduced_Size2,  !number of equations             |
     *                 Times_BounSolv, !number of right sizes           |
     *                 Gauss_Matrix,   !left matrix                     |
     *                 nDimen4,        !left matrix 1st dimension       |
     *                 IPIV,           !wrk array (number of equations) |
     *                 Incident_Field, !right sides matrix              |
     *                 nDimen4,        !right sides matrix 1st dimension|
     *                 ifail)          !exit status                    !|
                                                                       !|
           if (ifail .eq. 0) Then !------------------------+            |
              do      i=1,Reduced_Size2 !===============+  |            |
              do      j=1,Times_BounSolv !===========+  |  |            |
                 Lambda(i,j) = Incident_Field(i,j)  !|  |  |            |
              enddo  !===============================+  |  |            |
              enddo  !==================================+  |            |
           endif !-----------------------------------------+            |
                                                                       !|
        endif !---------------------------------------------------------+

        return
        end
