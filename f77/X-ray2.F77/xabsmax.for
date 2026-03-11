        subroutine xabsmax (x0, N_Total_Max, N_Total, UnitCoef, ! input
     *                      xabs, TERang, TER_angle)            ! output
 
        Integer    N_Total_Max, N_Total, i
        Complex*8  x0(N_Total_Max+1)
        Real*4     UnitCoef,            ! angle to radians
     *             xabs,                ! abs(x0_max)
     *             TERang,              ! critical TER angle (in radians)
     *             TER_angle,           ! critical TER angle (in UnitCoef)
     *             x
 
        xabs = 0.0
        do i=1,N_Total+1 !========================+
          x = abs(x0(i))                         !|
          if (x .gt. xabs) xabs = x              !|
        enddo  !==================================+

        TER_Angle = sqrt(xabs) / UnitCoef
        if (Abs(xabs) .lt. (1.E-24)) xabs = 1.   ! disable normalization
        TERang    = sqrt(xabs)
        return
        end
