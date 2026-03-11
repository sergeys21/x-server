        Subroutine      GetUname (UnitCoef, UnitName)
c------------------------------------------------------------
c This program returns the name of units based on the value of
c Units Coefficient (the conversion factor from Units to Radians)
c
c                      Author: Sergey Stepanov
c------------------------------------------------------------
        Real            UnitCoef,
     *                  pi, AngDegree, AngMinute,
     *                  AngSec, AngMrad, AngMurad
        Character       Unitname*(*)

        pi        = 4.*atan(1.)
        AngDegree = 2.*pi/360.
        AngMinute = AngDegree/60.
        AngSec    = AngDegree/3600.
        AngMrad   = 1.E-03
        AngMurad  = 1.E-06

        if     (abs((UnitCoef-AngDegree)/AngDegree) .lt. 1.E-5) Then !--+
          UnitName = 'degr.'                                           !|
        elseif (abs((UnitCoef-AngMinute)/AngMinute) .lt. 1.E-5) Then !--|
          UnitName = 'min.'                                            !|
        elseif (abs((UnitCoef-AngSec)/AngSec)       .lt. 1.E-5) Then !--|
          UnitName = 'sec.'                                            !|
        elseif (abs((UnitCoef-AngMrad)/AngMrad)     .lt. 1.E-5) Then !--|
          UnitName = 'mrad.'                                           !|
        elseif (abs((UnitCoef-AngMurad)/AngMurad)   .lt. 1.E-5) Then !--|
          UnitName = 'urad.'                                           !|
        else  !---------------------------------------------------------|
          UnitName = ' '                                               !|
        endif  !--------------------------------------------------------+
        return
        end
