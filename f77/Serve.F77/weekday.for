        Integer Function WeekDay (date, month, year,
     *                            txdate, ENGday, RUSday)
c-----------------------------------------------------
c This subroutine determines week day by given date.
c-----------------------------------------------------
c weekday         - day number within week (integer).
c date,month,year - given month day, month number and year (integer).
c                   If all nulls are passed, then current date is taken.
c txdate*10       - date in text form .
c ENGday*12       - English name of week day.
c RUSday*12       - Russian name of week day.
c
c                Author: S.Stepanov
c-----------------------------------------------------
        Integer         date, month, year
c       Character       txdate*10,  ENGday*12,  RUSday*12
        Character       txdate*(*), ENGday*(*), RUSday*(*)

        Integer         startyear, delyear,
     *                  days, i
        Character       ENGweek(7)*12,
     *                  RUSweek(7)*12

        Integer          monthdays(12)
        Common /mondays/ monthdays

c Origin of date counting: January 1, 1973, Monday:
        Data    startyear       /1973/

c Names of week days:
        Data    ENGweek         /'Monday   ','Tuesday ',
     *                           'Wednesday','Thursday',
     *                           'Friday   ','Saturday',
     *                           'Sunday   '/
        Data    RUSweek         /'Понедельник','Вторник',
     *                           'Среда      ','Четверг',
     *                           'Пятница    ','Суббота',
     *                           'Воскресенье'/

        txdate = ' '
        ENGday = ' '
        RUSday = ' '
        weekday = 1
        if (Len(txdate) .lt. 10) Then !-------------------------------+
           write (*,*) 'Weekday: too short datestamp field passed'   !|
           return                                                    !|
        endif  !------------------------------------------------------+
        if (Len(ENGday) .lt. 12 .OR. Len(RUSday) .lt. 12) Then !------+
           write (*,*) 'Weekday: too short EN/RU day field passed'   !|
           return                                                    !|
        endif  !------------------------------------------------------+

c Number of days in each month:
        monthdays(1)    =       31
        monthdays(2)    =       28
        monthdays(3)    =       31
        monthdays(4)    =       30
        monthdays(5)    =       31
        monthdays(6)    =       30
        monthdays(7)    =       31
        monthdays(8)    =       31
        monthdays(9)    =       30
        monthdays(10)   =       31
        monthdays(11)   =       30
        monthdays(12)   =       31

c Current date:
        if (date.eq.0 .and. month.eq.0 .and. year.eq.0) Then !---+
           Call GetDat (year,month,date)                        !|
        endif !--------------------------------------------------+

c If it is leap year, then February has 29 days:
        if (year.eq.4*(year/4)) monthdays(2)=29

c Number of years since origin:
        delyear = year - startyear

c Number of days since origin with account for extra days in leap years:
        days = 365*delyear + delyear/4

c  Number of days since the beginning of the year:
        if (month.ne.1) Then  !--------------+
          do i=1,month-1  !===============+  |
            days = days + monthdays(i)   !|  |
          enddo  !========================+  |
        endif  !-----------------------------+
        days = days+date-1

c Find the week day index:
        weekday = days-7*(days/7)+1
        if (weekday.lt.1 .or. weekday.gt.7) stop 'WEEKDAY: fatal day'
c The date in text form:
        write   (txdate(1:2), '(i2)')  date
        write   (txdate(4:5), '(i2)')  month
        write   (txdate(7:10),'(i4)')  year
        txdate(3:3) = '.'
        txdate(6:6) = '.'
        if (txdate(4:4).eq.' ') txdate(4:4)='0'

        ENGday = ENGweek(weekday)
        RUSday = RUSweek(weekday)

        return
        end
