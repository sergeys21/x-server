        Subroutine      DateStam (txdate,engday)

        Character       txdate*10, engday*12
        Character       rusday*12
        Integer         iday, month, iyear, iweek

        Integer         WeekDay
        External        WeekDay

        iday  = 0
        month = 0
        iyear = 0
        iweek = WeekDay (iday,month,iyear,txdate,engday,rusday)
        return
        end
