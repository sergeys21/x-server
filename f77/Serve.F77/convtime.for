        Subroutine Convert_Time (hh,mm,ss,timestamp)
        Integer*4       hh,mm,ss
        Character*9     timestamp
        write   (timestamp,1) hh,mm,ss
  1     format  (i3,2(':',i2))
        if (timestamp(5:5).eq.' ') timestamp(5:5)='0'
        if (timestamp(8:8).eq.' ') timestamp(8:8)='0'
        return
        end

