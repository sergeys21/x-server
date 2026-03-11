        Subroutine Expired (day_last, month_last, year_last)
c-----------------------------------------------------
c         Check for the given expiration date
c-----------------------------------------------------
c date,month,year - given day, month and year (all are integers)
c
c                   Author: S.Stepanov
c-----------------------------------------------------
        Integer         l, j
        Integer*2       day_last,  month_last,  year_last,
     *                  day_today, month_today, year_today,
     *                  idate(12), y2k/2000/
        Character       buffer*80

        Character       getOS*8
        External        getOS

        Logical         IsLate
        External        IsLate

c Get current date:
        Call    GetDat (year_today, month_today, day_today)
c The year is returned in the format like: 1998.
c However, after 2000 it might be something like: 0002
        if (year_today.lt.1900) year_today = year_today+y2k

        if (IsLate(day_last,  month_last,  year_last,
     *             day_today, month_today, year_today)) Then !-+
c       write (*,*) 'Current date:',                          !|
c    *              day_today, month_today, year_today        !|
        goto 2                                                !|
        endif  !-----------------------------------------------+

        if (getOS() .ne. 'windows') goto 1 !--------------------------+ return
                                                                     !|
c Now, test the dates of system files:                               !v
        Call  GetEnv32 ('WINDIR',buffer)
        l = Len_Trim (buffer)
        if (l.eq.0) goto 1
        do j=1,l  !================================+
          if (buffer(j:j).eq.'\') buffer(j:j)='/' !|
        enddo  !===================================+
        if (buffer(l:l).ne.'\' .AND.
     *      buffer(l:l).ne.'/') Then  !--+
          l = l+1                       !|
c         buffer(l:l) = '\'             !|
          buffer(l:l) = '/'             !|
        endif  !-------------------------+

c ...user.dat:
        Call GetFileDate_NN (buffer(1:l)//'user.dat',idate)
        day_today   = idate(3)
        month_today = idate(2)
        year_today  = idate(1)
        if (IsLate(day_last,  month_last,  year_last,
     *             day_today, month_today, year_today)) Then !-+
c       write (*,*) 'USER.dat date:',                         !|
c    *              day_today, month_today, year_today        !|
        goto 2                                                !|
        endif  !-----------------------------------------------+

c ...system.dat:
        Call GetFileDate_NN (buffer(1:l)//'system.dat',idate)
        day_today   = idate(3)
        month_today = idate(2)
        year_today  = idate(1)
        if (IsLate(day_last,  month_last,  year_last,
     *             day_today, month_today, year_today)) Then !-+
c       write (*,*) 'SYSTEM.dat date:',                       !|
c    *              day_today, month_today, year_today        !|
        goto 2                                                !|
        endif  !-----------------------------------------------+

c ...win386.swp:
        Call GetFileDate_NN (buffer(1:l)//'win386.swp',idate)
        day_today   = idate(3)
        month_today = idate(2)
        year_today  = idate(1)
        if (IsLate(day_last,  month_last,  year_last,
     *             day_today, month_today, year_today)) Then !-+
c       write (*,*) 'WIN386.SWP date:',                       !|
c    *              day_today, month_today, year_today        !|
        goto 2                                                !|
        endif  !-----------------------------------------------+
  1     return
  2     Call exit_quiet()
        end

c ==================================================================

        Logical Function IsLate(day_last,  month_last,  year_last,
     *                          day_today, month_today, year_today)

        Integer*2       day_last,  month_last,  year_last,
     *                  day_today, month_today, year_today

        Islate = .False.

        if (year_today.gt.year_last)         Then  !------+
          goto 1                                         !|
        elseif (year_today.eq.year_last)     Then  !------|
          if (month_today.gt.month_last)     Then  !--+   |
            goto 1                                   !|   |
          elseif (month_today.eq.month_last) Then  !--|   |
            if (day_today.ge.day_last)  goto 1       !|   |
          endif  !------------------------------------+   |
        endif  !------------------------------------------+
        return

  1     continue
        IsLate = .True.
        return
        end
