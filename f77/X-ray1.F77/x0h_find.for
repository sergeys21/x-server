        Subroutine X0h_find (Code, abc, isyngony, wave,
     *                       mode_search, hkl_base, hkl_limits,
     *                       Q_limits, QB_limits, PRC_min,
     *                       nfoundMax, nfound, hkl_found,
     *                       Q_found, QB_found, PRC_found, ifail)
c -------------------------------------------------------
c    Search for Bragg reflections satisfying given conditions:
c -------------------------------------------------------
c+=========================================+
        Include         'x0h_size.inc'    !|
c+=========================================+
        Integer         hkl_base(4), hkl_limits(4,2), inmax,
     *                  mode_search, nfoundMax, nfound,
     *                  hkl_found(4,nfoundMax), ifail,
     *                  ind(4), hkl_step(3), ill(4,2),
     *                  isyngony, icall, i, iwarn, isico

        Real*8          qb8, q8, d8

        Real            abc(6), wave, Q_limits(2), QB_limits(2),
     *                  PRC_min, Q_found(nfoundMax),
     *                  QB_found(nfoundMax), PRC_found(nfoundMax),
     *                  qb, qb1, qb2, q, q1, q2, q3, q4,
     *                  d, s, prc

        Character       code*20

        Character       txt(20)*80
        Common  /msg/   txt

        Real*8          Bragan8, DAngle
        External        Bragan8, DAngle
        Real            Xrh_Rel
        External        Xrh_Rel

c -------------------------------------------------------
        prc = 0.                        !keep GNU Fortran happy
        i = kcompMax                    !to suppress unused paramters warnings
        i = kolMax                      !to suppress unused paramters warnings
        i = ncodMax                     !to suppress unused paramters warnings
        i = natmMax                     !to suppress unused paramters warnings
        i = nwavMax                     !to suppress unused paramters warnings
        i = nedgMax                     !to suppress unused paramters warnings
        i = nabcMax                     !to suppress unused paramters warnings
        i = nedgTb                      !to suppress unused paramters warnings
        i = maxX0hWarnLines             !to suppress unused paramters warnings
c -------------------------------------------------------
        nfound = 0
        if (isyngony .ne. 2) Then !---+
          inmax = 3                  !|
        else !------------------------+
          inmax = 4                  !|
        endif  !----------------------+

        if (hkl_base(1)    .eq. 0  .and.
     *      hkl_base(2)    .eq. 0  .and.
     *      hkl_base(inmax).eq. 0)  Then !----------------+
          txt(1) = 'X0h_find: no base plane specified'   !|
          ifail = 1                                      !|
          goto 100                                       !|
        endif  !------------------------------------------+

        if (mode_search .lt. 1 .or. mode_search .gt. 3) Then !--+
          write (txt,1) mode_search                            !|
  1       format('X0h_find: incorrect search mode =',i3/
     *    'Valid modes are 1 to 3')                            !|
          ifail = 2                                            !|
          goto 100                                             !|
        endif  !------------------------------------------------+


        if (hkl_limits(1,2) .ge. hkl_limits(1,1)) Then !--------+
          hkl_step(1) = 1                                      !|
        else !--------------------------------------------------+
          hkl_step(1) =-1                                      !|
        endif  !------------------------------------------------+
        if (hkl_limits(2,2) .ge. hkl_limits(2,1)) Then !--------+
          hkl_step(2) = 1                                      !|
        else !--------------------------------------------------+
          hkl_step(2) =-1                                      !|
        endif  !------------------------------------------------+
        if (hkl_limits(inmax,2).ge.hkl_limits(inmax,1)) Then !--+
          hkl_step(3) = 1                                      !|
        else !--------------------------------------------------+
          hkl_step(3) =-1                                      !|
        endif  !------------------------------------------------+

c Min/max for Bragg indices scanning:
        do i = 1,inmax  !====================================+
          ill(i,1) = Min0(hkl_limits(i,1),hkl_limits(i,2))  !|
          ill(i,2) = Max0(hkl_limits(i,1),hkl_limits(i,2))  !|
        enddo  !=============================================+
c Min/max for Bragg angles:
        qb1 = Amin1(QB_limits(1),QB_limits(2))
        qb2 = Amax1(QB_limits(1),QB_limits(2))
c Min/Max for inter-planar angles:
        q1 = Amin1(Q_limits(1),Q_limits(2))
        q2 = Amax1(Q_limits(1),Q_limits(2))

c------------------------------------- Search loop:
        ind(1)     = hkl_limits(1,1)
  651   continue   !<------------------------------------------+
        ind(2)     = hkl_limits(2,1)                          !|
  652   continue   !<----------------------------------------+ |
        ind(inmax) = hkl_limits(inmax,1)                    !| |
  653   continue   !<--------------------------------------+ | |
        if (isyngony.eq.2) Then !---+                      | | |
          ind(3)=-(ind(1)+ind(2))  !|                      ^ ^ ^
        endif  !--------------------+                      | | |
                                                          !| | |
        if (ind(1)    .eq. 0  .and.                       !| | |
     *      ind(2)    .eq. 0  .and.                       !| | |
     *      ind(inmax).eq. 0)  goto 654 !--(0,0,0)-->----+ | | |
                                                        !| | | |
        iwarn = 0               !take QB=0 if l>2d      !| | | |
        isico = 0               !skip calc of sin/cos   !| | | |
        qb8 = Bragan8 (wave,ind,inmax,iwarn,            !| | | |
     *               d8,isico,s,s,s,s,s,s)              !| | | |
        qb = Sngl(qb8)                                  !| | | |
        d  = Sngl(d8)                                   !| | | |
c If no Bragg reflection, go to next:                    | | | |
        if (Abs(qb).lt.1.E-37)          goto 654  !-->---+ | | |
c If Bragg angle is not in range, go to next:            v ^ ^ ^
        if (qb.lt.qb1 .or. qb.gt.qb2)   goto 654  !-->---+ | | |
                                                        !| | | |
        q8 = DAngle (ind,hkl_base,inmax)                !v ^ ^ ^
        q  = Sngl(q8)                                   !| | | |
c If reflection is from the same plane:                  | | | |
c cc Commented by Sergey on July 20, 2001:               | | | |
c cc    if (q.eq.0.)                    goto 654  !-->---+ | | |
                                                        !| | | |
        if     (mode_search .eq. 1) Then  !---------+    | | | |
c From QB-Theta1 to QB-Theta2:                      |    | | | |
          q3 = Amin1(qb-q1,qb-q2)                  !|    v ^ ^ ^
          q4 = Amax1(qb-q1,qb-q2)                  !|    | | | |
          if (q.lt.q3 .or. q.gt.q4)     goto 654 !--+->--+ | | |
        elseif (mode_search .eq. 2) Then  !---------+    | | | |
c From Theta1 to QB-Theta2:                         |    | | | |
          q3 = Amin1(Q_limits(1),qb-Q_limits(2))   !|    v ^ ^ ^
          q4 = Amax1(Q_limits(1),qb-Q_limits(2))   !|    | | | |
          if (q.lt.q3 .or. q.gt.q4)     goto 654 !--+->--+ | | |
        elseif (mode_search .eq. 3) Then  !---------+    | | | |
c From Theta1 to Theta2:                            |    | | | |
          if (q.lt.q1 .or. q.gt.q2)     goto 654 !--+->--+ | | |
        endif  !------------------------------------+    | | | |
                                                        !v ^ ^ ^
        if (PRC_min .gt. 0.) Then  !-----------------+   | | | |
c This presumes that X0h1 has been called and we     |   | | | |
c already have abc(6), wave, code, isyngony:         |   | | | |
          icall = 5                                 !|   | | | |
          prc = 100.*Xrh_Rel(abc, wave, code,       !|   | | | |
     *                       ind(1), ind(2),        !|   | | | |
     *                       ind(inmax),            !|   | | | |
     *                       isyngony, icall,       !|   | | | |
     *                       ifail)                 !|   v | | |
          if (ifail.ne.0) goto 99  !-----------------+---+-+-+-+---+
          if (prc.lt.PRC_min)    goto 654  !---------+->-+ ^ ^ ^   v
        endif  !-------------------------------------+   | | | |
                                                        !v | | |
        nfound = nfound+1                               !| | | |
        if (nfound .le. nfoundMax) Then  !---+           | | | |
          do i=1,inmax  !==================+ |           | | | |
            hkl_found(i,nfound) = ind(i)  !| |           | | | |
          enddo  !=========================+ |           | | | |
          Q_found(nfound)   = q             !|           | | | |
          QB_found(nfound)  = qb            !|           | | | |
          if (PRC_min .gt. 0.)              !|           | | | |
     *    PRC_found(nfound) = prc           !|           | | | |
        endif !------------------------------+          !| | | |
                                                        !| | | |
  654   continue   !<------------------------------------+ | | |
        ind(inmax) = ind(inmax)+hkl_step(3)               !| | |
        if (ind(inmax) .ge. ill(inmax,1) .and.            !| | |
     *      ind(inmax) .le. ill(inmax,2))      goto 653  !-+ | |
        ind(2) = ind(2)+hkl_step(2)                         !| |
        if (ind(2)     .ge. ill(2,1)     .and.              !| |
     *      ind(2)     .le. ill(2,2))          goto 652  !---+ |
        ind(1) = ind(1)+hkl_step(1)                           !|
        if (ind(1)     .ge. ill(1,1)     .and.                !|
     *      ind(1)     .le. ill(1,2))          goto 651  !-----+

  100   continue
        return

  99    continue
        ifail = iabs(ifail)
        goto 100
        end
