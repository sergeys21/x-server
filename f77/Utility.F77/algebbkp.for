c This subroutines Backup and Restore Common blocks for
c Algebra, Algebra2 and Algebra3

c =====================================================
        Subroutine      LatticeBackup ()
c-----------------------------------------------------
        Integer lbev4
        Real    a4(6), vol4
        Common /back1/ a4, vol4, lbev4

        Real    at4(6), g4(3,3), gt4(3,3)
        Common  /back2/ at4, g4, gt4

        Integer lbev8
        Real*8  vol8
        Real*8  at8(6), g8(3,3), gt8(3,3)
        Common /Dback2/ at8, g8, gt8, vol8, lbev8

        Integer lbev4_bkp, lbev8_bkp, i, j, k
        Real    BKP4(1+6+6+9+9)
        Real*8  BKP8(1+6+9+9)
        Common  /LatticeBKP/    BKP8, BKP4, lbev4_bkp, lbev8_bkp

        if (lbev4 .ne. 0) Then !----------------------+
          BKP4(1) = vol4                             !|
          do i=1,6 !============+                     |
            BKP4(1+i) = a4(i)  !|                     |
            BKP4(7+i) = at4(i) !|                     |
          enddo !===============+                    !|
          k = 13                                     !|
          do i=1,3  !===================+             |
            do j=1,3  !===============+ |             |
              BKP4(k+j)   = g4(i,j)  !| |             |
              BKP4(k+9+j) = gt4(i,j) !| |             |
            enddo  !==================+ |             |
            k = k+3                    !|             |
          enddo  !======================+             |
          lbev4_bkp = lbev4                          !|
          if (lbev8 .ne. 0) Then !----------------+   |
            BKP8(1) = vol8                       !|   |
            do i=1,6  !============+              |   |
              BKP8(1+i) = at8(i)  !|              |   |
            enddo  !===============+              |   |
            k = 7                                !|   |
            do i=1,3  !=====================+     |   |
              do j=1,3  !================+  |     |   |
                BKP8(k+j)   = g8(i,j)   !|  |     |   |
                BKP8(k+9+j) = gt8(i,j)  !|  |     |   |
              enddo  !===================+  |     |   |
              k = k+3                      !|     |   |
            enddo  !========================+     |   |
            lbev8_bkp = lbev8                    !|   |
          endif !---------------------------------+   |
        endif !---------------------------------------+
        return
        end


c =====================================================
        Subroutine      LatticeRestore ()
c-----------------------------------------------------
        Integer lbev4
        Real    a4(6), vol4
        Common /back1/ a4, vol4, lbev4

        Real    at4(6), g4(3,3), gt4(3,3)
        Common  /back2/ at4, g4, gt4

        Integer lbev8
        Real*8  vol8
        Real*8  at8(6), g8(3,3), gt8(3,3)
        Common /Dback2/ at8, g8, gt8, vol8, lbev8

        Integer lbev4_bkp, lbev8_bkp, i, j, k
        Real    BKP4(1+6+6+9+9)
        Real*8  BKP8(1+6+9+9)
        Common  /LatticeBKP/    BKP8, BKP4, lbev4_bkp, lbev8_bkp

        if (lbev4_bkp .ne. 0) Then !------------------+
          vol4 = BKP4(1)                             !|
          do i=1,6 !============+                     |
            a4(i)  = BKP4(1+i) !|                     |
            at4(i) = BKP4(7+i) !|                     |
          enddo !===============+                    !|
          k = 13                                     !|
          do i=1,3  !===================+             |
            do j=1,3  !===============+ |             |
              g4(i,j)  = BKP4(k+j)   !| |             |
              gt4(i,j) = BKP4(k+9+j) !| |             |
            enddo  !==================+ |             |
            k = k+3                    !|             |
          enddo  !======================+             |
          lbev4_bkp = lbev4                          !|
          if (lbev8 .ne. 0) Then !----------------+   |
            vol8 = BKP8(1)                       !|   |
            do i=1,6  !============+              |   |
              at8(i) = BKP8(1+i)  !|              |   |
            enddo  !===============+              |   |
            k = 7                                !|   |
            do i=1,3  !=====================+     |   |
              do j=1,3  !================+  |     |   |
                g8(i,j)  = BKP8(k+j)    !|  |     |   |
                gt8(i,j) = BKP8(k+9+j)  !|  |     |   |
              enddo  !===================+  |     |   |
              k = k+3                      !|     |   |
            enddo  !========================+     |   |
            lbev8 = lbev8_bkp                    !|   |
          else  !---------------------------------|   |
            lbev8 = 0   !erase & recalc backlattic|   |
          endif !---------------------------------+   |
        else  !---------------------------------------|
          stop 'LatticeRestore: no prev.lattice'     !|
        endif !---------------------------------------+
        return
        end
