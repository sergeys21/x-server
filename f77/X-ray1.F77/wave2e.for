        Real*8 Function wave2energy8 (wave8)

        Real*8    w2e
c https://www.physics.nist.gov/MajResFac/SURF/SURF/schwingersample.html
        Parameter (w2e=12.3984186)      !E=w2e/wave, NIST web site

        Real*8  wave8

        if (wave8 .gt. 0.) Then !----+
          wave2energy8 = w2e/wave8  !|
        else !-----------------------|
          wave2energy8 = 0.         !|
        endif  !---------------------+
        return
        end

c =====================================================================

        Real*4 Function wave2energy (wave)

        Real*8    w2e
c https://www.physics.nist.gov/MajResFac/SURF/SURF/schwingersample.html
        Parameter (w2e=12.3984186)      !E=w2e/wave, NIST web site

        Real*4  wave

        if (wave .gt. 0.) Then !----------+
          wave2energy = Sngl(w2e/wave)   !|
        else !----------------------------|
          wave2energy = 0.               !|
        endif !---------------------------+
        return
        end
