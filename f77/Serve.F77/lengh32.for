        integer function   lengh (txt)
c------------------------------------------------
c        Evaluation of text string length
c------------------------------------------------
        character txt*(*)

        lengh = Len_Trim(txt)

        return
        end
