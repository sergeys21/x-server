        subroutine      Case    (txt,icase)
c-----------------------------------------------
c This routine transforms the text string "txt"
c into the upper or lower case
c
c icase=0 - transform "txt" into the UPPER case
c icase=1 - transform "txt" into the LOWER case
c
c                 By Sergey Stepanov
c-----------------------------------------------
        integer         icase
        character       txt*(*)

c Transform ENGLISH letters:
        call    CaseLat (txt,icase)
c Transform RUSSIAN letters (in the CP866 coding):
        call    CaseRus (txt,icase)
        return
        end
