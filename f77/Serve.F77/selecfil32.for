c ATTENTION: this a dummy emulator for For32!

        Subroutine    Selecfil (Dir_Mask,
     *                          Selected_Files,
     *                          Selected_Indexes,
     *                          Max_Selected,
     *                          Num_Selected,
     *                          Multi_Flag,
     *                          Nx,Ny,Ncolumns,Nlines)
c-----------------------------------------------------
c This subroutine provides file choice on the screen menu
c with the help of arrow keys and the "Ins" and "Del" keys.
c
c               Author: S.Stepanov
c-----------------------------------------------------
        integer         Max_Selected, Num_Selected,
     *                  Multi_Flag, Nx, Ny,
     *                  Ncolumns, Nlines,
     *                  Selected_Indexes(Max_Selected)
        character       Dir_Mask*(*),
     *                  Selected_Files(Max_Selected)*(*)

        integer         i
        character       byte

        byte = Dir_Mask         !to suppress compiler complains about unused vars
        i    = Multi_Flag       !to suppress compiler complains about unused vars
        i    = Nx               !to suppress compiler complains about unused vars
        i    = Ny               !to suppress compiler complains about unused vars
        i    = Ncolumns         !to suppress compiler complains about unused vars
        i    = Nlines           !to suppress compiler complains about unused vars

        Num_Selected = 0
        do i=1,Max_Selected  !==================+
          Selected_Files(i)   = ' '            !|
          Selected_indexes(i) = 0              !|
        enddo  !================================+
        return
        end
