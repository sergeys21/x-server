The programs and files in this distribution constitute a backend of 
X-ray Server at https://x-server.gmca.aps.anl.gov

They are licensed under GNU General Public License v3.0 (GPL-3.0).

Below are the steps to make these programs compiled and working on your 
computer.

Install needed software: GNU Fortran, CodeBlocks IDE and LAPACK. For example,
on Ubuntu, Debian and derivative operating systems the installation command 
will be:
$ sudo apt install gfortran codeblocks liblapack-dev

Add X0H environment (e.g. for ~/.bashrc):
$ export X0H=<path_to_x0h_folder>
$ export PATH=$PATH:$X0H

Start CodeBlocks:
$ codeblocks

In CodeBlocks:
 - At startup choose "GNU FORTRAN Compiler" as default compiler
 - Go to Settings -> Compiler -> Linker -> add "lapack" library
 - Go to File -> Open -> open f77/CodeBlocks/_xserver/xserver.workspace
 - Go to Build -> Build Workspace
 
Alternatively (not using CodeBlocks):
  cd f77/CodeBlocks/_xserver
  make

In either case (CodeBlocks or make),Install compiled executables into the x0h directory:

$ f77/CodeBlocks/xserver_install_linux.pl

The main programs are (for documentation and details see the X-ray Server 
webpage at https://x-server.gmca.aps.anl.gov):
x0h_web   - calculates X-ray scattering factors
x0p_web   - provides search for Bragg planes
gid_slm7  - calculates X-ray Bragg diffraction from crystals with optional surface multilayers
ter_slm5  - calculates X-ray specular reflection from mirrors and multilayers
brl_build - builds X-ray multiple Bragg diffraction geometry for brls12
brls12    - calculates X-ray multiple Bragg diffraction from crystals
trds_99   - calculates X-ray diffuse scattering from interface roughness in mirrors and multilayers
mag_sl98  - calculates X-ray resonant scattering from multilayers
mag_sl99  - calculates X-ray resonant scattering from multilayers (simplified model for hard X-rays)

Programs x0h_web, x0p_web, and brl_build take input through their command 
line in the form: key1=value key2=vallue ...
Run them without parameters to see all available keys.

The rest of the programs need an input file to be prepared and specified
in the command line. See input files templates per each program in the 
x0h directory.

The remaining programs (grd2dat, grd2gnu, grd2grd, cofu_dec) are auxiliary.

- Sergey Stepanov, 2026/03/03
