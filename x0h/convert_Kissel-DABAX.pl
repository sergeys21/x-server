#!/usr/bin/perl -w
#--------------------------------------------------------------
# This script converts data from f1f2_asf_Kissel.dat which contains data
# in the DABAX project form: 
# {E_KeV, g", g', f', (f+df1), df2 } => into E_eV, (f+df1), df2
#
#  #########################################################
#  # Something is wrong here:                              #
#  # ------------------------                              #
#  # According to f1f2_asf_Kissel.dat, g"=-df2, but in few #
#  # cases it gives negative df2 for La, Ac, and Th, i.e.  #
#  # NEGATIVE absorption at some energies. Those data      #
#  # points are skipped (replaced with $nodata).           #
#  #########################################################
#
# The Kissel data are based on these two publications:
#
#
# [1] Elastic Scattering of Gamma-Rays and X-Rays by Atoms, by
# P. P. Kane, Lynn Kissel, R. H. Pratt and S. C. Roy, Physics Reports
# Vol. 140, 75-159 (1986).
#
# [2] Rayleigh Scattering - Elastic Photon Scattering by Bound
# Electrons, by Lynn Kissel and R. H. Pratt, in Atomic Inner-Shell
# Physics, edited by Bernd Crasemann (Plenum Publishing: New York,
# 1985).
#
#
#			Author: Sergey Stepanov
#
# Version-1:    2021/06/25
# Version-2:    2026/05/05 adapted to extra column is the DABAX file.
#--------------------------------------------------------------
  use strict;
  use warnings 'all';
# use warnings::unused;

  select STDOUT; $|=1;							#set unbuffered output
  select STDOUT; $|=1;							#set unbuffered output

  my ($X0H, $script, $status);
  my ($inp, $out, $lbl, $nodata, $i, $crystal);
  my (@lines, $str, @words, $nwords);

  $script = $0 =~ s/^.*[\/\\]//gr;					#remove everything before last "/" or "\"

  $inp = 'f1f2_asf_Kissel.dat';						#DABAX-format input file
  if (! -e $inp) {
     $str = 'wget -O '.$inp.' https://ftp.esrf.fr/pub/scisoft/xop2.3/DabaxFiles/'.$inp;
     $status = system($str);
     if ($status) {die '*** '.$script.': failed downloading '.$inp.' with:'."\n".$str."\n";}
  }
  $out  = 'kissel.dat';							#X0H-formatted output file
  $nodata = -.999900;

  if (defined $ENV{X0H})	{$X0H=$ENV{X0H};}
  else				{$X0H=dirname(abs_path($0));}
  $X0H =~ s/\\/\//g;
  if    (-e $X0H.'/support_subs.pl') {do $X0H.'/support_subs.pl';}
  else  {die $0.': cannot find support_subs.pl for X0H='.$X0H."\n";}

  if (! -e $inp) {
     die '*** '.$script.': No input file ['.$inp.'] in the current directory.'."\n";
  }
  
  $status = &read_entire_file($inp,\@lines);
  if ($status) {
     die '*** '.$script.': Error reading input file ['.$inp.'].'."\n";
  }

### Open output file for writing:
  open (DAT,'>',$out) || die 'Cannot open '.$out."\n";			#overwrite DAT file
  select DAT; $|=1;							#set unbuffered output
  $i = 0;
  foreach $str (@lines) {
     $i++;
     $lbl = "#UO       E (keV)          g''           g'           f'           f1";
     if ($str =~ /\Q${lbl}\E/)    {print DAT '#UO       E (eV)           f1            f2',"\n"; next;}
     $lbl = "2:g''[=-f''CL=-f2]  3:g' 4:f'  5:f1";
     if ($str =~ /\Q${lbl}\E/)        {next;}
     $lbl = "PhotonEnergy[KeV]  -f2  g'  f'  f1";
     if ($str =~ /\Q${lbl}\E/)        {next;}
     $lbl = "#UO  *** END OF DATA ***";
     if ($str =~ /\Q${lbl}\E/)        {next;}
     if ($str =~ /^\#UO\s+INPUT/)     {next;}
     if ($str =~ /^\#UO\s+OUTPUT/)    {next;}
     if ($str =~ /^\#UO\s+TITLE /)    {next;}
     if ($str =~ /^\#UO\s+# COLUMNS/) {next;}
     if ($str =~ /^\#UO\s+# UNIQUE /) {next;}
     if ($str =~ /^\#UO\s+REMARK/)    {next;}
     if ($str =~ /^\#UO\s+ASF:/)      {next;}
     if ($str =~ /^\#UO\s+LABEL/)     {next;}
     if ($str =~ /^\#N /)             {next;}
     if ($str =~ /^\#L /)             {next;}
     if ($str =~ /^\#UF1ADD /)        {next;}
     if ($str =~ /^\#S /) {
        $crystal = $str;
        $crystal =~ s/^\#S\s+//;
     }
     if ($str =~ /^\#/ || length($str) == 0) {print DAT $str,"\n"; next;}
     $str =~ s/^\s+//;
     @words = split(/\s+/,$str);
     $nwords = @words;
     if ($nwords != 6) {
        die '*** '.$script.': incorrect record ['.$str.'] in data file ['.$inp.'] -- nWords='.$nwords."\n";
     }
     if (! &is_float($words[0]) || ! &is_float($words[1]) || ! &is_float($words[4]) || ! &is_float($words[5])) {
        die '*** '.$script.': incorrect record ['.$str.'] in data file ['.$inp.'] -- non-float data. Line='.$i.', crystal='.$crystal."\n";
     }
     elsif ($words[5] < 0) {
        print STDOUT '!!! '.$script.': incorrect record ['.$str.'] in data file ['.$inp.'] -- negative absorption. Line='.$i.', crystal='.$crystal."\n";
     }
     elsif ($words[5] != -$words[1]) {
        die '*** '.$script.': incorrect record ['.$str.'] in data file ['.$inp.'] -- g2 # -f2. Line='.$i.', crystal='.$crystal."\n";
     }
     if ($words[0] <= 0) {next;}
     if ($words[5] <= 0) {$words[5] = $nodata;}
     printf DAT ' %14.5f   %-13.6g   %-13.6g'."\n", 1000*$words[0],$words[4],$words[5];    #from KeV to eV
  }
  close(DAT);
  print STDOUT "\n".$script.': All done.'."\n";
  exit 0;
