#!/usr/bin/perl -w
#--------------------------------------------------------------
# This script converts data from multiple files mac[NN].csv of 
# the Github project:
#
#           https://github.com/usnistgov/FFAST.jl
#
# which contains Chantler's X-ray scattering data and converts
# it into a form acceptable by X0H and similar to the Henke 
# tables of f1, f2.
#
# The Chantler data are based on these two publications:
#
# [1] C.T. Chantler, J. Phys. Chem. Ref. Data 29(4), 597-1048 (2000). 
# "Detailed Tabulation of Atomic Form Factors, Photoelectric Absorption 
# and Scattering Cross Section, and Mass Attenuation Coefficients in the 
# Vicinity of Absorption Edges in the Soft X-Ray (Z = 30-36, Z = 60-89, 
# E = 0.1 keV-10 keV), Addressing Convergence Issues of Earlier Work."
# 
# [2] C.T. Chantler, J. Phys. Chem. Ref. Data 24, 71 (1995). 
# "Theoretical Form Factor, Attenuation and Scattering Tabulation for 
# Z=1-92 from E=1-10 eV to E=0.4-1.0 MeV."
#
# The details are available at NIST:
#
# Chantler, C., Olsen, K., Dragoset, R., Kishore, A., Kotochigova, S. and Zucker, D. (2003), 
# X-Ray Form Factor, Attenuation and Scattering Tables (version 2.0), 
# http://physics.nist.gov/ffast,
# 
# https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables
# 
# Web input form:
# https://physics.nist.gov/PhysRefData/FFast/html/form.html
#
#
#			Author: Sergey Stepanov
#
# Version-1:    2021/06/22
#--------------------------------------------------------------
  use strict;
  use warnings 'all';
# use warnings::unused;

  select STDOUT; $|=1;							#set unbuffered output
  select STDOUT; $|=1;							#set unbuffered output

  my ($X0H, $script, @elements, $code, $status);
  my ($zmin, $zmax, $z, $path, $csv, $out, $j);
  my (@lines, @words, $nwords, $str, $nodata);

  $script = $0 =~ s/^.*[\/\\]//gr;					#remove everything before last "/" or "\"

  $zmin = 1;
  $zmax = 92;
  $path = 'data/';							#filenames are mac[NN].csv. Place them into this directory
  $out  = 'chantler.dat';						#output file
  $nodata = -.999900;

  if (defined $ENV{X0H})	{$X0H=$ENV{X0H};}
  else				{$X0H=dirname(abs_path($0));}
  $X0H =~ s/\\/\//g;
  if    (-e $X0H.'/support_subs.pl') {do $X0H.'/support_subs.pl';}
  else  {die $0.': cannot find support_subs.pl for X0H='.$X0H."\n";}

  if (&import_elements($X0H.'/elements.dat',\@elements)) {exit;}

  if (! -e $path) {
     die '*** '.$script.': No directory ['.$path.'] with the Chantler data mac[NN].csv'."\n";
  }
  
### Open output file for writing:
  open (DAT,'>',$out) || die 'Cannot open '.$out."\n";			#overwrite DAT file
  select DAT; $|=1;							#set unbuffered output
  print DAT '#F f1f2_Chantler.dat'."\n"
  .'#UT Elastic Photon-Atom Scattering, anomalous scattering factors.'."\n"
  .'#C  This file has been created using '.$script.' on '.scalar(localtime)."\n"
  .'#C  It is a part of X-ray Server, https://x-server.gmca.aps.anl.gov'."\n"
  .'#UD '."\n"
  .'#UD The scattering values in this file have been obtained from:'."\n"
  .'#UD https://github.com/usnistgov/FFAST.jl/tree/master/data/mac[NN].csv'."\n"
  .'#UD and a WWW interface to these data is at:'."\n"
  .'#UD https://physics.nist.gov/PhysRefData/FFast/html/form.html'."\n"
  .'#UD '."\n"
  .'#UD Reference'."\n"
  .'#UD   [1] C.T. Chantler, J. Phys. Chem. Ref. Data 29(4), 597-1048 (2000).'."\n"
  .'#UD   "Detailed Tabulation of Atomic Form Factors, Photoelectric Absorption'."\n"
  .'#UD   and Scattering Cross Section, and Mass Attenuation Coefficients in the'."\n"
  .'#UD   Vicinity of Absorption Edges in the Soft X-Ray (Z = 30-36, Z = 60-89,'."\n"
  .'#UD   E = 0.1 keV-10 keV), Addressing Convergence Issues of Earlier Work."'."\n"
  .'#UD   '."\n"
  .'#UD   [2] C.T. Chantler, J. Phys. Chem. Ref. Data 24, 71 (1995). '."\n"
  .'#UD   "Theoretical Form Factor, Attenuation and Scattering Tabulation for '."\n"
  .'#UD   Z=1-92 from E=1-10 eV to E=0.4-1.0 MeV."'."\n";

  for ($z=$zmin; $z<=$zmax; $z++) {
     $j = $z-1;
     if ($elements[$j]{'z'} != $z) {
        die '*** '.$script.': Unexpected z='.$elements[$j]{'z'}.' in the elements array, expecting '.$z."\n";
     }
     $csv = $path.'mac['.$z.'].csv';
     if (! -e $csv) {
        die '*** '.$script.': No file ['.$csv.'] found'."\n";
     }
     $code = $elements[$j]{'code'};
     print DAT '#S '.$z.'  '.$code."\n"
              .'#UB '.$nodata."\n"
              .'#UF1ADD 0.0'."\n"
              .'#L PhotonEnergy[eV]  f1  f2'."\n"
              .'#UO  '.$code.', Z='.$z.', (Energy (eV),f1,f2)'."\n";
     print STDOUT $script.': Opening '.$csv."\n";
     $status = &read_entire_file($csv,\@lines);
     if ($status) {
         die '*** '.$script.': Error reading data file ['.$csv.'].'."\n";
     }
     foreach $str (@lines) {
        if ($str =~ /^E,f1,f2,/) {next;}
        if (length($str) == 0) {next;}
        @words = split(/,/,$str);
        $nwords = @words;
        if ($nwords != 8) {
           die '*** '.$script.': incorrect record ['.$str.'] in data file ['.$csv.'] -- nWords='.$nwords."\n";
        }
        if (! &is_float($words[0]) || ! &is_float($words[1]) || ! &is_float($words[2])) {
           die '*** '.$script.': incorrect record ['.$str.'] in data file ['.$csv.'] -- non-float data'."\n";
        }
        if ($words[2] <= 0) {$words[2] = $nodata;}
        printf DAT ' %12.5f   %-13.6g   %-13.6g'."\n", 1000*$words[0],$words[1],$words[2];    #from KeV to eV
     }
  }
  close(DAT);
  print STDOUT "\n".$script.': All done.'."\n";
  exit 0;
