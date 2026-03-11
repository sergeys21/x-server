#!/usr/bin/perl -w
#--------------------------------------------------------------
# This script reads data from the DABAX file f1f2_Chantler.dat
# (https://ftp.esrf.fr/pub/scisoft/xop2.3/DabaxFiles/f1f2_Chantler.dat)
# which contains Chantler's X-ray scattering data and converts
# it into a form acceptable by X0H and similar to the Henke 
# tables of f1, f2.
#
# The Chantler data are  based on these two publications:
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
# Version-1:    2026/03/10
#--------------------------------------------------------------
  use strict;
# use warnings::unused;
  use File::Basename 'dirname';
  use Cwd 'abs_path';

  my ($X0H, $script, $status);
  my ($path, $inp, $out, $l, $nodata);
  my (@lines, $nlines, $line, @words, $nwords, $w);

  if (defined $ENV{X0H})	{$X0H=$ENV{X0H};}
  else				{$X0H=dirname(abs_path($0));}
  $X0H =~ s/\\/\//g;
  if    (-e $X0H.'/support_subs.pl') {do $X0H.'/support_subs.pl';}
  else  {die $0.': cannot find support_subs.pl for X0H='.$X0H."\n";}

  select STDOUT; $|=1;							#set unbuffered output
  select STDOUT; $|=1;							#set unbuffered output

  $script = $0 =~ s/^.*[\/\\]//gr;					#remove everything before last "/" or "\"

  $path = './';								#where the input and output files should be
  $inp  = 'f1f2_Chantler.dat';						#input file
  $out  = 'chantler.dat';						#output file
  $nodata = -.999900;

  $status = &read_entire_file($path.$inp,\@lines);
  if ($status) {die '*** '.$script.': failed reading '.$path.$inp."\n\n";}
  $nlines = @lines;

### Open output file for writing (may overwrite):
  open (DAT,'>',$path.$out) || die '*** '.$script.': cannot open '.$path.$out."\n";
  select DAT; $|=1;							#set unbuffered output
  print DAT '#C This file has been converted from the DABAX file '.$inp."\n"
           .'#C into the X-ray Server format using '.$script.' on '.localtime()."\n";
  for ($l=0; $l<$nlines; $l++) {
     $line = $lines[$l];
     if ($line =~ /^\#/) {
        if    ($line =~ /^\#L PhotonEnergy/) {print DAT '#L PhotonEnergy[eV]  f1  f2'."\n";}
        elsif ($line =~ /^\#UO\s+E\s+/)      {print DAT '#UO   E             f1          f2'."\n";}
        elsif ($line =~ /^\#UO\s+keV\s+/)    {print DAT '#UO   eV            e/atom      e/atom'."\n";}
        else                                 {print DAT $line."\n";}
     } else {
        $line =~ s/^\s+//;
        @words = split(/\s+/,$line);
        $nwords = @words;
        if ($nwords != 7) {die '*** '.$script.': unexpected number of words='.$nwords.' in line '.($l+1).': ['.$line.']'."\n";}
        foreach $w (@words) {
           if (! &is_float($w)) {die '*** '.$script.': unexpected word=['.$w.'] in line '.($l+1).': ['.$line.']'."\n";}
        }
        if ($words[2] <= 0) {$words[2] = $nodata;}
        printf DAT ' %12.5f   %-13.6g   %-13.6g'."\n", 1000*$words[0],abs($words[1]),$words[2];	#from KeV to eV
     }
  }
  close(DAT);
  exit 0;

