#!/usr/bin/perl -w
  use warnings 'all';
# use warnings::unused;
  use strict;
  use File::Copy 'move';
  use File::Path 'rmtree';
  use File::Basename 'dirname';
  use Cwd 'abs_path';

  select STDERR; $|=1;				# Set unbuffered output
  select STDOUT; $|=1;				# Set unbuffered output

  my ($dest_x0h, $dirs_prefix, $dirs_postfix, $codeblocks, @dirs);
  my (@dirs_x0h, $ndirs, $i, @exefiles, $exefile, $script);

  print STDOUT 'This script installs X-server executables for linux'."\n";

  $codeblocks = 'CodeBlocks';
  if (defined $ENV{X0H}) {
     $dest_x0h = $ENV{X0H};
     $script = abs_path($0);
     $dirs_prefix = dirname($script);
     @dirs = split(/\n/,`find $dirs_prefix -name xserver.workspace`);
     if ($dirs[0] !~ $codeblocks) {die '*** Cannot find '.$codeblocks.' folder under '.$dirs_prefix."\n";}
     $dirs[0] =~ s/${codeblocks}.*$/$codeblocks/;
     $dirs_prefix = $dirs[0];
     print STDOUT 'dir="'.$dirs_prefix.'"'."\n";
  } else {
     die '*** Please define X0H environment pointing to the x0h folder.'."\n";
  }
  print STDOUT '       X0H directory: '.$dest_x0h."\n";
  print STDOUT 'Codeblocks directory: '.$dirs_prefix."\n";
  if (! -e $dest_x0h)    {die '*** X0H folder "'.$dest_x0h.'" does not exist.';}
  if (! -e $dirs_prefix) {die '*** Codeblocks folder "'.$dirs_prefix.'" does not exist.';}

  @dirs_x0h = qw(
brl_build
brls12
cofu_dec
gid_slm7
grd2dat
grd2gnu
grd2grd
mag_sl98
mag_sl99
ter_slm5
trds_99
x0h_web
x0p_web
);

  print STDOUT 'cd '.$dest_x0h."\n";
  chdir $dest_x0h;
  $ndirs = @dirs_x0h;
  for ($i=0; $i<$ndirs; $i++) {
     $dirs_x0h[$i] = $dirs_prefix.'/'.$dirs_x0h[$i].'/bin/';
     if    (-e $dirs_x0h[$i].'Release') {$dirs_postfix = 'Release';}
     elsif (-e $dirs_x0h[$i].'Debug')   {$dirs_postfix = 'Debug';}
     else {
        print STDOUT '!!! Neither '.$dirs_x0h[$i].'Release nor '.$dirs_x0h[$i].'Debug found!'."\n";
        undef $dirs_postfix;
     }
     if (defined $dirs_postfix) {
        $dirs_x0h[$i] .= $dirs_postfix;
        @exefiles = glob($dirs_x0h[$i].'/*');
        foreach $exefile (@exefiles) {
           print STDOUT 'Moving '.$exefile.' to '.$dest_x0h.'/'."\n";
           move($exefile, $dest_x0h.'/');
        }
        rmtree([$dirs_x0h[$i]]);
     } else {
        print STDOUT '!!! No compiled executables in '.$dirs_x0h[$i]."\n";
     }
  }

  print STDOUT "\n\n".'Done installing X-Server executables.'."\n";
