#!/usr/bin/perl -w
#--------------------------------------------------------------
# These are common subs for  most of other scripts
#
#		Author: Sergey Stepanov
#
# Version-1:  2026/03/03 derived from common_subs.pl
#--------------------------------------------------------------
# List of functions:
# ^^^^^^^^^^^^^^^^^
#  1.                        $sub_name = &sub_name([$sub_index);					#1=current sub, 2=above sub, 3=above-above...
#  2.                          $status = &try_load_module($module,[\@functions,$silent]);
#  3.                          $status = &try_load_lib($libdir,[$silent]);
#  4.                                    &try_do($file);						#version of "do $file" with error checking (Linux-specific)
#  5.                                    &try_do2($file);						#version of "do $file" with error checking (little less strict than above, but not Linux-specific)
#  6.                         $funcref = &function_exists($funcname);					#after call: if (defined $funcref) {$funcref->(@args);}
#  7.                        $continue = &ask_to_continue([$msg,$prompt,'noexit',$Yeslbl,$Nolbl]);	#if 'noexit' specified, the function returns 1=continue,0=exit
#  8.                         $consent = &ask_to_overwrite($file);
#  9.                                    &press_key_continue([$msg,$prompt]);
# 10.                                    &press_key_exit([$msg,$exitcode,'info']);
# 11.                                    &zenity_info ($message,['nowait']);
# 12.                                    &zenity_error ($message,[$stderr,'nowait']);
# 13.                                    &zenity_error_exit ($message,[$stderr,'nowait']);
# 14.                                    &zenity_proc (\$message,$type,[$stderr,$nowait]);		#$type=info|error
# 15.                           $yesno = &zenity_question ($message,[$Yeslbl,$Nolbl]);			#returns: 1=OK/Yes, 0=Cancel/No
# 16.                     $choiceIndex = &zenity_list($message,\@choices,[$default]);                    #returns index of choice or -1 if cancel or none selected
# 17.                                    &zenity_bring_on_top($title);					#for Linux only, needs wmctrl
# 18.                                    &disable_interrupts([$message]);
# 19.                         $radians = &asinSS($val);							#uses asin from PDL or Math::Trig, if available
# 20.                         $max_val = &maxval(@array);						#replaces List::Util qw(min max);
# 21.                         $min_val = &minval(@array);						#replaces List::Util qw(min max);
# 22.                         $max_val = &maxvalRef(\@array);						#replaces List::Util qw(min max);
# 23.                         $min_val = &minvalRef(\@array);						#replaces List::Util qw(min max);
# 24.                           $i_min = &findIndexMinRef(\@array);
# 25.                           $i_max = &findIndexMaxRef(\@array);
# 26.                           $y_max = &findMaxRef(\@array);
# 27.                  ($y_min,$y_max) = &findMinMaxRef(\@array,[$npts]);
# 28.                ($rmsVal,$rmsDev) = &findRMS(\@array);
# 29.                 ($y_max,$pt_max,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_peak]) = &findPeak($npts,@array);
# 30.                 ($y_max,$pt_max,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_peak]) = &findNoisyPeak($npts,@array);
# 31. ($y_max,$y_min,$pt_max,$pt_min,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_peak])  = &findPeakRef($npts,\@array,[$bkgmode,'smooth|noisy']);
# 32. ($y_max,$y_min,$pt_max,$pt_min,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_valey]) = &findValleyRef($npts,\@array,[$bkgmode,'smooth|noisy']);
# 33.             ($x_centroid,[$sum]) = &findCentroidRef($npts,\@arrx,\@arry,[$discriminator]);
# 34.             ($x_centroid,[$sum]) = &findCentroidFixRef($npts,$x0,$dx,\@arry,[$discriminator]);
# 35.                                    &smoothArray(\@array);
# 36.  ($pt_min,$y_min,$pt_max,$y_max) = &arrayStatsRef(\@array);
# 37.  ($pt_min,$y_min,$pt_max,$y_max) = &arrayStats($npts,@array);
# 38.                     ($name,$ext) = &analyzeFilename($file);
# 39.                      $scriptName = &makeScriptName([1|'short'|'noext']);					#optionally remove .pl
# 40.                       $age_hours = &getFileAgeHours($file);						#returns -1 if file does not exist.
# 41.                           $yesno = &isMounted($mountPoint);
# 42.                 ($status,[$msg]) = &read_data_file ($file,\@xdata,\@ydata,[$xcol,$ycol]); 		#if xcol=0, fills @xdata with pt number
# 43.                           $epoch = &timestamp2epoch($timestamp);						#stamp is YYYY-MM-DD HH:MM[:SS] or YYYY/MM/DD HH:MM[:SS]
# 44.                       $timestamp = &epoch2timestamp([$epoch,{/-},'nosec']);				#stamp is YYYY-MM-DD HH:MM:SS or YYYY/MM/DD HH:MM:SS; if no $epoch, reports current time
# 45.                          $sorted = &sort_by_x('incr|decr',\@x,@y1,[\@y2,\@y3,...]);			#1=sorted,0=did_not_need,-1=error
# 46.                            $npts = &remove_duplicates($filter,$step,\@x,@y1,[\@y2,\@y3,...]);		#the arrays must be sorted over x and 'quasi-regular'
# 47.                            $path = &Which($prg);								#returns path to $prg or undef if not found. On win32 tries .exe,.bat,.com
# 48.             ($dirname,$filename) = &split_filename_dirname($filespec);
# 49.                         $radians = &degrees2rad($degrees);
# 50.                         $degrees = &rad2degrees($radians);
# 51.                           $index = &find_first_index($element,\@array);					#returns -1 if nothing found
# 52.                           @index = &find_all_index($element,\@array);					#returns empty list if nothing found
# 53.                           $index = &find_closest_index($val,\@array);					#finds index of element closest to given value
# 54.                          $status = &execute_external_script(\@cmdline);					#returns error if cannot execute
# 55.                                    &abort_external_script(;$pid);						#uses global $EXTERNAQL_PID if $pid not passed
# 56.                          $status = &zenity_progress($time,$header);					#returns 1 if process was interrupted
# 57.                            $bool = &is_float($str);							#returns 1/0 (True/False if the value is any float incl. exponents)
# 58.                            $bool = &is_decimal($str);							#returns 1/0 (True/False if the value is like xxx.yyy)
# 59.                            $bool = &is_int($str);								#returns 1/0 (True/False)
# 60.                            $bool = &is_mounted($dir,[$refdir]);						#returns 1/0/-1 (True/False/NofFound)
# 61.                            $bool = &is_mountpoint($dir);							#returns 1/0/-1 (True/False/NofFound)
# 62.                       $bitstring = &dec2bin($intword,['r[everse]','t[rim]'|$nbits]);			#convert integer into bits string; specify 'reverse' for bits parsing
# 63.                         $intword = &bin2dec($bitstring);							#convert bits string into integer
# 64.                          $status = &read_cmdline_flags($prg,\@flags);					#read cmdline flags of a program. The program must support the "-h" option
# 65.                          $xcoord = &pt2coord($xpt,\@xarray);                                               #convert non-integer pt into coordinate
# 66.                            $host = &gethost(['short');
# 67.                                    &sort_hash_array(\@hasharray,$key);
# 68.                          $status = &read_entire_file($file,\@lines,[$lock]);				#returns 1=error, 0=OK
# 69.                             $msg = &call_trace([$LUN,'noprint|no|0]);					#function calling trace (by default prints to STDOUT)
# 70.                          $status = &get_login_environment();
# 71.                          $status = &ziplist($zipfile,\%hash);						#Fills %hash as key=filename, value=length
# 72.                          $status = &zipextract($zipfile,$file,\@lines);					#Fills @lines with file content
# 73.                          $status = &import_elements(\@elements);						#fills hash array of the Periodic Table elements (name, code, z, a)
# 74.                       $userAgent = &getLatestUserAgent();

#--------------------------------------------------------------
# use ENV;							#load the ENV
# use Time::Local;						#needed by timestamp2epoch
# use Carp qw(shortmess);					#may be needed for stack trace in report_exit() and report_only()
# use Tk;							#may be needed by zenity_progress (MS Windows only)
# use Tk::ProgressBar;						#may be needed by zenity_progress (MS Windows only)
# use POSIX ":sys_wait_h";					#may be used by execute_external_script and abort_external_script

  use warnings 'all';
  use warnings::unused;
  use strict;
# no strict "vars";
# no strict "subs";
# no strict "refs";

#=================================================================

sub sub_name(;$) {			#optional level: 1=current caller, 2=above, 3=above-above...
  my ($j, $script);
# The caller($i) function contains the calls stack where for each function $j we can decode:
#      0         1      2        3         4         5          6          7        8       9        10
# ($package,$filename,$line,$subroutine,$hasargs,$wantarray,$evaltext,$is_require,$hints,$bitmask,$hinthash) = caller($i);
# Note that all components 0 to 10 may be defined.
# We are interested in the component $subroutine(3) for the caller of this sub_name() function(1).
# In other words, caller(0) provides data for sub_name() itself
# while caller(1) provides data for the function calling sub_name().
  $j = shift(@_); 
  if (!defined $j || $j < 1) {$j = 1;}
  if (defined caller($j)) {
     $script = (caller($j))[3];
     $script = (split(/::/,$script))[-1];	#remove the prefix, e.g.:   main::setEpicsRetry => setEpicsRetry
  } else {
     $script = &makeScriptName();
  }
# print STDOUT $script."\n";
  return $script;
}

#=================================================================

#sub try_load_module($module,[\@functions,$silent]) {
sub try_load_module($;\@) {
  my ($script, $module, $funcref, $silent, $i, $e);

  $module  = shift(@_);
  $funcref = shift(@_);
  $silent  = shift(@_);		if (!defined $silent) {$silent = 0;}

  select STDOUT; $|=1;					#set unbuffered output

  if (defined $funcref && @$funcref > 0) {
     for ($i=0; $i<@$funcref; $i++) {
        if ($i == 0) {$module .= ' qw/'.@$funcref[$i];}
        else         {$module .= ' '.@$funcref[$i];}
     }
     $module .= '/';
  }
# print STDOUT 'use '.$module."\n";

  eval('use '.$module);
  ### $@ is set if the string to be eval-ed did not compile or if Perl code
  ### executed during evaluation die()d. In these cases the value of $@ is
  ### the compile error, or the argument to die.
  if ($@) {
     if (! $silent) {
        $e = $@;
        $script = &makeScriptName();
        print STDOUT "\n".'*** '.$script.': module "'.$module.'" failed to load.'."\n".$e."\n\n";
     }
     return 1;
  } else {
     return 0;
  }
}

#=================================================================

#sub try_load_lib($libdir,[$silent]) {
sub try_load_lib($;$) {
  my ($script, $libdir, $silent, $e);

  $libdir = shift(@_);
  $silent = shift(@_);		if (!defined $silent) {$silent = 0;}

  select STDOUT; $|=1;					#set unbuffered output

# print STDOUT 'use lib "'.$libdir.'"'."\n";

  eval('use lib "'.$libdir.'"');
  ### $@ is set if the string to be eval-ed did not compile or if Perl code
  ### executed during evaluation die()d. In these cases the value of $@ is
  ### the compile error, or the argument to die.
  if ($@) {
     if (! $silent) {
        $e = $@;
        $script = &makeScriptName();
        print STDOUT "\n".'*** '.$script.': failed to load library from dir='.$libdir."\n".$e."\n\n";
        return 1;
     }
  } else {
     return 0;
  }
}

#=================================================================

#sub try_do($file) {                                                    # Linux-specific!
sub try_do($) {
  my ($file, $st);
  $file = shift(@_);
  if (-e $file) {
     $st = eval `cat $file`;                                            #"cat" is a Linux command.
     if ($@)              {die 'Failed parsing '.$file.': '.$@;}        #'$@' is error message
     elsif ($!)           {die 'Could not eval '.$file.': '.$!;}        #'$!' is error message
#    elsif (!defined $st) {die 'Failed to eval '.$file.': '.$!;}        #this does not work when $file contains subs: then $st is undef
#    elsif (! $st)        {die 'Failed to run '.$file;}                 #this does not work when $file contains subs: then $st is undef
  } else                  {die 'Did not find '.$file;}
}

#=================================================================

#sub try_do2($file) {                                                   #this is less strict than try_do(), but not Linux-specific
sub try_do2($) {
  my ($file, $st);
  $file = shift(@_);
  if (-e $file) {
     $st = do $file;
     if ($@)              {die 'Failed parsing '.$file.': '.$@;}        #'$@' is error message
     elsif ($!)           {die 'Could not do '.$file.': '.$!;}          #'$!' is error message
#    elsif (!defined $st) {die 'Failed to do '.$file.': '.$!;}          #this does not work when $file contains subs: then $st is undef
#    elsif (! $st)        {die 'Failed to run '.$file;}                 #this does not work when $file contains subs: then $st is undef
  } else                  {die 'Did not find '.$file;}
}

#=================================================================

sub function_exists ($) {
  ### This function determines if function exists.
  ### An example use:
  ### $funcref = &function_exists($funcname);
  ### if (defined $funcref) {$funcref->(@args);}	#the way to call function by name
  ### else                  {die $funcname.' is undefined';}
  my $funcname = shift;
  if (defined &{$funcname}) {return \&{$funcname};}
  else                      {return;}
}

#=================================================================

#sub ask_to_continue ([$msg,$prompt,'noexit',$Yeslbl,$Nolbl]) {
sub ask_to_continue (;$$$$$) {

  our ($NONINTERACTIVE, $USE_ZENITY_FOR_EXIT);
  my ($msg, $prompt, $txt, $noexit, $Yeslbl, $Nolbl);
  $msg    = shift(@_);
  $prompt = shift(@_);		if (!defined $prompt) {$prompt = 'Continue?';}
  $noexit = shift(@_);		if (!defined $noexit) {$noexit = 'exit';}
  $Yeslbl = shift(@_);
  $Nolbl  = shift(@_);
 
  if (!defined $NONINTERACTIVE)      {$NONINTERACTIVE=0;}
  if (!defined $USE_ZENITY_FOR_EXIT) {$USE_ZENITY_FOR_EXIT=0;}

  if ($USE_ZENITY_FOR_EXIT && $NONINTERACTIVE == 0) {
     if	(defined $msg) {$txt = $msg."\n\n".$prompt;}
     else              {$txt = $prompt;}
     $_ = &zenity_question ($txt,$Yeslbl,$Nolbl);			#returns 1 or 0
  } else {
     if (defined $msg)    {print STDOUT "\n".$msg."\n\n";}
     if (defined $prompt) {print STDOUT $prompt;}
     else                 {print STDOUT 'Continue [n]?_';}
     if ($NONINTERACTIVE == 0) {$_=uc <STDIN>;  chomp $_;}
     else                      {$_='Y'; print STDOUT 'Y';}
  }
  if ($_ ne 'Y' && $_ ne '1') {
     print STDOUT ' ....... Exiting..'."\n\n";
     if ($noexit ne 'noexit') {
        &clearAllCaMonitors();
        exit;
     } else {
        return 0;
     }
  } else {
     print STDOUT ' ....... Continuing..'."\n\n";
     return 1;
  }
}

#=================================================================

#sub ask_to_overwrite ($file) {
sub ask_to_overwrite ($) {

  our ($NONINTERACTIVE, $USE_ZENITY_FOR_EXIT);
  my ($file, $yesno);
  if (!defined $NONINTERACTIVE)      {$NONINTERACTIVE      = 0;}
  if (!defined $USE_ZENITY_FOR_EXIT) {$USE_ZENITY_FOR_EXIT = 0;}
  $file = shift(@_);
  if (-e $file) {			#file exists
    if (-f $file) {			#file is file
      if (-w $file) {			#file is writable
         if ($USE_ZENITY_FOR_EXIT) {
            $yesno = &zenity_question('File "'.$file.'" exists. Overwrite?');
            return $yesno;
         } else {
            print STDOUT 'File "'.$file.'" exists. Overwrite [n]?_';
            if ($NONINTERACTIVE == 0) {$_=uc <STDIN>;  chomp $_;}
            else                      {$_='Y'; print STDOUT 'Y';}
            print STDOUT "\n";
            if ($_ eq 'Y') {return 1;}
            else           {return 0;}
         }
      } else {
        return 0;			#file is writable, NO overwrite
      }
    } else {
      return 0;				#file is directory, NO overwrite
    }
  } else {
    return 1;				#no file, OK to create new file
  }
}

#=================================================================

#sub press_key_continue ([$msg,$prompt]) {
sub press_key_continue (;$$) {

  our $NONINTERACTIVE;
  my $msg    = shift(@_);
  my $prompt = shift(@_);
  if (!defined $NONINTERACTIVE) {$NONINTERACTIVE=0;}

  if (defined $msg && length($msg) > 0) {print STDOUT "\n".$msg."\n\n";}
  if ($NONINTERACTIVE == 0) {
     if (defined $prompt) {print STDOUT  $prompt;}
     else                 {print STDOUT  'Hit <ENTER> to continue...';}
     $_ = <STDIN>;
     print STDOUT "\n";
  }
}

#=================================================================

#sub press_key_exit ([$msg,$exitcode,'info|error',$stderr]) {
sub press_key_exit (;$$$$) {
  our ($NONINTERACTIVE, $USE_ZENITY_FOR_EXIT);
  my ($msg, $exitcode, $icontype, $stderr, $LUN);
  $msg      = shift(@_);
  $exitcode = shift(@_);
  $icontype = shift(@_);					#error/info
  $stderr   = shift(@_);					#1/0
  if (!defined $NONINTERACTIVE)      {$NONINTERACTIVE = 0;}
  if (!defined $USE_ZENITY_FOR_EXIT) {$USE_ZENITY_FOR_EXIT = 0;}
  if (!defined $exitcode)	     {$exitcode = 0;}
  if ($exitcode ne '1')		     {$exitcode = 0;}
  if (!defined $icontype || $icontype !~ /^info$/i) {$icontype = 'error';}
  if (!defined $stderr   || $stderr   != 1)         {$stderr   = 0;}

  &clearAllCaMonitors();

  if ($stderr) {$LUN = \*STDERR;}
  else         {$LUN = \*STDOUT;}

  if ($USE_ZENITY_FOR_EXIT) {
     ### Graphical mode
     if (!defined $msg || length($msg) == 0) {$msg = 'Click OK to exit...';}
     if ($icontype =~ /^info$/i) {&zenity_info ($msg);}
     else                        {&zenity_error ($msg,$stderr);}
  }  else  {
     ### Text mode
     if (defined $msg && length($msg) > 0) {print $LUN "\n".$msg."\n";}
     if (!$NONINTERACTIVE) {
        print $LUN  "\n".'Hit <ENTER> to exit...';
        $_ = <STDIN>;
     }
     print $LUN "\n";
  }
  exit $exitcode;
}

#=================================================================

#sub zenity_info ($message,['nowait']) {
sub zenity_info ($;$) {

### Open graphical window to report a notification
### zenity --info --text='a-a\nb-b' --title='Test' --no-wrap --ellipsize

  my ($msg, $nowait, $stderr);
  $msg    = shift(@_);
  $nowait = shift(@_);
  $stderr = 0;					#output to STDOUT
  &zenity_proc(\$msg,'info',$stderr,$nowait);
}

#=================================================================

#sub zenity_error ($message,[$stderr,$nowait]) {
sub zenity_error ($;$$) {

### Open graphical window to report an error
### zenity --error --text='a-a\nb-b' --title='Test' --no-wrap --ellipsize

  my ($msg, $stderr, $nowait);
  $msg    = shift(@_);
  $stderr = shift(@_);
  $nowait = shift(@_);
  &zenity_proc(\$msg,'error',$stderr,$nowait);
}

#=================================================================

#sub zenity_error_exit ($message,[$stderr,'nowait']) {
sub zenity_error_exit ($;$$) {

### Open graphical window to report an error and exits
### zenity --error --text='a-a\nb-b' --title='Test' --no-wrap --ellipsize

  our ($IGNORE_ERRORS);
  my ($msg, $stderr, $nowait);
  $msg    = shift(@_);
  $stderr = shift(@_);
  $nowait = shift(@_);
  &zenity_proc(\$msg,'error',$stderr,$nowait);
  if (!defined $IGNORE_ERRORS) {$IGNORE_ERRORS = 0;}
  if (!$IGNORE_ERRORS) {
     &clearAllCaMonitors();
     exit 1;
  }
}

#=================================================================

#sub zenity_proc (\$message,$type,[$stderr,'wait|nowait']) {
sub zenity_proc (\$$;$$) {

### Open graphical window to report an error or post an information
### ($type can be 'error' and 'info')

  my ($script, $msg_r, $msg, $type, $stderr, $wait);
  my ($LUN, $prefix, $title, $path, @cmd, $string);

  $script = &makeScriptName();

  $msg_r  = shift(@_);                                                  #pointer to message
  $type   = shift(@_);                                                  #'error' or 'info'
  $stderr = shift(@_);  if (!defined $stderr) {$stderr = 0;}		#0=stdout or 1=stderr
  $wait   = shift(@_);  if (!defined $wait)   {$wait = 'wait';}		#'nowait' or any meaning 'wait'

  $msg = $$msg_r;
  $msg =~ s/\"/\\\"/g;
  $msg =~ s/\'/`/g;
# $msg =~ s/\&/and/g;
  $msg =~ s/\&/&amp;/g;			#2024.06: use HTML notation
# $msg =~ s/</\\\</g;
# $msg =~ s/</\</g;
  $msg =~ s/</&lt;/g;			#2024.06: use HTML notation
# $msg =~ s/>/\\\>/g;
# $msg =~ s/>/\>/g;
  $msg =~ s/>/&gt;/g;			#2024.06: use HTML notation
  if ($^O !~ /Win/i) {
     $msg =~ s/\(/\\\(/g;
     $msg =~ s/\)/\\\)/g;
  }
  $title = $script.'-'.$type;
  $title =~ s/ /_/g;
  @cmd = ('--'.$type,
          '--title='."'".$title."'",
          '--text='."'".$msg."'");
  $path = '/usr/bin/zenity';
  if ($^O !~ /Win/i && -f $path) {					#Linux with zenity installed
### zenity --error --text='bla-bla-aaaaaaaaaaaaaaa ee\rbla' --title='Test' --no-wrap --ellipsize --window-icon='/usr/share/pixmaps/applet-critical.png'
     unshift @cmd, $path;						#prepend with 'path'
     push @cmd, '--no-wrap';						#this helps to reduce excessive white space below text
     push @cmd, '--ellipsize';						#this fixes the high window size with long texts
#    push @cmd, '2>/dev/null';						#redirect stderr
#    push @cmd, '&>/dev/null';						#redirect both stdout and stderr
     push @cmd, ' >/dev/null', '2>&1';					#redirect both stdout and stderr (preferred)
     if ($wait =~ /^nowait$/i) {
        push @cmd, '&';
#       unshift (@cmd,'screen','-d','-m');
        unshift (@cmd,'nohup');
     }
     $string = join(' ',@cmd);						#merge from array
     &zenity_bring_on_top($title);
     system($string);							#only string form allows no-wait
  } else {								#either no zenity or no Linux
     print STDERR '*** '.$script.': no zenity available!'."\n";
  }

  $LUN = \*STDOUT;
  if ($stderr) {$LUN = \*STDERR;}
  $prefix = '';
  if ($type =~ /^error$/) {$prefix = '*** ';}
  print $LUN "\n".$prefix.$script.':'."\n".$$msg_r."\n\n";
}

#=================================================================

#sub zenity_question ($message,[$Yeslbl,$Nolbl]) {			#return 1/0=yes/no
sub zenity_question ($;$$) {

### Open graphical window to query a choice of 2 options
### zenity --question --text='a-a\nb-b' --ok-label='OK' --cancel-label='Cancel' --title='Test' --no-wrap --ellipsize

  my ($script, $msg, $Yeslbl, $Nolbl, $title, $path, $yesno);

  $script = &makeScriptName();

  $msg    = shift(@_);	if (!defined $msg) {$msg = '?';}
  $Yeslbl = shift(@_);	if (!defined $Yeslbl) {$Yeslbl = 'OK';}
  $Nolbl  = shift(@_);	if (!defined $Nolbl)  {$Nolbl  = 'Cancel';}

  $title = $script.'-question';
  $title =~ s/ /_/g;
  $path   = '/usr/bin/zenity';
  if ($^O !~ /Win/i && -f $path) {			                #Linux with zenity installed
     ### xmessage -buttons Yes:1,No:0 'bla-bla-bla'
     ### zenity --question --no-wrap --ellipsize --text='bla-bla-bla\rbla' --title='Test' --ok-label='Continue' --cancel-label='Dismiss' --window-icon='/usr/share/pixmaps/applet-critical.png'
     &zenity_bring_on_top($title);
     $yesno = system($path.' --question'
                          .' --no-wrap'
                          .' --ellipsize'
                          .' --ok-label='    ."'" .$Yeslbl ."'"
                          .' --cancel-label='."'" .$Nolbl  ."'"
                          .' --title='       ."'" .$title  ."'"
                          .' --text='        ."'" .$msg    ."'"
                          .' 2>/dev/null');
     if ($yesno ne '0') {$yesno = 0;}		                #Cancel returns '1'
     else               {$yesno = 1;}		                #OK returns '0'
  } else {
     print STDERR '*** '.$script.': no zenity available!'."\n";
  }
  ### Replace \r by carriage return:
  $msg =~ s/\\\\n/\n/g;
  $msg =~ s/\\n/\n/g;
  $msg =~ s/\\r/\n/g;
  print STDOUT "\n".$script.":\n".$msg.'  => ';
  if ($yesno == 1) {print STDOUT $Yeslbl."\n";}
  else             {print STDOUT $Nolbl."\n";}
  return $yesno;
}

#=================================================================

#sub zenity_list($message,\@choices,[$default]) {                       #returns index of choice or -1 if cancel or none selected
sub zenity_list ($\@;$) {

### Open graphical window to request a single choice from a list of options
### zenity --list --radiolist --text='bla-bla-bla' --title='Test' --column='Select' --column='Option' TRUE 'aaa' FALSE 'bbb' FALSE 'ccc'

  my ($script, $subname, $msg, $choices_ref, $default);
  my ($str, $path, $title, $nchoices);
  my ($i, $j, $result);

  $script = &makeScriptName();
  $subname = (split(/::/,(caller(0))[3]))[-1];

  $msg         = shift(@_);
  $choices_ref = shift(@_);     $nchoices = @$choices_ref;
  $default     = shift(@_);     if (!defined $default) {$default = -1;}

  for ($i=0; $i<$nchoices; $i++) {
     for ($j=0; $j<$nchoices; $j++) {
        if ($i != $j && $$choices_ref[$i] eq $$choices_ref[$j]) {
           print STDOUT '*** '.$subname.': duplicate choice "'.$$choices_ref[$i].'"'."\n";
           return -1;
        }
     }     
  }  

  $title = $script.'-question';
  $title =~ s/ /_/g;
  $path   = '/usr/bin/zenity';
  if ($^O !~ /Win/i && -f $path) {			                #Linux with zenity installed
     ### zenity --list --radiolist --text='bla-bla-bla' --title='Test' --column='Select' --column='Option' TRUE 'aaa' FALSE 'bbb' FALSE 'ccc'
     &zenity_bring_on_top($title);
     $str = $path
          . ' --list'
          . ' --radiolist'
          . ' --text='        ."'".$msg."'"
          . ' --title='       ."'".$title."'"
          . ' --column='      ."'Select'"
          . ' --column='      ."'Option'";
     for ($i=0; $i<$nchoices; $i++) {
        if ($i == $default) {$str .= ' TRUE';}
        else                {$str .= ' FALSE';}
        $str .= " '".$$choices_ref[$i]."'";
     }
     $str .= ' 2>/dev/null';
#    print STDOUT $str,"\n";
     $result = `$str`;
     chomp($result);
#    print STDOUT $result,"\n";
     for ($i=0; $i<$nchoices; $i++) {
        $str = $$choices_ref[$i];
        chomp($str);
        if ($result eq $str) {return $i;}
     }
     return -1;
  } else {
     print STDERR '*** '.$script.': no zenity available!'."\n";
     return -1;
  }
}


#=================================================================

#sub zenity_bring_on_top ($title) {
sub zenity_bring_on_top ($) {

### Open graphical window to report an error

  my ($title, $wmctrl, $bash, $sleep);

  $title = shift;
  $wmctrl = '/usr/bin/wmctrl';
  $bash   = '/usr/bin/bash';
  $sleep  = '/bin/sleep';
  if (-e $wmctrl && -e $bash && -e $sleep) {
     system($bash.' -c "'.$sleep.' 2 && '.$wmctrl.' -r '.$title.' -b toggle,above" &');
  }
  return;
}
#=================================================================

#sub disable_interrupts ([$message]) {
sub disable_interrupts  (;$) {
### Disable all interrupts so that the interrupt
### routine is not called several times
  my ($script, $msg);
  $msg = shift(@_);
  if (!defined $msg) {
     $script = &makeScriptName();
     $msg = $script.': further interrupts are disabled.';
  }
  print STDOUT "\n".$msg."\n\n";
  $SIG{HUP}  = 'IGNORE';	#(1) hangup (see signum.h)
  $SIG{INT}  = 'IGNORE';	#(2) interrupt (see signum.h)
  $SIG{QUIT} = 'IGNORE';	#(3) quit (see signum.h)
# $SIG{KILL} = 'IGNORE';	#(9) this is what is sent by kill -9 -- unblockable
  $SIG{TERM} = 'IGNORE';	#(15) terminate -- sent by kill without paramaters
# $SIG{STOP} = 'IGNORE';	#(19) stop -- unblockable
}

#=================================================================

sub asinSS ($) {
### This sub uses standard asin from PDL or Math::Trig, if they are available
  local $a = $_[0];
  if    ($a > 1.)  {$a = 1.;}
  elsif ($a < -1.) {$a = -1.;}
  if (defined &asin) {asin($a);}
  else               {atan2($a,sqrt(1-$a*$a));}
}

#=================================================================

#sub maxval (@array) {
sub maxval (@) {
### This is a home-made replacement for List::Util qw(min max);
  my ($script, $npts, $max_, $i);

  $script = (split(/::/,(caller(0))[3]))[-1];

  $npts = @_;
  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0;
  }
  $max_ = $_[0];
  foreach $i (1 .. $#_) {
     if (defined $_[$i]) {$max_ = $_[$i] if $_[$i] > $max_;}
  }
  return $max_;
}

#=================================================================

#sub minval (@array) {
sub minval (@) {
### This is a home-made replacement for List::Util qw(min max);
  my ($script, $npts, $min_, $i);

  $script = (split(/::/,(caller(0))[3]))[-1];

  $npts = @_;
  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0;
  }
  $min_ = $_[0];
  foreach $i (1 .. $#_) {
     if (defined $_[$i]) {$min_ = $_[$i] if $_[$i] < $min_;}
  }
  return $min_;
}

#=================================================================

#sub maxvalRef (\@array) {
sub maxvalRef (\@) {
### This is a home-made replacement for List::Util qw(min max);
  my ($script, $npts, $arrayRef, $element, $max_);
  $script   = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
  $npts     = @$arrayRef;
  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0;
  }
  $max_ = @$arrayRef[0];
  foreach $element (@$arrayRef) {
     if (defined $element) {$max_ = $element if $element > $max_;}
  }
  return $max_;
}

#=================================================================

#sub minvalRef (\@array) {
sub minvalRef (\@) {
### This is a home-made replacement for List::Util qw(min max);
  my ($script, $npts, $arrayRef, $element, $min_);
  $script   = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
  $npts     = @$arrayRef;
  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0;
  }
  $min_ = @$arrayRef[0];
  foreach $element (@$arrayRef) {
     if (defined $element) {$min_ = $element if $element < $min_;}
  }
  return $min_;
}

#=================================================================

#sub findIndexMinRef (\@array) {
sub findIndexMinRef (\@) {
### This script finds index of the array element with min value
  my ($script, $npts, $arrayRef, $a_, $pt, $array_min, $index_min);
  $script   = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
  $npts     = @$arrayRef;
  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0;
  }
  $array_min = @$arrayRef[0];
  $index_min = 0;
  for ($pt=1; $pt < $npts; $pt++) {
     $a_ = @$arrayRef[$pt];
     if ($a_ < $array_min) {$array_min=$a_; $index_min=$pt;}
  }
  return $index_min;
}

#=================================================================

#sub findIndexMaxRef (\@array) {
sub findIndexMaxRef (\@) {
### This script finds index of the array element with max value
  my ($script, $npts, $arrayRef, $a_, $pt, $array_max, $index_max);
  $script   = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
  $npts     = @$arrayRef;
  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0;
  }
  $array_max = @$arrayRef[0];
  $index_max = 0;
  for ($pt=1; $pt < $npts; $pt++) {
     $a_ = @$arrayRef[$pt];
     if ($a_ > $array_max) {$array_max=$a_; $index_max=$pt;}
  }
  return $index_max;
}

#=================================================================

#sub findMaxRef (\@array) {
sub findMaxRef (\@) {
### This subroutine returns ($array_max)
  my ($script, $arrayRef, $array_max);
  $script   = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
### Find maximum:
  $array_max = &maxvalRef($arrayRef);
  return $array_max;
}

#=================================================================

#sub findMinMaxRef (\@array,[$npts]) {
sub findMinMaxRef (\@;$) {
### This subroutine returns ($array_min,$array_max)
  my ($script, $npts, $nmax, $arrayRef, $array_max, $array_min, $pt, $a_);
  $script   = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
  $npts     = shift(@_);
  $nmax     = @$arrayRef;
  if ($nmax <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0;
  }
  if (!defined $npts || $npts > $nmax) {$npts = $nmax;}

### Find maximum/minimum:
  $array_max = @$arrayRef[0];		#array maximum value
  $array_min = @$arrayRef[0];		#array minimum value
  for ($pt=1; $pt < $npts; $pt++) {
     $a_ = @$arrayRef[$pt];
     if ($a_ > $array_max) {$array_max = $a_;}
     if ($a_ < $array_min) {$array_min = $a_;}
  }
  return $array_min,$array_max;
}

#=================================================================

# sub findRMS(\@array) {		#returns ($rmsVal,$rmsDev)

sub findRMS (\@) {
### This subroutine returns ($rmsVal,$rmsDev)
  my ($script, $npts, $arrayRef, $rmsVal, $rmsDev, $i, $a_);
  $script   = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
  $npts     = @$arrayRef;
  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0;
  }

### Find mean RMS:
  $rmsVal = 0;
  for ($i=0; $i<$npts; $i++) {
     $a_ = @$arrayRef[$i];
     $rmsVal += $a_*$a_;
  }
  $rmsVal = sqrt($rmsVal/$npts);

### Find RMS deviation:
  $rmsDev = 0;
  for ($i=0; $i<$npts; $i++) {
     $a_ = @$arrayRef[$i] - $rmsVal;
     $rmsDev += $a_*$a_;
  }
  $rmsDev = sqrt($rmsDev/$npts);

  return $rmsVal,$rmsDev;
}

#=================================================================

#sub findPeak ($npts, @array) {
sub findPeak ($@) {
### This subroutine returns:
### ($array_max,$pt_max,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_peak])
### where $x_left,$x_right,$x_peak are respective interpolated value (i.e.
### they are a higher precision than $pt_left,$pt_right,$pt_FWHM.
  my ($script, $npts, @array, $narg, $array_max, $pt_max, $pt_FWHM);
  my ($pt_left, $pt_right, $x_left, $x_right, $x_peak);
  my ($half_intensity, $pt, $dy);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $narg   = @_;
  $npts   = shift(@_);

  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0,-1;
  }
  elsif ($npts > $narg-1) {
     printf STDERR "\n".'!!! '.$script.': received %d array elements instead of %d declared!'."\n\n",($narg-1),$npts;
     $npts = $narg-1;
  }

### Find maximum
  $pt_max    = -1;				# array index at maximum value
  $pt_FWHM   = -1;				# peak width in pts (if exists)
  $array_max = -1.E37;				# array maximum value
  for ($pt=0; $pt < $npts; $pt++) {
     $array[$pt] = shift(@_);
     if ($array[$pt] > $array_max) {
              $array_max = $array[$pt];
              $pt_max    = $pt;
     }
  }
  $half_intensity = 0.5*$array_max;
### Climb down to the left of the peak:
  $pt_left = $pt_max;
  while (($array[$pt_left] > $half_intensity) && ($pt_left > 0)) {
     $pt_left--;
  }
### Climb down to the right of the peak:
  $pt_right = $pt_max;
  while (($array[$pt_right] > $half_intensity) && ($pt_right < $npts-1)) {
     $pt_right++;
  }
### If peak found:
  if (($pt_right < $npts-1) && ($pt_left > 0)) {
     $pt_FWHM = $pt_right - $pt_left + 1;
     $x_left = $pt_left;
     $dy = $array[$pt_left+1] - $array[$pt_left];
     if ($dy != 0.) {$x_left += ($half_intensity - $array[$pt_left]) / $dy;}
     $x_right = $pt_right;
     $dy = $array[$pt_right-1] - $array[$pt_right];
     if ($dy != 0.) {$x_right -= ($half_intensity - $array[$pt_right]) / $dy;}
     $x_peak = 0.5 * ($x_left + $x_right);
  } else {
     $pt_FWHM = -1;
     $x_left  = $pt_left;
     $x_right = $pt_right;
     $x_peak  = $pt_max;
  }
  return $array_max,$pt_max,$pt_FWHM,$pt_left,$pt_right,$x_left,$x_right,$x_peak;
}

#=================================================================

#sub findNoisyPeak ($npts, @array) {
sub findNoisyPeak ($@) {
### This subroutine returns:
### ($array_max,$pt_max,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_peak])
### where $x_left,$x_right,$x_peak are respective interpolated value (i.e.
### they are a higher precision than $pt_left,$pt_right,$pt_FWHM.
### +---------------------------------------+
### |ATTENTION: It does not smooth anything!|
### +---------------------------------------+
### All the difference between "normal" and "noisy" peak
### finding is that in "normal" case we climb down the peak
### to find FWHM, while in "noisy" we climb up from the tails.
###
  my ($script, $npts, @array, $narg, $array_max, $pt_max, $pt_FWHM);
  my ($pt_left, $pt_right, $x_left, $x_right, $x_peak);
  my ($half_intensity, $pt, $dy);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $narg   = @_;
  $npts   = shift(@_);

  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0,-1,0,0,0,0,0;
  }
  elsif ($npts > $narg-1) {
     printf STDERR "\n".'!!! '.$script.': received %d array elements instead of %d declared!'."\n\n",($narg-1),$npts;
     $npts = $narg-1;
  }

  $array_max = -1.E37;			#array maximum value
  $pt_max    = -1;			#array index at maximum value
  $pt_FWHM   = -1;			#peak width in pts (if exists)
### Find maximum:
  for ($pt=0; $pt < $npts; $pt++) {
     $array[$pt] = shift(@_);
     if ($array[$pt] > $array_max) {
        $array_max = $array[$pt];
        $pt_max    = $pt;
     }
  }
  $half_intensity = 0.5*$array_max;
### Climb up on the peak from curve start:
  $pt_left = 0;
  while (($array[$pt_left] < $half_intensity) && ($pt_left < $pt_max)) {
     $pt_left++;
  }
### Climb up on the peak from curve end:
  $pt_right = $npts-1;
  while (($array[$pt_right] < $half_intensity) && ($pt_right > $pt_max)) {
     $pt_right--;
  }
### If peak found:
  if (($pt_right < $npts-1) && ($pt_left > 0)) {
     $pt_FWHM = $pt_right - $pt_left + 1;
     $x_left = $pt_left;
     $dy = $array[$pt_left] - $array[$pt_left-1];
     if ($dy != 0.) {$x_left += ($half_intensity - $array[$pt_left]) / $dy;}
     $x_right = $pt_right;
     $dy = $array[$pt_right] - $array[$pt_right+1];
     if ($dy != 0.) {$x_right -= ($half_intensity - $array[$pt_right]) / $dy;}
     $x_peak = 0.5 * ($x_left + $x_right);
  } else {
     $pt_FWHM = -1;
     $x_left  = $pt_left;
     $x_right = $pt_right;
     $x_peak  = $pt_max;
  }
  return $array_max,$pt_max,$pt_FWHM,$pt_left,$pt_right,$x_left,$x_right,$x_peak;
}

#=================================================================

#sub findPeakRef ($npts, \@array, [$bkgmode,'smooth|noisy']) {
sub findPeakRef ($\@;$$) {
### This subroutine returns:
### ($array_max,$array_min,$pt_max,$pt_min,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_peak])
### where $x_left,$x_right,$x_peak are respective interpolated value (i.e.
### they are a higher precision than $pt_left,$pt_right,$pt_FWHM).
### +---------------------------------------+
### |ATTENTION: It does not smooth anything!|
### +---------------------------------------+
### All the difference between "normal" and "noisy" peak
### finding is that in "normal" case we climb down the peak
### to find FWHM, while in "noisy" we climb up from the tails.
###
  my ($script, $npts, $arrayRef, $bkgmode, $searchmode);
  my ($array_max, $array_min, $pt_max, $pt_min);
  my ($pt_FWHM, $array_span, $a_, $half_intensity);
  my ($pt_left, $pt_right, $x_left, $x_right, $x_peak, $pt, $dy);

  $script     = (split(/::/,(caller(0))[3]))[-1];
  $npts       = shift(@_);
  $arrayRef   = shift(@_);
  $bkgmode    = shift(@_);               		#zero(default) or auto - differs from findValley
  $searchmode = shift(@_);               		#smooth(default) or noisy
  if (!defined $bkgmode)    {$bkgmode    = 'zero';}
  if (!defined $searchmode) {$searchmode = 'smooth';}

  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0,-1,-1,0,0,0,0,0;
  }

  $array_max = -1.E37;				        #array maximum value
  $array_min = +1.E37;				        #array minimum value
  $pt_max    = -1;					#array index at maximum value
  $pt_min    = -1;					#array index at minimum value
  $pt_FWHM   = -1;					#peak width in pts (if exists)

### Find maximum/minimum:
  for ($pt=0; $pt < $npts; $pt++) {
     $a_ = @$arrayRef[$pt];
     if ($a_ > $array_max) {$array_max = $a_; $pt_max = $pt;}
     if ($a_ < $array_min) {$array_min = $a_; $pt_min = $pt;}
  }
  if ($array_max < 0) {
     print STDOUT '*** '.$script.': all intensity is negative (max='.$array_max.')'."\n";
     print STDOUT '*** '.$script.': possibly incorrect background data.'."\n";
  }
  elsif ($array_min < 0) {
     print STDOUT '!!! '.$script.': some intensity is negative (min='.$array_min.')'."\n";
     print STDOUT '!!! '.$script.': possibly incorrect background data.'."\n";
  }

  if    ($bkgmode =~ /^zero$/i) {
     if ($array_max >= 0) {$half_intensity = 0.5*$array_max;}
     else {
        $pt_left  = 0;
        $pt_right = $npts-1;
        $x_left = $pt_left;
        $x_right = $pt_right;
        $x_peak = $pt_max;
        return $array_max,$array_min,$pt_max,$pt_min,$pt_FWHM,$pt_left,$pt_right,$x_left,$x_right,$x_peak;
     }
  }
  elsif ($bkgmode =~ /^auto$/i) {
     $array_span = $array_max - $array_min;
     $half_intensity = 0.5*$array_span + $array_min;
  } else {
     &press_key_exit('*** '.$script.': incorrect bkgmode='.$bkgmode,1);
  }

  if ($searchmode =~ /^smooth$/i) {
     ### Climb up on the peak from curve start:
     $pt_left = 0;
     while ((@$arrayRef[$pt_left] < $half_intensity) && ($pt_left < $pt_max)) {$pt_left++;}
     ### Climb up on the peak from curve end:
     $pt_right = $npts-1;
     while ((@$arrayRef[$pt_right] < $half_intensity) && ($pt_right > $pt_max)) {$pt_right--;}
  }          
  elsif ($searchmode =~ /^noisy$/i) {    
     ### Climb up on the peak from curve start:
     $pt_left = 0;
     while ((@$arrayRef[$pt_left]  < $half_intensity) && ($pt_left < $pt_max)) {$pt_left++;}
     ### Climb up on the peak from curve end:
     $pt_right = $npts-1;
     while ((@$arrayRef[$pt_right] < $half_intensity) && ($pt_right > $pt_max)) {$pt_right--;}
  }
  else {
     &press_key_exit('*** '.$script.': incorrect searchmode='.$searchmode,1);
  }
### If peak found:
  if (($pt_right < $npts-1) && ($pt_left > 0)) {
     $pt_FWHM = $pt_right - $pt_left + 1;
     $x_left = $pt_left;
     $dy = @$arrayRef[$pt_left] - @$arrayRef[$pt_left-1];
     if ($dy != 0.) {$x_left += ($half_intensity - @$arrayRef[$pt_left]) / $dy;}
     $x_right = $pt_right;
     $dy = @$arrayRef[$pt_right] - @$arrayRef[$pt_right+1];
     if ($dy != 0.) {$x_right -= ($half_intensity - @$arrayRef[$pt_right]) / $dy;}
     $x_peak = 0.5 * ($x_left + $x_right);
#    print STDOUT $script.":\n";
#    print STDOUT '   Half='.$half_intensity."\n";
#    print STDOUT '    Left='.@$arrayRef[$pt_left].'   '.@$arrayRef[$pt_left-1]."\n";
#    print STDOUT '   Right='.@$arrayRef[$pt_right].'  '.@$arrayRef[$pt_right+1]."\n";
#    print STDOUT '    pt_left='.$pt_left.'   x_left='.$x_left."\n";
#    print STDOUT '   pt_right='.$pt_right.'  x_right='.$x_right."\n";
  } else {
     $x_left  = $pt_left;
     $x_right = $pt_right;
     $x_peak  = $pt_max;
  }
  return $array_max,$array_min,$pt_max,$pt_min,$pt_FWHM,$pt_left,$pt_right,$x_left,$x_right,$x_peak;
}

#=================================================================

#sub findValleyRef ($npts, \@array, [$bkgmode,'smooth|noisy']) {
sub findValleyRef ($\@;$$) {
### This subroutine returns:
### ($array_max,$array_min,$pt_max,$pt_min,$pt_FWHM,$pt_left,$pt_right,[$x_left,$x_right,$x_valley])
### where $x_left,$x_right,$x_valley are respective interpolated value (i.e.
### they are a higher precision than $pt_left,$pt_right,$pt_FWHM).
### +---------------------------------------+
### |ATTENTION: It does not smooth anything!|
### +---------------------------------------+
### All the difference between "normal" and "noisy" valley
### finding is that in "normal" case we climb up from the valley
### to find FWHM, while in "noisy" we climb down from the tails.
###
  my ($script, $npts, $arrayRef, $bkgmode, $searchmode);
  my ($array_max, $array_min, $pt_max, $pt_min);
  my ($pt_FWHM, $pt, $array_span, $a_, $half_valley);
  my ($pt_left, $pt_right, $x_left, $x_right, $x_valley, $dy);

  $script     = (split(/::/,(caller(0))[3]))[-1];
  $npts       = shift(@_);
  $arrayRef   = shift(@_);
  $bkgmode    = shift(@_);               		#zero or auto(default) - differs from findPeak
  $searchmode = shift(@_);               		#smooth(default) or noisy
  if (!defined $bkgmode)    {$bkgmode    = 'auto';}
  if (!defined $searchmode) {$searchmode = 'smooth';}

  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0,-1,-1,0,0,0,0,0;
  }

  $array_max = -1.E37;				        #array maximum value
  $array_min = +1.E37;				        #array minimum value
  $pt_max    = -1;					#array index at maximum value
  $pt_min    = -1;					#array index at minimum value
  $pt_FWHM   = -1;					#valley width in pts (if exists)

### Find maximum/minimum:
  for ($pt=0; $pt < $npts; $pt++) {
     $a_ = @$arrayRef[$pt];
     if ($a_ > $array_max) {$array_max = $a_; $pt_max = $pt;}
     if ($a_ < $array_min) {$array_min = $a_; $pt_min = $pt;}
  }
  if ($array_max < 0) {
     print STDOUT '*** '.$script.': all intensity is negative (max='.$array_max.')'."\n";
     print STDOUT '*** '.$script.': possibly incorrect background data.'."\n";
  }
  elsif ($array_min < 0) {
     print STDOUT '!!! '.$script.': some intensity is negative (min='.$array_min.')'."\n";
     print STDOUT '!!! '.$script.': possibly incorrect background data.'."\n";
  }

  if    ($bkgmode =~ /^zero$/i) {
     if ($array_min >= 0) {$half_valley = 2.0*$array_min;}
     else                 {$half_valley = 0.5*$array_min;}
  }
  elsif ($bkgmode =~ /^auto$/i) {
     $array_span = $array_max - $array_min;
     $half_valley = 0.5*$array_span + $array_min;
  } else {
     &press_key_exit('*** '.$script.': incorrect bkgmode='.$bkgmode,1);
  }

  if ($searchmode =~ /^smooth$/i) {
     ### Climb up to the left of the valley:
     $pt_left = $pt_max;
     while ((@$arrayRef[$pt_left] < $half_valley) && ($pt_left > 0)) {$pt_left--;}
     ### Climb up to the right of the valley:
     $pt_right = $pt_max;
     while ((@$arrayRef[$pt_right] < $half_valley) && ($pt_right < $npts-1)) {$pt_right++;}
  }          
  elsif ($searchmode =~ /^noisy$/i) {    
     ### Climb down on the valley from curve start:
     $pt_left = 0;
     while ((@$arrayRef[$pt_left]  > $half_valley) && ($pt_left < $pt_max)) {$pt_left++;}
     ### Climb down on the valley from curve end:
     $pt_right = $npts-1;
     while ((@$arrayRef[$pt_right] > $half_valley) && ($pt_right > $pt_max)) {$pt_right--;}
  }
  else {
     &press_key_exit('*** '.$script.': incorrect searchmode='.$searchmode,1);
  }

### If valley found:
  if (($pt_right < $npts-1) && ($pt_left > 0)) {
     $pt_FWHM = $pt_right - $pt_left + 1;
     $x_left = $pt_left;
     $dy = @$arrayRef[$pt_left] - @$arrayRef[$pt_left-1];
     if ($dy != 0.) {$x_left += ($half_valley - @$arrayRef[$pt_left]) / $dy;}
     $x_right = $pt_right;
     $dy = @$arrayRef[$pt_right] - @$arrayRef[$pt_right+1];
     if ($dy != 0.) {$x_right -= ($half_valley - @$arrayRef[$pt_right]) / $dy;}
     $x_valley = 0.5 * ($x_left + $x_right);
#    print STDOUT $script.":\n";
#    print STDOUT '   Half='.$half_valley."\n";
#    print STDOUT '    Left='.@$arrayRef[$pt_left].'   '.@$arrayRef[$pt_left-1]."\n";
#    print STDOUT '   Right='.@$arrayRef[$pt_right].'  '.@$arrayRef[$pt_right+1]."\n";
#    print STDOUT '    pt_left='.$pt_left.'   x_left='.$x_left."\n";
#    print STDOUT '   pt_right='.$pt_right.'  x_right='.$x_right."\n";
  } else {
     $x_left   = $pt_left;
     $x_right  = $pt_right;
     $x_valley = $pt_min;
  }
  return $array_max,$array_min,$pt_max,$pt_min,$pt_FWHM,$pt_left,$pt_right,$x_left,$x_right,$x_valley;
}

#=================================================================

#($x_centroid,[$sum]) = &findCentroidRef($npts,\@arrx,\@arry,[$discriminator]);
#sub findCentroidRef ($npts, \@arrx, \@arry, [$discriminator]) {
sub findCentroidRef ($\@\@;$) {
### This subroutine returns $x_centroid
### ATTENTION: It does not smooth anything!
  my ($script, $npts, $arrxRef, $arryRef, $discriminator);
  my ($x_centroid, $sumxy, $sumx, $sumy, $i_warn);
  my ($x_, $y_, $pt, $discr, $y_max);
  my ($x1_, $y1_, $sum, $percent);

  $script        = (split(/::/,(caller(0))[3]))[-1];
  $npts          = shift(@_);
  $arrxRef       = shift(@_);
  $arryRef       = shift(@_);
  $discriminator = shift(@_);
  if (!defined $discriminator) {$discriminator=0;}

  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0;
  }

  $x_centroid = -1.E37;
  $sumxy      = 0.0;
  $sumx       = 0.0;
  $sumy       = 0.0;
  $sum        = 0.0;				#integral intensity above discriminator
  $i_warn     = 0;

### Find discriminated intensity level:
  if ($discriminator > 0) {
     $y_max =  -1.E37;
     for ($pt=0; $pt < $npts; $pt++) {
        $y_ = @$arryRef[$pt];
        if ($y_ > $y_max) {$y_max = $y_;}
     }
     $discr = $y_max * $discriminator;
  } else {
     $discr = 0.;
  }

### Find sums:
  for ($pt=0; $pt < $npts; $pt++) {
     $x_ = @$arrxRef[$pt];
     $y_ = @$arryRef[$pt]- $discr;
     if ($y_ < 0) {$i_warn++;  $y_ = 0;}
     $sumx += $x_;
     $sumy  += $y_;
     $sumxy += $x_ * $y_;
     if ($pt > 0) {$sum += 0.5*($y_ + $y1_) * ($x_ - $x1_);}
     $x1_ = $x_;
     $y1_ = $y_;
  }
  if ($npts > 0) {
     if ($sumy != 0.0) {$x_centroid = $sumxy/$sumy;}    #real centroid
     else              {$x_centroid = $sumx/$npts;}     #average x
  }

  if ($i_warn > 0) {
     $percent = (100.*$i_warn) / $npts;
     printf STDOUT '!!! '.$script.': intensity was '.$i_warn.' of '.$npts
                                  .' times (%1.0f%s) below discrimination level ('
                                  .$discriminator.') -- ignored'."\n",
                                   $percent,'%';
   }

  return $x_centroid, $sum;
}

#=================================================================

# ($x_centroid,[$sum]) = &findCentroidFixRef($npts,$x0,$dx,\@arry,[$discriminator]);
#sub findCentroidFixRef($npts,$x0,$dx,\@arry, [$discriminator]) {
sub findCentroidFixRef ($$$\@;$) {
### This subroutine returns $x_centroid for an equally spaced array
### ATTENTION: It does not smooth anything!
  my ($script, $npts, $x0, $dx, $arryRef, $discriminator);
  my ($x_centroid, $sumxy, $sumx, $sumy, $i_warn);
  my ($pt, $x_, $y_, $discr, $y_max);
  my ($y1_, $sum, $percent);

  $script        = (split(/::/,(caller(0))[3]))[-1];
  $npts          = shift(@_);
  $x0            = shift(@_);
  $dx            = shift(@_);
  $arryRef       = shift(@_);
  $discriminator = shift(@_);
  if (!defined $discriminator) {$discriminator=0;}

  if ($npts <= 0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0;
  }

  $x_centroid = -1.E37;
  $sumxy      = 0.0;
  $sumx       = 0.0;
  $sumy       = 0.0;
  $sum        = 0.0;
  $i_warn     = 0;

### Find discriminated intensity level:
  if ($discriminator > 0) {
     $y_max =  -1.E37;
     for ($pt=0; $pt < $npts; $pt++) {
        $y_ = @$arryRef[$pt];
        if ($y_ > $y_max) {$y_max = $y_;}
     }
     $discr = $y_max * $discriminator;
  } else {
     $discr = 0.;
  }

### Find sums:
  for ($pt=0; $pt < $npts; $pt++) {
     $x_ = $x0 + $dx*$pt;
     $y_ = @$arryRef[$pt] - $discr;
     if ($y_ < 0) {$i_warn++;  $y_ = 0;}
     $sumx += $x_;
     $sumy  += $y_;
     $sumxy += $x_ * $y_;
     if ($pt > 0) {$sum += 0.5*($y_ + $y1_) * $dx;}
     $y1_ = $y_;
  }
  if ($npts > 0) {
     if ($sumy != 0.0) {$x_centroid = $sumxy/$sumy;}    #real centroid
     else              {$x_centroid = $sumx/$npts;}     #average x
  }

  if ($i_warn > 0) {
     $percent = (100.*$i_warn) / $npts;
     printf STDOUT '!!! '.$script.': intensity was '.$i_warn.' of '.$npts
                                  .' times (%1.0f%s) below discrimination level ('
                                  .$discriminator.') -- ignored'."\n",
                                   $percent,'%';
  }

  return $x_centroid, $sum;
}

#=================================================================

#sub smoothArray (\@array) {
sub smoothArray (\@) {
# May be used by any curve analysis program.
# This routine does weighted array averaging over 7 nearest points.
# The algorithm is due to Alexander Chuzo.

  my ($arref, @w, @q, $d, $i, $j, $k, $nPts);

  $arref = shift(@_);
  @w = (0.005,0.025,0.625,1.,0.625,0.025,0.005);	# weights
  @q = ();						#cached arrx
  $d = 0.0;						#divider

  no strict 'refs';
  $nPts = @$arref;
  if ($nPts < 4) {return;}		    		#nothing to smooth

  for($j=-3; $j<=3; $j++) {                 		#since we are overwriting the array
    if ($j >= 0) {$q[3+$j] = @$arref[$j];}  		#we need to make a cached copy
    else         {$q[3+$j] = 0.0;}
  }

  for ($i=0; $i<$nPts; $i++) {
     @$arref[$i] = 0.0;
     $d          = 0.0;
     for($j=$i-3; $j<=$i+3; $j++) {
       if (($j >= 0) && ($j < $nPts)) {
          $k           = $j-($i-3);
          @$arref[$i] += $w[$k] * $q[$k];
          $d          += $w[$k];
       }
     }
     if ($d != 0.) {@$arref[$i] = @$arref[$i]/$d;}
     else          {@$arref[$i] = 0.0;}

     for ($k=0; $k<6; $k++) {$q[$k]=$q[$k+1];} 		#shift cached copy of arr2 by one element
     if ($i<$nPts-4) {$q[6] = @$arref[$i+4];}  		#add next element to the cache (if exists)
     else            {$q[6] = 0.0;}
  }
  use strict 'refs';
}

#=================================================================

#sub arrayStatsRef (\@array) {
sub arrayStatsRef (\@) {
### This subroutine returns array statistics: ($pt_min,$y_min,$pt_max,$y_max)

  my ($script, $arrayRef, $npts, $y_min, $y_max, $pt_min, $pt_max, $pt);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $arrayRef = shift(@_);
  $npts     = @$arrayRef;

  if ($npts <=0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0,0,0,0;
  }

  $y_min    = +1.E37;				#array minimum value
  $y_max    = -1.E37;				#array maximum value
  $pt_min   = -1;				#array index at minimum value
  $pt_max   = -1;				#array index at maximum value
### Find maximum
  for ($pt=0; $pt < $npts; $pt++) {
     if (@$arrayRef[$pt] > $y_max) {$y_max=@$arrayRef[$pt]; $pt_max=$pt;}
     if (@$arrayRef[$pt] < $y_min) {$y_min=@$arrayRef[$pt]; $pt_min=$pt;}
  }
  return $pt_min,$y_min,$pt_max,$y_max;
}

#=================================================================

#sub arrayStats ($npts, @array) {
sub arrayStats ($@) {
# This subroutine returns array statistics: ($pt_min,$y_min,$pt_max,$y_max)

  my ($script, $npts, @array, $narg);
  my ($y_min, $y_max, $pt_min, $pt_max, $pt, $pt1);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $narg   = @_;
  $npts   = shift(@_);

  if ($npts <=0) {
     print STDERR "\n".'!!! '.$script.': received zero-length array!'."\n\n";
     return 0, 0, 0, 0;
  }
  elsif ($npts > $narg-1) {
     printf STDERR "\n".'!!! '.$script.': received %d array elements instead of %d declared!'."\n\n",($narg-1),$npts;
     $npts = $narg-1;
  }

  $y_min  = +1.E37;				#array minimum value
  $y_max  = -1.E37;				#array maximum value
  $pt_min = -1;					#array index at minimum value
  $pt_max = -1;					#array index at maximum value
### Find maximum
  for ($pt=0; $pt<$npts; $pt++) {
     $array[$pt] = shift(@_);
     if (!defined $array[$pt]) {
        $pt1 = $pt+1;
        print STDOUT '*** '.$script.': undefined array for pt='.$pt1.' of '.$npts."\n";
     }
     else {
        if ($array[$pt] > $y_max) {$y_max=$array[$pt]; $pt_max=$pt;}
        if ($array[$pt] < $y_min) {$y_min=$array[$pt]; $pt_min=$pt;}
     }
  }
  return $pt_min,$y_min,$pt_max,$y_max;
}

#=================================================================

#sub &analyzeFilename ($file) {			#returns ($name,$ext)
sub analyzeFilename ($) {
### This subroutine returns name and extension of a given file

  my ($file, $name, $ext, @a, $n);

  $file = shift(@_);
  $name = $file;
  $ext  = '';
  @a = split /\./, $file;
  $n = @a;
  if ($n > 1) {
     $ext  = $a[$n-1];
     $name =~ s/\.$ext$//;
  }
  return $name, $ext;
}

#=================================================================

#sub makeScriptName ([1|'short'|'noext']) {		#returns truncated script name (no path)
sub makeScriptName (;$) {
  my ($script, $short);
  $short = shift(@_); if (!defined $short) {$short=0;}
# $script = $0;
# $script =~ s/\\/\//g;				#change \ to /
# $script =~ s/^.*\///g;			#remove everything before last "/"
  $script = $0 =~ s/^.*[\/\\]//gr;		#same as above 3 in 1 line: remove everything before last "/" or "\"

  if ($short eq '1'        || 
      $short =~ /^short$/i || 
      $short =~ /^noext/i) {$script =~ s/\.pl$//i;}
  return $script;
}

#=================================================================

#sub getFileAgeHours ($file) {  	#returns $age or -1 if file does not exist
sub getFileAgeHours($) {

  my ($script, $file, $age_hours);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $file = shift(@_);
  if (-f $file) {				#if file exists:
     $age_hours = 24 * (-M $file);
     if ($age_hours < 0) {
        print STDOUT '!!! '.$script.': file "'.$file.'" has modification'
                                   .' timestamp in future'."\n";
        return 0;
     } else {
        return $age_hours;
     }
  } else {
     return -1;					#file does not exist
  }
}

#=================================================================

#sub isMounted ($mountPoint) {	#returns 1(yes) is mount point is mounted; otherwise 0(no)
sub isMounted ($) {
  my ($script, $mountPoint, $procMountFile, $handle, $record, @words, $nwords);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $mountPoint = shift(@_);
  if (!defined $mountPoint) {
     print STDOUT '!!! '.$script.': unspecified mountpoint'."\n";
     return 1;
  }
  if ($^O !~ /Win/i) {
### Proc filesystem which tracks current mounts on Linux:
     $procMountFile = '/proc/mounts';
     $handle = open (MOUNTLIST, '< '.$procMountFile);
     if (!defined $handle) {
        print STDOUT '!!! '.$script.': cannot open "'.$procMountFile.'" -- '
                    .'presuming mountpoint='.$mountPoint.' is mounted'."\n";
        return 1;
     }
     while (<MOUNTLIST>) {
        $record = $_;
        chomp $record;				# trail CR
        @words = split " ", $record;
        $nwords = @words;
        if ($nwords >=2) {
           if ($words[1] =~ /^$mountPoint$/){
### The mount is found
              close(MOUNTLIST);
              return 1;
           }
        }
     }
### The mount is NOT found
     close(MOUNTLIST);
     return 0;
  }
### Microsoft Windows:
  return 1;
}

#=================================================================

#($status,$msg) = &read_data_file ($file,\@xdata,\@ydata,[$xcol,$ycol]);
sub read_data_file ($\@\@;$$) {
  my ($script, $file, $xref, $yref, $xcol, $ycol, $maxcol, $ncol);
  my ($handle, $line, $npts, $string, @arrline, $msg);
  $script = (split(/::/,(caller(0))[3]))[-1];
  $file   = shift(@_);
  $xref   = shift(@_);
  $yref   = shift(@_);
  $xcol   = shift(@_); if (!defined $xcol) {$xcol=1;}
  $ycol   = shift(@_); if (!defined $ycol) {$ycol=2;}
  if ($xcol == $ycol) {
     $msg = '*** '.$script.': Xcol='.$xcol.' and Ycol='.$ycol.' are the same';
     return 1,$msg;
  }
  $maxcol = &maxval(($xcol,$ycol));

### Open data file:
  if (-f $file) {
     print STDOUT $script.': Opening file="'.$file.'", xcol='.$xcol.', ycol='.$ycol."\n";
     $handle = open (DAT, '< '.$file);
     if (!defined $handle) {
        $msg = '*** '.$script.': Error opening data file='.$file;
        return 1,$msg;
     }
#    select DAT; $|=1;
  } else {
     $msg = '*** '.$script.': data file: '.$file.' does not exist';
     return 1,$msg;
  }
### Read the file ignoring the comment lines:
  $line = 0;
  $npts = 0;
  while (<DAT>) {						#begin loop over lines in file
     $line++;
     if (!defined $_) {next;}
     $string = $_;						#read line
     chomp $string;						#trail CR
     $string =~ s/\#.*$//g;					#trail comments started by "#" ("."=any char, "*"=any times, "$"=EOL)
     $string =~ s/,/ /g;					#replace comas by spaces
     $string =~ s/\t/ /g;					#replace tabs by spaces
     $string =~ s/ +/ /g;					#compress spaces ("+" means match 1 or more times)
     $string =~ s/^ *//g;					#remove beginning-of-the-line spaces
     $string =~ s/ *$//g;					#remove end-of-the-line spaces
     if (length($string) == 0) {next;}
     @arrline = split ' ', $string;  				#split the line into words
     $ncol = @arrline;
     if ($ncol < $maxcol) {
        $msg = '*** '.$script.': data file '.$file.' has '.$ncol.' < '.$maxcol.' columns at line='.$line;
        return 1,$msg;
     }
     if ($xcol > 0) {$$xref[$npts] = $arrline[$xcol-1];}	#x-data
     else           {$$xref[$npts] = $npts;}            	#x-data
     $$yref[$npts] = $arrline[$ycol-1];				#y-data
     $npts++;
  }								#end of loop over lines in data file
  close(DAT);							#close data file
  return 0;
}	##read_data_file


#=================================================================

#$epoch = &timestamp2epoch($timestamp); 		#expected: YYYY-MM-DD HH:MM[:SS]
sub timestamp2epoch ($) {
  my ($script, $timestamp, $g1, $g2, $epoch);
  my (@ct, $el);					#$sec,$min,$hour,$mday,$mon,$year

  $script = (split(/::/,(caller(0))[3]))[-1];
  $timestamp = shift(@_);
  $epoch = 0;
  if ((length($timestamp) != 19 && 			#length(YYYY-MM-DD HH:MM:SS)
       length($timestamp) != 16) || 			#length(YYYY-MM-DD HH:MM)
       $timestamp !~ /[-\/_]/    || 			#e.g. YYYY-MM-DD or YYYY/MM/DD or YYYY_MM_DD
       $timestamp !~ ':'         ||			#e.g. HH:MM
       $timestamp !~ ' ') {				#should contain space between date & time
     print STDOUT '!!! '.$script.': incorrect time stamp ['.$timestamp.']'."\n";
     return $epoch;
  }
  ($g1,$g2) = split ' ', $timestamp;			#split into date & time
  ($ct[5], $ct[4], $ct[3]) = split /[-\/_]/, $g1;	#($year, $month, $mday)
  ($ct[2], $ct[1], $ct[0]) = split ':',      $g2;	#($hour, $min, $sec)
  if (!defined $ct[0]) {$ct[0] = '00';}
  foreach $el (@ct) {
     if (length($el) > 4 || length($el) < 2 || $el !~ /^\d+$/) {
        print STDOUT '!!! '.$script.': incorrect time stamp ['.$timestamp.']'."\n";
        return $epoch;
     }
  }
  $ct[5] -= 1900;
  $ct[4] -= 1;
  if (!defined &function_exists('timelocal')) {
     if (&try_load_module('Time::Local')) {return 0;}
  }
  $epoch = timelocal($ct[0],$ct[1],$ct[2],$ct[3],$ct[4],$ct[5]);	#use Time::Local
  return $epoch;
}


#=================================================================

#$timestamp = &epoch2timestamp([$epoch,{/-},'nosec']); 		#returns: [YYYY-MM-DD HH:MM:SS] or [YYYY/MM/DD HH:MM:SS]
sub epoch2timestamp (;$$$) {
  my ($script, $epoch, $sep, $nos, $timestamp);
  my (@ct);						#$sec,$min,$hour,$mday,$mon,$year
                                                        # 0    1     2     3     4    5

  $script = (split(/::/,(caller(0))[3]))[-1];
  $epoch  = shift(@_);
  $sep    = shift(@_);  if (!defined $sep || $sep ne '/')    {$sep='-';}
  $nos    = shift(@_);  if (!defined $nos || $nos !~/^nos/i) {$nos='seconds';}

  if (defined $epoch)  {@ct = localtime($epoch);}
  else                 {@ct = localtime();}
  
  if ($nos eq 'seconds') {
     $timestamp = sprintf '%04d%s%02d%s%02d %02d:%02d:%02d',($ct[5]+1900),$sep,($ct[4]+1),$sep,$ct[3], $ct[2],$ct[1],$ct[0];
  } else {
     $timestamp = sprintf '%04d%s%02d%s%02d %02d:%02d',     ($ct[5]+1900),$sep,($ct[4]+1),$sep,$ct[3], $ct[2],$ct[1];
  }

  return $timestamp;
}


#=================================================================


# $status = &sort_by_x('incr|decr',\@x,@y1,[\@y2,\@y3,...]);
sub sort_by_x ($\@\@;\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@) {
  my ($script, $type, $xref, @yref, $npts, $pt, $j, $j1, $tmp, $nSort, $y);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $type = shift(@_);
  $xref = shift(@_);
  @yref = @_;				#the rest: array of pointers to y-arrays
  $npts = @$xref;			#array length
  if ($type !~ /^incr/i && $type !~ /^decr/) {
     print STDOUT '*** '.$script.': incorrect sorting type='.$type."\n";
     return -1;
  }
  if ($npts < 2) {return 0;}		#no sorting required
  $nSort = 0;
  if ($type =~ /^incr/i) {
#    print STDOUT $script.': sorting data in increasing order...'."\n";
     for ($pt=0; $pt<($npts-1); $pt++) {
        for ($j=0; $j<($npts-1-$pt); $j++) {
           $j1 = $j + 1;
           if (@$xref[$j1] < @$xref[$j]) {
              $nSort++;
              $tmp        = @$xref[$j];
              @$xref[$j]  = @$xref[$j1];
              @$xref[$j1] = $tmp;
              foreach $y (@yref) {
                 $tmp     = @$y[$j];
                 @$y[$j]  = @$y[$j1];
                 @$y[$j1] = $tmp;
              }
           }
        }
     }
  } else {
#    print STDOUT $script.': sorting data in decreasing order...'."\n";
     for ($pt=0; $pt<($npts-1); $pt++) {
        for ($j=0; $j<($npts-1-$pt); $j++) {
           $j1 = $j + 1;
           if (@$xref[$j1] > @$xref[$j]) {
              $nSort++;
              $tmp        = @$xref[$j];
              @$xref[$j]  = @$xref[$j1];
              @$xref[$j1] = $tmp;
              foreach $y (@yref) {
                 $tmp     = @$y[$j];
                 @$y[$j]  = @$y[$j1];
                 @$y[$j1] = $tmp;
              }
           }
        }
     }
  }
  if ($nSort > 0) {
#    print STDOUT '!!! '.$script.': data has been re-ordered.'."\n";
     return 1;
  } else {
#    print STDOUT '!!! '.$script.': data sorting was not required.'."\n";
     return 0;
  }
}


#=================================================================


# $npts = &remove_duplicates($filter,$step,\@x,@y1,[\@y2,\@y3,...]);
sub remove_duplicates($$\@\@;\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@) {

### The arrays must be sorted over x and 'quasi-regular', i.e.
### $x(last)-$x(first)/($npts-1) is about an average step

  my ($script, $fltr, $step, $xref, @yref, $npts, $pt, $y, $diff, $nDup);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $fltr = shift(@_);
  $step = shift(@_);
  $xref = shift(@_);
  @yref = @_;
  $npts = @$xref;
  if ($npts < 2) {return $npts;}			#no duplicates here

  if (!defined $step || $step =~ /^auto$/i) {
     $step = (@$xref[$npts-1] - @$xref[0])/($npts-1);	#presumes 'quasi-regular' data spacing
  }
  if ($fltr <=0 || $fltr > 0.5) {
     print STDOUT '*** '.$script.': duplicates filter='.$fltr.' is not in range [0.0-0.5]'."\n";
     return $npts;
  }
# print STDOUT $script.': filter='.$fltr.'  step='.$step."\n";
  $fltr = $fltr*abs($step);				#points spaced below 20% of step are duplcates

### Remove duplicates:
  $nDup   = 0;
  $pt     = 1;
  while ($pt < $npts) {
     $diff = @$xref[$pt] - @$xref[$pt-1];
     if (abs($diff) < $fltr) {
        @$xref[$pt-1] = 0.5*(@$xref[$pt] + @$xref[$pt-1]);
#       print STDOUT $script.': filter*step='.$fltr.'  diff='.$diff.'  x(average)='.@$xref[$pt-1]."\n";
        splice(@$xref, $pt, 1);				#compress x array
        foreach $y (@yref) {
           @$y[$pt-1] = 0.5*(@$y[$pt] + @$y[$pt-1]);
           splice(@$y, $pt, 1);				#compress data arrays
        }
        $pt--;
        $npts--;
        $nDup++;
     }
     $pt++;
  }
  if ($nDup > 0) {
     print STDOUT '!!! '.$script.': '.$nDup.' duplicate data points removed'
                        .' ('.$npts.' total remaining)'."\n";
  }
  return $npts;
}


#=================================================================

#sub Which ($prg) {	#returns path to $prg or undef if not found. On win32 tries .exe,.bat,.com
sub Which ($) {

  my ($program, $path, @parr, $file, $os, @ext, $ex, $filewin);
  $program = shift;				#get the passed value
  return if (not defined($program));		#return if nothing is provided
  if ($^O !~ /Win/i) {$os = 'unix';}
  else               {$os = 'win32';}
  @ext = ('.exe', '.bat', '.com');
# use ENV;					#load the ENV
  $path = $ENV{PATH};				#load the path
  $path =~ s/\\/\//g;				#replace all \ by /
  if ($os eq 'unix') {$path =~ s/:/;/g;}	#replace all ":"  by ";"
  else               {$path = './;'.$path;}	#on Windows, also add current directory
  $path =~ s/\/;/;/g;				#replace all "/;" by ";"
  @parr = split(/;/,$path);			#now make an array (split by ";"_
  foreach (@parr) {				#loop and find if the file is in one of the paths
    $file = $_ . '/' . $program;		#concatenate the file
    #print STDOUT '  Which: checking '.$file."\n";
    return $file if ((-e $file) && (-f $file));	#return the path if it was found
    if ($os eq 'win32') {
      foreach $ex (@ext) {			#loop over possible extensions
        if ($file !~ /${ex}$/i) {
          $filewin = $file.$ex;
          return $filewin if ((-e $filewin) && (-f $filewin));
        }
      }
    }
  }
  return;
}


#=================================================================

#sub split_filename_dirname ($filespec) {	#returns ($dirname,$filename)
sub split_filename_dirname ($) {
  my ($filespec, $dirname, $filename);
  $filespec = shift(@_);
  if ($filespec =~ m|[/\\]$|) {			#if ends by '\' or '/'
     ($dirname, $filename) = ($filespec, '');
  } else {
    ($dirname, $filename) = ($filespec =~ m|^(.*[/\\])([^/\\]+?)$|);	#split by  '\' or '/'
    if (!defined $dirname && $^O =~ /Win/i) {				#if unsuccessful & Windows
       ($dirname, $filename) = ($filespec =~ m|^([a-zA-Z]\:)(.*)$|);	#split by drive name (e.g. 'c:')
    }
    if (!defined $dirname) {($dirname, $filename) = ('',$filespec);}
  }
  return ($dirname, $filename);
}


#=================================================================

# $radians = &degrees2rad($degrees);		#this is a replacement to deg2rad of Math::Trig
sub degrees2rad($) {
  my ($degrees, $radians, $pi);
  $degrees = shift(@_);
  $pi = 4.*atan2(1.,1.);
  if (!defined $pi) {$pi = 3.1416;}
  $radians = $degrees*$pi/180.;
  return $radians
}


#=================================================================


# $degrees = &rad2degrees($radians);		#this is a replacement to rad2deg of Math::Trig
sub rad2degrees($) {
  my ($degrees, $radians, $pi);
  $radians = shift(@_);
  $pi = 4.*atan2(1.,1.);
  if (!defined $pi) {$pi = 3.1416;}
  $degrees = $radians*180./$pi;
  return $radians
}

#=================================================================


#sub find_first_index($element,\@array) {
sub find_first_index($\@) {			#returns index of element in array or -1
  my ($element, $arrayRef, $ind, $narr);

  $element  = shift(@_);
  $arrayRef = shift(@_);

  $narr = @$arrayRef;
  $ind = -1;
  do {$ind++;} while ($ind < $narr && @$arrayRef[$ind] ne $element);
  if ($ind < $narr) {return $ind;} else {return -1;}
}


#=================================================================

#sub find_all_index($element,\@array) {
sub find_all_index($\@) {			#returns index of element in array or -1
  my ($element, $arrayRef, @ind, $ilast);

  $element  = shift(@_);
  $arrayRef = shift(@_);

  $ilast = @$arrayRef-1;
  @ind = grep{@$arrayRef[$_]=~/^$element$/} 0..$ilast;
  return @ind;
}


#=================================================================


#sub find_closest_index($val,\@array) {
sub find_closest_index($\@) {			#returns index of element in array closest to val
  my ($val, $arrayRef, $diff, $d, $ind, $narr, $i);

  $val      = shift(@_);
  $arrayRef = shift(@_);

  $narr = @$arrayRef;
  $diff = 1.e+100;
  $ind = 0;
  for ($i=0; $i<$narr; $i++) {
     $d = abs(@$arrayRef[$i] - $val);
     if ($d < $diff) {$diff=$d; $ind=$i;}
  }
  return $ind;
}


#=================================================================

# sub execute_external_script(\@cmd) {		#perl, -W, $script, $arg1, ...
sub execute_external_script(\@) {
  our ($EXTERNAL_PID, $EXTERNAL_PRG);
  my ($cmdref, $script, $pid, $kid, $WNOHANG);
  my ($code, $arg, @subpid, $not_perl);
  my ($cmd, $msg, @processes, $p_, @words, $nw);

  $cmdref = shift(@_);
  undef $EXTERNAL_PID;
  undef $EXTERNAL_PRG;

### WNOHANG is a constant defined in POSIX:
# use POSIX ":sys_wait_h";
### This constant is actually "1". So instead
### of using POSIX we define a variable:
  $WNOHANG = 1;

  $script = &makeScriptName();

  $pid = fork;
  if (!defined $pid) {
     print STDOUT '*** '.$script.': cannot fork to execute '."@$cmdref\n";
     return 1;
  }
  elsif ($pid == 0) {
     ### Client process: execute external program
     print STDOUT $script.':  ==> Calling:'."\n@$cmdref\n";
     exec(@$cmdref);							#this will completely replace the child with called program, but retain PID
  }
  else {								#monitor the script until it exits (also prepare data for signal handler)
     ### We are in parent process
     ### 1. Find the script name
     foreach $arg (@$cmdref) {						#find the name of executed script
        if ($arg =~ /\.pl$/) {
           $EXTERNAL_PRG = $arg;
           $EXTERNAL_PRG =~ s/^.*\///g;
           last;
        }
     }
     if (!defined $EXTERNAL_PRG) {
        $EXTERNAL_PRG = @$cmdref[0];
        $not_perl = 1;
     }
     ### 2. Find the script PID
     if (@$cmdref[0] =~ /perl$/ || 					#direct call of perl script
         @$cmdref[0] =~ /perl\.exe$/i ||                                #direct call of perl script
         defined $not_perl) {	                                        #we are calling not perl
        $EXTERNAL_PID = $pid;						#this is for aborting by signal handler
     } else {								#script is started in terminal, find the PID of the script
        if ($^O !~ /Win/i) {
           @subpid = `/usr/bin/pgrep -P $pid`;				#find all subprocesses for given parent
        } else {
           @subpid = ();
           $cmd = 'wmic process where name="perl.exe" get ProcessId,CommandLine /format:csv '
           . '| grep -iv "grep\|zenity" '
           . '| grep -i '.$EXTERNAL_PRG;
           $msg = `$cmd`;
           $msg =~ s/\'//g;						#remove all "'"
           @processes = split(/\n/, $msg);
           foreach $p_ (@processes) {
              @words = split(',', $p_);
              $nw = @words;
              if ($words[$nw-1] =~ /^\d+$/) {
                 push(@subpid,$words[$nw-1]);
              } else {
                 print $script.': Unexpected PID='.$words[$nw-1].' for process='.$p_."\n";
              }
           }
        }
        $nw = @subpid;
        if ($nw > 0) {
           $EXTERNAL_PID = $subpid[0];					#take the first found only. Is it right?
           if ($nw != 1) {
              print STDOUT '!!!'.$script.': found '.$nw.' child PIDs for script '.$EXTERNAL_PRG."\n";
           }
        } else {
           $EXTERNAL_PID = $pid;
        }
     }
     ### 3. Wait for child to exit
     do {						#wait the script to finish
        select(undef,undef,undef,0.05);
        $kid = waitpid($pid, $WNOHANG);			#wait for particular child
     } while ($kid == 0);				#this should be "0" while process runs
     $code = $? >> 8;
     undef $EXTERNAL_PID;
     undef $EXTERNAL_PRG;
     if ($code) {
        print STDOUT '!!!'.$script.': execution of ['."@$cmdref".'] finished with status='.$code."\n";
     }
     return $code;
  }
}

#=================================================================

#sub abort_external_script (;$pid) {
sub abort_external_script(;$) {
  our ($EXTERNAL_PID, $EXTERNAL_PRG);
  my  ($script, $pid, $kid, $nkilled, $code);
  my  ($wait_, $wait_step, $wait_max, $WNOHANG);

### WNOHANG is a constant defined in POSIX:
# use POSIX ":sys_wait_h";
### This constant is actually "1". So instead
### of using POSIX we define a variable:
  $WNOHANG = 1;

  $pid = shift(@_);   if (defined $pid) {$EXTERNAL_PID = $pid;}

  if (!defined $EXTERNAL_PID) {return 0;}				#nothing to do
  else                        {$pid = $EXTERNAL_PID;}
  if (!defined $EXTERNAL_PRG ) {$EXTERNAL_PRG = 'external script';}

  $script = &makeScriptName();
  printf STDOUT '!!! '.$script.': killing [%s] with PID=%s'."\n", $EXTERNAL_PRG, $EXTERNAL_PID;

  $nkilled = kill('SIGINT',$EXTERNAL_PID);			#sends to main process
# printf STDOUT '!!! '.$script.': killed %d processes.'."\n", $nkilled;

### The script we are trying to kill may have its own signal handler and
### it may caLL other scripts which it may try to kill on exit. Therefore
# we need to wait until it finishes:
  $wait_     = 0;
  $wait_step = 0.05;
  $wait_max  = 30;
  do {								#wait the script to finish
     select(undef,undef,undef,$wait_step);
     $wait_ += $wait_step;
     $kid = waitpid($pid, $WNOHANG);				#wait for particular child
  } while ($kid == 0 && $wait_ < $wait_max);			#this should be "0" while process runs
  $code = $? >> 8;
  undef $EXTERNAL_PID;
  undef $EXTERNAL_PRG;
  if ($code) {
     print STDOUT '!!!!'.$script.': execution of '.$EXTERNAL_PRG.' finished with status='.$code."\n";
  }
  return $code;
}

#===========================================================================

sub zenity_progress($$) {			# $seconds,$msg

### This routine displays countdown window. On Linux it
### uses zenity, but on Windows pipes do not take command
### line argumnents and then this script requires
### Tk and Tk::ProgressBar to display the progress

  our ($zenity_pid, $mw, $mw_progress_label);
  our ($mw_progress_bar, $mw_progress_abort);
  my ($script, $waitmax, $wait_, $msg, $zenity, @zenargs);
  my ($start_time, $elapsed, $gui, $percent_remaining);

  $script = (split(/::/,(caller(0))[3]))[-1];

  $waitmax = shift(@_);		if ($waitmax < 1) {return;}
  $msg     = shift(@_);		$msg =~ s/[\n\r]/ /g;

  $start_time = time();
  if ($^O !~ /Win/i) {$zenity = &Which('zenity');}
  else               {$zenity = &Which('zenity.exe'); undef $zenity;} #cygwin zenity works, but requires X11

  if (defined $zenity) {
     ### zenity --progress --text="bla-bla-aaaaaaaaaaaaaaa ee\rbla" --title="Beamstop-IN" --percentage=99
     $wait_ = sprintf $msg.'. Remaining time %ds', $waitmax;
     @zenargs = ('--progress'
                ,'--title=Progress'
                ,'--text="'.$wait_.'"'
                ,'--percentage=99'
                );
     if ($^O !~ /Win/i) {
        ### The canonical version does not suppress Gtk errors, even with '2>/dev/null'.
        ### Therefore, we must use Open3, which is a bidirectional pipe catching both
        ### zenity STDOUT and zenity STDERR.
        # push @zenargs, '2>/dev/null';
        # $zenity_pid = open(ZPIPE, '|-', $zenity, @zenargs);
        use IPC::Open3;
        $zenity_pid = open3(\*ZPIPE, \*CHLD_OUT, \*CHLD_ERR, $zenity, @zenargs);
     } else {
        $zenity_pid = open(ZPIPE, '|-', join(' ',$zenity,@zenargs));		#only string works on Windows
     }
     if (defined $zenity_pid) {                                              #if pipe to zenity opened OK
        select ZPIPE; $|=1;							#set unbuffered output
        $SIG{PIPE} = sub {our $zenity_pid;
                          undef $zenity_pid;
                          close(CHLD_OUT);
                          close(CHLD_ERR);
                         };
        do {
           select(undef,undef,undef,0.5);
           $elapsed = time() - $start_time;
           $wait_ = sprintf $msg.'. Remaining time %ds', ($waitmax-$elapsed);
           $percent_remaining = sprintf '%d', 100*($waitmax-$elapsed)/$waitmax;
           if ($percent_remaining > 99) {$percent_remaining=99;}
           if (waitpid($zenity_pid,1) !=0 ) {undef $zenity_pid;}		#no wait
           if (defined $zenity_pid) {
              print ZPIPE '# '.$wait_."\n".$percent_remaining."\n";
           } else {
              print STDOUT $wait_."\n";
           }
        } while ($elapsed < $waitmax && defined $zenity_pid);
        if (defined $zenity_pid) {
           kill('KILL',$zenity_pid);
           close(ZPIPE);
           return 0;
        } else {
           $script = &makeScriptName();
           $wait_ = sprintf '%ds', ($waitmax-$elapsed);
           print STDOUT '!!! '.$script.': Progress aborted '.$wait_.' before completion.'."\n\n";
           return 1;
        }
     } else {                                                                #pipe to zenity did not open
        do {
           select(undef,undef,undef,0.5);
           $elapsed = time() - $start_time;
           $wait_ = sprintf $msg.'. Remaining time %ds', ($waitmax-$elapsed);
           print STDOUT $wait_."\n";
        } while ($elapsed < $waitmax);
     }
  } else {                                                                   #there is no zenity app on the system
#----------------------------------------------#
     if (!defined &function_exists('Exists')) {
        # use Tk;
        # use Tk::ProgressBar;
        if (&try_load_module('Tk'))              {die 'cannot load Tk';}
        if (&try_load_module('Tk::ProgressBar')) {die 'cannot load Tk::ProgressBar';}
     }
#----------------------------------------------#
     if (! Exists($mw)) {
        $mw = MainWindow->new(-title => 'Progress');
        $mw->geometry( "600x150" );
        $gui = 0;
     } else {
        $gui = 1;
     }
     $zenity_pid = 1;
     $percent_remaining = 100;
     $wait_ = sprintf $msg.'. Remaining time %ds', $waitmax;
     $mw_progress_label = $mw->Label(-text =>$wait_,
                          )->pack(-padx=>10,-pady=>10);

     $mw_progress_bar = $mw->ProgressBar(
                       -width => 30,
                       -from => 0,
                       -to => 100,
                       -blocks => 10,
                       -gap => 0,
                       -colors => [0, 'blue'],
                       -variable => \$percent_remaining,
                       )->pack(-fill => 'x',-padx=>10,-pady=>10);

     $mw_progress_abort = $mw->Button(-text => 'Abort',
                             -command=> sub {
                                $mw_progress_label->destroy();
                                $mw_progress_abort->destroy();
                                $mw_progress_bar->destroy();
                                if (!$gui) {$mw->destroy();}
                                undef $zenity_pid;
                             },
                          )->pack(-padx=>10,-pady=>10);
     do {
        select(undef,undef,undef,0.5);
        $elapsed = time() - $start_time;
        $wait_ = sprintf $msg.'. Remaining time %ds', ($waitmax-$elapsed);
        $percent_remaining = sprintf '%d', 100*($waitmax-$elapsed)/$waitmax;
        $mw_progress_label->configure(-text=>$wait_);
        $mw->update;
#       print STDOUT $wait_."\n";
     } while ($elapsed < $waitmax && defined $zenity_pid);
     if (defined $zenity_pid) {
        $mw_progress_label->destroy();
        $mw_progress_abort->destroy();
        $mw_progress_bar->destroy();
        if (!$gui) {$mw->destroy();}
        undef $zenity_pid;
        return 0;
     } else {
        $script = &makeScriptName();
        $wait_ = sprintf '%ds', ($waitmax-$elapsed);
        print STDOUT '!!! '.$script.': Progress aborted '.$wait_.' before completion.'."\n\n";
        return 1;
     }
     no strict "subs";
     if (!$gui) {MainLoop;}
     use strict "subs";
  }
}

#===========================================================================

sub is_float($) {					# returns 1/0 (True/False)
# use Data::Types;					# may conflict with this module (same names)
  my $str = shift(@_);
  if ($str =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {return 1;}
  else                                                            {return 0;}
}

#===========================================================================

sub is_decimal($) {					# returns 1/0 (True/False)
# use Data::Types;					# may conflict with this module (same names)
  my $str = shift(@_);
### Similar to is_float(), but allows fixed-point decimals only (no exponents)
  if ($str =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)$/) {return 1;}
  else                                          {return 0;}
}

#===========================================================================

sub is_int($) {						# returns 1/0 (True/False)
# use Data::Types;					# may conflict with this module (same names)
  my $str = shift(@_);
  if ($str =~ /^[+-]?\d+$/) {return 1;}
  else                      {return 0;}
}

#===========================================================================

sub is_mounted ($;$) {
  use Cwd 'abs_path';
  my ($dir, $ref, @stat1, @stat2);
# my ($dirlink);
  $dir = shift(@_);
  $ref = shift(@_);
  if (!defined $ref) {$ref='/';}
  if (! -e $dir) {return -1;}
  if (! -e $ref) {return -1;}
# $dirlink = readlink($dir);				#this does not if link is relative because it will return relative path and stat will fail
# if (defined $dirlink) {$dir = $dirlink;}
  $dir   = abs_path($dir);
  $ref   = abs_path($ref);
# ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks)=stat($filename);
  @stat1 = stat($ref);
  @stat2 = stat($dir);
  if ($stat2[0] != $stat1[0]) {return 1;}		#device numbers do not coincide
  else                        {return 0;}
}

#===========================================================================

sub is_mountpoint ($) {
  use Cwd 'abs_path';
  my ($dir, @stat1, @stat2);
# my ($dirlink);
  $dir = shift(@_);
  if (! -e $dir) {return -1;}
# $dirlink = readlink($dir);				#this does not if link is relative because it will return relative path and stat will fail
# if (defined $dirlink) {$dir = $dirlink;}
# ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks)=stat($filename);
  $dir   = abs_path($dir);
  @stat1 = stat($dir.'/.');
  @stat2 = stat($dir.'/..');
  if ($stat2[0] != $stat1[0]) {return 1;}		#device numbers do not coincide
  else                        {return 0;}
}

#===========================================================================

#sub dec2bin ($decimal,['r[everse]','t[rim]'|$nbits]) {
sub dec2bin ($;$$) {
  ### http://www.geos.ed.ac.uk/~bmg/software/Perl%20Books/cookbook/ch02_05.htm
  my ($script, $str, $decimal, $reverse, $nbits);
  $script = (split(/::/,(caller(0))[3]))[-1];
  $decimal = shift @_;					#this must be integer
  $reverse = shift @_;					#for later decoding to bits
  $nbits   = shift @_;					#16,24,32, or 'trim'
  if (!defined $reverse) {$reverse = 'n';}
  if (!defined $nbits) {
     if ($reverse =~ /^r/i) {$nbits = 32;}		#normally 24
     else                   {$nbits = 'trim';} 
  }
  if (! &is_int($decimal)) {
     print STDOUT '!!! '.$script.': non-integer number='.$decimal."\n";
  }
# $str = sprintf("%032b", $decimal);			#sprintf and unpack do same
  $str = unpack("B32", pack("N", $decimal));		#sprintf and unpack do same
  if ($nbits =~ /^t/i) {
     $str =~ s/^0+(?=\d)//;				#trim leading zeros
  } else {
     if (&is_int($nbits)) {$str = substr($str,length($str)-$nbits,$nbits);}
     else {print STDOUT '*** '.$script.': unexpected nbits='.$nbits."\n";}
  }
  ### Reversing is required if we want to decode bits: normally when we  
  ### print a number as a binary, higher bits are on the left, but for 
  ### decoding we want lower bits to be on the left (the first characters
  ### of the string):
  if ($reverse =~ /^r/i) {$str = reverse $str;}		#reverse the order of bits
  return $str;
}

#===========================================================================

#sub bin2dec ($binstr) {
sub bin2dec ($) {
  my ($binstr, $decimal);
  $binstr = shift @_;
  ### These two methods convert to unsigned int32:
  ### (they seem to do the same thing)
  ### http://www.geos.ed.ac.uk/~bmg/software/Perl%20Books/cookbook/ch02_05.htm
  $decimal = unpack("N", pack("B32", substr("0" x 32 . $binstr, -32)));
# $decimal = oct('0b' . shift);
  ### Convert from unsigned int32 to signed int32
  $decimal = ($decimal & 0x80000000) ? -((~$decimal & 0xffffffff) + 1) : $decimal; 
  return $decimal;
}

#===========================================================================

sub uint16_int16 ($) {
  ### Convert from unsigned int16 to signed int16
  my $v = shift;
  return ($v & 0x8000) ? -((~$v & 0xffff) + 1) : $v; 
}

#===========================================================================

sub uint32_int32 ($) {
  ### Convert from unsigned int32 to signed int32
  my $v = shift;
  return ($v & 0x80000000) ? -((~$v & 0xffffffff) + 1) : $v; 
}

#=================================================================

# sub read_cmdline_flags($prg,\@flags) {		#fills @flags, returns status
sub read_cmdline_flags($\@) {
  my ($prg, $fRef, $buf, $i);

  $prg  = shift(@_);
  $fRef = shift(@_);
  @$fRef = ();
  if (!defined &Which($prg)) {return 1;}
  $buf = `${prg} -h`;
  if ($buf =~ /\Q*** \E/ || length($buf) == 0) {return 1;}
  $buf =~ s/\n/ /g;
  $buf =~ s/\s+/ /g;
  $buf =~ s/Syntax: //i;
  $buf =~ s/\s*${prg}\s*//i;
# print STDOUT $buf,"\n";
  @$fRef = split(/\s+/,$buf);
  for ($i=0; $i<@$fRef; $i++) {
     @$fRef[$i] =~ s/^\[//;
     @$fRef[$i] =~ s/\]$//;
     if (@$fRef[$i] =~ /help/) {splice(@$fRef, $i, 1); $i--;}
     if (@$fRef[$i] =~ /\?/)   {splice(@$fRef, $i, 1); $i--;}
  } 
  return 0;
}

#=================================================================

# $xcoord = &pt2coord($xpt,\@xarray);
#sub pt2coord($xpt,\@xarray) {                          #returns $xcoord
sub pt2coord($\@) {                                     #returns $xcoord
# This subroutine converts non-integer (interpolated) point of
# array @xarray into respective coordinate.
  my ($xpt, $xarr_ref, $pmax, $pt0, $pt1, $dx, $xcoord);
  $xpt      = shift(@_);
  $xarr_ref = shift(@_);
 
  $pmax = @$xarr_ref - 1;

  $pt0 = int($xpt);
  if    ($pt0 < 0)     {$pt0 = 0;}
  elsif ($pt0 > $pmax) {$pt0 = $pmax;}
 
  $dx = $xpt - $pt0;
 
  if ($dx >= 0) {$pt1 = $pt0+1;}
  else          {$pt1 = $pt0-1;}
  if    ($pt1 < 0)     {$pt1 = 0;}
  elsif ($pt1 > $pmax) {$pt1 = $pmax;}
  
  $dx = abs($dx);
  $xcoord = @$xarr_ref[$pt0]*(1-$dx) + @$xarr_ref[$pt1]*($dx);
  return $xcoord;
}

#===========================================================================

sub gethost(;$) {					# returns full or short hostname
  my ($host, $shortopt);
# my (@parts);

  $shortopt = shift(@_);
  if (!defined $shortopt)     {$shortopt = 0;}
  if ($shortopt =~ /^short/i) {$shortopt = 1;}
  if ($shortopt ne '1')       {$shortopt = 0;}
  $host = $ENV{HOSTNAME};
  if (!defined $host) {$host=`hostname`; $host=~s/\n//g;}
  if (!defined $host) {$host='unknown';}
  $host = lc($host);
  if ($shortopt) {
     $host  =~ s/\..*//;					#remove everything after first "."
#    @parts = split(/\./,$host);
#    $host  = $parts[0];
  }
  return $host;
}

#==========================================================================

#sub sort_hash_array(\@hasharray,$key) {
sub sort_hash_array(\@$) {
  my ($script, $arref, $key, @sorted);
  $script = (split(/::/,(caller(0))[3]))[-1];
  $arref  = shift(@_);
  $key    = shift(@_);
  if (defined $$arref[0]->{$key}) {
     @sorted = sort { $a->{$key} <=> $b->{$key} } @$arref;
     @$arref = @sorted;
  } else {
     print STDOUT '*** '.$script.': No hash for key='.$key.' in the array. Skipping sorting.'."\n";
  }
}

#==========================================================================

# sub read_entire_file($file,\@lines,[$lock]) {
sub read_entire_file($\@;$) {
# use Fcntl;					                #not used (defines O_RDONLY constant (O_RDONLY=0,O_WRONLY=1,O_RDWR=2))
  my ($script, $file, $linesref, $lock, $buffer, $status);

  $script = (split(/::/,(caller(0))[3]))[-1];
  $file     = shift(@_);
  $linesref = shift(@_);
  $lock     = shift(@_);	if (!defined $lock) {$lock=0;}

  @$linesref = ();
  $status = sysopen(FH, $file, 0);		                #0 is used instead of O_RDONLY defined in "use Fcntl";
  if (!defined $status) {print STDOUT '*** '.$script.': Could not open '.$file."\n"; return 1;}
  if ($lock) {flock(FH, 2);}
  sysread(FH, $buffer, -s FH);
  close(FH);
  @$linesref = split(/\r?\n/,$buffer);
  return 0;
}

#==========================================================================

sub call_trace(;$$) {
  my ($i, $j, $LUN, $OUT, $msg, $scr);

  $LUN = shift(@_);		if (!defined $LUN) {$LUN = \*STDOUT;}
  $OUT = shift(@_);		if (!defined $OUT) {$OUT = 1;}

  if ($OUT =~ /no/i || $OUT eq '0') {$OUT=0;} else {$OUT=1;}

  $scr = &makeScriptName();
  $msg = 'Function calling trace:'."\n";
  $i = 0;
  do {
     $i++;
     $j = ( caller($i) )[3];
     if (defined $j) {$msg .= '  - '.$j."\n";}
     else            {$msg .= '  - '.$scr."\n";}
  } while (defined $j);
  if ($OUT) {print $LUN $msg;}
  return $msg;
}

#=================================================================

sub get_login_environment {
  my ($shell, $env, @lines, $line, $key, $val);
  if ($^O =~ /^Win/i) {return 1;}				#this function does not work on Windows (no getpwuid)
  $shell = shift || (getpwuid($<))[8];				#getpwuid does not exist on Windows
  $env = `$shell -lc env`;
  @lines = split(/\n/,$env);
  foreach $line (@lines) {
     if ($line =~ /=/) {
        ($key,$val) = split("=",$line,2);
        if ($key !~ /BASH_FUNC_module/ && length($key)>0) {$ENV{$key} = $val;}
     }
  }
  return 0;
}
  
#=================================================================

#sub ziplist($zipfile,\%hash) {			#Fills %hash as key=filename, value=length
sub ziplist($\%) {

  my ($scr, $zipfile, $hashref, $prg, @lines, $line, @words);

  $scr = (split(/::/,(caller(0))[3]))[-1];
  $zipfile = shift(@_);
  $hashref = shift(@_);

  if (! -e $zipfile) {print STDERR '*** '.$scr.': zipfile ['.$zipfile.'] does not exist!'."\n"; return -1;}
  $prg = &Which('unzip');
  if (! defined $prg) { print STDERR '*** '.$scr.': no unzip program found!'."\n"; return -1;}
  
  @lines = split(/\n/,`$prg -l $zipfile 2>/dev/null`);
  foreach $line (@lines) {
     if ($line =~ /^Archive/ || $line =~ /^\s*Length/ || $line =~ /^--/ || $line =~ /files?$/) {next;}
     $line =~ s/^\s*//;
     @words = split(/\s+/,$line);
     if (@words != 4) {next;}
     $$hashref{lc($words[3])} = $words[0];
  }
  return 0;
}
  
#=================================================================

#sub zipextract($zipfile,$file,\@lines) {			#Fills @lines with file content
sub zipextract($\%) {

  my ($scr, $zipfile, $file, $linesref, $prg, %hash);

  $scr = (split(/::/,(caller(0))[3]))[-1];
  $zipfile  = shift(@_);
  $file     = shift(@_);
  $linesref = shift(@_);

  @$linesref = ();
  if (! -e $zipfile) {print STDERR '*** '.$scr.': zipfile ['.$zipfile.'] does not exist!'."\n"; return -1;}
  $prg = &Which('unzip');
  if (! defined $prg) {print STDERR '*** '.$scr.': no unzip program found!'."\n"; return -1;}

  if (&ziplist($zipfile,\%hash)) {return -1;}
  if (!defined $hash{$file}) {print STDERR '*** '.$scr.': no ['.$file.'] in zip ['.$zipfile.']!'."\n"; return -1;}
# if ($hash{$file} eq 0) {print STDERR '*** '.$scr.': ['.$file.'] in zip ['.$zipfile.'] is empty!'."\n"; return -1;}
  
  ### -p=extract files to pipe, no messages, -C=match filenames case-insensitively
  @$linesref = split(/\n/,`$prg -pC $zipfile $file 2>/dev/null`);
  return 0;
}

#=================================================================

#sub import_elements($file,\@elements) {			#fills elements data and returns status
sub import_elements($\@) {

# This script imports elements data from elements.dat into a hash array.
# The keywords in the array are ('name', 'code', 'z', 'a').

  my ($script, $elref, $file, @lines, @words);
  my ($status, $str, $nwords, $msg);
        
  $script = &sub_name();
  $file   = shift(@_);
  $elref  = shift(@_);
  @$elref = ();

  if (! -e $file) {
     $msg = '*** '.$script.':  Cannot find data file ['.$file.'].';
     print STDOUT $msg."\n";
     return -1;
  }
  print STDOUT $script.': Opening elements file: '.$file."\n";
  $status = &read_entire_file($file,\@lines);
  if ($status) {
     $msg = '*** '.$script.': Error reading data file ['.$file.'].';
     print STDOUT $msg."\n";
     return -1;
  }
  foreach $str (@lines) {
     if ($str =~ /^\#/) {next;}
     $str =~ s/\s+$//;					#remove end-of-line spaces
     $str =~ s/^\s+//;					#remove beg-of-line spaces
     if (length($str) == 0) {next;}
     @words = split(/\s+/,$str);
     $nwords = @words;
     if ($nwords != 4) {
        $msg = '*** '.$script.': incorrect record ['.$str.'] in data file ['.$file.'] -- nWords='.$nwords;
        print STDOUT $msg."\n";
        return -1;
     }
     if (! &is_int($words[2])) {
        $msg = '*** '.$script.': incorrect record ['.$str.'] in data file ['.$file.'] -- non-int Z='.$words[2];
        print STDOUT $msg."\n";
        return -1;
     }
     $words[3] =~ s/[\(\)]//g;
     if (! &is_float($words[3])) {
        $msg = '*** '.$script.': incorrect record ['.$str.'] in data file ['.$file.'] -- non-float a='.$words[3];
        print STDOUT $msg."\n";
        return -1;
     }
     push @$elref, {name=>$words[0], code=>$words[1], z=>$words[2], a=>$words[3]};
  }
  @$elref = sort{$a->{'z'} <=> $b->{'z'}} @$elref;	#sort over Z
  return 0;
}

#=================================================================

sub getLatestUserAgent() {

  my ($ver, $ua);
### Check at https://www.whatismybrowser.com/detect/what-is-my-user-agent/
  $ver = '150.0';
  $ua = sprintf 'Mozilla/5.0 (X11; Linux x86_64; rv:%s) Gecko/20100101 Firefox/%s',$ver,$ver;
  return $ua;
}

