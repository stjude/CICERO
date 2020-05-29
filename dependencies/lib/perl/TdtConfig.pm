#!/usr/bin/env perl
# Helper functions for using config files in shell script format
package TdtConfig;

use strict;
use warnings;
use Carp;
use File::Path qw/ make_path /;

# Helper function for reading config
# Parameters:
# 1. Category (ignored if option 2 is a file)
# 2. Config name or file
# Returns hash ref of var => val pairs
sub readConfig {
  my $cat = shift; 
  my $conf = shift; 

  my $file = &findConfig($cat, $conf);
  chomp $file;
  confess "Config ($cat, $conf) not found" unless(-f $file);
  open(IN, $file) or die "Could not open $file";
  my $out = {};
  while(my $line = <IN>) {
    chomp $line;
    next unless($line);
    next if(substr($line, 0, 1) eq '#');
    my($var, $val) = split(/\t/, $line, 2);
    $out->{$var} = $val if($var);
  }
  close(IN);
  return $out;
}

# Helper function for writing config in sh format
# - $file: config file to read (in sh format)
# - $config: configuration as a hash ref of var => val pairs
# - $desc: optional description to write to the file
sub writeConfig {
  my($file, $config, $desc) = @_;
  open(OUT, ">$file") or die "Could not open $file";
  print OUT "##SIP $desc\n" if($desc);
  while(my ($var, $val) = each (%$config) ) {
    print OUT "$var\t$val\n";
  }
  close(OUT);
}

# This searches for a file named <conf>.config.txt in the category configuration
# directory under the root configuration directory.
# The root configuration directory used will be:
# The directory in the SJ_CONFIGS environment variable
sub findConfig {
	my $cat = shift; 
	my $conf = shift; 

	if($conf && substr($conf, 0, 1) eq "/"){
		if (-f "$conf"){
			my $abs_path = `readlink -f $conf`; 
			return $abs_path; 
		}
		else{
			die "No config file at given absolute path: $conf"; 
		}
	}
	my $sj = $ENV{"SJ_CONFIGS"}; 
	if (!$sj){
		confess "Configuration directory not set"; 	
	}
	my $dir = "$sj/$cat";
	if(! $conf){
		return $dir; 
	}
	my $file = "$dir/$conf.config.txt";
	if (-f "$file"){
		$file = `readlink -f $file`; 
		chomp $file;
		return $file; 
	}  
	else{
		my $alias = "$dir/aliases.txt";
		my $new_conf; 
		if (-f "$alias"){
			$new_conf = `awk -v alias=$conf '\$1 == alias { print \$2; exit }' $alias`;
			chomp $new_conf;
		}
		if ($new_conf){
			$file = "$dir/$new_conf.config.txt"; 
			if (-f "$file"){
				$file = `readlink -f $file`; 
				chomp $file; 
				return $file; 
			}
		} 
	}
	
	return "";  
}

1;
