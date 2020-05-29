package MiscUtils;
# mne 9/2010

use strict;
use Exporter;

use Carp qw(cluck confess carp);
use FileUtils qw(read_simple_file);
use List::Util qw(min max sum);

@MiscUtils::ISA = qw(Exporter);
@MiscUtils::EXPORT_OK = qw(
                           load_reorder
average
median
decile
unique_list
unique_ordered_list
unique_list_hashes
dump_die
log_message
split_list
trim_flanking_whitespace
seconds_to_time
longest
build_argv_list
unquote
get_hash_option
get_core_count
shell_cmd
interval_overlap
float_trim
get_java_version
module_dependencies
mean_variance_standard_deviation
                          );


sub load_reorder {
  my (%options) = @_;
  my $jobs = $options{"-jobs"} || die "-jobs";
  my $load_regexp = $options{"-bucket-regexp"} || die "-bucket-regexp";
  my $VERBOSE = 0;

  #
  # create a job object containing each "bucket"
  #
  my %pending_jobs;
  my $job_id = 0;
  foreach my $job_name (@{$jobs}) {
    my %job;
    $job{name} = $job_name;
    $job_name =~ /$load_regexp/ || die;
    $job{buckets} = [ $1 ];
    my $id = ++$job_id;
    $job{id} = $id;
    $pending_jobs{$id} = \%job;
  }

  my @sorted;
  my %load;
  while (%pending_jobs) {
    my $choice;

    printf STDERR "current load: %s\n", join " ", map {$_ . "=" . $load{$_}} sort {$a <=> $b} keys %load if $VERBOSE;
    foreach my $id (keys %pending_jobs) {
      my $job = $pending_jobs{$id};
      my $buckets = $job->{buckets};
      my $any_load;
      foreach (@{$buckets}) {
	if ($load{$_}) {
	  $any_load = 1;
	  last;
	}
      }

      if (!$any_load) {
	$choice = $id;
	last;
      }
    }

    unless ($choice) {
      #
      # no job hits only unloaded buckets.
      # try to pick a job which LEAST overlaps previous jobs.
      #

      my $best_load;
      my $best_id;

      foreach my $id (keys %pending_jobs) {
	my $job = $pending_jobs{$id};
	my $buckets = $job->{buckets};

	my %ub = map {$_, 1} @{$buckets};
	my $current_load = 0;
	foreach (keys %ub) {
	  $current_load += $load{$_} || 0;
	}
	foreach (@{$buckets}) {
	  $current_load++;
	}

	if (not(defined($best_load)) or $current_load < $best_load) {
	  printf STDERR "new_load:%d\n", $current_load if $VERBOSE;
	  $best_load = $current_load;
	  $best_id = $id;
	}
      }

      $choice = $best_id;
    }

    die unless $choice;

    my $job = $pending_jobs{$choice};
    foreach (@{$job->{buckets}}) {
      print STDERR "add load to $_\n" if $VERBOSE;
      $load{$_}++;
    }
    push @sorted, $job;
    delete $pending_jobs{$choice};
  }

  my @names = map {$_->{name}} @sorted;

  return \@names;
}

sub average {
  # return the average (mean) value of a list
  # argument: ref to list

  #    eval(join "+",@_) / scalar @_;
  # spiffy, but expensive
  my $total = 0;
  if (ref $_[0]) {
    cluck "division by zero" unless @{$_[0]};
    foreach (@{$_[0]}) {
      $total += $_;
    }
    return ($total / scalar @{$_[0]});
  } else {
    foreach (@_) {
      $total += $_;
    }
    cluck "division by zero" unless @_;
    return ($total / scalar @_);
  }
}

sub median {
  # return the median value of a list
  if (ref $_[0] eq "ARRAY") {
    return (sort {$a <=> $b} @{$_[0]})[scalar @{$_[0]} / 2];
  } else {
    return (sort {$a <=> $b} @_)[scalar @_ / 2];
  }
}

sub unique_ordered_list {
  # remove duplicates from array while preserving sort order
  my ($array_in, %options) = @_;
  my @out;
  my %saw;
  foreach my $thing (@{$array_in}) {
    next if $saw{$thing};
    $saw{$thing} = 1;
    push @out, $thing;
  }
#  return wantarray ? @out : \@out;
  # causes problems if caller is a list context, e.g. if creating a hash
  return \@out;
}

sub unique_list {
  # remove duplicates from array while preserving sort order
  my ($array_in, %options) = @_;
  my %v = map {$_, 1} @{$array_in};
  my @results;
  if ($options{"-sort"}) {
    @results = sort keys %v;
  } else {
    @results = keys %v;
  }
  return \@results;
}

sub dump_die {
  # dump a hash row and exit or warn
  my ($row, $message, $warn) = @_;
  $message = "" unless $message;
  my $label = $warn ? "WARNING" : "ERROR";
  if ($warn) {
    my $old = $Carp::Verbose;
    $Carp::Verbose = 1;
    cluck sprintf "%s: %s\n", $label, $message;
    $Carp::Verbose = $old;
  } else {
    printf STDERR "%s: %s\n", $label, $message;
  }
  foreach (sort keys %{$row}) {
    my $thing = $row->{$_};
    my $v;
    if (!defined($thing)) {
      $v = '[undef]';
    } elsif (ref $thing eq "ARRAY") {
      $v = sprintf '[ARRAYREF] %s', join ",", @{$thing};
    } elsif (ref $thing eq "HASH") {
      $v = "[HASHREF]";
    } elsif (ref $thing eq "SCALAR") {
      $v = "[SCALAR] " . $$thing;
    } elsif (ref $thing) {
      $v = sprintf "[ref to %s]", ref $thing;
    } else {
      $v = $thing;
    }
    printf STDERR "  %s: %s\n", $_, $v;
  }
  confess unless $warn;
}

sub log_message {
  my ($msg, %options) = @_;
  chomp $msg;
  printf STDERR "%s: %s\n", scalar(localtime), $msg;
}

sub split_list {
  # returns list of references to smaller lists of specified size.
  # args:
  #  - reference to a long list
  #  - maximum number of elements in each returned sublist
  # 
  # Returns list of references to smaller lists.
  # 9/2004: work better with fractional values
  my ($big_list_ref, $step_size) = @_;
  my $total_size = scalar @{$big_list_ref};

  my $ptr;
  my $end;
  my $max = $total_size - 1;
  my @results;
  my @slice;
  my $last_start;
  my $start = 0;
  for ($ptr = 0; $ptr < $total_size; $ptr += $step_size) {
    # $ptr is the virtual index, very possibly fractional.
    $start = int($ptr);
    # beginning slice array index in integer space
    $start = $last_start + 1 if defined($last_start) and $start <= $last_start;
    # ensure we're always moving forward through the array and
    # never repeat an index.  This can happen if we're passed a fractional value
    # that's too small.
    last if $start >= $total_size;
    # after adjustments, might be past the end already
    $end = int($ptr + $step_size - 1);
    # ending slice array index in integer space
    $end = $start if $end < $start;
    # ensure slice is at least one element long
    $end = $max if $end > $max;
    # end bounds checking
    @slice = @{$big_list_ref}[$start .. $end];
#    printf STDERR "DEBUG: start=%s end=%s: %s\n", $start, $end, join ",", @slice;
    push @results, [ @slice ];
    $last_start = $start if not(defined($last_start)) or $start > $last_start;
    # record last start position
  }
  
  return wantarray ? @results : \@results;
}

sub trim_flanking_whitespace {
  my ($thing) = @_;
  $thing =~ s/^\s+//;
  $thing =~ s/\s+$//;
  return $thing;
}

sub seconds_to_time {
  # convert a time in seconds to days, hours, minutes, seconds
  my ($secs, %options) = @_;

  my $second = 1;
  my $minute = $second * 60;
  my $hour = $minute * 60;
  my $day = $hour * 24;
  
  my $days = 0;
  my $hours = 0;
  my $minutes = 0;
  my $seconds = 0;

  my $count;

  my @set = (
	     [$day, \$days, "day"],
	     [$hour, \$hours, "hour"],
	     [$minute, \$minutes, "minute"],
	     [$second, \$seconds, "second"]
	    );
  
  my ($unit, $ref, $label);

  foreach (@set) {
    ($unit, $ref) = @{$_};
    $$ref = int($secs / $unit);
    if ($$ref) {
      $secs = $secs % $unit;
    }
  }

  if ($options{"-string"}) {
    my @things;
    my $long = $options{"-long"};
    foreach (@set) {
      ($unit, $ref, $label) = @{$_};
      if ($$ref) {
	if ($long) {
	  push @things, sprintf '%d %s%s', $$ref, $label, $$ref == 1 ? "" : "s";
	} else {
	  push @things, sprintf '%d%s', $$ref, substr($label,0,1);
	}
      }
    }
    return join ",", @things;
  } else {
    return ($days, $hours, $minutes, $seconds);
  }
}

sub longest {
  my ($thing) = @_;
  my $longest = 0;
  if (ref $thing eq "ARRAY") {
    foreach my $s (@{$thing}) {
      my $l = length($s);
      if ($l > $longest) {
	$longest = $l;
      }
    }
  } else {
    die;
  }
  return $longest;
}

sub unique_list_hashes {
  my ($rows, $key, %options) = @_;
  my %v = map {$_->{$key}, 1} @{$rows};
  return [sort keys %v];
}

sub build_argv_list {
  # e.g.
  # build_argv_list("-flags" => \%FLAGS, "-single" => "file", "-set" => "files");
  my (%options) = @_;
  my $flags = $options{"-flags"} || die "-flags";
  my $f_single = $options{"-single"} || die;
  my $f_set = $options{"-set"} || die;
  my $f_glob = $options{"-glob"};

  my @set;
  push @set, $flags->{$f_single} if $flags->{$f_single};
  if ($f_glob and my $p = $flags->{$f_glob}) {
    my @hits = glob($p);
    if (@hits) {
      push @set, @hits;
    } else {
      confess "no match for $p";
    }
  }


  if (my $fn = $flags->{$f_set}) {
    my $list = read_simple_file($fn);
    push @set, @{$list};
  }
  confess sprintf 'specify -%s [file] | -%s [listfile] | glob [pattern]', $f_single, $f_set unless @set;

  return \@set;
}

sub unquote {
  my ($string) = @_;
  $string =~ s/^([\"\'])(.*)\1$/$2/;
  return $string;
}

sub get_hash_option {
  my ($hash, @params) = @_;
  my $required = 1;
  my @results;
  foreach my $param (@params) {
    dump_die($hash, "can't find $param") if $required and not(exists $hash->{$param});
    push @results, $hash->{$param};
  }
  return wantarray ? @results : $results[0];
}

sub get_core_count {
  my $count = 0;
  my $verbose = 0;

  printf STDERR "cpuinfo debug start:\n" if $verbose;
  if (open(TMPCC, "/proc/cpuinfo")) {
    while (<TMPCC>) {
      $count++ if /^processor/;
      print STDERR $_ if $verbose;
    }
    close TMPCC;
  }
  printf STDERR "cpuinfo debug end\n" if $verbose;
  return $count;
}

sub shell_cmd {
  my (%options) = @_;
  my $cmd = $options{"-cmd"} || confess "-cmd";
  my $log = $options{"-log"};
  my $return_exit = $options{"-return-exit-code"};
  $log = 1 unless defined $log;
  log_message($cmd) if $log;
  my $result;
  if ($options{"-backtick"}) {
    $result = `$cmd`;
    die "$cmd failed with $?" if $?;
    chomp $result;
  } else {
    system $cmd;
    my $code = $?;
    if ($return_exit) {
      $result = $code;
    } elsif ($code) {
      die "$cmd failed with $code" if $code;
    } else {
      $result = 1;
    }
  }
  return $result;
}

sub interval_overlap {
  my ($s1, $e1, $s2, $e2) = @_;
  my $overlap = 0;
  if ($s1 > $e2 or $e1 < $s2) {
    # no overlap: 
    # - 1st interval starts after 2nd interval ends
    # - 1st interval ends before 2nd interval starts

#  } elsif (0 and $s1 == $s2 and $e1 == $e2) {
    # identical
#    $overlap = $e1 - $s1 + 1;
#  } elsif (0 and $s1 >= $s2 and $e1 <= $e2 and $s1 <= $s2 and $e1 <= $e2) {
    # interval 1 contained within interval 2
#    print STDERR "1 in 2\n";
#    $overlap = $e1 - $s1 + 1;
#  } elsif (0 and $s2 >= $s1 and $e2 <= $e1 and $s2 <= $s1 and $e2 <= $e1) {
    # interval 2 contained within interval 1
#    print STDERR "2 in 1\n";
#    $overlap = $e2 - $s2 + 1;
  } else {
    # overlap
#    die "$s1 $e1 $s2 $e2";
    my $trim_start = max($s1, $s2);
    my $trim_end = min($e1, $e2);
    $overlap = $trim_end - $trim_start + 1;
  }
  return $overlap;
}

sub decile {
  my ($array, $decile) = @_;
  die "first argument must be arrayref" unless ref $_[0] eq "ARRAY";
  die "decile must be fractional value > 0 and <= 1" unless $decile > 0 and $decile <= 1;
  my @s = sort {$a <=> $b} @{$_[0]};
  my $idx = int(@s * $decile);
  $idx = @s - 1 if $idx >= @s;
  return $s[$idx];
}

sub float_trim {
  my ($v_in, %options) = @_;

  my $sprintf;
  my $precision = $options{"-precision"};
  if ($precision) {
    # optionally trim floating point precision to X places; DESTRUCTIVE!
    $sprintf = '%.' . $precision . 'f';
  }
  my $v_out = $v_in;

  if (/^([\-\+]?)(\d+\.\d+)$/) {
    my ($sign, $float) = ($1, $2);

    if ($precision) {
      $float =~ /\.(\d+)$/ || die;
      if (length $1 > $precision) {
	$float = sprintf $sprintf, $float;
      }
    }

    if ($float =~ s/^0\./\./) {
      # e.g. 0.38 => .38
#      die "$v_in $float";
    }

    if ($float =~ s/0+$//) {
      # e.g. 2.078000 => 2.078
#      print STDERR "tail $v_in $float\n";
    }

    if ($float eq ".") {
      # 0.0 => . => 0
      $float = 0;
    }

    $v_out = $sign . $float;

#    printf STDERR "in:$v_in sign:$sign value:$float out:$v_out\n";

  }

  return $v_out;
}

sub get_java_version {
  my (%options) = @_;
  my $info = `java -version 2>&1`;
  my $minimum = $options{"-require-minimum"};
  my $exact = $options{"-require-exact"};

  my $result;
  if ($info =~ /java version "((\d+\.\d+)\S+)"/) {
    my ($full, $major) = ($1, $2);
    $result = $major;

    if ($minimum) {
      die sprintf "minimum java version %s required, your loaded version is %s",
	$minimum, $major if $major < $minimum;
    } elsif ($exact) {
      die sprintf "exact java version %s required, your loaded version is %s",
	$exact, $major unless $major == $exact;
    }
  } else {
    die "can't detect java version from $info";
  }
  return $result;
}

sub module_dependencies {

  my ($module, $path);

  my %src;

  while (($module, $path) = each %INC) {
    # HACKY: local system policies rather than a global policy
    if ($module =~ /^\// or
	$path =~ /^\/hpcf\/apps\/perl\/install\//) {
      # ignore:
      # - module name starts with a slash (?)
      # - system/builtin module, ignore
    } elsif ($path =~ /cluster_code\/trunk\/([^\/]+)\//) {
      # local SJ compbio package
      $src{$1}{$module} = 1;
    } else {
      $src{third_party}{$module} = 1;
    }
  }

  printf STDERR "detected module dependencies:\n";
  foreach my $class (sort keys %src) {
    printf STDERR "  %s:\n", $class;
    foreach my $m (sort keys %{$src{$class}}) {
      printf STDERR "    %s\n", $m;
    }
  }
}

sub mean_variance_standard_deviation {
  # given a [reference to a] list, return mean, variance and standard deviation
  my $listref = shift;

  my $mean = average($listref);
#  my @squared_variances = map {abs($mean - $_) ** 2} @{$listref};
  my @squared_variances = map {($mean - $_) ** 2} @{$listref};
  # abs _required_: perl's exponentiation is weird w/negative nonvariables;
  # ie "print -2 ** 2" returns -4 (!).
  # Oops; this is an operator precedence problem:
  # "print (-2) ** 2" works correctly.  "Never mind."

  my $variance = sum(@squared_variances) / @{$listref};
  # average of the sum of squared differences from the mean (whew)
  my $standard_deviation = sqrt($variance);

  return ($mean, $variance, $standard_deviation);
}


1;
