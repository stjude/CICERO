package Counter;

use MethodMaker qw(
		   total
		   count
		   delta
		   mod
		   datestamp
                   start_time
                   ticks
		  );
use strict;
use Configurable;

#use misc;
use MiscUtils qw(average seconds_to_time);

@Counter::ISA = qw(Configurable Exporter);
@Counter::EXPORT_OK = qw(generate_counter);

use constant TICKS_TO_TRACK => 10;

sub new {
  my ($type, $all, %options) = @_;
  my $self = {};
  bless $self, $type;

  $self->ticks([]);
  $self->start_time(time);
  if ($all) {
    if (ref $all eq "ARRAY") {
      $self->total(scalar @{$all});
    } elsif (ref $all eq "HASH") {
      $self->total(scalar keys %{$all});
    } elsif (ref $all) {
      die "huh?";
    } else {
      $self->total($all);
    }
  }
  $self->count(0);
  $self->configure(%options);

  return $self;
}

sub next {
  my ($self, $thing, %options) = @_;

  my $count = $self->count();
  if (my $increment = $options{"-increment"}) {
    $count += $increment;
    $self->count($count);
  } else {
    $self->count(++$count);
  }

  my $mod = $self->mod();
  return if $mod and $count % $mod != 0;

  my $total = $self->total();
  my $delta = $self->delta();

  my $frag = "";
  my $percent;
  if ($total) {
    $percent = $count * 100 / $total;
    $frag = sprintf '/%d (%.1f%%)', $total, $percent;
    if ($delta) {
      my $chunk = int($total * ($delta / 100));
#      print "$count: $delta $chunk\n";
      return unless $count % $chunk == 0;
    }
  } else {
    die "delta doesn't work without total" if $delta;
  }

  if ($options{"-percent"}) {
    return $percent;
  } else {
    printf STDERR "%s: ", scalar localtime if $self->datestamp;

    printf STDERR "Processing %s \#%d%s",
      ($thing || ""),
	$count,
	  $frag;

    if (1 and defined($percent)) {
      my $elapsed_secs = time - $self->start_time();

      my $ticks = $self->ticks;
      my $tick = {
		  time => time(),
		  percent => $percent,
		 };
      shift @{$ticks} if scalar @{$ticks} >= TICKS_TO_TRACK;
      push @{$ticks}, $tick;

      my $total_time;
      if (@{$ticks} > 1) {
	my @estimated_total_times;
	for (my $i = 1; $i < scalar @{$ticks}; $i++) {
	  my $diff_time = $ticks->[$i]->{time} - $ticks->[$i - 1]->{time};
	  # raw elapsed time from last tick to this one
	  my $diff_pct = $ticks->[$i]->{percent} - $ticks->[$i - 1]->{percent};
	  # total percent progress made between last tick and this one

  	  # diff_elapsed:2 diff_pct:3.33333333333333 

          # ($diff_pct / 100) = $diff_time / $total
	  # $diff_pct * $total = 100 * $diff_time
	  # total = (100 * $diff_time) / $diff_pct
	  
	  # 3.3333 / 100 = 2 / $total
	  # 3.3333 * $total = 100 * 2
	  # total = (100 * 2) / 3.333
	  # total = 60

	  my $est_total = ($diff_time * 100) / $diff_pct;
	  # the amount of time the total would take
	  # based ONLY on most recent progress
	  push @estimated_total_times, $est_total;

#	  printf STDERR "diff_elapsed:%d diff_pct:%s est_total:%s\n", $diff_time, $diff_pct, $est_total;
	}

	my $est_total_based_on_recent = average(\@estimated_total_times);
	my $est_total_based_on_progress = $elapsed_secs * (100 / $percent);

	foreach ($est_total_based_on_recent,
		 $est_total_based_on_progress) {
	  $_ = $elapsed_secs if $_ < $elapsed_secs;
	}

	my $weight_to_progress = 0.80;

	$total_time = ($est_total_based_on_progress * $weight_to_progress) +
	  ($est_total_based_on_recent * (1 - $weight_to_progress));

#	printf STDERR "total_est_recent:%d total_est_progress:%d weighted:%d\n", $est_total_based_on_recent, $est_total_based_on_progress, $total_time;

      } else {
	$total_time = $elapsed_secs * (100 / $percent);
      }

      my $remaining = $total_time - $elapsed_secs;

      printf STDERR " elapsed:%s est_total:%s remaining:%s",
      seconds_to_time($elapsed_secs, "-string" => 1) || "0",
      seconds_to_time($total_time, "-string" => 1) || "??",
      seconds_to_time($remaining, "-string" => 1) || "??";
    }
    print STDERR "\n";
  }
}

sub generate_counter {
  my (%options) = @_;
  my $count = 0;
  my $total = $options{"-total"} || die;

  my $coderef;

  my $code = '$coderef = sub {
$count++;
';

  if (my $mod = $options{"-mod"}) {
    $code .= sprintf 'return undef unless $count %% %d == 0;' . "\n", $mod;
  }

  if (my $percent = $options{"-percent"}) {
    # return percentage of total progress
    $code .= 'return $count * 100 / $total;';
  } else {
    # just return raw count
    $code .= 'return $count;';
  }

  $code .= "\n};";
  eval $code;
  return $coderef;

  print $code;

}

1;
