package Cluster;
# cluster job submission wrapper
# MNE 5/2012
#
# TO DO: memory configuration, etc.

use strict;

use File::Basename;

use Configurable;
use ClusterLogFile;

use Carp qw(confess);
#use DevelopmentPath;

@Cluster::ISA = qw(Configurable Exporter);
@Cluster::EXPORT_OK = qw(run_now_or_later wait_for_jobs run_jobs);

@Cluster::ARGV_RAW = @main::ARGV;
# copy of raw ARGV before any parsing/modification

use MethodMaker qw(
		    command
		    outfile

		    app
		    project
		    save_output
                    misc

memory_reserve_mb
memory_limit_mb
memory

node_class
tracking_dir
force

create_bsub_file
create_job_file

debug
queue

core_limit

processors
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->save_output(1);
#  $self->node_class("idataplex");
  # 3/2013: disable in an attempt to avoid very long queue times

  $self->create_bsub_file(1);
  $self->create_job_file(1);
  $self->core_limit(0);
  # by default suppress creation of coredumps

  # alternatives:
  # node0: "blade"
  # lm*: "large_mem"
  # verna: "-m verna" (NOT -R)
  $self->configure(%options);
  return $self;
}

sub guess_app {
  my ($self) = @_;
  my $cmd = $self->command;
  if ($cmd and !$self->app) {
    # try to guess application profile from command line
    if ($cmd =~ /[\/ ]java/) {
#      $self->app("jdk6u30");
# 5/2019: disable esp. as Java 1.8+ now typically required
    } elsif ($cmd =~ /samtools/) {
      $self->app("samtools");
    }
  }
}

sub get_tracking_file {
  my ($self, %options) = @_;
  my $suffix = $options{"-suffix"} || die "-suffix";
  my $outfile = $self->outfile || confess "need -outfile";
  my $track_dir = $self->tracking_dir();
  my $result;
  if ($track_dir) {
    $result = sprintf '%s/%s.cluster.%s', $track_dir, basename($outfile), $suffix;
  } else {
    $result = sprintf '%s.cluster.%s', $outfile, $suffix;
  }
  return $result;
}

sub run {
  my ($self, %options) = @_;
  my $project = $self->project || die "need -project";
  my $no_outfile = $options{"-no-outfile"};
  my ($outfile, $cluster_out_file, $cluster_err_file, $job_file);

  my $needed = 1;
  my $verbose = 0;

  if ($no_outfile) {
    $cluster_out_file = basename($0) . ".cluster.out";
    $cluster_err_file = basename($0) . ".cluster.err";
  } else {
    $outfile = $self->outfile || die "need -outfile";
    $cluster_out_file = $self->get_tracking_file("-suffix" => "out");
    $cluster_err_file = $self->get_tracking_file("-suffix" => "err");
    $job_file = $self->get_tracking_file("-suffix" => "job");

    if (-s $cluster_out_file) {
      my $clf = new ClusterLogFile("-cluster" => $self);
      if (! $clf->ok()) {
	my $msg = $clf->get_error_summary();
	if (0) {
	  printf STDERR "last run had errors (%s): unlinking %s and %s\n", $msg, $outfile, $cluster_out_file;
	  unlink($outfile, $cluster_out_file);
	} else {
	  printf STDERR "last run had errors (%s): renaming %s\n", $msg, join ",", $outfile, $cluster_out_file, $cluster_err_file;
	  foreach my $obs ($outfile, $cluster_out_file, $cluster_err_file) {
	    if (-f $obs) {
	      my $new = $obs . ".crash";
	      rename($obs, $new) || die;
	    }
	    die if -f $obs;
	  }
	}
      }
    }

#  if (-s $outfile) {
    if (-e $outfile and not($self->force)) {
      # might be zero length, e.g. variant detection on low-coverage reference
#    print STDERR "already done\n";
      $needed = 0;
    }

    if ($needed and -s $job_file) {
      # job isn't finished but a job was submitted to the cluster;
      # see if it's running or dead
      open(JTMP, $job_file) || die;
      my $job_id = <JTMP>;
      chomp $job_id;
      close JTMP;

      open(JCHECK, "bjobs $job_id|") || die;
      my $first = <JCHECK>;
      chomp $first;
      my @lines;
      while (<JCHECK>) {
	chomp;
	push @lines, $_;
      }
      close JCHECK;
      if ($first =~ /JOBID/) {
	my @h = split /\s+/, $first;
	my %row;
	@row{@h} = split /\s+/, $lines[0];
	my $status = $row{STAT};

	printf STDERR "job status for %s: %s, ", $outfile, $status if $verbose;
	if ($status eq "EXIT") {
	  print STDERR "resubmitting\n" if $verbose;
	} elsif ($status eq "RUN") {
	  # job is running
	  print STDERR "not resubmitting\n" if $verbose;
	  $needed = 0;
	} elsif ($status eq "DONE") {
	  printf STDERR "resubmitting (outfile $outfile deleted??)\n" if $verbose;
	  # - user deleted outfile, so resubmit OK
	  # - outfile is not actually produced by job, which could lead
	  #   to redundant job submission  :/
	} elsif ($status eq "PEND") {
	  print STDERR "not resubmitting\n" if $verbose;
	  $needed = 0;
	} elsif ($status eq "ZOMBI") {
	  # saw this happen for some jobs where state was UNKWN,
	  # after killing them status changed to ZOMBI
	  print STDERR "resubmitting\n" if $verbose;
	} else {
	  die "unhandled status $status, fix me";
	  # pending, etc.
	}
      }
    }
  }

  if ($needed) {
    unlink($cluster_out_file, $cluster_err_file);
    # make sure log files are new (not appended)

#    $self->guess_app();
# 3/2014: disabled: not compatible w/clinical cluster

    my $command = $self->command || die "need command";
    if (0) {
      print STDERR "DEBUG CMD\n";
      $command = "sleep 10";
    }

    my $cmd = "bsub";
    $cmd .= " -P $project";
    $cmd .= sprintf " -q %s", $self->queue if $self->queue;
    $cmd .= sprintf " -C %d", $self->core_limit();

    my $processors = $self->processors() || 1;
    if ($processors > 1) {
      $cmd .= sprintf ' -n %d -R "span[hosts=1]"', $processors;
    }

    if (my $reserve = $self->memory_reserve_mb() || $self->memory()) {
	my $limit = $self->memory_limit_mb() || $self->memory() || die "need memory limit";

	if ($processors > 1) {
	  die "don't know how to divide $reserve" unless $reserve =~ /^\d+$/;
	  $reserve /= $processors;
	}

#	my $memory_config = sprintf ' -R "rusage[mem=%d]" -M %d -v %d',	$reserve, $limit, $limit;
	my $memory_config = sprintf ' -R "rusage[mem=%d]"', $reserve;
	# 10/4/2016:
	# - when requesting multiple cores, divide memory usage by core count
	# - don't use "-v" anymore

	$cmd .= $memory_config;
    }

    if (my $nc = $self->node_class) {
	$cmd .= " -R " . $nc;
    }

    if ($self->save_output) {
      $cmd .= " -o " . $cluster_out_file;
      $cmd .= " -e " . $cluster_err_file;
    } elsif (0) {
      # suppress email
      $cmd .= " -o /dev/null";
      $cmd .= " -e /dev/null";
    }

    $cmd .= " -app " . $self->app() if $self->app();
    $cmd .= " " . $self->misc() if $self->misc();

    $command = '"' . $command . '"' if $command =~ /\|/;
    $cmd .= " " . $command;

    # bsub -app jdk6u30 -P JZ_matrix -o test.out -e test.err env java -cp /home/medmonso/lib/bambino.jar Ace2.ReadReport -bam /nfs_exports/genomes/1/PCGP/BucketRaw/SJINF/SJINF013_G-TB-03-0007.bam -passthrough job_SJINF013.tab -min-mapq 0

    if ($self->create_bsub_file) {
      my $bsub_file = $self->get_tracking_file("-suffix" => "bsub");
      open(WSUB, ">" . $bsub_file) || die "can't write to $bsub_file";
      printf WSUB "#!/bin/sh\n";
      printf WSUB "%s\n", $cmd;
      close WSUB;
    }

#    die $cmd if ($self->processors || 1) > 1;

    if ($self->debug()) {
      die "DEBUG, not submitting $cmd";
    }

    open(CSUB, $cmd . " 2>&1|") || die;
    while (<CSUB>) {
      print;
      # Job <906000> is submitted to queue <pcgp>.
      if ($self->create_job_file() and /Job <(\d+)>/) {
	my $job_id = $1;
	open(JSAVE, ">" . $job_file) || die "can't write to $job_file";
	printf JSAVE "%d\n", $job_id;
	close JSAVE;
      }
    }
  }
}

sub done {
  my ($self, %options) = @_;
  my $done = 0;
  my $require_outfile = $options{"-no-outfile"} ? 0 : 1;

  if ($require_outfile ? -s $self->outfile() : 1) {
    # outfile exists, but be sure job is actually done in case
    # program is writing to final output filename  :/
    my $jobfile = $self->get_tracking_file("-suffix" => "job");
    if (-s $jobfile) {
      open(JTMP, $jobfile) || die "can't open $jobfile";
      my $job_id = <JTMP>;
      chomp $job_id;
      close JTMP;

      my $info = $self->get_job_info("-job" => $job_id);
      my $status = $info->{STAT};
#      printf STDERR "status: %s\n", $status;
      if (
	  $status eq "DONE"
	  or $status eq "EXIT"
	  or $status eq ""
	  # if not found by bsub
	 ) {
	$done = 1;
      } elsif ($status eq "RUN" or $status eq "PEND") {
	$done = 0;
      } else {
	die "unhandled job status $status";
      }
    } else {
      # no job file, only outfile
      $done = 1;
    }
  }
  return $done;
}

sub get_job_info {
  my ($self, %options) = @_;
  my $job_id = $options{"-job"} || die "-job";

  open(JCHECK, "bjobs $job_id|") || die;
  my $first = <JCHECK>;
  chomp $first if defined $first;
  my @lines;
  while (<JCHECK>) {
    chomp;
    push @lines, $_;
  }
  close JCHECK;
  my %row;
  if ($first =~ /JOBID/) {
    my @h = split /\s+/, $first;
    @row{@h} = split /\s+/, $lines[0];
  }
  return \%row;
}

sub run_now_or_later {
  # STATIC, exportable
  my (%options) = @_;
  my $now = grep {$_ eq "-now"} @Cluster::ARGV_RAW;

  unless ($now) {
    my @params = map {/\s/ ? '"' . $_ . '"' : $_} @Cluster::ARGV_RAW;
    # quote parameters with whitespace to protect from shell
    push @params, "-now";
#    push @params, "-devel-path" if DevelopmentPath::is_development();

    my $cmd = join " ", $0, @params;

    my $outfile = $options{"-outfile"};
    my $memory = $options{"-memory"} || die "-memory";
    my $project = $options{"-project"} || "PCGP";

    my $c = new Cluster(
      "-command" => $cmd,
      "-memory" => $memory,
      "-outfile" => $outfile,
      "-project" => $project,
#      "-debug" => 1,
#      "-save_output" => 0,
      "-create_bsub_file" => 0,
      "-create_job_file" => 0
	);
    $c->run(
      "-no-outfile" => (defined $outfile ? 0 : 1),
	);

  }

  return $now;
}

sub wait_for_jobs {
  # STATIC
  my (%options) = @_;
  my $jobs = $options{"-jobs"} || die;
  my $sleep = $options{"-sleep"} || 60;
  my $min_age_before_considered_crashed = 190;

  while (1) {
    my $all_done = 1;
    my %done_timestamp;
    my $count_done = 0;
    foreach my $job (@{$jobs}) {
      if (-s $job->outfile) {
	# finished
	$count_done++;
      } elsif ($job->done("-no-outfile" => 1)) {
	# even if job has successfully finished, there may be a delay
	# before the output file is visible on this node (NFS issue?)
	$all_done = 0;
	my $outfile = $job->outfile();
	if (my $time = $done_timestamp{$outfile}) {
	  my $age = time - $time;
	  if ($age < $min_age_before_considered_crashed) {
	    # wait for file to appear on filesystem
	    printf STDERR "waiting for %s for %d...\n", $outfile, $age;
	  } else {
	    # waited too long
	    confess(sprintf "ERROR: job terminated for %s but no outfile!", $job->outfile);
	  }
	} else {
	  printf STDERR "job done but no outfile %s, waiting a little more...\n", $outfile;
	  $done_timestamp{$outfile} = time;
	}
      } else {
	# job still in system, wait
#	print STDERR "job not done\n";
	$all_done = 0;
      }
    }
    if ($all_done) {
      last;
    } else {
      my $total = scalar @{$jobs};
      printf STDERR "have %d/%d (%.1f%%), waiting for %d...\n",
	$count_done, $total,
	  ($count_done * 100 / $total),
	    $sleep;
      sleep $sleep;
    }
  }
}

sub run_jobs {
  # STATIC
  my (%options) = @_;
  my $jobs = $options{"-jobs"} || die;
  my $project = $options{"-project"} || "PCGP";
  # hack
  my $ram = $options{"-memory"} || die "-memory";

  my @jobs;
  foreach my $job (@{$jobs}) {
    my $c = new Cluster(%{$job},
			"-memory" => $ram,
			"-project" => $project,
		       );
    $c->run();
    push @jobs, $c;
  }

  wait_for_jobs("-jobs" => \@jobs,
		"-sleep" => $options{"-sleep"});
}


1;
