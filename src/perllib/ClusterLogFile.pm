package ClusterLogFile;
# simple LSF cluster log file parsing
# mne 10/2012

use strict;
use Carp qw(confess);

use Configurable;

@ClusterLogFile::ISA = qw(Configurable Exporter);
@ClusterLogFile::EXPORT_OK = qw();

use MethodMaker qw(
	cluster
        filename

	ok
	errors
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->parse();
  return $self;
}

sub parse {
  my ($self, %options) = @_;
  if (my $c = $self->cluster()) {
    $self->filename($c->get_tracking_file("-suffix" => "out"));
  }
  my $filename = $self->filename() || die "-cluster or -filename";
  my $ok;
  my %errors;

  if (-e $filename) {
    open(CLOG, $filename) || confess "can't open $filename";
    while (<CLOG>) {
      chomp;
      last if /^The output \(if any\) follows:/;
      if (/job killed/) {
	$errors{killed} = $_;
	$ok = 0;
      } elsif (/^Exited with exit code \w+\./) {
	$errors{exited} = $_;
	$ok = 0;
      } elsif (/^Successfully completed/) {
	$ok = 1;
      }
    }
    close CLOG;
    die "can't determine cluster job outcome in $filename" unless defined($ok);
  } else {
    printf STDERR "log file %s doesn't exist, assuming OK\n", $filename;
    $self->ok(1);
  }
  $self->ok($ok);
  $self->errors(\%errors);
}

sub get_error_summary {
  my ($self) = @_;
  my $errors = $self->errors();
  return join " ", map {$_ . ": " . $errors->{$_}} sort keys %{$errors};
}

1;
