package WorkingFile;
# write to a temporary file until completely written, then
# rename to a permanent file.  Prevents the "complete" file
# from existing until it's complete.

use strict;
use Carp qw(confess);
use FileHandle;

#use misc;
use FileUtils qw(aggressive_rename);
use Configurable;

@WorkingFile::ISA = qw(Configurable);

use MethodMaker qw(
		   complete_name
		   tempfile
		   fh
		   compress
                   level

                   skip_if_exists_mode
		   tabix
		  );


sub new {
  my ($type, $complete_name, %options) = @_;
  die "final filename required" unless $complete_name;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->complete_name($complete_name);
  return $self;
}

sub outfile_setup {
  my ($self) = @_;
  my $complete_name = $self->complete_name();
  if (my $type = $self->compress()) {
    if ($type eq "bz2") {
      $complete_name .= ".bz2" unless $complete_name =~ /\.bz2$/i;
    } else {
      $complete_name .= ".gz" unless $complete_name =~ /\.gz$/i;
    }
    $self->complete_name($complete_name);
  }

  my $fn = $self->complete_name() . ".tmp";
  # don't use tmpnam() and copy to target because file might be very large...
  # if the tempfile is in same location we can just rename it, which is 
  # much faster and uses no temporary disk space.
  $self->tempfile($fn);

  return $fn;
}

sub output_filehandle {
  my ($self, %options) = @_;
  my $fh = $self->fh;
  unless ($fh) {
    my $fn = $self->outfile_setup();
    $fh = new FileHandle();
    my $cmd;
    if (my $type = $self->compress) {
      my $level = $self->level || 1;
      if ($type eq "bz2") {
	$cmd = sprintf "|bzip2 -%d >%s", $level, $fn;
      } elsif ($type eq "bgzip") {
	$cmd = sprintf "|bgzip >%s", $fn;
      } elsif ($type eq "gz" or $type eq "1") {
	$cmd = sprintf "|gzip -%d >%s", $level, $fn;
      } else {
	die "unhandled compression type $type";
      }
    } elsif ($options{"-rw"}) {
      $cmd = ($options{"-cmd"} || "") . "+>$fn";
    } else {
      $cmd = ($options{"-cmd"} || "") . ">$fn";
    }
#    printf STDERR "cmd: %s\n", $cmd;
    $fh->open($cmd) || confess "can't open $cmd: $!";
    die "ERROR" if $?;
    $self->fh($fh);
  }
  return $fh;
}

sub finish {
  my ($self) = @_;
  my $fh = $self->fh;
  my $tempfile = $self->tempfile || confess "wtf, no tempfile in WF finish!";
  my $complete = $self->complete_name;
  if ($fh) {
    $fh->close || confess "close failed: $!";
    die "ERROR" if $?;
  }
  aggressive_rename($tempfile, $complete,
		    "-skip-if-exists" => $self->skip_if_exists_mode());

#  rename($tempfile, $complete) || die "rename failed";
  unlink $tempfile;
  if (my $tt = $self->tabix) {
    my $prefix;
    if ($tt eq "vcf") {
      # ok
    } else {
      die "-tabix must be vcf";
    }

    my $cmd = sprintf 'tabix -p %s %s', $tt, $complete;
    system $cmd;
    die "$cmd exited with $?" if $?;
  }

}

sub DESTROY {
  if (my $fn = $_[0]->tempfile) {
    unlink $fn if -f $fn;
  }
}

1;
