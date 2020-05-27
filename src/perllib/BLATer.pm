package BLATer;
# BLAT wrapper

use strict;
use Carp qw(confess);

use Bio::SearchIO;

use Configurable;
use Exporter;
use TemporaryFileWrangler;
use WorkingFile;

@BLATer::ISA = qw(Configurable Exporter);
@BLATer::EXPORT_OK = qw();

use MethodMaker qw(
tfw
tempfile_base

null_mode
verbose

minScore
stepSize
tileSize
		  );

my $VERBOSE = 0;

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  my $tfw = new TemporaryFileWrangler();
#  $tfw->verbose(1);
  $self->tempfile_base($tfw->get_tempfile("-glob" => 1));
  $self->tfw($tfw);
  $self->configure(%options);
  return $self;
}

sub blat {
  my ($self, %options) = @_;

  my $parser;

  if ($self->null_mode) {
    $parser = Bio::SearchIO->new();
    # not sure this will fly
  } else {
    my $db_fa = $self->fa_setup(($options{"-database"} || die "-database"), "database");
    my $query_fa = $self->fa_setup(($options{"-query"} || die "-query"), "query");
    my $outfile = $self->get_tempfile("-append" => ".blat");

    printf STDERR "db=%s query=%s out=%s\n", $db_fa, $query_fa, $outfile if $VERBOSE;
    if ($ENV{DEBUG_BLAT_DB_FASTA}) {
      foreach my $ref (["query", $query_fa],
		       ["db", $db_fa]) {
	my ($label, $file) = @{$ref};
	printf STDERR "%s FASTA:\n", $label;
	open(DEBUGTMP, $file) || die;
	while(<DEBUGTMP>) {
	  print STDERR $_;
	}
      }

    }

    my @cmd = "blat";
    push @cmd, sprintf "-minScore=%d", $self->minScore if $self->minScore();
    push @cmd, sprintf "-stepSize=%d", $self->stepSize if $self->stepSize();
    push @cmd, sprintf "-tileSize=%d", $self->tileSize if $self->tileSize();
    push @cmd, "-out=pslx";
    if ($options{"-protein"}) {
      push @cmd, "-t=prot";
      push @cmd, "-q=prot";
    }
    push @cmd, $db_fa;
    push @cmd, $query_fa;
    push @cmd, $outfile;
    push @cmd, "> /dev/null" unless $VERBOSE;

    my $cmd = join " ", @cmd;

    printf STDERR "running %s\n", $cmd if $VERBOSE;
    system($cmd);

    die sprintf "error running %s, exit %d", $cmd, $? if $?;
    die "where is outfile $outfile" unless -f $outfile;

    if ($self->verbose) {
      printf STDERR "BLAT output:\n";
      open(BLATTMP, $outfile) || die;
      while (<BLATTMP>) {
	print STDERR $_;
      } 
      close BLATTMP;
    }

    $parser = Bio::SearchIO->new(-file => $outfile, -format => 'psl');
    # no specific pslx parser??
  }

  return $parser;
}

sub fa_setup {
  my ($self, $thing, $tag) = @_;
  my $fa_file;
  die unless $tag;
  if (ref $thing) {
    # hash of sequences
    $fa_file = $self->get_tempfile("-append" => sprintf ".%s.fa", $tag);
    write_fasta_set("-file" => $fa_file, "-reads" => $thing);
  } else {
    if (-f $thing) {
      # pre-built filename
      $fa_file = $thing;
    } else {
      # single sequence string
      $fa_file = $self->get_tempfile("-append" => sprintf ".%s.fa", $tag);
      my %reads;
      $reads{query} = $thing;
      write_fasta_set("-file" => $fa_file, "-reads" => \%reads);
    }
  }
  return $fa_file;
}

sub write_fasta_set {
  # hack: should probably figure out how to do w/SeqIO but this is so simple...
  my (%options) = @_;
  my $file = $options{"-file"} || die;
  my $reads = $options{"-reads"} || die;

  my $wf = new WorkingFile($file);
  my $fh = $wf->output_filehandle();

  foreach my $id (sort keys %{$reads}) {
    my $sequence = $reads->{$id};
    # QUESTION: mask *?  how does blat deal with these?
    printf $fh ">%s\n%s\n", $id, $sequence;
  }
  $wf->finish();
}

sub get_tempfile {
  my ($self, %options) = @_;
  my $append = $options{"-append"} || die "-append";
  my $base = $self->tempfile_base() || die;
  return $base . $append;
  # glob mode will clean up these varieties
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
=head1 LICENCE AND COPYRIGHT
Copyright 2019 St. Jude Children's Research Hospital 

Licensed under a modified version of the Apache License, Version 2.0
(the "License") for academic research use only; you may not use this
file except in compliance with the License. To inquire about commercial
use, please contact the St. Jude Office of Technology Licensing at
scott.elmer@stjude.org.

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
