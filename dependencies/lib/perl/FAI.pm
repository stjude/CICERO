package FAI;
# read .fai FASTA index files
# see also http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm#Indexed_Fasta_Files
# (this version apparently not installed, don't want to disrupt)
# 
# see also ~/wc/compbio/bambino/trunk/src/main/java/Ace2/FAIIndexRecord.java
#
# MNE 8/2014

use strict;

use Carp qw(confess);
use POSIX qw(SEEK_SET);

use Configurable;
use Exporter;
use MiscUtils qw(dump_die);
use ReferenceNameMapper;
use TdtConfig;
use CacheManager;

@FAI::ISA = qw(Configurable Exporter);
@FAI::EXPORT_OK = qw(get_fai);

use MethodMaker qw(
fasta
index

sequence_cache

chunk_ref
last_chunk_setup
start_base
rnm

buffered_id
buffered_start
buffered_end
buffered_chunk
buffered_buffer_length

chunk_cache_chromosomes
genome
limit_chrom_cache
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->buffered_buffer_length(1000);
  $self->configure(%options);
  $self->reset_cache();
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $fasta = $self->fasta();
  unless ($fasta) {
    my $genome = $self->genome() || die "-fasta or -genome";
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $self->fasta($fasta = $config_genome->{FASTA} || die);
  }

  die "need -fasta or -genome" unless $fasta;
  confess "where is $fasta" unless -s $fasta;
  my $fai = $fasta . ".fai";
  confess "where is $fai" unless -s $fai;
  open(FAI, $fai) || die;
  my $rnm = new ReferenceNameMapper();
  $self->rnm($rnm);
  my @headers = qw(
      sequence_id
      sequence_length
      file_sequence_offset
      nt_per_line
      bytes_per_line
  );

  my %index;
  while (<FAI>) {
    chomp;
    my %row;
    my @f = split /\t/, $_;
    die unless @f == @headers;
    @row{@headers} = @f;
    my $id = $row{sequence_id} || die;
    $rnm->add_name($id);
#    print STDERR "index $id\n";
    $index{$id} = \%row;
  }
  close FAI;
  $self->index(\%index);
}

sub reset_cache {
  my ($self) = @_;
  my $limit = $self->limit_chrom_cache() || 10000;
  # by default, set very high value so effectively no limit
  my $cm = new CacheManager("-cache_limit" => $limit);
#  $cm->verbose(1);
  $self->sequence_cache($cm);
}

sub get_irow {
  my ($self, %options) = @_;
  my $id_raw = $options{"-id"} || die "-id";
  my $id_idx = $self->find_name($id_raw) || confess "no hit for $id_raw";
  my $irow = $self->index->{$id_idx} || die;
  return $irow;
}

sub get_sequence {
  my ($self, %options) = @_;

  my $irow = $self->get_irow(%options) || die "can't find index entry";
  my $expected_len = $irow->{sequence_length} || die;

  my $seq_ref;
  my $sid = $irow->{sequence_id} || die;

  if ($seq_ref = $self->sequence_cache->get($sid)) {
    # cached
  } else {
#    dump_die(\%options, "fasta load", 1);
    printf STDERR "loading FASTA for %s...", $sid;
    my $offset = $irow->{file_sequence_offset} || die;
    open(FASTA, $self->fasta) || die;
    seek(FASTA, $offset, SEEK_SET) || die "seek failed";
    my $seq = "";
    while (<FASTA>) {
      last if /^>/;
#    print;
      chomp;
      $seq .= $_;
    }
    close FASTA;
    $seq_ref = \$seq;
    $self->sequence_cache->put($sid, $seq_ref);
    print STDERR "done\n";
  }
  die "length mismatch" unless length($$seq_ref) == $expected_len;

  return $seq_ref;
}

sub get_chunk {
  my ($self, %options) = @_;
  my $base_number = $options{"-start"} || die "-start";
  # 1-based
  my $length = $options{"-length"};
  if (not($length) and my $end = $options{"-end"}) {
    $length = ($end - $base_number) + 1;
  }
  confess "specify -length or -end" unless $length;

  my $seq;

  if ($self->chunk_cache_chromosomes()) {
    my $seq_ref = $self->get_sequence(%options) || die;
    $seq = substr($$seq_ref, $base_number - 1, $length);
  } else {
    my $irow = $self->get_irow(%options) || die "can't find index entry";
    my $nt_per_line = $irow->{nt_per_line} || die;
    my $bytes_per_line = $irow->{bytes_per_line} || die;

    my $line_index = int(($base_number - 1) / $nt_per_line);
    my $line_offset = ($base_number - 1) % $nt_per_line;

    my $offset = $irow->{file_sequence_offset} || die;
    $offset += $line_index * $bytes_per_line;

    open(FASTA, $self->fasta) || die "can't open " . $self->fasta;
    seek(FASTA, $offset, SEEK_SET) || die "seek failed";
    $seq = <FASTA>;
    chomp $seq;
    $seq = substr($seq, $line_offset);

    while (1) {
      my $slen = length($seq);
      if ($slen == $length) {
	last;
      } elsif ($slen > $length) {
	my $extra = $slen - $length;
	$seq = substr($seq, 0, $length);
	last;
      } else {
	my $line = <FASTA>;
	chomp $line;
	if ($line) {
	  if ($line =~ /^>/) {
	    last;
	  } else {
	    $seq .= $line;
	  }
	} else {
	  last;
	}
      }
    }
  }

  return $seq;
}

sub sanity_check {
  # compare user-provided reference sequence vs. actual
  my ($self, %options) = @_;
  my $user_seq = $options{"-sequence"} || die "-sequence";
  my $actual_seq = $self->get_chunk(%options, "-length" => length($user_seq));

  foreach ($user_seq, $actual_seq) {
    $_ = uc($_);
  }

  return $user_seq eq $actual_seq ? 1 : 0;
  # TO DO: masking, etc.?
}

sub chunk_setup {
  my ($self, %options) = @_;
  if (my $last = $self->last_chunk_setup()) {
    # 
    my $all_agree = 1;
#    dump_die(\%options, "this", 1);
    while (my ($k, $v) = each %options) {
      if (!(exists $last->{$k} and $options{$k} eq $last->{$k})) {
	$all_agree = 0;
	last;
      }
    }
    if ($all_agree) {
      # cache hit
#      printf STDERR "FAI chunk cache hit\n";
      return;
    }
  }
  $self->last_chunk_setup(\%options);

  my $chunk = $self->get_chunk(%options);
  my $start = $options{"-start"} || die "-start";

  $self->chunk_ref(\$chunk);
  $self->start_base($start);
}

sub get_chunked_base {
  my ($self, %options) = @_;
  my $start = $options{"-start"} || die "-start";
  my $length = $options{"-length"} || die "-length";
  my $chunk_ref = $self->chunk_ref() || die;
  my $start_base = $self->start_base();
  die unless defined $start_base;

  return substr($$chunk_ref, $start - $start_base, $length);
}

sub get_chunk_buffered {
  my ($self, %options) = @_;
  my $id = $options{"-id"} || die "-id";
  my $start = $options{"-start"} || die "-start";
  my $length = $options{"-length"} || die "-length";
  my $start_base = $self->start_base();
  my $buffered_start = $self->buffered_start();

  my $cached;
  if ($buffered_start and ($self->buffered_id() || "") eq $id) {
    # previously queried this chromosome
    my $qend = $start + ($length - 1);
    if ($start >= $buffered_start and $qend <= $self->buffered_end) {
#      print STDERR "cached\n";
#      die "cached";
    } else {
#      print STDERR "new query\n";
    }
  }

  unless ($cached) {
    my $blen = $self->buffered_buffer_length || die;
    $buffered_start = $start - $blen;
    my $bend = $options{"-end"} || $start + ($length - 1);
    $bend += $blen;
    # FIX ME: trim end too!!
    $buffered_start = 1 if $buffered_start < 1;
#    print STDERR "query $buffered_start $bend\n";
    
    my %o = %options;
    $o{"-start"} = $buffered_start;
    $o{"-end"} = $bend;
    delete $o{"-length"};
    my $chunk = $self->get_chunk(%o);
    $self->buffered_chunk(\$chunk);
    $self->buffered_start($buffered_start);
    $self->buffered_end($bend);
    $self->buffered_id($id);
  }

  my $chunk_ref = $self->buffered_chunk() || die;

  return substr($$chunk_ref, $start - $buffered_start, $length);
  
}

sub find_name {
  my ($self, $name) = @_;
  return $self->rnm->find_name($name);
}

sub generate_chunk_query {
  my ($self) = @_;
  my $chunk_ref = $self->chunk_ref() || die;
  my $start_base = $self->start_base();

  my $code = sub {
    # args:
    #  (1) query basse number (1-based)
    #  (2) length
    return substr($$chunk_ref, $_[0] - $start_base, $_[1]);
  };
  return $code;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
