package RefFlat2Genome;
# map refFlat gene models to genomic bases
# MNE 8/2014

use strict;
use Carp qw(confess);

use List::Util qw(min max);

use Configurable;
use Exporter;

use Bio::Tools::CodonTable;

use DelimitedFile;
use RefFlatFile;
use GenomeUtils qw(complement);
use MiscUtils qw(dump_die);
use FAI;

use constant STATUS_OK_PERFECT => 1;
# perfect match between genome-translated AA and canonical AA
use constant STATUS_OK_TRUNCATED => 2;
# genome-translated AA matches canonical AA, but complete translation
# coordinates are not provided (e.g. CRLF2).  Per JZ this can
# legitimately happen because full genomic sequence is not always
# available.

# TO DO: type for case that's truncated and also has a few mismatches?

use constant STATUS_ERROR_NO_AA_DATA => -1;
use constant STATUS_ERROR_AA_MISMATCH_MINOR => -2;
# same AA length, but a small number of disagreements
use constant STATUS_ERROR_AA_MISMATCH_MAJOR => -3;
# same AA length, but a higher level of mismatches
use constant STATUS_ERROR_AA_LENGTH_MISMATCH_SINGLE_TRANSCRIPT => -4;
use constant STATUS_ERROR_AA_LENGTH_MISMATCH_MULTI_TRANSCRIPT => -5;
use constant STATUS_ERROR_CDS_SYNC => -6;
# might not be an error, could be NR_, e.g. NR_000005

use constant MIN_AA_MISMATCH_FREQUENCY_FOR_MAJOR_MISMATCH => .04;
# level of mismatches between my (genomic-based) AA translation
# and the canonical AA from the refSeq to consider problematic
# enough that the mapping should not be used.

use constant MAX_AA_MISMATCHES_FOR_MINOR_MISMATCH => 2;
# threshold to consider minor/major mismatch

@RefFlat2Genome::ISA = qw(Configurable Exporter);
@RefFlat2Genome::EXPORT_OK = qw();

use MethodMaker qw(
	refflat
        refgene2protein
fasta

        by_accession
        r2p
        r2gene
        fai
verbose

cache_genome_mappings
infer_aa
rff
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);

  my $fasta = $self->fasta() || die "-fasta";
  $self->fai(new FAI("-fasta" => $fasta));

  return $self;
}

sub parse {
  my ($self, %options) = @_;

  my $fn_rf = $options{"-refflat"} || die "-refflat";

  my $rf = new RefFlatFile();
  $rf->canonical_references_only(1);
  $rf->parse_file(
		  "-refflat" => $fn_rf,
		  "-type" => "refgene"
		 );
  $self->rff($rf);
  my %index;
  foreach my $row (@{$rf->rows}) {
    my $key = $row->{name} || die;
    next if $key =~ /^NR_/ and not $self->infer_aa();
    # ignore non-coding RNAs
    push @{$index{$key}}, $row;
  }
  $self->by_accession(\%index);

  my $r2p_fn = $options{"-refgene2protein"} || die "-refgene2protein";
  printf STDERR "loading %s...", $r2p_fn;
  my $df = new DelimitedFile("-file" => $r2p_fn,
			     "-headers" => 1,
			    );
  my %r2p;
  my %r2gene;
  while (my $row = $df->get_hash()) {
    # to do: can this logic
    my $acc = $row->{accession} || die "no accession field";
    die "duplicate" if $r2p{$acc};
    $r2p{$acc} = $row->{protein};

    my $gene = $row->{gene} || die "no gene field";
    $r2gene{$acc} = $gene;
  }
  $self->r2p(\%r2p);
  $self->r2gene(\%r2gene);
  print STDERR "done\n";
}

sub find {
  # by default this returns an arrayref since a transcript may be
  # mapped to more than one position
  my ($self, %options) = @_;

  my $acc = $options{"-accession"} || confess "-accession";
  my $set_raw = $self->by_accession->{$acc} || die "can't find accession $acc";

  my @set;
  if ($self->cache_genome_mappings) {
    # save genome mappings in refFlat hashes.
    # this is memory intensive!
    @set = @{$set_raw};
  } else {
    foreach my $raw (@{$set_raw}) {
      # don't return raw hashes, instead create copies where genome
      # mapping info will be added.  This is memory intensive -- if we
      # cache process memory will grow substantially if we cache results.
      die if $raw->{is_genome_mapped};
      # should NOT be saved in raw data!
      push @set, { %{$raw} };
    }
  }

  if (my $reference = $options{"-reference"}) {
    # optional: filter results to those mapped on specifed chrom only
    my @filtered;

    foreach my $hit (@set) {
      push @filtered, $hit if chrom_compare($reference, $hit->{chrom} || die);
    }
    @set = @filtered;
  }

  $self->genome_map(%options,
		    "-set" => \@set);
  # map on demand and save to temporary hashes only

  return \@set;
}

sub genome_map {
  my ($self, %options) = @_;
  my $acc = $options{"-accession"} || die;
#  my $set = $self->by_accession->{$acc} || die "can't find $acc";
  my $set = $options{"-set"} || die "-set";
  # temporary copies for user results
  my $aa = $self->r2p->{$acc};
  my $infer_aa = $self->infer_aa;

  if (not($infer_aa) and not($aa)) {
    printf STDERR "WARNING: no AA available for %s, can't verify translation!\n", $acc;
    # e.g. 
    # NM_000658.2: This RefSeq was permanently suppressed because currently there is insufficient support for the transcript and the protein.
  }

  my $verbose = $self->verbose;
  my $map_count = scalar @{$set};

  foreach my $row (@{$set}) {
    next if $row->{is_genome_mapped};

    my $seq_ref = $self->fai->get_sequence("-id" => $row->{chrom});

    my $cds_start = $row->{cdsStart} + 1;
    # convert from interbase to in-base
    my $exon_count = scalar @{$row->{exons}};

    if ($verbose) {
      printf STDERR "processing %s on %s\n", $acc, $row->{chrom};
      printf STDERR "CDS: %d-%d\n", @{$row}{qw(cdsStart cdsEnd)};
      printf STDERR "tx: %d-%d\n", @{$row}{qw(txStart txEnd)};

      foreach my $ref (@{$row->{exons}}) {
	printf STDERR "raw exon: %d-%d\n", $ref->{start}, $ref->{end};
      }
    }

    my $ei_check = 0;
    foreach my $exon (@{$row->{exons}}) {
      printf STDERR "exon: %d-%d\n", $exon->{start}, $exon->{end} if $verbose;
      last if $cds_start >= $exon->{start} && $cds_start <= $exon->{end};
      # possible edge case here??
      $ei_check++;
    }
    # index of exon containing codon start

#    die "start exon != 0 for $acc: $ei" unless $ei == 0;
    # NM_000015

    if ($ei_check >= $exon_count) {
      # sanity check
      printf STDERR "ERROR: fail to find exon with cds_start for %s!\n", $acc;
      $row->{aa_trans_status} = STATUS_ERROR_CDS_SYNC;
      # may never happen now that we're skipping NR_ sequences
      printf STDERR "AA translation status for %s: %d\n", $acc, $row->{aa_trans_status};
      # hack
      next;
    }

    my @map;
    my $in_cds = 0;
    my $ei = 0;
    my $exon = $row->{exons}->[$ei];
    my $ref_base_num = $exon->{start};
    # already adjusted to 1-based

    while (1) {
      $in_cds = 1 if $ref_base_num == $cds_start;

      my %info;
      $info{ref_base_num} = $ref_base_num;
      $info{exon_number} = $ei + 1;
      $info{in_cds} = $in_cds;

#      printf STDERR "saved entry at %d\n", $ref_base_num;
      push @map, \%info;

      if ($ref_base_num >= $exon->{end}) {
	$ei++;
	if ($ei >= $exon_count) {
	  printf STDERR "out of exons, stopping translation...\n";
	  last;
	} else {
	  $exon = $row->{exons}->[$ei];
	  $ref_base_num = $exon->{start};
	}
      } else {
	$ref_base_num++;
      }
    }

    my $is_rc;
    if ($row->{strand} eq "+") {
      # OK
    } elsif ($row->{strand} eq "-") {
      #
      # transcript on - strand: reorient data
      #
      $is_rc = 1;
      my $cds_end = $row->{cdsEnd} || die;
      my @filtered;
      my $in_cds = 0;
      my %saw_exon;
      my $new_exno = 0;
      foreach my $entry (reverse @map) {
        # reverse base ordering
	$in_cds = 1 if $entry->{ref_base_num} == $cds_end;
	$entry->{in_cds} = $in_cds;
	my $old_exno = $entry->{exon_number} || die;
	unless ($saw_exon{$old_exno}) {
	  $saw_exon{$old_exno} = 1;
	  $new_exno++;
	}
	$entry->{exon_number} = $new_exno;
	push @filtered, $entry;
      }
      @map = @filtered;

    } else {
      die;
    }

    my $cds_start_i;
    for ($cds_start_i = 0; $cds_start_i < @map; $cds_start_i++) {
      last if $map[$cds_start_i]->{in_cds};
    }

    my $cds_phase_start_i = $cds_start_i % 3;
    # index to start fake AA translation in so that we are
    # in phase with the coding region.
    # these simulated AA sequences can be used for SV in-frame checks.

    my $start_codon_number = - int($cds_start_i / 3);
    # includes leading negative (and 0) codon numbers
    # for 5' UTR.

    # assign codon numbers:
    my $nt_count = 0;
    my $sanity_checked;
    my $coding_base_count = 0;
    for (my $i = 0; $i < @map; $i++) {
      my $entry = $map[$i];
      if ($i >= $cds_phase_start_i) {
	# in phase
	$entry->{codon_number} = int($nt_count / 3) + 1 + $start_codon_number;
	$entry->{codon_base} = ($nt_count % 3) + 1;
	$nt_count++;

	if ($entry->{in_cds}) {
	  $entry->{coding_base_number} = ++$coding_base_count;
	  # http://www.hgvs.org/mutnomen/recs-DNA.html
	  # - there is no nucleotide 0
	  # - nucleotide 1 is the A of the ATG-translation initiation codon 
	  unless ($sanity_checked) {
	    die "sanity fail, CDS start codon # != 1!" unless $entry->{codon_number} == 1;
	    $sanity_checked = 1;
	  }
	}
      }
    }

#  printf STDERR "CDS: start=%d end=%d\n", $cds_start, $cds_end;
#  printf STDERR "AA=%s\n", $aa;

    my $ct = new Bio::Tools::CodonTable();

    my @codon;
    my @translated;
    my $stopped;
    my $map_i = 0;
    foreach my $entry (@map) {
      # 5' -> 3', in transcript order
      unless ($entry->{in_cds}) {
	# CDS translation hasn't started yet
	$entry->{in_utr} = 5;
	printf STDERR "5' UTR base %d\n", $entry->{ref_base_num} if $verbose;
#	next;
      }
      if ($stopped) {
	# CDS translation finished
	$entry->{in_utr} = 3;
	$entry->{coding_base_number} = "";
	printf STDERR "3' UTR base %d\n", $entry->{ref_base_num} if $verbose;
#	next;
      }

      unless ($map_i++ >= $cds_phase_start_i) {
	# don't start translating codons until we reach the
	# first base that's in-frame with the start codon
	printf STDERR "skipping leading out-of-frame 5' UTR base\n";
	next;
      }

      push @codon, $entry;
      if (@codon == 3) {
	# enough bases for a full codon, translate
	my $nt = "";
	my @bases;
	my @coding_base_numbers;
	foreach my $e (@codon) {
	  my $bn = $e->{ref_base_num};
	  my $ref_base = substr($$seq_ref, $bn - 1, 1);
	  $e->{ref_base_genome} = $ref_base;
	  $ref_base = complement($ref_base) if $is_rc;
	  # already reversed above, just needs complement
	  $e->{ref_base_transcript} = $ref_base;
	  # corrected for strand
	  $nt .= $ref_base;
	  push @bases, $bn;
	  push @coding_base_numbers, $e->{coding_base_number} || "";
	}
	my $a = $ct->translate($nt);
	foreach my $e (@codon) {
	  $e->{AA} = $a;
	}

	printf STDERR "%s (%s) => %s => %s (%d)\n",
	join(",", @bases),
	join(",", @coding_base_numbers),
	$nt, $a, $codon[0]->{codon_number} if $verbose;

	if ($entry->{in_cds}) {
	  # coding sequence
	  $stopped = 1 if $ct->is_ter_codon($nt);
	  push @translated, $a unless $stopped;
	}

	@codon = ();
      }
    }

    #
    #  generate fake translations for 5' UTR.
    #  these may be used in in-frame detection.
    #


    if (0) {
      foreach my $entry (@map) {
	dump_die($entry, "debug", 1);
      }
    }

    $row->{codon_map} = \@map;

    my $aa_trans = join "", @translated;
#    printf STDERR "generated AA for %s: %s\n", $acc, $aa_trans;

    if (not($aa) and $self->infer_aa()) {
      # if we don't have the formal AA sequence (e.g. from refSeq),
      # infer from the genome mapping.  Might be fine but could also
      # be garbage, no way to verify (QC not available)
      $aa = $aa_trans;
      $self->r2p->{$acc} = $aa;
      # save prediction for later use e.g. by sv_inframe.pl
    }


    my $aa_trans_status;
    if ($aa) {
      if ($aa eq $aa_trans) {
	$aa_trans_status = STATUS_OK_PERFECT;
      } else {
	printf STDERR "ERROR: AA translation sync failure for %s!:\n  expected:\n%s\n  got:\n%s\n", $acc, $aa, $aa_trans;
	if (length($aa) == length($aa_trans)) {
	  #
	  # length matches: might just be a subtle disagreement
          #
	  my $mismatches = 0;
	  my $len = length($aa);
	  for (my $i = 0; $i < $len; $i++) {
	    $mismatches++ unless substr($aa, $i, 1) eq substr($aa_trans, $i, 1);
	  }

	  my $mm_freq = $mismatches / length($aa);

	  printf STDERR "mismatch count for %s: len %d mismatches %d freq %.04f\n", $acc, length($aa), $mismatches, $mm_freq;

	  if ($mm_freq >= MIN_AA_MISMATCH_FREQUENCY_FOR_MAJOR_MISMATCH) {
	    $aa_trans_status = STATUS_ERROR_AA_MISMATCH_MAJOR;
	  } else {
	    # acceptable level of mismatches
	    # not sure why these happens; genome revisions between
	    # sequence used for NM_ and major genome build?
	    $aa_trans_status = STATUS_ERROR_AA_MISMATCH_MINOR;
	  }
	} else {
	  #
	  # length doesn't match
	  #
	  if (index($aa, $aa_trans) != -1) {
	    # mapping only available for a subset of the canonical
	    # sequence: this can happen sometimes
	    $aa_trans_status = STATUS_OK_TRUNCATED;
	  }

	  if (not($aa_trans_status) and length($aa_trans) < length($aa)) {
	    my $mismatches = 0;
	    my $len = length($aa_trans);
	    my $clean = 0;
	    for (my $i = 0; $i < $len; $i++) {
	      $mismatches++ unless substr($aa, $i, 1) eq substr($aa_trans, $i, 1);
	      $clean++ unless $mismatches;
	    }
	    printf STDERR "acc:%s length:%d clean:%d mismatches:%d\n", $acc, $len, $clean, $mismatches;
	  }

	  unless ($aa_trans_status) {
	    $aa_trans_status = $map_count == 1 ? STATUS_ERROR_AA_LENGTH_MISMATCH_SINGLE_TRANSCRIPT : STATUS_ERROR_AA_LENGTH_MISMATCH_MULTI_TRANSCRIPT;
	    # for debugging purposes track accessions with single mappings
	    # separately.  Seems more critical if the only available mapping
	    # fails than of one of several fail.
	    printf STDERR "ERROR: AA translation length mismatch for %s! expected:%d got:%d strand=%s code=%d\n", $acc, length($aa), length($aa_trans), $row->{strand}, $aa_trans_status;
	  }
	}
      }
    } else {
      $aa_trans_status = STATUS_ERROR_NO_AA_DATA;
      # no data, can't compare
    }
    $row->{aa_trans_status} = $aa_trans_status;
    printf STDERR "AA translation status for %s: %d\n", $acc, $aa_trans_status;

    $row->{is_genome_mapped} = 1;
  }
}

sub chrom_compare {
  my ($chr_a, $chr_b) = @_;
#  printf STDERR "compare %s:", join " ", $chr_a, $chr_b;
  foreach ($chr_a, $chr_b) {
    $_ = lc($_);
    s/^chr//;
  }
  my $c = $chr_a eq $chr_b;
#  printf STDERR "%d\n", $c || 0;
  return $c;
}

sub get_aa {
  my ($self, %options) = @_;
  my $acc = $options{"-accession"} || die;
  return $self->r2p->{$acc};
}

sub is_intronic {
  # does the specified site fall within the transcript but not in an exon?
  my ($self, %options) = @_;
  my $hit = $options{"-entry"} || die;
  my $pos = $options{"-base-number"} || die;
  my $map = $hit->{codon_map} || die;
  my %all_bases;
  # currently all bases; option for CDS only?
  foreach my $entry (@{$map}) {
    my $bn = $entry->{ref_base_num} || die;
    $all_bases{$bn} = 1;
  }
  my @all_b = keys %all_bases;
  my $min = min(@all_b);
  my $max = max(@all_b);
  # compensate for strand

  my $result;
  if ($pos >= $min and $pos <= $max) {
    $result = $all_bases{$pos} ? 0 : 1;
  } else {
    # outside of transcript mapping
    $result = undef;
  }
  return $result;
}

sub get_gene {
  my ($self, $nm) = @_;
  return $self->r2gene->{$nm};
}

sub get_reference_sequence {
  my ($self, $id) = @_;
  return $self->fai->get_sequence("-id" => $id);
}

sub get_utr5_aa_from_map {
  my ($self, %options) = @_;
  my $hit = $options{"-hit"} || die "-hit";
  my $map = $hit->{codon_map} || die;
  my %codon2aa;
  foreach my $ref (@{$map}) {
    next if $ref->{in_utr} == 3;
    my $aa = $ref->{AA} || next;
    my $cnum = $ref->{codon_number};
    $codon2aa{$cnum} = $aa;
  }

  return join "", map {$codon2aa{$_}} sort {$a <=> $b} keys %codon2aa;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
