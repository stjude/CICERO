package GeneAnnotation;
# find HUGO gene annotations by location
#
# TO DO:
#  - use official HUGO data?
#
# MNE 8/2013

use strict;
use Carp qw(confess);

use MiscUtils qw(dump_die);

@GeneAnnotation::ISA = qw(Configurable Exporter);
@GeneAnnotation::EXPORT_OK = qw();

use Configurable;
use Exporter;

use TdtConfig;
use BucketMap;
use GenomeUtils qw(cook_chromosome_name);
use GeneListCollapser;
#use RefFlatFile qw(uniquify_accessions);
# ?? not sure why this isn't working
use RefFlatFile;


use MethodMaker qw(
all_genes

style

gene_exon_region_dir

genome
bucket_map
is_initialized

results_genes
results_rows
results_genes_genomic_order

single_gene_mode

ignore_non_coding

search_start
search_end

intronic_debug
refgene_flatfile

verbose
is_interbase
strand

nm2rows
clean_sharp
unique_refflat_accessions
		  );


use DBTools qw(
		get_dbi_hg18
		get_dbi_hg19
		get_dbi_mm9
		selectall_hashref
);

use constant CHUNK_SIZE => 100000;

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->ignore_non_coding(1);
  # 9/4/2013:
  # there may be contentions between coding and noncoding genes.
  # e.g. X:70463693-70463693 hits both ZMYM3 (NM_) BCYRN1 (NR_, noncoding).
  # A more sophisticated approach would be to only eliminate the noncoding
  # entries if contention is required.  However, I am on vacation so
  # all this is gravy.
  $self->clean_sharp(1);
  $self->configure(%options);
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  return if $self->is_initialized();

  my $style = $self->style() || die "specify -style [refgene|gene_exon_region|refgene_flatfile]";
  my $rows;

  my $VERBOSE = $self->verbose();

  my $is_interbase;

  if ($style eq "refgene") {
    $is_interbase = 1;
    my $genome = $self->genome() || die "specify -genome";
    my $dbi;
    if ($genome eq "hg19" or $genome eq "GRCh37-lite") {
      $dbi = get_dbi_hg19();
    } elsif ($genome eq "hg18") {
      $dbi = get_dbi_hg18();
    } elsif ($genome eq "mm9" or $genome eq "MGSCv37") {
      $dbi = get_dbi_mm9();
    } else {
      die "unknown database for genome $genome";
    }
    $rows = selectall_hashref($dbi, "select chrom,txStart,txEnd,name2,name from refGene");
  } elsif ($style eq "gene_exon_region") {
    $is_interbase = 0;
    my $dir = $self->gene_exon_region_dir() || confess "-gene_exon_region_dir";
    my @files = glob($dir . "/*region.txt");
    die unless @files;
    $rows = [];
    foreach my $fn (@files) {
      open(GERTMP, $fn) || die;
      my @ger;
      while (<GERTMP>) {
	chomp;
	my @f = split /\t/, $_;
	die unless @f == 8;
	my %r;
	@r{qw(
name2
exon_number
name
chrom
txStart
txEnd
unknown
source
)} = @f;
	# - an accession may be mapped to multiple chroms
	# - it may also be mapped multiple times WITHIN a chrom!
	#   e.g. NM_001098844
	# - NEED TO SORT BY ACCESSION -> START POSITION -> EXON NUMBER,
	#   resetting range when necessary
	push @ger, \%r;
      }
      close GERTMP;
#      die scalar @ger;

      my $last_acc = "bogus";
      my $last_exno = 0;
      my $current_row;

      foreach my $row (sort {
	$a->{name} cmp $b->{name} ||
	    $a->{txStart} <=> $b->{txStart} ||
	    $a->{exon_number} <=> $b->{exon_number}
		       } @ger) {
	printf STDERR "ROW: %s %s %s %s %d-%d\n", @{$row}{qw(name name2 exon_number chrom txStart txEnd)} if $VERBOSE;
	if ($row->{name} ne $last_acc or
#	    $row->{exon_number} < $last_exno) {
	    # FAIL: NR_036199 has 1 exon and is mapped twice on the same chrom
	    $row->{exon_number} <= $last_exno) {
	  $current_row = $row;
	  $current_row->{raw_rows} = [ { %{$row} } ];
	  push @{$rows}, $current_row;
	} else {
	  # same ref
	  push @{$current_row->{raw_rows}}, { %{$row} };
	  $current_row->{txEnd} = $row->{txEnd};
	  die "gurgle" if $current_row->{txStart} > $current_row->{txEnd};
	}

	$last_acc = $row->{name};
	$last_exno = $row->{exon_number};
      }

      if (0) {
	foreach my $r (@{$rows}) {
	  printf "FINAL %s %s %s %d-%d\n",
	  @{$r}{qw(name name2 chrom txStart txEnd)};
	}
      }
    }
  } elsif ($style eq "refgene_flatfile") {
    $is_interbase = 1;
    $rows = [];
    my $fn = $self->refgene_flatfile();
    unless ($fn) {
      my $genome = $self->genome || die "no genome or file";
      my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
      $fn = $config_genome->{REFSEQ_REFFLAT} || die;
    }

    my $strand_filter = $self->strand();

    my @headers_refgene = qw(
		      bin
		      name
		      chrom
		      strand
		      txStart
		      txEnd
		      cdsStart
		      cdsEnd
		      exonCount
		      exonStarts
		      exonEnds
		      score
		      name2
		      cdsStartStat
		      cdsEndStat
		      exonFrames
		   );

    my @headers_refseq = qw(
		      name2
		      name
		      chrom
		      strand
		      txStart
		      txEnd
		      cdsStart
		      cdsEnd
		      exonCount
		      exonStarts
		      exonEnds
		   );
    # different field count, gene symbol in different column

    if (ref $fn and ref $fn eq "ARRAY") {
      # already-parsed rows, e.g. from RefFlatFile.pm
      $rows = $fn;
    } else {
      open(RGTMP, $fn) || die "can't open $fn";
      my @headers;
      my $first = 1;
      my $skipped = 0;
      while (<RGTMP>) {
	chomp;
	if ($first) {
	  $first = 0;
	  next if /^\#/;
	  # /nfs_exports/genomes/1/Homo_sapiens/Hg19/mRNA/RefSeq/refFlat.txt
	}
	my @f = split /\t/, $_;
	unless (@headers) {
	  if (@f == 16) {
	    # refGene
	    @headers = @headers_refgene;
	  } elsif (@f == 11) {
	    # UCSC? / refSeq "sharp"
	    @headers = @headers_refseq;
	  } else {
	    die sprintf "field count mismatch! found %d expected %d", scalar(@f), scalar(@headers_refgene);
	  }
	}
	
	my %r;
	@r{@headers} = @f;
#	dump_die(\%r);

	my $usable = 1;

	if ($strand_filter) {
	  my $strand = $r{strand} || die;
	  $usable = 0 unless $strand eq $strand_filter;
	}

	if ($usable) {
	  push @{$rows}, \%r;
	} else {
	  $skipped++;
	}
      }
      close RGTMP;
      printf STDERR "skipped because on wrong strand: %d\n", $skipped if $skipped;
    }
  } else {
    die "unimplemented style $style";
  }

  RefFlatFile::uniquify_accessions($rows) if $self->unique_refflat_accessions();

  die "interbase not defined" unless defined $is_interbase;
  $self->is_interbase($is_interbase);

  my %map;
  my %all_genes;
  my %nm2rows;
  my %saw;

  #
  #  hash gene locations:
  #
  my %prefix2coding = (
    "NM" => 1,
    "NR" => 0,
    "NG" => 0,
      );

  my $genome = $self->genome;
  my $clean_sharp = $self->clean_sharp();

  foreach my $row (@{$rows}) {
    unless ($row->{name2}) {
      $row->{name2} = $row->{name} || "unknown";
      # normally this should never happen, however some sources e.g.
      # COMBINED_REFFLAT mix data from different formats, some of which
      # are missing this field.
    }

    if ($clean_sharp and $row->{name2}) {
      # remove refFlat "sharp" uniqueness suffixes
      $row->{name2} =~ s/_loc\w+$//;
    }

    my $chr_raw = $row->{chrom} || die;

    my $chr = cook_chromosome_name($chr_raw,
				   "-genome" => $genome,
				   "-return-unknown" => 1) || die;

#    printf STDERR "range: %s (raw=%s) %s %s %d-%d\n", $chr, $chr_raw, $row->{name}, $row->{name2}, $row->{txStart}, $row->{txEnd};

    my $map = $map{$chr};
    unless ($map) {
      $map = $map{$chr} = new BucketMap("-chunk" => CHUNK_SIZE);
    }

    if ($self->ignore_non_coding) {
      my $refseq = $row->{name};
      my $is_coding;
      if ($refseq =~ /^([A-Z]+)_/) {
	my $prefix = $1;
	$is_coding = $prefix2coding{$prefix};
	die "coding prefix not defined for $prefix" unless defined $is_coding;
      } else {
	# 10/2020: CICERO's refFlat file may now have non-refSeq records
	# appended from e.g. ENSEMBL
	printf STDERR "WARNING: accession %s not in refSeq format, assuming non-coding\n", $refseq;
      }
      next unless $is_coding;
    }

#    my $key = join "_", @{$row}{qw(name2 chrom txStart txEnd)};
#    next if $saw{$key};
    # 9/2014: this suppresses some legitimate isoforms, so disabling:
    # 
    #
    # CRLF2   NM_001012288    chrX    -       1314893 1331616 1314893 1325338 5      1314893,1317418,1321271,1325325,1331448, 1315014,1317581,1321405,1325492,1331616,
    # CRLF2   NM_022148       chrX    -       1314893 1331616 1314893 1331527 6      1314893,1317418,1321271,1325325,1327698,1331448, 1315014,1317581,1321405,1325492,1327801,1331616,
    #
    # however I can't remember why these were excluded in the first place.
#    $saw{$key} = 1;

    printf STDERR "final range %s\n", join ",", $row->{name}, $row->{name2}, $row->{chrom}, $row->{txStart}, $row->{txEnd} if $VERBOSE;

    push @{$all_genes{$row->{name2}}}, $row if $row->{name2};
    push @{$nm2rows{$row->{name}}}, $row;

    $map->add_range(
      "-start" => $row->{txStart},
      "-end" => $row->{txEnd},
      "-value" => $row
	);
  }
  $self->bucket_map(\%map);
  $self->all_genes(\%all_genes);
  $self->nm2rows(\%nm2rows);
  $self->is_initialized(1);
}

sub find {
  # find by genomic location
  my ($self, %options) = @_;
  $self->setup();

  my $is_interbase = $self->is_interbase();

  my $ref_name = $options{"-reference"} || die "-reference";
  my $single_gene_mode = $self->single_gene_mode();

  my $chr = cook_chromosome_name(
				 $ref_name,
				 "-genome" => $self->genome,
				 "-return-unknown" => 1) || die;
  my @usable;
  my %genes;

  my $map = $self->bucket_map()->{$chr};

  if ($map) {
    my @hits;

    my $start = $options{"-start"} || die "-start";
    my $end = $options{"-end"} || die "-end";

    $self->search_start($start);
    $self->search_end($end);

#    push @hits, @{$map->fuzzy_find("-site" => $start)};
#    push @hits, @{$map->fuzzy_find("-site" => $end)} unless $start == $end;
    # FAIL: doesn't query properly for large intervals!!
    if ($start == $end) {
      push @hits, @{$map->fuzzy_find("-site" => $start)};
    } else {
      die "end position $end is not greater than start position $start" unless $end > $start;
      push @hits, @{$map->fuzzy_find("-start" => $start, "-end" => $end)};
    }

    foreach my $hit (@hits) {
      next if $hit->{txEnd} < $start;
#      next if $hit->{txStart} > $end;

      my $spos = $hit->{txStart};
      $spos++ if $is_interbase;
      # - query coordinates are in-base
      # - refFlat coordinates are interbase and need correction
      # - GENE_EXON_REGION already uses in-base
      next if $spos > $end;

      push @usable, $hit;
    }
    %genes = map {$_->{name2}, 1} @usable;
  } else {
    printf STDERR "can't find map for $chr, ref=$ref_name\n";
    # OK for MT
  }

  if ($single_gene_mode) {
    my $glc = new GeneListCollapser();
    $glc->collapse("-hash" => \%genes, "-single" => 1);
  }

  $self->results_rows(\@usable);
  $self->results_genes([sort keys %genes]);

  foreach my $r (@usable) {
    die if $r->{txStart} > $r->{txEnd};
    # QC
  }

  #
  #  build genomically-sorted unique list:
  #
  my @sorted = sort {
    $a->{txStart} <=> $b->{txStart} ||
	$a->{txEnd} <=> $b->{txEnd} ||
	  $a->{name2} cmp $b->{name2}
	    # sort by coordinate, then symbol name
	    # e.g. chr13.56170652-63822116
	    # PRR20D,PRR20A,PRR20C,PRR20B,PRR20E,PCDH17,DIAPH3,TDRD3,MIR3169,PCDH20
	    # without this PRR20* genes can be reported in arbitrary order

  } @usable;


  my %saw;
  my @genes_genomic;
  # unique genes in interval, sorted genomically
  foreach my $r (@sorted) {
    my $g = $r->{name2} || dump_die($r, "no name2 data");
    next if $saw{$g};
#    printf STDERR "%d\n", $r->{txStart};
    push @genes_genomic, $g;
    $saw{$g} = 1;
  }
  $self->results_genes_genomic_order(\@genes_genomic);

  return scalar @usable;
}

sub is_valid_gene {
  my ($self, %options) = @_;
  $self->setup();
  my $g = $options{"-gene"} || die "-gene";
  return $self->all_genes()->{$g};
}

sub is_intronic {
  my ($self) = @_;
  my $rr = $self->results_rows() || die;
  my $result;

  my $start = $self->search_start() || die;
  my $end = $self->search_end() || die;

  my %intronic;
  my %not_intronic;

  my $VERBOSE = $self->intronic_debug() || 0;

  printf STDERR "starting search for %d-%d\n", $start, $end if $VERBOSE;

  foreach my $row (@{$rr}) {
    my $transcript_id = $row->{name} || die;
    my $exons = $row->{raw_rows} || die;

    if ($VERBOSE) {
      foreach my $exon (@{$exons}) {
	dump_die($exon, "exon", 1);
      }
    }

    if (@{$exons} > 1) {
      # skip for single-exon entries, e.g. NR_036199
      my $last_start;
      my $intronic;
      for (my $i=0; $i < @{$exons} -1; $i++) {
	die if $last_start and $last_start > $exons->[$i]->{txStart};
	$last_start = $exons->[$i]->{txStart};

	my $intron_start = ($exons->[$i]->{txEnd} || die) + 1;
	my $intron_end = ($exons->[$i + 1]->{txStart} || die) - 1;
	printf STDERR "intron for %s: %d-%d\n", $transcript_id, $intron_start, $intron_end if $VERBOSE;

	$intronic = 1 if $start >= $intron_start and $end <= $intron_end;
      }

      if ($intronic) {
	$intronic{$transcript_id} = 1;
      } else {
	$not_intronic{$transcript_id} = 1;
      }
    }
  }

  if (0) {
    die sprintf "some but not all intronic (intronic:%s not:%s)",
    join(",", sort keys %intronic),
    join(",", sort keys %not_intronic)
	if $VERBOSE and %intronic and %not_intronic;
  }
  
  return (%intronic and not(%not_intronic)) ? 1 : 0;
}

sub find_gene {
  #
  # find by gene symbol
  #
  my ($self, $gene) = @_;
  $self->setup();
  return $self->all_genes->{$gene};
}

sub find_accession {
  my ($self, $nm, %options) = @_;
  $self->setup();
  my $hits = $self->nm2rows->{$nm};
  if ($options{"-symbol"}) {
    my %genes;
    foreach my $hit (@{$hits}) {
      $genes{$hit->{name2} || die} = 1;
    }
    if (keys %genes == 1) {
      return (keys %genes)[0];
    } else {
      confess "ambiguous gene symbols in refflat!";
    }
  } else {
    return $hits;
  }
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
