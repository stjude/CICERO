package RefFlatFile;

use strict;
use Exporter;
use Carp qw(confess cluck);

use Digest::MD5 qw(md5_hex);
use Configurable;

@RefFlatFile::ISA = qw(Configurable Exporter);
@RefFlatFile::EXPORT_OK = qw(uniquify_accessions);

use List::Util qw(min max);

use Bio::Tools::CodonTable;

use Reporter;
use MiscUtils qw(dump_die);
use FAI;
use GenomeUtils qw(complement);
use CacheManager;
use DebugRAM qw(debug_ram);
use Set::IntSpan;

use MethodMaker qw(
rows
max_entries_per_source
canonical_references_only
missing_genes_ok
preserve_interbase
idx_acc
idx_gene
cache_base_translations
translation_cache
cache_stats
cache_limit
strip_sharp_annotations
cm
unique_refflat_accessions
		  );

use constant REFGENE_FIELDS => qw(
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
# 16 fields

# use the same annotation labels as Cicero does
# (try to keep the vocabulary down):
use constant CICERO_INTRON => "intron";
use constant CICERO_CODING => "coding";
use constant CICERO_INTERGENIC => "intergenic";
use constant CICERO_UTR_5 => "5utr";
use constant CICERO_UTR_3 => "3utr";

@RefFlatFile::EXPORT_OK = qw(
CICERO_INTRON
CICERO_CODING
CICERO_INTERGENIC
CICERO_UTR_5
CICERO_UTR_3
);

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->cache_stats({qw(hits 0 misses 0)});
  $self->cache_limit(2000);
  # chews up a lot of RAM, flush periodically.
  # works much better if input data is sorted genomically or by transcript.
  $self->reset_cache();
  $self->rows([]);
  $self->configure(%options);
  return $self;
}

sub reset_cache {
  my ($self) = @_;
  $self->translation_cache({});
}

sub parse_file {
  my ($self, %options) = @_;
  my $type = $options{"-type"} || die "-type";
  my $file = $options{"-refflat"} || die "-refflat";
  # main
  my $t2g_helper = $options{"-index"};
  # helper file to map transcript IDs to gene symbols
  # (handling depends on format)
  my $max_entries = $self->max_entries_per_source();
  my $canonical_only = $self->canonical_references_only();
  my $missing_genes_ok = $self->missing_genes_ok();
  my $strip_sharp = $self->strip_sharp_annotations();

  my $is_refgene = $type eq "refgene";
  my $is_ensembl = $type eq "ensembl";
  my $is_aceview = $type eq "aceview";
  my $is_ucsc = $type eq "refflat";
  my $fuzzy_lookup = $options{"-fuzzy"};
  my $no_junctions = $options{"-no-junctions"};
  my $no_md5 = $options{"-no-md5"};
  my $skip_nr = $options{"-skip-nr"};
  my $lite = $options{"-lite"};

  my %priority = (
    "refgene" => 1,
    "ensembl" => 2,
    "refflat" => 3,
    "aceview" => 4
      );
  my $priority = $priority{$type} || die;

  debug_ram("refflat parse", "start") if $ENV{DEBUG_RAM};

  #
  #  map transcript IDs to genes:
  #
  my %t2g;
  if ($t2g_helper) {
    printf STDERR "parsing %s...\n", $t2g_helper;
    open(LOOKUP, $t2g_helper) || die;
    if ($is_refgene or $is_ensembl) {
      while (<LOOKUP>) {
	chomp;
	my ($t, $g) = split /\t/, $_;
	$t2g{$t} = $g;
      }
    } elsif ($is_ucsc) {
      while (<LOOKUP>) {
	chomp;
	my ($t, $g) = (split /\t/, $_)[0,4];
	$t =~ s/\.\d+$// if $fuzzy_lookup;
	$t2g{$t} = $g;
      }
    } else {
      die "unhandled $type";
    }
  }

  #
  #  parse refflat:
  #
#  printf STDERR "parsing %s...\n", $file;
  my $rows = $self->rows || die;
  open(RF, $file) || die "can't open $file";
  my $count = 0;
  my @headers = REFGENE_FIELDS;
  my $has_headers;
  while (<RF>) {
    chomp;
#    printf STDERR "line=%s\n", $_;
    if ($count == 0 and /^#/ and $is_refgene) {
      s/^#//;
      my @f = split /\t/, $_;
      @headers = @f;
      next;
    }
    $count++;
    if ($max_entries and $count > $max_entries) {
      printf STDERR "DEBUG, quitting at %d rows\n", $max_entries;
      last;
    }

    my @f = split /\t/, $_;
    my %r;
#    @r{REFGENE_FIELDS()} = @f;

    while (@headers > @f) {
      pop @headers;
    }

    @r{@headers} = @f;
    next if $skip_nr and $r{name} =~ /^NR_/;

    $r{priority} = $priority;

    if ($canonical_only) {
      my $chr = $r{chrom} || die;
      die unless $chr =~ /^chr(\w+)/;
      my $num = $1;
      unless ($num =~ /^\d+$/ or
	      $num eq "X" or
	      $num eq "Y" or
	      $num eq "M" or 
	      $num eq "MT") {
	next;
      }
    }

    if ($is_aceview) {
      $r{gene} = (split /\./, $r{name})[0];
    } elsif ($is_refgene) {
      if (exists $r{name2}) {
	$r{gene} = $r{name2};
      } elsif (exists $r{geneName}) {
	$r{gene} = $r{geneName};
      } elsif (@f == 11) {
	$r{gene} = $r{bin};
	# hack for refFlat.txt
      } else {
	die "can't identify gene field!";
      }
    } else {
      my $t = $r{name};
      $t =~ s/\.\d+$// if $fuzzy_lookup;
      $r{gene} = $t2g{$t};
    }

    unless ($r{gene}) {
      my $msg =  "no gene for " . $r{name};
      if ($missing_genes_ok) {
	printf STDERR "WARNING: %s\n", $msg;
      } else {
	die $msg;
      }
    }

    if ($strip_sharp) {
      # remove refFlat "sharp" annotations for records mapped
      # to multiple places
      $r{gene} =~ s/_loc([A-Z]|Par)$// if $r{gene} and $r{gene} !~ /^_loc/;
      $r{name} =~ s/\-loc([A-Z]|Par)$// if $r{name};
    }

    my @starts = split /,/, $r{exonStarts};
    my @ends = split /,/, $r{exonEnds};
    die unless @starts == @ends;
    die "exon count doesn't match starts for " . $r{gene} unless $r{exonCount} == scalar @starts;

    if ($lite) {
      # - tweak to save RAM for some use cases
      # - may break some features of this module!
      delete $r{exonStarts};
      delete $r{exonEnds};
      # for code that just uses parsed exons array
      delete $r{priority};
      delete $r{exonCount};
      delete $r{bin};
#      dump_die(\%r);
    }

    unless ($self->preserve_interbase()) {
      # coordinates are interbase:
      # convert starts to in-base
      foreach (@starts) {
	$_++;
      }
    }

    # generate exons:
    my @exons;
    for (my $i = 0; $i < @starts; $i++) {
      my %exon;
      $exon{start} = $starts[$i];
      $exon{end} = $ends[$i];
      push @exons, \%exon;
    }
    $r{exons} = \@exons;

    unless ($no_junctions) {
      # generate junctions:
      my @junctions;
      for (my $i = 0; $i < (@starts - 1); $i++) {
	my %junction;
	$junction{start} = $ends[$i];
	$junction{end} = $starts[$i + 1];
	push @junctions, \%junction;
      }
      $r{junctions} = \@junctions;
    }

    unless ($no_md5) {
      my @md5 = ($r{chrom}, $r{strand});
      push @md5, map {$_->{start}, $_->{end}} @exons;
      $r{md5} = md5_hex(@md5);
      # for duplicate-checking porpoises
    }
    push @{$rows}, \%r;
  }
  close RF;

  uniquify_accessions($rows) if $self->unique_refflat_accessions();

  debug_ram("refflat parse", "end") if $ENV{DEBUG_RAM};
}

sub get_reporter {
  my ($self, %options) = @_;
  my $rpt = new Reporter(
    %options,
    "-delimiter" => "\t",
    "-labels" => [ REFGENE_FIELDS ]
      );
  $rpt->headers_done(1);
  # no header line for this format
  return $rpt;
}

sub get_annotation_for_position {
  # given a refFlat format row and 1-based base number,
  # return CICERO-style region annotation (coding, intron, etc.)
  my ($self, %options) = @_;
  my $rf_row = $options{"-row"} || die "-row";
  my $base_num = $options{"-base"} || die "-base";
  my $extended = $options{"-extended"};
  my $intergenic_ok = $options{"-intergenic-ok"};

  my $tx_s = $rf_row->{txStart} || die;
  my $tx_e = $rf_row->{txEnd} || die;
  my $cds_s = $rf_row->{cdsStart} || die;
  my $cds_e = $rf_row->{cdsEnd} || die;
  my $strand = $rf_row->{strand} || die;

  if ($options{"-noncoding-synthesize-coding"} and
      $rf_row->{name} =~ /^NR_/ and
      $cds_s == $cds_e) {
    # fudge coding regions for non-coding exons so exonic regions
    # will be reported as such rather than 3' UTR
    # (make default?)
    $cds_s = $tx_s;
    $cds_e = $tx_e;
  }

  my $type;
  foreach ($tx_s, $cds_s) {
    $_++;
    # convert from interbase to in-base
  }

#  dump_die(\%options, "debug", 1);
#  dump_die($rf_row, "row debug", 1);

  my $feature_number = 0;
  my $exon_count;
  if ($base_num >= $tx_s and $base_num <= $tx_e) {
    # inside transcript (expected if we got this far)
    my @starts = split /,/, ($rf_row->{exonStarts} || die);
    my @ends = split /,/, ($rf_row->{exonEnds} || die);
    die unless @starts == @ends;
    $exon_count = scalar @starts;
    my $in_exon;
    for (my $i = 0; $i < @starts; $i++) {
      my $start = $starts[$i] + 1;
      # convert from interbase to in-base
      my $end = $ends[$i];
      $feature_number++;
      # exon/intron #

      if ($base_num >= $start and $base_num <= $end) {
	$in_exon = 1;
	last;
      }

      if ($i < @starts - 1) {
	my $is = $end + 1;
	my $ie = $starts[$i + 1];
	last if $base_num >= $is and $base_num <= $ie;
	# in intron: quit to record feature #
      }
    }

    if ($in_exon) {
      if ($base_num >= $cds_s and $base_num <= $cds_e) {
	# in CDS
	$type = CICERO_CODING;
      } else {
	# UTR
	if ($strand eq "+") {
	  $type = $base_num < $cds_s ? CICERO_UTR_5 : CICERO_UTR_3;
	} elsif ($strand eq "-") {
	  $type = $base_num > $cds_e ? CICERO_UTR_5 : CICERO_UTR_3;
	} else {
	  die;
	}
      }
    } else {
      $type = CICERO_INTRON;
    }
  } elsif ($intergenic_ok) {
    $type = CICERO_INTERGENIC;
  } else {
    dump_die($rf_row, "WTF: site not in tx");
  }

  if ($extended) {
    if ($strand eq "-") {
      $feature_number = $exon_count + 1 - $feature_number;
      $feature_number-- if $type eq CICERO_INTRON;
    }
    return ($type, $feature_number, $strand);
  } else {
    return $type;
  }
}

sub get_base_translations {
  # get base-by-base annotations for entire transcript between txStart/txEnd.
  # - category: utr5/utr3/exon/intron
  # - flag: is_coding
  # - exon/intron numbers (strand-adjusted)
  # useful for e.g. variant annotation.
  #
  # TO DO: track utr3/utr5 separately, and guarantee an exon entry
  # if it is partially/entirely noncoding?

  my ($self, %options) = @_;
  my $exons_only = $options{"-exons-only"};
  my $rf_row = $options{"-row"} || die "-row";
#  die "need preserve_interbase" unless $self->preserve_interbase;

  # build map of exonic bases
  # for entire txStart-txEnd:
  # - if exonic, coding if in CDS, otherwise UTR3/5 (strand)
  # - else intron
  # assign exon/intron #s as we go
  # - compare results with e.g. annovar exon #s, etc.
  # - run with same db as Annovar uses?

  my $cache_base_translations = $self->cache_base_translations();

  my $CACHE_NEW = 1;
  my $cache;
  my $cm;
  if ($CACHE_NEW) {
    $cm = $self->cm();
    unless ($cm) {
      $cm = new CacheManager("-cache_limit" => ($self->cache_limit || die));
      $self->cm($cm);
    }
  } else {
    $cache = $self->translation_cache();
  }

  my $acc = $rf_row->{name} || die;
#  printf STDERR "get_base_translations() for %s\n", $acc;

  my $chrom = $rf_row->{chrom} || die;
  my $tx_s = $rf_row->{txStart} || die;
  my $tx_e = $rf_row->{txEnd} || die;
  my $cds_s = $rf_row->{cdsStart} || die;
  my $cds_e = $rf_row->{cdsEnd} || die;
  my $strand = $rf_row->{strand} || die;

  my $cache_key;
  if ($cache_base_translations) {
    my $verbose = $ENV{RFF_CACHE_DEBUG};
    my @o;
    foreach (sort keys %options) {
      push @o, $_, $options{$_};
    }
    $cache_key = join ".", $acc, md5_hex(@o);
    # careful: an accession might be mapped to multiple locations/strands

    my $hit;
    if ($CACHE_NEW) {
      $hit = $cm->get($cache_key);
    } else {
      $hit = $cache->{$cache_key};
    }
    my $stats = $self->cache_stats();
    my $tag = $hit ? "hits" : "misses";
    printf STDERR "miss for $acc\n" if not($hit) and $verbose;
    $stats->{$tag}++;

    if ($CACHE_NEW) {
      my $cache_count = $cm->get_count();
      printf STDERR "RFF cache: entries:%d hits:%d misses:%d detail:%s\n", $cache_count, $stats->{hits}, $stats->{misses}, join " ", map {(split /\./, $_)[0]} sort $cm->get_keys() if $verbose;
    }  else {
      my $cache_count = scalar keys %{$cache};
      printf STDERR "RFF cache: entries:%d hits:%d misses:%d detail:%s\n", $cache_count, $stats->{hits}, $stats->{misses}, join " ", map {(split /\./, $_)[0]} sort keys %{$cache} if $verbose;
      $self->reset_cache() if $cache_count > $self->cache_limit();
    }
    # TO DO:
    # - smarter cache pruning?  maybe delete the longest-ago accessed entry
    #   using a sequence #?
    # - detect thrashing?
 
    return $hit if $hit;
  }

  my $type;
  foreach ($tx_s, $cds_s) {
    $_++;
    # convert from interbase to in-base
  }

  #
  # build map of exonic bases:
  #
  my @starts = split /,/, ($rf_row->{exonStarts} || die);
  my @ends = split /,/, ($rf_row->{exonEnds} || die);
  die unless @starts == @ends;
  my $exon_count = $rf_row->{exonCount};
  die unless $exon_count == @starts;

  my $exon_number = 0;
  my $intron_number = 0;

  my %exon;
  my %intron;
  my %splice_edges;

  for (my $i = 0; $i < @starts; $i++) {
    my $start = $starts[$i] + 1;
    # convert from interbase to in-base
    my $end = $ends[$i];

    die "txStart != 1st start" if $i == 0 and $start != $tx_s;

    $exon_number++;
    foreach ($start .. $end) {
      $exon{$_} = $exon_number;
    }
    
    $splice_edges{$start} = 1 unless $i == 0;
    $splice_edges{$end} = 1 unless $i == @starts - 1;

    unless ($exons_only) {
      if ($i < @starts - 1) {
	# build intron
	my $i_start = $end + 1;
	my $i_end = $starts[$i + 1];
	$intron_number++;
	foreach ($i_start .. $i_end) {
	  $intron{$_} = $intron_number;
	}
      }
    }
  }

  if ($strand eq "-") {
    #
    #  reverse exon/intron numbers for transcripts on -
    #
    my $esub = $exon_count + 1;
    foreach (values %exon) {
      $_ = $esub - $_;
    }
    foreach (values %intron) {
      $_ = $exon_count - $_;
    }
  }

  my @rows;

  for (my $base_num = $tx_s; $base_num <= $tx_e; $base_num++) {
    unless ($exons_only) {
      die unless exists($exon{$base_num}) or exists($intron{$base_num});
    }
    my $fnum;
    if ($fnum = $exon{$base_num}) {
      if ($base_num >= $cds_s and $base_num <= $cds_e) {
	# in CDS
	$type = CICERO_CODING;
      } else {
	# UTR
	if ($strand eq "+") {
	  $type = $base_num < $cds_s ? CICERO_UTR_5 : CICERO_UTR_3;
	} elsif ($strand eq "-") {
	  $type = $base_num > $cds_e ? CICERO_UTR_5 : CICERO_UTR_3;
	} else {
	  die;
	}
      }
    } elsif ($exons_only) {
      next;
    } else {
      $type = CICERO_INTRON;
      $fnum = $intron{$base_num};
    }
    die unless $fnum;

#    printf "%d: %s num:%d\n", $base_num, $type, $fnum;

    my %r;
    $r{base_number} = $base_num;
    $r{feature} = $type;
    $r{feature_number} = $fnum;
    $r{is_splice_edge} = 1 if $splice_edges{$base_num};
    push @rows, \%r;
  }

  if ($options{"-generate-codons"}) {
    my $fasta = $options{"-fasta"} || die "-fasta";
    my $fa = new FAI("-fasta" => $fasta);
    my $USE_CHUNK = 1;

    if ($USE_CHUNK) {
      my $cds_start = ($rf_row->{cdsStart} || die) + 1;
      my $cds_end = $rf_row->{cdsEnd} || die;
#      dump_die($rf_row, "cds ends before start: $cds_end $cds_start") if $cds_end < $cds_start;
      if ($cds_end < $cds_start) {
#	dump_die($rf_row, "cds ends before start: $cds_end $cds_start", 1);
	# can happen for NR_ entries, harmless
	my $result = $options{"-hash"} ? {} : [];
	if ($cache_base_translations) {
	  if ($CACHE_NEW) {
	    $cm->put($cache_key, $result);
	  } else {
	    $cache->{$cache_key} = $result;
	  }
	}
	return $result;
      }
      $fa->chunk_setup("-id" => $chrom,
		       "-start" => $cds_start,
		       "-end" => $cds_end);
    }

    my $strand = $rf_row->{strand} || die;
    my $direction;
    my $idx;
    if ($strand eq "+") {
      $idx = 0;
      $direction = 1;
    } elsif ($strand eq "-") {
      $direction = -1;
      $idx = $#rows;
    } else {
      die;
    }

    my @codon_queue;
    my $codon_number = 1;
    my $ct = new Bio::Tools::CodonTable();
#    my $max = scalar @rows;

    while (1) {
#      print "$idx\n";
      last if $idx < 0 or $idx >= @rows;
      if ($rows[$idx]->{feature} eq CICERO_CODING) {
#	dump_die($rows[$idx], "debug", 1);

	my $base;
	if ($USE_CHUNK) {
	  $base = $fa->get_chunked_base(
					"-start" => $rows[$idx]->{base_number},
					"-length" => 1
					);
	} else {
	  $base = $fa->get_chunk("-start" => $rows[$idx]->{base_number},
				  "-id" => $chrom,
				  "-length" => 1);
	}
	$base = complement($base) if $strand eq "-";

	push @codon_queue, $rows[$idx];
	$rows[$idx]->{codon_base} = $base;
	$rows[$idx]->{codon_number} = $codon_number;

	if (@codon_queue == 3) {
	  my $seq = join "", map {$_->{codon_base}} @codon_queue;
	  my $code = $ct->translate($seq) || die "can't translate $seq";
	  my $cn = 0;
	  foreach (@codon_queue) {
	    $_->{codon_base_number} = ++$cn;
	    $_->{codon_code} = $code;
#	    $_->{is_split_codon_exon_boundary} = 0;
	    # waste of memory
	  }

	  for (my $i = 0; $i < 2; $i++) {
	    my $p1 = $codon_queue[$i]->{base_number};
	    my $p2 = $codon_queue[$i + 1]->{base_number};
	    my $dist = abs($p1 - $p2);
	    if ($dist > 1) {
	      $codon_queue[$i]->{is_split_codon_exon_boundary} = 1;
	      $codon_queue[$i + 1]->{is_split_codon_exon_boundary} = 1;
	    }
	  }

	  $codon_number++;

	  @codon_queue = ();
	}
      }

      $idx += $direction;
    }
  }

  my $result;
  if ($options{"-hash"}) {
    my %hash;
    foreach my $r (@rows) {
      $hash{$r->{base_number} || die} = $r;
    }
    $result = \%hash;
  } else {
    $result = \@rows;
  }

  if ($cache_base_translations) {
    if ($CACHE_NEW) {
      $cm->put($cache_key, $result);
    } else {
      $cache->{$cache_key} = $result;
    }
  }
  return $result;
}

sub get_feature_summary {
  #
  # build summary of feature locations
  # 
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die;
  my $lines = $self->get_base_translations("-row" => $row);

  my %touched;
  foreach my $r (@{$lines}) {
    my $feature = $r->{feature} || die;
    my $bn = $r->{base_number} || die;
#    my $fno = $r->{feature_number} || dump_die($r, "no feature number for $accession");
    my $fno = $r->{feature_number};

    my $fname = $feature;
    $fname = "exon" if $fname eq "coding";

    $touched{$fname}{$fno}{$bn} = 1;
  }

  my @lines;
  foreach my $fname (sort keys %touched) {
    foreach my $fno (sort {$a <=> $b} keys %{$touched{$fname}}) {
      my @bases = keys %{$touched{$fname}{$fno}};

      my %r;
      $r{feature} = $fname;
      $r{feature_number} = $fno;
      $r{start} = min(@bases);
      $r{end} = max(@bases);

      push @lines, \%r;
    }
  }

  return \@lines;
}

sub find_by_accession {
  my ($self, $acc) = @_;
  my $idx = $self->idx_acc();
  unless ($idx) {
    my $rows = $self->rows() || die;
    $idx = $self->idx_acc({});
    foreach my $r (@{$rows}) {
      my $acc = $r->{name} || die;
      push @{$idx->{$acc}}, $r;
    }
  }
  $acc =~ s/\.\d+$//;
  # refflat files have unversioned accessions
  return $idx->{$acc};
}

sub get_rows {
  my ($self) = @_;
  return $self->rows;
}

sub find_by_gene {
  # TO DO: add GeneSymbolMapper, in case query symbol is older?
  my ($self, $gene) = @_;
  my $idx = $self->idx_gene();
  unless ($idx) {
    my $rows = $self->rows() || die;
    $idx = $self->idx_gene({});
    foreach my $r (@{$rows}) {
      my $gene = $r->{gene} || dump_die($r, "no gene field");
      push @{$idx->{$gene}}, $r;
    }
  }
  return $idx->{$gene};
}

sub find_gene_interval {
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die;
  my $buffer = $options{"-buffer"};
  die "-buffer" unless defined $buffer;
  my $single = $options{"-single"};
  confess "-single [0|1]" unless defined $single;
  my $buffer_upstream = $options{"-buffer-upstream"};
  my $unique_overlapped_exons = $options{"-unique-overlapped-exons"};

  my $wanted = $self->find_by_gene($gene);
  die "no rows for $gene" unless $wanted;

  if ($single) {
    # attempt to create a single interval spanning all mappings.
    # Good for creating a superset of isoform, however may not
    # be what you want if there are multiple mappings far away on
    # the same chrom.

    my %bucket;
    # bucket results by chrom and strand; some genes are mapped to
    # multiple chroms (e.g. DUX4), some to multiple strands within
    # a chrom (e.g. SSX2)
    foreach my $r (@{$wanted}) {
      my $key = join "_", @{$r}{qw(chrom strand)};
      push @{$bucket{$key}}, $r;
    }

    my @results;
    foreach my $key (sort keys %bucket) {
      # separate rows by chrom for e.g. DUX4
      my (@starts, @ends, $chr, $strand);
      my $set_exons = new Set::IntSpan();
      foreach my $r (@{$bucket{$key}}) {
	$chr = $r->{chrom} || die;
	$strand = $r->{strand} || die;
	# fixed for set

	if ($unique_overlapped_exons) {
	  my $exons = $r->{exons} || die;
	  foreach my $exon (@{$exons}) {
	    my $start = $exon->{start} || die;
	    my $end = $exon->{end} || die;
	    die unless $start < $end;
	    my $range = join "-", $start, $end;
	    printf STDERR "exon %s %s\n", $r->{name}, $range;
	    $set_exons = union $set_exons $range;
	  }
	} else {
	  push @starts, $r->{txStart} || die;
	  push @ends, $r->{txEnd} || die;
	}
      }

      my ($start, $end);
      if ($unique_overlapped_exons) {
#	printf STDERR "merged: %s %s\n", $chr, $set_exons;
	my @runs = split /,/, run_list $set_exons;
	foreach my $run (@runs) {
	  my @f = split /\-/, $run;
	  die unless @f == 2;
	  push @results, [ $chr, @f ];
	}
      } else {
	$start = min(@starts) - $buffer;
	$end = max(@ends) + $buffer;
	if ($buffer_upstream) {
#	printf STDERR "before %s\n", join ",", $gene, $chr, $start, $end, $strand;
	  if ($strand eq "+") {
	    # translation begins at genomic mapping start
	    $start -= $buffer_upstream;
	  } elsif ($strand eq "-") {
	    # translation begins at genomic mapping end
	    $end += $buffer_upstream;
	  }
#	printf STDERR "after %s\n", join ",", $gene, $chr, $start, $end, $strand;
	}
	push @results, [ $chr, $start, $end ];
      }
    }

    if (@results > 1) {
      if ($options{"-multi-ok"}) {
	return \@results;
      } else {
	die "multiple chroms for $gene, specify -multi-ok to return arrayref of results";
      }      
    } else {
      return @{$results[0]};
    }
  } else {
    #
    # separate interval for each mapping (safer)
    #
    # TO DO:
    # - skip perfect duplicates?  could happen if e.g. skipped exon
    # - merge overlaps?
    my @intervals;
    foreach my $r (@{$wanted}) {
      my $chrom = $r->{chrom} || die;
      my $start = $r->{txStart} || die;
      my $end = $r->{txEnd} || die;
      push @intervals, [$chrom, $start - $buffer, $end + $buffer];
      die "not implemented" if $buffer_upstream;
    }
    return \@intervals;
  }
}

sub get_strand_for_accession {
  my ($self, $acc) = @_;
  my $result;
  if (my $rows = $self->find_by_accession($acc)) {
    my %strand = map {$_->{strand}, 1} @{$rows};
    die unless scalar keys %strand == 1;
    ($result) = keys %strand;
  }
  return $result;
}

sub uniquify_accessions {
  # STATIC / exported
#  cluck "start uniquify";
  my ($rows) = @_;

  my %counts;
  foreach my $r (@{$rows}) {
    $counts{$r->{name}}++;
  }

  my $idx = 0;
  my %idx2letter = map {$idx++, $_} "A" .. "Z";
  # could done something with ord() but don't think that's PC these days

  my %counter;
  foreach my $r (@{$rows}) {
    my $nm = $r->{name};
    if ($counts{$nm} > 1) {
      my $count = $counter{$nm}++;
      my $new = sprintf '%s-loc%s', $nm, $idx2letter{$count} || die;
      printf STDERR "%s => %s\n", $nm, $new;
      $r->{name} = $new;
    }
    $counts{$r->{name}}++;
  }
}

sub get_genes {
  my ($self) = @_;
  $self->find_by_gene("bogus");
  # init index
  return [ keys %{$self->idx_gene} ];
}

sub get_neighbor_exon_distances {
  # given an intron position in a transcript, return distances
  # to nearest exons
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die;
  my $pos = $options{"-base"} || die;
  dump_die($row, "Debug", 1);

  my @es = split /,/, $row->{"exonStarts"} || die;
  foreach (@es) {
    $_++;
    # convert from interbase to in-base
  }
  my @ee = split /,/, $row->{"exonEnds"} || die;

  die unless @es == @ee;
  my $end = @es - 1;
  my %results;
  my $strand = $row->{strand} || die;

  for (my $i = 0; $i < $end; $i++) {
    # examine exon pairs
    die "out of order" unless $es[$i] < $ee[$i];

    my $end_this = $ee[$i];
    my $start_next = $es[$i + 1];

    if ($pos > $end_this and $pos < $start_next) {
      # found exons bounding intronic site

      my ($k1, $k2);
      if ($strand eq "+") {
	$k1 = "downstream";
	# downstream to next exon
	$k2 = "upstream";
	# upstream to previous exon
      } elsif ($strand eq "-") {
	$k1 = "upstream";
	$k2 = "downstream";
      } else {
	die;
      }

      $results{$k1} = {
	"boundary" => $start_next,
	"distance" => $start_next - $pos
      };

      $results{$k2} = {
	"boundary" => $end_this,
	"distance" => $pos - $end_this
      };

    }
#    printf STDERR "  %d-%d\n", $es[$i + 1], $ee[$i + 1];
  }
  die "can't find site" unless %results;

  return \%results;

}

1;
