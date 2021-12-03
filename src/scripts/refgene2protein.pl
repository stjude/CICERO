#!/bin/env perl
# extract features from blob of GenBank-formatted refGene records:
# ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.gbff.gz
# MNE 8/2014
#
# TO DO: swarm option with split preprocessed files?

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use Bio::SeqIO;

use MiscUtils qw(dump_die log_message);
use FileUtils qw(universal_open read_simple_file);
use File::Basename;
use Digest::MD5 qw(md5_hex);

use Reporter;
use DelimitedFile qw(df_bucket_by_header);
use EUtilities;

my %FLAGS;
my @FEATURES;

use constant MIN_AA_MISMATCH_FREQUENCY_FOR_MAJOR_MISMATCH => .04;
# copied from RefFlat2Genome.pm

use constant TRANSCRIPT_VARIANT_NUMBER_NON_NUMERIC => 9999;
# if transcipt number is not available, or alphanumeric

my $NUMERIC_TV_CONVERT_ALPHA = 1;
# for "transcript_variant_numeric", convert single-letter identifiers
# like a, b, c to coresponding 1, 2, 3, to preserve sortability.
# e.g. KRAS,
# https://www.ncbi.nlm.nih.gov/nuccore/NM_001369786.1

use constant REPORT_HEADERS => qw(
				   accession
				   version
				   gi
				   gene
				   protein
				   transcript_variant
				   transcript_variant_numeric
				   accession_versioned
				);

my @GBFF;
GetOptions(\%FLAGS,
	   "-rna-gbff=s" => \@GBFF,
	   # download from NCBI:
	   # ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.?.rna.gbff.gz
	   # where ? is 1-X
	   "-glob-gbff",
	   "-single-gbff=s",

	   "-out=s",

	   "-fetch-old-genbank=s",
	   # determine older transcript revisions needed and retrieve
	   # from NCBI

	   "-debug-eu=s",
	   "-preprocess-only",

	   "-generate-fasta",
	   "-cds-only",
	   # FASTA option

	   "-nm-only",
	  );

#push @GBFF, glob("*.gbff.*") if $FLAGS{"glob-gbff"};
push @GBFF, glob("*.gbff.gz") if $FLAGS{"glob-gbff"};
push @GBFF, $FLAGS{"single-gbff"} if $FLAGS{"single-gbff"};
# include suffix so outfile is excluded if present

if (@GBFF) {
  my $fasta_mode = $FLAGS{"generate-fasta"};

  my $of;
  if ($of = $FLAGS{out}) {
  } else {
    my $suffix = $fasta_mode ? ".fa" : ".refgene2protein.tab";
    if (@GBFF == 1) {
      $of = basename($GBFF[0]) . $suffix;
    } elsif ($fasta_mode) {
      $of = "refseq.fa";
    } else {
      $of = "refgene2protein.tab";
    }
  }

  if ($FLAGS{"generate-fasta"}) {
    generate_fasta("-in" => \@GBFF,
		   "-out" => $of);
  } else {
    generate_refgene2protein(
			     "-in" => \@GBFF,
			     "-out" => $of,
			    );
  }
} elsif (my $fn = $FLAGS{"fetch-old-genbank"}) {
  fetch_old_genbank($fn);
} else {
  die "?\n";
}

sub generate_refgene2protein {
  my %options = @_;
  my @infiles;
  my $in = $options{"-in"} || die "-in";
  my $outfile = $options{"-out"} || die "-out";
  if (-f $in) {
    @infiles = $in;
  } elsif (ref $in) {
    @infiles = @{$in};
  } else {
    die;
  }
  die unless @infiles;

  my $preprocess = $FLAGS{"single-gbff"} ? 0 : 1;

  if ($preprocess) {
    # parser is extremely slow, so digest input to just refGene records.
    # even parsing the ~50k NM_ accessions is very slow but MUCH
    # faster than having to parse everything.
    my $md5 = md5_hex(map {$_, -s $_} @infiles);
    my $f_cooked = sprintf "refgene.%s.txt.gz", $md5;
    my $f_split_template = sprintf "refgene_split.%s.txt.%%03d.gz", $md5;

    unless (-s $f_cooked) {
      my $wf = new WorkingFile($f_cooked, "-compress" => "gz");
      my $fh_out = $wf->output_filehandle();

      my $wf_chunk;
      my $fh_chunk;
      my $chunk_num = 0;
      my $nm_count = 0;
      my $chunk_size = 1000;

      my $nm_only = $FLAGS{"nm-only"};

      foreach my $infile (@infiles) {
	printf STDERR "preprocessing %s...\n", $infile;
	my $fh = universal_open($infile);
	my $usable;
	while (<$fh>) {
	  if (/^LOCUS\s+(\w+)/) {
#	    $usable = $1 =~ /^NM_/;
	    if ($nm_only) {
	      $usable = $1 =~ /^NM_/;
	    } else {
	      $usable = $1 =~ /^N[MR]_/;
	    }
	    # 12/2019: also include NR_ as these can have transcript
	    # variant numbers, which may be used by VEP+ to assist
	    # isoform filtering

	    if ($usable) {
	      if ($nm_count++ % $chunk_size == 0) {
		$wf_chunk->finish() if $wf_chunk;
		my $cfn = sprintf $f_split_template, $chunk_num++;
#		print STDERR " $cfn\n";
		$wf_chunk = new WorkingFile($cfn, "-compress" => "gz");
		$fh_chunk = $wf_chunk->output_filehandle();
	      }
#	      print "$nm_count\n";
	    }

	  }
	  if ($usable) {
	    print $fh_out $_;
	    print $fh_chunk $_;
	  }
	}
      }
      $wf->finish();
      $wf_chunk->finish();
    }

    @infiles = $f_cooked;
    if ($FLAGS{"preprocess-only"}) {
      printf STDERR "preprocessing finished\n";
      exit(0);
    }
  }

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [ REPORT_HEADERS ],
			 "-auto_qc" => 1,
			);

  my $skipped = 0;
  my $records = 0;
  my $ping = 1000;
  my $wrote = 0;
  foreach my $blob (@infiles) {
    log_message("parsing $blob...");
    my $fh = universal_open($blob);
    my $stream = Bio::SeqIO->new("-fh" => $fh,
				 -format => 'GenBank');
    while (my $record = $stream->next_seq()) {
      # Bio::SeqIO::genbank

      my $info = parse_one_genbank($record);
      $info->{protein} = "" unless exists $info->{protein};
      # won't be present for NR_

      my $acc = $info->{accession} || die;
      my $protein = $info->{protein} || "";
      # not present for NR_
      my $usable = $acc =~ /^N[MR]_/ ? 1 : 0;
      # exclude predicted XM_ records which don't appear in refFlat.txt
      die "refseq accession $acc has no protein sequence" if $acc =~ /^NM_/ and !$protein;
      # sanity check

      if ($usable) {
	$rpt->end_row($info);
	$wrote++;
      } else {
	$skipped++;
      }

      if (++$records % $ping == 0) {
	log_message(sprintf "raw:%d wrote:%d at:%s", $records, $wrote, $info->{accession});
      }

    }
  }
  $rpt->finish;

  printf STDERR "skipped %d records\n", $skipped;
}

sub fetch_old_genbank {
  my ($fn) = @_;
  my $ids = find_older_ids_to_query($fn);
  # find older transcript versions needed

  my $eu = new EUtilities();
  my $result_files = $eu->fetch_genbank("-ids" => $ids);
  # fetch old records from NCBI, and cache

  if (my $filter = $FLAGS{"debug-eu"}) {
    $result_files = [ grep {/$filter/} @{$result_files} ];
    die unless @{$result_files};
  }

  my $parsed = "refgene2protein_old.tab";
  if (-s $parsed) {
    # cache since this is slow
    printf STDERR "using cached %s\n", $parsed;
  } else {
    generate_refgene2protein(
			     "-in" => $result_files,
			     "-out" => $parsed
			    );
  }

  # compare refSeq (latest) with older versions to see which actually
  # have different AA annotations:
  my $latest = df_bucket_by_header("-file" => $fn, "-column" => "accession");
  my $older = df_bucket_by_header("-file" => $parsed, "-column" => "accession");

  my $rpt = new Reporter(
			 "-file" => "refgene2protein_merged.tab",
			 "-delimiter" => "\t",
			 "-labels" => [ REPORT_HEADERS, "note" ],
			 "-auto_qc" => 1,
			);

  foreach my $nm (sort keys %{$latest}) {
    my $latest_set = $latest->{$nm};
    die unless @{$latest_set} == 1;
    my $latest_ref = $latest_set->[0];
    $latest_ref->{note} = "latest";
    $rpt->end_row($latest_ref);

    if (my $older_set = $older->{$nm}) {
      # older revisions are available
      my %saw;
      my $latest_aa = $latest_ref->{protein};
      $saw{$latest_aa} = 1;

      my @sorted = sort {$b->{version} <=> $a->{version}} @{$older_set};
      # traverse versions backward from refSeq

      foreach my $ref (@sorted) {
	my $aa = $ref->{protein} || die;
	unless ($saw{$aa}) {
	  # ignore versions with AAs we've handled already
	  my $tag;
	  if (length($aa) == length($latest_aa)) {
	    # refFlat.txt mapping verification process already tolerates
	    # minor differences
	    my $disagreements = 0;
	    my $len = length $aa;
	    for (my $i = 0; $i < $len; $i++) {
	      $disagreements++ if substr($aa,$i,1) ne substr($latest_aa,$i,1);
	    }
	    my $fraction = $disagreements / $len;
	    if ($fraction >= MIN_AA_MISMATCH_FREQUENCY_FOR_MAJOR_MISMATCH) {
	      # significant difference, e.g. NM_000366, with significant
	      # changes late in the sequence
	      # in the last portion
	      $tag = sprintf "major difference (%.03f)", $fraction;
	    } else {
	      # minor disagreement
	      # keep everything: possibly try ALL revisions with different
	      # AAs in RefFlat2Genome.pm and keep the best one:
	      # a perfect mapping is better than one with 1-4% disagreement.
	      $tag = sprintf "minor difference (%.03f)", $fraction;
	    }
	  } else {
	    # always save: very likely to introduce compatibility problems
	    $tag = "length change";
	  }

	  $ref->{note} = $tag || die;
	  $rpt->end_row($ref);
	  $saw{$aa} = 1;
	}
      }
    }
  }
  $rpt->finish();

}

sub find_older_ids_to_query {
  my ($fn) = @_;

  my $df = new DelimitedFile(
			     "-file" => $fn,
			     "-headers" => 1,
			    );

  my @needed;
  while (my $row = $df->get_hash()) {
    my $acc = $row->{accession} || die;
    if ($acc =~ /^NM_/) {
      my $version = $row->{version} || die;
      for (my $v = $version - 1; $v >= 1; $v--) {
	push @needed, join ".", $acc, $v;
      }
    }
  }
  return \@needed;
}

sub parse_one_genbank {
  # Bio::SeqIO::genbank
  my ($record) = @_;
  my $verbose = 0;

  my $accession = $record->accession_number();
  my $version = $record->version();
  my $desc = $record->desc() || die "no desc for $accession";

  my $cds_count = 0;
#  my $transcript_variant = $desc =~ /transcript variant (\d+)/ ? $1 : 0;
#  my $transcript_variant = $desc =~ /transcript variant ([\w\-]+)/ ? $1 : 0;
  my $transcript_variant = $desc =~ /transcript variant ([^,]+)/ ? $1 : "";
  # - some genes e.g. KRAS use letters instead of numbers.
  # - sometimes these can be names, e.g. EGFRvIII
  # - may contain dashes, e.g. NM_001323302: transcript variant JNK1-a1
  # - or quote characters: ... (PTCH1), transcript variant 1a', mRNA

  my $transcript_variant_numeric;
  if ($transcript_variant =~ /^\d+$/) {
    # natively a number
    $transcript_variant_numeric = $transcript_variant;
  } elsif ($NUMERIC_TV_CONVERT_ALPHA and $transcript_variant =~ /^[A-Z]$/i) {
    # convert single letter ID to numeric equivalent starting at 1,
    # to preserve sortability of the set
    $transcript_variant_numeric = 1 + ord(uc($transcript_variant)) - ord("A");
  } else {
    $transcript_variant_numeric = TRANSCRIPT_VARIANT_NUMBER_NON_NUMERIC;
  }

  my %r;
  $r{accession} = $accession;
  $r{version} = $version;
  $r{accession_versioned} = join ".", $accession, $version;

#  $r{gi} = $record->primary_id;
  # fail: this is the GI of the refSeq record, not the protein gi
  $r{gi} = "";
  # protein gi might not be present
  $r{transcript_variant} = $transcript_variant;
  $r{transcript_variant_numeric} = $transcript_variant_numeric;

  print STDERR "\n" if $verbose > 1;
  printf STDERR "starting %s.%d...\n", $accession, $version if $verbose;

  my %genes;
  foreach my $feature ($record->get_SeqFeatures()) {
    printf STDERR "  %s\n", $feature->primary_tag() if $verbose > 1;
    my $ptag = $feature->primary_tag();
    if ($ptag eq "gene") {
      my @values = $feature->get_tag_values($ptag);
      die if @values > 1;
      $genes{$values[0]} = 1;
    } elsif ($ptag eq "CDS") {
      $cds_count++;
      foreach my $tag ($feature->get_all_tags()) {
	printf STDERR "tag:%s values:%s\n", $tag, join ",", $feature->get_tag_values($tag) if $verbose > 1;
	if ($tag eq "translation") {
	  my @values = $feature->get_tag_values($tag);
	  die unless @values == 1;
	  $r{protein} = $values[0];
	} elsif ($tag eq "db_xref") {
	  my @values = $feature->get_tag_values($tag);
	  foreach my $v (@values) {
	    my @f = split /:/, $v;
#	    die sprintf "unexpected field count in $v" unless @f == 2;
	    $r{gi} = $f[1] if $f[0] eq "GI";
	  }
	}
      }
    }
  }
  $r{gene} = join ",", sort keys %genes;

  die "no info!" unless %r;
#  die sprintf "CDS count %d for %s", $cds_count, $accession if $cds_count != 1;
  die "multiple CDS" if $cds_count > 1;
  # will need to return multiple records if this happens

  return \%r;
}

sub generate_fasta {
  my %options = @_;
  my @infiles;
  my $in = $options{"-in"} || die "-in";
  my $outfile = $options{"-out"} || die "-out";
  my $cds_only = $FLAGS{"cds-only"};

  if (-f $in) {
    @infiles = $in;
  } elsif (ref $in) {
    @infiles = @{$in};
  } else {
    die;
  }
  die unless @infiles;

  my $wf = new WorkingFile($outfile);
  my $fh_out = $wf->output_filehandle();

  my $records = 0;
  my $ping = 1000;
  my $wrote = 0;
  foreach my $blob (@infiles) {
    log_message("parsing $blob...");
    my $fh = universal_open($blob);
    my $stream = Bio::SeqIO->new("-fh" => $fh,
				 -format => 'GenBank');
    while (my $record = $stream->next_seq()) {
      # Bio::SeqIO::genbank

      my $info = parse_genbank_for_fasta($record);
      my $acc = $info->{accession} || die;
      my $gene = $info->{gene} || die;
      my $map = $info->{map};
      # some records don't have this, e.g. NR_037872.1

#      printf $fh_out ">%s %s %s\n%s\n", $acc, $gene, $info->{desc}, ($info->{seq} || die);
      # what was the problem with this format?  Maybe BLAST requires
      # a certain defline format, etc.?
#      printf $fh_out ">%s %s\n%s\n", $acc, $info->{desc}, ($info->{seq} || die);

      my $sequence = $info->{seq} || die "no sequence";
      if ($cds_only) {
	# extract coding region only
	my @cds = split /,/, $info->{CDS};
	die "multi CDS" if @cds > 1;
	my ($start, $end) = split /\-/, $cds[0];
	my $len = ($end - $start) + 1;
	my $cds = substr($sequence, $start - 1, $len);
	$sequence = $cds;
      }

      printf $fh_out ">%s %s /gene=%s",
	$acc, $info->{desc}, $gene;
      printf $fh_out " /map=%s", $map if $map;
      printf $fh_out "\n%s\n", $sequence;
      $wrote++;

      if (++$records % $ping == 0) {
	log_message(sprintf "raw:%d wrote:%d at:%s", $records, $wrote, $info->{accession});
      }
    }
  }
  $wf->finish;

}

sub parse_genbank_for_fasta {
  # Bio::SeqIO::genbank
  my ($record) = @_;
  my $verbose = 0;

  my $accession = $record->accession_number();
  my $version = $record->version();
  my $desc = $record->desc() || die "no desc for $accession";

  my %r;
  $r{accession} = join ".", $accession, $version;
  $r{seq} = $record->seq();
  $r{desc} = $desc;

  print STDERR "\n" if $verbose > 1;
  printf STDERR "starting %s.%d...\n", $accession, $version if $verbose;

  my %genes;
  my %map;
  my %cds;
  foreach my $feature ($record->get_SeqFeatures()) {
    printf STDERR "  %s\n", $feature->primary_tag() if $verbose > 1;
    my $ptag = $feature->primary_tag();
    if ($ptag eq "gene") {
      my @values = $feature->get_tag_values($ptag);
      die "multiple gene tag values" if @values > 1;
      $genes{$values[0]} = 1;
    } elsif ($ptag eq "source") {
      if ($feature->has_tag("map")) {
	my @values = $feature->get_tag_values("map");
	foreach my $v (@values) {
	  $map{$v} = 1;
	}
      }
    } elsif ($ptag eq "CDS") {
      my $range = join "-", $feature->start(), $feature->end();
      $cds{$range} = 1;
    }
  }
  $r{gene} = join ",", sort keys %genes;
  # horrible hack, are multiple symbols ever an issue for these records?
  $r{map} = join ",", sort keys %map;
  # similarly for mapping
  $r{CDS} = join ",", sort keys %cds;

  die "no info!" unless %r;
#  die sprintf "CDS count %d for %s", $cds_count, $accession if $cds_count != 1;
  return \%r;
}


__DATA__

=pod

=head1 INTRODUCTION

This script builds the "refgene2protein.tab" file used by
sv_inframe.pl.  It requires BioPerl and various modules
distributed with Cicero.

=head1 USAGE

1. create a new directory and cd into it.

2. download all "rna.gbff.gz" files from ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/.  As of this writing there are 7 files, for example human.1.rna.gbff.gz, human.2.rna.gbff.gz, etc.

3. run "refgene2protein.pl -glob-gbff".  This will parse the GenBank
records and create the refgene2protein.tab file.
