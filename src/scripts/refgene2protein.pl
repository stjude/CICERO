#!/bin/env perl
# extract features from blob of GenBank-formatted refGene records:
# ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.rna.gbff.gz
# MNE 8/2014

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use Bio::SeqIO;

use MiscUtils qw(dump_die);
use FileUtils qw(universal_open read_simple_file);

use Reporter;
use DelimitedFile qw(df_bucket_by_header);
use EUtilities;

my %FLAGS;
my @FEATURES;

use constant MIN_AA_MISMATCH_FREQUENCY_FOR_MAJOR_MISMATCH => .04;
# copied from RefFlat2Genome.pm

use constant REPORT_HEADERS => qw(
				   accession
				   version
				   gi
				   gene
				   protein
				);

my @GBFF;
GetOptions(\%FLAGS,
	   "-rna-gbff=s" => \@GBFF,
	   # download from NCBI:
	   # ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.?.rna.bff.gz
	   # where ? is 1-X
	   "-glob-gbff",

	   "-out=s",

	   "-fetch-old-genbank=s",
	   # determine older transcript revisions needed and retrieve
	   # from NCBI

	   "-debug-eu=s",

	  );

push @GBFF, glob("*.gbff.*") if $FLAGS{"glob-gbff"};

if (@GBFF) {
  generate_refgene2protein(
			   "-in" => \@GBFF,
			   "-out" => ($FLAGS{out} || "refgene2protein.tab")
			  );
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

  my $rpt = new Reporter(
			 "-file" => $outfile,
			 "-delimiter" => "\t",
			 "-labels" => [ REPORT_HEADERS ],
			 "-auto_qc" => 1,
			);

  my $skipped_non_nm = 0;

  foreach my $blob (@infiles) {
    printf STDERR "parsing %s...\n", $blob;
    my $fh = universal_open($blob);
    my $stream = Bio::SeqIO->new("-fh" => $fh,
				 -format => 'GenBank');
    while (my $record = $stream->next_seq()) {
      # Bio::SeqIO::genbank
      my $info = parse_one_genbank($record);
      if ($info->{protein}) {
	# only protein-coding (i.e. not NR_)
	if ($info->{accession} =~ /NM_/) {
	  # exclude predicted XM_ records which don't appear in refFlat.txt
	  $rpt->end_row($info);
	} else {
	  $skipped_non_nm++;
	}
      }
    }
  }
  $rpt->finish;

  printf STDERR "skipped %d non-NM records\n", $skipped_non_nm;
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

  my $cds_count = 0;
  my %r;
  $r{accession} = $accession;
  $r{version} = $version;
#  $r{gi} = $record->primary_id;
  # fail: this is the GI of the refSeq record, not the protein gi
  $r{gi} = "";
  # protein gi might not be present

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
