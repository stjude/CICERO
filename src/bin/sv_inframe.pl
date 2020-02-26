#!/usr/bin/env perl
# determine whether fusion events are in-frame
# MNE 8/2014
#
# prerequisites:
#
#     module load blat
#     cbload cicero
#
# TO DO:
# - use config system for params! (warn for rnaseq though)
#
# - QC: find cases where geneA anchoring fails but CICERO says coding
#   - sync problems due to AA translation inconsistencies,
#     missing refSeqs, etc.
# - QC: independent 3-frame BLAT check (see "try #2" comment section below)
# - track transcript-specific annotations by hash and report by index
#   like other fields.  Truly "global" annotations are fewer.
#   also deal with duplicates, e.g. suppressed_nm.tab.frame.tab
# - ability to manually build up a list of files to process on the command
#   line (avoiding need for listfile)

use strict;
use warnings;

use 5.10.1;
# many 3rd-party modules installed only for this version

use Getopt::Long;
use Carp qw(confess cluck);

use File::Basename;
use List::Util qw(min max sum);
use Algorithm::Combinatorics qw(combinations);
use Bio::Tools::CodonTable;
use Set::IntSpan;

use MiscUtils qw(dump_die);
use DelimitedFile;
# use Reporter;
use GeneAnnotation;
use WorkingFile;
use BLATer;
use RefFlat2Genome;
use FileUtils qw(read_simple_file);
use CommandLineRebuilder;
use Cluster;
use GenomeUtils qw(reverse_complement);
use AnnotationField;
use RefFlatFile;
use TdtConfig;
use Counter;

#my $MIN_AA_FOR_GENE_A = 10;
#my $MIN_AA_FOR_GENE_B = 10;
my $MIN_TUPLE_AA = 10;
# insufficient, see multiple_matching_frames.tab
# tricky case, intra-chrom fusion, how to fix??

my $GENE_FLANK_SEARCH_NT = 1000;

my $BLAT_AA_MIN_IDENTICAL = 10;
my $BLAT_AA_MAX_MISMATCHES = 1;
my $BLAT_AA_MAX_GAPS = 0;

my $CODON_SEEK_MAX_TRIES = 10;
my $CODON_SEEK_MAX_TRIES_BEYOND_TRANSCRIPT = 1000;

my $MIN_AA_DOWNSTREAM_OF_M_IN_UTR_CHECK = 10;

my $CLUSTER_RAM = 10000;

my $UTR_TUPLE_SEARCH_WIGGLE_CODONS = 3;
# set to 0 to disable

my $NOT_INFRAME_IF_ENTIRELY_UTR3 = 1;
# if set, never consider in-frame if the search sequence upstream of
# PosA or downstream of PosB is entirely in the 3' UTR

my $REANNOTATE_INFRAME_GENES = 1;
# for in-frame results, reannotate geneA and geneB columns
# if required

#my $ENABLE_SYNTHETIC_AA_INTRON = 1;
my $ENABLE_SYNTHETIC_AA_INTRON = 0;
# unfinished, disable

use constant SV_MATCH_CODE_FRAMESHIFT => 0;
use constant SV_MATCH_CODE_INFRAME => 1;
# GeneA and GeneB are in-frame
use constant SV_MATCH_CODE_UTR5_GENEB_COMPLETE_CDS => 2;
# GeneB portion contains canonical CDS start
use constant SV_MATCH_CODE_POTENTIAL_UTR5_FUSION => 3;

use constant STOP_CODON_CODE => "*";

use constant CICERO_FEATURE_INTERGENIC => "intergenic";

my $VERBOSE = 0;

my $FULL_LENGTH_MIN_GENE_A_ANCHOR = 10;
#my $FULL_LENGTH_MIN_BLAT_IDENTITY_TO_MERGE = 0.93;
my $FULL_LENGTH_MIN_BLAT_IDENTITY_TO_MERGE = 0.95;
my $FULL_LENGTH_START_COMPARE_LENGTH = 10;

my $FULL_LENGTH_DIGEST_EXCLUDE_INTERSTITIAL_STOP = 1;
# when creating summaries of full-length sequences, discard those
# predicted to have an interstitial stop codon between geneA and geneB.
# Note that this may eliminate ALL results for some fusions.

my $PRESORT = 1;
my $PROGRESS_PING = 10;
use constant DEFAULT_CHROM_CACHE_LIMIT => 2;
# only if pre-sort enabled
# a value of 2 likely quickly prunes still-needed chroms so there
# will be some repetitive re-loading, however process memory usage
# will be kept down which is one of our goals

my $CACHE_TRANSCRIPT_MAPPINGS = 1;
use constant DEFAULT_TRANSCRIPT_CACHE_LIMIT => 25;

my %KNOWN_FEATURES = map {$_, 1} qw(
                                     intron
                                     coding
                                     intergenic
                                     5utr
                                     3utr
                                  );

my @params = (
              "-single=s",
              "-list=s",

              "-genome=s",

              "-refflat=s",
              # refFlat.txt used in CICERO run, e.g.
              # /nfs_exports/genomes/1/Homo_sapiens/Hg19/mRNA/RefSeq/refFlat.txt
              "-fasta=s",

              "-refgene2protein=s",
              # mapping of refseq accession to protein,
              # see refgene2protein.pl

              "-save-tempfiles",

              "-test-single=s",
              "-test-all",

              "-verbose-r2g",
              "-verbose" => \$VERBOSE,

              "-sj-tag-outfile",

              "-extract-example=s",

              "-cache",
              # don't overwrite outfiles if already exist
              "-force",

              "-no-blat",

              "-ignore-duplicate-outfile",
              # debug only

              "-cache-genome-mappings=i" => \$CACHE_TRANSCRIPT_MAPPINGS,

              "-hack-extract-inframe",

              "-force-refflat-map",
              # accept refFlat.txt mappings even if they disagree
              # with refSeq AA.  DEBUGGING ONLY!

              "-utr-frame-tests",

              "-fq",
              "-cluster",

              "-restrict-a-nm=s",
              "-restrict-b-nm=s",
              # debug
              "-restrict-junction=s",

              "-junction-mode",

              "-infer-aa",
              "-blat-min-score=i",
              "-junction-flanking=i",

              "-not-inframe-if-utr3=i" => \$NOT_INFRAME_IF_ENTIRELY_UTR3,

              "-report-coding-base-number",

              "-patch-coding-base-number",

              "-hack2",

              "-generate-qpos",

              "-min-aa-tuple=i" => \$MIN_TUPLE_AA,
              "-min-aa-blat=i" => \$BLAT_AA_MIN_IDENTICAL,

	      "-enable-full-length",


	      "-digest-full-length=s",
	      "-prep-excerpt=s",
	      "-run-digest",
	      "-fusions=s",

	      "-frame-debug=s",

	      "-crest",

	      "-cache-whole-chrom=i",
	      "-limit-chrom-cache=i",
	      "-limit-transcript-cache=i",

	      "-presort=i" => \$PRESORT,
	      "-ping=i" => \$PROGRESS_PING,
	     );

my %FLAGS;
GetOptions(\%FLAGS, @params);

$FLAGS{"report-coding-base-number"} = 1;
# 12/2015: revisiting code for other purposes so also enable this

my $JUNCTION_FLANKING_SEQUENCE = $FLAGS{"junction-flanking"} || 60;
# when generating contig sequence from a junction, # of nt to
# extract on either side the junction.
# 60 = 20 codons of AA for anchoring on each side; enough?

my $BLAT_MIN_SCORE = exists $FLAGS{"blat-min-score"} ? $FLAGS{"blat-min-score"} : 20;
# blat -minScore override of 20 required to detect 10-AA hit for GSHSMRYFYT in
#
# >peptide_frame_geneA
# RAGLSPSSPPGSHSMRYFYTA
#
# vs.
#
# >NM_005514
# MLVMAPRTVLLLLSAALALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDRTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEPSSQSTVPIVGIVAGLAVLAVVVIGAVVAAVMCRRKSSGGKGGSYSQAACSDSAQGSDVSLTA

my $clr = new CommandLineRebuilder("-parameters" => \@params, "-flags" => \%FLAGS);
#$clr->exclude_parameter("-single");
$clr->exclude_parameter("-list");
$clr->exclude_parameter("-cluster");


if (my $acc = $FLAGS{"test-single"}) {
  test_single($acc);
  exit(0);
} elsif ($FLAGS{"test-all"}) {
  test_all();
  exit(0);
} elsif ($FLAGS{"extract-example"}) {
  extract_example();
  exit(0);
} elsif ($FLAGS{"hack-extract-inframe"}) {
  extract_inframe_hack();
  exit(0);
} elsif ($FLAGS{hack2}){
  hack2();
  exit(0);
} elsif ($FLAGS{"utr-frame-tests"}) {
  utr_frame_tests();
  exit(0);
} elsif ($FLAGS{cluster}) {
  # submit jobs to cluster
  cluster_jobs();
  exit(0);
} elsif ($FLAGS{"digest-full-length"}) {
  digest_full_length();
  exit(0);
} elsif ($FLAGS{"prep-excerpt"}) {
  prep_excerpt_files();
  exit(0);
} elsif ($FLAGS{"frame-debug"}) {
  frame_debug();
  exit(0);
}

if ($FLAGS{"patch-coding-base-number"}) {
  patch_coding_base_numbers();
} elsif ($FLAGS{"single"} or $FLAGS{"list"}) {
  process_svs();
} else {
  die "specify -single [CICERO/CREST file] or -list [file containing list of filenames]\n";
}

sub get_infiles {
  my $infiles;
  if (my $fn = $FLAGS{"single"}) {
    $infiles = [ $fn ];
  } elsif (my $lf = $FLAGS{"list"}) {
    $infiles = read_simple_file($lf);
  } else {
    die "specify -single FILE or -list LISTFILE";
  }
  return $infiles;
}

sub process_svs {
  printf STDERR "configuration:\n";
  printf STDERR "  pre-sorting data: ";
  if ($PRESORT) {
    $FLAGS{"limit-chrom-cache"} = DEFAULT_CHROM_CACHE_LIMIT unless exists $FLAGS{"limit-chrom-cache"};
    printf STDERR "yes\n";
    printf STDERR "    cache whole chromosomes: %s\n", $FLAGS{"cache-whole-chrom"} ? "yes" : "no";
    printf STDERR "    chromosome cache limit: %d\n", $FLAGS{"limit-chrom-cache"};
  } else {
    print STDERR "no\n";
  }

  print STDERR "  cache transcript mappings?: ";
  if ($CACHE_TRANSCRIPT_MAPPINGS) {
    $FLAGS{"limit-transcript-cache"} = DEFAULT_TRANSCRIPT_CACHE_LIMIT unless exists $FLAGS{"limit-transcript-cache"};
    printf STDERR "yes (accession cache limit=%d)\n", $FLAGS{"limit-transcript-cache"};
  } else {
    print STDERR "no\n";
  }

  my $infiles = get_infiles();

  my $r2g = get_r2g();

  my $refgene_fn = get_refflat_file();
  my $infer_aa = $FLAGS{"infer-aa"};

  my $ga = new GeneAnnotation(
                                "-style" => "refgene_flatfile",
                                "-refgene_flatfile" => $refgene_fn,
#                              "-verbose" => 1
                               );
  # quickly find transcript IDs from regions

  if ($infer_aa) {
    $ga->ignore_non_coding(0);
    # required for e.g. Ensembl transcripts, where we don't
    # have the AA for comparison
  }

  my $junction_mode = $FLAGS{"junction-mode"};
  if ($junction_mode and $PRESORT) {
    printf STDERR "junction mode: disabling presort option\n";
    $PRESORT = 0;
  }

  my @sv_extra_fields;
  if ($junction_mode) {
    # add columns for generated CICERO-format fields
    push @sv_extra_fields, qw(
                               chrA
                               posA
                               ortA
                               featureA

                               chrB
                               posB
                               ortB
                               featureB

                               sv_ort
                               contig
                           );
  }

  push @sv_extra_fields, (
                         "sv_any_inframe",
                         "sv_any_inframe_canonical",

                         "sv_inframe",
                         #                                           "sv_inframe_aa",
 # AA predicted from contig sequence, in-frame for GeneA transcript

                         "sv_refseqA",
                         "sv_refseqA_codon",
                         "sv_refseqA_exon",
                         "sv_refseqA_anchor_type",

                         "sv_refseqB",
                         "sv_refseqB_codon",
                         "sv_refseqB_exon",
                         "sv_refseqB_anchor_type",
                         "sv_posB_adjusted",
                         "sv_AA",
                         "sv_desc",

                         "sv_processing_exception",
                         "sv_general_info",
                        );

  if ($junction_mode) {
    push @sv_extra_fields, (
                            "sv_refseqA_inframe",
#                            "sv_refseqA_inframe_transcripts",
                            "sv_refseqB_inframe",
#                            "sv_refseqB_inframe_transcripts",
                           );
    # on second thought don't report the transcripts for now:
    # if geneA status is inferred from results, the list may not
    # be comprehensive, because it relies on geneB being synced first.
  }

  if ($junction_mode) {
    # CICERO input already has these annotations.
    # leave separate code block, may want to optionally add in the
    # future for comparisons
    push @sv_extra_fields, (
                            "sv_refseqA_annot",
                            "sv_refseqB_annot",
                           );
  }

  if ($FLAGS{"report-coding-base-number"}) {
    push @sv_extra_fields, "sv_refseqA_coding_base_number";
    push @sv_extra_fields, "sv_refseqB_coding_base_number";
    push @sv_extra_fields, "sv_refseqB_last_coding_base_number";
  }

  push @sv_extra_fields, (
                          "sv_refseqA_AA_index",
                          # last codon in geneA
                          "sv_refseqB_AA_index",
                          # first codon in geneB
                          "sv_refseqA_contig_index",
                          # contig index of last base of last codon of
                          # geneA anchor
                          "sv_refseqB_contig_index",
                          # contig index of 1st base of 1st codon of
                          # geneB anchor
                          "sv_interstitial_AA",
                          # AA sequence between geneA and geneB mapping
                          "sv_frame_index",
                         );
  # 2/2016

  unique_outfile_check($infiles);

  my $enable_full_length = $FLAGS{"enable-full-length"};

  foreach my $cicero (@{$infiles}) {
    printf STDERR "%s: processing file %s...\n", scalar(localtime()), $cicero;
    my $start_time = time;
    my $df = new DelimitedFile(
                               "-file" => $cicero,
                               "-headers" => 1,
                              );
    my $outfile = get_outfile($cicero);

    if ($FLAGS{cache} and -s $outfile and not($FLAGS{force})) {
      printf STDERR "skipping %s, exists\n", $outfile;
      next;
    }


    my $rpt = $df->get_reporter(
                                "-file" => $outfile,
                                "-extra" => \@sv_extra_fields
                               );
    $rpt->auto_qc(1);

    my $rpt_full_length;
    my %nm2protein;
    my $f_out_full_length = $outfile . ".full_length.tab";
    if ($enable_full_length) {
      my $f_r2p = get_r2p_file();
      my $df = new DelimitedFile("-file" => $f_r2p,
				 "-headers" => 1,
				);
      while (my $r = $df->get_hash()) {
	$nm2protein{$r->{accession}} = $r->{protein};
      }

      $rpt_full_length = new Reporter(
				      "-file" => $f_out_full_length,
				      "-delimiter" => "\t",
				      "-labels" => [
						    "sample",
						    "geneA",
						    "geneB",
						    "gene_a_nm",
						    "gene_b_nm",
						    "contig",
						    "contig_frame",

						    "AA_gene_a_full",
						    "AA_gene_b_full",

						    "AA_gene_a",
						    "AA_interstitial",
						    "AA_gene_b",
						    "QC_span",
						    "QC_pass",
						    "full_length_sequence",
						    "full_length_sequence_length",
						    "full_length_sequence_contains_entire_gene_b",
						    "notes",
						   ]
				     );
    }

    # while (my $row = $df->next("-ref" => 1)) {  # headerless

    my @rows_raw;
    while (my $row = $df->get_hash()) {
      push @rows_raw, $row;
    }
    my @rows_process = @rows_raw;
    my $c = new Counter(\@rows_raw, "-mod" => $PROGRESS_PING);

    if ($PRESORT) {
      foreach my $h (qw(chrA chrB geneA geneB)) {
	dump_die($rows_raw[0], "where is $h") unless exists $rows_raw[0]->{$h};
      }

      @rows_process = sort {$a->{chrA} cmp $b->{chrA} ||
			      $a->{chrB} cmp $b->{chrB} ||
				$a->{geneA} cmp $b->{geneA} ||
				  $a->{geneB} cmp $b->{geneB}}
	@rows_raw;

#      foreach my $r (@rows_process) {
#	printf "%s\n", join " ", @{$r}{qw(chrA chrB geneA geneB)};
#      }
    }

    foreach my $row (@rows_process) {
      $c->next(join ".", @{$row}{qw(chrA chrB geneA geneB)});
      if ($FLAGS{crest}) {
	# experiment: attempt to use CREST-format input
	$row->{qposA} = ($row->{posA} || dump_die($row, "no CREST qposA field"));
	# in CICERO, qposA/qposB are the breakpoint locations in the CONTIG.
	# 3346068


	$row->{featureA} = $row->{featureB} = "coding";
	# HACK

	# need:
	# annotateA / annotateB (just used as backup if featureA empty?)

      }

      my $junction_peptide;
      if ($junction_mode) {
        if (my $r = $FLAGS{"restrict-junction"}) {
          next unless $row->{junction} eq $r;
        }

        junction_setup($row, $ga, $r2g);
        $junction_peptide = $row->{junction_peptide};
      }

      dump_die($row, "starting processing", 1) if $VERBOSE;

      my $contig = get_contig($row);
      #  my $gene_a = $row->{GeneA} || die;
      #  my $gene_b = $row->{GeneB} || die;
      # these may NOT be reliable for our purposes:
      # e.g. CICERO output refers to KMT2A but refFlat file refers to MLL.

      my $annot_global = new AnnotationField("-unique" => 1);
      my $annot_trouble = new AnnotationField("-unique" => 1);
      # processing problem(s) that might affect the results

      my $is_cicero = is_cicero($row);

      if ($is_cicero) {
        my $orientation = $row->{sv_ort} || die;
        if ($orientation ne ">") {
          foreach my $f (@sv_extra_fields) {
            $row->{$f} = "" unless defined $row->{$f};
          }
#          dump_die($row);
          $row->{sv_processing_exception} = "unhandled orientation $orientation";
#          $rpt->end_row($row);
          next;
          # hacktacular: we "should" if/then block this but indents are
          # already severe  :(
        }
      }

      my $frames = get_frames($contig);
      # possible AA frames for this sequence

      my $frames_masked = mask_frames_for_gene_a(
                                                 "-row" => $row,
                                                 "-frames" => $frames
                                                );

      my $b_nms = find_nms_from_ga(
                                   "-ga" => $ga,
                                   "-row" => $row,
                                   "-which" => "B"
                                  );

      # find all transcripts touched by breakpoint B
      if (my $restrict = $FLAGS{"restrict-b-nm"}) {
        die "restrict-b-nm failed" unless grep {$_ eq $restrict} @{$b_nms};
        $b_nms = [ $restrict ];
      }

      my @pairs;
      # pairs of comparisons of GeneB transcripts vs. GeneA transcripts
      my @geneb_utr_matches;

      if (@{$b_nms}) {
        $annot_global->add("geneB_transcripts", $b_nms);
      } else {
        $annot_trouble->add("geneB_no_transcripts");
      }

      my $b_chrom = get_chrom($row, "B") || die;
      my $b_bn = get_base_number($row, "B") || die;
      my $b_ort = get_orientation($row, "B") || die;

      # my $b_seek_dir;
      # # seek direction to go downstream in the transcript mapping for GeneB
      # if ($b_ort eq "+") {
      #         $b_seek_dir = 1;
      # } elsif ($b_ort eq "-") {
      #         $b_seek_dir = -1;
      # } else {
      #         die;
      # }
      my $b_feature = get_feature($row, "B") || die "can't find feature B";

      my $b_global_inframe = "";
      my %b_global_inframe_transcripts;
      # separate global reporting of whether any transcript is in-frame.
      # not sure if this will be good enough; have to see.
      # can't simply add a comparison field to @pairs below because
      # evaluation of geneA is done even when geneB is not in frame,
      # which in the current workflow stops processing altogether.

      #
      #  Anchor gene B first, then check mapping to A.
      #  This is done because some events only require a mapping to B.
      #
      foreach my $b_nm (@{$b_nms}) {
        #
        # find the codon containing the breakpoint site
        # - for geneA, breakpoint is the last mapped nt before the break
        # - however site may not be perfectly accurate (blat)
        #

        printf STDERR "new B: %s\n", $b_nm if $VERBOSE;
        my $set_passed = get_transcript_mappings(
                                                 "-end" => "B",
                                                 "-accession" => $b_nm,
                                                 "-row" => $row,
                                                 "-r2g" => $r2g,
                                                 "-annot-global" => $annot_global,
                                                 "-annot-trouble" => $annot_trouble
                                                );

        #
        #  we now know which transcript mapping touches the GeneB breakpoint
        #  and exactly where the breakpoint falls in the mapping
        #
        my %gene_b_inframe;
        my $processing_exception;
        my %frame_blat;

        if (@{$set_passed}) {
          my ($hit, $codon_index, $b_codon_number, $b_final_pos, $b_coding_basenum) = @{$set_passed->[0]};

#          dump_die($hit);

          my $b_aa = $r2g->get_aa("-accession" => $b_nm);
          printf STDERR "GeneB AA: %s\n", ($b_aa || die) if $VERBOSE;

          #
          #  - perform single BLAT of geneB AA vs all 3 frames
          #  - cache results by subject sequence and use in frame_search()
          #  - this is much more efficient than running an
          #    individual BLAT for each frame.
          #  - results only used if tuple lookup fails (e.g. CDS mismatch)
          #
          my $blat = get_blat();
          $blat->verbose(1) if $VERBOSE;

          my $parser = $blat->blat(
                                   "-query" => {
                                                $b_nm => ($r2g->get_aa("-accession" => $b_nm) || die),

                                               },
#                                   "-database" => $frames,
                                   "-database" => $frames_masked,
                                   "-protein" => 1,
                                  );
          my $result = $parser->next_result;
          # one object per query sequence (only one query seq)
          if ($result) {
            while (my $hit = $result->next_hit()) {
              # hits from this query to a database sequence
              # (Bio::Search::Hit::HitI)
              push @{$frame_blat{$hit->name()}}, $hit;
            }
          }

          if ($junction_peptide) {
            my $junction_peptide_found;
            foreach my $fid (keys %{$frames}) {
              my $frame = $frames->{$fid};
              if (index($frame, $junction_peptide) != -1) {
                $junction_peptide_found = 1;
              }
            }

            $annot_trouble->add("ERROR_no_frames_match_user_peptide")
              unless ($junction_peptide_found);
            # serious problem: why no match?  transcript sync issue?
          }

          foreach my $fid (keys %{$frames}) {
            my $frame = $frames->{$fid};

            next if $junction_peptide and
              index($frame, $junction_peptide) == -1;
            # match required if peptide specified

            my ($idx, $fragment, $frame_error, $utr_codons, $blat_mode) = frame_search(
                                                "-row" => $row,
                                                "-r2g" => $r2g,
                                                "-frame" => $frame,
                                                "-hit" => $hit,
                                                "-index" => $codon_index,
                                                "-length" => $MIN_TUPLE_AA,
                                                "-direction" => 1,
                                                "-blat" => $frame_blat{$fid},
                                                  "-annot-global" => $annot_global,
                                                "-annot-trouble" => $annot_trouble,
                                                "-end" => "B",
                                               );
            # find unique anchoring in AA.  The fragment size will
            # be expanded for uniqueness if necessary.
            # the formal codon number is taken from the RefFlat2Genome mapping,
            # however the string lookup is required to separate the
            # GeneA portion from the GeneB portion for lookup for GeneB.
            #
            # the fragment will NOT include the codon at the breakpoint
            # as this might not match Cicero's contig sequence (see
            # comments in get_aa_fragment()), and will be 1 AA away in the
            # appropriate direction.

	    my $fragment_extended = $fragment;

	    if ($utr_codons) {
	      $annot_global->add("geneb_utr_codons", $utr_codons);
	      # debug for evaluation/testing

	      my ($idx2, $fragment2, $frame_error2, $utr_codons2, $blat_mode2) = frame_search(
						"-row" => $row,
						"-r2g" => $r2g,
						"-frame" => $frame,
						"-hit" => $hit,
						"-index" => $codon_index,
						"-length" => $MIN_TUPLE_AA + abs($b_codon_number) + 20,
						"-direction" => 1,
						"-blat" => $frame_blat{$fid},
  					        "-annot-global" => $annot_global,
						"-annot-trouble" => $annot_trouble,
						"-end" => "B",
					       );
	      $fragment_extended = $fragment2;
	      # find a longer anchor that definitely extends into the
	      # start of geneB.  This is critical for generating full-length
	      # sequences where geneB lands in 5' UTR.
	    }

	    if ($frame_error) {
	      # special case: no upstream fragment!  e.g. in 1st codon,
	      # how to handle??
	      $annot_trouble->add($frame_error, $b_nm);
	      # will be reported multiple times; good example for
	      # converting to hash model for this field
	    } elsif ($idx != -1) {
	      #
	      #  this frame is in-frame with GeneB
	      #
	      $b_global_inframe = 1;
	      $b_global_inframe_transcripts{$b_nm} = 1;
	      # global flag indicating at least one transcript for
	      # geneB is in-frame.

              unless ($blat_mode) {
                if (length($fragment) > $MIN_TUPLE_AA) {
                  $annot_global->add("geneB_expanded_tuple_required", $fragment);
                  printf STDERR "final expanded fragment: %s\n", $fragment;
                } elsif (length($fragment) < $MIN_TUPLE_AA) {
                  $annot_global->add("geneB_truncated_tuple_required", $fragment);
                  printf STDERR "final truncated fragment: %s\n", $fragment;
                }
              }

              my $frame_raw = $frames->{$fid};

              my ($upstream, $upstream_full, $downstream);

              if ($idx == 0) {
                # match is at start of sequence (!)
                # mapping site ambiguity?
                # ambiguous sequence?
                $upstream = $upstream_full = "";
                $downstream = $frame_raw;
              } else {
                # get portion of AA contig involving geneA:
                #
#                $upstream = substr($frame_raw, 0, $idx - 1);
                # -1: don't include codon at breakpoint which may touch
                # both transcripts and be altered, which may interfere
                # with mapping to geneA.
                $upstream = substr($frame_raw, 0, $idx);
                # revise:
                # - skipping a base introduces a gap for perfect SVs
                # - BLAT is used as a backup
                $upstream_full = substr($frame_raw, 0, $idx);
                $downstream = substr($frame_raw, $idx);
              }

              if ($VERBOSE) {
                printf STDERR "upstream: %s idx=%d\n", $upstream, $idx;
                printf STDERR "upstream_full: %s\n", $upstream_full;
                printf STDERR "downstream: %s\n", $downstream;
              }

              #
              #  special checks for GeneB which do not require
              #  anchoring of GeneA:
              #
              my $map = $hit->{codon_map} || die;
              my %cm = map {$_->{ref_base_num}, $_} @{$map};
              my $cm_min = min(keys %cm);
              my $cm_max = max(keys %cm);

              my $last_cn;
              my $utr_call = SV_MATCH_CODE_FRAMESHIFT;
              # if none of below checks pass, consider out of frame

              my $max_seek = length($downstream);
              # only search the sequence present in contig downstream

              my $b_seek_dir;
              # seek direction to go downstream in the transcript mapping for GeneB
              if ($hit->{strand} eq "+") {
                $b_seek_dir = 1;
              } elsif ($hit->{strand} eq "-") {
                $b_seek_dir = -1;
              } else {
                die;
              }

              my $j = $b_final_pos;
              for (my $try = 0; $try < $max_seek; $try++) {
                my $e = $cm{$j};
                #                dump_die($e);
                if (defined $e->{codon_number}) {
                  my $cn = $e->{codon_number};
                  die "direction sanity check failed: last_CN=$last_cn this_CN=$cn" if defined($last_cn) and $cn < $last_cn;
                  # sanity check: make sure we're proceeding
                  # down the GeneB transcript
                  $last_cn = $cn;
                  if (($e->{AA} || "") eq "M" and $cn == 1) {
                    # found canonical start codon for GeneB
                    # from Yongjin's email 9/3/2014:
                    # "Search start codon M in s2, if the start codon M is the
                    #  same with geneB, then mark it as 2
                    #  (5UTR fusion, like P2RY8-CRLF2)"
                    $utr_call = SV_MATCH_CODE_UTR5_GENEB_COMPLETE_CDS;
                    last;
                  }
                }
                $j += $b_seek_dir;
                last if $j < $cm_min or $j > $cm_max;
              }

#              printf STDERR "debug: downstream=%s\n", $downstream;

              if (not($utr_call) and
                  find_qualifying_coding_frame(
                                               "-frame" => $frame_raw,
                                               "-upstream" => $upstream_full,
                                               "-downstream" => $downstream
                                              )) {
                # "if M is within geneB and there is some AAs
                # matches after M then return 3 (potential 5UTR
                # fusion, like ZEB2-CXCR4)"
                $utr_call = SV_MATCH_CODE_POTENTIAL_UTR5_FUSION;
              }
              # QUESTION: what if M occurs late in the contig sequence?
              # is M + 9 AA really bad if AA may continue past contig?

              if (1) {
                my %info;
                $info{gene_a_nm} = "";
                $info{gene_a_codon_number} = "";
                $info{gene_a_anchor_type} = "";
                # n/a
                $info{gene_b_nm} = $b_nm;
                $info{gene_b_final_pos} = $b_final_pos;
                $info{gene_b_anchor_type} = $blat_mode ? "blat" : "tuple";
                $info{gene_b_codon_number} = $b_codon_number;
                $info{inframe_call} = $utr_call;

                push @geneb_utr_matches, \%info;
              }

	      #
	      #  compare GeneB frame with GeneA transcripts:
	      #
	      my $info = {};
	      $info->{gene_b_nm} = $b_nm;
	      $info->{gene_b_aa_downstream} = $fragment;
	      $info->{gene_b_aa_downstream_extended} = $fragment_extended;
	      $info->{gene_b_final_pos} = $b_final_pos;
	      $info->{gene_b_anchor_type} = $blat_mode ? "blat" : "tuple";
	      $info->{frame_id} = $fid;
	      $info->{index} = $idx;
	      $info->{gene_b_index} = $idx;
	      $info->{frame_raw} = $frame_raw;
	      $info->{frame_index} = frame_offset($fid);

              printf STDERR "raw frame in-frame for GeneB: %s\n", $frame_raw if $VERBOSE;

              $info->{frame_upstream} = $upstream;
              # the portion of the frame in-frame for GeneB upstream
              # of the GeneB match (i.e. containing GeneA)
              $info->{gene_b_codon_number} = $b_codon_number;
              $info->{gene_b_coding_basenum} = $b_coding_basenum;


              $gene_b_inframe{$fid} = $info;
              # QC / sanity

              #
              #  compare frame aligned with GeneB's AA translation to
              #  GeneA's AA translation to determine whether in-frame
              #

              # - find GeneA transcript(s)
              my $a_nms = find_nms_from_ga(
                                           "-ga" => $ga,
                                           "-row" => $row,
                                           "-which" => "A"
                                          );
#              die join "\n", @{$a_nms};

              if (my $restrict = $FLAGS{"restrict-a-nm"}) {
                die "restrict-a-nm failed" unless grep {$_ eq $restrict} @{$a_nms};
                $a_nms = [ $restrict ];
              }

              $annot_trouble->add("geneA_no_transcript_IDs_found") unless @{$a_nms};

              foreach my $a_nm (@{$a_nms}) {

                if ($VERBOSE) {
                  printf STDERR "pair: B=%s A=%s fid=%s\n", $b_nm, $a_nm, $fid;
                  printf STDERR "geneB AA: %s\n", $r2g->get_aa("-accession" => $b_nm);
                  printf STDERR "upstream: %s\n", $upstream;
                }

                my %pair_info = %{$info};
                $pair_info{gene_a_nm} = $a_nm;

                my ($a_codon_number, $is_inframe, $a_anchor_type, $a_upstream_idx_end);

#                foreach my $try_utr5 (0, 1) {
                foreach my $try_utr5 (0) {
                  # disable 5' UTR check
                  ($a_codon_number,
                   $is_inframe,
                   $a_anchor_type,
                   $a_upstream_idx_end) = match_a(
                                   "-row" => $row,
                                   "-upstream" => $upstream,
                                   "-nm" => $a_nm,
                                   "-r2g" => $r2g,
                                   "-try-utr5" => $try_utr5,
                                   "-annot-trouble" => $annot_trouble,
                                   "-annot-global" => $annot_global,
                                  );
                  last if $is_inframe;
                }

                if (length($is_inframe) and $is_inframe == SV_MATCH_CODE_INFRAME) {
                  die "undef A idx end" unless defined $a_upstream_idx_end;

                  my $after_a = substr($upstream_full, $a_upstream_idx_end);
                  my $fail;

                  my @check = $after_a;
                  push @check, $fragment if $utr_codons;
                  # for 5' UTR types, a stop may also appear in
                  # geneB lookup start, so search the fragment too.
                  #
                  # FIX ME: what if additional stops downstream but
                  # before canonical start??

#                  foreach ($after_a, $fragment) {
                  foreach (@check) {
                    next unless $_;
                    $fail = 1 if index($_, STOP_CODON_CODE) != -1;
                  }

                  if ($fail) {
                    $is_inframe = SV_MATCH_CODE_FRAMESHIFT;
                    # not exactly
                    $annot_global->add("not_inframe_due_to_intermediate_stop_codon");
                  }
                }

                #
                #  generate exon # annotation, if possible:
                #
                add_exon_annotation($row, \%pair_info, $r2g);
                # FIX ME: may 'splode unless both genes anchorable

                $pair_info{inframe_call} = $is_inframe;
                $pair_info{gene_a_codon_number} = $a_codon_number;
                $pair_info{gene_a_anchor_type} = $a_anchor_type;
                if (defined $a_upstream_idx_end) {
                  my $ai = $a_upstream_idx_end - 1;
                  # change to be the last position of gene A,
                  # similar to gene_b_index (first position of gene B)
                  $pair_info{gene_a_index} = $ai;
                  my $bi = $pair_info{gene_b_index};
                  my $aa_between = "";
                  if ($bi > $ai) {
                    $aa_between = substr($frame_raw, $ai + 1, $bi - $ai - 1);
                  }
                  $pair_info{interstitial_AA} = $aa_between;

                  my $fo = frame_offset($fid);

                  my $contig_index_a = ($ai * 3) + $fo + 2;
                  # contig index of last base of last codon of geneA anchor

                  my $contig_index_b = ($bi * 3) + $fo;
                  # contig index of 1st base of 1st codon of geneB anchor

                  $pair_info{gene_a_contig_index} = $contig_index_a;
                  $pair_info{gene_b_contig_index} = $contig_index_b;
                }

                if ($pair_info{gene_a_anchor_type} and
                    $pair_info{gene_b_anchor_type} and
                    not($is_inframe)) {
#                  dump_die(\%pair_info, "WTF: 2 anchors but no inframe call?");
                  # => stop codon in AA sequence, see above
                }

                foreach my $f (qw(inframe_call gene_a_codon_number gene_a_anchor_type)) {
                  die "$f not defined" unless defined $pair_info{$f};
                }

                push @pairs, \%pair_info;
              }
            } else { # $idx != -1
              $b_global_inframe = 0 unless length($b_global_inframe);
              # evaluated but not found: init result code if necessary
            }
          }  # $fid
        } else {
          # refSeq IDs identified for GeneB, but no fully-mapped CDS matches
          # - gene model mapping problem?
          $processing_exception = 1;
          $annot_trouble->add("geneB_breakpoint_not_mapped_to_CDS", $b_nm);
        }

        unless ($processing_exception) {
          #
          # post-processing QC for GeneB:
          #
          my $hit_frame_count = scalar keys %gene_b_inframe;
          if ($hit_frame_count > 1) {
            $annot_trouble->add("ERROR_geneB_multiple_matching_frames", $b_nm);
          } elsif ($hit_frame_count == 0) {
            #
            # unable to identify frame compatible with GeneB excerpt
            #

            if (%frame_blat) {
              # WTF: successful BLAT of GeneB AA vs. frames
#              die "no frame alignment for GeneA, but blat hit, CHECK ME";
              $annot_trouble->add("geneB_anchor_fail_but_found_blat_hit", $b_nm);
              # possibly subtle mismatch between
              #   (a) geneA AA derived from reference sequence
              #   (b) geneA AA in refSeq record
              #   (c) geneA AA segment derived from Cicero contig sequence
              # and/or breakpoint site problems
              # see no_frames_but_blat_hit.tab + slides
              #
              # => hopefully this code now no longer reachable since
              #    BLAT results are used in frame_search()
            } else {
              # no BLAT hit
              # not sure if we ever want to lower minScore,
              # e.g. no_geneA_frame_match_blat_minscore19.tab;
              # in this case a hit of 10 with 1 mismatch, nothing further
              printf STDERR "GeneB: no frames found and no BLAT hit\n";
              $annot_trouble->add("geneB_no_inframe_contig_found");
            }
          }
        }
      }                                # $b_nm

      my $sv_refseqA_inframe;
      my %sv_refseqA_inframe_transcripts;
      if ($junction_mode and $df->headers->{junction_peptide}) {
        #
        #  in junction mode with user-specified peptide,
        #  separate reporting of GeneA in-frame status.
        #
        foreach my $p (@pairs) {
          if ($p->{inframe_call} == SV_MATCH_CODE_INFRAME) {
            $sv_refseqA_inframe = 1;
            # derivable from existing results
            $sv_refseqA_inframe_transcripts{$p->{gene_a_nm}} = 1;
          }
        }

        unless ($sv_refseqA_inframe) {
          # evaluation needed
          my @hits;
          foreach my $fid (keys %{$frames}) {
            my $frame = $frames->{$fid};
            push @hits, $frame unless index($frame, $junction_peptide) == -1;
          }
          if (@hits) {
            # might not be any matches if user peptide doesn't match
            # any of the 3 frames (rare)
            die "expected 1 hit, got " . scalar @hits unless @hits == 1;

            my $frame = $hits[0];
            # frame in-frame with user-specified peptide

            $sv_refseqA_inframe = genea_frame_check(
                                                    "-frame" => $hits[0],
                                                    "-row" => $row,
                                                    "-ga" => $ga,
                                                    "-r2g" => $r2g,
                                                    "-annot-global" => $annot_global,
                                                    "-annot-trouble" => $annot_trouble,
                                                    "-results" => \%sv_refseqA_inframe_transcripts
                                                   );
          }

        }
      }

#      dump_die($row, "> 1 pairs evaluated") if @pairs > 1;
      # looking for example

      #
      #  report:
      #
      my $sv_any_inframe = "";
      my $sv_any_inframe_canonical = "";
      my $sv_inframe = "";
      my $sv_aa = "";
      my $sv_refseqA = "";
      my $sv_refseqB = "";
      my $sv_refseqA_codon = "";
      my $sv_refseqB_codon = "";
      my $sv_refseqA_exon = "";
      my $sv_refseqB_exon = "";
      my $sv_refseqB_coding_base_number = "";
      my $sv_refseqA_anchor = "";
      my $sv_refseqB_anchor = "";
      my $sv_posB_adjusted = "";
      my $sv_inframe_desc = "";
      my $sv_refseqA_AA_index = "";
      my $sv_refseqB_AA_index = "";
      my $sv_interstitial_AA = "";
      my $sv_refseqA_contig_index = "";
      my $sv_refseqB_contig_index = "";
      my $sv_frame_index = "";

      my @eval = @pairs;
      my $report_utr = 1;
      foreach my $pair (@eval) {
        $report_utr = 0 if $pair->{inframe_call};
        # if paired inframe eval
      }
      push @eval, @geneb_utr_matches if $report_utr;
      # if any paired comparisons are in-frame, no need to
      # report checks for geneB alone

      if (@eval) {
        #
        #  sometimes geneB is anchored to more than one putative frame.
        #  If one of these is in-frame with geneA and one isn't,
        #  just keep the in-frame result.
        #
        my %saw;
        foreach my $pair (@eval) {
          my $nm_a = $pair->{gene_a_nm} || next;
          my $nm_b = $pair->{gene_b_nm} || next;
          my $inframe = $pair->{inframe_call};
          next unless length($inframe);
          $saw{$nm_a}{$nm_b}{$inframe} = 1;
        }

        my %prune;
        foreach my $nm_a (keys %saw) {
          foreach my $nm_b (keys %{$saw{$nm_a}}) {
            my $types = $saw{$nm_a}{$nm_b};
            $prune{$nm_a}{$nm_b} = 1 if $types->{SV_MATCH_CODE_FRAMESHIFT()} and $types->{SV_MATCH_CODE_INFRAME()};
          }
        }

        if (%prune) {
          my @filtered;
          foreach my $pair (@eval) {
            my $usable = 1;
            my $nm_a = $pair->{gene_a_nm} || next;
            my $nm_b = $pair->{gene_b_nm} || next;
            my $inframe = $pair->{inframe_call};
            if ($prune{$nm_a}{$nm_b} and $inframe == SV_MATCH_CODE_FRAMESHIFT) {
              $usable = 0;
              printf STDERR "prune %s %s\n", $nm_a, $nm_b;
            }
            push @filtered, $pair if $usable;
          }
          @eval = @filtered;
        }

      }

      if (@eval) {
        #
        # some evaluations were made
        #
        my %inframe_symbols;

        foreach my $pair (@eval) {
          my $call = $pair->{inframe_call};
          if (length($call)) {
            # actual call (non-blank)
            unless (length($sv_any_inframe)) {
              # initialize output only if at least one call is available
              $sv_any_inframe = 0;
              $sv_any_inframe_canonical = 0;
            }
            if ($call) {
              $sv_any_inframe = 1;

              my @inframe_ends;
              if ($call eq SV_MATCH_CODE_INFRAME) {
                # both ends in-frame
                $sv_any_inframe_canonical = 1;
                @inframe_ends = qw(a b);
              } elsif ($call eq SV_MATCH_CODE_UTR5_GENEB_COMPLETE_CDS) {
                # geneB only
                @inframe_ends = qw(b);
              }

              foreach my $end (@inframe_ends) {
                my $nm = $pair->{"gene_" . $end . "_nm"} || die;
                my $gene = $ga->find_accession($nm, "-symbol" => 1) || die;
                $inframe_symbols{$end}{$gene} = 1;
              }

            }
          }
        }

        if ($REANNOTATE_INFRAME_GENES) {
          foreach my $end (qw(a b)) {
            if ($inframe_symbols{$end}) {
              my @syms = keys %{$inframe_symbols{$end}};
              if (@syms == 1) {
                # exactly one symbol in in-frame on this end
                my $key = sprintf 'gene%s', uc($end);
                my $current = $row->{$key};
                if ($current and $current ne $syms[0]) {
                  $annot_global->add(sprintf 'reannotated %s from %s to %s',
                                     $key, $current, $syms[0]);
                  $row->{$key} = $syms[0];
                }
              }
            }
          }
        }

        $sv_aa = join ",", map {$_->{frame_raw} || ""} @eval;
        $sv_refseqA = join ",", map {$_->{gene_a_nm}} @eval;
        $sv_refseqB = join ",", map {$_->{gene_b_nm}} @eval;
        $sv_refseqA_codon = join ",", map {$_->{gene_a_codon_number}} @eval;
        $sv_refseqB_codon = join ",", map {$_->{gene_b_codon_number}} @eval;
        $sv_refseqA_exon = join ",", map {$_->{gene_a_exon_number} || ""} @eval;
        $sv_refseqB_exon = join ",", map {$_->{gene_b_exon_number} || ""} @eval;
#        $sv_refseqB_coding_base_number = join ",", map {$_->{gene_b_coding_basenum} || ""} @eval;
        $sv_refseqA_anchor = join ",", map {$_->{gene_a_anchor_type}} @eval;
        $sv_refseqB_anchor = join ",", map {$_->{gene_b_anchor_type}} @eval;
        $sv_inframe = join ",", map {$_->{inframe_call}} @eval;
        $sv_posB_adjusted = join ",", map {$_->{gene_b_final_pos}} @eval;
#        $sv_refseqA_AA_index = join ",", map {$_->{gene_a_index}} @eval;
        $sv_refseqA_AA_index = build_cdl(\@eval, "gene_a_index");
        $sv_refseqB_AA_index = build_cdl(\@eval, "gene_b_index");
        $sv_interstitial_AA = build_cdl(\@eval, "interstitial_AA");
        $sv_refseqA_contig_index = build_cdl(\@eval, "gene_a_contig_index");
        $sv_refseqB_contig_index = build_cdl(\@eval, "gene_b_contig_index");
        $sv_inframe_desc = join ";", map {$_->{inframe_desc} || ""} @eval;
        $sv_frame_index = build_cdl(\@eval, "frame_index");
      } else {
        #
        # no paired evaluations made.
        #
        # 12/2015: report accessions even in these cases
        #
        my $nms_a = find_nms_from_ga(
                                     "-ga" => $ga,
                                     "-row" => $row,
                                     "-which" => "A"
                                     );
        my $nms_b = find_nms_from_ga(
                                     "-ga" => $ga,
                                     "-row" => $row,
                                     "-which" => "B"
                                     );

        my $any_intergenic_salvaged;

        $any_intergenic_salvaged = 1 if salvage_intergenic_genes(
                                 "-row" => $row,
                                 "-which" => "A",
                                 "-nms" => $nms_a,
                                 "-rff" => $r2g->rff(),
                                );
        $any_intergenic_salvaged = 1 if salvage_intergenic_genes(
                                 "-row" => $row,
                                 "-which" => "B",
                                 "-nms" => $nms_b,
                                 "-rff" => $r2g->rff(),
                                );
        # for intergenic events, find nearby refseqs

        if ((@{$nms_a} and !@{$nms_b}) or
            (@{$nms_b} and !@{$nms_a})) {
          # we have a refSeq on one end but not the other (e.g. IGH)
          # create a dummy entry for the empty set so we will at
          # least report the refSeq for the other gene.
          push @{$nms_a}, "" unless @{$nms_a};
          push @{$nms_b}, "" unless @{$nms_b};
        }

        my @eval;
        # generate a pairing for each A/B refGene as if processing
        # proceeded normally

        foreach my $nm_a (@{$nms_a}) {
          foreach my $nm_b (@{$nms_b}) {
            my %r;
            $r{gene_a_nm} = $nm_a;
            $r{gene_b_nm} = $nm_b;
            push @eval, \%r;
          }
        }
         $sv_refseqA = join ",", map {$_->{gene_a_nm}} @eval;
         $sv_refseqB = join ",", map {$_->{gene_b_nm}} @eval;

        my @blank = map {""} @eval;

        foreach my $v (
                       \$sv_inframe,
                       \$sv_aa,
                       \$sv_refseqA_codon,
                       \$sv_refseqB_codon,
                       \$sv_refseqA_exon,
                       \$sv_refseqB_exon,
                       \$sv_refseqB_coding_base_number,
                       \$sv_refseqA_anchor,
                       \$sv_refseqB_anchor,
                       \$sv_posB_adjusted,
                       \$sv_refseqA_AA_index,
                       \$sv_refseqB_AA_index,
                       \$sv_interstitial_AA,
                       \$sv_refseqA_contig_index,
                       \$sv_refseqB_contig_index,
                       \$sv_frame_index,
                      ) {
          $$v = join ",", @blank;
        }
        $sv_inframe_desc = join ";", @blank;

        if ($any_intergenic_salvaged) {
          # call event not in frame for intergenic events.
          # HOWEVER: what if before 5' UTR?
          $sv_inframe = join ",", map {0} @eval;
        }
      }

      $row->{sv_general_info} = $annot_global->get_field();
      $row->{sv_processing_exception} = $annot_trouble->get_field();

      $row->{sv_AA} = $sv_aa;
      $row->{sv_any_inframe} = $sv_any_inframe;
      $row->{sv_any_inframe_canonical} = $sv_any_inframe_canonical;
      $row->{sv_inframe} = $sv_inframe;
      $row->{sv_refseqA} = $sv_refseqA;
      $row->{sv_refseqB} = $sv_refseqB;
      $row->{sv_refseqA_codon} = $sv_refseqA_codon;
      $row->{sv_refseqB_codon} = $sv_refseqB_codon;
      $row->{sv_refseqA_exon} = $sv_refseqA_exon;
      $row->{sv_refseqB_exon} = $sv_refseqB_exon;
      $row->{sv_refseqB_coding_base_number} = $sv_refseqB_coding_base_number;
      $row->{sv_refseqA_anchor_type} = $sv_refseqA_anchor;
      $row->{sv_refseqB_anchor_type} = $sv_refseqB_anchor;
      $row->{sv_posB_adjusted} = $sv_posB_adjusted;
#      $row->{sv_inframe_desc} = $sv_inframe_desc;
      $row->{sv_desc} = $sv_inframe_desc;
      # we now attempt to always report

      $row->{sv_refseqA_inframe} = $sv_refseqA_inframe;
      $row->{sv_refseqA_inframe_transcripts} = join ",", sort keys %sv_refseqA_inframe_transcripts;
      $row->{sv_refseqB_inframe} = $b_global_inframe;
      $row->{sv_refseqB_inframe_transcripts} = join ",", sort keys %b_global_inframe_transcripts;
      $row->{sv_refseqA_AA_index} = $sv_refseqA_AA_index;
      $row->{sv_refseqB_AA_index} = $sv_refseqB_AA_index;
      $row->{sv_interstitial_AA} = $sv_interstitial_AA;
      $row->{sv_frame_index} = $sv_frame_index;
      $row->{sv_refseqA_contig_index} = $sv_refseqA_contig_index;
      $row->{sv_refseqB_contig_index} = $sv_refseqB_contig_index;
      foreach my $end (qw(A B)) {
        my $key = sprintf 'sv_refseq%s_annot', $end;
        $row->{$key} = gene_feature_annot(
                                          "-row" => $row,
                                          "-which" => $end,
                                          "-ga" => $ga
                                         );
      }

      patch_row_coding_bases("-row" => $row,
                             "-r2g" => $r2g);
#      $rpt->end_row($row);

      if ($rpt_full_length) {
	report_full_length(
			   "-row" => $row,
			   "-eval" => \@eval,
			   "-nm2protein" => \%nm2protein,
			   "-rpt" => $rpt_full_length,
			   "-r2g" => $r2g,
			  );
	# may have to be updated to deal with pre-sorting if necessary
      }
    }

    foreach my $r (@rows_raw) {
      # write rows in ORIGINAL order regardless of whether they were
      # processed in a sorted order
      $rpt->end_row($r);
    }
    $rpt->finish();

    if ($rpt_full_length) {
      $rpt_full_length->finish();
      if ($FLAGS{"run-digest"}) {
	$FLAGS{"digest-full-length"} = $f_out_full_length;
	digest_full_length();
      }
    }

    my $end_time = time;
    printf STDERR "processing time: %d\n", $end_time - $start_time if $VERBOSE;

  }

}

sub get_frames {
  my ($contig) = @_;

  my $ct = new Bio::Tools::CodonTable();

  my @nt = split //, $contig;

  my %frames;
  foreach my $offset (0 .. 2) {
    my $string = substr($contig, $offset);
    my @codons = $string =~ /(...)/g;
    my $id = sprintf "frame_offset_%d", $offset;
    my $sequence = join "", map {$ct->translate($_)} @codons;
    printf STDERR "frame %s: %s\n", $id, $sequence if $VERBOSE;
    $frames{$id} = $sequence;
  }
  return \%frames;
}

sub find_nms_from_ga {
  # find accession from position via GeneAnnotation.pm
  my (%options) = @_;
  my $ga = $options{"-ga"} || die "-ga";
  my $row = $options{"-row"};
  my $which = $options{"-which"};
  my $strand_mode = $options{"-strand"};
  my $rows_mode = $options{"-return-rows"};
  my $junction_strand = $row->{junction_strand};

  my $ref_name;
  unless ($ref_name = $options{"-chrom"}) {
    die unless $row and $which;
    $ref_name = get_chrom($row, $which) || die;
  }

  my $pos;
  unless ($pos = $options{"-position"}) {
    $pos = $row->{"pos" . $which} || $row->{"Pos" . $which} || dump_die($row, "can't find pos A");
    #  note that geneB breakpoint may NOT be at the same breakpoint
    #  indicated by the contig sequence, see example.tab.
    #  mapping ambiguity??
    #  - breakpoint appears to be at the soft clip site
  }

  confess "need ref_name" unless $ref_name;
  confess "need pos" unless $pos;

  my %nm;

#  foreach my $try (0, 1) {
  foreach my $try (0) {
    # expanded search disabled until we find an example where it's needed
    my $seek_start = $pos;
    my $seek_end = $pos;
    if ($try > 0) {
      # first attempt might fail if site is not in CDS.
      # Expand search in case fusion event lands in e.g. 5' UTR.
      $seek_start -= $GENE_FLANK_SEARCH_NT;
      $seek_end += $GENE_FLANK_SEARCH_NT;
    }

    if ($ga->find(
                  "-reference" => $ref_name,
                  "-start" => $seek_start,
                  "-end" => $seek_end
                 )) {
      my $hits = $ga->results_rows();

      if ($junction_strand) {
        my @filtered;
        foreach my $h (@{$hits}) {
          push @filtered, $h if $h->{strand} eq $junction_strand;
          # if user specified junction strand orientation,
          # only process isoforms w/that strand
        }
        $hits = \@filtered;
      }

      if ($strand_mode) {
        %nm = map {($_->{strand} || die "no strand"), 1} @{$hits};
      } else {
        %nm = map {($_->{name} || die "no acc"), $_} @{$hits};
      }
      die "TEST ME: found gene after expanding search" if $try;
      last;
    } else {
#      dump_die($row, "can't find gene hit for site try=$try $ref_name $seek_start $seek_end", 1);
    }
  }

  if ($rows_mode) {
    return [values %nm];
  } else {
    return [sort keys %nm];
  }
}

sub get_blat {
  my $blat = new BLATer();
  $blat->null_mode(1) if $FLAGS{"no-blat"};
  $blat->minScore($BLAT_MIN_SCORE);
  if ($FLAGS{"save-tempfiles"}) {
    printf STDERR "saving tempfiles\n";
    $blat->tfw->auto_unlink(0);
  }
  return $blat;
}

sub test_single {
  #
  # test single sequence:
  #
  my ($acc) = @_;
  my $r2g = new RefFlat2Genome(
                               "-fasta" => get_fasta_file()
                              );
  $r2g->parse(
              "-refflat" => get_refflat_file(),
              "-refgene2protein" => get_r2p_file()
             );

  $r2g->verbose(1);


#    my $acc = "NM_005933";
    # mapped to +, OK

#    my $acc = "NM_004529";
    # mapped to -, OK

#    my $acc = "NM_000014";
    # single mismatch in AA translation: ???
    # - doesn't seem to involve a split exon

#    my $acc = "NM_000658";
    # no AA translation available from parsed refFlat dump
    # => suppressed sequence

    # my $acc = "NM_001022";
    # protein length translation error (!)
    # => convert cds_start to in-base before exon start check
    # => check to make sure start exon index is valid

#    my $acc = "NM_001075";
    # protein length translation error #2 (!)
    # => don't understand this translation at all, mine looks correct

#    my $acc = "NM_000363";
    # ERROR: AA translation length mismatch for NM_000363! expected:210 got:4 code=-3

#  $acc = "NM_001987";
  # ETV6, +

  $r2g->find("-accession" => $acc);
}

sub test_all {
  #
  #  verify all translations:
  #
  my $r2g = new RefFlat2Genome(
                               "-fasta" => get_fasta_file()
                              );
  $r2g->parse(
              "-refflat" => ($FLAGS{refgene} || die),
              "-refgene2protein" => get_r2p_file()
             );

  my $set = $r2g->by_accession();
  printf "total:%d\n", scalar keys %{$set};
  foreach my $id (sort keys %{$set}) {
    printf STDERR "processing %s...\n", $id;
    $r2g->find("-accession" => $id);
  }
}

sub get_aa_fragment {
  my %options = @_;
  my $hit = $options{"-hit"} || die;
  my $index = $options{"-index"};
  die unless defined $index;
  my $length = $options{"-length"} || die;
  my $direction = $options{"-direction"} || die;

  my $map = $hit->{codon_map} || die;
  my %saw_codons;
  my %utr_codons;
  my @aa;

  my $i = $index;
  my $utr_wiggle_enabled = $UTR_TUPLE_SEARCH_WIGGLE_CODONS;
  my $all_utr3 = 1;

  my $SKIP_EDGE_BASE = 0;
  # Now disabled.  The original reason for this was to prevent altered
  # codons at SV boundaries from breaking tuple lookup check.  OTOH:
  # - blat is now used as a backup, so this is no longer fatal
  # - this behavior guarantees a gap even if the codon boundaries
  #   are PERFECT on both edges of the SV.

  my $build_length = $length;
  $build_length++ if $SKIP_EDGE_BASE;

  while (1) {
#    last if $map->[$i]->{in_utr};
    # ran out of sequence!
    # we may need to go into the UTR, see e.g. geneA_no_upstream_fragment.tab
    # where the SV occurs in 1st codon!

    my $cnum = $map->[$i]->{codon_number};
    last unless defined $cnum;
    # straggler nucleotides at extreme 5' / 3' will not
    # be mapped to a codon # (insufficient)

    if (my $utr_type = $map->[$i]->{in_utr}) {
#      dump_die($map->[$i]);
      $utr_codons{$cnum} = 1;
      if ($utr_type == 3) {
      } elsif ($utr_type == 5) {
        $all_utr3 = 0;
      } else {
        die "unknown UTR code";
      }
    } else {
      $utr_wiggle_enabled = 0;
      $all_utr3 = 0;
    }

#    printf STDERR "at %d, cnum=%d\n", $i, $cnum;
    unless ($saw_codons{$cnum}) {
      my $aa = $map->[$i]->{AA};
      unless ($aa) {
        if ($i <= 1 or (@{$map} - $i <= 2)) {
          # if at extreme end of 5' or 3' UTR, it's possible to
          # have orphaned bases which are insufficient to form a
          # fake codon.  Not an error, just no data
          last;
        } else {
          # shouldn't happen
          dump_die($map->[$i], "no AA entry and not at UTR end! " . scalar @{$map} . " " . $i);
        }
      }

      if ($utr_wiggle_enabled) {
        # a leading UTR codon which may be troublesome if posA/B mapping
        # is not accurate
        $utr_wiggle_enabled--;
        # toss codon
      } else {
        if ($direction > 0) {
          push @aa, $aa;
        } else {
          unshift @aa, $aa;
        }
        last if @aa == $build_length;
      }

      $saw_codons{$cnum} = 1;
    }

    $i += $direction;
    last if $i < 0 or $i >= @{$map};
  }

  printf STDERR "WARNING: too-short AA %s\n", join("", @aa) unless @aa >= $length;
  # might happen if:
  # - breakpoint too close to start of transcript
  # - non-unique AA, extension required and hit start of transcript

  if ($SKIP_EDGE_BASE) {
    #
    # exclude the breakpoint codon as this might change in CICERO contig seq:
    # see e.g. JZ_SJRB051_D_example2_out_of_frame.tab
    # searching for fragment SWITFLLAKG in frame DC*RSMMYCLHSSANWKGHVNLYI*HNPAVRYLLK*ILHWC*KFLGSHFY*LKWMLTKIVQENL*TKTKETNGEKYLKYF
    # searching for fragment SWITFLLAKG in frame LLKKYDVLFALFSKLERTCELIYLTQPSSSISTEINSALVLKVSWITFLLAKVDADENCTREFMNKNKRDKWRKIFEVL
    #   => FAIL!
    # searching for fragment SWITFLLAKG in frame TVEEV*CIVCTLQQIGKDM*TYIFDTTQQFDIY*NKFCIGAKSFLDHIFIS*SGC*RKLYKRIYEQKQKRQMAKNI*STF

    # codon for this site spans exons:
    # 48923156,48923157,48923158 => AAA => K (202)
    # >>>48923159<<<,48934153,48934154 => GGG => G (203)
    # 48934155,48934156,48934157 => GAA => E (204)
    #
    #   => the codon containing the breakpoint should not be considered during
    #      the string lookup phase because Cicero's contig sequence might
    #      not resolve to the same codon!
    if ($direction == -1) {
      pop @aa;
    } elsif ($direction == 1) {
      shift @aa;
    } else {
      die;
    }
  }

  my $fragment = join "", @aa;

  return ($fragment, scalar(keys %utr_codons), $all_utr3);
}

sub get_feature {
  my ($row, $which) = @_;
  my @keys;
  push @keys, "feature" . $which;
  push @keys, "annotate" . $which;
  my $feature;
  foreach (@keys) {
    $feature = $row->{$_};
    last if $feature;
  }
  return $feature;
}

sub get_chrom {
  my ($row, $which) = @_;

  my @keys = "chr" . $which;
  push @keys, "Chr" . $which;
  push @keys, "#chr" . $which if $which eq "A";
  # CREST

  my $ref_name;
  foreach (@keys) {
    $ref_name = $row->{$_};
    last if $ref_name;
  }
  dump_die($row, "can't find chr") unless $ref_name;
  return $ref_name;
}

sub get_base_number {
  my ($row, $which) = @_;
  my $pos = $row->{"pos" . $which} || $row->{"Pos" . $which} || dump_die($row, "can't find pos A");
  return $pos;
}

sub frame_search {
  my (%options) = @_;
  my $frame = $options{"-frame"} || die;
  my $hit = $options{"-hit"} || die;
  my $codon_index = $options{"-index"};
  die unless defined $codon_index;
  my $length = $options{"-length"} || die;
  my $direction = $options{"-direction"} || die;
  my $end = $options{"-end"} || die;
  my $annot_trouble = $options{"-annot-trouble"} || new AnnotationField("-unique" => 1);
  my $annot_global = $options{"-annot-global"} || new AnnotationField("-unique" => 1);

#  dump_die($hit, $end);

  my ($idx, $fragment, $frame_error, $utr_codons, $all_utr3);
  my ($this_frag_length, $last_frag_length, $utr3_fail);

  while (1) {
    ($fragment, $utr_codons, $all_utr3) = get_aa_fragment(
                                               "-hit" => $hit,
                                               "-index" => $codon_index,
                                               "-length" => $length,
                                               "-direction" => $direction
                                              );

    $utr3_fail = 1 if $all_utr3 and $NOT_INFRAME_IF_ENTIRELY_UTR3;

    if ($utr3_fail) {
#      $frame_error = sprintf '%sstream_AA_entirely_in_utr3', $direction == 1 ? "down" : "up";
      # not really an error so don't report one, since this will
      # result in an inframe call of null rather than 0.
      # A value of 0 is appropriate since we know we don't want
      # inframe calls to these regions.
      $idx = -1;
      last;
    }

    unless ($fragment) {
      # geneA_no_upstream_fragment.tab
      $frame_error = sprintf "no_%sstream_AA", $direction == 1 ? "down" : "up";
      $idx = -1;
      last;
    }

    if (0) {
      printf STDERR "DEBUG: disabling frag lookup!!!\n";
      $fragment = "xyzzy1234";
    }

    printf STDERR "searching for fragment %s in frame %s\n", $fragment, $frame if $VERBOSE;
    my $this_frag_length = length($fragment);
    if (defined($last_frag_length) and $this_frag_length <= $last_frag_length) {
      $frame_error = "short_upstream_AA_and_cannot_expand";
      $idx = -1;
      last;
    }
    $last_frag_length = $this_frag_length;

    $idx = index($frame, $fragment);
    if ($idx == -1) {
      # not found, quit
      last;
    } else {
      # verify sequence is unique in frame.
      # fragment must be unique to properly anchor match position in contig.
      # if not, expand length of search AA fragment
      my $idx2 = index($frame, $fragment, $idx + 1);
      if ($idx2 == -1) {
        # fragment is unique, stop
        last;
      } else {
        printf STDERR "WARNING: multiple hits in AA for fragment %s in %s; expanding fragment\n", $fragment, $frame;
        $length++;
      }
    }
  }

  my $blat_mode = 0;
#    if ($all_utr3 and $NOT_INFRAME_IF_ENTIRELY_UTR3) {

  my $blat = $options{"-blat"};
  if ($idx == -1 and $blat and not($utr3_fail)) {
    # try blat if not found, unless we're in 3' UTR and don't
    # want inframe calls there
    die "multiple hits!" if @{$blat} > 1;
    my $bhit = $blat->[0];

    my @hsp;
    my %hsp_problems;
    while (my $hsp = $bhit->next_hsp) {
      if (hsp_filter($hsp, \%hsp_problems)) {
        printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
          $hsp->score,
            $hsp->strand,
              $hsp->range("query"),
                $hsp->range("hit"),
                  $hsp->num_identical(),
                    $hsp->frac_identical("query"),
                      $hsp->length("query"),
                        $hsp->length("hit"),
                          $hsp->length("total"),
                            $hsp->query_string(),
                              $hsp->hit_string();

        push @hsp, $hsp;
      }
    }

    if (@hsp > 1) {
#      $hsp_problems{multiple_filtered_hsp} = 1;
      my @sorted = sort {$b->score <=> $a->score} @hsp;
      @hsp = $sorted[0];
      # if multiple hits, use highest-scoring one
      $annot_global->add(sprintf "gene%s_multiple_blat_hits=%s", $end, $hit->{name});
    }

    if (@hsp == 1) {
      my $hsp = $hsp[0];
      $blat_mode = 1;
      $idx = ($hsp->range("hit"))[0] - 1;
      $fragment = substr($frame, $idx);
      if (length($fragment) > $MIN_TUPLE_AA) {
        # trim the fragment to size expected by ordinary tuple check.
        # Leaving it full length may interfere with other checks, e.g.
        # for a stop in 5' UTR.
        $fragment = substr($fragment, 0, $MIN_TUPLE_AA);
      }
      # TO DO: ANCHOR IN UTR5 ALSO?
    } elsif (%hsp_problems) {
      $annot_trouble->add(sprintf("gene%s_blat_issues", $end),
                          [ sort keys %hsp_problems ]);
    }
  }

  if ($ENABLE_SYNTHETIC_AA_INTRON and
      $idx == -1 and
      not($utr3_fail)
     ) {
    # TO DO:
    # - determine if site is actually in intron
    # - generate fake protein sequence
    # - index and blat lookup (CAN'T use blat info above as this
    #   was for refGene only)
    ($idx, $fragment) = intronic_frame_search(%options);
  }

  return ($idx, $fragment, $frame_error, $utr_codons, $blat_mode);
}


sub get_outfile {
  my ($infile) = @_;
  my $outfile;
  if ($FLAGS{"sj-tag-outfile"}) {
    if ($infile =~ /(SJ\w+[0-9]+\w+)/) {
      $outfile = sprintf "%s_%s.frame.tab", $1, basename($infile);
    } else {
      printf STDERR "WARNING: can't find SJ tag in %s\n", $infile;
      $outfile = basename($infile) . ".frame.tab";
    }
  } elsif ($FLAGS{fq}) {
    $outfile = $infile . ".frame.tab";
  } else {
    $outfile = basename($infile) . ".frame.tab";
  }
  return $outfile;
}

sub extract_example {
  # extract last site worked on from error log
  my $log = $FLAGS{"extract-example"};
  open(IN, $log) || die;
  my $fn;
  my ($pos_a, $pos_b);
  while (<IN>) {
    if (/processing file (.*)\.\.\./) {
      $fn = $1;
    } elsif (/^\s+posA: (\d+)/) {
      $pos_a = $1;
    } elsif (/^\s+posB: (\d+)/) {
      $pos_b = $1;
    }
  }

  die unless $fn and $pos_a and $pos_b;

  my $outfile = "example.tab";

  my $cmd = sprintf "egrep -i 'posa|%d.*%d' %s > %s\n", $pos_a, $pos_b, $fn, $outfile;
  die $cmd;



  die $log;
  die;
}

sub get_contig {
  my ($row) = @_;
  my $contig = $row->{contig} || $row->{ConsensusSeq} || die "row does not contain contig or ConsensusSeq field";
  return $contig;
}

sub is_cicero {
  # is this row from CICERO?
  my ($row) = @_;
  return exists $row->{sv_ort};
}

sub extract_inframe_hack {
  # pull out SVs Yongjin has identified as in-frame for testing
  my $lines = read_simple_file("inframe_fusions.txt");
  my %track;
  foreach my $line (@{$lines}) {
    my @f = split /\t/, $line;
    die unless @f == 3;
    my ($sample, $g1, $g2) = @f;
    my $key = join "_", $sample, $g1, $g2;
    $track{$key} = 1;
  }

  my $infile = "PHL5.final_fusions.txt";
  my $df = new DelimitedFile("-file" => $infile,
                             "-headers" => 1,
                             );
  my $rpt = $df->get_reporter(
                              "-file" => $infile . ".inframe.tab",
                             );

  my @delete;

  my %matches;

  while (my $row = $df->get_hash()) {
#    dump_die($row);
    my $usable;
    my $key = join "_", @{$row}{qw(sample geneA geneB)};
    my $key2 = join "_", @{$row}{qw(sample geneB geneA)};

    if ($track{$key}) {
      $rpt->end_row($row);
      push @{$matches{$key}}, $row;
      push @delete, $key;
    } elsif (0 and $track{$key2}) {
      # appropriate??
      $rpt->end_row($row);
      delete $track{$key2};
    }
  }
  $rpt->finish();

  foreach my $key (keys %matches) {
    my $set = $matches{$key};
    my %fa = map {$_->{featureA}, 1} @{$set};
    if (keys %fa == 1) {
      printf STDERR "single type for %s: %s\n", $key, join ",", keys %fa;
    }
  }

  foreach my $key (@delete) {
    # only delete at finish, may be multiple entries
    delete $track{$key};
  }

  printf STDERR "NOT FOUND:\n";
  foreach (sort keys %track) {
    print STDERR "  $_\n";
  }


}

sub get_orientation {
  my ($row, $type) = @_;
  return $row->{"ort" . $type};
}

sub is_utr5 {
  my ($feature) = @_;
  die "unknown feature code $feature" unless $KNOWN_FEATURES{$feature};
  return $feature eq "5utr";
}

sub find_qualifying_coding_frame {
  my (%options) = @_;
  my $upstream = $options{"-upstream"};
  confess "-upstream" unless defined $upstream;
  # might be blank in some (wack?) cases
  my $downstream = $options{"-downstream"} || die "-downstream";
  my $frame;
  if ($options{"-no-sanity"}) {
    $frame = $upstream . $downstream;
  } else {
    $frame = $options{"-frame"} || die;
    die "mismatch" unless $upstream . $downstream eq $frame;
  }
  my $downstream_idx = index($frame, $downstream);
  die if $downstream_idx == -1;
#  die $downstream_idx;

  my $pattern = sprintf "(M[^\\*]{%d,})", $MIN_AA_DOWNSTREAM_OF_M_IN_UTR_CHECK;
  my $ok = 0;

  while ($frame =~ /$pattern/g) {
    my $pos = pos($frame);
    my $end_idx = $pos - 1;
    # last base index
    my $site_idx = $pos - length($1);
    printf STDERR "frame:%s found %s at %d, end=%d ds_idx=%d\n", $frame, $1, $site_idx, $end_idx, $downstream_idx if $VERBOSE;
#    $ok = 1 if $end_idx >= $downstream_idx;
    # only pass if the end of the match is in GeneB
    $ok = 1 if $end_idx >= ($downstream_idx - 1);
    # -1:
    #   - the stop codon may or may not be present in the sequence
    #     (i.e. AA may continue beyond end of contig)
    #   - if present, its coordinate is not reflected in the match position
    #     above, so adjust with -1.  This is for an edge case where
    #     the first codon in geneB is a stop for an event starting in GeneA.

  }
  return $ok;
}

sub utr_frame_tests {
  my @tests;

  push @tests, {
                "-desc" => "stop is 1st codon in GeneB",
                "-upstream" => "XXXXXXXXXXM1234567890",
                "-downstream" => "*YYYYYYYYYYYYYYYYYYYY",
                "-expected-usable" => 1
               };

  push @tests, {
                "-desc" => "downstream at border",
                "-upstream" => "YYYYYYYYYY",
                "-downstream" => "M1234567890*",
                "-expected-usable" => 1
               };

  push @tests, {
                "-desc" => "split, ends in B",
                "-upstream" => "YYYYYYM123456789",
                "-downstream" => "0*YYYYYYYYYY",
                "-expected-usable" => 1
               };

  push @tests, {
                "-desc" => "ends in B, no stop",
                "-upstream" => "XXXXXXXXXX",
                "-downstream" => "YYYYYYYYYYM12345678901",
                "-expected-usable" => 1
               };

  push @tests, {
                "-desc" => "size OK, but entirely in geneA",
                "-upstream" => "M1234567890*",
                "-downstream" => "YYYYYYYYYYYYYYYYYYYY",
                "-expected-usable" => 0
               };

  push @tests, {
                "-desc" => "2 matches, only 1 long enough",
                "-upstream" => "M12345*XXXXXXXXX",
                "-downstream" => "YYYM12345678901",
                "-expected-usable" => 1
               };

  push @tests, {
                "-desc" => "too short, and in upstream",
                "-upstream" => "M123456789*",
                "-downstream" => "YYYYYYYYYYYYYYYYYYYY",
                "-expected-usable" => 0
               };
  # to add:
  # - example w/multiple matches but only one good one

  foreach my $test (@tests) {
    my ($status) = find_qualifying_coding_frame(
                                                %{$test},
                                                "-no-sanity" => 1
                                               );
    my $exp = $test->{"-expected-usable"};
    dump_die($test, "no expected usable status") unless defined $exp;
    printf "%s: expected %d, got %d\n", $test->{"-desc"}, $exp, $status;
    die "sanity fail" unless $exp == $status;
  }
}

sub cluster_jobs {
  my $infiles = get_infiles();
  die "no files" unless @{$infiles};

  unique_outfile_check($infiles);
  foreach my $infile (@{$infiles}) {
    my $outfile = get_outfile($infile);

    my $cmd = $clr->get_command_line("-single" => $infile);
    my $c = new Cluster();
    $c->node_class("");
    $c->memory_reserve_mb($CLUSTER_RAM);
    $c->memory_limit_mb($CLUSTER_RAM);
    $c->outfile($outfile);
    $c->project("SV_inframe_check");
    $c->command($cmd);
    $c->run();
  }
}

sub unique_outfile_check {
  my ($infiles) = @_;

  my %saw_out;
  foreach my $cicero (@{$infiles}) {
    # ensure outfiles are unique before starting run
    # (don't want to crash in the middle of a long list)
    my $outfile = get_outfile($cicero);
    die "duplicate outfile $outfile; consider using -fq?" if $saw_out{$outfile} and not($FLAGS{"ignore-duplicate-outfile"});
    $saw_out{$outfile} = 1;
  }
}

sub hsp_filter {
  my ($hsp, $hsp_problems) = @_;
  $hsp_problems = {} unless $hsp_problems;

  die "hit not on + strand" unless $hsp->strand() == 1;
  my $usable = 1;

  if ($hsp->gaps() > $BLAT_AA_MAX_GAPS) {
    $usable = 0;
    $hsp_problems->{too_many_gaps} = 1;
  }

  my $mismatches = $hsp->length("query") - $hsp->num_identical();
  if ($mismatches > $BLAT_AA_MAX_MISMATCHES) {
    $usable = 0;
    $hsp_problems->{too_many_mismatches} = 1;
  }

  unless ($hsp->num_identical() >= $BLAT_AA_MIN_IDENTICAL) {
    $usable = 0;
    $hsp_problems->{insufficient_identity} = 1;
  }
  return $usable;
}

sub add_exon_annotation {
  #
  #  add exon annotation (affected exon number, total exons)
  #  for in-frame events
  #
  my ($row, $pair_info, $r2g) = @_;

  #
  #  identify exon for both geneA and geneB:
  #
  foreach my $end_code (qw(A B)) {
    my $key = sprintf 'gene_%s_nm', lc($end_code);
    my $nm = $pair_info->{$key} || die;
    my $chrom = get_chrom($row, $end_code) || die;
    my $bn = get_base_number($row, $end_code) || die;

    my $set_raw = $r2g->find(
                             "-accession" => $nm,
                             "-reference" => $chrom,
                            );
    # a transcript may be mapped to multiple chroms and/or positions

    my @set_passed;
    foreach my $hit (@{$set_raw}) {
      if ($hit->{aa_trans_status} > 0 or
          # validated mapping
          $hit->{aa_trans_status} == RefFlat2Genome::STATUS_ERROR_AA_MISMATCH_MINOR
         ) {
        # usable

        my $map = $hit->{codon_map} || die;

        my $seek_dir;
        # genomic seek direction to anchor the mapping to an exon,
        # if necessary.
        my $strand = $hit->{strand} || die;
        if ($end_code eq "A") {
          # for GeneA, seek upstream
          $seek_dir = $strand eq "+" ? -1 : 1;
        } elsif ($end_code eq "B") {
          # for GeneB, seek downstream
          $seek_dir = $strand eq "+" ? 1 : -1;
        } else {
          die;
        }

        my %cm;
        # index genomic positions
        my %exons;
        foreach my $entry (@{$map}) {
          $cm{$entry->{ref_base_num}} = $entry;
          $exons{$entry->{exon_number}} = 1 if $entry->{exon_number};
        }

        $pair_info->{sprintf 'gene_%s_exon_count', lc($end_code)} = max(keys %exons);

        my $final_pos;

        my $cm_min = min(keys %cm);
        my $cm_max = max(keys %cm);
        my $max_tries = $cm_max - $cm_min;

        for (my $try = 0; $try <= $max_tries; $try++) {
          # find the nearest exon
          my $search_pos = $bn + ($try * $seek_dir);
          printf STDERR "code %s, try %d, pos %d\n", $end_code, $try, $search_pos if $VERBOSE;
          if (my $entry = $cm{$search_pos}) {
            # hit a mapped base
            dump_die($entry) unless $entry->{exon_number};
            $final_pos = $search_pos;
            last;
          }
        }

        my $exno = -1;
        if ($final_pos) {
          my $entry = $cm{$final_pos};
          $exno = $entry->{exon_number} || die;
        }
        my $key = sprintf 'gene_%s_exon_number', lc($end_code);
        $pair_info->{$key} = $exno;
      }
    }
  }

  #
  #  generate "pretty" event descriptions:
  #
  #   KMT2A:NM_005933[1,6]-AFF1:NM_001166693[8,21]
  #    (KMT2A e6 -> AFF1 e8)
  #
  my $genea_exon = $pair_info->{gene_a_exon_number} || 0;
  my $geneb_exon = $pair_info->{gene_b_exon_number};
  my $geneb_end = $pair_info->{gene_b_exon_count};

  my $genea_nm = $pair_info->{gene_a_nm} || die;
  my $geneb_nm = $pair_info->{gene_b_nm} || die;

  my $genea_gene = $r2g->get_gene($genea_nm) || $genea_nm;
  my $geneb_gene = $r2g->get_gene($geneb_nm) || $geneb_nm;
  # might not be an annotation if using non-refSeq library

  my $desc = sprintf '%s:%s[%d,%d]-%s:%s[%d,%d]',
    $genea_gene, $genea_nm, 1, $genea_exon,
      $geneb_gene, $geneb_nm, $geneb_exon, $geneb_end;
  $pair_info->{inframe_desc} = $desc;
}

sub junction_setup {
  #
  # synthesize data required for SV analysis from a junction.
  #
  my ($row, $ga, $r2g) = @_;

  my $broken;

  my $strand;
  # TO DO:
  # - option to accept a user-specified strand, if specified, don't
  #   derive ourselves from transcripts

  #
  #  derive strand from matching transcripts:
  #
  my $j = parse_junction($row->{junction} || dump_die($row, "no junction field"));
  my %strands;
  if (my $strand = $row->{junction_strand}) {
    # use strand specified by user hint
    $strands{$strand} = 1;
    die unless $strand eq "+" or $strand eq "-";
  } else {
    # autodetect
    foreach my $type (qw(start end)) {
      my $chr = $j->{$type}{ref_name};
      my $pos = $j->{$type}{base_number};

      my $strands = find_nms_from_ga(
                                     "-ga" => $ga,
                                     "-chrom" => $chr,
                                     "-position" => $pos,
                                     "-strand" => 1
                                    );
      foreach (@{$strands}) {
        $strands{$_} = 1;
      }
    }
  }

  my @strands = keys %strands;

  if (@strands == 0) {
    $broken = "no strand info found";
  } elsif (@strands > 1) {
    $broken = "matching transcripts hit multiple strands";
  } else {
    $strand = $strands[0];
  }

  my $sv_ort;
  my $contig;
  if ($broken) {
    $sv_ort = "?: $broken";
    $contig = "error";
  } else {
    $sv_ort = ">";
    my $chrom_a = $j->{start}{ref_name} || die;
    my $chrom_b = $j->{end}{ref_name} || die;

    my $pos_start = $j->{start}{base_number} || die;
    my $pos_end = $j->{end}{base_number} || die;

#    printf "%s\n", join ".", $chrom, $pos_start, $pos_end;

    my $rs_a = $r2g->get_reference_sequence($chrom_a);
    my $rs_b = $r2g->get_reference_sequence($chrom_b);
    my $before = substr($$rs_a, $pos_start - $JUNCTION_FLANKING_SEQUENCE, $JUNCTION_FLANKING_SEQUENCE);
    my $after = substr($$rs_b, $pos_end - 1, $JUNCTION_FLANKING_SEQUENCE);

    $contig = $before . $after;
    my ($pos_a, $pos_b);
    if ($strand eq "+") {
      $pos_a = $pos_start;
      $pos_b = $pos_end;
    } elsif ($strand eq "-") {
      $contig = reverse_complement($contig);
      $pos_a = $pos_end;
      $pos_b = $pos_start;
    } else {
      die "unhandled strand value";
    }

    $row->{chrA} = $chrom_a;
    $row->{chrB} = $chrom_b;

    $row->{ortA} = $row->{ortB} = $strand;

    $row->{featureA} = $row->{featureB} = "unknown";
    # since we don't know at this point, mark unknown to allow
    # extra wiggle room during searching if it's in an intron

    $row->{posA} = $pos_a;
    $row->{posB} = $pos_b;
    $row->{contig} = $contig;
  }
  $row->{contig} = $contig;
  $row->{sv_ort} = $sv_ort;
  die unless $row->{sv_ort};

}

sub parse_junction {
  # FIX ME: move into shared utility lib?
  my ($j) = @_;
  my @f = split /,/, $j;
  die unless @f == 2;
  my %result;
  foreach my $ref (
                   [ "start", $f[0] ],
                   [ "end", $f[1] ]
                  ) {
    my ($tag, $spec) = @{$ref};
    my @g = split /:/, $spec;
    die "$spec must be :-delimited and have 3 fields" unless @g == 3;
    my ($ref_name, $base_number, $strand) = @g;
    $result{$tag} = {
                     ref_name => $ref_name,
                     base_number => $base_number,
                     strand => $strand
                    };
  }
  return \%result;
}

sub genea_frame_check {
  #
  # check geneA frame status independently
  #
  my (%options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $frame_full = $options{"-frame"} || die "-frame";
  my $ga = $options{"-ga"} || die "-ga";
  my $r2g = $options{"-r2g"} || die "-r2g";
  my $annot_global = $options{"-annot-global"} || die "-annot-global";
  my $annot_trouble = $options{"-annot-trouble"} || die;
  my $transcripts_a_ref = $options{"-results"} || die;

  my $a_frame_fragment = substr($frame_full, 0, int(length($frame_full) * .55));
  # excerpt of the frame for geneA only.  Since the contig sequence
  # for the junction was simulated by us, we know the 5' portion
  # is from geneA.  This trimming is required otherwise we may
  # inappropriately match to the geneB portion of the sequence.
  # Allow a little extra to ensure the entire sequence through the
  # junction site is present in the sequence.

  my $a_nms = find_nms_from_ga(
                               "-ga" => $ga,
                               "-row" => $row,
                               "-which" => "A",
                              );

  my (%in_frame, %not_in_frame);

  foreach my $a_nm (@{$a_nms}) {
    my $set_passed = get_transcript_mappings(
                                             "-end" => "A",
                                             "-accession" => $a_nm,
                                             "-row" => $row,
                                             "-r2g" => $r2g,
                                             "-annot-global" => $annot_global,
                                             "-annot-trouble" => $annot_trouble
                                            );

    my %frame_blat;

    if (@{$set_passed}) {
      my ($hit, $codon_index, $codon_number, $final_pos) = @{$set_passed->[0]};

      my $a_aa = $r2g->get_aa("-accession" => $a_nm) || die;
      printf STDERR "A AA: %s\n", $a_aa;

      my %frames;
      my $fid = "peptide_frame_geneA";
      # excerpt of frame for geneA only
      $frames{$fid} = $a_frame_fragment;

      my $blat = get_blat();
      $blat->verbose(1);

      my $parser = $blat->blat(
                               "-query" => {
                                            $a_nm => $a_aa
                                           },
                               "-database" => \%frames,
                               "-protein" => 1,
                              );
      my $result = $parser->next_result;
      # one object per query sequence (only one query seq)
      if ($result) {
        while (my $hit = $result->next_hit()) {
          # hits from this query to a database sequence
          # (Bio::Search::Hit::HitI)
          push @{$frame_blat{$hit->name()}}, $hit;
        }
      }

      my ($idx, $fragment, $frame_error, $utr_codons, $blat_mode) = frame_search(
                                         "-frame" => $a_frame_fragment,
                                         "-hit" => $hit,
                                         "-index" => $codon_index,
                                         "-length" => $MIN_TUPLE_AA,
                                         "-direction" => -1,
                                         # upstream back into geneA

                                         "-blat" => $frame_blat{$fid},
                                         "-end" => "A",
                                                                                );
      # find unique anchoring in AA.  The fragment size will
      # be expanded for uniqueness if necessary.
      # the formal codon number is taken from the RefFlat2Genome mapping,
      # however the string lookup is required to separate the
      # GeneA portion from the GeneB portion for lookup for GeneB.
      #
      # the fragment will NOT include the codon at the breakpoint
      # as this might not match Cicero's contig sequence (see
      # comments in get_aa_fragment()), and will be 1 AA away in the
      # appropriate direction.

      if ($frame_error) {
        # special case: no upstream fragment!  e.g. in 1st codon,
        # how to handle??
        $annot_trouble->add($frame_error, $a_nm);
        # will be reported multiple times; good example for
        # converting to hash model for this field
      } elsif ($idx == -1) {
        # not in-frame
        $not_in_frame{$a_nm} = 1;
      } else {
        #
        #  in-frame
        #
        $in_frame{$a_nm} = 1;
#        $annot_global->add("TEST ME: geneA blat required for inframe call") if $blat_mode;
        $transcripts_a_ref->{$a_nm} = 1;
        #        last;
        # don't stop since we need to identify all inframe accessions
      }
    }
  }

  my $a_inframe;
  if (%in_frame) {
    $a_inframe = 1;
  } elsif (%not_in_frame) {
    $a_inframe = 0;
  } else {
    $a_inframe = "";
  }
  return $a_inframe;
}

sub get_transcript_mappings {
  my (%options) = @_;
  my $nm  = $options{"-accession"} || die;
  my $r2g = $options{"-r2g"} || die;
  my $row = $options{"-row"} || die;
  my $end = $options{"-end"} || die;
  die "end must be A or B" unless $end eq "A" or $end eq "B";
  my $chrom = get_chrom($row, $end) || die;
  my $bn = get_base_number($row, $end) || die;

  my $ort = get_orientation($row, $end) || die;
  my $feature = get_feature($row, $end) || die;
  my $annot_trouble = $options{"-annot-trouble"} || die;
  my $annot_global = $options{"-annot-global"} || die;

  my $set_raw = $r2g->find(
                           "-accession" => $nm,
                           "-reference" => $chrom,
                          );
  # a transcript may be mapped to multiple chroms and/or positions on
  # the chrom.

  my @set_passed;
  printf STDERR "raw set count: %d\n", scalar @{$set_raw} if $VERBOSE;
  foreach my $hit (@{$set_raw}) {
    unless ($hit->{aa_trans_status} > 0) {
      # only use transcripts where AA mapping to genome has been validated
      $annot_trouble->add(sprintf("gene%s_CDS_mapping_exception", $end),
                          sprintf "%s:%s", $nm, $hit->{aa_trans_status});
      if ($hit->{aa_trans_status} == RefFlat2Genome::STATUS_ERROR_AA_MISMATCH_MINOR) {
        # occasional small mismatches, e.g. NM_000014 which has
        # a single-AA mismatch.  Continue anyway, though this could
        # potentially lead to AA matching problems.  Failure to anchor
        # to GeneA will still be caught and reported anyway.
      } elsif ($FLAGS{"force-refflat-map"}) {
        printf STDERR "debug: tolerating broken CDS mapping\n";
      } else {
        # don't process other exception types, more serious
        next;
      }
    }

    if ($ort ne $hit->{strand}) {
      $annot_global->add(sprintf("gene%s_transcript_strand_disagree", $end),
                         $hit->{name});
      # might be nothing: maybe a transcript unrelated to the SV call?
    }
    my $seek_dir;
    my $strand = $hit->{strand};
    die unless $strand eq "+" or $strand eq "-";
    if ($end eq "B") {
      # for geneB, search downstream in the transcript mapping
      $seek_dir = $strand eq "+" ? 1 : -1;
    } elsif ($end eq "A") {
      # for geneA, search upstream in the transcript mapping
      $seek_dir = $strand eq "+" ? -1 : 1;
    } else {
      die;
    }

    #          dump_die($hit);

    my $map = $hit->{codon_map} || die;
    my %cm;
    # index because we may need to seek
    my $ei = 0;
    foreach my $entry (@{$map}) {
      $cm{$entry->{ref_base_num}} = $entry;
      $entry->{entry_index} = $ei++;
#      dump_die($entry, "debug", 1);
    }

    my $max_tries;
    my $cm_min = min(keys %cm);
    my $cm_max = max(keys %cm);
    if ($feature eq "coding") {
      # if annotated as coding, require the site be pretty close
      $max_tries = $CODON_SEEK_MAX_TRIES;
    } elsif ($bn >= $cm_min and $bn <= $cm_max) {
      # site is within transcript: allow broader search to
      # account for long introns
      $max_tries = $cm_max - $cm_min;
    } else {
      # allow some flexibility outside of the transcript if we want to
      # e.g. extend beyond the UTR.  We should be much stricter though
      # to prevent accidentally hitting multiple transcripts mapped to
      # the same chrom.  For example chr15.82729181 hits within one
      # mapping for NM_001291420 but is also within the genomic
      # mapping length of another:
      #
      # GOLGA6L9        NM_001291420    chr15   +       82722987        82731586       82726274 82728611        8       82722987,82724030,82725017,82725762,82726233,82727068,82728166,82728478,        82723099,82724150,82725098,82725850,82726802,82727129,82728269,82731586,
      # GOLGA6L9        NM_001291420    chr15   +       82804768        83108111       83102799 83105136        8       82804768,82805815,83101542,83102287,83102758,83103593,83104691,83105003,        82804880,82805935,83101623,83102375,83103327,83103654,83104794,83108111,
      #
      # without this restriction there would be a match to the 2nd transcript
      # as well because it is within the range of $cm_min - $cm_max.
      #
      $max_tries = $CODON_SEEK_MAX_TRIES_BEYOND_TRANSCRIPT;
      printf STDERR "site %s.%d outside of transcript mapping for %s (%d-%d), limiting search to %d nt\n", $chrom, $bn, $hit->{name}, $cm_min, $cm_max, $max_tries;
    }

    my $found;
    for (my $try = 0; $try <= $max_tries; $try++) {
      # require that GeneB touch a base handled by parsed and
      # verified exon mapping (including simulated UTR codons).
      # CICERO blat coordinates may not be accurate, in which
      # case we need to seek downstream in geneB transcript.
      # (if we were anchoring to GeneA first, we'd seek
      # upstream.)
      my $search_pos = $bn + ($try * $seek_dir);
#      printf STDERR "try %d, pos %d\n", $try, $search_pos;

      if (my $entry = $cm{$search_pos}) {
        # hit a mapped base
        push @set_passed, [ $hit, $entry->{entry_index}, $entry->{codon_number}, $search_pos, $entry->{coding_base_number} ];
        $found = 1;
        last;
        # found site
      }
    }

    if (not($found) and $r2g->is_intronic(
                                          "-entry" => $hit,
                                          "-base-number" => $bn
                                         )) {
      # we have a good transcript->genome mapping, but target
      # site was not found.  Look for evidence site is intronic.
      $annot_global->add(sprintf("gene%s_intronic", $end), $nm);
      # record in case e.g. this is the only transcript
    }
  }

  if (@set_passed > 1) {
    printf STDERR "ERROR: site touched more than once in same NM_!\n";
    foreach my $r (@set_passed) {
      printf STDERR "index:%d codon_number:%d search_pos:%d\n", @{$r}[1,2,3];
      dump_die($r->[0], "debug", 1);
    }
    die;
  }

#  printf STDERR "found index:%d codon_number:%d search_pos:%d\n",
#    @{$set_passed[0]}[1,2,3] if @set_passed;

  return \@set_passed;
}

sub mask_frames_for_gene_a {
  #
  # when performing the initial anchoring of geneB, blat may be required
  # in some cases.  When using BLAT we want to match against the geneB
  # portion of the frame only.  This is necessary e.g. for junctions
  # in the same gene, where the geneA portion is in-frame but the geneB
  # portion is not: blat may match to the geneA portion which is inappropriate.
  # As a workaround, mask the geneA portion of the frame AA before blat.
  #
  my (%options) = @_;
  my $row = $options{"-row"} || die;
  my $frames = $options{"-frames"} || die;

  my %masked;
  foreach my $fid (keys %{$frames}) {
    my $raw = $frames->{$fid};
    my $mask_until;
    if ($FLAGS{"junction-mode"}) {
      # our generated contig sequence, 50-50 GeneA/B
      $mask_until = int(length($raw) * 0.45);
    } else {
      #
      # CICERO SV:
      #
# From: Li, Yongjin
# Sent: Thursday, August 28, 2014 10:30 AM
# To: Edmonson, Michael
# Subject: RE: CICERO in-frame checks
#
# The columns qposA and qposB are breakpoints in the contig. Again, there might be some base shift. In most cases, the breakpoint should be the same. But for fusions with insertion, like EPOR fusion, C11orf95-RELA and C11orf95-MAML2, the breakpoints qposA and qposB are different.
#
# You may use (qposB+10) to extract contig sequence mapped to geneB, and match the AA sequence to geneB.

      my $qposa = $row->{qposA};
      unless ($qposa) {
        if ($FLAGS{"generate-qpos"}) {
          $qposa = int((length($raw) * 3) / 2);
          # nucleotide contig coordinates
        } else {
          die "no qposA field";
        }
      }

      $mask_until = int((($qposa - 10) / 3) - 1);
      # position minus recommended wiggle room, converted to AA index,
      # -1 for good measure
#      die sprintf "mask range error: qposA=%d %s", $qposa, join ",", $mask_until, length($raw) if $mask_until <= 0 or $mask_until >= length $raw;
    }
    my $masked = $raw;
    for (my $i = 0; $i < $mask_until; $i++) {
      substr($masked, $i, 1) = "X";
      # mask sequence rather than trimming it to preserve existing code,
      # which uses hit index to determine codon number
    }
    printf STDERR "masked frame %s: %s\n", $fid, $masked if $VERBOSE;
    $masked{$fid} = $masked;
  }
  return \%masked;
}

sub gene_feature_annot {
  my (%options) = @_;
  my $row = $options{"-row"} || die;
  my $which = $options{"-which"} || die;
  my $ga = $options{"-ga"} || die;

  my $rows = find_nms_from_ga(
                              "-ga" => $ga,
                              "-row" => $row,
                              "-which" => $which,
                              "-return-rows" => 1
                             );

  my %all_types;
  my $rf = new RefFlatFile();

  if (@{$rows}) {
    my $base_num = get_base_number($row, $which) || die;
    foreach my $rf_row (@{$rows}) {
      my $type = $rf->get_annotation_for_position(
                                                  "-row" => $rf_row,
                                                  "-base" => $base_num,
                                                 ) || die;
      my $acc = $rf_row->{name} || die;
      $all_types{$type}{$acc} = 1;
    }
  }

  $all_types{intergenic}{no_transcripts_found} = 1 unless %all_types;
  # no rows (or rows with usable strand)


  my @things;
  foreach my $type (sort keys %all_types) {
    my @entries = sort keys %{$all_types{$type}};
    push @things, join "=", $type, join ",", @entries;
  }
  my $result = join ";", @things;
  return $result;
}

sub get_config_entry {
  my ($param, $ckey) = @_;

  my $fn = $FLAGS{$param};
  unless ($fn) {
    my $genome = $FLAGS{genome} || die sprintf "specify -%s FILE or -genome VERSION", $param;
    my $config_genome = TdtConfig::readConfig('genome', $genome) || die "can't find config for $genome";
    $fn = $config_genome->{$ckey} || die "no config entry for $ckey";
  }
  confess "no $param file" unless $fn;
  confess "$fn?" unless -s $fn;
  return $fn;
}

sub get_refflat_file {
  return get_config_entry("refflat", "REFSEQ_REFFLAT");
}

sub get_fasta_file {
  return get_config_entry("fasta", "FASTA");
}

sub get_r2p_file {
  return get_config_entry("refgene2protein", "REFGENE2PROTEIN");
}

sub patch_coding_base_numbers {
  my $infiles = get_infiles();
  my $r2g = get_r2g();
  foreach my $infile (@{$infiles}) {
    my $outfile = ($FLAGS{fq} ? $infile : basename($infile)) . ".patch.tab";

    my $df = new DelimitedFile("-file" => $infile,
                               "-headers" => 1,
                              );
    my $rpt = $df->get_reporter(
                                "-file" => $outfile,
                                "-extra" => [
                                             qw(
                                               sv_refseqA_coding_base_number
                                               sv_refseqB_coding_base_number
                                               sv_refseqB_last_coding_base_number
                                              )
                                            ]
                               );

    # while (my $row = $df->next("-ref" => 1)) {  # headerless
    while (my $row = $df->get_hash()) {
      patch_row_coding_bases("-row" => $row,
                             "-r2g" => $r2g);
      $rpt->end_row($row);
    }
    $rpt->finish();
  }
}

sub patch_row_coding_bases {
  # add coding base # annotations to SVs in an output row.
  # implemented this way so can be used to patch fields
  # from existing reports (vs. requiring a re-run)
  my (%options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $r2g = $options{"-r2g"} || die "-r2g";

  foreach my $end (qw(A B)) {
    my $f = "sv_refseq" . $end;
    die "no field $f" unless exists $row->{$f};
    my @basenum;
    my @basenum_cds_end;

    my $field = $row->{$f};
    my @d = $field =~ /,/g;
    my $count = @d + 1;
    # might be trailing blank fields, e.g. if reporting on geneB event only
    #
    # sv_refseqA: NM_005517,
    # sv_refseqB: NM_005966,NM_005966

    foreach my $acc (split /,/, $field, $count) {
      # may be blank
      my $basenum = "";
      my $basenum_end = "";
      if ($acc) {
        # find nearest coding base for each accession
        # - main: RefGene2Genome.pm
        # - possible backup method: RefFlatFile::get_base_translations
        #   (simpler range-based model, could be modified to add cds #)
        ($basenum, $basenum_end) = get_nearest_cds_base(
                                        "-end" => $end,
                                        "-accession" => $acc,
                                        "-row" => $row,
                                        "-r2g" => $r2g,
                                       );
      }
      push @basenum, $basenum || "";
      push @basenum_cds_end, $basenum_end || "";
    }

    my $k = sprintf "sv_refseq%s_coding_base_number", $end;
    $row->{$k} = join ",", @basenum;

    $k = sprintf "sv_refseq%s_last_coding_base_number", $end;
    $row->{$k} = join ",", @basenum_cds_end;
  }
}

sub get_r2g {
  my $refgene_fn = get_refflat_file();
  my $r2p_fn = get_r2p_file();

  my $cwc = $FLAGS{"cache-whole-chrom"} || 0;
  my $lcc = $FLAGS{"limit-chrom-cache"} || 0;
  my $ltc = $FLAGS{"limit-transcript-cache"};

  my $r2g = new RefFlat2Genome(
                               "-fasta" => get_fasta_file(),
			       "-cache_whole_chrom" => $cwc,
			       "-limit_chrom_cache" => $lcc,
			       "-limit_transcript_cache" => $ltc,
                              );
  $r2g->cache_genome_mappings(1) if $CACHE_TRANSCRIPT_MAPPINGS;

  $r2g->verbose(1) if $FLAGS{"verbose-r2g"};

  my $infer_aa = $FLAGS{"infer-aa"};
  $r2g->infer_aa(1) if $infer_aa;

  $r2g->parse(
              "-refflat" => $refgene_fn,
              "-refgene2protein" => $r2p_fn
             );
  return $r2g;
}


sub get_nearest_cds_base {
  my (%options) = @_;
  my $nm  = $options{"-accession"} || die;
  my $r2g = $options{"-r2g"} || die;
  my $row = $options{"-row"} || die;
  my $end = $options{"-end"} || die;
  die "end must be A or B" unless $end eq "A" or $end eq "B";
  my $chrom = get_chrom($row, $end) || die;
  my $bn = get_base_number($row, $end) || die;

  my $ort = get_orientation($row, $end) || die;
  my $feature = get_feature($row, $end) || die;

  my $set_raw = $r2g->find(
                           "-accession" => $nm,
                           "-reference" => $chrom,
                          );
  # a transcript may be mapped to multiple chroms and/or positions on
  # the chrom.

  my @set_passed;
  foreach my $hit (@{$set_raw}) {
    next unless $hit->{aa_trans_status} > 0;
    # only use transcripts where AA mapping to genome has been validated
    next if $ort ne $hit->{strand};
    # do transcripts ever legitimately doubly-map in opposite orientations??

    my $map = $hit->{codon_map} || die;
    my %cm;
    # index
    foreach my $entry (@{$map}) {
      $cm{$entry->{ref_base_num}} = $entry;
    }

    my $max_tries;
    my $cm_min = min(keys %cm);
    my $cm_max = max(keys %cm);
    if ($feature eq "coding") {
      # if annotated as coding, require the site be pretty close
      $max_tries = $CODON_SEEK_MAX_TRIES;
    } elsif ($bn >= $cm_min and $bn <= $cm_max) {
      # site is within transcript: allow broader search to
      # account for long introns
      $max_tries = $cm_max - $cm_min;
    } else {
      # allow some flexibility outside of the transcript if we want to
      # e.g. extend beyond the UTR.  We should be much stricter though
      # to prevent accidentally hitting multiple transcripts mapped to
      # the same chrom.  For example chr15.82729181 hits within one
      # mapping for NM_001291420 but is also within the genomic
      # mapping length of another:
      #
      # GOLGA6L9        NM_001291420    chr15   +       82722987        82731586       82726274 82728611        8       82722987,82724030,82725017,82725762,82726233,82727068,82728166,82728478,        82723099,82724150,82725098,82725850,82726802,82727129,82728269,82731586,
      # GOLGA6L9        NM_001291420    chr15   +       82804768        83108111       83102799 83105136        8       82804768,82805815,83101542,83102287,83102758,83103593,83104691,83105003,        82804880,82805935,83101623,83102375,83103327,83103654,83104794,83108111,
      #
      # without this restriction there would be a match to the 2nd transcript
      # as well because it is within the range of $cm_min - $cm_max.
      #
      $max_tries = $CODON_SEEK_MAX_TRIES_BEYOND_TRANSCRIPT;
#      printf STDERR "site %s.%d outside of transcript mapping for %s (%d-%d), limiting search to %d nt\n", $chrom, $bn, $hit->{name}, $cm_min, $cm_max, $max_tries;
    }

    my $found_entry;
    for (my $try = 0; $try <= $max_tries; $try++) {
      # require that GeneB touch a base handled by parsed and
      # verified exon mapping (including simulated UTR codons).
      # CICERO blat coordinates may not be accurate, in which
      # case we need to seek downstream in geneB transcript.
      # (if we were anchoring to GeneA first, we'd seek
      # upstream.)
      foreach my $seek_dir (1, -1) {
        my $search_pos = $bn + ($try * $seek_dir);
        #      printf STDERR "try %d, pos %d\n", $try, $search_pos;

        last if $found_entry = $cm{$search_pos};
        # hit a mapped base
      }
      last if $found_entry;
    }
    if ($found_entry) {
      my $strand = $hit->{strand};
      my @genomic = sort {$a <=> $b} keys %cm;
      my ($search_start, $search_dir, $search_end);
      my $last_coding_base_number;
      if ($strand eq "+") {
        # stop codon should be first coding entry from the end
        $search_start = $#genomic;
        $search_dir = -1;
        $search_end = 0;
      } elsif ($strand eq "-") {
        # stop should be first coding entry from the start
        $search_start = 0;
        $search_dir = 1;
        $search_end = @genomic;
      } else {
        die;
      }

      for (my $i = $search_start; $i != $search_end; $i += $search_dir) {
        my $e = $cm{$genomic[$i]};
#        dump_die($e, "debug", 1);
        if ($e->{in_cds} and !$e->{in_utr}) {
          $last_coding_base_number = $e->{coding_base_number};
          printf STDERR "last coding base number in %s/%s: %d\n", $hit->{name}, $strand, $last_coding_base_number if $VERBOSE;
          last;
        }
      }
      $found_entry->{last_coding_base_number} = $last_coding_base_number || die "can't find last coding base";
      # hack

      $found_entry->{distance} = abs($bn - ($found_entry->{ref_base_num} || dump_die($found_entry, "no ref_base_num in found entry")));
      push @set_passed, $found_entry;
      dump_die($hit, "saving hit", 1) if $VERBOSE;
    }
  }

  my $coding_base_number;
  my $last_coding_base_number;

  @set_passed = sort {$a->{distance} <=> $b->{distance}} @set_passed;
  # if multiple mappings for the same accession match, use the one
  # where the result is closest to the original target site
  # /home/medmonso/work/yongjin/bugs/2016_06_07_inframe/crash.txt

  if (@set_passed) {
    $coding_base_number = $set_passed[0]->{coding_base_number};
    $last_coding_base_number = $set_passed[0]->{last_coding_base_number};
  }
  return ($coding_base_number, $last_coding_base_number);
}

sub salvage_intergenic_genes {
  my (%options) = @_;
  my $row = $options{"-row"} || die;
  my $end = $options{"-which"} || die;
  die unless $end eq "A" or $end eq "B";
  my $nms = $options{"-nms"} || die;
  my $rff = $options{"-rff"} || die;

  my $any_added = 0;

  unless ($nms and @{$nms}) {
    # only if list is not already populated
    my $feature = $row->{"feature" . $end};
    if ($feature eq CICERO_FEATURE_INTERGENIC) {
      my $chrom = get_chrom($row, $end);
      my $gene = $row->{"gene" . $end};
      my $pos = $row->{"pos" . $end};
      my $rows = $rff->find_by_gene($gene);
      if ($rows and (@{$rows})) {
        my %distance;
        foreach my $r (@{$rows}) {
          if (chrom_compare($chrom, $r->{chrom} || die)) {
            # mapping to same chromosome as SV
            my $tx_start = $r->{txStart} || die;
            my $tx_end = $r->{txEnd} || die;
            my $distance;
            if ($pos < $tx_start) {
              $distance = $tx_start - $pos;
            } else {
              $distance = $pos - $tx_end;
            }
            push @{$distance{$distance}}, $r unless $r->{name} =~ /NR_/;
            # skip non-coding
          }
        }
        if(%distance) {
          my ($closest) = sort {$a <=> $b} keys %distance;
          foreach my $r (@{$distance{$closest}}) {
            my $nm = $r->{name} || die;
            push @{$nms}, $nm;
            $any_added = 1;
          }
        }
      }
    }
  }

  return $any_added;
}

sub chrom_compare {
  my ($c1, $c2) = @_;
  foreach ($c1, $c2) {
    s/^chr//;
  }
  return $c1 eq $c2;
}

sub intronic_frame_search {
  my (%options) = @_;
  my $frame = $options{"-frame"} || die;
  my $hit = $options{"-hit"} || die;
  my $codon_index = $options{"-index"};
  die unless defined $codon_index;
  my $length = $options{"-length"} || die;
  my $direction = $options{"-direction"} || die;
  my $end = $options{"-end"} || die;
  my $annot_global = $options{"-annot-global"} || confess("need -annot-global");
  my $row = $options{"-row"} || confess("-row");
  my $r2g = $options{"-r2g"} || confess("-r2g");

  my $map = $hit->{codon_map} || die;
  my @ref_bases = map {$_->{ref_base_num} || die} @{$map};

  my %map;
  foreach my $r (@{$map}) {
    $map{$r->{ref_base_num} || die} = $r;
  }

  my $chrom = get_chrom($row, $end) || die;
  my $sv_pos = get_base_number($row, $end) || die;
  my $genome = $r2g->fai->get_sequence("-id" => $chrom);

  my ($idx, $fragment);

  if ($r2g->is_intronic(
                        "-entry" => $hit,
                        "-base-number" => $sv_pos,
                       )) {
    my $strand = $map->[1]->{ref_base_num} > $map->[0]->{ref_base_num} ? "+" : "-";
    dump_die($map->[$codon_index], "debug", 1);
#    die join ",", $map->[0]->{ref_base_num}, $map->[1]->{ref_base_num};

    my $exon_edge = $map->[$codon_index] || die;
    my $last_codon_base = $exon_edge->{codon_base} || die;
    my $last_codon_number = $exon_edge->{codon_number} || die;

    my $exon_edge_base = $exon_edge->{ref_base_num};
    my $genomic_direction = $exon_edge_base < $sv_pos ? 1 : -1;
    my $reverse_direction = $genomic_direction == 1 ? -1 : 1;
    my $pos = $exon_edge_base + $genomic_direction;
    if ($genomic_direction < 0) {
      die unless $strand eq "+";
    } else {
      die unless $strand eq "-";
    }

    # TO DO:
    # - clone edge entry plus a few?

    #
    # move from exon edge towards target (SV site)
    #

    my %map;

    my %encode;
    my %codon2aa;

    while (1) {
      my %r;
      $r{ref_base_num} = $pos;
      $r{ref_base_genome} = substr($$genome, $pos - 1, 1);
      my $tbase = $r{ref_base_genome};
      if ($strand eq "+") {
      } else {
        $tbase = reverse_complement($tbase);
      }
      $r{ref_base_transcript} = $tbase;

      my $codon_base = $last_codon_base - 1;
      # since we are always moving back from exon edge,
      # codon base (and fake codon number) will always decrease
      if ($codon_base < 1) {
        $codon_base = 3;
        $last_codon_number--;
      }
      $last_codon_base = $codon_base;
      $r{codon_base} = $codon_base;
      $r{codon_number} = $last_codon_number;

      $encode{$last_codon_number}{$codon_base} = $tbase;

      printf STDERR "%s\n", join " ", @r{qw(ref_base_num ref_base_genome ref_base_transcript codon_base codon_number)};

      $map{$pos} = \%r;

      last if ($pos == $sv_pos);
      $pos += $genomic_direction;
    }

    my $ct = new Bio::Tools::CodonTable();
    foreach my $codon_number (sort {$a <=> $b} keys %encode) {
      if (scalar(keys %{$encode{$codon_number}}) == 3) {
        my $codon = join "", @{$encode{$codon_number}}{qw(1 2 3)};
        $codon2aa{$codon_number} = $ct->translate($codon);
      }
    }

    if (1) {
      print STDERR "dump fake codons:\n";
      foreach my $pos (sort keys %map) {
        my $ref = $map{$pos};
        printf STDERR "%s\n", join " ",
          $ref->{ref_base_num},
            $ref->{codon_number},
              $codon2aa{$ref->{codon_number}},
                $ref->{ref_base_genome},
                  $ref->{ref_base_transcript},
                    $ref->{codon_base};
          }
    }

    # start at SV site
    # IGNORE that codon
    # continue towards known edge gathering codons
    my $sv_codon_number = $map{$sv_pos}{codon_number} || die;
    $pos = $sv_pos;
    my @aa;

    my %saw;
    while (1) {
      my $entry = $map{$pos} || last;
      my $codon_number = $entry->{codon_number};
      unless ($codon_number == $sv_codon_number) {
        # skip codon number associated with site as it may contain
        # some sequence from other end of SV and so a different AA
        unless ($saw{$codon_number}) {
          push @aa, $codon2aa{$codon_number} || die;
          last if @aa == $length;
          $saw{$codon_number} = 1;
        }
      }
      $pos += $reverse_direction;
    }

    $fragment = join "", @aa;
    $idx = index($frame, $fragment);
  }

  if (defined($idx) and $idx != -1) {
    die $frame;
    $annot_global->add(sprintf "gene%s_intronic_anchor", $end);
  }

  return ($idx, $fragment);
}

sub hack2 {
  my $r2g = get_r2g();
  my $genome = $r2g->fai->get_sequence("-id" => 2);

  my $edge = 29446394 - 20;
  my $posB = 29447684 - 20;

  my $len = $posB - $edge;

  my $chunk = substr($$genome, $edge, $len);
  $chunk = reverse_complement($chunk);
  # to form AA sequence at exon edge

  my $frames = get_frames($chunk);
  foreach my $id (sort keys %{$frames}) {
    my $frame = $frames->{$id};
    printf STDERR "%s (%d): %s\n", $id, length($frame), $frame;
  }
}

sub frame_offset {
  my ($fid) = @_;
  $fid =~ /frame_offset_(\d+)$/ || die;
  my $offset = $1;
  die unless $offset >= 0 and $offset <= 2;
  return $offset;
}

sub match_a {
  #
  # look for upstream fragment in geneA
  #
  my (%options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $a_nm = $options{"-nm"} || die "-nm";
  my $upstream = $options{"-upstream"};
  die "-upstream" unless defined $upstream;
  # may be blank
  my $r2g = $options{"-r2g"} || die "-r2g";
  my $annot_trouble = $options{"-annot-trouble"} || die;
  my $annot_global = $options{"-annot-global"} || die;
  my $infer_aa = $FLAGS{"infer-aa"};
  my $try_utr5 = $options{"-try-utr5"};

  my $a_codon_number = "";
  my $is_inframe = "";
  # 0 is reserved for when we have enough information
  # for a reasonable guess
  my $a_anchor_type = "";
  my $a_upstream_idx_end;
  # index in upstream sequence (geneA) *after* match

  my $a_aa;
  if ($try_utr5) {
    my $annot_null = new AnnotationField("-unique" => 1);
    my $set_passed = get_transcript_mappings(
                                             "-end" => "A",
                                             "-accession" => $a_nm,
                                             "-row" => $row,
                                             "-r2g" => $r2g,
                                             "-annot-global" => $annot_null,
                                             "-annot-trouble" => $annot_null
                                            );
    if (@{$set_passed}) {
      my ($hit, $codon_index, $codon_number, $final_pos) = @{$set_passed->[0]};
      $a_aa = $r2g->get_utr5_aa_from_map("-hit" => $hit);
    } else {
      printf STDERR "ERROR: no passed AA for %s\n", $a_nm;
      $a_aa = $r2g->get_aa(
                           "-accession" => $a_nm,
                          );
    }
  } else {
    $a_aa = $r2g->get_aa(
                         "-accession" => $a_nm,
                        );
#    die "WTF" unless $a_aa;
  }
  printf STDERR "UTR5:%d AA:%s\n", $try_utr5, ($a_aa || "");
  # may not be a refgene2protein entry, e.g. NM_033557 (suppressed)

  if (not($a_aa) and $infer_aa) {
    # non-refSeq transcripts where we don't know the AA
    # sequence ahead of time.  Full parsing is required
    # to initialize these.
    my $annot_null = new AnnotationField("-unique" => 1);
    my $set_passed = get_transcript_mappings(
                                             "-end" => "A",
                                             "-accession" => $a_nm,
                                             "-row" => $row,
                                             "-r2g" => $r2g,
                                             "-annot-global" => $annot_null,
                                             "-annot-trouble" => $annot_null
                                            );
    $a_aa = $r2g->get_aa("-accession" => $a_nm);
  }

  unless ($a_aa) {
    $annot_trouble->add("geneA_no_refseq_AA_possibly_suppressed", $a_nm);
    return ($a_codon_number, $is_inframe, $a_anchor_type, $a_upstream_idx_end);
  }
  printf STDERR "geneA AA: %s\n", $a_aa if $VERBOSE;

  my $idx = index($a_aa, $upstream);

  if ($idx != -1) {
    # entire upstream sequence before geneB match
    # matches geneA's AA
    $a_codon_number = $idx + 1 + length($upstream);
    $is_inframe = SV_MATCH_CODE_INFRAME;
    $a_anchor_type = "perfect";
    $a_upstream_idx_end = length($upstream);
  }

  unless ($is_inframe) {
    #
    #  try #2:
    #
    # downstream AA fragment does not perfectly match.
    # try tuple-based lookup.
    # fail: try other methods?  e.g. if GeneB has some UTR
    # or might legitmately be out of frame
    # - or: BLAT
    # - BLAT might be a good idea: subtle differences in AA
    #   might cause lookup to fail even if event is in-frame
    # ALTERNATIVE: why not just blat geneB AA vs OTHER FRAMES?
    # - don't we expect GeneB to hit one of the frames?
    #   (qc check)
    # - putting it another way: if GeneB can't be blat'd
    # against one of the 3 frames, isn't that an error?
    # - the risk of continuing to do it the way we've been
    #   doing is that any failure is treated as the event being
    #   out of frame!
    my $us_len = length $upstream;

    if ($us_len < $MIN_TUPLE_AA) {
      $annot_trouble->add("upstream_AA_too_short_to_anchor_to_geneA");
    }

    for (my $i = $us_len - $MIN_TUPLE_AA;
         $i >= 0; $i--) {
      my $tuple = substr($upstream, $i, $MIN_TUPLE_AA);
      my $idx = index($a_aa, $tuple);
      printf STDERR "geneA search tuple: %s, idx=%d\n", $tuple, $idx if $VERBOSE;
      if ($idx != -1) {
        # hit
        my $idx2 = index($a_aa, $tuple, $idx + 1);
        if ($idx2 != -1) {
          # AMBIGUOUS TUPLE, NEEDS WORK,
          # but can be tricky, see
          # gene_b_ambiguous_tuple_tricky_run_of_S.tab where
          # extending the tuple fails beyond the S run  :(
          # extend tuple as above?
          # this is important because we want to report
          # the affected codon # unambiguously
          printf STDERR "ambiguous tuple %s!!\n", $tuple;
          $annot_global->add("ERROR_geneA_ambiguous_tuple_anchor", $a_nm);
          # record problem but continue processing;
          # maybe a later tuple will be be okay
        } else {
          $is_inframe = SV_MATCH_CODE_INFRAME;
          $a_codon_number = $idx + 1 + $MIN_TUPLE_AA;
          # first match to canonical AA in GeneA near breakpoint
          # (approximate codon at breakpoint site)
          $a_anchor_type = "tuple";
          $a_upstream_idx_end = $i + $MIN_TUPLE_AA;
          last;
        }
      }
    }
  }

  unless ($is_inframe) {
    # try BLATing the upstream sequence vs. the in-frame frame.
    # a small mismatch in the downstream sequence vs. the
    # AA for geneA may disrupt the tuple lookup process.

    my $blat = get_blat();
    # QUESTION:
    # should $frame_raw be trimmed after 1st stop codon?
    # (i.e. not really in frame if event aligns after a stop)

    my $parser = $blat->blat(
                             "-query" => {
                                          "contig_upstream" => $upstream
                                         },
                             "-database" => {
                                             "geneA_AA" => $a_aa
                                            },
                             "-protein" => 1,
                            );
    # blat sequence downstream of in-frame GeneA hit vs. GeneB AA
    my $result = $parser->next_result();
    # one result object per query sequence

    if ($result) {
      my @hits;
      my @hsp;
      my %hsp_problems;

      while (my $hit = $result->next_hit()) {
        # hits from this query to a database sequence
        # (Bio::Search::Hit::HitI)
        push @hits, $hit;
        while (my $hsp = $hit->next_hsp) {
          printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
            $hsp->score,
              $hsp->strand,
                $hsp->range("query"),
                  $hsp->range("hit"),
                    $hsp->num_identical(),
                      $hsp->frac_identical("query"),
                        $hsp->length("query"),
                          $hsp->length("hit"),
                            $hsp->length("total"),
                              $hsp->query_string(),
                                $hsp->hit_string();

          push @hsp, $hsp if hsp_filter($hsp, \%hsp_problems);
        }
      }

      my @blat_problems;

      push @blat_problems, "multiple_hits" if @hits > 1;
      push @blat_problems, "multiple_hsp" if @hsp > 1;
      # could be e.g. an ambiguous tuple, which would
      # also fail above

      if (@blat_problems) {
        # serious issue: ambiguity, etc.
        $annot_trouble->add("ERROR_geneA_blat_exception", \@blat_problems);
        $a_codon_number = "";
        $a_anchor_type = "";
      } elsif (@hsp) {
        # acceptable match
        $is_inframe = SV_MATCH_CODE_INFRAME;
        $a_codon_number = ($hsp[0]->range("hit"))[0];
        # point at which blat hit starts, which might
        # exclude e.g. some inserted sequence
        $a_anchor_type = "blat";

        my @query_range = $hsp[0]->range("query");
        $a_upstream_idx_end = $query_range[1];
      } else {
        # not handled
        $annot_trouble->add("ERROR_geneA_inframe_failed_but_has_blat_hit", $a_nm);
        # not specifically diagnosed
        $annot_trouble->add("geneA_blat_issues", [ sort keys %hsp_problems ]) if %hsp_problems;
        $a_codon_number = "";
        $a_anchor_type = "";
      }
    } else {
      # no BLAT hit, so we now have some confidence
      # this is out of frame
      # QUESTION: test low complexity cases like
      # gene_b_ambiguous_tuple_tricky_run_of_S.tab.frame.tab
      # where tuple is a run of S?
      $is_inframe = SV_MATCH_CODE_FRAMESHIFT;
      $a_codon_number = "";
      $a_anchor_type = "";
    }
  }                                # blat

  if ($try_utr5) {
    # codon numbering won't work if based on UTR5 matching
    $a_codon_number = "";
  }

#  die join ",", $a_codon_number, $is_inframe, $a_anchor_type, $a_upstream_idx_end;

  return ($a_codon_number, $is_inframe, $a_anchor_type, $a_upstream_idx_end);

}

sub build_cdl {
  my ($rows, $field) = @_;
  my @v;
  foreach my $r (@{$rows}) {
    my $v = $r->{$field};
    $v = "" unless defined $v;
    push @v, $v;
  }
  return join ",", @v;
}

sub report_full_length {
  my %options = @_;
  my $row_src = $options{"-row"} || die;
#  dump_die($row_src, "Debug", 1);

  my $eval = $options{"-eval"} || die;
  my $rpt = $options{"-rpt"} || die;
  my $nm2protein = $options{"-nm2protein"} || die;
  my $r2g = $options{"-r2g"} || die;

  foreach my $row (@{$eval}) {
    dump_die($row, "debug", 1);
    next unless $row->{inframe_call};
    # only applicable for in-frame calls

    my @notes;

    my %r;
    foreach my $f (qw(sample contig)) {
      # fields from CICERO row
      $r{$f} = $row_src->{$f} || dump_die($row_src, "no data for $f");
    }

    my $fail;
    foreach my $f (qw(gene_a_nm gene_b_nm frame_raw)) {
      # fields from transcript-specific row
#      $r{$f} = $row->{$f} || dump_die($row, "no data for $f");
      $r{$f} = $row->{$f};
      $fail = 1 unless $r{$f};
    }
    next if $fail;
    # geneA info might be missing for e.g. in-frame call code 2

    $r{contig_frame} = $row->{frame_raw};

    $r{geneA} = $r2g->get_gene($row->{gene_a_nm}) || $row->{gene_a_nm};
    $r{geneB} = $r2g->get_gene($row->{gene_b_nm}) || $row->{gene_b_nm};
    # gene symbols for this transcript pairing

    my $geneA_protein = $nm2protein->{$row->{gene_a_nm}} || die;
    my $geneB_protein = $nm2protein->{$row->{gene_b_nm}} || die;
    $r{AA_gene_a_full} = $geneA_protein;
    $r{AA_gene_b_full} = $geneB_protein;

    my $interstitial = $row->{interstitial_AA};
#    dump_die($row, "interstitial AA, test me!!") if $interstitial;
    $r{AA_interstitial} = $interstitial;
    push @notes, "contains_interstitial_sequence" if $interstitial;

#    my $geneA_upstream = $row->{"frame_upstream"} || die;
    # large fragment from contig.  Note this includes geneA sequence
    # PLUS ANY INTERSTITIAL SEQUENCE.
    my $idx_geneA_end = $row->{"gene_a_index"};
    my $idx_geneB_start = $row->{"gene_b_index"};
    die unless defined $idx_geneA_end and defined $idx_geneB_start;
    my $geneA_upstream = substr($row->{"frame_upstream"}, 0, $idx_geneA_end + 1);
    # the upstream sequence, EXCLUDING any interstitial AA

    my $geneB_downstream = $row->{"gene_b_aa_downstream"} || die;
    # small fragment: anchor to larger piece?
    my $geneB_downstream_extended = $row->{gene_b_aa_downstream_extended};
    # only longer if site is in 5' UTR

    my $geneA_upstream_anchor = $geneA_upstream;

    # require anchor matches exactly once.  however the contig may go
    # out of sync upstream, possibly because it's genomic sequence
    # and so intronic regions won't match mRNA from refFlat record.
    if (unique_in_string($geneA_protein, $geneA_upstream_anchor)) {
      # success
    } elsif (my $chunk = unique_in_substring($geneA_protein, $geneA_upstream_anchor)) {
      # a substring uniquely matches geneA, record trimming required
      $geneA_upstream_anchor = $chunk;
      push @notes, sprintf "geneA_anchor_trim_bases=%d", index($geneA_upstream, $chunk);
    }

    # if (length($geneA_upstream_anchor) < $FULL_LENGTH_MIN_GENE_A_ANCHOR) {
    unless (unique_in_string($geneA_protein, $geneA_upstream_anchor)) {
      my $msg = "ERROR: can't anchor to geneA, geneA_upstream=$geneA_upstream geneA=$geneA_protein";
      push @notes, $msg;
      foreach my $f (qw(
    			 AA_gene_a
    			 AA_gene_b
    			 full_length_sequence
    			 full_length_sequence_length
    			 full_length_sequence_contains_entire_gene_b
    			 QC_span
    		      )) {
    	$r{$f} = "";
      }
      $r{notes} = join ",", @notes;

      $rpt->end_row(\%r);
      next;
    }

    $geneA_protein =~ /^(.*$geneA_upstream_anchor)/i || die;
    my $upstream_full = $1;
    $r{AA_gene_a} = $upstream_full;

#    die "geneB anchor fail for $geneB_downstream in $geneB_protein" unless unique_in_string($geneB_protein, $geneB_downstream);
    # short so might not be unique?

    my @ahits = $geneB_protein =~ /${geneB_downstream}/ig;
#    printf STDERR "ahits=%d\n", scalar @ahits;

    my $full_sequence;
    if (unique_in_string($geneB_protein, $geneB_downstream)) {
      $geneB_protein =~ /($geneB_downstream.*)$/i || die;
      my $downstream_full = $1;

      $r{AA_gene_b} = $downstream_full;

      $full_sequence = join "", $upstream_full, $interstitial, $downstream_full;
    } elsif (my $substr = unique_in_substring($geneB_protein, $geneB_downstream)) {
      # the downstream sequence doesn't perfectly match geneB, but a substring
      # does: possibly leading 5' UTR or intronic sequence, or anchoring
      # requiring blat rather than a substring
      $geneB_protein =~ /($substr.*)$/i || die;
      my $downstream_full = $1;
      $r{AA_gene_b} = $downstream_full;

      my $idx = index($geneB_downstream, $substr) || die;
      $interstitial .= substr($geneB_downstream, 0, $idx);
      $r{AA_interstitial} = $interstitial;

      $full_sequence = join "", $upstream_full, $interstitial, $downstream_full;
    } elsif ($row->{gene_b_codon_number} < 1) {
      #
      # geneB site is upstream of start / possibly 5' UTR
      #
      my $feature_b = $row_src->{featureB} || die;
      if ($feature_b eq "5utr") {
	push @notes, "geneB_5_prime_UTR";
      } elsif ($feature_b eq "intron") {
	push @notes, "geneB_intron";
      } else {
	push @notes, "geneB_potential_intron_or_5_prime_UTR";
	push @notes, sprintf "cicero_geneB_feature=%s", $feature_b;
      }

      my $idx_up = index($row->{frame_raw}, $geneA_upstream_anchor);
#      dump_die($row, "can't anchor geneA $geneA_upstream_anchor idx=$idx_up") unless $idx_up == 0;
      die if $idx_up == -1;

      my $frame_down = substr($row->{frame_raw}, $idx_up + length($geneA_upstream_anchor));
      # downstream portion of the frame AFTER geneA, including any interstitial
      # sequence, and presumably at some point the start of geneB

      my $dl = length($frame_down);

      for (my $i = 0; $i < $dl; $i++) {
	my $chunk = substr($frame_down, $i);
	$chunk = substr($chunk, 0, $MIN_TUPLE_AA) if length($chunk) > $MIN_TUPLE_AA;
	last if length($chunk) < 7;
	# hack

	if (unique_in_string($geneB_protein, $chunk)) {
#	  printf STDERR "unique hit for $chunk\n";
	  my $interstitial = "";
	  $interstitial = substr($frame_down, 0, $i) if ($i > 0);

	  $geneB_protein =~ /($chunk.*)$/i || die;
	  my $downstream_full = $1;
	  $r{AA_gene_b} = $downstream_full;

#	  printf STDERR "FINAL:\n%s\n", join "\n", $upstream_full, $interstitial, $downstream_full;

	  $r{AA_interstitial} = $interstitial;
	  $full_sequence = join "", $upstream_full, $interstitial, $downstream_full;
	  last;
	}
      }

      unless ($full_sequence) {
	push @notes, "geneB_direct_frame_coding_anchor_failed";
	# can't anchor to geneB coding sequence directly from frame sequence
	my $b_start_frag = substr($geneB_protein, 0, $MIN_TUPLE_AA);
	my $bi = index($geneB_downstream_extended, $b_start_frag);
	die "geneB start anchor fail, can't find $b_start_frag in $geneB_downstream_extended" unless $bi >= 0;
	# ensure geneB coding start is found in downstream extended sequence
	die unless unique_in_string($geneB_downstream_extended, $b_start_frag);

	my $interstitial_pt1 = "";
	for (my $i=0; $i < length($frame_down); $i++) {
	  my $chunk = substr($frame_down, $i);
	  if (index($geneB_downstream_extended, $chunk) == 0) {
	    # done
	    last;
	  } else {
	    $interstitial_pt1 .= substr($frame_down, $i, 1);
	  }
	}

	my $interstitial_pt2 = substr($geneB_downstream_extended, 0, $bi);

	$r{AA_interstitial} = join "", $interstitial_pt1, $interstitial_pt2;

	$full_sequence = join "", $upstream_full, $interstitial_pt1, $interstitial_pt2, $geneB_protein;
	$r{AA_gene_b} = $geneB_protein;

	dump_die($row, sprintf("geneB UTR5 anchor, frame:%s geneB_down:%s down_ext:%s frame_down:%s geneB_full:%s codon:%s interstitial1:%s interstitial2:%s final_full:%s",
		 $row->{frame_raw},
		 $geneB_downstream,
		 $geneB_downstream_extended,
		 $frame_down,
		 $geneB_protein,
                 $row->{gene_b_codon_number},
		 $interstitial_pt1,
		 $interstitial_pt2,
			       $full_sequence
			      ), 1
		);
      }
#    } elsif (($geneB_protein =~ /${geneB_downstream}/ig) > 1) {
# WTF?  why doesn't this work??
    } elsif (@ahits > 1) {
      #
      # anchoring sequence is ambiguous in geneB
      #
      my $target_codon = $row->{gene_b_codon_number} || die;
      # codon number nearest to geneB position

      #
      # find index closest to the annotated codon #:
      #
      my $last = -1;
      my %distance;
      while (1) {
	my $idx = index($geneB_protein, $geneB_downstream, $last + 1);
	last if $idx == -1;
#	print "idx $idx\n";
	$distance{$idx} = abs($idx - $target_codon);
	$last = $idx;
      }
      my ($idx_nearest) = sort {$distance{$a} <=> $distance{$b}} keys %distance
;
      my $downstream_full = substr($geneB_protein, $idx_nearest);
      $r{AA_gene_b} = $downstream_full;

      $full_sequence = join "", $upstream_full, $interstitial, $downstream_full;
#      die join "\n", $row->{frame_raw}, $upstream_full, $interstitial, $downstream_full;
    } else {
      dump_die($row, "geneB anchor fail for $geneB_downstream in $geneB_protein");
    }

    $r{full_length_sequence} = $full_sequence;
    $r{full_length_sequence_length} = length $full_sequence;

    #
    # QC check:
    #

#    push @notes, "full_raw_frame_not_found" unless unique_in_string($full_sequence, $row->{frame_raw});
    # this is a desirable sanity check to make sure the CICERO
    # predicted contig actually appears in the final fusion sequence
    # (correctly incorporates interstitial sequence, etc.)  HOWEVER
    # I'm not sure if this is always reasonable: what if frame
    # contains genomic sequence that isn't in the refFlat?  Or if
    # there are SNVs affecting the protein sequence?

    my $flank = 10;
    my $qc_start = $idx_geneA_end - ($flank - 1);
    my $qc_end = $idx_geneB_start + ($flank - 1);

#    my $cnum = $row->{gene_b_codon_number};
#    if ($cnum < 0) {
#      $qc_end += abs($cnum) + 10;
      # if contig includes 5' UTR, extended until we actually hit coding
      # sequence of geneB
#    }

    my $qc_chunk = substr($row->{frame_raw}, $qc_start, ($qc_end - $qc_start) + 1);
    $r{QC_span} = $qc_chunk;

    my $qc_pass = 1;

    unless (unique_in_string($full_sequence, $qc_chunk)) {
      $qc_pass = 0;
      push @notes, "QC_span_frame_not_found";
    }
    # while this contains less geneA/geneB flanking sequence than the full frame,
    # no guarantee event might not happen e.g. near an exon edge and so still not work

    if ($r{AA_interstitial} =~ /\*/) {
      $qc_pass = 0;
      push @notes, "QC_interstitial_termination";
    }

    foreach my $f (qw(AA_gene_a AA_gene_b)) {
      dump_die(\%r, "stop found in gene sequence") if $r{$f} =~ /\*/;
    }

    if (my $is = $r{AA_interstitial}) {
      if (length($is) >= 5) {
	my $ok = 0;
	my $frame = $row->{frame_raw};
	my $il = length($is);
	for (my $i = $il; $i > 0; $i--) {
	  my $chunk = substr($is, 0, $i);
	  my $idx = index($frame, $chunk);
#	  printf STDERR "test %s: %d\n", $chunk, $idx;
	  if ($idx != -1) {
	    $ok = 1;
	    last;
	  }
	}
	dump_die(\%r, "fragment of final interstitial not found in raw frame") unless $ok;
      }
    }
    $r{QC_pass} = $qc_pass;

    $r{full_length_sequence_contains_entire_gene_b} = $r{AA_gene_b} eq $r{AA_gene_b_full} ? 1: 0;

    $r{notes} = join ",", @notes;
    $rpt->end_row(\%r);

  }
}

sub unique_in_string {
  my ($full, $substring) = @_;
  $substring =~ s/\*/\\\*/g;
  # prevent stop codon from being interpreted as regex character
  my @hits = $full =~ /$substring/ig;
  return @hits == 1;
}

sub digest_full_length {
  #
  # - reannotate input report recording reasons to keep/reject?
  # - record "best" per sample? (just the longest??)
  # - record "best" overall
  # - record if same event observed in > 1 sample
  #
  my ($infile) = $FLAGS{"digest-full-length"};
  my $df = new DelimitedFile("-file" => $infile,
			     "-headers" => 1,
			     );

  my %all;
  my %sample2protein;
  my %protein2sample;
  my %protein2geneA;
  my %protein2geneB;
  my %fusions;
  while (my $row = $df->get_hash()) {

    next if $FULL_LENGTH_DIGEST_EXCLUDE_INTERSTITIAL_STOP and
      $row->{AA_interstitial} =~ /\*/;

    my $seq = $row->{full_length_sequence} || next;
    # some intra-fusion events can't (yet) produce these
    $all{$seq} = 1;
    my $sample = $row->{sample} || die;
    my $geneA = $row->{geneA} || die;
    my $geneB = $row->{geneB} || die;
    my $fusion = join "_", $geneA, $geneB;
    $fusions{$fusion} = 1;

    my $geneA_nm = $row->{gene_a_nm} || die;
    my $geneB_nm = $row->{gene_b_nm} || die;

    $sample2protein{$sample}{$seq} = join "_", $geneA_nm, $geneB_nm;
    # might be > 1 NM pairing, just save one

    $protein2sample{$seq}{$sample} = 1;
    $protein2geneA{$seq}{$geneA} = 1;
    $protein2geneB{$seq}{$geneB} = 1;
  }
  die "multiple fusions in input file" if scalar keys %fusions > 1;
  my ($fusion) = keys %fusions;

  printf STDERR "merge samples...\n";
  my %best_in_sample;
  foreach my $sample (sort keys %sample2protein) {
    my ($best, $worst_frac) = find_best_protein($sample2protein{$sample});
#    $best_in_sample{$sample} = $best;
    $best_in_sample{$best} = $worst_frac;
  }

  printf STDERR "merge global...\n";
  my ($best_overall, $best_overall_worst_frac) = find_best_protein(\%all);

  my %starts = map {substr($_, 0, $FULL_LENGTH_START_COMPARE_LENGTH), 1} keys %all;
  my $unique_start_sequences = scalar keys %starts;

  my $rpt = new Reporter(
			 "-file" => $infile . ".digest.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   fusion
					   sample
					   sample_count
					   geneA
					   geneB
					   full_length
					   length
					   unique_start_sequences
					   best_in_sample
					   best_in_sample_blat_identity
					   best_overall
					   best_overall_blat_identity
					)
				      ],
			 "-auto_qc" => 1,
			);

  foreach my $full (sort keys %all) {
    my %r;
    $r{fusion} = $fusion;
    $r{full_length} = $full;
    $r{length} = length $full;
    $r{sample} = join ",", sort keys %{$protein2sample{$full}};
    $r{sample_count} = scalar keys %{$protein2sample{$full}};
    $r{geneA} = join ",", sort keys %{$protein2geneA{$full}};
    $r{geneB} = join ",", sort keys %{$protein2geneB{$full}};
    $r{unique_start_sequences} = $unique_start_sequences;

    my $best_in_sample = 0;
    my $best_in_sample_blat_identity = "";
    if (my $worst_frac = $best_in_sample{$full}) {
      $best_in_sample = 1;
      $best_in_sample_blat_identity = $worst_frac;
    }
    $r{best_in_sample} = $best_in_sample;
    $r{best_in_sample_blat_identity} = $best_in_sample_blat_identity;

    my $is_best_overall = 0;
    my $best_overall_blat_identity = "";

    if ($full eq $best_overall) {
      $is_best_overall = 1;
      $best_overall_blat_identity = $best_overall_worst_frac;
    }
    $r{best_overall} = $is_best_overall;
    $r{best_overall_blat_identity} = $best_overall_blat_identity;

    $rpt->end_row(\%r);
  }

  $rpt->finish();


}

sub find_best_protein {
  #
  #  eliminate sequences that can be (near-)perfectly aligned to longer
  #  sequences:
  #
  my ($proteins_hash) = @_;

  my %full = %{$proteins_hash};
  # copy as we'll be modifying

  my $raw_count = scalar keys %full;

  my $worst_identity = 100;
  # in case there is only 1 sequence
  my $done = 0;
  while (!$done) {
    printf STDERR "checking %d sequences\n", scalar keys %full;
    $done = 1;
    my @full_by_len = sort {length($a) <=> length($b)} keys %full;
  IDX:
    for (my $idx = 0; $idx < @full_by_len - 1; $idx++) {
      my @idx_rest = ($idx + 1 .. @full_by_len - 1);
      my @rest = @full_by_len[@idx_rest];
#      my %rest = map {"idx" . $_, $full_by_len[$_]} @idx_rest;
      my %rest = map {$full{$_}, $_} @rest;

      my $query_seq = $full_by_len[$idx];
      my $qlen = length($query_seq);

      printf STDERR "compare %s/%s with %d: %s\n", $full{$query_seq}, $query_seq, scalar(@rest), join ", ", @rest;
#      dump_die(\%rest);
      my $blat = get_blat();
      $blat->verbose(1) if $VERBOSE;

      my $parser = $blat->blat(
			       "-query" => {
					    "query" => $query_seq,
					   },
			       #				   "-database" => $frames,
			       "-database" => \%rest,
			       "-protein" => 1,
			       "-dump-fasta" => 1,
			      );
      my $result = $parser->next_result;
      # one object per query sequence (only one query seq)
      if ($result) {
	my %hit2identity;
	# bucket results by target: the best result is not
	# necessarily the first one!

	while (my $hit = $result->next_hit()) {
	  # hits from this query to a database sequence
	  # (Bio::Search::Hit::HitI)

	  my $target_id = $hit->name();

	  printf STDERR "hits to %s, qlen=%d\n", $target_id, $qlen;
	  my @hsp;
	  while (my $hsp = $hit->next_hsp) {
	    printf STDERR "    score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	      $hsp->score,
		$hsp->strand,
		  $hsp->range("query"),
		    $hsp->range("hit"),
		      $hsp->num_identical(),
			$hsp->frac_identical("query"),
			  $hsp->length("query"),
			    $hsp->length("hit"),
			      $hsp->length("total"),
				$hsp->query_string(),
				  $hsp->hit_string();
	    push @hsp, $hsp;
	  }

	  my $num_identical;
	  if (@hsp > 1) {

	    my $iter = combinations(\@hsp, 2);
	    printf STDERR "raw: %s\n", join " ", @hsp;

	    # eliminate any hsp whose query range is a subset of another hsp:
	    my %bad;
	    my $hsp_seq = 0;
	    my $track_key = "__forever_unclean";
	    foreach my $hsp (@hsp) {
	      $hsp->{$track_key} = $hsp_seq++;
	    }

	    while (my $c = $iter->next()) {
	      my ($hsp1, $hsp2) = @{$c};
	      my $set1 = new Set::IntSpan join "-", $hsp1->range("query");
	      my $set2 = new Set::IntSpan join "-", $hsp2->range("query");
	      $bad{$hsp1->{$track_key}} = 1 if subset $set1 $set2;
	      $bad{$hsp2->{$track_key}} = 1 if superset $set1 $set2;
	    }
	    @hsp = grep {!$bad{$_->{$track_key}}} @hsp;

	    $num_identical = sum map {$_->num_identical()} @hsp;
	    # imperfect: what about partial overlaps???

	  } else {
	    $num_identical = $hsp[0]->num_identical();
	  }

	  my $frac_identical = $num_identical / $qlen;
	  # consider the full length of the query rather than just the
	  # aligned portion in the hsp(s)
	  die if $hit2identity{$target_id};
	  $hit2identity{$target_id} = $frac_identical;

	}  # next_hit()

	dump_die(\%hit2identity, "debug", 1);

	my @best = sort {$hit2identity{$b} <=> $hit2identity{$a}} keys %hit2identity;
	my $best_id = $best[0];
	my $frac_identical = $hit2identity{$best_id};
	if ($frac_identical >= $FULL_LENGTH_MIN_BLAT_IDENTITY_TO_MERGE) {
	  if (not(defined $worst_identity) or $frac_identical < $worst_identity) {
	    $worst_identity = $frac_identical;
	  }
	  $done = 0;
	  printf STDERR "deleting %s/%s, blat %.3f%% match to %s/%s\n", $full{$query_seq}, $query_seq, $frac_identical * 100, $full{$rest{$best_id}}, $rest{$best_id};
	  delete $full{$query_seq};
	  last IDX;
	} else {
	  printf STDERR "reject merge with %f identity\n", $frac_identical;
	}

      } else {
	die "no blat result for query";
      }
    }
  }

  die "returns array" unless wantarray();

  if (keys %full == 1) {
    # success
    return ((keys %full)[0], $worst_identity);
  } else {
    # couldn't collapse to one sequence with required strictness,
    # just use longest
    cluck sprintf "not exactly 1 left: orig=%d final=%d worst=%f", $raw_count, scalar keys %full, $worst_identity;
    $worst_identity = "merge_fail";

    my @sorted = sort {length($b) <=> length($a)} keys %full;

    return($sorted[0], "merge_fail");
  }

}

sub prep_excerpt_files {
  # extract input data for target fusions from a cicero output file
  my $f_all = $FLAGS{"prep-excerpt"} || die;

  # get target fusions first and limit all-rows loading to just those
  # with matching genes?

  #
  #  load list of fusions we're searching for and build list of all
  #  geneA's of interest:
  #
  my $f_fusions = $FLAGS{fusions} || die "-fusions";
  my $df_fusions = new DelimitedFile(
				     "-file" => $f_fusions,
				     "-headers" => 1,
				    );
  my %all_wanted_gene_a;
  my @search_fusions;
  while (my $row = $df_fusions->get_hash()) {
    push @search_fusions, $row;
    my $fusion = $row->{fusion} || dump_die($row, "fusion");
    my @f = split /_/, $fusion;
    die unless @f == 2;
    my ($gene_a_list, $gene_b_list) = @f;

    foreach my $g (split /,/, $gene_a_list) {
      $all_wanted_gene_a{$g} = 1;
    }
  }

  #
  #  load fusions to search against, saving only those of potential interest.
  #  much more efficient esp. if we eventually scan a large set of result
  #  files, etc.:
  #
  my @raw_rows;
  my $df_all = new DelimitedFile(
				 "-file" => $f_all,
				 "-headers" => 1,
				);
  while (my $row = $df_all->get_hash()) {
    my $geneA = $row->{geneA} || die;
    my @gene_a = split /,/, $row->{geneA};
    my ($usable) = grep {$all_wanted_gene_a{$_}} @gene_a;
    # not necessarily a match for both, but this dramatically reduces
    # the search space
    push @raw_rows, $row if $usable;
  }

  my @h_out;
  foreach (@{$df_all->headers_raw()}) {
    last if $_ eq "frame";
    push @h_out, $_;
  }

  my $rpt_summary = new Reporter(
			 "-file" => basename($f_fusions) . ".summary.tab",
			 "-delimiter" => "\t",
			 "-labels" => [
				       qw(
					   fusion
					   wanted
					   wanted_count

					   missing
					   missing_count
					   all_missing

					   additional
					   additional_count

					   missing_adjusted
					   missing_count_adjusted
					   additional_adjusted
					   additional_count_adjusted
					)
				      ],
			 "-auto_qc" => 1,
			);


  foreach my $row (@search_fusions) {
    #
    # search for each query fusion:
    #
    my $fusion = $row->{fusion} || dump_die($row, "fusion");
    printf STDERR "%s\n", $fusion;
    my @f = split /_/, $fusion;
    die unless @f == 2;
    my ($gene_from, $gene_to) = @f;

    my %wanted_samples;
    # all samples user is looking for
    my %needed_samples;
    # tracker: how many still left
    my %other_samples;
    # tracker: additional samples found
    if (my $sl = $row->{sample_list}) {
      %wanted_samples = map {$_, 1} split /,/, $sl;
      %needed_samples = %wanted_samples;
    }

    my $outfile = $fusion . ".input.tab";
    my $rpt = $df_all->get_reporter(
				"-file" => $outfile,
				"-auto_qc" => 1,
				"-extra" => [ "fusion_name" ],
				"-extra-prepend" => 1,
			       );
    $rpt->labels(\@h_out);

  # while (my $row = $df->next("-ref" => 1)) {  # headerless
    my $found = 0;
    foreach my $row (@raw_rows) {
      my %gene_a = map {$_, 1} split /,/, $row->{geneA};
      my %gene_b = map {$_, 1} split /,/, $row->{geneB};
      my @frame_calls = split /,/, $row->{frame};
      my $frame_usable = grep {$_ and $_ > 0} @frame_calls;

      if ($gene_a{$gene_from} and $gene_b{$gene_to}) {

	if ($frame_usable) {
	  $row->{fusion_name} = $fusion;
	  $rpt->end_row($row);
	  $found++;
	}

	# always perform these checks so we don't consider a record
	# missing if it's just not usable because it's not called in-frame
	my $sample = $row->{sample} || die "-sample";
	if (%wanted_samples) {
	  # sample tracking active
	  if ($wanted_samples{$sample}) {
	    delete $needed_samples{$sample};
	  } else {
	    $other_samples{$sample} = 1;
	  }
	}

      }
    }

    my %r;
    $r{fusion} = $fusion;
    $r{wanted} = join ",", sort keys %wanted_samples;
    $r{wanted_count} = scalar keys %wanted_samples;
    $r{missing} = join ",", sort keys %needed_samples;
    $r{missing_count} = scalar keys %needed_samples;

    $r{additional} = join ",", sort keys %other_samples;
    $r{additional_count} = scalar keys %other_samples;

    my %other_samples_adjusted = %other_samples;
    my %missing_adjusted;
    if ($r{missing_count}) {
      if ($r{additional_count}) {
	my %other_trimmed;
	foreach my $s_other (keys %other_samples) {
	  my $s = $s_other;
	  $s =~ s/_.*//;
	  $other_trimmed{$s} = $s_other;
	}

	foreach my $s_needed (keys %needed_samples) {
	  my $s = $s_needed;
	  $s =~ s/_.*//;
	  if (my $other = $other_trimmed{$s}) {
	    # some data from this sample observed, so don't consider missing
	    delete $other_samples_adjusted{$other};
	  } else {
	    # still can't find
	    $missing_adjusted{$s_needed} = 1;
	  }
	}
      } else {
	%missing_adjusted = %needed_samples;
      }
    }

    $r{missing_adjusted} = join ",", sort keys %missing_adjusted;
    $r{missing_count_adjusted} = scalar keys %missing_adjusted;
    $r{additional_adjusted} = join ",", sort keys %other_samples_adjusted;
    $r{additional_count_adjusted} = scalar keys %other_samples_adjusted;

    my $all_missing = 0;
    if ($found) {
      $rpt->finish();
    } else {
      printf STDERR "ERROR: no matches for $gene_from $gene_to\n";
      $all_missing = 1;
      # no matching samples found at all
    }
    $r{all_missing} = $all_missing;

    $rpt_summary->end_row(\%r);
  }

  $rpt_summary->finish();

}


sub unique_in_substring {
  my ($full_sequence, $downstream) = @_;

  my $dl = length($downstream);
  my $unique;
  for (my $i = 1; $i < $dl; $i++) {
    # start at 1 rather than 0 because full-length has been tried already
    my $chunk = substr($downstream, $i);
    printf STDERR "idx:%d chunk:%s full:%s\n", $i, $chunk, $full_sequence;

    if (unique_in_string($full_sequence, $chunk)) {
      $unique = $chunk;
      last;
    }
  }
  return $unique;
}

sub frame_debug {
  my $sequence = $FLAGS{"frame-debug"} || die "-frame-debug";
  get_frames($sequence);
}

