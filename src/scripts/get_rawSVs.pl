#!/usr/bin/env perl

## Exit codes:
## 1: Invalid arguments
## 2: Failed to combine internal events file (cat)
## 3: Failed to combine fusion events file (cat)
## 4: Failed to create main annotation directory (mkdir)
## 5: Failed to create temporary annotation directory (mkdir)
## 6: Breakpoint 1 had invalid orientation
## 7: Breakpoint 2 had invalid orientation
## 8: Failed to create Fasta file for assembly (cat)
## 9: Failed to run blat on first breakpoint fasta
## 10: Failed to run blat on second breakpoint fasta

use strict;
use warnings;

use Carp;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path;
use Cwd qw[abs_path];
use List::Util qw[min max];
use File::Glob ':globally';

use DelimitedFile;
use File::Temp qw/ tempdir /;

use CiceroSCValidator qw($lowqual_cutoff LEFT_CLIP RIGHT_CLIP);
use CiceroUtil qw(prepare_reads_file parse_range rev_comp
	is_PCR_dup read_fa_file get_discordant_reads get_sclip_reads normalizeChromosomeName exist_multiplename_checking);

require CiceroExtTools;

use TdtConfig;

use Transcript;
use Gene;
use GeneModel;

use constant BADFUSION_DISTANCE_CUTOFF => 200;

my $debug = 0;

# input/output
my ($genome, $ref_genome, $gene_model_file);
my ($out_dir, $gene_info_file, $junction_file, $known_itd_file);
my ($input_bam, $sample);
my ($all_output, $internal, $DNA) = (0, 0, 0);
my ($max_num_hits);
my ($excluded_gene_file, $excluded_fusion_file, $complex_region_file, $excluded_chroms, $gold_gene_file);
my ( $help, $man, $version, $usage );

if(@ARGV == 0){
	#TODO: get the correct usage string
	print STDERR "Usage: $0 -g <genome> -i <bam_file> -o <out_dir> -f <gene_info>\n";
	exit 1;
}

my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,	
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'  => \$genome,
	'ref_genome=s'  => \$ref_genome,
	'genemodel=s'		=> \$gene_model_file,
	'max_num_hits=i'	=> \$max_num_hits,
	'internal!' => \$internal,
	'all!' => \$all_output,
	'DNA!' => \$DNA,
	'f|gene_info_file=s' => \$gene_info_file,
	'known_itd_file=s'	=> \$known_itd_file,
	'j|junction_file=s' => \$junction_file,
	's|sample=s'		=> \$sample,
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
);

my $conf;
# Load configuration values for the genome build
if (&TdtConfig::findConfig("genome", $genome)){
	$conf = &TdtConfig::readConfig("genome", $genome);
}
else{
	croak("no config");
}

$ref_genome = $conf->{FASTA} unless($ref_genome && -f $ref_genome);
$gene_model_file = $conf->{'REFSEQ_REFFLAT'} unless($gene_model_file);
print STDERR "gene_model_file: $gene_model_file\n";
$excluded_gene_file = $conf->{EXCLUDED_GENES} unless($excluded_gene_file);
$excluded_fusion_file = $conf->{EXCLUDED_FUSIONS} unless($excluded_fusion_file);
print STDERR "excluded_fusion_file: $excluded_fusion_file\n";
$known_itd_file = $conf->{KNOWN_ITD_FILE} unless($known_itd_file);
print STDERR "KNOWN_ITD_FILE: ", $known_itd_file, "\n";
$gold_gene_file = $conf->{CICERO_GOLD_GENE_LIST_FILE};
print STDERR "CICERO_GOLD_GENE_LIST_FILE:", $gold_gene_file, "\n";
$excluded_chroms = $conf->{EXCLD_CHR} unless($excluded_chroms);
$complex_region_file = $conf->{COMPLEX_REGIONS} unless($complex_region_file);

# Load Cicero-specific configuration settings
$conf = &TdtConfig::readConfig('app', 'cicero');
$max_num_hits = $conf->{MAX_NUM_HITS} unless($max_num_hits);

# Assume sample name is the bam prefix
$sample = basename($input_bam, ".bam") unless($sample);

my @excluded_chroms = split(/,/,$excluded_chroms);

# Load list of complex regions to be used later to determine if
# breakpoints are located within a "complex" region
my @complex_regions;
#if ($complex_region_file && -s $complex_region_file){
	open (my $CRF, $complex_region_file) or die "Cannot open $complex_region_file: $!";
	while(<$CRF>){
		chomp;
		next if(/Start/);
		my ($name, $chr, $start, $end) = split(/\t/);
		my $cr = {
			name => $name,
			chr => $chr,
			start => $start,
			end => $end
			};
		push @complex_regions, $cr;
	}
	close($CRF);
	my %cr_hash = map { $_->{name} => 1 } @complex_regions;
#}


# Combine all of the individual results from Cicero.pl
my $unfiltered_file = "$out_dir/unfiltered.fusion.txt";
if($internal) {
	$unfiltered_file = "$out_dir/unfiltered.internal.txt";
	unless (-s $unfiltered_file){
		`cat $out_dir/*/unfiltered.internal.txt > $unfiltered_file`;
		if ($?){
			my $err = $!;
			print STDERR "Error combining internal events: $err\n";
			exit 2;
		}		
	}
}
else{
	unless (-s $unfiltered_file){
		`cat $out_dir/*/unfiltered.fusion.txt > $unfiltered_file`;
		if ($?){
			my $err = $!;
			print STDERR "Error combining fusion events: $err\n";
			exit 3;
		}
	}
}
print "Unfiltered fusions file: $unfiltered_file\n";

if (! $gene_info_file || ! -e $gene_info_file){
	my $out_prefix = basename($input_bam, ".bam");
	$gene_info_file = "$out_prefix.gene_info.txt";
	$gene_info_file = File::Spec->catfile($out_dir, $gene_info_file);
}

my %gene_info = ();
print STDERR "\ngene_info_file: $gene_info_file\n" if($debug);
open my $GI, "$gene_info_file" or die "cannot open < $gene_info_file: $!";
while(<$GI>){
	chomp;
	my ($name, $gRange, $strand, $mRNA_length, $cnt, $sc_cutoff) = split(/\t/);
	$gene_info{$name} = $sc_cutoff;
}
close $GI;


# Read the gene model from a refflat file
my $gm_format = "REFFLAT";
my $gm = GeneModel->new if($gene_model_file);
$gm->from_file($gene_model_file, $gm_format);

croak "Specify a genome: $ref_genome" unless (-f $ref_genome);

# Load the BAM file
my $sam_d = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $ref_genome);
my @seq_ids = $sam_d->seq_ids;

my %excluded = ();
open(my $EXC, "$excluded_gene_file") or die "Cannot open $excluded_gene_file: $!";
while(<$EXC>){
	my $line = $_;
	chomp($line);
	$excluded{$line} = 1;
}
close($EXC);

my %gold_genes = ();
if($gold_gene_file && -e $gold_gene_file){
	open(my $GGF, "$gold_gene_file");
	while(<$GGF>){
		my $line = $_;
		chomp($line);
		$gold_genes{$line} = 1;
	}
	close($GGF);
}

my %bad_fusions;
my $df = new DelimitedFile(
	"-file" => $excluded_fusion_file,
	"-headers" => 1,
);

while (my $row = $df->get_hash()) {
	my ($chrA, $posA, $chrB, $posB) = ($row->{chrA}, $row->{posA}, $row->{chrB}, $row->{posB});
	$bad_fusions{$chrA.":".$posA.":".$chrB.":".$posB} = 1;	      
}

sub is_bad_fusion{
	my ($chrA, $posA, $chrB, $posB) = @_;
	$chrA = "chr".$chrA unless($chrA =~ /chr/);
	$chrB = "chr".$chrB unless($chrB =~ /chr/);
	foreach my $xx (keys %bad_fusions) {
		my ($chr1, $pos1, $chr2, $pos2) = split(":",$xx);
		#print STDERR " badfusionlist |".$chr1.":".$pos1.":".$chr2.":".$pos2."\n";
		return 1 if($chrA eq $chr1 && $chrB eq $chr2 &&
			    abs($pos1 - $posA) < BADFUSION_DISTANCE_CUTOFF && abs($pos2 - $posB) < BADFUSION_DISTANCE_CUTOFF);# cutoff is based on the cutoff of merging GTEx false positive fusions from CICERO running
		return 1 if($chrA eq $chr2 && $chrB eq $chr1 &&
			    abs($pos2 - $posA) < BADFUSION_DISTANCE_CUTOFF && abs($pos1 - $posB) < BADFUSION_DISTANCE_CUTOFF);
	}
	return 0;
}

my %known_ITDs = ();
open(my $ITD_F, $known_itd_file) or die "Cannot open $known_itd_file: $!";
while(<$ITD_F>){
	chomp;
	my ($gene, $chr, $start, $end) = split(/\t/,$_);
	$known_ITDs{$gene} = [$start, $end];
}


my @raw_SVs = ();
my %gene_recurrance = ();
my %contig_recurrance = ();
my %genepairs = ();

### end of annotate.pl head

open(my $UNF, "$unfiltered_file") or die "Cannot open $unfiltered_file: $!";
print STDERR "Reading unfiltered file\n";
while(my $line = <$UNF>){
	chomp($line);
	my @fields = split("\t", $line);
	my ($gene1, $gene2) = ($fields[1], $fields[2]);
	my $cutoff = $fields[4];
	my $qseq = $fields[17];

	if($gene1 eq "NA"){
			my $crA = in_complex_region($fields[8], $fields[5]);
			if($crA){
				$gene1 = $crA;
			}
			else{
				$gene1 = join(":", $fields[8], $fields[5]);
			}
	}
	if($gene2 eq "NA"){
			my $crB = in_complex_region($fields[19], $fields[6]);
			if($crB){
				$gene2 = $crB;
			}
			else{
				$gene2 = join(":", $fields[19], $fields[6]);
			}
	}

	my @genes1 = split(/,|\|/, $gene1);
	my @genes2 = split(/,|\|/, $gene2);
	my $bad_gene = 0;

	# to remove genes with multiple potential partners.
	foreach my $g1 (@genes1) {
		$bad_gene = 1 if(exists($excluded{$g1}));
		foreach my $g2 (@genes2){
			$bad_gene = 1 if(exists($excluded{$g2}));
			next if ($g1 eq $g2);
			my $gene_pair = ($g1 le $g2) ? join(":",$g1,$g2) : join(":",$g2,$g1);
			next if(exists($genepairs{$gene_pair})); 
			$genepairs{$gene_pair} = 1;
			if(exists($gene_recurrance{$g1})){$gene_recurrance{$g1}++;}
			else{$gene_recurrance{$g1} = 1;}
			if(exists($gene_recurrance{$g2})){$gene_recurrance{$g2}++;}
			else{$gene_recurrance{$g2} = 1;}
		}
	}
	next if($bad_gene);

	my $bad_evidence = 0;
	my $first_bp = {
		reads_num => $fields[3],
		gene => $gene1,
		tpos => $fields[5],
		ort => $fields[7],
		tname => $fields[8],
		qstart => $fields[11],
		qend => $fields[12],
		qstrand => $fields[13],
		matches => $fields[14],
		percent => $fields[15],
		repeat => $fields[16]
	};
	$first_bp->{clip} = $first_bp->{qstrand}*$first_bp->{ort};
	my $qpos = $first_bp->{ort} > 0 ? $first_bp->{qend} : $first_bp->{qstart};
	if($qpos > 30 && $qpos < length($qseq) - 30){
		my $junc_seq = substr($qseq, $qpos - 30, 60);
		$bad_evidence += 2 if(low_complexity($junc_seq));
	}

	my ($gap, $same_gene) = ($fields[28], $fields[29]);
	next if(($same_gene && !$internal) || (!$same_gene && $internal));

	$bad_evidence += 1 if($gap >5);
	$bad_evidence += 1 if($first_bp->{matches} < 50);
	$bad_evidence += 1 if($first_bp->{repeat} > 0.9);
	$bad_evidence += 1 if($first_bp->{percent} < 0.95);

	my $second_bp = {
		reads_num => 0,
		gene => $gene2,
		tpos => $fields[6],
		ort => $fields[18],
		tname => $fields[19],
		qstart => $fields[22],
		qend => $fields[23],
		qstrand => $fields[24],
		matches => $fields[25],
		percent => $fields[26],
		repeat => $fields[27]
	};
	$second_bp->{clip} = $second_bp->{qstrand}*$second_bp->{ort};

        # If either breakpoint is in an excluded chromosome, skip it
	next if(is_bad_chrom($first_bp->{tname}) || is_bad_chrom($second_bp->{tname}));

	# Check if both breakpoints are in complex regions, if so skip
	my($crA, $crB) = (in_complex_region($first_bp->{tname}, $first_bp->{tpos}), in_complex_region($second_bp->{tname}, $second_bp->{tpos}));
	next if($crA && $crB);
	$first_bp->{gene} = $crA if($crA);
	$second_bp->{gene} = $crB if($crB);

	# Check if the reference uses 'chr' prefixes
	$first_bp->{tname} = normalizeChromosomeName($seq_ids[0], $first_bp->{tname});
	$second_bp->{tname} = normalizeChromosomeName($seq_ids[0], $second_bp->{tname});

	# Ensure that the breakpoint chromosome names match
	next if(!$internal && is_bad_fusion($first_bp->{tname}, $first_bp->{tpos}, $second_bp->{tname}, $second_bp->{tpos}));

	# Determine the variant type: CTX, Internal_inv, Interal_splicing, Internal_dup, ITX, read_through, DEL, INS
	my $type = get_type($first_bp, $second_bp, $same_gene);
	print STDERR "type: $type\n" if ($debug);
	# If we're not processing combined results (-all) and this is an Internal event
	if(!$all_output && $type =~ m/Internal/){
		if(%gold_genes){
			 next unless ($type eq 'Internal_dup' && exist_multiplename_checking(\%gold_genes, $first_bp->{gene}));
		}
		else{
			next unless ($type eq 'Internal_dup' && is_good_ITD($first_bp, $second_bp));
		}
		unless(is_good_ITD($first_bp, $second_bp)){
			next if(!$DNA && ($first_bp->{reads_num} < 5*$cutoff || $first_bp->{reads_num} < 10));
		}

		next if($type eq 'Internal_dup' && is_bad_fusion($first_bp->{tname}, $first_bp->{tpos}, $second_bp->{tname}, $second_bp->{tpos}));
	}

	my $tmp_SV = {
		first_bp => $first_bp,
		second_bp => $second_bp,
		};
	# If our SV is not a duplicate (same chromosome and within 10bp of the breakpoints)
	if($tmp_SV && ! is_dup_raw_SV(\@raw_SVs, $tmp_SV)){
		push @raw_SVs, $tmp_SV;

		# If we've seen this sequence string, increment
		# recurrance count
		if(exists($contig_recurrance{$qseq})){
			$contig_recurrance{$qseq}++;
		}else{ # Check the reverse complement as well.
			$qseq = rev_comp($qseq);
			if(exists($contig_recurrance{$qseq})){
				$contig_recurrance{$qseq}++;
			}else{$contig_recurrance{$qseq} = 1;}
		}
		$tmp_SV->{junc_seq} = $qseq;
		$tmp_SV->{type} = $type;
	}
}
close($UNF);

if($junction_file){
  print STDERR "Reading junction file\n";
  open(my $JUNC, "$junction_file");
  while(my $line = <$JUNC>){
	chomp($line);
	# chrX:1212 is the hg38 coordinate prefix for CRLF2 (chrX:1,190,449-1,212,815(GRCh38/hg38))
	# chrX:1331 is the hg19 coordinate prefix for CRLF2 (chrX:1,314,869-1,331,616(GRCh37/hg19))
	next unless($line =~ m/novel/ || $line =~ m/chrX:1212/ || $line =~ m/chrX:1331/);
	my @fields = split("\t",$line);
	my ($junction, $gene, $qc_flanking, $qc_perfect_reads, $qc_clean_reads) = @fields[0,3,5,8,9];
	unless($line =~ m/chrX:1212/ || $line =~ m/chrX:1331/){
		next if($qc_perfect_reads < 2 || $qc_flanking < 5);
		next if($qc_perfect_reads + $qc_clean_reads < 5);
	}
	my ($chr1, $pos1, $strand1, $chr2, $pos2, $strand2) = split(/:|,/, $junction);
	next if (abs($pos1 - $pos2) < 10000);
	next if(is_bad_chrom($chr1));
	next if(is_bad_fusion($chr1, $pos1, $chr2, $pos2));

	my($crA, $crB) = (in_complex_region($chr1, $pos1), in_complex_region($chr2, $pos2));
	next if($crA && $crB);

	my $cutoff = -1;
	my ($gene1,$gene2)= ("NA", "NA");

	foreach("+", "-"){
		my $gm_tree = $gm->sub_model($chr1, $_);
		last if(!defined($gm_tree));
		my ($cutoffA, $cutoffB) = (-1, -1);
		if($crA){
			$gene1 = $crA;
		}else{
			my @tmp = $gm_tree->intersect([$pos1 - 5000, $pos1 + 5000]);
			foreach my $g (@tmp){
				$g=$g->val;
				$cutoffA = $gene_info{$g->name} if($gene_info{$g->name} > $cutoffA);
				$gene1 = $gene1 eq "NA" ? $g->name : $gene1.",".$g->name;
			}
		}
		if($crB){
			$gene2 = $crB;
		}
		else{
			my @tmp = $gm_tree->intersect([$pos2 - 5000, $pos2 + 5000]);
			foreach my $g (@tmp){
				$g=$g->val;
				$cutoffB = $gene_info{$g->name} if($gene_info{$g->name} > $cutoffB);
				$gene2 = $gene2 eq "NA" ? $g->name : $gene2.",".$g->name;
			}
		}
		$cutoff = ($cutoffA + $cutoffB)/2  if($cutoffA + $cutoffB > 2*$cutoff);
	}
	next if($qc_perfect_reads < $cutoff);

	# Ensure that the breakpoint chromosome names match
	$chr1 = normalizeChromosomeName($seq_ids[0], $chr1);
	$chr2 = normalizeChromosomeName($seq_ids[0], $chr2);

	if($cutoff == -1){
		my $bg_reads1 =  count_coverage($sam_d, $chr1, $pos1);
		my $bg_reads2 =  count_coverage($sam_d, $chr1, $pos2);
		next if($bg_reads1*$bg_reads2 == 0 || ($qc_flanking/$bg_reads1 < 0.01 && $qc_flanking/$bg_reads2 < 0.01));
	}

	my $first_bp = {
		clip => RIGHT_CLIP,
		reads_num => $qc_perfect_reads,
		gene => $gene1,
		tname => $chr1,
		tpos => $pos1,
		qstrand => $strand1 eq "+" ? 1 : -1,
	};
	$first_bp->{ort} = $first_bp->{clip}*$first_bp->{qstrand};

	my $second_bp = {
		clip => LEFT_CLIP,
		reads_num => $qc_perfect_reads,
		gene => $gene2,
		tname => $chr2,
		tpos => $pos2,
		qstrand => $strand2 eq "+" ? 1 : -1,
	};
	my $same_gene = same_gene($first_bp->{gene}, $second_bp->{gene});
	$second_bp->{ort} = $second_bp->{clip}*$second_bp->{qstrand};
	my $type = get_type($first_bp, $second_bp, $same_gene);
	next if($type eq 'Internal_splicing');

	my $tmp_SV = {
		first_bp => $first_bp,
		second_bp => $second_bp,
		type => "Junction",
		};
	if($tmp_SV && ! is_dup_raw_SV(\@raw_SVs, $tmp_SV)){
		push @raw_SVs, $tmp_SV;
	}
  }
  close($JUNC);
}

my $New_excluded_file = "$out_dir/excluded.new.txt";
if ($internal){
	$New_excluded_file = "$out_dir/excluded.new.internal.txt";
}
open(my $NEXC, ">$New_excluded_file");
foreach my $g (sort { $gene_recurrance{$b} <=> $gene_recurrance{$a} } keys %gene_recurrance) {
	last if($gene_recurrance{$g} < $max_num_hits*10);
	next if($g eq "NA" || $g eq "IGH" || $g eq "TCRA" || $g eq "TCRB" || $g eq "TARP");
	next if($cr_hash{$g});
#	$excluded{$g} = 1;#Tian to rescue IGH-CRLF2 when CRLF2 has many partner genes, will not add CRLF2 to excluded
	print STDERR "Adding $g to excluded with recurrance: ".$gene_recurrance{$g}."\n";
	print $NEXC join("\t",$g, $gene_recurrance{$g}),"\n";
}
close($NEXC);

`mkdir  -p $out_dir/tmp_anno`;
if ($?){
	my $err = $!;
	print STDERR "Error creating annotation directory: $err\n";
	exit 4;
}
my $annotation_dir = tempdir(DIR => "$out_dir/tmp_anno");
`mkdir -p $annotation_dir`;
if ($?){
	my $err = $!;
	print STDERR "Error creating temp directory for annotation: $err\n";
	exit 5;
}
print STDERR "Annotation Dir: $annotation_dir\n" if($debug);

### write SVs to file

my $raw_file = "$out_dir/raw.fusion.txt";
if($internal) {
	$raw_file = "$out_dir/raw.internal.txt";
}

open(hFo, ">$raw_file");
foreach my $raw_SV (@raw_SVs){
	my ($first_bp, $second_bp) = ($raw_SV->{first_bp}, $raw_SV->{second_bp});
	my @row = ($first_bp->{gene}, $second_bp->{gene},
				$first_bp->{reads_num}, $second_bp->{reads_num},
				$first_bp->{tpos}, $second_bp->{tpos},
				$first_bp->{ort}, $first_bp->{tname}, $first_bp->{qstart}, $first_bp->{qend},
				$first_bp->{qstrand}, $first_bp->{matches}, $first_bp->{percent}, $first_bp->{repeat}, 
				$first_bp->{clip},
				$second_bp->{ort}, $second_bp->{tname}, $second_bp->{qstart}, $second_bp->{qend},
				$second_bp->{qstrand}, $second_bp->{matches}, $second_bp->{percent}, $second_bp->{repeat}, 
				$second_bp->{clip},
				$raw_SV->{junc_seq},$raw_SV->{type}
				);
	my $out_string = join("\t", map { defined ? $_ : '' } @row);

	print hFo $out_string, "\n";
}
close(hFo);




### start of annotate.pl tail

sub is_good_ITD {
	my($bp1, $bp2) = @_;
	my @genes = split(/,|\|/, $bp1->{gene});

	foreach my $g1 (@genes) {
		return 1 if(%known_ITDs && exists($known_ITDs{$g1}) &&
                	$bp1->{tpos} > $known_ITDs{$g1}[0] &&
                	$bp1->{tpos} < $known_ITDs{$g1}[1] &&
                	$bp2->{tpos} > $known_ITDs{$g1}[0] &&
                	$bp2->{tpos} < $known_ITDs{$g1}[1]);
	}

	return 0;
}

sub is_dup_raw_SV {
	my($r_SVs, $sv) = @_;
	foreach my $s (@{$r_SVs}) {
		return 1
		if( abs($s->{first_bp}->{tpos} - $sv->{first_bp}->{tpos}) < 10 &&
		    abs($s->{second_bp}->{tpos} - $sv->{second_bp}->{tpos}) < 10 &&
			$s->{first_bp}->{tname} eq $sv->{first_bp}->{tname} &&
			$s->{second_bp}->{tname} eq $sv->{second_bp}->{tname}
		);

		if( abs($s->{first_bp}->{tpos} - $sv->{second_bp}->{tpos}) < 10 &&
		    abs($s->{second_bp}->{tpos} - $sv->{first_bp}->{tpos}) < 10 &&
			$s->{first_bp}->{tname} eq $sv->{second_bp}->{tname} &&
			$s->{second_bp}->{tname} eq $sv->{first_bp}->{tname}
		){
			$s->{second_bp} = $sv->{first_bp};
			return 1;
		}
	}
	return 0;
}

sub count_coverage {
	my ($sam, $chr, $pos) = @_;
	my $seg = $sam->segment(-seq_id => $chr, -start => $pos, -end => $pos);
	return 0 unless $seg;
	my $n = 0;
	my $itr = $seg->features(-iterator => 1);
	while( my $a = $itr->next_seq) {
		next unless($a->start && $a->end); #why unmapped reads here?
		$n++;
	}
	return $n;
}


sub in_complex_region{
	my ($chr, $pos) = @_;
	my $full_chr = ($chr =~ m/chr/) ? $chr : "chr$chr";
	foreach my $cr (@complex_regions){
		return $cr->{name} if($cr->{chr} eq $full_chr && $pos > $cr->{start} && $pos < $cr->{end});
	}
	return 0;
}

sub is_bad_chrom{
	my $chr = shift;
	foreach my $bad_chr (@excluded_chroms){
		return 1 if($chr =~ /$bad_chr/i);
	}
	return 0;
}

sub low_complexity{
	my $sequence = shift;
	my $max_run_nt = 20;
	my $mask_seq = $sequence;
	$mask_seq =~ s/((.+)\2{9,})/'N' x length $1/eg;
	return 1 if $mask_seq =~ /(N{$max_run_nt,})/;

	my $len= length($mask_seq);
	my $seg_len = 25;
	for (my $i=0; $i<$len-$seg_len; $i++){
		my $sub_seq = substr $mask_seq, $i, $seg_len;
		my $n = @{[$sub_seq =~ /(N)/g]};
		return 1 if($n>20);
	}
	return 0;
}

sub get_type {

	my ($first_bp, $second_bp, $same_gene, $transcript_strand) = @_;
	my $debug = 0;#Tian
	print STDERR "=== get_type ===", join("\t", $first_bp->{gene}, $second_bp->{gene}, $second_bp->{tname}, $first_bp->{tname}), "\n" if($debug);
	return "CTX" if($second_bp->{tname} ne $first_bp->{tname});

	my $clip1 = $first_bp->{qstrand}*$first_bp->{ort};
	if($same_gene) { # internal events
		return "Internal_inv" if($second_bp->{qstrand}*$first_bp->{qstrand} < 0);
		print STDERR "clip1: ", $clip1, "\ttpos1: ", $first_bp->{tpos}, "\ttpos2:", $second_bp->{tpos}, "\n" if($debug);
		return 'Internal_splicing' if($clip1*($first_bp->{tpos} -  $second_bp->{tpos} - 10) < 0);
		return 'Internal_dup';
	}
	else {
		return 'ITX' if($second_bp->{qstrand}*$first_bp->{qstrand} < 0);
		if($clip1 * ($first_bp->{tpos} - $second_bp->{tpos} - 10) < 0){
			return "read_through" if(abs($first_bp->{tpos} - $second_bp->{tpos}) < 100000);
			return 'DEL';
		}
		return 'INS';
	}
	return 'undef';
}

sub same_gene{

	my ($gene1, $gene2) = @_;
	my @g1_names = split(/,|\|/,$gene1);
	my @g2_names = split(/,|\|/,$gene2);
	my $same_gene=0;
	foreach my $g1 (@g1_names){
		foreach my $g2 (@g2_names){
			return 0 if($g1 eq "NA" && $g2 eq "NA");
			return 1 if($g1 eq $g2);
		}
	}
	return 0;
}

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
