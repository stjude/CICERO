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

use constant POSITION_RESCUE_MAX_EXON_DISTANCE => 15;
# maximum distance from rescue position to an exon/splice boundary
use constant POSITION_RESCUE_MAX_INTRON_DISTANCE => 100000;
# maximum adjusted position distance, should be approximate max intron size

my $debug = 0;

my ($blat_server, $blat_port, $dir_2bit);
my $blat_client_options = ' -out=psl -nohead > /dev/null 2>&1';
my $cap3_options = " -o 25 -z 2 -h 60 > /dev/null";

my $rmdup = 0;
my $paired = 1;

# input/output
my ($genome, $ref_genome, $gene_model_file);
my ($out_dir, $known_itd_file, $known_fusion_file);
my ($input_bam, $raw_file, $read_len);
my ($all_output, $internal, $DNA) = (0, 0, 0);
my ($min_hit_len, $max_num_hits, $min_fusion_distance);
$min_hit_len = 25;
my $sc_shift = 10;
my ($excluded_gene_file, $complex_region_file, $gold_gene_file);
my ($help, $man, $version, $usage );

if(@ARGV == 0){
	#TODO: get the correct usage string
	print STDERR "Usage: $0 -g <genome> -i <bam_file> -r <raw_file> -o <out_dir> -f <gene_info>\n";
	exit 1;
}

my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,
	'r|raw_file=s'	=> \$raw_file,
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'  => \$genome,
	'ref_genome=s'  => \$ref_genome,
	'genemodel=s'		=> \$gene_model_file,
	'excluded-genes=s'	=> \$excluded_gene_file,
	'excluded-fusions=s'	=> \$excluded_fusion_file,

	'blatserver' =>	\$blat_server,
	'blatport=s'		=> \$blat_port,
	'min_hit_len=i'		=> \$min_hit_len,
	'max_num_hits=i'	=> \$max_num_hits,
	'c|cluster=i'   => \$sc_shift,
	'paired!'		=> \$paired,
	'rmdup!'		=> \$rmdup,
	'l|read_len=i'	=> \$read_len,
	'internal!' => \$internal,
	'all!' => \$all_output,
	'DNA!' => \$DNA,
	'known_itd_file=s'	=> \$known_itd_file,
	'known_fusion_file=s'	=> \$known_fusion_file,
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
$blat_server = $conf->{BLAT_HOST} unless($blat_server);
$blat_port = $conf->{BLAT_PORT} unless($blat_port);
$dir_2bit = '/';
$gene_model_file = $conf->{'REFSEQ_REFFLAT'} unless($gene_model_file);
print STDERR "gene_model_file: $gene_model_file\n";
$excluded_gene_file = $conf->{EXCLUDED_GENES} unless($excluded_gene_file);
$known_itd_file = $conf->{KNOWN_ITD_FILE} unless($known_itd_file);
print STDERR "KNOWN_ITD_FILE: ", $known_itd_file, "\n";
$known_fusion_file = $conf->{KNOWN_FUSIONS} unless($known_fusion_file);
$gold_gene_file = $conf->{CICERO_GOLD_GENE_LIST_FILE};
print STDERR "CICERO_GOLD_GENE_LIST_FILE:", $gold_gene_file, "\n";
$complex_region_file = $conf->{COMPLEX_REGIONS} unless($complex_region_file);

# Load Cicero-specific configuration settings
$conf = &TdtConfig::readConfig('app', 'cicero');
$min_hit_len = $conf->{MIN_HIT_LEN} unless($min_hit_len);
$max_num_hits = $conf->{MAX_NUM_HITS} unless($max_num_hits);
$min_fusion_distance = $conf->{MIN_FUSION_DIST} unless($min_fusion_distance);

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

# Intialize CAP3 and BLAT utilities
# These variables will be global
my $assembler = Assembler->new(
	-PRG => "cap3",
	-OPTIONS => $cap3_options
);

my $mapper = Mapper->new(
	-PRG => join(' ', ("gfClient", $blat_server, $blat_port)),
	-OPTIONS => $blat_client_options,
	-BIT2_DIR => $dir_2bit,
	-MIN_HIT_LEN => $min_hit_len,
	-MAX_NUM_HITS => $max_num_hits,
	-MIN_FS_DIST => $min_fusion_distance
);

# Read the gene model from a refflat file
my $gm_format = "REFFLAT";
my $gm = GeneModel->new if($gene_model_file);
$gm->from_file($gene_model_file, $gm_format);

croak "Specify a genome: $ref_genome" unless (-f $ref_genome);

# Load the BAM file
my $sam_d = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $ref_genome);
my @seq_ids = $sam_d->seq_ids;
my $validator = CiceroSCValidator->new();
$validator->remove_validator('strand_validator') if(!$paired);

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

my %known_ITDs = ();
open(my $ITD_F, $known_itd_file) or die "Cannot open $known_itd_file: $!";
while(<$ITD_F>){
	chomp;
	my ($gene, $chr, $start, $end) = split(/\t/,$_);
	$known_ITDs{$gene} = [$start, $end];
}

# Load gene pairs that are enhancer activated
# involving T-cell receptors and immunoglobulin loci
my %enhancer_activated_genes = ();
my %known_fusion_partners = ();
if($known_fusion_file && -e $known_fusion_file){
   open(my $KFF, $known_fusion_file);
   while(<$KFF>){
	chomp;
	my ($gene1, $gene2) = split(/\t/,$_);
	next if exists($enhancer_activated_genes{$gene1});
	next if exists($enhancer_activated_genes{$gene2});
	$enhancer_activated_genes{$gene2} = $gene1 if($gene1 =~ m/^IG.$/ || $gene1 =~ m/^TR.$/);
	$enhancer_activated_genes{$gene1} = $gene2 if($gene2 =~ m/^IG.$/ || $gene2 =~ m/^TR.$/);
	$known_fusion_partners{$gene1}{$gene2} = 1;
	$known_fusion_partners{$gene2}{$gene1} = 1;
   }
   close($KFF);
}

my @raw_SVs = ();

### end of annotate.pl head


###Load split raw SV file

#my $raw_file = "$out_dir/raw.fusion.txt";
#if($internal) {
#	$raw_file = "$out_dir/raw.internal.txt";
#}

open my $RSV, '<',  "$raw_file";
chomp(my @lines = <$RSV>);
close $RSV;

### Convert lines to SV objects

foreach my $line (@lines){

	my @fields = split("\t", $line);
	my $first_bp = {
		reads_num => $fields[2],
		gene => $fields[0],
		tpos => $fields[4],
		ort => $fields[6],
		tname => $fields[7],
		qstart => $fields[8],
		qend => $fields[9],
		qstrand => $fields[10],
		matches => $fields[11],
		percent => $fields[12],
		repeat => $fields[13],
		clip => $fields[14]
	};

	my $second_bp = {
		reads_num => $fields[3],
		gene => $fields[1],
		tpos => $fields[5],
		ort => $fields[15],
		tname => $fields[16],
		qstart => $fields[17],
		qend => $fields[18],
		qstrand => $fields[19],
		matches => $fields[20],
		percent => $fields[21],
		repeat => $fields[22],
		clip => $fields[23]
	};

	my $sv = {
		first_bp => $first_bp,
		second_bp => $second_bp,
		};

	$sv->{junc_seq} = $fields[24];
	$sv->{type} = $fields[25];

	push @raw_SVs, $sv;
}


my $tmp_name = $raw_file;
$tmp_name =~ s/.txt//;

`mkdir  -p $tmp_name`;
if ($?){
	my $err = $!;
	print STDERR "Error creating annotation directory: $err\n";
	exit 4;
}

my $annotation_dir = tempdir(DIR => "$tmp_name");

`mkdir -p $annotation_dir`;
if ($?){
	my $err = $!;
	print STDERR "Error creating temp directory for annotation: $err\n";
	exit 5;
}


print STDERR "Processing raw SVs\n";
my @annotated_SVs;
foreach my $sv (@raw_SVs){
#print STDERR "$sv->{type}\n";
	next if($sv->{type} eq "Internal_splicing");

	my ($first_bp, $second_bp, $contigSeq) = ($sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq});
	my ($gene1, $gene2) = ($first_bp->{gene}, $second_bp->{gene});
	my @genes1 = split(/,|\|/, $gene1);
	my @genes2 = split(/,|\|/, $gene2);
	my $bad_gene = 0;
	print STDERR "xxx\n" if(abs($sv->{second_bp}->{tpos} - 170818803)<10 || abs($sv->{first_bp}->{tpos} - 170818803)<10);
	foreach my $g1 (@genes1) {
		if(exists($excluded{$g1})) {$bad_gene = 1; last;}
	}
	next if($bad_gene);
	foreach my $g2 (@genes2){
		if(exists($excluded{$g2})) {$bad_gene = 1; last;}
	}
	next if($bad_gene);
	print STDERR "next if($contigSeq &&  > $max_num_hits)\n" if(abs($sv->{second_bp}->{tpos} - 170818803)<10 || abs($sv->{first_bp}->{tpos} - 170818803)<10);

	my $bp1_site = join("_", $first_bp->{tname}, $first_bp->{tpos}, $first_bp->{clip});
	my $bp2_site = join("_", $second_bp->{tname}, $second_bp->{tpos}, $second_bp->{clip});
	$bp1_site = normalizeChromosomeName($seq_ids[0], $bp1_site);
	$bp2_site = normalizeChromosomeName($seq_ids[0], $bp2_site);


	my $start_run = time();
	print STDERR "\nstart to quantify the fusion... ", join(" ", $sv->{first_bp}->{tname}, $sv->{first_bp}->{tpos}, $sv->{second_bp}->{tname}, $sv->{second_bp}->{tpos}), "\n" if(abs($sv->{second_bp}->{tpos} - 170818803)<10 || abs($sv->{first_bp}->{tpos} - 170818803)<10);
	print STDERR "\nstart to quantify the fusion... ", join(" ", $sv->{first_bp}->{tname}, $sv->{first_bp}->{tpos}, $sv->{second_bp}->{tname}, $sv->{second_bp}->{tpos}), "\n" if($debug);
	my @quantified_SVs = quantification(-SAM => $sam_d,
		 	-GeneModel => $gm,
		 	-VALIDATOR => $validator,
		 	-PAIRED => $paired,
		 	-SV => $sv,
		 	-ANNO_DIR => $annotation_dir
		   );
 	my $end_run = time();
	my $run_time = $end_run - $start_run;
	print STDERR "Quantification run time: $run_time\n" if($debug);

	foreach my $quantified_SV (@quantified_SVs){
		my $annotated_SV = annotate($gm, $quantified_SV) if($quantified_SV);
		print STDERR "annotated_SV: $annotated_SV\n" if($debug);
		next unless($annotated_SV);
		my ($first_bp, $second_bp, $type) = ($annotated_SV->{first_bp}, $annotated_SV->{second_bp}, $annotated_SV->{type});
		if(!$all_output && $type =~ m/Internal/){
			if($type eq 'Internal_dup'){
				next unless (is_good_ITD($first_bp, $second_bp) || exist_multiplename_checking(\%gold_genes, $first_bp->{gene}));
			}
			elsif(!$DNA){
				next if($type eq 'Internal_splicing');
			        next unless( $first_bp->{feature} =~ m/coding/ ||
				    $second_bp->{feature} =~ m/coding/);
			}
		}
		push @annotated_SVs, $annotated_SV;
	}
}

### write SVs to file 

my $annotated_file = $raw_file;
$annotated_file =~ s/raw/quantified/;

open(hFo, ">$annotated_file");
foreach my $annotated_SV (@annotated_SVs){
	my ($first_bp, $second_bp) = ($annotated_SV->{first_bp}, $annotated_SV->{second_bp});
	my @row = ($first_bp->{gene}, $second_bp->{gene},
				$first_bp->{reads_num}, $second_bp->{reads_num},
				$first_bp->{tpos}, $second_bp->{tpos},
				$first_bp->{ort}, $first_bp->{tname}, $first_bp->{qpos}, $first_bp->{qstart}, $first_bp->{qend},
				$first_bp->{qstrand}, $first_bp->{ts_strand}, $first_bp->{matches}, $first_bp->{percent}, $first_bp->{repeat}, 
				$first_bp->{clip}, $first_bp->{area}, $first_bp->{feature}, $first_bp->{annotate_score},
				$second_bp->{ort}, $second_bp->{tname}, $second_bp->{qpos}, $second_bp->{qstart}, $second_bp->{qend},
				$second_bp->{qstrand}, $second_bp->{ts_strand}, $second_bp->{matches}, $second_bp->{percent}, $second_bp->{repeat}, 
				$second_bp->{clip}, $second_bp->{area}, $second_bp->{feature}, $second_bp->{annotate_score},
				$annotated_SV->{junc_seq},$annotated_SV->{type},$annotated_SV->{ort}
				);
	my $out_string = join("\t", map { defined ? $_ : '' } @row);

	print hFo $out_string, "\n";
}
close(hFo);

rmtree(["$annotation_dir"]);

print "annotated_SVs: ", scalar @annotated_SVs, "\n";

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

sub in_complex_region{
	my ($chr, $pos) = @_;
	my $full_chr = ($chr =~ m/chr/) ? $chr : "chr$chr";
	foreach my $cr (@complex_regions){
		return $cr->{name} if($cr->{chr} eq $full_chr && $pos > $cr->{start} && $pos < $cr->{end});
	}
	return 0;
}

sub annotate {
	my %args = @_;
	my ($gm, $SV) = @_;
	my $debug = 0;
	print "\n=== annotate SV ===\n" if($debug);
	my ($first_bp, $second_bp, $qseq) = ($SV->{first_bp}, $SV->{second_bp}, $SV->{junc_seq});
	my($annotated_first_bp, $annotated_second_bp) = ($first_bp, $second_bp);
	my($crA, $crB) = (in_complex_region($first_bp->{tname}, $first_bp->{tpos}), in_complex_region($second_bp->{tname}, $second_bp->{tpos}));
	if($crA){
		$annotated_first_bp->{gene} = $crA;
		$annotated_first_bp->{feature} = '5utr';
		$annotated_first_bp->{annotate_score} = 0.5;
		$annotated_second_bp = annotate_enhancer_gene_bp($second_bp) unless($crB);
	}

	if($crB){
		$annotated_second_bp->{gene} = $crB;
		$annotated_second_bp->{feature} = '5utr';
		$annotated_second_bp->{annotate_score} = 0.5;
		$annotated_first_bp = annotate_enhancer_gene_bp($first_bp) unless($crA);
	}

	print STDERR "\nannotation... --- junc_seq ", $qseq, "\n" if($debug);
	print STDERR "clip info: ", $first_bp->{clip}, "\t", $second_bp->{clip}, "\n" if($debug);
	print STDERR "\nstart bp1 annotation .....\n" if($debug);
	$annotated_first_bp = annotate_bp($first_bp) unless($annotated_first_bp->{feature});
	print STDERR "bp1 annotation: ", join("\t", $annotated_first_bp->{feature}, $annotated_first_bp->{annotate_score},
		$annotated_first_bp->{gene}), "\n" if($debug);

	print STDERR "\nstart bp2 annotation .....\n" if($debug);
	$annotated_second_bp = annotate_bp($second_bp) unless($annotated_second_bp->{feature});
	print STDERR "bp2 annotation: ", join("\t", $annotated_second_bp->{feature}, $annotated_second_bp->{annotate_score},
		$annotated_second_bp->{gene}), "\n" if($debug);

	my $qseq_ort = sign($annotated_first_bp->{annotate_score} + $annotated_second_bp->{annotate_score}, 'd');
	if($qseq_ort < 0){
		$qseq = rev_comp($qseq);
		$annotated_first_bp->{qstrand} = -1*$annotated_first_bp->{qstrand};
		$annotated_first_bp->{ort} = -1*$annotated_first_bp->{ort};
		$annotated_second_bp->{qstrand} = -1*$annotated_second_bp->{qstrand};
		$annotated_second_bp->{ort} = -1*$annotated_second_bp->{ort};
		$annotated_first_bp->{qpos} = length($qseq) - $annotated_first_bp->{qpos} + 1;
		$annotated_second_bp->{qpos} = length($qseq) - $annotated_second_bp->{qpos} + 1;
	}
	print STDERR "\nfirst bp ort -- ", $first_bp->{ort}, "\n" if($debug);

	my $annotated_SV;
	print STDERR "\njunc_seq ", $qseq, "\nfinished annotation\n\n" if($debug);
	if($first_bp->{ort} < 0){

		$annotated_SV = {
			junc_seq => $qseq,
			first_bp => $annotated_second_bp,
			second_bp => $annotated_first_bp,
			ort => '>',
			};
	}
	else {
		$annotated_SV = {
			junc_seq => $qseq,
			first_bp => $annotated_first_bp,
			second_bp => $annotated_second_bp,
			ort => '>',
			};
	}
	my $same_gene = same_gene($annotated_first_bp->{gene}, $annotated_second_bp->{gene});

	my $type = ($qseq_ort ==0) ? get_type($annotated_first_bp, $annotated_second_bp, $same_gene) : get_type($annotated_first_bp, $annotated_second_bp, $same_gene, $annotated_first_bp->{qstrand});
	return if($type eq "Internal_splicing");
	print STDERR "type = $type\n" if($debug);

	if($qseq_ort == 0){
		$annotated_SV->{ort} = '?';
	}
	$annotated_SV->{type} = $type;
	return $annotated_SV;
} #end of annotate

sub sign{

	my $a = shift;
	my $b = shift;

	if($b eq 'c'){
		return '+' if($a > 0);	
		return '-' if($a < 0);	
		return '=' if($a == 0);	
	}

	if($b eq 'd'){
		return 1 if($a > 0);	
		return -1 if($a < 0);	
		return 0 if($a == 0);	
	}
}

sub annotate_enhancer_gene_bp{

	my $bp = shift;
	my $debug = 0;
	my $chr = ($bp->{tname} =~ m/chr/) ? $bp->{tname} : "chr".$bp->{tname};
	my $tpos = $bp->{tpos};
	my $strand = ($bp->{qstrand} > 0) ? '+' : '-';
	print STDERR "\n=== annotating bp at ", $chr, ":", $tpos, "\t$strand ===\n" if($debug);
	my $qseq_ort = ($bp->{ort} > 0) ? '+' : '-';
	my $rev_strand = ($bp->{qstrand} > 0) ? '-' : '+';

	my $extend_size = 1000000;
	
		my ($start, $end) = ($tpos - $extend_size, $tpos + $extend_size);
		my $gm_tree = $gm->sub_model($chr, $strand);
		return if(!defined($gm_tree));
		my @tmp = $gm_tree->intersect([$start, $end]);
		foreach my $g (@tmp){
			$g=$g->val;
			next unless(exists($enhancer_activated_genes{$g->name}));
			my ($tmp_feature, $tmp_score);
			print STDERR "gene at $strand is: ", join("\t", $g->name, $g->start, $g->end),"\n" if($debug);
			my $check_point = ($qseq_ort eq $strand) ? ($tpos - 10) : ($tpos + 10);
			$tmp_feature = $g->get_feature($chr, $check_point, $strand);
			print STDERR "$tmp_feature = g->get_feature($chr, $check_point, $strand)\n" if($debug);
			$tmp_score = 1 if($tmp_feature eq 'coding');
			$tmp_score = 0.8 if($tmp_feature =~ m/utr/);
			$tmp_score = 0.5 if($tmp_feature eq 'intron');
			$tmp_score = 0.1 if($tmp_feature eq 'intergenic');
	
			$bp->{annotate_score} = $tmp_score;
			$bp->{feature} = $tmp_feature;
			$bp->{ts_strand} = $bp->{qstrand};
			$bp->{gene} = $g->name;
			return $bp if($bp->{annotate_score} == 1);
		}

		print STDERR "gm_tree = gm->sub_model($chr, $rev_strand)\n" if($debug);
		$gm_tree = $gm->sub_model($chr, $rev_strand);
		@tmp = $gm_tree->intersect([$start, $end]);
		# return gene orientation, qseq_ort, annotation ...
		foreach my $g (@tmp){
			$g=$g->val;
			next unless(exists($enhancer_activated_genes{$g->name}));
			my ($tmp_feature, $tmp_score);
			print STDERR "gene at $rev_strand is: ", join("\t", $g->name, $g->start, $g->end),"\n" if($debug);
			my $check_point = ($qseq_ort eq $strand) ? ($tpos - 10) : ($tpos + 10);
			$tmp_feature = $g->get_feature($chr, $check_point, $rev_strand);
			print STDERR "$tmp_feature = g->get_feature($chr, $check_point, $rev_strand)\n" if($debug);
			$tmp_score = -1 if($tmp_feature eq 'coding');
			$tmp_score = -0.8 if($tmp_feature =~ m/utr/);
			$tmp_score = -0.5 if($tmp_feature eq 'intron');
			$tmp_score = -0.1 if($tmp_feature eq 'intergenic');
			$bp->{annotate_score} = $tmp_score;
			$bp->{feature} = $tmp_feature;
			$bp->{ts_strand} = -1*$bp->{qstrand};
			$bp->{gene} = $g->name;
			return $bp if(abs($bp->{annotate_score}) == 1);
		}
	return $bp;
}

sub annotate_bp{

	my $bp = shift;
	my $debug = 0;
	my $chr = ($bp->{tname} =~ m/chr/) ? $bp->{tname} : "chr".$bp->{tname};
	my $tpos = $bp->{tpos};
	print STDERR "\n=== annotating bp at ", $chr, ":", $tpos, " ===\n" if($debug);
	my $strand = ($bp->{qstrand} > 0) ? '+' : '-';
	my $qseq_ort = ($bp->{ort} > 0) ? '+' : '-';
	my $rev_strand = ($bp->{qstrand} > 0) ? '-' : '+';
	$bp->{annotate_score} = 0;
	$bp->{feature} = 'intergenic';
	$bp->{ts_strand} = 0;
	$bp->{gene} = 'NA';
	my $dist = 40000;

	#priority order: same_direction exon > diff_direction exon > same_direction utr > diff_direction utr
	#   > same_direction intron > diff_direction intron > min_distance intergenic (for all same/diff direction intergenic genes)
	foreach my $extend_size (10, 5000, 10000, 40000){
	
		my ($start, $end) = ($tpos - $extend_size, $tpos + $extend_size);
		my $gm_tree = $gm->sub_model($chr, $strand);
		return if(!defined($gm_tree));
		my @tmp = $gm_tree->intersect([$start, $end]);
		foreach my $g (@tmp){
			$g=$g->val;
			my ($tmp_feature, $tmp_score);
			print STDERR "gene at $strand is: ", join("\t", $g->name, $g->start, $g->end),"\n" if($debug);
			my $check_point = ($qseq_ort eq $strand) ? ($tpos - 10) : ($tpos + 10);
			$tmp_feature = $g->get_feature($chr, $check_point, $strand);
			print STDERR "$tmp_feature = g->get_feature($chr, $check_point, $strand)\n" if($debug);
			$tmp_score = 1 if($tmp_feature eq 'coding');
			$tmp_score = 0.8 if($tmp_feature =~ m/utr/);
			$tmp_score = 0.5 if($tmp_feature eq 'intron');
			$tmp_score = 0.1 if($tmp_feature eq 'intergenic');
			my $tmp_dist = (abs($g->start - $tpos) < abs($g->end - $tpos)) ? abs($g->start - $tpos) : abs($g->end - $tpos);
	
			return $bp if($bp->{annotate_score} == 1);

                        if($tmp_feature eq 'intergenic')
                        {
                                if($tmp_dist < $dist && abs($bp->{annotate_score}) < 0.11){#Tian if saved feature is +/- intergenic, update the annotation if dist is smaller
                                        $bp->{annotate_score} = $tmp_score;
                                        $bp->{feature} = $tmp_feature;
                                        $bp->{ts_strand} = $bp->{qstrand};
                                        $bp->{gene} = $g->name;
                                        $dist = $tmp_dist;
                                        #print STDERR "update1 ", $tmp_score, " ", $tmp_feature, " ", $bp->{ts_strand}, " ", $bp->{gene}, " ", $dist, "\n";
                                }
                        }
                        else{
                                if($tmp_score > abs($bp->{annotate_score})){
                                        $bp->{annotate_score} = $tmp_score;
                                        $bp->{feature} = $tmp_feature;
                                        $bp->{ts_strand} = $bp->{qstrand};
                                        $bp->{gene} = $g->name;
                                        #print STDERR "update1 ", $tmp_score, " ", $tmp_feature, " ", $bp->{ts_strand}, " ", $bp->{gene}, " ", $dist, "\n";
                                }
                        }
		}

		print STDERR "gm_tree = gm->sub_model($chr, $rev_strand)\n" if($debug);
		$gm_tree = $gm->sub_model($chr, $rev_strand);
		@tmp = $gm_tree->intersect([$start, $end]);
		# return gene orientation, qseq_ort, annotation ...
		foreach my $g (@tmp){
			$g=$g->val;
			my ($tmp_feature, $tmp_score);
			print STDERR "gene at $rev_strand is: ", join("\t", $g->name, $g->start, $g->end),"\n" if($debug);
			my $check_point = ($qseq_ort eq $strand) ? ($tpos - 10) : ($tpos + 10);
			$tmp_feature = $g->get_feature($chr, $check_point, $rev_strand);
			print STDERR "$tmp_feature = g->get_feature($chr, $check_point, $rev_strand)\n" if($debug);
			$tmp_score = -1 if($tmp_feature eq 'coding');
			$tmp_score = -0.8 if($tmp_feature =~ m/utr/);
			$tmp_score = -0.5 if($tmp_feature eq 'intron');
			$tmp_score = -0.1 if($tmp_feature eq 'intergenic');
			my $tmp_dist = (abs($g->start - $tpos) < abs($g->end - $tpos)) ? abs($g->start - $tpos) : abs($g->end - $tpos);
			
			return $bp if($bp->{annotate_score} == -1);

                        if($tmp_feature eq 'intergenic')
                        {
                                if($tmp_dist < $dist && abs($bp->{annotate_score}) < 0.11){#Tian if saved feature is +/- intergenic, update the annotation if dist is smaller
                                        $bp->{annotate_score} = $tmp_score;
                                        $bp->{feature} = $tmp_feature;
                                        $bp->{ts_strand} = -1*$bp->{qstrand};
                                        $bp->{gene} = $g->name;
                                        $dist = $tmp_dist;
                                        #print STDERR "update2 ", $tmp_score, " ", $tmp_feature, " ", $bp->{ts_strand}, " ", $bp->{gene}, " ", $dist, "\n";
                                }
                        }
                        else{
                                if(abs($tmp_score) > abs($bp->{annotate_score})){
                                        $bp->{annotate_score} = $tmp_score;
                                        $bp->{feature} = $tmp_feature;
                                        $bp->{ts_strand} = -1*$bp->{qstrand};
                                        $bp->{gene} = $g->name;
                                        #print STDERR "update2 ", $tmp_score, " ", $tmp_feature, " ", $bp->{ts_strand}, " ", $bp->{gene}, " ", $dist, "\n";
                                }
                        }
		}
	}
	return $bp;
}

sub quantification {
	my $debug = 0;
	my %args = @_;
	my ($gm, $sam, $validator, $paired, $SV, $anno_dir) =
		($args{-GeneModel}, $args{-SAM}, $args{-VALIDATOR}, $args{-PAIRED}, $args{-SV}, $args{-ANNO_DIR});
	my ($bp1, $bp2) = ($SV->{first_bp}, $SV->{second_bp});
	my ($chr1, $pos1, $start1, $end1) = ($bp1->{tname}, $bp1->{tpos}, $bp1->{ort}, $bp1->{tstart}, $bp1->{tend});
	my ($chr2, $pos2, $start2, $end2) = ($bp2->{tname}, $bp2->{tpos}, $bp2->{ort}, $bp2->{tstart}, $bp2->{tend});
	$debug = 0 if(abs($pos1 - 170818803)<10 || abs($pos2 - 170818803)<10);
	print STDERR "xxx\n" if(abs($pos1 - 170818803)<10 || abs($pos2 - 170818803)<10);
	my $fixSC1 = $bp1->{reads_num} < 10 ? 1 : 0;
	my $fixSC2 = $bp2->{reads_num} < 10 ? 1 : 0;

	# right clip or left clip?
	my $clip1;
	if($bp1->{ort} != 1 && $bp1->{ort} != -1){
		print STDERR "bp1->ort ", $bp1->{ort}, " error!\n";
		exit 6;
	}
	$clip1 = $bp1->{ort}*$bp1->{qstrand};

	# right clip or left clip?
	my $clip2;
	if($bp2->{ort} != 1 && $bp2->{ort} != -1){
		print STDERR "bp2->ort ", $bp2->{ort}, " error!\n";
		exit 7;
	}
	$clip2 = $bp2->{ort}*$bp2->{qstrand};

	my $gap_size = 0;
	if($chr1 eq $chr2){
		$gap_size = abs($pos2-$pos1) if(($pos1 < $pos2 && $clip1 == RIGHT_CLIP && $clip2 == LEFT_CLIP) ||
		  ($pos1 > $pos2 && $clip2 == RIGHT_CLIP && $clip1 == LEFT_CLIP));
	}

	my $rmdup=1;
	my $clip1x = $clip1 + 1;
	my $fa_file1 = "$anno_dir/".join(".", $chr1,$pos1, ($clip1+1), "fa");
	my $output_mate = 1;
	$output_mate = 0 if($SV->{type} eq "Internal_dup");	
	prepare_reads_file(
			-OUT => $fa_file1,
		        -SAM => $sam,
			-CHR =>$chr1,
			-POS => $pos1,
			-CLIP => $clip1,
			-VALIDATOR => $validator,
			-PAIRED => $paired,
			-RMDUP => $rmdup,
			-MIN_SC => 1,
			-SC_SHIFT => $sc_shift,
			-MIN_SC_LEN => 3,
			-GAP_SIZE => $gap_size,
			-FIXSC => $fixSC1,
			-UNMAPPED_CUTOFF => 1000,
			-MATE => $output_mate
	        	) unless(-s $fa_file1);
	print STDERR "fa_file1: *$fa_file1*\n" if($debug);

	my $fa_file2 = "$anno_dir/".join(".", $chr2, $pos2, ($clip2+1), "fa");
	prepare_reads_file(
			-OUT => $fa_file2,
		        -SAM => $sam,
			-CHR =>$chr2,
			-POS => $pos2,
			-CLIP => $clip2,
			-VALIDATOR => $validator,
			-PAIRED => $paired,
			-RMDUP => $rmdup,
			-MIN_SC => 1,
			-SC_SHIFT => $sc_shift,
			-MIN_SC_LEN => 3,
			-GAP_SIZE => $gap_size,
			-FIXSC => $fixSC2,
			-UNMAPPED_CUTOFF => 1000,
			-MATE => $output_mate
	        	) unless(-s $fa_file2);
	print STDERR "fa_file2: *$fa_file2*\n" if($debug);
	return unless((-f $fa_file1 && -s $fa_file1) || (-f $fa_file2 && -s $fa_file2));

	my $fa_file = "$anno_dir/reads.$chr1.$pos1.$chr2.$pos2.fa";
	if($fa_file1 eq $fa_file2){
		#`cat $fa_file1 > $fa_file`;
		$fa_file = $fa_file1;
		#`cat $fa_file1.qual > $fa_file.qual` if(-s "$fa_file1.qual");
	}
	else {
		unlink $fa_file if(-s $fa_file);
		#unlink "$fa_file.qual" if(-s "$fa_file.qual");
		my $arg = "";
		$arg .= " $fa_file1 " if (-f $fa_file1 && -s $fa_file1);
		$arg .= " $fa_file2 " if (-f $fa_file2 && -s $fa_file2);
		if ($arg ne ""){
			`cat $arg >> $fa_file`;
			if ($?){
				my $err = $!;
				print STDERR "Error creating fasta file: $err\n";
				exit 8;
			}
		}
		#`cat $fa_file1 >> $fa_file` if(-f $fa_file1 && -s $fa_file1);
		#`cat $fa_file2 >> $fa_file` if(-f $fa_file2 && -s $fa_file2);
		#`cat $fa_file1.qual >> $fa_file.qual` if(-f "$fa_file1.qual" && -s "$fa_file1.qual");
		#`cat $fa_file2.qual >> $fa_file.qual` if(-f "$fa_file2.qual" && -s "$fa_file2.qual");
	}
	print STDERR "to do assembly ...\n" if($debug);
	my($contig_file, $sclip_count, $contig_reads) = $assembler->run($fa_file);

	my @mappings;
	print STDERR "start mapping ... $contig_file\n" if($debug && -s $contig_file);
	print STDERR join("\t", $chr1, $pos1, $clip1, $read_len), "\n" if($debug);
	my $ref_chr1 = normalizeChromosomeName($seq_ids[0], $chr1);
	push @mappings, $mapper->run(-QUERY => $contig_file, -scChr => $ref_chr1, -scSite=>$pos1, -CLIP=>$clip1, -READ_LEN => $read_len, "-blat-adjust-sc" => $gm) if(-s $contig_file);
	print STDERR "number of mapping: ", scalar @mappings, "\n" if($debug);
	my $ref_chr2 = normalizeChromosomeName($seq_ids[0], $chr2);
	push @mappings, $mapper->run(-QUERY => $contig_file, -scChr => $ref_chr2, -scSite=>$pos2, -CLIP=>$clip2, -READ_LEN => $read_len, "-blat-adjust-sc" => $gm) if(-s $contig_file);
	push @mappings, $mapper->run(-QUERY => $contig_file, -scChr => $ref_chr2, -scSite=>$pos2, -CLIP=>$clip2, -READ_LEN => $read_len, "-blat-adjust-sc" => $gm) if(($SV->{type} eq 'Internal_dup' || !@mappings) && -s $contig_file);

	my @qSVs;
	foreach my $sv (@mappings){
		print STDERR "\n***mapping of new contig: ", $sv->{junc_seq}, "\n" if($debug);

		my ($first_bp, $second_bp, $qseq) = ($sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq});
		my ($ortA, $chrA, $tstartA, $tendA, $qstartA, $qendA, $qstrandA, $matchesA, $percentA, $repeatA) =
		   ($first_bp->{ort}, $first_bp->{tname}, $first_bp->{tstart}, $first_bp->{tend}, $first_bp->{qstart}, $first_bp->{qend}, $first_bp->{qstrand}, $first_bp->{matches}, $first_bp->{percent}, $first_bp->{repeat});
		my ($ortB ,$chrB, $tstartB, $tendB, $qstartB, $qendB, $qstrandB, $matchesB, $percentB, $repeatB) =
		   ($second_bp->{ort}, $second_bp->{tname}, $second_bp->{tstart}, $second_bp->{tend}, $second_bp->{qstart}, $second_bp->{qend}, $second_bp->{qstrand}, $second_bp->{matches}, $second_bp->{percent}, $second_bp->{repeat});
		if($bp1->{tname} =~ m/chr/) {
			if ($chrA !~ m/^chr/){
			  $chrA = "chr".$chrA;
			}
			if ($chrB !~ m/^chr/){
			  $chrB = "chr".$chrB;
			}
		}	
		print STDERR "first_bp: ",  join("\t", $ortA, $chrA, $tstartA, $tendA, $qstartA, $qendA, $qstrandA, $matchesA, $repeatA), "\n" if($debug);
		print STDERR "second_bp: ", join("\t", $ortB, $chrB, $tstartB, $tendB, $qstartB, $qendB, $qstrandB, $matchesB, $repeatB), "\n" if($debug);
		my ($qposA, $qposB) = ($ortA > 0) ? ($qendA, $qstartB) : ($qstartA, $qendB);
		my ($clipA, $clipB) = ($ortA*$qstrandA, $ortB*$qstrandB);
		my $tposA = ($clipA > 0) ? $tendA : $tstartA;
		my $tposB = ($clipB > 0) ? $tendB : $tstartB;

		print STDERR "first_bp: ", join("\t", $ortA, $chrA, $tposA, $qstrandA), "\n" if($debug);
		print STDERR "second_bp: ", join("\t", $ortB, $chrB, $tposB, $qstrandB), "\n" if($debug);
		print STDERR "bp1: ", join("\t", $bp1->{ort}, $bp1->{tname}, $bp1->{tpos}, $bp1->{qstrand}), "\n" if($debug);
		print STDERR "bp2: ", join("\t", $bp2->{ort}, $bp2->{tname}, $bp2->{tpos}, $bp2->{qstrand}), "\n" if($debug);
	
		#next unless(($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA)<50 &&
		#	    $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB)<50) ||
		#	    ($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA)<50 &&
                #            $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB)<50));
		# to do alignment

		my $usable;
		if (($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA)<50 &&
			     $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB)<50) ||
			    ($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA)<50 &&
			     $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB)<50)) {
		  # passes original requirement, OK
		  $usable = 1;
		} else {
		  #
		  # rescue site if both ends are near a splice site
		  #
		  # - chrA/tposA/chrB/tposB are the blat align-adjusted
		  #   positions
		  $usable = breakpoint_annotation_qc($chrA, $tposA, $chrB, $tposB, $bp1, $bp2);
		};
		next unless $usable;

		my $tmp_ctg_file = "$anno_dir/reads.$chr1.$pos1.$chr2.$pos2.fa.tmp.contig";
		open(my $CTG, ">$tmp_ctg_file");
		print $CTG ">ctg\n$qseq\n";
		close($CTG);

		my ($psl_file1, $psl_file2) = ("$anno_dir/bp1.psl", "$anno_dir/bp2.psl",);
		
		unlink $psl_file1 if(-f $psl_file1); unlink $psl_file2 if(-f $psl_file2);
		if (-s $fa_file1){
			`blat -noHead -maxIntron=5 $tmp_ctg_file $fa_file1 $psl_file1`;
			if ($?){
				my $err = $!;
				print STDERR "Error running blat: $err\n";
				print STDERR "File: $fa_file1\n";
				exit 9;
			}
		}
		if (-s $fa_file2){
			`blat -noHead -maxIntron=5 $tmp_ctg_file $fa_file2 $psl_file2`;
			if ($?){
				my $err = $!;
				print STDERR "Error running blat: $err\n";
				print STDERR "File: $fa_file2\n";
				exit 10;
			}
		}
		my ($readsA, $areaA, $readsB, $areaB) = (0,1,0,1);
		my $shift_bases = 5;
		if($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA)<50 &&
		   $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB)<50){
				#$tposA = $bp1->{tpos};
				#$tposB = $bp2->{tpos};
				($readsA, $areaA) = get_junc_reads($psl_file1, $qposA, $ortA, $shift_bases) if(-f $psl_file1);
				($readsB, $areaB) = get_junc_reads($psl_file2, $qposB, $ortB, $shift_bases) if(-f $psl_file2);
		}
		elsif($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA)<50 &&
                   $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB)<50){
				#$tposA = $bp2->{tpos};
				#$tposB = $bp1->{tpos};
				($readsA, $areaA) = get_junc_reads($psl_file2, $qposA, $ortA, $shift_bases) if(-f $psl_file2);
				($readsB, $areaB) = get_junc_reads($psl_file1, $qposB, $ortB, $shift_bases) if(-f $psl_file1);
		}
		#intron size usually <100kb
		elsif($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA)<100000 &&
                   $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB)<50){
                                ($readsA, $areaA) = get_junc_reads($psl_file1, $qposA, $ortA, $shift_bases) if(-f $psl_file1);
                                ($readsB, $areaB) = get_junc_reads($psl_file2, $qposB, $ortB, $shift_bases) if(-f $psl_file2);
                }
		elsif($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA)<50 &&
                   $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB)<100000){
                                ($readsA, $areaA) = get_junc_reads($psl_file1, $qposA, $ortA, $shift_bases) if(-f $psl_file1);
                                ($readsB, $areaB) = get_junc_reads($psl_file2, $qposB, $ortB, $shift_bases) if(-f $psl_file2);
                }
                elsif($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA)<100000 &&
                   $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB)<50){
                                ($readsA, $areaA) = get_junc_reads($psl_file2, $qposA, $ortA, $shift_bases) if(-f $psl_file2);
                                ($readsB, $areaB) = get_junc_reads($psl_file1, $qposB, $ortB, $shift_bases) if(-f $psl_file1);
                }
		elsif($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA)<50 &&
                   $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB)<100000){
                                ($readsA, $areaA) = get_junc_reads($psl_file2, $qposA, $ortA, $shift_bases) if(-f $psl_file2);
                                ($readsB, $areaB) = get_junc_reads($psl_file1, $qposB, $ortB, $shift_bases) if(-f $psl_file1);
                }

	my $selected_bp1 = {
		clip => $clipA,
		ort => $ortA,
		tname => $chrA,
		tpos => $tposA,
		tstart => $tstartA,
		tend => $tendA,
		qpos => $qposA,
		qstart => $qstartA,
		qend => $qendA,
		qstrand => $qstrandA,
		matches => $matchesA,
		percent => $percentA,
		repeat => $repeatA,
		reads_num => $readsA,
		area => $areaA
	};

	my $selected_bp2 = {
		clip => $clipB,
		ort => $ortB,
		tname => $chrB,
		tpos => $tposB,
		tstart => $tstartB,
		tend => $tendB,
		qpos => $qposB,
		qstart => $qstartB,
		qend => $qendB,
		qstrand => $qstrandB,
		matches => $matchesB,
		percent => $percentB,
		repeat => $repeatB,
		reads_num => $readsB,
		area => $areaB,
	};

		my $tmp_SV = {
			junc_seq => $qseq,
			first_bp => $selected_bp1,
			second_bp => $selected_bp2
			};
		
		push @qSVs, $tmp_SV if($selected_bp1->{tpos} && $selected_bp2->{tpos});
	}
	return @qSVs;
}

sub get_junc_reads{

	my ($psl_file, $bp, $ort, $cutoff) = @_;
	my %junc_reads = ();
	my $coverage = 0;
	my $debug = 0;#Tiam
	open(hFi, $psl_file);
	my %supports = ();
	while(<hFi>){
		my $line = $_;
		chomp($line);
		my @fields = split(/\t/, $line);
		my ($matches, $qstrand, $qname, $qstart, $qend, $tstart, $tend) = ($fields[0], $fields[8], $fields[9], $fields[11], $fields[12], $fields[15], $fields[16]);
		next unless($matches > $min_hit_len);
		my $percent = $matches/($qend - $qstart);
		next if($percent <= 0.95);
		next unless($tend > $bp + $cutoff && $tstart < $bp - $cutoff);
		$junc_reads{$qname} = 1;
		$qstrand = ($qstrand eq '+') ? 1 : -1;
		my $support_len = ($ort > 0) ? ($tend - $bp) : ($bp - $tstart);
		print STDERR "$qname: support_len = ($ort > 0) ? ($tend - $bp) : ($bp - $tstart)\n" if($debug);
		print STDERR "coverage = $coverage + $support_len\n" if($debug);
		if(! exists($supports{$support_len})){
			$supports{$support_len} = 1;
		}
		else{
			next if($supports{$support_len} == 2);
			$supports{$support_len}++;
		}
		$coverage = $coverage + $support_len;
	}
	close(hFi);
	my @rtn = (scalar (keys %junc_reads), $coverage);
	return @rtn;
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

sub breakpoint_annotation_qc {
  my ($chrA, $tposA, $chrB, $tposB, $bp1, $bp2) = @_;

  my $pass_max_distance = (
			   ($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA) < POSITION_RESCUE_MAX_INTRON_DISTANCE &&
			    $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB) < POSITION_RESCUE_MAX_INTRON_DISTANCE) ||
			   ($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA) < POSITION_RESCUE_MAX_INTRON_DISTANCE &&
			    $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB) < POSITION_RESCUE_MAX_INTRON_DISTANCE)
			  );

  my $usable_count = 0;
  my @status;

  if ($pass_max_distance) {
    #
    #  only perform gene checks if shift distance is acceptable
    #
    foreach my $set (
		     [$chrA, $tposA],
		     [$chrB, $tposB]
		    ) {
      my ($chr, $pos) = @{$set};

      $chr = "chr" . $chr unless $chr =~ /^chr/;
      # ensure UCSC-style prefix is present to match refFlat

      my $usable_reason;
      my $usable_transcript;

      unless ($chr =~ /M/) {
      SEARCH:
	foreach my $strand (qw(+ -)) {
	  my $tree = $gm->sub_model($chr, $strand);
	  if (defined $tree) {
	    foreach my $model ($tree->intersect([$pos, $pos])) {
	      my $gene_ref = $model->val();
	      foreach my $transcript (@{$gene_ref->transcripts}) {
		foreach my $er (@{$transcript->exons}) {
		  my ($start, $end) = @{$er};
		  $start++;
		  # convert UCSC interbase to 1-based start/end
		  foreach my $edge ($start, $end) {
		    my $distance = abs($pos - $edge);
		    if ($distance <= POSITION_RESCUE_MAX_EXON_DISTANCE) {
		      $usable_transcript = $transcript;
		      $usable_reason = "near_exon_edge";
		      last SEARCH;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      my @info;
      if ($usable_reason) {
	push @info, $usable_reason, $usable_transcript->name(), $usable_transcript->refseq_id();
	$usable_count++;
      } else {
	push @info, "far_from_splice";
      }

      push @status, join ".", @info, $chr, $pos;
    }
  } else {
    push @status, join ".", "shift_beyond_max_intron_distance", $chrA, $tposA, $chrB, $tposB, $bp1->{tname}, $bp1->{tpos}, $bp2->{tname}, $bp2->{tpos};
  }

  my $final_call = $usable_count == 2 ? 1 : 0;
  printf STDERR "rescue:%d %s\n", $final_call, join " ", @status;

  return $final_call;
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
