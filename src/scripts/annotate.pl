#!/usr/bin/env perl 

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
	is_PCR_dup read_fa_file get_discordant_reads get_sclip_reads normalizeChromosomeName);

require CiceroExtTools;

use TdtConfig; 

use Transcript;
use Gene;
use GeneModel;

my $debug = 0;

my $out_header = join("\t", "sample", "geneA", "chrA", "posA", "ortA", "featureA", "geneB", "chrB", "posB", "ortB", "featureB", 
		 	"sv_ort", "readsA", "readsB", "matchA", "matchB", "repeatA", "repeatB", , "coverageA", "coverageB",
			 "ratioA", "ratioB", "qposA", "qposB", "total_readsA", "total_readsB", "contig", "type");


my ($blat_server, $blat_port, $dir_2bit);
my $blat_client_options = ' -out=psl -nohead > /dev/null 2>&1';
my $cap3_options = " -o 25 -z 2 -h 60 > /dev/null";

my $rmdup = 0;
my $paired = 1;

# input/output
my ($genome, $ref_genome, $gene_model_file, $header);
my ($config_file, $out_dir, $gene_info_file, $junction_file, $known_itd_file, $known_fusion_file);
my ($input_bam, $read_len, $sample);
my ($all_output, $internal, $DNA) = (0, 0, 0);
my ($min_hit_len, $max_num_hits, $min_fusion_distance);
$min_hit_len = 25;
my $sc_shift = 10; 
my ($blacklist_gene_file, $blacklist_fusion_file, $complex_region_file, $excluded_chroms, $gold_gene_file);
my ( $help, $man, $version, $usage );

if(@ARGV == 0){
	#TODO: get the correct usage string
	print STDERR "Usage: $0 -g <genome> -i <bam_file> -o <out_dir> -f <gene_info>\n"; 
	exit 1; 
}

my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,	
	'c|config_file=s'	=> \$config_file,
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'  => \$genome,
	'header'    => \$header,
	'ref_genome=s'  => \$ref_genome,
	'genemodel=s'		=> \$gene_model_file,
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
	'f|gene_info_file=s' => \$gene_info_file,
	'known_itd_file=s'	=> \$known_itd_file,
	'known_fusion_file=s'	=> \$known_fusion_file,
	'j|junction_file=s' => \$junction_file,
	's|sample=s'		=> \$sample,
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
);

if ($header){
	print STDOUT $out_header; 
	exit;
}

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
$blacklist_gene_file = $conf->{BLACKLIST_GENES} unless($blacklist_gene_file);
$blacklist_fusion_file = $conf->{BLACKLIST_FUSIONS} unless($blacklist_fusion_file);
print STDERR "blacklist_fusion_file: $blacklist_fusion_file\n";
$known_itd_file = $conf->{KNOWN_ITD_FILE} unless($known_itd_file);
print STDERR "KNOWN_ITD_FILE: ", $known_itd_file, "\n";
$known_fusion_file = $conf->{KNOWN_FUSIONS} unless($known_fusion_file);
$gold_gene_file = $conf->{CLINCLS_GOLD_GENE_LIST_FILE};
print STDERR "CLINCLS_GOLD_GENE_LIST_FILE:", $gold_gene_file, "\n";
$excluded_chroms = $conf->{EXCLD_CHR} unless($excluded_chroms);
$complex_region_file = $conf->{COMPLEX_REGIONS} unless($complex_region_file);

# Load Cicero-specific configuration settings
$conf = &TdtConfig::readConfig('app', 'cicero'); 
$min_hit_len = $conf->{MIN_HIT_LEN} unless($min_hit_len);
$max_num_hits = $conf->{MAX_NUM_HITS} unless($max_num_hits);
$min_fusion_distance = $conf->{MIN_FUSION_DIST} unless($min_fusion_distance);

# Assume sample name is the bam prefix
$sample = basename($input_bam, ".bam") unless($sample);

my @excluded_chroms = split(/,/,$excluded_chroms);

# Load list of complex regions to be used later to determine if
# breakpoints are located within a "complex" region
my @complex_regions;
#if ($complex_region_file && -s $complex_region_file){
	open (my $CRF, $complex_region_file);
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
			exit $err;
		}		
	}
}
else{
	unless (-s $unfiltered_file){
		`cat $out_dir/*/unfiltered.fusion.txt > $unfiltered_file`;
		if ($?){
			my $err = $!;
			print STDERR "Error combining fusion events: $err\n"; 
			exit $err;
		}
	}
}
print "Unfiltered fusions file: $unfiltered_file\n"; 

if (! $gene_info_file || ! -e $gene_info_file){
	my $out_prefix = basename($input_bam, ".bam");
	$gene_info_file = "$out_prefix.gene_info.txt";
	$gene_info_file = File::Spec->catfile($out_dir, $gene_info_file);
}
my $out_file = $unfiltered_file;
$out_file =~ s/unfiltered/annotated/;
print STDERR "unfiltered results: $unfiltered_file\n" if($debug);

my %gene_info = ();
print STDERR "\ngene_info_file: $gene_info_file\n" if($debug);
open my $GI, "$gene_info_file" or die "cannot open < $gene_info_file: $!";
while(<$GI>){
	chomp;
	my ($name, $gRange, $strand, $mRNA_length, $cnt, $sc_cutoff) = split(/\t/);
	$gene_info{$name} = $sc_cutoff;
}
close $GI;

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

my %blacklist = ();
open(my $BLK, "$blacklist_gene_file");
while(<$BLK>){
	my $line = $_;
	chomp($line);
	$blacklist{$line} = 1;
}
close($BLK);

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
if(!$internal){
	my $df = new DelimitedFile(
	       "-file" => $blacklist_fusion_file,
	       "-headers" => 1,
	      );

	while (my $row = $df->get_hash()) {
	      my @gene_a = split(",", $row->{geneA});
	      my @gene_b = split(",", $row->{geneB});
	      my ($chrA, $posA, $chrB, $posB) = ($row->{chrA}, $row->{posA}, $row->{chrB}, $row->{posB});
	      foreach my $g1 (@gene_a){
			$g1 = join(":", $chrA, $posA) if($g1 eq 'NA');
		foreach my $g2 (@gene_b){
			$g2 = join(":", $chrB, $posB) if($g2 eq 'NA');
			my $gene_pair = ($g1 le $g2) ? join(":",$g1,$g2) : join(":",$g2,$g1);
			$bad_fusions{$gene_pair} = [$chrA, $posA, $chrB, $posB];
			}
		}
	}
}

sub is_bad_fusion{

	my ($chrA, $posA, $chrB, $posB) = @_;
	$chrA = "chr".$chrA unless($chrA =~ /chr/);
	$chrB = "chr".$chrB unless($chrB =~ /chr/);
	foreach my $xx (values %bad_fusions) {
		return 1 if($chrA eq @{$xx}[0] && $chrB eq @{$xx}[2] &&
			    abs(@{$xx}[1] - $posA) < 10000 && abs(@{$xx}[3] - $posB) < 10000);
		return 1 if($chrA eq @{$xx}[2] && $chrB eq @{$xx}[0] &&
			    abs(@{$xx}[3] - $posA) < 10000 && abs(@{$xx}[1] - $posB) < 10000);
	}
	return 0;
}

my %known_ITDs = ();
open(my $ITD_F, $known_itd_file);
while(<$ITD_F>){
	chomp;
	my ($gene, $chr, $start, $end) = split(/\t/,$_);
	$known_ITDs{$gene} = [$start, $end];
}

# Load gene pairs that are enhancer activated
# involving T-cell receptors and immunoglobulin loci
my %enhancer_activated_genes = ();
if($known_fusion_file && -e $known_fusion_file){
   open(my $KFF, $known_fusion_file);
   while(<$KFF>){
	chomp;
	my ($gene1, $gene2) = split(/\t/,$_);
	next if exists($enhancer_activated_genes{$gene1});
	next if exists($enhancer_activated_genes{$gene2});
	$enhancer_activated_genes{$gene2} = $gene1 if($gene1 =~ m/^IG.$/ || $gene1 =~ m/^TR.$/);
	$enhancer_activated_genes{$gene1} = $gene2 if($gene2 =~ m/^IG.$/ || $gene2 =~ m/^TR.$/);
   }
   close($KFF);
}

my @raw_SVs = ();
my %gene_recurrance = ();
my %contig_recurrance = ();
my %genepairs = (); 
my %breakpoint_sites = (); 
open(my $UNF, "$unfiltered_file");
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
	my $bad_fusion = 0;
	foreach my $g1 (@genes1) {
		$bad_gene = 1 if(exists($blacklist{$g1}));
		foreach my $g2 (@genes2){
			$bad_gene = 1 if(exists($blacklist{$g2}));
			next if ($g1 eq $g2);
			my $gene_pair = ($g1 le $g2) ? join(":",$g1,$g2) : join(":",$g2,$g1);
			$bad_fusion = 1 if(exists($bad_fusions{$gene_pair})); 
			last if($bad_fusion);
			next if(exists($genepairs{$gene_pair})); 
			$genepairs{$gene_pair} = 1;
			if(exists($gene_recurrance{$g1})){$gene_recurrance{$g1}++;}
			else{$gene_recurrance{$g1} = 1;}
			if(exists($gene_recurrance{$g2})){$gene_recurrance{$g2}++;}
			else{$gene_recurrance{$g2} = 1;}
		}
	}
	next if($bad_fusion || $bad_gene);

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
			 next unless ($type eq 'Internal_dup' && exists($gold_genes{$first_bp->{gene}}));
		}
		else{
			next unless ($type eq 'Internal_dup' && is_good_ITD($first_bp, $second_bp));
		}
		unless(is_good_ITD($first_bp, $second_bp)){
			next if(!$DNA && ($first_bp->{reads_num} < 5*$cutoff || $first_bp->{reads_num} < 10));
		}
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
	my $bad_fusion = 0;

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
		my $gene_pair = ($gene1 le $gene2) ? join(":",$gene1,$gene2) : join(":",$gene2,$gene1);
		$bad_fusion = 1 if(exists($bad_fusions{$gene_pair})); 
		last if($bad_fusion);
	}
	next if($bad_fusion);
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

my $New_blacklist_file = "$out_dir/blacklist.new.txt";
if ($internal){
	$New_blacklist_file = "$out_dir/blacklist.new.internal.txt";
}
open(my $NBLK, ">$New_blacklist_file");
foreach my $g (sort { $gene_recurrance{$b} <=> $gene_recurrance{$a} } keys %gene_recurrance) {
	last if($gene_recurrance{$g} < $max_num_hits*10);
	next if($g eq "NA" || $g eq "IGH" || $g eq "TCRA" || $g eq "TCRB" || $g eq "TARP");
	#next if($cr_hash{$g});
	$blacklist{$g} = 1;
	print STDERR "Adding $g to blacklist with recurrance: ".$gene_recurrance{$g}."\n"; 
	print $NBLK join("\t",$g, $gene_recurrance{$g}),"\n";
}
close($NBLK);

print STDERR "out file is: $out_file\nnumber of SVs: ", scalar @raw_SVs, "\n" if($debug);
`mkdir  -p $out_dir/tmp_anno`;
if ($?){
	my $err = $!;
	print STDERR "Error creating annotation directory: $err\n"; 
	exit $err;
}
my $annotation_dir = tempdir(DIR => "$out_dir/tmp_anno");
`mkdir -p $annotation_dir`;
if ($?){
	my $err = $!;
	print STDERR "Error creating temp directory for annotation: $err\n"; 
	exit $err;
}
print STDERR "Annotation Dir: $annotation_dir\n" if($debug); 

print STDERR "Reading cover files\n"; 
# Loop over the soft clip files again
my @cover_files = <$out_dir/*.cover>;
foreach my $fn (@cover_files) {
    my @path = split(/\//, $fn);
	open(my $IN, "$fn");
	while(<$IN>){
		chomp;
		my $line = $_;
		chomp($line);

		my ($chr, $pos, $clip, $sc_cover, $cover, $psc, $nsc, $pn, $nn) = split(/\t/,$line);
		$chr = normalizeChromosomeName($seq_ids[0], $chr);

		$clip = RIGHT_CLIP if($clip eq "+");
		$clip = LEFT_CLIP if($clip eq "-");
		# If this is a new site, add it to the list of breakpoints
		my $site= $chr."_".$pos."_".$clip;
		if(not exists($breakpoint_sites{$site})){
			$breakpoint_sites{$site} = $line;
		}
		else{ # Look for breakpoints already added that are within +/- 5bp
			my $sc_cover0 = 0;
			for(my $s=-5; $s<=5; $s++){
				my $tmp_pos = $pos + $s;
				my $tmp_site= $chr."_".$tmp_pos."_".$clip;
				# If the new site is within +/- 5bp and its coverage is greater than the prior site
				if(exists($breakpoint_sites{$tmp_site}) && $sc_cover > $sc_cover0){
					$sc_cover0 = $sc_cover;
					$breakpoint_sites{$site} = $line;
				}
			}
		}
	}
	close($IN);
}

print STDERR "Processing raw SVs\n"; 
my @annotated_SVs;
foreach my $sv (@raw_SVs){
	next if($sv->{type} eq "Internal_splicing");

	my ($first_bp, $second_bp, $contigSeq) = ($sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq});
	my ($gene1, $gene2) = ($first_bp->{gene}, $second_bp->{gene});
	my @genes1 = split(/,|\|/, $gene1);
	my @genes2 = split(/,|\|/, $gene2);
	my $bad_gene = 0;
	print STDERR "xxx\n" if(abs($sv->{second_bp}->{tpos} - 170818803)<10 || abs($sv->{first_bp}->{tpos} - 170818803)<10);
	foreach my $g1 (@genes1) {
		if(exists($blacklist{$g1})) {$bad_gene = 1; last;}
	}
	next if($bad_gene);
	foreach my $g2 (@genes2){
		if(exists($blacklist{$g2})) {$bad_gene = 1; last;}
	}
	next if($bad_gene);
	print STDERR "next if($contigSeq && ", $contig_recurrance{$contigSeq}," > $max_num_hits)\n" if(abs($sv->{second_bp}->{tpos} - 170818803)<10 || abs($sv->{first_bp}->{tpos} - 170818803)<10); 

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
				next unless (is_good_ITD($first_bp, $second_bp) || exists($gold_genes{$first_bp->{gene}}));
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
print "annotated_SVs: ", scalar @annotated_SVs, "\n"; 

my @uniq_SVs;
foreach my $sv (@annotated_SVs){
	print STDERR "xxx\n" if(abs($sv->{second_bp}->{tpos} - 170818803)<10 || abs($sv->{first_bp}->{tpos} - 170818803)<10);
	my ($bp1, $bp2, $qseq) = ($sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq});
	if(exists($blacklist{$bp1->{gene}}) || exists($blacklist{$bp2->{gene}})){
		print STDERR "Removing duplicate: gene1: ".$bp1->{gene}." gene2: ".$bp2->{gene}."\n";
		next;
	}
	if($sv && ! is_dup_SV(\@uniq_SVs, $sv)){
		push @uniq_SVs, $sv;
	}
	else{
		print STDERR "Duplicate SV: gene1: ".$bp1->{gene}." gene2: ".$bp2->{gene}."\n";
	}
}

print "unique_SVs ", scalar @uniq_SVs, "\n"; 

open(hFo, ">$out_file");
print hFo $out_header, "\n";
foreach my $sv (@uniq_SVs){
	my ($bp1, $bp2, $qseq, $type) = ($sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq}, $sv->{type});
	my ($geneA, $geneB) = ($bp1->{gene}, $bp2->{gene});
	my $ratioA = (exists($gene_info{$geneA}) && $gene_info{$geneA} > 0 && $bp1->{feature} ne 'intergenic') ? 
		     (($bp1->{reads_num}+0.01)/80)/$gene_info{$geneA} : ($bp1->{reads_num}+0.01)/(count_coverage($sam_d, $bp1->{tname}, $bp1->{tpos}) + 1);
	my $ratioB = (exists($gene_info{$geneB}) && $gene_info{$geneB} > 0 && $bp2->{feature} ne 'intergenic') ? 
		     (($bp2->{reads_num}+0.01)/80)/$gene_info{$geneB} : ($bp2->{reads_num}+0.01)/(count_coverage($sam_d, $bp2->{tname}, $bp2->{tpos}) + 1);
	$ratioA = 1 if($ratioA > 1);  $ratioB = 1 if($ratioB > 1);

	$bp1->{qstrand} =  ($bp1->{qstrand}>0) ? '+' : '-';
	$bp2->{qstrand} =  ($bp2->{qstrand}>0) ? '+' : '-';

	my $bp1_site = $bp1->{tname}."_".$bp1->{tpos}."_". $bp1->{clip};
	my $bp2_site = $bp2->{tname}."_".$bp2->{tpos}."_". $bp2->{clip};

	my ($mafA,$mafB) = ($ratioA, $ratioB);
	my ($pmafA,$nmafA,$pmafB,$nmafB) = (0, 0, 0, 0);
	my ($pscA, $nscA, $pnA, $nnA, $pscB, $nscB, $pnB, $nnB) = 0;

	if(exists($breakpoint_sites{$bp1_site}) && $breakpoint_sites{$bp1_site} ne "1"){
		my @bp1_fields = split(/\t/,$breakpoint_sites{$bp1_site});
		($pscA, $nscA, $pnA, $nnA) = @bp1_fields[5,6,7,8];
	}
	else{
	    for (my $s = -5; $s<=5; $s++){
		my $tmp_pos = $bp1->{tpos}+$s;
		$bp1_site = $bp1->{tname}."_".$tmp_pos."_". $bp1->{clip};

	       if(exists($breakpoint_sites{$bp1_site}) && $breakpoint_sites{$bp1_site} ne "1"){
		 my @bp1_fields = split(/\t/,$breakpoint_sites{$bp1_site});
		 ($pscA, $nscA, $pnA, $nnA) = @bp1_fields[5,6,7,8];
		 last;
	       }
	    }
	}
	$pmafA = $pscA / $pnA if($pnA); $nmafA = $nscA / $nnA if($nnA);
	$mafA = ($pmafA > $nmafA) ? $pmafA : $nmafA if($pnA || $nnA);
	$mafA = 1 if($mafA > 1);

	if(exists($breakpoint_sites{$bp2_site}) && $breakpoint_sites{$bp2_site} ne "1"){
		my @bp2_fields = split(/\t/,$breakpoint_sites{$bp2_site});
		($pscB, $nscB, $pnB, $nnB) = @bp2_fields[5,6,7,8];
	}
	else{
	    for (my $s = -5; $s<=5; $s++){
		my $tmp_pos = $bp2->{tpos}+$s;
		$bp2_site = $bp1->{tname}."_".$tmp_pos."_". $bp2->{clip};
	        if(exists($breakpoint_sites{$bp2_site}) && $breakpoint_sites{$bp2_site} ne "1"){
		   my @bp2_fields = split(/\t/,$breakpoint_sites{$bp2_site});
		   ($pscB, $nscB, $pnB, $nnB) = @bp2_fields[5,6,7,8];
		   last;
	   	}
	    }
	}
	$pmafB = $pscB / $pnB if($pnB); $nmafB = $nscB / $nnB if($nnB);
	$mafB = ($pmafB > $nmafB) ? $pmafB : $nmafB if($pnB || $nnB);
	$mafB = 1 if($mafB > 1);

	my ($total_readsA, $total_readsB) = (0,0);
	if($pnA && $nnA){
		$total_readsA = $pnA + $nnA;
	}
	else {
		$total_readsA = count_coverage($sam_d, $bp1->{tname}, $bp1->{tpos});
	}
	$total_readsA = $bp1->{reads_num} if($bp1->{reads_num} > $total_readsA);

	if($pnB && $nnB){
		$total_readsB = $pnB + $nnB;
	}
	else {
		$total_readsB = count_coverage($sam_d, $bp2->{tname}, $bp2->{tpos});
	}
	$total_readsB = $bp2->{reads_num} if($bp2->{reads_num} > $total_readsB);

	unless($seq_ids[0] =~ m/chr/) {$bp1->{tname} = "chr".$bp1->{tname}; $bp2->{tname} = "chr".$bp2->{tname};}
	my $out_string = join("\t", $sample, $bp1->{gene}, $bp1->{tname}, $bp1->{tpos}, $bp1->{qstrand}, $bp1->{feature}, 
				$bp2->{gene}, $bp2->{tname}, $bp2->{tpos}, $bp2->{qstrand}, $bp2->{feature}, $sv->{ort}, 
				$bp1->{reads_num}, $bp2->{reads_num}, $bp1->{matches}, $bp2->{matches}, sprintf("%.2f", $bp1->{repeat}), 
				sprintf("%.2f", $bp2->{repeat}), $bp1->{area}, $bp2->{area}, sprintf("%.2f", $mafA), sprintf("%.2f", $mafB),
				 $bp1->{qpos}, $bp2->{qpos}, $total_readsA, $total_readsB, $qseq, $type);
	print hFo $out_string, "\n";
}	
close(hFo);
rmtree(["$annotation_dir"]);

sub is_good_ITD {
	my($bp1, $bp2) = @_;
	my $gene = $bp1->{gene};
	return 1 if(%known_ITDs && exists($known_ITDs{$gene}) && 
	   	 $bp1->{tpos} > $known_ITDs{$gene}[0] && 
	   	 $bp1->{tpos} < $known_ITDs{$gene}[1] && 
	   	 $bp2->{tpos} > $known_ITDs{$gene}[0] && 
	   	 $bp2->{tpos} < $known_ITDs{$gene}[1]);
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

sub is_dup_SV {
	my($r_SVs, $sv) = @_;
	foreach my $s (@{$r_SVs}) {
		my $more_reads = ($s->{first_bp}->{reads_num} + $s->{second_bp}->{reads_num} >= $sv->{first_bp}->{reads_num} + $sv->{second_bp}->{reads_num}) ? 1 : 0;
		my $longer_contig = ($s->{first_bp}->{matches} + $s->{second_bp}->{matches} >= $sv->{first_bp}->{matches} + $sv->{second_bp}->{matches}) ? 1 : 0;
		return 1
		if( 	($more_reads || $longer_contig) &&
			abs($s->{first_bp}->{tpos} - $sv->{first_bp}->{tpos}) < 10 &&
			abs($s->{second_bp}->{tpos} - $sv->{second_bp}->{tpos}) < 10 &&
			$s->{first_bp}->{tname} eq $sv->{first_bp}->{tname} &&
			$s->{second_bp}->{tname} eq $sv->{second_bp}->{tname});
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

sub is_bad_chrom{
	my $chr = shift;
	foreach my $bad_chr (@excluded_chroms){
		return 1 if($chr =~ /$bad_chr/i);
	}
	return 0;
}

sub count_genes {
	my ($chr,$start,$end, $transcript_strand) = @_;
	my $cnt = 0;
	my $debug = 0;
	$transcript_strand = $transcript_strand > 0 ? "+" : "-" if($transcript_strand);
	foreach my $strand( "+", "-" ) {
		next if($transcript_strand && $strand ne $transcript_strand);
		my $full_chr = ($chr=~/chr/) ? $chr : "chr".$chr;
		my $tree = $gm->sub_model($full_chr, $strand);
		next if(!$tree);
		print STDERR "my \@tmp = $tree->intersect([$start, $end])\n" if($debug);
		if($start > $end) {my $tmp = $start; $start = $end; $end = $tmp;}
		my @tmp = $tree->intersect([$start, $end]);
		print STDERR $strand, "\tnumber of genes: ", scalar @tmp, "\n" if($debug);
		foreach my $tnode (@tmp) {
			my $g=$tnode->val;
			my $gRange = $chr.":".$g->start."-".$g->end;
			print STDERR $g->name, "\tgRange = $gRange\tstart:$start\tend:$end\n" if($debug);
			return 0 if($g->start <= $start + 10000 && $g->end >= $end - 10000);
			if($g->start > $start && $g->end < $end && $g->get_cds_length > 300){
				# to solve the overlapping gene problem.
				my $overlap = 0;
				foreach my $t (@tmp){
					my $tg = $t->val;
					next if($tg->name eq $g->name);
					$overlap = 1 if($tg->start < $g->start && $tg->end > $g->end);
				}
				$cnt++ unless($overlap);
			}
		}
	}
	return $cnt;
}

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
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
	
			if($tmp_score > $bp->{annotate_score} || ($tmp_dist < $dist && $bp->{annotate_score} == 0.1)){
				$bp->{annotate_score} = $tmp_score;
				$bp->{feature} = $tmp_feature;
				$bp->{ts_strand} = $bp->{qstrand};
				$bp->{gene} = $g->name;
				$dist = $tmp_dist if($bp->{annotate_score} == 0.1);
			}
			return $bp if($bp->{annotate_score} == 1); 
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
			if(abs($tmp_score) > abs($bp->{annotate_score}) || ($tmp_dist < $dist && abs($bp->{annotate_score}) == 0.1)) {
				$bp->{annotate_score} = $tmp_score;
				$bp->{feature} = $tmp_feature;
				$bp->{ts_strand} = -1*$bp->{qstrand};
				$bp->{gene} = $g->name;
				$dist = $tmp_dist if(abs($bp->{annotate_score}) == 0.1);
			}
			return $bp if(abs($bp->{annotate_score}) == 1); 
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
		print STDERR "bp1->ort ", $bp1->{ort}, " error!\n" if($debug);
		exit;
	}
	$clip1 = $bp1->{ort}*$bp1->{qstrand};

	# right clip or left clip?
	my $clip2;
	if($bp2->{ort} != 1 && $bp2->{ort} != -1){
		print STDERR "bp2->ort ", $bp2->{ort}, " error!\n" if($debug);
		exit;
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
	prepare_reads_file(-OUT => $fa_file1,
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
	prepare_reads_file(-OUT => $fa_file2,
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
			  exit $err;
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
	push @mappings, $mapper->run(-QUERY => $contig_file, -scChr => $ref_chr1, -scSite=>$pos1, -CLIP=>$clip1, -READ_LEN => $read_len) if(-s $contig_file);
	print STDERR "number of mapping: ", scalar @mappings, "\n" if($debug);
	my $ref_chr2 = normalizeChromosomeName($seq_ids[0], $chr2);
	push @mappings, $mapper->run(-QUERY => $contig_file, -scChr => $ref_chr2, -scSite=>$pos2, -CLIP=>$clip2, -READ_LEN => $read_len) if(-s $contig_file);
	push @mappings, $mapper->run(-QUERY => $contig_file, -scChr => $ref_chr2, -scSite=>$pos2, -CLIP=>$clip2, -READ_LEN => $read_len) if(($SV->{type} eq 'Internal_dup' || !@mappings) && -s $contig_file);

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
	
		next unless(($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA)<50 &&
			    $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB)<50) || 
			    ($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA)<50 &&
                            $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB)<50));
		# to do alignment
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
				exit $err;
			}
		}
		if (-s $fa_file2){
			`blat -noHead -maxIntron=5 $tmp_ctg_file $fa_file2 $psl_file2`;
			if ($?){
				my $err = $!;
				print STDERR "Error running blat: $err\n"; 
				exit $err;
			}
		}
		my ($readsA, $areaA, $readsB, $areaB) = (0,1,0,1);
		my $shift_bases = 5;
		if($chrA eq $bp1->{tname} && abs($bp1->{tpos} - $tposA)<50 &&
		   $bp2->{tname} eq $chrB && abs($bp2->{tpos} - $tposB)<50){ 
				$tposA = $bp1->{tpos};
				$tposB = $bp2->{tpos};
				($readsA, $areaA) = get_junc_reads($psl_file1, $qposA, $ortA, $shift_bases) if(-f $psl_file1);
				($readsB, $areaB) = get_junc_reads($psl_file2, $qposB, $ortB, $shift_bases) if(-f $psl_file2);
		}
		elsif($bp2->{tname} eq $chrA && abs($bp2->{tpos} - $tposA)<50 &&
                   $bp1->{tname} eq $chrB && abs($bp1->{tpos} - $tposB)<50){
				$tposA = $bp2->{tpos};
				$tposB = $bp1->{tpos};
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
	my $debug = 0;
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

sub get_genes {

	my ($gm, $chr, $start, $end) = @_;
	my ($f_tree, $r_tree) = ($gm->sub_model($chr, "+"), $gm->sub_model($chr, "-"));
	my (%genes, @f_genes, @r_genes);
	push @f_genes, $f_tree->intersect([$start, $end]); 
	push @r_genes, $r_tree->intersect([$start, $end]); 
	if(!@f_genes && !@r_genes){ 
		push @f_genes, $f_tree->intersect([$start-5000, $end+5000]); 
		push @r_genes, $r_tree->intersect([$start-5000, $end+5000]); 
	}
	if(!@f_genes && !@r_genes){ 
		push @f_genes, $f_tree->intersect([$start-10000, $end+10000]); 
		push @r_genes, $r_tree->intersect([$start-10000, $end+10000]); 
	}
	return (\@f_genes, \@r_genes);
}

sub get_gene_name {
	my ($gm, $chr, $start, $end) = @_;
	my ($f_tree, $r_tree) = ($gm->sub_model($chr, "+"), $gm->sub_model($chr, "-"));
	my (%genes, @f_genes, @r_genes);
	push @f_genes, $f_tree->intersect([$start, $end]); 
	push @r_genes, $r_tree->intersect([$start, $end]); 
	if(!@f_genes && !@r_genes){ 
		push @f_genes, $f_tree->intersect([$start-5000, $end+5000]); 
		push @r_genes, $r_tree->intersect([$start-5000, $end+5000]); 
	}
	if(!@f_genes && !@r_genes){ 
		push @f_genes, $f_tree->intersect([$start-10000, $end+10000]); 
		push @r_genes, $r_tree->intersect([$start-10000, $end+10000]); 
	}
	if(!@f_genes && !@r_genes){ 
		push @f_genes, $f_tree->intersect([$start-40000, $end+40000]); 
		push @r_genes, $r_tree->intersect([$start-40000, $end+40000]); 
	}
	if(!@f_genes && !@r_genes) {
		$genes{'NA'} = 0;
		return \%genes;
	};

	foreach my $g (@f_genes){
		my @gene_names = split(/,|\|/,$g->val->name);
		foreach my $g1 (@gene_names){
			$genes{$g1} = 1;
		} 
	}	

	foreach my $g (@r_genes){
		my @gene_names = split(/,|\|/,$g->val->name);
		foreach my $g1 (@gene_names){
			$genes{$g1} = -1;
		} 
	}	
	return \%genes;
}

sub get_type {

	my ($first_bp, $second_bp, $same_gene, $transcript_strand) = @_;
	my $debug = 0;
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
