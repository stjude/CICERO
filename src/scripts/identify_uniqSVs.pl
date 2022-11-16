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

my $out_header = join("\t", "sample", "geneA", "chrA", "posA", "ortA", "featureA", "geneB", "chrB", "posB", "ortB", "featureB",
		 	"sv_ort", "readsA", "readsB", "matchA", "matchB", "repeatA", "repeatB", , "coverageA", "coverageB",
			 "ratioA", "ratioB", "qposA", "qposB", "total_readsA", "total_readsB", "contig", "type");

# input/output
my ($genome, $ref_genome, $header);
my ($out_dir,$gene_info_file, $known_fusion_file);
my ($input_bam, $annotated_file, $sample);
my ($excluded_gene_file);
my ($help, $man, $version, $usage );

if(@ARGV == 0){
	#TODO: get the correct usage string
	print STDERR "Usage: $0 -g <genome> -i <bam_file> -o <out_dir> -f <gene_info>\n";
	exit 1;
}

my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,	
	'a|annotated_file=s'	=> \$annotated_file,
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'  => \$genome,
	'header'    => \$header,
	'ref_genome=s'  => \$ref_genome,
	'excluded-genes=s'	=> \$excluded_gene_file,

	'f|gene_info_file=s' => \$gene_info_file,
	'known_fusion_file=s'	=> \$known_fusion_file,
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

$excluded_gene_file = $conf->{EXCLUDED_GENES} unless($excluded_gene_file);
$known_fusion_file = $conf->{KNOWN_FUSIONS} unless($known_fusion_file);

# Load Cicero-specific configuration settings
$conf = &TdtConfig::readConfig('app', 'cicero');

# Assume sample name is the bam prefix
$sample = basename($input_bam, ".bam") unless($sample);

if (! $gene_info_file || ! -e $gene_info_file){
	my $out_prefix = basename($input_bam, ".bam");
	$gene_info_file = "$out_prefix.gene_info.txt";
	$gene_info_file = File::Spec->catfile($out_dir, $gene_info_file);
}

my $uniq_file = $annotated_file;
$uniq_file =~ s/quantified/annotated/;

my %gene_info = ();
print STDERR "\ngene_info_file: $gene_info_file\n" if($debug);
open my $GI, "$gene_info_file" or die "cannot open < $gene_info_file: $!";
while(<$GI>){
	chomp;
	my ($name, $gRange, $strand, $mRNA_length, $cnt, $sc_cutoff) = split(/\t/);
	$gene_info{$name} = $sc_cutoff;
}
close $GI;

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

my %breakpoint_sites = ();

### end of annotate.pl head

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

open my $ASV, '<',  "$annotated_file";
chomp(my @lines = <$ASV>);
close $ASV;

my @annotated_SVs;
foreach my $line (@lines){

	my @fields = split("\t", $line);
	my $first_bp = {
		reads_num => $fields[2],
		gene => $fields[0],
		tpos => $fields[4],
		ort => $fields[6],
		tname => $fields[7],
		qpos => $fields[8],
		qstart => $fields[9],
		qend => $fields[10],
		qstrand => $fields[11],
		ts_strand => $fields[12],
		matches => $fields[13],
		percent => $fields[14],
		repeat => $fields[15],
		clip => $fields[16],
		area => $fields[17],
		feature => $fields[18],
		annotate_score => $fields[19]
	};

	my $second_bp = {
		reads_num => $fields[3],
		gene => $fields[1],
		tpos => $fields[5],
		ort => $fields[20],
		tname => $fields[21],
		qpos => $fields[22],
		qstart => $fields[23],
		qend => $fields[24],
		qstrand => $fields[25],
		ts_strand => $fields[26],
		matches => $fields[27],
		percent => $fields[28],
		repeat => $fields[29],
		clip => $fields[30],
		area => $fields[31],
		feature => $fields[32],
		annotate_score => $fields[33]
	};

	my $sv = {
		first_bp => $first_bp,
		second_bp => $second_bp,
		};

	$sv->{junc_seq} = $fields[34];
	$sv->{type} = $fields[35];
	$sv->{ort} = $fields[36];

	push @annotated_SVs, $sv;
}

print STDERR "Processing annotated SVs\n";
my @uniq_SVs;
foreach my $sv (@annotated_SVs){
	print STDERR "xxx\n" if(abs($sv->{second_bp}->{tpos} - 170818803)<10 || abs($sv->{first_bp}->{tpos} - 170818803)<10);
	my ($bp1, $bp2, $qseq) = ($sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq});
	if(exist_multiplename_checking(\%excluded, $bp1->{gene}) || exist_multiplename_checking(\%excluded, $bp2->{gene})){
		# If the highly recurrent fusion doesn't involve known partners, remove it.
		if(!(exist_multiplename_pair_checking(\%known_fusion_partners, $bp1->{gene}, $bp2->{gene}))){
			print STDERR "Removing duplicate: gene1: ".$bp1->{gene}." gene2: ".$bp2->{gene}."\n";
			next;
		}
	}
	if($sv) {	#&& ! is_dup_SV(\@uniq_SVs, $sv)){
		push @uniq_SVs, $sv;
	}
	else{
		print STDERR "Duplicate SV: gene1: ".$bp1->{gene}." gene2: ".$bp2->{gene}."\n";
	}
}

print "unique_SVs ", scalar @uniq_SVs, "\n";

my @out_strings = ();
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
	push @out_strings, $out_string
}


#print hFo $out_string, "\n";

open(hFo, ">$uniq_file");
print hFo $out_header, "\n";
print hFo join ("\n", @out_strings);
print hFo "\n";
close(hFo);

### start of annotate.pl tail

sub exist_multiplename_pair_checking {
	my %fusionpartner_list = %{(shift)};
	my $targetgene1 = shift;#e.g. targetgene UBTF,MIR6782
	my $targetgene2 = shift;#e.g. targetgene UBTF,MIR6782
	my @genes1 = split(/,|\|/, $targetgene1);
	my @genes2 = split(/,|\|/, $targetgene2);

	foreach my $g1 (@genes1) {
		foreach my $g2 (@genes2){
			return 1 if(exists($known_fusion_partners{$g1}{$g2}));
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
		if( 	($more_reads || $longer_contig) &&
			abs($s->{first_bp}->{tpos} - $sv->{first_bp}->{tpos}) < 10 &&
			abs($s->{second_bp}->{tpos} - $sv->{second_bp}->{tpos}) < 10 &&
			$s->{first_bp}->{tname} eq $sv->{first_bp}->{tname} &&
			$s->{second_bp}->{tname} eq $sv->{second_bp}->{tname}) {
			return 1;
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
