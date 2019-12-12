#!/usr/bin/env perl

use warnings; 
use strict;

use Carp;
use Getopt::Long;
use English;
use Pod::Usage;
use Data::Dumper;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants;
use File::Temp qw/ tempfile tempdir /;
use File::Spec;
use File::Basename;
use Cwd;
use List::MoreUtils qw/ uniq /;
use lib dirname($0);
my $script_dir = dirname($0);
#custom packages
#use SCValidator qw($lowqual_cutoff $min_percent_id $min_percent_hq LEFT_CLIP RIGHT_CLIP);
use CiceroUtil qw(parse_range is_PCR_dup);
use TdtConfig; 

use constant FQ_BASE_NUMBER => 33;

use GeneModel;
use Gene;
my ($blacklist_file, $excluded_chroms, $genome, $ref_genome, $gene_model_file);
my $rmdup = 0;
my $out_dir;
my $read_len;

my ($input_bam, $sample);
my $paired = 1;
my ( $help, $man, $version, $usage );
my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,	
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'  => \$genome,
	'genemodel=s'		=> \$gene_model_file,
	'paired!'		=> \$paired,
	'rmdup!'		=> \$rmdup,
	'l|read_len=i'	=> \$read_len,
	's|sample=s'		=> \$sample,
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );

if($input_bam) {
	$input_bam = File::Spec->rel2abs($input_bam);
	$sample = basename($input_bam, ".bam") unless($sample);
}
else{
	croak "You must specify the input bam file or sample name";
}

my $input_base = fileparse($input_bam, '.bam');

#setup output dir and workind directory
$out_dir = getcwd if(!$out_dir);

chomp $genome; 
my $conf = &TdtConfig::findConfig("genome", $genome);
if ($conf){
	$conf = &TdtConfig::readConfig("genome", $genome); 	
	$gene_model_file = $conf->{'REFSEQ_REFFLAT'} unless($gene_model_file);
	$ref_genome = $conf->{'FASTA'}; 
	$blacklist_file = $conf->{'BLACKLIST_GENES'}; 	
	$excluded_chroms = $conf->{EXCLD_CHR} unless($excluded_chroms);
}
else{
	croak "Unknown genome name: $genome, $conf\n";
}

my $gm_format = "REFFLAT";
print "GMF: $gene_model_file\n";
my $gm = GeneModel->new if($gene_model_file);
$gm->from_file($gene_model_file, $gm_format);
my $sam = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $ref_genome);

my $out_prefix = $input_base;
my @path = split /\//, $input_bam;
my $bam_file = $path[-1];
my @bam_file=split /-/, $bam_file;
$sample=$bam_file[0] unless (defined $sample);

`mkdir -p $out_dir`;

my %blacklist = ();
if($blacklist_file && -s $blacklist_file){
open(my $BLK, $blacklist_file);
while(<$BLK>){
	my $line = $_;
	chomp($line);
	$blacklist{$line} = 1;
}
close($BLK);
}

my @chroms = $sam->seq_ids;
my @excluded_chroms = split(/,/,$excluded_chroms);
my $gene_info_file;
$gene_info_file = "$out_prefix.gene_info.txt";
$gene_info_file = File::Spec->catfile($out_dir, $gene_info_file);
open my $GI, ">$gene_info_file";
foreach my $chr (@chroms){
	next if(is_bad_chrom($chr));
	my $full_chr = ($chr=~m/chr/) ? $chr : "chr".$chr;
	my ($start, $end) = (1, $sam->length($chr));

	foreach my $strand( "+", "-" ) {
		my $tree = $gm->sub_model($full_chr, $strand);
		next unless($tree);
		my @tmp = $tree->intersect([$start, $end]);
		#print "number of genes: ", scalar @tmp, " at ", $full_chr, " ", $strand, "\n";
		foreach my $tnode (@tmp) {
			my $g=$tnode->val;
			#$gRange =~ s/chr//;
			my $mRNA_length = $g->get_mRNA_length;
			my $gRange = $chr.":".$g->start."-".$g->end;
			my $cnt = 0;
			$sam->fetch($gRange, 
				sub {
					my $a = shift;
					return if( ($a->flag & 0x0400) || ($a->flag & 0x0004) ); #PCR dumplicate or unmapped
					$cnt++;
				}) unless(exists($blacklist{$g->name}));

			my $sc_cutoff = ($read_len-20)*$cnt/(100*$mRNA_length);
			my $tmp_gene = {
				name        => $g->name, 
			   	range       => $gRange, 
			   	strand      => $strand, 
				mRNA_length => $mRNA_length, 
				reads_cnt   => $cnt,
				sc_cutoff   => $sc_cutoff
			  };
			#print STDERR join("\t", $tmp_gene->{name}, $tmp_gene->{range}, $tmp_gene->{strand}, $tmp_gene->{mRNA_length}, $tmp_gene->{reads_cnt}, $tmp_gene->{sc_cutoff}), "\n";
			print $GI join("\t", $tmp_gene->{name}, $tmp_gene->{range}, $tmp_gene->{strand}, $tmp_gene->{mRNA_length}, $tmp_gene->{reads_cnt}, $tmp_gene->{sc_cutoff}), "\n";
		}
	}
} #end for chr
close $GI;
undef $gm;

sub calc_sc_cutoff{
	my ($sam, $range, $mRNA_length) = @_;
	my $cnt=0;
	$sam->fetch($range, 
		sub {
			my $a = shift;
			return if( ($a->flag & 0x0400) || ($a->flag & 0x0004) ); #PCR dumplicate or unmapped
			$cnt++;
		});
		my $RPK_value = 1000*$cnt/$mRNA_length;
		my $sc_cutoff = ($read_len-20)*$cnt/(100*$mRNA_length);
		#print join("\t", $range, $mRNA_length, $cnt, $RPK_value, $sc_cutoff), "\n";
		return $sc_cutoff;
}

sub is_bad_chrom{
	my $chr = shift;
	foreach my $bad_chr (@excluded_chroms){
		return 1 if($chr =~ /$bad_chr/i);
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
