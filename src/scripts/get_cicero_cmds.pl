#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long;
use English;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use TdtConfig; 

if (@ARGV == 0){
	pod2usage(1); 
	exit 1; 
}

my $out_prefix;
my $queue = $ENV{"AFC_DEFAULT_QUEUE"}; 

my ( $help, $man, $version, $usage );
my $bam_file;
my $genome; 
my $read_length = 100;
my $output_dir;
my $cluster_arg = 10;
my $optionOK = GetOptions(
	'out_prefix=s'	=> \$out_prefix,
	'q|queue=s'		=> \$queue,
	'h|help|?'		=> \$help,
	'i|bam_d=s'	=> \$bam_file,
	'genome=s'	=> \$genome,
	'l|read_len=i'	=> \$read_length,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
	'o|outdir=s'	=> \$output_dir,
    'c|cluster=i'   => \$cluster_arg,
);
if( !$bam_file) {
	croak "you must provide an input bam file to run the program\n";
}

my ($genome_file, $gene_model_file);

my $conf = &TdtConfig::findConfig("genome", $genome); 
if ($conf){
	$conf = &TdtConfig::readConfig("genome", $genome); 	
	$gene_model_file = $conf->{'REFSEQ_REFFLAT'};
	$genome_file = $conf->{'FASTA'}; 	
}
else{
	croak "Unknown genome name: $genome\n";
}

$out_prefix = fileparse($bam_file, '.bam') if (!$out_prefix); 

my @files = <$output_dir/$out_prefix.*.SC>;
foreach my $file (@files){
	my $out = $file; 
	$out =~ s/.SC$//;
	mkdir $out;  
	my $cmd = "Cicero.pl -genome $genome -i $bam_file -o $out -l $read_length -f $file -c $cluster_arg";
	print $cmd."\n";
}
