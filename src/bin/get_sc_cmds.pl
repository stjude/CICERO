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

my $sample;
my $queue = $ENV{"AFC_DEFAULT_QUEUE"}; 

my ( $help, $man, $version, $usage );
my $prefix = '';
my $bam_file;
my $genome; 
my $read_length = 100;
my $output_dir;
my $optionOK = GetOptions(
	's|sample=s'	=> \$sample,
	'q|queue=s'		=> \$queue,
	'h|help|?'		=> \$help,
	'i|bam_d=s'	=> \$bam_file,
	'genome=s'	=> \$genome,
	'l|read_len=i'	=> \$read_length,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
	'o|outdir=s'	=> \$output_dir,
);
if( !$bam_file) {
	croak "you must provide sample names or input bam files to run the program\n";
}

my ($genome_file, $gene_model_file, $excluded_chr);

my $conf = &TdtConfig::findConfig("genome", $genome); 
if ($conf){
	$conf = &TdtConfig::readConfig("genome", $genome); 	
	$gene_model_file = $conf->{'REFSEQ_REFFLAT'};
	$genome_file = $conf->{'FASTA'}; 
	$excluded_chr = $conf->{'EXCLD_CHR'}; 	
}
else{
	croak "Unknown genome name: $genome\n";
}

my @chrs = split(/,/, $excluded_chr); 

open my $chrFile, "<", $conf->{'CHR_LENGTHS'};
while (my $chr = <$chrFile>){
	($chr, my $o) = split(/\s/, $chr);
	my $skip = 0; 
	foreach my $bad (@chrs){
		$skip = 1 if ($chr =~ /.*$bad.*/i); 
	}
	next if ($skip);  
	my $cmd = "extract_range.pl --ref_genome $genome_file -i $bam_file -o $output_dir -r $chr -l $read_length -m 2 -min_sc_len 3";
	print $cmd."\n";
}

=head1 LICENCE AND COPYRIGHT
Copyright 2019 St. Jude Children's Research Hospital 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
