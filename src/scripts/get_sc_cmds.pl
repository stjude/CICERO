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
use CiceroUtil qw(normalizeChromosomeName);

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
my $sc_shift = 3;
my $output_dir;
my $excludes_file;
my $optionOK = GetOptions(
	's|sample=s'	=> \$sample,
	'q|queue=s'		=> \$queue,
	'h|help|?'		=> \$help,
	'i|bam_d=s'	=> \$bam_file,
	'genome=s'	=> \$genome,
	'l|read_len=i'	=> \$read_length,
	'c|cluster=i'   => \$sc_shift,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
	'o|outdir=s'	=> \$output_dir,
	'e|exclude=s'   => \$excludes_file
);
if( !$bam_file) {
	croak "you must provide sample names or input bam files to run the program\n";
}

my ($genome_file, $gene_model_file, $excluded_chr, $chr_lengths);

my $conf = &TdtConfig::findConfig("genome", $genome); 
if ($conf){
	$conf = &TdtConfig::readConfig("genome", $genome); 	
	$gene_model_file = $conf->{'REFSEQ_REFFLAT'};
	$genome_file = $conf->{'FASTA'}; 
	$excluded_chr = $conf->{'EXCLD_CHR'};
	$excludes_file = $conf->{'EXCLUDED_REGIONS'}; 	
}
else{
	croak "Unknown genome name: $genome\n";
}

my $sam_d = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $genome_file);
my @seq_ids = $sam_d->seq_ids;


my @chrs = split(/,/, $excluded_chr); 
my %regions; 
if($excludes_file){
	readRegions($excludes_file);
}


open my $chrFile, "<", $conf->{'CHR_LENGTHS'};
while (my $chr = <$chrFile>){
	($chr, my $len) = split(/\s/, $chr);
	$chr = normalizeChromosomeName($seq_ids[0], $chr); 
	
	my $skip = 0; 
	foreach my $bad (@chrs){
		$skip = 1 if ($chr =~ /.*$bad.*/i); 
	}
	next if ($skip);  
	if($excludes_file){
		my $s = 0; # Intialize to start of the chromosome to search for regions.
		foreach my $start (sort {$a <=> $b} keys %{$regions{$chr}}){
			my $end = $regions{$chr}{$start};
			# If the region is 0 length, skip it
			if ($start == $end){
				$s = $end;
				next;
			} 
			# If the region start is equal to current start position, 
			# we need to skip to the end of this region. 
			if ($s == $start){
				$s = $end + 1;
			}
			else{ 
				# If the start of the region is not the current start,  
				# we need to create a region from the current start to the
				# beginning of this region.
				my $region_end = $start - 1;
				if ($s < $region_end){
					my $cmd = "extract_range.pl --ref_genome $genome_file -i $bam_file -o $output_dir -r $chr:$s-$region_end -l $read_length -m 2 -min_sc_len 3 -c $sc_shift";
					print $cmd."\n";
				}
				# The new current start is the start of this region.
				$s = $end + 1; 
			}
		}
		# If our last excluded region didn't extend to the end
		# of the chromosome, add a region. 
		if ($s < $len){
 			my $cmd = "extract_range.pl --ref_genome $genome_file -i $bam_file -o $output_dir -r $chr:$s-$len -l $read_length -m 2 -min_sc_len 3 -c $sc_shift";
			print $cmd."\n";
		} 
	}
	else{
		my $cmd = "extract_range.pl --ref_genome $genome_file -i $bam_file -o $output_dir -r $chr -l $read_length -m 2 -min_sc_len 3 -c $sc_shift";
		print $cmd."\n";
	}
}

sub readRegions{
	my ($file) = @_;
	open(my $fh, "<", $file); 
	while(my $line = <$fh>){
		chomp $line; 
		my ($chrom, $start, $end, $type) = split("\t", $line);
		$chrom = normalizeChromosomeName($seq_ids[0], $chrom);
		$regions{$chrom}{$start} = $end; 
	}
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
