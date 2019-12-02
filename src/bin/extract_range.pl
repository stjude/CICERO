#!/usr/bin/env perl

use strict;
use warnings;

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
#custom packages
use CiceroSCValidator qw($lowqual_cutoff $min_percent_hq LEFT_CLIP RIGHT_CLIP);
use CiceroUtil qw(parse_range is_PCR_dup);

use constant FQ_BASE_NUMBER => 33;

my $debug = 0;
my $rmdup = 0;
my ($out_dir, $ref_genome);
my $min_sc_reads = 2;
my($read_len, $min_sc_len) = (100, 10);
my $sc_shift = 3;

# input/output
my ($out_prefix, $range, $input_bam );
my $out_suffix = "sclip.txt";
my $paired = 1;
my ( $help, $man, $version, $usage );
my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,	
	'o|out_dir=s'	=> \$out_dir,
	'ref_genome=s'  => \$ref_genome,
	'p=s'			=> \$out_prefix,
	'l|read_len=i'	=> \$read_len,
	'paired!'		=> \$paired,
	'rmdup!'		=> \$rmdup,
	'min_sc_len=i'	=>	\$min_sc_len,
	'm|min_sc_reads=i'	=>	\$min_sc_reads,
	'lq_cutoff=i'	=> \$lowqual_cutoff,
	'min_pct_hq=i'	=> \$min_percent_hq,
	'r|range=s'		=> \$range,
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );

my $start_dir = getcwd;
if($input_bam) {
	$input_bam = File::Spec->rel2abs($input_bam);
}
else{
	croak "You must specify the input bam file or sample name";
}
croak "You must provide the reference genome in fasta format!" if(!$ref_genome);

my $input_base = fileparse($input_bam, '.bam');

#setup output dir and working directory
$out_dir = getcwd if(!$out_dir);
mkdir $out_dir if(!-e $out_dir || ! -d $out_dir);

# figure out output prefix
$out_prefix = $input_base if(!$out_prefix);

my $sam = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $ref_genome);
my ($chr, $start, $end)=split /:|-/,$range;
print "range is ", join("\t", $chr, $start, $end), "\nmin_sc_reads: $min_sc_reads\tmin_sc_len:$min_sc_len\n", if($debug);

my $output_file;
my $validator = CiceroSCValidator->new();
if(!$paired) {
	$validator->remove_validator("strand_validator");
}

my $tmp = $chr;
$tmp = $tmp . ".$start" if($start);
$tmp = $tmp . ".$end" if($end);
$output_file = join('.', $out_prefix, $tmp, "cover");
$output_file = File::Spec->catfile($out_dir, $output_file);
print "output_file: $output_file\n" if($debug);
open my $SC_file, ">$output_file";
#my ($cov) = $sam->features(-type => 'coverage', -flags=>{DUPLICATE=>0}, -seq_id=> $chr,
       #       		        -start => $start, -end => $end);
        #my @c_d = $cov->coverage;
my($pcover, $ncover) = extract_range_sclip(
	-SAM => $sam, 
	-RANGE => $range, 
	-VALIDATOR => $validator);

foreach my $p (keys(%{$pcover})) {
	#my $c = ($rmdup ? scalar(@{$pcover->{$p}}) : $pcover->{$p});
	my @c = @{$pcover->{$p}};
	my $tc = count_coverage($sam, $chr, $p);
	print $SC_file join("\t", $chr, $p, "+", @c), "\t", $tc->{"+"},"\t", $tc->{"-"}, "\n" if($c[0]>=$min_sc_reads);
}
foreach my $p (keys(%{$ncover})) {
	#my $c = ($rmdup ? scalar(@{$ncover->{$p}}) : $ncover->{$p});
	my @c = @{$ncover->{$p}};
	my $tc = count_coverage($sam, $chr, $p);
	print $SC_file join("\t", $chr, $p, "-", @c), "\t", $tc->{"+"},"\t", $tc->{"-"}, "\n" if($c[0]>=$min_sc_reads);
}
exit(0);

#sub count_coverage {
#	my ($sam, $chr, $pos, $clip) = @_;
#	my ($c) = $sam->features(-type => 'coverage', -seq_id=> $chr, 
#		-start => $pos, -end => $pos);
#	my @c_d = $c->coverage;
#	return $c_d[0];
#}

sub count_coverage {
	my ($sam, $chr, $pos) = @_;
	my $seg = $sam->segment(-seq_id => $chr, -start => $pos, -end => $pos);
	#print STDERR "seg: ", $seg, "\n";
	my ($pos_N, $neg_N) = (0, 0);
	if($seg){
		my $itr = $seg->features(-iterator => 1);
		while( my $a = $itr->next_seq) {
			next unless($a->start && $a->end);
			if($a->strand > 0){
				$pos_N++;
			}
			else{
				$neg_N++;
			}
		}
	}
	my %coverage=();
	$coverage{"+"} = $pos_N;
	$coverage{"-"} = $neg_N;
	return \%coverage;
}

sub extract_range_sclip {
	my %arg = @_;
	my $debug = 0;
	my $sam = $arg{-SAM} || croak "missing -SAM";
	my $range = $arg{-RANGE} || croak "missing -RANGE";
#	my $output_file = $arg{-OUTPUT} || croak "missing -OUTPUT";
	my $validator = $arg{-VALIDATOR} || croak "missing -VALIDATOR";
	my ($chr, $start, $end) = split /[:|-]/, $range;

	my (%right_cover, %left_cover, %right_count, %left_count, %neg_right_count, %neg_left_count, %plus_right_count, %plus_left_count);
	my (%right_sclip_lens, %left_sclip_lens);
	my (%long_right_cover, %long_left_cover);
	print STDERR "$range\n" if($debug);
	$sam->fetch($range, 
		sub {
			my $a = shift;
			my $cigar_str = $a->cigar_str;
			my @cigar_array = @{$a->cigar_array};
			my $cigar_length = scalar @cigar_array;

			return if($a->cigar_str !~ m/S/);
			return if($paired && !$a->mpos); 
			return if($a->flag & 0x0400); #PCR duplicate
			return unless($a->start);
			#return unless(abs($a->end - 29446394) < 10);
			#return unless(abs($a->start - 73800850) < 10);
			#return unless($a->qname =~ m/C0FVUACXX120313:3:2302:5930:71041/);
			#return unless($a->qname =~ m/16389/);
			if($paired && !$a->proper_pair){ #paired but mate is not mapped 
				$validator->remove_validator("strand_validator");
			}
			#print STDERR "\n1\t", join("\t", $a->qname, $a->start, $a->flag, $cigar_str), "\n" if($debug);
			print STDERR join("\t", $a->qname, $cigar_array[0]->[0], $cigar_array[0]->[1]), "\n" if($debug);
			my ($left_validated, $right_validated);
			$left_validated = $validator->validate($a, LEFT_CLIP) if($cigar_array[0]->[0] eq 'S');
			$right_validated = $validator->validate($a, RIGHT_CLIP) if($cigar_array[-1]->[0] eq 'S');
			if($cigar_array[0]->[0] eq 'S' && $left_validated){
				my $ort = '-';
				my (@sclip_lens, @sclip_sites);
				#push @sclip_sites, [$cigar_array[0]->[1], $a->start];
				push @sclip_lens, $cigar_array[0]->[1];
				push @sclip_sites, $a->start;
				# to solve very short base match problem
				if($cigar_length >= 4 && $cigar_array[0]->[1] >= $min_sc_len-3 && # S>=17 and M<=3 and N
					$cigar_array[1]->[0] eq 'M' && $cigar_array[1]->[1] <= 5 &&
				   	$cigar_array[2]->[0] eq 'N'){
					push @sclip_lens, $cigar_array[0]->[1] + $cigar_array[1]->[1];
					push @sclip_sites, $a->start + $cigar_array[1]->[1] + $cigar_array[2]->[1];
				}
			#print STDERR "2\t", join("\t", $a->qname, $a->start, $a->flag, $cigar_str, $cigar_array[0]->[1]), "\n" if($debug);

			    for(my $i=0;$i<=$#sclip_sites;$i++){ # one or two soft-clip sites

				my ($pos, $sclip_len) = ($sclip_sites[$i], $sclip_lens[$i]);
			print STDERR join("\t", $a->qname, $a->start, $a->flag, $cigar_str, $cigar_array[0]->[1], $pos, $sclip_len, "-"), "\n" if($debug);
			#print STDERR "xxxx\t", join("\t", $a->qname, $a->start, $a->flag, $cigar_str, $cigar_array[0]->[1], $pos, $sclip_len, "+"), "\n" if($debug);
#				print "($pos < $start || $pos > $end)\n";
				next if($sclip_len<$min_sc_len);
				if(exists $left_sclip_lens{$pos}) {
					$left_sclip_lens{$pos} .= "+$sclip_len";
				}
				else{
					$left_sclip_lens{$pos} = $sclip_len;
				}
				if(exists $left_cover{$pos}) {
					$left_cover{$pos} += $sclip_len;
					$left_count{$pos} += 1;
				}
				else{
					$left_cover{$pos} = $sclip_len;
					$left_count{$pos} = 1;
				}

				if($a->strand>0){
					if(exists $plus_left_count{$pos}) {
						$plus_left_count{$pos} += 1;
					}
					else{
						$plus_left_count{$pos} = 1;
					}
				}
				else{
					if(exists $neg_left_count{$pos}) {
						$neg_left_count{$pos} += 1;
					}
					else{
						$neg_left_count{$pos} = 1;
					}
				}

				next if($sclip_len<20);
				if(exists $long_left_cover{$pos}) { # the number of reads with sc length longer than 20bp
					$long_left_cover{$pos} += $sclip_len;
				}
				else{
					$long_left_cover{$pos} = $sclip_len;
				}
			    } # end for
			} # end if

			if($cigar_array[-1]->[0] eq 'S' &&  $right_validated){
				my (@sclip_lens, @sclip_sites);
				my $ort = '+';
				push @sclip_lens, $cigar_array[-1]->[1];
				push @sclip_sites, $a->end;
				#my $read_len = length($a->query->dna);
			print STDERR "3\t", join("\t", $a->qname, $a->start, $a->flag, $cigar_str), "\n" if($debug);
				# to solve very short base match problem
				if($cigar_length >=4 && $cigar_array[-1]->[1] >= $min_sc_len-3 && 
					$cigar_array[-2]->[0] eq 'M' && $cigar_array[-2]->[1] <= 5
					&& $cigar_array[-3]->[0] eq 'N'){
					print STDERR "wrong mapping at:", join("\t", $a->seq_id, $ort, $a->start, $a->end),"\n" if($debug); 
					print STDERR "CIGAR:", join("\t", $cigar_array[-4]->[0],$cigar_array[-4]->[1], 
						 $cigar_array[-3]->[0], $cigar_array[-3]->[1],
						 $cigar_array[-2]->[0], $cigar_array[-2]->[1], 
						 $cigar_array[-1]->[0], $cigar_array[-1]->[1]),"\n" if($debug); 
					push @sclip_lens, $cigar_array[-1]->[1] + $cigar_array[-2]->[1];
					push @sclip_sites, $a->end - $cigar_array[-2]->[1] - $cigar_array[-3]->[1];
				}

			    for(my $i=0;$i<=$#sclip_sites;$i++){

				my ($pos, $sclip_len) = ($sclip_sites[$i], $sclip_lens[$i]);
			print STDERR "yyyy\t", join("\t", $a->qname, $a->start, $a->flag, $cigar_str, $cigar_array[0]->[1], $pos, $sclip_len, "-"), "\n" if($debug);
				#next if($pos < $start || $pos > $end);
				next if($sclip_len<$min_sc_len);
				if(exists $right_sclip_lens{$pos}) {
					$right_sclip_lens{$pos} .= "+$sclip_len";
				}
				else{
					$right_sclip_lens{$pos} = $sclip_len;
				}

				if(exists $right_cover{$pos}) {
					$right_cover{$pos} += $sclip_len;
					print STDERR "right_cover{$pos} += $sclip_len\n" if($debug);
					$right_count{$pos} += 1;
				}
				else{
					$right_cover{$pos} = $sclip_len;
					print STDERR "right_cover{$pos} = $sclip_len\n" if($debug);
					$right_count{$pos} = 1;
				}

				if($a->strand > 0){
					if(exists $plus_right_count{$pos}) {
						$plus_right_count{$pos} += 1;
					}
					else{
						$plus_right_count{$pos} = 1;
					}
				}
				else{
					if(exists $neg_right_count{$pos}) {
						$neg_right_count{$pos} += 1;
					}
					else{
						$neg_right_count{$pos} = 1;
					}
				}

				next if($sclip_len<20);
				if(exists $long_right_cover{$pos}) { # the number of reads with sc length longer than 20bp
					$long_right_cover{$pos} += $sclip_len;
				}
				else{
					$long_right_cover{$pos} = $sclip_len;
				}
			    } # end for
			} #end if
			if($paired) { #paired but mate is not mapped, add back
				$validator->add_validator("strand_validator"); 
			}
			return;
		}
	);

	#combine soft-clip cluster
	my (%right_added, %left_added, %right_cover_cluster, %left_cover_cluster);
	print STDERR "number of left sites: ", scalar (keys %left_cover), "\n" if($debug);
	foreach my $p (sort{ $a <=> $b} keys %left_cover){
		my $exist = 0;
		for(my $i=$p-$sc_shift; $i<=$p+$sc_shift; $i++){
			if(exists($left_cover_cluster{$i})) {$exist = 1; last;}
		}
		next if($exist);
		
		my $cover = 0; # depth of soft-clipping reads
		my $max_p = $p; # the position with the max number of soft-clip reads
		for(my $j=$p; $j<=$p+$sc_shift; $j++){
			next unless(exists($left_cover{$j}));
			print STDERR "p = $p\tmax_p = $j if(", $left_cover{$j}," > ", $left_cover{$max_p},")\n" if($debug);
			# to select the sc-site with the max coverage
			$max_p = $j if( $left_cover{$j} > $left_cover{$max_p});
			print STDERR "cover += left_cover{$j} if(exists(left_cover{$j}))\n" if($debug);
			print STDERR "$cover += ", $left_cover{$j}, "\n" if($debug);
			next unless(exists($left_cover{$j}));
			$cover += $left_cover{$j};
		}
		next if($cover < $min_sc_reads);

		#the combined site +- 10 region should have at least $min_sc_reads
		my ($cover1,$count1,$neg_count1, $plus_count1, $long_sc_reads) = (0, 0, 0, 0, 0);
		for(my $k=$max_p-$sc_shift; $k<=$max_p+$sc_shift; $k++){
			$cover1 += $left_cover{$k} if(exists($left_cover{$k}));
			$count1 += $left_count{$k} if(exists($left_count{$k}));
			$neg_count1 += $neg_left_count{$k} if(exists($neg_left_count{$k}));
			$plus_count1 += $plus_left_count{$k} if(exists($plus_left_count{$k}));
			$long_sc_reads += $long_left_cover{$k} if(exists($long_left_cover{$k}));
		}
		next if($count1 < $min_sc_reads || $long_sc_reads < 20);
		print STDERR "max_p = $max_p\tcover = $cover\tlong_sc_reads:$long_sc_reads\tcover1:$cover1\n" if($debug);
		$left_cover_cluster{$max_p} = [$count1, $cover1, $plus_count1, $neg_count1];
		print STDERR "left_cover_cluster{$max_p} = $cover1\n\n" if($debug);
	}
	
	print STDERR "number of right sites: ", scalar keys(%right_cover), "\n", join(" x ", keys(%right_cover)),"\n"  if($debug);
	foreach my $p (sort{ $a <=> $b} keys %right_cover){
		my $exist = 0;
		for(my $i=$p-$sc_shift; $i<=$p+$sc_shift; $i++){
			if(exists($right_cover_cluster{$i})) {$exist = 1; last;}
		}
		next if($exist);
		
		my $cover = 0; my $max_p = $p;
		for(my $j=$p; $j<=$p+$sc_shift; $j++){
			# to select the sc-site with the max coverage
			next unless(exists($right_cover{$j}));
			print STDERR "p = $p\tmax_p = $j if(", $right_cover{$j}," > ", $right_cover{$max_p},")\n" if($debug);
			$max_p = $j if($right_cover{$j} > $right_cover{$max_p});
			print STDERR "cover += right_cover{$j} if(exists(right_cover{$j}))\n" if($debug);
			print STDERR "$cover += ", $right_cover{$j}, "\n" if($debug);
			$cover += $right_cover{$j} if(exists($right_cover{$j}));
		}
		print STDERR "max_p = $max_p\tcover = $cover\n" if($debug);
		next if($cover < $min_sc_reads);
		 my ($cover1,$count1,$long_sc_reads, $plus_count1, $neg_count1) = (0, 0, 0, 0, 0);
		 my @tmp_sclip_lens = ();
		 for(my $k=$max_p-$sc_shift; $k<=$max_p+$sc_shift; $k++){
			$cover1 += $right_cover{$k} if(exists($right_cover{$k}));
			$count1 += $right_count{$k} if(exists($right_count{$k}));
			$plus_count1 += $plus_right_count{$k} if(exists($plus_right_count{$k}));
			$neg_count1 += $neg_right_count{$k} if(exists($neg_right_count{$k}));
			#push @tmp_sclip_lens, split("+", $right_sclip_lens{$k}) if(exists($right_sclip_lens{$k}));
			print STDERR "cover1 += right_cover{$k} if(exists(right_cover{$k}))\n" if($debug);
			print STDERR "$cover1 += ", $right_cover{$k}, "\n" if($debug && exists($right_cover{$k}));
			$long_sc_reads += $long_right_cover{$k} if(exists($long_right_cover{$k}));
		}
		next if($count1 < $min_sc_reads || $long_sc_reads < 20);
		$right_cover_cluster{$max_p} = [$count1, $cover1, $plus_count1, $neg_count1];
		print STDERR "right_cover_cluster{$max_p} = $cover1\n\n" if($debug);
	}
	
	return(\%right_cover_cluster, \%left_cover_cluster);
}
=head1 NAME

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
