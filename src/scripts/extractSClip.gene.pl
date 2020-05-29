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
use List::Util;
use List::MoreUtils qw/ uniq /;
use TdtConfig; 

#custom packages
use CiceroSCValidator qw($lowqual_cutoff $min_percent_id $min_percent_hq LEFT_CLIP RIGHT_CLIP);
use constant FQ_BASE_NUMBER => 33;

if (@ARGV < 1){
	pod2usage(1);
	exit 1;
}	

use GeneModel;
use Gene;
my ($genome, $ref_genome, $gene_model_file);
my ($read_len, $min_sc_reads, $min_sc_len, $sc_shift) = (100, 2, 10, 3);
my $gm_format = "REFFLAT";

my $rmdup = 0;
my $min_percent_id = 90;
my ($out_dir, $gene_info_file);

# input/output
my ($range, $input_bam, $sample);
my $paired = 1;
my $DNA;
my ( $help, $man, $version, $usage );
my $optionOK = GetOptions(
	'i|in|input=s'	=> \$input_bam,	
	'f|gene_info=s'	=> \$gene_info_file,
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'  => \$genome,
	's=s'			=> \$sample,
	'paired!'		=> \$paired,
	'rmdup!'		=> \$rmdup,
	'm|min_sc_reads=i'	=> \$min_sc_reads,
	'l|read_len=i'	=> \$read_len,
	'lq_cutoff=i'	=> \$lowqual_cutoff,
	'min_pct_id=i'	=> \$min_percent_id,
	'min_pct_hq=i'	=> \$min_percent_hq,
	'DNA!'		=> \$DNA,
	'r|range=s'		=> \$range,
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
chomp $genome; 
my $conf = &TdtConfig::findConfig("genome", $genome);
if ($conf){
	$conf = &TdtConfig::readConfig("genome", $genome); 	
	$ref_genome = $conf->{'FASTA'}; 	
	$gene_model_file = $conf->{'REFSEQ_REFFLAT'} unless($gene_model_file);  
}
else{
	croak "Unknown genome name: $ref_genome, $conf\n";
}
croak "You need specify the input gene model file" unless ($gene_model_file);

#my $gm = GeneModelCicero->new if($gene_model_file);
my $gm = GeneModel->new if($gene_model_file);
$gm->from_file($gene_model_file, $gm_format);

my $sam = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $ref_genome);

my @gene_info;
open my $GI, "$gene_info_file";
while(my $line = <$GI>){
   chomp($line);
   my @fields = split(/\t/, $line);
   my $tmp_gene = {
		name        => $fields[0], 
	   	range       => $fields[1], 
	   	strand      => $fields[2], 
		mRNA_length => $fields[3], 
		reads_count => $fields[4], 
		sc_cutoff   => $fields[5]
 	};
   push @gene_info, $tmp_gene;
}
close $GI;

#sub count_coverage {

my $validator = CiceroSCValidator->new();
if(!$paired) {
	$validator->remove_validator("strand_validator");
}

`mkdir -p $out_dir`;
my $cover_file = "$out_dir/local.cover";
my $SC_file = "$out_dir/local.SC";
open my $SC, ">$cover_file";
open my $SC2, ">$SC_file";
foreach my $tmp_gene (@gene_info){
	my ($gRange, $sc_cutoff, $r_cnt, $mRNA_length) = ($tmp_gene->{range}, $tmp_gene->{sc_cutoff}/2, $tmp_gene->{reads_count}, $tmp_gene->{mRNA_length});
	next if($r_cnt == 0);
	my ($chr, $start, $end) = split /[:|-]/, $gRange;

	my($pcover, $ncover) = extract_range_sclip(
		-SAM => $sam, 
		-RANGE => $gRange, 
		-VALIDATOR => $validator);

	my $ratio = -1;
	foreach my $p (keys(%{$pcover})) {
		#my $c = ($rmdup ? scalar(@{$pcover->{$p}}) : $pcover->{$p});
		my @c = @{$pcover->{$p}};
		my $tc = count_coverage($sam, $chr, $p);
		print $SC join("\t", $chr, $p, "+", @c), "\t", $tc->{"+"},"\t", $tc->{"-"}, "\n" if($c[0]>=$min_sc_reads);
		if($DNA){
			my $p_ratio = $c[2]/($tc->{"+"} + 0.001);
			my $n_ratio = $c[3]/($tc->{"-"} + 0.001);
			$ratio =  ($p_ratio > $n_ratio) ? $p_ratio : $n_ratio;
		}
		else{
			$ratio = ($c[1]/($read_len-20))/($r_cnt*$read_len/(4*$mRNA_length));
		}
		print $SC2 join("\t", $chr, $p, "+", $c[0], sprintf("%.2f",$sc_cutoff), sprintf("%.3f",$ratio), $c[1]), "\n";
	}
	foreach my $p (keys(%{$ncover})) {
		#my $c = ($rmdup ? scalar(@{$ncover->{$p}}) : $ncover->{$p});
		my @c = @{$ncover->{$p}};
		my $tc = count_coverage($sam, $chr, $p);
		print $SC join("\t", $chr, $p, "-", @c), "\t", $tc->{"+"},"\t", $tc->{"-"}, "\n" if($c[0]>=$min_sc_reads);
		if($DNA){
			my $p_ratio = $c[2]/($tc->{"+"} + 0.001);
			my $n_ratio = $c[3]/($tc->{"-"} + 0.001);
			$ratio =  ($p_ratio > $n_ratio) ? $p_ratio : $n_ratio;
		}
		else{
			$ratio = ($c[1]/($read_len-20))/($r_cnt*$read_len/(4*$mRNA_length));
		}
		print $SC2 join("\t", $chr, $p, "-", $c[0], sprintf("%.2f",$sc_cutoff), sprintf("%.3f", $ratio), $c[1]), "\n";
	}
}
close $SC;
close $SC2;

exit(0);

sub low_complexity{
	my $sequence = shift;
        my $max_single_nt = 0.8 * length($sequence);
        my $max_run_nt = 0.6 * length($sequence);

        return 1 if @{[$sequence =~ /(A)/g]} > $max_single_nt;
        return 1 if @{[$sequence =~ /(C)/g]} > $max_single_nt;
        return 1 if @{[$sequence =~ /(T)/g]} > $max_single_nt;
        return 1 if @{[$sequence =~ /(G)/g]} > $max_single_nt;
        #return 1 if $sequence =~ /(A{$max_run_nt,})/;

	my $mask_seq = $sequence;
	$mask_seq =~ s/((.+)\2{3,})/'N' x length $1/eg;
        return 1 if @{[$sequence =~ /(N)/g]} > $max_single_nt;
	return 1 if $mask_seq =~ /(N{$max_run_nt,})/;
	return 0;
}

sub calc_sc_cutoff{
	my ($sam, $range, $mRNA_length) = @_;
	my $cnt=0;
	$sam->fetch($range, 
		sub {
			my $a = shift;
			return if( ($a->flag & 0x0400) || ($a->flag & 0x0004) ); #PCR duplicate or unmapped
			$cnt++;
		});
#		print "gene: $g->name, gRange: $gRange\ncnt: $cnt\nmRNA_length: $mRNA_length\n";
		my $RPK_value = 1000*$cnt/$mRNA_length;
		my $sc_cutoff = ($read_len-20)*$cnt/(100*$mRNA_length);
		print join("\t", $range, $mRNA_length, $cnt, $RPK_value, $sc_cutoff), "\n";
		return $sc_cutoff;
		#return [$gRange, $mRNA_length, $cnt, $RPK_value, $sc_cutoff];
}


sub count_coverage {
	my ($sam, $chr, $pos) = @_;
	my $seg = $sam->segment(-seq_id => $chr, -start => $pos, -end => $pos);
	#print STDERR "seg: ", $seg, "\n";
	return 0 unless $seg;
	my ($pos_N, $neg_N) = (0, 0);
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
			#return unless(abs($a->start - 28608251) < 20);
			#return unless($a->qname =~ m/C0FVUACXX120313:3:2302:5930:71041/);
			#return unless($a->qname =~ m/16389/);
			if($paired && !$a->proper_pair){ #paired but mate is not mapped 
				$validator->remove_validator("strand_validator");
			}
			#print STDERR "\n1\t", join("\t", $a->qname, $a->start, $a->flag, $cigar_str), "\n" if($debug);
			my ($left_validated, $right_validated);
			$left_validated = $validator->validate($a, LEFT_CLIP) if($cigar_array[0]->[0] eq 'S');
			print STDERR join("\t", $a->qname, $a->start, $cigar_array[0]->[0], $cigar_array[0]->[1], $left_validated), "\n" if($debug);
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

