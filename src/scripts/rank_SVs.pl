#!/usr/bin/env perl 

use strict;
use warnings; 

use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Path;
use Cwd qw[abs_path getcwd];
use List::Util qw[min max];
use TdtConfig; 

use DelimitedFile;
use File::Temp qw/ tempdir /;

my ($annotated_file, $out_file);
my ($min_reads_cnt, $max_repeat_score, $min_ratio) = (2,0.7,0.01);
my $read_len;
my $genome;# = "GRCh37-lite";
my $optionOK = GetOptions(
	'i|input=s'	=> \$annotated_file,
	'o|out_file=s'	=> \$out_file,
	'l|read_len=s'	=> \$read_len,
	'mr|min_ratio=s'	=> \$min_ratio,
	'genome=s'  => \$genome,
);

my $min_match_len = $read_len*0.5;
my $min_coverage = $read_len*0.5;

my $out_dir = dirname($annotated_file);
$out_file = "$out_dir/final_fusions.txt" unless($out_file);

my ($known_fusion_file, $known_itd_file);
my $conf; 
if (&TdtConfig::findConfig("genome", $genome)){
	$conf = &TdtConfig::readConfig("genome", $genome); 
}
else{
	croak("no config");
}

$known_fusion_file = $conf->{KNOWN_FUSIONS} unless($known_fusion_file);
print STDERR "known_fusion_file: $known_fusion_file\n";
$known_itd_file = $conf->{KNOWN_ITD_FILE} unless($known_itd_file);

my %known_fusions=();
my %known_fusion_partners=();

my $df = new DelimitedFile(
	       "-file" => $annotated_file,
	       "-headers" => 1,
	      );

open(my $KF, "$known_fusion_file");
while(<$KF>){

	chomp;
	my ($fg1, $fg2) = split(/\t/, $_);
	if($fg1 gt $fg2){
		my $tmp = $fg1;
		$fg1 = $fg2;
		$fg2 = $tmp;
	}
	my $fusion = $fg1.":".$fg2;
	next if(exists($known_fusions{$fusion}));
	$known_fusions{$fusion} = 1;
	if(exists($known_fusion_partners{$fg1})){
		$known_fusion_partners{$fg1} ++;
	}
	else {
		$known_fusion_partners{$fg1} = 1;
	}

	if(exists($known_fusion_partners{$fg2})){
		$known_fusion_partners{$fg2} ++;
	}
	else {
		$known_fusion_partners{$fg2} = 1;
	}
}
close($KF);

my %known_ITDs = ();
open(my $ITD_F, $known_itd_file);
while(<$ITD_F>){
	chomp;
	my ($gene, $chr, $start, $end) = split(/\t/,$_);
	$known_ITDs{$gene} = [$start, $end];
}
close($ITD_F);

my %all_svs = ();
my %sv_seqs = ();
my @uniq_SVs;
my $N_col = 0; # number of columns
while (my $row = $df->get_hash()) {

	next if($row->{readsA} eq 'readsA');
	my $first_bp = {
		reads_num => $row->{readsA},
		gene => $row->{geneA},
		tname => $row->{chrA},
		tpos => $row->{posA},
		qstrand => $row->{ortA},
		feature => $row->{featureA},
		matches => $row->{matchA},
		repeat => $row->{repeatA},
		area => $row->{coverageA},
		expression_ratio => $row->{ratioA1},
		total_reads => $row->{total_readsA},
		qpos => $row->{qposA},
		maf => $row->{ratioA},
		sv_refseq => $row->{sv_refseqA},
		sv_refseq_codon => $row->{sv_refseqA_codon},
		sv_refseq_exon => $row->{sv_refseqA_exon},
		sv_refseq_anchor_type => $row->{sv_refseqA_anchor_type},
		sv_refseq_coding_base_number => $row->{sv_refseqA_coding_base_number},
		sv_refseq_last_coding_base_number => $row->{sv_refseqA_last_coding_base_number},
		sv_refseq_AA_index => $row->{sv_refseqA_AA_index},
		sv_refseq_contig_index => $row->{sv_refseqA_contig_index},
	};

	foreach my $k (keys %{$first_bp}) {
		$first_bp->{$k} = '' unless(defined($first_bp->{$k}));
	}

	
	my $second_bp = {
		reads_num => $row->{readsB},
		gene => $row->{geneB},
		tname => $row->{chrB},
		tpos => $row->{posB},
		qstrand => $row->{ortB},
		feature => $row->{featureB},
		matches => $row->{matchB},
		repeat => $row->{repeatB},
		area => $row->{coverageB},
		expression_ratio => $row->{ratioB1},
		total_reads => $row->{total_readsB},
		qpos => $row->{qposB},
		maf => $row->{ratioB},
		sv_refseq => $row->{sv_refseqB},
		sv_refseq_codon => $row->{sv_refseqB_codon},
		sv_refseq_exon => $row->{sv_refseqB_exon},
		sv_refseq_anchor_type => $row->{sv_refseqB_anchor_type},
		sv_refseq_coding_base_number => $row->{sv_refseqB_coding_base_number},
		sv_refseq_last_coding_base_number => $row->{sv_refseqB_last_coding_base_number},
		sv_refseq_AA_index => $row->{sv_refseqB_AA_index},
		sv_refseq_contig_index => $row->{sv_refseqB_contig_index},
	};

	if($row->{sv_posB_adjusted} =~ m/,/){
		my @sv_posB_adjusted = split(/,/, $row->{sv_posB_adjusted});
		$second_bp->{tpos} = $sv_posB_adjusted[0] if(@sv_posB_adjusted && $sv_posB_adjusted[0] =~ m/^\d+$/);
	}
=pod
	my $second_bp = {
		reads_num => $fields[13],
		gene => $gene2,
		tname => $fields[7],
		tpos => $fields[8],
		qstrand => $fields[9],
		feature => $fields[10],
		matches => $fields[15],
		repeat => $fields[17],
		area => $fields[19],
		maf => $fields[21],
		expression_ratio => $fields[25],
		total_reads => $fields[27],
		qpos => $fields[23],
	};
=cut
	my ($sample, $sv_ort, $qseq, $type, $inframe) = ($row->{sample}, $row->{sv_ort}, $row->{contig}, $row->{type}, $row->{sv_inframe});

	my $tmp_SV = {
		junc_seq => $qseq,
		first_bp => $first_bp,
		second_bp => $second_bp,
		ort => $sv_ort,
		type => $type,
		sample => $sample,
		frame => $inframe,
		sv_AA => $row->{sv_AA},
		sv_desc => $row->{sv_desc},
		sv_processing_exception => $row->{sv_processing_exception},
		sv_general_info => $row->{sv_general_info},
		sv_interstitial_AA => $row->{sv_interstitial_AA},
		sv_frame_index => $row->{sv_frame_index}
		};

	push @uniq_SVs, scoring($tmp_SV) unless(is_dup_SV(\@uniq_SVs, $tmp_SV));
}

my @sorted_SVs = sort{   $a->{rating} cmp $b->{rating} || 
			$b->{score} <=> $a->{score}
		}@uniq_SVs;
print STDERR "number of SVs: ", scalar @uniq_SVs, "\n";

my $out_header = join("\t", "sample", "geneA", "chrA", "posA", "ortA", "featureA", "geneB", "chrB", "posB", "ortB", 
			"featureB", "sv_ort", "readsA", "readsB", "matchA", "matchB", "repeatA", "repeatB","coverageA", 
			"coverageB", "ratioA", "ratioB", "qposA", "qposB", "total_readsA", "total_readsB", "contig", "type");
$out_header .="\tscore\trating\tmedal\tfunctional effect";
$out_header = join("\t", $out_header, "frame", "sv_refseqA", "sv_refseqA_codon", "sv_refseqA_exon", "sv_refseqA_anchor_type", "sv_refseqA_coding_base_number", "sv_refseqA_last_coding_base_number", "sv_refseqA_AA_index", "sv_refseqA_contig_index"); 
$out_header = join("\t", $out_header, "sv_refseqB", "sv_refseqB_codon", "sv_refseqB_exon", "sv_refseqB_anchor_type", "sv_refseqB_coding_base_number", "sv_refseqB_last_coding_base_number", "sv_refseqB_AA_index", "sv_refseqB_contig_index"); 
$out_header = join("\t", $out_header, "sv_AA", "sv_desc", "sv_processing_exception", "sv_general_info", "sv_interstitial_AA","sv_frame_index"); 
my $k=0;
my @HQ_SVs;
my %HQ_genes=();
foreach my $sv (@sorted_SVs){

	my ($sample, $bp1, $bp2, $qseq, $type, $score) = 
	($sv->{sample}, $sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq}, $sv->{type}, $sv->{score});

	if($sv->{rating} eq 'HQ'){
		push @HQ_SVs, $sv; 
		$k++;
		next;
	}

	if($sv->{type} eq 'read_through'){
	 	push @HQ_SVs, $sv if($sv->{rating} eq 'RT');
		next;
	}

	last if($sv->{rating} eq 'bad');
	if ($k <= 10 && $bp1->{reads_num} > $min_reads_cnt && $bp2->{reads_num} > $min_reads_cnt){
			$sv->{rating} = 'HQ';
	}
	elsif($k < 10 || ($bp1->{reads_num} <= $min_reads_cnt || $bp2->{reads_num} <= $min_reads_cnt)){
		$sv->{rating} = 'HQ' if($sv->{frame} && $sv->{frame} =~ /1|2/);
	}
	else{
		$sv->{rating} = 'LQ';
	}

	if($sv->{medal} > 1 && $sv->{frame} && $sv->{frame} =~ /1|2/){
		$sv->{rating} = 'HQ' if($sv->{rating} eq 'LQ');
		$sv->{rating} = 'LQ' if($sv->{rating} eq 'bad' && $sv->{type} ne 'read_through');
	}
	push @HQ_SVs, $sv;
	$k++;
}
print STDERR "number of remaining SVs: ", scalar @HQ_SVs, "\n";

my @final_SVs = sort{   $a->{rating} cmp $b->{rating} || 
			$b->{score} <=> $a->{score}
		}@HQ_SVs;

open(my $hFo, ">$out_file");
print $hFo $out_header, "\n";
foreach my $sv (@final_SVs){

	my ($sample, $bp1, $bp2, $qseq, $type, $score) = 
	($sv->{sample}, $sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq}, $sv->{type}, $sv->{score});

	my $type2 = 'other';

	if($sv->{type} =~ /Internal_dup/){
		$type2 = 'ITD';
	}
	elsif($sv->{frame} =~ /1/){
		$type2 = 'Fusion';
	}
	elsif($sv->{frame} =~ /2/){
		$type2 = 'upTSS';
	}
	else{
		$type2 = 'NLoss' if ($bp1->{feature} ne 'coding' && $bp2->{feature} eq 'coding');
		$type2 = 'CLoss' if ($bp1->{feature} eq 'coding' && $bp2->{feature} ne 'coding');
		$type2 = 'Fusion' if ($bp1->{feature} eq 'coding' && $bp2->{feature} eq 'coding');
	}

	my $out_string = join("\t", $sample, $bp1->{gene}, $bp1->{tname}, $bp1->{tpos}, $bp1->{qstrand}, $bp1->{feature}, 
	   $bp2->{gene}, $bp2->{tname}, $bp2->{tpos}, $bp2->{qstrand}, $bp2->{feature}, $sv->{ort}, $bp1->{reads_num}, $bp2->{reads_num}, 
	   $bp1->{matches}, $bp2->{matches}, sprintf("%.3f", $bp1->{repeat}), sprintf("%.3f", $bp2->{repeat}), $bp1->{area}, $bp2->{area},
	   sprintf("%.3f",$bp1->{maf}), sprintf("%.3f", $bp2->{maf}), $bp1->{qpos}, $bp2->{qpos}, $bp1->{total_reads}, $bp2->{total_reads}, $qseq, $type);

	$out_string .= "\t".join("\t", sprintf("%.2f", $sv->{score}), $sv->{rating}, $sv->{medal}, $type2);
	$out_string = join("\t", $out_string, $sv->{frame}, $bp1->{sv_refseq}, $bp1->{sv_refseq_codon},$bp1->{sv_refseq_exon},$bp1->{sv_refseq_anchor_type},
			$bp1->{sv_refseq_coding_base_number},$bp1->{sv_refseq_last_coding_base_number}, $bp1->{sv_refseq_AA_index}, $bp1->{sv_refseq_contig_index}); 
	$out_string = join("\t", $out_string, $bp2->{sv_refseq}, $bp2->{sv_refseq_codon},$bp2->{sv_refseq_exon},$bp2->{sv_refseq_anchor_type},
			$bp2->{sv_refseq_coding_base_number},$bp2->{sv_refseq_last_coding_base_number}, $bp2->{sv_refseq_AA_index}, $bp2->{sv_refseq_contig_index}); 
	$out_string = join("\t", $out_string, $sv->{sv_AA}, $sv->{sv_desc}, $sv->{sv_processing_exception}, $sv->{sv_general_info}, $sv->{sv_interstitial_AA}, $sv->{sv_frame_index}); 
	print $hFo "$out_string\n";
}
close($hFo);

sub scoring {
	my $sv = shift;
	my ($sample, $bp1, $bp2, $qseq, $type) = 
	($sv->{sample}, $sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq}, $sv->{type});

	my($fg1, $fg2) = ($bp1->{gene}, $bp2->{gene});
	my $fusion = ($fg1 lt $fg2) ? "$fg1:$fg2" : "$fg2:$fg1";
	my $debug = 0;
	my $medal = 0;
	my $rating = 'LQ';
	$medal = 1 if(exists($known_fusion_partners{$fg1}));
	$medal = 1 if(exists($known_fusion_partners{$fg2}) && $sv->{ort} eq "?");
	$medal = 2 if(exists($known_fusion_partners{$fg2}) && $sv->{ort} eq ">");
	$medal = 3 if(exists($known_fusion_partners{$fg1}) && exists($known_fusion_partners{$fg2}) && $sv->{type} !~ /Internal/);
	$medal = 4 if(exists($known_fusions{$fusion}));
	
	if(exists($known_ITDs{$fg1}) && $type eq 'Internal_dup' &&
	   ($bp1->{tpos} - $known_ITDs{$fg1}[0]) * ($bp1->{tpos} < $known_ITDs{$fg1}[1]) > 0 &&
	   ($bp2->{tpos} - $known_ITDs{$fg1}[0]) * ($bp2->{tpos} < $known_ITDs{$fg1}[1]) > 0){
			$rating = 'HQ'; $medal = 4;
	}

	my ($matchA, $matchB, $ratioA, $ratioB) = ($bp1->{matches}, $bp2->{matches}, $bp1->{maf}, $bp2->{maf});

	$ratioA = ($ratioA < $min_ratio) ? exp($min_ratio/$ratioA) : 1 if($ratioA > 0);
	$ratioB = ($ratioB < $min_ratio) ? exp($min_ratio/$ratioB)  : 1 if($ratioB > 0);


	$matchA = ($matchA < $min_match_len) ? exp(0.5*($matchA - $min_match_len)) : 1;
	$matchB = ($matchB < $min_match_len) ? exp(0.5*($matchB - $min_match_len)) : 1;

	my $featureA = 0.5; $featureA = 0.8 if($bp1->{feature} =~ m/intron/);
	$featureA = 0.9 if($bp1->{feature} =~ m/utr/); $featureA = 1 if($bp1->{feature} =~ m/coding/);
	my $featureB = 0.5; $featureB = 0.8 if($bp2->{feature} =~ m/intron/);
	$featureB = 0.9 if($bp2->{feature} =~ m/utr/); $featureB = 1 if($bp2->{feature} =~ m/coding/);
	my $ort = 1; 
	$ort = 2 if($sv->{ort} eq ">");
	my $frame = 1; 
	$frame = 2 if($sv->{frame} && $sv->{frame} =~ m/1|2/);

	my $scoreA = $bp1->{area}*$ratioA*$matchB*(1-$bp2->{repeat});
	my $scoreB = $bp2->{area}*$ratioB*$matchA*(1-$bp1->{repeat});
	my $score = 0.5*($scoreA+$scoreB)*$ort*$frame;
	$score *= 0.1 if( $frame < 2 && $type eq 'DEL' || $type eq 'read_through');
	$score *= 0.1 if( $frame < 2 && $type eq 'Internal_dup');
	$score *= 0.1 if(!exists($known_ITDs{$fg1}) && $type eq 'Internal_dup');
	$score *= 0.1 if($bp1->{tname} eq $bp2->{tname} && abs($bp1->{tpos} - $bp2->{tpos}) < 40000 && $bp1->{qstrand} ne $bp2->{qstrand});
	$score *= 0.1 if($bp1->{tname} eq $bp2->{tname} && abs($bp1->{tpos} - $bp2->{tpos}) < 10000 && $bp1->{qstrand} ne $bp2->{qstrand});

	$rating = 'RT'  if($sv->{type} eq 'read_through');
	$rating = 'bad' if($bp1->{reads_num} + $bp2->{reads_num} <= $min_reads_cnt);
	$rating = 'bad' if($bp1->{repeat} > $max_repeat_score &&  $bp2->{repeat} > $max_repeat_score);
	$rating = 'bad' if($bp1->{maf} < $min_ratio && $bp2->{maf} < $min_ratio);
	$rating = 'bad' if($score < 1);
	$rating = 'HQ'  if($medal == 4);
	$sv->{score} = $score;
	$sv->{rating} = $rating;
	$sv->{medal} = $medal;
	print STDERR "fusion is $fusion\tscore = $score\trating is $rating\tmedal is $medal\n" if($debug);
	return $sv;
}	

sub is_dup_SV {
	my($r_SVs, $sv) = @_;
	foreach my $s (@{$r_SVs}) {
		my $more_reads = ($s->{first_bp}->{reads_num} + $s->{second_bp}->{reads_num} > $sv->{first_bp}->{reads_num} + $sv->{second_bp}->{reads_num}) ? 1 : 0;
		return 1
		if( 	$more_reads &&
			abs($s->{first_bp}->{tpos} - $sv->{first_bp}->{tpos}) < 10 &&
			abs($s->{second_bp}->{tpos} - $sv->{second_bp}->{tpos}) < 10 &&
			$s->{first_bp}->{tname} eq $sv->{first_bp}->{tname} &&
			$s->{second_bp}->{tname} eq $sv->{second_bp}->{tname});
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
