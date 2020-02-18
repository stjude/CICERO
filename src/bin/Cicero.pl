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
use Bio::SearchIO;
use Bio::SeqIO;
use File::Temp qw/ tempfile tempdir /;
use File::Spec;
use File::Path;
use File::Copy;
use File::Basename;
use List::MoreUtils qw/ uniq /;
use Cwd;

use TdtConfig;

use CiceroSCValidator qw($min_percent_hq LEFT_CLIP RIGHT_CLIP);
use CiceroUtil qw(print_out prepare_reads_file parse_range rev_comp is_PCR_dup 
	          read_fa_file get_discordant_reads get_sclip_reads);
require CiceroExtTools;

use Transcript;
use Gene;
use GeneModel;

my $debug = 0;
# input/output
my ($out_dir, $input_bam, $sclip_file);
#my $out_suffix = "predSV.txt";
my ($sample, $genome, $ref_genome, $gene_model_file, $gm);
my $gm_format = "REFFLAT";

# external programs related variable
# 1. cap3 related variables
my $cap3_options = " -o 25 -z 2 -h 60 -y 10 > /dev/null 2>&1";

# 2 blat related variables, using blat server and blat standalone
my $blat_client_exe = "gfClient";
my $blat_client_options = '-out=psl -nohead > /dev/null 2>&1';
my ($blat_server, $blat_port, $dir_2bit); 
my ($paired, $rmtmp, $rmdup) = (1, 1, 1);
# other options 
my ($read_len, $min_fusion_distance);
my ($min_sclip_len, $min_hit_len) = (5, 25);
my ($min_sclip_reads, $expression_ratio_cutoff, $max_num_hits) = (2, 0.01, 3);

if(@ARGV == 0){
	#TODO: get the correct usage string
	print STDERR "Usage: $0 -g <genome> -i <bam_file> -o <out_dir> -f <sclip_file>\n";
	exit 1; 
}

# common help options
my ( $help, $man, $version, $usage );
my $optionOK = GetOptions(
	# input/output 
	's|sample=s'	=> \$sample,
	'i|input_bam=s'	=> \$input_bam,
	'ratio=f'	=> \$expression_ratio_cutoff,
	'f|sclipfile=s'	=> \$sclip_file,
	'o|out_dir=s'	=> \$out_dir,
	'l|read_len=i'	=> \$read_len,

	'genome=s'	=> \$genome,
	# external programs location and options
	'cap3opt=s'		=> \$cap3_options,
	'blatclientopt=s'	=> \$blat_client_options,
	'blatserver=s'	=> \$blat_server,
	'blatport=i'	=> \$blat_port,
	'2bitdir=s'		=> \$dir_2bit,
	#RNAseq support
	'genemodel=s'		=> \$gene_model_file,
	'gmformat=s'		=> \$gm_format,
	#other related parameters
	'min_fs_dist=i'	=> \$min_fusion_distance,
	'm|min_sclip_reads=i'		=> \$min_sclip_reads,
	'min_sclip_len=i'	=> \$min_sclip_len,
	'min_hit_len=i'		=> \$min_hit_len,
	'max_num_hits=i'	=> \$max_num_hits,
	'paired!'		=> \$paired,
	'rmtmp=i'		=> \$rmtmp,
	'rmdup=i'		=> \$rmdup,

	# common help parameters
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'version'		=> \$version,
	'debug'			=> \$debug,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );

# figure out input file
if(!$input_bam) {
	croak "You need specify input bam file(s)";
}

if(!$sclip_file) {
	croak "You need to specify the softclipping file";
}
$sclip_file = File::Spec->rel2abs($sclip_file);

my $conf = &TdtConfig::readConfig("genome", $genome);
$ref_genome = $conf->{'FASTA'} unless ($ref_genome && -e $ref_genome);  
$blat_server = $conf->{'BLAT_HOST'} unless ($blat_server); 
$blat_port = $conf->{'BLAT_PORT'} unless ($blat_port); 
$gene_model_file = $conf->{'REFSEQ_REFFLAT'} unless ($gene_model_file && -e $gene_model_file); 
$dir_2bit = '/';

if(!$ref_genome) {
	croak "You need to specify the reference genome";
}

$conf = &TdtConfig::readConfig('app', 'cicero'); 
$min_hit_len = $conf->{MIN_HIT_LEN} unless($min_hit_len);
$max_num_hits = $conf->{MAX_NUM_HITS} unless($max_num_hits);
$min_fusion_distance = $conf->{MIN_FUSION_DIST} unless($min_fusion_distance);
$min_sclip_reads = $conf->{MIN_SC_READS} unless($min_sclip_reads);
$min_sclip_len = $conf->{MIN_SC_LEN} unless($min_sclip_len);

croak "You need specify the input gene model file" unless ($gene_model_file);
if($gene_model_file) {
	$gm = GeneModel->new if($gene_model_file);
	$gm->from_file($gene_model_file, $gm_format);
}
# set up the external programs and validators
# Those variable will be global
my $assembler = Assembler->new( 
	-PRG => 'cap3',
	-OPTIONS => $cap3_options
);

my $mapper = Mapper->new(
	-PRG => join(' ', ($blat_client_exe, $blat_server, $blat_port)),
	-OPTIONS => $blat_client_options,
	-BIT2_DIR => $dir_2bit,
	-MIN_HIT_LEN => $min_hit_len,
	-MAX_NUM_HITS => $max_num_hits,
	-MIN_FS_DIST => $min_fusion_distance,
);

my $input_base;
$input_bam = File::Spec->rel2abs($input_bam);
$sample = basename($input_bam, ".bam") unless($sample);

my $validator = CiceroSCValidator->new();
$validator->remove_validator('strand_validator') if(!$paired);

#setup output and working directory
$out_dir = getcwd if(!$out_dir);
mkdir $out_dir if(!-e $out_dir || ! -d $out_dir);
my $tmp_dir = "$out_dir/tmp"; 
mkdir $tmp_dir if(!-e $tmp_dir || ! -d $tmp_dir);  

my $sam_d = Bio::DB::Sam->new( -bam => $input_bam, -fasta => $ref_genome);

system("touch $out_dir/unfiltered.fusion.txt"); 
system("touch $out_dir/unfiltered.internal.txt"); 

# the softclip file is sorted, so no need to re-sort it
open my $SCLIP, "<$sclip_file" or croak "can't open $sclip_file:$OS_ERROR";
while( my $line = <$SCLIP> ) {
	chomp $line;
	my ($chr, $pos, $ort, $sc_cnt, $sc_cutoff, $cover, $expression_ratio) = split /\t/, $line;
	next if($sc_cnt < $min_sclip_reads);
	if($sc_cutoff > 0){ 
		#print STDERR "next if($expression_ratio < $expression_ratio_cutoff)\n";
		#next if(defined $expression_ratio && $expression_ratio < $expression_ratio_cutoff);
		#next if($sc_cnt<$sc_cutoff || $expression_ratio < $expression_ratio_cutoff);
	}else{
		my $cnt = count_coverage($sam_d, $chr, $pos);
		next if(!$cnt);
		$sc_cutoff = (-1)*$cnt/20;
		next unless(($cnt < 2000 && $sc_cnt > abs($sc_cutoff)/10) || $sc_cnt > abs($sc_cutoff));
	}
	my $clip = ($ort eq '-')? LEFT_CLIP: RIGHT_CLIP;
	print STDERR "detect_SV($sam_d, $chr, $pos, $sc_cnt, $sc_cutoff, $clip)\n" if($debug);
	detect_SV($sam_d, $chr, $pos, $sc_cnt, $sc_cutoff, $clip);
}
close $SCLIP;

sub count_coverage {
	my ($sam, $chr, $pos) = @_;
	my $seg = $sam->segment(-seq_id => $chr, -start => $pos, -end => $pos);
	return 0 unless $seg;
	my $n = 0;
	my $itr = $seg->features(-iterator => 1);
	while( my $a = $itr->next_seq) {
		next unless($a->start && $a->end);
		next if(abs($a->start - $a->end) > 2*$read_len);
		$n++;
	}
	return $n;
}

sub detect_SV{
	my ($sam_d, $chr, $sc_site, $sc_cover, $sc_cutoff, $clip) = @_;
	my $debug = 0;
	print "\nsam_d: $sam_d, out_dir:$out_dir\nchr: $chr, sc_site: $sc_site, sc_cover: $sc_cover, sc_cutoff: $sc_cutoff, clip: $clip\n" if($debug);
	my ($tid) = $sam_d->header->parse_region($chr);
	my $chr_len = $sam_d->header->target_len->[$tid];

	my $unmapped_cutoff = 6*$sc_cutoff;
	my $fixSC = 0;

	# to fix the inaccurate soft-clipping
	$fixSC = 1 if($sc_cover < 10 && abs($sc_cutoff) < 10); 
	my $fa_file = "$tmp_dir/". join(".", $chr, $sc_site, ($clip+1), "fa");
	my $now_string = localtime;  # e.g., "Thu Oct 13 04:54:34 1994"

	print STDERR "prepare_reads_file:", join("\n", $fa_file, $sam_d, $chr, $sc_site, $clip, $validator, $paired, $rmdup, $min_sclip_reads, $min_sclip_len, $unmapped_cutoff),"\n" if($debug);
	prepare_reads_file(-OUT => $fa_file,
		           -SAM => $sam_d,
			   -CHR =>$chr, 
			   -POS => $sc_site, 
		   	-CLIP => $clip, 
		   	-VALIDATOR => $validator,
		   	-PAIRED => $paired,
		   	-RMDUP => $rmdup,
		   	-MIN_SC => $min_sclip_reads,
			-MIN_SC_LEN => $min_sclip_len,
			-UNMAPPED_CUTOFF => $unmapped_cutoff,
			-FIXSC => $fixSC,
	        	);

	return unless(-f $fa_file && -s $fa_file);
	print STDERR "start to assemble reads: at $fa_file", "\n" if($debug);
	my($contig_file, $sclip_count, $contig_reads) = $assembler->run($fa_file); 
	if(not -s $contig_file){ 
		unlink glob "$fa_file*" if($rmtmp && $fa_file =~ /fa/);
		return; }

	print STDERR "finished assembly and start to map contigs...", "\n" if($debug);
	my @mappings;
	print STDERR "start mapping ... $contig_file\nsc_site: $sc_site\tclip: $clip\tmin_hit_len: $min_hit_len\n" if($debug && -s $contig_file);
	my $ref_chr = $chr; $ref_chr =~ s/chr//;
	push @mappings, $mapper->run(-QUERY=>$contig_file, -scChr => $ref_chr, -scSite=>$sc_site, -CLIP=>$clip, -READ_LEN=>$read_len) if(-s $contig_file);

	my %num_of_mappings = ();
	foreach my $sv (@mappings){
		# count the number of mappings for each contig
		my $ctg = $sv->{contig_name};
		if(exists($num_of_mappings{$ctg})) {$num_of_mappings{$ctg}++;}
		else{$num_of_mappings{$ctg} = 1;}
	}

	print STDERR "rm $fa_file*\n" if($debug && $rmtmp && $fa_file =~ /fa/);
	unlink glob "$fa_file*" if($rmtmp && $fa_file =~ /fa/);
	my $n_m = scalar @mappings;
	print STDERR "\nnumber of mappings: $n_m\tmin_hit_len: $min_hit_len\nmappings: @mappings\n" if($debug);

	my $unfiltered_fusion_file = "$out_dir/unfiltered.fusion.txt";
	my $unfiltered_internal_events_file = "$out_dir/unfiltered.internal.txt";
	open(my $UFF, ">>$unfiltered_fusion_file");
	open(my $UIF, ">>$unfiltered_internal_events_file");

		# If we make it this far, create the output files before exiting. 
		# Nothing has failed at this point, there are simply no results.
		return if($n_m == 0);
		foreach my $sv (@mappings){
			
			my ($bp1, $bp2, $qseq, $qname) = ($sv->{first_bp}, $sv->{second_bp}, $sv->{junc_seq}, $sv->{contig_name});
			if($chr =~ m/chr/ && $bp1->{tname} !~ m/chr/) {$bp1->{tname} = "chr".$bp1->{tname};}
			if($chr =~ m/chr/ && $bp2->{tname} !~ m/chr/) {$bp2->{tname} = "chr".$bp2->{tname};}

			$sv->{gap} = $bp1->{ort} > 0 ? $bp2->{qstart} - $bp1->{qend} : $bp1->{qstart} - $bp2->{qend};
			print STDERR "bp1: ", join("\t",$bp1->{tname}, $bp1->{qstrand}, $bp1->{tend}, $bp1->{tstart}, $bp1->{ort}, $bp1->{qend}, $bp1->{qstart}), "\n" if($debug);
			print STDERR "bp2: ", join("\t",$bp2->{tname}, $bp2->{qstrand}, $bp2->{tend}, $bp2->{tstart}, $bp2->{ort}, $bp2->{qend}, $bp2->{qstart}), "\n" if($debug);
	
=pos
			if($num_of_mappings{$qname} >= $max_num_hits){
				my $reads_num2 = second_sc_chk($sam_d, $bp2);
				print STDERR "reads_num2: $reads_num2\tn_m: $n_m\n" if($debug);
				next if($reads_num2 == 0);
			}
=cut
			my ($pos1, $pos2);
			$pos1 = ($bp1->{qstrand}*$bp1->{ort} > 0)? $bp1->{tend}: $bp1->{tstart};
			$bp1->{tpos} = $pos1;	
				
			$pos2 = ($bp2->{qstrand}*$bp2->{ort} > 0) ? $bp2->{tend}: $bp2->{tstart};
			$bp2->{tpos} = $pos2;	

			my $g1_chr = ($bp1->{tname} =~ m/chr/) ? $bp1->{tname} : "chr".$bp1->{tname};
			my $g2_chr = ($bp2->{tname} =~ m/chr/) ? $bp2->{tname} : "chr".$bp2->{tname};
			next if($g2_chr =~ m/M/ || $g1_chr =~ m/M/);
			my ($g1_start, $g1_end, $g1_strand) = ( $bp1->{tstart}, $bp1->{tend}, $bp1->{qstrand});
			my ($g2_start, $g2_end, $g2_strand) = ($bp2->{tstart}, $bp2->{tend}, $bp2->{qstrand});
			print STDERR join("\t", $g1_chr, $g2_chr, $g1_strand, $g2_strand, $sc_cover), "\n" if($debug);

			my (@genes, $gene1, @genes2, $gene2);
			my ($f_tree1, $r_tree1) = ($gm->sub_model($g1_chr, "+"), $gm->sub_model($g1_chr, "-"));
			my ($f_tree2, $r_tree2) = ($gm->sub_model($g2_chr, "+"), $gm->sub_model($g2_chr, "-"));
			push @genes, $f_tree1->intersect([$g1_start, $g1_end]) if(defined($f_tree1));
			push @genes, $r_tree1->intersect([$g1_start, $g1_end]) if(defined($r_tree1));
			push @genes2, $f_tree2->intersect([$g2_start, $g2_end]) if(defined($f_tree2)); 
			push @genes2, $r_tree2->intersect([$g2_start, $g2_end]) if(defined($r_tree2));

			my $same_gene = 0;
			$same_gene = 1 if((!@genes || !@genes2) && $g1_chr eq $g2_chr && 
					  (abs($pos1-$g2_start) < 1000 || abs($pos1-$g2_end) < 1000 ||
					   abs($pos2-$g1_start) < 1000 || abs($pos2-$g1_end) < 1000));
			$same_gene = 1 if((!@genes && !@genes2) && $g1_chr eq $g2_chr && 
					   (abs($pos1-$g2_start) < 5000 || abs($pos1-$g2_end) < 5000 ||
					   abs($pos2-$g1_start) < 5000 || abs($pos2-$g1_end) < 5000 ||
					   (abs($pos1-$pos2) < $min_fusion_distance/10 && $bp1->{qstrand} eq $bp2->{qstrand})));

				if(!@genes){ 
					push @genes, $f_tree1->intersect([$g1_start-5000, $g1_end+5000]) if(defined($f_tree1)); 
					push @genes, $r_tree1->intersect([$g1_start-5000, $g1_end+5000]) if(defined($r_tree1)); 
				}
				if(!@genes2 && !$same_gene){ 
					push @genes2, $f_tree2->intersect([$g2_start-5000, $g2_end+5000]) if(defined($f_tree2)); 
					push @genes2, $r_tree2->intersect([$g2_start-5000, $g2_end+5000]) if(defined($r_tree2)); 
				}
				if(!@genes){ 
					push @genes, $f_tree1->intersect([$g1_start-10000, $g1_end+10000]) if(defined($f_tree1)); 
					push @genes, $r_tree1->intersect([$g1_start-10000, $g1_end+10000]) if(defined($r_tree1)); 
				}
				if(!@genes2 && !$same_gene){ 
					push @genes2, $f_tree2->intersect([$g2_start-10000, $g2_end+10000]) if(defined($f_tree2)); 
					push @genes2, $r_tree2->intersect([$g2_start-10000, $g2_end+10000]) if(defined($r_tree2)); 
				}
			if(!@genes){ $gene1="NA"; }
			else{ $gene1=$genes[0]->val->name; }
			if(scalar @genes > 1){ for(my $i=1; $i<=$#genes; $i++) {$gene1=$gene1.",".$genes[$i]->val->name} }

			if(!$same_gene){ # to deal with overlapping genes
				foreach my $g1 (@genes){
					foreach my $g2 (@genes2){
				  		$same_gene = 1 if($g1 eq $g2 && $g1 ne 'NA');
				  		next if($g1_chr ne $g2_chr);
						$same_gene = 1 if($g1->val->start - 1000 < $g2->val->start && $g1->val->end + 1000 > $g2->val->end);
						$same_gene = 1 if($g2->val->start - 1000 < $g1->val->start && $g2->val->end + 1000 > $g1->val->end);
					}
				}
			}

			if(!@genes2){
				$gene2 = $same_gene ? $gene1 : "NA";
			}
			else{ $gene2=$genes2[0]->val->name; }
			if(scalar @genes2 > 1){ for(my $i=1; $i<=$#genes2; $i++) {$gene2=$gene2.",".$genes2[$i]->val->name} }

			my $out_string = join("\t", $bp1->{ort}, $g1_chr, $bp1->{tstart}, $bp1->{tend}, $bp1->{qstart}, $bp1->{qend}, 
				$bp1->{qstrand}, $bp1->{matches}, $bp1->{percent}, $bp1->{repeat}, $qseq, $bp2->{ort}, $g2_chr, 
				$bp2->{tstart}, $bp2->{tend}, $bp2->{qstart}, $bp2->{qend}, $bp2->{qstrand}, $bp2->{matches}, $bp2->{percent}, 
				$bp2->{repeat}, $sv->{gap}, $same_gene); 
			print STDERR join("\t", $sample, $gene1, $gene2, $sc_cover, $sc_cutoff, $pos1, $pos2, $out_string), "\n" if($debug); 
			if($same_gene) {print $UIF join("\t", $sample, $gene1, $gene2, $sc_cover, $sc_cutoff, $pos1, $pos2, $out_string), "\n";}
 			else { print $UFF join("\t", $sample, $gene1, $gene2, $sc_cover, $sc_cutoff, $pos1, $pos2, $out_string), "\n"; }
		}
		close($UFF);
		close($UIF);
}

sub second_sc_chk {
	my $sam = shift;
	my $bp = shift;
	my @reads;
	my $debug = 0;
	my ($chr, $pos) = ($bp->{tname},$bp->{tpos});
	my $clip = $bp->{qstrand}*$bp->{ort};
	print STDERR "bp_chk = get_sclip_reads(-SAM => $sam, -VALIDATOR => $validator, -CHR =>$chr, -POS => $pos, -CLIP => $clip, -MIN_SC_LEN => 3)\n" if($debug);
	my $bp_chk = get_sclip_reads(-SAM => $sam, -VALIDATOR => $validator, -CHR =>$chr, -POS => $pos, -CLIP => $clip, -MIN_SC_LEN => 3);
	push @reads, @{$bp_chk->{reads}} if(defined $bp_chk);
	return scalar @reads;
}

=head1 VERSION
This documentation refers to Cicero.pl version 0.1.7.

=head1 USAGE
	
	This program depends on several things that need to be installed and/or
	specified.  The program uses BioPerl and Bio::DB::Sam module to parse 
	the files and bam files.  Also it uses Blat software suites to do genome
	mapping and alignment.  To make the program efficient, it also requires
	a blat server setup.  And the program uses CAP3 assembler.

=head1 AUTHOR
Yongjin Li 
Liqing Tian (Liqing.Tian@STJUDE.ORG)

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

