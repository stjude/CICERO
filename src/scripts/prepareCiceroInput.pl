#!/usr/bin/env perl 

use strict;
use warnings; 

use locale; 

use Carp;
use Getopt::Long;
use English;
use Pod::Usage;

use TdtConfig; 

my ($genome, $out_prefix, $gene_info_file, $out_dir, $read_len); 
my $expression_ratio_cutoff = 0.01;
my $split_every_n = 5000;
my $out_suffix = "sclip.txt";
my ( $help, $man, $version, $usage );
my $min_sclip_reads = 2;
my $optionOK = GetOptions(
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'	=> \$genome,
	'p=s'			=> \$out_prefix,
	'l=i'			=> \$read_len,
	's|split_every=i'	=> \$split_every_n,
	'ratio=f'		=> \$expression_ratio_cutoff,
	'm|min_sclip_reads=i'		=> \$min_sclip_reads,
	'f|gene_info_file=s' => \$gene_info_file,
	'h|help|?'		=> \$help,
	'man'			=> \$man,
	'usage'			=> \$usage,
	'v|version'		=> \$version,
);

pod2usage(-verbose=>2) if($man or $usage);
pod2usage(1) if($help or $version );

mkdir $out_dir if(!-d $out_dir);

my $conf; 
if (&TdtConfig::findConfig("genome", $genome)){
	$conf = &TdtConfig::readConfig("genome", $genome); 
}
else{
	croak "no config"; 
}

my $excluded_list_file = $conf->{'EXCLUDED_GENES'};
open(my $EXC, "<", "$excluded_list_file");
my %excluded_list = ();
print STDERR "Parsing excluded list of genes $excluded_list_file\n"; 
while(<$EXC>){
	my $line = $_;
	chomp($line);
	$excluded_list{$line} = 1;
}
close($EXC);

my $breakpoints_file = $conf->{'KNOWN_BREAKPOINTS'};
open(my $BRK, "<", "$breakpoints_file");
my %breakpoints = (); 
print STDERR "Parsing known break points $breakpoints_file\n"; 
while(<$BRK>){
	my $line = $_;
	chomp($line);
	my ($chr, $pos) = split(/\t/, $line); 
	$breakpoints{$chr}{$pos} = 1;
}
close($BRK);

my %gene_info;
print STDERR "Parsing gene info file\n"; 
print "\ngene_info_file: $gene_info_file\n";
open my $GI, "$gene_info_file";
while(<$GI>){
	chomp;
	my ($gene_name, $gRange, $strand, $mRNA_length, $cnt, $sc_cutoff) = split(/\t/);
	my ($chr, $gStart, $gEnd) = split /[:|-]/, $gRange;
	push @{$gene_info{$chr}}, [$gene_name, $gStart, $gEnd, $mRNA_length, $cnt, $sc_cutoff];
}
close $GI;

my $SC_file = "$out_dir/$out_prefix.SC.txt";

`ls $out_dir/$out_prefix*.cover | xargs cat > $SC_file`;
if ($?){
	my $err = $!;
	print STDERR "Error combining cover files: $err\n"; 
	exit $err;
}

print STDERR $SC_file, "\n";
open my $SCI, "$SC_file";

my @intra_SCs = ();
my @intergenic_SCs = ();
print STDERR "Parsing Soft Clips file: $SC_file\n"; 
while(<$SCI>){
	chomp;
	my ($chr, $site, $strand, $sc_cnt, $cover, $psc, $nsc, $pn, $nn) = split(/\t/);
	my $SC_MAF = $psc/($pn+0.01) > $nsc/($nn+0.01) ? $psc/($pn+0.01) : $nsc/($nn+0.01);
	next if($sc_cnt < $min_sclip_reads && (! exists $breakpoints{$chr}{$site} && ! exists $breakpoints{'chr'.$chr}{$site}));
	next unless (exists($gene_info{$chr}));
	my @gInfo = @{$gene_info{$chr}};
	my $intra = 0;
	foreach my $g (@gInfo){
		my ($gene_name, $gStart, $gEnd, $mRNA_length, $r_cnt, $sc_cutoff) =  @{$g};
		next if($r_cnt == 0);
		next if ($gStart > $site || $gEnd < $site);
		$intra = 1;
		last if ($sc_cutoff < 0.01);

		# to remove excluded list genes
		my $exclude = 0;
		my @g_names = split(/,|\|/,$gene_name);
		foreach my $gname (@g_names){
			$exclude = 1 if(exists($excluded_list{$gname}));
		}
		last if($exclude==1);

		# Compute an adjusted gene expression ratio
		my $expression_ratio = ($cover/($read_len-20))/($r_cnt*$read_len/(2*$mRNA_length));
		if($expression_ratio > $expression_ratio_cutoff || $SC_MAF > 0.05){
			push @intra_SCs, [$chr, $site, $strand, $sc_cnt, $sc_cutoff, $expression_ratio, $cover] if($chr);
			last;
		}
	}
	next if($sc_cnt < 5 || $intra);
	push @intra_SCs, [$chr, $site, $strand, $sc_cnt, -1, 0, $cover];
}
close $SCI;

print "start to write out seperate files...\n";
my $idx_bin=0;
my $n_SC_reads = 0; 
my $IN;
my $SC_info_file = "$out_dir/$out_prefix.$idx_bin.SC";
print "SC_info_file: $SC_info_file\n";
foreach my $sc (@intra_SCs){
		my ($chr, $site, $strand, $cnt, $sc_cutoff, $expression_ratio, $cover) = @{$sc};
		if($n_SC_reads > $split_every_n || $idx_bin == 0){

			$idx_bin++;
			$SC_info_file = "$out_dir/$out_prefix.$idx_bin.SC";
			print "SC_info_file: $SC_info_file\n";
			open $IN, ">$SC_info_file";
			print $IN join("\t", $chr, $site, $strand, $cnt, $sc_cutoff, $cover, $expression_ratio), "\n" if($idx_bin>0);
			$n_SC_reads = 0;
		}
		else{
			print $IN join("\t", $chr, $site, $strand, $cnt, $sc_cutoff, $cover, $expression_ratio), "\n";
			$n_SC_reads += $cnt; 
		}
} #end for chr
close $IN;

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
