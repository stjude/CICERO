#!/usr/bin/env perl 

use strict;
use warnings; 

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
my $min_sclip_reads = 2;#Liqing
my $optionOK = GetOptions(
	'o|out_dir=s'	=> \$out_dir,
	'genome=s'	=> \$genome,
	'p=s'			=> \$out_prefix,
	'l=i'			=> \$read_len,
	's|split_every=i'	=> \$split_every_n,
	'ratio=f'		=> \$expression_ratio_cutoff,
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

my $black_list_file = $conf->{'BLACKLIST_GENES'};
open(my $BLK, "<", "$black_list_file");
my %black_list = ();
print STDERR "Parsing black list of genes $black_list_file\n"; 
while(<$BLK>){
	my $line = $_;
	chomp($line);
	$black_list{$line} = 1;
}
close($BLK);

my %gene_info;
#$gene_info_file = "$out_dir/gInfo.txt";
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
`cat $out_dir/$out_prefix*.cover > $SC_file`;
print STDERR $SC_file, "\n";
open my $SCI, "$SC_file";

my @intra_SCs = ();
my @intergenic_SCs = ();
print STDERR "Parsing Soft Clips file: $SC_file\n"; 
while(<$SCI>){
	chomp;
	my ($chr, $site, $strand, $sc_cnt, $cover, $psc, $nsc, $pn, $nn) = split(/\t/);
	my $SC_MAF = $psc/($pn+0.01) > $nsc/($nn+0.01) ? $psc/($pn+0.01) : $nsc/($nn+0.01);
	next if($sc_cnt < $min_sclip_reads);#Liqing
	next unless (exists($gene_info{$chr}));#Liqing
	my @gInfo = @{$gene_info{$chr}};
	my $intra = 0;
	foreach my $g (@gInfo){
		my ($gene_name, $gStart, $gEnd, $mRNA_length, $r_cnt, $sc_cutoff) =  @{$g};
		next if($r_cnt == 0);
		next if ($gStart > $site || $gEnd < $site);
		$intra = 1;
		last if ($sc_cutoff < 0.01);

		# to remove black list genes
		my $blk = 0;
		my @g_names = split(/,|\|/,$gene_name);
		foreach my $gname (@g_names){
			$blk = 1 if(exists($black_list{$gname}));
		}
		last if($blk==1);

		my $expression_ratio = ($cover/($read_len-20))/($r_cnt*$read_len/(2*$mRNA_length));
		if($expression_ratio > $expression_ratio_cutoff || $SC_MAF > 0.05){
		#if($cnt > $sc_cutoff){
			push @intra_SCs, [$chr, $site, $strand, $sc_cnt, $sc_cutoff, $expression_ratio, $cover] if($chr);
			last;
		}
	}
	next if($sc_cnt < 5 || $intra);
	#push @intergenic_SCs, [$chr, $site, $strand, $sc_cnt, -1, 0, $cover] if($intra==0 && $sc_cnt > 5);#Liqing
	push @intra_SCs, [$chr, $site, $strand, $sc_cnt, -1, 0, $cover];#Liqing
}
close $SCI;
#print "number of SCs: ", scalar @intra_SCs, "\t", scalar @intergenic_SCs, "\n";#Liqing

#push @intra_SCs, @intergenic_SCs;#Liqing

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

