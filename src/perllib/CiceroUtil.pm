package CiceroUtil;

use strict;
use Carp;
use English;
use Bio::SeqIO;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants;
use File::Spec;
use File::Path qw(remove_tree);
use File::Temp qw/ tempfile tempdir /;
use Cwd;
use CiceroSCValidator qw(LEFT_CLIP RIGHT_CLIP);
require Exporter;
our @ISA = ('Exporter');
our @EXPORT = qw(print_out is_new_pair prepare_reads_file parse_range
	is_PCR_dup read_fa_file read_config_file
	get_sclip_reads get_discordant_reads get_junction_reads rev_comp get_mate_reads);

sub rev_comp {
    my $str = shift;
    $str = reverse $str;
    $str =~ tr/ACTG/TGAC/;
    return $str;
}

sub is_garbage{

	my $seq = shift;
	my $max_Ns = 5;
	my $len = length($seq);
	my $n = @{[$seq =~ /(N)/g]};
	return 1 if($seq =~ /(N{$max_Ns,})/ || $n > 0.3*$len);
	return 0;

}

# something about binormial distribution
sub choose {
    my ($n,$k) = @_;
    my ($result,$j) = (1,1);
	    
    return 0 if $k>$n||$k<0;
    $k = ($n - $k) if ($n - $k) <$k;
	    
    while ($j <= $k ) {
        $result *= $n--;
        $result /= $j++;
    }
    return $result;
}

sub dbinorm { #desity funtion 
	my($n, $k, $p) = @_;
	return unless($n > 0 && $k >= 0 && $n >= $k && $p >= 0 && $p <=1);
	my $prob = log(choose($n, $k)) + ($n-$k)*log(1-$p) + $k*log($p);
	return exp($prob);
}
sub find_smallest_cover {
	my ($c, $p, $crit) = @_;
	my @s_cover;
	for( my $i = 1; $i < $c; $i++) {
		my $n = 1;
		while(1) {
			my $tmp = 0;
			for( my $j = 0; $j <= $i && $j < $n; $j++) {
				$tmp += dbinorm($n, $j, $p);
			}
			if( $tmp < $crit) {
				$s_cover[$i] = $n - 1;
				last;
			}
			$n++;
		}
	}
	return \@s_cover;
}

sub read_config_file {
	my $config_file = shift;
	my %config_values;
	open CFG, $config_file or die "Error: Unable to open $config_file\n";
	while (<CFG>)
	{
		chomp;
		my ($key, $value) = split(/\s+/, $_);
		next unless(defined $key && defined $value);
		$config_values{$key} = $value;
	}
	close CFG;
	return \%config_values;
}

sub read_fa_file {
	my $file = shift;
	my $in = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
	my %seqs;
	while( my $seq=$in->next_seq()) {
		$seqs{$seq->display_id} = $seq->seq;
	}
	return \%seqs;
}

# rmdup: remove PCR duplication related code
# you can replace it using your own criteria 
sub is_new_pair {
    my ($a, $sclip_len, $pairs) = @_;
    foreach my $p (@{$pairs}) {
        return undef if($a->start == $p->[0] && $a->end == $p->[1] &&
            $a->mate_start == $p->[2] && $a->mate_end == $p->[3] && $sclip_len == $p->[4]);
    }
    return 1;
}

# check PCR duplicate
sub is_PCR_dup {
	my ($a, $pairs, $sclip_len) = @_;
	my ($s, $e, $ms, $me ) = ($a->start, $a->end, $a->mate_start ? $a->mate_start : 0,
		$a->mate_end ? $a->mate_end : 0);
	return 1 if($a->flag & RFLAGS->{DUPLICATE});
    foreach my $p (@{$pairs}) {
		return 1 if($s == $p->[0] && $e == $p->[1] &&
			$ms == $p->[2] && $me == $p->[3] && $sclip_len == $p->[4]);
	}
	return;
}

# parse the range of input, format is: chr:start-end
# start and end is optional
sub parse_range {
	my $range = shift;
	my ($chr, $start, $end);
	my @field = split /:/, $range;
	$chr = $field[0];
#	$chr = substr($chr, 3) if($chr =~ /^chr/);
	if(@field > 1) {
		@field = split /-/, $field[1];
		$start = $field[0];
		$end = $field[1] if($field[1]);
		if($start !~ m/^\d+$/ or $end !~ m/^\d+$/) {
			croak "wrong range format, need to be: chr:start-end";
		}
	}
	return ($chr, $start, $end);
}

sub get_junction_reads {
	my %args = @_;
	my ($sam, $chr, $site, $extend_length, $clip, $gap_size) = 
		($args{-SAM}, $args{-CHR}, $args{-site}, $args{-extend_length}, $args{-CLIP}, $args{-GAP_SIZE});

	my $debug = 0; 
	print STDERR "=== get_junction_reads  === $chr:$site\t$clip\t$gap_size\n" if($debug);
	$extend_length = 15 unless(defined $extend_length);
	my ($start, $end) = ($site - $extend_length, $site + $extend_length);
	my $range = "$chr:$start-$end";

	print STDERR "range: $range\n" if($debug);
	my @reads;
	#my $strand_validator = $paired ? 1 : 0;
	$sam->fetch($range, sub { #to find discordant reads
		my $a = shift;
		my $is_paired = $a->paired;
		my $is_proper = $a->proper_pair; 
		my $is_unmapped = $a->unmapped;
		my $mate_is_unmapped = $a->munmapped;
		my $reversed = $a->reversed; 

		return unless($is_paired);
		return if($is_unmapped || $mate_is_unmapped);
		return if($a->flag & RFLAGS->{DUPLICATE}); 
		return unless($a->cigar_str =~ /N/);
		my @cigar_array = @{$a->cigar_array};

		#return if($debug && $a->qname ne 'HWI-ST988:130:D1TFEACXX:3:1111:6222:25769');
		print $a->qname, "\t", $a->start, "\t", $a->cigar_str, "\n"  if($debug);
		my $n = scalar @cigar_array;
		return if($n < 3);
		my ($gap_start_site, $gap_end_site) = ($a->start, $a->start);
		my $junc = 0;
 		for(my $i=0; $i<$n; $i++){
			print "--- ", join("\t", $cigar_array[$i]->[0], $i, $cigar_array[$i]->[1], $cigar_array[$i-1]->[1], $cigar_array[$i+1]->[1]), "\n"  if($debug); 
			if($cigar_array[$i]->[0] ne 'N' || ($i>0 && $cigar_array[$i-1]->[1] <= 3) || $cigar_array[$i+1]->[1] <= 3){
				$gap_start_site += $cigar_array[$i]->[1]; 
				print join("\t", $clip, $gap_start_site, $gap_end_site, $start, $end), "---\n" if($debug);
				next;
			}
			else{
				last unless($cigar_array[$i-1]->[0] eq 'M' && $cigar_array[$i-1]->[1] > 10 && 
					    $cigar_array[$i+1]->[0] eq 'M' && $cigar_array[$i+1]->[1] > 10 &&
					    abs($cigar_array[$i]->[1] - $gap_size) < 10);
				$gap_end_site = $gap_start_site + $cigar_array[$i]->[1];
				print join("\t", $gap_size, $gap_start_site, $gap_end_site, $start, $end), "\n" if($debug);
				if($clip > 0 && $gap_start_site > $start && $gap_start_site < $end){
					print "right-clip and junc = 1\n" if($debug);
					$junc = 1; last;
				}
				if($clip < 0 && $gap_end_site > $start && $gap_end_site < $end){
					print "left-clip and junc = 1\n" if($debug);
					$junc = 1; last;
				}
				$gap_start_site += $cigar_array[$i]->[1];
			}
		}

		return unless($junc);
		print join("\t", $a->qname, $a->start, $a->cigar_str), "\n" if($debug);

		my ($qname, @qscore) = ($a->qname, $a->qscore);
		my $strand = $reversed ? -1: 1; # mate reads strand
		push @reads, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			cigar => $a->cigar_str,
			qual => \@qscore,
			strand => $strand,
		};
	}); # @readsA
	print STDERR "number of junc reads: ", scalar @reads, "\n" if($debug);
	return @reads;
}

sub get_unmapped_reads {

	my %args = @_;
	my ($sam, $chr, $start, $end, $strand, $validator, $paired, $rmdup) = 
		($args{-SAM}, $args{-CHR}, $args{-START}, $args{-END}, $args{-STRAND}, 
		$args{-VALIDATOR}, $args{-PAIRED}, $args{-RMDUP} );
	my @rtn;
	my $debug = 0;
	print "\n=== get_unmapped_reads ===\n" if($debug);
	my $strand_validator = $paired ? 1 : 0;
	my $range = $chr;
	my $range2 = '';
	my ($start2, $end2);
	print STDERR "start: $start\tend: $end\n" if($debug);
	if($end - $start <= 100){
		$start -= 100;
		$end += 100;
	}
	$range = $chr . ":" . $start . "-" . $end if($start && $end);

	if($end - $start > 800){
		$end2 = $start + 400;
		$range = $chr . ":" . $start . "-" . $end2;
		$start2 = $end - 400;
		$range2 = $chr . ":" . $start2 . "-" . $end;
	} 

	print STDERR "range1: $range\n" if($debug);
	$sam->fetch($range, sub {
		my $a = shift;
		return if($a->flag & RFLAGS->{DUPLICATE}); 
		my $is_paired = $a->paired;
		my $is_proper = $a->proper_pair; 
		my $is_unmapped = $a->unmapped;
		my $mate_is_unmapped = $a->munmapped;
		my $reversed = $a->reversed; 
		my $mate_reversed = $a->mreversed;

		return unless ($is_paired && $is_unmapped && !$mate_is_unmapped);
		my $this_strand = $mate_reversed ? -1: 1; # mate reads strand
		return unless($this_strand == $strand);
		return unless($a->mate_start >= $start && $a->mate_start <= $end);
		return unless(is_high_quality($a));
		return if(is_garbage($a->query->dna));

		my $qscore = $a->qscore;
		print join("\t", $a->qname, $a->flag, $a->query->dna, $is_unmapped, $a->unmapped, $a->paired, $a->mate_start, $this_strand), "\n" if($debug);
		push @rtn, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			qual => $qscore,
		};
	}
	);
	print STDERR "number of reads from range1 is ", scalar @rtn, "\n" if($debug);

	print STDERR "range2: $range2\n" if($debug);
	$sam->fetch($range2, sub {
		my $a = shift;
		return if($a->flag & RFLAGS->{DUPLICATE}); 
		my $is_paired = $a->paired;
		my $is_unmapped = $a->unmapped;
		my $mate_is_unmapped = $a->munmapped;
		my $mate_reversed = $a->mreversed;

		return unless ($is_paired && $is_unmapped && !$mate_is_unmapped);
		my $this_strand = $mate_reversed ? -1: 1; # mate reads strand
		return if($this_strand != $strand);

		return unless($a->mate_start >= $start && $a->mate_start <= $end);
		return unless(is_high_quality($a));
		return if(is_garbage($a->query->dna));

		my @qscore = $a->qscore;
		print join("\t", $a->qname, $a->flag, $a->query->dna, $is_unmapped, $a->unmapped, $a->paired, $a->mate_start, $this_strand), "\n" if($debug);
		push @rtn, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			qual => \@qscore,
		};
	}
	) if($range2 ne '');
	print STDERR "number of unmapped reads is ", scalar @rtn, "\n" if($debug);
	return @rtn;
}

sub is_high_quality {
#	print STDERR "quality validator\n";
	my $a = shift;
	my ($lowqual_cutoff, $min_percent_hq) = (20, 80);
#	my $cigar_str = $a->cigar_str;
	my @cigar_array = @{$a->cigar_array};
	my @qual = $a->qscore;
	my $n_hq = 0;
	for( my $i = 0; $i <= $#qual; $i++) {
		$n_hq++ if($qual[$i] >= $lowqual_cutoff);
	}
	return ($n_hq * 100 > $#qual * $min_percent_hq) ? 1 : undef;
}

sub get_mate_reads {

	my %args = @_;
	my ($sam, $chr, $start, $end, $array_ref) = 
		($args{-SAM}, $args{-CHR}, $args{-START}, $args{-END}, $args{-READS});

	my $debug = 0;
	print "=== get_mate_reads ===\n" if($debug);
	my @rtn;
	my $range = $chr;
	my $range2 = '';
	my ($start2, $end2);
	print STDERR "start: $start\tend: $end\n" if($debug);
	if($end - $start <= 100){
		$start -= 100;
		$end += 100;
	}
	$range = $chr . ":" . $start . "-" . $end if($start && $end);
=pod
	my $range_cutoff = 800;
	if($end - $start > $range_cutoff){
		$end2 = $start + $range_cutoff/2;
		$range = $chr . ":" . $start . "-" . $end2;
		$start2 = $end - $range_cutoff/2;
		$range2 = $chr . ":" . $start2 . "-" . $end;
	} 
	print STDERR "range1: $range\trange2: $range2\n" if($debug);
=cut
	my @reads = @{$array_ref};
	my %reads_hash = ();
	foreach my $r (@reads){
		my $qname = $r->{name};
		my $end = ($r->{sam_flag} & 0x0040)? 1 : 2;
		$reads_hash{$qname} = $end;
	}

	my %mate = ();

	$sam->fetch($range, sub {
		my $a = shift;

		my $is_paired = $a->paired;
		my $is_unmapped = $a->unmapped;
		return unless($is_paired);
		return if($is_unmapped);

		my $qname = $a->qname;
		return unless (exists($reads_hash{$qname}));

		my $end = ($a->flag & 0x0040)? 1 : 2;
		return if($end == $reads_hash{$qname});
		return if(is_garbage($a->query->dna));
	 
		#	print STDERR "found mate pair for ", $r->{name}, ", it is", $r_name, "\n";
		my @qscore = $a->qscore;
		print STDERR $a->qname, "\n" if($debug);
		push @rtn, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			qual => \@qscore,
		};
		$mate{$a->qname} = 1;
	});
	print STDERR "finished range1\n" if($debug);

	$sam->fetch($range2, sub {
		my $a = shift;

		my $is_paired = $a->paired;
		my $is_unmapped = $a->unmapped;
		return unless($is_paired);
		return if($is_unmapped);

		my $qname = $a->qname;
		return unless (exists($reads_hash{$qname}));

		my $end = ($a->flag & 0x0040)? 1 : 2;
		return if($end == $reads_hash{$qname});
		return if(is_garbage($a->query->dna));
	 
		#	print STDERR "found mate pair for ", $r->{name}, ", it is", $r_name, "\n";
		my @qscore = $a->qscore;
		print STDERR $a->qname, "\n" if($debug);
		push @rtn, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			qual => \@qscore,
		} unless(exists($mate{$a->qname}));
	}) if($range2 ne '');
	return @rtn;
}

sub get_discordant_reads {
	my %args = @_;
	my ($sam, $rangeA, $rangeB) = 
		($args{-SAM}, $args{-rangeA}, $args{-rangeB});

	my ($chrA, $leftA, $rightA) = split /[:|-]/, $rangeA;
	my ($chrB, $leftB, $rightB) = split /[:|-]/, $rangeB;

	my %readsA = ();
	my %readsB = ();
	my (@readsA, @readsB);
	$sam->fetch($rangeA, sub { #to find discordant reads
		my $a = shift;
		my $is_paired = $a->paired;
		my $is_proper = $a->proper_pair; 
		my $is_unmapped = $a->unmapped;
		my $mate_is_unmapped = $a->munmapped;
		my $reversed = $a->reversed; 

		return unless($is_paired);
		return if($is_unmapped || $mate_is_unmapped);
		return unless($a->start >= $leftA && $a->end <= $rightA);	
		return if($a->flag & RFLAGS->{DUPLICATE}); 

		my ($qname, $qscore) = ($a->qname, $a->qscore);
		my $strand = $reversed ? -1: 1; # mate reads strand
		$readsA{$qname} = $strand;
		push @readsA, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			qual => $qscore,
			strand => $strand,
		};
	}); # @readsA

	$sam->fetch($rangeB, sub { #to find discordant reads
		my $a = shift;
		my $is_paired = $a->paired;
		my $is_proper = $a->proper_pair; 
		my $is_unmapped = $a->unmapped;
		my $mate_is_unmapped = $a->munmapped;
		my $reversed = $a->reversed; 

		return unless($is_paired);
		return if($is_unmapped || $mate_is_unmapped);
		return unless($a->start >= $leftB && $a->end <= $rightB);	
		return if($a->flag & RFLAGS->{DUPLICATE}); 

		my ($qname, $qscore) = ($a->qname, $a->qscore);
		my $strand = ($a->flag & RFLAGS->{REVERSED}) ? -1: 1; # mate reads strand
#		return unless($a->start == 54247864);
		return unless(exists($readsA{$qname}));
		$readsB{$qname} = $strand;
		push @readsB, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			qual => $qscore,
			strand => $strand,
		};
	}); # @readsB

	foreach my $r (@readsA){
		push @readsB, $r if(exists($readsB{$r->{name}}));
	}
	return @readsB;
}

sub get_del_reads {
	my %args = @_;
	my ($sam, $chr, $pos, $validator, $del_len) = 
	   ($args{-SAM}, $args{-CHR}, $args{-POS}, $args{-VALIDATOR}, $args{-DEL_LEN});

	my $debug = 0;
	print "== get_del_reads  ==\n" if($debug);
	my @rtn;
	my $range = $chr;
	$range = $chr . ":" . $pos . "-" . $pos;
	print STDERR "range: $range\n" if($debug);

	my $del_tag = $del_len."N" if($del_len);
	$sam->fetch($range, sub { #to find reads with one end mapped, the other end unmapped
		my $a = shift;
		return if($a->flag & RFLAGS->{DUPLICATE}); 
		return unless($a->start && $a->end);	
#		return unless($a->start == 54247864);
		my $cigar_str = $a->cigar_str;

		return if($cigar_str !~ m/$del_tag/);
		return if ($a->start > $pos || $a->end < $pos);

	});
}

sub get_sclip_reads {
	my %args = @_;
	my ($sam, $chr, $sc_site, $clip, $min_len_cutoff, $sc_shift, $validator, $paired, $rmdup, $fixSC) = 
		($args{-SAM}, $args{-CHR}, $args{-POS}, $args{-CLIP}, $args{-MIN_SC_LEN}, $args{-SC_SHIFT},
		 $args{-VALIDATOR}, $args{-PAIRED}, $args{-RMDUP}, $args{-FIXSC});

	my $debug = 0;
	print STDERR "== get_sclip_reads  ==\n" if($debug);
	my @rtn;
	my $mt_start0=1e20;
	my $mt_end0=0;
	my $max_sclip_len = 0;
	$min_len_cutoff = 1 unless $min_len_cutoff;
	$sc_shift = 3 unless $sc_shift;
	my $range = $chr;
	my ($start, $end) = ($sc_site - $sc_shift, $sc_site + $sc_shift);
	$range = $chr . ":" . $start . "-" . $end if($start && $end);
	print STDERR "range: $range\tfixSC: $fixSC\n" if($debug);

	my @pairs;
	my $strand_validator = $paired ? 1 : 0;
	my ($plus, $minus) = (0, 0);

	my %sc_reads = ();
	my %sc_pos = ();
	$sam->fetch($range, sub { #to find reads with one end mapped, the other end unmapped
		my $a = shift;
		return if($a->flag & RFLAGS->{DUPLICATE}); 
		return unless($a->start && $a->end);	
#		return unless($a->start == 54247864);
		my $cigar_str = $a->cigar_str;

		return if($cigar_str !~ m/S/);
		#return unless ($a->qname =~ /HWI-ST671:118:D05A3ACXX:2:2306:7105:92919/);
		print STDERR "\n getSC\t", join("\t", $a->qname, $a->cigar_str, $a->start, $a->end, $a->proper_pair, "clip:".$clip), "\n" if($debug);
		my @cigar_array = @{$a->cigar_array};

        	my ($seq, $qual);
		my ($left_validated, $right_validated, $boundary_softclip, $sclip_len, $pos) = (0, 0, 0, 0, 0);

		$pos = $a->start if($clip == LEFT_CLIP && $cigar_array[0]->[0] eq 'S');
		$pos = $a->end if($clip == RIGHT_CLIP && $cigar_array[-1]->[0] eq 'S');
		$boundary_softclip = 1 if($pos >= $start && $pos <= $end); # the softclipping position is in range
		print STDERR $a->qname, "\tafter boundary check and boundary_softclip is $boundary_softclip\npos:\t$pos\nstart:\t$start\nend:\t$end\n" if($debug);

		if($fixSC){
			return if($a->end < $start or $a->start > $end); # the softclipping position is not in range
			return if($boundary_softclip == 0 && scalar @cigar_array < 4);
		}
		else{
			return if($pos > $a->end or $pos < $a->start); 
			return if($boundary_softclip == 0);
		}

		if($clip == LEFT_CLIP && $cigar_array[0]->[0] eq 'S'){
            		$sclip_len = $cigar_array[0]->[1];
			return if($sclip_len < 3);
		#	print "left pos = ", $a->start, "\n" if($debug);
			$left_validated = 1 if($validator && $validator->validate($a, $clip));
		}

		if($clip == RIGHT_CLIP && $cigar_array[-1]->[0] eq 'S'){
         		$sclip_len = $cigar_array[-1]->[1]; 
			return if($sclip_len < 3);
			#print "right pos = ", $a->end, "\n" if($debug);
			$right_validated = 1 if($validator && $validator->validate($a, $clip));
		}
		if($left_validated + $right_validated == 0){
			print STDERR "not valid SC\n" if($debug); 
			return;
		} # wrong clip signal

		if(is_garbage($a->query->dna)) {print STDERR $a->qname," is garbage\n" if($debug); return;}
		print STDERR "$left_validated + $right_validated == 0\npos: $pos\n", "a->start: ", $a->start, "\ta->end: ", $a->end, "\t", $a->query->dna,"\n" if($debug);
		 
		if(!$boundary_softclip && scalar @cigar_array >= 4){ # to check unmapped soft clip sites

			if($left_validated && $cigar_array[1]->[1] <= 6 &&
			   $cigar_array[1]->[0] eq 'M' && $cigar_array[2]->[0] eq 'N' && $cigar_array[2]->[1] > 100){
				$pos = $a->start + $cigar_array[1]->[1] - $cigar_array[2]->[1];
				$sclip_len = $cigar_array[0]->[1] + $cigar_array[1]->[1];
			}

			if($right_validated && $cigar_array[-2]->[1] <= 6 && 
			   $cigar_array[-2]->[0] eq 'M' && $cigar_array[-3]->[0] eq 'N' && $cigar_array[-3]->[1] > 100){
				$pos = $a->end - $cigar_array[-2]->[1] - $cigar_array[-3]->[1];
            			$sclip_len = $cigar_array[-1]->[1] + $cigar_array[-2]->[1]; 
			}
			return unless($pos >= $start && $pos <= $end); # adjusted softclipping position is in range
		}

		#return if($sclip_len < $min_len_cutoff); 
		#print STDERR join("\t", $sclip_len, $max_sclip_len, $min_len_cutoff), "\n";
		$max_sclip_len = $sclip_len if($max_sclip_len < $sclip_len);
		return if(@pairs > 0 && $rmdup && is_PCR_dup($a, \@pairs, $sclip_len)); #it's a PCR duplicate
		push @pairs, [$a->start, $a->end, $a->mate_start ? $a->mate_start : 0, 
				$a->mate_end ? $a->mate_end : 0, $sclip_len] if($rmdup);

		my $mt_pos = $a->mate_start;
		($a->flag & RFLAGS->{M_REVERSED}) ? $minus++: $plus++; # mate reads strand
		my $mate_start = $a->mate_start;
		$mt_start0 = $mate_start if($mt_start0 >= $mate_start);
		$mt_end0 = $mate_start if($mt_end0 <= $mate_start);
		my @qscore = $a->qscore;
		#return if (exists($sc_reads{$a->qname})); # both ends of the read is soft-clipped
		print STDERR $a->qname, " is a valid sc read\n" if($debug);
		$sc_reads{$a->qname} = 1;
		push @rtn, {
			name => $a->qname,
			sam_flag => $a->flag,
			seq => $a->query->dna,
			qual => \@qscore,
			sc_pos => $pos
		};

		# to find the number of sc reads at each site
		if(exists($sc_pos{$pos})) {
			$sc_pos{$pos}++;
		}
		else{
			$sc_pos{$pos} = 1;
		}
		}
	);
	print STDERR "-=finished sam=-\nnumber of reads: ", scalar @rtn, "\n" if($debug);

	print STDERR "max_sclip_len: $max_sclip_len\n" if($debug);
	return undef if($max_sclip_len < $min_len_cutoff);
	
	# to find the site with the highest number of soft-clip reads
	my ($pos, $max) =(0, 0);
	foreach my $m (keys %sc_pos) {
		if($sc_pos{$m}>$max){
			$max=$sc_pos{$m};
			$pos=$m;
		}
	}

	my $strand = $plus > $minus ? 1 : -1;
	#return @rtn;
	my %SC = (
		reads		=> \@rtn,
		pos		=> $pos,
		mt_start	=> $mt_start0,
		mt_end		=> $mt_end0,
		mt_strand	=> $strand # strand of mate of SC reads
	);
	return \%SC;
}

sub prepare_reads_file{

	my %args = @_;
	my ($out_file, $sam, $chr, $sc_site, $clip, $validator, $paired, $rmdup, $min_SC, $min_sclip_len, $unmapped_cutoff, $gap_size, $fixSC, $sc_shift, $output_mate) = 
		($args{-OUT}, $args{-SAM}, $args{-CHR}, $args{-POS}, $args{-CLIP}, $args{-VALIDATOR}, $args{-PAIRED}, $args{-RMDUP}, 
		 $args{-MIN_SC}, $args{-MIN_SC_LEN}, $args{-UNMAPPED_CUTOFF}, $args{-GAP_SIZE}, $args{-FIXSC}, $args{-SC_SHIFT}, $args{-MATE});

		my $debug = 0;
		print STDERR "\n=== prepare_reads_file ===\n$out_file\n" if($debug);
		$min_sclip_len = 20 unless(defined $min_sclip_len);
		$output_mate = 1 unless(defined $output_mate);
		$sc_shift = 3 unless $sc_shift;
		print STDERR "to get soft-clip reads $chr:$sc_site, $clip, $min_sclip_len ...\n" if($debug);
		my $sc = get_sclip_reads(-SAM => $sam,
					   -CHR =>$chr, 
					   -POS => $sc_site, 
				   	-CLIP => $clip, 
				   	-SC_SHIFT=> $sc_shift,
				   	-MIN_SC_LEN=> $min_sclip_len,
				   	-VALIDATOR => $validator,
				   	-PAIRED => $paired,
				   	-RMDUP => $rmdup,
				   	-FIXSC => $fixSC,
			        	);
		print STDERR "get_sclip_reads done!\n" if($debug);
		my (@reads, $mt_strand, $mt_start, $mt_end);
		if(defined $sc) {
			push @reads, @{$sc->{reads}};
			$mt_strand = $sc->{mt_strand};
			$mt_start = $sc->{mt_start};
			$mt_end = $sc->{mt_end};
			print STDERR "info for mate reads: ", join("\t", $mt_strand, $mt_start, $mt_end, $sc_site), "\n" if($debug);
			$mt_start = $sc_site + 10 if($mt_start > $sc_site + 10);
			$mt_end = $sc_site - 10 if($mt_end < $sc_site - 10);
			print STDERR "got ", scalar @{$sc->{reads}}, " SC reads\n" if($debug);
			print STDERR "info for unmapped reads: ", join("\t", $mt_strand, $mt_start, $mt_end), "\n" if($debug);
		}

		print STDERR "to get junction reads $chr:$sc_site, $clip, *$gap_size* ...\n" if($debug && $gap_size);
		push @reads, get_junction_reads(-SAM => $sam, -CHR => $chr, -site => $sc_site, -CLIP=>$clip, -GAP_SIZE=>$gap_size)  if($gap_size && $gap_size > 10000);
		my $n_sc = scalar @reads;
		print STDERR "got ", $n_sc, " SC or junction reads\n" if($debug);
		return if($n_sc < $min_SC);

		if ($n_sc >= 2000 || !defined($sc)){
			print STDERR "print_out(-READS=>\@reads, -OUT_FILE=>$out_file, -FORMAT=>'fa')\n" if($debug);
			print_out(-READS=>\@reads, -OUT_FILE=>$out_file, -FORMAT=>'fa');
			return 1;
		}
		my @unmapped_reads;
		@unmapped_reads = get_unmapped_reads(-SAM => $sam,
					   -CHR =>$chr, 
					   -START => $mt_start, 
					   -END => $mt_end, 
				   	-STRAND => $mt_strand, 
				   	-VALIDATOR => $validator,
				   	-PAIRED => $paired,
				   	-RMDUP => $rmdup,
			) if($n_sc <= $unmapped_cutoff);
		print STDERR "number of unmapped reads: ",  scalar @unmapped_reads, "\n" if($debug);
		my $n_um = scalar @unmapped_reads;
		if($n_um < $n_sc*50 && $n_um < 500){
			push @reads, @unmapped_reads;
		}
		else{
			print STDERR "too many unmapped at ", join("\t", $chr, $sc_site, $mt_start, $mt_end, $mt_strand, $paired, $rmdup, $n_sc, $n_um), "\n";
		}

		if($output_mate){
		print STDERR "to get mate reads $chr:$mt_start-$mt_end ...\n" if($debug);
		my @mate_reads;
		@mate_reads = get_mate_reads(-SAM => $sam,
					   -CHR =>$chr, 
					   -START => $mt_start, 
					   -END => $mt_end, 
 					   -READS => \@reads,
				) if($paired && scalar @reads < 2000);
		
		print STDERR "number of mate reads: ",  scalar @mate_reads, "\n" if($debug);
		push @reads, @mate_reads;
		}
		print STDERR "print_out(-READS=>\@reads, -OUT_FILE=>$out_file, -FORMAT=>'fa')\n" if($debug);
		print_out(-READS=>\@reads, -OUT_FILE=>$out_file, -FORMAT=>'fa', -RMDUP => $rmdup);
		return 1;
}

sub print_out{

	my %args = @_;
        my ($array_ref, $out_file, $format, $rmdup) = ($args{-READS}, $args{-OUT_FILE}, $args{-FORMAT}, $args{-RMDUP});
        my @reads = @{$array_ref};

	my $nr = scalar @reads;
	my $rand_cutoff = 1;
	$rand_cutoff = 2000/$nr if($nr > 2000);

	open(my $fa_fh, ">$out_file");
	#my $out_file_Q = $out_file.".qual";
	#open(my $qual_fh, ">$out_file_Q");
	my (%end1_reads, %end2_reads);
	foreach my $r (@reads){
		#next unless(exists($end1_reads{$r->{name}}) && exists($end2_reads{$r->{name}}));
		#print STDERR "testing ... \n";

		#next if($rmdup && $r->{sam_flag} & 0x0400); #PCR duplicates
		next if($r->{sam_flag} & 0x0200); #QC failure

		if($rand_cutoff < 1){
			my $rrr = rand();
			next if(rand() > $rand_cutoff);
		}

		my $r_ort = ($r->{sam_flag} & 0x0040)? 'f': 'r'; # read orientation
		if($r_ort eq 'f'){
			next if(exists($end1_reads{$r->{name}}));
			$end1_reads{$r->{name}} = 1;
		}
		else{
			next if(exists($end2_reads{$r->{name}}));
			$end2_reads{$r->{name}} = 1;
		}
		my $seq = $r->{seq};
		my @qual = @{$r->{qual}};
		if($r->{sam_flag} & 0x0010){
			$seq = rev_comp($r->{seq});
			@qual = reverse(@{$r->{qual}});
			#$score_str = $align->qstring;
		}
		if($format eq 'fq'){
			$r_ort = ($r->{sam_flag} & 0x0040)? '1': '2'; # read orientation
			#print $fa_fh "@", $r->{name}, "/", $r_ort,"\n+\n", $seq, "\n", $align->qstring,
			
		}
		$r_ort = ($r->{sam_flag} & 0x0040)? 'f': 'r'; # read orientation
		print $fa_fh ">", $r->{name}, ".", $r_ort,"\n", $seq, "\n" if($format eq 'fa');
		#print $qual_fh ">", $r->{name}, ".", $r_ort,"\n", join(" ", @qual), "\n";
	}
	close($fa_fh); 
	#close($qual_fh);
}

my $start_dir;
BEGIN {
	$start_dir = getcwd;
}
END {
	chdir($start_dir);
#    remove_tree($work_dir) ;
}

1;
