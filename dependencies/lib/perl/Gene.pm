package Gene;
use strict;
use Transcript;
use Carp;
use Data::Dumper;

# a light weight gene structure is used here

my @Gene_slots;
BEGIN {
	@Gene_slots = qw(NAME CHR START END STRAND EXONS TRANSCRIPTS TYPE);
}
use enum @Gene_slots;

my %attribute = (
    	name         => NAME,
	chr	 	=> CHR,
    	start        => START,
    	end          => END,
	strand	     => STRAND,
	transcripts  => TRANSCRIPTS,
	type	     => TYPE,
);

#using an array instead of a hash for the node
sub _accessor {
    my $index = shift;
    return sub {
        my $self = shift;
        return undef unless $self;
        if (@_) {
          $self->[$index] = shift;
        }
        return $self->[$index];
    };
}

while(my($at, $idx) = each %attribute) {
    no strict 'refs';
    *$at = _accessor($idx);
}

sub new {
    my $class = shift;
    my $obj = [];
	$obj->[TRANSCRIPTS] = [];
    if (@_) {
		my %arg = @_;
        	$obj->[NAME]        = $arg{-NAME} if($arg{-NAME});
        	$obj->[CHR]         = $arg{-CHR} if($arg{-CHR});
		$obj->[START]       = $arg{-START} if($arg{-START});
		$obj->[END]         = $arg{-END} if($arg{-END});
		$obj->[STRAND]      = $arg{-STRAND} if($arg{-STRAND});
		$obj->[TRANSCRIPTS] = $arg{-TRANSCRIPTS} if($arg{-TRANSCRIPTS});
    }
    return bless $obj, $class;
}

sub add_transcript {
	my ($self, $fea) = @_;
	croak "You must add a Transcript type into a gene" 
		 unless ($fea->isa('Transcript'));
	if($self->[STRAND] && $self->[STRAND] ne $fea->strand) {
		croak "The transcript has different orientation with the gene";
	}
	if($self->[CHR] && $self->[CHR] ne $fea->chr) {
		croak "The transcript is on different chr with the gene";
	}
#	if($self->[TYPE] && $fea->type ne $fea->type) {
#		croak "The type of the transcript are different from the gene";
#	}
	$self->[STRAND] = $fea->strand;
	$self->[CHR] = $fea->chr;
#	$self->[TYPE] = $fea->type;

	push @{$self->[TRANSCRIPTS]}, $fea;

	if(defined $self->[NAME]) { 
		my @names = split(/,/,$self->[NAME]);
		my $exist = 0;
		foreach my $name (@names){
			if($name eq $fea->name) {$exist = 1; last;}
		}
		$self->[NAME] .= ",".$fea->name unless($exist);
		#$self->[NAME] = ($self->[NAME] eq $fea->name) ? $self->[NAME] : $self->[NAME].",".$fea->name;
	}
	else {
		$self->[NAME] = $fea->name; 
	}
	
	#update the start and end of the gene
	$self->[START] = $fea->start if(!$self->[START] || $self->[START] > $fea->start);
	$self->[END]   = $fea->end if(!$self->[END] || $self->[END] < $fea->end);
}

sub get_start {
	my ($self, $pos, $ext) = @_;
	my $rtn = $pos;
	#print STDERR "=== get_start: ", $pos, " ===\n";
	foreach my $t (@{$self->[TRANSCRIPTS]}) {
		my $tmp = $t->get_start($pos, $ext);
		print STDERR "get_start: ", $tmp, "\n";
		$rtn = $tmp if($tmp < $rtn);
	}
	return $rtn;
}

sub get_end {
	my ($self, $pos, $ext) = @_;
	my $rtn = $pos;
	#print STDERR "get_end: ", $pos, "\n";
	foreach my $t (@{$self->[TRANSCRIPTS]}) {
		my $tmp = $t->get_end($pos, $ext);
		#print STDERR "get_end: ", $tmp, "\n";
		$rtn = $tmp if($tmp > $rtn);
	}
	return $rtn;
}

sub get_feature {
	my ($self, $chr, $pos, $strand) = @_;

	my $feature;
	my ($score,$tmp_score) = (0,0);
	foreach my $t ( @{$self->transcripts} ) {
		my ($tmp_feature, $tmp_score);
		$tmp_feature = $t->get_feature($chr, $pos, $strand);
		if(!$tmp_feature){
			$tmp_score = 0.1;
			$tmp_feature = "intergenic";
			#print STDERR "tmp_feature = intergenic\n";
		}
		$tmp_score = 1 if($tmp_feature eq 'coding');
		$tmp_score = 0.8 if($tmp_feature =~ /utr/);
		$tmp_score = 0.5 if($tmp_feature eq 'intron');
		$tmp_score = 0.1 if($tmp_feature eq 'intergenic');
		if($tmp_score > $score) {
			$score = $tmp_score;
			$feature = $tmp_feature;
		}
	}
	return $feature;
}

sub get_mRNA_length {
	my $self = shift;
	my $rtn = 0;
	foreach my $t ( @{$self->transcripts} ) {
	#foreach my $t ( @{$g->[TRANSCRIPTS]} ) {
		#print join("\t", $t->name, $t->refseq_id, $t->chr, $t->start, $t->end, $t->strand, $t->cds_start, $t->cds_end, $t->type), "\n";
		#print "length: ",$t->get_length, "\n";
		$rtn = $t->get_length if($rtn < $t->get_length);
		#last;
	}
	return $rtn;
}

sub get_cds_length {
	my $self = shift;
	my $rtn = 0;
	foreach my $t (@{$self->[TRANSCRIPTS]}) {
		next if($t->refseq_id =~ m/NR/);
		#print STDERR "\nname\t".$t->name,"\n";
		#print STDERR "refseq_id\t".$t->refseq_id,"\n";
		#print STDERR "chr\t".$t->chr,"\n";
		#print STDERR "start\t".$t->start,"\n";
		#print STDERR "end\t".$t->end,"\n";
		#print STDERR "strand\t".$t->strand,"\n";
		#print STDERR "cds_start\t".$t->cds_start,"\n";
		#print STDERR "cds_end\t".$t->cds_end,"\n";
		#print STDERR "type\t".$t->type,"\n";
		my $cds_len = abs($t->cds_end - $t->cds_start);
		$rtn = $cds_len if($rtn < $cds_len);
	}
	return $rtn;
}


sub overlap {
	my ($self, $fea) = @_;

	if(ref($fea) eq 'ARRAY') { 
		foreach my $t ( @{$self->[TRANSCRIPTS]} ) {
			return 1 if($t->overlap($fea));
		}
	}
	elsif($fea->isa('Transcript')) {
		return if($fea->strand &&  $self->[STRAND] ne $fea->strand );
		return if($fea->chr && $self->[CHR] ne $fea->chr) ;
		#return if($fea->type && $self->[TYPE] ne $fea->type);
		foreach my $e ( @{$fea->exons}) {
			foreach my $t ( @{$self->[TRANSCRIPTS]} ) {
				return 1 if($t->overlap($e));
			}
		}
	}
	else {
		croak "Not implemented overlap";
	}
	return 0;
}

1;
