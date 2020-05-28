package Transcript;
use strict;
use Carp;
use Data::Dumper;

# we are going to use a light weight Transcript model here 
my @Transcript_slots;
BEGIN {
	@Transcript_slots = qw(NAME REFSEQ_ID CHR START END STRAND CDS_START CDS_END EXONS TYPE);
}
use enum @Transcript_slots;

my %attribute = (
    	name         	=> NAME,
    	refseq_id    	=> REFSEQ_ID,
	chr	 	=> CHR,
    	start        	=> START,
    	end          	=> END,
	strand		 => STRAND,
    	cds_start    	=> CDS_START,
	cds_end      	=> CDS_END,
	exons        	=> EXONS,
	type		=> TYPE,
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

    if (@_) {
	my %arg = @_;
        $obj->[NAME]      = $arg{-NAME}      if($arg{-NAME});
	$obj->[REFSEQ_ID] = $arg{-REFSEQ_ID} if($arg{-REFSEQ_ID});
	$obj->[CHR]       = $arg{-CHR}       if($arg{-CHR});
	$obj->[START]     = $arg{-START}     if($arg{-START});
	$obj->[END]       = $arg{-END}       if($arg{-END});
	$obj->[STRAND]    = $arg{-STRAND}    if($arg{-STRAND});
	$obj->[CDS_START] = $arg{-CDS_START} if($arg{-CDS_START});
	$obj->[CDS_END]   = $arg{-CDS_END}   if($arg{-CDS_END});
	$obj->[EXONS]     = $arg{-EXONS}     if($arg{-EXONS});
	$obj->[TYPE]      = $arg{-TYPE}      if($arg{-TYPE});
    }
    return bless $obj, $class;
}

sub get_feature {

	my ($self, $chr, $bp, $strand) = @_;
	#print STDERR "self $self, chr $chr, bp $bp, strand $strand\n";

	return 'intergenic' if($bp > $self->end && $strand eq '+');
	return 'intergenic' if($bp < $self->start && $strand eq '-');
	my @tmp = @{$self->[EXONS]};
	#print STDERR "--strand $strand\n";
	foreach my $e (@tmp) {
		if($strand eq '+' && $bp > $e->[0] && $bp < $e->[1]){
			if ($self->strand > 0){
				return '5utr' if($bp <= $self->cds_start);
				return '3utr' if($bp > $self->cds_end);
			}
			if ($self->strand < 0){
				return '3utr' if($bp <= $self->cds_start);
				return '5utr' if($bp > $self->cds_end);
			}
			return 'coding';
		}
		if($strand eq '-' && $bp > $e->[0] && $bp < $e->[1]){
			if ($self->strand > 0){
				return '5utr' if($bp <= $self->cds_start);
				return '3utr' if($bp > $self->cds_end);
			}
			if ($self->strand < 0){
				return '3utr' if($bp <= $self->cds_start);
				return '5utr' if($bp > $self->cds_end);
			}
			return 'coding';
		}
	}
	return 'intron';
}

sub get_start {
	my ($self, $pos, $ext) = @_;
	my @tmp;
	foreach my $e( @{$self->[EXONS]} ) {
		if($e->[1] < $pos) {
			push @tmp, $e;
			next;
		}
		last;
	}
	my $len = 0;
	while(scalar @tmp > 0) {
		my $e = pop @tmp;
		if($e->[1] >= $pos) {
			my $l = $pos - $e->[0];
			if($l + $len < $ext) {
				$len = $l;
				next;
			}
			return $pos - $ext;
		}
		if($e->[1] - $e->[0] + 1 + $len < $ext) {
			$len += ($e->[1] - $e->[0] + 1);
			next;
		}
		return ($e->[1] - $ext + $len);
	}
	return $self->start;
}

sub get_end {
	my ($self, $pos, $ext) = @_;
	my @tmp = @{$self->[EXONS]};
	my $len = 0;
	while(scalar @tmp > 0) {
		my $e = shift @tmp;
		next if($e->[1] < $pos);
		if($e->[0] <= $pos ) {
			return $pos + $ext if($e->[1] - $pos >= $ext);
			$len = $e->[1] - $pos;
			next;
		}
		if($e->[1] - $e->[0] + 1 + $len < $ext) {
			$len += ($e->[1] - $e->[0] + 1);
			next;
		}
		return ($e->[0] + $ext - $len);
	}
	return $self->end;
}


sub get_length{

	my $self=shift;
	my @tmp = @{$self->[EXONS]};
	my $len = 0;
	foreach my $e (@tmp) {
		$len += ($e->[1] - $e->[0] + 1);
	}
	return $len;
}

sub overlap {
	my $self = shift;
	my $range = shift;
	croak "Range must be a ref of array" unless(ref($range) eq 'ARRAY');

	foreach my $e ( @{$self->[EXONS]} ) {
		return 1 if($e->[0] <= $range->[1] && $e->[1] >= $range->[0]);
	}
	return;
}

1;
