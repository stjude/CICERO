package GeneModel;
use strict;
use Carp;
use Data::Dumper;
use Tree::Interval;
use Tree::DAG_Node;
use English;
use Gene;

sub new {
    my $class = shift;
    my $obj = {};
	$obj->{FMODEL} = {};
	$obj->{RMODEL} = {};
    if (@_) {
		my %arg = @_;
		$obj->{FMODEL} = $arg{-FMODEL} if($arg{-FMODEL});
		$obj->{RMODEL} = $arg{-RMODEL} if($arg{-RMODEL});
    }
    return bless $obj, $class;
}

sub get_all_chr {
	my $self = shift;
	return keys %{$self->{FMODEL}};
}

sub add_gene {
	my ($self, $gene) = @_;
	if($gene->strand eq '+') {
		$self->{FMODEL}->{$gene->chr}->insert([$gene->start, $gene->end], $gene);
	}
	else {
		$self->{FMODEL}->{$gene->chr}->insert([$gene->start, $gene->end], $gene);
	}
}

sub remove_gene {
	my $self = shift;
	my $gene = shift;
	if($gene->strand eq '+') {
		$self->{FMODEL}->{$gene->chr}->remove([$gene->start, $gene->end], $gene);
	}
	else {
		$self->{RMODEL}->{$gene->chr}->remove([$gene->start, $gene->end], $gene);
	}
}

sub n_genes {
	my $self = shift;
	my $n = 0;
	foreach my $c (keys %{$self->{FMODEL}}) {
		$n += $self->{FMODEL}->{$c}->size;
	}
	foreach my $c (keys %{$self->{RMODEL}}) {
		$n += $self->{RMODEL}->{$c}->size;
	}
	return $n;
}

sub sub_model {
	my ($self, $chr, $strand) = @_;
	if($strand eq '+' ) {
		return $self->{FMODEL}->{$chr};
	}
	else {
		return $self->{RMODEL}->{$chr};
	}
}
sub look_up_gene {
	my ($self, $chr, $gene) = @_;
	if($gene->isa('Gene')) {
		if($gene->strand eq '+') {
			return $self->{FMODEL}->{$chr}->lookup([$gene->start, $gene->end]);
		}
		else {
			return $self->{RMODEL}->{$chr}->lookup([$gene->start, $gene->end]);
		}
	}
	elsif(ref($gene) eq 'ARRAY') {
		return (
			$self->{FMODEL}->{$chr}->lookup([$gene->start, $gene->end]),
			$self->{RMODEL}->{$chr}->lookup([$gene->start, $gene->end]) );
	}
	else {
		croak "Not implemented look_up_gene for $gene";
	}
}

sub overlap_genes {
	my ($self, $chr,  $fea, $ort) = @_;
	my @genes;

	if(ref($fea) eq 'ARRAY') {
		@genes = $self->_overlap_array($chr, $fea, $ort);
	}
	elsif($fea->isa('Transcript')) {
		foreach my $e (@{$fea->exons}) {
			push @genes, $self->_overlap_array($chr, $e,  
				$fea->strand ? $fea->strand : undef);
		}
	}
	else {
		croak "Not implemented overlap for $fea";
	}
	my @rtn;
	foreach my $g (@genes) {
		push @rtn, $g if($g->overlap($fea));
	}
	return @rtn;
}

sub _overlap_array {
	my($self, $chr, $a, $ort) = @_;
	if($ort) {
		return $ort == '+' ? $self->{FMODEL}->{$chr}->intersect($a) :
			$self->{RMODEL}->{$chr}->intersect($a);
	}
	else {
		return ($self->{FMODEL}->{$chr}->intersect($a), $self->{RMODEL}->{$chr}->intersect($a));
	}
}

sub model {
	my ($self, $chr, $strand, $val) = @_;
	croak "Please specify chr and strand for GeneModel" 
		unless $chr && $strand;
	my $model = $strand eq '+' ? $self->{FMODEL} : $self->{RMODEL};

	if($val) {
		croak "The model for each chrom must be a Tree::Interval" 
			unless $val->isa('Tree::Interval');
		$model->{$chr} = $val;
	}
	return $model->{$chr};
}

sub from_file {
	my ($self, $file, $format) = @_;
	croak "$format is not implemented yet!" unless uc($format) eq "REFFLAT";
	#print STDERR "Processing gene model file, please wait ... ";
	my $transcripts = _refFlat2Transcripts($file);
	my $model = _Transcript2GeneModel($transcripts);
	foreach my $c (keys %{$model}) {
		foreach my $s (keys %{$model->{$c}}) {
			if($s eq '+') {
				$self->{FMODEL}->{$c} = $model->{$c}{$s};
			}
			else {
				$self->{RMODEL}->{$c} = $model->{$c}{$s};
			}
		}
	}
	#print STDERR "Done.\n";
}

sub _refFlat2Transcripts {
	my $file = shift;
	my %transcripts;

	open my $FILE, "<$file" or croak "Can't open $file: $OS_ERROR";
	while( my $line = <$FILE> ) {
		next if($line =~ /^#/); #annotation line
		#next unless($line =~ 'PTX1');
		my ($name, $refseq_id, $chr, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, 
		          $exonStarts, $exonEnds) = split(/\s+/, $line);
		my $type = $cdsEnd == $cdsStart ? "NR" : "NM";
		my @starts = split /,/, $exonStarts;
		my @ends = split /,/, $exonEnds;
		if(scalar @starts != scalar @ends) {
			warn "Number of exon starts is not same as exon ends!";
			next;
		}
		$transcripts{$chr} = {} unless (exists $transcripts{$chr});
		$transcripts{$chr}->{$strand} = [] unless (exists $transcripts{$chr}->{$strand});
		my @exons;
		for( my $i = 0; $i < $exonCount; $i++) {
			push @exons, [$starts[$i], $ends[$i]];
		}

		#$obj->[CDS_END]   = $arg{-CDS_END}   if($arg{-CDS_END});
		my $t = Transcript->new(-START => $starts[0], -END => $ends[$#ends], -CHR => $chr, 
			-STRAND => $strand eq '+' ? 1 : -1, -NAME => $name, -REFSEQ_ID => $refseq_id,
			-EXONS => \@exons, -TYPE => $type, -CDS_START => $cdsStart, -CDS_END => $cdsEnd);
		push @{$transcripts{$chr}->{$strand}}, $t if($t->start && $t->end);
		#print STDERR join("\t", $t->name, $t->start, $t->end), "\n";
	}
	#sort the list of transcripts
	foreach my $c (keys %transcripts) {
		#print STDERR "chrom $c\n";
		foreach my $s (keys %{$transcripts{$c}}) {
			#print STDERR "strand $s\n";
			my @tmp = @{$transcripts{$c}->{$s}};
			@tmp = sort { $a->start <=> $b->start || $a->end <=> $b->end } @tmp;
			$transcripts{$c}->{$s} = \@tmp;
		}
	}
	return \%transcripts;
}

sub _Transcript2GeneModel {
	my $transcripts = shift;
	my %genemodel;

	foreach my $c (keys %{$transcripts}) {
		$genemodel{$c} = {};	
		foreach my $s (keys %{$transcripts->{$c}}) {
			my $tree = $genemodel{$c}->{$s} = Tree::Interval->new();
			foreach my $t (@{$transcripts->{$c}{$s}}) {
				my @ol = $tree->intersect([$t->start, $t->end]);
				if(scalar @ol == 0 ) {
					my $gene = Gene->new();
					$gene->add_transcript($t);
					$tree->insert([$t->start, $t->end], $gene);
					next;
				}
				my @real_ol;
				foreach my $g (@ol) {
					push @real_ol, $g->val if($g->val->overlap($t));
				}
				if(scalar @real_ol == 0) { #gene in gene?
#					print STDERR "Gene in Gene: ", $t->name, " in ", $ol[0]->val->name, " \n";
					my $gene = Gene->new();
					$gene->add_transcript($t);
					$tree->insert([$t->start, $t->end], $gene);
					next;
				}
				if(scalar @real_ol == 1 ) { # transcript belongs to gene
					my $g = $real_ol[0];
#					print STDERR "Same Gene group: ", $t->name, " in ", $real_ol[0]->name, " \n";
					if($t->start < $g->start || $t->end > $g->end) {
						$tree->delete([$g->start, $g->end], $g);
						$g->add_transcript($t);
						$tree->insert([$g->start, $g->end], $g);
						next;
					}
					$g->add_transcript($t);
					next;
				}
				# a transcript make more genes same
#				print STDERR "Join Gene group: ", $t->name, " join ", $real_ol[0]->name, "and", $real_ol[1]->name, " \n";
				my $gene = Gene->new();
				$gene->add_transcript($t);
				foreach my $g (@real_ol) {
					$tree->delete([$g->start, $g->end], $g);
					foreach my $tp (@{$g->transcripts}) {
						$gene->add_transcript($tp);
					}
				}
				$tree->insert([$gene->start, $gene->end], $gene);
			}
		}
	}
	return \%genemodel;
}

sub _print_tree {
	return;
	my $tree = shift;
	my $t = Tree::DAG_Node->lol_to_tree($tree->root->as_lol);
	local $, = "\n";
	print @{$t->draw_ascii_tree}, "\n";
}
1;
