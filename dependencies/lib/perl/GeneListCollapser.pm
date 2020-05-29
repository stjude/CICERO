package GeneListCollapser;
# attempts to prune a gene list to its primary component
# used in COSMIC parsing, etc.
# MNE 8/2013

use strict;
use Configurable;
use Exporter;

@GeneListCollapser::ISA = qw(Configurable Exporter);
@GeneListCollapser::EXPORT_OK = qw();

use MethodMaker qw(
single_entry_mode
verbose
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub collapse {
  my ($self, %options) = @_;

  my $genes = $options{"-hash"} || die "-hash";

  my @before = sort keys %{$genes};
  
  $self->filter_component_gene($genes);
  $self->filter_class($genes, "-orf" => 1);
  $self->filter_class($genes, "-loc" => 1);
  $self->filter_class($genes, "-ensg" => 1);

  if ($self->single_entry_mode() || $options{"-single"}) {
    my @g = sort keys %{$genes};
    if (@g > 1) {
      printf STDERR "still multiple genes (%s), stripping to ", join ",", @g if $self->verbose();
      for (my $i = 1; $i < @g; $i++) {
	delete $genes->{$g[$i]};
      }
      printf STDERR "%s\n", keys %{$genes} if $self->verbose();
    }
  }

  if (@before < scalar keys %{$genes}) {
    printf STDERR "GeneListCollapser: collapsed %s to %s\n",
    join(",", @before), join(",", sort keys %{$genes});
  }



  return $genes;
}

sub filter_component_gene {
  # given e.g. C11orf52,HSPB2-C11ORF52, delete HSPB2-C11ORF52
  my ($self, $genes) = @_;
  my @all_genes = keys %{$genes};
  my @reject;
  foreach my $g (@all_genes) {
    next if $g =~ /\-/;
    foreach my $other (@all_genes) {
      next if $g eq $other;
      next unless $other =~ /\-/;
      my $reject;
      foreach my $component (split /\-/, $other) {
	push @reject, $other if lc($component) eq lc($g);
      }
    }
  }
  if (@reject) {
    my @before = sort keys %{$genes};
    delete @{$genes}{@reject};
    printf STDERR "gene cleanup (component): before:%s after:%s\n",
    join(",", @before),
    join(",", sort keys %{$genes}) if $self->verbose();
  }
}

sub filter_class {
  # given C4orf39,TRIM61, delete C4orf39
  # given LOC100128398,ZNF606, delete LOC100128398
  my ($self, $genes, %options) = @_;
  my @all_genes = keys %{$genes};
  my @icky;
  my $label;
  if ($options{"-orf"}) {
    $label = "ORF";
    @icky = grep {/^\w\d+orf\d+$/i} @all_genes;
  } elsif ($options{"-loc"}) {
    $label = "LOC";
    @icky = grep {/^LOC/} @all_genes;
  } elsif ($options{"-ensg"}) {
    $label = "ENSG";
    @icky = grep {/^ENSG/} @all_genes;
  } else {
    die;
  }

  if (@icky and @icky < @all_genes) {
    my @before = sort keys %{$genes};
    delete @{$genes}{@icky};
    printf STDERR "gene cleanup (%s): before:%s after:%s\n",
    $label,
    join(",", @before),
    join(",", sort keys %{$genes}) if $self->verbose();
  }
}



1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
