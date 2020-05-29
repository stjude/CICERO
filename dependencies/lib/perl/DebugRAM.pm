package DebugRAM;
# simple RAM usage debugging tool
# function call rather than OO to simplify code insertion
# MNE 2/2020

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@DebugRAM::ISA = qw(Configurable Exporter);
@DebugRAM::EXPORT_OK = qw(
debug_ram
);

sub debug_ram {
  # EXPORTED
  my ($label, $type) = @_;
  die "usage: label, start|end" unless $label and $type;

  # get RAM usage:
  my $cmd = sprintf 'ps u';
  open(PSTMP, sprintf '/bin/ps u %d|', $$) || die;
  my $hl = <PSTMP>;
  chomp $hl;
  my $dl = <PSTMP>;
  close PSTMP;

  my @h = split /\s+/, $hl;
  my @d = split /\s+/, $dl;
  # HACK: breaks for command portion, but should work for earlier fields
  my %info;
  @info{@h} = @d;

  if ($type eq "start") {
    $DebugRAM::TRACK{$label} = \%info;
  } elsif ($type eq "end") {
    my $start = $DebugRAM::TRACK{$label} || die "no tracking info for $label";
    printf STDERR "RAM change for %s: RSS:%d VSZ:%d\n",
      $label,
	$info{RSS} - $start->{RSS},
	  $info{VSZ} - $start->{VSZ};

  } else {
    die "unhandled type $type";
  }
}



1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
