package SortedColumnStreamer;
# given a DelimitedFile with grouped/sorted data, parse a block of
# rows until value changes, also die if that value is seen again
# later.  For example read all rows for a gene, die if name seen again
# (meaning an unsorted file, a fatal error in this use case).  Saves
# RAM when parsing huge but pre-sorted files.
# MNE 1/2020

use strict;
use Exporter;
use Carp qw(confess);

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@SortedColumnStreamer::ISA = qw(Configurable Exporter);
@SortedColumnStreamer::EXPORT_OK = qw();

use MethodMaker qw(
		    df
		    file

		    field
		    queue
		    die_if_unsorted
		    track_saw
hp
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->die_if_unsorted(1);
  $self->configure(%options);
  confess "specify -hp 0|1" unless defined $self->hp();
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  $self->track_saw({});
  unless ($self->df) {
    my $df;
    my $infile = $self->file() || die "-file";
    if ($self->hp) {
      $df = new DelimitedFileHP("-file" => $infile);
      $df->prepare_query("-fields" => [ $self->field() ]);
    } else {
      $df = new DelimitedFile("-file" => $infile,
			      "-headers" => 1,
			     );
    }
    $self->df($df);
  }
}

sub next_set {
  my ($self, %options) = @_;
  my $df = $self->df() || die;
  my $field = $self->field() || die "-field";
  my @queue;
  my $old_queue = $self->queue();
  my $hp_mode = $self->hp();

  while (1) {
    my $row;
    if ($old_queue) {
      $row = $old_queue;
      $old_queue = undef;
      $self->queue(0);
    } elsif ($hp_mode) {
      last unless $row = $df->next_row();
    } else {
      last unless $row = $df->get_hash();
    }

    my $value;
    if ($hp_mode) {
      ($value) = @{$df->get_query()};
    } else {
      die "no field $field" unless exists $row->{$field};
      $value = $row->{$field};
    }

    my $change;
    if (@queue) {
      my $queue_value = $hp_mode ? $df->get_value($field, $queue[0]) : $queue[0]->{$field};
      $change = 1 if $value ne $queue_value;
    }

    if ($change) {
      # new value
      $self->queue($row);
      last;
    } else {
      # same value or start of file
      push @queue, $row;
    }
  }

  if (@queue) {
    if ($self->die_if_unsorted()) {
      # sanity check file ordering
      my $queue_value = $hp_mode ? $df->get_value($field, $queue[0]) : $queue[0]->{$field};
      my $saw = $self->track_saw() || die;
      confess sprintf "ERROR: input file is not sorted by field %s, already processed row set for %s", $self->field, $queue_value if $saw->{$queue_value};
      $saw->{$queue_value} = 1;
    }
    return \@queue;
  } else {
    return undef;
  }

}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
