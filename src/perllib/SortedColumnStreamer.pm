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
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->die_if_unsorted(1);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  $self->track_saw({});
  unless ($self->df) {
    my $df = new DelimitedFile("-file" => ($self->file() || die "-file"),
			       "-headers" => 1,
			     );
    $self->df($df)
  }
}

sub next_set {
  my ($self, %options) = @_;
  my $df = $self->df() || die;
  my $field = $self->field() || die "-field";
  my @queue;
  my $old_queue = $self->queue();

  while (1) {
    my $row;
    if ($old_queue) {
      $row = $old_queue;
      $old_queue = undef;
      $self->queue(0);
    } else {
      last unless $row = $df->get_hash();
    }

    die "no field $field" unless exists $row->{$field};
    my $value = $row->{$field};
    if (@queue and $value ne $queue[0]->{$field}) {
      # new value
      $self->queue($row);
      # TO DO: add out-of-order tracking here
      last;
    } else {
      # same value or start of file
      push @queue, $row;
    }
  }

  if (@queue) {
    if ($self->die_if_unsorted()) {
      # sanity check file ordering
      my $v_this = $queue[0]->{$field};
      my $saw = $self->track_saw || die;
      confess sprintf "ERROR: input file is not sorted by field %s, already processed row set for %s", $self->field, $v_this if $saw->{$v_this};
      $saw->{$v_this} = 1;
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
