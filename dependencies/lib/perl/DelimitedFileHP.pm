package DelimitedFileHP;
# higher-performance delimited file parser:
# avoid parsing into hashes by using arrayref
# MNE 9/2017

use strict;
use Exporter;
use Carp qw(confess);
use File::Basename;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use FileHandle;
use WorkingFile;
use FileUtils qw(universal_open);

@DelimitedFileHP::ISA = qw(Configurable Exporter);
@DelimitedFileHP::EXPORT_OK = qw();

use MethodMaker qw(
	file
delimiter
fh_in
fh_out
wf_out

headers
headers_extra
headers_out

h2i
current_line
header_exclude_regexp
header_exclude_field
header_indices_wanted

header_out_rename

query_indexes
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->delimiter("\t");
  $self->header_exclude_field({});
  $self->header_exclude_regexp([]);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  unless ($self->fh_in()) {
    my $fh_in = universal_open($self->file() || die "-file");
    $self->fh_in($fh_in);
  }
  my $headers = $self->next_row();
  $self->headers($headers);
  my $i = 0;
  my %h2i;
  for (my $i = 0; $i < @{$headers}; $i++) {
    my $h = $headers->[$i];
    die "duplicate header $h" if exists $h2i{$h};
    $h2i{$h} = $i;
  }
  $self->h2i(\%h2i);
}

sub next_row {
  my ($self) = @_;
  my $fh = $self->fh_in();
  if (my $line = <$fh>) {
    chomp $line;
    my $l = [ split $self->delimiter(), $line, -1 ];
    return $self->current_line($l);
  } else {
    $self->current_line(undef);
    return undef;
  }
}

sub get_hash {
  # convert current line to hash.  SLOW, so defeats main purpose of module
  my ($self) = @_;
  my %h;
  @h{@{$self->headers}} = @{$self->current_line};
  return \%h;
}

sub get_value {
  # get value of named field from given row
  my ($self, $field, $line) = @_;
  $line = $self->current_line() unless $line;
  my $idx = $self->h2i()->{$field};
  confess(sprintf "can't find a column named %s: available fields=%s", $field, join ",", sort keys %{$self->h2i()}) unless defined $idx;
  return $line->[$idx];
}

sub write_init {
  # prepare to write output
  my ($self, %options) = @_;
  my $outfile = $options{"-file"};
  unless ($outfile) {
    my $suffix = $options{"-suffix"} || ".filtered.tab";
    die "$suffix needs leading period" unless $suffix =~ /^\./;
    $outfile = basename($self->file) . $suffix;
  }
  my $fh;
  if ($options{"-stdout"} or $outfile eq "-") {
    $fh = *main::STDOUT;
  }  else {
    my $wf = new WorkingFile(
			     $outfile,
			     "-compress" => $options{"-compress"},
			    );
    # don't just pass through %options because it may contain params
    # unknown to WorkingFile and cause method not found crashes
    $self->wf_out($wf);
    $fh = $wf->output_filehandle();
  }

  my $hre = $self->header_exclude_regexp();
  my $hef = $self->header_exclude_field();
  my @headers;
  my $headers = $self->headers;
  if (@{$hre} or %{$hef}) {
    my @header_indices_wanted;
    for (my $i = 0; $i < @{$headers}; $i++) {
      my $usable = 1;
      foreach my $re (@{$hre}) {
	if ($headers->[$i] =~ /$re/) {
	  $usable = 0;
	  last;
	}
      }
      $usable = 0 if $hef and $hef->{$headers->[$i]};

      if ($usable) {
	push @headers, $headers->[$i];
	push @header_indices_wanted, $i;
      }
    }
    $self->header_indices_wanted(\@header_indices_wanted);
  } else {
    @headers = @{$self->headers};
  }

  if (my $rename = $self->header_out_rename) {
    foreach (@headers) {
      $_ = $rename->{$_} if $rename->{$_};
    }
  }

  push @headers, @{$self->headers_extra()} if $self->headers_extra();
  printf $fh "%s\n", join $self->delimiter(), @headers;
  $self->headers_out(\@headers);
  $self->fh_out($fh);
}

sub write_row {
  my ($self, %options) = @_;
  my $fh = $self->fh_out();
  my $line = $options{"-line"} || $self->current_line();
  my @extra;
  if (my $hex = $self->headers_extra) {
    my $extra = $options{"-extra"} || die "-extra";
    foreach my $h (@{$hex}) {
      dump_die($extra, "where is $h") unless exists $extra->{$h};
      push @extra, $extra->{$h};
    }
  }

  if (my $hiw = $self->header_indices_wanted) {
    printf $fh "%s\n", join $self->delimiter, @{$line}[@{$hiw}], @extra;
  } else {
    printf $fh "%s\n", join $self->delimiter, @{$line}, @extra;
  }
}

sub write_finish {
  # finish and close writer
  my ($self) = @_;
  $self->wf_out->finish() if $self->wf_out();
  # close if writing to file
}

sub add_header_exclude_regexp {
  my ($self, $regexp) = @_;
  push @{$self->header_exclude_regexp}, $regexp;
}

sub add_header_exclude_field {
  my ($self, $field) = @_;
  $self->header_exclude_field->{$field} = 1;
}

sub has_header {
  my ($self, $field) = @_;
  return exists $self->h2i->{$field};
}

sub prepare_query {
  # build index lookup list for a set of fields (SQL prepare-ish)
  my ($self, %options) = @_;
  my $fields = $options{"-fields"} || die "-fields";
  my @idx;
  foreach my $f (@{$fields}) {
    my $idx = $self->h2i()->{$f};
    unless (defined $idx) {
      printf STDERR "ERROR: no field %s, available:\n";
      foreach $_ (sort keys (%{$self->h2i()})) {
	printf STDERR "  %s\n", $_;
      }
      die;
    }
      
    push @idx, $idx;
  }
  $self->query_indexes(\@idx);
}

sub get_query {
  my ($self) = @_;
  return [ @{$self->current_line()}[@{$self->query_indexes()}] ];
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
