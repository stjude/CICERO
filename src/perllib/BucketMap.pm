package BucketMap;
# approximate hashing

use strict;
use Configurable;

use Carp qw(confess);

@BucketMap::ISA = qw(Configurable Exporter);
@BucketMap::EXPORT_OK = qw();

use MethodMaker qw(
		    chunk
		    hash
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->hash({});
  $self->configure(%options);
  die "need -chunk" unless $self->chunk();
  return $self;
}

sub get_key {
  my ($self, $raw) = @_;
  return int($raw / $self->chunk);
}

sub add_keyed {
  my ($self, $key, $value) = @_;
  $self->hash->{$key}->{$value} = $value;
  # quasi-unique
}

sub add_range {
  my ($self, %options) = @_;
  my $start = $options{"-start"};
  die "-start" unless defined $start;
  my $key_start = $self->get_key($start);
  my $key_end = $self->get_key($options{"-end"} || die);
  for (my $i = $key_start; $i <= $key_end; $i++) {
    $self->add_keyed($i, $options{"-value"} || die);
  }
}

sub fuzzy_find {
  my ($self, %options) = @_;
  my @keys;
  if (my $site = $options{"-site"}) {
    # single
    my $key = $self->get_key($site);
    push @keys, $key if defined $key;
  } elsif ($options{"-start"}) {
    my $key_start = $self->get_key($options{"-start"} || die);
    my $key_end = $self->get_key($options{"-end"} || die);
    @keys = ($key_start .. $key_end);
  } else {
    die;
  }

  my @hits;
  foreach my $key (@keys) {
    my $hits = $self->hash->{$key};
    push @hits, values %{$hits} if $hits;
  }

  return \@hits;
}

1;

=head1 LICENCE AND COPYRIGHT
Copyright 2019 St. Jude Children's Research Hospital 

Licensed under a modified version of the Apache License, Version 2.0
(the "License") for academic research use only; you may not use this
file except in compliance with the License. To inquire about commercial
use, please contact the St. Jude Office of Technology Licensing at
scott.elmer@stjude.org.

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
