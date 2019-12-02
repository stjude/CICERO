package CacheManager;
# manage a small hash cache
# MNE 11/2018
# TO DO:
# - prune by oldest access
# - prune by smallest hit count?
# - total reset?

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@CacheManager::ISA = qw(Configurable Exporter);
@CacheManager::EXPORT_OK = qw();

use MethodMaker qw(
	cache
	cache_stats
	cache_last_access
	cache_limit
	request_count
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->reset();
  return $self;
}

sub reset {
  my ($self) = @_;
  $self->cache({});
  $self->cache_stats({});
  $self->cache_last_access({});
  $self->request_count(0);
}

sub track_and_prune {
  my ($self, $key) = @_;
  my $rn = $self->request_count() + 1;
  $self->request_count($rn);
  my $cache_last_access = $self->cache_last_access();
  $cache_last_access->{$key} = $rn;

  my $limit = $self->cache_limit || die "-cache_limit";
  my $cache = $self->cache();
  if (scalar keys %{$cache} > $limit) {
    my @sorted = sort {$cache_last_access->{$a} <=> $cache_last_access->{$b}} keys %{$cache};
    my $key_prune = $sorted[0];
#    printf STDERR "pruning %s\n", $key_prune;

    delete $cache->{$key_prune};
    delete $cache_last_access->{$key_prune};
  }
}

sub get_count {
  my ($self) = @_;
  return scalar keys %{$self->cache};
}

sub get_keys {
  my ($self) = @_;
  return keys %{$self->cache};
}

sub get {
  my ($self, $key) = @_;
  my $hit = $self->cache()->{$key};
  my $tag = $hit ? "hits" : "misses";
  $self->cache_stats()->{$tag}++;
  $self->track_and_prune($key);
  return $hit;
}

sub put {
  my ($self, $key, $value) = @_;
  $self->track_and_prune($key);
  $self->cache()->{$key} = $value;
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
=head1 LICENCE AND COPYRIGHT
Copyright 2019 St. Jude Children's Research Hospital 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
