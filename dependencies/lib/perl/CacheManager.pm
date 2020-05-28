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

pruned
verbose
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
  my ($self, $key, $type) = @_;
  my $rn = $self->request_count() + 1;
  $self->request_count($rn);
  my $cache_last_access = $self->cache_last_access();
  $cache_last_access->{$key} = $rn;

  $self->cache_debug($key, $type) if $self->verbose();

  my $limit = $self->cache_limit || die "-cache_limit";
  my $cache = $self->cache();
  if (scalar keys %{$cache} > $limit) {
    my @sorted = sort {$cache_last_access->{$a} <=> $cache_last_access->{$b}} keys %{$cache};
    my ($key_prune, @keep) = @sorted;
    printf STDERR "CacheManager: pruning %s, keeping %s\n", $key_prune, join ",", @keep if $self->verbose();

    delete $cache->{$key_prune};
    delete $cache_last_access->{$key_prune};

    $self->pruned($key_prune);
  } else {
    $self->pruned(0);
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
  $self->track_and_prune($key, "get");
  return $hit;
}

sub put {
  my ($self, $key, $value) = @_;
#  $self->track_and_prune($key, "put");
  $self->cache()->{$key} = $value;
  $self->track_and_prune($key, "put");
}

sub cache_debug {
  my ($self, $key, $type) = @_;
  my $stats = $self->cache_stats() || die;
  printf STDERR "cache debug for %s %s\n", $type, $key;
  printf STDERR "  summary:\n";
  foreach my $k (keys %{$stats}) {
    printf STDERR "    %s: %d\n", $k, $stats->{$k};
  }

  my $cache_last_access = $self->cache_last_access();
  printf STDERR "  last access:\n";
  foreach my $k (sort {$cache_last_access->{$b} <=> $cache_last_access->{$a}} keys %{$cache_last_access}) {
    printf STDERR "    %s: %s\n", $k, $cache_last_access->{$k};
  }
}


1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
