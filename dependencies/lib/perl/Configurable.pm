package Configurable;
# simple package configuration
# MNE

use strict;

#sub import {
#  my $code = sprintf '{package %s; no strict qw(vars); push @ISA, qw(Configurable);}', scalar caller();
#  eval $code;
#  print "$code\n";
#}

sub configure {
  my ($self, %options) = @_;
  while (my ($key, $value) = each %options) {
    $key =~ s/^-//;
    $self->$key($value);
  }
}

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub override_option {
  # override one or more object values with values in provided hashref
  my ($self, $option_ref, @options) = @_;
  my (@values, $option, $key);
  foreach $option (@options) {
    $key = '-' . $option;
    push @values, exists $option_ref->{$key} ? $option_ref->{$key} : $self->$option();
  }
  return @options == 1 ? $values[0] : @values;
}

1;
