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
