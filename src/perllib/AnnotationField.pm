package AnnotationField;
# standardization/QC for catchall annotation fields
# MNE 9/2014

use strict;
use Configurable;
use Exporter;

@AnnotationField::ISA = qw(Configurable Exporter);
@AnnotationField::EXPORT_OK = qw();

use MethodMaker qw(
        annotation_delimiter
        annotation_tag_value_delimiter
	annotation_value_delimiter

data
unique
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->annotation_delimiter(";");
  # main separator for annotation entries.  Use ";" as values may
  # contain arrays.
  $self->annotation_tag_value_delimiter("=");
  # separates annotation tags from values
  $self->annotation_value_delimiter(",");
  # delimiter for when tag value is an array
  $self->data([]);
  $self->configure(%options);
  return $self;
}

sub add {
  #
  # TO DO: check values for illegal use of reserved delimiter chars
  #
  my ($self, $key, $thing) = @_;
  my $annotation_listref = $self->data || die;
  my $value;
  if (defined $thing) {
    my $type = ref $thing;
    if ($type eq "ARRAY") {
      $value = join $self->annotation_value_delimiter(), @{$thing};
    } elsif ($type) {
      die;
    } else {
      $value = $thing;
    }
  }
  if (defined $value) {
    # tag=value
    push @{$annotation_listref}, join $self->annotation_tag_value_delimiter(), $key, $value;
  } else {
    # tag only
    push @{$annotation_listref}, $key;
  }
}

sub get_field {
  my ($self) = @_;
  my $data = $self->data || die;
  my @out;
  if ($self->unique) { 
    my %saw;
    foreach (@{$data}) {
      unless ($saw{$_}) {
	push @out, $_;
	$saw{$_} = 1;
      }
    }
  } else {
    @out = @{$data};
  }

  return join $self->annotation_delimiter(), @out;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
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
