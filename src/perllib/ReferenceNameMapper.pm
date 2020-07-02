package ReferenceNameMapper;
# try to connect different versions of reference sequence names

use strict;
use Configurable;
use Exporter;

@ReferenceNameMapper::ISA = qw(Configurable Exporter);
@ReferenceNameMapper::EXPORT_OK = qw();

use MethodMaker qw(
	db_names
	db_names_uc
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->db_names({});
  $self->db_names_uc({});
  $self->configure(%options);
  return $self;
}

sub add_name {
  my ($self, $name) = @_;
  $self->db_names->{$name} = 1;
  $self->db_names_uc->{uc($name)} = 1;
}

sub find_name {
  my ($self, $query) = @_;
  my $db_names = $self->db_names();
  my $db_names_uc = $self->db_names_uc();
  
  my @queries;
  push @queries, $query;
  my $check = $query;
  $check = uc($check);
  $check =~ s/^chr//i;
  if ($check eq "MT") {
    push @queries, "M";
  } elsif ($check eq "M") {
    push @queries, "MT";
  }

  my @try;
  foreach my $q (@queries) {
    push @try, $q;
    push @try, "chr" . $q;
    my $q2 = uc($query);
    $q2 =~ s/^chr//i;
    push @try, $q2;
  }

  my $result;

  foreach my $q (@try) {
    $result = $q if $db_names->{$q};
    last if defined $result;
  }

  unless (defined $result) {
    foreach my $q (@try) {
      $result = $q if $db_names_uc->{uc($q)};
      last if defined $result;
    }
  }

  return $result;
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
