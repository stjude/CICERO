package MethodMaker;
# generate simple variable accessor methods for a package
# MNE

use strict;

sub import {
  my ($package, @names) = @_;

  my $buffer = sprintf "\{ package %s;\n", scalar caller();
  foreach (@names) {
    $buffer .= sprintf '
sub %s {
  return(defined $_[1] ? $_[0]->{"%s"} = $_[1] : $_[0]->{"%s"});
}
', $_, $_, $_;
  }
  $buffer .= "\n\}\n";

#print STDERR $buffer;
  eval $buffer;
  
  die "yikes: $buffer" if $@;
}

1;

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
