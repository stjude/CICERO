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
