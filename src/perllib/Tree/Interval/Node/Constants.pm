package Tree::Interval::Node::Constants;

use strict;
use Carp;
use vars qw( $VERSION @EXPORT );

$VERSION = '0.1';

require Exporter;
*import = \&Exporter::import;

my @Node_slots;
my @Node_colors;

BEGIN { 
    @Node_slots  = qw(PARENT LEFT RIGHT COLOR KEY VAL MAX INTERVAL); 
    @Node_colors = qw(RED BLACK);
}

@EXPORT = (@Node_colors, map {"_$_"} @Node_slots);

use enum @Node_colors;
use enum @Node_slots;

# enum doesn't allow symbols to start with "_", but we want them 
foreach my $s (@Node_slots) {
    no strict 'refs';
    *{"_$s"} = \&$s;
    delete $Tree::RB::Node::Constants::{$s};
} 

1;
