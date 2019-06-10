package Tree::Interval::Node;

use strict;
use Carp;
use Tree::Interval::Node::Constants;
use vars qw( $VERSION @EXPORT_OK );

require Exporter;
*import    = \&Exporter::import;
@EXPORT_OK = qw[set_color color_of parent_of left_of right_of];

$VERSION = '0.1';

# key and interval is the same thing
my %attribute = (
    key    => _KEY,
    val    => _VAL,
    color  => _COLOR,
    parent => _PARENT,
    left   => _LEFT,
    right  => _RIGHT,
	max    => _MAX,
	interval => _INTERVAL,
);

#using an array instead of a hash for the node
sub _accessor {
    my $index = shift;
    return sub {
        my $self = shift;
		return undef unless $self;
        if (@_) {
          $self->[$index] = shift;
        }
        return $self->[$index];
    };
}

while(my($at, $idx) = each %attribute) {
    no strict 'refs';
    *$at = _accessor($idx);
}

sub new {
    my $class = shift;
    my $obj = [];

    if (@_) {
        $obj->[_KEY] = $obj->[_INTERVAL] = shift;
        $obj->[_VAL] = shift;
    }
    return bless $obj, $class;
}

sub left_most {
    my $self = shift;
    while ($self->[_LEFT]) {
        $self = $self->[_LEFT];
    }
    return $self;
}

sub right_most {
    my $self = shift;
    while ($self->[_RIGHT]) {
        $self = $self->[_RIGHT];
    }
    return $self;
}

#find left_most leaf
sub leaf {
    my $self = shift;
    while (my $any_child = $self->[_LEFT] || $self->[_RIGHT]) {
        $self = $any_child;
    }
    return $self;
}

sub successor {
    my $self = shift;
    if ($self->[_RIGHT]) {
        return $self->[_RIGHT]->left_most;
    }
    my $parent = $self->[_PARENT];
    while ($parent && $parent->[_RIGHT] && $self == $parent->[_RIGHT]) {
        $self = $parent;
        $parent = $parent->[_PARENT];
    }
    return $parent;
}

sub predecessor {
    my $self = shift;
    if ($self->[_LEFT]) {
        return $self->[_LEFT]->right_most;
    }
    my $parent = $self->[_PARENT];
    while ($parent && $parent->[_LEFT] && $self == $parent->[_LEFT]) {
        $self = $parent;
        $parent = $parent->[_PARENT];
    }
    return $parent;
}

sub as_lol {
    my $self = shift;
    my $node = shift || $self;
    my $aref;
    push @$aref,
         $node->[_LEFT]
           ? $self->as_lol($node->[_LEFT])
           : '*';
    push @$aref,
         $node->[_RIGHT]
           ? $self->as_lol($node->[_RIGHT])
           : '*';
    my $color = ($node->[_COLOR] == RED ? 'R' : 'B');
    no warnings 'uninitialized';
    push @$aref, "$color:[$node->[_KEY][0],$node->[_KEY][1]]:$node->[_MAX]";
    return $aref;
}

sub strip {
    my $self = shift;
    my $callback = shift;

    my $x = $self;
    while($x) {
        my $leaf = $x->leaf;
        $x = $leaf->[_PARENT];

        # detach $leaf from the (sub)tree
        no warnings "uninitialized";
        if($leaf == $x->[_LEFT]) {
            undef $x->[_LEFT];
        }
        else {
            undef $x->[_RIGHT];
        }
        undef $leaf->[_PARENT];
        if($callback) {
            $callback->($leaf);
        }

        if(!$x->[_LEFT] && !$x->[_RIGHT]) {
            $x = $x->[_PARENT];
        }
    }
}

sub DESTROY { $_[0]->strip; }

# Null aware accessors to assist with rebalancings during insertion and deletion
#
# A weird case of Java to the rescue!
# These are inspired by http://www.javaresearch.org/source/jdk142/java/util/TreeMap.java.html
# which was found via http://en.wikipedia.org/wiki/Red-black_tree#Implementations

# do wen need it? as we have accessors already
sub set_color {
    my ($node, $color) = @_;
    if($node) {
        $node->[_COLOR] = $color; 
    }
}

sub color_of {
    $_[0] ? $_[0]->[_COLOR] : BLACK;
}

sub parent_of {
    $_[0] ? $_[0]->[_PARENT] : undef;
}

sub left_of {
    $_[0] ? $_[0]->[_LEFT] : undef;
}

sub right_of {
    $_[0] ? $_[0]->[_RIGHT] : undef;
}

sub _overlap {
	my ($a, $b) = @_;
	return 1 if($a->[0] <= $b->[1] && $a->[1] >= $b->[0]);
	return undef;
}	

sub intersect {
	my $x = shift;
	my $interval = shift;
	return if(!$x);
#	print $x->val->name, "\t", $x->key->[0], "\t", $x->key->[1], "\n";
	my @rtn;
	if(_overlap($x->interval, $interval)) {
		push @rtn, $x;
	}
#	my $y = $x->parent;
	if($x->left && $x->left->max >= $interval->[0] ) { # && (!$y || _overlap($interval, [$y->interval->[0], $x->left->max]))) {
		push @rtn, $x->left->intersect($interval);
	}
	push @rtn, $x->right->intersect($interval) if($x->right && _overlap($interval, [$x->interval->[0], $x->right->max]));
	return @rtn;
}

1;

