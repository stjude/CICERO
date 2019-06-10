package Tree::Interval;

use strict;
use Carp;

use Tree::Interval::Node qw[set_color color_of parent_of left_of right_of];
use Tree::Interval::Node::Constants;
use Tree::DAG_Node;
use vars qw( $VERSION @EXPORT_OK );
$VERSION = 0.1;

require Exporter;
*import    = \&Exporter::import;
@EXPORT_OK = qw[LUEQUAL LUGTEQ LULTEQ LUGREAT LULESS LUNEXT LUPREV];

use enum qw{ LUEQUAL  LUGTEQ  LULTEQ LUGREAT LULESS LUNEXT LUPREV };

# object slots
use enum qw{ ROOT CMP SIZE };

sub new {
    my ($class, $cmp) = @_;
    my $obj = [];
    $obj->[SIZE] = 0;
    if($cmp) {
        ref $cmp eq 'CODE'
          or croak('Invalid arg: codref expected');
        $obj->[CMP] = $cmp;
    }
	else {
		# default compare
		$obj->[CMP] = sub { $_[0]->[0] <=> $_[1]->[0] || $_[0]->[1] <=> $_[1]->[1]};
	}
    return bless $obj => $class;
}


sub DESTROY { $_[0]->[ROOT]->DESTROY if $_[0]->[ROOT] }

sub root { $_[0]->[ROOT] }
sub size { $_[0]->[SIZE] }

sub left_most {
    my $self = shift;
    return undef unless $self->[ROOT];
    return $self->[ROOT]->left_most;
}

sub right_most {
    my $self = shift;
    return undef unless $self->[ROOT];
    return $self->[ROOT]->right_most;
}

# return all the intervals intersect with the interval
sub intersect {
	my $self = shift;
	my $interval = shift;
	return my @tmp  unless $self->[ROOT];
	return $self->[ROOT]->intersect($interval);
}

sub lookup {
    my $self = shift;
    my $key  = shift;
    defined $key
      or croak("Can't use undefined value as key");
    my $mode = shift || LUEQUAL;
    my $cmp = $self->[CMP];

    my $y;
    my $x = $self->[ROOT]
      or return;
    my $next_child;
    while($x) {
        $y = $x;
        if($cmp ? $cmp->($key, $x->[_KEY]) == 0
                : $key eq $x->[_KEY]) {
            # found it!
            if($mode == LUGREAT || $mode == LUNEXT) {
                $x = $x->successor;
            }
            elsif($mode == LULESS || $mode == LUPREV) {
                $x = $x->predecessor;
            }
            return wantarray
              ? ($x->[_VAL], $x)
              : $x->[_VAL];
        }
        if($cmp ? $cmp->($key, $x->[_KEY]) < 0
                : $key lt $x->[_KEY]) {
            $next_child = _LEFT;
        }
        else {
            $next_child = _RIGHT;
        }
        $x = $x->[$next_child];
    }
    # Didn't find it :(
    if($mode == LUGTEQ || $mode == LUGREAT) {
        if($next_child == _LEFT) {
            return wantarray ? ($y->[_VAL], $y) : $y->[_VAL];
        }
        else {
            my $next = $y->successor
              or return;
            return wantarray ? ($next->[_VAL], $next) : $next->[_VAL];
        }
    }
    elsif($mode == LULTEQ || $mode == LULESS) {
        if($next_child == _RIGHT) {
            return wantarray ? ($y->[_VAL], $y) : $y->[_VAL];
        }
        else {
            my $next = $y->predecessor
              or return;
            return wantarray ? ($next->[_VAL], $next) : $next->[_VAL];
        }
    }
    return;
}

sub insert {
    my $self = shift;
    my $key_or_node = shift;
    defined $key_or_node
      or croak("Can't use undefined value as key or node");
    my $val = shift;

    my $cmp = $self->[CMP];
    my $z = (ref $key_or_node eq 'Tree::Interval::Node')
              ? $key_or_node
              : Tree::Interval::Node->new($key_or_node => $val);

    my $y;
    my $x = $self->[ROOT];
    while($x) {
        $y = $x;
        # Handle case of inserting node with duplicate key.
        if($cmp ? $cmp->($z->[_KEY], $x->[_KEY]) == 0
                : $z->[_KEY] eq $x->[_KEY])
        {
			warn "The same key (range) is already in the tree
				it will be replaced!";
            my $old_val = $x->[_VAL];
            $x->[_VAL] = $z->[_VAL];
            return $old_val;
        }

        if($cmp ? $cmp->($z->[_KEY], $x->[_KEY]) < 0
                : $z->[_KEY] lt $x->[_KEY])
        {
            $x = $x->[_LEFT];
        }
        else {
            $x = $x->[_RIGHT];
        }
    }
    # insert new node
    $z->[_PARENT] = $y;
    if(not defined $y) {
        $self->[ROOT] = $z;
    }
    else {
        if($cmp ? $cmp->($z->[_KEY], $y->[_KEY]) < 0
                : $z->[_KEY] lt $y->[_KEY])
        {
            $y->[_LEFT] = $z;
        }
        else {
            $y->[_RIGHT] = $z;
        }
    }
	_update_max($z);
    $self->_fix_after_insertion($z);
    $self->[SIZE]++;
}

sub _fix_after_insertion {
    my $self = shift;
    my $x = shift or croak('Missing arg: node');

    $x->[_COLOR] = RED;
    while($x != $self->[ROOT] && $x->[_PARENT][_COLOR] == RED) {
        my ($child, $rotate1, $rotate2);
        if(($x->[_PARENT] || 0) == ($x->[_PARENT][_PARENT][_LEFT] || 0)) {
            ($child, $rotate1, $rotate2) = (_RIGHT, '_left_rotate', '_right_rotate');
        }
        else {
            ($child, $rotate1, $rotate2) = (_LEFT, '_right_rotate', '_left_rotate');
        }
        my $y = $x->[_PARENT][_PARENT][$child];

        if($y && $y->[_COLOR] == RED) {
            $x->[_PARENT][_COLOR] = BLACK;
            $y->[_COLOR] = BLACK;
            $x->[_PARENT][_PARENT][_COLOR] = RED;
            $x = $x->[_PARENT][_PARENT];
        }
        else {
            if($x == ($x->[_PARENT][$child] || 0)) {
                $x = $x->[_PARENT];
                $self->$rotate1($x);
            }
            $x->[_PARENT][_COLOR] = BLACK;
            $x->[_PARENT][_PARENT][_COLOR] = RED;
            $self->$rotate2($x->[_PARENT][_PARENT]);
        }
    }
    $self->[ROOT][_COLOR] = BLACK;
}

sub delete {
    my ($self, $key_or_node) = @_;
    defined $key_or_node
      or croak("Can't use undefined value as key or node");

    my $z = (ref $key_or_node eq 'Tree::Interval::Node')
              ? $key_or_node
              : ($self->lookup($key_or_node))[1];
    return unless $z;

    my $y;
    if($z->[_LEFT] && $z->[_RIGHT]) {
        # (Notes kindly provided by Christopher Gurnee)
        # When deleting a node 'z' which has two children from a binary search tree, the
        # typical algorithm is to delete the successor node 'y' instead (which is
        # guaranteed to have at most one child), and then to overwrite the key/values of
        # node 'z' (which is still in the tree) with the key/values (which we don't want
        # to lose) from the now-deleted successor node 'y'.

        # Since we need to return the deleted item, it's not good enough to overwrite the
        # key/values of node 'z' with those of node 'y'. Instead we swap them so we can
        # return the deleted values.

        $y = $z->predecessor;
        ($z->[_KEY], $y->[_KEY]) = ($y->[_KEY], $z->[_KEY]);
        ($z->[_VAL], $y->[_VAL]) = ($y->[_VAL], $z->[_VAL]);
        ($z->[_INTERVAL], $y->[_INTERVAL]) = ($y->[_INTERVAL], $z->[_INTERVAL]);
    }
    else {
        $y = $z;
    }
	
    # splice out $y
    my $x = $y->[_LEFT] || $y->[_RIGHT];
    if(defined $x) {
        $x->[_PARENT] = $y->[_PARENT];
        if(! defined $y->[_PARENT]) {
            $self->[ROOT] = $x;
        }
        elsif($y == $y->[_PARENT][_LEFT]) {
            $y->[_PARENT][_LEFT] = $x;
        }
        else {
            $y->[_PARENT][_RIGHT] = $x;
        }
        # Null out links so they are OK to use by _fix_after_deletion
        delete @{$y}[_PARENT, _LEFT, _RIGHT];
		_update_max($x);
        # Fix replacement
        if($y->[_COLOR] == BLACK) {
            $self->_fix_after_deletion($x);
        }
    }
    elsif(! defined $y->[_PARENT]) {
        # return if we are the only node
        delete $self->[ROOT];
    }
    else {
        # No children. Use self as phantom replacement and unlink
        if($y->[_COLOR] == BLACK) {
            $self->_fix_after_deletion($y);
        }
        if(defined $y->[_PARENT]) {
            no warnings 'uninitialized';
            if($y == $y->[_PARENT][_LEFT]) {
                delete $y->[_PARENT][_LEFT];
            }
            elsif($y == $y->[_PARENT][_RIGHT]) {
                delete $y->[_PARENT][_RIGHT];
            }
			my $tmp = $y->[_PARENT];
            delete $y->[_PARENT];
			_update_max($tmp);
        }
    }
    $self->[SIZE]--;
    return $y;
}

sub _fix_after_deletion {
    my $self = shift;
    my $x = shift or croak('Missing arg: node');

    while($x != $self->[ROOT] && color_of($x) == BLACK) {
        my ($child1, $child2, $rotate1, $rotate2);
        no warnings 'uninitialized';
        if($x == left_of(parent_of($x))) {
            ($child1,    $child2,   $rotate1,       $rotate2) =
            (\&right_of, \&left_of, '_left_rotate', '_right_rotate');
        }
        else {
            ($child1,   $child2,    $rotate1,        $rotate2) =
            (\&left_of, \&right_of, '_right_rotate', '_left_rotate');
        }
        use warnings;

        my $w = $child1->(parent_of($x));
        if(color_of($w) == RED) {
            set_color($w, BLACK);
            set_color(parent_of($x), RED);
            $self->$rotate1(parent_of($x));
            $w = right_of(parent_of($x));
        }
        if(color_of($child2->($w)) == BLACK &&
           color_of($child1->($w)) == BLACK) {
            set_color($w, RED);
            $x = parent_of($x);
        }
        else {
            if(color_of($child1->($w)) == BLACK) {
                set_color($child2->($w), BLACK);
                set_color($w, RED);
                $self->$rotate2($w);
                $w = $child1->(parent_of($x));
            }
            set_color($w, color_of(parent_of($x)));
            set_color(parent_of($x), BLACK);
            set_color($child1->($w), BLACK);
            $self->$rotate1(parent_of($x));
            $x = $self->[ROOT];
        }
    }
    set_color($x, BLACK);
}


sub _max3 {
	my ($a, $b, $c) = @_;
	my $min = ($a || $b || $c) - 1;
	$a = $a || $min;
	$b = $b || $min;
	$c = $c || $min;
	return $a if($a >= $b && $a >= $c);
	return $b if($b >= $a && $b >= $c);
	return $c if($c >= $a && $c >= $b);
}

sub _max {
	my $x = shift;
	my $tmp = _max3($x->[_INTERVAL][1], 
		$x->[_LEFT] ? $x->[_LEFT][_MAX] : undef, 
		$x->[_RIGHT] ? $x->[_RIGHT][_MAX] : undef);
	$x->[_MAX] = $tmp;
	return $tmp;
}

sub _update_max {
	my $x = shift;
	#update the max field for each parent node
	_max($x);
	my $k = $x->[_PARENT];
	while($k) {
		my $tmp = $k->[_MAX];
		_max($k);
		last if($tmp == $k->[_MAX]); # no need to update further
		$k = $k->[_PARENT];
	}
}

sub _left_rotate {
    my $self = shift;
    my $x = shift or croak('Missing arg: node');

    my $y = $x->[_RIGHT]
      or return;
    $x->[_RIGHT] = $y->[_LEFT];
    if($y->[_LEFT]) {
        $y->[_LEFT]->[_PARENT] = $x;
    }
    $y->[_PARENT] = $x->[_PARENT];
    if(not defined $x->[_PARENT]) {
        $self->[ROOT] = $y;
    }
    else {
        $x == $x->[_PARENT]->[_LEFT]
          ? $x->[_PARENT]->[_LEFT]  = $y
          : $x->[_PARENT]->[_RIGHT] = $y;
    }
    $y->[_LEFT]   = $x;
    $x->[_PARENT] = $y;
	#update max
	_max($x); _max($y);
}

sub _right_rotate {
    my $self = shift;
    my $y = shift or croak('Missing arg: node');

    my $x = $y->[_LEFT]
      or return;
    $y->[_LEFT] = $x->[_RIGHT];
    if($x->[_RIGHT]) {
        $x->[_RIGHT]->[_PARENT] = $y
    }
    $x->[_PARENT] = $y->[_PARENT];
    if(not defined $y->[_PARENT]) {
        $self->[ROOT] = $x;
    }
    else {
        $y == $y->[_PARENT]->[_RIGHT]
          ? $y->[_PARENT]->[_RIGHT] = $x
          : $y->[_PARENT]->[_LEFT]  = $x;
    }
    $x->[_RIGHT] = $y;
    $y->[_PARENT] = $x;
	_max($y); _max($x);
}

1; 

# Magic true value required at end of module
__END__

=head1 NAME

Tree::RB - Perl implementation of the Red/Black tree, a type of balanced binary search tree. 


=head1 VERSION

This document describes Tree::RB version 0.1


=head1 SYNOPSIS

    use Tree::RB;

    my $tree = Tree::RB->new;
    $tree->put('France'  => 'Paris');
    $tree->put('England' => 'London');
    $tree->put('Hungary' => 'Budapest');
    $tree->put('Ireland' => 'Dublin');
    $tree->put('Egypt'   => 'Cairo');
    $tree->put('Germany' => 'Berlin');

    $tree->put('Alaska' => 'Anchorage'); # D'oh!
    $tree->delete('Alaska');

    print $tree->get('Ireland'); # 'Dublin'

    print $tree->min->key; # 'Egypt' 
    print $tree->max->key; # 'Ireland' 
    print $tree->size; # 6

    # print items, ordered by key
    my $it = $tree->iter;

    while(my $node = $it->next) {
        sprintf "key = %s, value = %s\n", $node->key, $node->val;
    }

    # print items in reverse order
    $it = $tree->rev_iter;

    while(my $node = $it->next) {
        sprintf "key = %s, value = %s\n", $node->key, $node->val;
    }

    # Hash interface
    tie my %capital, 'Tree::RB';

    # or do this to store items in descending order 
    tie my %capital, 'Tree::RB', sub { $_[1] cmp $_[0] };

    $capital{'France'}  = 'Paris';
    $capital{'England'} = 'London';
    $capital{'Hungary'} = 'Budapest';
    $capital{'Ireland'} = 'Dublin';
    $capital{'Egypt'}   = 'Cairo';
    $capital{'Germany'} = 'Berlin';

    # print items in order
    while(my ($key, $val) = each %capital) {
        printf "key = $key, value = $val\n";
    }

=head1 DESCRIPTION

This is a Perl implementation of the Red/Black tree, a type of balanced binary search tree. 

A tied hash interface is also provided to allow ordered hashes to be used.

See the Wikipedia article at L<http://en.wikipedia.org/wiki/Red-black_tree> for more on Red/Black trees.


=head1 INTERFACE

=head2 new([CODEREF])

Creates and returns a new tree. If a reference to a subroutine is passed to
new(), the subroutine will be used to override the tree's default lexical
ordering and provide a user a defined ordering. 

This subroutine should be just like a comparator subroutine used with L<sort>, 
except that it doesn't do the $a, $b trick.

For example, to get a case insensitive ordering

    my $tree = Tree::RB->new(sub { lc $_[0] cmp lc $_[1]});
    $tree->put('Wall'  => 'Larry');
    $tree->put('Smith' => 'Agent');
    $tree->put('mouse' => 'micky');
    $tree->put('duck'  => 'donald');

    my $it = $tree->iter;

    while(my $node = $it->next) {
        sprintf "key = %s, value = %s\n", $node->key, $node->val;
    }

=head2 resort(CODEREF)

Changes the ordering of nodes within the tree. The new ordering is
specified by a comparator subroutine which must be passed to resort().

See L</new> for further information about the comparator.

=head2 size()

Returns the number of nodes in the tree.

=head2 root()

Returns the root node of the tree. This will either be undef
if no nodes have been added to the tree, or a L<Tree::RB::Node> object.
See the L<Tree::RB::Node> manual page for details on the Node object.

=head2 min()

Returns the node with the minimal key.

=head2 max()

Returns the node with the maximal key.

=head2 lookup(KEY, [MODE])

When called in scalar context, lookup(KEY) returns the value
associated with KEY.

When called in list context, lookup(KEY) returns a list whose first
element is the value associated with KEY, and whose second element
is the node containing the key/value.

An optional MODE parameter can be passed to lookup() to influence
which key is returned.

The values of MODE are constants that are exported on demand by
Tree::RB

    use Tree::RB qw[LUEQUAL LUGTEQ LULTEQ LUGREAT LULESS LUNEXT LUPREV];

=over

=item LUEQUAL

Returns node exactly matching the key.

=item LUGTEQ

Returns the node exactly matching the specified key, 
if this is not found then the next node that is greater than the specified key is returned.

=item LULTEQ

Returns the node exactly matching the specified key, 
if this is not found then the next node that is less than the specified key is returned.

=item LUGREAT

Returns the node that is just greater than the specified key - not equal to. 
This mode is similar to LUNEXT except that the specified key need not exist in the tree.

=item LULESS

Returns the node that is just less than the specified key - not equal to. 
This mode is similar to LUPREV except that the specified key need not exist in the tree.

=item LUNEXT

Looks for the key specified, if not found returns C<undef>. 
If the node is found returns the next node that is greater than 
the one found (or C<undef> if there is no next node). 

This can be used to step through the tree in order.

=item LUPREV

Looks for the key specified, if not found returns C<undef>. 
If the node is found returns the previous node that is less than 
the one found (or C<undef> if there is no previous node). 

This can be used to step through the tree in reverse order.

=back

=head2 get(KEY)

get() is an alias for lookup().

=head2 iter([KEY])

Returns an iterator object that can be used to traverse the tree in order.

The iterator object supports a 'next' method that returns the next node in the
tree or undef if all of the nodes have been visited.

See the synopsis for an example.

If a key is supplied, the iterator returned will traverse the tree in order starting from
the node with key greater than or equal to the specified key.

    $it = $tree->iter('France');
    my $node = $it->next;
    print $node->key; # -> 'France'

=head2 rev_iter([KEY])

Returns an iterator object that can be used to traverse the tree in reverse order.

If a key is supplied, the iterator returned will traverse the tree in order starting from
the node with key less than or equal to the specified key.

    $it = $tree->rev_iter('France');
    my $node = $it->next;
    print $node->key; # -> 'England'

=head2 hseek(KEY, [{-reverse => 1|0}])

For tied hashes, determines the next entry to be returned by each.

    tie my %capital, 'Tree::RB';

    $capital{'France'}  = 'Paris';
    $capital{'England'} = 'London';
    $capital{'Hungary'} = 'Budapest';
    $capital{'Ireland'} = 'Dublin';
    $capital{'Egypt'}   = 'Cairo';
    $capital{'Germany'} = 'Berlin';
    tied(%capital)->hseek('Germany');

    ($key, $val) = each %capital;
    print "$key, $val"; # -> Germany, Berlin 

The direction of iteration can be reversed by passing a hashref with key '-reverse' and value 1
to hseek after or instead of KEY, e.g. to iterate over the hash in reverse order:

    tied(%capital)->hseek({-reverse => 1});
    $key = each %capital;
    print $key; # -> Ireland 

The following calls are equivalent

    tied(%capital)->hseek('Germany', {-reverse => 1});
    tied(%capital)->hseek({-key => 'Germany', -reverse => 1});

=head2 put(KEY, VALUE)

Adds a new node to the tree. 

The first argument is the key of the node, the second is its value. 

If a node with that key already exists, its value is replaced with 
the given value and the old value is returned. Otherwise, undef is returned.

=head2 delete(KEY)

If the tree has a node with the specified key, that node is
deleted from the tree and returned, otherwise C<undef> is returned.


=head1 DEPENDENCIES

L<enum>


=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

Please report any bugs or feature requests to
C<bug-tree-rb@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.


=head1 AUTHOR

Arun Prasad  C<< <arunbear@cpan.org> >>

Some documentation has been borrowed from Benjamin Holzman's L<Tree::RedBlack>
and Damian Ivereigh's libredblack (L<http://libredblack.sourceforge.net/>).

=head1 ACKNOWLEDGEMENTS

Thanks for bug reports go to Anton Petrusevich, Wes Thompson and Christopher Gurnee.

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2007, Arun Prasad C<< <arunbear@cpan.org> >>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.
