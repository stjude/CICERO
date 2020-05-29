package CommandLineRebuilder;

# my $clr = new CommandLineRebuilder("-parameters" => \@options,
#				   "-flags" => \%FLAGS);

use strict;
use Carp qw(confess);

use Configurable;

@CommandLineRebuilder::ISA = qw(Configurable Exporter);
@CommandLineRebuilder::EXPORT_OK = qw();

use MethodMaker qw(
	parameters
        flags
exclude_list
include_list
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->exclude_list({});
  $self->include_list({});
  $self->configure(%options);
  return $self;
}

sub get_command_line {
  my ($self, %options) = @_;
  my $flags = $self->flags();
  my $exclude = $self->exclude_list();
  my $include = $self->include_list();
  my @params;
  my $parameters = $self->parameters();

  for (my $i = 0; $i < @{$parameters}; $i++) {
    # for each possible parameter
    my $param_name = $parameters->[$i];
#    printf STDERR "debug param %s\n", $param_name;
    die "param name is a ref, WTF?" if ref $param_name;

    if ($param_name =~ /\=/) {
      # parameter accepting a value
      my $name = $param_name;
      die unless $name =~ /\=\w$/;
      $name =~ s/\=\w//;
      my $key = $name;
      $key =~ s/^\-// || die;
      if ($exclude->{$name}) {
	# don't pass through this parameter from environment
	my $v = $options{$name};
	push @params, $name, $v if defined $v;
	# custom value specified here (e.g. a single file to process)
	if (($i < scalar(@{$parameters}) - 1) and ref $parameters->[$i + 1]) {
	  # argument to param is an arrayref, skip
	  $i++;
	}
      } else {
	my $v;
	if (exists $flags->{$key}) {
	  # single value placed into ARGV
	  $v = $flags->{$key};
	} elsif (exists $options{$name}) {
	  # user-specified value
	  $v = $options{$name};
	  if (($i < scalar(@{$parameters}) - 1) and ref $parameters->[$i + 1]) {
	    # argument to param is an arrayref, skip
	    $i++;
	  }

	} elsif ($i < scalar(@{$parameters}) - 1) {
	  if (ref $parameters->[$i + 1]) {
	    my $set = $parameters->[++$i];
	    if (ref $set eq "ARRAY") {
	      foreach (@{$set}) {
		push @params, $name, $_;
	      }
	    } elsif (ref $set eq "SCALAR") {
#	      push @params, $name, $$set;
	      $v = $$set;
	      # pass defined check below
	    } else {
	      die "unhandled ref type " . ref $set;
	    }
	  }
	}
	push @params, $name, quote_check($v) if defined $v;
      }
    } else {
      # flag-style parameter
      my $key = $param_name;
      $key =~ s/^\-//;
      unless ($exclude->{$param_name}) {
	push @params, $param_name if $flags->{$key} or
	    $options{$param_name} or
	    $include->{$param_name};
      }
      if (($i < scalar(@{$parameters}) - 1) and ref $parameters->[$i + 1]) {
	# argument to param is an arrayref, skip
	$i++;
      }

    }
  }
#  push @params, "-devel-path" if DevelopmentPath::is_development();
  my $cmd = join " ", $0, @params;
  return $cmd;
}

sub exclude_parameter {
  my ($self, $exclude) = @_;
  die "param requires leading -" unless $exclude =~ /^\-/;
  $self->exclude_list->{$exclude} = 1;
}

sub include_parameter {
  my ($self, $include) = @_;
  die "param requires leading -" unless $include =~ /^\-/;
  $self->include_list->{$include} = 1;
}

sub get_values_for {
  # return values for a single flag or parameter.
  # Since I present have a head cold, duplicate-y code from above.
  my ($self, %options) = @_;
  my $wanted = $options{"-param"} || die;
  $wanted =~ s/^\-//;

  my $parameters = $self->parameters();
  my $flags = $self->flags();

  my $result;
  my $found;

  for (my $i = 0; $i < @{$parameters}; $i++) {
    my $param_name = $parameters->[$i];
    my $key;
    my $v;
    if ($param_name =~ /\=/) {
      # parameter accepting a value
      my $name = $param_name;
      $name =~ s/\=\w$// || die;
      $key = $name;
      $key =~ s/^\-// || die;
      if (exists $flags->{$key}) {
	# single value placed into ARGV
	$v = $flags->{$key};
      } elsif (exists $options{$name}) {
	# user-specified value
	$v = $options{$name};
      } elsif ($i < scalar(@{$parameters}) - 1) {
	# lookahead to reference variable
	if (ref $parameters->[$i + 1]) {
	  my $set = $parameters->[++$i];
	  if (ref $set eq "ARRAY") {
	    $v = [ @{$set} ];
	  } elsif (ref $set eq "SCALAR") {
	    $v = $$set;
	  } else {
	    die "unhandled ref type " . ref $set;
	  }
	}
      }
    } else {
      # flag-style parameter
      $key = $param_name;
      $key =~ s/^\-//;
      $v = $flags->{$key} if exists $flags->{$key};
    }

    if ($key eq $wanted) {
      $found = 1;

      if (defined $v) {
	if (my $type = ref $v) {
	  die unless $type eq "ARRAY";
	  if (@{$v}) {
	    # only return if array actually contains data
	    $result = $v;
	  }
	} else {
	  # scalar
	  $result = [ $v ];
	}
      }
      last;
    }
  }
  
  confess("parameter $wanted not found") unless $found;

  return $result;

}

sub quote_check {
  my ($v) = @_;
  $v = "'" . $v . "'" if $v =~ /\W/;
  # quote values if necessary to protect from shell, e.g.
  # e.g. 
  #  hg19_pos(1-based)
  # becomes
  #  'hg19_pos(1-based)'
  #
  # ....otherwise bash will complain
  return $v;
}


1;
