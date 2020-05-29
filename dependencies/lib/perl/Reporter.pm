# abstract report generation
# TO DO:
# - auto-QC mode to verify output has same # of columns in each row?

package Reporter;

use strict;
use Configurable;
use WorkingFile;
use Carp qw(confess);
use MiscUtils qw(dump_die);
use FileUtils qw(universal_open);

use MethodMaker qw (
		    labels
                    print_labels
		    delimiter
		    headers_done
		    html
		    wrote
		    title
		    refs
		    wf
		    fh
		    file
                    compress
                    html_headers_done

                    html_start_hook
		    start_hook
		    end_hook
		    tr
		    html_strip
                    allow_undef

auto_qc
comment_header
rename_headers
		   );

@Reporter::ISA = qw(Configurable);

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);

  my $file = $options{"-file"} || "";
  if ($file eq "-") {
    $self->fh(*main::STDOUT);
  } elsif ($file) {
    my $wf = new WorkingFile($file,
			     "-compress" => $options{"-compress"}
			    );
    $self->wf($wf);
    $self->fh($wf->output_filehandle);
  } elsif (my $fh = $options{"-fh"}) {
    $self->fh($fh);
  }
  return $self;
}

sub add {
  # add columns/cells to current row
  my $self = shift;
  push @{$self->{"queue"}}, @_;
}

sub end_row {
  # finish the row/line
  my ($self, @stuff) = @_;
  my $queue = $self->{"queue"};
  my $labels = $self->labels();
  my $old_fh;
  my $fh = $self->fh;
  if ($fh) {
    $old_fh = select();
    select $fh;
  }
  my $allow_undef = $self->allow_undef;

  if (@stuff) {
    if (ref $stuff[0] eq "HASH") {
      # values in hash
      my $i = 0;
      my $label;
      my $value;
      my $html_strip = $self->html_strip;

      if (0) {
	printf STDERR "labels: %s\n", join ",", @{$labels};
	dump_die($stuff[0], "debug row", 1);
      }
#      die join ",", @{$labels};
      for ($i=0; $i <= $#$labels; $i++) {
	$label = $labels->[$i];

	unless (exists $stuff[0]->{$label}) {
	  if ($allow_undef) {
	    $stuff[0]->{$label} = "";
	  } else {
	    dump_die($stuff[0], "no value for $label");
	    # die "undef value for $label" unless defined $stuff[0]->{$label};
	  }
	}

	$value = $stuff[0]->{$label};
	$value =~ s/<[^>]+>//g if $html_strip;
	push @{$queue}, $value;
      }
    } else {
      push @{$queue}, @stuff;
    }
  }

  my $html = $self->html();
  my $delimiter = $self->delimiter();
  my $refs = $self->refs;
  my $tr_tag = $self->tr || '<tr>';

  $self->write_headers() unless $self->headers_done();

  if ($delimiter) {
    if ($refs) {
      printf "%s\n", join $delimiter, map {$$_} @{$refs};
    } else {
      printf "%s\n", join $delimiter, @{$queue};
    }
  } elsif ($html) {
    printf "%s%s</tr>\n", $tr_tag, join " ", map {"<td>" . (defined $_ ? $_ : "")} @{$queue};
    if ($self->wrote(($self->wrote || 0) + 1) % 100 == 0) {
      print "\n</table>\n";
      $self->headers_done(0);
    }
  } else {
    die "no headers!" unless $labels;
    my @mapped;
    for (my $i=0; $i < @{$queue}; $i++) {
      push @mapped, sprintf '%s%s', 
	($i <= $#$labels ? $labels->[$i] . ":" : ""), $queue->[$i];
    }
    printf "%s\n", join " ", @mapped;
  }

  @{$queue} = ();

  select $old_fh if $old_fh;
}

sub run_hook {
  # STATIC
  my ($self, $callback) = @_;
  &$callback($self, $self->fh()) if $callback;
}

sub finish {
  my $fh = $_[0]->fh;
  my $old_fh;
  if ($fh) {
    $old_fh = select();
    select $fh;
  }

  if ($_[0]->html) {
    print "</table>\n";
  }

  $_[0]->run_hook($_[0]->end_hook);

#  printf STDERR "finish rpt %s html=%s\n", $_[0], $_[0]->html;

  if (($_[0]->html || 0) == 2) {
    print "</body></html>\n" if $_[0]->html() == 2;
  }

  if ($_[0]->wf) {
    $_[0]->wf->finish();
  } else {
    close $fh if $fh;
  }
  select $old_fh if $fh;
  
  if ($_[0]->auto_qc()) {
    $_[0]->run_auto_qc();
  } elsif (my $file = $_[0]->file) {
    printf STDERR "NOTE: auto_qc not enabled for %s\n", $file;
  }
}

sub bind {
  # bind column names to variable references.
  # sets headers for report.
  my ($self, @list) = @_;
  my (@headers, @refs);
  while (@list) {
    my ($header, $ref) = splice(@list, 0, 2);
    die "bogus ref for $header" unless ref $ref and ref $ref eq "SCALAR";
    push @headers, $header;
    push @refs, $ref;
  }
  $self->labels(\@headers);
  $self->refs(\@refs);
}

sub end_row_hash_optimized {
  # finish the row/line
  my ($self, $ref) = @_;
  my $labels = $self->labels();
  my @out = @{$ref}{@{$labels}};
  my $delimiter = $self->delimiter();

  unless ($self->headers_done) {
    my $use = $self->print_labels || $labels;
    printf "%s\n", join $delimiter, @{$use} if @{$use};
    $self->headers_done(1);
  }

  printf "%s\n", join $delimiter, @out;
}

sub write_headers {
  my ($self) = @_;
  my $html = $self->html();
  my $delimiter = $self->delimiter();
  my $labels = $self->labels();

  unless ($self->headers_done()) {
    my $old_fh;
    if (my $fh = $self->fh) {
      $old_fh = select();
      select $fh;
    }

    $self->run_hook($self->start_hook);
    if ($html) {
      if ($html == 2 and !$self->html_headers_done) {
	my $title = $self->title || "";
	printf '<html><head><title>%s</title></head><body>', $title || "Report";
	printf '<h1>%s</h1>', $title if $title;
	$self->run_hook($self->html_start_hook);
	$self->html_headers_done(1);
	# only write headers once
      }
      print '<table border=1>';
    }
    $self->run_hook($self->start_hook);
    if ($labels) {
      my $use = $self->print_labels || $labels;
      my @h;
      if (my $rename = $self->rename_headers()) {
	@h = map {exists $rename->{$_} ? $rename->{$_} : $_} @{$use};
      } else {
	@h = @{$use};
      }

      if ($delimiter) {
	printf "%s%s\n",
	($self->comment_header ? "#" : ""),
	join $delimiter, @h;
      } elsif ($html) {
	printf "<tr>%s</tr>\n", join " ", map {"<th>" . $_} @h;
      }
    }
    $self->headers_done(1);

    select $old_fh if $old_fh;
  }
}

sub run_auto_qc {
  my ($self) = @_;
  my $file = $self->file || return;
  # no file in STDOUT mode

  $file .= ".gz" if $self->compress() and not($file) =~ /\.gz$/;

  my $delimiter = $self->delimiter || die;

  my $fh = universal_open($file) || die "can't open $file";
  my %counts;
  my %example;
  while (<$fh>) {
    chomp;
    my @f = split /$delimiter/, $_, -1;
    # -1: don't strip trailing empty fields
    my $count = scalar @f;
    $counts{$count}++;
    $example{$count} = $_;
  }
  close $fh;

  if (scalar keys %counts > 1) {
    # FAIL
    printf STDERR "ERROR in %s: multiple column counts!:\n", $file;
    foreach my $count (sort {$a <=> $b} keys %counts) {
      printf STDERR "  cols=%d row_count=%d example=%s\n",
      $count, $counts{$count}, $example{$count};
    }
    die;
  }
}


1;
