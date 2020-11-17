package DelimitedFile;

use strict;
use FileHandle;
use Carp qw(confess);

use Reporter;
#use misc;
# rude export behavior, I think universal_open() is the only function
use FileUtils qw(universal_open);

@DelimitedFile::ISA = qw(Configurable Exporter);
@DelimitedFile::EXPORT_OK = qw(df_bucket_by_header);

use Configurable;
use Exporter;


use constant HEADER_DETECT_READ_LINE_COUNT => 30;

#use constant DELIMITERS => ("\t", ",", '\s+');
#use constant DELIMITERS => ("\t", ",", '\s');
use constant DELIMITERS => ("\t", ",");
# 2/2005:
#  - change from "\s+" to "\s" (in case of null column)
#  - oops: disallow \s+ altogether...if tab-delimited file contains text fields,
#          the spaces could be interpreted as "delimiters".
#          Is delim of \s+ *ever* correct??
#         
#    for tab-d

use MethodMaker qw(
		   delimiter
                   quote_character

		   fh
		   file
		   headers
		   headers_raw
		   headers_orig
		   autotrim
		   last_line
		   field_count
		   skip_until
                   unquote
                   skipped_lines
                   headers_uppercase
skip_comments
skip_double_comments
save_last_line

header_line
skip_duplicate_header_lines
uncomment_header
double_comments
		  );

sub new {
  # -file ...
  # -headers ...
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->save_last_line(1);
  $self->skip_duplicate_header_lines(1);
  $self->double_comments([]);
  $self->configure(%options);
#  open my $fh, $options{"-file"};
# ??? http://www.perl.com/pub/1999/11/sins.html ???
  $self->field_count(-1);
  $self->open_file() unless $self->fh();
  $self->header_setup() if $self->headers();
  return $self;
}

sub open_file {
  # (re-)open file
  my ($self) = @_;
  my $file = $self->file || confess "need -file";
  #  printf STDERR "file: %s\n", $options{"-file"};
  my $fh;
  if ($file =~ /\.gz$/i) {
    $fh = universal_open($file);
  } else {
    $fh = new FileHandle();
    $fh->open($file) || confess "can't open file $file: $?";
  }
  $self->fh($fh);
}

sub header_setup {
  my ($self, %options) = @_;

  my $pattern = $self->skip_until();
  my @headers;
  if ($self->headers == 3) {
    # autodetect header line (and delimiter)
    my ($i, $line);
    my $fh = $self->fh;
    my %delims;
    my %dcount;
    my %first_line;
    for ($i=0; $i < HEADER_DETECT_READ_LINE_COUNT; $i++) {
      $line = <$fh>;
      chomp $line;
      $line =~ s/\x0d$//;
      
      last unless $line;
      my $delim = get_best_delimiter($line);
      if ($delim) {
	# valid delimiter detected
	$delims{$delim}++;
	my @fields = split($delim, $line);
	my $fcount = scalar @fields;
	$dcount{$delim}{$fcount}++;
	unless (defined $first_line{$delim}{$fcount}) {
	  $first_line{$delim}{$fcount} = [$i, \@fields];
	}
      }
    }
    
    my $best_delim = (sort {$delims{$b} <=> $delims{$a}} keys %delims)[0];
    my $best_count = (sort {$dcount{$best_delim}{$b} <=> $dcount{$best_delim}{$a}}
		      keys %{$dcount{$best_delim}})[0];
    my ($header_line_number, $header_ref) = @{$first_line{$best_delim}{$best_count}};
    $self->delimiter($best_delim);

    @headers = @{$header_ref};
    $fh = $self->open_file;
    # re-open file and move pointer past header
    my @skipped_lines;
    for ($i=0; $i <= $header_line_number; $i++) {
      $line = <$fh>;
      $line =~ s/\x0d$//;
      push @skipped_lines, $line;
    }
    $self->skipped_lines(\@skipped_lines);
  } elsif ($self->headers() == 4) {
    # passed-in headers for headerless file
    my $headers_raw = $self->headers_raw || die "mode 4 requires headers_raw";
    @headers = @{$headers_raw};
  } else {
    # autodetect
    while (1) {
      @headers = $self->next();
      last if $pattern ? $self->last_line =~ /$pattern/ : 1;
      # skip leading junk lines until we see pattern
    }
  }
  my %map;

  if ($self->headers == 2) {
    # strip quotes from header names
    foreach (@headers) {
      s/^\"//;
      s/\"$//;
    }
  }
  
  if ($self->headers_uppercase) {
    foreach (@headers) {
      $_ = uc($_);
    }
  }
  if ($self->uncomment_header) {
    $headers[0] =~ s/^\#//;
  }


  $self->headers_raw([ @headers ]);
  $self->field_count(scalar @headers);
  my @headers_cooked;
  for (my $i=0; $i < @headers; $i++) {
#    printf STDERR "%s\n", $headers[$i];
    if (exists $map{$headers[$i]}) {
      # duplicate header
      for (my $suffix = 2; ; $suffix++) {
	my $label = $headers[$i] . $suffix;
	unless (exists $map{$label}) {
	  $map{$label} = $i;
	  printf STDERR "Duplicate header \"%s\", using \"%s\" instead.\n", $headers[$i], $label;
	  push @headers_cooked, $label;
	  last;
	}
      }
    } else {
      $map{$headers[$i]} = $i;
      push @headers_cooked, $headers[$i];
    }
  }
  $self->header_line($self->last_line);
  # needs work for other than basic parsing
  $self->headers_orig($self->headers_raw);
  $self->headers_raw(\@headers_cooked);
  $self->headers(\%map);
}

sub set_manual_headers {
  my ($self, $headers) = @_;
  my %map;
  for (my $i = 0; $i < @{$headers}; $i++) {
    $map{$headers->[$i]} = $i;
  }
  $self->headers(\%map);
  $self->headers_raw($headers);
}

sub get {
  # get the named columns from the next entry
  # requires -headers
  my ($self, @columns) = @_;
  my $map = $self->headers();
  my @results;
  if (my @fields = $self->next()) {
    foreach (@columns) {
      if (exists $map->{$_}) {
	push @results, $fields[$map->{$_}];
      } else {
	die sprintf "unknown column $_ (%s)", join ",", keys %{$map};
      }
    }
    return @results;
  } else {
    return ();
  }
}

sub get_hash {
  # return next row as a hash
  # requires -headers
  my ($self, %options) = @_;
  my $unquote = $options{"-unquote"} || $self->unquote();

  my $fields;
  
  while (1) {
    $fields = $self->next("-ref" => 1);
#  die scalar @{$fields};
    if ($self->skip_duplicate_header_lines) {
      last unless $fields;
      my $hl = $self->header_line;
      if ($hl and $self->last_line eq $hl) {
	printf STDERR "DelimitedFile: skipping duplicate header line\n";
      } else {
	last;
      }
    } else {
      last;
    }
  }

  if ($fields and @{$fields}) {
    my %hash;
#    while (my ($label, $index) = each %{$self->headers}) {
#      $hash{$label} = $fields[$index];
#    }
    if ($unquote) {
      foreach (@{$fields}) {
	s/^\"//;
	s/\"$//;
      }
    }

    confess "no headers" unless $self->headers_raw;
    @hash{@{$self->headers_raw}} = @{$fields};
    return \%hash;
  } else {
    return undef;
  }
}

sub get_all {
  # get the named columns from the next entry
  # requires -headers
  my ($self, @columns) = @_;
  my $map = $self->headers();

  foreach (@columns) {
    die sprintf "unknown column $_ (%s)", join ",", keys %{$map}
      unless exists $map->{$_};
  }
  my @map = map { $map->{$_} } @columns;
  my (@fields, @results);
  my $fields;

  if (@columns == 1) {
    while (@fields = $self->next()) {
      push @results, $fields[$map[0]];
    }
  } elsif (0) {
    while (@fields = $self->next()) {
      push @results, [ map { $fields[$_] } @map ];
    }
  } else {
    while ($fields = $self->next("-ref" => 1)) {
      push @results, [ @{$fields}[@map] ];
    }
  }
  return \@results;
}

sub get_all_optimized {
  #
  # fast 'n' loose
  #
  my ($self, %options) = @_;
  my $map = $self->headers();
  my $fh = $self->fh();
  my $delimiter = $self->delimiter();
  my $hash = $options{"-hash"};
  my $trim_whitespace = $options{"-trim-whitespace"};

  my $columns;
  if ($options{"-all"}) {
    # return all columns
    $columns = $self->headers_raw();
  } else {
    $columns = $options{"-columns"} || die "specify -columns";
  }

  my @map;
  if ($options{"-index"}) {
    # column indices
    @map = @{$columns};
  } else {
    # column names
    foreach (@{$columns}) {
      die sprintf "unknown column $_ (%s)", join ",", keys %{$map}
	unless exists $map->{$_};
    }
    @map = map { $map->{$_} } @{$columns};
  }

  my @results;
  my $line;
  my $ref;

  if ($options{"-hash"}) {
    while ($line = <$fh>) {
      chomp $line;
      $delimiter = get_best_delimiter($line) unless $delimiter;
      $line =~ s/\r$//;
      my %r;
      @r{@{$columns}} = (split(/$delimiter/, $line))[@map];
      if ($trim_whitespace) {
	foreach (values %r) {
	  if (defined $_) {
	    s/^\s+//;
	    s/\s+$//;
	  }
	}
      }
      push @results, \%r;
    }
  } else {
    while ($line = <$fh>) {
      chomp $line;
      $delimiter = get_best_delimiter($line) unless $delimiter;
      $line =~ s/\r$//;
      $ref = [(split(/$delimiter/, $line))[@map]];
      if ($trim_whitespace) {
	foreach (@{$ref}) {
	  s/^\s+//;
	  s/\s+$//;
	}
      }
      push @results, $ref;
    }
  }

  return \@results;
}


sub next {
  my ($self, %options) = @_;
  my $fh = $self->fh;
  my $line;
  if ($options{"-line"}) {
    # parse given line
    $line = $options{"-line"};
  } else {
    $line = <$fh>;
  }
  if ($line) {
    my $delimiter = $self->delimiter();
    unless ($delimiter) {
      if ($self->skip_comments()) {
	while ($line =~ /^#/) {
	  $line = <$fh>;
	}
      } elsif ($self->skip_double_comments()) {
	while ($line =~ /^##/) {
	  push @{$self->double_comments()}, $line;
	  $line = <$fh>;
	}
      }
      $self->delimiter($delimiter = get_best_delimiter($line));
    }

    my $unquote = $self->unquote();

    if ($unquote) {
      my $quote_char = $self->quote_character();
      $self->quote_character($quote_char = get_best_quote($line)) unless $quote_char;
      $delimiter = $quote_char . $delimiter . $quote_char;
      # HACK: if values are quoted, add the quotation character to
      # the delimiter.  This prevents delimiter characters within
      # the quoted string from breaking the split() below.
      #
      # e.g. in Affymetrix annotation file HC_G110_annot.csv,
      # column "Annotation Date" contains value "Mar 31, 2005"
      # so splitting on "," wrongly splits the date into 2 fields!
      #
      # "better" approach would mask quoted regions rather than 
      # using split(), but hey.
    }

    my @fields;

    if ($self->autotrim) {
      while (1) {
	last if index($line, $delimiter) >= 0;
	$line = <$fh>;
	return () unless defined $line;
      }
    }

    chomp $line;
    $line =~ s/\x0d$//;
    # remove trailing carriage returns (msdos)

    $self->last_line($line);
#    return split /$delimiter/, $line;
#    return split /$delimiter/, $line, 11;
    my $fc = $self->field_count();
    @fields = split /$delimiter/, $line, $fc;
    # don't strip trailing empty fields
    
    my $have = scalar @fields;
    if ($have < $fc) {
      # there are fewer cells in this row than in the "headers" row.
      # pad with empty (but not undef) values.
      my $i;
      for ($i=$have + 1; $i <= $fc; $i++) {
	$fields[$i - 1] = "";
      }
    } elsif ($fc > 0 and $have > $fc) {
      printf STDERR "WTF: have %d cells, expected %d\n", $have, $fc;
    }

    if ($unquote) {
      foreach (@fields) {
	s/^\"//;
	s/\"$//;
      }
    }
    
    return $options{"-ref"} ? \@fields : @fields;
  } else {
    $self->last_line("");
    return ();
  }
}

sub get_best_quote {
  # STATIC
  # quoted values: get most likely quotation character
  my ($line) = @_;
  my $best_char;
  my $best_count = 0;
  foreach my $char ("'", '"') {
    my @hits = $line =~ /$char/g;
    if (@hits and @hits > $best_count) {
      $best_char = $char;
      $best_count = scalar @hits;
    }
  }
  return $best_char;
}

sub get_best_delimiter {
  # STATIC
  my ($line) = @_;
  chomp $line;

  my ($best_d, $c, @fields);
  my $best_c = 0;
  my $all_valid;
  my $d;
  foreach $d (DELIMITERS) {
    @fields = split /$d/, $line;
    $c = scalar @fields;
    $all_valid = 1;

    if (0) {
      # 2/05: don't disqualify a delimiter if a column is null
      foreach (@fields) {
	#      $all_valid = 0, last unless /\w/;
	$all_valid = 0, last unless /\S/;
      }
    }
#    printf STDERR "%s: %d av=%d\n", $d, $c, $all_valid;

    if ($all_valid and $c > 1 and $c > $best_c) {
      $best_c = $c;
      $best_d = $d;
    }
  }
  
  confess "no delimiter detected!" unless $best_d;
  
  return $best_d;
}

sub column_sort {
  # sort file by some columns
  my ($self, %options) = @_;

  my $columns = $options{"-columns"} || die "-columns";
  my $headers = $self->headers() || die "column titles required";
  
  my @rows;
  my $row;
  while ($row = $self->next("-ref" => 1)) {
    push @rows, $row;
  }

  my @i;
  my $index;
  foreach (@{$columns}) {
    $index = $headers->{$_};
    confess "can't find column $_!" unless defined $index;
    push @i, $index;
  }
  
  my $results = multi_column_sort("-rows" => \@rows,
				  "-columns" => \@i);

  if (my $out_fn = $options{"-write"}) {
    # write results to file

    my $delimiter = $self->delimiter();

    open(DF_OUT, ">$out_fn") || die "can't write to $out_fn";

    if (my $h = $self->headers_raw) {
      printf DF_OUT "%s\n", join $delimiter, @{$h};
    }

    foreach (@{$results}) {
      printf DF_OUT "%s\n", join $delimiter, @{$_};
    }
    close DF_OUT;
  }

  return $results;
}

sub get_reporter {
  # return new reporter, cloning current report format
  my ($self, %options) = @_;
  my @labels = @{$self->headers_raw};
  if (my $map = $options{"-rename"}) {
    # rename column headers.
    # note it's the responsibility of the application to actually
    # copy the values in each row.
    @labels = map {$map->{$_} || $_} @labels;
    delete $options{"-rename"};
  }

  if (my $fields = $options{"-extra"}) {
    my @f;
    if ($options{"-clobber"}) {
      my %old = map {$_, 1} @labels;
      foreach my $f (@{$fields}) {
	push @f, $f unless $old{$f};
      }
      delete $options{"-clobber"};
    } else {
      @f = @{$fields};
    }

    if (@f) {
      if ($options{"-extra-prepend"}) {
	unshift @labels, @f;
	delete $options{"-extra-prepend"};
      } else {
	push @labels, @f;
      }
    }
    delete $options{"-extra"};
  }
  if (my $fields = $options{"-remove"}) {
    foreach my $f (@{$fields}) {
      @labels = grep {$_ ne $f} @labels;
    }
    delete $options{"-remove"};
  }

  my $delimiter = $options{"-delimiter"} || $self->delimiter;

  return new Reporter(
		      %options,
		      "-delimiter" => $options{"-html"} ? undef : $delimiter,
		      "-labels" => \@labels,
		     );
}

sub df_bucket_by_header {
  # STATIC utility routine
  my (%options) = @_;
  my $file = $options{"-file"} || die "-file";
  my $column = $options{"-column"} || die "-column";

  my $df = new DelimitedFile("-file" => $file,
			     "-headers" => 1,
			     );
  my %bucket;
  while (my $row = $df->get_hash()) {
    die unless exists $row->{$column};
    my $value = $row->{$column};
    push @{$bucket{$value}}, $row;
    # TO DO: modes to specify unique bucketing, etc.
  }
  return \%bucket;
}

sub check_headers_unique {
  my ($self, %options) = @_;
  my %count;
  foreach my $h (@{$self->headers_orig}) {
    $count{$h}++;
  }

  my @broken;
  foreach (sort keys %count) {
    push @broken, $_ if $count{$_} > 1;
  }

  if (@broken and $options{"-die"}) {
    die sprintf "ERROR: some column labels are not unique (%s).  Unique labels are required to guarantee information is read from the correct column.\n", join ", ", @broken;
  }

  return @broken ? 0 : 1;
}


1;
