package FileUtils;
# mne 3/02

use strict;

use Digest::MD5 qw(md5_hex);
use FileHandle;
use File::Find;
use Cwd 'abs_path';
use Carp qw(cluck confess);
use File::Path;
use File::Basename;
use POSIX qw(uname);
use Compress::Zlib;

use Configurable;
#use misc;
# replaced find_program() with find_binary(), I think that's it

use constant MAX_FILE_LENGTH => 255;
# HACK: maximum length of a filename on OS
# ...is this available programatically? POSIX::FILENAME_MAX returns 1024,
# maybe that's the max path length??

@FileUtils::ISA = qw(Configurable Exporter);
@FileUtils::EXPORT_OK = qw(
                           delete_directory
                           decompress_file
                           digest_cachefile
                           write_simple_file
                           read_simple_file
                           sort_files_by_date
                           get_local_arch_path
                           get_local_arch_binary
get_df
mkpath_cd

                           decompress_gzip
                           mtime
test_magic
find_binary
count_file_lines
universal_open
aggressive_rename
newer_than
line_delimiter_qc
open_filehandle

md5_file
is_offline
today_file
tab_cat
                          );

sub delete_directory {
  # recursively delete a directory (all contents and subdirectories).
  my ($target_dir) = @_;

  print STDERR "delete_directory is obsolete!  USE File::Path::rmtree!\n";

  if (-d $target_dir) {
    my %dir_depths;
    my $depth;
    my $dir_separator = '/';
    # hack, not portable!
    my @sd;

    my $first_pass = sub {
      if (-d $_) {
	if ($_ ne '.' and $_ ne '..') {
	  @sd = split /$dir_separator/, $File::Find::name;
	  $depth = scalar @sd;
	  push @{$dir_depths{$depth}}, $File::Find::name;
	}
      } else {
	# file, symlink, etc
	unlink $File::Find::name;
      }
    };

  
    find($first_pass, $target_dir);
    # delete files and get subdirectory tree 

    foreach $depth (sort {$b <=> $a} keys %dir_depths) {
      #    printf "%d:\n%s\n", $depth, join "\n", @{$dir_depths{$depth}};
      foreach (@{$dir_depths{$depth}}) {
	rmdir($_) || print STDERR "ack: can't rmdir $_\n";
      }
    }

    if (0) {
      print "debug, finding beneath $target_dir...\n";
      system "find $target_dir -ls";
    }

    if (cwd() eq $target_dir) {
      cluck "can't delete current directory!\n";
    }

    rmdir($target_dir);

    return -d $target_dir ? undef : 1;
    # return status
  } else {
    return 1;
  }
}

sub decompress_file {
  my ($file) = @_;
  my $binary;
  if ($file =~ /\.gz$/i) {
    $binary = "gzip";
  } elsif ($file =~ /\.bz2$/i) {
    $binary = "bzip2";
  } else {
    die;
  }
#  my $cmd = sprintf '%s -d %s', (find_program($binary) || die), $file;
  my $cmd = sprintf '%s -d %s', (find_binary($binary) || die), $file;
  system $cmd;
  return $? ? undef : 1;
}

sub digest_cachefile {
  # would have used md5_base64 but that can include "/" (pathname interference)
  my ($thing, %options) = @_;
  my $basedir = $options{"-base"} || ".md5_cache";

  my ($md5, $basename);
  if (ref $thing) {
    # array
    $md5 = md5_hex(@{$thing});
    $basename = join "", @{$thing};
  } else {
    $md5 = md5_hex($thing);
    $basename = $thing;
  }

  $basename = substr($basename, 0, MAX_FILE_LENGTH) if length($basename) > MAX_FILE_LENGTH;

  my $fn = join "/",
    $basedir,
      substr($md5, -2),
	substr($md5, -4),
	  $md5,
	    $basename;

  my $dir = dirname($fn);
  my $old_umask = umask();
  umask(0);
  unless (-d $dir) {
    mkpath($dir,0,0777) || die "can't mkpath $dir";
    # create world-readable/writable, ignoring umask
  }
  umask($old_umask);
  
  return $fn;
}

sub write_simple_file {
  # to do: hash
  my ($thing, $filename, %options) = @_;
  open(WRTMP, ">" . $filename) || die "can't write to $filename";
  if (ref $thing eq "ARRAY") {
    foreach (@{$thing}) {
      print WRTMP "$_\n";
    }
  } elsif (ref $thing eq "SCALAR") {
    print WRTMP $$thing;
  } else {
    die "unknown ref type";
  }
  close WRTMP;
}

sub read_simple_file {
  # TO DO: trim trailing \r (DOS)
  my ($filename, %options) = @_;
  my $tokenize = $options{"-tokenize"};
  confess "no filename" unless $filename;
  open(RDTMP, $filename) || confess "can't open $filename: $!";
  my $result;
  if ($options{"-raw"}) {
    local $/ = undef;
    my $blob = <RDTMP>;
    close RDTMP;
    $result = \$blob;
  } elsif (my $col = $options{"-hash"}) {
    # return results bucketed by a column
    my $single = $options{"-single"};
    require DelimitedFile;
#    printf STDERR "loading $filename...\n";
    my $df = new DelimitedFile("-file" => $filename,
			       "-headers" => 1);
    die "unknown field $col" unless defined $df->headers()->{$col};
    my %rows;
    if ($options{"-array"}) {
      while (my $row = $df->get_hash()) {
	my $key = $row->{$col} || die "no data for $col";
	push @{$rows{$key}}, $row;
      }
    } else {
      while (my $row = $df->get_hash()) {
	my $key = $row->{$col} || die "no data for $col";
	if ($single) {
	  if (exists $row->{$single}) {
	    $rows{$key} = $row->{$single};
	  } else {
	    die "no entry for $single";
	  }
	} else {
	  $rows{$key} = $row;
	} 
      }
    }
    $result = \%rows;
  } elsif ($options{"-as-hash"}) {
    # return arrayref of hashed rows
    require DelimitedFile;
    my $df = new DelimitedFile("-file" => $filename,
			       "-headers" => 1);
    my @rows;
    while (my $row = $df->get_hash()) {
      push @rows, $row;
    }
    if ($options{"-return-headers"}) {
      return (\@rows, $df->headers_raw());
    } else {
      $result = \@rows;
    }
  } elsif ($options{"-hash1"}) {
    # single-column file, return as hash
    $result = {};
    while (<RDTMP>) {
      chomp;
      s/\r$//;
      # DOS line ending
      if ($tokenize) {
	s/^\s+//;
	s/\s+$//;
	my @tokens = split /\s+/, $_;
	foreach (@tokens) {
	  $result->{$_} = 1;
	}
      } else {
	$result->{$_} = 1;
      }
    }
  } else {
    my @results;
    while (<RDTMP>) {
      chomp;
      s/\r$//;
      # DOS line ending
      push @results, $_;
    }
    close RDTMP;
    $result = \@results;
  }
  return $result;
}

sub sort_files_by_date {
  my (%options) = @_;
  my $files = $options{"-files"} || die "-files";
  my $type = $options{"-type"} || "newest";
  
  my %mtime = map {$_, (stat($_))[9]} @{$files};

  my @results;
  if ($type eq "newest") {
    @results = sort {$mtime{$b} <=> $mtime{$a}} @{$files};
  } else {
    die;
  }

  if ($options{"-single"}) {
    return $results[0];
  } else {
    return \@results;
  }
}

sub get_local_arch_path {
  #
  # get architecture-dependent local path
  #
  # FIX ME: use $ENV{MACHTYPE} instead if available??
  #
  my (%options) = @_;
  my ($sysname, $nodename, $release, $version, $machine) = uname();
  my $root = $options{"-root"} || $ENV{HOME} || "";
  my $path = sprintf "%s/local/%s_%s/", $root, ($sysname || ""), $machine;
  unless (-d $path) {
    $path = sprintf "%s/local/%s_%s/", "/h1/edmonsom", ($sysname || ""), $machine;
  }
  return $path;
}

sub get_local_arch_binary {
  my ($program, %options) = @_;
  return get_local_arch_path(%options) . "/bin/" . $program;
}

sub decompress_gzip {
  # uncompress a gzip-format file.
  my (%options) = @_;
  my $file = $options{"-file"} || die "need -file";
  my $gz = gzopen($file, "rb") || die;

  my $out_fn = $file;
  $out_fn =~ s/\.gz$//i;
  die if $out_fn eq $file;

  my $wf = new WorkingFile($out_fn);
  my $fh = $wf->output_filehandle;

  my $buffer;
  my $status;
  while (($status = $gz->gzread($buffer)) > 0) {
    print $fh $buffer;
  }
  $gz->gzclose;
  $wf->finish;

  if ($status == 0) {
    # file ended normally
    return 1;
  } else {
    die "gz error!";
  }
}

sub get_df {
  #
  # system call, bleah;
  # Filesys::DiskUsage doesn't seem to work on this machine
  #
  my ($dir) = @_;
  open(DF, sprintf 'df -k %s|', $dir) || die;
  my $header = <DF>;
  my $l1 = <DF>;
  my $l2 = <DF>;
  close DF;
  chomp $header;
  chomp $l1;
  chomp $l2;
  $header =~ s/Mounted on/Mounted_on/ || die;
  my %df;
  my @labels = split /\s+/, $header;
  @df{@labels} = split /\s+/, $l2;
  $df{$labels[0]} = $l1;
  return \%df;
}

sub mkpath_cd {
  # chdir into a directory, mkpath'ing if it doesn't exist.
  my ($dir) = @_;
  unless (-d $dir) {
    mkpath($dir) || die "mkpath $dir";
  }
  chdir($dir) || die "chdir $dir";
}

sub mtime {
  my ($fn) = @_;
  return (stat($fn))[9];
}

sub test_magic {
  #
  # test files for consistency with known magic numbers
  #
  my (%options) = @_;
  my $ok;
  my $format = $options{"-format"} || die "-format";
  my $file = $options{"-file"} || die "-file";

  open(MAGICTMP, $file) || confess "can't open $file";
  my $buf;
  read(MAGICTMP, $buf, 2);
  close MAGICTMP;

  my $hex = unpack "H*", $buf;

  if ($format eq "gzip") {
    $ok = $hex eq "1f8b" ? 1 : 0;
  } else {
    die "unhandled format $format";
  }

  return $ok;
}

sub find_binary {
  my ($name, %options) = @_;
  my $result;
  foreach my $dir (split /:/, $ENV{PATH} || "") {
    my $prog = $dir . "/" . $name;
    if (-f $prog and -s $prog) {
      if (-x $prog) {
	$result = $prog;
	last;
      } else {
	printf STDERR "warning: %s is not executable\n", $prog;
	$result = $prog;
	last;
	# accept anyway: in ctest, this works even though not executable
      }
    }
  }

  my %bin2module;
  $bin2module{"blat"} = "blat";
  $bin2module{"blastall"} = "blast";
  $bin2module{"gfClient"} = "blat";
  $bin2module{"tabix"} = "tabix";
  $bin2module{"bgzip"} = "tabix";
  $bin2module{"bcftools"} = "bcftools";
  $bin2module{"indelmatch"} = "indelmatch";
  $bin2module{"variant_effect_predictor.pl"} = "vep";

  my %cbbin2module;
  foreach (qw(snv_high_20_tn.sh snv_low_tn.sh snv_germline.sh)) {
    $cbbin2module{$_} = "snv";
  }
  $cbbin2module{"indel_query.pl"} = "variants2matrix";
  $cbbin2module{"pmid2abstract.pl"} = "common-scripts-internal";

  if (($options{"-die"} or $options{"-warn"}) and !$result) {
    my $msg = sprintf "ERROR: can't find %s on PATH; missing required module load or app profile?", $name;
    my $module = $bin2module{$name};
    $msg .= sprintf ' Try running this first: "module load %s"', $module if $module;

    my $cb = $cbbin2module{$name};
    $msg .= sprintf ' Try running this first: "cbload %s"', $cb if $cb;

    if ($options{"-die"}) {
      confess $msg;
    } else {
      printf STDERR "%s\n", $msg;
    }
  }

  return $result;
}

sub count_file_lines {
  my ($fn) = @_;
  open(LCTMP, $fn) || die;
  my $count = 0;
  while (<LCTMP>) {
    $count++;
  }
  close LCTMP;
  return $count;
}

sub universal_open {
  # transparently open regular or compressed files.
  # copied from rude misc.pm
  my ($fn) = @_;
  my $fmt;
  if ($fn =~ /\.bz2$/i) {
    # compressed w/bzip2
    $fmt = find_binary("bzip2") . ' -dc %s|';
  } elsif ($fn =~ /\.g?z$/i) {
    # compress/gzip
    $fmt = find_binary("gzip") . ' -dc %s|';
  } else {
    $fmt = '%s';
  }
  my $fh = new FileHandle();
  $fh->open(sprintf $fmt, $fn) || die "universal_open of $fn failed: $!";
  return $fh;
}

sub aggressive_rename {
  # if regular rename() fails, try unlink/rename instead.
  # (copied from rude misc.pm)
  my ($from, $to, %options) = @_;
  confess "WTF: no from file $from" unless $from and -f $from;
  confess "WTF: no to file" unless $to;
  my $verbose = 1;
  my $errors_fatal = 1;
  $verbose = $options{"-verbose"} if exists $options{"-verbose"};
  $errors_fatal = $options{"-fatal"} if exists $options{"-fatal"};
  my $skip_if_exists = $options{"-skip-if-exists"};
  my $error;

  if ($skip_if_exists and -s $to) {
    unlink $from;
    return 1;
  } elsif (rename($from, $to)) {
    return 1;
  } else {
#    print STDERR "rename $from => $to failed!\n";
    if (-f $to) {
      # target filename exists
      if (unlink $to) {
	if (rename($from, $to)) {
	  print STDERR "successful unlink-style rename for $to\n" if $verbose;
	} else {
	  $error = "doh! rename of $from failed after $to unlink!";
	}
      } else {
	$error = "argh, can't unlink $to";
      }
    } else {
      $error = "rename failed and $to doesn't exist";
    }
  }

  if ($error) {
    if ($errors_fatal) {
      confess $error;
    } else {
      print STDERR "$error\n" if $verbose;
      return undef;
    }
  } else {
    return 1;
  }
}

sub newer_than {
  # is file1 newer (modified more recently than) file2?
  my ($f1, $f2) = @_;
  confess "undefined f1" unless defined $f1;
  confess "undefined f2" unless defined $f2;
  my $f1_mt = -f $f1 ? (stat($f1))[9] : 0;
  my $f2_mt = -f $f2 ? (stat($f2))[9] : 0;
  return $f1_mt > $f2_mt ? 1 : 0;
}

sub line_delimiter_qc {
  my (%options) = @_;
  my $file = $options{"-file"} || die "-file";
  my $delimiter = $options{"-delimiter"} || die "-delimiter";

  my $fh = universal_open($file) || die "can't open $file: $!";

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
  $fh->close();

  if (scalar keys %counts > 1) {
    # FAIL
    printf STDERR "ERROR in %s: multiple column counts!:\n", $file;
    foreach my $count (sort {$a <=> $b} keys %counts) {
      printf STDERR "  cols=%d row_count=%d example=%s\n",
      $count, $counts{$count}, $example{$count};
    }
    confess;
  }
}

sub open_filehandle {
  my ($fn) = @_;
  my $fh = new FileHandle();
  $fh->open($fn) || die "can't open $fn: $!";
  return $fh;
}

sub md5_file {
  my ($fn) = @_;
#  printf STDERR "MD5 %s\n", $fn;
  open(MD5FNTMP, $fn) || confess "can't open $fn";

  my $md5 = new Digest::MD5();
  $md5->addfile(*MD5FNTMP);
  close MD5FNTMP;
  return $md5->hexdigest();
}

sub is_offline {
  my ($fn_raw, %options) = @_;
  my $fn = abs_path($fn_raw);
  # resolve to final file if symlink
  my $type = $options{"-type"} || 0;

  my %offline;
  $offline{0} = 1;
  $offline{2048} = 1;
  delete $offline{2048} if $type == 2;
  # continue if desired in this case, file may be recallable on access

  my $chunk = `stat $fn`;
  my $is_offline = 0;
  if ($chunk =~ /Blocks: (\d+)/) {
    $is_offline = 1 if $offline{$1};
  } else {
    die "can't find Blocks: in stat for $fn";
  }
  return $is_offline;
}

sub today_file {
  my ($template) = @_;
  die "$template must contain %s" unless $template =~ /%s/;
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
#my $outfile = $FLAGS{out} || basename($f_ranges) . ".patched.txt";
  my $date = sprintf "%d_%02d_%02d", 1900+$year, $mon + 1, $mday;
  return sprintf $template, $date;
}

sub tab_cat {
  my (%options) = @_;
  my $infiles = $options{"-files"} || die "-files";
  my $outfile = $options{"-out"} || die "-out";
  my $remove_dups = $options{"-remove-duplicates"};

  my $wf = new WorkingFile($outfile);
  my $fh_out = $wf->output_filehandle();
  my $first_headers;

  my $empty_files = 0;
  my $header_first;
  my %saw;
  foreach my $f_in (@{$infiles}) {
    unless (-s $f_in) {
      $empty_files++;
      next;
    }
    printf STDERR "%s...\n", $f_in;

    my $fh = new FileHandle();
    $fh->open($f_in) || die "can't open $f_in: $!";
    my $header_this = <$fh>;
    if ($header_first) {
      die "header mismatch" unless $header_this eq $header_first;
    } else {
      $header_first = $header_this;
      print $fh_out $header_this;
    }

    if ($remove_dups) {
      while (<$fh>) {
	print $fh_out $_ unless $saw{$_};
	$saw{$_} = 1;
      }
    } else {
      while (<$fh>) {
	print $fh_out $_;
      }
    }
  }
  $wf->finish();
  printf STDERR "skipped %d empty files\n", $empty_files if $empty_files;
}

1;
