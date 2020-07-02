package TabixFile;
# tabix-indexed file

use strict;
use FileHandle;
use Exporter;
use Carp qw(confess cluck);
use File::Copy;

use Vcf;
use Configurable;
use ReferenceNameMapper;
use TemporaryFileWrangler;
use WorkingFile;

use FileUtils qw(find_binary universal_open);
use MiscUtils qw(dump_die split_list);

use constant MAX_INTERVALS_PER_QUERY => 200;
#use constant MAX_INTERVALS_PER_QUERY => 500;
# how many is too many?

@TabixFile::ISA = qw(Configurable Exporter);
@TabixFile::EXPORT_OK = qw();

use MethodMaker qw(
	file
        binary_name
        preset
        need_headers
        wiggle_bases
rnm
warned_chr

all_header_lines
headers
warn_missing_reference
index
f_chr
f_start
f_end

sort_intervals
indel_wiggle_bases
tfw
verbose

suppress_duplicate_intervals
suppress_duplicate_output

vcf2tab_mode
vcf
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->binary_name("tabix");
  $self->suppress_duplicate_intervals(1);
  # e.g. dbNSFP batch tabix for 2 variants at same site (chr4.106190863.G.[AT])
  # querying the same interval twice returns duplicate records, which
  # cause problems with annotation as ambiguity is detected
  $self->suppress_duplicate_output(1);
  # in case intervals overlap, don't want to return same records twice
  $self->warn_missing_reference(1);
  $self->warned_chr({});
  $self->need_headers(1);
  $self->sort_intervals(1);
  $self->tfw(new TemporaryFileWrangler());
#  $self->verbose(1);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $file = $self->file() || die "-file";
  my $exe = $self->binary_name();
  find_binary($exe, "-die" => 1);
  # ensure tabix on path
  # TO DO: get official tabix perl module installed
  confess "where is $file" unless -s $file;

  my $tbi = $file . ".tbi";
  if ($self->index() and not(-s $tbi)) {
    my $f_chr = $self->f_chr || die "f_chr";
    my $f_start = $self->f_start || die "f_start";
    my $f_end = $self->f_end;
    my $fh = universal_open($file);
    my $first = <$fh>;
    chomp $first;
    my @f = split /\t/, $first;
    my $i = 0;
    my %name2num;
    foreach (@f) {
      $name2num{$_} = ++$i;
    }
    my $n_chr = $name2num{$f_chr} || $name2num{"#" . $f_chr} || die "no fnum for $f_chr";
    my $n_start = $name2num{$f_start} || die "no fnum for $f_start";
    my $n_end;
    if ($f_end) {
      $n_end = $name2num{$f_end} || die;
    } else {
      $n_end = $n_start;
    }

    my $cmd = sprintf "tabix -s %d -b %d -e %d %s",
      $n_chr, $n_start, $n_end, $file;
    system $cmd;
    die "error running $cmd: $!" if $?;
    die unless -s $tbi;
  }

  my $fh = new FileHandle();
  $fh->open(sprintf "%s -l %s|", $exe, $file) || die;

  my $rnm = new ReferenceNameMapper();
  # index reference names used in tabix file
  while (<$fh>) {
    chomp;
    $rnm->add_name($_);
  }
  $self->rnm($rnm);
  $fh->close() || die "tabix call error";
}

sub query {
  my ($self, %options) = @_;

  #
  # init raw (user format) list of intervals
  # 
  my @intervals_raw;
  my $variants_raw;
  if (my $is = $options{"-intervals"}) {
    # pre-build list of query variants (unadjusted for fuzzy search)
    push @intervals_raw, @{$is};
  } elsif ($variants_raw = $options{"-variants"}) {
    # array of Variant.pm instances
  } else {
    # single manually-specified variant
    my $user_chr = $options{"-reference"} || die "-reference | -intervals";
#  my $user_pos = $options{"-pos"} || die "-pos";
    my $user_pos = $options{"-pos"};
    # might be querying entire chrom
    my $user_end = $options{"-end"};

    if ($user_pos) {
      my $end = $user_end || $user_pos;
      push @intervals_raw, sprintf '%s:%d-%d', $user_chr, $user_pos, $end;
    } else {
      push @intervals_raw, $user_chr;
    }
  }

  my $hash_mode = $options{"-hash"};
  my $vcf_mode = $options{"-vcf-parser"};
  my $vcf2tab_mode = $self->vcf2tab_mode();
  my $rnm = $self->rnm() || die;

  #
  #  initialize intervals:
  #  - adjust reference names to those used by tabix
  #  - adjust query positions if requested (e.g. for indel overlap)
  #  - store in sortable arrayref
  #
  my @intervals_parsed;
  foreach my $ir (@intervals_raw) {
    my @f = split /:/, $ir;
    my ($user_chr, $user_start, $user_end);
    $user_chr = $f[0];
    if (@f == 2) {
      @f = split /\-/, $f[1];
      die unless @f == 2;
      ($user_start, $user_end) = @f;
    }
    
    my $db_chr = $rnm->find_name($user_chr);
    if ($db_chr) {
      my ($db_start, $db_end);
      my $interval;
      if ($user_start) {
	my $wiggle = $self->wiggle_bases || 0;
	$db_start = $user_start - $wiggle;
	$db_end = $user_end + $wiggle;
	$db_start = 1 if $db_start < 1;
#	$interval = sprintf '%s:%d-%d', $db_chr, $db_start, $db_end;
	# interval corrected for names used by database
      }
      push @intervals_parsed, [ $db_chr, $db_start, $db_end ];
    } else {
      $self->chrom_warn($user_chr);
    }
  }

  if ($variants_raw) {
    my $indel_wiggle_bases = $self->indel_wiggle_bases();
    die "specify indel_wiggle_bases" unless defined $indel_wiggle_bases;
    # require be specified (even 0) so we don't get

    my %saw;
    foreach my $v (@{$variants_raw}) {
      my $user_chr = $v->reference_name;
      my $db_chr = $rnm->find_name($user_chr);
      if ($db_chr) {
	my $wiggle = 0;
	$wiggle = $indel_wiggle_bases if $v->is_insertion or $v->is_deletion or $v->is_complex;
	my $user_start = $v->start;
	my $user_end = $v->end;
	my $db_start = $user_start - $wiggle;
	my $db_end = $user_end + $wiggle;
	# TO DO: apply wiggle to indels only??
	$db_start = 1 if $db_start < 1;

	my $key = join ".", $db_chr, $db_start, $db_end;
	push @intervals_parsed, [ $db_chr, $db_start, $db_end ] unless $saw{$key};
	$saw{$key} = 1;
      } else {
	$self->chrom_warn($user_chr);
      }
    }
  }
  
  if ($self->sort_intervals) {
    # sort intervals by reference name and position
    # for (significantly) better performance
    @intervals_parsed = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @intervals_parsed;
  }

  my @intervals_final;
  # final list in tabix format
  my %is;
  my $suppress_duplicate_intervals = $self->suppress_duplicate_intervals();
  foreach my $i (@intervals_parsed) {
    my $interval;
    if (defined $i->[1]) {
      $interval = sprintf '%s:%d-%d', @{$i};
    } else {
      $interval = $i->[0];
      # full chromosome query
    }
    if ($suppress_duplicate_intervals and $is{$interval}) {
#      die "duplicate interval $interval";
    } else {
      push @intervals_final, $interval;
      $is{$interval} = 1;
    }
  }

  # die "no intervals to query!" unless @intervals_final;
  return [] unless @intervals_final;
  # could happen for single row query, e.g. attempt to add
  # NHLBI annotation on chrY (no data)

  my $interval_sets;
  # sets of intervals that can be handled by a single tabix call
  if (@intervals_final <= MAX_INTERVALS_PER_QUERY) {
    $interval_sets = [ \@intervals_final ];
  } else {
    $interval_sets = split_list(\@intervals_final, MAX_INTERVALS_PER_QUERY);
  }

#  die "multiple chunks required, VCF code can only handle one!" if $vcf_mode and @{$interval_sets} > 1;
  # rewrite to use tempfile??

  my $results;
  my $need_headers = $self->need_headers();

  if ($vcf_mode and
      @{$interval_sets} == 1 and
      @{$interval_sets->[0]} == 1 and
      ($self->file =~ /\.vcf$/i or $self->file =~ /\.gz/)
      ) {
    # single interval: hand off to VCF parser.
    # make sure filename is well-formed because parser looks at file
    # extension to determine whether to gunzip, otherwise it crashes
    # (e.g. on dnanexus where filenames change)
    my $vcf = Vcf->new(
      file => $self->file,
      region => $interval_sets->[0]->[0]
	);
    $vcf->parse_header();
    return $vcf;
  } else {
    my $temp_vcf_fn;
    my $temp_vcf_wf;
    my $temp_vcf_fh;

    if ($vcf_mode or $vcf2tab_mode) {
      $temp_vcf_fn = $self->tfw->get_tempfile("-append" => "vcf_temp");
#      printf STDERR "tempfile: %s\n", $temp_vcf_fn;
      $temp_vcf_wf = new WorkingFile($temp_vcf_fn);
      $temp_vcf_fh = $temp_vcf_wf->output_filehandle();
      $need_headers = 1;
    }

    my $verbose = $self->verbose;

    foreach my $intervals_cooked (@{$interval_sets}) {
      printf STDERR "querying %d intervals: %s\n", scalar (@{$intervals_cooked}), join " ", @{$intervals_cooked} if $verbose;
      my @cmd = $self->binary_name;
      push @cmd, "-p", $self->preset() if $self->preset();
      push @cmd, "-h" if $need_headers;
      push @cmd, $self->file || die;
      push @cmd, join " ", @{$intervals_cooked};
      my $cmd = join " ", @cmd;
      printf STDERR "tabix: %s\n", $cmd if $verbose;

      my $fh = new FileHandle();
      $fh->open($cmd . "|") || die "error opening $cmd: $!";

      my @header_lines;
      my @hit_lines;
      while (<$fh>) {
	print $temp_vcf_fh $_ if $vcf_mode or $vcf2tab_mode;
	# in multi-interval VCF mode, copy raw data to tempfile

	chomp;
	if (/^#/) {
	  push @header_lines, $_;
	} else {
	  push @hit_lines, $_;
	}
      }
      $fh->close();

      if ($need_headers) {
	die "where are headers " unless @header_lines;
	$self->all_header_lines([@header_lines]);
	my $header_line = pop @header_lines;

	my @headers = split /\t/, $header_line;
	$headers[0] =~ s/^#//;
	$self->headers(\@headers);
	$self->need_headers(0);
	$need_headers = 0;
	# only parse this with first request
      }

      if ($vcf_mode or $vcf2tab_mode) {
	# already handled
      } elsif ($hash_mode) {
	$results = [] unless $results;
	my $headers = $self->headers() || die;
	foreach my $l (@hit_lines) {
#	printf STDERR "raw: %s\n", $l;
#	my @f = split /\t/, $l;
	  my @f = split /\t/, $l, scalar @{$headers};
	  unless (@f == @{$headers}) {
	    printf STDERR "expected %d, got %d\n", scalar(@{$headers}), scalar @f;
	    for (my $i = 0; $i < @f; $i++) {
	      printf STDERR "  %s: %s\n", $headers->[$i], $f[$i];
	    }
	    die;
	  }

	  my %r;
	  @r{@{$headers}} = @f;
	  push @{$results}, \%r;
	}
      } else {
	$results = [] unless $results;
	push @{$results}, @hit_lines;
      }
    }

    if ($vcf_mode) {
      $temp_vcf_wf->finish();
      if (0) {
	print STDERR "temp VCF: $temp_vcf_fn\n";
	copy($temp_vcf_fn, "debug_filtered.vcf");
      }
      my $vcf = Vcf->new(file => $temp_vcf_fn);
      $vcf->parse_header();
      return $vcf;
    } elsif ($vcf2tab_mode) {
      $temp_vcf_wf->finish();
      my $cmd = sprintf 'vcf2tab.pl -file %s -stdout -sj-post -no-insertion-adjustment -quiet %s|', $temp_vcf_fn, $vcf2tab_mode;
      if (0) {
	printf STDERR "cmd: %s\n", $cmd;
	sleep 1000;
      }
      my $fh_v2t = new FileHandle();
      $fh_v2t->open($cmd) || die "can't open $cmd";
      my $first = <$fh_v2t>;
      chomp $first;
      my @h = split /\t/, $first;
      $self->headers(\@h);

      while (<$fh_v2t>) {
	chomp;
	my %r;
	@r{@h} = split /\t/, $_, -1; 
	push @{$results}, \%r;
      }
      $fh_v2t->close();

    }
  }

  if ($self->suppress_duplicate_output and $results and @{$results}) {
    my $headers = $self->headers;
    my %saw;
    my @out;
    foreach my $r (@{$results}) {
      my $line = join "_", @{$r}{@{$headers}};
      if ($saw{$line}) {
#	die "duplicate";
      } else {
	push @out, $r;
	$saw{$line} = 1;
      }
    }
    $results = \@out;
  }

  return $results;
}

sub get_ref_name {
  my ($self, $user_chr) = @_;
  return $self->rnm->find_name($user_chr);
}

sub chrom_warn {
  my ($self, $user_chr) = @_;
  if (not($self->warned_chr->{$user_chr})) {
    # fix me: 
    cluck sprintf "WARNING: tabix file %s doesn't appear to contain ref seq %s\n", $self->file, $user_chr if $self->warn_missing_reference();
    $self->warned_chr->{$user_chr} = 1;
  }
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
