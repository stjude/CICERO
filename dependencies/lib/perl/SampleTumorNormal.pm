package SampleTumorNormal;
# determine tumor/normal status using either sample naming library
# or hacky pattern recognition
# MNE 10/2015

use strict;
use File::Basename;
use File::Copy;

use Configurable;
use Exporter;
use SampleName;
use MiscUtils qw(dump_die);

@SampleTumorNormal::ISA = qw(Configurable Exporter);

use constant TYPE_NORMAL => "N";
use constant TYPE_TUMOR => "T";

my %CODE2TN = (
  "C" => TYPE_TUMOR,
  # cell line (assuming tumor??)

  "G" => TYPE_NORMAL,
  "H" => TYPE_NORMAL,
  "I" => TYPE_NORMAL,
  # germline

  "D" => TYPE_TUMOR,
  "E" => TYPE_TUMOR,
  "F" => TYPE_TUMOR,
  # primary diagnosis

  "O" => TYPE_TUMOR,
  # other tumor

  "R" => TYPE_TUMOR,
  "S" => TYPE_TUMOR,
  "T" => TYPE_TUMOR,
  # relapse

  "M" => TYPE_TUMOR,
  # metastasis
  "A" => TYPE_TUMOR,
  # autopsy

  "X" => TYPE_TUMOR,
  "Y" => TYPE_TUMOR,
  "Z" => TYPE_TUMOR,
  # xenograft


    );
# http://hc-wiki.stjude.org/display/compbio/Sample+Naming
# how to handle C (cell line)?  Not sure this is guaranteed tumor



use MethodMaker qw(
		    type
		    type_full
		    is_tumor
		    is_normal
		    tn
subject
disease
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub parse {
  my ($self, $thing) = @_;
  $thing = basename($thing);
  # if given fully-qualified BAM file

  my $type = "";
  my $parse_ok = 1;
  my $subject = "";
  my $disease = "";
  my $type_full = "";
  if (my %info = SampleName::parse($thing, "WARN")) {
    # parse OK
    $type = $info{type};
    $subject = $info{subject};
#  } elsif ($thing =~ /^(SJ([A-Z]+)\d+)_([A-Z])/) {
  } elsif ($thing =~ /^(SJ([A-Z\d]+[A-Z])\d+)_([A-Z])(\d*)/) {
    # SJAML040538_D2.bam
    # SJE2A012_D.bam
#    ($subject, $disease, $type) = ($1, $2, $3);
    my $num;
    ($subject, $disease, $type, $num) = ($1, $2, $3, $4);
    $type_full = $num ? $type . $num : $type;
  } elsif ($thing =~ /^Bin([A-Z])\d+_/) {
    # MiSeq bin BAM
    $type = $1;
  } else {
    $parse_ok = 0;
  }
  $self->subject($subject);
  $self->disease($disease);

  my $is_tumor = 0;
  my $is_normal = 0;
  my $tn = "";

  if ($type) {
    $tn = $CODE2TN{$type};
    if ($tn eq TYPE_NORMAL) {
      $is_normal = 1;
    } elsif ($tn eq TYPE_TUMOR) {
      $is_tumor = 1;
    } else {
      die "unhandled type $type";
    }
  }
  $self->type($type);
  $self->type_full($type_full);
  $self->tn($tn);
  $self->is_tumor($is_tumor);
  $self->is_normal($is_normal);

  return $parse_ok;
}

sub is_valid_disease_code {
  my ($self, $code) = @_;
  return $CODE2TN{$code};
}

sub code_is_tumor {
  my ($self, $code) = @_;
  $code =~ s/\d+$//;
  my $tn = $CODE2TN{$code} || die "unknown code $code";
  return $tn eq TYPE_TUMOR;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
