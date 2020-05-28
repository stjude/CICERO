package SampleName;

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use Data::Compare;
use Carp qw(confess cluck);
# use Log::Log4perl qw(:easy);

$Carp::MaxArgLen = 0;   # Show full length of function argument. Default is 64.
$Carp::MaxArgNums = 0;  # Show all arguments. Default is 8.


#    parse(<name>,<error behavior>)
# 
#usages  
#    parse("SJTALL001_D-TB-01-0284","WARN")
#    parse("SJTALL001_D-TB-01-0284","DIE")
sub parse {
  my ($sampname, $fb) = @_;
  my  $failbehavior = uc($fb);
  my %h = ();
  my $sid;
  if ($sampname =~ /^SJ([A-Za-z0-9]*[A-Za-z]\d?)(\d{3})_([A-Z])(-(.+))?$/) {
    # type 1
    %h = ("subject" => "SJ" . $1 .$2,
          "disease" => $1,
# ADT: 1/8/2015 - This behavior is inconsistent with the documentation and is not the desired functionality. 
# The samplename library should return the generic type code rather than the specific letter for V1 sample names. 
#    "type" => $3,  # actual type is wanted here for S, H, etc.
          "type" => convertType($3),
          "index" => getIndex($3),
          "sid" => undef,
          "number" => $2,
          "sample" => "SJ".$1.$2."_".$3,
          "barcode" => undef,
          "conventionVersion" => 1
         );
    if ($h{"index"} < 1) {
      if ($failbehavior eq "WARN") {
    warn "ERROR: Unknown sample type: ".$h{"type"}."\n";
    return my %error;
      } elsif ($failbehavior eq "DIE") {
    die "FATAL: Unknown sample type: ".$h{"type"}."\n";
      } else {
    die "FATAL: Unknown ERROR behavior: $failbehavior\nERROR: Unknown sample type: ".$h{"type"}."\n";
      }
    }
    if (defined($4)) {$sid = $5;}
  } elsif ($sampname =~ /^SJ([A-Za-z0-9]*[A-Za-z]\d?)(\d{6})_([A-Z])(\d+)(-(.+))?$/) {
    # type 2
    if ($4 =~ /^0/) {
      if ($failbehavior eq "WARN") {
    warn "ERROR: Padded index:$sampname\n";
    return my %error;
      } elsif ($failbehavior eq "DIE") {
    die "FATAL: Padded index:$sampname\n";
      } else {
    die "FATAL: Unknown ERROR behavior: $failbehavior\n Padded index:$sampname\n";
      }
    }
    %h = ("subject" => "SJ" . $2,
          "disease" => $1,
          "type" => $3,
          "index" => $4,
          "sid" => undef,
          "number" => $2,
          "sample" => "SJ".$1.$2."_".$3.$4,
          "barcode" => undef,
          "conventionVersion" => 2
         );
    if (checkType($3) == 0) {
      if ($failbehavior eq "WARN") {
    warn "ERROR: Unknown sample type: ".$h{"type"}."\n";
    return my %error;
      } elsif ($failbehavior eq "DIE") {
    die "FATAL: Unknown sample type: ".$h{"type"}."\n";
      } else {
    die "FATAL: Unknown ERROR behavior: $failbehavior\nERROR: Unknown sample type: ".$h{"type"}."\n";
      }
    }
    if (defined($5)) {$sid = $6;}
  } else {
    if ($failbehavior eq "WARN") {
      warn "ERROR: Incorrectly formatted sample name: $sampname\n";
      return %h;
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: Incorrectly formatted sample name: $sampname\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\nIncorrectly formatted sample name: $sampname\n";
    }
  }
  if (defined($sid)) {
    $h{"sid"} =  $sid;
  }
    $h{"barcode"} = $sampname;
  return %h;
}


sub parseBamName {
  my ($filename, $fb) = @_;
  $filename =~ s/\_\d\d\d\d\d\d\d\d\d\.bam// || $filename =~ s/\_\d\d\d\d\d\d\d\d\dM\.bam// || $filename =~ s/\.bam//;
  return parse($filename, $fb);
}


#    parseDisease(<name>,<error behavior>)
# 
#usages  
#    parseDisease("SJTALL001_D-TB-01-0284","WARN")
#    parseDisease("SJTALL001_D-TB-01-0284","DIE")
#
# Convenience method to return the disease only
sub parseDisease {
  my %h = parse(@_);
  return $h{"disease"};
}


sub build {
  my %sampname = (%{shift(@_)});
  my $fb = shift(@_);
  my $failbehavior = uc($fb);
  my $full = 0;  # SAMPLE name only
  my $name;
  if (defined($fb)) { $failbehavior = uc($fb); } else { $failbehavior = 'DIE'; }
  
  if (checkHash(\%sampname, $failbehavior)  == 0) { return ""; }
  if (defined($sampname{"sid"})) { $full = 1; }  # full BARCODE

  if ($sampname{"conventionVersion"} == 1) {
    $name = buildV1SampleFromHash(\%sampname, $failbehavior, $full);
  } else {
    $name = buildV2SampleFromHash(\%sampname,  $full);
  }
  return $name;
}



sub checkHash {
  my %sampname = (%{shift(@_)});
  my $failbehavior = shift(@_);
  my $disease;
  my $number;


  # must have "conventionVersion" of 1 or 2
  if (!(defined($sampname{"conventionVersion"}) && checkConventionVersion($sampname{"conventionVersion"}))) {   
    if ($failbehavior eq "WARN") {
      warn "ERROR:  no or incorrect conventionVersion provided.  Should be 1 or 2.\n";
      return 0;
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: no or incorrect conventionVersion provided. Should be 1 or 2.\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\nno incorrect conventionVersion provided.  Should be 1 or 2.\n";
    }
  }
  
  # must have valid "subject"
  if (defined($sampname{"subject"})) {
    my $badsubject = 1;
    if ($sampname{"conventionVersion"} == 1) {
      if ($sampname{"subject"} =~ /^SJ([A-Za-z0-9]+)(\d\d\d)$/) { $badsubject = 0; $disease = $1; $number = $2;}
    } else {
      if ($sampname{"subject"} =~ /^SJ(\d\d\d\d\d\d)$/) { $badsubject = 0;  $number = $1;}
    }
    if ($badsubject) {
      if ($failbehavior eq "WARN") {
    warn "ERROR: bad subject value provided: ".$sampname{"subject"}."\n";
    return 0;
      } elsif ($failbehavior eq "DIE") {
    die "FATAL: bad subject value provided: ".$sampname{"subject"}."\n";
      } else {
    die "FATAL: Unknown ERROR behavior: $failbehavior\nbad subject value provided: ".$sampname{"subject"}."\n";
      } 
    } 
  } else {
    if ($failbehavior eq "WARN") {
      warn "ERROR: no subject value provided in hash\n";
    return 0;
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: no subject value provided in hash\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\nno subject value provided in hash\n";
    }
  }

  # subject disease must match disease
  if ($sampname{"conventionVersion"} == 1) {
    if (defined($sampname{"disease"})) {
      if ($sampname{"disease"} ne $disease) {
    if ($failbehavior eq "WARN") {
      warn "ERROR: provided disease does not match subject id: ".$sampname{"disease"}."\n";
      return 0;
    } elsif ($failbehavior eq "DIE") {
      die "FATAL:  provided disease does not match subject id: ".$sampname{"disease"}."\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\n provided disease does not match subject id: ".$sampname{"disease"}."\n";
    }  
      }
    } else {
      if ($failbehavior eq "WARN") {
    warn "ERROR: no disease value provided in hash\n";
    return 0;
      } elsif ($failbehavior eq "DIE") {
    die "FATAL: no disease value provided in hash\n";
      } else {
    die "FATAL: Unknown ERROR behavior: $failbehavior\nno disease value provided in hash\n";
      }
    }
  }
  
  # subject number must match number
  if (defined($sampname{"number"})) {
    if ($sampname{"number"} ne $number) {
      if ($failbehavior eq "WARN") {
    warn "ERROR: provided number does not match subject id: ".$sampname{"number"}."\n";
    return 0;
      } elsif ($failbehavior eq "DIE") {
    die "FATAL: provided number does not match subject id: ".$sampname{"number"}."\n";
      } else {
    die "FATAL: Unknown ERROR behavior: $failbehavior\n provided number does not match subject id: ".$sampname{"number"}."\n";
      }  
    }
  } else {
    if ($failbehavior eq "WARN") {
      warn "ERROR: no number value provided in hash\n";
      return 0;
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: no number value provided in hash\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\nno number value provided in hash\n";
    }
  }
  
  # must have valid type
  if (checkType($sampname{"type"}) == 0) {  
    if ($failbehavior eq "WARN") {
      warn "ERROR: bad type value provided: ".$sampname{"type"}."\n";
      return 0;
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: bad type value provided: ".$sampname{"type"}."\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\nbad type value provided: ".$sampname{"type"}."\n";
    }  
  }

  # must have index provided
  if ((!(defined($sampname{"index"}))) || $sampname{"index"} =~ /\D+/ || $sampname{"index"} < 1) {
    if ($failbehavior eq "WARN") {
      warn "ERROR: incorrect or no index provided\n";
      return 0;
    } elsif ($failbehavior eq "DIE") {
      die "FATAL:  incorrect or no index provided\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\nincorrect or no index provided\n";
    } 
  }
  return 1;
}





sub pairToString {
  my ($casename, $controlname, $fb) = @_;
  my  $failbehavior = uc($fb);

  # in case they call with just on sample or something
  if (scalar @_ < 2) {
    die "Improper usage of pairToString(<casename>,<controlname>,<DIE|WARN>)\n";  # fb is not defined.
  }

  if ($casename eq $controlname) {
    if ($failbehavior eq "WARN") {
      warn "ERROR: case and control have the same name $casename\n";
      return ("");
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: case and control have the same name $casename\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\n case and control have the same name $casename\n";
    } 
  }

  my %case = parse($casename, $fb);
  my %control = parse($controlname, $fb);
  my $pair;

  my %empty;
  if (Compare(\%empty,\%case)) { return (""); }
  if (Compare(\%empty,\%control)) { return (""); }
  
  # both v1
     # same disease
     # different disease
  # both v2
     # same disease
     # different disease
  # mix versions 
     # same disease
     # different disease
  
  if ($case{"conventionVersion"} == 1 && $control{"conventionVersion"} == 1) {
    if ($case{"disease"} eq $control{"disease"}) {
      if ($case{"number"} eq $control{"number"}) {
        $pair = $case{"subject"} . "_" . extendType($case{"type"},$case{"index"},$fb) . "_" . extendType($control{"type"},$control{"index"},$fb);
      } else {  # different donors from same disease
        $pair = $case{"subject"} . "_" . extendType($case{"type"},$case{"index"},$fb)  . "_" . $control{"subject"} . "_" . extendType($control{"type"},$control{"index"},$fb);
      }
    } else {   # mixed diseases
      $pair = $case{"subject"} . "_" . extendType($case{"type"},$case{"index"},$fb)  . "_" . $control{"subject"} . "_" . extendType($control{"type"},$control{"index"},$fb);
    }
  } elsif ($case{"conventionVersion"} == 2 && $control{"conventionVersion"} == 2) {
    if ($case{"disease"} eq $control{"disease"}) {
      if ($case{"number"} eq $control{"number"}){
        $pair = "SJ"  . $case{"disease"} . $case{"number"} . "_" .   $case{"type"} . $case{"index"} . "_" . $control{"type"} . $control{"index"};
      } else { # Different donors from same disease
        $pair = "SJ"  . $case{"disease"} . $case{"number"} . "_" . $case{"type"} . $case{"index"}. "_" . "SJ"  . $control{"disease"} . $control{"number"}  . "_" . $control{"type"} . $control{"index"};
      }
    } else {  # mixed diseases
      $pair = "SJ"  . $case{"disease"} . $case{"number"} . "_" . $case{"type"} . $case{"index"}. "_" . "SJ"  . $control{"disease"} . $control{"number"}  . "_" . $control{"type"} . $control{"index"};
    }
  } elsif ($case{"conventionVersion"} == 1 && $control{"conventionVersion"} == 2) {    ## mixed versions
      $pair = $case{"subject"} . "_" . extendType($case{"type"},$case{"index"},$fb)  . "_" . "SJ"  . $control{"disease"} . $control{"number"} . "_" .   $control{"type"} . $control{"index"} ;
  } elsif ($case{"conventionVersion"} == 2 && $control{"conventionVersion"} == 1) {    ## mixed versions
    $pair = "SJ"  . $case{"disease"} . $case{"number"} . "_" .   $case{"type"} . $case{"index"} . "_". $control{"subject"} . "_" . extendType($control{"type"},$control{"index"},$fb);
  }
  return $pair;
}


sub stringToPair {
  my ($pair, $fb)  = @_;
  my  $failbehavior = uc($fb);
  my $caseV;
  my $controlV;
  my $case = "";
  my $control = "";
  my $shortform;
  
  # check for only expected chars
  if (!($pair =~ /^[A-Z0-9\_]+$/)) {
    if ($failbehavior eq "WARN") {
      warn "ERROR: incorrectly formatted pair name, $pair\n";
      return ("","");
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: incorrectly formatted pair name, $pair\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\n incorrectly formatted pair name, $pair\n";
    }  
  }

  # check both v1, same subject
  my @x = split(/_/, $pair);
  if (scalar @x == 3) {
    $shortform = 1;
  } elsif (scalar @x == 4) {
    $shortform = 0;
    # check that it should not be shortform
    if ($x[0] eq $x[2]) {
      if ($failbehavior eq "WARN") {
        warn "ERROR: $x[0] eq $x[2], should use short form\n";
        return ("","");
      } elsif ($failbehavior eq "DIE") {
    #die "FATAL:  $x[0] eq $x[2], should use short form\n";
    warn "ERROR:  $x[0] eq $x[2], should use short form\n";
      } else {
        #die "FATAL: Unknown ERROR behavior: $failbehavior\n $x[0] eq $x[2], should use short form\n";
        warn "ERROR: Unknown ERROR behavior: $failbehavior\n $x[0] eq $x[2], should use short form\n";
      }
    }
  } else {
    if ($failbehavior eq "WARN") {
      warn "ERROR: incorrectly formatted pair name, $pair, wrong number of parts\n";
      return ("","");
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: incorrectly formatted pair name, $pair, wrong number of parts\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\n incorrectly formatted pair name, $pair, wrong number of parts\n";
    }  
  }

  if ($shortform) {
    my $subject  = $x[0];
    my $typeCase = $x[1];
    my $typeControl = $x[2];
    $case = $subject . "_" . $typeCase;
    $control = $subject . "_" . $typeControl;
  } else {
    $case = $x[0] . "_" . $x[1];
    $control = $x[2] . "_" . $x[3];
  }
  
  # make sure the built-up samples parse well
  my %empty;
  my %parseCase = parse($case, $failbehavior);
  if (Compare(\%empty,\%parseCase)) { return ("",""); }
  my %parseControl = parse($control, $failbehavior);
  if ( Compare(\%empty,\%parseControl)) { return ("",""); }

  return ($case, $control);
}


sub getIndex {
  my ($type) = @_;
  my $i;
  if ($type eq "Y" || $type eq "S" || $type eq "E" || $type eq "B" || $type eq "H") {
    $i = 2;
  } elsif ($type eq "Z"  || $type eq "F" || $type eq "T" || $type eq "I") {
    $i = 3;
  } elsif ($type eq "G" || $type eq "D" || $type eq "X" || $type eq "A" || $type eq "M" || $type eq "R" || $type eq "C") {
    $i = 1;
  } else {
    $i = 0;
  }
  return $i;
}

sub checkType {
  my ($type) = @_;
  if ($type eq "G" || $type eq "D" || $type eq "X" || $type eq "A" || $type eq "M" || $type eq "R" || $type eq "C" || $type eq "O") {
    return 1;
  }
  return 0;
}


sub checkConventionVersion {
  my ($x) = @_;
  if ($x == 1 || $x == 2) {
    return 1;
  }
  return 0;
}


sub convertType {
  my ($type) = @_;
  if ($type eq "Y" || $type eq "Z") {
    $type ='X';
  } elsif ($type eq "S" || $type eq "T") {
    $type = 'R';
  } elsif ($type eq "E" || $type eq "F" ) {
    $type = 'D';
  } elsif ($type eq "B" ) {
    $type = 'A';
  } elsif ($type eq "H" || $type eq "I") {
    $type = 'G';
  }  
  return $type;
}


sub buildV1SampleFromHash {
  my %sampname = (%{shift(@_)});
  my $failbehavior = shift(@_);
  my $isbarcode = shift(@_);
  my $newtype = extendType($sampname{"type"}, $sampname{"index"}, $failbehavior);
  if ($newtype eq "") { return ""; }
  my $n = $sampname{"subject"} . "_" . $newtype;
  if ($isbarcode) {
    $n .=  "-" . $sampname{"sid"};
  }
  return $n;
}

sub buildV2SampleFromHash {
  my %sampname = (%{shift(@_)});
  my $isbarcode = shift(@_);
  my $n = "SJ"  . $sampname{"disease"} . $sampname{"number"} . "_" . $sampname{"type"} . $sampname{"index"};
  if ($isbarcode) {
    $n .=  "-" . $sampname{"sid"};
  }
  return $n;
}


sub extendType {
  my ($type, $index, $failbehavior) = @_;
  my $badindex = 0;
  if ($type eq 'X') {
    if ($index == 2) {
      $type = "Y";
    } elsif  ($index == 3) {
      $type = "Z";
    } elsif  ($index == 1) {
      ## leave type as is
    } else {
      $badindex = 1;
    }
  } elsif ($type eq 'R') {
    if ($index == 2) {
      $type = "S";
    } elsif ($index == 3) {
      $type = "T";
    } elsif  ($index == 1) {
      ## leave type as is
    } else {
      $badindex = 1;
    } 
  } elsif ($type eq 'D') {
    if ($index == 2) {
      $type = "E";
    } elsif  ($index == 3) {
      $type = "F";
    } elsif  ($index == 1) {
      ## leave type as is
    } else {
      $badindex = 1;
    }
  } elsif ($type eq 'A') {
    if ($index == 2) {
      $type = "B";
    } elsif  ($index == 1) {
      ## leave type as is
    } else {
      $badindex = 1;
    }
  } elsif ($type eq 'G') {
    if ($index == 2) {
      $type = "H";
    } elsif  ($index == 3) {
      $type = "I";
    } elsif  ($index == 1) {
      ## leave type as is
    } else {
      $badindex = 1;
    }
  }

  if ($badindex) {
    if ($failbehavior eq "WARN") {
      warn "ERROR: bad index ($index) for type ($type)\n";
      return "";
    } elsif ($failbehavior eq "DIE") {
      die "FATAL: bad index ($index) for type ($type)\n";
    } else {
      die "FATAL: Unknown ERROR behavior: $failbehavior\nbad index ($index) for type ($type)\n";
    } 
  }
  return $type;
}

sub getFullType {
  my ($type2, $fb) = @_;
  my  $failbehavior = uc($fb);
  my $ret;

  my $type = &convertType($type2);

  if ($type eq "G") {
    $ret = "GERMLINE";
  } elsif ($type eq "D") {
    $ret = "DIAGNOSIS";
  } elsif ($type eq "X") {
    $ret = "XENOGRAFT";
  } elsif ($type eq "A") {
    $ret = "AUTOPSY"
  } elsif ($type eq "M") {
    $ret = "METASTATIC";
  } elsif ($type eq "R") {
    $ret = "RELAPSE"
  } elsif ($type eq "C") {
    $ret = "CELL_LINE";
  } else {
    if ($failbehavior eq "DIE") {
      die "FATAL: Unknown single sample type code: $type\n";
    } else {
      warn "ERROR: Unknown single sample type code: $type\n";
      return "";
    }
  }
  return $ret;
}


1;
