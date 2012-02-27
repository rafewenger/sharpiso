#!/usr/bin/perl
# test and compare different versions of isodual3D
# updated by R. Wenger Feb. 25, 2012

use strict;

my $testdir = "data";

my @proglist = @ARGV;
my @input_options;
my $fastflag = 0;
my $veryfastflag = 0;
my $alloptions = 0;
my $found_difference = 0;
my $isoval_offset = 0;
my $outfile = "temp.off";
my $outfile0 = "temp0.off";

# isodual3D arguments which take an input value/string.
my @isodual3D_options = ( "-subsample",  "-position" );

while (scalar(@proglist) > 0 &&
       $proglist[0] =~ /-.*/) {

  my $new_option = shift(@proglist);

  if ($new_option eq "-fast") { 
    $fastflag = 1; 
    next;
  };

  if ($new_option eq "-veryfast") { 
    $veryfastflag = 1; 
    next;
  };

  if ($new_option eq "-offset") {
    $isoval_offset = shift(@proglist);
    next;
  }

  if ($new_option eq "-alloptions") {
    $alloptions = 1;
    next;
  }

  push(@input_options, $new_option);

  if (scalar(@proglist) > 0) {
    if (is_isodual3D_option("$new_option")) {
      $new_option = shift(@proglist);
      push(@input_options, $new_option);
    }
  }
}

if (scalar(@proglist) < 1) { usage_error(); };

my %testdata;

$testdata{cube_A20}{fname} = "cube3D.A20x.nrrd";
$testdata{cube_A20}{isovalue} = [ 4.9, 5, 5.5 ];

$testdata{cube_B20}{fname} = "cube3D.B20x.nrrd";
$testdata{cube_B20}{isovalue} = [ 5, 5.5];

$testdata{twocubes_A21}{fname} = "twocubes3D.A21x.nrrd";
$testdata{twocubes_A21}{isovalue} = [ 4, 4.5 ];

$testdata{annulus_A31}{fname} = "annulus3D.A31x.nrrd";
$testdata{annulus_A31}{isovalue} = [ 4, 4.5 ];

$testdata{annulus_B31}{fname} = "annulus3D.B31x.nrrd";
$testdata{annulus_B31}{isovalue} = [ 4, 4.5 ];

$testdata{annulus_C31}{fname} = "annulus3D.C31x.nrrd";
$testdata{annulus_C31}{isovalue} = [ 4, 4.5 ];


my $prog0 = shift(@proglist);

if (scalar(@proglist) == 0) {

  usage_error();
}
else {

  if ($alloptions) {
    compare_executables_all_options();
  } 
  else {
    compare_executables(@input_options);
  }

  report_on_differences();
}


if (-e "$outfile0") { system "rm $outfile0"; };
if (-e "$outfile") { system "rm $outfile"; }

# ***********************************

sub usage_error {
 print "Usage: testisodual3D.pl [OPTIONS] {prog1} {prog2} [{prog3} ...]";
 exit(10);
}

# compare executables with all possible options
sub compare_executables_all_options {

  compare_executables("-position centroid");
  print "\n";
  compare_executables("-position cube_center");
  print "\n";
  compare_executables("-position gradC");
  print "\n";
  compare_executables("-position gradCS");
  print "\n";
  compare_executables("-position gradN");
  print "\n";
  compare_executables("-position gradNS");
  print "\n";
}

sub compare_executables {

  my @option_list = @_;

  foreach my $tdata (keys %testdata) {

    my $tfile = "$testdir"."/"."$testdata{$tdata}{fname}";
    my @isovalue_list = @{$testdata{$tdata}{isovalue}};
    foreach my $isoval (@isovalue_list) {

      $isoval = $isoval + $isoval_offset;
      run_isodual3D($prog0, $tfile, "$outfile0", \@option_list, $isoval);
      foreach my $ijkdual (@proglist) {
        run_isodual3D($ijkdual, $tfile, "$outfile", \@option_list, $isoval);
        diff_files("$outfile0", "$outfile");
      }
    }
  }
}

sub run_isodual3D {

  scalar(@_) > 4 ||
    die "Error in sub run_isodual3D. Requires at least 5 parameters.\n";

  my $ijkdual = $_[0];
  my $input_filename = $_[1];
  my $output_filename = $_[2];
  my @option_list = @{$_[3]};
  my $isoval1 = $_[4];

  my $command_line;
  $command_line = 
    "./$ijkdual @option_list -s -o $output_filename $isoval1 $input_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program ijkdual abnormally terminated.\n";

}

sub diff_files {

  scalar(@_) == 2 ||
    die "Error in sub run_isodual3D. Requires exactly 2 parameters.\n";

  my $flag = system("diff --brief $_[0] $_[1] > /dev/null");
  if ($flag != 0) { 
    print "*** Output .off files differ. ***\n\n";
    $found_difference = 1; 
  };

}

sub report_on_differences {

  print "\n";
  if ($found_difference) {
    print "*** Differences detected. ***";
  }
  else {
    print "No differences detected.";
  }

}

# return true (1) if $_[0] is an isodual3D option which takes an 
#   input value/string
sub is_isodual3D_option {

  scalar(@_) == 1 ||
    die "Error in sub is_isodual3D_option. Requires exactly 1 parameter.\n";

  foreach my $option (@isodual3D_options) {
    if ($_[0] eq "$option") {
      return 1;
    }
  };

  return 0;
}


