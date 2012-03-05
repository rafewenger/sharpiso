#!/usr/bin/perl
# test and compare different versions of anisograd
# updated by R. Wenger Mar 5, 2012

use strict;

my $testdir = "cylinder";

my @proglist = @ARGV;
my @input_options;
my $fastflag = 0;
my $veryfastflag = 0;
my $alloptions = 0;
my $found_difference = 0;
my $outfile = "temp.nrrd";
my $outfile0 = "temp0.nrrd";

# anisograd arguments which take an input value/string.
my @anisograd_options = ( );

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

  if ($new_option eq "-alloptions") {
    $alloptions = 1;
    next;
  }

  push(@input_options, $new_option);

  if (scalar(@proglist) > 0) {
    if (is_anisograd_option("$new_option")) {
      $new_option = shift(@proglist);
      push(@input_options, $new_option);
    }
  }
}

if (scalar(@proglist) < 1) { usage_error(); };

my %testdata;

$testdata{curveA_11x}{fname} = "curveA.11x.nrrd";



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
 print "Usage: testanisograd.pl [OPTIONS] {prog1} {prog2} [{prog3} ...]";
 exit(10);
}

# compare executables with all possible options
sub compare_executables_all_options {

  compare_executables("-position centroid");
  print "\n";
}

sub compare_executables {

  my @option_list = @_;

  foreach my $tdata (keys %testdata) {

    my $tfile = "$testdir"."/"."$testdata{$tdata}{fname}";
    run_anisograd($prog0, $tfile, "$outfile0", \@option_list);
    foreach my $anisograd (@proglist) {
      run_anisograd($anisograd, $tfile, "$outfile", \@option_list);
      diff_files("$outfile0", "$outfile");
    }
  }
}

sub run_anisograd {

  scalar(@_) > 3 ||
    die "Error in sub run_anisograd. Requires at least 4 parameters.\n";

  my $anisograd = $_[0];
  my $input_filename = $_[1];
  my $output_filename = $_[2];
  my @option_list = @{$_[3]};

  my $command_line;
  $command_line = 
    "./$anisograd @option_list $input_filename $output_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program anisograd abnormally terminated.\n";

}

sub diff_files {

  scalar(@_) == 2 ||
    die "Error in sub run_anisograd. Requires exactly 2 parameters.\n";

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

# return true (1) if $_[0] is an anisograd option which takes an 
#   input value/string
sub is_anisograd_option {

  scalar(@_) == 1 ||
    die "Error in sub is_anisograd_option. Requires exactly 1 parameter.\n";

  foreach my $option (@anisograd_options) {
    if ($_[0] eq "$option") {
      return 1;
    }
  };

  return 0;
}


