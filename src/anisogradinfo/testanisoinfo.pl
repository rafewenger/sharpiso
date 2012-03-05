#!/usr/bin/perl
# test and compare different versions of anisogradinfo
# created by R. Wenger: Mar. 5, 2012

use strict;

my ($prog0, $prog1);
my ($found_difference);

my $temp0 = "temp0.out";
my $temp1 = "temp1.out";

parse_command_line(@ARGV);

compare_executables("cylinder/curveA.11x.nrrd", "-vertex 871 -k -c -m -N");

# *********************************************

sub usage_error {
  print "Usage: testanisoinfo.pl {prog1} {prog2}\n";
  exit(10);
}

sub parse_command_line {

  my @command_line = @_;

  if (scalar(@_) != 2) { usage_error(); }

  $prog0 = $command_line[0];
  $prog1 = $command_line[1];
}

sub compare_executables {

  (scalar(@_) == 2) ||
    die "Error in compare_executables.  Requires 2 parameters.\n";

  my $scalar_input = shift;
  my $options = shift;

  run_anisogradinfo($prog0, $scalar_input, $temp0, $options);
  run_anisogradinfo($prog1, $scalar_input, $temp1, $options);
  diff_files($temp0, $temp1);
}

sub run_anisogradinfo {

  (scalar(@_) == 4) ||
    die "Error in run_anisogradinfo.  Requires 4 parameters.\n";

  my $prog = shift;
  my $input_scalar = shift;
  my $output_filename = shift;
  my $options = shift;

  my $command_line =
    "./$prog $options  $input_scalar > $output_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program ijksharpinfo abnormally terminated.\n";
}

sub diff_files {

  scalar(@_) == 2 ||
    die "Error in sub diff_files. Requires exactly 2 parameters.\n";

  my $flag = system("diff --brief $_[0] $_[1] > /dev/null");
  if ($flag != 0) { 
    print "*** Output differs. ***\n\n";
    $found_difference = 1; 
  };

}
