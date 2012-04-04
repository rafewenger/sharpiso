#!/usr/bin/perl
# test and compare different versions of sharpinfo
# created by R. Wenger: Dec. 3, 2011

use strict;

my ($prog0, $prog1);
my ($found_difference);

my $temp0 = "temp0.out";
my $temp1 = "temp1.out";

parse_command_line(@ARGV);

compare_executables("data/corner3D.B.0.5x.nrrd", "data/corner3D.B.0.5x.grad.nrrd", "-listg -cube 86 -isovalue 3.2");
compare_executables("data/corner3D.B.0.5x.nrrd", "data/corner3D.B.0.5x.grad.nrrd", "-listg -cube 86 -gradCS -isovalue 3.2");
compare_executables("data/corner3D.B.0.5x.nrrd", "data/corner3D.B.0.5x.grad.nrrd", "-listg -cube 86 -gradNS -isovalue 3.2");
compare_executables("data/corner3D.B.0.5x.nrrd", "data/corner3D.B.0.5x.grad.nrrd", "-listg -cube 86 -gradES -isovalue 3.2");
compare_executables("data/corner3D.B.0.5x.nrrd", "data/corner3D.B.0.5x.grad.nrrd", "-listg -cube 86 -gradEC -isovalue 3.2");

# *********************************************

sub usage_error {
  print "Usage: testsharpinfo.pl {prog1} {prog2}\n";
  exit(10);
}

sub parse_command_line {

  my @command_line = @_;

  if (scalar(@_) != 2) { usage_error(); }

  $prog0 = $command_line[0];
  $prog1 = $command_line[1];
}

sub compare_executables {

  (scalar(@_) == 3) ||
    die "Error in compare_executables.  Requires 3 parameters.\n";

  my $scalar_input = shift;
  my $gradient_input = shift;
  my $options = shift;

  run_sharpinfo($prog0, $scalar_input, $gradient_input, $temp0, $options);
  run_sharpinfo($prog1, $scalar_input, $gradient_input, $temp1, $options);
  diff_files($temp0, $temp1);
}

sub run_sharpinfo {

  (scalar(@_) == 5) ||
    die "Error in run_sharpinfo.  Requires 5 parameters.\n";

  my $prog = shift;
  my $input_scalar = shift;
  my $input_gradient = shift;
  my $output_filename = shift;
  my $options = shift;

  my $command_line =
    "./$prog $options  $input_scalar $input_gradient > $output_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program ijksharpinfo abnormally terminated.\n";
}

sub diff_files {

  scalar(@_) == 2 ||
    die "Error in sub run_ijkdual. Requires exactly 2 parameters.\n";

  my $flag = system("diff --brief $_[0] $_[1] > /dev/null");
  if ($flag != 0) { 
    print "*** Output differs. ***\n\n";
    $found_difference = 1; 
  };

}
