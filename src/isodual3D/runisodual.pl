#!/usr/bin/perl
# run isodual3D with different edge intersection calculations
# updated by R. Wenger Apr. 12, 2012

use strict;

my $flag_reposition = 0;
my @input_list = @ARGV;

my @isodual3D_options;
while (scalar(@input_list) > 0 &&
       $input_list[0] =~ /-.*/) {

  my $new_option = shift(@input_list);

  if ("$new_option" eq "-reposition") 
    { $flag_reposition = 1; }

  push(@isodual3D_options, $new_option);
}

if (scalar(@input_list) != 2) { usage_error(); }

my $isovalue = $input_list[0];
my $infile = $input_list[1];

run_isodual3D($isovalue, $infile, "centroid");
run_isodual3D($isovalue, $infile, "gradES");
run_isodual3D($isovalue, $infile, "gradEC");
run_isodual3D($isovalue, $infile, "gradC");
run_isodual3D($isovalue, $infile, "gradIE");
run_isodual3D($isovalue, $infile, "gradIES");
run_isodual3D($isovalue, $infile, "gradCD");
run_isodual3D($isovalue, $infile, "gradN");
run_isodual3D($isovalue, $infile, "gradNS");
run_isodual3D($isovalue, $infile, "gradNIE");
run_isodual3D($isovalue, $infile, "gradNIES");

print "\n\n";
system("fgrep \"degree 1 or 3\" \*.count");

# ***********************************

sub usage_error {
 print "Usage: runisodual3D.pl {isovalue} {scalar input file}]";
 exit(10);
}

sub run_isodual3D {

  scalar(@_) == 3 ||
    die "Error in sub runisodual. Requires 4 parameters.\n";

  my $isovalue = $_[0];
  my $input_filename = $_[1];
  my $position_option = $_[2];

  my $prefix = "$position_option";
  if ($flag_reposition)
    { $prefix = "$prefix" . ".R"; }

  my $output_filename = "$prefix" . ".off";
  my $line_filename = "$prefix" . ".line";
  my $count_filename = "$prefix" . ".count";

  my $command_line;
  if (-e "./isodual3D") {
    # Use version in current directory
    $command_line = 
      "./isodual3D @isodual3D_options -trimesh -position $position_option -s -o $output_filename $isovalue $input_filename";
  }
  else {
    # Use version in bin
    $command_line = 
      "isodual3D @isodual3D_options -trimesh -position $position_option -s -o $output_filename $isovalue $input_filename";
  }

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program isodual3D abnormally terminated.\n";
  $command_line = "findedge 140 $output_filename";
  # print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program findedge abnormally terminated.\n";

  $command_line = "findEdgeCount $line_filename > $count_filename";
  # print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program findedge abnormally terminated.\n";
}



