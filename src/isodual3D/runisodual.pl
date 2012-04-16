#!/usr/bin/perl
# run isodual3D with different edge intersection calculations
# updated by R. Wenger Apr. 12, 2012

use strict;

if (scalar(@ARGV) != 2) { usage_error(); }

my $isovalue = $ARGV[0];
my $infile = $ARGV[1];

run_isodual3D($isovalue, $infile, "gradES");
run_isodual3D($isovalue, $infile, "gradEC");
run_isodual3D($isovalue, $infile, "gradIE");
run_isodual3D($isovalue, $infile, "gradNIE");
run_isodual3D($isovalue, $infile, "gradNS");

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

  my $output_filename = "$position_option" . ".off";
  my $line_filename = "$position_option" . ".line";
  my $count_filename = "$position_option" . ".count";

  my $command_line;
  $command_line = 
    "./isodual3D -trimesh -position $position_option -s -o $output_filename $isovalue $input_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program isodual3D abnormally terminated.\n";

  $command_line = "findedge 140 $output_filename";
  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program findedge abnormally terminated.\n";

  $command_line = "findEdgeCount $line_filename > $count_filename";
  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program findedge abnormally terminated.\n";
}



