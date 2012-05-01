#!/usr/bin/perl
# run isodual3D with different edge intersection calculations
# updated by R. Wenger Apr. 12, 2012

use strict;

my $flag_reposition = 0;
my @input_list = @ARGV;
my $allpositions = 0;
my @count_files;

# isodual3D arguments which take an input value/string.
my @isodual3D_options = ( "-subsample",  "-position", "-round" );

my @input_options;
while (scalar(@input_list) > 0 &&
       $input_list[0] =~ /-.*/) {

  my $new_option = shift(@input_list);

  if ("$new_option" eq "-reposition") 
    { $flag_reposition = 1; }

  if ($new_option eq "-allpos" || $new_option eq "-allpositions") {
    $allpositions = 1;
    next;
  }

  push(@input_options, $new_option);

  if (scalar(@input_list) > 0) {
    if (is_isodual3D_option("$new_option")) {
      $new_option = shift(@input_list);
      push(@input_options, $new_option);
    }
  }

}

if (scalar(@input_list) != 2) { usage_error(); }

my $isovalue = $input_list[0];
my $infile = $input_list[1];

if ($allpositions) {
  run_isodual3D_position($isovalue, $infile, "centroid");
  run_isodual3D_position($isovalue, $infile, "gradES");
  run_isodual3D_position($isovalue, $infile, "gradEC");
  run_isodual3D_position($isovalue, $infile, "gradC");
  run_isodual3D_position($isovalue, $infile, "gradIE");
  run_isodual3D_position($isovalue, $infile, "gradIES");
  run_isodual3D_position($isovalue, $infile, "gradCD");
  run_isodual3D_position($isovalue, $infile, "gradCDdup");
  run_isodual3D_position($isovalue, $infile, "gradN");
  run_isodual3D_position($isovalue, $infile, "gradNS");
  run_isodual3D_position($isovalue, $infile, "gradNIE");
  run_isodual3D_position($isovalue, $infile, "gradNIES");

  print "\n\n";
}
else {
  run_isodual3D($isovalue, $infile, "out.off", \@input_options);
}

output_count();

# ***********************************

sub usage_error {
 print "Usage: runisodual3D.pl {isovalue} {scalar input file}]";
 exit(10);
}

sub run_isodual3D_position {

  scalar(@_) == 3 ||
    die "Error in sub runisodual. Requires 3 parameters.\n";

  my $isovalue = $_[0];
  my $input_filename = $_[1];
  my $position_option = $_[2];

  my $prefix = "$position_option";
  if ($flag_reposition)
    { $prefix = "$prefix" . ".R"; }

  my $output_filename = "$prefix" . ".off";
  my @option_list = ("-position $position_option", @input_options);
  run_isodual3D($isovalue, $input_filename, $output_filename, \@option_list);
}

sub run_isodual3D {

  scalar(@_) == 4 ||
    die "Error in sub run_isodual. Requires 4 parameters.\n";

  my $isovalue = $_[0];
  my $input_filename = $_[1];
  my $output_filename = $_[2];
  my @option_list = @{$_[3]};

  my $line_filename = "$output_filename";
  $line_filename =~ s/.off/.line/;
  my $count_filename = "$output_filename";
  $count_filename =~ s/.off/.count/;

  my $command_line =
      "isodual3D @option_list -trimesh -s -o $output_filename $isovalue $input_filename";

  if (-e "./isodual3D") {
    # Use version in current directory
    $command_line = "./" . "$command_line";
  }
  # else use version in bin

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program isodual3D abnormally terminated.\n";
  $command_line = "findedge 140 $output_filename >& /dev/null";
  # print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program findedge abnormally terminated.\n";

  $command_line = "findEdgeCount $line_filename > $count_filename";
  # print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program findedge abnormally terminated.\n";

  push(@count_files, $count_filename);
}

sub output_count {

  foreach my $count_filename (@count_files) {
    system("fgrep -H \"degree 1 or 3\" $count_filename");
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

