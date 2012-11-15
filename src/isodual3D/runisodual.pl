#!/usr/bin/perl
# run isodual3D with different edge intersection calculations
# updated by R. Wenger Apr. 12, 2012

use strict;
use File::Basename;

my $flag_reposition = 0;
my @input_list = @ARGV;
my $allpositions = 0;
my @count_files;
my $annulus_flag = 0;
my $flange_flag = 0;
my $twocubes_flag = 0;
my $predefined_data = 0;
my $data_dir = "data";
my $print_command_line = 1;
my $print_count_filename = 1;

# isodual3D arguments which take an input value/string.
my @isodual3D_options = ( "-subsample",  "-position", "-pos", "-gradient",
                          "-round", "-max_dist", "-snap_dist" );

my @input_options=();
while (scalar(@input_list) > 0 &&
       $input_list[0] =~ /-.*/) {

  my $new_option = shift(@input_list);

  if ("$new_option" eq "-reposition") 
    { $flag_reposition = 1; }

  if ($new_option eq "-allpos" || $new_option eq "-allpositions") {
    $allpositions = 1;
    next;
  }

  if ("$new_option" eq "-annulus") {
    $annulus_flag = 1;
    next;
  }

  if ("$new_option" eq "-flange") {
    $flange_flag = 1;
    next;
  }

  if ("$new_option" eq "-twocubes") {
    $twocubes_flag = 1;
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

if ($annulus_flag || $twocubes_flag || $flange_flag) {
  $predefined_data = 1;
}

if ($predefined_data) {
  $print_command_line = 0;
  $print_count_filename = 0;

  if ($annulus_flag) 
    { run_isodual3D_on_annulus(\@input_options); }

  if ($flange_flag) 
    { run_isodual3D_on_flange(\@input_options); }

  if ($twocubes_flag) 
    { run_isodual3D_on_twocubes(\@input_options); }
}
else {

  if (scalar(@input_list) != 2) { usage_error(); }

  my $isovalue = $input_list[0];
  my $infile = $input_list[1];

  if ($allpositions) {
    run_isodual3D_all_positions($isovalue, $infile);
    output_count();
  }
  else
    {
      $print_command_line = 1;
      $print_count_filename = 0;
      compute_and_output_count
        ($isovalue, $infile, "out.off", \@input_options);
    }
}


# ***********************************
# run_isodual3D SUBROUTINES
# ***********************************

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
      "isodual3D @option_list -trimesh -o $output_filename $isovalue $input_filename";

  if (-e "./isodual3D") {
    # Use version in current directory
    $command_line = "./" . "$command_line";
  }
  # else use version in bin

  if ($print_command_line) 
    { print "$command_line\n"; }

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

sub compute_and_output_count {

  scalar(@_) == 4 ||
    die "Error in sub compute_and_output_count. Requires 4 parameters.\n";

  my $isovalue = $_[0];
  my $input_pathname = $_[1];
  my $output_filename = $_[2];
  my @option_list = @{$_[3]};

  my $input_filename = fileparse("$input_pathname");
   
  run_isodual3D(@_);
  print "$input_filename:$isovalue ";
  output_count();
  @count_files = ();
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

sub run_isodual3D_all_positions {

  scalar(@_) == 2 ||
    die "Error in sub run_isodual_all_positions. Requires 2 parameters.\n";

  my $isovalue = $_[0];
  my $infile = $_[1];

  run_isodual3D_position($isovalue, $infile, "centroid");
  run_isodual3D_position($isovalue, $infile, "gradES");
  run_isodual3D_position($isovalue, $infile, "gradEC");
  run_isodual3D_position($isovalue, $infile, "gradCD");
  run_isodual3D_position($isovalue, $infile, "gradCDdup");
  run_isodual3D_position($isovalue, $infile, "gradC");
  run_isodual3D_position($isovalue, $infile, "gradIE");
  run_isodual3D_position($isovalue, $infile, "gradIEDir");
  run_isodual3D_position($isovalue, $infile, "gradIES");
  run_isodual3D_position($isovalue, $infile, "gradN");
  run_isodual3D_position($isovalue, $infile, "gradNS");
  run_isodual3D_position($isovalue, $infile, "gradNIE");
  run_isodual3D_position($isovalue, $infile, "gradNIES");

  print "\n\n";
}

sub run_isodual3D_on_annulus {

  scalar(@_) == 1 ||
    die "Error in sub run_isodual3D_on_annulus. Requires 1 parameter.\n";

  my @option_list = @{$_[0]};

  my @annulus_files = 
    ( "annulus3D.A31x.nrrd", "annulus3D.B31x.nrrd", "annulus3D.C31x.nrrd",
      "annulus3D.D31x.nrrd", "annulus3D.E31x.nrrd", "annulus3D.F31x.nrrd");

  my @isovalue_list = ( 4.0, 4.01, 4.1, 4.2, 4.3, 4.4, 4.5 );

  foreach my $filename (@annulus_files) {
    foreach my $isovalue (@isovalue_list) {
      compute_and_output_count
        ($isovalue, "$data_dir/$filename", "out.off", \@option_list);
    }
  }

}

sub run_isodual3D_on_flange {

  scalar(@_) == 1 ||
    die "Error in sub run_isodual3D_on_flange. Requires 1 parameter.\n";

  my @option_list = @{$_[0]};

  my @flange_files = 
    ( "flange3D.A61x.nrrd", "flange3D.B61x.nrrd", "flange3D.C61x.nrrd",
      "flange3D.D61x.nrrd", "flange3D.E61x.nrrd", "flange3D.F61x.nrrd");

  my @isovalue_list = ( 3.0, 3.01, 3.1, 3.2, 3.3, 3.4, 3.5 );

  foreach my $filename (@flange_files) {
    foreach my $isovalue (@isovalue_list) {
      compute_and_output_count
        ($isovalue, "$data_dir/$filename", "out.off", \@option_list);
    }
  }

}

sub run_isodual3D_on_twocubes {

  scalar(@_) == 1 ||
    die "Error in sub run_isodual3D_on_twocubes. Requires 1 parameter.\n";

  my @option_list = @{$_[0]};

  my @twocubes_files = 
    ( "twocubes3D.A21x.nrrd", "twocubes3D.B21x.nrrd", "twocubes3D.C21x.nrrd" );

  my @isovalue_list = ( 4.0, 4.01, 4.1, 4.2, 4.3, 4.4, 4.5 );

  foreach my $filename (@twocubes_files) {
    foreach my $isovalue (@isovalue_list) {
      compute_and_output_count
        ($isovalue, "$data_dir/$filename", "out.off", \@option_list);
    }
  }

}


# ***********************************
# MISC SUBROUTINES
# ***********************************

sub output_count {

  foreach my $count_filename (@count_files) {

    if ($print_count_filename) 
      { system("fgrep -H \"degree 1 or 3\" $count_filename"); }
    else
      { system("fgrep \"degree 1 or 3\" $count_filename"); }
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

sub usage_error {
 print "Usage: runisodual3D.pl [OPTIONS] {isovalue} {scalar input file}]\n";
 print "OPTIONS: [-allpos] [isodual3D options]\n";
 exit(10);
}

