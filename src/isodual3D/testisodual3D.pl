#!/usr/bin/perl
# test and compare different versions of isodual3D
# updated by R. Wenger Feb. 25, 2012

use strict;

my $testdir = "data";

my @proglist = @ARGV;
my @input_options;
my @input_options0;
my $fastflag = 0;
my $veryfastflag = 0;
my $veryveryfastflag = 0;
my $alloptions = 0;
my $found_difference = 0;
my $isoval_offset = 0;
my $outfile = "temp.off";
my $outfile0 = "temp0.off";
my $min_diff_flag = 0;
my $min_diff = 0.0;

my %data_flag;
my $use_all_data = 0;

# isodual3D arguments which take an input value/string.
my @isodual3D_options = ( "-subsample",  "-position", "-round", 
                          "-max_dist", "-max_eigen" );

while (scalar(@proglist) > 0 &&
       $proglist[0] =~ /^\-.*/) {

  my $new_option = shift(@proglist);

  if ($new_option eq "-fast") { 
    $fastflag = 1; 
    next;
  };

  if ($new_option eq "-veryfast") { 
    $veryfastflag = 1; 
    next;
  };

  if ($new_option eq "-veryveryfast") { 
    $veryveryfastflag = 1; 
    $veryfastflag = 1;
    next;
  };

  if ($new_option eq "-flange") {
    $data_flag{flange} = 1;
    next;
  }

  if ($new_option eq "-annulus") {
    $data_flag{annulus} = 1;
    next;
  }

  if ($new_option eq "-cubes") {
    $data_flag{cubes} = 1;
    next;
  }

  if ($new_option eq "-twocubes") {
    $data_flag{twocubes} = 1;
    next;
  }

  if ($new_option eq "-offset") {
    $isoval_offset = shift(@proglist);
    next;
  }

  if ($new_option eq "-alloptions") {
    $alloptions = 1;
    next;
  }

  if ($new_option eq "-min_diff") {
    $min_diff_flag = 1;
    $min_diff = shift(@proglist);
    next;
  }

  if ($new_option eq "-max_eigen0") {
    push(@input_options0, "-max_eigen");
    $new_option = shift(@proglist);
    push(@input_options0, $new_option);
    next;
  }

  if ($new_option eq "-lindstrom_fast0") {
    push(@input_options0, "-lindstrom_fast");
    next;
  }

  if ($new_option eq "-single_isov0") {
    push(@input_options0, "-single_isov");
    next;
  }

  if ($new_option eq "-no_check_disk0") {
    push(@input_options0, "-no_check_disk");
    next;
  }

  push(@input_options, $new_option);
  push(@input_options0, $new_option);

  if (scalar(@proglist) > 0) {
    if (is_isodual3D_option("$new_option")) {
      $new_option = shift(@proglist);
      push(@input_options, $new_option);
      push(@input_options0, $new_option);
    }
  }
}

if (scalar(@proglist) < 1) { usage_error(); };

my %testdata;

if ((scalar keys %data_flag) == 0) {
  # use all data
  $use_all_data = 1;
}


if ($use_all_data) {

  if (!$veryveryfastflag) {
    $testdata{cube_A20}{fname} = "cube3D.A20x.nrrd";
    $testdata{cube_A20}{isovalue} = [ 4.9, 5, 5.5 ];

    $testdata{cube_B20}{fname} = "cube3D.B20x.nrrd";
    $testdata{cube_B20}{isovalue} = [ 5, 5.5];


  }
}



if ($use_all_data || defined($data_flag{cubes})) {

  if (!$veryfastflag) {

    if (!$veryveryfastflag) {

      $testdata{tcube110A}{fname} = "tcube110A.31x.nrrd";
      $testdata{tcube110A}{isovalue} = [ 4.9, 5, 5.5 ];

      $testdata{tcube111A}{fname} = "tcube111A.31x.nrrd";
      $testdata{tcube111A}{isovalue} = [ 4.9, 5, 5.5 ];

      $testdata{tcube211A}{fname} = "tcube211A.31x.nrrd";
      $testdata{tcube211A}{isovalue} = [ 4.9, 5, 5.5 ];
    }
  }
}


if ($use_all_data || defined($data_flag{twocubes})) {

  if (!$veryfastflag) {
    $testdata{twocubes_A21}{fname} = "twocubes3D.A21x.nrrd";
    $testdata{twocubes_A21}{isovalue} = [ 4, 4.5 ];
  }

  $testdata{twocubes_B21}{fname} = "twocubes3D.B21x.nrrd";
  $testdata{twocubes_B21}{isovalue} = [ 4, 4.5 ];
}

if ($use_all_data || defined($data_flag{annulus})) {

  if (!$veryveryfastflag) {

    if (!$veryfastflag) {
      $testdata{annulus_A31}{fname} = "annulus3D.A31x.nrrd";
      $testdata{annulus_A31}{isovalue} = [ 4, 4.5 ];

      $testdata{annulus_B31}{fname} = "annulus3D.B31x.nrrd";
      $testdata{annulus_B31}{isovalue} = [ 4, 4.1, 4.5 ];

      $testdata{annulus_C31}{fname} = "annulus3D.C31x.nrrd";
      $testdata{annulus_C31}{isovalue} = [ 4, 4.5 ];

      $testdata{annulus_D31}{fname} = "annulus3D.D31x.nrrd";
      $testdata{annulus_D31}{isovalue} = [ 4, 4.5 ];

      $testdata{annulus_E31}{fname} = "annulus3D.E31x.nrrd";
      $testdata{annulus_E31}{isovalue} = [ 4, 4.5 ];
    }

    $testdata{annulus_F31}{fname} = "annulus3D.F31x.nrrd";
    $testdata{annulus_F31}{isovalue} = [ 4, 4.5 ];
  }
}

if ($use_all_data || defined($data_flag{flange})) {

  if (!$veryfastflag) {
    $testdata{flange_A61}{fname} = "flange3D.A61x.nrrd";
    $testdata{flange_A61}{isovalue} = [ 4, 4.5 ];

    $testdata{flange_B61}{fname} = "flange3D.B61x.nrrd";
    $testdata{flange_B61}{isovalue} = [ 4, 4.5 ];

    $testdata{flange_C61}{fname} = "flange3D.C61x.nrrd";
    $testdata{flange_C61}{isovalue} = [ 4, 4.5 ];

    $testdata{flange_D61}{fname} = "flange3D.D61x.nrrd";
    $testdata{flange_D61}{isovalue} = [ 4, 4.5 ];

    $testdata{flange_E61}{fname} = "flange3D.E61x.nrrd";
    $testdata{flange_E61}{isovalue} = [ 4, 4.5 ];
  }

  $testdata{flange_F61}{fname} = "flange3D.F61x.nrrd";
  $testdata{flange_F61}{isovalue} = [ 4, 4.5, 8 ];
}


my $prog0 = shift(@proglist);

if (scalar(@proglist) == 0) {

  run_isodual3D_count_sharp_edges(@input_options);
}
else {

  if ($alloptions) {
    compare_executables_all_options();
  } 
  else {
    compare_executables("@input_options0", "@input_options");
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

# run isodual3D and count number of sharp edges
sub run_isodual3D_count_sharp_edges {

  my @option_list = ("-trimesh", "@_");

  foreach my $tdata (keys %testdata) {

    my $tfile = "$testdir"."/"."$testdata{$tdata}{fname}";
    my @isovalue_list = @{$testdata{$tdata}{isovalue}};
    foreach my $isoval (@isovalue_list) {

      $isoval = $isoval + $isoval_offset;
      run_isodual3D($prog0, $tfile, "$outfile0", \@option_list, $isoval);

      count_sharp_edges("$outfile0");
    }
  }

}

# compare executables with all possible options
sub compare_executables_all_options {

  compare_executables("-position centroid @input_options");
  print "\n";
  compare_executables("-position cube_center @input_options");
  print "\n";
  compare_executables("-position gradC @input_options");
  print "\n";
  compare_executables("-position gradCS @input_options");
  print "\n";
  compare_executables("-position gradN @input_options");
  print "\n";
  compare_executables("-position gradNS @input_options");
  print "\n";
  compare_executables("-position gradES @input_options");
  print "\n";
  compare_executables("-position gradEC @input_options");
  print "\n";
  compare_executables("-position gradIE @input_options");
  print "\n";
  compare_executables("-position gradIEDir @input_options");
  print "\n";
  compare_executables("-position gradCD @input_options");
  print "\n";
  compare_executables("-position gradCDdup @input_options");
  print "\n";
  compare_executables("-position gradCD -sep_pos @input_options");
  print "\n";
  compare_executables("-position gradCD -sep_neg @input_options");
  print "\n";
  compare_executables("-position gradCD -resolve_ambig @input_options");
  print "\n";
  compare_executables("-position centroid -trimesh @input_options");
  print "\n";
  compare_executables("-position gradEC -trimesh @input_options");
  print "\n";
  compare_executables("-position gradNS -trimesh @input_options");
  print "\n";
  compare_executables("-position centroid -uniform_trimesh @input_options");
  print "\n";
  compare_executables("-position gradEC -uniform_trimesh @input_options");
  print "\n";
}

sub compare_executables {

  my @option_list0 = @_[0];
  my @option_list;

  if (defined @_[1]) 
    { @option_list = @_[1]; }
  else
    { @option_list = @option_list0; }

  foreach my $tdata (keys %testdata) {

    my $tfile = "$testdir"."/"."$testdata{$tdata}{fname}";
    my @isovalue_list = @{$testdata{$tdata}{isovalue}};
    foreach my $isoval (@isovalue_list) {

      $isoval = $isoval + $isoval_offset;
      run_isodual3D($prog0, $tfile, "$outfile0", \@option_list0, $isoval);
      foreach my $isodual (@proglist) {
        run_isodual3D($isodual, $tfile, "$outfile", \@option_list, $isoval);
        diff_files("$outfile0", "$outfile");
      }
    }
  }
}

sub run_isodual3D {

  scalar(@_) > 4 ||
    die "Error in sub run_isodual3D. Requires at least 5 parameters.\n";

  my $isodual = $_[0];
  my $input_filename = $_[1];
  my $output_filename = $_[2];
  my @option_list = @{$_[3]};
  my $isoval1 = $_[4];

  my $command_line;
  $command_line = 
    "./$isodual @option_list -s -o $output_filename $isoval1 $input_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program isodual3D abnormally terminated.\n";

}

# count number of sharp edges
sub count_sharp_edges {

  scalar(@_) != 0 ||
    die "Error in sub count_sharp_edges. Requires one parameter.\n";

  my $output_filename = $_[0];

  my $line_filename = $output_filename;
  $line_filename =~ s/.off/.line/;
  my $count_filename = $output_filename;
  $count_filename =~ s/.off/.count/;

  my $command_line = "findedge 140 $output_filename";
  # print "$command_line\n";
  system("$command_line >& /dev/null") == 0 ||
    die "Program findedge abnormally terminated.\n";

  $command_line = "findEdgeCount $line_filename > $count_filename";
  # print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program findedge abnormally terminated.\n";

  system("fgrep \"degree 1 or 3\" $count_filename");

  $command_line = "ijkmeshinfo -terse -manifold $output_filename";
  my $return_val = system("$command_line");
  $return_val = $return_val >> 8;
  ($return_val == 0 || $return_val == 1) ||
    die "Program ijkmeshinfo abnormally terminated.\n";
}

sub diff_files {

  scalar(@_) == 2 ||
    die "Error in sub run_isodual3D. Requires exactly 2 parameters.\n";

  my $flag = system("diff --brief $_[0] $_[1] > /dev/null");
  if ($flag != 0) { 
    print "*** Output .off files differ. ***\n";
    $found_difference = 1; 

    if ($min_diff_flag) 
      { system("ijkmeshdiff -terse -min_diff $min_diff $_[0] $_[1]"); }
    print "\n";
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


