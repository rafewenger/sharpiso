#!/usr/bin/perl
# test and compare different versions of mergesharp
# updated by R. Wenger April 16, 2013

use strict;

my $testdir = "data";

my @proglist = @ARGV;
my @input_options;
my @input_options0;
my $fastflag = 0;
my $veryfastflag = 0;
my $alloptions = 0;
my $found_difference = 0;
my $isoval_offset = 0;
my $outfile = "temp.off";
my $outfile0 = "temp0.off";
my $min_diff_flag = 0;
my $min_diff = 0.0;
my $flag_info = 0;
my $flag_sort = 0;
my $random_flag = 0;
my $random_seed = 0;
my $bzero_flag = 0;
my $axis_size = 20;      # axis size for random test

my %data_flag;
my $use_all_data = 0;


# mergesharp arguments which take an input value/string.
my @mergesharp_options = ( "-subsample",  "-position", "-round", 
                           "-max_dist", "-max_grad_dist",
                           "-max_eigen", "-merge_linf_th",
                           "-min_triangle_angle", "-min_normal_angle");

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

  if ($new_option eq "-hermite") {
    $data_flag{hermite} = 1;
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

  if ($new_option eq "-sort") { 
    $flag_sort = 1; 
    next;
  };

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

  if ($new_option eq "-info") {
    $flag_info = 1;
  }

  if ($new_option eq "-random") {
    $random_flag = 1;
    $random_seed = shift(@proglist);
    next;
  }

  if ($new_option eq "-axis_size") {
    $axis_size = shift(@proglist);
    next;
  }

  if ($new_option eq "-bzero") {
    $bzero_flag = 1;
    next;
  }

  push(@input_options, $new_option);
  push(@input_options0, $new_option);

  if (scalar(@proglist) > 0) {
    if (is_mergesharp_option("$new_option")) {
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

  if (!$veryfastflag) {
    $testdata{cube_A}{fname} = "cube.A.nrrd";
    $testdata{cube_A}{isovalue} = [ 4.9, 5, 5.5 ];

    $testdata{cube_B}{fname} = "cube.B.nrrd";
    $testdata{cube_B}{isovalue} = [ 5, 5.5];
  }
}



if ($use_all_data || defined($data_flag{cubes})) {

  if (!$veryfastflag) {
    $testdata{cube_dir111_A}{fname} = "cube.dir111.A.nrrd";
    $testdata{cube_dir111_A}{isovalue} = [ 4.9, 5, 5.5 ];

    $testdata{cube_dir111_B}{fname} = "cube.dir111.B.nrrd";
    $testdata{cube_dir111_B}{isovalue} = [ 4.9, 5, 5.5 ];

    $testdata{cube_dir321_A}{fname} = "cube.dir111.A.nrrd";
    $testdata{cube_dir321_A}{isovalue} = [ 4.9, 5, 5.5 ];

    $testdata{cube_dir321_A}{fname} = "cube.dir111.B.nrrd";
    $testdata{cube_dir321_A}{isovalue} = [ 4.9, 5, 5.5 ];

    $testdata{cube_dir321_A}{fname} = "cube.dir321.A.nrrd";
    $testdata{cube_dir321_A}{isovalue} = [ 4.9, 5, 5.5 ];
  }
  else {
    $testdata{cube_dir321_A}{fname} = "cube.dir321.A.nrrd";
    $testdata{cube_dir321_A}{isovalue} = [ 4.9, 5 ];
  }


}


if ($use_all_data || defined($data_flag{twocubes})) {

  if (!$veryfastflag) {
    $testdata{twocubes_A}{fname} = "twocubes.A.nrrd";
    $testdata{twocubes_A}{isovalue} = [ 9.8, 10, 10.1 ];

    $testdata{twocubes_dir111_A}{fname} = "twocubes.dir111.A.nrrd";
    $testdata{twocubes_dir111_A}{isovalue} = [ 9.8, 10, 10.1 ];
  }
  else {
    $testdata{twocubes_dir111_A}{fname} = "twocubes.dir111.A.nrrd";
    $testdata{twocubes_dir111_A}{isovalue} = [ 9.8, 10 ];
  }

}

if ($use_all_data || defined($data_flag{annulus})) {

  if (!$veryfastflag) {

    $testdata{annulus_dir100_A}{fname} = "annulus.dir100.A.nrrd";
    $testdata{annulus_dir100_A}{isovalue} = [ 4.9, 5 ];
  }

  $testdata{annulus_dir100_A}{fname} = "annulus.dir111.A.nrrd";
  $testdata{annulus_dir100_A}{isovalue} = [ 4.9, 5 ];
}

if ($use_all_data || defined($data_flag{flange})) {

  if (!$veryfastflag) {
    $testdata{flange_dir100_A}{fname} = "flange.dir100.A.nrrd";
    $testdata{flange_dir100_A}{isovalue} = [ 4.9, 5 ];
  }

  $testdata{flange_dir111_A}{fname} = "flange.dir111.A.nrrd";
  $testdata{flange_dir111_A}{isovalue} = [ 4.9, 5 ];
}

if (defined($data_flag{hermite})) {

    $testdata{weld}{fname} = "weld.nhdr";
    $testdata{weld}{normalFile} = "weld-normals.off";
    $testdata{weld}{isovalue} = [ 0.5 ];
}


# test data parameters for -random
my %random_test;

$random_test{seed} = $random_seed;
$random_test{isovalue} = 4;
$random_test{axis_size} = $axis_size;
$random_test{num_random} = 10;
$random_test{nrrd_filename} = "rtest.nrrd";


my $prog0 = shift(@proglist);

if (scalar(@proglist) == 0) {

  if ($random_flag) {
    run_mergesharp_on_random(@input_options);
  }
  else {
    run_mergesharp_count_sharp_edges(@input_options);
  }
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
 print "Usage: testmergesharp.pl [OPTIONS] {prog1} {prog2} [{prog3} ...]";
 exit(10);
}

# run mergesharp and count number of sharp edges
sub run_mergesharp_count_sharp_edges {

  my @option_list = ("-trimesh", "@_");

  foreach my $tdata (keys %testdata) {

    my $tfile = "$testdir"."/"."$testdata{$tdata}{fname}";
    my @isovalue_list = @{$testdata{$tdata}{isovalue}};
    my @option_listB = @option_list;
    if (exists $testdata{$tdata}{normalFile}) {
      my $nfile = "$testdir"."/"."$testdata{$tdata}{normalFile}";
      push(@option_listB, "-normal", "$nfile");
    };

    foreach my $isoval (@isovalue_list) {

      $isoval = $isoval + $isoval_offset;
      run_mergesharp($prog0, $tfile, "$outfile0", \@option_listB, $isoval);

      count_sharp_edges("$outfile0");
    }
  }

}

# run mergesharp on random sharp data
# count sharp edges and check manifold conditions.
sub run_mergesharp_on_random {

  my @option_list = ("-trimesh", "@_");

  my $num_random = $random_test{num_random};
  my $isoval = $random_test{isovalue} + $isoval_offset;
  my $nrrd_filename = $random_test{nrrd_filename};
  my $axis_size = $random_test{axis_size};
  my $randompos_seed = $random_test{seed};
  my $randomdir_seed = $randompos_seed+10000;

  for (my $i = 0; $i < $num_random; $i++) {

    my $command_line =
      "ijkgenscalar -grad -dim 3 -asize $axis_size -randompos $randompos_seed -randomdir $randomdir_seed -n 5 -field cube -dir \"1 1 1\" -side_dir \"1 0 0\" -s";
    if ($bzero_flag) {
      $command_line = $command_line . " -bzero";
    }
    $command_line = $command_line . " $nrrd_filename";
    print "$command_line\n";

    system("$command_line") == 0 ||
      die "Program ijkgenscalar abnormally terminated.\n";

    run_mergesharp
      ($prog0, $nrrd_filename, "$outfile0", \@option_list, $isoval);
    count_sharp_edges("$outfile0");

    $randompos_seed++;
    $randomdir_seed++;
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
  compare_executables("-position gradNS -interpolate_edgeI -dist2centroid @input_options");
  print "\n";
  compare_executables("-position gradNS -sharp_edgeI -dist2centroid @input_options");
  print "\n";
  compare_executables("-position gradNS -dist2center @input_options");
  print "\n";
  compare_executables("-position gradCD -sep_pos @input_options");
  print "\n";
  compare_executables("-position gradCD -sep_neg @input_options");
  print "\n";
  compare_executables("-position gradCD -no_merge_sharp @input_options");
  print "\n";
  compare_executables("-position gradCD -no_merge_sharp -resolve_ambig @input_options");
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
    my @option_listB0 = @option_list0;
    my @option_listB = @option_list;
    if (exists $testdata{$tdata}{normalFile}) {
      my $nfile = "$testdir"."/"."$testdata{$tdata}{normalFile}";
      push(@option_listB0, "-normal", "$nfile");
      push(@option_listB, "-normal", "$nfile");
    };

    foreach my $isoval (@isovalue_list) {

      $isoval = $isoval + $isoval_offset;
      run_mergesharp($prog0, $tfile, "$outfile0", \@option_listB0, $isoval);
      if ($flag_sort) { system("ijkmesh sort $outfile0"); }
      foreach my $mergesharp (@proglist) {
        run_mergesharp($mergesharp, $tfile, "$outfile", \@option_listB, $isoval);
        if ($flag_sort) { system("ijkmesh sort $outfile"); }
        diff_files("$outfile0", "$outfile");
      }
    }
  }
}

sub run_mergesharp {

  scalar(@_) > 4 ||
    die "Error in sub run_mergesharp. Requires at least 5 parameters.\n";

  my $mergesharp = $_[0];
  my $input_filename = $_[1];
  my $output_filename = $_[2];
  my @option_list = @{$_[3]};
  my $isoval1 = $_[4];

  my $command_line;
  $command_line = "./$mergesharp @option_list";
  if (!$flag_info) 
    { $command_line = "$command_line" . " -s"; }
  $command_line = 
    "$command_line" . " -o $output_filename $isoval1 $input_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program mergesharp abnormally terminated.\n";

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

  my $command_line = "findsharp 140 $output_filename";
  # print "$command_line\n";
  system("$command_line >& /dev/null") == 0 ||
    die "Program findsharp abnormally terminated.\n";

  $command_line = "countdegree $line_filename > $count_filename";
  # print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program countdegree abnormally terminated.\n";

  system("fgrep \"degree 1 or 3\" $count_filename");

  $command_line = "ijkmeshinfo -report_deep -terse -manifold $output_filename";
  my $return_val = system("$command_line");
  $return_val = $return_val >> 8;
  ($return_val == 0 || $return_val == 1) ||
    die "Program ijkmeshinfo abnormally terminated.\n";
}

sub diff_files {

  scalar(@_) == 2 ||
    die "Error in sub run_mergesharp. Requires exactly 2 parameters.\n";

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

# return true (1) if $_[0] is an mergesharp option which takes an 
#   input value/string
sub is_mergesharp_option {

  scalar(@_) == 1 ||
    die "Error in sub is_mergesharp_option. Requires exactly 1 parameter.\n";

  foreach my $option (@mergesharp_options) {
    if ($_[0] eq "$option") {
      return 1;
    }
  };

  return 0;
}


