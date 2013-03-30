#!/usr/bin/perl
# test and compare different versions of grad2hermite
# updated by R. Wenger March 12, 2013

use strict;

my $testdir = "data";

my @proglist = @ARGV;
my @input_options;
my $fastflag = 0;
my $veryfastflag = 0;
my $veryveryfastflag = 0;
my $found_difference = 0;
my $isoval_offset = 0;
my $outfile = "temp.off";
my $outfile0 = "temp0.off";

my %data_flag;
my $use_all_data = 0;


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

  if ($new_option eq "-hermite") {
    $data_flag{hermite} = 1;
    next;
  }

  if ($new_option eq "-offset") {
    $isoval_offset = shift(@proglist);
    next;
  }

  push(@input_options, $new_option);
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
  run_grad2hermite_on_data();
}
else {
  compare_executables();
  report_on_differences();
}


if (-e "$outfile0") { system "rm $outfile0"; };
if (-e "$outfile") { system "rm $outfile"; }

# ***********************************

sub usage_error {
 print "Usage: testgrad2hermite.pl [OPTIONS] {prog1} {prog2} [{prog3} ...]";
 exit(10);
}

# run grad2hermite
sub run_grad2hermite_on_data {

  my @option_list = @input_options;

  foreach my $tdata (keys %testdata) {

    my $tfile = "$testdir"."/"."$testdata{$tdata}{fname}";
    my @isovalue_list = @{$testdata{$tdata}{isovalue}};

    foreach my $isoval (@isovalue_list) {
      run_grad2hermite($prog0, $tfile, "$outfile0", \@option_list, $isoval);


      my $command_line = "runisodual.pl -normal $outfile0 $isoval $tfile";
      print "$command_line\n";
      system("$command_line") == 0 ||
        die "Program runisodual.pl abnormally terminated.\n";
    }
  }

}


sub compare_executables {

  my @option_list = @input_options;

  foreach my $tdata (keys %testdata) {

    my $tfile = "$testdir"."/"."$testdata{$tdata}{fname}";
    my @isovalue_list = @{$testdata{$tdata}{isovalue}};

    foreach my $isoval (@isovalue_list) {

      run_grad2hermite($prog0, $tfile, "$outfile0", \@option_list, $isoval);
      foreach my $isodual (@proglist) {
        run_grad2hermite($isodual, $tfile, "$outfile", \@option_list, $isoval);
        diff_files("$outfile0", "$outfile");
      }
    }
  }
}

sub run_grad2hermite {

  scalar(@_) == 5 ||
    die "Error in sub run_isodual3D. Requires 4 parameters.\n";

  my $grad2hermite = $_[0];
  my $input_filename = $_[1];
  my $output_filename = $_[2];
  my @option_list = @{$_[3]};
  my $isoval1 = $_[4];

  my $gradient_filename = $input_filename;
  $gradient_filename =~ s/.nrrd/.grad.nrrd/;

  my $command_line;
  $command_line = "./$grad2hermite";
  $command_line = 
    "$command_line" . " @option_list $isoval1 $input_filename $gradient_filename $output_filename";

  print "$command_line\n";
  system("$command_line") == 0 ||
    die "Program isodual3D abnormally terminated.\n";

}

sub diff_files {

  scalar(@_) == 2 ||
    die "Error in sub run_grad2hermite. Requires exactly 2 parameters.\n";

  my $flag = system("diff --brief $_[0] $_[1] > /dev/null");
  if ($flag != 0) { 
    print "*** Output .off files differ. ***\n";
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

