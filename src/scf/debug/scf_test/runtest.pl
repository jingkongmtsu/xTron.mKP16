#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use File::Basename;
use File::Copy;

#####################################################
#  this program is to run our calculation against
#  the reference data
#####################################################

# we need to pass in some args
my $num_args = $#ARGV + 1;
if ($num_args != 3 && $num_args != 4) {
	print "\nUsage: runtest.pl xtron_path scf_method basis_name (restart)\n";
	exit;
}

# read in binary path, method and basis and check
my $XTRON  = $ARGV[0];
my $METHOD = $ARGV[1];
my $BASIS  = $ARGV[2];

# do we do restart?
my $redo   = 0;
if ($num_args == 4) {
	my $val = $ARGV[3];
	if ($val eq "restart") {
		$redo = 1;
	}else{
		print "the third argument is not correct!\n";
		exit;
	}
}

# check that whether the xtron exist?
my $binary = "$XTRON"."/xtron.exe";
if (-e "$binary") {
	print "the xtron.exe is found at the path of $XTRON\n";
}else{
	die "we can not locate the binary xtron.exe through the path you provided: $XTRON\n";
}

# check method
$METHOD =~ tr/a-z/A-Z/;
if ($METHOD eq "HF") {
	print "checking HF method\n";
}elsif ($METHOD eq "B3LYP") {
	print "checking B3LYP method\n";
}else {
	print "the input method is $METHOD\n";
	print "We only support HF/B3LYP for testing\n";
	die "the input method is not supported\n";
}

# check basis set
my $basname; 
if ($BASIS eq "aug-cc-pvtz") {
	$basname = "$BASIS";
}elsif ($BASIS eq "6-31G**") {
	$basname = "6-31Gss";
}elsif ($BASIS eq "g3large") {
	$basname = "g3large";
}else{
	print "the input basis sets is $BASIS\n";
	print "the supported basis set are list here below: \n";
	print "Dunning-type: aug-cc-pvtz\n";
	print "Pople-type: 6-31G** and g3large\n";
	die "the input basis set name is not supported\n";
}

# whether this is DFT method?
my $is_DFT = 0;
if ($METHOD eq "B3LYP") {
	$is_DFT = 1;
}

# creteria comparing with reference
# for DFT, because we have difference on setting up
# the grids etc. between GAMESS and ours, therefore
# for multi-atoms system we expect the absolute energy
# difference should increase
# so we set creteria to be 0.0001
my $creteria = 0.000001;
if ($is_DFT == 1) {
	$creteria = 0.0001;
}

# form a list of all of input files
my @inputs  = glob("./input/*.in");

# set up a scratch dir
my $scr = "our_scr";
if ($redo == 0) {
	if (-e "$scr") {
		system("rm -rf $scr");
	}
	mkdir ("$scr", 0777) || die "can not creat the scr dir!!\n";
}

####################################
#     now generate the input       #
#################################### 
my $in;
foreach $in(@inputs) { 

	# set the job 
	my $name = fileparse($in,".in");

	# set input file
	my $input = "$scr/"."$name".".in";
	my $out   = "$scr/$name".".log";

	# let's see whether we need to generate the input file and do calculation
	my $doCal = 1;
	if ($redo == 1) {
		if (-f "$out") {
			$doCal = 0;
		}
	}

	# write the file
	if ($doCal == 1) {
		open (WRITE, ">$input") || die "can not open the $input file for writing!!!\n";

		# read in the input
		open (READ, "$in") || die "can not open the $in file for reading!!!\n";
		while(<READ>) {
			print (WRITE  $_);
		}
		close(READ);
		print (WRITE  "\n");

		# basis set file
		my $bas = "bas/"."$basname".".bas";
		open (READ, "$bas") || die "can not open the $bas file for reading!!!\n";
		while(<READ>) {
			print (WRITE $_);
		}
		close(READ);

		# now it's keyword
		print (WRITE  "\n");
		print (WRITE  "%xcfunc\n");
		print (WRITE  "name $METHOD\n");
		print (WRITE  "%end\n");
		print (WRITE  "\n");

		print (WRITE  "%scf\n");
		print (WRITE  "max_scf_cycles = 50\n");
		# for some special cases, core guess does not lead to a correct
		# scf result, so we drop it
		print (WRITE  "scf_guess = core\n");
		print (WRITE  "scf_convergence  = 1.0E-6\n");
		print (WRITE  "%end\n");	
		print (WRITE  "\n");

		print (WRITE  "%scfconv\n");
		print (WRITE  "scf_algorithm  ediis\n");
		print (WRITE  "%end\n");
		print (WRITE  "\n");

		print (WRITE  "%xcints\n");
		print (WRITE  "grid_points  128  302\n");
		print (WRITE  "%end\n");
		print (WRITE  "\n");

		print (WRITE  "%gints\n");
		print (WRITE  "shellpair_threshold_gints4d  1.0E-14\n");
		print (WRITE  "shellpair_threshold_gints2d  1.0E-14\n");
		print (WRITE  "gints4d_threshold            1.0E-14\n");
		print (WRITE  "gints2d_threshold            1.0E-14\n");
		print (WRITE  "cs_threshold                 1.0E-14\n");
		print (WRITE  "gints4deriv_threshold        1.0E-14\n");
		print (WRITE  "gints2deriv_threshold        1.0E-14\n");
		print (WRITE  "%end\n");
		print (WRITE  "\n");

		print (WRITE  "%scfintscontroller\n");
		print (WRITE  "integral_controller_option  NO_CONTROLLER\n");
		print (WRITE  "%end\n");
		print (WRITE  "\n");
		close(WRITE);

		# now run the calculation
		# you may have different way to do it
		my $out  = "$scr/$name".".log";
		system("$binary", "$input", "$out") && die "running the calculation $input failed\n";
	}

	# get the energy from reference data
	my $ref = "gamess_";
	if ("$METHOD" eq "HF" || "$METHOD" eq "hf") {
		$ref = "$ref"."hf"."_";
	}elsif ("$METHOD" eq "B3LYP") {
		$ref = "$ref"."B3LYP"."_";
	}
	$ref = "$ref"."$basname".".txt";
	my $E0 = 0.0;
	open (READ, "$ref") || die "can not open the $ref file for reading!!!\n";
	my $get_it = 0;
	while(<READ>) {
		my $line = $_;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		my @tmp = split(/\s+/, $line);
		my $t   = $tmp[0];
		if ($t eq $name) {
			$E0 = $tmp[1];
			$get_it = 1;
			last;
		}
	}
	close(READ);
	if ($get_it == 0) {
		die "failed to get the standard data from reference data\n";
	}

	# now let's compare
	# we note that the total energy value is always the last value
	my $E = 0.0;
	open (READ, "$out") || die "can not open the $out file for reading!!!\n";
	$get_it = 0;
	while(<READ>) {
		if ($_ =~ /@@@@/) {
			my $line = $_;
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;
			my @tmp = split(/\s+/, $line);
			$E = $tmp[$#tmp];
			$get_it = 1;
		}
	}
	close(READ);
	if ($get_it == 0) {
		die "failed to get the result data from output\n";
	}

	# finally let's see 
	my $diff = abs($E-$E0);
	#print "our energy is $E", " standard energy is $E0", " diff is $diff\n";
	if ($diff>$creteria) {
		printf ("%-25s is failed for difference %-16.10f \n", $name, $diff);
	}else{
		printf ("%-25s is pass with difference %-16.10f!! \n", $name, $diff);
	}
}

