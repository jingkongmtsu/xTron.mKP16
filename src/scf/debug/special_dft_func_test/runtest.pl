#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use File::Basename;
use File::Copy;

#####################################################
#  this program is to run our calculation so that 
#  to test current implementation against reference
#####################################################

# we need to pass in some args
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
	print "\nUsage: runtest.pl scf_method basis_name \n";
	exit;
}

# read in method and basis and check
my $METHOD = $ARGV[0];
my $BASIS  = $ARGV[1];

# check method
$METHOD =~ tr/a-z/A-Z/;
if ($METHOD eq "KP14") {
	print "checking KP14 method\n";
}else {
	print "the input method is $METHOD\n";
	print "We only support KP14 for testing\n";
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
if ($METHOD eq "KP14") {
	$is_DFT = 1;
}

# creteria comparing with reference
my $creteria = 0.00001;

# form a list of all of input files
my @inputs  = glob("./input/*.in");

# set up a scratch dir
my $scr = "test_scr";
if (-e "$scr") {
	system("rm -rf $scr");
}
mkdir ("$scr", 0777) || die "can not creat the scr dir!!\n";

# set up the reference data
my $ref = "ref_";
if ("$METHOD" eq "B05") {
	$ref = "$ref"."B05"."_";
}elsif ("$METHOD" eq "KP14") {
	$ref = "$ref"."KP14"."_";
}
$ref = "$ref"."$basname".".txt";

####################################
#     now generate the input       #
#################################### 
my $in;
foreach $in(@inputs) { 

	# set the job 
	my $name = fileparse($in,".in");

	# set input file
	my $input = "$scr/"."$name".".in";
	open (WRITE, ">$input") || die "can not open the $input file for writing!!!\n";

	# read in the input, make up the first geom
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
	print (WRITE  "name HF\n");
	print (WRITE  "%end\n");
	print (WRITE  "\n");

	print (WRITE  "%scf\n");
	print (WRITE  "max_scf_cycles = 200\n");
	print (WRITE  "scf_convergence  = 1.0E-6\n");
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
	print (WRITE  "\n\n\n");

	# read in the input
	open (READ, "$in") || die "can not open the $in file for reading!!!\n";
	while(<READ>) {
		print (WRITE  $_);
	}
	close(READ);
	print (WRITE  "\n");

	# basis set file
	$bas = "bas/"."$basname".".bas";
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
	print (WRITE  "max_scf_cycles = 200\n");
	print (WRITE  "scf_convergence  = 1.0E-5\n");
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
	system("./SCF", "$input", "$out") && die "running the calculation $input failed\n";

	# now we have the output, get the energy
	open (READ, "$out") || die "can not open the $out file for parsing result data!!!\n";
	my $E = 0.0;
	my $getit = 0;
	while(<READ>) {
		if ($_ =~ /scf converged energy is:/) {
			my $line = $_;
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;
			my @tmp = split(/\s+/, $line);
			$E = $tmp[8];
			$getit = $getit + 1;
		}
	}
	close(READ);
	if ($getit != 2) {
		print "something wrong with SCF convergence for the case $name\n";
		die "please check the SCF convergence!!!\n";
	}

	# now get the reference
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

	# finally let's see 
	my $diff = abs($E-$E0);
	#print "our energy is $E", " standard energy is $E0", " diff is $diff\n";
	if ($diff>$creteria) {
		printf ("%-25s is failed for difference %-16.10f \n", $name, $diff);
	}else{
		printf ("%-25s is pass with difference %-16.10f!! \n", $name, $diff);
	}
}

