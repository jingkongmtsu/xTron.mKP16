#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use File::Basename;
use File::Copy;

#########################################################
# if you use the matrix diffComp function to print out
# the difference value for every matrix element, here
# you can use this function to collect the output error
# data and print out them into more clear way
#########################################################
#
# set up the input file and error creteria
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
	print "\nUsage: matrix_error_analysis.pl input_file_name error_threshold_value \n";
	exit;
}
my $input    = $ARGV[0];
my $creteria = $ARGV[1];

# whether we omit all of difference value less than the threshold value?
my $we_omit_small_val  = 1;
my $n_large_error_term = 0;  # error that larger than the criteria
my $largestError       = 0;  # this is the largest abs error
my $largestErrorRate   = 0; # corresponding largest abs error rate

# also set an array to store the largest error per each percentage
my @larErrorArray;

# input file
for(my $i=0; $i<=9; $i = $i+1) {

	# set the error base
	my $error_down = 0.1*$i;
	my $error_up   = 0.1*($i+1);
	if ($i == 0) {
		$error_down = 0.01;
	}

	# set the le
	my $le = 0;

	# let's print out  the range
	print "\n\n####################################\n";
	print "# error range is between $error_down and $error_up\n";
	print "####################################\n";

	# read in the input
	open (READ, "$input") || die "can not open the $input file for reading!!!\n";
	while(<READ>) {
		if ($_ =~ /abs difference is/) {
			my @tmp = split /\s+/, $_;
			my $error = abs($tmp[12]);
			my $abs_error = $tmp[8];
			if ($we_omit_small_val == 1) {
				if ($abs_error < $creteria) {
					next;
				}
			}
			if ($error>$error_down && $error<=$error_up) {
				print $_;
				$n_large_error_term += 1;
				if ($abs_error>$le) {
					$le = $abs_error;
				}
			}
		}
	}
	close(READ);

	# now push the le
	push @larErrorArray, $le;
}

# let's print out the largest range
print "\n\n####################################\n";
print "# error range > 1.0 \n";
print "####################################\n";

my $lemax = 0.0;
# read in the input
open (READ, "$input") || die "can not open the $input file for reading!!!\n";
while(<READ>) {
	if ($_ =~ /abs difference is/) {
		my @tmp = split /\s+/, $_;
		my $error = abs($tmp[12]);
		my $abs_error = $tmp[8];
		if ($we_omit_small_val == 1) {
			if ($abs_error < $creteria) {
				next;
			}
		}
		if ($error>=1.0) {
			print $_;
			$n_large_error_term += 1;
			if ($abs_error>$lemax) {
				$lemax = $abs_error;
			}
		}
		if ($abs_error>$largestError) {
			$largestError = $abs_error; 
			$largestErrorRate = $error; 
		}
	}
}
close(READ);

# then let's grasp some useful information
my $do_check = 1;
my $nbas = 0;
my $MAE;
my $RMSD;
if($do_check == 1) {
	open (READ, "$input") || die "can not open the $input file for reading!!!\n";
	while(<READ>) {
		if ($_ =~ /maximum error/) {
			my @tmp = split /\s+/, $_;
			my $t = $tmp[3];
			$MAE  = $t;
		}
		if ($_ =~ /RMSD error/) {
			my @tmp = split /\s+/, $_;
			my $t = $tmp[3];
			$RMSD  = $t;
		}
		if ($_ =~ /number of basis sets/) {
			my @tmp = split /\s+/, $_;
			my $t = $tmp[4];
			$nbas  = $t;
		}
	}
	close(READ);
}

print "maximum error is ", $MAE, "\n";
print "RMSD is ", $RMSD, "\n";
print "the largest absolute error is ", $largestError, " and the corresponding rate is ", $largestErrorRate, "\n";
print "total number of large error terms is ", $n_large_error_term, "\n";
my $n = ($nbas*($nbas+1))/2;
if ($n > 0) {
	my $rate = $n_large_error_term/$n*100.0;
	print "nbas", $nbas,"\n";
	print "total number of large error terms rate in percentage ", $rate, "\n";
}
print "the largest error for each percentage\n";
my $i = 0;
my $e;
foreach $e(@larErrorArray) {
	my $error_down = 0.1*$i;
	my $error_up   = 0.1*($i+1);
	if ($i == 0) {
		$error_down = 0.01;
	}
	print "error range is between $error_down and $error_up : ", $e, "\n";
	$i = $i+1;
}
print "error range is larger than 1.0 : ", $lemax, "\n";

