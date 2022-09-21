#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

######################################################
# This program is used to generate the basis set     #
# order, which will be used in integral generation   #
# fenglai liu                                        #
# Jan. 2012                                          #
######################################################

# max L to generate basis set order
my $maxL = 20;

# the type of basis set order, right now two types are 
# available: LIBINT_ORDER and TESTING_ORDER
# in default, we generate basis set order with LIBINT_ORDER
my $basis_set_order = "LIBINT_ORDER";

#analyze the command line parameter
if ("$#ARGV" == 0) {
	$maxL = $ARGV[0];
}elsif("$#ARGV" == 1) {
	$maxL = $ARGV[0];
	$basis_set_order = $ARGV[1];
} else {
	print "We use default parameters\n";
	print "Default maxL is 20, and basis set order we follow libint program\n";
}

# test the parameter
if ($maxL < 0) {
	die "The maximum L provided should not be lower than 0\n";
}

# prepare file writing
my $file = "basisSetOrder.txt";
open (WRITE, ">$file") || die "can not open the $file file for writing!!!\n";

#generating the basis set order
my $i;
my $j;
my $k;
my $am;
my $nx; 
my $ny; 
my $nz; 
if ($basis_set_order eq "LIBINT_ORDER") {
	for($am=0; $am<= $maxL; $am++) {
		for($i=0; $i<=$am; $i++) {
			$nx = $am - $i;
			for($j=0; $j<=$i; $j++) {
				$ny = $i-$j; 
				$nz = $j;   
				if ($ny == 0 && $nz == 0) {
					print (WRITE "$nx,  $ny,  $nz,  /*L = $nx*/\n");
				}else{
					print (WRITE "$nx,  $ny,  $nz,\n");
				}
			}
		}
	}
}elsif ($basis_set_order eq "TESTING_ORDER"){
	for($am=0; $am<= $maxL; $am++) {
		for($i=0; $i<=$am; $i++) {
			$nz = $i;
			for($j=0; $j<=$am-$i; $j++) {
				$ny = $j;
				$nx = $am - $nz - $ny;
				if ($nx == 0 && $ny == 0) {
					print (WRITE "$nx,  $ny,  $nz,  /*L = $nz*/\n");
				}else{
					print (WRITE "$nx,  $ny,  $nz,\n");
				}
			}   
		}   
	}   
} else {
	die ("basis set order is not correct");
}
close(WRITE);


