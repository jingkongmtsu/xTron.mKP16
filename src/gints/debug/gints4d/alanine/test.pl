#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use File::Basename;
use File::Copy;

# setting the geometry
my $start = 12;
my $final = 13;

# set the basis set file
my $basis = "basis1.txt";

# set input file
my $file = "../large.in";
open (WRITE, ">$file") || die "can not open the $file file for writing!!!\n";

# global section
print (WRITE  "%global_infor\n");
print (WRITE "multi_threads  true\n");
print (WRITE "print_timing   true\n");
print (WRITE  "%end\n");
print (WRITE  "\n");

# now generate the input 
my $d;
for ($d=$start; $d<=$final; $d = $d + 1) {

	# now let's molecule section
	my $xyz  = "ala-"."$d".".xyz";
	open (READ, "$xyz") || die "can not open the $xyz file for writing!!!\n";
	while(<READ>) {
		print (WRITE   "$_");
	}
	close(READ);
	print (WRITE  "\n\n\n");
	open (READ, "$basis") || die "can not open the $basis file for writing!!!\n";
	while(<READ>) {
		print (WRITE   "$_");
	}
	close(READ);
	print (WRITE  "\n\n\n");

	# gints section
	print (WRITE  "%gints\n");
	print (WRITE "gints4d_threshold             1.0E-10\n");
	print (WRITE "shellpair_threshold_gints4d   1.0E-10\n");
	print (WRITE  "%end\n");
	print (WRITE  "\n\n\n");
}
close(WRITE);

