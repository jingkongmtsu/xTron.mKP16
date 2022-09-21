#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use File::Basename;
use File::Copy;

sub getAtomic;

#####################################################
#  this program is to run GAMESS for reference data 
#  used to compare our data
#####################################################

# we need to pass in some args
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
	print "Usage: runtest.pl scf_method basis_name\n";
	exit;
}

# read in method and basis and check
my $METHOD = $ARGV[0];
my $BASIS  = $ARGV[1];

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
	print "basis set is $BASIS\n";
	$basname = "aug-cc-pvtz";
}elsif ($BASIS eq "6-31G**") {
	print "basis set is $BASIS\n";
	$basname = "6-31Gss";
}else{
	print "the input basis sets is $BASIS\n";
	print "the supported basis set are list here below: \n";
	print "Dunning-type: aug-cc-pvtz\n";
	print "Pople-type: 6-31G**\n";
	die "the input basis set name is not supported\n";
}

# form a list of all of input files
my @inputs  = glob("./input/*.in");

# set up a scratch dir
# we only do it when restart is not performing
my $scr = "hf_scr";
if ($METHOD eq "B3LYP") {
	$scr = "b3lyp_scr";
}
$scr = "$scr"."$basname";
if (-e "$scr") {
	system("rm -rf $scr");
}
mkdir ("$scr", 0777) || die "can not creat the scr dir!!\n";

####################################
#     now generate the input       #
#################################### 
my $in;
foreach $in(@inputs) { 

	# set the job 
	my $name = fileparse($in,".in");
	my $job  = $name;

	# is it open shell? etc.
	my $charge = 0;
	my $mult   = 1;
	my @coords;
	my @atoms;
	open (READ, "$in") || die "can not open the $in file for reading geometry data!!!\n";
	while(<READ>) {

		# this is the starting of the geometry section
		if ($_ =~ /molecule/ && $_ =~ /%/) {

			# begin to read on line
			my $line = <READ>;
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;

			my @tmp = split(/\s+/, $line);
			$charge = int($tmp[0]);
			$mult   = int($tmp[1]);

			# now read in the coord
			while(<READ>) {

				# break out the loop
				if ($_ =~ /end/) {
					last;
				}

				# read in the coordinates
				my $line = $_;
				$line =~ s/^\s+//;
				$line =~ s/\s+$//;
				my @tmp = split(/\s+/, $line);
				push(@atoms,  $tmp[0]);
				push(@coords, $tmp[1]);
				push(@coords, $tmp[2]);
				push(@coords, $tmp[3]);
			}

			# we should break out
			last;
		}
	}
	close(READ);

	# set input file
	my $input = "$scr/"."$name".".inp";
	open (WRITE, ">$input") || die "can not open the gamess $input file for writing!!!\n";

	# scf field
	my $scftyp = "RHF";
	if ($mult != 1) {
		$scftyp = "UHF";
	}

	# DFT type
	# now it's non
	my $dftype = "NONE";
	if ($METHOD eq "B3LYP") {
		$dftype   = "B3LYPV1R";
	}

	# run type
	my $runtyp = "ENERGY";

	# whether the 
	my $useSphere = 0;
	if ($BASIS eq "aug-cc-pvtz") {
		$useSphere = 1;
	}

	# keywords section
	print (WRITE  " \$CONTRL SCFTYP=$scftyp DFTTYP=$dftype RUNTYP=$runtyp \$END\n");
	print (WRITE  " \$CONTRL ICHARG=$charge MULT=$mult MAXIT=100 NOSYM=1 \$END\n");
	if ($useSphere == 1) {
		print (WRITE  " \$CONTRL ISPHER=1 \$END\n");
	}
	print (WRITE  " \$SYSTEM MWORDS=200 \$END\n");
	if ($BASIS eq "aug-cc-pvtz") {
		print (WRITE  " \$BASIS  GBASIS=ACCT  \$END\n");
	}elsif ($BASIS eq "6-31G**") {
		print (WRITE  " \$BASIS  GBASIS=N31 NGAUSS=6 \$END\n");
		print (WRITE  " \$BASIS  NDFUNC=1 NPFUNC=1 POLAR=POPN31 \$END\n");
	}else{
		print "the input basis sets is $BASIS\n";
		print "the supported basis set are list here below\n";
		print "Dunning-type: aug-cc-pvtz\n";
		print "Pople-type: 6-31G**\n";
		die "the input basis set is not supported\n";
	}
	print (WRITE  " \$GUESS  GUESS=HUCKEL   \$END\n");
	print (WRITE  " \$SCF DIRSCF=.TRUE. FDIFF=.FALSE. DIIS=.TRUE. EXTRAP=.FALSE. \$END\n");
	if ($METHOD ne "HF") {
		print (WRITE  " \$DFT NRAD=128 NLEB=302 \$END\n");
	}
	print (WRITE  "\n");

	# now let's molecule section
	print (WRITE  " \$DATA\n");
	print (WRITE  "calculation for $job\n");
	print (WRITE  "C1 \n");
	my $atom;
	my $atom_index = 0;
	foreach $atom(@atoms) { 
		my $atomic = getAtomic($atom);
		my $x = $coords[3*$atom_index+0];
		my $y = $coords[3*$atom_index+1];
		my $z = $coords[3*$atom_index+2];
		my $line = "$atom"."  "."$atomic"."  "."$x"."  "."$y"."  "."$z";
		print (WRITE   "$line\n");
		$atom_index = $atom_index + 1;
	}
	print (WRITE  " \$END\n");
	print (WRITE  "\n");
	close(WRITE);
}

####################################
#     now run the calculation      #
#     and parse output             #
####################################
chdir ("$scr") || die "we can not change to the dir of $scr!!\n";
@inputs  = glob("*.inp");
foreach $in(@inputs) { 

	# run gamess to get the output file
	# you may have different way to do it
	system("/home/fenglai/bin/gms", "$in") && die "gamess running for the calculation $in failed\n";

	# now we have the output, get the energy
	my $name = "$in";
	chop($name);
	chop($name);
	chop($name);
	chop($name);
	my $out  = "$name".".log";
	open (READ, "$out") || die "can not open the $out file for parsing result gamess data!!!\n";
	my $E = 0.0;
	while(<READ>) {
		if ($_ =~ /FINAL/) {
			my $line = $_;
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;
			my @tmp = split(/\s+/, $line);
			$E = $tmp[4];
			last;
		}
	}
	close(READ);

	# print out the result
	printf STDOUT ("%-25s  %-16.10f \n", $name, $E);
}
chdir ("../") || die "we can not change to the upper level dir !!\n";
 
##########################
#  utility function      #
##########################
sub getAtomic {

	my $name = $_[0];
	if ($name eq  "H" || $name eq "h") {
		return "1.0";
	}elsif ($name eq  "He" || $name eq "he") {
		return "2.0";
	}elsif ($name eq  "Li" || $name eq "li") {
		return "3.0";
	}elsif ($name eq  "Be" || $name eq "be") {
		return "4.0";
	}elsif ($name eq  "B" || $name eq "b") {
		return "5.0";
	}elsif ($name eq  "C" || $name eq "c") {
		return "6.0";
	}elsif ($name eq  "N" || $name eq "n") {
		return "7.0";
	}elsif ($name eq  "O" || $name eq "o") {
		return "8.0";
	}elsif ($name eq  "F" || $name eq "f") {
		return "9.0";
	}elsif ($name eq  "Ne" || $name eq "ne") {
		return "10.0";
	}elsif ($name eq  "Na" || $name eq "na") {
		return "11.0";
	}elsif ($name eq  "Mg" || $name eq "mg") {
		return "12.0";
	}elsif ($name eq  "Al" || $name eq "al") {
		return "13.0";
	}elsif ($name eq  "Si" || $name eq "si") {
		return "14.0";
	}elsif ($name eq  "P" || $name eq "p") {
		return "15.0";
	}elsif ($name eq  "S" || $name eq "s") {
		return "16.0";
	}elsif ($name eq  "Cl" || $name eq "cl") {
		return "17.0";
	}elsif ($name eq  "Ar" || $name eq "ar") {
		return "18.0";
	}elsif ($name eq  "K"  || $name eq "k") {
		return "19.0"; 
	}elsif ($name eq  "Ca" || $name eq "ca") {
		return "20.0";
	}else{
		print "$name\n";
		die "I did not recognize this atom symbol\n";
	}
}
