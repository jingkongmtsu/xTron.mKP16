#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

#########################################################
#  this perl file is going to explore the maximum number
#  of shells contained in the basis set data, as well as
#  maximum number of contractions for a given shell
#########################################################

# set the basis set style 
my $GAUSSIAN_STYLE = 1;
my $NWCHEM_STYLE   = 2;

# get the files list
my @files = glob("*.bas");
#my @files = glob("basis_lib/*");

# these are the limits we may used in the program
# let's try to print out the basis set larger than this limits
my $nPrimsLimit  = 20;
my $nShellsLimit = 30;

# if the basis set file has G, H etc. angular momentum, print it out
my $ang_limit    = 3;

# set the global compare result
my $maxnshells = 0;
my $maxnprims  = 0;

# print out the details for each basis set file
my $f;
foreach $f(@files) {

	# set the data
	my $maxPrim   = 0;
	my $maxNShell = 0;

	# firstly, we need to judge which style the basis set 
	# data file it is?
	my $style = 0;
	open (READ, "$f") || die "can not open the basis set data file $f for reading!!!\n";
	while(<READ>) {
		my $line = $_;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		if ($line =~ /^BASIS/i ) {
			$style = $NWCHEM_STYLE;
			last;
		}
		if ($line =~ /^\*\*\*\*/ ) {
			$style = $GAUSSIAN_STYLE;
			last;
		}
	}
	close(READ);

	# check that whether we have got the style?
	if ($style == 0) {
		print "We failed to identify the basis set data style, only Gaussian/NWChem are supported\n";
		exit(1);
	}

	# parse the nwchem style of file
	if ($style == $NWCHEM_STYLE) {

		# set the maximum ang etc. for this file
		my $maxAng    = 0;
		my $maxShells = 0;
		my $maxPrims  = 0;

		# read in data
		open (READ, "$f") || die "can not open the basis set data file $f for reading!!!\n";
		while(<READ>) {

			# we omit all of comment lines
			if ($_ =~ /#/ ) {
				next;
			}

			# get a section of atom shell data
			if ($_ =~ /^BASIS/i ) {
				my $nshells = 0;
				my $nprim   = 0;

				# for the first shell data, we need to read it's shelltype
				my $line = <READ>;
				$line =~ s/^\s+//;
				$line =~ s/\s+$//;
				my @tmp = split /\s+/, $line;
				my $shelltype = $tmp[1];

				# now let's read in the shell data, until the next shell title
				# or the end of the atom shell section
				while(<READ>) {

					# shall we stop here?
					if ($_ =~ /^END/i ) {
						last;
					}

					# trim the data string
					# split the data according to space
					my $s = $_;
					$s =~ s/^\s+//;
					$s =~ s/\s+$//;
					my @tmp = split /\s+/, $s;

					# the title of each shell starts with atom symbol
					# and the shell type, which is what we needed
					if ($tmp[0] !~ /^[0-9,.]/) {
						$shelltype = $tmp[1];

						# also let's compare number of primitives
						if ($nprim>$maxPrim) {
							$maxPrim  = $nprim;
						}
						if ($nprim>$maxPrims) {
							$maxPrims = $nprim;
						}
						$nprim = 0;

						# let's go to next data section
						next;
					}

					# let's count the number of shells
					# we only do it when we read the first line of primitive function data,
					# that is to say; nprim = 0
					#
					# for the SP shell, they do not support the 
					# general format of shell data; therefore for 
					# SP shell the nshells are 1
					if ($nprim == 0) {
						if ($shelltype eq "SP" || $shelltype eq "sp") {
							$nshells += 1;
						}else{
							$nshells += scalar(@tmp) - 1;
						}

						# let's see the shell type
						# we are going to derive the maxAng
						my $ang = 0;
						if ($shelltype eq "S" || $shelltype eq "s") {
							$ang = 0;
						}elsif ($shelltype eq "SP" || $shelltype eq "sp") {
							$ang = 1;
						}elsif ($shelltype eq "P"  || $shelltype eq "p") {
							$ang = 1;
						}elsif ($shelltype eq "D"  || $shelltype eq "d") {
							$ang = 2;
						}elsif ($shelltype eq "F"  || $shelltype eq "f") {
							$ang = 3;
						}elsif ($shelltype eq "G"  || $shelltype eq "g") {
							$ang = 4;
						}elsif ($shelltype eq "H"  || $shelltype eq "h") {
							$ang = 5;
						}elsif ($shelltype eq "I"  || $shelltype eq "i") {
							$ang = 6;
						}else{
							print "$_\n";
							print "the shell type is not supported, the value is ", $shelltype, "\n";
							die "illegal shell type\n";
						}

						# set the maxAng
						if ($maxAng<$ang) {
							$maxAng = $ang;
						}
					}

					# counting the primitive functions
					$nprim = $nprim + 1;
				}

				# at the end of this atom shell data, let's see how many
				# shells we counted
				if ($nshells>$maxNShell) {
					$maxNShell = $nshells;
				}
				if ($nshells>$maxShells) {
					$maxShells = $nshells;
				}
			}
		}
		close(READ);

		# print out the max information
		print "*************************************************************\n";
		print "basis set file name: ", $f, "\n";
		print "maximum number of shells     in atom: ", $maxShells, "\n";
		print "maximum number of primitives in atom: ", $maxPrims, "\n";
		print "highest angular momentum is ", $maxAng, "\n";
		print "*************************************************************\n\n";
	}

	# parse the Gaussian style of file
	if ($style == $GAUSSIAN_STYLE) {

		# set the maximum ang for this file
		my $maxAng    = 0;
		my $maxShells = 0;
		my $maxPrims  = 0;

		# read in data
		open (READ, "$f") || die "can not open the basis set data file $f for reading!!!\n";
		while(<READ>) {

			# get a section of atom shell data
			my $line = $_;
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;
			my @tmp = split /\s+/, $line;

			# let's see whether this is the atom symbol with a "0"
			my $nshells = 0;
			if (scalar(@tmp) == 2 && $tmp[1] eq "0") {
				while(<READ>) {

					# we break when the shell sign met
					if ($_ =~ /\Q****\E/) {
						last;
					}

					# now let's read the shell information
					my $l = $_;
					$l =~ s/^\s+//;
					$l =~ s/\s+$//;
					my @tmp1 = split /\s+/, $l;

					# now set the data
					if ($tmp1[0] !~ /^[0-9,.]/) {
						$nshells += 1;
						my $shelltype = $tmp1[0];
						my $nprim = $tmp1[1];

						# also let's compare number of primitives
						if ($nprim>$maxPrim) {
							$maxPrim = $nprim;
						}
						if ($nprim>$maxPrims) {
							$maxPrims = $nprim;
						}

						# analyze the shell type
						my $ang = 0;
						if ($shelltype eq "S" || $shelltype eq "s") {
							$ang = 0;
						}elsif ($shelltype eq "SP" || $shelltype eq "sp") {
							$ang = 1;
						}elsif ($shelltype eq "P"  || $shelltype eq "p") {
							$ang = 1;
						}elsif ($shelltype eq "D"  || $shelltype eq "d") {
							$ang = 2;
						}elsif ($shelltype eq "F"  || $shelltype eq "f") {
							$ang = 3;
						}elsif ($shelltype eq "G"  || $shelltype eq "g") {
							$ang = 4;
						}elsif ($shelltype eq "H"  || $shelltype eq "h") {
							$ang = 5;
						}elsif ($shelltype eq "I"  || $shelltype eq "i") {
							$ang = 6;
						}else{
							print "$_\n";
							print "the shell type is not supported, the value is ", $shelltype, "\n";
							die "illegal shell type\n";
						}

						# set the maxAng
						if ($maxAng<$ang) {
							$maxAng = $ang;
						}
					}
				}
			}

			# at the end of this atom shell data, let's see how many shells
			if ($nshells>$maxNShell) {
				$maxNShell = $nshells;
			}
			if ($nshells>$maxShells) {
				$maxShells = $nshells;
			}
		}
		close(READ);

		# print out the max information
		print "*************************************************************\n";
		print "basis set file name: ", $f, "\n";
		print "maximum number of shells     in atom: ", $maxShells, "\n";
		print "maximum number of primitives in atom: ", $maxPrims, "\n";
		print "highest angular momentum is ", $maxAng, "\n";
		print "*************************************************************\n\n";
	}

	# now let's print out
	if ($maxNShell > $nShellsLimit) {
		print "basis set file name: ", $f, "\n";
		print "maximum number of shells for an atom: ", $maxNShell, "\n";
	}
	if ($maxPrim > $nPrimsLimit) {
		print "basis set file name: ", $f, "\n";
		print "maximum number of primitive functions for a shell: ", $maxPrim, "\n\n";
	}

	# make a record
	if ($maxnshells < $maxNShell) {
		$maxnshells = $maxNShell;
	}
	if ($maxnprims < $maxPrim) {
		$maxnprims = $maxPrim;
	}
}	

# print out the global data
print "maximum number of shells for an atom among all of basis set files: ", $maxnshells, "\n";
print "maximum number of primitive functions for a shell among all of basis set files: ", $maxnprims, "\n";

