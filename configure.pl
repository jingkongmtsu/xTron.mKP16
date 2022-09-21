#
# note:
#
# we currently disable the checking for GSL
# since the SCF do not use GSL library anymore
#
#!/usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use Cwd;
use Env;

#
# this program is used to instruct cmake to compiling the whole program
# author: fenglai liu
#

#################################
#   the parameters set here     #
#################################

# compilers
# we do not have C compiler inside
# All are C++ plus fortran
my $cxx_compiler = "g++";
my $for_compiler = "gfortran";
my $disable_compiler_warning = 0;

# we may need the library from intel
# if we use the Intel compiler
my $intel_path = "";
if (defined $ENV{"INTEL_ROOT"} && length $ENV{"INTEL_ROOT"} > 0) {
	$intel_path = $ENV{"INTEL_ROOT"};
}

# build dir and source dir
# source dir is always default current dir
my $current_dir= getcwd();
my $source_dir = $current_dir; 
my $build_dir  = "$current_dir"."/"."build";

# build type
my $build_type = "Release";

# blas/lapack vendor
my $blas_vendor = "netlib";

# netlib_root
my $netlib_root = "";
if (defined $ENV{"NETLIB_ROOT"} && length $ENV{"NETLIB_ROOT"} > 0) {
	$netlib_root = $ENV{"NETLIB_ROOT"};
}

# mkl_root
my $mkl_root = "";
if (defined $ENV{"MKLROOT"} && length $ENV{"MKLROOT"} > 0) {
	$mkl_root = $ENV{"MKLROOT"};
}
my $mkl_major_version  = -1;
my $use_sequential_mkl =  1;

# check TBBROOT
my $tbb_root = "";
if (defined $ENV{"TBBROOT"} && length $ENV{"TBBROOT"} > 0) {
	$tbb_root = $ENV{"TBBROOT"};
}

# we need to have the boost dir to be set
my $boost_dir = "";
if (defined $ENV{"BOOST_ROOT"} && length $ENV{"BOOST_ROOT"} > 0) {
	$boost_dir = $ENV{"BOOST_ROOT"};
}

# we need to have the gsl dir to be set
my $gsl_dir = "";
if (defined $ENV{"GSLROOT"} && length $ENV{"GSLROOT"} > 0) {
	$gsl_dir = $ENV{"GSLROOT"};
}

# now set the building target
# for debugging modules, we have "DEBUG" + module name,
# for example; target = DSCF means debugging scf module
my $target;

# whether we should compile with integral derivatives code?
# in default we do not, this is because it's too time consuming
# we will check the target defined by the user, then see
# whether we bring up the integral derivatives code
# if the value is 0, we do not do deriv
# if the value is 1, we will do it 
my $with_gints_deriv = 0;

# specify the external integral library path
# and to see whether we use the external integral library
# 0 is not, which is default;
# 1 is yes to use it. and you must define the integral lib path
my $integral_lib_path;
my $use_external_int_lib = 0;

# do you want to use the GDM?
my $use_gdm = 1;
my $gdm_lib_path;

# setting for using Intel Phi co-processor
# in default it's not enabled
my $use_mic = 0;
my $check_local_mic_infor = -1;

#------------------------------------------------------------------------------
# print out help message if no argument is assigned 
#------------------------------------------------------------------------------
sub help_print {
	print "\n";
	print "Here below is the options to build this program.\n";
	print "Command argument format is like: compiler=intel\n";
	print "\n";
	print "cxx_compiler: (default choice is gnu)\n";
	print " * gnu      -- GNU g++ compilers \n";
	print " * intel    -- Intel icpc compilers \n";
	print "\n";
	print "for_compiler: (default choice is gnu)\n";
	print " * gnu      -- GNU gfortran compilers \n";
	print " * intel    -- Intel ifort compilers \n";
	print "\n";
	print "warning_off : whether to disable the warning message of compiler(default choice is false)\n";
	print "\n";
	print "build_dir   : binary file building directory\n";
	print "            default is build folder under current DIR\n";
	print "\n";
	print "source_dir  : source code directory\n";
	print "            default is current DIR\n";
	print "\n";
	print "build_type  : building options for compilers, default is release with optimization on code\n";
	print " * debug    -- compiling with debug choice\n";
	print " * release  -- compiling with optimization for release version\n";
	print "\n";
	print "debug_module: building binary for debugging the given module, see the list below\n";
	print " * gints    -- debug gints module  (main function is in gints/debug/main.cpp)\n";
	print " * xcints   -- debug xcints module (main function is in xcints/debug/main.cpp)\n";
	print " * scf      -- debug scf module    (main function is in scf/debug/main.cpp)\n";
	print " * result   -- debug result module (main function is in result/debug/main.cpp)\n";
	print "\n";
	print "release_module: building binary for real calculation, all of main functions under src/main folder\n";
	print " * scf      -- building code for single point calculation (main function is in main/scf/sp.cpp)\n";
	print "\n";
	print "blas_vendor : BLAS/LAPACK library provider (default is mkl)\n";
	print " * mkl      -- MKL library from intel(env. MKLROOT should be defined)\n";
	print " * netlib   -- library from http://www.netlib.org/\n";
	print "\n";
	print "mkl_type    : if you use mkl library, do you use parallel one or sequential one(default is sequential)\n";
	print " * seq      -- represents the sequential library of MKL\n";
	print " * par      -- represents the parallel version of MKL\n";
	print "\n";
	print "boost_dir   : BOOST dir provided by the user, in default it's empty\n";
	print "              this value will be used by cmake in searching boost library\n";
	print "\n";
	print "intel_dir   : Intel dir provided by the user, in default it's empty\n";
	print "              this value will be used by cmake in searching Intel library\n";
	print "\n";
	print "gsl_dir     : GNU GSL library dir provided by the user, in default it's empty\n";
	print "              this value will be used by cmake in searching gnu gsl library\n";
	print "\n";
	print "int_lib_path: because the analytical integral library compilation costs most of the time\n";
	print "              therefore you can pre-compile it into library and use it directly here so\n";
	print "              that to save compilation time. This keyword specify the path of integral\n";
	print "              library. The analytical integral includes ERI and ERI derivatives etc. \n";
	print "\n";
	print "gdm_lib_path: if you want to use GDM algorithm for handle SCF convergence, please\n";
	print "              specify the library path here\n";
	print "\n";
	exit 0;
}

#------------------------------------------------------------------------------
# parsing the input parameters
# please make sure that the command line we parsed matches the above
# information in help print function 
#------------------------------------------------------------------------------
sub para_parsing {
	my $cml_para;
	foreach $cml_para(@ARGV) {

		# separate the keyword and value 
		# for parameter we have to use 
		# "=" to connect the keyword and its value
		# the only exception is help parameter,
		# which has no parameter
		my $key;
		my $val;
		if ($cml_para =~ /=/) {
			my @para = split /=/, $cml_para;
			my $len  = @para;
			if ($len == 2) {
				$key = $para[0];
				$val = $para[1];
			}else{
				print "Unkown command argument is given: ", $cml_para, "\n";
				die "Invalid command argument is given\n";
			}
		}else {

			# is it help message?
			if ($cml_para eq "help") {
				help_print();
			}	

			# now print out error message
			print "Unkown command argument is given: ", $cml_para, "\n";
			print "The correct format is like below:\n";
			print "./configure.pl compiler=gnu\n";
			die "Invalid command argument is given\n";
		}

		# now let's set the arguments
		if ($key eq "cxx_compiler") {

			# it's compiler setting
			if ($val eq "gnu") {
				$cxx_compiler = "g++";
			}elsif ($val eq "intel") {
				$cxx_compiler = "icpc";
			}else{
				print "Your compiler option is: ", $val, "\n";
				die "Invalid C++ compiler option is given\n";
			}

		}elsif ($key eq "for_compiler") {

			# it's compiler setting
			if ($val eq "gnu") {
				$for_compiler = "gfortran";
			}elsif ($val eq "intel") {
				print "**********************************************************\n";
				print "                         WARNING!!!                       \n";
				print "currently we have trouble in using the ifort              \n";
				print "it seems for PSTS functional, no matter whether we link\n   ";
				print "with ifcore or ifcoremt, the XC matrix result for DFT\n     ";
				print "is just like a mess. However, it does not have any problem\n";
				print "when we switch to GFORTRAN. Therefore if you choose ifort \n";
				print "we will report an error here, and please choose gfortran  \n";
				print "instead\n";
				print "\n";
				print "We note that the testing ifort version is 17.0.2. We do not\n";
				print "know other version results, but in case we just disable all\n";
				print "of use of ifort\n";
				print "**********************************************************\n\n";
				$for_compiler = "ifort";
				die "Invalid fortran compiler option is given\n";
			}else{
				print "Your compiler option is: ", $val, "\n";
				die "Invalid fortran compiler option is given\n";
			}

		}elsif ($key eq "build_dir") {

				$build_dir = $val;

		}elsif ($key eq "source_dir") {

			# now let's try to get the source dir
			$source_dir = $val;

		}elsif ($key eq "build_type") {

			# now let's try to get the build dir
			if ($val eq "debug") {
				$build_type = "Debug";
			}elsif($val eq "release") {
				$build_type = "Release";
			}else{
				print $val, "\n";
				die "Invalid build type given\n";
			}

		}elsif ($key eq "blas_vendor") {

			# now let's try to get the blas_vendor
			$blas_vendor = $val;

		}elsif ($key eq "mkl_type") {

			# check the value
			if ($val eq "seq") {
				$use_sequential_mkl =  1;
			}elsif($val eq "par") {
				$use_sequential_mkl =  0;
			}else{
				print $val, "\n";
				die "Invalid mkl_type value given, only value of \"seq\" or \"par\" are allowed\n";
			}

		}elsif ($key eq "boost_dir") {

			# now let's try to get the boost_dir
			$boost_dir = $val;

		}elsif ($key eq "intel_dir") {

			# now let's try to get the intel path
			$intel_path = $val;

		}elsif ($key eq "gsl_dir") {

			# now let's try to get the gsl_dir
			$gsl_dir = $val;

		}elsif ($key eq "debug_module") {

			# we need to check that whether target has already
			# been defined?
			if (defined $target) {
				print "the target has already been set before define debug_module\n";
				print "we will overwrite the original definition with debug_module\n";
			}

			# below are all debugging modules
			$target = "D";
			if ($val eq "gints") {
				$target = "$target"."GINTS";
			}elsif ($val eq "xcints") {
				$target = "$target"."XCINTS";
			}elsif ($val eq "scf") {
				$target = "$target"."SCF";
			}elsif ($val eq "result") {
				$target = "$target"."RESULT";
			}else{
				die "Invalid target $val for debug module\n";
			}

		}elsif ($key eq "release_module") {

			# we need to check that whether target has already
			# been defined?
			if (defined $target) {
				print "the target has already been set before define release_module\n";
				print "we will overwrite the original definition with release_module\n";
			}

			# below are all debugging modules
			if ($val eq "scf") {
				$target = "SCF";
			}else{
				die "Invalid target $target for release module\n";
			}

		}elsif ($key eq "int_lib_path") {

			# let's read in integral library path
			# we will check it later
			$integral_lib_path = $val;

		}elsif ($key eq "gdm_lib_path") {

			# let's read in gdm library path
			# we will check it later
			$gdm_lib_path = $val;

		}elsif ($key eq "warning_off") {

			if ($val eq "true" || $val eq "TRUE") {
				$disable_compiler_warning = 1;
			}elsif ($val eq "false" || $val eq "FALSE") {
				$disable_compiler_warning = 0;
			}else {
				print "Invalid value for keyword warning_off, only true or false allowed\n";
				die "This value for this keyword from command line is not recognized\n";
			}

		}else {
			print "Invalid keyword: $key\n";
			die "This keyword input from command line is not recognized\n";
		}
	}

	# finally we need to check
	# wether we have fined the target?
	if (not defined $target) {
		print "you did not define any target module for compiling\n";
		print "therefore we choose the default target is SCF \n";
		$target = "SCF";
	}

	# check the MIC card related information
	if ($use_mic == 1) {

		# if mic_is_local not defined, we need to issue an error here
		if ($check_local_mic_infor < 0) {
			print "since you choose to use the MIC card(use_mic is true), you need to specify whether\n";
			print "the MIC card is locally installed in this computer(set mic_is_local to true or false).\n";
			print "If mic_is_local is true, we will perform additional checking for the MIC card information.\n";
			print "For more information, please see the description about mic_is_local in print mode\n";
			die "the setting of mic_is_local is undefined, please use ./configure.pl help to get more information\n";
		}
	}
}

#------------------------------------------------------------------------------
# check the mkl version
#------------------------------------------------------------------------------
sub mkl_inquire {

	# check the mkl_root
	if (length($mkl_root) == 0) {
		print "You did not define the MKLROOT ENV\n";
		die "Can not find the MKL\n";
	}

	# now let's try to get mkl's version number
	# so that later in mkl we can generate 
	# mkl library
	my $file1="$mkl_root"."/"."include"."/"."mkl.h";
	my $file2="$mkl_root"."/"."include"."/"."mkl_version.h";
	my @files= ("$file1", "$file2");
	my $f;
	my $get_it = 0;
	for $f(@files) {

		# the file may not exist
		if (! -e "$f") {
			next;
		}

		# try to catch the information
		open(READ, "$f") || die "can not open the $f through the given path $mkl_root!\n";
		while(<READ>) {
			if ($_ =~ /define/ && $_ =~ /__INTEL_MKL__/) {
				my @tmp = split(/\s+/, $_);
				my $l = @tmp;
				if ($l == 3) {
					$mkl_major_version = $tmp[2];
					$get_it = 1;
				}
			}
		}
		close(READ);

		# if we get it then we break out
		if ($get_it == 1) {
			last;
		}
	}

	# now let's check
	my $ver = $mkl_major_version;
	if ($ver != 11 && $ver != 2017 && $ver != 2018) {
		die "currently we only support mkl in major version of 11, 2017 and 2018\n";
	}

	# return the version
	return $mkl_major_version;
}

#------------------------------------------------------------------------------
# check the tbb library
#------------------------------------------------------------------------------
sub tbb_inquire {

	# check the tbb_root
	if (length($tbb_root) == 0) {
		print "You did not define the TBBROOT ENV\n";
		print "TBB library is a necessary part used in this program\n";
		die "Please download it from https://www.threadingbuildingblocks.org\n";
	}
}

#------------------------------------------------------------------------------
# check the netlib blas+lapack library
#------------------------------------------------------------------------------
sub netlib_inquire {

	# check the netlib_root
	if (length($netlib_root) == 0) {
		print "You did not define the NETLIB_ROOT ENV\n";
		print "netlib library provides BLAS+LAPACK for this program\n";
		die "Please download it from http://www.netlib.org/blas/ and http://www.netlib.org/lapack/\n";
	}
}

#------------------------------------------------------------------------------
# check that whether we need to compile/link against analytical integral deriv
#------------------------------------------------------------------------------
sub integral_deriv_check {

	# now let's check the module
	# we list the modules that will need the 
	# integral derivatives code
	if ($target eq "DGINTS") {
		$with_gints_deriv = 1;
	}
}

#------------------------------------------------------------------------------
# check the external integral library path if it's set
# this function must be called after integral_deriv_check
#------------------------------------------------------------------------------
sub external_integral_lib_check {

	if (not defined $integral_lib_path) {
		$use_external_int_lib = 0;
	}else{

		# issue warning if you choose to use the external integral library
		print "**********************************************************\n";
		print "                         WARNING!!!                       \n";
		print "You have chosen to use external gints engine library\n";
		print "Be caution that the library provided may result in\n";
		print "unsolved reference to compiler libraries. For example,\n";
		print "if you use gints engine library compiled by Intel compiler\n";
		print "but you do not link Intel library in linking, you may have\n";
		print "unsolved reference in the linking stage\n";
		print "**********************************************************\n\n";

		# set the flag to use external integral lib
		$use_external_int_lib = 1;

		# now let's check the path
		# it's already defined
		if (! -d "$integral_lib_path") {
			print "you want to use the pre-compiled integral library such as ERI\n";
			print "however, the library path you pass in is not exist\n";
			die "invalid integral path of $integral_lib_path\n";
		}

		# now let's check several examples
		# such as ERI, nai etc.
		my $i;
		for($i=0; $i<4; $i++) {

			# set name
			my $name = "twobodyoverlap";
			if ($i == 1) {
				$name = "kinetic";
			}
			if ($i == 2) {
				$name = "nai";
			}
			if ($i == 3) {
				$name = "eri";
			}

			# get the lib name
			my $libname = "libgints_engine_"."$name"."_d0.a";
			$libname = "$integral_lib_path"."/"."$libname";
			if (! -e "$libname") {
				print "you want to use the pre-compiled integral library such as $libname\n";
				print "however, this library does not exist\n";
				die "can not find the corresponding integral library\n";
			}
		}

		# now check integral derivatives
		if ($with_gints_deriv == 1) {
			for($i=0; $i<2; $i++) {

				# set name
				my $order = "d1";
				if ($i == 1) {
					$order = "d2";
				}

				# get the lib name
				my $libname = "libgints_engine_eri_"."$order".".a";
				$libname = "$integral_lib_path"."/"."$libname";
				if (! -e "$libname") {
					print "you want to use the pre-compiled integral library such as $libname\n";
					print "however, this library does not exist\n";
					die "can not find the corresponding integral library\n";
				}
			}
		}

		# finally let's write the external integral library
		# path into the env, and we will get it later from cmake
		$ENV{GINTS_ENGINE_LIB_PATH} = "$integral_lib_path";
	}
}

#------------------------------------------------------------------------------
# check the gdm library path if it's set
#------------------------------------------------------------------------------
sub gdm_lib_check {

	if (not defined $gdm_lib_path) {
		$use_gdm = 0;

		# issue warning if you choose to use the external integral library
		print "**********************************************************\n";
		print "                         WARNING!!!                       \n";
		print "GDM function is disabled since you choose not to link\n";
		print "with GDM library\n";
		print "**********************************************************\n\n";

	}else{

		# issue warning if you choose to use the external integral library
		print "**********************************************************\n";
		print "                         WARNING!!!                       \n";
		print "You have chosen to use external GDM library\n";
		print "Be caution that the library provided may result in\n";
		print "unsolved reference to compiler libraries. For example,\n";
		print "if you use the library compiled by Intel compiler\n";
		print "but you do not link Intel library in linking, you may have\n";
		print "unsolved reference in the linking stage\n";
		print "**********************************************************\n\n";

		# set the flag 
		$use_gdm = 1;

		# now let's check the path
		# it's already defined
		if (! -d "$gdm_lib_path") {
			print "you want to use the pre-compiled GDM library\n";
			print "however, the library path you pass in is not exist\n";
			die "invalid integral path of $gdm_lib_path\n";
		}

		# get the lib name
		my $libname = "$gdm_lib_path"."/libgdm.a";
		if (! -e "$libname") {
			print "you want to use the pre-compiled gdm library $libname\n";
			print "however, this library does not exist\n";
			die "can not find the corresponding gdm library\n";
		}

		# finally let's write the external gdm library
		# path into the env, and we will get it later from cmake
		$ENV{GDM_LIB_PATH} = "$gdm_lib_path";
	}
}

#------------------------------------------------------------------------------
# check the boost library
#------------------------------------------------------------------------------
sub boost_inquire {

	# check the boost_root
	if (length($boost_dir) == 0) {
		print "**********************************************************\n";
		print "                         WARNING!!!                       \n";
		print "You did not define the BOOST_ROOT ENV\n";
		print "boost library is a necessary part used in this program\n\n";
		print "You may use the libraries already installed on /use/lib\n";
		print "therefore it does not cause error here, but you may get\n";
		print "a problem later running cmake. If so, please download\n";
		print "the boost library from http://www.boost.org/\n";
		print "**********************************************************\n\n";
	}else{
		# if we get it from the input, we write it into the system env
		my $val = $ENV{BOOST_ROOT};
		if (not defined $val) {
			$ENV{BOOST_ROOT} = "$boost_dir";
		}
	}
}

#------------------------------------------------------------------------------
# check the Intel compiler library
#------------------------------------------------------------------------------
sub intel_inquire {

	# check the intel root path
	if (length($intel_path) == 0) {
		print "You did not define the INTEL_ROOT ENV\n";
		print "You can define it through the command line parameter, see the output of help option\n";
		print "INTEL_ROOT is needed if you use the Intel compiler, we also use the Intel libraries\n";
		print "Please define the Intel root path, like /opt/intel\n";
		die "Intel root path is not defined\n";
	}

	# now let's check that whether the imf library exist
	if (! -f "$intel_path/lib/intel64/libimf.so") {
		print "We are trying to catch the libimf.so in the $intel_path/lib/intel64\n";
		print "However we did not get it, it should be over there unless the path is wrong. Please check it\n";
		die "Intel libimf.so(math library) is not found\n";
	}

	# if we get it from the input, we write it into the system env
	print "INTEL_ROOT is defined as $intel_path\n";
	my $val = $ENV{INTEL_ROOT};
	if (not defined $val) {
		$ENV{INTEL_ROOT} = "$intel_path";
	}
}

#------------------------------------------------------------------------------
# check the gnu gsl library
#------------------------------------------------------------------------------
sub gsl_inquire {

	# check the gsl_root
	if (length($gsl_dir) == 0) {
		print "**********************************************************\n";
		print "                         WARNING!!!                       \n";
		print "You did not define the GSLROOT ENV\n";
		print "GNU GSL library is a necessary part used in this program\n\n";
		print "You may use the libraries already installed on /use/lib\n";
		print "therefore it does not cause error here, but you may get\n";
		print "a problem later running cmake. If so, please download\n";
		print "the gsl library from https://www.gnu.org/software/gsl\n";
		print "**********************************************************\n\n";
	}else{
		# if we get it from the input, we write it into the system env
		my $val = $ENV{GSLROOT};
		if (not defined $val) {
			$ENV{GSLROOT} = "$gsl_dir";
		}
	}
}

#------------------------------------------------------------------------------
# check the MIC coprocessor use
#------------------------------------------------------------------------------
sub check_mic_infor {

	# currently if you want to compile the source code also on MIC card,
	# you can only use the Intel compiler
	if ($cxx_compiler ne "icpc") {
		print "if you want to compile code also for Intel Phi co-processor, you need to use\n";
		print "Intel C++ compiler icpc\n";
		die "please set the cxx_compiler to intel\n";
	}
}

#------------------------------------------------------------------------------
# print out the user settings
#------------------------------------------------------------------------------
sub print_settings {

	print "**********************************************************\n";
	print "*                  USER OPTIONS LIST                     *\n";
	print "**********************************************************\n";
	print "Source code directory is  : ", $source_dir, "\n";
	print "Building directory is     : ", $build_dir, "\n";
	print "Building type is          : ", $build_type, "\n";
	print "BLAS/LAPACK vendor is     : ", $blas_vendor, "\n";
	if ($blas_vendor eq "mkl") {
		print "MKL's root is             : ", $mkl_root, "\n";
		print "MKL's major version is    : ", $mkl_major_version, "\n";
		if ($use_sequential_mkl == 1) {
			print "MKL is using the sequential type library", "\n";
		}else{
			print "MKL is using the parallel type library", "\n";
		}
	}
	print "TBB's root is             : ", $tbb_root, "\n";
	if (length($boost_dir) != 0) {
		print "BOOST dir is              : ", $boost_dir, "\n";
	}
	if (length($intel_path) != 0) {
		print "INTEL dir is              : ", $intel_path, "\n";
	}
	if (length($gsl_dir) != 0) {
		print "GNU GSL dir is            : ", $gsl_dir, "\n";
	}
	print "target module is          : ", $target, "\n";
	if ($with_gints_deriv == 0) {
		print "with gints derivatives    : ", "NO", "\n";
	}else{
		print "with gints derivatives    : ", "YES", "\n";
	}
	if ($use_external_int_lib == 0) {
		print "use external integral lib?: ", "NO", "\n";
	}else{
		print "use external integral lib?: ", "YES", "\n";
		print "external integral lib path: ", "$integral_lib_path", "\n";
	}
	if ($use_gdm == 0) {
		print "use external GDM lib?     : ", "NO", "\n";
	}else{
		print "use external GDM lib?     : ", "YES", "\n";
		print "external GDM lib path     : ", "$gdm_lib_path", "\n";
	}
	print "\n\n\n";
}

#------------------------------------------------------------------------------
# prepare work for running cmake
#------------------------------------------------------------------------------
sub prepare_cmake {

   # check that whether we have CMakeLists.txt in current folder
	my $cmake_file = "$source_dir"."/"."CMakeLists.txt";
	if (! -e "$cmake_file" ) {
		print $cmake_file;
		die "In current folder we can not find the CMakeLists.txt to run CMAKE\n";
	}

	# if use mkl library, we need to know mkl version first
	# also check all other possible blas+lapack libraries choice
	if ($blas_vendor eq "mkl") {
		$mkl_major_version = mkl_inquire(); 
	}elsif ($blas_vendor eq "netlib") {
		netlib_inquire(); 
	}

	# checking TBB
	tbb_inquire();

	# checking boost
	boost_inquire();

	# if we use Intel compiler(C++), we check 
	# the Intel library usage
	if ($cxx_compiler eq "icpc") {
		intel_inquire();
	}

	# checking gnu gsl
	# disabled here
	gsl_inquire();

	# check that wether we go with gints derivatives code?
	integral_deriv_check(); 

	# check the external integral library information if it's set
	external_integral_lib_check(); 

	# check the GDM library
	gdm_lib_check(); 

   # prepare the new build dir
	# this is the last step
	if (-d "$build_dir") {
		print "**********************************************************\n";
		print "                         WARNING!!!                       \n";
		print "The building directory exists, we will remove it\n";
		print "**********************************************************\n\n";
		system("rm -rf $build_dir");
	}
	mkdir($build_dir, 0777) || die "can not creat the build dir\n";
}

#------------------------------------------------------------------------------
# pass these things to cmake to finish building
#------------------------------------------------------------------------------
sub run_cmake {

	# enter the build dir
	chdir($build_dir);
	my $cmake_command = "cmake";

   # add in compiler setting
	$cmake_command = "$cmake_command"." -D CMAKE_CXX_COMPILER=$cxx_compiler";
	$cmake_command = "$cmake_command"." -D CMAKE_Fortran_COMPILER=$for_compiler";
	$cmake_command = "$cmake_command"." -D COMPILER_WARNING_OFF=$disable_compiler_warning";

   # add in the building type
	$cmake_command = "$cmake_command"." -D CMAKE_BUILD_TYPE=$build_type";

   # add in the BLAS/LAPACK vendor information
	$cmake_command = "$cmake_command"." -D BLAS_VENDOR=$blas_vendor";
	if (length($mkl_root) > 0 && $blas_vendor eq "mkl") {
		$cmake_command = "$cmake_command"." -D MKLROOT=$mkl_root";
		$cmake_command = "$cmake_command"." -D MKL_MAJOR_VERSION=$mkl_major_version";
		$cmake_command = "$cmake_command"." -D USE_SEQUENTIAL_MKL_LIB=$use_sequential_mkl";
	}

	# add in target module information
	$cmake_command = "$cmake_command"." -D TARGET_MODULE=$target";

	# check whether we go with gints derivatives code?
	$cmake_command = "$cmake_command"." -D WITH_GINTS_DERIV=$with_gints_deriv";

	# check whether we go with external gints engine library?
	$cmake_command = "$cmake_command"." -D WITH_EXTERNAL_GINTS_ENGINE_LIB=$use_external_int_lib";

	# check whether we go with GDM library?
	$cmake_command = "$cmake_command"." -D WITH_GDM_LIB=$use_gdm";

	# are you going to use MIC card?
	if ($use_mic == 1) {
		$cmake_command = "$cmake_command"." -D WITH_INTEL_PHI_COPROC=1";
	}else{
		$cmake_command = "$cmake_command"." -D WITH_INTEL_PHI_COPROC=0";
	}

   #finally, add in source dir
	$cmake_command = "$cmake_command"." $source_dir";
	#print "$cmake_command\n";
	my $status = system("$cmake_command");
	if ($status != 0) {
		print "The return status from system command in perl is not 0\n";
		print "Please check that whether you have install the cmake command??\n";
		print "or any libraries is not properly installed etc. so cmake returns error information??\n";
	}
}

#################################
# @@@@ main script for run      #
#################################

#
# now let's parse the input choice
#
para_parsing();

#
# redirect the output to the file of configure.log
#
if (-e "configure.log") {
	system ("rm -rf configure.log");
}
open (STDOUT, "| tee -ai configure.log");

#
# put the full command line to the output
#
my $commandline = join " ", $0, @ARGV;
print "**********************************************************\n";
print "*              FULL COMMAND LINE OPTIONS                 *\n";
print "**********************************************************\n";
print "perl $commandline", "\n";
print "\n";

#
# now do prepare work before calling cmake
#
prepare_cmake();

#
# prepare the user's setting
#
print_settings(); 

#
# now running cmake
#
run_cmake();



