/**
 * \file   libgen.h
 * \brief  general data used in the whole program
 * \author Fenglai Liu and Jing Kong
 *
 */

#ifndef LIBGEN_H
#define LIBGEN_H
#include<cstddef>

/**
 * re-define the types here
 * we do not use unsigned type for int, since the common int
 * should be good enough
 */
//
// integer
//
typedef size_t        UInt;
typedef int           Int;
typedef long long     LInt;

//
// floating 
//
#ifdef WITH_SINGLE_PRECISION
typedef float         Double;
#else
typedef double        Double;
#endif
typedef long double   LDouble;


/**
 * general constant used in the program
 */
#define  MINUS_HALF   -0.5E0
#define  MINUS_ONE    -1.0E0
#define  MINUS_TWO    -2.0E0
#define  MINUS_THREE  -3.0E0
#define  MINUS_FOUR   -4.0E0
#define  MINUS_FIVE   -5.0E0
#define  MINUS_SIX    -6.0E0
#define  MINUS_SEVEN  -7.0E0
#define  MINUS_EIGHT  -8.0E0
#define  MINUS_NINE   -9.0E0
#define  MINUS_TEN    -10.0E0
#define  ZERO          0.0E0
#define  ONE           1.0E0
#define  TWO           2.0E0
#define  THREE         3.0E0
#define  FOUR          4.0E0
#define  FIVE          5.0E0
#define  SIX           6.0E0
#define  SEVEN         7.0E0
#define  EIGHT         8.0E0
#define  NINE          9.0E0
#define  TEN           10.0E0
#define  HALF          0.5E0

///
/// define a large number used in the program
/// this should be usually good enough
/// 
/// on the other hand, this number does not exceed the 
/// range of float limit, therefore it's still safe
/// to perform arithmetic operations
///
/// so we choose the float point limit to define it
///
#define LARGE_FLOAT_NUMBER    1.0E+38

///
/// physical and math constants
///
#define  BOHR_TO_ANGSTROM  0.5291772109217E0
#define  ANGSTROM_TO_BOHR  1.8897261245590E0
#ifndef PI
#define  PI                3.14159265358979323846E0
#endif
#define  DEGREE_TO_RADIAN  0.01745329251994329576E0


///
/// constant used to measure the difference between two Double value
/// since the double accuracy is around 15-16, so we set it as 1.0E-14 
/// therefore we take its lower bound
///
#define  THRESHOLD_MATH  1.0E-14

///
/// in general, let's define the maximum number of 
/// molecule sections permitted in the program
///
///
#define MAX_MOLECULE_SECTION_NUMBER  10000

///
/// now this section defines the processor type
/// currently we can load the work to three types of processors: CPU, INTEL_Phi and GPGPU
///
#define  CPU          1
#define  INTEL_PHI    2
#define  GPGPU        3

///
/// define the type of shell, used universally in the whole program
///  
/// TYPE_NORM refers to the basis set data in "normal" order, which is the original definition
/// of basis set data(could be all spherical, like ccpvdz, or all Cartesian, like 6-31G)
///
/// TYPE_CART refers to the basis set data all in "Cartesian" order. No matter how the original
/// basis set data is defined, here the basis set are all as Cartesian type
///
#define  TYPE_NORM    0
#define  TYPE_CART    1

///
/// here defines the memory allocation/free status for using offload module 
///
#define ALLOC   alloc_if(1)
#define REUSE   alloc_if(0)
#define FREE    free_if(1)
#define RETAIN  free_if(0)

///
/// for using the SIMD operations, we have a requirement to align the data to the cache line
/// current cache line for CPU is 64 bit (Intel Xeon Phi and scalable processors etc.)
///
#define CPU_SIMD_ALIGNMENT         64

///
/// for the SIMD operations, let's define the basic loop count. currently for utilizing the 
/// AVX-512 instruction set, the foundamental loop count is 8
///
#define SIMD_LOOP_COUNT            8

#endif

