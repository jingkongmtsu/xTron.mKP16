/**
 * \file   scalablevec.h
 * 
 * All of vector form of data uses Intel TBB's scalable allocator,
 * and you can also choose the use other vectors, too. For STL one,
 * We did not test the STL allocator in the multi-threading env.
 * So be careful to use it!!
 * 
 * For TBB, we use the allocator since version 4.2.
 *
 * STL allocator is only "THREAD SAFE" in multi-threading env. since C++11
 * according to what C++11 manual says. So if you use with_std_alloc, 
 * be careful to check the version of allocator!!!
 *
 * In C++11, std::allocator is thread safe. From its definition:
 * 20.6.9.1/6: Remark: the storage is obtained by calling ::operator new(std::size_t)
 * and from the definition of ::operator new:
 * 18.6.1.4: The library versions of operator new and operator delete, 
 * user replacement versions of global operator new and operator delete, 
 * and the C standard library functions calloc, malloc, realloc, and free 
 * shall not introduce data races (1.10) as a result of concurrent calls 
 * from different threads.
 *
 * The Intel TBB scalable allocator is able to allocate memory
 * in scalable way for multi-threading env. For single threading
 * env. scalable allocator behaves in the similar way of STL one.
 * For more information, please refer to the TBB's document.
 *
 * Comparing with STL allocator, for large memory chunk, there's 
 * no memory overhead for TBB allocator; however, it's 3-4 times
 * slower than the STL one(so as the deallocation process). For 
 * some memory allocation between large and small, it's said that
 * the memory overhead may be significant. But so far we do not
 * have data on hand to show it.
 *
 * see the post of https://software.intel.com/en-us/forums/topic/516219
 * for discussion on the difference between STL allocator 
 * and TBB one, as we showed above.
 *
 * \author Fenglai Liu and Jing Kong
 *
 */

#ifndef SCALABLEVEC_H
#define SCALABLEVEC_H
#include "libgen.h"
#include<vector>
#include "tbb/scalable_allocator.h"

//
// Double type vector is used everywhere
// however, UInt vector and Int vector is useful in some particular space,
// for example; may be in DFT part
// so we also declare them
//
typedef std::vector<Double,tbb::scalable_allocator<Double> >   DoubleVec;
typedef std::vector<UInt,tbb::scalable_allocator<UInt> >       UIntVec;
typedef std::vector<Int,tbb::scalable_allocator<Int> >         IntVec;

//
// some times we may need the exact double precision vectors for calculation
// for example, in the functional derivatives calculation of DFT part;
// the fortran code defines all of inputs as REAL*8
// this is necessary for accurate calculation
// so we add the case here below
//
typedef std::vector<double,tbb::scalable_allocator<double> >   DoublePrecVec;

#endif

