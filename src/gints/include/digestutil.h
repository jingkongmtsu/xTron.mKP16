/**
 * \file   digestutil.h
 * \brief  providing some utility functions for digestion
 * \author Fenglai Liu
 */
#ifndef DIGESTUTIL_H
#define DIGESTUTIL_H
#include "libgen.h"

namespace matrix {
	class Mtrx;
}

namespace shell {
	class MolShell;
}

namespace digestutil {

	using namespace shell;
	using namespace matrix;

	///
	/// in the case of symmetrical digestion, we note that the result block could 
	/// be in upper part, or lower part; this is same for both J and K situation.
	/// now we will update the result between upper and lower part for the result
	/// matrix
	///
	/// on the other hand, the updation(add transpose) will make the diagonal
	/// block 2 times larger, we need to cancel this side effect
	///
	void postProcessJKMatrix(const MolShell& s, Mtrx& K);

	/**
	 * this function performs basis set normalization scaling 
	 * for the given integral I(nbra1,nbra2,nket1,nket2)
	 * such function could be used in forming the integral boundary
	 * like Cauchy–Schwarz screening for ERI
	 *
	 * \param LBra1Min, LBra1Max : angular momentum for bra1 dimension
	 * \param LBra2Min, LBra2Max : angular momentum for bra2 dimension
	 * \param LKet1Min, LKet1Max : angular momentum for ket1 dimension
	 * \param LKet2Min, LKet2Max : angular momentum for ket2 dimension
	 * \param bra1Scale          : basis set normalization vector for bra1 dimension
	 * \param bra2Scale          : basis set normalization vector for bra2 dimension
	 * \param ket1Scale          : basis set normalization vector for ket1 dimension
	 * \param ket2Scale          : basis set normalization vector for ket2 dimension
	 * \return        I          : four index integral array
	 */
	void scale4BodyInts(const UInt& LBra1Min, const UInt& LBra1Max, 
			const UInt& LBra2Min, const UInt& LBra2Max, 
			const UInt& LKet1Min, const UInt& LKet1Max, 
			const UInt& LKet2Min, const UInt& LKet2Max, 
			const Double* bra1Scale, const Double* bra2Scale, 
			const Double* ket1Scale, const Double* ket2Scale,
			Double* I); 
	
};

#endif
