/**
 * \file    scfconvsol.h
 * \brief   solutions class provided to solve all of SCF convergence problems
 * \author  Fenglai Liu 
 */
#ifndef SCFCONVSOL_H
#define SCFCONVSOL_H
#include "libgen.h"
#include "scfmacro.h"

namespace scfparam {
	class SCFParam;
}

namespace globalinfor {
	class GlobalInfor;
}

namespace scfconvsol {

	using namespace scfmacro;
	using namespace scfparam;
	using namespace globalinfor;

	/**
	 *
	 * for the SCF convergence progress, each SCF iteration we may have "bad status"
	 * (known as problem for scf, defined in scfmacro.h), and for each of these
	 * problem we may have a couple of solutions for solving it until good status coming 
	 * out(the solutions are also defined in scfmacro.h).
	 *
	 * This class is to arrange the possible solutions(user can define them) for all
	 * of problems encounter in SCF.
	 *
	 * The idea for "problem-solution" pick up is designed only for one given SCF iteration.
	 * Here we does not assume that the problem for this SCF iteration may affect the 
	 * following SCF iterations. 
	 *
	 * This may be the case. If we have the SCF convergence problem and we solve it
	 * in the current SCF iteration, such problem will not come back again. However,
	 * we may also see the case that the same problem (or similar problem) will occur
	 * again and again for the SCF iterations.
	 *
	 * So in all of cases, for each SCF iteration we just treat the given problem as 
	 * a new one and we explore the solution space for an answer. Therefore, you need
	 * to reset the index in the end of SCF interation.
	 */
	class SCFConvSol {

		private:

			bool printSCFConvSol;                   ///< whether to print out the content of this class?
			UInt spaceVecRepeatIndex;               ///< solution index for spaceVecRepeat
			UInt currentScfIndexNotInsideIndex;     ///< solution index for currentScfIndexNotInside
			UInt scfOscillationIndex;               ///< solution index for scfOscillation

			///
			/// solutions provided for soliing the problem that the DIIS/EDIIS etc. space
			/// vector repeating by comparing with the history archive
			///
			UInt spaceVecRepeat[MAXSOL];      
			
			///
			/// solutions provided by user for solving the problem that
			/// the current scf index is not inside the result scf vector space
			///
			UInt currentScfIndexNotInside[MAXSOL];   

			///
			/// solution provided to solve the oscillation problem
			///
			UInt scfOscillation[MAXSOL];

		public:

			///
			/// contructor: to initialize the solution arrays from user defined input file
			///
			SCFConvSol(const GlobalInfor& globInfor, const SCFParam& par);

			///
			/// destructor
			///
			~SCFConvSol() { };

			///
			/// debug printing function
			///
			void print() const;

			///
			/// reset the referring index of the solution arrays
			/// use it when the scf convergence algorithm is applied,
			/// this is for the next round of SCF iteration
			///
			void resetSolArrayIndex() {
				spaceVecRepeatIndex           = 0;  
				currentScfIndexNotInsideIndex = 0; 
				scfOscillationIndex           = 0;      
			};

			///
			/// for the given problem identified, we will return the solution
			/// from corresponding solution array
			///
			/// we use the index value (defined for each array) to target
			/// the solution value, each time we refer the solution value
			/// we will increase the correspinding index value until we
			/// we run out of solutions inside the array
			///
			/// if in this case, we will just return the null solution - which
			/// means no solution anymore
			///
			/// we note that such problem-solution pick up should be for a 
			/// given SCF iteration. If  the iteration is done (in the end of scfconv),
			/// we should reset the solution index. 
			///
			UInt getCurrentSolution(UInt problem); 

			///
			/// whether print out content of this class?
			///
			bool doPrint() const { return printSCFConvSol; };
	};

}

#endif

