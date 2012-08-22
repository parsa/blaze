//=================================================================================================
/*!
//  \file blaze/math/solvers/GaussianElimination.h
//  \brief Header file for the GaussianElimination class
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_SOLVERS_GAUSSIANELIMINATION_H_
#define _BLAZE_MATH_SOLVERS_GAUSSIANELIMINATION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/problems/LSE.h>
#include <blaze/math/solvers/Solver.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of the Gaussian elimination direct linear system solver.
// \ingroup lse_solvers
//
// TODO: description
// TODO: Problem formulation: \f$ A \cdot x + b = 0 \f$ !!
*/
class GaussianElimination : public Solver
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit GaussianElimination();
   //@}
   //**********************************************************************************************

   //**Solver functions****************************************************************************
   /*!\name Solver functions */
   //@{
   inline bool solve( LSE& lse );
          bool solve( const CMatMxN& A, const VecN& b, VecN& x );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MatMxN A_;  //!< TODO
   VecN   b_;  //!< TODO
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  SOLVER FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief TODO
//
// \param lse TODO
// \return TODO
// \exception std::invalid_argument Invalid matrix size.
// \exception std::invalid_argument Invalid right-hand side vector size.
//
// TODO: description
// TODO: Problem formulation: \f$ A \cdot x + b = 0 \f$ !!
*/
inline bool GaussianElimination::solve( LSE& lse )
{
   return solve( lse.A_, lse.b_, lse.x_ );
}
//*************************************************************************************************

} // namespace blaze

#endif
