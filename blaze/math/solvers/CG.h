//=================================================================================================
/*!
//  \file blaze/math/solvers/CG.h
//  \brief Header file for the conjugate gradient solver
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

#ifndef _BLAZE_MATH_SOLVERS_CG_H_
#define _BLAZE_MATH_SOLVERS_CG_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/CompressedMatrix.h>
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
/*!\brief A conjugate gradient solver.
// \ingroup lse_solvers
//
// TODO: description
// TODO: Problem formulation: \f$ A \cdot x + b = 0 \f$ !!
*/
class CG : public Solver
{
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit CG();
   //@}
   //**********************************************************************************************

   //**Solver functions****************************************************************************
   /*!\name Solver functions */
   //@{
   bool solve( LSE& lse );
   bool solve( const CMatMxN& A, const VecN& b, VecN& x );
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VecN r_;  //!< TODO
   VecN d_;  //!< TODO
   VecN h_;  //!< TODO
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
//
// TODO: description
// TODO: Problem formulation: \f$ A \cdot x + b = 0 \f$ !!
*/
inline bool CG::solve( LSE& lse ) {
   return solve( lse.A_, lse.b_, lse.x_ );
}
//*************************************************************************************************

} // namespace blaze

#endif
