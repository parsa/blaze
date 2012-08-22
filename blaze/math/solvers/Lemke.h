//=================================================================================================
/*!
//  \file blaze/math/solvers/Lemke.h
//  \brief The Lemke pivoting algorithm for solving LCPs.
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

#ifndef _BLAZE_MATH_SOLVERS_LEMKE_H_
#define _BLAZE_MATH_SOLVERS_LEMKE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iosfwd>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/problems/LCP.h>
#include <blaze/math/solvers/Solver.h>
#include <blaze/system/Precision.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The Lemke pivoting algorithm for solving LCPs.
// \ingroup complementarity_solvers
//
// TODO
*/
class Lemke : public Solver
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit Lemke();
   //@}
   //**********************************************************************************************

   //**Solver functions****************************************************************************
   /*!\name Solver functions */
   //@{
   bool solve( LCP& lcp );
   bool solve( LCP& lcp, const VecN& d );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void printTableau( std::ostream& os ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   bool isComponentwiseNonnegative( const VecN& v ) const;
   bool isComponentwisePositive   ( const VecN& v ) const;
   bool isLexicographicallyLess   ( size_t i1, real f1, size_t i2, real f2 ) const;
   bool isLexicographicallyGreater( size_t i1, real f1, size_t i2, real f2 ) const;
   void pivot                     ( size_t block, size_t drive );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   DynamicVector<ptrdiff_t> basics_;     //!< TODO
   DynamicVector<ptrdiff_t> nonbasics_;  //!< TODO
   MatMxN                   M_;          //!< TODO
   MatMxN                   Q_;          //!< TODO
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************

} // namespace blaze

#endif
