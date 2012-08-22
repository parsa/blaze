//=================================================================================================
/*!
//  \file blaze/math/problems/MixedLCP.h
//  \brief A data structure for mixed linear complementarity problems
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

#ifndef _BLAZE_MATH_PROBLEMS_MIXEDLCP_H_
#define _BLAZE_MATH_PROBLEMS_MIXEDLCP_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Functions.h>
#include <blaze/math/Infinity.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief A mixed linear complementarity problem (MLCP) data structure.
// \ingroup math
//
// The LCP class represent a mixed linear complementarity problem of the form

                          \f[
                          \left(\begin{array}{*{2}{c}}
                          A_{11} & A_{12} \\
                          A_{21} & A_{22} \\
                          \end{array}\right) \cdot

                          \left(\begin{array}{*{1}{c}}
                          x_{1} \\
                          x_{2} \\
                          \end{array}\right) +

                          \left(\begin{array}{*{1}{c}}
                          b_{1} \\
                          b_{2} \\
                          \end{array}\right) \leq

                          \left(\begin{array}{*{1}{c}}
                          0 \\
                          0 \\
                          \end{array}\right) \quad\perp\quad

                          \left(\begin{array}{*{1}{c}}
                          x_{1} \\
                          x_{2} \\
                          \end{array}\right) \geq

                          \left(\begin{array}{*{1}{c}}
                          0 \\
                          0 \\
                          \end{array}\right)
                          \f]\n
*/
struct MixedLCP
{
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t size       ()               const;
   inline size_t equations  ()               const;
   inline size_t constraints()               const;
   inline void   project    ( size_t index );
   inline real   lbound     ( size_t index ) const;
   inline real   ubound     ( size_t index ) const;
   inline real   residual   ( size_t index ) const;
   inline real   residual   ()               const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   CMatMxN A11_;  //!< The upper left part of the system matrix \f$ A_{11} \f$.
   CMatMxN A12_;  //!< The upper right part of the system matrix \f$ A_{12} \f$.
   CMatMxN A21_;  //!< The lower left part of the system matrix \f$ A_{21} \f$.
   CMatMxN A22_;  //!< The lower right part of the system matrix \f$ A_{22} \f$.
   VecN    b1_;   //!< The upper part of the right-hand side vector \f$ b_{1} \f$.
   VecN    b2_;   //!< The lower part of the right-hand side vector \f$ b_{2} \f$.
   VecN    x1_;   //!< The upper part of the vector of unknowns \f$ x_{1} \f$.
   VecN    x2_;   //!< The lower part of the vector of unknowns \f$ x_{2} \f$.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the size of the mixed linear complementarity problem.
//
// \return The actual size of the MLCP.
*/
inline size_t MixedLCP::size() const
{
   return x1_.size() + x2_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of equations of the MLCP.
//
// \return The number of equations.
*/
inline size_t MixedLCP::equations() const
{
   return x1_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of constraints of the MLCP.
//
// \return The number of constraints.
*/
inline size_t MixedLCP::constraints() const
{
   return x2_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Projects the unknown at the given index on the solution range.
//
// \param index Access index. The index has to be in the range \f$ [0..size) \f$.
// \return void
*/
inline void MixedLCP::project( size_t index )
{
   if( index < x1_.size() )
      return;

   index -= x1_.size();
   x2_[index] = max( real(0), x2_[index] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the lower bound of the unknown at the given index.
//
// \param index Access index. The index has to be in the range \f$ [0..size) \f$.
// \return void
*/
inline real MixedLCP::lbound( size_t index ) const
{
   if( index < x1_.size() )
      return -inf;
   else return real(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the upper bound of the unknown at the given index.
//
// \param index Access index. The index has to be in the range \f$ [0..size) \f$.
// \return void
*/
inline real MixedLCP::ubound( size_t /*index*/ ) const
{
   return inf;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the residual of the unknown at the given index.
//
// \param index Access index. The index has to be in the range \f$ [0..size) \f$.
// \return The residual at index \a index.
*/
inline real MixedLCP::residual( size_t index ) const
{
   // Calculating the LSE residual by Ax+b
   if( index < x1_.size() )
      return ( A11_ * x1_ )[index] + ( A12_, x2_ )[index];

   index -= x1_.size();

   // Calculating the LCP residual by max( x - xmax, min( x - xmin, Ax+b ) )
   return min( x2_[index],
               ( A21_ * x1_ )[index] + ( A22_ * x2_ )[index] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculates the maximum norm of the residual of the mixed LCP.
//
// \return The maximum norm of the global residual of the MLCP.
*/
inline real MixedLCP::residual() const
{
   real rmax( 0 );

   for( size_t i=0; i<size(); ++i )
      rmax = max( rmax, std::fabs( residual( i ) ) );

   return rmax;
}
//*************************************************************************************************

} // namespace blaze

#endif
