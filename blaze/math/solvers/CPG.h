//=================================================================================================
/*!
//  \file blaze/math/solvers/CPG.h
//  \brief Implementation of the conjugate projected gradient algorithm
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

#ifndef _BLAZE_MATH_SOLVERS_CPG_H_
#define _BLAZE_MATH_SOLVERS_CPG_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/problems/BoxLCP.h>
#include <blaze/math/problems/ContactLCP.h>
#include <blaze/math/problems/LCP.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/solvers/Solver.h>
#include <blaze/util/Assert.h>
#include <blaze/util/ColorMacros.h>
#include <blaze/util/logging/DebugSection.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of the conjugate projected gradient algorithm.
// \ingroup complementarity_solvers
//
// TODO: description of the CPG solver
// TODO: capabilities of the CPG solver (which LCP problems, etc)
// TODO: known issues of the CPG solver
*/
class CPG : public Solver
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit CPG();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename CP > bool solve( CP& cp );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VecN                r_;         //!< TODO
   VecN                w_;         //!< TODO
   VecN                p_;         //!< TODO
   DynamicVector<int>  activity_;  //!< TODO
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
/*!\brief TODO
//
// \param cp TODO
// \return void
//
// TODO
*/
template< typename CP >  // Type of the complementarity problem
bool CPG::solve( CP& cp )
{
   const size_t n( cp.size() );
   const CMatMxN& A( cp.A_ );
   const VecN&    b( cp.b_ );

   bool converged( false );
   VecN& x( cp.x_ );
   size_t activeSetChanges( 0 );
   real alpha( 0 ), alpha_nom( 0 ), alpha_denom( 1 );
   real beta( 0 ),  beta_nom( 0 ),  beta_denom( 0 );
   real tmp( 0 );

   BLAZE_INTERNAL_ASSERT( isSymmetric( A ), "The CPG solver requires that the system matrix is symmetric" );

   // Allocating helper data
   r_.resize( n, false );
   w_.resize( n, false );
   p_.resize( n, false );
   activity_.resize( n, false );

   // Determining activity and project initial solution to feasible region
   for( size_t i=0; i<n; ++i ) {
      if( x[i] <= cp.lbound( i ) ) {
         x[i] = cp.lbound( i );
         ++activeSetChanges;
         activity_[i] = -1;
      }
      else if( x[i] >= cp.ubound( i ) ) {
         x[i] = cp.ubound( i );
         ++activeSetChanges;
         activity_[i] = 1;
      }
      else {
         activity_[i] = 0;
      }
   }

   // Computing the initial residual
   lastPrecision_ = cp.residual();
   if( lastPrecision_ < threshold_ )
      converged = true;

   // Choosing the initial values such that the descent direction conjugation process is disabled
   p_ = real(0);
   w_ = real(0);

   size_t it( 0 );
   for( ; !converged && it < maxIterations_; ++it )
   {
      // Computing the steepest descent direction
      r_ = -( A*x + b );

      // Projecting the gradient and the previous descent direction
      beta_nom = real(0);
      beta_denom = alpha_denom;

      for( size_t i=0; i<n; ++i ) {
         tmp = r_[i];
         if( activity_[i] == -1 ) {
            tmp   = max( tmp,   0 );
            p_[i] = max( p_[i], 0 );
         }
         else if( activity_[i] == 1 ) {
            tmp   = min( tmp,   0 );
            p_[i] = min( p_[i], 0 );
         }
         beta_nom += w_[i] * tmp;
         w_[i] = tmp;
      }

      if( beta_denom == 0 ) {
         // No conjugation can be performed, fallback to steepest descent
         beta = 0;
      }
      else {
         beta = -beta_nom / beta_denom;
      }

      BLAZE_INTERNAL_ASSERT( !isnan( beta ), "Conjugation coefficient is nan" );

      // Choosing the next descent direction conjugated to all previous directions
      p_ = w_ + beta * p_;

      // Finding the minimum along the descent direction p
      alpha_nom   = trans(r_) * p_;
      alpha_denom = trans(p_) * A * p_;

      if( alpha_denom == 0 )
         // In case p^T A p is zero, no reduction of the objective function can be obtained
         // in direction p no matter which alpha we choose
         alpha = 0;
      else
         alpha = alpha_nom / alpha_denom;

      if( alpha == 0 ) {
         if( beta == 0 ) {
            // p is the steepest descent direction since beta = 0 but we still cannot
            // make any progress along p => minimum
            break;
         }
         else {
            // Retry with steepest descent direction
            continue;
         }
      }

      // Descending along p and projecting
      activeSetChanges = 0;

      for( size_t i=0; i<n; ++i ) {
         if( activity_[i] != 0 && p_[i] == real(0) ) {
            // In case the bounds are depending on the unknowns this ensures that
            // the unknowns stay at the bounds
            if( activity_[i] == -1 )
               x[i] = cp.lbound( i );
            else
               x[i] = cp.ubound( i );
         }
         else {
            x[i] += alpha * p_[i];

            if( x[i] <= cp.lbound( i ) ) {
               x[i] = cp.lbound( i );
               if( activity_[i] != -1 )
                  ++activeSetChanges;
               activity_[i] = -1;
            }
            else if( x[i] >= cp.ubound( i ) ) {
               x[i] = cp.ubound( i );
               if( activity_[i] != +1 )
                  ++activeSetChanges;
               activity_[i] = +1;
            }
            else {
               if( activity_[i] != 0 )
                  ++activeSetChanges;
               activity_[i] = 0;
            }
         }
      }

      // Computing the residual (TODO we should improve this)
      lastPrecision_ = cp.residual();
      if( lastPrecision_ < threshold_ )
         converged = true;
   }

   BLAZE_LOG_DEBUG_SECTION( log ) {
      if( converged )
         log << "      Solved the quadratic program in " << it << " CPG iterations.";
      else
         log << BLAZE_YELLOW << "      WARNING: Did not solve the quadratic program within accuracy. (" << lastPrecision_ << ")" << BLAZE_OLDCOLOR;
   }

   lastIterations_ = it;

   return converged;
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPLICIT TEMPLATE INSTANTIATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if !defined(_MSC_VER)
extern template bool CPG::solve<LCP>( LCP& );
extern template bool CPG::solve<BoxLCP>( BoxLCP& );
extern template bool CPG::solve<ContactLCP>( ContactLCP& );
#endif
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
