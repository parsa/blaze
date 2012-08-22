//*************************************************************************************************
/*!
//  \file src/math/solvers/CG.cpp
//  \brief Source file for the conjugate gradient solver
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
//*************************************************************************************************


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <stdexcept>
#include <blaze/math/Functions.h>
#include <blaze/math/solvers/CG.h>
#include <blaze/util/Assert.h>
#include <blaze/util/ColorMacros.h>
#include <blaze/util/logging/DebugSection.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for the conjugate gradient solver.
*/
CG::CG()
   : r_()  // TODO
   , d_()  // TODO
   , h_()  // TODO
{}
//*************************************************************************************************




//=================================================================================================
//
//  SOLVER FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief TODO
//
// \param A TODO
// \param b TODO
// \param x TODO
// \return void
// \exception std::invalid_argument System matrix is not square.
// \exception std::invalid_argument System matrix is not symmetric.
// \exception std::invalid_argument Invalid right-hand side vector size.
//
// TODO: description
// TODO: Problem formulation: \f$ A \cdot x + b = 0 \f$ !!
*/
bool CG::solve( const CMatMxN& A, const VecN& b, VecN& x )
{
   const size_t n( b.size() );
   bool converged( false );
   real alpha, beta, delta;

   if( A.rows() != A.columns() )
      throw std::invalid_argument( "System matrix is not square" );

   if( !A.isSymmetric() )
      throw std::invalid_argument( "System matrix is not symmetric" );

   if( A.rows() != b.size() )
      throw std::invalid_argument( "Invalid right-hand side vector size" );

   // Allocating helper data
   r_.resize( n, false );
   d_.resize( n, false );
   h_.resize( n, false );

   // Preparing the vector of unknowns
   x.resize( n, false );
   x.reset();

   // Computing the initial residual
   r_ = A * x + b;

   // Initial convergence test
   lastPrecision_ = 0;
   for( size_t i=0; i<n; ++i ) {
      lastPrecision_ = max( lastPrecision_, std::fabs( r_[i] ) );
   }

   if( lastPrecision_ < threshold_ )
      converged = true;

   delta = trans(r_) * r_;

   d_ = -r_;

   // Performing the CG iterations
   size_t it( 0 );

   for( ; !converged && it<maxIterations_; ++it )
   {
      h_ = A * d_;

      alpha = delta / ( trans(d_) * h_ );

      x  += alpha * d_;
      r_ += alpha * h_;

      lastPrecision_ = 0;
      for( size_t i=0; i<n; ++i ) {
         lastPrecision_ = max( lastPrecision_, std::fabs( r_[i] ) );
      }

      if( lastPrecision_ < threshold_ ) {
         converged = true;
         break;
      }

      beta = trans(r_) * r_;

      d_ = ( beta / delta ) * d_ - r_;

      delta = beta;
   }

   BLAZE_LOG_DEBUG_SECTION( log ) {
      if( converged )
         log << "      Solved the linear system in " << it << " CG iterations.";
      else
         log << BLAZE_YELLOW << "      WARNING: Did not solve the linear system within accuracy. (" << lastPrecision_ << ")" << BLAZE_OLDCOLOR;
   }

   lastIterations_ = it;

   return converged;
}
//*************************************************************************************************

} // namespace blaze
