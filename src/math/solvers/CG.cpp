//*************************************************************************************************
/*!
//  \file src/math/solvers/CG.cpp
//  \brief Source file for the conjugate gradient solver
//
//  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
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

   if( !isSymmetric( A ) )
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
