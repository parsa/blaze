//*************************************************************************************************
/*!
//  \file src/math/solvers/GaussianElimination.cpp
//  \brief Source file for the GaussianElimination class
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

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <blaze/math/Accuracy.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/solvers/GaussianElimination.h>
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
/*!\brief The default constructor for the GaussianElimination class.
*/
GaussianElimination::GaussianElimination()
   : A_()  // TODO
   , b_()  // TODO
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
// \return TODO
// \exception std::invalid_argument Invalid matrix size.
// \exception std::invalid_argument Invalid right-hand side vector size.
//
// TODO: description
// TODO: problem formulation \f$ A \cdot x + b = 0 \f$ !!
*/
bool GaussianElimination::solve( const CMatMxN& A, const VecN& b, VecN& x )
{
   if( A.rows() != A.columns() )
      throw std::invalid_argument( "Invalid matrix size" );

   if( A.rows() != b.size() )
      throw std::invalid_argument( "Invalid right-hand side vector size" );

   const size_t n( b.size() );

   // Allocating helper data
   A_ =  A;
   b_ = -b;
   x.resize( n, false );

   size_t pi, pj;
   lastPrecision_ = real(0);

   // Initializing the pivot vector
   DynamicVector<size_t> p( n );
   for( size_t j=0; j<n; ++j ) {
      p[j] = j;
   }

   // Performing the Gaussian elimination
   for( size_t j=0; j<n; ++j )
   {
      size_t max( j );
      real max_val( std::fabs( A_(p[max],j) ) );

      // Partial search for pivot
      for( size_t i=j+1; i<n; ++i ) {
         if( std::fabs( A_(p[i],j) ) > max_val ) {
            max = i;
            max_val = std::fabs( A_(p[max],j) );
         }
      }

      // Swapping rows such the pivot lies on the diagonal
      std::swap( p[max], p[j] );
      pj = p[j];

      if( !isDefault( A_(pj,j) ) )
      {
         // Eliminating the column below the diagonal
         for( size_t i=j+1; i<n; ++i )
         {
            pi = p[i];
            const real f = A_(pi,j) / A_(pj,j);

            reset( A_(pi,j) );

            for( size_t k=j+1; k<n; ++k ) {
               A_(pi,k) -= A_(pj,k) * f;
            }

            b_[pi] -= b_[pj] * f;
         }
      }
      else {
         // Asserting that the column is zero below the diagonal
         for( size_t i=j+1; i<n; ++i ) {
            BLAZE_INTERNAL_ASSERT( isDefault( A_(p[i],j) ), "Fatal error in Gaussian elimination" );
         }
      }
   }

   // Performing the backward substitution
   for( size_t i=n-1; i<n; --i )
   {
      pi = p[i];
      real rhs = b_[pi];

      for( size_t j=i+1; j<n; ++j ) {
         rhs -= x[j] * A_(pi,j);
      }

      if( std::fabs( A_(pi,i) ) > accuracy ) {
         x[i] = rhs / A_(pi,i);
      }
      else {
         // This will introduce errors in the solution
         reset( x[i] );
         lastPrecision_ = max( lastPrecision_, std::fabs( rhs ) );
      }
   }

   BLAZE_LOG_DEBUG_SECTION( log ) {
      if( lastPrecision_ < threshold_ )
         log << "      Solved the linear system using Gaussian elimination.";
      else
         log << BLAZE_YELLOW << "      WARNING: Did not solve the linear system within accuracy. (" << lastPrecision_ << ")" << BLAZE_OLDCOLOR;
   }

   lastIterations_ = 1;

   return lastPrecision_ < threshold_;
}
//*************************************************************************************************

} // namespace blaze
