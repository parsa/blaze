//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatVarExpr.h
//  \brief Header file for the sparse matrix variance expression
//
//  Copyright (C) 2012-2019 Klaus Iglberger - All Rights Reserved
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
//=================================================================================================

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATVAREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATVAREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Exception.h>
#include <blaze/math/expressions/ScalarExpandExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Computes the variance for the given sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the variance computation.
// \return The variance of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the <a href="https://en.wikipedia.org/wiki/Variance">variance</a> for
// the given sparse matrix \a sm. Both the non-zero and zero elements of the sparse matrix are
// taken into account. Example:

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int> A{ { 1, 3, 2 }
                          , { 2, 6, 4 }
                          , { 9, 6, 3 } };

   const double v = var( A );  // Results in 6.5
   \endcode

// In case the size of the given matrix is smaller than 2, a \a std::invalid_argument is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) var( const SparseMatrix<MT,SO>& sm )
{
   using BT = UnderlyingBuiltin_t<MT>;

   const size_t n( size( ~sm ) );

   if( n < 2UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix" );
   }

   const auto m( expand( mean( ~sm ), rows( ~sm ), columns( ~sm ) ) );

   return sum( map( (~sm) - m, Pow2() ) ) * inv( BT( n-1UL ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the row-/columnwise variance function for the given sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the variance computation.
// \return The row-/columnwise variance of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the row-/columnwise
// <a href="https://en.wikipedia.org/wiki/Variance">variance</a> for the given sparse matrix
// \a sm. In case \a RF is set to \a rowwise, the function returns a column vector containing
// the variance of each row of \a sm. In case \a RF is set to \a columnwise, the function
// returns a row vector containing the variance of each column of \a dm. Both the non-zero
// and zero elements of the sparse matrix are taken into account. Example:

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int> A{ { 1, 3, 2 }
                          , { 2, 6, 4 }
                          , { 9, 6, 3 } };

   const auto v1 = var<rowwise>( A );     // Results in ( 1  4  9 )
   const auto v2 = var<columnwise>( A );  // Results in ( 19  3  1 )
   \endcode

// In case \a RF is set to \a rowwise and the number of columns of the given matrix is smaller
// than 2 or in case \a RF is set to \a columnwise and the number of rows of the given matrix is
// smaller than 2, a \a std::invalid_argument is thrown.
*/
template< bool RF      // Reduction flag
        , typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) var( const SparseMatrix<MT,SO>& sm )
{
   using BT = UnderlyingBuiltin_t<MT>;

   const size_t n( RF == rowwise ? columns( ~sm ) : rows( ~sm ) );

   if( n < 2UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix" );
   }

   const auto m( expand( mean<RF>( ~sm ), n ) );

   return sum<RF>( map( (~sm) - m, Pow2() ) ) * inv( BT( n-1UL ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
