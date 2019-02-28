//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecReduceExpr.h
//  \brief Header file for the sparse vector reduce expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECREDUCEEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECREDUCEEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/functors/Add.h>
#include <blaze/math/functors/Max.h>
#include <blaze/math/functors/Min.h>
#include <blaze/math/functors/Mult.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Performs a custom reduction operation on the given sparse vector.
// \ingroup sparse_vector
//
// \param sv The given sparse vector for the reduction computation.
// \param op The reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the given sparse vector \a sv by means of the
// given reduction operation \a op:

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization

   const double totalsum1 = reduce( a, blaze::Add() );
   const double totalsum2 = reduce( a, []( double a, double b ){ return a + b; } );
   \endcode

// As demonstrated in the example it is possible to pass any binary callable as custom reduction
// operation. See \ref custom_operations for a detailed overview of the possibilities of custom
// operations.
//
// Please note that the evaluation order of the reduction operation is unspecified. Thus the
// behavior is non-deterministic if \a op is not associative or not commutative. Also, the
// operation is undefined if the given reduction operation modifies the values.
*/
template< typename VT    // Type of the sparse vector
        , bool TF        // Transpose flag
        , typename OP >  // Type of the reduction operation
inline decltype(auto) reduce( const SparseVector<VT,TF>& sv, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using CT = CompositeType_t<VT>;
   using ET = ElementType_t<VT>;

   if( (~sv).size() == 0UL ) return ET{};

   CT tmp( ~sv );

   BLAZE_INTERNAL_ASSERT( tmp.size() == (~sv).size(), "Invalid vector size" );

   const auto end( tmp.end() );
   auto element( tmp.begin() );

   if( element == end ) return ET{};

   ET redux( element->value() );
   ++element;

   for( ; element!=end; ++element ) {
      redux = op( redux, element->value() );
   }

   return redux;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given sparse vector by means of addition.
// \ingroup sparse_vector
//
// \param sv The given sparse vector for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the given sparse vector \a sv by means of
// addition:

   \code
   blaze::CompressedVector<int> a{ 1, 2, 3, 4 };

   const int totalsum = sum( a );  // Results in 10
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) sum( const SparseVector<VT,TF>& sv )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~sv, Add() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduces the given sparse vector by means of multiplication.
// \ingroup sparse_vector
//
// \param sv The given sparse vector for the reduction operation.
// \return The result of the reduction operation.
//
// This function reduces the non-zero elements of the given sparse vector \a sv by means of
// multiplication:

   \code
   blaze::CompressedVector<int> a{ 1, 2, 3, 4 };
   const int totalprod = prod( a );  // Results in 24
   \endcode

// Please note that the evaluation order of the reduction operation is unspecified.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) prod( const SparseVector<VT,TF>& sv )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~sv, Mult() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the smallest element of the sparse vector.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \return The smallest sparse vector element.
//
// This function returns the smallest non-zero element of the given sparse vector. This function
// can only be used for element types that support the smaller-than relationship. In case the
// given vector currently has a size of 0, the returned value is the default value (e.g. 0 in
// case of fundamental data types).
//
// \note In case the sparse vector is not completely filled, the implicit zero elements are NOT
// taken into account. Example: the following compressed vector has only 2 non-zero elements.
// However, the minimum of this vector is 1:

   \code
   blaze::CompressedVector<int> a{ 1, 0, 3, 0 };
   const int totalmin = min( a );  // Results in 1
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) min( const SparseVector<VT,TF>& sv )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~sv, Min() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the largest element of the sparse vector.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \return The largest sparse vector element.
//
// This function returns the largest non-zero element of the given sparse vector. This function
// can only be used for element types that support the smaller-than relationship. In case the
// given vector currently has a size of 0, the returned value is the default value (e.g. 0 in
// case of fundamental data types).
//
// \note In case the sparse vector is not completely filled, the implicit zero elements are NOT
// taken into account. Example: the following compressed vector has only 2 non-zero elements.
// However, the maximum of this vector is -1:

   \code
   blaze::CompressedVector<int> a{ -1, 0, -3, 0 };
   const int totalmin = max( a );  // Results in -1
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) max( const SparseVector<VT,TF>& sv )
{
   BLAZE_FUNCTION_TRACE;

   return reduce( ~sv, Max() );
}
//*************************************************************************************************

} // namespace blaze

#endif
