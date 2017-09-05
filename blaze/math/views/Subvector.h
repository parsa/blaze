//=================================================================================================
/*!
//  \file blaze/math/views/Subvector.h
//  \brief Header file for the implementation of the Subvector view
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_SUBVECTOR_H_
#define _BLAZE_MATH_VIEWS_SUBVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/CrossExpr.h>
#include <blaze/math/expressions/VecEvalExpr.h>
#include <blaze/math/expressions/VecMapExpr.h>
#include <blaze/math/expressions/VecScalarDivExpr.h>
#include <blaze/math/expressions/VecScalarMultExpr.h>
#include <blaze/math/expressions/VecSerialExpr.h>
#include <blaze/math/expressions/Vector.h>
#include <blaze/math/expressions/VecTransExpr.h>
#include <blaze/math/expressions/VecVecAddExpr.h>
#include <blaze/math/expressions/VecVecDivExpr.h>
#include <blaze/math/expressions/VecVecMapExpr.h>
#include <blaze/math/expressions/VecVecMultExpr.h>
#include <blaze/math/expressions/VecVecSubExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/views/subvector/BaseTemplate.h>
#include <blaze/math/views/subvector/Dense.h>
#include <blaze/math/views/subvector/Sparse.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup subvector
//
// \param vector The vector containing the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given vector.
// The following example demonstrates the creation of a dense and sparse subvector:

   \code
   blaze::DynamicVector<double,blaze::columnVector> d;
   blaze::CompressedVector<int,blaze::rowVector> s;
   // ... Resizing and initialization

   // Creating a dense subvector of size 8, starting from index 4
   auto dsv = subvector<4UL,8UL>( d );

   // Creating a sparse subvector of size 7, starting from index 5
   auto ssv = subvector<5UL,7UL>( s );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// Please note that this function creates an unaligned dense or sparse subvector. For instance,
// the creation of the dense subvector is equivalent to the following function call:

   \code
   auto dsv = subvector<unaligned,4UL,8UL>( d );
   \endcode

// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions. However, especially in case of dense subvectors this may
// result in considerable performance improvements. In order to create an aligned subvector the
// following function call has to be used:

   \code
   auto dsv = subvector<aligned,4UL,8UL>( d );
   \endcode

// Note however that in this case the given \a index and \a size are subject to additional checks
// to guarantee proper alignment.
*/
template< size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline decltype(auto) subvector( Vector<VT,TF>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned,I,N>( ~vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given constant vector.
// \ingroup subvector
//
// \param vector The constant vector containing the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given constant
// vector. The following example demonstrates the creation of a dense and sparse subvector:

   \code
   const blaze::DynamicVector<double,blaze::columnVector> d( ... );
   const blaze::CompressedVector<int,blaze::rowVector> s( ... );
   // ... Resizing and initialization

   // Creating a dense subvector of size 8, starting from index 4
   auto dsv = subvector<4UL,8UL>( d );

   // Creating a sparse subvector of size 7, starting from index 5
   auto ssv = subvector<5UL,7UL>( s );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// Please note that this function creates an unaligned dense or sparse subvector. For instance,
// the creation of the dense subvector is equivalent to the following function call:

   \code
   auto dsv = subvector<unaligned,4UL,8UL>( d );
   \endcode

// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions. However, especially in case of dense subvectors this may
// result in considerable performance improvements. In order to create an aligned subvector the
// following function call has to be used:

   \code
   auto dsv = subvector<aligned,4UL,8UL>( d );
   \endcode

// Note however that in this case the given \a index and \a size are subject to additional checks
// to guarantee proper alignment.
*/
template< size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline decltype(auto) subvector( const Vector<VT,TF>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned,I,N>( ~vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given temporary vector.
// \ingroup subvector
//
// \param vector The temporary vector containing the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given
// temporary vector. In case the subvector is not properly specified (i.e. if the specified
// first index is greater than the total size of the given vector or the subvector is specified
// beyond the size of the vector) a \a std::invalid_argument exception is thrown.
*/
template< size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline decltype(auto) subvector( Vector<VT,TF>&& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned,I,N>( ~vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup subvector
//
// \param vector The vector containing the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given dense or sparse vector, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned subvector:

   \code
   blaze::DynamicVector<double,blaze::columnVector> d;
   blaze::CompressedVector<int,blaze::rowVector> s;
   // ... Resizing and initialization

   // Creating an aligned dense subvector of size 8 starting from index 4
   auto dsv = subvector<aligned,4UL,8UL>( d );

   // Creating an unaligned subvector of size 7 starting from index 3
   auto ssv = subvector<unaligned,3UL,7UL>( s );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions and the given \a index is subject to additional checks to
// guarantee proper alignment. However, especially in case of dense subvectors this may result
// in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of the subvector must be aligned. The following source code gives some examples
// for a double precision dynamic vector, assuming that AVX is available, which packs 4 \c double
// values into a SIMD vector:

   \code
   blaze::DynamicVector<double,blaze::columnVector> d( 17UL );
   // ... Resizing and initialization

   // OK: Starts at the beginning, i.e. the first element is aligned
   auto dsv1 = subvector<aligned,0UL,13UL>( d );

   // OK: Start index is a multiple of 4, i.e. the first element is aligned
   auto dsv2 = subvector<aligned,4UL,7UL>( d );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   auto dsv3 = subvector<aligned,8UL,9UL>( d );

   // Error: Start index is not a multiple of 4, i.e. the first element is not aligned
   auto dsv4 = subvector<aligned,5UL,8UL>( d );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline Subvector<VT,AF,I,N> subvector( Vector<VT,TF>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector<VT,AF,I,N>( ~vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given constant vector.
// \ingroup subvector
//
// \param vector The constant vector containing the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given constant dense or sparse vector, based on the specified alignment flag \a AF. The
// following example demonstrates the creation of both an aligned and unaligned subvector:

   \code
   const blaze::DynamicVector<double,blaze::columnVector> d( ... );
   const blaze::CompressedVector<int,blaze::rowVector> s( ... );

   // Creating an aligned dense subvector of size 8 starting from index 4
   auto dsv = subvector<aligned,4UL,8UL>( d );

   // Creating an unaligned subvector of size 7 starting from index 3
   auto ssv = subvector<unaligned,3UL,7UL>( s );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions and the given \a index is subject to additional checks to
// guarantee proper alignment. However, especially in case of dense subvectors this may result
// in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of the subvector must be aligned. The following source code gives some examples
// for a double precision dynamic vector, assuming that AVX is available, which packs 4 \c double
// values into a SIMD vector:

   \code
   const blaze::DynamicVector<double,blaze::columnVector> d( ... );

   // OK: Starts at the beginning, i.e. the first element is aligned
   auto dsv1 = subvector<aligned,0UL,13UL>( d );

   // OK: Start index is a multiple of 4, i.e. the first element is aligned
   auto dsv2 = subvector<aligned,4UL,7UL>( d );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   auto dsv3 = subvector<aligned,8UL,9UL>( d );

   // Error: Start index is not a multiple of 4, i.e. the first element is not aligned
   auto dsv4 = subvector<aligned,5UL,8UL>( d );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline const Subvector<const VT,AF,I,N>
   subvector( const Vector<VT,TF>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector<const VT,AF,I,N>( ~vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given temporary vector.
// \ingroup subvector
//
// \param vector The temporary vector containing the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given temporary dense or sparse vector, based on the specified alignment flag \a AF. In
// case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of
// the vector) or any alignment restrictions are violated, a \a std::invalid_argument exception
// is thrown.
*/
template< bool AF      // Alignment flag
        , size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline Subvector<VT,AF,I,N> subvector( Vector<VT,TF>&& vector )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector<VT,AF,I,N>( ~vector );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup subvector
//
// \param vector The vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given vector.
// The following example demonstrates the creation of a dense and sparse subvector:

   \code
   blaze::DynamicVector<double,blaze::columnVector> d;
   blaze::CompressedVector<int,blaze::rowVector> s;
   // ... Resizing and initialization

   // Creating a dense subvector of size 8, starting from index 4
   auto dsv = subvector( d, 4UL, 8UL );

   // Creating a sparse subvector of size 7, starting from index 5
   auto ssv = subvector( s, 5UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// Please note that this function creates an unaligned dense or sparse subvector. For instance,
// the creation of the dense subvector is equivalent to the following function call:

   \code
   auto dsv = subvector<unaligned>( d, 4UL, 8UL );
   \endcode

// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions. However, especially in case of dense subvectors this may
// result in considerable performance improvements. In order to create an aligned subvector the
// following function call has to be used:

   \code
   auto dsv = subvector<aligned>( d, 4UL, 8UL );
   \endcode

// Note however that in this case the given \a index and \a size are subject to additional checks
// to guarantee proper alignment.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline decltype(auto) subvector( Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given constant vector.
// \ingroup subvector
//
// \param vector The constant vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given constant
// vector. The following example demonstrates the creation of a dense and sparse subvector:

   \code
   const blaze::DynamicVector<double,blaze::columnVector> d( ... );
   const blaze::CompressedVector<int,blaze::rowVector> s( ... );
   // ... Resizing and initialization

   // Creating a dense subvector of size 8, starting from index 4
   auto dsv = subvector( d, 4UL, 8UL );

   // Creating a sparse subvector of size 7, starting from index 5
   auto ssv = subvector( s, 5UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// Please note that this function creates an unaligned dense or sparse subvector. For instance,
// the creation of the dense subvector is equivalent to the following function call:

   \code
   auto dsv = subvector<unaligned>( d, 4UL, 8UL );
   \endcode

// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions. However, especially in case of dense subvectors this may
// result in considerable performance improvements. In order to create an aligned subvector the
// following function call has to be used:

   \code
   auto dsv = subvector<aligned>( d, 4UL, 8UL );
   \endcode

// Note however that in this case the given \a index and \a size are subject to additional checks
// to guarantee proper alignment.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline decltype(auto) subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given temporary vector.
// \ingroup subvector
//
// \param vector The temporary vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given
// temporary vector. In case the subvector is not properly specified (i.e. if the specified
// first index is greater than the total size of the given vector or the subvector is specified
// beyond the size of the vector) a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline decltype(auto) subvector( Vector<VT,TF>&& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<unaligned>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given vector.
// \ingroup subvector
//
// \param vector The vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given dense or sparse vector, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned subvector:

   \code
   blaze::DynamicVector<double,blaze::columnVector> d;
   blaze::CompressedVector<int,blaze::rowVector> s;
   // ... Resizing and initialization

   // Creating an aligned dense subvector of size 8 starting from index 4
   auto dsv = subvector<aligned>( d, 4UL, 8UL );

   // Creating an unaligned subvector of size 7 starting from index 3
   auto ssv = subvector<unaligned>( s, 3UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions and the given \a index is subject to additional checks to
// guarantee proper alignment. However, especially in case of dense subvectors this may result
// in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of the subvector must be aligned. The following source code gives some examples
// for a double precision dynamic vector, assuming that AVX is available, which packs 4 \c double
// values into a SIMD vector:

   \code
   blaze::DynamicVector<double,blaze::columnVector> d( 17UL );
   // ... Resizing and initialization

   // OK: Starts at the beginning, i.e. the first element is aligned
   auto dsv1 = subvector<aligned>( d, 0UL, 13UL );

   // OK: Start index is a multiple of 4, i.e. the first element is aligned
   auto dsv2 = subvector<aligned>( d, 4UL, 7UL );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   auto dsv3 = subvector<aligned>( d, 8UL, 9UL );

   // Error: Start index is not a multiple of 4, i.e. the first element is not aligned
   auto dsv4 = subvector<aligned>( d, 5UL, 8UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline Subvector<VT,AF> subvector( Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector<VT,AF>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given constant vector.
// \ingroup subvector
//
// \param vector The constant vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given constant dense or sparse vector, based on the specified alignment flag \a AF. The
// following example demonstrates the creation of both an aligned and unaligned subvector:

   \code
   const blaze::DynamicVector<double,blaze::columnVector> d( ... );
   const blaze::CompressedVector<int,blaze::rowVector> s( ... );

   // Creating an aligned dense subvector of size 8 starting from index 4
   auto dsv = subvector<aligned>( d, 4UL, 8UL );

   // Creating an unaligned subvector of size 7 starting from index 3
   auto ssv = subvector<unaligned>( s, 3UL, 7UL );
   \endcode

// In case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of the
// vector) a \a std::invalid_argument exception is thrown.
//
// In contrast to unaligned subvectors, which provide full flexibility, aligned subvectors pose
// additional alignment restrictions and the given \a index is subject to additional checks to
// guarantee proper alignment. However, especially in case of dense subvectors this may result
// in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). In order to be properly aligned the
// first element of the subvector must be aligned. The following source code gives some examples
// for a double precision dynamic vector, assuming that AVX is available, which packs 4 \c double
// values into a SIMD vector:

   \code
   const blaze::DynamicVector<double,blaze::columnVector> d( ... );

   // OK: Starts at the beginning, i.e. the first element is aligned
   auto dsv1 = subvector<aligned>( d, 0UL, 13UL );

   // OK: Start index is a multiple of 4, i.e. the first element is aligned
   auto dsv2 = subvector<aligned>( d, 4UL, 7UL );

   // OK: The start index is a multiple of 4 and the subvector includes the last element
   auto dsv3 = subvector<aligned>( d, 8UL, 9UL );

   // Error: Start index is not a multiple of 4, i.e. the first element is not aligned
   auto dsv4 = subvector<aligned>( d, 5UL, 8UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline const Subvector<const VT,AF>
   subvector( const Vector<VT,TF>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector<const VT,AF>( ~vector, index, size );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific subvector of the given temporary vector.
// \ingroup subvector
//
// \param vector The temporary vector containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specific subvector of the vector.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing an aligned or unaligned subvector of the
// given temporary dense or sparse vector, based on the specified alignment flag \a AF. In
// case the subvector is not properly specified (i.e. if the specified first index is greater
// than the total size of the given vector or the subvector is specified beyond the size of
// the vector) or any alignment restrictions are violated, a \a std::invalid_argument exception
// is thrown.
*/
template< bool AF      // Alignment flag
        , typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline Subvector<VT,AF> subvector( Vector<VT,TF>&& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector<VT,AF>( ~vector, index, size );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector addition.
// \ingroup subvector
//
// \param vector The constant vector/vector addition.
// \return View on the specified subvector of the addition.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector addition.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecAddExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF,I,N>( (~vector).leftOperand() ) +
          subvector<AF,I,N>( (~vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector addition.
// \ingroup subvector
//
// \param vector The constant vector/vector addition.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the addition.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector addition.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecAddExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) +
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector subtraction.
// \ingroup subvector
//
// \param vector The constant vector/vector subtraction.
// \return View on the specified subvector of the subtraction.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector subtraction.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecSubExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF,I,N>( (~vector).leftOperand() ) -
          subvector<AF,I,N>( (~vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector subtraction.
// \ingroup subvector
//
// \param vector The constant vector/vector subtraction.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the subtraction.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector subtraction.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecSubExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) -
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector multiplication.
// \ingroup subvector
//
// \param vector The constant vector/vector multiplication.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector multiplication.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF,I,N>( (~vector).leftOperand() ) *
          subvector<AF,I,N>( (~vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector multiplication.
// \ingroup subvector
//
// \param vector The constant vector/vector multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector multiplication.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecMultExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) *
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector division.
// \ingroup subvector
//
// \param vector The constant vector/vector division.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector division.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecDivExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF,I,N>( (~vector).leftOperand() ) /
          subvector<AF,I,N>( (~vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector division.
// \ingroup subvector
//
// \param vector The constant vector/vector division.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector division.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecDivExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand() , index, size ) /
          subvector<AF>( (~vector).rightOperand(), index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector cross product.
// \ingroup subvector
//
// \param vector The constant vector/vector cross product.
// \return View on the specified subvector of the cross product.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector cross product.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const CrossExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector< VectorType_<VT>, AF, I, N >( ~vector );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/vector cross product.
// \ingroup subvector
//
// \param vector The constant vector/vector cross product.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the cross product.
//
// This function returns an expression representing the specified subvector of the given
// vector/vector cross product.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const CrossExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return Subvector< VectorType_<VT>, AF >( ~vector, index, size );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar multiplication.
// \ingroup subvector
//
// \param vector The constant vector/scalar multiplication.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar multiplication.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecScalarMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF,I,N>( (~vector).leftOperand() ) * (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar multiplication.
// \ingroup subvector
//
// \param vector The constant vector/scalar multiplication.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the multiplication.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar multiplication.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecScalarMultExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand(), index, size ) * (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar division.
// \ingroup subvector
//
// \param vector The constant vector/scalar division.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar division.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecScalarDivExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF,I,N>( (~vector).leftOperand() ) / (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector/scalar division.
// \ingroup subvector
//
// \param vector The constant vector/scalar division.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the division.
//
// This function returns an expression representing the specified subvector of the given
// vector/scalar division.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecScalarDivExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~vector).leftOperand(), index, size ) / (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given unary vector map operation.
// \ingroup subvector
//
// \param vector The constant unary vector map operation.
// \return View on the specified subvector of the unary map operation.
//
// This function returns an expression representing the specified subvector of the given unary
// vector map operation.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecMapExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return map( subvector<AF,I,N>( (~vector).operand() ), (~vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given unary vector map operation.
// \ingroup subvector
//
// \param vector The constant unary vector map operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the unary map operation.
//
// This function returns an expression representing the specified subvector of the given unary
// vector map operation.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecMapExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return map( subvector<AF>( (~vector).operand(), index, size ), (~vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given binary vector map operation.
// \ingroup subvector
//
// \param vector The constant binary vector map operation.
// \return View on the specified subvector of the binary map operation.
//
// This function returns an expression representing the specified subvector of the given binary
// vector map operation.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecMapExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return map( subvector<AF,I,N>( (~vector).leftOperand() ),
               subvector<AF,I,N>( (~vector).rightOperand() ),
               (~vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given binary vector map operation.
// \ingroup subvector
//
// \param vector The constant binary vector map operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the binary map operation.
//
// This function returns an expression representing the specified subvector of the given binary
// vector map operation.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecVecMapExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return map( subvector<AF>( (~vector).leftOperand() , index, size ),
               subvector<AF>( (~vector).rightOperand(), index, size ),
               (~vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector evaluation operation.
// \ingroup subvector
//
// \param vector The constant vector evaluation operation.
// \return View on the specified subvector of the evaluation operation.
//
// This function returns an expression representing the specified subvector of the given vector
// evaluation operation.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecEvalExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return eval( subvector<AF,I,N>( (~vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector evaluation operation.
// \ingroup subvector
//
// \param vector The constant vector evaluation operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the evaluation operation.
//
// This function returns an expression representing the specified subvector of the given vector
// evaluation operation.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecEvalExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return eval( subvector<AF>( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector serialization operation.
// \ingroup subvector
//
// \param vector The constant vector serialization operation.
// \return View on the specified subvector of the serialization operation.
//
// This function returns an expression representing the specified subvector of the given vector
// serialization operation.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecSerialExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return serial( subvector<AF,I,N>( (~vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector serialization operation.
// \ingroup subvector
//
// \param vector The constant vector serialization operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the serialization operation.
//
// This function returns an expression representing the specified subvector of the given vector
// serialization operation.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecSerialExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return serial( subvector<AF>( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector transpose operation.
// \ingroup subvector
//
// \param vector The constant vector transpose operation.
// \return View on the specified subvector of the transpose operation.
//
// This function returns an expression representing the specified subvector of the given vector
// transpose operation.
*/
template< bool AF        // Alignment flag
        , size_t I       // Index of the first subvector element
        , size_t N       // Size of the subvector
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecTransExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return trans( subvector<AF,I,N>( (~vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given vector transpose operation.
// \ingroup subvector
//
// \param vector The constant vector transpose operation.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the transpose operation.
//
// This function returns an expression representing the specified subvector of the given vector
// transpose operation.
*/
template< bool AF        // Alignment flag
        , typename VT >  // Vector base type of the expression
inline decltype(auto) subvector( const VecTransExpr<VT>& vector, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   return trans( subvector<AF>( (~vector).operand(), index, size ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of another subvector.
// \ingroup subvector
//
// \param sv The constant subvector.
// \return View on the specified subvector of the other subvector.
//
// This function returns an expression representing the specified subvector of the given subvector.
*/
template< bool AF1     // Required alignment flag
        , size_t I     // Index of the first subvector element
        , size_t N     // Size of the subvector
        , typename VT  // Type of the dense vector
        , bool AF2 >   // Present alignment flag
inline const Subvector<VT,AF1>
   subvector( const Subvector<VT,AF2>& sv )
{
   BLAZE_FUNCTION_TRACE;

   if( I + N > sv.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }

   return Subvector<VT,AF1>( sv.operand(), sv.offset() + I, N );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of another subvector.
// \ingroup subvector
//
// \param sv The constant subvector.
// \return View on the specified subvector of the other subvector.
//
// This function returns an expression representing the specified subvector of the given subvector.
*/
template< bool AF1     // Required alignment flag
        , size_t I     // Required subvector offset
        , size_t N     // Required size of the subvector
        , typename VT  // Type of the dense vector
        , bool AF2     // Present alignment flag
        , size_t I2    // Present subvector offset
        , size_t N2 >  // Present size of the subvector
inline const Subvector<VT,AF1,I+I2,N>
   subvector( const Subvector<VT,AF2,I2,N2>& sv )
{
   BLAZE_FUNCTION_TRACE;

   if( I + N > sv.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }

   return Subvector<VT,AF1,I+I2,N>( sv.operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of another subvector.
// \ingroup subvector
//
// \param sv The constant subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \return View on the specified subvector of the other subvector.
//
// This function returns an expression representing the specified subvector of the given subvector.
*/
template< bool AF1         // Required alignment flag
        , typename VT      // Type of the dense vector
        , bool AF2         // Present alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline const Subvector<VT,AF1>
   subvector( const Subvector<VT,AF2,SAs...>& sv, size_t index, size_t size )
{
   BLAZE_FUNCTION_TRACE;

   if( index + size > sv.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }

   return Subvector<VT,AF1>( sv.operand(), sv.offset() + index, size );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Subvector operators */
//@{
template< typename VT, bool AF, size_t... SAs >
inline void reset( Subvector<VT,AF,SAs...>& sv );

template< typename VT, bool AF, size_t... SAs >
inline void reset( Subvector<VT,AF,SAs...>&& sv );

template< typename VT, bool AF, size_t... SAs >
inline void clear( Subvector<VT,AF,SAs...>& sv );

template< typename VT, bool AF, size_t... SAs >
inline void clear( Subvector<VT,AF,SAs...>&& sv );

template< bool RF, typename VT, bool AF, size_t... SAs >
inline bool isDefault( const Subvector<VT,AF,SAs...>& sv );

template< typename VT, bool AF, size_t... SAs >
inline bool isIntact( const Subvector<VT,AF,SAs...>& sv ) noexcept;

template< typename VT, bool AF, size_t... SAs, bool TF >
inline bool isSame( const Subvector<VT,AF,SAs...>& a, const Vector<VT,TF>& b ) noexcept;

template< typename VT, bool TF, bool AF, size_t... SAs >
inline bool isSame( const Vector<VT,TF>& a, const Subvector<VT,AF,SAs...>& b ) noexcept;

template< typename VT1, bool AF1, size_t... SAs1, typename VT2, bool AF2, size_t... SAs2 >
inline bool isSame( const Subvector<VT1,AF1,SAs1...>& a, const Subvector<VT2,AF2,SAs2...>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given subvector.
// \ingroup subvector
//
// \param sv The subvector to be resetted.
// \return void
*/
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline void reset( Subvector<VT,AF,SAs...>& sv )
{
   sv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given temporary subvector.
// \ingroup subvector
//
// \param sv The temporary subvector to be resetted.
// \return void
*/
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline void reset( Subvector<VT,AF,SAs...>&& sv )
{
   sv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given subvector.
// \ingroup subvector
//
// \param sv The subvector to be cleared.
// \return void
*/
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline void clear( Subvector<VT,AF,SAs...>& sv )
{
   sv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given temporary subvector.
// \ingroup subvector
//
// \param sv The temporary subvector to be cleared.
// \return void
*/
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline void clear( Subvector<VT,AF,SAs...>&& sv )
{
   sv.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense subvector is in default state.
// \ingroup subvector
//
// \param sv The dense subvector to be tested for its default state.
// \return \a true in case the given dense subvector is component-wise zero, \a false otherwise.
//
// This function checks whether the dense row is in default state. For instance, in case the
// subvector is instantiated for a vector of built-in integral or floating point data type,
// the function returns \a true in case all subvector elements are 0 and \a false in case any
// subvector element is not 0.
*/
template< bool RF          // Relaxation flag
        , typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline bool isDefault_backend( const DenseSubvector<VT,AF,SAs...>& sv )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<sv.size(); ++i )
      if( !isDefault<RF>( sv[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse subvector is in default state.
// \ingroup subvector
//
// \param row The sparse subvector to be tested for its default state.
// \return \a true in case the given subvector is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse subvector is in default state. For instance, in case
// the subvector is instantiated for a vector of built-in integral or floating point data type,
// the function returns \a true in case all subvector elements are 0 and \a false in case any
// element is not 0.
*/
template< bool RF          // Relaxation flag
        , typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline bool isDefault_backend( const SparseSubvector<VT,AF,SAs...>& sv )
{
   using blaze::isDefault;

   for( const auto& element : ~sv )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given subvector is in default state.
// \ingroup subvector
//
// \param sv The subvector to be tested for its default state.
// \return \a true in case the given subvector is component-wise zero, \a false otherwise.
//
// This function checks whether the subvector is in default state. For instance, in case the
// subvector is instantiated for a vector of built-in integral or floating point data type,
// the function returns \a true in case all subvector elements are 0 and \a false in case any
// subvector element is not 0. The following example demonstrates the use of the \a isDefault
// function:

   \code
   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isDefault( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode
*/
template< bool RF          // Relaxation flag
        , typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline bool isDefault( const Subvector<VT,AF,SAs...>& sv )
{
   return isDefault_backend<RF>( ~sv );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given subvector vector are intact.
// \ingroup subvector
//
// \param sv The subvector to be tested.
// \return \a true in case the given subvector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the subvector are intact, i.e. if its state
// is valid. In case the invariants are intact, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isIntact( subvector( v, 10UL, 20UL ) ) ) { ... }
   \endcode
*/
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline bool isIntact( const Subvector<VT,AF,SAs...>& sv ) noexcept
{
   return ( sv.offset() + sv.size() <= sv.operand().size() &&
            isIntact( sv.operand() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given vector and subvector represent the same observable state.
// \ingroup subvector
//
// \param a The subvector to be tested for its state.
// \param b The vector to be tested for its state.
// \return \a true in case the subvector and vector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given subvector refers to the entire
// range of the given vector and by that represents the same observable state. In this case,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename VT    // Type of the vector
        , bool AF        // Alignment flag
        , size_t... SAs  // Compile time subvector arguments
        , bool TF >      // Transpose flag
inline bool isSame( const Subvector<VT,AF,SAs...>& a, const Vector<VT,TF>& b ) noexcept
{
   return ( isSame( a.operand(), ~b ) && ( a.size() == (~b).size() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given vector and subvector represent the same observable state.
// \ingroup subvector
//
// \param a The vector to be tested for its state.
// \param b The subvector to be tested for its state.
// \return \a true in case the vector and subvector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given subvector refers to the entire
// range of the given vector and by that represents the same observable state. In this case,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename VT      // Type of the vector
        , bool TF          // Transpose flag
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline bool isSame( const Vector<VT,TF>& a, const Subvector<VT,AF,SAs...>& b ) noexcept
{
   return ( isSame( ~a, b.operand() ) && ( (~a).size() == b.size() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given subvectors represent the same observable state.
// \ingroup subvector
//
// \param a The first subvector to be tested for its state.
// \param b The second subvector to be tested for its state.
// \return \a true in case the two subvectors share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given subvectors refer to exactly the
// same range of the same vector. In case both subvectors represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename VT1      // Type of the vector of the left-hand side subvector
        , bool AF1          // Alignment flag of the left-hand side subvector
        , size_t... SAs1    // Compile time subvector arguments of the left-hand side subvector
        , typename VT2      // Type of the vector of the right-hand side subvector
        , bool AF2          // Alignment flag of the right-hand side subvector
        , size_t... SAs2 >  // Compile time subvector arguments of the right-hand side subvector
inline bool isSame( const Subvector<VT1,AF1,SAs1...>& a, const Subvector<VT2,AF2,SAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand() ) &&
            ( a.offset() == b.offset() ) &&
            ( a.size() == b.size() ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a subvector.
// \ingroup subvector
//
// \param lhs The target left-hand side subvector.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the vector
        , bool AF         // Alignment flag
        , size_t... SAs   // Compile time subvector arguments
        , typename VT2    // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryAssign( const Subvector<VT1,AF,SAs...>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.operand(), ~rhs, lhs.offset() + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a subvector.
// \ingroup subvector
//
// \param lhs The target left-hand side subvector.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the vector
        , bool AF         // Alignment flag
        , size_t... SAs   // Compile time subvector arguments
        , typename VT2    // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const Subvector<VT1,AF,SAs...>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.operand(), ~rhs, lhs.offset() + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a subvector.
// \ingroup subvector
//
// \param lhs The target left-hand side subvector.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the vector
        , bool AF         // Alignment flag
        , size_t... SAs   // Compile time subvector arguments
        , typename VT2    // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool trySubAssign( const Subvector<VT1,AF,SAs...>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.operand(), ~rhs, lhs.offset() + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a subvector.
// \ingroup subvector
//
// \param lhs The target left-hand side subvector.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the vector
        , bool AF         // Alignment flag
        , size_t... SAs   // Compile time subvector arguments
        , typename VT2    // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryMultAssign( const Subvector<VT1,AF,SAs...>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.operand(), ~rhs, lhs.offset() + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a subvector.
// \ingroup subvector
//
// \param lhs The target left-hand side subvector.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1    // Type of the vector
        , bool AF         // Alignment flag
        , size_t... SAs   // Compile time subvector arguments
        , typename VT2    // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryDivAssign( const Subvector<VT1,AF,SAs...>& lhs, const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryDivAssign( lhs.operand(), ~rhs, lhs.offset() + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given subvector.
// \ingroup subvector
//
// \param sv The subvector to be derestricted.
// \return Subvector without access restrictions.
//
// This function removes all restrictions on the data access to the given subvector. It returns a
// subvector that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline decltype(auto) derestrict( Subvector<VT,AF,SAs...>& sv )
{
   return subvector( derestrict( sv.operand() ), sv.offset(), sv.size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary subvector.
// \ingroup subvector
//
// \param sv The temporary subvector to be derestricted.
// \return Subvector without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary subvector. It
// returns a subvector that does provide the same interface but does not have any restrictions on
// the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT      // Type of the vector
        , bool AF          // Alignment flag
        , size_t... SAs >  // Compile time subvector arguments
inline decltype(auto) derestrict( Subvector<VT,AF,SAs...>&& sv )
{
   return subvector( derestrict( sv.operand() ), sv.offset(), sv.size() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESTRICTED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF, size_t... SAs >
struct IsRestricted< Subvector<VT,AF,SAs...> >
   : public BoolConstant< IsRestricted<VT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF, size_t... SAs >
struct HasConstDataAccess< DenseSubvector<VT,AF,SAs...> >
   : public BoolConstant< HasConstDataAccess<VT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF, size_t... SAs >
struct HasMutableDataAccess< DenseSubvector<VT,AF,SAs...> >
   : public BoolConstant< HasMutableDataAccess<VT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool AF, size_t... SAs >
struct IsAligned< AlignedSubvector<VT,AF,SAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
