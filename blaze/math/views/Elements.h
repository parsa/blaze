//=================================================================================================
/*!
//  \file blaze/math/views/Elements.h
//  \brief Header file for the implementation of the Elements view
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_VIEWS_ELEMENTS_H_
#define _BLAZE_MATH_VIEWS_ELEMENTS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
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
#include <blaze/math/InitializerList.h>
#include <blaze/math/IntegerSequence.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/elements/BaseTemplate.h>
#include <blaze/math/views/elements/Dense.h>
#include <blaze/math/views/elements/Sparse.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/PtrdiffT.h>
#include <blaze/util/SmallVector.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/TypeList.h>
#include <blaze/util/Types.h>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a selection of elements of the given vector.
// \ingroup elements
//
// \param vector The vector containing the elements.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given vector.

   \code
   using blaze::rowVector;

   blaze::DynamicVector<double,rowVector> d;
   blaze::CompressedVector<double,rowVector> s;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd element of the dense vector d
   auto elements1 = elements<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd element of the sparse vector s
   auto elements2 = elements<4UL,2UL>( S );
   \endcode

// By default, the provided element indices are checked at runtime. In case any element is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of elements in the given vector) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto elements1 = elements<1UL,3UL>( d, unchecked );
   auto elements2 = elements<4UL,2UL>( s, unchecked );
   \endcode
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( Vector<VT,TF>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Elements_<VT,I,Is...>;
   return ReturnType( ~vector, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of elements of the given constant vector.
// \ingroup elements
//
// \param vector The constant vector containing the elements.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given constant
// vector.

   \code
   using blaze::rowVector;

   const blaze::DynamicVector<double,rowVector> d( ... );
   const blaze::CompressedVector<double,rowVector> s( ... );

   // Creating a view on the 1st and 3rd element of the dense vector d
   auto elements1 = elements<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd element of the sparse vector s
   auto elements2 = elements<4UL,2UL>( S );
   \endcode

// By default, the provided element indices are checked at runtime. In case any element is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of elements in the given vector) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto elements1 = elements<1UL,3UL>( d, unchecked );
   auto elements2 = elements<4UL,2UL>( s, unchecked );
   \endcode
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( const Vector<VT,TF>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Elements_<const VT,I,Is...>;
   return ReturnType( ~vector, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of elements of the given temporary vector.
// \ingroup elements
//
// \param vector The temporary vector containing the elements.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing an selection of elements of the given temporary
// vector. In case any element is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of the elements in the given vector) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( Vector<VT,TF>&& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Elements_<VT,I,Is...>;
   return ReturnType( ~vector, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of elements of the given vector.
// \ingroup elements
//
// \param vector The vector containing the elements.
// \param indices Pointer to the first index of the selected elements.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given vector.

   \code
   using blaze::rowVector;

   blaze::DynamicVector<double,rowVector> d;
   blaze::CompressedVector<double,rowVector> s;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd element of the dense vector d
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto elements1 = elements( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd element of the sparse vector s
   const std::array<size_t,2UL> indices2{ 4UL, 2UL };
   auto elements2 = elements( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided element indices are checked at runtime. In case any element is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of elements in the given vector) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto elements1 = elements( d, indices1.data(), indices1.size(), unchecked );
   auto elements2 = elements( s, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( Vector<VT,TF>& vector, const T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Elements_<VT>;
   return ReturnType( ~vector, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of elements of the given constant vector.
// \ingroup elements
//
// \param vector The constant vector containing the elements.
// \param indices Pointer to the first index of the selected elements.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given constant
// vector.

   \code
   using blaze::rowVector;

   const blaze::DynamicVector<double,rowVector> d( ... );
   const blaze::CompressedVector<double,rowVector> s( ... );

   // Creating a view on the 1st and 3rd element of the dense vector d
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto elements1 = elements( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd element of the sparse vector s
   const std::array<size_t,2UL> indices2{ 4UL, 2UL };
   auto elements2 = elements( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided element indices are checked at runtime. In case any element is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of elements in the given vector) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto elements1 = elements( d, indices1.data(), indices1.size(), unchecked );
   auto elements2 = elements( s, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( const Vector<VT,TF>& vector, const T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Elements_<const VT>;
   return ReturnType( ~vector, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of elements of the given temporary vector.
// \ingroup elements
//
// \param vector The temporary vector containing the elements.
// \param indices Pointer to the first index of the selected elements.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given temporary
// vector. In case any element is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of elements in the given vector) a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( Vector<VT,TF>&& vector, const T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Elements_<VT>;
   return ReturnType( ~vector, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given vector.
// \ingroup elements
//
// \param vector The vector containing the elements.
// \param indices The sequence of element indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given vector.
// In case any element is not properly specified (i.e. if any specified index is greater than
// or equal to the total number of elements in the given vector) a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the vector
        , size_t... Is        // Element indices
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( VT&& vector, index_sequence<Is...> indices, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   UNUSED_PARAMETER( indices );

   return elements<Is...>( std::forward<VT>( vector ), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given vector.
// \ingroup elements
//
// \param vector The vector containing the elements.
// \param indices The list of element indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given vector.
// In case any element is not properly specified (i.e. if any specified index is greater than
// or equal to the total number of elements in the given vector) a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the vector
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( VT&& vector, initializer_list<T> indices, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( std::forward<VT>( vector ), indices.begin(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given vector.
// \ingroup elements
//
// \param vector The vector containing the elements.
// \param indices The array of element indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given vector.
// In case any element is not properly specified (i.e. if any specified index is greater than
// or equal to the total number of elements in the given vector) a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the vector
        , typename T          // Type of the element indices
        , size_t N            // Number of indices
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( VT&& vector, const std::array<T,N>& indices, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( std::forward<VT>( vector ), indices.data(), N, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given vector.
// \ingroup elements
//
// \param vector The vector containing the elements.
// \param indices The vector of element indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given vector.
// In case any element is not properly specified (i.e. if any specified index is greater than
// or equal to the total number of elements in the given vector) a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the vector
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( VT&& vector, const std::vector<T>& indices, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( std::forward<VT>( vector ), indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given vector.
// \ingroup elements
//
// \param vector The vector containing the elements.
// \param indices The vector of element indices.
// \param args Optional arguments.
// \return View on the specified elements of the vector.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing a selection of elements of the given vector.
// In case any element is not properly specified (i.e. if any specified index is greater than
// or equal to the total number of elements in the given vector) a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the vector
        , typename T          // Type of the element indices
        , size_t N            // Number of preallocated elements
        , typename... REAs >  // Optional arguments
inline decltype(auto) elements( VT&& vector, const SmallVector<T,N>& indices, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( std::forward<VT>( vector ), indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector/vector addition.
// \ingroup elements
//
// \param vector The constant vector/vector addition.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the addition.
//
// This function returns an expression representing the specified selection of elements on the
// given vector/vector addition.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecVecAddExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<CEAs...>( (~vector).leftOperand(), args... ) +
          elements<CEAs...>( (~vector).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector/vector subtraction.
// \ingroup elements
//
// \param vector The constant vector/vector subtraction.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the subtraction.
//
// This function returns an expression representing the specified selection of elements on the
// given vector/vector subtraction.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecVecSubExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<CEAs...>( (~vector).leftOperand(), args... ) -
          elements<CEAs...>( (~vector).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector/vector multiplication.
// \ingroup elements
//
// \param vector The constant vector/vector multiplication.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the multiplication.
//
// This function returns an expression representing the specified selection of elements on the
// given vector/vector multiplication.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecVecMultExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<CEAs...>( (~vector).leftOperand(), args... ) *
          elements<CEAs...>( (~vector).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector/vector division.
// \ingroup elements
//
// \param vector The constant vector/vector division.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the division.
//
// This function returns an expression representing the specified selection of elements on the
// given vector/vector division.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecVecDivExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<CEAs...>( (~vector).leftOperand(), args... ) /
          elements<CEAs...>( (~vector).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector/vector cross product.
// \ingroup elements
//
// \param vector The constant vector/vector cross product.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the cross product.
//
// This function returns an expression representing the specified selection of elements on the
// given vector/vector cross product.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const CrossExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Elements_< VectorType_t<VT>, CEAs... >;
   return ReturnType( ~vector, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector/scalar multiplication.
// \ingroup elements
//
// \param vector The constant vector/scalar multiplication.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the multiplication.
//
// This function returns an expression representing the specified selection of elements on the
// given vector/scalar multiplication.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecScalarMultExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<CEAs...>( (~vector).leftOperand(), args... ) * (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector/scalar division.
// \ingroup elements
//
// \param vector The constant vector/scalar division.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the division.
//
// This function returns an expression representing the specified selection of elements on the
// given vector/scalar division.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecScalarDivExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<CEAs...>( (~vector).leftOperand(), args... ) / (~vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given unary vector map operation.
// \ingroup elements
//
// \param vector The constant unary vector map operation.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the unary map operation.
//
// This function returns an expression representing the specified selection of elements on the
// given unary vector map operation.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecMapExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( elements<CEAs...>( (~vector).operand(), args... ), (~vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given binary vector map operation.
// \ingroup elements
//
// \param vector The constant binary vector map operation.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the binary map operation.
//
// This function returns an expression representing the specified selection of elements on the
// given binary vector map operation.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecVecMapExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( elements<CEAs...>( (~vector).leftOperand(), args... ),
               elements<CEAs...>( (~vector).rightOperand(), args... ),
               (~vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector evaluation operation.
// \ingroup elements
//
// \param vector The constant vector evaluation operation.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the evaluation operation.
//
// This function returns an expression representing the specified selection of elements on the
// given vector evaluation operation.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecEvalExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( elements<CEAs...>( (~vector).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector serialization operation.
// \ingroup elements
//
// \param vector The constant vector serialization operation.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the serialization operation.
//
// This function returns an expression representing the specified selection of elements on the
// given vector serialization operation.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecSerialExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( elements<CEAs...>( (~vector).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on the given vector transpose operation.
// \ingroup elements
//
// \param vector The constant vector transpose operation.
// \param args The runtime element arguments.
// \return View on the specified selection of elements on the transpose operation.
//
// This function returns an expression representing the specified selection of elements on the
// given vector transpose operation.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const VecTransExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return trans( elements<CEAs...>( (~vector).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another element selection.
// \ingroup elements
//
// \param e The given element selection.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
//
// This function returns an expression representing the specified selection of elements on the
// given element selection.
*/
template< size_t I1           // First required element index
        , size_t... Is1       // Remaining required element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t I2           // First present element index
        , size_t... Is2       // Remaining present element indices
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( Elements<VT,TF,DF,I2,Is2...>& e, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return elements< indices[I1], indices[Is1]... >( e.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another constant element selection.
// \ingroup elements
//
// \param e The given constant element selection.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
//
// This function returns an expression representing the specified selection of elements on the
// given constant element selection.
*/
template< size_t I1           // First required element index
        , size_t... Is1       // Remaining required element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t I2           // First present element index
        , size_t... Is2       // Remaining present element indices
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( const Elements<VT,TF,DF,I2,Is2...>& e, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return elements< indices[I1], indices[Is1]... >( e.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another temporary element selection.
// \ingroup elements
//
// \param e The given temporary element selection.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
//
// This function returns an expression representing the specified selection of elements on the
// given temporary element selection.
*/
template< size_t I1           // First required element index
        , size_t... Is1       // Remaining required element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t I2           // First present element index
        , size_t... Is2       // Remaining present element indices
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( Elements<VT,TF,DF,I2,Is2...>&& e, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return elements< indices[I1], indices[Is1]... >( e.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another element selection.
// \ingroup elements
//
// \param e The given element selection.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified selection of elements on the
// given element selection.
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( Elements<VT,TF,DF>& e, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( e.size() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
         }
      }
   }

   decltype(auto) indices( e.idces() );
   return elements( e.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another constant element selection.
// \ingroup elements
//
// \param e The given constant element selection.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified selection of elements on the
// given constant element selection.
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( const Elements<VT,TF,DF>& e, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( e.size() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
         }
      }
   }

   decltype(auto) indices( e.idces() );
   return elements( e.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another temporary element selection.
// \ingroup elements
//
// \param e The given temporary element selection.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified selection of elements on the
// given temporary element selection.
*/
template< size_t I            // First element index
        , size_t... Is        // Remaining element indices
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( Elements<VT,TF,DF>&& e, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( e.size() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
         }
      }
   }

   decltype(auto) indices( e.idces() );
   return elements( e.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another element selection.
// \ingroup elements
//
// \param e The given element selection.
// \param indices The container of element indices.
// \param n The total number of indices.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified selection of elements on the
// given element selection.
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( Elements<VT,TF,DF,CEAs...>& e,
                                const T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( e.size() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
         }
      }
   }

   decltype(auto) oldIndices( e.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return elements( e.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another constant element selection.
// \ingroup elements
//
// \param e The given constant element selection.
// \param indices The container of element indices.
// \param n The total number of indices.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified selection of elements on the
// given constant element selection.
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( const Elements<VT,TF,DF,CEAs...>& e,
                                const T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( e.size() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
         }
      }
   }

   decltype(auto) oldIndices( e.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return elements( e.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements on another temporary element selection.
// \ingroup elements
//
// \param e The given temporary element selection.
// \param indices The container of element indices.
// \param n The total number of indices.
// \param args The optional element arguments.
// \return View on the specified selection of elements on the other element selection.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified selection of elements on the
// given temporary element selection.
*/
template< typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename T          // Type of the element indices
        , typename... REAs >  // Optional element arguments
inline decltype(auto) elements( Elements<VT,TF,DF,CEAs...>&& e,
                                const T* indices, size_t n, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<REAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( e.size() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
         }
      }
   }

   decltype(auto) oldIndices( e.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return elements( e.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (SUBVECTOR)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given element selection.
// \ingroup elements
//
// \param e The selection of elements containing the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the element selection.
//
// This function returns an expression representing the specified subvector of the given element
// selection.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first subvector element
        , size_t N            // Size of the subvector
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename... RSAs >  // Optional arguments
inline decltype(auto) subvector( Elements<VT,TF,DF,CEAs...>& e, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( e, make_shifted_index_sequence<I,N>(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given constant element selection.
// \ingroup elements
//
// \param e The constant selection of elements containing the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the element selection.
//
// This function returns an expression representing the specified subvector of the given constant
// element selection.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first subvector element
        , size_t N            // Size of the subvector
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename... RSAs >  // Optional arguments
inline decltype(auto) subvector( const Elements<VT,TF,DF,CEAs...>& e, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( e, make_shifted_index_sequence<I,N>(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given temporary element selection.
// \ingroup elements
//
// \param e The temporary selection of elements containing the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the element selection.
//
// This function returns an expression representing the specified subvector of the given temporary
// element selection.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first subvector element
        , size_t N            // Size of the subvector
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename... RSAs >  // Optional arguments
inline decltype(auto) subvector( Elements<VT,TF,DF,CEAs...>&& e, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( e, make_shifted_index_sequence<I,N>(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given element selection.
// \ingroup elements
//
// \param e The selection of elements containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the element selection.
//
// This function returns an expression representing the specified subvector of the given element
// selection.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename... RSAs >  // Optional arguments
inline decltype(auto)
   subvector( Elements<VT,TF,DF,CEAs...>& e, size_t index, size_t size, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   SmallVector<size_t,128UL> indices( size );
   std::iota( indices.begin(), indices.end(), index );

   return elements( e, indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given constant element selection.
// \ingroup elements
//
// \param e The constant selection of elements containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the element selection.
//
// This function returns an expression representing the specified subvector of the given constant
// element selection.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename... RSAs >  // Optional arguments
inline decltype(auto)
   subvector( const Elements<VT,TF,DF,CEAs...>& e, size_t index, size_t size, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   SmallVector<size_t,128UL> indices( size );
   std::iota( indices.begin(), indices.end(), index );

   return elements( e, indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given temporary element selection.
// \ingroup elements
//
// \param e The temporary selection of elements containing the subvector.
// \param index The index of the first element of the subvector.
// \param size The size of the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the element selection.
//
// This function returns an expression representing the specified subvector of the given temporary
// element selection.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename VT         // Type of the vector
        , bool TF             // Transpose flag
        , bool DF             // Density flag
        , size_t... CEAs      // Compile time element arguments
        , typename... RSAs >  // Optional arguments
inline decltype(auto)
   subvector( Elements<VT,TF,DF,CEAs...>&& e, size_t index, size_t size, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   SmallVector<size_t,128UL> indices( size );
   std::iota( indices.begin(), indices.end(), index );

   return elements( std::move( e ), indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ELEMENTS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given selection of elements.
// \ingroup elements
//
// \param e The selection of elements to be resetted.
// \return void
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs >  // Compile time element arguments
inline void reset( Elements<VT,TF,DF,CEAs...>& e )
{
   e.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary selection of elements.
// \ingroup elements
//
// \param e The temporary selection of elements to be resetted.
// \return void
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs >  // Compile time element arguments
inline void reset( Elements<VT,TF,DF,CEAs...>&& e )
{
   e.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given selection of elements.
// \ingroup elements
//
// \param e The selection of elements to be cleared.
// \return void
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs >  // Compile time element arguments
inline void clear( Elements<VT,TF,DF,CEAs...>& e )
{
   e.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary selection of elements.
// \ingroup elements
//
// \param e The temporary selection of elements to be cleared.
// \return void
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs >  // Compile time element arguments
inline void clear( Elements<VT,TF,DF,CEAs...>&& e )
{
   e.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense element selection is in default state.
// \ingroup elements
//
// \param e The dense element selection to be tested for its default state.
// \return \a true in case the given element selection is component-wise zero, \a false otherwise.
//
// This function checks whether the dense element selection is in default state. For instance, in
// case the dense element selection is instantiated for a vector of built-in integral or floating
// point data type, the function returns \a true in case all elements are 0 and \a false in case
// any element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isDefault( elements( v, { 5UL, 10UL, 15UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( elements( v, { 5UL, 10UL, 15UL } ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CEAs >  // Compile time element arguments
inline bool isDefault( const Elements<VT,TF,true,CEAs...>& e )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<e.size(); ++i )
      if( !isDefault<RF>( e[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse element selection is in default state.
// \ingroup elements
//
// \param e The sparse element selection to be tested for its default state.
// \return \a true in case the given element selection is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse element selection is in default state. For instance, in
// case the sparse element selection is instantiated for a vector of built-in integral or floating
// point data type, the function returns \a true in case all elements are 0 and \a false in case
// any element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isDefault( elements( v, { 5UL, 10UL, 15UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( elements( v, { 5UL, 10UL, 15UL } ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename VT       // Type of the sparse vector
        , bool TF           // Transpose flag
        , size_t... CEAs >  // Compile time element arguments
inline bool isDefault( const Elements<VT,TF,false,CEAs...>& e )
{
   using blaze::isDefault;

   for( const auto& element : ~e )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given selection of elements are intact.
// \ingroup elements
//
// \param e The selection of elements to be tested.
// \return \a true in case the given selection's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the selection of elements are intact, i.e. if
// its state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicVector<int,rowVector> v;
   // ... Resizing and initialization
   if( isIntact( elements( v, { 5UL, 10UL, 15UL } ) ) ) { ... }
   \endcode
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs >  // Compile time element arguments
inline bool isIntact( const Elements<VT,TF,DF,CEAs...>& e ) noexcept
{
   return ( e.size() <= e.operand().size() && isIntact( e.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given element selection and vector represent the same observable state.
// \ingroup elements
//
// \param a The element selection to be tested for its state.
// \param b The vector to be tested for its state.
// \return \a true in case the element selection and vector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given selection of elements refers to the
// entire range of the given vector in ascending and consecutive order and by that represents
// the same observable state. In this case, the function returns \a true, otherwise it returns
// \a false.
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs >  // Compile time element arguments
inline bool isSame( const Elements<VT,TF,DF,CEAs...>& a, const Vector<VT,TF>& b ) noexcept
{
   if( !isSame( a.operand(), ~b ) || ( a.size() != (~b).size() ) )
      return false;

   decltype(auto) indices( a.idces() );
   for( size_t i=0UL; i<a.size(); ++i ) {
      if( indices[i] != i )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given vector and element selection represent the same observable state.
// \ingroup elements
//
// \param a The vector to be tested for its state.
// \param b The element selection to be tested for its state.
// \return \a true in case the vector and element selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given selection of elements refers to the
// entire range of the given vector in ascending and consecutive order and by that represents
// the same observable state. In this case, the function returns \a true, otherwise it returns
// \a false.
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs >  // Compile time element arguments
inline bool isSame( const Vector<VT,TF>& a, const Elements<VT,TF,DF,CEAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given element selection and subvector represent the same observable state.
// \ingroup elements
//
// \param a The element selection to be tested for its state.
// \param b The subvector to be tested for its state.
// \return \a true in case the element selection and vector share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given selection of elements refers to the
// entire range of the given subvector in ascending and consecutive order and by that represents
// the same observable state. In this case, the function returns \a true, otherwise it returns
// \a false.
*/
template< typename VT1      // Type of the left-hand side vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename VT2      // Type of the right-hand side vector
        , AlignmentFlag AF  // Alignment flag
        , size_t... CSAs >  // Compile time subvector arguments
inline bool isSame( const Elements<VT1,TF,DF,CEAs...>& a, const Subvector<VT2,AF,TF,DF,CSAs...>& b ) noexcept
{
   if( !isSame( a.operand(), b.operand() ) || ( a.size() != b.size() ) )
      return false;

   decltype(auto) indices( a.idces() );
   for( size_t i=0UL; i<a.size(); ++i ) {
      if( indices[i] != b.offset()+i )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given subvector and element selection represent the same observable state.
// \ingroup elements
//
// \param a The subvector to be tested for its state.
// \param b The selection of elements to be tested for its state.
// \return \a true in case the vector and element selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given selection of elements refers to the
// entire range of the given subvector in ascending and consecutive order and by that represents
// the same observable state. In this case, the function returns \a true, otherwise it returns
// \a false.
*/
template< typename VT1      // Type of the left-hand side vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time subvector arguments
        , typename VT2      // Type of the right-hand side vector
        , size_t... CEAs >  // Compile time element arguments
inline bool isSame( const Subvector<VT1,AF,TF,DF,CSAs...>& a, const Elements<VT2,TF,DF,CEAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given element selections represent the same observable state.
// \ingroup elements
//
// \param a The first selection of elements to be tested for its state.
// \param b The second selection of elements to be tested for its state.
// \return \a true in case the two element selections share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given element selections refer to exactly
// the same range of the same vector. In case both selections represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename VT1       // Type of the vector of the left-hand side selection of elements
        , bool TF1           // Transpose flag of the left-hand side selection of elements
        , bool DF1           // Density flag of the left-hand side selection of elements
        , size_t... CEAs1    // Compile time element arguments of the left-hand side selection of elements
        , typename VT2       // Type of the vector of the right-hand side selection of elements
        , bool TF2           // Transpose flag of the right-hand side selection of elements
        , bool DF2           // Density flag of the right-hand side selection of elements
        , size_t... CEAs2 >  // Compile time element arguments of the right-hand side selection of elements
inline bool isSame( const Elements<VT1,TF1,DF1,CEAs1...>& a,
                    const Elements<VT2,TF2,DF2,CEAs2...>& b ) noexcept
{
   if( !isSame( a.operand(), b.operand() ) || a.size() != b.size() )
      return false;

   decltype(auto) indices1( a.idces() );
   decltype(auto) indices2( b.idces() );

   return std::equal( indices1.begin(), indices1.end(), indices2.begin() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a selection of elements.
// \ingroup elements
//
// \param e The target selection of elements.
// \param index The index of the element to be set.
// \param value The value to be set to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename ET >     // Type of the element
inline bool trySet( const Elements<VT,TF,DF,CEAs...>& e, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < e.size(), "Invalid vector access index" );

   return trySet( e.operand(), e.idx(index), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a selection of elements.
// \ingroup elements
//
// \param e The target selection of elements.
// \param index The index of the element to be modified.
// \param value The value to be added to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename ET >     // Type of the element
inline bool tryAdd( const Elements<VT,TF,DF,CEAs...>& e, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < e.size(), "Invalid vector access index" );

   return tryAdd( e.operand(), e.idx(index), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a selection of elements.
// \ingroup elements
//
// \param e The target selection of elements.
// \param index The index of the element to be modified.
// \param value The value to be subtracted from the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename ET >     // Type of the element
inline bool trySub( const Elements<VT,TF,DF,CEAs...>& e, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < e.size(), "Invalid vector access index" );

   return trySub( e.operand(), e.idx(index), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a selection of elements.
// \ingroup elements
//
// \param e The target selection of elements.
// \param index The index of the element to be modified.
// \param value The factor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename ET >     // Type of the element
inline bool tryMult( const Elements<VT,TF,DF,CEAs...>& e, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < e.size(), "Invalid vector access index" );

   return tryMult( e.operand(), e.idx(index), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a selection of elements.
// \ingroup elements
//
// \param e The target selection of elements.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The factor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT     // Type of the vector
        , bool TF         // Transpose flag
        , bool DF         // Density flag
        , size_t... CEAs  // Compile time element arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Elements<VT,TF,DF,CEAs...>& e, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (~e).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (~e).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryMult( e.operand(), e.idx(i), value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a selection of elements.
// \ingroup elements
//
// \param e The target selection of elements.
// \param index The index of the element to be modified.
// \param value The divisor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename ET >     // Type of the element
inline bool tryDiv( const Elements<VT,TF,DF,CEAs...>& e, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < e.size(), "Invalid vector access index" );

   return tryDiv( e.operand(), e.idx(index), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a selection of elements.
// \ingroup elements
//
// \param e The target selection of elements.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The divisor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT     // Type of the vector
        , bool TF         // Transpose flag
        , bool DF         // Density flag
        , size_t... CEAs  // Compile time element arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Elements<VT,TF,DF,CEAs...>& e, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (~e).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (~e).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryDiv( e.operand(), e.idx(i), value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a selection of elements.
// \ingroup elements
//
// \param lhs The target left-hand side selection of elements.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1      // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename VT2 >    // Type of the right-hand side vector
inline bool tryAssign( const Elements<VT1,TF,DF,CEAs...>& lhs,
                       const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (~rhs).size() <= lhs.size(), "Invalid vector size" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySet( lhs.operand(), lhs.idx(i+index), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a selection of
//        elements.
// \ingroup elements
//
// \param lhs The target left-hand side selection of elements.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1      // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename VT2 >    // Type of the right-hand side vector
inline bool tryAddAssign( const Elements<VT1,TF,DF,CEAs...>& lhs,
                          const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (~rhs).size() <= lhs.size(), "Invalid vector size" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryAdd( lhs.operand(), lhs.idx(i+index), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a selection
//        of elements.
// \ingroup elements
//
// \param lhs The target left-hand side selection of elements.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1      // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename VT2 >    // Type of the right-hand side vector
inline bool trySubAssign( const Elements<VT1,TF,DF,CEAs...>& lhs,
                          const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (~rhs).size() <= lhs.size(), "Invalid vector size" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySub( lhs.operand(), lhs.idx(i+index), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a selection
//        of elements.
// \ingroup elements
//
// \param lhs The target left-hand side selection of elements.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1      // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename VT2 >    // Type of the right-hand side vector
inline bool tryMultAssign( const Elements<VT1,TF,DF,CEAs...>& lhs,
                           const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (~rhs).size() <= lhs.size(), "Invalid vector size" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryMult( lhs.operand(), lhs.idx(i+index), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a selection of
//        elements.
// \ingroup elements
//
// \param lhs The target left-hand side selection of elements.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1      // Type of the vector
        , bool TF           // Transpose flag
        , bool DF           // Density flag
        , size_t... CEAs    // Compile time element arguments
        , typename VT2 >    // Type of the right-hand side vector
inline bool tryDivAssign( const Elements<VT1,TF,DF,CEAs...>& lhs,
                          const Vector<VT2,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (~rhs).size() <= lhs.size(), "Invalid vector size" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryDiv( lhs.operand(), lhs.idx(i+index), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given element selection.
// \ingroup elements
//
// \param e The element selection to be derestricted.
// \return Element selection without access restrictions.
//
// This function removes all restrictions on the data access to the given element selection.
// It returns an element selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT     // Type of the vector
        , bool TF         // Transpose flag
        , bool DF         // Density flag
        , size_t I        // First element index
        , size_t... Is >  // Remaining element indices
inline decltype(auto) derestrict( Elements<VT,TF,DF,I,Is...>& e )
{
   return elements<I,Is...>( derestrict( e.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary element selection.
// \ingroup elements
//
// \param e The temporary element selection to be derestricted.
// \return Element selection without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary element
// selection. It returns an element selection that does provide the same interface but does
// not have any restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT     // Type of the vector
        , bool TF         // Transpose flag
        , bool DF         // Density flag
        , size_t I        // First element index
        , size_t... Is >  // Remaining element indices
inline decltype(auto) derestrict( Elements<VT,TF,DF,I,Is...>&& e )
{
   return elements<I,Is...>( derestrict( e.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given element selection.
// \ingroup elements
//
// \param e The element selection to be derestricted.
// \return Element selection without access restrictions.
//
// This function removes all restrictions on the data access to the given element selection.
// It returns an element selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the vector
        , bool TF      // Transpose flag
        , bool DF >    // Density flag
inline decltype(auto) derestrict( Elements<VT,TF,DF>& e )
{
   decltype(auto) indices( e.idces() );
   return elements( derestrict( e.operand() ), indices.data(), indices.size(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary element selection.
// \ingroup elements
//
// \param e The temporary element selection to be derestricted.
// \return Element selection without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary element
// selection. It returns an element selection that does provide the same interface but does not
// have any restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename VT  // Type of the vector
        , bool TF      // Transpose flag
        , bool DF >    // Density flag
inline decltype(auto) derestrict( Elements<VT,TF,DF>&& e )
{
   decltype(auto) indices( e.idces() );
   return elements( derestrict( e.operand() ), indices.data(), indices.size(), unchecked );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF, bool DF, size_t I, size_t... Is >
struct Size< Elements<VT,TF,DF,I,Is...>, 0UL >
   : public PtrdiffT<1UL+sizeof...(Is)>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAXSIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF, bool DF, size_t I, size_t... Is >
struct MaxSize< Elements<VT,TF,DF,I,Is...>, 0UL >
   : public PtrdiffT<1UL+sizeof...(Is)>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESTRICTED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF, bool DF, size_t... CEAs >
struct IsRestricted< Elements<VT,TF,DF,CEAs...> >
   : public IsRestricted<VT>
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
template< typename VT, bool TF, size_t... CEAs >
struct HasConstDataAccess< Elements<VT,TF,true,CEAs...> >
   : public HasConstDataAccess<VT>
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
template< typename VT, bool TF, size_t... CEAs >
struct HasMutableDataAccess< Elements<VT,TF,true,CEAs...> >
   : public HasMutableDataAccess<VT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ELEMENTSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, bool TF, bool DF, size_t... CEAs1, size_t... CEAs2 >
struct ElementsTrait< Elements<VT,TF,DF,CEAs1...>, CEAs2... >
{
   using Type = ElementsTrait_t< ResultType_t< Elements<VT,TF,DF,CEAs1...> >, CEAs2... >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
