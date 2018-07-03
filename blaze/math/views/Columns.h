//=================================================================================================
/*!
//  \file blaze/math/views/Columns.h
//  \brief Header file for the implementation of the Columns view
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

#ifndef _BLAZE_MATH_VIEWS_COLUMNS_H_
#define _BLAZE_MATH_VIEWS_COLUMNS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <array>
#include <vector>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DeclExpr.h>
#include <blaze/math/expressions/MatEvalExpr.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/expressions/MatMatMapExpr.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/expressions/MatMatSubExpr.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/expressions/MatSerialExpr.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/TVecMatMultExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/IntegerSequence.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/ColumnsTrait.h>
#include <blaze/math/traits/RowsTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/Forward.h>
#include <blaze/math/views/column/ColumnData.h>
#include <blaze/math/views/columns/BaseTemplate.h>
#include <blaze/math/views/columns/Dense.h>
#include <blaze/math/views/columns/Sparse.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/PtrdiffT.h>
#include <blaze/util/SmallVector.h>
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
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   auto columns2 = columns<4UL,2UL>( S );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns<1UL,3UL>( D, unchecked );
   auto columns2 = columns<4UL,2UL>( S, unchecked );
   \endcode
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT,I,Is...>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given constant matrix.
// \ingroup columns
//
// \param matrix The constant matrix containing the columns.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   const blaze::DynamicMatrix<double,columnMajor> D( ... );
   const blaze::CompressedMatrix<double,columnMajor> S( ... );

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   auto columns2 = columns<4UL,2UL>( S );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns<1UL,3UL>( D, unchecked );
   auto columns2 = columns<4UL,2UL>( S, unchecked );
   \endcode
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const Matrix<MT,SO>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Columns_<const MT,I,Is...>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given temporary matrix.
// \ingroup columns
//
// \param matrix The temporary matrix containing the columns.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given temporary
// matrix. In case any column is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>&& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT,I,Is...>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto columns1 = columns( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   const std::array<size_t,2uL> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, indices1.data(), indices1.size(), unchecked );
   auto columns2 = columns( S, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>& matrix, const T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT>;
   return ReturnType( ~matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given constant matrix.
// \ingroup columns
//
// \param matrix The constant matrix containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto columns1 = columns( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   const std::array<size_t,2uL> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, indices1.data(), indices1.size(), unchecked );
   auto columns2 = columns( S, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const Matrix<MT,SO>& matrix, const T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Columns_<const MT>;
   return ReturnType( ~matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given temporary matrix.
// \ingroup columns
//
// \param matrix The temporary matrix containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given temporary
// matrix. In case any column is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>&& matrix, const T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT>;
   return ReturnType( ~matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The sequence of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.
// In case any column is not properly specified (i.e. if any specified index is greater than or
// equal to the total number of columns in the given matrix) a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT         // Type of the matrix
        , size_t... Is        // Column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, index_sequence<Is...> indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   UNUSED_PARAMETER( indices );

   return columns<Is...>( std::forward<MT>( matrix ), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The list of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.
// In case any column is not properly specified (i.e. if any specified index is greater than or
// equal to the total number of columns in the given matrix) a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, initializer_list<T> indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.begin(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The array of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.
// In case any column is not properly specified (i.e. if any specified index is greater than or
// equal to the total number of columns in the given matrix) a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , size_t N            // Number of indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, const std::array<T,N>& indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.data(), N, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The vector of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.
// In case any column is not properly specified (i.e. if any specified index is greater than or
// equal to the total number of columns in the given matrix) a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, const std::vector<T>& indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The vector of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.
// In case any column is not properly specified (i.e. if any specified index is greater than or
// equal to the total number of columns in the given matrix) a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , size_t N            // Number of preallocated elements
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, const SmallVector<T,N>& indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.data(), indices.size(), args... );
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
/*!\brief Creating a view on a selection of columns on the given matrix/matrix addition.
// \ingroup columns
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the addition.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/matrix addition.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatMatAddExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (~matrix).leftOperand(), args... ) +
          columns<CCAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/matrix subtraction.
// \ingroup columns
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the subtraction.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/matrix subtraction.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatMatSubExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (~matrix).leftOperand(), args... ) -
          columns<CCAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given Schur product.
// \ingroup columns
//
// \param matrix The constant Schur product.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the subtraction.
//
// This function returns an expression representing the specified selection of columns on the
// given Schur product.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const SchurExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (~matrix).leftOperand(), args... ) %
          columns<CCAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/matrix multiplication.
// \ingroup columns
//
// \param matrix The constant matrix/matrix multiplication.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/matrix multiplication.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatMatMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * columns<CCAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given outer product.
// \ingroup columns
//
// \param matrix The constant outer product.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the outer product.
//
// This function returns an expression representing the specified selection of columns on the
// given outer product.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const VecTVecMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * elements<CCAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/scalar multiplication.
// \ingroup columns
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/scalar multiplication.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatScalarMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (~matrix).leftOperand(), args... ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/scalar division.
// \ingroup columns
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the division.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/scalar division.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatScalarDivExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (~matrix).leftOperand(), args... ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given unary matrix map operation.
// \ingroup columns
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given unary matrix map operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( columns<CCAs...>( (~matrix).operand(), args... ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given binary matrix map operation.
// \ingroup columns
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given binary matrix map operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatMatMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( columns<CCAs...>( (~matrix).leftOperand(), args... ),
               columns<CCAs...>( (~matrix).rightOperand(), args... ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix evaluation operation.
// \ingroup columns
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix evaluation operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatEvalExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( columns<CCAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix serialization operation.
// \ingroup columns
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix serialization operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatSerialExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( columns<CCAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix declaration operation.
// \ingroup columns
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix declaration operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const DeclExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (~matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix transpose operation.
// \ingroup columns
//
// \param matrix The constant matrix transpose operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix transpose operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , typename = EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) > >
inline decltype(auto) columns( const MatTransExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return trans( rows<CCAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given column selection.
// \ingroup columns
//
// \param c The selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
//
// This function returns an expression representing the specified columns of the given column
// selection.
*/
template< size_t I1           // First required column index
        , size_t... Is1       // Remaining required column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First present column index
        , size_t... Is2       // Remaining present column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( Columns<MT,SO,DF,SF,I2,Is2...>& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return columns< indices[I1], indices[Is1]... >( c.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given constant column selection.
// \ingroup columns
//
// \param c The constant selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
//
// This function returns an expression representing the specified columns of the given constant
// column selection.
*/
template< size_t I1           // First required column index
        , size_t... Is1       // Remaining required column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First present column index
        , size_t... Is2       // Remaining present column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const Columns<MT,SO,DF,SF,I2,Is2...>& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return columns< indices[I1], indices[Is1]... >( c.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given temporary column selection.
// \ingroup columns
//
// \param c The temporary selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
//
// This function returns an expression representing the specified columns of the given temporary
// column selection.
*/
template< size_t I1           // First required column index
        , size_t... Is1       // Remaining required column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First present column index
        , size_t... Is2       // Remaining present column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( Columns<MT,SO,DF,SF,I2,Is2...>&& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return columns< indices[I1], indices[Is1]... >( c.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given column selection.
// \ingroup columns
//
// \param c The selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified columns of the given column
// selection.
*/
template< size_t I            // First required column index
        , size_t... Is        // Remaining required column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs      // Compile time column arguments
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( Columns<MT,SO,DF,SF,CCAs...>& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( c.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   decltype(auto) indices( c.idces() );
   return columns( c.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given constant column selection.
// \ingroup columns
//
// \param c The constant selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified columns of the given constant
// column selection.
*/
template< size_t I            // First required column index
        , size_t... Is        // Remaining required column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs      // Compile time column arguments
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const Columns<MT,SO,DF,SF,CCAs...>& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( c.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   decltype(auto) indices( c.idces() );
   return columns( c.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given temporary column selection.
// \ingroup columns
//
// \param c The temporary selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified columns of the given temporary
// column selection.
*/
template< size_t I            // First required column index
        , size_t... Is        // Remaining required column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs      // Compile time column arguments
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( Columns<MT,SO,DF,SF,CCAs...>&& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( c.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   decltype(auto) indices( c.idces() );
   return columns( c.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given column selection.
// \ingroup columns
//
// \param c The selection of columns containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given column
// selection.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs      // Compile time column arguments
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( Columns<MT,SO,DF,SF,CCAs...>& c, const T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( c.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   decltype(auto) oldIndices( c.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return columns( c.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given constant column selection.
// \ingroup columns
//
// \param c The constant selection of columns containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given constant
// column selection.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs      // Compile time column arguments
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const Columns<MT,SO,DF,SF,CCAs...>& c, const T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( c.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   decltype(auto) oldIndices( c.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return columns( c.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given temporary column selection.
// \ingroup columns
//
// \param c The temporary selection of columns containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given temporary
// column selection.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs      // Compile time column arguments
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( Columns<MT,SO,DF,SF,CCAs...>&& c, const T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( c.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   decltype(auto) oldIndices( c.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return columns( c.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ELEMENTS)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given vector/matrix multiplication.
// \ingroup columns
//
// \param vector The constant vector/matrix multiplication.
// \param args The runtime element arguments.
// \return View on the specified elements of the multiplication.
//
// This function returns an expression representing the specified elements of the given
// vector/matrix multiplication.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const TVecMatMultExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return (~vector).leftOperand() * columns<CEAs...>( (~vector).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ROW)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given column selection.
// \ingroup columns
//
// \param columns The selection of columns containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the column selection.
//
// This function returns an expression representing the specified row of the given column selection.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Columns<MT,SO,DF,SF,I,Is...>& columns, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<I,Is...>( row<CRAs...>( columns.operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given constant column selection.
// \ingroup columns
//
// \param columns The constant selection of columns containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the column selection.
//
// This function returns an expression representing the specified row of the given constant
// column selection.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const Columns<MT,SO,DF,SF,I,Is...>& columns, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<I,Is...>( row<CRAs...>( columns.operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given temporary column selection.
// \ingroup columns
//
// \param columns The temporary selection of columns containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the column selection.
//
// This function returns an expression representing the specified row of the given temporary
// column selection.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Columns<MT,SO,DF,SF,I,Is...>&& columns, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<I,Is...>( row<CRAs...>( columns.operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given column selection.
// \ingroup columns
//
// \param columns The selection of columns containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the column selection.
//
// This function returns an expression representing the specified row of the given column selection.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Columns<MT,SO,DF,SF>& columns, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( row<CRAs...>( columns.operand(), args... ), columns.idces() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given constant column selection.
// \ingroup columns
//
// \param columns The constant selection of columns containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the column selection.
//
// This function returns an expression representing the specified row of the given constant
// column selection.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const Columns<MT,SO,DF,SF>& columns, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( row<CRAs...>( columns.operand(), args... ), columns.idces() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given temporary column selection.
// \ingroup columns
//
// \param columns The temporary selection of columns containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the column selection.
//
// This function returns an expression representing the specified row of the given temporary
// column selection.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Columns<MT,SO,DF,SF>&& columns, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( row<CRAs...>( columns.operand(), args... ), columns.idces() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (COLUMN)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given column selection.
// \ingroup columns
//
// \param columns The selection of columns containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the column selection.
//
// This function returns an expression representing the specified column of the given column
// selection.
*/
template< size_t I1           // Column index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First column index
        , size_t... Is        // Remaining column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Columns<MT,SO,DF,SF,I2,Is...>& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is... };
   return column<indices[I1]>( columns.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given constant column selection.
// \ingroup columns
//
// \param columns The constant selection of columns containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the column selection.
//
// This function returns an expression representing the specified column of the given constant
// column selection.
*/
template< size_t I1           // Column index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First column index
        , size_t... Is        // Remaining column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const Columns<MT,SO,DF,SF,I2,Is...>& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is... };
   return column<indices[I1]>( columns.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given temporary column selection.
// \ingroup columns
//
// \param columns The temporary selection of columns containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the column selection.
//
// This function returns an expression representing the specified column of the given temporary
// column selection.
*/
template< size_t I1           // Column index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First column index
        , size_t... Is        // Remaining column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Columns<MT,SO,DF,SF,I2,Is...>&& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is... };
   return column<indices[I1]>( columns.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given column selection.
// \ingroup columns
//
// \param columns The selection of columns containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given column
// selection.
*/
template< size_t... CCAs1     // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs2     // Compile time column arguments
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Columns<MT,SO,DF,SF,CCAs2...>& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const ColumnData<CCAs1...> cd( args... );
   decltype(auto) indices( columns.idces() );

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      if( indices.size() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return column( columns.operand(), indices[cd.column()], args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given constant column selection.
// \ingroup columns
//
// \param columns The constant selection of columns containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given constant
// column selection.
*/
template< size_t... CCAs1     // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs2     // Compile time column arguments
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const Columns<MT,SO,DF,SF,CCAs2...>& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const ColumnData<CCAs1...> cd( args... );
   decltype(auto) indices( columns.idces() );

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      if( indices.size() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return column( columns.operand(), indices[cd.column()], args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given temporary column selection.
// \ingroup columns
//
// \param columns The temporary selection of columns containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given temporary
// column selection.
*/
template< size_t... CCAs1     // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CCAs2     // Compile time column arguments
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Columns<MT,SO,DF,SF,CCAs2...>&& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const ColumnData<CCAs1...> cd( args... );
   decltype(auto) indices( columns.idces() );

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      if( indices.size() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return column( columns.operand(), indices[cd.column()], args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given column selection.
// \ingroup columns
//
// \param columns The column selection to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void reset( Columns<MT,SO,DF,SF,CCAs...>& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary column selection.
// \ingroup columns
//
// \param columns The temporary column selection to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void reset( Columns<MT,SO,DF,SF,CCAs...>&& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column of the given column selection.
// \ingroup columns
//
// \param columns The column selection to be resetted.
// \param i The index of the column to be resetted.
// \return void
//
// This function resets the values in the specified column of the given column selection to their
// default value. Note that the capacity of the column remains unchanged.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void reset( Columns<MT,SO,DF,SF,CCAs...>& columns, size_t i )
{
   columns.reset( i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given column selection.
// \ingroup columns
//
// \param columns The column selection to be cleared.
// \return void
//
// Clearing a column selection is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void clear( Columns<MT,SO,DF,SF,CCAs...>& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary column selection.
// \ingroup columns
//
// \param columns The column selection to be cleared.
// \return void
//
// Clearing a column selection is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void clear( Columns<MT,SO,DF,SF,CCAs...>&& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense column selection is in default state.
// \ingroup columns
//
// \param columns The dense column selection to be tested for its default state.
// \return \a true in case the given dense column selection is component-wise zero, \a false otherwise.
//
// This function checks whether the dense column selection is in default state. For instance, in
// case the column selection is instantiated for a built-in integral or floating point data type,
// the function returns \a true in case all column elements are 0 and \a false in case any column
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool isDefault( const Columns<MT,SO,true,SF,CCAs...>& columns )
{
   using blaze::isDefault;

   if( SO == false ) {
      for( size_t i=0UL; i<columns.rows(); ++i )
         for( size_t j=0UL; j<columns.columns(); ++j )
            if( !isDefault<RF>( columns(i,j) ) )
               return false;
   }
   else {
      for( size_t j=0UL; j<columns.columns(); ++j )
         for( size_t i=0UL; i<columns.rows(); ++i )
            if( !isDefault<RF>( columns(i,j) ) )
               return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse column selection is in default state.
// \ingroup columns
//
// \param columns The sparse column selection to be tested for its default state.
// \return \a true in case the given sparse column selection is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse column selection is in default state. For instance, in
// case the column selection is instantiated for a built-in integral or floating point data type,
// the function returns \a true in case all column elements are 0 and \a false in case any column
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool isDefault( const Columns<MT,SO,false,SF,CCAs...>& columns )
{
   using blaze::isDefault;

   for( size_t j=0UL; j<columns.columns(); ++j ) {
      for( auto element=columns.cbegin(j); element!=columns.cend(j); ++element )
         if( !isDefault<RF>( element->value() ) ) return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given column selection are intact.
// \ingroup columns
//
// \param columns The column selection to be tested.
// \return \a true in case the given column selection's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the column selection are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isIntact( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool isIntact( const Columns<MT,SO,DF,SF,CCAs...>& columns ) noexcept
{
   return ( columns.rows() == columns.operand().rows() &&
            columns.columns() <= columns.operand().columns() &&
            isIntact( columns.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given column selection and matrix represent the same observable state.
// \ingroup columns
//
// \param a The column selection to be tested for its state.
// \param b The matrix to be tested for its state.
// \return \a true in case the column selection and matrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to all columns
// of the given matrix in ascending and consecutive order and by that represents the same observable
// state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT     // Type of the matrix
        , bool SO1        // Storage order of the left-hand side column selection
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool isSame( const Columns<MT,SO1,DF,SF,CCAs...>& a, const Matrix<MT,SO2>& b ) noexcept
{
   if( !isSame( a.operand(), ~b ) || ( a.rows() != (~b).rows() ) || ( a.columns() != (~b).columns() ) )
      return false;

   decltype(auto) indices( a.idces() );
   for( size_t j=0UL; j<a.columns(); ++j ) {
      if( indices[j] != j )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given matrix and column selection represent the same observable state.
// \ingroup columns
//
// \param a The matrix to be tested for its state.
// \param b The column selection to be tested for its state.
// \return \a true in case the matrix and column selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to all columns
// of the given matrix in ascending and consecutive order and by that represents the same observable
// state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO1          // Storage order of the left-hand side matrix
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs    // Compile time column arguments
        , bool SO2 >        // Storage order of the right-hand side column selection
inline bool isSame( const Matrix<MT,SO1>& a, const Columns<MT,SO2,DF,SF,CCAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given column selection and submatrix represent the same observable state.
// \ingroup columns
//
// \param a The column selection to be tested for its state.
// \param b The submatrix to be tested for its state.
// \return \a true in case the column selection and submatrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to same columns
// as the given submatrix in ascending and consecutive order and by that represents the same
// observable state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO1          // Storage order of the left-hand side column selection
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs    // Compile time column arguments
        , AlignmentFlag AF  // Alignment flag
        , bool SO2          // Storage order of the right-hand side submatrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isSame( const Columns<MT,SO1,DF,SF,CCAs...>& a, const Submatrix<MT,AF,SO2,DF,CSAs...>& b ) noexcept
{
   if( !isSame( a.operand(), b.operand() ) || ( a.rows() != (~b).rows() ) || ( a.columns() != (~b).columns() ) )
      return false;

   decltype(auto) indices( a.idces() );
   for( size_t j=0UL; j<a.columns(); ++j ) {
      if( indices[j] != b.column()+j )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given submatrix and column selection represent the same observable state.
// \ingroup columns
//
// \param a The submatrix to be tested for its state.
// \param b The column selection to be tested for its state.
// \return \a true in case the submatrix and column selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to same columns
// as the given submatrix in ascending and consecutive order and by that represents the same
// observable state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO1          // Storage order of the left-hand side submatrix
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , bool SO2          // Storage order of the right-hand side column selection
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool isSame( const Submatrix<MT,AF,SO1,DF,CSAs...>& a, const Columns<MT,SO2,DF,SF,CCAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given column selections represent the same observable state.
// \ingroup columns
//
// \param a The first selection of columns to be tested for its state.
// \param b The second selection of columns to be tested for its state.
// \return \a true in case the two column selections share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given column selections refer to exactly
// the same range of the same matrix. In case both selections represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT        // Type of the matrix
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF            // Symmetry flag
        , size_t... CCAs1    // Compile time column arguments of the left-hand side column selection
        , size_t... CCAs2 >  // Compile time column arguments of the right-hand side column selection
inline bool isSame( const Columns<MT,SO,DF,SF,CCAs1...>& a,
                    const Columns<MT,SO,DF,SF,CCAs2...>& b ) noexcept
{
   if( !isSame( a.operand(), b.operand() ) || a.rows() != b.rows() || a.columns() != b.columns() )
      return false;

   decltype(auto) indices1( a.idces() );
   decltype(auto) indices2( b.idces() );

   return std::equal( indices1.begin(), indices1.end(), indices2.begin() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense column selection.
// \ingroup columns
//
// \param c The dense column selection to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given dense column selection by means of the specified matrix type
// or matrix inversion algorithm \c IF (see the InversionFlag documentation):

   \code
   invert<asLower>( A );     // Inversion of a lower triangular matrix
   invert<asUniUpper>( A );  // Inversion of an upper unitriangular matrix
   invert<byLU>( A );        // Inversion by means of an LU decomposition
   invert<byLLH>( A );       // Inversion by means of a Cholesky decomposition
   ...
   \endcode

// The matrix inversion fails if ...
//
//  - ... the given column selection is not a square matrix;
//  - ... the given column selection is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a c may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline DisableIf_t< HasMutableDataAccess_v<MT> > invert( Columns<MT,SO,true,SF,CCAs...>& c )
{
   using RT = ResultType_t< Columns<MT,SO,true,SF,CCAs...> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION  ( RT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( RT );

   RT tmp( c );
   invert<IF>( tmp );
   c = tmp;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be set.
// \param j The column index of the element to be set.
// \param value The value to be set to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool trySet( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return trySet( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be added to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryAdd( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryAdd( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be subtracted from the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool trySub( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return trySub( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The factor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryMult( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryMult( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The factor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (~c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (~c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (~c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (~c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryMult( c.operand(), row, c.idx(j), m, n, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The divisor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryDiv( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryDiv( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The divisor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (~c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (~c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (~c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (~c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryDiv( c.operand(), row, c.idx(j), m, n, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a column vector to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                       const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryAssign( lhs.operand(), ~rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a row vector to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                       const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySet( lhs.operand(), row, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to the band of a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be assigned.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                       const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySet( lhs.operand(), row+i, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a matrix to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the matrix
        , bool SO1        // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool tryAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                       const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(~rhs).columns(); ++j ) {
      if( !tryAssign( lhs.operand(), blaze::column( ~rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a column vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAddAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryAddAssign( lhs.operand(), ~rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAddAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryAdd( lhs.operand(), row, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to the band of a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be added.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryAdd( lhs.operand(), row+i, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the matrix
        , bool SO1        // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool tryAddAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(~rhs).columns(); ++j ) {
      if( !tryAddAssign( lhs.operand(), blaze::column( ~rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a column vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool trySubAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return trySubAssign( lhs.operand(), ~rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool trySubAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySub( lhs.operand(), row, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to the band of a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be subtracted.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool trySubAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySub( lhs.operand(), row+i, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the matrix
        , bool SO1        // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool trySubAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(~rhs).columns(); ++j ) {
      if( !trySubAssign( lhs.operand(), blaze::column( ~rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a column vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be multiplied.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryMultAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                           const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryMultAssign( lhs.operand(), ~rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a row vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be multiplied.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryMultAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                           const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryMult( lhs.operand(), row, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to the band
//        of a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be multiplied.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryMultAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                           const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryMult( lhs.operand(), row+i, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the Schur product assignment of a matrix to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix for the Schur product.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the matrix
        , bool SO1        // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool trySchurAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(~rhs).columns(); ++j ) {
      if( !tryMultAssign( lhs.operand(), blaze::column( ~rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a column vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector divisor.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryDivAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryDivAssign( lhs.operand(), ~rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector divisor.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryDivAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryDiv( lhs.operand(), row, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to the band of a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector divisor.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryDivAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryDiv( lhs.operand(), row+i, lhs.idx( column+i ), (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given column selection.
// \ingroup columns
//
// \param c The column selection to be derestricted.
// \return Column selection without access restrictions.
//
// This function removes all restrictions on the data access to the given column selection.
// It returns a column selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t I        // First element index
        , size_t... Is >  // Remaining element indices
inline decltype(auto) derestrict( Columns<MT,SO,DF,SF,I,Is...>& c )
{
   return columns<I,Is...>( derestrict( c.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary column selection.
// \ingroup columns
//
// \param c The temporary column selection to be derestricted.
// \return Column selection without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary column
// selection. It returns a column selection that does provide the same interface but does not
// have any restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t I        // First element index
        , size_t... Is >  // Remaining element indices
inline decltype(auto) derestrict( Columns<MT,SO,DF,SF,I,Is...>&& c )
{
   return columns<I,Is...>( derestrict( c.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given column selection.
// \ingroup columns
//
// \param c The column selection to be derestricted.
// \return Column selection without access restrictions.
//
// This function removes all restrictions on the data access to the given column selection.
// It returns a column selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Columns<MT,SO,DF,SF>& c )
{
   decltype(auto) indices( c.idces() );
   return columns( derestrict( c.operand() ), indices.data(), indices.size(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary column selection.
// \ingroup columns
//
// \param r The temporary column selection to be derestricted.
// \return Column selection without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary column selection.
// It returns a column selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Columns<MT,SO,DF,SF>&& r )
{
   decltype(auto) indices( r.idces() );
   return columns( derestrict( r.operand() ), indices.data(), indices.size(), unchecked );
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs >
struct Size< Columns<MT,SO,DF,SF,CCAs...>, 0UL >
   : public Size<MT,0UL>
{};

template< typename MT, bool SO, bool DF, bool SF, size_t I, size_t... Is >
struct Size< Columns<MT,SO,DF,SF,I,Is...>, 1UL >
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs >
struct MaxSize< Columns<MT,SO,DF,SF,CCAs...>, 0UL >
   : public MaxSize<MT,0UL>
{};

template< typename MT, bool SO, bool DF, bool SF, size_t I, size_t... Is >
struct MaxSize< Columns<MT,SO,DF,SF,I,Is...>, 1UL >
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs >
struct IsRestricted< Columns<MT,SO,DF,SF,CCAs...> >
   : public IsRestricted<MT>
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
template< typename MT, bool SO, bool SF, size_t... CCAs >
struct HasConstDataAccess< Columns<MT,SO,true,SF,CCAs...> >
   : public HasConstDataAccess<MT>
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
template< typename MT, bool SO, bool SF, size_t... CCAs >
struct HasMutableDataAccess< Columns<MT,SO,true,SF,CCAs...> >
   : public HasMutableDataAccess<MT>
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
template< typename MT, bool SO, bool SF, size_t... CCAs >
struct IsAligned< Columns<MT,SO,true,SF,CCAs...> >
   : public IsAligned<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs, size_t... CRAs >
struct RowsTrait< Columns<MT,SO,DF,SF,CCAs...>, CRAs... >
{
   using Type = RowsTrait_t< ResultType_t< Columns<MT,SO,DF,SF,CCAs...> >, CRAs... >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs1, size_t... CCAs2 >
struct ColumnTrait< Columns<MT,SO,DF,SF,CCAs1...>, CCAs2... >
{
   using Type = ColumnTrait_t< ResultType_t< Columns<MT,SO,DF,SF,CCAs1...> >, CCAs2... >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs1, size_t... CCAs2 >
struct ColumnsTrait< Columns<MT,SO,DF,SF,CCAs1...>, CCAs2... >
{
   using Type = ColumnsTrait_t< ResultType_t< Columns<MT,SO,DF,SF,CCAs1...> >, CCAs2... >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs, ptrdiff_t... CBAs >
struct BandTrait< Columns<MT,SO,DF,SF,CCAs...>, CBAs... >
{
   using Type = BandTrait_t< ResultType_t< Columns<MT,SO,DF,SF,CCAs...> >, CBAs... >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
