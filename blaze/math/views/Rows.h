//=================================================================================================
/*!
//  \file blaze/math/views/Rows.h
//  \brief Header file for the implementation of the Rows view
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

#ifndef _BLAZE_MATH_VIEWS_ROWS_H_
#define _BLAZE_MATH_VIEWS_ROWS_H_


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
#include <blaze/math/expressions/MatVecMultExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/IntegerSequence.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/Forward.h>
#include <blaze/math/views/row/RowData.h>
#include <blaze/math/views/rows/BaseTemplate.h>
#include <blaze/math/views/rows/Dense.h>
#include <blaze/math/views/rows/Sparse.h>
#include <blaze/util/Assert.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
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
/*!\brief Creating a view on a selection of rows of the given matrix.
// \ingroup rows
//
// \param matrix The matrix containing the rows.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd row of the dense matrix D
   auto rows1 = rows<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd row of the sparse matrix S
   auto rows2 = rows<4UL,2UL>( S );
   \endcode

// By default, the provided row indices are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than or equal to the total number of rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto rows1 = rows<1UL,3UL>( D, unchecked );
   auto rows2 = rows<4UL,2UL>( S, unchecked );
   \endcode
*/
template< size_t I            // First row index
        , size_t... Is        // Remaining row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( Matrix<MT,SO>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Rows_<MT,I,Is...>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of rows of the given constant matrix.
// \ingroup rows
//
// \param matrix The constant matrix containing the rows.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the 1st and 3rd row of the dense matrix D
   auto rows1 = rows<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd row of the sparse matrix S
   auto rows2 = rows<4UL,2UL>( S );
   \endcode

// By default, the provided row indices are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than or equal to the total number of rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto rows1 = rows<1UL,3UL>( D, unchecked );
   auto rows2 = rows<4UL,2UL>( S, unchecked );
   \endcode
*/
template< size_t I            // First row index
        , size_t... Is        // Remaining row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( const Matrix<MT,SO>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Rows_<const MT,I,Is...>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of rows of the given temporary matrix.
// \ingroup rows
//
// \param matrix The temporary matrix containing the rows.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given temporary
// matrix. In case any row is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of rows in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I            // First row index
        , size_t... Is        // Remaining row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( Matrix<MT,SO>&& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Rows_<MT,I,Is...>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of rows of the given matrix.
// \ingroup rows
//
// \param matrix The matrix containing the rows.
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd row of the dense matrix D
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto rows1 = rows( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd row of the sparse matrix S
   const std::array<size_t,2uL> indices2{ 4UL, 2UL };
   auto rows2 = rows( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided row indices are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than or equal to the total number of rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto rows1 = rows( D, indices1.data(), indices1.size(), unchecked );
   auto rows2 = rows( S, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( Matrix<MT,SO>& matrix, const T* indices, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Rows_<MT>;
   return ReturnType( ~matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of rows of the given constant matrix.
// \ingroup rows
//
// \param matrix The constant matrix containing the rows.
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd row of the dense matrix D
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto rows1 = rows( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd row of the sparse matrix S
   const std::array<size_t,2uL> indices2{ 4UL, 2UL };
   auto rows2 = rows( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided row indices are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than or equal to the total number of rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto rows1 = rows( D, indices1.data(), indices1.size(), unchecked );
   auto rows2 = rows( S, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( const Matrix<MT,SO>& matrix, const T* indices, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Rows_<const MT>;
   return ReturnType( ~matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of rows of the given temporary matrix.
// \ingroup rows
//
// \param matrix The temporary matrix containing the rows.
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args Optional arguments.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given temporary
// matrix. In case any row is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of rows in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( Matrix<MT,SO>&& matrix, const T* indices, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Rows_<MT>;
   return ReturnType( ~matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows of the given matrix.
// \ingroup rows
//
// \param matrix The matrix containing the rows.
// \param indices The sequence of row indices.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given matrix. In
// case any row is not properly specified (i.e. if any specified index is greater than or equal
// to the total number of rows in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< typename MT         // Type of the matrix
        , size_t... Is        // Row indices
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( MT&& matrix, index_sequence<Is...> indices, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   UNUSED_PARAMETER( indices );

   return rows<Is...>( std::forward<MT>( matrix ), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows of the given matrix.
// \ingroup rows
//
// \param matrix The matrix containing the rows.
// \param indices The list of row indices.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given matrix. In
// case any row is not properly specified (i.e. if any specified index is greater than or equal
// to the total number of rows in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( MT&& matrix, initializer_list<T> indices, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows( std::forward<MT>( matrix ), indices.begin(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows of the given matrix.
// \ingroup rows
//
// \param matrix The matrix containing the rows.
// \param indices The array of row indices.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given matrix. In
// case any row is not properly specified (i.e. if any specified index is greater than or equal
// to the total number of rows in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the row indices
        , size_t N            // Number of indices
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( MT&& matrix, const std::array<T,N>& indices, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows( std::forward<MT>( matrix ), indices.data(), N, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows of the given matrix.
// \ingroup rows
//
// \param matrix The matrix containing the rows.
// \param indices The vector of row indices.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given matrix. In
// case any row is not properly specified (i.e. if any specified index is greater than or equal
// to the total number of rows in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( MT&& matrix, const std::vector<T>& indices, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows( std::forward<MT>( matrix ), indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows of the given matrix.
// \ingroup rows
//
// \param matrix The matrix containing the rows.
// \param indices The vector of row indices.
// \param args Optional arguments.
// \return View on the specified rows of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing a selection of rows of the given matrix. In
// case any row is not properly specified (i.e. if any specified index is greater than or equal
// to the total number of rows in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the row indices
        , size_t N            // Number of preallocated elements
        , typename... RRAs >  // Optional arguments
inline decltype(auto) rows( MT&& matrix, const SmallVector<T,N>& indices, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows( std::forward<MT>( matrix ), indices.data(), indices.size(), args... );
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
/*!\brief Creating a view on a selection of rows on the given matrix/matrix addition.
// \ingroup rows
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the addition.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix/matrix addition.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatMatAddExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CRAs...>( (~matrix).leftOperand(), args... ) +
          rows<CRAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix/matrix subtraction.
// \ingroup rows
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the subtraction.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix/matrix subtraction.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatMatSubExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CRAs...>( (~matrix).leftOperand(), args... ) -
          rows<CRAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given Schur product.
// \ingroup rows
//
// \param matrix The constant Schur product.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the subtraction.
//
// This function returns an expression representing the specified selection of rows on the given
// Schur product.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const SchurExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CRAs...>( (~matrix).leftOperand(), args... ) %
          rows<CRAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix/matrix multiplication.
// \ingroup rows
//
// \param matrix The constant matrix/matrix multiplication.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix/matrix multiplication.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatMatMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CRAs...>( (~matrix).leftOperand(), args... ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given outer product.
// \ingroup rows
//
// \param matrix The constant outer product.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the outer product.
//
// This function returns an expression representing the specified selection of rows on the given
// outer product.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const VecTVecMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<CRAs...>( (~matrix).leftOperand(), args... ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix/scalar multiplication.
// \ingroup rows
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix/scalar multiplication.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatScalarMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CRAs...>( (~matrix).leftOperand(), args... ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix/scalar division.
// \ingroup rows
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the division.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix/scalar division.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatScalarDivExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CRAs...>( (~matrix).leftOperand(), args... ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given unary matrix map operation.
// \ingroup rows
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// unary matrix map operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatMapExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( rows<CRAs...>( (~matrix).operand(), args... ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given binary matrix map operation.
// \ingroup rows
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// binary matrix map operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatMatMapExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( rows<CRAs...>( (~matrix).leftOperand(), args... ),
               rows<CRAs...>( (~matrix).rightOperand(), args... ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix evaluation operation.
// \ingroup rows
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix evaluation operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatEvalExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( rows<CRAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix serialization operation.
// \ingroup rows
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix serialization operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatSerialExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( rows<CRAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix declaration operation.
// \ingroup rows
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix declaration operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const DeclExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CRAs...>( (~matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of rows on the given matrix transpose operation.
// \ingroup rows
//
// \param matrix The constant matrix transpose operation.
// \param args The runtime row arguments.
// \return View on the specified selection of rows on the multiplication.
//
// This function returns an expression representing the specified selection of rows on the given
// matrix transpose operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , typename... RRAs  // Runtime row arguments
        , typename = EnableIf_t< ( sizeof...( CRAs ) + sizeof...( RRAs ) > 0UL ) > >
inline decltype(auto) rows( const MatTransExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return trans( columns<CRAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given row selection.
// \ingroup rows
//
// \param r The selection of rows containing the rows.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
//
// This function returns an expression representing the specified rows of the given row selection.
*/
template< size_t I1           // First required row index
        , size_t... Is1       // Remaining required row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First present row index
        , size_t... Is2       // Remaining present row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( Rows<MT,SO,DF,SF,I2,Is2...>& r, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return rows< indices[I1], indices[Is1]... >( r.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given constant row selection.
// \ingroup rows
//
// \param r The constant selection of rows containing the rows.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
//
// This function returns an expression representing the specified rows of the given constant row
// selection.
*/
template< size_t I1           // First required row index
        , size_t... Is1       // Remaining required row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First present row index
        , size_t... Is2       // Remaining present row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( const Rows<MT,SO,DF,SF,I2,Is2...>& r, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return rows< indices[I1], indices[Is1]... >( r.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given temporary row selection.
// \ingroup rows
//
// \param r The temporary selection of rows containing the rows.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
//
// This function returns an expression representing the specified rows of the given temporary row
// selection.
*/
template< size_t I1           // First required row index
        , size_t... Is1       // Remaining required row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First present row index
        , size_t... Is2       // Remaining present row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( Rows<MT,SO,DF,SF,I2,Is2...>&& r, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is2... };
   return rows< indices[I1], indices[Is1]... >( r.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given row selection.
// \ingroup rows
//
// \param r The selection of rows containing the rows.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified rows of the given row selection.
*/
template< size_t I            // First required row index
        , size_t... Is        // Remaining required row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs      // Compile time row arguments
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( Rows<MT,SO,DF,SF,CRAs...>& r, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( r.rows() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   decltype(auto) indices( r.idces() );
   return rows( r.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given constant row selection.
// \ingroup rows
//
// \param r The constant selection of rows containing the rows.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified rows of the given constant row
// selection.
*/
template< size_t I            // First required row index
        , size_t... Is        // Remaining required row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs      // Compile time row arguments
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( const Rows<MT,SO,DF,SF,CRAs...>& r, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( r.rows() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   decltype(auto) indices( r.idces() );
   return rows( r.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given temporary row selection.
// \ingroup rows
//
// \param r The temporary selection of rows containing the rows.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified rows of the given temporary row
// selection.
*/
template< size_t I            // First required row index
        , size_t... Is        // Remaining required row indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs      // Compile time row arguments
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( Rows<MT,SO,DF,SF,CRAs...>&& r, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( r.rows() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   decltype(auto) indices( r.idces() );
   return rows( r.operand(), { indices[I], indices[Is]... }, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given row selection.
// \ingroup rows
//
// \param r The selection of rows containing the rows.
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given row selection.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs      // Compile time row arguments
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( Rows<MT,SO,DF,SF,CRAs...>& r, const T* indices, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( r.rows() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   decltype(auto) oldIndices( r.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return rows( r.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given constant row selection.
// \ingroup rows
//
// \param r The constant selection of rows containing the rows.
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given constant row
// selection.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs      // Compile time row arguments
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( const Rows<MT,SO,DF,SF,CRAs...>& r, const T* indices, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( r.rows() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   decltype(auto) oldIndices( r.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return rows( r.operand(), newIndices.data(), newIndices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given temporary row selection.
// \ingroup rows
//
// \param r The temporary selection of rows containing the rows.
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args The optional row arguments.
// \return View on the specified rows of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given temporary row
// selection.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs      // Compile time row arguments
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( Rows<MT,SO,DF,SF,CRAs...>&& r, const T* indices, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( r.rows() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   decltype(auto) oldIndices( r.idces() );
   SmallVector<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( oldIndices[indices[i]] );
   }

   return rows( r.operand(), newIndices.data(), newIndices.size(), args... );
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
/*!\brief Creating a view on a selection of elements of the given matrix/vector multiplication.
// \ingroup rows
//
// \param vector The constant matrix/vector multiplication.
// \param args The runtime element arguments.
// \return View on the specified elements of the multiplication.
//
// This function returns an expression representing the specified elements of the given
// matrix/vector multiplication.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const MatVecMultExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return rows<CEAs...>( (~vector).leftOperand(), args... ) * (~vector).rightOperand();
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
/*!\brief Creating a view on a specific row of the given row selection.
// \ingroup rows
//
// \param rows The selection of rows containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the row selection.
//
// This function returns an expression representing the specified row of the given row selection.
*/
template< size_t I1           // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First row index
        , size_t... Is        // Remaining row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Rows<MT,SO,DF,SF,I2,Is...>& rows, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is... };
   return row<indices[I1]>( rows.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given constant row selection.
// \ingroup rows
//
// \param rows The constant selection of rows containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the row selection.
//
// This function returns an expression representing the specified row of the given constant row
// selection.
*/
template< size_t I1           // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First row index
        , size_t... Is        // Remaining row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const Rows<MT,SO,DF,SF,I2,Is...>& rows, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is... };
   return row<indices[I1]>( rows.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given temporary row selection.
// \ingroup rows
//
// \param rows The temporary selection of rows containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the row selection.
//
// This function returns an expression representing the specified row of the given temporary row
// selection.
*/
template< size_t I1           // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I2           // First row index
        , size_t... Is        // Remaining row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Rows<MT,SO,DF,SF,I2,Is...>&& rows, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   static constexpr size_t indices[] = { I2, Is... };
   return row<indices[I1]>( rows.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given row selection.
// \ingroup rows
//
// \param rows The selection of rows containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given row selection.
*/
template< size_t... CRAs1     // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs2     // Compile time row arguments
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Rows<MT,SO,DF,SF,CRAs2...>& rows, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const RowData<CRAs1...> rd( args... );
   decltype(auto) indices( rows.idces() );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      if( indices.size() <= rd.row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return row( rows.operand(), indices[rd.row()], args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given constant row selection.
// \ingroup rows
//
// \param rows The constant selection of rows containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given constant row
// selection.
*/
template< size_t... CRAs1     // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs2     // Compile time row arguments
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const Rows<MT,SO,DF,SF,CRAs2...>& rows, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const RowData<CRAs1...> rd( args... );
   decltype(auto) indices( rows.idces() );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      if( indices.size() <= rd.row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return row( rows.operand(), indices[rd.row()], args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given temporary row selection.
// \ingroup rows
//
// \param rows The temporary selection of rows containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the row selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given temporary row
// selection.
*/
template< size_t... CRAs1     // Compile time row arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t... CRAs2     // Compile time row arguments
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( Rows<MT,SO,DF,SF,CRAs2...>&& rows, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const RowData<CRAs1...> rd( args... );
   decltype(auto) indices( rows.idces() );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      if( indices.size() <= rd.row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return row( rows.operand(), indices[rd.row()], args... );
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
/*!\brief Creating a view on a specific column of the given row selection.
// \ingroup rows
//
// \param rows The selection of rows containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the row selection.
//
// This function returns an expression representing the specified column of the given row selection.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I            // First row index
        , size_t... Is        // Remaining row indices
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Rows<MT,SO,DF,SF,I,Is...>& rows, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<I,Is...>( column<CCAs...>( rows.operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given constant row selection.
// \ingroup rows
//
// \param rows The constant selection of rows containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the row selection.
//
// This function returns an expression representing the specified column of the given constant
// row selection.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I            // First row index
        , size_t... Is        // Remaining row indices
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const Rows<MT,SO,DF,SF,I,Is...>& rows, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<I,Is...>( column<CCAs...>( rows.operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given temporary row selection.
// \ingroup rows
//
// \param rows The temporary selection of rows containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the row selection.
//
// This function returns an expression representing the specified column of the given temporary
// row selection.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , size_t I            // First row index
        , size_t... Is        // Remaining row indices
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Rows<MT,SO,DF,SF,I,Is...>&& rows, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements<I,Is...>( column<CCAs...>( rows.operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given row selection.
// \ingroup rows
//
// \param rows The selection of rows containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the row selection.
//
// This function returns an expression representing the specified column of the given row selection.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Rows<MT,SO,DF,SF>& rows, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( column<CCAs...>( rows.operand(), args... ), rows.idces() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given constant row selection.
// \ingroup rows
//
// \param rows The constant selection of rows containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the row selection.
//
// This function returns an expression representing the specified column of the given constant
// row selection.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const Rows<MT,SO,DF,SF>& rows, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( column<CCAs...>( rows.operand(), args... ), rows.idces() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given temporary row selection.
// \ingroup rows
//
// \param rows The temporary selection of rows containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the row selection.
//
// This function returns an expression representing the specified column of the given temporary
// row selection.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( Rows<MT,SO,DF,SF>&& rows, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( column<CCAs...>( rows.operand(), args... ), rows.idces() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given row selection.
// \ingroup rows
//
// \param rows The row selection to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void reset( Rows<MT,SO,DF,SF,CRAs...>& rows )
{
   rows.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary row selection.
// \ingroup rows
//
// \param rows The temporary row selection to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void reset( Rows<MT,SO,DF,SF,CRAs...>&& rows )
{
   rows.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row of the given row selection.
// \ingroup rows
//
// \param rows The row selection to be resetted.
// \param i The index of the row to be resetted.
// \return void
//
// This function resets the values in the specified row of the given row selection to their
// default value. Note that the capacity of the row remains unchanged.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void reset( Rows<MT,SO,DF,SF,CRAs...>& rows, size_t i )
{
   rows.reset( i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given row selection.
// \ingroup rows
//
// \param rows The row selection to be cleared.
// \return void
//
// Clearing a row selection is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void clear( Rows<MT,SO,DF,SF,CRAs...>& rows )
{
   rows.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary row selection.
// \ingroup rows
//
// \param rows The row selection to be cleared.
// \return void
//
// Clearing a row selection is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void clear( Rows<MT,SO,DF,SF,CRAs...>&& rows )
{
   rows.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense row selection is in default state.
// \ingroup rows
//
// \param rows The dense row selection to be tested for its default state.
// \return \a true in case the given dense row selection is component-wise zero, \a false otherwise.
//
// This function checks whether the dense row selection is in default state. For instance, in
// case the row selection is instantiated for a built-in integral or floating point data type,
// the function returns \a true in case all row elements are 0 and \a false in case any row
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( rows( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( rows( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isDefault( const Rows<MT,SO,true,SF,CRAs...>& rows )
{
   using blaze::isDefault;

   if( SO == true ) {
      for( size_t i=0UL; i<rows.rows(); ++i )
         for( size_t j=0UL; j<rows.columns(); ++j )
            if( !isDefault<RF>( rows(i,j) ) )
               return false;
   }
   else {
      for( size_t j=0UL; j<rows.columns(); ++j )
         for( size_t i=0UL; i<rows.rows(); ++i )
            if( !isDefault<RF>( rows(i,j) ) )
               return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse row selection is in default state.
// \ingroup rows
//
// \param rows The sparse row selection to be tested for its default state.
// \return \a true in case the given sparse row selection is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse row selection is in default state. For instance, in
// case the row selection is instantiated for a built-in integral or floating point data type,
// the function returns \a true in case all row elements are 0 and \a false in case any row
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( rows( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( rows( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isDefault( const Rows<MT,SO,false,SF,CRAs...>& rows )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<rows.rows(); ++i ) {
      for( auto element=rows.cbegin(i); element!=rows.cend(i); ++element )
         if( !isDefault<RF>( element->value() ) ) return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given row selection are intact.
// \ingroup rows
//
// \param rows The row selection to be tested.
// \return \a true in case the given row selection's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the row selection are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::DynamicMatrix<double,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( rows( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isIntact( const Rows<MT,SO,DF,SF,CRAs...>& rows ) noexcept
{
   return ( rows.rows() <= rows.operand().rows() &&
            rows.columns() == rows.operand().columns() &&
            isIntact( rows.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given row selection and matrix represent the same observable state.
// \ingroup rows
//
// \param a The row selection to be tested for its state.
// \param b The matrix to be tested for its state.
// \return \a true in case the row selection and matrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given row selection refers to all rows of
// the given matrix in ascending and consecutive order and by that represents the same observable
// state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT     // Type of the matrix
        , bool SO1        // Storage order of the left-hand side row selection
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CRAs  // Compile time row arguments
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool isSame( const Rows<MT,SO1,DF,SF,CRAs...>& a, const Matrix<MT,SO2>& b ) noexcept
{
   if( !isSame( a.operand(), ~b ) || ( a.rows() != (~b).rows() ) || ( a.columns() != (~b).columns() ) )
      return false;

   decltype(auto) indices( a.idces() );
   for( size_t i=0UL; i<a.rows(); ++i ) {
      if( indices[i] != i )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given matrix and row selection represent the same observable state.
// \ingroup rows
//
// \param a The matrix to be tested for its state.
// \param b The row selection to be tested for its state.
// \return \a true in case the matrix and row selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given row selection refers to all rows of
// the given matrix in ascending and consecutive order and by that represents the same observable
// state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO1          // Storage order of the left-hand side matrix
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs    // Compile time row arguments
        , bool SO2 >        // Storage order of the right-hand side row selection
inline bool isSame( const Matrix<MT,SO1>& a, const Rows<MT,SO2,DF,SF,CRAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given row selection and submatrix represent the same observable state.
// \ingroup rows
//
// \param a The row selection to be tested for its state.
// \param b The submatrix to be tested for its state.
// \return \a true in case the row selection and submatrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given row selection refers to same rows
// as the given submatrix in ascending and consecutive order and by that represents the same
// observable state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO1          // Storage order of the left-hand side row selection
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs    // Compile time row arguments
        , AlignmentFlag AF  // Alignment flag
        , bool SO2          // Storage order of the right-hand side submatrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isSame( const Rows<MT,SO1,DF,SF,CRAs...>& a, const Submatrix<MT,AF,SO2,DF,CSAs...>& b ) noexcept
{
   if( !isSame( a.operand(), b.operand() ) || ( a.rows() != (~b).rows() ) || ( a.columns() != (~b).columns() ) )
      return false;

   decltype(auto) indices( a.idces() );
   for( size_t i=0UL; i<a.rows(); ++i ) {
      if( indices[i] != b.row()+i )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given submatrix and row selection represent the same observable state.
// \ingroup rows
//
// \param a The submatrix to be tested for its state.
// \param b The row selection to be tested for its state.
// \return \a true in case the submatrix and row selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given row selection refers to same rows
// as the given submatrix in ascending and consecutive order and by that represents the same
// observable state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , bool SO1          // Storage order of the left-hand side submatrix
        , bool DF           // Density flag
        , size_t... CSAs    // Compile time submatrix arguments
        , bool SO2          // Storage order of the right-hand side row selection
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isSame( const Submatrix<MT,AF,SO1,DF,CSAs...>& a, const Rows<MT,SO2,DF,SF,CRAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given row selections represent the same observable state.
// \ingroup rows
//
// \param a The first selection of rows to be tested for its state.
// \param b The second selection of rows to be tested for its state.
// \return \a true in case the two row selections share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given row selections refer to exactly
// the same range of the same matrix. In case both selections represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT        // Type of the matrix
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF            // Symmetry flag
        , size_t... CRAs1    // Compile time row arguments of the left-hand side row selection
        , size_t... CRAs2 >  // Compile time row arguments of the right-hand side row selection
inline bool isSame( const Rows<MT,SO,DF,SF,CRAs1...>& a,
                    const Rows<MT,SO,DF,SF,CRAs2...>& b ) noexcept
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
/*!\brief In-place inversion of the given dense row selection.
// \ingroup rows
//
// \param r The dense row selection to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given dense row selection by means of the specified matrix type or
// matrix inversion algorithm \c IF (see the InversionFlag documentation):

   \code
   invert<asLower>( A );     // Inversion of a lower triangular matrix
   invert<asUniUpper>( A );  // Inversion of an upper unitriangular matrix
   invert<byLU>( A );        // Inversion by means of an LU decomposition
   invert<byLLH>( A );       // Inversion by means of a Cholesky decomposition
   ...
   \endcode

// The matrix inversion fails if ...
//
//  - ... the given row selection is not a square matrix;
//  - ... the given row selection is singular and not invertible.
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
// exception \a r may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline DisableIf_t< HasMutableDataAccess_v<MT> > invert( Rows<MT,SO,true,SF,CRAs...>& r )
{
   using RT = ResultType_t< Rows<MT,SO,true,SF,CRAs...> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION  ( RT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( RT );

   RT tmp( r );
   invert<IF>( tmp );
   r = tmp;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a row selection.
// \ingroup rows
//
// \param r The target row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool trySet( const Rows<MT,SO,DF,SF,CRAs...>& r, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < r.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < r.columns(), "Invalid column access index" );

   return trySet( r.operand(), r.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a row selection.
// \ingroup rows
//
// \param r The target row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryAdd( const Rows<MT,SO,DF,SF,CRAs...>& r, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < r.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < r.columns(), "Invalid column access index" );

   return tryAdd( r.operand(), r.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a row selection.
// \ingroup rows
//
// \param r The target row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool trySub( const Rows<MT,SO,DF,SF,CRAs...>& r, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < r.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < r.columns(), "Invalid column access index" );

   return trySub( r.operand(), r.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a row selection.
// \ingroup rows
//
// \param r The target row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryMult( const Rows<MT,SO,DF,SF,CRAs...>& r, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < r.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < r.columns(), "Invalid column access index" );

   return tryMult( r.operand(), r.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a row selection.
// \ingroup rows
//
// \param r The target row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Rows<MT,SO,DF,SF,CRAs...>& r, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (~r).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (~r).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (~r).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (~r).columns(), "Invalid number of columns" );

   const size_t iend( row + m );

   for( size_t i=row; i<iend; ++i ) {
      if( !tryMult( r.operand(), r.idx(i), column, m, n, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a row selection.
// \ingroup rows
//
// \param r The target row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryDiv( const Rows<MT,SO,DF,SF,CRAs...>& r, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < r.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < r.columns(), "Invalid column access index" );

   return tryDiv( r.operand(), r.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a row selection.
// \ingroup rows
//
// \param r The target row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Rows<MT,SO,DF,SF,CRAs...>& r, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (~r).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (~r).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (~r).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (~r).columns(), "Invalid number of columns" );

   const size_t iend( row + m );

   for( size_t i=row; i<iend; ++i ) {
      if( !tryDiv( r.operand(), r.idx(i), column, m, n, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a column vector to a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                       const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySet( lhs.operand(), lhs.idx( row+i ), column, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a row vector to a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                       const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   return tryAssign( lhs.operand(), ~rhs, lhs.idx( row ), column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to the band of a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                       const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySet( lhs.operand(), lhs.idx( row+i ), column+i, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a matrix to a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool tryAssign( const Rows<MT1,SO1,DF,SF,CRAs...>& lhs,
                       const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).rows(); ++i ) {
      if( !tryAssign( lhs.operand(), blaze::row( ~rhs, i, unchecked ), lhs.idx( row+i ), column ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a column vector to a row
//        selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAddAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryAdd( lhs.operand(), lhs.idx( row+i ), column, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a row vector to a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAddAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   return tryAddAssign( lhs.operand(), ~rhs, lhs.idx( row ), column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to the band of a
//        row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryAdd( lhs.operand(), lhs.idx( row+i ), column+i, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool tryAddAssign( const Rows<MT1,SO1,DF,SF,CRAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).rows(); ++i ) {
      if( !tryAddAssign( lhs.operand(), blaze::row( ~rhs, i, unchecked ), lhs.idx( row+i ), column ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a column vector to a row
//        selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool trySubAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySub( lhs.operand(), lhs.idx( row+i ), column, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a row vector to a row
//        selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool trySubAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   return trySubAssign( lhs.operand(), ~rhs, lhs.idx( row ), column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to the band of
//        a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool trySubAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !trySub( lhs.operand(), lhs.idx( row+i ), column+i, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool trySubAssign( const Rows<MT1,SO1,DF,SF,CRAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).rows(); ++i ) {
      if( !trySubAssign( lhs.operand(), blaze::row( ~rhs, i, unchecked ), lhs.idx( row+i ), column ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a column vector to a
//        row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryMultAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                           const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryMult( lhs.operand(), lhs.idx( row+i ), column, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a row vector to a row
//        selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryMultAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                           const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   return tryMultAssign( lhs.operand(), ~rhs, lhs.idx( row ), column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to the band
//        of a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryMultAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                           const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryMult( lhs.operand(), lhs.idx( row+i ), column+i, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the Schur product assignment of a matrix to a row
//        selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename MT2    // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline bool trySchurAssign( const Rows<MT1,SO1,DF,SF,CRAs...>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).rows(); ++i ) {
      if( !tryMultAssign( lhs.operand(), blaze::row( ~rhs, i, unchecked ), lhs.idx( row+i ), column ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a column vector to a row
//        selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryDivAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryDiv( lhs.operand(), lhs.idx( row+i ), column, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a row vector to a row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryDivAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   return tryDivAssign( lhs.operand(), ~rhs, lhs.idx( row ), column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to the band of a
//        row selection.
// \ingroup rows
//
// \param lhs The target left-hand side row selection.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT     // Type of the right-hand side vector
        , bool TF >       // Transpose flag of the right-hand side vector
inline bool tryDivAssign( const Rows<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   UNUSED_PARAMETER( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(~rhs).size(); ++i ) {
      if( !tryDiv( lhs.operand(), lhs.idx( row+i ), column+i, (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given row selection.
// \ingroup rows
//
// \param r The row selection to be derestricted.
// \return Row selection without access restrictions.
//
// This function removes all restrictions on the data access to the given row selection. It
// returns a row selection that does provide the same interface but does not have any restrictions
// on the data access.\n
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
inline decltype(auto) derestrict( Rows<MT,SO,DF,SF,I,Is...>& r )
{
   return rows<I,Is...>( derestrict( r.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary row selection.
// \ingroup rows
//
// \param r The temporary row selection to be derestricted.
// \return Row selection without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary row selection.
// It returns a row selection that does provide the same interface but does not have any
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
inline decltype(auto) derestrict( Rows<MT,SO,DF,SF,I,Is...>&& r )
{
   return rows<I,Is...>( derestrict( r.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given row selection.
// \ingroup rows
//
// \param r The row selection to be derestricted.
// \return Row selection without access restrictions.
//
// This function removes all restrictions on the data access to the given row selection. It
// returns a row selection that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Rows<MT,SO,DF,SF>& r )
{
   decltype(auto) indices( r.idces() );
   return rows( derestrict( r.operand() ), indices.data(), indices.size(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary row selection.
// \ingroup rows
//
// \param r The temporary row selection to be derestricted.
// \return Row selection without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary row selection.
// It returns a row selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Rows<MT,SO,DF,SF>&& r )
{
   decltype(auto) indices( r.idces() );
   return rows( derestrict( r.operand() ), indices.data(), indices.size(), unchecked );
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
template< typename MT, bool SO, bool DF, bool SF, size_t I, size_t... Is >
struct Size< Rows<MT,SO,DF,SF,I,Is...>, 0UL >
   : public PtrdiffT<1UL+sizeof...(Is)>
{};

template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct Size< Rows<MT,SO,DF,SF,CRAs...>, 1UL >
   : public Size<MT,1UL>
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
template< typename MT, bool SO, bool DF, bool SF, size_t I, size_t... Is >
struct MaxSize< Rows<MT,SO,DF,SF,I,Is...>, 0UL >
   : public PtrdiffT<1UL+sizeof...(Is)>
{};

template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct MaxSize< Rows<MT,SO,DF,SF,CRAs...>, 1UL >
   : public MaxSize<MT,1UL>
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct IsRestricted< Rows<MT,SO,DF,SF,CRAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct HasConstDataAccess< Rows<MT,SO,true,SF,CRAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct HasMutableDataAccess< Rows<MT,SO,true,SF,CRAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct IsAligned< Rows<MT,SO,true,SF,CRAs...> >
   : public IsAligned<MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
