//=================================================================================================
/*!
//  \file blaze/math/views/Row.h
//  \brief Header file for the implementation of the Row view
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

#ifndef _BLAZE_MATH_VIEWS_ROW_H_
#define _BLAZE_MATH_VIEWS_ROW_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

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
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsOpposedView.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSubmatrix.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/row/BaseTemplate.h>
#include <blaze/math/views/row/Dense.h>
#include <blaze/math/views/row/Sparse.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/TrueType.h>
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
/*!\brief Creating a view on a specific row of the given matrix.
// \ingroup row
//
// \param matrix The matrix containing the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row<3UL>( D );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row<4UL>( S );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row<3UL>( D, unchecked );
   auto row4 = row<4UL>( S, unchecked );
   \endcode
*/
template< size_t I            // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT,I>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given constant matrix.
// \ingroup row
//
// \param matrix The constant matrix containing the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row<3UL>( D );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row<4UL>( S );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row<3UL>( D, unchecked );
   auto row4 = row<4UL>( S, unchecked );
   \endcode
*/
template< size_t I            // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const Matrix<MT,SO>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Row_<const MT,I>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given temporary matrix.
// \ingroup row
//
// \param matrix The temporary matrix containing the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given temporary
// matrix. In case the row is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the rows in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I            // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>&& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT,I>;
   return ReturnType( ~matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given matrix.
// \ingroup row
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row( D, 3UL );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row( S, 4UL );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row( D, 3UL, unchecked );
   auto row4 = row( S, 4UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT>;
   return ReturnType( ~matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given constant matrix.
// \ingroup row
//
// \param matrix The constant matrix containing the row.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row( D, 3UL );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row( S, 4UL );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row( D, 3UL, unchecked );
   auto row4 = row( S, 4UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const Matrix<MT,SO>& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Row_<const MT>;
   return ReturnType( ~matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given temporary matrix.
// \ingroup row
//
// \param matrix The temporary matrix containing the row.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given temporary
// matrix. In case the row is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the rows in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>&& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT>;
   return ReturnType( ~matrix, index, args... );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/matrix addition.
// \ingroup row
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime row arguments.
// \return View on the specified row of the addition.
//
// This function returns an expression representing the specified row of the given matrix/matrix
// addition.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatAddExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (~matrix).leftOperand(), args... ) +
          row<CRAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/matrix subtraction.
// \ingroup row
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime row arguments.
// \return View on the specified row of the subtraction.
//
// This function returns an expression representing the specified row of the given matrix/matrix
// subtraction.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatSubExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (~matrix).leftOperand(), args... ) -
          row<CRAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given Schur product.
// \ingroup row
//
// \param matrix The constant Schur product.
// \param args The runtime row arguments.
// \return View on the specified row of the Schur product.
//
// This function returns an expression representing the specified row of the given Schur product.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const SchurExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (~matrix).leftOperand(), args... ) *
          row<CRAs...>( (~matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/matrix multiplication.
// \ingroup row
//
// \param matrix The constant matrix/matrix multiplication.
// \param args The runtime row arguments
// \return View on the specified row of the multiplication.
//
// This function returns an expression representing the specified row of the given matrix/matrix
// multiplication.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (~matrix).leftOperand(), args... ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given outer product.
// \ingroup row
//
// \param matrix The constant outer product.
// \param args Optional row arguments.
// \return View on the specified row of the outer product.
//
// This function returns an expression representing the specified row of the given outer product.
*/
template< size_t I            // Row index
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const VecTVecMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   UNUSED_PARAMETER( args... );

   if( !Contains< TypeList<RRAs...>, Unchecked >::value ) {
      if( (~matrix).rows() <= I ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return (~matrix).leftOperand()[I] * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given outer product.
// \ingroup row
//
// \param matrix The constant outer product.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the outer product.
//
// This function returns an expression representing the specified row of the given outer product.
*/
template< typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const VecTVecMultExpr<MT>& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   UNUSED_PARAMETER( args... );

   if( !Contains< TypeList<RRAs...>, Unchecked >::value ) {
      if( (~matrix).rows() <= index ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return (~matrix).leftOperand()[index] * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/scalar multiplication.
// \ingroup row
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime row arguments
// \return View on the specified row of the multiplication.
//
// This function returns an expression representing the specified row of the given matrix/scalar
// multiplication.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatScalarMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (~matrix).leftOperand(), args... ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/scalar division.
// \ingroup row
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime row arguments
// \return View on the specified row of the division.
//
// This function returns an expression representing the specified row of the given matrix/scalar
// division.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatScalarDivExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (~matrix).leftOperand(), args... ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given unary matrix map operation.
// \ingroup row
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime row arguments
// \return View on the specified row of the unary map operation.
//
// This function returns an expression representing the specified row of the given unary matrix
// map operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMapExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( row<CRAs...>( (~matrix).operand(), args... ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given binary matrix map operation.
// \ingroup row
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime row arguments
// \return View on the specified row of the binary map operation.
//
// This function returns an expression representing the specified row of the given binary matrix
// map operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatMapExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( row<CRAs...>( (~matrix).leftOperand(), args... ),
               row<CRAs...>( (~matrix).rightOperand(), args... ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix evaluation operation.
// \ingroup row
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime row arguments
// \return View on the specified row of the evaluation operation.
//
// This function returns an expression representing the specified row of the given matrix
// evaluation operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatEvalExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( row<CRAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix serialization operation.
// \ingroup row
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime row arguments
// \return View on the specified row of the serialization operation.
//
// This function returns an expression representing the specified row of the given matrix
// serialization operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatSerialExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( row<CRAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix declaration operation.
// \ingroup row
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime row arguments
// \return View on the specified row of the declaration operation.
//
// This function returns an expression representing the specified row of the given matrix
// declaration operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const DeclExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (~matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix transpose operation.
// \ingroup row
//
// \param matrix The constant matrix transpose operation.
// \param args The runtime row arguments
// \return View on the specified row of the transpose operation.
//
// This function returns an expression representing the specified row of the given matrix
// transpose operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatTransExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return trans( column<CRAs...>( (~matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROW OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given row.
// \ingroup row
//
// \param row The row to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void reset( Row<MT,SO,DF,SF,CRAs...>& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary row.
// \ingroup row
//
// \param row The temporary row to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void reset( Row<MT,SO,DF,SF,CRAs...>&& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given row.
// \ingroup row
//
// \param row The row to be cleared.
// \return void
//
// Clearing a row is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void clear( Row<MT,SO,DF,SF,CRAs...>& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary row.
// \ingroup row
//
// \param row The temporary row to be cleared.
// \return void
//
// Clearing a row is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void clear( Row<MT,SO,DF,SF,CRAs...>&& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense row is in default state.
// \ingroup row
//
// \param row The dense row to be tested for its default state.
// \return \a true in case the given dense row is component-wise zero, \a false otherwise.
//
// This function checks whether the dense row is in default state. For instance, in case the
// row is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all row elements are 0 and \a false in case any row element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( row( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isDefault( const Row<MT,SO,true,SF,CRAs...>& row )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<row.size(); ++i )
      if( !isDefault<RF>( row[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse row is in default state.
// \ingroup row
//
// \param row The sparse row to be tested for its default state.
// \return \a true in case the given row is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse row is in default state. For instance, in case the
// row is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all row elements are 0 and \a false in case any row element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( row( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< bool RF           // Relaxation flag
        , typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isDefault( const Row<MT,SO,false,SF,CRAs...>& row )
{
   using blaze::isDefault;

   for( const auto& element : row )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given row are intact.
// \ingroup row
//
// \param row The row to be tested.
// \return \a true in case the given row's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the row are intact, i.e. if its state is valid.
// In case the invariants are intact, the function returns \a true, else it will return \a false.
// The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isIntact( const Row<MT,SO,DF,SF,CRAs...>& row ) noexcept
{
   return ( row.row() < row.operand().rows() &&
            isIntact( row.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two regular rows.
// \ingroup row
//
// \param a The first row to be tested for its state.
// \param b The second row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of two regular
// rows. In case both rows represent the same observable state, the function returns \a true,
// otherwise it returns \a false.
*/
template< typename MT1       // Type of the matrix of the left-hand side row
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF1           // Symmetry flag of the left-hand side row
        , size_t... CRAs1    // Compile time row arguments of the left-hand side row
        , typename MT2       // Type of the matrix of the right-hand side row
        , bool SF2           // Symmetry flag of the right-hand side row
        , size_t... CRAs2 >  // Compile time row arguments of the right-hand side row
inline DisableIf_< Or< IsSubmatrix<MT1>, IsSubmatrix<MT2> >, bool >
   isSame_backend( const Row<MT1,SO,DF,SF1,CRAs1...>& a,
                   const Row<MT2,SO,DF,SF2,CRAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand() ) && ( a.row() == b.row() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the left row being a row on a submatrix.
// \ingroup row
//
// \param a The first row to be tested for its state.
// \param b The second row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the left row
// being a row on a submatrix. In case both rows represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1       // Type of the submatrix of the left-hand side row
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF1           // Symmetry flag of the left-hand side row
        , size_t... CRAs1    // Compile time row arguments of the left-hand side row
        , typename MT2       // Type of the matrix of the right-hand side row
        , bool SF2           // Symmetry flag of the right-hand side row
        , size_t... CRAs2 >  // Compile time row arguments of the right-hand side row
inline EnableIf_< And< IsSubmatrix<MT1>, Not< IsSubmatrix<MT2> > >, bool >
   isSame_backend( const Row<MT1,SO,DF,SF1,CRAs1...>& a,
                   const Row<MT2,SO,DF,SF2,CRAs2...>& b ) noexcept
{
   return ( isSame( a.operand().operand(), b.operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() + a.operand().row() == b.row() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the right row being a row on a submatrix.
// \ingroup row
//
// \param a The first row to be tested for its state.
// \param b The second row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the right row
// being a row on a submatrix. In case both rows represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1       // Type of the matrix of the left-hand side row
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF1           // Symmetry flag of the left-hand side row
        , size_t... CRAs1    // Compile time row arguments of the left-hand side row
        , typename MT2       // Type of the submatrix of the right-hand side row
        , bool SF2           // Symmetry flag of the right-hand side row
        , size_t... CRAs2 >  // Compile time row arguments of the right-hand side row
inline EnableIf_< And< Not< IsSubmatrix<MT1> >, IsSubmatrix<MT2> >, bool >
   isSame_backend( const Row<MT1,SO,DF,SF1,CRAs1...>& a,
                   const Row<MT2,SO,DF,SF2,CRAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() == b.row() + b.operand().row() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two rows on submatrices.
// \ingroup row
//
// \param a The first row to be tested for its state.
// \param b The second row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of both rows
// being rows on submatrices. In case both rows represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1       // Type of the submatrix of the left-hand side row
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF1           // Symmetry flag of the left-hand side row
        , size_t... CRAs1    // Compile time row arguments of the left-hand side row
        , typename MT2       // Type of the submatrix of the right-hand side row
        , bool SF2           // Symmetry flag of the right-hand side row
        , size_t... CRAs2 >  // Compile time row arguments of the right-hand side row
inline EnableIf_< And< IsSubmatrix<MT1>, IsSubmatrix<MT2> >, bool >
   isSame_backend( const Row<MT1,SO,DF,SF1,CRAs1...>& a,
                   const Row<MT2,SO,DF,SF2,CRAs2...>& b ) noexcept
{
   return ( isSame( a.operand().operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() + a.operand().row() == b.row() + b.operand().row() ) &&
            ( a.operand().column() == b.operand().column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given rows represent the same observable state.
// \ingroup row
//
// \param a The first row to be tested for its state.
// \param b The second row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This overload of the isSame() function tests if the two given rows refer to exactly the same
// range of the same matrix. In case both rows represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1       // Type of the matrix of the left-hand side row
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF1           // Symmetry flag of the left-hand side row
        , size_t... CRAs1    // Compile time row arguments of the left-hand side row
        , typename MT2       // Type of the matrix of the right-hand side row
        , bool SF2           // Symmetry flag of the right-hand side row
        , size_t... CRAs2 >  // Compile time row arguments of the right-hand side row
inline bool isSame( const Row<MT1,SO,DF,SF1,CRAs1...>& a,
                    const Row<MT2,SO,DF,SF2,CRAs2...>& b ) noexcept
{
   return isSame_backend( a, b );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
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
inline bool tryAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                       const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.operand(), ~rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
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
inline bool tryAddAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.operand(), ~rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
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
inline bool trySubAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.operand(), ~rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
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
inline bool tryMultAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                           const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.operand(), ~rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
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
inline bool tryDivAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryDivAssign( lhs.operand(), ~rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given row.
// \ingroup row
//
// \param r The row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given row. It returns a row
// object that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Row index
inline decltype(auto) derestrict( Row<MT,SO,DF,SF,I>& r )
{
   return row<I>( derestrict( r.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary row.
// \ingroup row
//
// \param r The temporary row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary row. It
// returns a row object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Row index
inline decltype(auto) derestrict( Row<MT,SO,DF,SF,I>&& r )
{
   return row<I>( derestrict( r.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given row.
// \ingroup row
//
// \param r The row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given row. It returns a row
// object that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Row<MT,SO,DF,SF>& r )
{
   return row( derestrict( r.operand() ), r.row(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary row.
// \ingroup row
//
// \param r The temporary row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary row. It
// returns a row object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Row<MT,SO,DF,SF>&& r )
{
   return row( derestrict( r.operand() ), r.row(), unchecked );
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct IsRestricted< Row<MT,SO,DF,SF,CRAs...> >
   : public BoolConstant< IsRestricted<MT>::value >
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
struct HasConstDataAccess< Row<MT,SO,true,SF,CRAs...> >
   : public BoolConstant< HasConstDataAccess<MT>::value >
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
struct HasMutableDataAccess< Row<MT,SO,true,SF,CRAs...> >
   : public BoolConstant< HasMutableDataAccess<MT>::value >
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
struct IsAligned< Row<MT,SO,true,SF,CRAs...> >
   : public BoolConstant< And< IsAligned<MT>, Or< IsRowMajorMatrix<MT>, IsSymmetric<MT> > >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct IsPadded< Row<MT,SO,true,SF,CRAs...> >
   : public BoolConstant< And< IsPadded<MT>, Or< IsRowMajorMatrix<MT>, IsSymmetric<MT> > >::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISOPPOSEDVIEW SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool DF, size_t... CRAs >
struct IsOpposedView< Row<MT,false,DF,false,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTORTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs, size_t... CSAs >
struct SubvectorTrait< Row<MT,SO,DF,SF,CRAs...>, CSAs... >
{
   using Type = SubvectorTrait_< ResultType_< Row<MT,SO,DF,SF,CRAs...> >, CSAs... >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
