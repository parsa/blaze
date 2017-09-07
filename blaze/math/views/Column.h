//=================================================================================================
/*!
//  \file blaze/math/views/Column.h
//  \brief Header file for the implementation of the Column view
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

#ifndef _BLAZE_MATH_VIEWS_COLUMN_H_
#define _BLAZE_MATH_VIEWS_COLUMN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

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
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsOpposedView.h>
#include <blaze/math/typetraits/IsSubmatrix.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/views/column/BaseTemplate.h>
#include <blaze/math/views/column/Dense.h>
#include <blaze/math/views/column/Sparse.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given matrix.
// \ingroup column
//
// \param matrix The matrix containing the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column<3UL>( D );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column<4UL>( S );
   \endcode

// In case the column is not properly specified (i.e. if the specified index is greater than or
// equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I     // Column index
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Column<MT,I> column( Matrix<MT,SO>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return Column<MT,I>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given constant matrix.
// \ingroup column
//
// \param matrix The constant matrix containing the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   const blaze::DynamicMatrix<double,columnMajor> D( ... );
   const blaze::CompressedMatrix<double,columnMajor> S( ... );

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column<3UL>( D );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column<4UL>( S );
   \endcode

// In case the column is not properly specified (i.e. if the specified index is greater than or
// equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I     // Column index
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline const Column<const MT,I> column( const Matrix<MT,SO>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return Column<const MT,I>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given temporary matrix.
// \ingroup column
//
// \param matrix The temporary matrix containing the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given temporary
// matrix. In case the column is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I     // Column index
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Column<MT,I> column( Matrix<MT,SO>&& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return Column<MT,I>( ~matrix );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given matrix.
// \ingroup column
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column( D, 3UL );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column( S, 4UL );
   \endcode

// In case the column is not properly specified (i.e. if the specified index is greater than or
// equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Column<MT> column( Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return Column<MT>( ~matrix, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given constant matrix.
// \ingroup column
//
// \param matrix The constant matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   const blaze::DynamicMatrix<double,columnMajor> D( ... );
   const blaze::CompressedMatrix<double,columnMajor> S( ... );

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column( D, 3UL );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column( S, 4UL );
   \endcode

// In case the column is not properly specified (i.e. if the specified index is greater than or
// equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline const Column<const MT> column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return Column<const MT>( ~matrix, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given temporary matrix.
// \ingroup column
//
// \param matrix The temporary matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given temporary
// matrix. In case the column is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline Column<MT> column( Matrix<MT,SO>&& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return Column<MT>( ~matrix, index );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix addition.
// \ingroup column
//
// \param matrix The constant matrix/matrix addition.
// \return View on the specified column of the addition.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// addition.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatAddExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return column<I>( (~matrix).leftOperand() ) + column<I>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix addition.
// \ingroup column
//
// \param matrix The constant matrix/matrix addition.
// \param index The index of the column.
// \return View on the specified column of the addition.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// addition.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatAddExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) + column( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix subtraction.
// \ingroup column
//
// \param matrix The constant matrix/matrix subtraction.
// \return View on the specified column of the subtraction.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// subtraction.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatSubExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return column<I>( (~matrix).leftOperand() ) - column<I>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix subtraction.
// \ingroup column
//
// \param matrix The constant matrix/matrix subtraction.
// \param index The index of the column.
// \return View on the specified column of the subtraction.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// subtraction.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatSubExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) - column( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given Schur product.
// \ingroup column
//
// \param matrix The constant Schur product.
// \return View on the specified column of the Schur product.
//
// This function returns an expression representing the specified column of the given Schur
// product.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const SchurExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return column<I>( (~matrix).leftOperand() ) * column<I>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given Schur product.
// \ingroup column
//
// \param matrix The constant Schur product.
// \param index The index of the column.
// \return View on the specified column of the Schur product.
//
// This function returns an expression representing the specified column of the given Schur
// product.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const SchurExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) * column( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix multiplication.
// \ingroup column
//
// \param matrix The constant matrix/matrix multiplication.
// \return View on the specified column of the multiplication.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// multiplication.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * column<I>( (~matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix multiplication.
// \ingroup column
//
// \param matrix The constant matrix/matrix multiplication.
// \param index The index of the column.
// \return View on the specified column of the multiplication.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// multiplication.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatMultExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * column( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given outer product.
// \ingroup column
//
// \param matrix The constant outer product.
// \return View on the specified column of the outer product.
//
// This function returns an expression representing the specified column of the given outer
// product.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const VecTVecMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * (~matrix).rightOperand()[I];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given outer product.
// \ingroup column
//
// \param matrix The constant outer product.
// \param index The index of the column.
// \return View on the specified column of the outer product.
//
// This function returns an expression representing the specified column of the given outer
// product.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const VecTVecMultExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * (~matrix).rightOperand()[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar multiplication.
// \ingroup column
//
// \param matrix The constant matrix/scalar multiplication.
// \return View on the specified column of the multiplication.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// multiplication.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return column<I>( (~matrix).leftOperand() ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar multiplication.
// \ingroup column
//
// \param matrix The constant matrix/scalar multiplication.
// \param index The index of the column.
// \return View on the specified column of the multiplication.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// multiplication.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatScalarMultExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar division.
// \ingroup column
//
// \param matrix The constant matrix/scalar division.
// \return View on the specified column of the division.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// division.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatScalarDivExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return column<I>( (~matrix).leftOperand() ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar division.
// \ingroup column
//
// \param matrix The constant matrix/scalar division.
// \param index The index of the column.
// \return View on the specified column of the division.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// division.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatScalarDivExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given unary matrix map operation.
// \ingroup column
//
// \param matrix The constant unary matrix map operation.
// \return View on the specified column of the unary map operation.
//
// This function returns an expression representing the specified column of the given unary
// matrix map operation.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMapExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return map( column<I>( (~matrix).operand() ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given unary matrix map operation.
// \ingroup column
//
// \param matrix The constant unary matrix map operation.
// \param index The index of the column.
// \return View on the specified column of the unary map operation.
//
// This function returns an expression representing the specified column of the given unary
// matrix map operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMapExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return map( column( (~matrix).operand(), index ), (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given binary matrix map operation.
// \ingroup column
//
// \param matrix The constant binary matrix map operation.
// \return View on the specified column of the binary map operation.
//
// This function returns an expression representing the specified column of the given binary
// matrix map operation.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatMapExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return map( column<I>( (~matrix).leftOperand() ),
               column<I>( (~matrix).rightOperand() ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given binary matrix map operation.
// \ingroup column
//
// \param matrix The constant binary matrix map operation.
// \param index The index of the column.
// \return View on the specified column of the binary map operation.
//
// This function returns an expression representing the specified column of the given binary
// matrix map operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatMatMapExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return map( column( (~matrix).leftOperand() , index ),
               column( (~matrix).rightOperand(), index ),
               (~matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix evaluation operation.
// \ingroup column
//
// \param matrix The constant matrix evaluation operation.
// \return View on the specified column of the evaluation operation.
//
// This function returns an expression representing the specified column of the given matrix
// evaluation operation.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatEvalExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return eval( column<I>( (~matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix evaluation operation.
// \ingroup column
//
// \param matrix The constant matrix evaluation operation.
// \param index The index of the column.
// \return View on the specified column of the evaluation operation.
//
// This function returns an expression representing the specified column of the given matrix
// evaluation operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatEvalExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return eval( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix serialization operation.
// \ingroup column
//
// \param matrix The constant matrix serialization operation.
// \return View on the specified column of the serialization operation.
//
// This function returns an expression representing the specified column of the given matrix
// serialization operation.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatSerialExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return serial( column<I>( (~matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix serialization operation.
// \ingroup column
//
// \param matrix The constant matrix serialization operation.
// \param index The index of the column.
// \return View on the specified column of the serialization operation.
//
// This function returns an expression representing the specified column of the given matrix
// serialization operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatSerialExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return serial( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix declaration operation.
// \ingroup column
//
// \param matrix The constant matrix declaration operation.
// \return View on the specified column of the declaration operation.
//
// This function returns an expression representing the specified column of the given matrix
// declaration operation.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const DeclExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return column<I>( (~matrix).operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix declaration operation.
// \ingroup column
//
// \param matrix The constant matrix declaration operation.
// \param index The index of the column.
// \return View on the specified column of the declaration operation.
//
// This function returns an expression representing the specified column of the given matrix
// declaration operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const DeclExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).operand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix transpose operation.
// \ingroup column
//
// \param matrix The constant matrix transpose operation.
// \return View on the specified column of the transpose operation.
//
// This function returns an expression representing the specified column of the given matrix
// transpose operation.
*/
template< size_t I       // Column index
        , typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatTransExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return trans( row<I>( (~matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix transpose operation.
// \ingroup column
//
// \param matrix The constant matrix transpose operation.
// \param index The index of the column.
// \return View on the specified column of the transpose operation.
//
// This function returns an expression representing the specified column of the given matrix
// transpose operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) column( const MatTransExpr<MT>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return trans( row( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMN OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given column.
// \ingroup column
//
// \param column The column to be resetted.
// \return void
*/
template< typename MT      // Type of the matrix
        , bool SO          // Storage order
        , bool DF          // Density flag
        , bool SF          // Symmetry flag
        , size_t... CAs >  // Compile time column arguments
inline void reset( ColumnImpl<MT,SO,DF,SF,CAs...>& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary column.
// \ingroup column
//
// \param column The temporary column to be resetted.
// \return void
*/
template< typename MT      // Type of the matrix
        , bool SO          // Storage order
        , bool DF          // Density flag
        , bool SF          // Symmetry flag
        , size_t... CAs >  // Compile time column arguments
inline void reset( ColumnImpl<MT,SO,DF,SF,CAs...>&& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given column.
// \ingroup column
//
// \param column The column to be cleared.
// \return void
//
// Clearing a column is equivalent to resetting it via the reset() function.
*/
template< typename MT      // Type of the matrix
        , bool SO          // Storage order
        , bool DF          // Density flag
        , bool SF          // Symmetry flag
        , size_t... CAs >  // Compile time column arguments
inline void clear( ColumnImpl<MT,SO,DF,SF,CAs...>& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary column.
// \ingroup column
//
// \param column The temporary column to be cleared.
// \return void
//
// Clearing a column is equivalent to resetting it via the reset() function.
*/
template< typename MT      // Type of the matrix
        , bool SO          // Storage order
        , bool DF          // Density flag
        , bool SF          // Symmetry flag
        , size_t... CAs >  // Compile time column arguments
inline void clear( ColumnImpl<MT,SO,DF,SF,CAs...>&& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense column is in default state.
// \ingroup column
//
// \param column The dense column to be tested for its default state.
// \return \a true in case the given dense column is component-wise zero, \a false otherwise.
//
// This function checks whether the dense column is in default state. For instance, in case the
// column is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all column elements are 0 and \a false in case any column element is not 0.
// The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( column( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< bool RF          // Relaxation flag
        , typename MT      // Type of the dense matrix
        , bool SO          // Storage order
        , bool SF          // Symmetry flag
        , size_t... CAs >  // Compile time column arguments
inline bool isDefault( const ColumnImpl<MT,SO,true,SF,CAs...>& column )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<column.size(); ++i )
      if( !isDefault<RF>( column[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse column is in default state.
// \ingroup column
//
// \param column The sparse column to be tested for its default state.
// \return \a true in case the given column is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse column is in default state. For instance, in case the
// column is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all column elements are 0 and \a false in case any column element is not 0.
// The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( column( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< bool RF          // Relaxation flag
        , typename MT      // Type of the sparse matrix
        , bool SO          // Storage order
        , bool SF          // Symmetry flag
        , size_t... CAs >  // Compile time column arguments
inline bool isDefault( const ColumnImpl<MT,SO,false,SF,CAs...>& column )
{
   using blaze::isDefault;

   for( const auto& element : column )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given column are intact.
// \ingroup column
//
// \param column The column to be tested.
// \return \a true in case the given column's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the column are intact, i.e. if its state is
// valid. In case the invariants are intact, the function returns \a true, else it will return
// \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isIntact( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT      // Type of the matrix
        , bool SO          // Storage order
        , bool DF          // Density flag
        , bool SF          // Symmetry flag
        , size_t... CAs >  // Compile time column arguments
inline bool isIntact( const ColumnImpl<MT,SO,DF,SF,CAs...>& column ) noexcept
{
   return ( column.column() < column.operand().columns() &&
            isIntact( column.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two regular columns.
// \ingroup column
//
// \param a The first column to be tested for its state.
// \param b The second column to be tested for its state.
// \return \a true in case the two columns share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of two
// regular columns. In case both columns represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1      // Type of the matrix of the left-hand side column
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF1          // Symmetry flag of the left-hand side column
        , size_t... CAs1    // Compile time column arguments of the left-hand side column
        , typename MT2      // Type of the matrix of the right-hand side column
        , bool SF2          // Symmetry flag of the right-hand side column
        , size_t... CAs2 >  // Compile time column arguments of the right-hand side column
inline DisableIf_< Or< IsSubmatrix<MT1>, IsSubmatrix<MT2> >, bool >
   isSame_backend( const ColumnImpl<MT1,SO,DF,SF1,CAs1...>& a,
                   const ColumnImpl<MT2,SO,DF,SF2,CAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand() ) && ( a.column() == b.column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the left column being a column on a submatrix.
// \ingroup column
//
// \param a The first column to be tested for its state.
// \param b The second column to be tested for its state.
// \return \a true in case the two columns share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the left
// column being a column on a submatrix. In case both columns represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1      // Type of the submatrix of the left-hand side column
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF1          // Symmetry flag of the left-hand side column
        , size_t... CAs1    // Compile time column arguments of the left-hand side column
        , typename MT2      // Type of the matrix of the right-hand side column
        , bool SF2          // Symmetry flag of the right-hand side column
        , size_t... CAs2 >  // Compile time column arguments of the right-hand side column
inline EnableIf_< And< IsSubmatrix<MT1>, Not< IsSubmatrix<MT2> > >, bool >
   isSame_backend( const ColumnImpl<MT1,SO,DF,SF1,CAs1...>& a,
                   const ColumnImpl<MT2,SO,DF,SF2,CAs2...>& b ) noexcept
{
   return ( isSame( a.operand().operand(), b.operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.column() + a.operand().column() == b.column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the right column being a column on a submatrix.
// \ingroup column
//
// \param a The first column to be tested for its state.
// \param b The second column to be tested for its state.
// \return \a true in case the two columns share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the right
// column being a column on a submatrix. In case both columns represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1      // Type of the matrix of the left-hand side column
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF1          // Symmetry flag of the left-hand side column
        , size_t... CAs1    // Compile time column arguments of the left-hand side column
        , typename MT2      // Type of the submatrix of the right-hand side column
        , bool SF2          // Symmetry flag of the right-hand side column
        , size_t... CAs2 >  // Compile time column arguments of the right-hand side column
inline EnableIf_< And< Not< IsSubmatrix<MT1> >, IsSubmatrix<MT2> >, bool >
   isSame_backend( const ColumnImpl<MT1,SO,DF,SF1,CAs1...>& a,
                   const ColumnImpl<MT2,SO,DF,SF2,CAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.column() == b.column() + b.operand().column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two columns on submatrices.
// \ingroup column
//
// \param a The first column to be tested for its state.
// \param b The second column to be tested for its state.
// \return \a true in case the two columns share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of both columns
// being columns on submatrices. In case both columns represent the same observable state, the
// function returns \a true, otherwise it returns \a false.
*/
template< typename MT1      // Type of the submatrix of the left-hand side column
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF1          // Symmetry flag of the left-hand side column
        , size_t... CAs1    // Compile time column arguments of the left-hand side column
        , typename MT2      // Type of the submatrix of the right-hand side column
        , bool SF2          // Symmetry flag of the right-hand side column
        , size_t... CAs2 >  // Compile time column arguments of the right-hand side column
inline EnableIf_< And< IsSubmatrix<MT1>, IsSubmatrix<MT2> >, bool >
   isSame_backend( const ColumnImpl<MT1,SO,DF,SF1,CAs1...>& a,
                   const ColumnImpl<MT2,SO,DF,SF2,CAs2...>& b ) noexcept
{
   return ( isSame( a.operand().operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.column() + a.operand().column() == b.column() + b.operand().column() ) &&
            ( a.operand().row() == b.operand().row() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given columns represent the same observable state.
// \ingroup column
//
// \param a The first column to be tested for its state.
// \param b The second column to be tested for its state.
// \return \a true in case the two columns share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given columns refer to exactly the
// same range of the same matrix. In case both columns represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1      // Type of the matrix of the left-hand side column
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF1          // Symmetry flag of the left-hand side column
        , size_t... CAs1    // Compile time column arguments of the left-hand side column
        , typename MT2      // Type of the matrix of the right-hand side column
        , bool SF2          // Symmetry flag of the right-hand side column
        , size_t... CAs2 >  // Compile time column arguments of the right-hand side column
inline bool isSame( const ColumnImpl<MT1,SO,DF,SF1,CAs1...>& a,
                    const ColumnImpl<MT2,SO,DF,SF2,CAs2...>& b ) noexcept
{
   return isSame_backend( a, b );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the matrix
        , bool SO        // Storage order
        , bool DF        // Density flag
        , bool SF        // Symmetry flag
        , size_t... CAs  // Compile time column arguments
        , typename VT >  // Type of the right-hand side vector
inline bool tryAssign( const ColumnImpl<MT,SO,DF,SF,CAs...>& lhs,
                       const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAssign( lhs.operand(), ~rhs, index, lhs.column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the matrix
        , bool SO        // Storage order
        , bool DF        // Density flag
        , bool SF        // Symmetry flag
        , size_t... CAs  // Compile time column arguments
        , typename VT >  // Type of the right-hand side vector
inline bool tryAddAssign( const ColumnImpl<MT,SO,SF,SF,CAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryAddAssign( lhs.operand(), ~rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the matrix
        , bool SO        // Storage order
        , bool DF        // Density flag
        , bool SF        // Symmetry flag
        , size_t... CAs  // Compile time column arguments
        , typename VT >  // Type of the right-hand side vector
inline bool trySubAssign( const ColumnImpl<MT,SO,DF,SF,CAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return trySubAssign( lhs.operand(), ~rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the matrix
        , bool SO        // Storage order
        , bool DF        // Density flag
        , bool SF        // Symmetry flag
        , size_t... CAs  // Compile time column arguments
        , typename VT >  // Type of the right-hand side vector
inline bool tryMultAssign( const ColumnImpl<MT,SO,DF,SF,CAs...>& lhs,
                           const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryMultAssign( lhs.operand(), ~rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the matrix
        , bool SO        // Storage order
        , bool DF        // Density flag
        , bool SF        // Symmetry flag
        , size_t... CAs  // Compile time column arguments
        , typename VT >  // Type of the right-hand side vector
inline bool tryDivAssign( const ColumnImpl<MT,SO,DF,SF,CAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.size() - index, "Invalid vector size" );

   return tryDivAssign( lhs.operand(), ~rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given column.
// \ingroup column
//
// \param c The column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given column. It returns a
// column object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Compile time column arguments
inline decltype(auto) derestrict( ColumnImpl<MT,SO,DF,SF,I>& c )
{
   return column<I>( derestrict( c.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary column.
// \ingroup column
//
// \param c The temporary column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary column. It
// returns a column object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Compile time column arguments
inline decltype(auto) derestrict( ColumnImpl<MT,SO,DF,SF,I>&& c )
{
   return column<I>( derestrict( c.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given column.
// \ingroup column
//
// \param c The column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given column. It returns a
// column object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( ColumnImpl<MT,SO,DF,SF>& c )
{
   return column( derestrict( c.operand() ), c.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary column.
// \ingroup column
//
// \param c The temporary column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary column. It
// returns a column object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( ColumnImpl<MT,SO,DF,SF>&& c )
{
   return column( derestrict( c.operand() ), c.column() );
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CAs >
struct IsRestricted< ColumnImpl<MT,SO,DF,SF,CAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CAs >
struct HasConstDataAccess< ColumnImpl<MT,SO,true,SF,CAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CAs >
struct HasMutableDataAccess< ColumnImpl<MT,SO,true,SF,CAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CAs >
struct IsAligned< ColumnImpl<MT,SO,true,SF,CAs...> >
   : public BoolConstant< And< IsAligned<MT>, Or< IsColumnMajorMatrix<MT>, IsSymmetric<MT> > >::value >
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
template< typename MT, bool SO, bool SF, size_t... CAs >
struct IsPadded< ColumnImpl<MT,SO,true,SF,CAs...> >
   : public BoolConstant< And< IsPadded<MT>, Or< IsColumnMajorMatrix<MT>, IsSymmetric<MT> > >::value >
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
template< typename MT, bool DF, size_t... CAs >
struct IsOpposedView< ColumnImpl<MT,false,DF,false,CAs...> >
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CAs, size_t... SAs >
struct SubvectorTrait< ColumnImpl<MT,SO,DF,SF,CAs...>, SAs... >
{
   using Type = SubvectorTrait_< ResultType_< ColumnImpl<MT,SO,DF,SF,CAs...> >, SAs... >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
