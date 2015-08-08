//=================================================================================================
/*!
//  \file blaze/math/views/Column.h
//  \brief Header file for all restructuring column functions
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
//=================================================================================================

#ifndef _BLAZE_MATH_VIEWS_COLUMN_H_
#define _BLAZE_MATH_VIEWS_COLUMN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/traits/ColumnExprTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsMatAbsExpr.h>
#include <blaze/math/typetraits/IsMatConjExpr.h>
#include <blaze/math/typetraits/IsMatEvalExpr.h>
#include <blaze/math/typetraits/IsMatImagExpr.h>
#include <blaze/math/typetraits/IsMatMatAddExpr.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsMatMatSubExpr.h>
#include <blaze/math/typetraits/IsMatRealExpr.h>
#include <blaze/math/typetraits/IsMatScalarDivExpr.h>
#include <blaze/math/typetraits/IsMatScalarMultExpr.h>
#include <blaze/math/typetraits/IsMatSerialExpr.h>
#include <blaze/math/typetraits/IsMatTransExpr.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/typetraits/IsVecTVecMultExpr.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given matrix.
// \ingroup views
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix.

   \code
   using blaze::columnMajor;

   typedef blaze::DynamicMatrix<double,columnMajor>     DenseMatrix;
   typedef blaze::CompressedMatrix<double,columnMajor>  SparseMatrix;

   DenseMatrix D;
   SparseMatrix S;
   // ... Resizing and initialization

   // Creating a view on the 3rd column of the dense matrix D
   blaze::DenseColumn<DenseMatrix> = column( D, 3UL );

   // Creating a view on the 4th column of the sparse matrix S
   blaze::SparseColumn<SparseMatrix> = column( S, 4UL );
   \endcode
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename DisableIf< Or< IsComputation<MT>, IsTransExpr<MT> >
                         , typename ColumnExprTrait<MT>::Type >::Type
   column( Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   typedef typename ColumnExprTrait<MT>::Type  ReturnType;
   return ReturnType( ~matrix, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given constant matrix.
// \ingroup views
//
// \param matrix The constant matrix containing the column.
// \param index The index of the column.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix.

   \code
   using blaze::columnMajor;

   typedef blaze::DynamicMatrix<double,columnMajor>     DenseMatrix;
   typedef blaze::CompressedMatrix<double,columnMajor>  SparseMatrix;

   DenseMatrix D;
   SparseMatrix S;
   // ... Resizing and initialization

   // Creating a view on the 3rd column of the dense matrix D
   blaze::DenseColumn<DenseMatrix> = column( D, 3UL );

   // Creating a view on the 4th column of the sparse matrix S
   blaze::SparseColumn<SparseMatrix> = column( S, 4UL );
   \endcode
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename DisableIf< Or< IsComputation<MT>, IsTransExpr<MT> >
                         , typename ColumnExprTrait<const MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   typedef typename ColumnExprTrait<const MT>::Type  ReturnType;
   return ReturnType( ~matrix, index );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix addition.
// \ingroup views
//
// \param matrix The constant matrix/matrix addition.
// \param index The index of the column.
// \return View on the specified column of the addition.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// addition.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatAddExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) + column( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix subtraction.
// \ingroup views
//
// \param matrix The constant matrix/matrix subtraction.
// \param index The index of the column.
// \return View on the specified column of the subtraction.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// subtraction.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatSubExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) - column( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix multiplication.
// \ingroup views
//
// \param matrix The constant matrix/matrix multiplication.
// \param index The index of the column.
// \return View on the specified column of the multiplication.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// multiplication.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatMultExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * column( (~matrix).rightOperand(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given outer product.
// \ingroup views
//
// \param matrix The constant outer product.
// \param index The index of the column.
// \return View on the specified column of the outer product.
//
// This function returns an expression representing the specified column of the given outer
// product.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsVecTVecMultExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return (~matrix).leftOperand() * (~matrix).rightOperand()[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar multiplication.
// \ingroup views
//
// \param matrix The constant matrix/scalar multiplication.
// \param index The index of the column.
// \return View on the specified column of the multiplication.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// multiplication.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarMultExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar division.
// \ingroup views
//
// \param matrix The constant matrix/scalar division.
// \param index The index of the column.
// \return View on the specified column of the division.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// division.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarDivExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return column( (~matrix).leftOperand(), index ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix \a abs operation.
// \ingroup views
//
// \param matrix The constant matrix \a abs operation.
// \param index The index of the column.
// \return View on the specified column of the \a abs operation.
//
// This function returns an expression representing the specified column of the given matrix
// \a abs operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatAbsExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return abs( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix \a conj operation.
// \ingroup views
//
// \param matrix The constant matrix \a conj operation.
// \param index The index of the column.
// \return View on the specified column of the \a conj operation.
//
// This function returns an expression representing the specified column of the given matrix
// \a conj operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatConjExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return conj( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix \a real operation.
// \ingroup views
//
// \param matrix The constant matrix \a real operation.
// \param index The index of the column.
// \return View on the specified column of the \a real operation.
//
// This function returns an expression representing the specified column of the given matrix
// \a real operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatRealExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return real( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix \a imag operation.
// \ingroup views
//
// \param matrix The constant matrix \a imag operation.
// \param index The index of the column.
// \return View on the specified column of the \a imag operation.
//
// This function returns an expression representing the specified column of the given matrix
// \a imag operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatImagExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return imag( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix evaluation operation.
// \ingroup views
//
// \param matrix The constant matrix evaluation operation.
// \param index The index of the column.
// \return View on the specified column of the evaluation operation.
//
// This function returns an expression representing the specified column of the given matrix
// evaluation operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatEvalExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return eval( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix serialization operation.
// \ingroup views
//
// \param matrix The constant matrix serialization operation.
// \param index The index of the column.
// \return View on the specified column of the serialization operation.
//
// This function returns an expression representing the specified column of the given matrix
// serialization operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatSerialExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return serial( column( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix transpose operation.
// \ingroup views
//
// \param matrix The constant matrix transpose operation.
// \param index The index of the column.
// \return View on the specified column of the transpose operation.
//
// This function returns an expression representing the specified column of the given matrix
// transpose operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatTransExpr<MT>, typename ColumnExprTrait<MT>::Type >::Type
   column( const Matrix<MT,SO>& matrix, size_t index )
{
   BLAZE_FUNCTION_TRACE;

   return trans( row( (~matrix).operand(), index ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
