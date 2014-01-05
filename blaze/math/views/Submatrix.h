//=================================================================================================
/*!
//  \file blaze/math/views/Submatrix.h
//  \brief Header file for all restructuring submatrix functions
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

#ifndef _BLAZE_MATH_VIEWS_SUBMATRIX_H_
#define _BLAZE_MATH_VIEWS_SUBMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/traits/SubmatrixExprTrait.h>
#include <blaze/math/typetraits/IsMatAbsExpr.h>
#include <blaze/math/typetraits/IsMatEvalExpr.h>
#include <blaze/math/typetraits/IsMatMatAddExpr.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsMatMatSubExpr.h>
#include <blaze/math/typetraits/IsMatScalarDivExpr.h>
#include <blaze/math/typetraits/IsMatScalarMultExpr.h>
#include <blaze/math/typetraits/IsMatTransExpr.h>
#include <blaze/math/typetraits/IsVecTVecMultExpr.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL RESTRUCTURING OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix addition.
// \ingroup views
//
// \param matrix The constant matrix/matrix addition.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the addition.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix addition.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatAddExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix( (~matrix).leftOperand() , row, column, m, n ) +
          submatrix( (~matrix).rightOperand(), row, column, m, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix subtraction.
// \ingroup views
//
// \param matrix The constant matrix/matrix subtraction.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the subtraction.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix subtraction.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatSubExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix( (~matrix).leftOperand() , row, column, m, n ) -
          submatrix( (~matrix).rightOperand(), row, column, m, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/matrix multiplication.
// \ingroup views
//
// \param matrix The constant matrix/matrix multiplication.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the multiplication.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/matrix multiplication.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatMultExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   typename MT::LeftOperand  left ( (~matrix).leftOperand()  );
   typename MT::RightOperand right( (~matrix).rightOperand() );

   return submatrix( left, row, 0UL, m, left.columns() ) * submatrix( right, 0UL, column, right.rows(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given outer product.
// \ingroup views
//
// \param matrix The constant outer product.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the outer product.
//
// This function returns an expression representing the specified submatrix of the given
// outer product.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsVecTVecMultExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return subvector( (~matrix).leftOperand(), row, m ) * subvector( (~matrix).rightOperand(), column, n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/scalar multiplication.
// \ingroup views
//
// \param matrix The constant matrix/scalar multiplication.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the multiplication.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/scalar multiplication.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarMultExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix( (~matrix).leftOperand(), row, column, m, n ) * (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix/scalar division.
// \ingroup views
//
// \param matrix The constant matrix/scalar division.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the division.
//
// This function returns an expression representing the specified submatrix of the given
// matrix/scalar division.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarDivExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix( (~matrix).leftOperand(), row, column, m, n ) / (~matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix abs operation.
// \ingroup views
//
// \param matrix The constant matrix abs operation.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the abs operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// abs operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatAbsExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return abs( submatrix( (~matrix).operand(), row, column, m, n ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix evaluation operation.
// \ingroup views
//
// \param matrix The constant matrix evaluation operation.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the evaluation operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// evaluation operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatEvalExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return eval( submatrix( (~matrix).operand(), row, column, m, n ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given matrix transpose operation.
// \ingroup views
//
// \param matrix The constant matrix transpose operation.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specified submatrix of the transpose operation.
//
// This function returns an expression representing the specified submatrix of the given matrix
// transpose operation.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatTransExpr<MT>, typename SubmatrixExprTrait<MT,false>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return trans( submatrix( (~matrix).operand(), column, row, n, m ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
