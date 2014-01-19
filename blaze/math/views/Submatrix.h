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
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsMatAbsExpr.h>
#include <blaze/math/typetraits/IsMatEvalExpr.h>
#include <blaze/math/typetraits/IsMatMatAddExpr.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/math/typetraits/IsMatMatSubExpr.h>
#include <blaze/math/typetraits/IsMatScalarDivExpr.h>
#include <blaze/math/typetraits/IsMatScalarMultExpr.h>
#include <blaze/math/typetraits/IsMatTransExpr.h>
#include <blaze/math/typetraits/IsTransExpr.h>
#include <blaze/math/typetraits/IsVecTVecMultExpr.h>
#include <blaze/math/views/AlignmentFlag.h>
#include <blaze/util/DisableIf.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup views
//
// \param matrix The matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given matrix.
// The following example demonstrates the creation of a dense and sparse submatrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>     DenseMatrix;
   typedef blaze::CompressedMatrix<int,columnMajor>  SparseMatrix;

   DenseMatrix  D;
   SparseMatrix S;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 16
   blaze::DenseSubmatrix<DenseMatrix> dsm = submatrix( D, 0UL, 16UL, 8UL, 4UL );

   // Creating a sparse submatrix of size 7x3, starting in row 2 and column 4
   blaze::SparseSubmatrix<SparseMatrix> ssm = submatrix( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// In case the submatrix is not properly specified (i.e. if the specified row or column is larger
// than the total number of rows or columns of the given matrix or the submatrix is specified
// beyond the number of rows or columns of the matrix) a \a std::invalid_argument exception is
// thrown.
//
// Please note that this function creates an unaligned dense or sparse submatrix. For instance,
// the creation of the dense submatrix is equivalent to the following three function calls:

   \code
   blaze::DenseSubmatrix<DenseMatrix>           dsm = submatrix<unaligned>( D, 0UL, 16UL, 8UL, 4UL );
   blaze::DenseSubmatrix<DenseMatrix,unaligned> dsm = submatrix           ( D, 0UL, 16UL, 8UL, 4UL );
   blaze::DenseSubmatrix<DenseMatrix,unaligned> dsm = submatrix<unaligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions. However, especially in case of dense submatrices this may
// result in considerable performance improvements. In order to create an aligned submatrix the
// following function call has to be used:

   \code
   blaze::DenseSubmatrix<DenseMatrix,aligned> dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// Note however that in this case the given arguments \a row, \a columns, \a m, and \a n are
// subject to additional checks to guarantee proper alignment.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename SubmatrixExprTrait<MT,unaligned>::Type
   submatrix( Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned>( ~matrix, row, column, m, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup views
//
// \param matrix The matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given matrix.
// The following example demonstrates the creation of a dense and sparse submatrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>     DenseMatrix;
   typedef blaze::CompressedMatrix<int,columnMajor>  SparseMatrix;

   DenseMatrix  D;
   SparseMatrix S;
   // ... Resizing and initialization

   // Creating a dense submatrix of size 8x4, starting in row 0 and column 16
   blaze::DenseSubmatrix<DenseMatrix> dsm = submatrix( D, 0UL, 16UL, 8UL, 4UL );

   // Creating a sparse submatrix of size 7x3, starting in row 2 and column 4
   blaze::SparseSubmatrix<SparseMatrix> ssm = submatrix( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// In case the submatrix is not properly specified (i.e. if the specified row or column is larger
// than the total number of rows or columns of the given matrix or the submatrix is specified
// beyond the number of rows or columns of the matrix) a \a std::invalid_argument exception is
// thrown.
//
// Please note that this function creates an unaligned dense or sparse submatrix. For instance,
// the creation of the dense submatrix is equivalent to the following three function calls:

   \code
   blaze::DenseSubmatrix<DenseMatrix>           dsm = submatrix<unaligned>( D, 0UL, 16UL, 8UL, 4UL );
   blaze::DenseSubmatrix<DenseMatrix,unaligned> dsm = submatrix           ( D, 0UL, 16UL, 8UL, 4UL );
   blaze::DenseSubmatrix<DenseMatrix,unaligned> dsm = submatrix<unaligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions. However, especially in case of dense submatrices this may
// result in considerable performance improvements. In order to create an aligned submatrix the
// following function call has to be used:

   \code
   blaze::DenseSubmatrix<DenseMatrix,aligned> dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );
   \endcode

// Note however that in this case the given arguments \a row, \a columns, \a m, and \a n are
// subject to additional checks to guarantee proper alignment.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename SubmatrixExprTrait<const MT,unaligned>::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<unaligned>( ~matrix, row, column, m, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup views
//
// \param matrix The matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given dense or sparse matrix, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned submatrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>     DenseMatrix;
   typedef blaze::CompressedMatrix<int,columnMajor>  SparseMatrix;

   DenseMatrix  D;
   SparseMatrix S;
   // ... Resizing and initialization

   // Creating an aligned dense submatrix of size 8x4, starting in row 0 and column 16
   blaze::DenseSubmatrix<DenseMatrix,aligned> dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );

   // Creating an unaligned sparse submatrix of size 7x3, starting in row 2 and column 4
   blaze::SparseSubmatrix<SparseMatrix,unaligned> ssm = submatrix<unaligned>( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// In case the submatrix is not properly specified (i.e. if the specified row or column is larger
// than the total number of rows or columns of the given matrix or the submatrix is specified
// beyond the number of rows or columns of the matrix) a \a std::invalid_argument exception is
// thrown.
//
// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions and the given \a row, \a column, \a m, and \a n arguments are
// subject to additional checks to guarantee proper alignment. However, especially in case of dense
// submatrices this may result in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). The following source code gives some
// examples for a double precision dense matrix, assuming that AVX is available, which packs 4
// \c double values into an intrinsic vector:

   \code
   using blaze::rowMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>      MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType,aligned>  SubmatrixType;

   MatrixType D( 13UL, 17UL );
   // ... Resizing and initialization

   // OK: Starts at position (0,0) and the number of rows and columns are a multiple of 4
   SubmatrixType dsm1 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 12UL );

   // OK: First row and column and the number of rows and columns are all a multiple of 4
   SubmatrixType dsm2 = submatrix<aligned>( D, 4UL, 12UL, 8UL, 16UL );

   // OK: First row and column are a multiple of 4 and the submatrix includes the last row and column
   SubmatrixType dsm3 = submatrix<aligned>( D, 4UL, 0UL, 9UL, 17UL );

   // Error: First row is not a multiple of 4
   SubmatrixType dsm4 = submatrix<aligned>( D, 2UL, 4UL, 12UL, 12UL );

   // Error: First column is not a multiple of 4
   SubmatrixType dsm5 = submatrix<aligned>( D, 0UL, 2UL, 8UL, 8UL );

   // Error: The number of rows is not a multiple of 4 and the submatrix does not include the last row
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 7UL, 8UL );

   // Error: The number of columns is not a multiple of 4 and the submatrix does not include the last column
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 11UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DisableIf< Or< IsComputation<MT>, IsTransExpr<MT> >
                         , typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   typedef typename SubmatrixExprTrait<MT,AF>::Type  ReturnType;
   return ReturnType( ~matrix, row, column, m, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific submatrix of the given matrix.
// \ingroup views
//
// \param matrix The matrix containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \return View on the specific submatrix of the matrix.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing an aligned or unaligned submatrix of the
// given dense or sparse matrix, based on the specified alignment flag \a AF. The following
// example demonstrates the creation of both an aligned and unaligned submatrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>     DenseMatrix;
   typedef blaze::CompressedMatrix<int,columnMajor>  SparseMatrix;

   DenseMatrix  D;
   SparseMatrix S;
   // ... Resizing and initialization

   // Creating an aligned dense submatrix of size 8x4, starting in row 0 and column 16
   blaze::DenseSubmatrix<DenseMatrix,aligned> dsm = submatrix<aligned>( D, 0UL, 16UL, 8UL, 4UL );

   // Creating an unaligned sparse submatrix of size 7x3, starting in row 2 and column 4
   blaze::SparseSubmatrix<SparseMatrix,unaligned> ssm = submatrix<unaligned>( S, 2UL, 4UL, 7UL, 3UL );
   \endcode

// In case the submatrix is not properly specified (i.e. if the specified row or column is larger
// than the total number of rows or columns of the given matrix or the submatrix is specified
// beyond the number of rows or columns of the matrix) a \a std::invalid_argument exception is
// thrown.
//
// In contrast to unaligned submatrices, which provide full flexibility, aligned submatrices pose
// additional alignment restrictions and the given \a row, \a column, \a m, and \a n arguments are
// subject to additional checks to guarantee proper alignment. However, especially in case of dense
// submatrices this may result in considerable performance improvements.
//
// The alignment restrictions refer to system dependent address restrictions for the used element
// type and the available vectorization mode (SSE, AVX, ...). The following source code gives some
// examples for a double precision dense matrix, assuming that AVX is available, which packs 4
// \c double values into an intrinsic vector:

   \code
   using blaze::rowMajor;

   typedef blaze::DynamicMatrix<double,rowMajor>      MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType,aligned>  SubmatrixType;

   MatrixType D( 13UL, 17UL );
   // ... Resizing and initialization

   // OK: Starts at position (0,0) and the number of rows and columns are a multiple of 4
   SubmatrixType dsm1 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 12UL );

   // OK: First row and column and the number of rows and columns are all a multiple of 4
   SubmatrixType dsm2 = submatrix<aligned>( D, 4UL, 12UL, 8UL, 16UL );

   // OK: First row and column are a multiple of 4 and the submatrix includes the last row and column
   SubmatrixType dsm3 = submatrix<aligned>( D, 4UL, 0UL, 9UL, 17UL );

   // Error: First row is not a multiple of 4
   SubmatrixType dsm4 = submatrix<aligned>( D, 2UL, 4UL, 12UL, 12UL );

   // Error: First column is not a multiple of 4
   SubmatrixType dsm5 = submatrix<aligned>( D, 0UL, 2UL, 8UL, 8UL );

   // Error: The number of rows is not a multiple of 4 and the submatrix does not include the last row
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 7UL, 8UL );

   // Error: The number of columns is not a multiple of 4 and the submatrix does not include the last column
   SubmatrixType dsm6 = submatrix<aligned>( D, 0UL, 0UL, 8UL, 11UL );
   \endcode

// In case any alignment restrictions are violated, a \a std::invalid_argument exception is thrown.
*/
template< bool AF      // Alignment flag
        , typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline typename DisableIf< Or< IsComputation<MT>, IsTransExpr<MT> >
                         , typename SubmatrixExprTrait<const MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   typedef typename SubmatrixExprTrait<const MT,AF>::Type  ReturnType;
   return ReturnType( ~matrix, row, column, m, n );
}
//*************************************************************************************************




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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatAddExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF>( (~matrix).leftOperand() , row, column, m, n ) +
          submatrix<AF>( (~matrix).rightOperand(), row, column, m, n );
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatSubExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF>( (~matrix).leftOperand() , row, column, m, n ) -
          submatrix<AF>( (~matrix).rightOperand(), row, column, m, n );
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatMatMultExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   typename MT::LeftOperand  left ( (~matrix).leftOperand()  );
   typename MT::RightOperand right( (~matrix).rightOperand() );

   return submatrix<AF>( left, row, 0UL, m, left.columns() ) *
          submatrix<AF>( right, 0UL, column, right.rows(), n );
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsVecTVecMultExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return subvector<AF>( (~matrix).leftOperand(), row, m ) *
          subvector<AF>( (~matrix).rightOperand(), column, n );
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarMultExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF>( (~matrix).leftOperand(), row, column, m, n ) * (~matrix).rightOperand();
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatScalarDivExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return submatrix<AF>( (~matrix).leftOperand(), row, column, m, n ) / (~matrix).rightOperand();
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatAbsExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return abs( submatrix<AF>( (~matrix).operand(), row, column, m, n ) );
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatEvalExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return eval( submatrix<AF>( (~matrix).operand(), row, column, m, n ) );
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
template< bool AF      // Alignment flag
        , typename MT  // Type of the matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsMatTransExpr<MT>, typename SubmatrixExprTrait<MT,AF>::Type >::Type
   submatrix( const Matrix<MT,SO>& matrix, size_t row, size_t column, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   return trans( submatrix<AF>( (~matrix).operand(), column, row, n, m ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
