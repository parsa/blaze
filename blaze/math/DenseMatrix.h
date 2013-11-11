//=================================================================================================
/*!
//  \file blaze/math/DenseMatrix.h
//  \brief Header file for the DenseMatrix CRTP base class
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

#ifndef _BLAZE_MATH_DENSEMATRIX_H_
#define _BLAZE_MATH_DENSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DMatAbsExpr.h>
#include <blaze/math/expressions/DMatDMatAddExpr.h>
#include <blaze/math/expressions/DMatDMatMultExpr.h>
#include <blaze/math/expressions/DMatDMatSubExpr.h>
#include <blaze/math/expressions/DMatDVecMultExpr.h>
#include <blaze/math/expressions/DMatEvalExpr.h>
#include <blaze/math/expressions/DMatScalarDivExpr.h>
#include <blaze/math/expressions/DMatScalarMultExpr.h>
#include <blaze/math/expressions/DMatSMatAddExpr.h>
#include <blaze/math/expressions/DMatSMatMultExpr.h>
#include <blaze/math/expressions/DMatSMatSubExpr.h>
#include <blaze/math/expressions/DMatSVecMultExpr.h>
#include <blaze/math/expressions/DMatTDMatAddExpr.h>
#include <blaze/math/expressions/DMatTDMatMultExpr.h>
#include <blaze/math/expressions/DMatTDMatSubExpr.h>
#include <blaze/math/expressions/DMatTransExpr.h>
#include <blaze/math/expressions/DMatTransposer.h>
#include <blaze/math/expressions/DMatTSMatAddExpr.h>
#include <blaze/math/expressions/DMatTSMatMultExpr.h>
#include <blaze/math/expressions/DMatTSMatSubExpr.h>
#include <blaze/math/expressions/DVecTDVecMultExpr.h>
#include <blaze/math/expressions/SMatDMatMultExpr.h>
#include <blaze/math/expressions/SMatDMatSubExpr.h>
#include <blaze/math/expressions/SMatTDMatMultExpr.h>
#include <blaze/math/expressions/SMatTDMatSubExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/TDMatDMatMultExpr.h>
#include <blaze/math/expressions/TDMatDVecMultExpr.h>
#include <blaze/math/expressions/TDMatSMatAddExpr.h>
#include <blaze/math/expressions/TDMatSMatMultExpr.h>
#include <blaze/math/expressions/TDMatSMatSubExpr.h>
#include <blaze/math/expressions/TDMatSVecMultExpr.h>
#include <blaze/math/expressions/TDMatTDMatMultExpr.h>
#include <blaze/math/expressions/TDMatTSMatMultExpr.h>
#include <blaze/math/expressions/TDVecDMatMultExpr.h>
#include <blaze/math/expressions/TDVecTDMatMultExpr.h>
#include <blaze/math/expressions/TSMatDMatMultExpr.h>
#include <blaze/math/expressions/TSMatDMatSubExpr.h>
#include <blaze/math/expressions/TSMatTDMatMultExpr.h>
#include <blaze/math/expressions/TSVecDMatMultExpr.h>
#include <blaze/math/expressions/TSVecTDMatMultExpr.h>
#include <blaze/math/Matrix.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseMatrix operators */
//@{
template< typename T1, typename T2 >
inline bool operator==( const DenseMatrix<T1,false>& lhs, const DenseMatrix<T2,false>& rhs );

template< typename T1, typename T2 >
inline bool operator==( const DenseMatrix<T1,true>& lhs, const DenseMatrix<T2,true>& rhs );

template< typename T1, typename T2, bool SO >
inline bool operator==( const DenseMatrix<T1,SO>& lhs, const DenseMatrix<T2,!SO>& rhs );

template< typename T1, typename T2, bool SO >
inline bool operator==( const DenseMatrix<T1,SO>& lhs, const SparseMatrix<T2,false>& rhs );

template< typename T1, typename T2, bool SO >
inline bool operator==( const DenseMatrix<T1,SO>& lhs, const SparseMatrix<T2,true>& rhs );

template< typename T1, bool SO1, typename T2, bool SO2 >
inline bool operator==( const SparseMatrix<T1,SO1>& lhs, const DenseMatrix<T2,SO2>& rhs );

template< typename T1, typename T2 >
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator==( const DenseMatrix<T1,false>& mat, T2 scalar );

template< typename T1, typename T2 >
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator==( const DenseMatrix<T1,true>& mat, T2 scalar );

template< typename T1, typename T2, bool SO >
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator==( T1 scalar, const DenseMatrix<T2,SO>& mat );

template< typename T1, bool SO1, typename T2, bool SO2 >
inline bool operator!=( const DenseMatrix<T1,SO1>& lhs, const DenseMatrix<T2,SO2>& rhs );

template< typename T1, bool SO1, typename T2, bool SO2 >
inline bool operator!=( const DenseMatrix<T1,SO1>& lhs, const SparseMatrix<T2,SO2>& rhs );

template< typename T1, bool SO1, typename T2, bool SO2 >
inline bool operator!=( const SparseMatrix<T1,SO1>& lhs, const DenseMatrix<T2,SO2>& rhs );

template< typename T1, typename T2, bool SO >
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator!=( const DenseMatrix<T1,SO>& mat, T2 scalar );

template< typename T1, typename T2, bool SO >
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator!=( T1 scalar, const DenseMatrix<T2,SO>& mat );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two rwo-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side matrix for the comparison.
// \param rhs The right-hand side matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side dense matrix
inline bool operator==( const DenseMatrix<T1,false>& lhs, const DenseMatrix<T2,false>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t i=0; i<A.rows(); ++i ) {
      for( size_t j=0; j<A.columns(); ++j ) {
         if( !equal( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two column-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side matrix for the comparison.
// \param rhs The right-hand side matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side dense matrix
inline bool operator==( const DenseMatrix<T1,true>& lhs, const DenseMatrix<T2,true>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t j=0; j<A.columns(); ++j ) {
      for( size_t i=0; i<A.rows(); ++i ) {
         if( !equal( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two dense matrices with different storage order.
// \ingroup dense_matrix
//
// \param lhs The left-hand side matrix for the comparison.
// \param rhs The right-hand side matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order
inline bool operator==( const DenseMatrix<T1,SO>& lhs, const DenseMatrix<T2,!SO>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   const size_t rows   ( A.rows() );
   const size_t columns( A.columns() );
   const size_t block  ( 16 );

   for( size_t ii=0; ii<rows; ii+=block ) {
      const size_t iend( ( rows < ii+block )?( rows ):( ii+block ) );
      for( size_t jj=0; jj<columns; jj+=block ) {
         const size_t jend( ( columns < jj+block )?( columns ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               if( !equal( A(i,j), B(i,j) ) ) return false;
            }
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a dense matrix and a row-major sparse matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side row-major sparse matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side sparse matrix
        , bool SO >    // Storage order of the left-hand side dense matrix
inline bool operator==( const DenseMatrix<T1,SO>& lhs, const SparseMatrix<T2,false>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;
   typedef typename RemoveReference<CT2>::Type::ConstIterator  ConstIterator;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the dense matrix and sparse matrix operand
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   size_t j( 0 );

   for( size_t i=0; i<B.rows(); ++i ) {
      j = 0;
      for( ConstIterator element=B.begin(i); element!=B.end(i); ++element, ++j ) {
         for( ; j<element->index(); ++j ) {
            if( !isDefault( A(i,j) ) ) return false;
         }
         if( !equal( element->value(), A(i,j) ) ) return false;
      }
      for( ; j<A.columns(); ++j ) {
         if( !isDefault( A(i,j) ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a dense matrix and a column-major sparse matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side column-major sparse matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side sparse matrix
        , bool SO >    // Storage order of the left-hand side dense matrix
inline bool operator==( const DenseMatrix<T1,SO>& lhs, const SparseMatrix<T2,true>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;
   typedef typename RemoveReference<CT2>::Type::ConstIterator  ConstIterator;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the dense matrix and sparse matrix operand
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   size_t i( 0 );

   for( size_t j=0; j<B.columns(); ++j ) {
      i = 0;
      for( ConstIterator element=B.begin(j); element!=B.end(j); ++element, ++i ) {
         for( ; i<element->index(); ++i ) {
            if( !isDefault( A(i,j) ) ) return false;
         }
         if( !equal( element->value(), A(i,j) ) ) return false;
      }
      for( ; i<A.rows(); ++i ) {
         if( !isDefault( A(i,j) ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a sparse matrix and a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side sparse matrix
        , bool SO1     // Storage order of the left-hand side sparse matrix
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO2 >   // Storage order of the right-hand side sparse matrix
inline bool operator==( const SparseMatrix<T1,SO1>& lhs, const DenseMatrix<T2,SO2>& rhs )
{
   return ( rhs == lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a row-major dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side row-major dense matrix for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
//
// If all values of the matrix are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side scalar
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator==( const DenseMatrix<T1,false>& mat, T2 scalar )
{
   typedef typename T1::CompositeType  CT1;

   // Evaluation of the dense matrix operand
   CT1 A( ~mat );

   // In order to compare the matrix and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   for( size_t i=0; i<A.rows(); ++i ) {
      for( size_t j=0; j<A.columns(); ++j ) {
         if( !equal( A(i,j), scalar ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a column-major dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side column-major dense matrix for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
//
// If all values of the matrix are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side scalar
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator==( const DenseMatrix<T1,true>& mat, T2 scalar )
{
   typedef typename T1::CompositeType  CT1;

   // Evaluation of the dense matrix operand
   CT1 A( ~mat );

   // In order to compare the matrix and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   for( size_t j=0; j<A.columns(); ++j ) {
      for( size_t i=0; i<A.rows(); ++i ) {
         if( !equal( A(i,j), scalar ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a scalar value and a dense matrix.
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the comparison.
// \param mat The right-hand side dense matrix for the comparison.
// \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
//
// If all values of the matrix are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsNumeric<T1>, bool >::Type
   operator==( T1 scalar, const DenseMatrix<T2,SO>& mat )
{
   return ( mat == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are not equal, \a false if they are equal.
*/
template< typename T1  // Type of the left-hand side dense matrix
        , bool SO1     // Storage order of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO2 >   // Storage order of the right-hand side dense matrix
inline bool operator!=( const DenseMatrix<T1,SO1>& lhs, const DenseMatrix<T2,SO2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a dense matrix and a sparse matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side sparse matrix for the comparison.
// \return \a true if the two matrices are not equal, \a false if they are equal.
*/
template< typename T1  // Type of the left-hand side dense matrix
        , bool SO1     // Storage order of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side sparse matrix
        , bool SO2 >   // Storage order of the right-hand side sparse matrix
inline bool operator!=( const DenseMatrix<T1,SO1>& lhs, const SparseMatrix<T2,SO2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a sparse matrix and a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are not equal, \a false if they are equal.
*/
template< typename T1  // Type of the left-hand side sparse matrix
        , bool SO1     // Storage order of the left-hand side sparse matrix
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO2 >   // Storage order right-hand side dense matrix
inline bool operator!=( const SparseMatrix<T1,SO1>& lhs, const DenseMatrix<T2,SO2>& rhs )
{
   return !( rhs == lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if at least one element of the matrix is different from the scalar, \a false if not.
//
// If one value of the matrix is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side scalar
        , bool SO >    // Storage order
inline typename EnableIf< IsNumeric<T2>, bool >::Type
   operator!=( const DenseMatrix<T1,SO>& mat, T2 scalar )
{
   return !( mat == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a scalar value and a dense matrix.
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the comparison.
// \param mat The right-hand side dense matrix for the comparison.
// \return \a true if at least one element of the matrix is different from the scalar, \a false if not.
//
// If one value of the matrix is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order
inline typename EnableIf< IsNumeric<T1>, bool >::Type
   operator!=( T1 scalar, const DenseMatrix<T2,SO>& mat )
{
   return !( mat == scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseMatrix functions */
//@{
template< typename MT, bool SO >
bool isnan( const DenseMatrix<MT,SO>& dm );

template< typename MT, bool SO >
bool isDiagonal( const DenseMatrix<MT,SO>& dm );

template< typename MT, bool SO >
bool isSymmetric( const DenseMatrix<MT,SO>& dm );

template< typename MT, bool SO >
const typename MT::ElementType min( const DenseVector<MT,SO>& dm );

template< typename MT, bool SO >
const typename MT::ElementType max( const DenseVector<MT,SO>& dm );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense matrix for not-a-number elements.
// \ingroup dense_matrix
//
// \param dm The matrix to be checked for not-a-number elements.
// \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
//
// This function checks the dense matrix for not-a-number (NaN) elements. If at least one
// element of the matrix is not-a-number, the function returns \a true, otherwise it returns
// \a false.

   \code
   blaze::DynamicMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isnan( A ) ) { ... }
   \endcode

// Note that this function only works for matrices with floating point elements. The attempt to
// use it for a matrix with a non-floating point element type results in a compile time error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
bool isnan( const DenseMatrix<MT,SO>& dm )
{
   typedef typename MT::CompositeType  CT;

   CT A( ~dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<A.columns(); ++j )
            if( isnan( A(i,j) ) ) return true;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( size_t i=0UL; i<A.rows(); ++i )
            if( isnan( A(i,j) ) ) return true;
      }
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the give dense matrix is diagonal.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is diagonal, \a false if not.
//
// This function tests whether the matrix is diagonal, i.e. if the non-diagonal elements are
// default elements. In case of integral or floating point data types, a diagonal matrix has
// the form

                        \f[\left(\begin{array}{*{5}{c}}
                        aa     & 0      & 0      & \cdots & 0  \\
                        0      & bb     & 0      & \cdots & 0  \\
                        0      & 0      & cc     & \cdots & 0  \\
                        \vdots & \vdots & \vdots & \ddots & 0  \\
                        0      & 0      & 0      & 0      & xx \\
                        \end{array}\right)\f]
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
bool isDiagonal( const DenseMatrix<MT,SO>& dm )
{
   const size_t rows   ( (~dm).rows()    );
   const size_t columns( (~dm).columns() );

   if( rows != columns ) return false;

   if( SO == rowMajor ) {
      for( size_t i=1UL; i<rows; ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !isDefault( (~dm)(i,j) ) || !isDefault( (~dm)(j,i) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=1UL; j<columns; ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            if( !isDefault( (~dm)(i,j) ) || !isDefault( (~dm)(j,i) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is symmetric.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is symmetric, \a false if not.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
bool isSymmetric( const DenseMatrix<MT,SO>& dm )
{
   const size_t rows   ( (~dm).rows()    );
   const size_t columns( (~dm).columns() );

   if( rows != columns ) return false;

   if( SO == rowMajor ) {
      for( size_t i=1UL; i<rows; ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !equal( (~dm)(i,j), (~dm)(j,i) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=1UL; j<columns; ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            if( !equal( (~dm)(i,j), (~dm)(j,i) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the smallest element of the dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The smallest dense matrix element.
//
// This function returns the smallest element of the given dense matrix. This function can
// only be used for element types that support the smaller-than relationship. In case the
// matrix currently has either 0 rows or 0 columns, the returned value is the default value
// (e.g. 0 in case of fundamental data types).
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
const typename MT::ElementType min( const DenseMatrix<MT,SO>& dm )
{
   using blaze::min;

   typedef typename MT::ElementType    ET;
   typedef typename MT::CompositeType  CT;

   CT A( ~dm );  // Evaluation of the dense matrix operand

   if( A.rows() == 0UL || A.columns() == 0UL ) return ET();

   ET minimum( A(0,0) );

   if( SO == rowMajor ) {
      for( size_t j=1UL; j<A.columns(); ++j )
         minimum = min( minimum, A(0UL,j) );
      for( size_t i=1UL; i<A.rows(); ++i )
         for( size_t j=0UL; j<A.columns(); ++j )
            minimum = min( minimum, A(i,j) );
   }
   else {
      for( size_t i=1UL; i<A.rows(); ++i )
         minimum = min( minimum, A(i,0UL) );
      for( size_t j=1UL; j<A.columns(); ++j )
         for( size_t i=0UL; i<A.rows(); ++i )
            minimum = min( minimum, A(i,j) );
   }

   return minimum;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the largest element of the dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix.
// \return The largest dense matrix element.
//
// This function returns the largest element of the given dense matrix. This function can
// only be used for element types that support the smaller-than relationship. In case the
// matrix currently has either 0 rows or 0 columns, the returned value is the default value
// (e.g. 0 in case of fundamental data types).
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Transpose flag
const typename MT::ElementType max( const DenseMatrix<MT,SO>& dm )
{
   using blaze::max;

   typedef typename MT::ElementType    ET;
   typedef typename MT::CompositeType  CT;

   CT A( ~dm );  // Evaluation of the dense matrix operand

   if( A.rows() == 0UL || A.columns() == 0UL ) return ET();

   ET maximum( A(0,0) );

   if( SO == rowMajor ) {
      for( size_t j=1UL; j<A.columns(); ++j )
         maximum = max( maximum, A(0UL,j) );
      for( size_t i=1UL; i<A.rows(); ++i )
         for( size_t j=0UL; j<A.columns(); ++j )
            maximum = max( maximum, A(i,j) );
   }
   else {
      for( size_t i=1UL; i<A.rows(); ++i )
         maximum = max( maximum, A(i,0UL) );
      for( size_t j=1UL; j<A.columns(); ++j )
         for( size_t i=0UL; i<A.rows(); ++i )
            maximum = max( maximum, A(i,j) );
   }

   return maximum;
}
//*************************************************************************************************

} // namespace blaze

#endif
