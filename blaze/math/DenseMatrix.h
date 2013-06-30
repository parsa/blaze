//=================================================================================================
/*!
//  \file blaze/math/DenseMatrix.h
//  \brief Header file for the DenseMatrix CRTP base class
//
//  Copyright (C) 2011 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. This library is free software; you can redistribute
//  it and/or modify it under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 3, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with a special
//  exception for linking and compiling against the Blaze library, the so-called "runtime
//  exception"; see the file COPYING. If not, see http://www.gnu.org/licenses/.
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

} // namespace blaze

#endif
