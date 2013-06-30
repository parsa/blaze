//=================================================================================================
/*!
//  \file blaze/math/SparseMatrix.h
//  \brief Header file for the SparseMatrix CRTP base class
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

#ifndef _BLAZE_MATH_SPARSEMATRIX_H_
#define _BLAZE_MATH_SPARSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DVecTSVecMultExpr.h>
#include <blaze/math/expressions/SMatAbsExpr.h>
#include <blaze/math/expressions/SMatDVecMultExpr.h>
#include <blaze/math/expressions/SMatEvalExpr.h>
#include <blaze/math/expressions/SMatScalarDivExpr.h>
#include <blaze/math/expressions/SMatScalarMultExpr.h>
#include <blaze/math/expressions/SMatSMatAddExpr.h>
#include <blaze/math/expressions/SMatSMatMultExpr.h>
#include <blaze/math/expressions/SMatSMatSubExpr.h>
#include <blaze/math/expressions/SMatSVecMultExpr.h>
#include <blaze/math/expressions/SMatTransExpr.h>
#include <blaze/math/expressions/SMatTSMatAddExpr.h>
#include <blaze/math/expressions/SMatTSMatMultExpr.h>
#include <blaze/math/expressions/SMatTSMatSubExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/SVecTDVecMultExpr.h>
#include <blaze/math/expressions/SVecTSVecMultExpr.h>
#include <blaze/math/expressions/TDVecSMatMultExpr.h>
#include <blaze/math/expressions/TDVecTSMatMultExpr.h>
#include <blaze/math/expressions/TSMatDVecMultExpr.h>
#include <blaze/math/expressions/TSMatSMatMultExpr.h>
#include <blaze/math/expressions/TSMatSMatSubExpr.h>
#include <blaze/math/expressions/TSMatSVecMultExpr.h>
#include <blaze/math/expressions/TSMatTSMatAddExpr.h>
#include <blaze/math/expressions/TSMatTSMatMultExpr.h>
#include <blaze/math/expressions/TSMatTSMatSubExpr.h>
#include <blaze/math/expressions/TSVecSMatMultExpr.h>
#include <blaze/math/expressions/TSVecTSMatMultExpr.h>
#include <blaze/math/Matrix.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseMatrix operators */
//@{
template< typename T1, typename T2, bool SO >
inline bool operator==( const SparseMatrix<T1,false>& lhs, const SparseMatrix<T2,false>& rhs );

template< typename T1, typename T2, bool SO >
inline bool operator==( const SparseMatrix<T1,true>& lhs, const SparseMatrix<T2,true>& rhs );

template< typename T1, typename T2, bool SO >
inline bool operator==( const SparseMatrix<T1,SO>& lhs, const SparseMatrix<T2,!SO>& rhs );

template< typename T1, bool SO1, typename T2, bool SO2 >
inline bool operator!=( const SparseMatrix<T1,SO1>& lhs, const SparseMatrix<T2,SO2>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two row-major sparse matrices.
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side sparse matrix for the comparison.
// \return \a true if the two sparse matrices are equal, \a false if not.
*/
template< typename T1    // Type of the left-hand side sparse matrix
        , typename T2 >  // Type of the right-hand side sparse matrix
inline bool operator==( const SparseMatrix<T1,false>& lhs, const SparseMatrix<T2,false>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;
   typedef typename RemoveReference<CT1>::Type::ConstIterator  LhsConstIterator;
   typedef typename RemoveReference<CT2>::Type::ConstIterator  RhsConstIterator;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two sparse matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t i=0; i<A.rows(); ++i )
   {
      const LhsConstIterator lend( A.end(i) );
      const RhsConstIterator rend( B.end(i) );

      LhsConstIterator lelem( A.begin(i) );
      RhsConstIterator relem( B.begin(i) );

      while( lelem < lend && relem < rend ) {
         if( lelem->index() < relem->index() && !isDefault( (lelem++)->value() ) )
            return false;
         else if( lelem->index() > relem->index() && !isDefault( (relem++)->value() ) )
            return false;
         else if( !equal( (lelem++)->value(), (relem++)->value() ) )
            return false;
      }

      while( lelem < lend ) {
         if( !isDefault( (lelem++)->value() ) )
            return false;
      }

      while( relem < rend ) {
         if( !isDefault( (relem++)->value() ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two column-major sparse matrices.
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side sparse matrix for the comparison.
// \return \a true if the two sparse matrices are equal, \a false if not.
*/
template< typename T1    // Type of the left-hand side sparse matrix
        , typename T2 >  // Type of the right-hand side sparse matrix
inline bool operator==( const SparseMatrix<T1,true>& lhs, const SparseMatrix<T2,true>& rhs )
{
   typedef typename T1::CompositeType  CT1;
   typedef typename T2::CompositeType  CT2;
   typedef typename RemoveReference<CT1>::Type::ConstIterator  LhsConstIterator;
   typedef typename RemoveReference<CT2>::Type::ConstIterator  RhsConstIterator;

   // Early exit in case the matrix sizes don't match
   if( (~lhs).rows() != (~rhs).rows() || (~lhs).columns() != (~rhs).columns() )
      return false;

   // Evaluation of the two sparse matrix operands
   CT1 A( ~lhs );
   CT2 B( ~rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t j=0; j<A.columns(); ++j )
   {
      const LhsConstIterator lend( A.end(j) );
      const RhsConstIterator rend( B.end(j) );

      LhsConstIterator lelem( A.begin(j) );
      RhsConstIterator relem( B.begin(j) );

      while( lelem < lend && relem < rend ) {
         if( lelem->index() < relem->index() && !isDefault( (lelem++)->value() ) )
            return false;
         else if( lelem->index() > relem->index() && !isDefault( (relem++)->value() ) )
            return false;
         else if( !equal( (lelem++)->value(), (relem++)->value() ) )
            return false;
      }

      while( lelem < lend ) {
         if( !isDefault( (lelem++)->value() ) )
            return false;
      }

      while( relem < rend ) {
         if( !isDefault( (relem++)->value() ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two sparse matrices with different storage order.
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side sparse matrix for the comparison.
// \return \a true if the two sparse matrices are equal, \a false if not.
*/
template< typename T1  // Type of the left-hand side sparse matrix
        , typename T2  // Type of the right-hand side sparse matrix
        , bool SO >    // Storage order
inline bool operator==( const SparseMatrix<T1,SO>& lhs, const SparseMatrix<T2,!SO>& rhs )
{
   const typename T2::TransposeType tmp( ~rhs );
   return ( ~lhs == tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two sparse matrices.
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side sparse matrix for the comparison.
// \return \a true if the two sparse matrices are not equal, \a false if they are equal.
*/
template< typename T1  // Type of the left-hand side sparse matrix
        , bool SO1     // Storage order of the left-hand side sparse matrix
        , typename T2  // Type of the right-hand side sparse matrix
        , bool SO2 >   // Storage order of the right-hand side sparse matrix
inline bool operator!=( const SparseMatrix<T1,SO1>& lhs, const SparseMatrix<T2,SO2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
