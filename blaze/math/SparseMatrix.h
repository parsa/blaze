//=================================================================================================
/*!
//  \file blaze/math/SparseMatrix.h
//  \brief Header file for the SparseMatrix CRTP base class
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
#include <blaze/math/shims/IsNaN.h>
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

      while( lelem != lend && relem != rend ) {
         if( lelem->index() < relem->index() && !isDefault( (lelem++)->value() ) )
            return false;
         else if( lelem->index() > relem->index() && !isDefault( (relem++)->value() ) )
            return false;
         else if( !equal( (lelem++)->value(), (relem++)->value() ) )
            return false;
      }

      while( lelem != lend ) {
         if( !isDefault( (lelem++)->value() ) )
            return false;
      }

      while( relem != rend ) {
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

      while( lelem != lend && relem != rend ) {
         if( lelem->index() < relem->index() && !isDefault( (lelem++)->value() ) )
            return false;
         else if( lelem->index() > relem->index() && !isDefault( (relem++)->value() ) )
            return false;
         else if( !equal( (lelem++)->value(), (relem++)->value() ) )
            return false;
      }

      while( lelem != lend ) {
         if( !isDefault( (lelem++)->value() ) )
            return false;
      }

      while( relem != rend ) {
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




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseMatrix functions */
//@{
template< typename MT, bool SO >
bool isnan( const SparseMatrix<MT,SO>& sm );

template< typename MT, bool SO >
bool isDiagonal( const SparseMatrix<MT,SO>& sm );

template< typename MT, bool SO >
bool isSymmetric( const SparseMatrix<MT,SO>& sm );

template< typename MT, bool SO >
const typename MT::ElementType min( const SparseVector<MT,SO>& sm );

template< typename MT, bool SO >
const typename MT::ElementType max( const SparseVector<MT,SO>& sm );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse matrix for not-a-number elements.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked for not-a-number elements.
// \return \a true if at least one element of the sparse matrix is not-a-number, \a false otherwise.
//
// This function checks the sparse matrix for not-a-number (NaN) elements. If at least one
// element of the matrix is not-a-number, the function returns \a true, otherwise it returns
// \a false.

   \code
   blaze::CompressedMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isnan( A ) ) { ... }
   \endcode

// Note that this function only works for matrices with floating point elements. The attempt to
// use it for a matrix with a non-floating point element type results in a compile time error.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
bool isnan( const SparseMatrix<MT,SO>& sm )
{
   typedef typename MT::CompositeType  CT;
   typedef typename RemoveReference<CT>::Type::ConstIterator  ConstIterator;

   CT A( ~sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( ConstIterator element=A.begin(i); element!=A.end(i); ++element )
            if( isnan( element->value() ) ) return true;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( ConstIterator element=A.begin(j); element!=A.end(j); ++element )
            if( isnan( element->value() ) ) return true;
      }
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the give sparse matrix is diagonal.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
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
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
bool isDiagonal( const SparseMatrix<MT,SO>& sm )
{
   typedef typename MT::ConstIterator  ConstIterator;

   const size_t rows   ( (~sm).rows()    );
   const size_t columns( (~sm).columns() );

   if( rows != columns ) return false;

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<rows; ++i ) {
         for( ConstIterator element=(~sm).begin(i); element!=(~sm).end(i); ++element )
            if( element->index() != i && !isDefault( element->value() ) )
               return false;
      }
   }
   else {
      for( size_t j=0UL; j<columns; ++j ) {
         for( ConstIterator element=(~sm).begin(j); element!=(~sm).end(j); ++element )
            if( element->index() != j && !isDefault( element->value() ) )
               return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is symmetric.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is symmetric, \a false if not.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
bool isSymmetric( const SparseMatrix<MT,SO>& sm )
{
   typedef typename MT::ConstIterator  ConstIterator;

   const size_t rows   ( (~sm).rows()    );
   const size_t columns( (~sm).columns() );

   if( rows != columns ) return false;

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<rows; ++i ) {
         for( ConstIterator element=(~sm).begin(i); element!=(~sm).end(i); ++element )
         {
            const size_t index( element->index() );

            if( isDefault( element->value() ) )
               continue;

            const ConstIterator pos( (~sm).lowerBound( index, i ) );
            if( pos == (~sm).end(index) || pos->index() != i || !equal( pos->value(), element->value() ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<columns; ++j ) {
         for( ConstIterator element=(~sm).begin(j); element!=(~sm).end(j); ++element )
         {
            const size_t index( element->index() );

            if( isDefault( element->value() ) )
               continue;

            const ConstIterator pos( (~sm).lowerBound( j, index ) );
            if( pos == (~sm).end(index) || pos->index() != j || !equal( pos->value(), element->value() ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the smallest element of the sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix.
// \return The smallest sparse matrix element.
//
// This function returns the smallest element of the given sparse matrix. This function can
// only be used for element types that support the smaller-than relationship. In case the
// matrix currently has either 0 rows or 0 columns, the returned value is the default value
// (e.g. 0 in case of fundamental data types).
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
const typename MT::ElementType min( const SparseMatrix<MT,SO>& sm )
{
   using blaze::min;

   typedef typename MT::ElementType    ET;
   typedef typename MT::CompositeType  CT;
   typedef typename RemoveReference<CT>::Type::ConstIterator  ConstIterator;

   CT A( ~sm );  // Evaluation of the sparse matrix operand

   const size_t nonzeros( A.nonZeros() );

   if( nonzeros == 0UL ) {
      return ET();
   }

   ET minimum = ET();
   if( nonzeros == A.rows() * A.columns() ) {
      minimum = A.begin( 0UL )->value();
   }

   const size_t index( ( SO == rowMajor )?( A.rows() ):( A.columns() ) );

   for( size_t i=0UL; i<index; ++i ) {
      const ConstIterator end( A.end( i ) );
      ConstIterator element( A.begin( i ) );
      for( ; element!=end; ++element )
         minimum = min( minimum, element->value() );
   }

   return minimum;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the largest element of the sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix.
// \return The largest sparse matrix element.
//
// This function returns the largest element of the given sparse matrix. This function can
// only be used for element types that support the smaller-than relationship. In case the
// matrix currently has either 0 rows or 0 columns, the returned value is the default value
// (e.g. 0 in case of fundamental data types).
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
const typename MT::ElementType max( const SparseMatrix<MT,SO>& sm )
{
   using blaze::max;

   typedef typename MT::ElementType    ET;
   typedef typename MT::CompositeType  CT;
   typedef typename RemoveReference<CT>::Type::ConstIterator  ConstIterator;

   CT A( ~sm );  // Evaluation of the sparse matrix operand

   const size_t nonzeros( A.nonZeros() );

   if( nonzeros == 0UL ) {
      return ET();
   }

   ET maximum = ET();
   if( nonzeros == A.rows() * A.columns() ) {
      maximum = A.begin( 0UL )->value();
   }

   const size_t index( ( SO == rowMajor )?( A.rows() ):( A.columns() ) );

   for( size_t i=0UL; i<index; ++i ) {
      const ConstIterator end( A.end( i ) );
      ConstIterator element( A.begin( i ) );
      for( ; element!=end; ++element )
         maximum = max( maximum, element->value() );
   }

   return maximum;
}
//*************************************************************************************************

} // namespace blaze

#endif
