//=================================================================================================
/*!
//  \file blaze/math/adaptors/LowerMatrix.h
//  \brief Header file for the implementation of a lower matrix adaptor
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

#ifndef _BLAZE_MATH_ADAPTORS_LOWERMATRIX_H_
#define _BLAZE_MATH_ADAPTORS_LOWERMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/adaptors/lowermatrix/BaseTemplate.h>
#include <blaze/math/adaptors/lowermatrix/Dense.h>
#include <blaze/math/adaptors/lowermatrix/Sparse.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/Forward.h>
#include <blaze/math/Functions.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DerestrictTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MathTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/Columns.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/Rows.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/Unused.h>
#include <blaze/util/valuetraits/IsTrue.h>


namespace blaze {

//=================================================================================================
//
//  LOWERMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name LowerMatrix operators */
//@{
template< typename MT, bool SO, bool DF >
inline void reset( LowerMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline void reset( LowerMatrix<MT,SO,DF>& m, size_t i );

template< typename MT, bool SO, bool DF >
inline void clear( LowerMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline bool isDefault( const LowerMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline bool isIntact( const LowerMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline void swap( LowerMatrix<MT,SO,DF>& a, LowerMatrix<MT,SO,DF>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given lower matrix.
// \ingroup lower_matrix
//
// \param m The lower matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( LowerMatrix<MT,SO,DF>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the specified row/column of the given lower matrix.
// \ingroup lower_matrix
//
// \param m The lower matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given lower matrix to their
// default value. In case the given matrix is a \a rowMajor matrix the function resets the values
// in row \a i, if it is a \a columnMajor matrix the function resets the values in column \a i.
// Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( LowerMatrix<MT,SO,DF>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given lower matrix.
// \ingroup lower_matrix
//
// \param m The lower matrix to be cleared.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void clear( LowerMatrix<MT,SO,DF>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given lower matrix is in default state.
// \ingroup lower_matrix
//
// \param m The lower matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
//
// This function checks whether the matrix is in default state. For instance, in case the
// matrix is instantiated for a built-in integral or floating point data type, the function
// returns \a true in case all matrix elements are 0 and \a false in case any matrix element
// is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;

   LowerMatrix< DynamicMatrix<int> > A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isDefault( const LowerMatrix<MT,SO,DF>& m )
{
   return isDefault( m.matrix_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given lower matrix are intact.
// \ingroup lower_matrix
//
// \param m The lower matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the lower matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;

   LowerMatrix< DynamicMatrix<int> > A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isIntact( const LowerMatrix<MT,SO,DF>& m )
{
   return ( isIntact( m.matrix_ ) && isLower( m.matrix_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
// \ingroup lower_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void swap( LowerMatrix<MT,SO,DF>& a, LowerMatrix<MT,SO,DF>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side dense vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side dense vector
inline bool tryAssign( const LowerMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.rows() - row, "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   if( column <= row )
      return true;

   const size_t iend( min( column - row, (~rhs).size() ) );

   for( size_t i=0UL; i<iend; ++i ) {
      if( !isDefault( (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side dense vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side dense vector
inline bool tryAssign( const LowerMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t ibegin( ( row < column )?( 0UL ):( row - column + 1UL ) );

   for( size_t i=ibegin; i<(~rhs).size(); ++i ) {
      if( !isDefault( (~rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse vector to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side sparse vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side sparse vector
inline bool tryAssign( const LowerMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.rows() - row, "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   typedef typename VT::ConstIterator  RhsIterator;

   if( column <= row )
      return true;

   const RhsIterator last( (~rhs).lowerBound( column - row ) );

   for( RhsIterator element=(~rhs).begin(); element!=last; ++element ) {
      if( !isDefault( element->value() ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse vector to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side sparse vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side sparse vector
inline bool tryAssign( const LowerMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   typedef typename VT::ConstIterator  RhsIterator;

   const RhsIterator last( (~rhs).end() );
   RhsIterator element( (~rhs).lowerBound( ( row < column )?( 0UL ):( row - column + 1UL ) ) );

   for( ; element!=last; ++element ) {
      if( !isDefault( element->value() ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense matrix to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side dense matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side dense matrix
inline bool tryAssign( const LowerMatrix<MT1,SO,DF>& lhs,
                       const DenseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( row + 1UL >= column + N )
      return true;

   const size_t iend( min( column + N - row - 1UL, M ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row + i >= column );
      const size_t jbegin( ( containsDiagonal )?( row + i - column + 1UL ):( 0UL ) );

      for( size_t j=jbegin; j<N; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense matrix to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side dense matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side dense matrix
inline bool tryAssign( const LowerMatrix<MT1,SO,DF>& lhs,
                       const DenseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( row + 1UL >= column + N )
      return true;

   const size_t jbegin( ( row < column )?( 0UL ):( row - column + 1UL ) );

   for( size_t j=jbegin; j<N; ++j )
   {
      const size_t iend( min( column + j - row, M ) );

      for( size_t i=0UL; i<iend; ++i ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side sparse matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline bool tryAssign( const LowerMatrix<MT1,SO,DF>& lhs,
                       const SparseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   typedef typename MT2::ConstIterator  RhsIterator;

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( row + 1UL >= column + N )
      return true;

   const size_t iend( min( column + N - row - 1UL, M ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row + i >= column );
      const size_t index( ( containsDiagonal )?( row + i - column + 1UL ):( 0UL ) );

      const RhsIterator last( (~rhs).end(i) );
      RhsIterator element( (~rhs).lowerBound( i, index ) );

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side sparse matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline bool tryAssign( const LowerMatrix<MT1,SO,DF>& lhs,
                       const SparseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).rows() <= lhs.rows() - row, "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( (~rhs).columns() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   typedef typename MT2::ConstIterator  RhsIterator;

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( row + 1UL >= column + N )
      return true;

   const size_t jbegin( ( row < column )?( 0UL ):( row - column + 1UL ) );

   for( size_t j=jbegin; j<N; ++j )
   {
      const size_t index( column + j - row );
      const RhsIterator last( (~rhs).lowerBound( min( index, M ), j ) );

      for( RhsIterator element=(~rhs).begin(j); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const LowerMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
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
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryAddAssign( const LowerMatrix<MT1,SO1,DF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
// \param rhs The right-hand side vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool trySubAssign( const LowerMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a lower matrix.
// \ingroup lower_matrix
//
// \param lhs The target left-hand side lower matrix.
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
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool trySubAssign( const LowerMatrix<MT1,SO1,DF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the instance without the access restrictions to the upper part.
// \ingroup lower_matrix
//
// \param m The lower matrix to be derestricted.
// \return Reference to the matrix without access restrictions.
//
// This function returns a reference to the given lower matrix instance that has no access
// restrictions to the upper part of the matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline MT& derestrict( LowerMatrix<MT,SO,DF>& m )
{
   return m.matrix_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct Rows< LowerMatrix<MT,SO,DF> > : public Rows<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct Columns< LowerMatrix<MT,SO,DF> > : public Columns<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSQUARE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsSquare< LowerMatrix<MT,SO,DF> > : public IsTrue<true>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsLower< LowerMatrix<MT,SO,DF> > : public IsTrue<true>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsAdaptor< LowerMatrix<MT,SO,DF> > : public IsTrue<true>
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
template< typename MT, bool SO, bool DF >
struct IsRestricted< LowerMatrix<MT,SO,DF> > : public IsTrue<true>
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
template< typename MT, bool SO >
struct HasConstDataAccess< LowerMatrix<MT,SO,true> > : public IsTrue<true>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsResizable< LowerMatrix<MT,SO,DF> > : public IsTrue< IsResizable<MT>::value >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REMOVEADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct RemoveAdaptor< LowerMatrix<MT,SO,DF> >
{
   typedef MT  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DERESTRICTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DerestrictTrait< LowerMatrix<MT,SO,DF> >
{
   typedef MT&  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct AddTrait< LowerMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename AddTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< StaticMatrix<T,M,N,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct AddTrait< LowerMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename AddTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< HybridMatrix<T,M,N,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct AddTrait< LowerMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   typedef typename AddTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< DynamicMatrix<T,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct AddTrait< LowerMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   typedef typename AddTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< CompressedMatrix<T,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct AddTrait< LowerMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct AddTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< LowerMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< HermitianMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< LowerMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef LowerMatrix< typename AddTrait<MT1,MT2>::Type >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct SubTrait< LowerMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename SubTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< StaticMatrix<T,M,N,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct SubTrait< LowerMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename SubTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< HybridMatrix<T,M,N,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SubTrait< LowerMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   typedef typename SubTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< DynamicMatrix<T,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SubTrait< LowerMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   typedef typename SubTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< CompressedMatrix<T,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct SubTrait< LowerMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct SubTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< LowerMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< HermitianMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< LowerMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef LowerMatrix< typename SubTrait<MT1,MT2>::Type >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, typename T >
struct MultTrait< LowerMatrix<MT,SO,DF>, T >
{
   typedef LowerMatrix< typename MultTrait<MT,T>::Type >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< T, LowerMatrix<MT,SO,DF> >
{
   typedef LowerMatrix< typename MultTrait<T,MT>::Type >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
};

template< typename MT, bool SO, bool DF, typename T, size_t N >
struct MultTrait< LowerMatrix<MT,SO,DF>, StaticVector<T,N,false> >
{
   typedef typename MultTrait< MT, StaticVector<T,N,false> >::Type  Type;
};

template< typename T, size_t N, typename MT, bool SO, bool DF >
struct MultTrait< StaticVector<T,N,true>, LowerMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< StaticVector<T,N,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, typename T, size_t N >
struct MultTrait< LowerMatrix<MT,SO,DF>, HybridVector<T,N,false> >
{
   typedef typename MultTrait< MT, HybridVector<T,N,false> >::Type  Type;
};

template< typename T, size_t N, typename MT, bool SO, bool DF >
struct MultTrait< HybridVector<T,N,true>, LowerMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< HybridVector<T,N,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, typename T >
struct MultTrait< LowerMatrix<MT,SO,DF>, DynamicVector<T,false> >
{
   typedef typename MultTrait< MT, DynamicVector<T,false> >::Type  Type;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< DynamicVector<T,true>, LowerMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< DynamicVector<T,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, typename T >
struct MultTrait< LowerMatrix<MT,SO,DF>, CompressedVector<T,false> >
{
   typedef typename MultTrait< MT, CompressedVector<T,false> >::Type  Type;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< CompressedVector<T,true>, LowerMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< CompressedVector<T,true>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< LowerMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename MultTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< StaticMatrix<T,M,N,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< LowerMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename MultTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< HybridMatrix<T,M,N,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct MultTrait< LowerMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   typedef typename MultTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< DynamicMatrix<T,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct MultTrait< LowerMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   typedef typename MultTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< CompressedMatrix<T,SO1>, LowerMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct MultTrait< LowerMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct MultTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< LowerMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< HermitianMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< LowerMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef LowerMatrix< typename MultTrait<MT1,MT2>::Type >  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, typename T >
struct DivTrait< LowerMatrix<MT,SO,DF>, T >
{
   typedef LowerMatrix< typename DivTrait<MT,T>::Type >  Type;
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MATHTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MathTrait< LowerMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef LowerMatrix< typename MathTrait<MT1,MT2>::HighType >  HighType;
   typedef LowerMatrix< typename MathTrait<MT1,MT2>::LowType  >  LowType;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct SubmatrixTrait< LowerMatrix<MT,SO,DF> >
{
   typedef typename SubmatrixTrait<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct RowTrait< LowerMatrix<MT,SO,DF> >
{
   typedef typename RowTrait<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct ColumnTrait< LowerMatrix<MT,SO,DF> >
{
   typedef typename ColumnTrait<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
