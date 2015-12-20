//=================================================================================================
/*!
//  \file blaze/math/adaptors/StrictlyUpperMatrix.h
//  \brief Header file for the implementation of a strictly upper triangular matrix adaptor
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

#ifndef _BLAZE_MATH_ADAPTORS_STRICTLYUPPERMATRIX_H_
#define _BLAZE_MATH_ADAPTORS_STRICTLYUPPERMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/adaptors/strictlyuppermatrix/BaseTemplate.h>
#include <blaze/math/adaptors/strictlyuppermatrix/Dense.h>
#include <blaze/math/adaptors/strictlyuppermatrix/Sparse.h>
#include <blaze/math/adaptors/uniuppermatrix/BaseTemplate.h>
#include <blaze/math/adaptors/uppermatrix/BaseTemplate.h>
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
#include <blaze/math/typetraits/Columns.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/Rows.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/Unused.h>
#include <blaze/util/valuetraits/IsTrue.h>


namespace blaze {

//=================================================================================================
//
//  STRICTLYUPPERMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name StrictlyUpperMatrix operators */
//@{
template< typename MT, bool SO, bool DF >
inline void reset( StrictlyUpperMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline void reset( StrictlyUpperMatrix<MT,SO,DF>& m, size_t i );

template< typename MT, bool SO, bool DF >
inline void clear( StrictlyUpperMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline bool isDefault( const StrictlyUpperMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline bool isIntact( const StrictlyUpperMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline void swap( StrictlyUpperMatrix<MT,SO,DF>& a, StrictlyUpperMatrix<MT,SO,DF>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given strictly upper matrix.
// \ingroup strictly_upper_matrix
//
// \param m The strictly upper matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( StrictlyUpperMatrix<MT,SO,DF>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the specified row/column of the given strictly upper matrix.
// \ingroup strictly_upper_matrix
//
// \param m The strictly upper matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given strictly upper matrix
// to their default value. In case the given matrix is a \a rowMajor matrix the function resets
// the values in row \a i, if it is a \a columnMajor matrix the function resets the values in
// column \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( StrictlyUpperMatrix<MT,SO,DF>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given strictly upper matrix.
// \ingroup strictly_upper_matrix
//
// \param m The strictly upper matrix to be cleared.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void clear( StrictlyUpperMatrix<MT,SO,DF>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given resizable strictly upper matrix is in default state.
// \ingroup strictly_upper_matrix
//
// \param m The strictly upper matrix to be tested for its default state.
// \return \a true in case the given matrix is in default state, \a false otherwise.
//
// This function checks whether the resizable strictly upper triangular matrix is in default
// state.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isDefault_backend( const StrictlyUpperMatrix<MT,SO,DF>& m, TrueType )
{
   return ( m.rows() == 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given fixed-size strictly upper matrix is in default state.
// \ingroup strictly_upper_matrix
//
// \param m The strictly upper matrix to be tested for its default state.
// \return \a true in case the given matrix is in default state, \a false otherwise.
//
// This function checks whether the fixed-size strictly upper triangular matrix is in default
// state.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isDefault_backend( const StrictlyUpperMatrix<MT,SO,DF>& m, FalseType )
{
   return isIdentity( m );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given strictly upper matrix is in default state.
// \ingroup strictly_upper_matrix
//
// \param m The strictly upper matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
//
// This function checks whether the strictly upper triangular matrix is in default state. The
// following example demonstrates the use of the \a isDefault function:

   \code
   using blaze::DynamicMatrix;
   using blaze::StrictlyUpperMatrix;
   using blaze::rowMajor;

   StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> > A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isDefault( const StrictlyUpperMatrix<MT,SO,DF>& m )
{
   return isDefault_backend( m, typename IsResizable<MT>::Type() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given strictly upper matrix are intact.
// \ingroup strictly_upper_matrix
//
// \param m The strictly upper matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the strictly upper matrix are intact, i.e.
// if its state is valid. In case the invariants are intact, the function returns \a true, else
// it will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   using blaze::DynamicMatrix;
   using blaze::StrictlyUpperMatrix;

   StrictlyUpperMatrix< DynamicMatrix<int> > A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isIntact( const StrictlyUpperMatrix<MT,SO,DF>& m )
{
   return m.isIntact();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
// \ingroup strictly_upper_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
// \exception no-throw guarantee.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void swap( StrictlyUpperMatrix<MT,SO,DF>& a, StrictlyUpperMatrix<MT,SO,DF>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.rows() - row, "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   const size_t ibegin( ( column <= row )?( 0UL ):( column - row ) );

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
/*!\brief Predict invariant violations by the assignment of a dense vector to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   if( row < column )
      return true;

   const size_t iend( min( row - column + 1UL, (~rhs).size() ) );

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
/*!\brief Predict invariant violations by the assignment of a sparse vector to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.rows() - row, "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   typedef typename VT::ConstIterator  RhsIterator;

   const RhsIterator last( (~rhs).end() );
   RhsIterator element( (~rhs).lowerBound( ( column <= row )?( 0UL ):( column - row ) ) );

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
/*!\brief Predict invariant violations by the assignment of a sparse vector to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( (~rhs).size() <= lhs.columns() - column, "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   typedef typename VT::ConstIterator  RhsIterator;

   if( row < column )
      return true;

   const RhsIterator last( (~rhs).lowerBound( row - column + 1UL ) );

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
/*!\brief Predict invariant violations by the assignment of a dense matrix to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT1,SO,DF>& lhs,
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

   if( column + 1UL >= row + M )
      return true;

   const size_t ibegin( ( column <= row )?( 0UL ):( column - row ) );

   for( size_t i=ibegin; i<M; ++i )
   {
      const size_t jend( min( row + i - column + 1UL, N ) );

      for( size_t j=0UL; j<jend; ++j ) {
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
/*!\brief Predict invariant violations by the assignment of a dense matrix to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT1,SO,DF>& lhs,
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

   if( column + 1UL >= row + M )
      return true;

   const size_t jend( min( row + M - column, N ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column + j >= row );
      const size_t ibegin( ( containsDiagonal )?( column + j - row ):( 0UL ) );

      for( size_t i=ibegin; i<M; ++i ) {
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
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT1,SO,DF>& lhs,
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

   if( column + 1UL >= row + M )
      return true;

   const size_t ibegin( ( column < row )?( 0UL ):( column - row ) );

   for( size_t i=ibegin; i<M; ++i )
   {
      const size_t index( row + i - column + 1UL );
      const RhsIterator last( (~rhs).lowerBound( i, min( index, N ) ) );

      for( RhsIterator element=(~rhs).begin(i); element!=last; ++element ) {
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
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAssign( const StrictlyUpperMatrix<MT1,SO,DF>& lhs,
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

   if( column + 1UL >= row + M )
      return true;

   const size_t jend( min( row + M - column, N ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column + j >= row );
      const size_t index( ( containsDiagonal )?( column + j - row ):( 0UL ) );

      const RhsIterator last( (~rhs).end(j) );
      RhsIterator element( (~rhs).lowerBound( index, j ) );

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
/*!\brief Predict invariant violations by the addition assignment of a vector to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAddAssign( const StrictlyUpperMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a strictly upper
//        matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool tryAddAssign( const StrictlyUpperMatrix<MT1,SO1,DF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a strictly
//        upper matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool trySubAssign( const StrictlyUpperMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a strictly
//        upper matrix.
// \ingroup strictly_upper_matrix
//
// \param lhs The target left-hand side strictly upper matrix.
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
inline bool trySubAssign( const StrictlyUpperMatrix<MT1,SO1,DF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the instance without the access restrictions to the lower part.
// \ingroup math_shims
//
// \param m The strictly upper matrix to be derestricted.
// \return Reference to the matrix without access restrictions.
//
// This function returns a reference to the given strictly upper matrix instance that has no
// access restrictions to the lower part of the matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline MT& derestrict( StrictlyUpperMatrix<MT,SO,DF>& m )
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
struct Rows< StrictlyUpperMatrix<MT,SO,DF> > : public Rows<MT>
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
struct Columns< StrictlyUpperMatrix<MT,SO,DF> > : public Columns<MT>
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
struct IsSquare< StrictlyUpperMatrix<MT,SO,DF> > : public IsTrue<true>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsStrictlyUpper< StrictlyUpperMatrix<MT,SO,DF> > : public IsTrue<true>
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
struct IsAdaptor< StrictlyUpperMatrix<MT,SO,DF> > : public IsTrue<true>
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
struct IsRestricted< StrictlyUpperMatrix<MT,SO,DF> > : public IsTrue<true>
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
struct HasConstDataAccess< StrictlyUpperMatrix<MT,SO,true> > : public IsTrue<true>
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
template< typename MT, bool SO, bool DF >
struct IsAligned< StrictlyUpperMatrix<MT,SO,DF> > : public IsTrue< IsAligned<MT>::value >
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
template< typename MT, bool SO, bool DF >
struct IsPadded< StrictlyUpperMatrix<MT,SO,DF> > : public IsTrue< IsPadded<MT>::value >
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
struct IsResizable< StrictlyUpperMatrix<MT,SO,DF> > : public IsTrue< IsResizable<MT>::value >
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
struct RemoveAdaptor< StrictlyUpperMatrix<MT,SO,DF> >
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
struct DerestrictTrait< StrictlyUpperMatrix<MT,SO,DF> >
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
struct AddTrait< StrictlyUpperMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename AddTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< StaticMatrix<T,M,N,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct AddTrait< StrictlyUpperMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename AddTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< HybridMatrix<T,M,N,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct AddTrait< StrictlyUpperMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   typedef typename AddTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< DynamicMatrix<T,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool AF, bool PF, bool SO2 >
struct AddTrait< StrictlyUpperMatrix<MT,SO1,DF>, CustomMatrix<T,AF,PF,SO2> >
{
   typedef typename AddTrait< MT, CustomMatrix<T,AF,PF,SO2> >::Type  Type;
};

template< typename T, bool AF, bool PF, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< CustomMatrix<T,AF,PF,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< CustomMatrix<T,AF,PF,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct AddTrait< StrictlyUpperMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   typedef typename AddTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct AddTrait< CompressedMatrix<T,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename AddTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct AddTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< HermitianMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< LowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UniLowerMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< UniLowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyLowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename AddTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UpperMatrix<MT2,SO2,DF2> >
{
   typedef UpperMatrix< typename AddTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< UpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef UpperMatrix< typename AddTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   typedef UniUpperMatrix< typename AddTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< UniUpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef UniUpperMatrix< typename AddTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct AddTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename AddTrait<MT1,MT2>::Type >  Type;
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
struct SubTrait< StrictlyUpperMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename SubTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< StaticMatrix<T,M,N,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct SubTrait< StrictlyUpperMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename SubTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< HybridMatrix<T,M,N,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SubTrait< StrictlyUpperMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   typedef typename SubTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< DynamicMatrix<T,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool AF, bool PF, bool SO2 >
struct SubTrait< StrictlyUpperMatrix<MT,SO1,DF>, CustomMatrix<T,AF,PF,SO2> >
{
   typedef typename SubTrait< MT, CustomMatrix<T,AF,PF,SO2> >::Type  Type;
};

template< typename T, bool AF, bool PF, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< CustomMatrix<T,AF,PF,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< CustomMatrix<T,AF,PF,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SubTrait< StrictlyUpperMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   typedef typename SubTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< CompressedMatrix<T,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename SubTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct SubTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< HermitianMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< LowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UniLowerMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniLowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyLowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename SubTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UpperMatrix<MT2,SO2,DF2> >
{
   typedef UpperMatrix< typename SubTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef UpperMatrix< typename SubTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   typedef UpperMatrix< typename SubTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef UniUpperMatrix< typename SubTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename SubTrait<MT1,MT2>::Type >  Type;
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
struct MultTrait< StrictlyUpperMatrix<MT,SO,DF>, T, typename EnableIf< IsNumeric<T> >::Type >
{
   typedef StrictlyUpperMatrix< typename MultTrait<MT,T>::Type >  Type;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< T, StrictlyUpperMatrix<MT,SO,DF>, typename EnableIf< IsNumeric<T> >::Type >
{
   typedef StrictlyUpperMatrix< typename MultTrait<T,MT>::Type >  Type;
};

template< typename MT, bool SO, bool DF, typename T, size_t N >
struct MultTrait< StrictlyUpperMatrix<MT,SO,DF>, StaticVector<T,N,false> >
{
   typedef typename MultTrait< MT, StaticVector<T,N,false> >::Type  Type;
};

template< typename T, size_t N, typename MT, bool SO, bool DF >
struct MultTrait< StaticVector<T,N,true>, StrictlyUpperMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< StaticVector<T,N,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, typename T, size_t N >
struct MultTrait< StrictlyUpperMatrix<MT,SO,DF>, HybridVector<T,N,false> >
{
   typedef typename MultTrait< MT, HybridVector<T,N,false> >::Type  Type;
};

template< typename T, size_t N, typename MT, bool SO, bool DF >
struct MultTrait< HybridVector<T,N,true>, StrictlyUpperMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< HybridVector<T,N,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, typename T >
struct MultTrait< StrictlyUpperMatrix<MT,SO,DF>, DynamicVector<T,false> >
{
   typedef typename MultTrait< MT, DynamicVector<T,false> >::Type  Type;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< DynamicVector<T,true>, StrictlyUpperMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< DynamicVector<T,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, typename T, bool AF, bool PF >
struct MultTrait< StrictlyUpperMatrix<MT,SO,DF>, CustomVector<T,AF,PF,false> >
{
   typedef typename MultTrait< MT, CustomVector<T,AF,PF,false> >::Type  Type;
};

template< typename T, bool AF, bool PF, typename MT, bool SO, bool DF >
struct MultTrait< CustomVector<T,AF,PF,true>, StrictlyUpperMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< CustomVector<T,AF,PF,true>, MT >::Type  Type;
};

template< typename MT, bool SO, bool DF, typename T >
struct MultTrait< StrictlyUpperMatrix<MT,SO,DF>, CompressedVector<T,false> >
{
   typedef typename MultTrait< MT, CompressedVector<T,false> >::Type  Type;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< CompressedVector<T,true>, StrictlyUpperMatrix<MT,SO,DF> >
{
   typedef typename MultTrait< CompressedVector<T,true>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< StrictlyUpperMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   typedef typename MultTrait< MT, StaticMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< StaticMatrix<T,M,N,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< StaticMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< StrictlyUpperMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   typedef typename MultTrait< MT, HybridMatrix<T,M,N,SO2> >::Type  Type;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< HybridMatrix<T,M,N,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< HybridMatrix<T,M,N,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct MultTrait< StrictlyUpperMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   typedef typename MultTrait< MT, DynamicMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< DynamicMatrix<T,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< DynamicMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool AF, bool PF, bool SO2 >
struct MultTrait< StrictlyUpperMatrix<MT,SO1,DF>, CustomMatrix<T,AF,PF,SO2> >
{
   typedef typename MultTrait< MT, CustomMatrix<T,AF,PF,SO2> >::Type  Type;
};

template< typename T, bool AF, bool PF, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< CustomMatrix<T,AF,PF,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< CustomMatrix<T,AF,PF,SO1>, MT >::Type  Type;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct MultTrait< StrictlyUpperMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   typedef typename MultTrait< MT, CompressedMatrix<T,SO2> >::Type  Type;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< CompressedMatrix<T,SO1>, StrictlyUpperMatrix<MT,SO2,DF> >
{
   typedef typename MultTrait< CompressedMatrix<T,SO1>, MT >::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct MultTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< HermitianMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< LowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UniLowerMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniLowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyLowerMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef typename MultTrait<MT1,MT2>::Type  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename MultTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename MultTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename MultTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename MultTrait<MT1,MT2>::Type >  Type;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename MultTrait<MT1,MT2>::Type >  Type;
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
struct DivTrait< StrictlyUpperMatrix<MT,SO,DF>, T, typename EnableIf< IsNumeric<T> >::Type >
{
   typedef StrictlyUpperMatrix< typename DivTrait<MT,T>::Type >  Type;
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
struct MathTrait< StrictlyUpperMatrix<MT1,SO1,DF1>, StrictlyUpperMatrix<MT2,SO2,DF2> >
{
   typedef StrictlyUpperMatrix< typename MathTrait<MT1,MT2>::HighType >  HighType;
   typedef StrictlyUpperMatrix< typename MathTrait<MT1,MT2>::LowType  >  LowType;
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
struct SubmatrixTrait< StrictlyUpperMatrix<MT,SO,DF> >
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
struct RowTrait< StrictlyUpperMatrix<MT,SO,DF> >
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
struct ColumnTrait< StrictlyUpperMatrix<MT,SO,DF> >
{
   typedef typename ColumnTrait<MT>::Type  Type;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
