//=================================================================================================
/*!
//  \file blaze/math/adaptors/UniUpperMatrix.h
//  \brief Header file for the implementation of a upper unitriangular matrix adaptor
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
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

#ifndef _BLAZE_MATH_ADAPTORS_UNIUPPERMATRIX_H_
#define _BLAZE_MATH_ADAPTORS_UNIUPPERMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/adaptors/uniuppermatrix/BaseTemplate.h>
#include <blaze/math/adaptors/uniuppermatrix/Dense.h>
#include <blaze/math/adaptors/uniuppermatrix/Sparse.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/BinaryMapTrait.h>
#include <blaze/math/traits/ColumnsTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DeclDiagTrait.h>
#include <blaze/math/traits/DeclHermTrait.h>
#include <blaze/math/traits/DeclLowTrait.h>
#include <blaze/math/traits/DeclSymTrait.h>
#include <blaze/math/traits/DeclUppTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowsTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/UnaryMapTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/YieldsIdentity.h>
#include <blaze/math/typetraits/YieldsUniUpper.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/TrueType.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/Unused.h>


namespace blaze {

//=================================================================================================
//
//  UNIUPPERMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name UniUpperMatrix operators */
//@{
template< typename MT, bool SO, bool DF >
inline void reset( UniUpperMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline void reset( UniUpperMatrix<MT,SO,DF>& m, size_t i );

template< typename MT, bool SO, bool DF >
inline void clear( UniUpperMatrix<MT,SO,DF>& m );

template< bool RF, typename MT, bool SO, bool DF >
inline bool isDefault( const UniUpperMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline bool isIntact( const UniUpperMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
inline void swap( UniUpperMatrix<MT,SO,DF>& a, UniUpperMatrix<MT,SO,DF>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given uniupper matrix.
// \ingroup uniupper_matrix
//
// \param m The uniupper matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( UniUpperMatrix<MT,SO,DF>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the specified row/column of the given uniupper matrix.
// \ingroup uniupper_matrix
//
// \param m The uniupper matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given uniupper matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( UniUpperMatrix<MT,SO,DF>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given uniupper matrix.
// \ingroup uniupper_matrix
//
// \param m The uniupper matrix to be cleared.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void clear( UniUpperMatrix<MT,SO,DF>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given resizable uniupper matrix is in default state.
// \ingroup uniupper_matrix
//
// \param m The uniupper matrix to be tested for its default state.
// \return \a true in case the given matrix is in default state, \a false otherwise.
//
// This function checks whether the resizable upper unitriangular matrix is in default state.
*/
template< bool RF      // Relaxation flag
        , typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isDefault_backend( const UniUpperMatrix<MT,SO,DF>& m, TrueType )
{
   return ( m.rows() == 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given fixed-size uniupper matrix is in default state.
// \ingroup uniupper_matrix
//
// \param m The uniupper matrix to be tested for its default state.
// \return \a true in case the given matrix is in default state, \a false otherwise.
//
// This function checks whether the fixed-size upper unitriangular matrix is in default state.
*/
template< bool RF      // Relaxation flag
        , typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isDefault_backend( const UniUpperMatrix<MT,SO,DF>& m, FalseType )
{
   return isIdentity<RF>( m );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given uniupper matrix is in default state.
// \ingroup uniupper_matrix
//
// \param m The uniupper matrix to be tested for its default state.
// \return \a true in case the given matrix is component-wise zero, \a false otherwise.
//
// This function checks whether the upper unitriangular matrix is in default state. The following
// example demonstrates the use of the \a isDefault function:

   \code
   using blaze::DynamicMatrix;
   using blaze::UniUpperMatrix;
   using blaze::rowMajor;

   UniUpperMatrix< DynamicMatrix<int,rowMajor> > A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( A ) ) { ... }
   \endcode
*/
template< bool RF      // Relaxation flag
        , typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isDefault( const UniUpperMatrix<MT,SO,DF>& m )
{
   return isDefault_backend<RF>( m, typename IsResizable<MT>::Type() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given uniupper matrix are intact.
// \ingroup uniupper_matrix
//
// \param m The uniupper matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the uniupper matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   using blaze::DynamicMatrix;
   using blaze::UniUpperMatrix;

   UniUpperMatrix< DynamicMatrix<int> > A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isIntact( const UniUpperMatrix<MT,SO,DF>& m )
{
   return m.isIntact();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
// \ingroup uniupper_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void swap( UniUpperMatrix<MT,SO,DF>& a, UniUpperMatrix<MT,SO,DF>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given uniupper dense matrix.
// \ingroup uniupper_matrix
//
// \param m The uniupper dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given uniupper dense matrix by means of the specified matrix inversion
// algorithm \c IF. The The inversion fails if the given matrix is singular and not invertible.
// In this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a m may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert( UniUpperMatrix<MT,SO,true>& m )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   if( IF == asLower || IF == asUniLower ) {
      BLAZE_INTERNAL_ASSERT( isIdentity( m ), "Violation of preconditions detected" );
      return;
   }

   constexpr InversionFlag flag( ( IF == byLU || IF == asGeneral || IF == asUpper || IF == asUniUpper)
                                 ? ( asUniUpper )
                                 : ( asDiagonal ) );

   invert<flag>( derestrict( m ) );

   BLAZE_INTERNAL_ASSERT( isIntact( m ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief LU decomposition of the given uniupper dense matrix.
// \ingroup upper_matrix
//
// \param A The uniupper matrix to be decomposed.
// \param L The resulting lower triangular matrix.
// \param U The resulting upper triangular matrix.
// \param P The resulting permutation matrix.
// \return void
//
// This function performs the dense matrix (P)LU decomposition of a uniupper n-by-n matrix. The
// resulting decomposition is written to the three distinct matrices \c L, \c U, and \c P, which
// are resized to the correct dimensions (if possible and necessary).
//
// \note The LU decomposition will never fail, even for singular matrices. However, in case of a
// singular matrix the resulting decomposition cannot be used for a matrix inversion or solving
// a linear system of equations.
*/
template< typename MT1, bool SO1, typename MT2, typename MT3, typename MT4, bool SO2 >
inline void lu( const UniUpperMatrix<MT1,SO1,true>& A, DenseMatrix<MT2,SO1>& L,
                DenseMatrix<MT3,SO1>& U, Matrix<MT4,SO2>& P )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( MT3 );

   using ET2 = ElementType_t<MT2>;
   using ET4 = ElementType_t<MT4>;

   const size_t n( (~A).rows() );

   decltype(auto) L2( derestrict( ~L ) );

   (~U) = A;

   resize( ~L, n, n );
   reset( L2 );

   resize( ~P, n, n );
   reset( ~P );

   for( size_t i=0UL; i<n; ++i ) {
      L2(i,i)   = ET2(1);
      (~P)(i,i) = ET4(1);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of an uniupper matrix.
// \ingroup matrix
//
// \param mat The target uniupper matrix.
// \param i The row index of the element to be set.
// \param j The column index of the element to be set.
// \param value The value to be set to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool trySet( const UniUpperMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < (~mat).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < (~mat).columns(), "Invalid column access index" );

   UNUSED_PARAMETER( mat );

   if( i < j )
      return true;
   else if( i == j )
      return isOne( value );
   else
      return isDefault( value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of an uniupper matrix.
// \ingroup matrix
//
// \param mat The target uniupper matrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be added to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool tryAdd( const UniUpperMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < (~mat).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < (~mat).columns(), "Invalid column access index" );

   UNUSED_PARAMETER( mat );

   if( i < j )
      return true;
   else
      return isDefault( value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of an uniupper matrix.
// \ingroup matrix
//
// \param mat The target uniupper matrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be subtracted from the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool trySub( const UniUpperMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   return tryAdd( mat, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of an uniupper matrix.
// \ingroup matrix
//
// \param mat The target uniupper matrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The factor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool tryMult( const UniUpperMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < (~mat).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < (~mat).columns(), "Invalid column access index" );

   UNUSED_PARAMETER( mat );

   return ( i != j || isOne( value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of an uniupper matrix.
// \ingroup matrix
//
// \param mat The target uniupper matrix.
// \param row The index of the first row of the range to be multiplied.
// \param column The index of the first column of the range to be multiplied.
// \param m The number of rows of the range to be multiplied.
// \param n The number of columns of the range to be multiplied.
// \param value The factor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const UniUpperMatrix<MT,SO,DF>& mat, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (~mat).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (~mat).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (~mat).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (~mat).columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( mat );

   return ( row >= column + n ) || ( column >= row + m ) || isOne( value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of an uniupper matrix.
// \ingroup matrix
//
// \param mat The target uniupper matrix.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The divisor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool tryDiv( const UniUpperMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   return tryMult( mat, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of an uniupper matrix.
// \ingroup matrix
//
// \param mat The target uniupper matrix.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The divisor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const UniUpperMatrix<MT,SO,DF>& mat, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (~mat).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (~mat).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (~mat).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (~mat).columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( mat );

   return ( row >= column + n ) || ( column >= row + m ) || isOne( value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   if( column >= row + (~rhs).size() )
      return true;

   const bool containsDiagonal( column >= row );
   const size_t ibegin( ( !containsDiagonal )?( 0UL ):( column - row + 1UL ) );

   if( containsDiagonal && !isOne( (~rhs)[column-row] ) )
      return false;

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
/*!\brief Predict invariant violations by the assignment of a dense vector to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   if( row < column )
      return true;

   const bool containsDiagonal( row < column + (~rhs).size() );
   const size_t iend( min( row - column, (~rhs).size() ) );

   for( size_t i=0UL; i<iend; ++i ) {
      if( !isDefault( (~rhs)[i] ) )
         return false;
   }

   if( containsDiagonal && !isOne( (~rhs)[iend] ) )
      return false;

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to the band of an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side dense vector to be assigned.
// \param band The index of the band the right-hand side vector is assigned to.
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
        , typename VT  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag of the right-hand side dense vector
inline bool tryAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const DenseVector<VT,TF>& rhs,
                       ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs, row, column );

   if( band == 0L ) {
      for( size_t i=0UL; i<(~rhs).size(); ++i ) {
         if( !isOne( (~rhs)[i] ) )
            return false;
      }
   }
   else if( band < 0L ) {
      for( size_t i=0UL; i<(~rhs).size(); ++i ) {
         if( !isDefault( (~rhs)[i] ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse vector to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename VT::ConstIterator;

   if( column >= row + (~rhs).size() )
      return true;

   const bool containsDiagonal( column >= row );
   const size_t index( ( containsDiagonal )?( column - row ):( 0UL ) );
   const RhsIterator last( (~rhs).end() );
   RhsIterator element( (~rhs).lowerBound( index ) );

   if( containsDiagonal ) {
      if( element == last || element->index() != index || !isOne( element->value() ) )
         return false;
      ++element;
   }

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
/*!\brief Predict invariant violations by the assignment of a sparse vector to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename VT::ConstIterator;

   if( row < column )
      return true;

   const bool containsDiagonal( row < column + (~rhs).size() );
   const size_t index( row - column );
   const RhsIterator last( (~rhs).lowerBound( index ) );

   if( containsDiagonal ) {
      if( last == (~rhs).end() || last->index() != index || !isOne( last->value() ) )
         return false;
   }

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
/*!\brief Predict invariant violations by the assignment of a sparse vector to the band of an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side sparse vector to be assigned.
// \param band The index of the band the right-hand side vector is assigned to.
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
        , typename VT  // Type of the right-hand side sparse vector
        , bool TF >    // Transpose flag of the right-hand side sparse vector
inline bool tryAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const SparseVector<VT,TF>& rhs,
                       ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs, row, column );

   if( band == 0L ) {
      if( (~rhs).nonZeros() != (~rhs).size() )
         return false;
      for( const auto& element : ~rhs ) {
         if( !isOne( element.value() ) )
            return false;
      }
   }
   else if( band < 0L ) {
      for( const auto& element : ~rhs ) {
         if( !isDefault( element.value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense matrix to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                       const DenseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
      return true;

   const size_t ibegin( ( column < row )?( 0UL ):( column - row ) );

   for( size_t i=ibegin; i<M; ++i )
   {
      const size_t jend( min( row + i - column, N ) );

      for( size_t j=0UL; j<jend; ++j ) {
         if( !isDefault( (~rhs)(i,j) ) )
            return false;
      }

      const bool containsDiagonal( row + i < column + N );

      if( containsDiagonal && !isOne( (~rhs)(i,jend) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense matrix to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                       const DenseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
      return true;

   const size_t jend( min( row + M - column, N ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column + j >= row );

      if( containsDiagonal && !isOne( (~rhs)(column+j-row,j) ) )
         return false;

      const size_t ibegin( ( containsDiagonal )?( column + j - row + 1UL ):( 0UL ) );

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
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                       const SparseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename MT2::ConstIterator;

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
      return true;

   const size_t ibegin( ( column < row )?( 0UL ):( column - row ) );

   for( size_t i=ibegin; i<M; ++i )
   {
      const bool containsDiagonal( row + i < column + N );

      const size_t index( row + i - column );
      const RhsIterator last( (~rhs).lowerBound( i, min( index, N ) ) );

      if( containsDiagonal ) {
         if( last == (~rhs).end(i) || ( last->index() != index ) || !isOne( last->value() ) )
            return false;
      }

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
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool tryAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                       const SparseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename MT2::ConstIterator;

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
      return true;

   const size_t jend( min( row + M - column, N ) );

   for( size_t j=0UL; j<jend; ++j )
   {
      const bool containsDiagonal( column + j >= row );
      const size_t index( ( containsDiagonal )?( column + j - row ):( 0UL ) );

      const RhsIterator last( (~rhs).end(j) );
      RhsIterator element( (~rhs).lowerBound( index, j ) );

      if( containsDiagonal ) {
         if( element == last || ( element->index() != index ) || !isOne( element->value() ) )
            return false;
         ++element;
      }

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
/*!\brief Predict invariant violations by the addition assignment of a dense vector to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side dense vector to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                          const DenseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

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
/*!\brief Predict invariant violations by the addition assignment of a dense vector to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side dense vector to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                          const DenseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

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
/*!\brief Predict invariant violations by the addition assignment of a dense vector to the band
//        of an uniupper matrix.
// \ingroup uniupper_matrix
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side dense vector to be added.
// \param band The index of the band the right-hand side vector is assigned to.
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
        , typename VT  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag of the right-hand side dense vector
inline bool tryAddAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const DenseVector<VT,TF>& rhs,
                          ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs, row, column );

   if( band <= 0L ) {
      for( size_t i=0UL; i<(~rhs).size(); ++i ) {
         if( !isDefault( (~rhs)[i] ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a sparse vector to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side sparse vector to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                          const SparseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename VT::ConstIterator;

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
/*!\brief Predict invariant violations by the addition assignment of a sparse vector to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side sparse vector to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                          const SparseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename VT::ConstIterator;

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
/*!\brief Predict invariant violations by the addition assignment of a sparse vector to the band
//        of an uniupper matrix.
// \ingroup uniupper_matrix
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side sparse vector to be added.
// \param band The index of the band the right-hand side vector is assigned to.
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
        , typename VT  // Type of the right-hand side sparse vector
        , bool TF >    // Transpose flag of the right-hand side sparse vector
inline bool tryAddAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const SparseVector<VT,TF>& rhs,
                          ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs, row, column );

   if( band <= 0L ) {
      for( const auto& element : ~rhs ) {
         if( !isDefault( element.value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a dense matrix to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side dense matrix to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                          const DenseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
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
/*!\brief Predict invariant violations by the addition assignment of a dense matrix to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side dense matrix to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                          const DenseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
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
/*!\brief Predict invariant violations by the addition assignment of a sparse matrix to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side sparse matrix to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                          const SparseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename MT2::ConstIterator;

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
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
/*!\brief Predict invariant violations by the addition assignment of a sparse matrix to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side sparse matrix to be added.
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
inline bool tryAddAssign( const UniUpperMatrix<MT1,SO,DF>& lhs,
                          const SparseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   using RhsIterator = typename MT2::ConstIterator;

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( column >= row + M )
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
/*!\brief Predict invariant violations by the subtraction assignment of a vector to an uniupper
//        matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool trySubAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAddAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to the band of
//        an uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side vector to be subtracted.
// \param band The index of the band the right-hand side vector is assigned to.
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
inline bool trySubAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const Vector<VT,TF>& rhs,
                          ptrdiff_t band, size_t row, size_t column )
{
   return tryAddAssign( lhs, ~rhs, band, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to an uniupper
//        matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
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
inline bool trySubAssign( const UniUpperMatrix<MT1,SO1,DF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAddAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side vector to be multiplied.
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
        , typename VT >  // Type of the right-hand side vector
inline bool tryMultAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                           const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );

   UNUSED_PARAMETER( lhs );

   return ( column < row || (~rhs).size() <= column - row || isOne( (~rhs)[column-row] ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to an
//        uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side vector to be multiplied.
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
        , typename VT >  // Type of the right-hand side vector
inline bool tryMultAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                           const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   return ( row < column || (~rhs).size() <= row - column || isOne( (~rhs)[row-column] ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a dense vector to the
//        band of an uniupper matrix.
// \ingroup uniupper_matrix
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side dense vector to be multiplied.
// \param band The index of the band the right-hand side vector is assigned to.
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
        , typename VT  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag of the right-hand side dense vector
inline bool tryMultAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const DenseVector<VT,TF>& rhs,
                           ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs, row, column );

   if( band == 0L ) {
      for( size_t i=0UL; i<(~rhs).size(); ++i ) {
         if( !isOne( (~rhs)[i] ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a sparse vector to
//        the band of an uniupper matrix.
// \ingroup uniupper_matrix
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side sparse vector to be multiplied.
// \param band The index of the band the right-hand side vector is assigned to.
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
        , typename VT  // Type of the right-hand side sparse vector
        , bool TF >    // Transpose flag of the right-hand side sparse vector
inline bool tryMultAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const SparseVector<VT,TF>& rhs,
                           ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).size() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs, row, column );

   if( band == 0L ) {
      if( (~rhs).nonZeros() != (~rhs).size() )
         return false;
      for( const auto& element : ~rhs ) {
         if( !isOne( element.value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the Schur product assignment of a matrix to an uniupper
//        matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side matrix for the Schur product.
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
inline bool trySchurAssign( const UniUpperMatrix<MT1,SO1,DF>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (~rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (~rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   UNUSED_PARAMETER( lhs );

   const size_t M( (~rhs).rows()    );
   const size_t N( (~rhs).columns() );

   if( ( row >= column + N ) || ( column >= row + M ) )
      return true;

   size_t i( row < column ? column - row : 0UL );
   size_t j( column < row ? row - column : 0UL );

   for( ; i<M && j<N; ++i, ++j )
   {
      if( !isOne( (~rhs)(i,j) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to an uniupper matrix.
// \ingroup uniupper_matrix
//
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side vector divisor.
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
inline bool tryDivAssign( const UniUpperMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryMultAssign( lhs, ~rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to the band of
//        an uniupper matrix.
// \ingroup uniupper_matrix
// \param lhs The target left-hand side uniupper matrix.
// \param rhs The right-hand side vector divisor.
// \param band The index of the band the right-hand side vector is assigned to.
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
inline bool tryDivAssign( const UniUpperMatrix<MT,SO,DF>& lhs, const Vector<VT,TF>& rhs,
                          ptrdiff_t band, size_t row, size_t column )
{
   return tryMultAssign( lhs, rhs, band, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the instance without the access restrictions to the lower part.
// \ingroup math_shims
//
// \param m The uniupper matrix to be derestricted.
// \return Reference to the matrix without access restrictions.
//
// This function returns a reference to the given uniupper matrix instance that has no access
// restrictions to the lower part of the matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline MT& derestrict( UniUpperMatrix<MT,SO,DF>& m )
{
   return m.matrix_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct Size< UniUpperMatrix<MT,SO,DF>, 0UL >
   : public Size<MT,0UL>
{};

template< typename MT, bool SO, bool DF >
struct Size< UniUpperMatrix<MT,SO,DF>, 1UL >
   : public Size<MT,1UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAXSIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct MaxSize< UniUpperMatrix<MT,SO,DF>, 0UL >
   : public MaxSize<MT,0UL>
{};

template< typename MT, bool SO, bool DF >
struct MaxSize< UniUpperMatrix<MT,SO,DF>, 1UL >
   : public MaxSize<MT,1UL>
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
struct IsSquare< UniUpperMatrix<MT,SO,DF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsUniUpper< UniUpperMatrix<MT,SO,DF> >
   : public TrueType
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
struct IsAdaptor< UniUpperMatrix<MT,SO,DF> >
   : public TrueType
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
struct IsRestricted< UniUpperMatrix<MT,SO,DF> >
   : public TrueType
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
struct HasConstDataAccess< UniUpperMatrix<MT,SO,true> >
   : public TrueType
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
struct IsAligned< UniUpperMatrix<MT,SO,DF> >
   : public IsAligned<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCONTIGUOUS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsContiguous< UniUpperMatrix<MT,SO,DF> >
   : public IsContiguous<MT>
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
struct IsPadded< UniUpperMatrix<MT,SO,DF> >
   : public IsPadded<MT>
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
struct IsResizable< UniUpperMatrix<MT,SO,DF> >
   : public IsResizable<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSHRINKABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsShrinkable< UniUpperMatrix<MT,SO,DF> >
   : public IsShrinkable<MT>
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
struct RemoveAdaptor< UniUpperMatrix<MT,SO,DF> >
{
   using Type = MT;
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
template< typename T1, typename T2 >
struct AddTraitEval1< T1, T2
                    , EnableIf_t< ( IsUniUpper_v<T1> && IsStrictlyUpper_v<T2> ) ||
                                  ( IsStrictlyUpper_v<T1> && IsUniUpper_v<T2> ) > >
{
   using Type = UniUpperMatrix< typename AddTraitEval2<T1,T2>::Type >;
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
struct SubTrait< UniUpperMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   using Type = SubTrait_t< MT, StaticMatrix<T,M,N,SO2> >;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< StaticMatrix<T,M,N,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = SubTrait_t< StaticMatrix<T,M,N,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct SubTrait< UniUpperMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   using Type = SubTrait_t< MT, HybridMatrix<T,M,N,SO2> >;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< HybridMatrix<T,M,N,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = SubTrait_t< HybridMatrix<T,M,N,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SubTrait< UniUpperMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   using Type = SubTrait_t< MT, DynamicMatrix<T,SO2> >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< DynamicMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = SubTrait_t< DynamicMatrix<T,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool AF, bool PF, bool SO2 >
struct SubTrait< UniUpperMatrix<MT,SO1,DF>, CustomMatrix<T,AF,PF,SO2> >
{
   using Type = SubTrait_t< MT, CustomMatrix<T,AF,PF,SO2> >;
};

template< typename T, bool AF, bool PF, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< CustomMatrix<T,AF,PF,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = SubTrait_t< CustomMatrix<T,AF,PF,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SubTrait< UniUpperMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   using Type = SubTrait_t< MT, CompressedMatrix<T,SO2> >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< CompressedMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = SubTrait_t< CompressedMatrix<T,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SubTrait< UniUpperMatrix<MT,SO1,DF>, IdentityMatrix<T,SO2> >
{
   using Type = StrictlyUpperMatrix< SubTrait_t< MT, IdentityMatrix<T,SO2> > >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SubTrait< IdentityMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = StrictlyUpperMatrix< SubTrait_t< IdentityMatrix<T,SO1>, MT > >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct SubTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< HermitianMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< LowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, UniLowerMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniLowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< StrictlyLowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = SubTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, UpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SubTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SubTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SubTrait< UniUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SubTrait_t<MT1,MT2> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SCHURTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct SchurTrait< UniUpperMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   using Type = UpperMatrix< SchurTrait_t< MT, StaticMatrix<T,M,N,SO2> > >;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SchurTrait< StaticMatrix<T,M,N,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = UpperMatrix< SchurTrait_t< StaticMatrix<T,M,N,SO1>, MT > >;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct SchurTrait< UniUpperMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   using Type = UpperMatrix< SchurTrait_t< MT, HybridMatrix<T,M,N,SO2> > >;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct SchurTrait< HybridMatrix<T,M,N,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = UpperMatrix< SchurTrait_t< HybridMatrix<T,M,N,SO1>, MT > >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SchurTrait< UniUpperMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   using Type = UpperMatrix< SchurTrait_t< MT, DynamicMatrix<T,SO2> > >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SchurTrait< DynamicMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = UpperMatrix< SchurTrait_t< DynamicMatrix<T,SO1>, MT > >;
};

template< typename MT, bool SO1, bool DF, typename T, bool AF, bool PF, bool SO2 >
struct SchurTrait< UniUpperMatrix<MT,SO1,DF>, CustomMatrix<T,AF,PF,SO2> >
{
   using Type = UpperMatrix< SchurTrait_t< MT, CustomMatrix<T,AF,PF,SO2> > >;
};

template< typename T, bool AF, bool PF, bool SO1, typename MT, bool SO2, bool DF >
struct SchurTrait< CustomMatrix<T,AF,PF,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = UpperMatrix< SchurTrait_t< CustomMatrix<T,AF,PF,SO1>, MT > >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SchurTrait< UniUpperMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   using Type = UpperMatrix< SchurTrait_t< MT, CompressedMatrix<T,SO2> > >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SchurTrait< CompressedMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = UpperMatrix< SchurTrait_t< CompressedMatrix<T,SO1>, MT > >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct SchurTrait< UniUpperMatrix<MT,SO1,DF>, IdentityMatrix<T,SO2> >
{
   using Type = IdentityMatrix< MultTrait_t< ElementType_t<MT>, T >, SO2 >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct SchurTrait< IdentityMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = IdentityMatrix< MultTrait_t< T, ElementType_t<MT> >, SO1 >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct SchurTrait< UniUpperMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   using Type = UpperMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UniUpperMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< HermitianMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UniUpperMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   using Type = DiagonalMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< LowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = DiagonalMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UniUpperMatrix<MT1,SO1,DF1>, UniLowerMatrix<MT2,SO2,DF2> >
{
   using Type = IdentityMatrix< MultTrait_t< ElementType_t<MT1>, ElementType_t<MT2> >, SO1 >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UniLowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = IdentityMatrix< MultTrait_t< ElementType_t<MT1>, ElementType_t<MT2> >, SO1 >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UniUpperMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   using Type = DiagonalMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< StrictlyLowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = DiagonalMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UniUpperMatrix<MT1,SO1,DF1>, UpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< SchurTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct SchurTrait< UniUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UniUpperMatrix< SchurTrait_t<MT1,MT2> >;
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
struct MultTrait< UniUpperMatrix<MT,SO,DF>, T, EnableIf_t< IsNumeric_v<T> > >
{
   using Type = UpperMatrix< MultTrait_t<MT,T> >;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< T, UniUpperMatrix<MT,SO,DF>, EnableIf_t< IsNumeric_v<T> > >
{
   using Type = UpperMatrix< MultTrait_t<T,MT> >;
};

template< typename MT, bool SO, bool DF, typename T, size_t N >
struct MultTrait< UniUpperMatrix<MT,SO,DF>, StaticVector<T,N,false> >
{
   using Type = MultTrait_t< MT, StaticVector<T,N,false> >;
};

template< typename T, size_t N, typename MT, bool SO, bool DF >
struct MultTrait< StaticVector<T,N,true>, UniUpperMatrix<MT,SO,DF> >
{
   using Type = MultTrait_t< StaticVector<T,N,true>, MT >;
};

template< typename MT, bool SO, bool DF, typename T, size_t N >
struct MultTrait< UniUpperMatrix<MT,SO,DF>, HybridVector<T,N,false> >
{
   using Type = MultTrait_t< MT, HybridVector<T,N,false> >;
};

template< typename T, size_t N, typename MT, bool SO, bool DF >
struct MultTrait< HybridVector<T,N,true>, UniUpperMatrix<MT,SO,DF> >
{
   using Type = MultTrait_t< HybridVector<T,N,true>, MT >;
};

template< typename MT, bool SO, bool DF, typename T >
struct MultTrait< UniUpperMatrix<MT,SO,DF>, DynamicVector<T,false> >
{
   using Type = MultTrait_t< MT, DynamicVector<T,false> >;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< DynamicVector<T,true>, UniUpperMatrix<MT,SO,DF> >
{
   using Type = MultTrait_t< DynamicVector<T,true>, MT >;
};

template< typename MT, bool SO, bool DF, typename T, bool AF, bool PF >
struct MultTrait< UniUpperMatrix<MT,SO,DF>, CustomVector<T,AF,PF,false> >
{
   using Type = MultTrait_t< MT, CustomVector<T,AF,PF,false> >;
};

template< typename T, bool AF, bool PF, typename MT, bool SO, bool DF >
struct MultTrait< CustomVector<T,AF,PF,true>, UniUpperMatrix<MT,SO,DF> >
{
   using Type = MultTrait_t< CustomVector<T,AF,PF,true>, MT >;
};

template< typename MT, bool SO, bool DF, typename T >
struct MultTrait< UniUpperMatrix<MT,SO,DF>, CompressedVector<T,false> >
{
   using Type = MultTrait_t< MT, CompressedVector<T,false> >;
};

template< typename T, typename MT, bool SO, bool DF >
struct MultTrait< CompressedVector<T,true>, UniUpperMatrix<MT,SO,DF> >
{
   using Type = MultTrait_t< CompressedVector<T,true>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< UniUpperMatrix<MT,SO1,DF>, StaticMatrix<T,M,N,SO2> >
{
   using Type = MultTrait_t< MT, StaticMatrix<T,M,N,SO2> >;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< StaticMatrix<T,M,N,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = MultTrait_t< StaticMatrix<T,M,N,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, size_t M, size_t N, bool SO2 >
struct MultTrait< UniUpperMatrix<MT,SO1,DF>, HybridMatrix<T,M,N,SO2> >
{
   using Type = MultTrait_t< MT, HybridMatrix<T,M,N,SO2> >;
};

template< typename T, size_t M, size_t N, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< HybridMatrix<T,M,N,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = MultTrait_t< HybridMatrix<T,M,N,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct MultTrait< UniUpperMatrix<MT,SO1,DF>, DynamicMatrix<T,SO2> >
{
   using Type = MultTrait_t< MT, DynamicMatrix<T,SO2> >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< DynamicMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = MultTrait_t< DynamicMatrix<T,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool AF, bool PF, bool SO2 >
struct MultTrait< UniUpperMatrix<MT,SO1,DF>, CustomMatrix<T,AF,PF,SO2> >
{
   using Type = MultTrait_t< MT, CustomMatrix<T,AF,PF,SO2> >;
};

template< typename T, bool AF, bool PF, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< CustomMatrix<T,AF,PF,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = MultTrait_t< CustomMatrix<T,AF,PF,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct MultTrait< UniUpperMatrix<MT,SO1,DF>, CompressedMatrix<T,SO2> >
{
   using Type = MultTrait_t< MT, CompressedMatrix<T,SO2> >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< CompressedMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = MultTrait_t< CompressedMatrix<T,SO1>, MT >;
};

template< typename MT, bool SO1, bool DF, typename T, bool SO2 >
struct MultTrait< UniUpperMatrix<MT,SO1,DF>, IdentityMatrix<T,SO2> >
{
   using Type = UniUpperMatrix< MultTrait_t< MT, IdentityMatrix<T,SO2> > >;
};

template< typename T, bool SO1, typename MT, bool SO2, bool DF >
struct MultTrait< IdentityMatrix<T,SO1>, UniUpperMatrix<MT,SO2,DF> >
{
   using Type = UniUpperMatrix< MultTrait_t< IdentityMatrix<T,SO1>, MT > >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2, bool NF >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, SymmetricMatrix<MT2,SO2,DF2,NF> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, bool NF, typename MT2, bool SO2, bool DF2 >
struct MultTrait< SymmetricMatrix<MT1,SO1,DF1,NF>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, HermitianMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< HermitianMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, LowerMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< LowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, UniLowerMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniLowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< StrictlyLowerMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = MultTrait_t<MT1,MT2>;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, UpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< MultTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UpperMatrix< MultTrait_t<MT1,MT2> >;
};

template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct MultTrait< UniUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UniUpperMatrix< MultTrait_t<MT1,MT2> >;
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
struct DivTrait< UniUpperMatrix<MT,SO,DF>, T, EnableIf_t< IsNumeric_v<T> > >
{
   using Type = UpperMatrix< DivTrait_t<MT,T> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UNARYMAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename OP >
struct UnaryMapTrait< MT, OP, EnableIf_t< YieldsUniUpper_v<OP,MT> &&
                                          !YieldsIdentity_v<OP,MT> > >
{
   using Type = UniUpperMatrix< UnaryMapTrait_t< RemoveAdaptor_t<MT>, OP > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BINARYMAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, typename OP >
struct BinaryMapTrait< MT1, MT2, OP, EnableIf_t< YieldsUniUpper_v<OP,MT1,MT2> &&
                                                 !YieldsIdentity_v<OP,MT1,MT2> > >
{
   using Type = UniUpperMatrix< BinaryMapTrait_t< RemoveAdaptor_t<MT1>, RemoveAdaptor_t<MT2>, OP > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSYMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclSymTrait< UniUpperMatrix<MT,SO,DF> >
{
   using Type = IdentityMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLHERMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclHermTrait< UniUpperMatrix<MT,SO,DF> >
{
   using Type = IdentityMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclLowTrait< UniUpperMatrix<MT,SO,DF> >
{
   using Type = IdentityMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclUppTrait< UniUpperMatrix<MT,SO,DF> >
{
   using Type = UniUpperMatrix<MT,SO,DF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLDIAGTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclDiagTrait< UniUpperMatrix<MT,SO,DF> >
{
   using Type = IdentityMatrix< ElementType_t<MT>, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct HighType< UniUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UniUpperMatrix< typename HighType<MT1,MT2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOWTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct LowType< UniUpperMatrix<MT1,SO1,DF1>, UniUpperMatrix<MT2,SO2,DF2> >
{
   using Type = UniUpperMatrix< typename LowType<MT1,MT2>::Type >;
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
template< typename MT, bool SO, bool DF, size_t... CSAs >
struct SubmatrixTrait< UniUpperMatrix<MT,SO,DF>, CSAs... >
{
   using Type = SubmatrixTrait_t<MT,CSAs...>;
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
template< typename MT, bool SO, bool DF, size_t... CRAs >
struct RowTrait< UniUpperMatrix<MT,SO,DF>, CRAs... >
{
   using Type = RowTrait_t<MT,CRAs...>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, size_t... CRAs >
struct RowsTrait< UniUpperMatrix<MT,SO,DF>, CRAs... >
{
   using Type = RowsTrait_t<MT,CRAs...>;
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
template< typename MT, bool SO, bool DF, size_t... CCAs >
struct ColumnTrait< UniUpperMatrix<MT,SO,DF>, CCAs... >
{
   using Type = ColumnTrait_t<MT,CCAs...>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, size_t... CCAs >
struct ColumnsTrait< UniUpperMatrix<MT,SO,DF>, CCAs... >
{
   using Type = ColumnsTrait_t<MT,CCAs...>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, ptrdiff_t... CBAs >
struct BandTrait< UniUpperMatrix<MT,SO,DF>, CBAs... >
{
   using Type = BandTrait_t<MT,CBAs...>;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
