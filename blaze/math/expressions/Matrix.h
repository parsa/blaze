//=================================================================================================
/*!
//  \file blaze/math/expressions/Matrix.h
//  \brief Header file for the Matrix base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATRIX_H_
#define _BLAZE_MATH_EXPRESSIONS_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Assert.h>
#include <blaze/util/logging/FunctionTrace.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup matrix Matrices
// \ingroup math
*/
/*!\brief Base class for matrices.
// \ingroup matrix
//
// The Matrix class is a base class for all dense and sparse matrix classes within the Blaze
// library. It provides an abstraction from the actual type of the matrix, but enables a
// conversion back to this type via the 'Curiously Recurring Template Pattern' (CRTP).
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order
struct Matrix
{
   //**Type definitions****************************************************************************
   typedef MT  MatrixType;  //!< Type of the matrix.
   //**********************************************************************************************

   //**Non-const conversion operator***************************************************************
   /*!\brief Conversion operator for non-constant matrices.
   //
   // \return Reference of the actual type of the matrix.
   */
   inline MatrixType& operator~() {
      return *static_cast<MatrixType*>( this );
   }
   //**********************************************************************************************

   //**Const conversion operator*******************************************************************
   /*!\brief Conversion operator for constant matrices.
   //
   // \return Constant reference of the actual type of the matrix.
   */
   inline const MatrixType& operator~() const {
      return *static_cast<const MatrixType*>( this );
   }
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Matrix global functions */
//@{
template< typename MT, bool SO >
inline size_t rows( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
inline size_t columns( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
inline size_t capacity( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
inline size_t capacity( const Matrix<MT,SO>& m, size_t i );

template< typename MT, bool SO >
inline size_t nonZeros( const Matrix<MT,SO>& m );

template< typename MT, bool SO >
inline size_t nonZeros( const Matrix<MT,SO>& m, size_t i );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void assign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void addAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void subAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline void multAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
inline bool isSame( const Matrix<MT1,SO1>& a, const Matrix<MT2,SO2>& b );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of rows of the matrix.
// \ingroup matrix
//
// \param m The given matrix.
// \return The number of rows of the matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t rows( const Matrix<MT,SO>& m )
{
   return (~m).rows();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
// \ingroup matrix
//
// \param m The given matrix.
// \return The number of columns of the matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t columns( const Matrix<MT,SO>& m )
{
   return (~m).columns();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
// \ingroup matrix
//
// \param m The given matrix.
// \return The capacity of the matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t capacity( const Matrix<MT,SO>& m )
{
   return (~m).capacity();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current capacity of the specified row/column.
// \ingroup matrix
//
// \param m The given matrix.
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// This function returns the current capacity of the specified row/column. In case the
// storage order is set to \a rowMajor the function returns the capacity of row \a i,
// in case the storage flag is set to \a columnMajor the function returns the capacity
// of column \a i.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t capacity( const Matrix<MT,SO>& m, size_t i )
{
   return (~m).capacity( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the total number of non-zero elements in the matrix
// \ingroup matrix
//
// \param m The given matrix.
// \return The number of non-zero elements in the dense matrix.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t nonZeros( const Matrix<MT,SO>& m )
{
   return (~m).nonZeros();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row/column.
// \ingroup matrix
//
// \param m The given matrix.
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// This function returns the current number of non-zero elements in the specified row/column.
// In case the storage order is set to \a rowMajor the function returns the number of non-zero
// elements in row \a i, in case the storage flag is set to \a columnMajor the function returns
// the number of non-zero elements in column \a i.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
inline size_t nonZeros( const Matrix<MT,SO>& m, size_t i )
{
   return (~m).nonZeros( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a matrix to a matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be assigned.
// \return void
//
// This function implements the default assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void assign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).columns(), "Invalid number of columns" );

   (~lhs).assign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a matrix to a matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be added.
// \return void
//
// This function implements the default addition assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void addAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).columns(), "Invalid number of columns" );

   (~lhs).addAssign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a matrix to matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be subtracted.
// \return void
//
// This function implements the default subtraction assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void subAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == (~rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).columns(), "Invalid number of columns" );

   (~lhs).subAssign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the multiplication assignment of a matrix to a matrix.
// \ingroup matrix
//
// \param lhs The target left-hand side matrix.
// \param rhs The right-hand side matrix to be multiplied.
// \return void
//
// This function implements the default multiplication assignment of a matrix to a matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void multAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (~lhs).columns() == (~rhs).rows(), "Invalid matrix sizes" );

   (~lhs).multAssign( ~rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the two given matrices represent the same observable state.
// \ingroup matrix
//
// \param a The first matrix to be tested for its state.
// \param b The second matrix to be tested for its state.
// \return \a true in case the two matrices share a state, \a false otherwise.
//
// The isSame function provides an abstract interface for testing if the two given matrices
// represent the same observable state. This happens for instance in case \c a and \c b refer
// to the same matrix or in case \c a and \c b are aliases for the same matrix. In case both
// matrices represent the same observable state, the function returns \a true, other it returns
// \a false.

   \code
   typedef blaze::DynamicMatrix<int>          MatrixType;
   typedef blaze::DenseSubmatrix<MatrixType>  SubmatrixType;

   MatrixType mat1( 4UL, 5UL );  // Setup of a 4x5 dynamic matrix
   MatrixType mat2( 4UL, 5UL );  // Setup of a second 4x5 dynamic matrix

   SubmatrixType sub1 = submatrix( mat1, 0UL, 0UL, 4UL, 5UL );  // Submatrix fully covering mat1
   SubmatrixType sub2 = submatrix( mat1, 1UL, 1UL, 2UL, 3UL );  // Submatrix partially covering mat1
   SubmatrixType sub3 = submatrix( mat1, 1UL, 1UL, 2UL, 3UL );  // Submatrix partially covering mat1

   isSame( mat1, mat1 );  // returns true since both objects refer to the same matrix
   isSame( mat1, mat2 );  // returns false since mat1 and mat2 are two different matrices
   isSame( mat1, sub1 );  // returns true since sub1 represents the same observable state as mat1
   isSame( mat1, sub3 );  // returns false since sub3 only covers part of mat1
   isSame( sub2, sub3 );  // returns true since sub1 and sub2 refer to exactly the same part of mat1
   isSame( sub1, sub3 );  // returns false since sub1 and sub3 refer to different parts of mat1
   \endcode
*/
template< typename MT1  // Type of the left-hand side matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool isSame( const Matrix<MT1,SO1>& a, const Matrix<MT2,SO2>& b )
{
   return ( IsSame<MT1,MT2>::value &&
            reinterpret_cast<const void*>( &a ) == reinterpret_cast<const void*>( &b ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
