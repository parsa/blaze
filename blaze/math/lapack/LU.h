//=================================================================================================
/*!
//  \file blaze/math/lapack/LU.h
//  \brief Header file for LAPACK LU decomposition functions
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

#ifndef _BLAZE_MATH_LAPACK_LU_H_
#define _BLAZE_MATH_LAPACK_LU_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/cast.hpp>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/constraints/Double.h>
#include <blaze/util/constraints/Float.h>
#include <blaze/util/Exception.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
extern "C" {

void sgetrf_( int* m, int* n, float*  a, int* lda, int* ipiv, int* info );
void dgetrf_( int* m, int* n, double* a, int* lda, int* ipiv, int* info );
void cgetrf_( int* m, int* n, float*  a, int* lda, int* ipiv, int* info );
void zgetrf_( int* m, int* n, double* a, int* lda, int* ipiv, int* info );

}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LAPACK LU DECOMPOSITION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LU decomposition functions */
//@{
inline void sgetrf( int* m, int* n, float* a, int* lda, int* ipiv, int* info );

inline void dgetrf( int* m, int* n, double* a, int* lda, int* ipiv, int* info );

inline void cgetrf( int* m, int* n, complex<float>* a, int* lda, int* ipiv, int* info );

inline void zgetrf( int* m, int* n, complex<double>* a, int* lda, int* ipiv, int* info );

template< typename MT, bool SO >
inline void sgetrf( DenseMatrix<MT,SO>& A, int* ipiv );

template< typename MT, bool SO >
inline void dgetrf( DenseMatrix<MT,SO>& A, int* ipiv );

template< typename MT, bool SO >
inline void cgetrf( DenseMatrix<MT,SO>& A, int* ipiv );

template< typename MT, bool SO >
inline void zgetrf( DenseMatrix<MT,SO>& A, int* ipiv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense single precision matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK sgetrf() function, which uses partial pivoting with row interchanges. The
// decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the sgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void sgetrf( int* m, int* n, float* a, int* lda, int* ipiv, int* info )
{
   sgetrf_( m, n, a, lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense double precision matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK dgetrf() function, which uses partial pivoting with row interchanges. The
// decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the dgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void dgetrf( int* m, int* n, double* a, int* lda, int* ipiv, int* info )
{
   dgetrf_( m, n, a, lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense single precision complex matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK cgetrf() function, which uses partial pivoting with row interchanges. The
// decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the cgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void cgetrf( int* m, int* n, complex<float>* a, int* lda, int* ipiv, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cgetrf_( m, n, reinterpret_cast<float*>( a ), lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense double precision complex matrix.
// \ingroup lapack
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param a Pointer to the first element of the matrix.
// \param lda The total number of elements between two rows/columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK zgetrf() function, which uses partial pivoting with row interchanges. The
// decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the zgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
inline void zgetrf( int* m, int* n, complex<double>* a, int* lda, int* ipiv, int* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   zgetrf_( m, n, reinterpret_cast<double*>( a ), lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense single precision matrix.
// \ingroup lapack
//
// \param A The matrix to be decomposed.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \return void
// \exception std::invalid_argument Decomposition of singular matrix failed.
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK sgetrf() function, which uses partial pivoting with row interchanges. Note
// that the function only works for general, non-adapted matrices with \c float element type. The
// attempt to call the function with adaptors or matrices of any other element type results in a
// compile time error!\n
//
// The decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed. The LU decomposition fails if \a A is a singular matrix, which cannot be inverted.
// In this case a \a std::std::invalid_argument exception is thrown.
//
// For more information on the sgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void sgetrf( DenseMatrix<MT,SO>& A, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT::ElementType );

   int m   ( boost::numeric_cast<int>( (~A).rows()    ) );
   int n   ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda ( boost::numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   sgetrf( &m, &n, (~A).data(), &lda, ipiv, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for LU decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Decomposition of singular matrix failed" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense double precision matrix.
// \ingroup lapack
//
// \param A The matrix to be decomposed.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \return void
// \exception std::invalid_argument Decomposition of singular matrix failed.
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK dgetrf() function, which uses partial pivoting with row interchanges. Note
// that the function only works for general, non-adapted matrices with \c double element type. The
// attempt to call the function with adaptors or matrices of any other element type results in a
// compile time error!\n
//
// The decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed. The LU decomposition fails if \a A is a singular matrix, which cannot be inverted.
// In this case a \a std::std::invalid_argument exception is thrown.
//
// For more information on the dgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void dgetrf( DenseMatrix<MT,SO>& A, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT::ElementType );

   int m   ( boost::numeric_cast<int>( (~A).rows()    ) );
   int n   ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda ( boost::numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   dgetrf( &m, &n, (~A).data(), &lda, ipiv, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for LU decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Decomposition of singular matrix failed" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense single precision complex matrix.
// \ingroup lapack
//
// \param A The matrix to be decomposed.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \return void
// \exception std::invalid_argument Decomposition of singular matrix failed.
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK cgetrf() function, which uses partial pivoting with row interchanges. Note
// that the function only works for general, non-adapted matrices with \c complex<float> element
// type. The attempt to call the function with adaptors or matrices of any other element type
// results in a compile time error!\n
//
// The decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed. The LU decomposition fails if \a A is a singular matrix, which cannot be inverted.
// In this case a \a std::std::invalid_argument exception is thrown.
//
// For more information on the cgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void cgetrf( DenseMatrix<MT,SO>& A, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_FLOAT_TYPE( typename MT::ElementType::value_type );

   int m   ( boost::numeric_cast<int>( (~A).rows()    ) );
   int n   ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda ( boost::numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   cgetrf( &m, &n, (~A).data(), &lda, ipiv, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for LU decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Decomposition of singular matrix failed" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense double precision complex matrix.
// \ingroup lapack
//
// \param A The matrix to be decomposed.
// \param ipiv Auxiliary array for the pivot indices; size >= min( M, N ).
// \return void
// \exception std::invalid_argument Decomposition of singular matrix failed.
//
// This function performs the dense matrix LU decomposition of a general \f$ M \times N \f$ matrix
// based on the LAPACK zgetrf() function, which uses partial pivoting with row interchanges. Note
// that the function only works for general, non-adapted matrices with \c complex<double> element
// type. The attempt to call the function with adaptors or matrices of any other element type
// results in a compile time error!\n
//
// The decomposition has the form

                          \f[ A = P \dot L \dot U, \f]\n

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: In case of a column-major
// matrix, \c L is stored in the lower part of \a A and \c U is stored in the upper part. The unit
// diagonal elements of \c L are not stored. In case \a A is a row-major matrix the result is
// transposed. The LU decomposition fails if \a A is a singular matrix, which cannot be inverted.
// In this case a \a std::std::invalid_argument exception is thrown.
//
// For more information on the zgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function does not provide any exception safety guarantee, i.e. in case an exception
// is thrown \a A may already have been modified.
// \note This function can only be used if the fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void zgetrf( DenseMatrix<MT,SO>& A, int* ipiv )
{
   using boost::numeric_cast;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( typename MT::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_DOUBLE_TYPE( typename MT::ElementType::value_type );

   int m   ( boost::numeric_cast<int>( (~A).rows()    ) );
   int n   ( boost::numeric_cast<int>( (~A).columns() ) );
   int lda ( boost::numeric_cast<int>( (~A).spacing() ) );
   int info( 0 );

   zgetrf( &m, &n, (~A).data(), &lda, ipiv, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for LU decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Decomposition of singular matrix failed" );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
